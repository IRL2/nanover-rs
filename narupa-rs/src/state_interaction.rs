use crate::simulation::{IMDInteraction, InteractionKind};
use crate::state_broadcaster::StateBroadcaster;
use prost_types::{value::Kind, Struct};
use std::sync::{Arc, Mutex};

pub fn read_forces(state_clone: &Arc<Mutex<StateBroadcaster>>) -> Vec<IMDInteraction> {
    let state_interactions: Vec<IMDInteraction> = {
        let state = state_clone.lock().unwrap();
        let interaction_iter = state.iter();
        interaction_iter
            .filter(|kv| kv.0.starts_with("interaction."))
            .map(|kv| {
                let value = kv.1;
                println!("{:?}", kv);
                read_state_interaction(value)
            })
            .filter_map(|interaction| interaction.ok())
            .collect()
    };
    state_interactions
}

fn read_state_interaction(state_interaction: &prost_types::Value) -> Result<IMDInteraction, ()> {
    // Extract the interaction content that is under several layers of enums.
    // Fail already if we cannot: it means the input does not have the expected format.
    let content = match &state_interaction.kind {
        Some(Kind::StructValue(inside)) => inside,
        _ => return Err(()),
    };

    let kind = get_kind(content, "interaction_type")?;
    let max_force = get_number_or_error(content, "max_force")?;
    let scale = get_number_or_error(content, "scale")?.unwrap_or(1.0);
    let particles = get_particles(content, "particles")?;
    let position = get_position(content, "position")?;
    println!("Interaction with atoms {particles:?}");

    Ok(IMDInteraction::new(
        position, particles, kind, max_force, scale,
    ))
}

fn get_number_or_error(content: &Struct, key: &str) -> Result<Option<f64>, ()> {
    let content_value = content.fields.get(key);
    match content_value {
        None => Ok(None),
        Some(inside) => match inside.kind {
            Some(Kind::NumberValue(value)) => Ok(Some(value)),
            _ => Err(()),
        },
    }
}

fn get_kind(content: &Struct, key: &str) -> Result<InteractionKind, ()> {
    let kind = content.fields.get(key).ok_or(())?;
    let kind = match &kind.kind {
        Some(Kind::StringValue(inside)) => inside,
        _ => return Err(()),
    };
    if kind == "gaussian" {
        Ok(InteractionKind::GAUSSIAN)
    } else if kind == "spring" {
        Ok(InteractionKind::HARMONIC)
    } else {
        Err(())
    }
}

fn get_particles(content: &Struct, key: &str) -> Result<Vec<usize>, ()> {
    let particles = content.fields.get(key).ok_or(())?;
    match &particles.kind {
        Some(Kind::ListValue(values)) => {
            Ok(values
                .values
                .iter()
                .filter_map(|v| v.kind.as_ref())
                // We just ignore invalid values
                .filter_map(|v| match v {
                    Kind::NumberValue(inner) => Some((*inner) as usize),
                    _ => None,
                })
                .collect())
        }
        _ => Err(()),
    }
}

fn get_position(content: &Struct, key: &str) -> Result<[f64; 3], ()> {
    let position = content.fields.get(key);
    let position: Vec<f64> = match position {
        Some(inside) => match &inside.kind {
            Some(Kind::ListValue(values)) => values
                .values
                .iter()
                .filter_map(|v| v.kind.as_ref())
                .filter_map(|v| match v {
                    Kind::NumberValue(inner) => Some(*inner),
                    _ => None,
                })
                .collect(),
            _ => return Err(()),
        },
        _ => return Err(()),
    };
    if position.len() != 3 {
        return Err(());
    }
    let position: [f64; 3] = position.try_into().unwrap();
    Ok(position)
}
