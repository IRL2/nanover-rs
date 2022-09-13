use crate::simulation::{IMDInteraction, InteractionKind};
use prost_types::value::Kind;


pub fn read_state_interaction(state_interaction: &prost_types::Value) -> Result<IMDInteraction, ()> {
    // Extract the interaction content that is under several layers of enums.
    // Fail already if we cannot: it means the input does not have the expected format.
    let content = match &state_interaction.kind {
        Some(Kind::StructValue(inside)) => inside,
        _ => return Err(()),
    };

    let kind = content.fields.get("interaction_type");
    if kind.is_none() {
        return Err(());
    };
    let kind = kind.unwrap();
    let kind = match &kind.kind {
        Some(Kind::StringValue(inside)) => inside,
        _ => return Err(()),
    };
    let kind = if kind == "gaussian" {
        InteractionKind::GAUSSIAN
    } else if kind == "spring" {
        InteractionKind::HARMONIC
    } else {
        return Err(());
    };

    let max_force = content.fields.get("max_force");
    let max_force = match max_force {
        None => None,
        Some(inside) => match inside.kind {
            Some(Kind::NumberValue(value)) => Some(value),
            _ => return Err(()),
        },
    };

    let scale = content.fields.get("scale");
    let scale = match scale {
        None => 1.0,
        Some(inside) => match inside.kind {
            Some(Kind::NumberValue(value)) => value,
            _ => return Err(()),
        },
    };

    let particles = content.fields.get("particles");
    let particles: Vec<usize> = match particles {
        Some(inside) => {
            match &inside.kind {
                Some(Kind::ListValue(values)) => {
                    values
                        .values
                        .iter()
                        .filter_map(|v| v.kind.as_ref())
                        // We just ignore invalid values
                        .filter_map(|v| match v {
                            Kind::NumberValue(inner) => Some((*inner) as usize),
                            _ => None,
                        })
                        .collect()
                }
                _ => return Err(()),
            }
        }
        _ => return Err(()),
    };

    let position = content.fields.get("position");
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

    Ok(IMDInteraction::new(
        position, particles, kind, max_force, scale,
    ))
}
