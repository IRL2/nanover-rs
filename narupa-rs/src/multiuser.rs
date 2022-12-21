use crate::broadcaster::Broadcaster;
use crate::proto::protocol::command::{CommandMessage, CommandReply};
use crate::proto::protocol::state::StateUpdate;
use crate::services::commands::Command;
use crate::state_broadcaster::StateBroadcaster;
use prost_types::value::Kind;
use prost_types::{ListValue, Struct, Value};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

pub struct RadialOrient {
    state: Arc<Mutex<StateBroadcaster>>,
}

impl RadialOrient {
    pub fn new(state: Arc<Mutex<StateBroadcaster>>) -> Self {
        RadialOrient { state }
    }
}

impl Command for RadialOrient {
    fn run(&self, input: CommandMessage) -> CommandReply {
        let default_radius = 1.0;
        let radius = extract_radius(&input, default_radius);
        let full_circle = std::f64::consts::PI * 2.0;
        let avatars = get_avatar_ids(&self.state);
        let count = avatars.len();
        let angles = (0..count).map(|i| i as f64 * full_circle / (count as f64));
        let updates = to_user_origin_update(avatars.iter().zip(angles).map(|(avatar, angle)| {
            let converted_angle = 0.5 * (-angle - full_circle / 4.0);
            (
                format!("user-origin.{avatar}"),
                [
                    (
                        "position",
                        vec![radius * angle.cos(), 0.0, radius * angle.sin()],
                    ),
                    (
                        "rotation",
                        vec![0.0, converted_angle.sin(), 0.0, converted_angle.cos()],
                    ),
                ],
            )
        }));
        self.state.lock().unwrap().send(updates).unwrap();
        CommandReply { result: None }
    }

    fn arguments(&self) -> Option<Struct> {
        Some(Struct {
            fields: BTreeMap::from([(
                String::from("radius"),
                Value {
                    kind: Some(Kind::NumberValue(1.0)),
                },
            )]),
        })
    }
}

fn get_avatar_ids(state: &Arc<Mutex<StateBroadcaster>>) -> Vec<String> {
    let state = state.lock().unwrap();
    state
        .iter()
        .filter_map(|(key, _)| {
            if key.starts_with("avatar.") {
                Some(String::from(key.split_once('.').unwrap().1))
            } else {
                None
            }
        })
        .collect()
}

fn extract_radius(input: &CommandMessage, default: f64) -> f64 {
    input.arguments.as_ref().map_or(default, |args| {
        args.fields.get("radius").map_or(default, |value| {
            if let Some(Kind::NumberValue(radius)) = value.kind {
                radius
            } else {
                default
            }
        })
    })
}

fn to_user_origin_update<'a>(
    input: impl Iterator<Item = (String, [(&'a str, Vec<f64>); 2])>,
) -> StateUpdate {
    StateUpdate {
        changed_keys: Some(Struct {
            fields: BTreeMap::from_iter(input.map(|(key, value)| {
                (
                    String::from(key),
                    orient_inner_struct_value(&value)
                )
            })),
        }),
    }
}

fn number_to_value(number: &f64) -> Value {
    Value{ kind: Some(Kind::NumberValue(*number)) }
}

fn list_of_numbers(values: &Vec<f64>) -> Value {
    Value {
        kind: Some(Kind::ListValue(ListValue {
            values: values
                .iter()
                .map(number_to_value)
                .collect(),
        })),
    }
}

fn orient_inner_struct_value(content: &[(&str, Vec<f64>)]) -> Value {
    Value {
        kind: Some(Kind::StructValue(Struct {
            fields: BTreeMap::from_iter(content.into_iter().map(|(key, value)| {
                (
                    String::from(*key),
                    list_of_numbers(value),
                )
            })),
        })),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(1, 1.0)]
    #[case(20, 4.0)]
    fn test_radial_orient_completeness(#[case] number_of_avatars: usize, #[case] radius: f64) {
        let state = Arc::new(Mutex::new(StateBroadcaster::new(None)));
        let mut avatars_ids = populate_avatars(&state, number_of_avatars);
        let command = RadialOrient::new(Arc::clone(&state));
        let arguments = build_arguments(radius);
        command.run(arguments);

        let mut orient_avatar_ids = get_avatar_ids_from_origins(&state);

        orient_avatar_ids.sort();
        avatars_ids.sort();

        assert_eq!(orient_avatar_ids, avatars_ids);
    }

    #[rstest]
    #[case(1, 1.0)]
    #[case(20, 4.0)]
    fn test_radial_orient_distance(#[case] number_of_avatars: usize, #[case] radius: f64) {
        let state = Arc::new(Mutex::new(StateBroadcaster::new(None)));
        populate_avatars(&state, number_of_avatars);
        let command = RadialOrient::new(Arc::clone(&state));
        let arguments = build_arguments(radius);
        command.run(arguments);

        let positions = get_positions_from_origins(&state);
        let mut distances = positions.iter().map(compute_distance);

        assert!(distances.all(|distance| distance == radius));
    }

    fn build_arguments(radius: f64) -> CommandMessage {
        CommandMessage {name: String::from("command/name"), arguments: Some(Struct{fields: BTreeMap::from([(String::from("radius"), number_to_value(&radius))])})}
    }

    fn populate_avatars(state: &Arc<Mutex<StateBroadcaster>>, number_of_avatars: usize) -> Vec<String> {
        let avatar_ids: Vec<String> = (0..number_of_avatars).map(|id| format!("{id}")).collect();
        let avatars = avatar_ids.iter().map(|id| (format!("avatar.{id}"), Value{kind: Some(Kind::StructValue(Struct{fields: BTreeMap::new()}))}));
        let update = StateUpdate { changed_keys: Some(Struct{fields: BTreeMap::from_iter(avatars)}) };
        state.lock().unwrap().send(update).unwrap();
        avatar_ids
    }

    fn get_avatar_ids_from_origins(state: &Arc<Mutex<StateBroadcaster>>) -> Vec<String> {
        state.lock().unwrap().iter().filter_map(|(key, _)| if key.starts_with("user-origin.") {Some(String::from(key.split_once('.').unwrap().1))} else {None}).collect()
    }

    fn get_positions_from_origins(state: &Arc<Mutex<StateBroadcaster>>) -> Vec<[f64; 3]> {
        state
            .lock()
            .unwrap()
            .iter()
            .filter_map(|(key, value)|
                if !key.starts_with("user-origin.") {None}
                else {extract_position(value)}
            )
            .collect()
    }

    fn extract_position(value: &Value) -> Option<[f64; 3]> {
        let Some(Kind::StructValue(ref inner_struct)) = value.kind else {
            return None;
        };
        let fields = &inner_struct.fields;
        let perhaps_positions_value = fields.get("positions");
        let Some(positions_value) = perhaps_positions_value else {
            return None;
        };
        let Some(Kind::ListValue(ref position_vector)) = positions_value.kind else {
            return None;
        };
        let perhaps_vector_numbers: Option<Vec<f64>> = position_vector
            .values
            .iter()
            .map(unpack_number)
            .collect();
        let Some(vector_numbers) = perhaps_vector_numbers else {
            return None;
        };
        vector_numbers
            .try_into()
            .ok()
    }

    fn unpack_number(value: &Value) -> Option<f64> {
        let Some(Kind::NumberValue(number)) = value.kind else {
            return None;
        };
        Some(number)
    }

    fn compute_distance(vector: &[f64; 3]) -> f64 {
        (vector[0].powf(2.0) + vector[1].powf(2.0) + vector[2].powf(2.0)).sqrt()
    }
}