use std::sync::{Arc, Mutex};
use crate::broadcaster::Broadcaster;
use crate::proto::protocol::state::StateUpdate;
use crate::proto::protocol::command::{CommandMessage, CommandReply};
use crate::state_broadcaster::StateBroadcaster;
use crate::services::commands::Command;
use prost_types::value::Kind;
use prost_types::{Struct, Value, ListValue};
use std::collections::BTreeMap;

pub struct RadialOrient {
    state: Arc<Mutex<StateBroadcaster>>,
}

impl RadialOrient {
    pub fn new(state: Arc<Mutex<StateBroadcaster>>) -> Self {
        RadialOrient{ state }
    }
}

impl Command for RadialOrient {
    fn run(&self, input: CommandMessage) -> CommandReply {
        let default_radius = 1.0;
        let radius = input.arguments
            .map_or(default_radius, |args| args
                .fields
                .get("radius")
                .map_or(default_radius, |value|
                    if let Some(Kind::NumberValue(radius)) = value.kind {
                        radius
                    } else {
                        default_radius
                    })
            );
        let full_circle = std::f64::consts::PI * 2.0;
        let avatars: Vec<String> = {
            let state = self.state.lock().unwrap();
            state
                .iter()
                .filter_map(|(key, _)| if key.starts_with("avatar.") {
                    Some(String::from(key.split_once('.').unwrap().1))
                } else {None})
                .collect()
        };
        let count = avatars.len();
        let angles = (0..count).map(|i| i as f64 * full_circle / (count as f64));
        let updates = StateUpdate{ changed_keys: Some(Struct {fields: BTreeMap::from_iter(avatars
            .iter()
            .zip(angles)
            .map(|(avatar, angle)| {
                let converted_angle = 0.5 * (-angle - full_circle / 4.0);
                (
                    format!("user-origin.{avatar}"),
                    [
                        ("position", vec![radius * angle.cos(), 0.0, radius * angle.sin()]),
                        ("rotation", vec![0.0, converted_angle.sin(), 0.0, converted_angle.cos()]),
                    ],
                )
            })
            .map(|(key, value)| {
                (
                    key,
                    Value{ kind: Some(Kind::StructValue(Struct{ fields: BTreeMap::from_iter(
                        value.into_iter().map(|(key, value)| {(
                            String::from(key),
                            Value{kind: Some(
                                Kind::ListValue(ListValue {
                                    values: value
                                        .iter()
                                        .map(|number| {
                                            Value {kind: Some(Kind::NumberValue(*number))}
                                        })
                                        .collect()
                                })
                            )}
                        )})
                    )}))},
                )
            })
        )})};
        self.state.lock().unwrap().send(updates).unwrap();
        CommandReply { result: None }
    }

    fn arguments(&self) -> Option<Struct> {
        None
    }
}