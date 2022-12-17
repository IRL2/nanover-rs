use std::collections::{HashMap, BTreeMap};
use std::sync::{Arc, Mutex};

use crate::playback::PlaybackOrder;
use crate::proto::protocol::command::command_server;
use crate::proto::protocol::command::{
    CommandMessage, CommandReply, GetCommandsReply, GetCommandsRequest,
};
use crate::broadcaster::Broadcaster;
use crate::proto::protocol::state::StateUpdate;
use crate::state_broadcaster::StateBroadcaster;
use prost::alloc::vec::Vec;
use prost_types::{Struct, Value, ListValue};
use prost_types::value::Kind;
use tokio::sync::mpsc::Sender;

pub use crate::proto::protocol::command::command_server::CommandServer;

pub trait Command: Send + Sync {
    fn run(&self, input: CommandMessage) -> CommandReply;
    fn arguments(&self) -> Option<Struct>;
}

pub struct PlaybackCommand {
    channel: Sender<PlaybackOrder>,
    order: PlaybackOrder
}

impl PlaybackCommand {
    pub fn new(channel: Sender<PlaybackOrder>, order: PlaybackOrder) -> Self {
        Self {channel, order}
    }
}

impl Command for PlaybackCommand {
    fn run(&self, _input: CommandMessage) -> CommandReply {
        self.channel.try_send(self.order).unwrap();
        CommandReply { result: None }
    }

    fn arguments(&self) -> Option<Struct> {
        None
    }
}

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

pub struct CommandService {
    commands: HashMap<String, Box<dyn Command>>,
}

impl CommandService {
    pub fn new(commands: HashMap<String, Box<dyn Command>>) -> Self {
        CommandService { commands }
    }
}

#[tonic::async_trait]
impl command_server::Command for CommandService {
    async fn get_commands(
        &self,
        _request: tonic::Request<GetCommandsRequest>,
    ) -> Result<tonic::Response<GetCommandsReply>, tonic::Status> {
        // TODO: Provide the command arguments
        let command_list: Vec<CommandMessage> = self
            .commands
            .iter()
            .map(|(name, command)| CommandMessage { name: name.clone(), arguments: command.arguments() })
            .collect();
        let reply = GetCommandsReply {
            commands: command_list,
        };
        Ok(tonic::Response::new(reply))
    }

    async fn run_command(
        &self,
        request: tonic::Request<CommandMessage>,
    ) -> Result<tonic::Response<CommandReply>, tonic::Status> {
        let request = request.into_inner();
        let name = request.name.as_str();
        println!("Received command named {name}");
        let command = self.commands.get(name);
        let Some(command) = command else {
            return Err(tonic::Status::invalid_argument("Not implemented yet"));
        };
        let response = command.run(request);
        Ok(tonic::Response::new(response))
    }
}
