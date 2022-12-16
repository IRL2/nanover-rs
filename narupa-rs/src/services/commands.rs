use std::collections::HashMap;

use crate::playback::PlaybackOrder;
use crate::proto::protocol::command::command_server;
use crate::proto::protocol::command::{
    CommandMessage, CommandReply, GetCommandsReply, GetCommandsRequest,
};
use prost::alloc::vec::Vec;
use tokio::sync::mpsc::Sender;

pub use crate::proto::protocol::command::command_server::CommandServer;

pub trait Command: Send + Sync {
    fn run(&self, input: CommandMessage) -> CommandReply;
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
            .keys()
            .map(|name| CommandMessage {name: name.clone(), arguments: None})
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
