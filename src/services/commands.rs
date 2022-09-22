use std::sync::{Arc, Mutex, mpsc::Sender};
use crate::playback::PlaybackOrder;
use crate::proto::protocol::command::command_server::Command;
use crate::proto::protocol::command::{
    CommandMessage, CommandReply, GetCommandsReply, GetCommandsRequest,
};
use prost::alloc::vec::Vec;

pub use crate::proto::protocol::command::command_server::CommandServer;

pub struct CommandService {
    channel: Arc<Mutex<Sender<PlaybackOrder>>>,
    //channel: tokio::sync::mpsc::Sender<PlaybackOrder>,
}

impl CommandService {
    pub fn new(channel: Sender<PlaybackOrder>) -> Self {
        CommandService{ channel: Arc::new(Mutex::new(channel)) }
    }
    /*
    pub fn new(channel: tokio::sync::mpsc::Sender<PlaybackOrder>) -> Self {
        CommandService { channel }
    }
    */
}

#[tonic::async_trait]
impl Command for CommandService {
    async fn get_commands(
        &self,
        _request: tonic::Request<GetCommandsRequest>,
    ) -> Result<tonic::Response<GetCommandsReply>, tonic::Status> {
        let command_list: Vec<CommandMessage> = Vec::new();
        let reply = GetCommandsReply {
            commands: command_list,
        };
        Ok(tonic::Response::new(reply))
    }

    async fn run_command(
        &self,
        request: tonic::Request<CommandMessage>,
    ) -> Result<tonic::Response<CommandReply>, tonic::Status> {
        let name = request.into_inner().name;
        let name = name.as_str();
        println!("Received command named {name}");
        match name {
            "playback/play" => {
                self
                    .channel
                    .lock()
                    .unwrap()
                    .send(PlaybackOrder::Play)
                    .unwrap();
                Ok(tonic::Response::new(CommandReply { result: None }))
            },
            "playback/pause" => {
                self
                    .channel
                    .lock()
                    .unwrap()
                    .send(PlaybackOrder::Pause)
                    .unwrap();
                Ok(tonic::Response::new(CommandReply { result: None }))
            },
            "playback/reset" => {
                self
                    .channel
                    .lock()
                    .unwrap()
                    .send(PlaybackOrder::Reset)
                    .unwrap();
                Ok(tonic::Response::new(CommandReply { result: None }))
            },
            _ => {
                Err(tonic::Status::invalid_argument("Not implemented yet"))
            }
        }
    }
}
