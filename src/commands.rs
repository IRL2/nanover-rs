use crate::proto::protocol::command::command_server::Command;
use crate::proto::protocol::command::{GetCommandsRequest, GetCommandsReply, CommandMessage, CommandReply};
use prost::alloc::vec::Vec;

pub struct CommandService {}

#[tonic::async_trait]
impl Command for CommandService {
    async fn get_commands(
        &self,
        _request: tonic::Request<GetCommandsRequest>,
    ) -> Result<tonic::Response<GetCommandsReply>, tonic::Status> {
        let command_list: Vec<CommandMessage> = Vec::new();
        let reply = GetCommandsReply {commands: command_list};
        Ok(tonic::Response::new(reply))
    }
    
    async fn run_command(
        &self,
        _request: tonic::Request<CommandMessage>,
    ) -> Result<tonic::Response<CommandReply>, tonic::Status> {
        Err(tonic::Status::invalid_argument("Not implemented yet"))
    }
}