use crate::proto::protocol::state::{
    StateUpdate,
    SubscribeStateUpdatesRequest,
    UpdateStateRequest, UpdateStateResponse,
    UpdateLocksRequest, UpdateLocksResponse,
};
use crate::proto::protocol::state::state_server::State;
use futures::Stream;
use tonic::{Response, Status};
use tokio::sync::mpsc;
use tokio_stream::{wrappers::ReceiverStream, StreamExt};
use std::{pin::Pin, time::Duration};


type ResponseStream = Pin<Box<dyn Stream<Item = Result<StateUpdate, tonic::Status>> + Send + Sync>>;

pub struct StateService {}

#[tonic::async_trait]
impl State for StateService {
    type SubscribeStateUpdatesStream = ResponseStream;

    async fn subscribe_state_updates(
        &self,
        _request: tonic::Request<SubscribeStateUpdatesRequest>,
    ) -> Result<tonic::Response<Self::SubscribeStateUpdatesStream>, tonic::Status> {
        let repeat = std::iter::repeat(StateUpdate {changed_keys: None});
        let mut stream = Box::pin(tokio_stream::iter(repeat).throttle(Duration::from_millis(200)));

        // spawn and channel are required if you want handle "disconnect" functionality
        // the `out_stream` will not be polled after client disconnect
        let (tx, rx) = mpsc::channel(128);
        tokio::spawn(async move {
            while let Some(item) = stream.next().await {
                match tx.send(Result::<_, Status>::Ok(item)).await {
                    Ok(_) => {
                        // item (server response) was queued to be send to client
                    }
                    Err(_item) => {
                        // output_stream was build from rx and both are dropped
                        break;
                    }
                }
            }
        });
        let output_stream = ReceiverStream::new(rx);
        Ok(Response::new(
            Box::pin(output_stream) as Self::SubscribeStateUpdatesStream
        ))
    }

    async fn update_state(
        &self,
        _request: tonic::Request<UpdateStateRequest>,
    ) -> Result<tonic::Response<UpdateStateResponse>, tonic::Status> {
        Ok(Response::new(UpdateStateResponse {success: true}))
    }

    async fn update_locks(
        &self,
        _request: tonic::Request<UpdateLocksRequest>,
    ) -> Result<tonic::Response<UpdateLocksResponse>, tonic::Status> {
        Ok(Response::new(UpdateLocksResponse {success: true}))
    }
}