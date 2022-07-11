use crate::proto::protocol::state::{
    StateUpdate,
    SubscribeStateUpdatesRequest,
    UpdateStateRequest, UpdateStateResponse,
    UpdateLocksRequest, UpdateLocksResponse,
};
use crate::proto::protocol::state::state_server::State;
use crate::state_broadcaster::StateBroadcaster;
use crate::broadcaster::{Broadcaster, BroadcastReceiver};
use futures::Stream;
use tonic::{Response, Status};
use tokio::sync::mpsc;
use tokio_stream::{wrappers::ReceiverStream, StreamExt};
use std::{pin::Pin, time::Duration, sync::{Arc, Mutex}};

pub use crate::proto::protocol::state::state_server::StateServer;


type ResponseStream = Pin<Box<dyn Stream<Item = Result<StateUpdate, tonic::Status>> + Send + Sync>>;

struct StateUpdateIterator {
    //frame_source: Arc<Mutex<FrameData>>
    update_source: Arc<Mutex<BroadcastReceiver<StateUpdate>>>
}

impl Iterator for StateUpdateIterator {
    type Item = Option<StateUpdate>;

    fn next(&mut self) -> Option<Self::Item> {
        let update = self.update_source.lock().unwrap().recv();
        Some(update)
    }
}

pub struct StateService {
    shared_state: Arc<Mutex<StateBroadcaster>>,
}

impl StateService {
    pub fn new(shared_state: Arc<Mutex<StateBroadcaster>>) -> Self {
        Self {shared_state}
    }
}

#[tonic::async_trait]
impl State for StateService {
    type SubscribeStateUpdatesStream = ResponseStream;

    async fn subscribe_state_updates(
        &self,
        request: tonic::Request<SubscribeStateUpdatesRequest>,
    ) -> Result<tonic::Response<Self::SubscribeStateUpdatesStream>, tonic::Status> {
        let update_source = {
            self.shared_state.lock().unwrap().get_rx()
        };
        let update_iterator = StateUpdateIterator {update_source};
        let interval = (request.into_inner().update_interval * 1000.0) as u64;
        let mut stream = Box::pin(tokio_stream::iter(update_iterator)
            .throttle(Duration::from_millis(interval)));

        // spawn and channel are required if you want handle "disconnect" functionality
        // the `out_stream` will not be polled after client disconnect
        let (tx, rx) = mpsc::channel(128);
        tokio::spawn(async move {
            while let Some(item) = stream.next().await {
                if let Some(update) = item {
                    match tx.send(Result::<_, Status>::Ok(update)).await {
                        Ok(_) => {
                            // item (server response) was queued to be send to client
                        }
                        Err(_) => {
                            // output_stream was build from rx and both are dropped
                            break;
                        }
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
        request: tonic::Request<UpdateStateRequest>,
    ) -> Result<tonic::Response<UpdateStateResponse>, tonic::Status> {
        if let Some(update) = request.into_inner().update {
            self.shared_state.lock().unwrap().send(update).unwrap();
        }
        Ok(Response::new(UpdateStateResponse {success: true}))
    }

    async fn update_locks(
        &self,
        _request: tonic::Request<UpdateLocksRequest>,
    ) -> Result<tonic::Response<UpdateLocksResponse>, tonic::Status> {
        Ok(Response::new(UpdateLocksResponse {success: true}))
    }
}