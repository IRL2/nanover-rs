use crate::broadcaster::{BroadcastReceiver, Broadcaster};
use narupa_proto::state::state_server::State;
use narupa_proto::state::{
    StateUpdate, SubscribeStateUpdatesRequest, UpdateLocksRequest, UpdateLocksResponse,
    UpdateStateRequest, UpdateStateResponse,
};
use crate::state_broadcaster::StateBroadcaster;
use futures::Stream;
use prost_types::value::Kind;
use std::collections::BTreeMap;
use std::{
    pin::Pin,
    sync::{Arc, Mutex},
    time::Duration,
};
use tokio::sync::mpsc;
use tokio_stream::{wrappers::ReceiverStream, StreamExt};
use tonic::{Response, Status};
use log::trace;

pub use narupa_proto::state::state_server::StateServer;

type ResponseStream = Pin<Box<dyn Stream<Item = Result<StateUpdate, tonic::Status>> + Send + Sync>>;

struct StateUpdateIterator {
    update_source: Arc<Mutex<BroadcastReceiver<StateUpdate>>>,
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
        Self { shared_state }
    }
}

#[tonic::async_trait]
impl State for StateService {
    type SubscribeStateUpdatesStream = ResponseStream;

    async fn subscribe_state_updates(
        &self,
        request: tonic::Request<SubscribeStateUpdatesRequest>,
    ) -> Result<tonic::Response<Self::SubscribeStateUpdatesStream>, tonic::Status> {
        let update_source = { self.shared_state.lock().unwrap().get_rx() };
        let update_iterator = StateUpdateIterator { update_source };
        let interval = (request.into_inner().update_interval * 1000.0) as u64;
        let mut stream =
            Box::pin(tokio_stream::iter(update_iterator).throttle(Duration::from_millis(interval)));

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
        let inner_request = request.into_inner();
        let token = inner_request.access_token;
        let can_update = match inner_request.update {
            None => {
                // There is nothing to update so it necesserally succeed.
                true
            }
            Some(update) => {
                let mut state = self.shared_state.lock().unwrap();
                let can_update = state.send_with_locks(update, &token).is_ok();
                can_update
            }
        };
        Ok(Response::new(UpdateStateResponse {
            success: can_update,
        }))
    }

    async fn update_locks(
        &self,
        request: tonic::Request<UpdateLocksRequest>,
    ) -> Result<tonic::Response<UpdateLocksResponse>, tonic::Status> {
        let inner_request = request.into_inner();
        let success_response = Ok(Response::new(UpdateLocksResponse { success: true }));

        trace!("update_locks {inner_request:?}");

        // This is an undocumented behavior in the python version:
        // the lock_keys attribute of the request can be None, but
        // this is not handled by the python version of te server.
        // Hre, we assume that this case is valid and results in a
        // no-op.
        if inner_request.lock_keys.is_none() {
            return success_response;
        };
        let lock_keys = inner_request.lock_keys.unwrap();

        let token = inner_request.access_token;
        let mut requested_updates: BTreeMap<String, Option<Duration>> = BTreeMap::new();
        for update in lock_keys.fields.into_iter() {
            let (key, value) = update;
            match value.kind {
                // Protobuf values can be empty without being null.
                // Here, we conflate a lack of value and a null value.
                None | Some(Kind::NullValue(_)) => requested_updates.insert(key.clone(), None),

                Some(Kind::NumberValue(seconds)) => {
                    requested_updates.insert(key.clone(), Some(Duration::from_secs_f64(seconds)))
                }
                Some(_) => {
                    return Err(tonic::Status::invalid_argument("Misformated lock update."));
                }
            };
        }

        let mut state = self.shared_state.lock().unwrap();
        match state.atomic_lock_updates(token, requested_updates) {
            Ok(()) => success_response,
            Err(()) => Ok(Response::new(UpdateLocksResponse { success: false })),
        }
    }
}
