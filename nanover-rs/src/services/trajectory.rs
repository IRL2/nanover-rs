use crate::broadcaster::{BroadcastReceiver, Broadcaster};
use crate::frame_broadcaster::FrameBroadcaster;
use futures::Stream;
use log::{debug, trace};
use nanover_proto::trajectory::{
    trajectory_service_server::TrajectoryService, GetFrameRequest, GetFrameResponse,
};
use std::sync::{Arc, Mutex};
use std::{pin::Pin, time::Duration};
use tokio::sync::mpsc;
use tokio_stream::{wrappers::ReceiverStream, StreamExt};
use tonic::{Response, Status};

pub use nanover_proto::trajectory::trajectory_service_server::TrajectoryServiceServer;

type ResponseStream = Pin<Box<dyn Stream<Item = Result<GetFrameResponse, Status>> + Send + Sync>>;

struct FrameResponseIterator {
    frame_source: Arc<Mutex<BroadcastReceiver<GetFrameResponse>>>,
}

impl Iterator for FrameResponseIterator {
    // This is a perpetual iterator: it always yields something.
    // When there is no response to send (the simulation did not send a new
    // frame yet), then we return Some(none) so we do not stop the iteration.
    // Stopping the itaretion is the responsability of the caller.
    type Item = Option<GetFrameResponse>;

    fn next(&mut self) -> Option<Self::Item> {
        let frame = self.frame_source.lock().unwrap().recv();
        Some(frame)
    }
}

pub struct Trajectory {
    frame_source: Arc<Mutex<FrameBroadcaster>>,
    cancel_rx: Arc<Mutex<tokio::sync::oneshot::Receiver<()>>>,
}

impl Trajectory {
    pub fn new(
        frame_source: Arc<Mutex<FrameBroadcaster>>,
        cancel_rx: tokio::sync::oneshot::Receiver<()>,
    ) -> Self {
        Self {
            frame_source,
            cancel_rx: Arc::new(Mutex::new(cancel_rx)),
        }
    }
}

#[tonic::async_trait]
impl TrajectoryService for Trajectory {
    type SubscribeLatestFramesStream = ResponseStream;

    async fn subscribe_latest_frames(
        &self,
        request: tonic::Request<GetFrameRequest>,
    ) -> Result<tonic::Response<Self::SubscribeLatestFramesStream>, tonic::Status> {
        debug!("New client subscribed to the frames!");
        let interval = (request.into_inner().frame_interval * 1000.0) as u64;
        let receiver = self.frame_source.lock().unwrap().get_rx();
        let responses = FrameResponseIterator {
            frame_source: receiver,
        };
        let mut stream = Box::pin(
            tokio_stream::iter(responses)
                // if frame_source did not have a frame to return, it sent
                // Some(None) instead. The filter_map ommits these useless
                // values.
                .filter_map(|item| item)
                .throttle(Duration::from_millis(interval)),
        );
        let (tx, rx) = mpsc::channel(128);
        let cancel_rx = Arc::clone(&self.cancel_rx);
        tokio::spawn(async move {
            while let Some(item) = stream.next().await {
                {
                    let mut cancel_lock = cancel_rx.lock().unwrap();
                    match cancel_lock.try_recv() {
                        Ok(_) | Err(tokio::sync::oneshot::error::TryRecvError::Closed) => {
                            cancel_lock.close();
                            break;
                        }
                        Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {}
                    }
                }
                if item.frame.is_some() {
                    match tx.send(Result::<_, Status>::Ok(item)).await {
                        Ok(_) => {
                            // item (server response) was queued to be send to client
                        }
                        Err(item) => {
                            trace!("Error in stream: {item}");
                            // output_stream was build from rx and both are dropped
                            break;
                        }
                    }
                }
            }
            debug!("Client disconnected from the frames.");
        });
        let output_stream = ReceiverStream::new(rx);
        Ok(Response::new(
            Box::pin(output_stream) as Self::SubscribeLatestFramesStream
        ))
    }
}
