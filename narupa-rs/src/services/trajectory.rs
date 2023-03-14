use crate::broadcaster::{BroadcastReceiver, Broadcaster};
use narupa_proto::frame::FrameData;
use crate::frame_broadcaster::FrameBroadcaster;
use narupa_proto::trajectory::{
    trajectory_service_server::TrajectoryService, GetFrameRequest, GetFrameResponse,
};
use futures::Stream;
use log::debug;
use std::sync::{Arc, Mutex};
use std::{pin::Pin, time::Duration};
use tokio::sync::mpsc;
use tokio_stream::{wrappers::ReceiverStream, StreamExt};
use tonic::{Response, Status};

pub use narupa_proto::trajectory::trajectory_service_server::TrajectoryServiceServer;

type ResponseStream = Pin<Box<dyn Stream<Item = Result<GetFrameResponse, Status>> + Send + Sync>>;

struct FrameResponseIterator {
    frame_source: Arc<Mutex<BroadcastReceiver<FrameData>>>,
    frame_index: u32,
}

impl Iterator for FrameResponseIterator {
    type Item = GetFrameResponse;

    fn next(&mut self) -> Option<Self::Item> {
        let frame = self.frame_source.lock().unwrap().recv();
        self.frame_index = self.frame_index.wrapping_add(1);
        Some(GetFrameResponse {
            frame,
            frame_index: self.frame_index.wrapping_sub(1),
        })
    }
}

pub struct Trajectory {
    frame_source: Arc<Mutex<FrameBroadcaster>>,
    cancel_rx: Arc<Mutex<tokio::sync::oneshot::Receiver<()>>>,
}

impl Trajectory {
    pub fn new(frame_source: Arc<Mutex<FrameBroadcaster>>, cancel_rx: tokio::sync::oneshot::Receiver<()>) -> Self {
        Self { frame_source, cancel_rx: Arc::new(Mutex::new(cancel_rx)) }
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
            frame_index: 0,
        };
        let mut stream =
            Box::pin(tokio_stream::iter(responses).throttle(Duration::from_millis(interval)));
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
                        Err(_item) => {
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
