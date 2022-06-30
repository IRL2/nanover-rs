use crate::proto::protocol::trajectory::{
    trajectory_service_server::TrajectoryService,
    GetFrameRequest,
    GetFrameResponse,
};
use crate::frame::FrameData;
use crate::frame_broadcaster::FrameBroadcaster;
use crate::broadcaster::{BroadcastReceiver, Broadcaster};
use futures::Stream;
use std::{error::Error, io::ErrorKind, net::ToSocketAddrs, pin::Pin, time::Duration};
use tonic::{transport::Server, Request, Response, Status, Streaming};
use tokio_stream::{wrappers::ReceiverStream, StreamExt};
use tokio::sync::mpsc;
use std::sync::{Arc, Mutex};

type ResponseStream = Pin<Box<dyn Stream<Item = Result<GetFrameResponse, Status>> + Send + Sync>>;

struct FrameResponseIterator {
    frame_source: Arc<Mutex<FrameData>>
}

impl Iterator for FrameResponseIterator {
    type Item = GetFrameResponse;

    fn next(&mut self) -> Option<Self::Item> {
        let frame = Some(self.frame_source.lock().unwrap().clone());
        Some(GetFrameResponse {frame: frame, frame_index: 0})
    }
}

pub struct Trajectory {
//    frame_source: Arc<Mutex<FrameData>>,
    frame_source: Arc<Mutex<FrameData>>,
}

impl Trajectory {
    pub fn new(frame_source: Arc<Mutex<FrameData>>) -> Self {
        Self {frame_source: frame_source}
    }
}

#[tonic::async_trait]
impl TrajectoryService for Trajectory {
    type SubscribeLatestFramesStream = ResponseStream;

    async fn subscribe_latest_frames(
        &self,
        _request: tonic::Request<GetFrameRequest>,
    ) -> Result<tonic::Response<Self::SubscribeLatestFramesStream>, tonic::Status> {
        println!("Hello there!");
        //let mut frame = FrameData::empty();
        //frame.insert_number_value("message.number", 2.3).unwrap();
        //let frame = self.frame_source.recv().await.unwrap().unwrap();
        //let repeat = std::iter::repeat(GetFrameResponse {frame_index: 1, frame: Some(frame)});
        //let mut stream = Box::pin(tokio_stream::iter(repeat).throttle(Duration::from_millis(200)));
        let receiver = Arc::clone(&self.frame_source);
        let responses = FrameResponseIterator {frame_source: receiver};
        let mut stream = Box::pin(tokio_stream::iter(responses).throttle(Duration::from_millis(200)));
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
            println!("\tclient disconnected");
        });
        let output_stream = ReceiverStream::new(rx);
        Ok(Response::new(
            Box::pin(output_stream) as Self::SubscribeLatestFramesStream
        ))
    }
}
