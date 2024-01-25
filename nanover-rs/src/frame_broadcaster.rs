use crate::broadcaster::{BroadcastSendError, Broadcaster, BroadcasterSignal, ReceiverVec};
use nanover_proto::frame::{FrameData, GetFrameResponse};
use nanover_proto::Mergeable;
use std::sync::mpsc::Sender;
use std::sync::{Arc, Mutex};

pub struct FrameBroadcaster {
    receivers: ReceiverVec<GetFrameResponse>,
    current: Mutex<GetFrameResponse>,
    signal_tx: Option<Sender<BroadcasterSignal>>,
    next_frame_index: u32,
}

impl FrameBroadcaster {
    pub fn send_frame(&mut self, frame: FrameData) -> Result<(), BroadcastSendError> {
        let response = GetFrameResponse {
            frame_index: self.next_frame_index,
            frame: Some(frame),
        };
        self.next_frame_index += 1;
        self.send(response)
    }

    pub fn send_reset_frame(&mut self, frame: FrameData) -> Result<(), BroadcastSendError> {
        let response = GetFrameResponse {
            frame_index: 0,
            frame: Some(frame),
        };
        self.next_frame_index = 1;
        self.send(response)
    }
}

impl FrameBroadcaster {
    pub fn new(base_frame: GetFrameResponse, signal_tx: Option<Sender<BroadcasterSignal>>) -> Self {
        let next_frame_index = base_frame.frame_index + 1;
        Self {
            receivers: Arc::new(Mutex::new(Vec::new())),
            current: Mutex::new(base_frame),
            signal_tx,
            next_frame_index,
        }
    }
}

impl Broadcaster for FrameBroadcaster {
    type Content = GetFrameResponse;

    fn get_receivers(&self) -> ReceiverVec<GetFrameResponse> {
        Arc::clone(&self.receivers)
    }

    fn get_current(&self) -> Self::Content {
        self.current.lock().unwrap().clone()
    }

    fn update_current(&mut self, other: &Self::Content) {
        self.current.lock().unwrap().merge(other);
    }

    fn get_signal_tx(&self) -> Option<Sender<BroadcasterSignal>> {
        self.signal_tx.as_ref().cloned()
    }
}
