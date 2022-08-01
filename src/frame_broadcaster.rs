use crate::broadcaster::{Broadcaster, Mergeable, ReceiverVec};
use crate::frame::FrameData;
use std::sync::{Arc, Mutex};

pub struct FrameBroadcaster {
    receivers: ReceiverVec<FrameData>,
    current: Mutex<FrameData>,
}

impl FrameBroadcaster {}

impl FrameBroadcaster {
    pub fn new(base_frame: FrameData) -> Self {
        Self {
            receivers: Arc::new(Mutex::new(Vec::new())),
            current: Mutex::new(base_frame),
        }
    }
}

impl Broadcaster for FrameBroadcaster {
    type Content = FrameData;

    fn get_receivers(&self) -> ReceiverVec<FrameData> {
        Arc::clone(&self.receivers)
    }

    fn get_current(&self) -> Self::Content {
        self.current.lock().unwrap().clone()
    }

    fn update_current(&mut self, other: &Self::Content) {
        self.current.lock().unwrap().merge(other);
    }
}
