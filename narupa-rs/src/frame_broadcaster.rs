use crate::broadcaster::{Broadcaster, BroadcasterSignal, Mergeable, ReceiverVec};
use crate::frame::FrameData;
use std::sync::mpsc::Sender;
use std::sync::{Arc, Mutex};

pub struct FrameBroadcaster {
    receivers: ReceiverVec<FrameData>,
    current: Mutex<FrameData>,
    signal_tx: Option<Sender<BroadcasterSignal>>,
}

impl FrameBroadcaster {}

impl FrameBroadcaster {
    pub fn new(base_frame: FrameData, signal_tx: Option<Sender<BroadcasterSignal>>) -> Self {
        Self {
            receivers: Arc::new(Mutex::new(Vec::new())),
            current: Mutex::new(base_frame),
            signal_tx,
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

    fn get_signal_tx(&self) -> Option<Sender<BroadcasterSignal>> {
        self.signal_tx.as_ref().cloned()
    }
}
