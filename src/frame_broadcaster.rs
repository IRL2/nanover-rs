use std::vec::Vec;
use std::sync::Arc;
use crate::frame::FrameData;
use crate::broadcaster::{Broadcaster, BroadcastReceiver};

pub struct FrameBroadcaster {
    receivers: Vec<Arc<BroadcastReceiver<FrameData>>>
}

impl FrameBroadcaster {}

impl Broadcaster<FrameData> for FrameBroadcaster {
    fn get_receivers(&mut self) -> &mut Vec<Arc<BroadcastReceiver<FrameData>>> {
        &mut self.receivers
    }
}