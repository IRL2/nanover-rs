use std::sync::{Arc, Mutex};

pub struct BroadcastReceiver<T> {
    content: Mutex<Option<T>>,
}

pub trait Broadcaster<T> {
    fn get_rx(&mut self) -> Arc<BroadcastReceiver<T>> {
        let receiver = Arc::new(BroadcastReceiver::new());
        let clone = Arc::clone(&receiver);
        self.get_receivers().push(receiver);
        clone
    }

    fn send(&self, _item: T) -> Result<(), ()> {
        Ok(())
    }

    fn get_receivers(&mut self) -> &mut Vec<Arc<BroadcastReceiver<T>>>;
}

impl<T> BroadcastReceiver<T> {
    fn new() -> Self {
        Self {content: Mutex::new(None)}
    }

    fn recv(&mut self) -> Result<T, ()> {
        Err(())
    }
}