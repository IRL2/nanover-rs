use std::sync::{Arc, Mutex};
use std::fmt::Debug;

pub type ReceiverVec<T> = Arc<Mutex<Vec<Arc<Mutex<BroadcastReceiver<T>>>>>>;

pub struct BroadcastReceiver<T> {
    content: Option<T>,
}

pub trait Broadcaster {
    type Content: Mergeable + Clone + Debug;

    fn get_rx(&mut self) -> Arc<Mutex<BroadcastReceiver<Self::Content>>> {
        let current = self.get_current();
        let receiver = Arc::new(Mutex::new(BroadcastReceiver::new(current)));
        let clone = Arc::clone(&receiver);
        self.get_receivers().lock().unwrap().push(receiver);
        clone
    }

    fn send(&mut self, item: Self::Content) -> Result<(), ()> {
        self.update_current(&item);
        let receivers = self.get_receivers();
        let receivers_locked = receivers.lock().unwrap();
        for receiver in &*receivers_locked {
            let cloned_receiver = Arc::clone(receiver);
            let mut cloned_locked = cloned_receiver.lock().unwrap();
            let mut content = &mut cloned_locked.content;
            match &mut content {
                None => *content = Some(item.clone()),
                Some(c) => c.merge(&item),
            }
        }
        Ok(())
    }

    fn get_receivers(&self) -> ReceiverVec<Self::Content>;
    fn get_current(&self) -> Self::Content;
    fn update_current(&mut self, other: &Self::Content);
}

impl<T: Debug> BroadcastReceiver<T> {
    pub fn new(current: T) -> Self {
        Self {content: Some(current)}
    }

    pub fn recv(&mut self) -> Option<T> {
        self.content.take()
    }
}

pub trait Mergeable {
    fn merge(&mut self, other: &Self);
}