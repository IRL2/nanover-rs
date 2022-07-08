use std::sync::{Arc, Mutex};

pub type ReceiverVec<T> = Arc<Mutex<Vec<Arc<Mutex<BroadcastReceiver<T>>>>>>;

pub struct BroadcastReceiver<T> {
    content: Option<T>,
}

pub trait Broadcaster {
    type Update;
    type Content: Mergeable<Self::Update> + Clone;

    fn get_rx(&mut self) -> Arc<Mutex<BroadcastReceiver<Self::Content>>> {
        let current = self.get_current();
        let receiver = Arc::new(Mutex::new(BroadcastReceiver::new(current.clone())));
        let clone = Arc::clone(&receiver);
        self.get_receivers().lock().unwrap().push(receiver);
        clone
    }

    fn send(&mut self, item: Self::Update) -> Result<(), ()> {
        self.update_current(&item);
        let receivers = self.get_receivers();
        let receivers_locked = receivers.lock().unwrap();
        for receiver in &*receivers_locked {
            let cloned_receiver = Arc::clone(receiver);
            let mut cloned_locked = cloned_receiver.lock().unwrap();
            let mut content = &mut cloned_locked.content;
            match &mut content {
                None => *content = Some(self.get_current().clone()),
                Some(c) => c.merge(&item),
            }
        }
        Ok(())
    }

    fn get_receivers(&self) -> ReceiverVec<Self::Content>;
    fn get_current(&self) -> Self::Content;
    fn update_current(&mut self, other: &Self::Update);
}

impl<T> BroadcastReceiver<T> {
    pub fn new(current: T) -> Self {
        Self {content: Some(current)}
    }

    pub fn recv(&mut self) -> Option<T> {
        self.content.take()
    }
}

pub trait Mergeable<T> {
    fn merge(&mut self, other: &T);
}