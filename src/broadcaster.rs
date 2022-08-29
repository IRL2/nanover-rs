use std::fmt::Debug;
use std::sync::{Arc, Mutex};

pub type ReceiverVec<T> = Arc<Mutex<Vec<Arc<Mutex<BroadcastReceiver<T>>>>>>;

/// Broadcast data to several consumers at various rates.
///
/// Data is streamed as a set of differences with the previous
/// state. The broadcaster accumutates these differences for each
/// consumer, and for a full state that is sent as base for new
/// consumers. To be able to do this, the data needs to implement
/// the `Mergeable` trait.
///
/// The broadcaster can have multiple producers and multiple
/// consumers. Producers send updates using the `send` method
/// while new consumers are created using the `get_rx` one.
///
/// Creating a new consumer with `get_rx` returns a `BroadcarstReveiver`.
/// Its `recv` method pops the accumulated updates for that consumer.
///
/// Implementors need to implement `get_receivers`, `get_current`, and
/// `update_current`. They also need to specify the type of the data with
/// the `Content` attribute.
pub trait Broadcaster {
    type Content: Mergeable + Clone + Debug;

    /// Produce a new consumer.
    ///
    /// The consumer is a `BroadcastReceiver` that holds its
    /// accumulated update.
    fn get_rx(&mut self) -> Arc<Mutex<BroadcastReceiver<Self::Content>>> {
        let current = self.get_current();
        let receiver = Arc::new(Mutex::new(BroadcastReceiver::new(current)));
        let clone = Arc::clone(&receiver);
        self.get_receivers().lock().unwrap().push(receiver);
        clone
    }

    /// Send an update to all the receivers.
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

    /// List the receivers.
    fn get_receivers(&self) -> ReceiverVec<Self::Content>;
    /// Get a copy of the full accumulated state.
    fn get_current(&self) -> Self::Content;
    /// Update the full accumulated state.
    fn update_current(&mut self, other: &Self::Content);
}

/// A consumer of a broadcast.
///
/// The receiver holds the accumulated update since the last time `recv`
/// was called.
pub struct BroadcastReceiver<T> {
    content: Option<T>,
}

impl<T: Debug> BroadcastReceiver<T> {
    pub fn new(current: T) -> Self {
        Self {
            content: Some(current),
        }
    }

    /// Pop the accumulated update.
    ///
    /// Returns `None` if no update have been sent since the last call.
    pub fn recv(&mut self) -> Option<T> {
        self.content.take()
    }
}

/// Mergeable data for a broadcaster.
pub trait Mergeable {
    /// Combine another instance with the current instance, in place.
    fn merge(&mut self, other: &Self);
}
