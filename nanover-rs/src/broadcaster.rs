use nanover_proto::Mergeable;
use std::fmt::Debug;
use std::sync::{mpsc::Sender, Arc, Mutex, Weak};
use std::time::Instant;
use thiserror::Error;

pub type ReceiverVec<T> = Arc<Mutex<Vec<Weak<Mutex<BroadcastReceiver<T>>>>>>;

#[derive(Debug, Error)]
#[error("Broadcaster cannot send message.")]
pub struct BroadcastSendError {}

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
/// When a broadcaster executes some actions, it sends signals over a
/// channel. These signals are of type `BroadcasterSignal`. The broadcaster
/// gives access to the channel's sender with the `get_signal_tx` method
/// that returns an option. Signals are only sent if there is a channel.
///
/// Implementors need to implement `get_receivers`, `get_current`,
/// `update_current`, and `get_signal_tx`. They also need to specify the type
/// of the data with the `Content` attribute.
pub trait Broadcaster {
    type Content: Mergeable + Clone + Debug;

    /// Produce a new consumer.
    ///
    /// The consumer is a `BroadcastReceiver` that holds its
    /// accumulated update.
    ///
    /// Send the `BroacasterSignal::NewReceiver` signal.
    fn get_rx(&mut self) -> Arc<Mutex<BroadcastReceiver<Self::Content>>> {
        self.send_broadaster_signal(BroadcasterSignal::NewReceiver(Instant::now()))
            .unwrap();
        let current = self.get_current();
        let receiver = Arc::new(Mutex::new(BroadcastReceiver::new(current)));
        let clone = Arc::downgrade(&receiver);
        self.get_receivers().lock().unwrap().push(clone);
        receiver
    }

    /// Send an update to all the receivers.
    ///
    /// Sends the `BroadcasterSignal::Send` signal every time. Sends also the
    /// `BroadcasterSignal::RemoveReceiver` for each receiver that gets removed.
    fn send(&mut self, item: Self::Content) -> Result<(), BroadcastSendError> {
        self.send_broadaster_signal(BroadcasterSignal::Send(Instant::now()))?;
        let mut to_remove: Vec<usize> = Vec::new();
        self.update_current(&item);
        let receivers = self.get_receivers();
        let mut receivers_locked = receivers.lock().unwrap();
        for (index, receiver) in receivers_locked.iter().enumerate() {
            match receiver.upgrade() {
                Some(cloned_receiver) => {
                    let mut cloned_locked = cloned_receiver.lock().unwrap();
                    let mut content = &mut cloned_locked.content;
                    match &mut content {
                        None => *content = Some(item.clone()),
                        Some(c) => c.merge(&item),
                    }
                }
                None => {
                    to_remove.push(index);
                }
            };
        }
        for index in to_remove.into_iter().rev() {
            self.send_broadaster_signal(BroadcasterSignal::RemoveReceiver(Instant::now()))?;
            receivers_locked.remove(index);
        }
        Ok(())
    }

    fn send_broadaster_signal(&self, signal: BroadcasterSignal) -> Result<(), BroadcastSendError> {
        if let Some(tx) = self.get_signal_tx() {
            tx.send(signal).map_err(|_| BroadcastSendError {})?;
        }
        Ok(())
    }

    /// List the receivers.
    fn get_receivers(&self) -> ReceiverVec<Self::Content>;
    /// Get a copy of the full accumulated state.
    fn get_current(&self) -> Self::Content;
    /// Update the full accumulated state.
    fn update_current(&mut self, other: &Self::Content);

    fn get_signal_tx(&self) -> Option<Sender<BroadcasterSignal>>;
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

/// Signal sent for some actions of a broadcaster.
pub enum BroadcasterSignal {
    /// Sent when an update is sent to the broadcaster.
    Send(Instant),
    /// Sent when a receiver is created.
    NewReceiver(Instant),
    /// Sent when a receiver is removed from the broadcaster.
    RemoveReceiver(Instant),
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;

    // We need a very simple type that implements Mergeable.
    // We cannot implement Mergeable on an existing type because
    // of the orphan rule, so we create a new type here.
    #[derive(Debug, Clone, Copy)]
    struct DummyData {
        data: usize,
    }

    impl DummyData {
        pub fn new(data: usize) -> Self {
            Self { data }
        }
    }

    impl Mergeable for DummyData {
        fn merge(&mut self, other: &Self) {
            self.data += other.data;
        }
    }

    struct DummyBroadcaster {
        receivers: ReceiverVec<DummyData>,
        current: DummyData,
        signal_tx: Option<Sender<BroadcasterSignal>>,
    }

    impl DummyBroadcaster {
        pub fn new(signal_tx: Option<Sender<BroadcasterSignal>>) -> Self {
            DummyBroadcaster {
                receivers: Arc::new(Mutex::new(Vec::new())),
                current: DummyData::new(0),
                signal_tx,
            }
        }
    }

    impl Broadcaster for DummyBroadcaster {
        type Content = DummyData;

        fn get_current(&self) -> Self::Content {
            self.current
        }

        fn get_receivers(&self) -> ReceiverVec<Self::Content> {
            Arc::clone(&self.receivers)
        }

        fn update_current(&mut self, other: &Self::Content) {
            self.current.merge(other);
        }

        fn get_signal_tx(&self) -> Option<Sender<BroadcasterSignal>> {
            self.signal_tx.clone()
        }
    }

    fn assert_num_receivers<T>(broadcaster: &T, expected: usize)
    where
        T: Broadcaster,
    {
        assert_eq!(broadcaster.get_receivers().lock().unwrap().len(), expected);
    }

    #[test]
    fn test_adding_receivers() {
        let mut broadcaster = DummyBroadcaster::new(None);

        // Initially there is no receiver in the broadcaster.
        assert_num_receivers(&broadcaster, 0);

        // We add some receivers.
        let rx0 = broadcaster.get_rx();
        let rx1 = broadcaster.get_rx();
        let rx2 = broadcaster.get_rx();
        assert_num_receivers(&broadcaster, 3);

        // We remove the receivers. This is only reflected after a send.
        // We first remove one receiver, then the 2 others. This makes sure
        // that the behaviour is correct regardless of the number of receivers
        // dropped in between two sends.
        drop(rx0);
        assert_num_receivers(&broadcaster, 3);
        broadcaster.send(DummyData::new(1)).unwrap();
        assert_num_receivers(&broadcaster, 2);
        drop(rx1);
        drop(rx2);
        broadcaster.send(DummyData::new(2)).unwrap();
        assert_num_receivers(&broadcaster, 0);
    }

    /// When sending an update, the `Send` broadcaster signal is sent.
    #[test]
    fn test_sending_send_signal() {
        let (signal_tx, signal_rx) = mpsc::channel();
        let mut broadcaster = DummyBroadcaster::new(Some(signal_tx));

        // Before we do anything, the channel is empty.
        assert!(signal_rx.try_recv().is_err());

        broadcaster.send(DummyData::new(0)).unwrap();
        broadcaster.send(DummyData::new(1)).unwrap();
        assert!(matches! {signal_rx.try_recv().unwrap(), BroadcasterSignal::Send(_)});
        assert!(matches! {signal_rx.try_recv().unwrap(), BroadcasterSignal::Send(_)});
        assert!(signal_rx.try_recv().is_err());
    }

    /// Adding and removing receivers send the appropriate broadcaster signals.
    #[test]
    fn test_sending_receivers_signal() {
        let (signal_tx, signal_rx) = mpsc::channel();
        let mut broadcaster = DummyBroadcaster::new(Some(signal_tx));

        // Before we do anything, the channel is empty.
        assert!(signal_rx.try_recv().is_err());

        let rx = broadcaster.get_rx();
        drop(rx);
        broadcaster.send(DummyData::new(0)).unwrap();

        assert!(matches! {signal_rx.try_recv().unwrap(), BroadcasterSignal::NewReceiver(_)});
        assert!(matches! {signal_rx.try_recv().unwrap(), BroadcasterSignal::Send(_)});
        assert!(matches! {signal_rx.try_recv().unwrap(), BroadcasterSignal::RemoveReceiver(_)});
    }
}
