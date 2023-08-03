use log::trace;
use prost_types::{value::Kind, Struct, Value};
use std::collections::{btree_map, BTreeMap};
use std::sync::mpsc::Sender;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::broadcaster::{Broadcaster, BroadcasterSignal, ReceiverVec};
use narupa_proto::state_update::StateUpdate;

#[derive(Debug)]
pub struct StateLock {
    token: String,
    timeout: Option<Instant>,
}

pub struct StateBroadcaster {
    state: BTreeMap<String, Value>,
    locks: BTreeMap<String, StateLock>,
    receivers: ReceiverVec<StateUpdate>,
    signal_tx: Option<Sender<BroadcasterSignal>>,
}

/// Broadcast state updates to multiple consumers.
///
/// The broadcaster maintains the accumulated state, accumulated updates
/// each consumer, and the state locks.
///
/// # Example
///
/// ```
/// use std::collections::BTreeMap;
/// use prost_types::value::Kind;
/// use narupa_rs::state_broadcaster::StateBroadcaster;
/// use narupa_rs::broadcaster::Broadcaster;
/// use narupa_proto::Mergeable;
/// use narupa_proto::state::StateUpdate;
///
/// // Create a state broadcaster with an empty state
/// let mut broadcaster = StateBroadcaster::new(None);
///
/// // Create a consumer to receive state updates
/// let mut receiver_A = broadcaster.get_rx();
/// // At this point the receiver's state is the initial empty one
/// // from the brand new broadcaster.
/// assert!(receiver_A.lock().unwrap().recv() == Some(broadcaster.get_current()));
/// // Because we poped the latest update, the receiver is now empty.
/// assert!(receiver_A.lock().unwrap().recv() == None);
///
/// // Send an update.
/// // An update is directly a `StateUpdate` from the narupa protobuf
/// // protocol. Building one by hand is unpleasant, though these updates
/// // are meant to be built by the clients.
/// let mut raw_update = BTreeMap::new();
/// let key = "hello".to_string();
/// let value = prost_types::Value{kind: Some(Kind::StringValue("world".to_string()))};
/// raw_update.insert(key, value);
/// let update = StateUpdate {changed_keys: Some(prost_types::Struct {fields: raw_update})};
/// // We send a copy so we can reuse the update
/// broadcaster.send(update.clone());
///
/// // Now the receiver has received the update
/// assert!(receiver_A.lock().unwrap().recv() == Some(update.clone()));
///
/// // We can create other receivers.
/// // Each receiver maintains its own accumulated update.
/// // A new receiver starts with the full accumulated update.
/// let receiver_B = broadcaster.get_rx();
///
/// // Send another update.
/// let mut raw_update = BTreeMap::new();
/// let key = "foo".to_string();
/// let value = prost_types::Value{kind: Some(Kind::StringValue("bar".to_string()))};
/// raw_update.insert(key, value);
/// let update_2 = StateUpdate {changed_keys: Some(prost_types::Struct {fields: raw_update})};
/// // We send a copy so we can reuse the update
/// broadcaster.send(update_2.clone());
///
/// // The first receiver only has the last update in memory.
/// assert!(receiver_A.lock().unwrap().recv() == Some(update_2.clone()));
///
/// // The second receiver has a merge of both the first and the second
/// // update since we did not pop the first update yet.
/// let mut merged_updates = update.clone();
/// merged_updates.merge(&update_2);
/// assert!(receiver_B.lock().unwrap().recv() == Some(merged_updates));
/// ```
impl StateBroadcaster {
    pub fn new(signal_tx: Option<Sender<BroadcasterSignal>>) -> Self {
        let state = BTreeMap::new();
        let locks = BTreeMap::new();
        let receivers = Arc::new(Mutex::new(Vec::new()));
        Self {
            state,
            receivers,
            locks,
            signal_tx,
        }
    }

    /// Iterate over the keys and values of the accumulated state.
    pub fn iter(&self) -> btree_map::Iter<String, Value> {
        self.state.iter()
    }

    /// Update one or more state locks.
    pub fn atomic_lock_updates(
        &mut self,
        token: String,
        requested_updates: BTreeMap<String, Option<Duration>>,
    ) -> Result<(), ()> {
        let now = Instant::now();

        trace!("Atomic lock update {requested_updates:?}");

        let (requested_removals, requested_adds): (
            Vec<(String, Option<Duration>)>,
            Vec<(String, Option<Duration>)>,
        ) = requested_updates
            .into_iter()
            .partition(|(_, value)| value.is_none());

        for (key, _) in requested_removals {
            match self.locks.get(&key) {
                Some(lock) if can_write(lock, &now, &token) => {
                    trace!("Remove lock {key}");
                    self.locks.remove(&key);
                }
                _ => {}
            }
        }

        let mut updates: BTreeMap<String, StateLock> = BTreeMap::new();
        for (key, duration) in requested_adds {
            match self.locks.get(&key) {
                Some(lock) if !can_write(lock, &now, &token) => {
                    // There is a lock and we do not own it.
                    // We cancel the whole update.
                    return Err(());
                }
                _ => {
                    trace!("Adding lock {key}");
                    let timeout = duration.map(|d| now + d);
                    updates.insert(
                        key,
                        StateLock {
                            token: token.clone(),
                            timeout,
                        },
                    );
                }
            };
        }
        self.locks.append(&mut updates);
        trace!("Current locks {:?}", self.locks);

        Ok(())
    }

    pub fn can_write_key(&self, key: &String, now: &Instant, token: &String) -> bool {
        match self.locks.get(key) {
            None => true,
            Some(lock) => can_write(lock, now, token),
        }
    }

    pub fn send_with_locks(&mut self, item: StateUpdate, token: &String) -> Result<(), ()> {
        let now = Instant::now();
        let can_update = match &item.changed_keys {
            None => true,
            Some(changed_keys) => changed_keys
                .fields
                .keys()
                .all(|key| self.can_write_key(key, &now, token)),
        };
        if can_update {
            self.send(item).unwrap();
            Ok(())
        } else {
            Err(())
        }
    }
}

fn has_lock_expired(lock: &StateLock, now: &Instant) -> bool {
    match lock.timeout {
        None => false,
        Some(timeout) if &timeout < now => false,
        Some(_) => true,
    }
}

fn can_write(lock: &StateLock, now: &Instant, token: &String) -> bool {
    &lock.token == token || has_lock_expired(lock, now)
}

impl Broadcaster for StateBroadcaster {
    type Content = StateUpdate;

    fn get_receivers(&self) -> ReceiverVec<Self::Content> {
        Arc::clone(&self.receivers)
    }

    fn get_current(&self) -> Self::Content {
        let structure = Struct {
            fields: self.state.clone(),
        };
        StateUpdate {
            changed_keys: Some(structure),
        }
    }

    fn update_current(&mut self, other: &Self::Content) {
        match &other.changed_keys {
            None => {}
            Some(ref changes) => {
                self.state.extend(changes.fields.clone());
                self.state.retain(|_, value| {
                    // A Null value in an update marks the key-value pair to be
                    // deleted in the shared state. So a Null value in the shared
                    // state is a key-value pair to be deleted. Here, we treat the
                    // absence of value and an actual NullValue as being the same.
                    !matches!(value.kind, None | Some(Kind::NullValue(_)))
                });
            }
        }
    }

    fn get_signal_tx(&self) -> Option<Sender<BroadcasterSignal>> {
        self.signal_tx.as_ref().cloned()
    }
}
