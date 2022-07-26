use std::sync::{Arc, Mutex};
use std::collections::{btree_map, BTreeMap};
use std::time::{Instant, Duration};
use prost_types::{Value, Struct, value::Kind};

use crate::proto::protocol::state::StateUpdate;
use crate::broadcaster::{Broadcaster, Mergeable, ReceiverVec};

pub struct StateLock {
    token: String,
    timeout: Option<Instant>,
}

pub struct StateBroadcaster {
    state: BTreeMap<String, Value>,
    locks: BTreeMap<String, StateLock>,
    receivers: ReceiverVec<StateUpdate>,
}

impl StateBroadcaster {
    pub fn new() -> Self {
        let state = BTreeMap::new();
        let locks = BTreeMap::new();
        let receivers = Arc::new(Mutex::new(Vec::new()));
        Self {state, receivers, locks}
    }

    pub fn iter(&self) -> btree_map::Iter<String, Value> {
        self.state.iter()
    }

    pub fn atomic_lock_updates(
            &mut self,
            token: String,
            requested_updates: BTreeMap<String, Option<Duration>>,
    ) -> Result<(), ()> {
        let now = Instant::now();

        let requested_removals: Vec<&String> = requested_updates.iter()
            .filter_map(|kv| match kv.1 {
                None => Some(kv.0),
                Some(_) => None,
            })
            .collect();

        for key in requested_removals {
            match self.locks.get(key) {
                Some(lock) if can_write(lock, &now, &token) => {
                    self.locks.remove(key);
                },
                _ => {},
            }
        }

        let mut updates: BTreeMap<String, StateLock> = BTreeMap::new();
        for (key, duration) in requested_updates {
            match self.locks.get(&key) {
                Some(lock) if !can_write(lock, &now, &token) => {
                    // There is a lock and we do not own it.
                    // We cancel the whole update.
                    return Err(());
                },
                _ => {
                    let timeout = match duration {
                        None => None,
                        Some(d) => Some(now + d),
                    };
                    updates.insert(key, StateLock {token: token.clone(), timeout})
                }
            };
        }
        self.locks.append(&mut updates);

        Ok(())
    }
    
    pub fn can_write_key(&self, key: &String, now: &Instant, token: &String) -> bool {
        match self.locks.get(key) {
            None => true,
            Some(lock) => can_write(lock, now, token),
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

impl Default for StateBroadcaster {
    fn default() -> Self {
        Self::new()
    }
}

impl Broadcaster for StateBroadcaster {
    type Content = StateUpdate;

    fn get_receivers(&self) -> ReceiverVec<Self::Content> {
        Arc::clone(&self.receivers)
    }

    fn get_current(&self) -> Self::Content {
        let structure = Struct{fields: self.state.clone()};
        StateUpdate{changed_keys: Some(structure)}
    }

    fn update_current(&mut self, other: &Self::Content) {
        match &other.changed_keys {
            None => {},
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
}

impl Mergeable for StateUpdate {
    fn merge(&mut self, other: &Self) {
        match (&mut self.changed_keys, &other.changed_keys) {
            (_, None) => {},
            (Some(ref mut current), Some(changed)) => {
                current.fields.extend(changed.clone().fields);
            },
            (None, Some(changed)) => {
                self.changed_keys = Some(changed.clone());
            }
        }
    }
}