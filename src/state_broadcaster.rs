use std::sync::{Arc, Mutex};
use std::collections::BTreeMap;
use prost_types::{Value, Struct};

use crate::proto::protocol::state::StateUpdate;
use crate::broadcaster::{Broadcaster, Mergeable, ReceiverVec};

pub struct StateBroadcaster {
    state: BTreeMap<String, Value>,
    receivers: ReceiverVec<StateUpdate>,
}

impl StateBroadcaster {
    pub fn new() -> Self {
        let state = BTreeMap::new();
        let receivers = Arc::new(Mutex::new(Vec::new()));
        Self {state, receivers}
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
                // TODO: prune keys when the value is nil
                // use self.state.retain
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