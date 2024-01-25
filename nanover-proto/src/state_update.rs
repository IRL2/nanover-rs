pub use crate::protocol::state::StateUpdate;
use crate::Mergeable;

impl Mergeable for StateUpdate {
    fn merge(&mut self, other: &Self) {
        match (&mut self.changed_keys, &other.changed_keys) {
            (_, None) => {}
            (Some(ref mut current), Some(changed)) => {
                current.fields.extend(changed.clone().fields);
            }
            (None, Some(changed)) => {
                self.changed_keys = Some(changed.clone());
            }
        }
    }
}
