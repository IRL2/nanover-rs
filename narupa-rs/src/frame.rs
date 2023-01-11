use narupa_proto::frame::FrameData;
use crate::broadcaster::Mergeable;

impl Mergeable for FrameData {
    fn merge(&mut self, other: &FrameData) {
        let cloned_other = other.clone();
        self.values.extend(cloned_other.values);
        self.arrays.extend(cloned_other.arrays);
    }
}