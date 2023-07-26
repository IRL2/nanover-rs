use crate::protocol::value_array::Values;
use crate::protocol::{FloatArray, IndexArray, StringArray, ValueArray};
use crate::Mergeable;
use prost_types::Value;
use std::collections::HashMap;
use pack_prost::ToProstValue;

pub use crate::protocol::trajectory::FrameData;
pub use crate::protocol::trajectory::GetFrameResponse;

#[derive(Debug)]
pub struct ExistingDataError {}

impl FrameData {
    pub fn empty() -> FrameData {
        let values: HashMap<::prost::alloc::string::String, Value> = HashMap::new();
        let arrays: HashMap<::prost::alloc::string::String, ValueArray> = HashMap::new();

        FrameData { values, arrays }
    }

    pub fn insert_number_value(&mut self, key: &str, value: f64) -> Result<(), ExistingDataError> {
        let serialed_value = value.into_prost_value();
        match self.values.insert(key.to_string(), serialed_value) {
            None => Ok(()),
            Some(_) => Err(ExistingDataError {}),
        }
    }

    pub fn insert_float_array(
        &mut self,
        key: &str,
        value: Vec<f32>,
    ) -> Result<(), ExistingDataError> {
        let float_array = Values::FloatValues(FloatArray { values: value });
        let value_array = ValueArray {
            values: Some(float_array),
        };
        match self.arrays.insert(key.to_string(), value_array) {
            None => Ok(()),
            Some(_) => Err(ExistingDataError {}),
        }
    }

    pub fn insert_index_array(
        &mut self,
        key: &str,
        value: Vec<u32>,
    ) -> Result<(), ExistingDataError> {
        let index_array = Values::IndexValues(IndexArray { values: value });
        let value_array = ValueArray {
            values: Some(index_array),
        };
        match self.arrays.insert(key.to_string(), value_array) {
            None => Ok(()),
            Some(_) => Err(ExistingDataError {}),
        }
    }

    pub fn insert_string_array(
        &mut self,
        key: &str,
        value: Vec<String>,
    ) -> Result<(), ExistingDataError> {
        let string_array = Values::StringValues(StringArray { values: value });
        let value_array = ValueArray {
            values: Some(string_array),
        };
        match self.arrays.insert(key.to_string(), value_array) {
            None => Ok(()),
            Some(_) => Err(ExistingDataError {}),
        }
    }
}

impl Mergeable for FrameData {
    fn merge(&mut self, other: &FrameData) {
        let cloned_other = other.clone();
        self.values.extend(cloned_other.values);
        self.arrays.extend(cloned_other.arrays);
    }
}

impl Mergeable for GetFrameResponse {
    fn merge(&mut self, other: &Self) {
        // We keep the index of the initial frame. However, if the new frame has
        // an index of 0, then it is a resetting frame: the existing frame is
        // erased and the frame index is set to 0 to convey this information.
        let other_frame_index = other.frame_index;
        if other_frame_index == 0 {
            self.frame_index = 0;
            self.frame = other.frame.clone();
        } else {
            match (&mut self.frame, &other.frame) {
                (_, None) => {}
                (None, Some(other_frame)) => {
                    let _ = self.frame.insert(other_frame.clone());
                }
                (Some(self_frame), Some(other_frame)) => {
                    self_frame.merge(other_frame);
                }
            }
        }
    }
}