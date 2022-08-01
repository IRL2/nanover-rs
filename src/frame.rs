use crate::broadcaster::Mergeable;
use crate::proto::protocol::value_array::Values;
use crate::proto::protocol::{FloatArray, IndexArray, ValueArray};
use prost_types::{value::Kind, Value};
use std::collections::HashMap;

pub use crate::proto::protocol::trajectory::FrameData;

impl FrameData {
    pub fn empty() -> FrameData {
        let values: HashMap<::prost::alloc::string::String, Value> = HashMap::new();
        let arrays: HashMap<::prost::alloc::string::String, ValueArray> = HashMap::new();

        FrameData { values, arrays }
    }

    pub fn insert_number_value(&mut self, key: &str, value: f64) -> Result<(), ()> {
        let serialed_value = Value {
            kind: Some(Kind::NumberValue(value)),
        };
        match self.values.insert(key.to_string(), serialed_value) {
            None => Ok(()),
            Some(_) => Err(()),
        }
    }

    pub fn insert_float_array(&mut self, key: &str, value: Vec<f32>) -> Result<(), ()> {
        let float_array = Values::FloatValues(FloatArray { values: value });
        let value_array = ValueArray {
            values: Some(float_array),
        };
        match self.arrays.insert(key.to_string(), value_array) {
            None => Ok(()),
            Some(_) => Err(()),
        }
    }

    pub fn insert_index_array(&mut self, key: &str, value: Vec<u32>) -> Result<(), ()> {
        let index_array = Values::IndexValues(IndexArray { values: value });
        let value_array = ValueArray {
            values: Some(index_array),
        };
        match self.arrays.insert(key.to_string(), value_array) {
            None => Ok(()),
            Some(_) => Err(()),
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
