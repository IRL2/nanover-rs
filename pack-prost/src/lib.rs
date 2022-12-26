use prost_types::{Value, ListValue, value::Kind};

pub trait UnPack<T> {
    fn unpack(self) -> Option<T>;
}

/// Build a Value from the object
/// 
/// ```
/// use pack_prost::Pack;
/// use prost_types::{Value, value::Kind};
/// let number: f64 = 42.0;
/// let value = number.pack();
/// assert_eq!(value, Value{ kind: Some(Kind::NumberValue(42.0)) });
/// ```
/// 
/// It is safe to pack numbers smaller than f64 as they can be safely
/// cast to f64.
/// 
/// ```
/// use pack_prost::Pack;
/// use prost_types::{Value, value::Kind};
/// let number: f32 = 42.0;
/// let value = number.pack();
/// assert_eq!(value, Value{ kind: Some(Kind::NumberValue(42.0)) });
/// ```
/// 
/// Text can also be packed.
/// 
/// ```
/// use pack_prost::Pack;
/// use prost_types::{Value, value::Kind};
/// let text: String = String::from("Hello!");
/// let value = text.pack();
/// assert_eq!(
///     value,
///     Value{ kind: Some(Kind::StringValue("Hello!".to_string())) }
/// );
/// ```
/// 
/// ```
/// use pack_prost::Pack;
/// use prost_types::{Value, value::Kind};
/// let text: &str = "Hello!";
/// let value = text.pack();
/// assert_eq!(
///     value,
///     Value{ kind: Some(Kind::StringValue("Hello!".to_string())) }
/// );
/// ```
/// 
/// Vectors can be packed as well if the type they contain can be packed.
/// 
/// ```
/// use pack_prost::Pack;
/// use prost_types::{Value, ListValue, value::Kind};
/// let vector: Vec<bool> = vec![true, true, false];
/// let value = vector.pack();
/// assert_eq!(
///     value,
///     Value{ kind: Some(Kind::ListValue(ListValue { values: vec![
///         Value{ kind: Some(Kind::BoolValue(true)) },
///         Value{ kind: Some(Kind::BoolValue(true)) },
///         Value{ kind: Some(Kind::BoolValue(false)) },
///     ] })) }
/// );
/// ```
pub trait Pack {
    fn pack(self) -> Value;
}

impl UnPack<f64> for &Value {
    /// Get a number out of a reference to a prost Value
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, value::Kind};
    /// let value = Value{ kind: Some(Kind::NumberValue(42.0)) };
    /// let unpacked: f64 = (&value).unpack().unwrap();
    /// assert_eq!(unpacked, 42.0);
    /// ```
    fn unpack(self) -> Option<f64> {
        let Some(Kind::NumberValue(number)) = self.kind else {
            return None;
        };
        Some(number)
    }
}

impl UnPack<f64> for Value {
    /// Get a number out of a prost Value
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, value::Kind};
    /// let value = Value{ kind: Some(Kind::NumberValue(42.0)) };
    /// let unpacked: f64 = value.unpack().unwrap();
    /// assert_eq!(unpacked, 42.0);
    /// ```
    fn unpack(self) -> Option<f64> {
        let Some(Kind::NumberValue(number)) = self.kind else {
            return None;
        };
        Some(number)
    }
}

impl UnPack<bool> for &Value {
    /// Get a boolean out of a reference to a prost Value
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, value::Kind};
    /// let value = Value{ kind: Some(Kind::BoolValue(true)) };
    /// let unpacked: bool = (&value).unpack().unwrap();
    /// assert_eq!(unpacked, true);
    /// ```
    fn unpack(self) -> Option<bool> {
        let Some(Kind::BoolValue(value)) = self.kind else {
            return None;
        };
        Some(value)
    }
}

impl UnPack<bool> for Value {
    /// Get a boolean out of a prost Value
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, value::Kind};
    /// let value = Value{ kind: Some(Kind::BoolValue(true)) };
    /// let unpacked: bool = value.unpack().unwrap();
    /// assert_eq!(unpacked, true);
    /// ```
    fn unpack(self) -> Option<bool> {
        let Some(Kind::BoolValue(value)) = self.kind else {
            return None;
        };
        Some(value)
    }
}

impl UnPack<String> for &Value {
    /// Get a string out of a reference to a prost Value
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, value::Kind};
    /// let value = Value{ kind: Some(Kind::StringValue(String::from("Hello"))) };
    /// let unpacked: String = (&value).unpack().unwrap();
    /// assert_eq!(unpacked, String::from("Hello"));
    /// ```
    fn unpack(self) -> Option<String> {
        let Some(Kind::StringValue(ref content)) = self.kind else {
            return None;
        };
        Some(content.clone())
    }
}

impl UnPack<String> for Value {
    /// Get a string out of a prost Value
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, value::Kind};
    /// let value = Value{ kind: Some(Kind::StringValue(String::from("Hello"))) };
    /// let unpacked: String = value.unpack().unwrap();
    /// assert_eq!(unpacked, String::from("Hello"));
    /// ```
    fn unpack(self) -> Option<String> {
        let Some(Kind::StringValue(content)) = self.kind else {
            return None;
        };
        Some(content)
    }
}

impl<'a, T> UnPack<Vec<T>> for &'a Value where &'a Value: UnPack<T> {
    /// Get a vector of homogeneous values out of a reference to a prost value.
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, ListValue, value::Kind};
    /// let value = Value { kind: Some(Kind::ListValue(ListValue{ values: vec![
    ///     Value { kind: Some(Kind::NumberValue(23.0)) },
    ///     Value { kind: Some(Kind::NumberValue(42.0)) },
    ///     Value { kind: Some(Kind::NumberValue(12.0)) },
    /// ]}))};
    /// let unpacked: Vec<f64> = (&value).unpack().unwrap();
    /// assert_eq!(unpacked, vec![23.0, 42.0, 12.0]);
    /// ```
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, ListValue, value::Kind};
    /// let value = Value { kind: Some(Kind::ListValue(ListValue{ values: vec![
    ///     Value { kind: Some(Kind::StringValue(String::from("one"))) },
    ///     Value { kind: Some(Kind::StringValue(String::from("two"))) },
    ///     Value { kind: Some(Kind::StringValue(String::from("three"))) },
    ///     Value { kind: Some(Kind::StringValue(String::from("four"))) },
    /// ]}))};
    /// let unpacked: Vec<String> = (&value).unpack().unwrap();
    /// assert_eq!(
    ///     unpacked,
    ///     vec![
    ///         String::from("one"),
    ///         String::from("two"),
    ///         String::from("three"),
    ///         String::from("four"),
    ///     ]
    /// );
    /// ```
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, ListValue, value::Kind};
    /// let value = Value { kind: Some(Kind::ListValue(ListValue{ values: vec![
    ///     Value { kind: Some(Kind::ListValue(ListValue{ values: vec![
    ///         Value { kind: Some(Kind::StringValue(String::from("one"))) },
    ///         Value { kind: Some(Kind::StringValue(String::from("two"))) },
    ///         Value { kind: Some(Kind::StringValue(String::from("three"))) },
    ///         Value { kind: Some(Kind::StringValue(String::from("four"))) },
    ///     ]}))},
    ///     Value { kind: Some(Kind::ListValue(ListValue{ values: vec![
    ///         Value { kind: Some(Kind::StringValue(String::from("five"))) },
    ///         Value { kind: Some(Kind::StringValue(String::from("six"))) },
    ///     ]}))},
    /// ]}))};
    /// let unpacked: Vec<Vec<String>> = (&value).unpack().unwrap();
    /// assert_eq!(
    ///     unpacked,
    ///     vec![
    ///         vec![
    ///             String::from("one"),
    ///             String::from("two"),
    ///             String::from("three"),
    ///             String::from("four"),
    ///         ],
    ///         vec![
    ///             String::from("five"),
    ///             String::from("six"),
    ///         ],
    ///     ]
    /// );
    /// ```
    /// 
    /// ```
    /// use pack_prost::UnPack;
    /// use prost_types::{Value, ListValue, value::Kind};
    /// let value = Value { kind: Some(Kind::ListValue(ListValue{ values: vec![
    ///     Value { kind: Some(Kind::NumberValue(23.0)) },
    ///     Value { kind: Some(Kind::StringValue(String::from("two"))) },
    ///     Value { kind: Some(Kind::BoolValue(false)) },
    /// ]}))};
    /// let unpacked: Option<Vec<f64>> = (&value).unpack();
    /// assert_eq!(unpacked, None);
    /// ```
    fn unpack(self) -> Option<Vec<T>> {
        let Some(Kind::ListValue(ref vector_value)) = self.kind else {
            return None;
        };
        vector_value.values.iter().map(|item| item.unpack()).collect()
    }
}

impl<T> UnPack<Vec<T>> for Value where Value: UnPack<T> {
    fn unpack(self) -> Option<Vec<T>> {
        let Some(Kind::ListValue(vector_value)) = self.kind else {
            return None;
        };
        vector_value.values.into_iter().map(|item| item.unpack()).collect()
    }
}

impl Pack for f64 {
    fn pack(self) -> Value {
        Value{ kind: Some(Kind::NumberValue(self)) }
    }
}

impl Pack for f32 {
    fn pack(self) -> Value {
        Value{ kind: Some(Kind::NumberValue(self as f64)) }
    }
}

impl Pack for String {
    fn pack(self) -> Value {
        Value { kind: Some(Kind::StringValue(self)) }
    }
}

impl Pack for &String {
    fn pack(self) -> Value {
        Value { kind: Some(Kind::StringValue(self.clone())) }
    }
}

impl Pack for &str {
    fn pack(self) -> Value {
        Value { kind: Some(Kind::StringValue(self.to_owned())) }
    }
}

impl Pack for bool {
    fn pack(self) -> Value {
        Value { kind: Some(Kind::BoolValue(self)) }
    }
}

impl<T> Pack for Vec<T> where T: Pack {
    fn pack(self) -> Value {
        let values = self.into_iter().map(|item| item.pack()).collect();
        Value { kind: Some(Kind::ListValue(ListValue{ values })) }
    }
}