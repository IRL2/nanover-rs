use prost_types::{value::Kind, ListValue, Value};

pub trait UnPack<T> {
    fn unpack(self) -> Option<T>;
}

/// Build a Value from the object
///
/// ```
/// use pack_prost::ToProstValue;
/// use prost_types::{Value, value::Kind};
/// let number: f64 = 42.0;
/// let value = number.into_prost_value();
/// assert_eq!(value, Value{ kind: Some(Kind::NumberValue(42.0)) });
/// ```
///
/// It is safe to pack numbers smaller than f64 as they can be safely
/// cast to f64.
///
/// ```
/// use pack_prost::ToProstValue;
/// use prost_types::{Value, value::Kind};
/// let number: f32 = 42.0;
/// let value = number.into_prost_value();
/// assert_eq!(value, Value{ kind: Some(Kind::NumberValue(42.0)) });
/// ```
///
/// Text can also be packed.
///
/// ```
/// use pack_prost::ToProstValue;
/// use prost_types::{Value, value::Kind};
/// let text: String = String::from("Hello!");
/// let value = text.into_prost_value();
/// assert_eq!(
///     value,
///     Value{ kind: Some(Kind::StringValue("Hello!".to_string())) }
/// );
/// ```
///
/// ```
/// use pack_prost::ToProstValue;
/// use prost_types::{Value, value::Kind};
/// let text: &str = "Hello!";
/// let value = text.into_prost_value();
/// assert_eq!(
///     value,
///     Value{ kind: Some(Kind::StringValue("Hello!".to_string())) }
/// );
/// ```
///
/// Vectors can be packed as well if the type they contain can be packed.
///
/// ```
/// use pack_prost::ToProstValue;
/// use prost_types::{Value, ListValue, value::Kind};
/// let vector: Vec<bool> = vec![true, true, false];
/// let value = vector.into_prost_value();
/// assert_eq!(
///     value,
///     Value{ kind: Some(Kind::ListValue(ListValue { values: vec![
///         Value{ kind: Some(Kind::BoolValue(true)) },
///         Value{ kind: Some(Kind::BoolValue(true)) },
///         Value{ kind: Some(Kind::BoolValue(false)) },
///     ] })) }
/// );
/// ```
pub trait ToProstValue {
    fn into_prost_value(self) -> Value;
    fn to_prost_value(&self) -> Value;
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

impl<'a, T> UnPack<Vec<T>> for &'a Value
where
    &'a Value: UnPack<T>,
{
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
        vector_value
            .values
            .iter()
            .map(|item| item.unpack())
            .collect()
    }
}

impl<T> UnPack<Vec<T>> for Value
where
    Value: UnPack<T>,
{
    fn unpack(self) -> Option<Vec<T>> {
        let Some(Kind::ListValue(vector_value)) = self.kind else {
            return None;
        };
        vector_value
            .values
            .into_iter()
            .map(|item| item.unpack())
            .collect()
    }
}

impl ToProstValue for f64 {
    fn into_prost_value(self) -> Value {
        Value {
            kind: Some(Kind::NumberValue(self)),
        }
    }

    fn to_prost_value(&self) -> Value {
        Value {
            kind: Some(Kind::NumberValue(*self)),
        }
    }
}

impl ToProstValue for f32 {
    fn into_prost_value(self) -> Value {
        Value {
            kind: Some(Kind::NumberValue(self as f64)),
        }
    }

    fn to_prost_value(&self) -> Value {
        Value {
            kind: Some(Kind::NumberValue(*self as f64)),
        }
    }
}

impl ToProstValue for String {
    fn into_prost_value(self) -> Value {
        Value {
            kind: Some(Kind::StringValue(self)),
        }
    }

    fn to_prost_value(&self) -> Value {
        Value {
            kind: Some(Kind::StringValue(self.clone())),
        }
    }
}

impl ToProstValue for &String {
    fn into_prost_value(self) -> Value {
        Value {
            kind: Some(Kind::StringValue(self.clone())),
        }
    }

    fn to_prost_value(&self) -> Value {
        Value {
            kind: Some(Kind::StringValue(self.to_string())),
        }
    }
}

impl ToProstValue for &str {
    fn into_prost_value(self) -> Value {
        Value {
            kind: Some(Kind::StringValue(self.to_owned())),
        }
    }

    fn to_prost_value(&self) -> Value {
        Value {
            kind: Some(Kind::StringValue(self.to_string())),
        }
    }
}

impl ToProstValue for bool {
    fn into_prost_value(self) -> Value {
        Value {
            kind: Some(Kind::BoolValue(self)),
        }
    }

    fn to_prost_value(&self) -> Value {
        Value {
            kind: Some(Kind::BoolValue(*self)),
        }
    }
}

impl<T> ToProstValue for Vec<T>
where
    T: ToProstValue,
{
    fn into_prost_value(self) -> Value {
        let values = self
            .into_iter()
            .map(|item| item.into_prost_value())
            .collect();
        Value {
            kind: Some(Kind::ListValue(ListValue { values })),
        }
    }

    fn to_prost_value(&self) -> Value {
        let values = self.iter().map(|item| item.to_prost_value()).collect();
        Value {
            kind: Some(Kind::ListValue(ListValue { values })),
        }
    }
}
