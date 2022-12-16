// Copied from https://stackoverflow.com/a/74550371
macro_rules! test_ressource {($fname:expr) => (
    concat!(env!("CARGO_MANIFEST_DIR"), $fname) // assumes Linux ('/')!
)}

pub(crate) use test_ressource;