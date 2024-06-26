use std::io;

/// When a PDB field is ill-formatted, this enum tells what field has the issue.
#[derive(Debug)]
pub enum FieldError {
    Serial,
    ResidueIdentifier,
    Position,
    Conect,
}

/// What when wrong when reading an ATOM line in a PDB file?
#[derive(Debug)]
pub enum FormatError {
    FieldFormat(FieldError),
    LineTooShort,
    Unexpected,
    InconsistentKey,
    // (expected, found)
    UnexpectedFieldNumber(usize, usize),
    MissingField(String),
    MissingLoopKeys,
    ContentOutOfData,
    MissformatedData,
}

#[derive(Debug)]
pub enum ReadError {
    FormatError(FormatError, usize),
    IOError(io::Error),
}
