use std::io;

/// When a PDB field is ill-formatted, this enum tells what field has the issue.
#[derive(Debug)]
pub enum FieldError {
    Serial,
    ResidueIdentifier,
    Position,
}

/// What when wrong when reading an ATOM line in a PDB file?
#[derive(Debug)]
pub enum FormatError {
    FieldFormat(FieldError),
    LineTooShort,
    Unexpected,
    InconsistentKey,
    UnexpectedFieldNumber,
    MissingField(String),
    MissingLoopKeys,
}

#[derive(Debug)]
pub enum ReadError {
    FormatError(FormatError, usize),
    IOError(io::Error),
}