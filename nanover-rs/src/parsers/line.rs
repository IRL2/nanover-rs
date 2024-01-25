use crate::parsers::Position;

/// The parsed fields from an ATOM line in a PDB file.
#[derive(Debug, PartialEq)]
pub struct PDBLine {
    pub serial: isize,
    pub atom_name: String,
    pub alternate: char,
    pub residue_name: String,
    pub chain_identifier: char,
    pub residue_identifier: isize,
    pub insertion_code: char,
    pub position: Position,
    pub element_symbol: Option<usize>,
}
