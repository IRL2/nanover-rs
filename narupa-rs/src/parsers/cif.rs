use std::{
    collections::HashMap,
    io::BufRead,
    str::FromStr,
};
use crate::parsers::chemistry::lookup_element_symbol;
use crate::parsers::errors::*;
use crate::parsers::molecular_system::MolecularSystem;
use crate::parsers::line::PDBLine;

#[derive(Debug)]
enum PDBXContext {
    /// We are out of any particular context
    /// (e.g. beginning of the file or after a loop)
    Idle,
    /// We just read loop_ or a a loop key.
    /// Next line can be a loop key.
    /// The argument is the name of the loop, it is optional as
    /// the name is unknown until we read the first key.
    LoopKey(Option<String>),
    /// We read a loop record, we cannot accept anything but new records.
    /// The argument is the name of the loop.
    Loop(String),
    // We read a data block with the name in argument.
    Data(String),
    // We read a save block with the name in argument.
    // Save(String)
}

pub fn read_cif<F>(input: F) -> Result<MolecularSystem, ReadError>
where
    F: BufRead,
{
    let mut context = PDBXContext::Idle;
    let mut loop_keys: Option<Vec<String>> = None;
    let mut atoms: Vec<PDBLine> = Vec::new();
    for (lineno, line) in input.lines().enumerate() {
        let line = line.map_err(ReadError::IOError)?;
        let mut tokens = line.split_ascii_whitespace();
        context = match context {
            PDBXContext::Idle => {
                match tokens.next() {
                    // We ignore empty lines.
                    None => PDBXContext::Idle,
                    Some(first) if first.starts_with("data_") => {
                        let Some((_, name)) = first.split_once('_') else {
                            return Err(ReadError::FormatError(FormatError::MissformatedData, lineno));
                        };
                        PDBXContext::Data(String::from(name))
                    },
                    _ => return Err(ReadError::FormatError(FormatError::ContentOutOfData, lineno)),
                }
            }
            PDBXContext::Data(name) => {
                match tokens.next() {
                    // We ignore empty lines.
                    None => PDBXContext::Data(name),
                    // We only read the first data block. If we encounter another one,
                    // it means we are done.
                    Some(first) if first.starts_with("data_") => break,
                    Some(first) if first == "loop_" => PDBXContext::LoopKey(None),
                    // Data items have a first token starting with '_' and the rest of the
                    // tokens are the value. We do not read any at the moment so we ignore
                    // these lines.
                    Some(first) if first.starts_with('_') => PDBXContext::Data(name),
                    Some(first) if first.starts_with('#') => PDBXContext::Data(name),
                    _ => {
                        return Err(ReadError::FormatError(FormatError::Unexpected, lineno));
                    }
                }
            }
            PDBXContext::LoopKey(perhaps_name) => {
                match (tokens.next(), perhaps_name.clone()) {
                    (None, _) => PDBXContext::LoopKey(perhaps_name), // We ignore empty lines
                    (Some(first), _) if first.starts_with('#') => PDBXContext::Idle,
                    (Some(first), Some(name)) if !first.starts_with(&name) => {
                        if first.starts_with('_') {
                            return Err(ReadError::FormatError(
                                FormatError::InconsistentKey,
                                lineno,
                            ));
                        }
                        let loop_keys = match loop_keys {
                            None => {
                                return Err(ReadError::FormatError(
                                    FormatError::MissingLoopKeys,
                                    lineno,
                                ))
                            }
                            Some(ref keys) => keys,
                        };
                        // We only store what we need.
                        if name == "_atom_site" {
                            atoms.push(
                                parse_cif_atom_line(&line, loop_keys)
                                    .map_err(|e| ReadError::FormatError(e, lineno))?,
                            );
                        }
                        PDBXContext::Loop(name)
                    }
                    (Some(first), _) => {
                        let mut name_parts = first.split('.');
                        let base_key = name_parts
                            .next()
                            .ok_or(ReadError::FormatError(FormatError::InconsistentKey, lineno))?;
                        let sub_key = String::from(name_parts.as_str());
                        if let Some(ref mut keys) = loop_keys {
                            keys.push(sub_key);
                        } else {
                            loop_keys = Some(vec![sub_key]);
                        };
                        PDBXContext::LoopKey(Some(String::from(base_key)))
                    }
                }
            }
            PDBXContext::Loop(name) => {
                let loop_keys = &loop_keys
                    .clone()
                    .ok_or(ReadError::FormatError(FormatError::MissingLoopKeys, lineno))?;
                // We only store what we need.
                if name == "_atom_site" {
                    atoms.push(
                        parse_cif_atom_line(&line, loop_keys)
                            .map_err(|e| ReadError::FormatError(e, lineno))?,
                    );
                }
                PDBXContext::Loop(name)
            }
        };
    }

    Ok(MolecularSystem::from(atoms)
        .add_intra_residue_bonds()
        .add_inter_residue_bonds()
    )
}

fn parse_cif_atom_line(line: &str, loop_keys: &Vec<String>) -> Result<PDBLine, FormatError> {
    let fields: Vec<String> = line.split_ascii_whitespace().map(String::from).collect();
    if fields.len() != loop_keys.len() {
        return Err(FormatError::UnexpectedFieldNumber);
    };
    let line: HashMap<String, String> =
        HashMap::from_iter(loop_keys.iter().zip(fields).map(|(k, v)| (k.clone(), v)));

    let x: f64 = extract_number(&line, "Cartn_x", FieldError::Position)?;
    let y: f64 = extract_number(&line, "Cartn_y", FieldError::Position)?;
    let z: f64 = extract_number(&line, "Cartn_z", FieldError::Position)?;
    let position = [x / 10.0, y / 10.0, z / 10.0]; // We use nanometers

    Ok(PDBLine {
        serial: extract_number(&line, "id", FieldError::Serial)?,
        atom_name: extract_string(&line, "auth_atom_id")?,
        alternate: extract_char_with_default(&line, "label_alt_id", ' '),
        residue_name: extract_string(&line, "auth_comp_id")?,
        chain_identifier: extract_char_with_default(&line, "auth_asym_id", ' '),
        residue_identifier: extract_number(&line, "auth_seq_id", FieldError::ResidueIdentifier)?,
        insertion_code: extract_char_with_default(&line, "pdbx_PDB_ins_code", ' '),
        position,
        element_symbol: lookup_element_symbol(
            line.get("type_symbol").unwrap_or(&String::from(" ")),
        ),
    })
}

fn extract_number<N>(
    line: &HashMap<String, String>,
    key: &str,
    error: FieldError,
) -> Result<N, FormatError>
where
    N: FromStr,
{
    line.get(key)
        .ok_or(FormatError::MissingField(String::from(key)))?
        .parse::<N>()
        .or(Err(FormatError::FieldFormat(error)))
}

fn extract_char_with_default(line: &HashMap<String, String>, key: &str, default: char) -> char {
    line.get(key)
        .unwrap_or(&String::from(""))
        .chars()
        .next()
        .unwrap_or(default)
}

fn extract_string(line: &HashMap<String, String>, key: &str) -> Result<String, FormatError> {
    Ok(line
        .get(key)
        .ok_or(FormatError::MissingField(String::from(key)))?
        .clone())
}