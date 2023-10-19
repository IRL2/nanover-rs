use crate::parsers::chemistry::lookup_element_symbol;
use crate::parsers::errors::*;
use crate::parsers::line::PDBLine;
use crate::parsers::molecular_system::MolecularSystem;
use std::{collections::HashMap, io::BufRead, str::FromStr};

#[derive(Debug)]
enum PDBXContext {
    /// We are out of any particular context
    /// (e.g. beginning of the file or after a loop)
    Idle,
    /// We just read loop_ or a a loop key.
    /// Next line can be a loop key.
    /// The first argument is the name of the data block we are in.
    /// The second argument is the name of the loop, it is optional as
    /// the name is unknown until we read the first key.
    LoopKey(String, Option<String>),
    /// We read a loop record, we cannot accept anything but new records.
    /// The first argument is the name of the data block we are in.
    /// The second argument is the name of the loop.
    Loop(String, String),
    // We read a data block with the name in argument.
    Data(String),
    // We read a field that spans multiple lines.
    // The argument is the name of the data block we are in.
    // The second is false until we start reading the content of the field.
    Field(String, bool),
    // We read a save block with the name in argument.
    // Save(String)
}

#[derive(PartialEq, Eq, Hash)]
struct BondedAtom {
    pub chain: String,
    pub residue_name: String,
    pub residue_num: isize,
    pub atom_name: String,
}

impl From<&PDBLine> for BondedAtom {
    fn from(value: &PDBLine) -> Self {
        BondedAtom {
            chain: String::from(value.chain_identifier),
            residue_name: value.residue_name.clone(),
            residue_num: value.residue_identifier,
            atom_name: value.atom_name.clone(),
        }
    }
}

type PreBond = (BondedAtom, BondedAtom);

pub fn read_cif<F>(input: F) -> Result<MolecularSystem, ReadError>
where
    F: BufRead,
{
    let mut context = PDBXContext::Idle;
    let mut loop_keys: Option<Vec<String>> = None;
    let mut atoms: Vec<PDBLine> = Vec::new();
    let mut bonds: Vec<PreBond> = Vec::new();
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
                            return Err(ReadError::FormatError(
                                FormatError::MissformatedData,
                                lineno,
                            ));
                        };
                        PDBXContext::Data(String::from(name))
                    }
                    _ => {
                        return Err(ReadError::FormatError(
                            FormatError::ContentOutOfData,
                            lineno,
                        ))
                    }
                }
            }
            PDBXContext::Data(name) => {
                match tokens.next() {
                    // We ignore empty lines.
                    None => PDBXContext::Data(name),
                    // We only read the first data block. If we encounter another one,
                    // it means we are done.
                    Some(first) if first.starts_with("data_") => break,
                    Some("loop_") => PDBXContext::LoopKey(name, None),
                    // Data items have a first token starting with '_' and the rest of the
                    // tokens are the value. We do not read any at the moment so we ignore
                    // these lines. However, the content of the item can span multiple lines.
                    // In that case, we need to switch context so we can ignore the full value.
                    Some(first) if first.starts_with('_') => {
                        if tokens.next().is_some() {
                            PDBXContext::Data(name)
                        } else {
                            PDBXContext::Field(name, false)
                        }
                    }
                    Some(first) if first.starts_with('#') => PDBXContext::Data(name),
                    _ => {
                        return Err(ReadError::FormatError(FormatError::Unexpected, lineno));
                    }
                }
            }
            PDBXContext::LoopKey(data_name, perhaps_name) => {
                match (tokens.next(), perhaps_name.clone()) {
                    (None, _) => PDBXContext::LoopKey(data_name, perhaps_name), // We ignore empty lines
                    (Some(first), _) if first.starts_with('#') => {
                        loop_keys = None;
                        PDBXContext::Data(data_name)
                    }
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
                        match name.as_str() {
                            "_atom_site" => {
                                atoms.push(
                                    parse_cif_atom_line(&line, loop_keys)
                                        .map_err(|e| ReadError::FormatError(e, lineno))?,
                                );
                            }
                            "_struct_conn" => {
                                if let Some(bond) = parse_cif_bond_line(&line, loop_keys)
                                    .map_err(|e| ReadError::FormatError(e, lineno))?
                                {
                                    bonds.push(bond);
                                };
                            }
                            _ => (),
                        }
                        PDBXContext::Loop(data_name, name)
                    }
                    (Some(first), _) => {
                        let (base_key, sub_key) = first
                            .split_once('.')
                            .ok_or(ReadError::FormatError(FormatError::InconsistentKey, lineno))?;
                        let sub_key = String::from(sub_key);
                        if let Some(ref mut keys) = loop_keys {
                            keys.push(sub_key);
                        } else {
                            loop_keys = Some(vec![sub_key]);
                        };
                        PDBXContext::LoopKey(data_name, Some(String::from(base_key)))
                    }
                }
            }
            PDBXContext::Loop(data_name, name) => {
                if line.starts_with('#') {
                    loop_keys = None;
                    PDBXContext::Data(data_name)
                } else {
                    let loop_keys = &loop_keys
                        .clone()
                        .ok_or(ReadError::FormatError(FormatError::MissingLoopKeys, lineno))?;
                    // We only store what we need.
                    match name.as_str() {
                        "_atom_site" => {
                            atoms.push(
                                parse_cif_atom_line(&line, loop_keys)
                                    .map_err(|e| ReadError::FormatError(e, lineno))?,
                            );
                        }
                        "_struct_conn" => {
                            if let Some(bond) = parse_cif_bond_line(&line, loop_keys)
                                .map_err(|e| ReadError::FormatError(e, lineno))?
                            {
                                bonds.push(bond);
                            };
                        }
                        _ => (),
                    }
                    PDBXContext::Loop(data_name, name)
                }
            }
            PDBXContext::Field(data_name, content_started) => {
                if content_started && line.starts_with(';') {
                    PDBXContext::Data(data_name)
                } else if !content_started && !line.starts_with(';') {
                    let trimmed_line = line.trim();
                    if !trimmed_line.starts_with('\'') || !trimmed_line.ends_with('\'') {
                        return Err(ReadError::FormatError(FormatError::Unexpected, lineno));
                    }
                    PDBXContext::Data(data_name)
                } else {
                    PDBXContext::Field(data_name, true)
                }
            }
        };
    }

    let bonds = match_bonds_to_atoms(&bonds, &atoms);
    let mut molecular_system = MolecularSystem::from(atoms);
    molecular_system.bonds.extend(bonds);

    Ok(molecular_system
        .add_intra_residue_bonds()
        .add_inter_residue_bonds())
}

fn match_bonds_to_atoms(originals: &[PreBond], atoms: &[PDBLine]) -> Vec<(usize, usize, f32)> {
    let atom_to_index: HashMap<BondedAtom, usize> = HashMap::from_iter(
        atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| (atom.into(), index)),
    );
    originals
        .iter()
        .filter_map(|bond| {
            let from = atom_to_index.get(&bond.0);
            let to = atom_to_index.get(&bond.1);
            from.zip(to).map(|pre| (*pre.0, *pre.1, 1.0))
        })
        .collect()
}

fn parse_cif_bond_line(
    line: &str,
    loop_keys: &Vec<String>,
) -> Result<Option<PreBond>, FormatError> {
    let fields: Vec<String> = line.split_ascii_whitespace().map(String::from).collect();
    if fields.len() != loop_keys.len() {
        return Err(FormatError::UnexpectedFieldNumber(
            loop_keys.len(),
            fields.len(),
        ));
    };
    let line: HashMap<String, String> =
        HashMap::from_iter(loop_keys.iter().zip(fields).map(|(k, v)| (k.clone(), v)));

    let from = BondedAtom {
        chain: String::from(extract_char_with_default(&line, "ptnr1_label_asym_id", ' ')),
        residue_name: extract_string(&line, "ptnr1_label_comp_id")?,
        residue_num: extract_number(&line, "ptnr1_label_seq_id", FieldError::Conect)?,
        atom_name: extract_string(&line, "ptnr1_label_atom_id")?,
    };
    let to = BondedAtom {
        chain: String::from(extract_char_with_default(&line, "ptnr2_label_asym_id", ' ')),
        residue_name: extract_string(&line, "ptnr2_label_comp_id")?,
        residue_num: extract_number(&line, "ptnr2_label_seq_id", FieldError::Conect)?,
        atom_name: extract_string(&line, "ptnr2_label_atom_id")?,
    };

    let bond_type = extract_string(&line, "conn_type_id")?.to_uppercase();
    // The bond type can be "covale", "disulf","hydrog", or "metalc".
    // Only the two first types are covalent bond we want in our topology.
    // Cif documentation: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_conn.conn_type_id.html
    if bond_type == "COVALE" || bond_type == "DISULF" {
        Ok(Some((from, to)))
    } else {
        Ok(None)
    }
}

fn parse_cif_atom_line(line: &str, loop_keys: &Vec<String>) -> Result<PDBLine, FormatError> {
    let fields: Vec<String> = line.split_ascii_whitespace().map(String::from).collect();
    if fields.len() != loop_keys.len() {
        return Err(FormatError::UnexpectedFieldNumber(
            loop_keys.len(),
            fields.len(),
        ));
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
        .ok_or_else(|| FormatError::MissingField(String::from(key)))?
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
        .ok_or_else(|| FormatError::MissingField(String::from(key)))?
        .clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_ressource;
    use std::fs::File;
    use std::io::BufReader;

    /// Read a file downloaded from the PDB.
    #[test]
    fn test_file_from_pdb() {
        let filepath = test_ressource!("/1bta.cif");
        let file = File::open(filepath).expect("Cound not open test file.");
        let buffer = BufReader::new(file);
        let molecular_system = read_cif(buffer).expect("Error when parsing the file.");
        assert_eq!(molecular_system.atom_count(), 1434);
    }

    #[test]
    fn test_conect_from_cif() {
        let filepath = test_ressource!("/mercaptoethanol.cif");
        let file = File::open(filepath).expect("Could not open test file.");
        let buffer = BufReader::new(file);
        let molecular_system = read_cif(buffer).expect("Error when parsing the file.");
        let reference = vec![
            (0, 1, 1.0),   // 1@C1 - 1@C2
            (0, 2, 1.0),   // 1@C1 - 1@O1
            (1, 3, 1.0),   // 1@C2 - 1@S2
            (3, 7, 1.0),   // 1@S2 - 2@S2
            (4, 5, 1.0),   // 2@C1 - 2@C2
            (4, 6, 1.0),   // 2@C1 - 2@O1
            (5, 7, 1.0),   // 2@C2 - 2@S2
            (8, 9, 1.0),   // 3@C1 - 3@C2
            (8, 13, 1.0),  // 3@C1 - 3@C6
            (9, 10, 1.0),  // 3@C2 - 3@C3
            (10, 11, 1.0), // 3@C3 - 3@C4
            (11, 12, 1.0), // 3@C4 - 3@C5
            (12, 13, 1.0), // 3@C5 - 3@C6
        ];
        assert_eq!(molecular_system.bonds, reference);
    }
}
