use std::{
    collections::HashMap,
    io::{self, BufRead},
    str::FromStr,
};
use components::get_bond_templates;

pub type Position = [f64; 3];
pub type Bond = [isize; 2];

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

/// The parsed fields from an ATOM line in a PDB file.
#[derive(Debug, PartialEq)]
struct PDBLine {
    serial: isize,
    atom_name: String,
    alternate: char,
    residue_name: String,
    chain_identifier: char,
    residue_identifier: isize,
    insertion_code: char,
    position: Position,
    element_symbol: Option<usize>,
}

pub struct MolecularSystem {
    pub names: Vec<String>,
    pub elements: Vec<Option<usize>>,
    pub positions: Vec<Position>,
    pub atom_resindex: Vec<usize>,
    pub resnames: Vec<String>,
    pub resids: Vec<isize>,
    pub bonds: Vec<(usize, usize)>,
}

impl MolecularSystem {
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn iter_residues(&self) -> ResidueIterator {
        ResidueIterator::new(&self)
    }
}

pub struct ResidueIterator<'a> {
    system: &'a MolecularSystem,
    particle_index: usize,
}

impl<'a> ResidueIterator<'a> {
    pub fn new(system: &'a MolecularSystem) -> Self {
        ResidueIterator { system, particle_index: 0 }
    }
}

impl<'a> Iterator for ResidueIterator<'a> {
    type Item = ResidueView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let number_of_atoms = self.system.atom_resindex.len();
        let max_index = number_of_atoms - 1;
        if self.particle_index >= max_index {
            return None;
        }
        let start = self.particle_index;
        let mut index = start;
        let residue_index = self.system.atom_resindex[start];
        let mut current_residue_index = residue_index;
        while index < max_index && current_residue_index == residue_index {
            index += 1;
            current_residue_index = self.system.atom_resindex[index];
        }
        self.particle_index = index;
        Some(ResidueView {
            system: self.system,
            start_index: start,
            next_index: index,
        })
    }
}

pub struct ResidueView<'a> {
    system: &'a MolecularSystem,
    start_index: usize,
    next_index: usize,
}

impl<'a> ResidueView<'a> {
    pub fn find_atom_position(&self, name: &str) -> Option<usize> {
        println!("{:?}", &self.system.names[self.start_index..self.next_index]);
        self.system
            .names[self.start_index..self.next_index]
            .iter()
            .position(|n| name.trim() == n.trim())
            .map(|position| self.start_index + position)
    }

    pub fn name(&self) -> &String {
        let residue_index = self.system.atom_resindex[self.start_index];
        &self.system.resnames[residue_index]
    }
}

pub fn read_pdb<F>(input: F) -> Result<MolecularSystem, ReadError>
where
    F: BufRead,
{
    let atoms = read_pdb_atoms(input)?;
    Ok(add_intra_residue_bonds(flatten_atoms(atoms)))
}

fn read_pdb_atoms<F>(input: F) -> Result<Vec<PDBLine>, ReadError>
where
    F: BufRead,
{
    let mut output = Vec::new();
    for (lineno, line) in input.lines().enumerate() {
        let line = line.map_err(ReadError::IOError)?;
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            let parsed_line =
                parse_pdb_atom_line(&line).map_err(|e| ReadError::FormatError(e, lineno))?;
            output.push(parsed_line);
        }
    }

    Ok(output)
}

fn flatten_atoms(atoms: Vec<PDBLine>) -> MolecularSystem {
    // TODO: Filter alternates
    // TODO: Filter insertions
    if atoms.is_empty() {
        return MolecularSystem {
            names: vec![],
            elements: vec![],
            positions: vec![],
            atom_resindex: vec![],
            resnames: vec![],
            resids: vec![],
            bonds: vec![],
        };
    }
    let mut names = Vec::new();
    let mut elements = Vec::new();
    let mut positions = Vec::new();
    let mut atom_resnames = Vec::new();
    let mut atom_resids = Vec::new();
    let mut atom_insertion_codes = Vec::new();
    let mut alternates = Vec::new();
    atoms.iter().for_each(|atom| {
        names.push(atom.atom_name.clone());
        elements.push(atom.element_symbol);
        positions.push(atom.position);
        atom_resnames.push(atom.residue_name.clone());
        atom_resids.push(atom.residue_identifier);
        atom_insertion_codes.push(atom.insertion_code);
        alternates.push(atom.alternate);
    });

    let mut resnames = vec![atom_resnames[0].clone()];
    let mut resids = vec![atom_resids[0]];
    let mut insertion_codes = vec![atom_insertion_codes[0]];
    let mut current_residue_index = 0;
    let mut atom_resindex = vec![current_residue_index];
    let mut previous = (
        (&atom_resids[0], atom_resnames[0].clone()),
        atom_insertion_codes[0],
    );
    let mut residue_iter = atom_resids
        .iter()
        .zip(atom_resnames)
        .zip(atom_insertion_codes);
    residue_iter.next(); // We already looked at the first residue.
    for residue in residue_iter {
        let ((resid, resname), insertion_code) = &residue;
        if residue != previous {
            current_residue_index += 1;
            resids.push(**resid);
            resnames.push(resname.clone());
            insertion_codes.push(*insertion_code);
        }
        atom_resindex.push(current_residue_index);
        previous = ((resid, resname.to_string()), *insertion_code);
    }

    MolecularSystem {
        names,
        elements,
        positions,
        atom_resindex,
        resnames,
        resids,
        bonds: vec![],
    }
}

fn add_intra_residue_bonds(input: MolecularSystem) -> MolecularSystem {
    let components = get_bond_templates();
    let mut bonds: Vec<(usize, usize)> = Vec::new();
    for residue in input.iter_residues() {
        let residue_name: &str = &residue.name();
        if let Some(templates) = components.get(residue_name) {
            for bond_template in templates {
                let from = residue.find_atom_position(bond_template.from.trim());
                let to = residue.find_atom_position(bond_template.to.trim());
                println!("Bond between {}{from:?} and {}{to:?}", bond_template.from, bond_template.to);
                if let (Some(from), Some(to)) = (from, to) {
                    bonds.push((from, to));
                }
            }
        }
    }
    MolecularSystem {
        names: input.names,
        elements: input.elements,
        positions: input.positions,
        atom_resindex: input.atom_resindex,
        resnames: input.resnames,
        resids: input.resids,
        bonds
    }
}

fn parse_pdb_atom_line(line: &str) -> Result<PDBLine, FormatError> {
    let element_symbol: Option<usize>;
    if line.len() < 54 {
        // A line should be 80 columns. However, nothing after the positions is
        // required. Therefore, it is OK for a line to stop after the positions.
        return Err(FormatError::LineTooShort);
    } else if line.len() < 78 {
        // Since we are looking at the length of the line, we can check if it is
        // long enough to contain an element symbol.
        element_symbol = None;
    } else {
        element_symbol = lookup_element_symbol(&line[76..78]);
    }
    // From this point, we know the line is long enough to retrieve all the fields.
    // We cannot cause a panic by slicing the line.

    let atom_name = String::from(&line[12..16]);
    let alternate = line.chars().nth(16).unwrap();
    let residue_name = String::from(&line[17..20]);
    let chain_identifier = line.chars().nth(21).unwrap();
    let insertion_code = line.chars().nth(27).unwrap();

    let serial = line[6..11]
        .trim()
        .parse()
        .map_err(|_| FormatError::FieldFormat(FieldError::Serial))?;
    let residue_identifier = line[22..26]
        .trim()
        .parse()
        .map_err(|_| FormatError::FieldFormat(FieldError::ResidueIdentifier))?;

    let x: Result<f64, _> = line[30..38].trim().parse();
    let y: Result<f64, _> = line[38..46].trim().parse();
    let z: Result<f64, _> = line[46..54].trim().parse();
    if x.is_err() || y.is_err() || z.is_err() {
        return Err(FormatError::FieldFormat(FieldError::Position));
    }
    // We use nanometers instead of angstrom.
    let position = [x.unwrap() / 10.0, y.unwrap() / 10.0, z.unwrap() / 10.0];

    Ok(PDBLine {
        serial,
        atom_name,
        alternate,
        residue_name,
        chain_identifier,
        residue_identifier,
        insertion_code,
        position,
        element_symbol,
    })
}

fn lookup_element_symbol(symbol: &str) -> Option<usize> {
    #[rustfmt::skip]
    let element_number = match symbol.trim() {
        "H"  =>  1,                                                                                                                                                                                                               "He" =>   2,
        "Li" =>  3, "Be" =>  4,                                                                                                                                  "B"  =>   5, "C"  =>   6, "N"  =>   7, "O"  =>   8, "F"  =>   9, "Ne" =>  10,
        "Na" => 11, "Mg" => 12,                                                                                                                                  "Al" =>  13, "Si" =>  14, "P"  =>  15, "S"  =>  16, "Cl" =>  17, "Ar" =>  18,
        "K"  => 19, "Ca" => 20, "Sc" => 21, "Ti" =>  22, "V"  =>  23, "Cr" =>  24, "Mn" =>  25, "Fe" =>  26, "Co" =>  27, "Ni" =>  28, "Cu" =>  29, "Zn" =>  30, "Ga" =>  31, "Ge" =>  32, "As" =>  33, "Se" =>  34, "Br" =>  35, "Kr" =>  36,
        "Rb" => 37, "Sr" => 38, "Y"  => 39, "Zr" =>  40, "Nb" =>  41, "Mo" =>  42, "Tc" =>  43, "Ru" =>  44, "Rh" =>  45, "Pd" =>  46, "Ag" =>  47, "Cd" =>  48, "In" =>  49, "Sn" =>  50, "Sb" =>  51, "Te" =>  52, "I"  =>  53, "Xe" =>  54,
        "Cs" => 55, "Ba" => 56,             "Hf" =>  72, "Ta" =>  73, "W"  =>  74, "Re" =>  75, "Os" =>  76, "Ir" =>  77, "Pt" =>  78, "Au" =>  79, "Hg" =>  80, "Tl" =>  81, "Pb" =>  82, "Bi" =>  83, "Po" =>  84, "At" =>  85, "Rn" =>  86,
        "Fr" => 87, "Ra" => 88,             "Rf" => 104, "Db" => 105, "Sg" => 106, "Bh" => 107, "Hs" => 108, "Mt" => 109, "Ds" => 110, "Rg" => 111, "Cn" => 112, "Nh" => 113, "Fl" => 114, "Mc" => 115, "Lv" => 116, "Ts" => 117, "Og" => 118,
        /* Lanthanides */                   "La" =>  57, "Ce" =>  58, "Pr" =>  59, "Nd" =>  60, "Pm" =>  61, "Sm" =>  62, "Eu" =>  63, "Gd" =>  64, "Tb" =>  65, "Dy" =>  66, "Ho" =>  67, "Er" =>  68, "Tm" =>  69, "Yb" =>  70, "Lu" =>  71,
        /* Actinides   */                   "Ac" =>  89, "Th" =>  90, "Pa" =>  91, "U"  =>  92, "Np" =>  93, "Pu" =>  94, "Am" =>  95, "Cm" =>  96, "Bk" =>  97, "Cf" =>  98, "Es" =>  99, "Fm" => 100, "Md" => 101, "No" => 102, "Lr" => 103,  
        _ => 0,        
    };
    if element_number == 0 {
        None
    } else {
        Some(element_number)
    }
}

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
    // TODO: Implement data and save blocks.
    // We read a data block with the name in argument.
    // Data(String)
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
                    // We do not read data blocks. If we encounter one, it
                    // means the rest of the file will be only data blocks
                    // until EOF. Therefore, there is no point to keep reading.
                    Some(first) if first.starts_with("data_") => break,
                    Some(first) if first == "loop_" => PDBXContext::LoopKey(None),
                    // Data items have a first token starting with '_' and the rest of the
                    // tokens are the value. We do not read any at the moment so we ignore
                    // these lines.
                    Some(first) if first.starts_with('_') => PDBXContext::Idle,
                    Some(first) if first.starts_with('#') => PDBXContext::Idle,
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

    Ok(flatten_atoms(atoms))
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case("H", Some(1))]
    #[case("C", Some(6))]
    #[case("  C ", Some(6))] // Are spaces trimmed correctly?
    #[case("not an element", None)]
    fn test_lookup_element_symbol(#[case] symbol: &str, #[case] expected_number: Option<usize>) {
        let number = lookup_element_symbol(symbol);
        assert_eq!(number, expected_number);
    }

    #[rstest]
    // Full line
    #[case("ATOM      1  N   GLY A   3      17.119   0.186  36.320  1.00 64.10           N  ")]
    // Truncated after the element symbol
    #[case("ATOM      1  N   GLY A   3      17.119   0.186  36.320  1.00 64.10           N")]
    fn test_parse_pdb_atom_line(#[case] line: &str) {
        let expected = PDBLine {
            serial: 1,
            atom_name: String::from(" N  "),
            alternate: ' ',
            residue_name: String::from("GLY"),
            chain_identifier: 'A',
            residue_identifier: 3,
            insertion_code: ' ',
            position: [1.7119, 0.0186, 3.6320],
            element_symbol: Some(7),
        };
        let atom = parse_pdb_atom_line(line).unwrap();
        assert_eq!(atom, expected);
    }
}
