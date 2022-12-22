use crate::parsers::chemistry::lookup_element_symbol;
use crate::parsers::errors::*;
use crate::parsers::line::PDBLine;
use crate::parsers::molecular_system::MolecularSystem;
use crate::parsers::Bond;
use std::collections::HashMap;
use std::io::BufRead;

pub fn read_pdb<F>(input: F) -> Result<MolecularSystem, ReadError>
where
    F: BufRead,
{
    let (atoms, bonds) = read_pdb_atoms(input)?;
    let mut atoms = MolecularSystem::from(atoms);
    atoms.bonds.extend(bonds);
    Ok(atoms.add_intra_residue_bonds().add_inter_residue_bonds())
}

fn read_pdb_atoms<F>(input: F) -> Result<(Vec<PDBLine>, Vec<Bond>), ReadError>
where
    F: BufRead,
{
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for (lineno, line) in input.lines().enumerate() {
        let line = line.map_err(ReadError::IOError)?;
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            let parsed_line =
                parse_pdb_atom_line(&line).map_err(|e| ReadError::FormatError(e, lineno))?;
            atoms.push(parsed_line);
        } else if line.starts_with("CONECT") {
            let line_bonds =
                parse_pdb_conect_line(&line).map_err(|e| ReadError::FormatError(e, lineno))?;
            bonds.extend(line_bonds);
        }
    }

    let bonds = match_bonds_to_atoms(&bonds, &atoms);

    Ok((atoms, bonds))
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

fn parse_pdb_conect_line(line: &str) -> Result<Vec<(isize, isize)>, FormatError> {
    let mut line_bonds = Vec::new();
    let base = line
        .get(6..11)
        .ok_or(FormatError::FieldFormat(FieldError::Conect))?;
    let base = base
        .trim()
        .parse()
        .map_err(|_| FormatError::FieldFormat(FieldError::Conect))?;
    for index in (11..31).step_by(5) {
        let other = line.get(index..(index + 5));
        let Some(other) = other else {break};
        let other = other
            .trim()
            .parse()
            .map_err(|_| FormatError::FieldFormat(FieldError::Conect))?;
        line_bonds.push((base, other));
    }
    Ok(line_bonds)
}

fn match_bonds_to_atoms(
    originals: &[(isize, isize)],
    atoms: &[PDBLine],
) -> Vec<(usize, usize, f32)> {
    let serial_to_index = atoms
        .iter()
        .enumerate()
        .map(|(index, atom)| (atom.serial, index));
    let serial_to_index: HashMap<isize, usize> = HashMap::from_iter(serial_to_index);
    originals
        .iter()
        .filter_map(|bond| {
            let from = serial_to_index.get(&bond.0);
            let to = serial_to_index.get(&bond.1);
            from.zip(to).map(|(from, to)| (*from, *to, 1.0))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::test_ressource;
    use rstest::rstest;
    use std::fs::File;
    use std::io::BufReader;

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

    #[test]
    fn test_conect() {
        let filepath = test_ressource!("/mercaptoethanol.pdb");
        let file = File::open(filepath).expect("Cound not open test file.");
        let buffer = BufReader::new(file);
        let molecular_system = read_pdb(buffer).expect("Error parsing the file.");
        let reference = vec![
            (0, 1, 1.0),
            (0, 2, 1.0), // 1287 - 0
            (1, 0, 1.0),
            (1, 3, 1.0), // 1288 - 1
            (2, 0, 1.0), // 1289 - 2
            (3, 1, 1.0),
            (3, 7, 1.0), // 1290 - 3
            (4, 5, 1.0),
            (4, 6, 1.0), // 1291 - 4
            (5, 4, 1.0),
            (5, 7, 1.0), // 1292 - 5
            (6, 4, 1.0), // 1293 - 6
            (7, 3, 1.0),
            (7, 5, 1.0), // 1294 - 7
            (8, 9, 1.0),
            (8, 13, 1.0), // 1295 - 8
            (9, 8, 1.0),
            (9, 10, 1.0), // 1296 - 9
            (10, 9, 1.0),
            (10, 11, 1.0), // 1297 - 10
            (11, 10, 1.0),
            (11, 12, 1.0), // 1298 - 11
            (12, 11, 1.0),
            (12, 13, 1.0), // 1299 - 12
            (13, 8, 1.0),
            (13, 12, 1.0),
        ];
        assert_eq!(molecular_system.bonds, reference);
    }
}
