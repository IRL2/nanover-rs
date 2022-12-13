use std::io::BufRead;
use crate::parsers::chemistry::lookup_element_symbol;
use crate::parsers::errors::*;
use crate::parsers::molecular_system::MolecularSystem;
use crate::parsers::line::PDBLine;


pub fn read_pdb<F>(input: F) -> Result<MolecularSystem, ReadError>
where
    F: BufRead,
{
    let atoms = read_pdb_atoms(input)?;
    Ok(MolecularSystem::from(atoms)
        .add_intra_residue_bonds()
        .add_inter_residue_bonds()
    )
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

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