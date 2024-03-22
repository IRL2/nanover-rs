use crate::parsers::chains::ChainIterator;
use crate::parsers::line::PDBLine;
use crate::parsers::residues::{ResidueIterator, ResidueView};
use crate::parsers::Position;
use components::{get_bond_templates, ResidueType};

#[derive(Default)]
pub struct MolecularSystem {
    pub names: Vec<String>,
    pub elements: Vec<Option<usize>>,
    pub positions: Vec<Position>,
    pub atom_resindex: Vec<usize>,
    pub resnames: Vec<String>,
    pub resids: Vec<isize>,
    pub residue_chain_index: Vec<usize>,
    pub chain_identifiers: Vec<String>,
    pub bonds: Vec<(usize, usize, f32)>,
}

impl MolecularSystem {
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn iter_residues(&self) -> ResidueIterator {
        ResidueIterator::new(self, 0, self.atom_count())
    }

    pub fn iter_chains(&self) -> ChainIterator {
        ChainIterator::new(self)
    }

    pub fn add_intra_residue_bonds(mut self) -> MolecularSystem {
        let components = get_bond_templates();
        let mut bonds: Vec<(usize, usize, f32)> = Vec::new();
        for residue in self.iter_residues() {
            let residue_name: &str = residue.name();
            let Some(templates) = components.get(residue_name) else {
                continue;
            };
            for bond_template in &templates.bonds {
                let from = residue.find_atom_position(bond_template.from.trim());
                let to = residue.find_atom_position(bond_template.to.trim());
                let order = bond_template.order;
                if let (Some(from), Some(to)) = (from, to) {
                    bonds.push((from, to, order));
                }
            }
        }
        self.bonds.append(&mut bonds);
        self
    }

    pub fn add_inter_residue_bonds(mut self) -> MolecularSystem {
        let components = get_bond_templates();
        let mut bonds: Vec<(usize, usize, f32)> = Vec::new();

        for chain in self.iter_chains() {
            let previous_residues = chain.iter_residues();
            let current_residues = chain.iter_residues().skip(1);
            for (previous, current) in previous_residues.zip(current_residues) {
                let template_previous = components.get(previous.name());
                let template_current = components.get(current.name());
                let (Some(template_previous), Some(template_current)) =
                    (template_previous, template_current)
                else {
                    continue;
                };
                let type_previous = &template_previous.residue_type;
                let type_current = &template_current.residue_type;
                if type_previous != type_current {
                    continue;
                }
                if type_previous == &ResidueType::Peptide {
                    bonds.append(&mut make_peptide_bond(&previous, &current))
                }
            }
        }
        self.bonds.append(&mut bonds);
        self
    }
}

impl From<Vec<PDBLine>> for MolecularSystem {
    fn from(atoms: Vec<PDBLine>) -> Self {
        // TODO: Filter alternates
        // TODO: Filter insertions
        if atoms.is_empty() {
            return MolecularSystem::default();
        }

        let mut names = Vec::new();
        let mut elements = Vec::new();
        let mut positions = Vec::new();
        let mut atom_resnames = Vec::new();
        let mut atom_resids = Vec::new();
        let mut atom_insertion_codes = Vec::new();
        let mut alternates = Vec::new();
        let mut atom_chain_identifiers = Vec::new();
        atoms.iter().for_each(|atom| {
            names.push(atom.atom_name.clone());
            elements.push(atom.element_symbol);
            positions.push(atom.position);
            atom_resnames.push(atom.residue_name.clone());
            atom_resids.push(atom.residue_identifier);
            atom_insertion_codes.push(atom.insertion_code);
            alternates.push(atom.alternate);
            atom_chain_identifiers.push(atom.chain_identifier);
        });

        let mut current_residue_index = 0;
        let mut atom_resindex = vec![current_residue_index];
        let mut resnames = vec![atom_resnames[0].clone()];
        let mut resids = vec![atom_resids[0]];
        let mut insertion_codes = vec![atom_insertion_codes[0]];

        let mut current_chain_index = 0;
        let mut residue_chain_index = vec![current_chain_index];
        let mut current_chain_identifier = atom_chain_identifiers[0];
        let mut chain_identifiers = vec![String::from(atom_chain_identifiers[0])];
        let mut residue_iter = atom_resids
            .iter()
            .zip(atom_resnames)
            .zip(atom_insertion_codes)
            .zip(atom_chain_identifiers);
        let mut previous = residue_iter.next().unwrap();
        for residue in residue_iter {
            let (((resid, resname), insertion_code), chain_identifier) = &residue;
            if residue != previous {
                current_residue_index += 1;
                resids.push(**resid);
                resnames.push(resname.clone());
                insertion_codes.push(*insertion_code);
                if chain_identifier != &current_chain_identifier {
                    current_chain_identifier = *chain_identifier;
                    current_chain_index += 1;
                    chain_identifiers.push(String::from(current_chain_identifier));
                }
                residue_chain_index.push(current_chain_index);
            }
            atom_resindex.push(current_residue_index);
            previous = residue;
        }

        MolecularSystem {
            names,
            elements,
            positions,
            atom_resindex,
            resnames,
            resids,
            residue_chain_index,
            chain_identifiers,
            bonds: vec![],
        }
    }
}

fn make_peptide_bond(nter: &ResidueView, cter: &ResidueView) -> Vec<(usize, usize, f32)> {
    let maybe_c_on_nter = nter.find_atom_position("C");
    let maybe_n_on_cter = cter.find_atom_position("N");
    match (maybe_c_on_nter, maybe_n_on_cter) {
        (Some(c_on_nter), Some(n_on_cter)) => {
            vec![(c_on_nter, n_on_cter, 1.0)]
        }
        _ => vec![],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::fixture;
    use rstest::rstest;

    /// Create a MolecularSystem with two three-peptides on separate chains.
    ///
    /// Each peptide contains three alanines. The bonds are not included as
    /// the fixture is meant to test if the bonds can be generated correctly.
    /// Note: the coordinates are not meant to be usable!
    #[fixture]
    fn two_peptides_without_bonds() -> MolecularSystem {
        let atoms_per_residue = 12;
        let number_of_residues = 6;
        let residues_per_chain = 3;
        let number_of_chains = 2;
        let names_for_one_residue = make_string_vector(vec![
            "N", "CA", "C", "O", "CB", "H", "H2", "H3", "HA", "HB1", "HB2", "HB3",
        ]);
        let elements_for_one_residue = vec![7, 6, 6, 8, 6, 1, 1, 1, 1, 1, 1, 1]
            .into_iter()
            .map(Some)
            .collect();
        let positions = std::iter::repeat([0.0; 3])
            .take(atoms_per_residue * number_of_residues)
            .collect();
        let resnames = std::iter::repeat(String::from("ALA"))
            .take(number_of_residues)
            .collect();
        let resids = (1..=number_of_residues as isize).collect();
        let chain_identifiers = make_string_vector(vec!["A", "B"]);
        MolecularSystem {
            names: tile(names_for_one_residue, number_of_residues),
            elements: tile(elements_for_one_residue, number_of_residues),
            positions,
            atom_resindex: repeat((0..number_of_residues).collect(), atoms_per_residue),
            resnames,
            resids,
            residue_chain_index: repeat((0..number_of_chains).collect(), residues_per_chain),
            chain_identifiers,
            bonds: vec![],
        }
    }

    #[rstest]
    fn test_intra_residue_bonds(two_peptides_without_bonds: MolecularSystem) {
        let number_of_residues = 6;
        let number_of_atoms_per_residue = 12;
        let reference_for_one_residue = vec![
            (0, 1, 1.0),  // N-CA
            (0, 5, 1.0),  // N-H
            (0, 6, 1.0),  // N-H2
            (1, 2, 1.0),  // CA-C
            (1, 4, 1.0),  // CA-CB
            (1, 8, 1.0),  // CA-HA
            (2, 3, 2.0),  // C-O
            (4, 9, 1.0),  // CB-HB1
            (4, 10, 1.0), // CB-HB2
            (4, 11, 1.0), // CB-HB3
            (0, 7, 1.0),  // N-H3 use an extra bond in components
        ];
        let mut reference = Vec::new();
        for residue_index in 0..number_of_residues {
            let offset = number_of_atoms_per_residue * residue_index;
            for bond in &reference_for_one_residue {
                reference.push((bond.0 + offset, bond.1 + offset, bond.2));
            }
        }

        let two_peptides = two_peptides_without_bonds.add_intra_residue_bonds();
        assert_eq!(two_peptides.bonds, reference);
    }

    #[rstest]
    fn test_inter_residue_bonds(two_peptides_without_bonds: MolecularSystem) {
        let reference = vec![
            // Chain A
            (2, 12, 1.0),
            (14, 24, 1.0),
            // Chain B
            (38, 48, 1.0),
            (50, 60, 1.0),
        ];

        let two_peptides = two_peptides_without_bonds.add_inter_residue_bonds();
        assert_eq!(two_peptides.bonds, reference);
    }

    #[rstest]
    fn test_residue_iterator(two_peptides_without_bonds: MolecularSystem) {
        let reference = vec![(0, 12), (12, 24), (24, 36), (36, 48), (48, 60), (60, 72)];
        let residues: Vec<(usize, usize)> = two_peptides_without_bonds
            .iter_residues()
            .map(|r| (r.start_index, r.next_index))
            .collect();
        assert_eq!(residues, reference);
    }

    #[rstest]
    fn test_chain_iterator(two_peptides_without_bonds: MolecularSystem) {
        let number_of_chains = 2;
        let number_of_residues = 6;
        let number_of_atoms_per_residue = 12;
        let chain_b_start_index =
            (number_of_residues / number_of_chains) * number_of_atoms_per_residue;
        let reference = vec![
            (0, chain_b_start_index),
            (chain_b_start_index, two_peptides_without_bonds.atom_count()),
        ];
        let chains: Vec<(usize, usize)> = two_peptides_without_bonds
            .iter_chains()
            .map(|c| (c.start_index, c.next_index))
            .collect();
        assert_eq!(chains, reference);
    }

    fn tile<T>(input: Vec<T>, number: usize) -> Vec<T>
    where
        T: Clone,
    {
        input
            .iter()
            .cycle()
            .take(input.len() * number)
            .cloned()
            .collect()
    }

    fn repeat<T>(input: Vec<T>, number: usize) -> Vec<T>
    where
        T: Clone,
    {
        let mut output = Vec::new();
        for element in input {
            for _ in 0..number {
                output.push(element.clone());
            }
        }
        output
    }

    fn make_string_vector(input: Vec<&str>) -> Vec<String> {
        input.iter().map(|element| String::from(*element)).collect()
    }
}
