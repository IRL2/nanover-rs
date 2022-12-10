use crate::parsers::chains::ChainIterator;
use crate::parsers::line::PDBLine;
use crate::parsers::residues::{ResidueIterator, ResidueView};
use crate::parsers::Position;
use components::{get_bond_templates, ResidueType};

pub struct MolecularSystem {
    pub names: Vec<String>,
    pub elements: Vec<Option<usize>>,
    pub positions: Vec<Position>,
    pub atom_resindex: Vec<usize>,
    pub resnames: Vec<String>,
    pub resids: Vec<isize>,
    pub residue_chain_index: Vec<usize>,
    pub chain_identifiers: Vec<String>,
    pub bonds: Vec<(usize, usize)>,
}

impl Default for MolecularSystem {
    fn default() -> Self {
        MolecularSystem {
            names: vec![],
            elements: vec![],
            positions: vec![],
            atom_resindex: vec![],
            resnames: vec![],
            resids: vec![],
            residue_chain_index: vec![],
            chain_identifiers: vec![],
            bonds: vec![],
        }
    }
}

impl MolecularSystem {
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn iter_residues(&self) -> ResidueIterator {
        ResidueIterator::new(&self, 0, self.atom_count())
    }

    pub fn iter_chains(&self) -> ChainIterator {
        ChainIterator::new(&self)
    }

    pub fn add_intra_residue_bonds(mut self) -> MolecularSystem {
        let components = get_bond_templates();
        let mut bonds: Vec<(usize, usize)> = Vec::new();
        for residue in self.iter_residues() {
            let residue_name: &str = &residue.name();
            let Some(templates) = components.get(residue_name) else {
                continue;
            };
            for bond_template in &templates.bonds {
                let from = residue.find_atom_position(bond_template.from.trim());
                let to = residue.find_atom_position(bond_template.to.trim());
                if let (Some(from), Some(to)) = (from, to) {
                    bonds.push((from, to));
                }
            }
        }
        self.bonds.append(&mut bonds);
        self
    }

    pub fn add_inter_residue_bonds(mut self) -> MolecularSystem {
        let components = get_bond_templates();
        let mut bonds: Vec<(usize, usize)> = Vec::new();

        for chain in self.iter_chains() {
            let previous_residues = chain.iter_residues();
            let current_residues = chain.iter_residues().skip(1);
            for (previous, current) in previous_residues.zip(current_residues) {
                let template_previous = components.get(previous.name());
                let template_current = components.get(current.name());
                let (Some(template_previous), Some(template_current)) = (template_previous, template_current) else {
                    continue;
                };
                let type_previous = &template_previous.residue_type;
                let type_current = &template_current.residue_type;
                if type_previous != type_current {
                    continue;
                }
                match type_previous {
                    ResidueType::Peptide => {
                        bonds.append(&mut make_peptide_bond(&previous, &current))
                    }
                    _ => {}
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

fn make_peptide_bond(nter: &ResidueView, cter: &ResidueView) -> Vec<(usize, usize)> {
    let maybe_c_on_nter = nter.find_atom_position("C");
    let maybe_n_on_cter = cter.find_atom_position("N");
    match (maybe_c_on_nter, maybe_n_on_cter) {
        (Some(c_on_nter), Some(n_on_cter)) => {
            vec![(c_on_nter, n_on_cter)]
        }
        _ => vec![],
    }
}
