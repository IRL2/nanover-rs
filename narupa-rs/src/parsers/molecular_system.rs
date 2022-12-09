use components::get_bond_templates;
use crate::parsers::Position;
use crate::parsers::line::PDBLine;

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

impl MolecularSystem {
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    pub fn iter_residues(&self) -> ResidueIterator {
        ResidueIterator::new(&self)
    }

    pub fn add_intra_residue_bonds(mut self) -> MolecularSystem {
        let components = get_bond_templates();
        let mut bonds: Vec<(usize, usize)> = Vec::new();
        for residue in self.iter_residues() {
            let residue_name: &str = &residue.name();
            if let Some(templates) = components.get(residue_name) {
                for bond_template in &templates.bonds {
                    let from = residue.find_atom_position(bond_template.from.trim());
                    let to = residue.find_atom_position(bond_template.to.trim());
                    if let (Some(from), Some(to)) = (from, to) {
                        bonds.push((from, to));
                    }
                }
            }
        }
        self.bonds.append(&mut bonds);
        self
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

impl From<Vec<PDBLine>> for MolecularSystem {
    fn from(atoms: Vec<PDBLine>) -> Self {
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
                residue_chain_index: vec![],
                chain_identifiers: vec![],
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
        let mut previous = (
            (
                (&atom_resids[0], atom_resnames[0].clone()),
                atom_insertion_codes[0],
            ),
            atom_chain_identifiers[0],            
        );
        let mut residue_iter = atom_resids
            .iter()
            .zip(atom_resnames)
            .zip(atom_insertion_codes)
            .zip(atom_chain_identifiers);
        residue_iter.next(); // We already looked at the first residue.
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
            previous = (((resid, resname.to_string()), *insertion_code), *chain_identifier);
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