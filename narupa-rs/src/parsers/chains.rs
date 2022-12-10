use crate::parsers::residues::ResidueIterator;
use crate::parsers::MolecularSystem;

pub struct ChainView<'a> {
    system: &'a MolecularSystem,
    start_index: usize,
    next_index: usize,
}

impl<'a> ChainView<'a> {
    pub fn iter_residues(&self) -> ResidueIterator {
        ResidueIterator::new(&self.system, self.start_index, self.next_index)
    }
}

pub struct ChainIterator<'a> {
    system: &'a MolecularSystem,
    particle_index: usize,
}

impl<'a> ChainIterator<'a> {
    pub fn new(system: &'a MolecularSystem) -> Self {
        Self {
            system,
            particle_index: 0,
        }
    }
}

impl<'a> Iterator for ChainIterator<'a> {
    type Item = ChainView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.particle_index >= self.system.atom_count() {
            return None;
        }

        let start = self.particle_index;
        let current_residue = self.system.atom_resindex[start];
        let current_chain = self.system.residue_chain_index[current_residue];
        let mut i = 0;
        for atom in self.particle_index..self.system.atom_count() {
            i = atom;
            let residue = self.system.atom_resindex[atom];
            let chain = self.system.residue_chain_index[residue];
            if residue != current_residue && chain != current_chain {
                break;
            }
        }
        self.particle_index = i + 1;

        Some(ChainView {
            system: self.system,
            start_index: start,
            next_index: self.particle_index,
        })
    }
}
