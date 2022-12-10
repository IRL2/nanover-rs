use crate::parsers::MolecularSystem;

pub struct ResidueIterator<'a> {
    system: &'a MolecularSystem,
    particle_index: usize,
    end: usize,
}

impl<'a> ResidueIterator<'a> {
    pub fn new(system: &'a MolecularSystem, start: usize, end: usize) -> Self {
        ResidueIterator {
            system,
            particle_index: start,
            end: end,
        }
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
        while index < max_index && index < self.end && current_residue_index == residue_index {
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
        self.system.names[self.start_index..self.next_index]
            .iter()
            .position(|n| name.trim() == n.trim())
            .map(|position| self.start_index + position)
    }

    pub fn name(&self) -> &str {
        let residue_index = self.system.atom_resindex[self.start_index];
        &self.system.resnames[residue_index]
    }
}
