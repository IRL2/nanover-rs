pub type Position = [f64; 3];
pub type Bond = (usize, usize, f32); // from, to, order

pub mod chains;
pub mod chemistry;
pub mod cif;
pub mod errors;
pub mod line;
pub mod molecular_system;
pub mod pdb;
pub mod residues;

pub use cif::read_cif;
pub use molecular_system::MolecularSystem;
pub use pdb::read_pdb;
