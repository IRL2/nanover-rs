pub type Position = [f64; 3];
pub type Bond = [isize; 2];

pub mod chemistry;
pub mod errors;
pub mod molecular_system;
pub mod line;
pub mod pdb;
pub mod cif;
pub mod residues;
pub mod chains;

pub use molecular_system::MolecularSystem;
pub use pdb::read_pdb;
pub use cif::read_cif;