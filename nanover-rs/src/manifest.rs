use std::fs::File;
use std::io::BufReader;
use thiserror::Error;

use crate::simulation::{OpenMMSimulation, XMLParsingError};

#[derive(Error, Debug)]
#[error("Something went wrong while reading the manifest.")]
pub struct ManifestReadError {}

#[derive(Error, Debug)]
pub enum LoadDefaultError {
    #[error("No default defined in the manifest.")]
    NoDefault,
    #[error("{0}")]
    LoadSimulationError(#[from] LoadSimulationError),
}

#[derive(Error, Debug)]
pub enum LoadNextError {
    #[error("No next simulation defined in the manifest.")]
    NoNext,
    #[error("{0}")]
    LoadSimulationError(#[from] LoadSimulationError),
}

#[derive(Error, Debug)]
pub enum LoadSimulationError {
    #[error("The simulation file cannot be openned: {0}")]
    CannotOpen(#[from] std::io::Error),
    #[error("Index {0} is not defined in the manifest.")]
    NoIndex(usize),
    #[error("{0}")]
    XMLParsingError(#[from] XMLParsingError),
}

pub type SimulationEntry = String;

enum ManifestContent {
    SingleBuffer(Vec<u8>),
    MultiplePath(Vec<SimulationEntry>),
}

pub struct Manifest {
    content: ManifestContent,
    current_index: Option<usize>,
}

impl Manifest {
    pub fn from_simulation_xml_paths<T>(files: Vec<T>) -> Self
    where
        T: ToString,
    {
        let entries = files.iter().map(|path| path.to_string()).collect();
        let content = ManifestContent::MultiplePath(entries);
        let current_index = None;
        Self {
            content,
            current_index,
        }
    }

    pub fn from_simulation_xml_bytes(bytes: &[u8]) -> Self {
        Self {
            content: ManifestContent::SingleBuffer(bytes.to_vec()),
            current_index: None,
        }
    }

    pub fn load_default(&mut self) -> Result<OpenMMSimulation, LoadDefaultError> {
        let default_index = self
            .get_default_index()
            .ok_or(LoadDefaultError::NoDefault)?;
        Ok(self.load_index(default_index)?)
    }

    pub fn load_index(&mut self, index: usize) -> Result<OpenMMSimulation, LoadSimulationError> {
        match &self.content {
            ManifestContent::SingleBuffer(_) if index != 0 => {
                Err(LoadSimulationError::NoIndex(index))
            }
            ManifestContent::SingleBuffer(bytes) => {
                let buffer = BufReader::new(bytes.as_slice());
                let simulation = OpenMMSimulation::from_xml(buffer)?;
                self.current_index = Some(index);
                Ok(simulation)
            }
            ManifestContent::MultiplePath(_) => {
                let path = self
                    .get_path_for_index(index)
                    .ok_or(LoadSimulationError::NoIndex(index))?
                    .clone();
                self.current_index = Some(index);
                let xml_file = File::open(path)?;
                let buffer = BufReader::new(xml_file);
                let simulation = OpenMMSimulation::from_xml(buffer)?;
                Ok(simulation)
            }
        }
    }

    pub fn load_next(&mut self) -> Result<OpenMMSimulation, LoadNextError> {
        let next_index = self.get_next_index().ok_or(LoadNextError::NoNext)?;
        Ok(self.load_index(next_index)?)
    }

    pub fn list_simulations(&self) -> Vec<SimulationEntry> {
        match &self.content {
            ManifestContent::SingleBuffer(_) => vec!["default".into()],
            ManifestContent::MultiplePath(content) => content.clone(),
        }
    }

    fn get_default_index(&self) -> Option<usize> {
        match &self.content {
            ManifestContent::SingleBuffer(_) => Some(0),
            ManifestContent::MultiplePath(content) if content.is_empty() => None,
            ManifestContent::MultiplePath(_) => Some(0),
        }
    }

    fn get_next_index(&self) -> Option<usize> {
        match (self.current_index, &self.content) {
            (_, ManifestContent::SingleBuffer(_)) => None,
            (_, ManifestContent::MultiplePath(content)) if content.is_empty() => None,
            (None, ManifestContent::MultiplePath(_)) => Some(0),
            (Some(index), ManifestContent::MultiplePath(content)) if index < content.len() - 1 => {
                Some(index + 1)
            }
            (Some(_), ManifestContent::MultiplePath(_)) => Some(0),
        }
    }

    fn get_path_for_index(&self, index: usize) -> Option<&SimulationEntry> {
        match &self.content {
            ManifestContent::SingleBuffer(_) => None,
            ManifestContent::MultiplePath(content) => content.get(index),
        }
    }
}
