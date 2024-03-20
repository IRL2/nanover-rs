use std::fs::File;
use std::io::BufReader;
use thiserror::Error;

use crate::{
    application::InputPath,
    recording::{RecordingReadError, ReplaySimulation},
    simulation::{OpenMMSimulation, Simulation, ToFrameData, XMLParsingError},
};

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
    #[error("The input file is not a NanoVer recording.")]
    NotARecording,
    #[error("Recording file version {0} is not supported.")]
    UnsuportedFormatVersion(u64),
    #[error("An error occured while reading a recording file: {0}")]
    ReadError(std::io::Error),
}

impl From<RecordingReadError> for LoadSimulationError {
    fn from(value: RecordingReadError) -> Self {
        match value {
            RecordingReadError::IOError(error) => Self::ReadError(error),
            RecordingReadError::NotARecording => Self::NotARecording,
            RecordingReadError::UnsuportedFormatVersion(version) => {
                Self::UnsuportedFormatVersion(version)
            }
        }
    }
}

enum ManifestContent {
    SingleBuffer(Vec<u8>),
    MultiplePath(Vec<InputPath>),
}

pub enum LoadedSimulation {
    OpenMM(OpenMMSimulation),
    Recording(ReplaySimulation),
}

impl Simulation for LoadedSimulation {
    fn step(&mut self, steps: i32) {
        match self {
            Self::OpenMM(simulation) => simulation.step(steps),
            Self::Recording(simulation) => simulation.step(steps),
        }
    }

    fn reset(&mut self) {
        match self {
            Self::OpenMM(simulation) => simulation.reset(),
            Self::Recording(simulation) => simulation.reset(),
        }
    }
}

impl ToFrameData for LoadedSimulation {
    fn to_framedata(
        &self,
        with_velocity: bool,
        with_forces: bool,
    ) -> nanover_proto::trajectory::FrameData {
        match self {
            Self::OpenMM(simulation) => simulation.to_framedata(with_velocity, with_forces),
            Self::Recording(simulation) => simulation.to_framedata(with_velocity, with_forces),
        }
    }

    fn to_topology_framedata(&self) -> nanover_proto::trajectory::FrameData {
        match self {
            Self::OpenMM(simulation) => simulation.to_topology_framedata(),
            Self::Recording(simulation) => simulation.to_topology_framedata(),
        }
    }
}

pub struct Manifest {
    content: ManifestContent,
    current_index: Option<usize>,
}

impl Manifest {
    pub fn from_simulation_input_paths(files: Vec<InputPath>) -> Self {
        let content = ManifestContent::MultiplePath(files);
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

    pub fn load_default(&mut self) -> Result<LoadedSimulation, LoadDefaultError> {
        let default_index = self
            .get_default_index()
            .ok_or(LoadDefaultError::NoDefault)?;
        Ok(self.load_index(default_index)?)
    }

    pub fn load_index(&mut self, index: usize) -> Result<LoadedSimulation, LoadSimulationError> {
        match &self.content {
            ManifestContent::SingleBuffer(_) if index != 0 => {
                Err(LoadSimulationError::NoIndex(index))
            }
            ManifestContent::SingleBuffer(bytes) => {
                let buffer = BufReader::new(bytes.as_slice());
                let simulation = OpenMMSimulation::from_xml(buffer)?;
                self.current_index = Some(index);
                Ok(LoadedSimulation::OpenMM(simulation))
            }
            ManifestContent::MultiplePath(_) => {
                self.current_index = Some(index);
                let entry = self
                    .get_path_for_index(index)
                    .ok_or(LoadSimulationError::NoIndex(index))?;
                match entry {
                    InputPath::OpenMM(path) => {
                        let xml_file = File::open(path)?;
                        let buffer = BufReader::new(xml_file);
                        let simulation = OpenMMSimulation::from_xml(buffer)?;
                        Ok(LoadedSimulation::OpenMM(simulation))
                    }
                    InputPath::Recording(path) => {
                        let simulation = path.try_into()?;
                        Ok(LoadedSimulation::Recording(simulation))
                    }
                }
            }
        }
    }

    pub fn load_next(&mut self) -> Result<LoadedSimulation, LoadNextError> {
        let next_index = self.get_next_index().ok_or(LoadNextError::NoNext)?;
        Ok(self.load_index(next_index)?)
    }

    pub fn list_simulations(&self) -> Vec<String> {
        match &self.content {
            ManifestContent::SingleBuffer(_) => vec!["default".into()],
            ManifestContent::MultiplePath(content) => content
                .iter()
                .map(|entry| match entry {
                    InputPath::OpenMM(path) => path.clone(),
                    InputPath::Recording(path) => path.to_string(),
                })
                .collect(),
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

    fn get_path_for_index(&self, index: usize) -> Option<&InputPath> {
        match &self.content {
            ManifestContent::SingleBuffer(_) => None,
            ManifestContent::MultiplePath(content) => content.get(index),
        }
    }
}
