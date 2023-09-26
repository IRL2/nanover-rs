use crate::parsers::{errors::ReadError, read_cif, read_pdb, MolecularSystem};
use quick_xml::events::{BytesEnd, BytesStart, Event};
use quick_xml::Writer;
use std::ffi::CString;
use std::io::{BufReader, Cursor};
use std::str;
use thiserror::Error;

#[derive(Debug)]
pub enum XMLTarget {
    System,
    Integrator,
}

impl TryFrom<&[u8]> for XMLTarget {
    type Error = ();

    fn try_from(value: &[u8]) -> Result<Self, Self::Error> {
        match value {
            b"System" => Ok(Self::System),
            b"Integrator" => Ok(Self::Integrator),
            _ => Err(()),
        }
    }
}

#[derive(Debug)]
pub enum ReadState {
    Unstarted,
    Ignore,
    CopyXML(XMLTarget),
    CopyStructure,
}

enum StructureType {
    None,
    Pdb,
    Pdbx,
}

pub(crate) struct PreSimulation {
    structure: Vec<u8>,
    system: Writer<Cursor<Vec<u8>>>,
    integrator: Writer<Cursor<Vec<u8>>>,
    structure_type: StructureType,
}

impl PreSimulation {
    pub fn new() -> Self {
        PreSimulation {
            structure: vec![],
            system: Writer::new(Cursor::new(Vec::new())),
            integrator: Writer::new(Cursor::new(Vec::new())),
            structure_type: StructureType::None,
        }
    }

    pub fn set_structure_type(&mut self, structure_type: &[u8]) -> Result<(), XMLParsingError> {
        self.structure_type = match structure_type {
            b"pdb" => StructureType::Pdb,
            b"pdbx" => StructureType::Pdbx,
            // Cannot happen because of the condition above
            _ => return Err(XMLParsingError::UnknownStructureType),
        };
        Ok(())
    }

    pub fn add_structure_line(&mut self, line: &[u8]) {
        self.structure.extend(line.iter());
    }

    pub fn start_xml_element(&mut self, target: &XMLTarget, element: &BytesStart) {
        let name: &str = std::str::from_utf8(element.name().into_inner()).unwrap();
        let mut elem = BytesStart::from_content(name, name.len());
        elem.extend_attributes(element.attributes().map(|attr| attr.unwrap()));
        let writer = self.choose_writer(target);
        writer.write_event(Event::Start(elem)).unwrap();
    }

    pub fn end_xml_element(&mut self, target: &XMLTarget, element: &BytesEnd) {
        let name: &str = std::str::from_utf8(element.name().into_inner()).unwrap();
        let writer = self.choose_writer(target);
        let elem = BytesEnd::new(name);
        writer.write_event(Event::End(elem)).unwrap();
    }

    pub fn empty_xml_element(&mut self, target: &XMLTarget, element: &BytesStart) {
        let name: &str = std::str::from_utf8(element.name().into_inner()).unwrap();
        let writer = self.choose_writer(target);
        let mut elem = BytesStart::from_content(name, name.len());
        elem.extend_attributes(element.attributes().map(|attr| attr.unwrap()));
        writer.write_event(Event::Empty(elem)).unwrap();
    }

    pub fn close_xml(&mut self, target: &XMLTarget) {
        let writer = self.choose_writer(target);
        writer.write_event(Event::Eof).unwrap();
    }

    pub fn finish(self) -> Result<(CString, CString, MolecularSystem), XMLParsingError> {
        let system_content =
            CString::new(str::from_utf8(&self.system.into_inner().into_inner()).unwrap()).unwrap();
        let integrator_content =
            CString::new(str::from_utf8(&self.integrator.into_inner().into_inner()).unwrap())
                .unwrap();
        let structure = match self.structure_type {
            StructureType::None => return Err(XMLParsingError::NoStructureFound),
            StructureType::Pdb => {
                let input = BufReader::new(Cursor::new(self.structure));
                read_pdb(input).map_err(XMLParsingError::PDBReadError)?
            }
            StructureType::Pdbx => {
                let input = BufReader::new(Cursor::new(self.structure));
                read_cif(input).map_err(XMLParsingError::PDBxReadError)?
            }
        };
        Ok((system_content, integrator_content, structure))
    }

    fn choose_writer(&mut self, target: &XMLTarget) -> &mut Writer<Cursor<Vec<u8>>> {
        match target {
            XMLTarget::System => &mut self.system,
            XMLTarget::Integrator => &mut self.integrator,
        }
    }
}

#[derive(Error, Debug)]
pub enum XMLParsingError {
    #[error("Not reading the expected file format.")]
    UnexpectedFileFormat,
    #[error("Error at position {1}: {0:?}.")]
    XMLError(quick_xml::Error, usize),
    #[error("Unexpected event at position {1}: {0}.")]
    LogicError(String, usize),
    #[error("The number of particles in the structure ({0}) does not match the system {1}.")]
    UnexpectedNumberOfParticles(usize, usize),
    #[error("Unrecognised structure type.")]
    UnknownStructureType,
    #[error("No structure found in the file.")]
    NoStructureFound,
    #[error("Error while reading the embedded PDB structure: {0:?}")]
    PDBReadError(ReadError),
    #[error("Error while reading the embedded PDBx structure: {0:?}.")]
    PDBxReadError(ReadError),
}
