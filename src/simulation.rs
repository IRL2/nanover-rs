extern crate openmm_sys;
extern crate pdbtbx;

use std::fs::File;
use std::io::prelude::*;
use std::str;
use std::io::{BufReader, Read, Cursor};
use std::ffi::CString;
use openmm_sys::{
    OpenMM_Force,
    OpenMM_System,
    OpenMM_System_create, OpenMM_System_destroy,
    OpenMM_System_addForce,
    OpenMM_System_addParticle,
    OpenMM_NonbondedForce_create,
    OpenMM_NonbondedForce_addParticle,
    OpenMM_Vec3,
    OpenMM_Vec3Array,
    OpenMM_Vec3Array_create, OpenMM_Vec3Array_destroy,
    OpenMM_Vec3Array_set,
    OpenMM_Integrator, OpenMM_Integrator_destroy,
    OpenMM_Integrator_step,
    OpenMM_VerletIntegrator_create,
    OpenMM_Context,
    OpenMM_Context_create,
    OpenMM_Context_setPositions,
    OpenMM_Context_destroy,
    OpenMM_Context_getState,
    OpenMM_State_destroy,
    OpenMM_State_DataType_OpenMM_State_Positions,
    OpenMM_State_getPositions,
    OpenMM_Vec3Array_getSize,
    OpenMM_Vec3Array_get,
    OpenMM_Vec3_scale,
    OpenMM_XmlSerializer_deserializeSystem,
    OpenMM_XmlSerializer_deserializeIntegrator,
    OpenMM_System_getNumParticles,
};
use quick_xml::Writer;
use quick_xml::events::{Event, BytesEnd, BytesStart};
use quick_xml::Reader;

use crate::frame::FrameData;

pub trait Simulation {
    fn step(&mut self, steps: i32);
}

pub trait ToFrameData {
    fn to_framedata(&self) -> FrameData;
}

pub struct TestSimulation {
    system: *mut OpenMM_System,
    init_pos: *mut OpenMM_Vec3Array,
    integrator: *mut OpenMM_Integrator,
    context: *mut OpenMM_Context,
}

impl TestSimulation {
    pub fn new() -> Self {
        let sim = unsafe {
            let system = OpenMM_System_create();
            let nonbond = OpenMM_NonbondedForce_create();
            OpenMM_System_addForce(system, nonbond as *mut OpenMM_Force);
            let init_pos = OpenMM_Vec3Array_create(3);

            for a in 0..3 {
                let pos = OpenMM_Vec3 {
                    x: 0.95 * a as f64,
                    y: 0.95 * a as f64,
                    z: 0.95 * a as f64,
                };

                OpenMM_Vec3Array_set(init_pos, a, pos);
                OpenMM_System_addParticle(system, 39.95);
                OpenMM_NonbondedForce_addParticle(nonbond, 0.0, 0.3350, 0.996);
            }

            let integrator = OpenMM_VerletIntegrator_create(0.004) as *mut OpenMM_Integrator;
            let context = OpenMM_Context_create(system, integrator);
            OpenMM_Context_setPositions(context, init_pos);
            Self {system, init_pos, integrator, context}
        };
        sim
    }
}

impl Drop for TestSimulation {
    fn drop(&mut self) {
        println!("Frop!");
        unsafe {
            OpenMM_Vec3Array_destroy(self.init_pos);
            OpenMM_Context_destroy(self.context);
            OpenMM_Integrator_destroy(self.integrator);
            OpenMM_System_destroy(self.system);
        }
    }
}

impl Simulation for TestSimulation {
    fn step(&mut self, steps: i32) {
        unsafe {
            OpenMM_Integrator_step(self.integrator, steps);
        }
    }
}

impl ToFrameData for TestSimulation {
    fn to_framedata(&self) -> FrameData {
        let mut positions = Vec::<f32>::new();
        unsafe {
            let state = OpenMM_Context_getState(
                self.context,
                OpenMM_State_DataType_OpenMM_State_Positions as i32,
                0,
            );
            let pos_state = OpenMM_State_getPositions(state);
            let particle_count = OpenMM_Vec3Array_getSize(pos_state);
            for i in 0..particle_count {
                let pos = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(pos_state, i), 1.0);
                positions.push(pos.x as f32);
                positions.push(pos.y as f32);
                positions.push(pos.z as f32);
            }
            OpenMM_State_destroy(state);
        }
        let mut frame = FrameData::empty();
        frame.insert_number_value("particle.count", (positions.len() / 3) as f64).unwrap();
        frame.insert_float_array("particle.positions", positions).unwrap();
        frame
    }
}

enum ReadState {
    Ignore,
    CopyXML,
    CopyStructure,
}

enum StructureType {
    None,
    PDB,
    PDBx,
}

pub struct XMLSimulation {
    system: *mut OpenMM_System,
    init_pos: *mut OpenMM_Vec3Array,
    integrator: *mut OpenMM_Integrator,
    context: *mut OpenMM_Context,
}

impl XMLSimulation {
    pub fn new<R: Read>(input: BufReader<R>) -> Self {
        let mut reader = Reader::from_reader(input);
        reader.trim_text(true);
        let mut buf = Vec::new();

        let mut read_state = ReadState::Ignore;
        let mut writer = Writer::new(Cursor::new(Vec::new()));

        let mut structure_type = StructureType::None;
        let mut structure_buffer: Vec<u8> = Vec::new();

        let mut system_content: Option<Vec<u8>> = None;
        let mut integrator_content: Option<Vec<u8>> = None;

        loop {
            match reader.read_event(&mut buf) {
                // Structure
                Ok(Event::Start(ref e)) if e.name() == b"pdbx" || e.name() == b"pdb" => {
                    let name = e.name();
                    let name_str = str::from_utf8(name).unwrap();
                    println!("Start {name_str}");
                    read_state = ReadState::CopyStructure;
                    structure_type = match name {
                        b"pdb" => StructureType::PDB,
                        b"pdbx" => StructureType::PDBx,
                        // Cannot happen because of the condition above
                        _ => panic!("Unrecognised structure type."),
                    };
                },
                Ok(Event::End(ref e)) if e.name() == b"pdbx" || e.name() == b"pdb" => {
                    let name = e.name();
                    let name_str = str::from_utf8(name).unwrap();
                    println!("End {name_str}");
                    read_state = ReadState::Ignore;
                },
                Ok(Event::Text(ref e)) => {
                    if let ReadState::CopyStructure = read_state {
                        structure_buffer.extend(e.iter());
                    }
                }

                // Structure and Integrator XML
                Ok(Event::Start(ref e)) => {
                    let name = e.name();
                    let name_str = str::from_utf8(name).unwrap();
                    if name == b"System" || name == b"Integrator" {
                        println!("Start {name_str}");
                        writer = Writer::new(Cursor::new(Vec::new()));
                        read_state = ReadState::CopyXML;
                    }
                    if let ReadState::CopyXML = read_state {
                        let mut elem = BytesStart::owned(name.to_vec(), name.len());
                        elem.extend_attributes(e.attributes().map(|attr| attr.unwrap()));
                        writer.write_event(Event::Start(elem)).unwrap();
                    }
                },
                Ok(Event::End(ref e)) => {
                    let name = e.name();
                    let name_str = str::from_utf8(name).unwrap();
                    if let ReadState::CopyXML = read_state {
                        let elem = BytesEnd::owned(name.to_vec());
                        writer.write_event(Event::End(elem)).unwrap();
                    }
                    if name == b"System" || name == b"Integrator" {
                        println!("End {name_str}");
                        read_state = ReadState::Ignore;

                        writer.write_event(Event::Eof).unwrap();

                        let to_write = writer;
                        writer = Writer::new(Cursor::new(Vec::new()));
                        if name == b"System" {
                            system_content = Some(to_write.into_inner().into_inner());
                        }
                        else if name == b"Integrator" {
                            integrator_content = Some(to_write.into_inner().into_inner());
                        }
                    }
                },
                Ok(Event::Empty(ref e)) => {
                    let name = e.name();
                    let name_str = str::from_utf8(name).unwrap();
                    if name == b"System" || name == b"Integrator" {
                        println!("Empty {name_str}");
                        read_state = ReadState::CopyXML;
                    }
                    if let ReadState::CopyXML = read_state {
                        let mut elem = BytesStart::owned(name.to_vec(), name.len());
                        elem.extend_attributes(e.attributes().map(|attr| attr.unwrap()));
                        writer.write_event(Event::Empty(elem)).unwrap();
                    }
                    if name == b"System" || name == b"Integrator" {
                        let to_write = writer;
                        writer = Writer::new(Cursor::new(Vec::new()));
                        if name == b"System" {
                            system_content = Some(to_write.into_inner().into_inner());
                        }
                        else if name == b"Integrator" {
                            integrator_content = Some(to_write.into_inner().into_inner());
                        }
                    }
                }

                // End
                Ok(Event::Eof) => {break},
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => {},
            }
        }

        let system_content = match system_content {
            None => panic!("No system found in the file."),
            Some(content) => CString::new(str::from_utf8(&content).unwrap()).unwrap(),
        };
        let integrator_content = match integrator_content {
            None => panic!("No integrator found in the file."),
            Some(content) => CString::new(str::from_utf8(&content).unwrap()).unwrap(),
        };
        let structure = match structure_type {
            StructureType::None => panic!("No structure found."),
            StructureType::PDB => {
                let (structure, _) = pdbtbx::open_pdb_raw(
                    BufReader::new(Cursor::new(structure_buffer)),
                    pdbtbx::Context::None,
                    pdbtbx::StrictnessLevel::Loose
                ).unwrap();
                structure
            },
            StructureType::PDBx => {
                let (structure, _) = pdbtbx::open_mmcif_raw(
                    str::from_utf8(&structure_buffer).unwrap(),
                    pdbtbx::StrictnessLevel::Loose
                ).unwrap();
                structure
            },
        };
        let n_atoms = structure.atom_count();
        println!("Particles in the structure: {n_atoms}");

        let sim = unsafe {
            println!("Entering the unsafe section");
            let init_pos = OpenMM_Vec3Array_create(n_atoms.try_into().unwrap());
            for (i, atom) in structure.atoms().enumerate() {
                // The input structure is in angstsroms, but we need to provide nm.
                let position = OpenMM_Vec3 {
                    x: atom.x() / 10.0,
                    y: atom.y() / 10.0,
                    z: atom.z() / 10.0,
                };
                OpenMM_Vec3Array_set(init_pos, i.try_into().unwrap(), position);
            }
            println!("Coordinate array built");

            //println!("{}", system_content.to_str().unwrap());
            {
                let mut file = File::create("foo.xml").unwrap();
                file.write_all(system_content.to_bytes()).unwrap();
            }
            let system = OpenMM_XmlSerializer_deserializeSystem(system_content.as_ptr());
            println!("System read");
            println!("Paticles in system: {}", OpenMM_System_getNumParticles(system));
            let integrator = OpenMM_XmlSerializer_deserializeIntegrator(integrator_content.as_ptr());
            println!("Integrator read");
            let context = OpenMM_Context_create(system, integrator);
            OpenMM_Context_setPositions(context, init_pos);
            Self {system, init_pos, integrator, context}
        };
        sim
    }
}

impl Drop for XMLSimulation {
    fn drop(&mut self) {
        println!("Frop!");
        unsafe {
            OpenMM_Vec3Array_destroy(self.init_pos);
            OpenMM_Context_destroy(self.context);
            OpenMM_Integrator_destroy(self.integrator);
            OpenMM_System_destroy(self.system);
        }
    }
}

impl Simulation for XMLSimulation {
    fn step(&mut self, steps: i32) {
        unsafe {
            OpenMM_Integrator_step(self.integrator, steps);
        }
    }
}

impl ToFrameData for XMLSimulation {
    fn to_framedata(&self) -> FrameData {
        let mut positions = Vec::<f32>::new();
        unsafe {
            let state = OpenMM_Context_getState(
                self.context,
                OpenMM_State_DataType_OpenMM_State_Positions as i32,
                0,
            );
            let pos_state = OpenMM_State_getPositions(state);
            let particle_count = OpenMM_Vec3Array_getSize(pos_state);
            for i in 0..particle_count {
                let pos = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(pos_state, i), 1.0);
                positions.push(pos.x as f32);
                positions.push(pos.y as f32);
                positions.push(pos.z as f32);
            }
            OpenMM_State_destroy(state);
        }
        let mut frame = FrameData::empty();
        frame.insert_number_value("particle.count", (positions.len() / 3) as f64).unwrap();
        frame.insert_float_array("particle.positions", positions).unwrap();
        frame
    }
}
