extern crate openmm_sys;
extern crate pdbtbx;

use std::fs::File;
use std::io::prelude::*;
use std::str;
use std::io::{BufReader, Read, Cursor};
use std::ffi::{CString, CStr};
use std::collections::{BTreeMap, HashSet};
use openmm_sys::{
    OpenMM_System,
    OpenMM_System_addForce,
    OpenMM_System_destroy,
    OpenMM_Vec3,
    OpenMM_Vec3Array,
    OpenMM_Vec3Array_create, OpenMM_Vec3Array_destroy,
    OpenMM_Vec3Array_set,
    OpenMM_Integrator, OpenMM_Integrator_destroy,
    OpenMM_Integrator_step,
    OpenMM_Context,
    OpenMM_Context_create,
    OpenMM_Context_setPositions,
    OpenMM_Context_destroy,
    OpenMM_Context_getState,
    OpenMM_State_destroy,
    OpenMM_State_DataType_OpenMM_State_Positions,
    OpenMM_State_getPositions,
    OpenMM_State_getPeriodicBoxVectors,
    OpenMM_Vec3Array_getSize,
    OpenMM_Vec3Array_get,
    OpenMM_Vec3_scale,
    OpenMM_XmlSerializer_deserializeSystem,
    OpenMM_XmlSerializer_deserializeIntegrator,
    OpenMM_System_getNumParticles,
    OpenMM_System_getParticleMass,
    OpenMM_Context_getPlatform,
    OpenMM_Platform_getName,
    OpenMM_Force,
    OpenMM_CustomExternalForce,
    OpenMM_CustomExternalForce_create,
    OpenMM_CustomExternalForce_destroy,
    OpenMM_CustomExternalForce_addPerParticleParameter,
    OpenMM_CustomExternalForce_addParticle,
    OpenMM_CustomExternalForce_setParticleParameters,
    OpenMM_DoubleArray_create,
    OpenMM_DoubleArray_destroy,
    OpenMM_DoubleArray_set,
};
use quick_xml::Writer;
use quick_xml::events::{Event, BytesEnd, BytesStart};
use quick_xml::Reader;

use crate::frame::FrameData;

#[derive(Debug)]
pub enum InteractionKind {
    GAUSSIAN,
    HARMONIC,
}

#[derive(Debug)]
pub struct IMDInteraction {
    position: [f64; 3],
    particles: Vec<usize>,
    pub kind: InteractionKind,
    max_force: Option<f64>,
}

impl IMDInteraction {
    pub fn new(position: [f64; 3], particles: Vec<usize>, kind: InteractionKind, max_force: Option<f64>) -> Self {
        Self {position, particles, kind, max_force}
    }
}

pub struct InteractionForce {
    pub selection: usize,
    pub force: [f64; 3],
}

pub struct Interaction {
    pub forces: Vec<InteractionForce>,
}

pub trait Simulation {
    fn step(&mut self, steps: i32);
}

pub trait ToFrameData {
    fn to_framedata(&self) -> FrameData;
    fn to_topology_framedata(&self) -> FrameData;
}

pub trait IMD {
    fn update_imd_forces(&mut self, interactions: Vec<Interaction>) -> Result<(), ()>;
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
    topology: pdbtbx::PDB,
    platform_name: String,
    imd_force: *mut OpenMM_CustomExternalForce,
    n_particles: usize,
    previous_particle_touched: HashSet<i32>,
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

            {
                let mut file = File::create("foo.xml").unwrap();
                file.write_all(system_content.to_bytes()).unwrap();
            }
            let system = OpenMM_XmlSerializer_deserializeSystem(system_content.as_ptr());
            println!("System read");
            let n_particles = OpenMM_System_getNumParticles(system);
            let imd_force = Self::add_imd_force(n_particles);
            let cast_as_force = imd_force as *mut OpenMM_Force;
            OpenMM_System_addForce(system, cast_as_force);
            let n_particles = n_particles as usize;
            if n_particles != n_atoms {
                panic!("The number of particles in the structure does not match the system.");
            }
            println!("Particles in system: {}", n_particles);
            let integrator = OpenMM_XmlSerializer_deserializeIntegrator(integrator_content.as_ptr());
            println!("Integrator read");
            let context = OpenMM_Context_create(system, integrator);
            OpenMM_Context_setPositions(context, init_pos);

            let platform = OpenMM_Context_getPlatform(context);
            let platform_name = CStr::from_ptr(OpenMM_Platform_getName(platform))
                .to_str().unwrap()
                .to_string();
            
            let previous_particle_touched = HashSet::new();

            Self {
                system,
                init_pos,
                integrator,
                context,
                topology: structure,
                platform_name,
                imd_force,
                n_particles,
                previous_particle_touched,
            }
        };
        sim
    }
    
    pub fn get_platform_name(&self) -> String {
        self.platform_name.clone()
    }

    unsafe fn add_imd_force(n_particles: i32) -> *mut OpenMM_CustomExternalForce {
        let energy_expression = CString::new("-fx * x - fy * y - fz * z").unwrap();
        let force = OpenMM_CustomExternalForce_create(energy_expression.into_raw() as *const i8);
        OpenMM_CustomExternalForce_addPerParticleParameter(
            force, CString::new("fx").unwrap().into_raw());
        OpenMM_CustomExternalForce_addPerParticleParameter(
            force, CString::new("fy").unwrap().into_raw());
        OpenMM_CustomExternalForce_addPerParticleParameter(
            force, CString::new("fz").unwrap().into_raw());
        for i in 0..n_particles {
            let zeros = OpenMM_DoubleArray_create(3);
            OpenMM_DoubleArray_set(zeros, 0, 0.0);
            OpenMM_DoubleArray_set(zeros, 1, 0.0);
            OpenMM_DoubleArray_set(zeros, 2, 0.0);
            OpenMM_CustomExternalForce_addParticle(force, i, zeros);
        }
        force
    }

    pub fn compute_forces(&self, imd_interaction: &Vec<IMDInteraction>) -> Vec<Interaction> {
        unsafe {
            let state = OpenMM_Context_getState(
                self.context,
                OpenMM_State_DataType_OpenMM_State_Positions as i32,
                0,
            );
            let pos_state = OpenMM_State_getPositions(state);
            let interactions = imd_interaction.iter().map(|imd| {
                let selection: Vec<i32> = imd.particles.iter()
                    // *Ignore* particle indices that are out of bound.
                    .filter(|p| **p < self.n_particles)
                    // Convert particle indices to i32 so we can use them
                    // in OpenMM methods. *Ignore* indices that do not fit.
                    .flat_map(|p| (*p).try_into())
                    .collect();
                let (com_sum, total_mass) = selection.iter()
                    // Retrieve the position and mass of the selected particles.
                    .map(|p| {
                        let position = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(pos_state, *p), 1.0);
                        let mass = OpenMM_System_getParticleMass(self.system, *p);
                        ([position.x, position.y, position.z], mass)
                    })
                    .fold(([0.0, 0.0, 0.0], 0.0), |acc, x| {
                        (
                            [
                                acc.0[0] + x.0[0] * x.1,
                                acc.0[1] + x.0[1] * x.1,
                                acc.0[2] + x.0[2] * x.1,
                            ],
                            acc.1 + x.1
                        )
                    });
                let com = [
                    com_sum[0] / total_mass,
                    com_sum[1] / total_mass,
                    com_sum[2] / total_mass,
                ];
                let c = imd.position;
                let diff = [
                    com[0] - c[0],
                    com[1] - c[1],
                    com[2] - c[2], 
                ];
                let sigma_sqr = 1.0;  // For now we use this as a constant. It is in the python version.
                let distance_sqr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                let gauss = (-distance_sqr / (2.0 * sigma_sqr)).exp();
                let com_force = [
                    -(diff[0] / sigma_sqr) * gauss,
                    -(diff[1] / sigma_sqr) * gauss,
                    -(diff[2] / sigma_sqr) * gauss,
                ];
                let force_per_particle = [
                    com_force[0] / self.n_particles as f64,
                    com_force[1] / self.n_particles as f64,
                    com_force[2] / self.n_particles as f64,
                ];
                Interaction {forces: selection.iter().map(|p| InteractionForce {selection: *p as usize, force: force_per_particle}).collect()}
            })
            .collect();

            OpenMM_State_destroy(state);

            interactions
        }
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
            OpenMM_CustomExternalForce_destroy(self.imd_force);
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
        let mut box_vectors = Vec::<f32>::new();
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

            let mut a = OpenMM_Vec3{x: 0.0, y: 0.0, z: 0.0};
            let mut b = OpenMM_Vec3{x: 0.0, y: 0.0, z: 0.0};
            let mut c = OpenMM_Vec3{x: 0.0, y: 0.0, z: 0.0};
            OpenMM_State_getPeriodicBoxVectors(state, &mut a, &mut b, &mut c);
            box_vectors.push(a.x as f32);
            box_vectors.push(a.y as f32);
            box_vectors.push(a.z as f32);
            box_vectors.push(b.x as f32);
            box_vectors.push(b.y as f32);
            box_vectors.push(b.z as f32);
            box_vectors.push(c.x as f32);
            box_vectors.push(c.y as f32);
            box_vectors.push(c.z as f32);

            OpenMM_State_destroy(state);
        }
        let mut frame = FrameData::empty();
        frame.insert_number_value("particle.count", (positions.len() / 3) as f64).unwrap();
        frame.insert_float_array("particle.positions", positions).unwrap();
        frame.insert_float_array("system.box.vectors", box_vectors).unwrap();

        frame
    }

    fn to_topology_framedata(&self) -> FrameData {
        let mut frame = FrameData::empty();

        let n_particles = self.topology.atom_count();
        frame.insert_number_value("particle.count", (n_particles) as f64).unwrap();

        let elements: Vec<u32> = self.topology.atoms()
            .map(|atom| {atom.atomic_number().unwrap_or(0).try_into().unwrap()}).collect();
        frame.insert_index_array("particle.elements", elements).unwrap();

        frame
    }
}

impl IMD for XMLSimulation {
    fn update_imd_forces(&mut self, interactions: Vec<Interaction>) -> Result<(), ()> {
        let mut forces = zeroed_out(&self.previous_particle_touched);
        let accumulated_forces = accumulate_forces(interactions);
        self.previous_particle_touched = HashSet::new();
        accumulated_forces.iter().for_each(|kv| {self.previous_particle_touched.insert(*kv.0);});

        forces.extend(accumulated_forces);
        for (index, force) in forces {
            if index as usize >= self.n_particles {
                return Err(());
            }
            unsafe {
                let force_array = OpenMM_DoubleArray_create(3);
                OpenMM_DoubleArray_set(force_array, 0, force[0]);
                OpenMM_DoubleArray_set(force_array, 1, force[1]);
                OpenMM_DoubleArray_set(force_array, 2, force[2]);
                OpenMM_CustomExternalForce_setParticleParameters(
                    self.imd_force, index, index, force_array);
            }
        }
        Ok(())
    }
}

fn accumulate_forces(interactions: Vec<Interaction>) -> BTreeMap<i32, [f64; 3]> {
    let mut btree: BTreeMap<i32, [f64;3]> = BTreeMap::new();
    for interaction in interactions {
        for particle in interaction.forces {
            let index: i32 = particle.selection
                .try_into()
                .expect("Particle index does not fit an i32.");
            btree.entry(index)
                .and_modify(|f| {
                    f[0] += particle.force[0];
                    f[1] += particle.force[1];
                    f[2] += particle.force[2];
                })
                .or_insert(particle.force);
        }
    }
    btree
}

fn zeroed_out(indices: &HashSet<i32>) -> BTreeMap<i32, [f64; 3]> {
    let mut btree: BTreeMap<i32, [f64; 3]> = BTreeMap::new();
    indices.iter().for_each(|i| {btree.insert(*i, [0.0, 0.0, 0.0]);});
    btree
}