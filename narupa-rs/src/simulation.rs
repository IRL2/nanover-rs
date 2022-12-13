extern crate openmm_sys;

use openmm_sys::{
    OpenMM_Context, OpenMM_Context_create, OpenMM_Context_destroy, OpenMM_Context_getPlatform,
    OpenMM_Context_getState, OpenMM_Context_setPositions, OpenMM_CustomExternalForce,
    OpenMM_CustomExternalForce_addParticle, OpenMM_CustomExternalForce_addPerParticleParameter,
    OpenMM_CustomExternalForce_create, OpenMM_CustomExternalForce_destroy,
    OpenMM_CustomExternalForce_setParticleParameters,
    OpenMM_CustomExternalForce_updateParametersInContext, OpenMM_DoubleArray_create,
    OpenMM_DoubleArray_set, OpenMM_Force, OpenMM_Integrator, OpenMM_Integrator_destroy,
    OpenMM_Integrator_step, OpenMM_Platform_getName, OpenMM_Platform_getNumPlatforms,
    OpenMM_Platform_loadPluginsFromDirectory, OpenMM_State_DataType_OpenMM_State_Positions,
    OpenMM_State_destroy, OpenMM_State_getPeriodicBoxVectors, OpenMM_State_getPositions,
    OpenMM_System, OpenMM_System_addForce, OpenMM_System_destroy, OpenMM_System_getNumParticles,
    OpenMM_System_getParticleMass, OpenMM_Vec3, OpenMM_Vec3Array, OpenMM_Vec3Array_create,
    OpenMM_Vec3Array_destroy, OpenMM_Vec3Array_get, OpenMM_Vec3Array_getSize, OpenMM_Vec3Array_set,
    OpenMM_Vec3_scale, OpenMM_XmlSerializer_deserializeIntegrator,
    OpenMM_XmlSerializer_deserializeSystem, OpenMM_State, OpenMM_Context_setState,
};
use quick_xml::events::{BytesEnd, BytesStart, Event};
use quick_xml::Reader;
use quick_xml::Writer;
use std::collections::{BTreeMap, HashSet};
use std::env;
use std::ffi::{CStr, CString};
use std::io::{BufReader, Cursor, Read};
use std::str;

use crate::frame::FrameData;
use crate::parsers::{read_pdb, read_cif, MolecularSystem};

type Coordinate = [f64; 3];
type CoordMap = BTreeMap<i32, Coordinate>;

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
    scale: f64,
}

impl IMDInteraction {
    pub fn new(
        position: [f64; 3],
        particles: Vec<usize>,
        kind: InteractionKind,
        max_force: Option<f64>,
        scale: f64,
    ) -> Self {
        Self {
            position,
            particles,
            kind,
            max_force,
            scale,
        }
    }
}

pub struct InteractionForce {
    pub selection: usize,
    pub force: Coordinate,
}

pub struct Interaction {
    pub forces: Vec<InteractionForce>,
}

pub trait Simulation {
    fn step(&mut self, steps: i32);
    fn reset(&mut self);
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
    Pdb,
    Pdbx,
}

pub struct XMLSimulation {
    system: *mut OpenMM_System,
    init_pos: *mut OpenMM_Vec3Array,
    integrator: *mut OpenMM_Integrator,
    context: *mut OpenMM_Context,
    topology: MolecularSystem,
    platform_name: String,
    imd_force: *mut OpenMM_CustomExternalForce,
    n_particles: usize,
    previous_particle_touched: HashSet<i32>,
    initial_state: *mut OpenMM_State,
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
                        b"pdb" => StructureType::Pdb,
                        b"pdbx" => StructureType::Pdbx,
                        // Cannot happen because of the condition above
                        _ => panic!("Unrecognised structure type."),
                    };
                }
                Ok(Event::End(ref e)) if e.name() == b"pdbx" || e.name() == b"pdb" => {
                    let name = e.name();
                    let name_str = str::from_utf8(name).unwrap();
                    println!("End {name_str}");
                    read_state = ReadState::Ignore;
                }
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
                }
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
                        } else if name == b"Integrator" {
                            integrator_content = Some(to_write.into_inner().into_inner());
                        }
                    }
                }
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
                        } else if name == b"Integrator" {
                            integrator_content = Some(to_write.into_inner().into_inner());
                        }
                    }
                }

                // End
                Ok(Event::Eof) => break,
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => {}
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
            StructureType::Pdb => {
                let input = BufReader::new(Cursor::new(structure_buffer));
                read_pdb(input).expect("Could not read the PDB.")
            }
            StructureType::Pdbx => {
                let input = BufReader::new(Cursor::new(structure_buffer));
                read_cif(input).expect("Could not read the PDBx.")
            }
        };
        let n_atoms = structure.atom_count();
        println!("Particles in the structure: {n_atoms}");

        let sim = unsafe {
            println!("Entering the unsafe section");

            println!("Loading plugins");
            match env::var("OPENMM_PLUGIN_DIR") {
                Ok(dirname) => {
                    let lib_directory = CString::new(dirname).unwrap();
                    OpenMM_Platform_loadPluginsFromDirectory(lib_directory.into_raw());
                }
                Err(_) => {
                    println!("No plugin to load, set OPENMM_PLUGIN_DIR");
                }
            }

            let n_platform = OpenMM_Platform_getNumPlatforms();
            println!("Number of platforms registered: {n_platform}");

            let init_pos = OpenMM_Vec3Array_create(n_atoms.try_into().unwrap());
            for (i, atom) in structure.positions.iter().enumerate() {
                // The input structure is in angstsroms, but we need to provide nm.
                let position = OpenMM_Vec3 {
                    x: atom[0],
                    y: atom[1],
                    z: atom[2],
                };
                OpenMM_Vec3Array_set(init_pos, i.try_into().unwrap(), position);
            }
            println!("Coordinate array built");

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
            let integrator =
                OpenMM_XmlSerializer_deserializeIntegrator(integrator_content.as_ptr());
            println!("Integrator read");
            let context = OpenMM_Context_create(system, integrator);
            OpenMM_Context_setPositions(context, init_pos);

            let initial_state = OpenMM_Context_getState(
                context, OpenMM_State_DataType_OpenMM_State_Positions as i32, 0);

            let platform = OpenMM_Context_getPlatform(context);
            let platform_name = CStr::from_ptr(OpenMM_Platform_getName(platform))
                .to_str()
                .unwrap()
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
                initial_state,
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
            force,
            CString::new("fx").unwrap().into_raw(),
        );
        OpenMM_CustomExternalForce_addPerParticleParameter(
            force,
            CString::new("fy").unwrap().into_raw(),
        );
        OpenMM_CustomExternalForce_addPerParticleParameter(
            force,
            CString::new("fz").unwrap().into_raw(),
        );
        for i in 0..n_particles {
            let zeros = OpenMM_DoubleArray_create(3);
            OpenMM_DoubleArray_set(zeros, 0, 0.0);
            OpenMM_DoubleArray_set(zeros, 1, 0.0);
            OpenMM_DoubleArray_set(zeros, 2, 0.0);
            OpenMM_CustomExternalForce_addParticle(force, i, zeros);
        }
        force
    }

    pub fn compute_forces(&self, imd_interaction: &[IMDInteraction]) -> Vec<Interaction> {
        unsafe {
            let state = OpenMM_Context_getState(
                self.context,
                OpenMM_State_DataType_OpenMM_State_Positions as i32,
                0,
            );
            let pos_state = OpenMM_State_getPositions(state);
            let interactions = imd_interaction
                .iter()
                .map(|imd| {
                    let max_force = imd.max_force.unwrap_or(f64::INFINITY);
                    let selection = filter_selection(&imd.particles, self.n_particles);
                    let masses = get_selection_masses_from_system(&selection, self.system);
                    let particle_positions = get_selection_positions_from_state_positions(&selection, pos_state);
                    let com = compute_com(&particle_positions, &masses);
                    let interaction_position = imd.position;
                    let diff = [
                        com[0] - interaction_position[0],
                        com[1] - interaction_position[1],
                        com[2] - interaction_position[2],
                    ];
                    let sigma = 1.0; // For now we use this as a constant. It is in the python version.
                    let com_force = match imd.kind {
                        InteractionKind::GAUSSIAN => compute_gaussian_force(diff, sigma),
                        InteractionKind::HARMONIC => compute_harmonic_force(diff, sigma),
                    };
                    build_interaction(
                        &com_force,
                        self.n_particles,
                        imd.scale,
                        &selection,
                        &masses,
                        max_force,
                    )
                })
                .collect();

            OpenMM_State_destroy(state);

            interactions
        }
    }

    pub fn reset_state(&mut self) {
        unsafe {
            OpenMM_Context_setState(self.context, self.initial_state);
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
            OpenMM_State_destroy(self.initial_state);
        }
    }
}

impl Simulation for XMLSimulation {
    fn step(&mut self, steps: i32) {
        unsafe {
            OpenMM_Integrator_step(self.integrator, steps);
        }
    }

    fn reset(&mut self) {
        self.reset_state();
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

            let mut a = OpenMM_Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            };
            let mut b = OpenMM_Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            };
            let mut c = OpenMM_Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            };
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
        frame
            .insert_number_value("particle.count", (positions.len() / 3) as f64)
            .unwrap();
        frame
            .insert_float_array("particle.positions", positions)
            .unwrap();
        frame
            .insert_float_array("system.box.vectors", box_vectors)
            .unwrap();

        frame
    }

    fn to_topology_framedata(&self) -> FrameData {
        let mut frame = FrameData::empty();

        let n_particles = self.topology.atom_count();
        frame
            .insert_number_value("particle.count", (n_particles) as f64)
            .unwrap();

        let elements: Vec<u32> = self
            .topology
            .elements
            .iter()
            .map(|e| e.unwrap_or(0) as u32)
            .collect();
        frame
            .insert_index_array("particle.elements", elements)
            .unwrap();
        
        let atom_resindex: Vec<u32> = self.topology.atom_resindex
            .iter()
            .map(|e| *e as u32)
            .collect();
        frame
            .insert_number_value("residue.count", atom_resindex.len() as f64)
            .unwrap();
        frame
            .insert_index_array("particle.residues", atom_resindex)
            .unwrap();
        frame
            .insert_string_array("residue.names", self.topology.resnames.clone())
            .unwrap();
        let resids: Vec<f32> = self.topology.resids.iter().map(|e| *e as f32).collect();
        frame
            .insert_float_array("residue.ids", resids)
            .unwrap();

        let residue_chain_index: Vec<u32> = self.topology.residue_chain_index
            .iter()
            .map(|e| *e as u32)
            .collect();
        frame
            .insert_number_value("chain.count", residue_chain_index.len() as f64)
            .unwrap();
        frame
            .insert_index_array("residue.chains", residue_chain_index)
            .unwrap();
        frame
            .insert_string_array("chain.names", self.topology.chain_identifiers.clone())
            .unwrap();

        let (bonds, bond_orders): (Vec<[u32; 2]>, Vec<f32>) = self.topology
            .bonds
            .iter()
            .map(|bond| ([bond.0 as u32, bond.1 as u32], bond.2))
            .unzip();
        let bonds: Vec<u32> = bonds.into_iter().flatten().collect();
        if !bonds.is_empty() {
            frame
                .insert_index_array("bond.pairs", bonds)
                .unwrap();
            frame.insert_float_array("bond.orders", bond_orders).unwrap();
        }

        frame
    }
}

impl IMD for XMLSimulation {
    fn update_imd_forces(&mut self, interactions: Vec<Interaction>) -> Result<(), ()> {
        let mut forces = zeroed_out(&self.previous_particle_touched);
        let accumulated_forces = accumulate_forces(interactions);
        self.previous_particle_touched = HashSet::new();
        accumulated_forces.iter().for_each(|kv| {
            self.previous_particle_touched.insert(*kv.0);
        });
        forces.extend(accumulated_forces);

        for (index, force) in &forces {
            if *index as usize >= self.n_particles {
                return Err(());
            }
            unsafe {
                let force_array = OpenMM_DoubleArray_create(3);
                OpenMM_DoubleArray_set(force_array, 0, force[0]);
                OpenMM_DoubleArray_set(force_array, 1, force[1]);
                OpenMM_DoubleArray_set(force_array, 2, force[2]);
                OpenMM_CustomExternalForce_setParticleParameters(
                    self.imd_force,
                    *index,
                    *index,
                    force_array,
                );
            }
        }

        unsafe {
            OpenMM_CustomExternalForce_updateParametersInContext(self.imd_force, self.context);
        }

        Ok(())
    }
}

fn compute_com(positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
    let (com_sum, total_mass) = positions
        .iter()
        .zip(masses)
        .fold(([0.0, 0.0, 0.0], 0.0), |acc, x| {
            (
                [
                    acc.0[0] + x.0[0] * x.1,
                    acc.0[1] + x.0[1] * x.1,
                    acc.0[2] + x.0[2] * x.1,
                ],
                acc.1 + x.1,
            )
        });
    [
        com_sum[0] / total_mass,
        com_sum[1] / total_mass,
        com_sum[2] / total_mass,
    ]
}

fn accumulate_forces(interactions: Vec<Interaction>) -> CoordMap {
    let mut btree: CoordMap = BTreeMap::new();
    for interaction in interactions {
        for particle in interaction.forces {
            let index: i32 = particle
                .selection
                .try_into()
                .expect("Particle index does not fit an i32.");
            btree
                .entry(index)
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

fn build_interaction(
    com_force: &[f64; 3],
    n_particles: usize,
    scale: f64,
    selection: &[i32],
    masses: &[f64],
    max_force: f64,
) -> Interaction {
    let force_per_particle = [
        scale * com_force[0] / n_particles as f64,
        scale * com_force[1] / n_particles as f64,
        scale * com_force[2] / n_particles as f64,
    ];
    Interaction {
        forces: selection
            .iter()
            .zip(masses)
            .map(|pm| InteractionForce {
                selection: *pm.0 as usize,
                force: [
                    (force_per_particle[0] * pm.1).clamp(-max_force, max_force),
                    (force_per_particle[1] * pm.1).clamp(-max_force, max_force),
                    (force_per_particle[2] * pm.1).clamp(-max_force, max_force),
                ],
            })
            .collect(),
    }
}

fn filter_selection(particles: &[usize], n_particles: usize) -> Vec<i32> {
    let selection: HashSet<i32> = particles
        .iter()
        // *Ignore* particle indices that are out of bound.
        .filter(|p| **p < n_particles)
        // Convert particle indices to i32 so we can use them
        // in OpenMM methods. *Ignore* indices that do not fit.
        .flat_map(|p| (*p).try_into())
        // By collecting into a HashSet, we remove duplicate indices.
        .collect();
    selection.into_iter().collect()
}

/// Create a map with the provided indices and arrays of zeros as values.
fn zeroed_out(indices: &HashSet<i32>) -> CoordMap {
    let mut btree: CoordMap = BTreeMap::new();
    indices.iter().for_each(|i| {
        btree.insert(*i, [0.0, 0.0, 0.0]);
    });
    btree
}

fn compute_gaussian_force(diff: Coordinate, sigma: f64) -> Coordinate {
    let sigma_sqr = sigma * sigma;
    let distance_sqr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
    let gauss = (-distance_sqr / (2.0 * sigma_sqr)).exp();
    [
        -(diff[0] / sigma_sqr) * gauss,
        -(diff[1] / sigma_sqr) * gauss,
        -(diff[2] / sigma_sqr) * gauss,
    ]
}

fn compute_harmonic_force(diff: Coordinate, k: f64) -> Coordinate {
    let factor = -2.0 * k;
    [factor * diff[0], factor * diff[1], factor * diff[2]]
}

unsafe fn get_selection_masses_from_system(
        selection: &[i32],
        system: *mut OpenMM_System,
    ) -> Vec<f64> {
    selection
        .iter()
        .map(|p| OpenMM_System_getParticleMass(system, *p))
        .collect()
}

unsafe fn get_selection_positions_from_state_positions(
        selection: &[i32],
        pos_state: *const OpenMM_Vec3Array,
    ) -> Vec<[f64; 3]> {
    selection
        .iter()
        .map(|p| {
            let position =
                OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(pos_state, *p), 1.0);
            [position.x, position.y, position.z]
        })
        .collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::collections::HashSet;

    const EXP_1: f64 = 0.6065306597126334;  // exp(-0.5)
    const UNIT: Coordinate = [0.5773502691896258; 3]; // [1, 1, 1] / |[1, 1, 1]|

    #[test]
    fn test_zeroed_out() {
        let request_indices: [i32; 4] = [2, 4, 7, 9];
        let set_indices: HashSet<i32> = HashSet::from(request_indices);
        let zeroed = zeroed_out(&set_indices);
        let keys: Vec<_> = zeroed.keys().cloned().collect();
        let values: Vec<_> = zeroed.values().cloned().collect();
        assert_eq!(keys, request_indices);
        assert_eq!(values, [[0.0; 3]; 4]);
    }

    #[rstest]
    #[case([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-EXP_1, 0.0, 0.0])]
    #[case([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [ EXP_1, 0.0, 0.0])]
    #[case([1.0, 3.0, 0.0], [1.0, 2.0, 0.0], [0.0, -EXP_1, 0.0])]
    #[case([1.0, 3.0, 3.0], [1.0, 3.0, 2.0], [0.0, 0.0, -EXP_1])]
    #[case(UNIT, [0.0, 0.0, 0.0], [-EXP_1 * UNIT[0]; 3])]
    fn test_gaussian_force(
            #[case] position: Coordinate,
            #[case] interaction_position: Coordinate, 
            #[case] expected_force: Coordinate,
        ) {
        let sigma = 1.0;
        let diff = [
            position[0] - interaction_position[0],
            position[1] - interaction_position[1],
            position[2] - interaction_position[2],
        ];
        let force = compute_gaussian_force(diff, sigma);
        assert_f64_near!(force[0], expected_force[0]);
        assert_f64_near!(force[1], expected_force[1]);
        assert_f64_near!(force[2], expected_force[2]);
    }

    #[rstest]
    #[case([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-2.0, 0.0, 0.0])]
    #[case([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [ 2.0, 0.0, 0.0])]
    #[case([1.0, 3.0, 0.0], [1.0, 2.0, 0.0], [0.0, -2.0, 0.0])]
    #[case([1.0, 3.0, 3.0], [1.0, 3.0, 2.0], [0.0, 0.0, -2.0])]
    #[case(UNIT, [0.0, 0.0, 0.0], [-2.0 * UNIT[0]; 3])]
    #[case([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [-2.0, -2.0, -2.0])]
    #[case([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [0.0, 0.0, 0.0])]
    #[case([-1.0, -1.0, -1.0], [0.0, 0.0, 0.0], [2.0, 2.0, 2.0])]
    fn test_harmonic_force(
        #[case] position: Coordinate,
        #[case] interaction_position: Coordinate,
        #[case] expected_force: Coordinate,
    ) {
        let k = 1.0;
        let diff = [
            position[0] - interaction_position[0],
            position[1] - interaction_position[1],
            position[2] - interaction_position[2],
        ];
        let force = compute_harmonic_force(diff, k);
        assert_f64_near!(force[0], expected_force[0]);
        assert_f64_near!(force[1], expected_force[1]);
        assert_f64_near!(force[2], expected_force[2]);
    }

    #[test]
    fn test_accumulate_forces() {
        let interactions: Vec<Interaction> = vec![
            Interaction {forces: vec![
                InteractionForce {
                    selection: 43,
                    force: [1.0, 2.0, 3.0],
                },
                InteractionForce {
                    selection: 12,
                    force: [2.1, 3.2, 4.3],
                }
            ]},
            Interaction {forces: vec![
                InteractionForce {
                    selection: 2,
                    force: [3.2, 4.3, 5.4],
                }
            ]},
            Interaction {forces: vec![
                InteractionForce {
                    selection: 43,
                    force: [4.3, 5.4, 6.5],
                }
            ]}
        ];
        let accumulated = accumulate_forces(interactions);
        let mut expected = BTreeMap::new();
        expected.insert(2, [3.2, 4.3, 5.4]);
        expected.insert(12, [2.1, 3.2, 4.3]);
        expected.insert(43, [5.3, 7.4, 9.5]);

        assert_eq!(accumulated, expected);
    }

    #[test]
    fn test_filter_selection() {
        let input: Vec<usize> = vec![4, 5, 6, 11, 10, 5, 9];
        let n_particles = 10;
        let mut result = filter_selection(&input, n_particles);
        result.sort();
        assert_eq!(result, vec![4, 5, 6, 9]);
    }
}