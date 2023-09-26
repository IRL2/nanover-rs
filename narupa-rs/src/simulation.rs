extern crate openmm_sys;
use thiserror::Error;

use openmm_sys::{OpenMM_Vec3Array, OpenMM_Vec3Array_get, OpenMM_Vec3_scale};
use std::collections::{BTreeMap, HashSet};

use narupa_proto::frame::FrameData;

type Coordinate = [f64; 3];
type CoordMap = BTreeMap<i32, Coordinate>;

#[derive(Debug, Default)]
pub enum InteractionKind {
    #[default]
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
    id: Option<String>,
}

impl IMDInteraction {
    pub fn new(
        position: [f64; 3],
        particles: Vec<usize>,
        kind: InteractionKind,
        max_force: Option<f64>,
        scale: f64,
        id: Option<String>,
    ) -> Self {
        Self {
            position,
            particles,
            kind,
            max_force,
            scale,
            id,
        }
    }
}

impl Default for IMDInteraction {
    fn default() -> Self {
        Self {
            position: [0.0, 0.0, 0.0],
            particles: vec![],
            kind: InteractionKind::default(),
            max_force: Some(20000.0),
            scale: 1.0,
            id: None,
        }
    }
}

#[derive(Debug)]
pub struct InteractionForce {
    pub selection: usize,
    pub force: Coordinate,
}

#[derive(Debug)]
pub struct Interaction {
    pub forces: Vec<InteractionForce>,
    pub energy: Option<f64>,
    pub id: Option<String>,
}

pub trait Simulation {
    fn step(&mut self, steps: i32);
    fn reset(&mut self);
}

pub trait ToFrameData {
    fn to_framedata(&self) -> FrameData;
    fn to_topology_framedata(&self) -> FrameData;
}

#[derive(Debug, Error)]
#[error("Particle with index {index} out of range for {particle_count} particles.")]
pub struct ParticleOutOfRange {
    index: i32,
    particle_count: usize,
}

impl ParticleOutOfRange {
    pub fn new(index: i32, particle_count: usize) -> Self {
        Self {
            index,
            particle_count,
        }
    }
}

pub trait IMD {
    fn update_imd_forces(&mut self, interactions: &[Interaction])
        -> Result<(), ParticleOutOfRange>;
}

fn compute_com(positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
    let (com_sum, total_mass) =
        positions
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

pub fn accumulate_forces(interactions: &[Interaction]) -> CoordMap {
    let mut btree: CoordMap = BTreeMap::new();
    for interaction in interactions {
        for particle in &interaction.forces {
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

#[allow(clippy::too_many_arguments)]
fn build_interaction(
    com_force: &[f64; 3],
    scale: f64,
    selection: &[i32],
    masses: &[f64],
    max_force: f64,
    energy: Option<f64>,
    id: Option<String>,
) -> Interaction {
    assert_eq!(selection.len(), masses.len());
    let n_particles = selection.len();
    let force_per_particle = [
        scale * com_force[0] / n_particles as f64,
        scale * com_force[1] / n_particles as f64,
        scale * com_force[2] / n_particles as f64,
    ];
    Interaction {
        forces: selection
            .iter()
            .zip(masses)
            .map(|(index, mass)| InteractionForce {
                selection: *index as usize,
                force: [
                    (force_per_particle[0] * mass).clamp(-max_force, max_force),
                    (force_per_particle[1] * mass).clamp(-max_force, max_force),
                    (force_per_particle[2] * mass).clamp(-max_force, max_force),
                ],
            })
            .collect(),
        energy: energy.map(|energy| masses.iter().map(|mass| mass * scale * energy).sum::<f64>()),
        id,
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
pub fn zeroed_out(indices: &HashSet<i32>) -> CoordMap {
    let mut btree: CoordMap = BTreeMap::new();
    indices.iter().for_each(|i| {
        btree.insert(*i, [0.0, 0.0, 0.0]);
    });
    btree
}

fn compute_gaussian_force(diff: Coordinate, sigma: f64) -> (Coordinate, f64) {
    let sigma_sqr = sigma * sigma;
    let distance_sqr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
    let gauss = (-distance_sqr / (2.0 * sigma_sqr)).exp();
    let force = [
        -(diff[0] / sigma_sqr) * gauss,
        -(diff[1] / sigma_sqr) * gauss,
        -(diff[2] / sigma_sqr) * gauss,
    ];
    let energy = -gauss;
    (force, energy)
}

fn compute_harmonic_force(diff: Coordinate, k: f64) -> (Coordinate, f64) {
    let factor = -2.0 * k;
    let distance_sqr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
    let force = [factor * diff[0], factor * diff[1], factor * diff[2]];
    let energy = k * distance_sqr;
    (force, energy)
}

unsafe fn get_selection_positions_from_state_positions(
    selection: &[i32],
    pos_state: *const OpenMM_Vec3Array,
) -> Vec<[f64; 3]> {
    selection
        .iter()
        .map(|p| {
            let position = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(pos_state, *p), 1.0);
            [position.x, position.y, position.z]
        })
        .collect()
}

fn get_masses_for_selection(selection: &[i32], masses: &[f64]) -> Vec<f64> {
    selection
        .iter()
        .map(|index| {
            let index: usize = *index as usize;
            masses[index]
        })
        .collect()
}

pub fn compute_forces_inner(
    positions: *const OpenMM_Vec3Array,
    masses: &[f64],
    n_particles: usize,
    imd_interactions: &[IMDInteraction],
) -> Vec<Interaction> {
    let interactions = imd_interactions
        .iter()
        .map(|imd| compute_forces_single(imd, positions, masses, n_particles))
        .collect();
    interactions
}

fn compute_forces_single(
    imd: &IMDInteraction,
    positions: *const OpenMM_Vec3Array,
    masses: &[f64],
    n_particles: usize,
) -> Interaction {
    let max_force = imd.max_force.unwrap_or(f64::INFINITY);
    let selection = filter_selection(&imd.particles, n_particles);
    let masses_selection = get_masses_for_selection(&selection, masses);
    let particle_positions =
        unsafe { get_selection_positions_from_state_positions(&selection, positions) };
    let com = compute_com(&particle_positions, &masses_selection);
    let interaction_position = imd.position;
    let diff = [
        com[0] - interaction_position[0],
        com[1] - interaction_position[1],
        com[2] - interaction_position[2],
    ];
    let sigma = 1.0; // For now we use this as a constant. It is in the python version.
    let (com_force, energy) = match imd.kind {
        InteractionKind::GAUSSIAN => compute_gaussian_force(diff, sigma),
        InteractionKind::HARMONIC => compute_harmonic_force(diff, sigma),
    };
    build_interaction(
        &com_force,
        imd.scale,
        &selection,
        &masses_selection,
        max_force,
        Some(energy),
        imd.id.clone(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use openmm_sys::{
        OpenMM_Vec3, OpenMM_Vec3Array, OpenMM_Vec3Array_create, OpenMM_Vec3Array_destroy,
        OpenMM_Vec3Array_set,
    };
    use rstest::rstest;
    use std::collections::HashSet;

    const EXP_1: f64 = 0.6065306597126334; // exp(-0.5)
    const EXP_3: f64 = 0.22313016014842982; // exp(-3.0/2.0)
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
    #[case([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-EXP_1, 0.0, 0.0], -EXP_1)]
    #[case([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [ EXP_1, 0.0, 0.0], -EXP_1)]
    #[case([1.0, 3.0, 0.0], [1.0, 2.0, 0.0], [0.0, -EXP_1, 0.0], -EXP_1)]
    #[case([1.0, 3.0, 3.0], [1.0, 3.0, 2.0], [0.0, 0.0, -EXP_1], -EXP_1)]
    #[case(UNIT, [0.0, 0.0, 0.0], [-EXP_1 * UNIT[0]; 3], -EXP_1)]
    fn test_gaussian_force(
        #[case] position: Coordinate,
        #[case] interaction_position: Coordinate,
        #[case] expected_force: Coordinate,
        #[case] expected_energy: f64,
    ) {
        let sigma = 1.0;
        let diff = [
            position[0] - interaction_position[0],
            position[1] - interaction_position[1],
            position[2] - interaction_position[2],
        ];
        let (force, energy) = compute_gaussian_force(diff, sigma);
        assert_f64_near!(force[0], expected_force[0]);
        assert_f64_near!(force[1], expected_force[1]);
        assert_f64_near!(force[2], expected_force[2]);
        assert_f64_near!(energy, expected_energy);
    }

    #[rstest]
    #[case([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-2.0, 0.0, 0.0], 1.0)]
    #[case([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [ 2.0, 0.0, 0.0], 1.0)]
    #[case([1.0, 3.0, 0.0], [1.0, 2.0, 0.0], [0.0, -2.0, 0.0], 1.0)]
    #[case([1.0, 3.0, 3.0], [1.0, 3.0, 2.0], [0.0, 0.0, -2.0], 1.0)]
    #[case(UNIT, [0.0, 0.0, 0.0], [-2.0 * UNIT[0]; 3], 1.0)]
    #[case([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [-2.0, -2.0, -2.0], 3.0)]
    #[case([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [0.0, 0.0, 0.0], 0.0)]
    #[case([-1.0, -1.0, -1.0], [0.0, 0.0, 0.0], [2.0, 2.0, 2.0], 3.0)]
    fn test_harmonic_force(
        #[case] position: Coordinate,
        #[case] interaction_position: Coordinate,
        #[case] expected_force: Coordinate,
        #[case] expected_energy: f64,
    ) {
        let k = 1.0;
        let diff = [
            position[0] - interaction_position[0],
            position[1] - interaction_position[1],
            position[2] - interaction_position[2],
        ];
        let (force, energy) = compute_harmonic_force(diff, k);
        assert_f64_near!(force[0], expected_force[0]);
        assert_f64_near!(force[1], expected_force[1]);
        assert_f64_near!(force[2], expected_force[2]);
        assert_f64_near!(energy, expected_energy);
    }

    #[test]
    fn test_accumulate_forces() {
        let interactions: Vec<Interaction> = vec![
            Interaction {
                forces: vec![
                    InteractionForce {
                        selection: 43,
                        force: [1.0, 2.0, 3.0],
                    },
                    InteractionForce {
                        selection: 12,
                        force: [2.1, 3.2, 4.3],
                    },
                ],
                energy: None,
                id: None,
            },
            Interaction {
                forces: vec![InteractionForce {
                    selection: 2,
                    force: [3.2, 4.3, 5.4],
                }],
                energy: None,
                id: None,
            },
            Interaction {
                forces: vec![InteractionForce {
                    selection: 43,
                    force: [4.3, 5.4, 6.5],
                }],
                energy: None,
                id: None,
            },
        ];
        let accumulated = accumulate_forces(&interactions);
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

    fn single_interaction() -> IMDInteraction {
        IMDInteraction {
            position: [0.0, 0.0, 0.0],
            particles: vec![1],
            ..Default::default()
        }
    }

    struct SampleParticles {
        positions: *mut OpenMM_Vec3Array,
        masses: Vec<f64>,
        n_particles: usize,
    }

    impl SampleParticles {
        unsafe fn new() -> Self {
            let n_particles: usize = 50;
            let positions = OpenMM_Vec3Array_create(n_particles as i32);
            (0..n_particles).for_each(|index| {
                let coord = index as f64;
                OpenMM_Vec3Array_set(
                    positions,
                    index as i32,
                    OpenMM_Vec3 {
                        x: coord,
                        y: coord,
                        z: coord,
                    },
                )
            });
            let masses = (0..n_particles).map(|i| i as f64 + 1.0).collect();
            Self {
                positions,
                masses,
                n_particles,
            }
        }
    }

    impl Drop for SampleParticles {
        fn drop(&mut self) {
            unsafe {
                OpenMM_Vec3Array_destroy(self.positions);
            }
        }
    }

    fn assert_interaction_force_near(expected: &InteractionForce, actual: &InteractionForce) {
        assert_eq!(expected.selection, actual.selection);
        expected
            .force
            .iter()
            .zip(actual.force)
            .for_each(|(expected, actual)| assert_f64_near!(*expected, actual));
    }

    #[rstest]
    #[case(-1.0)]
    #[case(0.0)]
    #[case(100.0)]
    /// Tests that the interaction force calculation gives the expected result on a single atom, at
    /// a particular position, with varying scale.
    fn test_interaction_force_single(#[case] scale: f64) {
        let mut single_interaction = single_interaction();
        single_interaction.scale = scale;
        let masses;
        let interaction = unsafe {
            let particles = SampleParticles::new();
            masses = particles.masses.clone();
            compute_forces_single(
                &single_interaction,
                particles.positions,
                &particles.masses,
                particles.n_particles,
            )
        };
        let expected_energy = -EXP_3 * scale * masses[single_interaction.particles[0]];
        let expected_interaction = Interaction {
            forces: single_interaction
                .particles
                .iter()
                .map(|index| InteractionForce {
                    selection: *index,
                    force: [-EXP_3 * scale * masses[*index]; 3],
                })
                .collect(),
            energy: Some(expected_energy),
            id: None,
        };

        assert_f64_near!(
            interaction.energy.unwrap(),
            expected_interaction.energy.unwrap()
        );
        assert!(expected_interaction.id.is_none());
        assert_eq!(interaction.forces.len(), expected_interaction.forces.len());
        interaction
            .forces
            .iter()
            .zip(expected_interaction.forces)
            .for_each(|(actual, expected)| assert_interaction_force_near(&expected, actual));
    }

    #[test]
    fn test_build_interaction() {
        let com_force = [1.0, 2.0, 3.0];
        let scale = 2.0;
        let selection = [3, 4, 5];
        let masses = vec![3.0, 4.0, 5.0];
        let max_force = 20000.0;
        let energy = 42.24;
        let id = None;

        let expected = Interaction {
            forces: selection
                .iter()
                .enumerate()
                .map(|(index, particle)| InteractionForce {
                    selection: *particle as usize,
                    force: [
                        scale * masses[index] * com_force[0] / selection.len() as f64,
                        scale * masses[index] * com_force[1] / selection.len() as f64,
                        scale * masses[index] * com_force[2] / selection.len() as f64,
                    ],
                })
                .collect(),
            energy: Some(masses.iter().map(|mass| scale * mass * energy).sum::<f64>()),
            id: id.clone(),
        };

        let actual = build_interaction(
            &com_force,
            scale,
            &selection,
            &masses,
            max_force,
            Some(energy),
            id,
        );
        expected
            .forces
            .iter()
            .zip(actual.forces)
            .for_each(|(e, a)| assert_interaction_force_near(e, &a));
    }
}
