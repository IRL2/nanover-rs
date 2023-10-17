use log::{error, info, warn};
use narupa_proto::trajectory::FrameData;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::{thread, time};
use tokio::sync::mpsc::{error::TryRecvError, Receiver};

use crate::broadcaster::BroadcastSendError;
use crate::frame_broadcaster::FrameBroadcaster;
use crate::manifest::{LoadDefaultError, LoadSimulationError, Manifest};
use crate::playback::{PlaybackOrder, PlaybackState};
use crate::simulation::{
    CoordMap, OpenMMSimulation, Simulation, ToFrameData, XMLParsingError, IMD,
};
use crate::state_broadcaster::StateBroadcaster;
use crate::state_interaction::read_forces;

/// A simulation with all the book keeping required for the thread.
struct TrackedSimulation {
    simulation: OpenMMSimulation,
    simulation_frame: usize,
    reset_counter: usize,
    simulation_counter: usize,
    user_forces: CoordMap,
    user_energies: f64,
}

impl TrackedSimulation {
    fn new(simulation: OpenMMSimulation, simulation_counter: usize) -> Self {
        let user_forces = BTreeMap::new();
        let user_energies = 0.0;
        let reset_counter = 0;
        let simulation_frame = 0;
        Self {
            simulation,
            simulation_frame,
            reset_counter,
            simulation_counter,
            user_forces,
            user_energies,
        }
    }

    fn reset(&mut self) {
        self.simulation.reset();
        self.reset_counter += 1;
    }

    fn step(&mut self, n_frames: usize) {
        self.simulation.step(n_frames as i32);
        self.simulation_frame += n_frames;
    }

    fn reset_counter(&self) -> usize {
        self.reset_counter
    }

    fn simulation_counter(&self) -> usize {
        self.simulation_counter
    }

    fn simulation(&self) -> &OpenMMSimulation {
        &self.simulation
    }

    fn simulation_mut(&mut self) -> &mut OpenMMSimulation {
        &mut self.simulation
    }

    fn user_forces(&self) -> &CoordMap {
        &self.user_forces
    }

    fn update_user_forces(&mut self, force_map: CoordMap) {
        self.user_forces = force_map;
    }

    fn user_energies(&self) -> f64 {
        self.user_energies
    }

    fn update_user_energies(&mut self, user_energies: f64) {
        self.user_energies = user_energies;
    }

    fn simulation_frame(&self) -> usize {
        self.simulation_frame
    }
}

fn next_stop(
    current_frame: usize,
    frame_interval: usize,
    force_interval: usize,
) -> (usize, bool, bool) {
    let next_frame_stop = frame_interval - current_frame % frame_interval;
    let next_force_stop = force_interval - current_frame % force_interval;
    match next_frame_stop.cmp(&next_force_stop) {
        Ordering::Less => (next_frame_stop, true, false),
        Ordering::Greater => (next_force_stop, false, true),
        Ordering::Equal => (next_frame_stop, true, true),
    }
}

fn next_frame_stop(current_frame: usize, frame_interval: usize) -> usize {
    frame_interval - current_frame % frame_interval
}

fn starting_frame(simulation: &TrackedSimulation) -> FrameData {
    let mut frame = simulation.simulation().to_topology_framedata();
    frame
        .insert_number_value(
            "system.simulation.counter",
            simulation.simulation_counter() as f64,
        )
        .unwrap();
    frame
        .insert_number_value("system.reset.counter", simulation.reset_counter() as f64)
        .unwrap();
    frame
}

fn apply_forces(
    state_clone: &Arc<Mutex<StateBroadcaster>>,
    simulation: &mut OpenMMSimulation,
    simulation_tx: std::sync::mpsc::Sender<usize>,
) -> (CoordMap, f64) {
    let state_interactions = read_forces(state_clone);
    let imd_interactions = simulation.compute_forces(&state_interactions);
    simulation_tx.send(imd_interactions.len()).unwrap();
    let forces = simulation.update_imd_forces(&imd_interactions).unwrap();
    let interactions = imd_interactions
        .into_iter()
        .filter_map(|interaction| interaction.energy)
        .sum();
    (forces, interactions)
}

fn add_force_map_to_frame(force_map: &CoordMap, frame: &mut FrameData) {
    let indices = force_map.keys().map(|index| *index as u32).collect();
    let sparse = force_map
        .values()
        .flatten()
        .map(|value| *value as f32)
        .collect();
    frame
        .insert_index_array("forces.user.index", indices)
        .unwrap();
    frame
        .insert_float_array("forces.user.sparse", sparse)
        .unwrap();
}

fn send_regular_frame(
    simulation: &TrackedSimulation,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
    with_velocities: bool,
    with_forces: bool,
) -> Result<(), BroadcastSendError> {
    let mut frame = simulation
        .simulation
        .to_framedata(with_velocities, with_forces);
    let energy_total = simulation.user_energies();
    frame
        .insert_number_value("energy.user.total", energy_total)
        .unwrap();
    frame
        .insert_number_value("system.reset.counter", simulation.reset_counter() as f64)
        .unwrap();
    add_force_map_to_frame(simulation.user_forces(), &mut frame);
    let mut source = sim_clone.lock().unwrap();
    source.send_frame(frame)
}

fn send_reset_frame(
    simulation: &TrackedSimulation,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
) -> Result<(), BroadcastSendError> {
    let frame = starting_frame(simulation);
    sim_clone.lock().unwrap().send_reset_frame(frame)
}

fn next_simulation_counter(maybe_simulation: Option<&TrackedSimulation>) -> usize {
    maybe_simulation
        .map(|simulation| simulation.simulation_counter() + 1)
        .unwrap_or_default()
}

#[allow(clippy::too_many_arguments)]
pub fn run_simulation_thread(
    mut simulations_manifest: Manifest,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
    state_clone: Arc<Mutex<StateBroadcaster>>,
    simulation_interval: u64,
    frame_interval: usize,
    force_interval: usize,
    verbose: bool,
    mut playback_rx: Receiver<PlaybackOrder>,
    simulation_tx: std::sync::mpsc::Sender<usize>,
    auto_reset: bool,
    run_on_start: bool,
    with_velocities: bool,
    with_forces: bool,
) -> Result<(), XMLParsingError> {
    let mut maybe_simulation: Option<TrackedSimulation> = match simulations_manifest.load_default()
    {
        Ok(simulation) => Some(TrackedSimulation::new(simulation, 0)),
        Err(LoadDefaultError::NoDefault) => None,
        Err(LoadDefaultError::LoadSimulationError(load_simulation_error)) => {
            match load_simulation_error {
                LoadSimulationError::CannotOpen(err) => {
                    error!("Cannot read simulation file: {err}");
                    None
                }
                LoadSimulationError::NoIndex(index) => {
                    error!("No simulation with index {index}.");
                    None
                }
                LoadSimulationError::XMLParsingError(error) => return Err(error),
            }
        }
    };

    tokio::task::spawn_blocking(move || {
        let mut playback_state = PlaybackState::new(run_on_start);
        let interval = Duration::from_millis(simulation_interval);
        if let Some(ref simulation) = maybe_simulation {
            if send_reset_frame(simulation, Arc::clone(&sim_clone)).is_err() {
                return;
            }
            info!("Platform: {}", simulation.simulation.get_platform_name());
            info!("Start simulating");
        } else {
            info!("No simulation loaded yes.");
        }
        info!("Simulation interval: {simulation_interval}");

        loop {
            let now = time::Instant::now();
            loop {
                let order_result = playback_rx.try_recv();
                match order_result {
                    Ok(PlaybackOrder::Reset) => {
                        if let Some(ref mut simulation) = maybe_simulation {
                            simulation.reset()
                        } else {
                            warn!("No simulation loaded, ignoring RESET command.");
                        };
                    }
                    Ok(PlaybackOrder::Step) => {
                        if let Some(ref mut simulation) = maybe_simulation {
                            let delta_frames =
                                next_frame_stop(simulation.simulation_frame(), frame_interval);
                            simulation.step(delta_frames);
                            if send_regular_frame(
                                simulation,
                                Arc::clone(&sim_clone),
                                with_velocities,
                                with_forces,
                            )
                            .is_err()
                            {
                                return;
                            }
                            playback_state.update(PlaybackOrder::Step);
                        } else {
                            warn!("No simulation loaded, ignoring STEP command.");
                        }
                    }
                    Ok(PlaybackOrder::Load(simulation_index)) => {
                        maybe_simulation = match simulations_manifest.load_index(simulation_index) {
                            Ok(new_simulation) => {
                                let simulation = TrackedSimulation::new(
                                    new_simulation,
                                    next_simulation_counter(maybe_simulation.as_ref()),
                                );
                                if send_reset_frame(&simulation, Arc::clone(&sim_clone)).is_err() {
                                    return;
                                }

                                Some(simulation)
                            }
                            Err(error) => {
                                error!("Could not load simulation with index {simulation_index}: {error}");
                                warn!("No new simulation loaded, keep using the previously loaded one if any.");
                                maybe_simulation
                            }
                        }
                    }
                    Ok(PlaybackOrder::Next) => {
                        maybe_simulation = match simulations_manifest.load_next() {
                            Ok(new_simulation) => {
                                let simulation = TrackedSimulation::new(
                                    new_simulation,
                                    next_simulation_counter(maybe_simulation.as_ref()),
                                );
                                if send_reset_frame(&simulation, Arc::clone(&sim_clone)).is_err() {
                                    return;
                                }

                                Some(simulation)
                            }
                            Err(error) => {
                                error!("Could not load the next simulation: {error}");
                                warn!("No new simulation loaded, keep using the previously loaded one if any.");
                                maybe_simulation
                            }
                        }
                    }
                    Ok(order) => playback_state.update(order),
                    // The queue of order is empty so we are done handling them.
                    Err(TryRecvError::Empty) => break,
                    // The server thread is done so we sould end the simulation.
                    Err(TryRecvError::Disconnected) => return,
                }
            }

            if let Some(ref mut simulation) = maybe_simulation {
                if playback_state.is_playing() {
                    let (delta_frames, do_frames, do_forces) = next_stop(
                        simulation.simulation_frame(),
                        frame_interval,
                        force_interval,
                    );
                    simulation.step(delta_frames);

                    if do_forces {
                        let (force_map, user_energies) = apply_forces(
                            &state_clone,
                            simulation.simulation_mut(),
                            simulation_tx.clone(),
                        );
                        simulation.update_user_forces(force_map);
                        simulation.update_user_energies(user_energies);
                    }

                    let system_energy = simulation.simulation.get_total_energy();

                    if do_frames
                        && send_regular_frame(
                            simulation,
                            Arc::clone(&sim_clone),
                            with_velocities,
                            with_forces,
                        )
                        .is_err()
                    {
                        return;
                    };

                    let elapsed = now.elapsed();
                    let time_left = match interval.checked_sub(elapsed) {
                        Some(d) => d,
                        None => Duration::from_millis(0),
                    };
                    if verbose {
                        info!(
                            "Simulation frame {}. Time to sleep {time_left:?}. Total energy {system_energy:.2} kJ/mol.",
                            simulation.simulation_frame(),
                        );
                    };
                    if auto_reset && !system_energy.is_finite() {
                        simulation.reset();
                    }
                    thread::sleep(time_left);
                }
            }
        }
    });

    Ok(())
}
