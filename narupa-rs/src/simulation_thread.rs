use log::{info, warn, error};
use std::cmp::Ordering;
use std::convert::TryInto;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::{thread, time};
use tokio::sync::mpsc::{error::TryRecvError, Receiver};

use crate::frame_broadcaster::FrameBroadcaster;
use crate::manifest::{Manifest, LoadDefaultError, LoadSimulationError};
use crate::playback::{PlaybackOrder, PlaybackState};
use crate::simulation::{OpenMMSimulation, Simulation, ToFrameData, XMLParsingError, IMD};
use crate::state_broadcaster::StateBroadcaster;
use crate::state_interaction::read_forces;

fn next_stop(current_frame: u64, frame_interval: u32, force_interval: u32) -> (i32, bool, bool) {
    let frame_interval: u64 = frame_interval as u64;
    let force_interval: u64 = force_interval as u64;
    let next_frame_stop = frame_interval - current_frame % frame_interval;
    let next_force_stop = force_interval - current_frame % force_interval;
    match next_frame_stop.cmp(&next_force_stop) {
        Ordering::Less => (next_frame_stop.try_into().unwrap(), true, false),
        Ordering::Greater => (next_force_stop.try_into().unwrap(), false, true),
        Ordering::Equal => (next_frame_stop.try_into().unwrap(), true, true),
    }
}

fn next_frame_stop(current_frame: u64, frame_interval: u32) -> i32 {
    let frame_interval: u64 = frame_interval as u64;
    (frame_interval - current_frame % frame_interval)
        .try_into()
        .unwrap()
}

fn apply_forces(
    state_clone: &Arc<Mutex<StateBroadcaster>>,
    simulation: &mut OpenMMSimulation,
    simulation_tx: std::sync::mpsc::Sender<usize>,
) {
    let state_interactions = read_forces(state_clone);
    let imd_interactions = simulation.compute_forces(&state_interactions);
    simulation_tx.send(imd_interactions.len()).unwrap();
    simulation.update_imd_forces(imd_interactions).unwrap();
}

pub fn run_simulation_thread(
    mut simulations_manifest: Manifest,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
    state_clone: Arc<Mutex<StateBroadcaster>>,
    simulation_interval: u64,
    frame_interval: u32,
    force_interval: u32,
    verbose: bool,
    mut playback_rx: Receiver<PlaybackOrder>,
    simulation_tx: std::sync::mpsc::Sender<usize>,
    auto_reset: bool,
) -> Result<(), XMLParsingError> {
    let mut maybe_simulation: Option<OpenMMSimulation> = match simulations_manifest.load_default() {
        Ok(simulation) => Some(simulation),
        Err(LoadDefaultError::NoDefault) => None,
        Err(LoadDefaultError::LoadSimulationError(load_simulation_error)) => {
            match load_simulation_error {
                LoadSimulationError::CannotOpen(_) => None,
                LoadSimulationError::NoIndex(_) => None,
                LoadSimulationError::XMLParsingError(error) => return Err(error),
            }
        }
    };

    tokio::task::spawn_blocking(move || {
        let mut playback_state = PlaybackState::new(true);
        let interval = Duration::from_millis(simulation_interval);
        if let Some(ref simulation) = maybe_simulation {
            if sim_clone
                .lock()
                .unwrap()
                .send_reset_frame(simulation.to_topology_framedata())
                .is_err()
            {
                return;
            }
            info!("Platform: {}", simulation.get_platform_name());
            info!("Start simulating");
        } else {
            info!("No simulation loaded yes.");
        }
        info!("Simulation interval: {simulation_interval}");

        let mut current_simulation_frame: u64 = 0;
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
                                next_frame_stop(current_simulation_frame, frame_interval);
                            simulation.step(delta_frames);
                            current_simulation_frame += delta_frames as u64;
                            let frame = simulation.to_framedata();
                            sim_clone.lock().unwrap().send_frame(frame).unwrap();
                            playback_state.update(PlaybackOrder::Step);
                        } else {
                            warn!("No simulation loaded, ignoring STEP command.");
                        }
                    }
                    Ok(PlaybackOrder::Load(simulation_index)) => {
                        maybe_simulation = match simulations_manifest.load_index(simulation_index) {
                            Ok(new_simulation) => Some(new_simulation),
                            Err(error) => {
                                error!("Could not load simulation with index {simulation_index}: {error}");
                                warn!("No new simulation loaded, keep using the previously loaded one if any.");
                                maybe_simulation
                            }
                        }
                    }
                    Ok(PlaybackOrder::Next) => {
                        maybe_simulation = match simulations_manifest.load_next() {
                            Ok(new_simulation) => Some(new_simulation),
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
                    let (delta_frames, do_frames, do_forces) =
                        next_stop(current_simulation_frame, frame_interval, force_interval);
                    simulation.step(delta_frames);
                    current_simulation_frame += delta_frames as u64;
                    if do_frames {
                        let frame = simulation.to_framedata();
                        let mut source = sim_clone.lock().unwrap();
                        if source.send_frame(frame).is_err() {
                            return;
                        };
                    }
                    if do_forces {
                        apply_forces(&state_clone, simulation, simulation_tx.clone());
                    }

                    let elapsed = now.elapsed();
                    let time_left = match interval.checked_sub(elapsed) {
                        Some(d) => d,
                        None => Duration::from_millis(0),
                    };
                    let energy = simulation.get_total_energy();
                    if verbose {
                        info!(
                            "Simulation frame {current_simulation_frame}. Time to sleep {time_left:?}. Total energy {energy:.2} kJ/mol."
                        );
                    };
                    if auto_reset && !energy.is_finite() {
                        simulation.reset();
                    }
                    thread::sleep(time_left);
                }
            }
        }
    });

    Ok(())
}
