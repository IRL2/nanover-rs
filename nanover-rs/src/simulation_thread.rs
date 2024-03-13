use log::{error, info, trace, warn};
use std::sync::{Arc, Mutex};
use std::time;
use tokio::sync::mpsc::{error::TryRecvError, Receiver};

use crate::frame_broadcaster::FrameBroadcaster;
use crate::manifest::{LoadDefaultError, LoadSimulationError, Manifest};
use crate::playback::{PlaybackOrder, PlaybackState};
use crate::simulation::XMLParsingError;
use crate::state_broadcaster::StateBroadcaster;
use crate::tracked_simulation::specific::{
    next_simulation_counter, send_reset_frame, SpecificSimulationTracked,
};
use crate::tracked_simulation::Configuration;

fn playback_loop(
    playback_rx: &mut Receiver<PlaybackOrder>,
    playback_state: &mut PlaybackState,
    simulations_manifest: &mut Manifest,
    mut maybe_simulation: Option<SpecificSimulationTracked>,
    sim_clone: &Arc<Mutex<FrameBroadcaster>>,
    configuration: &Configuration,
) -> (Option<SpecificSimulationTracked>, bool) {
    let mut keep_going = true;
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
                    simulation.step();
                    if simulation
                        .send_regular_frame(
                            Arc::clone(sim_clone),
                            configuration.with_velocities,
                            configuration.with_forces,
                        )
                        .is_err()
                    {
                        keep_going = false;
                        break;
                    }
                    playback_state.update(PlaybackOrder::Step);
                } else {
                    warn!("No simulation loaded, ignoring STEP command.");
                }
            }
            Ok(PlaybackOrder::Load(simulation_index)) => {
                maybe_simulation = match simulations_manifest.load_index(simulation_index) {
                    Ok(new_simulation) => {
                        let simulation = SpecificSimulationTracked::new(
                            new_simulation,
                            next_simulation_counter(maybe_simulation.as_ref()),
                            configuration,
                        );
                        if send_reset_frame(&simulation, Arc::clone(sim_clone)).is_err() {
                            keep_going = false;
                            break;
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
                        let simulation = SpecificSimulationTracked::new(
                            new_simulation,
                            next_simulation_counter(maybe_simulation.as_ref()),
                            configuration,
                        );
                        if send_reset_frame(&simulation, Arc::clone(sim_clone)).is_err() {
                            keep_going = false;
                            break;
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
            Err(TryRecvError::Disconnected) => {
                keep_going = false;
                break;
            }
        }
    }

    (maybe_simulation, keep_going)
}

// TODO: Redo error handling
fn load_initial_simulation(
    simulations_manifest: &mut Manifest,
    configuration: &Configuration,
) -> Result<Option<SpecificSimulationTracked>, XMLParsingError> {
    match simulations_manifest.load_default() {
        Ok(simulation) => Ok(Some(SpecificSimulationTracked::new(
            simulation,
            0,
            configuration,
        ))),
        Err(LoadDefaultError::NoDefault) => Ok(None),
        Err(LoadDefaultError::LoadSimulationError(load_simulation_error)) => {
            match load_simulation_error {
                LoadSimulationError::CannotOpen(err) => {
                    error!("Cannot read simulation file: {err}");
                    Ok(None)
                }
                LoadSimulationError::NoIndex(index) => {
                    error!("No simulation with index {index}.");
                    Ok(None)
                }
                LoadSimulationError::XMLParsingError(error) => Err(error),
                LoadSimulationError::ReadError(error) => {
                    error!("An error occured while reading the recording: {error}.");
                    Ok(None)
                }
                LoadSimulationError::UnsuportedFormatVersion(version) => {
                    error!("The file use a version of the recording format that is unsuported ({version}).");
                    Ok(None)
                }
                LoadSimulationError::NotARecording => {
                    error!("The file is not a recording.");
                    Ok(None)
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run_simulation_thread(
    mut cancel_rx: tokio::sync::oneshot::Receiver<()>,
    mut simulations_manifest: Manifest,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
    state_clone: Arc<Mutex<StateBroadcaster>>,
    mut playback_rx: Receiver<PlaybackOrder>,
    simulation_tx: std::sync::mpsc::Sender<usize>,
    run_on_start: bool,
    configuration: Configuration,
) -> Result<(), XMLParsingError> {
    let mut maybe_simulation = load_initial_simulation(&mut simulations_manifest, &configuration)?;

    tokio::task::spawn_blocking(move || {
        let mut playback_state = PlaybackState::new(run_on_start);
        if let Some(ref simulation) = maybe_simulation {
            if send_reset_frame(simulation, Arc::clone(&sim_clone)).is_err() {
                return;
            }
            if let SpecificSimulationTracked::OpenMM(openmm_simulation) = &simulation {
                info!("Platform: {}", openmm_simulation.get_platform_name());
            }
            info!("Start simulating");
        } else {
            info!("No simulation loaded yes.");
        }
        info!(
            "Simulation interval: {} ms",
            configuration.simulation_interval.as_millis()
        );

        loop {
            let now = time::Instant::now();
            match cancel_rx.try_recv() {
                Ok(_) | Err(tokio::sync::oneshot::error::TryRecvError::Closed) => {
                    trace!("Simulation loop ended.");
                    break;
                }
                Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {
                    trace!("Tick simulation loop");
                }
            };
            let keep_going;
            (maybe_simulation, keep_going) = playback_loop(
                &mut playback_rx,
                &mut playback_state,
                &mut simulations_manifest,
                maybe_simulation,
                &sim_clone,
                &configuration,
            );
            if !keep_going {
                return;
            }

            if let Some(ref mut simulation) = maybe_simulation {
                match simulation {
                    SpecificSimulationTracked::Recording(simulation) => {
                        if playback_state.is_playing() {
                            simulation.simulation_loop_iteration(&sim_clone, &now, &configuration)
                        }
                    }
                    SpecificSimulationTracked::OpenMM(simulation) => {
                        if playback_state.is_playing() {
                            simulation.simulation_loop_iteration(
                                &sim_clone,
                                &state_clone,
                                &simulation_tx,
                                &now,
                                &configuration,
                            )
                        }
                    }
                }
            }
        }
    });

    Ok(())
}
