use std::cmp::Ordering;
use std::convert::TryInto;
use std::fs::File;
use std::io::{self, BufReader};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::{thread, time};
use tokio::sync::mpsc::{error::TryRecvError, Receiver};
use thiserror::Error;
use log::info;

use crate::broadcaster::Broadcaster;
use crate::frame_broadcaster::FrameBroadcaster;
use crate::playback::{PlaybackOrder, PlaybackState};
use crate::simulation::{Simulation, ToFrameData, OpenMMSimulation, IMD, XMLParsingError};
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

#[derive(Error, Debug)]
pub enum SimulationSetupError {
    #[error("Cannot open the input file: {0}")]
    InputFileIOError(#[from] io::Error),
    #[error("Cannot parse the input file: {0}")]
    CannotParse(#[from] XMLParsingError),
}

pub fn run_simulation_thread(
    xml_path: String,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
    state_clone: Arc<Mutex<StateBroadcaster>>,
    simulation_interval: u64,
    frame_interval: u32,
    force_interval: u32,
    verbose: bool,
    mut playback_rx: Receiver<PlaybackOrder>,
    simulation_tx: std::sync::mpsc::Sender<usize>,
    auto_reset: bool,
) -> Result<(), SimulationSetupError> {
    let file = File::open(xml_path)?;
    let file_buffer = BufReader::new(file);
    let mut simulation = OpenMMSimulation::from_xml(file_buffer)?;

    tokio::task::spawn_blocking(move || {
        let mut playback_state = PlaybackState::new(true);
        let interval = Duration::from_millis(simulation_interval);
        {
            if let Err(_) = sim_clone
                .lock()
                .unwrap()
                .send(simulation.to_topology_framedata()) {
                    return;
                }
        }
        info!("Platform: {}", simulation.get_platform_name());
        info!("Simulation interval: {simulation_interval}");
        info!("Start simulating");
        let mut current_simulation_frame: u64 = 0;
        loop {
            let now = time::Instant::now();
            loop {
                let order_result = playback_rx.try_recv();
                match order_result {
                    Ok(PlaybackOrder::Reset) => simulation.reset(),
                    Ok(PlaybackOrder::Step) => {
                        let delta_frames =
                            next_frame_stop(current_simulation_frame, frame_interval);
                        simulation.step(delta_frames);
                        current_simulation_frame += delta_frames as u64;
                        let frame = simulation.to_framedata();
                        sim_clone.lock().unwrap().send(frame).unwrap();
                        playback_state.update(PlaybackOrder::Step);
                    }
                    Ok(order) => playback_state.update(order),
                    // The queue of order is empty so we are done handling them.
                    Err(TryRecvError::Empty) => break,
                    // The server thread is done so we sould end the simulation.
                    Err(TryRecvError::Disconnected) => return,
                }
            }

            if playback_state.is_playing() {
                let (delta_frames, do_frames, do_forces) =
                    next_stop(current_simulation_frame, frame_interval, force_interval);
                simulation.step(delta_frames);
                current_simulation_frame += delta_frames as u64;
                if do_frames {
                    let frame = simulation.to_framedata();
                    let mut source = sim_clone.lock().unwrap();
                    if let Err(_) = source.send(frame) {return};
                }
                if do_forces {
                    apply_forces(&state_clone, &mut simulation, simulation_tx.clone());
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
    });

    Ok(())
}
