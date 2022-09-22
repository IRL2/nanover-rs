use std::convert::TryInto;
use std::fs::File;
use std::io::BufReader;
use std::time::Duration;
use std::{thread, time};
use std::sync::{Arc, Mutex, mpsc::Receiver};
use std::sync::mpsc::TryRecvError;
use std::cmp::Ordering;

use crate::playback::{PlaybackOrder, PlaybackState};
use crate::simulation::{XMLSimulation, ToFrameData, Simulation, IMD};
use crate::frame_broadcaster::FrameBroadcaster;
use crate::state_broadcaster::StateBroadcaster;
use crate::state_interaction::read_forces;
use crate::broadcaster::Broadcaster;

fn next_stop(current_frame: u64, frame_interval: u64, force_interval: u64) -> (i32, bool, bool) {
    let next_frame_stop = frame_interval - current_frame % frame_interval;
    let next_force_stop = force_interval - current_frame % force_interval;
    match next_frame_stop.cmp(&next_force_stop) {
        Ordering::Less => (next_frame_stop.try_into().unwrap(), true, false),
        Ordering::Greater => (next_force_stop.try_into().unwrap(), false, true),
        Ordering::Equal => (next_frame_stop.try_into().unwrap(), true, true),
    }
}

fn apply_forces(state_clone: &Arc<Mutex<StateBroadcaster>>, simulation: &mut XMLSimulation) {
    let state_interactions = read_forces(state_clone);
    let imd_interactions = simulation.compute_forces(&state_interactions);
    simulation.update_imd_forces(imd_interactions).unwrap();
}

pub fn run_simulation_thread(
        xml_path: String,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        state_clone: Arc<Mutex<StateBroadcaster>>,
        simulation_interval: u64,
        frame_interval: i32,
        force_interval: i32,
        verbose: bool,
        playback_rx: Receiver<PlaybackOrder>,
) {
    tokio::task::spawn_blocking(move || {
        // TODO: check if there isn't a throttled iterator, otherwise write one.
        let file = File::open(xml_path).unwrap();
        let file_buffer = BufReader::new(file);
        let mut simulation = XMLSimulation::new(file_buffer);
        let mut playback_state = PlaybackState::new(true);
        let interval = Duration::from_millis(simulation_interval);
        {
            sim_clone
                .lock()
                .unwrap()
                .send(simulation.to_topology_framedata())
                .unwrap();
        }
        println!("Platform: {}", simulation.get_platform_name());
        println!("Simulation interval: {}", simulation_interval);
        println!("Start simulating");
        let mut current_simulation_frame: u64 = 0;
        for _i in 0.. {
            let now = time::Instant::now();
            loop {
                let order_result = playback_rx.try_recv();
                match order_result {
                    Ok(PlaybackOrder::Reset) => simulation.reset(),
                    Ok(order) => playback_state.update(order),
                    // The queue of order is empty so we are done handling them.
                    Err(TryRecvError::Empty) => break,
                    // The server thread is done so we sould end the simulation.
                    Err(TryRecvError::Disconnected) => return,
                }
            }

            if playback_state.is_playing() {
                let (delta_frames, do_frames, do_forces) = next_stop(
                    current_simulation_frame,
                    frame_interval as u64,
                    force_interval as u64,
                );
                simulation.step(delta_frames);
                current_simulation_frame += delta_frames as u64;
                if do_frames {
                    let frame = simulation.to_framedata();
                    let mut source = sim_clone.lock().unwrap();
                    source.send(frame).unwrap();
                }
                if do_forces {
                    apply_forces(&state_clone, &mut simulation);
                }

                let elapsed = now.elapsed();
                let time_left = match interval.checked_sub(elapsed) {
                    Some(d) => d,
                    None => Duration::from_millis(0),
                };
                if verbose {
                    println!("Simulation frame {current_simulation_frame}. Time to sleep {time_left:?}");
                };
                thread::sleep(time_left);
            }
        }
    });
}