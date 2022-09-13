use std::convert::TryInto;
use std::fs::File;
use std::io::BufReader;
use std::time::Duration;
use std::{thread, time};
use std::sync::{Arc, Mutex};
use std::cmp::Ordering;
use crate::simulation::{IMDInteraction, XMLSimulation, ToFrameData, Simulation, IMD};
use crate::frame_broadcaster::FrameBroadcaster;
use crate::state_broadcaster::StateBroadcaster;
use crate::state_interaction::read_state_interaction;
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

pub fn run_simulation_thread(
        xml_path: String,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        state_clone: Arc<Mutex<StateBroadcaster>>,
        simulation_interval: u64,
        frame_interval: i32,
        force_interval: i32,
        verbose: bool,
) {
    tokio::task::spawn_blocking(move || {
        // TODO: check if there isn't a throttled iterator, otherwise write one.
        let file = File::open(xml_path).unwrap();
        let file_buffer = BufReader::new(file);
        let mut simulation = XMLSimulation::new(file_buffer);
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
            let (delta_frames, do_frames, do_forces) = next_stop(
                current_simulation_frame,
                frame_interval as u64,
                force_interval as u64,
            );
            let now = time::Instant::now();
            simulation.step(delta_frames);
            current_simulation_frame += delta_frames as u64;
            if do_frames {
                let frame = simulation.to_framedata();
                let mut source = sim_clone.lock().unwrap();
                source.send(frame).unwrap();
            }
            if do_forces {
                let state_interactions: Vec<IMDInteraction> = {
                    let state = state_clone.lock().unwrap();
                    let interaction_iter = state.iter();
                    interaction_iter
                        .filter(|kv| kv.0.starts_with("interaction."))
                        .map(|kv| {
                            let value = kv.1;
                            read_state_interaction(value)
                        })
                        .filter(|result| result.is_ok())
                        .filter_map(|interaction| interaction.ok())
                        .collect()
                };
                let imd_interactions = simulation.compute_forces(&state_interactions);
                simulation.update_imd_forces(imd_interactions).unwrap();
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
    });
}