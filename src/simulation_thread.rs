use std::fs::File;
use std::io::BufReader;
use std::time::Duration;
use std::{thread, time};
use std::sync::{Arc, Mutex};
use crate::simulation::{IMDInteraction, XMLSimulation, ToFrameData, Simulation, IMD};
use crate::frame_broadcaster::FrameBroadcaster;
use crate::state_broadcaster::StateBroadcaster;
use crate::state_interaction::read_state_interaction;
use crate::broadcaster::Broadcaster;

pub fn run_simulation_thread(
        xml_path: String,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        state_clone: Arc<Mutex<StateBroadcaster>>,
        simulation_interval: u64,
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
        for _i in 0.. {
            let now = time::Instant::now();
            //println!("{i}");
            simulation.step(10);
            {
                let frame = simulation.to_framedata();
                let mut source = sim_clone.lock().unwrap();
                source.send(frame).unwrap();
            }
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
                    //.inspect(|interaction| {
                    //    println!("{interaction:?}");
                    //})
                    .collect()
            };
            let imd_interactions = simulation.compute_forces(&state_interactions);

            simulation.update_imd_forces(imd_interactions).unwrap();
            let elapsed = now.elapsed();
            let time_left = match interval.checked_sub(elapsed) {
                Some(d) => d,
                None => Duration::from_millis(0),
            };
            //println!("Time to sleep {time_left:?}");
            thread::sleep(time_left);
        }
    });
}