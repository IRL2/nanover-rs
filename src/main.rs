extern crate clap;

use narupa_rs::broadcaster::Broadcaster;
use narupa_rs::frame::FrameData;
use narupa_rs::frame_broadcaster::FrameBroadcaster;
use narupa_rs::services::commands::{CommandServer, CommandService};
use narupa_rs::services::state::{StateServer, StateService};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::simulation::{
    IMDInteraction, InteractionKind, Simulation, ToFrameData, XMLSimulation, IMD,
};
use narupa_rs::state_broadcaster::StateBroadcaster;
use prost_types::value::Kind;
use std::fs::File;
use std::io::BufReader;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::{thread, time};
use tonic::transport::Server;

use clap::Parser;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(value_parser, default_value = "17-ala.xml")]
    input_xml_path: String,
    #[clap(short, long, value_parser, default_value_t = 38801)]
    port: usize,
    #[clap(short, long, value_parser, default_value_t = 30)]
    simulation_fps: usize,
}

fn read_state_interaction(state_interaction: &prost_types::Value) -> Result<IMDInteraction, ()> {
    // Extract the interaction content that is under several layers of enums.
    // Fail already if we cannot: it means the input does not have the expected format.
    let content = match &state_interaction.kind {
        Some(Kind::StructValue(inside)) => inside,
        _ => return Err(()),
    };

    let kind = content.fields.get("interaction_type");
    if kind.is_none() {
        return Err(());
    };
    let kind = kind.unwrap();
    let kind = match &kind.kind {
        Some(Kind::StringValue(inside)) => inside,
        _ => return Err(()),
    };
    let kind = if kind == "gaussian" {
        InteractionKind::GAUSSIAN
    } else if kind == "spring" {
        InteractionKind::HARMONIC
    } else {
        return Err(());
    };

    let max_force = content.fields.get("max_force");
    let max_force = match max_force {
        None => None,
        Some(inside) => match inside.kind {
            Some(Kind::NumberValue(value)) => Some(value),
            _ => return Err(()),
        },
    };

    let scale = content.fields.get("scale");
    let scale = match scale {
        None => 1.0,
        Some(inside) => match inside.kind {
            Some(Kind::NumberValue(value)) => value,
            _ => return Err(()),
        },
    };

    let particles = content.fields.get("particles");
    let particles: Vec<usize> = match particles {
        Some(inside) => {
            match &inside.kind {
                Some(Kind::ListValue(values)) => {
                    values
                        .values
                        .iter()
                        .filter_map(|v| v.kind.as_ref())
                        // We just ignore invalid values
                        .filter_map(|v| match v {
                            Kind::NumberValue(inner) => Some((*inner) as usize),
                            _ => None,
                        })
                        .collect()
                }
                _ => return Err(()),
            }
        }
        _ => return Err(()),
    };

    let position = content.fields.get("position");
    let position: Vec<f64> = match position {
        Some(inside) => match &inside.kind {
            Some(Kind::ListValue(values)) => values
                .values
                .iter()
                .filter_map(|v| v.kind.as_ref())
                .filter_map(|v| match v {
                    Kind::NumberValue(inner) => Some(*inner),
                    _ => None,
                })
                .collect(),
            _ => return Err(()),
        },
        _ => return Err(()),
    };
    if position.len() != 3 {
        return Err(());
    }
    let position: [f64; 3] = position.try_into().unwrap();

    Ok(IMDInteraction::new(
        position, particles, kind, max_force, scale,
    ))
}

fn run_simulation_thread(
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

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let xml_path = cli.input_xml_path;
    let port = cli.port;
    let simulation_interval = ((1.0 / (cli.simulation_fps) as f64) * 1000.0) as u64;

    // We have 2 separate threads: one runs the simulation, and the other one
    // runs the GRPC server. Here, we setup how the two threads talk
    // to each other.
    // TODO: actually implement the command service
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(FrameBroadcaster::new(empty_frame)));
    let shared_state = Arc::new(Mutex::new(StateBroadcaster::new()));

    // Run the simulation thread.
    let sim_clone = Arc::clone(&frame_source);
    let state_clone = Arc::clone(&shared_state);
    run_simulation_thread(xml_path, sim_clone, state_clone, simulation_interval);

    // Run the GRPC server on the main thread.
    println!("Let's go!");
    let address = format!("[::]:{port}");
    let server = Trajectory::new(Arc::clone(&frame_source));
    let command_service = CommandService {};
    let state_service = StateService::new(Arc::clone(&shared_state));
    Server::builder()
        .add_service(TrajectoryServiceServer::new(server))
        .add_service(CommandServer::new(command_service))
        .add_service(StateServer::new(state_service))
        .serve(address.to_socket_addrs().unwrap().next().unwrap())
        .await
        .unwrap();
    Ok(())
}
