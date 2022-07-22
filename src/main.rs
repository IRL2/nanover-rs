use narupa_rs::frame::FrameData;
use narupa_rs::broadcaster::Broadcaster;
use narupa_rs::frame_broadcaster::FrameBroadcaster;
use narupa_rs::state_broadcaster::StateBroadcaster;
use narupa_rs::simulation::{
    Simulation,
    ToFrameData,
    XMLSimulation,
    IMD,
    IMDInteraction,
    InteractionKind,
};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::services::commands::{CommandService, CommandServer};
use narupa_rs::services::state::{StateService, StateServer};
use std::{time, thread};
use std::time::Duration;
use tonic::transport::Server;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::BufReader;
use prost_types::value::Kind;


fn read_state_interaction(state_interaction: &prost_types::Value) -> Result<IMDInteraction, ()> {
    // Extract the interaction content that is under several layers of enums.
    // Fail already if we cannot: it means the input does not have the expected format.
    let content = match &state_interaction.kind {
        Some(Kind::StructValue(inside)) => {inside},
        _ => return Err(()),
    };

    let kind = content.fields.get("interaction_type");
    if kind.is_none() {return Err(())};
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
        Some(inside) => {
            match inside.kind {
                Some(Kind::NumberValue(value)) => Some(value),
                _ => return Err(()),
            }
        }
    };

    let scale = content.fields.get("scale");
    let scale = match scale {
        None => 1.0,
        Some(inside) => {
            match inside.kind {
                Some(Kind::NumberValue(value)) => value,
                _ => return Err(()),
            }
        }
    };

    let particles = content.fields.get("particles");
    let particles: Vec<usize> = match particles {
        Some(inside) => {
            match &inside.kind {
                Some(Kind::ListValue(values)) => {
                    values.values
                        .iter()
                        .filter_map(|v| v.kind.as_ref())
                        // We just ignore invalid values
                        .filter_map(|v| match v {
                            Kind::NumberValue(inner) => Some((*inner) as usize),
                            _ => None,
                        })
                        .collect()
                },
                _ => return Err(()),
            }
        }
        _ => return Err(()),
    };

    let position = content.fields.get("position");
    let position: Vec<f64> = match position {
        Some(inside) => match &inside.kind {
            Some(Kind::ListValue(values)) => {
                values.values
                    .iter()
                    .filter_map(|v| v.kind.as_ref())
                    .filter_map(|v| match v {
                        Kind::NumberValue(inner) => Some(*inner),
                        _ => None,
                    })
                    .collect()
            },
            _ => return Err(()),
        },
        _ => return Err(()),
    };
    if position.len() != 3 {return Err(())}
    let position: [f64; 3] = position.try_into().unwrap();

    Ok(IMDInteraction::new(position, particles, kind, max_force, scale))
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // We have 2 separate threads: one runs the simulation, and the other one
    // runs the GRPC server. Here, we setup how the two threads talk
    // to each other.
    // TODO: actually implement the state service
    // TODO: actually implement the command service
    let empty_frame = FrameData::empty();
    //let frame_source = Arc::new(Mutex::new(empty_frame));
    let frame_source = Arc::new(Mutex::new(FrameBroadcaster::new(empty_frame)));

    let shared_state = Arc::new(Mutex::new(StateBroadcaster::new()));

    // Run the simulation thread.
    // TODO: build the simulation from an input file and
    // provide the file as a CLI argument
    let sim_clone = Arc::clone(&frame_source);
    let state_clone = Arc::clone(&shared_state);
    tokio::task::spawn_blocking(move || {
        // TODO: check if there isn't a throttled iterator, otherwise write one.
        // TODO: make the throttling interval a CLI argument
        //let mut simulation = TestSimulation::new();
        let file = File::open("17-ala.xml").unwrap();
        let file_buffer = BufReader::new(file);
        let mut simulation = XMLSimulation::new(file_buffer);
        let interval = Duration::from_millis(33);
        {
            sim_clone.lock().unwrap()
                .send(simulation.to_topology_framedata()).unwrap();
        }
        println!("Platform: {}", simulation.get_platform_name());
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

    // Run the GRPC server on the main thread.
    // TODO: make the address and port CLI arguments
    println!("Let's go!");
    let server = Trajectory::new(Arc::clone(&frame_source));
    let command_service = CommandService {};
    let state_service = StateService::new(Arc::clone(&shared_state));
    Server::builder()
        .add_service(TrajectoryServiceServer::new(server))
        .add_service(CommandServer::new(command_service))
        .add_service(StateServer::new(state_service))
        .serve("[::]:38801".to_socket_addrs().unwrap().next().unwrap())
        .await
        .unwrap();
    Ok(())
}