use narupa_rs::frame::FrameData;
use narupa_rs::simulation::{Simulation, ToFrameData, TestSimulation};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::services::commands::{CommandService, CommandServer};
use narupa_rs::services::state::{StateService, StateServer};
use std::{time, thread};
use std::time::Duration;
use tonic::transport::Server;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};


#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // We have 2 separate threads: one runs the simulation, and the other one
    //runs the GRPC server. Here, we setup how the two threads talk
    // to each other.
    // TODO: actually implement the state service
    // TODO: actually implement the command service
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(empty_frame));

    // Run the simulation thread.
    // TODO: build the simulation from an input file and
    // provide the file as a CLI argument
    let sim_clone = Arc::clone(&frame_source);
    tokio::task::spawn_blocking(move || {
        // TODO: check if there isn't a throttled iterator, otherwise write one.
        // TODO: make the throttling interval a CLI argument
        let mut simulation = TestSimulation::new();
        let interval = Duration::from_millis(200);
        for i in 0.. {
            let now = time::Instant::now();
            println!("{i}");
            simulation.step(10);
            {
                let frame = simulation.to_framedata();
                let mut source = sim_clone.lock().unwrap();
                *source = frame;
            }
            let elapsed = now.elapsed();
            let time_left = interval - elapsed;
            thread::sleep(time_left);
        }
    });

    // Run the GRPC server on the main thread.
    // TODO: make the address and port CLI arguments
    println!("Let's go!");
    let server = Trajectory::new(Arc::clone(&frame_source));
    let command_service = CommandService {};
    let state_service = StateService {};
    Server::builder()
        .add_service(TrajectoryServiceServer::new(server))
        .add_service(CommandServer::new(command_service))
        .add_service(StateServer::new(state_service))
        .serve("[::]:50051".to_socket_addrs().unwrap().next().unwrap())
        .await
        .unwrap();
    Ok(())
}