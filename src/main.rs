extern crate clap;

use narupa_rs::frame::FrameData;
use narupa_rs::frame_broadcaster::FrameBroadcaster;
use narupa_rs::services::commands::{CommandServer, CommandService};
use narupa_rs::services::state::{StateServer, StateService};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::state_broadcaster::StateBroadcaster;
use narupa_rs::simulation_thread::run_simulation_thread;
use std::convert::TryInto;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};
use tonic::transport::Server;

use clap::Parser;

/// A Narupa IMD server.
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// The path to the Narupa XML file describing the simulation to run.
    #[clap(value_parser, default_value = "17-ala.xml")]
    input_xml_path: String,
    /// Port the server will listen.
    #[clap(short, long, value_parser, default_value_t = 38801)]
    port: usize,
    /// Throtle the simulation at this rate.
    #[clap(short, long, value_parser, default_value_t = 30)]
    simulation_fps: usize,
    /// Sends a frame every STEPS dynamics steps.
    #[clap(short, long, value_parser, default_value_t = 5)]
    frame_interval: u32,
    /// Display simulation advancement.
    #[clap(short, long, value_parser, default_value_t = false)]
    verbose: bool,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read the user arguments.
    let cli = Cli::parse();
    let xml_path = cli.input_xml_path;
    let port = cli.port;
    let simulation_interval = ((1.0 / (cli.simulation_fps) as f64) * 1000.0) as u64;
    let frame_interval: i32 = cli
        .frame_interval
        .try_into()
        .expect("Invalid frame interval value.");
    let verbose = cli.verbose;

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
    run_simulation_thread(
        xml_path,
        sim_clone,
        state_clone,
        simulation_interval,
        frame_interval,
        verbose,
    );

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
