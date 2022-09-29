extern crate clap;

use narupa_rs::broadcaster::BroadcasterSignal;
use narupa_rs::frame::FrameData;
use narupa_rs::frame_broadcaster::FrameBroadcaster;
use narupa_rs::services::commands::{CommandServer, CommandService};
use narupa_rs::services::state::{StateServer, StateService};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::state_broadcaster::StateBroadcaster;
use narupa_rs::simulation_thread::run_simulation_thread;
use narupa_rs::playback::PlaybackOrder;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use tokio::sync::mpsc::{self, Sender, Receiver};
use tonic::transport::Server;

use std::fs::File;
use std::io::Write;

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
    #[clap(short='f', long, value_parser, default_value_t = 5)]
    frame_interval: u32,
    /// Update the interactions every STEPS dynamics steps.
    #[clap(short='i', long, value_parser, default_value_t = 10)]
    force_interval: u32,
    /// Display simulation advancement.
    #[clap(short, long, value_parser, default_value_t = false)]
    verbose: bool,
    #[clap(long, value_parser)]
    statistics: Option<String>,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read the user arguments.
    let cli = Cli::parse();
    let xml_path = cli.input_xml_path;
    let port = cli.port;
    let simulation_interval = ((1.0 / (cli.simulation_fps) as f64) * 1000.0) as u64;
    let frame_interval = cli .frame_interval;
    let force_interval = cli .force_interval;
    let verbose = cli.verbose;
    let statistics_file = cli
        .statistics
        .map_or(
            None,
            |path| Some(
                File::create(path)
                    .expect("Cannot open statistics file.")
            )
        );

    // We have 2 separate threads: one runs the simulation, and the other one
    // runs the GRPC server. Here, we setup how the two threads talk
    // to each other.
    let (frame_tx, frame_rx) = std::sync::mpsc::channel();
    let (state_tx, _state_rx) = std::sync::mpsc::channel();
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(FrameBroadcaster::new(empty_frame, Some(frame_tx))));
    let shared_state = Arc::new(Mutex::new(StateBroadcaster::new(Some(state_tx))));

    let (playback_tx, playback_rx): (Sender<PlaybackOrder>, Receiver<PlaybackOrder>) = mpsc::channel(100);

    if let Some(mut output) = statistics_file {
        tokio::task::spawn_blocking(move || {
            let start = Instant::now();
            loop {
                let frame_signal = frame_rx.try_recv();
                match frame_signal {
                    Ok(BroadcasterSignal::Send(instant)) => {
                        write!(output, "{:?}\n", instant.duration_since(start)).unwrap();
                    },
                    Err(std::sync::mpsc::TryRecvError::Empty) => (),
                    Err(std::sync::mpsc::TryRecvError::Disconnected) => break,
                    Ok(_) => (), 
                };
            };
        });
    };

    // Run the simulation thread.
    let sim_clone = Arc::clone(&frame_source);
    let state_clone = Arc::clone(&shared_state);
    run_simulation_thread(
        xml_path,
        sim_clone,
        state_clone,
        simulation_interval,
        frame_interval,
        force_interval,
        verbose,
        playback_rx,
    );

    // Run the GRPC server on the main thread.
    println!("Let's go!");
    let address = format!("[::]:{port}");
    let server = Trajectory::new(Arc::clone(&frame_source));
    let command_service = CommandService::new(playback_tx);
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
