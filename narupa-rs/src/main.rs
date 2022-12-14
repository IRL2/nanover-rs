extern crate clap;

use narupa_rs::frame::FrameData;
use narupa_rs::frame_broadcaster::FrameBroadcaster;
use narupa_rs::services::commands::{CommandServer, CommandService};
use narupa_rs::services::state::{StateServer, StateService};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::state_broadcaster::StateBroadcaster;
use narupa_rs::simulation_thread::run_simulation_thread;
use narupa_rs::observer_thread::run_observer_thread;
use narupa_rs::playback::PlaybackOrder;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::time::Duration;
use tokio::net::UdpSocket;
use tokio::sync::mpsc::{self, Sender, Receiver};
use tokio::time;
use tonic::transport::Server;
use network_interface::NetworkInterface;
use network_interface::NetworkInterfaceConfig;


use clap::Parser;

/// A Narupa IMD server.
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// The path to the Narupa XML file describing the simulation to run.
    #[clap(value_parser, default_value = "17-ala.xml")]
    input_xml_path: String,
    /// IP address to bind.
    #[clap(short, long, value_parser, default_value = "0.0.0.0")]
    address: String,
    /// Port the server will listen.
    #[clap(short, long, value_parser, default_value_t = 38801)]
    port: usize,
    /// Throtle the simulation at this rate.
    #[clap(short, long, value_parser, default_value_t = 30.0)]
    simulation_fps: f64,
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
    #[clap(long, value_parser, default_value_t = 4.0)]
    statistics_fps: f64,
    #[clap(short, long, value_parser, default_value = "Narupa-RS iMD Server")]
    name: String,
}

async fn serve_essd(name: String, port: usize) {
    let mut interval = time::interval(Duration::from_secs_f32(0.5));
    let id = uuid::Uuid::new_v4();

    let socket = UdpSocket::bind("0.0.0.0:0").await.unwrap();
    socket.set_broadcast(true).unwrap();
    loop {
        interval.tick().await;
        let network_interfaces = NetworkInterface::show().unwrap();
        for interface in network_interfaces.iter() {
            let Some(address) = interface.addr else {continue};
            let Some(broadcast_address) = address.broadcast() else {continue};

            let server_address = address.ip();
            let message = format!("{{\"name\": \"{name}\", \"address\": \"{server_address}\", \"port\": {port}, \"id\": \"{id}\", \"essd_version\": \"1.0.0\", \"services\": {{\"imd\": {port}, \"trajectory\": {port}}}}}");
            let message = message.as_bytes();
            socket.send_to(message, format!("{broadcast_address}:54545")).await.unwrap();
        }
    }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read the user arguments.
    let cli = Cli::parse();
    let xml_path = cli.input_xml_path;
    let address = format!("{}:{}", cli.address, cli.port);
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
    let statistics_interval = ((1.0 / cli.statistics_fps) * 1000.0) as u64;

    // We have 3 separate threads: one runs the simulation, one
    // runs the GRPC server, and one observes what is happening
    // to provide some statistics. Here, we setup how the threads talk
    // to each other.
    let (frame_tx, frame_rx) = std::sync::mpsc::channel();
    let (state_tx, state_rx) = std::sync::mpsc::channel();
    let (simulation_tx, simulation_rx) = std::sync::mpsc::channel();
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(FrameBroadcaster::new(empty_frame, Some(frame_tx))));
    let shared_state = Arc::new(Mutex::new(StateBroadcaster::new(Some(state_tx))));
    let (playback_tx, playback_rx): (Sender<PlaybackOrder>, Receiver<PlaybackOrder>) = mpsc::channel(100);

    // Observe what is happening for statistics.
    if let Some(output) = statistics_file {
        run_observer_thread(output, statistics_interval, frame_rx, state_rx, simulation_rx);
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
        simulation_tx,
    );

    // Advertise the server with ESSD
    println!("Advertise the server with ESSD");
    tokio::task::spawn(serve_essd(cli.name, cli.port));

    // Run the GRPC server on the main thread.
    println!("Let's go!");
    let socket_address = address.to_socket_addrs().unwrap().next().unwrap();
    println!("Listening to {socket_address}");
    let server = Trajectory::new(Arc::clone(&frame_source));
    let command_service = CommandService::new(playback_tx);
    let state_service = StateService::new(Arc::clone(&shared_state));
    Server::builder()
        .add_service(TrajectoryServiceServer::new(server))
        .add_service(CommandServer::new(command_service))
        .add_service(StateServer::new(state_service))
        .serve(socket_address)
        .await
        .unwrap();
    Ok(())
}
