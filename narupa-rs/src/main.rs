extern crate clap;

use futures::TryFutureExt;
use narupa_rs::essd::serve_essd;
use narupa_proto::frame::FrameData;
use narupa_rs::frame_broadcaster::FrameBroadcaster;
use narupa_rs::multiuser::RadialOrient;
use narupa_rs::observer_thread::run_observer_thread;
use narupa_rs::playback::PlaybackCommand;
use narupa_rs::playback::PlaybackOrder;
use narupa_rs::services::commands::{Command, CommandServer, CommandService};
use narupa_rs::services::state::{StateServer, StateService};
use narupa_rs::services::trajectory::{Trajectory, TrajectoryServiceServer};
use narupa_rs::simulation_thread::XMLBuffer;
use narupa_rs::simulation_thread::run_simulation_thread;
use narupa_rs::state_broadcaster::StateBroadcaster;
use std::collections::HashMap;
use std::fs::File;
use std::net::IpAddr;
use std::net::SocketAddr;
use std::error::Error as StdError;
use std::sync::{Arc, Mutex};
use std::process::ExitCode;
use std::io::BufReader;
use tokio::sync::mpsc::{self, Receiver, Sender};
use tonic::transport::Server;
use thiserror::Error;
use log::{error, info};
use env_logger::Builder;
use log::LevelFilter;

use clap::Parser;

#[derive(Error, Debug)]
#[error("The statistic file cannot be open.")]
struct CannotOpenStatisticFile;

/// A Narupa IMD server.
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// The path to the Narupa XML file describing the simulation to run.
    #[clap(value_parser)]
    input_xml_path: Option<String>,
    /// IP address to bind.
    #[clap(short, long, value_parser, default_value = "0.0.0.0")]
    address: IpAddr,
    /// Port the server will listen.
    #[clap(short, long, value_parser, default_value_t = 38801)]
    port: u16,
    /// Throtle the simulation at this rate.
    #[clap(short, long, value_parser, default_value_t = 30.0)]
    simulation_fps: f64,
    /// Sends a frame every STEPS dynamics steps.
    #[clap(short = 'f', long, value_parser, default_value_t = 5)]
    frame_interval: u32,
    /// Show the simulation progression and some performance data.
    #[clap(long, value_parser, default_value_t = false)]
    progression: bool,
    /// Update the interactions every STEPS dynamics steps.
    #[clap(short = 'i', long, value_parser, default_value_t = 10)]
    force_interval: u32,
    /// Display more information about what the software does.
    #[clap(short, long, value_parser, default_value_t = false)]
    verbose: bool,
    /// Be very verbose about what the software does.
    #[clap(short, long, value_parser, default_value_t = false)]
    trace: bool,
    #[clap(long, value_parser)]
    statistics: Option<String>,
    #[clap(long, value_parser, default_value_t = 4.0)]
    statistics_fps: f64,
    /// Server name to advertise for autoconnect.
    #[clap(short, long, value_parser, default_value = "Narupa-RS iMD Server")]
    name: String,
}

async fn main_to_wrap(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {

    let (cancel_tx, cancel_rx) = tokio::sync::oneshot::channel();
    tokio::spawn(async move {
        tokio::signal::ctrl_c().await.unwrap();
        // Your handler here
        info!("Closing the server. Goodbye!");
        cancel_tx.send(()).unwrap();
    });

    // Read the user arguments.
    let xml_path = cli.input_xml_path;
    let simulation_interval = ((1.0 / cli.simulation_fps) * 1000.0) as u64;
    let frame_interval = cli.frame_interval;
    let force_interval = cli.force_interval;
    let verbose = cli.progression;
    let statistics_file = cli
        .statistics
        .map(File::create)
        .transpose()
        .map_err(|_| CannotOpenStatisticFile)?;
    let statistics_interval = ((1.0 / cli.statistics_fps) * 1000.0) as u64;
    let socket_address = SocketAddr::new(cli.address, cli.port);

    // We have 3 separate threads: one runs the simulation, one
    // runs the GRPC server, and one observes what is happening
    // to provide some statistics. Here, we setup how the threads talk
    // to each other.
    let (frame_tx, frame_rx) = std::sync::mpsc::channel();
    let (state_tx, state_rx) = std::sync::mpsc::channel();
    let (simulation_tx, simulation_rx) = std::sync::mpsc::channel();
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(FrameBroadcaster::new(
        empty_frame,
        Some(frame_tx),
    )));
    let shared_state = Arc::new(Mutex::new(StateBroadcaster::new(Some(state_tx))));
    let (playback_tx, playback_rx): (Sender<PlaybackOrder>, Receiver<PlaybackOrder>) =
        mpsc::channel(100);

    let mut commands: HashMap<String, Box<dyn Command>> = HashMap::new();
    commands.insert(
        "playback/play".into(),
        Box::new(PlaybackCommand::new(
            playback_tx.clone(),
            PlaybackOrder::Play,
        )),
    );
    commands.insert(
        "playback/pause".into(),
        Box::new(PlaybackCommand::new(
            playback_tx.clone(),
            PlaybackOrder::Pause,
        )),
    );
    commands.insert(
        "playback/reset".into(),
        Box::new(PlaybackCommand::new(
            playback_tx.clone(),
            PlaybackOrder::Reset,
        )),
    );
    commands.insert(
        "playback/step".into(),
        Box::new(PlaybackCommand::new(
            playback_tx.clone(),
            PlaybackOrder::Step,
        )),
    );
    commands.insert(
        "multiuser/radially-orient-origins".into(),
        Box::new(RadialOrient::new(Arc::clone(&shared_state))),
    );

    // Observe what is happening for statistics.
    if let Some(output) = statistics_file {
        run_observer_thread(
            output,
            statistics_interval,
            frame_rx,
            state_rx,
            simulation_rx,
        );
    };

    // Run the simulation thread.
    let sim_clone = Arc::clone(&frame_source);
    let state_clone = Arc::clone(&shared_state);
    let xml_buffer = if let Some(path) = xml_path {
        let xml_file = File::open(path)?;
        XMLBuffer::FileBuffer(BufReader::new(xml_file))
    } else {
        let bytes = include_bytes!("../17-ala.xml");
        XMLBuffer::BytesBuffer(BufReader::new(bytes))
    };
    run_simulation_thread(
        xml_buffer,
        sim_clone,
        state_clone,
        simulation_interval,
        frame_interval,
        force_interval,
        verbose,
        playback_rx,
        simulation_tx,
        true,
    )?;

    // Advertise the server with ESSD
    info!("Advertise the server with ESSD");
    tokio::task::spawn(serve_essd(cli.name, cli.port));


    // Run the GRPC server on the main thread.
    info!("Listening to {socket_address}");
    let server = Trajectory::new(Arc::clone(&frame_source));
    let command_service = CommandService::new(commands);
    let state_service = StateService::new(Arc::clone(&shared_state));
    tokio::task::spawn(
        Server::builder()
            .add_service(TrajectoryServiceServer::new(server))
            .add_service(CommandServer::new(command_service))
            .add_service(StateServer::new(state_service))
            .serve_with_shutdown(socket_address, cancel_rx.unwrap_or_else(|_| ()))
    )
    .await??;

    Ok(())
}

#[tokio::main]
async fn main() -> ExitCode {
    let cli = Cli::parse();

    if std::env::var("RUST_LOG").is_ok() {
        env_logger::init();
    } else {
        let mut verbosity_level = LevelFilter::Info;
        if cli.verbose {verbosity_level = LevelFilter::Debug};
        if cli.trace {verbosity_level = LevelFilter::Trace};

        let mut builder = Builder::new();
        builder
            .filter_module("narupa_rs", verbosity_level)
            .format_target(false)
            .init();
    }

    let run_status = main_to_wrap(cli).await;
    let Err(ref error) = run_status else {
        return ExitCode::SUCCESS;
    };

    // The Display trait from tonic's errors is not very expressive for the
    // end user. We need to dig out the underlying error.
    let maybe_transport_error = error.downcast_ref::<tonic::transport::Error>();
    let maybe_source_error = match maybe_transport_error {
        None => None,
        Some(transport_error) => transport_error.source(),
    };
    let maybe_hyper_error = match maybe_source_error {
        None => None,
        Some(source_error) => source_error.downcast_ref::<hyper::Error>(),
    };
    if let Some(hyper_error) = maybe_hyper_error {
        error!("{hyper_error}");
        return ExitCode::FAILURE;
    };

    error!("{error}");
    ExitCode::FAILURE
}