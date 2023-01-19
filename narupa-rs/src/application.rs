extern crate clap;

use futures::TryFutureExt;
use log::{error, info};
use narupa_proto::frame::FrameData;
use crate::essd::serve_essd;
use crate::frame_broadcaster::FrameBroadcaster;
use crate::multiuser::RadialOrient;
use crate::observer_thread::run_observer_thread;
use crate::playback::PlaybackCommand;
use crate::playback::PlaybackOrder;
use crate::services::commands::{Command, CommandServer, CommandService};
use crate::services::state::{StateServer, StateService};
use crate::services::trajectory::{Trajectory, TrajectoryServiceServer};
use crate::simulation::XMLParsingError;
use crate::simulation_thread::run_simulation_thread;
use crate::simulation_thread::XMLBuffer;
use crate::state_broadcaster::StateBroadcaster;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::net::IpAddr;
use std::net::SocketAddr;
use std::sync::{Arc, Mutex};
use thiserror::Error;
use tokio::sync::mpsc::{self, Receiver, Sender};
use tonic::transport::Server;

use clap::Parser;

#[derive(Error, Debug)]
#[error("The statistic file cannot be open.")]
struct CannotOpenStatisticFile;
unsafe impl Send for CannotOpenStatisticFile {}
impl From<CannotOpenStatisticFile> for Box<dyn std::error::Error + Send> {
    fn from(_value: CannotOpenStatisticFile) -> Self {
        Box::new(CannotOpenStatisticFile)
    }
}

/// A Narupa IMD server.
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct Cli {
    /// The path to the Narupa XML file describing the simulation to run.
    #[clap(value_parser)]
    pub input_xml_path: Option<String>,
    /// IP address to bind.
    #[clap(short, long, value_parser, default_value = "0.0.0.0")]
    pub address: IpAddr,
    /// Port the server will listen.
    #[clap(short, long, value_parser, default_value_t = 38801)]
    pub port: u16,
    /// Throtle the simulation at this rate.
    #[clap(short, long, value_parser, default_value_t = 30.0)]
    pub simulation_fps: f64,
    /// Sends a frame every STEPS dynamics steps.
    #[clap(short = 'f', long, value_parser, default_value_t = 5)]
    pub frame_interval: u32,
    /// Show the simulation progression and some performance data.
    #[clap(long, value_parser, default_value_t = false)]
    pub progression: bool,
    /// Update the interactions every STEPS dynamics steps.
    #[clap(short = 'i', long, value_parser, default_value_t = 10)]
    pub force_interval: u32,
    /// Display more information about what the software does.
    #[clap(short, long, value_parser, default_value_t = false)]
    pub verbose: bool,
    /// Be very verbose about what the software does.
    #[clap(short, long, value_parser, default_value_t = false)]
    pub trace: bool,
    #[clap(long, value_parser)]
    pub statistics: Option<String>,
    #[clap(long, value_parser, default_value_t = 4.0)]
    pub statistics_fps: f64,
    /// Server name to advertise for autoconnect.
    #[clap(short, long, value_parser, default_value = "Narupa-RS iMD Server")]
    pub name: String,
}

impl Default for Cli {
    fn default() -> Self {
        Cli {
            input_xml_path: None,
            address: IpAddr::from([0, 0, 0, 0]),
            port: 38801,
            simulation_fps: 30.0,
            frame_interval: 5,
            force_interval: 10,
            progression: false,
            verbose: false,
            trace: false,
            statistics: None,
            statistics_fps: 4.0,
            name: "Narupa-RS iMD Server".to_owned(),
        }
    }
}

#[derive(Error, Debug)]
pub enum AppError {
    #[error("Cannot open the input file.")]
    CannotOpenInputFile(#[from] std::io::Error),
    #[error("Cannot parse input file.")]
    CannotParseInputFile(#[from] XMLParsingError),
    #[error("Cannot open statistics file.")]
    CannotOpenStatisticFile,
    #[error("Server cannot establish connection.")]
    TransportError(#[from] tonic::transport::Error),
    #[error("Internal server error.")]
    JoinError(#[from] tokio::task::JoinError),
}

impl From<CannotOpenStatisticFile> for AppError {
    fn from(_value: CannotOpenStatisticFile) -> Self {
        AppError::CannotOpenStatisticFile
    }
}

pub async fn main_to_wrap(cli: Cli, cancel_rx: tokio::sync::oneshot::Receiver<()>) -> Result<(), AppError> {
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
            .serve_with_shutdown(socket_address, cancel_rx.unwrap_or_else(|_| ())),
    )
    .await??;

    Ok(())
}