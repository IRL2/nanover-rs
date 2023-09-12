extern crate clap;

use crate::broadcaster::BroadcastReceiver;
use crate::broadcaster::Broadcaster;
use crate::essd::serve_essd;
use crate::frame_broadcaster::FrameBroadcaster;
use crate::manifest::Manifest;
use crate::multiuser::RadialOrient;
use crate::observer_thread::run_observer_thread;
use crate::playback::ListSimulations;
use crate::playback::LoadCommand;
use crate::playback::PlaybackCommand;
use crate::playback::PlaybackOrder;
use crate::services::commands::{Command, CommandServer, CommandService};
use crate::services::state::{StateServer, StateService};
use crate::services::trajectory::{Trajectory, TrajectoryServiceServer};
use crate::simulation::XMLParsingError;
use crate::simulation_thread::run_simulation_thread;
use crate::state_broadcaster::StateBroadcaster;
use futures::TryFutureExt;
use indexmap::IndexMap;
use log::{debug, error, info, trace};
use narupa_proto::frame::FrameData;
use narupa_proto::trajectory::GetFrameResponse;
use prost::Message;
use std::fs::File;
use std::net::IpAddr;
use std::net::SocketAddr;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::time::Instant;
use thiserror::Error;
use tokio::io::AsyncWriteExt;
use tokio::sync::mpsc::{self, Receiver, Sender};
use tokio::sync::oneshot;
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

pub struct CancellationReceivers {
    server: tokio::sync::oneshot::Receiver<()>,
    trajectory: tokio::sync::oneshot::Receiver<()>,
    state: tokio::sync::oneshot::Receiver<()>,
    traj_service: tokio::sync::oneshot::Receiver<()>,
    state_service: tokio::sync::oneshot::Receiver<()>,
}

impl CancellationReceivers {
    #[allow(clippy::type_complexity)]
    pub fn unpack(
        self,
    ) -> (
        oneshot::Receiver<()>,
        oneshot::Receiver<()>,
        oneshot::Receiver<()>,
        oneshot::Receiver<()>,
        oneshot::Receiver<()>,
    ) {
        (
            self.server,
            self.trajectory,
            self.state,
            self.traj_service,
            self.state_service,
        )
    }
}

#[derive(Debug, Error)]
#[error("Not all cancellation requests could be sent.")]
pub struct CancellationError {}

pub struct CancellationSenders {
    server: tokio::sync::oneshot::Sender<()>,
    trajectory: tokio::sync::oneshot::Sender<()>,
    state: tokio::sync::oneshot::Sender<()>,
    traj_service: tokio::sync::oneshot::Sender<()>,
    state_service: tokio::sync::oneshot::Sender<()>,
}

impl CancellationSenders {
    pub fn send(self) -> Result<(), CancellationError> {
        self.server.send(()).map_err(|_| CancellationError {})?;
        self.trajectory.send(()).map_err(|_| CancellationError {})?;
        self.state.send(()).map_err(|_| CancellationError {})?;
        self.traj_service
            .send(())
            .map_err(|_| CancellationError {})?;
        self.state_service
            .send(())
            .map_err(|_| CancellationError {})?;
        Ok(())
    }
}

pub fn cancellation_channels() -> (CancellationSenders, CancellationReceivers) {
    let (server_tx, server_rx) = tokio::sync::oneshot::channel();
    let (trajectory_tx, trajectory_rx) = tokio::sync::oneshot::channel();
    let (state_tx, state_rx) = tokio::sync::oneshot::channel();
    let (traj_service_tx, traj_service_rx) = tokio::sync::oneshot::channel();
    let (state_service_tx, state_service_rx) = tokio::sync::oneshot::channel();
    (
        CancellationSenders {
            server: server_tx,
            trajectory: trajectory_tx,
            state: state_tx,
            traj_service: traj_service_tx,
            state_service: state_service_tx,
        },
        CancellationReceivers {
            server: server_rx,
            trajectory: trajectory_rx,
            state: state_rx,
            traj_service: traj_service_rx,
            state_service: state_service_rx,
        },
    )
}

/// A Narupa IMD server.
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct Cli {
    /// The path to the Narupa XML file describing the simulation to run.
    #[clap(value_parser)]
    pub input_xml_path: Vec<String>,
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
    /// Record the trajecory
    #[clap(long, value_parser)]
    pub trajectory: Option<String>,
    /// Record the updates from the shared state
    #[clap(long, value_parser)]
    pub state: Option<String>,
    /// Start the simulation paused
    #[clap(long, value_parser, default_value_t = false)]
    pub start_paused: bool,
}

impl Default for Cli {
    fn default() -> Self {
        Cli {
            input_xml_path: Vec::new(),
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
            trajectory: None,
            state: None,
            start_paused: false,
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
    #[error("Cannot open the trajectory file.")]
    CannotOpenTrajectoryFile(std::io::Error),
    #[error("Cannot open the state file.")]
    CannotOpenStateFile(std::io::Error),
}

impl From<CannotOpenStatisticFile> for AppError {
    fn from(_value: CannotOpenStatisticFile) -> Self {
        AppError::CannotOpenStatisticFile
    }
}

async fn record_broadcaster<T>(
    start: Instant,
    receiver: Arc<Mutex<BroadcastReceiver<T>>>,
    mut output: tokio::fs::File,
    mut cancel_rx: tokio::sync::oneshot::Receiver<()>,
    id: usize,
) -> std::io::Result<()>
where
    T: Message,
{
    // The file starts with a MAGIC_NUMBER to indicate it is the expected file
    // format. This is followed by a version number. This was introduced with
    // version 2, so files that do not start with the magic number may be
    // version 1.
    const MAGIC_NUMBER: u64 = 6661355757386708963;
    const FORMAT_VERSION: u64 = 2;
    output.write_all(&MAGIC_NUMBER.to_le_bytes()).await?;
    output.write_all(&FORMAT_VERSION.to_le_bytes()).await?;

    let duration = Duration::from_millis(33);
    let mut interval = tokio::time::interval(duration);
    loop {
        match cancel_rx.try_recv() {
            Ok(_) | Err(tokio::sync::oneshot::error::TryRecvError::Closed) => break,
            Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {}
        };
        interval.tick().await;
        let Some(frame) = receiver.lock().unwrap().recv() else {
            continue;
        };
        let timestamp = Instant::now().saturating_duration_since(start).as_micros();
        let now_bytes = timestamp.to_le_bytes();
        let mut buffer: Vec<u8> = Vec::new();
        frame.encode(&mut buffer)?;
        // println!("Frame: {buffer:?}");
        let frame_byte_size = (buffer.len() as u64).to_le_bytes();
        output.write_all(&now_bytes).await?;
        output.write_all(&frame_byte_size).await?;
        output.write_all(&buffer).await?;
        trace!("Record {id}");
    }

    println!("Terminate recorder {id}");
    Ok(())
}

pub async fn main_to_wrap(cli: Cli, cancel_rx: CancellationReceivers) -> Result<(), AppError> {
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
    let frame_index = 0;
    let empty_frame = GetFrameResponse {
        frame_index,
        frame: Some(FrameData::empty()),
    };
    let frame_source = Arc::new(Mutex::new(FrameBroadcaster::new(
        empty_frame,
        Some(frame_tx),
    )));
    let shared_state = Arc::new(Mutex::new(StateBroadcaster::new(Some(state_tx))));
    let (playback_tx, playback_rx): (Sender<PlaybackOrder>, Receiver<PlaybackOrder>) =
        mpsc::channel(100);

    let simulation_manifest = if !xml_path.is_empty() {
        if xml_path.len() == 1 {
            info!("Running {}", xml_path.get(0).unwrap());
        } else {
            info! {"Running a queue of {} simulations.", xml_path.len()};
            xml_path.iter().for_each(|path| debug!("* {path}"));
        }
        Manifest::from_simulation_xml_paths(xml_path)
    } else {
        let bytes = include_bytes!("../17-ala.xml");
        info!("Running the demo simulation.");
        Manifest::from_simulation_xml_bytes(bytes)
    };

    let mut commands: IndexMap<String, Box<dyn Command>> = IndexMap::new();
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
        "playback/next".into(),
        Box::new(PlaybackCommand::new(
            playback_tx.clone(),
            PlaybackOrder::Next,
        )),
    );
    commands.insert(
        "playback/load".into(),
        Box::new(LoadCommand::new(playback_tx.clone())),
    );
    commands.insert(
        "playback/list".into(),
        Box::new(ListSimulations::new(simulation_manifest.list_simulations())),
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

    let (
        cancel_server_rx,
        cancel_trajectory_rx,
        cancel_state_rx,
        cancel_traj_serv_rx,
        cancel_state_serv_rx,
    ) = cancel_rx.unpack();
    let syncronous_start = Instant::now();
    if let Some(path) = cli.trajectory {
        let file = tokio::fs::File::create(path)
            .await
            .map_err(AppError::CannotOpenTrajectoryFile)?;
        let receiver = frame_source.lock().unwrap().get_rx();
        tokio::spawn(record_broadcaster(
            syncronous_start,
            receiver,
            file,
            cancel_trajectory_rx,
            0,
        ));
    }
    if let Some(path) = cli.state {
        let file = tokio::fs::File::create(path)
            .await
            .map_err(AppError::CannotOpenStateFile)?;
        let receiver = shared_state.lock().unwrap().get_rx();
        tokio::spawn(record_broadcaster(
            syncronous_start,
            receiver,
            file,
            cancel_state_rx,
            1,
        ));
    }

    // Run the simulation thread.
    let sim_clone = Arc::clone(&frame_source);
    let state_clone = Arc::clone(&shared_state);

    run_simulation_thread(
        simulation_manifest,
        sim_clone,
        state_clone,
        simulation_interval,
        frame_interval,
        force_interval,
        verbose,
        playback_rx,
        simulation_tx,
        true,
        !cli.start_paused,
    )?;

    // Advertise the server with ESSD
    info!("Advertise the server with ESSD");
    tokio::task::spawn(serve_essd(cli.name, cli.port));

    // Run the GRPC server on the main thread.
    info!("Listening to {socket_address}");
    let server = Trajectory::new(Arc::clone(&frame_source), cancel_traj_serv_rx);
    let command_service = CommandService::new(commands);
    let state_service = StateService::new(Arc::clone(&shared_state), cancel_state_serv_rx);
    tokio::task::spawn(
        Server::builder()
            .add_service(TrajectoryServiceServer::new(server))
            .add_service(CommandServer::new(command_service))
            .add_service(StateServer::new(state_service))
            .serve_with_shutdown(socket_address, cancel_server_rx.unwrap_or_else(|_| ())),
    )
    .await??;

    Ok(())
}
