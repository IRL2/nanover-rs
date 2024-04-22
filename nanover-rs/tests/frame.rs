use std::collections::BTreeMap;
use std::collections::HashMap;
use std::net::{IpAddr, SocketAddr};
use std::time::{Duration, Instant};

use log::info;
use nanover_rs::application::InputPath;
use prost_types::Struct;
use serde::Deserialize;
use tokio::net::TcpListener;
use tokio::net::UdpSocket;
use tokio::time::timeout;
use tokio_stream::StreamExt;
use tonic::transport::Channel;
use tonic::Streaming;
use uuid::Uuid;

use nanover_proto::command::command_client::CommandClient;
use nanover_proto::command::CommandMessage;
use nanover_proto::trajectory::{
    trajectory_service_client::TrajectoryServiceClient, GetFrameRequest,
};
use nanover_proto::trajectory::{FrameData, GetFrameResponse};
use nanover_rs::application::{
    cancellation_channels, main_to_wrap, AppError, CancellationSenders, Cli,
};
use nanover_rs::test_ressource;
use pack_prost::{ToProstValue, UnPack};

use test_log::test;

struct Server {
    handle: tokio::task::JoinHandle<Result<(), AppError>>,
    cancel_tx: Option<CancellationSenders>,
    name: String,
    port: u16,
}

impl Server {
    async fn new(arguments: Cli) -> Self {
        let name = arguments.name.clone();
        let (cancel_tx, cancel_rx) = cancellation_channels();
        let socket_address = SocketAddr::new(arguments.address, arguments.port);
        let listener = TcpListener::bind(socket_address).await.unwrap();
        let address = listener.local_addr().unwrap();
        let handle = tokio::spawn(main_to_wrap(arguments, cancel_rx, Some(listener)));
        Server {
            handle,
            cancel_tx: Some(cancel_tx),
            name,
            port: address.port(),
        }
    }

    pub fn stop(&mut self) {
        // We should be able to stop the server several time. Calling this
        // method means we want the server to be stopped, not that we want
        // that specific call to stop the server. Therefore, we can ignore
        // the send failing or the transmitter being None.
        if let Some(tx) = self.cancel_tx.take() {
            tx.send().ok();
        };
    }

    pub fn close(mut self) -> tokio::task::JoinHandle<Result<(), AppError>> {
        self.stop();
        self.handle
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn port(&self) -> u16 {
        self.port
    }
}

struct Client {
    trajectory: TrajectoryServiceClient<Channel>,
    command: CommandClient<Channel>,
}

impl Client {
    async fn new(port: u16) -> Result<Self, tonic::transport::Error> {
        let endpoint = format!("http://127.0.0.1:{port}");
        info!("Client connecting to {endpoint}");
        let trajectory = TrajectoryServiceClient::connect(endpoint.clone()).await?;
        let command = CommandClient::connect(endpoint.clone()).await?;
        Ok(Self {
            trajectory,
            command,
        })
    }

    async fn reset(&mut self) {
        let reset_request = CommandMessage {
            name: "playback/reset".into(),
            arguments: None,
        };
        self.command.run_command(reset_request).await.unwrap();
    }

    async fn step(&mut self) {
        let reset_request = CommandMessage {
            name: "playback/step".into(),
            arguments: None,
        };
        self.command.run_command(reset_request).await.unwrap();
    }

    async fn next(&mut self) {
        let next_request = CommandMessage {
            name: "playback/next".into(),
            arguments: None,
        };
        self.command.run_command(next_request).await.unwrap();
    }

    async fn load(&mut self, index: usize) {
        let arguments = Struct {
            fields: BTreeMap::from([("index".into(), (index as f64).to_prost_value())]),
        };

        let load_request = CommandMessage {
            name: "playback/load".into(),
            arguments: Some(arguments),
        };
        self.command.run_command(load_request).await.unwrap();
    }

    async fn step_frame(&mut self, frames: &mut Streaming<GetFrameResponse>) -> GetFrameResponse {
        self.step().await;
        self.next_frame(frames).await
    }

    async fn next_frame(&mut self, frames: &mut Streaming<GetFrameResponse>) -> GetFrameResponse {
        let frame_response = timeout(Duration::from_secs(1), frames.next())
            .await
            .unwrap() // Panic if we reached the timeout
            .unwrap() // Panic if we reached the end of the stream
            .unwrap(); // Panic if there was a transport error
        let frame = frame_response.frame.as_ref().unwrap();
        info!("Frame ({}): {frame:?}", frame_response.frame_index);
        frame_response
    }
}

async fn create_server_client_pair() -> Result<(Server, Client), tonic::transport::Error> {
    let essd_name = format!("Test NanoVer server {}", Uuid::new_v4());
    let port: u16 = 0;
    let arguments = Cli {
        port,
        // Some tests need to identify this server in the ESSD records.
        name: essd_name,
        // Some tests need to switch simulation so we load several.
        input_xml_path: vec![
            InputPath::OpenMM(test_ressource!("17-ala.xml").into()),
            InputPath::OpenMM(test_ressource!("buckyballs.xml").into()),
            InputPath::OpenMM(test_ressource!("helicene.xml").into()),
        ],
        // Some tests need to access the actual first frame so we start paused.
        start_paused: true,
        ..Default::default()
    };
    let server = Server::new(arguments).await;
    info!("Server started on port {}", server.port);
    tokio::time::sleep(Duration::from_micros(100)).await;
    let client = Client::new(server.port()).await?;

    Ok((server, client))
}

fn get_number(frame: &FrameData, key: &str) -> f64 {
    frame.values.get(key).unwrap().unpack().unwrap()
}

fn has_value(frame: &FrameData, key: &str) -> bool {
    frame.values.get(key).is_some()
}

pub const ESSD_DEFAULT_PORT: u16 = 54545;
const MAXIMUM_MESSAGE_SIZE: usize = 1024;

#[derive(Deserialize, PartialEq, Debug, Clone)]
pub struct ServiceHub {
    name: String,
    address: IpAddr,
    id: String,
    essd_version: String,
    services: HashMap<String, u16>,
}

impl ServiceHub {
    pub fn new(
        name: String,
        address: IpAddr,
        id: String,
        essd_version: String,
        services: HashMap<String, u16>,
    ) -> Self {
        Self {
            name,
            address,
            id,
            essd_version,
            services,
        }
    }

    pub fn trajectory(&self) -> Option<u16> {
        self.services.get("trajectory").copied()
    }

    pub fn multiplayer(&self) -> Option<u16> {
        self.services.get("multiplayer").copied()
    }

    pub fn address(&self) -> IpAddr {
        self.address
    }

    pub fn name(&self) -> String {
        self.name.clone()
    }

    pub fn id(&self) -> String {
        self.id.clone()
    }

    pub fn essd_version(&self) -> String {
        self.essd_version.clone()
    }

    pub fn services(&self) -> HashMap<String, u16> {
        self.services.clone()
    }
}

pub async fn find_servers(
    port: u16,
    search_time: Duration,
) -> std::io::Result<HashMap<String, ServiceHub>> {
    // We need a tokio socket to work with async, that socket needs the "reuse
    // port" flag set which needs a socket2 socket. Socket2's socket are tricky
    // to build, standard library's socket are easier to create.
    let address = format!("0.0.0.0:{port}").parse::<SocketAddr>().unwrap();
    let std_socket = std::net::UdpSocket::bind(address)?;
    let flag_socket: socket2::Socket = std_socket.into();
    flag_socket.set_nonblocking(true)?;
    flag_socket.set_broadcast(true)?;

    // set_reuse_port is not available on Windows
    // #[cfg(not(target_os = "windows"))]
    // flag_socket.set_reuse_port(true)?;

    let socket = UdpSocket::from_std(flag_socket.into())?;

    let mut buffer = [0u8; MAXIMUM_MESSAGE_SIZE];
    let mut servers = HashMap::new();
    let start_time = Instant::now();
    while start_time.elapsed() < search_time {
        let n_bytes = timeout(Duration::from_millis(100), socket.recv(&mut buffer))
            .await
            .unwrap_or(Ok(0))?;
        if n_bytes > 0 {
            let hub: ServiceHub = serde_json::from_slice(&buffer[..n_bytes])?;
            servers.insert(hub.id.clone(), hub);
        }
    }
    Ok(servers)
}

#[test(tokio::test)]
async fn test_simulation_counter() {
    let (server, mut client) = create_server_client_pair().await.unwrap();

    let request = GetFrameRequest::default();
    let mut frames = client
        .trajectory
        .subscribe_latest_frames(request)
        .await
        .unwrap()
        .into_inner();

    // We should not need to call playback/step to get a frame here. The first frame is sent before
    // we start looking for blackback orders.
    // TODO: Fix this so client.next_frame is enough here.
    let mut frame_response = client.step_frame(&mut frames).await;
    let mut frame = frame_response.frame.as_ref().unwrap();
    //We sometimes get the frame too early. In that case, we get the initial frame the broadcaster
    //was built with as the frame of interest was not sent yet. If this is the case, we want to
    //ignore that frame.
    if frame.arrays.is_empty() && frame.values.is_empty() {
        frame_response = client.next_frame(&mut frames).await;
        frame = frame_response.frame.as_ref().unwrap();
    };

    // This is a frame from a new simulation, the frame_index should be 0 to tell the client to
    // clear the aggregated frame.
    assert_eq!(frame_response.frame_index, 0);
    // This is the first simulation loaded by the server. The counter must be set to 0.
    let simulation_counter: f64 = get_number(frame, "system.simulation.counter");
    let simulation_counter = simulation_counter as usize;
    assert_eq!(simulation_counter, 0);
    // We have not reset yet, so the reset counter must be 0.
    let reset_counter: f64 = get_number(frame, "system.reset.counter");
    let reset_counter = reset_counter as usize;
    assert_eq!(reset_counter, 0);

    // The first frame is split into topology and coodinates.
    let frame_response = client.next_frame(&mut frames).await;
    assert_eq!(frame_response.frame_index, 1);

    // RESETTING
    client.reset().await;
    let frame_response = client.step_frame(&mut frames).await;
    let frame = frame_response.frame.as_ref().unwrap();
    assert_eq!(frame_response.frame_index, 2);
    // A reset should not affect the simulation counter
    assert!(!has_value(frame, "system.simulation,counter"));
    // The reset counter should be incremented, though.
    let reset_counter: f64 = get_number(frame, "system.reset.counter");
    let reset_counter = reset_counter as usize;
    assert_eq!(reset_counter, 1);

    // LOADING NEXT SIMULATION
    client.next().await;
    let frame_response = client.next_frame(&mut frames).await;
    let frame = frame_response.frame.as_ref().unwrap();
    // We loaded a new simulation so the frame index is 0 as we need to reset the accumulated
    // frame.
    assert_eq!(frame_response.frame_index, 0);
    // We loaded a new simulation so the simulation counter must be upgraded.
    let simulation_counter: f64 = get_number(frame, "system.simulation.counter");
    let simulation_counter = simulation_counter as usize;
    assert_eq!(simulation_counter, 1);
    // Loading a simulation resets the reset counter.
    let reset_counter: f64 = get_number(frame, "system.reset.counter");
    let reset_counter = reset_counter as usize;
    assert_eq!(reset_counter, 0);

    // LOADING A SIMULATION WITH playback/load
    client.load(2).await;
    let frame_response = client.next_frame(&mut frames).await;
    let frame = frame_response.frame.as_ref().unwrap();
    // We loaded a new simulation so the frame index is 0 as we need to reset the accumulated
    // frame.
    assert_eq!(frame_response.frame_index, 0);
    // We loaded a new simulation so the simulation counter must be upgraded.
    let simulation_counter: f64 = get_number(frame, "system.simulation.counter");
    let simulation_counter = simulation_counter as usize;
    assert_eq!(simulation_counter, 2);
    // Loading a simulation resets the reset counter.
    let reset_counter: f64 = get_number(frame, "system.reset.counter");
    let reset_counter = reset_counter as usize;
    assert_eq!(reset_counter, 0);

    // This prevents the test from hanging, but I do not know why.
    // TODO: fix the server so this step is not necessary.
    client.step().await;

    server.close().await.unwrap().unwrap();
}

async fn is_server_in_essd(name: &str) -> bool {
    let available_servers = find_servers(ESSD_DEFAULT_PORT, Duration::from_secs(1))
        .await
        .unwrap();
    info!("Servers from ESSD: {available_servers:?}");
    available_servers.values().any(|hub| hub.name() == name)
}

#[test(tokio::test)]
async fn test_essd_stop() {
    info!("TEST_ESSD_STOP");
    let (server, _client) = create_server_client_pair().await.unwrap();
    let server_name = server.name().to_owned();

    // At this point, the ESSD server should be running.
    info!("Requesting ESSD before closing the server.");
    let server_is_found = is_server_in_essd(&server_name).await;
    assert!(server_is_found);

    // We stop the server and we give it some time to finish gracefully.
    server.close().await.unwrap().unwrap();
    tokio::time::sleep(Duration::from_millis(500)).await;

    // Now, we should not receive anything from ESSD anymore.
    info!("Requesting ESSD after closing the server.");
    let server_is_found = is_server_in_essd(&server_name).await;
    assert!(!server_is_found);
}
