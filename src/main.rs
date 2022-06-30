use narupa_rs::frame::FrameData;
use narupa_rs::simulation::{Simulation, ToFrameData, TestSimulation};
use narupa_rs::proto::protocol::trajectory::trajectory_service_server::TrajectoryServiceServer;
use narupa_rs::trajectory::Trajectory;
use std::{time, thread};
use std::time::Duration;
use tonic::transport::Server;
use std::net::ToSocketAddrs;
use std::sync::{Arc, Mutex};


#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(empty_frame));
    let sim_clone = Arc::clone(&frame_source);
    tokio::task::spawn_blocking(move || {
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
    println!("Let's go!");
    let server = Trajectory::new(Arc::clone(&frame_source));
    Server::builder()
        .add_service(TrajectoryServiceServer::new(server))
        .serve("[::1]:50051".to_socket_addrs().unwrap().next().unwrap())
        .await
        .unwrap();
    Ok(())
}