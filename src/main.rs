use narupa_rs::frame::FrameData;
use narupa_rs::simulation::{Simulation, ToFrameData, TestSimulation};
use narupa_rs::proto::protocol::trajectory::trajectory_service_server::TrajectoryServiceServer;
use narupa_rs::trajectory::Trajectory;
use tokio::time;
use std::time::Duration;
use tonic::transport::Server;
use std::net::ToSocketAddrs;
use tokio::sync::mpsc;
use std::sync::{Arc, Mutex};

/*
fn main() {
    let mut simulation = TestSimulation::new();
    println!("Hoy!");
    println!("Hello, world!");
    simulation.step(10);
    println!("{:?}", simulation.to_framedata());
}
*/

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let empty_frame = FrameData::empty();
    let frame_source = Arc::new(Mutex::new(empty_frame));
    let sim_clone = Arc::clone(&frame_source);
    tokio::task::spawn_blocking(move || {
        let mut simulation = TestSimulation::new();
        let mut interval = time::interval(Duration::from_millis(33));
        for i in 0.. {
            println!("{i}");
            simulation.step(10);
            let frame = simulation.to_framedata();
            let mut source = sim_clone.lock().unwrap();
            *source = frame;
            //println!("{:?}", frame);
            //interval.tick().await;
            /*
            match tx.blocking_send(Result::<_, ()>::Ok(frame)) {
                Ok(_) => {
                    // item (server response) was queued to be send to client
                }
                Err(_item) => {
                    // output_stream was build from rx and both are dropped
                    break;
                }
            }
            */
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
