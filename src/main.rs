use narupa_rs::simulation::{Simulation, TestSimulation};
use narupa_rs::frame::FrameData;

fn main() {
    let mut frame = FrameData::empty();
    frame.insert_number_value("particle.count", 10.0);


    let mut simulation = TestSimulation::new();
    println!("Hoy!");
    println!("Hello, world!");
    simulation.step(10);
    println!("end")
}
