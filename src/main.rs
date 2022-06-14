use narupa_rs::simulation::{Simulation, ToFrameData, TestSimulation};

fn main() {
    let mut simulation = TestSimulation::new();
    println!("Hoy!");
    println!("Hello, world!");
    simulation.step(10);
    println!("{:?}", simulation.to_framedata());
}
