use narupa_rs::simulation::{Simulation, TestSimulation};


fn main() {
    let mut simulation = TestSimulation::new();
    println!("Hoy!");
    println!("Hello, world!");
    simulation.step(10);
    println!("end")
}
