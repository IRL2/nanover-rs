use crate::simulation::{Simulation, ToFrameData};

pub struct ReplaySimulation {}

impl Simulation for ReplaySimulation {
    fn step(&mut self, steps: i32) {
        unimplemented!();
    }

    fn reset(&mut self) {
        unimplemented!();
    }
}

impl ToFrameData for ReplaySimulation {
    fn to_framedata(
        &self,
        with_velocity: bool,
        with_forces: bool,
    ) -> nanover_proto::trajectory::FrameData {
        unimplemented!();
    }

    fn to_topology_framedata(&self) -> nanover_proto::trajectory::FrameData {
        unimplemented!();
    }
}
