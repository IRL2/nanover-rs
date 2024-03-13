use std::{
    sync::{Arc, Mutex},
    time::Duration,
};

use nanover_proto::trajectory::FrameData;

use super::{openmm::TrackedOpenMMSimulation, replay::TrackedReplaySimulation};
use crate::{
    broadcaster::BroadcastSendError, frame_broadcaster::FrameBroadcaster,
    manifest::LoadedSimulation, simulation::ToFrameData,
};

pub struct Configuration {
    pub simulation_interval: Duration,
    pub frame_interval: usize,
    pub force_interval: usize,
    pub with_velocities: bool,
    pub with_forces: bool,
    pub auto_reset: bool,
    pub verbose: bool,
}

pub trait TrackedSimulation {
    fn step(&mut self);
    fn reset(&mut self);
    fn reset_counter(&self) -> usize;
    fn simulation_counter(&self) -> usize;
    fn send_regular_frame(
        &self,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        with_velocities: bool,
        with_forces: bool,
    ) -> Result<(), BroadcastSendError>;
}

pub enum SpecificSimulationTracked {
    OpenMM(TrackedOpenMMSimulation),
    Recording(TrackedReplaySimulation),
}

impl SpecificSimulationTracked {
    pub fn new(
        simulation: LoadedSimulation,
        simulation_counter: usize,
        configuration: &Configuration,
    ) -> Self {
        match simulation {
            LoadedSimulation::OpenMM(sim) => {
                let tracked = TrackedOpenMMSimulation::new(
                    sim,
                    simulation_counter,
                    configuration.frame_interval,
                );
                Self::OpenMM(tracked)
            }
            LoadedSimulation::Recording(sim) => {
                let tracked = TrackedReplaySimulation::new(sim, simulation_counter);
                Self::Recording(tracked)
            }
        }
    }

    pub fn step(&mut self) {
        match self {
            Self::OpenMM(simulation) => simulation.step(),
            Self::Recording(simulation) => simulation.step(),
        }
    }

    pub fn reset(&mut self) {
        match self {
            Self::OpenMM(simulation) => simulation.reset(),
            Self::Recording(simulation) => simulation.reset(),
        }
    }

    pub fn reset_counter(&self) -> usize {
        match self {
            Self::OpenMM(simulation) => simulation.reset_counter(),
            Self::Recording(simulation) => simulation.reset_counter(),
        }
    }

    pub fn simulation_counter(&self) -> usize {
        match self {
            Self::OpenMM(simulation) => simulation.simulation_counter(),
            Self::Recording(simulation) => simulation.simulation_counter(),
        }
    }

    pub fn send_regular_frame(
        &self,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        with_velocities: bool,
        with_forces: bool,
    ) -> Result<(), BroadcastSendError> {
        match self {
            Self::OpenMM(simulation) => {
                simulation.send_regular_frame(sim_clone, with_velocities, with_forces)
            }
            Self::Recording(simulation) => {
                simulation.send_regular_frame(sim_clone, with_velocities, with_forces)
            }
        }
    }
}

impl ToFrameData for SpecificSimulationTracked {
    fn to_framedata(&self, with_velocity: bool, with_forces: bool) -> FrameData {
        match self {
            Self::OpenMM(simulation) => simulation
                .simulation
                .to_framedata(with_velocity, with_forces),
            Self::Recording(simulation) => simulation
                .simulation
                .to_framedata(with_velocity, with_forces),
        }
    }

    fn to_topology_framedata(&self) -> FrameData {
        match self {
            Self::OpenMM(simulation) => simulation.simulation.to_topology_framedata(),
            Self::Recording(simulation) => simulation.simulation.to_topology_framedata(),
        }
    }
}

fn starting_frame(simulation: &SpecificSimulationTracked) -> FrameData {
    let mut frame = simulation.to_topology_framedata();
    frame
        .insert_number_value(
            "system.simulation.counter",
            simulation.simulation_counter() as f64,
        )
        .unwrap();
    frame
        .insert_number_value("system.reset.counter", simulation.reset_counter() as f64)
        .unwrap();
    frame
}

pub fn send_reset_frame(
    simulation: &SpecificSimulationTracked,
    sim_clone: Arc<Mutex<FrameBroadcaster>>,
) -> Result<(), BroadcastSendError> {
    let frame = starting_frame(simulation);
    sim_clone.lock().unwrap().send_reset_frame(frame)
}

pub fn next_simulation_counter(maybe_simulation: Option<&SpecificSimulationTracked>) -> usize {
    maybe_simulation
        .map(|simulation| simulation.simulation_counter() + 1)
        .unwrap_or_default()
}
