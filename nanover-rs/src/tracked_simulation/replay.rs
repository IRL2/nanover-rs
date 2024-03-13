use std::{
    sync::{Arc, Mutex},
    thread,
    time::{Duration, Instant},
};

use super::{specific::TrackedSimulation, Configuration};
use crate::{
    broadcaster::BroadcastSendError,
    frame_broadcaster::FrameBroadcaster,
    recording::ReplaySimulation,
    simulation::{Simulation, ToFrameData},
};

const DEFAULT_DELAY: f32 = 1.0 / 30.0;

pub struct TrackedReplaySimulation {
    pub simulation: ReplaySimulation,
    reset_counter: usize,
    simulation_counter: usize,
}

impl TrackedReplaySimulation {
    pub fn new(simulation: ReplaySimulation, simulation_counter: usize) -> Self {
        let reset_counter = 0;
        Self {
            simulation,
            reset_counter,
            simulation_counter,
        }
    }

    pub fn delay_to_next_frame(&self) -> Option<u128> {
        self.simulation.delay_to_next_frame()
    }

    pub fn read_next_frame(&mut self) {
        self.simulation.next_frame();
    }

    pub fn simulation_loop_iteration(
        &mut self,
        sim_clone: &Arc<Mutex<FrameBroadcaster>>,
        now: &Instant,
        configuration: &Configuration,
    ) {
        self.send_regular_frame(
            sim_clone.clone(),
            configuration.with_velocities,
            configuration.with_forces,
        )
        .expect("Cannot send replay frame");
        let delay_to_next_frame = self.delay_to_next_frame();
        self.read_next_frame();
        if let Some(delay) = delay_to_next_frame {
            let elapsed = now.elapsed();
            let time_left = match Duration::from_micros(delay as u64).checked_sub(elapsed) {
                Some(d) => d,
                None => Duration::from_micros(0),
            };
            thread::sleep(time_left);
        } else {
            thread::sleep(Duration::from_secs_f32(DEFAULT_DELAY));
        }
    }
}

impl TrackedSimulation for TrackedReplaySimulation {
    fn step(&mut self) {
        unimplemented!();
    }

    fn reset(&mut self) {
        self.simulation.reset();
        self.reset_counter += 1;
    }

    fn reset_counter(&self) -> usize {
        self.reset_counter
    }

    fn simulation_counter(&self) -> usize {
        self.simulation_counter
    }

    fn send_regular_frame(
        &self,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        with_velocities: bool,
        with_forces: bool,
    ) -> Result<(), BroadcastSendError> {
        let frame = self.simulation.to_framedata(with_velocities, with_forces);
        sim_clone.lock().unwrap().send_frame(frame)
    }
}
