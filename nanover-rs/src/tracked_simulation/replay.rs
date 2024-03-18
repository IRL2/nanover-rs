use std::{
    sync::{Arc, Mutex},
    thread,
    time::{Duration, Instant},
};

use nanover_proto::trajectory::GetFrameResponse;

use super::{specific::TrackedSimulation, Configuration};
use crate::{
    broadcaster::BroadcastSendError,
    frame_broadcaster::FrameBroadcaster,
    recording::{RecordPair, ReplaySimulation},
    simulation::{Simulation, ToFrameData},
};

const DEFAULT_DELAY: u128 = ((1.0 / 30.0) * 1_000_000.0) as u128;

pub struct TrackedReplaySimulation {
    pub simulation: ReplaySimulation,
    reset_counter: usize,
    simulation_counter: usize,
    current_time: u128,
    reached_end: bool,
}

impl TrackedReplaySimulation {
    pub fn new(simulation: ReplaySimulation, simulation_counter: usize) -> Self {
        let reset_counter = 0;
        let current_time = 0;
        let reached_end = false;
        Self {
            simulation,
            reset_counter,
            simulation_counter,
            current_time,
            reached_end,
        }
    }

    fn time_next_record(&self) -> Option<u128> {
        self.simulation.time_next_record()
    }

    fn time_current_record(&self) -> Option<u128> {
        self.simulation.time_current_record()
    }

    pub fn read_next_frame(&mut self) {
        self.simulation.next_frame();
    }

    fn last_frame_read(&self) -> Option<&RecordPair<GetFrameResponse>> {
        self.simulation.last_frame_read()
    }

    pub fn simulation_loop_iteration(
        &mut self,
        sim_clone: &Arc<Mutex<FrameBroadcaster>>,
        now: &Instant,
        configuration: &Configuration,
    ) {
        if self.reached_end {
            thread::sleep(Duration::from_micros(DEFAULT_DELAY as u64));
            return;
        }
        let Some(next_time) = self.time_current_record() else {
            thread::sleep(Duration::from_micros(DEFAULT_DELAY as u64));
            return;
        };

        if self.current_time > next_time {
            self.send_regular_frame(
                sim_clone.clone(),
                configuration.with_velocities,
                configuration.with_forces,
            )
            .expect("Cannot send replay frame");
            if self
                .last_frame_read()
                .as_ref()
                .and_then(|pair| pair.next.as_ref())
                .is_none()
            {
                self.reached_end = true;
            } else {
                self.read_next_frame();
            }
        };

        let Some(next_time) = self.time_next_record() else {
            // There are no more record to read so we just stall, but not too long because we still
            // need to deal with the playback orders.
            thread::sleep(Duration::from_micros(DEFAULT_DELAY as u64));
            return;
        };

        let delay_to_next_frame = next_time
            .checked_sub(self.current_time)
            .unwrap_or(DEFAULT_DELAY);
        let delay = delay_to_next_frame.min(DEFAULT_DELAY);
        let elapsed = now.elapsed();
        let time_left = Duration::from_micros(delay as u64).saturating_sub(elapsed);
        thread::sleep(time_left);
        self.current_time += delay;
    }
}

impl TrackedSimulation for TrackedReplaySimulation {
    fn step(&mut self) {
        self.simulation.step(1);
    }

    fn reset(&mut self) {
        self.simulation.reset();
        self.current_time = 0;
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
