use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
    thread,
    time::{Duration, Instant},
};

use nanover_proto::{state_update::StateUpdate, trajectory::GetFrameResponse, Mergeable};

use super::{specific::TrackedSimulation, Configuration};
use crate::{
    broadcaster::{BroadcastSendError, Broadcaster},
    frame_broadcaster::FrameBroadcaster,
    recording::{RecordPair, ReplaySimulation},
    simulation::Simulation,
    state_broadcaster::StateBroadcaster,
};

const DEFAULT_DELAY: u128 = ((1.0 / 30.0) * 1_000_000.0) as u128;
const SERVER_LOCK_TOKEN: &str = "server";

pub struct TrackedReplaySimulation {
    pub simulation: ReplaySimulation,
    aggregated_state: StateUpdate,
    reset_counter: usize,
    simulation_counter_base: usize,
    last_original_simulation_counter_sent: Option<usize>,
    current_time: u128,
    reached_end_frames: bool,
    reached_end_state: bool,
}

impl TrackedReplaySimulation {
    pub fn new(simulation: ReplaySimulation, simulation_counter_base: usize) -> Self {
        let reset_counter = 0;
        let current_time = 0;
        let reached_end_frames = false;
        let reached_end_state = false;
        let last_original_simulation_counter_sent = None;
        let aggregated_state = StateUpdate::default();
        Self {
            simulation,
            aggregated_state,
            reset_counter,
            simulation_counter_base,
            last_original_simulation_counter_sent,
            current_time,
            reached_end_frames,
            reached_end_state,
        }
    }

    fn time_next_record(&self) -> Option<u128> {
        self.simulation.time_next_record()
    }

    fn time_current_frame(&self) -> Option<u128> {
        self.simulation.time_current_frame()
    }

    fn time_current_state(&self) -> Option<u128> {
        self.simulation.time_current_state()
    }

    pub fn read_next_frame(&mut self) {
        self.simulation.next_frame();
    }

    pub fn read_next_state(&mut self) {
        self.simulation.next_state();
    }

    fn last_frame_read(&self) -> Option<&RecordPair<GetFrameResponse>> {
        self.simulation.last_frame_read()
    }

    fn last_state_read(&self) -> Option<&RecordPair<StateUpdate>> {
        self.simulation.last_state_read()
    }

    pub fn simulation_loop_iteration(
        &mut self,
        sim_clone: &Arc<Mutex<FrameBroadcaster>>,
        state_clone: &Arc<Mutex<StateBroadcaster>>,
        now: &Instant,
        configuration: &Configuration,
    ) {
        state_clone
            .lock()
            .unwrap()
            .atomic_lock_updates(
                SERVER_LOCK_TOKEN.into(),
                BTreeMap::from([("scene".into(), Some(Duration::from_secs(1)))]),
            )
            .expect("Could not lock the scene");
        let reached_end_recording = self.reached_end_frames && self.reached_end_state;
        if reached_end_recording {
            thread::sleep(Duration::from_micros(DEFAULT_DELAY as u64));
            return;
        }

        if let Some(next_frame_time) = self.time_current_frame() {
            if self.current_time > next_frame_time {
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
                    self.reached_end_frames = true;
                } else {
                    self.read_next_frame();
                }
            }
        };

        if let Some(next_state_time) = self.time_current_state() {
            if self.current_time > next_state_time {
                self.send_state_update(state_clone.clone());
                if self
                    .last_state_read()
                    .as_ref()
                    .and_then(|pair| pair.next.as_ref())
                    .is_none()
                {
                    self.reached_end_state = true;
                } else {
                    self.read_next_state();
                }
            }
        }

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

    fn send_state_update(&mut self, state_clone: Arc<Mutex<StateBroadcaster>>) {
        let state_update = self
            .simulation
            .last_state_read()
            .cloned()
            .unwrap_or_default()
            .current
            .as_record();
        self.aggregated_state.merge(&state_update);
        state_clone
            .lock()
            .unwrap()
            .send(state_update)
            .expect("Cound not send state update");
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
        self.last_original_simulation_counter_sent
            .unwrap_or_default()
    }

    fn send_regular_frame(
        &mut self,
        sim_clone: Arc<Mutex<FrameBroadcaster>>,
        _with_velocities: bool,
        _with_forces: bool,
    ) -> Result<(), BroadcastSendError> {
        let pair = self.last_frame_read().cloned().unwrap_or_default();
        let frame_response = pair.current.as_record();
        let frame_index = frame_response.frame_index;
        let mut frame = frame_response.frame.unwrap_or_default();
        let mut simulation_counter: Option<usize> = None;
        frame
            .values
            .entry("system.simulation.counter".to_string())
            .and_modify(|counter| {
                let counter_base = self.simulation_counter_base as f64;
                if let Some(prost_types::value::Kind::NumberValue(ref mut counter_from_record)) =
                    counter.kind
                {
                    let simulation_counter_to_send = *counter_from_record + counter_base;
                    simulation_counter = Some(simulation_counter_to_send as usize);
                    *counter_from_record = simulation_counter_to_send;
                };
            });
        if simulation_counter.is_some() {
            self.last_original_simulation_counter_sent = simulation_counter;
        }

        if frame_index == 0 {
            sim_clone.lock().unwrap().send_reset_frame(frame)
        } else {
            sim_clone.lock().unwrap().send_frame(frame)
        }
    }

    fn clear_state(&mut self, state_clone: Arc<Mutex<StateBroadcaster>>) {
        let Some(changed_keys) = &self.aggregated_state.changed_keys else {
            return;
        };
        let keys_to_remove = BTreeMap::from_iter(changed_keys.fields.keys().map(|key| {
            (
                key.clone(),
                prost_types::Value {
                    kind: Some(prost_types::value::Kind::NullValue(0)),
                },
            )
        }));
        let changed_keys = Some(prost_types::Struct {
            fields: keys_to_remove,
        });
        let update = StateUpdate { changed_keys };
        {
            let mut locked_state = state_clone.lock().unwrap();
            locked_state
                .send(update)
                .expect("Could not send reset update.");
            locked_state
                .atomic_lock_updates(
                    SERVER_LOCK_TOKEN.into(),
                    BTreeMap::from([("scene".into(), None)]),
                )
                .expect("Could not release the lock.");
        }
        self.aggregated_state = StateUpdate::default();
    }
}
