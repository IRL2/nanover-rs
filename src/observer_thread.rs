use std::sync::mpsc::Receiver;
use std::time::{Duration, Instant};
use std::fs::File;
use std::io::Write;
use std::thread;

use crate::broadcaster::BroadcasterSignal;

pub fn run_observer_thread(
        mut output_file: File,
        interval_microseconds: u64,
        frame_rx: Receiver<BroadcasterSignal>,
        simulation_rx: Receiver<usize>,
    ) {
    tokio::task::spawn_blocking(move || {
        let interval = Duration::from_millis(interval_microseconds);
        let mut previous = Instant::now();
        let mut keep_running = true;
        while keep_running {
            let now = Instant::now();
            let mut mean_fps = AverageAccumulator::new();
            loop {
                let frame_signal = frame_rx.try_recv();
                match frame_signal {
                    Ok(BroadcasterSignal::Send(instant)) => {
                        let frame_interval_in_seconds = instant.duration_since(previous).as_secs_f64();
                        let fps = 1.0 / frame_interval_in_seconds;
                        mean_fps.push(fps);
                        previous = instant;
                    },
                    Err(std::sync::mpsc::TryRecvError::Empty) => break,
                    Err(std::sync::mpsc::TryRecvError::Disconnected) => keep_running = false,
                    Ok(_) => (),
                };
            };
            let mut max_interactions = 0;
            loop {
                let num_interactions = simulation_rx.try_recv();
                match num_interactions {
                    Ok(num) => if num > max_interactions {max_interactions = num},
                    Err(std::sync::mpsc::TryRecvError::Empty) => break,
                    Err(std::sync::mpsc::TryRecvError::Disconnected) => keep_running = false,
                }
            }
            write!(output_file, "{:?}\t{:?}\t{}\n", now, mean_fps.average(), max_interactions).unwrap();
            let elapsed = now.elapsed();
            let time_left = match interval.checked_sub(elapsed) {
                Some(d) => d,
                None => Duration::from_millis(0),
            };
            thread::sleep(time_left);
        };
    });
}

struct AverageAccumulator {
    avg: f64,
    count: usize,
}

impl AverageAccumulator {
    pub fn new() -> Self {
        Self {avg: 0.0, count: 0}
    }

    pub fn push(&mut self, value: f64) {
        let n = self.count as f64;
        self.avg = (value + n * self.avg) / (n  + 1.0);
        self.count += 1;
    }

    pub fn average(&self) -> f64 {
        if self.count > 0 {self.avg} else {f64::NAN}
    }
}

impl Default for AverageAccumulator {
    fn default() -> Self {
        Self::new()    
    }
}