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
        _state_rx: Receiver<BroadcasterSignal>,
    ) {
    tokio::task::spawn_blocking(move || {
        let interval = Duration::from_millis(interval_microseconds);
        let start = Instant::now();
        let mut keep_running = true;
        while keep_running {
            let now = Instant::now();
            loop {
                let frame_signal = frame_rx.try_recv();
                match frame_signal {
                    Ok(BroadcasterSignal::Send(instant)) => {
                        write!(output_file, "{:?}\t{:?}\n", now, instant.duration_since(start)).unwrap();
                    },
                    Err(std::sync::mpsc::TryRecvError::Empty) => break,
                    Err(std::sync::mpsc::TryRecvError::Disconnected) => keep_running = false,
                    Ok(_) => (),
                };
            };
            let elapsed = now.elapsed();
            let time_left = match interval.checked_sub(elapsed) {
                Some(d) => d,
                None => Duration::from_millis(0),
            };
            thread::sleep(time_left);
        };
    });
}