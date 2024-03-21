use log::trace;
use std::{
    fs::File,
    io::{self, BufReader, Read, Seek},
    marker::PhantomData,
};
use thiserror::Error;

use nanover_proto::{
    state_update::StateUpdate,
    trajectory::{FrameData, GetFrameResponse},
    Mergeable,
};

use crate::{
    application::RecordingPath,
    simulation::{Simulation, ToFrameData},
};

const MAGIC_NUMBER: u64 = 6661355757386708963;
const FORMAT_VERSION: u64 = 2;

#[derive(Clone, Debug)]
pub struct TimedRecord<T> {
    record: T,
    timestamp: u128,
}

#[derive(Clone, Debug)]
pub struct RecordPair<T> {
    pub current: TimedRecord<T>,
    pub next: Option<TimedRecord<T>>,
}

impl<T> TimedRecord<T>
where
    T: Clone,
{
    pub fn record(&self) -> &T {
        &self.record
    }

    pub fn as_record(self) -> T {
        self.record
    }

    pub fn timestamp(&self) -> u128 {
        self.timestamp
    }
}

impl<T> Default for TimedRecord<T>
where
    T: Default,
{
    fn default() -> Self {
        Self {
            record: T::default(),
            timestamp: 0,
        }
    }
}

impl<T> Default for RecordPair<T>
where
    T: Default,
{
    fn default() -> Self {
        Self {
            current: TimedRecord::default(),
            next: None,
        }
    }
}

struct RecordingFile<T>
where
    T: Default,
{
    source: BufReader<File>,
    first_record_position: u64,
    last_read: RecordPair<T>,
    aggregate: T,
    _phantom: PhantomData<T>,
}

#[derive(Error, Debug)]
pub enum RecordingReadError {
    #[error("Error while reading the file.")]
    IOError(#[from] std::io::Error),
    #[error("The input file is not a NanoVer recording.")]
    NotARecording,
    #[error("File version {0} is not supported.")]
    UnsuportedFormatVersion(u64),
}

impl<T> RecordingFile<T>
where
    T: prost::Message + Default + Mergeable + Clone,
{
    fn try_new(source: File) -> Result<Self, RecordingReadError> {
        let mut buffered_source = BufReader::new(source);
        let magic_number = read_u64(&mut buffered_source)?;
        if magic_number != MAGIC_NUMBER {
            return Err(RecordingReadError::NotARecording);
        }
        let format_version = read_u64(&mut buffered_source)?;
        if format_version != FORMAT_VERSION {
            return Err(RecordingReadError::UnsuportedFormatVersion(format_version));
        }
        let first_record_position = buffered_source.stream_position()?;
        let first_read = Self::first_read(&mut buffered_source);
        Ok(Self {
            source: buffered_source,
            first_record_position,
            last_read: first_read,
            aggregate: T::default(),
            _phantom: PhantomData,
        })
    }

    fn first_read(source: &mut BufReader<File>) -> RecordPair<T> {
        let next_record = read_one_frame(source).ok();
        RecordPair {
            current: TimedRecord {
                record: T::default(),
                timestamp: 0,
            },
            next: next_record,
        }
    }

    fn current_record(&self) -> T {
        self.last_read.current.record.clone()
    }

    pub fn last_read(&self) -> &RecordPair<T> {
        &self.last_read
    }

    fn reset(&mut self) {
        self.source
            .seek(io::SeekFrom::Start(self.first_record_position))
            .expect("Seeking failed.");
        self.last_read = Self::first_read(&mut self.source);
        self.aggregate = T::default();
    }

    fn next_frame_pair(&mut self) -> std::io::Result<RecordPair<T>> {
        match &self.last_read.next {
            None => Ok(self.last_read.clone()),
            Some(record) => {
                Mergeable::merge(&mut self.aggregate, &record.record);
                let next = read_one_frame(&mut self.source).ok();
                let pair = RecordPair {
                    current: TimedRecord {
                        record: self.aggregate.clone(),
                        timestamp: record.timestamp(),
                    },
                    next,
                };
                self.last_read = pair.clone();
                Ok(pair)
            }
        }
    }
}

pub struct ReplaySimulation {
    frame_source: Option<RecordingFile<GetFrameResponse>>,
    state_source: Option<RecordingFile<StateUpdate>>,
}

impl TryFrom<&RecordingPath> for ReplaySimulation {
    type Error = RecordingReadError;

    fn try_from(value: &RecordingPath) -> Result<Self, Self::Error> {
        let trajectory_file = value.trajectory.as_ref().map(File::open).transpose()?;
        let state_file = value.state.as_ref().map(File::open).transpose()?;
        Self::try_new(trajectory_file, state_file)
    }
}

impl ReplaySimulation {
    pub fn try_new(
        frame_source: Option<File>,
        state_source: Option<File>,
    ) -> Result<Self, RecordingReadError> {
        let frame_source = frame_source
            .map(RecordingFile::<GetFrameResponse>::try_new)
            .transpose()?;
        let state_source = state_source
            .map(RecordingFile::<StateUpdate>::try_new)
            .transpose()?;
        Ok(Self {
            frame_source,
            state_source,
        })
    }

    pub fn last_frame_read(&self) -> Option<&RecordPair<GetFrameResponse>> {
        Some(self.frame_source.as_ref()?.last_read())
    }

    pub fn last_state_read(&self) -> Option<&RecordPair<StateUpdate>> {
        Some(self.state_source.as_ref()?.last_read())
    }

    pub fn time_next_frame(&self) -> Option<u128> {
        Some(
            self.frame_source
                .as_ref()?
                .last_read
                .next
                .as_ref()?
                .timestamp(),
        )
    }

    pub fn time_next_state(&self) -> Option<u128> {
        Some(
            self.state_source
                .as_ref()?
                .last_read
                .next
                .as_ref()?
                .timestamp(),
        )
    }

    pub fn time_next_record(&self) -> Option<u128> {
        match (self.time_next_frame(), self.time_next_state()) {
            (None, None) => None,
            (Some(frame), Some(state)) => {
                trace!("next: frame: {frame}, state: {state}");
                Some(frame.min(state))
            }
            (Some(frame), None) => Some(frame),
            (None, Some(state)) => Some(state),
        }
    }

    pub fn time_current_frame(&self) -> Option<u128> {
        Some(self.frame_source.as_ref()?.last_read.current.timestamp())
    }

    pub fn time_current_state(&self) -> Option<u128> {
        Some(self.state_source.as_ref()?.last_read.current.timestamp())
    }

    pub fn next_frame(&mut self) {
        if let Some(ref mut source) = self.frame_source {
            source.next_frame_pair().expect("Could not load next frame");
        };
    }

    pub fn next_state(&mut self) {
        if let Some(ref mut source) = self.state_source {
            source
                .next_frame_pair()
                .expect("Could not load next state update");
        };
    }
}

impl Simulation for ReplaySimulation {
    fn step(&mut self, steps: i32) {
        let steps = steps.max(0);
        for _ in 0..steps {
            self.next_frame();
        }
    }

    fn reset(&mut self) {
        if let Some(ref mut source) = self.frame_source {
            source.reset()
        };
        if let Some(ref mut source) = self.state_source {
            source.reset()
        };
    }
}

impl ToFrameData for ReplaySimulation {
    fn to_framedata(&self, _with_velocity: bool, _with_forces: bool) -> FrameData {
        self.to_topology_framedata()
    }

    fn to_topology_framedata(&self) -> FrameData {
        self.frame_source
            .as_ref()
            .and_then(|source| source.current_record().frame)
            .unwrap_or_default()
    }
}

fn read_u128(file: &mut BufReader<File>) -> io::Result<u128> {
    let mut buffer = [0u8; 16];
    file.read_exact(&mut buffer)?;
    Ok(u128::from_le_bytes(buffer))
}

pub(crate) fn read_u64(file: &mut BufReader<File>) -> io::Result<u64> {
    let mut buffer = [0u8; 8];
    file.read_exact(&mut buffer)?;
    Ok(u64::from_le_bytes(buffer))
}

pub(crate) fn read_one_frame<T>(file: &mut BufReader<File>) -> io::Result<TimedRecord<T>>
where
    T: prost::Message + Default,
{
    let timestamp = read_u128(file)?;
    let size = read_u64(file)?;
    let mut update_bytes = Vec::new();
    file.take(size).read_to_end(&mut update_bytes)?;
    let record = T::decode(update_bytes.as_slice())?;
    Ok(TimedRecord { timestamp, record })
}
