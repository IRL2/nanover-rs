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

use crate::simulation::{Simulation, ToFrameData};

const MAGIC_NUMBER: u64 = 6661355757386708963;
const FORMAT_VERSION: u64 = 2;

#[derive(Clone)]
pub struct TimedRecord<T> {
    record: T,
    timestamp: u128,
}

#[derive(Clone)]
pub struct RecordPair<T> {
    pub current: TimedRecord<T>,
    pub next: Option<TimedRecord<T>>,
}

impl<T> TimedRecord<T>
where
    T: Clone,
{
    fn record(&self) -> T {
        self.record.clone()
    }

    fn timestamp(&self) -> u128 {
        self.timestamp
    }
}

#[derive(Clone)]
enum LastRead<T> {
    NoMoreToRead(Option<RecordPair<T>>),
    Read(RecordPair<T>),
}

struct RecordingFile<T> {
    source: BufReader<File>,
    first_record_position: u64,
    last_read: LastRead<T>,
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

    fn first_read(source: &mut BufReader<File>) -> LastRead<T> {
        let next_record = read_one_frame(source).ok();
        LastRead::Read(RecordPair {
            current: TimedRecord {
                record: T::default(),
                timestamp: 0,
            },
            next: next_record,
        })
    }

    fn current_record(&self) -> T {
        match &self.last_read {
            LastRead::Read(record) => record.current.record.clone(),
            LastRead::NoMoreToRead(Some(record)) => record.current.record.clone(),
            LastRead::NoMoreToRead(None) => T::default(),
        }
    }

    fn delay_to_next_record(&self) -> Option<u128> {
        unimplemented!()
    }

    fn reset(&mut self) {
        self.source
            .seek(io::SeekFrom::Start(self.first_record_position))
            .expect("Seeking failed.");
        self.last_read = Self::first_read(&mut self.source);
        self.aggregate = T::default();
    }

    fn next_frame_pair(&mut self) -> std::io::Result<LastRead<T>> {
        let last_read: TimedRecord<T> = match &self.last_read {
            LastRead::NoMoreToRead(record) => return Ok(LastRead::NoMoreToRead(record.clone())),
            LastRead::Read(pair) => match &pair.next {
                Some(record) => record.clone(),
                None => return Ok(LastRead::NoMoreToRead(Some(pair.clone()))),
            },
        };

        Mergeable::merge(&mut self.aggregate, &last_read.record());
        let next = read_one_frame(&mut self.source).ok();
        let pair = RecordPair {
            current: TimedRecord {
                record: self.aggregate.clone(),
                timestamp: last_read.timestamp(),
            },
            next,
        };

        self.last_read = LastRead::Read(pair.clone());

        Ok(LastRead::Read(pair))
    }

    fn seek(&mut self, time: u128) -> std::io::Result<LastRead<T>> {
        let start_time = match &self.last_read {
            LastRead::NoMoreToRead(Some(pair)) => {
                return Ok(LastRead::NoMoreToRead(Some(pair.clone())))
            }
            LastRead::NoMoreToRead(None) => {
                return Ok(LastRead::NoMoreToRead(None));
            }
            LastRead::Read(pair) => pair.current.timestamp(),
        };
        if time < start_time {
            self.reset();
        }
        loop {
            match &self.last_read {
                LastRead::NoMoreToRead(_) => break,
                LastRead::Read(ref record) => match record.next {
                    Some(ref next) => {
                        if next.timestamp() > time {
                            break;
                        }
                    }
                    None => break,
                },
            };

            self.next_frame_pair()?;
        }
        Ok(self.last_read.clone())
    }
}

pub struct ReplaySimulation {
    frame_source: Option<RecordingFile<GetFrameResponse>>,
    state_source: Option<RecordingFile<StateUpdate>>,
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

    fn seek_frame(&mut self, time: u128) -> std::io::Result<Option<TimedRecord<GetFrameResponse>>> {
        let Some(ref mut source) = self.frame_source else {
            return Ok(None);
        };
        match source.seek(time)? {
            LastRead::NoMoreToRead(None) => Ok(None),
            LastRead::NoMoreToRead(Some(pair)) => Ok(Some(pair.current)),
            LastRead::Read(pair) => Ok(Some(pair.current)),
        }
    }

    fn seek_state(&mut self, time: u128) -> std::io::Result<Option<TimedRecord<StateUpdate>>> {
        let Some(ref mut source) = self.state_source else {
            return Ok(None);
        };
        match source.seek(time)? {
            LastRead::NoMoreToRead(None) => Ok(None),
            LastRead::NoMoreToRead(Some(pair)) => Ok(Some(pair.current)),
            LastRead::Read(pair) => Ok(Some(pair.current)),
        }
    }

    pub fn seek(
        &mut self,
        time: u128,
    ) -> std::io::Result<(
        Option<TimedRecord<GetFrameResponse>>,
        Option<TimedRecord<StateUpdate>>,
    )> {
        let frame = self.seek_frame(time)?;
        let state = self.seek_state(time)?;
        Ok((frame, state))
    }

    pub fn delay_to_next_frame(&self) -> Option<u128> {
        self.frame_source.as_ref()?.delay_to_next_record()
    }

    pub fn next_frame(&mut self) {
        if let Some(ref mut source) = self.frame_source {
            source.next_frame_pair().expect("Could not load next frame");
        };
    }
}

impl Simulation for ReplaySimulation {
    fn step(&mut self, _steps: i32) {
        unimplemented!();
    }

    fn reset(&mut self) {
        if let Some(ref mut source) = self.frame_source {
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
