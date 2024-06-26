use std::{collections::BTreeMap, path::PathBuf};

use crate::services::commands::Command;
use nanover_proto::command::{CommandMessage, CommandReply};
use pack_prost::{ToProstValue, UnPack};
use prost_types::{value::Kind, Struct, Value};
use tokio::sync::mpsc::Sender;

#[derive(Debug, Clone, Copy)]
pub enum PlaybackOrder {
    Play,
    Pause,
    Reset,
    Step,
    Load(usize),
    Next,
}

pub struct PlaybackState {
    playing: bool,
}

impl PlaybackState {
    pub fn new(start_playing: bool) -> Self {
        PlaybackState {
            playing: start_playing,
        }
    }

    pub fn is_playing(&self) -> bool {
        self.playing
    }

    pub fn update(&mut self, order: PlaybackOrder) {
        match order {
            PlaybackOrder::Play => self.playing = true,
            PlaybackOrder::Pause => self.playing = false,
            // Stepping implies to pause the trajectory
            PlaybackOrder::Step => self.playing = false,
            // Not our responsability here
            PlaybackOrder::Reset | PlaybackOrder::Load(_) | PlaybackOrder::Next => {}
        }
    }
}

pub struct PlaybackCommand {
    channel: Sender<PlaybackOrder>,
    order: PlaybackOrder,
}

impl PlaybackCommand {
    pub fn new(channel: Sender<PlaybackOrder>, order: PlaybackOrder) -> Self {
        Self { channel, order }
    }
}

impl Command for PlaybackCommand {
    fn run(&self, _input: CommandMessage) -> CommandReply {
        self.channel.try_send(self.order).unwrap();
        CommandReply { result: None }
    }

    fn arguments(&self) -> Option<Struct> {
        None
    }
}

pub struct LoadCommand {
    channel: Sender<PlaybackOrder>,
}

impl LoadCommand {
    pub fn new(channel: Sender<PlaybackOrder>) -> Self {
        Self { channel }
    }
}

impl Command for LoadCommand {
    fn run(&self, input: CommandMessage) -> CommandReply {
        let Some(argument_struct) = input.arguments else {
            // The client did not provide any arguments
            // TODO: return an error to the client
            return CommandReply { result: None };
        };
        let arguments = argument_struct.fields;
        let Some(simulation_index_value) = arguments.get("index") else {
            // The client did not provide the index of the simulation they want
            // to load.
            // TODO: return an error to the client
            return CommandReply { result: None };
        };
        let maybe_simulation_index: Option<f64> = simulation_index_value.unpack();
        let Some(maybe_simulation_index) = maybe_simulation_index else {
            // The simulation index provided by the user is not a number.
            // TODO: return an error to the client
            return CommandReply { result: None };
        };
        let maybe_simulation_index = maybe_simulation_index.trunc() as i64;
        let Ok(simulation_index) = TryInto::<usize>::try_into(maybe_simulation_index) else {
            // The simulation index is negative or too large.
            // TODO: return an error to the client
            return CommandReply { result: None };
        };

        self.channel
            .try_send(PlaybackOrder::Load(simulation_index))
            .unwrap();
        // TODO: indicate to the client that the command is valid. We do not
        // know if it succeeded, though.
        CommandReply { result: None }
    }

    fn arguments(&self) -> Option<Struct> {
        let arguments = BTreeMap::from([(
            "index".into(),
            Value {
                kind: Some(Kind::NullValue(0)),
            },
        )]);
        Some(Struct { fields: arguments })
    }
}

pub struct ListSimulations {
    simulations: Vec<String>,
}

impl ListSimulations {
    pub fn new(simulations: Vec<String>) -> Self {
        Self {
            simulations: reduce_paths(&simulations),
        }
    }
}

impl Command for ListSimulations {
    fn run(&self, _input: CommandMessage) -> CommandReply {
        let simulation_list = self.simulations.to_prost_value();
        let result = BTreeMap::from([("simulations".into(), simulation_list)]);
        CommandReply {
            result: Some(Struct { fields: result }),
        }
    }

    fn arguments(&self) -> Option<Struct> {
        None
    }
}

/// Remove the common part at the beginning of a series of paths.
fn reduce_paths(paths: &[String]) -> Vec<String> {
    // When there is only one path, it makes no sense to return the part that
    // differ amongst paths. Would we do that, we would return an empty path as
    // it is identical to itself. Instead, we return the file name if it can be
    // obtained, or the full path if something prevents us from getting a file
    // name (i.e. path is the root or ends with ..).
    if paths.len() == 1 {
        let path: PathBuf = paths[0].clone().into();
        return if let Some(file_name) = path.file_name() {
            vec![file_name.to_string_lossy().to_string()]
        } else {
            vec![paths[0].clone()]
        };
    }

    // TODO: If all the paths are identical, then we return a vector of empty
    // strings. We need to figure out what to do instead. The solution may also
    // solve the the edge case of a single path being provided.
    let mut iter = paths.iter();
    let Some(first) = iter.next() else {
        return Vec::new();
    };
    let reference: PathBuf = first.into();
    let mut reference: Vec<_> = reference.components().collect();
    for path in iter {
        let pathbuf = Into::<PathBuf>::into(path);
        let components = pathbuf.components();
        reference = reference
            .into_iter()
            .zip(components)
            .take_while(|(a, b)| *a == *b)
            .map(|(a, _b)| a)
            .collect();
    }
    let n_components_to_remove = reference.len();
    paths
        .iter()
        .map(|p| cleave_path(Into::<PathBuf>::into(p), n_components_to_remove))
        .collect()
}

fn cleave_path(path: PathBuf, n_components_to_remove: usize) -> String {
    let mut out = PathBuf::new();
    path.components()
        .skip(n_components_to_remove)
        .for_each(|component| out.push(component));
    out.to_string_lossy().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_playback_init(#[case] starting_play: bool) {
        let state = PlaybackState::new(starting_play);
        assert_eq!(state.is_playing(), starting_play);
    }

    #[rstest]
    #[case(true, PlaybackOrder::Play, true)]
    #[case(false, PlaybackOrder::Play, true)]
    #[case(true, PlaybackOrder::Pause, false)]
    #[case(false, PlaybackOrder::Pause, false)]
    #[case(true, PlaybackOrder::Reset, true)]
    #[case(false, PlaybackOrder::Reset, false)]
    #[case(true, PlaybackOrder::Step, false)]
    #[case(false, PlaybackOrder::Step, false)]
    fn test_playback_update(
        #[case] previous_play: bool,
        #[case] order: PlaybackOrder,
        #[case] expected_play: bool,
    ) {
        let mut state = PlaybackState::new(previous_play);
        state.update(order);
        assert_eq!(state.is_playing(), expected_play);
    }
}
