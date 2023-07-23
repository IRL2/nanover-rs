use std::collections::BTreeMap;

use narupa_proto::command::{CommandMessage, CommandReply};
use crate::services::commands::Command;
use prost_types::{Struct, Value, value::Kind};
use tokio::sync::mpsc::Sender;
use pack_prost::{UnPack, ToProstValue};

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

        self.channel.try_send(PlaybackOrder::Load(simulation_index)).unwrap();
        // TODO: indicate to the client that the command is valid. We do not
        // know if it succeeded, though.
        CommandReply { result: None }
    }

    fn arguments(&self) -> Option<Struct> {
        let arguments = BTreeMap::from([("index".into(), Value { kind: Some(Kind::NullValue(0)) })]);
        Some(Struct { fields: arguments })
    }
}

pub struct ListSimulations {
    simulations: Vec<String>,
}

impl ListSimulations {
    pub fn new(simulations: Vec<String>) -> Self {
        Self { simulations }
    }
}

impl Command for ListSimulations {
    fn run(&self, _input: CommandMessage) -> CommandReply {
        let simulation_list = self.simulations.to_prost_value();
        let result = BTreeMap::from([("simulations".into(), simulation_list)]);
        CommandReply { result: Some(Struct { fields: result }) }
    }

    fn arguments(&self) -> Option<Struct> {
        None
    }
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
