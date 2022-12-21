use tokio::sync::mpsc::Sender;
use crate::services::commands::Command;
use crate::proto::protocol::command::{CommandMessage, CommandReply};
use prost_types::Struct;

#[derive(Debug, Clone, Copy)]
pub enum PlaybackOrder {
    Play,
    Pause,
    Reset,
    Step,
}

pub struct PlaybackState {
    playing: bool,
}

impl PlaybackState {
    pub fn new(start_playing: bool) -> Self {
        PlaybackState { playing: start_playing }
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
            PlaybackOrder::Reset => {},
        }
    }
}

pub struct PlaybackCommand {
    channel: Sender<PlaybackOrder>,
    order: PlaybackOrder
}

impl PlaybackCommand {
    pub fn new(channel: Sender<PlaybackOrder>, order: PlaybackOrder) -> Self {
        Self {channel, order}
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