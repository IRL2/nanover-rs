pub mod protocol {
    tonic::include_proto!("nanover.protocol");
    pub mod trajectory {
        tonic::include_proto!("nanover.protocol.trajectory");
    }
    pub mod state {
        tonic::include_proto!("nanover.protocol.state");
    }
    pub mod command {
        tonic::include_proto!("nanover.protocol.command");
    }
}

/// Mergeable data for a broadcaster.
pub trait Mergeable {
    /// Combine another instance with the current instance, in place.
    fn merge(&mut self, other: &Self);
}

pub mod frame;
pub mod state_update;
pub use protocol::*;
