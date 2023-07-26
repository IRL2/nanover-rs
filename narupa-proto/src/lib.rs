pub mod protocol {
    tonic::include_proto!("narupa.protocol");
    pub mod trajectory {
        tonic::include_proto!("narupa.protocol.trajectory");
    }
    pub mod state {
        tonic::include_proto!("narupa.protocol.state");
    }
    pub mod command {
        tonic::include_proto!("narupa.protocol.command");
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