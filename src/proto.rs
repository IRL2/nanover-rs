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
