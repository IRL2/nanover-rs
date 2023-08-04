fn main() -> Result<(), Box<dyn std::error::Error>> {
    tonic_build::configure().build_server(true).compile(
        &[
            "proto/narupa/protocol/array.proto",
            // "proto/narupa/protocol/address.proto",
            "proto/narupa/protocol/trajectory/frame.proto",
            "proto/narupa/protocol/trajectory/get_frame.proto",
            "proto/narupa/protocol/trajectory/trajectory_service.proto",
            "proto/narupa/protocol/state/state_service.proto",
            // "proto/narupa/protocol/instance/connect_to_trajectory.proto",
            // "proto/narupa/protocol/instance/instance_service.proto",
            // "proto/narupa/protocol/instance/representation_service.proto",
            "proto/narupa/protocol/command/command_service.proto",
        ],
        &[
            "proto",
            "../../../Anaconda3/envs/narupa/Lib/site-packages/grpc_tools/_proto",
        ],
    )?;

    Ok(())
}
