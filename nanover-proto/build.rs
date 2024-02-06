fn main() -> Result<(), Box<dyn std::error::Error>> {
    tonic_build::configure().build_server(true).compile(
        &[
            "proto/nanover/protocol/array.proto",
            // "proto/nanover/protocol/address.proto",
            "proto/nanover/protocol/trajectory/frame.proto",
            "proto/nanover/protocol/trajectory/get_frame.proto",
            "proto/nanover/protocol/trajectory/trajectory_service.proto",
            "proto/nanover/protocol/state/state_service.proto",
            // "proto/nanover/protocol/instance/connect_to_trajectory.proto",
            // "proto/nanover/protocol/instance/instance_service.proto",
            // "proto/nanover/protocol/instance/representation_service.proto",
            "proto/nanover/protocol/command/command_service.proto",
        ],
        &["proto"],
    )?;

    Ok(())
}
