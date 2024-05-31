# nanover-rs --- A NanoVer-OpenMM server written in Rust

## Build instructions

NanoVer-RS is written in Rust and therefore requires the Rust tool chain to build. Instructions to install the tool chain can be found at <https://www.rust-lang.org/tools/install>.

The program uses gRPC and needs the protobuf compiler, protoc, to compile.

OpenMM is compiled as part of building NanoVer-RS. It requires Python, cmake, and doxygen.

On linux, the GTK3 headers are needed to compile `nanover-gui`.

The command to build a release binary is `cargo build --release`. One would then need to copy the OpenMM library and plugin files from the path indicated during the compilation.

NanoVer-RS uses a modified version of [openmm-sys](https://docs.rs/crate/openmm-sys/latest) to bind to OpenMM. The documentation for this library describes how to use alternative builds of OpenMM, including the one from conda, for faster compilation.

## Runtime controls with environment variables

The `OPENMM_PLUGIN_DIR` environment variable allows to use OpenMM plugins located elsewhere than in a `lib` directory next to the executable.

NanoVer-RS uses [env_logger](https://docs.rs/env_logger/latest/env_logger/) to manage logs. This library uses the `RUST_LOG` environment variable to control the log level. For example, `RUST_LOG=trace` would set the log level to the most verbose.

## Dependencies

### Linux

Compiling the project requires to install some system packages. The project requires a C and C++ compiler, the GTK3 header files, and the `protoc` protobuf compiler.

On Fedora, these requirements can be installed with the following commands:

```bash
sudo dnf groupinstall "Development Tools" "Development Libraries" "C Development Tools and Libraries"
sudo dnf install gtk3-devel protobuf-compiler protobuf-devel
```
