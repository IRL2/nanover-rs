# narups-rs --- A NanoVer-OpenMM server written in Rust

## Dependencies

### Linux

Compiling the project requires to install some system packages. The project requires a C and C++ compiler, the GTK3 header files, and the `protoc` protobuf compiler.

On Fedora, these requirements can be installed with the following commands:

```bash
sudo dnf groupinstall "Development Tools" "Development Libraries" "C Development Tools and Libraries"
sudo dnf install gtk3-devel protobuf-compiler protobuf-devel
```
