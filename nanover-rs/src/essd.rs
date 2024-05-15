use log::trace;
use network_interface::NetworkInterface;
use network_interface::NetworkInterfaceConfig;
use std::net::SocketAddr;
use std::time::Duration;
use tokio::net::UdpSocket;
use tokio::time;

pub async fn serve_essd(
    name: String,
    port: u16,
    mut cancel_rx: tokio::sync::oneshot::Receiver<()>,
) -> std::io::Result<()> {
    let mut interval = time::interval(Duration::from_secs_f32(0.5));
    let id = uuid::Uuid::new_v4();

    // We need a tokio socket to work with async, that socket needs the "reuse
    // port" flag set which needs a socket2 socket. Socket2's socket are tricky
    // to build, standard library's socket are easier to create.
    let address = format!("0.0.0.0:{port}").parse::<SocketAddr>().unwrap();
    let std_socket = std::net::UdpSocket::bind(address)?;
    let flag_socket: socket2::Socket = std_socket.into();
    flag_socket.set_nonblocking(true)?;
    flag_socket.set_broadcast(true)?;

    // These options are not available on windows
    #[cfg(not(target_os = "windows"))]
    {
        flag_socket.set_reuse_port(true)?;
        flag_socket.set_reuse_address(true)?;
    }

    let socket = UdpSocket::bind("0.0.0.0:0").await.unwrap();
    socket.set_broadcast(true).unwrap();
    loop {
        match cancel_rx.try_recv() {
            Ok(_) | Err(tokio::sync::oneshot::error::TryRecvError::Closed) => {
                trace!("ESSD received cancellation order.");
                break;
            }
            Err(tokio::sync::oneshot::error::TryRecvError::Empty) => {}
        };

        interval.tick().await;
        let network_interfaces = NetworkInterface::show().unwrap();
        for interface in network_interfaces.iter() {
            for address in &interface.addr {
                let Some(broadcast_address) = address.broadcast() else {
                    continue;
                };
                let server_address = address.ip();
                let message = format!("{{\"name\": \"{name}\", \"address\": \"{server_address}\", \"port\": {port}, \"id\": \"{id}\", \"essd_version\": \"1.0.0\", \"services\": {{\"imd\": {port}, \"trajectory\": {port}, \"multiplayer\": {port}}}}}");
                let message = message.as_bytes();
                socket
                    .send_to(message, format!("{broadcast_address}:54545"))
                    .await
                    .unwrap();
            }
        }
    }
    Ok(())
}
