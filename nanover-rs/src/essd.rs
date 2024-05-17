use log::trace;
use log::warn;
use network_interface::NetworkInterface;
use network_interface::NetworkInterfaceConfig;
use std::net::IpAddr;
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
    // port" flag set which needs a socket2 socket.
    let address = format!("0.0.0.0:0").parse::<SocketAddr>().unwrap();
    let socket = UdpSocket::bind(address).await.unwrap();
    let flag_socket = socket2::SockRef::from(&socket);
    flag_socket.set_nonblocking(true)?;
    flag_socket.set_broadcast(true)?;

    // These options are not available on windows
    #[cfg(not(target_os = "windows"))]
    {
        flag_socket.set_reuse_port(true)?;
        flag_socket.set_reuse_address(true)?;
    }

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
                // UDP broadcast does not work with IPv6
                if matches!(broadcast_address, IpAddr::V6(_)) {
                    continue;
                };
                let server_address = address.ip();
                let message = format!("{{\"name\": \"{name}\", \"address\": \"{server_address}\", \"port\": {port}, \"id\": \"{id}\", \"essd_version\": \"1.0.0\", \"services\": {{\"imd\": {port}, \"trajectory\": {port}, \"multiplayer\": {port}}}}}");
                let message = message.as_bytes();
                let result = socket
                    .send_to(message, format!("{broadcast_address}:54545"))
                    .await;
                if let Err(error) = result {
                    warn!("Failed to broadcast ESSD to {broadcast_address}:54545 ({error})")
                }
            }
        }
    }
    Ok(())
}
