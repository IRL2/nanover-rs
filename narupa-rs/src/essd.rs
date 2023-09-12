use network_interface::NetworkInterface;
use network_interface::NetworkInterfaceConfig;
use std::time::Duration;
use tokio::net::UdpSocket;
use tokio::time;

pub async fn serve_essd(name: String, port: u16) {
    let mut interval = time::interval(Duration::from_secs_f32(0.5));
    let id = uuid::Uuid::new_v4();

    let socket = UdpSocket::bind("0.0.0.0:0").await.unwrap();
    socket.set_broadcast(true).unwrap();
    loop {
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
}
