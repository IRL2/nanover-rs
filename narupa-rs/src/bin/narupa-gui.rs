use eframe::egui;
use log::{LevelFilter, SetLoggerError, debug, trace};
use narupa_proto::command::{command_client::CommandClient, GetCommandsRequest, CommandMessage};
use narupa_rs::application::{main_to_wrap, AppError, Cli, cancellation_channels, CancellationSenders};
use tonic::transport::Channel;
use std::{sync::Mutex, net::{SocketAddr, Ipv4Addr, IpAddr}, collections::BTreeMap, marker::PhantomData, str::FromStr};

fn main() {
    init_logging().expect("Could not setup logging.");
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "Narupa-RS server",
        native_options,
        Box::new(|cc| Box::new(MyEguiApp::new(cc))),
    );
}

static LOG_VECTOR: Mutex<Vec<(log::Level, String)>> = Mutex::new(Vec::new());
static UI_LOGGER: UILogger = UILogger;

struct UILogger;

impl log::Log for UILogger {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.level() <= log::Level::Trace
    }

    fn log(&self, record: &log::Record) {
        if self.enabled(record.metadata()) && record.target().starts_with("narupa") {
            let mut logs = LOG_VECTOR.lock().unwrap();
            logs.push((record.level(), record.args().to_string()));
        }
    }

    fn flush(&self) {
        let mut logs = LOG_VECTOR.lock().unwrap();
        logs.clear();
    }
}

fn init_logging() -> Result<(), SetLoggerError> {
    log::set_logger(&UI_LOGGER).map(|()| log::set_max_level(LevelFilter::Trace))
}

struct Client {
    command: CommandClient<Channel>,
}

impl Client {
    pub fn connect(address: &SocketAddr, runtime_handle: &tokio::runtime::Handle) -> Result<Self, tonic::transport::Error> {
        let endpoint = format!("http://{address}");
        let command = runtime_handle.block_on(CommandClient::connect(endpoint))?;
        Ok(Client { command })
    }

    pub fn get_command_list(&mut self, runtime_handle: &tokio::runtime::Handle) -> Result<Vec<String>, tonic::Status> {
        let request = GetCommandsRequest{};
        let reply = runtime_handle.block_on(self.command.get_commands(request))?.into_inner();
        Ok(reply.commands.into_iter().map(|message| message.name).collect())
    }

    pub fn run_command(&mut self, name: String, runtime_handle: &tokio::runtime::Handle) -> Result<(), tonic::Status> {
        let request = CommandMessage{ name, arguments: None };
        let _reply = runtime_handle.block_on(self.command.run_command(request))?;
        Ok(())
    }
}

struct Server {
    handle: tokio::task::JoinHandle<Result<(), AppError>>,
    cancel_tx: Option<CancellationSenders>,
}

impl Server {
    pub fn new(arguments: Cli, runtime_handle: &tokio::runtime::Handle) -> Self {
        let (cancel_tx, cancel_rx) = cancellation_channels();
        let handle =
            runtime_handle.spawn(main_to_wrap(arguments, cancel_rx));
        Server {
            handle,
            cancel_tx: Some(cancel_tx),
        }
    }

    pub fn is_running(&self) -> bool {
        !self.handle.is_finished()
    }

    pub fn stop(&mut self) {
        // We should be able to stop the server several time. Calling this
        // method means we want the server to be stopped, not that we want
        // that specific call to stop the server. Therefore, we can ignore
        // the send failing or the transmitter being None.
        if let Some(tx) = self.cancel_tx.take() {
            tx.send().ok();
        };
    }

    pub fn close(mut self) -> tokio::task::JoinHandle<Result<(), AppError>> {
        self.stop();
        self.handle
    }
}

struct NumericField<NumType> where NumType: FromStr {
    label: String,
    default: String,
    raw: String,
    _phantom: PhantomData<NumType>,
}

impl<NumType> NumericField<NumType> where NumType: FromStr {
    fn new(label: impl ToString, default: impl ToString) -> Self {
        Self {
            label: label.to_string(),
            default: default.to_string(),
            raw: default.to_string(),
            _phantom: PhantomData,
        }
    }

    fn convert(&self) -> Result<NumType, <NumType as FromStr>::Err> {
        self.raw.parse()
    }

    fn is_valid(&self) -> bool {
        self.convert().is_ok()
    }

    fn widget(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let text_color;
            if self.is_valid() {
                ui.label(&self.label);
                text_color = None;
            } else {
                ui.label(egui::RichText::new(&self.label).color(egui::Color32::RED));
                text_color = Some(egui::Color32::RED);
            }
            egui::TextEdit::singleline(&mut self.raw).text_color_opt(text_color).show(ui);
            if ui.button("Set to default").clicked() {
                self.raw = self.default.clone();
            }
        });
    }
}

struct MyEguiApp {
    runtime: tokio::runtime::Runtime,
    reference: Cli,
    input_type: InputSelection,
    input_path: Option<String>,
    server: Option<Server>,
    client: Option<Client>,
    error: Option<String>,
    log_level: LevelFilter,
    show_progression: bool,
    server_name: String,
    port: NumericField<u16>,
    simulation_fps: NumericField<f64>,
    frame_interval: NumericField<u32>,
    force_interval: NumericField<u32>,
    record_statistics: bool,
    statistics: Option<String>,
    statistics_fps: NumericField<f64>,
    record_trajectory: bool,
    trajectory: Option<String>,
    record_state: bool,
    state: Option<String>,
}

impl Default for MyEguiApp {
    fn default() -> Self {
        let reference = Cli::default();
        let runtime = tokio::runtime::Builder::new_multi_thread()
                .enable_all()
                .build()
                .unwrap();
        let _enter = runtime.enter();
        Self {
            runtime,
            reference: Cli::default(),
            input_type: InputSelection::DefaultInput,
            input_path: None,
            server: None,
            client: None,
            error: None,
            log_level: LevelFilter::Info,
            show_progression: reference.progression,
            server_name: reference.name,
            port: NumericField::new("Port", reference.port),
            simulation_fps: NumericField::new("Simulation FPS", reference.simulation_fps),
            frame_interval: NumericField::new("Frame interval", reference.frame_interval),
            force_interval: NumericField::new("Force interval", reference.force_interval),
            record_statistics: false,
            statistics: reference.statistics,
            statistics_fps: NumericField::new("Statistics FPS", reference.statistics_fps),
            record_trajectory: false,
            trajectory: None,
            record_state: false,
            state: None,
        }
    }
}

impl MyEguiApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        // Customize egui here with cc.egui_ctx.set_fonts and cc.egui_ctx.set_visuals.
        // Restore app state using cc.storage (requires the "persistence" feature).
        // Use the cc.gl (a glow::Context) to create graphics shaders and buffers that you can use
        // for e.g. egui::PaintCallback.
        Self::default()
    }

    fn input_selection(&mut self, ui: &mut egui::Ui) {
        let mut file_picked = if let Some(ref path) = self.input_path {
            path.clone()
        } else {
            "".to_owned()
        };
        ui.vertical(|ui| {
            ui.horizontal(|ui| {
                ui.radio_value(
                    &mut self.input_type,
                    InputSelection::DefaultInput,
                    "Demonstration input",
                )
            });
            ui.horizontal(|ui| {
                ui.radio_value(
                    &mut self.input_type,
                    InputSelection::FileInput,
                    "File input",
                );
                if ui.text_edit_singleline(&mut file_picked).changed() {
                    self.input_path = Some(file_picked);
                };
                if ui.button("Select files").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("Narupa XML", &["xml"])
                        .pick_file()
                    {
                        self.input_path = Some(path.display().to_string());
                        self.input_type = InputSelection::FileInput;
                    };
                };
            });
        });
    }

    fn run_button(&mut self, ui: &mut egui::Ui) {
        let ready = match (&self.input_type, &self.input_path) {
            (InputSelection::DefaultInput, _) => true,
            (InputSelection::FileInput, Some(_)) => true,
            _ => false,
        };
        let ready = ready
            && self.port.is_valid()
            && self.simulation_parameters_are_valid()
            && self.recording_paramaters_are_valid();
        let text = match &self.input_type {
            InputSelection::DefaultInput => "Run demonstration input!",
            InputSelection::FileInput => "Run the selected file!",
        };
        let button = egui::widgets::Button::new(text);
        if ui.add_enabled(ready, button).clicked() {
            self.start_server();
        }
    }

    fn stop_button(&mut self, ui: &mut egui::Ui) {
        if ui.button("Stop server").clicked() {
            self.stop_server();
        }
    }

    fn error_message(&mut self, ui: &mut egui::Ui) {
        self.collect_error();
        let Some(ref error) = self.error else {return};
        let error_text = egui::RichText::new(error)
            .color(egui::Color32::RED)
            .strong();
        ui.label(error_text);
    }

    fn verbosity_selector(&mut self, ui: &mut egui::Ui, progression: bool) {
        egui::CollapsingHeader::new("Verbosity")
            .show(ui, |ui| {
                ui.vertical(|ui| {
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut self.log_level, LevelFilter::Info, "Normal");
                        ui.radio_value(&mut self.log_level, LevelFilter::Debug, "Verbose");
                        ui.radio_value(&mut self.log_level, LevelFilter::Trace, "Super verbose");
                    });
                    if progression {
                        ui.checkbox(&mut self.show_progression, "Show simulation progression");
                    };
                })
            });
    }

    fn network_parameters(&mut self, ui: &mut egui::Ui) {
        let mut header = egui::RichText::new("Network");
        let port_is_valid = self.port.is_valid();
        if !port_is_valid {
            header = header.color(egui::Color32::RED);
        }
        egui::CollapsingHeader::new(header)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Server name");
                    ui.text_edit_singleline(&mut self.server_name);
                    if ui.button("Set to default").clicked() {
                        self.server_name = self.reference.name.clone();
                    }
                });
                self.port.widget(ui);
            });
    }

    fn simulation_parameters(&mut self, ui: &mut egui::Ui) {
        let mut header = egui::RichText::new("Simulation");
        if !self.simulation_parameters_are_valid() {
            header = header.color(egui::Color32::RED);
        };
        egui::CollapsingHeader::new(header)
            .show(ui, |ui| {
                self.simulation_fps.widget(ui);
                self.frame_interval.widget(ui);
                self.force_interval.widget(ui);
            });
    }

    fn recording_paramaters(&mut self, ui: &mut egui::Ui) {
        let mut statistics_path = if let Some(ref path) = self.statistics {
            path.clone()
        } else {
            "".to_string()
        };
        let mut trajectory_path = if let Some(ref path) = self.trajectory {
            path.clone()
        } else {
            "".to_string()
        };
        let mut state_path = if let Some(ref path) = self.state {
            path.clone()
        } else {
            "".to_string()
        };

        let mut header = egui::RichText::new("Recording");
        if !self.recording_paramaters_are_valid() {
            header = header.color(egui::Color32::RED);
        };
        egui::CollapsingHeader::new(header)
            .show(ui, |ui| {
                ui.checkbox(&mut self.record_statistics, "Record statistics");
                ui.horizontal(|ui| {
                    let text_field = egui::TextEdit::singleline(&mut statistics_path);
                    let label = if self.statistics_is_valid() {
                        egui::Label::new("Statistics file:")
                    } else {
                        egui::Label::new(egui::RichText::new("Statistics file:").color(egui::Color32::RED))
                    };
                    ui.add_enabled(self.record_statistics, label);
                    if ui.add_enabled(self.record_statistics, text_field).changed() {
                        self.statistics = Some(statistics_path);
                    };
                    if ui.add(egui::Button::new("Select file")).clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .save_file()
                        {
                            self.statistics = Some(path.display().to_string());
                            self.record_statistics = true;
                        };
                    }
                });
                self.statistics_fps.widget(ui);
                ui.checkbox(&mut self.record_trajectory, "Record trajectory");
                ui.horizontal(|ui| {
                    let text_field = egui::TextEdit::singleline(&mut trajectory_path);
                    let label = if self.trajectory_is_valid() {
                        egui::Label::new("Trajectory file:")
                    } else {
                        egui::Label::new(egui::RichText::new("Trajectory file:").color(egui::Color32::RED))
                    };
                    ui.add_enabled(self.record_trajectory, label);
                    if ui.add_enabled(self.record_trajectory, text_field).changed() {
                        self.trajectory = Some(trajectory_path);
                    };
                    if ui.add(egui::Button::new("Select file")).clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .save_file()
                        {
                            self.trajectory = Some(path.display().to_string());
                            self.record_trajectory = true;
                        };
                    }
                });
                ui.checkbox(&mut self.record_state, "Record shared state");
                ui.horizontal(|ui| {
                    let text_field = egui::TextEdit::singleline(&mut state_path);
                    let label = if self.state_is_valid() {
                        egui::Label::new("Shared state file:")
                    } else {
                        egui::Label::new(egui::RichText::new("Shared state file:").color(egui::Color32::RED))
                    };
                    ui.add_enabled(self.record_state, label);
                    if ui.add_enabled(self.record_state, text_field).changed() {
                        self.state = Some(state_path);
                    };
                    if ui.add(egui::Button::new("Select file")).clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .save_file()
                        {
                            self.state = Some(path.display().to_string());
                            self.record_state = true;
                        };
                    }
                });
            });
    }

    fn command_buttons(&mut self, ui: &mut egui::Ui) {
        let known_commands = BTreeMap::from([
            ("playback/play", "Play"),
            ("playback/pause", "Pause"),
            ("playback/step", "Step"),
            ("playback/reset", "Reset"),
            ("multiuser/radially-orient-origins", "Radially orient origins"),
        ]);
        egui::CollapsingHeader::new("Commands").show(ui, |ui| {
            let Some(commands) = self.get_command_list() else {
                return;
            };
            if commands.is_empty() {
                return;
            };
            commands.chunks(4).for_each(|row| {
                ui.horizontal(|ui| {
                    row.into_iter().for_each(|command| {
                        let button_label = *known_commands.get(command.as_str()).unwrap_or(&command.as_str());
                        if ui.button(button_label).clicked() {
                            self.run_client_command(command.to_string());
                        }});
                });
            });
        });
    }

    fn log_window(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical()
            .stick_to_bottom(true)
            .auto_shrink([true, true])
            .max_height(ui.available_height())
            .show(ui, |ui| {
                let logs = LOG_VECTOR.lock().unwrap();
                logs.iter().for_each(|(level, message)| {
                    if level <= &self.log_level {
                        ui.label(format!("[{level}] {message}"));
                    }
                });
            });
    }
}

impl MyEguiApp {
    fn simulation_parameters_are_valid(&self) -> bool {
        self.simulation_fps.is_valid()
        && self.frame_interval.is_valid()
        && self.force_interval.is_valid()
    }

    fn statistics_is_valid(&self) -> bool {
        if self.record_statistics {
            let Some(ref path) = self.statistics else {
                return false;
            };
            !path.is_empty()
        } else {
            true
        }
    }

    fn trajectory_is_valid(&self) -> bool {
        if self.record_trajectory {
            let Some(ref path) = self.trajectory else {
                return false;
            };
            !path.is_empty()
        } else {
            true
        }
    }

    fn state_is_valid(&self) -> bool {
        if self.record_state {
            let Some(ref path) = self.state else {
                return false;
            };
            !path.is_empty()
        } else {
            true
        }
    }

    fn recording_paramaters_are_valid(&self) -> bool {
        (!self.record_statistics || (self.statistics_fps.is_valid() && self.statistics_is_valid()))
        && (!self.record_trajectory || self.trajectory_is_valid())
        && (!self.record_state || self.state_is_valid())
    }
}

impl MyEguiApp {
    fn build_arguments(&self) -> Result<Cli, ()> {
        let port = self.port.convert().map_err(|_| ())?;
        let simulation_fps = self.simulation_fps.convert().map_err(|_| ())?;
        let frame_interval = self.frame_interval.convert().map_err(|_| ())?;
        let force_interval = self.force_interval.convert().map_err(|_| ())?;

        let mut arguments = Cli::default();
        arguments.progression = self.show_progression;
        arguments.port = port;
        arguments.name = self.server_name.clone();
        if let InputSelection::FileInput = self.input_type {
            arguments.input_xml_path = self
                .input_path
                .as_ref()
                .map(|p| vec![p.clone()])
                .unwrap_or_else(|| Vec::new());
        };
        arguments.simulation_fps = simulation_fps;
        arguments.frame_interval = frame_interval;
        arguments.force_interval = force_interval;

        if self.record_statistics {
            let Some(ref statistics) = self.statistics else {
                return Err(());
            };
            arguments.statistics = Some(statistics.clone());
        } else {
            arguments.statistics = None;
        };

        if self.record_trajectory {
            let Some(ref trajectory) = self.trajectory else {
                return Err(());
            };
            arguments.trajectory = Some(trajectory.clone());
        } else {
            arguments.trajectory = None;
        }

        Ok(arguments)
    }

    fn start_server(&mut self) {
        log::logger().flush();
        self.clear_error();

        let Ok(arguments) = self.build_arguments() else {
            self.error = Some("An invalid argument was provided. Check red fields.".to_string());
            return;
        };

        let runtime_handle = self.runtime.handle();
        let server = Server::new(arguments, runtime_handle);
        self.server = Some(server);
    }

    fn stop_server(&mut self) {
        self.client = None;
        if let Some(mut server) = self.server.take() {
            server.stop();
        };
    }

    fn is_running(&self) -> bool {
        let Some(ref server) = self.server else {
            return false;
        };
        server.is_running()
    }

    fn is_idle(&self) -> bool {
        !self.is_running()
    }

    fn server_has_issue(&self) -> bool {
        self.server.is_some() && self.is_idle()
    }

    fn try_connect_client(&mut self) {
        if self.is_running() && self.client.is_none() {
            let address = SocketAddr::new(IpAddr::from(Ipv4Addr::new(127, 0, 0, 1)), self.port.convert().unwrap());
            let client = Client::connect(&address, self.runtime.handle());
            match client {
                Ok(client) => {
                    self.client = Some(client);
                    debug!("Client connected. Some features will be activated.");
                }
                Err(status) => {
                    trace!("Cannot connect a client: {:?}", status);
                }
            }
        }
    }

    fn get_command_list(&mut self) -> Option<Vec<String>> {
        if self.client.is_none() {
            self.try_connect_client();
        };
        let Some(ref mut client) = self.client else {
            return None;
        };
        match client.get_command_list(self.runtime.handle()) {
            Ok(commands) => {
                trace!("Commands: {commands:?}");
                Some(commands)
            }
            Err(status) => {
                trace!("Could not get the list of commands: {:?}", status);
                None
            }
        }
    }

    fn run_client_command(&mut self, name: String) {
        if self.client.is_none() {
            self.try_connect_client();
        };
        let Some(ref mut client) = self.client else {
            return;
        };
        match client.run_command(name.clone(), self.runtime.handle()) {
            Ok(_) => debug!("Successfully ran {name}"),
            Err(status) => debug!("Error while running {name}: {status:?}."),
        };
    }
}

impl MyEguiApp {
    fn collect_error(&mut self) {
        if self.server_has_issue() {
            let Some(server) = self.server.take() else {return};
            let Err(error) = self.runtime.block_on(server.close()).unwrap() else {return};
            self.error = Some(format!("{error}"));
        }
    }

    fn clear_error(&mut self) {
        self.error = None;
    }
}

impl eframe::App for MyEguiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Narupa server");
            self.error_message(ui);
            if self.is_idle() {
                self.input_selection(ui);
                self.verbosity_selector(ui, true);
                self.network_parameters(ui);
                self.simulation_parameters(ui);
                self.recording_paramaters(ui);
                self.run_button(ui);
            } else {
                ui.label("Server is running.");
                self.verbosity_selector(ui, false);
                self.command_buttons(ui);
                self.stop_button(ui);
                self.log_window(ui);
                ctx.request_repaint();
            }
        });
    }
}

#[derive(PartialEq)]
enum InputSelection {
    DefaultInput,
    FileInput,
}
