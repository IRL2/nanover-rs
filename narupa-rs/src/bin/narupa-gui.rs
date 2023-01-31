use eframe::egui;
use log::{LevelFilter, SetLoggerError};
use narupa_rs::application::{main_to_wrap, AppError, Cli};
use std::{sync::Mutex, num::{ParseIntError, ParseFloatError}};
use tokio::runtime::Runtime;

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
        metadata.level() <= log::Level::Debug
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
    log::set_logger(&UI_LOGGER).map(|()| log::set_max_level(LevelFilter::Debug))
}

struct Server {
    thread: std::thread::JoinHandle<Result<(), AppError>>,
    cancel_tx: Option<tokio::sync::oneshot::Sender<()>>,
}

impl Server {
    pub fn new(arguments: Cli) -> Self {
        let (cancel_tx, cancel_rx) = tokio::sync::oneshot::channel();
        let runtime = Runtime::new().expect("Unable to create Runtime");
        let _enter = runtime.enter();
        let handle =
            std::thread::spawn(move || runtime.block_on(main_to_wrap(arguments, cancel_rx)));
        Server {
            thread: handle,
            cancel_tx: Some(cancel_tx),
        }
    }

    pub fn is_running(&self) -> bool {
        !self.thread.is_finished()
    }

    pub fn stop(&mut self) {
        // We should be able to stop the server several time. Calling this
        // method means we want the server to be stopped, not that we want
        // that specific call to stop the server. Therefore, we can ignore
        // the send failing or the transmitter being None.
        if let Some(tx) = self.cancel_tx.take() {
            tx.send(()).ok();
        };
    }

    pub fn close(mut self) -> Result<(), AppError> {
        self.stop();
        self.thread.join().unwrap()
    }
}

struct MyEguiApp {
    reference: Cli,
    input_type: InputSelection,
    input_path: Option<String>,
    server: Option<Server>,
    error: Option<String>,
    log_level: LevelFilter,
    show_progression: bool,
    server_name: String,
    port: String,
    simulation_fps: String,
    frame_interval: String,
    force_interval: String,
    record_statistics: bool,
    statistics: Option<String>,
    statistics_fps: String,
    record_trajectory: bool,
    trajectory: Option<String>,
}

impl Default for MyEguiApp {
    fn default() -> Self {
        let reference = Cli::default();
        Self {
            reference: Cli::default(),
            input_type: InputSelection::DefaultInput,
            input_path: None,
            server: None,
            error: None,
            log_level: LevelFilter::Info,
            show_progression: reference.progression,
            server_name: reference.name,
            port: format!("{}", reference.port),
            simulation_fps: format!("{}", reference.simulation_fps),
            frame_interval: format!("{}", reference.frame_interval),
            force_interval: format!("{}", reference.force_interval),
            record_statistics: false,
            statistics: reference.statistics,
            statistics_fps: format!("{}", reference.statistics_fps),
            record_trajectory: false,
            trajectory: None,
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
            && self.port_is_valid()
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
        let port_is_valid = self.port_is_valid();
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
                ui.horizontal(|ui| {
                    if port_is_valid {
                        ui.label("Port");
                    } else {
                        ui.label(egui::RichText::new("Port").color(egui::Color32::RED));
                    }
                    let text_color = if port_is_valid {None} else {Some(egui::Color32::RED)};
                    egui::TextEdit::singleline(&mut self.port).text_color_opt(text_color).show(ui);
                    if ui.button("Set to default").clicked() {
                        self.port = self.reference.port.to_string();
                    }
                });
            });
    }

    fn simulation_parameters(&mut self, ui: &mut egui::Ui) {
        let mut header = egui::RichText::new("Simulation");
        if !self.simulation_parameters_are_valid() {
            header = header.color(egui::Color32::RED);
        };
        egui::CollapsingHeader::new(header)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    let text_color;
                    if self.simulation_fps_is_valid() {
                        ui.label("Simulation FPS");
                        text_color = None;
                    } else {
                        ui.label(egui::RichText::new("Simulation FPS").color(egui::Color32::RED));
                        text_color = Some(egui::Color32::RED);
                    }
                    egui::TextEdit::singleline(&mut self.simulation_fps).text_color_opt(text_color).show(ui);
                    if ui.button("Set to default").clicked() {
                        self.simulation_fps = self.reference.simulation_fps.to_string();
                    }
                });

                ui.horizontal(|ui| {
                    let text_color;
                    if self.frame_interval_is_valid() {
                        ui.label("Frame interval");
                        text_color = None;
                    } else {
                        ui.label(egui::RichText::new("Frame interval").color(egui::Color32::RED));
                        text_color = Some(egui::Color32::RED);
                    }
                    egui::TextEdit::singleline(&mut self.frame_interval).text_color_opt(text_color).show(ui);
                    if ui.button("Set to default").clicked() {
                        self.frame_interval = self.reference.frame_interval.to_string();
                    }
                });

                ui.horizontal(|ui| {
                    let text_color;
                    if self.force_interval_is_valid() {
                        ui.label("Force interval");
                        text_color = None;
                    } else {
                        ui.label(egui::RichText::new("Force interval").color(egui::Color32::RED));
                        text_color = Some(egui::Color32::RED);
                    }
                    egui::TextEdit::singleline(&mut self.force_interval).text_color_opt(text_color).show(ui);
                    if ui.button("Set to default").clicked() {
                        self.force_interval = self.reference.force_interval.to_string();
                    }
                })
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
                ui.horizontal(|ui| {
                    let text_color;
                    let label = if self.statistics_fps_is_valid() {
                        text_color = None;
                        egui::Label::new("Statistics FPS")
                    } else {
                        text_color = Some(egui::Color32::RED);
                        egui::Label::new(egui::RichText::new("Statistics FPS:").color(egui::Color32::RED))
                    };
                    let text_edit = egui::TextEdit::singleline(&mut self.statistics_fps).text_color_opt(text_color);
                    ui.add_enabled(self.record_statistics, label);
                    ui.add_enabled(self.record_statistics, text_edit);
                    if ui.add_enabled(self.record_statistics, egui::Button::new("Set to default")).clicked() {
                        self.statistics_fps = format!("{}", self.reference.statistics_fps);
                    }
                });
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
    fn convert_port(&self) -> Result<u16, ParseIntError> {
        self.port.parse()
    }

    fn port_is_valid(&self) -> bool {
        self.convert_port().is_ok()
    }

    fn convert_simulation_fps(&self) -> Result<f64, ParseFloatError> {
        self.simulation_fps.parse()
    }

    fn simulation_fps_is_valid(&self) -> bool {
        self.convert_simulation_fps().is_ok()
    }
    
    fn convert_frame_interval(&self) -> Result<u32, ParseIntError> {
        self.frame_interval.parse()
    }

    fn frame_interval_is_valid(&self) -> bool {
        self.convert_frame_interval().is_ok()
    }

    fn convert_force_interval(&self) -> Result<u32, ParseIntError> {
        self.force_interval.parse()
    }

    fn force_interval_is_valid(&self) -> bool {
        self.convert_force_interval().is_ok()
    }

    fn simulation_parameters_are_valid(&self) -> bool {
        self.simulation_fps_is_valid()
        && self.frame_interval_is_valid()
        && self.force_interval_is_valid()
    }

    fn convert_statistics_fps(&self) -> Result<f64, ParseFloatError> {
        self.statistics_fps.parse()
    }

    fn statistics_fps_is_valid(&self) -> bool {
        self.convert_statistics_fps().is_ok()
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
    fn recording_paramaters_are_valid(&self) -> bool {
        (!self.record_statistics || (self.statistics_fps_is_valid() && self.statistics_is_valid()))
        &&
        (!self.record_trajectory || self.trajectory_is_valid())
    }
}

impl MyEguiApp {
    fn build_arguments(&self) -> Result<Cli, ()> {
        let port = self.convert_port().map_err(|_| ())?;
        let simulation_fps = self.convert_simulation_fps().map_err(|_| ())?;
        let frame_interval = self.convert_frame_interval().map_err(|_| ())?;
        let force_interval = self.convert_force_interval().map_err(|_| ())?;

        let mut arguments = Cli::default();
        arguments.progression = self.show_progression;
        arguments.port = port;
        arguments.name = self.server_name.clone();
        if let InputSelection::FileInput = self.input_type {
            arguments.input_xml_path = self.input_path.clone()
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

        let server = Server::new(arguments);
        self.server = Some(server)
    }

    fn stop_server(&mut self) {
        let Some(mut server) = self.server.take() else {
            return;
        };
        server.stop();
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
}

impl MyEguiApp {
    fn collect_error(&mut self) {
        if self.server_has_issue() {
            let Some(server) = self.server.take() else {return};
            let Err(error) = server.close() else {return};
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
