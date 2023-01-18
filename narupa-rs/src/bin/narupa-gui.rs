use eframe::egui;
use log::LevelFilter;
use tokio::runtime::Runtime;
use narupa_rs::application::{main_to_wrap, Cli, AppError};
use env_logger::Builder;

fn main() {
    let mut builder = Builder::new();
    builder
        .filter_module("narupa_rs", LevelFilter::Debug)
        .filter_module("narupa_gui", LevelFilter::Debug)
        .format_target(false)
        .init();
    let native_options = eframe::NativeOptions::default();
    eframe::run_native("Narupa-RS server", native_options, Box::new(|cc| Box::new(MyEguiApp::new(cc))));
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
        let handle = std::thread::spawn(move || runtime.block_on(main_to_wrap(arguments, cancel_rx)));
        Server {thread: handle, cancel_tx: Some(cancel_tx)}
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
    input_type: InputSelection,
    input_path: Option<String>,
    server: Option<Server>,
    error: Option<String>,
}

impl Default for MyEguiApp {
    fn default() -> Self {
        Self {
            input_type: InputSelection::DefaultInput,
            input_path: None,
            server: None,
            error: None,
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
        let mut file_picked = if let Some(ref path) = self.input_path {path.clone()} else {"".to_owned()};
        ui.vertical(|ui| {
            ui.horizontal(|ui| ui.radio_value(&mut self.input_type, InputSelection::DefaultInput, "Demonstration input"));
            ui.horizontal(|ui| {
                ui.radio_value(&mut self.input_type, InputSelection::FileInput, "File input");
                if ui.text_edit_singleline(&mut file_picked).changed() {
                    self.input_path = Some(file_picked);
                };
                if ui.button("Browse files").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("Narupa XML", &["xml"]).pick_file() {
                        self.input_path = Some(path.display().to_string());
                        println!("File picked: {:?}", path);
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
        let error_text = egui::RichText::new(error).color(egui::Color32::RED).strong();
        ui.label(error_text);
    }

    fn start_server(&mut self) {
        self.clear_error();
        let mut arguments = Cli::default();
        if let InputSelection::FileInput = self.input_type {
            arguments.input_xml_path = self.input_path.clone()
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
                self.run_button(ui);
            } else {
                ui.label("Server is running.");
                self.stop_button(ui);
            }
        });
    }
}

#[derive(PartialEq)]
enum InputSelection {
    DefaultInput,
    FileInput,
}