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

struct MyEguiApp {
    input_type: InputSelection,
    input_path: Option<String>,
    server_thread: Option<std::thread::JoinHandle<Result<(), AppError>>>,
}

impl Default for MyEguiApp {
    fn default() -> Self {
        Self {
            input_type: InputSelection::DefaultInput,
            input_path: None,
            server_thread: None,
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
        let button = egui::widgets::Button::new("Run!");
        if ui.add_enabled(ready, button).clicked() {
            self.start_server();
        }
    }

    fn start_server(&mut self) {
        let mut arguments = Cli::default();
        if let InputSelection::FileInput = self.input_type {
            arguments.input_xml_path = self.input_path.clone()
        };
        let runtime = Runtime::new().expect("Unable to create Runtime");
        let _enter = runtime.enter();
        let handle = std::thread::spawn(move || runtime.block_on(main_to_wrap(arguments)));
        self.server_thread = Some(handle)
    }

    fn is_running(&self) -> bool {
        let Some(ref thread_handle) = self.server_thread else {
            return false;
        };
        !thread_handle.is_finished()
    }

    fn is_idle(&self) -> bool {
        !self.is_running()
    }
}

impl eframe::App for MyEguiApp {
   fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Narupa server");
            if self.is_idle() {
                self.input_selection(ui);
                self.run_button(ui);
            } else {
                ui.label("Server is running.");
            }
        });
    }
}

#[derive(PartialEq)]
enum InputSelection {
    DefaultInput,
    FileInput,
}