use std::process::ExitCode;
use log::LevelFilter;
use env_logger::Builder;
use clap::Parser;
use tokio::runtime::Runtime;
use std::error::Error;
use log::{error, info};
use narupa_rs::application::{main_to_wrap, Cli, AppError, cancellation_channels};

fn main() -> ExitCode {
    let cli = Cli::parse();

    if std::env::var("RUST_LOG").is_ok() {
        env_logger::init();
    } else {
        let mut verbosity_level = LevelFilter::Info;
        if cli.verbose {
            verbosity_level = LevelFilter::Debug
        };
        if cli.trace {
            verbosity_level = LevelFilter::Trace
        };

        let mut builder = Builder::new();
        builder
            .filter_module("narupa_rs", verbosity_level)
            .filter_module("narupa_cli", verbosity_level)
            .format_target(false)
            .init();
    }

    let runtime = Runtime::new().expect("Unable to create Runtime");
    let _enter = runtime.enter();
    
    let (cancel_tx, cancel_rx) = cancellation_channels();
    tokio::spawn(async move {
        tokio::signal::ctrl_c().await.unwrap();
        // Your handler here
        info!("Closing the server. Goodbye!");
        cancel_tx.send().unwrap();
    });

    let run_status: Result<(), AppError> =
        std::thread::scope(|scope| scope.spawn(|| runtime.block_on(main_to_wrap(cli, cancel_rx))).join())
            .unwrap();
    let Err(ref error) = run_status else {
        return ExitCode::SUCCESS;
    };

    // The Display trait from tonic's errors is not very expressive for the
    // end user. We need to dig out the underlying error.
    let maybe_transport_error = if let AppError::TransportError(transport_error) = error {
        Some(transport_error)
    } else {
        None
    };
    //let maybe_transport_error = error.downcast_ref::<tonic::transport::Error>();
    let maybe_source_error = match maybe_transport_error {
        None => None,
        Some(transport_error) => transport_error.source(),
    };
    let maybe_hyper_error = match maybe_source_error {
        None => None,
        Some(source_error) => source_error.downcast_ref::<hyper::Error>(),
    };
    if let Some(hyper_error) = maybe_hyper_error {
        error!("{hyper_error}");
        return ExitCode::FAILURE;
    };

    error!("{error}");
    ExitCode::FAILURE
}
