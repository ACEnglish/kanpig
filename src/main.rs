extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use kanpig::commands::{Commands, KanpigCommand};

/// Entrypoint for the commands
#[derive(Parser, Clone, Debug)]
#[command(name = "kanpig")]
#[command(about = "Kmer ANalysis of PIleups for Genotyping")]
#[command(author = "ACEnglish", version)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

fn setup_logging(args: &impl KanpigCommand) {
    let level = if args.debug() {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };

    pretty_env_logger::formatted_timed_builder()
        .filter_level(level)
        .init();

    info!("params: {:#?}", args);
    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::GT(mut cmd) => {
            setup_logging(&cmd);
            cmd.run()
        }
        Commands::Plup(mut cmd) => {
            setup_logging(&cmd);
            cmd.run();
        }
        Commands::Trio(mut cmd) => {
            setup_logging(&cmd);
            cmd.run();
        }
        Commands::Mosaic(mut cmd) => {
            setup_logging(&cmd);
            cmd.run();
        }
    }
}
