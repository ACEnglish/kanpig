extern crate pretty_env_logger;

#[macro_use]
extern crate log;

mod genotype_main;
mod kplib;
mod plup_main;
use crate::{genotype_main::genotype_main, plup_main::plup_main};
use clap::Parser;
use kplib::{Cli, Commands, KanpigParams};

fn setup_logging(args: &(impl KanpigParams + std::fmt::Debug)) {
    let level = if args.trace() {
        log::LevelFilter::Trace
    } else if args.debug() {
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
        Commands::Gt(args) => {
            setup_logging(&args);
            genotype_main(args);
        }
        Commands::Plup(args) => {
            setup_logging(&args);
            plup_main(args)
        }
    };
}
