extern crate pretty_env_logger;

#[macro_use]
extern crate log;

mod kplib;
use clap::Parser;
use kplib::{genotyper_main, plup_main, Cli, Commands};

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Gt(args) => genotyper_main(args),
        Commands::Plup(args) => plup_main(args),
    }
}
