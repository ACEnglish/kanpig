extern crate pretty_env_logger;

#[macro_use]
extern crate log;

mod kplib;
use clap::Parser;
use kplib::{plup_main, genotyper_main, Cli, Commands, GTArgs, PlupArgs};

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Gt(args) => genotyper_main(args),
        Commands::Plup(args) => plup_main(args),
    }
}
