extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;

mod cli;

use crate::cli::ArgParser;

fn main() {
    pretty_env_logger::formatted_timed_builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = ArgParser::parse();

    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }
    println!("{:?}", args);
}
