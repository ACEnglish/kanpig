extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use chrono::Local;
use clap::Parser;
use kanpig::commands::{Commands, KanpigCommand};
use log::LevelFilter;
use pretty_env_logger::env_logger::fmt::{Color, Formatter, WriteStyle};
use std::io::Write;

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
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    pretty_env_logger::formatted_timed_builder()
        .filter_level(level)
        .format(|buf: &mut Formatter, record| {
            let ts = Local::now().format("%Y-%m-%dT%H:%M:%S");
            let mut level_style = buf.style();

            match record.level() {
                log::Level::Error => level_style.set_color(Color::Red).set_bold(true),
                log::Level::Warn => level_style.set_color(Color::Yellow).set_bold(true),
                log::Level::Info => level_style.set_color(Color::Green),
                log::Level::Debug => level_style.set_color(Color::Blue),
                log::Level::Trace => level_style.set_color(Color::Cyan),
            };

            let target = record.module_path().unwrap_or_default();

            writeln!(
                buf,
                " {} {:<5} {:<25} > {}",
                ts,
                level_style.value(record.level()),
                target,
                record.args()
            )
        })
        .write_style(WriteStyle::Auto)
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
