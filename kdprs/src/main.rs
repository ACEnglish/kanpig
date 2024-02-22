extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;

mod cli;
mod kmer;
mod similarity;

use crate::cli::ArgParser;
use crate::kmer::seq_to_kmer;
use crate::similarity::{cosinesim, weighted_cosinesim};

fn main() {
    pretty_env_logger::formatted_timed_builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = ArgParser::parse();

    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let a = seq_to_kmer("ACATACAATACAACATACATACCATGGACACAGTA".into(), 3);
    let b = seq_to_kmer("ACATAGAATTAGACATACATACCATGGACACAGTA".into(), 3);
    println!("{:?}", a);
    println!("{:?}", b);
    println!("{:?}", cosinesim(&a, &b)); // 0.8952744
    println!("{:?}", weighted_cosinesim(&a, &b)); // 0.8952744
    println!("{:?}", args);

}
