extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use noodles_vcf::{self as vcf};

mod bedparser;
mod chunker;
mod cli;
mod comparisons;
mod kmer;
mod regions;
mod similarity;

use crate::chunker::filter_entry;
use crate::cli::ArgParser;
use crate::kmer::seq_to_kmer;
use crate::regions::{build_region_tree, ContigMap};
use crate::similarity::{cosinesim, weighted_cosinesim};
use std::path::PathBuf;

/// Remove contigs from keep that aren't in reduce
/// Works in-place
fn check_vcf_contigs(keep: &mut ContigMap, reduce: &ContigMap) {
    let mut removed = vec![];
    for i in keep.keys() {
        if !reduce.contains_key(i) {
            removed.push(i.clone());
        }
    }
    if !removed.is_empty() {
        warn!(
            "{} --input contigs do not have a --vcf contig. Removing",
            removed.len()
        );
        for i in removed {
            keep.shift_remove(&i);
        }
    }
}

fn main() {
    pretty_env_logger::formatted_timed_builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = ArgParser::parse();

    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let mut input_vcf = vcf::reader::Builder::default()
        .build_from_path(args.io.input)
        .expect("Unable to parse vcf");
    let input_header = input_vcf.read_header().expect("Unable to parse header");
    let mut m_contigs = input_header.contigs().clone();

    let tree = build_region_tree(&m_contigs, args.io.regions);

    let mut m_input = chunker::VCFIter::new(input_vcf, input_header, tree, args.kd.clone());
    let mut cnt = 0;
    for entry in &mut m_input {
        cnt += 1;
    }
    println!("parsed {} entries", cnt);
}
