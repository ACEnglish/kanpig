extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use noodles_vcf::{self as vcf};

mod cli;
mod kmer;
mod chunker;
mod regions;
mod similarity;
mod comparisons;

use crate::cli::ArgParser;
use crate::kmer::seq_to_kmer;
use crate::similarity::{cosinesim, weighted_cosinesim};
use crate::regions::{build_region_tree, ContigMap};

/// Remove contigs from keep that aren't in reduce
/// Works in-place
fn check_vcf_contigs(keep: &mut ContigMap, reduce: &ContigMap){
    let mut removed = vec![];
    for i in keep.keys() {
        if !reduce.contains_key(i) {
            removed.push(i.clone());
        }
    }
    if removed.len() != 0 {
        warn!("{} --input contigs do not have a --vcf contig. Removing", removed.len());
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

    let a = seq_to_kmer("ACATACAATACAACATACATACCATGGACACAGTA".into(), 3);
    let b = seq_to_kmer("ACATAGAATTAGACATACATACCATGGACACAGTA".into(), 3);
    println!("{:?}", a);
    println!("{:?}", b);
    println!("{:?}", cosinesim(&a, &b)); // 0.8952744
    println!("{:?}", weighted_cosinesim(&a, &b)); // 0.8952744
    println!("{:?}", args);

    // Must be an input
    let mut input_vcf = vcf::reader::Builder::default()
        .build_from_path(args.io.input)
        .expect("Unable to parse vcf");
    let input_header = input_vcf.read_header().expect("Unable to parse header");
    let mut m_contigs = input_header.contigs().clone();
    
    let mut base_vcf = match args.io.vcf {
        Some(i) => {
            let mut fh = vcf::reader::Builder::default().build_from_path(i).expect("Unable to parse vcf");
            let header = fh.read_header().expect("Unable to parse header");
            check_vcf_contigs(&mut m_contigs, header.contigs());
            Some(fh)
        },
        None => None
    };

    // Get all of its contigs, that's all we're parsing

    // Validate that all of the
    build_region_tree(&m_contigs, 1);
}
