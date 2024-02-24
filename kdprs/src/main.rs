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

use crate::chunker::VcfChunker;
use crate::cli::ArgParser;
use crate::regions::{build_region_tree, ContigMap};

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
    info!("starting kdp");
    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let mut input_vcf = vcf::reader::Builder::default()
        .build_from_path(args.io.input)
        .expect("Unable to parse vcf");
    let input_header = input_vcf.read_header().expect("Unable to parse header");
    let m_contigs = input_header.contigs().clone();

    let tree = build_region_tree(&m_contigs, args.io.bed);

    //let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();
    let mut m_input = VcfChunker::new(input_vcf, input_header, tree, args.kd.clone());
    let mut cnt = 0;
    for chunk in &mut m_input {
        cnt += 1;
        //println!("{}", cnt);
        //for i in chunk {
            //println!("{} {}", i.chromosome(), i.position());
        //}
    }
    //if let Ok(report) = guard.report().build() { println!("report: {:?}", &report); };
    println!("parsed {} entries", cnt);
    info!("finished kdp");
}
