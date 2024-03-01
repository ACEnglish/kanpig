extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use noodles_vcf::{
    self as vcf,
    record::genotypes::{sample::Value, Genotypes},
};
//use rayon::prelude::*;
use crossbeam_channel::{unbounded, Receiver, Sender};
use std::fs::File;
use std::io::BufWriter;
use std::thread;

mod bamparser;
mod bedparser;
mod chunker;
mod cli;
mod haplotype;
mod kmeans;
mod kmer;
mod metrics;
mod pathscore;
mod pileup;
mod regions;
mod vargraph;
mod vcf_traits;

use crate::bamparser::BamParser;
use crate::chunker::VcfChunker;
use crate::cli::ArgParser;
use crate::pathscore::PathScore;
use crate::regions::build_region_tree;
use crate::vargraph::Variants;

// What are the channels going to be transferring?
type InputType = (ArgParser, Vec<vcf::Record>);
type OutputType = (Variants, PathScore, PathScore);

fn main() {
    pretty_env_logger::formatted_timed_builder()
        .filter_level(log::LevelFilter::Info)
        .init();

    let mut args = ArgParser::parse();
    info!("starting");
    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let mut input_vcf = vcf::reader::Builder::default()
        .build_from_path(args.io.input.clone())
        .expect("Unable to parse vcf");
    let input_header = input_vcf.read_header().expect("Unable to parse header");

    // Ensure sample is correctly setup
    if args.io.sample.is_none() {
        if input_header.sample_names().is_empty() {
            error!("--input contains no samples");
            std::process::exit(1);
        }
        args.io.sample = Some(input_header.sample_names()[0].clone());
        info!("Setting sample to {}", args.io.sample.as_ref().unwrap());
    } else if !input_header
        .sample_names()
        .contains(args.io.sample.as_ref().unwrap())
    {
        error!(
            "--sample {} not in --input {}",
            args.io.sample.unwrap(),
            args.io.input.display()
        );
        std::process::exit(1);
    }
    info!("params: {:?}", args);

    let m_contigs = input_header.contigs().clone();
    let tree = build_region_tree(&m_contigs, &args.io.bed);

    let out_buf = BufWriter::new(File::create(&args.io.out).expect("Error Creating Kmers"));
    let mut writer = vcf::Writer::new(out_buf);
    writer.write_header(&input_header);

    let mut m_input = VcfChunker::new(input_vcf, input_header.clone(), tree, args.kd.clone());

    // Create channels for communication between threads
    let (sender, receiver): (Sender<Option<InputType>>, Receiver<Option<InputType>>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    // Spawn worker threads
    for _ in 0..args.io.threads {
        let receiver = receiver.clone();
        let result_sender = result_sender.clone();
        thread::spawn(move || {
            for item in receiver.into_iter().flatten() {
                let result = process_item(item);
                result_sender.send(result).unwrap();
            }
        });
    }

    // Send items to worker threads
    let mut num_chunks = 0;
    for i in &mut m_input {
        sender.send(Some((args.clone(), i))).unwrap();
        num_chunks += 1;
    }

    // Signal worker threads to exit
    for _ in 0..args.io.threads {
        sender.send(None).unwrap();
    }

    // Collect results from worker threads
    for _ in 0..num_chunks {
        let (m_graph, p1, p2) = result_receiver.recv().unwrap();

        for var_idx in m_graph.node_indices {
            let mut cur_var = match &m_graph.graph.node_weight(var_idx).unwrap().entry {
                Some(var) => var.clone(),
                None => {
                    continue;
                }
            };
            //https://docs.rs/noodles-vcf/0.49.0/noodles_vcf/record/genotypes/struct.Genotypes.html
            let gt = match (p1.path.contains(&var_idx), p2.path.contains(&var_idx)) {
                (true, true) => "1|1",
                (true, false) => "1|0",
                (false, true) => "0|1",
                (false, false) => "0|0",
            };
            let keys = "GT:GQ".parse().unwrap();
            let genotypes = Genotypes::new(
                keys,
                vec![vec![Some(Value::from(gt)), Some(Value::from(13))]],
            );

            *cur_var.genotypes_mut() = genotypes.clone();
            let _ = writer.write_record(&input_header, &cur_var);
        }
    }

    info!("finished");
}

fn process_item(chunk_in: InputType) -> OutputType {
    let (args, chunk) = chunk_in;
    let mut m_bam = BamParser::new(args.io.bam, args.io.reference, args.kd.clone());
    let m_graph = Variants::new(chunk, args.kd.kmer);
    //println!("Analyzing {:?}", m_graph);
    let (h1, h2) = m_bam.find_haps(&m_graph.chrom, m_graph.start, m_graph.end);
    // If the haplotypes are identical, we only need to traverse once
    // And it holds the overall coverage
    let p1 = m_graph.apply_coverage(&h1, &args.kd);
    //println!("H1 Best Path {:?}", p1);
    let p2 = m_graph.apply_coverage(&h2, &args.kd);
    //println!("H2 Best Path {:?}", p2);
    (m_graph, p1, p2)
}
