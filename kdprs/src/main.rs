extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use noodles_vcf::{
    self as vcf,
    record::genotypes::{sample::Value, Genotypes},
};
use std::fs::File;
use std::io::BufWriter;

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
use crate::regions::build_region_tree;
use crate::vargraph::Variants;

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
    let tree = build_region_tree(&m_contigs, args.io.bed);

    //let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();

    let out_buf = BufWriter::new(File::create(&args.io.out).expect("Error Creating Kmers"));
    let mut writer = vcf::Writer::new(out_buf);
    writer.write_header(&input_header);

    let mut m_input = VcfChunker::new(input_vcf, input_header.clone(), tree, args.kd.clone());
    let mut m_bam = BamParser::new(args.io.bam, args.io.reference, args.kd.clone());
    /*
     * threading strategy: m_input should be windowed and each is sent to a thread
     * Because m_bam.open(), each thread will then have its own file handler to the bam/reference
     * The information I'm missing is how many variants go to each thread. I could send them
     * independentally, but that'll make a lot of opens. And since I don't know.. wait, I do know
     * how many parts there are because of the region tree. So the tree total size / threads is how
     * approximately how many parts I'll have. I might end up with some idle threads. But lets not
     * over optimize, yet. iter.chunks() works, I guess. But I'll have to 'magic number' the size
     */
    //
    for chunk in &mut m_input {
        let m_graph = Variants::new(chunk, args.kd.kmer);
        println!("Analyzing {:?}", m_graph);
        let (h1, h2) = m_bam.find_haps(&m_graph.chrom, m_graph.start, m_graph.end);
        let p1 = m_graph.apply_coverage(&h1, &args.kd);
        let p2 = m_graph.apply_coverage(&h2, &args.kd);
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
            // Should be handling
            let _ = writer.write_record(&input_header, &cur_var);
            /*if var_idx == 0 || var_idx == m_graph.graph.node_count() - 1 {
                continue
            }
            gt = var in p1.path + | + var in p2.path
            m_graph.nodes[var_idx].variant.unwrap()
            if var in p1
            */
        }
        //println!("Graph: {:?}", m_graph);

        // nodes that have None for both coverage are 0/0
        // nodes that have None for half shouldn't happen...? They will, actually.
        // So I need to change the sim to just a Vec and I'll push the coverages on
        // If it has one, then ... well I guess I don't know
        // Maybe I shouldn't be applying coverage but returning the best path
        // Then I can m_graph.pull_variants(path1, path2, coverage)
    }

    //if let Ok(report) = guard.report().build() { println!("report: {:?}", &report); };

    info!("finished");
}
