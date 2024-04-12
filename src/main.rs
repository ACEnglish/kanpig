extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use crossbeam_channel::{select, unbounded, Receiver, Sender};
use indicatif::{ProgressBar, ProgressStyle};
use noodles_vcf::{self as vcf};
use std::thread;
mod kanpig;

use kanpig::{
    build_region_tree, cluster_haplotypes, ArgParser, BamParser, GenotypeAnno, Variants,
    VcfChunker, VcfWriter,
};

type InputType = (ArgParser, Vec<vcf::Record>);
type OutputType = Vec<GenotypeAnno>;

fn main() {
    let args = ArgParser::parse();
    let level = if args.io.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    pretty_env_logger::formatted_timed_builder()
        .filter_level(level)
        .init();

    info!("starting");
    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }
    info!("params: {:#?}", args);

    let mut input_vcf = vcf::reader::Builder::default()
        .build_from_path(args.io.input.clone())
        .expect("Unable to parse vcf");
    let input_header = input_vcf.read_header().expect("Unable to parse vcf header");

    let m_contigs = input_header.contigs().clone();
    let tree = build_region_tree(&m_contigs, &args.io.bed);

    let mut writer = VcfWriter::new(&args.io.out, input_header.clone(), &args.io.sample);

    // We send the writer to the reader so that we can pipe filtered variants forward
    let mut m_input = VcfChunker::new(
        input_vcf,
        input_header.clone(),
        tree,
        args.kd.clone(),
        &mut writer,
    );

    // Create channels for communication between threads
    let (sender, receiver): (Sender<Option<InputType>>, Receiver<Option<InputType>>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    info!("spawning {} threads", args.io.threads);
    for _ in 0..args.io.threads {
        let receiver = receiver.clone();
        let result_sender = result_sender.clone();
        thread::spawn(move || {
            for (m_args, chunk) in receiver.into_iter().flatten() {
                let mut m_graph = Variants::new(chunk, m_args.kd.kmer);
                let mut m_bam =
                    BamParser::new(m_args.io.bam, m_args.io.reference, m_args.kd.clone());
                let (haps, coverage) = m_bam.find_haps(&m_graph.chrom, m_graph.start, m_graph.end);
                let (h1, h2) = cluster_haplotypes(haps, coverage, &m_args.kd);
                debug!("Hap1 out {:?}", h1);
                debug!("Hap2 out {:?}", h2);
                let p1 = m_graph.apply_coverage(&h1, &m_args.kd);
                let p2 = m_graph.apply_coverage(&h2, &m_args.kd);
                result_sender
                    .send(m_graph.to_annotated(&p1, &p2, coverage))
                    .unwrap();
            }
        });
    }

    // Send items to worker threads
    let mut num_chunks: u64 = 0;
    info!("parsing input");
    for i in &mut m_input {
        sender.send(Some((args.clone(), i))).unwrap();
        num_chunks += 1;
    }

    if num_chunks == 0 {
        error!("No variants to be analyzed");
        std::process::exit(1);
    }

    // Signal worker threads to exit
    for _ in 0..args.io.threads {
        sender.send(None).unwrap();
    }

    info!("collecting output");
    let sty =
        ProgressStyle::with_template(" [{elapsed_precise}] {bar:45.cyan/blue} > {pos:>7}/{len:7}")
            .unwrap()
            .progress_chars("##-");
    let pbar = ProgressBar::new(num_chunks);
    pbar.set_style(sty.clone());

    let mut phase_group: i32 = 0;
    loop {
        select! {
            recv(result_receiver) -> result => {
                match result {
                    Ok(annotated_entries) => {
                        for entry in annotated_entries {
                            writer.anno_write(entry, phase_group);
                        }
                        phase_group += 1;
                        pbar.inc(1);
                        if phase_group as u64 == num_chunks {
                            break;
                        }
                    },
                    Err(e) => {
                        debug!("Problem {:?}", e);
                        break;
                    }
                }
            },
        }
    }
    pbar.finish();

    info!("finished");
}
