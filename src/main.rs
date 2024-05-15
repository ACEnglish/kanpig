extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use crossbeam_channel::{select, unbounded, Receiver, Sender};
use indicatif::{ProgressBar, ProgressStyle};
use noodles_vcf::{self as vcf};
use std::thread;
mod kplib;

use kplib::{
    build_region_tree, diploid_haplotypes, haploid_haplotypes, ArgParser, BamParser, GenotypeAnno,
    PathScore, Ploidy, PloidyRegions, Variants, VcfChunker, VcfWriter,
};

type InputType = Vec<vcf::Record>;
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

    // So this object... will have a bed file of chrom\tstart\tend\tploidy
    // if its 0, we will put everything as ./.
    // if its 1, we'll make 1 or 0
    // default is diploid e.g. 0/1
    // And we'll only allow 0,1
    // So this means there will be a make and female bed.
    // female bed will set chrY to 0
    // male bed will chrY 1 and non-par chrX 1
    // After you build this thing, test it by seeing if you can get Zero to mask out variants
    let ploidy = PloidyRegions::new(&args.io.ploidy_bed);

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
        let m_args = args.clone();
        let receiver = receiver.clone();
        let result_sender = result_sender.clone();
        let m_ploidy = ploidy.clone();
        thread::spawn(move || {
            let mut m_bam = BamParser::new(m_args.io.bam, m_args.io.reference, m_args.kd.clone());
            for chunk in receiver.into_iter().flatten() {
                let mut m_graph = Variants::new(chunk, m_args.kd.kmer, m_args.kd.maxhom);

                let ploidy = m_ploidy.get_ploidy(&m_graph.chrom, m_graph.start);
                // For zero, we don't have to waste time going into the bam
                if ploidy == Ploidy::Zero {
                    result_sender
                        .send(m_graph.take_annotated(&[], 0, &ploidy))
                        .unwrap();
                    continue;
                }

                let (haps, coverage) = m_bam.find_haps(&m_graph.chrom, m_graph.start, m_graph.end);

                let haps = match ploidy {
                    Ploidy::Haploid => haploid_haplotypes(haps, coverage, &m_args.kd),
                    _ => diploid_haplotypes(haps, coverage, &m_args.kd),
                    // and then eventually this could allow a --ploidy flag to branch to
                    // polyploid_haplotypes
                };

                let paths: Vec<PathScore> = haps
                    .iter()
                    .map(|h| m_graph.apply_coverage(h, &m_args.kd))
                    .collect();

                result_sender
                    .send(m_graph.take_annotated(&paths, coverage, &ploidy))
                    .unwrap();
            }
        });
    }

    // Send items to worker threads
    let mut num_chunks: u64 = 0;
    info!("parsing input");
    for i in &mut m_input {
        sender.send(Some(i)).unwrap();
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
        ProgressStyle::with_template(" [{elapsed_precise}] {bar:44.cyan/blue} > {pos} completed")
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
    info!("genotype counts: {:#?}", writer.gtcounts);
    info!("finished");
}
