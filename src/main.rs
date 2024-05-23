extern crate pretty_env_logger;

#[macro_use]
extern crate log;

use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use indicatif::{ProgressBar, ProgressStyle};
use noodles_vcf::{self as vcf};
use peak_alloc::PeakAlloc;
use std::sync::{Arc, Mutex};
use std::thread;
use std::thread::JoinHandle;

mod kplib;

use kplib::{
    build_region_tree, diploid_haplotypes, haploid_haplotypes, ArgParser, BamParser, GenotypeAnno,
    PathScore, Ploidy, PloidyRegions, Variants, VcfChunker, VcfWriter,
};

type InputType = Option<Vec<vcf::Record>>;
type OutputType = Option<Vec<GenotypeAnno>>;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

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
    info!("params: {:#?}", args);
    if !args.validate() {
        error!("please fix arguments");
        std::process::exit(1);
    }

    let mut input_vcf = vcf::reader::Builder::default()
        .build_from_path(args.io.input.clone())
        .expect("Unable to parse vcf");
    let input_header = input_vcf.read_header().expect("Unable to parse vcf header");

    let m_contigs = input_header.contigs().clone();
    let tree = build_region_tree(&m_contigs, &args.io.bed);

    let ploidy = PloidyRegions::new(&args.io.ploidy_bed);

    // Create channels for communication between threads
    let (task_sender, task_receiver): (Sender<InputType>, Receiver<InputType>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    info!("spawning {} threads", args.io.threads);
    let task_handles: Vec<JoinHandle<()>> = (0..args.io.threads)
        .map(|_| {
            let m_args = args.clone();
            let m_receiver = task_receiver.clone();
            let m_result_sender = result_sender.clone();
            let m_ploidy = ploidy.clone();

            thread::spawn(move || {
                let mut m_bam =
                    BamParser::new(m_args.io.bam, m_args.io.reference, m_args.kd.clone());
                loop {
                    match m_receiver.recv() {
                        Ok(None) | Err(_) => break,
                        Ok(Some(chunk)) => {
                            let mut m_graph =
                                Variants::new(chunk, m_args.kd.kmer, m_args.kd.maxhom);

                            let ploidy = m_ploidy.get_ploidy(&m_graph.chrom, m_graph.start);
                            // For zero, we don't have to waste time going into the bam
                            if ploidy == Ploidy::Zero {
                                m_result_sender
                                    .send(Some(m_graph.take_annotated(&[], 0, &ploidy)))
                                    .unwrap();
                                continue;
                            }

                            let (haps, coverage) =
                                m_bam.find_haps(&m_graph.chrom, m_graph.start, m_graph.end);

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

                            m_result_sender
                                .send(Some(m_graph.take_annotated(&paths, coverage, &ploidy)))
                                .unwrap();
                        }
                    }
                }
                // This should give a result
            })
        })
        .collect();

    //Before we start the workers, we'll start the writer
    let writer = Arc::new(Mutex::new(VcfWriter::new(
        &args.io.out,
        input_header.clone(),
        &args.io.sample,
    )));

    let cloned_writer = writer.clone();
    let num_variants = Arc::new(Mutex::new(0));
    let cloned_num_variants = num_variants.clone();

    let write_handler = std::thread::spawn(move || {
        let sty = ProgressStyle::with_template(
            " [{elapsed_precise}] {bar:44.cyan/blue} > {pos} completed",
        )
        .unwrap()
        .progress_chars("##-");
        let mut pbar: Option<ProgressBar> = None;
        let mut phase_group: i32 = 0;
        let mut completed_variants: u64 = 0;
        let mut m_writer = cloned_writer.lock().unwrap();

        loop {
            match result_receiver.recv() {
                Ok(None) | Err(_) => {
                    pbar.expect("I actually shouldn't be expecting the bar")
                        .finish();
                    break;
                }
                Ok(Some(result)) => {
                    let mut rsize: u64 = 0;
                    for entry in result {
                        m_writer.anno_write(entry, phase_group);
                        rsize += 1;
                    }

                    if let Some(ref mut bar) = pbar {
                        bar.inc(rsize);
                    } else {
                        completed_variants += rsize;
                        // check if the reader is finished so we can setup the pbar
                        let value_guard = cloned_num_variants.lock().unwrap();
                        if *value_guard != 0 {
                            let t_bar = ProgressBar::new(*value_guard).with_style(sty.clone());
                            t_bar.inc(completed_variants);
                            pbar = Some(t_bar);
                        }
                    }
                    phase_group += 1;
                }
            }
        }
    });

    info!("building variant graphs");
    let mut m_input = VcfChunker::new(
        input_vcf,
        input_header.clone(),
        tree,
        args.kd.clone(),
        result_sender.clone(),
    );

    // Send items to worker threads
    let mut hang_tracker: i64 = 0;
    for i in &mut m_input {
        task_sender.send(Some(i)).unwrap();

        // The reader can get way ahead of the tasks, so we monitor memory usage
        // and let threads catch up
        let mut num_waits = 0;
        while PEAK_ALLOC.current_usage_as_gb() >= args.io.mem && num_waits < 10 {
            std::thread::sleep(std::time::Duration::from_millis(100));
            warn!(
                "throttling vcf reading with memory @ {}",
                PEAK_ALLOC.current_usage_as_gb()
            );
            num_waits += 1;
        }
        if num_waits != 0 {
            hang_tracker += 1;
        } else {
            hang_tracker -= 1;
        }
        if hang_tracker >= 10 {
            warn!("memory seems pretty full, consider setting a higher --mem");
            hang_tracker = 0;
        }
        hang_tracker = hang_tracker.max(-10);
    }

    if m_input.chunk_count == 0 {
        error!("No variants to be analyzed");
        std::process::exit(1);
    }

    // Signal worker threads to exit
    for _ in 0..args.io.threads {
        task_sender.send(None).unwrap();
    }

    // We now know how many variants will be parsed and can turn on the bar
    {
        let mut value_guard = num_variants.lock().unwrap();
        *value_guard = m_input.call_count + m_input.skip_count;
    }

    info!("genotyping");
    for handle in task_handles {
        handle.join().unwrap();
    }

    // There will be no more results made
    result_sender.send(None).unwrap();

    // Wait on the writer
    write_handler.join().unwrap();
    let writer = writer.lock().unwrap();
    info!("genotype counts: {:#?}", writer.gtcounts);
    info!("finished");
}
