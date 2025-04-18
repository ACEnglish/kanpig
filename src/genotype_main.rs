use crossbeam_channel::{unbounded, Receiver, Sender};
use indicatif::{ProgressBar, ProgressStyle};
use noodles_vcf::{self as vcf};
use rust_htslib::faidx;
use std::{
    sync::{Arc, Mutex},
    thread::{self, JoinHandle},
};

use crate::kplib::{
    build_region_tree, BamParser, GTArgs, GenotypeAnno, IOParams, PathScore, Ploidy, PloidyRegions,
    PlupParser, ReadParser, Variants, VcfChunker, VcfWriter,
};

type InputType = Option<Vec<vcf::variant::RecordBuf>>;
type OutputType = Option<Vec<GenotypeAnno>>;

fn hp_sorter(a: &Option<u8>, b: &Option<u8>) -> std::cmp::Ordering {
    match (a, b) {
        // If both are Some, reverse order
        (Some(va), Some(vb)) => vb.cmp(va),

        // If one is None and the other is Some(1), None comes first
        (Some(1), None) => std::cmp::Ordering::Greater,
        (None, Some(1)) => std::cmp::Ordering::Less,

        // If one is None and the other is Some(2), None comes last
        (Some(2), None) => std::cmp::Ordering::Less,
        (None, Some(2)) => std::cmp::Ordering::Greater,

        // If both are None, do nothing
        (None, None) => std::cmp::Ordering::Equal,

        // Default fallback (not strictly needed with the cases above)
        _ => std::cmp::Ordering::Equal,
    }
}

fn write_thread(
    result_receiver: Receiver<OutputType>,
    wt_io: IOParams,
    wt_header: vcf::Header,
    wt_num_variants: Arc<Mutex<u64>>,
) {
    let mut m_writer = VcfWriter::new(&wt_io.out, wt_header.clone(), &wt_io.sample);

    let mut pbar: Option<ProgressBar> = None;
    let sty =
        ProgressStyle::with_template(" [{elapsed_precise}] {bar:44.cyan/blue} > {pos} completed")
            .unwrap()
            .progress_chars("・🐷🥫");

    let mut completed_variants: u64 = 0;
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
                    m_writer.anno_write(entry);
                    rsize += 1;
                }

                if let Some(ref mut bar) = pbar {
                    bar.inc(rsize);
                } else {
                    completed_variants += rsize;
                    // check if the reader is finished so we can setup the pbar
                    let value = *wt_num_variants.lock().unwrap();
                    if value != 0 {
                        let t_bar = ProgressBar::new(value).with_style(sty.clone());
                        t_bar.inc(completed_variants);
                        pbar = Some(t_bar);
                    }
                }
            }
        }
    }
    if m_writer.iupac_fixed {
        warn!("Some IUPAC codes in REF sequences have been fixed in output");
    }
    info!("genotype counts: {:#?}", m_writer.gtcounts);
}

fn task_thread(
    m_args: GTArgs,
    m_receiver: Receiver<InputType>,
    m_result_sender: Sender<OutputType>,
    m_ploidy: PloidyRegions,
) {
    let reference = faidx::Reader::from_path(&m_args.io.reference).unwrap();
    let mut m_reads: Box<dyn ReadParser> =
        match m_args.io.reads.file_name().and_then(|name| name.to_str()) {
            Some(name) if name.ends_with(".plup.gz") => Box::new(PlupParser::new(
                m_args.io.reads,
                reference,
                m_args.kd.clone(),
            )),
            _ => Box::new(BamParser::new(
                m_args.io.reads,
                m_args.io.reference,
                reference,
                m_args.kd.clone(),
            )),
        };

    loop {
        match m_receiver.recv() {
            Ok(None) | Err(_) => break,
            Ok(Some(chunk)) => {
                let mut m_graph = Variants::new(chunk, m_args.kd.kmer, m_args.kd.maxhom);

                let ploidy = m_ploidy.get_ploidy(&m_graph.chrom, m_graph.start);
                // For zero, we don't have to waste time going into the bam
                if ploidy == Ploidy::Zero {
                    m_result_sender
                        .send(Some(m_graph.take_annotated(&[], 0, &ploidy)))
                        .unwrap();
                    continue;
                }

                let (haps, coverage) =
                    m_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                let haps = ploidy.cluster(haps, coverage, &m_args.kd);

                // Only need to build the full graph sometimes
                let should_build = !haps.is_empty()
                    && !m_args.kd.one_to_one
                    && m_graph.node_indices.len() <= (m_args.kd.maxnodes + 2);
                m_graph.build(should_build);

                let mut paths: Vec<PathScore> = haps
                    .iter()
                    .map(|h| m_graph.apply_coverage(h, &m_args.kd))
                    .filter(|p| *p != PathScore::default())
                    .collect();
                paths.sort_by(|a, b| hp_sorter(&a.hp, &b.hp));

                // Sort paths based on their HP if set
                m_result_sender
                    .send(Some(m_graph.take_annotated(&paths, coverage, &ploidy)))
                    .unwrap();
            }
        }
    }
    // This should give a result
}

pub fn genotype_main(args: GTArgs) {
    let mut input_vcf = vcf::io::reader::Builder::default()
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
                task_thread(m_args, m_receiver, m_result_sender, m_ploidy);
            })
        })
        .collect();

    // Before we start the workers, we'll start the writer
    // This is the semaphore for the progress bar that communicates between main and writer
    let num_variants = Arc::new(Mutex::new(0));

    let wt_io = args.io.clone();
    let wt_header = input_header.clone();
    let wt_num_variants = num_variants.clone();

    let write_handler = thread::spawn(move || {
        write_thread(result_receiver, wt_io, wt_header.clone(), wt_num_variants);
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
    for i in &mut m_input {
        task_sender.send(Some(i)).unwrap();
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
        info!("genotyping {} variants", value_guard);
    }

    for handle in task_handles {
        handle.join().unwrap();
    }

    // There will be no more results made
    result_sender.send(None).unwrap();

    // Wait on the writer
    write_handler.join().unwrap();
    info!("finished");
}
