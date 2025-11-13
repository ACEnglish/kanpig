use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use noodles_vcf::{self as vcf};
use std::{
    path::PathBuf,
    sync::{Arc, Mutex},
    thread::{self, JoinHandle},
};

use crate::{
    commands::KanpigCommand,
    file_validators,
    kplib::{
        build_region_tree, germ_genotyper::Genotyper, hp_sorter, open_reads, open_writer_thread,
        ChannelInput, ChannelOutput, GraphParams, PathScore, Ploidy, PloidyRegions, Variants,
        VcfChunker,
    },
};
fn task_thread(
    m_args: GermCommand,
    m_receiver: Receiver<ChannelInput>,
    m_result_sender: Sender<ChannelOutput>,
    m_ploidy: PloidyRegions,
) {
    let mut m_reads = open_reads(
        m_args.io.reads,
        m_args.io.reference,
        m_args.io.sample.expect("Sample should have been set"),
        0, // First sample is index 0 in the HaplotypeMeta vectros
        1, // One total sample will be opened (for HaplotypeMeta)
        &m_args.graph,
    );
    let genotyper = Genotyper::with_optional_config(m_args.gqconfig);
    loop {
        match m_receiver.recv() {
            Ok(None) | Err(_) => break,
            Ok(Some(chunk)) => {
                let mut m_graph = Variants::new(chunk, m_args.graph.kmer);
                debug!(
                    "Chunk {:?}:{:?}-{:?} w/ {}",
                    m_graph.chrom,
                    m_graph.start,
                    m_graph.end,
                    m_graph.node_indices.len() - 2
                );
                let ploidy = m_ploidy.get_ploidy(&m_graph.chrom, m_graph.start);
                // For zero, we don't have to waste time going into the bam
                if ploidy == Ploidy::Zero {
                    m_result_sender
                        .send(m_graph.take_annotated(vec![&[]], vec![0], vec![&ploidy], &genotyper))
                        .unwrap();
                    continue;
                }

                let (haps, coverage) =
                    m_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                let haps = ploidy.cluster(
                    haps,
                    coverage,
                    0,
                    m_args.hps_weight,
                    m_args.hapsim,
                    m_args.ab,
                    &m_args.graph,
                );

                // Only need to build the full graph sometimes
                let should_build = !haps.is_empty()
                    && !m_args.graph.one_to_one
                    && m_graph.node_indices.len() <= (m_args.graph.maxnodes + 2);
                m_graph.build(should_build);

                let mut paths: Vec<PathScore> = haps
                    .iter()
                    .map(|h| m_graph.apply_haplotype(h, &m_args.graph))
                    .filter(|p| *p != PathScore::default())
                    .collect();
                paths.sort_by(|a, b| hp_sorter(&a.meta.hp[0], &b.meta.hp[0]));
                m_result_sender
                    .send(m_graph.take_annotated(
                        vec![&paths],
                        vec![coverage],
                        vec![&ploidy],
                        &genotyper,
                    ))
                    .unwrap();
            }
        }
    }
    // This should give a result
}
#[derive(Parser, Debug, Clone)]
pub struct GermCommand {
    #[command(flatten)]
    pub io: IOParams,

    #[command(flatten)]
    pub graph: GraphParams,

    /// Clustering weight for haplotagged reads (off=0.0, full=1.0)
    #[arg(long, default_value_t = 1.0, help_heading = "Genotyping")]
    pub hps_weight: f32,

    /// Collapse haplotypes of similar size (off=1)
    #[arg(long, default_value_t = 1.0, help_heading = "Genotyping")]
    pub hapsim: f32,

    /// Minimum allele balance for compound het lower VAF (off=0)
    #[arg(long, default_value_t = 0.0, help_heading = "Genotyping")]
    pub ab: f32,

    /// GQ Config
    #[arg(long, help_heading = "Genotyping")]
    pub gqconfig: Option<PathBuf>,
}

#[derive(clap::Args, Clone, Debug)]
pub struct IOParams {
    /// VCF to genotype
    #[arg(short, long, help_heading = "I/O")]
    pub input: PathBuf,

    /// Reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(short, long, help_heading = "I/O")]
    pub reads: PathBuf,

    /// Reference genome
    #[arg(short = 'f', long, help_heading = "I/O")]
    pub reference: PathBuf,

    /// Output VCF (unsorted, uncompressed) [default: stdout]
    #[arg(short, long, help_heading = "I/O")]
    pub out: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 1, help_heading = "I/O")]
    pub threads: usize,

    /// Output VCF sample name
    #[arg(long, help_heading = "I/O")]
    pub sample: Option<String>,

    /// Bed file of non-diploid regions
    #[arg(long, help_heading = "I/O")]
    pub ploidy_bed: Option<PathBuf>,

    /// Regions to analyze
    #[arg(long, help_heading = "I/O")]
    pub bed: Option<PathBuf>,

    /// Verbose logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,
}

impl KanpigCommand for GermCommand {
    fn debug(&self) -> bool {
        self.io.debug
    }
    /// Validate command line arguments
    fn validate(&self) -> bool {
        let mut is_ok = true;

        is_ok &= file_validators::validate_file(&self.io.input, "--input");
        is_ok &= file_validators::validate_reads(&self.io.reads, &self.graph);
        is_ok &= file_validators::validate_reference(&self.io.reference);

        if let Some(bed_file) = &self.io.bed {
            is_ok &= file_validators::validate_file(bed_file, "--bed");
        }

        if self.graph.sizemin < 10 {
            warn!("--sizemin is recommended to be at least 10");
        }

        if self.graph.kmer >= 8 {
            warn!("--kmer above 8 becomes memory intensive");
        }

        if self.graph.kmer < 1 {
            error!("--kmer must be at least 1");
            is_ok = false;
        }

        if self.graph.sizemin < self.graph.kmer.into() {
            error!("--sizemin must be ≥ --kmer");
            is_ok = false;
        }

        if self.graph.sizesim < 0.0 || self.graph.sizesim > 1.0 {
            error!("--sizesim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.graph.seqsim < 0.0 || self.graph.seqsim > 1.0 {
            error!("--seqsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.hapsim < 0.0 || self.hapsim > 1.0 {
            error!("--hapsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.graph.maxpaths < 1 {
            error!("--maxpaths must be at least 1");
            is_ok = false;
        }

        if self.io.threads < 1 {
            error!("--threads must be at least 1");
            is_ok = false;
        }

        is_ok
    }

    fn run(&mut self) {
        let mut input_vcf = vcf::io::reader::Builder::default()
            .build_from_path(self.io.input.clone())
            .expect("Unable to parse vcf");

        let input_header = input_vcf.read_header().expect("Unable to parse vcf header");

        if self.io.sample.is_none() {
            if input_header.sample_names().is_empty() {
                error!("--input contains no samples. --sample name must be provided");
                std::process::exit(1);
            }
            let samp_name = input_header.sample_names()[0].clone();
            info!("Setting sample to {}", samp_name);
            self.io.sample = Some(samp_name);
        }

        let m_contigs = input_header.contigs().clone();

        let tree = build_region_tree(&m_contigs, &self.io.bed);

        let ploidy = PloidyRegions::new(&self.io.ploidy_bed);

        // Create channels for communication between threads
        let (task_sender, task_receiver): (Sender<ChannelInput>, Receiver<ChannelInput>) =
            unbounded();
        let (result_sender, result_receiver): (Sender<ChannelOutput>, Receiver<ChannelOutput>) =
            unbounded();

        info!("spawning {} threads", self.io.threads);
        let task_handles: Vec<JoinHandle<()>> = (0..self.io.threads)
            .map(|_| {
                let m_args = self.clone();
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

        let write_handler = open_writer_thread(
            result_receiver,
            self.io.out.clone(),
            vec![self
                .io
                .sample
                .clone()
                .expect("Sample should have already been set")],
            input_header.clone(),
            num_variants.clone(),
        );

        info!("building variant graphs");
        let mut m_input = VcfChunker::new(
            input_vcf,
            input_header.clone(),
            tree,
            &self.io.reference,
            self.graph.clone(),
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
        for _ in 0..self.io.threads {
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
}
