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
        build_region_tree, open_reads, open_writer_thread, ChannelInput, ChannelOutput, KDParams,
        PathScore, Ploidy, PloidyRegions, Variants, VcfChunker,
    },
};

fn task_thread(
    m_args: TrioCommand,
    m_receiver: Receiver<ChannelInput>,
    m_result_sender: Sender<ChannelOutput>,
    m_ploidy: PloidyRegions,
) {
    let mut pro_reads = open_reads(
        m_args.io.proband.clone(),
        m_args.io.reference.clone(),
        m_args.io.sample.expect("Sample should have been set"),
        0, // First sample is index 0 in the HaplotypeMeta vectros
        3, // One total sample will be opened (for HaplotypeMeta)
        &m_args.kd,
    );
    let mut mat_reads = open_reads(
        m_args.io.mother.clone(),
        m_args.io.reference.clone(),
        "Mat".to_string(),
        1, // First sample is index 0 in the HaplotypeMeta vectros
        3, // One total sample will be opened (for HaplotypeMeta)
        &m_args.kd,
    );
    let mut pat_reads = open_reads(
        m_args.io.father.clone(),
        m_args.io.reference.clone(),
        "Pat".to_string(),
        2, // First sample is index 0 in the HaplotypeMeta vectros
        3, // One total sample will be opened (for HaplotypeMeta)
        &m_args.kd,
    );

    loop {
        match m_receiver.recv() {
            Ok(None) | Err(_) => break,
            Ok(Some(chunk)) => {
                let mut m_graph = Variants::new(chunk, m_args.kd.kmer, m_args.kd.maxhom);

                let ploidy = m_ploidy.get_ploidy(&m_graph.chrom, m_graph.start);
                // For zero, we don't have to waste time going into the bam
                if ploidy == Ploidy::Zero {
                    m_result_sender
                        .send(m_graph.take_annotated(vec![&[]], vec![0], vec![&ploidy]))
                        .unwrap();
                    continue;
                }

                let (pro_haps, pro_coverage) =
                    pro_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                let pro_haps = ploidy.cluster(pro_haps, pro_coverage, 0, &m_args.kd);

                let (mat_haps, mat_coverage) =
                    mat_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                let mat_haps = ploidy.cluster(mat_haps, mat_coverage, 1, &m_args.kd);

                let (pat_haps, pat_coverage) =
                    pat_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                let pat_haps = ploidy.cluster(pat_haps, pat_coverage, 2, &m_args.kd);

                // Only need to build the full graph sometimes
                //let should_build = !pro_haps.is_empty()
                //&& !m_args.kd.one_to_one
                //&& m_graph.node_indices.len() <= (m_args.kd.maxnodes + 2);
                m_graph.build(true);

                // Haplotype to paths
                let paths: Vec<PathScore> = pro_haps
                    .into_iter()
                    .chain(mat_haps)
                    .chain(pat_haps)
                    .map(|h| m_graph.apply_haplotype(&h, &m_args.kd))
                    .filter(|p| *p != PathScore::default())
                    .collect();

                // This is weird and should maybe be done by apply haplotype?
                // or maybe in the separation below
                //paths.sort_by(|a, b| hp_sorter(&a.meta.hp[0], &b.meta.hp[0]));

                // Separate paths back out to the samples
                let num_samples = 3;
                let mut separated_paths: Vec<Vec<PathScore>> = vec![Vec::new(); num_samples];
                for path in paths {
                    for (bit, s_paths) in separated_paths.iter_mut().enumerate().take(num_samples) {
                        if (path.meta.samples_flag & (1 << bit)) != 0 {
                            s_paths.push(path.clone());
                        }
                    }
                }

                let separated_paths: Vec<&[PathScore]> =
                    separated_paths.iter().map(|bin| bin.as_slice()).collect();

                m_result_sender
                    .send(m_graph.take_annotated(
                        separated_paths,
                        vec![pro_coverage, mat_coverage, pat_coverage],
                        vec![&ploidy, &ploidy, &ploidy],
                    ))
                    .unwrap();
            }
        }
    }
    // This should give a result
}
#[derive(Parser, Debug, Clone)]
pub struct TrioCommand {
    #[command(flatten)]
    pub io: IOParams,

    #[command(flatten)]
    pub kd: KDParams,
}

#[derive(clap::Args, Clone, Debug)]
pub struct IOParams {
    /// VCF to genotype
    #[arg(short, long, help_heading = "I/O")]
    pub input: PathBuf,

    /// Proband reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(long, help_heading = "I/O")]
    pub proband: PathBuf,

    /// Maternal reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(long, help_heading = "I/O")]
    pub mother: PathBuf,

    /// Paternal reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(long, help_heading = "I/O")]
    pub father: PathBuf,

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

impl KanpigCommand for TrioCommand {
    fn debug(&self) -> bool {
        self.io.debug
    }
    /// Validate command line arguments
    fn validate(&self) -> bool {
        let mut is_ok = true;

        is_ok &= file_validators::validate_file(&self.io.input, "--input");
        is_ok &= file_validators::validate_reads(&self.io.proband, &self.kd);
        is_ok &= file_validators::validate_reads(&self.io.mother, &self.kd);
        is_ok &= file_validators::validate_reads(&self.io.father, &self.kd);
        is_ok &= file_validators::validate_reference(&self.io.reference);

        if let Some(bed_file) = &self.io.bed {
            is_ok &= file_validators::validate_file(bed_file, "--bed");
        }

        if self.kd.sizemin < 10 {
            warn!("--sizemin is recommended to be at least 10");
        }

        if self.kd.kmer >= 8 {
            warn!("--kmer above 8 becomes memory intensive");
        }

        if self.kd.kmer < 1 {
            error!("--kmer must be at least 1");
            is_ok = false;
        }

        if self.kd.sizemin < self.kd.kmer.into() {
            error!("--sizemin must be ≥ --kmer");
            is_ok = false;
        }

        if self.kd.sizesim < 0.0 || self.kd.sizesim > 1.0 {
            error!("--sizesim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.kd.seqsim < 0.0 || self.kd.seqsim > 1.0 {
            error!("--seqsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.kd.hapsim < 0.0 || self.kd.hapsim > 1.0 {
            error!("--hapsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.kd.maxpaths < 1 {
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
            vec![
                self.io
                    .sample
                    .clone()
                    .expect("Sample should have already been set"),
                "Mat".to_string(),
                "Pat".to_string(),
            ],
            input_header.clone(),
            num_variants.clone(),
        );

        info!("building variant graphs");
        let mut m_input = VcfChunker::new(
            input_vcf,
            input_header.clone(),
            tree,
            self.kd.clone(),
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
