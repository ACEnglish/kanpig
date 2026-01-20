use clap::{Parser, ValueEnum};
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
        build_region_tree,
        cluster::collapse_haplotypes,
        find_subintervals,
        germ_genotyper::Genotyper,
        open_reads, open_writer_thread,
        polycluster::{self, ToPolyCluParams},
        readparsers::collect_read_data,
        trio_genotyper::trio_genotyper,
        ChannelInput, ChannelOutput, GraphParams, PathScore, Ploidy, PloidyRegions, VariantGraph,
        VcfChunker,
    },
};

fn task_thread(
    m_args: TrioCommand,
    m_receiver: Receiver<ChannelInput>,
    m_result_sender: Sender<ChannelOutput>,
    m_ploidy: Vec<PloidyRegions>,
) {
    let pro_reads = open_reads(
        m_args.io.proband.clone(),
        m_args.io.reference.clone(),
        m_args.io.proband_sample.clone(),
        0, // First sample is index 0 in the HaplotypeMeta vectros
        3, // Three total samples will be opened (for HaplotypeMeta)
        &m_args.graph,
    );

    let pat_reads = open_reads(
        m_args.io.father.clone(),
        m_args.io.reference.clone(),
        m_args.io.father_sample.clone(),
        1,
        3,
        &m_args.graph,
    );

    let mat_reads = open_reads(
        m_args.io.mother.clone(),
        m_args.io.reference.clone(),
        m_args.io.mother_sample.clone(),
        2,
        3,
        &m_args.graph,
    );

    // These need to be pulled out so we can do the polyclustering
    // on both TrioCommand and MosaicCommand
    let pclu_params = m_args.to_polyclu_params();
    let germ_genotyper = Genotyper::from_config_file(m_args.gqconfig);
    let mut reads = vec![pro_reads, pat_reads, mat_reads];

    loop {
        match m_receiver.recv() {
            Ok(None) | Err(_) => break,
            Ok(Some(chunk)) => {
                let mut m_graph = VariantGraph::new(chunk, m_args.graph.kmer);

                let ploidy_owned: Vec<Ploidy> = m_ploidy
                    .iter()
                    .map(|p| p.get_ploidy(&m_graph.chrom, m_graph.start))
                    .collect();
                let ploidy: Vec<&Ploidy> = ploidy_owned.iter().collect();

                let mut pileup_data = collect_read_data(&mut reads, &m_graph);

                let subintv = find_subintervals(&pileup_data.reads, m_args.graph.neighdist);
                for (sub_start, sub_end) in subintv.into_iter() {
                    let Some(mut subgraph) = m_graph.make_subgraph(sub_start, sub_end) else {
                        continue;
                    };
                    // This is the block we're reusing in the different parts
                    let (haplos, coverages, ref_coverage) =
                        pileup_data.subset_to_interval(sub_start, sub_end, m_args.graph.kmer);

                    let n_haps = haplos.len();
                    if n_haps <= m_args.graph.mincoverage || n_haps > m_args.graph.maxcoverage {
                        debug!(
                            "Region skipped extreme coverage {}x @ {}:{}-{}",
                            n_haps, subgraph.chrom, sub_start, sub_end
                        );
                        m_result_sender
                            .send(subgraph.take_annotated(
                                vec![&[], &[], &[]],
                                coverages.to_vec(),
                                &ploidy,
                                &germ_genotyper,
                            ))
                            .unwrap();
                        continue;
                    }

                    let cluster_result = polycluster::perform_clustering(&haplos, &pclu_params, 3);

                    let read_counts = polycluster::count_reads(
                        cluster_result.k,
                        &ref_coverage,
                        &cluster_result.assignments,
                        &haplos,
                    );

                    debug!("Read Counts:\n {:?}", read_counts);
                    let gts = trio_genotyper(&read_counts, &cluster_result.quality);

                    let clustered_haps = collapse_haplotypes(cluster_result, haplos, gts);

                    let should_build = !clustered_haps.is_empty()
                        && !m_args.graph.one_to_one
                        && subgraph.node_indices.len() <= (m_args.graph.maxnodes + 2);
                    subgraph.build(should_build);

                    let paths: Vec<PathScore> = clustered_haps
                        .into_iter()
                        .map(|h| subgraph.apply_haplotype(&h, &m_args.graph))
                        .filter(|p| *p != PathScore::default())
                        .collect();

                    let separated_paths = polycluster::separate_paths_by_sample(paths, 3);

                    m_result_sender
                        .send(subgraph.take_annotated(
                            separated_paths.iter().map(|bin| bin.as_slice()).collect(),
                            coverages.to_vec(),
                            &ploidy,
                            &germ_genotyper,
                        ))
                        .unwrap();
                }

                // Remaining refcovered
                let remaining_variants =
                    m_graph.take_refcovered(pileup_data.coverages, &ploidy, &germ_genotyper);
                // Guard because I don't know what happens to variant-less graphs
                if remaining_variants.is_some() {
                    m_result_sender.send(remaining_variants).unwrap();
                }
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
    pub graph: GraphParams,

    /// Minimum number of reads in a cluster
    #[arg(long, default_value_t = 5, help_heading = "Genotyping")]
    pub msmin: usize,

    /// Max clusters
    #[arg(long, default_value_t = 5, help_heading = "Genotyping")]
    pub maxclust: usize,

    /// Clustering weight for haplotagged reads (off=0.0, full=1.0)
    #[arg(long, default_value_t = 0.25, help_heading = "Genotyping")]
    pub hps_weight: f32,

    /// Clustering weight for haplotype lengths (off=0.0, full=1.0)
    #[arg(long, default_value_t = 0.25, help_heading = "Genotyping")]
    pub len_weight: f32,

    /// Only cluster on haplotype lengths
    #[arg(
        long,
        default_value_t = false,
        help_heading = "Genotyping",
        hide = true
    )]
    pub lengthonly: bool,

    /// GQ Config
    #[arg(long, help_heading = "Genotyping")]
    pub gqconfig: Option<PathBuf>,
}

impl ToPolyCluParams for TrioCommand {
    fn to_polyclu_params(&self) -> polycluster::PolyCluParams {
        polycluster::PolyCluParams {
            msmin: self.msmin,
            maxclust: self.maxclust,
            hps_weight: self.hps_weight,
            len_weight: self.len_weight,
            lengthonly: self.lengthonly,
            minkfreq: self.graph.minkfreq,
            ..Default::default() // Fill remaining fields with defaults
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
#[clap(rename_all = "UPPER")]
pub enum Karyotype {
    XX,
    XY,
}

#[derive(clap::Args, Clone, Debug)]
pub struct IOParams {
    /// VCF to genotype
    #[arg(short, long, help_heading = "I/O")]
    pub input: PathBuf,

    /// Proband reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(long, help_heading = "I/O")]
    pub proband: PathBuf,

    /// Paternal reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(long, help_heading = "I/O")]
    pub father: PathBuf,

    /// Maternal reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(long, help_heading = "I/O")]
    pub mother: PathBuf,

    /// Reference genome
    #[arg(short = 'f', long, help_heading = "I/O")]
    pub reference: PathBuf,

    /// Output VCF (unsorted, uncompressed) [default: stdout]
    #[arg(short, long, help_heading = "I/O")]
    pub out: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 1, help_heading = "I/O")]
    pub threads: usize,

    /// Output RNAMES file
    #[arg(long, help_heading = "I/O")]
    pub rnames: Option<PathBuf>,

    /// Output VCF proband sample name
    #[arg(long, default_value = "PRO", help_heading = "I/O")]
    pub proband_sample: String,

    /// Output VCF paternal sample name
    #[arg(long, default_value = "PAT", help_heading = "I/O")]
    pub father_sample: String,

    /// Output VCF maternal sample name
    #[arg(long, default_value = "MAT", help_heading = "I/O")]
    pub mother_sample: String,

    /// Bed file of XY karyotype
    #[arg(long = "XYploidy-bed", help_heading = "I/O")]
    pub xyploidy_bed: Option<PathBuf>,

    /// Bed file of XX karyotype
    #[arg(long = "XXploidy-bed", help_heading = "I/O")]
    pub xxploidy_bed: Option<PathBuf>,

    /// Proband karyotype
    #[arg(long, help_heading = "I/O")]
    pub karyotype: Option<Karyotype>,

    /// Regions to analyze
    #[arg(long, help_heading = "I/O")]
    pub bed: Option<PathBuf>,

    /// Minimum coverage of a region to analyze
    #[arg(long, default_value_t = 1, help_heading = "I/O")]
    pub min_coverage: usize,

    /// Maximum coverage of a region to analyze
    #[arg(long, default_value_t = 1000, help_heading = "I/O")]
    pub max_coverage: usize,

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
        is_ok &= file_validators::validate_reads(&self.io.proband, &self.graph);
        is_ok &= file_validators::validate_reads(&self.io.mother, &self.graph);
        is_ok &= file_validators::validate_reads(&self.io.father, &self.graph);
        is_ok &= file_validators::validate_reference(&self.io.reference);

        if let Some(bed_file) = &self.io.bed {
            is_ok &= file_validators::validate_file(bed_file, "--bed");
        }

        is_ok &= match (
            &self.io.xxploidy_bed,
            &self.io.xyploidy_bed,
            &self.io.karyotype,
        ) {
            (Some(b1), Some(b2), Some(_k)) => {
                file_validators::validate_file(b1, "--XXploidy-bed")
                    && file_validators::validate_file(b2, "--XYploidy-bed")
            }
            (None, None, None) => true,
            _ => {
                warn!(
                    "--karyotype, --XXploidy-bed, and --XYploidy-bed must all or none be specified"
                );
                false
            }
        };

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

        if self.graph.maxpaths < 1 {
            error!("--maxpaths must be at least 1");
            is_ok = false;
        }

        if self.io.threads < 1 {
            error!("--threads must be at least 1");
            is_ok = false;
        }

        if self.maxclust < 4 {
            error!("--maxclust must be at least 4");
            is_ok = false;
        }

        is_ok
    }

    fn run(&mut self) {
        let mut input_vcf = vcf::io::reader::Builder::default()
            .build_from_path(self.io.input.clone())
            .expect("Unable to parse vcf");

        let input_header = input_vcf.read_header().expect("Unable to parse vcf header");

        info!(
            "Setting samples to {}, {}, {}",
            self.io.proband_sample, self.io.father_sample, self.io.mother_sample
        );

        let m_contigs = input_header.contigs().clone();

        let tree = build_region_tree(&m_contigs, &self.io.bed);

        let xy_ploidy = PloidyRegions::new(&self.io.xyploidy_bed);
        let xx_ploidy = PloidyRegions::new(&self.io.xxploidy_bed);
        let pro_ploidy = match self.io.karyotype {
            Some(Karyotype::XY) => xy_ploidy.clone(),
            Some(Karyotype::XX) => xx_ploidy.clone(),
            _ => PloidyRegions::new(&None),
        };

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
                let m_pro_ploidy = pro_ploidy.clone();
                let m_xy_ploidy = xy_ploidy.clone();
                let m_xx_ploidy = xx_ploidy.clone();

                thread::spawn(move || {
                    task_thread(
                        m_args,
                        m_receiver,
                        m_result_sender,
                        vec![m_pro_ploidy, m_xy_ploidy, m_xx_ploidy],
                    );
                })
            })
            .collect();

        // Before we start the workers, we'll start the writer
        // This is the semaphore for the progress bar that communicates between main and writer
        let num_variants = Arc::new(Mutex::new(0));

        let write_handler = open_writer_thread(
            result_receiver,
            self.io.out.clone(),
            self.io.rnames.clone(),
            vec![
                self.io.proband_sample.clone(),
                self.io.father_sample.clone(),
                self.io.mother_sample.clone(),
            ],
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
