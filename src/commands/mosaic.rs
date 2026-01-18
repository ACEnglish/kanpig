use clap::{ArgAction, Parser};
use crossbeam_channel::{unbounded, Receiver, Sender};
use ndarray::Axis;
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
        annotator::FiltFlags,
        build_region_tree,
        cluster::collapse_haplotypes,
        germ_genotyper::{GTstate, GenotypeMode, Genotyper, GenotyperConfig},
        mosaic_genotyper::{GenotypeHypothesis, MosaicGenotyper},
        open_reads, open_writer_thread,
        readparsers::collect_read_data,
        polycluster::{self, ToPolyCluParams},
        ChannelInput, ChannelOutput, GraphParams, PathScore, Ploidy, PloidyRegions, ReadParser,
        VariantGraph, VcfChunker,
    },
};

fn separate_paths_by_vaf(
    paths: Vec<Vec<PathScore>>,
    gts: GenotypeHypothesis,
) -> (Vec<Vec<PathScore>>, Vec<Vec<PathScore>>) {
    let mut germ = vec![Vec::new(); paths.len()];
    let mut soma = vec![Vec::new(); paths.len()];

    for (idx, sample_paths) in paths.into_iter().enumerate() {
        for path in sample_paths {
            if gts.germline_alleles.contains(&path.meta.id) {
                germ[idx].push(path);
            } else {
                soma[idx].push(path);
            }
        }
    }

    (germ, soma)
}

fn task_thread(
    m_args: MosaicCommand,
    m_receiver: Receiver<ChannelInput>,
    m_result_sender: Sender<ChannelOutput>,
    m_ploidy: PloidyRegions,
) {
    // For each bam (mainly one) open_reads add to list
    let n_samples = m_args.io.reads.len();
    let mut reads: Vec<Box<dyn ReadParser>> = m_args
        .io
        .reads
        .iter()
        .zip(m_args.io.sample.iter())
        .enumerate()
        .map(|(idx, (filename, sample))| {
            open_reads(
                filename.clone(),
                m_args.io.reference.clone(),
                sample.clone(),
                idx,
                n_samples,
                &m_args.graph,
            )
        })
        .collect();

    let pclu_params = m_args.to_polyclu_params();
    let genotyper = MosaicGenotyper::with_params(
        0.001, // This works fine enough, not worth it to tweak, probably
        m_args.alpha,
        m_args.beta,
        m_args.soma_vaf,
        1,
    );
    let cfg = GenotyperConfig {
        mode: GenotypeMode::Bino,
        ..Default::default()
    };
    let germ_genotyper = Genotyper::from_config(cfg);

    loop {
        match m_receiver.recv() {
            Ok(None) | Err(_) => break,
            Ok(Some(chunk)) => {
                let mut m_graph = VariantGraph::new(chunk, m_args.graph.kmer);

                let ploidy = m_ploidy.get_ploidy(&m_graph.chrom, m_graph.start);
                // For zero, we don't have to waste time going into the bam
                if ploidy == Ploidy::Zero {
                    m_result_sender
                        .send(m_graph.take_annotated(
                            vec![&[]; n_samples],
                            vec![0; n_samples],
                            vec![&ploidy; n_samples],
                            &germ_genotyper,
                        ))
                        .unwrap();
                    continue;
                }

                let pileup_data = collect_read_data(&mut reads, &m_graph);

                let n_haps = pileup_data.haplos.len();
                if n_haps <= m_args.graph.mincoverage || n_haps > m_args.graph.maxcoverage {
                    debug!(
                        "Region skipped extreme coverage {}x @ {}:{}-{}",
                        n_haps, m_graph.chrom, m_graph.start, m_graph.end
                    );
                    m_result_sender
                        .send(m_graph.take_annotated(
                            vec![&[]; n_samples],
                            pileup_data.coverages.to_vec(),
                            vec![&ploidy; n_samples],
                            &germ_genotyper,
                        ))
                        .unwrap();
                    continue;
                }

                let cluster_result =
                    polycluster::perform_clustering(&pileup_data.haplos, &pclu_params, n_samples);

                let read_counts = polycluster::count_reads(
                    cluster_result.k,
                    &pileup_data.ref_coverage,
                    &cluster_result.assignments,
                    &pileup_data.haplos,
                );

                debug!("Read Counts:\n {:?}", read_counts);
                let allele_support: Vec<u32> = read_counts
                    .axis_iter(Axis(0)) // Iterate over cols
                    .map(|row| row.sum() as u32)
                    .collect();
                let gts = genotyper.genotype(&allele_support);

                debug!("GTs; {:#?}", gts);
                let gts = match gts {
                    Some(g) => g,
                    None => {
                        m_result_sender
                            .send(m_graph.take_annotated(
                                vec![&[]; n_samples],
                                pileup_data.coverages.to_vec(),
                                vec![&ploidy; n_samples],
                                &germ_genotyper,
                            ))
                            .unwrap();
                        continue;
                    }
                };

                let clustered_haps = collapse_haplotypes(
                    cluster_result,
                    pileup_data.haplos,
                    vec![gts.genotype.observed_alleles.clone(); n_samples],
                );

                debug!("Haps: {:#?}", clustered_haps);

                let should_build = !clustered_haps.is_empty()
                    && !m_args.graph.one_to_one
                    && m_graph.node_indices.len() <= (m_args.graph.maxnodes + 2);
                m_graph.build(should_build);

                // An id is put on the pathscore.meta so we can still tie it
                // back to the gt.genotype.observed_alleles
                let paths: Vec<PathScore> = clustered_haps
                    .into_iter()
                    .map(|h| m_graph.apply_haplotype(&h, &m_args.graph))
                    .filter(|p| *p != PathScore::default())
                    .collect();

                // Now, for each sample, separate the germline from the somatic
                let separated_paths = polycluster::separate_paths_by_sample(paths, n_samples);
                let (germline_paths, somatic_paths) =
                    separate_paths_by_vaf(separated_paths, gts.genotype);

                // And then m_graph.take_annotated on the germline paths
                let mut send_back = m_graph.take_annotated(
                    germline_paths.iter().map(|bin| bin.as_slice()).collect(),
                    pileup_data.coverages.to_vec(),
                    vec![&ploidy; n_samples],
                    &germ_genotyper,
                );

                // Before updating the somatic in place
                if let Some(ref mut send_back) = send_back {
                    // VCF entry
                    for (_record, annos) in &mut *send_back {
                        // Sample columns
                        for (sample_idx, (anno, sample_paths)) in
                            annos.iter_mut().zip(&somatic_paths).enumerate()
                        {
                            let mut any_update = false;
                            // Check all paths this sample had reads in
                            for path in sample_paths {
                                if path.path.contains(&anno.var_idx) {
                                    // This sample will need to be updated
                                    any_update = true;
                                    if !anno.gt.contains("1") {
                                        anno.gq = gts.quality_score.round() as i32;
                                    }
                                    anno.rnames.extend(path.meta.rnames.clone());
                                    anno.filt |= FiltFlags::SOMATIC;

                                    // TODO: wrong for haploid regions
                                    *anno.ad[1].get_or_insert(0) +=
                                        path.meta.coverage[sample_idx] as i32;
                                    // Logic here is that the germline reference allele
                                    // coverage is inflated during GenotypeAnno calculation
                                    // so we're moving it to the alternate as part of the
                                    // somatic support
                                    *anno.ad[0].get_or_insert(0) = anno.ad[0]
                                        .unwrap_or(0)
                                        .saturating_sub(path.meta.coverage[sample_idx] as i32);
                                }
                            }

                            if any_update {
                                // Redo the genotyping with the new coverage
                                if let (Some(rcov), Some(acov)) = (anno.ad[0], anno.ad[1]) {
                                    let ngt = germ_genotyper.genotype(rcov as u64, 0, acov as u64);
                                    anno.sq = ngt.sq as i32;
                                    anno.gq = ngt.gq.round() as i32;
                                    anno.gt_state = ngt.state;
                                    anno.gt = match anno.gt_state {
                                        GTstate::Ref => "0/0".to_string(),
                                        GTstate::Het => "0/1".to_string(),
                                        GTstate::Hom => "1/1".to_string(),
                                        GTstate::Non => "./.".to_string(),
                                    };
                                }
                            }
                        }
                    }
                }

                m_result_sender.send(send_back).unwrap();
            }
        }
    }
    // This should give a result
}

#[derive(Parser, Debug, Clone)]
pub struct MosaicCommand {
    #[command(flatten)]
    pub io: IOParams,

    #[command(flatten)]
    pub graph: GraphParams,

    /// Minimum number of reads in a cluster
    #[arg(long, default_value_t = 2, help_heading = "Genotyping")]
    pub msmin: usize,

    /// Max clusters
    #[arg(long, default_value_t = 8, help_heading = "Genotyping")]
    pub maxclust: usize,

    /// Clustering weight for haplotagged reads (off=0.0, full=1.0)
    #[arg(long, default_value_t = 0.25, help_heading = "Genotyping")]
    pub hps_weight: f32,

    /// Clustering weight for haplotype lengths (off=0.0, full=1.0)
    #[arg(long, default_value_t = 0.25, help_heading = "Genotyping")]
    pub len_weight: f32,

    /// Minimum haplotype size difference for K estimation
    #[arg(long, default_value_t = 2, help_heading = "Genotyping")]
    pub bandwidth: usize,

    /// VAF prior alpha
    #[arg(long, default_value_t = 1.0, help_heading = "Genotyping")]
    pub alpha: f64,

    /// VAF prior beta
    #[arg(long, default_value_t = 15.0, help_heading = "Genotyping")]
    pub beta: f64,

    /// Max somatic VAF
    #[arg(long, default_value_t = 0.20, help_heading = "Genotyping")]
    pub soma_vaf: f64,
}

impl ToPolyCluParams for MosaicCommand {
    fn to_polyclu_params(&self) -> polycluster::PolyCluParams {
        polycluster::PolyCluParams {
            msmin: self.msmin,
            maxclust: self.maxclust,
            hps_weight: self.hps_weight,
            len_weight: self.len_weight,
            minkfreq: self.graph.minkfreq,
            bandwidth: Some(self.bandwidth as f64),
            ..Default::default() // Fill remaining fields with defaults
        }
    }
}

#[derive(clap::Args, Clone, Debug)]
pub struct IOParams {
    /// VCF to genotype
    #[arg(short, long, help_heading = "I/O")]
    pub input: PathBuf,

    /// Reads to genotype (indexed .bam, .cram, or .plup.gz; can be specified multiple times).
    #[arg(long, help_heading = "I/O", action = ArgAction::Append)]
    pub reads: Vec<PathBuf>,

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

    /// Output VCF sample names (one per `--reads`; can be specified multiple times)
    #[arg(long, default_value = "SAMPLE", help_heading = "I/O", action = ArgAction::Append)]
    pub sample: Vec<String>,

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

impl KanpigCommand for MosaicCommand {
    fn debug(&self) -> bool {
        self.io.debug
    }
    /// Validate command line arguments
    fn validate(&self) -> bool {
        let mut is_ok = true;

        // Per-Bam
        is_ok &= file_validators::validate_file(&self.io.input, "--input");
        if self.io.reads.is_empty() {
            warn!("No `--reads` provided");
            is_ok = false;
        }
        for i in self.io.reads.iter() {
            is_ok &= file_validators::validate_reads(i, &self.graph);
        }

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

        if self.graph.maxpaths < 1 {
            error!("--maxpaths must be at least 1");
            is_ok = false;
        }

        if self.io.threads < 1 {
            error!("--threads must be at least 1");
            is_ok = false;
        }

        if self.maxclust < 2 {
            error!("--maxclust must be at least 2");
            is_ok = false;
        }

        is_ok
    }

    fn run(&mut self) {
        let mut input_vcf = vcf::io::reader::Builder::default()
            .build_from_path(self.io.input.clone())
            .expect("Unable to parse vcf");

        let input_header = input_vcf.read_header().expect("Unable to parse vcf header");

        // Validate --reads / --samples
        if self.io.sample.len() != self.io.reads.len() {
            if self.io.sample.is_empty() {
                self.io.sample = (0..self.io.reads.len())
                    .map(|i| format!("SAMPLE{}", i))
                    .collect();
                info!("Setting samples to {}", self.io.sample.join(", "));
            } else {
                error!(
                    "Expected one --sample for each --read, got {} and {}",
                    self.io.sample.len(),
                    self.io.reads.len()
                );
                std::process::exit(1);
            }
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
            self.io.rnames.clone(),
            self.io.sample.clone(),
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
