use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use ndarray::Array2;
use noodles_vcf::{self as vcf};
use rand::SeedableRng;
use std::{
    path::PathBuf,
    sync::{Arc, Mutex},
    thread::{self, JoinHandle},
};

use crate::{
    commands::KanpigCommand,
    file_validators,
    kplib::{
        build_region_tree, metrics, open_reads, open_writer_thread, ChannelInput, ChannelOutput,
        Haplotype, KDParams, PathScore, Ploidy, PloidyRegions, Variants, VcfChunker,
    },
};

fn select_k_by_ratio(losses: &[(usize, f32, f64)]) -> usize {
    if losses.len() < 3 {
        return losses.first().map(|(k, _, _)| *k).unwrap_or(1);
    }

    let mut ratios = Vec::new();
    for i in 1..losses.len() {
        let prev = losses[i - 1].1;
        let curr = losses[i].1;
        ratios.push(prev / curr);
    }

    for i in 1..ratios.len() {
        if ratios[i - 1] / ratios[i] > 2.0 {
            // Ratio dropped significantly (tune the 2.0 threshold if needed)
            return losses[i].0;
        }
    }

    // Fall back: choose last K
    losses.last().unwrap().0
}

// Input: Vec<(k: usize, loss: f64, silhouette: f64)>
fn select_optimal_k(metrics: &[(usize, f32, f64)]) -> usize {
    if metrics.len() < 3 {
        return metrics.first().map(|(k, _, _)| *k).unwrap_or(1);
    }

    // Compute loss deltas and silhouette differences
    let mut loss_deltas = Vec::new();
    for i in 1..metrics.len() {
        let delta = metrics[i - 1].1 - metrics[i].1;
        loss_deltas.push(delta);
    }

    // Normalize loss deltas (to detect where diminishing returns begin)
    let max_delta = loss_deltas[0];
    let elbow_k = loss_deltas
        .iter()
        .enumerate()
        .find(|(_, delta)| **delta < 0.25 * max_delta) // 25% cutoff
        .map(|(i, _)| metrics[i + 1].0) // i+1 because delta is from i to i+1
        .unwrap_or(metrics.last().unwrap().0);

    // Also find the max silhouette score's K
    let max_sil_k = metrics
        .iter()
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
        .map(|(k, _, _)| *k)
        .unwrap_or(elbow_k);

    // Combine: pick the smaller of elbow and max silhouette to avoid overfitting
    elbow_k.min(max_sil_k)
}

fn is_valid_k(
    counts: &std::collections::HashMap<usize, usize>,
    total: usize,
    min_reads: usize,
) -> (bool, usize) {
    //let min_allowed = (total as f32 * min_frac).ceil() as usize;
    let num_small_clusters = counts.values().filter(|&&v| v < min_reads).count();
    (num_small_clusters == 0, num_small_clusters)
}

fn task_thread(
    m_args: TrioCommand,
    m_receiver: Receiver<ChannelInput>,
    m_result_sender: Sender<ChannelOutput>,
    m_ploidy: PloidyRegions,
) {
    let mut pro_reads = open_reads(
        m_args.io.proband.clone(),
        m_args.io.reference.clone(),
        m_args.io.proband_sample.clone(),
        0, // First sample is index 0 in the HaplotypeMeta vectros
        3, // One total sample will be opened (for HaplotypeMeta)
        &m_args.kd,
    );
    let mut mat_reads = open_reads(
        m_args.io.mother.clone(),
        m_args.io.reference.clone(),
        m_args.io.mother_sample.clone(),
        1, // First sample is index 0 in the HaplotypeMeta vectros
        3, // One total sample will be opened (for HaplotypeMeta)
        &m_args.kd,
    );
    let mut pat_reads = open_reads(
        m_args.io.father.clone(),
        m_args.io.reference.clone(),
        m_args.io.father_sample.clone(),
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
                //let pro_haps = ploidy.cluster(pro_haps, pro_coverage, 0, &m_args.kd);

                let (mat_haps, mat_coverage) =
                    mat_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                //let mat_haps = ploidy.cluster(mat_haps, mat_coverage, 1, &m_args.kd);

                let (pat_haps, pat_coverage) =
                    pat_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                //let pat_haps = ploidy.cluster(pat_haps, pat_coverage, 2, &m_args.kd);

                let haplos: Vec<Haplotype> = pro_haps
                    .into_iter()
                    .chain(mat_haps)
                    .chain(pat_haps)
                    .collect();

                let dist: Array2<f32> =
                    Array2::from_shape_fn((haplos.len(), haplos.len()), |(i, j)| {
                        1.0 - (metrics::seqsim(
                            &haplos[i].kfeat,
                            &haplos[j].kfeat,
                            m_args.kd.minkfreq as f32,
                        ))
                    });

                let results: Vec<(usize, f32, f64)> = (1..(6.min(haplos.len())))
                    .map(|i| {
                        let k = i as usize;
                        let mut medoids = kmedoids::random_initialization(
                            haplos.len(),
                            k,
                            &mut rand::rngs::StdRng::seed_from_u64(21),
                        );

                        let (loss, assignments, _, _): (f32, _, _, _) =
                            kmedoids::fasterpam(&dist.view(), &mut medoids, 100);
                        let (sil2, _): (f64, _) =
                            kmedoids::medoid_silhouette(&dist, &medoids, false);
                        let (sil, _): (f64, _) = kmedoids::silhouette(&dist, &assignments, false);

                        let mut counts = std::collections::HashMap::<usize, usize>::new();
                        for item in assignments {
                            *counts.entry(item).or_insert(0) += 1;
                        }

                        let (is_valid, num_small_clusters) = is_valid_k(&counts, haplos.len(), 3);
                        debug!(
                            "K {}: Loss = {}, Sil = {}, Sil2 = {}, Small = {} {}",
                            k, loss, sil, sil2, is_valid, num_small_clusters
                        );
                        debug!("Assign: {:#?}", counts);

                        // TODO: This is good for this example, but we could just have a spurious
                        // read.. So you got to come up with another way to look at this
                        //let penalty = num_small_clusters as f64 / haplos.len as f64;
                        //let adjusted_silhouette = sil * (1.0 - penalty);
                        (k, loss, sil)
                    })
                    .collect();

                let opt = select_optimal_k(results.as_slice());
                let rbe = select_k_by_ratio(results.as_slice());
                debug!("Selected {} from optimal {} from RBE", opt, rbe);
                continue;
                //let pat_haps = ploidy.cluster(all_haps, pat_coverage, 2, &m_args.kd);
                // Only need to build the full graph sometimes
                //let should_build = !pro_haps.is_empty()
                //&& !m_args.kd.one_to_one
                //&& m_graph.node_indices.len() <= (m_args.kd.maxnodes + 2);
                m_graph.build(true);

                // Haplotypes to PathScores
                let paths: Vec<PathScore> = pro_haps
                    .into_iter()
                    .chain(mat_haps)
                    .chain(pat_haps)
                    .map(|h| m_graph.apply_haplotype(&h, &m_args.kd))
                    .filter(|p| *p != PathScore::default())
                    .collect();

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

                // I don't like this, maybe refactor take_annotated?
                //bin.sort_by(|a, b| hp_sorter(&a.meta.hp[0], &b.meta.hp[0]);
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

    /// Output VCF proband sample name
    #[arg(long, default_value = "PRO", help_heading = "I/O")]
    pub proband_sample: String,

    /// Output VCF maternal sample name
    #[arg(long, default_value = "MAT", help_heading = "I/O")]
    pub mother_sample: String,

    /// Output VCF paternal sample name
    #[arg(long, default_value = "PAT", help_heading = "I/O")]
    pub father_sample: String,

    // XYploidy_bed
    // XXploidy_bed
    // proband_karyotype XY or XX, which will then just point to whatever ploidy bed
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

        info!(
            "Setting samples to {}, {}, {}",
            self.io.proband_sample, self.io.mother_sample, self.io.father_sample
        );

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
                self.io.proband_sample.clone(),
                self.io.mother_sample.clone(),
                self.io.father_sample.clone(),
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
