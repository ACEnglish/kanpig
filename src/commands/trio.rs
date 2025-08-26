use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use ndarray::{Array, Array2};
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
        build_region_tree, hp_sorter, metrics, open_reads, open_writer_thread, trio_genotyper,
        ChannelInput, ChannelOutput, Haplotype, KDParams, MeanShift, PathScore, Ploidy,
        PloidyRegions, Variants, VcfChunker,
    },
};

fn cluster_quality(k: usize, quality_scores: &[f64], labels: &[usize]) -> Vec<f64> {
    // Find the maximum label to determine the size needed
    let mut scores = vec![0.0; k + 1]; // +2 because we need indices 0 through max_label+1

    // Reference cluster is always 0 and quality of 1
    scores[0] = 1.0;

    // Get unique labels and calculate mean for each cluster
    let mut unique_labels: Vec<usize> = labels.to_vec();
    unique_labels.sort_unstable();
    unique_labels.dedup();

    for cluster in unique_labels {
        let cluster_scores: Vec<f64> = quality_scores
            .iter()
            .zip(labels.iter())
            .filter(|(_, &label)| label == cluster)
            .map(|(&score, _)| score)
            .collect();

        if !cluster_scores.is_empty() {
            let mean = cluster_scores.iter().sum::<f64>() / cluster_scores.len() as f64;
            scores[cluster + 1] = mean;
        }
    }

    scores
}

/// Count reads in each cluster
fn count_reads(k: usize, ref_coverage: &[usize; 3], assignments: &[usize], haplos: &[Haplotype]) -> Array2<usize> {
    let mut read_counts = Array::<usize, _>::zeros((k + 1, 3));
    read_counts[[0, 0]] = ref_coverage[0];
    read_counts[[0, 1]] = ref_coverage[1];
    read_counts[[0, 2]] = ref_coverage[2];
    for (&label, hap) in assignments.iter().zip(haplos.iter()) {
        let cluster_idx = label + 1;
        let sample_idx = hap.meta.samples_flag.trailing_zeros() as usize;
        read_counts[[cluster_idx, sample_idx]] += 1;
    }
    read_counts
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

    let mut pat_reads = open_reads(
        m_args.io.father.clone(),
        m_args.io.reference.clone(),
        m_args.io.father_sample.clone(),
        1,
        3,
        &m_args.kd,
    );

    let mut mat_reads = open_reads(
        m_args.io.mother.clone(),
        m_args.io.reference.clone(),
        m_args.io.mother_sample.clone(),
        2,
        3,
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

                let (pat_haps, pat_coverage) =
                    pat_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                //let pat_haps = ploidy.cluster(pat_haps, pat_coverage, 1, &m_args.kd);

                let (mat_haps, mat_coverage) =
                    mat_reads.find_pileups(&m_graph.chrom, m_graph.start, m_graph.end);
                //let mat_haps = ploidy.cluster(mat_haps, mat_coverage, 2, &m_args.kd);

                let ref_coverage = [
                    pro_coverage as usize - pro_haps.len(),
                    pat_coverage as usize - pat_haps.len(),
                    mat_coverage as usize - mat_haps.len(),
                ];

                let haplos: Vec<Haplotype> = pro_haps
                    .into_iter()
                    .chain(pat_haps)
                    .chain(mat_haps)
                    .collect();

                // TODO: is_empty checks on the input haplotypes.
                if haplos.len() <= 1 {
                    continue;
                }

                // MeanShift to determine K
                let sizes: Vec<f64> = haplos.iter().map(|x| x.size as f64).collect();
                // TODO PARAM and dynamic figure this out
                let mut ms = MeanShift::new().min_size(m_args.minreads);
                let ms_result = ms.fit(&sizes);

                // KMedoid Clustering
                let k = ms_result.cluster_centers.len();
                debug!("Setting K to {:?}", k);
                // Use MeanShift centers to start the medoids
                let mut medoids = ms_result.medoids.clone();

                let (assignments, quality) = if m_args.lengthonly {
                    (ms_result.labels, vec![1.0; k + 1])
                } else {
                    // Kmedoid Clustering
                    let dist: Array2<f32> =
                        Array2::from_shape_fn((haplos.len(), haplos.len()), |(i, j)| {
                            let mut dist: f32 = 1.0
                                - (metrics::seqsim(
                                    &haplos[i].kfeat,
                                    &haplos[j].kfeat,
                                    m_args.kd.minkfreq as f32,
                                ));
                            // if same sample and different hp, hps_weight penalty
                            let i_samp = haplos[i].meta.samples_flag;
                            let j_samp = haplos[j].meta.samples_flag;
                            if i_samp == j_samp {
                                match (
                                    haplos[i].meta.hp[i_samp.trailing_zeros() as usize],
                                    haplos[j].meta.hp[j_samp.trailing_zeros() as usize],
                                ) {
                                    (Some(group_i), Some(group_j)) if group_i != group_j => {
                                        //TODO PARAM
                                        dist *= 1.25;
                                    }
                                    _ => (),
                                }
                            }
                            if ms_result.labels[i] != ms_result.labels[j] {
                                //TODO PARAM
                                dist *= 1.25;
                            }
                            dist
                        });

                    let (_loss, assignments, _, _): (f32, _, _, _) =
                        kmedoids::fasterpam(&dist.view(), &mut medoids, 100);
                    let (_, quality): (f64, Vec<f64>) =
                        kmedoids::medoid_silhouette(&dist, &medoids, true);
                    let quality = cluster_quality(k, &quality, &assignments);
                    (assignments, quality)
                };

                let read_counts = count_reads(k, &ref_coverage, &assignments, &haplos);

                debug!("Read Counts:\n {:?}", read_counts);
                // TODO: Should be using this GQ?
                let (gts, _gqs) = trio_genotyper(&read_counts, &quality);

                let mut clustered_haps: Vec<Haplotype> =
                    medoids.iter().map(|i| haplos[*i].clear_clone()).collect();

                let mut hp_cnt = Array::<u16, _>::zeros((k, 3, 2));

                // Collapse Haplotypes
                assignments
                    .into_iter()
                    .zip(haplos)
                    .for_each(|(cluster_idx, m_hap)| {
                        // Sample index inside the HaplotypeMeta
                        let idx = m_hap.meta.samples_flag.trailing_zeros() as usize;
                        // Only apply reads to the clustered_hap if it goes together
                        if gts[idx].contains(&(cluster_idx + 1)) {
                            let k_hap = &mut clustered_haps[cluster_idx];
                            k_hap.meta.coverage[idx] += 1;
                            k_hap.meta.ps[idx] = k_hap.meta.ps[idx].or(m_hap.meta.ps[idx]);
                            k_hap.meta.hp[idx] = k_hap.meta.hp[idx].or(m_hap.meta.hp[idx]);
                            k_hap.meta.samples_flag |= m_hap.meta.samples_flag;

                            if let Some(val) = m_hap.meta.hp[idx] {
                                hp_cnt[[cluster_idx, idx, val as usize - 1]] += 1;
                            }
                        }
                        // TODO: Reassignment of reads assigned to unused clusters?
                    });

                // HP tag for GT order
                for (i, m_hap) in clustered_haps.iter_mut().enumerate() {
                    for j in 0..3 {
                        if m_hap.meta.hp[j].is_some() {
                            let max_idx: u8 = hp_cnt
                                .slice(ndarray::s![i, j, ..])
                                .iter()
                                .cloned()
                                .enumerate()
                                .max_by_key(|&(_, val)| val)
                                .map(|(idx, _)| idx)
                                .unwrap_or(1)
                                .try_into()
                                .unwrap();
                            m_hap.meta.hp[j] = Some(max_idx + 1);
                        }
                    }
                }

                let should_build = !clustered_haps.is_empty()
                    && !m_args.kd.one_to_one
                    && m_graph.node_indices.len() <= (m_args.kd.maxnodes + 2);
                m_graph.build(should_build);

                // Haplotypes to PathScores
                let paths: Vec<PathScore> = clustered_haps
                    .into_iter()
                    .map(|h| m_graph.apply_haplotype(&h, &m_args.kd))
                    .filter(|p| *p != PathScore::default())
                    .collect();

                // Separate paths back out to the samples
                let num_samples = 3;
                let mut separated_paths: Vec<Vec<PathScore>> = vec![Vec::new(); num_samples];
                for path in paths {
                    for (bit, s_paths) in separated_paths.iter_mut().enumerate().take(num_samples) {
                        if (path.meta.samples_flag & (1 << bit)) != 0 {
                            let mut p = path.clone();
                            p.meta.samples_flag = bit;
                            s_paths.push(p);
                        }
                    }
                }
                // HP Ordering
                let separated_paths: Vec<&[PathScore]> = separated_paths
                    .iter_mut()
                    .map(|bin| {
                        bin.sort_by(|a, b| {
                            hp_sorter(
                                &a.meta.hp[a.meta.samples_flag],
                                &b.meta.hp[b.meta.samples_flag],
                            )
                        });
                        bin.as_slice()
                    })
                    .collect();

                m_result_sender
                    .send(m_graph.take_annotated(
                        separated_paths,
                        vec![pro_coverage, pat_coverage, mat_coverage],
                        vec![&ploidy, &ploidy, &ploidy], // TODO: set this up for each
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

    /// Only cluster on haplotype lengths
    #[arg(long, default_value_t = false, help_heading = "Trio")]
    pub lengthonly: bool,

    /// Minimum number of reads in a cluster
    #[arg(long, default_value_t = 3, help_heading = "Trio")]
    pub minreads: usize,
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

    /// Output VCF proband sample name
    #[arg(long, default_value = "PRO", help_heading = "I/O")]
    pub proband_sample: String,

    /// Output VCF paternal sample name
    #[arg(long, default_value = "PAT", help_heading = "I/O")]
    pub father_sample: String,

    /// Output VCF maternal sample name
    #[arg(long, default_value = "MAT", help_heading = "I/O")]
    pub mother_sample: String,

    // TODO: XYploidy_bed
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
            self.io.proband_sample, self.io.father_sample, self.io.mother_sample
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
