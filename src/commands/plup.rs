use crate::{
    commands::KanpigCommand,
    file_validators,
    kplib::{open_bam, pileup::ReadPileup, SequenceMeta},
};
use clap::Parser;
use crossbeam_channel::{unbounded, Receiver, Sender};
use serde::{Deserialize, Serialize};

use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{self, IndexedReader, Read},
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    thread::{self, JoinHandle},
};

type InputType = Option<(String, u64, u64)>;
type OutputType = Option<Vec<ReadPileup>>;

/// Processes a specified region in a BAM file, filtering reads based on user-defined parameters and returning the results.
///
/// # Parameters
/// - `reader`: A mutable reference to an `IndexedReader` for reading the BAM file. The reader must be initialized and have a valid index loaded.
/// - `chrom`: A string representing the chromosome or reference sequence name to query.
/// - `start`: The start position (inclusive) of the region to fetch, in 0-based coordinates.
/// - `end`: The end position (exclusive) of the region to fetch, in 0-based coordinates.
/// - `params`: A reference to a `PlupCommand` struct containing user-defined filtering criteria, including:
///     - `mapq`: Minimum mapping quality required for reads to be included.
///     - `mapflag`: Bitwise flags for filtering reads based on their SAM flag values.
///     - `sizemin`: Minimum size threshold for reads to be included in the pileup.
///     - `sizemax`: Maximum size threshold for reads to be included in the pileup.
///
/// # Returns
/// - `OutputType`: A collection of processed reads from the specified region that meet the filtering criteria, organized as defined in `OutputType`.
///
/// # Panics
/// This function panics if the `fetch` operation on the BAM reader fails, which can occur if the specified region is invalid or if there is an issue with the BAM file or its index.
///
/// # Example
/// ```rust
/// let mut reader = IndexedReader::from_path("example.bam").unwrap();
/// let params = PlupCommand {
///     mapq: 30,
///     mapflag: 0,
///     sizemin: 50,
///     sizemax: 500,
/// };
/// let chrom = String::from("chr1");
/// let start = 100_000;
/// let end = 200_000;
///
/// let result = process_bam_region(&mut reader, &chrom, start, end, &params);
/// // Process the result...
/// ```
fn process_bam_region(
    reader: &mut IndexedReader,
    chrom: &String,
    start: u64,
    end: u64,
    params: &PlupCommand,
) -> OutputType {
    reader
        .fetch((chrom, start, end))
        .expect("Failed to fetch region");

    let mut ret = Vec::with_capacity((params.chunk_size as usize) * 2000);
    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r.expect("Failed to parse record");
        if !record.seq().is_empty()
            && record.mapq() >= params.mapq
            && (record.flags() & params.mapflag) == 0
            && record.reference_start().unsigned_abs() >= start
            && record.reference_start().unsigned_abs() < end
        {
            ret.push(ReadPileup::new_record(
                chrom.clone(),
                &record,
                params.sizemin,
                params.sizemax,
                SequenceMeta::new(0, 1),
            ));
        }
    }
    Some(ret)
}

/// Splits the reference sequences in a BAM file into regions of a specified size.
///
/// # Parameters
/// - `bam_path`: A reference to a `PathBuf` pointing to the BAM file. The BAM file must be indexed.
/// - `chunk_size`: The size of each region, in base pairs. The last region for a reference sequence may be smaller if the reference length is not a multiple of `chunk_size`.
///
/// # Returns
/// - `Vec<(String, u64, u64)>`: A vector of tuples where each tuple contains:
///   - The name of the reference sequence (`String`).
///   - The start position (inclusive) of the region (`u64`), in 0-based coordinates.
///   - The end position (exclusive) of the region (`u64`), in 0-based coordinates.
///
/// # Panics
/// - This function panics if the BAM file cannot be opened or if its header contains invalid UTF-8.
/// - Panics if any reference sequence's length cannot be determined.
///
/// # Example
/// ```rust
/// use std::path::PathBuf;
///
/// let bam_path = PathBuf::from("example.bam");
/// let chunk_size = 1_000_000; // Split regions into 1 Mb chunks
///
/// let regions = split_into_regions(&bam_path, chunk_size);
/// for (chrom, start, end) in regions {
///     println!("{}:{}-{}", chrom, start, end);
/// }
/// ```
fn split_into_regions(bam_path: &PathBuf, chunk_size: usize) -> Vec<(String, u64, u64)> {
    let reader = open_bam(bam_path).expect("BAM already checked");
    let header = reader.header().to_owned();

    (0..header.target_count())
        .flat_map(|tid| {
            let target_name = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in target name");
            let target_len = header.target_len(tid).expect("Failed to get target length") as usize;

            (0..target_len).step_by(chunk_size).map(move |start| {
                let end = usize::min(start + chunk_size, target_len);
                (target_name.clone(), start as u64, end as u64)
            })
        })
        .collect()
}

#[derive(Parser, Serialize, Deserialize, Debug, Clone)]
pub struct PlupCommand {
    /// Input BAM/CRAM file
    #[arg(short, long)]
    pub bam: PathBuf,

    /// Reference file for CRAMs
    #[arg(short, long)]
    pub reference: Option<PathBuf>,

    /// Output plup (unsorted, uncompressed) [default: stdout]
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,

    /// Minimum size of variant to index
    #[arg(long, default_value_t = 50)]
    pub sizemin: u32,

    /// Maximum size of variant to index
    #[arg(long, default_value_t = 10000)]
    pub sizemax: u32,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Ignore alignments matching flag
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,

    /// Chunksize in Mbp
    #[arg(long, default_value_t = 25)]
    pub chunk_size: u64,

    /// Verbose logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,
}

impl KanpigCommand for PlupCommand {
    fn debug(&self) -> bool {
        self.debug
    }

    fn validate(&self) -> bool {
        let mut is_ok = true;

        is_ok &= file_validators::validate_bam(&self.bam);

        if let Some(ref_path) = &self.reference {
            is_ok &= file_validators::validate_reference(ref_path);
        }

        if self.sizemin < 20 {
            warn!("--sizemin is recommended to be at least 20.");
        }

        is_ok
    }

    /// Runner for kanpig plup command
    fn run(&mut self) {
        let regions = split_into_regions(&self.bam, (self.chunk_size as usize) * 1000000);
        let num_regions = regions.len() as u64;
        info!("{} regions to process", num_regions);

        // Create channels for communication between threads
        let (task_sender, task_receiver): (Sender<InputType>, Receiver<InputType>) = unbounded();
        let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) =
            unbounded();

        let write_handler = {
            let m_args = self.clone();
            let sty = ProgressStyle::with_template(
                " [{elapsed_precise}] {bar:44.cyan/blue} > {pos} completed",
            )
            .unwrap()
            .progress_chars("・🐷🥫");
            thread::spawn(move || {
                let mut writer: Box<dyn Write> = match m_args.output {
                    Some(ref path) => {
                        let m_page = page_size::get() * 1000;
                        let file = File::create(path).expect("Error Creating Output File");
                        Box::new(BufWriter::with_capacity(m_page, file))
                    }
                    None => Box::new(BufWriter::new(std::io::stdout())),
                };
                // Header
                let serialized = serde_json::to_string(&m_args).expect("Error writing header");
                let prefixed = format!("# {}\n", serialized);
                let _ = writer.write_all(prefixed.as_bytes());

                let mut n_reads = 0;
                let pbar = ProgressBar::new(num_regions).with_style(sty);
                pbar.inc(0);
                loop {
                    match result_receiver.recv() {
                        Ok(None) | Err(_) => break,
                        Ok(Some(readplups)) => {
                            for read in readplups {
                                writeln!(writer, "{}", read).expect("Error writing to output file");
                                n_reads += 1;
                            }
                            pbar.inc(1);
                        }
                    }
                }
                pbar.finish();
                info!("processed {} reads", n_reads);
            })
        };

        info!("spawning {} threads", self.threads);
        let task_handles: Vec<JoinHandle<()>> = (0..self.threads)
            .map(|_| {
                let m_args = self.clone();
                let m_receiver = task_receiver.clone();
                let m_result_sender = result_sender.clone();
                thread::spawn(move || {
                    let mut m_bam = open_bam(&m_args.bam).expect("BAM already checked");
                    if let Some(ref ref_name) = m_args.reference {
                        let _ = m_bam.set_reference(ref_name);
                    }
                    loop {
                        match m_receiver.recv() {
                            Ok(None) | Err(_) => break,
                            Ok(Some(chunk)) => {
                                let _ = m_result_sender.send(process_bam_region(
                                    &mut m_bam, &chunk.0, chunk.1, chunk.2, &m_args,
                                ));
                            }
                        }
                    }
                })
            })
            .collect();

        // Send items to worker threads
        for i in regions {
            task_sender.send(Some(i)).unwrap();
        }

        // Signal worker threads to exit
        for _ in 0..self.threads {
            task_sender.send(None).unwrap();
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
