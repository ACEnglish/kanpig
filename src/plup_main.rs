use crate::kplib::{PlupArgs, ReadPileup};
use crossbeam_channel::{unbounded, Receiver, Sender};
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
/// - `params`: A reference to a `PlupArgs` struct containing user-defined filtering criteria, including:
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
/// let params = PlupArgs {
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
    params: &PlupArgs,
) -> OutputType {
    reader
        .fetch((chrom, start, end))
        .expect("Failed to fetch region");

    Some(
        reader
            .records()
            .filter_map(|r| {
                r.ok().filter(|rec| {
                    !rec.seq().is_empty()
                        && rec.mapq() >= params.mapq
                        && (rec.flags() & params.mapflag) == 0
                        && rec.reference_start().unsigned_abs() >= start
                        && rec.reference_start().unsigned_abs() < end
                })
            })
            .map(|record| ReadPileup::new(record, params.sizemin, params.sizemax))
            .collect(),
    )
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
    let reader = IndexedReader::from_path(bam_path).expect("Failed to open BAM file");
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

pub fn plup_main(args: PlupArgs) {
    let regions = split_into_regions(&args.bam, (args.chunk_size as usize) * 1000000);
    let num_regions = regions.len() as u64;
    info!("{} regions to process", num_regions);

    // Create channels for communication between threads
    let (task_sender, task_receiver): (Sender<InputType>, Receiver<InputType>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    let write_handler = {
        let m_args = args.clone();
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

            let bam = bam::Reader::from_path(m_args.bam).expect("Error opening BAM file");
            let header_main = bam::Header::from_template(bam.header());
            let header = bam::HeaderView::from_header(&header_main);
            let mut n_reads = 0;
            let pbar = ProgressBar::new(num_regions).with_style(sty);
            pbar.inc(0);
            loop {
                match result_receiver.recv() {
                    Ok(None) | Err(_) => break,
                    Ok(Some(readplups)) => {
                        for read in readplups {
                            let chrom = std::str::from_utf8(header.tid2name(read.chrom as u32))
                                .expect("Unable to lookup tid");
                            writeln!(writer, "{}", read.to_string(chrom))
                                .expect("Error writing to output file");
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

    info!("spawning {} threads", args.threads);
    let task_handles: Vec<JoinHandle<()>> = (0..args.threads)
        .map(|_| {
            let m_args = args.clone();
            let m_receiver = task_receiver.clone();
            let m_result_sender = result_sender.clone();
            thread::spawn(move || {
                let mut m_bam =
                    IndexedReader::from_path(&m_args.bam).expect("Failed to open BAM file");
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
    for _ in 0..args.threads {
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
