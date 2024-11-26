use crate::kplib::PlupArgs;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read, Record};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{mpsc, Arc};
use std::thread;

type ProcessedRead = (i32, i64, i64, String);

fn process_read(record: Record, sizemin: u32, sizemax: u32) -> ProcessedRead {
    let chrom = record.tid();
    let start = record.reference_start();
    let end = record.reference_end();

    // Build the CIGAR operations string
    let mut cigar_ops = Vec::new();
    let mut read_offset = 0;
    let mut align_offset = 0;

    for cigar in record.cigar().iter() {
        match cigar.char() {
            // Operations that consume read sequence and affect the offset
            'I' if sizemin <= cigar.len() && cigar.len() <= sizemax => {
                // Get the inserted sequence
                let sequence = String::from_utf8_lossy(
                    &record.seq().as_bytes()[read_offset..read_offset + cigar.len() as usize],
                )
                .to_string();
                cigar_ops.push(format!("{}:{}", align_offset, sequence));
                read_offset += cigar.len() as usize;
            }
            'D' if sizemin <= cigar.len() && cigar.len() <= sizemax => {
                // Record deletion without adjusting the offset
                cigar_ops.push(format!("{}:{}", align_offset, cigar.len()));
                align_offset += cigar.len() as usize;
            }
            'M' | 'X' | '=' => {
                read_offset += cigar.len() as usize;
                align_offset += cigar.len() as usize;
            }
            'I' | 'S' => {
                read_offset += cigar.len() as usize;
            }
            'D' => {
                align_offset += cigar.len() as usize;
            }
            // Handle 'H' and 'P' explicitly to ignore them
            'H' | 'P' => {}
            _ => {
                error!("What is this code?: {}", cigar.char());
            }
        }
    }

    if cigar_ops.is_empty() {
        (chrom, start, end, ".".to_string())
    } else {
        (chrom, start, end, cigar_ops.join(","))
    }
}

pub fn plup_main(args: PlupArgs) {
    let mut bam = bam::Reader::from_path(&args.bam).expect("Error opening BAM file");

    if let Some(ref_name) = args.reference {
        let _ = bam.set_reference(ref_name);
    }

    let (tx, rx): (mpsc::Sender<ProcessedRead>, mpsc::Receiver<ProcessedRead>) = mpsc::channel();

    // Spawn a thread for writing to ensure reads are written in order
    let writer_thread = {
        thread::spawn(move || {
            let mut writer: Box<dyn Write> = match args.output {
                Some(ref path) => {
                    let m_page = page_size::get() * 1000;
                    let file = File::create(path).expect("Error Creating Output File");
                    Box::new(BufWriter::with_capacity(m_page, file))
                }
                None => Box::new(BufWriter::new(std::io::stdout())),
            };

            let lbam = bam::Reader::from_path(args.bam).expect("Error opening BAM file");
            let header_main = bam::Header::from_template(lbam.header());
            let header = bam::HeaderView::from_header(&header_main);
            for result in rx {
                let chrom = std::str::from_utf8(header.tid2name(result.0 as u32)).expect("is okay");
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    chrom, result.1, result.2, result.3
                )
                .expect("Error writing to output file");
            }
        })
    };

    // Configure a custom thread pool
    let pool = ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("Failed to build thread pool");

    pool.install(|| {
        let mut buffer = Vec::with_capacity(args.batch_size);
        let tx = Arc::new(tx.clone()); // Cloneable sender for parallel tasks

        for record in bam.records().filter_map(|result| {
            match result {
                Ok(record) => {
                    // Filter out records that are empty or have the unwanted flags
                    if !record.seq().is_empty() && (record.flags() & args.mapflag) == 0 {
                        Some(record) // Only keep the valid records
                    } else {
                        None // Filter out invalid records
                    }
                }
                Err(_) => None, // Filter out any errors
            }
        }) {
            buffer.push(record);

            // When the buffer is full, process it in parallel
            if buffer.len() == args.batch_size {
                let chunk = std::mem::take(&mut buffer); // Take the batch
                let tx = Arc::clone(&tx);
                rayon::spawn(move || {
                    for record in chunk {
                        let _ = tx.send(process_read(record, args.sizemin, args.sizemax));
                    }
                });
            }
        }

        // Process any remaining records in the buffer
        if !buffer.is_empty() {
            let tx = Arc::clone(&tx);
            rayon::spawn(move || {
                for record in buffer {
                    let _ = tx.send(process_read(record, args.sizemin, args.sizemax));
                }
            });
        }
    });

    // Close the writer thread
    drop(tx);
    writer_thread.join().expect("Error joining writer thread");
}
