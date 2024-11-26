use crate::kplib::PlupArgs;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{self, Read, Record};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{mpsc, Arc};
use std::thread;

type ProcessedRead = Option<(i32, i64, i64, String)>;
fn process_read(record: Record) -> ProcessedRead {
    //if record.seq().is_empty() || (record.flags() & 3840) != 0 {
        //return None;
    //}
    let chrom = record.tid();
    let start = record.pos();
    let end = start + record.seq_len() as i64;
    
    // Build the CIGAR operations string
    let mut cigar_ops = Vec::new();
    let mut offset = 0;

    for cigar in record.cigar().iter() {
            match cigar.char() {
                // Operations that consume read sequence and affect the offset
                'I' if cigar.len() >= 50 => {
                    // Get the inserted sequence
                    let sequence = String::from_utf8_lossy(
                        &record.seq().as_bytes()[offset..offset + cigar.len() as usize],
                    )
                    .to_string();
                    cigar_ops.push(format!("{}:{}", offset, sequence));
                    offset += cigar.len() as usize;
                }
                'D' if cigar.len() >= 50 => {
                    // Record deletion without adjusting the offset
                    cigar_ops.push(format!("{}:{}", offset, cigar.len()));
                }
                'M' | 'X' | '=' | 'S' | 'I' => {
                    offset += cigar.len() as usize;
                }
                // Handle 'H' and 'P' explicitly to ignore them
                'H' | 'P' => {}
                _ => {}
            }
    }

    if cigar_ops.is_empty() {
        Some((chrom, start, end, ".".to_string()))
    } else {
        Some((chrom, start, end, cigar_ops.join(",")))
    }
}

pub fn plup_main(args: PlupArgs) {
    let mut bam = bam::Reader::from_path(&args.input).expect("Error opening BAM file");
    let BATCH_SIZE = 100000; // 100k ~3.5G RSS
    let (tx, rx): (mpsc::Sender<ProcessedRead>, mpsc::Receiver<ProcessedRead>) = mpsc::channel();

    // Spawn a thread for writing to ensure reads are written in order
    let writer_thread = {
        thread::spawn(move || {
            let output_file = File::create(args.output).expect("Error creating output file");
            let mut writer = BufWriter::new(output_file);
            let lbam = bam::Reader::from_path(args.input).expect("Error opening BAM file");
            let header_main = bam::Header::from_template(lbam.header());
            let header = bam::HeaderView::from_header(&header_main);
            for result in rx {
                let line = result.expect("Error processing read");
                let chrom = std::str::from_utf8(header.tid2name(line.0 as u32)).expect("is okay");
                writeln!(writer, "{}\t{}\t{}\t{}", chrom, line.1, line.2, line.3)
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
            let mut buffer = Vec::with_capacity(BATCH_SIZE);
            let tx = Arc::new(tx.clone()); // Cloneable sender for parallel tasks

            for record in bam.records()
                .filter_map(|result| {
                    match result {
                        Ok(record) => {
                            // Filter out records that are empty or have the unwanted flags
                            if !record.seq().is_empty() && (record.flags() & 3840) == 0 {
                                Some(record) // Only keep the valid records
                            } else {
                                None // Filter out invalid records
                            }
                        }
                        Err(_) => None, // Filter out any errors
                    }
                })
            {
                buffer.push(record);

                // When the buffer is full, process it in parallel
                if buffer.len() == BATCH_SIZE {
                    let chunk = std::mem::take(&mut buffer); // Take the batch
                    let tx = Arc::clone(&tx);
                    rayon::spawn(move || {
                        for record in chunk {
                            if let Some(data) = process_read(record) {
                                tx.send(Some(data)).expect("Error sending data to channel");
                            }
                        }
                    });
                }
            }

            // Process any remaining records in the buffer
            if !buffer.is_empty() {
                let tx = Arc::clone(&tx);
                rayon::spawn(move || {
                    for record in buffer {
                        if let Some(data) = process_read(record) {
                            tx.send(Some(data)).expect("Error sending data to channel");
                        }
                    }
                });
            }
        });

    // Close the writer thread
    drop(tx);
    writer_thread.join().expect("Error joining writer thread");
}
