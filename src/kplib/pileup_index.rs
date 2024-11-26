use crate::kplib::PlupArgs;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bgzf::Writer;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;

type ProcessedRead = Option<(i32, i64, i64, String)>;
fn process_read(record: Record) -> ProcessedRead {
    let chrom = record.tid();
    let start = record.pos();
    let end = start + record.seq_len() as i64;

    // Build the CIGAR operations string
    let mut cigar_ops = Vec::new();
    let mut offset = 0;

    for cigar in record.cigar().iter() {
        match cigar.char() {
            'M' | 'X' | '=' | 'S' | 'H' => {
                offset += cigar.len() as usize;
            }
            'I' if cigar.len() >= 50 => {
                // Get the inserted sequence
                let sequence = String::from_utf8_lossy(
                    &record.seq().as_bytes()[offset..offset + cigar.len() as usize],
                )
                .to_string();
                cigar_ops.push(format!("{}:{}", offset, sequence));
            }
            'D' if cigar.len() >= 50 => {
                cigar_ops.push(format!("{}:{}", offset, cigar.len()));
                offset += cigar.len() as usize;
            }
            _ => {}
        }
    }

    if cigar_ops.is_empty() {
        None
    } else {
        Some((chrom, start, end, cigar_ops.join(",")))
    }
}

pub fn plup_main(args: PlupArgs) {
    let mut bam = bam::Reader::from_path(&args.input).expect("Error opening BAM file");
    let (tx, rx): (mpsc::Sender<ProcessedRead>, mpsc::Receiver<ProcessedRead>) = mpsc::channel();

    // Spawn a thread for writing to ensure reads are written in order
    let writer_thread = {
        thread::spawn(move || {
            let mut writer = rust_htslib::bgzf::Writer::from_path(args.output).expect("Error creating BGZF writer");
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

    // Process BAM reads in parallel using the custom thread pool
    pool.install(|| {
        bam.records()
            .par_bridge()
            .for_each_with(tx.clone(), |tx, result| {
                if let Ok(record) = result {
                    if let Some(data) = process_read(record) {
                        tx.send(Some(data)).expect("Error sending data to channel");
                    }
                }
            });
    });

    // Close the writer thread
    drop(tx);
    writer_thread.join().expect("Error joining writer thread");
}
