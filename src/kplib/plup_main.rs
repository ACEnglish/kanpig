use crate::kplib::{PlupArgs, ReadPileup};
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{mpsc, Arc};
use std::thread;

pub fn plup_main(args: PlupArgs) {
    let mut bam = bam::Reader::from_path(&args.bam).expect("Error opening BAM file");

    if let Some(ref_name) = args.reference {
        let _ = bam.set_reference(ref_name);
    }

    let (tx, rx): (mpsc::Sender<ReadPileup>, mpsc::Receiver<ReadPileup>) = mpsc::channel();

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
                // Need to pass the chrom in
                let chrom =
                    std::str::from_utf8(header.tid2name(result.chrom as u32)).expect("is okay");
                writeln!(writer, "{}", result.to_string(chrom))
                    .expect("Error writing to output file");
            }
        })
    };

    let pool = ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .expect("Failed to build thread pool");

    pool.install(|| {
        let mut buffer = Vec::with_capacity(args.batch_size);
        let tx = Arc::new(tx.clone());

        for record in bam.records().filter_map(|result| match result {
            Ok(record) => {
                if record.seq().is_empty()
                    || record.mapq() < args.mapq
                    || (record.flags() & args.mapflag) != 0
                {
                    None
                } else {
                    Some(record)
                }
            }
            Err(_) => None,
        }) {
            buffer.push(record);

            if buffer.len() == args.batch_size {
                let chunk = std::mem::take(&mut buffer);
                let tx = Arc::clone(&tx);
                rayon::spawn(move || {
                    for record in chunk {
                        let _ = tx.send(ReadPileup::new(record, args.sizemin, args.sizemax));
                    }
                });
            }
        }

        // Process any remaining records in the buffer
        if !buffer.is_empty() {
            let tx = Arc::clone(&tx);
            rayon::spawn(move || {
                for record in buffer {
                    let _ = tx.send(ReadPileup::new(record, args.sizemin, args.sizemax));
                }
            });
        }
    });

    // Close the writer thread
    drop(tx);
    writer_thread.join().expect("Error joining writer thread");
}
