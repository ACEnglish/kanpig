use crate::kplib::{PlupArgs, ReadPileup};
use crossbeam_channel::{unbounded, Receiver, Sender};
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{self, IndexedReader, Read},
};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::thread;
use std::thread::JoinHandle;

type InputType = Option<(String, u64, u64)>;
type OutputType = Option<Vec<ReadPileup>>;

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
                        && (rec.reference_start().unsigned_abs() >= start)
                        && (rec.reference_start().unsigned_abs() < end)
                })
            })
            .map(|record| ReadPileup::new(record, params.sizemin, params.sizemax))
            .collect()
    )
}

fn split_into_regions(bam_path: &PathBuf, chunk_size: usize) -> Vec<(String, u64, u64)> {
    let reader = IndexedReader::from_path(bam_path).expect("Failed to open BAM file");
    let header = reader.header().to_owned();

    (0..header.target_count())
        .flat_map(|tid| {
            let target_name = String::from_utf8(header.tid2name(tid as u32).to_vec())
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
    info!("{} regions to process", regions.len());

    // Create channels for communication between threads
    let (task_sender, task_receiver): (Sender<InputType>, Receiver<InputType>) = unbounded();
    let (result_sender, result_receiver): (Sender<OutputType>, Receiver<OutputType>) = unbounded();

    let write_handler = {
        let m_args = args.clone();
        thread::spawn(move || {
            let mut writer: Box<dyn Write> = match m_args.output {
                Some(ref path) => {
                    let m_page = page_size::get() * 1000;
                    let file = File::create(path).expect("Error Creating Output File");
                    Box::new(BufWriter::with_capacity(m_page, file))
                }
                None => Box::new(BufWriter::new(std::io::stdout())),
            };
            let lbam = bam::Reader::from_path(m_args.bam).expect("Error opening BAM file");
            let header_main = bam::Header::from_template(lbam.header());
            let header = bam::HeaderView::from_header(&header_main);
            let mut n_reads = 0;
            loop {
                match result_receiver.recv() {
                    Ok(None) | Err(_) => break,
                    Ok(Some(readplups)) => {
                        for read in readplups {
                            let chrom = std::str::from_utf8(header.tid2name(read.chrom as u32))
                                .expect("is okay");
                            writeln!(writer, "{}", read.to_string(chrom))
                                .expect("Error writing to output file");
                            n_reads += 1;
                        }
                    }
                }
            }
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
