use std::io::{self, Read};
use crate::kplib::PasteArgs;
use rust_htslib::bgzf;

pub struct PlainVCF<R> {
    reader: R,
    buffer: Vec<u8>,
    position: usize,
    leftover_line: Option<String>,
}

impl<R: Read> PlainVCF<R> {
    /// Creates a new PlainVCF wrapper around a reader
    pub fn new(reader: R) -> Self {
        PlainVCF {
            reader,
            buffer: Vec::new(),
            position: 0,
            leftover_line: None,
        }
    }

    /// Reads the next line from the reader
    pub fn read_line(&mut self) -> io::Result<Option<String>> {
        if let Some(line) = self.leftover_line.take() {
            return Ok(Some(line));
        }

        let mut line = Vec::new();

        loop {
            // Check if there is data in the buffer to process
            while self.position < self.buffer.len() {
                let byte = self.buffer[self.position];
                self.position += 1;

                if byte == b'\n' {
                    // End of line
                    return String::from_utf8(line)
                        .map(Some)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));
                } else {
                    line.push(byte);
                }
            }

            // If we reach here, the buffer is exhausted, so we refill it
            self.buffer.clear();
            self.position = 0;
            let mut temp_buf = [0; 1024];
            let bytes_read = self.reader.read(&mut temp_buf)?;

            if bytes_read == 0 {
                // End of file
                if line.is_empty() {
                    return Ok(None);
                } else {
                    return String::from_utf8(line)
                        .map(Some)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));
                }
            }

            self.buffer.extend_from_slice(&temp_buf[..bytes_read]);
        }
    }

    /// Reads all header lines starting with "##"
    pub fn get_header(&mut self) -> io::Result<Vec<String>> {
        let mut headers = Vec::new();

        while let Some(line) = self.read_line()? {
            if line.starts_with("##") {
                headers.push(line);
            } else {
                self.leftover_line = Some(line); // Store the first non-header line
                break;
            }
        }

        Ok(headers)
    }

    /// Gets the 10th+ columns (1-based index) of the next line, tab-delimited
    pub fn get_sample(&mut self) -> io::Result<Option<String>> {
        if let Some(line) = self.read_line()? {
            let mut columns = line.split('\t');

            // Skip the first 9 columns
            for _ in 0..9 {
                if columns.next().is_none() {
                    return Ok(Some(String::new()));
                }
            }

            // Collect the remaining columns as a single string
            let sample_columns: String = columns.collect::<Vec<&str>>().join("\t");
            return Ok(Some(sample_columns));
        }
        Ok(None)
    }
}

// Example usage
pub fn paste_main(_args: PasteArgs) -> Result<(), Box<dyn std::error::Error>> {
    let p1 = "/Users/english/code/truvari/repo_utils/test_files/variants/input1.vcf.gz".to_string();
    let data1 = b"##fileformat=VCFv4.2\n##source=MyTool\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2\nchr1\t12345\trs123456\tA\tT\t100\tPASS\t.\tGT\t0|1\t1|1\n";
    let data2 = b"##fileformat=VCFv4.2\n##source=MyTool\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE3\tSAMPLE4\nchr1\t12345\trs123456\tA\tT\t100\tPASS\t.\tGT\t0|1\t1|1\n";

    // let reader =  &data1[..];
    let reader = bgzf::Reader::from_path(&p1)?;
    let mut plain_vcf1 = PlainVCF::new(reader);

    //let reader2 = &data2[..];
    let reader2 = bgzf::Reader::from_path(&p1)?;
    let mut plain_vcf2 = PlainVCF::new(reader2);

    let headers = plain_vcf1.get_header()?;
    for header in headers {
        println!("{}", header);
    }
    let _ = plain_vcf2.get_header()?;

    while let Some(line) = plain_vcf1.read_line()? {
        let Some(sample) = plain_vcf2.get_sample()? else { panic!("bad line at {}", line) };

        print!("{}\t{}", line, sample);
        print!("\n");
    }

    Ok(())
}

