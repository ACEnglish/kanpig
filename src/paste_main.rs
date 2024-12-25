use crate::kplib::PasteArgs;
use rust_htslib::bgzf;
use std::{
    fs::File,
    io::{self, BufWriter, Read, Write},
};

struct PlainVCF<R> {
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
pub fn paste_main(args: PasteArgs) -> Result<(), Box<dyn std::error::Error>> {
    // Should be able to make bgzip writers, also depending..
    let mut output: Box<dyn Write> = match args.output {
        Some(ref path) => {
            let m_page = page_size::get() * 1000;
            let file = File::create(path).expect("Error creating output file");
            Box::new(BufWriter::with_capacity(m_page, file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut readers: Vec<_> = args
        .inputs
        .iter()
        .map(|path| {
            let reader = bgzf::Reader::from_path(path).expect("Unable to read");
            PlainVCF::new(reader)
        })
        .collect();

    output.write_all(readers[0].get_header()?.join("\n").as_bytes())?;
    output.write_all(b"\n")?;

    for i in &mut readers[1..] {
        let _ = i.get_header()?;
    }

    while let Some(line) = readers[0].read_line()? {
        output.write_all(line.as_bytes())?;
        for i in &mut readers[1..] {
            let Some(sample) = i.get_sample()? else {
                panic!("bad line at {}", line)
            };
            output.write_all(b"\t")?;
            output.write_all(sample.as_bytes())?;
        }
        output.write_all(b"\n")?;
    }

    Ok(())
}
