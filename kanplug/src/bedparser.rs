use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

type FileHandler = io::Result<io::Lines<io::BufReader<File>>>;

pub fn read_lines<P>(filename: P) -> FileHandler
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub struct BedParser {
    /* Read tab delimited bed files while ensuring entries have start < end.
    It also ensures entries are sorted .. except in the case of chrA to chrB
    back to chrA*/
    file: std::path::PathBuf,
    prev_chrom: String,
    prev_start: u64,
}

impl BedParser {
    pub fn new(path: &Path) -> Self {
        Self {
            file: path.to_path_buf(),
            prev_chrom: String::new(),
            prev_start: 0,
        }
    }

    pub fn parse(&mut self) -> Vec<(String, u64, u64)> {
        if let Ok(lines) = read_lines(&mut self.file) {
            lines
                .flatten()
                .map(|line| {
                    let collection: Vec<&str> = line.split('\t').collect();
                    if collection.len() < 3 {
                        error!("malformed bed line: {}", line);
                        std::process::exit(1);
                    }
                    let chrom = collection[0].to_string();
                    let m_start = collection[1].parse::<u64>().unwrap();
                    let m_stop = collection[2].parse::<u64>().unwrap();

                    if chrom != self.prev_chrom {
                        self.prev_chrom = chrom.clone();
                        self.prev_start = 0;
                    }

                    if m_stop <= m_start {
                        error!("malformed bed line: stop <= start {}", line);
                        std::process::exit(1);
                    }
                    if m_start < self.prev_start {
                        error!(
                            "bed file unordered `sort -k3n -k1,2n` offending line {}",
                            line
                        );
                        std::process::exit(1);
                    }
                    (chrom, m_start, m_stop)
                })
                .collect()
        } else {
            panic!("unable to read bed file");
        }
    }
}
