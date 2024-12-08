use std::{
    fs::File,
    io::{self, BufRead},
    path::Path,
};

type FileHandler = io::Result<io::Lines<io::BufReader<File>>>;

fn read_lines<P>(filename: P) -> FileHandler
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

#[derive(Debug)]
pub struct BedEntry {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub data: Option<Vec<String>>,
}

impl BedParser {
    pub fn new(path: &Path) -> Self {
        Self {
            file: path.to_path_buf(),
            prev_chrom: String::new(),
            prev_start: 0,
        }
    }

    pub fn parse(&mut self) -> Vec<BedEntry> {
        if let Ok(lines) = read_lines(&mut self.file) {
            lines
                .map_while(Result::ok)
                .map(|line| {
                    let collection: Vec<&str> = line.split('\t').collect();
                    if collection.len() < 3 {
                        error!("malformed bed line: {}", line);
                        std::process::exit(1);
                    }
                    let chrom = collection[0].to_string();
                    let start = collection[1].parse::<u64>().unwrap();
                    let end = collection[2].parse::<u64>().unwrap();
                    let data = if collection.len() >= 4 {
                        Some(collection[3..].iter().map(|x| x.to_string()).collect())
                    } else {
                        None
                    };

                    if chrom != self.prev_chrom {
                        self.prev_chrom.clone_from(&chrom);
                        self.prev_start = 0;
                    }

                    if end <= start {
                        error!("malformed bed line: stop <= start {}", line);
                        std::process::exit(1);
                    }
                    if start < self.prev_start {
                        error!(
                            "bed file unordered `sort -k3n -k1,2n` offending line {}",
                            line
                        );
                        std::process::exit(1);
                    }
                    BedEntry {
                        chrom,
                        start,
                        end,
                        data,
                    }
                })
                .collect()
        } else {
            panic!("unable to read bed file");
        }
    }
}
