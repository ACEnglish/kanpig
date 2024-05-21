use indexmap::IndexMap;
use noodles_vcf::header::record::value::{map::contig::Name, map::Contig, Map};
use std::collections::HashMap;
use std::collections::VecDeque;

pub type ContigMap = IndexMap<Name, Map<Contig>>;
pub type Regions = HashMap<String, VecDeque<(u64, u64)>>;
use crate::kplib::BedParser;

/// create a HashMap with keys of chromsome names and
/// values a list of start, end positions with regions
/// which should be ananlyzed
pub fn build_region_tree(
    vcf_contigs: &ContigMap,
    includebed: &Option<std::path::PathBuf>,
) -> Regions {
    let mut m_contigs = HashMap::new();
    for (k, v) in vcf_contigs {
        let name = k.to_string();
        let length = match v.length() {
            Some(l) => l,
            None => {
                panic!("--vcf contig header missing length");
            }
        };
        m_contigs.insert(name, VecDeque::from([(0, length as u64)]));
    }

    if includebed.is_none() {
        return m_contigs;
    }

    // Parse the Bed Lines and return them as the IndexMap.
    let mut ret = HashMap::new();
    let mut m_parser = BedParser::new(&includebed.clone().unwrap());
    let mut prev_chrom = String::new();
    let mut prev_start: u64 = 0;

    for entry in m_parser
        .parse()
        .into_iter()
        .filter(|entry| m_contigs.contains_key(&entry.chrom))
    {
        if entry.chrom != prev_chrom {
            prev_chrom = entry.chrom.clone();
            prev_start = 0;
        }

        if entry.end <= entry.start {
            error!("malformed bed line: end <= start @ {:?}", entry);
            std::process::exit(1);
        }

        if entry.end < prev_start {
            error!("bed file unordered `sort -k3n -k1,2n` @ {:?}", entry);
            std::process::exit(1);
        }

        ret.entry(entry.chrom)
            .or_insert_with(VecDeque::new)
            .push_back((entry.start, entry.end));
    }

    ret
}
