use indexmap::IndexMap;
use noodles_vcf::header::record::value::{map::contig::Name, map::Contig, Map};

pub type ContigMap = IndexMap<Name, Map<Contig>>;
pub type Regions = IndexMap<String, Vec<(u64, u64)>>;
use crate::bedparser::BedParser;

pub fn build_region_tree(
    vcf_contigs: &ContigMap,
    includebed: Option<std::path::PathBuf>,
) -> Regions {
    let mut m_contigs = IndexMap::new();
    for (k, v) in vcf_contigs {
        let name = k.to_string();
        let length = match v.length() {
            Some(l) => l,
            None => {
                panic!("--vcf contig header missing length");
            }
        };
        m_contigs.insert(name, [(0, length as u64)].to_vec());
    }

    if includebed.is_none() {
        return m_contigs;
    }

    // Parse the Bed Lines and return them as the IndexMap.
    let mut ret = IndexMap::new();
    let mut m_parser = BedParser::new(&includebed.unwrap());
    for (chrom, m_start, m_stop) in m_parser.parse().into_iter() {
        if !m_contigs.contains_key(&chrom) {
            continue;
        }
        if !ret.contains_key(&chrom) {
            ret.insert(chrom.clone(), Vec::<(u64, u64)>::new());
        }

        if let Some(val) = ret.get_mut(&chrom) {
            val.push((m_start, m_stop))
        };
    }

    ret
}
