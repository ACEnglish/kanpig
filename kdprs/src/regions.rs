use noodles_vcf::header::record::value::{map::Contig, Map, map::contig::Name};
use indexmap::IndexMap;

pub type ContigMap = IndexMap<Name, Map<Contig>>;
// Now I need a bed parser


pub fn build_region_tree(vcfA_contigs: &ContigMap, includebed: u8) -> IndexMap<String, Vec<(u64, u64)>> {

    let mut ret = IndexMap::new(); // chrom : [[start, end]]
    for (k, v) in vcfA_contigs {
        let name = k.to_string();
        let length = match v.length() {
            Some(l) => l,
            None => {
                panic!("--vcf contig header missing length");
                std::process::exit(1);
            }
        };
        ret.insert(name, [(0, length as u64)].to_vec());
    }

    ret
}
