use noodles_vcf::header::record::value::{map::Contig, Map, map::contig::Name};
use indexmap::IndexMap;

pub type ContigMap = IndexMap<Name, Map<Contig>>;

pub fn build_region_tree(vcfA_contigs: &ContigMap, includebed: u8) {
    println!("{:?}", vcfA_contigs);
}
