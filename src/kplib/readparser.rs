use crate::kplib::{PileupSet, ReadsMap};

pub trait ReadParser {
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (ReadsMap, PileupSet, u64);
}
