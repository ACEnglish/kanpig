use crate::kplib::{pileup::ReadPileup, CoverageTrack};

///
/// Leverage the CoverageTrack to find sub-intervals of pileups that
/// are new neighdist chunks
///
pub fn find_subintervals(reads: &Vec<ReadPileup>, neighdist: u64) -> Vec<(u64, u64)> {
    let mut coords = vec![];
    for read in reads {
        for pileup in &read.pileups {
            coords.push((pileup.position - neighdist, pileup.end + neighdist));
        }
    }
    let mut cov_track = CoverageTrack::new(Some(coords), 0);
    cov_track.intervals.merge_overlaps();

    cov_track
        .intervals
        .into_iter()
        .map(|i| (i.start, i.stop))
        .collect()
}
