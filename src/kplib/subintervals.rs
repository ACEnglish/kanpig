
///
/// Leverage the CoverageTrack to find sub-intervals of pileups that
/// are new neighdist chunks
///
pub fn find_subintervals(reads: Vec<ReadData>, neighdist: u64) -> Vec<(u64, u64)> {
    let coords = vec![];
    for read in reads {
        for pileup in read.pileups {
            coords.push((pileup.position - neighdist, pileup.position + neighdist));
        }
    }
    let mut cov_track = CoverageTrack(Some(coords), 0);
    cov_track.intervals.merge_overlaps(); 
    let ret = vec![];
    for intv in cov_track.intervals {
        ret.push((intv.begin, intv.end));
    }

    ret
}
