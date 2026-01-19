use rust_lapper::{Interval, Lapper};

#[derive(Clone)]
pub struct CoverageTrack {
    pub intervals: Lapper<u64, ()>,
    pub buffer: u64,
}

impl CoverageTrack {
    pub fn new(reads: Option<Vec<(u64, u64)>>, buffer: u64) -> Self {
        let intervals: Vec<Interval<u64, ()>> = reads
            .unwrap_or_default()
            .into_iter()
            .map(|(start, stop)| Interval {
                start,
                stop,
                val: (),
            })
            .collect();

        Self {
            intervals: Lapper::new(intervals),
            buffer,
        }
    }

    pub fn count_spanning_reads(&self, mut query_start: u64, mut query_end: u64) -> u64 {
        query_start -= self.buffer;
        query_end += self.buffer;
        self.intervals
            .find(query_start, query_end)
            .filter(|read| read.start <= query_start && query_end <= read.stop)
            .count() as u64
    }
}
