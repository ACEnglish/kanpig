use crate::kplib::{GenotypeAnno, KDParams, KdpVcf, Ploidy, Regions};
use crossbeam_channel::Sender;
use noodles_vcf::{self as vcf};
use petgraph::graph::NodeIndex;
use std::collections::VecDeque;
use std::io::BufRead;

/// Takes a vcf and filtering parameters to create in iterable which will
/// return chunks of variants in the same neighborhood
pub struct VcfChunker<R: BufRead> {
    pub m_vcf: vcf::io::Reader<R>,
    pub m_header: vcf::Header,
    regions: Regions,
    params: KDParams,
    // Variables for tracking chunks
    cur_chrom: String,
    cur_end: u64,
    // When iterating, we will encounter a variant that no longer
    // fits in the current chunk. We need to hold on to it for the
    // next chunk
    hold_entry: Option<vcf::Record>,
    pub chunk_count: u64,
    pub call_count: u64,
    pub skip_count: u64,
    result_sender: Sender<Option<Vec<GenotypeAnno>>>,
}

impl<R: BufRead> VcfChunker<R> {
    pub fn new(
        m_vcf: vcf::io::Reader<R>,
        m_header: vcf::Header,
        regions: Regions,
        params: KDParams,
        result_sender: Sender<Option<Vec<GenotypeAnno>>>,
    ) -> Self {
        Self {
            m_vcf,
            m_header,
            regions,
            params,
            cur_chrom: String::new(),
            cur_end: 0,
            hold_entry: None,
            chunk_count: 0,
            call_count: 0,
            skip_count: 0,
            result_sender,
        }
    }

    /// Checks if entry passes all parameter conditions including
    /// within --bed regions, passing, and within expected size
    fn filter_entry(&mut self, entry: &vcf::Record) -> bool {
        if self.params.passonly & entry.is_filtered() {
            return false;
        }

        if !entry.valid_alt() {
            return false;
        }

        let size = entry.size();
        if self.params.sizemin > size || self.params.sizemax < size {
            return false;
        }

        // Is it inside a region
        let mut default = VecDeque::new();
        let m_coords = self
            .regions
            .get_mut(&entry.reference_sequence_name().to_string())
            .unwrap_or(&mut default);

        if m_coords.is_empty() {
            return false;
        }

        let (var_up, var_dn) = entry.boundaries();

        let (mut reg_up, mut reg_dn) = (0, 0);
        while let Some((coord_up, coord_dn)) = m_coords.front() {
            if var_up > *coord_dn {
                m_coords.pop_front();
                continue;
            }
            reg_up = *coord_up;
            reg_dn = *coord_dn;
            break;
        }

        if !(reg_up <= var_up && var_dn <= reg_dn) {
            return false;
        }

        // Need to make sure its sequence resolved
        // We'll make an exception for <DEL> in the future? (ref fetch, though..)
        true
    }

    /// Return the next vcf entry which passes parameter conditions
    fn get_next_entry(&mut self) -> Option<vcf::Record> {
        let mut entry = vcf::Record::default();

        loop {
            match self.m_vcf.read_record(&mut entry) {
                Ok(0) => return None,
                Err(e) => {
                    error!("skipping invalid VCF entry {:?}", e);
                    continue;
                }
                Ok(_) => {
                    if self.filter_entry(&entry) {
                        return Some(entry);
                    } else {
                        self.skip_count += 1;
                        let _ = self.result_sender.send(Some(vec![GenotypeAnno::new(
                            entry.clone(),
                            &NodeIndex::new(0),
                            &[],
                            0,
                            &Ploidy::Zero,
                        )]));
                    }
                }
            }
        }
    }

    /// Checks if this variant is within params.chunksize distance of last
    /// seen variant in this chunk
    /// If we wanted to be TR aware, when checking new_chunk, we don't just look at
    /// cur_end but also the TR catalog. We want to chunk all TR changes together
    /// regardless of their distance.
    fn entry_in_chunk(&mut self, entry: &vcf::Record) -> bool {
        let check_chrom = entry.reference_sequence_name().to_string();
        let new_chrom = !self.cur_chrom.is_empty() && check_chrom != self.cur_chrom;

        let (start, end) = entry.boundaries();
        let new_chunk = self.cur_end != 0 && self.cur_end + self.params.chunksize < start;

        self.cur_chrom = check_chrom;
        self.cur_end = if new_chrom {
            end
        } else {
            std::cmp::max(self.cur_end, end)
        };

        !(new_chrom || new_chunk)
    }
}

impl<R: BufRead> Iterator for VcfChunker<R> {
    type Item = Vec<vcf::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut ret = self.hold_entry.take().into_iter().collect::<Vec<_>>();

        while let Some(entry) = self.get_next_entry() {
            self.call_count += 1;
            if self.entry_in_chunk(&entry) {
                ret.push(entry);
            } else {
                self.hold_entry = Some(entry);
                self.chunk_count += 1;
                return Some(ret);
            }
        }

        if !ret.is_empty() {
            self.chunk_count += 1;
            Some(ret)
        } else {
            info!(
                "{} variants in {} chunks",
                self.call_count, self.chunk_count
            );
            info!("{} variants skipped", self.skip_count);
            None
        }
    }
}
