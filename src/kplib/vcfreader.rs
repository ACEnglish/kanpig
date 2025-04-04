use crate::kplib::{GenotypeAnno, KDParams, KdpVcf, Ploidy, Regions};
use crossbeam_channel::Sender;
use noodles_vcf::{self as vcf, variant::RecordBuf};
use petgraph::graph::NodeIndex;
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
    hold_entry: Option<RecordBuf>,
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
    fn filter_entry(&mut self, entry: &RecordBuf) -> bool {
        if self.params.passonly & entry.is_filtered(&self.m_header) {
            return false;
        }

        if !entry.valid_alt() {
            return false;
        }

        let size = entry.size() as u32;
        if self.params.sizemin > size || self.params.sizemax < size {
            return false;
        }

        // Is it inside a region
        let Some(m_coords) = self.regions.get_mut(entry.reference_sequence_name()) else {
            return false;
        };

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
    fn get_next_entry(&mut self) -> Option<RecordBuf> {
        let mut entry = RecordBuf::default();

        loop {
            match self.m_vcf.read_record_buf(&self.m_header, &mut entry) {
                Ok(0) => return None,
                Err(e) => {
                    error!("skipping invalid VCF entry {:?}", e);
                    continue;
                }
                Ok(_) => {
                    if self.filter_entry(&entry) {
                        // Clear samples early
                        *entry.samples_mut() = vcf::variant::record_buf::Samples::default();
                        return Some(entry);
                    } else {
                        self.skip_count += 1;
                        let _ = self.result_sender.send(Some(vec![GenotypeAnno::new(
                            entry.clone(),
                            &NodeIndex::new(0),
                            &[],
                            0,
                            &Ploidy::Zero,
                            0,
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
    fn entry_in_chunk(&mut self, entry: &RecordBuf) -> bool {
        let check_chrom = entry.reference_sequence_name().to_string();
        let new_chrom = !self.cur_chrom.is_empty() && check_chrom != self.cur_chrom;

        let (start, end) = entry.boundaries();
        let new_chunk = self.cur_end != 0 && self.cur_end + self.params.neighdist < start;

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
    type Item = Vec<RecordBuf>;

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
