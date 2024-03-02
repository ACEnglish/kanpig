use std::collections::VecDeque;
use std::io::BufRead;

use noodles_vcf::{self as vcf};

use crate::kanplug::{KDParams, Regions, KdpVcf};

/// Takes a vcf and filtering parameters to create in iterable which will
/// return chunks of variants in the same neighborhood
pub struct VcfChunker<R: BufRead> {
    pub m_vcf: vcf::reader::Reader<R>,
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
    chunk_count: u64,
    call_count: u64,
}

impl<R: BufRead> VcfChunker<R> {
    pub fn new(
        m_vcf: vcf::reader::Reader<R>,
        m_header: vcf::Header,
        regions: Regions,
        params: KDParams,
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
        }
    }

    /// Checks if entry passes all parameter conditions including
    /// within --bed regions, passing, and within expected size
    fn filter_entry(&mut self, entry: &vcf::Record) -> bool {
        // Is it inside a region
        let mut default = VecDeque::new();
        let m_coords = self
            .regions
            .get_mut(&entry.chromosome().to_string())
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

        if self.params.passonly & entry.is_filtered() {
            return false;
        }

        let size = entry.size();
        if self.params.sizemin > size || self.params.sizemax < size {
            return false;
        }

        // Need to make sure its sequence resolved
        // We'll make an exception for <DEL> in the future? (ref fetch, though..)
        true
    }

    /// Return the next vcf entry which passes parameter conditions
    fn get_next_entry(&mut self) -> Option<vcf::Record> {
        //let mut entry = vcf::Record::default();
        let mut entry = vcf::Record::default();

        loop {
            match self.m_vcf.read_record(&self.m_header, &mut entry) {
                Ok(0) => return None,
                Err(e) => {
                    error!("skipping invalid VCF entry {:?}", e);
                    continue;
                }
                Ok(_) if self.filter_entry(&entry) => return Some(entry),
                _ => continue,
            }
        }
    }

    /// Checks if this variant is within params.chunksize distance of last
    /// seen variant in this chunk
    fn entry_in_chunk(&mut self, entry: &vcf::Record) -> bool {
        let check_chrom = entry.chromosome().to_string();
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
                "{} chunks of {} variants",
                self.chunk_count, self.call_count
            );
            None
        }
    }
}
