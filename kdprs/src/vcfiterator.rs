use std::io::BufRead;
use std::collections::VecDeque;

use noodles_vcf::{self as vcf};

use crate::cli::KDParams;
use crate::comparisons;
use crate::regions::Regions;


pub struct VCFIter<R: BufRead> {
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

}

impl<R: BufRead> VCFIter<R> {
    pub fn new(
        m_vcf: vcf::reader::Reader<R>,
        m_header: vcf::Header,
        regions: Regions,
        params: KDParams,
    ) -> Self {
        Self {
            m_vcf:m_vcf,
            m_header:m_header,
            regions:regions,
            params:params,
            cur_chrom:String::new(),
            cur_end:0,
            hold_entry:None,
        }
    }

    fn filter_entry(&mut self, entry: &vcf::Record) -> bool {
        if self.params.passonly & comparisons::entry_is_filtered(entry) {
            return false;
        }

        let size = comparisons::entry_size(entry);
        if self.params.sizemin > size || self.params.sizemax < size {
            return false;
        }

        let chrom = entry.chromosome().to_string();
        let m_coords: &mut VecDeque<(u64, u64)> = match self.regions.get_mut(&chrom) {
            Some(i) => i,
            None => return false,
        };

        if m_coords.is_empty() {
            return false;
        }

        let (var_up, var_dn) = comparisons::entry_boundaries(entry, false);

        let mut reg_up = 0;
        let mut reg_dn = 0;
        while let Some((coord_up, coord_dn)) = m_coords.front() {
            if var_up > *coord_dn {
                m_coords.pop_front();
                continue;
            }
            reg_up = *coord_up;
            reg_dn = *coord_dn;
            break;
        }

        if var_up < reg_up {
            return false;
        }

        if !(reg_up <= var_up && var_dn <= reg_dn) {
            return false;
        }

        true
    }

    fn get_next_entry(&mut self) -> Option<vcf::Record> {
        let mut entry = vcf::Record::default();

        loop {
            match self.m_vcf.read_record(&self.m_header, &mut entry) {
                Ok(0) => return None,
                Err(e) => {
                    error!("skipping invalid VCF entry {:?}", e);
                    continue
                },
                Ok(_) if self.filter_entry(&entry) => return Some(entry),
                _ => continue,
            }
        }
    }
    
    fn entry_in_chunk(&mut self, entry: &vcf::Record) -> bool {
        let check_chrom = entry.chromosome().to_string();
        let new_chrom = self.cur_chrom.len() != 0 && check_chrom != self.cur_chrom;
        let (start, end) = comparisons::entry_boundaries(entry, false);
        let new_chunk = self.cur_end != 0 && self.cur_end + self.params.chunksize < start;

        // Update where we're looking
        self.cur_chrom = check_chrom;
        self.cur_end = std::cmp::max(self.cur_end, end);

        !(new_chunk || new_chrom)
    }
}

impl<R: BufRead> Iterator for VCFIter<R> {
    type Item = Vec<vcf::Record>;
    
    fn next(&mut self) -> Option<Self::Item> {
        let mut ret = match &self.hold_entry {
            Some(entry) => {
                vec![entry.clone()]
            },
            None => vec![],
        };
        self.hold_entry = None;

        loop {
            match self.get_next_entry() {
                Some(entry) => {
                    if self.entry_in_chunk(&entry) {
                        ret.push(entry);
                    } else {
                        self.hold_entry = Some(entry);
                        return Some(ret);
                    }
                },
                None => break,
            };
        }

        match ret.is_empty() {
            false => Some(ret),
            true => None,
        }
    }
}
