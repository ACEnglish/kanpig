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
}

impl<R: BufRead> VCFIter<R> {
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
}

impl<R: BufRead> Iterator for VCFIter<R> {
    type Item = vcf::Record;

    fn next(&mut self) -> Option<Self::Item> {
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
}
