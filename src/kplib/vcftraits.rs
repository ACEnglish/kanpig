use crate::kplib::seq_to_kmer;
use noodles_vcf::{
    variant::record::AlternateBases, variant::record::Filters, variant::RecordBuf, Header,
};
use std::str::FromStr;

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum Svtype {
    Ins,
    Del,
    Dup,
    Inv,
    Snp,
    Unk,
    //Repl should be one
}

impl FromStr for Svtype {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "DEL" => Ok(Svtype::Del),
            "INS" => Ok(Svtype::Ins),
            "DUP" => Ok(Svtype::Dup),
            "INV" => Ok(Svtype::Inv),
            "SNP" => Ok(Svtype::Snp),
            _ => Ok(Svtype::Unk),
        }
    }
}

pub trait KdpVcf {
    fn to_kfeat(&self, kmer: u8, maxhom: usize) -> (Vec<f32>, i64);
    fn boundaries(&self) -> (u64, u64);
    fn size(&self) -> u64;
    fn is_filtered(&self, header: &Header) -> bool;
    fn valid_alt(&self) -> bool;
    fn get_alt(&self) -> &str;
}

impl KdpVcf for RecordBuf {
    /// Convert variant sequence to Kfeat
    fn to_kfeat(&self, kmer: u8, maxhom: usize) -> (Vec<f32>, i64) {
        let ref_seq = self.reference_bases().to_string();
        let alt_seq = self.get_alt();

        let size = alt_seq.len() as i64 - ref_seq.len() as i64;

        let m_ref = seq_to_kmer(&ref_seq.as_bytes()[1..], kmer, false, maxhom);
        let m_alt = seq_to_kmer(&alt_seq.as_bytes()[1..], kmer, false, maxhom);

        let m_ret: Vec<_> = m_alt
            .iter()
            .zip(m_ref.iter())
            .map(|(&x, &y)| (x - y))
            .collect();

        (m_ret, size)
    }

    /// start and end positions of an entry
    fn boundaries(&self) -> (u64, u64) {
        let start: u64 = u64::try_from(usize::from(self.variant_start().unwrap())).unwrap() - 1;
        let end: u64 = start + self.reference_bases().len() as u64;
        (start, end)
    }

    /// grab entry's length from either SVLEN field or infer it from the REF ALT fields
    fn size(&self) -> u64 {
        let r_len: u64 = self.reference_bases().len() as u64;
        let a_len: u64 = self.get_alt().len() as u64;

        if r_len == a_len {
            if r_len == 1 {
                return 0;
            } else {
                return r_len;
            }
        }

        r_len.abs_diff(a_len)
    }

    /// checks if an entry's FILTER is '.' or PASS, true if it is filtered
    fn is_filtered(&self, header: &Header) -> bool {
        !(self.filters().is_empty()
            || self.filters().iter(header).any(|res| match res {
                Ok(s) => s == "PASS",
                Err(_) => false,
            }))
    }

    /// Alternate sequence isn't '.' or '*' or bnd or symbolic
    fn valid_alt(&self) -> bool {
        let alt = self.get_alt();
        alt != "." && alt != "*" && !alt.contains(':') && !alt.contains('<')
    }

    /// Returns the first alternate allele or a blank string with '.' if there isn't any
    fn get_alt(&self) -> &str {
        let alts = self.alternate_bases();
        match alts.len() {
            0 => ".",
            _ => alts.iter().next().expect("I just checked").unwrap(),
            //.to_string(), // I don't like all this String when str should be simplier
        }
    }
}
