use crate::kanplug::seq_to_kmer;
use noodles_vcf::{
    self as vcf, record::alternate_bases::allele, record::info::field, record::Filters,
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

/// Convert vcf::Record to kfeat
pub trait KdpVcf {
    fn to_kfeat(&self, kmer: u8) -> (Vec<f32>, i64);
    fn boundaries(&self) -> (u64, u64);
    fn size(&self) -> u64;
    fn is_filtered(&self) -> bool;
    fn variant_type(&self) -> Svtype;
}

impl KdpVcf for vcf::Record {
    /// Convert variant sequence to Kfeat
    fn to_kfeat(&self, kmer: u8) -> (Vec<f32>, i64) {
        let ref_seq = self.reference_bases().to_string();
        let alt_seq = self
            .alternate_bases()
            .first()
            .expect("Can only work on sequence resolved variants")
            .to_string();

        let size = alt_seq.len() as i64 - ref_seq.len() as i64;

        let m_ref = seq_to_kmer(&ref_seq.as_bytes()[1..], kmer, false);
        let m_alt = seq_to_kmer(&alt_seq.as_bytes()[1..], kmer, false);

        let m_ret: Vec<_> = m_alt
            .iter()
            .zip(m_ref.iter())
            .map(|(&x, &y)| (x - y))
            .collect();

        (m_ret, size)
    }

    /// start and end positions of an entry
    fn boundaries(&self) -> (u64, u64) {
        let start: u64 = u64::try_from(usize::from(self.position())).unwrap() - 1;
        let end: u64 = u64::try_from(usize::from(self.end().expect("No Variant End"))).unwrap();
        (start, end)
    }

    /// grab entry's length from either SVLEN field or infer it from the REF ALT fields
    fn size(&self) -> u64 {
        let svlen = self
            .info()
            .get(&field::Key::from_str("SVLEN").expect("No SVLEN INFO"));
        if let Some(Some(field::Value::Integer(svlen))) = svlen {
            return svlen.unsigned_abs() as u64;
        } else if let Some(Some(field::Value::Array(field::value::Array::Integer(svlen)))) = svlen {
            return svlen[0].expect("Bad SVLEN").unsigned_abs() as u64;
        }

        let r_len: u64 = self.reference_bases().len() as u64;
        let a_len: u64 = match self.alternate_bases().first() {
            Some(allele::Allele::Bases(alt)) => alt.len() as u64,
            Some(allele::Allele::Symbol(_alt)) => {
                let (start, end) = self.boundaries();
                start.abs_diff(end) + 1
            }
            _ => 0,
        };

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
    fn is_filtered(&self) -> bool {
        match &self.filters() {
            Some(map) => **map != Filters::Pass,
            None => false,
        }
    }

    /// return the Svtype of a vcf entry
    fn variant_type(&self) -> Svtype {
        match self
            .info()
            .get(&field::Key::from_str("SVTYPE").expect("Unable to make key"))
        {
            // INFO/SVTYPE
            Some(Some(field::Value::String(svtype))) => svtype.parse().expect("Bad SVTYPE"),
            Some(Some(field::Value::Array(field::value::Array::String(svtype)))) => {
                svtype[0].clone().expect("Bad SVTYPE").parse().unwrap()
            }
            // Direct from REF/ALT
            _ => match self.alternate_bases().first() {
                Some(allele::Allele::Bases(alt)) => {
                    let asz = alt.len();
                    let rsz = self.reference_bases().len();
                    if asz > rsz {
                        Svtype::Ins
                    } else if asz < rsz {
                        Svtype::Del
                    } else if asz == 1 {
                        Svtype::Snp
                    } else {
                        Svtype::Unk
                    }
                }
                Some(allele::Allele::Symbol(alt)) => {
                    Svtype::from_str(&alt.to_string()).expect("Bad Symbolic Alt")
                }
                _ => Svtype::Unk,
            },
        }
    }
}
