use crate::kplib::seq_to_kmer;
use noodles_vcf::{
    self as vcf,
    variant::record::info::field::value::Value,
    variant::record::info::field::key::Key,
};
use std::cmp::Ordering;
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
    fn to_kfeat(&self, kmer: u8, maxhom: usize) -> (Vec<f32>, i64);
    fn boundaries(&self) -> (u64, u64);
    fn size(&self) -> u64;
    fn is_filtered(&self) -> bool;
    fn variant_type(&self) -> Svtype;
    fn is_symbolic(&self) -> bool;
    fn is_bnd(&self) -> bool;
}

impl KdpVcf for vcf::Record {
    /// Convert variant sequence to Kfeat
    fn to_kfeat(&self, kmer: u8, maxhom: usize) -> (Vec<f32>, i64) {
        let ref_seq = self.reference_bases().to_string();
        let alt_seq = self
            .alternate_bases()
            .first()
            .expect("Can only work on sequence resolved variants")
            .to_string();

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
        let start: u64 = u64::try_from(usize::from(self.position())).unwrap() - 1;
        let end: u64 = u64::try_from(usize::from(self.end().expect("No Variant End"))).unwrap();
        (start, end)
    }

    /// grab entry's length from either SVLEN field or infer it from the REF ALT fields
    fn size(&self) -> u64 {
        let svlen = self
            .info()
            .get(&Key::from_str("SVLEN").unwrap_or_else(|_| panic!("No SVLEN INFO")));

        if let Some(Some(Value::Integer(svlen))) = svlen {
            return svlen.unsigned_abs() as u64;
        } else if let Some(Some(Value::Array(field::value::Array::Integer(svlen)))) = svlen {
            return svlen
                .first()
                .unwrap_or_else(|| panic!("Bad SVLEN"))
                .unwrap()
                .unsigned_abs() as u64;
        }

        let r_len: u64 = self.reference_bases().len() as u64;
        let a_len: u64 = if self.is_symbolic() {
                let (start, end) = self.boundaries();
                start.abs_diff(end) + 1
            } else {
                match self.alternate_bases().first() {
                    Some(alt) => alt.len(),
                    None => 0
                }
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
            .get(&Key::from_str("SVTYPE").expect("Unable to make key"))
        {
            // INFO/SVTYPE
            Some(Some(Value::String(svtype))) => svtype.parse().expect("Bad SVTYPE"),
            Some(Some(Value::Array(field::value::Array::String(svtype)))) => svtype
                .first()
                .cloned()
                .unwrap_or_else(|| panic!("Bad SVTYPE"))
                .expect("parsed")
                .parse()
                .unwrap(),
            // Direct from REF/ALT
            _ => match self.alternate_bases().first() {
                Some(allele::Allele::Bases(alt)) => {
                    match alt.len().cmp(&self.reference_bases().len()) {
                        Ordering::Greater => Svtype::Ins,
                        Ordering::Less => Svtype::Del,
                        Ordering::Equal if alt.len() == 1 => Svtype::Snp,
                        _ => Svtype::Unk,
                    }
                }
                Some(allele::Allele::Symbol(alt)) => Svtype::from_str(&alt.to_string())
                    .unwrap_or_else(|_| panic!("Bad Symbolic Alt")),
                _ => Svtype::Unk,
            },
        }
    }

    /// Checks if its a symbolic allele e.g. <DEL>
    /// Returns false if its a monozygotic reference
    fn is_symbolic(&self) -> bool {
        match self.alternate_bases().first() {
            Some(alt) => alt.contains('<'),
            None => false
        }
    }

    fn is_bnd(&self) -> bool {
        match self.alternate_bases().first() {
            Some(alt) => (alt.contains('[') || alt.contains(']')) && alt.contains(':'),
            None => false
        }
    }
}
