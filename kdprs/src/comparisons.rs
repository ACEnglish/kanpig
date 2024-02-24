use noodles_vcf::{
    self as vcf, record::alternate_bases::allele, record::info::field, record::Filters,
};

use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
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

pub fn entry_boundaries(entry: &vcf::Record, ins_inflate: bool) -> (u64, u64) {
    let mut start: u64 = u64::try_from(usize::from(entry.position())).unwrap() - 1;
    let mut end: u64 = u64::try_from(usize::from(entry.end().expect("No Variant End"))).unwrap();
    if ins_inflate & (entry_variant_type(entry) == Svtype::Ins) {
        let size = entry_size(entry);
        start -= size / 2;
        end += size / 2;
    }
    (start, end)
}

pub fn entry_size(entry: &vcf::Record) -> u64 {
    let svlen = entry
        .info()
        .get(&field::Key::from_str("SVLEN").expect("No SVLEN INFO"));
    if let Some(Some(field::Value::Integer(svlen))) = svlen {
        return svlen.unsigned_abs() as u64;
    } else if let Some(Some(field::Value::Array(field::value::Array::Integer(svlen)))) = svlen {
        return svlen[0].expect("Bad SVLEN").unsigned_abs() as u64;
    }

    let r_len: u64 = entry.reference_bases().len() as u64;
    let a_len: u64 = match entry.alternate_bases().first() {
        Some(allele::Allele::Bases(alt)) => alt.len() as u64,
        Some(allele::Allele::Symbol(_alt)) => {
            let (start, end) = entry_boundaries(entry, false);
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

    r_len.abs_diff(a_len) as u64
}

pub fn sizesim(size_a: usize, size_b: usize) -> f32 {
    if ((size_a == 0) || (size_b == 0)) && size_a == size_b {
        return 1.0;
    }
    std::cmp::max(std::cmp::min(size_a, size_b), 1) as f32
        / std::cmp::max(std::cmp::max(size_a, size_b), 1) as f32
}

pub fn overlaps(s1: u64, e1: u64, s2: u64, e2: u64) -> bool {
    std::cmp::max(s1, s2) < std::cmp::min(e1, e2)
}

pub fn entry_is_filtered(entry: &vcf::Record) -> bool {
    // Filter is None or PASS not in filter
    match entry.filters() {
        Some(map) => *map != Filters::Pass,
        None => false,
    }
}

pub fn entry_variant_type(entry: &vcf::Record) -> Svtype {
    match entry
        .info()
        .get(&field::Key::from_str("SVTYPE").expect("Unable to make key"))
    {
        // INFO/SVTYPE
        Some(Some(field::Value::String(svtype))) => svtype.parse().expect("Bad SVTYPE"),
        Some(Some(field::Value::Array(field::value::Array::String(svtype)))) => {
            svtype[0].clone().expect("Bad SVTYPE").parse().unwrap()
        }
        // Direct from REF/ALT
        _ => match entry.alternate_bases().first() {
            Some(allele::Allele::Bases(alt)) => {
                let asz = alt.len();
                let rsz = entry.reference_bases().len();
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

// keep? Maybe use this instead of cosine..?
/*
pub fn seqsim(seq_a: &String, seq_b: &String) -> f32 {
    let align_res = edlibAlignRs(
        seq_a.as_bytes(),
        seq_b.as_bytes(),
        &EdlibAlignConfigRs::default(),
    );
    let totlen: f32 = (seq_a.len() + seq_b.len()) as f32;
    (totlen - align_res.editDistance as f32) / totlen
}

// keep? Maybe use this instead of cosine.. ?
pub fn unroll_compare(a_seq: &String, b_seq: &String, p: usize, up: bool) -> f32 {
    let b_len = b_seq.len();
    let f = p % b_len;
    let position = b_len - f; // I'm worried about signs here
    if position >= b_len {
        return 0.0; // should never be called unless Symbolic alts are present, in which case we
                    // can't compare
    }
    // If up, a_seq is upstream of b_seq
    let rolled = match up {
        true => format!("{}{}", &b_seq[position..], &b_seq[..position]),
        false => format!("{}{}", &b_seq[..position], &b_seq[position..]),
    };
    seqsim(a_seq, &rolled)
}
*/
