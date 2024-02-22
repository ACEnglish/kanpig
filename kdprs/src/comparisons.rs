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

pub fn coords_within(
    qstart: usize,
    qend: usize,
    rstart: usize,
    rend: usize,
    end_within: bool,
) -> bool {
    let ending = if end_within {
        qend <= rend
    } else {
        qend < rend
    };
    (qstart >= rstart) & ending
}

pub fn entry_boundaries(entry: &vcf::Record, ins_inflate: bool) -> (usize, usize) {
    let mut start = usize::from(entry.position()) - 1;
    let mut end = usize::from(entry.end().expect("No Variant End"));
    if ins_inflate & (entry_variant_type(entry) == Svtype::Ins) {
        let size = entry_size(entry);
        start -= size / 2;
        end += size / 2;
    }
    (start, end)
}

pub fn entry_size(entry: &vcf::Record) -> usize {
    let svlen = entry
        .info()
        .get(&field::Key::from_str("SVLEN").expect("No SVLEN INFO"));
    if let Some(Some(field::Value::Integer(svlen))) = svlen {
        return svlen.unsigned_abs() as usize;
    } else if let Some(Some(field::Value::Array(field::value::Array::Integer(svlen)))) = svlen {
        return svlen[0].expect("Bad SVLEN").unsigned_abs() as usize;
    }

    let r_len = entry.reference_bases().len();
    let a_len = match entry.alternate_bases().first() {
        Some(allele::Allele::Bases(alt)) => alt.len(),
        Some(allele::Allele::Symbol(_alt)) => {
            let (start, end) = entry_boundaries(entry, false);
            start.abs_diff(end) + 1 // I don't understand why I have to add 1 to match the
                                    // python code.
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


pub fn sizesim(size_a: usize, size_b: usize) -> (f32, isize) {
    if ((size_a == 0) || (size_b == 0)) && size_a == size_b {
        return (1.0, 0);
    }
    let pct = std::cmp::max(std::cmp::min(size_a, size_b), 1) as f32
        / std::cmp::max(std::cmp::max(size_a, size_b), 1) as f32;
    let diff = size_a as isize - size_b as isize;
    (pct, diff)
}


pub fn entry_within(entry: &vcf::Record, rstart: usize, rend: usize) -> bool {
    let (qstart, qend) = entry_boundaries(entry, false);
    let end_within = entry_variant_type(entry) != Svtype::Ins;
    coords_within(qstart, qend, rstart, rend, end_within)
}


pub fn overlaps(s1: usize, e1: usize, s2: usize, e2: usize) -> bool {
    std::cmp::max(s1, s2) < std::cmp::min(e1, e2)
}

pub fn entry_is_filtered(entry: &vcf::Record) -> bool {
    // Filter is None or PASS not in filter
    match entry.filters() {
        Some(map) => *map != Filters::Pass,
        None => false,
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
