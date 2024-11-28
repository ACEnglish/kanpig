use crate::kplib::BedParser;
use rust_lapper::{Interval, Lapper};
use std::{collections::HashMap, str::FromStr};

#[derive(PartialEq, Debug)]
pub enum Ploidy {
    Zero,
    Haploid,
    Diploid,
    Polyploid,
    Unset,
}

impl FromStr for Ploidy {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim() {
            "0" => Ok(Ploidy::Zero),
            "1" => Ok(Ploidy::Haploid),
            _ => Ok(Ploidy::Unset),
        }
    }
}

impl Ploidy {
    pub fn from_value(value: u64) -> Self {
        match value {
            0 => Ploidy::Zero,
            1 => Ploidy::Haploid,
            2 => Ploidy::Diploid,
            3 => Ploidy::Polyploid,
            _ => Ploidy::Unset,
        }
    }
}

pub type Iv = Interval<u64, u64>;
pub type IvLookup = HashMap<String, Lapper<u64, u64>>;

#[derive(Clone)]
pub struct PloidyRegions {
    intervals: IvLookup,
}

impl PloidyRegions {
    pub fn new(path: &Option<std::path::PathBuf>) -> Self {
        if path.is_none() {
            return PloidyRegions {
                intervals: IvLookup::new(),
            };
        }

        let mut m_parser = BedParser::new(&path.clone().unwrap());
        let mut hold_intv = HashMap::<String, Vec<(u64, u64, u64)>>::new();
        for entry in m_parser.parse().into_iter() {
            if entry.data.is_none() {
                error!("No ploidy (0 or 1) in entry {:?}", entry);
                std::process::exit(1);
            }
            let m_ploy: Ploidy = entry.data.unwrap()[0].parse().unwrap();
            hold_intv
                .entry(entry.chrom)
                .or_default()
                .push((entry.start, entry.end, m_ploy as u64));
        }
        // Convert the original HashMap into a new HashMap<String, Lapper>
        let intervals: IvLookup = hold_intv
            .into_iter()
            .map(|(key, values)| {
                let ivs: Vec<Iv> = values
                    .into_iter()
                    .map(|(a, b, c)| Iv {
                        start: a,
                        stop: b,
                        val: c,
                    })
                    .collect();
                (key, Lapper::new(ivs))
            })
            .collect();

        PloidyRegions { intervals }
    }

    pub fn get_ploidy(&self, chrom: &String, start: u64) -> Ploidy {
        if let Some(lapper) = self.intervals.get(chrom) {
            match lapper.find(start, start + 1).next() {
                Some(i) => Ploidy::from_value(i.val),
                None => Ploidy::Unset,
            }
        } else {
            Ploidy::Unset
        }
    }
}
