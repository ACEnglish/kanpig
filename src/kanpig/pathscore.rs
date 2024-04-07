use crate::kanpig::{metrics, Haplotype, KDParams, VarNode};
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    pub sizesim: f32,
    pub seqsim: f32,
    pub coverage: Option<u64>,
    pub path: Vec<NodeIndex>,
    pub align_pct: f32, // percent of the haplotype used 
}

impl Eq for PathScore {}

impl PartialEq for PathScore {
    fn eq(&self, other: &Self) -> bool {
        self.sizesim == other.sizesim && self.seqsim == other.seqsim
    }
}

impl Ord for PathScore {
    fn cmp(&self, other: &Self) -> Ordering {
        match self
            .sizesim
            .partial_cmp(&other.sizesim)
            .unwrap_or(Ordering::Equal)
        {
            Ordering::Equal => self
                .seqsim
                .partial_cmp(&other.seqsim)
                .unwrap_or(Ordering::Equal),
            other_ordering => other_ordering,
        }
    }
}

impl PartialOrd for PathScore {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Default for PathScore {
    fn default() -> PathScore {
        PathScore {
            path: vec![],
            sizesim: 0.0,
            seqsim: 0.0,
            coverage: None,
            align_pct: 0.0,
        }
    }
}

impl PathScore {
    pub fn new(
        graph: &DiGraph<VarNode, ()>,
        path: Vec<NodeIndex>,
        targets: &Vec<Haplotype>,
        target_size: i64,
        params: &KDParams,
    ) -> Self {
        let path_size: i64 = path
            .iter()
            .filter_map(|&node_index| graph.node_weight(node_index))
            .map(|x| x.size)
            .sum();

        let path_k: Vec<f32> = path
            .iter()
            .filter_map(|&node_index| graph.node_weight(node_index))
            .map(|x| x.kfeat.as_ref())
            .fold(
                vec![0f32; 4_usize.pow(params.kmer.into())],
                |acc: Vec<f32>, other: &Vec<f32>| {
                    acc.iter().zip(other).map(|(x, y)| x + y).collect()
                },
            );

        // Return the partials in order from all to least
        for hap_parts in targets {
            if path_size.signum() != hap_parts.size.signum() {
                continue;
            }

            let sizesim = metrics::sizesim(path_size.unsigned_abs(), hap_parts.size.unsigned_abs());
            debug!("szsim: {}", sizesim);
            if sizesim < params.sizesim {
                continue;
            }

            let seqsim = metrics::seqsim(&path_k, &hap_parts.kfeat, params.minkfreq as f32);
            debug!("sqsim: {}", seqsim);
            if seqsim < params.seqsim {
                continue;
            }
            // Stop on the first one that matches
            // This has weird tying implications, I think
            return PathScore {
                path,
                sizesim,
                seqsim,
                coverage: None,
                align_pct: (hap_parts.size.unsigned_abs() as f32 / target_size.unsigned_abs() as f32).abs()
            };
        }
        PathScore::default()
    }
}
