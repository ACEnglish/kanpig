use crate::kanpig::{metrics, Haplotype, KDParams, VarNode};
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    pub sizesim: f32,
    pub seqsim: f32,
    pub coverage: Option<u64>,
    pub path: Vec<NodeIndex>,
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
        }
    }
}

impl PathScore {
    pub fn new(
        graph: &DiGraph<VarNode, ()>,
        path: Vec<NodeIndex>,
        target: &Haplotype,
        params: &KDParams,
    ) -> Self {
        let path_size: i64 = path
            .iter()
            .filter_map(|&node_index| graph.node_weight(node_index))
            .map(|x| x.size)
            .sum();

        if path_size.signum() != target.size.signum() {
            return PathScore::default();
        }

        let sizesim = metrics::sizesim(path_size.unsigned_abs(), target.size.unsigned_abs());
        debug!("szsim: {}", sizesim);
        if sizesim < params.sizesim {
            return PathScore::default();
        }

        let path_k: Vec<f32> = path
            .iter()
            .filter_map(|&node_index| graph.node_weight(node_index))
            .map(|x| x.kfeat.as_ref())
            .fold(
                vec![0f32; 4_usize.pow(params.kmer.into())],
                |acc: Vec<f32>, other: &Vec<f32>| {
                    acc.iter()
                        .zip(other)
                        .map(|(x, y)| x + y) // This line sums corresponding elements
                        .collect()
                },
            );

        let seqsim = metrics::seqsim(&path_k, &target.kfeat, params.minkfreq as f32);
        debug!("sqsim: {}", seqsim);
        if seqsim < params.seqsim {
            return PathScore::default();
        }

        PathScore {
            path,
            sizesim,
            seqsim,
            coverage: None,
        }
    }
}
