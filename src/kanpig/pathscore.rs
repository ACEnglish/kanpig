use crate::kanpig::{Haplotype, metrics, VarNode, KDParams};
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    pub path: Vec<NodeIndex>,
    pub sizesim: f32,
    pub seqsim: f32,
    pub coverage: Option<u64>,
}

impl Eq for PathScore {}

impl PartialEq for PathScore {
    fn eq(&self, other: &Self) -> bool {
        self.sizesim == other.sizesim && self.seqsim == other.seqsim
    }
}

impl Ord for PathScore {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.sizesim < other.sizesim {
            Ordering::Less
        } else if self.sizesim > other.sizesim {
            Ordering::Greater
        } else {
            self.seqsim
                .partial_cmp(&other.seqsim)
                .unwrap_or(Ordering::Equal)
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
        // This is essentially find best path and find all paths where
        // I'm going to sum the path VarNodes and then do size/cos similarity
        // PathScores will be comparable
        let path_size: i64 = path
            .iter()
            .filter_map(|&node_index| graph.node_weight(node_index))
            .map(|x| x.size)
            .sum();

        if path_size.signum() != target.size.signum() {
            return PathScore {
                path,
                sizesim: -1.0, // ew
                seqsim: -1.0,
                coverage: None,
            };
        }

        let mut sizesim = metrics::sizesim(path_size.unsigned_abs(), target.size.unsigned_abs());
        // No need for seqsim because sizesim is alredy a failure
        if sizesim < params.sizesim {
            return PathScore {
                path,
                sizesim: -1.0,
                seqsim: -1.0,
                coverage: None,
            };
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

        let mut seqsim = metrics::seqsim(&path_k, &target.kfeat, params.minkfreq as f32);

        if seqsim < params.seqsim {
            seqsim = -1.0;
            sizesim = -1.0;
        }

        PathScore {
            path,
            sizesim,
            seqsim,
            coverage: None,
        }
    }
}
