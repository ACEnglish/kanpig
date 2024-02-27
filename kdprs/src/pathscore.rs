use crate::cli::KDParams;
use crate::haplotype::Haplotype;
use crate::metrics::{self as metrics};
use crate::vargraph::VarNode;
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    path: Vec<NodeIndex>,
    sizesim: f32,
    cossim: f32,
}

impl Eq for PathScore {}

impl PartialEq for PathScore {
    fn eq(&self, other: &Self) -> bool {
        self.sizesim == other.sizesim && self.cossim == other.cossim
    }
}

impl Ord for PathScore {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.sizesim < other.sizesim {
            Ordering::Less
        } else if self.sizesim > other.sizesim {
            Ordering::Greater
        } else {
            self.cossim
                .partial_cmp(&other.cossim)
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
            cossim: 0.0,
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
                sizesim: 0.0,
                cossim: 0.0,
            };
        }

        let sizesim = metrics::sizesim(path_size.unsigned_abs(), target.size.unsigned_abs());
        // No need for cossim because sizesim is alredy a failure
        if sizesim < params.pctsize {
            return PathScore {
                path,
                sizesim,
                cossim: 0.0,
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

        // I might have messedup these signs
        let cossim = if std::cmp::max(target.size.unsigned_abs(), path_size.unsigned_abs())
            < params.wcoslen
        {
            metrics::weighted_cosinesim(&path_k, &target.kfeat)
        } else {
            metrics::cosinesim(&path_k, &target.kfeat)
        };

        PathScore {
            path,
            sizesim,
            cossim,
        }
    }
}
