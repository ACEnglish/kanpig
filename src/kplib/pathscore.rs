use crate::kplib::{metrics, Haplotype, KDParams, VarNode};
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    pub score: f32, // Score(P) = S((SS + SZ) / 2) − λ ⋅ ∣L(P)−E∣
    pub sizesim: f32,
    pub seqsim: f32,
    pub coverage: Option<u64>,
    pub path: Vec<NodeIndex>,
    pub full_target: bool, // Does this path use partial
    pub is_ref: bool,      // This path tried to use a reference allele
}

impl Eq for PathScore {}

impl PartialEq for PathScore {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
    }
}

impl Ord for PathScore {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .partial_cmp(&other.score)
            .unwrap_or(Ordering::Equal)
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
            score: 0.0,
            path: vec![],
            sizesim: 0.0,
            seqsim: 0.0,
            coverage: None,
            full_target: false,
            is_ref: true,
        }
    }
}

impl PathScore {
    pub fn new(
        graph: &DiGraph<VarNode, ()>,
        path: Vec<NodeIndex>,
        path_size: i64,
        targets: &[Haplotype],
        params: &KDParams,
    ) -> Self {
        let mut path_k: Option<Vec<f32>> = None;
        let mut best_path = PathScore {
            is_ref: false,
            ..Default::default()
        };
        // Return the partials in order from all to least
        for hap_parts in targets {
            if path_size.signum() != hap_parts.size.signum() {
                continue;
            }

            let sizesim = metrics::sizesim(path_size.unsigned_abs(), hap_parts.size.unsigned_abs());
            //debug!("szsim: {}", sizesim);
            if sizesim < params.sizesim {
                continue;
            }

            if path_k.is_none() {
                // only make if it is ever needed
                path_k = Some(
                    path.iter()
                        .filter_map(|&node_index| graph.node_weight(node_index))
                        .map(|x| x.kfeat.as_ref())
                        .fold(
                            vec![0f32; 4_usize.pow(params.kmer.into())],
                            |acc: Vec<f32>, other: &Vec<f32>| {
                                acc.iter().zip(other).map(|(x, y)| x + y).collect()
                            },
                        ),
                );
            }

            let seqsim = metrics::seqsim(
                path_k.as_ref().unwrap(),
                &hap_parts.kfeat,
                params.minkfreq as f32,
            );

            if seqsim < params.seqsim {
                continue;
            }

            let score = ((seqsim + sizesim) / 2.0)
                - (params.gpenalty * hap_parts.parts.len().abs_diff(path.len()) as f32)
                - (params.fpenalty * hap_parts.partial as f32);

            if score > best_path.score {
                best_path = PathScore {
                    score,
                    path: path.clone(),
                    sizesim,
                    seqsim,
                    coverage: None,
                    full_target: hap_parts.partial == 0,
                    is_ref: false,
                };
            }
        }

        best_path
    }
}
