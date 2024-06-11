use crate::kplib::{metrics, Haplotype, KDParams, VarNode};
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    pub sizesim: f32,
    pub seqsim: f32,
    pub coverage: Option<u64>,
    pub path: Vec<NodeIndex>,
    pub full_target: bool, // is this path against the full target
    pub is_ref: bool,      // This path tried to use a reference allele
}

impl Eq for PathScore {}

impl PartialEq for PathScore {
    fn eq(&self, other: &Self) -> bool {
        self.full_target == other.full_target
            && self.sizesim == other.sizesim
            && self.seqsim == other.seqsim
    }
}

impl Ord for PathScore {
    // Sort by mean of size and sequence
    fn cmp(&self, other: &Self) -> Ordering {
        match self
            .full_target
            .partial_cmp(&other.full_target)
            .unwrap_or(Ordering::Equal)
        {
            Ordering::Equal => {
                let m_score = (self.sizesim + self.seqsim) / 2.0;
                let o_score = (other.sizesim + other.seqsim) / 2.0;
                m_score.partial_cmp(&o_score).unwrap_or(Ordering::Equal)
            }
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

        // Return the partials in order from all to least
        for (i, hap_parts) in targets.iter().enumerate() {
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
            //debug!("sqsim: {}", seqsim);
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
                full_target: i == 0,
                is_ref: false,
            };
        }
        PathScore {
            is_ref: false,
            ..Default::default()
        }
        //PathScore::default()
    }
}
