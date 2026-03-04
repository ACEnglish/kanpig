use crate::kplib::{metrics, vargraph::VarNode, GraphParams, Haplotype, KmerVec, SequenceMeta};
use petgraph::graph::{DiGraph, NodeIndex};
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct PathScore {
    pub score: f32, // Score(P) = S((SS + SZ) / 2) − λ ⋅ ∣L(P)−E∣
    #[allow(dead_code)]
    pub sizesim: f32,
    pub seqsim: f32,
    pub path: Vec<NodeIndex>,
    pub full_target: bool, // Does this path use partial
    pub meta: SequenceMeta,
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
            full_target: false,
            meta: SequenceMeta::default(),
        }
    }
}

impl PathScore {
    pub fn new(
        graph: &DiGraph<VarNode, ()>,
        path: Vec<NodeIndex>,
        path_size: i64,
        targets: &[Haplotype],
        params: &GraphParams,
        target: &Haplotype,
    ) -> Self {
        let mut path_k: Option<KmerVec> = None;
        let mut best_path = PathScore {
            meta: target.meta.clone(),
            ..Default::default()
        };
        // Return the partials in order from all to least
        for hap_parts in targets {
            if path_size.signum() != hap_parts.size.signum() {
                continue;
            }

            let sizesim = metrics::sizesim(path_size.unsigned_abs(), hap_parts.size.unsigned_abs());

            if sizesim < params.sizesim {
                continue;
            }

            if path_k.is_none() {
                let mut kv = KmerVec::blank();
                for node in path.iter().filter_map(|&n| graph.node_weight(n)) {
                    kv += &node.kmers;
                }
                path_k = Some(kv);
            }

            let pk = path_k.as_ref().unwrap();

            let fine_sim = pk.fine_similarity(&hap_parts.kmers);
            if fine_sim < params.seqsim {
                continue;
            }

            let coarse_sim = pk.coarse_similarity(&hap_parts.kmers);
            let seqsim = (fine_sim * coarse_sim).sqrt(); // Geometric Mean

            let mut score =
                ((seqsim + sizesim) / 2.0) - (params.fpenalty * hap_parts.partial as f32);

            if params.squish {
                score -= params.gpenalty * (path.len() - 1) as f32;
            } else {
                score -= params.gpenalty * hap_parts.n.abs_diff(path.len() as u64) as f32;
            }

            if score > best_path.score {
                best_path = PathScore {
                    score,
                    path: path.clone(),
                    sizesim,
                    seqsim,
                    full_target: hap_parts.partial == 0,
                    meta: target.meta.clone(),
                };
            }
        }

        best_path
    }
}
