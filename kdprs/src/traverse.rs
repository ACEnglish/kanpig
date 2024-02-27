use crate::cli::KDParams;
use crate::haplotype::Haplotype;
use crate::metrics::{self as metrics};
use crate::vargraph::VarNode;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct PathScore {
    path: Vec<NodeIndex>,
    sizesim: f32,
    cossim: f32,
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
        let cossim = if std::cmp::max(target.size.unsigned_abs(), path_size.unsigned_abs()) < params.wcoslen
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

/// Recursive depth first search of a VarGraph that's guided by size similarity
/// Search stops after maxpaths have been checked
/// Assumes NodeIndex 0 is src node and NodeIndex -1 is snk
/// Returns a PathScore just has the NodeIndex of nodes covred
pub fn find_path(
    graph: &DiGraph<VarNode, ()>,
    target: &Haplotype,
    params: &KDParams,
    mut npaths: u64,
    cur_len: i64,
    i_cur_node: Option<NodeIndex>,
    mut i_path: Option<Vec<NodeIndex>>,
    i_best_path: Option<PathScore>,
) -> (Option<PathScore>, u64) {
    // First call, setup the snk
    // Either way, unwrap the i_*
    let (path, mut best_path, cur_len, cur_node) = match i_cur_node {
        Some(node) => {
            i_path.as_mut().expect("How?").push(node);
            (i_path.unwrap(), i_best_path.unwrap(), cur_len, node)
        }
        None => (
            vec![],
            PathScore {
                path: vec![],
                sizesim: 0.0,
                cossim: 0.0,
            },
            0,
            NodeIndex::new(0),
        ),
    };

    let cur_len = cur_len + graph.node_weight(cur_node).unwrap().size;

    // Order next nodes by how close they get us to the haplotype's length
    let mut diffs: Vec<_> = graph
        .edges(cur_node)
        .map(|edge| {
            let target_node = edge.target();
            let size = graph.node_weight(target_node).unwrap().size;
            ((target.size - (cur_len + size)).abs(), target_node)
        })
        .collect();
    diffs.sort_by_key(|&(len_diff, _)| len_diff);
    for (_, next_node) in diffs {
        // stop case - snk node
        if next_node.index() == graph.node_count() - 1 {
            // Let the PathScore have the path
            let n_best_path = PathScore::new(graph, path.clone(), target, params);
            best_path = if n_best_path > best_path {
                n_best_path
            } else {
                best_path
            };
            npaths += 1;
        } else {
            let n_path = path.clone();
            // Best path if we go to the next node
            // I'm not getting the max path update, its going in, but not coming out
            (best_path, npaths) = match find_path(
                graph,
                target,
                params,
                npaths,
                cur_len,
                Some(next_node),
                Some(n_path),
                Some(best_path.clone()),
            ) {
                (Some(new_best), mp) => {
                    if new_best > best_path {
                        (new_best, mp)
                    } else {
                        (best_path, mp)
                    }
                }
                (None, mp) => (best_path, mp),
            };
        }
        if npaths > params.maxpaths {
            break;
        }
    }

    (Some(best_path), npaths)
}
