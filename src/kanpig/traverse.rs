/// Approaches for applying Haplotypes to a VarGraph
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;

use crate::kanpig::{Haplotype, KDParams, PathScore, VarNode};

type PathResult = (Option<PathScore>, u64);

/// Recursive depth first search of a VarGraph that's guided by size similarity
/// Search stops after maxpaths have been checked
/// Assumes NodeIndex 0 is src node and NodeIndex -1 is snk
/// Returns a PathScore
pub fn brute_force_find_path(
    graph: &DiGraph<VarNode, ()>,
    target: &Haplotype,
    params: &KDParams,
    mut npaths: u64,
    cur_len: i64,
    i_cur_node: Option<NodeIndex>,
    mut i_path: Option<Vec<NodeIndex>>,
    i_best_path: Option<PathScore>,
) -> PathResult {
    let (path, mut best_path, cur_len, cur_node) = match i_cur_node {
        Some(node) => {
            i_path.as_mut().unwrap().push(node);
            (i_path.unwrap(), i_best_path.unwrap(), cur_len, node)
        }
        None => (vec![], PathScore::default(), 0, NodeIndex::new(0)),
    };

    let cur_len = cur_len + graph.node_weight(cur_node).unwrap().size;
    // Order next nodes by how close they get us to the haplotype's length
    let mut diffs: Vec<_> = graph
        .edges(cur_node)
        .map(|edge| {
            let target_node = edge.target();
            let size = graph.node_weight(target_node).unwrap().size;
            ((target.size.abs_diff(cur_len + size)), target_node)
        })
        .collect();
    diffs.sort_by_key(|&(len_diff, _)| len_diff);

    for (_, next_node) in diffs {
        if next_node.index() == graph.node_count() - 1 {
            let n_best_path = PathScore::new(graph, path.clone(), target, params);
            best_path = best_path.max(n_best_path);
            npaths += 1;
        } else {
            let (new_best, mp) = brute_force_find_path(
                graph,
                target,
                params,
                npaths,
                cur_len,
                Some(next_node),
                Some(path.clone()),
                Some(best_path.clone()),
            );

            if let Some(new_best) = new_best {
                best_path = best_path.max(new_best);
            }
            npaths = mp;
        }
        if npaths > params.maxpaths {
            break;
        }
    }

    (Some(best_path), npaths)
}

/// 1-to-1 : Before we go crazy searching for a path,
/// Lets just see if there is some single VarNode that
/// matches to the haplotype pretty well
pub fn one_to_one(
    graph: &DiGraph<VarNode, ()>,
    target: &Haplotype,
    params: &KDParams,
) -> Option<PathResult> {
    // we were asked not to do this
    if params.skip_one {
        return None;
    }
    let mut candidates: Vec<_> = graph
        .node_indices()
        .filter_map(|target_node| {
            let candidate = PathScore::new(graph, vec![target_node], target, params);
            if candidate.seqsim > 0.0 {
                Some(candidate)
            } else {
                None
            }
            //(node.size >= size_range_lower) & (node.size <= size_range_upper) & (
        })
        .collect();

    if candidates.is_empty() {
        None
    } else {
        candidates.sort();
        Some((candidates.iter().max().cloned(), 0))
    }
}
