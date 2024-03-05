/// Approaches for applying Haplotypes to a VarGraph
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::prelude::EdgeIndex;
use petgraph::visit::{Dfs, EdgeRef};
use std::collections::HashSet;
use std::iter::FromIterator;

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
    skip_edges: &Vec<EdgeIndex>,
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
        .edges(cur_node) // remove
        .filter(|edge| !skip_edges.contains(&edge.id()))
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
                skip_edges,
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

/// Return nodes which that have a 1-to-1 match to the haplotype
pub fn get_one_to_one(
    graph: &DiGraph<VarNode, ()>,
    target: &Haplotype,
    params: &KDParams,
) -> Vec<PathScore> {
    graph
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
        .collect()
}

/// Remove edges from graph that lead to paths which never pass through the kept nodes
pub fn prune_graph(
    graph: &DiGraph<VarNode, ()>,
    kept_paths: &[PathScore],
    source: &NodeIndex,
    sink: &NodeIndex,
) -> Vec<EdgeIndex> {
    let mut visited = HashSet::new();
    let mut dfs = Dfs::new(&graph, *sink);
    let kept_nodes: HashSet<NodeIndex> = HashSet::from_iter(
        kept_paths
            .iter()
            .flat_map(|i| i.path.clone())
            .collect::<Vec<NodeIndex>>(),
    );

    while let Some(node) = dfs.next(&graph) {
        visited.insert(node);
    }

    //edges_to_remove
    graph
        .edges_directed(*source, petgraph::Direction::Outgoing)
        .filter(|e| !kept_nodes.contains(&e.target()))
        .map(|e| e.id())
        .collect::<Vec<_>>()
}
