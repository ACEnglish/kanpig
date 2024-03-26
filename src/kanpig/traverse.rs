/// Approaches for applying Haplotypes to a VarGraph
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::prelude::EdgeIndex;
use petgraph::visit::{Dfs, EdgeRef};
use std::collections::HashSet;
use std::iter::FromIterator;

use crate::kanpig::{Haplotype, KDParams, PathScore, VarNode};

#[derive(Clone)]
pub struct PathNodeState {
    pub size: i64,
    pub dist: u64,
    pub node: NodeIndex,
    pub path: Vec<NodeIndex>,
}

/// Search of a VarGraph that's guided by size similarity
/// Search stops after maxpaths have been checked
/// Assumes NodeIndex 0 is src node and NodeIndex -1 is snk
/// Returns a PathScore
pub fn brute_force_find_path(
    graph: &DiGraph<VarNode, ()>,
    target: &Haplotype,
    params: &KDParams,
    skip_edges: &[EdgeIndex],
) -> PathScore {
    let start_path = PathNodeState {
        size: 0,
        dist: target.size.unsigned_abs(), // this is for sorting
        node: NodeIndex::new(0),
        path: vec![],
    };
    let mut stack: Vec<PathNodeState> = vec![start_path];
    let mut best_path = PathScore::default();
    let mut npaths = 0;
    let snk_node = NodeIndex::new(graph.node_count() - 1);

    while let Some(cur_path) = stack.pop() {
        let mut any_push = false;

        // Throw all of cur_node's neighbors on the stack
        // Except snk_node, which is an indicator that the
        // current path has ended
        for next_node in graph.edges(cur_path.node).filter_map(|edge| {
            if skip_edges.contains(&edge.id()) {
                None
            } else {
                Some(edge.target())
            }
        }) {
            if next_node == snk_node {
                best_path =
                    best_path.max(PathScore::new(graph, cur_path.path.clone(), target, params));
                debug!("best path {:?}", best_path);
                npaths += 1;
            } else {
                any_push = true;
                let nsize = cur_path.size + graph.node_weight(next_node).unwrap().size;
                let mut npath = cur_path.path.clone();
                npath.push(next_node);
                stack.push(PathNodeState {
                    size: nsize,
                    dist: target.size.abs_diff(nsize),
                    node: next_node,
                    path: npath,
                });
            }
        }

        // Only need to sort when we've added to the stack
        if any_push {
            stack.sort_by_key(|node| std::cmp::Reverse(node.dist));
        }

        if npaths > params.maxpaths {
            break;
        }
    }

    best_path
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
