use crate::cli::KDParams;
use crate::haplotype::Haplotype;
use crate::pathscore::PathScore;
use crate::vargraph::VarNode;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;

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
            ((target.size - (cur_len + size)).abs(), target_node)
        })
        .collect();
    diffs.sort_by_key(|&(len_diff, _)| len_diff);

    for (_, next_node) in diffs {
        if next_node.index() == graph.node_count() - 1 {
            // Let the PathScore have the path
            let n_best_path = PathScore::new(graph, path.clone(), target, params);
            best_path = best_path.max(n_best_path);
            npaths += 1;
        } else {
            let n_path = path.clone();
            // Best path if we go to the next node
            // I'm not getting the max path update, its going in, but not coming out
            let (new_best, mp) = find_path(
                graph,
                target,
                params,
                npaths,
                cur_len,
                Some(next_node),
                Some(n_path),
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
