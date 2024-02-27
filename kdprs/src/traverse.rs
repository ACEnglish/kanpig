use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use crate::cli::KDParams;
use crate::vargraph::VarNode;
use crate::haplotype::Haplotype;

struct PathScore {
    path: Vec<NodeIndex>,
    sizesim: f32,
    cossim: f32,
}

impl PathScore {
    pub fn new(graph: DiGraph<VarNode, ()>, path: Vec<NodeIndex>, target: Haplotype, params: KDParams) -> Self {
        // This is essentially find best path and find all paths where
        // I'm going to sum the path VarNodes and then do size/cos similarity
        // PathScores will be comparable
        let path_size = [graph.nodes[i].size for i in path].sum();
        if path_size.signum() != target.size.signum() {
            return PathScore { path, sizesim: 0.0, cossim: 0.0 };
        }
        
        let sizesim =  sizesimilarity(path_size.abs(), target.size.abs());
        // No need for cossim because sizesim is alredy a failure
        if sizesim < params.pctsize {
            return PathScore { path, sizesim, cossim: 0.0 };
        }

        let path_k = [graph.nodes[i].kfeat for i in path].sum();
        let cossim = if wcoslen < std::cmp::max(target.size.abs(), path_size.abs()) {
            weighted_cosinesime(path_k, target.kfeat)
        } else {
            cosinesim(path_k, target.kfeat)
        };

        PathScore { path, sizesim, cossim }
    }
}

/// Recursive depth first search of a VarGraph that's guided by size similarity
/// Search stops after maxpaths have been checked
/// Assumes NodeIndex 0 is src node and NodeIndex -1 is snk
/// Returns a PathScore just has the NodeIndex of nodes covred
pub fn find_path(
    graph: DiGraph<VarNode, ()>,
    target: Haplotype,
    mut maxpaths: u64,
    cur_len: i64,
    i_cur_node: Option<NodeIndex>,
    mut i_path: Option<Vec<NodeIndex>>,
    i_best_path: Option<PathScore>,
) -> Option<PathScore> {
    // First call, setup the snk
    let (mut path, mut best_path, mut cur_len) = match i_cur_node {
        Some(node) => {
            i_path.as_mut().expect("How?").push(node);
            (i_path, i_best_path, cur_len)
        }
        None => {
            (Some(vec![NodeIndex::new(0)]), Some(PathScore { path: vec![], sizesim:0.0, cossim:0.0 }), 0)
        }
    };

    let cur_len = cur_len + graph.node_weight(*path.as_ref().expect("How?").last().unwrap()).unwrap().size;

    // Order next nodes by how close they get us to the haplotype's length
    let mut diffs: Vec<_> = graph
        .edges(*path.as_ref().expect("How?").last().unwrap())
        .map(|edge| {
            let target_node = edge.target();
            let size = graph.node_weight(target_node).unwrap().size;
            (
                (target.size - (cur_len + size)).abs(),
                target_node,
            )
        })
        .collect();
    diffs.sort_by_key(|&(len_diff, _)| len_diff);

    for (_, next_node) in diffs {
        // I don't know if I need this, anymore
        if next_node.index() == graph.node_count() // stop case - snk node
        {
            // Let the PathScore have the path
            let n_best_path = Some(PathScore::new(graph, path.clone().unwrap(), target, params));
            best_path = if best_path > n_best_path {
                n_best_path
            } else {
                best_path
            };
            maxpaths -= 1;
        } else {
            let n_path = path.clone();
            // Best path if we go to the next node
            // I'm not getting the max path update, its going in, but not coming out
            let n_best_path = find_path(graph, target, maxpaths, cur_len, i_cur_node, n_path, best_path);

            best_path = if best_path > n_best_path {
                n_best_path
            } else {
                best_path
            };

            maxpaths -= 1;
            if maxpaths <= 0 {
                break;
            }
        }
        if maxpaths <= 0 {
            break;
        }
    }

    best_path
}
