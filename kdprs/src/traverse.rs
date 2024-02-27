use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use crate::vargraph::VarNode;
use crate::haplotype::Haplotype;

/// Recursive depth first search of a VarGraph that's guided by size similarity
/// Search stops after maxpaths have been checked
/// Assumes NodeIndex 0 is src node and NodeIndex -1 is snk
pub fn find_path(
    graph: DiGraph<VarNode, ()>,
    target: Haplotype,
    mut maxpaths: u64,
    cur_len: i64,
    i_cur_node: Option<NodeIndex>,
    mut i_path: Option<Vec<NodeIndex>>,
    i_best_path: Option<Vec<NodeIndex>>,
) -> Vec<NodeIndex> {
    // First call, setup the snk
    let (mut path, mut best_path, mut cur_len) = match i_cur_node {
        Some(node) => {
            i_path.as_mut().expect("How?").push(node);
            (i_path, i_best_path, cur_len)
        }
        None => {
            (Some(vec![NodeIndex::new(0)]), Some(vec![]), 0)
        }
    };

    let cur_len = cur_len + graph.node_weight(*path.as_ref().expect("How?").last().unwrap()).unwrap().size;

    // I'm questioning this sort
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
        if next_node.index() == graph.node_count()
            && *path.as_ref().expect("How?").last().unwrap() == NodeIndex::new(0)
        {
            best_path = None; //get_better_path(path, haplotype, params);
            maxpaths -= 1;
        } else {
            let n_path = path.clone();
            for sub_path in Vec::<usize>::new() {
                //find_path(awe....) {
                best_path = None; //get_better_path(sub_path, best_path, params);
                maxpaths -= 1;
                if maxpaths <= 0 {
                    break;
                }
            }
        }
        if maxpaths <= 0 {
            break;
        }
    }

    best_path.unwrap()
}
// Sum the VarNodes in path by size
// ensure is same sign and  over the kdparams, otherwise we aren't including it
// Sum the VarNodes by
// check if this path is better
// if path > best_path:
//  path = best_path
// most -= 1;
