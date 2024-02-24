use crate::comparisons::{entry_boundaries, overlaps};
use crate::kmer::var_to_kfeat;
use itertools::Itertools;
use noodles_vcf::{self as vcf};
use petgraph::graph::{DiGraph, NodeIndex};

#[derive(Debug)]
pub struct VarNode {
    name: String,
    start: u64,
    end: u64,
    size: i64,
    kfeat: Option<Vec<f32>>,
    entry: Option<vcf::Record>,
}

impl VarNode {
    pub fn new(entry: vcf::Record, kmer: u8) -> Self {
        let name = "".to_string(); // Want to make a hash for these names, I think. For debugging.
        let (start, end) = entry_boundaries(&entry, false);
        let (kfeat, size) = var_to_kfeat(&entry, kmer);
        Self {
            name,
            start,
            end,
            size,
            kfeat: Some(kfeat),
            entry: Some(entry),
        }
    }

    pub fn new_anchor(name: &str) -> Self {
        Self {
            name: name.to_string(),
            start: 0,
            end: 0,
            size: 0,
            kfeat: None,
            entry: None,
        }
    }
}

pub fn vars_to_graph(
    mut variants: Vec<vcf::Record>,
    kmer: u8,
) -> (Vec<NodeIndex>, DiGraph<VarNode, ()>) {
    let mut graph = DiGraph::new();

    let node_indices: Vec<NodeIndex<_>> = variants
        .drain(..) // drain lets us move the entry without a clone
        .map(|entry| graph.add_node(VarNode::new(entry, kmer)))
        .collect();

    for pair in node_indices.iter().combinations(2) {
        if let [Some(up_node), Some(dn_node)] = [graph.node_weight(*pair[0]), graph.node_weight(*pair[1])] {
            if !overlaps(up_node.start, up_node.end, dn_node.start, dn_node.end) {
                graph.add_edge(*pair[0], *pair[1], ());
            }
        } else {
            panic!("Fetched a node that doesn't exist");
        }
    }

    // Do I want to have the src/sink placed somewhere special?
    // Like, can I node_indices[-2], [-1] if I need it?
    let src_node = graph.add_node(VarNode::new_anchor("src"));
    let snk_node = graph.add_node(VarNode::new_anchor("snk"));

    // Add edges between 'src' node and all other nodes
    for node_index in &node_indices {
        graph.add_edge(src_node, *node_index, ());
    }

    // Add edges between all nodes and 'snk' node
    for node_index in &node_indices {
        graph.add_edge(*node_index, snk_node, ());
    }

    (node_indices, graph)
}
