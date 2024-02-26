use crate::metrics::overlaps;
use crate::vcf_traits::KdpVcf;
use itertools::Itertools;
use noodles_vcf::{self as vcf};
use petgraph::graph::{DiGraph, NodeIndex};

/// Every --input variant is placed inside a node is turned into a graph.
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
        // Want to make a hash for these names for debugging, I think.
        let name = "".to_string();
        let (start, end) = entry.boundaries();
        let (kfeat, size) = entry.to_kfeat(kmer);
        Self {
            name,
            start,
            end,
            size,
            kfeat: Some(kfeat),
            entry: Some(entry),
        }
    }

    /// For the 'src' and 'snk' nodes, just need the name
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

pub struct Variants {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub node_indices: Vec<NodeIndex>,
    pub graph: DiGraph<VarNode, ()>,
}

/// Build a graph of all variants in a chunk.
/// Assumes variants are ordered by position (small to large)
/// Variants will have edges to every downstream variant that it does not overlap
/// The graph has an upstream 'src' node that point to every variant node
/// The graph has a dnstream 'snk' node that is pointed to by every variant node and 'src'
impl Variants {
    pub fn new(mut variants: Vec<vcf::Record>, kmer: u8) -> Self {
        if variants.is_empty() {
            panic!("Cannot create a graph from no variants");
        }

        let mut graph = DiGraph::new();

        let (chrom, start, end) = Variants::get_region(&variants);

        let mut node_indices: Vec<NodeIndex<_>> = variants
            .drain(..) // drain lets us move the entry without a clone
            .map(|entry| graph.add_node(VarNode::new(entry, kmer)))
            .collect();

        for pair in node_indices.iter().combinations(2) {
            if let [Some(up_node), Some(dn_node)] =
                [graph.node_weight(*pair[0]), graph.node_weight(*pair[1])]
            {
                if !overlaps(up_node.start, up_node.end, dn_node.start, dn_node.end) {
                    graph.add_edge(*pair[0], *pair[1], ());
                }
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
        node_indices.push(src_node);

        // Add edges between all nodes and 'snk' node
        for node_index in &node_indices {
            graph.add_edge(*node_index, snk_node, ());
        }
        node_indices.push(snk_node);

        Variants {
            chrom,
            start,
            end,
            node_indices,
            graph,
        }
    }

    fn get_region(entries: &Vec<vcf::Record>) -> (String, u64, u64) {
        let chrom = entries[0].chromosome().to_string();
        let mut min_start = u64::MAX;
        let mut max_end = 0;

        for e in entries {
            let (start, end) = e.boundaries();
            if start < min_start {
                min_start = start;
            }
            if end > max_end {
                max_end = end;
            }
        }

        (chrom, min_start, max_end)
    }
}
