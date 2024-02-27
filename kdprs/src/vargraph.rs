use crate::cli::KDParams;
use crate::haplotype::Haplotype;
use crate::metrics::overlaps;
use crate::pathscore::PathScore;
use crate::vcf_traits::KdpVcf;
use itertools::Itertools;
use noodles_vcf::{self as vcf};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;


/// Every --input variant is placed inside a node is turned into a graph.
#[derive(Debug)]
pub struct VarNode {
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub size: i64,
    pub kfeat: Vec<f32>,
    pub entry: Option<vcf::Record>,
    pub coverage: (Option<u64>, Option<u64>),
    pub cossim: (Option<f32>, Option<f32>),
    pub sizesim: (Option<f32>, Option<f32>),
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
            kfeat,
            entry: Some(entry),
            coverage: (None, None),
            cossim: (None, None),
            sizesim: (None, None),
        }
    }

    /// For the 'src' and 'snk' nodes, just need the name
    pub fn new_anchor(name: &str, kmer: u8) -> Self {
        Self {
            name: name.to_string(),
            start: 0,
            end: 0,
            size: 0,
            kfeat: vec![0f32; 4_usize.pow(kmer.into())],
            entry: None,
            coverage: (None, None),
            cossim: (None, None),
            sizesim: (None, None),
        }
    }
}

#[derive(Debug)]
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

        let mut node_indices = Vec::<NodeIndex<_>>::with_capacity(variants.len() + 2);
        node_indices.push(graph.add_node(VarNode::new_anchor("src", kmer)));

        node_indices.append(
            &mut variants
                .drain(..) // drain lets us move the entry without a clone
                .map(|entry| graph.add_node(VarNode::new(entry, kmer)))
                .collect(),
        );

        node_indices.push(graph.add_node(VarNode::new_anchor("snk", kmer)));

        for pair in node_indices.iter().combinations(2) {
            if let [Some(up_node), Some(dn_node)] =
                [graph.node_weight(*pair[0]), graph.node_weight(*pair[1])]
            {
                if !overlaps(up_node.start, up_node.end, dn_node.start, dn_node.end) {
                    graph.add_edge(*pair[0], *pair[1], ());
                }
            }
        }

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

    // Find the path through this graph that best fits
    // the haplotype push coverage onto the VarNodes
    pub fn apply_coverage(&self, hap: &Haplotype, params: &KDParams) -> Option<PathScore> {
        find_path(&self.graph, hap, params, 0, 0, None, None, None).0
    }
}

/// Recursive depth first search of a VarGraph that's guided by size similarity
/// Search stops after maxpaths have been checked
/// Assumes NodeIndex 0 is src node and NodeIndex -1 is snk
/// Returns a PathScore just has the NodeIndex of nodes covred

fn find_path(
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
