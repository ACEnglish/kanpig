use crate::kplib::{
    metrics::overlaps, traverse::brute_force_find_path, traverse::get_one_to_one,
    vcftraits::KdpVcf, ChannelOutput, GenotypeAnno, GraphParams, Haplotype, PathScore, Ploidy,
};
use itertools::Itertools;
use noodles_vcf::variant::RecordBuf;
use petgraph::graph::{DiGraph, NodeIndex};

#[derive(Debug)]
pub struct VarNode {
    pub start: u64,
    pub end: u64,
    pub size: i64,
    pub entry: Option<RecordBuf>,
    pub kfeat: Vec<f32>,
}

impl VarNode {
    pub fn new(entry: RecordBuf, kmer: u8) -> Self {
        // Want to make a hash for these names for debugging, I think.
        let (start, end) = entry.boundaries();
        let (kfeat, size) = entry.to_kfeat(kmer);
        Self {
            start,
            end,
            size,
            entry: Some(entry),
            kfeat,
        }
    }

    pub fn new_anchor(kmer: u8) -> Self {
        Self {
            start: 0,
            end: 0,
            size: 0,
            entry: None,
            kfeat: vec![0f32; 4_usize.pow(kmer.into())],
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
    pub fn new(mut variants: Vec<RecordBuf>, kmer: u8) -> Self {
        if variants.is_empty() {
            panic!("Cannot create a graph from no variants");
        }

        let mut graph = DiGraph::new();

        let (chrom, start, end) = Variants::get_region(&variants);
        let mut node_indices = Vec::<NodeIndex<_>>::with_capacity(variants.len() + 2);
        node_indices.push(graph.add_node(VarNode::new_anchor(kmer)));

        node_indices.append(
            &mut variants
                .drain(..) // drain lets us move the entry without a clone
                .map(|entry| graph.add_node(VarNode::new(entry, kmer)))
                .collect(),
        );

        node_indices.push(graph.add_node(VarNode::new_anchor(kmer)));

        Self {
            chrom,
            start,
            end,
            node_indices,
            graph,
        }
    }

    /// Build the edges in the graph
    /// Only needs to be done fully if there are pileups to consider
    pub fn build(&mut self, full: bool) {
        if full {
            for pair in self.node_indices.iter().combinations(2) {
                if let [Some(up_node), Some(dn_node)] = [
                    self.graph.node_weight(*pair[0]),
                    self.graph.node_weight(*pair[1]),
                ] {
                    if !overlaps(up_node.start, up_node.end, dn_node.start, dn_node.end) {
                        self.graph.add_edge(*pair[0], *pair[1], ());
                    }
                }
            }
        } else {
            let src_node = NodeIndex::new(0);
            let snk_node = NodeIndex::new(self.graph.node_count() - 1);
            self.graph.add_edge(src_node, snk_node, ());
        }
    }

    /// Again, TR aware, we need to set the bounds for doing the pileup
    /// to the TR boundaries.
    fn get_region(entries: &[RecordBuf]) -> (String, u64, u64) {
        let chrom = entries[0].reference_sequence_name().to_string();

        let (min_start, max_end) = entries.iter().fold((u64::MAX, 0), |acc, e| {
            let (start, end) = e.boundaries();
            (acc.0.min(start), acc.1.max(end))
        });

        (chrom, min_start, max_end)
    }

    // Find the path through this graph that best fits
    // the haplotype push coverage onto the VarNodes
    pub fn apply_haplotype(&self, hap: &Haplotype, params: &GraphParams) -> PathScore {
        // if there are no variants in the hap, we don't want to apply the coverage.
        if params.one_to_one || (self.node_indices.len() - 2) > params.maxnodes {
            get_one_to_one(&self.graph, hap, params)
                .iter()
                .max()
                .cloned()
                .unwrap_or_else(PathScore::default)
        } else {
            brute_force_find_path(&self.graph, hap, params)
        }
    }
    /// Transform the graph back into annotated variants
    /// Note that this will take the entries out of the graph's VarNodes
    pub fn take_annotated(
        &mut self,
        paths: Vec<&[PathScore]>,
        coverages: Vec<u64>,
        ploidy: Vec<&Ploidy>,
    ) -> ChannelOutput {
        self.node_indices
            .iter_mut()
            .filter_map(|var_idx| {
                self.graph
                    .node_weight_mut(*var_idx)
                    .unwrap()
                    .entry
                    .take()
                    .map(|entry| {
                        let mut annos = Vec::with_capacity(coverages.len());
                        for (i, (&cov, &ploid)) in coverages.iter().zip(ploidy.iter()).enumerate() {
                            annos.push(GenotypeAnno::new(
                                var_idx, paths[i], cov, ploid, self.start, i,
                            ));
                        }
                        Some((entry, annos))
                    })
            })
            .collect::<ChannelOutput>()
    }
}
