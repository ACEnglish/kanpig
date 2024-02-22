/*
 * This is going to hold methods that will
 * Convert a vcf entry into a Variant Node
 * Create a graph from a set of vcf entries
*/
pub struct VariantNode {
    pub kfeat: Vec<u32>,
    pub size: u64,
}

/// from_entry would be better
pub fn var_to_node(variants: vcf::Record, kmer: u8) -> VariantNode {

}
