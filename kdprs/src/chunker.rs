pub struct Filter {
    pub sizemin: u64,
    pub sizemax: u64,
}

/*
 * build_region_tree
 * merge_region_tree_overlaps
 * filter (region, size, pass, resolved) on each file
 *
 * file_zipper to put them together
 * chunker to make the units for parsing
 * so chunk/zip can be one.. except that sometimes its 1 vcf and sometimes 2. So we want to keep
 * them separate
 *
 * regions = build_region_tree
 * file1 = filter(vcf, filter settings, regions)
 * file2 = filter(vcf, filter settings, regions)
 * zipped = file_zipper([file1, file2])
 * chunks = chunker(zipped)
 */
