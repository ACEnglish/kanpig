"""
Approximate phasing of a comparison VCF to match a phased baseline VCF using kmer featurization distances

Todo: 
Validate Inputs e.g. (need tbi)
Instead of max_nodes, I can probably do a max-paths
I think an A* could be faster if I can use the sizediff of current path from hap-path as the heuristic
Figure out if I can do something about the size-separator. Like, it worked good for 'isn't random', but
    What I really need is "are the same"
"""
import sys
import logging
import truvari

import kdp

def pull_variants(graph, used, h1_min_path, h1, h2_min_path, h2, chunk_id, sample=0):
    """
    Update the variants in the two paths and return them
    """
    ret_entries = []
    for node in used:
        g1 = int(node in h1_min_path.path)
        g2 = int(node in h2_min_path.path)
        v = graph.nodes[node]["variant"]
        v.samples[sample]["GT"] = (g1, g2)
        v.samples[sample].phased = True
        v.samples[sample]["PG"] = chunk_id
        v.samples[sample]["SZ"] = (round(h1_min_path.sizesim, 3) if g1 else None,
                                   round(h2_min_path.sizesim, 3) if g2 else None)
        v.samples[sample]["CS"] = (round(h1_min_path.cossim, 3) if g1 else None,
                                   round(h2_min_path.cossim, 3) if g2 else None)
        v.samples[sample]["AD"] = (h1.coverage, h2.coverage)
        ret_entries.append(v)
    return ret_entries

def get_bounds(cnk):
    """
    Need min/max position of variants
    """
    mstart = sys.maxsize
    mend = 0
    for i in cnk:
        mstart = min(mstart, i.start)
        mend = max(mend, i.stop)
    return mstart, mend

def phase_region(up_variants, hap1, hap2, params, chunk_id=None, sample=0):
    """
    Phase the variants from a region over two haplotypes
    Returns a list of variants
    """
    ret_entries = []
    graph, unused_vars = kdp.vars_to_graph(up_variants, params.kmer)
    unused_cnt = len(unused_vars)
    for entry in unused_vars:
        entry.samples[sample]['GT'] = (0, 0)
        ret_entries.append(entry)

    # No changes to apply
    logging.debug('h1 %s', hap1)
    if hap1.n:
        h1_paths = kdp.find_hap_paths(graph, hap1, params)
    else:
        h1_paths = []
    h1_min_path = kdp.get_best_path(h1_paths, params)
    logging.debug('chose h1 %s', h1_min_path)

    logging.debug('h2 %s', hap2)
    if hap2.n:
        h2_paths = kdp.find_hap_paths(graph, hap2, params)
    else:
        h2_paths = []
    h2_min_path = kdp.get_best_path(h2_paths, params)
    logging.debug('chose h2 %s', h2_min_path)

    used = set(h1_min_path.path) | set(h2_min_path.path)
    ret_entries.extend(pull_variants(
        graph, used, h1_min_path, hap1, h2_min_path, hap2, chunk_id, sample))

    unused_cnt = 0
    for node in set(graph.nodes) - used:
        if node in ['src', 'snk']:
            continue
        v = graph.nodes[node]["variant"]
        v.samples[0]['GT'] = (0, 0)
        ret_entries.append(v)
        unused_cnt += 1
    logging.debug("Cleared %d variants' genotypes", unused_cnt)
    return ret_entries


def kdp_job_vcf(chunk, params, header=None, sample=0):
    """
    Phase a chunk of variants
    """
    chunk_dict, chunk_id = chunk
    ret = []

    comp_entries = chunk_dict['comp']
    if len(comp_entries) == 0:
        return ret

    for entry in comp_entries:  # TODO just masking?
        entry.samples[sample]["GT"] = (None, None)
    
    # No base variants, assuming all comp is covered
    base_entries = chunk_dict['base']
    if len(base_entries) == 0:
        for entry in comp_entries:
            entry.samples[sample]["GT"] = (0, 0)
            ret.append(entry)
        return ret

    hap1, hap2 = kdp.vcf_haps(base_entries, params.kmer)

    # gross but don't have to keep passing around the header
    for _ in comp_entries:
        _.translate(header)

    chunk_id = str(chunk_id)
    for entry in phase_region(comp_entries, hap1, hap2, params, chunk_id, sample):
        ret.append(entry)
    return ret

def kdp_job_bam(chunk, bam, reference, params, header=None, sample=0):
    """
    Phase a chunk of variants
    """
    chunk_dict, chunk_id = chunk
    ret = []

    comp_entries = chunk_dict['comp']
    if len(comp_entries) == 0:
        return ret

    for entry in comp_entries:  # TODO just masking?
        entry.samples[sample]["GT"] = (None, None)
    
    chrom = comp_entries[0].chrom
    n_tries = 1 # turn this off, for now. I think its slow - good proof of concept, though
    while n_tries and comp_entries:
        start, end = get_bounds(comp_entries)
        refseq = reference.fetch(chrom, start - params.chunksize, end + params.chunksize)
        hap1, hap2 = kdp.bam_haps(bam, refseq, chrom, start, end, params)
        n_tries -= 1
        # Neither hap describes a variant, possibly due to bad boundaries from a giant variant
        # First attempt is to just remove the largest variant and then retry
        # Long term solution is that we collect the Haplotypes over the region and record
        # with them the start/end they span. Then, we work to apply those variants over the sub-graph 
        # within that span
        # Like, this whole approach we're currently using with letting the chunk window dictate the variants
        # isn't great. The reads should speak for themselves and then get placed onto the chunk
        # I think I found it https://web.eecs.umich.edu/~dkoutra/papers/18-HashAlign-PAKDD.pdf
        if hap1.n == 0 and hap2.n == 0: 
            # Remove the largest comp_entry
            largest_idx = 0
            for i in range(len(comp_entries)):
                if truvari.entry_size(comp_entries[i]) > truvari.entry_size(comp_entries[largest_idx]):
                    largst_idx = i
            ret.append(comp_entries.pop(largest_idx))
        else:
            n_tries = 0

    if not comp_entries:
        return ret

    # gross but don't have to keep passing around the header
    for _ in comp_entries:
        _.translate(header)

    chunk_id = str(chunk_id)
    for entry in phase_region(comp_entries, hap1, hap2, params, chunk_id, sample):
        ret.append(entry)
    return ret
