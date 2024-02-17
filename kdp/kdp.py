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

import kdp


def phased_k(variants, kmer):
    """
    """
    h1_k = kdp.seq_to_kmer("", kmer)
    h2_k = kdp.seq_to_kmer("", kmer)
    h1_s = 0
    h2_s = 0
    h1_n = 0
    h2_n = 0
    for v in variants:
        k, sz = kdp.var_to_kfeat(v, kmer)
        if v.samples[0]['GT'][0] == 1:
            h1_k += k
            h1_s += sz
            h1_n += 1
        if len(v.samples[0]['GT']) > 1 and v.samples[0]['GT'][1] == 1:
            h2_k += k
            h2_s += sz
            h2_n += 1
    return h1_k, h1_s, h1_n, h2_k, h2_s, h2_n


def pull_variants(graph, used, h1_min_path, h2_min_path, chunk_id, sample=0):
    """
    Update the variants in the two paths and return them
    """
    ret_entries = []
    for node in used:
        g1 = int(node in h1_min_path.path)
        g2 = int(node in h2_min_path.path)
        v = graph.nodes[node]["variant"]
        v.samples[sample]['GT'] = (g1, g2)
        v.samples[sample].phased = True
        v.samples[sample]["PG"] = chunk_id
        v.samples[sample]['SZ'] = (round(h1_min_path.sizesim, 3) if g1 else None,
                                   round(h2_min_path.sizesim, 3) if g2 else None)
        v.samples[sample]['CS'] = (round(h1_min_path.cossim, 3) if g1 else None,
                                   round(h2_min_path.cossim, 3) if g2 else None)
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

def phase_region(up_variants, p_variants, pg=False, chunk_id=None, kmer=3, min_cos=0.90, min_size=0.90, sample=0,
                 max_paths=10000):
    """
    Phase the variants from a region over the reference and haplotypes
    Returns a list of variants
    """
    ret_entries = []
    hap1_k, hap1_size, hap1_n, hap2_k, hap2_size, hap2_n = phased_k(p_variants, kmer)
    graph, unused_vars = kdp.vars_to_graph(up_variants, kmer)
    unused_cnt = len(unused_vars)
    for entry in unused_vars:
        entry.samples[sample]['GT'] = (0, 0)
        ret_entries.append(entry)

    logging.info(get_bounds(up_variants))

    # TODO: I need to exclude haps without any variants. It causes spurious FPs
    # in 'balanced' events
    if hap1_n: # Just test this
        h1_paths = kdp.find_hap_paths(graph,
                                    hap1_k,
                                    hap1_size,
                                    min_size,
                                    max_paths)
    else:
        h1_paths = []
    h1_min_path = kdp.get_best_path(h1_paths,
                                    min_cos=min_cos,
                                    min_size=min_size)

    if hap2_n: # Just test this
        h2_paths = kdp.find_hap_paths(graph,
                                      hap2_k,
                                      hap2_size,
                                      min_size,
                                      max_paths)
    else:
        h2_paths = []
    h2_min_path = kdp.get_best_path(h2_paths,
                                    min_cos=min_cos,
                                    min_size=min_size)

    used = set(h1_min_path.path) | set(h2_min_path.path)
    m_chunk_id = chunk_id
    if pg:
        m_chunk_id += f'.0'
    ret_entries.extend(pull_variants(
        graph, used, h1_min_path, h2_min_path, m_chunk_id, sample))

    # while multiple-phase groups is on and we've been applying variants somewhere
    p_cnt = 0
    while pg and (h1_min_path.path or h2_min_path.path):
        p_cnt += 1
        m_chunk_id = chunk_id + f'.{p_cnt}'
        h1_min_path = kdp.get_best_path(h1_paths,
                                        exclude=used,
                                        min_cos=min_cos,
                                        min_size=min_size)
        h2_min_path = kdp.get_best_path(h2_paths,
                                        exclude=used,
                                        min_cos=min_cos,
                                        min_size=min_size)
        new_used = set(h1_min_path.path) | set(h2_min_path.path)
        ret_entries.extend(pull_variants(
            graph, new_used, h1_min_path, h2_min_path, m_chunk_id, sample))
        # More variants to work on
        used |= new_used

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


def kdfp_job(chunk, max_paths=10000, phase_groups=False, header=None, kmer=3, min_cos=0.90, min_size=0.90, sample=0):
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

    base_entries = chunk_dict['base']
    if len(base_entries) == 0:
        for entry in comp_entries:
            entry.samples[sample]["GT"] = (0, 0)
            ret.append(entry)
        return ret

    # gross but don't have to keep passing around the header
    for _ in comp_entries:
        _.translate(header)

    chunk_id = str(chunk_id)
    for entry in phase_region(comp_entries, base_entries, phase_groups, chunk_id, kmer=kmer, min_cos=min_cos,
                              min_size=min_size, sample=sample, max_paths=max_paths):
        ret.append(entry)
    return ret
