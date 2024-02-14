"""
Approximate phasing of a comparison VCF to match a phased baseline VCF using kmer featurization distances

Todo: 
Validate Inputs e.g. (need tbi)
Instead of max_nodes, I can probably do a max-paths
I think an A* could be faster if I can use the sizediff of current path from hap-path as the heuristic
Figure out if I can do something about the size-separator. Like, it worked good for 'isn't random', but
    What I really need is "are the same"
"""
import re
import sys
import logging
import argparse
import itertools
import multiprocessing
from functools import partial, total_ordering
from collections import defaultdict
from dataclasses import dataclass, field

import pysam
import truvari

import numpy as np
import networkx as nx


def parse_args(args):
    """
    """
    parser = argparse.ArgumentParser(prog="kfdphase", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--base", type=str, required=True,
                        help="Phased VCF")
    parser.add_argument("-c", "--comp", type=str, required=True,
                        help="VCF to phase")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output VCF (stdout)")
    parser.add_argument("-r", "--regions", type=str, default=None,
                        help="Bed filename or comma-separated list of chrom:start-end regions to process")
    parser.add_argument("--kmer", type=int, default=3,
                        help="Kmer size(%(default)s)")
    parser.add_argument("-s", "--sample", type=str, default=None,
                        help="Name of sample to apply genotypes (first)")
    parser.add_argument("--passonly", action="store_true",
                        help="Only phase passing variants")
    parser.add_argument("--sizemin", type=int, default=0,
                        help="Minimum variant size (%(default)s)")
    parser.add_argument("--sizemax", type=int, default=50000,
                        help="Maximum variant size (%(default)s)")
    parser.add_argument("--maxpaths", type=int, default=10000,
                        help="Stop region processing after trying N paths (%(default)s)")
    parser.add_argument("--cossim", type=float, default=0.90,
                        help="Minimum cosine similarity (%(default)s)")
    parser.add_argument("--pctsize", type=truvari.restricted_float, default=0.90,
                        help="Minimum size similarity between a path and haplotype (%(default)s)")
    parser.add_argument("--pg", action="store_true",
                        help="Allow multiple phase groups (%(default)s)")
    parser.add_argument("--chunksize", type=truvari.restricted_int, default=500,
                        help="Max reference distance to phase calls (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args(args)
    return args


def encode_nuc(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'G':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'G':
        return 3
    return 0


def generate_kmers(sequence, kmer=3):
    result = 0
    for pos, i in enumerate(sequence[:kmer]):
        result += encode_nuc(i) << ((kmer - pos - 1) * 2)
    yield result
    mask = (4**(kmer-1)) - 1
    for i in sequence[1:]:
        nuc = encode_nuc(i)
        result = ((result & mask) << 2) + nuc
        yield result


def make_kmer(seq, kmer_len=6):
    """
    Make the kmer array of all kmers and those over min_freq
    """
    ret = np.zeros(4**kmer_len)
    for i in generate_kmers(seq.upper(), kmer_len):
        ret[i] += 1
    return ret


def get_size_dist(size):
    """
    Get given size's distribution's mean, standard deviation 
    This helps cosine similarity of longer sequences be more meaningful for 
    separating two somewhat similar sequences from two random sequences. However
    most 'partial' SVs are somewhat similar and what we're interested in here is
    an approximation of sequence similarity
    """
    def regress(x, coef, inter):
        powers_of_x = np.array([x**i for i in range(0, len(coef))]).T
        return np.dot(powers_of_x, coef) + inter
    mean = regress(size,
                   [0.00000000e+00, 3.43352323e-04, -7.11063537e-08,
                       7.23861424e-12, -2.77881017e-16],
                   0.03197948132428041)
    std = regress(size,
                  [0.00000000e+00, -6.15350480e-06,
                      7.67528930e-10, -3.54563342e-14],
                  0.025929011228045223)
    return mean, std


def cosinesim(a, b, size):
    """
    How many standard deviations is the score from random sequence given the size (turned off)
    """
    score = abs(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    #mean, std = get_size_dist(size)
    #print( (score - mean) / std )
    return score


def make_kfeat(entry, kmer=3):
    """
    Make the kmer featurization of this variant
    """
    alt = make_kmer(entry.alts[0], kmer)
    ref = make_kmer(entry.ref, kmer)
    szdiff = len(entry.alts[0]) - len(entry.ref)
    return alt - ref, szdiff


def vars_to_graph(variants, kmer=3):
    """
    For a sorted set of variants, make a graph
    Returns the digraph (and variants not used?)
    """
    keep_vars = []
    unused_vars = []
    for entry in variants:
        k, s = make_kfeat(entry, kmer)
        if k.sum() != 0:
            keep_vars.append((truvari.entry_to_hash(entry), entry, k, s))
        else:
            unused_vars.append(entry)
    graph = nx.DiGraph()
    graph.add_node('src', size=0)
    graph.add_node('snk', size=0)
    for key, v, k, s in keep_vars:
        graph.add_node(key, variant=v, kfeat=k, size=s)
        graph.add_edge('src', key)
        graph.add_edge(key, 'snk')

    # link the variants
    for i in range(len(keep_vars) - 1):
        up_key, up_var, _, _ = keep_vars[i]
        up = truvari.entry_boundaries(up_var)
        for j in range(i + 1, len(keep_vars)):
            dn_key, dn_var, _, _ = keep_vars[j]
            dn = truvari.entry_boundaries(dn_var)
            if not truvari.overlaps(up[0], up[1], dn[0], dn[1]):
                graph.add_edge(up_key, dn_key)
    return graph, unused_vars


def phased_k(variants, kmer):
    """
    """
    h1_k = make_kmer("", kmer)
    h2_k = make_kmer("", kmer)
    h1_s = 0
    h2_s = 0
    for v in variants:
        k, sz = make_kfeat(v, kmer)
        if v.samples[0]['GT'][0] == 1:
            h1_k += k
            h1_s += sz
        if len(v.samples[0]['GT']) > 1 and v.samples[0]['GT'][1] == 1:
            h2_k += k
            h2_s += sz
    return h1_k, h1_s, h2_k, h2_s


@total_ordering
@dataclass
class PhasePath():
    """
    Holds the path/similarities
    """
    cossim: float = 0
    sizesim: float = 0
    path: list = field(default_factory=list)

    def __lt__(self, other):
        # Trues are always worth more
        if round(self.sizesim, 4) == round(other.sizesim, 4):
            return self.cossim < other.cossim
        return self.sizesim < other.sizesim

    def __eq__(self, other):
        return self.sizesim == other.sizesim and self.cossim == other.cossim


def graph_phase_paths(graph, hap1_difference, hap1_size, hap2_difference, hap2_size, max_paths=10000):
    """
    This will return the paths and we'll let phase_region do the editing/deciding
    So this will return list of PhasePath
    """
    ret = []
    for cur_hap_diff, cur_hap_size in [(hap1_difference, hap1_size), (hap2_difference, hap2_size)]:
        cur_paths = []
        for path in dfs(graph, cur_hap_size):
            m_k = np.copy(graph.nodes[path[0]]['kfeat'])
            m_s = graph.nodes[path[0]]['size']
            for node in path[1:]:
                m_k += graph.nodes[node]['kfeat']
                m_s += graph.nodes[node]['size']

            m_dist = cosinesim(m_k, cur_hap_diff, m_s)
            # ensure same sign (same net effect of deletion/insertion)
            if (cur_hap_size ^ m_s) < 0:
                m_sz = 0
            else:
                m_sz, _ = truvari.sizesim(abs(cur_hap_size), abs(m_s))

            cur_paths.append(PhasePath(m_dist, m_sz, path))
            max_paths -= 1
            if max_paths <= 0:
                break
        ret.append(cur_paths)
    return ret

def dfs(g, target, cur_node=None, cur_len=0, path=None):
    """
    Yield paths with DFS with traversal guided by length difference from the target
    """
    if not cur_node:
        cur_node = 'src'
        path = []
    else:
        path.append(cur_node)

    cur_len += g.nodes[cur_node]['size']
    
    diffs = sorted([(abs(target - (cur_len + g.nodes[n]['size'])), n) for _, n in g.out_edges(cur_node)])
    for len_diff, next_node in diffs:
        if next_node == 'snk' and cur_node != 'src':
            yield list(path)
        else:
            n_path = list(path)
            for sub_path in dfs(g, target, next_node, cur_len, n_path):
                yield sub_path

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
        v.samples[sample]['SZ'] = (round(h1_min_path.sizesim, 3), round(h2_min_path.sizesim, 3))
        v.samples[sample]['CS'] = (round(h1_min_path.cossim, 3), round(h2_min_path.sizesim, 3))
        ret_entries.append(v)
    return ret_entries


def get_best_path(paths, exclude=None, min_cos=0.90, min_size=0.90):
    """
    Returns the best path
    """
    to_analyze = filter(lambda tup: tup.sizesim >= min_size and tup.cossim >= min_cos, paths)
    to_analyze = sorted(to_analyze, reverse=True)
    for path in to_analyze:
        # Don't allow paths under size similarity
        # Or those with previously used variants
        if (exclude is None or not set(path.path).intersection(exclude)):
            return path
    return PhasePath()


def phase_region(up_variants, p_variants, pg=False, chunk_id=None, kmer=3, min_cos=0.90, min_size=0.90, sample=0,
                 max_paths=10000):
    """
    Phase the variants from a region over the reference and haplotypes
    Returns a list of variants
    """
    ret_entries = []
    hap1_difference, hap1_size, hap2_difference, hap2_size = phased_k(p_variants, kmer)
    graph, unused_vars = vars_to_graph(up_variants, kmer)
    unused_cnt = len(unused_vars)
    for entry in unused_vars:
        entry.samples[sample]['GT'] = (0, 0)
        ret_entries.append(entry)

    # num_paths_approx = math.factorial(len(graph.nodes) - 2) + 1
    
    all_paths = graph_phase_paths(
        graph, hap1_difference, hap1_size, hap2_difference, hap2_size, max_paths)
    h1_min_path = get_best_path(all_paths[0], min_cos=min_cos, min_size=min_size)
    h2_min_path = get_best_path(all_paths[1], min_cos=min_cos, min_size=min_size)

    used = set(h1_min_path.path) | set(h2_min_path.path)
    m_chunk_id = chunk_id
    if pg:
        m_chunk_id += f'.0'
    ret_entries.extend(pull_variants(graph, used, h1_min_path, h2_min_path, m_chunk_id, sample))

    # while multiple-phase groups is on and we've been applying variants somewhere
    p_cnt = 0
    while pg and (h1_min_path.path or h2_min_path.path):
        p_cnt += 1
        m_chunk_id = chunk_id + f'.{p_cnt}'
        h1_min_path = get_best_path(0, all_paths, used, min_cos=min_cos, min_size=min_size)
        h2_min_path = get_best_path(1, all_paths, used, min_cos=min_cos, min_size=min_size)
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


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    truvari.setup_logging(args.debug)
    logging.info("Starting")

    base = pysam.VariantFile(args.base)
    comp = pysam.VariantFile(args.comp)

    matcher = truvari.Matcher()
    matcher.params.passonly = args.passonly
    matcher.params.sizefilt = args.sizemin
    matcher.params.sizemin = args.sizemin
    matcher.params.sizemax = args.sizemax
    matcher.params.chunksize = args.chunksize

    # v4.2.0 <=
    #regions = truvari.RegionVCFIterator(base, comp,
                                        #args.regions,
                                        #matcher.params.sizemax)
    #regions.merge_overlaps()
    #base_i = regions.iterate(base)
    #comp_i = regions.iterate(comp)
    
    # v4.2.1
    region_tree = truvari.build_region_tree(base, comp, args.regions)
    truvari.merge_region_tree_overlaps(region_tree)
    base_i = truvari.region_filter(base, region_tree)
    comp_i = truvari.region_filter(comp, region_tree)

    chunks = truvari.chunker(matcher, ('base', base_i), ('comp', comp_i))

    out_name = args.output
    if out_name.endswith(".vcf.gz"):
        out_name = args.output[:-len(".gz")]

    header = comp.header.copy()
    header.add_line(('##FORMAT=<ID=SZ,Number=R,Type=Float,'
                    'Description="Size similarity of phase group">'))
    header.add_line(('##FORMAT=<ID=CS,Number=R,Type=Float,'
                    'Description="Cosine similarity of phase group">'))
    header.add_line(('##FORMAT=<ID=PG,Number=1,Type=String,'
                    'Description="Phase group id from kdp">'))
    args.sample = 0 if args.sample is None else args.sample
    out_vcf = pysam.VariantFile(out_name, 'w', header=header)
    task = partial(kdfp_job,
                   max_paths=args.maxpaths,
                   phase_groups=args.pg,
                   header=header,
                   kmer=args.kmer,
                   min_cos=args.cossim,
                   min_size=args.pctsize,
                   sample=args.sample)

    for variant in itertools.chain.from_iterable(map(task, chunks)):
        out_vcf.write(variant)
    out_vcf.close()

    if args.output.endswith(".gz"):
        truvari.compress_index_vcf(out_name)
    logging.info("Finished")
