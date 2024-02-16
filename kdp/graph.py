import numpy as np
import networkx as nx
from functools import total_ordering
from dataclasses import dataclass, field

import truvari

import kdp


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


def vars_to_graph(variants, kmer=3):
    """
    For a sorted set of variants, make a graph
    Returns the digraph (and variants not used?)
    """
    keep_vars = []
    unused_vars = []
    for entry in variants:
        k, s = kdp.make_kfeat(entry, kmer)
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


def get_best_path(paths, exclude=None, min_cos=0.90, min_size=0.90):
    """
    Returns the best path
    """
    to_analyze = filter(lambda tup: tup.sizesim >=
                        min_size and tup.cossim >= min_cos, paths)
    to_analyze = sorted(to_analyze, reverse=True)
    for path in to_analyze:
        # Don't allow paths under size similarity
        # Or those with previously used variants
        if (exclude is None or not set(path.path).intersection(exclude)):
            return path
    return PhasePath()

def graph_phase_paths_orig(graph, hap1_difference, hap1_size, hap2_difference, hap2_size, max_paths=10000):
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

            m_dist = kdp.cosinesim(m_k, cur_hap_diff, m_s)
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



def graph_phase_paths(graph, cur_hap_diff, cur_hap_size, max_paths=10000):
    """
    This will return the paths and we'll let phase_region do the editing/deciding
    So this will return list of PhasePath
    """
    ret = []
    for path in dfs(graph, cur_hap_size):
        m_s = graph.nodes[path[0]]['size']
        for node in path[1:]:
            m_s += graph.nodes[node]['size']

        # ensure same sign (same net effect of deletion/insertion)
        if (cur_hap_size ^ m_s) < 0:
            m_sz = 0
        else:
            m_sz, _ = truvari.sizesim(abs(cur_hap_size), abs(m_s))
        # This is all I've done different
        if m_sz < min_size:
            continue

        m_k = np.copy(graph.nodes[path[0]]['kfeat'])
        for node in path[1:]:
            m_k += graph.nodes[node]['kfeat']

        m_dist = kdp.cosinesim(m_k, cur_hap_diff, m_s)
        ret.append(PhasePath(m_dist, m_sz, path))
        max_paths -= 1
        if max_paths <= 0:
            break
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

    diffs = sorted([(abs(target - (cur_len + g.nodes[n]['size'])), n)
                   for _, n in g.out_edges(cur_node)])
    for len_diff, next_node in diffs:
        if next_node == 'snk' and cur_node != 'src':
            yield list(path)
        else:
            n_path = list(path)
            for sub_path in dfs(g, target, next_node, cur_len, n_path):
                yield sub_path
