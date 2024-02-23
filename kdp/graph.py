import logging
import numpy as np
import networkx as nx
from functools import total_ordering
from dataclasses import dataclass, field

import truvari

import kdp

@total_ordering
@dataclass(order=False)
class PhasePath():
    """
    Holds the path/similarities
    """
    sizesim: float = 0
    cossim: float = 0
    path: list = field(default_factory=list)

    def __lt__(self, other):
        # if the size similarity is very close, then try to pick the more similar one
        # What I really need to be doing is trying to take into account the length of the path, as well
        #if abs(self.sizesim - other.sizesim) < 0.005:
            #return self.cossim < other.cossim
        # This weighing is questionable
        #return (self.sizesim / len(self.path)) < (other.sizesim / len(other.path))
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
        k, s = kdp.var_to_kfeat(entry, kmer)
        if k.sum() != 0:
            keep_vars.append((truvari.entry_to_hash(entry), entry, k, s))
        else:
            logging.debug("Removing %s", entry)
            unused_vars.append(entry)
    graph = nx.DiGraph()
    graph.add_node('src', size=0)
    graph.add_node('snk', size=0)
    for key, v, k, s in keep_vars:
        logging.debug("%s %s", key, str(v))
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


def get_best_path(paths, params, exclude=None):
    """
    Returns the best path, this used to do more work when we tried to track pg
    """
    to_analyze = filter(lambda tup: tup.sizesim >=
                        params.pctsize and tup.cossim >= params.cossim, paths)
    to_analyze = sorted(to_analyze, reverse=True)
    for path in to_analyze:
        # Don't allow paths under size similarity
        # Or those with previously used variants
        if (exclude is None or not set(path.path).intersection(exclude)):
            return path
    return PhasePath()

def find_hap_paths_align(graph, hap, params):
    """
    Alignment based path finding
    Returns a list of PhasePaths (actually should just return 1 - the best)
    1. build a matrix that has the similarity of every node against every hap
    2. Every one that matches above the threshold we keep
    3. Every doubly used hap need to pick a single best
    4. We then need to start at the snk and go to the first match.
    5. if the 
    """
    pass

def find_hap_paths(graph, hap, params):
    """
    This will return the paths and we'll let phase_region do the editing/deciding
    So this will return list of PhasePath
    """
    ret = []
    max_paths = params.maxpaths
    for path in dfs(graph, hap.size):
        if max_paths <= 0:
            break
        max_paths -= 1
        m_s = graph.nodes[path[0]]['size']
        for node in path[1:]:
            m_s += graph.nodes[node]['size']

        # ensure same sign (same net effect of deletion/insertion)
        if (hap.size ^ m_s) < 0:
            m_sz = 0
        else:
            m_sz, _ = truvari.sizesim(abs(hap.size), abs(m_s))

        # No need to waste time on this variant
        if m_sz < params.pctsize:
            continue

        m_k = np.copy(graph.nodes[path[0]]['kfeat'])
        for node in path[1:]:
            m_k += graph.nodes[node]['kfeat']

        if abs(m_s) < params.wcoslen:
            m_dist = kdp.weighted_cosinesim(m_k, hap.kfeat)
        else:
            m_dist = kdp.cosinesim(m_k, hap.kfeat)

        ret.append(PhasePath(m_dist, m_sz, path))
        logging.debug('found %s', ret[-1])
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
