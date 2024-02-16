
__version__ = '0.1.0-dev'

from kdp.kdp import *

from kdp.similarity import (
        cosinesim,
)

from kdp.kmer import (
    make_kfeat,
    make_kmer,
)

from kdp.graph import (
    PhasePath,
    vars_to_graph,
    get_best_path,
    graph_phase_paths,
)
