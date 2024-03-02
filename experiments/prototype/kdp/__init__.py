
__version__ = '0.1.0-dev'

from kdp.kdp import (
        kdp_job_vcf,
        kdp_job_bam,
)

from kdp.similarity import (
    cosinesim,
    weighted_cosinesim,
)

from kdp.kmer import (
    var_to_kfeat,
    seq_to_kmer,
)

from kdp.graph import (
    PhasePath,
    vars_to_graph,
    get_best_path,
    find_hap_paths,
)

from kdp.haps import (
    Haplotype,
    vcf_haps,
    bam_haps,
)

from kdp.cli import (
        IOParams,
        KDParams,
        parse_args,
)
