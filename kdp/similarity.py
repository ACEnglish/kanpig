import numpy as np

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
    return float(score) # have to cast for pysam/VCF
    #return simsimd.cosine(a, b) wasn't faster, didn't work :(


