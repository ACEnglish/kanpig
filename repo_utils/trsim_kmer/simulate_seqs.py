import sys
import json
import math
import random
import argparse

from collections import Counter, defaultdict

import kanpig
import truvari

def parse_args():
    parser = argparse.ArgumentParser(
        description="Tandem repeat simulation parameters")

    parser.add_argument("--num-experiments", type=int, default=1000)
    parser.add_argument("--motif-length", type=int, default=3)
    parser.add_argument("--n-repeats", type=int, default=40)
    parser.add_argument("--tr-mutation-rate", type=float, default=0.02)
    parser.add_argument("--gc-bias", type=float, default=0.4)
    parser.add_argument("--db-mutation-rate", type=float, default=0.05)
    parser.add_argument("--db-max-motif-expcon", type=int, default=5)
    parser.add_argument("--query-mutation-rate", type=float, default=0.02)
    parser.add_argument("--query-max-motif-expcon", type=int, default=2)
    parser.add_argument("--kmer-size", type=str, default="16,4")
    parser.add_argument("--dual", action="store_true")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--brief", action="store_true")
    parser.add_argument("--write-sim", type=str, default=None)
    args = parser.parse_args()

    if args.seed is None:
        args.seed = random.randint(0, 2**32 - 1)
    
    args.kmer_size = tuple(map(int, args.kmer_size.split(',')))
    return args


def simulate_tandem_repeat_noisy(unit_length=3, n_repeats=40, mutation_rate=0.02, gc_bias=0.4):
    at = (1 - gc_bias) / 2
    gc = gc_bias / 2
    unit = ''.join(random.choices(['A', 'T', 'C', 'G'],
                                  weights=[at, at, gc, gc],
                                  k=unit_length))
    # Each copy of the unit can drift slightly
    copies = [mutate_dna(unit, mutation_rate) for _ in range(n_repeats)]
    ret = ''.join(copies)

    return ret, unit


def simulate_targets(original_seq, motif, n, base_rate=0.05, motif_max=5):
    ret = []
    for i in range(n):
        new_seq = original_seq
        if random.choice([True, False]):
            trim = random.randint(0, motif_max) * len(motif)
            if trim != 0:
                new_seq = new_seq[:-trim]
        else:
            new_seq = new_seq + (motif * random.randint(0, motif_max))
        new_seq = mutate_dna(new_seq, base_rate)
        while new_seq in ret:
            new_seq = mutate_dna(new_seq, 0.005)
        ret.append(new_seq)
    return ret


def mutate_dna(sequence, mutation_rate=0.05):
    bases = "ACGT"
    sequence = list(sequence)

    for i in range(len(sequence)):
        if random.random() < mutation_rate:
            # Pick a different base than the current one
            sequence[i] = random.choice([b for b in bases if b != sequence[i]])

    return "".join(sequence)

def score(query_vec, target_vec, dual=False):
    a = query_vec.fine_similarity(target_vec)
    if not dual:
        return a
    b = query_vec.coarse_similarity(target_vec)
    return math.sqrt(a * b)

if __name__ == '__main__':
    args = parse_args()
    if not args.brief:
        print(json.dumps(vars(args), indent=2))

    if args.write_sim:
        args.write_sim = open(args.write_sim, 'w')
        args.write_sim.write("cansim_q1\tcansim_q2\tseqsim_q1\tseqsim_q2\n")

    random.seed(args.seed)

    invalid = 0  # Both queries shouldn't hit the same target
    odd_balls = 0  # Queries should hit a different target than 0, but don't
    tests = 0  # valid tests run
    same = 0  # tests where both queries hit the same target
    same_correct = 0  # tests where both queries hit the same correct target

    for _ in range(args.num_experiments):
        original_seq, motif = simulate_tandem_repeat_noisy(args.motif_length,
                                                           args.n_repeats,
                                                           args.tr_mutation_rate,
                                                           args.gc_bias)
        if args.debug:
            print(f"Simulated {len(original_seq)}bp TR with motif {motif}:")
            print(original_seq)

        target_seqs = simulate_targets(original_seq,
                                       motif,
                                       10,
                                       args.db_mutation_rate,
                                       args.db_max_motif_expcon)
        if args.debug:
            print("Targets:")
            print(">" + "\n>".join(target_seqs))

        target_vecs = [kanpig.KmerVec(
            t, args.kmer_size) for t in target_seqs]

        # Query is a less noisy version of target_0
        query_seq1, query_seq2 = simulate_targets(target_seqs[0],
                                                  motif,
                                                  2,
                                                  args.query_mutation_rate,
                                                  args.query_max_motif_expcon)
        if args.debug:
            print("Queries")
            print(">" + query_seq1)
            print(">" + query_seq2)

        query_vec1 = kanpig.KmerVec(query_seq1, args.kmer_size)
        query_vec2 = kanpig.KmerVec(query_seq2, args.kmer_size)

        # Establish queries' most similar target
        base_sim1 = truvari.seqsim(query_seq1, target_seqs[0])
        base_sim2 = truvari.seqsim(query_seq2, target_seqs[0])
        best_idx1 = 0
        best_idx2 = 0
        for idx, t in enumerate(target_seqs[1:]):
            sim1 = truvari.seqsim(query_seq1, t)
            if base_sim1 < sim1:
                base_sim1 = sim1
                best_idx1 = idx + 1

            sim2 = truvari.seqsim(query_seq2, t)
            if base_sim2 < sim2:
                base_sim2 = sim2
                best_idx2 = idx + 1

        if best_idx1 != best_idx2:
            invalid += 1
            if args.debug:
                print("Invalid. Queries should hit same target")
            continue
        elif args.debug:
            print(
                f"Best target is {best_idx1} ({base_sim1:.4f}, {base_sim2:.4f})")

        # Do the queries hit the best target
        can_sim1 = score(query_vec1, target_vecs[0], args.dual)
        can_sim2 = score(query_vec2, target_vecs[0], args.dual)
        best_can_idx1 = 0
        best_can_idx2 = 0

        for idx, t in enumerate(target_vecs[1:]):
            sim1 = score(query_vec1, t, args.dual)
            if can_sim1 < sim1:
                can_sim1 = sim1
                best_can_idx1 = idx + 1

            sim2 = score(query_vec2, t, args.dual)
            if can_sim2 < sim2:
                can_sim2 = sim2
                best_can_idx2 = idx + 1

        if args.debug:
            print(
                f"Queries hit {best_can_idx1} ({can_sim1:.4f}) & {best_can_idx2} ({can_sim2:.4f})")

        # Sometimes the classic seqsim will choose a different target from what the queries were
        # simulated from. But then cansim will still find the original target. Odd.
        if best_can_idx1 == best_can_idx2 and best_idx1 != 0 and best_can_idx1 == 0:
            odd_balls += 1
            continue

        # Just the same/correct ones to keep it less noisy
        if args.write_sim and best_can_idx1 == best_can_idx2 == best_idx1:
            args.write_sim.write(
                f"{can_sim1}\t{can_sim2}\t{base_sim1}\t{base_sim2}\n")

        tests += 1
        same += best_can_idx1 == best_can_idx2
        same_correct += best_can_idx1 == best_can_idx2 == best_idx1
        # same_wrong_both += (best_can_idx1 == best_can_idx2) and best_can_idx1 != best_idx1 and best_can_idx2 != best_idx1
        # same_wrong_one += (best_can_idx1 == best_can_idx2) and (best_can_idx1 != best_idx1 ^ best_can_idx2 != best_idx1)
        # different += best_can_idx1 != best_can_idx2

    if args.brief:
        print(invalid, odd_balls, tests, same, same_correct)
        exit(0)

    print('-' * 10)
    width = len(str(tests)) + 1
    print("Invalid: ", f"{invalid:>{width}}")
    print("Odd:     ", f"{odd_balls:>{width}}")
    print("Tests:   ", f"{tests:>{width}}")
    if tests:
        print("Same:    ", f"{same:>{width}}",
              f"{round(same / tests * 100, 1)}%")
        print("&Corr:   ", f"{same_correct:>{width}}",
              f"{round(same_correct / tests * 100, 1)}%")
    else:
        print("Same:    ", f"{0:>{width}}", f"0")
        print("&Corr:   ", f"{0:>{width}}", f"0")
