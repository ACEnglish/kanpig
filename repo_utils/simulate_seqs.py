import sys
import json
import random
import argparse

from collections import Counter, defaultdict

import kanpig
import truvari

def parse_args():
    parser = argparse.ArgumentParser(description="Tandem repeat simulation parameters")
    
    parser.add_argument("--num-experiments", type=int, default=1000)
    parser.add_argument("--motif-length", type=int, default=3)
    parser.add_argument("--n-repeats", type=int, default=40)
    parser.add_argument("--tr-mutation-rate", type=float, default=0.02)
    parser.add_argument("--gc-bias", type=float, default=0.4)
    parser.add_argument("--db-mutation-rate", type=float, default=0.05)
    parser.add_argument("--db-max-motif-expcon", type=int, default=5)
    parser.add_argument("--query-mutation-rate", type=float, default=0.02)
    parser.add_argument("--query-max-motif-expcon", type=int, default=2)
    parser.add_argument("--mink", type=int, default=1)
    parser.add_argument("--kmer-size", type=int, default=4)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    if args.seed is None:
        args.seed = random.randint(0, 2**32 - 1)

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


if __name__ == '__main__':
    args = parse_args()
    print(json.dumps(vars(args), indent=2))
    random.seed(args.seed)

    tests = 0
    same = 0
    same_correct = 0

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

        target_vecs = [kanpig.seq_to_kmer(t, args.kmer_size, False) for t in target_seqs]

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

        query_vec1 = kanpig.seq_to_kmer(query_seq1, args.kmer_size, False)
        query_vec2 = kanpig.seq_to_kmer(query_seq2, args.kmer_size, False)

        # I should be checking that the querys' seqsim is most similar to target 0
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
            if args.debug:
                print("Invalid. Queries shouldn't hit same target")
            continue
        elif args.debug:
            print(f"Best target is {best_idx1} ({base_sim1:.4f}, {base_sim2:.4f})")
        
        # Do the queries hit the best target
        base_sim1 = kanpig.cansim(query_vec1, target_vecs[0], args.mink)
        base_sim2 = kanpig.cansim(query_vec2, target_vecs[0], args.mink)
        best_can_idx1 = 0
        best_can_idx2 = 0

        for idx, t in enumerate(target_vecs[1:]):
            sim1 = kanpig.cansim(query_vec1, t, args.mink)
            if base_sim1 < sim1:
                base_sim1 = sim1
                best_can_idx1 = idx + 1

            sim2 = kanpig.cansim(query_vec1, t, args.mink)
            if base_sim2 < sim1:
                base_sim1 = sim1
                best_can_idx2 = idx + 1

        if args.debug:
            print(f"Queries hit {best_can_idx1} ({base_sim1:.4f}) & {best_can_idx2} ({base_sim2:.4f})")

        tests += 1
        same += best_can_idx1 == best_can_idx2
        same_correct += best_can_idx1 == best_can_idx2 == best_idx1
        #same_wrong_both += (best_can_idx1 == best_can_idx2) and best_can_idx1 != best_idx1 and best_can_idx2 != best_idx1
        #same_wrong_one += (best_can_idx1 == best_can_idx2) and (best_can_idx1 != best_idx1 ^ best_can_idx2 != best_idx1)
        #different += best_can_idx1 != best_can_idx2
    
    print('-' * 10)
    width = len(str(tests)) + 1
    print("Tests:", f"{tests:>{width}}")
    if tests:
        print("Same: ", f"{same:>{width}}", f"{round(same / tests * 100, 1)}%")
        print("&Corr:", f"{same_correct:>{width}}", f"{round(same_correct / tests * 100, 1)}%")
                

"""Notes
python simulate_seqs.py --num-experiments 1 --debug  --n-repeats 20 --tr-mutation-rate 0.05  --motif-length 6 --mink 1 --kmer-size 7 --seed 674768669

This one is weird. Kanpig approach hits 0 the best, seqsim on 3
"""
