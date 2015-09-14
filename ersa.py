#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from ersa.ersa_LL import Background, Relation, estimate_relation
from ersa.parser import get_pair_dict
from time import time
from sys import stdout
from argparse import ArgumentParser


"""Default parameter values"""
h = 10                  # in cM

def get_args():
    p = ArgumentParser(description="estimate combined number of generations between pairs of individuals")
    p.add_argument("matchfile", help="input match file")
    p.add_argument("-a", "--alpha", help="significance level (default: %(default).2f)",
                   type=float, default=0.05)
    p.add_argument("-c", help="number of autosomes (default: %(default)d for humans)",
                   type=int, default=22)
    p.add_argument("-d", "--dmax", help="max combined number of generations to test (default: %(default)d)",
                   type = int, default=10)
    p.add_argument("-l", help="mean number of segments shared in the population (default: %(default).1f)",
                   type=float, default=13.73)
    p.add_argument("-o", "--ofile", help="direct output to OFILE")
    p.add_argument("-r", help="expected number of recombination events per haploid genome per generation (default %(default).1f for humans)",
                   type=float, default=35.2548101)
    p.add_argument("-t", help="min segment length (in cM) to include in comparisons (default %(default).1f)",
                   type=float, default=2.5)
    p.add_argument("-th", "--theta", help="mean shared segment length (in cM) in the population (default %(default).3f)",
                   type=float, default=3.197036753)

    args = p.parse_args()
    return args

def main():
    args = get_args()

    start_time = time()

    print("--- Reading match file ---")

    pair_dict = get_pair_dict(args.matchfile, args.t, h)

    h0 = Background(args.t, args.theta, args.l)
    ha = Relation(args.c, args.r, args.t, args.theta, args.l)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))
    print()

    print("--- Solving ---")

    if args.ofile:
        output_file = open(args.ofile, "w")
    else:
        output_file = stdout

    print("Individual_1\tIndividual_2\tNumber_Meioses\tRelatedness\tLLn\tLlr\tlower_d_CI\tupper_d_CI",file=output_file)
    
    for pair, (n, s) in pair_dict.items():
        null_LL, max_alt_LL, d, reject, lower_d, upper_d = estimate_relation(n, s, h0, ha, args.dmax, args.alpha)
        pair1, pair2 = pair.split(':')
        print(pair1, '\t', pair2, '\t', d, '\t', reject, '\t', null_LL, '\t', max_alt_LL, '\t', lower_d, '\t', upper_d, file=output_file)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))

if __name__ == '__main__':
    main()
