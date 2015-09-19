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
from ersa.dbmanager import DbManager, Database


def get_args():
    p = ArgumentParser(description="estimate combined number of generations between pairs of individuals")
    p.add_argument("matchfile", help="input match file")
    p.add_argument("-a", "--alpha", help="significance level (default: %(default).2f)",
                   type=float, default=0.05)
    p.add_argument("-c", help="number of autosomes (default: %(default)d for humans)",
                   type=int, default=22)
    p.add_argument("-ci", help="generate confidence intervals",
                   action='store_true')
    p.add_argument("-d", "--dmax", help="max combined number of generations to test (default: %(default)d)",
                   type=int, default=10)
    p.add_argument("-l", help="mean number of segments shared in the population (default: %(default).1f)",
                   type=float, default=13.73)
    p.add_argument("-r", help="expected number of recombination events per haploid genome per generation (default %(default).1f for humans)",
                   type=float, default=35.2548101)
    p.add_argument("-t", help="min segment length (in cM) to include in comparisons (default %(default).1f)",
                   type=float, default=2.5)
    p.add_argument("-u", "--user", help="filter input file to only look at USER",
                   type=str)
    p.add_argument("-th", "--theta", help="mean shared segment length (in cM) in the population (default %(default).3f)",
                   type=float, default=3.197036753)


    group = p.add_mutually_exclusive_group()
    group.add_argument("-D", help="direct output to database D")
    group.add_argument("-o", "--ofile", help="direct output to OFILE")

    args = p.parse_args()
    return args


def gen_estimates(args, h0, ha, pair_dict):
    for pair, seg_list in pair_dict.items():
        dob = (None, None)  # TODO get dob from file
        s = [seg.length for seg in seg_list]
        n = len(s)
        est = estimate_relation(pair, dob, n, s, h0, ha, args.dmax, args.alpha, args.ci)
        pair1, pair2 = pair.split(':')
        d_est = est.d if est.reject else "NA"
        yield d_est, est, n, pair1, pair2, s, seg_list


def main():
    args = get_args()

    start_time = time()

    print("--- Reading match file ---")

    pair_dict = get_pair_dict(args.matchfile, args.t, args.user)

    h0 = Background(args.t, args.theta, args.l)
    ha = Relation(args.c, args.r, args.t, args.theta, args.l)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))
    print()

    print("--- Solving ---")

    output_file = None
    if not args.D:
        output_file = open(args.ofile, "w") if args.ofile else stdout
        print("{:<20} {:<20} {:<10} {:<10} {:>10} {:>10} {:>10}"
              .format("Indv_1", "Indv_2", "Rel_est1", "Rel_est2", "d_est", "N_seg", "Tot_cM"),
              file=output_file)

    if args.D:
        n_pairs = len(pair_dict)
        print("processing {:,} pairs..".format(n_pairs))
        ests, seg_lists = [], []
        for d_est, est, n, pair1, pair2, s, seg_list in gen_estimates(args, h0, ha, pair_dict):
            ests.append(est)
            seg_lists.append(seg_list)
        print("pushing results to database..")
        with DbManager(args.D) as db:
            db.insert(ests, seg_lists)
    else:
        for d_est, est, n, pair1, pair2, s, seg_list in gen_estimates(args, h0, ha, pair_dict):
            if est.rel_est is None:
                rel_est = ("NA", "NA")
            else:
                rel_est = est.rel_est
            print("{:<20} {:<20} {:10} {:10} {:>10} {:10} {:10,.2f}"
                  .format(pair1, pair2, rel_est[0], rel_est[1], d_est, n, sum(s)),
                  file=output_file)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))

if __name__ == '__main__':
    main()
