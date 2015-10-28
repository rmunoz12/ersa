#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from .ersa_LL import Background, Relation, estimate_relation
from .parser import get_pair_dict
from time import time
from sys import stdout
from argparse import ArgumentParser
from .dbmanager import DbManager


def get_args():
    p = ArgumentParser(description="estimate combined number of generations between pairs of individuals")
    p.add_argument("matchfile", help="input match file")
    p.add_argument("-a", "--alpha", help="significance level (default: %(default).2f)",
                   type=float, default=0.05)
    p.add_argument("--avuncular-adj", help="apply the adjustment to Na from Li et al. (2014) for avuncular (a=2, d=3) relationships",
                   action="store_true")
    p.add_argument("-c", help="number of autosomes (default: %(default)d for humans)",
                   type=int, default=22)
    p.add_argument("-ci", help="generate confidence intervals",
                   action='store_true')
    p.add_argument("-d", "--dmax", help="max combined number of generations to test (default: %(default)d)",
                   type=int, default=10)
    p.add_argument("--first_deg_adj", help="Include adjustments for first-degree relationships",
                   action="store_true")
    p.add_argument("-H", help="input matchfile contains an extra column at the end of each line with haploscores (discarded by ersa)",
                   action='store_true')
    p.add_argument("-l", help="mean number of segments shared in the population (default: %(default).1f)",
                   type=float, default=13.73)
    p.add_argument("--merge-segs", help="merge segments that are on the same chromosome and <= MERGE-SEGS bp apart (default No merge)",
                   type=int, default=-1)
    p.add_argument("--nomask", help="disable genomic region masking",
                   action="store_true")
    p.add_argument("-r", help="expected number of recombination events per haploid genome per generation (default %(default).1f for humans)",
                   type=float, default=35.2548101)
    p.add_argument("-t", help="min segment length (in cM) to include in comparisons (default %(default).1f)",
                   type=float, default=2.5)
    p.add_argument("-u", "--user", help="filter input file to only look at USER",
                   type=str)
    p.add_argument("-th", "--theta", help="mean shared segment length (in cM) in the population (default %(default).3f)",
                   type=float, default=3.197036753)
    p.add_argument("--skip-soft-delete", help="Assume the database is empty, don't soft-delete before inserting new data",
                   action='store_true', default=False)

    group = p.add_mutually_exclusive_group()
    group.add_argument("-D", help="direct output to database D")
    group.add_argument("-o", "--ofile", help="direct output to OFILE")

    group2 = p.add_mutually_exclusive_group()
    group2.add_argument("--insig-threshold", help="Threshold (cM) minimum to keep insignificant results (default: off)",
                        type=float, default=None)
    group2.add_argument("--keep-insig-by-seg", help="Keep insignificant results that have at least <first value> segments shared of the specified size <second value> (default: off)",
                        type=float, default=None, nargs=2)
    group2.add_argument("--keep-insignificant", help="push insignificant results to the database where d_est is NULL (default: discard below INSIG-THRESHOLD)",
                        action='store_true')

    args = p.parse_args()
    return args


def gen_estimates(args, h0, ha, pair_dict):
    """
    Returns
    -------
    (est, seg_list) : (Estimate, list[ersa.parser.SharedSegment])
        Tuple of estimate results and corresponding segment list.
    """
    for pair, seg_list in pair_dict.items():
        dob = (None, None)  # TODO get dob from file
        s = [seg.length for seg in seg_list]
        n = len(s)
        est = estimate_relation(pair, dob, n, s, h0, ha, args.dmax, args.alpha, args.ci)
        yield est, seg_list


def print_LLs(alts):
    print("{", end='')
    for i in range(len(alts)):
        alt = alts[i]
        print("\"{}\": {:>8,.3f}".format(str(alt[0] - 1), round(alt[2], 3)), end='')
        if i < len(alts) - 1:
            print(",", end=' ')
    print("}")

def main():
    args = get_args()

    start_time = time()

    print("--- Reading match file ---")

    pair_dict = get_pair_dict(args.matchfile, args.t, args.user, args.H, args.nomask, args.merge_segs)

    h0 = Background(args.t, args.theta, args.l)
    ha = Relation(args.c, args.r, args.t, args.theta, args.l,
                  args.first_deg_adj, args.nomask, args.avuncular_adj)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))
    print()

    print("--- Solving ---")

    if args.D:
        n_pairs = len(pair_dict)
        print("processing {:,} pairs..".format(n_pairs))
        ests, seg_lists = [], []
        total_segs = 0
        for est, seg_list in gen_estimates(args, h0, ha, pair_dict):
            keep = False
            if args.keep_insig_by_seg:
                n_needed = args.keep_insig_by_seg[0]
                l_needed = args.keep_insig_by_seg[1]
                count = sum(i.length > l_needed for i in seg_list)
                keep = True if count > n_needed else False
            # 'reject' => H0 is rejected, this pair is significant.
            if est.reject or args.keep_insignificant or \
                    (args.insig_threshold and est.cm >= args.insig_threshold) or \
                    (args.keep_insig_by_seg and keep):
                ests.append(est)
                seg_lists.append(seg_list)
                total_segs += len(seg_list)
        print("pushing results from '{}' to database... " \
              "({} pairs, {} segments)".format(args.matchfile, len(ests), total_segs))
        with DbManager(args.D, skip_soft_delete=args.skip_soft_delete) as db:
            db.insert(ests, seg_lists)
    else:
        output_file = open(args.ofile, "w") if args.ofile else stdout
        print("{:<20} {:<20} {:<25} {:>10} {:>10} {:>10} {:>8} {:>8} {}"
              .format("Indv_1", "Indv_2", "Rel_est1", "d_est", "N_seg", "Tot_cM", "Max-P", "NullLL", "LLs"),
              file=output_file)
        for est, seg_list in gen_estimates(args, h0, ha, pair_dict):
            d_est = est.d if est.reject else "NA"
            s = est.cm
            if est.rel_est is None:
                rel_est = ("NA", "NA")
            else:
                rel_est = est.rel_est
            print("{:<20} {:<20} {:25} {:>10} {:10} {:10,.2f} {:8,.3f} {:8,.3f}"
                  .format(est.indv1, est.indv2, rel_est[0], d_est, len(seg_list), s, est.p, est.null_LL),
                  file=output_file, end="  ")
            print_LLs(est.alts)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))
