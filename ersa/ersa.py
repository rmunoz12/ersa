#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from argparse import ArgumentParser
from datetime import datetime
import json
from time import time
from tqdm import tqdm
from .dbmanager import DbManager
from .ersa_LL import Background, Relation, estimate_relation
from .json_manager import StreamJSON
from .parser import get_pair_dict


def get_args():
    p = ArgumentParser(description="estimate combined number of generations between pairs "
                                   "of individuals")
    p.add_argument("matchfile", help="input match file")
    p.add_argument("-a", "--alpha", help="significance level (default: %(default).2f)",
                   type=float, default=0.05)
    p.add_argument("--avuncular-adj", help="apply the adjustment to Na from Li et al. (2014) "
                                           "for avuncular (a=2, d=3) relationships",
                   action="store_true")
    p.add_argument("-c", help="number of autosomes (default: %(default)d for humans)",
                   type=int, default=22)
    p.add_argument("-ci", help="generate confidence intervals",
                   action='store_true')
    p.add_argument("-d", "--dmax", help="max combined number of generations to test "
                                        "(default: %(default)d)",
                   type=int, default=10)
    p.add_argument("--first_deg_adj", help="Include adjustments for first-degree relationships",
                   action="store_true")
    p.add_argument("-H", help="input matchfile contains an extra column at the end of each "
                              "line with haploscores (discarded by ersa)",
                   action='store_true')
    p.add_argument("-l", help="mean number of segments shared in the population "
                              "(default: %(default).1f)",
                   type=float, default=13.73)
    p.add_argument("--merge-segs", help="merge segments that are on the same chromosome "
                                        "and <= MERGE-SEGS bp apart (default No merge)",
                   type=int, default=-1)
    p.add_argument("--nomask", help="disable genomic region masking",
                   action="store_true")
    p.add_argument("--progress", help="show progress bar during relation estimation",
                   action="store_true")
    p.add_argument("-r", help="expected number of recombination events per haploid genome "
                              "per generation (default %(default).1f for humans)",
                   type=float, default=35.2548101)
    p.add_argument("-t", help="min segment length (in cM) to include in comparisons "
                              "(default %(default).1f)",
                   type=float, default=2.5)
    p.add_argument("-u", "--user", help="filter input file to only look at USER",
                   type=str)
    p.add_argument("-th", "--theta", help="mean shared segment length (in cM) in the population "
                                          "(default %(default).3f)",
                   type=float, default=3.197036753)
    p.add_argument("--skip-soft-delete", help="Assume the database is empty, don't "
                                              "soft-delete before inserting new data",
                   action='store_true', default=False)

    group = p.add_mutually_exclusive_group()
    group.add_argument("-D", help="direct output to database D")
    group.add_argument("-o", "--ofile", help="direct output to a JSON file, OFILE",
                       type=str, default=None)

    group2 = p.add_mutually_exclusive_group()
    group2.add_argument("--insig-threshold", help="Threshold (cM) minimum to keep "
                                                  "insignificant results (default: off)",
                        type=float, default=None)
    group2.add_argument("--keep-insig-by-seg", help="Keep insignificant results that have "
                                                    "at least <first value> segments shared "
                                                    "of the specified size <second value> "
                                                    "(default: off)",
                        type=float, default=None, nargs=2)
    group2.add_argument("--keep-insignificant", help="push insignificant results to the "
                                                     "database where d_est is NULL "
                                                     "(default: discard below INSIG-THRESHOLD)",
                        action='store_true')

    args = p.parse_args()
    return args


def gen_estimates(args, h0, ha, pair_dict):
    """
    Generate estimates for a set of pairs.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments. The following variables
        are used.

        args.progress : bool | None
            If True, show the progress bar.

        args.dmax : int
            Maximum degree of relation to test.

        args.alpha : float
            Significance level for likelihood ratio test.

        args.ci : bool | None
            Controls whether confidence intervals are calculated.

    h0 : Background
        Population background likelihood calculator.

    ha : Relation
        Recent relationship likelihood calculator.

    pair_dict : dict[str: list[SharedSegments]]
        Map of pairs to segments shared by each pair.

    Yields
    -------
    (est, seg_list) : (Estimate, list[ersa.parser.SharedSegment])
        Tuple of estimate results and corresponding segment list.
    """
    g = pair_dict.items()
    if args.progress:
        g = tqdm(g, total=len(pair_dict))
    for pair, seg_list in g:
        dob = (None, None)  # TODO get dob from file
        s = [seg.length for seg in seg_list]
        n = len(s)
        est = estimate_relation(pair, dob, n, s, h0, ha, args.dmax, args.alpha, args.ci)
        yield est, seg_list


def filter_estimates(args, g):
    """
    Filter insignificant estimates generated by g according
    to command-line arguments. Significant estimates are never
    filtered out. Specific cases filtered are described under the
    args parameter.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments. The following variables control
        which insignificant results are filtered out.

        args.keep_insignificant : bool | None
            If True, keep all insignificant results.

        args.insig_threshold : float | None
            If not None and total cM shared is >= args.insig_threshold
            then the insignificant result is kept.

        args.keep_insig_by_seg : (int, float) | None
            If not None, then the first value controls the number of
            "long" segments needed for the insignificant result to
            be kept, where "long" segments are defined by the second
            value. Let args.keep_insig_by_seg = (n, l), then the
            insignificant result is kept if the number of shared
            segments with length greater than l is greater than n.

    g : T -> generator[(Estimate, list[ersa.parser.SharedSegment])]
        Function that yields ancestry estimates.

    Yields
    -------
    (est, seg_list) : (Estimate, list[ersa.parser.SharedSegment])
        Tuple of estimate results and corresponding segment list.
    """
    for est, seg_list in g:
        keep = False
        if args.keep_insig_by_seg:
            n_needed = args.keep_insig_by_seg[0]
            l_needed = args.keep_insig_by_seg[1]
            count = sum(i.length > l_needed for i in seg_list)
            keep = True if count > n_needed else False
        if est.reject or args.keep_insignificant or \
                (args.insig_threshold and est.cm >= args.insig_threshold) or \
                (args.keep_insig_by_seg and keep):
            yield est, seg_list


def estimates_to_dicts(g):
    """
    Generate a dictionary for each ancestry estimate.

    Parameters
    ----------
    g : T -> generator[(Estimate, list[ersa.parser.SharedSegment])]
        Function that yields ancestry estimates.

    Yields
    -------
    line : dict[str, T]
        Dictionary form of an estimate.
    """
    for est, seg_list in g:
        total_bp = 0
        for seg in seg_list:
            total_bp += seg.bpEnd - seg.bpStart + 1
        line = {"pair": est.pair, "indv1": est.indv1, "indv2": est.indv2,
                "d_est": est.d, "rel_est1": est.rel_est1, "rel_est2": est.rel_est2,
                "n": len(est.s), "total_cM": est.cm, "total_bp": total_bp,
                "LLs": est.LLs, "na": (len(est.s) - est.np),
                "created_date": str(datetime.utcnow()),
                "segments": [{'chromosome': seg.chrom, 'length': seg.length,
                              'bp_start': seg.bpStart, 'bp_end': seg.bpEnd}
                             for seg in seg_list],
                "max_model_p": est.p, "reject": est.reject}
        yield line


def main():
    args = get_args()
    start_time = time()

    print("--- Reading match file ---")

    pair_dict = get_pair_dict(args.matchfile, args.t, args.user,
                              args.H, args.nomask, args.merge_segs)

    h0 = Background(args.t, args.theta, args.l)
    ha = Relation(args.c, args.r, args.t, args.theta, args.l,
                  args.first_deg_adj, args.nomask, args.avuncular_adj)

    print("--- {} seconds ---\n".format(round(time() - start_time, 3)))

    print("--- Solving ---")

    print("processing {:,} pairs..".format(len(pair_dict)))
    ests, seg_lists = [], []
    total_segs = 0
    g = gen_estimates(args, h0, ha, pair_dict)
    g = filter_estimates(args, g)
    if args.D:
        for est, seg_list in g:
            ests.append(est)
            seg_lists.append(seg_list)
            total_segs += len(seg_list)
        print("pushing results from '{}'... "
              "({} pairs, {} segments)".format(args.matchfile, len(ests), total_segs))
        with DbManager(args.D, skip_soft_delete=args.skip_soft_delete) as db:
            db.insert(ests, seg_lists)
    else:
        est_stream = StreamJSON(estimates_to_dicts, g=g)
        if args.ofile:
            with open(args.ofile, 'w') as out:
                json.dump(est_stream, out)
        else:
            print(json.dumps(est_stream, indent=4))

    print("--- {} seconds ---".format(round(time() - start_time, 3)))
