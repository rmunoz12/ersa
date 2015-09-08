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
from sys import argv


"""Default parameter values"""
MAX_D = 40

t = 2.5                 # in cM
h = 10                  # in cM
theta = 3.197036753     # in cM
lambda_ = 13.73         #
r = 35.2548101          # ~for humans
c = 22                  # human autosomes


def main():
    argc = len(argv)
    assert 1 < argc < 3  # TODO add command line options
    path = argv[1]

    start_time = time()

    print("--- Reading match file ---")

    pair_dict = get_pair_dict(path, t, h)

    h0 = Background(t, theta, lambda_)
    ha = Relation(c, r, t, theta, lambda_)

    print("--- {} seconds ---".format(round(time() - start_time, 3)))
    print()

    print("--- Solving ---")

    output_file = open("test.out", "w")
    for pair, (n, s) in pair_dict.items():
        null_LL, max_alt_LL, d, reject = estimate_relation(pair, n, s, h0, ha, MAX_D)
        print("Pair: {}\td: {}\tReject: {}\tLLn: {}\tLLr: {}".
              format(pair, d, reject, null_LL, max_alt_LL,),
              file=output_file)

    output_file.close()
    print("--- {} seconds ---".format(round(time() - start_time, 3)))

if __name__ == '__main__':
    main()
