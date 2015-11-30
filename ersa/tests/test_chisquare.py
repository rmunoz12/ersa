"""Unit Tests for ersa/chisquare.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from math import log
from random import random
from ersa.chisquare import *

class Test_chisquare:

    num_iter = 40
    max_d = 10

    def test_likelihood_ratio_test(self):
        """
        For num_iter iterations, generate random values for La and Ln with La >= Ln, and perform the likelihood ratio
        test with random df in range [0, max_df] and alpha [0, 1]. Confirm that accepted/rejected tests match the
        output from LL_ratio_test()
        """
        for i in range(0, self.num_iter):
            La = random()
            Ln = random()*La

            LLn = log(Ln)
            LLa = log(La)
            alpha = random()

            ratio = -2 * LLn + 2 * LLa
            p = 1 - (1 - exp(-ratio / 2))
            reject = True if p < alpha else False
            assert LL_ratio_test(LLa, LLn, alpha) == reject


    def test_likelihood_ratio_CI(self):
        """
        For num_iter iterations, generate a set of maximum log likelihoods for each d in range [0, max_d], corresponding
        to one pair. Confirm that d values in the range returned by likelihood_ratio_CI reject the null hypothesis.
        """
        for i in range(0, self.num_iter):
            alts = []
            global_max_LL = 0
            for d in range(1, self.max_d + 1):
                max_LL = log(random())
                alts.append((d, 0, max_LL))
                if max_LL > global_max_LL:
                    global_max_LL = max_LL

            alpha = random()
            lower_d, upper_d = likelihood_ratio_CI(alts, global_max_LL, alpha)
            for alt in alts:
                if not LL_ratio_test(global_max_LL, alt[2], alpha):
                    assert alt[0] >= lower_d and alt[0] <= upper_d
