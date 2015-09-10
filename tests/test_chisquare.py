from ersa.chisquare import *
from scipy.stats import chi2
from random import random, randint
from math import log

class Test_chisquare:

    def test_likelihood_ratio_test(self):
        for i in range(0,40):
            La = random()
            Ln = random()*La

            LLn = log(Ln)
            LLa = log(La)
            df = randint(0,30)
            alpha = random()

            ratio = -2 * LLn + 2 * LLa
            p = 1 - chi2.cdf(ratio, df)
            reject = True if p < alpha else False
            assert LL_ratio_test(LLa, LLn, df, alpha) == reject


    def test_likelihood_ratio_CI(self):
        max_d = 10
        for i in range(0,40):
            alts = []
            global_max_LL = 0
            for d in range(1, max_d + 1):
                max_LL = log(random())
                alts.append((d, 0, max_LL))
                if max_LL > global_max_LL:
                    global_max_LL = max_LL

            df = randint(0, 30)
            alpha = random()
            lower_d, upper_d = likelihood_ratio_CI(alts, global_max_LL, df, alpha)
            for alt in alts:
                if not LL_ratio_test(global_max_LL, alt[2], df, alpha):
                    assert alt[0] >= lower_d and alt[0] <= upper_d
