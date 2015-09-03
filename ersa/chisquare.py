from scipy import stats
import math


def likelihoodRatioTest(Lr, Ln, df=2):
    ratio = -2*math.log(Lr) + 2*math.log(Ln)
    return 1-stats.chi2.cdf(ratio, df)