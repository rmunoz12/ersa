"""Likelihood ratio test"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from scipy import stats


def LL_ratio_test(LLr, LLn, df=2, alpha=0.05):
    """
    Perform a likelihood ratio test of LLr (alternative) and LLn (null)
    and return whether to accept (False) or reject (True) the null.

    Parameters
    ----------
    LLr: log-likelihood of alternative
    LLn: log-likelihood of null
    df: degrees of freedom for the ratio test
    alpha: confidence level
    """
    ratio = -2 * LLn + 2 * LLr 
    p = 1 - stats.chi2.cdf(ratio, df)
    return True if p < alpha else False


def likelihood_ratio_CI(alts, max_alt_LL, df=2, alpha=0.05):
    """
    For a given set of alternative models and the log-likelihood of the
    maximum model, returns the lower and upper estimates for d using
    a chi-square approximation for the likelihood ratio test with
    df degrees of freedom and an alpha confidence level.
    """
    lower_d, upper_d = None, None
    for alt in alts:
        if not LL_ratio_test(max_alt_LL, alt[2], df, alpha):
            d = alt[0]
            if not lower_d:
                lower_d, upper_d = d, d
            elif d < lower_d:
                lower_d = d
            elif d > upper_d:
                upper_d = d
    return lower_d, upper_d
