"""Likelihood ratio test"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from scipy import stats

def test_LL_ratio(LLr, LLn, df=2, alpha=0.05):
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
    ratio = -2 * LLr + 2 * LLn
    p = 1 - stats.chi2.cdf(ratio, df)
    return True if p < alpha else False

