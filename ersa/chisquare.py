"""Likelihood ratio test"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from scipy.stats import chi2


def LL_ratio_test(LLr, LLn, alpha=0.05):
    """
    Perform a likelihood ratio test of LLr (alternative) and LLn (null)
    using two degrees of freedom and return whether to accept or reject
    the null.

    Parameters
    ----------
    LLr : float
        Log-likelihood of alternative.

    LLn : float
        Log-likelihood of null.

    alpha : float
        Confidence level.

    Returns
    -------
    bool
        Accept (True) or reject (False) the null hypothesis.
    """
    p = LL_ratio_p(LLr, LLn, 2)
    return True if p < alpha else False


def LL_ratio_p(LLr, LLn, df=2):
    """
    Return p-value for likelihood ratio test using two degrees of
    freedom comparing LLr (alternative) and LLn (null).

    Parameters
    ----------
    LLr : float
        log-likelihood of alternative

    LLn : float
        log-likelihood of null

    Returns
    -------
    p : float
        p-value
    """
    ratio = -2 * LLn + 2 * LLr
    p = 1 - chi2.cdf(ratio, df)
    return p


def likelihood_ratio_CI(alts, max_alt_LL, alpha=0.05):
    """
    For a given set of alternative models and the log-likelihood of the
    maximum model, returns the lower and upper estimates for d using
    a chi-square approximation for the likelihood ratio test with
    df degrees of freedom and an alpha confidence level.

    Parameters
    ----------
    alts : (int, int , float)
        Tuple containing three elements corresponding to, in order,
        d value, number of segments in population (not used here),
        and the log-likelihood.

    max_alt_LL : float
        Log-likelihood of the maximum model

    alpha : float
        Confidence level.

    Returns
    -------
    (lower_d, upper_d) : (int, int)
        Lower and upper estimates for the degree of relation, d.
    """
    lower_d, upper_d = None, None
    for alt in alts:
        if not LL_ratio_test(max_alt_LL, alt[2], alpha):
            d = alt[0]
            if not lower_d:
                lower_d, upper_d = d, d
            elif d < lower_d:
                lower_d = d
            elif d > upper_d:
                upper_d = d
    return lower_d, upper_d
