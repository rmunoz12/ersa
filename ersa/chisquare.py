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
    ratio = -2 * LLr + 2 * LLn  # TODO flip in order to correct notation
    p = 1 - stats.chi2.cdf(ratio, df)
    return True if p < alpha else False

#confidence interval for d given alpha
def likelihood_ratio_CI(alts, max_alt_LL, df=2, alpha=0.05):

    #threshold from inverse cdf
    # chi_thresh = stats.chi2.ppf(1-alpha,df)

    lower_d = 999
    upper_d = 0
    for alt in alts:
        # ratio = -2 * alt[2] + 2 * null_likelihood
        # if ratio < chi_thresh:
        #     if alt[0] < lower_d:
        #         lower_d = alt[0]
        #     elif alt[0] > upper_d:
        #         upper_d = alt[0]
        if not LL_ratio_test(alt[2], max_alt_LL):
            d = alt[0]
            if d < lower_d:
                lower_d = d
            if d > upper_d:
                upper_d = d
    return lower_d, upper_d
