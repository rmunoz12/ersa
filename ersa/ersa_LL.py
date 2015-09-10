"""Classes implementing Lp and Lr."""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from math import exp, log, factorial
from operator import itemgetter
from ersa.chisquare import LL_ratio_test, likelihood_ratio_CI


class Background:
    """
    Class Background
    ----------------
    Represents a calculator for the null hypothesis, Lp.

    Created with empirical parameters, t, theta, lambba_.

    Provides the likelihood function, L(n, s), and log-likelihood
    function, LL(n, s).


    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    h0 = Background(2.5, 3.12, 15)
    h0.LL(n, s)


    Notes
    -----
    Lp(n,s|t) = Np(n|t) * Sp(s|t)

        Np(n|t) = Poisson distribution with mean equal to sample mean
                  of the number of segments shared in the population (mu)
        Sp(s|t) = Product over i in s of Fp(i|t)
        Fp(i|t) = (e^(i-t)/theta)/theta

        n = number of shared segments
        s = list of shared segment lengths (in cM)
        t = minimum threshold for segment lengths (in cM)
        h = maximum threshold for segment lengths (in cM)
        lambda_ = the sample mean of the number of segments shared
                  in the population
        theta = mean shared segment length (in cm) in the population
                for all segments of size >t and <h
    """
    def __init__(self, t, theta, lambda_):
        self.t = t
        self.theta = theta
        self.lambda_ = lambda_

    def _Fp(self, i):
        assert i >= self.t
        prob = exp(-(i - self.t) / (self.theta - self.t)) / (self.theta - self.t)
        return log(prob)

    def _Sp(self, s):
        result = 0
        for i in s:
            result += self._Fp(i)
        return result

    def _Np(self, n):
        l_prob = n * log(self.lambda_) - self.lambda_ - log(factorial(n))
        return l_prob

    def LL(self, n, s):
        ret = 0
        ret += self._Np(n)
        ret += self._Sp(s)
        return ret


class Relation(Background):
    """
    Class Relation
    --------------
    Represents a calculator for the alternative hypothesis, Lr.

    Created with empirical parameters c, r, t, theta, and lambda_. The
    latter three empirical parameters are described in Class Background.

    Provides the a function returning the maximum log-likelihood over
    all possible np, MLL(n, s).  Create this class for each d evaluated.


    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    ha = Relation(22, 35.3, 2.5, 3.12, 14)
    ha.MLL(n, s)


    See also
    -------
    Class Background


    Notes
    -----
    Lr = La(na, sa | d, a, t) * Lp(np, sp | t)

        Lp = follows description of the null hypothesis (see
             Class Background)
        La(na, sa | d, a, t) = Na(n | d, a, t) * Sa(sa | d, t)
        Sa(s | d, t) = product over i in s of Fa(i | t)
        Fa(i | d, t) = (e^(-d(i-t)/100)) / (100/d0)
        Na(n | d, a, t) = <eq. 8 in Huff et. al. 2011>

        a = number of ancestors shared (set to 2 based on supp. mat.)
        d = combined number of generations separating the individuals
            from their ancestor(s)
        n = number of shared segments
        na = number of shared segments inherited from recent ancestors
        np = number of shared segments shared due to the population background
        s = list of shared segment lengths (in cM)
        sa and sp = mutually exclusive subsets of s, with sa equal to a list
                    of segment lengths inherited from recent ancestors and
                    sp equal to a list of segment lengths shared due to the
                    background

        r = expected number of recombination events per haploid genome
            per generation
        c = number of autosomes
        p(t) = probability that a shared segment is longer than t
    """
    def __init__(self, c, r, t, theta, lambda_):
        super(Relation, self).__init__(t, theta, lambda_)
        self.c = c
        self.r = r
        self.a = 2  # see Huff et al 2011 supplemental material

    def _Fa(self, i, d):
        assert i >= self.t
        l_prob = (-d * (i-self.t) / 100) - log(100 / d)
        return l_prob

    def _Sa(self, s, d):
        result = 0
        for i in s:
            result += self._Fa(i, d)
        return result

    def _p(self, d):
        return exp((-d * self.t) / 100)

    def _Na(self, n, d):
        lambda_ = (self.a * (self.r * d + self.c) * self._p(d)) / (2 ** (d - 1))
        l_prob = n * log(lambda_) - lambda_ - log(factorial(n))
        return l_prob

    def _LLr(self, np, na, s, d):
        """
        Compute Lr given d and np.

        Requires s to be pre-sorted from smallest to largest.
        """
        result = 0
        result += self._Np(np)
        result += self._Na(na, d)
        result += self._Sp(s[:np])
        result += self._Sa(s[np:], d)
        return result

    def MLL(self, n, s, d):
        """
        For a given d, return the maximum log-likelihood (MLL).
        Requires s to be sorted smallest to largest.
        """
        max_mll, max_np = None, None
        for np in range(n + 1):
            mll = self._LLr(np, n - np, s, d)
            if max_mll is None:
                max_mll, max_np = mll, np
            elif max_mll < mll:
                max_mll, max_np = mll, np
        return max_np, max_mll


def estimate_relation(pair, n, s, h0, ha, max_d):
    """
    Tests a pair of individuals for a relation and
    returns relevant parameters.

    Requires s to a pre-sorted list of shared segments,
    from smallest to largest.  h0 and ha must be Background
    Relation objects.
    """
    assert isinstance(h0, Background)
    assert isinstance(ha, Relation)
    for i in range(1, len(s)):
        assert s[i - 1] <= s[i]

    null_LL = h0.LL(n, s)

    alts = []
    for d in range(1, max_d + 1):
        alt_np, alt_MLL = ha.MLL(n, s, d)
        alts.append((d, alt_np, alt_MLL))
    max_alt = max(alts, key=itemgetter(2))
    d = max_alt[0] - 1
    max_LL = max_alt[2]

    reject = LL_ratio_test(max_LL, null_LL)
    lower_d, upper_d = 0, 0
    if reject:
        lower_d, upper_d = likelihood_ratio_CI(alts, max_LL)
        lower_d -= 1
        upper_d -= 1

    # subtract one from d based on a = 2
    return null_LL, max_LL, d, reject, lower_d, upper_d
