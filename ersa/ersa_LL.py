"""Classes implementing Lp and Lr."""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from math import exp, log
from scipy.stats import poisson
import pylab as pl


class Background:
    """
    Class Background
    ----------------
    Represents the null hypothesis, Lp.

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
        assert i > self.t
        return (exp(-i - self.t)/self.theta) / self.theta

    def _Sp(self, s):
        result = 1
        for i in s:
            result *= self._Fp(i)
        return result

    def _Np(self, n):
        return poisson.pmf(n, self.lambda_)

    def L(self, n, s):
        return self._Np(n) * self._Sp(s)

    def LL(self, n, s):
        ret = 0
        ret += log(self._Np(n))
        ret += log(self._Sp(s))
        return ret



class Relation(Background):
    """
    Class Relation
    --------------
    Represents the alternative hypothesis, Lr.

    Created with empirical parameters c, r, t, theta, and lambda_, and the
    chosen d value. The latter three empirical parameters are described in
    Class Background.

    Provides the a function returning the maximum log-likelihood over
    all possible np, MLL(n, s).  Create this class for each d evaluated.


    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    ha = Relation(22, 35.3, 2.5, 3.12, 14)
    ha.MLL(n, s)


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
    def __init__(self, c, r, d, t, theta, lambda_):
        super(Relation, self).__init__(t, theta, lambda_)
        self.c = c
        self.r = r
        self.d = d
        self.a = 2  # see Huff et al 2011 supplemental material

    def _Fa(self, i):
        assert i > self.t
        return exp(-self.d * (i - self.t) / 100) / (100 / self.d)

    def _Sa(self, s):
        result = 1
        for i in s:
            result *= self._Fa(i)
        return result

    def _p(self):
        return exp((-self.d * self.t) / 100)

    def _Na(self, n):
        lambda_ = (-self.a * (self.r * self.d + self.c) * self._p()) / (2 ** (self.d - 1))
        return poisson.pmf(n, lambda_)

    def _La(self, na, sa):
        return self._Na(na) * self._Sa(sa)  # changed n on RHS of eq. 5 to na

    # s must be sorted smallest to largest
    def _MLr(self, np, na, s):
        result = 1
        result *= self._Np(np)
        result *= self._Na(na)
        result *= self._Sp(s[:np + 1])
        result *= self._Sa(s[np+1:])
        return result

    def MLL(self, n, s):
        mll_dict = {}
        s_sorted = sorted(s)
        for np in range(n + 1):
            mlr = self._MLr(np, n - np, s_sorted)
            mll = log(mlr)
            mll_dict[np] = mll
        max_key = max(mll_dict.keys(), key=lambda k: mll_dict[k])
        return max_key, mll_dict[max_key]


def main():
    t = 2.5       # in cM
    h = 10        # in cM
    theta = 3.12  # in cM
    lambda_ = 14  # estimated from Fig. 2
    r = 35.3      # ~for humans
    c = 22        # human autosomes
    pass


if __name__ == '__main__':
    main()
