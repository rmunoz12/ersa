# Classes implementing Lp, Lr, and MLr

from math import exp, log
from scipy.stats import poisson
import pylab as pl

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
class Background:
    def __init__(self, t, theta, lambda_):
        self.t = t
        self.theta = theta
        self.lambda_ = lambda_

    def _Fp(self, i):
        assert(i > self.t)
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


"""
Class Relation
--------------
Represents the alternative hypothesis, Lr.

Created with empirical parameters c, r, t, theta, and lambda_. The
latter three parameters are described in Class Background.

Provides the a function returning the maximum log-likelihood over
d and np, MLL(n, s)


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
class Relation(Background):
    def __init__(self, c, r, t, theta, lambda_):
        super(Relation, self).__init__(t, theta, lambda_)
        self.c = c
        self.r = r
        self.a = 2  # see Huff et al 2011 supplemental material

    def _Fa(self, i, d):
        assert(i > self.t)
        return exp(-d * (i - self.t) / 100) / (100 / d)

    def _Sa(self, s, d):
        result = 1
        for i in s:
            result *= self._Fa(i, d)
        return result

    def _p(self, d):
        return exp((-d * self.t) / 100)

    def _Na(self, n, d):
        lambda_ = (-self.a * (self.r * d + self.c) * self._p(d)) / (2 ** (d - 1))
        return poisson.pmf(n, lambda_)

    def _La(self, na, sa, d):
        return self._Na(na, d) * self._Sa(sa, d)  # changed n on RHS of eq. 5 to na

    # s must be sorted smallest to largest
    def _MLr(self, np, na, s, d):
        result = 1
        result *= self._Np(np)
        result *= self._Na(na, d)
        result *= self._Sp(s[:np + 1])
        result *= self._Sa(s[np+1:], d)
        return result

    def MLL(self, n, s):
        mll_dict = {}
        s_sorted = sorted(s)
        for d in range(40):
            mll = 0
            for np in range(n + 1):
                mlr = self._MLr(np, n - np, s_sorted, d)
                mll += log(mlr)
            mll_dict[d] = mll
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
