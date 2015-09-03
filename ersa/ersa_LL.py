# Classes implementing Lp, Lr, and MLr

from math import exp, log
from scipy.stats import poisson
import pylab as pl

"""
Class Background
--------------
Represents the null hypothesis, Lp.

Created with empirical parameters, t, theta, lambba_. Provides
the likelihood function, L(n, s), and log-likelihood function,
LL(n, s).

Example
-------
h0 = Background(2.5, 3.12, 15

Notes
-----
Lp(n,s|t) = Np(n|t) * Sp(s|t)

    Np(n|t) = Poisson distribution with mean equal to sample mean
              of the number of segments shared in the population (mu)
    Sp(s|t) = Product over i in s of Fp(i|t)
    Fp(i|t) = (e^(i-t)/theta)/theta

    t = 2.5 cM (empirical)
    h = 10 cM (empirical)
    lambda_ = the sample mean of the number of segments shared in the population (empirical)
    theta = mean shared segment length in the population for
            all segments of size >t and <ha, equal to 3.12 cM in
            Huff et. al. 2011 (empirical)
"""
class Background:
    def __init__(self, t, theta, lambda_):
        self.t = t
        self.theta = theta
        self.mu = lambda_

    def _Fp(self, i):
        assert(i > self.t)
        return (exp(-i - self.t)/self.theta) / self.theta

    def _Sp(self, s):
        result = 1
        for i in s:
            result *= self._Fp(i)
        return result

    def _Np(self, n):
        return poisson.pmf(n, self.mu)

    def L(self, n, s):
        return self._Np(n) * self._Sp(s)

    def LL(self, n, s):
        ret = 0
        ret += log(self._Np(n))
        ret += log(self._Sp(s))
        return ret


def main():
    t = 2.5
    theta = 3.12
    lambda_ = 15


if __name__ == '__main__':
    main()