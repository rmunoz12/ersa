"""Unit Tests for ersa_LL.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from ersa.ersa_LL import *
import pytest
from math import log
from scipy.stats import poisson

class TestBackground:
    t = 1
    theta = 2
    lambda_ = 3
    B = Background(t, theta, lambda_)

    def test_init(self):
        assert self.B.t == self.t
        assert self.B.theta == self.theta
        assert self.B.lambda_ == self.lambda_

    def test_Fp(self):
        with pytest.raises(AssertionError):
            self.B._Fp(self.t - 1)
        assert self.B._Fp(self.t) == log(1 / (self.theta - self.t))
        assert self.B._Fp(self.t + 1) == log(exp(-1 / (self.theta - self.t)) / (self.theta - self.t))

    def test_Sp(self):
        assert self.B._Sp([5]) == self.B._Fp(5)
        assert self.B._Sp([5, 10]) == self.B._Fp(5) + self.B._Fp(10)

    def test_Np(self):
        for n in range(10):
            Np_obs = self.B._Np(n)
            Np_exp = log(poisson.pmf(n, self.lambda_))
            assert abs(Np_exp - Np_obs) < (10 ** -8)

    def test_LL(self):
        s = [4, 5, 6, 7, 8, 10]
        n = len(s)
        assert self.B.LL(n, s) == self.B._Np(n) + self.B._Sp(s)


class TestRelation(Relation):
    def __init__(self):
        c = 1
        r = 2
        t = 3
        theta = 4
        lambda_ = 5
        super(TestRelation, self).__init__(c, r, t, theta, lambda_)

    def test_Fa(self):
        with pytest.raises(AssertionError):
            d = 5
            self._Fa(self.t - 1, d)
        assert self._Fa(5, 100) == -(5 - self.t)
        assert self._Fa(6, 1) == -(6 - self.t) / 100 - log(100)

    def test_Sa(self):
        d = 4
        assert self._Sa([5], d) == self._Fa(5, d)
        assert self._Sa([5, 10], d) == self._Fa(5, d) + self._Fa(10, d)

    def test_p(self):
        assert self._p(100) == exp(-self.t)

    def test_Na(self):
        for n in range(5):
            for d in range(1, 5):
                lambda_ = (self.a * (self.r * d + self.c) * self._p(d)) / (2 ** (d - 1))
                Na_obs = self._Na(n, d)
                Na_exp = log(poisson.pmf(n, lambda_))
                assert Na_obs == Na_exp

    def test_MLr(self):
        s = sorted([10, 8, 6, 4, 3])
        n = len(s)
        d = 3
        for np in range(n + 1):
            na = n - np
            mlr_obs = self._MLr(np, na, s, d)
            mlr_exp = self._Np(np) + self._Na(na, d) + self._Sp(s[:np]) + self._Sa(s[np:], d)
            assert mlr_obs == mlr_exp

    def test_MLL(self):
        pass  # TODO simple test cases / input files


def test_estimate_relation():
    pass  # TODO simple test cases / input files

