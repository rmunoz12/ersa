"""Unit Tests for ersa_LL.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from ersa.ersa_LL import *
from ersa.parser import get_pair_dict
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


class TestRelation():
    c = 1
    r = 2
    t = 3
    theta = 4
    lambda_ = 5
    R = Relation(c, r, t, theta, lambda_)

    def test_Fa(self):
        with pytest.raises(AssertionError):
            d = 5
            self.R._Fa(self.t - 1, d)
        assert self.R._Fa(5, 100) == -(5 - self.R.t)
        assert self.R._Fa(6, 1) == -(6 - self.R.t) / 100 - log(100)

    def test_Sa(self):
        d = 4
        assert self.R._Sa([5], d) == self.R._Fa(5, d)
        assert self.R._Sa([5, 10], d) == self.R._Fa(5, d) + self.R._Fa(10, d)

    def test_p(self):
        assert self.R._p(100) == exp(-self.R.t)

    def test_Na(self):
        for n in range(5):
            for d in range(1, 5):
                lambda_ = (self.R.a * (self.r * d + self.c) * self.R._p(d)) / (2 ** (d - 1))
                Na_obs = self.R._Na(n, d)
                Na_exp = log(poisson.pmf(n, lambda_))
                assert Na_obs == Na_exp

    def test_MLr(self):
        s = sorted([10, 8, 6, 4, 3])
        n = len(s)
        d = 3
        for np in range(n + 1):
            na = n - np
            mlr_obs = self.R._MLr(np, na, s, d)
            mlr_exp = self.R._Np(np) + self.R._Na(na, d) + self.R._Sp(s[:np]) + self.R._Sa(s[np:], d)
            assert mlr_obs == mlr_exp

    def test_MLL(self):
        max_np, max_mll = self.R.MLL(1, [4], 2)
        assert max_np == 1
        assert max_mll == -9.099384755487144

        max_np, max_mll = self.R.MLL(5, [10, 5, 3], 8)
        assert max_np == 5
        assert max_mll == -10.949250206207347


def test_estimate_relation():
    """
    Uses hypothetical data to check:

    Pair            LLn         LLr         Reject  d   low_d   upp_d
    TestA:TestB     -78.0734    -30.5848    True    7   6       9
    TestB:TestC     -78.0734    -30.5848    True    7   7       9
    """

    MAX_D = 10

    t = 2.5                 # in cM
    h = 10                  # in cM
    theta = 3.197036753     # in cM
    lambda_ = 13.73         #
    r = 35.2548101          # ~for humans
    c = 22                  # human autosomes

    pair_dict = get_pair_dict('tests/test_data/test_LL.match', t, h)
    h0 = Background(t, theta, lambda_)
    ha = Relation(c, r, t, theta, lambda_)
    for pair, (n, s) in pair_dict.items():
        LLn, LLr, d, reject, low_d, upp_d = estimate_relation(pair, n, s, h0, ha, MAX_D)
        if pair == 'TestA:TestB':
            assert LLn == -78.07341166722982
            assert LLr == -30.584814485465674
            assert reject
            assert d == 7
            assert low_d == 6
            assert upp_d == 9
        if pair == 'TestB:TestC':
            assert LLn == -31.873042343537115
            assert LLr == -15.445177211112602
            assert reject
            assert d == 9
            assert low_d == 7
            assert upp_d == 9
