from ersa.ersa_LL import *
import pytest
from math import log

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

