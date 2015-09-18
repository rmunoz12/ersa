"""Unit Tests for ersa/parser.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from ersa.dbhandler import *
from ersa.parser import get_pair_dict
from ersa.ersa_LL import Background, Relation, estimate_relation
import pytest


def get_test_data():
    MAX_D = 10
    t = 2.5                 # in cM
    h = 10                  # in cM
    theta = 3.197036753     # in cM
    lambda_ = 13.73         #
    r = 35.2548101          # ~for humans
    c = 22                  # human autosomes
    alpha = 0.05
    pair_dict = get_pair_dict('tests/test_data/test_LL.match', t)
    h0 = Background(t, theta, lambda_)
    ha = Relation(c, r, t, theta, lambda_)
    dob = (None, None)
    for pair, seg_list in pair_dict.items():
        s = [seg.length for seg in seg_list]
        n = len(s)
        est = estimate_relation(pair, dob, n, s, h0, ha, MAX_D, alpha)
        yield est, seg_list


def test_insert():
    db = DBhandler("sqlite:///", shared_pool=False)
    for e, s in get_test_data():
        db.insert(e, s)
    assert db.session.query(Result).count() == 2
    assert db.session.query(Likelihood).count() == 20
    assert db.session.query(Segment).count() == 10

    # ensure that duplicate pair combinations are
    # replaced with new data
    for e, s in get_test_data():
        db.insert(e, s)
    assert db.session.query(Result).count() == 2
    assert db.session.query(Likelihood).count() == 20
    assert db.session.query(Segment).count() == 10

