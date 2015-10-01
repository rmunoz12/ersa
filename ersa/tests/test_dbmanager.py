"""Unit Tests for ersa/parser.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from ersa.dbmanager import *
from ersa.parser import get_pair_dict
from ersa.ersa_LL import Background, Relation, estimate_relation


def get_test_data():
    MAX_D = 10
    t = 2.5                 # in cM
    h = 10                  # in cM
    theta = 3.197036753     # in cM
    lambda_ = 13.73         #
    r = 35.2548101          # ~for humans
    c = 22                  # human autosomes
    alpha = 0.05
    pair_dict = get_pair_dict('ersa/tests/test_data/test_LL.match', t)
    h0 = Background(t, theta, lambda_)
    ha = Relation(c, r, t, theta, lambda_)
    dob = (None, None)
    for pair, seg_list in pair_dict.items():
        s = [seg.length for seg in seg_list]
        n = len(s)
        est = estimate_relation(pair, dob, n, s, h0, ha, MAX_D, alpha)
        yield est, seg_list


def test_insert():
    with DbManager("sqlite:///", shared_pool=False) as db:
        ests, segs = [], []
        for e, s in get_test_data():
            ests.append(e)
            segs.append(s)

        db.insert(ests, segs)

        q = select([Result.__table__])
        res = db.conn.execute(q)
        count = 0
        for row in res:
            count += 1
        assert count == 2

        q = select([Segment.__table__])
        res = db.conn.execute(q)
        count = 0
        for row in res:
            count += 1
        assert count == 10


def test_soft_delete():
    with DbManager("sqlite:///", shared_pool=False) as db:
        pairs = ['TestA:TestB', 'TestC:TestB']
        ests, segs = [], []
        for e, s in get_test_data():
            ests.append(e)
            segs.append(s)

        db.insert(ests, segs)
        assert db.soft_delete(pairs) == 2

        db.insert(ests, segs)
        assert db.soft_delete(pairs) == 2

def test_delete():
    with DbManager("sqlite:///", shared_pool=False) as db:
        pairs = ['TestA:TestB', 'TestC:TestB']
        ests, segs = [], []
        for e, s in get_test_data():
            ests.append(e)
            segs.append(s)

        db.insert(ests, segs)
        db.soft_delete(pairs)
        n_deleted = db.delete()
        assert n_deleted['r'] == 2
        assert n_deleted['s'] == 10

        db.insert(ests, segs)
        db.insert(ests, segs)
        n_deleted = db.delete()
        assert n_deleted['r'] == 2
        assert n_deleted['s'] == 10

        # 2 results in db prior to this block
        db.insert(ests, segs)  # 4 results
        db.insert(ests, segs)  # 6 results
        db.insert(ests, segs)  # 8 results
        n_deleted = db.delete()
        assert n_deleted['r'] == 6
        assert n_deleted['s'] == 30
