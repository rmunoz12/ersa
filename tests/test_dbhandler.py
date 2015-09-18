"""Unit Tests for ersa/parser.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from ersa.dbhandler import *
import pytest


def test_insert():
    alts = []
    for d in range(9):
        alts.append((d, 3, -10.123 + d))
    s = [2.5 + x for x in range(1, 10)]  # Todo change to Shared Segments
    e = Estimate("TestA:TestB", (1950, 2010), 2, True, -20.3, -8.23, 1, 4, alts, s)

    insert("sqlite://", e, s)

