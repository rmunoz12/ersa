"""Unit Tests for ersa/parser.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


from ersa.parser import *
import pytest


def test_SharedSegment():
    with pytest.raises(AssertionError):
        SharedSegment(0)
    for n in [14, 16]:
        with pytest.raises(AssertionError):
            SharedSegment(n * [1])


def test_read_matchfile():
    pass

def test_get_pair_dict():
    pass
