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
    with pytest.raises(FileNotFoundError):
        read_matchfile("")
    path = "tests/test_data/test_LL.match"
    s_list = read_matchfile(path)
    assert len(s_list) == 14


def test_get_pair_dict():
    path = "tests/test_data/test_LL.match"
    pair_dict = get_pair_dict(path, 2.5)
    assert pair_dict['TestA:TestB'][0] == 7

    pair_dict = get_pair_dict(path, 2.5, "TestA")
    with pytest.raises(KeyError):
        assert pair_dict['TestB:TestC'][0] == 9

