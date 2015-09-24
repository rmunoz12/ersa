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
    path = "ersa/tests/test_data/test_LL.match"
    s_list = read_matchfile(path)
    len = 0
    for s in s_list:
        len += 1
    assert len == 14

    path = "ersa/tests/test_data/test_LL_haploscores.match"
    s_list = read_matchfile(path, haploscores=True)
    len = 0
    for s in s_list:
        len += 1
    assert len == 14

def test_get_pair_dict():
    path = "ersa/tests/test_data/test_LL.match"
    pair_dict = get_pair_dict(path, 2.5)
    assert len(pair_dict['TestA:TestB']) == 7

    pair_dict = get_pair_dict(path, 2.5, "TestA")
    with pytest.raises(KeyError):
        assert len(pair_dict['TestB:TestC']) == 9

