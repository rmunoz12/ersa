"""Unit Tests for ersa/mask.py"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

import pytest
from ersa.mask import *
from ersa.parser import SharedSegment


def test_mask_input_segs():
    t = 2.5
    b = 1 * 10 ** 6
    # Mask at chrome 9
    ml = 38293483  # start: 38,293,483
    mh = 72605261  # end:   72,605,261
    mask_cm = 8.15

    l = 6.29
    s_params = ["0", "user1", "0", "user2",
                "9", str(ml), str(mh),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0

    x = ml - b
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), "72605261",
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0

    y = mh + b
    l = mask_cm + 2
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0  # adjust length below t = 2.5

    l = mask_cm + 2.5
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 1
    assert res[0].length == l - mask_cm

    x = ml + 50
    y = mh + 50
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0

    x = ml - 50
    y = mh - 50
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0

    x = ml + 50
    y = mh + b + 1
    l = 6.29
    ratio = (y - mh) / (y - x)
    adj_len = round(l * ratio, 2)
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0  # adjust length below t = 2.5

    x = ml + 50
    y = mh + 20 * b
    l = 10.29
    ratio = (y - mh) / (y - x)
    adj_len = round(l * ratio, 2)
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 1
    assert res[0].length == adj_len

    x = ml - b
    y = mh - 50
    l = 6.29
    ratio = (ml - x) / (y - x)
    adj_len = round(l * ratio, 2)
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 0  # adjust length below t = 2.5

    x = ml - 20 * b
    y = mh - 50
    l = 10.29
    ratio = (ml - x) / (y - x)
    adj_len = round(l * ratio, 2)
    s_params = ["0", "user1", "0", "user2",
                "9", str(x), str(y),
                "abc1", "abc2", "100",
                str(l), "cM",
                "0", "0", "0"]
    s = SharedSegment(s_params)
    res = mask_input_segs([s], t)
    assert len(res) == 1
    assert res[0].length == adj_len


def test_total_masked():
    m = total_masked()
    m = round(m, 3)
    assert m == 119.92
