from __future__ import annotations

import numpy as np
import pybmad
import pytest


def test_slice_read_returns_numpy_array():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    sliced = c[1:4]
    assert isinstance(sliced, np.ndarray)
    np.testing.assert_array_equal(sliced, [1.0, 2.0, 3.0])


def test_slice_read_is_zero_copy():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    sliced = c[1:4]
    sliced[0] = 99.0
    assert c[1] == 99.0


def test_slice_step():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(c[::2], [0.0, 2.0, 4.0])


def test_slice_negative_bounds():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(c[-2:], [3.0, 4.0])


def test_slice_assignment():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    c[1:4] = [10.0, 20.0, 30.0]
    assert list(c) == [0.0, 10.0, 20.0, 30.0, 4.0]


def test_slice_assignment_length_mismatch_raises():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    with pytest.raises(ValueError, match="different sizes"):
        c[1:4] = [1.0, 2.0]


def test_view_slice_is_zero_copy():
    c = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    view = c.view()
    sliced = view[1:3]
    assert isinstance(sliced, np.ndarray)
    sliced[0] = -7.0
    assert c[1] == -7.0


def test_bool_slice_returns_list():
    b = pybmad.BoolAlloc1D()
    b.resize(4)
    for i, val in enumerate([True, False, True, False]):
        b[i] = val
    sliced = b[1:3]
    assert isinstance(sliced, list)
    assert sliced == [False, True]
