from __future__ import annotations

import numpy as np
import pytest
from pybmad import AllEncompassingStruct


@pytest.fixture  # scope="function"
def aes():
    return AllEncompassingStruct()


def test_array_1d_values_shape_dtype(aes: AllEncompassingStruct):
    arr = np.asarray(aes.real_dp_1d)
    assert arr.shape == (3,)
    assert arr.dtype == np.float64
    np.testing.assert_array_equal(arr, [1.0, 1.0, 1.0])


def test_array_1d_is_zero_copy(aes: AllEncompassingStruct):
    view = aes.real_dp_1d
    arr = np.asarray(view)
    arr[0] = 9.0
    assert view[0] == 9.0  # same view observes the write
    assert aes.real_dp_1d[0] == 9.0  # persisted into Fortran memory


def test_array_copy_is_independent(aes: AllEncompassingStruct):
    view = aes.real_dp_1d
    copy = np.array(view, copy=True)
    copy[0] = -5.0
    assert view[0] == 1.0


def test_array_2d_shape_and_fortran_order(aes: AllEncompassingStruct):
    arr = np.asarray(aes.real_dp_2d)
    assert arr.shape == (3, 4)
    assert arr.flags["F_CONTIGUOUS"]
    np.testing.assert_array_equal(arr, np.full((3, 4), 2.0))


def test_array_2d_is_zero_copy(aes: AllEncompassingStruct):
    arr = np.asarray(aes.real_dp_2d)
    arr[1, 2] = 7.0
    assert np.asarray(aes.real_dp_2d)[1, 2] == 7.0


def test_array_3d_is_zero_copy(aes: AllEncompassingStruct):
    arr = np.asarray(aes.real_dp_3d)
    arr[1, 2, 0] = 7.0
    assert np.asarray(aes.real_dp_3d)[1, 2, 0] == 7.0


def test_array_3d_shape_and_values(aes: AllEncompassingStruct):
    arr = np.asarray(aes.real_dp_3d)
    assert arr.shape == (3, 4, 5)
    np.testing.assert_array_equal(arr, np.full((3, 4, 5), 3.0))


def test_array_int_dtype(aes: AllEncompassingStruct):
    arr = np.asarray(aes.int_2d)
    assert arr.shape == (3, 4)
    assert arr.dtype == np.int32
    np.testing.assert_array_equal(arr, np.full((3, 4), 2))


def test_array_complex_dtype(aes: AllEncompassingStruct):
    arr = np.asarray(aes.complex_dp_1d)
    assert arr.dtype == np.complex128
    np.testing.assert_array_equal(arr, [2 + 2j, 2 + 2j, 2 + 2j])


def test_bool_array_has_no_zero_copy_protocol(aes: AllEncompassingStruct):
    # Fortran logical is 4 bytes with undocumented (?) false-y values, I think.
    # __array__ is intentionally not exposed for boolean arrays and falls back
    # to a listified result.
    assert not hasattr(aes.logical_1d, "__array__")
    assert hasattr(aes.real_dp_1d, "__array__")
    # TODO: stubgen bindings show __iter__ returning 'object'?
    assert list(aes.logical_1d) == [True, True, True]
    np.testing.assert_array_equal(np.asarray(aes.logical_1d), np.asarray([True, True, True]))
