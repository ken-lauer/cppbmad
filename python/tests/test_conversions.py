from __future__ import annotations

import pybmad
import pytest


def test_alloc_from_list():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    assert len(a) == 3
    assert list(a) == [1.0, 2.0, 3.0]


def test_int_alloc_from_list():
    a = pybmad.IntAlloc1D([1, 2, 3])
    assert list(a) == [1, 2, 3]


def test_alloc_from_view_is_deep_copy():
    src = pybmad.RealAlloc1D([1.0, 2.0])
    copy = pybmad.RealAlloc1D(src.view())
    assert list(copy) == [1.0, 2.0]
    copy[0] = 99.0
    assert src[0] == 1.0  # independent storage


def test_view_from_alloc_shares_memory():
    a = pybmad.RealAlloc1D([1.0, 2.0])
    view = pybmad.RealArray1D(a)
    assert list(view) == [1.0, 2.0]
    view[0] = 9.0
    assert a[0] == 9.0  # non-owning view aliases the allocator


def test_alloc_implicitly_converts_to_view_at_routine_boundary():
    # test_real_array's parameters are views; passing allocators directly
    # exercises the registered Alloc -> View implicit conversion.
    arr_inout = pybmad.RealAlloc1D([0.0, 0.0])
    res = pybmad.test.test_real_array(
        # TODO: type annotations say this is invalid, but it's implicitly
        # converted with nanobind by way of nb::implicitly_convertible
        arr_in=pybmad.RealAlloc1D([1.1, 2.2]),
        arr_inout=arr_inout,
    )
    assert list(res.arr_out) == pytest.approx([1.1, 2.2])
    assert list(arr_inout) == pytest.approx([1.0, 1.0])  # inout write reaches the allocator
