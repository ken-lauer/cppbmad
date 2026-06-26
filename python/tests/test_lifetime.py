"""Lifetime / keep-alive coverage for containers, views, numpy, and nested proxies.

The pattern is: build a parent, take a child that should keep the parent alive,
drop the parent (+ ``gc.collect()``), then assert the child still yields correct
data. A broken keep-alive would surface as wrong values or a crash.
"""

from __future__ import annotations

import gc
import sys

import numpy as np
import pybmad
import pytest
from pybmad import AllEncompassingStruct, TestSubStruct

# --------------------------------------------------------------------------
# Refcount increments (safe: no use-after-free even if keep-alive is broken)
# --------------------------------------------------------------------------


def test_view_increments_allocator_refcount():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    rc_before = sys.getrefcount(a)
    v = a.view()
    assert sys.getrefcount(a) > rc_before
    del v


def test_iter_increments_allocator_refcount():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    rc_before = sys.getrefcount(a)
    it = iter(a)
    assert sys.getrefcount(a) > rc_before
    del it


def test_numpy_array_holds_owner_reference():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    rc_before = sys.getrefcount(a)
    arr = np.asarray(a)
    assert sys.getrefcount(a) > rc_before  # ndarray owns a reference to `a`
    del arr


def test_slice_view_increments_allocator_refcount():
    a = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    rc_before = sys.getrefcount(a)
    s = a[1:4]
    assert sys.getrefcount(a) > rc_before
    del s


# --------------------------------------------------------------------------
# Survival after the parent is dropped
# --------------------------------------------------------------------------


def test_view_survives_allocator_deletion():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    v = a.view()
    del a
    gc.collect()
    assert list(v) == [1.0, 2.0, 3.0]


def test_numpy_array_survives_allocator_deletion():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    arr = np.asarray(a)
    del a
    gc.collect()
    assert arr.tolist() == [1.0, 2.0, 3.0]
    arr[0] = 42.0  # writing into still-valid backing memory
    assert arr[0] == 42.0


def test_slice_view_survives_allocator_deletion():
    a = pybmad.RealAlloc1D([0.0, 1.0, 2.0, 3.0, 4.0])
    s = a[1:4]
    del a
    gc.collect()
    assert s.tolist() == [1.0, 2.0, 3.0]


def test_explicit_view_from_alloc_survives_allocator_deletion():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    v = pybmad.RealArray1D(a)  # keep_alive<1, 2> ties `a` to the view
    del a
    gc.collect()
    assert list(v) == [1.0, 2.0, 3.0]


def test_iterator_survives_allocator_deletion():
    a = pybmad.RealAlloc1D([1.0, 2.0, 3.0])
    it = iter(a)
    del a
    gc.collect()
    assert list(it) == [1.0, 2.0, 3.0]


# --------------------------------------------------------------------------
# Nested struct-member proxy chains
# --------------------------------------------------------------------------


def test_member_proxy_survives_parent_deletion():
    aes = AllEncompassingStruct()
    aes.type_0d.sr.bbb = 5
    sub = aes.type_0d  # proxy into aes's memory; must keep aes alive
    del aes
    gc.collect()
    assert sub.sr.bbb == 5


def test_deep_member_proxy_survives_parent_deletion():
    aes = AllEncompassingStruct()
    aes.type_0d.sr.bbb = 7
    sr = aes.type_0d.sr  # two levels deep
    del aes
    gc.collect()
    assert sr.bbb == 7


def test_struct_array_element_survives_parent_deletion():
    aes = AllEncompassingStruct()
    aes.type_1d[0].sr.bbb = 9
    elem = aes.type_1d[0]  # element proxy -> array -> aes chain
    del aes
    gc.collect()
    assert elem.sr.bbb == 9


def test_type_alloc_element_survives_container_deletion():
    arr = TestSubStruct.new_array1d(3)
    arr[0].sr.bbb = 3
    elem = arr[0]
    del arr
    gc.collect()
    assert elem.sr.bbb == 3


def test_repeated_access_does_not_corrupt():
    aes = AllEncompassingStruct()
    aes.type_0d.sr.bbb = 11
    # churn through many transient proxies; the underlying value must be stable
    for _ in range(1000):
        assert aes.type_0d.sr.bbb == 11
    gc.collect()
    assert aes.type_0d.sr.bbb == 11


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", *sys.argv[1:], __file__]))
