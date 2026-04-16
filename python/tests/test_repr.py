from __future__ import annotations

import re

import pytest

from pybmad import (
    AllEncompassingStruct,
    CoordStruct,
    RealAlloc1D,
    SplineStruct,
    TestSubStruct,
)


@pytest.fixture
def spline():
    s = SplineStruct()
    s.x0 = 1.0
    s.y0 = 2.0
    return s


@pytest.fixture
def coord():
    return CoordStruct()


@pytest.fixture
def aes():
    a = AllEncompassingStruct()
    a.real_dp_0d = 3.14
    a.int_0d = 42
    return a


class TestStructRepr:
    """Test __repr__ output of proxy struct objects."""

    def test_spline_repr_format(self, spline):
        r = repr(spline)
        assert r.startswith("SplineStruct(0x")
        assert r.endswith(")")
        assert "x0=" in r
        assert "y0=" in r

    def test_coord_repr_format(self, coord):
        r = repr(coord)
        assert r.startswith("CoordStruct(0x")
        assert "vec=" in r

    def test_aes_repr_format(self, aes):
        r = repr(aes)
        assert r.startswith("AllEncompassingStruct(0x")
        assert "real_dp_0d=" in r
        assert "int_0d=" in r

    def test_sub_struct_repr_format(self):
        sub = TestSubStruct()
        r = repr(sub)
        assert r.startswith("TestSubStruct(0x")
        assert "sr=" in r

    def test_repr_contains_hex_address(self, spline):
        r = repr(spline)
        # Should contain a hex address after the class name
        assert re.search(r"0x[0-9a-fA-F]+", r) is not None

    def test_repr_reflects_values(self):
        s = SplineStruct()
        s.x0 = 42.0
        r = repr(s)
        assert "42" in r

    def test_different_instances_different_addresses(self):
        s1 = SplineStruct()
        s2 = SplineStruct()
        r1 = repr(s1)
        r2 = repr(s2)
        # Extract hex addresses
        addr1 = re.search(r"0x([0-9a-fA-F]+)", r1).group(1)
        addr2 = re.search(r"0x([0-9a-fA-F]+)", r2).group(1)
        assert addr1 != addr2


class TestContainerRepr:
    """Test __repr__ output of array container objects."""

    def test_real_alloc_1d_repr(self):
        c = RealAlloc1D()
        c.resize(3)
        view = c.view()
        view[0] = 1.0
        view[1] = 2.0
        view[2] = 3.0
        r = repr(c)
        # Should be a meaningful representation
        assert len(r) > 0

    def test_real_alloc_1d_empty_repr(self):
        c = RealAlloc1D()
        r = repr(c)
        assert len(r) > 0
