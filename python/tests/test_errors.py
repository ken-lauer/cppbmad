from __future__ import annotations

import pytest
from pybmad import (
    AllEncompassingStruct,
    CharacterAlloc1D,
    EleStructAlloc1D,
    IntAlloc1D,
    RealAlloc1D,
    TestSubStructAlloc1D,
)


class TestArrayIndexErrors:
    """Test that out-of-bounds access raises IndexError."""

    def test_real_alloc_1d_positive_oob(self):
        c = RealAlloc1D()
        c.resize(3)
        with pytest.raises(IndexError):
            c[3]

    def test_real_alloc_1d_negative_oob(self):
        c = RealAlloc1D()
        c.resize(3)
        with pytest.raises(IndexError):
            c[-4]

    def test_int_alloc_1d_positive_oob(self):
        c = IntAlloc1D()
        c.resize(3)
        with pytest.raises(IndexError):
            c[3]

    def test_int_alloc_1d_negative_oob(self):
        c = IntAlloc1D()
        c.resize(3)
        with pytest.raises(IndexError):
            c[-4]

    def test_character_alloc_1d_positive_oob(self):
        c = CharacterAlloc1D()
        c.resize(2, 10)  # n=2 strings, str_len=10
        with pytest.raises(IndexError):
            c[2]

    def test_character_alloc_1d_negative_oob(self):
        c = CharacterAlloc1D()
        c.resize(2, 10)
        with pytest.raises(IndexError):
            c[-3]


class TestStructArrayIndexErrors:
    """Test that out-of-bounds struct array access raises IndexError."""

    def test_ele_struct_alloc_positive_oob(self):
        arr = EleStructAlloc1D()
        arr.resize(2)
        with pytest.raises(IndexError):
            arr[2]

    def test_ele_struct_alloc_negative_oob(self):
        arr = EleStructAlloc1D()
        arr.resize(2)
        with pytest.raises(IndexError):
            arr[-3]

    def test_sub_struct_alloc_positive_oob(self):
        arr = TestSubStructAlloc1D()
        arr.resize(3)
        with pytest.raises(IndexError):
            arr[3]


class TestNegativeIndexing:
    """Verify that negative indexing works correctly (Python convention)."""

    def test_real_alloc_1d_negative_index(self):
        c = RealAlloc1D()
        c.resize(3)
        view = c.view()
        c[0] = 10.0
        view[1] = 20.0
        c[2] = 30.0
        assert c[-1] == pytest.approx(30.0)
        assert c[-2] == pytest.approx(20.0)
        assert c[-3] == pytest.approx(10.0)

    def test_int_alloc_1d_negative_index(self):
        c = IntAlloc1D()
        c.resize(3)
        view = c.view()
        view[0] = 1
        view[1] = 2
        c[2] = 3
        assert c[-1] == 3
        assert c[-3] == 1

    def test_ele_struct_alloc_negative_index(self):
        arr = EleStructAlloc1D()
        arr.resize(3)
        arr[0].name = "first"
        arr[2].name = "last"
        assert arr[-1].name == "last"
        assert arr[-3].name == "first"


class TestEmptyContainerAccess:
    """Test accessing elements of empty/unallocated containers."""

    def test_real_alloc_1d_empty(self):
        c = RealAlloc1D()
        assert len(c) == 0
        with pytest.raises(IndexError):
            c[0]

    def test_int_alloc_1d_empty(self):
        c = IntAlloc1D()
        assert len(c) == 0
        with pytest.raises(IndexError):
            c[0]

    def test_ele_struct_alloc_empty(self):
        arr = EleStructAlloc1D()
        assert len(arr) == 0
        with pytest.raises(IndexError):
            arr[0]


class TestSetitemErrors:
    """Test that out-of-bounds writes raise IndexError."""

    def test_real_alloc_1d_set_oob(self):
        c = RealAlloc1D()
        c.resize(3)
        with pytest.raises(IndexError):
            c[3] = 99.0

    def test_character_alloc_1d_set_oob(self):
        c = CharacterAlloc1D()
        c.resize(2, 10)
        with pytest.raises(IndexError):
            c[2] = "oops"


class TestTypeErrors:
    """Test that type mismatches raise appropriate errors."""

    def test_set_real_with_string(self):
        aes = AllEncompassingStruct()
        with pytest.raises(TypeError):
            aes.real_dp_0d = "not a number"

    def test_set_int_with_string(self):
        aes = AllEncompassingStruct()
        with pytest.raises(TypeError):
            aes.int_0d = "not an int"
