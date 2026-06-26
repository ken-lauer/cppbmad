from __future__ import annotations

import sys

import numpy.testing
import pytest
from pybmad import AllEncompassingStruct, TestSubStruct


@pytest.fixture
def aes():
    """Fixture to provide a fresh instance of AllEncompassingStruct."""
    return AllEncompassingStruct()


# -----------------------------------------------------------------------------
# 1. Real Attributes (rp/dp)
# -----------------------------------------------------------------------------
REAL_ATTRIBUTES = [
    # Real single precision
    pytest.param("real_rp_0d", "scalar", id="real_rp_0d"),
    pytest.param("real_rp_0d_ptr", "scalar", id="real_rp_0d_ptr"),
    pytest.param("real_rp_1d", "list", id="real_rp_1d"),
    pytest.param("real_rp_1d_alloc", "alloc", id="real_rp_1d_alloc"),
    pytest.param("real_rp_1d_ptr", "list", id="real_rp_1d_ptr"),
    pytest.param("real_rp_2d", "list", id="real_rp_2d"),
    pytest.param("real_rp_2d_alloc", "list", id="real_rp_2d_alloc"),
    pytest.param("real_rp_2d_ptr", "list", id="real_rp_2d_ptr"),
    pytest.param("real_rp_3d", "list", id="real_rp_3d"),
    pytest.param("real_rp_3d_alloc", "list", id="real_rp_3d_alloc"),
    pytest.param("real_rp_3d_ptr", "list", id="real_rp_3d_ptr"),
    # Real double precision
    pytest.param("real_dp_0d", "scalar", id="real_dp_0d"),
    pytest.param("real_dp_0d_ptr", "scalar", id="real_dp_0d_ptr"),
    pytest.param("real_dp_1d", "list", id="real_dp_1d"),
    pytest.param("real_dp_1d_alloc", "alloc", id="real_dp_1d_alloc"),
    pytest.param("real_dp_1d_ptr", "list", id="real_dp_1d_ptr"),
    pytest.param("real_dp_2d", "list", id="real_dp_2d"),
    pytest.param("real_dp_2d_alloc", "list", id="real_dp_2d_alloc"),
    pytest.param("real_dp_2d_ptr", "list", id="real_dp_2d_ptr"),
    pytest.param("real_dp_3d", "list", id="real_dp_3d"),
    pytest.param("real_dp_3d_alloc", "list", id="real_dp_3d_alloc"),
    pytest.param("real_dp_3d_ptr", "list", id="real_dp_3d_ptr"),
]


@pytest.mark.parametrize(("attr_name", "kind"), REAL_ATTRIBUTES)
def test_aes_real_attributes(aes: AllEncompassingStruct, attr_name, kind):
    """Test access to real (float) attributes."""
    val = getattr(aes, attr_name)

    print(kind, val, attr_name)  # noqa: T201
    if kind == "scalar":
        if val is not None and "0d" in attr_name:
            numpy.testing.assert_allclose(val, 0.0)
    elif attr_name.endswith(("_alloc", "_ptr")):
        assert len(val) == 0
    else:
        assert len(val)
        # if "1d" in attr_name:
        #     numpy.testing.assert_allclose(val[0], 1.0)
        # if "2d" in attr_name:
        #     numpy.testing.assert_allclose(val[0][0], 2.0)


# -----------------------------------------------------------------------------
# 2. Complex Attributes
# -----------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("attr", "expected"),
    [
        pytest.param("complex_dp_0d", 1.0 + 2j),
        pytest.param("complex_dp_1d", 2.0 + 2j),
        pytest.param("complex_dp_2d", 3.0 + 2j),
        pytest.param("complex_dp_3d", 4.0 + 2j),
    ],
)
def test_aes_complex(aes: AllEncompassingStruct, attr: str, expected):
    actual = getattr(aes, attr)
    numpy.testing.assert_allclose(actual=actual, desired=expected)


COMPLEX_ATTRIBUTES = [
    pytest.param("complex_dp_1d_alloc", "alloc", id="complex_dp_1d_alloc"),
    pytest.param("complex_dp_1d_ptr", "list", id="complex_dp_1d_ptr"),
    pytest.param("complex_dp_2d_alloc", "list", id="complex_dp_2d_alloc"),
    pytest.param("complex_dp_2d_ptr", "list", id="complex_dp_2d_ptr"),
    pytest.param("complex_dp_3d_alloc", "list", id="complex_dp_3d_alloc"),
    pytest.param("complex_dp_3d_ptr", "list", id="complex_dp_3d_ptr"),
]


@pytest.mark.parametrize(("attr_name", "kind"), COMPLEX_ATTRIBUTES)
def test_aes_complex_attributes(aes: AllEncompassingStruct, attr_name, kind):
    """Test access to complex attributes."""
    val = getattr(aes, attr_name)

    if kind == "scalar":
        if val is not None:
            assert isinstance(val, complex)
    else:
        assert len(val) >= 0


# -----------------------------------------------------------------------------
# 3. Integer Attributes (int and int8)
# -----------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("attr", "expected"),
    [
        pytest.param("int_0d", 0.0),
        pytest.param("int_1d", 1.0),
        pytest.param("int_2d", 2.0),
        pytest.param("int_3d", 3.0),
        pytest.param("int8_0d", 0.0),
        pytest.param("int8_1d", 1.0),
        pytest.param("int8_2d", 2.0),
        pytest.param("int8_3d", 3.0),
    ],
)
def test_aes_int(aes: AllEncompassingStruct, attr: str, expected):
    actual = getattr(aes, attr)
    numpy.testing.assert_allclose(actual=actual, desired=expected)


INT_ATTRIBUTES = [
    # Standard Int
    pytest.param("int_0d", "scalar", id="int_0d"),
    pytest.param("int_0d_ptr", "scalar", id="int_0d_ptr"),
    pytest.param("int_1d", "list", id="int_1d"),
    pytest.param("int_1d_alloc", "alloc", id="int_1d_alloc"),
    pytest.param("int_1d_ptr", "list", id="int_1d_ptr"),
    pytest.param("int_2d", "list", id="int_2d"),
    pytest.param("int_2d_alloc", "list", id="int_2d_alloc"),
    pytest.param("int_2d_ptr", "list", id="int_2d_ptr"),
    pytest.param("int_3d", "list", id="int_3d"),
    pytest.param("int_3d_alloc", "list", id="int_3d_alloc"),
    pytest.param("int_3d_ptr", "list", id="int_3d_ptr"),
    # Int8
    pytest.param("int8_0d", "scalar", id="int8_0d"),
    pytest.param("int8_0d_ptr", "scalar", id="int8_0d_ptr"),
    pytest.param("int8_1d", "list", id="int8_1d"),
    pytest.param("int8_1d_alloc", "alloc", id="int8_1d_alloc"),
    pytest.param("int8_1d_ptr", "list", id="int8_1d_ptr"),
    pytest.param("int8_2d", "list", id="int8_2d"),
    pytest.param("int8_2d_alloc", "list", id="int8_2d_alloc"),
    pytest.param("int8_2d_ptr", "list", id="int8_2d_ptr"),
    pytest.param("int8_3d", "list", id="int8_3d"),
    pytest.param("int8_3d_alloc", "list", id="int8_3d_alloc"),
    pytest.param("int8_3d_ptr", "list", id="int8_3d_ptr"),
]


@pytest.mark.parametrize(("attr_name", "kind"), INT_ATTRIBUTES)
def test_aes_int_attributes(aes: AllEncompassingStruct, attr_name, kind):
    """Test access to integer attributes."""
    val = getattr(aes, attr_name)

    if kind == "scalar":
        assert isinstance(val, (int, type(None)))
    elif kind == "alloc":
        assert len(val) == 0
        if hasattr(val, "resize"):
            val.resize(10)
            assert len(val) == 10
    else:
        assert len(val) >= 0


# -----------------------------------------------------------------------------
# 4. Logical (Boolean) Attributes
# -----------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("attr", "expected"),
    [
        pytest.param("logical_0d", True),
        pytest.param("logical_1d", True),
        pytest.param(
            "logical_2d",
            True,
            marks=pytest.mark.xfail(reason="no 2d logical array support"),
        ),
        pytest.param(
            "logical_3d",
            True,
            marks=pytest.mark.xfail(reason="no 3d logical array support"),
        ),
    ],
)
def test_aes_logical(aes: AllEncompassingStruct, attr: str, expected):
    actual = getattr(aes, attr)
    numpy.testing.assert_allclose(actual=actual, desired=expected)


LOGICAL_ATTRIBUTES = [
    pytest.param("logical_0d", "scalar", id="logical_0d"),
    pytest.param("logical_0d_ptr", "scalar", id="logical_0d_ptr"),
    pytest.param("logical_1d", "list", id="logical_1d"),
    pytest.param("logical_1d_ptr", "list", id="logical_1d_ptr"),
    pytest.param("logical_2d", "list", id="logical_2d"),
    pytest.param("logical_3d", "list", id="logical_3d"),
]


@pytest.mark.parametrize(("attr_name", "kind"), LOGICAL_ATTRIBUTES)
def test_aes_logical_attributes(aes: AllEncompassingStruct, attr_name, kind):
    """Test access to logical (bool) attributes."""
    val = getattr(aes, attr_name)

    if kind == "scalar":
        if val is not None:
            assert isinstance(val, bool)
    else:
        assert len(val) >= 0


# -----------------------------------------------------------------------------
# 5. Type (Struct) Attributes
# -----------------------------------------------------------------------------
TYPE_ATTRIBUTES = [
    pytest.param("type_0d", "scalar", id="type_0d"),
    pytest.param("type_0d_ptr", "ptr", id="type_0d_ptr"),
    pytest.param("type_1d", "list", id="type_1d"),
    pytest.param("type_1d_alloc", "alloc", id="type_1d_alloc"),
    pytest.param("type_1d_ptr", "ptr_list", id="type_1d_ptr"),
    pytest.param("type_2d", "list", id="type_2d"),
    pytest.param("type_2d_alloc", "list", id="type_2d_alloc"),
    pytest.param("type_2d_ptr", "list", id="type_2d_ptr"),
    pytest.param("type_3d", "list", id="type_3d"),
    pytest.param("type_3d_alloc", "list", id="type_3d_alloc"),
    pytest.param("type_3d_ptr", "list", id="type_3d_ptr"),
]


@pytest.mark.parametrize(("attr_name", "kind"), TYPE_ATTRIBUTES)
def test_aes_type_attributes(aes: AllEncompassingStruct, attr_name: str, kind: str):
    """Test access to struct-type attributes."""
    val = getattr(aes, attr_name)

    if kind == "ptr":
        # Per snippet: assert aes.type_0d_ptr is None
        assert val is None
    elif kind == "alloc":
        assert len(val) == 0
    elif kind == "list":
        if "alloc" in attr_name:
            assert len(val) == 0
        else:
            assert len(val) >= 0
    elif kind == "ptr_list":
        assert len(val) >= 0


def test_struct_alloc_access(aes: AllEncompassingStruct):
    aes = AllEncompassingStruct()
    assert len(aes.int8_1d_alloc) == 0
    aes.int8_1d_alloc.resize(10)
    assert len(aes.int8_1d_alloc) == 10


def test_derived_type_scalar_setter():
    aes = AllEncompassingStruct()
    sub = TestSubStruct()
    sub.sr.bbb = 77
    aes.type_0d = sub
    assert aes.type_0d.sr.bbb == 77


def test_derived_type_scalar_setter_is_deep_copy():
    aes = AllEncompassingStruct()
    sub = TestSubStruct()
    sub.sr.bbb = 77
    aes.type_0d = sub
    sub.sr.bbb = 99
    assert aes.type_0d.sr.bbb == 77  # copied, not aliased


def test_optional_scalar_pointer_reads_none_when_null():
    aes = AllEncompassingStruct()
    assert aes.real_rp_0d_ptr is None
    assert aes.int_0d_ptr is None
    assert aes.type_0d_ptr is None


def test_setting_null_optional_scalar_pointer_is_noop():
    # The Fortran setter only writes when the pointer is associated, so
    # assigning to a null pointer member is silently ignored.
    aes = AllEncompassingStruct()
    aes.real_rp_0d_ptr = 3.5
    assert aes.real_rp_0d_ptr is None


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", *sys.argv[1:], __file__]))
