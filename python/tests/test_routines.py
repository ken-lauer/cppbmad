from __future__ import annotations

import sys

import pybmad
import pytest


def container_create(container_cls, values: list):
    arr = container_cls()
    n = len(values)
    if values:
        arr.resize(0, n)
        for i, val in enumerate(values):
            arr[i] = val
    return arr


def get_optional_ids(*val):
    """Generates clean names for tests based on optional arg presence."""
    print(val)  # noqa: T201
    in_stat = "InOpt" if val[0] else "NoInOpt"
    inout_stat = "InOutOpt" if val[1] else "NoInOutOpt"
    return f"{in_stat}-{inout_stat}"


optionals = pytest.mark.parametrize(
    ("use_in_opt", "use_inout_opt"),
    [
        pytest.param(False, False, id="no-optional"),
        pytest.param(True, False, id="InOpt"),
        pytest.param(False, True, id="InOutOpt"),
        pytest.param(True, True, id="InOpt-InOutOpt"),
    ],
)


@optionals
def test_integer_scalar(use_in_opt: bool, use_inout_opt: bool):
    """Test integer scalar permutations."""
    val_in = 10
    val_inout = 20

    val_in_opt = 5 if use_in_opt else None
    val_inout_opt = 50 if use_inout_opt else None

    res = pybmad.test_integer_scalar(
        val_in=val_in, val_inout=val_inout, val_in_opt=val_in_opt, val_inout_opt=val_inout_opt
    )

    expected_out = 10
    if use_in_opt:
        expected_out += 5
    assert res.val_out == expected_out

    # TODO: this needs special handling on the python end (immutable)
    # 2. InOut Modification (+1)
    # expected_inout_final = 21
    # assert val_inout == expected_inout_final

    # 3. Status
    expected_status = [1 if use_in_opt else 0, 1 if use_inout_opt else 0]
    # list() conversion required for c_int array
    assert list(res.opt_status) == expected_status

    # TODO: this needs special handling on the python end (immutable)
    # 4. Optional InOut (+1)
    # if use_inout_opt:
    #     assert val_inout_opt == 51


@optionals
def test_integer_array(use_in_opt: bool, use_inout_opt: bool):
    """Test integer array (Alloc1D) wrappers."""
    input_data = [1, 2, 3]

    arr_in = container_create(pybmad.IntAlloc1D, input_data)
    arr_inout = container_create(pybmad.IntAlloc1D, [10, 20, 30])

    kw_args = {"arr_in": arr_in, "arr_inout": arr_inout}

    if use_in_opt:
        arr_in_opt = container_create(pybmad.IntAlloc1D, [5, 5, 5])
        kw_args["arr_in_opt"] = arr_in_opt

    if use_inout_opt:
        arr_inout_opt = container_create(pybmad.IntAlloc1D, [100, 100, 100])
        kw_args["arr_inout_opt"] = arr_inout_opt

    res = pybmad.test_integer_array(**kw_args)

    expected_out_vals = [x + (5 if use_in_opt else 0) for x in input_data]
    assert list(res.arr_out) == expected_out_vals

    # arr_inout modified in place due to C++ reference wrapper
    assert list(arr_inout) == [11, 21, 31]

    expected_stat = [1 if use_in_opt else 0, 1 if use_inout_opt else 0]
    assert list(res.opt_status) == expected_stat

    if use_inout_opt:
        assert list(kw_args["arr_inout_opt"]) == [101, 101, 101]


@optionals
def test_integer8_array(use_in_opt: bool, use_inout_opt: bool):
    """Test 64-bit integer allocations."""
    input_data = [2**33, 2**34]

    arr_in = container_create(pybmad.Int8Alloc1D, input_data)
    arr_inout = container_create(pybmad.Int8Alloc1D, [10, 20])

    kw_args = {"arr_in": arr_in, "arr_inout": arr_inout}

    if use_in_opt:
        kw_args["arr_in_opt"] = container_create(pybmad.Int8Alloc1D, [1, 1])

    if use_inout_opt:
        kw_args["arr_inout_opt"] = container_create(pybmad.Int8Alloc1D, [5, 5])

    res = pybmad.test_integer8_array(**kw_args)

    expected_add = 1 if use_in_opt else 0
    assert list(res.arr_out) == [x + expected_add for x in input_data]

    assert list(arr_inout) == [11, 21]


@optionals
def test_real_scalar(use_in_opt: bool, use_inout_opt: bool):
    val_in = 1.5
    val_inout = 2.5

    kw = {
        "val_in": val_in,
        "val_inout": val_inout,
        "val_in_opt": 0.5 if use_in_opt else None,
        "val_inout_opt": 9.5 if use_inout_opt else None,
    }

    res = pybmad.test_real_scalar(**kw)

    expected_out = 1.5 + (0.5 if use_in_opt else 0.0)
    assert res.val_out == pytest.approx(expected_out)
    # TODO: immutable pybmad
    # assert val_inout == pytest.approx(3.5)  # 2.5 + 1.0


@optionals
def test_real_array(use_in_opt: bool, use_inout_opt: bool):
    """Test RealAlloc1D."""
    arr_in = container_create(pybmad.RealAlloc1D, [1.1, 2.2])
    arr_inout = container_create(pybmad.RealAlloc1D, [0.0, 0.0])

    kw = {"arr_in": arr_in, "arr_inout": arr_inout}

    if use_in_opt:
        kw["arr_in_opt"] = container_create(pybmad.RealAlloc1D, [10.0, 10.0])
    if use_inout_opt:
        kw["arr_inout_opt"] = container_create(pybmad.RealAlloc1D, [-1.0, -1.0])

    res = pybmad.test_real_array(**kw)

    out_vals = list(res.arr_out)
    expected_add = 10.0 if use_in_opt else 0.0
    assert out_vals == pytest.approx([1.1 + expected_add, 2.2 + expected_add])

    inout_vals = list(arr_inout)
    assert inout_vals == pytest.approx([1.0, 1.0])  # 0 + 1.0

    if use_inout_opt:
        opt_inout_vals = list(kw["arr_inout_opt"])
        assert opt_inout_vals == pytest.approx([0.0, 0.0])  # -1.0 + 1.0


@optionals
def test_real16_scalar(use_in_opt: bool, use_inout_opt: bool):
    """Testing extended precision scalar mapping."""
    res = pybmad.test_real16_scalar(
        val_in=100.0,
        val_inout=200.0,
        val_in_opt=50.0 if use_in_opt else None,
        val_inout_opt=300.0 if use_inout_opt else None,
    )

    expect_out = 100.0 + (50.0 if use_in_opt else 0.0)
    assert res.val_out == pytest.approx(expect_out)
    # TODO immutable
    # assert val_inout == pytest.approx(201.0)


@optionals
@pytest.mark.xfail(reason="float128 isn't supported on clang/macos", condition=sys.platform == "darwin")
def test_real16_array(use_in_opt: bool, use_inout_opt: bool):
    """Note: mapped to Real16Alloc1D in types description provided."""
    arr_in = container_create(pybmad.Real16Alloc1D, [1.0])
    arr_inout = container_create(pybmad.Real16Alloc1D, [2.0])

    kw = {"arr_in": arr_in, "arr_inout": arr_inout}
    if use_in_opt:
        kw["arr_in_opt"] = container_create(pybmad.Real16Alloc1D, [0.5])
    if use_inout_opt:
        kw["arr_inout_opt"] = container_create(pybmad.Real16Alloc1D, [3.0])

    res = pybmad.test_real16_array(**kw)

    out = res.arr_out[0]
    expected = 1.0 + (0.5 if use_in_opt else 0.0)
    assert out == pytest.approx(expected)

    inout = arr_inout[0]
    assert inout == pytest.approx(3.0)  # 2.0 + 1.0


@optionals
def test_complex_array(use_in_opt: bool, use_inout_opt: bool):
    val_in = [complex(1, 2), complex(3, 4)]
    arr_in = container_create(pybmad.ComplexAlloc1D, val_in)
    arr_inout = container_create(pybmad.ComplexAlloc1D, [complex(0, 0), complex(0, 0)])
    arr_inout_opt = container_create(pybmad.ComplexAlloc1D, [complex(0, 0), complex(0, 0)])

    kw = {"arr_in": arr_in, "arr_inout": arr_inout}
    if use_in_opt:
        kw["arr_in_opt"] = container_create(pybmad.ComplexAlloc1D, [complex(1, 1), complex(1, 1)])
    if use_inout_opt:
        kw["arr_inout_opt"] = arr_inout_opt

    res = pybmad.test_complex_array(**kw)

    # Fortran: arr_out + (if opt: arr_in_opt)
    add_val = complex(1, 1) if use_in_opt else complex(0, 0)

    assert res.arr_out[0] == pytest.approx(val_in[0] + add_val)
    assert res.arr_out[1] == pytest.approx(val_in[1] + add_val)

    # InOut (added (1.0, 1.0))
    assert arr_inout[0] == pytest.approx(complex(1, 1))

    if use_inout_opt:
        assert all(v == pytest.approx(complex(1, 1)) for v in arr_inout_opt)


@optionals
def test_logical_scalar(use_in_opt: bool, use_inout_opt: bool):
    # Fortran: Out = In. If InOpt present: Out = (Out .eqv. InOpt)
    # Fortran: InOut = NOT InOut.

    kw = {
        "val_in": True,
        "val_inout": True,
        "val_in_opt": False if use_in_opt else None,
        "val_inout_opt": False if use_inout_opt else None,
    }

    res = pybmad.test_logical_scalar(**kw)

    expected_out = True
    if use_in_opt:
        # T .eqv. F -> F
        expected_out = False

    assert res.val_out == expected_out
    # TODO immutable
    # assert val_inout is False  # NOT True


@optionals
def test_logical_array(use_in_opt: bool, use_inout_opt: bool):
    # Logic: Out = In; InOut = NOT InOut
    arr_in = container_create(pybmad.BoolAlloc1D, [True, False])
    arr_inout = container_create(pybmad.BoolAlloc1D, [True, False])
    arr_inout_opt = container_create(pybmad.BoolAlloc1D, [False, True])

    kw = {"arr_in": arr_in, "arr_inout": arr_inout}
    if use_in_opt:
        # Pass opposite to verify EQV
        kw["arr_in_opt"] = container_create(pybmad.BoolAlloc1D, [False, True])
    if use_inout_opt:
        # Pass opposite to verify EQV
        kw["arr_inout_opt"] = arr_inout_opt

    res = pybmad.test_logical_array(**kw)

    if use_in_opt:
        # [T, F] eqv [F, T] -> [F, F]
        assert list(res.arr_out) == [False, False]
    else:
        assert list(res.arr_out) == [True, False]

    assert list(arr_inout) == [False, True]
    if use_inout_opt:
        # [T, F] eqv [F, T] -> [F, F]
        assert list(arr_inout_opt) == [True, False]


@optionals
def test_character_scalar(use_in_opt: bool, use_inout_opt: bool):
    val_in = "foo"
    val_inout = "bar"

    kw = {
        "val_in": val_in,
        "val_inout": val_inout,
        "val_in_opt": "123" if use_in_opt else None,
        "val_inout_opt": "baz" if use_inout_opt else None,
    }

    res = pybmad.test_character_scalar(**kw)

    expected_out = "foo"
    if use_in_opt:
        expected_out += "123"

    assert res.val_out.strip() == expected_out
    # TODO: inout strings not really a thing
    # assert val_inout.strip() == "bar_mod"


@optionals
def test_bunch_struct_scalar(use_in_opt: bool, use_inout_opt: bool):
    """
    Test passing structs (passed proxy class).
    Fortran logic:
       Output = Input
       InOut%ix_ele += 1
       Output%charge_tot += OptInput%charge_tot (if present)
    """

    s_in = pybmad.BunchStruct()
    s_in.ix_ele = 10
    s_in.charge_tot = 100.0

    s_inout = pybmad.BunchStruct()
    s_inout.ix_ele = 5
    s_inout.charge_tot = 50.0

    kw = {"val_in": s_in, "val_inout": s_inout}

    if use_in_opt:
        s_opt = pybmad.BunchStruct()
        s_opt.ix_ele = 999
        s_opt.charge_tot = 1000.0
        kw["val_in_opt"] = s_opt

    if use_inout_opt:
        s_io_opt = pybmad.BunchStruct()
        s_io_opt.ix_ele = 100
        kw["val_inout_opt"] = s_io_opt

    res = pybmad.test_bunch_struct_scalar(**kw)

    # Check InOut Modification in place
    # ix_ele + 1 -> 6
    assert s_inout.ix_ele == 6
    assert s_inout.charge_tot == pytest.approx(51.0)  # + 1.0

    # Check Output copy/modify
    # Copy of in (ix=10, chg=100)
    expected_ix = 10
    expected_charge = 100.0

    if use_in_opt:
        expected_ix += 999
        expected_charge += 1000.0

    assert res.val_out.ix_ele == expected_ix
    assert res.val_out.charge_tot == pytest.approx(expected_charge)

    assert res.opt_status[0] == (1 if use_in_opt else 0)

    if use_inout_opt:
        assert kw["val_inout_opt"].ix_ele == 101
