#include "pybmad/generated/SimUtils_routines_b.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyBinXCenter python_bin_x_center(int ix_bin, double bin1_x_min, double bin_delta) {
  auto _result = SimUtils::bin_x_center(ix_bin, bin1_x_min, bin_delta);
  auto py_result{PyBinXCenter{_result, ix_bin}};
  return py_result;
}
PyBitSet python_bit_set(int word, int pos, bool set_to_1) {
  SimUtils::bit_set(word, pos, set_to_1);
  auto py_result{PyBitSet{word}};
  return py_result;
}

void init_SimUtils_routines_b(nb::module_ &m) {
  nb::class_<SimUtils::BicubicCmplxEval>(m, "BicubicCmplxEval", "bicubic_cmplx_eval return type")
      .def_ro("df_dx", &SimUtils::BicubicCmplxEval::df_dx)
      .def_ro("df_dy", &SimUtils::BicubicCmplxEval::df_dy)
      .def_ro("f_val", &SimUtils::BicubicCmplxEval::f_val)
      .def("__len__", [](const SimUtils::BicubicCmplxEval &) { return 3; })
      .def("__getitem__", [](const SimUtils::BicubicCmplxEval &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.df_dx);
        if (i == 1)
          return nb::cast(s.df_dy);
        if (i == 2)
          return nb::cast(s.f_val);
        throw nb::index_error();
      });
  m.def(
      "bicubic_cmplx_eval",
      &SimUtils::bicubic_cmplx_eval,
      nb::arg("x_norm"),
      nb::arg("y_norm"),
      nb::arg("bi_coef"),
      R"""(Routine to evaluate a bicubic interpolating complex function.

Use the routine bicubic_interpolation_cmplx_coefs to generate bi_coef.

Note: In the equations below, the four points of the grid box being interpolated range
from (x0, y0) to (x0+dx, y0+dy).

Parameters
----------
x_norm : float
    x_norm = (x - x0) / dx

y_norm : float
    y_norm = (y - y0) / dy

bi_coef : BicubicCmplxCoefStruct
    Coefficients.

Returns
-------
f_val : complex
    Value of f.

df_dx : complex, optional
    Normalized first derivative: True df/dx = df_dx * dx

df_dy : complex, optional
    Normalized first derivative: True df/dy = df_dy * dy
)"""
  );
  m.def(
      "bin_index",
      &SimUtils::bin_index,
      nb::arg("x"),
      nb::arg("bin1_x_min"),
      nb::arg("bin_delta"),
      R"""(Helper function to locate the appropriate histogram bin index.

Parameters
----------
x : float
    Input value to bin.

bin1_x_min : float
    Minimum value of bin with index 1.

bin_delta : float
    Bin width.

Returns
-------
ix_bin : int
    Index of bin x is in.
)"""
  );
  nb::class_<PyBinXCenter>(m, "BinXCenter", "bin_x_center return type")
      .def_ro("x_center", &PyBinXCenter::x_center)
      .def_ro("ix_bin", &PyBinXCenter::ix_bin)
      .def("__len__", [](const PyBinXCenter &) { return 2; })
      .def("__getitem__", [](const PyBinXCenter &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.x_center);
        if (i == 1)
          return nb::cast(s.ix_bin);
        throw nb::index_error();
      });
  m.def(
      "bin_x_center",
      &python_bin_x_center,
      nb::arg("ix_bin"),
      nb::arg("bin1_x_min"),
      nb::arg("bin_delta"),
      R"""(Helper function to locate the center of a histogram bin.

Parameters
----------
ix_bin : int
    Index of bin under question.

bin1_x_min : float
    Minimum value of bin with index 1.

bin_delta : float
    Bin width.

Returns
-------
ix_bin : int
    Index of bin under question.
)"""
  );
  nb::class_<PyBitSet>(m, "BitSet", "bit_set return type")
      .def_ro("word", &PyBitSet::word)
      .def("__len__", [](const PyBitSet &) { return 1; })
      .def("__getitem__", [](const PyBitSet &s, int i) -> nb::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return nb::cast(s.word);
        throw nb::index_error();
      });
  m.def(
      "bit_set",
      &python_bit_set,
      nb::arg("word"),
      nb::arg("pos"),
      nb::arg("set_to_1"),
      R"""(Routine to set a bit in a word.

Parameters
----------
word : int
    Input word
    This parameter is an input/output and is modified in-place.
    As an output, word: Word with bit set.

pos : int
    position to set.

set_to_1 : bool
    If True then bit is set to 1. If False bit is set to 0.

Returns
-------
word : int
    Input word
    This parameter is an input/output and is modified in-place.
    As an output, word: Word with bit set.
)"""
  );
  nb::class_<SimUtils::BracketIndexForSpline>(
      m,
      "BracketIndexForSpline",
      "bracket_index_for_spline return type"
  )
      .def_ro("ix0", &SimUtils::BracketIndexForSpline::ix0)
      .def_ro("ok", &SimUtils::BracketIndexForSpline::ok)
      .def("__len__", [](const SimUtils::BracketIndexForSpline &) { return 2; })
      .def("__getitem__", [](const SimUtils::BracketIndexForSpline &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ix0);
        if (i == 1)
          return nb::cast(s.ok);
        throw nb::index_error();
      });
  m.def(
      "bracket_index_for_spline",
      &SimUtils::bracket_index_for_spline,
      nb::arg("x_knot"),
      nb::arg("x"),
      nb::arg("strict") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to find which interval to use for evaluating a spline.
If strict = False (default), x is in range if
      x_knot(1) - (x_knot(2) - x_knot(1)) < x < x_knot(n) + (x_knot(n) - x_knot(n-1))
If stric = True, x is in range if
      x_knot(1) <= x <= x_knot(n)
where n = size(x_knot)

Parameters
----------
x_knot : 1D array of float
    Array of x values.

x : float
    Evaluation point.

strict : bool, optional
    Default is False. Determines acceptible range.

print_err : bool, optional
    Default is True. Print error message if out of range?

Returns
-------
ix0 : int
    If ok = True, x is in the interval [x_knot(ix0), x_knot(ix0+1)]

ok : bool
    True if x is in range. False otherwise.
)"""
  );
}
