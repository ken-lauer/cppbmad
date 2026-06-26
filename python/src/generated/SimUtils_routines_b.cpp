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
  nb::class_<SimUtils::BicubicEval>(m, "BicubicEval", "bicubic_eval return type")
      .def_ro("df_dx", &SimUtils::BicubicEval::df_dx)
      .def_ro("df_dy", &SimUtils::BicubicEval::df_dy)
      .def_ro("f_val", &SimUtils::BicubicEval::f_val)
      .def("__len__", [](const SimUtils::BicubicEval &) { return 3; })
      .def("__getitem__", [](const SimUtils::BicubicEval &s, int i) -> nb::object {
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
      "bicubic_eval",
      &SimUtils::bicubic_eval,
      nb::arg("x_norm"),
      nb::arg("y_norm"),
      nb::arg("bi_coef"),
      R"""(Routine to evaluate a bicubic interpolating function.

Use the routine bicubic_interpolation_coefs to generate bi_coef.

Note: In the equations below, the four points of the grid box being interpolated range
from (x0, y0) to (x0+dx, y0+dy).

Parameters
----------
x_norm : float
    x_norm = (x - x0) / dx

y_norm : float
    y_norm = (y - y0) / dy

bi_coef : BicubicCoefStruct
    Coefficients.

Returns
-------
f_val : float
    Value of f.

df_dx : float, optional
    Normalized first derivative: True df/dx = df_dx * dx

df_dy : float, optional
    Normalized first derivative: True df/dy = df_dy * dy
)"""
  );
  m.def(
      "bicubic_interpolation_cmplx_coefs",
      &SimUtils::bicubic_interpolation_cmplx_coefs,
      nb::arg("field_at_box"),
      R"""(Routine to compute the complex coefficients for bicubic interpolation.

Use the routine bicubic_cmplx_eval to evaluate the interpolation function.

Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy

Parameters
----------
field_at_box : CmplxFieldAt2dBoxStruct
    Field and normalized derivatives at the 4 grid points.

Returns
-------
bi_coef : BicubicCmplxCoefStruct
    Coefficients.
)"""
  );
  m.def(
      "bicubic_interpolation_coefs",
      &SimUtils::bicubic_interpolation_coefs,
      nb::arg("field_at_box"),
      R"""(Routine to compute the coefficients for bicubic interpolation.

Use the routine bicubic_eval to evaluate the interpolation function.

Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy

Parameters
----------
field_at_box : FieldAt2dBoxStruct
    Field and normalized derivatives at the 4 grid points.

Returns
-------
bi_coef : BicubicCoefStruct
    Coefficients.
)"""
  );
  m.def(
      "bin_2d",
      &SimUtils::bin_2d,
      nb::arg("data1"),
      nb::arg("data2"),
      nb::arg("weight") = nb::none(),
      nb::arg("min1") = nb::none(),
      nb::arg("max1") = nb::none(),
      nb::arg("min2") = nb::none(),
      nb::arg("max2") = nb::none(),
      nb::arg("n_bins1") = nb::none(),
      nb::arg("n_bins2") = nb::none(),
      R"""(         result (bin_data)

Similiar to bin(...), but for two dimensions.

Parameters
----------
data1 : 1D array of float
    data to bin in dimension 1

data2 : 1D array of float
    data to bin in dimension 2

weight : 1D array of float, optional
    1D weights for each data. Default: 1

min1 : float, optional
    minimum considered for data1. Default: minval(data1)

max1 : float, optional
    maximum considered for data1. Default: maxval(data1)

min2 : float, optional
    minimum considered for data2. Default: minval(data2)

max2 : float, optional
    maximum considered for data2. Default: maxva2(data1)

n_bins1 : int, optional
    number of bins for dimension 1. Default: 2*size(data1)^(1/6)

n_bins2 : int, optional
    number of bins for dimension 2. Default: 2*size(data2)^(1/6)

Returns
-------
bin_data : GeneralBinStruct
    binned data, with .dim==2
)"""
  );
  m.def(
      "bin_data",
      &SimUtils::bin_data,
      nb::arg("data"),
      nb::arg("weight") = nb::none(),
      nb::arg("min") = nb::none(),
      nb::arg("max") = nb::none(),
      nb::arg("n_bins") = nb::none(),
      R"""(Bin centers are at [ x_min + 1/2*delta, x_min + 3/2*delta, ..., x_max - 1/2*delta ]
with delta  = (max-min)/n_bins

Parameters
----------
data : 1D array of float
    1D data to bin

weight : 1D array of float, optional
    1D weights for each data. Default: 1

min : float, optional
    minimum considered. Default: minval(data)

max : float, optional
    maximum considered. Default: maxval(data)

n_bins : int, optional
    number of bins. Default: 2*size(data)^(1/3) (Rice's rule)

Returns
-------
binned_data : BinStruct
    binned data.
)"""
  );
  m.def(
      "bin_data_density",
      &SimUtils::bin_data_density,
      nb::arg("bin_data"),
      nb::arg("x"),
      nb::arg("order") = nb::none(),
      R"""(Calculate the density of binned data at an arbitrary location x.
Zero is returned if x is out of bounds of the binned data.

Parameters
----------
bin_data : BinStruct
    binned data

x : float
    position to query

order : int, optional
    interpolation order: 0 or 1 only. Default: 1

Returns
-------
r : float
    density at x
)"""
  );
  m.def(
      "bin_data_density_2d",
      &SimUtils::bin_data_density_2d,
      nb::arg("bin_data"),
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("order") = nb::none(),
      R"""(Calculate the density of 2D binned data at an arbitrary location (x, y)
Zero is returned if x or y is out of bounds of the binned data

Parameters
----------
bin_data : GeneralBinStruct
    binned data. Must have .dim == 2

x : float
    x position to query

y : float
    y position to query

order : int, optional
    interpolation order: 0 or 1 only Default: 1

Returns
-------
r0 : float
    density at (x,y)
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
