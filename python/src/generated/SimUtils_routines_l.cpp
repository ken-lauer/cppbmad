#include "pybmad/generated/SimUtils_routines_l.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_l(nb::module_ &m) {
  m.def(
      "linear_fit",
      &SimUtils::linear_fit,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("n_data"),
      nb::arg("a"),
      nb::arg("b"),
      nb::arg("sig_a"),
      nb::arg("sig_b"),
      R"""(Wrapper for Fortran routine linear_fit

Parameters
----------
x : 1D array of float

y : 1D array of float

n_data : int

a : float

b : float

sig_a : float

sig_b : float
)"""
  );
  m.def(
      "linear_fit_2d",
      &SimUtils::linear_fit_2d,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("z"),
      R"""(Wrapper for Fortran routine linear_fit_2d

Parameters
----------
x : 1D array of float
    Array of x-values.

y : 1D array of float
    Array of y-values.

z : 1D array of float
    Array of z-values

Returns
-------
coef : 1D array of float (shape: 3)
    Coefficients of the linear fit
)"""
  );
  m.def(
      "logic_str",
      &SimUtils::logic_str,
      nb::arg("logic"),
      R"""(Wrapper for Fortran routine logic_str

Parameters
----------
logic : bool

Returns
-------
str : str
)"""
  );
  m.def(
      "lunget",
      &SimUtils::lunget,
      R"""(Wrapper for Fortran routine lunget

Parameters
----------
)"""
  );
}
