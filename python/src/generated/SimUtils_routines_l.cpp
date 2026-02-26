#include "pybmad/generated/SimUtils_routines_l.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_l(py::module &m) {
  m.def(
      "linear_fit",
      &SimUtils::linear_fit,
      py::arg("x"),
      py::arg("y"),
      py::arg("n_data"),
      py::arg("a"),
      py::arg("b"),
      py::arg("sig_a"),
      py::arg("sig_b"),
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("logic"),
      py::arg("str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine logic_str

Parameters
----------
logic : bool

str : str
)"""
  );
  m.def(
      "lunget",
      &SimUtils::lunget,
      py::arg("func_retval__"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine lunget

Parameters
----------
)"""
  );
}
