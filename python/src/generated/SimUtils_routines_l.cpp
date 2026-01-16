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
      R"""(Parameters
----------
x : 
y : 
n_data : 
a : 
b : 
sig_a : 
sig_b : 
)"""
  );
  m.def(
      "linear_fit_2d",
      &SimUtils::linear_fit_2d,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      R"""(Parameters
----------
x : float
    Array of x-values.
y : float
    Array of y-values.
z : float
    Array of z-values
coef : float
    Coefficients of the linear fit
)"""
  );
  m.def(
      "logic_str",
      &SimUtils::logic_str,
      py::arg("logic"),
      py::arg("str"),
      R"""(Parameters
----------
logic : 
str : 
)"""
  );
  m.def(
      "lunget",
      &SimUtils::lunget,
      py::arg("func_retval__"),
      R"""(Parameters
----------
lunget : 
)"""
  );
}
