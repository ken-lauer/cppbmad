#include "pybmad/generated/Bmad_routines_y.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_y(py::module &m) {
  m.def(
      "ylafun",
      &Bmad::ylafun,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine ylafun

Parameters
----------
x : float

y : float

z : float

Returns
-------
res : float
)"""
  );
}
