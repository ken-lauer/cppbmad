#include "pybmad/generated/SimUtils_routines_h.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_h(py::module &m) {
  m.def(
      "hanhan",
      &SimUtils::hanhan,
      py::arg("N"),
      py::arg("hh"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine hanhan

Parameters
----------
N : int

hh : 1D array of float
)"""
  );
}
