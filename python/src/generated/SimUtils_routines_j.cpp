#include "pybmad/generated/SimUtils_routines_j.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_j(py::module &m) {
  m.def(
      "j_bessel",
      &SimUtils::j_bessel,
      py::arg("m"),
      py::arg("arg"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine j_bessel

Parameters
----------
m : int

arg : float
    Bessel argument.

Returns
-------
j_bes : float
    Bessel value.
)"""
  );
}
