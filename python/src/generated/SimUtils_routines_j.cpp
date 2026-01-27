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
      py::arg("j_bes"),
      R"""(Wrapper for Fortran routine j_bessel

Parameters
----------
m : int

arg : float

j_bes : float

Returns
-------
m : int

arg : float

j_bes : float
)"""
  );
}
