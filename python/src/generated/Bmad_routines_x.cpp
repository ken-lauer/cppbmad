#include "pybmad/generated/Bmad_routines_x.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_x(py::module &m) {
  m.def(
      "xlafun",
      &Bmad::xlafun,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::arg("res"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine xlafun

Parameters
----------
x : float

y : float

z : float

res : float
)"""
  );
  m.def(
      "xraylib_nist_compound",
      &Bmad::xraylib_nist_compound,
      py::arg("name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function xraylib_nist_compound (name) result (indx)

Routine to return the xraylib index for a given NIST compound.
Taken from file xraylib/include/xraylib-nist_compounds.h

Parameters
----------
name : str
    Name of compound

Returns
-------
indx : int
    Compound index. -1 if not found.
)"""
  );
}
