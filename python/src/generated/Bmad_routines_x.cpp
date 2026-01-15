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
      R"""(Parameters
----------
x : 
y : 
z : 
res : 
)"""
  );
  m.def(
      "xraylib_nist_compound",
      &Bmad::xraylib_nist_compound,
      py::arg("name"),
      R"""(Function xraylib_nist_compound (name) result (indx)

Routine to return the xraylib index for a given NIST compound.
Taken from file xraylib/include/xraylib-nist_compounds.h

Parameters
----------
name : unknown
    Name of compound

Returns
-------
indx : int
    Compound index. -1 if not found.
)"""
  );
}
