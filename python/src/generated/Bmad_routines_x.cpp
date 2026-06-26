#include "pybmad/generated/Bmad_routines_x.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_x(nb::module_ &m) {
  m.def(
      "xlafun",
      &Bmad::xlafun,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("z"),
      R"""(Wrapper for Fortran routine xlafun

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
  m.def(
      "xraylib_nist_compound",
      &Bmad::xraylib_nist_compound,
      nb::arg("name"),
      R"""(Routine to return the xraylib index for a given NIST compound.
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
