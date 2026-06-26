#include "pybmad/generated/Bmad_routines_y.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_y(nb::module_ &m) {
  m.def(
      "ylafun",
      &Bmad::ylafun,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("z"),
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
