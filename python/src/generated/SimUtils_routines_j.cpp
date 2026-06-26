#include "pybmad/generated/SimUtils_routines_j.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_j(nb::module_ &m) {
  m.def(
      "j_bessel",
      &SimUtils::j_bessel,
      nb::arg("m"),
      nb::arg("arg"),
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
