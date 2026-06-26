#include "pybmad/generated/SimUtils_routines_x.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_x(nb::module_ &m) {
  m.def(
      "x0_radiation_length",
      &SimUtils::x0_radiation_length,
      nb::arg("species"),
      R"""(Routine to return the X0 raidation length for atomes.

Parameters
----------
Species : int
    Species ID.

Returns
-------
x0 : float
    Radiation length in kg/m^2. Set to real_garbage$ if species is not atomic or has atomic index greater than
    92.
)"""
  );
}
