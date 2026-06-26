#include "pybmad/generated/SimUtils_routines_h.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_h(nb::module_ &m) {
  m.def(
      "hanhan",
      &SimUtils::hanhan,
      nb::arg("N"),
      nb::arg("hh"),
      R"""(Wrapper for Fortran routine hanhan

Parameters
----------
N : int

hh : 1D array of float
)"""
  );
}
