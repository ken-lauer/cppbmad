#include "pybmad/generated/bsim_routines_l.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_l(nb::module_ &m) {
  m.def(
      "logical_to_python",
      &bsim::logical_to_python,
      nb::arg("logic"),
      R"""(Wrapper for Fortran routine logical_to_python

Parameters
----------
logic : bool

Returns
-------
string : str
)"""
  );
}
