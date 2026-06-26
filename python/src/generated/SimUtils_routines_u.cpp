#include "pybmad/generated/SimUtils_routines_u.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_u(nb::module_ &m) {
  m.def(
      "upcase_string",
      &SimUtils::upcase_string,
      nb::arg("string"),
      R"""(Wrapper for Fortran routine upcase_string

Parameters
----------
string : str
)"""
  );
}
