#include "pybmad/generated/SimUtils_routines_u.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_u(py::module &m) {
  m.def(
      "upcase_string",
      &SimUtils::upcase_string,
      py::arg("string"),
      R"""(Wrapper for Fortran routine upcase_string

Parameters
----------
string : str
)"""
  );
}
