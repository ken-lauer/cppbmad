#include "pybmad/generated/bsim_routines_l.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_l(py::module &m) {
  m.def(
      "logical_to_python",
      &bsim::logical_to_python,
      py::arg("logic"),
      py::arg("string"),
      R"""(Wrapper for Fortran routine logical_to_python

Parameters
----------
logic : 
string : 
)"""
  );
}
