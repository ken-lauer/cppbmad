#include "pybmad/generated/bsim_routines_h.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_h(py::module &m) {
  m.def(
      "hom_voltage",
      &bsim::hom_voltage,
      py::arg("lr_wake"),
      py::arg("voltage"),
      R"""(Wrapper for Fortran routine hom_voltage

Parameters
----------
lr_wake : 
voltage : 
)"""
  );
}
