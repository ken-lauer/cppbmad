#include "pybmad/generated/bsim_routines_h.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_h(py::module &m) {
  m.def(
      "hom_voltage",
      &bsim::hom_voltage,
      py::arg("lr_wake"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine hom_voltage

Parameters
----------
lr_wake : WakeLrModeStruct

Returns
-------
voltage : float
)"""
  );
}
