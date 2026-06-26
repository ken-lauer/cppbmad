#include "pybmad/generated/bsim_routines_h.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_h(nb::module_ &m) {
  m.def(
      "hom_voltage",
      &bsim::hom_voltage,
      nb::arg("lr_wake"),
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
