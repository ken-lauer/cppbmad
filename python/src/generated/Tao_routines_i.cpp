#include "pybmad/generated/Tao_routines_i.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Tao_routines_i(nb::module_ &m) {
  m.def(
      "integrate_max",
      &Tao::integrate_max,
      nb::arg("ix_start"),
      nb::arg("ix_ele"),
      nb::arg("datum_value"),
      nb::arg("ix_m"),
      nb::arg("branch"),
      nb::arg("vec"),
      nb::arg("datum"),
      R"""(Wrapper for Fortran routine integrate_max

Parameters
----------
ix_start : int

ix_ele : int

datum_value : float

ix_m : int

branch : BranchStruct

vec : 1D array of float

datum : TaoDataStruct
)"""
  );
  m.def(
      "integrate_min",
      &Tao::integrate_min,
      nb::arg("ix_start"),
      nb::arg("ix_ele"),
      nb::arg("datum_value"),
      nb::arg("ix_m"),
      nb::arg("branch"),
      nb::arg("vec"),
      nb::arg("datum"),
      R"""(Wrapper for Fortran routine integrate_min

Parameters
----------
ix_start : int

ix_ele : int

datum_value : float

ix_m : int

branch : BranchStruct

vec : 1D array of float

datum : TaoDataStruct
)"""
  );
}
