#include "pybmad/generated/Tao_routines_i.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Tao_routines_i(py::module &m) {
  m.def(
      "integrate_max",
      &Tao::integrate_max,
      py::arg("ix_start"),
      py::arg("ix_ele"),
      py::arg("datum_value"),
      py::arg("ix_m"),
      py::arg("branch"),
      py::arg("vec"),
      py::arg("datum"),
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
      py::arg("ix_start"),
      py::arg("ix_ele"),
      py::arg("datum_value"),
      py::arg("ix_m"),
      py::arg("branch"),
      py::arg("vec"),
      py::arg("datum"),
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
