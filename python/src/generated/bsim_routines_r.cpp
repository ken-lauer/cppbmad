#include "pybmad/generated/bsim_routines_r.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_r(nb::module_ &m) {
  m.def(
      "rf_cav_names",
      &bsim::rf_cav_names,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine rf_cav_names

Parameters
----------
lat : LatStruct
)"""
  );
}
