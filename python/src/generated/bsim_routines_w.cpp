#include "pybmad/generated/bsim_routines_w.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_w(nb::module_ &m) {
  m.def(
      "write_bunch_by_bunch_info",
      &bsim::write_bunch_by_bunch_info,
      nb::arg("lat"),
      nb::arg("bbu_beam"),
      nb::arg("bbu_param"),
      nb::arg("this_stage"),
      R"""(Wrapper for Fortran routine write_bunch_by_bunch_info

Parameters
----------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

this_stage : BbuStageStruct
)"""
  );
}
