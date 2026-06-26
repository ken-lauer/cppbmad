#include "pybmad/generated/bsim_routines_c.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_c(nb::module_ &m) {
  m.def(
      "check_rf_freq",
      &bsim::check_rf_freq,
      nb::arg("lat"),
      nb::arg("fb"),
      R"""(Wrapper for Fortran routine check_rf_freq

Parameters
----------
lat : LatStruct

fb : float
)"""
  );
  m.def(
      "count_lines_in_file",
      &bsim::count_lines_in_file,
      nb::arg("file_name"),
      R"""(Wrapper for Fortran routine count_lines_in_file

Parameters
----------
file_name : str

Returns
-------
lines : int
)"""
  );
}
