#include "pybmad/generated/bsim_routines_c.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_c(py::module &m) {
  m.def(
      "check_rf_freq",
      &bsim::check_rf_freq,
      py::arg("lat"),
      py::arg("fb"),
      R"""(Wrapper for Fortran routine check_rf_freq

Parameters
----------
lat : LatStruct

fb : float

Returns
-------
lat : LatStruct

fb : float
)"""
  );
  m.def(
      "count_lines_in_file",
      &bsim::count_lines_in_file,
      py::arg("file_name"),
      R"""(Wrapper for Fortran routine count_lines_in_file

Parameters
----------
file_name : character

lines : int

Returns
-------
file_name : character

lines : int
)"""
  );
}
