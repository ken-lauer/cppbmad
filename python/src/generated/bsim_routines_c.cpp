#include "pybmad/generated/bsim_routines_c.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_c(py::module& m) {
  m.def(
      "check_rf_freq",
      &bsim::check_rf_freq,
      py::arg("lat"),
      py::arg("fb"),
      R"""(Parameters
----------
lat : 
fb : 
)""");
  m.def(
      "count_lines_in_file",
      &bsim::count_lines_in_file,
      py::arg("file_name"),
      R"""(Parameters
----------
file_name : 
lines : 
)""");
}
