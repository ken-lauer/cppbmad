#include "pybmad/generated/bsim_routines_w.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers

void init_bsim_routines_w(py::module& m) {
  m.def(
      "write_bunch_by_bunch_info",
      &bsim::write_bunch_by_bunch_info,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("this_stage"),
      R"""(Parameters
----------
lat : 
bbu_beam : 
bbu_param : 
this_stage : 
)""");
}
