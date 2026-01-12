#include "pybmad/generated/bsim_routines_r.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers

void init_bsim_routines_r(py::module& m) {
  m.def(
      "rf_cav_names",
      &bsim::rf_cav_names,
      py::arg("lat"),
      R"""(Parameters
----------
lat : 
)""");
}
