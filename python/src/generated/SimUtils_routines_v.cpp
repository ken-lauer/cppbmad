#include "pybmad/generated/SimUtils_routines_v.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers

void init_SimUtils_routines_v(py::module& m) {
  m.def(
      "virtual_memory_usage",
      &SimUtils::virtual_memory_usage,
      R"""(Parameters
----------
usage : 
)""");
}
