#include "pybmad/generated/SimUtils_routines_v.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_v(py::module &m) {
  m.def(
      "virtual_memory_usage",
      &SimUtils::virtual_memory_usage,
      R"""(Parameters
----------
usage : 
)"""
  );
}
