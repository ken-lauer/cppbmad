#include "pybmad/generated/SimUtils_routines_v.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_v(py::module &m) {
  m.def(
      "virtual_memory_usage",
      &SimUtils::virtual_memory_usage,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine virtual_memory_usage

Parameters
----------

Returns
-------
usage : int
)"""
  );
}
