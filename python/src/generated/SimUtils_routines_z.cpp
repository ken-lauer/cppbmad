#include "pybmad/generated/SimUtils_routines_z.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_z(py::module &m) {
  m.def(
      "zig_table_init",
      &SimUtils::zig_table_init,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine zig_table_init ()

Private routine to initialize the Ziggurat lookup tables on first use.
Based on Marsaglia & Tsang (2000), "The Ziggurat Method for Generating
Random Variables", Journal of Statistical Software 5(8).
)"""
  );
}
