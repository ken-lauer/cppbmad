#include "pybmad/generated/SimUtils_routines_z.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_z(nb::module_ &m) {
  m.def(
      "zig_table_init",
      &SimUtils::zig_table_init,
      R"""(Private routine to initialize the Ziggurat lookup tables on first use.
Based on Marsaglia & Tsang (2000), "The Ziggurat Method for Generating
Random Variables", Journal of Statistical Software 5(8).
)"""
  );
}
