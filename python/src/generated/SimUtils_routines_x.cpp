#include "pybmad/generated/SimUtils_routines_x.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_x(py::module& m) {
  m.def(
      "x0_radiation_length",
      &SimUtils::x0_radiation_length,
      py::arg("species"),
      R"""(Function x0_radiation_length(species) result (x0)

Routine to return the X0 raidation length for atomes.

Parameters
----------
Species : int
    Species ID.

Returns
-------
x0 : float
    Radiation length in kg/m^2. Set to real_garbage$ if species is not atomic or has atomic index greater than
    92.
)""");
}
