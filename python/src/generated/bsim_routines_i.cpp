#include "pybmad/generated/bsim_routines_i.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers

void init_bsim_routines_i(py::module& m) {
  m.def(
      "insert_phase_trombone",
      &bsim::insert_phase_trombone,
      py::arg("branch"),
      R"""(Parameters
----------
branch : BranchStruct
    Lattice branch.
    This parameter is an input/output and is modified in-place. As an output: Lattice branch with trumbone at
    branch.ele(1).
)""");
}
