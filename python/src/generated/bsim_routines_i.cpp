#include "pybmad/generated/bsim_routines_i.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_i(nb::module_ &m) {
  m.def(
      "insert_phase_trombone",
      &bsim::insert_phase_trombone,
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine insert_phase_trombone

Parameters
----------
branch : BranchStruct
    Lattice branch.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Lattice branch with trumbone at branch.ele(1).
)"""
  );
}
