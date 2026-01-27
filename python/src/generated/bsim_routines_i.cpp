#include "pybmad/generated/bsim_routines_i.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_i(py::module &m) {
  m.def(
      "insert_phase_trombone",
      &bsim::insert_phase_trombone,
      py::arg("branch"),
      R"""(Wrapper for Fortran routine insert_phase_trombone

Parameters
----------
branch : BranchStruct
    Lattice branch.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Lattice branch with trumbone at branch.ele(1).

Returns
-------
branch : BranchStruct
    Lattice branch.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Lattice branch with trumbone at branch.ele(1).
)"""
  );
}
