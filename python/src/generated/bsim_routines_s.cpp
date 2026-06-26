#include "pybmad/generated/bsim_routines_s.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_s(nb::module_ &m) {
  m.def(
      "set_tune_3d",
      &bsim::set_tune_3d,
      nb::arg("branch"),
      nb::arg("target_tunes"),
      nb::arg("mask") = nb::none(),
      nb::arg("use_phase_trombone") = nb::none(),
      nb::arg("z_tune_set") = nb::none(),
      nb::arg("group_knobs") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Wrapper for Fortran routine set_tune_3d

Parameters
----------
branch : BranchStruct
    This parameter is an input/output and is modified in-place.
    As an output, branch: with adjusted quads and RF to match desired tunes.

target_tunes : 1D array of float (shape: 3)
    tunes for a, b, z modes (rad/2pi). Must include integer part.

mask : str, optional

use_phase_trombone : bool, optional
    Default False. If true, use a match element in phase trombone mode to adjust the tunes. The match element
    must be the first element in the lattice. Use insert_phase_trombone to insert one.

z_tune_set : bool, optional
    Default True. If false, do not try to set the synch tune.

group_knobs : 1D array of str (shape: 2), optional
    If set non-blank, use these group elements for tuning.

print_err : bool, optional
    Print error message if there is a problem? Default is True.

Returns
-------
everything_ok : bool
    Returns true or false if set was successful.
)"""
  );
}
