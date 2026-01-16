#include "pybmad/generated/bsim_routines_s.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_s(py::module &m) {
  m.def(
      "set_tune_3d",
      &bsim::set_tune_3d,
      py::arg("branch"),
      py::arg("target_tunes"),
      py::arg("mask") = py::none(),
      py::arg("use_phase_trombone") = py::none(),
      py::arg("z_tune_set") = py::none(),
      py::arg("group_knobs") = py::none(),
      py::arg("print_err") = py::none(),
      R"""(Parameters
----------
branch : BranchStruct
    This parameter is an input/output and is modified in-place. As an output: with adjusted quads and RF to
    match desired tunes.
target_tunes : float
    tunes for a, b, z modes (rad/2pi). Must include integer part.
mask : 
use_phase_trombone : bool, optional
    Default False. If true, use a match element in phase trombone mode to adjust the tunes.
z_tune_set : bool, optional
    Default True. If false, do not try to set the synch tune.
group_knobs : unknown, optional
    If set non-blank, use these group elements for tuning.
print_err : bool, optional
    Print error message if there is a problem? Default is True.
everything_ok : bool
    Returns true or false if set was successful.
)"""
  );
}
