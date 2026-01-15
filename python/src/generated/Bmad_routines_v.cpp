#include "pybmad/generated/Bmad_routines_v.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_v(py::module& m) {
  m.def(
      "valid_field_calc",
      &Bmad::valid_field_calc,
      py::arg("ele"),
      py::arg("field_calc"),
      py::arg("is_valid"),
      R"""(Parameters
----------
ele : EleStruct
    Lattice element.
field_calc : int
    bmad_standard$, etc.
is_valid : 
)""");
  m.def(
      "valid_fringe_type",
      &Bmad::valid_fringe_type,
      py::arg("ele"),
      py::arg("fringe_type"),
      py::arg("is_valid"),
      R"""(Parameters
----------
ele : EleStruct
    Lattice element.
fringe_type : int
    bmad_standard$, etc.
is_valid : 
)""");
  m.def(
      "valid_mat6_calc_method",
      &Bmad::valid_mat6_calc_method,
      py::arg("ele"),
      py::arg("species"),
      py::arg("mat6_calc_method"),
      py::arg("is_valid"),
      R"""(Parameters
----------
ele : EleStruct
    Lattice element.
species : 
    Type of particle being tracked. electron$, etc. or not_set$
mat6_calc_method : int
    bmad_standard$, etc.
is_valid : 
)""");
  m.def(
      "valid_spin_tracking_method",
      &Bmad::valid_spin_tracking_method,
      py::arg("ele"),
      py::arg("spin_tracking_method"),
      py::arg("is_valid"),
      R"""(Parameters
----------
ele : EleStruct
    Lattice element.
spin_tracking_method : int
    bmad_standard$, etc.
is_valid : 
)""");
  m.def(
      "valid_tracking_method",
      &Bmad::valid_tracking_method,
      py::arg("ele"),
      py::arg("species"),
      py::arg("tracking_method"),
      py::arg("is_valid"),
      R"""(Parameters
----------
ele : EleStruct
    Lattice element.
species : 
    Type of particle being tracked. electron$, etc. or not_set$
tracking_method : int
    bmad_standard$, etc.
is_valid : 
)""");
  m.def(
      "value_of_attribute",
      &Bmad::value_of_attribute,
      py::arg("ele"),
      py::arg("attrib_name"),
      py::arg("err_print_flag") = py::none(),
      py::arg("err_value") = py::none(),
      py::arg("value"),
      R"""(Parameters
----------
ele : EleStruct
    After this routine finishes Ptr_attrib will point to a variable within this element.
attrib_name : unknown
    Name of attribute. Must be uppercase. For example: "HKICK".
err_flag : bool
    Set True if attribtute not found. False otherwise.
err_print_flag : bool, optional
    If present and True then print an error message if there is an  error.
err_value : float, optional
    Value to set value argument if there is an error. Default is 0.
value : 
)""");
  m.def(
      "value_to_line",
      &Bmad::value_to_line,
      py::arg("line"),
      py::arg("value"),
      py::arg("str"),
      py::arg("typ"),
      py::arg("ignore_if_zero") = py::none(),
      py::arg("use_comma") = py::none(),
      R"""(Parameters
----------
line : 
value : 
str : 
typ : 
ignore_if_zero : 
use_comma : 
)""");
  m.def(
      "vec_to_polar",
      &Bmad::vec_to_polar,
      py::arg("vec"),
      py::arg("phase") = py::none(),
      py::arg("polar"),
      R"""(Parameters
----------
vec : float
    unitary spin vector
phase : float, optional
    Phase of the spinor, if not given then set to zero
polar : 
)""");
  m.def(
      "vec_to_spinor",
      &Bmad::vec_to_spinor,
      py::arg("vec"),
      py::arg("phase") = py::none(),
      py::arg("spinor"),
      R"""(Parameters
----------
vec : float
    Spin vector in cartesian coordinates
phase : float
    Phase of the spinor, if not given then set to zero
spinor : 
)""");
  m.def(
      "verify_valid_name",
      &Bmad::verify_valid_name,
      py::arg("name"),
      py::arg("ix_name"),
      py::arg("pure_name") = py::none(),
      py::arg("include_wild") = py::none(),
      R"""(Function verify_valid_name (name, ix_name, pure_name, include_wild) result (is_valid)

Routine to check if a name is well formed. Examples:
  "0>>Q0"                           -- Invalid (will only be valid after lattice expansion).
  "Q1##1"                           -- Invalid (double hash not accepted).
  "Q2A_C.\7#"                       -- Pure name (no "[", "]", "(", ")", "%" characters present).
  "Q3[GRID_FIELD(1)%FIELD_SCALE]"   -- Valid but not a pure name.
  "RFCAVITY::*"                     -- Valid if include_wild = True.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.


Parameters
----------
name : unknown
    Name(1:ix_name) is the string to check.
ix_name : int
    Number of characters in the name.
pure_name : bool, optional
    If True, reject names that contain "[", "]", "(", ")", "." characters. Default is False.
include_wild : bool, optional
    Name can include wild card characters and additionally type prefixes like "QUAD::". Default is False.

Returns
-------
is_valid : bool
    True if name is well formed. False otherwise.
)""");
}
