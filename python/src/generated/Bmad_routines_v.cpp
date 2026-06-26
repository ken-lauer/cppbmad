#include "pybmad/generated/Bmad_routines_v.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_v(nb::module_ &m) {
  m.def(
      "valid_field_calc",
      &Bmad::valid_field_calc,
      nb::arg("ele"),
      nb::arg("field_calc"),
      R"""(Wrapper for Fortran routine valid_field_calc

Parameters
----------
ele : EleStruct
    Lattice element.

field_calc : int
    bmad_standard$, etc.

Returns
-------
is_valid : bool
    True if a valid method. False otherwise.
)"""
  );
  m.def(
      "valid_fringe_type",
      &Bmad::valid_fringe_type,
      nb::arg("ele"),
      nb::arg("fringe_type"),
      R"""(Wrapper for Fortran routine valid_fringe_type

Parameters
----------
ele : EleStruct
    Lattice element.

fringe_type : int
    bmad_standard$, etc.

Returns
-------
is_valid : bool
    True if a valid method. False otherwise.
)"""
  );
  m.def(
      "valid_mat6_calc_method",
      &Bmad::valid_mat6_calc_method,
      nb::arg("ele"),
      nb::arg("species"),
      nb::arg("mat6_calc_method"),
      R"""(Wrapper for Fortran routine valid_mat6_calc_method

Parameters
----------
ele : EleStruct
    Lattice element.

species : int
    Type of particle being tracked. electron$, etc. or not_set$

mat6_calc_method : int
    bmad_standard$, etc.

Returns
-------
is_valid : bool
    True if a valid method. False otherwise.
)"""
  );
  m.def(
      "valid_spin_tracking_method",
      &Bmad::valid_spin_tracking_method,
      nb::arg("ele"),
      nb::arg("spin_tracking_method"),
      R"""(Wrapper for Fortran routine valid_spin_tracking_method

Parameters
----------
ele : EleStruct
    Lattice element.

spin_tracking_method : int
    bmad_standard$, etc.

Returns
-------
is_valid : bool
    True if a valid method. False otherwise.
)"""
  );
  m.def(
      "valid_tracking_method",
      &Bmad::valid_tracking_method,
      nb::arg("ele"),
      nb::arg("species"),
      nb::arg("tracking_method"),
      R"""(Wrapper for Fortran routine valid_tracking_method

Parameters
----------
ele : EleStruct
    Lattice element.

species : int
    Type of particle being tracked. electron$, etc. or not_set$

tracking_method : int
    bmad_standard$, etc.

Returns
-------
is_valid : bool
    True if a valid method. False otherwise.
)"""
  );
  nb::class_<Bmad::ValueOfAttribute>(m, "ValueOfAttribute", "value_of_attribute return type")
      .def_ro("err_flag", &Bmad::ValueOfAttribute::err_flag)
      .def_ro("value", &Bmad::ValueOfAttribute::value)
      .def("__len__", [](const Bmad::ValueOfAttribute &) { return 2; })
      .def("__getitem__", [](const Bmad::ValueOfAttribute &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "value_of_attribute",
      &Bmad::value_of_attribute,
      nb::arg("ele"),
      nb::arg("attrib_name"),
      nb::arg("err_print_flag") = nb::none(),
      nb::arg("err_value") = nb::none(),
      R"""(Wrapper for Fortran routine value_of_attribute

Parameters
----------
ele : EleStruct
    After this routine finishes Ptr_attrib will point to a variable within this element.

attrib_name : str
    Name of attribute. Must be uppercase. For example: "HKICK".

err_print_flag : bool, optional
    If present and True then print an error message if there is an  error.

err_value : float, optional
    Value to set value argument if there is an error. Default is 0.

Returns
-------
value : float
    Value of the attribute. Set to err_value if not found.

err_flag : bool, optional
    Set True if attribtute not found. False otherwise.
)"""
  );
  m.def(
      "value_to_line",
      &Bmad::value_to_line,
      nb::arg("line"),
      nb::arg("value"),
      nb::arg("str"),
      nb::arg("typ"),
      nb::arg("ignore_if_zero") = nb::none(),
      nb::arg("use_comma") = nb::none(),
      R"""(Wrapper for Fortran routine value_to_line

Parameters
----------
line : str

value : float

str : str

typ : str

ignore_if_zero : bool, optional

use_comma : bool, optional
)"""
  );
  m.def(
      "vec_to_polar",
      &Bmad::vec_to_polar,
      nb::arg("vec"),
      nb::arg("phase") = nb::none(),
      R"""(Wrapper for Fortran routine vec_to_polar

Parameters
----------
vec : 1D array of float (shape: 3)
    unitary spin vector

phase : float, optional
    Phase of the spinor, if not given then set to zero

Returns
-------
polar : SpinPolarStruct
)"""
  );
  m.def(
      "vec_to_spinor",
      &Bmad::vec_to_spinor,
      nb::arg("vec"),
      nb::arg("phase") = nb::none(),
      R"""(Wrapper for Fortran routine vec_to_spinor

Parameters
----------
vec : 1D array of float (shape: 3)
    Spin vector in cartesian coordinates

phase : float, optional
    Phase of the spinor, if not given then set to zero

Returns
-------
spinor : 1D array of complex (shape: 2)
    Spinor.
)"""
  );
  m.def(
      "verify_valid_name",
      &Bmad::verify_valid_name,
      nb::arg("name"),
      nb::arg("ix_name"),
      nb::arg("pure_name") = nb::none(),
      nb::arg("include_wild") = nb::none(),
      R"""(Routine to check if a name is well formed. Examples:
  "0>>Q0"                           -- Invalid (will only be valid after lattice expansion).
  "Q1##1"                           -- Invalid (double hash not accepted).
  "Q2A_C.\7#"                       -- Pure name (no "[", "]", "(", ")", "%" characters present).
  "Q3[GRID_FIELD(1)%FIELD_SCALE]"   -- Valid but not a pure name.
  "RFCAVITY::*"                     -- Valid if include_wild = True.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
name : str
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
)"""
  );
}
