#include "pybmad/generated/Bmad_routines_v.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_v(py::module &m) {
  m.def(
      "valid_field_calc",
      &Bmad::valid_field_calc,
      py::arg("ele"),
      py::arg("field_calc"),
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
      py::arg("ele"),
      py::arg("fringe_type"),
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
      py::arg("ele"),
      py::arg("species"),
      py::arg("mat6_calc_method"),
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
      py::arg("ele"),
      py::arg("spin_tracking_method"),
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
      py::arg("ele"),
      py::arg("species"),
      py::arg("tracking_method"),
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
  py::class_<Bmad::ValueOfAttribute, std::unique_ptr<Bmad::ValueOfAttribute>>(
      m,
      "ValueOfAttribute",
      "value_of_attribute return type"
  )
      .def_readonly("err_flag", &Bmad::ValueOfAttribute::err_flag)
      .def_readonly("value", &Bmad::ValueOfAttribute::value)
      .def("__len__", [](const Bmad::ValueOfAttribute &) { return 2; })
      .def("__getitem__", [](const Bmad::ValueOfAttribute &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "value_of_attribute",
      &Bmad::value_of_attribute,
      py::arg("ele"),
      py::arg("attrib_name"),
      py::arg("err_print_flag") = py::none(),
      py::arg("err_value") = py::none(),
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
      py::arg("line"),
      py::arg("value"),
      py::arg("str"),
      py::arg("typ"),
      py::arg("ignore_if_zero") = py::none(),
      py::arg("use_comma") = py::none(),
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
      py::arg("vec"),
      py::arg("phase") = py::none(),
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
      py::arg("vec"),
      py::arg("phase") = py::none(),
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
