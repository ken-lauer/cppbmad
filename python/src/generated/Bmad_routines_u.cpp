#include "pybmad/generated/Bmad_routines_u.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_u(py::module &m) {
  m.def(
      "update_ele_from_fibre",
      &Bmad::update_ele_from_fibre,
      py::arg("ele"),
      R"""(Subroutine update_ele_from_fibre (ele)

Routine to update a bmad lattice element when the associated PTC fibre has been modified.
Remember to call lattice_bookkeeper after calling this routine.

Parameters
----------
ele : EleStruct
    Element with corresponding ele.ptc_fibre fibre.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Modified element.

Returns
-------
ele : EleStruct
    Element with corresponding ele.ptc_fibre fibre.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Modified element.
)"""
  );
  m.def(
      "update_fibre_from_ele",
      &Bmad::update_fibre_from_ele,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine update_fibre_from_ele

Parameters
----------
ele : EleStruct
    Element with corresponding PTC fibre.

Returns
-------
survey_needed : bool
    Set True if a call to survey will be needed. Calling survey is avoided in this routine to save time if
    multiple elements are being updated.
)"""
  );
  m.def(
      "update_floor_angles",
      &Bmad::update_floor_angles,
      py::arg("floor"),
      py::arg("floor0") = py::none(),
      R"""(Wrapper for Fortran routine update_floor_angles

Parameters
----------
floor : FloorPositionStruct
    Position with input w matrix.
    This parameter is an input/output and is modified in-place.
    As an output, floor: Position with output angles.

floor0 : FloorPositionStruct, optional
    Reference position. There are two solutions related by: [theta, phi, psi] & [pi+theta, pi-phi, pi+psi] If
    floor0 is present, choose the solution "nearest" the angles in floor0.

Returns
-------
floor : FloorPositionStruct
    Position with input w matrix.
    This parameter is an input/output and is modified in-place.
    As an output, floor: Position with output angles.
)"""
  );
}
