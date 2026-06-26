#include "pybmad/generated/Bmad_routines_u.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_u(nb::module_ &m) {
  m.def(
      "update_ele_from_fibre",
      &Bmad::update_ele_from_fibre,
      nb::arg("ele"),
      R"""(Routine to update a bmad lattice element when the associated PTC fibre has been modified.
Remember to call lattice_bookkeeper after calling this routine.

Parameters
----------
ele : EleStruct
    Element with corresponding ele.ptc_fibre fibre.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Modified element.
)"""
  );
  m.def(
      "update_fibre_from_ele",
      &Bmad::update_fibre_from_ele,
      nb::arg("ele"),
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
      [](FloorPositionStruct &floor, FloorPositionStruct *floor0) {
        auto fn = static_cast<void (*)(FloorPositionStruct &, optional_ref<FloorPositionStruct>)>(
            &Bmad::update_floor_angles
        );
        return fn(floor, ptr_to_opt_ref(floor0));
      },
      nb::arg("floor"),
      nb::arg("floor0") = nb::none(),
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
)"""
  );
}
