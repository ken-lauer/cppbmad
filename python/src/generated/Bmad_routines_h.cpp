#include "pybmad/generated/Bmad_routines_h.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_h(py::module &m) {
  m.def(
      "hard_multipole_edge_kick",
      &Bmad::hard_multipole_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine hard_multipole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

Routine to track through the hard edge field of a multipole.
The dipole component is ignored and only quadrupole and higher multipoles are included.

This routine handles elements of type:
  sad_mult, sbend, quadrupole, sextupole

For sad_mult elements, ele%a_pole and ele%b_pole ae used for the multipole values.
For the other elements, k1 or k2 is used and it is assumed that we are in the element
frame of reference so tilt = 0.

Parameters
----------
ele : EleStruct
    Element with fringe.

param : LatParamStruct
    Tracking parameters.

particle_at : int
    Either first_track_edge$ or second_track_edge$.

orbit : CoordStruct
    Starting coordinates.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coordinates.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the fringe.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix including the fringe.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.
)"""
  );
  m.def(
      "has_attribute",
      &Bmad::has_attribute,
      py::arg("ele"),
      py::arg("attrib"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine has_attribute

Parameters
----------
ele : EleStruct

attrib : str

Returns
-------
has_it : bool
)"""
  );
  m.def(
      "has_curvature",
      &Bmad::has_curvature,
      py::arg("phot_ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function has_curvature (phot_ele) result (curved)

Routine to determine if a surface is potentially curved or is flat.

Parameters
----------
phot_ele : PhotonElementStruct
    From ele.photon

Returns
-------
curved : bool
    Set True if phot_eleace is curved.
)"""
  );
  m.def(
      "has_orientation_attributes",
      &Bmad::has_orientation_attributes,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function has_orientation_attributes (ele) result (has_attribs)

Routine to determine whether an element has orientation attributes like x_offset, etc.
Also see: has_attribute function.

Parameters
----------
ele : EleStruct
    Lattice element.

Returns
-------
has_attribs : bool
    True if ele has orientation attributes. False otherwise.
)"""
  );
  m.def(
      "hdf5_write_beam",
      &Bmad::hdf5_write_beam,
      py::arg("file_name"),
      py::arg("bunches"),
      py::arg("append"),
      py::arg("error"),
      py::arg("lat") = py::none(),
      py::arg("alive_only") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine hdf5_write_beam

Parameters
----------
file_name : str

bunches : 1D array of BunchStruct

append : bool

error : bool

lat : LatStruct, optional

alive_only : bool, optional
)"""
  );
  m.def(
      "hdf5_write_grid_field",
      &Bmad::hdf5_write_grid_field,
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("g_field"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine hdf5_write_grid_field

Parameters
----------
file_name : str

ele : EleStruct

g_field : 1D array of GridFieldStruct

err_flag : bool
)"""
  );
  m.def(
      "hwang_bend_edge_kick",
      &Bmad::hwang_bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orb"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine hwang_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

Subroutine to track through the edge field of an sbend using a 2nd order map.
Adapted from:
  Hwang and S. Y. Lee,
  "Dipole Fringe Field Thin Map for Compact Synchrotrons",
  Phys. Rev. ST Accel. Beams, 12, 122401, (2015).
See the Bmad manual for details.

Parameters
----------
ele : EleStruct
    SBend element.

param : LatParamStruct
    Rel charge.

particle_at : int
    first_track_edge$, or second_track_edge$

orb : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place.
    As an output, orb: Coords after tracking.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the edge.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix including the edge.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.
)"""
  );
}
