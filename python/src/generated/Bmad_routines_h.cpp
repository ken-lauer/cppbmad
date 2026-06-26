#include "pybmad/generated/Bmad_routines_h.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_h(nb::module_ &m) {
  m.def(
      "hard_multipole_edge_kick",
      &Bmad::hard_multipole_edge_kick,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("particle_at"),
      nb::arg("orbit"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Routine to track through the hard edge field of a multipole.
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
      nb::arg("ele"),
      nb::arg("attrib"),
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
      nb::arg("phot_ele"),
      R"""(Routine to determine if a surface is potentially curved or is flat.

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
      nb::arg("ele"),
      R"""(Routine to determine whether an element has orientation attributes like x_offset, etc.
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
      "hdf5_read_beam",
      [](std::string file_name,
         BeamStruct &beam,
         bool error,
         EleStruct *ele,
         PmdHeaderStruct *pmd_header,
         std::optional<bool> print_mom_shift_warning,
         std::optional<bool> conserve_momentum) {
        auto fn = static_cast<
            void (*)(std::string, BeamStruct &, bool, optional_ref<EleStruct>, optional_ref<PmdHeaderStruct>, std::optional<bool>, std::optional<bool>)>(
            &Bmad::hdf5_read_beam
        );
        return fn(
            file_name,
            beam,
            error,
            ptr_to_opt_ref(ele),
            ptr_to_opt_ref(pmd_header),
            print_mom_shift_warning,
            conserve_momentum
        );
      },
      nb::arg("file_name"),
      nb::arg("beam"),
      nb::arg("error"),
      nb::arg("ele") = nb::none(),
      nb::arg("pmd_header") = nb::none(),
      nb::arg("print_mom_shift_warning") = nb::none(),
      nb::arg("conserve_momentum") = nb::none(),
      R"""(Wrapper for Fortran routine hdf5_read_beam

Parameters
----------
file_name : str

beam : BeamStruct

error : bool

ele : EleStruct, optional

pmd_header : PmdHeaderStruct, optional

print_mom_shift_warning : bool, optional

conserve_momentum : bool, optional
)"""
  );
  m.def(
      "hdf5_read_grid_field",
      [](std::string file_name,
         EleStruct &ele,
         GridFieldStructArray1D g_field,
         bool err_flag,
         PmdHeaderStruct *pmd_header,
         std::optional<bool> combine) {
        auto fn = static_cast<
            void (*)(std::string, EleStruct &, GridFieldStructArray1D, bool, optional_ref<PmdHeaderStruct>, std::optional<bool>)>(
            &Bmad::hdf5_read_grid_field
        );
        return fn(file_name, ele, g_field, err_flag, ptr_to_opt_ref(pmd_header), combine);
      },
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("g_field"),
      nb::arg("err_flag"),
      nb::arg("pmd_header") = nb::none(),
      nb::arg("combine") = nb::none(),
      R"""(Wrapper for Fortran routine hdf5_read_grid_field

Parameters
----------
file_name : str

ele : EleStruct

g_field : 1D array of GridFieldStruct

err_flag : bool

pmd_header : PmdHeaderStruct, optional

combine : bool, optional
)"""
  );
  m.def(
      "hdf5_write_beam",
      [](std::string file_name,
         BunchStructArray1D bunches,
         bool append,
         bool error,
         LatStruct *lat,
         std::optional<bool> alive_only) {
        auto fn = static_cast<
            void (*)(std::string, BunchStructArray1D, bool, bool, optional_ref<LatStruct>, std::optional<bool>)>(
            &Bmad::hdf5_write_beam
        );
        return fn(file_name, bunches, append, error, ptr_to_opt_ref(lat), alive_only);
      },
      nb::arg("file_name"),
      nb::arg("bunches"),
      nb::arg("append"),
      nb::arg("error"),
      nb::arg("lat") = nb::none(),
      nb::arg("alive_only") = nb::none(),
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
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("g_field"),
      nb::arg("err_flag"),
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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("particle_at"),
      nb::arg("orb"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Subroutine to track through the edge field of an sbend using a 2nd order map.
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
