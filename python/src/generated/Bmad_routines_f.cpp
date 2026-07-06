#include "pybmad/generated/Bmad_routines_f.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyFibreToEle python_fibre_to_ele(
    Fibre &ptc_fibre,
    BranchStruct &branch,
    int ix_ele,
    std::optional<bool> from_mad = std::nullopt
) {
  auto _result = Bmad::fibre_to_ele(ptc_fibre, branch, ix_ele, from_mad);
  auto py_result{PyFibreToEle{_result, ix_ele}};
  return py_result;
}

void init_Bmad_routines_f(nb::module_ &m) {
  m.def(
      "fft1",
      &Bmad::fft1,
      nb::arg("a"),
      nb::arg("b"),
      nb::arg("n"),
      nb::arg("isn"),
      R"""(Wrapper for Fortran routine fft1

Parameters
----------
a : 1D array of float

b : 1D array of float

n : int

isn : int

Returns
-------
ierr : int
)"""
  );
  nb::class_<PyFibreToEle>(m, "FibreToEle", "fibre_to_ele return type")
      .def_ro("err_flag", &PyFibreToEle::err_flag)
      .def_ro("ix_ele", &PyFibreToEle::ix_ele)
      .def("__len__", [](const PyFibreToEle &) { return 2; })
      .def("__getitem__", [](const PyFibreToEle &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.ix_ele);
        throw nb::index_error();
      });
  m.def(
      "fibre_to_ele",
      &python_fibre_to_ele,
      nb::arg("ptc_fibre"),
      nb::arg("branch"),
      nb::arg("ix_ele"),
      nb::arg("from_mad") = nb::none(),
      R"""(Wrapper for Fortran routine fibre_to_ele

Parameters
----------
ptc_fibre : Fibre
    PTC fibre.

branch : BranchStruct
    branch containing elements.

ix_ele : int
    Index in ele(:) array of element last used.
    This parameter is an input/output and is modified in-place.
    As an output, ix_ele: Index to element created (upper index if more than one created).

from_mad : bool, optional
    If True, ignore PTC specific parameters like integrator_order. Default is False. True is used when the
    fibre has been created via MAD. In this case, the PTC specific parameters may not have good values.

Returns
-------
ix_ele : int
    Index in ele(:) array of element last used.
    This parameter is an input/output and is modified in-place.
    As an output, ix_ele: Index to element created (upper index if more than one created).

err_flag : bool
    Set true if there is an error. False otherwise. To do: lcavity energy change !? open or closed geometry?
    Energy patch
)"""
  );
  m.def(
      "field_attribute_free",
      &Bmad::field_attribute_free,
      nb::arg("ele"),
      nb::arg("attrib_name"),
      R"""(Routine to check if a field attribute is free to vary.

Field attributes are either normalized (EG K2 of a sextupole) or unnormalized (EG B2_GRADIENT of a sextupole).
Whether normalized or unnormalized attributes are free to vary will depend on the setting  of ele%field_master.

Generally, this routine should not be called directly. Use the routine attribute_free instead.

Parameters
----------
ele : EleStruct
    Element containing the attribute

attrib_name : str
    Name of the field attribute. Assumed upper case.

Returns
-------
free : bool
    Is the attribute free to vary? If the attribute is not recognized, free = True will be returned.
)"""
  );
  m.def(
      "finalize_reflectivity_table",
      &Bmad::finalize_reflectivity_table,
      nb::arg("table"),
      nb::arg("in_degrees"),
      R"""(Routine to finalize the construction of the reflectivity tables for a surface.

Parameters
----------
table : PhotonReflectTableStruct
    Surface tables to be finalized.
    This parameter is an input/output and is modified in-place.
    As an output, table: Finalized surface tables.

in_degrees : bool
    Table angles in degrees?
)"""
  );
  nb::class_<Bmad::FindElementEnds>(m, "FindElementEnds", "find_element_ends return type")
      .def_ro("ele1", &Bmad::FindElementEnds::ele1)
      .def_ro("ele2", &Bmad::FindElementEnds::ele2)
      .def("__len__", [](const Bmad::FindElementEnds &) { return 2; })
      .def("__getitem__", [](const Bmad::FindElementEnds &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ele1);
        if (i == 1)
          return nb::cast(s.ele2);
        throw nb::index_error();
      });
  m.def(
      "find_element_ends",
      &Bmad::find_element_ends,
      nb::arg("ele"),
      nb::arg("ix_multipass") = nb::none(),
      R"""(Wrapper for Fortran routine find_element_ends

Parameters
----------
ele : EleStruct
    Element to find the ends for.

ix_multipass : int, optional
    Which multipass pass to follow. Default is 1. This is ignored if there is no multipass elements.

Returns
-------
ele1 : EleStruct, optional
    Pointer to element just before ele.

ele2 : EleStruct, optional
    Pointer to ele itself or the last sub-element within ele. Note: ele1 and ele2 will be nullified if ele is
    in the lord part of the lattice and does not have any slaves. Note: For an element in the tracking part of
    the lattice: ele1.ix_ele = ele.ix_ele - 1 ele2        => ele Exception: For Beginning element (index 0),
    ele1 => ele
)"""
  );
  m.def(
      "find_fwhm",
      &Bmad::find_fwhm,
      nb::arg("bound"),
      nb::arg("args"),
      R"""(Finds the full width at half max of psi(t).  fwhm * c_light / TwoRtTwoLnTwo is taken as the bunch length.

Steps followed:
  Find value for p(0) that normalizes the solution to dpsi/dt.
  Find max value of p(t) for the value of p(0) found in the previous step.
  Find find tlower, tlower < 0, such that p(tlower) = pmax/2.
  Find find tupper, tupper > 0, such that p(tupper) = pmax/2.
  fwhm is tupper-tlower

Parameters
----------
bound : float
    -bound and +bound is integration bound.

args : 1D array of float (shape: 1:8)
    Parameters and constants of dpsi/dt.  See comments of psi_prime for details.

Returns
-------
fwhm : float
    Full width at half max of psi(t)
)"""
  );
  nb::class_<Bmad::FindMatchingFieldmap>(
      m,
      "FindMatchingFieldmap",
      "find_matching_fieldmap return type"
  )
      .def_ro("match_ele", &Bmad::FindMatchingFieldmap::match_ele)
      .def_ro("ix_field", &Bmad::FindMatchingFieldmap::ix_field)
      .def("__len__", [](const Bmad::FindMatchingFieldmap &) { return 2; })
      .def("__getitem__", [](const Bmad::FindMatchingFieldmap &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.match_ele);
        if (i == 1)
          return nb::cast(s.ix_field);
        throw nb::index_error();
      });
  m.def(
      "find_matching_fieldmap",
      &Bmad::find_matching_fieldmap,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("fm_type"),
      nb::arg("ignore_slaves") = nb::none(),
      R"""(Wrapper for Fortran routine find_matching_fieldmap

Parameters
----------
file_name : str
    File name associated with field to match to.

ele : EleStruct
    Element holding the field to be matched.

fm_type : int
    Type of fieldmap: cartesian_map$, cylindircal_map$, or gen_grad_map$, grid_field$

ignore_slaves : bool, optional
    If True, ignore any multipass slaves. Default is False.

Returns
-------
ix_field : int
    index of field. For example: matching field => match_ele.cartesian_map(ix_field) Set to -1 if no match
    found.

match_ele : EleStruct, optional
    Pointer to element with matched field. Nullified if no match found.
)"""
  );
  m.def(
      "find_normalization",
      &Bmad::find_normalization,
      nb::arg("bound"),
      nb::arg("p0"),
      nb::arg("args"),
      R"""(Finds value for boundary condition psi(0) that results in integral
of psi(t) from -bound to +bound to be 1.0.  This is done with the secant method.
Repeadedly calls integrate_psi with different values for psi(0).

Parameters
----------
bound : float
    -bound and +bound are integration boundaries

p0 : float
    Boundary condition psi(0)

args : 1D array of float (shape: 1:8)
    Parameters and constants of DEQ.  See psi_prime comments for details.

Returns
-------
pnrml : float
    Value for psi(0) that results in integral of psi(t) from -bound to +bound being equal to 1.0
)"""
  );
  nb::class_<Bmad::FloorAnglesToWMat>(m, "FloorAnglesToWMat", "floor_angles_to_w_mat return type")
      .def_ro("w_mat", &Bmad::FloorAnglesToWMat::w_mat)
      .def_ro("w_mat_inv", &Bmad::FloorAnglesToWMat::w_mat_inv)
      .def("__len__", [](const Bmad::FloorAnglesToWMat &) { return 2; })
      .def("__getitem__", [](const Bmad::FloorAnglesToWMat &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.w_mat);
        if (i == 1)
          return nb::cast(s.w_mat_inv);
        throw nb::index_error();
      });
  m.def(
      "floor_angles_to_w_mat",
      &Bmad::floor_angles_to_w_mat,
      nb::arg("theta"),
      nb::arg("phi"),
      nb::arg("psi"),
      R"""(Wrapper for Fortran routine floor_angles_to_w_mat

Parameters
----------
theta : float
    Azimuth angle.

phi : float
    Pitch angle.

psi : float
    Roll angle.

Returns
-------
w_mat : 2D array of float (shape: 3,3), optional
    Orientation matrix.

w_mat_inv : 2D array of float (shape: 3,3), optional
    Inverse Orientation matrix.
)"""
  );
  nb::class_<Bmad::FloorWMatToAngles>(m, "FloorWMatToAngles", "floor_w_mat_to_angles return type")
      .def_ro("theta", &Bmad::FloorWMatToAngles::theta)
      .def_ro("phi", &Bmad::FloorWMatToAngles::phi)
      .def_ro("psi", &Bmad::FloorWMatToAngles::psi)
      .def("__len__", [](const Bmad::FloorWMatToAngles &) { return 3; })
      .def("__getitem__", [](const Bmad::FloorWMatToAngles &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.theta);
        if (i == 1)
          return nb::cast(s.phi);
        if (i == 2)
          return nb::cast(s.psi);
        throw nb::index_error();
      });
  m.def(
      "floor_w_mat_to_angles",
      [](FixedArray2D<Real, 3, 3> w_mat, FloorPositionStruct *floor0) {
        auto fn = static_cast<Bmad::FloorWMatToAngles (*)(
            FixedArray2D<Real, 3, 3>,
            optional_ref<FloorPositionStruct>
        )>(&Bmad::floor_w_mat_to_angles);
        return fn(w_mat, ptr_to_opt_ref(floor0));
      },
      nb::arg("w_mat"),
      nb::arg("floor0") = nb::none(),
      R"""(Wrapper for Fortran routine floor_w_mat_to_angles

Parameters
----------
w_mat : 2D array of float (shape: 3,3)
    Orientation matrix.

floor0 : FloorPositionStruct, optional
    There are two solutions related by: [theta, phi, psi] & [pi+theta, pi-phi, pi+psi] If floor0 is present,
    choose the solution "nearest" the angles in floor0.

Returns
-------
theta : float
    Azimuth angle.

phi : float
    Pitch angle.

psi : float
    Roll angle.
)"""
  );
  m.def(
      "form_complex_taylor",
      &Bmad::form_complex_taylor,
      nb::arg("re_taylor"),
      nb::arg("im_taylor"),
      R"""(Subroutine to form a complex taylor from two taylor series representing
  the real and imaginary parts

Parameters
----------
re_taylor : TaylorStruct
    Real part

im_taylor : TaylorStruct
    Imaginary part

Returns
-------
complex_taylor : ComplexTaylorStruct
    combined complex taylor
)"""
  );
  nb::class_<Bmad::FormDigestedBmadFileName>(
      m,
      "FormDigestedBmadFileName",
      "form_digested_bmad_file_name return type"
  )
      .def_ro("digested_file", &Bmad::FormDigestedBmadFileName::digested_file)
      .def_ro("full_lat_file", &Bmad::FormDigestedBmadFileName::full_lat_file)
      .def("__len__", [](const Bmad::FormDigestedBmadFileName &) { return 2; })
      .def("__getitem__", [](const Bmad::FormDigestedBmadFileName &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.digested_file);
        if (i == 1)
          return nb::cast(s.full_lat_file);
        throw nb::index_error();
      });
  m.def(
      "form_digested_bmad_file_name",
      &Bmad::form_digested_bmad_file_name,
      nb::arg("lat_file"),
      nb::arg("use_line") = nb::none(),
      R"""(Subroutine to form the standard name of the Bmad digested file.
The standard digested file name has the suffix added to the file name:
    suffix = '.digested' + bmad_inc_version$
Exception: If the use_line argument is present and not blank, the suffix will be:
    suffix = '.' + use_line + '.digested' + bmad_inc_version$

Parameters
----------
lat_file : str
    Input lattice file name.

use_line : str, optional
    Line used for lattice expansion. If not present or blank, the line used is the one that was specified in
    the lattice file.

Returns
-------
digested_file : str
    Name of the digested file.

full_lat_file : str, optional
    Input lattice file name with full directory. Can be used for error messages.
)"""
  );
  m.def(
      "fringe_here",
      &Bmad::fringe_here,
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("particle_at"),
      R"""(Wrapper for Fortran routine fringe_here

Parameters
----------
ele : EleStruct
    Lattice element.

orbit : CoordStruct
    Particle position.

particle_at : int
    Either first_track_edge$ or second_track_edge$.

Returns
-------
is_here : bool
    True if there is a fringe. False if not.
)"""
  );
}
