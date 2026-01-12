#include "pybmad/generated/Bmad_routines_f.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyFringeHere python_fringe_here(
    EleProxy& ele,
    CoordProxy& orbit,
    int particle_at,
    bool is_here) {
  Bmad::fringe_here(ele, orbit, particle_at, is_here);
  auto py_result{PyFringeHere{is_here}};
  return py_result;
}

void init_Bmad_routines_f(py::module& m) {
  m.def(
      "fft1",
      &Bmad::fft1,
      py::arg("a"),
      py::arg("b"),
      py::arg("n"),
      py::arg("isn"),
      R"""(Parameters
  ----------
  a : 
  b : 
  n : 
  isn : 
  ierr : 
  )""");
  m.def(
      "field_attribute_free",
      &Bmad::field_attribute_free,
      py::arg("ele"),
      py::arg("attrib_name"),
      R"""(Function field_attribute_free (ele, attrib_name) result (free)

  Routine to check if a field attribute is free to vary.

  Field attributes are either normalized (EG K2 of a sextupole) or unnormalized (EG B2_GRADIENT of a sextupole).
  Whether normalized or unnormalized attributes are free to vary will depend on the setting  of ele%field_master.

  Generally, this routine should not be called directly. Use the routine attribute_free instead.

  Parameters
  ----------
  ele : EleStruct
      Element containing the attribute
  attrib_name : unknown
      Name of the field attribute. Assumed upper case.

  Returns
  -------
  free : bool
      Is the attribute free to vary? If the attribute is not recognized, free = True will be returned.
  )""");
  m.def(
      "finalize_reflectivity_table",
      &Bmad::finalize_reflectivity_table,
      py::arg("table"),
      py::arg("in_degrees"),
      R"""(Subroutine finalize_reflectivity_table (table, in_degrees)

  Routine to finalize the construction of the reflectivity tables for a surface.

  Parameters
  ----------
  table : PhotonReflectTableStruct
      Surface tables to be finalized.
      This parameter is an input/output and is modified in-place. As an output: Finalized surface tables.
  in_degrees : bool
      Table angles in degrees?
  )""");
  py::class_<Bmad::FindElementEnds, std::unique_ptr<Bmad::FindElementEnds>>(
      m, "FindElementEnds", "Fortran routine find_element_ends return value")
      .def_readonly("ele1", &Bmad::FindElementEnds::ele1)
      .def_readonly("ele2", &Bmad::FindElementEnds::ele2)
      .def("__len__", [](const Bmad::FindElementEnds&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::FindElementEnds& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.ele1);
            if (i == 1)
              return py::cast(s.ele2);
            throw py::index_error();
          });
  m.def(
      "find_element_ends",
      &Bmad::find_element_ends,
      py::arg("ele"),
      py::arg("ix_multipass") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to find the ends for.
  ele1 : EleStruct
      Pointer to element just before ele.
  ele2 : EleStruct
      Pointer to ele itself or the last sub-element within ele. Note: ele1 and ele2 will be nullified if ele is
      in the lord part of the lattice and does not have any slaves. Note: For an element in the tracking part of
      the lattice: ele1.ix_ele = ele.ix_ele - 1 ele2        => ele Exception: For Beginning element (index 0),
      ele1 => ele
  ix_multipass : int, optional
      Which multipass pass to follow. Default is 1. This is ignored if there is no multipass elements.
  )""");
  m.def(
      "find_fwhm",
      &Bmad::find_fwhm,
      py::arg("bound"),
      py::arg("args"),
      R"""(Subroutine find_fwhm(bound,args,fwhm)

  Finds the full width at half max of psi(t).  fwhm * c_light / TwoRtTwoLnTwo is taken as the bunch length.

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
  args : float
      Parameters and constants of dpsi/dt.  See comments of psi_prime for details.

  Returns
  -------
  fwhm : float
      Full width at half max of psi(t)
  )""");
  py::class_<
      Bmad::FindMatchingFieldmap,
      std::unique_ptr<Bmad::FindMatchingFieldmap>>(
      m,
      "FindMatchingFieldmap",
      "Fortran routine find_matching_fieldmap return value")
      .def_readonly("match_ele", &Bmad::FindMatchingFieldmap::match_ele)
      .def_readonly("ix_field", &Bmad::FindMatchingFieldmap::ix_field)
      .def("__len__", [](const Bmad::FindMatchingFieldmap&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::FindMatchingFieldmap& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.match_ele);
            if (i == 1)
              return py::cast(s.ix_field);
            throw py::index_error();
          });
  m.def(
      "find_matching_fieldmap",
      &Bmad::find_matching_fieldmap,
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("fm_type"),
      py::arg("ignore_slaves") = py::none(),
      R"""(Parameters
  ----------
  file_name : unknown
      File name associated with field to match to.
  ele : EleStruct
      Element holding the field to be matched.
  fm_type : int
      Type of fieldmap: cartesian_map$, cylindircal_map$, or gen_grad_map$, grid_field$
  match_ele : EleStruct
      Pointer to element with matched field. Nullified if no match found.
  ix_field : int
      index of field. For example: matching field => match_ele.cartesian_map(ix_field) Set to -1 if no match
      found.
  ignore_slaves : bool, optional
      If True, ignore any multipass slaves. Default is False.
  )""");
  m.def(
      "find_normalization",
      &Bmad::find_normalization,
      py::arg("bound"),
      py::arg("p0"),
      py::arg("args"),
      R"""(Subroutine find_normalization(bound,p0,args,pnrml)

  Finds value for boundary condition psi(0) that results in integral
  of psi(t) from -bound to +bound to be 1.0.  This is done with the secant method.
  Repeadedly calls integrate_psi with different values for psi(0).

  Parameters
  ----------
  bound : float
      -bound and +bound are integration boundaries
  p0 : float
      Boundary condition psi(0)
  args : float
      Parameters and constants of DEQ.  See psi_prime comments for details.

  Returns
  -------
  pnrml : float
      Value for psi(0) that results in integral of psi(t) from -bound to +bound being equal to 1.0
  )""");
  py::class_<Bmad::FloorAnglesToWMat, std::unique_ptr<Bmad::FloorAnglesToWMat>>(
      m,
      "FloorAnglesToWMat",
      "Fortran routine floor_angles_to_w_mat return value")
      .def_readonly("w_mat", &Bmad::FloorAnglesToWMat::w_mat)
      .def_readonly("w_mat_inv", &Bmad::FloorAnglesToWMat::w_mat_inv)
      .def("__len__", [](const Bmad::FloorAnglesToWMat&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::FloorAnglesToWMat& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.w_mat);
            if (i == 1)
              return py::cast(s.w_mat_inv);
            throw py::index_error();
          });
  m.def(
      "floor_angles_to_w_mat",
      &Bmad::floor_angles_to_w_mat,
      py::arg("theta"),
      py::arg("phi"),
      py::arg("psi"),
      R"""(Parameters
  ----------
  theta : float
      Azimuth angle.
  phi : float
      Pitch angle.
  psi : float
      Roll angle.
  w_mat : float
      Orientation matrix.
  w_mat_inv : float
      Inverse Orientation matrix.
  )""");
  py::class_<Bmad::FloorWMatToAngles, std::unique_ptr<Bmad::FloorWMatToAngles>>(
      m,
      "FloorWMatToAngles",
      "Fortran routine floor_w_mat_to_angles return value")
      .def_readonly("theta", &Bmad::FloorWMatToAngles::theta)
      .def_readonly("phi", &Bmad::FloorWMatToAngles::phi)
      .def_readonly("psi", &Bmad::FloorWMatToAngles::psi)
      .def("__len__", [](const Bmad::FloorWMatToAngles&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::FloorWMatToAngles& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.theta);
            if (i == 1)
              return py::cast(s.phi);
            if (i == 2)
              return py::cast(s.psi);
            throw py::index_error();
          });
  m.def(
      "floor_w_mat_to_angles",
      &Bmad::floor_w_mat_to_angles,
      py::arg("w_mat"),
      py::arg("floor0") = py::none(),
      R"""(Parameters
  ----------
  w_mat : float
      Orientation matrix.
  theta : float
      Azimuth angle.
  phi : float
      Pitch angle.
  psi : float
      Roll angle.
  floor0 : FloorPositionStruct, optional
      There are two solutions related by: [theta, phi, psi] & [pi+theta, pi-phi, pi+psi] If floor0 is present,
      choose the solution "nearest" the angles in floor0.
  )""");
  m.def(
      "form_complex_taylor",
      &Bmad::form_complex_taylor,
      py::arg("re_taylor"),
      py::arg("im_taylor"),
      R"""(Subroutine form_complex_taylor (re_taylor, im_taylor, complex_taylor)

  Subroutine to form a complex taylor from two taylor series representing
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
  )""");
  py::class_<
      Bmad::FormDigestedBmadFileName,
      std::unique_ptr<Bmad::FormDigestedBmadFileName>>(
      m,
      "FormDigestedBmadFileName",
      "Fortran routine form_digested_bmad_file_name return value")
      .def_readonly(
          "digested_file", &Bmad::FormDigestedBmadFileName::digested_file)
      .def_readonly(
          "full_lat_file", &Bmad::FormDigestedBmadFileName::full_lat_file)
      .def("__len__", [](const Bmad::FormDigestedBmadFileName&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::FormDigestedBmadFileName& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.digested_file);
            if (i == 1)
              return py::cast(s.full_lat_file);
            throw py::index_error();
          });
  m.def(
      "form_digested_bmad_file_name",
      &Bmad::form_digested_bmad_file_name,
      py::arg("lat_file"),
      py::arg("use_line") = py::none(),
      R"""(Subroutine form_digested_bmad_file_name (lat_file, digested_file, full_lat_file, use_line)

  Subroutine to form the standard name of the Bmad digested file.
  The standard digested file name has the suffix added to the file name:
      suffix = '.digested' + bmad_inc_version$
  Exception: If the use_line argument is present and not blank, the suffix will be:
      suffix = '.' + use_line + '.digested' + bmad_inc_version$


  Parameters
  ----------
  lat_file : unknown
      Input lattice file name.
  use_line : unknown, optional
      Line used for lattice expansion. If not present or blank, the line used is the one that was specified in
      the lattice file.

  Returns
  -------
  digested_file : unknown
      Name of the digested file.
  full_lat_file : unknown
      Input lattice file name with full directory. Can be used for error messages.
  )""");
  py::class_<PyFringeHere, std::unique_ptr<PyFringeHere>>(
      m, "FringeHere", "Fortran routine fringe_here return value")
      .def_readonly("is_here", &PyFringeHere::is_here)
      .def("__len__", [](const PyFringeHere&) { return 1; })
      .def("__getitem__", [](const PyFringeHere& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_here);
        throw py::index_error();
      });
  m.def(
      "fringe_here",
      &python_fringe_here,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("particle_at"),
      py::arg("is_here"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  orbit : CoordStruct
      Particle position.
  particle_at : int
      Either first_track_edge$ or second_track_edge$.
  is_here : 
  )""");
}
