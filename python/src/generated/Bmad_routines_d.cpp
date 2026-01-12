#include "pybmad/generated/Bmad_routines_d.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyDampingMatrixD {
  double gamma;
  double g_tot;
  double B0;
  double B1;
  double delta;
  int species;
};
PyDampingMatrixD python_damping_matrix_d(
    double gamma,
    double g_tot,
    double B0,
    double B1,
    double delta,
    int species,
    FixedArray2D<Real, 6, 6> mat) {
  Bmad::damping_matrix_d(gamma, g_tot, B0, B1, delta, species, mat);
  auto py_result{PyDampingMatrixD{gamma, g_tot, B0, B1, delta, species}};
  return py_result;
}
struct PyDefaultTrackingSpecies {
  int species;
};
PyDefaultTrackingSpecies python_default_tracking_species(
    LatParamProxy& param,
    int species) {
  Bmad::default_tracking_species(param, species);
  auto py_result{PyDefaultTrackingSpecies{species}};
  return py_result;
}
struct PyDiffractionPlateOrMaskHitSpot {
  int ix_section;
};
PyDiffractionPlateOrMaskHitSpot python_diffraction_plate_or_mask_hit_spot(
    EleProxy& ele,
    CoordProxy& orbit,
    int ix_section) {
  Bmad::diffraction_plate_or_mask_hit_spot(ele, orbit, ix_section);
  auto py_result{PyDiffractionPlateOrMaskHitSpot{ix_section}};
  return py_result;
}
struct PyDiffusionMatrixB {
  double gamma;
  double g_tot;
  int species;
};
PyDiffusionMatrixB python_diffusion_matrix_b(
    double gamma,
    double g_tot,
    int species,
    FixedArray2D<Real, 6, 6> mat) {
  Bmad::diffusion_matrix_b(gamma, g_tot, species, mat);
  auto py_result{PyDiffusionMatrixB{gamma, g_tot, species}};
  return py_result;
}

struct PyDistanceToAperture {
  bool no_aperture_here;
  double dist;
};
PyDistanceToAperture python_distance_to_aperture(
    CoordProxy& orbit,
    int particle_at,
    EleProxy& ele,
    double dist) {
  auto _result = Bmad::distance_to_aperture(orbit, particle_at, ele, dist);
  auto py_result{PyDistanceToAperture{_result, dist}};
  return py_result;
}
struct PyDpcGivenDe {
  double pc_old;
  double mass;
  double dE;
  double dpc;
};
PyDpcGivenDe python_dpc_given_de(
    double pc_old,
    double mass,
    double dE,
    double dpc) {
  Bmad::dpc_given_de(pc_old, mass, dE, dpc);
  auto py_result{PyDpcGivenDe{pc_old, mass, dE, dpc}};
  return py_result;
}

void init_Bmad_routines_d(py::module& m) {
  m.def(
      "damping_matrix_d",
      &python_damping_matrix_d,
      py::arg("gamma"),
      py::arg("g_tot"),
      py::arg("B0"),
      py::arg("B1"),
      py::arg("delta"),
      py::arg("species"),
      py::arg("mat"),
      R"""(Parameters
  ----------
  gamma : 
  g_tot : 
  B0 : 
  B1 : 
  delta : 
  species : 
  mat : 
  )""");
  py::class_<PyDampingMatrixD, std::unique_ptr<PyDampingMatrixD>>(
      m, "DampingMatrixD", "Fortran routine damping_matrix_d return value")
      .def_readonly("gamma", &PyDampingMatrixD::gamma)
      .def_readonly("g_tot", &PyDampingMatrixD::g_tot)
      .def_readonly("B0", &PyDampingMatrixD::B0)
      .def_readonly("B1", &PyDampingMatrixD::B1)
      .def_readonly("delta", &PyDampingMatrixD::delta)
      .def_readonly("species", &PyDampingMatrixD::species)
      .def("__len__", [](const PyDampingMatrixD&) { return 6; })
      .def("__getitem__", [](const PyDampingMatrixD& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.gamma);
        if (i == 1)
          return py::cast(s.g_tot);
        if (i == 2)
          return py::cast(s.B0);
        if (i == 3)
          return py::cast(s.B1);
        if (i == 4)
          return py::cast(s.delta);
        if (i == 5)
          return py::cast(s.species);
        throw py::index_error();
      });
  m.def(
      "deallocate_ele_pointers",
      &Bmad::deallocate_ele_pointers,
      py::arg("ele"),
      py::arg("nullify_only") = py::none(),
      py::arg("nullify_branch") = py::none(),
      py::arg("dealloc_poles") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element with pointers.
      This parameter is an input/output and is modified in-place. As an output: Element with deallocated
      pointers.
  nullify_only : bool, optional
      If present and True: Nullify & do not deallocate.
  nullify_branch : bool, optional
      Nullify ele.branch? Default is True.
  dealloc_poles : bool, optional
      Dealloc ele.a/b_pole, ele.a/b_pole_elec? Default is True.
  )""");
  m.def(
      "deallocate_expression_tree",
      &Bmad::deallocate_expression_tree,
      py::arg("tree"),
      R"""(Subroutine deallocate_expression_tree(tree)

  Routine to deallocate an expression tree.

  Parameters
  ----------
  tree : ExpressionTreeStruct
      Tree to deallocate.
      This parameter is an input/output and is modified in-place. As an output: Deallocated tree.
  )""");
  m.def(
      "deallocate_lat_pointers",
      &Bmad::deallocate_lat_pointers,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Lat with pointers.
      This parameter is an input/output and is modified in-place. As an output: Lat with deallocated pointers.
  )""");
  m.def(
      "default_tracking_species",
      &python_default_tracking_species,
      py::arg("param"),
      py::arg("species"),
      R"""(Parameters
  ----------
  param : LatParamStruct
      Parameters for a lattice branch.
  species : 
  )""");
  py::class_<
      PyDefaultTrackingSpecies,
      std::unique_ptr<PyDefaultTrackingSpecies>>(
      m,
      "DefaultTrackingSpecies",
      "Fortran routine default_tracking_species return value")
      .def_readonly("species", &PyDefaultTrackingSpecies::species)
      .def("__len__", [](const PyDefaultTrackingSpecies&) { return 1; })
      .def(
          "__getitem__",
          [](const PyDefaultTrackingSpecies& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.species);
            throw py::index_error();
          });
  m.def(
      "detector_pixel_pt",
      &Bmad::detector_pixel_pt,
      py::arg("orbit"),
      py::arg("ele"),
      R"""(Function detector_pixel_pt (orbit, ele) result (ix_pix)

  Routine to return the pixel a particle is hitting.

  Parameters
  ----------
  orbit : CoordStruct
      Orbit at surface.
  ele : EleStruct
      Detector element.

  Returns
  -------
  ix_pix : int
      index of ele.photon.pixel.pt(:,:) the particle is in.
  )""");
  m.def(
      "diffraction_plate_or_mask_hit_spot",
      &python_diffraction_plate_or_mask_hit_spot,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("ix_section"),
      R"""(Parameters
  ----------
  ele : EleStruct
      diffraction_plate or mask element.
  orbit : CoordStruct
      particle position.
  ix_section : 
  )""");
  py::class_<
      PyDiffractionPlateOrMaskHitSpot,
      std::unique_ptr<PyDiffractionPlateOrMaskHitSpot>>(
      m,
      "DiffractionPlateOrMaskHitSpot",
      "Fortran routine diffraction_plate_or_mask_hit_spot return value")
      .def_readonly("ix_section", &PyDiffractionPlateOrMaskHitSpot::ix_section)
      .def("__len__", [](const PyDiffractionPlateOrMaskHitSpot&) { return 1; })
      .def(
          "__getitem__",
          [](const PyDiffractionPlateOrMaskHitSpot& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.ix_section);
            throw py::index_error();
          });
  m.def(
      "diffusion_matrix_b",
      &python_diffusion_matrix_b,
      py::arg("gamma"),
      py::arg("g_tot"),
      py::arg("species"),
      py::arg("mat"),
      R"""(Parameters
  ----------
  gamma : 
  g_tot : 
  species : 
  mat : 
  )""");
  py::class_<PyDiffusionMatrixB, std::unique_ptr<PyDiffusionMatrixB>>(
      m, "DiffusionMatrixB", "Fortran routine diffusion_matrix_b return value")
      .def_readonly("gamma", &PyDiffusionMatrixB::gamma)
      .def_readonly("g_tot", &PyDiffusionMatrixB::g_tot)
      .def_readonly("species", &PyDiffusionMatrixB::species)
      .def("__len__", [](const PyDiffusionMatrixB&) { return 3; })
      .def("__getitem__", [](const PyDiffusionMatrixB& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.gamma);
        if (i == 1)
          return py::cast(s.g_tot);
        if (i == 2)
          return py::cast(s.species);
        throw py::index_error();
      });
  m.def(
      "distance_to_aperture",
      &python_distance_to_aperture,
      py::arg("orbit"),
      py::arg("particle_at"),
      py::arg("ele"),
      py::arg("dist"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle position.
  particle_at : int
      first_track_edge$, second_track_edge$, or in_between$
  ele : EleStruct
      Element containing aperture.
  no_aperture_here : bool
      True if aperture does not exist at the longitudinal location of the particle.
  dist : 
  )""");
  py::class_<PyDistanceToAperture, std::unique_ptr<PyDistanceToAperture>>(
      m,
      "DistanceToAperture",
      "Fortran routine distance_to_aperture return value")
      .def_readonly("no_aperture_here", &PyDistanceToAperture::no_aperture_here)
      .def_readonly("dist", &PyDistanceToAperture::dist)
      .def("__len__", [](const PyDistanceToAperture&) { return 2; })
      .def(
          "__getitem__",
          [](const PyDistanceToAperture& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.no_aperture_here);
            if (i == 1)
              return py::cast(s.dist);
            throw py::index_error();
          });
  m.def(
      "do_mode_flip",
      &Bmad::do_mode_flip,
      py::arg("ele"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Starting Element
      This parameter is an input/output and is modified in-place. As an output: Flipped element
  err_flag : bool
      Set True if there is an error. False otherwise.
  )""");
  m.def(
      "dpc_given_de",
      &python_dpc_given_de,
      py::arg("pc_old"),
      py::arg("mass"),
      py::arg("dE"),
      py::arg("dpc"),
      R"""(Parameters
  ----------
  pc_old : 
  mass : 
  dE : 
  dpc : 
  )""");
  py::class_<PyDpcGivenDe, std::unique_ptr<PyDpcGivenDe>>(
      m, "DpcGivenDe", "Fortran routine dpc_given_de return value")
      .def_readonly("pc_old", &PyDpcGivenDe::pc_old)
      .def_readonly("mass", &PyDpcGivenDe::mass)
      .def_readonly("dE", &PyDpcGivenDe::dE)
      .def_readonly("dpc", &PyDpcGivenDe::dpc)
      .def("__len__", [](const PyDpcGivenDe&) { return 4; })
      .def("__getitem__", [](const PyDpcGivenDe& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.pc_old);
        if (i == 1)
          return py::cast(s.mass);
        if (i == 2)
          return py::cast(s.dE);
        if (i == 3)
          return py::cast(s.dpc);
        throw py::index_error();
      });
  m.def(
      "drift_and_pipe_track_methods_adjustment",
      &Bmad::drift_and_pipe_track_methods_adjustment,
      py::arg("lat"),
      R"""(Subroutine drift_and_pipe_track_methods_adjustment(lat)

  Drift and pipe elements can be used in both photon and non-photon lines.
  A problem occures if, for example, a lattice file with both photon and
  non-photon lines contains a line like:
    drift::*[tracking_method] = taylor
  So this routine resets drift and pipe tracking_method and mat6_calc_method
  parameters in photon lines to bmad_standard if needed.

  Parameters
  ----------
  lat : LatStruct
      Lattice
      This parameter is an input/output and is modified in-place. As an output: Lattice with tracking methods
      adjusted if needed.
  )""");
  m.def(
      "drift_multipass_name_correction",
      &Bmad::drift_multipass_name_correction,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : 
  )""");
  m.def(
      "drift_orbit_time",
      &Bmad::drift_orbit_time,
      py::arg("orbit"),
      py::arg("beta0"),
      py::arg("delta_s") = py::none(),
      py::arg("delta_t") = py::none(),
      R"""(Subroutine drift_orbit_time(orbit, beta0, delta_s, delta_t)

  Simple routine to drift a particle orbit in time-based coordinates by a distance delta_s
    or a time delta_t
    If the particle has zero longitudinal velocity, then the particle is not drifted
    and a warning is printed.

  Parameters
  ----------
  orbit : CoordStruct
      particle orbit in time-based coordinates.
  beta0 : float
      reference velocity v/c.
  delta_s : float, optional
      s-coordinate distance to drift particle.
  delta_t : float, optional
      -coordinate distancet to drift particle.
  )""");
  m.def(
      "drift_particle_to_s",
      &Bmad::drift_particle_to_s,
      py::arg("p"),
      py::arg("s"),
      py::arg("branch"),
      R"""(Subroutine drift_particle_to_s (p, s, branch)

  Drift a particle to a given s-coordinate

  Parameters
  ----------
  p : CoordStruct
      Init particle position.
      This parameter is an input/output and is modified in-place. As an output: Final particle position.
  s : float
      Target s coordinate.
  branch : BranchStruct
      Branch being tracked through.
  )""");
  m.def(
      "drift_particle_to_t",
      &Bmad::drift_particle_to_t,
      py::arg("p"),
      py::arg("t"),
      py::arg("branch"),
      R"""(Subroutine drift_particle_to_t (p, t, branch)

  Drift a particle to a given t-coordinate

  Parameters
  ----------
  p : CoordStruct
      Init particle position.
      This parameter is an input/output and is modified in-place. As an output: Final particle position.
  t : float
      Target t coordinate.
  branch : BranchStruct
      Lattice branch being tracked through.
  )""");
  m.def(
      "dspline_len",
      &Bmad::dspline_len,
      py::arg("s_chord0"),
      py::arg("s_chord1"),
      py::arg("spline"),
      py::arg("dtheta_ref") = py::none(),
      R"""(Function dspline_len (s_chord0, s_chord1, spline, dtheta_ref) result (dlen)

  Routine to calculate the difference in length between the spline curve length and a referece line.
  Referece line is centroid chord (referece system of the spline) rotated by dtheta_ref.

  Parameters
  ----------
  s_chord0 : float
      Start position along centroid chord.
  s_chord1 : float
      Stop position along central_chord.
  spline : SplineStruct
      Spline of x-position as a function of s.
  dtheta_ref : float, optional
      angle to rotate the reference line from the centroid chord. Default is 0.

  Returns
  -------
  dlen : float
      L_spline - L_chord
  )""");
  m.def(
      "dynamic_aperture_point",
      &Bmad::dynamic_aperture_point,
      py::arg("branch"),
      py::arg("ele0"),
      py::arg("orb0"),
      py::arg("theta_xy"),
      py::arg("ap_param"),
      py::arg("check_xy_init") = py::none(),
      R"""(Subroutine dynamic_aperture_point (branch, ele0, orb0, theta_xy, ap_param, ap_point, check_xy_init)

  Subroutine to determine one dynamic aperture point by tracking.
  This routine works by determining where on a radial line y = const * x the aperture is.
  Here x and y are deviations from the reference orbit.

  Parameters
  ----------
  branch : BranchStruct
      Lattice branch to track through.
  ele0 : EleStruct
      Lattice element at start of tracking
  orb0 : CoordStruct
      reference orbit at the start of tracking.
  theta_xy : float
      Angle of radial line (in radians) in x-y space. Angle is "normalized" by .x_init, .y_init.
  ap_param : ApertureParamStruct
      Structure holding the input data:
  check_xy_init : bool, optional
      If True, do not check that aperture_param.x_init and .y_init are non-zero. Default is True.

  Returns
  -------
  ap_point : AperturePointStruct
  )""");
  m.def(
      "dynamic_aperture_scan",
      &Bmad::dynamic_aperture_scan,
      py::arg("aperture_param"),
      py::arg("pz_start"),
      py::arg("lat"),
      py::arg("print_timing") = py::none(),
      R"""(Subroutine dynamic_aperture_scan(aperture_scan, aperture_param, pz_start, lat, print_timing)

  Routine to do a set of dynamic aperture scans.

  Parameters
  ----------
  aperture_param : ApertureParamStruct
      Scan parameters
  lat : LatStruct
      Lattice.
  pz_start : float
      Starting phase space pz values.
  print_timing : bool, optional
      If True print info on calculation time. Default is True.

  Returns
  -------
  aperture_scan : ApertureScanStruct
      Set of scans. One for each pz_start(:).
  )""");
}
