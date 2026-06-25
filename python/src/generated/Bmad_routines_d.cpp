#include "pybmad/generated/Bmad_routines_d.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_d(py::module &m) {
  m.def(
      "damping_matrix_d",
      &Bmad::damping_matrix_d,
      py::arg("gamma"),
      py::arg("g_tot"),
      py::arg("B0"),
      py::arg("B1"),
      py::arg("delta"),
      py::arg("species"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine damping_matrix_d

Parameters
----------
gamma : float

g_tot : float

B0 : float

B1 : float

delta : float

species : int

Returns
-------
mat : 2D array of float (shape: 6,6)
)"""
  );
  m.def(
      "deallocate_ele_pointers",
      &Bmad::deallocate_ele_pointers,
      py::arg("ele"),
      py::arg("nullify_only") = py::none(),
      py::arg("nullify_branch") = py::none(),
      py::arg("dealloc_poles") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine deallocate_ele_pointers

Parameters
----------
ele : EleStruct
    Element with pointers.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with deallocated pointers.

nullify_only : bool, optional
    If present and True: Nullify & do not deallocate.

nullify_branch : bool, optional
    Nullify ele.branch? Default is True.

dealloc_poles : bool, optional
    Dealloc ele.a/b_pole, ele.a/b_pole_elec? Default is True.
)"""
  );
  m.def(
      "deallocate_expression_tree",
      &Bmad::deallocate_expression_tree,
      py::arg("tree"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine deallocate_expression_tree(tree)

Routine to deallocate an expression tree.

Parameters
----------
tree : ExpressionTreeStruct
    Tree to deallocate.
    This parameter is an input/output and is modified in-place.
    As an output, tree: Deallocated tree.
)"""
  );
  m.def(
      "deallocate_lat_pointers",
      &Bmad::deallocate_lat_pointers,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine deallocate_lat_pointers

Parameters
----------
lat : LatStruct
    Lat with pointers.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lat with deallocated pointers.
)"""
  );
  m.def(
      "default_tracking_species",
      &Bmad::default_tracking_species,
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine default_tracking_species

Parameters
----------
param : LatParamStruct
    Parameters for a lattice branch.

Returns
-------
species : int
    Default species to be used for tracking.
)"""
  );
  m.def(
      "detector_pixel_pt",
      &Bmad::detector_pixel_pt,
      py::arg("orbit"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
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
ix_pix : 1D array of int (shape: 2)
    index of ele.photon.pixel.pt(:,:) the particle is in.
)"""
  );
  m.def(
      "diffraction_plate_or_mask_hit_spot",
      &Bmad::diffraction_plate_or_mask_hit_spot,
      py::arg("ele"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine diffraction_plate_or_mask_hit_spot

Parameters
----------
ele : EleStruct
    diffraction_plate or mask element.

orbit : CoordStruct
    particle position.

Returns
-------
ix_section : int
    integer, Set to index of the section where the particle is. Set to zero if the photon is outside all clear
    areas.
)"""
  );
  m.def(
      "diffusion_matrix_b",
      &Bmad::diffusion_matrix_b,
      py::arg("gamma"),
      py::arg("g_tot"),
      py::arg("species"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine diffusion_matrix_b

Parameters
----------
gamma : float

g_tot : float

species : int

Returns
-------
mat : 2D array of float (shape: 6,6)
)"""
  );
  py::class_<Bmad::DistanceToAperture, std::unique_ptr<Bmad::DistanceToAperture>>(
      m,
      "DistanceToAperture",
      "distance_to_aperture return type"
  )
      .def_readonly("no_aperture_here", &Bmad::DistanceToAperture::no_aperture_here)
      .def_readonly("dist", &Bmad::DistanceToAperture::dist)
      .def("__len__", [](const Bmad::DistanceToAperture &) { return 2; })
      .def("__getitem__", [](const Bmad::DistanceToAperture &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.no_aperture_here);
        if (i == 1)
          return py::cast(s.dist);
        throw py::index_error();
      });
  m.def(
      "distance_to_aperture",
      &Bmad::distance_to_aperture,
      py::arg("orbit"),
      py::arg("particle_at"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine distance_to_aperture

Parameters
----------
orbit : CoordStruct
    Particle position.

particle_at : int
    first_track_edge$, second_track_edge$, or in_between$

ele : EleStruct
    Element containing aperture.

Returns
-------
no_aperture_here : bool
    True if aperture does not exist at the longitudinal location of the particle.

dist : float
    Normalized distance of the particle from the aperture.
)"""
  );
  m.def(
      "do_mode_flip",
      &Bmad::do_mode_flip,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine do_mode_flip

Parameters
----------
ele : EleStruct
    Starting Element
    This parameter is an input/output and is modified in-place.
    As an output, ele: Flipped element

Returns
-------
err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "dpc_given_de",
      &Bmad::dpc_given_de,
      py::arg("pc_old"),
      py::arg("mass"),
      py::arg("dE"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine dpc_given_de

Parameters
----------
pc_old : float

mass : float

dE : float

Returns
-------
dpc : float
)"""
  );
  m.def(
      "drift_and_pipe_track_methods_adjustment",
      &Bmad::drift_and_pipe_track_methods_adjustment,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
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
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with tracking methods adjusted if needed.
)"""
  );
  m.def(
      "drift_multipass_name_correction",
      &Bmad::drift_multipass_name_correction,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine drift_multipass_name_correction

Parameters
----------
lat : LatStruct
)"""
  );
  m.def(
      "drift_orbit_time",
      &Bmad::drift_orbit_time,
      py::arg("orbit"),
      py::arg("beta0"),
      py::arg("delta_s") = py::none(),
      py::arg("delta_t") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
)"""
  );
  m.def(
      "drift_particle_to_s",
      &Bmad::drift_particle_to_s,
      py::arg("p"),
      py::arg("s"),
      py::arg("branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine drift_particle_to_s (p, s, branch)

Drift a particle to a given s-coordinate

Parameters
----------
p : CoordStruct
    Init particle position.
    This parameter is an input/output and is modified in-place.
    As an output, p: Final particle position.

s : float
    Target s coordinate.

branch : BranchStruct
    Branch being tracked through.
)"""
  );
  m.def(
      "drift_particle_to_t",
      &Bmad::drift_particle_to_t,
      py::arg("p"),
      py::arg("t"),
      py::arg("branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine drift_particle_to_t (p, t, branch)

Drift a particle to a given t-coordinate

Parameters
----------
p : CoordStruct
    Init particle position.
    This parameter is an input/output and is modified in-place.
    As an output, p: Final particle position.

t : float
    Target t coordinate.

branch : BranchStruct
    Lattice branch being tracked through.
)"""
  );
  m.def(
      "dspline_len",
      &Bmad::dspline_len,
      py::arg("s_chord0"),
      py::arg("s_chord1"),
      py::arg("spline"),
      py::arg("dtheta_ref") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
)"""
  );
  m.def(
      "dynamic_aperture_point",
      &Bmad::dynamic_aperture_point,
      py::arg("branch"),
      py::arg("ele0"),
      py::arg("orb0"),
      py::arg("theta_xy"),
      py::arg("ap_param"),
      py::arg("check_xy_init") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
)"""
  );
  m.def(
      "dynamic_aperture_scan",
      &Bmad::dynamic_aperture_scan,
      py::arg("aperture_param"),
      py::arg("pz_start"),
      py::arg("lat"),
      py::arg("print_timing") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine dynamic_aperture_scan(aperture_scan, aperture_param, pz_start, lat, print_timing)

Routine to do a set of dynamic aperture scans.

Parameters
----------
aperture_param : ApertureParamStruct
    Scan parameters

pz_start : 1D array of float
    Starting phase space pz values.

lat : LatStruct
    Lattice.

print_timing : bool, optional
    If True print info on calculation time. Default is True.

Returns
-------
aperture_scan : 1D array of ApertureScanStruct
    Set of scans. One for each pz_start(:).
)"""
  );
}
