#include "pybmad/generated/Bmad_routines_t.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyTargetMinMaxCalc python_target_min_max_calc(
    FixedArray1D<Real, 3> r_corner1,
    FixedArray1D<Real, 3> r_corner2,
    double y_min,
    double y_max,
    double phi_min,
    double phi_max,
    std::optional<bool> initial = std::nullopt
) {
  Bmad::target_min_max_calc(r_corner1, r_corner2, y_min, y_max, phi_min, phi_max, initial);
  auto py_result{PyTargetMinMaxCalc{y_min, y_max, phi_min, phi_max}};
  return py_result;
}
PyTrack1TimeRungeKutta python_track1_time_runge_kutta(
    CoordStruct &orbit,
    EleStruct &ele,
    LatParamStruct &param,
    std::optional<double> t_end = std::nullopt,
    std::optional<double> dt_step = std::nullopt
) {
  auto _result = Bmad::track1_time_runge_kutta(orbit, ele, param, t_end, make_opt_ref(dt_step));
  auto py_result{PyTrack1TimeRungeKutta{_result, dt_step}};
  return py_result;
}
PyTrackADrift python_track_a_drift(
    CoordStruct &orb,
    double length,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<int> ele_orientation = std::nullopt,
    std::optional<bool> include_ref_motion = std::nullopt,
    std::optional<double> time = std::nullopt
) {
  Bmad::track_a_drift(
      orb,
      length,
      mat6,
      make_matrix,
      ele_orientation,
      include_ref_motion,
      make_opt_ref(time)
  );
  auto py_result{PyTrackADrift{time}};
  return py_result;
}
PyTypeTaylors python_type_taylors(
    TaylorStructArray1D bmad_taylor,
    std::optional<int> max_order = std::nullopt,
    optional_ref<CharacterAlloc1D> lines = std::nullopt,
    std::optional<int> n_lines = std::nullopt,
    std::optional<int> file_id = std::nullopt,
    std::optional<std::string> out_style = std::nullopt,
    std::optional<bool> clean = std::nullopt,
    std::optional<std::string> out_var_suffix = std::nullopt,
    std::optional<bool> append = std::nullopt
) {
  Bmad::type_taylors(
      bmad_taylor,
      max_order,
      lines,
      make_opt_ref(n_lines),
      file_id,
      out_style,
      clean,
      out_var_suffix,
      append
  );
  auto py_result{PyTypeTaylors{n_lines}};
  return py_result;
}

void init_Bmad_routines_t(py::module &m) {
  py::class_<Bmad::T6ToB123, std::unique_ptr<Bmad::T6ToB123>>(
      m,
      "T6ToB123",
      "t6_to_b123 return type"
  )
      .def_readonly("B1", &Bmad::T6ToB123::B1)
      .def_readonly("B2", &Bmad::T6ToB123::B2)
      .def_readonly("B3", &Bmad::T6ToB123::B3)
      .def_readonly("err_flag", &Bmad::T6ToB123::err_flag)
      .def("__len__", [](const Bmad::T6ToB123 &) { return 4; })
      .def("__getitem__", [](const Bmad::T6ToB123 &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.B1);
        if (i == 1)
          return py::cast(s.B2);
        if (i == 2)
          return py::cast(s.B3);
        if (i == 3)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "t6_to_b123",
      &Bmad::t6_to_b123,
      py::arg("t6"),
      py::arg("abz_tunes"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine t6_to_B123(N, abz_tunes, B1, B2, B3, err_flag)

This decomposes the one-turn matrix according to Equation 56 from
"Alternative approach to general coupled linear optics" by A. Wolski. PRSTAB.

Note that a sigma matrix can be assembeled from:  sigma = B1*emit_a + B2*emit_b + B3*emit_c

Parameters
----------
t6 : 2D array of float (shape: 6,6)
    1-turn transfer matrix.  RF assumed to be on.

abz_tunes : 1D array of float (shape: 3)
    a-mode and b-mode tunes.  Used to order eigensystem.

Returns
-------
B1 : 2D array of float (shape: 6,6)
    Beta matrix associated with a-mode.

B2 : 2D array of float (shape: 6,6)
    Beta matrix associated with b-mode.

B3 : 2D array of float (shape: 6,6)
    Beta matrix associated with c-mode.

err_flag : bool
    Set True if there is an error. False otherwise
)"""
  );
  m.def(
      "taper_mag_strengths",
      &Bmad::taper_mag_strengths,
      py::arg("lat"),
      py::arg("ref_lat") = py::none(),
      py::arg("except_") = py::none(),
      py::arg("err_flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine taper_mag_strengths

Parameters
----------
lat : LatStruct
    Lattice to vary.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with magnet strengths varied.

ref_lat : LatStruct, optional
    Reference lattice. If not present, lat will be used as the ref.

except : str, optional
    List of elements not to vary.

err_flag : bool, optional
)"""
  );
  py::class_<PyTargetMinMaxCalc, std::unique_ptr<PyTargetMinMaxCalc>>(
      m,
      "TargetMinMaxCalc",
      "target_min_max_calc return type"
  )
      .def_readonly("y_min", &PyTargetMinMaxCalc::y_min)
      .def_readonly("y_max", &PyTargetMinMaxCalc::y_max)
      .def_readonly("phi_min", &PyTargetMinMaxCalc::phi_min)
      .def_readonly("phi_max", &PyTargetMinMaxCalc::phi_max)
      .def("__len__", [](const PyTargetMinMaxCalc &) { return 4; })
      .def("__getitem__", [](const PyTargetMinMaxCalc &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.y_min);
        if (i == 1)
          return py::cast(s.y_max);
        if (i == 2)
          return py::cast(s.phi_min);
        if (i == 3)
          return py::cast(s.phi_max);
        throw py::index_error();
      });
  m.def(
      "target_min_max_calc",
      &python_target_min_max_calc,
      py::arg("r_corner1"),
      py::arg("r_corner2"),
      py::arg("y_min"),
      py::arg("y_max"),
      py::arg("phi_min"),
      py::arg("phi_max"),
      py::arg("initial") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine target_min_max_calc (r_corner1, r_corner2, y_min, y_max, phi_min, phi_max, initial)

Routine to calculate the min/max values for (y, phi).
min/max values are cumulative.

Parameters
----------
r_corner1 : 1D array of float (shape: 3)
    In target coords: A corner of the target. Must be normalized to 1.

r_corner2 : 1D array of float (shape: 3)
    In target coords: Adjacent corner of the target. Must be normalized to 1.

y_min : float
    min/max values. Only needed if initial = False.

y_max : float
    min/max values. Only needed if initial = False.

phi_min : float
    min/max values. Only needed if initial = False.

phi_max : float
    min/max values. Only needed if initial = False.

initial : bool, optional
    If present and True then this is the first edge for computation.

Returns
-------
y_min : float
    min/max values. Only needed if initial = False.

y_max : float
    min/max values. Only needed if initial = False.

phi_min : float
    min/max values. Only needed if initial = False.

phi_max : float
    min/max values. Only needed if initial = False.
)"""
  );
  py::class_<Bmad::TargetRotMats, std::unique_ptr<Bmad::TargetRotMats>>(
      m,
      "TargetRotMats",
      "target_rot_mats return type"
  )
      .def_readonly("w_to_target", &Bmad::TargetRotMats::w_to_target)
      .def_readonly("w_to_ele", &Bmad::TargetRotMats::w_to_ele)
      .def("__len__", [](const Bmad::TargetRotMats &) { return 2; })
      .def("__getitem__", [](const Bmad::TargetRotMats &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.w_to_target);
        if (i == 1)
          return py::cast(s.w_to_ele);
        throw py::index_error();
      });
  m.def(
      "target_rot_mats",
      &Bmad::target_rot_mats,
      py::arg("r_center"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine target_rot_mats (r_center, w_to_target, w_to_ele)

Routine to calculate the rotation matrices between ele coords and "target" coords.
By definition, in target coords r_center = [0, 0, 1].

Parameters
----------
r_center : 1D array of float (shape: 3)
    In lab coords: Center of target relative to phton emission point.

Returns
-------
w_to_target : 2D array of float (shape: 3,3)
    Rotation matrix from ele to target coords.

w_to_ele : 2D array of float (shape: 3,3)
    Rotation matrix from target to ele coords.
)"""
  );
  m.def(
      "taylor_equal_taylor",
      &Bmad::taylor_equal_taylor,
      py::arg("taylor1"),
      py::arg("taylor2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine taylor_equal_taylor

Parameters
----------
taylor1 : TaylorStruct

taylor2 : TaylorStruct
)"""
  );
  m.def(
      "taylor_inverse",
      &Bmad::taylor_inverse,
      py::arg("taylor_in"),
      py::arg("taylor_inv"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine taylor_inverse (taylor_in, taylor_inv, err)

Subroutine to invert a taylor map. Since the inverse map is truncated, it is not exact.

Parameters
----------
taylor_in : 1D array of TaylorStruct
    Input taylor map.

taylor_inv : 1D array of TaylorStruct
    Inverted taylor map.

Returns
-------
err : bool, optional
    Set True if there is no inverse. If not present then print an error message.
)"""
  );
  m.def(
      "taylor_propagate1",
      &Bmad::taylor_propagate1,
      py::arg("orb_taylor"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("ref_in") = py::none(),
      py::arg("spin_taylor") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine taylor_propagate1 (orb_taylor, ele, param, err_flag, ref_in, spin_taylor)

Subroutine to track (symplectic integration) a orbital map, and optionally a spin map, through an element.
The spin tracking is only done if spin_taylor is present and bmad_com%spin_tracking_on = T.
The alternative routine, if ele has a taylor map, is concat_taylor.

This routine will fail if there is no corresponding ptc fibre for this
element. In general, the transfer_map_calc routine should be used instead.

Parameters
----------
orb_taylor : 1D array of TaylorStruct
    Map to be tracked
    This parameter is an input/output and is modified in-place.
    As an output, orb_taylor: Map through element.

ele : EleStruct
    Element to track through

param : LatParamStruct

ref_in : CoordStruct, optional
    Particle to be tracked. Must be present if the particle to be tracked is not the reference particle or if
    the direction of propagation is backwards.

spin_taylor : 1D array of TaylorStruct, optional
    Spin map to be tracked
    This parameter is an input/output and is modified in-place.
    As an output, spin_taylor: Tracked spin map.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "taylor_to_mad_map",
      &Bmad::taylor_to_mad_map,
      py::arg("taylor"),
      py::arg("energy"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine taylor_to_mad_map (taylor, energy, map)

Subroutine to convert a Taylor map to a mad order 2 map.
If any of the Taylor terms have order greater than 2 they are ignored.

Parameters
----------
taylor : 1D array of TaylorStruct
    Taylor map.

energy : MadEnergyStruct
    Energy numbers.

Returns
-------
map : MadMapStruct
    Order 2 map.
)"""
  );
  m.def(
      "taylors_equal_taylors",
      &Bmad::taylors_equal_taylors,
      py::arg("taylor1"),
      py::arg("taylor2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine taylors_equal_taylors

Parameters
----------
taylor1 : 1D array of TaylorStruct

taylor2 : 1D array of TaylorStruct
)"""
  );
  m.def(
      "tilt_coords",
      &Bmad::tilt_coords,
      py::arg("tilt_val"),
      py::arg("coord"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tilt_coords

Parameters
----------
tilt_val : float
    Tilt value (could be the roll value for a bend)

coord : 1D array of float
    Coordinates of particle before rotation.
    This parameter is an input/output and is modified in-place.
    As an output, coord: Coordinates of particle after rotation.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before tilt.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix transfer matrix after tilt applied.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "tilt_coords_photon",
      &Bmad::tilt_coords_photon,
      py::arg("tilt_val"),
      py::arg("coord"),
      py::arg("w_mat") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tilt_coords_photon

Parameters
----------
tilt_val : float
    Tilt value (could be the roll value for a bend)

coord : 1D array of float
    Coordinates of particle before rotation.
    This parameter is an input/output and is modified in-place.
    As an output, coord: Coordinates of particle after rotation.

w_mat : 2D array of float (shape: 3,3), optional
    Rotation matrix before tilt.
    This parameter is an input/output and is modified in-place.
    As an output, w_mat: Rotation matrix after tilt.
)"""
  );
  m.def(
      "tilt_mat6",
      &Bmad::tilt_mat6,
      py::arg("mat6"),
      py::arg("tilt"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tilt_mat6

Parameters
----------
mat6 : 2D array of float (shape: 6,6)
    Untilted matrix.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Tilted matrix.

tilt : float
    Tilt angle.
)"""
  );
  py::class_<Bmad::ToEtaReading, std::unique_ptr<Bmad::ToEtaReading>>(
      m,
      "ToEtaReading",
      "to_eta_reading return type"
  )
      .def_readonly("reading", &Bmad::ToEtaReading::reading)
      .def_readonly("err", &Bmad::ToEtaReading::err)
      .def("__len__", [](const Bmad::ToEtaReading &) { return 2; })
      .def("__getitem__", [](const Bmad::ToEtaReading &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.reading);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "to_eta_reading",
      &Bmad::to_eta_reading,
      py::arg("eta_actual"),
      py::arg("ele"),
      py::arg("axis"),
      py::arg("add_noise"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine to_eta_reading (eta, ele, axis, add_noise, reading, err)

Compute the measured dispersion reading given the true dispersion and the
monitor offsets, noise, etc.

This routine will only give a nonzero reading for Bmad markers,
monitors, and instruments.

Parameters
----------
eta_actual : 1D array of float
    Actual (eta_x, eta_y) dispersion.

ele : EleStruct
    Element where the orbit is measured.

axis : int
    x_plane$ or y_plane$

add_noise : bool
    If True add noise to the reading

Returns
-------
reading : float
    BPM reading

err : bool
    Set True if there is an error. False otherwise.
)"""
  );
  py::class_<Bmad::ToFieldmapCoords, std::unique_ptr<Bmad::ToFieldmapCoords>>(
      m,
      "ToFieldmapCoords",
      "to_fieldmap_coords return type"
  )
      .def_readonly("x", &Bmad::ToFieldmapCoords::x)
      .def_readonly("y", &Bmad::ToFieldmapCoords::y)
      .def_readonly("z", &Bmad::ToFieldmapCoords::z)
      .def_readonly("cos_ang", &Bmad::ToFieldmapCoords::cos_ang)
      .def_readonly("sin_ang", &Bmad::ToFieldmapCoords::sin_ang)
      .def_readonly("err_flag", &Bmad::ToFieldmapCoords::err_flag)
      .def("__len__", [](const Bmad::ToFieldmapCoords &) { return 6; })
      .def("__getitem__", [](const Bmad::ToFieldmapCoords &s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.x);
        if (i == 1)
          return py::cast(s.y);
        if (i == 2)
          return py::cast(s.z);
        if (i == 3)
          return py::cast(s.cos_ang);
        if (i == 4)
          return py::cast(s.sin_ang);
        if (i == 5)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "to_fieldmap_coords",
      &Bmad::to_fieldmap_coords,
      py::arg("ele"),
      py::arg("local_orb"),
      py::arg("s_body"),
      py::arg("ele_anchor_pt"),
      py::arg("r0"),
      py::arg("curved_ref_frame"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine to_fieldmap_coords (ele, local_orb, s_body, ele_anchor_pt, r0, curved_ref_frame,
                                                              x, y, z, cos_ang, sin_ang, err_flag)

Routine to return the (x,y,s) position relative to a field map.

Parameters
----------
ele : EleStruct
    Element being tracked through.

local_orb : CoordStruct
    Particle orbit. Must be in local element coordinates.

s_body : float
    Longitudinal position relative to the entrance end of the element.

ele_anchor_pt : int
    anchor point of the field map (anchor_beginning$, anchor_center$, or anchor_end$).

r0 : 1D array of float (shape: 3)
    origin point of the fieldmap.

curved_ref_frame : bool
    If the element is a bend: Does the field map follow the bend reference coords?

Returns
-------
x : float
    Coords relative to the field map.

y : float
    Coords relative to the field map.

z : float
    Coords relative to the field map.

cos_ang : float
    cos and sin of coordinate rotation angle.

sin_ang : float
    cos and sin of coordinate rotation angle.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  py::class_<Bmad::ToOrbitReading, std::unique_ptr<Bmad::ToOrbitReading>>(
      m,
      "ToOrbitReading",
      "to_orbit_reading return type"
  )
      .def_readonly("reading", &Bmad::ToOrbitReading::reading)
      .def_readonly("err", &Bmad::ToOrbitReading::err)
      .def("__len__", [](const Bmad::ToOrbitReading &) { return 2; })
      .def("__getitem__", [](const Bmad::ToOrbitReading &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.reading);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "to_orbit_reading",
      &Bmad::to_orbit_reading,
      py::arg("orb"),
      py::arg("ele"),
      py::arg("axis"),
      py::arg("add_noise"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine to_orbit_reading (orb, ele, axis, add_noise, reading, err)

Calculate the measured reading on a bpm given the actual orbit and the
BPM's offsets, noise, etc.

This routine will only give a nonzero reading for Bmad markers,
monitors, and instruments.

Parameters
----------
orb : CoordStruct
    Orbit position at BPM.

ele : EleStruct
    Element where the orbit is measured.

axis : int
    x_plane$ or y_plane$

add_noise : bool
    If True add noise to the reading

Returns
-------
reading : float
    BPM reading

err : bool
    Set True if there is no valid reading. For example, if ele.is_on = False.
)"""
  );
  py::class_<Bmad::ToPhaseAndCouplingReading, std::unique_ptr<Bmad::ToPhaseAndCouplingReading>>(
      m,
      "ToPhaseAndCouplingReading",
      "to_phase_and_coupling_reading return type"
  )
      .def_readonly("reading", &Bmad::ToPhaseAndCouplingReading::reading)
      .def_readonly("err", &Bmad::ToPhaseAndCouplingReading::err)
      .def("__len__", [](const Bmad::ToPhaseAndCouplingReading &) { return 2; })
      .def("__getitem__", [](const Bmad::ToPhaseAndCouplingReading &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.reading);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "to_phase_and_coupling_reading",
      &Bmad::to_phase_and_coupling_reading,
      py::arg("ele"),
      py::arg("add_noise"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine to_phase_and_coupling_reading (ele, add_noise, reading, err)

Find the measured coupling values given the actual ones

This routine will only give a nonzero reading for Bmad markers,
monitors, and instruments.

Parameters
----------
ele : EleStruct
    Element where phase is measured.

add_noise : bool
    If True add noise to the reading

Returns
-------
reading : BpmPhaseCouplingStruct
    K and Cbar coupling parameters

err : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "to_photon_angle_coords",
      &Bmad::to_photon_angle_coords,
      py::arg("orb_in"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function to_photon_angle_coords (orb_in, ele) result (orb_out)

Routine to convert from standard photon coords to "angle" coords defined as:
      x, angle_x, y, angle_y, z, E-E_ref

Parameters
----------
orb_in : CoordStruct
    orbit in standard photon coords.

ele : EleStruct
    Reference element (generally the detector element.)

Returns
-------
orb_out : CoordStruct
    Transformed coordinates.
)"""
  );
  m.def(
      "to_surface_coords",
      &Bmad::to_surface_coords,
      py::arg("lab_orbit"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine to_surface_coords

Parameters
----------
lab_orbit : CoordStruct
    Photon position in laboratory coords.

ele : EleStruct
    Detector element.

Returns
-------
surface_orbit : CoordStruct
    Photon position in element body coordinates.
)"""
  );
  m.def(
      "touschek_lifetime",
      &Bmad::touschek_lifetime,
      py::arg("mode"),
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine touschek_lifetime(mode, Tl, lat)

Calculates the touschek lifetime for a lattice by calling touschek_rate1
for each element.
The loss rate at each element is averaged over one turn to obtain the lifetime.

This function assumes that the twiss parameters and closed orbit have
been calculated, and that mode has been populated.

This subroutine assumes a fixed momentum aperture.  The loss rate at each element
uses the same momentum aperture, mode%pz_aperture.

A common way to call this function is to first populate mode using
radiation integrals.  If an ideal lattice is used, the vertical
emittance must also be set to a reasonable value.  If the vertical
emittance is due only to quantum excitation, then it will likely be
several orders of magnitude smaller than any real physical situation, in which
case the integral in this function will have problems converging.

In addition to setting mode, also set lat%param%n_part to the number of particles
per bunch.

Parameters
----------
mode : NormalModesStruct
    beam properties

lat : LatStruct
    Accelerator Lattice

Returns
-------
Tl : float
    Touschek lifetime in seconds
)"""
  );
  m.def(
      "touschek_rate1",
      &Bmad::touschek_rate1,
      py::arg("mode"),
      py::arg("lat"),
      py::arg("ix") = py::none(),
      py::arg("s") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine touschek_rate1(mode, rate, lat, ix, s)

Calculates the touschek rate at the location specified by s or ix
This calculation is based on Piwinski 1998 "The Touschek Effect In
Strong Focusing Storage Rings".  This is the most general case, equation
31.

This function uses twiss_and_track_at_s to determine the Twiss parameters
at the location s or element index ix.

A common way to call this function is to first populate mode using
radiation integrals.  If an ideal lattice is used, the vertical
emittance must also be set to a reasonable value.  If the vertical
emittance is due only to quantum excitation, then it will likely be
several orders of magnitude smaller than any real physical situation, in which
case the integral in this function will have problems converging.
Additionally, mode%pz_aperture needs to be set to the momentum aperture.

In addition to setting mode, also set lat%param%n_part to the number of particles
per bunch.

IMPORTANT NOTE: If the lattice type is a circular lattice, then
                mode%a%emittance and mode%b%emittance are assumed to
                contain the normalized emittences.  If lattice geometry is
                open, the emittances are assumed to be
                unnormalized.

IMPORTANT NOTE: The output of this subroutine is the loss rate assuming
                that two particles are lost per collision, one with too
                much energy, and one with too little energy.  This agrees
                with Piwinski's original derivation, which assumes that the
                positive energy aperture is equal in magnitude to the
                negative energy aperture.  If you are studying an
                accelerator with a non-symmetric energy aperture, then
                this subroutine should be called twice, once with the positive
                aperture, and once with the negative aperture, and rate from
                each call should be halved and summed.

Parameters
----------
mode : NormalModesStruct
    beam properties

lat : LatStruct
    Lattice

ix : int, optional
    element index (either s or ix must be specified)

s : float, optional
    location in meters (either s or ix must be specified)

Returns
-------
rate : float
    Touschek rate, in units particle per second, assuming two particles per event.
)"""
  );
  m.def(
      "touschek_rate1_zap",
      &Bmad::touschek_rate1_zap,
      py::arg("mode"),
      py::arg("rate"),
      py::arg("lat"),
      py::arg("ix") = py::none(),
      py::arg("s") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine touschek_rate1_zap

Parameters
----------
mode : NormalModesStruct

rate : float

lat : LatStruct

ix : int, optional

s : float, optional
)"""
  );
  py::class_<Bmad::Track1, std::unique_ptr<Bmad::Track1>>(m, "Track1", "track1 return type")
      .def_readonly("end_orb", &Bmad::Track1::end_orb)
      .def_readonly("err_flag", &Bmad::Track1::err_flag)
      .def("__len__", [](const Bmad::Track1 &) { return 2; })
      .def("__getitem__", [](const Bmad::Track1 &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.end_orb);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "track1",
      &Bmad::track1,
      py::arg("start_orb"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("track") = py::none(),
      py::arg("ignore_radiation") = py::none(),
      py::arg("make_map1") = py::none(),
      py::arg("init_to_edge") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1

Parameters
----------
start_orb : CoordStruct
    Starting position.

ele : EleStruct
    Element to track through.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Modified if make_map1 is True.

param : LatParamStruct
    Reference particle info.

track : TrackStruct, optional
    Structure holding existing track.
    This parameter is an input/output and is modified in-place.
    As an output, track: Structure holding the track information if the

ignore_radiation : bool, optional
    If present and True then do not include radiation effects along with space charge effects.

make_map1 : bool, optional
    Make ele.mat6 and ele.spin_q components? Default is false.

init_to_edge : bool, optional
    Default is True. If True then force the tracked particle to begin at the element's edge. See above. Do not
    use this argument unless you know what you are doing.

Returns
-------
end_orb : CoordStruct
    End position.

err_flag : bool, optional
    Set true if there is an error. False otherwise. Note: The particle getting lost (EG hitting an aperture)
    is *not* an error. An error is something like start_orb not being properly initialized.
)"""
  );
  m.def(
      "track1_beam",
      &Bmad::track1_beam,
      py::arg("beam"),
      py::arg("ele"),
      py::arg("centroid") = py::none(),
      py::arg("direction") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_beam (beam, ele, err, centroid, direction)

Subroutine to track a beam of particles through an element.

Parameters
----------
beam : BeamStruct
    Starting beam position.
    This parameter is an input/output and is modified in-place.
    As an output, beam: Ending beam position.

ele : EleStruct
    element to track through.

centroid : 1D array of CoordStruct, optional
    Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before beam tracking by
    tracking a single particle.

direction : int, optional
    +1 (default) -> Track forward, -1 -> Track backwards.

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost for a CSR calc.
)"""
  );
  py::class_<Bmad::Track1Bmad, std::unique_ptr<Bmad::Track1Bmad>>(
      m,
      "Track1Bmad",
      "track1_bmad return type"
  )
      .def_readonly("err_flag", &Bmad::Track1Bmad::err_flag)
      .def_readonly("track", &Bmad::Track1Bmad::track)
      .def("__len__", [](const Bmad::Track1Bmad &) { return 2; })
      .def("__getitem__", [](const Bmad::Track1Bmad &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.track);
        throw py::index_error();
      });
  m.def(
      "track1_bmad",
      &Bmad::track1_bmad,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_bmad

Parameters
----------
orbit : CoordStruct
    Starting position
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Element

param : LatParamStruct

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix propagated through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
err_flag : bool, optional
    Set true if there is an error. False otherwise.

track : TrackStruct, optional
    Structure holding the track information if the lattice element does tracking step-by-step. See track1 for
    more details.
)"""
  );
  m.def(
      "track1_bmad_photon",
      &Bmad::track1_bmad_photon,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_bmad_photon

Parameters
----------
orbit : CoordStruct
    Starting position
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position

ele : EleStruct
    Element

param : LatParamStruct

Returns
-------
err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  m.def(
      "track1_bunch",
      &Bmad::track1_bunch,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("centroid") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("bunch_track") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_bunch (bunch, ele, err, centroid, direction, bunch_track)

Subroutine to track a bunch of particles through an element.

Parameters
----------
bunch : BunchStruct
    Starting bunch position.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position.

ele : EleStruct
    element to track through.

centroid : 1D array of CoordStruct, optional
    Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before beam tracking by
    tracking a single particle.

direction : int, optional
    +1 (default) -> Track forward, -1 -> Track backwards.

bunch_track : BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: Track information appended to track.

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost for a CSR calc.
)"""
  );
  m.def(
      "track1_bunch_csr",
      &Bmad::track1_bunch_csr,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("centroid"),
      py::arg("s_start") = py::none(),
      py::arg("s_end") = py::none(),
      py::arg("bunch_track") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_bunch_csr (bunch, ele, centroid, err, s_start, s_end, bunch_track)

Routine to track a bunch of particles through an element with csr radiation effects.

Parameters
----------
bunch : BunchStruct
    Starting bunch position.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position.

ele : EleStruct
    The element to track through. Must be part of a lattice.

centroid : 1D array of CoordStruct
    coord_struct, Approximate beam centroid orbit for the lattice branch. Calculate this before beam tracking
    by tracking a single particle.

s_start : float, optional
    Starting position relative to ele. Default = 0

s_end : float, optional
    Ending position. Default is ele length.

bunch_track : BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: track information if the tracking method does

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost.
)"""
  );
  m.def(
      "track1_bunch_csr3d",
      &Bmad::track1_bunch_csr3d,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("centroid"),
      py::arg("s_start") = py::none(),
      py::arg("s_end") = py::none(),
      py::arg("bunch_track") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_bunch_csr3d (bunch, ele, centroid, err, bunch_track)

EXPERIMENTAL. NOT CURRENTLY OPERATIONAL!

Routine to track a bunch of particles through an element using
steady-state 3D CSR.

Parameters
----------
bunch : BunchStruct
    Starting bunch position.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position.

ele : EleStruct
    The element to track through. Must be part of a lattice.

centroid : 1D array of CoordStruct
    coord_struct, Approximate beam centroid orbit for the lattice branch. Calculate this before beam tracking
    by tracking a single particle.

s_start : float, optional
    Starting position relative to ele. Default = 0

s_end : float, optional
    Ending position. Default is ele length.

bunch_track : BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: track information if the tracking method does

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost.

Notes
-----
The core routines are from the OpenCSR package developed at:
https://github.com/ChristopherMayes/OpenCSR
)"""
  );
  m.def(
      "track1_bunch_hom",
      &Bmad::track1_bunch_hom,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("direction") = py::none(),
      py::arg("bunch_track") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_bunch_hom (bunch, ele, direction, bunch_track)

Subroutine to track a bunch of particles through an element including wakefields.

Parameters
----------
bunch : BunchStruct
    Starting bunch position.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position.

ele : EleStruct
    The element to track through.

direction : int, optional
    +1 (default) -> Track forward, -1 -> Track backwards.

bunch_track : BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: Track information appended to track.
)"""
  );
  m.def(
      "track1_bunch_space_charge",
      &Bmad::track1_bunch_space_charge,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("track_to_same_s") = py::none(),
      py::arg("bunch_track") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_bunch_space_charge

Parameters
----------
bunch : BunchStruct
    Starting bunch position.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position.

ele : EleStruct
    Element to track through. Must be part of a lattice.

track_to_same_s : bool, optional
    Default is True. If True, drift particles to all have the same s-position.

bunch_track : BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: track information if the tracking method does

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost for a CSR calc.
)"""
  );
  m.def(
      "track1_crystal",
      &Bmad::track1_crystal,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_crystal (ele, param, orbit)

Routine to track diffraction from a crystal.

Parameters
----------
ele : EleStruct
    Element tracking through.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_diffraction_plate_or_mask",
      &Bmad::track1_diffraction_plate_or_mask,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_diffraction_plate_or_mask (ele, param, orbit)

Routine to track through diffraction plate and mask elements.

Parameters
----------
ele : EleStruct
    Diffraction plate or mask element.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_high_energy_space_charge",
      &Bmad::track1_high_energy_space_charge,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_high_energy_space_charge (ele, param, orbit)

Routine to apply the ultra-relative space charge kick to a particle at the end of an element.
The routine setup_high_energy_space_charge_calc must be called initially before any tracking is done.
This routine assumes a Gaussian bunch and is only valid with relativistic particles where the
effect of the space charge is small.

Parameters
----------
ele : EleStruct
    Element tracked through.

param : LatParamStruct

orbit : CoordStruct
    Starting position
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position
)"""
  );
  m.def(
      "track1_lens",
      &Bmad::track1_lens,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_lens (ele, param, orbit)

Routine to track through a lens.

Parameters
----------
ele : EleStruct
    Element tracking through.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_linear",
      &Bmad::track1_linear,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_linear

Parameters
----------
orbit : CoordStruct
    Starting position
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position

ele : EleStruct
    Element

param : LatParamStruct
)"""
  );
  m.def(
      "track1_lr_wake",
      &Bmad::track1_lr_wake,
      py::arg("bunch"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_lr_wake (bunch, ele)

Subroutine to put in the long-range wakes for particle tracking.

Parameters
----------
bunch : BunchStruct
    Bunch to track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Kicked bunch.

ele : EleStruct
    Element with wakes.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with updated wake amplitudes.
)"""
  );
  m.def(
      "track1_mad",
      &Bmad::track1_mad,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_mad (orbit, ele, param)

Subroutine to track through an element using a 2nd order transfer map.
Note: If map does not exist then one will be created.

Parameters
----------
orbit : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coords.

ele : EleStruct
    Element to track through.

param : LatParamStruct
    Lattice parameters.
)"""
  );
  m.def(
      "track1_mirror",
      &Bmad::track1_mirror,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_mirror (ele, param, orbit)

Routine to track reflection from a mirror.

Parameters
----------
ele : EleStruct
    Element tracking through.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_mosaic_crystal",
      &Bmad::track1_mosaic_crystal,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_mosaic_crystal (ele, param, orbit)

Routine to track diffraction from a crystal.

Parameters
----------
ele : EleStruct
    Element tracking through.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_multilayer_mirror",
      &Bmad::track1_multilayer_mirror,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_multilayer_mirror (ele, param, orbit)

Routine to track reflection from a multilayer_mirror.
Basic equations are from Kohn, "On the Theory of Reflectivity of an X-Ray Multilayer Mirror".

Parameters
----------
ele : EleStruct
    Element tracking through.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_radiation",
      &Bmad::track1_radiation,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("edge"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_radiation (orbit, ele, edge)

Subroutine to apply a kick to a particle to account for radiation dampling and/or fluctuations.

For tracking through a given element, this routine should be called initially when
the particle is at the entrance end and at the end when the particle is at the exit end, when
the orbit is with respect to laboratory (not element body) coordinates.
That is, each time this routine is called it applies half the radiation kick for the entire element.

Note: This routine is called by track1.

Parameters
----------
orbit : CoordStruct
    Particle position before radiation applied.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Particle position after radiation has been applied.

ele : EleStruct
    Element generating radiation.

edge : int
    Where the particle is: start_edge$ or end_edge$.
)"""
  );
  m.def(
      "track1_radiation_center",
      &Bmad::track1_radiation_center,
      py::arg("orbit"),
      py::arg("ele1"),
      py::arg("ele2"),
      py::arg("rad_damp") = py::none(),
      py::arg("rad_fluct") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_radiation_center (orbit, ele1, ele2, rad_damp, rad_fluct)

Used for elements that have been split in half: This routine applies a kick to a particle
to account for radiation dampling and/or fluctuations.

Also see: track1_radiation.

Parameters
----------
orbit : CoordStruct
    Particle at center of element before radiation applied.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Particle position after radiation has been applied.

ele1 : EleStruct
    First half of the split element.

ele2 : EleStruct
    Second half of the split element.

rad_damp : bool, optional
    If present, override setting of bmad_com.radiation_damping_on.

rad_fluct : bool, optional
    If present, override setting of bmad_com.radiation_fluctuations_on.
)"""
  );
  py::class_<Bmad::Track1RungeKutta, std::unique_ptr<Bmad::Track1RungeKutta>>(
      m,
      "Track1RungeKutta",
      "track1_runge_kutta return type"
  )
      .def_readonly("err_flag", &Bmad::Track1RungeKutta::err_flag)
      .def_readonly("track", &Bmad::Track1RungeKutta::track)
      .def("__len__", [](const Bmad::Track1RungeKutta &) { return 2; })
      .def("__getitem__", [](const Bmad::Track1RungeKutta &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.track);
        throw py::index_error();
      });
  m.def(
      "track1_runge_kutta",
      &Bmad::track1_runge_kutta,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_runge_kutta

Parameters
----------
orbit : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coords.

ele : EleStruct
    Ele_struct

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix propagated through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.

track : TrackStruct, optional
    Structure holding the track information.
)"""
  );
  m.def(
      "track1_sample",
      &Bmad::track1_sample,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_sample (ele, param, orbit)

Routine to track reflection from a sample element.

Parameters
----------
ele : EleStruct
    Element tracking through.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords to be transformed
    This parameter is an input/output and is modified in-place.
    As an output, orbit: final phase-space coords
)"""
  );
  m.def(
      "track1_spin",
      &Bmad::track1_spin,
      py::arg("start_orb"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("end_orb"),
      py::arg("make_quaternion") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_spin

Parameters
----------
start_orb : CoordStruct
    Starting coords.

ele : EleStruct
    Element to track through.

param : LatParamStruct
    Beam parameters.

end_orb : CoordStruct
    Ending coords.

make_quaternion : bool, optional
    If present and true then calculate the 1st order spin map which is represented as a quaternion.
)"""
  );
  m.def(
      "track1_spin_integration",
      &Bmad::track1_spin_integration,
      py::arg("start_orb"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("end_orb"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_spin_integration

Parameters
----------
start_orb : CoordStruct
    Starting coords.

ele : EleStruct
    Element to track through.

param : LatParamStruct
    Beam parameters.

end_orb : CoordStruct
    Ending coords.
)"""
  );
  m.def(
      "track1_spin_taylor",
      &Bmad::track1_spin_taylor,
      py::arg("start_orb"),
      py::arg("ele"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_spin_taylor

Parameters
----------
start_orb : CoordStruct
    Starting coords.

ele : EleStruct
    Element to track through.

param : LatParamStruct
    Beam parameters.

Returns
-------
end_orb : CoordStruct
)"""
  );
  m.def(
      "track1_sr_wake",
      &Bmad::track1_sr_wake,
      py::arg("bunch"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track1_sr_wake (bunch, ele)

Subroutine to apply the short range wake fields to a bunch.

Parameters
----------
bunch : BunchStruct
    Bunch of particles.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Bunch with wakefields applied to the particles.

ele : EleStruct
    Element with wakefields.
)"""
  );
  m.def(
      "track1_symp_lie_ptc",
      &Bmad::track1_symp_lie_ptc,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_symp_lie_ptc

Parameters
----------
orbit : CoordStruct
    Starting position
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position

ele : EleStruct
    Element

param : LatParamStruct

Returns
-------
track : TrackStruct, optional
    Structure holding the track information.
)"""
  );
  m.def(
      "track1_taylor",
      &Bmad::track1_taylor,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("taylor") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_taylor

Parameters
----------
orbit : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coords.

ele : EleStruct
    Element to track through.

taylor : 1D array of TaylorStruct (shape: 6), optional
    Alternative map to use instead of ele.taylor.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  py::class_<PyTrack1TimeRungeKutta, std::unique_ptr<PyTrack1TimeRungeKutta>>(
      m,
      "Track1TimeRungeKutta",
      "track1_time_runge_kutta return type"
  )
      .def_readonly("err_flag", &PyTrack1TimeRungeKutta::err_flag)
      .def_readonly("track", &PyTrack1TimeRungeKutta::track)
      .def_readonly("dt_step", &PyTrack1TimeRungeKutta::dt_step)
      .def("__len__", [](const PyTrack1TimeRungeKutta &) { return 3; })
      .def("__getitem__", [](const PyTrack1TimeRungeKutta &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.track);
        if (i == 2)
          return py::cast(s.dt_step);
        throw py::index_error();
      });
  m.def(
      "track1_time_runge_kutta",
      &python_track1_time_runge_kutta,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("t_end") = py::none(),
      py::arg("dt_step") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track1_time_runge_kutta

Parameters
----------
orbit : CoordStruct
    starting position, z-based coords
    This parameter is an input/output and is modified in-place.
    As an output, orbit: end position, z-based coords

ele : EleStruct
    element

param : LatParamStruct
    lattice parameters

t_end : float, optional
    If present, maximum time to which the particle will be tracked. Used for tracking with given time steps.
    The time orb.t at which tracking stops may be less than this if the particle gets to the end of the
    element

dt_step : float, optional
    If positive, next RK time step to take. This overrides bmad_com.init_ds_adaptive_tracking. Used by
    track_bunch_time.
    This parameter is an input/output and is modified in-place.
    As an output, dt_step: Next RK time step that this tracker would take based on the error tolerance.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise

track : TrackStruct, optional
    Contains array of the step-by-step particle trajectory along with the field at these positions. When
    tracking through multiple elements, the trajectory in an element is appended to the existing trajectory.
    To reset: Set track.n_pt = -1.

dt_step : float, optional
    If positive, next RK time step to take. This overrides bmad_com.init_ds_adaptive_tracking. Used by
    track_bunch_time.
    This parameter is an input/output and is modified in-place.
    As an output, dt_step: Next RK time step that this tracker would take based on the error tolerance.
)"""
  );
  py::class_<Bmad::TrackABeambeam, std::unique_ptr<Bmad::TrackABeambeam>>(
      m,
      "TrackABeambeam",
      "track_a_beambeam return type"
  )
      .def_readonly("track", &Bmad::TrackABeambeam::track)
      .def_readonly("mat6", &Bmad::TrackABeambeam::mat6)
      .def("__len__", [](const Bmad::TrackABeambeam &) { return 2; })
      .def("__getitem__", [](const Bmad::TrackABeambeam &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.track);
        if (i == 1)
          return py::cast(s.mat6);
        throw py::index_error();
      });
  m.def(
      "track_a_beambeam",
      &Bmad::track_a_beambeam,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_beambeam

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Beambeam element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
track : TrackStruct, optional
    Structure holding the track information if the lattice element does tracking step-by-step. See track1 for
    more details.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_bend",
      &Bmad::track_a_bend,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_bend

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Bend element.

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix to the element end.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "track_a_bend_photon",
      &Bmad::track_a_bend_photon,
      py::arg("orb"),
      py::arg("ele"),
      py::arg("length"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_a_bend_photon (orb, ele, length)

Routine to track a photon through a dipole bend.
The photon is traveling in a straight line but the reference frame
is curved in a circular shape.

Parameters
----------
orb : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orb: End position.

ele : EleStruct
    Bend element.

length : float
    length to track.
)"""
  );
  m.def(
      "track_a_capillary",
      &Bmad::track_a_capillary,
      py::arg("orb"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_a_capillary (orb, ele)

Routine to track through a capillary.

Parameters
----------
orb : CoordStruct
    Input photon coordinates.
    This parameter is an input/output and is modified in-place.
    As an output, orb: Output photon coordinates.

ele : EleStruct
    Capillary element
)"""
  );
  m.def(
      "track_a_converter",
      &Bmad::track_a_converter,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_converter

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    converter element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_crab_cavity",
      &Bmad::track_a_crab_cavity,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_crab_cavity

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    crab_cavity element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  py::class_<PyTrackADrift, std::unique_ptr<PyTrackADrift>>(
      m,
      "TrackADrift",
      "track_a_drift return type"
  )
      .def_readonly("time", &PyTrackADrift::time)
      .def("__len__", [](const PyTrackADrift &) { return 1; })
      .def("__getitem__", [](const PyTrackADrift &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.time);
        throw py::index_error();
      });
  m.def(
      "track_a_drift",
      &python_track_a_drift,
      py::arg("orb"),
      py::arg("length"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::arg("ele_orientation") = py::none(),
      py::arg("include_ref_motion") = py::none(),
      py::arg("time") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_drift

Parameters
----------
orb : CoordStruct
    Orbit at start of the drift.
    This parameter is an input/output and is modified in-place.
    As an output, orb: Orbit at end of the drift.

length : float
    Length to drift through in body coordinates. --    If orb.direction = 1, positive length is in +z
    direction and vice versa.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the drift.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix including the drift.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

ele_orientation : int, optional
    Element orientation. Default is orb.direction.

include_ref_motion : bool, optional
    Include effect of the motion of the reference particle? Default is True. False is basically only used by
    offset_particle. Additionally, if False, orb.s is not changed.

time : float, optional
    Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
    This parameter is an input/output and is modified in-place.
    As an output, time: Updated time.

Returns
-------
time : float, optional
    Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
    This parameter is an input/output and is modified in-place.
    As an output, time: Updated time.
)"""
  );
  m.def(
      "track_a_drift_photon",
      &Bmad::track_a_drift_photon,
      py::arg("orb"),
      py::arg("length"),
      py::arg("phase_relative_to_ref"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_drift_photon

Parameters
----------
orb : CoordStruct
    Orbit at start of the drift.
    This parameter is an input/output and is modified in-place.
    As an output, orb: Orbit at end of the drift

length : float
    Longitudinal length to drift through.

phase_relative_to_ref : bool
    If true then E field phase shift is relative to ref particle.
)"""
  );
  m.def(
      "track_a_foil",
      &Bmad::track_a_foil,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_foil

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    foil element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_gkicker",
      &Bmad::track_a_gkicker,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_gkicker

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Gkicker

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "track_a_lcavity",
      &Bmad::track_a_lcavity,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_lcavity

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Thick multipole element.

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "track_a_lcavity_old",
      &Bmad::track_a_lcavity_old,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_lcavity_old

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Thick multipole element.

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "track_a_mask",
      &Bmad::track_a_mask,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_mask

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Mask element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_match",
      &Bmad::track_a_match,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("err_flag") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_match

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Match element.

param : LatParamStruct
    Lattice parameters.

err_flag : bool, optional

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  py::class_<Bmad::TrackAPatch, std::unique_ptr<Bmad::TrackAPatch>>(
      m,
      "TrackAPatch",
      "track_a_patch return type"
  )
      .def_readonly("s_ent", &Bmad::TrackAPatch::s_ent)
      .def_readonly("ds_ref", &Bmad::TrackAPatch::ds_ref)
      .def_readonly("mat6", &Bmad::TrackAPatch::mat6)
      .def("__len__", [](const Bmad::TrackAPatch &) { return 3; })
      .def("__getitem__", [](const Bmad::TrackAPatch &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.s_ent);
        if (i == 1)
          return py::cast(s.ds_ref);
        if (i == 2)
          return py::cast(s.mat6);
        throw py::index_error();
      });
  m.def(
      "track_a_patch",
      &Bmad::track_a_patch,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("drift_to_exit") = py::none(),
      py::arg("track_spin") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_patch

Parameters
----------
ele : EleStruct
    patch element.

orbit : CoordStruct
    Starting phase space coords
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Coords after applying a patch transformation.

drift_to_exit : bool, optional
    If False then do not drift the particle from beginning to end face. Also do not correct for a reference
    energy shift. Default is True.

track_spin : bool, optional
    If True rotate the spin vector appropriately. If ele.spin_tracking_method = symp_lie_ptc -> default =
    True. Else default = False.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
s_ent : float, optional
    Longitudinal coordinate of the initial particle position in the frame of reference of the face where the
    particle exits. For a patch with positive z_offset and all other attributes zero, s_ent = -z_offset.

ds_ref : float, optional
    Distance reference particle travels from entrance to exit.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_patch_photon",
      &Bmad::track_a_patch_photon,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("drift_to_exit") = py::none(),
      py::arg("use_z_pos") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_a_patch_photon (ele, orbit, drift_to_exit, use_z_pos)

Routine to track through a patch element with a photon.
The steps for tracking are:
  1) Transform from entrance to exit coordinates.
  2) Drift particle from the entrance to the exit coordinants.

Parameters
----------
ele : EleStruct
    patch element.

orbit : CoordStruct
    Starting phase space coords
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Coords after applying a patch transformation.

drift_to_exit : bool, optional
    If False then do not drift the particle from start to ending faces. Default is True.

use_z_pos : bool, optional
    If present and True, use orbit.vec(5) as the true z-position relative to the start of the element instead
    of assuming that the particle is at the patch edge.
)"""
  );
  m.def(
      "track_a_pickup",
      &Bmad::track_a_pickup,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("err_flag") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_pickup

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Pickup element.

param : LatParamStruct
    Lattice parameters.

err_flag : bool, optional

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_quadrupole",
      &Bmad::track_a_quadrupole,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_quadrupole

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Quadrupole element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_rfcavity",
      &Bmad::track_a_rfcavity,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_rfcavity

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    rfcavity element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_sad_mult",
      &Bmad::track_a_sad_mult,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_sad_mult

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Sad_mult element.

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the sad_mult.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "track_a_sol_quad",
      &Bmad::track_a_sol_quad,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_sol_quad

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Sol_quad or solenoid element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  m.def(
      "track_a_thick_multipole",
      &Bmad::track_a_thick_multipole,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_thick_multipole

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Thick multipole element.

param : LatParamStruct
    Lattice parameters.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "track_a_wiggler",
      &Bmad::track_a_wiggler,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("make_matrix") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_wiggler

Parameters
----------
orbit : CoordStruct
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: End position.

ele : EleStruct
    Wiggler element.

param : LatParamStruct
    Lattice parameters.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix through the element.
)"""
  );
  py::class_<Bmad::TrackAZeroLengthElement, std::unique_ptr<Bmad::TrackAZeroLengthElement>>(
      m,
      "TrackAZeroLengthElement",
      "track_a_zero_length_element return type"
  )
      .def_readonly("err_flag", &Bmad::TrackAZeroLengthElement::err_flag)
      .def_readonly("track", &Bmad::TrackAZeroLengthElement::track)
      .def("__len__", [](const Bmad::TrackAZeroLengthElement &) { return 2; })
      .def("__getitem__", [](const Bmad::TrackAZeroLengthElement &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.track);
        throw py::index_error();
      });
  m.def(
      "track_a_zero_length_element",
      &Bmad::track_a_zero_length_element,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_a_zero_length_element

Parameters
----------
orbit : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coords.

ele : EleStruct
    Element tracked through.

param : LatParamStruct
    Lattice parameters.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.

track : TrackStruct, optional
    Structure holding the track information.
)"""
  );
  py::class_<Bmad::TrackAll, std::unique_ptr<Bmad::TrackAll>>(
      m,
      "TrackAll",
      "track_all return type"
  )
      .def_readonly("track_state", &Bmad::TrackAll::track_state)
      .def_readonly("err_flag", &Bmad::TrackAll::err_flag)
      .def_readonly("orbit0", &Bmad::TrackAll::orbit0)
      .def("__len__", [](const Bmad::TrackAll &) { return 3; })
      .def("__getitem__", [](const Bmad::TrackAll &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.track_state);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.orbit0);
        throw py::index_error();
      });
  m.def(
      "track_all",
      &Bmad::track_all,
      py::arg("lat"),
      py::arg("orbit"),
      py::arg("ix_branch") = py::none(),
      py::arg("init_lost") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_all

Parameters
----------
lat : LatStruct
    Lat to track through.

orbit : 1D array of CoordStruct
    orbit(0) is the starting coordinates for tracking. If not allocated, the zero orbit will be used.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Orbit array.

ix_branch : int, optional
    Index of branch to track. Default is 0 (main branch).

init_lost : bool, optional
    Default if False. If True, initialize orbit(N) terms that are not tracked through due to particle loss.

Returns
-------
track_state : int, optional
    Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.

err_flag : bool, optional
    Set true if particle lost or error. False otherwise

orbit0 : 1D array of CoordStruct, optional
    Orbit array for branch 0. Used to fill in the orbit at lord elemenets. Only needed when orbit(:) is not
    the orbit for branch 0.
)"""
  );
  m.def(
      "track_beam",
      &Bmad::track_beam,
      py::arg("lat"),
      py::arg("beam"),
      py::arg("ele1") = py::none(),
      py::arg("ele2") = py::none(),
      py::arg("centroid") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("bunch_tracks") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_beam (lat, beam, ele1, ele2, err, centroid, direction, bunch_tracks)

Subroutine to track a beam of particles from the end of
ele1 Through to the end of ele2. Both must be in the same lattice branch.

Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.

Parameters
----------
lat : LatStruct
    Lattice to track through.

beam : BeamStruct
    Beam at end of element ix1.
    This parameter is an input/output and is modified in-place.
    As an output, beam: Beam at end of element ix2.

ele1 : EleStruct, optional
    Starting element (this element is NOT tracked through). Default is lat.ele(0).

ele2 : EleStruct, optional
    Ending element. Default is lat.ele(lat.n_ele_track).

centroid : 1D array of CoordStruct, optional
    Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before beam tracking by
    tracking a single particle.

direction : int, optional
    +1 (default) -> Track forward, -1 -> Track backwards.

bunch_tracks : 1D array of BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_tracks: track information if the tracking method does

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost for a CSR calc.
)"""
  );
  m.def(
      "track_bunch",
      &Bmad::track_bunch,
      py::arg("lat"),
      py::arg("bunch"),
      py::arg("ele1") = py::none(),
      py::arg("ele2") = py::none(),
      py::arg("centroid") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("bunch_track") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_bunch (lat, bunch, ele1, ele2, err, centroid, direction, bunch_track)

Subroutine to track a particle bunch from the end of ele1 Through to the end of ele2.
Both must be in the same lattice branch.
With forward tracking, if ele2 is at or before ele1, the tracking will "wrap" around
the ends of the lattice.

Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.

Parameters
----------
lat : LatStruct
    Lattice to track through.

bunch : BunchStruct
    Bunch at end of element ix1.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Bunch at end of element ix2.

ele1 : EleStruct, optional
    Starting element (this element is NOT tracked through). Default is lat.ele(0).

ele2 : EleStruct, optional
    Ending element. Default is lat.ele(lat.n_ele_track).

centroid : 1D array of CoordStruct, optional
    Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before bunch tracking by
    tracking a single particle.

direction : int, optional
    +1 (default) -> Track forward, -1 -> Track backwards.

bunch_track : BunchTrackStruct, optional
    Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: track information if the tracking method does

Returns
-------
err : bool
    Set true if there is an error. EG: Too many particles lost for a CSR calc.
)"""
  );
  m.def(
      "track_bunch_time",
      &Bmad::track_bunch_time,
      py::arg("bunch"),
      py::arg("branch"),
      py::arg("t_end"),
      py::arg("s_end"),
      py::arg("dt_step") = py::none(),
      py::arg("extra_field") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_bunch_time

Parameters
----------
bunch : BunchStruct
    Coordinates must be time-coords in element body frame.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Coordinates will be time-coords in element body frame.

branch : BranchStruct
    Lattice branch being tracked through.

t_end : float
    Ending time.

s_end : float
    Ending s-position.

dt_step : 1D array of float, optional
    Initial step to take for each particle. Overrides bmad_com.init_ds_adaptive_tracking.
    This parameter is an input/output and is modified in-place.
    As an output, dt_step: Next RK time step that this tracker would take based on the error tolerance.

extra_field : 1D array of EmFieldStruct, optional
    Per particle static field to be added to the lattice element field. Eg used with space charge.
)"""
  );
  m.def(
      "track_bunch_to_s",
      &Bmad::track_bunch_to_s,
      py::arg("bunch"),
      py::arg("s"),
      py::arg("branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_bunch_to_s (bunch, s, branch)

Drift a bunch of particles to the same s coordinate

Parameters
----------
bunch : BunchStruct
    Input bunch position in s-based coordinate.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Output bunch position in s-based coordinate. Particles will be at the same s
    coordinate

s : float
    Target coordinate.

branch : BranchStruct
    Branch being tracked through.
)"""
  );
  m.def(
      "track_bunch_to_t",
      &Bmad::track_bunch_to_t,
      py::arg("bunch"),
      py::arg("t_target"),
      py::arg("branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_bunch_to_t (bunch, t_target, branch)

Drift a bunch of particles to the same t coordinate

Parameters
----------
bunch : BunchStruct
    Input bunch position in s-based coordinate.
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Output bunch position in s-based coordinate. Particles will be at the same t
    coordinate

t_target : float
    Target t coordinate.

branch : BranchStruct
    Lattice branch being tracked through.
)"""
  );
  m.def(
      "track_complex_taylor",
      &Bmad::track_complex_taylor,
      py::arg("start_orb"),
      py::arg("complex_taylor"),
      py::arg("end_orb"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine track_complex_taylor (start_orb, complex_taylor, end_orb)

Subroutine to track using a complex_taylor map.

Parameters
----------
start_orb : 1D array of complex
    Starting coords.

complex_taylor : 1D array of ComplexTaylorStruct
    complex_taylor map.

end_orb : 1D array of complex
    Ending coords.
)"""
  );
  py::class_<Bmad::TrackFromSToS, std::unique_ptr<Bmad::TrackFromSToS>>(
      m,
      "TrackFromSToS",
      "track_from_s_to_s return type"
  )
      .def_readonly("orbit_end", &Bmad::TrackFromSToS::orbit_end)
      .def_readonly("all_orb", &Bmad::TrackFromSToS::all_orb)
      .def_readonly("track_state", &Bmad::TrackFromSToS::track_state)
      .def("__len__", [](const Bmad::TrackFromSToS &) { return 3; })
      .def("__getitem__", [](const Bmad::TrackFromSToS &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.orbit_end);
        if (i == 1)
          return py::cast(s.all_orb);
        if (i == 2)
          return py::cast(s.track_state);
        throw py::index_error();
      });
  m.def(
      "track_from_s_to_s",
      &Bmad::track_from_s_to_s,
      py::arg("lat"),
      py::arg("s_start"),
      py::arg("s_end"),
      py::arg("orbit_start"),
      py::arg("ix_branch") = py::none(),
      py::arg("ix_ele_end") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_from_s_to_s

Parameters
----------
lat : LatStruct
    Lattice to track through

s_start : float
    Starting s-position.

s_end : float
    Ending s-position. If <= s_start then will wrap

orbit_start : CoordStruct
    Starting coordinates.

ix_branch : int, optional
    Lattice branch index. Default is 0 (main branch).

ix_ele_end : int, optional
    If present, ignore s_end and track to in between ix_ele_end and ix_ele_end+1

Returns
-------
orbit_end : CoordStruct
    Ending coordinates.

all_orb : 1D array of CoordStruct, optional
    If present then the orbit at the exit ends of the elements tracked through will be recorded in this
    structure.

track_state : int, optional
    Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.
)"""
  );
  m.def(
      "track_many",
      &Bmad::track_many,
      py::arg("lat"),
      py::arg("orbit"),
      py::arg("ix_start"),
      py::arg("ix_end"),
      py::arg("direction"),
      py::arg("ix_branch") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_many

Parameters
----------
lat : LatStruct
    Lat to track through.

orbit : 1D array of CoordStruct
    Coordinates at start of tracking.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Orbit.

ix_start : int
    Start index (See Note).

ix_end : int
    End index (See Note).

direction : int
    Direction to track. = +1 -> Track forward (+s) = -1 -> Track backward (-s)

ix_branch : int, optional
    Branch to track. Default is 0 (main lattice).

Returns
-------
track_state : int, optional
    Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.
)"""
  );
  m.def(
      "track_to_surface",
      &Bmad::track_to_surface,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("param"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine track_to_surface

Parameters
----------
ele : EleStruct
    Element

orbit : CoordStruct
    Coordinates in the element coordinate frame
    This parameter is an input/output and is modified in-place.
    As an output, orbit: At surface in local surface coordinate frame

param : LatParamStruct
    Branch parameters.

Returns
-------
w_surface : 2D array of float (shape: 3,3)
    real(rp), rotation matrix to transform to surface coords.
)"""
  );
  py::class_<Bmad::TrackUntilDead, std::unique_ptr<Bmad::TrackUntilDead>>(
      m,
      "TrackUntilDead",
      "track_until_dead return type"
  )
      .def_readonly("end_orb", &Bmad::TrackUntilDead::end_orb)
      .def_readonly("track", &Bmad::TrackUntilDead::track)
      .def("__len__", [](const Bmad::TrackUntilDead &) { return 2; })
      .def("__getitem__", [](const Bmad::TrackUntilDead &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.end_orb);
        if (i == 1)
          return py::cast(s.track);
        throw py::index_error();
      });
  m.def(
      "track_until_dead",
      &Bmad::track_until_dead,
      py::arg("start_orb"),
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(subroutine track_until_dead (start_orb, lat, end_orb, track)

Subroutine to track a particle arbitrarily through a lattice, forwards or backwards,
  until it is lost or exits the lattice.

  The starting element is located using start_orb%s.

Parameters
----------
start_orb : CoordStruct
    Starting coords.

lat : LatStruct
    lattice that contains and element at start_orb.s

Returns
-------
end_orb : CoordStruct
    final coords

track : TrackStruct, optional
    (optional)
)"""
  );
  py::class_<Bmad::TrackingRadMapSetup, std::unique_ptr<Bmad::TrackingRadMapSetup>>(
      m,
      "TrackingRadMapSetup",
      "tracking_rad_map_setup return type"
  )
      .def_readonly("rad_map", &Bmad::TrackingRadMapSetup::rad_map)
      .def_readonly("err_flag", &Bmad::TrackingRadMapSetup::err_flag)
      .def("__len__", [](const Bmad::TrackingRadMapSetup &) { return 2; })
      .def("__getitem__", [](const Bmad::TrackingRadMapSetup &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.rad_map);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "tracking_rad_map_setup",
      &Bmad::tracking_rad_map_setup,
      py::arg("ele"),
      py::arg("tollerance"),
      py::arg("ref_edge"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tracking_rad_map_setup

Parameters
----------
ele : EleStruct
    Element to setup. Matrices will be with respect to the map reference orbit.

tollerance : float
    Tolerance used for the computation.

ref_edge : int
    Edge that the matrices are referenced to. upstream_end$ or downstream_end$.

Returns
-------
rad_map : RadMapStruct
    Structure holding the matrices.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "transfer_ac_kick",
      &Bmad::transfer_ac_kick,
      py::arg("ac_in"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_ac_kick

Parameters
----------
ac_in : AcKickerStruct
    Input

Returns
-------
ac_out : AcKickerStruct, optional
    Gets set equal to ac_in
)"""
  );
  m.def(
      "transfer_branch",
      &Bmad::transfer_branch,
      py::arg("branch1"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_branch

Parameters
----------
branch1 : BranchStruct

Returns
-------
branch2 : BranchStruct
)"""
  );
  m.def(
      "transfer_branch_parameters",
      &Bmad::transfer_branch_parameters,
      py::arg("branch_in"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_branch_parameters

Parameters
----------
branch_in : BranchStruct
    Input branch.

Returns
-------
branch_out : BranchStruct
    Output branch with parameters set.
)"""
  );
  m.def(
      "transfer_branches",
      &Bmad::transfer_branches,
      py::arg("branch1"),
      py::arg("branch2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_branches

Parameters
----------
branch1 : 1D array of BranchStruct

branch2 : 1D array of BranchStruct
)"""
  );
  m.def(
      "transfer_ele",
      &Bmad::transfer_ele,
      py::arg("ele1"),
      py::arg("nullify_pointers") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_ele

Parameters
----------
ele1 : EleStruct

nullify_pointers : bool, optional
    If present and True then nullify the pointers in ele2 except for the ele2.lat and ele2.lord pointers. This
    gives a "bare bones" copy where one does not have to worry about deallocating allocated structure
    components later.

Returns
-------
ele2 : EleStruct
)"""
  );
  m.def(
      "transfer_ele_taylor",
      &Bmad::transfer_ele_taylor,
      py::arg("ele_in"),
      py::arg("taylor_order") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_ele_taylor

Parameters
----------
ele_in : EleStruct
    Element with the Taylor map.

taylor_order : int, optional
    Order to truncate the Taylor map at.

Returns
-------
ele_out : EleStruct
    Element receiving the Taylor map truncated to order taylor_order.
)"""
  );
  m.def(
      "transfer_eles",
      &Bmad::transfer_eles,
      py::arg("ele1"),
      py::arg("ele2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_eles

Parameters
----------
ele1 : 1D array of EleStruct

ele2 : 1D array of EleStruct
)"""
  );
  m.def(
      "transfer_fieldmap",
      &Bmad::transfer_fieldmap,
      py::arg("ele_in"),
      py::arg("who"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_fieldmap

Parameters
----------
ele_in : EleStruct
    Input element.

who : int
    Possibilities are: all$, cartesian_map$, cylindrical_map$, or grid_field$

Returns
-------
ele_out : EleStruct
    Output element.
)"""
  );
  m.def(
      "transfer_fixer_params",
      &Bmad::transfer_fixer_params,
      py::arg("fixer"),
      py::arg("to_stored"),
      py::arg("orbit") = py::none(),
      py::arg("who") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function transfer_fixer_params(fixer, to_stored, orbit, who) result (is_ok)

Set parameters of fixer.

Parameters
----------
fixer : EleStruct
    Fixer element to set.

to_stored : bool
    If False, set real Twiss from stored. If True, set stored Twiss from real.

orbit : CoordStruct, optional
    Used for 'phase_space' transfers. Used for input if to_stored = True.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Used for 'phase_space' transfers. Used for output if to_stored = False.

who : str, optional
    Who to set. Possibilities are: Groups: 'all', ' ' (default and same as 'all') Note: This excludes all
    'start' sets., 'twiss', 'a_twiss', 'b_twiss', 'cmat', 'x_dispersion', 'y_dispersion', 'dispersion',
    'chromatic', 'orbit', 'phase_space', 'spin', 'x_plane', 'y_plane', 'z_plane', 'start', 'start_spin',
    'start_phase_space', Individula Parameters: 'x', 'px', 'cmat_11', etc.

Returns
-------
is_ok : bool
    logical
)"""
  );
  m.def(
      "transfer_lat",
      &Bmad::transfer_lat,
      py::arg("lat1"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_lat

Parameters
----------
lat1 : LatStruct

Returns
-------
lat2 : LatStruct
)"""
  );
  m.def(
      "transfer_lat_parameters",
      &Bmad::transfer_lat_parameters,
      py::arg("lat_in"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_lat_parameters

Parameters
----------
lat_in : LatStruct
    Input lat.

Returns
-------
lat_out : LatStruct
    Output lat with parameters set.
)"""
  );
  m.def(
      "transfer_map_calc",
      &Bmad::transfer_map_calc,
      py::arg("lat"),
      py::arg("orb_map"),
      py::arg("ix1") = py::none(),
      py::arg("ix2") = py::none(),
      py::arg("ref_orb") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("one_turn") = py::none(),
      py::arg("unit_start") = py::none(),
      py::arg("concat_if_possible") = py::none(),
      py::arg("spin_map") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_map_calc

Parameters
----------
lat : LatStruct
    Lattice used in the calculation.

orb_map : 1D array of TaylorStruct
    Initial map (used when unit_start = False)
    This parameter is an input/output and is modified in-place.
    As an output, orb_map: Transfer map.

ix1 : int, optional
    Element start index for the calculation. Default is 0.

ix2 : int, optional
    Element end index for the calculation. Default is lat.n_ele_track.

ref_orb : CoordStruct, optional
    Reference orbit/particle at s1 around which the map is made. This arg is needed if: unit_start = True or
    particle is not the same as the reference particle of the lattice.

ix_branch : int, optional
    Lattice branch index. Default is 0.

one_turn : bool, optional
    If present and True, and if ix1 = ix2, and the lattice is circular, then construct the one-turn map from
    ix1 back to ix1. Default = False.

unit_start : bool, optional
    If present and False then orb_map will be used as the starting map instead of the unit map. Default = True

concat_if_possible : bool, optional
    If present and True then use map concatenation rather than tracking if a map is present for a given
    lattice element. See above. Default is False.

spin_map : 1D array of TaylorStruct, optional
    Input quaternion spin map. Output only computed if bmad_com.spin_tracking_on = T
    This parameter is an input/output and is modified in-place.
    As an output, spin_map: Quaternion spin map.

Returns
-------
err_flag : bool
    Set True if problem like number overflow, etc.
)"""
  );
  py::class_<Bmad::TransferMapFromSToS, std::unique_ptr<Bmad::TransferMapFromSToS>>(
      m,
      "TransferMapFromSToS",
      "transfer_map_from_s_to_s return type"
  )
      .def_readonly("ref_orb_out", &Bmad::TransferMapFromSToS::ref_orb_out)
      .def_readonly("err_flag", &Bmad::TransferMapFromSToS::err_flag)
      .def("__len__", [](const Bmad::TransferMapFromSToS &) { return 2; })
      .def("__getitem__", [](const Bmad::TransferMapFromSToS &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ref_orb_out);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "transfer_map_from_s_to_s",
      &Bmad::transfer_map_from_s_to_s,
      py::arg("lat"),
      py::arg("t_map"),
      py::arg("s1") = py::none(),
      py::arg("s2") = py::none(),
      py::arg("ref_orb_in") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("one_turn") = py::none(),
      py::arg("unit_start") = py::none(),
      py::arg("concat_if_possible") = py::none(),
      py::arg("spin_map") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine transfer_map_from_s_to_s (lat, t_map, s1, s2, ref_orb_in, ref_orb_out, ix_branch,
                                         one_turn, unit_start, err_flag, concat_if_possible, spin_map)

Subroutine to calculate the transfer map between longitudinal positions s1 to s2.

If s2 < s1 and lat%param%geometry is closed$ then the
calculation will 'wrap around' the lattice end.
For example, if s1 = 900 and s2 = 10 then the t_map is the map from
element 900 to the lattice end plus from 0 through 10.

If s2 < s1 and lat%param%geometry is open$ then the inverse of the forward map of s2 -> s1 is computed.

If s2 = s1 then you get the unit map except if one_turn = True and the lattice is circular.

Parameters
----------
lat : LatStruct
    Lattice used in the calculation.

t_map : 1D array of TaylorStruct
    Initial map (used when unit_start = False)
    This parameter is an input/output and is modified in-place.
    As an output, t_map: Transfer map.

s1 : float, optional
    Element start position for the calculation. Default is 0.

s2 : float, optional
    Element end position for the calculation. Default is lat.param.total_length.

ref_orb_in : CoordStruct, optional
    Reference orbit/particle at s1 around which the map is made. This arg is needed if: unit_start = True or
    particle is not the same as the reference particle of the lattice.

ix_branch : int, optional
    Lattice branch index. Default is 0 (main branch).

one_turn : bool, optional
    If present and True, and s1 = s2, and the lattice is circular: Construct the one-turn map from s1 back to
    s1. Otherwise t_map is unchanged or the unit map if unit_start = T. Default = False.

unit_start : bool, optional
    If present and False then t_map will be used as the starting map instead of the unit map. Default = True

concat_if_possible : bool, optional
    If present and True then use map concatenation rather than tracking if a map is present for a given
    lattice element. See above. Default is False.

spin_map : 1D array of TaylorStruct, optional
    Initial spin map.
    This parameter is an input/output and is modified in-place.
    As an output, spin_map: Final spin map. Only computed if bmad_com.spin_tracking_on = T.

Returns
-------
ref_orb_out : CoordStruct, optional
    Ending coordinates of the reference orbit. This is also the actual orbit of particle

err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  m.def(
      "transfer_mat2_from_twiss",
      &Bmad::transfer_mat2_from_twiss,
      py::arg("twiss1"),
      py::arg("twiss2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_mat2_from_twiss

Parameters
----------
twiss1 : TwissStruct
    Twiss parameters at the initial point.

twiss2 : TwissStruct
    Twiss parameters at the end point.

Returns
-------
mat : 2D array of float (shape: 2,2)
    Transfer matrix between the two points.
)"""
  );
  m.def(
      "transfer_mat_from_twiss",
      &Bmad::transfer_mat_from_twiss,
      py::arg("ele1"),
      py::arg("ele2"),
      py::arg("orb1"),
      py::arg("orb2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_mat_from_twiss

Parameters
----------
ele1 : EleStruct
    Element with twiss parameters for the starting point.

ele2 : EleStruct
    Element with twiss parameters for the ending point.

orb1 : 1D array of float (shape: 6)
    Reference orbit at ele1 (affects m(i,6) dispersion terms).

orb2 : 1D array of float (shape: 6)
    Reference orbit at ele2 (affects m(i,6) dispersion terms).

Returns
-------
m : 2D array of float (shape: 6,6)
    Transfer matrix between the two points.
)"""
  );
  m.def(
      "transfer_matrix_calc",
      &Bmad::transfer_matrix_calc,
      py::arg("lat"),
      py::arg("xfer_mat"),
      py::arg("xfer_vec") = py::none(),
      py::arg("ix1") = py::none(),
      py::arg("ix2") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("one_turn") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_matrix_calc

Parameters
----------
lat : LatStruct
    Lattice used in the calculation.

xfer_mat : 2D array of float (shape: 6,6)

xfer_vec : 1D array of float (shape: 6), optional

ix1 : int, optional
    Element start index for the calculation. Default is 0.

ix2 : int, optional
    Element end index for the calculation. Defaults: If ix1 is not present: ix2 = lat.n_ele_track If ix1 is
    present and lattice is closed: Calculate the one-turn matrix from ix1 back to ix1.

ix_branch : int, optional
    Branch index. Default is 0.

one_turn : bool, optional
    If present and True, and ix1 = ix2, and the lattice is closed: Construct the one-turn matrix from ix1 back
    to ix1. If False, (the default), and ix1 = ix2, mat6 is the unit matrix.
)"""
  );
  m.def(
      "transfer_twiss",
      &Bmad::transfer_twiss,
      py::arg("ele_in"),
      py::arg("reverse") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_twiss

Parameters
----------
ele_in : EleStruct
    Element with existing Twiss parameters.

reverse : bool, optional
    Reverse alpha and coupling as if particle is going in the reversed direction? Default is False.

Returns
-------
ele_out : EleStruct
    Element receiving the Twiss parameters.
)"""
  );
  m.def(
      "transfer_wake",
      &Bmad::transfer_wake,
      py::arg("wake_in"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine transfer_wake

Parameters
----------
wake_in : WakeStruct
    Input wake.

Returns
-------
wake_out : WakeStruct, optional
    Output wake.
)"""
  );
  m.def(
      "truncate_complex_taylor_to_order",
      &Bmad::truncate_complex_taylor_to_order,
      py::arg("complex_taylor_in"),
      py::arg("order"),
      py::arg("complex_taylor_out"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine truncate_complex_taylor_to_order (complex_taylor_in, order, complex_taylor_out)

Subroutine to throw out all terms in a complex_taylor map that are above a certain order.

Parameters
----------
complex_taylor_in : 1D array of ComplexTaylorStruct
    Input complex_taylor map.

order : int
    Order above which terms are dropped.

complex_taylor_out : 1D array of ComplexTaylorStruct
    Truncated complex_taylor map.
)"""
  );
  py::class_<Bmad::Twiss1Propagate, std::unique_ptr<Bmad::Twiss1Propagate>>(
      m,
      "Twiss1Propagate",
      "twiss1_propagate return type"
  )
      .def_readonly("twiss2", &Bmad::Twiss1Propagate::twiss2)
      .def_readonly("err", &Bmad::Twiss1Propagate::err)
      .def("__len__", [](const Bmad::Twiss1Propagate &) { return 2; })
      .def("__getitem__", [](const Bmad::Twiss1Propagate &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.twiss2);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "twiss1_propagate",
      &Bmad::twiss1_propagate,
      py::arg("twiss1"),
      py::arg("mat2"),
      py::arg("ele_key"),
      py::arg("length"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss1_propagate

Parameters
----------
twiss1 : TwissStruct
    Input Twiss parameters.

mat2 : 2D array of float (shape: 2,2)
    The transfer matrix.

ele_key : int
    quadrupole$, etc.

length : float
    Determines whether the phase is increasing or decreasing.

Returns
-------
twiss2 : TwissStruct
    Output Twiss parameters.

err : bool
    Set True if there is an error, false otherwise.
)"""
  );
  m.def(
      "twiss3_at_start",
      &Bmad::twiss3_at_start,
      py::arg("lat"),
      py::arg("err_flag"),
      py::arg("ix_branch") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss3_at_start (lat, error, ix_branch, tune3)

Subroutine to calculate the 3D twiss parameters of the three modes of the full 6D 1-turn transfer matrix.
This routine is for lattices with closed geometries. For open lattices see: twiss3_from_twiss2.

Note: The rf must be on for this calculation.

Parameters
----------
lat : LatStruct
    Lattice with

ix_branch : int, optional
    Branch index. 0 = default.

Returns
-------
tune3 : 1D array of float (shape: 3), optional
    Normal mode tunes
)"""
  );
  m.def(
      "twiss3_from_twiss2",
      &Bmad::twiss3_from_twiss2,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss3_from_twiss2 (ele)

Routine to calculate the 3D Twiss parameters given the 2D transverse Twiss parameters and some
longitudinal parameters.
Also see: twiss3_at_start

Parameters
----------
ele : EleStruct
    Lattice element at which the calculation is made.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element
)"""
  );
  m.def(
      "twiss3_propagate1",
      &Bmad::twiss3_propagate1,
      py::arg("ele1"),
      py::arg("ele2"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss3_propagate1 (ele1, ele2, err_flag)

Subroutine to propagate the twiss parameters using all three normal modes.
Subroutine from original mode3_mod.
)"""
  );
  m.def(
      "twiss3_propagate_all",
      &Bmad::twiss3_propagate_all,
      py::arg("lat"),
      py::arg("ix_branch") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss3_propagate_all (lat, ix_branch)

Subroutine to propagate the twiss parameters using all three normal modes.
Subroutine from original mode3_mod.

Parameters
----------
lat : LatStruct
    Lattice

ix_branch : int, optional
    : Branch index. 0 = default.
)"""
  );
  m.def(
      "twiss_and_track",
      py::overload_cast<
          LatStruct &,
          CoordArrayStructAlloc1D,
          std::optional<bool>,
          std::optional<bool>,
          std::optional<bool>>(&Bmad::twiss_and_track),
      py::arg("lat"),
      py::arg("orb_array"),
      py::arg("print_err") = py::none(),
      py::arg("calc_chrom") = py::none(),
      py::arg("use_particle_start") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss_and_track

This routine is an overloaded name for:
  Subroutine twiss_and_track_branch (lat, orb, status, ix_branch, print_err, calc_chrom, orb_start, use_particle_start)
  Subroutine twiss_and_track_all (lat, orb_array, status, print_err, calc_chrom, use_particle_start)

Routine to calculate the twiss parameters, transport matrices and orbit.

The essential difference between these two procedures is that
twiss_and_track_branch only does the main branch while twiss_and_track_all
does everything but the photon_fork elements.

Note: This is not necessarily the fastest way to do things since this
routine does the entire calculation from scratch.

For a circular ring: If the RF is on, the computed orbit will be the 6D closed orbit.
If the RF is off, the 4D transverse closed orbit using orbi(0)%vec(6) is computed.

For an open lattice, the orbit will be computed using orb(0) as
starting conditions.

If there is a problem the status argument settings are: in_stop_band$,
unstable$, non_symplectic$, in_stop_band$, non_symplectic$, xfer_mat_clac_failure$,
twiss_propagate_failure$, no_complete_orbit$, or no_closed_orbit$. Note: in_stop_band$, unstable$,
and non_symplectic$ refer to the 1-turn matrix which is computed with closed lattices.
For an open geometry branch, status = no_complete_orbit$ is for
where the particle is lost in tracking. A negative sign is used to differentiate an
error occuring in the first call to twiss_at_start from the second call to twiss_at_start.

If there is a problem in an open geometry branch, status argument setting is -N where N is the element
where the particle was lost in tracking (negative numbers are used here to avoid confusion with ok$
which is mapped to 1.

Parameters
----------
lat : LatStruct
    lattice.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lat with computed twiss parameters.

orb_array : 1D array of CoordArrayStruct
    Array of orbit arrays.
    As an output, orb_array: Used as the starting point for an open lattice.

print_err : bool, optional
    Default is True. If False, suppress error messages.

calc_chrom : bool, optional
    Default is False. If True, calculate the chromatic functions.

use_particle_start : bool, optional
    Default is False. If True use branch.particle_start for the starting orbit. Do not use both this and
    orb_start.

Returns
-------
status : int, optional
    Set ok$ if everything is OK and set to something else otherwise. See above for more details.
)"""
  );
  m.def(
      "twiss_and_track_at_s",
      &Bmad::twiss_and_track_at_s,
      py::arg("lat"),
      py::arg("s"),
      py::arg("ele_at_s") = py::none(),
      py::arg("orb") = py::none(),
      py::arg("orb_at_s") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("use_last") = py::none(),
      py::arg("compute_floor_coords") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss_and_track_at_s (lat, s, ele_at_s, orb, orb_at_s, ix_branch, err, use_last, compute_floor_coords)

Subroutine to return the twiss parameters and particle orbit at a
given longitudinal position.

When calculating the Twiss parameters, this routine assumes
that the lattice elements already contain the Twiss parameters calculated
for the ends of the elements.

Additionally, the orbit at the ends of the elements (contained in orb(:)) must be
precomputed when orb_at_s is present.

Precomputation of Twiss and orbit at the element ends may be done with the twiss_and_track routine.

See also:
  twiss_and_track_from_s_to_s
  twiss_and_track_intra_ele

Parameters
----------
lat : LatStruct
    Lattice.

s : float
    Longitudinal position. If s is negative the the position is taken to be lat.param.total_length - s.

ele_at_s : EleStruct, optional
    If the use_last argument is True, ele_at_s is taken to contain valid Twiss parameters stored from a
    previous call to this routine.
    This parameter is an input/output and is modified in-place.
    As an output, ele_at_s: Element structure holding the Twiss parameters.

orb : 1D array of CoordStruct, optional
    Orbit through the Lattice.

orb_at_s : CoordStruct, optional
    If the use_last argument is True, orb_at_s is taken to contain the valid orbit stored from a previous
    call.
    This parameter is an input/output and is modified in-place.
    As an output, orb_at_s: Particle position at the position s.

ix_branch : int, optional
    Branch index, Default is 0 (main lattice).

use_last : bool, optional
    If present and True, and if ele_at_s.s < s, then use ele_at_s and orb_at_s as the starting point for the
    present calculation. This can speed things up when the present s-position is in the middle of a long
    complicated element and the tracking (EG: Runge-Kutta) is slow.

compute_floor_coords : bool, optional
    If present and True then the global "floor" coordinates (without misalignments) will be calculated and put
    in ele_at_s.floor.

Returns
-------
err : bool, optional
    Set True if there is a problem in the calculation, False otherwise.
)"""
  );
  m.def(
      "twiss_and_track",
      py::overload_cast<
          LatStruct &,
          CoordStructAlloc1D,
          std::optional<int>,
          std::optional<bool>,
          std::optional<bool>,
          optional_ref<CoordStruct>,
          std::optional<bool>>(&Bmad::twiss_and_track),
      py::arg("lat"),
      py::arg("orb"),
      py::arg("ix_branch") = py::none(),
      py::arg("print_err") = py::none(),
      py::arg("calc_chrom") = py::none(),
      py::arg("orb_start") = py::none(),
      py::arg("use_particle_start") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine twiss_and_track

This routine is an overloaded name for:
  Subroutine twiss_and_track_branch (lat, orb, status, ix_branch, print_err, calc_chrom, orb_start, use_particle_start)
  Subroutine twiss_and_track_all (lat, orb_array, status, print_err, calc_chrom, use_particle_start)

Routine to calculate the twiss parameters, transport matrices and orbit.

The essential difference between these two procedures is that
twiss_and_track_branch only does the main branch while twiss_and_track_all
does everything but the photon_fork elements.

Note: This is not necessarily the fastest way to do things since this
routine does the entire calculation from scratch.

For a circular ring: If the RF is on, the computed orbit will be the 6D closed orbit.
If the RF is off, the 4D transverse closed orbit using orbi(0)%vec(6) is computed.

For an open lattice, the orbit will be computed using orb(0) as
starting conditions.

If there is a problem the status argument settings are: in_stop_band$,
unstable$, non_symplectic$, in_stop_band$, non_symplectic$, xfer_mat_clac_failure$,
twiss_propagate_failure$, no_complete_orbit$, or no_closed_orbit$. Note: in_stop_band$, unstable$,
and non_symplectic$ refer to the 1-turn matrix which is computed with closed lattices.
For an open geometry branch, status = no_complete_orbit$ is for
where the particle is lost in tracking. A negative sign is used to differentiate an
error occuring in the first call to twiss_at_start from the second call to twiss_at_start.

If there is a problem in an open geometry branch, status argument setting is -N where N is the element
where the particle was lost in tracking (negative numbers are used here to avoid confusion with ok$
which is mapped to 1.

Parameters
----------
lat : LatStruct
    lattice.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lat with computed twiss parameters.

orb : 1D array of CoordStruct
    Orbit to be computed
    As an output, orb: Initial conditions to be used for an open geometry lattices.
    As an output, orb: Energy at which the closed orbit is computed.
    This parameter is an input/output and is modified in-place.
    As an output, orb: Computed orbit.

ix_branch : int, optional
    Branch to track.

print_err : bool, optional
    Default is True. If False, suppress error messages.

calc_chrom : bool, optional
    Default is False. If True, calculate the chromatic functions.

orb_start : CoordStruct, optional
    If present, use this as the starting orbit.

use_particle_start : bool, optional
    Default is False. If True use branch.particle_start for the starting orbit. Do not use both this and
    orb_start.

Returns
-------
status : int, optional
    Set ok$ if everything is OK and set to something else otherwise. See above for more details.
)"""
  );
  py::class_<Bmad::TwissAndTrackFromSToS, std::unique_ptr<Bmad::TwissAndTrackFromSToS>>(
      m,
      "TwissAndTrackFromSToS",
      "twiss_and_track_from_s_to_s return type"
  )
      .def_readonly("orbit_end", &Bmad::TwissAndTrackFromSToS::orbit_end)
      .def_readonly("ele_end", &Bmad::TwissAndTrackFromSToS::ele_end)
      .def_readonly("err", &Bmad::TwissAndTrackFromSToS::err)
      .def("__len__", [](const Bmad::TwissAndTrackFromSToS &) { return 3; })
      .def("__getitem__", [](const Bmad::TwissAndTrackFromSToS &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.orbit_end);
        if (i == 1)
          return py::cast(s.ele_end);
        if (i == 2)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "twiss_and_track_from_s_to_s",
      &Bmad::twiss_and_track_from_s_to_s,
      py::arg("branch"),
      py::arg("orbit_start"),
      py::arg("s_end"),
      py::arg("ele_start") = py::none(),
      py::arg("compute_floor_coords") = py::none(),
      py::arg("compute_twiss") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_and_track_from_s_to_s

Parameters
----------
branch : BranchStruct
    Lattice branch to track through.

orbit_start : CoordStruct
    Starting phase space coordinates at s_start.

s_end : float
    Ending position.

ele_start : EleStruct, optional
    Holds the starting parameters at s_start.

compute_floor_coords : bool, optional
    If present and True then the global "floor" coordinates will be calculated and put in ele_end.floor.

compute_twiss : bool, optional
    Default True. If False, to save a little time, do not compute Twiss parameters.

Returns
-------
orbit_end : CoordStruct
    End phase space coordinates.

ele_end : EleStruct, optional
    Holds the ending Twiss parameters and the transfer matrix. If present then the ele_start argument must
    also be present.

err : bool, optional
    Set True if there is a problem like the particle gets lost in tracking
)"""
  );
  py::class_<Bmad::TwissAndTrackIntraEle, std::unique_ptr<Bmad::TwissAndTrackIntraEle>>(
      m,
      "TwissAndTrackIntraEle",
      "twiss_and_track_intra_ele return type"
  )
      .def_readonly("orbit_end", &Bmad::TwissAndTrackIntraEle::orbit_end)
      .def_readonly("err", &Bmad::TwissAndTrackIntraEle::err)
      .def("__len__", [](const Bmad::TwissAndTrackIntraEle &) { return 2; })
      .def("__getitem__", [](const Bmad::TwissAndTrackIntraEle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.orbit_end);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "twiss_and_track_intra_ele",
      &Bmad::twiss_and_track_intra_ele,
      py::arg("ele"),
      py::arg("param"),
      py::arg("l_start"),
      py::arg("l_end"),
      py::arg("track_upstream_end"),
      py::arg("track_downstream_end"),
      py::arg("orbit_start") = py::none(),
      py::arg("ele_start") = py::none(),
      py::arg("ele_end") = py::none(),
      py::arg("compute_floor_coords") = py::none(),
      py::arg("compute_twiss") = py::none(),
      py::arg("reuse_ele_end") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_and_track_intra_ele

Parameters
----------
ele : EleStruct
    Element to track through.

param : LatParamStruct

l_start : float
    Start position measured from the beginning of the element.

l_end : float
    Stop position measured from the beginning of the element.

track_upstream_end : bool
    If True then entrance effects are included in the tracking. But only if l_start = 0 and
    orbit_start.location /= inside$.

track_downstream_end : bool
    If True then exit effects are included in the tracking but only if l_end = ele.value(l$) (within
    bmad_com.significant_length tol)

orbit_start : CoordStruct, optional
    Starting phase space coordinates at l_start.

ele_start : EleStruct, optional
    Holds the starting Twiss parameters at l_start.

ele_end : EleStruct, optional
    If reuse_ele_end is set True then reuse ele_end from trancking instead of recomputing ele_end from
    scratch. This can save time.
    This parameter is an input/output and is modified in-place.
    As an output, ele_end: Holds the ending Twiss parameters at l_end (except for photons).

compute_floor_coords : bool, optional
    If present and True then the global "floor" coordinates (without misalignments) will be calculated and put
    in ele_end.floor.

compute_twiss : bool, optional
    Default True. If False, to save a little time, do not compute Twiss parameters. Also if ele_start is not
    present, no Twiss parameters are computed.

reuse_ele_end : bool, optional
    If present and True, and if ele_end has the correct lonigitudianal length and key type, reuse ele_end from
    trancking instead of recomputing ele_end from scratch. This can save time.

Returns
-------
orbit_end : CoordStruct, optional
    End phase space coordinates. If present then the orbit_start argument must also be present.

err : bool, optional
    Set True if there is a problem like the particle gets lost in tracking
)"""
  );
  py::class_<Bmad::TwissAtElement, std::unique_ptr<Bmad::TwissAtElement>>(
      m,
      "TwissAtElement",
      "twiss_at_element return type"
  )
      .def_readonly("start", &Bmad::TwissAtElement::start)
      .def_readonly("end", &Bmad::TwissAtElement::end)
      .def_readonly("average", &Bmad::TwissAtElement::average)
      .def("__len__", [](const Bmad::TwissAtElement &) { return 3; })
      .def("__getitem__", [](const Bmad::TwissAtElement &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.start);
        if (i == 1)
          return py::cast(s.end);
        if (i == 2)
          return py::cast(s.average);
        throw py::index_error();
      });
  m.def(
      "twiss_at_element",
      &Bmad::twiss_at_element,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_at_element

Parameters
----------
ele : EleStruct
    Element to be averaged

Returns
-------
start : EleStruct, optional
    Twiss and s at start of element.

end : EleStruct, optional
    Twiss and s at end of element.

average : EleStruct, optional
    Average Twiss and s of element.
)"""
  );
  m.def(
      "twiss_at_start",
      &Bmad::twiss_at_start,
      py::arg("lat"),
      py::arg("ix_branch") = py::none(),
      py::arg("type_out") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_at_start

Parameters
----------
lat : LatStruct
    Lat
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with twiss parameters computed.

ix_branch : int, optional
    Branch to use. Default is 0 (main branch).

type_out : bool, optional
    If True (the default), print an error message If the 1-turn matrix is unstable.

Returns
-------
status : int, optional
    Calculation status: ok$, in_stop_band$, unstable$, or non_symplectic$
)"""
  );
  py::class_<Bmad::TwissFromTracking, std::unique_ptr<Bmad::TwissFromTracking>>(
      m,
      "TwissFromTracking",
      "twiss_from_tracking return type"
  )
      .def_readonly("symp_err", &Bmad::TwissFromTracking::symp_err)
      .def_readonly("err_flag", &Bmad::TwissFromTracking::err_flag)
      .def("__len__", [](const Bmad::TwissFromTracking &) { return 2; })
      .def("__getitem__", [](const Bmad::TwissFromTracking &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.symp_err);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "twiss_from_tracking",
      &Bmad::twiss_from_tracking,
      py::arg("lat"),
      py::arg("ref_orb0"),
      py::arg("d_orb") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_from_tracking

Parameters
----------
lat : LatStruct
    Lat to track through.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Structure holding the Twiss parameters.

ref_orb0 : CoordStruct
    Reference orbit at lat.ele(0).

d_orb : 1D array of float, optional
    Vector of offsets to use. If not present or zero bmad_com.d_orb(:) will be used.

Returns
-------
symp_err : float
    A measure of how symplectic the constructed matrices were before symplecitification. mat_symp_check for
    more details.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "twiss_propagate1",
      &Bmad::twiss_propagate1,
      py::arg("ele1"),
      py::arg("ele2"),
      py::arg("forward") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_propagate1

Parameters
----------
ele1 : EleStruct
    Element holding the starting Twiss parameters for forwards propagation.
    This parameter is an input/output and is modified in-place.
    As an output, ele1: Element for the ending Twiss parameters for backwards propagation.

ele2 : EleStruct
    Element holding the transfer matrix and, if backwards propagation, the starting Twiss.
    This parameter is an input/output and is modified in-place.
    As an output, ele2: Element for the ending Twiss parameters for forward propagation.

forward : bool, optional
    Default is True. If false, propagate the Twiss backwards.

Returns
-------
err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "twiss_propagate_all",
      &Bmad::twiss_propagate_all,
      py::arg("lat"),
      py::arg("ix_branch") = py::none(),
      py::arg("ie_start") = py::none(),
      py::arg("ie_end") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_propagate_all

Parameters
----------
lat : LatStruct
    lattice.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with parameters computed for the branch.

ix_branch : int, optional
    Branch index. Default is 0 (main lattice).

ie_start : int, optional
    Starting element index. Default is 0. Note: The first element at which the Twiss parameters are calculated
    is ie_start+1.

ie_end : int, optional
    Ending element index, Default is branch.n_ele_track.

Returns
-------
err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "twiss_to_1_turn_mat",
      &Bmad::twiss_to_1_turn_mat,
      py::arg("twiss"),
      py::arg("phi"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine twiss_to_1_turn_mat

Parameters
----------
twiss : TwissStruct
    Structure holding the Twiss parameters. .beta .alpha

phi : float
    Tune in radians.

Returns
-------
mat2 : 2D array of float (shape: 2,2)
    1-turn matrix.
)"""
  );
  py::class_<Bmad::TypeComplexTaylors, std::unique_ptr<Bmad::TypeComplexTaylors>>(
      m,
      "TypeComplexTaylors",
      "type_complex_taylors return type"
  )
      .def_readonly("lines", &Bmad::TypeComplexTaylors::lines)
      .def_readonly("n_lines", &Bmad::TypeComplexTaylors::n_lines)
      .def("__len__", [](const Bmad::TypeComplexTaylors &) { return 2; })
      .def("__getitem__", [](const Bmad::TypeComplexTaylors &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.lines);
        if (i == 1)
          return py::cast(s.n_lines);
        throw py::index_error();
      });
  m.def(
      "type_complex_taylors",
      &Bmad::type_complex_taylors,
      py::arg("complex_taylor"),
      py::arg("max_order") = py::none(),
      py::arg("file_id") = py::none(),
      py::arg("out_type") = py::none(),
      py::arg("clean") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine type_complex_taylors (complex_taylor, max_order, lines, n_lines, file_id, out_type, clean)

Subroutine to print or put in a string array a Bmad taylor map.
If the lines(:) argument is not present, the element information is printed to the terminal.

Parameters
----------
complex_taylor : 1D array of ComplexTaylorStruct
    Array of taylors.

max_order : int, optional
    Maximum order to print.

file_id : int, optional
    If present, write output to a file with handle file_id.

out_type : str, optional
    Determins the string to be used for the output type column. '' (default)  -> '1', '2', '3', etc. 'PHASE'
    -> 'X', 'Px, 'Y', 'Py', 'Z', 'Pz' 'SPIN'        -> 'S1', 'Sx', 'Sy', 'Sz'  ! If size = 4: quaternion
    representation 'SPIN'        -> 'Sx', 'Sy', 'Sz'  ! If size = 3: 'NONE'        -> No out column Anything
    else -> Use this for the output column.

clean : bool, optional
    If True then do not include terms whose coefficients are negligible. Default is false

Returns
-------
lines : 1D array of str, optional
    : Character array to hold the output. If not present, the information is printed to the terminal. Char
    width should be 120 or above for out_type = 'PHASE' but can be less for other out_types.

n_lines : int, optional
    Number of lines in lines(:) that hold valid output. n_lines must be present if lines(:) is.
)"""
  );
  m.def(
      "type_coord",
      &Bmad::type_coord,
      py::arg("coord"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine type_coord

Parameters
----------
coord : CoordStruct
    Coordinate
)"""
  );
  py::class_<Bmad::TypeEle, std::unique_ptr<Bmad::TypeEle>>(m, "TypeEle", "type_ele return type")
      .def_readonly("lines", &Bmad::TypeEle::lines)
      .def_readonly("n_lines", &Bmad::TypeEle::n_lines)
      .def("__len__", [](const Bmad::TypeEle &) { return 2; })
      .def("__getitem__", [](const Bmad::TypeEle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.lines);
        if (i == 1)
          return py::cast(s.n_lines);
        throw py::index_error();
      });
  m.def(
      "type_ele",
      &Bmad::type_ele,
      py::arg("ele"),
      py::arg("type_zero_attrib") = py::none(),
      py::arg("type_mat6") = py::none(),
      py::arg("type_taylor") = py::none(),
      py::arg("twiss_out") = py::none(),
      py::arg("type_control") = py::none(),
      py::arg("type_wake") = py::none(),
      py::arg("type_floor_coords") = py::none(),
      py::arg("type_field") = py::none(),
      py::arg("type_wall") = py::none(),
      py::arg("type_rad_kick") = py::none(),
      py::arg("type_internal") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine type_ele

Parameters
----------
ele : EleStruct
    Element

type_zero_attrib : bool, optional
    If False then surpress printing of real attributes whose value is 0 or switch attributes that have their
    default value. Default is False.

type_mat6 : int, optional
    = 0   => Do not type ele.mat6 = 4   => Type 4X4 xy submatrix = 6   => Type full 6x6 matrix (Default)

type_taylor : bool, optional
    Print out taylor map terms? If ele.taylor is not allocated then this is ignored. Default is False.

twiss_out : int, optional
    Print the Twiss parameters at the element end? = 0         => Do not print the Twiss parameters = radians$
    => Print Twiss, phi in radians (Default). = degrees$  => Print Twiss, phi in degrees. = cycles$   => Print
    Twiss, phi in radians/2pi.

type_control : int, optional
    Print control status? If ele.branch.lat is not associated cannot print status info. = no$      => One line
    of info. = short$   => Almost all info except long knot point lists are truncated (default). = all$     =>
    Everything.

type_wake : bool, optional
    If True then print the long-range and short-range wakes information. If False then just print how many
    terms the wake has. Default is True. If ele.wake is not allocated then this is ignored.

type_floor_coords : bool, optional
    Default is False. If True then print the global ("floor") coordinates at the exit end of the element.

type_field : int, optional
    Print field maps? = no$      => One line of info (default). = short$   => Header info. No tables. = all$
    => Everything.

type_wall : bool, optional
    Default is False. If True, print wall info.

type_rad_kick : bool, optional
    Default is False. If True, print synch rad kick info.

type_internal : bool, optional
    Default is False. If True, print some internal parameters.

Returns
-------
lines : 1D array of str, optional
    Character array to hold the output. If not present, the information is printed to the terminal.

n_lines : int, optional
    Number of lines in lines(:) that hold valid output. n_lines must be present if lines(:) is.
)"""
  );
  m.def(
      "type_end_stuff",
      &Bmad::type_end_stuff,
      py::arg("li"),
      py::arg("nl"),
      py::arg("lines") = py::none(),
      py::arg("n_lines") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine type_end_stuff

Parameters
----------
li : 1D array of str

nl : int

lines : 1D array of str, optional

n_lines : int, optional
)"""
  );
  m.def(
      "type_expression_tree",
      &Bmad::type_expression_tree,
      py::arg("tree"),
      py::arg("indent") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine type_expression_tree (tree, indent)

Routine to print an expression tree in tree form.
Good for debugging.

Parameters
----------
tree : ExpressionTreeStruct
    Tree to print.

indent : int, optional
    Initial indent. Default is zero.
)"""
  );
  py::class_<Bmad::TypePtcFibre, std::unique_ptr<Bmad::TypePtcFibre>>(
      m,
      "TypePtcFibre",
      "type_ptc_fibre return type"
  )
      .def_readonly("lines", &Bmad::TypePtcFibre::lines)
      .def_readonly("n_lines", &Bmad::TypePtcFibre::n_lines)
      .def("__len__", [](const Bmad::TypePtcFibre &) { return 2; })
      .def("__getitem__", [](const Bmad::TypePtcFibre &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.lines);
        if (i == 1)
          return py::cast(s.n_lines);
        throw py::index_error();
      });
  m.def(
      "type_ptc_fibre",
      &Bmad::type_ptc_fibre,
      py::arg("ptc_fibre"),
      py::arg("print_coords") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine type_ptc_fibre (ptc_fibre, print_coords, lines, n_lines)

Routine to put information on a PTC fibre element into a string array.
If "lines" is not present, the information will be printed to the screen.

Parameters
----------
ptc_fibre : Fibre
    Fibre to type info of.

print_coords : bool, optional
    If True then print coordinate and patch information. Default is True.

Returns
-------
lines : 1D array of str, optional
    Character array to hold the output.

n_lines : int, optional
    Number of lines used in lines(:)
)"""
  );
  m.def(
      "type_ptc_layout",
      &Bmad::type_ptc_layout,
      py::arg("lay"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine type_ptc_layout (lay)

Subroutine to print the global information in a layout
)"""
  );
  py::class_<PyTypeTaylors, std::unique_ptr<PyTypeTaylors>>(
      m,
      "TypeTaylors",
      "type_taylors return type"
  )
      .def_readonly("n_lines", &PyTypeTaylors::n_lines)
      .def("__len__", [](const PyTypeTaylors &) { return 1; })
      .def("__getitem__", [](const PyTypeTaylors &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.n_lines);
        throw py::index_error();
      });
  m.def(
      "type_taylors",
      &python_type_taylors,
      py::arg("bmad_taylor"),
      py::arg("max_order") = py::none(),
      py::arg("lines") = py::none(),
      py::arg("n_lines") = py::none(),
      py::arg("file_id") = py::none(),
      py::arg("out_style") = py::none(),
      py::arg("clean") = py::none(),
      py::arg("out_var_suffix") = py::none(),
      py::arg("append") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine type_taylors

Parameters
----------
bmad_taylor : 1D array of TaylorStruct
    Array of taylors.

max_order : int, optional
    Maximum order to print.

lines : 1D array of str, optional
    Used with append = True. Output will start at n_lines+1
    This parameter is an input/output and is modified in-place.
    As an output, lines: Array to hold the output.

n_lines : int, optional
    Used with append = True. Output will start at n_lines+1.
    This parameter is an input/output and is modified in-place.
    As an output, n_lines: Number of lines in lines(

file_id : int, optional
    If present, write output to a file with handle file_id.

out_style : str, optional
    Determins the string to be used for the output type column. '' (default) -> 'X', 'Px, 'Y', 'Py', 'Z', 'Pz'
    If size(bmad_taylor) = 6 -> 'S1', 'Sx', 'Sy', 'Sz' (Spin quaternions) If size(bmad_taylor) = 4 -> '1',
    '2', '3', etc. Otherwiase 'NUMBER'     -> '1', '2', '3', etc. 'BMAD'       -> Bmad lattice file format.

clean : bool, optional
    If True then do not include terms whose coefficients are negligible. Default is false.

out_var_suffix : str, optional
    If out_style = 'SCIBMAD', out_var_suffix is used as the suffix of the variable holding the taylor map.
    Default is "z". For example, if "z" is the suffix then: Descriptor = "d_z", orbital map name = "v_z", ref
    orbit name = v0_z, and spin map name = "q_z".

append : bool, optional
    Default is False. If True, n_lines on input is the number of existing lines in lines(:) to save.

Returns
-------
n_lines : int, optional
    Used with append = True. Output will start at n_lines+1.
    This parameter is an input/output and is modified in-place.
    As an output, n_lines: Number of lines in lines(
)"""
  );
  m.def(
      "type_twiss",
      &Bmad::type_twiss,
      py::arg("ele"),
      py::arg("frequency_units") = py::none(),
      py::arg("compact_format") = py::none(),
      py::arg("lines") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine type_twiss

Parameters
----------
ele : EleStruct
    Element containing the Twiss parameters.

frequency_units : int, optional
    Units for phi: = radians$  => Type Twiss, use radians for phi (Default). = degrees$  => Type Twiss, use
    degrees for phi. = cycles$   => Type Twiss, use cycles (1 = 2pi) units.

compact_format : bool, optional
    If present and True then use a compact output form.

lines : 1D array of str, optional
    : Character array to hold the output. The string length should be at least 120 characters. 13 lines are
    needed for the verbose form. If not present, the information is printed to the terminal.

Returns
-------
n_lines : int, optional
    Number of lines in lines(:) that hold valid output. n_lines must be present if lines(:) is.
)"""
  );
}
