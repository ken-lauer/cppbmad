#include "pybmad/generated/Bmad_routines_o.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyOdeintBmadTime python_odeint_bmad_time(
    CoordStruct &orb,
    EleStruct &ele,
    LatParamStruct &param,
    int t_dir,
    double rf_time,
    TrackStruct *track = nullptr,
    std::optional<double> t_end = std::nullopt,
    EmFieldStruct *extra_field = nullptr
) {
  auto _result = Bmad::odeint_bmad_time(
      orb,
      ele,
      param,
      t_dir,
      rf_time,
      ptr_to_opt_ref(track),
      t_end,
      ptr_to_opt_ref(extra_field)
  );
  auto py_result{PyOdeintBmadTime{_result, rf_time}};
  return py_result;
}
PyOffsetParticle python_offset_particle(
    EleStruct &ele,
    bool set,
    CoordStruct &coord,
    std::optional<bool> set_tilt = std::nullopt,
    std::optional<bool> set_hvkicks = std::nullopt,
    std::optional<int> drift_to_edge = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> set_spin = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<double> time = std::nullopt
) {
  auto _result = Bmad::offset_particle(
      ele,
      set,
      coord,
      set_tilt,
      set_hvkicks,
      drift_to_edge,
      s_pos,
      set_spin,
      mat6,
      make_matrix,
      make_opt_ref(time)
  );
  auto py_result{PyOffsetParticle{_result, time}};
  return py_result;
}

void init_Bmad_routines_o(nb::module_ &m) {
  nb::class_<Bmad::OdeintBmad>(m, "OdeintBmad", "odeint_bmad return type")
      .def_ro("err_flag", &Bmad::OdeintBmad::err_flag)
      .def_ro("track", &Bmad::OdeintBmad::track)
      .def("__len__", [](const Bmad::OdeintBmad &) { return 2; })
      .def("__getitem__", [](const Bmad::OdeintBmad &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.track);
        throw nb::index_error();
      });
  m.def(
      "odeint_bmad",
      &Bmad::odeint_bmad,
      nb::arg("orbit"),
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("s1_body"),
      nb::arg("s2_body"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
Recipes.  See the NR book for more details.

Notice that this routine has an two tolerances:
  bmad_com%rel_tol_adaptive_tracking
  bmad_com%abs_tol_adaptive_tracking

Note: For elements where the reference energy is not constant (lcavity, etc.), and
with elements where the reference particle does not follow the reference trajectory (wigglers for example),
the calculation of z is "off" while the particle is inside the element. At the ends there is no problem.

Parameters
----------
orbit : CoordStruct
    Starting coords: (x, px, y, py, z, delta) in element body coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coords

ele : EleStruct
    Element to track through.

param : LatParamStruct
    Lattice parameters.

s1_body : float
    Starting point relative to physical entrance.

s2_body : float
    Ending point relative physical entrance.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix propagated through the element.

make_matrix : bool, optional
    If True then make the 6x6 transfer matrix.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise. Note: a particle getting lost, for example hitting an
    aperture, is *not* an error.

track : TrackStruct, optional
    Structure holding the track information.
)"""
  );
  nb::class_<PyOdeintBmadTime>(m, "OdeintBmadTime", "odeint_bmad_time return type")
      .def_ro("err_flag", &PyOdeintBmadTime::err_flag)
      .def_ro("dt_step", &PyOdeintBmadTime::dt_step)
      .def_ro("rf_time", &PyOdeintBmadTime::rf_time)
      .def("__len__", [](const PyOdeintBmadTime &) { return 3; })
      .def("__getitem__", [](const PyOdeintBmadTime &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.dt_step);
        if (i == 2)
          return nb::cast(s.rf_time);
        throw nb::index_error();
      });
  m.def(
      "odeint_bmad_time",
      &python_odeint_bmad_time,
      nb::arg("orb"),
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("t_dir"),
      nb::arg("rf_time"),
      nb::arg("track") = nb::none(),
      nb::arg("t_end") = nb::none(),
      nb::arg("extra_field") = nb::none(),
      R"""(Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
Recipes.  See the NR book for more details.

Tracking is done until the particle is lost or exits the element.

Parameters
----------
orb : CoordStruct
    Starting coords: (x, px, y, py, s, ps) [t-based]
    This parameter is an input/output and is modified in-place.
    As an output, orb: Ending coords

ele : EleStruct
    Element to track through.

param : LatParamStruct
    Beam parameters.

t_dir : int
    Direction of time travel = +/-1. Can be negative for patches. Will be -1 if element has a negative length.

rf_time : float
    Time relative to RF clock.
    This parameter is an input/output and is modified in-place.
    As an output, rf_time: Updated time.

track : TrackStruct, optional
    Structure holding the track information.

t_end : float, optional
    If present, maximum time to which the particle will be tracked. Used for tracking with given time steps.
    The time orb.t at which tracking stops may be less than this if the particle gets to the end of the
    element

extra_field : EmFieldStruct, optional
    Static field to be added to the element field. Eg used with space charge.

Returns
-------
rf_time : float
    Time relative to RF clock.
    This parameter is an input/output and is modified in-place.
    As an output, rf_time: Updated time.

err_flag : bool
    Set True if there is an error. False otherwise.

dt_step : float, optional
    Next RK time step that this tracker would take based on the error tolerance. Used by track_bunch_time.
)"""
  );
  nb::class_<PyOffsetParticle>(m, "OffsetParticle", "offset_particle return type")
      .def_ro("s_out", &PyOffsetParticle::s_out)
      .def_ro("spin_qrot", &PyOffsetParticle::spin_qrot)
      .def_ro("time", &PyOffsetParticle::time)
      .def("__len__", [](const PyOffsetParticle &) { return 3; })
      .def("__getitem__", [](const PyOffsetParticle &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.s_out);
        if (i == 1)
          return nb::cast(s.spin_qrot);
        if (i == 2)
          return nb::cast(s.time);
        throw nb::index_error();
      });
  m.def(
      "offset_particle",
      &python_offset_particle,
      nb::arg("ele"),
      nb::arg("set"),
      nb::arg("coord"),
      nb::arg("set_tilt") = nb::none(),
      nb::arg("set_hvkicks") = nb::none(),
      nb::arg("drift_to_edge") = nb::none(),
      nb::arg("s_pos") = nb::none(),
      nb::arg("set_spin") = nb::none(),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      nb::arg("time") = nb::none(),
      R"""(Wrapper for Fortran routine offset_particle

Parameters
----------
ele : EleStruct
    Element

set : bool
    T (= set$)   -> Translate from lab coords to the local element coords. F (= unset$) -> Translate back from
    element to lab coords.

coord : CoordStruct

set_tilt : bool, optional
    Default is True. T -> Rotate using ele.value(tilt$) and ele.value(roll$) for sbends. F -> Do not rotate

set_hvkicks : bool, optional
    Default is True. T -> Apply 1/2 any hkick or vkick.

drift_to_edge : int, optional
    no$             -> Do not propagate (drift) particle. no$ is default if s_pos is present. upstream_end$
    -> Propagate to upsteam edge. This is default if set = set$ and s_pos is not present. downstream_end$ ->
    Propagate to downsteam edge. This is default if set = unset$ and s_pos is not present. Note: "edge" is
    body edge if set = set$ and is laboratory (nominal non-misaligned) edge if set = unset$

s_pos : float, optional
    Longitudinal particle position: If set = set$: Relative to upstream end (in lab coords). If set = unset$:
    Relative to entrance end (in body coords).

set_spin : bool, optional
    Default if False. Rotate spin coordinates? Also bmad_com.spin_tracking_on must be T to rotate.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before off setting.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix transfer matrix after offsets applied.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

time : float, optional
    Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
    This parameter is an input/output and is modified in-place.
    As an output, time: Updated time.

Returns
-------
s_out : float, optional
    Longitudinal particle position. If set = set$: Relative to entrance end (in body coords). If set = unset$:
    Relative to upstream end (in lab coords).

spin_qrot : 1D array of float (shape: 0:3), optional
    Spin rotation quaternion

time : float, optional
    Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
    This parameter is an input/output and is modified in-place.
    As an output, time: Updated time.
)"""
  );
  m.def(
      "offset_photon",
      &Bmad::offset_photon,
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("set"),
      nb::arg("offset_position_only") = nb::none(),
      nb::arg("rot_mat") = nb::none(),
      R"""(Wrapper for Fortran routine offset_photon

Parameters
----------
ele : EleStruct
    Element

orbit : CoordStruct
    Coordinates of the particle.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Coordinates of particle.

set : bool
    T (= set$)   -> Translate from lab coords to element coords. F (= unset$) -> Translate from element coords
    to lab coords.

offset_position_only : bool, optional
    If present and True, only offset the position coordinates.

rot_mat : 2D array of float (shape: 3,3), optional
    Rotation matrix from starting coords to ending coords.
)"""
  );
  m.def(
      "one_turn_mat_at_ele",
      &Bmad::one_turn_mat_at_ele,
      nb::arg("ele"),
      nb::arg("phi_a"),
      nb::arg("phi_b"),
      R"""(Wrapper for Fortran routine one_turn_mat_at_ele

Parameters
----------
ele : EleStruct
    Reference element.

phi_a : float
    "a" mode tune in radians.

phi_b : float
    "b" mode tune in radians.

Returns
-------
mat4 : 2D array of float (shape: 4,4)
    1-Turn coupled matrix.
)"""
  );
  nb::class_<Bmad::OpenBinaryFile>(m, "OpenBinaryFile", "open_binary_file return type")
      .def_ro("iu", &Bmad::OpenBinaryFile::iu)
      .def_ro("iver", &Bmad::OpenBinaryFile::iver)
      .def_ro("is_ok", &Bmad::OpenBinaryFile::is_ok)
      .def("__len__", [](const Bmad::OpenBinaryFile &) { return 3; })
      .def("__getitem__", [](const Bmad::OpenBinaryFile &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.iu);
        if (i == 1)
          return nb::cast(s.iver);
        if (i == 2)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "open_binary_file",
      &Bmad::open_binary_file,
      nb::arg("file_name"),
      nb::arg("action"),
      nb::arg("r_name"),
      R"""(Routine to open a binary file for reading or writing.

Parameters
----------
file_name : str
    File to create.

action : str
    'READ' or 'WRITE'

r_name : str
    Calling routine name for error messages.

Returns
-------
iu : int
    Unit number of opened file.

iver : int
    Version number if action = 'READ'

is_ok : bool
    Open OK?
)"""
  );
  nb::class_<Bmad::OrbitAmplitudeCalc>(m, "OrbitAmplitudeCalc", "orbit_amplitude_calc return type")
      .def_ro("amp_a", &Bmad::OrbitAmplitudeCalc::amp_a)
      .def_ro("amp_b", &Bmad::OrbitAmplitudeCalc::amp_b)
      .def_ro("amp_na", &Bmad::OrbitAmplitudeCalc::amp_na)
      .def_ro("amp_nb", &Bmad::OrbitAmplitudeCalc::amp_nb)
      .def("__len__", [](const Bmad::OrbitAmplitudeCalc &) { return 4; })
      .def("__getitem__", [](const Bmad::OrbitAmplitudeCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.amp_a);
        if (i == 1)
          return nb::cast(s.amp_b);
        if (i == 2)
          return nb::cast(s.amp_na);
        if (i == 3)
          return nb::cast(s.amp_nb);
        throw nb::index_error();
      });
  m.def(
      "orbit_amplitude_calc",
      &Bmad::orbit_amplitude_calc,
      nb::arg("ele"),
      nb::arg("orb"),
      R"""(Wrapper for Fortran routine orbit_amplitude_calc

Parameters
----------
ele : EleStruct
    Element holding the Twiss parameters, dispersion and coupling info.

orb : CoordStruct
    Orbit coordinates at the exit end of ele.

Returns
-------
amp_a : float, optional
    a-mode amplitude

amp_b : float, optional
    b-mode amplitude

amp_na : float, optional
    a-mode, energy normalized, amplitude.

amp_nb : float, optional
    b-mode, energy normalized, amplitude.
)"""
  );
  m.def(
      "orbit_reference_energy_correction",
      &Bmad::orbit_reference_energy_correction,
      nb::arg("orbit"),
      nb::arg("p0c_new"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Wrapper for Fortran routine orbit_reference_energy_correction

Parameters
----------
orbit : CoordStruct
    Coordinates to correct.

p0c_new : float
    New reference momentum.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before correction.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix transfer matrix including correction.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "orbit_to_floor_phase_space",
      &Bmad::orbit_to_floor_phase_space,
      nb::arg("orbit"),
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine orbit_to_floor_phase_space

Parameters
----------
orbit : CoordStruct
    Particle orbit in local (not element) coordinates.

ele : EleStruct
    Lattice element particle is in.

Returns
-------
floor_phase_space : 1D array of float (shape: 6)
    Floor phase space
)"""
  );
  m.def(
      "orbit_to_local_curvilinear",
      &Bmad::orbit_to_local_curvilinear,
      nb::arg("orbit"),
      nb::arg("ele"),
      nb::arg("z_direction") = nb::none(),
      nb::arg("relative_to") = nb::none(),
      R"""(Wrapper for Fortran routine orbit_to_local_curvilinear

Parameters
----------
orbit : CoordStruct
    Particle orbit in laboratory (not body) coordinates.

ele : EleStruct
    Lattice element particle is in.

z_direction : int, optional
    Set to +1 or -1.  Z-direction of particle velocity relative to element z-axis. Default is ele.orientation
    * orbit.direction.

relative_to : int, optional
    not_set$ (default), upstream_end$, downstream_end$. If not_set$ then origin is at the entrance end.

Returns
-------
local_position : FloorPositionStruct
    Position in local coordinates.
)"""
  );
  nb::class_<Bmad::OrbitTooLarge>(m, "OrbitTooLarge", "orbit_too_large return type")
      .def_ro("param", &Bmad::OrbitTooLarge::param)
      .def_ro("is_too_large", &Bmad::OrbitTooLarge::is_too_large)
      .def("__len__", [](const Bmad::OrbitTooLarge &) { return 2; })
      .def("__getitem__", [](const Bmad::OrbitTooLarge &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.param);
        if (i == 1)
          return nb::cast(s.is_too_large);
        throw nb::index_error();
      });
  m.def(
      "orbit_too_large",
      &Bmad::orbit_too_large,
      nb::arg("orbit"),
      nb::arg("check_momentum") = nb::none(),
      R"""(Wrapper for Fortran routine orbit_too_large

Parameters
----------
orbit : CoordStruct
    Particle orbit.

check_momentum : bool, optional
    If True (default) check the momentum.

Returns
-------
is_too_large : bool
    True if orbit is too large. False otherwise.

param : LatParamStruct, optional
)"""
  );
  nb::class_<Bmad::OrderEvecsByNSimilarity>(
      m,
      "OrderEvecsByNSimilarity",
      "order_evecs_by_n_similarity return type"
  )
      .def_ro("evec", &Bmad::OrderEvecsByNSimilarity::evec)
      .def_ro("err_flag", &Bmad::OrderEvecsByNSimilarity::err_flag)
      .def("__len__", [](const Bmad::OrderEvecsByNSimilarity &) { return 2; })
      .def("__getitem__", [](const Bmad::OrderEvecsByNSimilarity &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.evec);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "order_evecs_by_n_similarity",
      &Bmad::order_evecs_by_n_similarity,
      nb::arg("eval"),
      nb::arg("mat_tunes"),
      nb::arg("Nmat"),
      R"""(This subroutine orderes the eigensystem such that Nmat.mat_symp_conj(N) is closest
to the identity.  Nmat is supplied externally.

Parameters
----------
eval : 1D array of complex (shape: 6)
    complex eigenvalues.

mat_tunes : 1D array of float (shape: 3)
    Three normal mode tunes, in radians.
    This parameter is an input/output and is modified in-place.
    As an output, mat_tunes: Ordered normal mode tunes, in radians.

Nmat : 2D array of float (shape: 6,6)
    Normalized, real eigen matrix from make_N.

Returns
-------
evec : 2D array of complex (shape: 6,6)
    complex eigenvectors arranged down columns.

err_flag : bool
    Set True if there is an error. False otherwise
)"""
  );
  m.def(
      "order_evecs_by_plane_dominance",
      &Bmad::order_evecs_by_plane_dominance,
      nb::arg("evec"),
      nb::arg("eval"),
      nb::arg("mat_tunes") = nb::none(),
      R"""(This subroutine orderes the eigensystem according to which modes dominate the horizontal,
vertical, and longitudinal planes.  This subroutine works well in machines
that are not strongly coupled.  In machines with strong coupling, where the relation
between the three eigenmodes a, b, c and the three lab coordinates x, y, z can change
through the machine, this subroutine will not provide consistent ordering.

Parameters
----------
evec : 2D array of complex (shape: 6,6)
    complex eigenvectors arranged down columns.
    This parameter is an input/output and is modified in-place.
    As an output, evec: Ordered complex eigenvectors.

eval : 1D array of complex (shape: 6)
    complex eigenvalues.
    This parameter is an input/output and is modified in-place.
    As an output, eval: Ordered complex eigenvalues.

mat_tunes : 1D array of float (shape: 3), optional
    Three normal mode tunes, in radians.
    This parameter is an input/output and is modified in-place.
    As an output, mat_tunes: Reordered same as evecs.
)"""
  );
  m.def(
      "order_evecs_by_tune",
      &Bmad::order_evecs_by_tune,
      nb::arg("evec"),
      nb::arg("eval"),
      nb::arg("mat_tunes"),
      nb::arg("abz_tunes"),
      R"""(This subroutine orders the eigensystem by matching the tunes of the eigensystem to
externally supplied tunes abz_tunes.  abz_tunes is in radians.

Parameters
----------
evec : 2D array of complex (shape: 6,6)
    complex eigenvectors arranged down columns.
    This parameter is an input/output and is modified in-place.
    As an output, evec: Ordered eigenvectors.

eval : 1D array of complex (shape: 6)
    complex eigenvalues.
    This parameter is an input/output and is modified in-place.
    As an output, eval: Ordered eigenvalues.

mat_tunes : 1D array of float (shape: 3)
    Three normal mode tunes, in radians.

abz_tunes : 1D array of float (shape: 3)
    Tunes to order eigensystem by.

Returns
-------
err_flag : bool
    Set to true if an error occured.
)"""
  );
  m.def(
      "order_particles_in_z",
      &Bmad::order_particles_in_z,
      nb::arg("bunch"),
      R"""(Routine to order the particles longitudinally in terms of decreasing %vec(5).
That is from large z (head of bunch) to small z.
Only live particles are ordered.

Parameters
----------
bunch : BunchStruct
    collection of particles.
)"""
  );
  m.def(
      "order_super_lord_slaves",
      &Bmad::order_super_lord_slaves,
      nb::arg("lat"),
      nb::arg("ix_lord"),
      R"""(Wrapper for Fortran routine order_super_lord_slaves

Parameters
----------
lat : LatStruct
    Lat.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lat with fixed controls.

ix_lord : int
    Index of lord element.
)"""
  );
  m.def(
      "osc_alloc_freespace_array",
      &Bmad::osc_alloc_freespace_array,
      nb::arg("nlo"),
      nb::arg("nhi"),
      nb::arg("npad"),
      R"""(Wrapper for Fortran routine osc_alloc_freespace_array

Parameters
----------
nlo : 1D array of int (shape: 3)

nhi : 1D array of int (shape: 3)

npad : 1D array of int (shape: 3)
)"""
  );
  m.def(
      "osc_alloc_image_array",
      &Bmad::osc_alloc_image_array,
      nb::arg("nlo"),
      nb::arg("nhi"),
      nb::arg("npad"),
      R"""(Wrapper for Fortran routine osc_alloc_image_array

Parameters
----------
nlo : 1D array of int (shape: 3)

nhi : 1D array of int (shape: 3)

npad : 1D array of int (shape: 3)
)"""
  );
  m.def(
      "osc_alloc_rectpipe_arrays",
      &Bmad::osc_alloc_rectpipe_arrays,
      nb::arg("nlo"),
      nb::arg("nhi"),
      nb::arg("npad"),
      R"""(Wrapper for Fortran routine osc_alloc_rectpipe_arrays

Parameters
----------
nlo : 1D array of int (shape: 3)

nhi : 1D array of int (shape: 3)

npad : 1D array of int (shape: 3)
)"""
  );
  m.def(
      "osc_getgrnpipe",
      &Bmad::osc_getgrnpipe,
      nb::arg("gam"),
      nb::arg("a"),
      nb::arg("b"),
      nb::arg("delta"),
      nb::arg("umin"),
      nb::arg("npad"),
      R"""(Wrapper for Fortran routine osc_getgrnpipe

Parameters
----------
gam : float

a : float

b : float

delta : 1D array of float (shape: 3)

umin : 1D array of float (shape: 3)

npad : 1D array of int (shape: 3)
)"""
  );
  m.def(
      "osc_read_rectpipe_grn",
      &Bmad::osc_read_rectpipe_grn,
      R"""(Wrapper for Fortran routine osc_read_rectpipe_grn
)"""
  );
  m.def(
      "osc_write_rectpipe_grn",
      &Bmad::osc_write_rectpipe_grn,
      nb::arg("apipe"),
      nb::arg("bpipe"),
      nb::arg("delta"),
      nb::arg("umin"),
      nb::arg("umax"),
      nb::arg("nlo"),
      nb::arg("nhi"),
      nb::arg("gamma"),
      R"""(Wrapper for Fortran routine osc_write_rectpipe_grn

Parameters
----------
apipe : float

bpipe : float

delta : 1D array of float (shape: 3)

umin : 1D array of float (shape: 3)

umax : 1D array of float (shape: 3)

nlo : 1D array of int (shape: 3)

nhi : 1D array of int (shape: 3)

gamma : float
)"""
  );
}
