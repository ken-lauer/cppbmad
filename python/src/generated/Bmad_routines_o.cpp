#include "pybmad/generated/Bmad_routines_o.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_o(py::module &m) {
  py::class_<Bmad::OdeintBmad, std::unique_ptr<Bmad::OdeintBmad>>(
      m,
      "OdeintBmad",
      "odeint_bmad return type"
  )
      .def_readonly("err_flag", &Bmad::OdeintBmad::err_flag)
      .def_readonly("track", &Bmad::OdeintBmad::track)
      .def("__len__", [](const Bmad::OdeintBmad &) { return 2; })
      .def("__getitem__", [](const Bmad::OdeintBmad &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.track);
        throw py::index_error();
      });
  m.def(
      "odeint_bmad",
      &Bmad::odeint_bmad,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("s1_body"),
      py::arg("s2_body"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine odeint_bmad (orbit, ele, param, s1_body, s2_body, err_flag, track, mat6, make_matrix)

Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
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
    This parameter is an input/output and is modified in-place. As an output: Ending coords
ele : EleStruct
    Element to track through.
param : LatParamStruct
    Lattice parameters.
s1_body : float
    Starting point relative to physical entrance.
s2_body : float
    Ending point relative physical entrance.
mat6 : float, optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix propagated
    through the element.
make_matrix : bool, optional
    If True then make the 6x6 transfer matrix.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise. Note: a particle getting lost, for example hitting an
    aperture, is *not* an error.
track : TrackStruct
    Structure holding the track information.
)"""
  );
  py::class_<Bmad::OdeintBmadTime, std::unique_ptr<Bmad::OdeintBmadTime>>(
      m,
      "OdeintBmadTime",
      "odeint_bmad_time return type"
  )
      .def_readonly("err_flag", &Bmad::OdeintBmadTime::err_flag)
      .def_readonly("dt_step", &Bmad::OdeintBmadTime::dt_step)
      .def("__len__", [](const Bmad::OdeintBmadTime &) { return 2; })
      .def("__getitem__", [](const Bmad::OdeintBmadTime &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.dt_step);
        throw py::index_error();
      });
  m.def(
      "odeint_bmad_time",
      &Bmad::odeint_bmad_time,
      py::arg("orb"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("t_dir"),
      py::arg("rf_time"),
      py::arg("track") = py::none(),
      py::arg("t_end") = py::none(),
      py::arg("extra_field") = py::none(),
      R"""(Subroutine odeint_bmad_time (orb, ele, param, t_dir, rf_time, err_flag, track, t_end, dt_step, extra_field)

Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
Recipes.  See the NR book for more details.

Tracking is done until the particle is lost or exits the element.

Parameters
----------
orb : CoordStruct
    Starting coords: (x, px, y, py, s, ps) [t-based]
    This parameter is an input/output and is modified in-place. As an output: Ending coords
ele : EleStruct
    Element to track through. .tracking_method -- Determines which subroutine to use to calculate the field.
    Note: BMAD does no supply em_field_custom. == custom$ then use em_field_custom /= custom$ then use
    em_field_standard
param : LatParamStruct
    Beam parameters.
t_dir : float
    Direction of time travel = +/-1. Can be negative for patches. Will be -1 if element has a negative length.
rf_time : float
    Time relative to RF clock.
    This parameter is an input/output and is modified in-place. As an output: Updated time.
track : TrackStruct, optional
    Structure holding the track information. .save_track   -- Logical: Set True if track is to be saved.
t_end : float, optional
    If present, maximum time to which the particle will be tracked. Used for tracking with given time steps.
    The time orb.t at which tracking stops may be less than this if the particle gets to the end of the
    element
extra_field : EmFieldStruct, optional
    Static field to be added to the element field. Eg used with space charge.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
dt_step : float
    Next RK time step that this tracker would take based on the error tolerance. Used by track_bunch_time.
)"""
  );
  py::class_<Bmad::OffsetParticle, std::unique_ptr<Bmad::OffsetParticle>>(
      m,
      "OffsetParticle",
      "offset_particle return type"
  )
      .def_readonly("s_out", &Bmad::OffsetParticle::s_out)
      .def_readonly("spin_qrot", &Bmad::OffsetParticle::spin_qrot)
      .def("__len__", [](const Bmad::OffsetParticle &) { return 2; })
      .def("__getitem__", [](const Bmad::OffsetParticle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.s_out);
        if (i == 1)
          return py::cast(s.spin_qrot);
        throw py::index_error();
      });
  m.def(
      "offset_particle",
      &Bmad::offset_particle,
      py::arg("ele"),
      py::arg("set"),
      py::arg("orbit"),
      py::arg("set_tilt") = py::none(),
      py::arg("set_hvkicks") = py::none(),
      py::arg("drift_to_edge") = py::none(),
      py::arg("s_pos") = py::none(),
      py::arg("set_spin") = py::none(),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::arg("time") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element
set : bool
    T (= set$)   -> Translate from lab coords to the local element coords. F (= unset$) -> Translate back from
    element to lab coords.
orbit : CoordStruct
    Coordinates of the particle.
    This parameter is an input/output and is modified in-place. As an output: Coordinates of particle.
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
s_out : float
    Longitudinal particle position. If set = set$: Relative to entrance end (in body coords). If set = unset$:
    Relative to upstream end (in lab coords).
set_spin : bool, optional
    Default if False. Rotate spin coordinates? Also bmad_com.spin_tracking_on must be T to rotate.
mat6 : float, optional
    Transfer matrix before off setting.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix transfer matrix
    after offsets applied.
make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
spin_qrot : float
    Spin rotation quaternion
time : float, optional
    Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
    This parameter is an input/output and is modified in-place. As an output: Updated time.
)"""
  );
  m.def(
      "offset_photon",
      &Bmad::offset_photon,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("set"),
      py::arg("offset_position_only") = py::none(),
      py::arg("rot_mat") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element
orbit : CoordStruct
    Coordinates of the particle.
    This parameter is an input/output and is modified in-place. As an output: Coordinates of particle.
set : bool
    T (= set$)   -> Translate from lab coords to element coords. F (= unset$) -> Translate from element coords
    to lab coords.
offset_position_only : bool, optional
    If present and True, only offset the position coordinates. -- Logical, optional: If present and True, only
    offset the position coordinates.
rot_mat : float, optional
    Rotation matrix from starting coords to ending coords.
)"""
  );
  m.def(
      "one_turn_mat_at_ele",
      &Bmad::one_turn_mat_at_ele,
      py::arg("ele"),
      py::arg("phi_a"),
      py::arg("phi_b"),
      R"""(Parameters
----------
ele : EleStruct
    Reference element. .a       -- "a" mode Twiss parameter structure. .b       -- "b" mode Twiss parameter
    structure. .c_mat   -- 2x2 C matrix. .gamma_c -- gamma associated with C matrix.
phi_a : float
    "a" mode tune in radians.
phi_b : float
    "b" mode tune in radians.
mat4 : float
    1-Turn coupled matrix.
)"""
  );
  py::class_<Bmad::OpenBinaryFile, std::unique_ptr<Bmad::OpenBinaryFile>>(
      m,
      "OpenBinaryFile",
      "open_binary_file return type"
  )
      .def_readonly("iu", &Bmad::OpenBinaryFile::iu)
      .def_readonly("iver", &Bmad::OpenBinaryFile::iver)
      .def_readonly("is_ok", &Bmad::OpenBinaryFile::is_ok)
      .def("__len__", [](const Bmad::OpenBinaryFile &) { return 3; })
      .def("__getitem__", [](const Bmad::OpenBinaryFile &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.iu);
        if (i == 1)
          return py::cast(s.iver);
        if (i == 2)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "open_binary_file",
      &Bmad::open_binary_file,
      py::arg("file_name"),
      py::arg("action"),
      py::arg("r_name"),
      R"""(Function open_binary_file (file_name, action, iu, r_name, iver) result (is_ok)

Routine to open a binary file for reading or writing.

Parameters
----------
file_name : unknown
    File to create.
action : unknown
    'READ' or 'WRITE'
r_name : unknown
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
  py::class_<Bmad::OrbitAmplitudeCalc, std::unique_ptr<Bmad::OrbitAmplitudeCalc>>(
      m,
      "OrbitAmplitudeCalc",
      "orbit_amplitude_calc return type"
  )
      .def_readonly("amp_a", &Bmad::OrbitAmplitudeCalc::amp_a)
      .def_readonly("amp_b", &Bmad::OrbitAmplitudeCalc::amp_b)
      .def_readonly("amp_na", &Bmad::OrbitAmplitudeCalc::amp_na)
      .def_readonly("amp_nb", &Bmad::OrbitAmplitudeCalc::amp_nb)
      .def("__len__", [](const Bmad::OrbitAmplitudeCalc &) { return 4; })
      .def("__getitem__", [](const Bmad::OrbitAmplitudeCalc &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.amp_a);
        if (i == 1)
          return py::cast(s.amp_b);
        if (i == 2)
          return py::cast(s.amp_na);
        if (i == 3)
          return py::cast(s.amp_nb);
        throw py::index_error();
      });
  m.def(
      "orbit_amplitude_calc",
      &Bmad::orbit_amplitude_calc,
      py::arg("ele"),
      py::arg("orb"),
      R"""(Parameters
----------
ele : EleStruct
    Element holding the Twiss parameters, dispersion and coupling info.
orb : CoordStruct
    Orbit coordinates at the exit end of ele.
amp_a : float
    a-mode amplitude
amp_b : float
    b-mode amplitude
amp_na : float
    a-mode, energy normalized, amplitude.
amp_nb : float
    b-mode, energy normalized, amplitude.
)"""
  );
  m.def(
      "orbit_reference_energy_correction",
      &Bmad::orbit_reference_energy_correction,
      py::arg("orbit"),
      py::arg("p0c_new"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Parameters
----------
orbit : CoordStruct
    Coordinates to correct.
p0c_new : float
    New reference momentum.
mat6 : float, optional
    Transfer matrix before correction.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix transfer matrix
    including correction.
make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "orbit_to_floor_phase_space",
      &Bmad::orbit_to_floor_phase_space,
      py::arg("orbit"),
      py::arg("ele"),
      R"""(Parameters
----------
orbit : CoordStruct
    Particle orbit in local (not element) coordinates.
ele : EleStruct
    Lattice element particle is in.
floor_phase_space : float
    Floor phase space
)"""
  );
  m.def(
      "orbit_to_local_curvilinear",
      &Bmad::orbit_to_local_curvilinear,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("z_direction") = py::none(),
      py::arg("relative_to") = py::none(),
      R"""(Parameters
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
local_position : FloorPositionStruct
    Position in local coordinates.
)"""
  );
  py::class_<Bmad::OrbitTooLarge, std::unique_ptr<Bmad::OrbitTooLarge>>(
      m,
      "OrbitTooLarge",
      "orbit_too_large return type"
  )
      .def_readonly("param", &Bmad::OrbitTooLarge::param)
      .def_readonly("is_too_large", &Bmad::OrbitTooLarge::is_too_large)
      .def("__len__", [](const Bmad::OrbitTooLarge &) { return 2; })
      .def("__getitem__", [](const Bmad::OrbitTooLarge &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.param);
        if (i == 1)
          return py::cast(s.is_too_large);
        throw py::index_error();
      });
  m.def(
      "orbit_too_large",
      &Bmad::orbit_too_large,
      py::arg("orbit"),
      py::arg("check_momentum") = py::none(),
      R"""(Parameters
----------
orbit : CoordStruct
    Particle orbit.
param : LatParamStruct
    .unstable_factor  -- Set if orbit is too large. Otherwise not set
check_momentum : bool, optional
    If True (default) check the momentum.
is_too_large : bool
    True if orbit is too large. False otherwise.
)"""
  );
  py::class_<Bmad::OrderEvecsByNSimilarity, std::unique_ptr<Bmad::OrderEvecsByNSimilarity>>(
      m,
      "OrderEvecsByNSimilarity",
      "order_evecs_by_n_similarity return type"
  )
      .def_readonly("evec", &Bmad::OrderEvecsByNSimilarity::evec)
      .def_readonly("err_flag", &Bmad::OrderEvecsByNSimilarity::err_flag)
      .def("__len__", [](const Bmad::OrderEvecsByNSimilarity &) { return 2; })
      .def("__getitem__", [](const Bmad::OrderEvecsByNSimilarity &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.evec);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "order_evecs_by_n_similarity",
      &Bmad::order_evecs_by_n_similarity,
      py::arg("eval"),
      py::arg("mat_tunes"),
      py::arg("Nmat"),
      R"""(Subroutine order_evecs_by_N_similarity(evec, eval, mat_tunes, Nmat, err_flag)

This subroutine orderes the eigensystem such that Nmat.mat_symp_conj(N) is closest
to the identity.  Nmat is supplied externally.

Parameters
----------
eval : complex
    complex eigenvalues.
evecr : complex
    complex eigenvectors arranged down columns.
mat_tunes : float
    Three normal mode tunes, in radians.
    This parameter is an input/output and is modified in-place. As an output: Ordered normal mode tunes, in
    radians.
Nmat : float
    Normalized, real eigen matrix from make_N.

Returns
-------
evec : complex
    complex eigenvectors arranged down columns.
err_flag : bool
    Set True if there is an error. False otherwise
)"""
  );
  m.def(
      "order_evecs_by_plane_dominance",
      &Bmad::order_evecs_by_plane_dominance,
      py::arg("evec"),
      py::arg("eval"),
      py::arg("mat_tunes") = py::none(),
      R"""(Subroutine order_evecs_by_plane_dominance(evec, eval, mat_tunes)

This subroutine orderes the eigensystem according to which modes dominate the horizontal,
vertical, and longitudinal planes.  This subroutine works well in machines
that are not strongly coupled.  In machines with strong coupling, where the relation
between the three eigenmodes a, b, c and the three lab coordinates x, y, z can change
through the machine, this subroutine will not provide consistent ordering.

Parameters
----------
eval : complex
    complex eigenvalues.
    This parameter is an input/output and is modified in-place. As an output: Ordered complex eigenvalues.
evec : complex
    complex eigenvectors arranged down columns.
    This parameter is an input/output and is modified in-place. As an output: Ordered complex eigenvectors.
mat_tunes : float, optional
    Three normal mode tunes, in radians.
    This parameter is an input/output and is modified in-place. As an output: Reordered same as evecs.
)"""
  );
  m.def(
      "order_evecs_by_tune",
      &Bmad::order_evecs_by_tune,
      py::arg("evec"),
      py::arg("eval"),
      py::arg("mat_tunes"),
      py::arg("abz_tunes"),
      R"""(Subroutine order_evecs_by_tune(evec, eval, mat_tunes, abz_tunes, err_flag)

This subroutine orders the eigensystem by matching the tunes of the eigensystem to
externally supplied tunes abz_tunes.  abz_tunes is in radians.

Parameters
----------
eval : complex
    complex eigenvalues.
    This parameter is an input/output and is modified in-place. As an output: Ordered eigenvalues.
evec : complex
    complex eigenvectors arranged down columns.
    This parameter is an input/output and is modified in-place. As an output: Ordered eigenvectors.
mat_tunes : float
    Three normal mode tunes, in radians.
abz_tunes : float
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
      py::arg("bunch"),
      R"""(Subroutine order_particles_in_z (bunch)

Routine to order the particles longitudinally in terms of decreasing %vec(5).
That is from large z (head of bunch) to small z.
Only live particles are ordered.

Parameters
----------
bunch : BunchStruct
    collection of particles. .particle(j).vec(5) -- Longitudinal position of j^th particle.
)"""
  );
  m.def(
      "order_super_lord_slaves",
      &Bmad::order_super_lord_slaves,
      py::arg("lat"),
      py::arg("ix_lord"),
      R"""(Parameters
----------
lat : LatStruct
    Lat with fixed controls.
ix_lord : int
    Index of lord element. Output
)"""
  );
  m.def(
      "osc_alloc_freespace_array",
      &Bmad::osc_alloc_freespace_array,
      py::arg("nlo"),
      py::arg("nhi"),
      py::arg("npad"),
      R"""(Parameters
----------
nlo : 
nhi : 
npad : 
)"""
  );
  m.def(
      "osc_alloc_image_array",
      &Bmad::osc_alloc_image_array,
      py::arg("nlo"),
      py::arg("nhi"),
      py::arg("npad"),
      R"""(Parameters
----------
nlo : 
nhi : 
npad : 
)"""
  );
  m.def(
      "osc_alloc_rectpipe_arrays",
      &Bmad::osc_alloc_rectpipe_arrays,
      py::arg("nlo"),
      py::arg("nhi"),
      py::arg("npad"),
      R"""(Parameters
----------
nlo : 
nhi : 
npad : 
)"""
  );
  m.def(
      "osc_getgrnpipe",
      &Bmad::osc_getgrnpipe,
      py::arg("gam"),
      py::arg("a"),
      py::arg("b"),
      py::arg("delta"),
      py::arg("umin"),
      py::arg("npad"),
      R"""(Parameters
----------
gam : 
a : 
b : 
delta : 
umin : 
npad : 
)"""
  );
  m.def("osc_read_rectpipe_grn", &Bmad::osc_read_rectpipe_grn, R"""()""");
  m.def(
      "osc_write_rectpipe_grn",
      &Bmad::osc_write_rectpipe_grn,
      py::arg("apipe"),
      py::arg("bpipe"),
      py::arg("delta"),
      py::arg("umin"),
      py::arg("umax"),
      py::arg("nlo"),
      py::arg("nhi"),
      py::arg("gamma"),
      R"""(Parameters
----------
apipe : 
bpipe : 
delta : 
umin : 
umax : 
nlo : 
nhi : 
gamma : 
)"""
  );
}
