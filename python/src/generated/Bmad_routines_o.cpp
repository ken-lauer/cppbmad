#include "pybmad/generated/Bmad_routines_o.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyOdeintBmadTime : public Bmad::OdeintBmadTime {
  double rf_time;
  PyOdeintBmadTime(Bmad::OdeintBmadTime _base, double rf_time)
      : Bmad::OdeintBmadTime(std::move(_base)), rf_time(rf_time) {}
};
PyOdeintBmadTime python_odeint_bmad_time(
    CoordProxy& orb,
    EleProxy& ele,
    LatParamProxy& param,
    int t_dir,
    double rf_time,
    optional_ref<TrackProxy> track = std::nullopt,
    std::optional<double> t_end = std::nullopt,
    optional_ref<EmFieldProxy> extra_field = std::nullopt) {
  auto _result = Bmad::odeint_bmad_time(
      orb, ele, param, t_dir, rf_time, track, t_end, extra_field);
  auto py_result{PyOdeintBmadTime{_result, rf_time}};
  return py_result;
}
struct PyOffsetParticle : public Bmad::OffsetParticle {
  std::optional<double> time;
  PyOffsetParticle(Bmad::OffsetParticle _base, std::optional<double> time)
      : Bmad::OffsetParticle(std::move(_base)), time(time) {}
};
PyOffsetParticle python_offset_particle(
    EleProxy& ele,
    bool set,
    CoordProxy& orbit,
    std::optional<bool> set_tilt = std::nullopt,
    std::optional<bool> set_hvkicks = std::nullopt,
    std::optional<int> drift_to_edge = std::nullopt,
    std::optional<double> s_pos = std::nullopt,
    std::optional<bool> set_spin = std::nullopt,
    std::optional<FixedArray2D<Real, 6, 6>> mat6 = std::nullopt,
    std::optional<bool> make_matrix = std::nullopt,
    std::optional<double> time = std::nullopt) {
  auto _result = Bmad::offset_particle(
      ele,
      set,
      orbit,
      set_tilt,
      set_hvkicks,
      drift_to_edge,
      s_pos,
      set_spin,
      mat6,
      make_matrix,
      make_opt_ref(time));
  auto py_result{PyOffsetParticle{_result, time}};
  return py_result;
}

struct PyOrbitTooLarge {
  LatParamProxy param;
  bool is_too_large;
};
PyOrbitTooLarge python_orbit_too_large(
    CoordProxy& orbit,
    std::optional<bool> check_momentum,
    bool is_too_large) {
  auto _result = Bmad::orbit_too_large(orbit, check_momentum, is_too_large);
  auto py_result{PyOrbitTooLarge{_result, is_too_large}};
  return py_result;
}
struct PyOscGetgrnpipe {
  double gam;
  double a;
  double b;
};
PyOscGetgrnpipe python_osc_getgrnpipe(
    double gam,
    double a,
    double b,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Int, 3> npad) {
  Bmad::osc_getgrnpipe(gam, a, b, delta, umin, npad);
  auto py_result{PyOscGetgrnpipe{gam, a, b}};
  return py_result;
}
struct PyOscWriteRectpipeGrn {
  double apipe;
  double bpipe;
  double gamma;
};
PyOscWriteRectpipeGrn python_osc_write_rectpipe_grn(
    double apipe,
    double bpipe,
    FixedArray1D<Real, 3> delta,
    FixedArray1D<Real, 3> umin,
    FixedArray1D<Real, 3> umax,
    FixedArray1D<Int, 3> nlo,
    FixedArray1D<Int, 3> nhi,
    double gamma) {
  Bmad::osc_write_rectpipe_grn(
      apipe, bpipe, delta, umin, umax, nlo, nhi, gamma);
  auto py_result{PyOscWriteRectpipeGrn{apipe, bpipe, gamma}};
  return py_result;
}

void init_Bmad_routines_o(py::module& m) {
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
  )""");
  py::class_<Bmad::OdeintBmad, std::unique_ptr<Bmad::OdeintBmad>>(
      m, "OdeintBmad", "Fortran routine odeint_bmad return value")
      .def_readonly("err_flag", &Bmad::OdeintBmad::err_flag)
      .def_readonly("track", &Bmad::OdeintBmad::track)
      .def("__len__", [](const Bmad::OdeintBmad&) { return 2; })
      .def("__getitem__", [](const Bmad::OdeintBmad& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.track);
        throw py::index_error();
      });
  m.def(
      "odeint_bmad_time",
      &python_odeint_bmad_time,
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
  )""");
  py::class_<PyOdeintBmadTime, std::unique_ptr<PyOdeintBmadTime>>(
      m, "OdeintBmadTime", "Fortran routine odeint_bmad_time return value")
      .def_readonly("err_flag", &PyOdeintBmadTime::err_flag)
      .def_readonly("dt_step", &PyOdeintBmadTime::dt_step)
      .def_readonly("rf_time", &PyOdeintBmadTime::rf_time)
      .def("__len__", [](const PyOdeintBmadTime&) { return 3; })
      .def("__getitem__", [](const PyOdeintBmadTime& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.dt_step);
        if (i == 2)
          return py::cast(s.rf_time);
        throw py::index_error();
      });
  m.def(
      "offset_particle",
      &python_offset_particle,
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
  )""");
  py::class_<PyOffsetParticle, std::unique_ptr<PyOffsetParticle>>(
      m, "OffsetParticle", "Fortran routine offset_particle return value")
      .def_readonly("s_out", &PyOffsetParticle::s_out)
      .def_readonly("spin_qrot", &PyOffsetParticle::spin_qrot)
      .def_readonly("time", &PyOffsetParticle::time)
      .def("__len__", [](const PyOffsetParticle&) { return 3; })
      .def("__getitem__", [](const PyOffsetParticle& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.s_out);
        if (i == 1)
          return py::cast(s.spin_qrot);
        if (i == 2)
          return py::cast(s.time);
        throw py::index_error();
      });
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
  )""");
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
  )""");
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
  )""");
  py::class_<Bmad::OpenBinaryFile, std::unique_ptr<Bmad::OpenBinaryFile>>(
      m, "OpenBinaryFile", "Fortran routine open_binary_file return value")
      .def_readonly("iu", &Bmad::OpenBinaryFile::iu)
      .def_readonly("iver", &Bmad::OpenBinaryFile::iver)
      .def_readonly("is_ok", &Bmad::OpenBinaryFile::is_ok)
      .def("__len__", [](const Bmad::OpenBinaryFile&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::OpenBinaryFile& s, int i) -> py::object {
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
  )""");
  py::class_<
      Bmad::OrbitAmplitudeCalc,
      std::unique_ptr<Bmad::OrbitAmplitudeCalc>>(
      m,
      "OrbitAmplitudeCalc",
      "Fortran routine orbit_amplitude_calc return value")
      .def_readonly("amp_a", &Bmad::OrbitAmplitudeCalc::amp_a)
      .def_readonly("amp_b", &Bmad::OrbitAmplitudeCalc::amp_b)
      .def_readonly("amp_na", &Bmad::OrbitAmplitudeCalc::amp_na)
      .def_readonly("amp_nb", &Bmad::OrbitAmplitudeCalc::amp_nb)
      .def("__len__", [](const Bmad::OrbitAmplitudeCalc&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::OrbitAmplitudeCalc& s, int i) -> py::object {
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
  )""");
  m.def(
      "orbit_to_floor_phase_space",
      &Bmad::orbit_to_floor_phase_space,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("floor_phase_space"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle orbit in local (not element) coordinates.
  ele : EleStruct
      Lattice element particle is in.
  floor_phase_space : 
  )""");
  m.def(
      "orbit_to_local_curvilinear",
      &Bmad::orbit_to_local_curvilinear,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("z_direction") = py::none(),
      py::arg("relative_to") = py::none(),
      py::arg("local_position"),
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
  local_position : 
  )""");
  m.def(
      "orbit_too_large",
      &python_orbit_too_large,
      py::arg("orbit"),
      py::arg("check_momentum") = py::none(),
      py::arg("is_too_large"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle orbit.
  param : LatParamStruct
      .unstable_factor  -- Set if orbit is too large. Otherwise not set
  check_momentum : bool, optional
      If True (default) check the momentum.
  is_too_large : 
  )""");
  py::class_<PyOrbitTooLarge, std::unique_ptr<PyOrbitTooLarge>>(
      m, "OrbitTooLarge", "Fortran routine orbit_too_large return value")
      .def_readonly("param", &PyOrbitTooLarge::param)
      .def_readonly("is_too_large", &PyOrbitTooLarge::is_too_large)
      .def("__len__", [](const PyOrbitTooLarge&) { return 2; })
      .def("__getitem__", [](const PyOrbitTooLarge& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.param);
        if (i == 1)
          return py::cast(s.is_too_large);
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
  )""");
  py::class_<
      Bmad::OrderEvecsByNSimilarity,
      std::unique_ptr<Bmad::OrderEvecsByNSimilarity>>(
      m,
      "OrderEvecsByNSimilarity",
      "Fortran routine order_evecs_by_n_similarity return value")
      .def_readonly("evec", &Bmad::OrderEvecsByNSimilarity::evec)
      .def_readonly("err_flag", &Bmad::OrderEvecsByNSimilarity::err_flag)
      .def("__len__", [](const Bmad::OrderEvecsByNSimilarity&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::OrderEvecsByNSimilarity& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.evec);
            if (i == 1)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
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
  )""");
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
  )""");
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
  )""");
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
  )""");
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
  )""");
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
  )""");
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
  )""");
  m.def(
      "osc_getgrnpipe",
      &python_osc_getgrnpipe,
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
  )""");
  py::class_<PyOscGetgrnpipe, std::unique_ptr<PyOscGetgrnpipe>>(
      m, "OscGetgrnpipe", "Fortran routine osc_getgrnpipe return value")
      .def_readonly("gam", &PyOscGetgrnpipe::gam)
      .def_readonly("a", &PyOscGetgrnpipe::a)
      .def_readonly("b", &PyOscGetgrnpipe::b)
      .def("__len__", [](const PyOscGetgrnpipe&) { return 3; })
      .def("__getitem__", [](const PyOscGetgrnpipe& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.gam);
        if (i == 1)
          return py::cast(s.a);
        if (i == 2)
          return py::cast(s.b);
        throw py::index_error();
      });
  m.def("osc_read_rectpipe_grn", &Bmad::osc_read_rectpipe_grn, R"""()""");
  m.def(
      "osc_write_rectpipe_grn",
      &python_osc_write_rectpipe_grn,
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
  )""");
  py::class_<PyOscWriteRectpipeGrn, std::unique_ptr<PyOscWriteRectpipeGrn>>(
      m,
      "OscWriteRectpipeGrn",
      "Fortran routine osc_write_rectpipe_grn return value")
      .def_readonly("apipe", &PyOscWriteRectpipeGrn::apipe)
      .def_readonly("bpipe", &PyOscWriteRectpipeGrn::bpipe)
      .def_readonly("gamma", &PyOscWriteRectpipeGrn::gamma)
      .def("__len__", [](const PyOscWriteRectpipeGrn&) { return 3; })
      .def(
          "__getitem__",
          [](const PyOscWriteRectpipeGrn& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.apipe);
            if (i == 1)
              return py::cast(s.bpipe);
            if (i == 2)
              return py::cast(s.gamma);
            throw py::index_error();
          });
}
