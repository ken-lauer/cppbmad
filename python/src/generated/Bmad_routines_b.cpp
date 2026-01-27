#include "pybmad/generated/Bmad_routines_b.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_b(py::module &m) {
  py::class_<Bmad::BbiKick, std::unique_ptr<Bmad::BbiKick>>(m, "BbiKick", "bbi_kick return type")
      .def_readonly("nk", &Bmad::BbiKick::nk)
      .def_readonly("dnk", &Bmad::BbiKick::dnk)
      .def("__len__", [](const Bmad::BbiKick &) { return 2; })
      .def("__getitem__", [](const Bmad::BbiKick &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.nk);
        if (i == 1)
          return py::cast(s.dnk);
        throw py::index_error();
      });
  m.def(
      "bbi_kick",
      &Bmad::bbi_kick,
      py::arg("x"),
      py::arg("y"),
      py::arg("sigma"),
      R"""(Wrapper for Fortran routine bbi_kick

Parameters
----------
x : float
    X coordinate.

y : float
    Y coordinate.

sigma : 1D array of float (shape: 2)
    Beam (x,y) sigmas.

Returns
-------
nk : 1D array of float (shape: 2)
    Normalized, dimensionless kick component. In terms of the the actual kick: nk = [kick_x / (xi_x * sigma_x
    / beta_x), kick_y / (xi_y * sigma_y / beta_y) nk = -4 * pi * [x/sigma_x, y/sigma_y] in the linear region

dnk : 2D array of float (shape: 2,2)
    derivatives of nk. EG: dnk(2,1) = dnk(2)/dy
)"""
  );
  m.def(
      "bbi_slice_calc",
      &Bmad::bbi_slice_calc,
      py::arg("ele"),
      py::arg("n_slice"),
      py::arg("z_slice"),
      R"""(Wrapper for Fortran routine bbi_slice_calc

Parameters
----------
ele : EleStruct
    beambeam element

n_slice : int
    Number of slices

z_slice : 1D array of float
    Array of slice positions 1:n_slice. zero padded for indexes greater than n_slice
)"""
  );
  m.def(
      "beam_envelope_ibs",
      &Bmad::beam_envelope_ibs,
      py::arg("sigma_mat"),
      py::arg("tail_cut"),
      py::arg("tau"),
      py::arg("energy"),
      py::arg("n_part"),
      py::arg("species"),
      R"""(Subroutine beam_envelope_ibs(sigma_mat, ibs_mat, tail_cut, tau, energy, n_part, species)

This is a sigma matrix based IBS calculation.
It takes the beam sigma matrix and returns a matrix with changes to the 2nd order
moments due to IBS.

Use ibs_mat to change the sigma matrix like this:
sigma_matrix_updated = sigma_matrix + ibs_mat*element_length
See subroutine transport_with_sr_and_ibs in this module.

Parameters
----------
sigma_mat : 2D array of float (shape: 6,6)
    beam sigma_matrix at element entrance

tail_cut : bool
    If true, then apply tail cut to coulomb logarithm.

tau : float
    horizontal betatron damping rate.  Needed if tail_cut is true.

energy : float
    beam energy in eV

n_part : float
    number of particles in the bunch

species : int
    Partical species.

Returns
-------
ibs_mat : 2D array of float (shape: 6,6)
    changes in 2nd order moments due to IBS are ibs_mat*element_length
)"""
  );
  m.def(
      "beam_equal_beam",
      &Bmad::beam_equal_beam,
      py::arg("beam1"),
      py::arg("beam2"),
      R"""(Wrapper for Fortran routine beam_equal_beam

Parameters
----------
beam1 : BeamStruct

beam2 : BeamStruct
)"""
  );
  py::class_<Bmad::BeamInitSetup, std::unique_ptr<Bmad::BeamInitSetup>>(
      m,
      "BeamInitSetup",
      "beam_init_setup return type"
  )
      .def_readonly("err_flag", &Bmad::BeamInitSetup::err_flag)
      .def_readonly("beam_init_set", &Bmad::BeamInitSetup::beam_init_set)
      .def("__len__", [](const Bmad::BeamInitSetup &) { return 2; })
      .def("__getitem__", [](const Bmad::BeamInitSetup &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.beam_init_set);
        throw py::index_error();
      });
  m.def(
      "beam_init_setup",
      &Bmad::beam_init_setup,
      py::arg("beam_init_in"),
      py::arg("ele"),
      py::arg("species"),
      py::arg("modes") = py::none(),
      R"""(Wrapper for Fortran routine beam_init_setup

Parameters
----------
beam_init_in : BeamInitStruct
    Input parameters

ele : EleStruct

species : int
    Beam particle species.

modes : NormalModesStruct, optional
    Normal mode parameters.

Returns
-------
err_flag : bool, optional
    Set true if there is an error. False otherwise.

beam_init_set : BeamInitStruct
    See above.
)"""
  );
  py::class_<Bmad::BeamTilts, std::unique_ptr<Bmad::BeamTilts>>(
      m,
      "BeamTilts",
      "beam_tilts return type"
  )
      .def_readonly("angle_xy", &Bmad::BeamTilts::angle_xy)
      .def_readonly("angle_xz", &Bmad::BeamTilts::angle_xz)
      .def_readonly("angle_yz", &Bmad::BeamTilts::angle_yz)
      .def_readonly("angle_xpz", &Bmad::BeamTilts::angle_xpz)
      .def_readonly("angle_ypz", &Bmad::BeamTilts::angle_ypz)
      .def("__len__", [](const Bmad::BeamTilts &) { return 5; })
      .def("__getitem__", [](const Bmad::BeamTilts &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.angle_xy);
        if (i == 1)
          return py::cast(s.angle_xz);
        if (i == 2)
          return py::cast(s.angle_yz);
        if (i == 3)
          return py::cast(s.angle_xpz);
        if (i == 4)
          return py::cast(s.angle_ypz);
        throw py::index_error();
      });
  m.def(
      "beam_tilts",
      &Bmad::beam_tilts,
      py::arg("S"),
      R"""(Subroutine beam_tilts(S, angle_xy, angle_xz, angle_yz, angle_xpz, angle_ypz)

Given a 6x6 matrix of second-order moments, this routine returns
the beam tilts.

angle_xy is obtained from the projection of the beam envelop into the
xy plane.  The angle is that between the major axis of the projected
beam envelope and the +x axis.  Positive angles are measured towards the
+y axis.

angle_xz is obtained from the projection of the beam envelop into the
xy plane.  The angle is that between the major axis of the projected beam envelope
and the +z axis.  Positive angles are measured towards the +x axis.

angle_yz is obtained from the projection of the beam envelop into the
yz plane.  The angle is that between the major axis of the projected beam envelope
and the +z axis.  Positive angles are measured towards the +y axis.

Parameters
----------
S : 2D array of float (shape: 6,6)
    matrix of second order moments of beam envelope

Returns
-------
angle_xy : float
    transverse tilt of beam envelope

angle_xz : float
    horizontal crabbing of beam envelope

angle_yz : float
    vertical crabbing of beam envelope

angle_xpz : float
    x-pz coupling

angle_ypz : float
    y-pz coupling
)"""
  );
  m.def(
      "beambeam_fibre_setup",
      &Bmad::beambeam_fibre_setup,
      py::arg("ele"),
      R"""(Subroutine beambeam_fibre_setup(ele, ptc_fibre)

Routine to setup a fibre to handle the beambeam interaction.

Parameters
----------
ele : EleStruct
    Bmad beambeam element.

Returns
-------
ptc_fibre : Fibre
    Corresponding PTC fibre.
)"""
  );
  m.def(
      "bend_edge_kick",
      &Bmad::bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orb"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::arg("track_spin") = py::none(),
      R"""(Subroutine bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix, track_spin)

Subroutine to track through the edge field of an sbend.
This routine is called by apply_element_edge_kick only.

Parameters
----------
ele : EleStruct
    SBend element.

param : LatParamStruct
    Rel charge.

particle_at : int
    first_track_edge$, or second_track_edge$.

orb : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place.
    As an output, orb: Coords after tracking.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before fringe.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix transfer matrix including fringe.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

track_spin : bool, optional
    If True then track the spin through the edge fields. Default: False.
)"""
  );
  m.def(
      "bend_exact_multipole_field",
      &Bmad::bend_exact_multipole_field,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::arg("local_ref_frame"),
      py::arg("calc_dfield") = py::none(),
      py::arg("calc_potential") = py::none(),
      R"""(Wrapper for Fortran routine bend_exact_multipole_field

Parameters
----------
ele : EleStruct
    Bend element.

param : LatParamStruct
    Lattice branch parameters.

orbit : CoordStruct
    particle position.

local_ref_frame : bool
    Is the particle position in the local element ref frame (as opposed to the lab frame)?

calc_dfield : bool, optional
    If present and True then calculate the field derivatives.

calc_potential : bool, optional
    Calc electric and magnetic potentials? Default is false.

Returns
-------
field : EmFieldStruct
    Field
)"""
  );
  m.def(
      "bend_length_has_been_set",
      &Bmad::bend_length_has_been_set,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine bend_length_has_been_set

Parameters
----------
ele : EleStruct
    Element to be checked.

Returns
-------
is_set : bool
    Note: will be set True for non-bend elements.
)"""
  );
  m.def(
      "bend_photon_e_rel_init",
      &Bmad::bend_photon_e_rel_init,
      py::arg("r_in") = py::none(),
      R"""(Function bend_photon_e_rel_init (r_in) result (E_rel)

Routine to convert a random number in the interval [0,1] to a photon energy.
The photon probability spectrum is:
  P(E_rel) = (3 / (5 * Pi)) * Integral_{E_rel}^{Infty} K_{5/3}(x) dx
Where
  P(E_rel)) = Probability of finding a photon at relative energy E_rel.
  E_rel     = Relative photon energy: E / E_crit, E_crit = Critical energy.
  K_{5/3}   = Modified Bessel function.

Notice that the P(E) is not the same as the distribution radiation energy since
the photons must be energy weighted.

There is a cut-off built into the calculation so that E_rel will be in the
range [0, 31.4]. The error in neglecting photons with E_rel > 31.4 translates
to neglecting one photon for every 10^15 generated.
If r_in is present:
  r_in = 0 => E_rel = 0
  r_in = 1 => E_rel = 31.4

Parameters
----------
r_in : float, optional
    Integrated probability in the range [0,1]. If not present, a random number will be used.

Returns
-------
E_rel : float
    Relative photon energy E/E_crit.
)"""
  );
  m.def(
      "bend_photon_energy_integ_prob",
      &Bmad::bend_photon_energy_integ_prob,
      py::arg("E_photon"),
      py::arg("g_bend"),
      py::arg("gamma"),
      R"""(Function bend_photon_energy_integ_prob (E_photon, g_bend, gamma) result (integ_prob)

Routine to find the integrated probability corresponding to emitting a photon
from a bend in the range [0, E_photon].

Parameters
----------
E_photon : float
    Photon energy.

g_bend : float
    1/rho bending strength.

gamma : float
    Relativistic gamma factor of generating charged particle.

Returns
-------
integ_prob : float
    Integrated probability. Will be in the range [0, 1].
)"""
  );
  m.def(
      "bend_photon_energy_normalized_probability",
      &Bmad::bend_photon_energy_normalized_probability,
      py::arg("E_rel"),
      R"""(Function bend_photon_energy_normalized_probability (E_rel) result (prob)

Routine to return the normalized probability that a photon will be emitted in a bend with energy
E_rel relative to the critical energy. The probability is normalized such that
  Integral[0,Infinity] dE_rel P(E_rel) = 1

Parameters
----------
E_rel : float
    Photon energy relative to the critical energy.

Returns
-------
prob : float
    Normalized probability.
)"""
  );
  m.def(
      "bend_photon_init",
      &Bmad::bend_photon_init,
      py::arg("g_bend_x"),
      py::arg("g_bend_y"),
      py::arg("gamma"),
      py::arg("E_min") = py::none(),
      py::arg("E_max") = py::none(),
      py::arg("E_integ_prob") = py::none(),
      py::arg("vert_angle_min") = py::none(),
      py::arg("vert_angle_max") = py::none(),
      py::arg("vert_angle_symmetric") = py::none(),
      py::arg("emit_probability") = py::none(),
      R"""(Subroutine bend_photon_init (g_bend_x, g_bend_y, gamma, orbit, E_min, E_max, E_integ_prob,
                                        vert_angle_min, vert_angle_max, vert_angle_symmetric, emit_probability)

Routine to initalize a photon for dipole bends and wigglers (but not undulators).
The photon is initialized using the standard formulas for bending radiation.

The energy of the photon is calculated in one of two ways:

  1) If E_integ_prob is present and non-negative, the photon energy E will be such that the integrated
      probability  [E_min, E] relative to the integrated probability in the range [E_min, E_max] is E_integ_prob.
      That is, E_integ_prob can be used to to give a set of photon energies equally spaced in terms of the
      integrated probability distribution.

  2) If E_integ_prob is not present, or is negative, the photon energy is chosen at random in
      the range [E_min, E_max].

An E_integ_prob of zero means that the generated photon will have energy E_min.
An E_integ_prob of one means that the generated photon will have energy E_max.

The photon's polarization, will have unit amplitude.

This routine assumes that the emitting charged particle is on-axis and moving in
the forward direction. To correct for the actual charged particle postion use the routine
  absolute_photon_position

Parameters
----------
g_bend_x : float
    Bending 1/rho component in horizontal plane.

g_bend_y : float
    Bending 1/rho component in vertical plane.

gamma : float
    Relativistic gamma factor of generating charged particle.

E_min : float, optional
    Minimum photon energy. Default is zero. Ignored if negative.

E_max : float, optional
    Maximum photon energy.  Default is Infinity. Ignored if negative. If non-positive then E_max will be taken
    to be Infinity.

E_integ_prob : float, optional
    , optional :: integrated energy probability. See above. If E_integ_prob is non-negative, it must be in the
    range [0, 1].

vert_angle_min : float, optional
    Minimum vertical angle to emit a photon. -pi/2 is used if argument not present or if argument is less than
    -pi/2.

vert_angle_max : float, optional
    Maximum vertical angle to emit a photon. pi/2 is used if argument not present or if argument is greater
    than pi/2.

vert_angle_symmetric : bool, optional
    Default is False. If True, photons will be emitted in the range [-vert_angle_max, -vert_angle_min] as well
    as the range [vert_angle_min, vert_angle_max]. In this case vert_angle_min/max must be positive.

emit_probability : float, optional
    Probability of emitting a photon in the range [E_min, E_max] or in the vertical angular range given. The
    probability is normalized so that the probability of emitting if no ranges are given is 1.

Returns
-------
orbit : CoordStruct
    Initialized photon.
)"""
  );
  m.def(
      "bend_photon_polarization_init",
      &Bmad::bend_photon_polarization_init,
      py::arg("g_bend_x"),
      py::arg("g_bend_y"),
      py::arg("E_rel"),
      py::arg("gamma_phi"),
      R"""(Subroutine bend_photon_polarization_init (g_bend_x, g_bend_y, E_rel, gamma_phi, orbit)

Routine to set a photon's polarization.
The photon's polarization will be either in the plane of the bend or out of the plane and
the magnitude will be 1.

Parameters
----------
g_bend_x : float
    Bending 1/rho component in horizontal plane.

g_bend_y : float
    Bending 1/rho component in vertical plane.

E_rel : float
    Relative photon energy E/E_crit.

gamma_phi : float
    gamma * phi where gamma is the beam relativistic factor and phi is the vertical photon angle (in radians).

Returns
-------
orbit : CoordStruct
    Photon coords
)"""
  );
  m.def(
      "bend_photon_vert_angle_init",
      &Bmad::bend_photon_vert_angle_init,
      py::arg("E_rel"),
      py::arg("gamma"),
      py::arg("r_in") = py::none(),
      py::arg("invert") = py::none(),
      R"""(Function bend_photon_vert_angle_init (E_rel, gamma, r_in, invert) result (phi)

Routine to convert an integrated probability to a vertical angle for emitting a photon from a bend.
The integrated probability is in the range [0,1] with 0 corresponding to a phi = -pi/2 and
integrated probability of 1 corresponding to phi = pi/2.

Parameters
----------
E_rel : float
    Relative photon energy E/E_crit.

gamma : float
    beam relativistic factor

r_in : float, optional
    Integrated probability in the range [0,1]. If not present, a random number will be used.

invert : bool, optional
    If True then take r_in as the inverse integrated probability with inverted probability = 1 - probability.
    This is useful to avoid round-off errors when for looking at the tail of the distribution where the
    integrated prob is very close to 1 and small deviations can have large effects. Default is False.

Returns
-------
phi : float
    The photon vertical emission angle (in radians). Note: phi is an increasing monotonic function of r_in.
)"""
  );
  py::class_<Bmad::BendShift, std::unique_ptr<Bmad::BendShift>>(
      m,
      "BendShift",
      "bend_shift return type"
  )
      .def_readonly("w_mat", &Bmad::BendShift::w_mat)
      .def_readonly("position2", &Bmad::BendShift::position2)
      .def("__len__", [](const Bmad::BendShift &) { return 2; })
      .def("__getitem__", [](const Bmad::BendShift &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.w_mat);
        if (i == 1)
          return py::cast(s.position2);
        throw py::index_error();
      });
  m.def(
      "bend_shift",
      &Bmad::bend_shift,
      py::arg("position1"),
      py::arg("g"),
      py::arg("delta_s"),
      py::arg("ref_tilt") = py::none(),
      R"""(Wrapper for Fortran routine bend_shift

Parameters
----------
position1 : FloorPositionStruct
    Position of particle in inital coordinate frame.

g : float
    Curvature (1/rho)

delta_s : float
    S-position of final frame relative to the initial frame.

ref_tilt : float, optional
    ref_tilt. Default: 0

Returns
-------
w_mat : 2D array of float (shape: 3,3), optional
    W matrix used in the transformation

position2 : FloorPositionStruct
    particle coordinates relative to the final frame.
)"""
  );
  m.def(
      "bend_vert_angle_integ_prob",
      &Bmad::bend_vert_angle_integ_prob,
      py::arg("vert_angle"),
      py::arg("E_rel"),
      py::arg("gamma"),
      R"""(Function bend_vert_angle_integ_prob (vert_angle, E_rel, gamma) result (integ_prob)

Routine to find the integrated probability corresponding to emitting a photon
from a bend and with relative energy E_rel in the vertical angle range [-pi/2, vert_angle/2].

Note: vert_angle is allowed to be out of the range [-pi/2, pi/2]. In this case, integ_prob
will be set to 0 or 1 as appropriate.

Parameters
----------
vert_angle : float
    Vertical angle.

E_rel : float
    Relative photon energy E/E_crit.

gamma : float
    Relativistic gamma factor of generating charged particle.

Returns
-------
integ_prob : float
    Integrated probability. Will be in the range [0, 1].
)"""
  );
  m.def(
      "bl_via_vlassov",
      &Bmad::bl_via_vlassov,
      py::arg("current"),
      py::arg("alpha"),
      py::arg("Energy"),
      py::arg("sigma_p"),
      py::arg("Vrf"),
      py::arg("omega"),
      py::arg("U0"),
      py::arg("circ"),
      py::arg("R"),
      py::arg("L"),
      R"""(Subroutine bl_via_vlassov(current,alpha,Energy,sigma_p,Vrf,omega,U0,circ,R,L,sigma_z)

This is a frontend for get_bl_from_fwhm from longitudinal_profile_mod.
See longitudinal_profile_mod for details.  In short, this implements a model of potential well distortion
based on the Vlassov equation which uses an effective Resistive, Inductive, and Capacitive impedance.

Parameters
----------
current : float
    Beam current in amps

alpha : float
    Momentum compaction

Energy : float
    beam energy

sigma_p : float
    energy spread

Vrf : float
    total RF voltage in Volts

omega : float
    rf frequency in radians/s

U0 : float
    energy loss per turn (eV)

circ : float
    circumpherence

R : float
    Resistive part of effective impedance

L : float
    Inductive part of effective impedance

Returns
-------
sigma_z : float
    Bunch length. FWHM/TwoRootTwoLogTwo from bunch profile
)"""
  );
  py::class_<Bmad::BmadParser, std::unique_ptr<Bmad::BmadParser>>(
      m,
      "BmadParser",
      "bmad_parser return type"
  )
      .def_readonly("lat", &Bmad::BmadParser::lat)
      .def_readonly("digested_read_ok", &Bmad::BmadParser::digested_read_ok)
      .def_readonly("err_flag", &Bmad::BmadParser::err_flag)
      .def_readonly("parse_lat", &Bmad::BmadParser::parse_lat)
      .def("__len__", [](const Bmad::BmadParser &) { return 4; })
      .def("__getitem__", [](const Bmad::BmadParser &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.lat);
        if (i == 1)
          return py::cast(s.digested_read_ok);
        if (i == 2)
          return py::cast(s.err_flag);
        if (i == 3)
          return py::cast(s.parse_lat);
        throw py::index_error();
      });
  m.def(
      "bmad_parser",
      &Bmad::bmad_parser,
      py::arg("lat_file"),
      py::arg("make_mats6") = py::none(),
      py::arg("use_line") = py::none(),
      R"""(Wrapper for Fortran routine bmad_parser

Parameters
----------
lat_file : character
    Name of the input file.

make_mats6 : bool, optional
    Compute the 6x6 transport matrices for the Elements? Default is True. Do not set False unless you know
    what you are doing.

use_line : character, optional
    If present and not blank, override the use statement in the lattice file and use use_line instead.

Returns
-------
lat : LatStruct
    Lat structure. See bmad_struct.f90 for more details.

digested_read_ok : bool, optional
    Set True if the digested file was successfully read. False otherwise.

err_flag : bool, optional
    Set true if there is an error, false otherwise. Note: err_flag does *not* include errors in lat_make_mat6
    since if there is a match element, there is an error raised since the Twiss parameters have not been set
    but this is expected.

parse_lat : LatStruct, optional
    List of elements used to construct the lattice. Useful if bmad_parser2 will be called. See bmad_parser2
    documentation.
)"""
  );
  m.def(
      "bmad_parser2",
      &Bmad::bmad_parser2,
      py::arg("lat_file"),
      py::arg("lat"),
      py::arg("orbit") = py::none(),
      py::arg("make_mats6") = py::none(),
      py::arg("err_flag") = py::none(),
      py::arg("parse_lat") = py::none(),
      R"""(Wrapper for Fortran routine bmad_parser2

Parameters
----------
lat_file : character
    Input file name.

lat : LatStruct
    lattice with existing layout.
    This parameter is an input/output and is modified in-place.
    As an output, lat: lattice with modifications.

orbit : 1D array of CoordStruct, optional
    closed orbit for when bmad_parser2 calls lat_make_mat6

make_mats6 : bool, optional
    Make the 6x6 transport matrices for then Elements? Default is True.

err_flag : bool, optional

parse_lat : LatStruct, optional
    Used by bmad_parser to pass to bmad_parser2 a list of elements that were defined in the lattice file but
    not used. This is useful in preventing errors being generated if group/overlay elements definded by
    lat_file refer to unused slaves in parse_lat.
)"""
  );
  m.def(
      "bmad_patch_parameters_to_ptc",
      &Bmad::bmad_patch_parameters_to_ptc,
      py::arg("ang"),
      py::arg("exi"),
      R"""(Wrapper for Fortran routine bmad_patch_parameters_to_ptc

Parameters
----------
ang : 1D array of float (shape: 3)

exi : 2D array of float (shape: 3,3)
)"""
  );
  m.def(
      "bp_set_ran_status",
      &Bmad::bp_set_ran_status,
      R"""(Wrapper for Fortran routine bp_set_ran_status
)"""
  );
  m.def(
      "branch_equal_branch",
      &Bmad::branch_equal_branch,
      py::arg("branch1"),
      py::arg("branch2"),
      R"""(Wrapper for Fortran routine branch_equal_branch

Parameters
----------
branch1 : BranchStruct

branch2 : BranchStruct
)"""
  );
  m.def(
      "branch_name",
      &Bmad::branch_name,
      py::arg("branch"),
      R"""(Wrapper for Fortran routine branch_name

Parameters
----------
branch : BranchStruct
    Lattice branch

Returns
-------
name : character
    Encoded name
)"""
  );
  m.def(
      "branch_to_ptc_m_u",
      &Bmad::branch_to_ptc_m_u,
      py::arg("branch"),
      R"""(Subroutine branch_to_ptc_m_u (branch)

Subroutine to create a PTC layout from a Bmad lattice branch.
Note: If lat_to_ptc_layout has already been setup, you should first do a
          call kill_ptc_layouts(lat)
This deallocates the pointers in PTC

Note: If not already done, before you call this routine you need to first call:
   call set_ptc (...)
[This is normally done in bmad_parser.]

Note: If a Bmad element is using a hard edge model (EG: RFcavity element), there
will be three corresponding PTC fibre elements: (drift, RF. drift) for example.
In this case, ele%ptc_fibre will be set to point to the last PTC fibre. That is the
exit end of ele will correspond to the exit end of ele%ptc_fibre.

Parameters
----------
branch : BranchStruct
    Input branch.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Pointers to generated layouts.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Pointer to PTC fibres
)"""
  );
  m.def(
      "bunch_equal_bunch",
      &Bmad::bunch_equal_bunch,
      py::arg("bunch1"),
      py::arg("bunch2"),
      R"""(Wrapper for Fortran routine bunch_equal_bunch

Parameters
----------
bunch1 : BunchStruct

bunch2 : BunchStruct
)"""
  );
}
