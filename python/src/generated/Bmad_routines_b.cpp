#include "pybmad/generated/Bmad_routines_b.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_b(nb::module_ &m) {
  nb::class_<Bmad::BbiKick>(m, "BbiKick", "bbi_kick return type")
      .def_ro("nk", &Bmad::BbiKick::nk)
      .def_ro("dnk", &Bmad::BbiKick::dnk)
      .def("__len__", [](const Bmad::BbiKick &) { return 2; })
      .def("__getitem__", [](const Bmad::BbiKick &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.nk);
        if (i == 1)
          return nb::cast(s.dnk);
        throw nb::index_error();
      });
  m.def(
      "bbi_kick",
      &Bmad::bbi_kick,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("sigma"),
      nb::arg("linear_kick") = nb::none(),
      R"""(Wrapper for Fortran routine bbi_kick

Parameters
----------
x : float
    X coordinate.

y : float
    Y coordinate.

sigma : 1D array of float (shape: 2)
    Beam (x,y) sigmas.

linear_kick : bool, optional
    Default False. If present and True, kick and dnk are computed using the extrapolated kick from the linear
    region.

Returns
-------
nk : 1D array of float (shape: 2)
    Normalized, dimensionless kick component. In terms of the the actual kick: nk = [kick_x / (xi_x * sigma_x
    / beta_x), kick_y / (xi_y * sigma_y / beta_y) nk = -4 * pi * [x/sigma_x, y/sigma_y] in the linear region

dnk : 2D array of float (shape: 2,2)
    derivatives of nk. EG: dnk(2,1) = nk(2)/dx
)"""
  );
  m.def(
      "bbi_slice_calc",
      &Bmad::bbi_slice_calc,
      nb::arg("ele"),
      nb::arg("n_slice"),
      nb::arg("z_slice"),
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
      nb::arg("sigma_mat"),
      nb::arg("tail_cut"),
      nb::arg("tau"),
      nb::arg("energy"),
      nb::arg("n_part"),
      nb::arg("species"),
      R"""(This is a sigma matrix based IBS calculation.
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
      nb::arg("beam1"),
      nb::arg("beam2"),
      R"""(Wrapper for Fortran routine beam_equal_beam

Parameters
----------
beam1 : BeamStruct

beam2 : BeamStruct
)"""
  );
  nb::class_<Bmad::BeamInitSetup>(m, "BeamInitSetup", "beam_init_setup return type")
      .def_ro("err_flag", &Bmad::BeamInitSetup::err_flag)
      .def_ro("beam_init_set", &Bmad::BeamInitSetup::beam_init_set)
      .def("__len__", [](const Bmad::BeamInitSetup &) { return 2; })
      .def("__getitem__", [](const Bmad::BeamInitSetup &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.beam_init_set);
        throw nb::index_error();
      });
  m.def(
      "beam_init_setup",
      [](BeamInitStruct &beam_init_in, EleStruct &ele, int species, NormalModesStruct *modes) {
        auto fn = static_cast<
            Bmad::
                BeamInitSetup (*)(BeamInitStruct &, EleStruct &, int, optional_ref<NormalModesStruct>)>(
            &Bmad::beam_init_setup
        );
        return fn(beam_init_in, ele, species, ptr_to_opt_ref(modes));
      },
      nb::arg("beam_init_in"),
      nb::arg("ele"),
      nb::arg("species"),
      nb::arg("modes") = nb::none(),
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
beam_init_set : BeamInitStruct
    See above.

err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  nb::class_<Bmad::BeamTilts>(m, "BeamTilts", "beam_tilts return type")
      .def_ro("angle_xy", &Bmad::BeamTilts::angle_xy)
      .def_ro("angle_xz", &Bmad::BeamTilts::angle_xz)
      .def_ro("angle_yz", &Bmad::BeamTilts::angle_yz)
      .def_ro("angle_xpz", &Bmad::BeamTilts::angle_xpz)
      .def_ro("angle_ypz", &Bmad::BeamTilts::angle_ypz)
      .def("__len__", [](const Bmad::BeamTilts &) { return 5; })
      .def("__getitem__", [](const Bmad::BeamTilts &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.angle_xy);
        if (i == 1)
          return nb::cast(s.angle_xz);
        if (i == 2)
          return nb::cast(s.angle_yz);
        if (i == 3)
          return nb::cast(s.angle_xpz);
        if (i == 4)
          return nb::cast(s.angle_ypz);
        throw nb::index_error();
      });
  m.def(
      "beam_tilts",
      &Bmad::beam_tilts,
      nb::arg("S"),
      R"""(Given a 6x6 matrix of second-order moments, this routine returns
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
      nb::arg("ele"),
      R"""(Routine to setup a fibre to handle the beambeam interaction.

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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("particle_at"),
      nb::arg("orb"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      nb::arg("track_spin") = nb::none(),
      R"""(Subroutine to track through the edge field of an sbend.
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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("orbit"),
      nb::arg("local_ref_frame"),
      nb::arg("calc_dfield") = nb::none(),
      nb::arg("calc_potential") = nb::none(),
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
      nb::arg("ele"),
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
      nb::arg("r_in") = nb::none(),
      R"""(Routine to convert a random number in the interval [0,1] to a photon energy.
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
      nb::arg("E_photon"),
      nb::arg("g_bend"),
      nb::arg("gamma"),
      R"""(Routine to find the integrated probability corresponding to emitting a photon
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
      nb::arg("E_rel"),
      R"""(Routine to return the normalized probability that a photon will be emitted in a bend with energy
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
      nb::arg("g_bend_x"),
      nb::arg("g_bend_y"),
      nb::arg("gamma"),
      nb::arg("E_min") = nb::none(),
      nb::arg("E_max") = nb::none(),
      nb::arg("E_integ_prob") = nb::none(),
      nb::arg("vert_angle_min") = nb::none(),
      nb::arg("vert_angle_max") = nb::none(),
      nb::arg("vert_angle_symmetric") = nb::none(),
      nb::arg("emit_probability") = nb::none(),
      R"""(                                        vert_angle_min, vert_angle_max, vert_angle_symmetric, emit_probability)

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
      nb::arg("g_bend_x"),
      nb::arg("g_bend_y"),
      nb::arg("E_rel"),
      nb::arg("gamma_phi"),
      R"""(Routine to set a photon's polarization.
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
      nb::arg("E_rel"),
      nb::arg("gamma"),
      nb::arg("r_in") = nb::none(),
      nb::arg("invert") = nb::none(),
      R"""(Routine to convert an integrated probability to a vertical angle for emitting a photon from a bend.
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
  nb::class_<Bmad::BendShift>(m, "BendShift", "bend_shift return type")
      .def_ro("w_mat", &Bmad::BendShift::w_mat)
      .def_ro("position2", &Bmad::BendShift::position2)
      .def("__len__", [](const Bmad::BendShift &) { return 2; })
      .def("__getitem__", [](const Bmad::BendShift &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.w_mat);
        if (i == 1)
          return nb::cast(s.position2);
        throw nb::index_error();
      });
  m.def(
      "bend_shift",
      &Bmad::bend_shift,
      nb::arg("position1"),
      nb::arg("g"),
      nb::arg("delta_s"),
      nb::arg("ref_tilt") = nb::none(),
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
position2 : FloorPositionStruct
    particle coordinates relative to the final frame.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix used in the transformation
)"""
  );
  m.def(
      "bend_vert_angle_integ_prob",
      &Bmad::bend_vert_angle_integ_prob,
      nb::arg("vert_angle"),
      nb::arg("E_rel"),
      nb::arg("gamma"),
      R"""(Routine to find the integrated probability corresponding to emitting a photon
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
      nb::arg("current"),
      nb::arg("alpha"),
      nb::arg("Energy"),
      nb::arg("sigma_p"),
      nb::arg("Vrf"),
      nb::arg("omega"),
      nb::arg("U0"),
      nb::arg("circ"),
      nb::arg("R"),
      nb::arg("L"),
      R"""(This is a frontend for get_bl_from_fwhm from longitudinal_profile_mod.
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
  nb::class_<Bmad::BmadParser>(m, "BmadParser", "bmad_parser return type")
      .def_ro("lat", &Bmad::BmadParser::lat)
      .def_ro("digested_read_ok", &Bmad::BmadParser::digested_read_ok)
      .def_ro("err_flag", &Bmad::BmadParser::err_flag)
      .def_ro("parse_lat", &Bmad::BmadParser::parse_lat)
      .def("__len__", [](const Bmad::BmadParser &) { return 4; })
      .def("__getitem__", [](const Bmad::BmadParser &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.lat);
        if (i == 1)
          return nb::cast(s.digested_read_ok);
        if (i == 2)
          return nb::cast(s.err_flag);
        if (i == 3)
          return nb::cast(s.parse_lat);
        throw nb::index_error();
      });
  m.def(
      "bmad_parser",
      &Bmad::bmad_parser,
      nb::arg("lat_file"),
      nb::arg("make_mats6") = nb::none(),
      nb::arg("use_line") = nb::none(),
      R"""(Wrapper for Fortran routine bmad_parser

Parameters
----------
lat_file : str
    Name of the input file.

make_mats6 : bool, optional
    Compute the 6x6 transport matrices for the Elements? Default is True. Do not set False unless you know
    what you are doing.

use_line : str, optional
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
      [](std::string lat_file,
         LatStruct &lat,
         std::optional<CoordStructArray1D> orbit,
         std::optional<bool> make_mats6,
         std::optional<bool> err_flag,
         LatStruct *parse_lat) {
        auto fn = static_cast<
            void (*)(std::string, LatStruct &, std::optional<CoordStructArray1D>, std::optional<bool>, std::optional<bool>, optional_ref<LatStruct>)>(
            &Bmad::bmad_parser2
        );
        return fn(lat_file, lat, orbit, make_mats6, err_flag, ptr_to_opt_ref(parse_lat));
      },
      nb::arg("lat_file"),
      nb::arg("lat"),
      nb::arg("orbit") = nb::none(),
      nb::arg("make_mats6") = nb::none(),
      nb::arg("err_flag") = nb::none(),
      nb::arg("parse_lat") = nb::none(),
      R"""(Wrapper for Fortran routine bmad_parser2

Parameters
----------
lat_file : str
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
      "bmad_parser_string_attribute_set",
      [](EleStruct &ele,
         std::string attrib_name,
         std::string delim,
         bool delim_found,
         ParserEleStruct *pele,
         std::optional<std::string> str_out) {
        auto fn = static_cast<
            void (*)(EleStruct &, std::string, std::string, bool, optional_ref<ParserEleStruct>, std::optional<std::string>)>(
            &Bmad::bmad_parser_string_attribute_set
        );
        return fn(ele, attrib_name, delim, delim_found, ptr_to_opt_ref(pele), str_out);
      },
      nb::arg("ele"),
      nb::arg("attrib_name"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("pele") = nb::none(),
      nb::arg("str_out") = nb::none(),
      R"""(Wrapper for Fortran routine bmad_parser_string_attribute_set

Parameters
----------
ele : EleStruct

attrib_name : str

delim : str

delim_found : bool

pele : ParserEleStruct, optional

str_out : str, optional
)"""
  );
  m.def(
      "bmad_patch_parameters_to_ptc",
      &Bmad::bmad_patch_parameters_to_ptc,
      nb::arg("ang"),
      nb::arg("exi"),
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
      nb::arg("branch1"),
      nb::arg("branch2"),
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
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine branch_name

Parameters
----------
branch : BranchStruct
    Lattice branch

Returns
-------
name : str
    Encoded name
)"""
  );
  m.def(
      "branch_to_ptc_m_u",
      &Bmad::branch_to_ptc_m_u,
      nb::arg("branch"),
      R"""(Subroutine to create a PTC layout from a Bmad lattice branch.
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
      nb::arg("bunch1"),
      nb::arg("bunch2"),
      R"""(Wrapper for Fortran routine bunch_equal_bunch

Parameters
----------
bunch1 : BunchStruct

bunch2 : BunchStruct
)"""
  );
}
