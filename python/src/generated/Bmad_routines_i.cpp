#include "pybmad/generated/Bmad_routines_i.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_i(nb::module_ &m) {
  m.def(
      "i_csr",
      &Bmad::i_csr,
      nb::arg("kick1"),
      nb::arg("i_bin"),
      nb::arg("csr"),
      R"""(Routine to calculate the CSR kick integral (at y = 0)

Parameters
----------
kick1 : CsrKick1Struct

i_bin : int
    Bin index.

csr : CsrStruct
)"""
  );
  m.def(
      "ibs1",
      &Bmad::ibs1,
      nb::arg("lat"),
      nb::arg("ibs_sim_params"),
      nb::arg("rates"),
      nb::arg("i") = nb::none(),
      nb::arg("s") = nb::none(),
      R"""(Calculates IBS growth rates at some location in a lattice.
The IBS rates are betatron growth rates.  That is, they are the rate of
change in sigma_x, sigma_y, and sigma_p.  The emittance growth
rate is twice the betatron growth rate.
1/T_emit = 2/T_betatron.
eg  emit(t) = emit_0 * exp(-2*t/T_betatron) because emit = sigma^2/beta

 Available IBS formulas (ibs_sim_params%formula):
   cimp - Completely Integrated Modified Piwinski
   bjmt - Bjorken-Mtingwa formulation general to bunched beams (time consuming)
   bane - Bane approximation of Bjorken-Mtingwa formulation
   mpzt - Modified Piwinski with Zotter's Integral
   mpxx - Modified Piwinski with a constant Coulomb log.
   kubo - Kubo and Oide's sigma matrix-based

Either i or s, but not both, must be specified.
)"""
  );
  m.def(
      "ibs_blowup1turn",
      &Bmad::ibs_blowup1turn,
      nb::arg("lat"),
      nb::arg("ibs_sim_params"),
      R"""(Updates beam emittances with effect of IBS for
one turn on the lattice.

Parameters
----------
lat : LatStruct
    lattice

ibs_sim_params : IbsSimParamStruct
    Parameters for calculation of IBS rates
)"""
  );
  nb::class_<Bmad::IbsDeltaCalc>(m, "IbsDeltaCalc", "ibs_delta_calc return type")
      .def_ro("delta_sigma_energy", &Bmad::IbsDeltaCalc::delta_sigma_energy)
      .def_ro("delta_emit_a", &Bmad::IbsDeltaCalc::delta_emit_a)
      .def_ro("delta_emit_b", &Bmad::IbsDeltaCalc::delta_emit_b)
      .def("__len__", [](const Bmad::IbsDeltaCalc &) { return 3; })
      .def("__getitem__", [](const Bmad::IbsDeltaCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.delta_sigma_energy);
        if (i == 1)
          return nb::cast(s.delta_emit_a);
        if (i == 2)
          return nb::cast(s.delta_emit_b);
        throw nb::index_error();
      });
  m.def(
      "ibs_delta_calc",
      &Bmad::ibs_delta_calc,
      nb::arg("lat"),
      nb::arg("ix"),
      nb::arg("ibs_sim_params"),
      nb::arg("sigma_mat") = nb::none(),
      R"""(Calculates change in energy spread and emittances due to IBS for a single element.

Parameters
----------
lat : LatStruct
    lattice for tracking

ix : int
    index of element to use: lat.ele(ix)

ibs_sim_params : IbsSimParamStruct
    parameters for calculation of IBS rates.

sigma_mat : 2D array of float (shape: 6,6), optional
    Beam's sigma matrix. Required for 'kubo' method.

Returns
-------
delta_sigma_energy : float, optional
    change in energy spread in eV

delta_emit_a : float, optional
    change in a-mode emittance (geometric)

delta_emit_b : float, optional
    change in b-mode emittance (geometric)
)"""
  );
  m.def(
      "ibs_equib_der",
      &Bmad::ibs_equib_der,
      nb::arg("lat"),
      nb::arg("ibs_sim_params"),
      nb::arg("inmode"),
      nb::arg("granularity"),
      R"""(Computes equilibrium beam sizes by calculating emittance growth rates from IBS growth rates.
Steps beam size through time till equilibrium is reached.

Parameters
----------
lat : LatStruct
    lattice for tracking

ibs_sim_params : IbsSimParamStruct
    parameters for IBS calculation

inmode : NormalModesStruct
    natural beam parameters

granularity : float
    Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter. Set to -1 to calculate
    element-by-element.

Returns
-------
ibsmode : NormalModesStruct
    beam parameters after IBS effects
)"""
  );
  m.def(
      "ibs_equib_rlx",
      &Bmad::ibs_equib_rlx,
      nb::arg("lat"),
      nb::arg("ibs_sim_params"),
      nb::arg("inmode"),
      nb::arg("ratio"),
      nb::arg("initial_blow_up"),
      nb::arg("granularity"),
      R"""(Iterates to equilibrium beam conditions using relaxation method

This method requires that the initial beam size be larger than the equilibrium beam size.
An initial_blow_up of 3 to 5 is a good place to start.

See ibs_rates subroutine for available IBS rate formulas.

Parameters
----------
lat : LatStruct
    lattice for tracking

ibs_sim_params : IbsSimParamStruct
    parameters for IBS calculation

inmode : NormalModesStruct
    natural beam parameters

ratio : float
    Ratio of vert_emit_coupling / vert_emit_total

initial_blow_up : 1D array of float (shape: 3)
    Factor multiplied to all thre bunch dimensions prior to starting iteration.

granularity : float
    Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter.

Returns
-------
ibsmode : NormalModesStruct
    beam parameters after IBS effects
)"""
  );
  m.def(
      "ibs_lifetime",
      &Bmad::ibs_lifetime,
      nb::arg("lat"),
      nb::arg("ibs_sim_params"),
      nb::arg("maxratio"),
      nb::arg("granularity"),
      R"""(This module computes the beam lifetime due to
the diffusion process according to equation 12
from page 129 of The Handbook for Accelerator
Physics and Engineering 2nd edition.

Parameters
----------
lat : LatStruct
    lattice for tracking.

ibs_sim_params : IbsSimParamStruct
    parameters for calculation of IBS rates.

maxratio : IbsMaxratioStruct
    Ax,y,p/sigma_x,y,p where Ax,y,p is the maximum sigma.  Note that this quantity is just the ratio, not the
    ratio squared.  For example, maxratio%Rx = 1.1 says that the maximum acceptable beamsize is 10% larger
    than the beamsize before IBS effects.

granularity : float
    Step size when slicing lattice.  -1 for element-by-element.

Returns
-------
lifetime : IbsLifetimeStruct
    structure returning IBS lifetimes
)"""
  );
  m.def(
      "ibs_matrix_c",
      &Bmad::ibs_matrix_c,
      nb::arg("sigma_mat"),
      nb::arg("tail_cut"),
      nb::arg("tau"),
      nb::arg("energy"),
      nb::arg("n_part"),
      nb::arg("species"),
      R"""(Wrapper for Fortran routine ibs_matrix_c

Parameters
----------
sigma_mat : 2D array of float (shape: 6,6)

tail_cut : bool

tau : float

energy : float

n_part : float

species : int

Returns
-------
ibs_mat : 2D array of float (shape: 6,6)
)"""
  );
  m.def(
      "ibs_rates1turn",
      &Bmad::ibs_rates1turn,
      nb::arg("lat"),
      nb::arg("ibs_sim_params"),
      nb::arg("granularity"),
      R"""(Calculates IBS risetimes for given lat
This is basically a front-end for the various formulas
available in this module of calculating IBS rates.

Parameters
----------
lat : LatStruct
    lattice for tracking.

ibs_sim_params : IbsSimParamStruct
    parameters for IBS calculation.

granularity : float
    slice length.  -1 for element-by-element.

Returns
-------
rates1turn : IbsStruct
    ibs rates for onr turn on the lattice.
)"""
  );
  m.def(
      "igfcoulombfun",
      &Bmad::igfcoulombfun,
      nb::arg("u"),
      nb::arg("v"),
      nb::arg("w"),
      nb::arg("gam"),
      nb::arg("dx"),
      nb::arg("dy"),
      nb::arg("dz"),
      R"""(Wrapper for Fortran routine igfcoulombfun

Parameters
----------
u : float

v : float

w : float

gam : float

dx : float

dy : float

dz : float

Returns
-------
res : float
)"""
  );
  m.def(
      "igfexfun",
      &Bmad::igfexfun,
      nb::arg("u"),
      nb::arg("v"),
      nb::arg("w"),
      nb::arg("gam"),
      nb::arg("dx"),
      nb::arg("dy"),
      nb::arg("dz"),
      R"""(Wrapper for Fortran routine igfexfun

Parameters
----------
u : float

v : float

w : float

gam : float

dx : float

dy : float

dz : float

Returns
-------
res : float
)"""
  );
  m.def(
      "igfeyfun",
      &Bmad::igfeyfun,
      nb::arg("u"),
      nb::arg("v"),
      nb::arg("w"),
      nb::arg("gam"),
      nb::arg("dx"),
      nb::arg("dy"),
      nb::arg("dz"),
      R"""(Wrapper for Fortran routine igfeyfun

Parameters
----------
u : float

v : float

w : float

gam : float

dx : float

dy : float

dz : float

Returns
-------
res : float
)"""
  );
  m.def(
      "igfezfun",
      &Bmad::igfezfun,
      nb::arg("u"),
      nb::arg("v"),
      nb::arg("w"),
      nb::arg("gam"),
      nb::arg("dx"),
      nb::arg("dy"),
      nb::arg("dz"),
      R"""(Wrapper for Fortran routine igfezfun

Parameters
----------
u : float

v : float

w : float

gam : float

dx : float

dy : float

dz : float

Returns
-------
res : float
)"""
  );
  m.def(
      "image_charge_kick_calc",
      &Bmad::image_charge_kick_calc,
      nb::arg("kick1"),
      nb::arg("csr"),
      R"""(Routine to calculate the image charge kick.

Parameters
----------
kick1 : CsrKick1Struct

csr : CsrStruct
)"""
  );
  m.def(
      "init_attribute_name1",
      &Bmad::init_attribute_name1,
      nb::arg("ix_key"),
      nb::arg("ix_attrib"),
      nb::arg("name"),
      nb::arg("attrib_state") = nb::none(),
      nb::arg("override") = nb::none(),
      R"""(Routine to initialize a single name in the element attribute name table.

Parameters
----------
ix_key : int
    Key index.

ix_attrib : int
    Attribute index.

name : str
    Attribute name. Should be uppercase if attrib_state = is_free$. Should contain non-uppercase characters if
    attrib_state = private$.

attrib_state : int, optional
    Class of attribute: does_not_exist$, is_free$, etc. Defaults to is_free$.

override : bool, optional
    Normally this routine throws an error if the [ix_key, ix_attrib] has been set previously. If override =
    True then the set is done and no error is generated.
)"""
  );
  m.def(
      "init_attribute_name_array",
      &Bmad::init_attribute_name_array,
      R"""(Private routine to initialize the attribute name array used by routines
in attribute_mod. Not meant for general use.
)"""
  );
  nb::class_<Bmad::InitBeamDistribution>(
      m,
      "InitBeamDistribution",
      "init_beam_distribution return type"
  )
      .def_ro("beam", &Bmad::InitBeamDistribution::beam)
      .def_ro("err_flag", &Bmad::InitBeamDistribution::err_flag)
      .def_ro("beam_init_set", &Bmad::InitBeamDistribution::beam_init_set)
      .def("__len__", [](const Bmad::InitBeamDistribution &) { return 3; })
      .def("__getitem__", [](const Bmad::InitBeamDistribution &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.beam);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.beam_init_set);
        throw nb::index_error();
      });
  m.def(
      "init_beam_distribution",
      [](EleStruct &ele,
         LatParamStruct &param,
         BeamInitStruct &beam_init,
         NormalModesStruct *modes,
         std::optional<bool> print_p0c_shift_warning,
         std::optional<bool> conserve_momentum) {
        auto fn = static_cast<Bmad::InitBeamDistribution (*)(
            EleStruct &,
            LatParamStruct &,
            BeamInitStruct &,
            optional_ref<NormalModesStruct>,
            std::optional<bool>,
            std::optional<bool>
        )>(&Bmad::init_beam_distribution);
        return fn(
            ele,
            param,
            beam_init,
            ptr_to_opt_ref(modes),
            print_p0c_shift_warning,
            conserve_momentum
        );
      },
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("beam_init"),
      nb::arg("modes") = nb::none(),
      nb::arg("print_p0c_shift_warning") = nb::none(),
      nb::arg("conserve_momentum") = nb::none(),
      R"""(                                                                    print_p0c_shift_warning, conserve_momentum)

Subroutine to initialize a beam of particles.
Initialization uses the downstream parameters of ele.

Note: This routine sets the random number generator according to the settings
in beam_int and at the end resets things to their initial state.

For more information on individual bunch initialization, see the
init_bunch_distribution routine.

Note: The optional "modes" argument generally is used to pass in normal mode parameters as
calculated from the lattice. If present, and if a parameter like beam_init%a_emit are
set negative, then the corresponding parameter in the modes structure is used.
If not present, a warning message is issued and the parameter is set to zero.
This is only used for parameters that cannot be negative.

Parameters
----------
ele : EleStruct
    element to initialize distribution at (downstream end).

param : LatParamStruct
    Lattice parameters

beam_init : BeamInitStruct
    Use "getf beam_init_struct" for more details.

modes : NormalModesStruct, optional
    Normal mode parameters. See above.

print_p0c_shift_warning : bool, optional
    Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

Returns
-------
beam : BeamStruct
    Structure with initialized particles.

err_flag : bool, optional
    Set true if there is an error, false otherwise.

beam_init_set : BeamInitStruct, optional
    Set to input beam_init with components like .a_emit set what is used in constructing the beam (which is
    different from beam_init.a_emit if this is set negative).
)"""
  );
  m.def(
      "init_bmad",
      &Bmad::init_bmad,
      R"""(Wrapper for Fortran routine init_bmad
)"""
  );
  m.def(
      "init_bmad_parser_common",
      [](LatStruct *lat) {
        auto fn = static_cast<void (*)(optional_ref<LatStruct>)>(&Bmad::init_bmad_parser_common);
        return fn(ptr_to_opt_ref(lat));
      },
      nb::arg("lat") = nb::none(),
      R"""(Wrapper for Fortran routine init_bmad_parser_common

Parameters
----------
lat : LatStruct, optional
)"""
  );
  nb::class_<Bmad::InitBunchDistribution>(
      m,
      "InitBunchDistribution",
      "init_bunch_distribution return type"
  )
      .def_ro("bunch", &Bmad::InitBunchDistribution::bunch)
      .def_ro("err_flag", &Bmad::InitBunchDistribution::err_flag)
      .def_ro("beam_init_used", &Bmad::InitBunchDistribution::beam_init_used)
      .def("__len__", [](const Bmad::InitBunchDistribution &) { return 3; })
      .def("__getitem__", [](const Bmad::InitBunchDistribution &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.bunch);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.beam_init_used);
        throw nb::index_error();
      });
  m.def(
      "init_bunch_distribution",
      [](EleStruct &ele,
         LatParamStruct &param,
         BeamInitStruct &beam_init,
         int ix_bunch,
         NormalModesStruct *modes,
         std::optional<bool> print_p0c_shift_warning,
         std::optional<bool> conserve_momentum) {
        auto fn = static_cast<Bmad::InitBunchDistribution (*)(
            EleStruct &,
            LatParamStruct &,
            BeamInitStruct &,
            int,
            optional_ref<NormalModesStruct>,
            std::optional<bool>,
            std::optional<bool>
        )>(&Bmad::init_bunch_distribution);
        return fn(
            ele,
            param,
            beam_init,
            ix_bunch,
            ptr_to_opt_ref(modes),
            print_p0c_shift_warning,
            conserve_momentum
        );
      },
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("beam_init"),
      nb::arg("ix_bunch"),
      nb::arg("modes") = nb::none(),
      nb::arg("print_p0c_shift_warning") = nb::none(),
      nb::arg("conserve_momentum") = nb::none(),
      R"""(                                                                         print_p0c_shift_warning, conserve_momentum)

Subroutine to initialize a distribution of particles of a bunch.
Initialization uses the downstream parameters of ele.

There are four distributions available:
  '', or 'ran_gauss' -- Random gaussian distribution.
  'ellipse'  -- concentric ellipses representing a Gaussian distribution
  'grid'     -- uniform rectangular grid
  'KV'       -- Kapchinsky-Vladimirsky distribution
See the Bmad manual for more information.

The distribution is matched to the Twiss parameters, centroid position, and Energy - z
correlation as specified. Coupling in the element ele is incorporated into the distribution.

Note: Except for the random number seed, the random number generator
parameters used for this routine are set from the beam_init argument.
That is, these parameters are independent of what is used everywhere else.

Note: Make sure: |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z

Note: The optional "modes" argument generally is used to pass in normal mode parameters as
calculated from the lattice. If present, and if a parameter like beam_init%a_emit are
set negative, then the corresponding parameter in the modes structure is used.
If not present, a warning message is issued and the parameter is set to zero.
This is only used for parameters that cannot be negative.

Note: To get good results, It is important to make sure that for
circular rings that beam_init%center is the correct closed orbit.
The closed orbit will shift if, for example, radiation damping is turned on.

Parameters
----------
ele : EleStruct
    element to initialize distribution at (downstream end).

param : LatParamStruct
    Lattice parameters

beam_init : BeamInitStruct
    Use "getf beam_init_struct" for more details.

ix_bunch : int
    Bunch index. 0 = bunch generated at time = 0.

modes : NormalModesStruct, optional
    Normal mode parameters. See above.

print_p0c_shift_warning : bool, optional
    Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

Returns
-------
bunch : BunchStruct
    Structure with initialized particles.

err_flag : bool, optional
    Set True if there is an error. False otherwise.

beam_init_used : BeamInitStruct, optional
    Set to input beam_init with components like .a_emit set what is used in constructing the beam (which can
    be different from beam_init.a_emit if this is set negative). If reading from a file, beam_init_used will
    equal beam_init.
)"""
  );
  m.def(
      "init_complex_taylor_series",
      &Bmad::init_complex_taylor_series,
      nb::arg("complex_taylor"),
      nb::arg("n_term"),
      nb::arg("save") = nb::none(),
      R"""(Subroutine to initialize a Bmad complex_taylor series (6 of these series make
a complex_taylor map). Note: This routine does not zero the structure. The calling
routine is responsible for setting all values.

Parameters
----------
complex_taylor : ComplexTaylorStruct
    Old structure.
    This parameter is an input/output and is modified in-place.
    As an output, complex_taylor: Initalized structure.

n_term : int
    Number of terms to allocate. n_term < 1 => complex_taylor.term pointer will be disassociated.

save : bool, optional
    If True then save any old terms when complex_taylor is resized. Default is False.
)"""
  );
  m.def(
      "init_coord",
      [](CoordStruct &orb,
         std::optional<FixedArray1D<Real, 6>> vec,
         EleStruct *ele,
         std::optional<int> element_end,
         std::optional<int> particle,
         std::optional<int> direction,
         std::optional<double> E_photon,
         std::optional<double> t_offset,
         std::optional<bool> shift_vec6,
         std::optional<FixedArray1D<Real, 3>> spin,
         std::optional<double> s_pos,
         std::optional<bool> random_on) {
        auto fn = static_cast<void (*)(
            CoordStruct &,
            std::optional<FixedArray1D<Real, 6>>,
            optional_ref<EleStruct>,
            std::optional<int>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            std::optional<double>,
            std::optional<bool>,
            std::optional<FixedArray1D<Real, 3>>,
            std::optional<double>,
            std::optional<bool>
        )>(&Bmad::init_coord);
        return fn(
            orb,
            vec,
            ptr_to_opt_ref(ele),
            element_end,
            particle,
            direction,
            E_photon,
            t_offset,
            shift_vec6,
            spin,
            s_pos,
            random_on
        );
      },
      nb::arg("orb"),
      nb::arg("vec") = nb::none(),
      nb::arg("ele") = nb::none(),
      nb::arg("element_end") = nb::none(),
      nb::arg("particle") = nb::none(),
      nb::arg("direction") = nb::none(),
      nb::arg("E_photon") = nb::none(),
      nb::arg("t_offset") = nb::none(),
      nb::arg("shift_vec6") = nb::none(),
      nb::arg("spin") = nb::none(),
      nb::arg("s_pos") = nb::none(),
      nb::arg("random_on") = nb::none(),
      R"""(Routine to initialize a coord_struct.

This routine is an overloaded name for:
  Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
  Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
  Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice),
or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6)
corresponding to the beginning_ele's value of ele%value(p0c_start$).

Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead
of the standard which is to set orb%t from orb%vec(5).

For photons:
  orb%vec(5) is set depending upon where the photon is relative to the element.
  If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle
      except if direction is set.

Parameters
----------
orb : CoordStruct
    Input orbit

vec : 1D array of float (shape: 6), optional
    Coordinate vector. If not present then taken to be zero.

ele : EleStruct, optional
    Particle is initialized to start at element_end of this ele.

element_end : int, optional
    upstream_end$, downstream_end$, inside$, or start_end$. Must be present if ele argument is present.
    start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1. Default is
    upstream_end$. Note: If ele is the beginning element (index zero), the setting of element_end will not
    matter.

particle : int, optional
    Particle type (electron$, etc.). If particle = not_set$ and orb_in is present, use orb_in.species instead.

E_photon : float, optional
    Photon energy if particle is a photon. Ignored otherwise.

t_offset : float, optional
    Offset of the reference time. This is non-zero when there are multiple bunches and the reference time for
    a particular particle is pegged to the time of the center of the bunch.

shift_vec6 : bool, optional
    If present and False, prevent the shift of orb.vec(6).

spin : 1D array of float (shape: 3), optional
    Particle spin. Taken to be zero if not present.

s_pos : float, optional
    Particle s-position. Only relavent if element_end = inside$.

random_on : bool, optional
    Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
    photon coords using a random number generator. If False, the photon coords will be centered within the
    distribution specified in the photon_init ele.
)"""
  );
  m.def(
      "init_coord",
      [](CoordStruct &orb_in,
         EleStruct *ele,
         std::optional<int> element_end,
         std::optional<int> particle,
         std::optional<int> direction,
         std::optional<double> E_photon,
         std::optional<double> t_offset,
         std::optional<bool> shift_vec6,
         std::optional<FixedArray1D<Real, 3>> spin,
         std::optional<double> s_pos,
         std::optional<bool> random_on) {
        auto fn = static_cast<CoordStruct (*)(
            CoordStruct &,
            optional_ref<EleStruct>,
            std::optional<int>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            std::optional<double>,
            std::optional<bool>,
            std::optional<FixedArray1D<Real, 3>>,
            std::optional<double>,
            std::optional<bool>
        )>(&Bmad::init_coord);
        return fn(
            orb_in,
            ptr_to_opt_ref(ele),
            element_end,
            particle,
            direction,
            E_photon,
            t_offset,
            shift_vec6,
            spin,
            s_pos,
            random_on
        );
      },
      nb::arg("orb_in"),
      nb::arg("ele") = nb::none(),
      nb::arg("element_end") = nb::none(),
      nb::arg("particle") = nb::none(),
      nb::arg("direction") = nb::none(),
      nb::arg("E_photon") = nb::none(),
      nb::arg("t_offset") = nb::none(),
      nb::arg("shift_vec6") = nb::none(),
      nb::arg("spin") = nb::none(),
      nb::arg("s_pos") = nb::none(),
      nb::arg("random_on") = nb::none(),
      R"""(Routine to initialize a coord_struct.

This routine is an overloaded name for:
  Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
  Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
  Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice),
or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6)
corresponding to the beginning_ele's value of ele%value(p0c_start$).

Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead
of the standard which is to set orb%t from orb%vec(5).

For photons:
  orb%vec(5) is set depending upon where the photon is relative to the element.
  If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle
      except if direction is set.

Parameters
----------
orb_in : CoordStruct
    Input orbit

ele : EleStruct, optional
    Particle is initialized to start at element_end of this ele.

element_end : int, optional
    upstream_end$, downstream_end$, inside$, or start_end$. Must be present if ele argument is present.
    start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1. Default is
    upstream_end$. Note: If ele is the beginning element (index zero), the setting of element_end will not
    matter.

particle : int, optional
    Particle type (electron$, etc.). If particle = not_set$ and orb_in is present, use orb_in.species instead.

E_photon : float, optional
    Photon energy if particle is a photon. Ignored otherwise.

t_offset : float, optional
    Offset of the reference time. This is non-zero when there are multiple bunches and the reference time for
    a particular particle is pegged to the time of the center of the bunch.

shift_vec6 : bool, optional
    If present and False, prevent the shift of orb.vec(6).

spin : 1D array of float (shape: 3), optional
    Particle spin. Taken to be zero if not present.

s_pos : float, optional
    Particle s-position. Only relavent if element_end = inside$.

random_on : bool, optional
    Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
    photon coords using a random number generator. If False, the photon coords will be centered within the
    distribution specified in the photon_init ele.

Returns
-------
orb_out : CoordStruct
    Initialized coordinate
)"""
  );
  m.def(
      "init_coord",
      [](CoordStruct &orb,
         EleStruct *ele,
         std::optional<int> element_end,
         std::optional<int> particle,
         std::optional<int> direction,
         std::optional<double> E_photon,
         std::optional<double> t_offset,
         std::optional<bool> shift_vec6,
         std::optional<FixedArray1D<Real, 3>> spin) {
        auto fn = static_cast<void (*)(
            CoordStruct &,
            optional_ref<EleStruct>,
            std::optional<int>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            std::optional<double>,
            std::optional<bool>,
            std::optional<FixedArray1D<Real, 3>>
        )>(&Bmad::init_coord);
        return fn(
            orb,
            ptr_to_opt_ref(ele),
            element_end,
            particle,
            direction,
            E_photon,
            t_offset,
            shift_vec6,
            spin
        );
      },
      nb::arg("orb"),
      nb::arg("ele") = nb::none(),
      nb::arg("element_end") = nb::none(),
      nb::arg("particle") = nb::none(),
      nb::arg("direction") = nb::none(),
      nb::arg("E_photon") = nb::none(),
      nb::arg("t_offset") = nb::none(),
      nb::arg("shift_vec6") = nb::none(),
      nb::arg("spin") = nb::none(),
      R"""(Routine to initialize a coord_struct.

This routine is an overloaded name for:
  Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
  Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
  Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice),
or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6)
corresponding to the beginning_ele's value of ele%value(p0c_start$).

Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead
of the standard which is to set orb%t from orb%vec(5).

For photons:
  orb%vec(5) is set depending upon where the photon is relative to the element.
  If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle
      except if direction is set.

Parameters
----------
orb : CoordStruct
    Input orbit

ele : EleStruct, optional
    Particle is initialized to start at element_end of this ele.

element_end : int, optional
    upstream_end$, downstream_end$, inside$, or start_end$. Must be present if ele argument is present.
    start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1. Default is
    upstream_end$. Note: If ele is the beginning element (index zero), the setting of element_end will not
    matter.

particle : int, optional
    Particle type (electron$, etc.). If particle = not_set$ and orb_in is present, use orb_in.species instead.

E_photon : float, optional
    Photon energy if particle is a photon. Ignored otherwise.

t_offset : float, optional
    Offset of the reference time. This is non-zero when there are multiple bunches and the reference time for
    a particular particle is pegged to the time of the center of the bunch.

shift_vec6 : bool, optional
    If present and False, prevent the shift of orb.vec(6).

spin : 1D array of float (shape: 3), optional
    Particle spin. Taken to be zero if not present.
)"""
  );
  m.def(
      "init_custom",
      &Bmad::init_custom,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine init_custom

Parameters
----------
lat : LatStruct
)"""
  );
  m.def(
      "init_ele",
      [](std::optional<int> key,
         std::optional<int> sub_key,
         std::optional<int> ix_ele,
         BranchStruct *branch) {
        auto fn = static_cast<EleStruct (*)(
            std::optional<int>,
            std::optional<int>,
            std::optional<int>,
            optional_ref<BranchStruct>
        )>(&Bmad::init_ele);
        return fn(key, sub_key, ix_ele, ptr_to_opt_ref(branch));
      },
      nb::arg("key") = nb::none(),
      nb::arg("sub_key") = nb::none(),
      nb::arg("ix_ele") = nb::none(),
      nb::arg("branch") = nb::none(),
      R"""(Wrapper for Fortran routine init_ele

Parameters
----------
key : int, optional
    Key to initialize to. EG: quadrupole$, etc.

sub_key : int, optional
    Sub-key to initialize to.

ix_ele : int, optional
    ix_ele index to initalize to. Default = -1.

branch : BranchStruct, optional
    Branch to point ele.branch and ele.ix_branch to. Otherwise ele.branch is nullified and ele.ix_branch = 0.

Returns
-------
ele : EleStruct
    Initialized element.
)"""
  );
  m.def(
      "init_em_taylor_series",
      &Bmad::init_em_taylor_series,
      nb::arg("em_taylor"),
      nb::arg("n_term"),
      nb::arg("save_old") = nb::none(),
      R"""(Subroutine to initialize a Bmad Em_taylor series (6 of these series make
a Em_taylor map). Note: This routine does not zero the structure. The calling
routine is responsible for setting all values.

Parameters
----------
em_taylor : EmTaylorStruct
    Old structure.
    This parameter is an input/output and is modified in-place.
    As an output, em_taylor: Initalized structure.

n_term : int
    Number of terms to allocate. n_term < 0 => em_taylor.term pointer will be disassociated.

save_old : bool, optional
    If True then save any old terms when em_taylor is resized. Default is False.
)"""
  );
  m.def(
      "init_fringe_info",
      [](EleStruct &ele, CoordStruct *orbit, std::optional<int> leng_sign) {
        auto fn = static_cast<
            FringeFieldInfoStruct (*)(EleStruct &, optional_ref<CoordStruct>, std::optional<int>)>(
            &Bmad::init_fringe_info
        );
        return fn(ele, ptr_to_opt_ref(orbit), leng_sign);
      },
      nb::arg("ele"),
      nb::arg("orbit") = nb::none(),
      nb::arg("leng_sign") = nb::none(),
      R"""(Wrapper for Fortran routine init_fringe_info

Parameters
----------
ele : EleStruct
    Lattice element associated with fringe_info.

orbit : CoordStruct, optional
    Particle position. Must be present for a full init. If not full init only fringe_info.has_fringe will be
    set.

leng_sign : int, optional
    Is element length positive (+1) or negative (-1)? Must be present if orbit is present.

Returns
-------
fringe_info : FringeFieldInfoStruct
    Fringe information.
)"""
  );
  m.def(
      "init_lat",
      &Bmad::init_lat,
      nb::arg("n") = nb::none(),
      nb::arg("init_beginning_ele") = nb::none(),
      R"""(Wrapper for Fortran routine init_lat

Parameters
----------
n : int, optional
    Upper bound lat.ele(0:) array is initialized to. Default is 10.

init_beginning_ele : bool, optional
    Init lat.ele(0)? Default is False.

Returns
-------
lat : LatStruct
    Initialized lat.
)"""
  );
  m.def(
      "init_multipole_cache",
      &Bmad::init_multipole_cache,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine init_multipole_cache

Parameters
----------
ele : EleStruct
    Element to init
    This parameter is an input/output and is modified in-place.
    As an output, ele: Initalized element.
)"""
  );
  m.def(
      "init_photon_from_a_photon_init_ele",
      &Bmad::init_photon_from_a_photon_init_ele,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("random_on") = nb::none(),
      R"""(Wrapper for Fortran routine init_photon_from_a_photon_init_ele

Parameters
----------
ele : EleStruct
    patch element.

param : LatParamStruct
    lat_param_struct.

random_on : bool, optional
    : Default is True. If False then use zero for all random numbers needed in the calc.

Returns
-------
orbit : CoordStruct
    Output photon coords.
)"""
  );
  nb::class_<Bmad::InitPhotonIntegProb>(
      m,
      "InitPhotonIntegProb",
      "init_photon_integ_prob return type"
  )
      .def_ro("E_photon", &Bmad::InitPhotonIntegProb::E_photon)
      .def_ro("integ_prob", &Bmad::InitPhotonIntegProb::integ_prob)
      .def("__len__", [](const Bmad::InitPhotonIntegProb &) { return 2; })
      .def("__getitem__", [](const Bmad::InitPhotonIntegProb &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.E_photon);
        if (i == 1)
          return nb::cast(s.integ_prob);
        throw nb::index_error();
      });
  m.def(
      "init_photon_integ_prob",
      &Bmad::init_photon_integ_prob,
      nb::arg("gamma"),
      nb::arg("g"),
      nb::arg("E_min"),
      nb::arg("E_max"),
      nb::arg("vert_angle_min") = nb::none(),
      nb::arg("vert_angle_max") = nb::none(),
      nb::arg("vert_angle_symmetric") = nb::none(),
      nb::arg("energy_integ_prob") = nb::none(),
      R"""(             vert_angle_max, vert_angle_symmetric, energy_integ_prob, E_photon) result (integ_prob)

Routine to calcuate the integrated probability of emitting a photon in a given vertical angle range
and in a given energy range

Parameters
----------
gamma : float
    Gamma factor of charged particle emitting photon.

g : float
    1/rho bending strength.

E_min : float
    Minimum photon energy.

E_max : float
    Maximum photon energy.

vert_angle_min : float, optional
    Lower bound of vertical angle range.

vert_angle_max : float, optional
    Upper bound of vertical angle range.

vert_angle_symmetric : bool, optional
    Use two symmetric ranges [-vert_angle_max, -vert_angle_min] and [vert_angle_min, vert_angle_max] instead
    of just [vert_angle_min, vert_angle_max]?

energy_integ_prob : float, optional
    If present, E_photon will be set to the photon energy such that the integrated probability of generating a
    photon in the given angle and energy range in the interval [E_min, E_photon] is energy_integ_prob. That
    is, energy_integ_prob = 0 => E_photon = E_min and energy_integ_prob = 1 => E_photon = E_max.

Returns
-------
integ_prob : float
    Integrated probablility of emitting a photon in given angle and energy range.

E_photon : float, optional
    See energy_integ_prob. E_photon must be present if energy_integ_prob is.
)"""
  );
  m.def(
      "init_spin_distribution",
      &Bmad::init_spin_distribution,
      nb::arg("beam_init"),
      nb::arg("ele"),
      R"""(Initializes a spin distribution according to beam_init%spin.

Parameters
----------
beam_init : BeamInitStruct
    Initialization parameters

Returns
-------
bunch : BunchStruct
    Bunch of particles. .particle(:).spin
)"""
  );
  m.def(
      "init_surface_segment",
      &Bmad::init_surface_segment,
      nb::arg("phot"),
      nb::arg("ix"),
      nb::arg("iy"),
      R"""(Routine to init the componentes in ele%photon%segmented%pt(ix,iy) for use with segmented surface calculations.

Parameters
----------
phot : PhotonElementStruct
    Surface structure.

ix : int
    index of grid point to init.

iy : int
    index of grid point to init.
)"""
  );
  m.def(
      "init_taylor_series",
      &Bmad::init_taylor_series,
      nb::arg("bmad_taylor"),
      nb::arg("n_term"),
      nb::arg("save_old") = nb::none(),
      R"""(Wrapper for Fortran routine init_taylor_series

Parameters
----------
bmad_taylor : TaylorStruct
    Old structure.
    This parameter is an input/output and is modified in-place.
    As an output, bmad_taylor: Initalized structure.

n_term : int
    Number of terms to allocate. n_term < 0 => bmad_taylor.term pointer will be disassociated.

save_old : bool, optional
    If True then save any old terms and ref orbit when bmad_taylor is resized. If False zero the ref orbit.
    Default is False.
)"""
  );
  m.def(
      "init_wake",
      &Bmad::init_wake,
      nb::arg("n_sr_long"),
      nb::arg("n_sr_trans"),
      nb::arg("n_sr_z"),
      nb::arg("n_lr_mode"),
      nb::arg("always_allocate") = nb::none(),
      R"""(Wrapper for Fortran routine init_wake

Parameters
----------
n_sr_long : int
    Number of terms: wake.sr.long.

n_sr_trans : int
    Number of terms: wake.sr.trans.

n_sr_z : int
    Number of terms: wake.sr.z.

n_lr_mode : int
    Number of terms: wake.lr.mode.

always_allocate : bool, optional
    If present and True then allways allocate wake even if n_lr_mode, etc. are all 0. Default is False.

Returns
-------
wake : WakeStruct, optional
    Initialized structure.
)"""
  );
  m.def(
      "insert_element",
      &Bmad::insert_element,
      nb::arg("lat"),
      nb::arg("insert_ele"),
      nb::arg("insert_index"),
      nb::arg("ix_branch") = nb::none(),
      nb::arg("orbit") = nb::none(),
      R"""(Wrapper for Fortran routine insert_element

Parameters
----------
lat : LatStruct
    lattice that will be modified
    This parameter is an input/output and is modified in-place.
    As an output, lat: lattice with new element inserted

insert_ele : EleStruct
    element to insert into the lat

insert_index : int

ix_branch : int, optional
    : branch index for the insertion. Default = 0.

orbit : 1D array of CoordStruct, optional
    orbit array to enlarge.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Enlarged orbit array.
)"""
  );
  m.def(
      "integrand_base",
      &Bmad::integrand_base,
      nb::arg("t"),
      nb::arg("args"),
      R"""(This vectorized private function is the integrand in equation 31 of Piwinski's paper.

This intetegrand has a sharp exponential decay, and so a change of variables from t to y where t=exp(y)
is applied.  This COV makes the integrand more evenly distributed over the domain of integration,
which makes it easier for qtrap to integrate.

The change of variables is done using integrand_base_cov, which is then integrated
using qtrap.

Parameters
----------
t : float
    Array of reals over which to evaluate the integrand.
)"""
  );
  m.def(
      "integrate_psi",
      &Bmad::integrate_psi,
      nb::arg("bound"),
      nb::arg("p0"),
      nb::arg("args"),
      R"""(Integrate psi(t) from -bound to +bound.  The integration is done in two parts.  First from 0 to -bound, then from
0 to +bound.

Parameters
----------
bound : float
    integration bound

p0 : float
    psi(0).  Boundary condition.

args : 1D array of float (shape: 1:8)
    Parameters and constants of DEQ.  See psi_prime comments for details.

Returns
-------
result : float
    Integral of psi from -bound to +bound.
)"""
  );
  m.def(
      "integrated_mats",
      &Bmad::integrated_mats,
      nb::arg("eles"),
      nb::arg("coos"),
      nb::arg("Lambda"),
      nb::arg("Theta"),
      nb::arg("Iota"),
      nb::arg("mode"),
      R"""(No docstring available.
)"""
  );
  m.def(
      "integration_timer",
      nb::overload_cast<EleStruct &, LatParamStruct &, CoordStruct &, CoordStruct &, double>(
          &Bmad::integration_timer
      ),
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("start"),
      nb::arg("orb_max"),
      nb::arg("tol"),
      R"""(Wrapper for Fortran routine integration_timer_ele

Parameters
----------
ele : EleStruct

param : LatParamStruct

start : CoordStruct

orb_max : CoordStruct

tol : float
)"""
  );
  m.def(
      "integration_timer",
      nb::overload_cast<Fibre &, FixedArray1D<Real, 6>, FixedArray1D<Real, 6>, double>(
          &Bmad::integration_timer
      ),
      nb::arg("a_fibre"),
      nb::arg("orbit"),
      nb::arg("orbit_max"),
      nb::arg("tol_dp"),
      R"""(Wrapper for Fortran routine integration_timer_fibre

Parameters
----------
a_fibre : Fibre

orbit : 1D array of float (shape: 6)

orbit_max : 1D array of float (shape: 6)

tol_dp : float
)"""
  );
  nb::class_<Bmad::InterpolateField>(m, "InterpolateField", "interpolate_field return type")
      .def_ro("E", &Bmad::InterpolateField::E)
      .def_ro("B", &Bmad::InterpolateField::B)
      .def("__len__", [](const Bmad::InterpolateField &) { return 2; })
      .def("__getitem__", [](const Bmad::InterpolateField &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.E);
        if (i == 1)
          return nb::cast(s.B);
        throw nb::index_error();
      });
  m.def(
      "interpolate_field",
      &Bmad::interpolate_field,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("z"),
      nb::arg("mesh3d"),
      R"""(Interpolate field on mesh

Parameters
----------
x : float
    coordinates to interpolate

y : float
    coordinates to interpolate

z : float
    coordinates to interpolate

mesh3d : Mesh3dStruct
    contains efield, bfield

Returns
-------
E : 1D array of float (shape: 3), optional
    interpolated electric field at x, y, z

B : 1D array of float (shape: 3), optional
    interpolated magnetic field at x, y, z
)"""
  );
  m.def(
      "ion_kick",
      &Bmad::ion_kick,
      nb::arg("orbit"),
      nb::arg("r_beam"),
      nb::arg("n_beam_part"),
      nb::arg("a_twiss"),
      nb::arg("b_twiss"),
      nb::arg("sig_ee"),
      R"""(Wrapper for Fortran routine ion_kick

Parameters
----------
orbit : CoordStruct
    Ion position.

r_beam : 1D array of float (shape: 2)
    Beam (x, y) position.

n_beam_part : float
    Number of beam particles.

a_twiss : TwissStruct
    Horizontal like beam twiss parameters.

b_twiss : TwissStruct
    vertical like beam twiss parameters.

sig_ee : float
    Sigma_E/E beam energy spread.

Returns
-------
kick : 1D array of float (shape: 3)
    (x, y, s) kick in m/sec.
)"""
  );
  m.def(
      "is_attribute",
      &Bmad::is_attribute,
      nb::arg("ix_attrib"),
      nb::arg("which"),
      R"""(Routine to determine if an attribute index corresponds to a control variable for overlys/groups.

Parameters
----------
ix_attrib : int
    Attribute index.

which : int
    control_var$, old_control_var$, all_control_var$, multipole$, elec_multipole$

Returns
-------
is_attrib : bool
    True if a control variable
)"""
  );
}
