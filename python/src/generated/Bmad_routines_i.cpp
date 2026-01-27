#include "pybmad/generated/Bmad_routines_i.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyInitAttributeName1 python_init_attribute_name1(
    bool is_ok,
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state = std::nullopt,
    std::optional<bool> override = std::nullopt
) {
  Bmad::init_attribute_name1(is_ok, ix_key, ix_attrib, name, attrib_state, override);
  auto py_result{PyInitAttributeName1{is_ok}};
  return py_result;
}

void init_Bmad_routines_i(py::module &m) {
  m.def(
      "ibs_matrix_c",
      &Bmad::ibs_matrix_c,
      py::arg("sigma_mat"),
      py::arg("tail_cut"),
      py::arg("tau"),
      py::arg("energy"),
      py::arg("n_part"),
      py::arg("species"),
      py::arg("ibs_mat"),
      R"""(Wrapper for Fortran routine ibs_matrix_c

Parameters
----------
sigma_mat : 2D array of float (shape: 6,6)

tail_cut : bool

tau : float

energy : float

n_part : float

species : int

ibs_mat : 2D array of float (shape: 6,6)
)"""
  );
  m.def(
      "igfcoulombfun",
      &Bmad::igfcoulombfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
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

res : float
)"""
  );
  m.def(
      "igfexfun",
      &Bmad::igfexfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
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

res : float
)"""
  );
  m.def(
      "igfeyfun",
      &Bmad::igfeyfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
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

res : float
)"""
  );
  m.def(
      "igfezfun",
      &Bmad::igfezfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
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

res : float
)"""
  );
  py::class_<PyInitAttributeName1, std::unique_ptr<PyInitAttributeName1>>(
      m,
      "InitAttributeName1",
      "init_attribute_name1 return type"
  )
      .def_readonly("is_ok", &PyInitAttributeName1::is_ok)
      .def("__len__", [](const PyInitAttributeName1 &) { return 1; })
      .def("__getitem__", [](const PyInitAttributeName1 &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "init_attribute_name1",
      &python_init_attribute_name1,
      py::arg("is_ok"),
      py::arg("ix_key"),
      py::arg("ix_attrib"),
      py::arg("name"),
      py::arg("attrib_state") = py::none(),
      py::arg("override") = py::none(),
      R"""(subroutine init_attribute_name1 (is_ok, ix_key, ix_attrib, name, attrib_state, override)

Routine to initialize a single name in the element attribute name table.

Parameters
----------
is_ok : bool
    Initial setting.
    This parameter is an input/output and is modified in-place.
    As an output, is_ok: Set False if there is a problem. Otherwise untouched.

ix_key : int
    Key index.

ix_attrib : int
    Attribute index.

name : character
    Attribute name. Should be uppercase if attrib_state = is_free$. Should contain non-uppercase characters if
    attrib_state = private$.

attrib_state : int, optional
    Class of attribute: does_not_exist$, is_free$, etc. Defaults to is_free$.

override : bool, optional
    Normally this routine throws an error if the [ix_key, ix_attrib] has been set previously. If override =
    True then the set is done and no error is generated.

Returns
-------
is_ok : bool
    Initial setting.
    This parameter is an input/output and is modified in-place.
    As an output, is_ok: Set False if there is a problem. Otherwise untouched.
)"""
  );
  m.def(
      "init_attribute_name_array",
      &Bmad::init_attribute_name_array,
      R"""(Subroutine init_attribute_name_array ()

Private routine to initialize the attribute name array used by routines
in attribute_mod. Not meant for general use.
)"""
  );
  py::class_<Bmad::InitBeamDistribution, std::unique_ptr<Bmad::InitBeamDistribution>>(
      m,
      "InitBeamDistribution",
      "init_beam_distribution return type"
  )
      .def_readonly("beam", &Bmad::InitBeamDistribution::beam)
      .def_readonly("err_flag", &Bmad::InitBeamDistribution::err_flag)
      .def_readonly("beam_init_set", &Bmad::InitBeamDistribution::beam_init_set)
      .def("__len__", [](const Bmad::InitBeamDistribution &) { return 3; })
      .def("__getitem__", [](const Bmad::InitBeamDistribution &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.beam);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.beam_init_set);
        throw py::index_error();
      });
  m.def(
      "init_beam_distribution",
      &Bmad::init_beam_distribution,
      py::arg("ele"),
      py::arg("param"),
      py::arg("beam_init"),
      py::arg("modes") = py::none(),
      py::arg("print_p0c_shift_warning") = py::none(),
      py::arg("conserve_momentum") = py::none(),
      R"""(Subroutine init_beam_distribution (ele, param, beam_init, beam, err_flag, modes, beam_init_set,
                                                                    print_p0c_shift_warning, conserve_momentum)

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
      &Bmad::init_bmad_parser_common,
      py::arg("lat") = py::none(),
      R"""(Wrapper for Fortran routine init_bmad_parser_common

Parameters
----------
lat : LatStruct, optional
)"""
  );
  py::class_<Bmad::InitBunchDistribution, std::unique_ptr<Bmad::InitBunchDistribution>>(
      m,
      "InitBunchDistribution",
      "init_bunch_distribution return type"
  )
      .def_readonly("bunch", &Bmad::InitBunchDistribution::bunch)
      .def_readonly("err_flag", &Bmad::InitBunchDistribution::err_flag)
      .def_readonly("beam_init_used", &Bmad::InitBunchDistribution::beam_init_used)
      .def("__len__", [](const Bmad::InitBunchDistribution &) { return 3; })
      .def("__getitem__", [](const Bmad::InitBunchDistribution &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.bunch);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.beam_init_used);
        throw py::index_error();
      });
  m.def(
      "init_bunch_distribution",
      &Bmad::init_bunch_distribution,
      py::arg("ele"),
      py::arg("param"),
      py::arg("beam_init"),
      py::arg("ix_bunch"),
      py::arg("modes") = py::none(),
      py::arg("print_p0c_shift_warning") = py::none(),
      py::arg("conserve_momentum") = py::none(),
      R"""(Subroutine init_bunch_distribution (ele, param, beam_init, ix_bunch, bunch, err_flag, modes, beam_init_used,
                                                                         print_p0c_shift_warning, conserve_momentum)

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
      py::arg("complex_taylor"),
      py::arg("n_term"),
      py::arg("save") = py::none(),
      R"""(Subroutine init_complex_taylor_series (complex_taylor, n_term, save)

Subroutine to initialize a Bmad complex_taylor series (6 of these series make
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
      py::overload_cast<
          CoordStruct &,
          FixedArray1D<Real, 6>,
          optional_ref<EleStruct>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<double>,
          std::optional<double>,
          std::optional<bool>,
          std::optional<FixedArray1D<Real, 3>>,
          std::optional<double>,
          std::optional<bool>>(&Bmad::init_coord),
      py::arg("orb"),
      py::arg("vec"),
      py::arg("ele") = py::none(),
      py::arg("element_end") = py::none(),
      py::arg("particle") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("E_photon") = py::none(),
      py::arg("t_offset") = py::none(),
      py::arg("shift_vec6") = py::none(),
      py::arg("spin") = py::none(),
      py::arg("s_pos") = py::none(),
      py::arg("random_on") = py::none(),
      R"""(Subroutine init_coord (...)

Routine to initialize a coord_struct.

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

vec : 1D array of float (shape: 6)
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

direction : int, optional
    +1 -> moving downstream +s direciton, -1 -> moving upstream. 0 -> Ignore. Default is to not change
    orb.direction except for photons which get set according to orb.vec(6).

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
      py::overload_cast<
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
          std::optional<bool>>(&Bmad::init_coord),
      py::arg("orb_in"),
      py::arg("ele") = py::none(),
      py::arg("element_end") = py::none(),
      py::arg("particle") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("E_photon") = py::none(),
      py::arg("t_offset") = py::none(),
      py::arg("shift_vec6") = py::none(),
      py::arg("spin") = py::none(),
      py::arg("s_pos") = py::none(),
      py::arg("random_on") = py::none(),
      R"""(Subroutine init_coord (...)

Routine to initialize a coord_struct.

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

direction : int, optional
    +1 -> moving downstream +s direciton, -1 -> moving upstream. 0 -> Ignore. Default is to not change
    orb.direction except for photons which get set according to orb.vec(6).

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
      py::overload_cast<
          CoordStruct &,
          optional_ref<EleStruct>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<double>,
          std::optional<double>,
          std::optional<bool>,
          std::optional<FixedArray1D<Real, 3>>>(&Bmad::init_coord),
      py::arg("orb"),
      py::arg("ele") = py::none(),
      py::arg("element_end") = py::none(),
      py::arg("particle") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("E_photon") = py::none(),
      py::arg("t_offset") = py::none(),
      py::arg("shift_vec6") = py::none(),
      py::arg("spin") = py::none(),
      R"""(Subroutine init_coord (...)

Routine to initialize a coord_struct.

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

direction : int, optional
    +1 -> moving downstream +s direciton, -1 -> moving upstream. 0 -> Ignore. Default is to not change
    orb.direction except for photons which get set according to orb.vec(6).

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
      py::arg("lat"),
      R"""(Wrapper for Fortran routine init_custom

Parameters
----------
lat : LatStruct
)"""
  );
  m.def(
      "init_ele",
      &Bmad::init_ele,
      py::arg("key") = py::none(),
      py::arg("sub_key") = py::none(),
      py::arg("ix_ele") = py::none(),
      py::arg("branch") = py::none(),
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
      py::arg("em_taylor"),
      py::arg("n_term"),
      py::arg("save_old") = py::none(),
      R"""(Subroutine init_em_taylor_series (em_taylor, n_term, save_old)

Subroutine to initialize a Bmad Em_taylor series (6 of these series make
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
      "init_lat",
      &Bmad::init_lat,
      py::arg("n") = py::none(),
      py::arg("init_beginning_ele") = py::none(),
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
      py::arg("ele"),
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
      py::arg("ele"),
      py::arg("param"),
      py::arg("random_on") = py::none(),
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
  py::class_<Bmad::InitPhotonIntegProb, std::unique_ptr<Bmad::InitPhotonIntegProb>>(
      m,
      "InitPhotonIntegProb",
      "init_photon_integ_prob return type"
  )
      .def_readonly("E_photon", &Bmad::InitPhotonIntegProb::E_photon)
      .def_readonly("integ_prob", &Bmad::InitPhotonIntegProb::integ_prob)
      .def("__len__", [](const Bmad::InitPhotonIntegProb &) { return 2; })
      .def("__getitem__", [](const Bmad::InitPhotonIntegProb &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.E_photon);
        if (i == 1)
          return py::cast(s.integ_prob);
        throw py::index_error();
      });
  m.def(
      "init_photon_integ_prob",
      &Bmad::init_photon_integ_prob,
      py::arg("gamma"),
      py::arg("g"),
      py::arg("E_min"),
      py::arg("E_max"),
      py::arg("vert_angle_min") = py::none(),
      py::arg("vert_angle_max") = py::none(),
      py::arg("vert_angle_symmetric") = py::none(),
      py::arg("energy_integ_prob") = py::none(),
      R"""(Function init_photon_integ_prob(gamma, g, E_min, E_max, vert_angle_min,
             vert_angle_max, vert_angle_symmetric, energy_integ_prob, E_photon) result (integ_prob)

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
E_photon : float, optional
    See energy_integ_prob. E_photon must be present if energy_integ_prob is.

integ_prob : float
    Integrated probablility of emitting a photon in given angle and energy range.
)"""
  );
  m.def(
      "init_spin_distribution",
      &Bmad::init_spin_distribution,
      py::arg("beam_init"),
      py::arg("ele"),
      R"""(Subroutine init_spin_distribution (beam_init, bunch, ele)

Initializes a spin distribution according to beam_init%spin.

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
      py::arg("phot"),
      py::arg("ix"),
      py::arg("iy"),
      R"""(Subroutine init_surface_segment (phot, ix, iy)

Routine to init the componentes in ele%photon%segmented%pt(ix,iy) for use with segmented surface calculations.

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
      py::arg("bmad_taylor"),
      py::arg("n_term"),
      py::arg("save_old") = py::none(),
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
      py::arg("n_sr_long"),
      py::arg("n_sr_trans"),
      py::arg("n_sr_z"),
      py::arg("n_lr_mode"),
      py::arg("always_allocate") = py::none(),
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
      py::arg("lat"),
      py::arg("insert_ele"),
      py::arg("ix_ele"),
      py::arg("ix_branch") = py::none(),
      py::arg("orbit") = py::none(),
      R"""(Wrapper for Fortran routine insert_element

Parameters
----------
lat : LatStruct
    lattice that will be modified
    This parameter is an input/output and is modified in-place.
    As an output, lat: lattice with new element inserted

insert_ele : EleStruct
    element to insert into the lat

ix_ele : int
    branch.ele(:) index where the new element is inserted.

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
      py::arg("t"),
      py::arg("args"),
      py::arg("func_retval__"),
      R"""(Function integrand_base(t)

This vectorized private function is the integrand in equation 31 of Piwinski's paper.

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
      py::arg("bound"),
      py::arg("p0"),
      py::arg("args"),
      R"""(Subroutine integrate_psi(bound,p0,args,result)

Integrate psi(t) from -bound to +bound.  The integration is done in two parts.  First from 0 to -bound, then from
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
      py::arg("eles"),
      py::arg("coos"),
      py::arg("Lambda"),
      py::arg("Theta"),
      py::arg("Iota"),
      py::arg("mode"),
      R"""(subroutine integrated_mats(eles,coos,Lambda,Theta,Iota,mode)
)"""
  );
  m.def(
      "integration_timer",
      py::overload_cast<EleStruct &, LatParamStruct &, CoordStruct &, CoordStruct &, double>(
          &Bmad::integration_timer
      ),
      py::arg("ele"),
      py::arg("param"),
      py::arg("start"),
      py::arg("orb_max"),
      py::arg("tol"),
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
      py::overload_cast<Fibre &, FixedArray1D<Real, 6>, FixedArray1D<Real, 6>, double>(
          &Bmad::integration_timer
      ),
      py::arg("a_fibre"),
      py::arg("orbit"),
      py::arg("orbit_max"),
      py::arg("tol_dp"),
      R"""(Wrapper for Fortran routine integration_timer_fibre

Parameters
----------
a_fibre : Fibre

orbit : 1D array of float (shape: 6)

orbit_max : 1D array of float (shape: 6)

tol_dp : float
)"""
  );
  m.def(
      "ion_kick",
      &Bmad::ion_kick,
      py::arg("orbit"),
      py::arg("r_beam"),
      py::arg("n_beam_part"),
      py::arg("a_twiss"),
      py::arg("b_twiss"),
      py::arg("sig_ee"),
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
      py::arg("ix_attrib"),
      py::arg("which"),
      R"""(Function is_attribute (ix_attrib, which) result (is_attrib)

Routine to determine if an attribute index corresponds to a control variable for overlys/groups.

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
