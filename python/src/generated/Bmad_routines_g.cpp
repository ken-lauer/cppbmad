#include "pybmad/generated/Bmad_routines_g.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_g(nb::module_ &m) {
  m.def(
      "g_bend_from_em_field",
      &Bmad::g_bend_from_em_field,
      nb::arg("b"),
      nb::arg("e"),
      nb::arg("orbit"),
      R"""(Routine to calculate the bending strength (1/bending_radius) for a given particle for a given field.
This will include the dipole bending field of an sbend.

Parameters
----------
B : 1D array of float (shape: 3)
    Magnetic field.

E : 1D array of float (shape: 3)
    Electric field

orbit : CoordStruct
    particle orbit

Returns
-------
g_bend : 1D array of float (shape: 3)
    bending strength vector.
)"""
  );
  nb::class_<Bmad::GBendingStrengthFromEmField>(
      m,
      "GBendingStrengthFromEmField",
      "g_bending_strength_from_em_field return type"
  )
      .def_ro("g", &Bmad::GBendingStrengthFromEmField::g)
      .def_ro("dg", &Bmad::GBendingStrengthFromEmField::dg)
      .def("__len__", [](const Bmad::GBendingStrengthFromEmField &) { return 2; })
      .def("__getitem__", [](const Bmad::GBendingStrengthFromEmField &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.g);
        if (i == 1)
          return nb::cast(s.dg);
        throw nb::index_error();
      });
  m.def(
      "g_bending_strength_from_em_field",
      &Bmad::g_bending_strength_from_em_field,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("s_rel"),
      nb::arg("orbit"),
      nb::arg("local_ref_frame"),
      R"""(Wrapper for Fortran routine g_bending_strength_from_em_field

Parameters
----------
ele : EleStruct
    Element being tracked thorugh.

param : LatParamStruct
    Lattice parameters.

s_rel : float
    Distance from the start of the element to the particle.

orbit : CoordStruct
    Particle position in lab (not element) frame.

local_ref_frame : bool
    Logical, If True then take the input coordinates and output g as being with respect to the frame of
    referene of the element (ignore misalignments).

Returns
-------
g : 1D array of float (shape: 3)
    g = (g_x, g_y, g_s) bending strength vector (|g| = 1/bend_radius).

dg : 2D array of float (shape: 3,3), optional
    dg(:)/dr gradient. Takes into account dg_x/dx in a bend due to curvilinear coords.
)"""
  );
  m.def(
      "g_integrals_calc",
      &Bmad::g_integrals_calc,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine g_integrals_calc

Parameters
----------
lat : LatStruct
    Lattice to integrate through.
)"""
  );
  m.def(
      "gamma_ref",
      &Bmad::gamma_ref,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine gamma_ref

Parameters
----------
ele : EleStruct
    Element to evaluate at.

Returns
-------
gamma : float
    Relativistic gamma factor Energy/mass*c^2.
)"""
  );
  m.def(
      "gen_grad1_to_gg_taylor",
      &Bmad::gen_grad1_to_gg_taylor,
      nb::arg("ele"),
      nb::arg("gen_grad"),
      nb::arg("iz"),
      R"""(Wrapper for Fortran routine gen_grad1_to_gg_taylor

Parameters
----------
ele : EleStruct
    Element containing the map.

gen_grad : GenGradMapStruct
    Gen_grad map.

iz : int
    z-plane index to evaluate.

Returns
-------
gg_taylor : 1D array of GgTaylorStruct (shape: 3)
    Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
)"""
  );
  m.def(
      "gen_grad_at_s_to_gg_taylor",
      &Bmad::gen_grad_at_s_to_gg_taylor,
      nb::arg("ele"),
      nb::arg("gen_grad"),
      nb::arg("s_pos"),
      R"""(Wrapper for Fortran routine gen_grad_at_s_to_gg_taylor

Parameters
----------
ele : EleStruct
    Element containing the map.

gen_grad : GenGradMapStruct
    Gen_grad map.

s_pos : float
    Position to evaluate gg_taylor at.

Returns
-------
gg_taylor : 1D array of GgTaylorStruct (shape: 3)
    Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
)"""
  );
  m.def(
      "gen_grad_field",
      &Bmad::gen_grad_field,
      nb::arg("deriv"),
      nb::arg("gg"),
      nb::arg("rho"),
      nb::arg("theta"),
      R"""(Wrapper for Fortran routine gen_grad_field

Parameters
----------
deriv : 1D array of float

gg : GenGrad1Struct

rho : float

theta : float

Returns
-------
field : 1D array of float (shape: 3)
)"""
  );
  m.def(
      "get_bl_from_fwhm",
      &Bmad::get_bl_from_fwhm,
      nb::arg("bound"),
      nb::arg("args"),
      R"""(Calculate bunch length as fwhm * c_light / TwoRtTwoLnTwo.
Where fwhm is full width at half max of solution to dpsi/dt.

Parameters
----------
bound : float
    -bound and +bound are lower and upper integration bound.

args : 1D array of float (shape: 1:8)
    Parameters and constants of dpsi/dt.  See comments of psi_prime for details.

Returns
-------
sigma : float
    Bunch length
)"""
  );
  m.def(
      "get_called_file",
      &Bmad::get_called_file,
      nb::arg("delim"),
      nb::arg("call_file"),
      nb::arg("err"),
      R"""(Wrapper for Fortran routine get_called_file

Parameters
----------
delim : str

call_file : str

err : bool
)"""
  );
  nb::class_<Bmad::GetEmitFromSigmaMat>(
      m,
      "GetEmitFromSigmaMat",
      "get_emit_from_sigma_mat return type"
  )
      .def_ro("normal", &Bmad::GetEmitFromSigmaMat::normal)
      .def_ro("err_flag", &Bmad::GetEmitFromSigmaMat::err_flag)
      .def("__len__", [](const Bmad::GetEmitFromSigmaMat &) { return 2; })
      .def("__getitem__", [](const Bmad::GetEmitFromSigmaMat &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.normal);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "get_emit_from_sigma_mat",
      &Bmad::get_emit_from_sigma_mat,
      nb::arg("sigma_mat"),
      nb::arg("Nmat") = nb::none(),
      R"""(Given a beam envelop sigma matrix sigma_mat, this returns the 3 normal mode
emittances.

The normal mode emittance of the sigma matrix are the eigenvalues of
sigma_mat . S

If Nmat is present, then the modes are ordered such that the eigensystem most
closely resembles Nmat.  If Nmat is not present, then the modes are ordered
according to which plane they dominate.

    / 0  1  0  0  0  0 \
    |-1  0  0  0  0  0 |
S = | 0  0  0  1  0  0 |
    | 0  0 -1  0  0  0 |
    | 0  0  0  0  0  1 |
    \ 0  0  0  0 -1  0 /

Parameters
----------
sigma_mat : 2D array of float (shape: 6,6)
    beam envelop sigma matrix

Nmat : 2D array of float (shape: 6,6), optional
    If present, then the emittanced will be ordered such that the eigensystem most closely resembles Nmat.

Returns
-------
normal : 1D array of float (shape: 3)
    normal mode emittances

err_flag : bool
    Set to true if something went wrong.  Otherwise set to false.
)"""
  );
  m.def(
      "get_list_of_names",
      &Bmad::get_list_of_names,
      nb::arg("ele"),
      nb::arg("err_str"),
      nb::arg("name_list"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
      R"""(This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  nb::class_<Bmad::GetNextWord>(m, "GetNextWord", "get_next_word return type")
      .def_ro("ix_word", &Bmad::GetNextWord::ix_word)
      .def_ro("delim", &Bmad::GetNextWord::delim)
      .def_ro("delim_found", &Bmad::GetNextWord::delim_found)
      .def_ro("err_flag", &Bmad::GetNextWord::err_flag)
      .def("__len__", [](const Bmad::GetNextWord &) { return 4; })
      .def("__getitem__", [](const Bmad::GetNextWord &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.ix_word);
        if (i == 1)
          return nb::cast(s.delim);
        if (i == 2)
          return nb::cast(s.delim_found);
        if (i == 3)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "get_next_word",
      &Bmad::get_next_word,
      nb::arg("word"),
      nb::arg("delim_list"),
      nb::arg("upper_case_word") = nb::none(),
      nb::arg("call_check") = nb::none(),
      R"""(Subroutine to get the next word from the input stream.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
word : str
    Word returned

delim_list : str
    List of valid delimiters

upper_case_word : bool, optional
    if True then convert word to upper case. Default is True.

call_check : bool, optional
    If present and True then check for 'call::<filename>' construct. Default is False.

Returns
-------
ix_word : int
    length of word argument

delim : str
    Actual delimiter found

delim_found : bool
    Set true if a delimiter found. A delimiter may not be found if the end of the line is reached first.

err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "get_overlay_group_names",
      [](EleStruct &ele,
         LatStruct &lat,
         ParserEleStruct &pele,
         std::string delim,
         bool delim_found,
         bool is_control_var_list,
         bool err_flag,
         CharacterAlloc1D *names_out) {
        auto fn = static_cast<
            void (*)(EleStruct &, LatStruct &, ParserEleStruct &, std::string, bool, bool, bool, optional_ref<CharacterAlloc1D>)>(
            &Bmad::get_overlay_group_names
        );
        return fn(
            ele,
            lat,
            pele,
            delim,
            delim_found,
            is_control_var_list,
            err_flag,
            ptr_to_opt_ref(names_out)
        );
      },
      nb::arg("ele"),
      nb::arg("lat"),
      nb::arg("pele"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("is_control_var_list"),
      nb::arg("err_flag"),
      nb::arg("names_out") = nb::none(),
      R"""(This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
is_control_var_list : bool
    If True then parsing "var = {...}" list. If False then parsing "group/overlay/girder = {...}" list.
)"""
  );
  m.def(
      "get_sequence_args",
      &Bmad::get_sequence_args,
      nb::arg("seq_name"),
      nb::arg("arg_list"),
      nb::arg("delim"),
      nb::arg("err_flag"),
      R"""(Subroutine to get the argument list for a replacement_line or a list.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  nb::class_<Bmad::GetSlaveList>(m, "GetSlaveList", "get_slave_list return type")
      .def_ro("slaves", &Bmad::GetSlaveList::slaves)
      .def_ro("n_slave", &Bmad::GetSlaveList::n_slave)
      .def("__len__", [](const Bmad::GetSlaveList &) { return 2; })
      .def("__getitem__", [](const Bmad::GetSlaveList &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.slaves);
        if (i == 1)
          return nb::cast(s.n_slave);
        throw nb::index_error();
      });
  m.def(
      "get_slave_list",
      &Bmad::get_slave_list,
      nb::arg("lord"),
      R"""(Wrapper for Fortran routine get_slave_list

Parameters
----------
lord : EleStruct
    The lord element.

Returns
-------
slaves : 1D array of ElePointerStruct
    : Array of slaves.

n_slave : int
    Number of slaves.
)"""
  );
  m.def(
      "get_switch",
      &Bmad::get_switch,
      nb::arg("name"),
      nb::arg("name_list"),
      nb::arg("switch_"),
      nb::arg("err"),
      nb::arg("ele"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      R"""(Wrapper for Fortran routine get_switch

Parameters
----------
name : str

name_list : 1D array of str

err : bool

ele : EleStruct

delim : str

delim_found : bool
)"""
  );
  m.def(
      "gg_taylor_equal_gg_taylor",
      &Bmad::gg_taylor_equal_gg_taylor,
      nb::arg("gg_taylor1"),
      nb::arg("gg_taylor2"),
      R"""(Wrapper for Fortran routine gg_taylor_equal_gg_taylor

Parameters
----------
gg_taylor1 : GgTaylorStruct

gg_taylor2 : GgTaylorStruct
)"""
  );
  m.def(
      "gg_taylors_equal_gg_taylors",
      &Bmad::gg_taylors_equal_gg_taylors,
      nb::arg("gg_taylor1"),
      nb::arg("gg_taylor2"),
      R"""(Wrapper for Fortran routine gg_taylors_equal_gg_taylors

Parameters
----------
gg_taylor1 : 1D array of GgTaylorStruct

gg_taylor2 : 1D array of GgTaylorStruct
)"""
  );
  m.def(
      "gpt_field_grid_scaling",
      &Bmad::gpt_field_grid_scaling,
      nb::arg("ele"),
      nb::arg("dimensions"),
      nb::arg("field_scale"),
      nb::arg("ref_time"),
      R"""(Wrapper for Fortran routine gpt_field_grid_scaling

Parameters
----------
ele : EleStruct

dimensions : int

field_scale : float

ref_time : float
)"""
  );
  m.def(
      "gpt_max_field_reference",
      &Bmad::gpt_max_field_reference,
      nb::arg("pt0"),
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine gpt_max_field_reference

Parameters
----------
pt0 : GridFieldPt1Struct

ele : EleStruct

Returns
-------
field_value : float
)"""
  );
  nb::class_<Bmad::GptToParticleBunch>(m, "GptToParticleBunch", "gpt_to_particle_bunch return type")
      .def_ro("bunch", &Bmad::GptToParticleBunch::bunch)
      .def_ro("err_flag", &Bmad::GptToParticleBunch::err_flag)
      .def("__len__", [](const Bmad::GptToParticleBunch &) { return 2; })
      .def("__getitem__", [](const Bmad::GptToParticleBunch &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.bunch);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "gpt_to_particle_bunch",
      &Bmad::gpt_to_particle_bunch,
      nb::arg("gpt_file"),
      nb::arg("ele"),
      R"""(Routine to initialize a bunch of particles from a GPT screen file.

Parameters
----------
gpt_file : str
    Name of GPT data file.

ele : EleStruct
    Lattice element whose downstream end coincident with the GPT screen.

Returns
-------
bunch : BunchStruct
    Particle bunch

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "gradient_shift_sr_wake",
      &Bmad::gradient_shift_sr_wake,
      nb::arg("ele"),
      nb::arg("param"),
      R"""(Wrapper for Fortran routine gradient_shift_sr_wake

Parameters
----------
ele : EleStruct
    Lcavity element.

param : LatParamStruct
    Lattice parameters

Returns
-------
grad_shift : float
    Shift in gradient
)"""
  );
  m.def(
      "grid_field_interpolate",
      &Bmad::grid_field_interpolate,
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("grid"),
      nb::arg("err_flag"),
      nb::arg("x1"),
      nb::arg("x2") = nb::none(),
      nb::arg("x3") = nb::none(),
      nb::arg("allow_s_out_of_bounds") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(                                                             allow_s_out_of_bounds, print_err)

Subroutine to interpolate the E and B fields on a rectilinear grid.

Parameters
----------
ele : EleStruct
    Element containing the grid.

orbit : CoordStruct
    Used for constructing an error message if the particle is out of bounds.

grid : GridFieldStruct
    Grid to interpolate.

err_flag : bool
    Set to true if there is an error. False otherwise.

x1 : float
    dimension 1 interpolation point.

x2 : float, optional
    dimension 2 interpolation point.

x3 : float, optional
    dimension 3 interpolation point.

allow_s_out_of_bounds : bool, optional
    allow s-coordinate grossly out of bounds to return zero field without an error. This is used when the
    field of one element overlaps the field of another. Default is False.

print_err : bool, optional
    print an error message if the particle is out of bounds? Default is True.
)"""
  );
}
