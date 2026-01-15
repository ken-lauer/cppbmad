#include "pybmad/generated/Bmad_routines_g.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_g(py::module &m) {
  m.def(
      "g_bend_from_em_field",
      &Bmad::g_bend_from_em_field,
      py::arg("b"),
      py::arg("e"),
      py::arg("orbit"),
      R"""(Function g_bend_from_em_field (B, E, orbit) result (g_bend)

Routine to calculate the bending strength (1/bending_radius) for a given particle for a given field.
This will include the dipole bending field of an sbend.

Parameters
----------
B : float
    Magnetic field.
E : float
    Electric field
orbit : CoordStruct
    particle orbit

Returns
-------
g_bend : float
    bending strength vector.
)"""
  );
  py::class_<Bmad::GBendingStrengthFromEmField, std::unique_ptr<Bmad::GBendingStrengthFromEmField>>(
      m,
      "GBendingStrengthFromEmField",
      "g_bending_strength_from_em_field return type"
  )
      .def_readonly("g", &Bmad::GBendingStrengthFromEmField::g)
      .def_readonly("dg", &Bmad::GBendingStrengthFromEmField::dg)
      .def("__len__", [](const Bmad::GBendingStrengthFromEmField &) { return 2; })
      .def("__getitem__", [](const Bmad::GBendingStrengthFromEmField &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.g);
        if (i == 1)
          return py::cast(s.dg);
        throw py::index_error();
      });
  m.def(
      "g_bending_strength_from_em_field",
      &Bmad::g_bending_strength_from_em_field,
      py::arg("ele"),
      py::arg("param"),
      py::arg("s_rel"),
      py::arg("orbit"),
      py::arg("local_ref_frame"),
      R"""(Parameters
----------
ele : EleStruct
    Element being tracked thorugh.
param : LatParamStruct
    Lattice parameters.
s_rel : float
    Distance from the start of the element to the particle.
orbit : CoordStruct
    Particle position in lab (not element) frame.
local_ref_frame : 
    Logical, If True then take the input coordinates and output g as being with respect to the frame of
    referene of the element (ignore misalignments).
g : float
    g = (g_x, g_y, g_s) bending strength vector (|g| = 1/bend_radius).
dg : float
    dg(:)/dr gradient. Takes into account dg_x/dx in a bend due to curvilinear coords.
)"""
  );
  m.def(
      "g_integrals_calc",
      &Bmad::g_integrals_calc,
      py::arg("lat"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice to integrate through.
)"""
  );
  m.def(
      "gamma_ref",
      &Bmad::gamma_ref,
      py::arg("ele"),
      py::arg("gamma"),
      R"""(Parameters
----------
ele : EleStruct
    Element to evaluate at.
gamma : 
)"""
  );
  m.def(
      "gen_grad1_to_em_taylor",
      &Bmad::gen_grad1_to_em_taylor,
      py::arg("ele"),
      py::arg("gen_grad"),
      py::arg("iz"),
      R"""(Parameters
----------
ele : unknown
    Element containing the map.
gen_grad : GenGradMapStruct
    Gen_grad map.
iz : int
    z-plane index to evaluate.
em_taylor : EmTaylorStruct
    Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
)"""
  );
  m.def(
      "gen_grad_at_s_to_em_taylor",
      &Bmad::gen_grad_at_s_to_em_taylor,
      py::arg("ele"),
      py::arg("gen_grad"),
      py::arg("s_pos"),
      R"""(Parameters
----------
ele : unknown
    Element containing the map.
gen_grad : GenGradMapStruct
    Gen_grad map.
s_pos : float
    Position to evaluate em_taylor at.
em_taylor : EmTaylorStruct
    Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
)"""
  );
  m.def(
      "gen_grad_field",
      &Bmad::gen_grad_field,
      py::arg("deriv"),
      py::arg("gg"),
      py::arg("rho"),
      py::arg("theta"),
      py::arg("field"),
      R"""(Parameters
----------
deriv : 
gg : 
rho : 
theta : 
field : 
)"""
  );
  m.def(
      "get_bl_from_fwhm",
      &Bmad::get_bl_from_fwhm,
      py::arg("bound"),
      py::arg("args"),
      R"""(Subroutine get_bl_from_fwhm(bound,args,sigma)

Calculate bunch length as fwhm * c_light / TwoRtTwoLnTwo.
Where fwhm is full width at half max of solution to dpsi/dt.

Parameters
----------
bound : float
    -bound and +bound are lower and upper integration bound.
args : float
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
      py::arg("delim"),
      py::arg("call_file"),
      py::arg("err"),
      R"""(Parameters
----------
delim : 
call_file : 
err : 
)"""
  );
  py::class_<Bmad::GetEmitFromSigmaMat, std::unique_ptr<Bmad::GetEmitFromSigmaMat>>(
      m,
      "GetEmitFromSigmaMat",
      "get_emit_from_sigma_mat return type"
  )
      .def_readonly("normal", &Bmad::GetEmitFromSigmaMat::normal)
      .def_readonly("err_flag", &Bmad::GetEmitFromSigmaMat::err_flag)
      .def("__len__", [](const Bmad::GetEmitFromSigmaMat &) { return 2; })
      .def("__getitem__", [](const Bmad::GetEmitFromSigmaMat &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.normal);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "get_emit_from_sigma_mat",
      &Bmad::get_emit_from_sigma_mat,
      py::arg("sigma_mat"),
      py::arg("Nmat") = py::none(),
      R"""(Subroutine get_emit_from_sigma_mat(sigma_mat, normal, Nmat, err_flag)

Given a beam envelop sigma matrix sigma_mat, this returns the 3 normal mode
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
sigma_mat : float
    beam envelop sigma matrix
Nmat : float, optional
    If present, then the emittanced will be ordered such that the eigensystem most closely resembles Nmat.

Returns
-------
normal : float
    normal mode emittances
err_flag : bool
    Set to true if something went wrong.  Otherwise set to false.
)"""
  );
  m.def(
      "get_next_word",
      &Bmad::get_next_word,
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim_list"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("upper_case_word") = py::none(),
      py::arg("call_check") = py::none(),
      py::arg("err_flag") = py::none(),
      R"""(Subroutine get_next_word (word, ix_word, delim_list, delim, delim_found, upper_case_word, call_check, err_flag)

Subroutine to get the next word from the input stream.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
word : unknown
    Word returned
delim_list : unknown
    List of valid delimiters
upper_case_word : bool, optional
    if True then convert word to upper case. Default is True.
call_check : bool, optional
    If present and True then check for 'call::<filename>' construct. Default is False. Output
ix_word : int
    length of word argument
delim : unknown
    Actual delimiter found
delim_found : bool
    Set true if a delimiter found. A delimiter may not be found if the end of the line is reached first.
err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "get_slave_list",
      &Bmad::get_slave_list,
      py::arg("lord"),
      py::arg("slaves"),
      R"""(Parameters
----------
lord : EleStruct
    The lord element.
slaves : ElePointerStruct
    : Array of slaves.
n_slave : int
    Number of slaves.
)"""
  );
  m.def(
      "gpt_field_grid_scaling",
      &Bmad::gpt_field_grid_scaling,
      py::arg("ele"),
      py::arg("dimensions"),
      py::arg("field_scale"),
      py::arg("ref_time"),
      R"""(Parameters
----------
ele : 
dimensions : 
field_scale : 
ref_time : 
)"""
  );
  m.def(
      "gpt_max_field_reference",
      &Bmad::gpt_max_field_reference,
      py::arg("pt0"),
      py::arg("ele"),
      py::arg("field_value"),
      R"""(Parameters
----------
pt0 : 
ele : 
field_value : 
)"""
  );
  py::class_<Bmad::GptToParticleBunch, std::unique_ptr<Bmad::GptToParticleBunch>>(
      m,
      "GptToParticleBunch",
      "gpt_to_particle_bunch return type"
  )
      .def_readonly("bunch", &Bmad::GptToParticleBunch::bunch)
      .def_readonly("err_flag", &Bmad::GptToParticleBunch::err_flag)
      .def("__len__", [](const Bmad::GptToParticleBunch &) { return 2; })
      .def("__getitem__", [](const Bmad::GptToParticleBunch &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.bunch);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "gpt_to_particle_bunch",
      &Bmad::gpt_to_particle_bunch,
      py::arg("gpt_file"),
      py::arg("ele"),
      R"""(Subroutine gpt_to_particle_bunch (gpt_file, ele, bunch, err_flag)

Routine to initialize a bunch of particles from a GPT screen file.


Parameters
----------
gpt_file : unknown
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
      py::arg("ele"),
      py::arg("param"),
      py::arg("grad_shift"),
      R"""(Parameters
----------
ele : EleStruct
    Lcavity element.
param : LatParamStruct
    Lattice parameters .n_part        -- Number of particles in a bunch .particle      -- Type of particle
grad_shift : 
)"""
  );
  m.def(
      "grid_field_interpolate",
      &Bmad::grid_field_interpolate,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("grid"),
      py::arg("err_flag"),
      py::arg("x1"),
      py::arg("x2") = py::none(),
      py::arg("x3") = py::none(),
      py::arg("allow_s_out_of_bounds") = py::none(),
      py::arg("print_err") = py::none(),
      R"""(Subroutine grid_field_interpolate (ele, orbit, grid, field, err_flag, x1, x2, x3, &
                                                             allow_s_out_of_bounds, print_err)

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

Returns
-------
field : GridFieldPtStruct
    Interpolated field (complex)
)"""
  );
}
