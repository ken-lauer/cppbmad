#include "pybmad/generated/Bmad_routines_g.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyGammaRef python_gamma_ref(EleProxy& ele, double gamma) {
  Bmad::gamma_ref(ele, gamma);
  auto py_result{PyGammaRef{gamma}};
  return py_result;
}
PyGenGradField python_gen_grad_field(
    RealAlloc1D& deriv,
    GenGrad1Proxy& gg,
    double rho,
    double theta,
    FixedArray1D<Real, 3> field) {
  Bmad::gen_grad_field(deriv, gg, rho, theta, field);
  auto py_result{PyGenGradField{rho, theta}};
  return py_result;
}
PyGetCalledFile python_get_called_file(
    std::string delim,
    std::string call_file,
    bool err) {
  Bmad::get_called_file(delim, call_file, err);
  auto py_result{PyGetCalledFile{delim, call_file, err}};
  return py_result;
}
PyGptFieldGridScaling python_gpt_field_grid_scaling(
    EleProxy& ele,
    int dimensions,
    double field_scale,
    double ref_time) {
  Bmad::gpt_field_grid_scaling(ele, dimensions, field_scale, ref_time);
  auto py_result{PyGptFieldGridScaling{dimensions, field_scale, ref_time}};
  return py_result;
}
PyGptMaxFieldReference python_gpt_max_field_reference(
    GridFieldPt1Proxy& pt0,
    EleProxy& ele,
    double field_value) {
  Bmad::gpt_max_field_reference(pt0, ele, field_value);
  auto py_result{PyGptMaxFieldReference{field_value}};
  return py_result;
}
PyGradientShiftSrWake python_gradient_shift_sr_wake(
    EleProxy& ele,
    LatParamProxy& param,
    double grad_shift) {
  Bmad::gradient_shift_sr_wake(ele, param, grad_shift);
  auto py_result{PyGradientShiftSrWake{grad_shift}};
  return py_result;
}

void init_Bmad_routines_g(py::module& m) {
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
  )""");
  py::class_<
      Bmad::GBendingStrengthFromEmField,
      std::unique_ptr<Bmad::GBendingStrengthFromEmField>>(
      m,
      "GBendingStrengthFromEmField",
      "g_bending_strength_from_em_field return type")
      .def_readonly("g", &Bmad::GBendingStrengthFromEmField::g)
      .def_readonly("dg", &Bmad::GBendingStrengthFromEmField::dg)
      .def(
          "__len__", [](const Bmad::GBendingStrengthFromEmField&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::GBendingStrengthFromEmField& s, int i) -> py::object {
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
  )""");
  m.def(
      "g_integrals_calc",
      &Bmad::g_integrals_calc,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Lattice to integrate through.
  )""");
  py::class_<PyGammaRef, std::unique_ptr<PyGammaRef>>(
      m, "GammaRef", "gamma_ref return type")
      .def_readonly("gamma", &PyGammaRef::gamma)
      .def("__len__", [](const PyGammaRef&) { return 1; })
      .def("__getitem__", [](const PyGammaRef& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.gamma);
        throw py::index_error();
      });
  m.def(
      "gamma_ref",
      &python_gamma_ref,
      py::arg("ele"),
      py::arg("gamma"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to evaluate at.
  gamma : 
  )""");
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
  )""");
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
  )""");
  py::class_<PyGenGradField, std::unique_ptr<PyGenGradField>>(
      m, "GenGradField", "gen_grad_field return type")
      .def_readonly("rho", &PyGenGradField::rho)
      .def_readonly("theta", &PyGenGradField::theta)
      .def("__len__", [](const PyGenGradField&) { return 2; })
      .def("__getitem__", [](const PyGenGradField& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.rho);
        if (i == 1)
          return py::cast(s.theta);
        throw py::index_error();
      });
  m.def(
      "gen_grad_field",
      &python_gen_grad_field,
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
  )""");
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
  )""");
  py::class_<PyGetCalledFile, std::unique_ptr<PyGetCalledFile>>(
      m, "GetCalledFile", "get_called_file return type")
      .def_readonly("delim", &PyGetCalledFile::delim)
      .def_readonly("call_file", &PyGetCalledFile::call_file)
      .def_readonly("err", &PyGetCalledFile::err)
      .def("__len__", [](const PyGetCalledFile&) { return 3; })
      .def("__getitem__", [](const PyGetCalledFile& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.call_file);
        if (i == 2)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "get_called_file",
      &python_get_called_file,
      py::arg("delim"),
      py::arg("call_file"),
      py::arg("err"),
      R"""(Parameters
  ----------
  delim : 
  call_file : 
  err : 
  )""");
  py::class_<
      Bmad::GetEmitFromSigmaMat,
      std::unique_ptr<Bmad::GetEmitFromSigmaMat>>(
      m, "GetEmitFromSigmaMat", "get_emit_from_sigma_mat return type")
      .def_readonly("normal", &Bmad::GetEmitFromSigmaMat::normal)
      .def_readonly("err_flag", &Bmad::GetEmitFromSigmaMat::err_flag)
      .def("__len__", [](const Bmad::GetEmitFromSigmaMat&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::GetEmitFromSigmaMat& s, int i) -> py::object {
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
  )""");
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
  )""");
  py::class_<Bmad::GetSlaveList, std::unique_ptr<Bmad::GetSlaveList>>(
      m, "GetSlaveList", "get_slave_list return type")
      .def_readonly("slaves", &Bmad::GetSlaveList::slaves)
      .def_readonly("n_slave", &Bmad::GetSlaveList::n_slave)
      .def("__len__", [](const Bmad::GetSlaveList&) { return 2; })
      .def("__getitem__", [](const Bmad::GetSlaveList& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.slaves);
        if (i == 1)
          return py::cast(s.n_slave);
        throw py::index_error();
      });
  m.def(
      "get_slave_list",
      &Bmad::get_slave_list,
      py::arg("lord"),
      R"""(Parameters
  ----------
  lord : EleStruct
      The lord element.
  slaves : ElePointerStruct
      : Array of slaves.
  n_slave : int
      Number of slaves.
  )""");
  py::class_<PyGptFieldGridScaling, std::unique_ptr<PyGptFieldGridScaling>>(
      m, "GptFieldGridScaling", "gpt_field_grid_scaling return type")
      .def_readonly("dimensions", &PyGptFieldGridScaling::dimensions)
      .def_readonly("field_scale", &PyGptFieldGridScaling::field_scale)
      .def_readonly("ref_time", &PyGptFieldGridScaling::ref_time)
      .def("__len__", [](const PyGptFieldGridScaling&) { return 3; })
      .def(
          "__getitem__",
          [](const PyGptFieldGridScaling& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.dimensions);
            if (i == 1)
              return py::cast(s.field_scale);
            if (i == 2)
              return py::cast(s.ref_time);
            throw py::index_error();
          });
  m.def(
      "gpt_field_grid_scaling",
      &python_gpt_field_grid_scaling,
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
  )""");
  py::class_<PyGptMaxFieldReference, std::unique_ptr<PyGptMaxFieldReference>>(
      m, "GptMaxFieldReference", "gpt_max_field_reference return type")
      .def_readonly("field_value", &PyGptMaxFieldReference::field_value)
      .def("__len__", [](const PyGptMaxFieldReference&) { return 1; })
      .def(
          "__getitem__",
          [](const PyGptMaxFieldReference& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.field_value);
            throw py::index_error();
          });
  m.def(
      "gpt_max_field_reference",
      &python_gpt_max_field_reference,
      py::arg("pt0"),
      py::arg("ele"),
      py::arg("field_value"),
      R"""(Parameters
  ----------
  pt0 : 
  ele : 
  field_value : 
  )""");
  py::class_<
      Bmad::GptToParticleBunch,
      std::unique_ptr<Bmad::GptToParticleBunch>>(
      m, "GptToParticleBunch", "gpt_to_particle_bunch return type")
      .def_readonly("bunch", &Bmad::GptToParticleBunch::bunch)
      .def_readonly("err_flag", &Bmad::GptToParticleBunch::err_flag)
      .def("__len__", [](const Bmad::GptToParticleBunch&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::GptToParticleBunch& s, int i) -> py::object {
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
  )""");
  py::class_<PyGradientShiftSrWake, std::unique_ptr<PyGradientShiftSrWake>>(
      m, "GradientShiftSrWake", "gradient_shift_sr_wake return type")
      .def_readonly("grad_shift", &PyGradientShiftSrWake::grad_shift)
      .def("__len__", [](const PyGradientShiftSrWake&) { return 1; })
      .def(
          "__getitem__",
          [](const PyGradientShiftSrWake& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.grad_shift);
            throw py::index_error();
          });
  m.def(
      "gradient_shift_sr_wake",
      &python_gradient_shift_sr_wake,
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
  )""");
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
  )""");
}
