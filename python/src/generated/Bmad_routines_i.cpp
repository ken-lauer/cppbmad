#include "pybmad/generated/Bmad_routines_i.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyIbsMatrixC python_ibs_matrix_c(
    FixedArray2D<Real, 6, 6> sigma_mat,
    bool tail_cut,
    double tau,
    double energy,
    double n_part,
    int species,
    FixedArray2D<Real, 6, 6> ibs_mat) {
  Bmad::ibs_matrix_c(
      sigma_mat, tail_cut, tau, energy, n_part, species, ibs_mat);
  auto py_result{PyIbsMatrixC{tail_cut, tau, energy, n_part, species}};
  return py_result;
}
PyIgfcoulombfun python_igfcoulombfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res) {
  Bmad::igfcoulombfun(u, v, w, gam, dx, dy, dz, res);
  auto py_result{PyIgfcoulombfun{u, v, w, gam, dx, dy, dz, res}};
  return py_result;
}
PyIgfexfun python_igfexfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res) {
  Bmad::igfexfun(u, v, w, gam, dx, dy, dz, res);
  auto py_result{PyIgfexfun{u, v, w, gam, dx, dy, dz, res}};
  return py_result;
}
PyIgfeyfun python_igfeyfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res) {
  Bmad::igfeyfun(u, v, w, gam, dx, dy, dz, res);
  auto py_result{PyIgfeyfun{u, v, w, gam, dx, dy, dz, res}};
  return py_result;
}
PyIgfezfun python_igfezfun(
    double u,
    double v,
    double w,
    double gam,
    double dx,
    double dy,
    double dz,
    double res) {
  Bmad::igfezfun(u, v, w, gam, dx, dy, dz, res);
  auto py_result{PyIgfezfun{u, v, w, gam, dx, dy, dz, res}};
  return py_result;
}
PyInitAttributeName1 python_init_attribute_name1(
    bool is_ok,
    int ix_key,
    int ix_attrib,
    std::string name,
    std::optional<int> attrib_state = std::nullopt,
    std::optional<bool> override = std::nullopt) {
  Bmad::init_attribute_name1(
      is_ok, ix_key, ix_attrib, name, attrib_state, override);
  auto py_result{PyInitAttributeName1{is_ok}};
  return py_result;
}
PyInitBeamDistribution python_init_beam_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    optional_ref<NormalModesProxy> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    std::optional<bool> conserve_momentum = std::nullopt) {
  auto _result = Bmad::init_beam_distribution(
      ele,
      param,
      beam_init,
      modes,
      print_p0c_shift_warning,
      make_opt_ref(conserve_momentum));
  auto py_result{PyInitBeamDistribution{_result, conserve_momentum}};
  return py_result;
}
PyInitBunchDistribution python_init_bunch_distribution(
    EleProxy& ele,
    LatParamProxy& param,
    BeamInitProxy& beam_init,
    int ix_bunch,
    optional_ref<NormalModesProxy> modes = std::nullopt,
    std::optional<bool> print_p0c_shift_warning = std::nullopt,
    std::optional<bool> conserve_momentum = std::nullopt) {
  auto _result = Bmad::init_bunch_distribution(
      ele,
      param,
      beam_init,
      ix_bunch,
      modes,
      print_p0c_shift_warning,
      make_opt_ref(conserve_momentum));
  auto py_result{PyInitBunchDistribution{_result, conserve_momentum}};
  return py_result;
}
PyInitSurfaceSegment python_init_surface_segment(
    PhotonElementProxy& phot,
    int ix,
    int iy) {
  Bmad::init_surface_segment(phot, ix, iy);
  auto py_result{PyInitSurfaceSegment{ix, iy}};
  return py_result;
}
PyIntegrandBase python_integrand_base(
    double t,
    RealAlloc1D& args,
    double func_retval__) {
  Bmad::integrand_base(t, args, func_retval__);
  auto py_result{PyIntegrandBase{func_retval__}};
  return py_result;
}
PyIntegrationTimerEle python_integration_timer_ele(
    EleProxy& ele,
    LatParamProxy& param,
    CoordProxy& start,
    CoordProxy& orb_max,
    double tol) {
  Bmad::integration_timer(ele, param, start, orb_max, tol);
  auto py_result{PyIntegrationTimerEle{tol}};
  return py_result;
}

void init_Bmad_routines_i(py::module& m) {
  py::class_<PyIbsMatrixC, std::unique_ptr<PyIbsMatrixC>>(
      m, "IbsMatrixC", "ibs_matrix_c return type")
      .def_readonly("tail_cut", &PyIbsMatrixC::tail_cut)
      .def_readonly("tau", &PyIbsMatrixC::tau)
      .def_readonly("energy", &PyIbsMatrixC::energy)
      .def_readonly("n_part", &PyIbsMatrixC::n_part)
      .def_readonly("species", &PyIbsMatrixC::species)
      .def("__len__", [](const PyIbsMatrixC&) { return 5; })
      .def("__getitem__", [](const PyIbsMatrixC& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.tail_cut);
        if (i == 1)
          return py::cast(s.tau);
        if (i == 2)
          return py::cast(s.energy);
        if (i == 3)
          return py::cast(s.n_part);
        if (i == 4)
          return py::cast(s.species);
        throw py::index_error();
      });
  m.def(
      "ibs_matrix_c",
      &python_ibs_matrix_c,
      py::arg("sigma_mat"),
      py::arg("tail_cut"),
      py::arg("tau"),
      py::arg("energy"),
      py::arg("n_part"),
      py::arg("species"),
      py::arg("ibs_mat"),
      R"""(Parameters
  ----------
  sigma_mat : 
  tail_cut : 
  tau : 
  energy : 
  n_part : 
  species : 
  ibs_mat : 
  )""");
  py::class_<PyIgfcoulombfun, std::unique_ptr<PyIgfcoulombfun>>(
      m, "Igfcoulombfun", "igfcoulombfun return type")
      .def_readonly("u", &PyIgfcoulombfun::u)
      .def_readonly("v", &PyIgfcoulombfun::v)
      .def_readonly("w", &PyIgfcoulombfun::w)
      .def_readonly("gam", &PyIgfcoulombfun::gam)
      .def_readonly("dx", &PyIgfcoulombfun::dx)
      .def_readonly("dy", &PyIgfcoulombfun::dy)
      .def_readonly("dz", &PyIgfcoulombfun::dz)
      .def_readonly("res", &PyIgfcoulombfun::res)
      .def("__len__", [](const PyIgfcoulombfun&) { return 8; })
      .def("__getitem__", [](const PyIgfcoulombfun& s, int i) -> py::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return py::cast(s.u);
        if (i == 1)
          return py::cast(s.v);
        if (i == 2)
          return py::cast(s.w);
        if (i == 3)
          return py::cast(s.gam);
        if (i == 4)
          return py::cast(s.dx);
        if (i == 5)
          return py::cast(s.dy);
        if (i == 6)
          return py::cast(s.dz);
        if (i == 7)
          return py::cast(s.res);
        throw py::index_error();
      });
  m.def(
      "igfcoulombfun",
      &python_igfcoulombfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
      R"""(Parameters
  ----------
  u : 
  v : 
  w : 
  gam : 
  dx : 
  dy : 
  dz : 
  res : 
  )""");
  py::class_<PyIgfexfun, std::unique_ptr<PyIgfexfun>>(
      m, "Igfexfun", "igfexfun return type")
      .def_readonly("u", &PyIgfexfun::u)
      .def_readonly("v", &PyIgfexfun::v)
      .def_readonly("w", &PyIgfexfun::w)
      .def_readonly("gam", &PyIgfexfun::gam)
      .def_readonly("dx", &PyIgfexfun::dx)
      .def_readonly("dy", &PyIgfexfun::dy)
      .def_readonly("dz", &PyIgfexfun::dz)
      .def_readonly("res", &PyIgfexfun::res)
      .def("__len__", [](const PyIgfexfun&) { return 8; })
      .def("__getitem__", [](const PyIgfexfun& s, int i) -> py::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return py::cast(s.u);
        if (i == 1)
          return py::cast(s.v);
        if (i == 2)
          return py::cast(s.w);
        if (i == 3)
          return py::cast(s.gam);
        if (i == 4)
          return py::cast(s.dx);
        if (i == 5)
          return py::cast(s.dy);
        if (i == 6)
          return py::cast(s.dz);
        if (i == 7)
          return py::cast(s.res);
        throw py::index_error();
      });
  m.def(
      "igfexfun",
      &python_igfexfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
      R"""(Parameters
  ----------
  u : 
  v : 
  w : 
  gam : 
  dx : 
  dy : 
  dz : 
  res : 
  )""");
  py::class_<PyIgfeyfun, std::unique_ptr<PyIgfeyfun>>(
      m, "Igfeyfun", "igfeyfun return type")
      .def_readonly("u", &PyIgfeyfun::u)
      .def_readonly("v", &PyIgfeyfun::v)
      .def_readonly("w", &PyIgfeyfun::w)
      .def_readonly("gam", &PyIgfeyfun::gam)
      .def_readonly("dx", &PyIgfeyfun::dx)
      .def_readonly("dy", &PyIgfeyfun::dy)
      .def_readonly("dz", &PyIgfeyfun::dz)
      .def_readonly("res", &PyIgfeyfun::res)
      .def("__len__", [](const PyIgfeyfun&) { return 8; })
      .def("__getitem__", [](const PyIgfeyfun& s, int i) -> py::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return py::cast(s.u);
        if (i == 1)
          return py::cast(s.v);
        if (i == 2)
          return py::cast(s.w);
        if (i == 3)
          return py::cast(s.gam);
        if (i == 4)
          return py::cast(s.dx);
        if (i == 5)
          return py::cast(s.dy);
        if (i == 6)
          return py::cast(s.dz);
        if (i == 7)
          return py::cast(s.res);
        throw py::index_error();
      });
  m.def(
      "igfeyfun",
      &python_igfeyfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
      R"""(Parameters
  ----------
  u : 
  v : 
  w : 
  gam : 
  dx : 
  dy : 
  dz : 
  res : 
  )""");
  py::class_<PyIgfezfun, std::unique_ptr<PyIgfezfun>>(
      m, "Igfezfun", "igfezfun return type")
      .def_readonly("u", &PyIgfezfun::u)
      .def_readonly("v", &PyIgfezfun::v)
      .def_readonly("w", &PyIgfezfun::w)
      .def_readonly("gam", &PyIgfezfun::gam)
      .def_readonly("dx", &PyIgfezfun::dx)
      .def_readonly("dy", &PyIgfezfun::dy)
      .def_readonly("dz", &PyIgfezfun::dz)
      .def_readonly("res", &PyIgfezfun::res)
      .def("__len__", [](const PyIgfezfun&) { return 8; })
      .def("__getitem__", [](const PyIgfezfun& s, int i) -> py::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return py::cast(s.u);
        if (i == 1)
          return py::cast(s.v);
        if (i == 2)
          return py::cast(s.w);
        if (i == 3)
          return py::cast(s.gam);
        if (i == 4)
          return py::cast(s.dx);
        if (i == 5)
          return py::cast(s.dy);
        if (i == 6)
          return py::cast(s.dz);
        if (i == 7)
          return py::cast(s.res);
        throw py::index_error();
      });
  m.def(
      "igfezfun",
      &python_igfezfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("dx"),
      py::arg("dy"),
      py::arg("dz"),
      py::arg("res"),
      R"""(Parameters
  ----------
  u : 
  v : 
  w : 
  gam : 
  dx : 
  dy : 
  dz : 
  res : 
  )""");
  py::class_<PyInitAttributeName1, std::unique_ptr<PyInitAttributeName1>>(
      m, "InitAttributeName1", "init_attribute_name1 return type")
      .def_readonly("is_ok", &PyInitAttributeName1::is_ok)
      .def("__len__", [](const PyInitAttributeName1&) { return 1; })
      .def(
          "__getitem__",
          [](const PyInitAttributeName1& s, int i) -> py::object {
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
      This parameter is an input/output and is modified in-place. As an output: Set False if there is a problem.
      Otherwise untouched.
  ix_key : int
      Key index.
  ix_attrib : int
      Attribute index.
  name : unknown
      Attribute name. Should be uppercase if attrib_state = is_free$. Should contain non-uppercase characters if
      attrib_state = private$.
  attrib_state : int, optional
      Class of attribute: does_not_exist$, is_free$, etc. Defaults to is_free$.
  override : bool, optional
      Normally this routine throws an error if the [ix_key, ix_attrib] has been set previously. If override =
      True then the set is done and no error is generated.
  )""");
  m.def(
      "init_attribute_name_array",
      &Bmad::init_attribute_name_array,
      R"""(Subroutine init_attribute_name_array ()

  Private routine to initialize the attribute name array used by routines
  in attribute_mod. Not meant for general use.

  )""");
  py::class_<PyInitBeamDistribution, std::unique_ptr<PyInitBeamDistribution>>(
      m, "InitBeamDistribution", "init_beam_distribution return type")
      .def_readonly("beam", &PyInitBeamDistribution::beam)
      .def_readonly("err_flag", &PyInitBeamDistribution::err_flag)
      .def_readonly("beam_init_set", &PyInitBeamDistribution::beam_init_set)
      .def_readonly(
          "conserve_momentum", &PyInitBeamDistribution::conserve_momentum)
      .def("__len__", [](const PyInitBeamDistribution&) { return 4; })
      .def(
          "__getitem__",
          [](const PyInitBeamDistribution& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.beam);
            if (i == 1)
              return py::cast(s.err_flag);
            if (i == 2)
              return py::cast(s.beam_init_set);
            if (i == 3)
              return py::cast(s.conserve_momentum);
            throw py::index_error();
          });
  m.def(
      "init_beam_distribution",
      &python_init_beam_distribution,
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
      Lattice parameters .particle              -- Type of particle.
  beam_init : BeamInitStruct
      Use "getf beam_init_struct" for more details.
  modes : NormalModesStruct, optional
      Normal mode parameters. See above.
  print_p0c_shift_warning : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
  shift_momentum : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

  Returns
  -------
  beam : BeamStruct
      Structure with initialized particles.
  err_flag : bool
      Set true if there is an error, false otherwise.
  beam_init_set : BeamInitStruct
      Set to input beam_init with components like .a_emit set what is used in constructing the beam (which is
      different from beam_init.a_emit if this is set negative).
  )""");
  m.def("init_bmad", &Bmad::init_bmad, R"""()""");
  m.def(
      "init_bmad_parser_common",
      &Bmad::init_bmad_parser_common,
      py::arg("lat") = py::none(),
      R"""(Parameters
  ----------
  lat : 
  )""");
  py::class_<PyInitBunchDistribution, std::unique_ptr<PyInitBunchDistribution>>(
      m, "InitBunchDistribution", "init_bunch_distribution return type")
      .def_readonly("bunch", &PyInitBunchDistribution::bunch)
      .def_readonly("err_flag", &PyInitBunchDistribution::err_flag)
      .def_readonly("beam_init_used", &PyInitBunchDistribution::beam_init_used)
      .def_readonly(
          "conserve_momentum", &PyInitBunchDistribution::conserve_momentum)
      .def("__len__", [](const PyInitBunchDistribution&) { return 4; })
      .def(
          "__getitem__",
          [](const PyInitBunchDistribution& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.bunch);
            if (i == 1)
              return py::cast(s.err_flag);
            if (i == 2)
              return py::cast(s.beam_init_used);
            if (i == 3)
              return py::cast(s.conserve_momentum);
            throw py::index_error();
          });
  m.def(
      "init_bunch_distribution",
      &python_init_bunch_distribution,
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
  shift_momentum : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

  Returns
  -------
  bunch : BunchStruct
      Structure with initialized particles.
  err_flag : bool
      Set True if there is an error. False otherwise.
  beam_init_used : BeamInitStruct
      Set to input beam_init with components like .a_emit set what is used in constructing the beam (which can
      be different from beam_init.a_emit if this is set negative). If reading from a file, beam_init_used will
      equal beam_init.
  )""");
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
      This parameter is an input/output and is modified in-place. As an output: Initalized structure.
  n_term : int
      Number of terms to allocate. n_term < 1 => complex_taylor.term pointer will be disassociated.
  save : bool, optional
      If True then save any old terms when complex_taylor is resized. Default is False.
  )""");
  m.def(
      "init_coord",
      py::overload_cast<
          CoordProxy&,
          FixedArray1D<Real, 6>,
          optional_ref<EleProxy>,
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
  orb_in : CoordStruct
      Input orbit.
  vec : float
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
  spin : float, optional
      Particle spin. Taken to be zero if not present.
  s_pos : float, optional
      Particle s-position. Only relavent if element_end = inside$.
  random_on : bool, optional
      Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
      photon coords using a random number generator. If False, the photon coords will be centered within the
      distribution specified in the photon_init ele.

  Returns
  -------
  orb : coord_struct
      Input orbit

  Notes
  -----
  Overloaded versions:
  )""");
  m.def(
      "init_coord",
      py::overload_cast<
          CoordProxy&,
          optional_ref<EleProxy>,
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
  orb_in : coord_struct
      Input orbit
  vec : float
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
  spin : float, optional
      Particle spin. Taken to be zero if not present.
  s_pos : float, optional
      Particle s-position. Only relavent if element_end = inside$.
  random_on : bool, optional
      Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
      photon coords using a random number generator. If False, the photon coords will be centered within the
      distribution specified in the photon_init ele.

  Returns
  -------
  orb : CoordStruct
      Initialized coordinate. Note: For photons, orb.vec(6) is computed as sqrt(1 - vec(2)^2 - vec(4)^2) if
      needed.
  orb_out : coord_struct
      Initialized coordinate

  Notes
  -----
  Overloaded versions:
  )""");
  m.def(
      "init_coord",
      py::overload_cast<
          CoordProxy&,
          optional_ref<EleProxy>,
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
  orb_in : CoordStruct
      Input orbit.
  vec : float
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
  spin : float, optional
      Particle spin. Taken to be zero if not present.
  s_pos : float, optional
      Particle s-position. Only relavent if element_end = inside$.
  random_on : bool, optional
      Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
      photon coords using a random number generator. If False, the photon coords will be centered within the
      distribution specified in the photon_init ele.

  Returns
  -------
  orb : coord_struct
      Input orbit

  Notes
  -----
  Overloaded versions:
  )""");
  m.def(
      "init_custom",
      &Bmad::init_custom,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : 
  )""");
  m.def(
      "init_ele",
      &Bmad::init_ele,
      py::arg("key") = py::none(),
      py::arg("sub_key") = py::none(),
      py::arg("ix_ele") = py::none(),
      py::arg("branch") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Initialized element.
  key : int, optional
      Key to initialize to. EG: quadrupole$, etc.
  sub_key : int, optional
      Sub-key to initialize to.
  ix_ele : int, optional
      ix_ele index to initalize to. Default = -1.
  branch : BranchStruct, optional
      Branch to point ele.branch and ele.ix_branch to. Otherwise ele.branch is nullified and ele.ix_branch = 0.
  )""");
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
      This parameter is an input/output and is modified in-place. As an output: Initalized structure.
  n_term : int
      Number of terms to allocate. n_term < 0 => em_taylor.term pointer will be disassociated.
  save_old : bool, optional
      If True then save any old terms when em_taylor is resized. Default is False.
  )""");
  m.def(
      "init_lat",
      &Bmad::init_lat,
      py::arg("n") = py::none(),
      py::arg("init_beginning_ele") = py::none(),
      R"""(Parameters
  ----------
  lat : LatStruct
      Initialized lat.
  n : int, optional
      Upper bound lat.ele(0:) array is initialized to. Default is 10.
  init_beginning_ele : bool, optional
      Init lat.ele(0)? Default is False.
  )""");
  m.def(
      "init_multipole_cache",
      &Bmad::init_multipole_cache,
      py::arg("ele"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to init
      This parameter is an input/output and is modified in-place. As an output: Initalized element.
  )""");
  m.def(
      "init_photon_from_a_photon_init_ele",
      &Bmad::init_photon_from_a_photon_init_ele,
      py::arg("ele"),
      py::arg("param"),
      py::arg("random_on") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      patch element.
  param : 
      lat_param_struct.
  orbit : CoordStruct
      Output photon coords.
  random_on : bool, optional
      : Default is True. If False then use zero for all random numbers needed in the calc.
  )""");
  py::class_<
      Bmad::InitPhotonIntegProb,
      std::unique_ptr<Bmad::InitPhotonIntegProb>>(
      m, "InitPhotonIntegProb", "init_photon_integ_prob return type")
      .def_readonly("E_photon", &Bmad::InitPhotonIntegProb::E_photon)
      .def_readonly("integ_prob", &Bmad::InitPhotonIntegProb::integ_prob)
      .def("__len__", [](const Bmad::InitPhotonIntegProb&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::InitPhotonIntegProb& s, int i) -> py::object {
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
  vert_angle_symmetric : float, optional
      Use two symmetric ranges [-vert_angle_max, -vert_angle_min] and [vert_angle_min, vert_angle_max] instead
      of just [vert_angle_min, vert_angle_max]?
  energy_integ_prob : float, optional
      If present, E_photon will be set to the photon energy such that the integrated probability of generating a
      photon in the given angle and energy range in the interval [E_min, E_photon] is energy_integ_prob. That
      is, energy_integ_prob = 0 => E_photon = E_min and energy_integ_prob = 1 => E_photon = E_max.

  Returns
  -------
  E_photon : float
      See energy_integ_prob. E_photon must be present if energy_integ_prob is.
  integ_prob : float
      Integrated probablility of emitting a photon in given angle and energy range.
  )""");
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
      Initialization parameters .spin(3)  -- (x, y, z) spin coordinates ele

  Returns
  -------
  bunch : BunchStruct
      Bunch of particles. .particle(:).spin
  )""");
  py::class_<PyInitSurfaceSegment, std::unique_ptr<PyInitSurfaceSegment>>(
      m, "InitSurfaceSegment", "init_surface_segment return type")
      .def_readonly("ix", &PyInitSurfaceSegment::ix)
      .def_readonly("iy", &PyInitSurfaceSegment::iy)
      .def("__len__", [](const PyInitSurfaceSegment&) { return 2; })
      .def(
          "__getitem__",
          [](const PyInitSurfaceSegment& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.ix);
            if (i == 1)
              return py::cast(s.iy);
            throw py::index_error();
          });
  m.def(
      "init_surface_segment",
      &python_init_surface_segment,
      py::arg("phot"),
      py::arg("ix"),
      py::arg("iy"),
      R"""(Subroutine init_surface_segment (phot, ix, iy)

  Routine to init the componentes in ele%photon%segmented%pt(ix,iy) for use with segmented surface calculations.

  Parameters
  ----------
  phot : unknown
      index of grid point to init.
  )""");
  m.def(
      "init_taylor_series",
      &Bmad::init_taylor_series,
      py::arg("bmad_taylor"),
      py::arg("n_term"),
      py::arg("save_old") = py::none(),
      R"""(Parameters
  ----------
  bmad_taylor : TaylorStruct
      Old structure.
      This parameter is an input/output and is modified in-place. As an output: Initalized structure.
  n_term : int
      Number of terms to allocate. n_term < 0 => bmad_taylor.term pointer will be disassociated.
  save_old : bool, optional
      If True then save any old terms and ref orbit when bmad_taylor is resized. If False zero the ref orbit.
      Default is False.
  )""");
  m.def(
      "init_wake",
      &Bmad::init_wake,
      py::arg("n_sr_long"),
      py::arg("n_sr_trans"),
      py::arg("n_sr_z"),
      py::arg("n_lr_mode"),
      py::arg("always_allocate") = py::none(),
      R"""(Parameters
  ----------
  wake : WakeStruct
      Initialized structure.
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
  )""");
  m.def(
      "insert_element",
      &Bmad::insert_element,
      py::arg("lat"),
      py::arg("insert_ele"),
      py::arg("ix_ele"),
      py::arg("ix_branch") = py::none(),
      py::arg("orbit") = py::none(),
      R"""(Parameters
  ----------
  lat : LatStruct
      lattice that will be modified
      This parameter is an input/output and is modified in-place. As an output: lattice with new element
      inserted
  insert_ele : EleStruct
      element to insert into the lat
  ix_ele : int
      branch.ele(:) index where the new element is inserted.
  ix_branch : int, optional
      : branch index for the insertion. Default = 0.
  orbit : CoordStruct, optional
      orbit array to enlarge.
      This parameter is an input/output and is modified in-place. As an output: Enlarged orbit array.
  )""");
  py::class_<PyIntegrandBase, std::unique_ptr<PyIntegrandBase>>(
      m, "IntegrandBase", "integrand_base return type")
      .def_readonly("func_retval__", &PyIntegrandBase::func_retval__)
      .def("__len__", [](const PyIntegrandBase&) { return 1; })
      .def("__getitem__", [](const PyIntegrandBase& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.func_retval__);
        throw py::index_error();
      });
  m.def(
      "integrand_base",
      &python_integrand_base,
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
      Array of reals over which to evaluate the integrand. <return value> -- REAL(rp): Array of reals containing
      values of integrand at t(:).
  )""");
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
  args : float
      Parameters and constants of DEQ.  See psi_prime comments for details.

  Returns
  -------
  result : float
      Integral of psi from -bound to +bound.
  )""");
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

  )""");
  py::class_<PyIntegrationTimerEle, std::unique_ptr<PyIntegrationTimerEle>>(
      m, "IntegrationTimerEle", "integration_timer_ele return type")
      .def_readonly("tol", &PyIntegrationTimerEle::tol)
      .def("__len__", [](const PyIntegrationTimerEle&) { return 1; })
      .def(
          "__getitem__",
          [](const PyIntegrationTimerEle& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.tol);
            throw py::index_error();
          });
  m.def(
      "integration_timer",
      py::overload_cast<
          EleProxy&,
          LatParamProxy&,
          CoordProxy&,
          CoordProxy&,
          double>(&python_integration_timer_ele),
      py::arg("ele"),
      py::arg("param"),
      py::arg("start"),
      py::arg("orb_max"),
      py::arg("tol"),
      R"""(Parameters
  ----------
  ele : 
  param : 
  start : 
  orb_max : 
  tol : 
  )""");
  m.def(
      "integration_timer",
      py::overload_cast<
          FibreRawStruct&,
          FixedArray1D<Real, 6>,
          FixedArray1D<Real, 6>,
          double>(&Bmad::integration_timer),
      py::arg("a_fibre"),
      py::arg("orbit"),
      py::arg("orbit_max"),
      py::arg("tol_dp"),
      R"""(Parameters
  ----------
  a_fibre : 
  orbit : 
  orbit_max : 
  tol_dp : 
  )""");
  m.def(
      "ion_kick",
      &Bmad::ion_kick,
      py::arg("orbit"),
      py::arg("r_beam"),
      py::arg("n_beam_part"),
      py::arg("a_twiss"),
      py::arg("b_twiss"),
      py::arg("sig_ee"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Ion position.
  r_beam : float
      Beam (x, y) position.
  n_beam_part : float
      Number of beam particles.
  a_twiss : TwissStruct
      Horizontal like beam twiss parameters.
  b_twiss : TwissStruct
      vertical like beam twiss parameters.
  sig_ee : float
      Sigma_E/E beam energy spread.
  kick : float
      (x, y, s) kick in m/sec.
  )""");
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
  )""");
}
