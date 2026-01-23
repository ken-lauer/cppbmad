#include "pybmad/generated/Bmad_routines_r.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyRadiationIntegrals python_radiation_integrals(
    LatStruct &lat,
    CoordStructArray1D orbit,
    std::optional<int> ix_cache = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
) {
  auto _result = Bmad::radiation_integrals(lat, orbit, make_opt_ref(ix_cache), ix_branch);
  auto py_result{PyRadiationIntegrals{_result, ix_cache}};
  return py_result;
}
PyReleaseRadIntCache python_release_rad_int_cache(int ix_cache) {
  Bmad::release_rad_int_cache(ix_cache);
  auto py_result{PyReleaseRadIntCache{ix_cache}};
  return py_result;
}

void init_Bmad_routines_r(py::module &m) {
  py::class_<Bmad::Rad1DampAndStocMats, std::unique_ptr<Bmad::Rad1DampAndStocMats>>(
      m,
      "Rad1DampAndStocMats",
      "rad1_damp_and_stoc_mats return type"
  )
      .def_readonly("rad_map", &Bmad::Rad1DampAndStocMats::rad_map)
      .def_readonly("err_flag", &Bmad::Rad1DampAndStocMats::err_flag)
      .def_readonly("rad_int1", &Bmad::Rad1DampAndStocMats::rad_int1)
      .def("__len__", [](const Bmad::Rad1DampAndStocMats &) { return 3; })
      .def("__getitem__", [](const Bmad::Rad1DampAndStocMats &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.rad_map);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.rad_int1);
        throw py::index_error();
      });
  m.def(
      "rad1_damp_and_stoc_mats",
      &Bmad::rad1_damp_and_stoc_mats,
      py::arg("ele"),
      py::arg("include_opening_angle"),
      py::arg("orb_in"),
      py::arg("orb_out"),
      py::arg("g2_tol"),
      py::arg("g3_tol"),
      py::arg("ele0") = py::none(),
      R"""(Subroutine rad1_damp_and_stoc_mats (ele, include_opening_angle, orb_in, orb_out, rad_map, g2_tol, g3_tol, err_flag, ele0, rad_int1)

  Routine to calculate the damping and stochastic matrices for a given lattice element.

  Parameters
  ----------
  ele : EleStruct
      Element under consideration.

  include_opening_angle : bool
      If True include the effect of the vertical opening angle of emitted radiation.

  Generally use True unless comparing against other codes. : 
  orb_in : CoordStruct
      Entrance orbit about which to compute the matrices.

  orb_out : CoordStruct
      Exit orbit.

  g2_tol : float
      Tollerance on g^2 per unit length (damping tolerance).

  g3_tol : float
      Tollerance on g^3 per unit length (stocastic tolerance).

  ele0 : EleStruct, optional
      Element before `ele`. Needed if and only if rad_int1 is present

  Returns
  -------
  rad_map : RadMapStruct
      Damping and stochastic matrices.

  %stoc_mat : 
      Variance matrix.

  err_flag : bool
      Set true if there is an error. False otherwise.

  rad_int1 : RadInt1Struct
      Radiation integrals
  )"""
  );
  py::class_<Bmad::RadDampAndStocMats, std::unique_ptr<Bmad::RadDampAndStocMats>>(
      m,
      "RadDampAndStocMats",
      "rad_damp_and_stoc_mats return type"
  )
      .def_readonly("rmap", &Bmad::RadDampAndStocMats::rmap)
      .def_readonly("mode", &Bmad::RadDampAndStocMats::mode)
      .def_readonly("xfer_nodamp_mat", &Bmad::RadDampAndStocMats::xfer_nodamp_mat)
      .def_readonly("err_flag", &Bmad::RadDampAndStocMats::err_flag)
      .def_readonly("rad_int_branch", &Bmad::RadDampAndStocMats::rad_int_branch)
      .def("__len__", [](const Bmad::RadDampAndStocMats &) { return 5; })
      .def("__getitem__", [](const Bmad::RadDampAndStocMats &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.rmap);
        if (i == 1)
          return py::cast(s.mode);
        if (i == 2)
          return py::cast(s.xfer_nodamp_mat);
        if (i == 3)
          return py::cast(s.err_flag);
        if (i == 4)
          return py::cast(s.rad_int_branch);
        throw py::index_error();
      });
  m.def(
      "rad_damp_and_stoc_mats",
      &Bmad::rad_damp_and_stoc_mats,
      py::arg("ele1"),
      py::arg("ele2"),
      py::arg("include_opening_angle"),
      py::arg("closed_orbit") = py::none(),
      R"""(Subroutine rad_damp_and_stoc_mats (ele1, ele2, include_opening_angle, rmap, mode, xfer_nodamp_mat, err_flag, closed_orbit, rad_int_branch)

  Routine to calculate the damping and stochastic variance matrices from exit end of ele1
  to the exit end of ele2. Use ele1 = ele2 to get 1-turn matrices.

  If ele2 is before ele1 the integration range if from ele1 to the branch end plus
  from the beginning to ele2.

  Note: The ele%mat6 matrices will be remade. By convention, these matrices
  do not include damping.

  Parameters
  ----------
  ele1 : EleStruct
      Start element of integration range.

  ele2 : EleStruct
      End element of integration range.

  include_opening_angle : bool
      If True include the effect of the vertical opening angle of emitted radiation.

  Generally use True unless comparing against other codes. : 
  closed_orbit : CoordStruct, optional
      Closed orbit. If not present this routine will calculate it.

  Returns
  -------
  rmap : RadMapStruct
      Damping and stochastic mats

  %stoc_mat : 
      stochastic variance matrix.

  mode : NormalModesStruct

  %dpz_damp : 
      Change in pz without RF.

  %pz_average : 
      Average pz due to damping.

  xfer_nodamp_mat : float
      Transfer matrix without damping.

  rad_int_branch : RadIntBranchStruct
      Array of element-by-element radiation integrals.

  err_flag : bool
      Set true if there is a problem.
  )"""
  );
  py::class_<Bmad::RadGIntegrals, std::unique_ptr<Bmad::RadGIntegrals>>(
      m,
      "RadGIntegrals",
      "rad_g_integrals return type"
  )
      .def_readonly("int_g", &Bmad::RadGIntegrals::int_g)
      .def_readonly("int_g3", &Bmad::RadGIntegrals::int_g3)
      .def("__len__", [](const Bmad::RadGIntegrals &) { return 2; })
      .def("__getitem__", [](const Bmad::RadGIntegrals &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.int_g);
        if (i == 1)
          return py::cast(s.int_g3);
        throw py::index_error();
      });
  m.def(
      "rad_g_integrals",
      &Bmad::rad_g_integrals,
      py::arg("ele"),
      py::arg("where"),
      py::arg("orb_in"),
      py::arg("orb_out"),
      py::arg("int_g2"),
      py::arg("g_tol"),
      py::arg("g2_tol"),
      py::arg("g3_tol"),
      R"""(Subroutine rad_g_integrals (ele, where, orb_in, orb_out, int_g, int_g2, int_g3, g_tol, g2_tol, g3_tol)

  Routine to calculate bending strength integrals (g(s) = 1/trajectory_bending_radius(s)) in
  laboratory coords.

  Parameters
  ----------
  ele : EleStruct
      Element under consideration.

  where : int
      What part of ele to integrate over.

  upstream$ -> 1st half of element : 
  downsteam$ -> 2nd half : 
  all$ -> everything. : 
  orb_in : CoordStruct
      Entrance orbit about which to compute the matrices.

  orb_out : CoordStruct
      Exit orbit.

  g_tol : float
      Tollerance on |g| per unit length.

  g2_tol : float
      Tollerance on g^2 per unit length.

  g3_tol : float
      Tollerance on g^3 per unit length.

  Returns
  -------
  int_g : float
      Integrals of (gx,gy) vector.

  gint_g2 : float
      integrals of |g|^2 and |g|^3.

  int_g3 : float
      integrals of |g|^2 and |g|^3.
  )"""
  );
  py::class_<PyRadiationIntegrals, std::unique_ptr<PyRadiationIntegrals>>(
      m,
      "RadiationIntegrals",
      "radiation_integrals return type"
  )
      .def_readonly("mode", &PyRadiationIntegrals::mode)
      .def_readonly("rad_int_by_ele", &PyRadiationIntegrals::rad_int_by_ele)
      .def_readonly("ix_cache", &PyRadiationIntegrals::ix_cache)
      .def("__len__", [](const PyRadiationIntegrals &) { return 3; })
      .def("__getitem__", [](const PyRadiationIntegrals &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.mode);
        if (i == 1)
          return py::cast(s.rad_int_by_ele);
        if (i == 2)
          return py::cast(s.ix_cache);
        throw py::index_error();
      });
  m.def(
      "radiation_integrals",
      &python_radiation_integrals,
      py::arg("lat"),
      py::arg("orbit"),
      py::arg("ix_cache") = py::none(),
      py::arg("ix_branch") = py::none(),
      R"""(Wrapper for Fortran routine radiation_integrals

  Parameters
  ----------
  lat : LatStruct
      Lattice to use. The calculation assumes that

  orbit : CoordStruct
      Closed orbit for the branch.

  mode : NormalModesStruct
      Parameters for the ("horizontal like") a-mode,

  ix_cache : int, optional
      Cache pointer.
      This parameter is an input/output and is modified in-place. As an output: Cache pointer. If ix_cache = 0
      at input then

  ix_branch : int, optional
      Lattice branch index. Default is 0.

  rad_int_by_ele : RadIntAllEleStruct
      Radiation integrals element by element.
  )"""
  );
  m.def(
      "radiation_map_setup",
      &Bmad::radiation_map_setup,
      py::arg("ele"),
      py::arg("ref_orbit_in") = py::none(),
      R"""(Subroutine radiation_map_setup (ele, err_flag, ref_orbit_in)

  Routine to calculate the radiation kick for a lattice element.

  Parameters
  ----------
  ele : EleStruct
      Element whose map is to be setup.
      This parameter is an input/output and is modified in-place. As an output: Element with map calculated.

  ref_orb : CoordStruct, optional
      If present, ignore ele_map.stale setting and make the map around this reference orbit.

  Returns
  -------
  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "ramper_slave_setup",
      &Bmad::ramper_slave_setup,
      py::arg("lat"),
      py::arg("force_setup") = py::none(),
      R"""(Wrapper for Fortran routine ramper_slave_setup

  Parameters
  ----------
  lat : LatStruct
      Lattice to be setup.
      This parameter is an input/output and is modified in-place. As an output: Lattice with ramper slaves
      setup.

  force_setup : bool, optional
      Default False. If True, do the setup even if lat.ramper_slave_bookkeeping = ok$. But the setup will never
      be done if lat.ramper_slave_bookkeeping = super_ok$.
  )"""
  );
  py::class_<Bmad::RamperValue, std::unique_ptr<Bmad::RamperValue>>(
      m,
      "RamperValue",
      "ramper_value return type"
  )
      .def_readonly("err_flag", &Bmad::RamperValue::err_flag)
      .def_readonly("value", &Bmad::RamperValue::value)
      .def("__len__", [](const Bmad::RamperValue &) { return 2; })
      .def("__getitem__", [](const Bmad::RamperValue &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "ramper_value",
      &Bmad::ramper_value,
      py::arg("ramper"),
      py::arg("r1"),
      R"""(Wrapper for Fortran routine ramper_value

  Parameters
  ----------
  ramper : EleStruct
      Ramper lord.

  r1 : ControlRamp1Struct
      Slave function.

  err_flag : bool
      Set True if there is an error, False otherwise.

  value : float
      Value of the slave function.
  )"""
  );
  m.def(
      "randomize_lr_wake_frequencies",
      &Bmad::randomize_lr_wake_frequencies,
      py::arg("ele"),
      R"""(Subroutine randomize_lr_wake_frequencies (ele, set_done)

  Routine to randomize the frequencies of the lr wake HOMs according to:
    freq = freq_in * (1 + lr_freq_spread) * rr)
  where rr is a Gaussian distributed random number with unit variance.

  Parameters
  ----------
  ele : EleStruct
      Element with wake. If no wake then nothing is done.
      This parameter is an input/output and is modified in-place. As an output: Element with wake frequencies
      set.

  %value : 
      Frequency.

  Returns
  -------
  %wake%lr%mode : 
      Set frequency.

  set_done : bool
      Set True if there where lr wakes to be set. False otherwise.
  )"""
  );
  m.def(
      "rchomp",
      &Bmad::rchomp,
      py::arg("rel"),
      py::arg("plc"),
      py::arg("out"),
      R"""(Wrapper for Fortran routine rchomp

  Parameters
  ----------
  rel : 
  plc : 
  out : 
  )"""
  );
  m.def(
      "re_allocate_eles",
      &Bmad::re_allocate_eles,
      py::arg("eles"),
      py::arg("n"),
      py::arg("save_old") = py::none(),
      py::arg("exact") = py::none(),
      R"""(Wrapper for Fortran routine re_allocate_eles

  Parameters
  ----------
  eles : ElePointerStruct
      Array of element pointers with possible old data.
      This parameter is an input/output and is modified in-place. As an output: Array of element pointers.

  n : int
      Array size to set.

  save_old : bool, optional
      If present and True then save the old data.

  exact : bool, optional
      If present and True then eles will have size = n If False (default), reallcation will not be done if eles
      is already large enough
  )"""
  );
  m.def(
      "re_allocate",
      py::overload_cast<Wall3dSectionStructAlloc1D, int, std::optional<bool>>(&Bmad::re_allocate),
      py::arg("section"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      R"""(Wrapper for Fortran routine re_allocate_wall3d_section_array

  Parameters
  ----------
  section : 
  n : 
  exact : 
  )"""
  );
  m.def(
      "re_allocate",
      py::overload_cast<Wall3dVertexStructAlloc1D, int, std::optional<bool>>(&Bmad::re_allocate),
      py::arg("v"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      R"""(Wrapper for Fortran routine re_allocate_wall3d_vertex_array

  Parameters
  ----------
  v : 
  n : 
  exact : 
  )"""
  );
  m.def(
      "re_associate_node_array",
      &Bmad::re_associate_node_array,
      py::arg("tree"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      R"""(Subroutine re_associate_node_array(tree, n, exact)

  Routine to resize the tree%node(:) array.

  Note: The data of the array is preserved but data at the end of the
  array will be lost if n is less than the original size of the array

  Parameters
  ----------
  tree : ExpressionTreeStruct

  n : int
      Size wanted.

  exact : bool, optional
      Default is False. If False, the size of the output array is permitted to be larger than n.
  )"""
  );
  m.def(
      "re_str",
      py::overload_cast<long double, std::string>(&Bmad::re_str),
      py::arg("rel"),
      py::arg("str_out"),
      R"""(Wrapper for Fortran routine re_str_qp

  Parameters
  ----------
  rel : 
  str_out : 
  )"""
  );
  m.def(
      "re_str",
      py::overload_cast<double, std::string>(&Bmad::re_str),
      py::arg("rel"),
      py::arg("str_out"),
      R"""(Wrapper for Fortran routine re_str_rp

  Parameters
  ----------
  rel : 
  str_out : 
  )"""
  );
  py::class_<Bmad::ReadBeamAscii, std::unique_ptr<Bmad::ReadBeamAscii>>(
      m,
      "ReadBeamAscii",
      "read_beam_ascii return type"
  )
      .def_readonly("beam", &Bmad::ReadBeamAscii::beam)
      .def_readonly("err_flag", &Bmad::ReadBeamAscii::err_flag)
      .def("__len__", [](const Bmad::ReadBeamAscii &) { return 2; })
      .def("__getitem__", [](const Bmad::ReadBeamAscii &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.beam);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "read_beam_ascii",
      &Bmad::read_beam_ascii,
      py::arg("file_name"),
      py::arg("beam_init"),
      R"""(Subroutine read_beam_ascii (file_name, beam, beam_init, err_flag, ele, print_mom_shift_warning, conserve_momentum)

  Subroutine to read in a beam definition file.
  If non_zero, the following components of beam_init are used to rescale the beam:
      %n_bunch
      %n_particle
      %charge_tot

  If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
  Otherwise the file is assumed to be ASCII.

  Parameters
  ----------
  iu : int
      File unit number

  file_name : unknown
      Name of beam file.

  beam_init : BeamInitStruct
      See above.

  ele : EleStruct, optional
      Element with reference energy, etc.

  print_mom_shift_warning : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

  shift_momentum : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

  Returns
  -------
  beam : BeamStruct
      Structure holding the beam information.

  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  py::class_<Bmad::ReadBeamFile, std::unique_ptr<Bmad::ReadBeamFile>>(
      m,
      "ReadBeamFile",
      "read_beam_file return type"
  )
      .def_readonly("beam", &Bmad::ReadBeamFile::beam)
      .def_readonly("err_flag", &Bmad::ReadBeamFile::err_flag)
      .def("__len__", [](const Bmad::ReadBeamFile &) { return 2; })
      .def("__getitem__", [](const Bmad::ReadBeamFile &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.beam);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "read_beam_file",
      &Bmad::read_beam_file,
      py::arg("file_name"),
      py::arg("beam_init"),
      py::arg("ele") = py::none(),
      py::arg("print_mom_shift_warning") = py::none(),
      py::arg("conserve_momentum") = py::none(),
      R"""(Subroutine read_beam_file (file_name, beam, beam_init, err_flag, ele, print_mom_shift_warning, conserve_momentum)

  Subroutine to read in a beam definition file.
  If non_zero, the following components of beam_init are used to rescale the beam:
      %n_bunch
      %n_particle
      %bunch_charge -> charge_tot
      %species

  If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
  Otherwise the file is assumed to be ASCII.

  Parameters
  ----------
  file_name : unknown
      Name of beam file.

  beam_init : BeamInitStruct
      See above.

  ele : EleStruct, optional
      Element with reference energy, etc.

  print_mom_shift_warning : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

  shift_momentum : bool, optional
      Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

  Returns
  -------
  beam : BeamStruct
      Structure holding the beam information.

  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "read_binary_cartesian_map",
      &Bmad::read_binary_cartesian_map,
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("cart_map"),
      py::arg("err_flag"),
      R"""(Subroutine read_binary_cartesian_map (file_name, ele, cart_map, err_flag)

  Routine to read a binary cartesian_map structure.

  Parameters
  ----------
  file_name : unknown
      File to create.

  ele : EleStruct
      Element associated with the map.

  Ouput: : 
  cart_map : 
      cartesian_map_struct, cartesian map.

  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "read_binary_cylindrical_map",
      &Bmad::read_binary_cylindrical_map,
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("cl_map"),
      py::arg("err_flag"),
      R"""(Subroutine read_binary_cylindrical_map (file_name, ele, cl_map, err_flag)

  Routine to read a binary cylindrical_map structure.

  Parameters
  ----------
  file_name : unknown
      File to create.

  ele : EleStruct
      Element associated with the map.

  Ouput: : 
  cl_map : 
      cylindrical_map_struct, cylindrical map.

  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "read_binary_grid_field",
      &Bmad::read_binary_grid_field,
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("g_field"),
      py::arg("err_flag"),
      R"""(Subroutine read_binary_grid_field (file_name, ele, g_field, err_flag)

  Routine to read a binary grid_field structure.

  Parameters
  ----------
  file_name : unknown
      File to create.

  ele : EleStruct
      Element associated with the map.

  Ouput: : 
  g_field : 
      grid_field_struct, cylindrical map.

  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "read_surface_reflection_file",
      &Bmad::read_surface_reflection_file,
      py::arg("file_name"),
      R"""(Subroutine read_surface_reflection_file (file_name, surface)

  Routine to read the reflection probability data for a given type of surface from a file.

  Parameters
  ----------
  file_name : unknown
      Name of the file.

  Returns
  -------
  surface : PhotonReflectSurfaceStruct
      Surface info.
  )"""
  );
  m.def(
      "reallocate_beam",
      &Bmad::reallocate_beam,
      py::arg("beam"),
      py::arg("n_bunch"),
      py::arg("n_particle") = py::none(),
      py::arg("extend") = py::none(),
      R"""(Wrapper for Fortran routine reallocate_beam

  Parameters
  ----------
  beam : BeamStruct
      Beam bunches are saved if save = True.
      This parameter is an input/output and is modified in-place. As an output: Allocated beam_struct structure.

  n_bunch : int
      Number of bunches.

  n_particle : int, optional
      Number of particles. Must be non-negative. If save = True then the number of particles in existing bunches
      will not be touched.

  extend : 
  )"""
  );
  m.def(
      "reallocate_bp_com_const",
      &Bmad::reallocate_bp_com_const,
      R"""(Wrapper for Fortran routine reallocate_bp_com_const
  )"""
  );
  m.def(
      "reallocate_bunch",
      &Bmad::reallocate_bunch,
      py::arg("n_particle"),
      py::arg("save") = py::none(),
      R"""(Wrapper for Fortran routine reallocate_bunch

  Parameters
  ----------
  bunch : BunchStruct
      Allocated bunch_struct structure.

  n_particle : int
      Number of particles. Must be non-negative.

  save : bool, optional
      If present and True then save the old bunch info.
  )"""
  );
  m.def(
      "reallocate_control",
      &Bmad::reallocate_control,
      py::arg("lat"),
      py::arg("n"),
      R"""(Wrapper for Fortran routine reallocate_control

  Parameters
  ----------
  lat : LatStruct
      Lattice.

  n : int
      Array size for lat.control(:) and lat.ic(:).
  )"""
  );
  m.def(
      "reallocate_coord",
      py::overload_cast<CoordArrayStructAlloc1D, LatStruct &>(&Bmad::reallocate_coord),
      py::arg("coord_array"),
      py::arg("lat"),
      R"""(Subroutine reallocate_coord (...)

  Routine to allocate or reallocate at allocatable coord_struct array.
  reallocate_coord is an overloaded name for:
    reallocate_coord_n (coord, n_coord)
    reallocate_coord_lat (coord, lat, ix_branch)

  Subroutine to allocate an allocatable coord_struct array to at least:
      coord(0:n_coord)                            if n_coord arg is used.
      coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

  The old coordinates are saved
  If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
  In any case, coord(n)%vec for n > 0 is set to zero.

  Parameters
  ----------
  coord : CoordStruct
      Allocatable array.
      This parameter is an input/output and is modified in-place. As an output: Allocated array.

  n_coord : int
      Minimum array upper bound wanted.

  lat : LatStruct
      Lattice

  ix_branch : int, optional
      Branch to use. Default is 0 (main branch).

  Notes
  -----
  Overloaded versions:
  )"""
  );
  m.def(
      "reallocate_coord",
      py::overload_cast<CoordStructAlloc1D, LatStruct &, std::optional<int>>(&Bmad::reallocate_coord
      ),
      py::arg("coord"),
      py::arg("lat"),
      py::arg("ix_branch") = py::none(),
      R"""(Subroutine reallocate_coord (...)

  Routine to allocate or reallocate at allocatable coord_struct array.
  reallocate_coord is an overloaded name for:
    reallocate_coord_n (coord, n_coord)
    reallocate_coord_lat (coord, lat, ix_branch)

  Subroutine to allocate an allocatable coord_struct array to at least:
      coord(0:n_coord)                            if n_coord arg is used.
      coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

  The old coordinates are saved
  If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
  In any case, coord(n)%vec for n > 0 is set to zero.

  Parameters
  ----------
  coord : CoordStruct
      Allocatable array.
      This parameter is an input/output and is modified in-place. As an output: Allocated array.

  n_coord : int
      Minimum array upper bound wanted.

  lat : LatStruct
      Lattice

  ix_branch : int, optional
      Branch to use. Default is 0 (main branch).

  Notes
  -----
  Overloaded versions:
  )"""
  );
  m.def(
      "reallocate_coord",
      py::overload_cast<CoordStructAlloc1D, int>(&Bmad::reallocate_coord),
      py::arg("coord"),
      py::arg("n_coord"),
      R"""(Subroutine reallocate_coord (...)

  Routine to allocate or reallocate at allocatable coord_struct array.
  reallocate_coord is an overloaded name for:
    reallocate_coord_n (coord, n_coord)
    reallocate_coord_lat (coord, lat, ix_branch)

  Subroutine to allocate an allocatable coord_struct array to at least:
      coord(0:n_coord)                            if n_coord arg is used.
      coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

  The old coordinates are saved
  If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
  In any case, coord(n)%vec for n > 0 is set to zero.

  Parameters
  ----------
  coord : CoordStruct
      Allocatable array.
      This parameter is an input/output and is modified in-place. As an output: Allocated array.

  n_coord : int
      Minimum array upper bound wanted.

  lat : LatStruct
      Lattice

  ix_branch : int, optional
      Branch to use. Default is 0 (main branch).

  Notes
  -----
  Overloaded versions:
  )"""
  );
  m.def(
      "reallocate_expression_stack",
      &Bmad::reallocate_expression_stack,
      py::arg("stack"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      R"""(Wrapper for Fortran routine reallocate_expression_stack

  Parameters
  ----------
  stack : unknown
      Existing stack array.
      This parameter is an input/output and is modified in-place. As an output: Resized stack.

  n : int
      Array size needed.

  exact : bool, optional
      If present and False then the size of the output array is permitted to be larger than n. Default is True.
  )"""
  );
  m.def(
      "rel_tracking_charge_to_mass",
      &Bmad::rel_tracking_charge_to_mass,
      py::arg("orbit"),
      py::arg("ref_species"),
      R"""(Wrapper for Fortran routine rel_tracking_charge_to_mass

  Parameters
  ----------
  orbit : CoordStruct
      Particle position structure.

  ref_species : int
      Reference species

  rel_charge : float
      Relative charge/mass
  )"""
  );
  m.def(
      "relative_mode_flip",
      &Bmad::relative_mode_flip,
      py::arg("ele1"),
      py::arg("ele2"),
      py::arg("func_retval__"),
      R"""(Wrapper for Fortran routine relative_mode_flip

  Parameters
  ----------
  ele1 : EleStruct
      Elements to compare.

  ele2 : EleStruct
      Elements to compare.

  relative_mode_flip : 
  )"""
  );
  py::class_<PyReleaseRadIntCache, std::unique_ptr<PyReleaseRadIntCache>>(
      m,
      "ReleaseRadIntCache",
      "release_rad_int_cache return type"
  )
      .def_readonly("ix_cache", &PyReleaseRadIntCache::ix_cache)
      .def("__len__", [](const PyReleaseRadIntCache &) { return 1; })
      .def("__getitem__", [](const PyReleaseRadIntCache &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.ix_cache);
        throw py::index_error();
      });
  m.def(
      "release_rad_int_cache",
      &python_release_rad_int_cache,
      py::arg("ix_cache"),
      R"""(Subroutine release_rad_int_cache (ix_cache)

  Subroutine to release the memory associated with caching wiggler values.
  See the radiation_integrals routine for further details.

  Parameters
  ----------
  ix_cache : int
      Cache number.
      This parameter is an input/output and is modified in-place. As an output: Cache number set to 0,
  )"""
  );
  m.def(
      "remove_constant_taylor",
      &Bmad::remove_constant_taylor,
      py::arg("taylor_in"),
      py::arg("taylor_out"),
      py::arg("c0"),
      py::arg("remove_higher_order_terms"),
      R"""(Subroutine remove_constant_taylor (taylor_in, taylor_out, c0, remove_higher_order_terms)

  Subroutine to remove the constant part of a taylor map.
  Optionally terms that are higher order than bmad_com%taylor_order can
  be removed.

  Note: It is assumed that taylor_out has been deallocated before the call to
  this routine. Calling this routine with the first two actual arguments the
  same is prohibited.

  Parameters
  ----------
  taylor_in : TaylorStruct
      Input taylor map.

  remove_higher_order_terms : bool
      If True then terms that are higher order than bmad_com.taylor_order are removed.

  Returns
  -------
  taylor_out : TaylorStruct
      Taylor with constant terms removed.

  c0 : float
      The constant part of the taylor map
  )"""
  );
  m.def(
      "remove_dead_from_bunch",
      &Bmad::remove_dead_from_bunch,
      py::arg("bunch_in"),
      R"""(Wrapper for Fortran routine remove_dead_from_bunch

  Parameters
  ----------
  bunch_in : BunchStruct
      Input bunch with alive and dead particles.

  bunch_out : BunchStruct
      Output bunch with only alive and pre_born particles. Note: bunch_out can be the same actual argument as
      bunch_in.
  )"""
  );
  m.def(
      "remove_eles_from_lat",
      &Bmad::remove_eles_from_lat,
      py::arg("lat"),
      py::arg("check_sanity") = py::none(),
      R"""(Wrapper for Fortran routine remove_eles_from_lat

  Parameters
  ----------
  lat : LatStruct
      Lattice to compress.
      This parameter is an input/output and is modified in-place. As an output: Compressed lattice.

  check_sanity : bool, optional
      If True (default) then call lat_sanity_check
  )"""
  );
  m.def(
      "remove_lord_slave_link",
      &Bmad::remove_lord_slave_link,
      py::arg("lord"),
      py::arg("slave"),
      R"""(Wrapper for Fortran routine remove_lord_slave_link

  Parameters
  ----------
  lord : EleStruct
      Lord element
      This parameter is an input/output and is modified in-place. As an output: Lord element with link info
      removed

  slave : EleStruct
      Slave element
      This parameter is an input/output and is modified in-place. As an output: Slave element with link info
      removed
  )"""
  );
  m.def(
      "reverse_lat",
      &Bmad::reverse_lat,
      py::arg("lat_in"),
      py::arg("track_antiparticle") = py::none(),
      R"""(Wrapper for Fortran routine reverse_lat

  Parameters
  ----------
  lat_in : LatStruct
      Input lattice to reverse.

  lat_rev : LatStruct
      Reversed lattice.

  track_antiparticle : bool, optional
      Set the particle species of the reversed lat to the anti-particle of lat_in? Default is True.
  )"""
  );
  m.def(
      "rf_coupler_kick",
      &Bmad::rf_coupler_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("phase"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(No longer in the codebase
  function rf_clock_setup (branch, n_rf_included, n_rf_excluded) result (ok)
    import
    implicit none
    type (branch_struct), target :: branch
    integer n_rf_included, n_rf_excluded
    logical ok
  end function

  Parameters
  ----------
  ele : EleStruct
      Element being tracked through

  param : LatParamStruct
      branch parameters.

  particle_at : int
      first_track_edge$, or second_track_edge$.

  phase : float
      phase of cavity

  orbit : CoordStruct
      Position before kick.
      This parameter is an input/output and is modified in-place. As an output: Position after kick.

  mat6 : float, optional
      Transfer matrix before the element.
      This parameter is an input/output and is modified in-place. As an output: Transfer matrix through the
      element.

  make_matrix : bool, optional
      Propagate the transfer matrix? Default is false.

  Returns
  -------
  ok
  )"""
  );
  m.def(
      "rf_is_on",
      &Bmad::rf_is_on,
      py::arg("branch"),
      py::arg("ix_ele1") = py::none(),
      py::arg("ix_ele2") = py::none(),
      R"""(Wrapper for Fortran routine rf_is_on

  Parameters
  ----------
  branch : BranchStruct
      Lattice branch to check.

  ix_ele1 : int, optional
      Start of range of elements to check. Default is 0.

  ix_ele2 : int, optional
      End of range of elements to check. Default is branch.n_ele_track.

  is_on : bool
      True if any rfcavity is powered. False otherwise.
  )"""
  );
  m.def(
      "rf_ref_time_offset",
      &Bmad::rf_ref_time_offset,
      py::arg("ele"),
      py::arg("ds") = py::none(),
      py::arg("time"),
      R"""(Wrapper for Fortran routine rf_ref_time_offset

  Parameters
  ----------
  ele : EleStruct
      RF Element being tracked through.

  ds : float, optional
      Distance of particle from start edge. Default is zero.

  time : float
      Offset time.
  )"""
  );
  m.def(
      "rfun",
      &Bmad::rfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("a"),
      py::arg("b"),
      py::arg("hz"),
      py::arg("i"),
      py::arg("j"),
      py::arg("res"),
      R"""(Wrapper for Fortran routine rfun

  Parameters
  ----------
  u : 
  v : 
  w : 
  gam : 
  a : 
  b : 
  hz : 
  i : 
  j : 
  res : 
  )"""
  );
  m.def(
      "rk_adaptive_time_step",
      &Bmad::rk_adaptive_time_step,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orb"),
      py::arg("t_dir"),
      py::arg("rf_time"),
      py::arg("dt_try"),
      py::arg("dt_did"),
      py::arg("dt_next"),
      py::arg("err_flag"),
      py::arg("extra_field") = py::none(),
      R"""(Wrapper for Fortran routine rk_adaptive_time_step

  Parameters
  ----------
  ele : 
  param : 
  orb : 
  t_dir : 
  rf_time : 
  dt_try : 
  dt_did : 
  dt_next : 
  err_flag : 
  extra_field : 
  )"""
  );
  m.def(
      "rk_time_step1",
      &Bmad::rk_time_step1,
      py::arg("ele"),
      py::arg("param"),
      py::arg("rf_time"),
      py::arg("orb"),
      py::arg("dt"),
      py::arg("new_orb"),
      py::arg("dr_dt") = py::none(),
      py::arg("err_flag"),
      py::arg("print_err") = py::none(),
      py::arg("extra_field") = py::none(),
      R"""(Wrapper for Fortran routine rk_time_step1

  Parameters
  ----------
  ele : 
  param : 
  rf_time : 
  orb : 
  dt : 
  new_orb : 
  r_err : 
  dr_dt : 
  err_flag : 
  print_err : 
  extra_field : 
  )"""
  );
  m.def(
      "rotate3",
      &Bmad::rotate3,
      py::arg("vec"),
      py::arg("angle"),
      py::arg("rvec"),
      R"""(Wrapper for Fortran routine rotate3

  Parameters
  ----------
  vec : 
  angle : 
  rvec : 
  )"""
  );
  m.def(
      "rotate_em_field",
      &Bmad::rotate_em_field,
      py::arg("field"),
      py::arg("w_mat"),
      py::arg("w_inv"),
      py::arg("calc_dfield") = py::none(),
      py::arg("calc_potential") = py::none(),
      R"""(Subroutine rotate_em_field (field, w_mat, w_inv, calc_dfield, calc_potential)

  Routine to transform the fields using the given rotation matrices.

  Parameters
  ----------
  field : EmFieldStruct
      E and B fields and derivatives.

  w_mat : float
      rotation matrix.

  w_inv : float
      rotation matrix inverse = transpose(w_mat)

  calc_dfield : bool, optional
      If present and True then rotate the field derivatives.

  calc_potential : bool, optional
      Rotate the magnetic vector potential? Default is false.
  )"""
  );
  m.def(
      "rotate_field_zx",
      &Bmad::rotate_field_zx,
      py::arg("field"),
      py::arg("theta"),
      R"""(Wrapper for Fortran routine rotate_field_zx

  Parameters
  ----------
  field : 
  theta : 
  )"""
  );
  m.def(
      "rotate_for_curved_surface",
      &Bmad::rotate_for_curved_surface,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("set"),
      py::arg("rot_mat"),
      R"""(Wrapper for Fortran routine rotate_for_curved_surface

  Parameters
  ----------
  ele : EleStruct
      reflecting element

  orbit : CoordStruct
      Photon position.

  set : bool
      True -> Transform body coords to local curved body coords.

  rot_mat : float
      When set = False, rotation matrix calculated from previous call with set = True.
      This parameter is an input/output and is modified in-place. As an output: When set = True, calculated
      rotation matrix.
  )"""
  );
  m.def(
      "rotate_spin",
      &Bmad::rotate_spin,
      py::arg("rot_vec"),
      py::arg("spin"),
      R"""(Wrapper for Fortran routine rotate_spin

  Parameters
  ----------
  rot_vec : float
      Rotation axis. Magnitude of rot_vec is the rotation angle.

  spin : float
      Initial coords.
      This parameter is an input/output and is modified in-place. As an output: Final coords.

  qrot : float
      : rotation quaternion.
  )"""
  );
  m.def(
      "rotate_spin_a_step",
      &Bmad::rotate_spin_a_step,
      py::arg("orbit"),
      py::arg("field"),
      py::arg("ele"),
      py::arg("ds"),
      R"""(Wrapper for Fortran routine rotate_spin_a_step

  Parameters
  ----------
  orbit : CoordStruct
      Initial orbit.
      This parameter is an input/output and is modified in-place. As an output: Orbit with rotated spin

  field : EmFieldStruct
      EM Field

  ele : 
      ele_struct, Element being tracked through.

  ds : float
      Longitudinal step in element body frame.
  )"""
  );
  m.def(
      "rotate_spin_given_field",
      &Bmad::rotate_spin_given_field,
      py::arg("orbit"),
      py::arg("sign_z_vel"),
      py::arg("BL") = py::none(),
      py::arg("EL") = py::none(),
      py::arg("qrot") = py::none(),
      R"""(Wrapper for Fortran routine rotate_spin_given_field

  Parameters
  ----------
  orbit : CoordStruct
      Initial orbit.
      This parameter is an input/output and is modified in-place. As an output: Orbit with rotated spin

  sign_z_vel : int
      +/- 1. Sign of direction of travel relative to the element.

  BL : float, optional
      Integrated field strength. Assumed zero if not present.

  EL : float, optional
      Integrated field strength. Assumed zero if not present.

  qrot : float, optional
      Initial rotation quaternion.
      This parameter is an input/output and is modified in-place. As an output: Rotation quaternion with
      rotation due to the field added in.
  )"""
  );
}
