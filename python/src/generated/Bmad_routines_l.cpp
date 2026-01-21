#include "pybmad/generated/Bmad_routines_l.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyLatEleLocator python_lat_ele_locator(
    std::string loc_str,
    LatStruct &lat,
    ElePointerStructAlloc1D eles,
    int n_loc,
    std::optional<bool> above_ubound_is_err = std::nullopt,
    std::optional<int> ix_dflt_branch = std::nullopt,
    std::optional<bool> order_by_index = std::nullopt,
    std::optional<bool> append_eles = std::nullopt
) {
  auto _result = Bmad::lat_ele_locator(
      loc_str,
      lat,
      eles,
      n_loc,
      above_ubound_is_err,
      ix_dflt_branch,
      order_by_index,
      append_eles
  );
  auto py_result{PyLatEleLocator{_result, n_loc}};
  return py_result;
}

void init_Bmad_routines_l(py::module &m) {
  m.def(
      "lafun",
      &Bmad::lafun,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::arg("res"),
      R"""(Parameters
  ----------
  x : 
  y : 
  z : 
  res : 
  )"""
  );
  m.def(
      "lat_compute_ref_energy_and_time",
      &Bmad::lat_compute_ref_energy_and_time,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Input lattice.
  err_flag : bool
      Set true if there is an error. False otherwise.
  )"""
  );
  py::class_<PyLatEleLocator, std::unique_ptr<PyLatEleLocator>>(
      m,
      "LatEleLocator",
      "lat_ele_locator return type"
  )
      .def_readonly("err", &PyLatEleLocator::err)
      .def_readonly("n_loc", &PyLatEleLocator::n_loc)
      .def("__len__", [](const PyLatEleLocator &) { return 2; })
      .def("__getitem__", [](const PyLatEleLocator &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err);
        if (i == 1)
          return py::cast(s.n_loc);
        throw py::index_error();
      });
  m.def(
      "lat_ele_locator",
      &python_lat_ele_locator,
      py::arg("loc_str"),
      py::arg("lat"),
      py::arg("eles"),
      py::arg("n_loc"),
      py::arg("above_ubound_is_err") = py::none(),
      py::arg("ix_dflt_branch") = py::none(),
      py::arg("order_by_index") = py::none(),
      py::arg("append_eles") = py::none(),
      R"""(Parameters
  ----------
  loc_str : unknown
      Element names or indexes. May be lower case.
  lat : LatStruct
      Lattice to search through.
  eles : ElePointerStruct
      If append_eles is True, save existing elements.
      This parameter is an input/output and is modified in-place. As an output: Array of matching elements.
  n_loc : int
      Number of existing elements. Used if append_eles is True.
      This parameter is an input/output and is modified in-place. As an output: Number of locations found.
  err : bool
      Set True if there is a decode error. Note: Not finding any matching element is not an error.
  above_ubound_is_err : unknown, optional
      Default is True. If the upper bound "e2" on an "e1:e2" range construct is an integer and above the maximum
      element index then treat this as an error?
  ix_dflt_branch : int, optional
      If present and not -1 then restrict search to specified branch. If not present or -1: Search all branches.
      Exception: For elements specified using
  order_by_index : bool, optional
      False is default. If True, order a component of loc_str like "quad::*" by element index instead of
      longitudinal s-position. Index ordering and s-position ordering
  append_eles : bool, optional
      Default is False. If True, found elements are appended to eles(:) array.
  )"""
  );
  m.def(
      "lat_equal_lat",
      &Bmad::lat_equal_lat,
      py::arg("lat_out"),
      py::arg("lat_in"),
      R"""(Parameters
  ----------
  lat_out : 
  lat_in : 
  )"""
  );
  m.def(
      "lat_geometry",
      &Bmad::lat_geometry,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      The lattice.
  )"""
  );
  m.def(
      "lat_make_mat6",
      &Bmad::lat_make_mat6,
      py::arg("lat"),
      py::arg("ix_ele") = py::none(),
      py::arg("ref_orb") = py::none(),
      py::arg("ix_branch") = py::none(),
      R"""(Parameters
  ----------
  lat : LatStruct
      Lat containing the elements.
  ix_ele : int, optional
      Index of the element. If not present
  ref_orb : CoordStruct, optional
      Coordinates of the reference orbit around which the matrix is calculated. If not present
  ix_branch : int, optional
      Branch index. Default is 0 (main lattice). -1 => All branches/all elements (ref_orb & ix_ele will be
      ignored).
  err_flag : bool
      True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "lat_sanity_check",
      &Bmad::lat_sanity_check,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Lattice to check
  err_flag : bool
      Set True if there is an error. False otherwise.
  )"""
  );
  m.def(
      "lat_to_ptc_layout",
      &Bmad::lat_to_ptc_layout,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Input lattice
  )"""
  );
  m.def(
      "lat_vec_equal_lat_vec",
      &Bmad::lat_vec_equal_lat_vec,
      py::arg("lat1"),
      py::arg("lat2"),
      R"""(Parameters
  ----------
  lat1 : 
  lat2 : 
  )"""
  );
  m.def(
      "lattice_bookkeeper",
      &Bmad::lattice_bookkeeper,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Lattice needing bookkeeping.
      This parameter is an input/output and is modified in-place. As an output: Lattice with bookkeeping done.
  err_flag : bool
      Set true if there is an error. False otherwise.
  )"""
  );
  m.def(
      "lcavity_rf_step_setup",
      &Bmad::lcavity_rf_step_setup,
      py::arg("ele"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lcavity element.
      This parameter is an input/output and is modified in-place. As an output: Element with ele.rf properly
      setup.
  )"""
  );
  m.def(
      "linear_bend_edge_kick",
      &Bmad::linear_bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orb"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine linear_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

  Subroutine to track through the edge field of an sbend.
  Apply only the first order kick, which is edge focusing.

  Parameters
  ----------
  ele : EleStruct
      SBend element.
  param : LatParamStruct
      Rel charge.
  particle_at : int
      first_track_edge$, or second_track_edge$,
  orb : CoordStruct
      Starting coords.
      This parameter is an input/output and is modified in-place. As an output: Coords after tracking.
  mat6 : float, optional
      Transfer matrix up to the edge.
      This parameter is an input/output and is modified in-place. As an output: Transfer matrix including the
      edge.
  make_matrix : float, optional
      Propagate the transfer matrix? Default is False.
  )"""
  );
  py::class_<Bmad::LinearCoef, std::unique_ptr<Bmad::LinearCoef>>(
      m,
      "LinearCoef",
      "linear_coef return type"
  )
      .def_readonly("err_flag", &Bmad::LinearCoef::err_flag)
      .def_readonly("coef", &Bmad::LinearCoef::coef)
      .def("__len__", [](const Bmad::LinearCoef &) { return 2; })
      .def("__getitem__", [](const Bmad::LinearCoef &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.coef);
        throw py::index_error();
      });
  m.def(
      "linear_coef",
      &Bmad::linear_coef,
      py::arg("stack"),
      R"""(Function linear_coef (stack, err_flag) result (coef)

  Routine to return the linear coefficient of a linear expression.

  Parameters
  ----------
  stack : ExpressionAtomStruct
      Expression stack.

  Returns
  -------
  err_flag : bool
      Set True if the expression is not linear
  coef : float
      Linear coefficient.
  )"""
  );
  m.def(
      "linear_to_spin_taylor",
      &Bmad::linear_to_spin_taylor,
      py::arg("q_map"),
      R"""(Parameters
  ----------
  q_map : float
      Linear quaternion map.
  spin_taylor : TaylorStruct
      Taylor map
  )"""
  );
  py::class_<Bmad::LoadParseLine, std::unique_ptr<Bmad::LoadParseLine>>(
      m,
      "LoadParseLine",
      "load_parse_line return type"
  )
      .def_readonly("end_of_file", &Bmad::LoadParseLine::end_of_file)
      .def_readonly("err_flag", &Bmad::LoadParseLine::err_flag)
      .def("__len__", [](const Bmad::LoadParseLine &) { return 2; })
      .def("__getitem__", [](const Bmad::LoadParseLine &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.end_of_file);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "load_parse_line",
      &Bmad::load_parse_line,
      py::arg("action"),
      py::arg("ix_start"),
      R"""(Subroutine load_parse_line (action, ix_start, end_of_file, err_flag)

  Subroutine to load characters from the input file.
  This subroutine is used by bmad_parser and bmad_parser2.
  This subroutine is not intended for general use.

  Parameters
  ----------
  action : unknown
      'continue', 'new_command', or 'init'
  ix_start : int
      Index in bp_com.parse_line string where to append stuff.

  Returns
  -------
  end_of_file : bool
      End of file reached?
  err_flag : bool
      Set True if there is an error. False otherwise
  bp_com%parse_line : 
      string to append to.
  )"""
  );
  m.def(
      "lord_edge_aligned",
      &Bmad::lord_edge_aligned,
      py::arg("slave"),
      py::arg("slave_edge"),
      py::arg("lord"),
      R"""(Parameters
  ----------
  slave : EleStruct
      Slave element.
  slave_edge : int
      End under consideration: entrance_end$, exit_end$, in_between$, etc.
  lord : EleStruct
      Lord element.
  is_aligned : int
      True if a lord edge is aligned with the slave edge. If slave_edge is not entrance_end$ nor exit_end$ then
      is_aligned is False.
  )"""
  );
  m.def(
      "low_energy_z_correction",
      &Bmad::low_energy_z_correction,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("ds"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Position before correction
  ele : unknown
      Element being tracked through.
  ds : float
      Longitudinal distance traveled by reference particle.
  mat6 : float, optional
      Transfer matrix before the multipole.
      This parameter is an input/output and is modified in-place. As an output: Transfer matrix transfer matrix
      including multipole.
  make_matrix : bool, optional
      Propagate the transfer matrix? Default is false.
  dz : float
      Change in z.
  )"""
  );
}
