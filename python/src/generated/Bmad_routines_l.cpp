#include "pybmad/generated/Bmad_routines_l.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
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

void init_Bmad_routines_l(nb::module_ &m) {
  m.def(
      "lafun",
      &Bmad::lafun,
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("z"),
      R"""(Wrapper for Fortran routine lafun

Parameters
----------
x : float

y : float

z : float

Returns
-------
res : float
)"""
  );
  m.def(
      "lat_compute_ref_energy_and_time",
      &Bmad::lat_compute_ref_energy_and_time,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine lat_compute_ref_energy_and_time

Parameters
----------
lat : LatStruct
    Input lattice.

Returns
-------
err_flag : bool
    Set true if there is an error. False otherwise.
)"""
  );
  nb::class_<PyLatEleLocator>(m, "LatEleLocator", "lat_ele_locator return type")
      .def_ro("err", &PyLatEleLocator::err)
      .def_ro("n_loc", &PyLatEleLocator::n_loc)
      .def("__len__", [](const PyLatEleLocator &) { return 2; })
      .def("__getitem__", [](const PyLatEleLocator &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.n_loc);
        throw nb::index_error();
      });
  m.def(
      "lat_ele_locator",
      &python_lat_ele_locator,
      nb::arg("loc_str"),
      nb::arg("lat"),
      nb::arg("eles"),
      nb::arg("n_loc"),
      nb::arg("above_ubound_is_err") = nb::none(),
      nb::arg("ix_dflt_branch") = nb::none(),
      nb::arg("order_by_index") = nb::none(),
      nb::arg("append_eles") = nb::none(),
      R"""(Wrapper for Fortran routine lat_ele_locator

Parameters
----------
loc_str : str
    Element names or indexes. May be lower case.

lat : LatStruct
    Lattice to search through.

eles : 1D array of ElePointerStruct
    If append_eles is True, save existing elements.
    This parameter is an input/output and is modified in-place.
    As an output, eles: Array of matching elements.

n_loc : int
    Number of existing elements. Used if append_eles is True.
    This parameter is an input/output and is modified in-place.
    As an output, n_loc: Number of locations found.

above_ubound_is_err : bool, optional
    Default is True. If the upper bound "e2" on an "e1:e2" range construct is an integer and above the maximum
    element index then treat this as an error? If False, treat e2 as the maximum element index.

ix_dflt_branch : int, optional
    If present and not -1 then restrict search to specified branch. If not present or -1: Search all branches.
    Exception: For elements specified using an integer index (EG: "43"), if ix_dflt_branch is not present or
    -1 use branch 0.

order_by_index : bool, optional
    False is default. If True, order a component of loc_str like "quad::*" by element index instead of
    longitudinal s-position. Index ordering and s-position ordering are different when there are super lords
    and super slaves.

append_eles : bool, optional
    Default is False. If True, found elements are appended to eles(:) array.

Returns
-------
n_loc : int
    Number of existing elements. Used if append_eles is True.
    This parameter is an input/output and is modified in-place.
    As an output, n_loc: Number of locations found.

err : bool, optional
    Set True if there is a decode error. Note: Not finding any matching element is not an error.
)"""
  );
  m.def(
      "lat_equal_lat",
      &Bmad::lat_equal_lat,
      nb::arg("lat_out"),
      nb::arg("lat_in"),
      R"""(Wrapper for Fortran routine lat_equal_lat

Parameters
----------
lat_out : LatStruct

lat_in : LatStruct
)"""
  );
  m.def(
      "lat_geometry",
      &Bmad::lat_geometry,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine lat_geometry

Parameters
----------
lat : LatStruct
    The lattice.
)"""
  );
  m.def(
      "lat_make_mat6",
      &Bmad::lat_make_mat6,
      nb::arg("lat"),
      nb::arg("ix_ele") = nb::none(),
      nb::arg("ref_orb") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      R"""(Wrapper for Fortran routine lat_make_mat6

Parameters
----------
lat : LatStruct
    Lat containing the elements.

ix_ele : int, optional
    Index of the element. If not present or negative, the matrices for all elements will be calculated.

ref_orb : 1D array of CoordStruct, optional
    Coordinates of the reference orbit around which the matrix is calculated. If not present then the
    referemce is taken to be the origin.

ix_branch : int, optional
    Branch index. Default is 0 (main lattice). -1 => All branches/all elements (ref_orb & ix_ele will be
    ignored).

Returns
-------
err_flag : bool, optional
    True if there is an error. False otherwise.
)"""
  );
  m.def(
      "lat_sanity_check",
      &Bmad::lat_sanity_check,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine lat_sanity_check

Parameters
----------
lat : LatStruct
    Lattice to check

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "lat_to_ptc_layout",
      &Bmad::lat_to_ptc_layout,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine lat_to_ptc_layout

Parameters
----------
lat : LatStruct
    Input lattice
)"""
  );
  m.def(
      "lat_vec_equal_lat_vec",
      &Bmad::lat_vec_equal_lat_vec,
      nb::arg("lat1"),
      nb::arg("lat2"),
      R"""(Wrapper for Fortran routine lat_vec_equal_lat_vec

Parameters
----------
lat1 : 1D array of LatStruct

lat2 : 1D array of LatStruct
)"""
  );
  m.def(
      "lattice_bookkeeper",
      &Bmad::lattice_bookkeeper,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine lattice_bookkeeper

Parameters
----------
lat : LatStruct
    Lattice needing bookkeeping.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with bookkeeping done.

Returns
-------
err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  m.def(
      "lcavity_rf_step_setup",
      &Bmad::lcavity_rf_step_setup,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine lcavity_rf_step_setup

Parameters
----------
ele : EleStruct
    Lcavity element.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with ele.rf properly setup.
)"""
  );
  m.def(
      "linear_bend_edge_kick",
      &Bmad::linear_bend_edge_kick,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("particle_at"),
      nb::arg("orb"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Subroutine to track through the edge field of an sbend.
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
    This parameter is an input/output and is modified in-place.
    As an output, orb: Coords after tracking.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the edge.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix including the edge.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.
)"""
  );
  nb::class_<Bmad::LinearCoef>(m, "LinearCoef", "linear_coef return type")
      .def_ro("err_flag", &Bmad::LinearCoef::err_flag)
      .def_ro("coef", &Bmad::LinearCoef::coef)
      .def("__len__", [](const Bmad::LinearCoef &) { return 2; })
      .def("__getitem__", [](const Bmad::LinearCoef &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.coef);
        throw nb::index_error();
      });
  m.def(
      "linear_coef",
      &Bmad::linear_coef,
      nb::arg("stack"),
      R"""(Routine to return the linear coefficient of a linear expression.

Parameters
----------
stack : 1D array of ExpressionAtomStruct
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
      nb::arg("q_map"),
      R"""(Wrapper for Fortran routine linear_to_spin_taylor

Parameters
----------
q_map : 2D array of float (shape: 0:3, 0:6)
    Linear quaternion map.

Returns
-------
spin_taylor : 1D array of TaylorStruct (shape: 0:3)
    Taylor map
)"""
  );
  nb::class_<Bmad::LoadParseLine>(m, "LoadParseLine", "load_parse_line return type")
      .def_ro("end_of_file", &Bmad::LoadParseLine::end_of_file)
      .def_ro("err_flag", &Bmad::LoadParseLine::err_flag)
      .def("__len__", [](const Bmad::LoadParseLine &) { return 2; })
      .def("__getitem__", [](const Bmad::LoadParseLine &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.end_of_file);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "load_parse_line",
      &Bmad::load_parse_line,
      nb::arg("action"),
      nb::arg("ix_start"),
      R"""(Subroutine to load characters from the input file.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
action : str
    'continue', 'new_command', or 'init'

ix_start : int
    Index in bp_com.parse_line string where to append stuff.

Returns
-------
end_of_file : bool
    End of file reached?

err_flag : bool, optional
    Set True if there is an error. False otherwise
)"""
  );
  m.def(
      "lord_edge_aligned",
      &Bmad::lord_edge_aligned,
      nb::arg("slave"),
      nb::arg("slave_edge"),
      nb::arg("lord"),
      R"""(Wrapper for Fortran routine lord_edge_aligned

Parameters
----------
slave : EleStruct
    Slave element.

slave_edge : int
    End under consideration: entrance_end$, exit_end$, in_between$, etc.

lord : EleStruct
    Lord element.

Returns
-------
is_aligned : bool
    True if a lord edge is aligned with the slave edge. If slave_edge is not entrance_end$ nor exit_end$ then
    is_aligned is False.
)"""
  );
  m.def(
      "low_energy_z_correction",
      &Bmad::low_energy_z_correction,
      nb::arg("orbit"),
      nb::arg("ele"),
      nb::arg("ds"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Wrapper for Fortran routine low_energy_z_correction

Parameters
----------
orbit : CoordStruct
    Position before correction

ele : EleStruct
    Element being tracked through.

ds : float
    Longitudinal distance traveled by reference particle.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the multipole.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix transfer matrix including multipole.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.

Returns
-------
dz : float
    Change in z.
)"""
  );
}
