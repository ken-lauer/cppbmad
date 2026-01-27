#include "pybmad/generated/Bmad_routines_e.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyExpectOneOf python_expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool is_ok
) {
  Bmad::expect_one_of(delim_list, check_input_delim, ele_name, delim, delim_found, is_ok);
  auto py_result{PyExpectOneOf{delim}};
  return py_result;
}

void init_Bmad_routines_e(py::module &m) {
  m.def(
      "e_accel_field",
      &Bmad::e_accel_field,
      py::arg("ele"),
      py::arg("voltage_or_gradient"),
      py::arg("bmad_standard_tracking") = py::none(),
      R"""(Wrapper for Fortran routine e_accel_field

Parameters
----------
ele : EleStruct
    Lcavity or rfcavity element.

voltage_or_gradient : int
    voltage$ or gradient$

bmad_standard_tracking : bool, optional
    Using bmad_standard tracking? Default is False.

Returns
-------
field : float
    Cavity field or gradient.
)"""
  );
  m.def(
      "e_crit_photon",
      &Bmad::e_crit_photon,
      py::arg("gamma"),
      py::arg("g_bend"),
      R"""(Function E_crit_photon (gamma, g_bend) result (E_crit)

Routine to calculate the photon critical energy in a bend.

Parameters
----------
gamma : float
    Gamma factor of charged particle emitting photon.

g_bend : float
    1/radius bending strength.

Returns
-------
E_crit : float
    Critical photon energy.
)"""
  );
  py::class_<Bmad::EigenDecomp6mat, std::unique_ptr<Bmad::EigenDecomp6mat>>(
      m,
      "EigenDecomp6mat",
      "eigen_decomp_6mat return type"
  )
      .def_readonly("eval", &Bmad::EigenDecomp6mat::eval)
      .def_readonly("evec", &Bmad::EigenDecomp6mat::evec)
      .def_readonly("err_flag", &Bmad::EigenDecomp6mat::err_flag)
      .def_readonly("tunes", &Bmad::EigenDecomp6mat::tunes)
      .def("__len__", [](const Bmad::EigenDecomp6mat &) { return 4; })
      .def("__getitem__", [](const Bmad::EigenDecomp6mat &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.eval);
        if (i == 1)
          return py::cast(s.evec);
        if (i == 2)
          return py::cast(s.err_flag);
        if (i == 3)
          return py::cast(s.tunes);
        throw py::index_error();
      });
  m.def(
      "eigen_decomp_6mat",
      &Bmad::eigen_decomp_6mat,
      py::arg("mat"),
      R"""(Subroutine eigen_decomp_6mat(mat, eval, evec, tunes, err_flag)

Compute eigenvalues and eigenvectors of a real 6x6 matrix.
The evals and evecs are in general complex.

Parameters
----------
mat : 2D array of float (shape: 6,6)
    6x6 real matrix.  Usually a transfer matrix or sigma matrix.

Returns
-------
eval : 1D array of complex (shape: 6)
    complex eigenvalues.

evec : 2D array of complex (shape: 6,6)
    complex eigenvectors arranged down columns.

err_flag : bool
    set to true if an error has occured.

tunes : 1D array of float (shape: 3), optional
    Mode tunes, in radians.
)"""
  );
  m.def(
      "ele_compute_ref_energy_and_time",
      &Bmad::ele_compute_ref_energy_and_time,
      py::arg("ele0"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("err_flag"),
      R"""(Wrapper for Fortran routine ele_compute_ref_energy_and_time

Parameters
----------
ele0 : EleStruct
    Previous element in lattice with starting energy and time values.

ele : EleStruct
    Lattice element
    This parameter is an input/output and is modified in-place.
    As an output, ele: Lattice element with reference energy and time.

param : LatParamStruct
    Lattice parameters.

err_flag : bool
    Set true if there is an error. False otherwise.
)"""
  );
  m.def(
      "ele_equal_ele",
      &Bmad::ele_equal_ele,
      py::arg("ele_out"),
      py::arg("ele_in"),
      R"""(Wrapper for Fortran routine ele_equal_ele

Parameters
----------
ele_out : EleStruct

ele_in : EleStruct
)"""
  );
  m.def(
      "ele_equals_ele",
      &Bmad::ele_equals_ele,
      py::arg("ele_in"),
      py::arg("update_nametable"),
      R"""(Subroutine ele_equals_ele (ele_out, ele_in, update_nametable)

Subroutine that is used to set an element equal to another.
Note: Use ele_equal_ele instead unless you know what you are doing.

Parameters
----------
ele_in : EleStruct
    Input element.

update_nametable : bool
    If true, update the nametable. If false, do not. Note: nametable updates can take time if this routine is
    called a many times. See remove_eles_from_lat as an example.

Returns
-------
ele_out : EleStruct
    Output element.
)"""
  );
  m.def(
      "ele_finalizer",
      &Bmad::ele_finalizer,
      py::arg("ele"),
      R"""(Subroutine ele_finalizer(ele)

Finalizer routine for ele_struct instances.
NOTE: Not currently used.

Parameters
----------
ele : EleStruct
    Element to cleanup.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with pointers deallocated as needed.
)"""
  );
  m.def(
      "ele_full_name",
      &Bmad::ele_full_name,
      py::arg("ele"),
      py::arg("template_") = py::none(),
      R"""(Wrapper for Fortran routine ele_full_name

Parameters
----------
ele : EleStruct
    Element in a lattice

Returns
-------
str : character
    : Name/location string.
)"""
  );
  m.def(
      "ele_geometry",
      &Bmad::ele_geometry,
      py::arg("floor_start"),
      py::arg("ele"),
      py::arg("len_scale") = py::none(),
      py::arg("ignore_patch_err") = py::none(),
      R"""(Wrapper for Fortran routine ele_geometry

Parameters
----------
floor_start : FloorPositionStruct
    Starting floor coordinates at upstream end. Not used for fiducial and girder elements.

ele : EleStruct
    Element to propagate the geometry through.

len_scale : float, optional
    factor to scale the length of the element. 1.0_rp => Output is geometry at end of element (default).
    0.5_rp => Output is geometry at center of element. -1.0_rp => Used to propagate geometry in reverse.

ignore_patch_err : bool, optional
    If present and True, ignore flexible patch errors. This is used by ele_compute_ref_energy_and_time to
    suppress unnecessary messages.

Returns
-------
floor_end : FloorPositionStruct, optional
    Output floor position. If not present then ele.floor will be used and ele.bookkeeping_state.floor_position
    will be set to ok$.
)"""
  );
  m.def(
      "ele_geometry_with_misalignments",
      &Bmad::ele_geometry_with_misalignments,
      py::arg("ele"),
      py::arg("len_scale") = py::none(),
      R"""(Wrapper for Fortran routine ele_geometry_with_misalignments

Parameters
----------
ele : EleStruct
    Lattice element under consideration.

len_scale : float, optional
    factor to scale the length of the element. 1.0_rp => Output is geometry at end of element (default).
    0.5_rp => Output is geometry at center of element. -1.0_rp => Used to propagate geometry in reverse.

Returns
-------
floor : FloorPositionStruct
    Floor position with misalignments
)"""
  );
  m.def(
      "ele_has_constant_ds_dt_ref",
      &Bmad::ele_has_constant_ds_dt_ref,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine ele_has_constant_ds_dt_ref

Parameters
----------
ele : EleStruct
    Element.

Returns
-------
is_const : bool
    True if reference velocity must be a constant.
)"""
  );
  m.def(
      "ele_has_nonzero_kick",
      &Bmad::ele_has_nonzero_kick,
      py::arg("ele"),
      py::arg("has_kick"),
      R"""(Wrapper for Fortran routine ele_has_nonzero_kick

Parameters
----------
ele : EleStruct
    Element with possible nonzero kicks.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with no kicks.

has_kick : bool
)"""
  );
  m.def(
      "ele_has_nonzero_offset",
      &Bmad::ele_has_nonzero_offset,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine ele_has_nonzero_offset

Parameters
----------
ele : EleStruct
    Element with possible nonzero offsets.

Returns
-------
has_offset : bool
    Set true is element has a non-zero offset.
)"""
  );
  m.def(
      "ele_is_monitor",
      &Bmad::ele_is_monitor,
      py::arg("ele"),
      py::arg("print_warning") = py::none(),
      R"""(Function ele_is_monitor (ele, print_warning) result (is_monitor)

Routine to check that an element is either a detector, instrument, monitor, or marker.
These are the elements where measurement errors can be defined.

Parameters
----------
ele : EleStruct
    Lattice element.

print_warning : bool, optional
    If True print a warning message if the element not a monitor like element. Default is True.

Returns
-------
is_monitor : bool
    Set True if the element is a monitor like element.
)"""
  );
  m.def(
      "ele_loc",
      &Bmad::ele_loc,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine ele_loc

Parameters
----------
ele : EleStruct
    Element to be identified

Returns
-------
loc : LatEleLocStruct
    Element identifier.
)"""
  );
  m.def(
      "ele_loc_name",
      &Bmad::ele_loc_name,
      py::arg("ele"),
      py::arg("show_branch0") = py::none(),
      py::arg("parens") = py::none(),
      R"""(Wrapper for Fortran routine ele_loc_name

Parameters
----------
ele : EleStruct
    Element in a lattice

show_branch0 : bool, optional
    Explicitly show branch for main lattice elements? Default is False.

parens : character, optional
    If present, enclose location string using the two characters supplied. Typically parens will be set to
    "()" or "[]".

Returns
-------
str : character
    Output string. Left justified.
)"""
  );
  py::class_<Bmad::EleMisalignmentLSCalc, std::unique_ptr<Bmad::EleMisalignmentLSCalc>>(
      m,
      "EleMisalignmentLSCalc",
      "ele_misalignment_l_s_calc return type"
  )
      .def_readonly("L_mis", &Bmad::EleMisalignmentLSCalc::L_mis)
      .def_readonly("S_mis", &Bmad::EleMisalignmentLSCalc::S_mis)
      .def("__len__", [](const Bmad::EleMisalignmentLSCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::EleMisalignmentLSCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.L_mis);
        if (i == 1)
          return py::cast(s.S_mis);
        throw py::index_error();
      });
  m.def(
      "ele_misalignment_l_s_calc",
      &Bmad::ele_misalignment_l_s_calc,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine ele_misalignment_l_s_calc

Parameters
----------
ele : EleStruct
    Element

Returns
-------
L_mis : 1D array of float (shape: 3)
    Misalignment vector relative to center of element

S_mis : 2D array of float (shape: 3,3)
    Misalignment matrix relative to center of element
)"""
  );
  m.def(
      "ele_nametable_index",
      &Bmad::ele_nametable_index,
      py::arg("ele"),
      R"""(Wrapper for Fortran routine ele_nametable_index

Parameters
----------
ele : EleStruct
    Element in a lattice.

Returns
-------
ix_nt : int
    Nametable index. lat.nametable.name(ix_nt) and lat.nametable.index(ix_nt) correspond with ele. Set to -1
    if ele is not a lattice element. For example, a slice_slave is not a lattice element.
)"""
  );
  m.def(
      "ele_order_calc",
      &Bmad::ele_order_calc,
      py::arg("lat"),
      R"""(Wrapper for Fortran routine ele_order_calc

Parameters
----------
lat : LatStruct
    Lattice to analyze.

Returns
-------
order : LatEleOrderStruct
    Structure holding the element order information.
)"""
  );
  m.def(
      "ele_reference_energy_correction",
      &Bmad::ele_reference_energy_correction,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("particle_at"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Wrapper for Fortran routine ele_reference_energy_correction

Parameters
----------
ele : EleStruct
    Element being tracked through.

orbit : CoordStruct
    Coordinates to correct.

particle_at : int
    first_track_edge$ (that is, entering the element), or second_track_edge$ (that is, leaving the element),
    or upstream_end$ (inherit ele.value(p0c_start$) ref), or downstream_end$ (inherit ele.value(p0c$)).
    inside$ (or anything else) -> Do nothing.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before correction.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix transfer matrix including correction.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "ele_rf_step_index",
      &Bmad::ele_rf_step_index,
      py::arg("E_ref"),
      py::arg("s_rel"),
      py::arg("ele"),
      R"""(Wrapper for Fortran routine ele_rf_step_index

Parameters
----------
E_ref : float
    Reference energy of step. If negative, ignore and use s_rel.

s_rel : float
    S-position relative to the beginning of the element

ele : EleStruct
    RF cavity.

Returns
-------
ix_step : int
    Corresponding index in the ele.rf.steps(:) array.
)"""
  );
  py::class_<Bmad::EleToFibre, std::unique_ptr<Bmad::EleToFibre>>(
      m,
      "EleToFibre",
      "ele_to_fibre return type"
  )
      .def_readonly("ptc_fibre", &Bmad::EleToFibre::ptc_fibre)
      .def_readonly("err_flag", &Bmad::EleToFibre::err_flag)
      .def("__len__", [](const Bmad::EleToFibre &) { return 2; })
      .def("__getitem__", [](const Bmad::EleToFibre &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ptc_fibre);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "ele_to_fibre",
      &Bmad::ele_to_fibre,
      py::arg("ele"),
      py::arg("use_offsets"),
      py::arg("integ_order") = py::none(),
      py::arg("steps") = py::none(),
      py::arg("for_layout") = py::none(),
      py::arg("ref_in") = py::none(),
      R"""(Wrapper for Fortran routine ele_to_fibre

Parameters
----------
ele : EleStruct
    Bmad element.

use_offsets : bool
    Does ptc_fibre include element offsets, pitches and tilt?

integ_order : int, optional
    Order for the sympletic integrator. Possibilities are: 2, 4, or 6 Overrides ele.value(integrator_order$).
    default = 2 (if not set with set_ptc).

steps : int, optional
    Number of integration steps. Overrides ele.value(ds_step$).

for_layout : bool, optional
    If True then fibre will be put in the PTC layout. Default is False.

ref_in : CoordStruct, optional
    Particle to be tracked. ref_particle$, electron$, etc. This argument should only be present when the fibre
    is not to be put in a layout.

Returns
-------
ptc_fibre : Fibre, optional
    PTC fibre element.

err_flag : bool
    Set True if setup OK. False otherwise.
)"""
  );
  m.def(
      "ele_to_ptc_magnetic_bn_an",
      &Bmad::ele_to_ptc_magnetic_bn_an,
      py::arg("ele"),
      py::arg("bn"),
      py::arg("an"),
      R"""(Subroutine ele_to_ptc_magnetic_bn_an (ele, bn, an, n_max)

Routine to compute the a(n) and b(n) magnetic multipole components of a magnet.
This is used to interface between eles and PTC fibres

Note: The multipole index uses the PTC convention of starting from 1 instead of zero.

Note: On the PTC side bn(1) is error field when creating a fibre but
is the total field when the fibre is being modified. This routine returns the error field.

Parameters
----------
ele : EleStruct
    Bmad Element.

bn : 1D array of float
    Normal multipole component.

an : 1D array of float
    Skew multipole component.

Returns
-------
n_max : int, optional
    Maximum non-zero multipole component. Set to zero if there are no multipoles.
)"""
  );
  m.def(
      "ele_to_spin_taylor",
      &Bmad::ele_to_spin_taylor,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orb0"),
      R"""(Wrapper for Fortran routine ele_to_spin_taylor

Parameters
----------
ele : EleStruct
    Lattice element.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with spin map.

param : LatParamStruct
    Branch parameters.

orb0 : CoordStruct
    Starting ref coords.
)"""
  );
  py::class_<Bmad::EleToTaylor, std::unique_ptr<Bmad::EleToTaylor>>(
      m,
      "EleToTaylor",
      "ele_to_taylor return type"
  )
      .def_readonly("orbital_taylor", &Bmad::EleToTaylor::orbital_taylor)
      .def_readonly("spin_taylor", &Bmad::EleToTaylor::spin_taylor)
      .def("__len__", [](const Bmad::EleToTaylor &) { return 2; })
      .def("__getitem__", [](const Bmad::EleToTaylor &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.orbital_taylor);
        if (i == 1)
          return py::cast(s.spin_taylor);
        throw py::index_error();
      });
  m.def(
      "ele_to_taylor",
      &Bmad::ele_to_taylor,
      py::arg("ele"),
      py::arg("orb0") = py::none(),
      py::arg("taylor_map_includes_offsets") = py::none(),
      py::arg("include_damping") = py::none(),
      R"""(Wrapper for Fortran routine ele_to_taylor

Parameters
----------
ele : EleStruct
    Element to construct map for.

orb0 : CoordStruct, optional
    Starting coords around which the Taylor map is evaluated. Default is the zero orbit.

taylor_map_includes_offsets : bool, optional
    If present then value overrides ele.taylor_map_includes_offsets.

include_damping : bool, optional
    Sets if radiation damping is included. Default is what is set in ptc_private.base_state.

Returns
-------
orbital_taylor : 1D array of TaylorStruct (shape: 6), optional
    Orbital taylor map. If not present then the map is put in ele.taylor.

spin_taylor : 1D array of TaylorStruct (shape: 0:3), optional
    Spin taylor map. If not present then the map is put in ele.spin_taylor.
)"""
  );
  m.def(
      "ele_unique_name",
      &Bmad::ele_unique_name,
      py::arg("ele"),
      py::arg("order"),
      R"""(Wrapper for Fortran routine ele_unique_name

Parameters
----------
ele : EleStruct
    Element to construct a unique name for.

order : LatEleOrderStruct
    Information on element ordering. Before calling this routine, use the routine ele_order_calc to compute
    this argument.

Returns
-------
unique_name : character
    Unique name that can can be used to identify ele. The simplist name will be constructed. For example, if
    the element name is unique, unique_name will be set to the element name.
)"""
  );
  m.def(
      "ele_value_has_changed",
      &Bmad::ele_value_has_changed,
      py::arg("ele"),
      py::arg("list"),
      py::arg("abs_tol"),
      py::arg("set_old"),
      R"""(Wrapper for Fortran routine ele_value_has_changed

Parameters
----------
ele : EleStruct
    Element under consideration.
    This parameter is an input/output and is modified in-place.
    As an output, ele: ele.old_value may be set depending upon setting of set_old

list : 1D array of int
    List of indexes of ele.value(:) array to check.

abs_tol : 1D array of float
    List of values such that if the change in parameter value is less than this it is not considered to have
    changed significantly.

set_old : bool
    If True then set ele.old_value(j) = ele.value(j) for j in list

Returns
-------
has_changed : bool
    Set True if a value has changed significantly.
)"""
  );
  m.def(
      "ele_vec_equal_ele_vec",
      &Bmad::ele_vec_equal_ele_vec,
      py::arg("ele1"),
      py::arg("ele2"),
      R"""(Wrapper for Fortran routine ele_vec_equal_ele_vec

Parameters
----------
ele1 : 1D array of EleStruct

ele2 : 1D array of EleStruct
)"""
  );
  py::class_<Bmad::ElecMultipoleField, std::unique_ptr<Bmad::ElecMultipoleField>>(
      m,
      "ElecMultipoleField",
      "elec_multipole_field return type"
  )
      .def_readonly("Ex", &Bmad::ElecMultipoleField::Ex)
      .def_readonly("Ey", &Bmad::ElecMultipoleField::Ey)
      .def_readonly("dE", &Bmad::ElecMultipoleField::dE)
      .def_readonly("compute_dE", &Bmad::ElecMultipoleField::compute_dE)
      .def("__len__", [](const Bmad::ElecMultipoleField &) { return 4; })
      .def("__getitem__", [](const Bmad::ElecMultipoleField &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.Ex);
        if (i == 1)
          return py::cast(s.Ey);
        if (i == 2)
          return py::cast(s.dE);
        if (i == 3)
          return py::cast(s.compute_dE);
        throw py::index_error();
      });
  m.def(
      "elec_multipole_field",
      &Bmad::elec_multipole_field,
      py::arg("a"),
      py::arg("b"),
      py::arg("n"),
      py::arg("coord"),
      R"""(Wrapper for Fortran routine elec_multipole_field

Parameters
----------
a : float
    Multipole skew component.

b : float
    Multipole normal component.

n : int
    Multipole order.

coord : CoordStruct

Returns
-------
Ex : float
    X field component

Ey : float
    Y field component.

dE : 2D array of float (shape: 2,2), optional
    Field derivatives: dfield(x,y)/d(x,y).

compute_dE : bool, optional
    If False, do not compute the field derivatives even if dE is present. Default is True.
)"""
  );
  py::class_<Bmad::ElementAtSBranch, std::unique_ptr<Bmad::ElementAtSBranch>>(
      m,
      "ElementAtSBranch",
      "element_at_s_branch return type"
  )
      .def_readonly("err_flag", &Bmad::ElementAtSBranch::err_flag)
      .def_readonly("s_eff", &Bmad::ElementAtSBranch::s_eff)
      .def_readonly("position", &Bmad::ElementAtSBranch::position)
      .def_readonly("ix_ele", &Bmad::ElementAtSBranch::ix_ele)
      .def("__len__", [](const Bmad::ElementAtSBranch &) { return 4; })
      .def("__getitem__", [](const Bmad::ElementAtSBranch &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.s_eff);
        if (i == 2)
          return py::cast(s.position);
        if (i == 3)
          return py::cast(s.ix_ele);
        throw py::index_error();
      });
  m.def(
      "element_at_s",
      py::overload_cast<BranchStruct &, double, bool, std::optional<bool>>(&Bmad::element_at_s),
      py::arg("branch"),
      py::arg("s"),
      py::arg("choose_max"),
      py::arg("print_err") = py::none(),
      R"""(Function element_at_s (...) result (ix_ele)

Function to return the index of the element at position s.

element_at_s is an overloaded name for:
  function element_at_s_lat (lat, s, choose_max, ix_branch, err_flag, s_eff, position, print_err) result (ix_ele)
  function element_at_s_branch (branch, s, choose_max, err_flag, s_eff, position, print_err) result (ix_ele)

The differnce between these two routine is that with element_at_s_lat, the branch is given by the lat
  and ix_ele arguments: branch = lat%branch(ix_ele). With element_at_s_branch, the branch is an argument.

Also see: pointer_to_element_at_s

ix_ele is choisen such that:
If choose_max = True:
    If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
    Else: branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
If choose_max = False:
    If s = branch%ele(0)%s: ix_ele = 0
    Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
    choose_max = True  => ix_ele = ix2
    choose_max = False => ix_ele = ix1

The setting of choose_max only makes a difference when s corresponds to an element boundary.

Note: For a circular lattice, s is evaluated at the effective s which
is modulo the branch length:
    s_eff = s - branch_length * floor(s/branch_length)

Note: If there are multiple elements that are at the given s position due to the presence of
an element with a negative length, which of the possible elements is actually chosen is ill-defined.

Parameters
----------
branch : BranchStruct
    Branch to use

s : float
    Longitudinal position.

choose_max : bool
    See above

print_err : bool, optional
    Print error message if there is an error? Default is True.

Returns
-------
err_flag : bool, optional
    Set True if s is out of bounds. False otherwise.

s_eff : float, optional
    Effective s. Equal to s with a open lattice. See above.

position : CoordStruct, optional
    Positional information.

ix_ele : int
    Index of element at s.
)"""
  );
  py::class_<Bmad::ElementAtSLat, std::unique_ptr<Bmad::ElementAtSLat>>(
      m,
      "ElementAtSLat",
      "element_at_s_lat return type"
  )
      .def_readonly("err_flag", &Bmad::ElementAtSLat::err_flag)
      .def_readonly("s_eff", &Bmad::ElementAtSLat::s_eff)
      .def_readonly("position", &Bmad::ElementAtSLat::position)
      .def_readonly("ix_ele", &Bmad::ElementAtSLat::ix_ele)
      .def("__len__", [](const Bmad::ElementAtSLat &) { return 4; })
      .def("__getitem__", [](const Bmad::ElementAtSLat &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.s_eff);
        if (i == 2)
          return py::cast(s.position);
        if (i == 3)
          return py::cast(s.ix_ele);
        throw py::index_error();
      });
  m.def(
      "element_at_s",
      py::overload_cast<LatStruct &, double, bool, std::optional<int>, std::optional<bool>>(
          &Bmad::element_at_s
      ),
      py::arg("lat"),
      py::arg("s"),
      py::arg("choose_max"),
      py::arg("ix_branch") = py::none(),
      py::arg("print_err") = py::none(),
      R"""(Function element_at_s (...) result (ix_ele)

Function to return the index of the element at position s.

element_at_s is an overloaded name for:
  function element_at_s_lat (lat, s, choose_max, ix_branch, err_flag, s_eff, position, print_err) result (ix_ele)
  function element_at_s_branch (branch, s, choose_max, err_flag, s_eff, position, print_err) result (ix_ele)

The differnce between these two routine is that with element_at_s_lat, the branch is given by the lat
  and ix_ele arguments: branch = lat%branch(ix_ele). With element_at_s_branch, the branch is an argument.

Also see: pointer_to_element_at_s

ix_ele is choisen such that:
If choose_max = True:
    If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
    Else: branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
If choose_max = False:
    If s = branch%ele(0)%s: ix_ele = 0
    Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
    choose_max = True  => ix_ele = ix2
    choose_max = False => ix_ele = ix1

The setting of choose_max only makes a difference when s corresponds to an element boundary.

Note: For a circular lattice, s is evaluated at the effective s which
is modulo the branch length:
    s_eff = s - branch_length * floor(s/branch_length)

Note: If there are multiple elements that are at the given s position due to the presence of
an element with a negative length, which of the possible elements is actually chosen is ill-defined.

Parameters
----------
lat : LatStruct
    Lattice of elements.

s : float
    Longitudinal position.

choose_max : bool
    See above

ix_branch : int, optional
    Branch index. Default is 0.

print_err : bool, optional
    Print error message if there is an error? Default is True.

Returns
-------
err_flag : bool, optional
    Set True if s is out of bounds. False otherwise.

s_eff : float, optional
    Effective s. Equal to s with a open lattice. See above.

position : CoordStruct, optional
    Positional information.

ix_ele : int
    Index of element at s.
)"""
  );
  m.def(
      "element_slice_iterator",
      &Bmad::element_slice_iterator,
      py::arg("ele"),
      py::arg("param"),
      py::arg("i_slice"),
      py::arg("n_slice_tot"),
      py::arg("sliced_ele"),
      py::arg("s_start") = py::none(),
      py::arg("s_end") = py::none(),
      R"""(Wrapper for Fortran routine element_slice_iterator

Parameters
----------
ele : EleStruct
    Element to slice and dice.

param : LatParamStruct
    Lattice parameters

i_slice : int
    Slice index

n_slice_tot : int
    Total number of slices.

sliced_ele : EleStruct

s_start : float, optional
    Starting edge of slice relative to beginning of element.

s_end : float, optional
    Ending edge of slice relative to beginning of element.
)"""
  );
  m.def(
      "ellipinc_test",
      &Bmad::ellipinc_test,
      R"""(Wrapper for Fortran routine ellipinc_test
)"""
  );
  py::class_<Bmad::EmFieldCalc, std::unique_ptr<Bmad::EmFieldCalc>>(
      m,
      "EmFieldCalc",
      "em_field_calc return type"
  )
      .def_readonly("field", &Bmad::EmFieldCalc::field)
      .def_readonly("err_flag", &Bmad::EmFieldCalc::err_flag)
      .def("__len__", [](const Bmad::EmFieldCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::EmFieldCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.field);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "em_field_calc",
      &Bmad::em_field_calc,
      py::arg("ele"),
      py::arg("param"),
      py::arg("s_pos"),
      py::arg("orbit"),
      py::arg("local_ref_frame"),
      py::arg("calc_dfield") = py::none(),
      py::arg("calc_potential") = py::none(),
      py::arg("use_overlap") = py::none(),
      py::arg("grid_allow_s_out_of_bounds") = py::none(),
      py::arg("rf_time") = py::none(),
      py::arg("used_eles") = py::none(),
      py::arg("print_err") = py::none(),
      py::arg("original_ele") = py::none(),
      R"""(Wrapper for Fortran routine em_field_calc

Parameters
----------
ele : EleStruct
    Lattice element.

param : LatParamStruct
    Lattice parameters.

s_pos : float
    Longitudinal position. If local_ref_frame = T: In Body coords relative to the entrance edge of the
    element. If local_ref_frame = F: In Lab coords relative to the upstream edge of the element.

orbit : CoordStruct
    Transverse coordinates.

local_ref_frame : bool
    Logical, If True then take the input coordinates and output fields as being with respect to the frame of
    referene of the element (ignore misalignments).

calc_dfield : bool, optional
    If present and True then calculate the field derivatives.

calc_potential : bool, optional
    Calc electric and magnetic potentials? Default is false. This is experimental and only implemented for
    wigglers at present.

use_overlap : bool, optional
    Add in overlap fields from other elements? Default is True.

grid_allow_s_out_of_bounds : bool, optional
    For grids, allow s-coordinate to be grossly out of bounds and return zero instead of an error? Default:
    False. Used internally for overlapping fields.

rf_time : float, optional
    Set the time relative to the RF clock. Normally this time is calculated using orbit.t or orbit.vec(5) but
    sometimes it is convenient to be able to override this. For example, time_runge_kutta uses this.

used_eles : 1D array of ElePointerStruct, optional
    For internal use only when this routine is called recursively. Used to prevent double counting when there
    is field overlap.

print_err : bool, optional
    Print an error message? Default is True. For example, if the particle is out of bounds when the field is
    defined on a grid.

original_ele : EleStruct, optional
    Used with recursive calls that pass the lord as the ele argument. In this case original_ele is the
    original ele argument.

Returns
-------
field : EmFieldStruct
    E and B fields and derivatives.

err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "em_field_derivatives",
      &Bmad::em_field_derivatives,
      py::arg("ele"),
      py::arg("param"),
      py::arg("s_pos"),
      py::arg("orbit"),
      py::arg("local_ref_frame"),
      py::arg("grid_allow_s_out_of_bounds") = py::none(),
      py::arg("rf_time") = py::none(),
      R"""(Subroutine em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, dfield, grid_allow_s_out_of_bounds, rf_time)

Routine to calculate field derivatives.
In theory this should be handled by em_filed_calc. In practice, em_field_calc is currently incomplete.

Parameters
----------
ele : EleStruct
    Element

param : LatParamStruct
    Lattice parameters.

s_pos : float
    Longitudinal position relative to the upstream edge of the element.

orbit : CoordStruct
    Transverse coordinates.

local_ref_frame : bool
    Logical, If True then take the input coordinates and output fields as being with respect to the frame of
    referene of the element (ignore misalignments).

grid_allow_s_out_of_bounds : bool, optional
    For grids, allow s-coordinate to be grossly out of bounds and return zero instead of an error? Default:
    False. Used internally for overlapping fields.

rf_time : float, optional
    RF clock time. If not present then the time will be calculated using the standard algorithm.

Returns
-------
dfield : EmFieldStruct
    E and B field derivatives. dfield.E and dfield.B are not touched.
)"""
  );
  m.def(
      "em_field_kick_vector_time",
      &Bmad::em_field_kick_vector_time,
      py::arg("ele"),
      py::arg("param"),
      py::arg("rf_time"),
      py::arg("orbit"),
      py::arg("err_flag"),
      py::arg("print_err") = py::none(),
      py::arg("extra_field") = py::none(),
      R"""(Subroutine em_field_kick_vector_time (ele, param, rf_time, orbit, dvec_dt, err_flag, print_err, extra_field))

Subroutine to convert particle coordinates from t-based to s-based system.

Parameters
----------
ele : EleStruct
    input particle

param : LatParamStruct
    Reference momentum. The sign indicates direction of p_s.

rf_time : float
    RF time.

orbit : CoordStruct
    in t-based system

err_flag : bool
    Set True if there is an error. False otherwise.

print_err : bool, optional
    Passed to em_field_calc

extra_field : EmFieldStruct, optional
    Static field to be added to the element field. Eg used with space charge.

Returns
-------
dvec_dt : 1D array of float (shape: 10)
    Derivatives.
)"""
  );
  m.def(
      "em_field_plus_em_field",
      &Bmad::em_field_plus_em_field,
      py::arg("field1"),
      py::arg("field2"),
      py::arg("field_tot"),
      R"""(Wrapper for Fortran routine em_field_plus_em_field

Parameters
----------
field1 : EmFieldStruct

field2 : EmFieldStruct

field_tot : EmFieldStruct
)"""
  );
  m.def(
      "em_taylor_equal_em_taylor",
      &Bmad::em_taylor_equal_em_taylor,
      py::arg("em_taylor1"),
      py::arg("em_taylor2"),
      R"""(Wrapper for Fortran routine em_taylor_equal_em_taylor

Parameters
----------
em_taylor1 : EmTaylorStruct

em_taylor2 : EmTaylorStruct
)"""
  );
  m.def(
      "em_taylors_equal_em_taylors",
      &Bmad::em_taylors_equal_em_taylors,
      py::arg("em_taylor1"),
      py::arg("em_taylor2"),
      R"""(Wrapper for Fortran routine em_taylors_equal_em_taylors

Parameters
----------
em_taylor1 : 1D array of EmTaylorStruct

em_taylor2 : 1D array of EmTaylorStruct
)"""
  );
  py::class_<Bmad::Emit6d, std::unique_ptr<Bmad::Emit6d>>(m, "Emit6d", "emit_6d return type")
      .def_readonly("mode", &Bmad::Emit6d::mode)
      .def_readonly("sigma_mat", &Bmad::Emit6d::sigma_mat)
      .def_readonly("rad_int_by_ele", &Bmad::Emit6d::rad_int_by_ele)
      .def("__len__", [](const Bmad::Emit6d &) { return 3; })
      .def("__getitem__", [](const Bmad::Emit6d &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.mode);
        if (i == 1)
          return py::cast(s.sigma_mat);
        if (i == 2)
          return py::cast(s.rad_int_by_ele);
        throw py::index_error();
      });
  m.def(
      "emit_6d",
      &Bmad::emit_6d,
      py::arg("ele_ref"),
      py::arg("include_opening_angle"),
      py::arg("closed_orbit") = py::none(),
      R"""(Subroutine emit_6d (ele_ref, include_opening_angle, mode, sigma_mat, closed_orbit, rad_int_by_ele)

Routine to calculate the three normal mode emittances, damping partition numbers, radiation integrals, etc.
Since the emattances, etc. are only an invariant in the limit of zero damping, the calculated
values will vary depending upon the reference element.

If the lattice geometry is open, only the radiation integrals is computed.

Parameters
----------
ele_ref : EleStruct
    Origin of the 1-turn maps used to evaluate the emittances.

include_opening_angle : bool
    If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
    comparing against other codes.

closed_orbit : 1D array of CoordStruct, optional
    Closed orbit. If not present this routine will calculate it.

Returns
-------
mode : NormalModesStruct
    Emittance and other info.

sigma_mat : 2D array of float (shape: 6,6)
    Sigma matrix.

rad_int_by_ele : RadIntAllEleStruct, optional
    Radiation integrals element-by-element.
)"""
  );
  m.def(
      "entering_element",
      &Bmad::entering_element,
      py::arg("orbit"),
      py::arg("particle_at"),
      R"""(Wrapper for Fortran routine entering_element

Parameters
----------
orbit : CoordStruct
    Particle orbit.

particle_at : int
    First_track_edge$ or second_track_edge$

Returns
-------
is_entering : bool
    Set True if particle is going from outside to inside and vice versa.
)"""
  );
  m.def(
      "envelope_radints",
      &Bmad::envelope_radints,
      py::arg("Lambda"),
      py::arg("Theta"),
      py::arg("Iota"),
      py::arg("alpha"),
      py::arg("emit"),
      R"""(subroutine envelope_radints(Lambda,Theta,Iota,alpha,emit)

Calculates damping decrement and emittance of the three
normal modes from the integrate diffusion, damping, and vertical
excitation matrices names Lambda, Theta, and Iota, respectively.
These three matrices are obtained from the subroutine integrated_mats.

The damping times can obtained from alpha using:
   tau = lattice_length/c_light/alpha
)"""
  );
  py::class_<Bmad::EnvelopeRadintsIbs, std::unique_ptr<Bmad::EnvelopeRadintsIbs>>(
      m,
      "EnvelopeRadintsIbs",
      "envelope_radints_ibs return type"
  )
      .def_readonly("alpha", &Bmad::EnvelopeRadintsIbs::alpha)
      .def_readonly("emit", &Bmad::EnvelopeRadintsIbs::emit)
      .def("__len__", [](const Bmad::EnvelopeRadintsIbs &) { return 2; })
      .def("__getitem__", [](const Bmad::EnvelopeRadintsIbs &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.alpha);
        if (i == 1)
          return py::cast(s.emit);
        throw py::index_error();
      });
  m.def(
      "envelope_radints_ibs",
      &Bmad::envelope_radints_ibs,
      py::arg("Lambda"),
      py::arg("Theta"),
      py::arg("Iota"),
      py::arg("eles"),
      py::arg("mode"),
      py::arg("tail_cut"),
      py::arg("npart"),
      py::arg("species"),
      R"""(subroutine envelope_radints_ibs(Lambda, Theta, Iota, eles, alpha, emit, mode, tail_cut, npart, species)

Calculates damping decrement and emittance of the three
normal modes by integrating the IBS, SR diffusion, and SR damping matrices.

The IBS depends on the envelope, and so this routine iterates to
locate the equilibrium beam envelope. This iterative process can fail to converge.

The damping times can obtained from alpha using:
   tau = lattice_length/c_light/alpha

alpha and emit are quantities for the three normal modes.
alpha and emit are ordered by plane dominance.

Only radiation from sbends and rbends is taken into account.
The one-turn transfer matrix at each element (slice) is obtained
by concatenating the individual element transfer matrices.

Parameters
----------
Lambda : 2D array of complex (shape: 6,6)
    Integrated damping matrix.

Theta : 2D array of complex (shape: 6,6)
    Integrated diffusion matrix.

Iota : 2D array of complex (shape: 6,6)
    Integrated vertical excitation matrix.

eles : 1D array of EleStruct
    array of element structures representing ring.

mode : NormalModesStruct
    normal_modes_struct

tail_cut : bool
    apply tail cut.

npart : float
    number of particles in bunch.

species : int
    Particle species.

Returns
-------
alpha : 1D array of float (shape: 3)
    Normal mode damping decrements.

emit : 1D array of float (shape: 3)
    Normal mode emittances.
)"""
  );
  m.def(
      "eq_ac_kicker",
      &Bmad::eq_ac_kicker,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_ac_kicker

Parameters
----------
f1 : AcKickerStruct

f2 : AcKickerStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_ac_kicker_freq",
      &Bmad::eq_ac_kicker_freq,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_ac_kicker_freq

Parameters
----------
f1 : AcKickerFreqStruct

f2 : AcKickerFreqStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_ac_kicker_time",
      &Bmad::eq_ac_kicker_time,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_ac_kicker_time

Parameters
----------
f1 : AcKickerTimeStruct

f2 : AcKickerTimeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_anormal_mode",
      &Bmad::eq_anormal_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_anormal_mode

Parameters
----------
f1 : AnormalModeStruct

f2 : AnormalModeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_aperture_param",
      &Bmad::eq_aperture_param,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_aperture_param

Parameters
----------
f1 : ApertureParamStruct

f2 : ApertureParamStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_aperture_point",
      &Bmad::eq_aperture_point,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_aperture_point

Parameters
----------
f1 : AperturePointStruct

f2 : AperturePointStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_aperture_scan",
      &Bmad::eq_aperture_scan,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_aperture_scan

Parameters
----------
f1 : ApertureScanStruct

f2 : ApertureScanStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_beam",
      &Bmad::eq_beam,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_beam

Parameters
----------
f1 : BeamStruct

f2 : BeamStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_beam_init",
      &Bmad::eq_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_beam_init

Parameters
----------
f1 : BeamInitStruct

f2 : BeamInitStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_bmad_common",
      &Bmad::eq_bmad_common,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_bmad_common

Parameters
----------
f1 : BmadCommonStruct

f2 : BmadCommonStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_bookkeeping_state",
      &Bmad::eq_bookkeeping_state,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_bookkeeping_state

Parameters
----------
f1 : BookkeepingStateStruct

f2 : BookkeepingStateStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_bpm_phase_coupling",
      &Bmad::eq_bpm_phase_coupling,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_bpm_phase_coupling

Parameters
----------
f1 : BpmPhaseCouplingStruct

f2 : BpmPhaseCouplingStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_branch",
      &Bmad::eq_branch,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_branch

Parameters
----------
f1 : BranchStruct

f2 : BranchStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_bunch",
      &Bmad::eq_bunch,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_bunch

Parameters
----------
f1 : BunchStruct

f2 : BunchStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_bunch_params",
      &Bmad::eq_bunch_params,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_bunch_params

Parameters
----------
f1 : BunchParamsStruct

f2 : BunchParamsStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_cartesian_map",
      &Bmad::eq_cartesian_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_cartesian_map

Parameters
----------
f1 : CartesianMapStruct

f2 : CartesianMapStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_cartesian_map_term",
      &Bmad::eq_cartesian_map_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_cartesian_map_term

Parameters
----------
f1 : CartesianMapTermStruct

f2 : CartesianMapTermStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_cartesian_map_term1",
      &Bmad::eq_cartesian_map_term1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_cartesian_map_term1

Parameters
----------
f1 : CartesianMapTerm1Struct

f2 : CartesianMapTerm1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_complex_taylor",
      &Bmad::eq_complex_taylor,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_complex_taylor

Parameters
----------
f1 : ComplexTaylorStruct

f2 : ComplexTaylorStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_complex_taylor_term",
      &Bmad::eq_complex_taylor_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_complex_taylor_term

Parameters
----------
f1 : ComplexTaylorTermStruct

f2 : ComplexTaylorTermStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_control",
      &Bmad::eq_control,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_control

Parameters
----------
f1 : ControlStruct

f2 : ControlStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_control_ramp1",
      &Bmad::eq_control_ramp1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_control_ramp1

Parameters
----------
f1 : ControlRamp1Struct

f2 : ControlRamp1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_control_var1",
      &Bmad::eq_control_var1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_control_var1

Parameters
----------
f1 : ControlVar1Struct

f2 : ControlVar1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_controller",
      &Bmad::eq_controller,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_controller

Parameters
----------
f1 : ControllerStruct

f2 : ControllerStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_coord",
      &Bmad::eq_coord,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_coord

Parameters
----------
f1 : CoordStruct

f2 : CoordStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_coord_array",
      &Bmad::eq_coord_array,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_coord_array

Parameters
----------
f1 : CoordArrayStruct

f2 : CoordArrayStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_cylindrical_map",
      &Bmad::eq_cylindrical_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_cylindrical_map

Parameters
----------
f1 : CylindricalMapStruct

f2 : CylindricalMapStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_cylindrical_map_term",
      &Bmad::eq_cylindrical_map_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_cylindrical_map_term

Parameters
----------
f1 : CylindricalMapTermStruct

f2 : CylindricalMapTermStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_cylindrical_map_term1",
      &Bmad::eq_cylindrical_map_term1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_cylindrical_map_term1

Parameters
----------
f1 : CylindricalMapTerm1Struct

f2 : CylindricalMapTerm1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_ele",
      &Bmad::eq_ele,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_ele

Parameters
----------
f1 : EleStruct

f2 : EleStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_ellipse_beam_init",
      &Bmad::eq_ellipse_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_ellipse_beam_init

Parameters
----------
f1 : EllipseBeamInitStruct

f2 : EllipseBeamInitStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_em_field",
      &Bmad::eq_em_field,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_em_field

Parameters
----------
f1 : EmFieldStruct

f2 : EmFieldStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_em_taylor",
      &Bmad::eq_em_taylor,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_em_taylor

Parameters
----------
f1 : EmTaylorStruct

f2 : EmTaylorStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_em_taylor_term",
      &Bmad::eq_em_taylor_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_em_taylor_term

Parameters
----------
f1 : EmTaylorTermStruct

f2 : EmTaylorTermStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_expression_atom",
      &Bmad::eq_expression_atom,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_expression_atom

Parameters
----------
f1 : ExpressionAtomStruct

f2 : ExpressionAtomStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_floor_position",
      &Bmad::eq_floor_position,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_floor_position

Parameters
----------
f1 : FloorPositionStruct

f2 : FloorPositionStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_gen_grad1",
      &Bmad::eq_gen_grad1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_gen_grad1

Parameters
----------
f1 : GenGrad1Struct

f2 : GenGrad1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_gen_grad_map",
      &Bmad::eq_gen_grad_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_gen_grad_map

Parameters
----------
f1 : GenGradMapStruct

f2 : GenGradMapStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_grid_beam_init",
      &Bmad::eq_grid_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_grid_beam_init

Parameters
----------
f1 : GridBeamInitStruct

f2 : GridBeamInitStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_grid_field",
      &Bmad::eq_grid_field,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_grid_field

Parameters
----------
f1 : GridFieldStruct

f2 : GridFieldStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_grid_field_pt",
      &Bmad::eq_grid_field_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_grid_field_pt

Parameters
----------
f1 : GridFieldPtStruct

f2 : GridFieldPtStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_grid_field_pt1",
      &Bmad::eq_grid_field_pt1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_grid_field_pt1

Parameters
----------
f1 : GridFieldPt1Struct

f2 : GridFieldPt1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_high_energy_space_charge",
      &Bmad::eq_high_energy_space_charge,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_high_energy_space_charge

Parameters
----------
f1 : HighEnergySpaceChargeStruct

f2 : HighEnergySpaceChargeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_interval1_coef",
      &Bmad::eq_interval1_coef,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_interval1_coef

Parameters
----------
f1 : Interval1CoefStruct

f2 : Interval1CoefStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_kv_beam_init",
      &Bmad::eq_kv_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_kv_beam_init

Parameters
----------
f1 : KvBeamInitStruct

f2 : KvBeamInitStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_lat",
      &Bmad::eq_lat,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_lat

Parameters
----------
f1 : LatStruct

f2 : LatStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_lat_ele_loc",
      &Bmad::eq_lat_ele_loc,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_lat_ele_loc

Parameters
----------
f1 : LatEleLocStruct

f2 : LatEleLocStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_lat_param",
      &Bmad::eq_lat_param,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_lat_param

Parameters
----------
f1 : LatParamStruct

f2 : LatParamStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_linac_normal_mode",
      &Bmad::eq_linac_normal_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_linac_normal_mode

Parameters
----------
f1 : LinacNormalModeStruct

f2 : LinacNormalModeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_mode3",
      &Bmad::eq_mode3,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_mode3

Parameters
----------
f1 : Mode3Struct

f2 : Mode3Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_mode_info",
      &Bmad::eq_mode_info,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_mode_info

Parameters
----------
f1 : ModeInfoStruct

f2 : ModeInfoStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_normal_modes",
      &Bmad::eq_normal_modes,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_normal_modes

Parameters
----------
f1 : NormalModesStruct

f2 : NormalModesStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_photon_element",
      &Bmad::eq_photon_element,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_photon_element

Parameters
----------
f1 : PhotonElementStruct

f2 : PhotonElementStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_photon_material",
      &Bmad::eq_photon_material,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_photon_material

Parameters
----------
f1 : PhotonMaterialStruct

f2 : PhotonMaterialStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_photon_reflect_surface",
      &Bmad::eq_photon_reflect_surface,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_photon_reflect_surface

Parameters
----------
f1 : PhotonReflectSurfaceStruct

f2 : PhotonReflectSurfaceStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_photon_reflect_table",
      &Bmad::eq_photon_reflect_table,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_photon_reflect_table

Parameters
----------
f1 : PhotonReflectTableStruct

f2 : PhotonReflectTableStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_photon_target",
      &Bmad::eq_photon_target,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_photon_target

Parameters
----------
f1 : PhotonTargetStruct

f2 : PhotonTargetStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_pixel_detec",
      &Bmad::eq_pixel_detec,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_pixel_detec

Parameters
----------
f1 : PixelDetecStruct

f2 : PixelDetecStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_pixel_pt",
      &Bmad::eq_pixel_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_pixel_pt

Parameters
----------
f1 : PixelPtStruct

f2 : PixelPtStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_pre_tracker",
      &Bmad::eq_pre_tracker,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_pre_tracker

Parameters
----------
f1 : PreTrackerStruct

f2 : PreTrackerStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_rad_int1",
      &Bmad::eq_rad_int1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_rad_int1

Parameters
----------
f1 : RadInt1Struct

f2 : RadInt1Struct

is_eq : bool
)"""
  );
  m.def(
      "eq_rad_int_all_ele",
      &Bmad::eq_rad_int_all_ele,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_rad_int_all_ele

Parameters
----------
f1 : RadIntAllEleStruct

f2 : RadIntAllEleStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_rad_int_branch",
      &Bmad::eq_rad_int_branch,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_rad_int_branch

Parameters
----------
f1 : RadIntBranchStruct

f2 : RadIntBranchStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_rad_map",
      &Bmad::eq_rad_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_rad_map

Parameters
----------
f1 : RadMapStruct

f2 : RadMapStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_rad_map_ele",
      &Bmad::eq_rad_map_ele,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_rad_map_ele

Parameters
----------
f1 : RadMapEleStruct

f2 : RadMapEleStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_ramper_lord",
      &Bmad::eq_ramper_lord,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_ramper_lord

Parameters
----------
f1 : RamperLordStruct

f2 : RamperLordStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_space_charge_common",
      &Bmad::eq_space_charge_common,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_space_charge_common

Parameters
----------
f1 : SpaceChargeCommonStruct

f2 : SpaceChargeCommonStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_spin_polar",
      &Bmad::eq_spin_polar,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_spin_polar

Parameters
----------
f1 : SpinPolarStruct

f2 : SpinPolarStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_spline",
      &Bmad::eq_spline,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_spline

Parameters
----------
f1 : SplineStruct

f2 : SplineStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_strong_beam",
      &Bmad::eq_strong_beam,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_strong_beam

Parameters
----------
f1 : StrongBeamStruct

f2 : StrongBeamStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_curvature",
      &Bmad::eq_surface_curvature,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_curvature

Parameters
----------
f1 : SurfaceCurvatureStruct

f2 : SurfaceCurvatureStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_displacement",
      &Bmad::eq_surface_displacement,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_displacement

Parameters
----------
f1 : SurfaceDisplacementStruct

f2 : SurfaceDisplacementStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_displacement_pt",
      &Bmad::eq_surface_displacement_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_displacement_pt

Parameters
----------
f1 : SurfaceDisplacementPtStruct

f2 : SurfaceDisplacementPtStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_h_misalign",
      &Bmad::eq_surface_h_misalign,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_h_misalign

Parameters
----------
f1 : SurfaceHMisalignStruct

f2 : SurfaceHMisalignStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_h_misalign_pt",
      &Bmad::eq_surface_h_misalign_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_h_misalign_pt

Parameters
----------
f1 : SurfaceHMisalignPtStruct

f2 : SurfaceHMisalignPtStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_segmented",
      &Bmad::eq_surface_segmented,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_segmented

Parameters
----------
f1 : SurfaceSegmentedStruct

f2 : SurfaceSegmentedStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_surface_segmented_pt",
      &Bmad::eq_surface_segmented_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_surface_segmented_pt

Parameters
----------
f1 : SurfaceSegmentedPtStruct

f2 : SurfaceSegmentedPtStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_target_point",
      &Bmad::eq_target_point,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_target_point

Parameters
----------
f1 : TargetPointStruct

f2 : TargetPointStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_taylor",
      &Bmad::eq_taylor,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_taylor

Parameters
----------
f1 : TaylorStruct

f2 : TaylorStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_taylor_term",
      &Bmad::eq_taylor_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_taylor_term

Parameters
----------
f1 : TaylorTermStruct

f2 : TaylorTermStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_track",
      &Bmad::eq_track,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_track

Parameters
----------
f1 : TrackStruct

f2 : TrackStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_track_point",
      &Bmad::eq_track_point,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_track_point

Parameters
----------
f1 : TrackPointStruct

f2 : TrackPointStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_twiss",
      &Bmad::eq_twiss,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_twiss

Parameters
----------
f1 : TwissStruct

f2 : TwissStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wake",
      &Bmad::eq_wake,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wake

Parameters
----------
f1 : WakeStruct

f2 : WakeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wake_lr",
      &Bmad::eq_wake_lr,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wake_lr

Parameters
----------
f1 : WakeLrStruct

f2 : WakeLrStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wake_lr_mode",
      &Bmad::eq_wake_lr_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wake_lr_mode

Parameters
----------
f1 : WakeLrModeStruct

f2 : WakeLrModeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wake_sr",
      &Bmad::eq_wake_sr,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wake_sr

Parameters
----------
f1 : WakeSrStruct

f2 : WakeSrStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wake_sr_mode",
      &Bmad::eq_wake_sr_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wake_sr_mode

Parameters
----------
f1 : WakeSrModeStruct

f2 : WakeSrModeStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wake_sr_z_long",
      &Bmad::eq_wake_sr_z_long,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wake_sr_z_long

Parameters
----------
f1 : WakeSrZLongStruct

f2 : WakeSrZLongStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wall3d",
      &Bmad::eq_wall3d,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wall3d

Parameters
----------
f1 : Wall3dStruct

f2 : Wall3dStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wall3d_section",
      &Bmad::eq_wall3d_section,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wall3d_section

Parameters
----------
f1 : Wall3dSectionStruct

f2 : Wall3dSectionStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_wall3d_vertex",
      &Bmad::eq_wall3d_vertex,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_wall3d_vertex

Parameters
----------
f1 : Wall3dVertexStruct

f2 : Wall3dVertexStruct

is_eq : bool
)"""
  );
  m.def(
      "eq_xy_disp",
      &Bmad::eq_xy_disp,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Wrapper for Fortran routine eq_xy_disp

Parameters
----------
f1 : XyDispStruct

f2 : XyDispStruct

is_eq : bool
)"""
  );
  m.def(
      "equal_sign_here",
      &Bmad::equal_sign_here,
      py::arg("ele"),
      py::arg("delim"),
      py::arg("is_here"),
      R"""(Wrapper for Fortran routine equal_sign_here

Parameters
----------
ele : EleStruct

delim : character

is_here : bool
)"""
  );
  m.def(
      "equivalent_taylor_attributes",
      &Bmad::equivalent_taylor_attributes,
      py::arg("ele_taylor"),
      py::arg("ele2"),
      R"""(Wrapper for Fortran routine equivalent_taylor_attributes

Parameters
----------
ele_taylor : EleStruct
    Element with a Taylor map

ele2 : EleStruct
    Element that might receive the Taylor map from ele_taylor.

Returns
-------
equiv : bool
    True if elements are equivalent.
)"""
  );
  m.def(
      "etdiv",
      &Bmad::etdiv,
      py::arg("A"),
      py::arg("B"),
      py::arg("C"),
      py::arg("D"),
      py::arg("E"),
      py::arg("F"),
      R"""(Wrapper for Fortran routine etdiv

Parameters
----------
A : float

B : float

C : float

D : float

E : float

F : float
)"""
  );
  py::class_<Bmad::EvaluateArrayIndex, std::unique_ptr<Bmad::EvaluateArrayIndex>>(
      m,
      "EvaluateArrayIndex",
      "evaluate_array_index return type"
  )
      .def_readonly("err_flag", &Bmad::EvaluateArrayIndex::err_flag)
      .def_readonly("word2", &Bmad::EvaluateArrayIndex::word2)
      .def_readonly("delim2", &Bmad::EvaluateArrayIndex::delim2)
      .def_readonly("this_index", &Bmad::EvaluateArrayIndex::this_index)
      .def("__len__", [](const Bmad::EvaluateArrayIndex &) { return 4; })
      .def("__getitem__", [](const Bmad::EvaluateArrayIndex &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.word2);
        if (i == 2)
          return py::cast(s.delim2);
        if (i == 3)
          return py::cast(s.this_index);
        throw py::index_error();
      });
  m.def(
      "evaluate_array_index",
      &Bmad::evaluate_array_index,
      py::arg("delim_list1"),
      py::arg("delim_list2"),
      R"""(Function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)

Function of evaluate the index of an array. Typically the text being parsed looks like:
     "5) = ..."         or
     "6).COMP = ..."

Parameters
----------
delim_list1 : character
    Delimitor after the integer. Normally ')'.

delim_list2 : character
    Delimitor list to mark the end of word2. Normally '='.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.

word2 : character
    Word found after delim1. Normally this should be blank.

delim2 : character
    Actual delimitor found after word2.

this_index : int
    Integer value
)"""
  );
  py::class_<Bmad::EvaluateLogical, std::unique_ptr<Bmad::EvaluateLogical>>(
      m,
      "EvaluateLogical",
      "evaluate_logical return type"
  )
      .def_readonly("iostat", &Bmad::EvaluateLogical::iostat)
      .def_readonly("this_logic", &Bmad::EvaluateLogical::this_logic)
      .def("__len__", [](const Bmad::EvaluateLogical &) { return 2; })
      .def("__getitem__", [](const Bmad::EvaluateLogical &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.iostat);
        if (i == 1)
          return py::cast(s.this_logic);
        throw py::index_error();
      });
  m.def(
      "evaluate_logical",
      &Bmad::evaluate_logical,
      py::arg("word"),
      R"""(Function evaluate_logical (word, iostat) result (this_logic)

Function of convert a string into a logical value.
Accepted possibilities are:
  .TRUE.  .FALSE.
   TRUE    FALSE
   T       F

Parameters
----------
word : character
    Input string.

Returns
-------
iostat : int
    Status: Returns 0 if conversion successful.

this_logic : bool
    Result.
)"""
  );
  m.def(
      "exact_bend_edge_kick",
      &Bmad::exact_bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orb"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine exact_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

Subroutine to track through the edge field of an sbend.
Uses routines adapted from PTC

Parameters
----------
ele : EleStruct
    SBend element.

param : LatParamStruct

particle_at : int
    first_track_edge$, or second_track_edge$.

orb : CoordStruct
    Coords after tracking.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the edge.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix through the edge.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.
)"""
  );
  m.def(
      "exp_bessi0",
      &Bmad::exp_bessi0,
      py::arg("t"),
      py::arg("B1"),
      py::arg("B2"),
      py::arg("func_retval__"),
      R"""(Function exp_bessi0(t, B1, B2)

This is essentially the Numercal Recipes bessi0 function multiplied by exp(-B1*t).

This overcomes an issue where exp(B2*t) may be huge and exp(-B1*t) may be small.
Evaluating exp(B2*t) may result in overflow, but exp((B2-B1)*t) has a moderate value.
Simplifying the algebra of B2-B1 suggests that is should always have a moderate magnitude.

Parameters
----------
t : float
    Scalar agrument to evaluate function at.

B1 : float
    Scalar value.  Eq. 33 from Piwinski's paper.

B2 : float
    Scalar value.  Eq. 34 from Piwinski's paper.
)"""
  );
  py::class_<PyExpectOneOf, std::unique_ptr<PyExpectOneOf>>(
      m,
      "ExpectOneOf",
      "expect_one_of return type"
  )
      .def_readonly("delim", &PyExpectOneOf::delim)
      .def("__len__", [](const PyExpectOneOf &) { return 1; })
      .def("__getitem__", [](const PyExpectOneOf &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.delim);
        throw py::index_error();
      });
  m.def(
      "expect_one_of",
      &python_expect_one_of,
      py::arg("delim_list"),
      py::arg("check_input_delim"),
      py::arg("ele_name"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("is_ok"),
      R"""(Function expect_one_of (delim_list, check_input_delim, ele_name, delim, delim_found) result (is_ok)

Routine to check either that the current delimitor or the next character in the parse stream is the
expected delimitor.
This routine is used for Bmad lattice file parsing and is not meant for general use.

Also see: expect_this

Parameters
----------
delim_list : character
    List of expected (valid) delimitors. If list contains a space character then no delimitor (indicating the
    end of the command) is a valid possibility.

check_input_delim : bool
    If True, then check if delim argument is in the delim_list. If False, check that the next character in the
    parse stream is an expected delimitor.

ele_name : character
    Lattice element under construction. Used for error messages.

delim : character
    Current delimitor that will be checked if check_input_delim = .true.
    This parameter is an input/output and is modified in-place.
    As an output, delim: Next delim if check_input_delim = False.

Returns
-------
delim : character
    Current delimitor that will be checked if check_input_delim = .true.
    This parameter is an input/output and is modified in-place.
    As an output, delim: Next delim if check_input_delim = False.
)"""
  );
  py::class_<Bmad::ExpectThis, std::unique_ptr<Bmad::ExpectThis>>(
      m,
      "ExpectThis",
      "expect_this return type"
  )
      .def_readonly("delim", &Bmad::ExpectThis::delim)
      .def_readonly("delim_found", &Bmad::ExpectThis::delim_found)
      .def("__len__", [](const Bmad::ExpectThis &) { return 2; })
      .def("__getitem__", [](const Bmad::ExpectThis &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        throw py::index_error();
      });
  m.def(
      "expect_this",
      &Bmad::expect_this,
      py::arg("expecting"),
      py::arg("check_delim"),
      py::arg("call_check"),
      py::arg("err_str"),
      py::arg("ele"),
      py::arg("is_ok"),
      R"""(Function expect_this (expecting, check_delim, call_check, err_str, ele, delim, delim_found) result (is_ok)

Checks that the next character or characters in the parse stream corresponds to the
characters in the expecting argument. For example, if expecting is ')={' these three characters
should be the next non-blank characters in the parse stream.

Also see: expect_one_of

Parameters
----------
expecting : character
    list of characters that are expected to be next in the parse stream.

check_delim : bool
    If True then use delim argument as first token to check. A blank character indicates end of command is
    expected.

call_check : bool
    If True then check for 'call::<filename>' construct.

err_str : character
    String used for error messages.

ele : EleStruct
    Element parameters being parsed.

Returns
-------
delim : character
    Final delim

delim_found : bool
    Is there a final delim (as opposed to end of command).
)"""
  );
  m.def(
      "expression_stack_to_string",
      &Bmad::expression_stack_to_string,
      py::arg("stack"),
      py::arg("polish") = py::none(),
      R"""(Function expression_stack_to_string (stack, polish) result (str)

Routine to convert an expression stack to a string

Parameters
----------
stack : 1D array of ExpressionAtomStruct
    arithmetic expression

polish : bool, optional
    logical, optional, Construct expression in reverse polish? Default is False.

Returns
-------
str : character
    : Expression in string form.
)"""
  );
  py::class_<Bmad::ExpressionStackValue, std::unique_ptr<Bmad::ExpressionStackValue>>(
      m,
      "ExpressionStackValue",
      "expression_stack_value return type"
  )
      .def_readonly("err_flag", &Bmad::ExpressionStackValue::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionStackValue::err_str)
      .def_readonly("value", &Bmad::ExpressionStackValue::value)
      .def("__len__", [](const Bmad::ExpressionStackValue &) { return 3; })
      .def("__getitem__", [](const Bmad::ExpressionStackValue &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.err_str);
        if (i == 2)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "expression_stack_value",
      &Bmad::expression_stack_value,
      py::arg("stack"),
      py::arg("var") = py::none(),
      py::arg("use_old") = py::none(),
      R"""(Function expression_stack_value (stack, err_flag, err_str, var, use_old) result (value)

Routine to evaluate a mathematical expression represented by an "expression stack".
Expression stacks are created by expression_string_to_stack.

Note: Stack elements with stack(i)%type == variable$ need to be evalauated before
calling this routine and the value placed in stack(i)%value.

Also see:
  expression_value
  expression_string_to_stack

Parameters
----------
stack : 1D array of ExpressionAtomStruct
    Expression to evaluate.

var : 1D array of ControlVar1Struct, optional
    Array of control variables. Used with Bmad controller elements.

use_old : bool, optional
    Use var.old_value? Must be present if var(:) is present.

Returns
-------
err_flag : bool
    True if there is an evaluation problem. False otherwise.

err_str : character
    Error string explaining error if there is one.

value : float
    Value of the expression.
)"""
  );
  py::class_<Bmad::ExpressionStringToStack, std::unique_ptr<Bmad::ExpressionStringToStack>>(
      m,
      "ExpressionStringToStack",
      "expression_string_to_stack return type"
  )
      .def_readonly("stack", &Bmad::ExpressionStringToStack::stack)
      .def_readonly("n_stack", &Bmad::ExpressionStringToStack::n_stack)
      .def_readonly("err_flag", &Bmad::ExpressionStringToStack::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionStringToStack::err_str)
      .def("__len__", [](const Bmad::ExpressionStringToStack &) { return 4; })
      .def("__getitem__", [](const Bmad::ExpressionStringToStack &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.stack);
        if (i == 1)
          return py::cast(s.n_stack);
        if (i == 2)
          return py::cast(s.err_flag);
        if (i == 3)
          return py::cast(s.err_str);
        throw py::index_error();
      });
  m.def(
      "expression_string_to_stack",
      &Bmad::expression_string_to_stack,
      py::arg("string"),
      R"""(Subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)

This routine creates an expression stack array which can be used
to evaluate an arithmethic expression.

Stack end elements not used are marked stack(i)%type = end_stack$

Stack elements with stack(i)%type = variable$ are elements that need
to be evaluated before calling expression_stack_value.

Also see:
  expression_value
  expression_stack_value

Parameters
----------
string : character
    Expression to be converted.

Returns
-------
stack : 1D array of ExpressionAtomStruct
    Expression evaluation stack.

n_stack : int
    number of "atoms" used by the expression

err_flag : bool
    Set True if there is an error (EG divide by 0).

err_str : character
    String describing the error.
)"""
  );
  py::class_<Bmad::ExpressionStringToTree, std::unique_ptr<Bmad::ExpressionStringToTree>>(
      m,
      "ExpressionStringToTree",
      "expression_string_to_tree return type"
  )
      .def_readonly("err_flag", &Bmad::ExpressionStringToTree::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionStringToTree::err_str)
      .def("__len__", [](const Bmad::ExpressionStringToTree &) { return 2; })
      .def("__getitem__", [](const Bmad::ExpressionStringToTree &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.err_str);
        throw py::index_error();
      });
  m.def(
      "expression_string_to_tree",
      &Bmad::expression_string_to_tree,
      py::arg("string"),
      py::arg("root_tree"),
      R"""(Subroutine expression_string_to_tree (string, root_tree, err_flag, err_str)

Routine to create an expression tree array which can be used
to evaluate an arithmethic expression.

Also see:
  expression_value
  expression_tree_value
  deallocate_expression_tree

Important! trees use pointers as opposed to allocatable arrays due to the ifort compiler not being able to
handle node%node(:) being an allocatable array. Thus deallocate_expression_tree must be called before
any tree instance goes out of scope.

Note types used:
  plus$, minus$, times$, divide$, power$, unary_minus$, unary_plus$
  constant$, numeric$, variable$, function$
  root$, parens$, func_parens$, square_brackets$, curly_brackets$
  arrow$, equal$, colon$, double_colon$, vertical_bar$, compound$

An expression string will be split on:
  Two character operators: "->", "::"
  operators: + -  * / ^ = : &
  brackets: [] () {}
  comma: ,

Root node name is "root" and is of type root$
Brackets in the expression string must be matched.
The corresponding tree node will have a name / type of:
  "[]" / square_brackets$,    "()" / parens$ or func_parens$,   "{}" / curley_brackets$

The root node, equal nodes, and all bracket nodes, will have an array of child nodes all of which will be comma nodes.

Parameters
----------
string : character
    Expression to be converted.

root_tree : ExpressionTreeStruct
    Only used when recursively called.

Returns
-------
err_flag : bool
    Set True if there is an error (EG divide by 0).

err_str : character
    String describing the error. Make length large to hold the expression.
)"""
  );
  m.def(
      "expression_tree_to_string",
      &Bmad::expression_tree_to_string,
      py::arg("tree"),
      py::arg("include_root") = py::none(),
      py::arg("n_node") = py::none(),
      py::arg("parent") = py::none(),
      R"""(Function expression_tree_to_string (tree, include_root, n_node, parent) result(str_out)

Routine to convert an expression tree to a expression string.

Parameters
----------
tree : ExpressionTreeStruct
    Root of tree to print.

include_root : bool, optional
    Default is True. If True, do not inculde in the output string the root node. Note: If the root node is of
    type root$, this node is always ignored.

n_node : int, optional
    Node index. parent.node(n_node) === tree. Internal use only. Used with recursive calls.

parent : ExpressionTreeStruct, optional
    Internal use only. Used with recusive calls.

Returns
-------
str_out : character
    Expression string.
)"""
  );
  py::class_<Bmad::ExpressionValue, std::unique_ptr<Bmad::ExpressionValue>>(
      m,
      "ExpressionValue",
      "expression_value return type"
  )
      .def_readonly("err_flag", &Bmad::ExpressionValue::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionValue::err_str)
      .def_readonly("value", &Bmad::ExpressionValue::value)
      .def("__len__", [](const Bmad::ExpressionValue &) { return 3; })
      .def("__getitem__", [](const Bmad::ExpressionValue &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.err_str);
        if (i == 2)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "expression_value",
      &Bmad::expression_value,
      py::arg("expression"),
      py::arg("var") = py::none(),
      py::arg("use_old") = py::none(),
      R"""(Function expression_value (expression, err_flag, err_str, var, use_old) result (value)

Routine to evaluate a mathematical expression encoded in a string.

Also see:
  expression_string_to_stack
  expression_stack_value

Parameters
----------
expression : character
    Expression string.

var : 1D array of ControlVar1Struct, optional
    Array of control variables. Used with Bmad controller elements.

use_old : bool, optional
    Use var.old_value? Must be present if var(:) is present.

Returns
-------
err_flag : bool
    True if there is an evaluation problem. False otherwise.

err_str : character, optional
    Error string explaining error if there is one.

value : float
    Value of the expression.
)"""
  );
}
