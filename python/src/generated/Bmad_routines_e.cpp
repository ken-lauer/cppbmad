#include "pybmad/generated/Bmad_routines_e.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyExpectOneOf python_expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string delim,
    bool delim_found
) {
  auto _result = Bmad::expect_one_of(delim_list, check_input_delim, ele_name, delim, delim_found);
  auto py_result{PyExpectOneOf{_result, delim}};
  return py_result;
}

void init_Bmad_routines_e(nb::module_ &m) {
  m.def(
      "e_accel_field",
      &Bmad::e_accel_field,
      nb::arg("ele"),
      nb::arg("voltage_or_gradient"),
      nb::arg("bmad_standard_tracking") = nb::none(),
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
      nb::arg("gamma"),
      nb::arg("g_bend"),
      R"""(Routine to calculate the photon critical energy in a bend.

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
  nb::class_<Bmad::EigenDecomp6mat>(m, "EigenDecomp6mat", "eigen_decomp_6mat return type")
      .def_ro("eval", &Bmad::EigenDecomp6mat::eval)
      .def_ro("evec", &Bmad::EigenDecomp6mat::evec)
      .def_ro("err_flag", &Bmad::EigenDecomp6mat::err_flag)
      .def_ro("tunes", &Bmad::EigenDecomp6mat::tunes)
      .def("__len__", [](const Bmad::EigenDecomp6mat &) { return 4; })
      .def("__getitem__", [](const Bmad::EigenDecomp6mat &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.eval);
        if (i == 1)
          return nb::cast(s.evec);
        if (i == 2)
          return nb::cast(s.err_flag);
        if (i == 3)
          return nb::cast(s.tunes);
        throw nb::index_error();
      });
  m.def(
      "eigen_decomp_6mat",
      &Bmad::eigen_decomp_6mat,
      nb::arg("mat"),
      R"""(Compute eigenvalues and eigenvectors of a real 6x6 matrix.
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
      nb::arg("ele0"),
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("err_flag"),
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
      nb::arg("ele_out"),
      nb::arg("ele_in"),
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
      nb::arg("ele_in"),
      nb::arg("update_nametable"),
      R"""(Subroutine that is used to set an element equal to another.
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
      nb::arg("ele"),
      R"""(Finalizer routine for ele_struct instances.
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
      nb::arg("ele"),
      nb::arg("template_") = nb::none(),
      R"""(Wrapper for Fortran routine ele_full_name

Parameters
----------
ele : EleStruct
    Element in a lattice

Returns
-------
str : str
    : Name/location string.
)"""
  );
  m.def(
      "ele_geometry",
      &Bmad::ele_geometry,
      nb::arg("floor_start"),
      nb::arg("ele"),
      nb::arg("len_scale") = nb::none(),
      nb::arg("ignore_patch_err") = nb::none(),
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
      nb::arg("ele"),
      nb::arg("len_scale") = nb::none(),
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
      nb::arg("ele"),
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
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine ele_has_nonzero_kick

Parameters
----------
ele : EleStruct
    Element with possible nonzero kicks.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with no kicks.

Returns
-------
has_kick : bool
)"""
  );
  m.def(
      "ele_has_nonzero_offset",
      &Bmad::ele_has_nonzero_offset,
      nb::arg("ele"),
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
      nb::arg("ele"),
      nb::arg("print_warning") = nb::none(),
      R"""(Routine to check that an element is either a detector, instrument, monitor, or marker.
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
      nb::arg("ele"),
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
      nb::arg("ele"),
      nb::arg("show_branch0") = nb::none(),
      nb::arg("parens") = nb::none(),
      R"""(Wrapper for Fortran routine ele_loc_name

Parameters
----------
ele : EleStruct
    Element in a lattice

show_branch0 : bool, optional
    Explicitly show branch for main lattice elements? Default is False.

parens : str, optional
    If present, enclose location string using the two characters supplied. Typically parens will be set to
    "()" or "[]".

Returns
-------
str : str
    Output string. Left justified.
)"""
  );
  nb::class_<Bmad::EleMisalignmentLSCalc>(
      m,
      "EleMisalignmentLSCalc",
      "ele_misalignment_l_s_calc return type"
  )
      .def_ro("L_mis", &Bmad::EleMisalignmentLSCalc::L_mis)
      .def_ro("S_mis", &Bmad::EleMisalignmentLSCalc::S_mis)
      .def("__len__", [](const Bmad::EleMisalignmentLSCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::EleMisalignmentLSCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.L_mis);
        if (i == 1)
          return nb::cast(s.S_mis);
        throw nb::index_error();
      });
  m.def(
      "ele_misalignment_l_s_calc",
      &Bmad::ele_misalignment_l_s_calc,
      nb::arg("ele"),
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
      nb::arg("ele"),
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
      nb::arg("lat"),
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
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("particle_at"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
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
      nb::arg("E_ref"),
      nb::arg("s_rel"),
      nb::arg("ele"),
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
  nb::class_<Bmad::EleToFibre>(m, "EleToFibre", "ele_to_fibre return type")
      .def_ro("ptc_fibre", &Bmad::EleToFibre::ptc_fibre)
      .def_ro("err_flag", &Bmad::EleToFibre::err_flag)
      .def("__len__", [](const Bmad::EleToFibre &) { return 2; })
      .def("__getitem__", [](const Bmad::EleToFibre &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ptc_fibre);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "ele_to_fibre",
      [](EleStruct &ele,
         bool use_offsets,
         std::optional<int> integ_order,
         std::optional<int> steps,
         std::optional<bool> for_layout,
         CoordStruct *ref_in) {
        auto fn = static_cast<Bmad::EleToFibre (*)(
            EleStruct &,
            bool,
            std::optional<int>,
            std::optional<int>,
            std::optional<bool>,
            optional_ref<CoordStruct>
        )>(&Bmad::ele_to_fibre);
        return fn(ele, use_offsets, integ_order, steps, for_layout, ptr_to_opt_ref(ref_in));
      },
      nb::arg("ele"),
      nb::arg("use_offsets"),
      nb::arg("integ_order") = nb::none(),
      nb::arg("steps") = nb::none(),
      nb::arg("for_layout") = nb::none(),
      nb::arg("ref_in") = nb::none(),
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
err_flag : bool
    Set True if setup OK. False otherwise.

ptc_fibre : Fibre, optional
    PTC fibre element.
)"""
  );
  m.def(
      "ele_to_ptc_magnetic_bn_an",
      &Bmad::ele_to_ptc_magnetic_bn_an,
      nb::arg("ele"),
      nb::arg("bn"),
      nb::arg("an"),
      R"""(Routine to compute the a(n) and b(n) magnetic multipole components of a magnet.
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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("orb0"),
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
  nb::class_<Bmad::EleToTaylor>(m, "EleToTaylor", "ele_to_taylor return type")
      .def_ro("orbital_taylor", &Bmad::EleToTaylor::orbital_taylor)
      .def_ro("spin_taylor", &Bmad::EleToTaylor::spin_taylor)
      .def("__len__", [](const Bmad::EleToTaylor &) { return 2; })
      .def("__getitem__", [](const Bmad::EleToTaylor &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.orbital_taylor);
        if (i == 1)
          return nb::cast(s.spin_taylor);
        throw nb::index_error();
      });
  m.def(
      "ele_to_taylor",
      [](EleStruct &ele,
         CoordStruct *orb0,
         std::optional<bool> taylor_map_includes_offsets,
         std::optional<bool> include_damping) {
        auto fn = static_cast<Bmad::EleToTaylor (*)(
            EleStruct &,
            optional_ref<CoordStruct>,
            std::optional<bool>,
            std::optional<bool>
        )>(&Bmad::ele_to_taylor);
        return fn(ele, ptr_to_opt_ref(orb0), taylor_map_includes_offsets, include_damping);
      },
      nb::arg("ele"),
      nb::arg("orb0") = nb::none(),
      nb::arg("taylor_map_includes_offsets") = nb::none(),
      nb::arg("include_damping") = nb::none(),
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
      nb::arg("ele"),
      nb::arg("order"),
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
unique_name : str
    Unique name that can can be used to identify ele. The simplist name will be constructed. For example, if
    the element name is unique, unique_name will be set to the element name.
)"""
  );
  m.def(
      "ele_value_has_changed",
      &Bmad::ele_value_has_changed,
      nb::arg("ele"),
      nb::arg("list"),
      nb::arg("abs_tol"),
      nb::arg("set_old"),
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
      nb::arg("ele1"),
      nb::arg("ele2"),
      R"""(Wrapper for Fortran routine ele_vec_equal_ele_vec

Parameters
----------
ele1 : 1D array of EleStruct

ele2 : 1D array of EleStruct
)"""
  );
  nb::class_<Bmad::ElecMultipoleField>(m, "ElecMultipoleField", "elec_multipole_field return type")
      .def_ro("Ex", &Bmad::ElecMultipoleField::Ex)
      .def_ro("Ey", &Bmad::ElecMultipoleField::Ey)
      .def_ro("dE", &Bmad::ElecMultipoleField::dE)
      .def_ro("compute_dE", &Bmad::ElecMultipoleField::compute_dE)
      .def("__len__", [](const Bmad::ElecMultipoleField &) { return 4; })
      .def("__getitem__", [](const Bmad::ElecMultipoleField &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.Ex);
        if (i == 1)
          return nb::cast(s.Ey);
        if (i == 2)
          return nb::cast(s.dE);
        if (i == 3)
          return nb::cast(s.compute_dE);
        throw nb::index_error();
      });
  m.def(
      "elec_multipole_field",
      &Bmad::elec_multipole_field,
      nb::arg("a"),
      nb::arg("b"),
      nb::arg("n"),
      nb::arg("coord"),
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
  nb::class_<Bmad::ElementAtSBranch>(m, "ElementAtSBranch", "element_at_s_branch return type")
      .def_ro("err_flag", &Bmad::ElementAtSBranch::err_flag)
      .def_ro("s_eff", &Bmad::ElementAtSBranch::s_eff)
      .def_ro("position", &Bmad::ElementAtSBranch::position)
      .def_ro("ix_ele", &Bmad::ElementAtSBranch::ix_ele)
      .def("__len__", [](const Bmad::ElementAtSBranch &) { return 4; })
      .def("__getitem__", [](const Bmad::ElementAtSBranch &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.s_eff);
        if (i == 2)
          return nb::cast(s.position);
        if (i == 3)
          return nb::cast(s.ix_ele);
        throw nb::index_error();
      });
  m.def(
      "element_at_s",
      nb::overload_cast<BranchStruct &, double, bool, std::optional<bool>>(&Bmad::element_at_s),
      nb::arg("branch"),
      nb::arg("s"),
      nb::arg("choose_max"),
      nb::arg("print_err") = nb::none(),
      R"""(Function to return the index of the element at position s.

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
ix_ele : int
    Index of element at s.

err_flag : bool, optional
    Set True if s is out of bounds. False otherwise.

s_eff : float, optional
    Effective s. Equal to s with a open lattice. See above.

position : CoordStruct, optional
    Positional information.
)"""
  );
  nb::class_<Bmad::ElementAtSLat>(m, "ElementAtSLat", "element_at_s_lat return type")
      .def_ro("err_flag", &Bmad::ElementAtSLat::err_flag)
      .def_ro("s_eff", &Bmad::ElementAtSLat::s_eff)
      .def_ro("position", &Bmad::ElementAtSLat::position)
      .def_ro("ix_ele", &Bmad::ElementAtSLat::ix_ele)
      .def("__len__", [](const Bmad::ElementAtSLat &) { return 4; })
      .def("__getitem__", [](const Bmad::ElementAtSLat &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.s_eff);
        if (i == 2)
          return nb::cast(s.position);
        if (i == 3)
          return nb::cast(s.ix_ele);
        throw nb::index_error();
      });
  m.def(
      "element_at_s",
      nb::overload_cast<LatStruct &, double, bool, std::optional<int>, std::optional<bool>>(
          &Bmad::element_at_s
      ),
      nb::arg("lat"),
      nb::arg("s"),
      nb::arg("choose_max"),
      nb::arg("ix_branch") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Function to return the index of the element at position s.

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
ix_ele : int
    Index of element at s.

err_flag : bool, optional
    Set True if s is out of bounds. False otherwise.

s_eff : float, optional
    Effective s. Equal to s with a open lattice. See above.

position : CoordStruct, optional
    Positional information.
)"""
  );
  m.def(
      "element_slice_iterator",
      &Bmad::element_slice_iterator,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("i_slice"),
      nb::arg("n_slice_tot"),
      nb::arg("sliced_ele"),
      nb::arg("s_start") = nb::none(),
      nb::arg("s_end") = nb::none(),
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
  nb::class_<Bmad::EmFieldCalc>(m, "EmFieldCalc", "em_field_calc return type")
      .def_ro("field", &Bmad::EmFieldCalc::field)
      .def_ro("err_flag", &Bmad::EmFieldCalc::err_flag)
      .def("__len__", [](const Bmad::EmFieldCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::EmFieldCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.field);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "em_field_calc",
      [](EleStruct &ele,
         LatParamStruct &param,
         double s_pos,
         CoordStruct &orbit,
         bool local_ref_frame,
         std::optional<bool> calc_dfield,
         std::optional<bool> calc_potential,
         std::optional<bool> use_overlap,
         std::optional<bool> grid_allow_s_out_of_bounds,
         std::optional<double> rf_time,
         std::optional<ElePointerStructAlloc1D> used_eles,
         std::optional<bool> print_err,
         EleStruct *original_ele) {
        auto fn = static_cast<Bmad::EmFieldCalc (*)(
            EleStruct &,
            LatParamStruct &,
            double,
            CoordStruct &,
            bool,
            std::optional<bool>,
            std::optional<bool>,
            std::optional<bool>,
            std::optional<bool>,
            std::optional<double>,
            std::optional<ElePointerStructAlloc1D>,
            std::optional<bool>,
            optional_ref<EleStruct>
        )>(&Bmad::em_field_calc);
        return fn(
            ele,
            param,
            s_pos,
            orbit,
            local_ref_frame,
            calc_dfield,
            calc_potential,
            use_overlap,
            grid_allow_s_out_of_bounds,
            rf_time,
            used_eles,
            print_err,
            ptr_to_opt_ref(original_ele)
        );
      },
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("s_pos"),
      nb::arg("orbit"),
      nb::arg("local_ref_frame"),
      nb::arg("calc_dfield") = nb::none(),
      nb::arg("calc_potential") = nb::none(),
      nb::arg("use_overlap") = nb::none(),
      nb::arg("grid_allow_s_out_of_bounds") = nb::none(),
      nb::arg("rf_time") = nb::none(),
      nb::arg("used_eles") = nb::none(),
      nb::arg("print_err") = nb::none(),
      nb::arg("original_ele") = nb::none(),
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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("s_pos"),
      nb::arg("orbit"),
      nb::arg("local_ref_frame"),
      nb::arg("grid_allow_s_out_of_bounds") = nb::none(),
      nb::arg("rf_time") = nb::none(),
      R"""(Routine to calculate field derivatives.
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
      [](EleStruct &ele,
         LatParamStruct &param,
         double rf_time,
         CoordStruct &orbit,
         bool err_flag,
         std::optional<bool> print_err,
         EmFieldStruct *extra_field) {
        auto fn = static_cast<FixedArray1D<Real, 10> (*)(
            EleStruct &,
            LatParamStruct &,
            double,
            CoordStruct &,
            bool,
            std::optional<bool>,
            optional_ref<EmFieldStruct>
        )>(&Bmad::em_field_kick_vector_time);
        return fn(ele, param, rf_time, orbit, err_flag, print_err, ptr_to_opt_ref(extra_field));
      },
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("rf_time"),
      nb::arg("orbit"),
      nb::arg("err_flag"),
      nb::arg("print_err") = nb::none(),
      nb::arg("extra_field") = nb::none(),
      R"""(Subroutine to convert particle coordinates from t-based to s-based system.

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
      nb::arg("field1"),
      nb::arg("field2"),
      R"""(Wrapper for Fortran routine em_field_plus_em_field

Parameters
----------
field1 : EmFieldStruct

field2 : EmFieldStruct

Returns
-------
field_tot : EmFieldStruct
)"""
  );
  nb::class_<Bmad::Emit6d>(m, "Emit6d", "emit_6d return type")
      .def_ro("mode", &Bmad::Emit6d::mode)
      .def_ro("sigma_mat", &Bmad::Emit6d::sigma_mat)
      .def_ro("rad_int_by_ele", &Bmad::Emit6d::rad_int_by_ele)
      .def("__len__", [](const Bmad::Emit6d &) { return 3; })
      .def("__getitem__", [](const Bmad::Emit6d &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.mode);
        if (i == 1)
          return nb::cast(s.sigma_mat);
        if (i == 2)
          return nb::cast(s.rad_int_by_ele);
        throw nb::index_error();
      });
  m.def(
      "emit_6d",
      &Bmad::emit_6d,
      nb::arg("ele_ref"),
      nb::arg("include_opening_angle"),
      nb::arg("closed_orbit") = nb::none(),
      R"""(Routine to calculate the three normal mode emittances, damping partition numbers, radiation integrals, etc.
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
      "energy_func",
      &Bmad::energy_func,
      nb::arg("integ_prob"),
      nb::arg("status"),
      R"""(Wrapper for Fortran routine energy_func

Parameters
----------
integ_prob : float

status : int

Returns
-------
de : float
)"""
  );
  m.def(
      "entering_element",
      &Bmad::entering_element,
      nb::arg("orbit"),
      nb::arg("particle_at"),
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
      nb::arg("Lambda"),
      nb::arg("Theta"),
      nb::arg("Iota"),
      nb::arg("alpha"),
      nb::arg("emit"),
      R"""(Calculates damping decrement and emittance of the three
normal modes from the integrate diffusion, damping, and vertical
excitation matrices names Lambda, Theta, and Iota, respectively.
These three matrices are obtained from the subroutine integrated_mats.

The damping times can obtained from alpha using:
   tau = lattice_length/c_light/alpha
)"""
  );
  nb::class_<Bmad::EnvelopeRadintsIbs>(m, "EnvelopeRadintsIbs", "envelope_radints_ibs return type")
      .def_ro("alpha", &Bmad::EnvelopeRadintsIbs::alpha)
      .def_ro("emit", &Bmad::EnvelopeRadintsIbs::emit)
      .def("__len__", [](const Bmad::EnvelopeRadintsIbs &) { return 2; })
      .def("__getitem__", [](const Bmad::EnvelopeRadintsIbs &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.alpha);
        if (i == 1)
          return nb::cast(s.emit);
        throw nb::index_error();
      });
  m.def(
      "envelope_radints_ibs",
      &Bmad::envelope_radints_ibs,
      nb::arg("Lambda"),
      nb::arg("Theta"),
      nb::arg("Iota"),
      nb::arg("eles"),
      nb::arg("mode"),
      nb::arg("tail_cut"),
      nb::arg("npart"),
      nb::arg("species"),
      R"""(Calculates damping decrement and emittance of the three
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
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_ac_kicker

Parameters
----------
f1 : AcKickerStruct

f2 : AcKickerStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_ac_kicker_freq",
      &Bmad::eq_ac_kicker_freq,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_ac_kicker_freq

Parameters
----------
f1 : AcKickerFreqStruct

f2 : AcKickerFreqStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_ac_kicker_time",
      &Bmad::eq_ac_kicker_time,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_ac_kicker_time

Parameters
----------
f1 : AcKickerTimeStruct

f2 : AcKickerTimeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_anormal_mode",
      &Bmad::eq_anormal_mode,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_anormal_mode

Parameters
----------
f1 : AnormalModeStruct

f2 : AnormalModeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_aperture_param",
      &Bmad::eq_aperture_param,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_aperture_param

Parameters
----------
f1 : ApertureParamStruct

f2 : ApertureParamStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_aperture_point",
      &Bmad::eq_aperture_point,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_aperture_point

Parameters
----------
f1 : AperturePointStruct

f2 : AperturePointStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_aperture_scan",
      &Bmad::eq_aperture_scan,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_aperture_scan

Parameters
----------
f1 : ApertureScanStruct

f2 : ApertureScanStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_beam",
      &Bmad::eq_beam,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_beam

Parameters
----------
f1 : BeamStruct

f2 : BeamStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_beam_init",
      &Bmad::eq_beam_init,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_beam_init

Parameters
----------
f1 : BeamInitStruct

f2 : BeamInitStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_bmad_common",
      &Bmad::eq_bmad_common,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_bmad_common

Parameters
----------
f1 : BmadCommonStruct

f2 : BmadCommonStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_bookkeeping_state",
      &Bmad::eq_bookkeeping_state,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_bookkeeping_state

Parameters
----------
f1 : BookkeepingStateStruct

f2 : BookkeepingStateStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_bpm_phase_coupling",
      &Bmad::eq_bpm_phase_coupling,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_bpm_phase_coupling

Parameters
----------
f1 : BpmPhaseCouplingStruct

f2 : BpmPhaseCouplingStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_branch",
      &Bmad::eq_branch,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_branch

Parameters
----------
f1 : BranchStruct

f2 : BranchStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_bunch",
      &Bmad::eq_bunch,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_bunch

Parameters
----------
f1 : BunchStruct

f2 : BunchStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_bunch_params",
      &Bmad::eq_bunch_params,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_bunch_params

Parameters
----------
f1 : BunchParamsStruct

f2 : BunchParamsStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_cartesian_map",
      &Bmad::eq_cartesian_map,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_cartesian_map

Parameters
----------
f1 : CartesianMapStruct

f2 : CartesianMapStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_cartesian_map_term",
      &Bmad::eq_cartesian_map_term,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_cartesian_map_term

Parameters
----------
f1 : CartesianMapTermStruct

f2 : CartesianMapTermStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_cartesian_map_term1",
      &Bmad::eq_cartesian_map_term1,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_cartesian_map_term1

Parameters
----------
f1 : CartesianMapTerm1Struct

f2 : CartesianMapTerm1Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_complex_taylor",
      &Bmad::eq_complex_taylor,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_complex_taylor

Parameters
----------
f1 : ComplexTaylorStruct

f2 : ComplexTaylorStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_complex_taylor_term",
      &Bmad::eq_complex_taylor_term,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_complex_taylor_term

Parameters
----------
f1 : ComplexTaylorTermStruct

f2 : ComplexTaylorTermStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_control",
      &Bmad::eq_control,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_control

Parameters
----------
f1 : ControlStruct

f2 : ControlStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_control_ramp1",
      &Bmad::eq_control_ramp1,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_control_ramp1

Parameters
----------
f1 : ControlRamp1Struct

f2 : ControlRamp1Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_control_var1",
      &Bmad::eq_control_var1,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_control_var1

Parameters
----------
f1 : ControlVar1Struct

f2 : ControlVar1Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_controller",
      &Bmad::eq_controller,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_controller

Parameters
----------
f1 : ControllerStruct

f2 : ControllerStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_coord",
      &Bmad::eq_coord,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_coord

Parameters
----------
f1 : CoordStruct

f2 : CoordStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_coord_array",
      &Bmad::eq_coord_array,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_coord_array

Parameters
----------
f1 : CoordArrayStruct

f2 : CoordArrayStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_cylindrical_map",
      &Bmad::eq_cylindrical_map,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_cylindrical_map

Parameters
----------
f1 : CylindricalMapStruct

f2 : CylindricalMapStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_cylindrical_map_term",
      &Bmad::eq_cylindrical_map_term,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_cylindrical_map_term

Parameters
----------
f1 : CylindricalMapTermStruct

f2 : CylindricalMapTermStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_cylindrical_map_term1",
      &Bmad::eq_cylindrical_map_term1,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_cylindrical_map_term1

Parameters
----------
f1 : CylindricalMapTerm1Struct

f2 : CylindricalMapTerm1Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_ele",
      &Bmad::eq_ele,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_ele

Parameters
----------
f1 : EleStruct

f2 : EleStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_ellipse_beam_init",
      &Bmad::eq_ellipse_beam_init,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_ellipse_beam_init

Parameters
----------
f1 : EllipseBeamInitStruct

f2 : EllipseBeamInitStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_em_field",
      &Bmad::eq_em_field,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_em_field

Parameters
----------
f1 : EmFieldStruct

f2 : EmFieldStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_expression_atom",
      &Bmad::eq_expression_atom,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_expression_atom

Parameters
----------
f1 : ExpressionAtomStruct

f2 : ExpressionAtomStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_floor_position",
      &Bmad::eq_floor_position,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_floor_position

Parameters
----------
f1 : FloorPositionStruct

f2 : FloorPositionStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_gen_grad_curve",
      &Bmad::eq_gen_grad_curve,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_gen_grad_curve

Parameters
----------
f1 : GenGradCurveStruct

f2 : GenGradCurveStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_gen_gradients",
      &Bmad::eq_gen_gradients,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_gen_gradients

Parameters
----------
f1 : GenGradientsStruct

f2 : GenGradientsStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_gg_taylor",
      &Bmad::eq_gg_taylor,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_gg_taylor

Parameters
----------
f1 : GgTaylorStruct

f2 : GgTaylorStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_gg_taylor_term",
      &Bmad::eq_gg_taylor_term,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_gg_taylor_term

Parameters
----------
f1 : GgTaylorTermStruct

f2 : GgTaylorTermStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_grid_beam_init",
      &Bmad::eq_grid_beam_init,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_grid_beam_init

Parameters
----------
f1 : GridBeamInitStruct

f2 : GridBeamInitStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_grid_field",
      &Bmad::eq_grid_field,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_grid_field

Parameters
----------
f1 : GridFieldStruct

f2 : GridFieldStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_grid_field_pt",
      &Bmad::eq_grid_field_pt,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_grid_field_pt

Parameters
----------
f1 : GridFieldPtStruct

f2 : GridFieldPtStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_grid_field_pt1",
      &Bmad::eq_grid_field_pt1,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_grid_field_pt1

Parameters
----------
f1 : GridFieldPt1Struct

f2 : GridFieldPt1Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_high_energy_space_charge",
      &Bmad::eq_high_energy_space_charge,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_high_energy_space_charge

Parameters
----------
f1 : HighEnergySpaceChargeStruct

f2 : HighEnergySpaceChargeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_interval1_coef",
      &Bmad::eq_interval1_coef,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_interval1_coef

Parameters
----------
f1 : Interval1CoefStruct

f2 : Interval1CoefStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_kv_beam_init",
      &Bmad::eq_kv_beam_init,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_kv_beam_init

Parameters
----------
f1 : KvBeamInitStruct

f2 : KvBeamInitStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_lat",
      &Bmad::eq_lat,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_lat

Parameters
----------
f1 : LatStruct

f2 : LatStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_lat_ele_loc",
      &Bmad::eq_lat_ele_loc,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_lat_ele_loc

Parameters
----------
f1 : LatEleLocStruct

f2 : LatEleLocStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_lat_param",
      &Bmad::eq_lat_param,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_lat_param

Parameters
----------
f1 : LatParamStruct

f2 : LatParamStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_linac_normal_mode",
      &Bmad::eq_linac_normal_mode,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_linac_normal_mode

Parameters
----------
f1 : LinacNormalModeStruct

f2 : LinacNormalModeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_mode3",
      &Bmad::eq_mode3,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_mode3

Parameters
----------
f1 : Mode3Struct

f2 : Mode3Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_mode_info",
      &Bmad::eq_mode_info,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_mode_info

Parameters
----------
f1 : ModeInfoStruct

f2 : ModeInfoStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_normal_modes",
      &Bmad::eq_normal_modes,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_normal_modes

Parameters
----------
f1 : NormalModesStruct

f2 : NormalModesStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_photon_element",
      &Bmad::eq_photon_element,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_photon_element

Parameters
----------
f1 : PhotonElementStruct

f2 : PhotonElementStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_photon_material",
      &Bmad::eq_photon_material,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_photon_material

Parameters
----------
f1 : PhotonMaterialStruct

f2 : PhotonMaterialStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_photon_reflect_surface",
      &Bmad::eq_photon_reflect_surface,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_photon_reflect_surface

Parameters
----------
f1 : PhotonReflectSurfaceStruct

f2 : PhotonReflectSurfaceStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_photon_reflect_table",
      &Bmad::eq_photon_reflect_table,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_photon_reflect_table

Parameters
----------
f1 : PhotonReflectTableStruct

f2 : PhotonReflectTableStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_photon_target",
      &Bmad::eq_photon_target,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_photon_target

Parameters
----------
f1 : PhotonTargetStruct

f2 : PhotonTargetStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_pixel_detec",
      &Bmad::eq_pixel_detec,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_pixel_detec

Parameters
----------
f1 : PixelDetecStruct

f2 : PixelDetecStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_pixel_pt",
      &Bmad::eq_pixel_pt,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_pixel_pt

Parameters
----------
f1 : PixelPtStruct

f2 : PixelPtStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_pre_tracker",
      &Bmad::eq_pre_tracker,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_pre_tracker

Parameters
----------
f1 : PreTrackerStruct

f2 : PreTrackerStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_rad_int1",
      &Bmad::eq_rad_int1,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_rad_int1

Parameters
----------
f1 : RadInt1Struct

f2 : RadInt1Struct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_rad_int_all_ele",
      &Bmad::eq_rad_int_all_ele,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_rad_int_all_ele

Parameters
----------
f1 : RadIntAllEleStruct

f2 : RadIntAllEleStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_rad_int_branch",
      &Bmad::eq_rad_int_branch,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_rad_int_branch

Parameters
----------
f1 : RadIntBranchStruct

f2 : RadIntBranchStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_rad_map",
      &Bmad::eq_rad_map,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_rad_map

Parameters
----------
f1 : RadMapStruct

f2 : RadMapStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_rad_map_ele",
      &Bmad::eq_rad_map_ele,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_rad_map_ele

Parameters
----------
f1 : RadMapEleStruct

f2 : RadMapEleStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_ramper_lord",
      &Bmad::eq_ramper_lord,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_ramper_lord

Parameters
----------
f1 : RamperLordStruct

f2 : RamperLordStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_space_charge_common",
      &Bmad::eq_space_charge_common,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_space_charge_common

Parameters
----------
f1 : SpaceChargeCommonStruct

f2 : SpaceChargeCommonStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_spin_polar",
      &Bmad::eq_spin_polar,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_spin_polar

Parameters
----------
f1 : SpinPolarStruct

f2 : SpinPolarStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_spline",
      &Bmad::eq_spline,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_spline

Parameters
----------
f1 : SplineStruct

f2 : SplineStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_strong_beam",
      &Bmad::eq_strong_beam,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_strong_beam

Parameters
----------
f1 : StrongBeamStruct

f2 : StrongBeamStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_curvature",
      &Bmad::eq_surface_curvature,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_curvature

Parameters
----------
f1 : SurfaceCurvatureStruct

f2 : SurfaceCurvatureStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_displacement",
      &Bmad::eq_surface_displacement,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_displacement

Parameters
----------
f1 : SurfaceDisplacementStruct

f2 : SurfaceDisplacementStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_displacement_pt",
      &Bmad::eq_surface_displacement_pt,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_displacement_pt

Parameters
----------
f1 : SurfaceDisplacementPtStruct

f2 : SurfaceDisplacementPtStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_h_misalign",
      &Bmad::eq_surface_h_misalign,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_h_misalign

Parameters
----------
f1 : SurfaceHMisalignStruct

f2 : SurfaceHMisalignStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_h_misalign_pt",
      &Bmad::eq_surface_h_misalign_pt,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_h_misalign_pt

Parameters
----------
f1 : SurfaceHMisalignPtStruct

f2 : SurfaceHMisalignPtStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_segmented",
      &Bmad::eq_surface_segmented,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_segmented

Parameters
----------
f1 : SurfaceSegmentedStruct

f2 : SurfaceSegmentedStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_surface_segmented_pt",
      &Bmad::eq_surface_segmented_pt,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_surface_segmented_pt

Parameters
----------
f1 : SurfaceSegmentedPtStruct

f2 : SurfaceSegmentedPtStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_target_point",
      &Bmad::eq_target_point,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_target_point

Parameters
----------
f1 : TargetPointStruct

f2 : TargetPointStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_taylor",
      &Bmad::eq_taylor,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_taylor

Parameters
----------
f1 : TaylorStruct

f2 : TaylorStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_taylor_term",
      &Bmad::eq_taylor_term,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_taylor_term

Parameters
----------
f1 : TaylorTermStruct

f2 : TaylorTermStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_track",
      &Bmad::eq_track,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_track

Parameters
----------
f1 : TrackStruct

f2 : TrackStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_track_point",
      &Bmad::eq_track_point,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_track_point

Parameters
----------
f1 : TrackPointStruct

f2 : TrackPointStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_twiss",
      &Bmad::eq_twiss,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_twiss

Parameters
----------
f1 : TwissStruct

f2 : TwissStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wake",
      &Bmad::eq_wake,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wake

Parameters
----------
f1 : WakeStruct

f2 : WakeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wake_lr",
      &Bmad::eq_wake_lr,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wake_lr

Parameters
----------
f1 : WakeLrStruct

f2 : WakeLrStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wake_lr_mode",
      &Bmad::eq_wake_lr_mode,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wake_lr_mode

Parameters
----------
f1 : WakeLrModeStruct

f2 : WakeLrModeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wake_sr",
      &Bmad::eq_wake_sr,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wake_sr

Parameters
----------
f1 : WakeSrStruct

f2 : WakeSrStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wake_sr_mode",
      &Bmad::eq_wake_sr_mode,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wake_sr_mode

Parameters
----------
f1 : WakeSrModeStruct

f2 : WakeSrModeStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wake_sr_z_long",
      &Bmad::eq_wake_sr_z_long,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wake_sr_z_long

Parameters
----------
f1 : WakeSrZLongStruct

f2 : WakeSrZLongStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wall3d",
      &Bmad::eq_wall3d,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wall3d

Parameters
----------
f1 : Wall3dStruct

f2 : Wall3dStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wall3d_section",
      &Bmad::eq_wall3d_section,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wall3d_section

Parameters
----------
f1 : Wall3dSectionStruct

f2 : Wall3dSectionStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_wall3d_vertex",
      &Bmad::eq_wall3d_vertex,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_wall3d_vertex

Parameters
----------
f1 : Wall3dVertexStruct

f2 : Wall3dVertexStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "eq_xy_disp",
      &Bmad::eq_xy_disp,
      nb::arg("f1"),
      nb::arg("f2"),
      R"""(Wrapper for Fortran routine eq_xy_disp

Parameters
----------
f1 : XyDispStruct

f2 : XyDispStruct

Returns
-------
is_eq : bool
)"""
  );
  m.def(
      "equal_sign_here",
      &Bmad::equal_sign_here,
      nb::arg("ele"),
      nb::arg("delim"),
      R"""(Wrapper for Fortran routine equal_sign_here

Parameters
----------
ele : EleStruct

delim : str

Returns
-------
is_here : bool
)"""
  );
  m.def(
      "equivalent_taylor_attributes",
      &Bmad::equivalent_taylor_attributes,
      nb::arg("ele_taylor"),
      nb::arg("ele2"),
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
      nb::arg("A"),
      nb::arg("B"),
      nb::arg("C"),
      nb::arg("D"),
      nb::arg("E"),
      nb::arg("F"),
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
  nb::class_<Bmad::EvaluateArrayIndex>(m, "EvaluateArrayIndex", "evaluate_array_index return type")
      .def_ro("err_flag", &Bmad::EvaluateArrayIndex::err_flag)
      .def_ro("word2", &Bmad::EvaluateArrayIndex::word2)
      .def_ro("delim2", &Bmad::EvaluateArrayIndex::delim2)
      .def_ro("this_index", &Bmad::EvaluateArrayIndex::this_index)
      .def("__len__", [](const Bmad::EvaluateArrayIndex &) { return 4; })
      .def("__getitem__", [](const Bmad::EvaluateArrayIndex &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.word2);
        if (i == 2)
          return nb::cast(s.delim2);
        if (i == 3)
          return nb::cast(s.this_index);
        throw nb::index_error();
      });
  m.def(
      "evaluate_array_index",
      &Bmad::evaluate_array_index,
      nb::arg("delim_list1"),
      nb::arg("delim_list2"),
      R"""(Function of evaluate the index of an array. Typically the text being parsed looks like:
     "5) = ..."         or
     "6).COMP = ..."

Parameters
----------
delim_list1 : str
    Delimitor after the integer. Normally ')'.

delim_list2 : str
    Delimitor list to mark the end of word2. Normally '='.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.

word2 : str
    Word found after delim1. Normally this should be blank.

delim2 : str
    Actual delimitor found after word2.

this_index : int
    Integer value
)"""
  );
  nb::class_<Bmad::EvaluateLogical>(m, "EvaluateLogical", "evaluate_logical return type")
      .def_ro("iostat", &Bmad::EvaluateLogical::iostat)
      .def_ro("this_logic", &Bmad::EvaluateLogical::this_logic)
      .def("__len__", [](const Bmad::EvaluateLogical &) { return 2; })
      .def("__getitem__", [](const Bmad::EvaluateLogical &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.iostat);
        if (i == 1)
          return nb::cast(s.this_logic);
        throw nb::index_error();
      });
  m.def(
      "evaluate_logical",
      &Bmad::evaluate_logical,
      nb::arg("word"),
      R"""(Function of convert a string into a logical value.
Accepted possibilities are:
  .TRUE.  .FALSE.
   TRUE    FALSE
   T       F

Parameters
----------
word : str
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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("particle_at"),
      nb::arg("orb"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Subroutine to track through the edge field of an sbend.
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
      nb::arg("t"),
      nb::arg("B1"),
      nb::arg("B2"),
      R"""(This is essentially the Numercal Recipes bessi0 function multiplied by exp(-B1*t).

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
  nb::class_<PyExpectOneOf>(m, "ExpectOneOf", "expect_one_of return type")
      .def_ro("is_ok", &PyExpectOneOf::is_ok)
      .def_ro("delim", &PyExpectOneOf::delim)
      .def("__len__", [](const PyExpectOneOf &) { return 2; })
      .def("__getitem__", [](const PyExpectOneOf &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.is_ok);
        if (i == 1)
          return nb::cast(s.delim);
        throw nb::index_error();
      });
  m.def(
      "expect_one_of",
      &python_expect_one_of,
      nb::arg("delim_list"),
      nb::arg("check_input_delim"),
      nb::arg("ele_name"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      R"""(Routine to check either that the current delimitor or the next character in the parse stream is the
expected delimitor.
This routine is used for Bmad lattice file parsing and is not meant for general use.

Also see: expect_this

Parameters
----------
delim_list : str
    List of expected (valid) delimitors. If list contains a space character then no delimitor (indicating the
    end of the command) is a valid possibility.

check_input_delim : bool
    If True, then check if delim argument is in the delim_list. If False, check that the next character in the
    parse stream is an expected delimitor.

ele_name : str
    Lattice element under construction. Used for error messages.

delim : str
    Current delimitor that will be checked if check_input_delim = .true.
    This parameter is an input/output and is modified in-place.
    As an output, delim: Next delim if check_input_delim = False.

Returns
-------
delim : str
    Current delimitor that will be checked if check_input_delim = .true.
    This parameter is an input/output and is modified in-place.
    As an output, delim: Next delim if check_input_delim = False.
)"""
  );
  nb::class_<Bmad::ExpectThis>(m, "ExpectThis", "expect_this return type")
      .def_ro("delim", &Bmad::ExpectThis::delim)
      .def_ro("delim_found", &Bmad::ExpectThis::delim_found)
      .def_ro("is_ok", &Bmad::ExpectThis::is_ok)
      .def("__len__", [](const Bmad::ExpectThis &) { return 3; })
      .def("__getitem__", [](const Bmad::ExpectThis &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.delim);
        if (i == 1)
          return nb::cast(s.delim_found);
        if (i == 2)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "expect_this",
      &Bmad::expect_this,
      nb::arg("expecting"),
      nb::arg("check_delim"),
      nb::arg("call_check"),
      nb::arg("err_str"),
      nb::arg("ele"),
      R"""(Checks that the next character or characters in the parse stream corresponds to the
characters in the expecting argument. For example, if expecting is ')={' these three characters
should be the next non-blank characters in the parse stream.

Also see: expect_one_of

Parameters
----------
expecting : str
    list of characters that are expected to be next in the parse stream.

check_delim : bool
    If True then use delim argument as first token to check. A blank character indicates end of command is
    expected.

call_check : bool
    If True then check for 'call::<filename>' construct.

err_str : str
    String used for error messages.

ele : EleStruct
    Element parameters being parsed.

Returns
-------
delim : str
    Final delim

delim_found : bool
    Is there a final delim (as opposed to end of command).
)"""
  );
  m.def(
      "expression_stack_to_string",
      &Bmad::expression_stack_to_string,
      nb::arg("stack"),
      nb::arg("polish") = nb::none(),
      R"""(Routine to convert an expression stack to a string

Parameters
----------
stack : 1D array of ExpressionAtomStruct
    arithmetic expression

polish : bool, optional
    logical, optional, Construct expression in reverse polish? Default is False.

Returns
-------
str : str
    : Expression in string form.
)"""
  );
  nb::class_<Bmad::ExpressionStackValue>(
      m,
      "ExpressionStackValue",
      "expression_stack_value return type"
  )
      .def_ro("err_flag", &Bmad::ExpressionStackValue::err_flag)
      .def_ro("err_str", &Bmad::ExpressionStackValue::err_str)
      .def_ro("value", &Bmad::ExpressionStackValue::value)
      .def("__len__", [](const Bmad::ExpressionStackValue &) { return 3; })
      .def("__getitem__", [](const Bmad::ExpressionStackValue &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.err_str);
        if (i == 2)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "expression_stack_value",
      &Bmad::expression_stack_value,
      nb::arg("stack"),
      nb::arg("var") = nb::none(),
      nb::arg("use_old") = nb::none(),
      R"""(Routine to evaluate a mathematical expression represented by an "expression stack".
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

err_str : str
    Error string explaining error if there is one.

value : float
    Value of the expression.
)"""
  );
  nb::class_<Bmad::ExpressionStringToStack>(
      m,
      "ExpressionStringToStack",
      "expression_string_to_stack return type"
  )
      .def_ro("stack", &Bmad::ExpressionStringToStack::stack)
      .def_ro("n_stack", &Bmad::ExpressionStringToStack::n_stack)
      .def_ro("err_flag", &Bmad::ExpressionStringToStack::err_flag)
      .def_ro("err_str", &Bmad::ExpressionStringToStack::err_str)
      .def("__len__", [](const Bmad::ExpressionStringToStack &) { return 4; })
      .def("__getitem__", [](const Bmad::ExpressionStringToStack &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.stack);
        if (i == 1)
          return nb::cast(s.n_stack);
        if (i == 2)
          return nb::cast(s.err_flag);
        if (i == 3)
          return nb::cast(s.err_str);
        throw nb::index_error();
      });
  m.def(
      "expression_string_to_stack",
      &Bmad::expression_string_to_stack,
      nb::arg("string"),
      R"""(This routine creates an expression stack array which can be used
to evaluate an arithmethic expression.

Stack end elements not used are marked stack(i)%type = end_stack$

Stack elements with stack(i)%type = variable$ are elements that need
to be evaluated before calling expression_stack_value.

Also see:
  expression_value
  expression_stack_value

Parameters
----------
string : str
    Expression to be converted.

Returns
-------
stack : 1D array of ExpressionAtomStruct
    Expression evaluation stack.

n_stack : int
    number of "atoms" used by the expression

err_flag : bool
    Set True if there is an error (EG divide by 0).

err_str : str
    String describing the error.
)"""
  );
  nb::class_<Bmad::ExpressionStringToTree>(
      m,
      "ExpressionStringToTree",
      "expression_string_to_tree return type"
  )
      .def_ro("err_flag", &Bmad::ExpressionStringToTree::err_flag)
      .def_ro("err_str", &Bmad::ExpressionStringToTree::err_str)
      .def("__len__", [](const Bmad::ExpressionStringToTree &) { return 2; })
      .def("__getitem__", [](const Bmad::ExpressionStringToTree &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.err_str);
        throw nb::index_error();
      });
  m.def(
      "expression_string_to_tree",
      &Bmad::expression_string_to_tree,
      nb::arg("string"),
      nb::arg("root_tree"),
      R"""(Routine to create an expression tree array which can be used
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
  operators: + -  * / ^ = : &
  brackets: [] () {}
  comma: ,

Reverse polish tree will be constructed for mathematical expressions.
Exception: Any expression with the "->" (arrow$) token is not considered mathematical.
The reason why "::" is still considered mathematical is due to Tao syntax.

Root node name is "root" and is of type root$
Brackets in the expression string must be matched.
The corresponding tree node will have a name / type of:
  "[]" / square_brackets$,    "()" / parens$ or func_parens$,   "{}" / curley_brackets$

The root node, equal nodes, and all bracket nodes, will have an array of child nodes all of which will be comma nodes.

Parameters
----------
string : str
    Expression to be converted.

root_tree : ExpressionTreeStruct
    Only used when recursively called.

Returns
-------
err_flag : bool
    Set True if there is an error (EG divide by 0).

err_str : str
    String describing the error. Make length large to hold the expression.
)"""
  );
  m.def(
      "expression_tree_to_string",
      [](ExpressionTreeStruct &tree,
         std::optional<bool> include_root,
         std::optional<int> n_node,
         ExpressionTreeStruct *parent) {
        auto fn = static_cast<std::string (*)(
            ExpressionTreeStruct &,
            std::optional<bool>,
            std::optional<int>,
            optional_ref<ExpressionTreeStruct>
        )>(&Bmad::expression_tree_to_string);
        return fn(tree, include_root, n_node, ptr_to_opt_ref(parent));
      },
      nb::arg("tree"),
      nb::arg("include_root") = nb::none(),
      nb::arg("n_node") = nb::none(),
      nb::arg("parent") = nb::none(),
      R"""(Routine to convert an expression tree to a expression string.

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
str_out : str
    Expression string.
)"""
  );
  nb::class_<Bmad::ExpressionValue>(m, "ExpressionValue", "expression_value return type")
      .def_ro("err_flag", &Bmad::ExpressionValue::err_flag)
      .def_ro("err_str", &Bmad::ExpressionValue::err_str)
      .def_ro("value", &Bmad::ExpressionValue::value)
      .def("__len__", [](const Bmad::ExpressionValue &) { return 3; })
      .def("__getitem__", [](const Bmad::ExpressionValue &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.err_str);
        if (i == 2)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "expression_value",
      &Bmad::expression_value,
      nb::arg("expression"),
      nb::arg("var") = nb::none(),
      nb::arg("use_old") = nb::none(),
      R"""(Routine to evaluate a mathematical expression encoded in a string.

Also see:
  expression_string_to_stack
  expression_stack_value

Parameters
----------
expression : str
    Expression string.

var : 1D array of ControlVar1Struct, optional
    Array of control variables. Used with Bmad controller elements.

use_old : bool, optional
    Use var.old_value? Must be present if var(:) is present.

Returns
-------
err_flag : bool
    True if there is an evaluation problem. False otherwise.

value : float
    Value of the expression.

err_str : str, optional
    Error string explaining error if there is one.
)"""
  );
}
