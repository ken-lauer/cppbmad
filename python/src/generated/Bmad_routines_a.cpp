#include "pybmad/generated/Bmad_routines_a.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_a(py::module& m) {
  py::class_<Bmad::AbMultipoleKick, std::unique_ptr<Bmad::AbMultipoleKick>>(
      m, "AbMultipoleKick", "ab_multipole_kick return type")
      .def_readonly("kx", &Bmad::AbMultipoleKick::kx)
      .def_readonly("ky", &Bmad::AbMultipoleKick::ky)
      .def_readonly("dk", &Bmad::AbMultipoleKick::dk)
      .def("__len__", [](const Bmad::AbMultipoleKick&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::AbMultipoleKick& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.kx);
            if (i == 1)
              return py::cast(s.ky);
            if (i == 2)
              return py::cast(s.dk);
            throw py::index_error();
          });
  m.def(
      "ab_multipole_kick",
      &Bmad::ab_multipole_kick,
      py::arg("a"),
      py::arg("b"),
      py::arg("n"),
      py::arg("ref_species"),
      py::arg("ele_orientation"),
      py::arg("coord"),
      py::arg("pole_type") = py::none(),
      py::arg("scale") = py::none(),
      R"""(Subroutine ab_multipole_kick (a, b, n, ref_species, ele_orientation, coord, kx, ky, dk, pole_type, scale)

Subroutine to put in the kick due to an ab_multipole.

Parameters
----------
a : float
    Multipole skew component.
b : float
    Multipole normal component.
n : float
    Multipole order.
ref_species : int
    Reference species.
ele_orientation : int
    Element orientation +1 = normal, -1 = reversed, 0 = Ignore orientation and tracking species (used with
    pole_type = magnetic$).
coord : CoordStruct
    Particle position and direction of travel.
pole_type : int, optional
    Type of multipole. magnetic$ (default) or electric$.
scale : float, optional
    Factor to scale the kicks. Default is 1. For pole_type = electric$, set scale to the longitudinal length
    of the field region.

Returns
-------
kx : float
    X kick.
ky : float
    Y kick.
dk : float
    Kick derivative: dkick(x,y)/d(x,y).
)""");
  m.def(
      "ab_multipole_kicks",
      &Bmad::ab_multipole_kicks,
      py::arg("an"),
      py::arg("bn"),
      py::arg("ix_pole_max"),
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("pole_type") = py::none(),
      py::arg("scale") = py::none(),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine ab_multipole_kicks (an, bn, ix_pole_max, ele, orbit, pole_type, scale, mat6, make_matrix)

Routine to put in the kick due to ab_multipole components in an element.
The kick will be corrected for the orientation of the element and the particle direction of travel.
Any difference between element p0c and orbit%p0c will be taken into account.

Also see the multipole_kicks routine.

Parameters
----------
an : float
    Skew multipole strengths.
bn : float
    Normal multipole strengths.
ix_pole_max : int
    Maximum pole index.
ele : EleStruct
    Lattice element containing the multipoles.
orbit : CoordStruct
    Particle position.
    This parameter is an input/output and is modified in-place. As an output: Kicked particle.
pole_type : int, optional
    Type of multipole. magnetic$ (default) or electric$.
scale : float, optional
    Factor to scale the kicks. Default is 1. For pole_type = electric$, set scale to the longitudinal length
    of the field region
mat6 : float, optional
    Transfer matrix before the multipole.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix transfer matrix
    including multipole.
make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)""");
  m.def(
      "absolute_photon_position",
      &Bmad::absolute_photon_position,
      py::arg("e_orb"),
      py::arg("photon_orb"),
      R"""(Subroutine absolute_photon_position (e_orb, photon_orb)

Routine to calculate the photon phase space coordinates given:
  1) The phase space coords of the emitting charged particle and
  2) The photon phase space coords relative to the emitting particle.
     The photon (x, y, z) position is ignored (it is assumed the photon is emitted at
     the charged particle position) and only the photon's (vx, vy, vz) velocity matters.

Parameters
----------
e_orb : CoordStruct
    charged particle position.
photon_orb : CoordStruct
    Photon position relative to e_orb.
    This parameter is an input/output and is modified in-place. As an output: Absolute photon position.
)""");
  m.def(
      "absolute_time_tracking",
      &Bmad::absolute_time_tracking,
      py::arg("ele"),
      py::arg("is_abs_time"),
      R"""(Parameters
----------
ele : EleStruct
    Element being tracked through.
is_abs_time : 
)""");
  m.def(
      "ac_kicker_amp",
      &Bmad::ac_kicker_amp,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("true_time") = py::none(),
      py::arg("ac_amp"),
      R"""(Parameters
----------
ele : EleStruct
    ac_kicker element.
orbit : CoordStruct
    Contains the time to evaluate the amplitude at.
true_time : float, optional
    The actual time. Normally this time is calculated using orbit.t or orbit.vec(5) but sometimes it is
    convenient to be able to override this. For example, time_runge_kutta uses this.
ac_amp : 
)""");
  py::class_<Bmad::ActionToXyz, std::unique_ptr<Bmad::ActionToXyz>>(
      m, "ActionToXyz", "action_to_xyz return type")
      .def_readonly("X", &Bmad::ActionToXyz::X)
      .def_readonly("err_flag", &Bmad::ActionToXyz::err_flag)
      .def("__len__", [](const Bmad::ActionToXyz&) { return 2; })
      .def("__getitem__", [](const Bmad::ActionToXyz& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.X);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "action_to_xyz",
      &Bmad::action_to_xyz,
      py::arg("ring"),
      py::arg("ix"),
      py::arg("J"),
      R"""(Subroutine action_to_xyz(ring, ix, J, X, err_flag)

Given the normal mode invariants and phases J of a particle, returns the canonical coordinates.

The J vector looks like:
J = (sqrt(2Ja)cos(phia), -sqrt(2Ja)sin(phia), sqrt(2Jb)cos(phib), -sqrt(2Jb)sin(phib), sqrt(2Jc)cos(phic), -sqrt(2Jc)sin(phic))

X is obtained from:
X = N . J
Where N is from the Eigen decomposition of the 1-turn transfer matrix.

Parameters
----------
ring : LatStruct
    lattice .a.tune   -- a-mode tune (horizontal-like) .b.tune   -- b-mode tune (vertical-like) .z.tune   --
    c-mode tune (synchrotron-like)
ix : int
    element index at which to calculate J
J : float
    Vector containing normal mode invariants and phases

Returns
-------
X : float
    canonical phase space coordinates of the particle
err_flag : bool
    Set to true on error.  Often means Eigen decomposition failed.
)""");
  m.def(
      "add_lattice_control_structs",
      &Bmad::add_lattice_control_structs,
      py::arg("ele"),
      py::arg("n_add_slave") = py::none(),
      py::arg("n_add_lord") = py::none(),
      py::arg("n_add_slave_field") = py::none(),
      py::arg("n_add_lord_field") = py::none(),
      py::arg("add_at_end") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Lord or slave element that needs extra control elements.
n_add_slave : int, optional
    Number of field slaves to add to lord. Default is zero.
n_add_lord : int, optional
    Number of field lords to add to slave. Default is zero.
n_add_slave_field : int, optional
    Number of field slaves to add to lord. Default is zero.
n_add_lord_field : int, optional
    Number of field lords to add to slave. Default is zero.
add_at_end : bool, optional
    Used when n_add_slave or n_add_slave_field is non-zero. If True then new space is added at the end of the
    array. If False then new space is added at the front of the array. Default is True.
)""");
  py::class_<Bmad::AddSuperimpose, std::unique_ptr<Bmad::AddSuperimpose>>(
      m, "AddSuperimpose", "add_superimpose return type")
      .def_readonly("err_flag", &Bmad::AddSuperimpose::err_flag)
      .def_readonly("super_ele_out", &Bmad::AddSuperimpose::super_ele_out)
      .def("__len__", [](const Bmad::AddSuperimpose&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::AddSuperimpose& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.super_ele_out);
            throw py::index_error();
          });
  m.def(
      "add_superimpose",
      &Bmad::add_superimpose,
      py::arg("lat"),
      py::arg("super_ele_in"),
      py::arg("ix_branch"),
      py::arg("save_null_drift") = py::none(),
      py::arg("create_jumbo_slave") = py::none(),
      py::arg("ix_insert") = py::none(),
      py::arg("mangle_slave_names") = py::none(),
      py::arg("wrap") = py::none(),
      R"""(Subroutine add_superimpose (lat, super_ele_in, ix_branch, err_flag, super_ele_out,
               save_null_drift, create_jumbo_slave, ix_insert, mangle_slave_names, wrap)

Routine to superimpose an element. If the element can be inserted
into the lat without making a super_lord element then this will be done.

Note: This routine, since it handles only one superposition, is not sufficient for
  superposition in a multipass region. For historical reasons, the extra code needed
  is buried in the parser_add_superimpose code. If you need to do multipass superpositions
  please contact David Sagan and this situation will be rectified.

Note: Bookkeeping like recalculating reference energies and recalculating transfer matrices
  is *not* done by this routine.

Parameters
----------
lat : LatStruct
    Lat to modify.
    This parameter is an input/output and is modified in-place. As an output: Modified lat.
super_ele_in : EleStruct
    Element to superimpose. .s               -- Position of end of element. Negative distances mean distance
    from the end.
ix_branch : int
    Branch index to put element.
save_null_drift : bool, optional
    Save a copy of a drift to be split as a null_ele? This is useful if further superpositions might use this
    drift as a reference element. After all superpositions are done, remove_eles_from_lat can be called to
    remove all null_eles. Default is False.
create_jumbo_slave : bool, optional
    Default is False. If True then super_slaves that are created that have super_ele_in as their super_lord
    are em_field elements.
ix_insert : int, optional
    If present and positive, and super_ele_in has zero length, use ix_insert as the index to insert
    super_ele_in at. ix_insert is useful when superposing next to another element that has zero or negative
    length (EG a patch) and you want to make sure that the superimposed element is on the correct side of the
    element.
mangle_slave_names : bool, optional
    If True (default), adjust slave names appropriately. Name mangeling can take time so bmad_parser will do
    this all at once at the end.
wrap : bool, optional
    If True (default), and if the superimposed element has an end that extends beyond the starting or ending
    edge of the lattice, wrap the element around the lattice so that the beginning portion of the element is
    at the lattice ending edge and the rest of the element is at the lattice start edge. If wrap = False, and
    the superimposed element has an end that extends beyound a lattice edge, extend the lattice to
    accommodate.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise
super_ele_out : EleStruct
    Pointer to the super element in the lattice.
)""");
  m.def(
      "add_this_multipass",
      &Bmad::add_this_multipass,
      py::arg("lat"),
      py::arg("m_slaves"),
      py::arg("lord_in") = py::none(),
      R"""(Parameters
----------
lat : 
m_slaves : 
lord_in : 
)""");
  m.def(
      "add_this_taylor_term",
      &Bmad::add_this_taylor_term,
      py::arg("ele"),
      py::arg("i_out"),
      py::arg("coef"),
      py::arg("expn"),
      R"""(Subroutine add_this_taylor_term (ele, i_out, coef, expn)

Subroutine used by bmad_parser and bmad_parser2 to parse the input file.
This subroutine is not intended for general use.

)""");
  m.def(
      "adjust_super_slave_names",
      &Bmad::adjust_super_slave_names,
      py::arg("lat"),
      py::arg("ix1_lord"),
      py::arg("ix2_lord"),
      py::arg("first_time") = py::none(),
      R"""(Subroutine adjust_super_slave_names (lat, ix1_lord, ix2_lord, first_time)

Routine to adjust the names of the slaves.
This routine is used by add_superimpose and is not meant for general use.

)""");
  m.def(
      "allocate_branch_array",
      &Bmad::allocate_branch_array,
      py::arg("lat"),
      py::arg("upper_bound"),
      R"""(Parameters
----------
lat : LatStruct
    .branch(:)  -- Branch array to be allocated.
upper_bound : int
    Desired upper bound.
)""");
  m.def(
      "allocate_lat_ele_array",
      &Bmad::allocate_lat_ele_array,
      py::arg("lat"),
      py::arg("upper_bound") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("do_ramper_slave_setup") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Lattice with element array. .branch(ix_branch).ele(:)  -- Element array to reallocate.
upper_bound : int, optional
    Optional desired upper bound. Default: 1.3*ubound(ele(:)) or 10 if ele is not allocated.
ix_branch : int, optional
    Branch index. Default is 0.
do_ramper_slave_setup : bool, optional
    Default False. If true, setup ramper slaves. Generally this needs to be done if reallocating with a fully
    formed lattice.
)""");
  m.def(
      "angle_between_polars",
      &Bmad::angle_between_polars,
      py::arg("polar1"),
      py::arg("polar2"),
      py::arg("angle"),
      R"""(Parameters
----------
polar1 : 
    (spin_polar_struct)
polar2 : 
    (spin_polar_struct)
angle : 
)""");
  m.def(
      "angle_to_canonical_coords",
      &Bmad::angle_to_canonical_coords,
      py::arg("orbit"),
      py::arg("coord_type") = py::none(),
      R"""(Parameters
----------
orbit : CoordStruct
    Orbit in angular coordinates.
    This parameter is an input/output and is modified in-place. As an output: Orbit in canonical coordinates.
coord_type : unknown, optional
    Angular coordinates type '' (default): (x, x' = dx/ds, y, y' = dy/ds, z, pz) 'ZGOUBI':     (x, x' = dx/ds,
    y, y' = dy/ds, dt = -z / (beta * c), pz)
)""");
  m.def(
      "aperture_bookkeeper",
      &Bmad::aperture_bookkeeper,
      py::arg("ele"),
      R"""(Subroutine aperture_bookkeeper (ele)

Routine to calculate aperture limits when ele%attribute_type is set to auto_aperture$

Parameters
----------
ele : EleStruct
    Element with aperture.
    This parameter is an input/output and is modified in-place. As an output: Element with apertures set.
)""");
  m.def(
      "apply_all_rampers",
      &Bmad::apply_all_rampers,
      py::arg("lat"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice.
    This parameter is an input/output and is modified in-place. As an output: Lattice with rampers applied.
err_flag : bool
    Set True if there is an error. False otherwise.
)""");
  m.def(
      "apply_energy_kick",
      &Bmad::apply_energy_kick,
      py::arg("dE"),
      py::arg("orbit"),
      py::arg("ddE_dr"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Parameters
----------
dE : float
    Energy change
orbit : CoordStruct
    Beginning coordinates
    This parameter is an input/output and is modified in-place. As an output: coordinates with added dE energy
    kick.
ddE_dr : 
    real(rp), Derivatives of dE [ddE_dx, ddE_dy].
mat6 : float, optional
    Transfer matrix before fringe.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix transfer matrix
    including energy kick.
make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)""");
  m.def(
      "apply_patch_to_ptc_fibre",
      &Bmad::apply_patch_to_ptc_fibre,
      py::arg("ele"),
      R"""(Subroutine apply_patch_to_ptc_fibre (ele)

Routine to take the patch parameters from a Bmad patch element and
transfer them to the associated PTC fibre.

Parameters
----------
ele : EleStruct
    Patch element. ele.ptc_fibre -- PTC Fibre which should be a marker.
)""");
  m.def(
      "apply_rampers_to_slave",
      &Bmad::apply_rampers_to_slave,
      py::arg("slave"),
      R"""(Parameters
----------
slave : EleStruct
    Element to apply ramper elements to.
err_flag : bool
    Set true if there is an error. False otherwise.
)""");
  m.def(
      "array_re_str",
      &Bmad::array_re_str,
      py::arg("arr"),
      py::arg("parens_in") = py::none(),
      py::arg("str_out"),
      R"""(Parameters
----------
arr : 
parens_in : 
str_out : 
)""");
  m.def(
      "astra_max_field_reference",
      &Bmad::astra_max_field_reference,
      py::arg("pt0"),
      py::arg("ele"),
      py::arg("field_value"),
      R"""(Parameters
----------
pt0 : 
ele : 
field_value : 
)""");
  m.def(
      "at_this_ele_end",
      &Bmad::at_this_ele_end,
      py::arg("now_at"),
      py::arg("where_at"),
      py::arg("is_at_this_end"),
      R"""(Parameters
----------
now_at : int
    Which end is under consideration: entrance_end$, exit_end$, surface$, or in_between$.
where_at : int
    Which ends have the aperture or fringe field: entrance_end$, exit_end$, continuous$, both_ends$,
    no_aperture$, surface$, wall_transition$.
is_at_this_end : 
)""");
  m.def(
      "attribute_bookkeeper",
      &Bmad::attribute_bookkeeper,
      py::arg("ele"),
      py::arg("force_bookkeeping") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element with attributes
    This parameter is an input/output and is modified in-place. As an output: Element with self-consistant
    attributes.
force_bookkeeping : bool, optional
    If present and True then force -- Logical, optional: If present and True then force attribute bookkeeping
    to be done independent of the state of ele.bookkeeping_stat.attributes. This will also cause
    attribute_bookkeeper to assume intelligent bookkeeping.
)""");
  py::class_<Bmad::AttributeFree1, std::unique_ptr<Bmad::AttributeFree1>>(
      m, "AttributeFree1", "attribute_free1 return type")
      .def_readonly("why_not_free", &Bmad::AttributeFree1::why_not_free)
      .def_readonly("free", &Bmad::AttributeFree1::free)
      .def("__len__", [](const Bmad::AttributeFree1&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::AttributeFree1& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.why_not_free);
            if (i == 1)
              return py::cast(s.free);
            throw py::index_error();
          });
  m.def(
      "attribute_free",
      py::overload_cast<
          int,
          std::string,
          LatProxy&,
          std::optional<bool>,
          std::optional<bool>,
          std::optional<bool>>(&Bmad::attribute_free),
      py::arg("ix_ele"),
      py::arg("attrib_name"),
      py::arg("lat"),
      py::arg("err_print_flag") = py::none(),
      py::arg("except_overlay") = py::none(),
      py::arg("dependent_attribs_free") = py::none(),
      R"""(Function attribute_free

Overloaded function for:
  Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag,
                               except_overlay, dependent_attribs_free, why_not_free) result (free)
  Function attribute_free2 (ele, attrib_name, err_print_flag,
                               except_overlay, dependent_attribs_free, why_not_free) result (free)
  Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag,
                               except_overlay, why_not_free) result (free)

Routine to check if an attribute is free to vary.

Attributes that cannot be changed directly include super_slave attributes (since
these attributes are controlled by their super_lords) and attributes that
are controlled by an overlay.

Also dependent variables such as the angle of a bend cannot be
  freely variable.

Parameters
----------
ix_ele : int
    Index of element in element array.
ix_branch : int
    Branch index of element.
ele : EleStruct
    Element containing the attribute
attrib_name : unknown
    Name of the attribute. Assumed upper case.
lat : LatStruct
    Lattice structure.
err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is not free.
except_overlay : bool, optional
    If present and True then an attribute that is controlled by an overlay will be treated as free. This is
    used by, for example, the create_overlay routine.
dependent_attribs_free : bool, optional
    If present and True then mark as free attributes that are dependent. For example, if ele.field_master = F,
    b1_field is dependent upon k1. Default is False. Use True when using intelligent bookkeeping.

Returns
-------
free : bool
    Set True if attribtute not found or attriubte cannot be changed directly.
why_not_free : int
    Possibilities are: field_master_dependent$  -> Dependent due to setting of ele.field_master. dependent$
    -> Not field_master_dependent$ but value is dependent upon the value of other attributes. does_not_exist$
    -> Attribute name is unrecognized or does not exist for the type of element. overlay_slave$           ->
    Attribute is controlled by an overlay lord. super_slave$             -> Attribute is controlled by
    element's super_lord. multipass_slave$         -> Attribute is controlled by element's multipass_lord.
)""");
  py::class_<Bmad::AttributeFree2, std::unique_ptr<Bmad::AttributeFree2>>(
      m, "AttributeFree2", "attribute_free2 return type")
      .def_readonly("why_not_free", &Bmad::AttributeFree2::why_not_free)
      .def_readonly("free", &Bmad::AttributeFree2::free)
      .def("__len__", [](const Bmad::AttributeFree2&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::AttributeFree2& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.why_not_free);
            if (i == 1)
              return py::cast(s.free);
            throw py::index_error();
          });
  m.def(
      "attribute_free",
      py::overload_cast<
          EleProxy&,
          std::string,
          std::optional<bool>,
          std::optional<bool>,
          std::optional<bool>>(&Bmad::attribute_free),
      py::arg("ele"),
      py::arg("attrib_name"),
      py::arg("err_print_flag") = py::none(),
      py::arg("except_overlay") = py::none(),
      py::arg("dependent_attribs_free") = py::none(),
      R"""(Function attribute_free

Overloaded function for:
  Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag,
                               except_overlay, dependent_attribs_free, why_not_free) result (free)
  Function attribute_free2 (ele, attrib_name, err_print_flag,
                               except_overlay, dependent_attribs_free, why_not_free) result (free)
  Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag,
                               except_overlay, why_not_free) result (free)

Routine to check if an attribute is free to vary.

Attributes that cannot be changed directly include super_slave attributes (since
these attributes are controlled by their super_lords) and attributes that
are controlled by an overlay.

Also dependent variables such as the angle of a bend cannot be
  freely variable.

Parameters
----------
ix_ele : int
    Index of element in element array.
ix_branch : int
    Branch index of element.
ele : EleStruct
    Element containing the attribute
attrib_name : unknown
    Name of the attribute. Assumed upper case.
lat : LatStruct
    Lattice structure.
err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is not free.
except_overlay : bool, optional
    If present and True then an attribute that is controlled by an overlay will be treated as free. This is
    used by, for example, the create_overlay routine.
dependent_attribs_free : bool, optional
    If present and True then mark as free attributes that are dependent. For example, if ele.field_master = F,
    b1_field is dependent upon k1. Default is False. Use True when using intelligent bookkeeping.

Returns
-------
free : bool
    Set True if attribtute not found or attriubte cannot be changed directly.
why_not_free : int
    Possibilities are: field_master_dependent$  -> Dependent due to setting of ele.field_master. dependent$
    -> Not field_master_dependent$ but value is dependent upon the value of other attributes. does_not_exist$
    -> Attribute name is unrecognized or does not exist for the type of element. overlay_slave$           ->
    Attribute is controlled by an overlay lord. super_slave$             -> Attribute is controlled by
    element's super_lord. multipass_slave$         -> Attribute is controlled by element's multipass_lord.
)""");
  py::class_<Bmad::AttributeFree3, std::unique_ptr<Bmad::AttributeFree3>>(
      m, "AttributeFree3", "attribute_free3 return type")
      .def_readonly("why_not_free", &Bmad::AttributeFree3::why_not_free)
      .def_readonly("free", &Bmad::AttributeFree3::free)
      .def("__len__", [](const Bmad::AttributeFree3&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::AttributeFree3& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.why_not_free);
            if (i == 1)
              return py::cast(s.free);
            throw py::index_error();
          });
  m.def(
      "attribute_free",
      py::overload_cast<
          int,
          int,
          std::string,
          LatProxy&,
          std::optional<bool>,
          std::optional<bool>,
          std::optional<bool>>(&Bmad::attribute_free),
      py::arg("ix_ele"),
      py::arg("ix_branch"),
      py::arg("attrib_name"),
      py::arg("lat"),
      py::arg("err_print_flag") = py::none(),
      py::arg("except_overlay") = py::none(),
      py::arg("dependent_attribs_free") = py::none(),
      R"""(Function attribute_free

Overloaded function for:
  Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag,
                               except_overlay, dependent_attribs_free, why_not_free) result (free)
  Function attribute_free2 (ele, attrib_name, err_print_flag,
                               except_overlay, dependent_attribs_free, why_not_free) result (free)
  Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag,
                               except_overlay, why_not_free) result (free)

Routine to check if an attribute is free to vary.

Attributes that cannot be changed directly include super_slave attributes (since
these attributes are controlled by their super_lords) and attributes that
are controlled by an overlay.

Also dependent variables such as the angle of a bend cannot be
  freely variable.

Parameters
----------
ix_ele : int
    Index of element in element array.
ix_branch : int
    Branch index of element.
ele : EleStruct
    Element containing the attribute
attrib_name : unknown
    Name of the attribute. Assumed upper case.
lat : LatStruct
    Lattice structure.
err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is not free.
except_overlay : bool, optional
    If present and True then an attribute that is controlled by an overlay will be treated as free. This is
    used by, for example, the create_overlay routine.
dependent_attribs_free : bool, optional
    If present and True then mark as free attributes that are dependent. For example, if ele.field_master = F,
    b1_field is dependent upon k1. Default is False. Use True when using intelligent bookkeeping.

Returns
-------
free : bool
    Set True if attribtute not found or attriubte cannot be changed directly.
why_not_free : int
    Possibilities are: field_master_dependent$  -> Dependent due to setting of ele.field_master. dependent$
    -> Not field_master_dependent$ but value is dependent upon the value of other attributes. does_not_exist$
    -> Attribute name is unrecognized or does not exist for the type of element. overlay_slave$           ->
    Attribute is controlled by an overlay lord. super_slave$             -> Attribute is controlled by
    element's super_lord. multipass_slave$         -> Attribute is controlled by element's multipass_lord.
)""");
  py::class_<Bmad::AttributeIndex1, std::unique_ptr<Bmad::AttributeIndex1>>(
      m, "AttributeIndex1", "attribute_index1 return type")
      .def_readonly("full_name", &Bmad::AttributeIndex1::full_name)
      .def_readonly("attrib_index", &Bmad::AttributeIndex1::attrib_index)
      .def("__len__", [](const Bmad::AttributeIndex1&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::AttributeIndex1& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.full_name);
            if (i == 1)
              return py::cast(s.attrib_index);
            throw py::index_error();
          });
  m.def(
      "attribute_index",
      py::overload_cast<
          EleProxy&,
          std::string,
          std::optional<bool>,
          std::optional<bool>>(&Bmad::attribute_index),
      py::arg("ele"),
      py::arg("name"),
      py::arg("can_abbreviate") = py::none(),
      py::arg("print_error") = py::none(),
      R"""(Function attribute_index (...) result (attrib_index)

Function to return the index of a attribute for a given BMAD element type
and the name of the attribute. Abbreviations are by default permitted but must be at
least 3 characters. Exception: overlay and group varialbe names may not
be abbreviated.

This routine is an overloaded name for:
  attribute_index1 (ele, name, full_name, can_abbreviate, print_error) result (attrib_index)
  attribute_index2 (key, name, full_name, can_abbreviate, print_error) result (attrib_index)

Note:
  If ele%key or key = 0 -> Entire name table will be searched.

See also:
  has_attribute
  attribute_info
  attribute_name

Parameters
----------
ele : EleStruct
    attribute_index will restrict the name search to valid attributes of the given element.
key : int
    Equivalent to ele.key.
name : unknown
    Attribute name. Must be uppercase.
can_abbreviate : bool, optional
    Can abbreviate names? Default is True.
print_error : bool, optional
    Default True. If false, do not print error message.

Returns
-------
full_name : unknown
    Non-abbreviated name.
attrib_index : int
    Index of the attribute. If the attribute name is not appropriate then 0 will be returned. Example: ele.key
    = sbend$ ix = attribute_index (ele, 'K1') Result: ix -> k1$

Notes
-----
Overloaded versions:
)""");
  py::class_<Bmad::AttributeIndex2, std::unique_ptr<Bmad::AttributeIndex2>>(
      m, "AttributeIndex2", "attribute_index2 return type")
      .def_readonly("full_name", &Bmad::AttributeIndex2::full_name)
      .def_readonly("attrib_index", &Bmad::AttributeIndex2::attrib_index)
      .def("__len__", [](const Bmad::AttributeIndex2&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::AttributeIndex2& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.full_name);
            if (i == 1)
              return py::cast(s.attrib_index);
            throw py::index_error();
          });
  m.def(
      "attribute_index",
      py::overload_cast<
          int,
          std::string,
          std::optional<bool>,
          std::optional<bool>>(&Bmad::attribute_index),
      py::arg("key"),
      py::arg("name"),
      py::arg("can_abbreviate") = py::none(),
      py::arg("print_error") = py::none(),
      R"""(Function attribute_index (...) result (attrib_index)

Function to return the index of a attribute for a given BMAD element type
and the name of the attribute. Abbreviations are by default permitted but must be at
least 3 characters. Exception: overlay and group varialbe names may not
be abbreviated.

This routine is an overloaded name for:
  attribute_index1 (ele, name, full_name, can_abbreviate, print_error) result (attrib_index)
  attribute_index2 (key, name, full_name, can_abbreviate, print_error) result (attrib_index)

Note:
  If ele%key or key = 0 -> Entire name table will be searched.

See also:
  has_attribute
  attribute_info
  attribute_name

Parameters
----------
ele : EleStruct
    attribute_index will restrict the name search to valid attributes of the given element.
key : int
    Equivalent to ele.key.
name : unknown
    Attribute name. Must be uppercase.
can_abbreviate : bool, optional
    Can abbreviate names? Default is True.
print_error : bool, optional
    Default True. If false, do not print error message.

Returns
-------
full_name : unknown
    Non-abbreviated name.
attrib_index : int
    Index of the attribute. If the attribute name is not appropriate then 0 will be returned. Example: ele.key
    = sbend$ ix = attribute_index (ele, 'K1') Result: ix -> k1$

Notes
-----
Overloaded versions:
)""");
  m.def(
      "attribute_name",
      py::overload_cast<int, int, std::optional<bool>>(&Bmad::attribute_name),
      py::arg("key"),
      py::arg("ix_att"),
      py::arg("show_private") = py::none(),
      R"""(Function attribute_name (...) result (attrib_name)

Function to return the name of an attribute for a particular type of
Bmad element.

This routine is an overloaded name for:
  attribute_name1 (ele, ix_att, show_private) result (attrib_name)
  attribute_name2 (key, ix_att, show_private) result (attrib_name)


Note: attribute_name (key, ix_att) is not able to handle overlay/group control variables.
Use attributge_name (ele, ix_att) is this is needed.

Parameters
----------
ele : EleStruct
    .key             -- Integer: Key name of element type (e.g. SBEND$, etc.)
key : int
    Key name of element type (e.g. sbend$, etc.)
ix_att : int
    Index of attribute (e.g. k1$)
show_private : bool, optional
    If False (default) return null_name$ for private attributes.

Returns
-------
attrib_name : unknown
    Name of attribute. First character is a "!" if there is a problem. Will always be upper case (even with
    private attributes). = "!BAD ELE KEY"           .key is invalid = "!BAD INDEX"             ix_att is
    invalid (out of range). = "!NULL" (null_name$)     ix_att does not correspond to an attribute or is
    private. Example: ele.key = sbend$ name = attribute_name (ele, k1$) Result: name -> "K1"

Notes
-----
Overloaded versions:
)""");
  m.def(
      "attribute_name",
      py::overload_cast<EleProxy&, int, std::optional<bool>>(
          &Bmad::attribute_name),
      py::arg("ele"),
      py::arg("ix_att"),
      py::arg("show_private") = py::none(),
      R"""(Function attribute_name (...) result (attrib_name)

Function to return the name of an attribute for a particular type of
Bmad element.

This routine is an overloaded name for:
  attribute_name1 (ele, ix_att, show_private) result (attrib_name)
  attribute_name2 (key, ix_att, show_private) result (attrib_name)


Note: attribute_name (key, ix_att) is not able to handle overlay/group control variables.
Use attributge_name (ele, ix_att) is this is needed.

Parameters
----------
ele : EleStruct
    .key             -- Integer: Key name of element type (e.g. SBEND$, etc.)
key : int
    Key name of element type (e.g. sbend$, etc.)
ix_att : int
    Index of attribute (e.g. k1$)
show_private : bool, optional
    If False (default) return null_name$ for private attributes.

Returns
-------
attrib_name : unknown
    Name of attribute. First character is a "!" if there is a problem. Will always be upper case (even with
    private attributes). = "!BAD ELE KEY"           .key is invalid = "!BAD INDEX"             ix_att is
    invalid (out of range). = "!NULL" (null_name$)     ix_att does not correspond to an attribute or is
    private. Example: ele.key = sbend$ name = attribute_name (ele, k1$) Result: name -> "K1"

Notes
-----
Overloaded versions:
)""");
  m.def(
      "attribute_type",
      &Bmad::attribute_type,
      py::arg("attrib_name"),
      py::arg("ele") = py::none(),
      R"""(Function attribute_type (attrib_name, ele) result (attrib_type)

Routine to return the logical type of an attribute.

A "switch" attribute is an attribute whose value corresponds to some string.
For example, the "COUPLER_AT" attirbute with value 1 corresponds to "ENTRANCE_END", etc.

A "struct" attribute is an attribute that is the name for a "structure". For example,
CARTESIAN_MAP is the name of the structure hoding a Cartesian map.

If attrib_name corresponds to a switch attribute, The routine switch_attrib_value_name can
be used to print the name corresponding to the attribute's value.

Note: The "storage type" of an attribute is different from the "logical type" returned by
this routine. For example, the logical type of attribute "n_slice" is integer. However, the
value of "n_slice" is stored as a real number in the ele_struct [in ele%value(n_slice$)].

Parameters
----------
attrib_name : unknown
    Name of the attribute. Must be upper case.
ele : EleStruct, optional
    Element associated with the attribute. Needed if attrib_name can correspond to an overlay or group
    variable.

Returns
-------
attrib_type : int
    Attribute type: is_string$, is_logical$, is_integer$, is_real$, is_switch$, is_struct$ or invalid_name$
    Note: An overlay or group variable will be marked invalid_name$ if ele is missing.
)""");
  m.def(
      "attribute_units",
      &Bmad::attribute_units,
      py::arg("attrib_name"),
      py::arg("unrecognized_units") = py::none(),
      R"""(Function attribute_units (attrib_name, unrecognized_units) result (attrib_units)

Routine to return the units associated with an attribute.
Example: attrib_units('P0C') -> 'eV'

Parameters
----------
attrib_name : unknown
    Name of the attribute. Must be upper case.
unrecognized_units : unknown, optional
    String to use if the attribute name is unrecognized. Note: Non-real attributes (EG: 'TRACKING_METHOD') are
    not recognized. Default is ""

Returns
-------
attrib_units : unknown
    Units associated with the attribute.
)""");
  m.def(
      "autoscale_phase_and_amp",
      &Bmad::autoscale_phase_and_amp,
      py::arg("ele"),
      py::arg("param"),
      py::arg("scale_phase") = py::none(),
      py::arg("scale_amp") = py::none(),
      py::arg("call_bookkeeper") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    RF element or e_gun.
    This parameter is an input/output and is modified in-place. As an output: element with phase and amplitude
    adjusted.
param : LatParamStruct
    lattice parameters
err_flag : 
    Logical, Set true if there is an error. False otherwise.
scale_phase : bool, optional
    Scale the phase? See above.
scale_amp : bool, optional
    Scale the amplitude? See above.
call_bookkeeper : bool, optional
    Call lattice_bookkeeper at end? Default is True.
)""");
  m.def(
      "average_twiss",
      &Bmad::average_twiss,
      py::arg("frac1"),
      py::arg("twiss1"),
      py::arg("twiss2"),
      py::arg("ave_twiss"),
      R"""(Parameters
----------
frac1 : float
    Fraction of twiss1 to use in the average.
twiss1 : TwissStruct
    Twiss parameters to average.
twiss2 : 
ave_twiss : 
)""");
}
