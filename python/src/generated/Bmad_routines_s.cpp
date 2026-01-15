#include "pybmad/generated/Bmad_routines_s.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_s(py::module &m) {
  m.def(
      "s_body_calc",
      &Bmad::s_body_calc,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("s_body"),
      R"""(Parameters
----------
orbit : CoordStruct
    Particle coordinates.
ele : EleStruct
    Lattice element
s_body : 
)"""
  );
  m.def(
      "s_calc",
      &Bmad::s_calc,
      py::arg("lat"),
      R"""(Parameters
----------
lat : LatStruct
)"""
  );
  m.def(
      "sad_mult_hard_bend_edge_kick",
      &Bmad::sad_mult_hard_bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine sad_mult_hard_bend_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

Routine to track through the hard edge bend fringe field for a bend or sad_mult element.
Only the bend field is taken into account here. Higher order multipolse must be handled elsewhere.

This routine assumes that the particle coordinates are with respect to the actual magnet face.
Thus finite e1/e2 must be taken into account by other routines.

SAD calls this the "linear" fringe even though it is nonlinear.

Parameters
----------
ele : EleStruct
    Element with fringe.
param : LatParamStruct
    Tracking parameters.
particle_at : int
    Either first_track_edge$ or second_track_edge$.
orbit : CoordStruct
    Starting coordinates.
    This parameter is an input/output and is modified in-place. As an output: Ending coordinates.
mat6 : float, optional
    Transfer matrix up to the fringe.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix including the
    fringe.
make_matrix : float, optional
    Propagate the transfer matrix? Default is False.
)"""
  );
  m.def(
      "sad_soft_bend_edge_kick",
      &Bmad::sad_soft_bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orb"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine sad_soft_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

Subroutine to track through the ("linear") bend soft edge field of an sbend or sad_mult.

Parameters
----------
ele : EleStruct
    SBend or sad_mult element.
param : LatParamStruct
particle_at : int
    first_track_edge$, or second_track_edge$.
orb : CoordStruct
    Starting coords.
    This parameter is an input/output and is modified in-place. As an output: Coords after tracking.
mat6 : float, optional
    Starting matrix
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix after fringe
    field
make_matrix : float, optional
    Propagate the transfer matrix? Default is False.
k0l : float, optional
    Used with sad_mult. If present, use this instead of ele.a_pole/.b_pole.
t0 : float, optional
    Used with sad_mult. If present, use this instead of ele.a_pole/.b_pole. Must be present if k0l is.
)"""
  );
  m.def(
      "save_a_beam_step",
      &Bmad::save_a_beam_step,
      py::arg("ele"),
      py::arg("beam"),
      py::arg("bunch_tracks") = py::none(),
      py::arg("s_body") = py::none(),
      py::arg("is_time_coords") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element being tracked through.
beam : BeamStruct
    Bunches in the beam whose parameters are to be saved.
bunch_tracks : BunchTrackStruct, optional
    Track with current bunch info appended on. This routine does nothing if this argument is not present.
s_body : float, optional
    Body s-position from beginning of element.
is_time_coords : bool, optional
    Default is False. If True, input beam is using time coordinates in which case there will be a conversion
    to s-coords before bunch_params are computed. Ouput:
)"""
  );
  m.def(
      "save_a_bunch_step",
      &Bmad::save_a_bunch_step,
      py::arg("ele"),
      py::arg("bunch"),
      py::arg("bunch_track") = py::none(),
      py::arg("s_body") = py::none(),
      py::arg("is_time_coords") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element being tracked through.
bunch : BunchStruct
    Bunch whose parameters are to be saved.
bunch_track : BunchTrackStruct, optional
    Track with current bunch info appended on. This routine does nothing if this argument is not present.
s_body : float, optional
    Body s-position from beginning of element.
is_time_coords : bool, optional
    Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
    to s-coords before bunch_params are computed. Ouput:
)"""
  );
  m.def(
      "save_a_step",
      &Bmad::save_a_step,
      py::arg("track"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("local_ref_frame"),
      py::arg("orb"),
      py::arg("s_rel"),
      py::arg("save_field") = py::none(),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::arg("rf_time") = py::none(),
      py::arg("strong_beam") = py::none(),
      R"""(Parameters
----------
track : TrackStruct
    Track with current position appended on.
ele : EleStruct
    Element being tracked through.
param : LatParamStruct
    Lattice parameters.
local_ref_frame : bool
    If True then input orb is with respect to body coordinates.
orb : CoordStruct
    trajectory at s with respect to element coordinates.
s_rel : float
    Longitudinal position wrt the element. If local_ref_frame = F: Lab coords. If local_ref_frame = T: body
    coords.
save_field : bool, optional
    Save electric and magnetic field values? Default is False.
mat6 : float, optional
    Matrix to store.
make_matrix : bool, optional
    Is mat6 a valid matrix? Default is False.
rf_time : float, optional
    RF clock time used for calculating the field.. If not present then the time will be calculated using the
    standard algorithm. This is only needed if save_field = True.
strong_beam : StrongBeambeamStruct, optional
    Strong beam info if tracking through a beambeam element. Ouput:
)"""
  );
  m.def(
      "sbend_body_with_k1_map",
      &Bmad::sbend_body_with_k1_map,
      py::arg("ele"),
      py::arg("dg"),
      py::arg("b1"),
      py::arg("param"),
      py::arg("n_step"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Sbend element.
dg : float
    Field error.
b1 : float
    b1 quadrupole strength * rel_charge_dir
param : LatParamStruct
    Branch parameters.
n_step : int
    Number of steps to divide the bend into. Only one step is taken by this routine.
orbit : CoordStruct
    Orbit at beginning of the bend.
    This parameter is an input/output and is modified in-place. As an output: Ending coordinates.
mat6 : float, optional
    Transfer matrix before element.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix with body added
    in.
make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "sc_adaptive_step",
      &Bmad::sc_adaptive_step,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("include_image"),
      py::arg("t_now"),
      py::arg("dt_step"),
      py::arg("sc_field"),
      R"""(Subroutine sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, dt_next)

Routine to track a bunch of particles with space charge for one step using
adaptive step size control and determine appropriate step size for the next step

Parameters
----------
bunch : BunchStruct
    Starting bunch position in t-based coordinates
    This parameter is an input/output and is modified in-place. As an output: Ending bunch position in t-based
    coordinates.
ele : EleStruct
    Nominal lattice element being tracked through.
include_image : bool
    Include image charge forces?
    This parameter is an input/output and is modified in-place. As an output: Set False if image charge calc
    no longer needed (Note
t_now : float
    Current time at the beginning of tracking
dt_step : float
    Initial SC time step to take
    This parameter is an input/output and is modified in-place. As an output: Step done.
sc_field : unknown
    : Array to hold space charge fields. Its length should be the number of particles.

Returns
-------
dt_next : float
    Next SC time step the tracker would take based on the error tolerance
)"""
  );
  m.def(
      "sc_step",
      &Bmad::sc_step,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("include_image"),
      py::arg("t_end"),
      py::arg("sc_field"),
      R"""(Subroutine sc_step(bunch, ele, include_image, t_end, n_emit)

Subroutine to track a bunch through a given time step with space charge

Parameters
----------
bunch : BunchStruct
    Starting bunch position in t-based coordinates
    This parameter is an input/output and is modified in-place. As an output: Ending bunch position in t-based
    coordinates after space charge kick.
ele : EleStruct
    Nominal element being tracked through.
include_image : bool
    Include image charge forces?
    This parameter is an input/output and is modified in-place. As an output: Set False if image charge calc
    no longer needed (Note
t_end : float
    Time at which the tracking ends.
sc_field : unknown
    : Array to hold space charge fields. Its length should be the number of particles.

Returns
-------
n_emit : int
    The number of particles emitted in this step.
)"""
  );
  m.def(
      "set_active_fixer",
      &Bmad::set_active_fixer,
      py::arg("fixer"),
      py::arg("turn_on") = py::none(),
      R"""(Subroutine set_active_fixer(fixer, turn_on, orbit)

Set the acvitive fixer element.
All other fixer/beginning_ele elements in the branch will be deactivated.

If turn_on is True (default), the fixer argument becomes the active fixer.
If turn_on is False, and fixer%is_on is also False, there is nothing to be done.
If turn_on is False, and fixer%is_on is True, turn this fixer off and turn on the beginning element.

Parameters
----------
fixer : EleStruct
    Fixer element to make active.
    This parameter is an input/output and is modified in-place. As an output: Element is now active.
turn_on : bool, optional
    If True (default), make this fixer the active element. If False, make the beginning element active.

Returns
-------
orbit : CoordStruct
    Load with stored fixer phase space and spin values.
)"""
  );
  m.def(
      "set_custom_attribute_name",
      &Bmad::set_custom_attribute_name,
      py::arg("custom_name"),
      py::arg("custom_index") = py::none(),
      R"""(Subroutine set_custom_attribute_name (custom_name, err_flag, custom_index)

Routine to add custom element attributes to the element attribute name table.

Parameters
----------
custom_name : unknown
    Name of the custom attribute. If prefixed by "<class>::" then the custom name will be set only for that
    element class. Example: "quadrupole::error" will set the alias custom namefor quadrupoles.
custom_index : int, optional
    Index used in assigning where in the ele_struct the custom attribute is put. If not present or 0 then the
    next unused slot is used.

Returns
-------
err_flag : bool
    Set True if an error. False otherwise.
)"""
  );
  py::class_<Bmad::SetEleAttribute, std::unique_ptr<Bmad::SetEleAttribute>>(
      m,
      "SetEleAttribute",
      "set_ele_attribute return type"
  )
      .def_readonly("err_flag", &Bmad::SetEleAttribute::err_flag)
      .def_readonly("err_id", &Bmad::SetEleAttribute::err_id)
      .def("__len__", [](const Bmad::SetEleAttribute &) { return 2; })
      .def("__getitem__", [](const Bmad::SetEleAttribute &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.err_id);
        throw py::index_error();
      });
  m.def(
      "set_ele_attribute",
      &Bmad::set_ele_attribute,
      py::arg("ele"),
      py::arg("set_string"),
      py::arg("err_print_flag") = py::none(),
      py::arg("set_lords") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element with attribute to set.
    This parameter is an input/output and is modified in-place. As an output: Element with attribute set.
set_string : unknown
    Attribute and value for set.
err_flag : bool
    Set True if there is an error, False otherwise.
err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is, for example, not free.
set_lords : bool, optional
    Default False. If True, set the super_lord(s) if the element is a super_slave.
err_id : int
    Set to an integer which identifies the error type. 0 = no error. The higher the error the further along
    the error was encountered.
)"""
  );
  m.def(
      "set_ele_defaults",
      &Bmad::set_ele_defaults,
      py::arg("ele"),
      py::arg("do_allocate") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element to init. .key          -- Type of element.
    This parameter is an input/output and is modified in-place. As an output: Initialized element.
do_allocate : bool, optional
    Do default allocation of element components? Default is True.
)"""
  );
  m.def(
      "set_ele_name",
      &Bmad::set_ele_name,
      py::arg("ele"),
      py::arg("name"),
      R"""(Parameters
----------
ele : EleStruct
    Element whose name is to be set.
    This parameter is an input/output and is modified in-place. As an output: Element with name set.
name : unknown
    Name to set.
)"""
  );
  m.def(
      "set_ele_real_attribute",
      &Bmad::set_ele_real_attribute,
      py::arg("ele"),
      py::arg("attrib_name"),
      py::arg("value"),
      py::arg("err_print_flag") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element with attribute to set.
    This parameter is an input/output and is modified in-place. As an output: Element with attribute set.
attrib_name : unknown
    Attribute name.
value : float
    value to set to.
err_flag : bool
    Set True if there is an error, False otherwise.
err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is, for example, not free.
)"""
  );
  py::class_<Bmad::SetEleStatusStale, std::unique_ptr<Bmad::SetEleStatusStale>>(
      m,
      "SetEleStatusStale",
      "set_ele_status_stale return type"
  )
      .def_readonly("ele", &Bmad::SetEleStatusStale::ele)
      .def_readonly("status_group", &Bmad::SetEleStatusStale::status_group)
      .def_readonly("set_slaves", &Bmad::SetEleStatusStale::set_slaves)
      .def("__len__", [](const Bmad::SetEleStatusStale &) { return 3; })
      .def("__getitem__", [](const Bmad::SetEleStatusStale &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ele);
        if (i == 1)
          return py::cast(s.status_group);
        if (i == 2)
          return py::cast(s.set_slaves);
        throw py::index_error();
      });
  m.def(
      "set_ele_status_stale",
      &Bmad::set_ele_status_stale,
      R"""(Parameters
----------
ele : EleStruct
    Element. .bookkeeping_state   -- Status block to set.
status_group : int
    Which flag groups to set. Possibilities are: attribute_group$, control_group$, floor_position_group$,
    s_position_group$, s_and_floor_position_group$, ref_energy_group$, or mat6_group$, all_groups$
set_slaves : bool
    If present and False then do not set the status for any slaves. Default is True.
)"""
  );
  m.def(
      "set_flags_for_changed_attribute",
      py::overload_cast<EleStruct &, int, std::optional<bool>>(
          &Bmad::set_flags_for_changed_attribute
      ),
      py::arg("ele"),
      py::arg("attrib"),
      py::arg("set_dependent") = py::none(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified for use with "intelligent" bookkeeping.
Also will do some dependent variable bookkeeping when a particular attribute has
been altered.

This routine should be called after the attribute has been set.

set_flags_for_changed_attribute is an overloaded name for:
  set_flags_for_changed_lat_attribute (lat, set_dependent)
  set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
  set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
  set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
  set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

The set_flags_for_changed_lat_attribute (lat) routine is used when one
does not know what has changed and wants a complete bookkeeping done.

NOTE: The attribute argument MUST be the component that was changed. For example:
    ele%value(x_offset$) = off_value
    call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
And NOT:
    call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

Parameters
----------
lat : LatStruct
    Lattice being modified.
    This parameter is an input/output and is modified in-place. As an output: Lattice with appropriate
    changes.
ele : 
    ele_struct, Element being modified.
real_attrib : float, optional
    Attribute that has been changed. For example: ele.value(hkick$). If not present then assume everything has
    potentially changed.
int_attrib : int
    Attribute that has been changed. For example: ele.mat6_calc_method.
logic_attrib : unknown
    ele.is_on.
all_attrib : AllPointerStruct
    Pointer to attribute.
set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "set_flags_for_changed_attribute",
      py::overload_cast<LatStruct &, std::optional<bool>>(&Bmad::set_flags_for_changed_attribute),
      py::arg("lat"),
      py::arg("set_dependent") = py::none(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified for use with "intelligent" bookkeeping.
Also will do some dependent variable bookkeeping when a particular attribute has
been altered.

This routine should be called after the attribute has been set.

set_flags_for_changed_attribute is an overloaded name for:
  set_flags_for_changed_lat_attribute (lat, set_dependent)
  set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
  set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
  set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
  set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

The set_flags_for_changed_lat_attribute (lat) routine is used when one
does not know what has changed and wants a complete bookkeeping done.

NOTE: The attribute argument MUST be the component that was changed. For example:
    ele%value(x_offset$) = off_value
    call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
And NOT:
    call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

Parameters
----------
lat : LatStruct
    Lattice being modified.
    This parameter is an input/output and is modified in-place. As an output: Lattice with appropriate
    changes.
ele : 
    ele_struct, Element being modified.
real_attrib : float, optional
    Attribute that has been changed. For example: ele.value(hkick$). If not present then assume everything has
    potentially changed.
int_attrib : int
    Attribute that has been changed. For example: ele.mat6_calc_method.
logic_attrib : unknown
    ele.is_on.
all_attrib : AllPointerStruct
    Pointer to attribute.
set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "set_flags_for_changed_attribute",
      py::overload_cast<EleStruct &, bool, std::optional<bool>>(
          &Bmad::set_flags_for_changed_attribute
      ),
      py::arg("ele"),
      py::arg("attrib"),
      py::arg("set_dependent") = py::none(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified for use with "intelligent" bookkeeping.
Also will do some dependent variable bookkeeping when a particular attribute has
been altered.

This routine should be called after the attribute has been set.

set_flags_for_changed_attribute is an overloaded name for:
  set_flags_for_changed_lat_attribute (lat, set_dependent)
  set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
  set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
  set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
  set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

The set_flags_for_changed_lat_attribute (lat) routine is used when one
does not know what has changed and wants a complete bookkeeping done.

NOTE: The attribute argument MUST be the component that was changed. For example:
    ele%value(x_offset$) = off_value
    call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
And NOT:
    call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

Parameters
----------
lat : LatStruct
    Lattice being modified.
    This parameter is an input/output and is modified in-place. As an output: Lattice with appropriate
    changes.
ele : 
    ele_struct, Element being modified.
real_attrib : float, optional
    Attribute that has been changed. For example: ele.value(hkick$). If not present then assume everything has
    potentially changed.
int_attrib : int
    Attribute that has been changed. For example: ele.mat6_calc_method.
logic_attrib : unknown
    ele.is_on.
all_attrib : AllPointerStruct
    Pointer to attribute.
set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "set_flags_for_changed_attribute",
      py::overload_cast<EleStruct &, std::optional<double>, std::optional<bool>>(
          &Bmad::set_flags_for_changed_attribute
      ),
      py::arg("ele"),
      py::arg("attrib") = py::none(),
      py::arg("set_dependent") = py::none(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified for use with "intelligent" bookkeeping.
Also will do some dependent variable bookkeeping when a particular attribute has
been altered.

This routine should be called after the attribute has been set.

set_flags_for_changed_attribute is an overloaded name for:
  set_flags_for_changed_lat_attribute (lat, set_dependent)
  set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
  set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
  set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
  set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

The set_flags_for_changed_lat_attribute (lat) routine is used when one
does not know what has changed and wants a complete bookkeeping done.

NOTE: The attribute argument MUST be the component that was changed. For example:
    ele%value(x_offset$) = off_value
    call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
And NOT:
    call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

Parameters
----------
lat : LatStruct
    Lattice being modified.
    This parameter is an input/output and is modified in-place. As an output: Lattice with appropriate
    changes.
ele : 
    ele_struct, Element being modified.
real_attrib : float, optional
    Attribute that has been changed. For example: ele.value(hkick$). If not present then assume everything has
    potentially changed.
int_attrib : int
    Attribute that has been changed. For example: ele.mat6_calc_method.
logic_attrib : unknown
    ele.is_on.
all_attrib : AllPointerStruct
    Pointer to attribute.
set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "set_fringe_on_off",
      &Bmad::set_fringe_on_off,
      py::arg("fringe_at"),
      py::arg("ele_end"),
      py::arg("on_or_off"),
      R"""(Parameters
----------
fringe_at : float
    Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
    This parameter is an input/output and is modified in-place. As an output: Modified fringe setting.
ele_end : int
    Element edge: entrance_end$ or exit_end$
on_or_off : int
    Turn on$ or off$
)"""
  );
  m.def(
      "set_lords_status_stale",
      &Bmad::set_lords_status_stale,
      py::arg("ele"),
      py::arg("stat_group"),
      py::arg("control_bookkeeping") = py::none(),
      py::arg("flag") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element
stat_group : int
    which status group to set. floor_position_group$, etc. See set_ele_status_stale for more details.
control_bookkeeping : bool, optional
    Call control_bookkeeper for each lord if needed? -- logical, optional: Call control_bookkeeper for each
    lord if needed? Default if False.
flag : int, optional
    Do not use. For coordinating recursion. ele.lat    -- Lat_struct: Lattice with status flags of lords of
    ele set.
)"""
  );
  m.def(
      "set_on_off",
      &Bmad::set_on_off,
      py::arg("key"),
      py::arg("lat"),
      py::arg("switch_"),
      py::arg("orb") = py::none(),
      py::arg("use_ref_orb") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("saved_values") = py::none(),
      py::arg("attribute") = py::none(),
      py::arg("set_val") = py::none(),
      R"""(Parameters
----------
key : int
    Class name of elements to be turned on or off. [quadrupole$, etc.]
lat : LatStruct
    lattice structure holding the elements.
    This parameter is an input/output and is modified in-place. As an output: Modified lattice.
switch : int
    on$            => Turn elements on. If saved_values argument is present, use this. If not present (only
    for logical attributes), set to True. off$           => Turn elements off (but will not store the present
    state). off_and_save$  => Save on/off state and then turn elements off. save_state$    => Save present
    on/off state. No turning on or off is done. restore_state$ => Restore saved on/off state from saved_values
    argument.
orb : CoordStruct, optional
    Needed for lat_make_mat6
use_ref_orb : bool, optional
    If present and true then use ele.map_ref_orb for the reference orbit for calculating .mat6. Default is
    false.
ix_branch : int, optional
    If present then only set for this lattice branch.
saved_values : float, optional
    Element-by element saved values of the component. Must be present if needed (EG if switch =
    restore_state$, etc.).
    This parameter is an input/output and is modified in-place. As an output: Saved values of the component.
attribute : unknown, optional
    Attribute to turn on/off. Eg: 'K2', 'MULTIPOLE_ON', etc. Default is 'IS_ON'. Must be upper case.
set_val : int, optional
    Value to set to. Overrides normal set value.
)"""
  );
  m.def(
      "set_orbit_to_zero",
      &Bmad::set_orbit_to_zero,
      py::arg("orbit"),
      py::arg("n1"),
      py::arg("n2"),
      py::arg("ix_noset") = py::none(),
      R"""(Parameters
----------
orbit : CoordStruct
    Array with particle positions in the range orbit(n1:n2) set to zero except for orbit(ix_noset).
n1 : int
    Lower bound of orbit(:) array subset.
n2 : int
    Upper bound of orbit(:) array subset.
ix_noset : int, optional
    If present then orbit(ix_noset) will not be zeroed.
)"""
  );
  m.def(
      "set_ptc",
      &Bmad::set_ptc,
      py::arg("e_tot") = py::none(),
      py::arg("particle") = py::none(),
      py::arg("taylor_order") = py::none(),
      py::arg("integ_order") = py::none(),
      py::arg("n_step") = py::none(),
      py::arg("no_cavity") = py::none(),
      py::arg("force_init") = py::none(),
      R"""(Parameters
----------
e_tot : float, optional
    Energy in eV.
particle : int, optional
    Type of particle: electron$, proton$, etc.
taylor_order : int, optional
    Maximum order of the taylor polynomials. 0 => Use default.
integ_order : int, optional
    Default Order for the drift-kick-drift sympletic integrator. Possibilities are: 2, 4, or 6 Default = 2
n_step : int, optional
    Default Number of integration steps. Default = 1
no_cavity : bool, optional
    No RF Cavity exists? Default = False. Corresponds to the nocavity option of the PTC init routine.
    no_cavity = .true. will turn any cavity into a drift.
force_init : bool, optional
    If present and True then force a PTC init.
)"""
  );
  m.def(
      "set_ptc_base_state",
      &Bmad::set_ptc_base_state,
      py::arg("component"),
      py::arg("set_val"),
      R"""(Parameters
----------
component : unknown
    Name of component. "TOTALPATH", "SPIN", "NOCAVITY", "TIME", etc. See the PTC internal_state structure for
    component names.
set_val : bool
    Value to set to. For TOTALPATH, True => 1, False => 0.
old_val : bool
    Old value.
)"""
  );
  m.def(
      "set_ptc_com_pointers",
      &Bmad::set_ptc_com_pointers,
      R"""(Subroutine set_ptc_com_pointers ()

Routine to set ptc_com pointers to PTC global variables.

)"""
  );
  m.def(
      "set_ptc_quiet",
      &Bmad::set_ptc_quiet,
      py::arg("channel"),
      py::arg("set"),
      py::arg("old_val"),
      R"""(Subroutine set_ptc_quiet (channel, set, old_val)

Routine to set the lielib_print(:) array or c_verbose logical to suppress informational messages
that can clutter the output from a program using PTC.

Note: Only suppress printing if ptc_com%print_info_messages = F.

Parameters
----------
channel : int
    Index in the lielib_print(:) array to set. 0 => c_verbose.
set : bool
    If set$ then set lielib_print(:). If unset$ then undo a previous set$.
old_val : int
    Old value needed for set = unset$.
    This parameter is an input/output and is modified in-place. As an output: Saved value for set = set$.
)"""
  );
  m.def(
      "set_ptc_verbose",
      &Bmad::set_ptc_verbose,
      py::arg("on"),
      R"""(Parameters
----------
on : 
)"""
  );
  m.def(
      "set_pwd_ele",
      &Bmad::set_pwd_ele,
      py::arg("lat"),
      py::arg("mode0"),
      py::arg("inductance"),
      R"""(Subroutine set_pwd_ele(lat,mode0,inductance)

Simulates the effect of potential well distortion by adjusting lat%ele(ix_pwd)%taylor(6)%term(2)%coef for an
element in the lattice.  This element will apply a pz kick based on the z coordinate.
Element is assumed to be at lat%ele(1).  The ibs_ring driver program
inserts a taylor element into lat%ele(1) if set to perform pwd calculations.

Parameters
----------
lat : LatStruct
    lattice
mode0 : NormalModesStruct
    .sig_z and .z.sige_e should be populated before calling this subroutine.
inductance : float
    An inductance-like parameter describing the distortion of the potential well.

Returns
-------
sigma : float
    Bunch length
)"""
  );
  m.def(
      "set_status_flags",
      &Bmad::set_status_flags,
      py::arg("stat"),
      R"""(Parameters
----------
bookkeeping_state : BookkeepingStateStruct
stat : int
    bookkeeping status. ok$, stale$, etc.
)"""
  );
  m.def(
      "set_tune",
      &Bmad::set_tune,
      py::arg("phi_a_set"),
      py::arg("phi_b_set"),
      py::arg("dk1"),
      py::arg("eles"),
      py::arg("branch"),
      py::arg("orb"),
      py::arg("print_err") = py::none(),
      py::arg("ok"),
      R"""(Parameters
----------
phi_a_set : float
    Horizontal set tune (radians)
phi_b_set : float
    Vertical set tune (radians)
dk1 : float
    Relative amount to vary a quad in tuning. The variation will be proportional to dk1. Those quads with a
    positive dk1(i) will be varied as one group and the quads with negative dk1(i) will be varied as another
    group. The routine choose_quads_for_set_tune can be used to calculate values for dk1.
eles : ElePointerStruct
    eles(i).ele points to quadrupole corresponding to dk1(i).
branch : BranchStruct
    Lattice branch to tune.
    This parameter is an input/output and is modified in-place. As an output: Q_tuned lattice branch
orb : CoordStruct
    If RF is off: Energy dE/E at which the tune is computed.
    This parameter is an input/output and is modified in-place. As an output: New closed orbit.
print_err : bool, optional
    Print error message if there is a problem? Default is True.
ok : 
)"""
  );
  m.def(
      "set_twiss",
      &Bmad::set_twiss,
      py::arg("branch"),
      py::arg("twiss_ele"),
      py::arg("ix_ele"),
      py::arg("match_deta_ds"),
      py::arg("err_flag"),
      py::arg("print_err") = py::none(),
      R"""(Parameters
----------
branch : BranchStruct
    Branch to modify.
twiss_ele : EleStruct
    Element with desired Twiss parameters.
ix_ele : int
    Match branch.ele(ix_ele) Twiss to twiss_ele.
match_deta_ds : bool
    If True, match deta_ds. If False, match etap.
err_flag : bool
    Set True if there is an error. False otherwise.
print_err : bool, optional
    Print an error message if there is an error? Default is True.
)"""
  );
  m.def(
      "set_z_tune",
      &Bmad::set_z_tune,
      py::arg("branch"),
      py::arg("z_tune"),
      py::arg("print_err") = py::none(),
      R"""(Parameters
----------
branch : BranchStruct
z_tune : float
    Longitudinal tune in radians (must be negative above transition).
ok : bool
    If present, returns true or false if set was successful. If not present, set_z_tune will bomb if tune
    could not be set. Notes: 1) The calculation assumes that Q_z < 1. 2) By convention a positive tune
    signifies a clockwise rotation in phase space so that the transverse tunes are positive. This means the
    longitudinal tune is negative above transition.
print_err : bool, optional
    Default is True. If False, suppress error messages
)"""
  );
  m.def(
      "settable_dep_var_bookkeeping",
      &Bmad::settable_dep_var_bookkeeping,
      py::arg("ele"),
      R"""(Subroutine settable_dep_var_bookkeeping (ele)

Subroutine to initialize dependent variables in an element.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

)"""
  );
  m.def(
      "setup_high_energy_space_charge_calc",
      &Bmad::setup_high_energy_space_charge_calc,
      py::arg("calc_on"),
      py::arg("branch"),
      py::arg("n_part"),
      py::arg("mode"),
      py::arg("beam_init") = py::none(),
      py::arg("closed_orb") = py::none(),
      R"""(Subroutine setup_high_energy_space_charge_calc (calc_on, branch, n_part, mode, beam_init, closed_orb)

Routine to initialize constants needed by the ultra relativistic space charge
tracking routine track1_high_energy_space_charge. This setup routine must be called if
the lattice or any of the other input parameters are changed.

Parameters used:
    a-mode emittance
    b-mode emittance
    sig_z bunch length
    sig_pz relative energy spread

Parameters
----------
calc_on : bool
    Turns on or off the space charge calculation.
branch : BranchStruct
    Lattice for tracking.
n_part : float
    Number of actual particles in a bunch. Used to compute the bunch charge.
mode : NormalModesStruct
    Structure holding the beam info. Will be combined with info in beam_init.
beam_init : BeamInitStruct, optional
    Structure holding beam info. Will be combined with info in mode.
closed_orb : CoordStruct, optional
    Closed orbit. If not present the closed orbit is taken to be zero.
)"""
  );
  m.def(
      "sigma_mat_ptc_to_bmad",
      &Bmad::sigma_mat_ptc_to_bmad,
      py::arg("sigma_mat_ptc"),
      py::arg("beta0"),
      R"""(Subroutine sigma_mat_ptc_to_bmad (sigma_mat_ptc, beta0, sigma_mat_bmad)

Routine to convert a PTC sigma matrix to a Bmad sigma matrix.
The conversion includes the conversion between Bmad and PTC time coordinate systems.

Since PTC uses delta_E/P0c and Bmad uses delta_P/P0c coordinates, and since
the relationship between delta_E and delta_P is nonlinear, this routine
simplifies the calculation and assumes that the particle beta is constant
over the range of particle energies.

Parameters
----------
sigma_mat_ptc : float
    PTC sigma matrix.
beta0 : float
    Reference particle velocity

Returns
-------
sigma_mat_bmad : float
    Bmad sigma matrix.
)"""
  );
  m.def(
      "significant_difference",
      &Bmad::significant_difference,
      py::arg("value1"),
      py::arg("value2"),
      py::arg("abs_tol") = py::none(),
      py::arg("rel_tol") = py::none(),
      py::arg("is_different"),
      R"""(Parameters
----------
value1 : float
    First value.
value2 : float
    Second value.
abs_tol : float, optional
    Absolute tolerance. Default is 0.
rel_tol : float, optional
    Relative tolerance. Default is 0.
is_different : 
)"""
  );
  m.def(
      "skip_ele_blender",
      &Bmad::skip_ele_blender,
      py::arg("ele"),
      py::arg("skip"),
      R"""(Parameters
----------
ele : 
skip : 
)"""
  );
  m.def(
      "slice_lattice",
      &Bmad::slice_lattice,
      py::arg("lat"),
      py::arg("ele_list"),
      py::arg("do_bookkeeping") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Lattice to slice.
    This parameter is an input/output and is modified in-place. As an output: Lattice with unwanted elements
    sliced out.
ele_list : unknown
    List of elements to retain. See the documentation for the lat_ele_locator routine for the syntax of the
    list.
error : bool
    Set True if there is an error Set False if not.
do_bookkeeping : bool, optional
    Default is True. If false, the calling routine is responsible for: * Modifying lat.particle_start if
    needed. * Calculating Twiss functions.
)"""
  );
  m.def(
      "soft_quadrupole_edge_kick",
      &Bmad::soft_quadrupole_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine soft_quadrupole_edge_kick (ele, param, particle_at, orbit, mat6, make_matrix)

Routine to add the SAD "linear" soft edge (for finite f1 or f2).
This routine assumes that the particle orbit has been rotated to the element reference frame.
This routine is called with sad_mult and quadrupole elements.

Parameters
----------
ele : EleStruct
    Element being tracked through
param : LatParamStruct
    Tracking parameters.
particle_at : int
    first_track_edge$, or second_track_edge$.
orbit : CoordStruct
    Position before kick.
    This parameter is an input/output and is modified in-place. As an output: Position after kick.
mat6 : float, optional
    Transfer matrix up to the edge.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix with edge kick
    added on.
make_matrix : float, optional
    Propagate the transfer matrix? Default is False.
)"""
  );
  m.def(
      "sol_quad_mat6_calc",
      &Bmad::sol_quad_mat6_calc,
      py::arg("ks_in"),
      py::arg("k1_in"),
      py::arg("tilt"),
      py::arg("length"),
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Parameters
----------
ks_in : 
k1_in : 
tilt : float
    quadrupole tilt.
length : float
    Sol_quad length.
ele : EleStruct
    Sol_quad element.
orbit : CoordStruct
    Orbit at beginning of the sol_quad.
mat6 : float, optional
    Transfer matrix up to the sol_quad.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix includeing the
    sol_quad.
make_matrix : bool, optional
    Extend the matrix?
)"""
  );
  m.def(
      "solve_psi_adaptive",
      &Bmad::solve_psi_adaptive,
      py::arg("t0"),
      py::arg("t1"),
      py::arg("p0"),
      py::arg("args"),
      R"""(Subroutine solve_psi_adaptive(t0,t1,p0,args,p1)

Solve dpsi/dt for psi(t1) using adaptive steps and method:
  "Implicit Bulirsch-Stoer method of Bader and Deuflhard."

The boundary condition p0 is psi(t0)

Parameters
----------
t0 : float
    initial time
t1 : float
    final time
p0 : float
    Boundary condition psi(t0)
args : float
    Parameters.  See psi_prime comments for details.

Returns
-------
p1 : float
    psi(t1)
)"""
  );
  m.def(
      "solve_psi_fixed_steps",
      &Bmad::solve_psi_fixed_steps,
      py::arg("t0"),
      py::arg("t1"),
      py::arg("p0"),
      py::arg("args"),
      py::arg("t"),
      py::arg("p"),
      R"""(Subroutine solve_psi_fixed_steps(t0,t1,p0,args,t,p)

Solve dpsi/dt for psi(t1) using fixed steps and method:
  "Implicit Bulirsch-Stoer method of Bader and Deuflhard."

The boundary condition p0 is psi(t0).

Number of steps is determined by SIZE(p).

Parameters
----------
t0 : float
    initial time
t1 : float
    final time
p0 : float
    Boundary condition psi(t0)
args : float
    Parameters.  See psi_prime comments for details.

Returns
-------
t : float
    Array of times from t0 to t1
p : float
    Array of psi evaluated at t(:)
)"""
  );
  m.def(
      "sort_complex_taylor_terms",
      &Bmad::sort_complex_taylor_terms,
      py::arg("complex_taylor_in"),
      R"""(subroutine sort_complex_taylor_terms (complex_taylor_in, complex_taylor_sorted)

Subroutine to sort the complex_taylor terms from "lowest" to "highest" of
a complex_taylor series.
This subroutine is needed because what comes out of PTC is not sorted.

Uses function complex_taylor_exponent_index to sort.

Note: complex_taylor_sorted needs to have been initialized.
Note: complex_taylor_sorted cannot be complex_taylor_in. That is it is not legal to write:
          call sort_complex_taylor_terms (this_complex_taylor, this_complex_taylor)

Parameters
----------
complex_taylor_in : ComplexTaylorStruct
    Unsorted complex_taylor series.

Returns
-------
complex_taylor_sorted : ComplexTaylorStruct
    Sorted complex_taylor series.
)"""
  );
  m.def(
      "spin_dn_dpz_from_mat8",
      &Bmad::spin_dn_dpz_from_mat8,
      py::arg("mat_1turn"),
      py::arg("dn_dpz_partial") = py::none(),
      py::arg("dn_dpz"),
      R"""(Parameters
----------
mat_1turn : float
    Spin-orbital matrix.
dn_dpz_partial : float, optional
    dn_dpz_partial(i,:) is dn_dpz with only one osccilation mode "excited". So dn_dpz_partial(1,:) represents
    a-mode excitation, etc.
error : bool
    Set True if there is an error. False otherwise.
dn_dpz : 
)"""
  );
  m.def(
      "spin_dn_dpz_from_qmap",
      &Bmad::spin_dn_dpz_from_qmap,
      py::arg("orb_mat"),
      py::arg("q_map"),
      py::arg("dn_dpz_partial"),
      py::arg("dn_dpz_partial2"),
      py::arg("n0") = py::none(),
      py::arg("dn_dpz"),
      R"""(Parameters
----------
orb_mat : float
    1-turn orbital matrix.
q_map : float
    1-turn spin linear quaternion map.
dn_dpz_partial : float
    ) is dn_dpz with only one osccilation mode "excited". So dn_dpz_partial(1,:) represents a-mode excitation,
    etc.
dn_dpz_partial2 : float
    ) is dn_dpz with only two osccilation modes "excited". So dn_dpz_partial(1,:) represents b-mode and c-mode
    excitation without the a-mode, etc.
error : bool
    Set True if there is an error. False otherwise.
n0 : float
    3,0).
dn_dpz : 
)"""
  );
  m.def(
      "spin_map1_normalize",
      &Bmad::spin_map1_normalize,
      py::arg("spin1"),
      R"""(Parameters
----------
spin1 : float
    Unnormalized spin map.
    This parameter is an input/output and is modified in-place. As an output: Normalized spin map.
)"""
  );
  py::class_<Bmad::SpinMat8ResonanceStrengths, std::unique_ptr<Bmad::SpinMat8ResonanceStrengths>>(
      m,
      "SpinMat8ResonanceStrengths",
      "spin_mat8_resonance_strengths return type"
  )
      .def_readonly("xi_sum", &Bmad::SpinMat8ResonanceStrengths::xi_sum)
      .def_readonly("xi_diff", &Bmad::SpinMat8ResonanceStrengths::xi_diff)
      .def("__len__", [](const Bmad::SpinMat8ResonanceStrengths &) { return 2; })
      .def("__getitem__", [](const Bmad::SpinMat8ResonanceStrengths &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.xi_sum);
        if (i == 1)
          return py::cast(s.xi_diff);
        throw py::index_error();
      });
  m.def(
      "spin_mat8_resonance_strengths",
      &Bmad::spin_mat8_resonance_strengths,
      py::arg("orb_evec"),
      py::arg("mat8"),
      R"""(Parameters
----------
orb_evec : complex
    Orbital eigenvector.
mat8 : float
    Spin/orbital matrix.
xi_sum : float
    Sum resonance strength.
xi_diff : float
    Difference resonance strength.
)"""
  );
  py::class_<Bmad::SpinMatToEigen, std::unique_ptr<Bmad::SpinMatToEigen>>(
      m,
      "SpinMatToEigen",
      "spin_mat_to_eigen return type"
  )
      .def_readonly("orb_eval", &Bmad::SpinMatToEigen::orb_eval)
      .def_readonly("orb_evec", &Bmad::SpinMatToEigen::orb_evec)
      .def_readonly("n0", &Bmad::SpinMatToEigen::n0)
      .def_readonly("spin_evec", &Bmad::SpinMatToEigen::spin_evec)
      .def_readonly("error", &Bmad::SpinMatToEigen::error)
      .def("__len__", [](const Bmad::SpinMatToEigen &) { return 5; })
      .def("__getitem__", [](const Bmad::SpinMatToEigen &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.orb_eval);
        if (i == 1)
          return py::cast(s.orb_evec);
        if (i == 2)
          return py::cast(s.n0);
        if (i == 3)
          return py::cast(s.spin_evec);
        if (i == 4)
          return py::cast(s.error);
        throw py::index_error();
      });
  m.def(
      "spin_mat_to_eigen",
      &Bmad::spin_mat_to_eigen,
      py::arg("orb_mat"),
      py::arg("spin_map"),
      R"""(Parameters
----------
orb_mat : float
    Orbital matrix.
spin_map : float
    Quaternion 0th & 1st order map.
orb_eval : complex
    Eigenvalues.
orb_evec : complex
    Orbital eigenvectors. orb_evec(j,:) is the j^th vector.
n0 : float
    n_0 invariant spin
spin_evec : complex
    Spin eigenvectors. spin_evec(j,:) is the j^th vector.
error : bool
    Set true if there is an error. False otherwise.
)"""
  );
  m.def(
      "spin_omega",
      &Bmad::spin_omega,
      py::arg("field"),
      py::arg("orbit"),
      py::arg("sign_z_vel"),
      py::arg("phase_space_coords") = py::none(),
      py::arg("omega"),
      R"""(Parameters
----------
field : 
orbit : 
sign_z_vel : 
phase_space_coords : 
omega : 
)"""
  );
  py::class_<Bmad::SpinQuatResonanceStrengths, std::unique_ptr<Bmad::SpinQuatResonanceStrengths>>(
      m,
      "SpinQuatResonanceStrengths",
      "spin_quat_resonance_strengths return type"
  )
      .def_readonly("xi_sum", &Bmad::SpinQuatResonanceStrengths::xi_sum)
      .def_readonly("xi_diff", &Bmad::SpinQuatResonanceStrengths::xi_diff)
      .def("__len__", [](const Bmad::SpinQuatResonanceStrengths &) { return 2; })
      .def("__getitem__", [](const Bmad::SpinQuatResonanceStrengths &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.xi_sum);
        if (i == 1)
          return py::cast(s.xi_diff);
        throw py::index_error();
      });
  m.def(
      "spin_quat_resonance_strengths",
      &Bmad::spin_quat_resonance_strengths,
      py::arg("orb_evec"),
      py::arg("spin_q"),
      R"""(Parameters
----------
orb_evec : complex
    Orbital eigenvector.
spin_q : float
    First order spin map.
xi_sum : float
    Sum resonance strength.
xi_diff : float
    Difference resonance strength.
)"""
  );
  m.def(
      "spin_taylor_to_linear",
      &Bmad::spin_taylor_to_linear,
      py::arg("spin_taylor"),
      py::arg("normalize"),
      py::arg("dref_orb"),
      py::arg("is_on"),
      py::arg("spin_map1"),
      R"""(Parameters
----------
spin_taylor : TaylorStruct
    Taylor spin map.
normalize : bool
    If True, normalize the linear map.
dref_orb : float
    Change in Reference orbit: output_map1_ref - input_taylor_ref.
is_on : bool
    Is map turned on? If not spin_map1 will be the unit map.
spin_map1 : 
)"""
  );
  m.def(
      "spinor_to_polar",
      &Bmad::spinor_to_polar,
      py::arg("spinor"),
      py::arg("polar"),
      R"""(Parameters
----------
spinor : complex
    Spinor
polar : 
)"""
  );
  m.def(
      "spinor_to_vec",
      &Bmad::spinor_to_vec,
      py::arg("spinor"),
      py::arg("vec"),
      R"""(Parameters
----------
spinor : complex
    Spinor Output
vec : 
)"""
  );
  m.def(
      "spline_fit_orbit",
      &Bmad::spline_fit_orbit,
      py::arg("start_orb"),
      py::arg("end_orb"),
      py::arg("spline_x"),
      py::arg("spline_y"),
      R"""(Parameters
----------
start_orb : CoordStruct
    Starting coords.
end_orb : CoordStruct
    Ending coords.
spline_x : float
    Spline coefs for the horizontal trajectory.
spline_y : float
    Spline coefs for vertical trajectory.
)"""
  );
  py::class_<Bmad::SplitLat, std::unique_ptr<Bmad::SplitLat>>(
      m,
      "SplitLat",
      "split_lat return type"
  )
      .def_readonly("ix_split", &Bmad::SplitLat::ix_split)
      .def_readonly("split_done", &Bmad::SplitLat::split_done)
      .def_readonly("err_flag", &Bmad::SplitLat::err_flag)
      .def("__len__", [](const Bmad::SplitLat &) { return 3; })
      .def("__getitem__", [](const Bmad::SplitLat &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ix_split);
        if (i == 1)
          return py::cast(s.split_done);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "split_lat",
      &Bmad::split_lat,
      py::arg("lat"),
      py::arg("s_split"),
      py::arg("ix_branch"),
      py::arg("add_suffix") = py::none(),
      py::arg("check_sanity") = py::none(),
      py::arg("save_null_drift") = py::none(),
      py::arg("choose_max") = py::none(),
      py::arg("ix_insert") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Original lat structure.
    This parameter is an input/output and is modified in-place. As an output: Modified lat structure.
s_split : float
    Position at which lat.branch(ix_branch) is to be split.
ix_branch : int
    Index of lat.branch(:) to use.
ix_split : int
    Index of element just before the s = s_split point.
split_done : bool
    True if lat was split.
add_suffix : bool, optional
    If True (default) add '#1' and '#2" suffixes to the split elements.
check_sanity : bool, optional
    If True (default) then call lat_sanity_check after the split to make sure everything is ok.
save_null_drift : bool, optional
    Save a copy of a drift to be split as a null_ele? This is useful when superpositions are done. See
    add_superimpose for more info. Default is False.
err_flag : bool
    Set true if there is an error, false otherwise.
choose_max : bool, optional
    If no splitting of an element is needed, that is, s_split is at an element boundary, there can be multiple
    possible values for ix_split if there exist zero length elements at the split point. If choose_max = True,
    ix_split will be chosen to be the maximum possible index and if choose_max = False ix_split will be chosen
    to be the minimal possible index. If s_split is not at an element boundary, the setting of choose_max is
    immaterial. If ix_insert is present, the default value of choose_max is set to give the closest element to
    ix_insert. If ix_insert is not present, the default value of choose_max is False.
ix_insert : int, optional
    Element index near the point to be split. ix_insert is useful in the case where there is a patch with a
    negative length which can create an ambiguity as to where to do the split In this case ix_insert will
    remove the ambiguity. Also useful to ensure where to split if there are elements with zero length nearby.
    Ignored if negative.
)"""
  );
  m.def(
      "sprint_spin_taylor_map",
      &Bmad::sprint_spin_taylor_map,
      py::arg("ele"),
      py::arg("start_orbit") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element to form map for.
    This parameter is an input/output and is modified in-place. As an output: Element with map.
start_orbit : float, optional
    Reference orbit for the map. Default is zero orbit.
)"""
  );
  m.def(
      "sr_longitudinal_wake_particle",
      &Bmad::sr_longitudinal_wake_particle,
      py::arg("ele"),
      py::arg("orbit"),
      R"""(Subroutine sr_longitudinal_wake_particle (ele, orbit)

Routine to apply the short-range wake longitudinal component kick to a particle and then add
to the existing longitudinal wake the contribution from the particle.

Parameters
----------
ele : EleStruct
    Element with wakes.
orbit : CoordStruct
    Particle coords.
    This parameter is an input/output and is modified in-place. As an output: coords after the kick.
)"""
  );
  m.def(
      "sr_transverse_wake_particle",
      &Bmad::sr_transverse_wake_particle,
      py::arg("ele"),
      py::arg("orbit"),
      R"""(Subroutine sr_transverse_wake_particle (ele, orbit)

Subroutine to apply the short-range wake transverse component of the kick to a particle and then add
to the existing transverse wake the contribution from the particle.

Parameters
----------
ele : EleStruct
    Element with wakes.
orbit : CoordStruct
    Starting particle coords.
    This parameter is an input/output and is modified in-place. As an output: Ending particle coords.
)"""
  );
  m.def(
      "sr_z_long_wake",
      &Bmad::sr_z_long_wake,
      py::arg("ele"),
      py::arg("bunch"),
      py::arg("z_ave"),
      R"""(Subroutine sr_z_long_wake (ele, bunch, z_ave)

Subroutine to apply the short-range z-wake kick to a particle.

Parameters
----------
ele : EleStruct
    Element with wake.
bunch : BunchStruct
    Bunch before wake applied.
z_ave : float
    Average z-position of all live particles.

Returns
-------
orbit : CoordStruct
    Ending particle coords.
)"""
  );
  m.def(
      "srdt_calc",
      &Bmad::srdt_calc,
      py::arg("lat"),
      py::arg("order"),
      py::arg("n_slices_gen_opt") = py::none(),
      py::arg("n_slices_sxt_opt") = py::none(),
      py::arg("per_ele_out") = py::none(),
      R"""(Subroutine srdt_calc(lat, srdt_sums, order, n_slices_gen_opt, n_slices_sxt_opt)

Calculate summation RDT terms up to order=1 or order=2 while slicing sextupoles
n_slices_sxt_opt times and all other elements n_slices_gen_opt times.

These formulas are documented in "The Sextupole Scheme for the Swiss Light Source (SLS): An Analytic Approach"
by Johan Bengtsson.  SLS Note 9/97.

The 2nd order formulas are documented in "Second-order driving terms due to sextupoles and
chromatic effects of quadrupoles" by Chun-xi Wang.  AOP-TN-2009-020.

Parameters
----------
lat : LatStruct
    lattice with Twiss parameters calculated.
order : int
    1 to calculate only first order terms.  2 to also calculate 2nd order terms.
n_slices_gen_opt : int, optional
    number of times to slice elements other than sextupoles.  Default is 10.
n_slices_sxt_opt : int, optional
    nubmer of times to slice sextupoles.  Default is 20.

Returns
-------
srdt_sums : SummationRdtStruct
    contains complex RDT strengths.
)"""
  );
  m.def(
      "srdt_lsq_solution",
      &Bmad::srdt_lsq_solution,
      py::arg("lat"),
      py::arg("var_indexes"),
      py::arg("ls_soln"),
      py::arg("n_slices_gen_opt") = py::none(),
      py::arg("n_slices_sxt_opt") = py::none(),
      py::arg("chrom_set_x_opt") = py::none(),
      py::arg("chrom_set_y_opt") = py::none(),
      py::arg("weight_in") = py::none(),
      R"""(Subroutine srdt_lsq_solution(lat, var_indexes, ls_soln, n_slices_gen_opt, n_slices_sxt_opt,
                                                    chrom_set_x_opt, chrom_set_y_opt, weight_in)

Given lat, finds K2 moments that set the chromaticity and zeros-out the real
and complex parts of the first order driving terms, that minimizes the sum of the squares
of the K2 moments.  i.e. the weakest sextupole scheme that sets chromaticity
and zeros out the first order terms.

Note:  This subroutine does not, in its present form, work well with knobs, overlays, or in lattices where
       multiple elements have the same name.

This subroutine assumes that Nsext > 18.

Parameters
----------
lat : LatStruct
    lattice with Twiss parameters calculated.
var_indexes : int
    indexes in lat.ele that are K2 variables.  Must be sorted smallest index to largest index.
n_slices_gen_opt : int, optional
    number of times to slice elements other than sextupoles.  Default is 10.
n_slices_sxt_opt : int, optional
    nubmer of times to slice sextupoles.  Default is 20.
chrom_set_x_opt : float, optional
    what to set x chromaticity to.  Default zero.
chrom_set_y_opt : float, optional
    what to set y chromaticity to.  Default zero.
weight_in : float, optional
    moment weights. Terms are: [wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002, wgt_h21000,
    wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200, If present, any terms equal to zero are given default
    values which is 1.0e4 for wgt_chrom_x and wgt_chrom_y and is 1.0 for everything else.

Returns
-------
ls_soln : float
    contains K2 for the indexes in var_indexes
)"""
  );
  m.def(
      "start_branch_at",
      &Bmad::start_branch_at,
      py::arg("lat"),
      py::arg("ele_start"),
      py::arg("move_end_marker"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice to modify.
    This parameter is an input/output and is modified in-place. As an output: Modified lattice.
ele_start : unknown
    Start element. Ele_start will identify the lattice branch to modify.
move_end_marker : bool
    If True then the end marker (if it is present) will be shifted like any other element. False means that
    the end marker will stay at the end.
error : bool
    Set True if there is an error Set False if not.
)"""
  );
  m.def(
      "stream_ele_end",
      &Bmad::stream_ele_end,
      py::arg("physical_end"),
      py::arg("ele_orientation"),
      py::arg("stream_end"),
      R"""(Parameters
----------
physical_end : int
    entrance_end$, exit_end$, surface$, etc.
ele_orientation : int
    Either 1 = Normal or -1 = element reversed.
stream_end : 
)"""
  );
  m.def(
      "string_attrib",
      &Bmad::string_attrib,
      py::arg("attrib_name"),
      py::arg("ele"),
      R"""(Subroutine string_attrib (attrib_name, ele, attrib_value)

Routine to return the value of a string attribute of a lattice element.
This routine is useful when attrib_name is specified by the program user.

For example:
  call string_attrib ('NAME', ele, attrib_value)  ! Will return attrib_value = ele%name

Parameters
----------
attrib_name : unknown
    Name of the type of element attribute.
ele : EleStruct
    Lattice element.

Returns
-------
attrib_value : unknown
    The string associated with the attribute.
)"""
  );
  py::class_<Bmad::StrongBeamSigmaCalc, std::unique_ptr<Bmad::StrongBeamSigmaCalc>>(
      m,
      "StrongBeamSigmaCalc",
      "strong_beam_sigma_calc return type"
  )
      .def_readonly("sigma", &Bmad::StrongBeamSigmaCalc::sigma)
      .def_readonly("bbi_const", &Bmad::StrongBeamSigmaCalc::bbi_const)
      .def_readonly("dsigma_ds", &Bmad::StrongBeamSigmaCalc::dsigma_ds)
      .def("__len__", [](const Bmad::StrongBeamSigmaCalc &) { return 3; })
      .def("__getitem__", [](const Bmad::StrongBeamSigmaCalc &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.sigma);
        if (i == 1)
          return py::cast(s.bbi_const);
        if (i == 2)
          return py::cast(s.dsigma_ds);
        throw py::index_error();
      });
  m.def(
      "strong_beam_sigma_calc",
      &Bmad::strong_beam_sigma_calc,
      py::arg("ele"),
      py::arg("s_pos"),
      R"""(Parameters
----------
ele : EleStruct
    Beambeam element.
s_pos : float
    Longitudinal position in lab coords of slice (used with hourglass effect correction).
sigma : float
    Strong beam x,y sigmas.
bbi_const : float
    BBI kick scale factor.
dsigma_ds : float
    sig_x and sig_y longitudinal derivatives.
)"""
  );
  m.def(
      "strong_beam_strength",
      &Bmad::strong_beam_strength,
      py::arg("ele"),
      py::arg("strength"),
      R"""(Parameters
----------
ele : EleStruct
    Beambeam element.
strength : 
)"""
  );
  m.def(
      "surface_grid_displacement",
      &Bmad::surface_grid_displacement,
      py::arg("ele"),
      py::arg("x"),
      py::arg("y"),
      py::arg("err_flag"),
      py::arg("z"),
      py::arg("dz_dxy") = py::none(),
      py::arg("extend_grid") = py::none(),
      R"""(Subroutine surface_grid_displacement (ele, x, y, err_flag, z, dz_dxy, extend_grid)

Routine to add in the z displacement defined by the grid

Parameters
----------
ele : EleStruct
    Element containing the grid x, y          -- real(rp): Photon coords at surface.
extend_grid : bool, optional
    If (x,y) past grid pretend (x,y) is at grid boundary. Default is False. Output
err_flag : bool
    Set True if there is a problem.
z : float
    surface height at (x, y).
dz_dxy : float, optional
    Surface slope at (x, y).
)"""
  );
  m.def(
      "symp_lie_bmad",
      &Bmad::symp_lie_bmad,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      py::arg("offset_ele") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element with transfer matrix
    This parameter is an input/output and is modified in-place. As an output: Element with transfer matrix.
param : LatParamStruct
    Parameters are needed for some elements.
orbit : CoordStruct
    Coordinates at the beginning of element.
    This parameter is an input/output and is modified in-place. As an output: Coordinates at the end of
    element.
track : TrackStruct
    Structure holding the track information. When tracking through multiple elements, the trajectory in an
    element is appended to the existing trajectory. To reset: Set track.n_pt = -1.
mat6 : float, optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place. As an output: Transfer matrix propagated
    through the element.
make_matrix : bool
    If True then make the 6x6 transfer matrix.
offset_ele : bool, optional
    Offset the element using ele.value(x_offset$), etc. Default is True.
)"""
  );
}
