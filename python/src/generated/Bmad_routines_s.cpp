#include "pybmad/generated/Bmad_routines_s.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyScAdaptiveStep python_sc_adaptive_step(
    BunchStruct &bunch,
    EleStruct &ele,
    bool include_image,
    double t_now,
    double dt_step,
    EmFieldStructArray1D sc_field
) {
  auto _result = Bmad::sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, sc_field);
  auto py_result{PyScAdaptiveStep{_result, include_image, dt_step}};
  return py_result;
}
PyScStep python_sc_step(
    BunchStruct &bunch,
    EleStruct &ele,
    bool include_image,
    double t_end,
    EmFieldStructArray1D sc_field
) {
  auto _result = Bmad::sc_step(bunch, ele, include_image, t_end, sc_field);
  auto py_result{PyScStep{_result, include_image}};
  return py_result;
}
PySetFringeOnOff python_set_fringe_on_off(double fringe_at, int ele_end, int on_or_off) {
  Bmad::set_fringe_on_off(fringe_at, ele_end, on_or_off);
  auto py_result{PySetFringeOnOff{fringe_at}};
  return py_result;
}
PySetPtcQuiet python_set_ptc_quiet(int channel, bool set, int old_val) {
  Bmad::set_ptc_quiet(channel, set, old_val);
  auto py_result{PySetPtcQuiet{old_val}};
  return py_result;
}

void init_Bmad_routines_s(py::module &m) {
  m.def(
      "s_body_calc",
      &Bmad::s_body_calc,
      py::arg("orbit"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine s_body_calc

Parameters
----------
orbit : CoordStruct
    Particle coordinates.

ele : EleStruct
    Lattice element

Returns
-------
s_body : float
    Body postion.
)"""
  );
  m.def(
      "s_calc",
      &Bmad::s_calc,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine s_calc

Parameters
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
      py::call_guard<py::gil_scoped_release>(),
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
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coordinates.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the fringe.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix including the fringe.

make_matrix : bool, optional
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
      py::call_guard<py::gil_scoped_release>(),
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
    This parameter is an input/output and is modified in-place.
    As an output, orb: Coords after tracking.

mat6 : 2D array of float (shape: 6,6), optional
    Starting matrix
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix after fringe field

make_matrix : bool, optional
    Propagate the transfer matrix? Default is False.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine save_a_beam_step

Parameters
----------
ele : EleStruct
    Element being tracked through.

beam : BeamStruct
    Bunches in the beam whose parameters are to be saved.

bunch_tracks : 1D array of BunchTrackStruct, optional
    Track up to now. If bunch_tracks.n_pt < 0, the structure will be reinitialized.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_tracks: Track with current bunch info appended on. This routine does nothing

s_body : float, optional
    Body s-position from beginning of element.

is_time_coords : bool, optional
    Default is False. If True, input beam is using time coordinates in which case there will be a conversion
    to s-coords before bunch_params are computed.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine save_a_bunch_step

Parameters
----------
ele : EleStruct
    Element being tracked through.

bunch : BunchStruct
    Bunch whose parameters are to be saved.

bunch_track : BunchTrackStruct, optional
    Track up to now. If bunch_track.n_pt < 0, the structure will be reinitialized.
    This parameter is an input/output and is modified in-place.
    As an output, bunch_track: Track with current bunch info appended on. This routine does nothing

s_body : float, optional
    Body s-position from beginning of element.

is_time_coords : bool, optional
    Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
    to s-coords before bunch_params are computed.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine save_a_step

Parameters
----------
track : TrackStruct
    Track up to now. If track.n_pt < 0, the structure will be reinitialized.
    This parameter is an input/output and is modified in-place.
    As an output, track: Track with current position appended on.

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

mat6 : 2D array of float (shape: 6,6), optional
    Matrix to store.

make_matrix : bool, optional
    Is mat6 a valid matrix? Default is False.

rf_time : float, optional
    RF clock time used for calculating the field.. If not present then the time will be calculated using the
    standard algorithm. This is only needed if save_field = True.

strong_beam : StrongBeamStruct, optional
    Strong beam info if tracking through a beambeam element.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sbend_body_with_k1_map

Parameters
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
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending coordinates.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix with body added in.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  py::class_<PyScAdaptiveStep, std::unique_ptr<PyScAdaptiveStep>>(
      m,
      "ScAdaptiveStep",
      "sc_adaptive_step return type"
  )
      .def_readonly("dt_next", &PyScAdaptiveStep::dt_next)
      .def_readonly("include_image", &PyScAdaptiveStep::include_image)
      .def_readonly("dt_step", &PyScAdaptiveStep::dt_step)
      .def("__len__", [](const PyScAdaptiveStep &) { return 3; })
      .def("__getitem__", [](const PyScAdaptiveStep &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.dt_next);
        if (i == 1)
          return py::cast(s.include_image);
        if (i == 2)
          return py::cast(s.dt_step);
        throw py::index_error();
      });
  m.def(
      "sc_adaptive_step",
      &python_sc_adaptive_step,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("include_image"),
      py::arg("t_now"),
      py::arg("dt_step"),
      py::arg("sc_field"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, dt_next)

Routine to track a bunch of particles with space charge for one step using
adaptive step size control and determine appropriate step size for the next step

Parameters
----------
bunch : BunchStruct
    Starting bunch position in t-based coordinates
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position in t-based coordinates.

ele : EleStruct
    Nominal lattice element being tracked through.

include_image : bool
    Include image charge forces?
    This parameter is an input/output and is modified in-place.
    As an output, include_image: Set False if image charge calc no longer needed (Note

t_now : float
    Current time at the beginning of tracking

dt_step : float
    Initial SC time step to take
    This parameter is an input/output and is modified in-place.
    As an output, dt_step: Step done.

sc_field : 1D array of EmFieldStruct
    : Array to hold space charge fields. Its length should be the number of particles.

Returns
-------
include_image : bool
    Include image charge forces?
    This parameter is an input/output and is modified in-place.
    As an output, include_image: Set False if image charge calc no longer needed (Note

dt_step : float
    Initial SC time step to take
    This parameter is an input/output and is modified in-place.
    As an output, dt_step: Step done.

dt_next : float
    Next SC time step the tracker would take based on the error tolerance
)"""
  );
  py::class_<PyScStep, std::unique_ptr<PyScStep>>(m, "ScStep", "sc_step return type")
      .def_readonly("n_emit", &PyScStep::n_emit)
      .def_readonly("include_image", &PyScStep::include_image)
      .def("__len__", [](const PyScStep &) { return 2; })
      .def("__getitem__", [](const PyScStep &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.n_emit);
        if (i == 1)
          return py::cast(s.include_image);
        throw py::index_error();
      });
  m.def(
      "sc_step",
      &python_sc_step,
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("include_image"),
      py::arg("t_end"),
      py::arg("sc_field"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine sc_step(bunch, ele, include_image, t_end, n_emit)

Subroutine to track a bunch through a given time step with space charge

Parameters
----------
bunch : BunchStruct
    Starting bunch position in t-based coordinates
    This parameter is an input/output and is modified in-place.
    As an output, bunch: Ending bunch position in t-based coordinates after space charge kick.

ele : EleStruct
    Nominal element being tracked through.

include_image : bool
    Include image charge forces?
    This parameter is an input/output and is modified in-place.
    As an output, include_image: Set False if image charge calc no longer needed (Note

t_end : float
    Time at which the tracking ends.

sc_field : 1D array of EmFieldStruct
    : Array to hold space charge fields. Its length should be the number of particles.

Returns
-------
include_image : bool
    Include image charge forces?
    This parameter is an input/output and is modified in-place.
    As an output, include_image: Set False if image charge calc no longer needed (Note

n_emit : int, optional
    The number of particles emitted in this step.
)"""
  );
  m.def(
      "set_active_fixer",
      &Bmad::set_active_fixer,
      py::arg("fixer"),
      py::arg("turn_on") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
    This parameter is an input/output and is modified in-place.
    As an output, fixer: Element is now active.

turn_on : bool, optional
    If True (default), make this fixer the active element. If False, make the beginning element active.

Returns
-------
orbit : CoordStruct, optional
    Load with stored fixer phase space and spin values.
)"""
  );
  m.def(
      "set_custom_attribute_name",
      &Bmad::set_custom_attribute_name,
      py::arg("custom_name"),
      py::arg("custom_index") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine set_custom_attribute_name (custom_name, err_flag, custom_index)

Routine to add custom element attributes to the element attribute name table.

Parameters
----------
custom_name : str
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ele_attribute

Parameters
----------
ele : EleStruct
    Element with attribute to set.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with attribute set.

set_string : str
    Attribute and value for set.

err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is, for example, not free.

set_lords : bool, optional
    Default False. If True, set the super_lord(s) if the element is a super_slave.

Returns
-------
err_flag : bool
    Set True if there is an error, False otherwise.

err_id : int, optional
    Set to an integer which identifies the error type. 0 = no error. The higher the error the further along
    the error was encountered.
)"""
  );
  m.def(
      "set_ele_defaults",
      &Bmad::set_ele_defaults,
      py::arg("ele"),
      py::arg("do_allocate") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ele_defaults

Parameters
----------
ele : EleStruct
    Element to init.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Initialized element.

do_allocate : bool, optional
    Do default allocation of element components? Default is True.
)"""
  );
  m.def(
      "set_ele_name",
      &Bmad::set_ele_name,
      py::arg("ele"),
      py::arg("name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ele_name

Parameters
----------
ele : EleStruct
    Element whose name is to be set.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with name set.

name : str
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ele_real_attribute

Parameters
----------
ele : EleStruct
    Element with attribute to set.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with attribute set.

attrib_name : str
    Attribute name.

value : float
    value to set to.

err_print_flag : bool, optional
    If present and False then suppress printing of an error message if attribute is, for example, not free.

Returns
-------
err_flag : bool
    Set True if there is an error, False otherwise.
)"""
  );
  m.def(
      "set_ele_status_stale",
      &Bmad::set_ele_status_stale,
      py::arg("ele"),
      py::arg("status_group"),
      py::arg("set_slaves") = py::none(),
      py::arg("old_eles") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ele_status_stale

Parameters
----------
ele : EleStruct
    Element to set.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element.

status_group : int
    Which flag groups to set. Possibilities are: attribute_group$, control_group$, floor_position_group$,
    s_position_group$, s_and_floor_position_group$, ref_energy_group$, or mat6_group$, all_groups$

set_slaves : bool, optional
    If present and False then do not set the status for any slaves. Default is True.

old_eles : 1D array of ElePointerStruct, optional
    List of elements already set. This argument is only used when this routine is called recursively.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified.
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
ele : EleStruct
    ele_struct, Element being modified.

set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.
)"""
  );
  m.def(
      "set_flags_for_changed_attribute",
      py::overload_cast<LatStruct &, std::optional<bool>>(&Bmad::set_flags_for_changed_attribute),
      py::arg("lat"),
      py::arg("set_dependent") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified.
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
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with appropriate changes.

set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified.
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
ele : EleStruct
    ele_struct, Element being modified.

set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine set_flags_for_changed_attribute (...)

Routine to mark an element or lattice as modified.
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
ele : EleStruct
    ele_struct, Element being modified.

set_dependent : bool, optional
    If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
    when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
    doing.
)"""
  );
  py::class_<PySetFringeOnOff, std::unique_ptr<PySetFringeOnOff>>(
      m,
      "SetFringeOnOff",
      "set_fringe_on_off return type"
  )
      .def_readonly("fringe_at", &PySetFringeOnOff::fringe_at)
      .def("__len__", [](const PySetFringeOnOff &) { return 1; })
      .def("__getitem__", [](const PySetFringeOnOff &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.fringe_at);
        throw py::index_error();
      });
  m.def(
      "set_fringe_on_off",
      &python_set_fringe_on_off,
      py::arg("fringe_at"),
      py::arg("ele_end"),
      py::arg("on_or_off"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_fringe_on_off

Parameters
----------
fringe_at : float
    Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
    This parameter is an input/output and is modified in-place.
    As an output, fringe_at: Modified fringe setting.

ele_end : int
    Element edge: entrance_end$ or exit_end$

on_or_off : int
    Turn on$ or off$

Returns
-------
fringe_at : float
    Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
    This parameter is an input/output and is modified in-place.
    As an output, fringe_at: Modified fringe setting.
)"""
  );
  m.def(
      "set_lords_status_stale",
      &Bmad::set_lords_status_stale,
      py::arg("ele"),
      py::arg("stat_group"),
      py::arg("control_bookkeeping") = py::none(),
      py::arg("flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_lords_status_stale

Parameters
----------
ele : EleStruct
    Element

stat_group : int
    which status group to set. floor_position_group$, etc. See set_ele_status_stale for more details.

control_bookkeeping : bool, optional
    Call control_bookkeeper for each lord if needed? Default if False.

flag : int, optional
    Do not use. For coordinating recursion.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_on_off

Parameters
----------
key : int
    Class name of elements to be turned on or off. [quadrupole$, etc.]

lat : LatStruct
    lattice structure holding the elements.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Modified lattice.

orb : 1D array of CoordStruct, optional
    Needed for lat_make_mat6

use_ref_orb : bool, optional
    If present and true then use ele.map_ref_orb for the reference orbit for calculating .mat6. Default is
    false.

ix_branch : int, optional
    If present then only set for this lattice branch.

saved_values : 1D array of float, optional
    Element-by element saved values of the component. Must be present if needed (EG if switch =
    restore_state$, etc.).
    This parameter is an input/output and is modified in-place.
    As an output, saved_values: Saved values of the component.

attribute : str, optional
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_orbit_to_zero

Parameters
----------
orbit : 1D array of CoordStruct
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ptc

Parameters
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ptc_base_state

Parameters
----------
component : str
    Name of component. "TOTALPATH", "SPIN", "NOCAVITY", "TIME", etc. See the PTC internal_state structure for
    component names.

set_val : bool
    Value to set to. For TOTALPATH, True => 1, False => 0.

Returns
-------
old_val : bool, optional
    Old value.
)"""
  );
  m.def(
      "set_ptc_com_pointers",
      &Bmad::set_ptc_com_pointers,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine set_ptc_com_pointers ()

Routine to set ptc_com pointers to PTC global variables.
)"""
  );
  py::class_<PySetPtcQuiet, std::unique_ptr<PySetPtcQuiet>>(
      m,
      "SetPtcQuiet",
      "set_ptc_quiet return type"
  )
      .def_readonly("old_val", &PySetPtcQuiet::old_val)
      .def("__len__", [](const PySetPtcQuiet &) { return 1; })
      .def("__getitem__", [](const PySetPtcQuiet &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.old_val);
        throw py::index_error();
      });
  m.def(
      "set_ptc_quiet",
      &python_set_ptc_quiet,
      py::arg("channel"),
      py::arg("set"),
      py::arg("old_val"),
      py::call_guard<py::gil_scoped_release>(),
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
    This parameter is an input/output and is modified in-place.
    As an output, old_val: Saved value for set = set$.

Returns
-------
old_val : int
    Old value needed for set = unset$.
    This parameter is an input/output and is modified in-place.
    As an output, old_val: Saved value for set = set$.
)"""
  );
  m.def(
      "set_ptc_verbose",
      &Bmad::set_ptc_verbose,
      py::arg("on"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_ptc_verbose

Parameters
----------
on : bool
)"""
  );
  m.def(
      "set_pwd_ele",
      &Bmad::set_pwd_ele,
      py::arg("lat"),
      py::arg("mode0"),
      py::arg("inductance"),
      py::call_guard<py::gil_scoped_release>(),
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
)"""
  );
  m.def(
      "set_status_flags",
      &Bmad::set_status_flags,
      py::arg("stat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_status_flags

Parameters
----------
stat : int
    bookkeeping status. ok$, stale$, etc.

Returns
-------
bookkeeping_state : BookkeepingStateStruct
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_tune

Parameters
----------
phi_a_set : float
    Horizontal set tune (radians)

phi_b_set : float
    Vertical set tune (radians)

dk1 : 1D array of float
    Relative amount to vary a quad in tuning. The variation will be proportional to dk1. Those quads with a
    positive dk1(i) will be varied as one group and the quads with negative dk1(i) will be varied as another
    group. The routine choose_quads_for_set_tune can be used to calculate values for dk1.

eles : 1D array of ElePointerStruct
    eles(i).ele points to quadrupole corresponding to dk1(i).

branch : BranchStruct
    Lattice branch to tune.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Q_tuned lattice branch

orb : 1D array of CoordStruct
    If RF is off: Energy dE/E at which the tune is computed.
    This parameter is an input/output and is modified in-place.
    As an output, orb: New closed orbit.

print_err : bool, optional
    Print error message if there is a problem? Default is True.

Returns
-------
ok : bool
    Set True if everything is ok. False otherwise.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_twiss

Parameters
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_z_tune

Parameters
----------
branch : BranchStruct

z_tune : float
    Longitudinal tune in radians (must be negative above transition).

print_err : bool, optional
    Default is True. If False, suppress error messages

Returns
-------
ok : bool, optional
    If present, returns true or false if set was successful. If not present, set_z_tune will bomb if tune
    could not be set.
)"""
  );
  m.def(
      "settable_dep_var_bookkeeping",
      &Bmad::settable_dep_var_bookkeeping,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
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
      py::call_guard<py::gil_scoped_release>(),
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

closed_orb : 1D array of CoordStruct, optional
    Closed orbit. If not present the closed orbit is taken to be zero.
)"""
  );
  m.def(
      "sigma_mat_ptc_to_bmad",
      &Bmad::sigma_mat_ptc_to_bmad,
      py::arg("sigma_mat_ptc"),
      py::arg("beta0"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine sigma_mat_ptc_to_bmad (sigma_mat_ptc, beta0, sigma_mat_bmad)

Routine to convert a PTC sigma matrix to a Bmad sigma matrix.
The conversion includes the conversion between Bmad and PTC time coordinate systems.

Since PTC uses delta_E/P0c and Bmad uses delta_P/P0c coordinates, and since
the relationship between delta_E and delta_P is nonlinear, this routine
simplifies the calculation and assumes that the particle beta is constant
over the range of particle energies.

Parameters
----------
sigma_mat_ptc : 2D array of float (shape: 6,6)
    PTC sigma matrix.

beta0 : float
    Reference particle velocity

Returns
-------
sigma_mat_bmad : 2D array of float (shape: 6,6)
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine significant_difference

Parameters
----------
value1 : float
    First value.

value2 : float
    Second value.

abs_tol : float, optional
    Absolute tolerance. Default is 0.

rel_tol : float, optional
    Relative tolerance. Default is 0.

Returns
-------
is_different : bool
    Set True if the difference is significant. False otherwise.
)"""
  );
  m.def(
      "skip_ele_blender",
      &Bmad::skip_ele_blender,
      py::arg("ele"),
      py::arg("skip"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine skip_ele_blender

Parameters
----------
ele : EleStruct

skip : bool
)"""
  );
  m.def(
      "slice_lattice",
      &Bmad::slice_lattice,
      py::arg("lat"),
      py::arg("ele_list"),
      py::arg("do_bookkeeping") = py::none(),
      py::arg("set_phase_zero") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine slice_lattice

Parameters
----------
lat : LatStruct
    Lattice to slice.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with unwanted elements sliced out.

ele_list : str
    List of elements to retain. See the documentation for the lat_ele_locator routine for the syntax of the
    list.

do_bookkeeping : bool, optional
    Default is True. If false, the calling routine is responsible for: * Modifying lat.particle_start if
    needed. * Calculating Twiss functions.

set_phase_zero : bool, optional
    Default is True. Set betatron phase to zero?

Returns
-------
error : bool
    Set True if there is an error Set False if not.
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
      py::call_guard<py::gil_scoped_release>(),
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
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Position after kick.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the edge.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix with edge kick added on.

make_matrix : bool, optional
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sol_quad_mat6_calc

Parameters
----------
ks_in : float

k1_in : float

tilt : float
    quadrupole tilt.

length : float
    Sol_quad length.

ele : EleStruct
    Sol_quad element.

orbit : CoordStruct
    Orbit at beginning of the sol_quad.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix up to the sol_quad.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix includeing the sol_quad.

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
      py::call_guard<py::gil_scoped_release>(),
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

args : 1D array of float (shape: 1:8)
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
      py::call_guard<py::gil_scoped_release>(),
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

args : 1D array of float (shape: 1:8)
    Parameters.  See psi_prime comments for details.

t : 1D array of float
    Array of times from t0 to t1

p : 1D array of float
    Array of psi evaluated at t(:)
)"""
  );
  m.def(
      "sort_complex_taylor_terms",
      &Bmad::sort_complex_taylor_terms,
      py::arg("complex_taylor_in"),
      py::call_guard<py::gil_scoped_release>(),
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
  py::class_<Bmad::SpinDnDpzFromMat8, std::unique_ptr<Bmad::SpinDnDpzFromMat8>>(
      m,
      "SpinDnDpzFromMat8",
      "spin_dn_dpz_from_mat8 return type"
  )
      .def_readonly("error", &Bmad::SpinDnDpzFromMat8::error)
      .def_readonly("dn_dpz", &Bmad::SpinDnDpzFromMat8::dn_dpz)
      .def("__len__", [](const Bmad::SpinDnDpzFromMat8 &) { return 2; })
      .def("__getitem__", [](const Bmad::SpinDnDpzFromMat8 &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.error);
        if (i == 1)
          return py::cast(s.dn_dpz);
        throw py::index_error();
      });
  m.def(
      "spin_dn_dpz_from_mat8",
      &Bmad::spin_dn_dpz_from_mat8,
      py::arg("mat_1turn"),
      py::arg("dn_dpz_partial") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_dn_dpz_from_mat8

Parameters
----------
mat_1turn : 2D array of float (shape: 8,8)
    Spin-orbital matrix.

dn_dpz_partial : 2D array of float (shape: 3,3), optional
    dn_dpz_partial(i,:) is dn_dpz with only one osccilation mode "excited". So dn_dpz_partial(1,:) represents
    a-mode excitation, etc.

Returns
-------
error : bool
    Set True if there is an error. False otherwise.

dn_dpz : 1D array of float (shape: 3)
    dn_dpz (l,n,m) coordinates.
)"""
  );
  py::class_<Bmad::SpinDnDpzFromQmap, std::unique_ptr<Bmad::SpinDnDpzFromQmap>>(
      m,
      "SpinDnDpzFromQmap",
      "spin_dn_dpz_from_qmap return type"
  )
      .def_readonly("error", &Bmad::SpinDnDpzFromQmap::error)
      .def_readonly("dn_dpz", &Bmad::SpinDnDpzFromQmap::dn_dpz)
      .def("__len__", [](const Bmad::SpinDnDpzFromQmap &) { return 2; })
      .def("__getitem__", [](const Bmad::SpinDnDpzFromQmap &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.error);
        if (i == 1)
          return py::cast(s.dn_dpz);
        throw py::index_error();
      });
  m.def(
      "spin_dn_dpz_from_qmap",
      &Bmad::spin_dn_dpz_from_qmap,
      py::arg("orb_mat"),
      py::arg("q_map"),
      py::arg("dn_dpz_partial"),
      py::arg("dn_dpz_partial2"),
      py::arg("n0") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_dn_dpz_from_qmap

Parameters
----------
orb_mat : 2D array of float (shape: 6,6)
    1-turn orbital matrix.

q_map : 2D array of float (shape: 0:3,0:6)
    1-turn spin linear quaternion map.

dn_dpz_partial : 2D array of float (shape: 3,3)
    ) is dn_dpz with only one osccilation mode "excited". So dn_dpz_partial(1,:) represents a-mode excitation,
    etc.

dn_dpz_partial2 : 2D array of float (shape: 3,3)
    ) is dn_dpz with only two osccilation modes "excited". So dn_dpz_partial(1,:) represents b-mode and c-mode
    excitation without the a-mode, etc.

n0 : 1D array of float (shape: 3), optional
    3,0).

Returns
-------
error : bool
    Set True if there is an error. False otherwise.

dn_dpz : 1D array of float (shape: 3)
    dn_dpz.
)"""
  );
  m.def(
      "spin_map1_normalize",
      &Bmad::spin_map1_normalize,
      py::arg("spin1"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_map1_normalize

Parameters
----------
spin1 : 2D array of float (shape: 0:3,0:6)
    Unnormalized spin map.
    This parameter is an input/output and is modified in-place.
    As an output, spin1: Normalized spin map.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_mat8_resonance_strengths

Parameters
----------
orb_evec : 1D array of complex (shape: 6)
    Orbital eigenvector.

mat8 : 2D array of float (shape: 6,6)
    Spin/orbital matrix.

Returns
-------
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_mat_to_eigen

Parameters
----------
orb_mat : 2D array of float (shape: 6,6)
    Orbital matrix.

spin_map : 2D array of float (shape: 0:3,0:6)
    Quaternion 0th & 1st order map.

Returns
-------
orb_eval : 1D array of complex (shape: 6)
    Eigenvalues.

orb_evec : 2D array of complex (shape: 6,6)
    Orbital eigenvectors. orb_evec(j,:) is the j^th vector.

n0 : 1D array of float (shape: 3)
    n_0 invariant spin

spin_evec : 2D array of complex (shape: 6,3)
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
      py::arg("omega"),
      py::arg("phase_space_coords") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_omega

Parameters
----------
field : EmFieldStruct

orbit : CoordStruct

sign_z_vel : int

omega : 1D array of float (shape: 3)

phase_space_coords : bool, optional
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_quat_resonance_strengths

Parameters
----------
orb_evec : 1D array of complex (shape: 6)
    Orbital eigenvector.

spin_q : 2D array of float (shape: 0:3,0:6)
    First order spin map.

Returns
-------
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spin_taylor_to_linear

Parameters
----------
spin_taylor : 1D array of TaylorStruct (shape: 0:3)
    Taylor spin map.

normalize : bool
    If True, normalize the linear map.

dref_orb : 1D array of float (shape: 6)
    Change in Reference orbit: output_map1_ref - input_taylor_ref.

is_on : bool
    Is map turned on? If not spin_map1 will be the unit map.

Returns
-------
spin_map1 : 2D array of float (shape: 0:3,0:6)
    First order spin map.
)"""
  );
  m.def(
      "spinor_to_polar",
      &Bmad::spinor_to_polar,
      py::arg("spinor"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spinor_to_polar

Parameters
----------
spinor : 1D array of complex (shape: 2)
    Spinor

Returns
-------
polar : SpinPolarStruct
    The resultant Unitary Vector in polar coordinates
)"""
  );
  m.def(
      "spinor_to_vec",
      &Bmad::spinor_to_vec,
      py::arg("spinor"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spinor_to_vec

Parameters
----------
spinor : 1D array of complex (shape: 2)
    Spinor

Returns
-------
vec : 1D array of float (shape: 3)
    spin vector in cartesian coordinates
)"""
  );
  m.def(
      "spline_fit_orbit",
      &Bmad::spline_fit_orbit,
      py::arg("start_orb"),
      py::arg("end_orb"),
      py::arg("spline_x"),
      py::arg("spline_y"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine spline_fit_orbit

Parameters
----------
start_orb : CoordStruct
    Starting coords.

end_orb : CoordStruct
    Ending coords.

spline_x : 1D array of float (shape: 0:3)
    Spline coefs for the horizontal trajectory.

spline_y : 1D array of float (shape: 0:3)
    Spline coefs for vertical trajectory.
)"""
  );
  m.def(
      "split_expression_string",
      &Bmad::split_expression_string,
      py::arg("expr"),
      py::arg("width"),
      py::arg("indent"),
      py::arg("break_str") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine split_expression_string (expr, width, indent, lines)

Routine to break an expression into a number of lines for a nicer display.
Used when printing expressions.

Parameters
----------
expr : str
    String containing the expression.

width : int
    Maximum width of split expression.

indent : int
    If positive: Number of spaces to indent for every line after the first. If negative: No indentation but
    first line is shortened by |indent|.

break_str : str, optional
    If present, only break lines at places where this string is.

Returns
-------
lines : 1D array of str
    Split expression.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine split_lat

Parameters
----------
lat : LatStruct
    Original lat structure.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Modified lat structure.

s_split : float
    Position at which lat.branch(ix_branch) is to be split.

ix_branch : int
    Index of lat.branch(:) to use.

add_suffix : bool, optional
    If True (default) add '#1' and '#2" suffixes to the split elements.

check_sanity : bool, optional
    If True (default) then call lat_sanity_check after the split to make sure everything is ok.

save_null_drift : bool, optional
    Save a copy of a drift to be split as a null_ele? This is useful when superpositions are done. See
    add_superimpose for more info. Default is False.

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

Returns
-------
ix_split : int
    Index of element just before the s = s_split point.

split_done : bool
    True if lat was split.

err_flag : bool, optional
    Set true if there is an error, false otherwise.
)"""
  );
  m.def(
      "sprint_spin_taylor_map",
      &Bmad::sprint_spin_taylor_map,
      py::arg("ele"),
      py::arg("start_orbit") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sprint_spin_taylor_map

Parameters
----------
ele : EleStruct
    Element to form map for.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with map.

start_orbit : 1D array of float (shape: 6), optional
    Reference orbit for the map. Default is zero orbit.
)"""
  );
  m.def(
      "sr_longitudinal_wake_particle",
      &Bmad::sr_longitudinal_wake_particle,
      py::arg("ele"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine sr_longitudinal_wake_particle (ele, orbit)

Routine to apply the short-range wake longitudinal component kick to a particle and then add
to the existing longitudinal wake the contribution from the particle.

Parameters
----------
ele : EleStruct
    Element with wakes.

orbit : CoordStruct
    Particle coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: coords after the kick.
)"""
  );
  m.def(
      "sr_transverse_wake_particle",
      &Bmad::sr_transverse_wake_particle,
      py::arg("ele"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine sr_transverse_wake_particle (ele, orbit)

Subroutine to apply the short-range wake transverse component of the kick to a particle and then add
to the existing transverse wake the contribution from the particle.

Parameters
----------
ele : EleStruct
    Element with wakes.

orbit : CoordStruct
    Starting particle coords.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending particle coords.
)"""
  );
  m.def(
      "sr_z_long_wake",
      &Bmad::sr_z_long_wake,
      py::arg("ele"),
      py::arg("bunch"),
      py::arg("z_ave"),
      py::call_guard<py::gil_scoped_release>(),
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
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("n_slices_gen_opt") = py::none(),
      py::arg("n_slices_sxt_opt") = py::none(),
      py::arg("chrom_set_x_opt") = py::none(),
      py::arg("chrom_set_y_opt") = py::none(),
      py::arg("weight_in") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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

var_indexes : 1D array of int
    indexes in lat.ele that are K2 variables.  Must be sorted smallest index to largest index.

n_slices_gen_opt : int, optional
    number of times to slice elements other than sextupoles.  Default is 10.

n_slices_sxt_opt : int, optional
    nubmer of times to slice sextupoles.  Default is 20.

chrom_set_x_opt : float, optional
    what to set x chromaticity to.  Default zero.

chrom_set_y_opt : float, optional
    what to set y chromaticity to.  Default zero.

weight_in : 1D array of float (shape: 10), optional
    moment weights. Terms are: [wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002, wgt_h21000,
    wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200, If present, any terms equal to zero are given default
    values which is 1.0e4 for wgt_chrom_x and wgt_chrom_y and is 1.0 for everything else.

Returns
-------
ls_soln : 1D array of float
    contains K2 for the indexes in var_indexes
)"""
  );
  m.def(
      "start_branch_at",
      &Bmad::start_branch_at,
      py::arg("lat"),
      py::arg("ele_start"),
      py::arg("move_end_marker"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine start_branch_at

Parameters
----------
lat : LatStruct
    Lattice to modify.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Modified lattice.

ele_start : str
    Start element. Ele_start will identify the lattice branch to modify.

move_end_marker : bool
    If True then the end marker (if it is present) will be shifted like any other element. False means that
    the end marker will stay at the end.

Returns
-------
error : bool
    Set True if there is an error Set False if not.
)"""
  );
  m.def(
      "stream_ele_end",
      &Bmad::stream_ele_end,
      py::arg("physical_end"),
      py::arg("ele_orientation"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine stream_ele_end

Parameters
----------
physical_end : int
    entrance_end$, exit_end$, surface$, etc.

ele_orientation : int
    Either 1 = Normal or -1 = element reversed.

Returns
-------
stream_end : int
    upstream_end$, downstream_end$, or set equal to physical_end if physical_end is neither entrance_end$ nor
    exit_end$
)"""
  );
  m.def(
      "string_attrib",
      &Bmad::string_attrib,
      py::arg("attrib_name"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine string_attrib (attrib_name, ele, attrib_value)

Routine to return the value of a string attribute of a lattice element.
This routine is useful when attrib_name is specified by the program user.

For example:
  call string_attrib ('NAME', ele, attrib_value)  ! Will return attrib_value = ele%name

Parameters
----------
attrib_name : str
    Name of the type of element attribute.

ele : EleStruct
    Lattice element.

Returns
-------
attrib_value : str
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine strong_beam_sigma_calc

Parameters
----------
ele : EleStruct
    Beambeam element.

s_pos : float
    Longitudinal position in lab coords of slice (used with hourglass effect correction).

Returns
-------
sigma : 1D array of float (shape: 2)
    Strong beam x,y sigmas.

bbi_const : float
    BBI kick scale factor.

dsigma_ds : 1D array of float (shape: 2)
    sig_x and sig_y longitudinal derivatives.
)"""
  );
  m.def(
      "strong_beam_strength",
      &Bmad::strong_beam_strength,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine strong_beam_strength

Parameters
----------
ele : EleStruct
    Beambeam element.

Returns
-------
strength : float
    Strong beam strength.
)"""
  );
  py::class_<Bmad::SurfaceGridDisplacement, std::unique_ptr<Bmad::SurfaceGridDisplacement>>(
      m,
      "SurfaceGridDisplacement",
      "surface_grid_displacement return type"
  )
      .def_readonly("err_flag", &Bmad::SurfaceGridDisplacement::err_flag)
      .def_readonly("z", &Bmad::SurfaceGridDisplacement::z)
      .def_readonly("dz_dxy", &Bmad::SurfaceGridDisplacement::dz_dxy)
      .def("__len__", [](const Bmad::SurfaceGridDisplacement &) { return 3; })
      .def("__getitem__", [](const Bmad::SurfaceGridDisplacement &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.z);
        if (i == 2)
          return py::cast(s.dz_dxy);
        throw py::index_error();
      });
  m.def(
      "surface_grid_displacement",
      &Bmad::surface_grid_displacement,
      py::arg("ele"),
      py::arg("x"),
      py::arg("y"),
      py::arg("extend_grid") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine surface_grid_displacement (ele, x, y, err_flag, z, dz_dxy, extend_grid)

Routine to add in the z displacement defined by the grid

Parameters
----------
ele : EleStruct
    Element containing the grid

x : float
    Photon coords at surface.

y : float
    Photon coords at surface.

extend_grid : bool, optional
    If (x,y) past grid pretend (x,y) is at grid boundary. Default is False.

Returns
-------
err_flag : bool
    Set True if there is a problem.

z : float
    surface height at (x, y).

dz_dxy : 1D array of float (shape: 2), optional
    Surface slope at (x, y).
)"""
  );
  py::class_<Bmad::SwitchAttribValueName, std::unique_ptr<Bmad::SwitchAttribValueName>>(
      m,
      "SwitchAttribValueName",
      "switch_attrib_value_name return type"
  )
      .def_readonly("is_default", &Bmad::SwitchAttribValueName::is_default)
      .def_readonly("name_list", &Bmad::SwitchAttribValueName::name_list)
      .def_readonly("attrib_val_name", &Bmad::SwitchAttribValueName::attrib_val_name)
      .def("__len__", [](const Bmad::SwitchAttribValueName &) { return 3; })
      .def("__getitem__", [](const Bmad::SwitchAttribValueName &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.is_default);
        if (i == 1)
          return py::cast(s.name_list);
        if (i == 2)
          return py::cast(s.attrib_val_name);
        throw py::index_error();
      });
  m.def(
      "switch_attrib_value_name",
      &Bmad::switch_attrib_value_name,
      py::arg("attrib_name"),
      py::arg("attrib_value"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function switch_attrib_value_name (attrib_name, attrib_value, ele,
                                    is_default, name_list) result (attrib_val_name)

Routine to return the name corresponding to the value of a given switch attribute.

This routine is for "switch" and "species" attributes. For example, the "aperture_type" attribute
can have value names of "Entrance_End", "Exit_End", etc.

Optionally, this routine can determine if the attribute value corresponds
to the default value. That is, the value that the attribute would have if
not specified in the lattice file.

Use the routine attribute_type to first test if the type of the attribute
corresponds to is_switch$.

Parameters
----------
attrib_name : str
    Name of the attribute. Must be upper case.

attrib_value : float
    Value of the attribute.

ele : EleStruct
    Lattice element that the attribute is contained in. Generally only needed to determine the default value.

Returns
-------
attrib_val_name : str
    Name corresponding to the value. Set to null_name$ if there is a problem.

is_default : bool, optional
    If True then the value of the attiribute corresponds to the default value. If this argument is present,
    the ele argument must also be present.

name_list : 1D array of str, optional
    List of names the switch can take. Deallocated if there is an error.
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
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine symp_lie_bmad

Parameters
----------
ele : EleStruct
    Element with transfer matrix
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with transfer matrix.

param : LatParamStruct
    Parameters are needed for some elements.

orbit : CoordStruct
    Coordinates at the beginning of element.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Coordinates at the end of element.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix propagated through the element.

make_matrix : bool, optional
    If True then make the 6x6 transfer matrix.

offset_ele : bool, optional
    Offset the element using ele.value(x_offset$), etc. Default is True.

Returns
-------
track : TrackStruct, optional
    Structure holding the track information. When tracking through multiple elements, the trajectory in an
    element is appended to the existing trajectory. To reset: Set track.n_pt = -1.
)"""
  );
}
