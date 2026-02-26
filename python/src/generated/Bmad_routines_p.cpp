#include "pybmad/generated/Bmad_routines_p.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_p(py::module &m) {
  m.def(
      "parse_cartesian_map",
      &Bmad::parse_cartesian_map,
      py::arg("ct_map"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parse_cartesian_map (ct_map, ele, lat, delim, delim_found, err_flag)

Subroutine to parse a "cartesian_map = {}" construct

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is private to bmad_parser_mod.
This must read in:
{type = ,
   dr = ,
   r0 = ,
   pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) )
   .
   .
   . ) },
)"""
  );
  m.def(
      "parse_cylindrical_map",
      &Bmad::parse_cylindrical_map,
      py::arg("cl_map"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine parse_cylindrical_map

Parameters
----------
cl_map : CylindricalMapStruct

ele : EleStruct

lat : LatStruct

delim : str

delim_found : bool

err_flag : bool
)"""
  );
  m.def(
      "parse_gen_grad_map",
      &Bmad::parse_gen_grad_map,
      py::arg("gg_map"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parse_gen_grad_map (gg_map, ele, lat, delim, delim_found, err_flag)

Subroutine to parse a "gen_grad_map = {}" construct
)"""
  );
  m.def(
      "parse_grid_field",
      &Bmad::parse_grid_field,
      py::arg("g_field"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine parse_grid_field

Parameters
----------
g_field : GridFieldStruct

ele : EleStruct

lat : LatStruct

delim : str

delim_found : bool

err_flag : bool
)"""
  );
  m.def(
      "parse_integer_list",
      &Bmad::parse_integer_list,
      py::arg("err_str"),
      py::arg("lat"),
      py::arg("int_array"),
      py::arg("exact_size"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("is_ok"),
      py::arg("open_delim") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_delim") = py::none(),
      py::arg("default_value") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parse_integer_list (err_str, lat, int_array, exact_size, delim, delim_found, open_delim,
                                      separator, close_delim, default_value) result (is_ok)

Routine to parse a list of integers of the form:
   open_delim integer_1 separator integer_2 . . . close_delim
Example:   "(1.2, 2.3, 4.4, 8.5)"

Similar to parse_integer_list2 except does not use allocatable array.
See parse_integer_list2 for more details
)"""
  );
  py::class_<Bmad::ParseIntegerList2, std::unique_ptr<Bmad::ParseIntegerList2>>(
      m,
      "ParseIntegerList2",
      "parse_integer_list2 return type"
  )
      .def_readonly("num_found", &Bmad::ParseIntegerList2::num_found)
      .def_readonly("delim", &Bmad::ParseIntegerList2::delim)
      .def_readonly("delim_found", &Bmad::ParseIntegerList2::delim_found)
      .def_readonly("is_ok", &Bmad::ParseIntegerList2::is_ok)
      .def("__len__", [](const Bmad::ParseIntegerList2 &) { return 4; })
      .def("__getitem__", [](const Bmad::ParseIntegerList2 &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.num_found);
        if (i == 1)
          return py::cast(s.delim);
        if (i == 2)
          return py::cast(s.delim_found);
        if (i == 3)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "parse_integer_list2",
      &Bmad::parse_integer_list2,
      py::arg("err_str"),
      py::arg("lat"),
      py::arg("int_array"),
      py::arg("num_expected") = py::none(),
      py::arg("open_delim") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_delim") = py::none(),
      py::arg("default_value") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parse_integer_list2 (err_str, lat, int_array, num_found, delim, delim_found, num_expected,
                                       open_delim, separator, close_delim, default_value) result (is_ok)

Routine to parse a list of integers of the form
   open_delim integer_1 separator integer_2 . . . close_delim
Example:   (1, 2, 4, 8)

Parameters
----------
err_str : str
    Error string to print if there is an error.

lat : LatStruct
    lattice

int_array : 1D array of int
    the array to be read in Optional:
    This parameter is an input/output and is modified in-place.
    As an output, int_array: Array of values.

Returns
-------
num_found : int
    number of elements.

delim : str
    Delimiter found where the parsing of the input line stops.

delim_found : bool
    Delimiter found? False if end of input command.

is_ok : bool
    Set True if everything is ok.
)"""
  );
  py::class_<Bmad::ParseRealList, std::unique_ptr<Bmad::ParseRealList>>(
      m,
      "ParseRealList",
      "parse_real_list return type"
  )
      .def_readonly("delim", &Bmad::ParseRealList::delim)
      .def_readonly("delim_found", &Bmad::ParseRealList::delim_found)
      .def_readonly("num_found", &Bmad::ParseRealList::num_found)
      .def("__len__", [](const Bmad::ParseRealList &) { return 3; })
      .def("__getitem__", [](const Bmad::ParseRealList &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.num_found);
        throw py::index_error();
      });
  m.def(
      "parse_real_list",
      &Bmad::parse_real_list,
      py::arg("lat"),
      py::arg("err_str"),
      py::arg("real_array"),
      py::arg("exact_size"),
      py::arg("is_ok"),
      py::arg("open_delim") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_delim") = py::none(),
      py::arg("default_value") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parse_real_list (lat, err_str, real_array, exact_size, delim, delim_found, open_delim,
                               separator, close_delim, default_value, num_found) result (is_ok)

Routine to parse a list of reals of the form:
   open_delim real_1 separator real_2 . . . close_delim
Example:   "(1.2, 2.3, 4.4, 8.5)"

Similar to parse_real_list2 except does not use allocatable array.
Also see: parse_real_matrix.

Parameters
----------
lat : LatStruct
    Lattice

err_str : str
    Error string to print if there is an error.

real_array : 1D array of float

exact_size : bool

open_delim : str, optional

separator : str, optional

close_delim : str, optional

default_value : float, optional

Returns
-------
delim : str

delim_found : bool

num_found : int, optional
)"""
  );
  py::class_<Bmad::ParseRealList2, std::unique_ptr<Bmad::ParseRealList2>>(
      m,
      "ParseRealList2",
      "parse_real_list2 return type"
  )
      .def_readonly("num_found", &Bmad::ParseRealList2::num_found)
      .def_readonly("delim", &Bmad::ParseRealList2::delim)
      .def_readonly("delim_found", &Bmad::ParseRealList2::delim_found)
      .def_readonly("is_ok", &Bmad::ParseRealList2::is_ok)
      .def("__len__", [](const Bmad::ParseRealList2 &) { return 4; })
      .def("__getitem__", [](const Bmad::ParseRealList2 &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.num_found);
        if (i == 1)
          return py::cast(s.delim);
        if (i == 2)
          return py::cast(s.delim_found);
        if (i == 3)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "parse_real_list2",
      &Bmad::parse_real_list2,
      py::arg("lat"),
      py::arg("err_str"),
      py::arg("real_array"),
      py::arg("num_expected") = py::none(),
      py::arg("open_brace") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_brace") = py::none(),
      py::arg("default_value") = py::none(),
      py::arg("single_value") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parse_real_list2 (lat, err_str, real_array, num_found, delim, delim_found, num_expected,
                           open_delim, separator, close_delim, default_value, single_value) result (is_ok)

Routine to parse a list of reals of the form:
   open_brace real_1 separator real_2 . . . close_brace
Example:   "(1.2, 2.3, 4.4, 8.5)"

Also see:
  pase_real_list
  parse_real_matrix.

Parameters
----------
lat : LatStruct
    lattice

err_str : str
    Error string to print if there is an error.

real_array : 1D array of float
    the array to be read in
    This parameter is an input/output and is modified in-place.
    As an output, real_array: Array of values

Returns
-------
num_found : int
    number of elements

delim : str
    Delimiter found where the parsing of the input line stops.

delim_found : bool
    Stopping delimiter found? False if end of input command.

is_ok : bool
    Set True if everything is ok
)"""
  );
  m.def(
      "parser_add_constant",
      &Bmad::parser_add_constant,
      py::arg("word"),
      py::arg("lat"),
      py::arg("redef_is_error"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine parser_add_constant

Parameters
----------
word : str

lat : LatStruct

redef_is_error : bool
)"""
  );
  m.def(
      "parser_call_check",
      &Bmad::parser_call_check,
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("call_found"),
      py::arg("err_flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_call_check(word, ix_word, delim, delim_found, call_found, err_flag))

Routine to check if there is a "call::XXX" construct in the input stream.
)"""
  );
  py::class_<Bmad::ParserFastComplexRead, std::unique_ptr<Bmad::ParserFastComplexRead>>(
      m,
      "ParserFastComplexRead",
      "parser_fast_complex_read return type"
  )
      .def_readonly("delim", &Bmad::ParserFastComplexRead::delim)
      .def_readonly("is_ok", &Bmad::ParserFastComplexRead::is_ok)
      .def("__len__", [](const Bmad::ParserFastComplexRead &) { return 2; })
      .def("__getitem__", [](const Bmad::ParserFastComplexRead &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "parser_fast_complex_read",
      &Bmad::parser_fast_complex_read,
      py::arg("cmplx_vec"),
      py::arg("ele"),
      py::arg("err_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parser_fast_complex_read (cmplx_vec, ele, delim, err_str)  result (is_ok)

Routine to read an array of complex numbers.

This routine assumes that the array values are pure numbers in the form "<re>" or "(<re> <im>)"
where <re> and <im> are real numbers (not expressions) and there are no commas except possibly
at the end of the array.

Parameters
----------
cmplx_vec : 1D array of complex
    Complex vector.

ele : EleStruct
    Lattice element associated with the array. Used for error messages.

err_str : str
    String used when printing error messages identifying where in the lattice file the error is occuring.

Returns
-------
delim : str
    Delimitor at end of array. Must be "," or "}"

is_ok : bool
    True if everything OK. False otherwise.
)"""
  );
  m.def(
      "parser_fast_integer_read",
      &Bmad::parser_fast_integer_read,
      py::arg("int_vec"),
      py::arg("ele"),
      py::arg("delim_wanted"),
      py::arg("err_str"),
      py::arg("is_ok"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parser_fast_integer_read (int_vec, ele, delim_wanted, err_str)  result (is_ok)
)"""
  );
  py::class_<Bmad::ParserFastRealRead, std::unique_ptr<Bmad::ParserFastRealRead>>(
      m,
      "ParserFastRealRead",
      "parser_fast_real_read return type"
  )
      .def_readonly("delim", &Bmad::ParserFastRealRead::delim)
      .def_readonly("n_real", &Bmad::ParserFastRealRead::n_real)
      .def_readonly("is_ok", &Bmad::ParserFastRealRead::is_ok)
      .def("__len__", [](const Bmad::ParserFastRealRead &) { return 3; })
      .def("__getitem__", [](const Bmad::ParserFastRealRead &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.n_real);
        if (i == 2)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "parser_fast_real_read",
      &Bmad::parser_fast_real_read,
      py::arg("real_vec"),
      py::arg("ele"),
      py::arg("end_delims"),
      py::arg("err_str"),
      py::arg("exact_size") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function parser_fast_real_read (real_vec, ele, end_delims, delim, err_str, exact_size, n_real)  result (is_ok)

Routine to read an array of real numbers.

This routine assumes that the array values are pure numbers in the form "<re1> <re2> ...,"
where <re1>, <re2>, etc. are real numbers (not expressions) and there are no commas except possibly,
at the end of the array.

Note: if end_delim is "," and next character is a delim but not ",", the next character is taken as the delim.

Parameters
----------
real_vec : 1D array of float
    Real vector.

ele : EleStruct
    Lattice element associated with the array. Used for error messages.

end_delims : str
    List of possible ending delimitors.

err_str : str
    String used when printing error messages identifying where in the lattice file the error is occuring.

exact_size : bool, optional
    If True (default), number of values must match real_vec size.

Returns
-------
delim : str
    Delimitor at end of array.

is_ok : bool
    True if everything OK. False otherwise.

n_real : int, optional
    Number of elements found.
)"""
  );
  m.def(
      "parser_file_stack",
      &Bmad::parser_file_stack,
      py::arg("how"),
      py::arg("file_name_in") = py::none(),
      py::arg("finished") = py::none(),
      py::arg("err") = py::none(),
      py::arg("open_file") = py::none(),
      py::arg("abort_on_open_error") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_file_stack (how, file_name_in, finished, err, open_file, abort_on_open_error)

Subroutine to keep track of the files that are opened for reading.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  m.def(
      "parser_get_integer",
      &Bmad::parser_get_integer,
      py::arg("int_val"),
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err"),
      py::arg("str1") = py::none(),
      py::arg("str2") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine parser_get_integer

Parameters
----------
int_val : int

word : str

ix_word : int

delim : str

delim_found : bool

err : bool

str1 : str, optional

str2 : str, optional
)"""
  );
  m.def(
      "parser_get_logical",
      &Bmad::parser_get_logical,
      py::arg("attrib_name"),
      py::arg("this_logic"),
      py::arg("ele_name"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine parser_get_logical

Parameters
----------
attrib_name : str

this_logic : bool

ele_name : str

delim : str

delim_found : bool

err : bool
)"""
  );
  m.def(
      "parser_identify_fork_to_element",
      &Bmad::parser_identify_fork_to_element,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_identify_fork_to_element (lat)

Routine to identify the elements the forks in a lattice are branching to.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  m.def(
      "parser_init_custom_elements",
      &Bmad::parser_init_custom_elements,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_init_custom_elements (lat)
)"""
  );
  m.def(
      "parser_print_line",
      &Bmad::parser_print_line,
      py::arg("lat"),
      py::arg("end_of_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_print_line(end_of_file)

This routine is called when a print statement is found in the lattice file.
)"""
  );
  m.def(
      "parser_read_lr_wake",
      &Bmad::parser_read_lr_wake,
      py::arg("ele"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_read_lr_wake (ele, delim, delim_found, err_flag)

Subroutine to read in a long-range wake field from an external file.
This subroutine is used by bmad_parser and bmad_parser2.

Parameters
----------
ele : EleStruct
    Element containing wake structure.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with wake information.
)"""
  );
  m.def(
      "parser_read_old_format_lr_wake",
      &Bmad::parser_read_old_format_lr_wake,
      py::arg("ele"),
      py::arg("lr_file_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_read_old_format_lr_wake (ele, lr_file_name)

Subroutine to read in a long-range wake field from an external file.
This subroutine is used by bmad_parser and bmad_parser2.

Parameters
----------
ele : EleStruct
    Element containing wake structure.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with wake information.

lr_file_name : str
    Name of long-range wake field file.
)"""
  );
  m.def(
      "parser_read_old_format_sr_wake",
      &Bmad::parser_read_old_format_sr_wake,
      py::arg("ele"),
      py::arg("sr_file_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_read_old_format_sr_wake (ele, sr_file_name)

Subroutine to read in a short-range wake field from an external file.
This subroutine is used by bmad_parser and bmad_parser2.

Parameters
----------
ele : EleStruct
    Element containing wake structure.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with wake information.

sr_file_name : str
    Name of short-range wake field file.
)"""
  );
  m.def(
      "parser_read_sr_wake",
      &Bmad::parser_read_sr_wake,
      py::arg("ele"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_read_sr_wake (ele, delim, delim_found, err_flag)

Subroutine to read in a short-range wake field.
This subroutine is used by bmad_parser and bmad_parser2.

Parameters
----------
ele : EleStruct
    Element containing wake structure.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with wake information.
)"""
  );
  m.def(
      "parser_transfer_control_struct",
      &Bmad::parser_transfer_control_struct,
      py::arg("con_in"),
      py::arg("lord"),
      py::arg("ix_var"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine parser_transfer_control_struct (con_in, con_out, lord, ix_var)

Routine to transfer the information from an input control_struct (which stores
the user input parameters) to a control_struct that will be stored in the lat%control
or lord%control%ramp for a ramper.

Parameters
----------
con_in : ControlStruct
    Input control structure.

lord : EleStruct
    Lord element associated with the control_struct.

ix_var : int
    If an expression stack evaluates to a constant, this routine will modify the expression stack to evaluate
    to the value of: lord.control.var(ix_var) * constant

Returns
-------
con_out : ControlStruct
    Output control structure.
)"""
  );
  m.def(
      "particle_in_global_frame",
      &Bmad::particle_in_global_frame,
      py::arg("orb"),
      py::arg("branch"),
      py::arg("in_time_coordinates") = py::none(),
      py::arg("in_body_frame") = py::none(),
      py::arg("w_mat_out") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function particle_in_global_frame (orb, in_time_coordinates, in_body_frame, w_mat_out) result (particle)

Returns the particle in global time coordinates given is coordinates orb in lattice lat.

Parameters
----------
orb : CoordStruct
    particle in s-coordinates

branch : BranchStruct
    branch that contains branch.ele(orb.ix_ele)

in_time_coordinates : bool, optional
    Default is false. If true, orb will taken as in time coordinates.

in_body_frame : bool, optional
    Default is true. If false, ele offsets will be ignored.

Returns
-------
particle : CoordStruct
    particle in global time coordinates
)"""
  );
  m.def(
      "particle_is_moving_backwards",
      &Bmad::particle_is_moving_backwards,
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine particle_is_moving_backwards

Parameters
----------
orbit : CoordStruct
    Particle coordinates

Returns
-------
is_moving_backwards : bool
    True if moving backward. False otherwise.
)"""
  );
  m.def(
      "particle_is_moving_forward",
      &Bmad::particle_is_moving_forward,
      py::arg("orbit"),
      py::arg("dir") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine particle_is_moving_forward

Parameters
----------
orbit : CoordStruct
    Particle coordinates

dir : int, optional
    +1 if tracking forward(default) or -1 to return True if tracking backwards.

Returns
-------
is_moving_forward : bool
    True if moving forward. False otherwise.
)"""
  );
  m.def(
      "particle_rf_time",
      &Bmad::particle_rf_time,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("reference_active_edge") = py::none(),
      py::arg("s_rel") = py::none(),
      py::arg("time_coords") = py::none(),
      py::arg("rf_freq") = py::none(),
      py::arg("abs_time") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine particle_rf_time

Parameters
----------
orbit : CoordStruct
    Particle coordinates

ele : EleStruct
    Element being tracked through.

reference_active_edge : bool, optional
    If True, and ele is a rfcavity or lcavity, use the active edge (edge of the region with non-zero field) as
    the reference point.

s_rel : float, optional
    Longitudinal position relative to the upstream edge of the element. Needed for relative time tracking when
    the particle is inside the element. Default is 0.

time_coords : bool, optional
    Default False. If True then orbit is using time based phase space coordinates.

rf_freq : float, optional
    If present, the returned time shifted by an integer multiple of 1/rf_freq to be in the range
    [-1/2*rf_freq, 1/2*rf_freq]. This is useful to avoid round-off errors.

abs_time : bool, optional
    If False (default) use setting of bmad_com.absolute_time_tracking. If True, use absolute time instead of
    relative time.

Returns
-------
time : float
    Current time.
)"""
  );
  m.def(
      "patch_flips_propagation_direction",
      &Bmad::patch_flips_propagation_direction,
      py::arg("x_pitch"),
      py::arg("y_pitch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine patch_flips_propagation_direction

Parameters
----------
x_pitch : float
    Rotaion around y-axis

y_pitch : float
    Rotation around x-axis.

Returns
-------
is_flip : bool
    True if patch does a flip
)"""
  );
  m.def(
      "patch_length",
      &Bmad::patch_length,
      py::arg("patch"),
      py::arg("ref_coords") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine patch_length

Parameters
----------
patch : EleStruct
    Patch element.

ref_coords : int, optional
    Reference coords to use. entrance_end$, exit_end$ Default is nint(patch.value(ref_coords$)).

Returns
-------
length : float
    Length of patch.
)"""
  );
  py::class_<
      Bmad::PhotonAbsorptionAndPhaseShift,
      std::unique_ptr<Bmad::PhotonAbsorptionAndPhaseShift>>(
      m,
      "PhotonAbsorptionAndPhaseShift",
      "photon_absorption_and_phase_shift return type"
  )
      .def_readonly("absorption", &Bmad::PhotonAbsorptionAndPhaseShift::absorption)
      .def_readonly("phase_shift", &Bmad::PhotonAbsorptionAndPhaseShift::phase_shift)
      .def_readonly("err_flag", &Bmad::PhotonAbsorptionAndPhaseShift::err_flag)
      .def("__len__", [](const Bmad::PhotonAbsorptionAndPhaseShift &) { return 3; })
      .def("__getitem__", [](const Bmad::PhotonAbsorptionAndPhaseShift &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.absorption);
        if (i == 1)
          return py::cast(s.phase_shift);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "photon_absorption_and_phase_shift",
      &Bmad::photon_absorption_and_phase_shift,
      py::arg("material"),
      py::arg("Energy"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_absorption_and_phase_shift (material, Energy, absorption, phase_shift, err_flag)

Routine to calcualte the absorption and phase shift values for a photon with a given
energy going through a particular material.

Parameters
----------
material : str
    Material name.

Energy : float
    Photon energy (eV).

Returns
-------
absorption : float
    E_field ~ Exp(-absorption * length)

phase_shift : float
    E_field Phase shift (radians) per unit length relative to vacuum.

err_flag : bool
    Set true if material not recognized.
)"""
  );
  py::class_<
      Bmad::PhotonAddToDetectorStatistics,
      std::unique_ptr<Bmad::PhotonAddToDetectorStatistics>>(
      m,
      "PhotonAddToDetectorStatistics",
      "photon_add_to_detector_statistics return type"
  )
      .def_readonly("ix_pt", &Bmad::PhotonAddToDetectorStatistics::ix_pt)
      .def_readonly("iy_pt", &Bmad::PhotonAddToDetectorStatistics::iy_pt)
      .def("__len__", [](const Bmad::PhotonAddToDetectorStatistics &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonAddToDetectorStatistics &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ix_pt);
        if (i == 1)
          return py::cast(s.iy_pt);
        throw py::index_error();
      });
  m.def(
      "photon_add_to_detector_statistics",
      &Bmad::photon_add_to_detector_statistics,
      py::arg("orbit0"),
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("pixel_pt") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_add_to_detector_statistics (orbit0, orbit, ele, ix_pt, iy_pt, pixel_pt)

Routine to add photon statistics to the appropriate pixel of a "detector" grid.

It is assumed that track_to_surface has been called so that the photon is at the
detector surface and that orbit%vec(1) and %vec(3) are in surface coords (needed for curved detectors).

Parameters
----------
orbit0 : CoordStruct
    Photon coords at beginning of lattice

orbit : CoordStruct
    Photon coords at the detector.

ele : EleStruct
    Element with grid.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with updatted grid.

pixel_pt : PixelPtStruct, optional
    If present then use this grid point instead of the grid point determined by the (x, y) coords of the
    photon

Returns
-------
ix_pt : int, optional
    Index of upgraded ele.photon.surface.grid.pt(:,:) point. These arguments are not set if the pixel_pt
    argument is present.

iy_pt : int, optional
    Index of upgraded ele.photon.surface.grid.pt(:,:) point. These arguments are not set if the pixel_pt
    argument is present.
)"""
  );
  py::class_<Bmad::PhotonReflection, std::unique_ptr<Bmad::PhotonReflection>>(
      m,
      "PhotonReflection",
      "photon_reflection return type"
  )
      .def_readonly("graze_angle_out", &Bmad::PhotonReflection::graze_angle_out)
      .def_readonly("phi_out", &Bmad::PhotonReflection::phi_out)
      .def("__len__", [](const Bmad::PhotonReflection &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonReflection &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.graze_angle_out);
        if (i == 1)
          return py::cast(s.phi_out);
        throw py::index_error();
      });
  m.def(
      "photon_reflection",
      &Bmad::photon_reflection,
      py::arg("graze_angle_in"),
      py::arg("energy"),
      py::arg("surface"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_reflection (graze_angle_in, energy, surface, graze_angle_out, phi_out)

Routine to reflect a photon from a surface including both diffuse and specular reflections.

Parameters
----------
graze_angle_in : float
    Incident grazing (not polar) angle in radians.

energy : float
    Photon energy in eV.

surface : PhotonReflectSurfaceStruct
    surface info

Returns
-------
graze_angle_out : float
    graze_angle in radians.

phi_out : float
    Azimuthal angle in radians.
)"""
  );
  m.def(
      "photon_reflection_std_surface_init",
      &Bmad::photon_reflection_std_surface_init,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_reflection_std_surface_init (surface)

Routine to initialize the standard proton reflection probability tables.
The standard tables are for 10 nm C film on Al substrate.
The surface roughness for diffuse scattering is 200 nm and the
the surface roughness correlation length is 5.5 um.

Returns
-------
surface : PhotonReflectSurfaceStruct
    photon_reflect_surface_struct
)"""
  );
  py::class_<Bmad::PhotonReflectivity, std::unique_ptr<Bmad::PhotonReflectivity>>(
      m,
      "PhotonReflectivity",
      "photon_reflectivity return type"
  )
      .def_readonly("p_reflect", &Bmad::PhotonReflectivity::p_reflect)
      .def_readonly("rel_p_specular", &Bmad::PhotonReflectivity::rel_p_specular)
      .def("__len__", [](const Bmad::PhotonReflectivity &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonReflectivity &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.p_reflect);
        if (i == 1)
          return py::cast(s.rel_p_specular);
        throw py::index_error();
      });
  m.def(
      "photon_reflectivity",
      &Bmad::photon_reflectivity,
      py::arg("angle"),
      py::arg("energy"),
      py::arg("surface"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_reflectivity (angle, energy, surface, p_reflect, rel_p_specular)

Routine to evaluate the photon reflectivity.
  probability of absorption          = 1 - p_reflect
  probability of reflection          = p_reflect
  probability of specular reflection = p_reflect * rel_p_specular
  probability of diffuse reflection  = p_reflect * (1 - rel_p_specular)

Use photon_reflection_std_surface_init or read_surface_reflection_file to get surface info.

Parameters
----------
angle : float
    Incident grazing angle in radians.

energy : float
    Photon energy in eV.

surface : PhotonReflectSurfaceStruct
    surface info

Returns
-------
p_reflect : float
    Reflection probability.

rel_p_specular : float
    Relative specular reflection probability.
)"""
  );
  m.def(
      "photon_target_corner_calc",
      &Bmad::photon_target_corner_calc,
      py::arg("aperture_ele"),
      py::arg("x_lim"),
      py::arg("y_lim"),
      py::arg("z_lim"),
      py::arg("source_ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_target_corner_calc (aperture_ele, x_lim, y_lim, z_lim, source_ele, corner)

Routine to calculate the corner coords in the source_ele ref frame.

Parameters
----------
aperture_ele : EleStruct
    Element containing the aperture

x_lim : float
    Transverse corner points in aperture_ele coord frame.

y_lim : float
    Transverse corner points in aperture_ele coord frame.

source_ele : EleStruct
    Photon source element.

Returns
-------
corner : TargetPointStruct
    Corner coords in source_ele ref frame.
)"""
  );
  m.def(
      "photon_target_setup",
      &Bmad::photon_target_setup,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine photon_target_setup (ele)

Routine to calculate and store the parmeters needed for photon targeting.
This routine is called by Bmad parsing routines and is not meant for general use.

Photon initialization with targeting is done by the routine init_photon_from_a_photon_init_ele
Which is called by init_coord.

Parameters
----------
ele : EleStruct
    Source element to setup. Element will be of type: sample, diffraction_plate or photon_init.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Source element with target parameters setup in ele.photon.target.
)"""
  );
  m.def(
      "photon_type",
      &Bmad::photon_type,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function photon_type (ele) result (e_type)

Routine to return the type of photon to be tracked: coherent$ or incoherent$.

Parameters
----------
ele : EleStruct
    Element being tracked through.

Returns
-------
e_type : int
    coherent$ or incoherent$
)"""
  );
  m.def(
      "physical_ele_end",
      &Bmad::physical_ele_end,
      py::arg("track_end"),
      py::arg("orbit"),
      py::arg("ele_orientation"),
      py::arg("return_stream_end") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine physical_ele_end

Parameters
----------
track_end : int
    first_track_edge$, second_track_edge$, surface$, or in_between$

orbit : CoordStruct
    Particle position.

ele_orientation : int
    Either 1 = Normal or -1 = element reversed.

return_stream_end : bool, optional
    If True return the stream end instead of the physical end. Default is False.

Returns
-------
physical_end : int
    Return_stream_end ->  Possibilities False             ->  entrance_end$, exit_end$, surface$, or
    in_between$ True              ->  upstream_end$, downstream_end$
)"""
  );
  m.def(
      "point_photon_emission",
      &Bmad::point_photon_emission,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::arg("direction"),
      py::arg("max_target_area"),
      py::arg("w_to_surface") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine point_photon_emission (ele, param, orbit, direction, max_target_area, w_to_surface)

Routine to emit a photon from a point that may be on a surface.
If there is a downstream target, the emission calc will take this into account.

Parameters
----------
ele : EleStruct
    Emitting element.

param : LatParamStruct
    lattice parameters.

orbit : CoordStruct
    phase-space coords of photon. --   Will be in curved surface coords if there is a curved surface.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Final phase-space coords

direction : int
    +1 -> Emit in forward +z direction, -1 -> emit backwards.

max_target_area : float
    Area of the solid angle photons may be emitted over. max_target_area is used for normalizing the photon
    field. generally will be equal to twopi or fourpi.

w_to_surface : 2D array of float (shape: 3,3), optional
    Rotation matrix for curved surface.
)"""
  );
  m.def(
      "pointer_to_branch",
      py::overload_cast<EleStruct &>(&Bmad::pointer_to_branch),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_branch

Routine to return a pointer to the lattice branch associated with a given name
or a given element.

This routine is an overloaded name for:
  pointer_to_branch_given_ele (ele) result (branch_ptr))
  pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

The lattice branch *associated* with a given element is not necessarily the
branch where the element is *located*. For example, all lords live in branch #0.
But the branch associated with a super_lord element is the branch of its slaves.

To get the branch where the element is located, simply use ele%ix_branch.

Note: Result is ambiguous if ele argument is associated with multiple branches
which can happen, for example, with overlay elements.

Parameters
----------
ele : EleStruct
    Element contained in the branch.

Returns
-------
branch_ptr : BranchStruct, optional
    Pointer to the branch. Nullified if there is no associated branch.
)"""
  );
  m.def(
      "pointer_to_branch",
      py::overload_cast<std::string, LatStruct &, std::optional<bool>, std::optional<int>>(
          &Bmad::pointer_to_branch
      ),
      py::arg("branch_name"),
      py::arg("lat"),
      py::arg("parameter_is_branch0") = py::none(),
      py::arg("blank_branch") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_branch

Routine to return a pointer to the lattice branch associated with a given name
or a given element.

This routine is an overloaded name for:
  pointer_to_branch_given_ele (ele) result (branch_ptr))
  pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

The lattice branch *associated* with a given element is not necessarily the
branch where the element is *located*. For example, all lords live in branch #0.
But the branch associated with a super_lord element is the branch of its slaves.

To get the branch where the element is located, simply use ele%ix_branch.

Note: Result is ambiguous if ele argument is associated with multiple branches
which can happen, for example, with overlay elements.

Parameters
----------
branch_name : str
    May be a branch name or a branch index.

lat : LatStruct
    Lattice to search.

parameter_is_branch0 : bool, optional
    If True, 'PARAMETER' is taken to be an alternative name for branch(0). Default is False.

blank_branch : int, optional
    Branch index if branch_name = ''. Default is blank is an error.

Returns
-------
branch_ptr : BranchStruct, optional
    Pointer to the branch. Nullified if there is no associated branch.
)"""
  );
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatStruct &, int, std::optional<int>>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("ix_ele"),
      py::arg("ix_branch") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_ele (...)

Routine to return a pointer to an element.
pointer_to_ele is an overloaded name for:
    Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
    Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
    Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
    Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
element in lat.

Note that using ele_name to locate an element is potentially dangerous if there
are multiple elements that have the same name. Better in this case is to use:
  lat_ele_locator

Also see:
  pointer_to_slave
  pointer_to_lord

Parameters
----------
lat : LatStruct
    Lattice.

ix_ele : int
    Index of element in lat.branch(ix_branch).

ix_branch : int, optional
    Index of the lat.branch(:) containing the element.

Returns
-------
ele_ptr : EleStruct, optional
    Pointer to the element. Nullified if no match or error.
)"""
  );
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatStruct &, LatEleLocStruct &>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("ele_loc"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_ele (...)

Routine to return a pointer to an element.
pointer_to_ele is an overloaded name for:
    Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
    Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
    Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
    Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
element in lat.

Note that using ele_name to locate an element is potentially dangerous if there
are multiple elements that have the same name. Better in this case is to use:
  lat_ele_locator

Also see:
  pointer_to_slave
  pointer_to_lord

Parameters
----------
lat : LatStruct
    Lattice.

ele_loc : LatEleLocStruct
    Location identification.

Returns
-------
ele_ptr : EleStruct, optional
    Pointer to the element. Nullified if no match or error.
)"""
  );
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatStruct &, std::string>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("ele_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_ele (...)

Routine to return a pointer to an element.
pointer_to_ele is an overloaded name for:
    Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
    Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
    Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
    Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
element in lat.

Note that using ele_name to locate an element is potentially dangerous if there
are multiple elements that have the same name. Better in this case is to use:
  lat_ele_locator

Also see:
  pointer_to_slave
  pointer_to_lord

Parameters
----------
lat : LatStruct
    Lattice.

ele_name : str
    Name or index of element.

Returns
-------
ele_ptr : EleStruct, optional
    Pointer to the element. Nullified if no match or error.
)"""
  );
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatStruct &, EleStruct &>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("foreign_ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_ele (...)

Routine to return a pointer to an element.
pointer_to_ele is an overloaded name for:
    Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
    Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
    Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
    Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
element in lat.

Note that using ele_name to locate an element is potentially dangerous if there
are multiple elements that have the same name. Better in this case is to use:
  lat_ele_locator

Also see:
  pointer_to_slave
  pointer_to_lord

Parameters
----------
lat : LatStruct
    Lattice.

foreign_ele : EleStruct
    Lattice element in another lattice.

Returns
-------
ele_ptr : EleStruct, optional
    Pointer to the element. Nullified if no match or error.
)"""
  );
  py::class_<Bmad::PointerToElementAtS, std::unique_ptr<Bmad::PointerToElementAtS>>(
      m,
      "PointerToElementAtS",
      "pointer_to_element_at_s return type"
  )
      .def_readonly("err_flag", &Bmad::PointerToElementAtS::err_flag)
      .def_readonly("s_eff", &Bmad::PointerToElementAtS::s_eff)
      .def_readonly("position", &Bmad::PointerToElementAtS::position)
      .def_readonly("ele", &Bmad::PointerToElementAtS::ele)
      .def("__len__", [](const Bmad::PointerToElementAtS &) { return 4; })
      .def("__getitem__", [](const Bmad::PointerToElementAtS &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.s_eff);
        if (i == 2)
          return py::cast(s.position);
        if (i == 3)
          return py::cast(s.ele);
        throw py::index_error();
      });
  m.def(
      "pointer_to_element_at_s",
      &Bmad::pointer_to_element_at_s,
      py::arg("branch"),
      py::arg("s"),
      py::arg("choose_max"),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_element_at_s (branch, s, choose_max, err_flag, s_eff, position) result (ele)

Function to return a pointer to the element at position s.
That is, return ele => branch%ele(ix_ele) such that:
If choose_max = True:
    If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
    Else: branch%ele(ix_ele)%s_strat <= s < branch%ele(ix_ele)%s
If choose_max = False:
    If s = branch%ele(0): ix_ele = 0
    Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
    choose_max = True  => ix_ele = ix2
    choose_max = False => ix_ele = ix1

Also see: element_at_s

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
    See above.

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

ele : EleStruct, optional
    Pointer to element at s.
)"""
  );
  m.def(
      "pointer_to_fibre",
      &Bmad::pointer_to_fibre,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_fibre

Parameters
----------
ele : EleStruct
    Bmad element

Returns
-------
assoc_fibre : Fibre, optional
    Pointer to the associated fibre.
)"""
  );
  py::class_<Bmad::PointerToFieldEle, std::unique_ptr<Bmad::PointerToFieldEle>>(
      m,
      "PointerToFieldEle",
      "pointer_to_field_ele return type"
  )
      .def_readonly("dz_offset", &Bmad::PointerToFieldEle::dz_offset)
      .def_readonly("field_ele", &Bmad::PointerToFieldEle::field_ele)
      .def("__len__", [](const Bmad::PointerToFieldEle &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToFieldEle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.dz_offset);
        if (i == 1)
          return py::cast(s.field_ele);
        throw py::index_error();
      });
  m.def(
      "pointer_to_field_ele",
      &Bmad::pointer_to_field_ele,
      py::arg("ele"),
      py::arg("ix_field_ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_field_ele

Parameters
----------
ele : EleStruct
    Element with sum number of associated field elements.

ix_field_ele : int
    Index of the field element to point to. This index runs from 1 to num_field_eles(ele).

Returns
-------
dz_offset : float, optional
    Longitudinal offset of ele upstream edge from the field ele pointed to.

field_ele : EleStruct, optional
    Pointer to the field element with index ix_field_ele. Will point to null if ix_field_ele is out of range.
)"""
  );
  py::class_<Bmad::PointerToGirder, std::unique_ptr<Bmad::PointerToGirder>>(
      m,
      "PointerToGirder",
      "pointer_to_girder return type"
  )
      .def_readonly("ix_slave_back", &Bmad::PointerToGirder::ix_slave_back)
      .def_readonly("girder", &Bmad::PointerToGirder::girder)
      .def("__len__", [](const Bmad::PointerToGirder &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToGirder &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ix_slave_back);
        if (i == 1)
          return py::cast(s.girder);
        throw py::index_error();
      });
  m.def(
      "pointer_to_girder",
      &Bmad::pointer_to_girder,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_girder

Parameters
----------
ele : EleStruct
    Element to check.

Returns
-------
ix_slave_back : int, optional
    Index back to ele. That is, pointer_to_slave(girder, ix_slave_back) will point back to ele. Set to -1 if
    no girder present

girder : EleStruct, optional
    : Pointer to the girder. Null if ele is not girder supported.
)"""
  );
  py::class_<Bmad::PointerToLord, std::unique_ptr<Bmad::PointerToLord>>(
      m,
      "PointerToLord",
      "pointer_to_lord return type"
  )
      .def_readonly("control", &Bmad::PointerToLord::control)
      .def_readonly("ix_slave_back", &Bmad::PointerToLord::ix_slave_back)
      .def_readonly("ix_control", &Bmad::PointerToLord::ix_control)
      .def_readonly("ix_ic", &Bmad::PointerToLord::ix_ic)
      .def_readonly("lord_ptr", &Bmad::PointerToLord::lord_ptr)
      .def("__len__", [](const Bmad::PointerToLord &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToLord &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.control);
        if (i == 1)
          return py::cast(s.ix_slave_back);
        if (i == 2)
          return py::cast(s.ix_control);
        if (i == 3)
          return py::cast(s.ix_ic);
        if (i == 4)
          return py::cast(s.lord_ptr);
        throw py::index_error();
      });
  m.def(
      "pointer_to_lord",
      &Bmad::pointer_to_lord,
      py::arg("slave"),
      py::arg("ix_lord"),
      py::arg("lord_type") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_lord

Parameters
----------
slave : EleStruct
    Slave element.

ix_lord : int
    Index of the lord.

lord_type : int, optional
    See above.

Returns
-------
control : ControlStruct, optional
    Pointer to control info for this lord/slave relationship. Nullified if there is an error.

ix_slave_back : int, optional
    Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) will point back to slave. Set
    to -1 if there is an error or the slave is a slice_slave.

ix_control : int, optional
    Index in lat.control(:) array the control argument is at. For ramper lord elements, ix_control is index
    for the lord.control.ramper(:) array.

ix_ic : int, optional
    Index of the lat.ic(:) element associated with the control argument.

lord_ptr : EleStruct, optional
    Pointer to the lord. Nullified if there is an error.
)"""
  );
  py::class_<Bmad::PointerToMultipassLord, std::unique_ptr<Bmad::PointerToMultipassLord>>(
      m,
      "PointerToMultipassLord",
      "pointer_to_multipass_lord return type"
  )
      .def_readonly("ix_pass", &Bmad::PointerToMultipassLord::ix_pass)
      .def_readonly("super_lord", &Bmad::PointerToMultipassLord::super_lord)
      .def_readonly("multi_lord", &Bmad::PointerToMultipassLord::multi_lord)
      .def("__len__", [](const Bmad::PointerToMultipassLord &) { return 3; })
      .def("__getitem__", [](const Bmad::PointerToMultipassLord &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ix_pass);
        if (i == 1)
          return py::cast(s.super_lord);
        if (i == 2)
          return py::cast(s.multi_lord);
        throw py::index_error();
      });
  m.def(
      "pointer_to_multipass_lord",
      &Bmad::pointer_to_multipass_lord,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_multipass_lord

Parameters
----------
ele : EleStruct
    Lattice element.

Returns
-------
ix_pass : int, optional
    Multipass turn number. Set to 0 if element is a multipass_lord. Set to -1 if element is not a
    multipass_slave.

super_lord : EleStruct, optional
    super_lord of the element. Set to NULL if ele is not a super_slave or super_lord. Note: if ele is a
    multipass_lord there are multiple possible super_lord slaves.

multi_lord : EleStruct, optional
    multipass_lord if there is one. Set to NULL if there is no multipass_lord.
)"""
  );
  m.def(
      "pointer_to_next_ele",
      &Bmad::pointer_to_next_ele,
      py::arg("this_ele"),
      py::arg("next_ele"),
      py::arg("offset") = py::none(),
      py::arg("skip_beginning") = py::none(),
      py::arg("follow_fork") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_next_ele

Parameters
----------
this_ele : EleStruct
    Starting element.

next_ele : EleStruct
    Element after this_ele (if offset = 1). Nullified if there is an error. EG bad this_ele.

offset : int, optional
    +1 -> return next element, +2 -> element after that, etc. Can be negative. Default = +1.

skip_beginning : bool, optional
    If True then skip beginning element #0 when wrapping around. Default is False.

follow_fork : bool, optional
    If True then fork at any fork element. Default is False.
)"""
  );
  py::class_<Bmad::PointerToSlave, std::unique_ptr<Bmad::PointerToSlave>>(
      m,
      "PointerToSlave",
      "pointer_to_slave return type"
  )
      .def_readonly("control", &Bmad::PointerToSlave::control)
      .def_readonly("ix_lord_back", &Bmad::PointerToSlave::ix_lord_back)
      .def_readonly("ix_control", &Bmad::PointerToSlave::ix_control)
      .def_readonly("ix_ic", &Bmad::PointerToSlave::ix_ic)
      .def_readonly("slave_ptr", &Bmad::PointerToSlave::slave_ptr)
      .def("__len__", [](const Bmad::PointerToSlave &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSlave &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.control);
        if (i == 1)
          return py::cast(s.ix_lord_back);
        if (i == 2)
          return py::cast(s.ix_control);
        if (i == 3)
          return py::cast(s.ix_ic);
        if (i == 4)
          return py::cast(s.slave_ptr);
        throw py::index_error();
      });
  m.def(
      "pointer_to_slave",
      &Bmad::pointer_to_slave,
      py::arg("lord"),
      py::arg("ix_slave"),
      py::arg("lord_type") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_slave (lord, ix_slave, control, lord_type, ix_lord_back, ix_control, ix_ic) result (slave_ptr)

Function to point to a slave of a lord.
Note: Ramper lords do not have any associated slaves (slaves are assigned dynamically at run time).

If lord_type = all$ (the default) the range for ix_slave is:
  1 to lord%n_slave                                 for "regular" slaves.
  lord%n_slave+1 to lord%n_slave+lord%n_slave_field for field overlap slaves.

If lord_type = field_lord$, only the field overlap slaves may be accessed and the range for ix_slave is:
  1 to lord%n_slave_field

Also see:
  pointer_to_lord
  pointer_to_super_lord
  pointer_to_ele
  num_lords

Parameters
----------
lord : EleStruct
    Lord element

ix_slave : int
    Index of the slave in the list of slaves controled by the lord..

lord_type : int, optional
    See above.

Returns
-------
control : ControlStruct, optional
    Pointer to control info for this lord/slave relationship. Nullified if there is an error.

ix_lord_back : int, optional
    Index back to the lord. That is, pointer_to_lord(slave_ptr, ix_lord_back) will point back to the lord. Set
    to -1 if there is an error.

ix_control : int, optional
    Index in lat.control(:) array the control argument is at.

ix_ic : int, optional
    Index of the lat.ic(:) element associated with the control argument.

slave_ptr : EleStruct, optional
    Pointer to the slave. Nullified if there is an error.
)"""
  );
  py::class_<Bmad::PointerToSuperLord, std::unique_ptr<Bmad::PointerToSuperLord>>(
      m,
      "PointerToSuperLord",
      "pointer_to_super_lord return type"
  )
      .def_readonly("control", &Bmad::PointerToSuperLord::control)
      .def_readonly("ix_slave_back", &Bmad::PointerToSuperLord::ix_slave_back)
      .def_readonly("ix_control", &Bmad::PointerToSuperLord::ix_control)
      .def_readonly("ix_ic", &Bmad::PointerToSuperLord::ix_ic)
      .def_readonly("lord_ptr", &Bmad::PointerToSuperLord::lord_ptr)
      .def("__len__", [](const Bmad::PointerToSuperLord &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSuperLord &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.control);
        if (i == 1)
          return py::cast(s.ix_slave_back);
        if (i == 2)
          return py::cast(s.ix_control);
        if (i == 3)
          return py::cast(s.ix_ic);
        if (i == 4)
          return py::cast(s.lord_ptr);
        throw py::index_error();
      });
  m.def(
      "pointer_to_super_lord",
      &Bmad::pointer_to_super_lord,
      py::arg("slave"),
      py::arg("lord_type") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_super_lord

Parameters
----------
slave : EleStruct
    Slave element.

lord_type : int, optional
    If present, only return a super_lord of this type.

Returns
-------
control : ControlStruct, optional
    Pointer to control info for this lord/slave relationship. Nullified if there is an error.

ix_slave_back : int, optional
    Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) will point back to slave. Set
    to -1 if there is an error or the slave is a slice_slave.

ix_control : int, optional
    Index in lat.control(:) array the control argument is at. For ramper lord elements, ix_control is index
    for the lord.control.ramper(:) array.

ix_ic : int, optional
    Index of the lat.ic(:) element associated with the control argument.

lord_ptr : EleStruct, optional
    Pointer to the lord.
)"""
  );
  py::class_<
      Bmad::PointerToSurfaceDisplacementPt,
      std::unique_ptr<Bmad::PointerToSurfaceDisplacementPt>>(
      m,
      "PointerToSurfaceDisplacementPt",
      "pointer_to_surface_displacement_pt return type"
  )
      .def_readonly("ix", &Bmad::PointerToSurfaceDisplacementPt::ix)
      .def_readonly("iy", &Bmad::PointerToSurfaceDisplacementPt::iy)
      .def_readonly("xx", &Bmad::PointerToSurfaceDisplacementPt::xx)
      .def_readonly("yy", &Bmad::PointerToSurfaceDisplacementPt::yy)
      .def_readonly("pt", &Bmad::PointerToSurfaceDisplacementPt::pt)
      .def("__len__", [](const Bmad::PointerToSurfaceDisplacementPt &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSurfaceDisplacementPt &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.ix);
        if (i == 1)
          return py::cast(s.iy);
        if (i == 2)
          return py::cast(s.xx);
        if (i == 3)
          return py::cast(s.yy);
        if (i == 4)
          return py::cast(s.pt);
        throw py::index_error();
      });
  m.def(
      "pointer_to_surface_displacement_pt",
      &Bmad::pointer_to_surface_displacement_pt,
      py::arg("ele"),
      py::arg("nearest"),
      py::arg("x"),
      py::arg("y"),
      py::arg("extend_grid") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_surface_displacement_pt (ele, nearest, x, y, ix, iy, extend_grid, xx, yy) result (pt)

Routine to point to the grid point struct associated with point (x,y).

Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.

Parameters
----------
ele : EleStruct
    Element containing the grid

nearest : bool
    If True, return pointer to nearest grid point. If False, return pointer to the grid point lower and left
    of (x,y).

x : float
    Photon position.

y : float
    Photon position.

extend_grid : bool, optional
    If (x,y) past grid pretend (x,y) is at grid boundary. Default is False.

Returns
-------
ix : int, optional
    Grid point index.

iy : int, optional
    Grid point index.

xx : float, optional
    Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
    the nearest grid boundary point.

yy : float, optional
    Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
    the nearest grid boundary point.

pt : SurfaceDisplacementPtStruct, optional
    Pointer to grid point. Will not be associated if (x,y) outside the grid.
)"""
  );
  py::class_<Bmad::PointerToSurfaceSegmentedPt, std::unique_ptr<Bmad::PointerToSurfaceSegmentedPt>>(
      m,
      "PointerToSurfaceSegmentedPt",
      "pointer_to_surface_segmented_pt return type"
  )
      .def_readonly("ix", &Bmad::PointerToSurfaceSegmentedPt::ix)
      .def_readonly("iy", &Bmad::PointerToSurfaceSegmentedPt::iy)
      .def_readonly("xx", &Bmad::PointerToSurfaceSegmentedPt::xx)
      .def_readonly("yy", &Bmad::PointerToSurfaceSegmentedPt::yy)
      .def_readonly("pt", &Bmad::PointerToSurfaceSegmentedPt::pt)
      .def("__len__", [](const Bmad::PointerToSurfaceSegmentedPt &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSurfaceSegmentedPt &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.ix);
        if (i == 1)
          return py::cast(s.iy);
        if (i == 2)
          return py::cast(s.xx);
        if (i == 3)
          return py::cast(s.yy);
        if (i == 4)
          return py::cast(s.pt);
        throw py::index_error();
      });
  m.def(
      "pointer_to_surface_segmented_pt",
      &Bmad::pointer_to_surface_segmented_pt,
      py::arg("ele"),
      py::arg("nearest"),
      py::arg("x"),
      py::arg("y"),
      py::arg("extend_grid") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_surface_segmented_pt (ele, nearest, x, y, ix, iy, extend_grid, xx, yy) result (pt)

Routine to point to the grid point struct associated with point (x,y).

Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.

Parameters
----------
ele : EleStruct
    Element containing the grid

nearest : bool
    If True, return pointer to nearest grid point. If False, return pointer to the grid point lower and left
    of (x,y).

x : float
    Photon position.

y : float
    Photon position.

extend_grid : bool, optional
    If (x,y) past grid pretend (x,y) is at grid boundary. Default is False.

Returns
-------
ix : int, optional
    Grid point index.

iy : int, optional
    Grid point index.

xx : float, optional
    Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
    the nearest grid boundary point.

yy : float, optional
    Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
    the nearest grid boundary point.

pt : SurfaceSegmentedPtStruct, optional
    Pointer to grid point. Will not be associated if (x,y) outside the grid.
)"""
  );
  py::class_<Bmad::PointerToWakeEle, std::unique_ptr<Bmad::PointerToWakeEle>>(
      m,
      "PointerToWakeEle",
      "pointer_to_wake_ele return type"
  )
      .def_readonly("delta_s", &Bmad::PointerToWakeEle::delta_s)
      .def_readonly("wake_ele", &Bmad::PointerToWakeEle::wake_ele)
      .def("__len__", [](const Bmad::PointerToWakeEle &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToWakeEle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.delta_s);
        if (i == 1)
          return py::cast(s.wake_ele);
        throw py::index_error();
      });
  m.def(
      "pointer_to_wake_ele",
      &Bmad::pointer_to_wake_ele,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine pointer_to_wake_ele

Parameters
----------
ele : EleStruct
    Lattice element.

Returns
-------
delta_s : float, optional
    distance of wake locaiton from beginning of ele.

wake_ele : EleStruct, optional
    Element having the associated wake. wake_ele will be nullified if there is no associated wake.
)"""
  );
  py::class_<Bmad::PointerToWall3d, std::unique_ptr<Bmad::PointerToWall3d>>(
      m,
      "PointerToWall3d",
      "pointer_to_wall3d return type"
  )
      .def_readonly("ds_offset", &Bmad::PointerToWall3d::ds_offset)
      .def_readonly("is_branch_wall", &Bmad::PointerToWall3d::is_branch_wall)
      .def_readonly("wall3d", &Bmad::PointerToWall3d::wall3d)
      .def("__len__", [](const Bmad::PointerToWall3d &) { return 3; })
      .def("__getitem__", [](const Bmad::PointerToWall3d &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ds_offset);
        if (i == 1)
          return py::cast(s.is_branch_wall);
        if (i == 2)
          return py::cast(s.wall3d);
        throw py::index_error();
      });
  m.def(
      "pointer_to_wall3d",
      &Bmad::pointer_to_wall3d,
      py::arg("ele"),
      py::arg("ix_wall") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pointer_to_wall3d (ele, ix_wall, ds_offset, is_branch_wall) result (wall3d)

Function to return a pointer to a wall3d structure associated
with a given lattice element.

Note: The wall associated with a the vacuum chamber is the branch%wall3d.

Parameters
----------
ele : EleStruct
    lattice element.

ix_wall : int, optional
    index in wall3d(:) array. Default is 1.

Returns
-------
ds_offset : float, optional
    Element offset: s(beginning of ele) - s(beginning of wall3d)

is_branch_wall : bool, optional
    Set True if wall3d points to branch.wall3d.

wall3d : Wall3dStruct, optional
    Pointer to the associated wall structure. Will be nullified if there is no associated wall.
)"""
  );
  m.def(
      "polar_to_spinor",
      &Bmad::polar_to_spinor,
      py::arg("polar"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine polar_to_spinor

Parameters
----------
polar : SpinPolarStruct
    includes polar phase

Returns
-------
spinor : 1D array of complex (shape: 2)
    Spinor
)"""
  );
  m.def(
      "polar_to_vec",
      &Bmad::polar_to_vec,
      py::arg("polar"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine polar_to_vec

Parameters
----------
polar : SpinPolarStruct
    Spin_polar_struct

Returns
-------
vec : 1D array of float (shape: 3)
    Real(3)
)"""
  );
  py::class_<Bmad::ProjectEmitToXyz, std::unique_ptr<Bmad::ProjectEmitToXyz>>(
      m,
      "ProjectEmitToXyz",
      "project_emit_to_xyz return type"
  )
      .def_readonly("sigma_x", &Bmad::ProjectEmitToXyz::sigma_x)
      .def_readonly("sigma_y", &Bmad::ProjectEmitToXyz::sigma_y)
      .def_readonly("sigma_z", &Bmad::ProjectEmitToXyz::sigma_z)
      .def("__len__", [](const Bmad::ProjectEmitToXyz &) { return 3; })
      .def("__getitem__", [](const Bmad::ProjectEmitToXyz &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.sigma_x);
        if (i == 1)
          return py::cast(s.sigma_y);
        if (i == 2)
          return py::cast(s.sigma_z);
        throw py::index_error();
      });
  m.def(
      "project_emit_to_xyz",
      &Bmad::project_emit_to_xyz,
      py::arg("ring"),
      py::arg("ix"),
      py::arg("mode"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine project_emit_to_xyz(ring, ix, mode, sigma_x, sigma_y, sigma_z)

Obtains the projected x, y, and z beamsizes by building the sigma matrix
from the normal mode emittances and 1-turn transfer matrix.
These projectes beamsize are what would be seen by instrumentation.

This method of projecting takes into account transverse and longitudinal coupling.

This method of obtaining the projected beam sizes is from "Alternitive approach to general
coupled linear optics" by Andrzej Wolski.

The normal mode emittances used to generate a beam envelop sigma matrix from the
1-turn transfer matrix.  The projected sizes are from the 1, 1 3, 3 and 5, 5 elements of
the sigma matrix.

Parameters
----------
ring : LatStruct
    the storage ring

ix : int
    element at which to make the projection

mode : NormalModesStruct
    normal mode emittances

Returns
-------
sigma_x : float
    projected horizontal beamsize

sigma_y : float
    projected vertical beamsize

sigma_z : float
    projected longitudinal beamsize
)"""
  );
  m.def(
      "psi_prime_sca",
      &Bmad::psi_prime_sca,
      py::arg("t"),
      py::arg("p"),
      py::arg("args"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine psi_prime_sca(t, p, dpdt, args)

This wraps the array-valued psi_prime function as a scalar.

See psi_prime comments for details.

Parameters
----------
t : float
    time relative to RF bucket

p : float
    psi(t)

args : 1D array of float (shape: 1:8)
    parameters and constants of DEQ

Returns
-------
dpdt : float
    dpsi_dt
)"""
  );
  m.def(
      "ptc_bookkeeper",
      &Bmad::ptc_bookkeeper,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine ptc_bookkeeper

Parameters
----------
lat : LatStruct
    Bmad lattice.
)"""
  );
  m.def(
      "ptc_calculate_tracking_step_size",
      &Bmad::ptc_calculate_tracking_step_size,
      py::arg("ptc_layout"),
      py::arg("kl_max"),
      py::arg("ds_max") = py::none(),
      py::arg("even_steps") = py::none(),
      py::arg("r_typical") = py::none(),
      py::arg("dx_tol_bend") = py::none(),
      py::arg("use_2nd_order") = py::none(),
      py::arg("crossover") = py::none(),
      py::arg("crossover_wiggler") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_calculate_tracking_step_size (ptc_layout, kl_max, ds_max,
                                even_steps, r_typical, dx_tol_bend, use_2nd_order)

Routine to calculate the optimum number of tracking steps and order
of the integrator (2, 4, or 6) for each fibre in a layout.

See the Bmad manual chapter on PTC for more details.

Parameters
----------
ptc_layout : Layout
    This parameter is an input/output and is modified in-place.
    As an output, ptc_layout: Lattice with the optimum number of tracking steps and integrator order.

kl_max : float
    Maximum K1*L per tracking step.

ds_max : float, optional
    Maximum ds for any step. Useful when including other physicas like space charge.

even_steps : 1D array of bool (shape: 2), optional
    Always use an even number of steps for a fibre? Useful if need to evaluate at the center of fibres.

r_typical : float, optional
    Typical transverse offset. Used for computing the effective contribution of K1*L due to sextupoles.

dx_tol_bend : float, optional
    Tolerable residual orbit in a bend.

use_2nd_order : bool, optional
    If present and True then force the use of 2nd order integrator.

crossover : 1D array of int (shape: 2), optional
    crossover points between orders for all elements except wigglers. Default is [4, 18].

crossover_wiggler : 1D array of int (shape: 2), optional
    crossover points for wigglers. Default is [30, 60].
)"""
  );
  py::class_<Bmad::PtcCheckForLostParticle, std::unique_ptr<Bmad::PtcCheckForLostParticle>>(
      m,
      "PtcCheckForLostParticle",
      "ptc_check_for_lost_particle return type"
  )
      .def_readonly("state", &Bmad::PtcCheckForLostParticle::state)
      .def_readonly("ptc_fibre", &Bmad::PtcCheckForLostParticle::ptc_fibre)
      .def("__len__", [](const Bmad::PtcCheckForLostParticle &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcCheckForLostParticle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.state);
        if (i == 1)
          return py::cast(s.ptc_fibre);
        throw py::index_error();
      });
  m.def(
      "ptc_check_for_lost_particle",
      &Bmad::ptc_check_for_lost_particle,
      py::arg("do_reset"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_check_for_lost_particle (state, ptc_fibre, do_reset)

Routine to check if a particle has been lost when tracking with PTC.

Parameters
----------
do_reset : bool
    If True then reset ptc flags.

Returns
-------
state : int
    Same as coord_struct.state. alive$, lost$, lost_neg_x$, etc.

ptc_fibre : Fibre, optional
    Pointer to fibre where particle lost. Nullified if particle alive.
)"""
  );
  m.def(
      "ptc_closed_orbit_calc",
      &Bmad::ptc_closed_orbit_calc,
      py::arg("branch"),
      py::arg("radiation_damping_on") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_closed_orbit_calc (branch, closed_orbit, radiation_damping_on)

Routine to calculate the closed orbit of a lattice branch using PTC.
This routine assumes the associated PTC layout has been crated
with lat_to_ptc_layout.

Parameters
----------
branch : BranchStruct
    Branch of a lattice.

radiation_damping_on : bool, optional
    If True, radiation dampling is included in the calculation. Default is the setting of
    bmad_com..radiation_damping_on.

Returns
-------
closed_orbit : 1D array of CoordStruct
    closed_orbit
)"""
  );
  py::class_<Bmad::PtcEmitCalc, std::unique_ptr<Bmad::PtcEmitCalc>>(
      m,
      "PtcEmitCalc",
      "ptc_emit_calc return type"
  )
      .def_readonly("norm_mode", &Bmad::PtcEmitCalc::norm_mode)
      .def_readonly("closed_orb", &Bmad::PtcEmitCalc::closed_orb)
      .def("__len__", [](const Bmad::PtcEmitCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcEmitCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.norm_mode);
        if (i == 1)
          return py::cast(s.closed_orb);
        throw py::index_error();
      });
  m.def(
      "ptc_emit_calc",
      &Bmad::ptc_emit_calc,
      py::arg("ele"),
      py::arg("sigma_mat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)

Routine to calculate emittances, etc.

Note: This routine calls the PTC init_all routine.

Parameters
----------
ele : EleStruct
    Element at which to evaluate the parameters.

Returns
-------
norm_mode : NormalModesStruct
    Normal_modes_struct %a%tune, %b%tune, %z%tune %a%alpha_damp, etc. %a%emittance, etc.

closed_orb : CoordStruct
    Closed orbit at ele (Bmad coordinates). Notice: This closed orbit is the closed orbit with radiation on.
)"""
  );
  m.def(
      "ptc_layouts_resplit",
      &Bmad::ptc_layouts_resplit,
      py::arg("dKL_max"),
      py::arg("l_max"),
      py::arg("l_max_drift_only"),
      py::arg("bend_dorb"),
      py::arg("sex_dx"),
      py::arg("even") = py::none(),
      py::arg("crossover") = py::none(),
      py::arg("crossover_wiggler") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_layouts_resplit (dKL_max, l_max, l_max_drift_only, bend_dorb, sex_dx,
                                                          even, crossover, crossover_wiggler)

Routine to resplit (that is, recalculate the number of integration steps for an element)
For the fibres in all layouts. After doing a resplit, the tune (and any other relavent
"adjustable" parameters) should be adjusted to the correct values.

Parameters
----------
dKL_max : float
    Maximum K1 * L quadrupole strength allowed for an integration step. Reasonable value would be something
    like 0.04.

l_max : float
    Maximum step length. Ignored if set to 0.

l_max_drift_only : bool
    If True then l_max is only used for splitting drifts.

bend_dorb : float
    Residual bend orbit error. With some integration methods a zero orbit at the start of the bend will not be
    zero at the end. In this case, bend_dorb sets a maximum allowable orbit deviation. If set to zero, this
    argument will be ignored. A resonable value is 10d-7. Note that the actual orbit deviation is not simply
    related to bend_dorb and can be larger. In any case, lowering bend_dorb (without making it zero) will
    lower the

sex_dx : float
    To split sextupoles, sex_dx is used as the reference position about which the quadrupole strength is
    calculated. This quadrupole strength is then used with dKL_max to calculate the number of integration
    steps. Set to zero to ignore.

even : bool, optional
    If True then each fibre  will have an even number of steps. If False then the number of steps will be odd.
    If not present then number of steps is not constrained to be even or odd.

crossover : 1D array of int (shape: 2), optional
    crossover(1) sets the maximum number of 2nd order integration steps to use. If the number of steps would
    exceed crossover(1) then integration is switched to 4th order. crossover(2) sets the maximum number of 4th
    order integration steps. If this number is exceeded, 6th order integration is used. Currently the default
    in PTC is [4, 18].

crossover_wiggler : 1D array of int (shape: 2), optional
    crossover for wiggler elements.
)"""
  );
  m.def(
      "ptc_one_turn_mat_and_closed_orbit_calc",
      &Bmad::ptc_one_turn_mat_and_closed_orbit_calc,
      py::arg("branch"),
      py::arg("pz") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_one_turn_mat_and_closed_orbit_calc (branch, pz)

Routine to compute the transfer matrices for the individual elements and closed orbit
for a lattice branch with closed geometry.

Note: PTC itself does not compute Twiss parameters. Use twiss_from_mat6 to compute this.

Parameters
----------
branch : BranchStruct
    Lattice branch.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Lattice branch containing the matrices.

pz : float, optional
    energy offset around which to calculate the matrices if there is no RF.
)"""
  );
  m.def(
      "ptc_ran_seed_put",
      &Bmad::ptc_ran_seed_put,
      py::arg("iseed"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine ptc_ran_seed_put

Parameters
----------
iseed : int
    0 -> Use system clock.
)"""
  );
  m.def(
      "ptc_set_rf_state_for_c_normal",
      &Bmad::ptc_set_rf_state_for_c_normal,
      py::arg("nocavity"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine ptc_set_rf_state_for_c_normal

Parameters
----------
nocavity : bool
    True -> RF is off and vice versa.
)"""
  );
  m.def(
      "ptc_set_taylor_order_if_needed",
      &Bmad::ptc_set_taylor_order_if_needed,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_set_taylor_order_if_needed()

Routine to see if the taylor_order for PTC needs to be set/changed.
For example, for a change in bmad_com%taylor_order.
)"""
  );
  py::class_<Bmad::PtcSpinCalc, std::unique_ptr<Bmad::PtcSpinCalc>>(
      m,
      "PtcSpinCalc",
      "ptc_spin_calc return type"
  )
      .def_readonly("norm_mode", &Bmad::PtcSpinCalc::norm_mode)
      .def_readonly("closed_orb", &Bmad::PtcSpinCalc::closed_orb)
      .def("__len__", [](const Bmad::PtcSpinCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcSpinCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.norm_mode);
        if (i == 1)
          return py::cast(s.closed_orb);
        throw py::index_error();
      });
  m.def(
      "ptc_spin_calc",
      &Bmad::ptc_spin_calc,
      py::arg("ele"),
      py::arg("sigma_mat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_spin_calc (ele, norm_mode, sigma_mat, closed_orb)

Routine to equilibrium polarizations, etc.

Parameters
----------
ele : EleStruct
    Element at which to evaluate the parameters.

Returns
-------
norm_mode : NormalModesStruct
    Normal_modes_struct %a%tune, %b%tune, %z%tune %a%alpha_damp, etc. %a%emittance, etc.

closed_orb : CoordStruct
    Closed orbit at ele (Bmad coordinates). Notice: This closed orbit is the closed orbit with radiation on.
)"""
  );
  py::class_<Bmad::PtcTrackAll, std::unique_ptr<Bmad::PtcTrackAll>>(
      m,
      "PtcTrackAll",
      "ptc_track_all return type"
  )
      .def_readonly("track_state", &Bmad::PtcTrackAll::track_state)
      .def_readonly("err_flag", &Bmad::PtcTrackAll::err_flag)
      .def("__len__", [](const Bmad::PtcTrackAll &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcTrackAll &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.track_state);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "ptc_track_all",
      &Bmad::ptc_track_all,
      py::arg("branch"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine ptc_track_all (branch, orbit, track_state, err_flag)

Routine to track from the start to the end of a lattice branch.

Parameters
----------
branch : BranchStruct
    Lat to track through.

orbit : 1D array of CoordStruct
    Coordinates at beginning of branch.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Orbit array.

Returns
-------
track_state : int, optional
    Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.

err_flag : bool, optional
    Set true if particle lost or error. False otherwise
)"""
  );
  m.def(
      "ptc_transfer_map_with_spin",
      &Bmad::ptc_transfer_map_with_spin,
      py::arg("branch"),
      py::arg("t_map"),
      py::arg("s_map"),
      py::arg("orb0"),
      py::arg("ix1") = py::none(),
      py::arg("ix2") = py::none(),
      py::arg("one_turn") = py::none(),
      py::arg("unit_start") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine ptc_transfer_map_with_spin

Parameters
----------
branch : BranchStruct
    Lattice branch used in the calculation.

t_map : 1D array of TaylorStruct (shape: 6)
    Initial orbital map (used when unit_start = False)
    This parameter is an input/output and is modified in-place.
    As an output, t_map: Orbital transfer map.

s_map : 1D array of TaylorStruct (shape: 4)
    Initial spin map (used when unit_start = False)
    This parameter is an input/output and is modified in-place.
    As an output, s_map: Quaternion spin transfer map.

orb0 : CoordStruct
    Initial orbit around which the map is made.

ix1 : int, optional
    Element start index for the calculation. Default is 0.

ix2 : int, optional
    Element end index for the calculation. Default is branch.n_ele_track.

one_turn : bool, optional
    If present and True, and if ix1 = ix2, and the lattice is circular, then construct the one-turn map from
    ix1 back to ix1. Default = False.

unit_start : bool, optional
    If present and False then t_map will be used as the starting map instead of the unit map. Default = True

Returns
-------
err_flag : bool
    Set True if problem like number overflow, etc.
)"""
  );
  m.def(
      "pwd_mat",
      &Bmad::pwd_mat,
      py::arg("lat"),
      py::arg("t6"),
      py::arg("inductance"),
      py::arg("sig_z"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function pwd_mat(t6, inductance, sig_z) result (t6_pwd)

Calculates potential well distortion as RF defocusing.  Calculates t6_pwd=t6.Mpwd,
where Mpwd is identity with 65 element proportional to the inductance.

Vpwd = -inductance * lat%param%n_part * e_charge * c_light**3 / SQRT(twopi) / sig_z**3 / omega_RF  !effective RF voltage from PWD
Mpwd(6,5) = omega_RF * Vpwd / c_light / lat%ele(0)%value(E_TOT$) * branch%ele(i)%value(l$) / lat%param%total_length

Parameters
----------
lat : LatStruct
    TYPE(lat_struct)

t6 : 2D array of float (shape: 6,6)
    1-turn transfer matrix

inductance : float
    Longitudinal inductance in Henrys.  Something on the order of nH.

sig_z : float
    Bunch length.

Returns
-------
t6_pwd : 2D array of float (shape: 6,6)
    1-turn transfer matrix with PWD defocusing applied
)"""
  );
}
