#include "pybmad/generated/Bmad_routines_p.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_p(nb::module_ &m) {
  m.def(
      "p_func",
      &Bmad::p_func,
      nb::arg("E_in"),
      R"""(Wrapper for Fortran routine p_func

Parameters
----------
E_in : float

Returns
-------
rr1 : float
)"""
  );
  m.def(
      "parse_cartesian_map",
      &Bmad::parse_cartesian_map,
      nb::arg("ct_map"),
      nb::arg("ele"),
      nb::arg("lat"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
      R"""(Subroutine to parse a "cartesian_map = {}" construct

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
      nb::arg("cl_map"),
      nb::arg("ele"),
      nb::arg("lat"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
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
      "parse_gen_gradients",
      &Bmad::parse_gen_gradients,
      nb::arg("gg_map"),
      nb::arg("ele"),
      nb::arg("lat"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
      R"""(Subroutine to parse a "gen_gradients = {}" construct (curved-coordinate
generalized gradients). Each curve holds one GG derivative tower a_n, b_n, or
b_s selected by "kind = <a|b|bs>" and harmonic "n = <int>".
)"""
  );
  m.def(
      "parse_grid_field",
      &Bmad::parse_grid_field,
      nb::arg("g_field"),
      nb::arg("ele"),
      nb::arg("lat"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
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
      nb::arg("err_str"),
      nb::arg("lat"),
      nb::arg("int_array"),
      nb::arg("exact_size"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("open_delim") = nb::none(),
      nb::arg("separator") = nb::none(),
      nb::arg("close_delim") = nb::none(),
      nb::arg("default_value") = nb::none(),
      R"""(                                      separator, close_delim, default_value) result (is_ok)

Routine to parse a list of integers of the form:
   open_delim integer_1 separator integer_2 . . . close_delim
Example:   "(1.2, 2.3, 4.4, 8.5)"

Similar to parse_integer_list2 except does not use allocatable array.
See parse_integer_list2 for more details
)"""
  );
  nb::class_<Bmad::ParseIntegerList2>(m, "ParseIntegerList2", "parse_integer_list2 return type")
      .def_ro("num_found", &Bmad::ParseIntegerList2::num_found)
      .def_ro("delim", &Bmad::ParseIntegerList2::delim)
      .def_ro("delim_found", &Bmad::ParseIntegerList2::delim_found)
      .def_ro("is_ok", &Bmad::ParseIntegerList2::is_ok)
      .def("__len__", [](const Bmad::ParseIntegerList2 &) { return 4; })
      .def("__getitem__", [](const Bmad::ParseIntegerList2 &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.num_found);
        if (i == 1)
          return nb::cast(s.delim);
        if (i == 2)
          return nb::cast(s.delim_found);
        if (i == 3)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "parse_integer_list2",
      &Bmad::parse_integer_list2,
      nb::arg("err_str"),
      nb::arg("lat"),
      nb::arg("int_array"),
      nb::arg("num_expected") = nb::none(),
      nb::arg("open_delim") = nb::none(),
      nb::arg("separator") = nb::none(),
      nb::arg("close_delim") = nb::none(),
      nb::arg("default_value") = nb::none(),
      R"""(                                       open_delim, separator, close_delim, default_value) result (is_ok)

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
  m.def(
      "parse_line_or_list",
      &Bmad::parse_line_or_list,
      nb::arg("sequence"),
      nb::arg("iseq_tot"),
      nb::arg("lat"),
      nb::arg("top_level"),
      R"""(Subroutine to parse a sequence.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  nb::class_<Bmad::ParseRealList>(m, "ParseRealList", "parse_real_list return type")
      .def_ro("delim", &Bmad::ParseRealList::delim)
      .def_ro("delim_found", &Bmad::ParseRealList::delim_found)
      .def_ro("num_found", &Bmad::ParseRealList::num_found)
      .def_ro("is_ok", &Bmad::ParseRealList::is_ok)
      .def("__len__", [](const Bmad::ParseRealList &) { return 4; })
      .def("__getitem__", [](const Bmad::ParseRealList &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.delim);
        if (i == 1)
          return nb::cast(s.delim_found);
        if (i == 2)
          return nb::cast(s.num_found);
        if (i == 3)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "parse_real_list",
      &Bmad::parse_real_list,
      nb::arg("lat"),
      nb::arg("err_str"),
      nb::arg("real_array"),
      nb::arg("exact_size"),
      nb::arg("open_delim") = nb::none(),
      nb::arg("separator") = nb::none(),
      nb::arg("close_delim") = nb::none(),
      nb::arg("default_value") = nb::none(),
      R"""(                               separator, close_delim, default_value, num_found) result (is_ok)

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
  nb::class_<Bmad::ParseRealList2>(m, "ParseRealList2", "parse_real_list2 return type")
      .def_ro("num_found", &Bmad::ParseRealList2::num_found)
      .def_ro("delim", &Bmad::ParseRealList2::delim)
      .def_ro("delim_found", &Bmad::ParseRealList2::delim_found)
      .def_ro("is_ok", &Bmad::ParseRealList2::is_ok)
      .def("__len__", [](const Bmad::ParseRealList2 &) { return 4; })
      .def("__getitem__", [](const Bmad::ParseRealList2 &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.num_found);
        if (i == 1)
          return nb::cast(s.delim);
        if (i == 2)
          return nb::cast(s.delim_found);
        if (i == 3)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "parse_real_list2",
      &Bmad::parse_real_list2,
      nb::arg("lat"),
      nb::arg("err_str"),
      nb::arg("real_array"),
      nb::arg("num_expected") = nb::none(),
      nb::arg("open_brace") = nb::none(),
      nb::arg("separator") = nb::none(),
      nb::arg("close_brace") = nb::none(),
      nb::arg("default_value") = nb::none(),
      nb::arg("single_value") = nb::none(),
      R"""(                           open_delim, separator, close_delim, default_value, single_value) result (is_ok)

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
      "parse_superimpose_command",
      &Bmad::parse_superimpose_command,
      nb::arg("lat"),
      nb::arg("ele"),
      nb::arg("pele"),
      nb::arg("delim"),
      R"""(No docstring available.
)"""
  );
  m.def(
      "parser2_add_superimpose",
      [](LatStruct &lat, EleStruct &super_ele_in, ParserEleStruct &pele, LatStruct *in_lat) {
        auto fn = static_cast<
            void (*)(LatStruct &, EleStruct &, ParserEleStruct &, optional_ref<LatStruct>)>(
            &Bmad::parser2_add_superimpose
        );
        return fn(lat, super_ele_in, pele, ptr_to_opt_ref(in_lat));
      },
      nb::arg("lat"),
      nb::arg("super_ele_in"),
      nb::arg("pele"),
      nb::arg("in_lat") = nb::none(),
      R"""(Wrapper for Fortran routine parser2_add_superimpose

Parameters
----------
lat : LatStruct

super_ele_in : EleStruct

pele : ParserEleStruct

in_lat : LatStruct, optional
)"""
  );
  m.def(
      "parser_add_branch",
      &Bmad::parser_add_branch,
      nb::arg("fork_ele"),
      nb::arg("lat"),
      nb::arg("sequence"),
      nb::arg("seq_name"),
      nb::arg("seq_indexx"),
      nb::arg("no_end_marker"),
      nb::arg("in_lat"),
      nb::arg("plat"),
      nb::arg("created_new_branch"),
      nb::arg("new_branch_name") = nb::none(),
      R"""(                                 seq_name, seq_indexx, no_end_marker, in_lat, plat, created_new_branch)

Subroutine to do line expansion.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  m.def(
      "parser_add_constant",
      &Bmad::parser_add_constant,
      nb::arg("word"),
      nb::arg("lat"),
      nb::arg("redef_is_error"),
      R"""(Wrapper for Fortran routine parser_add_constant

Parameters
----------
word : str

lat : LatStruct

redef_is_error : bool
)"""
  );
  nb::class_<Bmad::ParserAddLords>(m, "ParserAddLords", "parser_add_lords return type")
      .def_ro("lat", &Bmad::ParserAddLords::lat)
      .def_ro("check_lat", &Bmad::ParserAddLords::check_lat)
      .def("__len__", [](const Bmad::ParserAddLords &) { return 2; })
      .def("__getitem__", [](const Bmad::ParserAddLords &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.lat);
        if (i == 1)
          return nb::cast(s.check_lat);
        throw nb::index_error();
      });
  m.def(
      "parser_add_lords",
      &Bmad::parser_add_lords,
      nb::arg("lord_lat"),
      nb::arg("n_ele_max"),
      nb::arg("plat"),
      R"""(Subroutine to add overlay, group, and girder lords to the lattice.
For overlays and groups: If multiple elements have the same name then
use all of them.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
lord_lat : LatStruct
    List of lord elements to add to lat.

n_ele_max : int
    lord elements in lord_lat are in range [1:n_ele_max].

plat : ParserLatStruct
    Extra info needed to place the lord elements

Returns
-------
lat : LatStruct
    Lattice to add lord elements to.

check_lat : LatStruct, optional
    If slave elements of a lord are not in lat but are in check_lat, do not issue error message about slave
    elements not found.
)"""
  );
  m.def(
      "parser_add_superimpose",
      &Bmad::parser_add_superimpose,
      nb::arg("branch"),
      nb::arg("super_ele_in"),
      nb::arg("pele"),
      nb::arg("in_lat"),
      nb::arg("plat"),
      R"""(Wrapper for Fortran routine parser_add_superimpose

Parameters
----------
branch : BranchStruct

super_ele_in : EleStruct

pele : ParserEleStruct

in_lat : LatStruct

plat : ParserLatStruct
)"""
  );
  m.def(
      "parser_call_check",
      &Bmad::parser_call_check,
      nb::arg("word"),
      nb::arg("ix_word"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("call_found"),
      nb::arg("err_flag") = nb::none(),
      R"""(Routine to check if there is a "call::XXX" construct in the input stream.
)"""
  );
  m.def(
      "parser_debug_print_info",
      &Bmad::parser_debug_print_info,
      nb::arg("lat"),
      nb::arg("debug_line"),
      nb::arg("sequence") = nb::none(),
      R"""(Subroutine to remove all null_ele elements.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  m.def(
      "parser_error",
      [](std::string what1,
         std::optional<std::string> what2,
         std::optional<std::string> what3,
         std::optional<std::string> what4,
         SeqStruct *seq,
         ParserEleStruct *pele,
         std::optional<bool> stop_here,
         std::optional<int> level,
         std::optional<FArray1D<Real>> r_array,
         std::optional<FArray1D<Int>> i_array) {
        auto fn = static_cast<void (*)(
            std::string,
            std::optional<std::string>,
            std::optional<std::string>,
            std::optional<std::string>,
            optional_ref<SeqStruct>,
            optional_ref<ParserEleStruct>,
            std::optional<bool>,
            std::optional<int>,
            std::optional<FArray1D<Real>>,
            std::optional<FArray1D<Int>>
        )>(&Bmad::parser_error);
        return fn(
            what1,
            what2,
            what3,
            what4,
            ptr_to_opt_ref(seq),
            ptr_to_opt_ref(pele),
            stop_here,
            level,
            r_array,
            i_array
        );
      },
      nb::arg("what1"),
      nb::arg("what2") = nb::none(),
      nb::arg("what3") = nb::none(),
      nb::arg("what4") = nb::none(),
      nb::arg("seq") = nb::none(),
      nb::arg("pele") = nb::none(),
      nb::arg("stop_here") = nb::none(),
      nb::arg("level") = nb::none(),
      nb::arg("r_array") = nb::none(),
      nb::arg("i_array") = nb::none(),
      R"""(Routine to print an error message generated when parsing a lattice.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Parameters
----------
what1 : str
    First line in error message.

what2 : str, optional
    Second line in error message.

what3 : str, optional
    Third line in error message.

what4 : str, optional
    Fourth line in error message.

seq : SeqStruct, optional
    Used when error is generated during reading of a lattice file. Contains information such as file name, and
    line number where error was detected.

pele : ParserEleStruct, optional
    Used when error is associated with a lattice element. Contains information on the lattice element.

stop_here : bool, optional
    If present and True then immediately stop. Used with severe errors.

level : int, optional
    Possibilities are:

r_array : 1D array of float, optional
    Real numbers to be encoded in error message. See out_io doc.

i_array : 1D array of int, optional
    Integer numbers to be encoded in error message. See out_io doc.
)"""
  );
  nb::class_<Bmad::ParserExpandLine>(m, "ParserExpandLine", "parser_expand_line return type")
      .def_ro("n_ele_expand", &Bmad::ParserExpandLine::n_ele_expand)
      .def_ro("expanded_line", &Bmad::ParserExpandLine::expanded_line)
      .def("__len__", [](const Bmad::ParserExpandLine &) { return 2; })
      .def("__getitem__", [](const Bmad::ParserExpandLine &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.n_ele_expand);
        if (i == 1)
          return nb::cast(s.expanded_line);
        throw nb::index_error();
      });
  m.def(
      "parser_expand_line",
      [](int i_lev,
         std::string line_name,
         SeqStructAlloc1D sequence,
         CharacterAlloc1D &seq_name,
         IntAlloc1D &seq_indexx,
         bool no_end_marker,
         LatStruct *lat,
         LatStruct *in_lat) {
        auto fn = static_cast<Bmad::ParserExpandLine (*)(
            int,
            std::string,
            SeqStructAlloc1D,
            CharacterAlloc1D &,
            IntAlloc1D &,
            bool,
            optional_ref<LatStruct>,
            optional_ref<LatStruct>
        )>(&Bmad::parser_expand_line);
        return fn(
            i_lev,
            line_name,
            sequence,
            seq_name,
            seq_indexx,
            no_end_marker,
            ptr_to_opt_ref(lat),
            ptr_to_opt_ref(in_lat)
        );
      },
      nb::arg("i_lev"),
      nb::arg("line_name"),
      nb::arg("sequence"),
      nb::arg("seq_name"),
      nb::arg("seq_indexx"),
      nb::arg("no_end_marker"),
      nb::arg("lat") = nb::none(),
      nb::arg("in_lat") = nb::none(),
      R"""(              seq_name, seq_indexx, no_end_marker, n_ele_expand, lat, in_lat, expanded_line)

Subroutine to do line expansion.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

Note: Either lat and in_lat must be present or expanded_line must be present.

Parameters
----------
i_lev : int
    Subsequence level. 1 => Root level.

line_name : str
    Root line to expand.

sequence : 1D array of SeqStruct
    Array of sequencies.

seq_name : 1D array of str
    Array of sequence names.

seq_indexx : 1D array of int
    Index array for the sequence names.

no_end_marker : bool
    Put a marker named "end" at the end of the branch?

lat : LatStruct, optional
    Lattice to put the expanded line
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with new line. Except if expanded_line is present.

in_lat : LatStruct, optional
    Lattice with array of defined elements.

Returns
-------
n_ele_expand : int
    Number of elements in the finished line.

expanded_line : 1D array of BaseLineEleStruct, optional
    If present, lat argument will be ignored and the expanded line will be put into expanded_line.
)"""
  );
  nb::class_<Bmad::ParserFastComplexRead>(
      m,
      "ParserFastComplexRead",
      "parser_fast_complex_read return type"
  )
      .def_ro("delim", &Bmad::ParserFastComplexRead::delim)
      .def_ro("is_ok", &Bmad::ParserFastComplexRead::is_ok)
      .def("__len__", [](const Bmad::ParserFastComplexRead &) { return 2; })
      .def("__getitem__", [](const Bmad::ParserFastComplexRead &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.delim);
        if (i == 1)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "parser_fast_complex_read",
      &Bmad::parser_fast_complex_read,
      nb::arg("cmplx_vec"),
      nb::arg("ele"),
      nb::arg("err_str"),
      R"""(Routine to read an array of complex numbers.

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
      nb::arg("int_vec"),
      nb::arg("ele"),
      nb::arg("delim_wanted"),
      nb::arg("err_str"),
      R"""(No docstring available.
)"""
  );
  nb::class_<Bmad::ParserFastRealRead>(m, "ParserFastRealRead", "parser_fast_real_read return type")
      .def_ro("delim", &Bmad::ParserFastRealRead::delim)
      .def_ro("n_real", &Bmad::ParserFastRealRead::n_real)
      .def_ro("is_ok", &Bmad::ParserFastRealRead::is_ok)
      .def("__len__", [](const Bmad::ParserFastRealRead &) { return 3; })
      .def("__getitem__", [](const Bmad::ParserFastRealRead &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.delim);
        if (i == 1)
          return nb::cast(s.n_real);
        if (i == 2)
          return nb::cast(s.is_ok);
        throw nb::index_error();
      });
  m.def(
      "parser_fast_real_read",
      &Bmad::parser_fast_real_read,
      nb::arg("real_vec"),
      nb::arg("ele"),
      nb::arg("end_delims"),
      nb::arg("err_str"),
      nb::arg("exact_size") = nb::none(),
      R"""(Routine to read an array of real numbers.

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
      nb::arg("how"),
      nb::arg("file_name_in") = nb::none(),
      nb::arg("finished") = nb::none(),
      nb::arg("err") = nb::none(),
      nb::arg("open_file") = nb::none(),
      nb::arg("abort_on_open_error") = nb::none(),
      R"""(Subroutine to keep track of the files that are opened for reading.
This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  m.def(
      "parser_get_integer",
      &Bmad::parser_get_integer,
      nb::arg("int_val"),
      nb::arg("word"),
      nb::arg("ix_word"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err"),
      nb::arg("str1") = nb::none(),
      nb::arg("str2") = nb::none(),
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
      nb::arg("attrib_name"),
      nb::arg("this_logic"),
      nb::arg("ele_name"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err"),
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
      nb::arg("lat"),
      R"""(Routine to identify the elements the forks in a lattice are branching to.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.
)"""
  );
  m.def(
      "parser_init_custom_elements",
      &Bmad::parser_init_custom_elements,
      nb::arg("lat"),
      R"""(No docstring available.
)"""
  );
  m.def(
      "parser_print_line",
      &Bmad::parser_print_line,
      nb::arg("lat"),
      nb::arg("end_of_file"),
      R"""(This routine is called when a print statement is found in the lattice file.
)"""
  );
  m.def(
      "parser_read_lr_wake",
      &Bmad::parser_read_lr_wake,
      nb::arg("ele"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
      R"""(Subroutine to read in a long-range wake field from an external file.
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
      nb::arg("ele"),
      nb::arg("lr_file_name"),
      R"""(Subroutine to read in a long-range wake field from an external file.
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
      nb::arg("ele"),
      nb::arg("sr_file_name"),
      R"""(Subroutine to read in a short-range wake field from an external file.
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
      nb::arg("ele"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("err_flag"),
      R"""(Subroutine to read in a short-range wake field.
This subroutine is used by bmad_parser and bmad_parser2.

Parameters
----------
ele : EleStruct
    Element containing wake structure.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with wake information.
)"""
  );
  nb::class_<Bmad::ParserSetAttribute>(m, "ParserSetAttribute", "parser_set_attribute return type")
      .def_ro("delim", &Bmad::ParserSetAttribute::delim)
      .def_ro("delim_found", &Bmad::ParserSetAttribute::delim_found)
      .def_ro("err_flag", &Bmad::ParserSetAttribute::err_flag)
      .def_ro("pele", &Bmad::ParserSetAttribute::pele)
      .def("__len__", [](const Bmad::ParserSetAttribute &) { return 4; })
      .def("__getitem__", [](const Bmad::ParserSetAttribute &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.delim);
        if (i == 1)
          return nb::cast(s.delim_found);
        if (i == 2)
          return nb::cast(s.err_flag);
        if (i == 3)
          return nb::cast(s.pele);
        throw nb::index_error();
      });
  m.def(
      "parser_set_attribute",
      &Bmad::parser_set_attribute,
      nb::arg("how"),
      nb::arg("ele"),
      nb::arg("check_free") = nb::none(),
      nb::arg("heterogeneous_ele_list") = nb::none(),
      nb::arg("set_field_master") = nb::none(),
      R"""(                                                                heterogeneous_ele_list, set_field_master)

Subroutine used by bmad_parser and bmad_parser2 to get the value of
an attribute from the input file and set the appropriate value in an element.

This subroutine is not intended for general use.

Parameters
----------
how : int
    Either def$ if the element is being construct from scratch or redef$ if the element has already been
    formed and this is part of a "ele_name[attrib_name] = value" construct.

ele : EleStruct
    Element whose attribute this is.

check_free : bool, optional
    If present and True then an error will be generated if the attribute is not free to vary. Used by
    bmad_parser2.

heterogeneous_ele_list : bool, optional
    If True (default = False), we are parsing something like something like "*[tracking_method] =
    runge_kutta". In this case, runge_kutta may not be valid for this ele but this is not an error.

set_field_master : bool, optional
    If True (the default) set ele.field_master = T if the attribute to be set is something like DB_FIELD. If
    this routine is being called post lattice parsing, setting ele.field_master is *not* wanted.

Returns
-------
delim : str
    Delimiter found where the parsing of the input line stops.

delim_found : bool
    Delimiter found? False if end of input command.

err_flag : bool
    Set True if there is a problem parsing the input.

pele : ParserEleStruct, optional
    Structure to hold additional information that cannot be stored in the ele argument.
)"""
  );
  m.def(
      "parser_transfer_control_struct",
      &Bmad::parser_transfer_control_struct,
      nb::arg("con_in"),
      nb::arg("lord"),
      nb::arg("ix_var"),
      R"""(Routine to transfer the information from an input control_struct (which stores
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
      nb::arg("orb"),
      nb::arg("branch"),
      nb::arg("in_time_coordinates") = nb::none(),
      nb::arg("in_body_frame") = nb::none(),
      nb::arg("w_mat_out") = nb::none(),
      R"""(Returns the particle in global time coordinates given is coordinates orb in lattice lat.

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
      nb::arg("orbit"),
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
      nb::arg("orbit"),
      nb::arg("dir") = nb::none(),
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
      nb::arg("orbit"),
      nb::arg("ele"),
      nb::arg("reference_active_edge") = nb::none(),
      nb::arg("s_rel") = nb::none(),
      nb::arg("time_coords") = nb::none(),
      nb::arg("rf_freq") = nb::none(),
      nb::arg("abs_time") = nb::none(),
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
      nb::arg("x_pitch"),
      nb::arg("y_pitch"),
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
      nb::arg("patch"),
      nb::arg("ref_coords") = nb::none(),
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
  nb::class_<Bmad::PhotonAbsorptionAndPhaseShift>(
      m,
      "PhotonAbsorptionAndPhaseShift",
      "photon_absorption_and_phase_shift return type"
  )
      .def_ro("absorption", &Bmad::PhotonAbsorptionAndPhaseShift::absorption)
      .def_ro("phase_shift", &Bmad::PhotonAbsorptionAndPhaseShift::phase_shift)
      .def_ro("err_flag", &Bmad::PhotonAbsorptionAndPhaseShift::err_flag)
      .def("__len__", [](const Bmad::PhotonAbsorptionAndPhaseShift &) { return 3; })
      .def("__getitem__", [](const Bmad::PhotonAbsorptionAndPhaseShift &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.absorption);
        if (i == 1)
          return nb::cast(s.phase_shift);
        if (i == 2)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "photon_absorption_and_phase_shift",
      &Bmad::photon_absorption_and_phase_shift,
      nb::arg("material"),
      nb::arg("Energy"),
      R"""(Routine to calcualte the absorption and phase shift values for a photon with a given
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
  nb::class_<Bmad::PhotonAddToDetectorStatistics>(
      m,
      "PhotonAddToDetectorStatistics",
      "photon_add_to_detector_statistics return type"
  )
      .def_ro("ix_pt", &Bmad::PhotonAddToDetectorStatistics::ix_pt)
      .def_ro("iy_pt", &Bmad::PhotonAddToDetectorStatistics::iy_pt)
      .def("__len__", [](const Bmad::PhotonAddToDetectorStatistics &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonAddToDetectorStatistics &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ix_pt);
        if (i == 1)
          return nb::cast(s.iy_pt);
        throw nb::index_error();
      });
  m.def(
      "photon_add_to_detector_statistics",
      [](CoordStruct &orbit0, CoordStruct &orbit, EleStruct &ele, PixelPtStruct *pixel_pt) {
        auto fn = static_cast<Bmad::PhotonAddToDetectorStatistics (*)(
            CoordStruct &,
            CoordStruct &,
            EleStruct &,
            optional_ref<PixelPtStruct>
        )>(&Bmad::photon_add_to_detector_statistics);
        return fn(orbit0, orbit, ele, ptr_to_opt_ref(pixel_pt));
      },
      nb::arg("orbit0"),
      nb::arg("orbit"),
      nb::arg("ele"),
      nb::arg("pixel_pt") = nb::none(),
      R"""(Routine to add photon statistics to the appropriate pixel of a "detector" grid.

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
  nb::class_<Bmad::PhotonDiffuseScattering>(
      m,
      "PhotonDiffuseScattering",
      "photon_diffuse_scattering return type"
  )
      .def_ro("graze_angle_out", &Bmad::PhotonDiffuseScattering::graze_angle_out)
      .def_ro("phi_out", &Bmad::PhotonDiffuseScattering::phi_out)
      .def_ro("diffuse_param", &Bmad::PhotonDiffuseScattering::diffuse_param)
      .def("__len__", [](const Bmad::PhotonDiffuseScattering &) { return 3; })
      .def("__getitem__", [](const Bmad::PhotonDiffuseScattering &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.graze_angle_out);
        if (i == 1)
          return nb::cast(s.phi_out);
        if (i == 2)
          return nb::cast(s.diffuse_param);
        throw nb::index_error();
      });
  m.def(
      "photon_diffuse_scattering",
      &Bmad::photon_diffuse_scattering,
      nb::arg("graze_angle_in"),
      nb::arg("energy"),
      nb::arg("surface"),
      R"""(Routine to simulate the diffuse scattering of photons. The outgoing angles are
choosen using the Dugan distribution.

Also see: photon_reflection.
Use photon_reflection_std_surface_init or read_surface_reflection_file to get surface info.

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

diffuse_param : DiffuseParamStruct, optional
    Internal parameters used in the calculation. This is used for diagnostics and is not used in standard
    simulations.
)"""
  );
  nb::class_<Bmad::PhotonHitFunc>(m, "PhotonHitFunc", "photon_hit_func return type")
      .def_ro("status", &Bmad::PhotonHitFunc::status)
      .def_ro("d_radius", &Bmad::PhotonHitFunc::d_radius)
      .def("__len__", [](const Bmad::PhotonHitFunc &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonHitFunc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.status);
        if (i == 1)
          return nb::cast(s.d_radius);
        throw nb::index_error();
      });
  m.def(
      "photon_hit_func",
      &Bmad::photon_hit_func,
      nb::arg("track_len"),
      R"""(Routine to be used as an argument in zbrent in the capillary_photon_hit_spot_calc.
Made a module procedure (not nested) to avoid a stack trampoline.

Parameters
----------
track_len : float
    Place to position the photon.

Returns
-------
status : int
    Not set.

d_radius : float
    r_photon - r_wall
)"""
  );
  m.def(
      "photon_read_spline",
      &Bmad::photon_read_spline,
      nb::arg("spline_dir"),
      R"""(Routine to initialize a photon using a set of spline fits.

Parameters
----------
spline_dir : str
    Root directory for the spline fits.

Returns
-------
splines : PhotonInitSplinesStruct
    Spline structure
)"""
  );
  nb::class_<Bmad::PhotonReflection>(m, "PhotonReflection", "photon_reflection return type")
      .def_ro("graze_angle_out", &Bmad::PhotonReflection::graze_angle_out)
      .def_ro("phi_out", &Bmad::PhotonReflection::phi_out)
      .def("__len__", [](const Bmad::PhotonReflection &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonReflection &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.graze_angle_out);
        if (i == 1)
          return nb::cast(s.phi_out);
        throw nb::index_error();
      });
  m.def(
      "photon_reflection",
      &Bmad::photon_reflection,
      nb::arg("graze_angle_in"),
      nb::arg("energy"),
      nb::arg("surface"),
      R"""(Routine to reflect a photon from a surface including both diffuse and specular reflections.

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
      R"""($OMP THREADPRIVATE(dr_d_param_ptr, dr_surface_ptr, dr_old_integral, dr_tot_integral, &
$OMP                dr_ran1, dr_ran2, dr_j)

 Subroutine photon_reflection_std_surface_init (surface)

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
  nb::class_<Bmad::PhotonReflectivity>(m, "PhotonReflectivity", "photon_reflectivity return type")
      .def_ro("p_reflect", &Bmad::PhotonReflectivity::p_reflect)
      .def_ro("rel_p_specular", &Bmad::PhotonReflectivity::rel_p_specular)
      .def("__len__", [](const Bmad::PhotonReflectivity &) { return 2; })
      .def("__getitem__", [](const Bmad::PhotonReflectivity &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.p_reflect);
        if (i == 1)
          return nb::cast(s.rel_p_specular);
        throw nb::index_error();
      });
  m.def(
      "photon_reflectivity",
      &Bmad::photon_reflectivity,
      nb::arg("angle"),
      nb::arg("energy"),
      nb::arg("surface"),
      R"""(Routine to evaluate the photon reflectivity.
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
      nb::arg("aperture_ele"),
      nb::arg("x_lim"),
      nb::arg("y_lim"),
      nb::arg("z_lim"),
      nb::arg("source_ele"),
      R"""(Routine to calculate the corner coords in the source_ele ref frame.

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
      nb::arg("ele"),
      R"""(Routine to calculate and store the parmeters needed for photon targeting.
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
      nb::arg("ele"),
      R"""(Routine to return the type of photon to be tracked: coherent$ or incoherent$.

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
      nb::arg("track_end"),
      nb::arg("orbit"),
      nb::arg("ele_orientation"),
      nb::arg("return_stream_end") = nb::none(),
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
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("orbit"),
      nb::arg("direction"),
      nb::arg("max_target_area"),
      nb::arg("w_to_surface") = nb::none(),
      R"""(Routine to emit a photon from a point that may be on a surface.
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
  nb::class_<Bmad::PointerToAttribute>(m, "PointerToAttribute", "pointer_to_attribute return type")
      .def_ro("a_ptr", &Bmad::PointerToAttribute::a_ptr)
      .def_ro("err_flag", &Bmad::PointerToAttribute::err_flag)
      .def_ro("ix_attrib", &Bmad::PointerToAttribute::ix_attrib)
      .def("__len__", [](const Bmad::PointerToAttribute &) { return 3; })
      .def("__getitem__", [](const Bmad::PointerToAttribute &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.a_ptr);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.ix_attrib);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_attribute",
      &Bmad::pointer_to_attribute,
      nb::arg("ele"),
      nb::arg("attrib_name"),
      nb::arg("do_allocation"),
      nb::arg("err_print_flag") = nb::none(),
      nb::arg("do_unlink") = nb::none(),
      R"""(Wrapper for Fortran routine pointer_to_attribute

Parameters
----------
ele : EleStruct
    After this routine finishes Ptr_attrib will point to a variable within this element.

attrib_name : str
    Name of attribute. Must be uppercase. For example: "HKICK".

do_allocation : bool
    If True then do an allocation if needed. EG: The multipole An and Bn arrays need to be allocated before
    their use.

err_print_flag : bool, optional
    If present and False then suppress printing of an error message on error.

do_unlink : bool, optional
    Default False. If True and applicable, unlink the structure containing the attribute. See above for
    details.

Returns
-------
a_ptr : AllPointerStruct
    Pointer to the attribute.

err_flag : bool
    Set True if attribtute not found. False otherwise.

ix_attrib : int, optional
    If applicable, this is the index to the attribute in the ele.value(:), ele.control.var(:), ele.a_pole(:)
    or ele.b_pole(:) arrays. Set to 0 if not in any of these arrays.
)"""
  );
  m.def(
      "pointer_to_branch",
      nb::overload_cast<EleStruct &>(&Bmad::pointer_to_branch),
      nb::arg("ele"),
      R"""(Routine to return a pointer to the lattice branch associated with a given name
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
      nb::overload_cast<std::string, LatStruct &, std::optional<bool>, std::optional<int>>(
          &Bmad::pointer_to_branch
      ),
      nb::arg("branch_name"),
      nb::arg("lat"),
      nb::arg("parameter_is_branch0") = nb::none(),
      nb::arg("blank_branch") = nb::none(),
      R"""(Routine to return a pointer to the lattice branch associated with a given name
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
      nb::overload_cast<LatStruct &, int, std::optional<int>>(&Bmad::pointer_to_ele),
      nb::arg("lat"),
      nb::arg("ix_ele"),
      nb::arg("ix_branch") = nb::none(),
      R"""(Routine to return a pointer to an element.
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
      nb::overload_cast<LatStruct &, LatEleLocStruct &>(&Bmad::pointer_to_ele),
      nb::arg("lat"),
      nb::arg("ele_loc"),
      R"""(Routine to return a pointer to an element.
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
      nb::overload_cast<LatStruct &, std::string>(&Bmad::pointer_to_ele),
      nb::arg("lat"),
      nb::arg("ele_name"),
      R"""(Routine to return a pointer to an element.
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
      nb::overload_cast<LatStruct &, EleStruct &>(&Bmad::pointer_to_ele),
      nb::arg("lat"),
      nb::arg("foreign_ele"),
      R"""(Routine to return a pointer to an element.
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
  nb::class_<Bmad::PointerToElementAtS>(
      m,
      "PointerToElementAtS",
      "pointer_to_element_at_s return type"
  )
      .def_ro("err_flag", &Bmad::PointerToElementAtS::err_flag)
      .def_ro("s_eff", &Bmad::PointerToElementAtS::s_eff)
      .def_ro("position", &Bmad::PointerToElementAtS::position)
      .def_ro("ele", &Bmad::PointerToElementAtS::ele)
      .def("__len__", [](const Bmad::PointerToElementAtS &) { return 4; })
      .def("__getitem__", [](const Bmad::PointerToElementAtS &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.s_eff);
        if (i == 2)
          return nb::cast(s.position);
        if (i == 3)
          return nb::cast(s.ele);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_element_at_s",
      &Bmad::pointer_to_element_at_s,
      nb::arg("branch"),
      nb::arg("s"),
      nb::arg("choose_max"),
      nb::arg("print_err") = nb::none(),
      R"""(Function to return a pointer to the element at position s.
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
      nb::arg("ele"),
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
  nb::class_<Bmad::PointerToFieldEle>(m, "PointerToFieldEle", "pointer_to_field_ele return type")
      .def_ro("dz_offset", &Bmad::PointerToFieldEle::dz_offset)
      .def_ro("field_ele", &Bmad::PointerToFieldEle::field_ele)
      .def("__len__", [](const Bmad::PointerToFieldEle &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToFieldEle &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.dz_offset);
        if (i == 1)
          return nb::cast(s.field_ele);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_field_ele",
      &Bmad::pointer_to_field_ele,
      nb::arg("ele"),
      nb::arg("ix_field_ele"),
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
  nb::class_<Bmad::PointerToGirder>(m, "PointerToGirder", "pointer_to_girder return type")
      .def_ro("ix_slave_back", &Bmad::PointerToGirder::ix_slave_back)
      .def_ro("girder", &Bmad::PointerToGirder::girder)
      .def("__len__", [](const Bmad::PointerToGirder &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToGirder &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ix_slave_back);
        if (i == 1)
          return nb::cast(s.girder);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_girder",
      &Bmad::pointer_to_girder,
      nb::arg("ele"),
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
  nb::class_<Bmad::PointerToIndexedAttribute>(
      m,
      "PointerToIndexedAttribute",
      "pointer_to_indexed_attribute return type"
  )
      .def_ro("a_ptr", &Bmad::PointerToIndexedAttribute::a_ptr)
      .def_ro("err_flag", &Bmad::PointerToIndexedAttribute::err_flag)
      .def("__len__", [](const Bmad::PointerToIndexedAttribute &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToIndexedAttribute &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.a_ptr);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_indexed_attribute",
      &Bmad::pointer_to_indexed_attribute,
      nb::arg("ele"),
      nb::arg("ix_attrib"),
      nb::arg("do_allocation"),
      nb::arg("err_print_flag") = nb::none(),
      R"""(Wrapper for Fortran routine pointer_to_indexed_attribute

Parameters
----------
ele : EleStruct
    After this routine finishes A_ptr will point to a variable within this element.

ix_attrib : int
    Integer, Attribute index.

do_allocation : bool
    If True then do an allocation if needed. EG: The multipole An and Bn arrays need to be allocated before
    their use.

err_print_flag : bool, optional
    If present and False then suppress printing of an error message on error.

Returns
-------
a_ptr : AllPointerStruct
    Pointer to the attribute.

err_flag : bool
    Set True if attribtute not found. False otherwise.
)"""
  );
  nb::class_<Bmad::PointerToLord>(m, "PointerToLord", "pointer_to_lord return type")
      .def_ro("control", &Bmad::PointerToLord::control)
      .def_ro("ix_slave_back", &Bmad::PointerToLord::ix_slave_back)
      .def_ro("ix_control", &Bmad::PointerToLord::ix_control)
      .def_ro("ix_ic", &Bmad::PointerToLord::ix_ic)
      .def_ro("lord_ptr", &Bmad::PointerToLord::lord_ptr)
      .def("__len__", [](const Bmad::PointerToLord &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToLord &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.control);
        if (i == 1)
          return nb::cast(s.ix_slave_back);
        if (i == 2)
          return nb::cast(s.ix_control);
        if (i == 3)
          return nb::cast(s.ix_ic);
        if (i == 4)
          return nb::cast(s.lord_ptr);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_lord",
      &Bmad::pointer_to_lord,
      nb::arg("slave"),
      nb::arg("ix_lord"),
      nb::arg("lord_type") = nb::none(),
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
  nb::class_<Bmad::PointerToMultipassLord>(
      m,
      "PointerToMultipassLord",
      "pointer_to_multipass_lord return type"
  )
      .def_ro("ix_pass", &Bmad::PointerToMultipassLord::ix_pass)
      .def_ro("super_lord", &Bmad::PointerToMultipassLord::super_lord)
      .def_ro("multi_lord", &Bmad::PointerToMultipassLord::multi_lord)
      .def("__len__", [](const Bmad::PointerToMultipassLord &) { return 3; })
      .def("__getitem__", [](const Bmad::PointerToMultipassLord &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.ix_pass);
        if (i == 1)
          return nb::cast(s.super_lord);
        if (i == 2)
          return nb::cast(s.multi_lord);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_multipass_lord",
      &Bmad::pointer_to_multipass_lord,
      nb::arg("ele"),
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
      nb::arg("this_ele"),
      nb::arg("offset") = nb::none(),
      nb::arg("skip_beginning") = nb::none(),
      nb::arg("follow_fork") = nb::none(),
      nb::arg("ix_multipass") = nb::none(),
      R"""(Wrapper for Fortran routine pointer_to_next_ele

Parameters
----------
this_ele : EleStruct
    Starting element.

offset : int, optional
    +1 -> return next element, +2 -> element after that, etc. Can be negative. Default = +1.

skip_beginning : bool, optional
    If True then skip beginning element #0 when wrapping around. Default is False.

follow_fork : bool, optional
    If True then fork at any fork element. Default is False.

ix_multipass : int, optional
    Default = 1. Used to select the multipass branch if this_ele is a multipass_lord.

Returns
-------
next_ele : EleStruct, optional
    Element after this_ele (if offset = 1). Nullified if there is an error. EG bad this_ele.
)"""
  );
  nb::class_<Bmad::PointerToSlave>(m, "PointerToSlave", "pointer_to_slave return type")
      .def_ro("control", &Bmad::PointerToSlave::control)
      .def_ro("ix_lord_back", &Bmad::PointerToSlave::ix_lord_back)
      .def_ro("ix_control", &Bmad::PointerToSlave::ix_control)
      .def_ro("ix_ic", &Bmad::PointerToSlave::ix_ic)
      .def_ro("slave_ptr", &Bmad::PointerToSlave::slave_ptr)
      .def("__len__", [](const Bmad::PointerToSlave &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSlave &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.control);
        if (i == 1)
          return nb::cast(s.ix_lord_back);
        if (i == 2)
          return nb::cast(s.ix_control);
        if (i == 3)
          return nb::cast(s.ix_ic);
        if (i == 4)
          return nb::cast(s.slave_ptr);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_slave",
      &Bmad::pointer_to_slave,
      nb::arg("lord"),
      nb::arg("ix_slave"),
      nb::arg("slave_type") = nb::none(),
      R"""(Function to point to a slave of a lord.
Note: Ramper lords do not have any associated slaves (slaves are assigned dynamically at run time).

If slave_type = all$ (the default) the range for ix_slave is:
  1 to lord%n_slave                                 for "regular" slaves.
  lord%n_slave+1 to lord%n_slave+lord%n_slave_field for field overlap slaves.

If slave_type = field_slave$, only the field overlap slaves may be accessed and the range for ix_slave is:
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

slave_type : int, optional
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
  nb::class_<Bmad::PointerToSuperLord>(m, "PointerToSuperLord", "pointer_to_super_lord return type")
      .def_ro("control", &Bmad::PointerToSuperLord::control)
      .def_ro("ix_slave_back", &Bmad::PointerToSuperLord::ix_slave_back)
      .def_ro("ix_control", &Bmad::PointerToSuperLord::ix_control)
      .def_ro("ix_ic", &Bmad::PointerToSuperLord::ix_ic)
      .def_ro("lord_ptr", &Bmad::PointerToSuperLord::lord_ptr)
      .def("__len__", [](const Bmad::PointerToSuperLord &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSuperLord &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.control);
        if (i == 1)
          return nb::cast(s.ix_slave_back);
        if (i == 2)
          return nb::cast(s.ix_control);
        if (i == 3)
          return nb::cast(s.ix_ic);
        if (i == 4)
          return nb::cast(s.lord_ptr);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_super_lord",
      &Bmad::pointer_to_super_lord,
      nb::arg("slave"),
      nb::arg("lord_type") = nb::none(),
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
  nb::class_<Bmad::PointerToSurfaceDisplacementPt>(
      m,
      "PointerToSurfaceDisplacementPt",
      "pointer_to_surface_displacement_pt return type"
  )
      .def_ro("ix", &Bmad::PointerToSurfaceDisplacementPt::ix)
      .def_ro("iy", &Bmad::PointerToSurfaceDisplacementPt::iy)
      .def_ro("xx", &Bmad::PointerToSurfaceDisplacementPt::xx)
      .def_ro("yy", &Bmad::PointerToSurfaceDisplacementPt::yy)
      .def_ro("pt", &Bmad::PointerToSurfaceDisplacementPt::pt)
      .def("__len__", [](const Bmad::PointerToSurfaceDisplacementPt &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSurfaceDisplacementPt &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.ix);
        if (i == 1)
          return nb::cast(s.iy);
        if (i == 2)
          return nb::cast(s.xx);
        if (i == 3)
          return nb::cast(s.yy);
        if (i == 4)
          return nb::cast(s.pt);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_surface_displacement_pt",
      &Bmad::pointer_to_surface_displacement_pt,
      nb::arg("ele"),
      nb::arg("nearest"),
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("extend_grid") = nb::none(),
      R"""(Routine to point to the grid point struct associated with point (x,y).

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
  nb::class_<Bmad::PointerToSurfaceSegmentedPt>(
      m,
      "PointerToSurfaceSegmentedPt",
      "pointer_to_surface_segmented_pt return type"
  )
      .def_ro("ix", &Bmad::PointerToSurfaceSegmentedPt::ix)
      .def_ro("iy", &Bmad::PointerToSurfaceSegmentedPt::iy)
      .def_ro("xx", &Bmad::PointerToSurfaceSegmentedPt::xx)
      .def_ro("yy", &Bmad::PointerToSurfaceSegmentedPt::yy)
      .def_ro("pt", &Bmad::PointerToSurfaceSegmentedPt::pt)
      .def("__len__", [](const Bmad::PointerToSurfaceSegmentedPt &) { return 5; })
      .def("__getitem__", [](const Bmad::PointerToSurfaceSegmentedPt &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.ix);
        if (i == 1)
          return nb::cast(s.iy);
        if (i == 2)
          return nb::cast(s.xx);
        if (i == 3)
          return nb::cast(s.yy);
        if (i == 4)
          return nb::cast(s.pt);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_surface_segmented_pt",
      &Bmad::pointer_to_surface_segmented_pt,
      nb::arg("ele"),
      nb::arg("nearest"),
      nb::arg("x"),
      nb::arg("y"),
      nb::arg("extend_grid") = nb::none(),
      R"""(Routine to point to the grid point struct associated with point (x,y).

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
  nb::class_<Bmad::PointerToWakeEle>(m, "PointerToWakeEle", "pointer_to_wake_ele return type")
      .def_ro("delta_s", &Bmad::PointerToWakeEle::delta_s)
      .def_ro("wake_ele", &Bmad::PointerToWakeEle::wake_ele)
      .def("__len__", [](const Bmad::PointerToWakeEle &) { return 2; })
      .def("__getitem__", [](const Bmad::PointerToWakeEle &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.delta_s);
        if (i == 1)
          return nb::cast(s.wake_ele);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_wake_ele",
      &Bmad::pointer_to_wake_ele,
      nb::arg("ele"),
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
  nb::class_<Bmad::PointerToWall3d>(m, "PointerToWall3d", "pointer_to_wall3d return type")
      .def_ro("ds_offset", &Bmad::PointerToWall3d::ds_offset)
      .def_ro("is_branch_wall", &Bmad::PointerToWall3d::is_branch_wall)
      .def_ro("wall3d", &Bmad::PointerToWall3d::wall3d)
      .def("__len__", [](const Bmad::PointerToWall3d &) { return 3; })
      .def("__getitem__", [](const Bmad::PointerToWall3d &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.ds_offset);
        if (i == 1)
          return nb::cast(s.is_branch_wall);
        if (i == 2)
          return nb::cast(s.wall3d);
        throw nb::index_error();
      });
  m.def(
      "pointer_to_wall3d",
      &Bmad::pointer_to_wall3d,
      nb::arg("ele"),
      nb::arg("ix_wall") = nb::none(),
      R"""(Function to return a pointer to a wall3d structure associated
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
  nb::class_<Bmad::PointersToAttribute>(
      m,
      "PointersToAttribute",
      "pointers_to_attribute return type"
  )
      .def_ro("ptr_array", &Bmad::PointersToAttribute::ptr_array)
      .def_ro("err_flag", &Bmad::PointersToAttribute::err_flag)
      .def_ro("eles", &Bmad::PointersToAttribute::eles)
      .def_ro("ix_attrib", &Bmad::PointersToAttribute::ix_attrib)
      .def("__len__", [](const Bmad::PointersToAttribute &) { return 4; })
      .def("__getitem__", [](const Bmad::PointersToAttribute &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.ptr_array);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.eles);
        if (i == 3)
          return nb::cast(s.ix_attrib);
        throw nb::index_error();
      });
  m.def(
      "pointers_to_attribute",
      &Bmad::pointers_to_attribute,
      nb::arg("lat"),
      nb::arg("ele_name"),
      nb::arg("attrib_name"),
      nb::arg("do_allocation"),
      nb::arg("err_print_flag") = nb::none(),
      nb::arg("do_unlink") = nb::none(),
      R"""(Wrapper for Fortran routine pointers_to_attribute

Parameters
----------
lat : LatStruct
    Lattice.

ele_name : str
    Element name. Must be uppercase

attrib_name : str
    Attribute name. Must be uppercase. For example: "HKICK".

do_allocation : bool
    If True then do an allocation if needed. EG: The multipole An and Bn arrays need to be allocated before
    their use.

err_print_flag : bool, optional
    If present and False then suppress printing of an error message on error.

do_unlink : bool, optional

Returns
-------
ptr_array : 1D array of AllPointerStruct
    Pointer to the attribute. Size of ptr_array will be set to 0 if there is a problem.

err_flag : bool
    Set True if attribtute not found.

eles : 1D array of ElePointerStruct, optional
    Array of element pointers. size(eles) = size(ptr_array). If there are no associated lattice elements (EG
    if ele_name = 'PARTICLE_START'), eles(i).ele will be null.

ix_attrib : int, optional
    If applicable then this is the index to the attribute in the ele.value(:) array.
)"""
  );
  m.def(
      "polar_to_spinor",
      &Bmad::polar_to_spinor,
      nb::arg("polar"),
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
      nb::arg("polar"),
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
  m.def(
      "print_mesh3d",
      &Bmad::print_mesh3d,
      nb::arg("mesh3d"),
      R"""(Wrapper for Fortran routine print_mesh3d

Parameters
----------
mesh3d : Mesh3dStruct
)"""
  );
  m.def(
      "prob_x_diffuse",
      &Bmad::prob_x_diffuse,
      nb::arg("x"),
      nb::arg("d_param"),
      nb::arg("surface"),
      R"""(Contained routine to calculate integrated probability distribution in x = sin(graze_angle_out).

Parameters
----------
x : float
    sin(graze_angle_out)
)"""
  );
  nb::class_<Bmad::ProjectEmitToXyz>(m, "ProjectEmitToXyz", "project_emit_to_xyz return type")
      .def_ro("sigma_x", &Bmad::ProjectEmitToXyz::sigma_x)
      .def_ro("sigma_y", &Bmad::ProjectEmitToXyz::sigma_y)
      .def_ro("sigma_z", &Bmad::ProjectEmitToXyz::sigma_z)
      .def("__len__", [](const Bmad::ProjectEmitToXyz &) { return 3; })
      .def("__getitem__", [](const Bmad::ProjectEmitToXyz &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.sigma_x);
        if (i == 1)
          return nb::cast(s.sigma_y);
        if (i == 2)
          return nb::cast(s.sigma_z);
        throw nb::index_error();
      });
  m.def(
      "project_emit_to_xyz",
      &Bmad::project_emit_to_xyz,
      nb::arg("ring"),
      nb::arg("ix"),
      nb::arg("mode"),
      R"""(Obtains the projected x, y, and z beamsizes by building the sigma matrix
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
      "propagate_part_way",
      &Bmad::propagate_part_way,
      nb::arg("orb_start"),
      nb::arg("param"),
      nb::arg("pt"),
      nb::arg("info"),
      nb::arg("z_here"),
      nb::arg("runt"),
      R"""(Wrapper for Fortran routine propagate_part_way

Parameters
----------
orb_start : CoordStruct

param : LatParamStruct

pt : RadIntTrackPointStruct

info : RadIntInfoStruct

z_here : float

runt : EleStruct
)"""
  );
  m.def(
      "psi_prime_sca",
      &Bmad::psi_prime_sca,
      nb::arg("t"),
      nb::arg("p"),
      nb::arg("args"),
      R"""(This wraps the array-valued psi_prime function as a scalar.

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
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine ptc_bookkeeper

Parameters
----------
lat : LatStruct
    Bmad lattice.
)"""
  );
  m.def(
      "ptc_calculate_tracking_step_size",
      [](Layout &ptc_layout,
         double kl_max,
         std::optional<double> ds_max,
         BoolAlloc1D *even_steps,
         std::optional<double> r_typical,
         std::optional<double> dx_tol_bend,
         std::optional<bool> use_2nd_order,
         std::optional<FixedArray1D<Int, 2>> crossover,
         std::optional<FixedArray1D<Int, 2>> crossover_wiggler) {
        auto fn = static_cast<void (*)(
            Layout &,
            double,
            std::optional<double>,
            optional_ref<BoolAlloc1D>,
            std::optional<double>,
            std::optional<double>,
            std::optional<bool>,
            std::optional<FixedArray1D<Int, 2>>,
            std::optional<FixedArray1D<Int, 2>>
        )>(&Bmad::ptc_calculate_tracking_step_size);
        return fn(
            ptc_layout,
            kl_max,
            ds_max,
            ptr_to_opt_ref(even_steps),
            r_typical,
            dx_tol_bend,
            use_2nd_order,
            crossover,
            crossover_wiggler
        );
      },
      nb::arg("ptc_layout"),
      nb::arg("kl_max"),
      nb::arg("ds_max") = nb::none(),
      nb::arg("even_steps") = nb::none(),
      nb::arg("r_typical") = nb::none(),
      nb::arg("dx_tol_bend") = nb::none(),
      nb::arg("use_2nd_order") = nb::none(),
      nb::arg("crossover") = nb::none(),
      nb::arg("crossover_wiggler") = nb::none(),
      R"""(                                even_steps, r_typical, dx_tol_bend, use_2nd_order)

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
  nb::class_<Bmad::PtcCheckForLostParticle>(
      m,
      "PtcCheckForLostParticle",
      "ptc_check_for_lost_particle return type"
  )
      .def_ro("state", &Bmad::PtcCheckForLostParticle::state)
      .def_ro("ptc_fibre", &Bmad::PtcCheckForLostParticle::ptc_fibre)
      .def("__len__", [](const Bmad::PtcCheckForLostParticle &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcCheckForLostParticle &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.state);
        if (i == 1)
          return nb::cast(s.ptc_fibre);
        throw nb::index_error();
      });
  m.def(
      "ptc_check_for_lost_particle",
      &Bmad::ptc_check_for_lost_particle,
      nb::arg("do_reset"),
      R"""(Routine to check if a particle has been lost when tracking with PTC.

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
      nb::arg("branch"),
      nb::arg("radiation_damping_on") = nb::none(),
      R"""(Routine to calculate the closed orbit of a lattice branch using PTC.
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
  nb::class_<Bmad::PtcEmitCalc>(m, "PtcEmitCalc", "ptc_emit_calc return type")
      .def_ro("norm_mode", &Bmad::PtcEmitCalc::norm_mode)
      .def_ro("closed_orb", &Bmad::PtcEmitCalc::closed_orb)
      .def("__len__", [](const Bmad::PtcEmitCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcEmitCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.norm_mode);
        if (i == 1)
          return nb::cast(s.closed_orb);
        throw nb::index_error();
      });
  m.def(
      "ptc_emit_calc",
      &Bmad::ptc_emit_calc,
      nb::arg("ele"),
      nb::arg("sigma_mat"),
      R"""(Routine to calculate emittances, etc.

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
      "ptc_kill_map_with_radiation",
      &Bmad::ptc_kill_map_with_radiation,
      nb::arg("rad_map"),
      R"""(Routine to kill a binary file containing a ptc_rad_map_struct map

Parameters
----------
rad_map : PtcRadMapStruct
    Map with radiation included.
    This parameter is an input/output and is modified in-place.
    As an output, rad_map: Deallocated map.
)"""
  );
  m.def(
      "ptc_layouts_resplit",
      &Bmad::ptc_layouts_resplit,
      nb::arg("dKL_max"),
      nb::arg("l_max"),
      nb::arg("l_max_drift_only"),
      nb::arg("bend_dorb"),
      nb::arg("sex_dx"),
      nb::arg("even") = nb::none(),
      nb::arg("crossover") = nb::none(),
      nb::arg("crossover_wiggler") = nb::none(),
      R"""(                                                          even, crossover, crossover_wiggler)

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
      "ptc_linear_isf_calc",
      &Bmad::ptc_linear_isf_calc,
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine ptc_linear_isf_calc

Parameters
----------
branch : BranchStruct
    Lattice branch to analyze.

Returns
-------
ele_isf : 1D array of LinearEleIsfStruct
    ISF at every element.
)"""
  );
  m.def(
      "ptc_one_turn_mat_and_closed_orbit_calc",
      &Bmad::ptc_one_turn_mat_and_closed_orbit_calc,
      nb::arg("branch"),
      nb::arg("pz") = nb::none(),
      R"""(Routine to compute the transfer matrices for the individual elements and closed orbit
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
      nb::arg("iseed"),
      R"""(Wrapper for Fortran routine ptc_ran_seed_put

Parameters
----------
iseed : int
    0 -> Use system clock.
)"""
  );
  nb::class_<Bmad::PtcReadFlatFile>(m, "PtcReadFlatFile", "ptc_read_flat_file return type")
      .def_ro("err_flag", &Bmad::PtcReadFlatFile::err_flag)
      .def_ro("lat", &Bmad::PtcReadFlatFile::lat)
      .def("__len__", [](const Bmad::PtcReadFlatFile &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcReadFlatFile &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.lat);
        throw nb::index_error();
      });
  m.def(
      "ptc_read_flat_file",
      &Bmad::ptc_read_flat_file,
      nb::arg("flat_file"),
      nb::arg("create_end_marker") = nb::none(),
      nb::arg("from_mad") = nb::none(),
      R"""(Wrapper for Fortran routine ptc_read_flat_file

Parameters
----------
flat_file : 1D array of str
    Name(s) of PTC flat file(s).

create_end_marker : bool, optional
    Put a marker element named END at the end of the lattice brances? Default is True.

from_mad : bool, optional
    If True, ignore PTC specific parameters like integrator_order. Default is False. True is used when the
    fibre has been created via MAD. In this case, the PTC specific parameters may not have good values.

Returns
-------
err_flag : bool
    Set True if there is a problem.

lat : LatStruct, optional
    If present then setup a Bmad lattice.
)"""
  );
  nb::class_<Bmad::PtcReadMapWithRadiation>(
      m,
      "PtcReadMapWithRadiation",
      "ptc_read_map_with_radiation return type"
  )
      .def_ro("rad_map", &Bmad::PtcReadMapWithRadiation::rad_map)
      .def_ro("err_flag", &Bmad::PtcReadMapWithRadiation::err_flag)
      .def("__len__", [](const Bmad::PtcReadMapWithRadiation &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcReadMapWithRadiation &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.rad_map);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "ptc_read_map_with_radiation",
      &Bmad::ptc_read_map_with_radiation,
      nb::arg("file_name") = nb::none(),
      nb::arg("file_unit") = nb::none(),
      R"""(Routine to read a binary file containing a ptc_rad_map_struct map

Either file_name or file_unit must be present but not both.
File_unit is used when there are multiple maps in a file.
If file_unit is present, it is the responsibility of the calling routine to open the file beforehand
and to close the file afterwards.

Parameters
----------
file_name : str, optional
    Name of binary file.

file_unit : int, optional
    File unit number read from.

Returns
-------
rad_map : PtcRadMapStruct
    Map with radiation included.

err_flag : bool
    Set True if there is a read error.
)"""
  );
  m.def(
      "ptc_set_rf_state_for_c_normal",
      &Bmad::ptc_set_rf_state_for_c_normal,
      nb::arg("nocavity"),
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
      R"""(Routine to see if the taylor_order for PTC needs to be set/changed.
For example, for a change in bmad_com%taylor_order.
)"""
  );
  nb::class_<Bmad::PtcSetupMapWithRadiation>(
      m,
      "PtcSetupMapWithRadiation",
      "ptc_setup_map_with_radiation return type"
  )
      .def_ro("rad_map", &Bmad::PtcSetupMapWithRadiation::rad_map)
      .def_ro("err_flag", &Bmad::PtcSetupMapWithRadiation::err_flag)
      .def("__len__", [](const Bmad::PtcSetupMapWithRadiation &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcSetupMapWithRadiation &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.rad_map);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "ptc_setup_map_with_radiation",
      [](EleStruct &ele1,
         EleStruct *ele2,
         std::optional<int> map_order,
         std::optional<bool> include_damping,
         std::optional<bool> create_symplectic_map,
         CoordStruct *orbit1) {
        auto fn = static_cast<Bmad::PtcSetupMapWithRadiation (*)(
            EleStruct &,
            optional_ref<EleStruct>,
            std::optional<int>,
            std::optional<bool>,
            std::optional<bool>,
            optional_ref<CoordStruct>
        )>(&Bmad::ptc_setup_map_with_radiation);
        return fn(
            ele1,
            ptr_to_opt_ref(ele2),
            map_order,
            include_damping,
            create_symplectic_map,
            ptr_to_opt_ref(orbit1)
        );
      },
      nb::arg("ele1"),
      nb::arg("ele2") = nb::none(),
      nb::arg("map_order") = nb::none(),
      nb::arg("include_damping") = nb::none(),
      nb::arg("create_symplectic_map") = nb::none(),
      nb::arg("orbit1") = nb::none(),
      R"""(                                                                         create_symplectic_map, orbit1, err_flag)

Routine to construct a map including radiation damping and excitation.
Note: The setting of bmad_com%radiation_damping_on will determine if damping is included in the map.

ele1/ele2 must have an associated PTC layout (which can be constructed by calling lat_to_ptc_layout).

To track after calling this routine track by calling ptc_track_with_radiation.
To cleanup memory after using, call ptc_kill_map_with_radiation.
To save a map call ptc_write_map_with_radiation.
To read a saved map call ptc_read_map_with_radiation.
To set the random number seed call: ptc_ran_seed_put.

Parameters
----------
ele1 : EleStruct
    The map starts at the exit end of ele1.

ele2 : EleStruct, optional
    The map ends at the exit end of ele2. If not present, the 1-turn map will be constructed.

map_order : int, optional
    Order of the map. If not present or less than 1, the currently set order is used.

include_damping : bool, optional
    If True (default), the map will be constructed with radiation damping included. If False, the map will not
    be constructed with radiation dampling included.

create_symplectic_map : bool, optional
    If False, create a Taylor map. If True (default), create a partially inverted map which can be
    symplecitally tracked.

orbit1 : CoordStruct, optional
    Orbit at ele1 about which the map is constructed. If not present then the orbit will be computed using PTC
    tracking.

Returns
-------
rad_map : PtcRadMapStruct
    Transport map.

err_flag : bool, optional
    Set True if there is an error such as not associated PTC layout.
)"""
  );
  nb::class_<Bmad::PtcSpinCalc>(m, "PtcSpinCalc", "ptc_spin_calc return type")
      .def_ro("norm_mode", &Bmad::PtcSpinCalc::norm_mode)
      .def_ro("closed_orb", &Bmad::PtcSpinCalc::closed_orb)
      .def("__len__", [](const Bmad::PtcSpinCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcSpinCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.norm_mode);
        if (i == 1)
          return nb::cast(s.closed_orb);
        throw nb::index_error();
      });
  m.def(
      "ptc_spin_calc",
      &Bmad::ptc_spin_calc,
      nb::arg("ele"),
      nb::arg("sigma_mat"),
      R"""(Routine to equilibrium polarizations, etc.

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
      "ptc_spin_matching_calc",
      &Bmad::ptc_spin_matching_calc,
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine ptc_spin_matching_calc

Parameters
----------
branch : BranchStruct
    Lattice branch to analyze.

Returns
-------
match_info : 1D array of SpinMatchingStruct
    G-matrix and other stuff. The array will be allocated by this routine.
)"""
  );
  nb::class_<Bmad::PtcTrackAll>(m, "PtcTrackAll", "ptc_track_all return type")
      .def_ro("track_state", &Bmad::PtcTrackAll::track_state)
      .def_ro("err_flag", &Bmad::PtcTrackAll::err_flag)
      .def("__len__", [](const Bmad::PtcTrackAll &) { return 2; })
      .def("__getitem__", [](const Bmad::PtcTrackAll &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.track_state);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "ptc_track_all",
      &Bmad::ptc_track_all,
      nb::arg("branch"),
      nb::arg("orbit"),
      R"""(Routine to track from the start to the end of a lattice branch.

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
      "ptc_track_map_with_radiation",
      &Bmad::ptc_track_map_with_radiation,
      nb::arg("orbit"),
      nb::arg("rad_map"),
      nb::arg("rad_damp") = nb::none(),
      nb::arg("rad_fluct") = nb::none(),
      R"""(Routine to track through a map that includes radiation.

NOTE! Tracking without damping when the map was made with radiation (and vice versa)
will not give good results. So avoid this situation unless testing.

To construct the map, use the routine ptc_setup_map_with_radiation.
To cleanup memory after using, call ptc_kill_map_with_radiation.
To save a map call ptc_write_map_with_radiation.
To read a saved map call ptc_read_map_with_radiation.
To set the random number seed call: ptc_ran_seed_put.

Parameters
----------
orbit : CoordStruct
    Starting orbit.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Ending orbit after tracking through the map.

rad_map : PtcRadMapStruct
    Map with radiation included.

rad_damp : bool, optional
    Override the setting of bmad_com.radiation_damping_on.

rad_fluct : bool, optional
    Override the setting of bmad_com.radiation_fluctuations_on
)"""
  );
  m.def(
      "ptc_transfer_map_with_spin",
      &Bmad::ptc_transfer_map_with_spin,
      nb::arg("branch"),
      nb::arg("t_map"),
      nb::arg("s_map"),
      nb::arg("orb0"),
      nb::arg("ix1") = nb::none(),
      nb::arg("ix2") = nb::none(),
      nb::arg("one_turn") = nb::none(),
      nb::arg("unit_start") = nb::none(),
      R"""(Wrapper for Fortran routine ptc_transfer_map_with_spin

Parameters
----------
branch : BranchStruct
    Lattice branch used in the calculation.

t_map : 1D array of TaylorStruct (shape: 6)

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
    If present and False then t_map_init will be used as the starting map instead of the unit map. Default =
    True

Returns
-------
err_flag : bool
    Set True if problem like number overflow, etc.
)"""
  );
  m.def(
      "ptc_write_map_with_radiation",
      &Bmad::ptc_write_map_with_radiation,
      nb::arg("rad_map"),
      nb::arg("file_name") = nb::none(),
      nb::arg("file_unit") = nb::none(),
      R"""(Routine to create or append to a binary file containing a ptc_rad_map_struct map.

Either file_name or file_unit must be present but not both.
If file_unit is present, it is the responsibility of the calling routine to open the file beforehand
and to close the file afterwards.

Parameters
----------
rad_map : PtcRadMapStruct
    Map with radiation included.

file_name : str, optional
    Name of binary file to create.

file_unit : int, optional
    File unit number to append to.
)"""
  );
  m.def(
      "ptwo",
      &Bmad::ptwo,
      nb::arg("sigma"),
      nb::arg("t"),
      nb::arg("phi"),
      nb::arg("d_param"),
      R"""(unnormalized two-dimensional probability distribution in x and phi
polar angles relative to surface normal
azimuthal angle relative to plane of incidence (plane of incoming ray and surface normal)
1/y suppressed

Private routine.
)"""
  );
  m.def(
      "pwd_mat",
      &Bmad::pwd_mat,
      nb::arg("lat"),
      nb::arg("t6"),
      nb::arg("inductance"),
      nb::arg("sig_z"),
      R"""(Calculates potential well distortion as RF defocusing.  Calculates t6_pwd=t6.Mpwd,
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
