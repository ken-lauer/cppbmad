#include "pybmad/generated/Tao_routines_t.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyTaoDrawCurveData python_tao_draw_curve_data(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoCurveStruct &curve,
    bool have_data
) {
  Tao::tao_draw_curve_data(plot, graph, curve, have_data);
  auto py_result{PyTaoDrawCurveData{have_data}};
  return py_result;
}
PyTaoDrawHistogramData python_tao_draw_histogram_data(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoCurveStruct &curve,
    bool have_data
) {
  Tao::tao_draw_histogram_data(plot, graph, curve, have_data);
  auto py_result{PyTaoDrawHistogramData{have_data}};
  return py_result;
}
PyTaoEleShapeInfo python_tao_ele_shape_info(
    int ix_uni,
    EleStruct &ele,
    TaoEleShapeStructArray1D ele_shapes,
    std::optional<int> ix_shape_min = std::nullopt
) {
  auto _result = Tao::tao_ele_shape_info(ix_uni, ele, ele_shapes, make_opt_ref(ix_shape_min));
  auto py_result{PyTaoEleShapeInfo{_result, ix_shape_min}};
  return py_result;
}
PyTaoNextSwitch python_tao_next_switch(
    std::string line,
    CharacterAlloc1D &switch_list,
    bool return_next_word,
    std::optional<bool> neg_num_not_switch = std::nullopt,
    std::optional<bool> print_err = std::nullopt
) {
  auto _result =
      Tao::tao_next_switch(line, switch_list, return_next_word, neg_num_not_switch, print_err);
  auto py_result{PyTaoNextSwitch{_result, line}};
  return py_result;
}
PyTaoNextWord python_tao_next_word(std::string line) {
  auto _result = Tao::tao_next_word(line);
  auto py_result{PyTaoNextWord{_result, line}};
  return py_result;
}
PyTaoPointerToEleShape python_tao_pointer_to_ele_shape(
    int ix_uni,
    EleStruct &ele,
    TaoEleShapeStructArray1D ele_shape,
    std::optional<int> ix_shape_min = std::nullopt
) {
  auto _result = Tao::tao_pointer_to_ele_shape(ix_uni, ele, ele_shape, make_opt_ref(ix_shape_min));
  auto py_result{PyTaoPointerToEleShape{_result, ix_shape_min}};
  return py_result;
}
PyTaoPointerToUniverseStr python_tao_pointer_to_universe_str(
    std::string string,
    std::optional<bool> neg2_to_default = std::nullopt
) {
  auto _result = Tao::tao_pointer_to_universe(string, neg2_to_default);
  auto py_result{PyTaoPointerToUniverseStr{_result, string}};
  return py_result;
}
PyTaoRemoveBlankCharacters python_tao_remove_blank_characters(std::string str) {
  Tao::tao_remove_blank_characters(str);
  auto py_result{PyTaoRemoveBlankCharacters{str}};
  return py_result;
}

void init_Tao_routines_t(py::module &m) {
  m.def(
      "tao_abort_command_file",
      &Tao::tao_abort_command_file,
      py::arg("force_abort") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_abort_command_file

Parameters
----------
force_abort : bool, optional
    : If present and True, ignore s.global.cmd_file_abort_on_error and abort any open command files.
)"""
  );
  m.def(
      "tao_add_to_normal_mode_h_array",
      &Tao::tao_add_to_normal_mode_h_array,
      py::arg("h_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_add_to_normal_mode_h_array(h_str, h_array)

Routine to add on to the "h(:)" array holding the list of normal form
resonance driving terms to calculate.
If h_str is already in the h_array(:) list, nothing is done.

Parameters
----------
h_str : str
    Resonance driving term ID. EG: "110000"

Returns
-------
h_array : 1D array of ResonanceHStruct
    Array of resonance driving terms.
)"""
  );
  m.def(
      "tao_alias_cmd",
      &Tao::tao_alias_cmd,
      py::arg("alias"),
      py::arg("string"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_alias_cmd

Parameters
----------
alias : str
    Name of the tao command file.

string : str
    Command file arguments.
)"""
  );
  m.def(
      "tao_allocate_data_array",
      &Tao::tao_allocate_data_array,
      py::arg("u"),
      py::arg("n_data"),
      py::arg("exact") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_allocate_data_array

Parameters
----------
u : TaoUniverseStruct

n_data : int

exact : bool, optional
)"""
  );
  m.def(
      "tao_allocate_v1_var",
      &Tao::tao_allocate_v1_var,
      py::arg("n_v1"),
      py::arg("save_old"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_allocate_v1_var

Parameters
----------
n_v1 : int

save_old : bool
)"""
  );
  m.def(
      "tao_allocate_var_array",
      &Tao::tao_allocate_var_array,
      py::arg("n_var"),
      py::arg("default_good_user"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_allocate_var_array (n_var, default_good_user)

Routine to increase the s%var(:) array size.

Parameters
----------
n_var : int
    Size of s.var(:) wanted.
)"""
  );
  m.def(
      "tao_beam_emit_calc",
      &Tao::tao_beam_emit_calc,
      py::arg("plane"),
      py::arg("emit_type"),
      py::arg("ele"),
      py::arg("bunch_params"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_beam_emit_calc

Parameters
----------
plane : int
    x_plane$ or y_plane$.

emit_type : int
    Either projected_emit$ or apparent_emit$

ele : EleStruct
    Element.

bunch_params : BunchParamsStruct
    Bunch sigma matrix

Returns
-------
emit : float
    emittance.
)"""
  );
  m.def(
      "tao_beam_track",
      &Tao::tao_beam_track,
      py::arg("u"),
      py::arg("tao_lat"),
      py::arg("ix_branch"),
      py::arg("beam"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_beam_track (u, tao_lat, ix_branch, beam, calc_ok)

Routine to track a a beam of particles.

Parameters
----------
u : TaoUniverseStruct
    Universe to track through.

tao_lat : TaoLatticeStruct
    Structure containing the lattice.

ix_branch : int
    Branch index to track through.

beam : BeamStruct
    Initial beam distribution
    This parameter is an input/output and is modified in-place.
    As an output, beam: Final beam distribution.

Returns
-------
calc_ok : bool
    Set True if there were no problems, False otherwise.
)"""
  );
  m.def(
      "tao_beam_track_endpoint",
      &Tao::tao_beam_track_endpoint,
      py::arg("ele_id"),
      py::arg("lat"),
      py::arg("branch_str"),
      py::arg("where"),
      py::arg("u"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_beam_track_endpoint

Parameters
----------
ele_id : str
    Name or index of the element.

lat : LatStruct
    Lattice.

branch_str : str
    Branch where the tracking is done. '' => Branch not specified.

where : str
    'TRACK_END', 'TRACK_START', etc.. Used for error messages.

u : TaoUniverseStruct
    Universe beam is being tracked in.

Returns
-------
ele : EleStruct, optional
    Pointer to the track endpoint element. Nullified if error. Note: blank ele_id is handled if "where"
    contains 'END' or 'START'
)"""
  );
  m.def(
      "tao_branch_index",
      &Tao::tao_branch_index,
      py::arg("ix_branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_branch_index

Parameters
----------
ix_branch : int
    Nominal branch number.

Returns
-------
ix_this : int
    Branch number.
)"""
  );
  m.def(
      "tao_calc_data_at_s_pts",
      &Tao::tao_calc_data_at_s_pts,
      py::arg("tao_lat"),
      py::arg("curve"),
      py::arg("comp_sign"),
      py::arg("good"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_calc_data_at_s_pts

Parameters
----------
tao_lat : TaoLatticeStruct

curve : TaoCurveStruct

comp_sign : float

good : 1D array of bool
)"""
  );
  m.def(
      "tao_call_cmd",
      &Tao::tao_call_cmd,
      py::arg("file_name"),
      py::arg("cmd_arg") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_call_cmd

Parameters
----------
file_name : str
    Name of the tao command file.

cmd_arg : 1D array of str, optional
    Command file arguments.
)"""
  );
  m.def(
      "tao_cbar_wave_anal",
      &Tao::tao_cbar_wave_anal,
      py::arg("plot"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_cbar_wave_anal (plot)
)"""
  );
  m.def(
      "tao_change_ele",
      &Tao::tao_change_ele,
      py::arg("ele_name"),
      py::arg("attrib_name"),
      py::arg("num_str"),
      py::arg("update"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_change_ele (ele_name, attrib_name, num_str, update, err_flag)

Routine to change a variable in the model lattice.

Parameters
----------
ele_name : str
    Name of variable or element.

attrib_name : str
    Attribute name of element.

num_str : str
    Change in value. A '@' signifies a absolute set. A 'd' signifies a set relative design.

Returns
-------
err_flag : bool
    logical, Set true if there is an error, false otherwise.
)"""
  );
  m.def(
      "tao_change_tune",
      &Tao::tao_change_tune,
      py::arg("branch_str"),
      py::arg("mask_str"),
      py::arg("print_list"),
      py::arg("dqa_str"),
      py::arg("dqb_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_change_tune (branch_str, mask_str, print_list, dqa_str, dqb_str, err_flag)

Parameters
----------
branch_str : str
    List of branches to apply tune set to.

mask_str : str
    List of quadrupoles to veto.

print_list : bool
    If True, print a list of elements varied and coefficients.

dqa_str : str
    Expression for dQa tune.

dqb_str : str
    Expression for dQb tune.

Returns
-------
err_flag : bool
    logical, Set true if there is an error, false otherwise.
)"""
  );
  m.def(
      "tao_change_var",
      &Tao::tao_change_var,
      py::arg("name"),
      py::arg("num_str"),
      py::arg("silent"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_change_var (name, num_str, silent, err_flag)

Routine to change a variable in the model lattice.

Parameters
----------
name : str
    Name of variable or element.

num_str : str
    Change in value. A '@' signifies a absolute set. A 'd' signifies a set relative design.

silent : bool
    If True then do not print any info.

Returns
-------
err_flag : bool
    logical, Set true if there is an error, false otherwise.
)"""
  );
  m.def(
      "tao_change_z_tune",
      &Tao::tao_change_z_tune,
      py::arg("branch_str"),
      py::arg("dq_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_change_z_tune (branch_str, dq_str, err_flag)

Parameters
----------
branch_str : str
    List of branches to apply tune set to.

dq_str : str
    Expression for dQc tune.

Returns
-------
err_flag : bool
    logical, Set true if there is an error, false otherwise.
)"""
  );
  m.def(
      "tao_chrom_calc_needed",
      &Tao::tao_chrom_calc_needed,
      py::arg("data_type"),
      py::arg("data_source"),
      py::arg("do_chrom"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_chrom_calc_needed

Parameters
----------
data_type : str

data_source : str

do_chrom : bool
)"""
  );
  m.def(
      "tao_clear_cmd",
      &Tao::tao_clear_cmd,
      py::arg("cmd_line"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_clear_cmd

Parameters
----------
cmd_line : str
    Should be set to 'maps'.
)"""
  );
  m.def(
      "tao_clip_cmd",
      &Tao::tao_clip_cmd,
      py::arg("gang"),
      py::arg("where"),
      py::arg("value1"),
      py::arg("value2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_clip_cmd

Parameters
----------
gang : bool
    Gang all data d1 arrays together.

where : str
    Graph() to clip. Eg: 'top:x'

value1 : float

value2 : float
)"""
  );
  m.def(
      "tao_close_command_file",
      &Tao::tao_close_command_file,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_close_command_file
)"""
  );
  m.def(
      "tao_cmd_history_record",
      &Tao::tao_cmd_history_record,
      py::arg("cmd"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_cmd_history_record (cmd)

Subroutine to record a cmd in the command history stack
)"""
  );
  m.def(
      "tao_cmd_split",
      &Tao::tao_cmd_split,
      py::arg("cmd_line"),
      py::arg("n_word"),
      py::arg("cmd_word"),
      py::arg("extra_words_is_error"),
      py::arg("separator") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_cmd_split (cmd_line, n_word, cmd_word, extra_words_is_error, err, separator)

This routine splits the command line into "words" (everything between separators).

Parameters
----------
cmd_line : str
    The command line.

n_word : int
    Maximum number of words to split command line into.

cmd_word : 1D array of str
    The individual words.

extra_words_is_error : bool
    are extra words allowed at the end? If True then err argument is set True. If False then cmd_word(n_word)
    will contain everything after the n_word-1 word.

separator : str, optional
    a list of characters that, besides a blank space, signify a word boundary.

Returns
-------
err : bool
    error in splitting words For example: separator = '-+' cmd_line = 'model-design' cmd_word(1) = 'model'
    cmd_word(2) = '-' cmd_word(3) = 'design'

Notes
-----
Anything between single or double quotes is treated as a single word.
Quoted words have quote marks removed.
Whitespace or a separator inside of "{}", "()", or "[]" is ignored.
Whitespace after or before a comma is ignored.
)"""
  );
  m.def(
      "tao_command",
      &Tao::tao_command,
      py::arg("command_line"),
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_command

Parameters
----------
command_line : str
    command line

err : bool

Returns
-------
err_is_fatal : bool
    Set True on non-recoverable error. False otherwise
)"""
  );
  m.def(
      "tao_constraint_type_name",
      &Tao::tao_constraint_type_name,
      py::arg("datum"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_constraint_type_name

Parameters
----------
datum : TaoDataStruct
    Datum

Returns
-------
datum_name : str
    Appropriate name.
)"""
  );
  m.def(
      "tao_control_tree_list",
      &Tao::tao_control_tree_list,
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_control_tree_list

Parameters
----------
ele : EleStruct
    Lattice element to start at.

Returns
-------
tree : 1D array of ElePointerStruct
    Array of elements.
)"""
  );
  m.def(
      "tao_count_strings",
      &Tao::tao_count_strings,
      py::arg("string"),
      py::arg("pattern"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_count_strings

Parameters
----------
string : str
    the string to look at

pattern : str
    the search pattern

Returns
-------
num : int
    number of occurances
)"""
  );
  m.def(
      "tao_create_plot_window",
      &Tao::tao_create_plot_window,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_create_plot_window ()

Subroutine to create the plot window.
This soubroutine knows not to create a second window if one already exists.
)"""
  );
  m.def(
      "tao_curve_beam_ellipse_setup",
      &Tao::tao_curve_beam_ellipse_setup,
      py::arg("curve"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_curve_beam_ellipse_setup

Parameters
----------
curve : TaoCurveStruct
)"""
  );
  m.def(
      "tao_curve_check_universe",
      &Tao::tao_curve_check_universe,
      py::arg("curve"),
      py::arg("uni"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_curve_check_universe (curve, uni) result (is_ok)

Routine to check if the universe associated with a curve exists and is on.

Parameters
----------
curve : TaoCurveStruct
    Curve to check.
    This parameter is an input/output and is modified in-place.
    As an output, curve: Curve.valid set to False if needed.

uni : TaoUniverseStruct
    Associated universe

Returns
-------
is_ok : bool
    Set True if associated universe exists and is on.
)"""
  );
  m.def(
      "tao_curve_data_setup",
      &Tao::tao_curve_data_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::arg("curve"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_curve_data_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct

curve : TaoCurveStruct
)"""
  );
  m.def(
      "tao_curve_datum_calc",
      &Tao::tao_curve_datum_calc,
      py::arg("eles"),
      py::arg("plot"),
      py::arg("curve"),
      py::arg("who"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_curve_datum_calc (eles, plot, curve, who)

Routine to calculate datum values.
The values are calculated at the end of each eles(:)%ele element.

Parameters
----------
eles : 1D array of ElePointerStruct
    Array of elements.

plot : TaoPlotStruct

curve : TaoCurveStruct
    This parameter is an input/output and is modified in-place.
    As an output, curve: Structure holding the datum values

who : str
    Where to put the data. Either: "SYMBOL" or "LINE".
)"""
  );
  m.def(
      "tao_curve_ele_ref",
      &Tao::tao_curve_ele_ref,
      py::arg("curve"),
      py::arg("point_to_ele_ref"),
      py::arg("ele_track"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_curve_ele_ref

Parameters
----------
curve : TaoCurveStruct
    Curve with ref ele.

point_to_ele_ref : bool

ele_track : EleStruct
)"""
  );
  m.def(
      "tao_curve_ix_uni",
      &Tao::tao_curve_ix_uni,
      py::arg("curve"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_curve_ix_uni

Parameters
----------
curve : TaoCurveStruct
    Curve.

Returns
-------
ix_uni : int
    Universe index.
)"""
  );
  m.def(
      "tao_curve_name",
      &Tao::tao_curve_name,
      py::arg("curve"),
      py::arg("use_region") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_curve_name

Parameters
----------
curve : TaoCurveStruct
    Curve

use_region : bool, optional
    If present and True then use the region name instead of the plot name. Region name is 'NULL_REGION' if
    there is no assocaited region.

Returns
-------
curve_name : str
    Appropriate name.
)"""
  );
  py::class_<Tao::TaoCurveRmsCalc, std::unique_ptr<Tao::TaoCurveRmsCalc>>(
      m,
      "TaoCurveRmsCalc",
      "tao_curve_rms_calc return type"
  )
      .def_readonly("rms", &Tao::TaoCurveRmsCalc::rms)
      .def_readonly("mean", &Tao::TaoCurveRmsCalc::mean)
      .def("__len__", [](const Tao::TaoCurveRmsCalc &) { return 2; })
      .def("__getitem__", [](const Tao::TaoCurveRmsCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.rms);
        if (i == 1)
          return py::cast(s.mean);
        throw py::index_error();
      });
  m.def(
      "tao_curve_rms_calc",
      &Tao::tao_curve_rms_calc,
      py::arg("curve"),
      py::arg("who"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_curve_rms_calc

Parameters
----------
curve : TaoCurveStruct
    Curve to analyze.

who : str
    "LINE" or "SYMBOL".

Returns
-------
rms : float
    RMS. -1 => Curve has no data.

mean : float
    Mean.
)"""
  );
  m.def(
      "tao_d2_d1_name",
      &Tao::tao_d2_d1_name,
      py::arg("d1"),
      py::arg("show_universe") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_d2_d1_name

Parameters
----------
d1 : TaoD1DataStruct
    Data array.

show_universe : bool, optional
    Show the datum's universe. Default is True.

Returns
-------
d2_d1_name : str
    Appropriate name.
)"""
  );
  m.def(
      "tao_d2_data_stuffit",
      &Tao::tao_d2_data_stuffit,
      py::arg("u"),
      py::arg("d2_name"),
      py::arg("n_d1_data"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_d2_data_stuffit

Parameters
----------
u : TaoUniverseStruct

d2_name : str

n_d1_data : int
)"""
  );
  m.def(
      "tao_data_check",
      &Tao::tao_data_check,
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_data_check

Parameters
----------
err : bool
)"""
  );
  m.def(
      "tao_data_coupling_init",
      &Tao::tao_data_coupling_init,
      py::arg("branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_data_coupling_init

Parameters
----------
branch : BranchStruct
    New lattice branch.
)"""
  );
  m.def(
      "tao_data_sanity_check",
      &Tao::tao_data_sanity_check,
      py::arg("datum"),
      py::arg("print_err"),
      py::arg("default_data_type"),
      py::arg("uni") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_data_sanity_check

Parameters
----------
datum : TaoDataStruct
    Datum to check.

print_err : bool
    Print error message if data is not valid?

default_data_type : str
    Default data type associated with the datum's d2 structure.

uni : TaoUniverseStruct, optional
    Universe to use instead of datum.d1.d2.ix_universe

Returns
-------
is_valid : bool
    True if internally consistent.
)"""
  );
  m.def(
      "tao_data_show_use",
      &Tao::tao_data_show_use,
      py::arg("d2_data"),
      py::arg("lines") = py::none(),
      py::arg("nl") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_data_show_use

Parameters
----------
d2_data : TaoD2DataStruct

lines : 1D array of str, optional

nl : int, optional
)"""
  );
  m.def(
      "tao_data_type_substitute",
      &Tao::tao_data_type_substitute,
      py::arg("template_"),
      py::arg("curve"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_data_type_substitute (template, str_out, curve, graph)

Routine substitute the appropriate data type string for instances of "#ref" and
"#comp" in template.

Additionally, if template does not have a "|" character,
the string "|" + component will be added at the end of str_out.

Parameters
----------
curve : TaoCurveStruct
    curve.ele_ref_name is substituted for all instances of "#ref".

graph : TaoGraphStruct

Returns
-------
str_out : str
    String with substitutions.
)"""
  );
  m.def(
      "tao_data_useit_plot_calc",
      &Tao::tao_data_useit_plot_calc,
      py::arg("curve"),
      py::arg("graph"),
      py::arg("data"),
      py::arg("check_s_position"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_data_useit_plot_calc (curve, graph, data, check_s_position, most_invalid)

Routine to set the data for plotting.

Parameters
----------
curve : TaoCurveStruct
    tao_curve_struct

graph : TaoGraphStruct
    tao_graph_struct

data : 1D array of TaoDataStruct

check_s_position : bool
    If present and True then veto data that does not have an s-position.

Returns
-------
most_invalid : str
    String documenting biggest invalid data problem.
)"""
  );
  m.def(
      "tao_datum_has_associated_ele",
      &Tao::tao_datum_has_associated_ele,
      py::arg("data_type"),
      py::arg("branch_geometry") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_datum_has_associated_ele

Parameters
----------
data_type : str
    Type of data.

branch_geometry : int, optional
    Geometry of the associated lattice branch. open$ or closed$.

Returns
-------
has_associated_ele : int
)"""
  );
  py::class_<Tao::TaoDatumIntegrate, std::unique_ptr<Tao::TaoDatumIntegrate>>(
      m,
      "TaoDatumIntegrate",
      "tao_datum_integrate return type"
  )
      .def_readonly("valid_value", &Tao::TaoDatumIntegrate::valid_value)
      .def_readonly("why_invalid", &Tao::TaoDatumIntegrate::why_invalid)
      .def_readonly("result", &Tao::TaoDatumIntegrate::result)
      .def("__len__", [](const Tao::TaoDatumIntegrate &) { return 3; })
      .def("__getitem__", [](const Tao::TaoDatumIntegrate &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.valid_value);
        if (i == 1)
          return py::cast(s.why_invalid);
        if (i == 2)
          return py::cast(s.result);
        throw py::index_error();
      });
  m.def(
      "tao_datum_integrate",
      &Tao::tao_datum_integrate,
      py::arg("datum"),
      py::arg("branch"),
      py::arg("s_pos"),
      py::arg("values"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_datum_integrate (datum, branch, s_pos, values, valid_value, why_invalid) result (result)

Routine to calculate the integral, rms, or average of an array of values associated with a datum.

Parameters
----------
datum : TaoDataStruct
    Datum under consideration.

branch : BranchStruct
    Associated lattice branch.

s_pos : 1D array of float
    Array of s-positions of the values.

values : 1D array of float
    Array of values.

Returns
-------
valid_value : bool
    Set false if, for example, all s_pos(:) are the same.

why_invalid : str
    Information string if there is a problem.

result : float
    Integral, rms, or average depending upon datum.merit_type.
)"""
  );
  m.def(
      "tao_datum_name",
      &Tao::tao_datum_name,
      py::arg("datum"),
      py::arg("show_universe") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_datum_name

Parameters
----------
datum : TaoDataStruct
    Datum

show_universe : bool, optional
    Show the datum's universe. Default is True.

Returns
-------
datum_name : str
    Appropriate name.
)"""
  );
  m.def(
      "tao_datum_s_position",
      &Tao::tao_datum_s_position,
      py::arg("datum"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_datum_s_position (datum, ele) result (s_pos)

Routine to calculate the longitudinal position associated with a datum.

Parameters
----------
datum : TaoDataStruct
    Datum under conideration.

ele : EleStruct
    Associated lattice element.

Returns
-------
s_pos : float
    Associated longitudinal position.
)"""
  );
  m.def(
      "tao_de_optimizer",
      &Tao::tao_de_optimizer,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_de_optimizer

Returns
-------
abort : bool
    Set True if an user stop signal detected.
)"""
  );
  m.def(
      "tao_deallocate_plot_cache",
      &Tao::tao_deallocate_plot_cache,
      py::arg("plot_cache"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_deallocate_plot_cache

Parameters
----------
plot_cache : 1D array of TaoPlotCacheStruct
)"""
  );
  m.def(
      "tao_deallocate_tree",
      &Tao::tao_deallocate_tree,
      py::arg("tree"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_deallocate_tree (tree)

Routine to deallocate tree%node(:) and everything below it

Parameters
----------
tree : TaoEvalNodeStruct
    Root of tree to deallocate.
    This parameter is an input/output and is modified in-place.
    As an output, tree: Deallocated tree.
)"""
  );
  m.def(
      "tao_destroy_plot_window",
      &Tao::tao_destroy_plot_window,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_destroy_plot_window
)"""
  );
  m.def(
      "tao_dmerit_calc",
      &Tao::tao_dmerit_calc,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_dmerit_calc ()
)"""
  );
  m.def(
      "tao_dmodel_dvar_calc",
      &Tao::tao_dmodel_dvar_calc,
      py::arg("force_calc"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_dModel_dVar_calc (force_calc, err_flag)

Subroutine to calculate the dModel_dVar derivative matrix.

Parameters
----------
force_calc : bool
    If true then force recalculation of the matrix. If False then only calculate matrix if it doesn't exist.

Returns
-------
err_flag : bool
    Set true if there is an error. False otherwise.
)"""
  );
  m.def(
      "tao_do_wire_scan",
      &Tao::tao_do_wire_scan,
      py::arg("ele"),
      py::arg("theta"),
      py::arg("beam"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_do_wire_scan (ele, wire_params, theta, beam) result (moment)

Returns the beam's second moment using the wire along the specified angle.
Keep in mind that the actual correlation axis is 90 degrees off of the
wire angle

This simulates a fast wire scanner that performs the scan over only one
bunch. Obviously, this isn't realistic. Any dynamic effects will not be
accounted for!

Parameters
----------
ele : EleStruct

theta : float
    wire angle wrt x axis (in degrees)

beam : BeamStruct
    contains the beam distribution

Returns
-------
moment : float
    second moment along axis specified by angle.
)"""
  );
  m.def(
      "tao_draw_beam_chamber_wall",
      &Tao::tao_draw_beam_chamber_wall,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(NOTE: THIS ROUTINE IS NOT CURRENTLY ACITVE (NOT CALLED BY ANY OTHER ROUTINE).

Subroutine tao_draw_beam_chamber_wall (plot, graph)

Routine to draw the beam chamber wall.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  py::class_<PyTaoDrawCurveData, std::unique_ptr<PyTaoDrawCurveData>>(
      m,
      "TaoDrawCurveData",
      "tao_draw_curve_data return type"
  )
      .def_readonly("have_data", &PyTaoDrawCurveData::have_data)
      .def("__len__", [](const PyTaoDrawCurveData &) { return 1; })
      .def("__getitem__", [](const PyTaoDrawCurveData &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.have_data);
        throw py::index_error();
      });
  m.def(
      "tao_draw_curve_data",
      &python_tao_draw_curve_data,
      py::arg("plot"),
      py::arg("graph"),
      py::arg("curve"),
      py::arg("have_data"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_curve_data (plot, graph, curve, have_data)

Routine to draw a graph with data and/or variable curves.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph containing the curve.

curve : TaoCurveStruct
    Curve to draw.

have_data : bool
    Intitial state.
    This parameter is an input/output and is modified in-place.
    As an output, have_data: Is there any data to plot? Set True if so.

Returns
-------
have_data : bool
    Intitial state.
    This parameter is an input/output and is modified in-place.
    As an output, have_data: Is there any data to plot? Set True if so.
)"""
  );
  m.def(
      "tao_draw_ele_for_floor_plan",
      &Tao::tao_draw_ele_for_floor_plan,
      py::arg("plot"),
      py::arg("graph"),
      py::arg("tao_lat"),
      py::arg("ele"),
      py::arg("ele_shape"),
      py::arg("label_name"),
      py::arg("offset1"),
      py::arg("offset2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_ele_for_floor_plan (plot, graph, tao_lat, ele, ele_shape, label_name, offset1, offset2)

Routine to draw one lattice element or one datum location for the floor plan graph.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.

tao_lat : TaoLatticeStruct
    Lattice containing the element.

ele : EleStruct
    Element to draw.

ele_shape : TaoEleShapeStruct
    Shape to draw from s.plot_page.floor_plan.ele_shape(:) array. Will be NULL if no associated shape for this
    element.

label_name : str
    Shape label.

offset1 : float
    Transverse distances used to scale the drawing of the element shape.

offset2 : float
    Transverse distances used to scale the drawing of the element shape.
)"""
  );
  m.def(
      "tao_draw_floor_plan",
      &Tao::tao_draw_floor_plan,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_floor_plan (plot, graph)

Routine to draw a floor plan graph.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  m.def(
      "tao_draw_graph_axes",
      &Tao::tao_draw_graph_axes,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_graph_axes (plot, graph)

Routine to draw a just the graph part of a data graph.
The calling routine takes care of drawing any curves.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  py::class_<PyTaoDrawHistogramData, std::unique_ptr<PyTaoDrawHistogramData>>(
      m,
      "TaoDrawHistogramData",
      "tao_draw_histogram_data return type"
  )
      .def_readonly("have_data", &PyTaoDrawHistogramData::have_data)
      .def("__len__", [](const PyTaoDrawHistogramData &) { return 1; })
      .def("__getitem__", [](const PyTaoDrawHistogramData &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.have_data);
        throw py::index_error();
      });
  m.def(
      "tao_draw_histogram_data",
      &python_tao_draw_histogram_data,
      py::arg("plot"),
      py::arg("graph"),
      py::arg("curve"),
      py::arg("have_data"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_histogram_data (plot, graph, curve, have_data)

Routine to draw a graph with data and/or variable histograms.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph containing the histogram.

curve : TaoCurveStruct
    Histogram to draw.

have_data : bool
    Intitial state.
    This parameter is an input/output and is modified in-place.
    As an output, have_data: Is there any data to plot? Set True if so.

Returns
-------
have_data : bool
    Intitial state.
    This parameter is an input/output and is modified in-place.
    As an output, have_data: Is there any data to plot? Set True if so.
)"""
  );
  m.def(
      "tao_draw_lat_layout",
      &Tao::tao_draw_lat_layout,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_lat_layout (plot, graph)

Routine to draw a lattice layout graph.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  m.def(
      "tao_draw_plots",
      &Tao::tao_draw_plots,
      py::arg("do_clear") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_draw_plots (do_clear)

Subroutine to draw the plots on the plot window.

Parameters
----------
do_clear : bool, optional
    If present and False then call qp_clear_page. This argument is used when drawing PS or GIF.
)"""
  );
  py::class_<
      Tao::TaoEleGeometryWithMisalignments,
      std::unique_ptr<Tao::TaoEleGeometryWithMisalignments>>(
      m,
      "TaoEleGeometryWithMisalignments",
      "tao_ele_geometry_with_misalignments return type"
  )
      .def_readonly("valid_value", &Tao::TaoEleGeometryWithMisalignments::valid_value)
      .def_readonly("why_invalid", &Tao::TaoEleGeometryWithMisalignments::why_invalid)
      .def_readonly("value", &Tao::TaoEleGeometryWithMisalignments::value)
      .def("__len__", [](const Tao::TaoEleGeometryWithMisalignments &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEleGeometryWithMisalignments &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.valid_value);
        if (i == 1)
          return py::cast(s.why_invalid);
        if (i == 2)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "tao_ele_geometry_with_misalignments",
      &Tao::tao_ele_geometry_with_misalignments,
      py::arg("datum"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_ele_geometry_with_misalignments (datum, ele, valid_value, why_invalid) result (value)

Routine to evaluate a floor position with misalignments at a given element.
This routine is private and not for general use.

Parameters
----------
datum : TaoDataStruct
    Datum info

ele : EleStruct
    Lattice element to evaluate at.

Returns
-------
valid_value : bool
    Was able to evalute the datum?

why_invalid : str
    If not valid, why not.

value : float
    Datum value.
)"""
  );
  py::class_<PyTaoEleShapeInfo, std::unique_ptr<PyTaoEleShapeInfo>>(
      m,
      "TaoEleShapeInfo",
      "tao_ele_shape_info return type"
  )
      .def_readonly("e_shape", &PyTaoEleShapeInfo::e_shape)
      .def_readonly("label_name", &PyTaoEleShapeInfo::label_name)
      .def_readonly("y1", &PyTaoEleShapeInfo::y1)
      .def_readonly("y2", &PyTaoEleShapeInfo::y2)
      .def_readonly("ix_shape_min", &PyTaoEleShapeInfo::ix_shape_min)
      .def("__len__", [](const PyTaoEleShapeInfo &) { return 5; })
      .def("__getitem__", [](const PyTaoEleShapeInfo &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.e_shape);
        if (i == 1)
          return py::cast(s.label_name);
        if (i == 2)
          return py::cast(s.y1);
        if (i == 3)
          return py::cast(s.y2);
        if (i == 4)
          return py::cast(s.ix_shape_min);
        throw py::index_error();
      });
  m.def(
      "tao_ele_shape_info",
      &python_tao_ele_shape_info,
      py::arg("ix_uni"),
      py::arg("ele"),
      py::arg("ele_shapes"),
      py::arg("ix_shape_min") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_ele_shape_info

Parameters
----------
ix_uni : int
    Universe index.

ele : EleStruct
    Lattice element.

ele_shapes : 1D array of TaoEleShapeStruct
    Array of shapes to search.

ix_shape_min : int, optional
    Index of minimum ele_shape(:) index to start search from. Default is 1.
    This parameter is an input/output and is modified in-place.
    As an output, ix_shape_min: Ele_shape(

Returns
-------
label_name : str
    Label name.

y1 : float
    shape transverse sizes.

y2 : float
    shape transverse sizes.

e_shape : TaoEleShapeStruct, optional
    element shape. Will be nullified if no associated shape.

ix_shape_min : int, optional
    Index of minimum ele_shape(:) index to start search from. Default is 1.
    This parameter is an input/output and is modified in-place.
    As an output, ix_shape_min: Ele_shape(
)"""
  );
  py::class_<Tao::TaoEvalFloorOrbit, std::unique_ptr<Tao::TaoEvalFloorOrbit>>(
      m,
      "TaoEvalFloorOrbit",
      "tao_eval_floor_orbit return type"
  )
      .def_readonly("valid_value", &Tao::TaoEvalFloorOrbit::valid_value)
      .def_readonly("why_invalid", &Tao::TaoEvalFloorOrbit::why_invalid)
      .def_readonly("value", &Tao::TaoEvalFloorOrbit::value)
      .def("__len__", [](const Tao::TaoEvalFloorOrbit &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvalFloorOrbit &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.valid_value);
        if (i == 1)
          return py::cast(s.why_invalid);
        if (i == 2)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "tao_eval_floor_orbit",
      &Tao::tao_eval_floor_orbit,
      py::arg("datum"),
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("bunch_params"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_eval_floor_orbit (datum, ele, orbit, bunch_params, valid_value, why_invalid) result (value)

Routine to evaluate a floor_orbit datum at a given element.
This routine is private and not for general use.

Parameters
----------
datum : TaoDataStruct
    Datum info

ele : EleStruct
    Lattice element to evaluate at.

orbit : CoordStruct
    Particle orbit at element.

bunch_params : BunchParamsStruct
    Bunch parameters at element.

Returns
-------
valid_value : bool
    Was able to evalute the datum?

why_invalid : str
    If not valid, why not.

value : float
    Datum value.
)"""
  );
  py::class_<Tao::TaoEvaluateADatum, std::unique_ptr<Tao::TaoEvaluateADatum>>(
      m,
      "TaoEvaluateADatum",
      "tao_evaluate_a_datum return type"
  )
      .def_readonly("datum_value", &Tao::TaoEvaluateADatum::datum_value)
      .def_readonly("valid_value", &Tao::TaoEvaluateADatum::valid_value)
      .def_readonly("why_invalid", &Tao::TaoEvaluateADatum::why_invalid)
      .def("__len__", [](const Tao::TaoEvaluateADatum &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvaluateADatum &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.datum_value);
        if (i == 1)
          return py::cast(s.valid_value);
        if (i == 2)
          return py::cast(s.why_invalid);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_a_datum",
      &Tao::tao_evaluate_a_datum,
      py::arg("datum"),
      py::arg("u"),
      py::arg("tao_lat"),
      py::arg("called_from_lat_calc") = py::none(),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_a_datum

Parameters
----------
datum : TaoDataStruct
    What type of datum

u : TaoUniverseStruct
    Which universe to use.

tao_lat : TaoLatticeStruct
    Lattice to use.

called_from_lat_calc : bool, optional
    Default is False. If true, prevents infinite loop of this routine calling tao_lattice_calc

print_err : bool, optional
    Default is True. If False, do not print an error message.

Returns
-------
datum_value : float
    Value of the datum.

valid_value : bool
    Set false when there is a problem. Set true otherwise.

why_invalid : str, optional
    Tells why datum value is invalid.
)"""
  );
  py::class_<Tao::TaoEvaluateDatumAtS, std::unique_ptr<Tao::TaoEvaluateDatumAtS>>(
      m,
      "TaoEvaluateDatumAtS",
      "tao_evaluate_datum_at_s return type"
  )
      .def_readonly("err_str", &Tao::TaoEvaluateDatumAtS::err_str)
      .def_readonly("bad_datum", &Tao::TaoEvaluateDatumAtS::bad_datum)
      .def_readonly("value", &Tao::TaoEvaluateDatumAtS::value)
      .def("__len__", [](const Tao::TaoEvaluateDatumAtS &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvaluateDatumAtS &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err_str);
        if (i == 1)
          return py::cast(s.bad_datum);
        if (i == 2)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_datum_at_s",
      &Tao::tao_evaluate_datum_at_s,
      py::arg("datum"),
      py::arg("tao_lat"),
      py::arg("ele"),
      py::arg("ele_ref"),
      py::arg("valid_value"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_evaluate_datum_at_s (datum, tao_lat, ele, ele_ref, valid_value, err_str, bad_datum) result(value)

Routine to evaluate a datum at a given s-position in the lattice

Parameters
----------
datum : TaoDataStruct
    Datum to evaluate.

tao_lat : TaoLatticeStruct

ele : EleStruct
    Evaluation element.

ele_ref : EleStruct
    Reference element.

valid_value : bool
    True if evaluation was sucessful. False if not.

Returns
-------
err_str : str
    Error string for printing an error message.

bad_datum : bool
    True -> datum is malformed. False -> Could evaluate or evaluation problem was not due to the datum itself
    (EG: the lattice was unstable).

value : float
    Datum value.
)"""
  );
  py::class_<Tao::TaoEvaluateElementParameters, std::unique_ptr<Tao::TaoEvaluateElementParameters>>(
      m,
      "TaoEvaluateElementParameters",
      "tao_evaluate_element_parameters return type"
  )
      .def_readonly("err", &Tao::TaoEvaluateElementParameters::err)
      .def_readonly("values", &Tao::TaoEvaluateElementParameters::values)
      .def_readonly("info", &Tao::TaoEvaluateElementParameters::info)
      .def("__len__", [](const Tao::TaoEvaluateElementParameters &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvaluateElementParameters &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.err);
        if (i == 1)
          return py::cast(s.values);
        if (i == 2)
          return py::cast(s.info);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_element_parameters",
      &Tao::tao_evaluate_element_parameters,
      py::arg("param_name"),
      py::arg("print_err"),
      py::arg("dflt_source"),
      py::arg("dflt_ele") = py::none(),
      py::arg("dflt_component") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::arg("eval_point") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_element_parameters

Parameters
----------
param_name : str
    parameter name.

print_err : bool
    Print error message?

dflt_source : str
    Default source

dflt_ele : EleStruct, optional
    Default element if not specified by param_name.

dflt_component : str, optional
    Default component

dflt_uni : int, optional
    Default universe to use.

eval_point : int, optional

Returns
-------
err : bool
    True if there is an error in syntax. False otherwise

values : 1D array of float
    Array of datum values.

info : 1D array of TaoExpressionInfoStruct, optional
)"""
  );
  py::class_<Tao::TaoEvaluateExpression, std::unique_ptr<Tao::TaoEvaluateExpression>>(
      m,
      "TaoEvaluateExpression",
      "tao_evaluate_expression return type"
  )
      .def_readonly("value", &Tao::TaoEvaluateExpression::value)
      .def_readonly("err_flag", &Tao::TaoEvaluateExpression::err_flag)
      .def_readonly("info", &Tao::TaoEvaluateExpression::info)
      .def_readonly("stack", &Tao::TaoEvaluateExpression::stack)
      .def("__len__", [](const Tao::TaoEvaluateExpression &) { return 4; })
      .def("__getitem__", [](const Tao::TaoEvaluateExpression &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.info);
        if (i == 3)
          return py::cast(s.stack);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_expression",
      &Tao::tao_evaluate_expression,
      py::arg("expression"),
      py::arg("n_size"),
      py::arg("use_good_user"),
      py::arg("print_err") = py::none(),
      py::arg("dflt_component") = py::none(),
      py::arg("dflt_source") = py::none(),
      py::arg("dflt_ele_ref") = py::none(),
      py::arg("dflt_ele_start") = py::none(),
      py::arg("dflt_ele") = py::none(),
      py::arg("dflt_dat_or_var_index") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::arg("dflt_eval_point") = py::none(),
      py::arg("dflt_s_offset") = py::none(),
      py::arg("dflt_orbit") = py::none(),
      py::arg("datum") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_expression

Parameters
----------
expression : str
    Arithmetic expression.

n_size : int
    Size of the value array. If the expression evaluates to a a scalar, each value in the value array will get
    this value. If n_size = 0, the natural size is determined by the expression itself.

use_good_user : bool
    Use the good_user logical in evaluating good(:)

print_err : bool, optional
    If False then supress evaluation error messages. This does not affect syntax error messages. Default is
    True.

dflt_component : str, optional
    Component to use if not specified in the expression. 'model' (default), 'base', or 'design'.

dflt_source : str, optional
    Default source ('lat', 'data', etc.). Default is ''.

dflt_ele_ref : EleStruct, optional
    Default reference element.

dflt_ele_start : EleStruct, optional
    Default start element for ranges.

dflt_ele : EleStruct, optional
    Default element to evaluate at.

dflt_dat_or_var_index : str, optional
    Default datum or variable index to use.

dflt_uni : int, optional
    Default universe to use. If 0 or not present, use viewed universe.

dflt_eval_point : int, optional
    Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

dflt_s_offset : float, optional
    Default offset of eval_point. Default = 0.

dflt_orbit : CoordStruct, optional
    Default orbit to evaluate at.

datum : TaoDataStruct, optional
    If present, check to see that the expression does not depend upon a datum that will be evaluated after
    this datum. If so, this is an error.

Returns
-------
value : 1D array of float
    Value of arithmetic expression.

err_flag : bool
    True on an error. EG: Invalid expression. A divide by zero is not an error but good(:) will be set to
    False.

info : 1D array of TaoExpressionInfoStruct, optional
    Is the value valid?, etc. Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
    orbit.x[23]|good_user is False.

stack : 1D array of TaoEvalNodeStruct, optional
    Array of nodes of variable names. This is useful to check what datums or variables are used in the
    expression.
)"""
  );
  py::class_<Tao::TaoEvaluateExpressionNew, std::unique_ptr<Tao::TaoEvaluateExpressionNew>>(
      m,
      "TaoEvaluateExpressionNew",
      "tao_evaluate_expression_new return type"
  )
      .def_readonly("value", &Tao::TaoEvaluateExpressionNew::value)
      .def_readonly("err_flag", &Tao::TaoEvaluateExpressionNew::err_flag)
      .def_readonly("info", &Tao::TaoEvaluateExpressionNew::info)
      .def_readonly("stack", &Tao::TaoEvaluateExpressionNew::stack)
      .def("__len__", [](const Tao::TaoEvaluateExpressionNew &) { return 4; })
      .def("__getitem__", [](const Tao::TaoEvaluateExpressionNew &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.info);
        if (i == 3)
          return py::cast(s.stack);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_expression_new",
      &Tao::tao_evaluate_expression_new,
      py::arg("expression"),
      py::arg("n_size"),
      py::arg("use_good_user"),
      py::arg("print_err") = py::none(),
      py::arg("dflt_component") = py::none(),
      py::arg("dflt_source") = py::none(),
      py::arg("dflt_ele_ref") = py::none(),
      py::arg("dflt_ele_start") = py::none(),
      py::arg("dflt_ele") = py::none(),
      py::arg("dflt_dat_or_var_index") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::arg("dflt_eval_point") = py::none(),
      py::arg("dflt_s_offset") = py::none(),
      py::arg("dflt_orbit") = py::none(),
      py::arg("datum") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_expression_new

Parameters
----------
expression : str
    Arithmetic expression.

n_size : int
    Size of the value array. If the expression evaluates to a a scalar, each value in the value array will get
    this value. If n_size = 0, the natural size is determined by the expression itself.

use_good_user : bool
    Use the good_user logical in evaluating good(:)

print_err : bool, optional
    If False then supress evaluation error messages. This does not affect syntax error messages. Default is
    True.

dflt_component : str, optional
    Component to use if not specified in the expression. 'model' (default), 'base', or 'design'.

dflt_source : str, optional
    Default source ('lat', 'data', etc.). Default is ''.

dflt_ele_ref : EleStruct, optional
    Default reference element.

dflt_ele_start : EleStruct, optional
    Default start element for ranges.

dflt_ele : EleStruct, optional
    Default element to evaluate at.

dflt_dat_or_var_index : str, optional
    Default datum or variable index to use.

dflt_uni : int, optional
    Default universe to use. If 0 or not present, use viewed universe.

dflt_eval_point : int, optional
    Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

dflt_s_offset : float, optional
    Default offset of eval_point. Default = 0.

dflt_orbit : CoordStruct, optional
    Default orbit to evaluate at.

datum : TaoDataStruct, optional
    If present, check to see that the expression does not depend upon a datum that will be evaluated after
    this datum. If so, this is an error.

Returns
-------
value : 1D array of float
    Value of arithmetic expression.

err_flag : bool
    True on an error. EG: Invalid expression. A divide by zero is not an error but good(:) will be set to
    False.

info : 1D array of TaoExpressionInfoStruct, optional
    Is the value valid?, etc. Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
    orbit.x[23]|good_user is False.

stack : 1D array of TaoEvalNodeStruct, optional
    Array of nodes of variable names. This is useful to check what datums or variables are used in the
    expression.
)"""
  );
  py::class_<Tao::TaoEvaluateExpressionOld, std::unique_ptr<Tao::TaoEvaluateExpressionOld>>(
      m,
      "TaoEvaluateExpressionOld",
      "tao_evaluate_expression_old return type"
  )
      .def_readonly("value", &Tao::TaoEvaluateExpressionOld::value)
      .def_readonly("err_flag", &Tao::TaoEvaluateExpressionOld::err_flag)
      .def_readonly("info", &Tao::TaoEvaluateExpressionOld::info)
      .def_readonly("stack", &Tao::TaoEvaluateExpressionOld::stack)
      .def("__len__", [](const Tao::TaoEvaluateExpressionOld &) { return 4; })
      .def("__getitem__", [](const Tao::TaoEvaluateExpressionOld &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.info);
        if (i == 3)
          return py::cast(s.stack);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_expression_old",
      &Tao::tao_evaluate_expression_old,
      py::arg("expression"),
      py::arg("n_size"),
      py::arg("use_good_user"),
      py::arg("print_err") = py::none(),
      py::arg("dflt_component") = py::none(),
      py::arg("dflt_source") = py::none(),
      py::arg("dflt_ele_ref") = py::none(),
      py::arg("dflt_ele_start") = py::none(),
      py::arg("dflt_ele") = py::none(),
      py::arg("dflt_dat_or_var_index") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::arg("dflt_eval_point") = py::none(),
      py::arg("dflt_s_offset") = py::none(),
      py::arg("dflt_orbit") = py::none(),
      py::arg("datum") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_expression_old

Parameters
----------
expression : str
    Arithmetic expression.

n_size : int
    Size of the value array. If the expression evaluates to a a scalar, each value in the value array will get
    this value. If n_size = 0, the natural size is determined by the expression itself.

use_good_user : bool
    Use the good_user logical in evaluating good(:)

print_err : bool, optional
    If False then supress evaluation error messages. This does not affect syntax error messages. Default is
    True.

dflt_component : str, optional
    Component to use if not specified in the expression. 'model' (default), 'base', or 'design'.

dflt_source : str, optional
    Default source ('lat', 'data', etc.). Default is ''.

dflt_ele_ref : EleStruct, optional
    Default reference element.

dflt_ele_start : EleStruct, optional
    Default start element for ranges.

dflt_ele : EleStruct, optional
    Default element to evaluate at.

dflt_dat_or_var_index : str, optional
    Default datum or variable index to use.

dflt_uni : int, optional
    Default universe to use. If 0 or not present, use viewed universe.

dflt_eval_point : int, optional
    Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

dflt_s_offset : float, optional
    Default offset of eval_point. Default = 0.

dflt_orbit : CoordStruct, optional
    Default orbit to evaluate at.

datum : TaoDataStruct, optional
    If present, check to see that the expression does not depend upon a datum that will be evaluated after
    this datum. If so, this is an error.

Returns
-------
value : 1D array of float
    Value of arithmetic expression.

err_flag : bool
    True on an error. EG: Invalid expression. A divide by zero is not an error but good(:) will be set to
    False.

info : 1D array of TaoExpressionInfoStruct, optional
    Is the value valid?, etc. Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
    orbit.x[23]|good_user is False.

stack : 1D array of TaoEvalNodeStruct, optional
    Array of nodes of variable names. This is useful to check what datums or variables are used in the
    expression.
)"""
  );
  py::class_<Tao::TaoEvaluateLatOrBeamData, std::unique_ptr<Tao::TaoEvaluateLatOrBeamData>>(
      m,
      "TaoEvaluateLatOrBeamData",
      "tao_evaluate_lat_or_beam_data return type"
  )
      .def_readonly("err", &Tao::TaoEvaluateLatOrBeamData::err)
      .def_readonly("values", &Tao::TaoEvaluateLatOrBeamData::values)
      .def("__len__", [](const Tao::TaoEvaluateLatOrBeamData &) { return 2; })
      .def("__getitem__", [](const Tao::TaoEvaluateLatOrBeamData &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err);
        if (i == 1)
          return py::cast(s.values);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_lat_or_beam_data",
      &Tao::tao_evaluate_lat_or_beam_data,
      py::arg("data_name"),
      py::arg("print_err"),
      py::arg("default_source"),
      py::arg("dflt_ele_ref") = py::none(),
      py::arg("dflt_ele_start") = py::none(),
      py::arg("dflt_ele") = py::none(),
      py::arg("dflt_component") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::arg("dflt_eval_point") = py::none(),
      py::arg("dflt_s_offset") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(! private tao_scratch_values_calc, tao_eval_floor_orbit, tao_ele_geometry_with_misalignments

 Subroutine tao_evaluate_lat_or_beam_data (err, data_name, values, print_err, default_source, default_source,
               dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni, dflt_eval_point, dflt_s_offset)

 Routine to evaluate data with a lat or beam source of the form:
     <universe>@lat::<data_type>[<ix_ele_start>&<ix_ele>]|<component>

Parameters
----------
data_name : str
    data name.

print_err : bool
    Print error message?

dflt_ele_ref : EleStruct, optional
    Default reference element.

dflt_ele_start : EleStruct, optional
    Default start element.

dflt_ele : EleStruct, optional
    Default element to evaluate at.

dflt_component : str, optional
    Default component: 'model' (default), 'base', or 'design'.

dflt_uni : int, optional
    Default universe to use.

dflt_eval_point : int, optional
    Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

dflt_s_offset : float, optional
    Default offset of eval_point. Default = 0.

Returns
-------
err : bool
    True if there is an error. False otherwise.

values : 1D array of float
    Array of datum valuse.
)"""
  );
  py::class_<Tao::TaoEvaluateStackOld, std::unique_ptr<Tao::TaoEvaluateStackOld>>(
      m,
      "TaoEvaluateStackOld",
      "tao_evaluate_stack_old return type"
  )
      .def_readonly("value", &Tao::TaoEvaluateStackOld::value)
      .def_readonly("err_flag", &Tao::TaoEvaluateStackOld::err_flag)
      .def("__len__", [](const Tao::TaoEvaluateStackOld &) { return 2; })
      .def("__getitem__", [](const Tao::TaoEvaluateStackOld &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_stack_old",
      &Tao::tao_evaluate_stack_old,
      py::arg("stack"),
      py::arg("n_size_in"),
      py::arg("use_good_user"),
      py::arg("print_err"),
      py::arg("expression"),
      py::arg("info_in") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_evaluate_stack_old (stack, n_size_in, use_good_user, value, info, err_flag, print_err, expression)

Routine to evaluate an expression stack.

Parameters
----------
stack : 1D array of TaoEvalNodeStruct
    Expression stack

n_size_in : int
    Desired array size. If the expression evaluates to a a scalar, each value in the value array will get this
    value. If n_size = 0, the natural size is determined by the expression itself.

use_good_user : bool
    Use the good_user logical in evaluating good(:)

print_err : bool
    If False then supress evaluation error messages. This does not affect syntax error messages. Default is
    True.

expression : str
    Original expression. Used for error messages.

Returns
-------
value : 1D array of float
    Value of arithmetic expression.

err_flag : bool
    True on error. False otherwise
)"""
  );
  py::class_<Tao::TaoEvaluateTree, std::unique_ptr<Tao::TaoEvaluateTree>>(
      m,
      "TaoEvaluateTree",
      "tao_evaluate_tree return type"
  )
      .def_readonly("value", &Tao::TaoEvaluateTree::value)
      .def_readonly("err_flag", &Tao::TaoEvaluateTree::err_flag)
      .def("__len__", [](const Tao::TaoEvaluateTree &) { return 2; })
      .def("__getitem__", [](const Tao::TaoEvaluateTree &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "tao_evaluate_tree",
      &Tao::tao_evaluate_tree,
      py::arg("tao_tree"),
      py::arg("n_size"),
      py::arg("use_good_user"),
      py::arg("print_err"),
      py::arg("expression"),
      py::arg("info_in") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_tree

Parameters
----------
tao_tree : TaoEvalNodeStruct
    Expression tree

n_size : int

use_good_user : bool
    Use the good_user logical in evaluating good(:)

print_err : bool
    If False then supress evaluation error messages. This does not affect syntax error messages. Default is
    True.

expression : str
    Original expression. Used for error messages.

info_in : 1D array of TaoExpressionInfoStruct, optional

Returns
-------
value : 1D array of float
    Value(s) of the arithmetic expression.

err_flag : bool
    True on error. False otherwise
)"""
  );
  m.def(
      "tao_evaluate_tune",
      &Tao::tao_evaluate_tune,
      py::arg("q_str"),
      py::arg("q0"),
      py::arg("delta_input"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_evaluate_tune

Parameters
----------
q_str : str
    String expression.

q0 : float
    Default to use if q_str evaluates to zero. Also used to set the integer part of the tune.

delta_input : bool
    If true then qa_str and qb_str are deltas from present tune.

Returns
-------
q_val : float
    Tune value. Set zero if there is an error.
)"""
  );
  m.def(
      "tao_expression_hash_substitute",
      &Tao::tao_expression_hash_substitute,
      py::arg("expression_in"),
      py::arg("eval_ele") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_expression_hash_substitute(expression_in, expression_out, eval_ele)

Routine to, in the expression, substitute the evaluation lattice element name in place
of hash ("#") characters. Care is taken to only do this where it makes sense.
For example, "Q1##3" where here "##3" means the third instance of Q1, does not qualify.

Specifically, a substitution will be done if the character before the hash and the
character after are one of:
  [,]-*+/:|@<>, or a blank character, or the beginning or end of the expression

Parameters
----------
expression_in : str
    Expression.

eval_ele : EleStruct, optional
    Evaluation element name to substitute in. If not present, expression will not be modified.

Returns
-------
expression_out : str
    Expression with substitutions made.
)"""
  );
  m.def(
      "tao_expression_tree_to_string",
      &Tao::tao_expression_tree_to_string,
      py::arg("tree"),
      py::arg("include_root") = py::none(),
      py::arg("n_node") = py::none(),
      py::arg("parent") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_expression_tree_to_string (tree, include_root, n_node, parent) result(str_out)

Routine to convert an expression tree to a expression string.

Parameters
----------
tree : TaoEvalNodeStruct
    Tree to print.

include_root : bool, optional
    Default is True. If True, do not inculde in the output string the root node. Note: If the root node is of
    type root$, this node is always ignored.

n_node : int, optional
    Internal use only. Used with recursive calls.

parent : TaoEvalNodeStruct, optional
    Internal use only. Used with recusive calls.

Returns
-------
str_out : str
    Expression string.
)"""
  );
  py::class_<Tao::TaoFindPlotRegion, std::unique_ptr<Tao::TaoFindPlotRegion>>(
      m,
      "TaoFindPlotRegion",
      "tao_find_plot_region return type"
  )
      .def_readonly("err", &Tao::TaoFindPlotRegion::err)
      .def_readonly("region", &Tao::TaoFindPlotRegion::region)
      .def("__len__", [](const Tao::TaoFindPlotRegion &) { return 2; })
      .def("__getitem__", [](const Tao::TaoFindPlotRegion &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err);
        if (i == 1)
          return py::cast(s.region);
        throw py::index_error();
      });
  m.def(
      "tao_find_plot_region",
      &Tao::tao_find_plot_region,
      py::arg("where"),
      py::arg("print_flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_find_plot_region

Parameters
----------
where : str
    Region name.

print_flag : bool, optional
    If present and False then surpress error messages. Default is True.

Returns
-------
err : bool
    Set True on error. False otherwise.

region : TaoPlotRegionStruct, optional
    Region found.
)"""
  );
  m.def(
      "tao_fixer",
      &Tao::tao_fixer,
      py::arg("switch_"),
      py::arg("word1"),
      py::arg("word2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_fixer

Parameters
----------
word1 : str
    First word of command.

word2 : str
    Secton word of command.
)"""
  );
  py::class_<Tao::TaoFloorToScreen, std::unique_ptr<Tao::TaoFloorToScreen>>(
      m,
      "TaoFloorToScreen",
      "tao_floor_to_screen return type"
  )
      .def_readonly("x_screen", &Tao::TaoFloorToScreen::x_screen)
      .def_readonly("y_screen", &Tao::TaoFloorToScreen::y_screen)
      .def("__len__", [](const Tao::TaoFloorToScreen &) { return 2; })
      .def("__getitem__", [](const Tao::TaoFloorToScreen &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.x_screen);
        if (i == 1)
          return py::cast(s.y_screen);
        throw py::index_error();
      });
  m.def(
      "tao_floor_to_screen",
      &Tao::tao_floor_to_screen,
      py::arg("graph"),
      py::arg("r_floor"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_floor_to_screen

Parameters
----------
graph : TaoGraphStruct
    Graph defining the projection plane.

r_floor : 1D array of float (shape: 3)

Returns
-------
x_screen : float
    x-coordinate of projected point.

y_screen : float
    y-coordinate of projected point.
)"""
  );
  m.def(
      "tao_floor_to_screen_coords",
      &Tao::tao_floor_to_screen_coords,
      py::arg("graph"),
      py::arg("floor"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_floor_to_screen_coords

Parameters
----------
graph : TaoGraphStruct
    Graph defining the projection plane.

floor : FloorPositionStruct
    3D coordinate.

Returns
-------
screen : FloorPositionStruct
    Projected point
)"""
  );
  m.def(
      "tao_geodesic_lm_optimizer",
      &Tao::tao_geodesic_lm_optimizer,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_geodesic_lm_optimizer (abort)

Routine to minimize the merit function by varying variables until
the "data" as calculated from the model matches the measured data.

This subroutine is a wrapper for the "geodesic"
Levenburg - Marquardt method.

Returns
-------
abort : bool
    Set True if an user stop signal detected.
)"""
  );
  py::class_<Tao::TaoGetData, std::unique_ptr<Tao::TaoGetData>>(
      m,
      "TaoGetData",
      "tao_get_data return type"
  )
      .def_readonly("data_value", &Tao::TaoGetData::data_value)
      .def_readonly("data_weight", &Tao::TaoGetData::data_weight)
      .def_readonly("data_meas_value", &Tao::TaoGetData::data_meas_value)
      .def_readonly("data_ix_dModel", &Tao::TaoGetData::data_ix_dModel)
      .def("__len__", [](const Tao::TaoGetData &) { return 4; })
      .def("__getitem__", [](const Tao::TaoGetData &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.data_value);
        if (i == 1)
          return py::cast(s.data_weight);
        if (i == 2)
          return py::cast(s.data_meas_value);
        if (i == 3)
          return py::cast(s.data_ix_dModel);
        throw py::index_error();
      });
  m.def(
      "tao_get_data",
      &Tao::tao_get_data,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_get_data (data_value, data_weight, data_meas_value, dat_ix_dModel)

Subroutine to get the values of the data used in optimization and put them
in an array. The data is ordered starting with the first universe

Returns
-------
data_value : 1D array of float, optional
    Data model values.

data_weight : 1D array of float, optional
    Data weights in the merit function.

data_meas_value : 1D array of float, optional
    Data values when the data was taken.

data_ix_dModel : 1D array of int, optional
    Data ix_dModel indices
)"""
  );
  py::class_<Tao::TaoGetOptVars, std::unique_ptr<Tao::TaoGetOptVars>>(
      m,
      "TaoGetOptVars",
      "tao_get_opt_vars return type"
  )
      .def_readonly("var_value", &Tao::TaoGetOptVars::var_value)
      .def_readonly("var_step", &Tao::TaoGetOptVars::var_step)
      .def_readonly("var_delta", &Tao::TaoGetOptVars::var_delta)
      .def_readonly("var_weight", &Tao::TaoGetOptVars::var_weight)
      .def_readonly("var_ix", &Tao::TaoGetOptVars::var_ix)
      .def_readonly("ignore_if_weight_is_zero", &Tao::TaoGetOptVars::ignore_if_weight_is_zero)
      .def_readonly("ignore_if_not_limited", &Tao::TaoGetOptVars::ignore_if_not_limited)
      .def("__len__", [](const Tao::TaoGetOptVars &) { return 7; })
      .def("__getitem__", [](const Tao::TaoGetOptVars &s, int i) -> py::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return py::cast(s.var_value);
        if (i == 1)
          return py::cast(s.var_step);
        if (i == 2)
          return py::cast(s.var_delta);
        if (i == 3)
          return py::cast(s.var_weight);
        if (i == 4)
          return py::cast(s.var_ix);
        if (i == 5)
          return py::cast(s.ignore_if_weight_is_zero);
        if (i == 6)
          return py::cast(s.ignore_if_not_limited);
        throw py::index_error();
      });
  m.def(
      "tao_get_opt_vars",
      &Tao::tao_get_opt_vars,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_get_opt_vars

Returns
-------
var_value : 1D array of float, optional
    Variable model values.

var_step : 1D array of float, optional
    Variable step sizes.

var_delta : 1D array of float, optional
    Variable Merit deltas.

var_weight : 1D array of float, optional
    Variable weights in the merit function.

var_ix : 1D array of int, optional
    Variable s.var(:) indexes

ignore_if_weight_is_zero : bool, optional
    If present and True then ignore all variables whose merit weight is zero.

ignore_if_not_limited : bool, optional
    If present and True then ignore all variables with limit constraint that are not limited.
)"""
  );
  m.def(
      "tao_get_user_input",
      &Tao::tao_get_user_input,
      py::arg("prompt_str") = py::none(),
      py::arg("wait_flag") = py::none(),
      py::arg("cmd_in") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in)

Subroutine to get the next Tao command. In order of precedence, input may come from:
  1) s%com%cmd string (if s%com%use_cmd_here is set to True).
     Used for recalling commands from the history stack.
  3) A command file.
  4) The cmd_in argument (if present). Used, for example, when interfacing with Python.
  5) The terminal.

Note: A saved command string is present if a prior input string contained multiple commands.
For example, the following string is read from a command file or terminal or passed via cmd_in:
        "show ele 1; set opti de; run"
Then cmd_out would be "show ele 1" and "set opti de; run" would be saved for the next call to this routine.

Note: In single character mode, the input precedence order is ignored and input is taken from the terminal.

Parameters
----------
prompt_str : str, optional
    Primpt string to print at terminal. If not present then s.global.prompt_string will be used.

wait_flag : bool, optional
    Used for single mode: Wait state for get_a_char call.

cmd_in : str, optional
    Command to be used in place getting user input.

Returns
-------
cmd_out : str
    Command from the user.
)"""
  );
  m.def(
      "tao_graph_controller_setup",
      &Tao::tao_graph_controller_setup,
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_controller_setup

Parameters
----------
graph : TaoGraphStruct
)"""
  );
  m.def(
      "tao_graph_data_setup",
      &Tao::tao_graph_data_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_data_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  m.def(
      "tao_graph_data_slice_setup",
      &Tao::tao_graph_data_slice_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_data_slice_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  m.def(
      "tao_graph_dynamic_aperture_setup",
      &Tao::tao_graph_dynamic_aperture_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_dynamic_aperture_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  m.def(
      "tao_graph_histogram_setup",
      &Tao::tao_graph_histogram_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_histogram_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  m.def(
      "tao_graph_name",
      &Tao::tao_graph_name,
      py::arg("graph"),
      py::arg("use_region") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_name

Parameters
----------
graph : TaoGraphStruct
    Graph

use_region : bool, optional
    If present and True then use the region name instead of the plot name. Region name is 'NULL_REGION' if
    there is no assocaited region.

Returns
-------
graph_name : str
    Appropriate name.
)"""
  );
  m.def(
      "tao_graph_phase_space_setup",
      &Tao::tao_graph_phase_space_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_phase_space_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  py::class_<Tao::TaoGraphSMinMaxCalc, std::unique_ptr<Tao::TaoGraphSMinMaxCalc>>(
      m,
      "TaoGraphSMinMaxCalc",
      "tao_graph_s_min_max_calc return type"
  )
      .def_readonly("s_min", &Tao::TaoGraphSMinMaxCalc::s_min)
      .def_readonly("s_max", &Tao::TaoGraphSMinMaxCalc::s_max)
      .def("__len__", [](const Tao::TaoGraphSMinMaxCalc &) { return 2; })
      .def("__getitem__", [](const Tao::TaoGraphSMinMaxCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.s_min);
        if (i == 1)
          return py::cast(s.s_max);
        throw py::index_error();
      });
  m.def(
      "tao_graph_s_min_max_calc",
      &Tao::tao_graph_s_min_max_calc,
      py::arg("graph"),
      py::arg("branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_graph_s_min_max_calc(graph, branch, s_min, s_max)

Routine to calculate min and max for a graph when plot%x_axis_type is set to "s".

Parameters
----------
graph : TaoGraphStruct
    Graph to calculate for.

branch : BranchStruct
    Associated lattice branch.

Returns
-------
s_min : float
    Graph min. May be negative with graph.allow_wrap_around = T.

s_max : float
    Graph max.
)"""
  );
  m.def(
      "tao_graph_setup",
      &Tao::tao_graph_setup,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_graph_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  py::class_<Tao::TaoHelp, std::unique_ptr<Tao::TaoHelp>>(m, "TaoHelp", "tao_help return type")
      .def_readonly("lines", &Tao::TaoHelp::lines)
      .def_readonly("n_lines", &Tao::TaoHelp::n_lines)
      .def("__len__", [](const Tao::TaoHelp &) { return 2; })
      .def("__getitem__", [](const Tao::TaoHelp &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.lines);
        if (i == 1)
          return py::cast(s.n_lines);
        throw py::index_error();
      });
  m.def(
      "tao_help",
      &Tao::tao_help,
      py::arg("what1"),
      py::arg("what2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_help

Parameters
----------
what1 : str
    command to query. EG: "show".

what2 : str
    subcommand to query. EG: "element".

Returns
-------
lines : 1D array of str, optional
    If present then the output will be put in this string array instead of printing to the terminal.

n_lines : int, optional
    Must be present if lines is present. Number of lines used in the lines(:) array.
)"""
  );
  m.def(
      "tao_init",
      &Tao::tao_init,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init

Returns
-------
err_flag : bool
    Set Treu if there is an error. False otherwise.
)"""
  );
  m.def(
      "tao_init_beam_in_universe",
      &Tao::tao_init_beam_in_universe,
      py::arg("u"),
      py::arg("beam_init"),
      py::arg("track_start"),
      py::arg("track_end"),
      py::arg("comb_ds_save"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init_beam_in_universe

Parameters
----------
u : TaoUniverseStruct

beam_init : BeamInitStruct

track_start : str

track_end : str

comb_ds_save : float
)"""
  );
  m.def(
      "tao_init_beams",
      &Tao::tao_init_beams,
      py::arg("init_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_init_beams (init_file)

Subroutine to initialize beam stuff.

Parameters
----------
init_file : str
    Tao initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  m.def(
      "tao_init_data",
      &Tao::tao_init_data,
      py::arg("data_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_init_data (data_file)

Subroutine to initialize the tao data structures.

Parameters
----------
data_file : str
    Tao data initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  m.def(
      "tao_init_data_end_stuff",
      &Tao::tao_init_data_end_stuff,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init_data_end_stuff
)"""
  );
  m.def(
      "tao_init_data_in_universe",
      &Tao::tao_init_data_in_universe,
      py::arg("u"),
      py::arg("n_d2_add"),
      py::arg("keep_existing_data") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init_data_in_universe

Parameters
----------
u : TaoUniverseStruct

n_d2_add : int

keep_existing_data : bool, optional
)"""
  );
  m.def(
      "tao_init_dynamic_aperture",
      &Tao::tao_init_dynamic_aperture,
      py::arg("init_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_init_dynamic_aperture (init_file)

Routine to initalize dynamic aperture simulations.

Parameters
----------
init_file : str
    File setting dynamic_aperture parameters.
)"""
  );
  py::class_<Tao::TaoInitFindElements, std::unique_ptr<Tao::TaoInitFindElements>>(
      m,
      "TaoInitFindElements",
      "tao_init_find_elements return type"
  )
      .def_readonly("eles", &Tao::TaoInitFindElements::eles)
      .def_readonly("found_one", &Tao::TaoInitFindElements::found_one)
      .def("__len__", [](const Tao::TaoInitFindElements &) { return 2; })
      .def("__getitem__", [](const Tao::TaoInitFindElements &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.eles);
        if (i == 1)
          return py::cast(s.found_one);
        throw py::index_error();
      });
  m.def(
      "tao_init_find_elements",
      &Tao::tao_init_find_elements,
      py::arg("u"),
      py::arg("search_string"),
      py::arg("attribute") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init_find_elements

Parameters
----------
u : TaoUniverseStruct
    Universe to search

search_string : str
    What to search for

attribute : str, optional
    Check that attribute of element is free to vary.

Returns
-------
eles : 1D array of ElePointerStruct
    List of matching elements. Size is zero if no elements found.

found_one : bool, optional
    Set True if a matching element is found. However: Not set if no matching element found.
)"""
  );
  m.def(
      "tao_init_global",
      &Tao::tao_init_global,
      py::arg("init_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_init_global (init_file)

Subroutine to initialize the tao global structures.

Parameters
----------
init_file : str
    Tao initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  m.def(
      "tao_init_lattice",
      &Tao::tao_init_lattice,
      py::arg("lat_file"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init_lattice

Parameters
----------
lat_file : str

err_flag : bool
)"""
  );
  m.def(
      "tao_init_plotting",
      &Tao::tao_init_plotting,
      py::arg("plot_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_init_plotting

Parameters
----------
plot_file : str
)"""
  );
  m.def(
      "tao_init_variables",
      &Tao::tao_init_variables,
      py::arg("var_file"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_init_variables (var_file)

Subroutine to initialize the tao variable structures.

Parameters
----------
var_file : str
    Tao variable initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  py::class_<Tao::TaoInjectBeam, std::unique_ptr<Tao::TaoInjectBeam>>(
      m,
      "TaoInjectBeam",
      "tao_inject_beam return type"
  )
      .def_readonly("beam", &Tao::TaoInjectBeam::beam)
      .def_readonly("init_ok", &Tao::TaoInjectBeam::init_ok)
      .def("__len__", [](const Tao::TaoInjectBeam &) { return 2; })
      .def("__getitem__", [](const Tao::TaoInjectBeam &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.beam);
        if (i == 1)
          return py::cast(s.init_ok);
        throw py::index_error();
      });
  m.def(
      "tao_inject_beam",
      &Tao::tao_inject_beam,
      py::arg("u"),
      py::arg("model"),
      py::arg("ix_branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_inject_beam (u, model, ix_branch, beam, init_ok)

This will initialize the beam for a given lattice branch.

Trying to inject a beam of one species into a branch with a different ref species
(example: electron bunch into photon branch) is problematical. To avoid problems, Tao
will set not inject (init_ok = False) if there is a mismatch.

Parameters
----------
u : TaoUniverseStruct
    Universe containing the lattice.

model : TaoLatticeStruct
    Universe parameters.

ix_branch : int
    Lattice branch index to inject into.

Returns
-------
beam : BeamStruct
    Initial beam.

init_ok : bool
    Set False if there are problems. True otherwise.
)"""
  );
  m.def(
      "tao_inject_particle",
      &Tao::tao_inject_particle,
      py::arg("u"),
      py::arg("model"),
      py::arg("ix_branch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_inject_particle

Parameters
----------
u : TaoUniverseStruct

model : TaoLatticeStruct

ix_branch : int
)"""
  );
  py::class_<Tao::TaoIsValidName, std::unique_ptr<Tao::TaoIsValidName>>(
      m,
      "TaoIsValidName",
      "tao_is_valid_name return type"
  )
      .def_readonly("why_invalid", &Tao::TaoIsValidName::why_invalid)
      .def_readonly("is_valid", &Tao::TaoIsValidName::is_valid)
      .def("__len__", [](const Tao::TaoIsValidName &) { return 2; })
      .def("__getitem__", [](const Tao::TaoIsValidName &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.why_invalid);
        if (i == 1)
          return py::cast(s.is_valid);
        throw py::index_error();
      });
  m.def(
      "tao_is_valid_name",
      &Tao::tao_is_valid_name,
      py::arg("name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_is_valid_name

Parameters
----------
name : str
    Name to be checked.

Returns
-------
why_invalid : str
    Why invalid description.

is_valid : bool
    True if valid. False otherwise.
)"""
  );
  m.def(
      "tao_json_cmd",
      &Tao::tao_json_cmd,
      py::arg("input_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_json_cmd

Parameters
----------
input_str : str
    What to show.
)"""
  );
  m.def(
      "tao_key_info_to_str",
      &Tao::tao_key_info_to_str,
      py::arg("ix_key"),
      py::arg("ix_min_key"),
      py::arg("ix_max_key"),
      py::arg("key_str"),
      py::arg("header_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_key_info_to_str

Parameters
----------
ix_key : int

ix_min_key : int

ix_max_key : int

key_str : str

header_str : str
)"""
  );
  m.def(
      "tao_lat_bookkeeper",
      &Tao::tao_lat_bookkeeper,
      py::arg("u"),
      py::arg("tao_lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lat_bookkeeper

Parameters
----------
u : TaoUniverseStruct

tao_lat : TaoLatticeStruct

Returns
-------
err_flag : bool
    Set True if there is a problem. False otherwise.
)"""
  );
  m.def(
      "tao_lat_emit_calc",
      &Tao::tao_lat_emit_calc,
      py::arg("plane"),
      py::arg("emit_type"),
      py::arg("ele"),
      py::arg("modes"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lat_emit_calc

Parameters
----------
plane : int
    x_plane$ or y_plane$.

emit_type : int
    Either projected_emit$ or apparent_emit$

ele : EleStruct
    Element holding the Twiss and coupling parameters.

modes : NormalModesStruct
    Structure holding the emittances

Returns
-------
emit : float
    emittance.
)"""
  );
  m.def(
      "tao_lat_sigma_calc_needed",
      &Tao::tao_lat_sigma_calc_needed,
      py::arg("data_type"),
      py::arg("data_source"),
      py::arg("do_lat_sigma"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lat_sigma_calc_needed

Parameters
----------
data_type : str

data_source : str

do_lat_sigma : bool
)"""
  );
  m.def(
      "tao_lat_sigma_track",
      &Tao::tao_lat_sigma_track,
      py::arg("tao_lat"),
      py::arg("ix_branch"),
      py::arg("print_err") = py::none(),
      py::arg("force_calc") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_lat_sigma_track (tao_lat, calc_ok, ix_branch, print_err, force_calc)

Routine to track the 6x6 sigma matrix through the lattice using the lattice linear transfer matrices.

Parameters
----------
tao_lat : TaoLatticeStruct
    Structure containing the lattice.

ix_branch : int
    Branch index to track through.

print_err : bool, optional
    Default is False. Print error messages if, eg, lattice is unstable?

force_calc : bool, optional
    Default is False. If True, force the calculation to be done.

Returns
-------
calc_ok : bool
    Set True if there were no problems, False otherwise.
)"""
  );
  m.def(
      "tao_lattice_branches_equal_tao_lattice_branches",
      &Tao::tao_lattice_branches_equal_tao_lattice_branches,
      py::arg("tlb1"),
      py::arg("tlb2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lattice_branches_equal_tao_lattice_branches

Parameters
----------
tlb1 : 1D array of TaoLatticeBranchStruct

tlb2 : 1D array of TaoLatticeBranchStruct
)"""
  );
  py::class_<Tao::TaoLatticeCalc, std::unique_ptr<Tao::TaoLatticeCalc>>(
      m,
      "TaoLatticeCalc",
      "tao_lattice_calc return type"
  )
      .def_readonly("calc_ok", &Tao::TaoLatticeCalc::calc_ok)
      .def_readonly("print_err", &Tao::TaoLatticeCalc::print_err)
      .def("__len__", [](const Tao::TaoLatticeCalc &) { return 2; })
      .def("__getitem__", [](const Tao::TaoLatticeCalc &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.calc_ok);
        if (i == 1)
          return py::cast(s.print_err);
        throw py::index_error();
      });
  m.def(
      "tao_lattice_calc",
      &Tao::tao_lattice_calc,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lattice_calc

Returns
-------
calc_ok : bool
    Set False if there was an error in the calculation like a particle was lost or a lat is unstable.

print_err : bool, optional
    Default True. If False, do not print error messages if, for example, the lattice is unstable.
)"""
  );
  m.def(
      "tao_lattice_equal_tao_lattice",
      &Tao::tao_lattice_equal_tao_lattice,
      py::arg("lat1"),
      py::arg("lat2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lattice_equal_tao_lattice

Parameters
----------
lat1 : TaoLatticeStruct

lat2 : TaoLatticeStruct
)"""
  );
  m.def(
      "tao_limit_calc",
      &Tao::tao_limit_calc,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_limit_calc

Returns
-------
limited : bool
    Set True if a variable is past a limit.
)"""
  );
  m.def(
      "tao_lm_optimizer",
      &Tao::tao_lm_optimizer,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_lm_optimizer (abort)

Routine to minimize the merit function by varying variables until
the "data" as calculated from the model matches the measured data.

This subroutine is a wrapper for the mrqmin routine of Numerical Recipes.
See the Numerical Recipes writeup for more details.
'lm' stands for Levenburg - Marquardt. Otherwise known as LMDIF.

Returns
-------
abort : bool
    Set True if an user stop signal detected.
)"""
  );
  m.def(
      "tao_lmdif_optimizer",
      &Tao::tao_lmdif_optimizer,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_lmdif_optimizer

Returns
-------
abort : bool
    Set True if an user stop signal detected or there is a problem with calculating the merit function.
)"""
  );
  m.def(
      "tao_load_this_datum",
      &Tao::tao_load_this_datum,
      py::arg("vec"),
      py::arg("ele_ref"),
      py::arg("ele_start"),
      py::arg("ele"),
      py::arg("datum_value"),
      py::arg("valid_value"),
      py::arg("datum"),
      py::arg("branch"),
      py::arg("why_invalid") = py::none(),
      py::arg("good") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_load_this_datum

Parameters
----------
vec : 1D array of float

ele_ref : EleStruct

ele_start : EleStruct

ele : EleStruct

datum_value : float

valid_value : bool

datum : TaoDataStruct

branch : BranchStruct

why_invalid : str, optional

good : 1D array of bool, optional
)"""
  );
  py::class_<Tao::TaoLocateAllElements, std::unique_ptr<Tao::TaoLocateAllElements>>(
      m,
      "TaoLocateAllElements",
      "tao_locate_all_elements return type"
  )
      .def_readonly("eles", &Tao::TaoLocateAllElements::eles)
      .def_readonly("err", &Tao::TaoLocateAllElements::err)
      .def("__len__", [](const Tao::TaoLocateAllElements &) { return 2; })
      .def("__getitem__", [](const Tao::TaoLocateAllElements &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.eles);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "tao_locate_all_elements",
      &Tao::tao_locate_all_elements,
      py::arg("ele_list"),
      py::arg("ignore_blank") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_locate_all_elements

Parameters
----------
ele_list : str
    String with element names using element list format.

ignore_blank : bool, optional
    If present and true then do nothing if ele_list is blank. otherwise a blank is treated as an error.

Returns
-------
eles : 1D array of ElePointerStruct
    : Array of elements in the model lat.

err : bool
    Set true on error.
)"""
  );
  py::class_<Tao::TaoLocateElements, std::unique_ptr<Tao::TaoLocateElements>>(
      m,
      "TaoLocateElements",
      "tao_locate_elements return type"
  )
      .def_readonly("eles", &Tao::TaoLocateElements::eles)
      .def_readonly("err", &Tao::TaoLocateElements::err)
      .def("__len__", [](const Tao::TaoLocateElements &) { return 2; })
      .def("__getitem__", [](const Tao::TaoLocateElements &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.eles);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "tao_locate_elements",
      &Tao::tao_locate_elements,
      py::arg("ele_list"),
      py::arg("ix_universe"),
      py::arg("lat_type") = py::none(),
      py::arg("ignore_blank") = py::none(),
      py::arg("err_stat_level") = py::none(),
      py::arg("above_ubound_is_err") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("multiple_eles_is_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_locate_elements

Parameters
----------
ele_list : str
    String with element names using element list format.

ix_universe : int
    Universe to search. -1 => search s.global.default_universe. -2 (all unis) => error. ix_universe is ignored
    if ele_list starts with a universe specifier "N@".

lat_type : int, optional
    model$ (default), design$, or base$.

ignore_blank : bool, optional
    If present and true then do nothing if ele_list is blank. otherwise treated as an error.

err_stat_level : int, optional
    Status level for error messages. If not present, print with level s_error$. Use s_nooutput$ to prevent
    printing.

above_ubound_is_err : bool, optional

ix_branch : int, optional
    If present and non-negative then use this as the branch index for elements specified using an integer
    index (EG: "43"). If -1 use the default branch, search all branches.

multiple_eles_is_err : bool, optional
    If present and True then matching to more than one element is an error.

Returns
-------
eles : 1D array of ElePointerStruct
    : Array of elements in the model lat.

err : bool
    Set true on error.
)"""
  );
  m.def(
      "tao_mark_lattice_ele",
      &Tao::tao_mark_lattice_ele,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_mark_lattice_ele

Parameters
----------
lat : LatStruct
    Input lattice
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with elements marked.
)"""
  );
  py::class_<Tao::TaoMerit, std::unique_ptr<Tao::TaoMerit>>(m, "TaoMerit", "tao_merit return type")
      .def_readonly("calc_ok", &Tao::TaoMerit::calc_ok)
      .def_readonly("this_merit", &Tao::TaoMerit::this_merit)
      .def("__len__", [](const Tao::TaoMerit &) { return 2; })
      .def("__getitem__", [](const Tao::TaoMerit &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.calc_ok);
        if (i == 1)
          return py::cast(s.this_merit);
        throw py::index_error();
      });
  m.def(
      "tao_merit",
      &Tao::tao_merit,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_merit

Returns
-------
this_merit : float
    Merit value.

calc_ok : bool, optional
    Set False if there was an error in the calculation like a particle was lost or a lat is unstable.
)"""
  );
  py::class_<PyTaoNextSwitch, std::unique_ptr<PyTaoNextSwitch>>(
      m,
      "TaoNextSwitch",
      "tao_next_switch return type"
  )
      .def_readonly("switch_", &PyTaoNextSwitch::switch_)
      .def_readonly("err", &PyTaoNextSwitch::err)
      .def_readonly("line", &PyTaoNextSwitch::line)
      .def("__len__", [](const PyTaoNextSwitch &) { return 3; })
      .def("__getitem__", [](const PyTaoNextSwitch &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.switch_);
        if (i == 1)
          return py::cast(s.err);
        if (i == 2)
          return py::cast(s.line);
        throw py::index_error();
      });
  m.def(
      "tao_next_switch",
      &python_tao_next_switch,
      py::arg("line"),
      py::arg("switch_list"),
      py::arg("return_next_word"),
      py::arg("neg_num_not_switch") = py::none(),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_next_switch (line, switch_list, return_next_word, switch, err, neg_num_not_switch, print_err)

Subroutine look at the next word on the command line and match this word to a list of "switches"
given by the switch_list argument.

If switch_list(1) starts with a "-" or "#" character, switches are assumed to start with this character.
If switch_list(1) starts with any other character, everything is considered to be a switch.

Switch abbreviations are permitted.

If return_next_word = True then, when a non-switch word is encountered, the switch argument
will be set to that word and that word will be removed from the line argument.

If return_next_word = False then, when a non-switch word is encountered, the switch argument
will be set to '' and the non-switch word will be left on the line argument.

If the first non-blank character in line is a single or double quote. The word returned will be the
substring from the initial quote mark to the next matching quote mark. The quote marks will be removed
from the returned switch argument.

Parameters
----------
line : str
    Command line
    This parameter is an input/output and is modified in-place.
    As an output, line: Command line with first word removed if

switch_list : 1D array of str
    List of valid switches.

return_next_word : bool
    See above.

neg_num_not_switch : bool, optional
    If present and True then a word like "-34" will be treated as a non-switch.

print_err : bool, optional
    Default is True. If False, do not print unknown switch error.

Returns
-------
line : str
    Command line
    This parameter is an input/output and is modified in-place.
    As an output, line: Command line with first word removed if

err : bool
    Set True if the next word begins with '-' but there is no match to anything in switch_list.
)"""
  );
  py::class_<PyTaoNextWord, std::unique_ptr<PyTaoNextWord>>(
      m,
      "TaoNextWord",
      "tao_next_word return type"
  )
      .def_readonly("word", &PyTaoNextWord::word)
      .def_readonly("line", &PyTaoNextWord::line)
      .def("__len__", [](const PyTaoNextWord &) { return 2; })
      .def("__getitem__", [](const PyTaoNextWord &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.word);
        if (i == 1)
          return py::cast(s.line);
        throw py::index_error();
      });
  m.def(
      "tao_next_word",
      &python_tao_next_word,
      py::arg("line"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_next_word (line, word)

Routine to return the next word in a line.

Words are delimited by a space character except if the space is within quotes.
Additionally, spaces within brackets "(...)", "{...}", and "[...]" are ignored.
Outer quote marks will be removed in the returned word.

Parameters
----------
line : str
    String to parse.
    This parameter is an input/output and is modified in-place.
    As an output, line: String with first word removed.

Returns
-------
line : str
    String to parse.
    This parameter is an input/output and is modified in-place.
    As an output, line: String with first word removed.

word : str
    First word of line.
)"""
  );
  m.def(
      "tao_one_turn_map_calc_needed",
      &Tao::tao_one_turn_map_calc_needed,
      py::arg("data_type"),
      py::arg("data_source"),
      py::arg("do_one_turn_map"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_one_turn_map_calc_needed

Parameters
----------
data_type : str

data_source : str

do_one_turn_map : bool
)"""
  );
  m.def(
      "tao_open_file",
      &Tao::tao_open_file,
      py::arg("file"),
      py::arg("file_name"),
      py::arg("error_severity"),
      py::arg("binary") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_open_file

Parameters
----------
file : str

file_name : str
    File name.

error_severity : int
    Severity level used in the error message. Possibilities are s_fatal$, etc. See out_io doc for more
    details. Use -1 to not print a message if file cannot be opened.

binary : bool, optional
    If present and True then open a binary file, Defaut is False.

Returns
-------
iunit : int
    Logical unit number. Set to 0 if file not openable.
)"""
  );
  py::class_<Tao::TaoOpenScratchFile, std::unique_ptr<Tao::TaoOpenScratchFile>>(
      m,
      "TaoOpenScratchFile",
      "tao_open_scratch_file return type"
  )
      .def_readonly("err", &Tao::TaoOpenScratchFile::err)
      .def_readonly("iu", &Tao::TaoOpenScratchFile::iu)
      .def("__len__", [](const Tao::TaoOpenScratchFile &) { return 2; })
      .def("__getitem__", [](const Tao::TaoOpenScratchFile &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err);
        if (i == 1)
          return py::cast(s.iu);
        throw py::index_error();
      });
  m.def(
      "tao_open_scratch_file",
      &Tao::tao_open_scratch_file,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_open_scratch_file

Returns
-------
err : bool
    Set True if there is an error. False otherwise.

iu : int
    File handle unit number.
)"""
  );
  m.def(
      "tao_optimization_status",
      &Tao::tao_optimization_status,
      py::arg("datum"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_optimization_status

Parameters
----------
datum : TaoDataStruct
    Datum to evaluate.

Returns
-------
why_str : str
    Optimization status of the datum.
)"""
  );
  m.def(
      "tao_orbit_beta_wave_anal",
      &Tao::tao_orbit_beta_wave_anal,
      py::arg("plot"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_orbit_beta_wave_anal (plot)
)"""
  );
  m.def(
      "tao_oreint_building_wall_pt",
      &Tao::tao_oreint_building_wall_pt,
      py::arg("pt_in"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_oreint_building_wall_pt

Parameters
----------
pt_in : TaoBuildingWallPointStruct
    Building wall point.

Returns
-------
pt_out : TaoBuildingWallPointStruct
    Building wall point with orientation params applied.
)"""
  );
  py::class_<Tao::TaoParamValueAtS, std::unique_ptr<Tao::TaoParamValueAtS>>(
      m,
      "TaoParamValueAtS",
      "tao_param_value_at_s return type"
  )
      .def_readonly("err_flag", &Tao::TaoParamValueAtS::err_flag)
      .def_readonly("why_invalid", &Tao::TaoParamValueAtS::why_invalid)
      .def_readonly("print_err", &Tao::TaoParamValueAtS::print_err)
      .def_readonly("bad_datum", &Tao::TaoParamValueAtS::bad_datum)
      .def_readonly("value", &Tao::TaoParamValueAtS::value)
      .def("__len__", [](const Tao::TaoParamValueAtS &) { return 5; })
      .def("__getitem__", [](const Tao::TaoParamValueAtS &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.why_invalid);
        if (i == 2)
          return py::cast(s.print_err);
        if (i == 3)
          return py::cast(s.bad_datum);
        if (i == 4)
          return py::cast(s.value);
        throw py::index_error();
      });
  m.def(
      "tao_param_value_at_s",
      &Tao::tao_param_value_at_s,
      py::arg("dat_name"),
      py::arg("ele_to_s"),
      py::arg("ele_here"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_param_value_at_s

Parameters
----------
dat_name : str

ele_to_s : EleStruct
    Element whose exit end is at the evaluation s-position.

ele_here : EleStruct
    Lattice element that overlaps the s-position ele.s.

orbit : CoordStruct
    Orbit at the evaluation s-position.

Returns
-------
err_flag : bool
    Set true if parameter cannot be evaluated.

value : float
    Parameter value.

why_invalid : str, optional
    Set if  err_flag = True to document why is there a problem.

print_err : bool, optional
    Print error message on error? Default is True.

bad_datum : bool, optional
    Data_type is malformed.
)"""
  );
  m.def(
      "tao_param_value_routine",
      &Tao::tao_param_value_routine,
      py::arg("str"),
      py::arg("use_good_user"),
      py::arg("saved_prefix"),
      py::arg("stack"),
      py::arg("err_flag"),
      py::arg("print_err"),
      py::arg("dflt_component") = py::none(),
      py::arg("dflt_source") = py::none(),
      py::arg("dflt_ele_ref") = py::none(),
      py::arg("dflt_ele_start") = py::none(),
      py::arg("dflt_ele") = py::none(),
      py::arg("dflt_dat_or_var_index") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::arg("dflt_eval_point") = py::none(),
      py::arg("dflt_s_offset") = py::none(),
      py::arg("dflt_orbit") = py::none(),
      py::arg("datum") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_param_value_routine

Parameters
----------
str : str

use_good_user : bool

saved_prefix : str

stack : TaoEvalNodeStruct

err_flag : bool

print_err : bool

dflt_component : str, optional

dflt_source : str, optional

dflt_ele_ref : EleStruct, optional

dflt_ele_start : EleStruct, optional

dflt_ele : EleStruct, optional

dflt_dat_or_var_index : str, optional

dflt_uni : int, optional

dflt_eval_point : int, optional

dflt_s_offset : float, optional

dflt_orbit : CoordStruct, optional

datum : TaoDataStruct, optional
)"""
  );
  m.def(
      "tao_parse_command_args",
      &Tao::tao_parse_command_args,
      py::arg("cmd_line") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_parse_command_args

Parameters
----------
cmd_line : str, optional

Returns
-------
error : bool
    Set True if there is an error. False otherwise.
)"""
  );
  py::class_<Tao::TaoParseElementParamStr, std::unique_ptr<Tao::TaoParseElementParamStr>>(
      m,
      "TaoParseElementParamStr",
      "tao_parse_element_param_str return type"
  )
      .def_readonly("err", &Tao::TaoParseElementParamStr::err)
      .def_readonly("uni", &Tao::TaoParseElementParamStr::uni)
      .def_readonly("element", &Tao::TaoParseElementParamStr::element)
      .def_readonly("parameter", &Tao::TaoParseElementParamStr::parameter)
      .def_readonly("where", &Tao::TaoParseElementParamStr::where)
      .def_readonly("component", &Tao::TaoParseElementParamStr::component)
      .def("__len__", [](const Tao::TaoParseElementParamStr &) { return 6; })
      .def("__getitem__", [](const Tao::TaoParseElementParamStr &s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.err);
        if (i == 1)
          return py::cast(s.uni);
        if (i == 2)
          return py::cast(s.element);
        if (i == 3)
          return py::cast(s.parameter);
        if (i == 4)
          return py::cast(s.where);
        if (i == 5)
          return py::cast(s.component);
        throw py::index_error();
      });
  m.def(
      "tao_parse_element_param_str",
      &Tao::tao_parse_element_param_str,
      py::arg("in_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_parse_element_param_str

Parameters
----------
in_str : str
    String specifying a parameter of an element or elements.

Returns
-------
err : bool
    Set True if there is a parse error. False otherwise.

uni : str
    Universe substring.

element : str
    Element name.

parameter : str
    Element parameter name.

where : int
    One of not_set$, anchor_beginning$, anchor_center$, or anchor_end$.

component : str
    One of "model", "design", or "base".
)"""
  );
  py::class_<Tao::TaoParticleDataValue, std::unique_ptr<Tao::TaoParticleDataValue>>(
      m,
      "TaoParticleDataValue",
      "tao_particle_data_value return type"
  )
      .def_readonly("value", &Tao::TaoParticleDataValue::value)
      .def_readonly("err", &Tao::TaoParticleDataValue::err)
      .def("__len__", [](const Tao::TaoParticleDataValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoParticleDataValue &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "tao_particle_data_value",
      &Tao::tao_particle_data_value,
      py::arg("data_type"),
      py::arg("p"),
      py::arg("ele"),
      py::arg("ix_bunch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_particle_data_value (data_type, p, value, err, ele, ix_bunch)

Routine to calculate the value array of a data_type for an array of particles.

Parameters
----------
data_type : str
    Type of data.

p : 1D array of CoordStruct
    coord_struct, Array of particles containing the data.

ele : EleStruct
    Needed for "Ja" evaluation.

ix_bunch : int
    Bunch index.

Returns
-------
value : 1D array of float
    Array of values.

err : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "tao_pause_cmd",
      &Tao::tao_pause_cmd,
      py::arg("time"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pause_cmd

Parameters
----------
time : float
    Time to pause in seconds.
)"""
  );
  m.def(
      "tao_phase_space_axis_index",
      &Tao::tao_phase_space_axis_index,
      py::arg("data_type"),
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_phase_space_axis_index (data_type, err) result (ix_axis)

Routine to calculate the phase space axis index for a given data type.

Parameters
----------
data_type : str
    Type of data.

err : bool
    Set True if there is an error.

Returns
-------
ix_axis : int
    Axis index.
)"""
  );
  m.def(
      "tao_phase_wave_anal",
      &Tao::tao_phase_wave_anal,
      py::arg("plot"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_phase_wave_anal (plot)
)"""
  );
  py::class_<Tao::TaoPickUniverse, std::unique_ptr<Tao::TaoPickUniverse>>(
      m,
      "TaoPickUniverse",
      "tao_pick_universe return type"
  )
      .def_readonly("name_out", &Tao::TaoPickUniverse::name_out)
      .def_readonly("picked", &Tao::TaoPickUniverse::picked)
      .def_readonly("err", &Tao::TaoPickUniverse::err)
      .def_readonly("ix_uni", &Tao::TaoPickUniverse::ix_uni)
      .def_readonly("explicit_uni", &Tao::TaoPickUniverse::explicit_uni)
      .def("__len__", [](const Tao::TaoPickUniverse &) { return 5; })
      .def("__getitem__", [](const Tao::TaoPickUniverse &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.name_out);
        if (i == 1)
          return py::cast(s.picked);
        if (i == 2)
          return py::cast(s.err);
        if (i == 3)
          return py::cast(s.ix_uni);
        if (i == 4)
          return py::cast(s.explicit_uni);
        throw py::index_error();
      });
  m.def(
      "tao_pick_universe",
      &Tao::tao_pick_universe,
      py::arg("name_in"),
      py::arg("dflt_uni") = py::none(),
      py::arg("pure_uni") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pick_universe

Parameters
----------
name_in : str
    data name with possible universe spec.

dflt_uni : int, optional
    Default universe to use. Set to -1 if explicit universe is required.

pure_uni : bool, optional
    Default is False. See above

Returns
-------
name_out : str
    name_in without any "n@" beginning.

picked : 1D array of bool
    Array showing picked universes. The array will be resized if necessary.

err : bool
    Set True if an error is detected.

ix_uni : int, optional
    Set to the picked universe with the highest index.

explicit_uni : bool, optional
    Set True if name_in has explicit universe "n@" specification.
)"""
  );
  m.def(
      "tao_pipe_cmd",
      &Tao::tao_pipe_cmd,
      py::arg("input_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pipe_cmd

Parameters
----------
input_str : str
    What to show.
)"""
  );
  m.def(
      "tao_place_cmd",
      &Tao::tao_place_cmd,
      py::arg("where"),
      py::arg("who"),
      py::arg("no_buffer") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_place_cmd

Parameters
----------
where : str
    Region where the plot goes. Eg: 'top'.

who : str
    Type of plot. Eg: 'orbit'.

no_buffer : bool, optional
    If present and True then prevents buffering in the case when s.global.external_plotting = T
)"""
  );
  m.def(
      "tao_plot_cmd",
      &Tao::tao_plot_cmd,
      py::arg("where"),
      py::arg("component"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_plot_cmd

Parameters
----------
where : str
    Region name to identify the plot to set.

component : str
    Who to plot. EG: 'meas - design'
)"""
  );
  m.def(
      "tao_plot_data",
      &Tao::tao_plot_data,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_plot_data (plot, graph)

Routine to draw a graph with data and/or variable curves.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  m.def(
      "tao_plot_histogram",
      &Tao::tao_plot_histogram,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_plot_histogram (plot, graph)

Routine to draw one graph for the histogram analysis plot.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  m.def(
      "tao_plot_key_table",
      &Tao::tao_plot_key_table,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_plot_key_table (plot, graph)

Routine to draw a key table graph.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  m.def(
      "tao_plot_setup",
      &Tao::tao_plot_setup,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_plot_setup
)"""
  );
  m.def(
      "tao_plot_struct_transfer",
      &Tao::tao_plot_struct_transfer,
      py::arg("plot_in"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_plot_struct_transfer

Parameters
----------
plot_in : TaoPlotStruct
    Input structure.

Returns
-------
plot_out : TaoPlotStruct
    Output struture.
)"""
  );
  m.def(
      "tao_plot_wave",
      &Tao::tao_plot_wave,
      py::arg("plot"),
      py::arg("graph"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_plot_wave (plot, graph)

Routine to draw one graph for the wave analysis plot.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  m.def(
      "tao_pointer_to_building_wall_shape",
      &Tao::tao_pointer_to_building_wall_shape,
      py::arg("wall_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pointer_to_building_wall_shape

Parameters
----------
wall_name : str
    Name of the wall.

Returns
-------
e_shape : TaoEleShapeStruct, optional
    Associated shape. Nullified if there is no associated shape.
)"""
  );
  m.def(
      "tao_pointer_to_datum",
      &Tao::tao_pointer_to_datum,
      py::arg("d1"),
      py::arg("ele_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pointer_to_datum

Parameters
----------
d1 : TaoD1DataStruct
    D1 data struct to search.

ele_name : str
    Name of lattice element to match to.

Returns
-------
datum_ptr : TaoDataStruct, optional
    Pointer to the matched datum. Will be null if no match found.
)"""
  );
  py::class_<Tao::TaoPointerToDatumEle, std::unique_ptr<Tao::TaoPointerToDatumEle>>(
      m,
      "TaoPointerToDatumEle",
      "tao_pointer_to_datum_ele return type"
  )
      .def_readonly("valid", &Tao::TaoPointerToDatumEle::valid)
      .def_readonly("why_invalid", &Tao::TaoPointerToDatumEle::why_invalid)
      .def_readonly("ele", &Tao::TaoPointerToDatumEle::ele)
      .def("__len__", [](const Tao::TaoPointerToDatumEle &) { return 3; })
      .def("__getitem__", [](const Tao::TaoPointerToDatumEle &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.valid);
        if (i == 1)
          return py::cast(s.why_invalid);
        if (i == 2)
          return py::cast(s.ele);
        throw py::index_error();
      });
  m.def(
      "tao_pointer_to_datum_ele",
      &Tao::tao_pointer_to_datum_ele,
      py::arg("lat"),
      py::arg("ele_name"),
      py::arg("ix_ele"),
      py::arg("datum"),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_pointer_to_datum_ele (lat, ele_name, ix_ele, datum, valid, why_invalid, print_err) result (ele)

Routine to see if an element index corresponds to an element with a definite
location such as an overlay or multipass element.

If the element is a super_lord then the super_slave element at the exit end
of the lord will be returned. Otherwise ix_loc will be set to ix_ele.

Parameters
----------
lat : LatStruct
    Lattice

ix_ele : int
    Index of element.

datum : TaoDataStruct
    Used for error messages and gives branch index.

print_err : bool, optional
    Default is True. If False, do not print an error message.

Returns
-------
valid : bool
    Set False if element does not have a definite location. Set True otherwise

why_invalid : str, optional
    Tells why datum value is invalid.

ele : EleStruct, optional
    : Pointer to the element. Set to NULL if not valid or no associated element.
)"""
  );
  py::class_<PyTaoPointerToEleShape, std::unique_ptr<PyTaoPointerToEleShape>>(
      m,
      "TaoPointerToEleShape",
      "tao_pointer_to_ele_shape return type"
  )
      .def_readonly("dat_var_name", &PyTaoPointerToEleShape::dat_var_name)
      .def_readonly("dat_var_value", &PyTaoPointerToEleShape::dat_var_value)
      .def_readonly("e_shape", &PyTaoPointerToEleShape::e_shape)
      .def_readonly("ix_shape_min", &PyTaoPointerToEleShape::ix_shape_min)
      .def("__len__", [](const PyTaoPointerToEleShape &) { return 4; })
      .def("__getitem__", [](const PyTaoPointerToEleShape &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.dat_var_name);
        if (i == 1)
          return py::cast(s.dat_var_value);
        if (i == 2)
          return py::cast(s.e_shape);
        if (i == 3)
          return py::cast(s.ix_shape_min);
        throw py::index_error();
      });
  m.def(
      "tao_pointer_to_ele_shape",
      &python_tao_pointer_to_ele_shape,
      py::arg("ix_uni"),
      py::arg("ele"),
      py::arg("ele_shape"),
      py::arg("ix_shape_min") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pointer_to_ele_shape

Parameters
----------
ix_uni : int
    Universe index.

ele : EleStruct
    Lattice element.

ele_shape : 1D array of TaoEleShapeStruct
    Array of shapes to search.

ix_shape_min : int, optional
    Index of minimum ele_shape(:) index to start search from. Default is 1.
    This parameter is an input/output and is modified in-place.
    As an output, ix_shape_min: Ele_shape(

Returns
-------
dat_var_name : str, optional
    Name of datum or variable associated with e_shape. Will be set to "" if there is no associated datum or
    variable.

dat_var_value : float, optional
    Value of datum or variable associated with e_shape. Will be set to zero if there is no associated datum or
    variable.

ix_shape_min : int, optional
    Index of minimum ele_shape(:) index to start search from. Default is 1.
    This parameter is an input/output and is modified in-place.
    As an output, ix_shape_min: Ele_shape(

e_shape : TaoEleShapeStruct, optional
    Associated shape. Nullified if there is no associated shape.
)"""
  );
  m.def(
      "tao_pointer_to_tao_lat",
      &Tao::tao_pointer_to_tao_lat,
      py::arg("u"),
      py::arg("lat_type") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pointer_to_tao_lat

Parameters
----------
u : TaoUniverseStruct
    Universe to work with

lat_type : int, optional
    model$ (default), design$, or base$.

Returns
-------
tao_lat : TaoLatticeStruct, optional
    Tao_lat pointer. Points to u.model, u.design, or u.base
)"""
  );
  m.def(
      "tao_pointer_to_universe",
      py::overload_cast<int, std::optional<bool>>(&Tao::tao_pointer_to_universe),
      py::arg("ix_uni"),
      py::arg("neg2_to_default") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_pointer_to_universe (...) result (u)

Routine to set a pointer to a universe.

This is an overloaded routine for the:
 tao_pointer_to_universe_int (ix_uni, neg2_to_default) result (u)
 tao_pointer_to_universe_str (string, neg2_to_default) result (u)

Note: With a string argument, this routine can only handle single universe picks.
That is, it cannot handlle something like "[1,3,4]@...". To handle multiple universe picks, use:
  tao_pointer_to_universes

Parameters
----------
ix_uni : int
    Index to the s.u(:) array If ix_uni is -1 -> u(s.global.default_universe) will be used.

neg2_to_default : bool, optional
    i_uni = -2 (all universes) maps to the default uni? Default if False.

Returns
-------
u : TaoUniverseStruct, optional
    Universe pointer. u will be nullified if there is an error and an error message will be printed.
)"""
  );
  py::class_<PyTaoPointerToUniverseStr, std::unique_ptr<PyTaoPointerToUniverseStr>>(
      m,
      "TaoPointerToUniverseStr",
      "tao_pointer_to_universe_str return type"
  )
      .def_readonly("u", &PyTaoPointerToUniverseStr::u)
      .def_readonly("string", &PyTaoPointerToUniverseStr::string)
      .def("__len__", [](const PyTaoPointerToUniverseStr &) { return 2; })
      .def("__getitem__", [](const PyTaoPointerToUniverseStr &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.u);
        if (i == 1)
          return py::cast(s.string);
        throw py::index_error();
      });
  m.def(
      "tao_pointer_to_universe",
      py::overload_cast<std::string, std::optional<bool>>(&python_tao_pointer_to_universe_str),
      py::arg("string"),
      py::arg("neg2_to_default") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_pointer_to_universe (...) result (u)

Routine to set a pointer to a universe.

This is an overloaded routine for the:
 tao_pointer_to_universe_int (ix_uni, neg2_to_default) result (u)
 tao_pointer_to_universe_str (string, neg2_to_default) result (u)

Note: With a string argument, this routine can only handle single universe picks.
That is, it cannot handlle something like "[1,3,4]@...". To handle multiple universe picks, use:
  tao_pointer_to_universes

Parameters
----------
string : str
    String in the form "<ix_uni>@..." or, if no "@" is present, u will point to the default universe.
    This parameter is an input/output and is modified in-place.
    As an output, string: String with universe prefix stripped off.

neg2_to_default : bool, optional
    i_uni = -2 (all universes) maps to the default uni? Default if False.

Returns
-------
string : str
    String in the form "<ix_uni>@..." or, if no "@" is present, u will point to the default universe.
    This parameter is an input/output and is modified in-place.
    As an output, string: String with universe prefix stripped off.

u : TaoUniverseStruct, optional
    Universe pointer. u will be nullified if there is an error and an error message will be printed.
)"""
  );
  py::class_<Tao::TaoPointerToUniverses, std::unique_ptr<Tao::TaoPointerToUniverses>>(
      m,
      "TaoPointerToUniverses",
      "tao_pointer_to_universes return type"
  )
      .def_readonly("unis", &Tao::TaoPointerToUniverses::unis)
      .def_readonly("err", &Tao::TaoPointerToUniverses::err)
      .def_readonly("name_out", &Tao::TaoPointerToUniverses::name_out)
      .def_readonly("explicit_uni", &Tao::TaoPointerToUniverses::explicit_uni)
      .def("__len__", [](const Tao::TaoPointerToUniverses &) { return 4; })
      .def("__getitem__", [](const Tao::TaoPointerToUniverses &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.unis);
        if (i == 1)
          return py::cast(s.err);
        if (i == 2)
          return py::cast(s.name_out);
        if (i == 3)
          return py::cast(s.explicit_uni);
        throw py::index_error();
      });
  m.def(
      "tao_pointer_to_universes",
      &Tao::tao_pointer_to_universes,
      py::arg("name_in"),
      py::arg("dflt_uni") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_pointer_to_universes

Parameters
----------
name_in : str
    data name with possible universe spec.

dflt_uni : int, optional
    Default universe to use. Set to -1 if explicit universe is required.

Returns
-------
unis : 1D array of TaoUniversePointerStruct
    Array of pointers to picked universes. The array will be resized if necessary.

err : bool
    Set True if an error is detected.

name_out : str, optional
    name_in without any "n@" beginning.

explicit_uni : bool, optional
    Set True if name_in has explicit universe "n@" specification.
)"""
  );
  m.def(
      "tao_pointer_to_var_in_lattice",
      &Tao::tao_pointer_to_var_in_lattice,
      py::arg("var"),
      py::arg("ix_uni"),
      py::arg("ele"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_pointer_to_var_in_lattice (var, ix_uni, ele, err)

Routine to add a pointer from a given Tao variable
to the appropriate variable in a lattice.

Parameters
----------
var : TaoVarStruct
    Structure has the info of where to point.

ix_uni : int
    the universe to use

Returns
-------
err : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "tao_pointer_to_var_in_lattice2",
      &Tao::tao_pointer_to_var_in_lattice2,
      py::arg("var"),
      py::arg("ix_uni"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_pointer_to_var_in_lattice2 (var, ix_uni, err)

Routine to add a pointer from a given Tao variable
to the appropriate variable in a lattice.

Parameters
----------
var : TaoVarStruct
    Structure has the info of where to point.

ix_uni : int
    the universe to use

Returns
-------
err : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "tao_print_command_line_info",
      &Tao::tao_print_command_line_info,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_print_command_line_info
)"""
  );
  m.def(
      "tao_ptc_normal_form",
      &Tao::tao_ptc_normal_form,
      py::arg("do_calc"),
      py::arg("tao_lat"),
      py::arg("ix_branch"),
      py::arg("rf_on") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_ptc_normal_form

Parameters
----------
do_calc : bool
    Set True to do the calculation.

tao_lat : TaoLatticeStruct
    Lattice to work on.

ix_branch : int
    Branch of lattice to work on.

rf_on : int, optional
    RF state for calculation. yes$, no$, or maybe$ (default) maybe$ means that RF state in branch is used.
)"""
  );
  m.def(
      "tao_python_cmd",
      &Tao::tao_python_cmd,
      py::arg("input_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_python_cmd

Parameters
----------
input_str : str
    What to show.
)"""
  );
  m.def(
      "tao_quiet_set",
      &Tao::tao_quiet_set,
      py::arg("set"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_quiet_set

Parameters
----------
set : str
    True is silent running is wanted.
)"""
  );
  m.def(
      "tao_rad_int_calc_needed",
      &Tao::tao_rad_int_calc_needed,
      py::arg("data_type"),
      py::arg("data_source"),
      py::arg("do_rad_int"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_rad_int_calc_needed

Parameters
----------
data_type : str

data_source : str

do_rad_int : bool
)"""
  );
  m.def(
      "tao_re_allocate_expression_info",
      &Tao::tao_re_allocate_expression_info,
      py::arg("info"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_re_allocate_expression_info

Parameters
----------
info : 1D array of TaoExpressionInfoStruct
    This parameter is an input/output and is modified in-place.
    As an output, info: Allocated array with size(re) >= n.

n : int
    Size wanted.

exact : bool, optional
    If present and False then the size of the output array is permitted to be larger than n. Default is True.
)"""
  );
  m.def(
      "tao_re_associate_node_array",
      &Tao::tao_re_associate_node_array,
      py::arg("tree"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_re_associate_node_array(tree, n, exact)

Routine to resize the tree%node(:) array.

Note: The data of the array is preserved but data at the end of the
array will be lost if n is less than the original size of the array

Parameters
----------
tree : TaoEvalNodeStruct

n : int
    Size wanted.

exact : bool, optional
    Default is False. If False, the size of the output array is permitted to be larger than n.
)"""
  );
  m.def(
      "tao_re_execute",
      &Tao::tao_re_execute,
      py::arg("string"),
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_re_exectue (string, err)

Subroutine to execute a previous command.
)"""
  );
  m.def(
      "tao_read_cmd",
      &Tao::tao_read_cmd,
      py::arg("which"),
      py::arg("unis"),
      py::arg("file"),
      py::arg("silent"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_read_cmd

Parameters
----------
which : str

unis : str
    Universes to apply to

file : str

silent : bool
    Silent
)"""
  );
  m.def(
      "tao_read_phase_space_index",
      &Tao::tao_read_phase_space_index,
      py::arg("name"),
      py::arg("ixc"),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_read_phase_space_index

Parameters
----------
name : str
    character array holding the index. Must be in the range 1-6.

ixc : int
    location within <name> to evaluate index.

print_err : bool, optional
    If present and False then do not print an error message

Returns
-------
ix_ps : int
    Index at <name>(<ixc>:<ixc>). Returns 0 if bad index.
)"""
  );
  m.def(
      "tao_regression_test",
      &Tao::tao_regression_test,
      py::arg("cmd_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_regression_test

Parameters
----------
cmd_str : str
)"""
  );
  py::class_<PyTaoRemoveBlankCharacters, std::unique_ptr<PyTaoRemoveBlankCharacters>>(
      m,
      "TaoRemoveBlankCharacters",
      "tao_remove_blank_characters return type"
  )
      .def_readonly("str", &PyTaoRemoveBlankCharacters::str)
      .def("__len__", [](const PyTaoRemoveBlankCharacters &) { return 1; })
      .def("__getitem__", [](const PyTaoRemoveBlankCharacters &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "tao_remove_blank_characters",
      &python_tao_remove_blank_characters,
      py::arg("str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_remove_blank_characters

Parameters
----------
str : str
    Input string.
    This parameter is an input/output and is modified in-place.
    As an output, str: String with blank characters removed.

Returns
-------
str : str
    Input string.
    This parameter is an input/output and is modified in-place.
    As an output, str: String with blank characters removed.
)"""
  );
  m.def(
      "tao_run_cmd",
      &Tao::tao_run_cmd,
      py::arg("which"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_run_cmd

Parameters
----------
which : str
    which optimizer to use.

Returns
-------
abort : bool
    Set True if the run was aborted by the user, an at minimum condition, a singular matrix condition, etc..
    False otherwise.
)"""
  );
  m.def(
      "tao_scale_cmd",
      &Tao::tao_scale_cmd,
      py::arg("where"),
      py::arg("y_min_in"),
      py::arg("y_max_in"),
      py::arg("axis") = py::none(),
      py::arg("include_wall") = py::none(),
      py::arg("gang") = py::none(),
      py::arg("exact") = py::none(),
      py::arg("turn_autoscale_off") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_scale_cmd (where, y_min_in, y_max_in, axis, include_wall, gang, exact, turn_autoscale_off)

Routine to scale a plot.
If y_min = y_max, the scales will be chosen to show all the data.

Parameters
----------
where : str
    Region to scale. Eg: "top:x"

y_min_in : float
    Plot y-axis min value.

y_max_in : float
    Plot y-axis max value.

axis : str, optional
    'y', 'y2', or '' (both). Default = ''.

include_wall : bool, optional
    Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
    include the building wall position will be included in determining the the scale.

gang : str, optional
    'gang', 'nogang', ''. Default = ''.

exact : bool, optional
    Exact plot y_max, y_min to correspond to y_min_in, y_max_in? Default is False. Only relavent when y_min_in
    /= y_max_in.

turn_autoscale_off : bool, optional
    If True (default) then turn off plot.autoscale_y logical for all plots that are scaled.
)"""
  );
  py::class_<Tao::TaoScaleGraph, std::unique_ptr<Tao::TaoScaleGraph>>(
      m,
      "TaoScaleGraph",
      "tao_scale_graph return type"
  )
      .def_readonly("y_range", &Tao::TaoScaleGraph::y_range)
      .def_readonly("y2_range", &Tao::TaoScaleGraph::y2_range)
      .def("__len__", [](const Tao::TaoScaleGraph &) { return 2; })
      .def("__getitem__", [](const Tao::TaoScaleGraph &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.y_range);
        if (i == 1)
          return py::cast(s.y2_range);
        throw py::index_error();
      });
  m.def(
      "tao_scale_graph",
      &Tao::tao_scale_graph,
      py::arg("graph"),
      py::arg("y_min"),
      py::arg("y_max"),
      py::arg("axis") = py::none(),
      py::arg("include_wall") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_scale_graph (graph, y_min, y_max, axis, include_wall, y_range, y2_range)

Routine to scale the y-axis and/or y2-axis of a graph
If y_min = y_max then autoscaling will be done and the particular value of y_min and y_max is ignored.
Note: y_min/y_max is ignored if scaling y2-axis and graph%y2_mirrors_y = T.

Parameters
----------
graph : TaoGraphStruct
    Graph with axis/axes to be scaled.
    This parameter is an input/output and is modified in-place.
    As an output, graph: Graph with scaled axis/axes.

y_min : float
    Axis [min, max] must cover [y_min, y_max] if not autoscaling.

y_max : float
    Axis [min, max] must cover [y_min, y_max] if not autoscaling.

axis : str, optional
    Axis to scale. ''   -> scale y and y2 (default). 'y'  -> scale y-axis. 'y2' -> scale y2-axis

include_wall : bool, optional
    Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
    include the building wall position will be included in determining the the scale.

Returns
-------
y_range : 1D array of float (shape: 2), optional
    Only used by tao_scale_plot when ganging graphs.

y2_range : 1D array of float (shape: 2), optional
    Only used by tao_scale_plot when ganging graphs.
)"""
  );
  m.def(
      "tao_scale_ping_data",
      &Tao::tao_scale_ping_data,
      py::arg("u"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_scale_ping_data

Parameters
----------
u : TaoUniverseStruct
)"""
  );
  m.def(
      "tao_scale_plot",
      &Tao::tao_scale_plot,
      py::arg("plot"),
      py::arg("y_min_in"),
      py::arg("y_max_in"),
      py::arg("axis") = py::none(),
      py::arg("include_wall") = py::none(),
      py::arg("gang") = py::none(),
      py::arg("skip_lat_layout") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_scale_plot (plot, y_min_in, y_max_in, axis, include_wall, gang, skip_lat_layout)

Routine to scale the y-axis and/or y2-axis of the graphs of the plot.
If y_min_in = y_max_in then autoscaling will be done and the particular value
of y_min_in and y_max_in is ignored.

Parameters
----------
plot : TaoPlotStruct
    Plot with graphs to be scaled.
    This parameter is an input/output and is modified in-place.
    As an output, plot: Plot with scaled graphs.

y_min_in : float
    Axis [min, max] must cover [y_min_in, y_max_in] if not autoscaling.

y_max_in : float
    Axis [min, max] must cover [y_min_in, y_max_in] if not autoscaling.

axis : str, optional
    Axis to scale. ''   -> scale y and y2 (default). 'y'  -> scale y-axis. 'y2' -> scale y2-axis

include_wall : bool, optional
    Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
    include the building wall position will be included in determining the the scale.

gang : str, optional
    If autoscale then make all graph y-axes the same and/or make all y2-axes the same? ''        -> (default)
    Use setting of plot.autoscale_gang_y 'gang'    -> Gang graphs. 'nogang'  -> Do not gang graphs.

skip_lat_layout : bool, optional
    If True, skip scaling any lat_layout graphs. Default is false.
)"""
  );
  m.def(
      "tao_scratch_values_calc",
      &Tao::tao_scratch_values_calc,
      py::arg("ele_ref"),
      py::arg("ele_start"),
      py::arg("ele"),
      py::arg("datum"),
      py::arg("branch"),
      py::arg("orbit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_scratch_values_calc

Parameters
----------
ele_ref : EleStruct

ele_start : EleStruct

ele : EleStruct

datum : TaoDataStruct

branch : BranchStruct

orbit : 1D array of CoordStruct
)"""
  );
  m.def(
      "tao_set_beam_cmd",
      &Tao::tao_set_beam_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::arg("branch_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_beam_cmd (who, value_str, branch_str)

Routine to set various beam parameters.

Parameters
----------
who : str
    which parameter to set.

value_str : str
    Value to set to.

branch_str : str
    Branch to use. '' => branch 0.
)"""
  );
  m.def(
      "tao_set_beam_init_cmd",
      &Tao::tao_set_beam_init_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::arg("branch_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_beam_init_cmd (who, value_str, branch_str)

Routine to set beam_init variables

Parameters
----------
who : str
    which beam_init variable to set

value_str : str
    Value to set to.

branch_str : str
    Branch to use. '' => branch 0
)"""
  );
  m.def(
      "tao_set_bmad_com_cmd",
      &Tao::tao_set_bmad_com_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_bmad_com_cmd (who, value_str)

Routine to set bmad_com variables

Parameters
----------
who : str
    which bmad_com variable to set

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_branch_cmd",
      &Tao::tao_set_branch_cmd,
      py::arg("branch_str"),
      py::arg("component_str"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_branch_cmd (branch_str, component_str, value_str)

Routine to set lattice branch values.

Parameters
----------
branch_str : str
    Which branch to set.

component_str : str
    Which branch parameter to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_calculate_cmd",
      &Tao::tao_set_calculate_cmd,
      py::arg("switch_") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_calculate_cmd (switch)

Toggles off lattice calc and plotting.
)"""
  );
  m.def(
      "tao_set_curve_cmd",
      &Tao::tao_set_curve_cmd,
      py::arg("curve_name"),
      py::arg("component"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_curve_cmd (curve_name, component, value_str)

Routine to set var values.

Parameters
----------
curve_name : str
    Which curve to set.

component : str
    Which component to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_curve_invalid",
      &Tao::tao_set_curve_invalid,
      py::arg("curve"),
      py::arg("why_invalid"),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_curve_invalid (curve, why_invalid, print_err)

Routine to set curve%valid to False.

Parameters
----------
curve : TaoCurveStruct
    Curve to set.
    This parameter is an input/output and is modified in-place.
    As an output, curve: Curve properly set.

why_invalid : str
    Invalid information.

print_err : bool, optional
    If present and True then also print an error message.
)"""
  );
  m.def(
      "tao_set_data_cmd",
      &Tao::tao_set_data_cmd,
      py::arg("who_str"),
      py::arg("value_str"),
      py::arg("silent") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_data_cmd (who_str, value_str, silent)

Routine to set data values.

Parameters
----------
who_str : str
    Which data component(s) to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_data_useit_opt",
      &Tao::tao_set_data_useit_opt,
      py::arg("data") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_set_data_useit_opt

Parameters
----------
data : 1D array of TaoDataStruct, optional
    Data to work on. Default is all data in all universes.
)"""
  );
  m.def(
      "tao_set_default_cmd",
      &Tao::tao_set_default_cmd,
      py::arg("who_str"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_default_cmd (who_str, value_str)

Routine to set default values.

Parameters
----------
who_str : str
    Which default component(s) to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_drawing_cmd",
      &Tao::tao_set_drawing_cmd,
      py::arg("drawing"),
      py::arg("component"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_drawing_cmd (drawing, component, value_str)

Routine to set floor_plan and lat_layout parameters.

Parameters
----------
drawing : TaoDrawingStruct
    s.plot_page.floor_plan or s.plot_page.lat_layout.

component : str
    Which shape component to set.

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_dynamic_aperture_cmd",
      &Tao::tao_set_dynamic_aperture_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_dynamic_aperture_cmd (who, value_str)

Sets dynamic aperture parameters.

Parameters
----------
who : str
    which parameter to set.

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_elements_cmd",
      &Tao::tao_set_elements_cmd,
      py::arg("ele_list"),
      py::arg("attribute"),
      py::arg("value"),
      py::arg("update"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_elements_cmd (ele_list, attribute, value, update)

Sets element parameters.

Parameters
----------
ele_list : str
    which elements.

attribute : str
    Attribute to set.

value : str
    Value to set.
)"""
  );
  m.def(
      "tao_set_floor_plan_axis_label",
      &Tao::tao_set_floor_plan_axis_label,
      py::arg("graph"),
      py::arg("axis_in"),
      py::arg("axis_out"),
      py::arg("which"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_set_floor_plan_axis_label

Parameters
----------
graph : TaoGraphStruct

axis_in : QpAxisStruct

axis_out : QpAxisStruct

which : str
)"""
  );
  m.def(
      "tao_set_geodesic_lm_cmd",
      &Tao::tao_set_geodesic_lm_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_geodesic_lm_cmd (who, value_str)

Routine to set geodesic_lm variables

Parameters
----------
who : str
    which geodesic_lm variable to set

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_global_cmd",
      &Tao::tao_set_global_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_global_cmd (who, value_str)

Routine to set global variables

Parameters
----------
who : str
    which global variable to set

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_graph_cmd",
      &Tao::tao_set_graph_cmd,
      py::arg("graph_name"),
      py::arg("component"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_graph_cmd (graph_name, component, value_str)

Routine to set var values.

Parameters
----------
graph_name : str
    Which graph to set.

component : str
    Which component to set.

value_str : str
    What value to set it to.
)"""
  );
  py::class_<Tao::TaoSetIntegerValue, std::unique_ptr<Tao::TaoSetIntegerValue>>(
      m,
      "TaoSetIntegerValue",
      "tao_set_integer_value return type"
  )
      .def_readonly("var", &Tao::TaoSetIntegerValue::var)
      .def_readonly("error", &Tao::TaoSetIntegerValue::error)
      .def("__len__", [](const Tao::TaoSetIntegerValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetIntegerValue &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.var);
        if (i == 1)
          return py::cast(s.error);
        throw py::index_error();
      });
  m.def(
      "tao_set_integer_value",
      &Tao::tao_set_integer_value,
      py::arg("var_str"),
      py::arg("value_str"),
      py::arg("min_val") = py::none(),
      py::arg("max_val") = py::none(),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_integer_value (var, var_str, value_str, error, min_val, max_val, print_err)

Subroutine to read and set the value of an integer varialbe.

If the value is out of the range [min_val, max_val] then an error message will
be generated and the variable will not be set.

Parameters
----------
var_str : str
    Used for error messages.

value_str : str
    String with encoded value.

min_val : int, optional
    Minimum value.

max_val : int, optional
    Maximum value.

print_err : bool, optional
    If True, print error message. Default is true

Returns
-------
var : int
    Variable to set.

error : bool
    Set True on an error. False otherwise.
)"""
  );
  m.def(
      "tao_set_invalid",
      &Tao::tao_set_invalid,
      py::arg("datum"),
      py::arg("message"),
      py::arg("exterminate") = py::none(),
      py::arg("err_level") = py::none(),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_set_invalid

Parameters
----------
datum : TaoDataStruct
    Bad datum.

message : str
    Error message.

exterminate : bool, optional
    Default is False. If True, set datum.exists to False so that Tao will ignore this datum from now on.

err_level : int, optional
    s_error$ (default), s_warn$, etc.

print_err : bool, optional
    Default is True. If False, do not print an error message.

Returns
-------
why_invalid : str, optional
    Set to message if present.
)"""
  );
  m.def(
      "tao_set_key_cmd",
      &Tao::tao_set_key_cmd,
      py::arg("key_str"),
      py::arg("cmd_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_key_cmd (key_str, cmd_str)

Associates a command with a key press for single mode.

Parameters
----------
key_str : str
    keyboard key.

cmd_str : str
    Command associated with key.
)"""
  );
  m.def(
      "tao_set_lattice_cmd",
      &Tao::tao_set_lattice_cmd,
      py::arg("dest_lat"),
      py::arg("source_lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_lattice_cmd (dest_lat, source_lat)

Sets a lattice equal to another. This will also update the data structs

Parameters
----------
dest_lat : str
    Maybe: 'model', 'design', or 'base' with optional '@n' at beginning to indicate the universe

source_lat : str
    Maybe: 'model', 'design', or 'base'
)"""
  );
  py::class_<Tao::TaoSetLogicalValue, std::unique_ptr<Tao::TaoSetLogicalValue>>(
      m,
      "TaoSetLogicalValue",
      "tao_set_logical_value return type"
  )
      .def_readonly("var", &Tao::TaoSetLogicalValue::var)
      .def_readonly("error", &Tao::TaoSetLogicalValue::error)
      .def("__len__", [](const Tao::TaoSetLogicalValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetLogicalValue &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.var);
        if (i == 1)
          return py::cast(s.error);
        throw py::index_error();
      });
  m.def(
      "tao_set_logical_value",
      &Tao::tao_set_logical_value,
      py::arg("var_str"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_logical_value (var, var_str, value_str, error)

Subroutine to read and set the value of an logical varialbe.

If the value is out of the range [min_val, max_val] then an error message will
be generated and the variable will not be set.

Parameters
----------
var_str : str
    Used for error messages.

value_str : str
    String with encoded value.

Returns
-------
var : bool
    Variable to set.

error : bool
    Set True on an error. False otherwise.
)"""
  );
  m.def(
      "tao_set_openmp_n_threads",
      &Tao::tao_set_openmp_n_threads,
      py::arg("n_threads"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_openmp_n_threads (n_threads)

Routine to set OpenMP thread count.  Errors if OpenMP is not available.

Parameters
----------
n_threads : int
    Number of threads.
)"""
  );
  m.def(
      "tao_set_opt_vars",
      &Tao::tao_set_opt_vars,
      py::arg("var_vec"),
      py::arg("print_limit_warning") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_set_opt_vars

Parameters
----------
var_vec : 1D array of float
    Vector of variables.

print_limit_warning : bool, optional
    Print a warning if the value is past the variable's limits. Default is True.
)"""
  );
  m.def(
      "tao_set_opti_de_param_cmd",
      &Tao::tao_set_opti_de_param_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_opti_de_param_cmd (who, value_str)

Routine to set opti_de_param variables

Parameters
----------
who : str
    which opti_de_param variable to set

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_particle_start_cmd",
      &Tao::tao_set_particle_start_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_particle_start_cmd (who, value_str)

Routine to set particle_start variables.

Parameters
----------
who : str
    which particle_start variable to set

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_plot_cmd",
      &Tao::tao_set_plot_cmd,
      py::arg("plot_name"),
      py::arg("component"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_plot_cmd (plot_name, component, value_str)

Routine to set plot parameters.

Parameters
----------
plot_name : str
    Which plot to set.

component : str
    Which component to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_plot_page_cmd",
      &Tao::tao_set_plot_page_cmd,
      py::arg("component"),
      py::arg("value_str"),
      py::arg("value_str2") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_plot_page_cmd (component, value_str, value_str2)

 Set various aspects of the plotting window

Parameters
----------
component : str
    Which component to set.

value_str : str
    What value to set to.

value_str2 : str, optional
    2nd value if component is an array.
)"""
  );
  m.def(
      "tao_set_ptc_com_cmd",
      &Tao::tao_set_ptc_com_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_ptc_com_cmd (who, value_str)

Routine to set ptc_com variables

Parameters
----------
who : str
    which ptc_com variable to set

value_str : str
    Value to set to.
)"""
  );
  py::class_<Tao::TaoSetQpAxisStruct, std::unique_ptr<Tao::TaoSetQpAxisStruct>>(
      m,
      "TaoSetQpAxisStruct",
      "tao_set_qp_axis_struct return type"
  )
      .def_readonly("error", &Tao::TaoSetQpAxisStruct::error)
      .def_readonly("ix_uni", &Tao::TaoSetQpAxisStruct::ix_uni)
      .def("__len__", [](const Tao::TaoSetQpAxisStruct &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetQpAxisStruct &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.error);
        if (i == 1)
          return py::cast(s.ix_uni);
        throw py::index_error();
      });
  m.def(
      "tao_set_qp_axis_struct",
      &Tao::tao_set_qp_axis_struct,
      py::arg("qp_axis_name"),
      py::arg("component"),
      py::arg("qp_axis"),
      py::arg("value"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_qp_axis_struct (qp_axis_name, component, qp_axis, value, error, ix_uni)

Routine to set qp_axis_names of a qp_axis_struct.

Parameters
----------
qp_axis_name : str
    qp_axis name. Used for error messages.

component : str
    qp_axis component name.

qp_axis : QpAxisStruct
    qp_axis_struct with component to modify
    This parameter is an input/output and is modified in-place.
    As an output, qp_axis: qp_axis_struct with changed component value.

value : str
    Component value.

Returns
-------
error : bool
    Set true if there is an error. False otherwise.

ix_uni : int, optional
    Tao universe number in case the value depends upon a parameter of a particular universe.
)"""
  );
  py::class_<Tao::TaoSetQpPointStruct, std::unique_ptr<Tao::TaoSetQpPointStruct>>(
      m,
      "TaoSetQpPointStruct",
      "tao_set_qp_point_struct return type"
  )
      .def_readonly("error", &Tao::TaoSetQpPointStruct::error)
      .def_readonly("ix_uni", &Tao::TaoSetQpPointStruct::ix_uni)
      .def("__len__", [](const Tao::TaoSetQpPointStruct &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetQpPointStruct &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.error);
        if (i == 1)
          return py::cast(s.ix_uni);
        throw py::index_error();
      });
  m.def(
      "tao_set_qp_point_struct",
      &Tao::tao_set_qp_point_struct,
      py::arg("qp_point_name"),
      py::arg("component"),
      py::arg("qp_point"),
      py::arg("value"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_qp_point_struct (qp_point_name, component, qp_point, value, error, ix_uni)

Routine to set qp_point_names of a qp_point_struct.

Parameters
----------
qp_point_name : str
    qp_point name. Used for error messages.

component : str
    qp_point component name.

qp_point : QpPointStruct
    qp_point_struct with component to modify
    This parameter is an input/output and is modified in-place.
    As an output, qp_point: qp_point_struct with changed component value.

value : str
    Component value.

Returns
-------
error : bool
    Set true if there is an error. False otherwise.

ix_uni : int, optional
    Tao universe number in case the value depends upon a parameter of a particular universe.
)"""
  );
  py::class_<Tao::TaoSetQpRectStruct, std::unique_ptr<Tao::TaoSetQpRectStruct>>(
      m,
      "TaoSetQpRectStruct",
      "tao_set_qp_rect_struct return type"
  )
      .def_readonly("error", &Tao::TaoSetQpRectStruct::error)
      .def_readonly("ix_uni", &Tao::TaoSetQpRectStruct::ix_uni)
      .def("__len__", [](const Tao::TaoSetQpRectStruct &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetQpRectStruct &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.error);
        if (i == 1)
          return py::cast(s.ix_uni);
        throw py::index_error();
      });
  m.def(
      "tao_set_qp_rect_struct",
      &Tao::tao_set_qp_rect_struct,
      py::arg("qp_rect_name"),
      py::arg("component"),
      py::arg("qp_rect"),
      py::arg("value"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_qp_rect_struct (qp_rect_name, component, qp_rect, value, error, ix_uni)

Routine to set qp_rect_names of a qp_rect_struct.

Parameters
----------
qp_rect_name : str
    qp_rect name. Used for error messages.

component : str
    qp_rect component name.

qp_rect : QpRectStruct
    qp_rect_struct with component to modify
    This parameter is an input/output and is modified in-place.
    As an output, qp_rect: qp_rect_struct with changed component value.

value : str
    Component value.

Returns
-------
error : bool
    Set true if there is an error. False otherwise.

ix_uni : int, optional
    Tao universe number in case the value depends upon a parameter of a particular universe.
)"""
  );
  m.def(
      "tao_set_ran_state_cmd",
      &Tao::tao_set_ran_state_cmd,
      py::arg("state_string"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_ran_state_cmd (state_string)

Sets the random number generator state.

Parameters
----------
state_string : str
    Encoded random number state.
)"""
  );
  py::class_<Tao::TaoSetRealValue, std::unique_ptr<Tao::TaoSetRealValue>>(
      m,
      "TaoSetRealValue",
      "tao_set_real_value return type"
  )
      .def_readonly("var", &Tao::TaoSetRealValue::var)
      .def_readonly("error", &Tao::TaoSetRealValue::error)
      .def("__len__", [](const Tao::TaoSetRealValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetRealValue &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.var);
        if (i == 1)
          return py::cast(s.error);
        throw py::index_error();
      });
  m.def(
      "tao_set_real_value",
      &Tao::tao_set_real_value,
      py::arg("var_str"),
      py::arg("value_str"),
      py::arg("min_val") = py::none(),
      py::arg("max_val") = py::none(),
      py::arg("dflt_uni") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_real_value (var, var_str, value_str, error, min_val, max_val, dflt_uni)

Subroutine to read and set the value of a real variable.

If the value is out of the range [min_val, max_val] then an error message will
be generated and the variable will not be set.

Parameters
----------
var_str : str
    Used for error messages.

value_str : str
    String with encoded value.

min_val : float, optional
    Minimum value.

max_val : float, optional
    Maximum value.

dflt_uni : int, optional
    Default universe used to evaluate parameters.

Returns
-------
var : float
    Variable to set.

error : bool
    Set True on an error. False otherwise.
)"""
  );
  m.def(
      "tao_set_region_cmd",
      &Tao::tao_set_region_cmd,
      py::arg("region_name"),
      py::arg("component"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_region_cmd (region_name, component, value_str)

Routine to set region parameters.

Parameters
----------
region_name : str
    Which region to set.

component : str
    Which component to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_space_charge_com_cmd",
      &Tao::tao_set_space_charge_com_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_space_charge_com_cmd (who, value_str)

Routine to set space_charge_com variables

Parameters
----------
who : str
    which space_charge_com variable to set

value_str : str
    Value to set to.
)"""
  );
  m.def(
      "tao_set_symbolic_number_cmd",
      &Tao::tao_set_symbolic_number_cmd,
      py::arg("sym_str"),
      py::arg("num_str") = py::none(),
      py::arg("val") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_symbolic_number_cmd (sym_str, num_str, val)

Associates a given symbol with a given number.
Note: Either num_str or val argument must be present.

Parameters
----------
sym_str : str
    Symbol.

num_str : str, optional
    Symbol value expression.

val : float, optional
    Value of symbol
)"""
  );
  m.def(
      "tao_set_tune_cmd",
      &Tao::tao_set_tune_cmd,
      py::arg("branch_str"),
      py::arg("mask_str"),
      py::arg("print_list"),
      py::arg("qa_str"),
      py::arg("qb_str"),
      py::arg("delta_input"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_tune_cmd (branch_str, mask_str, print_list, qa_str, qb_str, delta_input)

Routine to set the transverse tunes.

Parameters
----------
branch_str : str
    List of branches to apply tune set to.

mask_str : str
    List of quadrupoles to veto.

print_list : bool
    If True, print a list of elements varied and coefficients.

qa_str : str
    Expression for Qa tune.

qb_str : str
    Expression for Qb tune.

delta_input : bool
    If true then qa_str and qb_str are deltas from present tune.
)"""
  );
  m.def(
      "tao_set_universe_cmd",
      &Tao::tao_set_universe_cmd,
      py::arg("uni"),
      py::arg("who"),
      py::arg("what"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_universe_cmd (uni, who, what)

Sets a universe on or off, or sets the recalculate or twiss_calc logicals, etc.

Parameters
----------
uni : str
    which universe; 0 => current viewed universe

who : str
    "on", "off", "recalculate", "dynamic_aperture_calc", "one_turn_map_calc", or "twiss_calc"

what : str
    "on" or "off" for who = "dynamic_aperture_calc", "one_turn_map_calc" or "twiss_calc".
)"""
  );
  m.def(
      "tao_set_var_cmd",
      &Tao::tao_set_var_cmd,
      py::arg("var_str"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_var_cmd (var_str, value_str)

Routine to set var values.

Parameters
----------
var_str : str
    Which var name to set.

value_str : str
    What value to set it to.
)"""
  );
  m.def(
      "tao_set_var_model_value",
      &Tao::tao_set_var_model_value,
      py::arg("var"),
      py::arg("value"),
      py::arg("print_limit_warning") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_set_var_model_value

Parameters
----------
var : TaoVarStruct
    Variable to set

value : float
    Value to set to

print_limit_warning : bool, optional
    Print a warning if the value is past the variable's limits. Default is True.
)"""
  );
  m.def(
      "tao_set_var_useit_opt",
      &Tao::tao_set_var_useit_opt,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_set_var_useit_opt
)"""
  );
  m.def(
      "tao_set_wave_cmd",
      &Tao::tao_set_wave_cmd,
      py::arg("who"),
      py::arg("value_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_wave_cmd (who, value_str, err)

Routine to set wave variables

Parameters
----------
who : str
    which wave variable to set

value_str : str
    Value to set to.

Returns
-------
err : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "tao_set_z_tune_cmd",
      &Tao::tao_set_z_tune_cmd,
      py::arg("branch_str"),
      py::arg("q_str"),
      py::arg("delta_input"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_set_z_tune_cmd (branch_str, q_str, delta_input)

Routine to set the z-tune.

Parameters
----------
branch_str : str
    List of branches to apply tune set to.

q_str : str
    Expression for Qc tune.

delta_input : bool
    If true then qa_str and qb_str are deltas from present tune.
)"""
  );
  m.def(
      "tao_setup_key_table",
      &Tao::tao_setup_key_table,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_setup_key_table
)"""
  );
  m.def(
      "tao_shape_init",
      &Tao::tao_shape_init,
      py::arg("shape"),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_shape_init

Parameters
----------
shape : TaoEleShapeStruct
    Shape

print_err : bool, optional
    If True then print an error message if there is a problem. Default is True.

Returns
-------
err : bool
    Set true if there is a problem translating the element class.
)"""
  );
  m.def(
      "tao_show_cmd",
      &Tao::tao_show_cmd,
      py::arg("what"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_show_cmd

Parameters
----------
what : str
    What to show.
)"""
  );
  m.def(
      "tao_show_constraints",
      &Tao::tao_show_constraints,
      py::arg("iunit"),
      py::arg("form"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_show_constraints (iunit, form)

Routine to show a list of datums and variables and how they contribute to the merit function.

Parameters
----------
iunit : int
    File unit to write to. 0 => print to the terminal.

form : str
    What to output: 'ALL'   -> All datums and variables. 'TOP10' -> Top datums and variables that contribute
    to the merit function.
)"""
  );
  py::class_<Tao::TaoShowThis, std::unique_ptr<Tao::TaoShowThis>>(
      m,
      "TaoShowThis",
      "tao_show_this return type"
  )
      .def_readonly("result_id", &Tao::TaoShowThis::result_id)
      .def_readonly("lines", &Tao::TaoShowThis::lines)
      .def_readonly("nl", &Tao::TaoShowThis::nl)
      .def("__len__", [](const Tao::TaoShowThis &) { return 3; })
      .def("__getitem__", [](const Tao::TaoShowThis &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.result_id);
        if (i == 1)
          return py::cast(s.lines);
        if (i == 2)
          return py::cast(s.nl);
        throw py::index_error();
      });
  m.def(
      "tao_show_this",
      &Tao::tao_show_this,
      py::arg("what"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_show_this

Parameters
----------
what : str
    What to show.

Returns
-------
result_id : str
    ID string idendifying what was shown. Used by the Tao GUI.

lines : 1D array of str
    Output.

nl : int
    Number of lines(:).
)"""
  );
  m.def(
      "tao_single_mode",
      &Tao::tao_single_mode,
      py::arg("char_"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_single_mode

Parameters
----------
)"""
  );
  m.def(
      "tao_single_track",
      &Tao::tao_single_track,
      py::arg("tao_lat"),
      py::arg("ix_branch"),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_single_track (tao_lat, calc_ok, ix_branch, print_err)

Routine to track a single particle and calculate lattice functions through a lattice.

Parameters
----------
tao_lat : TaoLatticeStruct
    Structure containing the lattice.

ix_branch : int
    Branch index to track through.

print_err : bool, optional
    Default False. Print error messages if, eg, lattice is unstable?

Returns
-------
calc_ok : bool
    Set True if there were no problems, False otherwise.
)"""
  );
  m.def(
      "tao_spin_matrices_calc_needed",
      &Tao::tao_spin_matrices_calc_needed,
      py::arg("data_type"),
      py::arg("data_source"),
      py::arg("do_calc"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_spin_matrices_calc_needed

Parameters
----------
data_type : str

data_source : str

do_calc : bool
)"""
  );
  m.def(
      "tao_spin_tracking_turn_on",
      &Tao::tao_spin_tracking_turn_on,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_spin_tracking_turn_on
)"""
  );
  py::class_<Tao::TaoSplitComponent, std::unique_ptr<Tao::TaoSplitComponent>>(
      m,
      "TaoSplitComponent",
      "tao_split_component return type"
  )
      .def_readonly("comp", &Tao::TaoSplitComponent::comp)
      .def_readonly("err", &Tao::TaoSplitComponent::err)
      .def("__len__", [](const Tao::TaoSplitComponent &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSplitComponent &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.comp);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "tao_split_component",
      &Tao::tao_split_component,
      py::arg("comp_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_split_component

Parameters
----------
comp_str : str
    Components. EG: 'meas - design'

Returns
-------
comp : 1D array of TaoDataVarComponentStruct
    Array of individual components.

err : bool
    Set True if there is an error, False otherwise.
)"""
  );
  m.def(
      "tao_srdt_calc_needed",
      &Tao::tao_srdt_calc_needed,
      py::arg("data_type"),
      py::arg("data_source"),
      py::arg("do_srdt"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_srdt_calc_needed

Parameters
----------
data_type : str

data_source : str

do_srdt : int
)"""
  );
  py::class_<Tao::TaoSubinUniNumber, std::unique_ptr<Tao::TaoSubinUniNumber>>(
      m,
      "TaoSubinUniNumber",
      "tao_subin_uni_number return type"
  )
      .def_readonly("name_out", &Tao::TaoSubinUniNumber::name_out)
      .def_readonly("ok", &Tao::TaoSubinUniNumber::ok)
      .def("__len__", [](const Tao::TaoSubinUniNumber &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSubinUniNumber &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.name_out);
        if (i == 1)
          return py::cast(s.ok);
        throw py::index_error();
      });
  m.def(
      "tao_subin_uni_number",
      &Tao::tao_subin_uni_number,
      py::arg("name_in"),
      py::arg("ix_uni"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_subin_uni_number

Parameters
----------
name_in : str
    Input name with "#" character

ix_uni : int
    Universe index.

Returns
-------
name_out : str
    Output name.

ok : bool
    False if multiple universes and no "#" in name_in. True otherwise.
)"""
  );
  m.def(
      "tao_svd_optimizer",
      &Tao::tao_svd_optimizer,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_svd_optimizer (abort)

Routine to minimize the merit function using svd.

Returns
-------
abort : bool
    Set True if svd step increases the merit function.
)"""
  );
  m.def(
      "tao_symbol_import_from_lat",
      &Tao::tao_symbol_import_from_lat,
      py::arg("lat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_symbol_import_from_lat

Parameters
----------
lat : LatStruct
)"""
  );
  m.def(
      "tao_taper_cmd",
      &Tao::tao_taper_cmd,
      py::arg("except_"),
      py::arg("uni_names"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_taper_cmd

Parameters
----------
except : str
    List of elements not to vary.

uni_names : str
    Universes to taper.
)"""
  );
  m.def(
      "tao_to_change_number",
      &Tao::tao_to_change_number,
      py::arg("num_str"),
      py::arg("n_size"),
      py::arg("change_number"),
      py::arg("abs_or_rel"),
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_to_change_number

Parameters
----------
num_str : str

n_size : int

change_number : 1D array of float

abs_or_rel : str

err : bool
)"""
  );
  m.def(
      "tao_to_int",
      &Tao::tao_to_int,
      py::arg("str"),
      py::arg("i_int"),
      py::arg("err"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_to_int (str, i_int, err)

Converts a string to an integer

If the string str is blank then i_int = 0
)"""
  );
  py::class_<Tao::TaoToPhaseAndCouplingReading, std::unique_ptr<Tao::TaoToPhaseAndCouplingReading>>(
      m,
      "TaoToPhaseAndCouplingReading",
      "tao_to_phase_and_coupling_reading return type"
  )
      .def_readonly("bpm_data", &Tao::TaoToPhaseAndCouplingReading::bpm_data)
      .def_readonly("valid_value", &Tao::TaoToPhaseAndCouplingReading::valid_value)
      .def("__len__", [](const Tao::TaoToPhaseAndCouplingReading &) { return 2; })
      .def("__getitem__", [](const Tao::TaoToPhaseAndCouplingReading &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.bpm_data);
        if (i == 1)
          return py::cast(s.valid_value);
        throw py::index_error();
      });
  m.def(
      "tao_to_phase_and_coupling_reading",
      &Tao::tao_to_phase_and_coupling_reading,
      py::arg("ele"),
      py::arg("why_invalid"),
      py::arg("datum"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)

Buffer routine for to_phase_and_coupling_reading.

Parameters
----------
ele : EleStruct
    The monitor.

Returns
-------
bpm_data : BpmPhaseCouplingStruct
    Monitor values

valid_value : bool
    Valid data value?
)"""
  );
  py::class_<Tao::TaoToReal, std::unique_ptr<Tao::TaoToReal>>(
      m,
      "TaoToReal",
      "tao_to_real return type"
  )
      .def_readonly("value", &Tao::TaoToReal::value)
      .def_readonly("err_flag", &Tao::TaoToReal::err_flag)
      .def("__len__", [](const Tao::TaoToReal &) { return 2; })
      .def("__getitem__", [](const Tao::TaoToReal &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.value);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "tao_to_real",
      &Tao::tao_to_real,
      py::arg("expression"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_to_real

Parameters
----------
expression : str
    arithmetic expression

Returns
-------
value : float
    Value of arithmetic expression.

err_flag : bool
    TRUE on error.
)"""
  );
  m.def(
      "tao_too_many_particles_lost",
      &Tao::tao_too_many_particles_lost,
      py::arg("beam"),
      py::arg("no_beam"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_too_many_particles_lost

Parameters
----------
beam : BeamStruct

no_beam : bool
)"""
  );
  m.def(
      "tao_top10_derivative_print",
      &Tao::tao_top10_derivative_print,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_top10_derivative_print ()

Routine to print out the top10 contributors to the merit function.
)"""
  );
  m.def(
      "tao_top10_merit_categories_print",
      &Tao::tao_top10_merit_categories_print,
      py::arg("iunit"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_top10_merit_categories_print (iunit)

Routine to print the top data and variable categories that contribute to
the merit function.

Parameters
----------
iunit : int
    File unit to write to. 0 => print to the terminal.
)"""
  );
  m.def(
      "tao_top_level",
      &Tao::tao_top_level,
      py::arg("command") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_top_level

Parameters
----------
command : str, optional
    Tao command string. If present, getting user input from the terminal is bypassed. This is used when
    interfacing to Python.

Returns
-------
errcode : int, optional
    Return error code: 0 => OK, Not 0 => Err.
)"""
  );
  py::class_<Tao::TaoTrackingEleIndex, std::unique_ptr<Tao::TaoTrackingEleIndex>>(
      m,
      "TaoTrackingEleIndex",
      "tao_tracking_ele_index return type"
  )
      .def_readonly("ix_branch", &Tao::TaoTrackingEleIndex::ix_branch)
      .def_readonly("ix_ele", &Tao::TaoTrackingEleIndex::ix_ele)
      .def("__len__", [](const Tao::TaoTrackingEleIndex &) { return 2; })
      .def("__getitem__", [](const Tao::TaoTrackingEleIndex &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ix_branch);
        if (i == 1)
          return py::cast(s.ix_ele);
        throw py::index_error();
      });
  m.def(
      "tao_tracking_ele_index",
      &Tao::tao_tracking_ele_index,
      py::arg("ele"),
      py::arg("datum"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_tracking_ele_index(ele, datum, ix_branch) result (ix_ele)

Routine to return the index in the tracking part of a lattice that corresponds to ele.

Parameters
----------
ele : EleStruct
    Lattice element.

datum : TaoDataStruct
    Datum

Returns
-------
ix_ele : int
    Element index associated with ele.

ix_branch : int, optional
    Lattice branch associated with element
)"""
  );
  m.def(
      "tao_turn_on_special_calcs_if_needed_for_plotting",
      &Tao::tao_turn_on_special_calcs_if_needed_for_plotting,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_turn_on_special_calcs_if_needed_for_plotting
)"""
  );
  m.def(
      "tao_type_expression_tree",
      &Tao::tao_type_expression_tree,
      py::arg("tree"),
      py::arg("indent") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_type_expression_tree (tree, indent)

Routine to print an expression tree in tree form.
Good for debugging.

Parameters
----------
tree : TaoEvalNodeStruct
    Tree to print.

indent : int, optional
    Initial indent. Default is zero.
)"""
  );
  m.def(
      "tao_uni_atsign_index",
      &Tao::tao_uni_atsign_index,
      py::arg("string"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function tao_uni_atsign_index(string) result (ix_amp)

Routine to return the index of an atsign ("@") character in a string if the atsign is
being used as a separator between a universe spec and the rest of the string.

For example:
  string = "[1:3]@orbit.x[5] => ix_amp = 6
  string = "orbit.x[5@0.2]   => ix_amp = 0 (no universe "@" present)

Parameters
----------
string : str
    String to parse

Returns
-------
ix_amp : int
    Index of universe "@". Set to zero if no universe "@" found.
)"""
  );
  m.def(
      "tao_universe_index",
      &Tao::tao_universe_index,
      py::arg("i_uni"),
      py::arg("neg2_to_default") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_universe_index

Parameters
----------
i_uni : int
    Nominal universe number.

neg2_to_default : bool, optional
    i_uni = -2 (all universes) maps to the default uni? Default if False.

Returns
-------
i_this_uni : int
    Universe number.
)"""
  );
  m.def(
      "tao_use_data",
      &Tao::tao_use_data,
      py::arg("action"),
      py::arg("data_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_use_data

Parameters
----------
action : str
    veto, use or restore

data_name : str
    the selected data
)"""
  );
  m.def(
      "tao_use_var",
      &Tao::tao_use_var,
      py::arg("action"),
      py::arg("var_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_use_var

Parameters
----------
action : str
    'use', 'veto', or 'restore'

var_name : str
    the selected variable name or all
)"""
  );
  m.def(
      "tao_user_is_terminating_optimization",
      &Tao::tao_user_is_terminating_optimization,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_user_is_terminating_optimization

Returns
-------
is_terminating : bool
    Set True of '.' is detected. False otherwise.
)"""
  );
  m.def(
      "tao_var1_name",
      &Tao::tao_var1_name,
      py::arg("var"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var1_name

Parameters
----------
var : TaoVarStruct
    Variable

Returns
-------
var1_name : str
    Appropriate name.
)"""
  );
  m.def(
      "tao_var_attrib_name",
      &Tao::tao_var_attrib_name,
      py::arg("var"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var_attrib_name

Parameters
----------
var : TaoVarStruct
    Variable

Returns
-------
var_attrib_name : str
    Attribute list.
)"""
  );
  m.def(
      "tao_var_check",
      &Tao::tao_var_check,
      py::arg("eles"),
      py::arg("attribute"),
      py::arg("silent"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var_check

Parameters
----------
eles : 1D array of ElePointerStruct
    Array of elements which have a changed attribute.

attribute : str
    Name of attribute changed.

silent : bool
    If True and the problem can be fixed, do not issue an error message.
)"""
  );
  m.def(
      "tao_var_repoint",
      &Tao::tao_var_repoint,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var_repoint
)"""
  );
  m.def(
      "tao_var_show_use",
      &Tao::tao_var_show_use,
      py::arg("v1_var"),
      py::arg("lines") = py::none(),
      py::arg("nl") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var_show_use

Parameters
----------
v1_var : TaoV1VarStruct
    tao_v1_var_struct

lines : 1D array of str, optional

nl : int, optional
)"""
  );
  m.def(
      "tao_var_target_calc",
      &Tao::tao_var_target_calc,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var_target_calc
)"""
  );
  m.def(
      "tao_var_useit_plot_calc",
      &Tao::tao_var_useit_plot_calc,
      py::arg("graph"),
      py::arg("var"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_var_useit_plot_calc

Parameters
----------
graph : TaoGraphStruct

var : 1D array of TaoVarStruct
)"""
  );
  m.def(
      "tao_var_write",
      &Tao::tao_var_write,
      py::arg("out_file"),
      py::arg("show_good_opt_only") = py::none(),
      py::arg("tao_format") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_var_write (out_file, show_good_opt_only, tao_format)

Routine to write the optimized variables. One file will be created for each universe.
The created file will have three sections:
  1) The variable values
  2) The list of constraints.
  3) A list of the top 10 constraints.
If out_file = '' the information will be dumped to the terminal.
In this case, only the variable values will be printed.

When tao_format = True, the output is in the form "set variable <name> = <value>"
so the file can be used as a Tao command file. If tao_format = False, the format
is suitable for inclusion in a Bmad lattice file.

Parameters
----------
out_file : str
    Name of output file. If blank. Ouput to the terminal.

show_good_opt_only : bool, optional
    Write only the variables used in the optimization? Default is False.

tao_format : bool, optional
    Output format. Default False. See above.
)"""
  );
  m.def(
      "tao_veto_vars_with_zero_dmodel",
      &Tao::tao_veto_vars_with_zero_dmodel,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_veto_vars_with_zero_dmodel ()

Routine to veto all variables with zero effect on data used in the merit function.
)"""
  );
  m.def(
      "tao_wave_analysis",
      &Tao::tao_wave_analysis,
      py::arg("plot"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_wave_analysis (plot)

Routine to do a wave anaysis.

Parameters
----------
plot : TaoPlotStruct
    Plot region setup by tao_wave_cmd.
    This parameter is an input/output and is modified in-place.
    As an output, plot: Plot with wave analysis curves.
)"""
  );
  m.def(
      "tao_wave_cmd",
      &Tao::tao_wave_cmd,
      py::arg("curve_name"),
      py::arg("plot_place"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_wave_cmd (curve_name, plot_place, err_flag)

Routine to do the initial setup for wave plotting.
The wave analysis is done by the routine tao_wave_analysis.

Parameters
----------
curve_name : str
    Character(*) curve for wave analysis.

plot_place : str
    Character(*) place on plot page to put the wave plot.
)"""
  );
  m.def(
      "tao_wave_fit",
      &Tao::tao_wave_fit,
      py::arg("curve"),
      py::arg("ix1"),
      py::arg("n_dat"),
      py::arg("coef"),
      py::arg("rms"),
      py::arg("f1"),
      py::arg("f2") = py::none(),
      py::arg("f3") = py::none(),
      py::arg("f4") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_wave_fit (curve, ix1, n_dat, coef, rms, f1, f2, f3, f4)

Routine for fitting the curve data to up to four functions using a least squares
SVD fit.

Parameters
----------
curve : TaoCurveStruct
    Curve containing the data.

ix1 : int
    Index of first point in the data array.

n_dat : int
    Number of data points.

coef : 1D array of float
    Fit coefficients.

rms : 1D array of float
    Variances with rms(n_func+1) = sqrt(chi^2/n_dat).

f1 : 1D array of float
    First fit function.

f2 : 1D array of float, optional
    Second fit function.

f3 : 1D array of float, optional
    third fit function.

f4 : 1D array of float, optional
    fourth fit function.
)"""
  );
  m.def(
      "tao_write_cmd",
      &Tao::tao_write_cmd,
      py::arg("what"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_write_cmd

Parameters
----------
what : str
    What to output. See the code for more details.
)"""
  );
  m.def(
      "tao_write_lines",
      &Tao::tao_write_lines,
      py::arg("iunit"),
      py::arg("line"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_write_lines (iunit, line)

Subroutine to write out a series of lines to a file or to the terminal.
It is assumed that any file has already been opened.

Parameters
----------
iunit : int
    File unit to write to. 0 => print to the terminal.

line : 1D array of str
    A series of lines.
)"""
  );
  m.def(
      "tao_x_axis_cmd",
      &Tao::tao_x_axis_cmd,
      py::arg("where"),
      py::arg("what"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_x_axis_cmd

Parameters
----------
where : str
    Region to axis. Eg: "top"

what : str
    "s" or "index"
)"""
  );
  m.def(
      "tao_x_scale_cmd",
      &Tao::tao_x_scale_cmd,
      py::arg("where"),
      py::arg("x_min_in"),
      py::arg("x_max_in"),
      py::arg("include_wall") = py::none(),
      py::arg("gang") = py::none(),
      py::arg("exact") = py::none(),
      py::arg("turn_autoscale_off") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_x_scale_cmd (where, x_min_in, x_max_in, err, include_wall, gang, exact, turn_autoscale_off)

Routine to scale a plot. If x_min = x_max
Then the scales will be chosen to show all the data.

Parameters
----------
where : str
    Region to scale. Eg: "top"

x_min_in : float
    Plot x-axis min value.

x_max_in : float
    Plot x-axis max value.

include_wall : bool, optional
    Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
    include the building wall position will be included in determining the the scale.

gang : str, optional
    'gang', 'nogang', ''. Default = ''.

exact : bool, optional
    Exact plot y_max, y_min to correspond to y_min_in, y_max_in? Default is False. Only relavent when y_min_in
    /= y_max_in.

turn_autoscale_off : bool, optional
    If True (default) then turn off plot.autoscale_x logical for all plots that are scaled.

Returns
-------
err : bool
    Set to True if the plot cannot be found. False otherwise.
)"""
  );
  m.def(
      "tao_x_scale_graph",
      &Tao::tao_x_scale_graph,
      py::arg("graph"),
      py::arg("x_min"),
      py::arg("x_max"),
      py::arg("include_wall") = py::none(),
      py::arg("have_scaled") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine tao_x_scale_graph

Parameters
----------
graph : TaoGraphStruct

x_min : float

x_max : float

include_wall : bool, optional

have_scaled : bool, optional
)"""
  );
  m.def(
      "tao_x_scale_plot",
      &Tao::tao_x_scale_plot,
      py::arg("plot"),
      py::arg("x_min_in"),
      py::arg("x_max_in"),
      py::arg("include_wall") = py::none(),
      py::arg("gang") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine tao_x_scale_plot (plot, x_min_in, x_max_in, include_wall, gang, have_scaled)

Routine to scale a plot. If x_min = x_max
Then the scales will be chosen to show all the data.

Parameters
----------
plot : TaoPlotStruct
    Plot to scale. Eg: "top"

x_min_in : float
    Plot x-axis min value.

x_max_in : float
    Plot x-axis max value.

include_wall : bool, optional
    Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
    include the building wall position will be included in determining the the scale.

gang : str, optional
    'gang', 'nogang', ''. Default = ''.

Returns
-------
have_scaled : bool, optional
    Has a graph been scaled?
)"""
  );
}
