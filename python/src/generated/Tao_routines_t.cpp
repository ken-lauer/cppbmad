#include "pybmad/generated/Tao_routines_t.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
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

void init_Tao_routines_t(nb::module_ &m) {
  m.def(
      "tao_abort_command_file",
      &Tao::tao_abort_command_file,
      nb::arg("force_abort") = nb::none(),
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
      nb::arg("h_str"),
      R"""(Routine to add on to the "h(:)" array holding the list of normal form
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
      nb::arg("alias"),
      nb::arg("string"),
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
      nb::arg("u"),
      nb::arg("n_data"),
      nb::arg("exact") = nb::none(),
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
      nb::arg("n_v1"),
      nb::arg("save_old"),
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
      nb::arg("n_var"),
      nb::arg("default_good_user"),
      R"""(Routine to increase the s%var(:) array size.

Parameters
----------
n_var : int
    Size of s.var(:) wanted.
)"""
  );
  m.def(
      "tao_beam_emit_calc",
      &Tao::tao_beam_emit_calc,
      nb::arg("plane"),
      nb::arg("emit_type"),
      nb::arg("ele"),
      nb::arg("bunch_params"),
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
      nb::arg("u"),
      nb::arg("tao_lat"),
      nb::arg("ix_branch"),
      nb::arg("beam"),
      R"""(Routine to track a a beam of particles.

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
      nb::arg("ele_id"),
      nb::arg("lat"),
      nb::arg("branch_str"),
      nb::arg("where"),
      nb::arg("u"),
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
      nb::arg("ix_branch"),
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
      nb::arg("tao_lat"),
      nb::arg("curve"),
      nb::arg("comp_sign"),
      nb::arg("good"),
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
      [](std::string file_name, CharacterAlloc1D *cmd_arg) {
        auto fn =
            static_cast<void (*)(std::string, optional_ref<CharacterAlloc1D>)>(&Tao::tao_call_cmd);
        return fn(file_name, ptr_to_opt_ref(cmd_arg));
      },
      nb::arg("file_name"),
      nb::arg("cmd_arg") = nb::none(),
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
      nb::arg("plot"),
      R"""(No docstring available.
)"""
  );
  m.def(
      "tao_change_ele",
      &Tao::tao_change_ele,
      nb::arg("ele_name"),
      nb::arg("attrib_name"),
      nb::arg("num_str"),
      nb::arg("update"),
      R"""(Routine to change a variable in the model lattice.

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
      nb::arg("branch_str"),
      nb::arg("mask_str"),
      nb::arg("print_list"),
      nb::arg("dqa_str"),
      nb::arg("dqb_str"),
      R"""(No docstring available.

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
      nb::arg("name"),
      nb::arg("num_str"),
      nb::arg("silent"),
      R"""(Routine to change a variable in the model lattice.

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
      nb::arg("branch_str"),
      nb::arg("dq_str"),
      R"""(No docstring available.

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
      nb::arg("data_type"),
      nb::arg("data_source"),
      R"""(Wrapper for Fortran routine tao_chrom_calc_needed

Parameters
----------
data_type : str

data_source : str

Returns
-------
do_chrom : bool
)"""
  );
  m.def(
      "tao_clear_cmd",
      &Tao::tao_clear_cmd,
      nb::arg("cmd_line"),
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
      nb::arg("gang"),
      nb::arg("where"),
      nb::arg("value1"),
      nb::arg("value2"),
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
      R"""(Wrapper for Fortran routine tao_close_command_file
)"""
  );
  m.def(
      "tao_cmd_history_record",
      &Tao::tao_cmd_history_record,
      nb::arg("cmd"),
      R"""(Subroutine to record a cmd in the command history stack
)"""
  );
  m.def(
      "tao_cmd_split",
      &Tao::tao_cmd_split,
      nb::arg("cmd_line"),
      nb::arg("n_word"),
      nb::arg("cmd_word"),
      nb::arg("extra_words_is_error"),
      nb::arg("separator") = nb::none(),
      R"""(This routine splits the command line into "words" (everything between separators).

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
      nb::arg("command_line"),
      nb::arg("err"),
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
      nb::arg("datum"),
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
      nb::arg("ele"),
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
      nb::arg("string"),
      nb::arg("pattern"),
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
      R"""(Subroutine to create the plot window.
This soubroutine knows not to create a second window if one already exists.
)"""
  );
  m.def(
      "tao_curve_beam_ellipse_setup",
      &Tao::tao_curve_beam_ellipse_setup,
      nb::arg("curve"),
      R"""(Wrapper for Fortran routine tao_curve_beam_ellipse_setup

Parameters
----------
curve : TaoCurveStruct
)"""
  );
  m.def(
      "tao_curve_check_universe",
      &Tao::tao_curve_check_universe,
      nb::arg("curve"),
      nb::arg("uni"),
      R"""(Routine to check if the universe associated with a curve exists and is on.

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
      nb::arg("plot"),
      nb::arg("graph"),
      nb::arg("curve"),
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
      nb::arg("eles"),
      nb::arg("plot"),
      nb::arg("curve"),
      nb::arg("who"),
      R"""(Routine to calculate datum values.
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
      nb::arg("curve"),
      nb::arg("point_to_ele_ref"),
      R"""(Wrapper for Fortran routine tao_curve_ele_ref

Parameters
----------
curve : TaoCurveStruct
    Curve with ref ele.

point_to_ele_ref : bool

Returns
-------
ele_track : EleStruct
)"""
  );
  m.def(
      "tao_curve_ix_uni",
      &Tao::tao_curve_ix_uni,
      nb::arg("curve"),
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
      nb::arg("curve"),
      nb::arg("use_region") = nb::none(),
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
  nb::class_<Tao::TaoCurveRmsCalc>(m, "TaoCurveRmsCalc", "tao_curve_rms_calc return type")
      .def_ro("rms", &Tao::TaoCurveRmsCalc::rms)
      .def_ro("mean", &Tao::TaoCurveRmsCalc::mean)
      .def("__len__", [](const Tao::TaoCurveRmsCalc &) { return 2; })
      .def("__getitem__", [](const Tao::TaoCurveRmsCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.rms);
        if (i == 1)
          return nb::cast(s.mean);
        throw nb::index_error();
      });
  m.def(
      "tao_curve_rms_calc",
      &Tao::tao_curve_rms_calc,
      nb::arg("curve"),
      nb::arg("who"),
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
      nb::arg("d1"),
      nb::arg("show_universe") = nb::none(),
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
      nb::arg("u"),
      nb::arg("d2_name"),
      nb::arg("n_d1_data"),
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
      nb::arg("err"),
      R"""(Wrapper for Fortran routine tao_data_check

Parameters
----------
err : bool
)"""
  );
  m.def(
      "tao_data_coupling_init",
      &Tao::tao_data_coupling_init,
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine tao_data_coupling_init

Parameters
----------
branch : BranchStruct
    New lattice branch.
)"""
  );
  m.def(
      "tao_data_sanity_check",
      [](TaoDataStruct &datum,
         bool print_err,
         std::string default_data_type,
         TaoUniverseStruct *uni) {
        auto fn = static_cast<
            bool (*)(TaoDataStruct &, bool, std::string, optional_ref<TaoUniverseStruct>)>(
            &Tao::tao_data_sanity_check
        );
        return fn(datum, print_err, default_data_type, ptr_to_opt_ref(uni));
      },
      nb::arg("datum"),
      nb::arg("print_err"),
      nb::arg("default_data_type"),
      nb::arg("uni") = nb::none(),
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
      [](TaoD2DataStruct &d2_data, CharacterAlloc1D *lines, std::optional<int> nl) {
        auto fn = static_cast<
            void (*)(TaoD2DataStruct &, optional_ref<CharacterAlloc1D>, std::optional<int>)>(
            &Tao::tao_data_show_use
        );
        return fn(d2_data, ptr_to_opt_ref(lines), nl);
      },
      nb::arg("d2_data"),
      nb::arg("lines") = nb::none(),
      nb::arg("nl") = nb::none(),
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
      nb::arg("template_"),
      nb::arg("curve"),
      nb::arg("graph"),
      R"""(Routine substitute the appropriate data type string for instances of "#ref" and
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
      nb::arg("curve"),
      nb::arg("graph"),
      nb::arg("data"),
      nb::arg("check_s_position"),
      R"""(Routine to set the data for plotting.

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
      nb::arg("data_type"),
      nb::arg("branch_geometry") = nb::none(),
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
  nb::class_<Tao::TaoDatumIntegrate>(m, "TaoDatumIntegrate", "tao_datum_integrate return type")
      .def_ro("valid_value", &Tao::TaoDatumIntegrate::valid_value)
      .def_ro("why_invalid", &Tao::TaoDatumIntegrate::why_invalid)
      .def_ro("result", &Tao::TaoDatumIntegrate::result)
      .def("__len__", [](const Tao::TaoDatumIntegrate &) { return 3; })
      .def("__getitem__", [](const Tao::TaoDatumIntegrate &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.valid_value);
        if (i == 1)
          return nb::cast(s.why_invalid);
        if (i == 2)
          return nb::cast(s.result);
        throw nb::index_error();
      });
  m.def(
      "tao_datum_integrate",
      &Tao::tao_datum_integrate,
      nb::arg("datum"),
      nb::arg("branch"),
      nb::arg("s_pos"),
      nb::arg("values"),
      R"""(Routine to calculate the integral, rms, or average of an array of values associated with a datum.

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
      nb::arg("datum"),
      nb::arg("show_universe") = nb::none(),
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
      nb::arg("datum"),
      nb::arg("ele"),
      R"""(Routine to calculate the longitudinal position associated with a datum.

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
      nb::arg("plot_cache"),
      R"""(Wrapper for Fortran routine tao_deallocate_plot_cache

Parameters
----------
plot_cache : 1D array of TaoPlotCacheStruct
)"""
  );
  m.def(
      "tao_deallocate_tree",
      &Tao::tao_deallocate_tree,
      nb::arg("tree"),
      R"""(Routine to deallocate tree%node(:) and everything below it

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
      R"""(Wrapper for Fortran routine tao_destroy_plot_window
)"""
  );
  m.def(
      "tao_dmerit_calc",
      &Tao::tao_dmerit_calc,
      R"""(No docstring available.
)"""
  );
  m.def(
      "tao_dmodel_dvar_calc",
      &Tao::tao_dmodel_dvar_calc,
      nb::arg("force_calc"),
      R"""(Subroutine to calculate the dModel_dVar derivative matrix.

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
      nb::arg("ele"),
      nb::arg("theta"),
      nb::arg("beam"),
      R"""(Returns the beam's second moment using the wire along the specified angle.
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
      nb::arg("plot"),
      nb::arg("graph"),
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
  nb::class_<PyTaoDrawCurveData>(m, "TaoDrawCurveData", "tao_draw_curve_data return type")
      .def_ro("have_data", &PyTaoDrawCurveData::have_data)
      .def("__len__", [](const PyTaoDrawCurveData &) { return 1; })
      .def("__getitem__", [](const PyTaoDrawCurveData &s, int i) -> nb::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return nb::cast(s.have_data);
        throw nb::index_error();
      });
  m.def(
      "tao_draw_curve_data",
      &python_tao_draw_curve_data,
      nb::arg("plot"),
      nb::arg("graph"),
      nb::arg("curve"),
      nb::arg("have_data"),
      R"""(Routine to draw a graph with data and/or variable curves.

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
      nb::arg("plot"),
      nb::arg("graph"),
      nb::arg("tao_lat"),
      nb::arg("ele"),
      nb::arg("ele_shape"),
      nb::arg("label_name"),
      nb::arg("offset1"),
      nb::arg("offset2"),
      R"""(Routine to draw one lattice element or one datum location for the floor plan graph.

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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw a floor plan graph.

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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw a just the graph part of a data graph.
The calling routine takes care of drawing any curves.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  nb::class_<PyTaoDrawHistogramData>(
      m,
      "TaoDrawHistogramData",
      "tao_draw_histogram_data return type"
  )
      .def_ro("have_data", &PyTaoDrawHistogramData::have_data)
      .def("__len__", [](const PyTaoDrawHistogramData &) { return 1; })
      .def("__getitem__", [](const PyTaoDrawHistogramData &s, int i) -> nb::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return nb::cast(s.have_data);
        throw nb::index_error();
      });
  m.def(
      "tao_draw_histogram_data",
      &python_tao_draw_histogram_data,
      nb::arg("plot"),
      nb::arg("graph"),
      nb::arg("curve"),
      nb::arg("have_data"),
      R"""(Routine to draw a graph with data and/or variable histograms.

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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw a lattice layout graph.

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
      nb::arg("do_clear") = nb::none(),
      R"""(Subroutine to draw the plots on the plot window.

Parameters
----------
do_clear : bool, optional
    If present and False then call qp_clear_page. This argument is used when drawing PS or GIF.
)"""
  );
  nb::class_<Tao::TaoEleGeometryWithMisalignments>(
      m,
      "TaoEleGeometryWithMisalignments",
      "tao_ele_geometry_with_misalignments return type"
  )
      .def_ro("valid_value", &Tao::TaoEleGeometryWithMisalignments::valid_value)
      .def_ro("why_invalid", &Tao::TaoEleGeometryWithMisalignments::why_invalid)
      .def_ro("value", &Tao::TaoEleGeometryWithMisalignments::value)
      .def("__len__", [](const Tao::TaoEleGeometryWithMisalignments &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEleGeometryWithMisalignments &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.valid_value);
        if (i == 1)
          return nb::cast(s.why_invalid);
        if (i == 2)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "tao_ele_geometry_with_misalignments",
      &Tao::tao_ele_geometry_with_misalignments,
      nb::arg("datum"),
      nb::arg("ele"),
      R"""(Routine to evaluate a floor position with misalignments at a given element.
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
  nb::class_<PyTaoEleShapeInfo>(m, "TaoEleShapeInfo", "tao_ele_shape_info return type")
      .def_ro("e_shape", &PyTaoEleShapeInfo::e_shape)
      .def_ro("label_name", &PyTaoEleShapeInfo::label_name)
      .def_ro("y1", &PyTaoEleShapeInfo::y1)
      .def_ro("y2", &PyTaoEleShapeInfo::y2)
      .def_ro("ix_shape_min", &PyTaoEleShapeInfo::ix_shape_min)
      .def("__len__", [](const PyTaoEleShapeInfo &) { return 5; })
      .def("__getitem__", [](const PyTaoEleShapeInfo &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.e_shape);
        if (i == 1)
          return nb::cast(s.label_name);
        if (i == 2)
          return nb::cast(s.y1);
        if (i == 3)
          return nb::cast(s.y2);
        if (i == 4)
          return nb::cast(s.ix_shape_min);
        throw nb::index_error();
      });
  m.def(
      "tao_ele_shape_info",
      &python_tao_ele_shape_info,
      nb::arg("ix_uni"),
      nb::arg("ele"),
      nb::arg("ele_shapes"),
      nb::arg("ix_shape_min") = nb::none(),
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
  m.def(
      "tao_ele_shape_input_to_struct",
      &Tao::tao_ele_shape_input_to_struct,
      nb::arg("shape_input"),
      nb::arg("namelist_name") = nb::none(),
      R"""(Wrapper for Fortran routine tao_ele_shape_input_to_struct

Parameters
----------
shape_input : TaoEleShapeInput

namelist_name : str, optional

Returns
-------
shape_struct : TaoEleShapeStruct
)"""
  );
  m.def(
      "tao_ele_shape_struct_to_input",
      &Tao::tao_ele_shape_struct_to_input,
      nb::arg("shape_struct"),
      R"""(Wrapper for Fortran routine tao_ele_shape_struct_to_input

Parameters
----------
shape_struct : TaoEleShapeStruct

Returns
-------
shape_input : TaoEleShapeInput
)"""
  );
  nb::class_<Tao::TaoEvalFloorOrbit>(m, "TaoEvalFloorOrbit", "tao_eval_floor_orbit return type")
      .def_ro("valid_value", &Tao::TaoEvalFloorOrbit::valid_value)
      .def_ro("why_invalid", &Tao::TaoEvalFloorOrbit::why_invalid)
      .def_ro("value", &Tao::TaoEvalFloorOrbit::value)
      .def("__len__", [](const Tao::TaoEvalFloorOrbit &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvalFloorOrbit &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.valid_value);
        if (i == 1)
          return nb::cast(s.why_invalid);
        if (i == 2)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "tao_eval_floor_orbit",
      &Tao::tao_eval_floor_orbit,
      nb::arg("datum"),
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("bunch_params"),
      R"""(Routine to evaluate a floor_orbit datum at a given element.
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
  nb::class_<Tao::TaoEvaluateADatum>(m, "TaoEvaluateADatum", "tao_evaluate_a_datum return type")
      .def_ro("datum_value", &Tao::TaoEvaluateADatum::datum_value)
      .def_ro("valid_value", &Tao::TaoEvaluateADatum::valid_value)
      .def_ro("why_invalid", &Tao::TaoEvaluateADatum::why_invalid)
      .def("__len__", [](const Tao::TaoEvaluateADatum &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvaluateADatum &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.datum_value);
        if (i == 1)
          return nb::cast(s.valid_value);
        if (i == 2)
          return nb::cast(s.why_invalid);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_a_datum",
      &Tao::tao_evaluate_a_datum,
      nb::arg("datum"),
      nb::arg("u"),
      nb::arg("tao_lat"),
      nb::arg("called_from_lat_calc") = nb::none(),
      nb::arg("print_err") = nb::none(),
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
  nb::class_<Tao::TaoEvaluateDatumAtS>(
      m,
      "TaoEvaluateDatumAtS",
      "tao_evaluate_datum_at_s return type"
  )
      .def_ro("err_str", &Tao::TaoEvaluateDatumAtS::err_str)
      .def_ro("bad_datum", &Tao::TaoEvaluateDatumAtS::bad_datum)
      .def_ro("value", &Tao::TaoEvaluateDatumAtS::value)
      .def("__len__", [](const Tao::TaoEvaluateDatumAtS &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvaluateDatumAtS &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.err_str);
        if (i == 1)
          return nb::cast(s.bad_datum);
        if (i == 2)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_datum_at_s",
      &Tao::tao_evaluate_datum_at_s,
      nb::arg("datum"),
      nb::arg("tao_lat"),
      nb::arg("ele"),
      nb::arg("ele_ref"),
      nb::arg("valid_value"),
      R"""(Routine to evaluate a datum at a given s-position in the lattice

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
  nb::class_<Tao::TaoEvaluateElementParameters>(
      m,
      "TaoEvaluateElementParameters",
      "tao_evaluate_element_parameters return type"
  )
      .def_ro("err", &Tao::TaoEvaluateElementParameters::err)
      .def_ro("values", &Tao::TaoEvaluateElementParameters::values)
      .def_ro("info", &Tao::TaoEvaluateElementParameters::info)
      .def("__len__", [](const Tao::TaoEvaluateElementParameters &) { return 3; })
      .def("__getitem__", [](const Tao::TaoEvaluateElementParameters &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.values);
        if (i == 2)
          return nb::cast(s.info);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_element_parameters",
      [](std::string param_name,
         bool print_err,
         std::string dflt_source,
         EleStruct *dflt_ele,
         std::optional<std::string> dflt_component,
         std::optional<int> dflt_uni,
         std::optional<int> eval_point) {
        auto fn = static_cast<Tao::TaoEvaluateElementParameters (*)(
            std::string,
            bool,
            std::string,
            optional_ref<EleStruct>,
            std::optional<std::string>,
            std::optional<int>,
            std::optional<int>
        )>(&Tao::tao_evaluate_element_parameters);
        return fn(
            param_name,
            print_err,
            dflt_source,
            ptr_to_opt_ref(dflt_ele),
            dflt_component,
            dflt_uni,
            eval_point
        );
      },
      nb::arg("param_name"),
      nb::arg("print_err"),
      nb::arg("dflt_source"),
      nb::arg("dflt_ele") = nb::none(),
      nb::arg("dflt_component") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("eval_point") = nb::none(),
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
  nb::class_<Tao::TaoEvaluateExpression>(
      m,
      "TaoEvaluateExpression",
      "tao_evaluate_expression return type"
  )
      .def_ro("value", &Tao::TaoEvaluateExpression::value)
      .def_ro("err_flag", &Tao::TaoEvaluateExpression::err_flag)
      .def_ro("info", &Tao::TaoEvaluateExpression::info)
      .def_ro("stack", &Tao::TaoEvaluateExpression::stack)
      .def("__len__", [](const Tao::TaoEvaluateExpression &) { return 4; })
      .def("__getitem__", [](const Tao::TaoEvaluateExpression &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.info);
        if (i == 3)
          return nb::cast(s.stack);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_expression",
      [](std::string expression,
         int n_size,
         bool use_good_user,
         std::optional<bool> print_err,
         std::optional<std::string> dflt_component,
         std::optional<std::string> dflt_source,
         EleStruct *dflt_ele_ref,
         EleStruct *dflt_ele_start,
         EleStruct *dflt_ele,
         std::optional<std::string> dflt_dat_or_var_index,
         std::optional<int> dflt_uni,
         std::optional<int> dflt_eval_point,
         std::optional<double> dflt_s_offset,
         CoordStruct *dflt_orbit,
         TaoDataStruct *datum) {
        auto fn = static_cast<Tao::TaoEvaluateExpression (*)(
            std::string,
            int,
            bool,
            std::optional<bool>,
            std::optional<std::string>,
            std::optional<std::string>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            std::optional<std::string>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            optional_ref<CoordStruct>,
            optional_ref<TaoDataStruct>
        )>(&Tao::tao_evaluate_expression);
        return fn(
            expression,
            n_size,
            use_good_user,
            print_err,
            dflt_component,
            dflt_source,
            ptr_to_opt_ref(dflt_ele_ref),
            ptr_to_opt_ref(dflt_ele_start),
            ptr_to_opt_ref(dflt_ele),
            dflt_dat_or_var_index,
            dflt_uni,
            dflt_eval_point,
            dflt_s_offset,
            ptr_to_opt_ref(dflt_orbit),
            ptr_to_opt_ref(datum)
        );
      },
      nb::arg("expression"),
      nb::arg("n_size"),
      nb::arg("use_good_user"),
      nb::arg("print_err") = nb::none(),
      nb::arg("dflt_component") = nb::none(),
      nb::arg("dflt_source") = nb::none(),
      nb::arg("dflt_ele_ref") = nb::none(),
      nb::arg("dflt_ele_start") = nb::none(),
      nb::arg("dflt_ele") = nb::none(),
      nb::arg("dflt_dat_or_var_index") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("dflt_eval_point") = nb::none(),
      nb::arg("dflt_s_offset") = nb::none(),
      nb::arg("dflt_orbit") = nb::none(),
      nb::arg("datum") = nb::none(),
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
  nb::class_<Tao::TaoEvaluateExpressionNew>(
      m,
      "TaoEvaluateExpressionNew",
      "tao_evaluate_expression_new return type"
  )
      .def_ro("value", &Tao::TaoEvaluateExpressionNew::value)
      .def_ro("err_flag", &Tao::TaoEvaluateExpressionNew::err_flag)
      .def_ro("info", &Tao::TaoEvaluateExpressionNew::info)
      .def_ro("stack", &Tao::TaoEvaluateExpressionNew::stack)
      .def("__len__", [](const Tao::TaoEvaluateExpressionNew &) { return 4; })
      .def("__getitem__", [](const Tao::TaoEvaluateExpressionNew &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.info);
        if (i == 3)
          return nb::cast(s.stack);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_expression_new",
      [](std::string expression,
         int n_size,
         bool use_good_user,
         std::optional<bool> print_err,
         std::optional<std::string> dflt_component,
         std::optional<std::string> dflt_source,
         EleStruct *dflt_ele_ref,
         EleStruct *dflt_ele_start,
         EleStruct *dflt_ele,
         std::optional<std::string> dflt_dat_or_var_index,
         std::optional<int> dflt_uni,
         std::optional<int> dflt_eval_point,
         std::optional<double> dflt_s_offset,
         CoordStruct *dflt_orbit,
         TaoDataStruct *datum) {
        auto fn = static_cast<Tao::TaoEvaluateExpressionNew (*)(
            std::string,
            int,
            bool,
            std::optional<bool>,
            std::optional<std::string>,
            std::optional<std::string>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            std::optional<std::string>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            optional_ref<CoordStruct>,
            optional_ref<TaoDataStruct>
        )>(&Tao::tao_evaluate_expression_new);
        return fn(
            expression,
            n_size,
            use_good_user,
            print_err,
            dflt_component,
            dflt_source,
            ptr_to_opt_ref(dflt_ele_ref),
            ptr_to_opt_ref(dflt_ele_start),
            ptr_to_opt_ref(dflt_ele),
            dflt_dat_or_var_index,
            dflt_uni,
            dflt_eval_point,
            dflt_s_offset,
            ptr_to_opt_ref(dflt_orbit),
            ptr_to_opt_ref(datum)
        );
      },
      nb::arg("expression"),
      nb::arg("n_size"),
      nb::arg("use_good_user"),
      nb::arg("print_err") = nb::none(),
      nb::arg("dflt_component") = nb::none(),
      nb::arg("dflt_source") = nb::none(),
      nb::arg("dflt_ele_ref") = nb::none(),
      nb::arg("dflt_ele_start") = nb::none(),
      nb::arg("dflt_ele") = nb::none(),
      nb::arg("dflt_dat_or_var_index") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("dflt_eval_point") = nb::none(),
      nb::arg("dflt_s_offset") = nb::none(),
      nb::arg("dflt_orbit") = nb::none(),
      nb::arg("datum") = nb::none(),
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
  nb::class_<Tao::TaoEvaluateExpressionOld>(
      m,
      "TaoEvaluateExpressionOld",
      "tao_evaluate_expression_old return type"
  )
      .def_ro("value", &Tao::TaoEvaluateExpressionOld::value)
      .def_ro("err_flag", &Tao::TaoEvaluateExpressionOld::err_flag)
      .def_ro("info", &Tao::TaoEvaluateExpressionOld::info)
      .def_ro("stack", &Tao::TaoEvaluateExpressionOld::stack)
      .def("__len__", [](const Tao::TaoEvaluateExpressionOld &) { return 4; })
      .def("__getitem__", [](const Tao::TaoEvaluateExpressionOld &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.info);
        if (i == 3)
          return nb::cast(s.stack);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_expression_old",
      [](std::string expression,
         int n_size,
         bool use_good_user,
         std::optional<bool> print_err,
         std::optional<std::string> dflt_component,
         std::optional<std::string> dflt_source,
         EleStruct *dflt_ele_ref,
         EleStruct *dflt_ele_start,
         EleStruct *dflt_ele,
         std::optional<std::string> dflt_dat_or_var_index,
         std::optional<int> dflt_uni,
         std::optional<int> dflt_eval_point,
         std::optional<double> dflt_s_offset,
         CoordStruct *dflt_orbit,
         TaoDataStruct *datum) {
        auto fn = static_cast<Tao::TaoEvaluateExpressionOld (*)(
            std::string,
            int,
            bool,
            std::optional<bool>,
            std::optional<std::string>,
            std::optional<std::string>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            std::optional<std::string>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            optional_ref<CoordStruct>,
            optional_ref<TaoDataStruct>
        )>(&Tao::tao_evaluate_expression_old);
        return fn(
            expression,
            n_size,
            use_good_user,
            print_err,
            dflt_component,
            dflt_source,
            ptr_to_opt_ref(dflt_ele_ref),
            ptr_to_opt_ref(dflt_ele_start),
            ptr_to_opt_ref(dflt_ele),
            dflt_dat_or_var_index,
            dflt_uni,
            dflt_eval_point,
            dflt_s_offset,
            ptr_to_opt_ref(dflt_orbit),
            ptr_to_opt_ref(datum)
        );
      },
      nb::arg("expression"),
      nb::arg("n_size"),
      nb::arg("use_good_user"),
      nb::arg("print_err") = nb::none(),
      nb::arg("dflt_component") = nb::none(),
      nb::arg("dflt_source") = nb::none(),
      nb::arg("dflt_ele_ref") = nb::none(),
      nb::arg("dflt_ele_start") = nb::none(),
      nb::arg("dflt_ele") = nb::none(),
      nb::arg("dflt_dat_or_var_index") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("dflt_eval_point") = nb::none(),
      nb::arg("dflt_s_offset") = nb::none(),
      nb::arg("dflt_orbit") = nb::none(),
      nb::arg("datum") = nb::none(),
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
  nb::class_<Tao::TaoEvaluateLatOrBeamData>(
      m,
      "TaoEvaluateLatOrBeamData",
      "tao_evaluate_lat_or_beam_data return type"
  )
      .def_ro("err", &Tao::TaoEvaluateLatOrBeamData::err)
      .def_ro("values", &Tao::TaoEvaluateLatOrBeamData::values)
      .def("__len__", [](const Tao::TaoEvaluateLatOrBeamData &) { return 2; })
      .def("__getitem__", [](const Tao::TaoEvaluateLatOrBeamData &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.values);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_lat_or_beam_data",
      [](std::string data_name,
         bool print_err,
         std::string default_source,
         EleStruct *dflt_ele_ref,
         EleStruct *dflt_ele_start,
         EleStruct *dflt_ele,
         std::optional<std::string> dflt_component,
         std::optional<int> dflt_uni,
         std::optional<int> dflt_eval_point,
         std::optional<double> dflt_s_offset) {
        auto fn = static_cast<Tao::TaoEvaluateLatOrBeamData (*)(
            std::string,
            bool,
            std::string,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            std::optional<std::string>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>
        )>(&Tao::tao_evaluate_lat_or_beam_data);
        return fn(
            data_name,
            print_err,
            default_source,
            ptr_to_opt_ref(dflt_ele_ref),
            ptr_to_opt_ref(dflt_ele_start),
            ptr_to_opt_ref(dflt_ele),
            dflt_component,
            dflt_uni,
            dflt_eval_point,
            dflt_s_offset
        );
      },
      nb::arg("data_name"),
      nb::arg("print_err"),
      nb::arg("default_source"),
      nb::arg("dflt_ele_ref") = nb::none(),
      nb::arg("dflt_ele_start") = nb::none(),
      nb::arg("dflt_ele") = nb::none(),
      nb::arg("dflt_component") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("dflt_eval_point") = nb::none(),
      nb::arg("dflt_s_offset") = nb::none(),
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
  nb::class_<Tao::TaoEvaluateStackOld>(
      m,
      "TaoEvaluateStackOld",
      "tao_evaluate_stack_old return type"
  )
      .def_ro("value", &Tao::TaoEvaluateStackOld::value)
      .def_ro("err_flag", &Tao::TaoEvaluateStackOld::err_flag)
      .def("__len__", [](const Tao::TaoEvaluateStackOld &) { return 2; })
      .def("__getitem__", [](const Tao::TaoEvaluateStackOld &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_stack_old",
      &Tao::tao_evaluate_stack_old,
      nb::arg("stack"),
      nb::arg("n_size_in"),
      nb::arg("use_good_user"),
      nb::arg("print_err"),
      nb::arg("expression"),
      nb::arg("info_in") = nb::none(),
      R"""(Routine to evaluate an expression stack.

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
  nb::class_<Tao::TaoEvaluateTree>(m, "TaoEvaluateTree", "tao_evaluate_tree return type")
      .def_ro("value", &Tao::TaoEvaluateTree::value)
      .def_ro("err_flag", &Tao::TaoEvaluateTree::err_flag)
      .def("__len__", [](const Tao::TaoEvaluateTree &) { return 2; })
      .def("__getitem__", [](const Tao::TaoEvaluateTree &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "tao_evaluate_tree",
      &Tao::tao_evaluate_tree,
      nb::arg("tao_tree"),
      nb::arg("n_size"),
      nb::arg("use_good_user"),
      nb::arg("print_err"),
      nb::arg("expression"),
      nb::arg("info_in") = nb::none(),
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
      nb::arg("q_str"),
      nb::arg("q0"),
      nb::arg("delta_input"),
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
      [](std::string expression_in, EleStruct *eval_ele) {
        auto fn = static_cast<std::string (*)(std::string, optional_ref<EleStruct>)>(
            &Tao::tao_expression_hash_substitute
        );
        return fn(expression_in, ptr_to_opt_ref(eval_ele));
      },
      nb::arg("expression_in"),
      nb::arg("eval_ele") = nb::none(),
      R"""(Routine to, in the expression, substitute the evaluation lattice element name in place
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
      [](TaoEvalNodeStruct &tree,
         std::optional<bool> include_root,
         std::optional<int> n_node,
         TaoEvalNodeStruct *parent) {
        auto fn = static_cast<std::string (*)(
            TaoEvalNodeStruct &,
            std::optional<bool>,
            std::optional<int>,
            optional_ref<TaoEvalNodeStruct>
        )>(&Tao::tao_expression_tree_to_string);
        return fn(tree, include_root, n_node, ptr_to_opt_ref(parent));
      },
      nb::arg("tree"),
      nb::arg("include_root") = nb::none(),
      nb::arg("n_node") = nb::none(),
      nb::arg("parent") = nb::none(),
      R"""(Routine to convert an expression tree to a expression string.

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
  nb::class_<Tao::TaoFindData>(m, "TaoFindData", "tao_find_data return type")
      .def_ro("err", &Tao::TaoFindData::err)
      .def_ro("d2_array", &Tao::TaoFindData::d2_array)
      .def_ro("d1_array", &Tao::TaoFindData::d1_array)
      .def_ro("d_array", &Tao::TaoFindData::d_array)
      .def_ro("re_array", &Tao::TaoFindData::re_array)
      .def_ro("log_array", &Tao::TaoFindData::log_array)
      .def_ro("str_array", &Tao::TaoFindData::str_array)
      .def_ro("int_array", &Tao::TaoFindData::int_array)
      .def_ro("component", &Tao::TaoFindData::component)
      .def("__len__", [](const Tao::TaoFindData &) { return 9; })
      .def("__getitem__", [](const Tao::TaoFindData &s, int i) -> nb::object {
        if (i < 0)
          i += 9;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.d2_array);
        if (i == 2)
          return nb::cast(s.d1_array);
        if (i == 3)
          return nb::cast(s.d_array);
        if (i == 4)
          return nb::cast(s.re_array);
        if (i == 5)
          return nb::cast(s.log_array);
        if (i == 6)
          return nb::cast(s.str_array);
        if (i == 7)
          return nb::cast(s.int_array);
        if (i == 8)
          return nb::cast(s.component);
        throw nb::index_error();
      });
  m.def(
      "tao_find_data",
      &Tao::tao_find_data,
      nb::arg("data_name"),
      nb::arg("ix_uni") = nb::none(),
      nb::arg("dflt_index") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Wrapper for Fortran routine tao_find_data

Parameters
----------
data_name : str
    The data name type. Eg: "3@orbit.x[2:5,10]|meas"

ix_uni : int, optional
    Index of default universe to use. If ix_uni = 0 then "viewed" universe will be used. Also, if not present
    then the "viewed" universe will be used.

dflt_index : str, optional
    If present and non-negative, and if no index is specified by the data_name argument, this index is used in
    the evaluation.

print_err : bool, optional
    Print error message if data is not found? Default is True.

Returns
-------
err : bool
    Err condition

d2_array : 1D array of TaoD2DataArrayStruct, optional
    Array of pointers to all the matching d2_data structure. Size(d2_array) = 0 if no structures found.

d1_array : 1D array of TaoD1DataArrayStruct, optional
    Array of pointers to all the matching d1_data structures. Size(d1_array) = 0 if no structures found.

d_array : 1D array of TaoDataArrayStruct, optional
    Array of pointers to all the matching tao_data_structs.  Size(d_array) = 0 if no structures found.

re_array : 1D array of TaoRealPointerStruct, optional
    Array of pointers to real component values.  Size(re_array) = 0 if no structures found.

log_array : 1D array of TaoLogicalArrayStruct, optional
    Array of pointers to logical component values.  Size(log_array) = 0 if no structures found.

str_array : 1D array of TaoStringArrayStruct, optional
    Array of pointers to character component values.  Size(str_array) = 0 if no structures found.

int_array : 1D array of TaoIntegerArrayStruct, optional
    Array of pointers to integer component values.  Size(int_array) = 0 if no structures found.

component : str, optional
    Name of the component. E.G: 'good_user' set to ' ' if no component present.
)"""
  );
  nb::class_<Tao::TaoFindPlotRegion>(m, "TaoFindPlotRegion", "tao_find_plot_region return type")
      .def_ro("err", &Tao::TaoFindPlotRegion::err)
      .def_ro("region", &Tao::TaoFindPlotRegion::region)
      .def("__len__", [](const Tao::TaoFindPlotRegion &) { return 2; })
      .def("__getitem__", [](const Tao::TaoFindPlotRegion &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.region);
        throw nb::index_error();
      });
  m.def(
      "tao_find_plot_region",
      &Tao::tao_find_plot_region,
      nb::arg("where"),
      nb::arg("print_flag") = nb::none(),
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
  nb::class_<Tao::TaoFindPlots>(m, "TaoFindPlots", "tao_find_plots return type")
      .def_ro("err", &Tao::TaoFindPlots::err)
      .def_ro("plot", &Tao::TaoFindPlots::plot)
      .def_ro("graph", &Tao::TaoFindPlots::graph)
      .def_ro("curve", &Tao::TaoFindPlots::curve)
      .def("__len__", [](const Tao::TaoFindPlots &) { return 4; })
      .def("__getitem__", [](const Tao::TaoFindPlots &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.plot);
        if (i == 2)
          return nb::cast(s.graph);
        if (i == 3)
          return nb::cast(s.curve);
        throw nb::index_error();
      });
  m.def(
      "tao_find_plots",
      &Tao::tao_find_plots,
      nb::arg("name"),
      nb::arg("where"),
      nb::arg("print_flag") = nb::none(),
      nb::arg("blank_means_all") = nb::none(),
      nb::arg("only_visible") = nb::none(),
      R"""(Wrapper for Fortran routine tao_find_plots

Parameters
----------
name : str
    Name of plot or region.

where : str
    Where to look: 'TEMPLATE', 'REGION', 'BOTH' For where = 'BOTH', if something is found in a plot region,
    then the templates will not be searched

print_flag : bool, optional
    If present and False then surpress error messages. Default is True.

blank_means_all : bool, optional
    If present and True then blank graph or curve fields get  interpreted as "*".

only_visible : bool, optional
    Default is True. If True and s.global.plot_on = True then only include visible plots. If False then plot
    visible setting is ignored.

Returns
-------
err : bool
    Set True on error. False otherwise.

plot : 1D array of TaoPlotArrayStruct, optional
    Array of plots. If error => size set to 0.

graph : 1D array of TaoGraphArrayStruct, optional
    Array of graphs. If error => size set to 0.

curve : 1D array of TaoCurveArrayStruct, optional
    Array of curves. If error => size set to 0.
)"""
  );
  nb::class_<Tao::TaoFindVar>(m, "TaoFindVar", "tao_find_var return type")
      .def_ro("err", &Tao::TaoFindVar::err)
      .def_ro("v1_array", &Tao::TaoFindVar::v1_array)
      .def_ro("v_array", &Tao::TaoFindVar::v_array)
      .def_ro("re_array", &Tao::TaoFindVar::re_array)
      .def_ro("log_array", &Tao::TaoFindVar::log_array)
      .def_ro("str_array", &Tao::TaoFindVar::str_array)
      .def_ro("component", &Tao::TaoFindVar::component)
      .def("__len__", [](const Tao::TaoFindVar &) { return 7; })
      .def("__getitem__", [](const Tao::TaoFindVar &s, int i) -> nb::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.v1_array);
        if (i == 2)
          return nb::cast(s.v_array);
        if (i == 3)
          return nb::cast(s.re_array);
        if (i == 4)
          return nb::cast(s.log_array);
        if (i == 5)
          return nb::cast(s.str_array);
        if (i == 6)
          return nb::cast(s.component);
        throw nb::index_error();
      });
  m.def(
      "tao_find_var",
      &Tao::tao_find_var,
      nb::arg("var_name"),
      nb::arg("print_err") = nb::none(),
      nb::arg("dflt_var_index") = nb::none(),
      R"""(Wrapper for Fortran routine tao_find_var

Parameters
----------
var_name : str
    Name of the variable.

print_err : bool, optional
    Print error message if data is not found? Default is True.

dflt_var_index : str, optional
    If present and "[...]" var selection substring is not present, then dflt_var_index will be used. [Do not
    include the brackets in this string.]

Returns
-------
err : bool
    err condition

v1_array : 1D array of TaoV1VarArrayStruct, optional
    Array of pointers to all the v1_var structures.

v_array : 1D array of TaoVarArrayStruct, optional
    Array of pointers to the variable data point.

re_array : 1D array of TaoRealPointerStruct, optional
    Array of pointers to the real component values.

log_array : 1D array of TaoLogicalArrayStruct, optional
    Array of pointers to logical component values.

str_array : 1D array of TaoStringArrayStruct, optional
    Array of pointers to character component values.

component : str, optional
    Name of the component. E.G: 'good_user' set to ' ' if no component present.
)"""
  );
  m.def(
      "tao_fixer",
      &Tao::tao_fixer,
      nb::arg("switch_"),
      nb::arg("word1"),
      nb::arg("word2"),
      R"""(Wrapper for Fortran routine tao_fixer

Parameters
----------
word1 : str
    First word of command.

word2 : str
    Secton word of command.
)"""
  );
  nb::class_<Tao::TaoFloorToScreen>(m, "TaoFloorToScreen", "tao_floor_to_screen return type")
      .def_ro("x_screen", &Tao::TaoFloorToScreen::x_screen)
      .def_ro("y_screen", &Tao::TaoFloorToScreen::y_screen)
      .def("__len__", [](const Tao::TaoFloorToScreen &) { return 2; })
      .def("__getitem__", [](const Tao::TaoFloorToScreen &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.x_screen);
        if (i == 1)
          return nb::cast(s.y_screen);
        throw nb::index_error();
      });
  m.def(
      "tao_floor_to_screen",
      &Tao::tao_floor_to_screen,
      nb::arg("graph"),
      nb::arg("r_floor"),
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
      nb::arg("graph"),
      nb::arg("floor"),
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
      R"""(Routine to minimize the merit function by varying variables until
the "data" as calculated from the model matches the measured data.

This subroutine is a wrapper for the "geodesic"
Levenburg - Marquardt method.

Returns
-------
abort : bool
    Set True if an user stop signal detected.
)"""
  );
  nb::class_<Tao::TaoGetData>(m, "TaoGetData", "tao_get_data return type")
      .def_ro("data_value", &Tao::TaoGetData::data_value)
      .def_ro("data_weight", &Tao::TaoGetData::data_weight)
      .def_ro("data_meas_value", &Tao::TaoGetData::data_meas_value)
      .def_ro("data_ix_dModel", &Tao::TaoGetData::data_ix_dModel)
      .def("__len__", [](const Tao::TaoGetData &) { return 4; })
      .def("__getitem__", [](const Tao::TaoGetData &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.data_value);
        if (i == 1)
          return nb::cast(s.data_weight);
        if (i == 2)
          return nb::cast(s.data_meas_value);
        if (i == 3)
          return nb::cast(s.data_ix_dModel);
        throw nb::index_error();
      });
  m.def(
      "tao_get_data",
      &Tao::tao_get_data,
      R"""(Subroutine to get the values of the data used in optimization and put them
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
  nb::class_<Tao::TaoGetOptVars>(m, "TaoGetOptVars", "tao_get_opt_vars return type")
      .def_ro("var_value", &Tao::TaoGetOptVars::var_value)
      .def_ro("var_step", &Tao::TaoGetOptVars::var_step)
      .def_ro("var_delta", &Tao::TaoGetOptVars::var_delta)
      .def_ro("var_weight", &Tao::TaoGetOptVars::var_weight)
      .def_ro("var_ix", &Tao::TaoGetOptVars::var_ix)
      .def_ro("ignore_if_weight_is_zero", &Tao::TaoGetOptVars::ignore_if_weight_is_zero)
      .def_ro("ignore_if_not_limited", &Tao::TaoGetOptVars::ignore_if_not_limited)
      .def("__len__", [](const Tao::TaoGetOptVars &) { return 7; })
      .def("__getitem__", [](const Tao::TaoGetOptVars &s, int i) -> nb::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return nb::cast(s.var_value);
        if (i == 1)
          return nb::cast(s.var_step);
        if (i == 2)
          return nb::cast(s.var_delta);
        if (i == 3)
          return nb::cast(s.var_weight);
        if (i == 4)
          return nb::cast(s.var_ix);
        if (i == 5)
          return nb::cast(s.ignore_if_weight_is_zero);
        if (i == 6)
          return nb::cast(s.ignore_if_not_limited);
        throw nb::index_error();
      });
  m.def(
      "tao_get_opt_vars",
      &Tao::tao_get_opt_vars,
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
      nb::arg("prompt_str") = nb::none(),
      nb::arg("wait_flag") = nb::none(),
      nb::arg("cmd_in") = nb::none(),
      R"""(Subroutine to get the next Tao command. In order of precedence, input may come from:
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
      nb::arg("graph"),
      R"""(Wrapper for Fortran routine tao_graph_controller_setup

Parameters
----------
graph : TaoGraphStruct
)"""
  );
  m.def(
      "tao_graph_data_setup",
      &Tao::tao_graph_data_setup,
      nb::arg("plot"),
      nb::arg("graph"),
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
      nb::arg("plot"),
      nb::arg("graph"),
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
      nb::arg("plot"),
      nb::arg("graph"),
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
      nb::arg("plot"),
      nb::arg("graph"),
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
      nb::arg("graph"),
      nb::arg("use_region") = nb::none(),
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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Wrapper for Fortran routine tao_graph_phase_space_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  nb::class_<Tao::TaoGraphSMinMaxCalc>(
      m,
      "TaoGraphSMinMaxCalc",
      "tao_graph_s_min_max_calc return type"
  )
      .def_ro("s_min", &Tao::TaoGraphSMinMaxCalc::s_min)
      .def_ro("s_max", &Tao::TaoGraphSMinMaxCalc::s_max)
      .def("__len__", [](const Tao::TaoGraphSMinMaxCalc &) { return 2; })
      .def("__getitem__", [](const Tao::TaoGraphSMinMaxCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.s_min);
        if (i == 1)
          return nb::cast(s.s_max);
        throw nb::index_error();
      });
  m.def(
      "tao_graph_s_min_max_calc",
      &Tao::tao_graph_s_min_max_calc,
      nb::arg("graph"),
      nb::arg("branch"),
      R"""(Routine to calculate min and max for a graph when plot%x_axis_type is set to "s".

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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Wrapper for Fortran routine tao_graph_setup

Parameters
----------
plot : TaoPlotStruct

graph : TaoGraphStruct
)"""
  );
  nb::class_<Tao::TaoHelp>(m, "TaoHelp", "tao_help return type")
      .def_ro("lines", &Tao::TaoHelp::lines)
      .def_ro("n_lines", &Tao::TaoHelp::n_lines)
      .def("__len__", [](const Tao::TaoHelp &) { return 2; })
      .def("__getitem__", [](const Tao::TaoHelp &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.lines);
        if (i == 1)
          return nb::cast(s.n_lines);
        throw nb::index_error();
      });
  m.def(
      "tao_help",
      &Tao::tao_help,
      nb::arg("what1"),
      nb::arg("what2"),
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
      nb::arg("u"),
      nb::arg("beam_init"),
      nb::arg("track_start"),
      nb::arg("track_end"),
      nb::arg("comb_ds_save"),
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
      nb::arg("init_file"),
      R"""(Subroutine to initialize beam stuff.

Parameters
----------
init_file : str
    Tao initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  m.def(
      "tao_init_data",
      &Tao::tao_init_data,
      nb::arg("data_file"),
      R"""(Subroutine to initialize the tao data structures.

Parameters
----------
data_file : str
    Tao data initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  m.def(
      "tao_init_data_end_stuff",
      &Tao::tao_init_data_end_stuff,
      R"""(Wrapper for Fortran routine tao_init_data_end_stuff
)"""
  );
  m.def(
      "tao_init_data_in_universe",
      &Tao::tao_init_data_in_universe,
      nb::arg("u"),
      nb::arg("n_d2_add"),
      nb::arg("keep_existing_data") = nb::none(),
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
      nb::arg("init_file"),
      R"""(Routine to initalize dynamic aperture simulations.

Parameters
----------
init_file : str
    File setting dynamic_aperture parameters.
)"""
  );
  nb::class_<Tao::TaoInitFindElements>(
      m,
      "TaoInitFindElements",
      "tao_init_find_elements return type"
  )
      .def_ro("eles", &Tao::TaoInitFindElements::eles)
      .def_ro("found_one", &Tao::TaoInitFindElements::found_one)
      .def("__len__", [](const Tao::TaoInitFindElements &) { return 2; })
      .def("__getitem__", [](const Tao::TaoInitFindElements &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.eles);
        if (i == 1)
          return nb::cast(s.found_one);
        throw nb::index_error();
      });
  m.def(
      "tao_init_find_elements",
      &Tao::tao_init_find_elements,
      nb::arg("u"),
      nb::arg("search_string"),
      nb::arg("attribute") = nb::none(),
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
      nb::arg("init_file"),
      R"""(Subroutine to initialize the tao global structures.

Parameters
----------
init_file : str
    Tao initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  m.def(
      "tao_init_lattice",
      &Tao::tao_init_lattice,
      nb::arg("lat_file"),
      nb::arg("err_flag"),
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
      nb::arg("plot_file"),
      R"""(Wrapper for Fortran routine tao_init_plotting

Parameters
----------
plot_file : str
)"""
  );
  m.def(
      "tao_init_variables",
      &Tao::tao_init_variables,
      nb::arg("var_file"),
      R"""(Subroutine to initialize the tao variable structures.

Parameters
----------
var_file : str
    Tao variable initialization file. If blank, there is no file so just use the defaults.
)"""
  );
  nb::class_<Tao::TaoInjectBeam>(m, "TaoInjectBeam", "tao_inject_beam return type")
      .def_ro("beam", &Tao::TaoInjectBeam::beam)
      .def_ro("init_ok", &Tao::TaoInjectBeam::init_ok)
      .def("__len__", [](const Tao::TaoInjectBeam &) { return 2; })
      .def("__getitem__", [](const Tao::TaoInjectBeam &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.beam);
        if (i == 1)
          return nb::cast(s.init_ok);
        throw nb::index_error();
      });
  m.def(
      "tao_inject_beam",
      &Tao::tao_inject_beam,
      nb::arg("u"),
      nb::arg("model"),
      nb::arg("ix_branch"),
      R"""(This will initialize the beam for a given lattice branch.

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
      nb::arg("u"),
      nb::arg("model"),
      nb::arg("ix_branch"),
      R"""(Wrapper for Fortran routine tao_inject_particle

Parameters
----------
u : TaoUniverseStruct

model : TaoLatticeStruct

ix_branch : int
)"""
  );
  nb::class_<Tao::TaoIsValidName>(m, "TaoIsValidName", "tao_is_valid_name return type")
      .def_ro("why_invalid", &Tao::TaoIsValidName::why_invalid)
      .def_ro("is_valid", &Tao::TaoIsValidName::is_valid)
      .def("__len__", [](const Tao::TaoIsValidName &) { return 2; })
      .def("__getitem__", [](const Tao::TaoIsValidName &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.why_invalid);
        if (i == 1)
          return nb::cast(s.is_valid);
        throw nb::index_error();
      });
  m.def(
      "tao_is_valid_name",
      &Tao::tao_is_valid_name,
      nb::arg("name"),
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
      nb::arg("input_str"),
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
      nb::arg("ix_key"),
      nb::arg("ix_min_key"),
      nb::arg("ix_max_key"),
      nb::arg("key_str"),
      nb::arg("header_str"),
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
      nb::arg("u"),
      nb::arg("tao_lat"),
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
      nb::arg("plane"),
      nb::arg("emit_type"),
      nb::arg("ele"),
      nb::arg("modes"),
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
      nb::arg("data_type"),
      nb::arg("data_source"),
      R"""(Wrapper for Fortran routine tao_lat_sigma_calc_needed

Parameters
----------
data_type : str

data_source : str

Returns
-------
do_lat_sigma : bool
)"""
  );
  m.def(
      "tao_lat_sigma_track",
      &Tao::tao_lat_sigma_track,
      nb::arg("tao_lat"),
      nb::arg("ix_branch"),
      nb::arg("print_err") = nb::none(),
      nb::arg("force_calc") = nb::none(),
      R"""(Routine to track the 6x6 sigma matrix through the lattice using the lattice linear transfer matrices.

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
      nb::arg("tlb1"),
      nb::arg("tlb2"),
      R"""(Wrapper for Fortran routine tao_lattice_branches_equal_tao_lattice_branches

Parameters
----------
tlb1 : 1D array of TaoLatticeBranchStruct

tlb2 : 1D array of TaoLatticeBranchStruct
)"""
  );
  nb::class_<Tao::TaoLatticeCalc>(m, "TaoLatticeCalc", "tao_lattice_calc return type")
      .def_ro("calc_ok", &Tao::TaoLatticeCalc::calc_ok)
      .def_ro("print_err", &Tao::TaoLatticeCalc::print_err)
      .def("__len__", [](const Tao::TaoLatticeCalc &) { return 2; })
      .def("__getitem__", [](const Tao::TaoLatticeCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.calc_ok);
        if (i == 1)
          return nb::cast(s.print_err);
        throw nb::index_error();
      });
  m.def(
      "tao_lattice_calc",
      &Tao::tao_lattice_calc,
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
      nb::arg("lat1"),
      nb::arg("lat2"),
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
      R"""(Routine to minimize the merit function by varying variables until
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
      R"""(Wrapper for Fortran routine tao_lmdif_optimizer

Returns
-------
abort : bool
    Set True if an user stop signal detected or there is a problem with calculating the merit function.
)"""
  );
  m.def(
      "tao_load_this_datum",
      [](FArray1D<Real> &vec,
         EleStruct &ele_ref,
         EleStruct &ele_start,
         EleStruct &ele,
         double datum_value,
         bool valid_value,
         TaoDataStruct &datum,
         BranchStruct &branch,
         std::optional<std::string> why_invalid,
         BoolAlloc1D *good) {
        auto fn = static_cast<void (*)(
            FArray1D<Real> &,
            EleStruct &,
            EleStruct &,
            EleStruct &,
            double,
            bool,
            TaoDataStruct &,
            BranchStruct &,
            std::optional<std::string>,
            optional_ref<BoolAlloc1D>
        )>(&Tao::tao_load_this_datum);
        return fn(
            vec,
            ele_ref,
            ele_start,
            ele,
            datum_value,
            valid_value,
            datum,
            branch,
            why_invalid,
            ptr_to_opt_ref(good)
        );
      },
      nb::arg("vec"),
      nb::arg("ele_ref"),
      nb::arg("ele_start"),
      nb::arg("ele"),
      nb::arg("datum_value"),
      nb::arg("valid_value"),
      nb::arg("datum"),
      nb::arg("branch"),
      nb::arg("why_invalid") = nb::none(),
      nb::arg("good") = nb::none(),
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
  nb::class_<Tao::TaoLocateAllElements>(
      m,
      "TaoLocateAllElements",
      "tao_locate_all_elements return type"
  )
      .def_ro("eles", &Tao::TaoLocateAllElements::eles)
      .def_ro("err", &Tao::TaoLocateAllElements::err)
      .def("__len__", [](const Tao::TaoLocateAllElements &) { return 2; })
      .def("__getitem__", [](const Tao::TaoLocateAllElements &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.eles);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "tao_locate_all_elements",
      &Tao::tao_locate_all_elements,
      nb::arg("ele_list"),
      nb::arg("ignore_blank") = nb::none(),
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
  nb::class_<Tao::TaoLocateElements>(m, "TaoLocateElements", "tao_locate_elements return type")
      .def_ro("eles", &Tao::TaoLocateElements::eles)
      .def_ro("err", &Tao::TaoLocateElements::err)
      .def("__len__", [](const Tao::TaoLocateElements &) { return 2; })
      .def("__getitem__", [](const Tao::TaoLocateElements &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.eles);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "tao_locate_elements",
      &Tao::tao_locate_elements,
      nb::arg("ele_list"),
      nb::arg("ix_universe"),
      nb::arg("lat_type") = nb::none(),
      nb::arg("ignore_blank") = nb::none(),
      nb::arg("err_stat_level") = nb::none(),
      nb::arg("above_ubound_is_err") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      nb::arg("multiple_eles_is_err") = nb::none(),
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
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine tao_mark_lattice_ele

Parameters
----------
lat : LatStruct
    Input lattice
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with elements marked.
)"""
  );
  nb::class_<Tao::TaoMerit>(m, "TaoMerit", "tao_merit return type")
      .def_ro("calc_ok", &Tao::TaoMerit::calc_ok)
      .def_ro("this_merit", &Tao::TaoMerit::this_merit)
      .def("__len__", [](const Tao::TaoMerit &) { return 2; })
      .def("__getitem__", [](const Tao::TaoMerit &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.calc_ok);
        if (i == 1)
          return nb::cast(s.this_merit);
        throw nb::index_error();
      });
  m.def(
      "tao_merit",
      &Tao::tao_merit,
      R"""(Wrapper for Fortran routine tao_merit

Returns
-------
this_merit : float
    Merit value.

calc_ok : bool, optional
    Set False if there was an error in the calculation like a particle was lost or a lat is unstable.
)"""
  );
  nb::class_<PyTaoNextSwitch>(m, "TaoNextSwitch", "tao_next_switch return type")
      .def_ro("switch_", &PyTaoNextSwitch::switch_)
      .def_ro("err", &PyTaoNextSwitch::err)
      .def_ro("line", &PyTaoNextSwitch::line)
      .def("__len__", [](const PyTaoNextSwitch &) { return 3; })
      .def("__getitem__", [](const PyTaoNextSwitch &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.switch_);
        if (i == 1)
          return nb::cast(s.err);
        if (i == 2)
          return nb::cast(s.line);
        throw nb::index_error();
      });
  m.def(
      "tao_next_switch",
      &python_tao_next_switch,
      nb::arg("line"),
      nb::arg("switch_list"),
      nb::arg("return_next_word"),
      nb::arg("neg_num_not_switch") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Subroutine look at the next word on the command line and match this word to a list of "switches"
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
  nb::class_<PyTaoNextWord>(m, "TaoNextWord", "tao_next_word return type")
      .def_ro("word", &PyTaoNextWord::word)
      .def_ro("line", &PyTaoNextWord::line)
      .def("__len__", [](const PyTaoNextWord &) { return 2; })
      .def("__getitem__", [](const PyTaoNextWord &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.word);
        if (i == 1)
          return nb::cast(s.line);
        throw nb::index_error();
      });
  m.def(
      "tao_next_word",
      &python_tao_next_word,
      nb::arg("line"),
      R"""(Routine to return the next word in a line.

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
      nb::arg("data_type"),
      nb::arg("data_source"),
      R"""(Wrapper for Fortran routine tao_one_turn_map_calc_needed

Parameters
----------
data_type : str

data_source : str

Returns
-------
do_one_turn_map : bool
)"""
  );
  m.def(
      "tao_open_file",
      &Tao::tao_open_file,
      nb::arg("file"),
      nb::arg("file_name"),
      nb::arg("error_severity"),
      nb::arg("binary") = nb::none(),
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
  nb::class_<Tao::TaoOpenScratchFile>(m, "TaoOpenScratchFile", "tao_open_scratch_file return type")
      .def_ro("err", &Tao::TaoOpenScratchFile::err)
      .def_ro("iu", &Tao::TaoOpenScratchFile::iu)
      .def("__len__", [](const Tao::TaoOpenScratchFile &) { return 2; })
      .def("__getitem__", [](const Tao::TaoOpenScratchFile &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.iu);
        throw nb::index_error();
      });
  m.def(
      "tao_open_scratch_file",
      &Tao::tao_open_scratch_file,
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
      nb::arg("datum"),
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
      nb::arg("plot"),
      R"""(No docstring available.
)"""
  );
  m.def(
      "tao_oreint_building_wall_pt",
      &Tao::tao_oreint_building_wall_pt,
      nb::arg("pt_in"),
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
  nb::class_<Tao::TaoParamValueAtS>(m, "TaoParamValueAtS", "tao_param_value_at_s return type")
      .def_ro("err_flag", &Tao::TaoParamValueAtS::err_flag)
      .def_ro("why_invalid", &Tao::TaoParamValueAtS::why_invalid)
      .def_ro("print_err", &Tao::TaoParamValueAtS::print_err)
      .def_ro("bad_datum", &Tao::TaoParamValueAtS::bad_datum)
      .def_ro("value", &Tao::TaoParamValueAtS::value)
      .def("__len__", [](const Tao::TaoParamValueAtS &) { return 5; })
      .def("__getitem__", [](const Tao::TaoParamValueAtS &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.why_invalid);
        if (i == 2)
          return nb::cast(s.print_err);
        if (i == 3)
          return nb::cast(s.bad_datum);
        if (i == 4)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "tao_param_value_at_s",
      &Tao::tao_param_value_at_s,
      nb::arg("dat_name"),
      nb::arg("ele_to_s"),
      nb::arg("ele_here"),
      nb::arg("orbit"),
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
      [](std::string str,
         bool use_good_user,
         std::string saved_prefix,
         TaoEvalNodeStruct &stack,
         bool err_flag,
         bool print_err,
         std::optional<std::string> dflt_component,
         std::optional<std::string> dflt_source,
         EleStruct *dflt_ele_ref,
         EleStruct *dflt_ele_start,
         EleStruct *dflt_ele,
         std::optional<std::string> dflt_dat_or_var_index,
         std::optional<int> dflt_uni,
         std::optional<int> dflt_eval_point,
         std::optional<double> dflt_s_offset,
         CoordStruct *dflt_orbit,
         TaoDataStruct *datum) {
        auto fn = static_cast<void (*)(
            std::string,
            bool,
            std::string,
            TaoEvalNodeStruct &,
            bool,
            bool,
            std::optional<std::string>,
            std::optional<std::string>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            optional_ref<EleStruct>,
            std::optional<std::string>,
            std::optional<int>,
            std::optional<int>,
            std::optional<double>,
            optional_ref<CoordStruct>,
            optional_ref<TaoDataStruct>
        )>(&Tao::tao_param_value_routine);
        return fn(
            str,
            use_good_user,
            saved_prefix,
            stack,
            err_flag,
            print_err,
            dflt_component,
            dflt_source,
            ptr_to_opt_ref(dflt_ele_ref),
            ptr_to_opt_ref(dflt_ele_start),
            ptr_to_opt_ref(dflt_ele),
            dflt_dat_or_var_index,
            dflt_uni,
            dflt_eval_point,
            dflt_s_offset,
            ptr_to_opt_ref(dflt_orbit),
            ptr_to_opt_ref(datum)
        );
      },
      nb::arg("str"),
      nb::arg("use_good_user"),
      nb::arg("saved_prefix"),
      nb::arg("stack"),
      nb::arg("err_flag"),
      nb::arg("print_err"),
      nb::arg("dflt_component") = nb::none(),
      nb::arg("dflt_source") = nb::none(),
      nb::arg("dflt_ele_ref") = nb::none(),
      nb::arg("dflt_ele_start") = nb::none(),
      nb::arg("dflt_ele") = nb::none(),
      nb::arg("dflt_dat_or_var_index") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("dflt_eval_point") = nb::none(),
      nb::arg("dflt_s_offset") = nb::none(),
      nb::arg("dflt_orbit") = nb::none(),
      nb::arg("datum") = nb::none(),
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
      nb::arg("cmd_line") = nb::none(),
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
  nb::class_<Tao::TaoParseElementParamStr>(
      m,
      "TaoParseElementParamStr",
      "tao_parse_element_param_str return type"
  )
      .def_ro("err", &Tao::TaoParseElementParamStr::err)
      .def_ro("uni", &Tao::TaoParseElementParamStr::uni)
      .def_ro("element", &Tao::TaoParseElementParamStr::element)
      .def_ro("parameter", &Tao::TaoParseElementParamStr::parameter)
      .def_ro("where", &Tao::TaoParseElementParamStr::where)
      .def_ro("component", &Tao::TaoParseElementParamStr::component)
      .def("__len__", [](const Tao::TaoParseElementParamStr &) { return 6; })
      .def("__getitem__", [](const Tao::TaoParseElementParamStr &s, int i) -> nb::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return nb::cast(s.err);
        if (i == 1)
          return nb::cast(s.uni);
        if (i == 2)
          return nb::cast(s.element);
        if (i == 3)
          return nb::cast(s.parameter);
        if (i == 4)
          return nb::cast(s.where);
        if (i == 5)
          return nb::cast(s.component);
        throw nb::index_error();
      });
  m.def(
      "tao_parse_element_param_str",
      &Tao::tao_parse_element_param_str,
      nb::arg("in_str"),
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
  nb::class_<Tao::TaoParticleDataValue>(
      m,
      "TaoParticleDataValue",
      "tao_particle_data_value return type"
  )
      .def_ro("value", &Tao::TaoParticleDataValue::value)
      .def_ro("err", &Tao::TaoParticleDataValue::err)
      .def("__len__", [](const Tao::TaoParticleDataValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoParticleDataValue &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "tao_particle_data_value",
      &Tao::tao_particle_data_value,
      nb::arg("data_type"),
      nb::arg("p"),
      nb::arg("ele"),
      nb::arg("ix_bunch"),
      R"""(Routine to calculate the value array of a data_type for an array of particles.

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
      nb::arg("time"),
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
      nb::arg("data_type"),
      nb::arg("err"),
      R"""(Routine to calculate the phase space axis index for a given data type.

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
      nb::arg("plot"),
      R"""(No docstring available.
)"""
  );
  nb::class_<Tao::TaoPickUniverse>(m, "TaoPickUniverse", "tao_pick_universe return type")
      .def_ro("name_out", &Tao::TaoPickUniverse::name_out)
      .def_ro("picked", &Tao::TaoPickUniverse::picked)
      .def_ro("err", &Tao::TaoPickUniverse::err)
      .def_ro("ix_uni", &Tao::TaoPickUniverse::ix_uni)
      .def_ro("explicit_uni", &Tao::TaoPickUniverse::explicit_uni)
      .def("__len__", [](const Tao::TaoPickUniverse &) { return 5; })
      .def("__getitem__", [](const Tao::TaoPickUniverse &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.name_out);
        if (i == 1)
          return nb::cast(s.picked);
        if (i == 2)
          return nb::cast(s.err);
        if (i == 3)
          return nb::cast(s.ix_uni);
        if (i == 4)
          return nb::cast(s.explicit_uni);
        throw nb::index_error();
      });
  m.def(
      "tao_pick_universe",
      &Tao::tao_pick_universe,
      nb::arg("name_in"),
      nb::arg("dflt_uni") = nb::none(),
      nb::arg("pure_uni") = nb::none(),
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
      nb::arg("input_str"),
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
      nb::arg("where"),
      nb::arg("who"),
      nb::arg("no_buffer") = nb::none(),
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
      nb::arg("where"),
      nb::arg("component"),
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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw a graph with data and/or variable curves.

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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw one graph for the histogram analysis plot.

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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw a key table graph.

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
      R"""(Wrapper for Fortran routine tao_plot_setup
)"""
  );
  m.def(
      "tao_plot_struct_transfer",
      &Tao::tao_plot_struct_transfer,
      nb::arg("plot_in"),
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
      nb::arg("plot"),
      nb::arg("graph"),
      R"""(Routine to draw one graph for the wave analysis plot.

Parameters
----------
plot : TaoPlotStruct
    Plot containing the graph.

graph : TaoGraphStruct
    Graph to plot.
)"""
  );
  nb::class_<Tao::TaoPointerToBranches>(
      m,
      "TaoPointerToBranches",
      "tao_pointer_to_branches return type"
  )
      .def_ro("branches", &Tao::TaoPointerToBranches::branches)
      .def_ro("unis", &Tao::TaoPointerToBranches::unis)
      .def_ro("err", &Tao::TaoPointerToBranches::err)
      .def("__len__", [](const Tao::TaoPointerToBranches &) { return 3; })
      .def("__getitem__", [](const Tao::TaoPointerToBranches &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.branches);
        if (i == 1)
          return nb::cast(s.unis);
        if (i == 2)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "tao_pointer_to_branches",
      &Tao::tao_pointer_to_branches,
      nb::arg("branch_str"),
      R"""(Wrapper for Fortran routine tao_pointer_to_branches

Parameters
----------
branch_str : str
    String specifying what branches to use.

Returns
-------
branches : 1D array of BranchPointerStruct
    Array of pointers to branches.

unis : 1D array of TaoUniversePointerStruct
    Array of corresponding universes.

err : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "tao_pointer_to_building_wall_shape",
      &Tao::tao_pointer_to_building_wall_shape,
      nb::arg("wall_name"),
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
      nb::arg("d1"),
      nb::arg("ele_name"),
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
  nb::class_<Tao::TaoPointerToDatumEle>(
      m,
      "TaoPointerToDatumEle",
      "tao_pointer_to_datum_ele return type"
  )
      .def_ro("valid", &Tao::TaoPointerToDatumEle::valid)
      .def_ro("why_invalid", &Tao::TaoPointerToDatumEle::why_invalid)
      .def_ro("ele", &Tao::TaoPointerToDatumEle::ele)
      .def("__len__", [](const Tao::TaoPointerToDatumEle &) { return 3; })
      .def("__getitem__", [](const Tao::TaoPointerToDatumEle &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.valid);
        if (i == 1)
          return nb::cast(s.why_invalid);
        if (i == 2)
          return nb::cast(s.ele);
        throw nb::index_error();
      });
  m.def(
      "tao_pointer_to_datum_ele",
      &Tao::tao_pointer_to_datum_ele,
      nb::arg("lat"),
      nb::arg("ele_name"),
      nb::arg("ix_ele"),
      nb::arg("datum"),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to see if an element index corresponds to an element with a definite
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
  nb::class_<PyTaoPointerToEleShape>(
      m,
      "TaoPointerToEleShape",
      "tao_pointer_to_ele_shape return type"
  )
      .def_ro("dat_var_name", &PyTaoPointerToEleShape::dat_var_name)
      .def_ro("dat_var_value", &PyTaoPointerToEleShape::dat_var_value)
      .def_ro("e_shape", &PyTaoPointerToEleShape::e_shape)
      .def_ro("ix_shape_min", &PyTaoPointerToEleShape::ix_shape_min)
      .def("__len__", [](const PyTaoPointerToEleShape &) { return 4; })
      .def("__getitem__", [](const PyTaoPointerToEleShape &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.dat_var_name);
        if (i == 1)
          return nb::cast(s.dat_var_value);
        if (i == 2)
          return nb::cast(s.e_shape);
        if (i == 3)
          return nb::cast(s.ix_shape_min);
        throw nb::index_error();
      });
  m.def(
      "tao_pointer_to_ele_shape",
      &python_tao_pointer_to_ele_shape,
      nb::arg("ix_uni"),
      nb::arg("ele"),
      nb::arg("ele_shape"),
      nb::arg("ix_shape_min") = nb::none(),
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
      nb::arg("u"),
      nb::arg("lat_type") = nb::none(),
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
      nb::overload_cast<int, std::optional<bool>>(&Tao::tao_pointer_to_universe),
      nb::arg("ix_uni"),
      nb::arg("neg2_to_default") = nb::none(),
      R"""(Routine to set a pointer to a universe.

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
  nb::class_<PyTaoPointerToUniverseStr>(
      m,
      "TaoPointerToUniverseStr",
      "tao_pointer_to_universe_str return type"
  )
      .def_ro("u", &PyTaoPointerToUniverseStr::u)
      .def_ro("string", &PyTaoPointerToUniverseStr::string)
      .def("__len__", [](const PyTaoPointerToUniverseStr &) { return 2; })
      .def("__getitem__", [](const PyTaoPointerToUniverseStr &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.u);
        if (i == 1)
          return nb::cast(s.string);
        throw nb::index_error();
      });
  m.def(
      "tao_pointer_to_universe",
      nb::overload_cast<std::string, std::optional<bool>>(&python_tao_pointer_to_universe_str),
      nb::arg("string"),
      nb::arg("neg2_to_default") = nb::none(),
      R"""(Routine to set a pointer to a universe.

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
  nb::class_<Tao::TaoPointerToUniverses>(
      m,
      "TaoPointerToUniverses",
      "tao_pointer_to_universes return type"
  )
      .def_ro("unis", &Tao::TaoPointerToUniverses::unis)
      .def_ro("err", &Tao::TaoPointerToUniverses::err)
      .def_ro("name_out", &Tao::TaoPointerToUniverses::name_out)
      .def_ro("explicit_uni", &Tao::TaoPointerToUniverses::explicit_uni)
      .def("__len__", [](const Tao::TaoPointerToUniverses &) { return 4; })
      .def("__getitem__", [](const Tao::TaoPointerToUniverses &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.unis);
        if (i == 1)
          return nb::cast(s.err);
        if (i == 2)
          return nb::cast(s.name_out);
        if (i == 3)
          return nb::cast(s.explicit_uni);
        throw nb::index_error();
      });
  m.def(
      "tao_pointer_to_universes",
      &Tao::tao_pointer_to_universes,
      nb::arg("name_in"),
      nb::arg("dflt_uni") = nb::none(),
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
      nb::arg("var"),
      nb::arg("ix_uni"),
      nb::arg("ele"),
      R"""(Routine to add a pointer from a given Tao variable
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
      nb::arg("var"),
      nb::arg("ix_uni"),
      R"""(Routine to add a pointer from a given Tao variable
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
      R"""(Wrapper for Fortran routine tao_print_command_line_info
)"""
  );
  m.def(
      "tao_print_vars",
      &Tao::tao_print_vars,
      nb::arg("iu"),
      nb::arg("ix_uni"),
      nb::arg("show_good_opt_only") = nb::none(),
      nb::arg("tao_format") = nb::none(),
      nb::arg("v_array") = nb::none(),
      R"""(Routine to print a list of set statements for the Bmad parameters controlled by the Tao variables.
The set statements are Bmad lattice format compatible.

When tao_format = True, the output is in the form "set variable <name> = <value>"
so the file can be used as a Tao command file. If tao_format = False, the format
is suitable for inclusion in a Bmad lattice file.

Parameters
----------
iu : int
    File unit number. 0 => print to the terminal.

ix_uni : int
    Universe index. If zero print slave parameters for all universes. If non-zero, only print set statements
    for slave parameters of this universe. Ignored if tao_format = True.

show_good_opt_only : bool, optional
    If True, only show slave parameters of variables used in optimization.

tao_format : bool, optional
    Output format. Default False. See above.

v_array : 1D array of TaoVarArrayStruct, optional
    Variable array. If present, restrict printing to parameters of these variables.
)"""
  );
  m.def(
      "tao_ptc_normal_form",
      &Tao::tao_ptc_normal_form,
      nb::arg("do_calc"),
      nb::arg("tao_lat"),
      nb::arg("ix_branch"),
      nb::arg("rf_on") = nb::none(),
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
      nb::arg("input_str"),
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
      nb::arg("set"),
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
      nb::arg("data_type"),
      nb::arg("data_source"),
      R"""(Wrapper for Fortran routine tao_rad_int_calc_needed

Parameters
----------
data_type : str

data_source : str

Returns
-------
do_rad_int : bool
)"""
  );
  m.def(
      "tao_re_allocate_expression_info",
      &Tao::tao_re_allocate_expression_info,
      nb::arg("info"),
      nb::arg("n"),
      nb::arg("exact") = nb::none(),
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
      nb::arg("tree"),
      nb::arg("n"),
      nb::arg("exact") = nb::none(),
      R"""(Routine to resize the tree%node(:) array.

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
      nb::arg("string"),
      nb::arg("err"),
      R"""(Subroutine to execute a previous command.
)"""
  );
  m.def(
      "tao_read_cmd",
      &Tao::tao_read_cmd,
      nb::arg("which"),
      nb::arg("unis"),
      nb::arg("file"),
      nb::arg("silent"),
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
      nb::arg("name"),
      nb::arg("ixc"),
      nb::arg("print_err") = nb::none(),
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
      nb::arg("cmd_str"),
      R"""(Wrapper for Fortran routine tao_regression_test

Parameters
----------
cmd_str : str
)"""
  );
  nb::class_<PyTaoRemoveBlankCharacters>(
      m,
      "TaoRemoveBlankCharacters",
      "tao_remove_blank_characters return type"
  )
      .def_ro("str", &PyTaoRemoveBlankCharacters::str)
      .def("__len__", [](const PyTaoRemoveBlankCharacters &) { return 1; })
      .def("__getitem__", [](const PyTaoRemoveBlankCharacters &s, int i) -> nb::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return nb::cast(s.str);
        throw nb::index_error();
      });
  m.def(
      "tao_remove_blank_characters",
      &python_tao_remove_blank_characters,
      nb::arg("str"),
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
      nb::arg("which"),
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
      nb::arg("where"),
      nb::arg("y_min_in"),
      nb::arg("y_max_in"),
      nb::arg("axis") = nb::none(),
      nb::arg("include_wall") = nb::none(),
      nb::arg("gang") = nb::none(),
      nb::arg("exact") = nb::none(),
      nb::arg("turn_autoscale_off") = nb::none(),
      R"""(Routine to scale a plot.
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
  nb::class_<Tao::TaoScaleGraph>(m, "TaoScaleGraph", "tao_scale_graph return type")
      .def_ro("y_range", &Tao::TaoScaleGraph::y_range)
      .def_ro("y2_range", &Tao::TaoScaleGraph::y2_range)
      .def("__len__", [](const Tao::TaoScaleGraph &) { return 2; })
      .def("__getitem__", [](const Tao::TaoScaleGraph &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.y_range);
        if (i == 1)
          return nb::cast(s.y2_range);
        throw nb::index_error();
      });
  m.def(
      "tao_scale_graph",
      &Tao::tao_scale_graph,
      nb::arg("graph"),
      nb::arg("y_min"),
      nb::arg("y_max"),
      nb::arg("axis") = nb::none(),
      nb::arg("include_wall") = nb::none(),
      R"""(Routine to scale the y-axis and/or y2-axis of a graph
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
      nb::arg("u"),
      R"""(Wrapper for Fortran routine tao_scale_ping_data

Parameters
----------
u : TaoUniverseStruct
)"""
  );
  m.def(
      "tao_scale_plot",
      &Tao::tao_scale_plot,
      nb::arg("plot"),
      nb::arg("y_min_in"),
      nb::arg("y_max_in"),
      nb::arg("axis") = nb::none(),
      nb::arg("include_wall") = nb::none(),
      nb::arg("gang") = nb::none(),
      nb::arg("skip_lat_layout") = nb::none(),
      R"""(Routine to scale the y-axis and/or y2-axis of the graphs of the plot.
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
      nb::arg("ele_ref"),
      nb::arg("ele_start"),
      nb::arg("ele"),
      nb::arg("datum"),
      nb::arg("branch"),
      nb::arg("orbit"),
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
      nb::arg("who"),
      nb::arg("value_str"),
      nb::arg("branch_str"),
      R"""(Routine to set various beam parameters.

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
      nb::arg("who"),
      nb::arg("value_str"),
      nb::arg("branch_str"),
      R"""(Routine to set beam_init variables

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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set bmad_com variables

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
      nb::arg("branch_str"),
      nb::arg("component_str"),
      nb::arg("value_str"),
      R"""(Routine to set lattice branch values.

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
      nb::arg("switch_") = nb::none(),
      R"""(Toggles off lattice calc and plotting.
)"""
  );
  m.def(
      "tao_set_curve_cmd",
      &Tao::tao_set_curve_cmd,
      nb::arg("curve_name"),
      nb::arg("component"),
      nb::arg("value_str"),
      R"""(Routine to set var values.

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
      nb::arg("curve"),
      nb::arg("why_invalid"),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to set curve%valid to False.

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
      nb::arg("who_str"),
      nb::arg("value_str"),
      nb::arg("silent") = nb::none(),
      R"""(Routine to set data values.

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
      nb::arg("data") = nb::none(),
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
      nb::arg("who_str"),
      nb::arg("value_str"),
      R"""(Routine to set default values.

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
      nb::arg("drawing"),
      nb::arg("component"),
      nb::arg("value_str"),
      R"""(Routine to set floor_plan and lat_layout parameters.

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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Sets dynamic aperture parameters.

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
      nb::arg("ele_list"),
      nb::arg("attribute"),
      nb::arg("value"),
      nb::arg("update"),
      R"""(Sets element parameters.

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
      "tao_set_flags_for_changed_attribute",
      [](TaoUniverseStruct &u,
         std::string ele_name,
         EleStruct *ele_ptr,
         AllPointerStruct *val_ptr,
         std::optional<std::string> who) {
        auto fn = static_cast<void (*)(
            TaoUniverseStruct &,
            std::string,
            optional_ref<EleStruct>,
            optional_ref<AllPointerStruct>,
            std::optional<std::string>
        )>(&Tao::tao_set_flags_for_changed_attribute);
        return fn(u, ele_name, ptr_to_opt_ref(ele_ptr), ptr_to_opt_ref(val_ptr), who);
      },
      nb::arg("u"),
      nb::arg("ele_name"),
      nb::arg("ele_ptr") = nb::none(),
      nb::arg("val_ptr") = nb::none(),
      nb::arg("who") = nb::none(),
      R"""(Wrapper for Fortran routine tao_set_flags_for_changed_attribute

Parameters
----------
u : TaoUniverseStruct
    Universe containing the lattice.

ele_name : str
    Associated "element" of the changed parameter.

ele_ptr : EleStruct, optional
    Pointer to the element. May be null, for example, if ele_name = "PARTICLE_START".

val_ptr : AllPointerStruct, optional
    optional: Pointer to the attribute that was changed. Must be present if ele_ptr is present.

who : str, optional
    Name of changed attribute. Only used with PARTICLE_START.
)"""
  );
  m.def(
      "tao_set_floor_plan_axis_label",
      &Tao::tao_set_floor_plan_axis_label,
      nb::arg("graph"),
      nb::arg("axis_in"),
      nb::arg("axis_out"),
      nb::arg("which"),
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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set geodesic_lm variables

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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set global variables

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
      nb::arg("graph_name"),
      nb::arg("component"),
      nb::arg("value_str"),
      R"""(Routine to set var values.

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
  nb::class_<Tao::TaoSetIntegerValue>(m, "TaoSetIntegerValue", "tao_set_integer_value return type")
      .def_ro("var", &Tao::TaoSetIntegerValue::var)
      .def_ro("error", &Tao::TaoSetIntegerValue::error)
      .def("__len__", [](const Tao::TaoSetIntegerValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetIntegerValue &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.var);
        if (i == 1)
          return nb::cast(s.error);
        throw nb::index_error();
      });
  m.def(
      "tao_set_integer_value",
      &Tao::tao_set_integer_value,
      nb::arg("var_str"),
      nb::arg("value_str"),
      nb::arg("min_val") = nb::none(),
      nb::arg("max_val") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Subroutine to read and set the value of an integer varialbe.

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
      nb::arg("datum"),
      nb::arg("message"),
      nb::arg("exterminate") = nb::none(),
      nb::arg("err_level") = nb::none(),
      nb::arg("print_err") = nb::none(),
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
      nb::arg("key_str"),
      nb::arg("cmd_str"),
      R"""(Associates a command with a key press for single mode.

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
      nb::arg("dest_lat"),
      nb::arg("source_lat"),
      R"""(Sets a lattice equal to another. This will also update the data structs

Parameters
----------
dest_lat : str
    Maybe: 'model', 'design', or 'base' with optional '@n' at beginning to indicate the universe

source_lat : str
    Maybe: 'model', 'design', or 'base'
)"""
  );
  nb::class_<Tao::TaoSetLogicalValue>(m, "TaoSetLogicalValue", "tao_set_logical_value return type")
      .def_ro("var", &Tao::TaoSetLogicalValue::var)
      .def_ro("error", &Tao::TaoSetLogicalValue::error)
      .def("__len__", [](const Tao::TaoSetLogicalValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetLogicalValue &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.var);
        if (i == 1)
          return nb::cast(s.error);
        throw nb::index_error();
      });
  m.def(
      "tao_set_logical_value",
      &Tao::tao_set_logical_value,
      nb::arg("var_str"),
      nb::arg("value_str"),
      R"""(Subroutine to read and set the value of an logical varialbe.

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
      nb::arg("n_threads"),
      R"""(Routine to set OpenMP thread count.  Errors if OpenMP is not available.

Parameters
----------
n_threads : int
    Number of threads.
)"""
  );
  m.def(
      "tao_set_opt_vars",
      &Tao::tao_set_opt_vars,
      nb::arg("var_vec"),
      nb::arg("print_limit_warning") = nb::none(),
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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set opti_de_param variables

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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set particle_start variables.

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
      nb::arg("plot_name"),
      nb::arg("component"),
      nb::arg("value_str"),
      R"""(Routine to set plot parameters.

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
      nb::arg("component"),
      nb::arg("value_str"),
      nb::arg("value_str2") = nb::none(),
      R"""(Set various aspects of the plotting window

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
      "tao_set_plotting",
      &Tao::tao_set_plotting,
      nb::arg("plot_input"),
      nb::arg("plot_page"),
      nb::arg("use_cmd_line_geom"),
      nb::arg("reverse") = nb::none(),
      R"""(Wrapper for Fortran routine tao_set_plotting

Parameters
----------
plot_input : TaoPlotPageInput

plot_page : TaoPlotPageStruct

use_cmd_line_geom : bool

reverse : bool, optional
)"""
  );
  m.def(
      "tao_set_ptc_com_cmd",
      &Tao::tao_set_ptc_com_cmd,
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set ptc_com variables

Parameters
----------
who : str
    which ptc_com variable to set

value_str : str
    Value to set to.
)"""
  );
  nb::class_<Tao::TaoSetQpAxisStruct>(m, "TaoSetQpAxisStruct", "tao_set_qp_axis_struct return type")
      .def_ro("error", &Tao::TaoSetQpAxisStruct::error)
      .def_ro("ix_uni", &Tao::TaoSetQpAxisStruct::ix_uni)
      .def("__len__", [](const Tao::TaoSetQpAxisStruct &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetQpAxisStruct &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.error);
        if (i == 1)
          return nb::cast(s.ix_uni);
        throw nb::index_error();
      });
  m.def(
      "tao_set_qp_axis_struct",
      &Tao::tao_set_qp_axis_struct,
      nb::arg("qp_axis_name"),
      nb::arg("component"),
      nb::arg("qp_axis"),
      nb::arg("value"),
      R"""(Routine to set qp_axis_names of a qp_axis_struct.

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
  nb::class_<Tao::TaoSetQpPointStruct>(
      m,
      "TaoSetQpPointStruct",
      "tao_set_qp_point_struct return type"
  )
      .def_ro("error", &Tao::TaoSetQpPointStruct::error)
      .def_ro("ix_uni", &Tao::TaoSetQpPointStruct::ix_uni)
      .def("__len__", [](const Tao::TaoSetQpPointStruct &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetQpPointStruct &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.error);
        if (i == 1)
          return nb::cast(s.ix_uni);
        throw nb::index_error();
      });
  m.def(
      "tao_set_qp_point_struct",
      &Tao::tao_set_qp_point_struct,
      nb::arg("qp_point_name"),
      nb::arg("component"),
      nb::arg("qp_point"),
      nb::arg("value"),
      R"""(Routine to set qp_point_names of a qp_point_struct.

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
  nb::class_<Tao::TaoSetQpRectStruct>(m, "TaoSetQpRectStruct", "tao_set_qp_rect_struct return type")
      .def_ro("error", &Tao::TaoSetQpRectStruct::error)
      .def_ro("ix_uni", &Tao::TaoSetQpRectStruct::ix_uni)
      .def("__len__", [](const Tao::TaoSetQpRectStruct &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetQpRectStruct &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.error);
        if (i == 1)
          return nb::cast(s.ix_uni);
        throw nb::index_error();
      });
  m.def(
      "tao_set_qp_rect_struct",
      &Tao::tao_set_qp_rect_struct,
      nb::arg("qp_rect_name"),
      nb::arg("component"),
      nb::arg("qp_rect"),
      nb::arg("value"),
      R"""(Routine to set qp_rect_names of a qp_rect_struct.

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
      nb::arg("state_string"),
      R"""(Sets the random number generator state.

Parameters
----------
state_string : str
    Encoded random number state.
)"""
  );
  nb::class_<Tao::TaoSetRealValue>(m, "TaoSetRealValue", "tao_set_real_value return type")
      .def_ro("var", &Tao::TaoSetRealValue::var)
      .def_ro("error", &Tao::TaoSetRealValue::error)
      .def("__len__", [](const Tao::TaoSetRealValue &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSetRealValue &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.var);
        if (i == 1)
          return nb::cast(s.error);
        throw nb::index_error();
      });
  m.def(
      "tao_set_real_value",
      &Tao::tao_set_real_value,
      nb::arg("var_str"),
      nb::arg("value_str"),
      nb::arg("min_val") = nb::none(),
      nb::arg("max_val") = nb::none(),
      nb::arg("dflt_uni") = nb::none(),
      R"""(Subroutine to read and set the value of a real variable.

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
      nb::arg("region_name"),
      nb::arg("component"),
      nb::arg("value_str"),
      R"""(Routine to set region parameters.

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
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set space_charge_com variables

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
      nb::arg("sym_str"),
      nb::arg("num_str") = nb::none(),
      nb::arg("val") = nb::none(),
      R"""(Associates a given symbol with a given number.
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
      nb::arg("branch_str"),
      nb::arg("mask_str"),
      nb::arg("print_list"),
      nb::arg("qa_str"),
      nb::arg("qb_str"),
      nb::arg("delta_input"),
      R"""(Routine to set the transverse tunes.

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
      nb::arg("uni"),
      nb::arg("who"),
      nb::arg("what"),
      R"""(Sets a universe on or off, or sets the recalculate or twiss_calc logicals, etc.

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
      nb::arg("var_str"),
      nb::arg("value_str"),
      R"""(Routine to set var values.

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
      nb::arg("var"),
      nb::arg("value"),
      nb::arg("print_limit_warning") = nb::none(),
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
      R"""(Wrapper for Fortran routine tao_set_var_useit_opt
)"""
  );
  m.def(
      "tao_set_wave_cmd",
      &Tao::tao_set_wave_cmd,
      nb::arg("who"),
      nb::arg("value_str"),
      R"""(Routine to set wave variables

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
      nb::arg("branch_str"),
      nb::arg("q_str"),
      nb::arg("delta_input"),
      R"""(Routine to set the z-tune.

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
      R"""(Wrapper for Fortran routine tao_setup_key_table
)"""
  );
  m.def(
      "tao_shape_init",
      &Tao::tao_shape_init,
      nb::arg("shape"),
      nb::arg("print_err") = nb::none(),
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
      nb::arg("what"),
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
      nb::arg("iunit"),
      nb::arg("form"),
      R"""(Routine to show a list of datums and variables and how they contribute to the merit function.

Parameters
----------
iunit : int
    File unit to write to. 0 => print to the terminal.

form : str
    What to output: 'ALL'   -> All datums and variables. 'TOP10' -> Top datums and variables that contribute
    to the merit function.
)"""
  );
  nb::class_<Tao::TaoShowThis>(m, "TaoShowThis", "tao_show_this return type")
      .def_ro("result_id", &Tao::TaoShowThis::result_id)
      .def_ro("lines", &Tao::TaoShowThis::lines)
      .def_ro("nl", &Tao::TaoShowThis::nl)
      .def("__len__", [](const Tao::TaoShowThis &) { return 3; })
      .def("__getitem__", [](const Tao::TaoShowThis &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.result_id);
        if (i == 1)
          return nb::cast(s.lines);
        if (i == 2)
          return nb::cast(s.nl);
        throw nb::index_error();
      });
  m.def(
      "tao_show_this",
      &Tao::tao_show_this,
      nb::arg("what"),
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
      nb::arg("char_"),
      R"""(Wrapper for Fortran routine tao_single_mode

Parameters
----------
)"""
  );
  m.def(
      "tao_single_track",
      &Tao::tao_single_track,
      nb::arg("tao_lat"),
      nb::arg("ix_branch"),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to track a single particle and calculate lattice functions through a lattice.

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
      nb::arg("data_type"),
      nb::arg("data_source"),
      R"""(Wrapper for Fortran routine tao_spin_matrices_calc_needed

Parameters
----------
data_type : str

data_source : str

Returns
-------
do_calc : bool
)"""
  );
  m.def(
      "tao_spin_tracking_turn_on",
      &Tao::tao_spin_tracking_turn_on,
      R"""(Wrapper for Fortran routine tao_spin_tracking_turn_on
)"""
  );
  nb::class_<Tao::TaoSplitComponent>(m, "TaoSplitComponent", "tao_split_component return type")
      .def_ro("comp", &Tao::TaoSplitComponent::comp)
      .def_ro("err", &Tao::TaoSplitComponent::err)
      .def("__len__", [](const Tao::TaoSplitComponent &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSplitComponent &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.comp);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "tao_split_component",
      &Tao::tao_split_component,
      nb::arg("comp_str"),
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
      nb::arg("data_type"),
      nb::arg("data_source"),
      R"""(Wrapper for Fortran routine tao_srdt_calc_needed

Parameters
----------
data_type : str

data_source : str

Returns
-------
do_srdt : int
)"""
  );
  nb::class_<Tao::TaoSubinUniNumber>(m, "TaoSubinUniNumber", "tao_subin_uni_number return type")
      .def_ro("name_out", &Tao::TaoSubinUniNumber::name_out)
      .def_ro("ok", &Tao::TaoSubinUniNumber::ok)
      .def("__len__", [](const Tao::TaoSubinUniNumber &) { return 2; })
      .def("__getitem__", [](const Tao::TaoSubinUniNumber &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.name_out);
        if (i == 1)
          return nb::cast(s.ok);
        throw nb::index_error();
      });
  m.def(
      "tao_subin_uni_number",
      &Tao::tao_subin_uni_number,
      nb::arg("name_in"),
      nb::arg("ix_uni"),
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
      R"""(Routine to minimize the merit function using svd.

Returns
-------
abort : bool
    Set True if svd step increases the merit function.
)"""
  );
  m.def(
      "tao_symbol_import_from_lat",
      &Tao::tao_symbol_import_from_lat,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine tao_symbol_import_from_lat

Parameters
----------
lat : LatStruct
)"""
  );
  m.def(
      "tao_taper_cmd",
      &Tao::tao_taper_cmd,
      nb::arg("except_"),
      nb::arg("uni_names"),
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
      nb::arg("num_str"),
      nb::arg("n_size"),
      nb::arg("change_number"),
      nb::arg("abs_or_rel"),
      nb::arg("err"),
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
      nb::arg("str"),
      nb::arg("i_int"),
      nb::arg("err"),
      R"""(Converts a string to an integer

If the string str is blank then i_int = 0
)"""
  );
  nb::class_<Tao::TaoToPhaseAndCouplingReading>(
      m,
      "TaoToPhaseAndCouplingReading",
      "tao_to_phase_and_coupling_reading return type"
  )
      .def_ro("bpm_data", &Tao::TaoToPhaseAndCouplingReading::bpm_data)
      .def_ro("valid_value", &Tao::TaoToPhaseAndCouplingReading::valid_value)
      .def("__len__", [](const Tao::TaoToPhaseAndCouplingReading &) { return 2; })
      .def("__getitem__", [](const Tao::TaoToPhaseAndCouplingReading &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.bpm_data);
        if (i == 1)
          return nb::cast(s.valid_value);
        throw nb::index_error();
      });
  m.def(
      "tao_to_phase_and_coupling_reading",
      &Tao::tao_to_phase_and_coupling_reading,
      nb::arg("ele"),
      nb::arg("why_invalid"),
      nb::arg("datum"),
      R"""(Buffer routine for to_phase_and_coupling_reading.

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
  nb::class_<Tao::TaoToReal>(m, "TaoToReal", "tao_to_real return type")
      .def_ro("value", &Tao::TaoToReal::value)
      .def_ro("err_flag", &Tao::TaoToReal::err_flag)
      .def("__len__", [](const Tao::TaoToReal &) { return 2; })
      .def("__getitem__", [](const Tao::TaoToReal &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.value);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "tao_to_real",
      &Tao::tao_to_real,
      nb::arg("expression"),
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
      "tao_to_top10",
      &Tao::tao_to_top10,
      nb::arg("top10"),
      nb::arg("value"),
      nb::arg("name"),
      nb::arg("c_index"),
      nb::arg("order"),
      R"""(Routine to order the largest contributors to the merit function in
a list. Call this routine for each contributor.

Note: Before first calling this routine set:
  top10(:)%valid = .false.

Parameters
----------
top10 : 1D array of TaoTop10Struct
    List of top contributors. Note that the list is not limited to 10 entries.

value : float
    value of the contributor.

name : str
    Name of the contributor..

c_index : int
    Index of the contributor.

order : str
    Ordering of the list. Possibilities are:
)"""
  );
  m.def(
      "tao_too_many_particles_lost",
      &Tao::tao_too_many_particles_lost,
      nb::arg("beam"),
      R"""(Wrapper for Fortran routine tao_too_many_particles_lost

Parameters
----------
beam : BeamStruct

Returns
-------
no_beam : bool
)"""
  );
  m.def(
      "tao_top10_derivative_print",
      &Tao::tao_top10_derivative_print,
      R"""(Routine to print out the top10 contributors to the merit function.
)"""
  );
  m.def(
      "tao_top10_merit_categories_print",
      &Tao::tao_top10_merit_categories_print,
      nb::arg("iunit"),
      R"""(Routine to print the top data and variable categories that contribute to
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
      nb::arg("command") = nb::none(),
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
  nb::class_<Tao::TaoTrackingEleIndex>(
      m,
      "TaoTrackingEleIndex",
      "tao_tracking_ele_index return type"
  )
      .def_ro("ix_branch", &Tao::TaoTrackingEleIndex::ix_branch)
      .def_ro("ix_ele", &Tao::TaoTrackingEleIndex::ix_ele)
      .def("__len__", [](const Tao::TaoTrackingEleIndex &) { return 2; })
      .def("__getitem__", [](const Tao::TaoTrackingEleIndex &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ix_branch);
        if (i == 1)
          return nb::cast(s.ix_ele);
        throw nb::index_error();
      });
  m.def(
      "tao_tracking_ele_index",
      &Tao::tao_tracking_ele_index,
      nb::arg("ele"),
      nb::arg("datum"),
      R"""(Routine to return the index in the tracking part of a lattice that corresponds to ele.

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
      R"""(Wrapper for Fortran routine tao_turn_on_special_calcs_if_needed_for_plotting
)"""
  );
  m.def(
      "tao_type_expression_tree",
      &Tao::tao_type_expression_tree,
      nb::arg("tree"),
      nb::arg("indent") = nb::none(),
      R"""(Routine to print an expression tree in tree form.
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
      nb::arg("string"),
      R"""(Routine to return the index of an atsign ("@") character in a string if the atsign is
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
      nb::arg("i_uni"),
      nb::arg("neg2_to_default") = nb::none(),
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
      nb::arg("action"),
      nb::arg("data_name"),
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
      nb::arg("action"),
      nb::arg("var_name"),
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
      nb::arg("var"),
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
      nb::arg("var"),
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
      nb::arg("eles"),
      nb::arg("attribute"),
      nb::arg("silent"),
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
      R"""(Wrapper for Fortran routine tao_var_repoint
)"""
  );
  m.def(
      "tao_var_show_use",
      [](TaoV1VarStruct &v1_var, CharacterAlloc1D *lines, std::optional<int> nl) {
        auto fn = static_cast<
            void (*)(TaoV1VarStruct &, optional_ref<CharacterAlloc1D>, std::optional<int>)>(
            &Tao::tao_var_show_use
        );
        return fn(v1_var, ptr_to_opt_ref(lines), nl);
      },
      nb::arg("v1_var"),
      nb::arg("lines") = nb::none(),
      nb::arg("nl") = nb::none(),
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
      R"""(Wrapper for Fortran routine tao_var_target_calc
)"""
  );
  m.def(
      "tao_var_useit_plot_calc",
      &Tao::tao_var_useit_plot_calc,
      nb::arg("graph"),
      nb::arg("var"),
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
      nb::arg("out_file"),
      nb::arg("show_good_opt_only") = nb::none(),
      nb::arg("tao_format") = nb::none(),
      R"""(Routine to write the optimized variables. One file will be created for each universe.
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
      R"""(Routine to veto all variables with zero effect on data used in the merit function.
)"""
  );
  m.def(
      "tao_wave_analysis",
      &Tao::tao_wave_analysis,
      nb::arg("plot"),
      R"""(Routine to do a wave anaysis.

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
      nb::arg("curve_name"),
      nb::arg("plot_place"),
      nb::arg("err_flag"),
      R"""(Routine to do the initial setup for wave plotting.
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
      nb::arg("curve"),
      nb::arg("ix1"),
      nb::arg("n_dat"),
      nb::arg("coef"),
      nb::arg("rms"),
      nb::arg("f1"),
      nb::arg("f2") = nb::none(),
      nb::arg("f3") = nb::none(),
      nb::arg("f4") = nb::none(),
      R"""(Routine for fitting the curve data to up to four functions using a least squares
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
      nb::arg("what"),
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
      nb::arg("iunit"),
      nb::arg("line"),
      R"""(Subroutine to write out a series of lines to a file or to the terminal.
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
      nb::arg("where"),
      nb::arg("what"),
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
      nb::arg("where"),
      nb::arg("x_min_in"),
      nb::arg("x_max_in"),
      nb::arg("include_wall") = nb::none(),
      nb::arg("gang") = nb::none(),
      nb::arg("exact") = nb::none(),
      nb::arg("turn_autoscale_off") = nb::none(),
      R"""(Routine to scale a plot. If x_min = x_max
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
      nb::arg("graph"),
      nb::arg("x_min"),
      nb::arg("x_max"),
      nb::arg("include_wall") = nb::none(),
      nb::arg("have_scaled") = nb::none(),
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
      nb::arg("plot"),
      nb::arg("x_min_in"),
      nb::arg("x_max_in"),
      nb::arg("include_wall") = nb::none(),
      nb::arg("gang") = nb::none(),
      R"""(Routine to scale a plot. If x_min = x_max
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
