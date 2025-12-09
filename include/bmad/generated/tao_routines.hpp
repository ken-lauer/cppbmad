#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.h"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"

using namespace Bmad;

namespace Tao {

// Skipped unusable routine avv:
// Routine in configuration skip list

// Skipped unusable routine callback:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine complex_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine do_loop_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_integrate_max(
    c_Int& ix_start /* 0D_NOT_integer */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Real& datum_value /* 0D_NOT_real */,
    c_Int& ix_m /* 0D_NOT_integer */,
    void* branch /* 0D_NOT_type */,
    void* vec /* 1D_ALLOC_real */,
    void* datum /* 0D_NOT_type */);
void integrate_max(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum);
extern "C" void fortran_integrate_min(
    c_Int& ix_start /* 0D_NOT_integer */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    c_Real& datum_value /* 0D_NOT_real */,
    c_Int& ix_m /* 0D_NOT_integer */,
    void* branch /* 0D_NOT_type */,
    void* vec /* 1D_ALLOC_real */,
    void* datum /* 0D_NOT_type */);
void integrate_min(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum);

// Skipped unusable routine jacobian:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_re_allocate_c_double(
    void* re /* 1D_ALLOC_real */,
    c_Int& n /* 0D_NOT_integer */,
    c_Bool* exact /* 0D_NOT_logical */,
    c_Real* init_val /* 0D_NOT_real */);
void re_allocate_c_double(
    RealAlloc1D& re,
    int n,
    std::optional<bool> exact = std::nullopt,
    optional_ref<double> init_val = std::nullopt);
extern "C" void fortran_tao_abort_command_file(
    c_Bool* force_abort /* 0D_NOT_logical */);
void tao_abort_command_file(std::optional<bool> force_abort = std::nullopt);
extern "C" void fortran_tao_add_to_normal_mode_h_array(
    c_Char h_str /* 0D_NOT_character */,
    void* h_array /* 1D_ALLOC_type */);
ResonanceHProxyAlloc1D tao_add_to_normal_mode_h_array(std::string h_str);
extern "C" void fortran_tao_alias_cmd(
    c_Char alias /* 0D_NOT_character */,
    c_Char string /* 0D_NOT_character */);
void tao_alias_cmd(std::string alias, std::string string);

// Skipped unusable routine tao_alias_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_allocate_data_array(
    void* u /* 0D_NOT_type */,
    c_Int& n_data /* 0D_NOT_integer */,
    c_Bool* exact /* 0D_NOT_logical */);
void tao_allocate_data_array(
    TaoUniverseProxy& u,
    int n_data,
    optional_ref<bool> exact = std::nullopt);
extern "C" void fortran_tao_allocate_v1_var(
    c_Int& n_v1 /* 0D_NOT_integer */,
    c_Bool& save_old /* 0D_NOT_logical */);
void tao_allocate_v1_var(int n_v1, bool save_old);
extern "C" void fortran_tao_allocate_var_array(
    c_Int& n_var /* 0D_NOT_integer */,
    c_Bool& default_good_user /* 0D_NOT_logical */);
void tao_allocate_var_array(int n_var, bool default_good_user);

// Skipped unusable routine tao_beam_branch_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_beam_emit_calc(
    c_Int& plane /* 0D_NOT_integer */,
    c_Int& emit_type /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    void* bunch_params /* 0D_NOT_type */,
    c_Real& emit /* 0D_NOT_real */);
void tao_beam_emit_calc(
    int plane,
    int emit_type,
    EleProxy& ele,
    BunchParamsProxy& bunch_params,
    double emit);

// Skipped unusable routine tao_beam_shake_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_beam_track(
    void* u /* 0D_NOT_type */,
    void* tao_lat /* 0D_NOT_type */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    void* beam /* 0D_NOT_type */,
    c_Bool& calc_ok /* 0D_NOT_logical */);
bool tao_beam_track(
    TaoUniverseProxy& u,
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    BeamProxy& beam);
extern "C" bool fortran_tao_beam_track_endpoint(
    c_Char ele_id /* 0D_NOT_character */,
    void* lat /* 0D_NOT_type */,
    c_Char branch_str /* 0D_NOT_character */,
    c_Char where /* 0D_NOT_character */,
    void* u /* 0D_NOT_type */,
    void* ele /* 0D_PTR_type */);
void tao_beam_track_endpoint(
    std::string ele_id,
    LatProxy& lat,
    std::string branch_str,
    std::string where,
    TaoUniverseProxy& u,
    EleProxy& ele);

// Skipped unusable routine tao_beam_uni_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_branch_index(
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Int& ix_this /* 0D_NOT_integer */);
void tao_branch_index(int ix_branch, int ix_this);

// Skipped unusable routine tao_building_wall_orientation_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_building_wall_point_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_building_wall_section_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_building_wall_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_c_command:
// Argument not defined: c_str (have: [])
// Argument not defined: tao_c_command (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_beam_track_element:
// Argument not defined: tao_c_get_beam_track_element (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_integer_array:
// Argument not defined: tao_c_get_integer_array (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_real_array:
// Argument not defined: tao_c_get_real_array (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_string_buffer:
// No matching docstring

// Skipped unusable routine tao_c_get_string_buffer_length:
// Argument not defined: tao_c_get_string_buffer_length (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_init_tao:
// Argument not defined: c_str (have: [])
// Argument not defined: tao_c_init_tao (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_integer_array_size:
// Argument not defined: tao_c_integer_array_size (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_interface_common_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_c_out_io_buffer_get_line:
// Argument not defined: n (have: [])
// Argument not defined: tao_c_out_io_buffer_get_line (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_out_io_buffer_num_lines:
// Argument not defined: tao_c_out_io_buffer_num_lines (have: [])
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_c_out_io_buffer_reset();
void tao_c_out_io_buffer_reset();

// Skipped unusable routine tao_c_real_array_size:
// Argument not defined: tao_c_real_array_size (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_string_size:
// Argument not defined: tao_c_string_size (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_calc_data_at_s_pts:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_call_cmd:
// Variable-sized in character array: cmd_arg(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_cbar_wave_anal:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
extern "C" void fortran_tao_change_ele(
    c_Char ele_name /* 0D_NOT_character */,
    c_Char attrib_name /* 0D_NOT_character */,
    c_Char num_str /* 0D_NOT_character */,
    c_Bool& update /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_change_ele(
    std::string ele_name,
    std::string attrib_name,
    std::string num_str,
    bool update);
extern "C" void fortran_tao_change_tune(
    c_Char branch_str /* 0D_NOT_character */,
    c_Char mask_str /* 0D_NOT_character */,
    c_Bool& print_list /* 0D_NOT_logical */,
    c_Char dqa_str /* 0D_NOT_character */,
    c_Char dqb_str /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_change_tune(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string dqa_str,
    std::string dqb_str);
extern "C" void fortran_tao_change_var(
    c_Char name /* 0D_NOT_character */,
    c_Char num_str /* 0D_NOT_character */,
    c_Bool& silent /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_change_var(std::string name, std::string num_str, bool silent);
extern "C" void fortran_tao_change_z_tune(
    c_Char branch_str /* 0D_NOT_character */,
    c_Char dq_str /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_change_z_tune(std::string branch_str, std::string dq_str);
extern "C" bool fortran_tao_chrom_calc_needed(
    c_Char data_type /* 0D_NOT_character */,
    c_Char data_source /* 0D_NOT_character */,
    c_Bool& do_chrom /* 0D_NOT_logical */);
void tao_chrom_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_chrom);
extern "C" void fortran_tao_clear_cmd(c_Char cmd_line /* 0D_NOT_character */);
void tao_clear_cmd(std::string cmd_line);
extern "C" void fortran_tao_clip_cmd(
    c_Bool& gang /* 0D_NOT_logical */,
    c_Char where /* 0D_NOT_character */,
    c_Real& value1 /* 0D_NOT_real */,
    c_Real& value2 /* 0D_NOT_real */);
void tao_clip_cmd(bool gang, std::string where, double value1, double value2);
extern "C" void fortran_tao_close_command_file();
void tao_close_command_file();

// Skipped unusable routine tao_cmd_end_calc:
// Module name unset? Internal error
extern "C" void fortran_tao_cmd_history_record(
    c_Char cmd /* 0D_NOT_character */);
void tao_cmd_history_record(std::string cmd);

// Skipped unusable routine tao_cmd_history_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_cmd_split:
// Variable-sized out character array: cmd_word(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_command(
    c_Char command_line /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool& err_is_fatal /* 0D_NOT_logical */);
bool tao_command(std::string command_line, bool err);

// Skipped unusable routine tao_command_file_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_common_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_constraint_type_name(
    void* datum /* 0D_NOT_type */,
    c_Char datum_name /* 0D_NOT_character */);
void tao_constraint_type_name(TaoDataProxy& datum, std::string datum_name);

// Skipped unusable routine tao_control_tree_list:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
extern "C" void fortran_tao_count_strings(
    c_Char string /* 0D_NOT_character */,
    c_Char pattern /* 0D_NOT_character */,
    c_Int& num /* 0D_NOT_integer */);
int tao_count_strings(std::string string, std::string pattern);
extern "C" void fortran_tao_create_plot_window();
void tao_create_plot_window();

// Skipped unusable routine tao_curve_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_curve_beam_ellipse_setup:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_check_universe:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_color_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_curve_data_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_datum_calc:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_ele_ref:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_curve_ix_uni:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_name:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_orbit_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_curve_rms_calc:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_curve_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_d1_data_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_d1_data_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_d1_data_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_d2_d1_name(
    void* d1 /* 0D_NOT_type */,
    c_Bool* show_universe /* 0D_NOT_logical */,
    c_Char d2_d1_name /* 0D_NOT_character */);
void tao_d2_d1_name(
    TaoD1DataProxy& d1,
    std::optional<bool> show_universe,
    std::string d2_d1_name);

// Skipped unusable routine tao_d2_data_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_d2_data_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_d2_data_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_d2_data_stuffit(
    void* u /* 0D_NOT_type */,
    c_Char d2_name /* 0D_NOT_character */,
    c_Int& n_d1_data /* 0D_NOT_integer */);
void tao_d2_data_stuffit(
    TaoUniverseProxy& u,
    std::string d2_name,
    int n_d1_data);

// Skipped unusable routine tao_data_array_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_data_check(c_Bool& err /* 0D_NOT_logical */);
void tao_data_check(bool err);
extern "C" void fortran_tao_data_coupling_init(void* branch /* 0D_NOT_type */);
void tao_data_coupling_init(BranchProxy& branch);
extern "C" bool fortran_tao_data_sanity_check(
    void* datum /* 0D_NOT_type */,
    c_Bool& print_err /* 0D_NOT_logical */,
    c_Char default_data_type /* 0D_NOT_character */,
    void* uni /* 0D_NOT_type */,
    c_Bool& is_valid /* 0D_NOT_logical */);
void tao_data_sanity_check(
    TaoDataProxy& datum,
    bool print_err,
    std::string default_data_type,
    optional_ref<TaoUniverseProxy> uni,
    bool is_valid);

// Skipped unusable routine tao_data_show_use:
// Variable-sized inout character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_data_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_data_type_substitute:
// Untranslated type: TaoCurveProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_data_useit_plot_calc:
// Untranslated type: TaoCurveProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_data_var_component_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_datum_has_associated_ele(
    c_Char data_type /* 0D_NOT_character */,
    c_Int* branch_geometry /* 0D_NOT_integer */,
    c_Int& has_associated_ele /* 0D_NOT_integer */);
void tao_datum_has_associated_ele(
    std::string data_type,
    std::optional<int> branch_geometry,
    int has_associated_ele);

// Skipped unusable routine tao_datum_input_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_datum_integrate(
    void* datum /* 0D_NOT_type */,
    void* branch /* 0D_NOT_type */,
    void* s_pos /* 1D_ALLOC_real */,
    void* values /* 1D_ALLOC_real */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Real& result /* 0D_NOT_real */);
struct TaoDatumIntegrate {
  bool valid_value;
  std::string why_invalid;
  double result;
};
TaoDatumIntegrate tao_datum_integrate(
    TaoDataProxy& datum,
    BranchProxy& branch,
    RealAlloc1D& s_pos,
    RealAlloc1D& values);
extern "C" bool fortran_tao_datum_name(
    void* datum /* 0D_NOT_type */,
    c_Bool* show_universe /* 0D_NOT_logical */,
    c_Char datum_name /* 0D_NOT_character */);
void tao_datum_name(
    TaoDataProxy& datum,
    std::optional<bool> show_universe,
    std::string datum_name);
extern "C" bool fortran_tao_datum_s_position(
    void* datum /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Real& s_pos /* 0D_NOT_real */);
double tao_datum_s_position(TaoDataProxy& datum, EleProxy& ele);
extern "C" void fortran_tao_de_optimizer(c_Bool& abort /* 0D_NOT_logical */);
bool tao_de_optimizer();
extern "C" void fortran_tao_deallocate_plot_cache(
    void* plot_cache /* 1D_ALLOC_type */);
void tao_deallocate_plot_cache(TaoPlotCacheProxyAlloc1D& plot_cache);

// Skipped unusable routine tao_deallocate_tree:
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)

// Skipped unusable routine tao_design_lat_input_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_destroy_plot_window();
void tao_destroy_plot_window();
extern "C" void fortran_tao_dmerit_calc();
void tao_dmerit_calc();
extern "C" void fortran_tao_dmodel_dvar_calc(
    c_Bool& force_calc /* 0D_NOT_logical */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_dmodel_dvar_calc(bool force_calc);
extern "C" bool fortran_tao_do_wire_scan(
    void* ele /* 0D_NOT_type */,
    c_Real& theta /* 0D_NOT_real */,
    void* beam /* 0D_NOT_type */,
    c_Real& moment /* 0D_NOT_real */);
double tao_do_wire_scan(EleProxy& ele, double theta, BeamProxy& beam);

// Skipped unusable routine tao_draw_beam_chamber_wall:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_draw_curve_data:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_draw_ele_for_floor_plan:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)
// Untranslated type: TaoEleShapeProxy (0D_PTR_type)

// Skipped unusable routine tao_draw_floor_plan:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_draw_graph_axes:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_draw_histogram_data:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_draw_lat_layout:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)
extern "C" void fortran_tao_draw_plots(c_Bool* do_clear /* 0D_NOT_logical */);
void tao_draw_plots(std::optional<bool> do_clear = std::nullopt);

// Skipped unusable routine tao_drawing_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_dynamic_aperture_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_ele_geometry_with_misalignments(
    void* datum /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Real& value /* 0D_NOT_real */);
struct TaoEleGeometryWithMisalignments {
  bool valid_value;
  std::string why_invalid;
  double value;
};
TaoEleGeometryWithMisalignments tao_ele_geometry_with_misalignments(
    TaoDataProxy& datum,
    EleProxy& ele);

// Skipped unusable routine tao_ele_pointer_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_ele_shape_info:
// Untranslated type: TaoEleShapeProxy (1D_ALLOC_type)
// Untranslated type: TaoEleShapeProxy (0D_PTR_type)

// Skipped unusable routine tao_ele_shape_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_ele_shape_input_to_struct:
// Untranslated type: TaoEleShapeInputProxy (0D_NOT_type)
// Untranslated type: TaoEleShapeProxy (0D_NOT_type)

// Skipped unusable routine tao_ele_shape_struct_to_input:
// Untranslated type: TaoEleShapeProxy (0D_NOT_type)
// Untranslated type: TaoEleShapeInputProxy (0D_NOT_type)

// Skipped unusable routine tao_ele_shape_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_eval_floor_orbit(
    void* datum /* 0D_NOT_type */,
    void* ele /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    void* bunch_params /* 0D_NOT_type */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Real& value /* 0D_NOT_real */);
struct TaoEvalFloorOrbit {
  bool valid_value;
  std::string why_invalid;
  double value;
};
TaoEvalFloorOrbit tao_eval_floor_orbit(
    TaoDataProxy& datum,
    EleProxy& ele,
    CoordProxy& orbit,
    BunchParamsProxy& bunch_params);

// Skipped unusable routine tao_eval_node_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_evaluate_a_datum(
    void* datum /* 0D_NOT_type */,
    void* u /* 0D_NOT_type */,
    void* tao_lat /* 0D_NOT_type */,
    c_Real& datum_value /* 0D_NOT_real */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Bool* called_from_lat_calc /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */);
struct TaoEvaluateADatum {
  double datum_value;
  bool valid_value;
  std::string why_invalid;
};
TaoEvaluateADatum tao_evaluate_a_datum(
    TaoDataProxy& datum,
    TaoUniverseProxy& u,
    TaoLatticeProxy& tao_lat,
    std::optional<bool> called_from_lat_calc = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" bool fortran_tao_evaluate_datum_at_s(
    void* datum /* 0D_NOT_type */,
    void* tao_lat /* 0D_NOT_type */,
    void* ele /* 0D_PTR_type */,
    void* ele_ref /* 0D_PTR_type */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    c_Char err_str /* 0D_NOT_character */,
    c_Bool& bad_datum /* 0D_NOT_logical */,
    c_Real& value /* 0D_NOT_real */);
struct TaoEvaluateDatumAtS {
  std::string err_str;
  bool bad_datum;
  double value;
};
TaoEvaluateDatumAtS tao_evaluate_datum_at_s(
    TaoDataProxy& datum,
    TaoLatticeProxy& tao_lat,
    EleProxy& ele,
    EleProxy& ele_ref,
    bool valid_value);

// Skipped unusable routine tao_evaluate_element_parameters:
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)

// Skipped unusable routine tao_evaluate_expression:
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)
// Untranslated type: TaoEvalNodeProxy (1D_ALLOC_type)

// Skipped unusable routine tao_evaluate_expression_new:
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)
// Untranslated type: TaoEvalNodeProxy (1D_ALLOC_type)

// Skipped unusable routine tao_evaluate_expression_old:
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)
// Untranslated type: TaoEvalNodeProxy (1D_ALLOC_type)
extern "C" void fortran_tao_evaluate_lat_or_beam_data(
    c_Bool& err /* 0D_NOT_logical */,
    c_Char data_name /* 0D_NOT_character */,
    void* values /* 1D_ALLOC_real */,
    c_Bool& print_err /* 0D_NOT_logical */,
    c_Char default_source /* 0D_NOT_character */,
    void* dflt_ele_ref /* 0D_PTR_type */,
    void* dflt_ele_start /* 0D_PTR_type */,
    void* dflt_ele /* 0D_PTR_type */,
    c_Char dflt_component /* 0D_NOT_character */,
    c_Int* dflt_uni /* 0D_NOT_integer */,
    c_Int* dflt_eval_point /* 0D_NOT_integer */,
    c_Real* dflt_s_offset /* 0D_NOT_real */);
struct TaoEvaluateLatOrBeamData {
  bool err;
  RealAlloc1D values;
};
TaoEvaluateLatOrBeamData tao_evaluate_lat_or_beam_data(
    std::string data_name,
    bool print_err,
    std::string default_source,
    optional_ref<EleProxy> dflt_ele_ref = std::nullopt,
    optional_ref<EleProxy> dflt_ele_start = std::nullopt,
    optional_ref<EleProxy> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt);

// Skipped unusable routine tao_evaluate_stack_old:
// Untranslated type: TaoEvalNodeProxy (1D_ALLOC_type)
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)

// Skipped unusable routine tao_evaluate_tree:
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)
extern "C" bool fortran_tao_evaluate_tune(
    c_Char q_str /* 0D_NOT_character */,
    c_Real& q0 /* 0D_NOT_real */,
    c_Bool& delta_input /* 0D_NOT_logical */,
    c_Real& q_val /* 0D_NOT_real */);
void tao_evaluate_tune(
    std::string q_str,
    double q0,
    bool delta_input,
    double q_val);
extern "C" void fortran_tao_expression_hash_substitute(
    c_Char expression_in /* 0D_NOT_character */,
    c_Char expression_out /* 0D_NOT_character */,
    void* eval_ele /* 0D_PTR_type */);
std::string tao_expression_hash_substitute(
    std::string expression_in,
    optional_ref<EleProxy> eval_ele = std::nullopt);

// Skipped unusable routine tao_expression_info_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_expression_tree_to_string:
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)

// Skipped unusable routine tao_find_data:
// Untranslated type: TaoD2DataArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoD1DataArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoDataArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoRealPointerProxy (1D_ALLOC_type)
// Untranslated type: TaoLogicalArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoStringArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoIntegerArrayProxy (1D_ALLOC_type)

// Skipped unusable routine tao_find_plot_region:
// Untranslated type: TaoPlotRegionProxy (0D_PTR_type)

// Skipped unusable routine tao_find_plots:
// Untranslated type: TaoPlotArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoGraphArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoCurveArrayProxy (1D_ALLOC_type)

// Skipped unusable routine tao_find_var:
// Untranslated type: TaoV1VarArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoVarArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoRealPointerProxy (1D_ALLOC_type)
// Untranslated type: TaoLogicalArrayProxy (1D_ALLOC_type)
// Untranslated type: TaoStringArrayProxy (1D_ALLOC_type)
extern "C" void fortran_tao_fixer(
    c_Char switch_ /* 0D_NOT_character */,
    c_Char word1 /* 0D_NOT_character */,
    c_Char word2 /* 0D_NOT_character */);
void tao_fixer(std::string switch_, std::string word1, std::string word2);

// Skipped unusable routine tao_floor_plan_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_floor_to_screen:
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_floor_to_screen_coords:
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_geo_lm_func:
// Routine in configuration skip list
extern "C" void fortran_tao_geodesic_lm_optimizer(
    c_Bool& abort /* 0D_NOT_logical */);
bool tao_geodesic_lm_optimizer();
extern "C" void fortran_tao_get_data(
    void* data_value /* 1D_ALLOC_real */,
    void* data_weight /* 1D_ALLOC_real */,
    void* data_meas_value /* 1D_ALLOC_real */,
    void* data_ix_dModel /* 1D_ALLOC_integer */);
struct TaoGetData {
  RealAlloc1D data_value;
  RealAlloc1D data_weight;
  RealAlloc1D data_meas_value;
  IntAlloc1D data_ix_dModel;
};
TaoGetData tao_get_data();
extern "C" void fortran_tao_get_opt_vars(
    void* var_value /* 1D_ALLOC_real */,
    void* var_step /* 1D_ALLOC_real */,
    void* var_delta /* 1D_ALLOC_real */,
    void* var_weight /* 1D_ALLOC_real */,
    void* var_ix /* 1D_ALLOC_integer */,
    c_Bool& ignore_if_weight_is_zero /* 0D_NOT_logical */,
    c_Bool& ignore_if_not_limited /* 0D_NOT_logical */);
struct TaoGetOptVars {
  RealAlloc1D var_value;
  RealAlloc1D var_step;
  RealAlloc1D var_delta;
  RealAlloc1D var_weight;
  IntAlloc1D var_ix;
  bool ignore_if_weight_is_zero;
  bool ignore_if_not_limited;
};
TaoGetOptVars tao_get_opt_vars();
extern "C" void fortran_tao_get_user_input(
    c_Char cmd_out /* 0D_NOT_character */,
    c_Char prompt_str /* 0D_NOT_character */,
    c_Bool* wait_flag /* 0D_NOT_logical */,
    c_Char cmd_in /* 0D_NOT_character */);
std::string tao_get_user_input(
    std::optional<std::string> prompt_str = std::nullopt,
    std::optional<bool> wait_flag = std::nullopt,
    std::optional<std::string> cmd_in = std::nullopt);

// Skipped unusable routine tao_global_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_graph_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_graph_controller_setup:
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_data_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_data_slice_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_dynamic_aperture_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_histogram_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_graph_name:
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_phase_space_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_s_min_max_calc:
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_setup:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_graph_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_help:
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_histogram_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_hook_branch_calc_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_command_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_curve_s_pt_def:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_data_sanity_check_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_draw_floor_plan_def:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_draw_graph_def:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_evaluate_a_datum_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_graph_postsetup_def:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_graph_setup_def:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_init1_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init2_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init_beam_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init_data_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init_global_def:
// Untranslated type: TaoGlobalProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_init_lattice_post_parse_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init_plotting_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init_read_lattice_info_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_init_var_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_lattice_calc_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_merit_data_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_merit_var_def:
// Untranslated type: TaoVarProxy (0D_NOT_type)

// Skipped unusable routine tao_hook_optimizer_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_parse_command_args_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_plot_setup_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_post_process_data_def:
// Routine in configuration skip list

// Skipped unusable routine tao_hook_show_cmd_def:
// Variable-sized inout character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_init(c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_init();
extern "C" void fortran_tao_init_beam_in_universe(
    void* u /* 0D_NOT_type */,
    void* beam_init /* 0D_NOT_type */,
    c_Char track_start /* 0D_NOT_character */,
    c_Char track_end /* 0D_NOT_character */,
    c_Real& comb_ds_save /* 0D_NOT_real */);
void tao_init_beam_in_universe(
    TaoUniverseProxy& u,
    BeamInitProxy& beam_init,
    std::string track_start,
    std::string track_end,
    double comb_ds_save);
extern "C" void fortran_tao_init_beams(c_Char init_file /* 0D_NOT_character */);
void tao_init_beams(std::string init_file);

// Skipped unusable routine tao_init_building_wall:
// Module name unset? Internal error
extern "C" void fortran_tao_init_data(c_Char data_file /* 0D_NOT_character */);
void tao_init_data(std::string data_file);
extern "C" void fortran_tao_init_data_end_stuff();
void tao_init_data_end_stuff();
extern "C" void fortran_tao_init_data_in_universe(
    void* u /* 0D_NOT_type */,
    c_Int& n_d2_add /* 0D_NOT_integer */,
    c_Bool* keep_existing_data /* 0D_NOT_logical */);
void tao_init_data_in_universe(
    TaoUniverseProxy& u,
    int n_d2_add,
    optional_ref<bool> keep_existing_data = std::nullopt);
extern "C" void fortran_tao_init_dynamic_aperture(
    c_Char init_file /* 0D_NOT_character */);
void tao_init_dynamic_aperture(std::string init_file);

// Skipped unusable routine tao_init_find_elements:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)
extern "C" void fortran_tao_init_global(
    c_Char init_file /* 0D_NOT_character */);
void tao_init_global(std::string init_file);
extern "C" void fortran_tao_init_lattice(
    c_Char lat_file /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void tao_init_lattice(std::string lat_file, bool err_flag);
extern "C" void fortran_tao_init_plotting(
    c_Char plot_file /* 0D_NOT_character */);
void tao_init_plotting(std::string plot_file);

// Skipped unusable routine tao_init_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_init_variables(
    c_Char var_file /* 0D_NOT_character */);
void tao_init_variables(std::string var_file);
extern "C" void fortran_tao_inject_beam(
    void* u /* 0D_NOT_type */,
    void* model /* 0D_NOT_type */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    void* beam /* 0D_NOT_type */,
    c_Bool& init_ok /* 0D_NOT_logical */);
struct TaoInjectBeam {
  BeamProxy beam;
  bool init_ok;
};
TaoInjectBeam tao_inject_beam(
    TaoUniverseProxy& u,
    TaoLatticeProxy& model,
    int ix_branch);
extern "C" void fortran_tao_inject_particle(
    void* u /* 0D_NOT_type */,
    void* model /* 0D_NOT_type */,
    c_Int& ix_branch /* 0D_NOT_integer */);
void tao_inject_particle(
    TaoUniverseProxy& u,
    TaoLatticeProxy& model,
    int ix_branch);

// Skipped unusable routine tao_integer_array_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_is_valid_name(
    c_Char name /* 0D_NOT_character */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Bool& is_valid /* 0D_NOT_logical */);
std::string tao_is_valid_name(std::string name, bool is_valid);
extern "C" void fortran_tao_json_cmd(c_Char input_str /* 0D_NOT_character */);
void tao_json_cmd(std::string input_str);
extern "C" void fortran_tao_key_info_to_str(
    c_Int& ix_key /* 0D_NOT_integer */,
    c_Int& ix_min_key /* 0D_NOT_integer */,
    c_Int& ix_max_key /* 0D_NOT_integer */,
    c_Char key_str /* 0D_NOT_character */,
    c_Char header_str /* 0D_NOT_character */);
void tao_key_info_to_str(
    int ix_key,
    int ix_min_key,
    int ix_max_key,
    std::string key_str,
    std::string header_str);

// Skipped unusable routine tao_key_input_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_lat_bookkeeper(
    void* u /* 0D_NOT_type */,
    void* tao_lat /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */);
bool tao_lat_bookkeeper(TaoUniverseProxy& u, TaoLatticeProxy& tao_lat);
extern "C" bool fortran_tao_lat_emit_calc(
    c_Int& plane /* 0D_NOT_integer */,
    c_Int& emit_type /* 0D_NOT_integer */,
    void* ele /* 0D_NOT_type */,
    void* modes /* 0D_NOT_type */,
    c_Real& emit /* 0D_NOT_real */);
void tao_lat_emit_calc(
    int plane,
    int emit_type,
    EleProxy& ele,
    NormalModesProxy& modes,
    double emit);
extern "C" bool fortran_tao_lat_sigma_calc_needed(
    c_Char data_type /* 0D_NOT_character */,
    c_Char data_source /* 0D_NOT_character */,
    c_Bool& do_lat_sigma /* 0D_NOT_logical */);
void tao_lat_sigma_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_lat_sigma);

// Skipped unusable routine tao_lat_sigma_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_lat_sigma_track(
    void* tao_lat /* 0D_NOT_type */,
    c_Bool& calc_ok /* 0D_NOT_logical */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_Bool* force_calc /* 0D_NOT_logical */);
bool tao_lat_sigma_track(
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> force_calc = std::nullopt);

// Skipped unusable routine tao_lattice_branch_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_lattice_branches_equal_tao_lattice_branches(
    void* tlb1 /* 1D_ALLOC_type */,
    void* tlb2 /* 1D_ALLOC_type */);
void tao_lattice_branches_equal_tao_lattice_branches(
    TaoLatticeBranchProxyAlloc1D& tlb1,
    TaoLatticeBranchProxyAlloc1D& tlb2);
extern "C" void fortran_tao_lattice_calc(
    c_Bool& calc_ok /* 0D_NOT_logical */,
    c_Bool& print_err /* 0D_NOT_logical */);
struct TaoLatticeCalc {
  bool calc_ok;
  bool print_err;
};
TaoLatticeCalc tao_lattice_calc();
extern "C" void fortran_tao_lattice_equal_tao_lattice(
    void* lat1 /* 0D_NOT_type */,
    void* lat2 /* 0D_NOT_type */);
void tao_lattice_equal_tao_lattice(
    TaoLatticeProxy& lat1,
    TaoLatticeProxy& lat2);

// Skipped unusable routine tao_lattice_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_limit_calc(c_Bool& limited /* 0D_NOT_logical */);
bool tao_limit_calc();
extern "C" void fortran_tao_lm_optimizer(c_Bool& abort /* 0D_NOT_logical */);
bool tao_lm_optimizer();
extern "C" void fortran_tao_lmdif_optimizer(c_Bool& abort /* 0D_NOT_logical */);
bool tao_lmdif_optimizer();
extern "C" void fortran_tao_load_this_datum(
    void* vec /* 1D_ALLOC_real */,
    void* ele_ref /* 0D_PTR_type */,
    void* ele_start /* 0D_PTR_type */,
    void* ele /* 0D_PTR_type */,
    c_Real& datum_value /* 0D_NOT_real */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    void* datum /* 0D_NOT_type */,
    void* branch /* 0D_NOT_type */,
    c_Char why_invalid /* 0D_NOT_character */,
    void* good /* 1D_ALLOC_logical */);
void tao_load_this_datum(
    RealAlloc1D& vec,
    EleProxy& ele_ref,
    EleProxy& ele_start,
    EleProxy& ele,
    double datum_value,
    bool valid_value,
    TaoDataProxy& datum,
    BranchProxy& branch,
    optional_ref<std::string> why_invalid = std::nullopt,
    optional_ref<BoolAlloc1D> good = std::nullopt);

// Skipped unusable routine tao_locate_all_elements:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine tao_locate_elements:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine tao_logical_array_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_mark_lattice_ele(void* lat /* 0D_NOT_type */);
void tao_mark_lattice_ele(LatProxy& lat);
extern "C" bool fortran_tao_merit(
    c_Bool& calc_ok /* 0D_NOT_logical */,
    c_Real& this_merit /* 0D_NOT_real */);
bool tao_merit(double this_merit);

// Skipped unusable routine tao_model_branch_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_model_element_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_mpi_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_mrq_func:
// Variable inout sized array: dy_da(:,:) 2D_NOT_real

// Skipped unusable routine tao_next_switch:
// Variable-sized in character array: switch_list(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_next_word(
    c_Char line /* 0D_NOT_character */,
    c_Char word /* 0D_NOT_character */);
std::string tao_next_word(std::string line);
extern "C" bool fortran_tao_one_turn_map_calc_needed(
    c_Char data_type /* 0D_NOT_character */,
    c_Char data_source /* 0D_NOT_character */,
    c_Bool& do_one_turn_map /* 0D_NOT_logical */);
void tao_one_turn_map_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_one_turn_map);
extern "C" void fortran_tao_open_file(
    c_Char file /* 0D_NOT_character */,
    c_Int& iunit /* 0D_NOT_integer */,
    c_Char file_name /* 0D_NOT_character */,
    c_Int& error_severity /* 0D_NOT_integer */,
    c_Bool* binary /* 0D_NOT_logical */);
int tao_open_file(
    std::string file,
    std::string file_name,
    int error_severity,
    std::optional<bool> binary = std::nullopt);
extern "C" bool fortran_tao_open_scratch_file(
    c_Bool& err /* 0D_NOT_logical */,
    c_Int& iu /* 0D_NOT_integer */);
bool tao_open_scratch_file(int iu);
extern "C" bool fortran_tao_optimization_status(
    void* datum /* 0D_NOT_type */,
    c_Char why_str /* 0D_NOT_character */);
void tao_optimization_status(TaoDataProxy& datum, std::string why_str);

// Skipped unusable routine tao_orbit_beta_wave_anal:
// Untranslated type: TaoPlotProxy (0D_NOT_type)

// Skipped unusable routine tao_oreint_building_wall_pt:
// Untranslated type: TaoBuildingWallPointProxy (0D_NOT_type)
// Untranslated type: TaoBuildingWallPointProxy (0D_NOT_type)
extern "C" bool fortran_tao_param_value_at_s(
    c_Char dat_name /* 0D_NOT_character */,
    void* ele_to_s /* 0D_NOT_type */,
    void* ele_here /* 0D_NOT_type */,
    void* orbit /* 0D_NOT_type */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Bool& print_err /* 0D_NOT_logical */,
    c_Bool& bad_datum /* 0D_NOT_logical */,
    c_Real& value /* 0D_NOT_real */);
struct TaoParamValueAtS {
  bool err_flag;
  std::string why_invalid;
  bool print_err;
  bool bad_datum;
};
TaoParamValueAtS tao_param_value_at_s(
    std::string dat_name,
    EleProxy& ele_to_s,
    EleProxy& ele_here,
    CoordProxy& orbit,
    double value);

// Skipped unusable routine tao_param_value_routine:
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)
extern "C" void fortran_tao_parse_command_args(
    c_Bool& error /* 0D_NOT_logical */,
    c_Char cmd_line /* 0D_NOT_character */);
bool tao_parse_command_args(optional_ref<std::string> cmd_line = std::nullopt);
extern "C" void fortran_tao_parse_element_param_str(
    c_Bool& err /* 0D_NOT_logical */,
    c_Char in_str /* 0D_NOT_character */,
    c_Char uni /* 0D_NOT_character */,
    c_Char element /* 0D_NOT_character */,
    c_Char parameter /* 0D_NOT_character */,
    c_Int& where /* 0D_NOT_integer */,
    c_Char component /* 0D_NOT_character */);
struct TaoParseElementParamStr {
  bool err;
  std::string uni;
  std::string element;
  std::string parameter;
  int where;
  std::string component;
};
TaoParseElementParamStr tao_parse_element_param_str(std::string in_str);
extern "C" void fortran_tao_particle_data_value(
    c_Char data_type /* 0D_NOT_character */,
    void* p /* 1D_ALLOC_type */,
    void* value /* 1D_ALLOC_real */,
    c_Bool& err /* 0D_NOT_logical */,
    void* ele /* 0D_NOT_type */,
    c_Int& ix_bunch /* 0D_NOT_integer */);
struct TaoParticleDataValue {
  RealAlloc1D value;
  bool err;
};
TaoParticleDataValue tao_particle_data_value(
    std::string data_type,
    CoordProxyAlloc1D& p,
    EleProxy& ele,
    int ix_bunch);
extern "C" void fortran_tao_pause_cmd(c_Real& time /* 0D_NOT_real */);
void tao_pause_cmd(double time);
extern "C" bool fortran_tao_phase_space_axis_index(
    c_Char data_type /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Int& ix_axis /* 0D_NOT_integer */);
int tao_phase_space_axis_index(std::string data_type, bool err);

// Skipped unusable routine tao_phase_wave_anal:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
extern "C" void fortran_tao_pick_universe(
    c_Char name_in /* 0D_NOT_character */,
    c_Char name_out /* 0D_NOT_character */,
    void* picked /* 1D_ALLOC_logical */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Int& ix_uni /* 0D_NOT_integer */,
    c_Bool& explicit_uni /* 0D_NOT_logical */,
    c_Int* dflt_uni /* 0D_NOT_integer */,
    c_Bool* pure_uni /* 0D_NOT_logical */);
struct TaoPickUniverse {
  std::string name_out;
  BoolAlloc1D picked;
  bool err;
  int ix_uni;
  bool explicit_uni;
};
TaoPickUniverse tao_pick_universe(
    std::string name_in,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<bool> pure_uni = std::nullopt);

// Skipped unusable routine tao_ping_scale_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_pipe_cmd(c_Char input_str /* 0D_NOT_character */);
void tao_pipe_cmd(std::string input_str);
extern "C" void fortran_tao_place_cmd(
    c_Char where /* 0D_NOT_character */,
    c_Char who /* 0D_NOT_character */,
    c_Bool* no_buffer /* 0D_NOT_logical */);
void tao_place_cmd(
    std::string where,
    std::string who,
    std::optional<bool> no_buffer = std::nullopt);

// Skipped unusable routine tao_place_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_plot_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_plot_cache_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_plot_cmd(
    c_Char where /* 0D_NOT_character */,
    c_Char component /* 0D_NOT_character */);
void tao_plot_cmd(std::string where, std::string component);

// Skipped unusable routine tao_plot_data:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_plot_histogram:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_plot_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_plot_key_table:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_plot_page_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_plot_page_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_plot_region_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_plot_setup();
void tao_plot_setup();

// Skipped unusable routine tao_plot_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_plot_struct_transfer:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoPlotProxy (0D_NOT_type)

// Skipped unusable routine tao_plot_wave:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_point_d1_to_data:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_point_v1_to_var:
// Untranslated type: TaoV1VarProxy (0D_NOT_type)
// Untranslated type: TaoVarProxy (1D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_pointer_to_branches:
// Untranslated type: BranchPointerProxy (1D_ALLOC_type)
// Untranslated type: TaoUniversePointerProxy (1D_ALLOC_type)

// Skipped unusable routine tao_pointer_to_building_wall_shape:
// Untranslated type: TaoEleShapeProxy (0D_PTR_type)
extern "C" bool fortran_tao_pointer_to_datum(
    void* d1 /* 0D_NOT_type */,
    c_Char ele_name /* 0D_NOT_character */,
    void* datum_ptr /* 0D_PTR_type */);
void tao_pointer_to_datum(
    TaoD1DataProxy& d1,
    std::string ele_name,
    TaoDataProxy& datum_ptr);
extern "C" bool fortran_tao_pointer_to_datum_ele(
    void* lat /* 0D_NOT_type */,
    c_Char ele_name /* 0D_NOT_character */,
    c_Int& ix_ele /* 0D_NOT_integer */,
    void* datum /* 0D_NOT_type */,
    c_Bool& valid /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Bool* print_err /* 0D_NOT_logical */,
    void* ele /* 0D_PTR_type */);
struct TaoPointerToDatumEle {
  bool valid;
  std::string why_invalid;
  EleProxy ele;
};
TaoPointerToDatumEle tao_pointer_to_datum_ele(
    LatProxy& lat,
    std::string ele_name,
    int ix_ele,
    TaoDataProxy& datum,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine tao_pointer_to_ele_shape:
// Untranslated type: TaoEleShapeProxy (1D_ALLOC_type)
// Untranslated type: TaoEleShapeProxy (0D_PTR_type)
extern "C" bool fortran_tao_pointer_to_tao_lat(
    void* u /* 0D_NOT_type */,
    c_Int* lat_type /* 0D_NOT_integer */,
    void* tao_lat /* 0D_PTR_type */);
void tao_pointer_to_tao_lat(
    TaoUniverseProxy& u,
    std::optional<int> lat_type,
    TaoLatticeProxy& tao_lat);

// Skipped unusable routine tao_pointer_to_universes:
// Untranslated type: TaoUniversePointerProxy (1D_ALLOC_type)

// Skipped unusable routine tao_pointer_to_var_in_lattice:
// Untranslated type: TaoVarProxy (0D_NOT_type)

// Skipped unusable routine tao_pointer_to_var_in_lattice2:
// Untranslated type: TaoVarProxy (0D_NOT_type)
extern "C" void fortran_tao_print_command_line_info();
void tao_print_command_line_info();

// Skipped unusable routine tao_print_vars:
// Untranslated type: TaoVarArrayProxy (1D_ALLOC_type)

// Skipped unusable routine tao_ptc_cmd:
// Module name unset? Internal error
extern "C" void fortran_tao_ptc_normal_form(
    c_Bool& do_calc /* 0D_NOT_logical */,
    void* tao_lat /* 0D_NOT_type */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Int* rf_on /* 0D_NOT_integer */);
void tao_ptc_normal_form(
    bool do_calc,
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<int> rf_on = std::nullopt);
extern "C" void fortran_tao_python_cmd(c_Char input_str /* 0D_NOT_character */);
void tao_python_cmd(std::string input_str);
extern "C" void fortran_tao_quiet_set(c_Char set /* 0D_NOT_character */);
void tao_quiet_set(std::string set);
extern "C" bool fortran_tao_rad_int_calc_needed(
    c_Char data_type /* 0D_NOT_character */,
    c_Char data_source /* 0D_NOT_character */,
    c_Bool& do_rad_int /* 0D_NOT_logical */);
void tao_rad_int_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_rad_int);

// Skipped unusable routine tao_re_allocate_expression_info:
// Untranslated type: TaoExpressionInfoProxy (1D_ALLOC_type)

// Skipped unusable routine tao_re_associate_node_array:
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)
extern "C" void fortran_tao_re_execute(
    c_Char string /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */);
void tao_re_execute(std::string string, bool err);
extern "C" void fortran_tao_read_cmd(
    c_Char which /* 0D_NOT_character */,
    c_Char unis /* 0D_NOT_character */,
    c_Char file /* 0D_NOT_character */,
    c_Bool& silent /* 0D_NOT_logical */);
void tao_read_cmd(
    std::string which,
    std::string unis,
    std::string file,
    bool silent);

// Skipped unusable routine tao_read_in_patterns:
// Module name unset? Internal error
extern "C" bool fortran_tao_read_phase_space_index(
    c_Char name /* 0D_NOT_character */,
    c_Int& ixc /* 0D_NOT_integer */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_Int& ix_ps /* 0D_NOT_integer */);
void tao_read_phase_space_index(
    std::string name,
    int ixc,
    std::optional<bool> print_err,
    int ix_ps);

// Skipped unusable routine tao_real_pointer_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_region_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_region_input_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_regression_test();
void tao_regression_test();
extern "C" void fortran_tao_remove_blank_characters(
    c_Char str /* 0D_NOT_character */);
void tao_remove_blank_characters(std::string str);
extern "C" void fortran_tao_run_cmd(
    c_Char which /* 0D_NOT_character */,
    c_Bool& abort /* 0D_NOT_logical */);
bool tao_run_cmd(std::string which);
extern "C" void fortran_tao_scale_cmd(
    c_Char where /* 0D_NOT_character */,
    c_Real& y_min_in /* 0D_NOT_real */,
    c_Real& y_max_in /* 0D_NOT_real */,
    c_Char axis /* 0D_NOT_character */,
    c_Bool* include_wall /* 0D_NOT_logical */,
    c_Char gang /* 0D_NOT_character */,
    c_Bool* exact /* 0D_NOT_logical */,
    c_Bool* turn_autoscale_off /* 0D_NOT_logical */);
void tao_scale_cmd(
    std::string where,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> exact = std::nullopt,
    std::optional<bool> turn_autoscale_off = std::nullopt);

// Skipped unusable routine tao_scale_graph:
// Untranslated type: TaoGraphProxy (0D_NOT_type)
extern "C" void fortran_tao_scale_ping_data(void* u /* 0D_NOT_type */);
void tao_scale_ping_data(TaoUniverseProxy& u);

// Skipped unusable routine tao_scale_plot:
// Untranslated type: TaoPlotProxy (0D_NOT_type)

// Skipped unusable routine tao_scratch_space_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_scratch_values_calc(
    void* ele_ref /* 0D_PTR_type */,
    void* ele_start /* 0D_PTR_type */,
    void* ele /* 0D_PTR_type */,
    void* datum /* 0D_NOT_type */,
    void* branch /* 0D_NOT_type */,
    void* orbit /* 1D_ALLOC_type */);
void tao_scratch_values_calc(
    EleProxy& ele_ref,
    EleProxy& ele_start,
    EleProxy& ele,
    TaoDataProxy& datum,
    BranchProxy& branch,
    CoordProxyAlloc1D& orbit);
extern "C" void fortran_tao_set_beam_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Char branch_str /* 0D_NOT_character */);
void tao_set_beam_cmd(
    std::string who,
    std::string value_str,
    std::string branch_str);
extern "C" void fortran_tao_set_beam_init_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Char branch_str /* 0D_NOT_character */);
void tao_set_beam_init_cmd(
    std::string who,
    std::string value_str,
    std::string branch_str);
extern "C" void fortran_tao_set_bmad_com_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_bmad_com_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_branch_cmd(
    c_Char branch_str /* 0D_NOT_character */,
    c_Char component_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_branch_cmd(
    std::string branch_str,
    std::string component_str,
    std::string value_str);
extern "C" void fortran_tao_set_calculate_cmd(
    c_Char switch_ /* 0D_NOT_character */);
void tao_set_calculate_cmd(optional_ref<std::string> switch_ = std::nullopt);
extern "C" void fortran_tao_set_curve_cmd(
    c_Char curve_name /* 0D_NOT_character */,
    c_Char component /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_curve_cmd(
    std::string curve_name,
    std::string component,
    std::string value_str);

// Skipped unusable routine tao_set_curve_invalid:
// Untranslated type: TaoCurveProxy (0D_NOT_type)
extern "C" void fortran_tao_set_data_cmd(
    c_Char who_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Bool* silent /* 0D_NOT_logical */);
void tao_set_data_cmd(
    std::string who_str,
    std::string value_str,
    optional_ref<bool> silent = std::nullopt);
extern "C" void fortran_tao_set_data_useit_opt(void* data /* 1D_ALLOC_type */);
void tao_set_data_useit_opt(
    optional_ref<TaoDataProxyAlloc1D> data = std::nullopt);
extern "C" void fortran_tao_set_default_cmd(
    c_Char who_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_default_cmd(std::string who_str, std::string value_str);

// Skipped unusable routine tao_set_drawing_cmd:
// Untranslated type: TaoDrawingProxy (0D_NOT_type)
extern "C" void fortran_tao_set_dynamic_aperture_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_dynamic_aperture_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_elements_cmd(
    c_Char ele_list /* 0D_NOT_character */,
    c_Char attribute /* 0D_NOT_character */,
    c_Char value /* 0D_NOT_character */,
    c_Bool& update /* 0D_NOT_logical */);
void tao_set_elements_cmd(
    std::string ele_list,
    std::string attribute,
    std::string value,
    bool update);

// Skipped unusable routine tao_set_flags_for_changed_attribute:
// Untranslated type: AllPointerProxy (0D_NOT_type)

// Skipped unusable routine tao_set_floor_plan_axis_label:
// Untranslated type: TaoGraphProxy (0D_NOT_type)
// Untranslated type: QpAxisProxy (0D_NOT_type)
// Untranslated type: QpAxisProxy (0D_NOT_type)
extern "C" void fortran_tao_set_geodesic_lm_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_geodesic_lm_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_global_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_global_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_graph_cmd(
    c_Char graph_name /* 0D_NOT_character */,
    c_Char component /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_graph_cmd(
    std::string graph_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_integer_value(
    c_Int& var /* 0D_NOT_integer */,
    c_Char var_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Bool& error /* 0D_NOT_logical */,
    c_Int* min_val /* 0D_NOT_integer */,
    c_Int* max_val /* 0D_NOT_integer */,
    c_Bool* print_err /* 0D_NOT_logical */);
struct TaoSetIntegerValue {
  int var;
  bool error;
};
TaoSetIntegerValue tao_set_integer_value(
    std::string var_str,
    std::string value_str,
    std::optional<int> min_val = std::nullopt,
    std::optional<int> max_val = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_set_invalid(
    void* datum /* 0D_NOT_type */,
    c_Char message /* 0D_NOT_character */,
    c_Char why_invalid /* 0D_NOT_character */,
    c_Bool* exterminate /* 0D_NOT_logical */,
    c_Int* err_level /* 0D_NOT_integer */,
    c_Bool* print_err /* 0D_NOT_logical */);
std::string tao_set_invalid(
    TaoDataProxy& datum,
    std::string message,
    std::optional<bool> exterminate = std::nullopt,
    std::optional<int> err_level = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_set_key_cmd(
    c_Char key_str /* 0D_NOT_character */,
    c_Char cmd_str /* 0D_NOT_character */);
void tao_set_key_cmd(std::string key_str, std::string cmd_str);
extern "C" void fortran_tao_set_lattice_cmd(
    c_Char dest_lat /* 0D_NOT_character */,
    c_Char source_lat /* 0D_NOT_character */);
void tao_set_lattice_cmd(std::string dest_lat, std::string source_lat);
extern "C" void fortran_tao_set_logical_value(
    c_Bool& var /* 0D_NOT_logical */,
    c_Char var_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Bool& error /* 0D_NOT_logical */);
struct TaoSetLogicalValue {
  bool var;
  bool error;
};
TaoSetLogicalValue tao_set_logical_value(
    std::string var_str,
    std::string value_str);
extern "C" void fortran_tao_set_openmp_n_threads(
    c_Int& n_threads /* 0D_NOT_integer */);
void tao_set_openmp_n_threads(int n_threads);
extern "C" void fortran_tao_set_opt_vars(
    void* var_vec /* 1D_ALLOC_real */,
    c_Bool* print_limit_warning /* 0D_NOT_logical */);
void tao_set_opt_vars(
    RealAlloc1D& var_vec,
    std::optional<bool> print_limit_warning = std::nullopt);
extern "C" void fortran_tao_set_opti_de_param_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_opti_de_param_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_particle_start_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_particle_start_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_plot_cmd(
    c_Char plot_name /* 0D_NOT_character */,
    c_Char component /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_plot_cmd(
    std::string plot_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_plot_page_cmd(
    c_Char component /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Char value_str2 /* 0D_NOT_character */);
void tao_set_plot_page_cmd(
    std::string component,
    std::string value_str,
    std::optional<std::string> value_str2 = std::nullopt);

// Skipped unusable routine tao_set_plotting:
// Untranslated type: TaoPlotPageInputProxy (0D_NOT_type)
// Untranslated type: TaoPlotPageProxy (0D_NOT_type)
extern "C" void fortran_tao_set_ptc_com_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_ptc_com_cmd(std::string who, std::string value_str);

// Skipped unusable routine tao_set_qp_axis_struct:
// Untranslated type: QpAxisProxy (0D_NOT_type)

// Skipped unusable routine tao_set_qp_point_struct:
// Untranslated type: QpPointProxy (0D_NOT_type)

// Skipped unusable routine tao_set_qp_rect_struct:
// Untranslated type: QpRectProxy (0D_NOT_type)
extern "C" void fortran_tao_set_ran_state_cmd(
    c_Char state_string /* 0D_NOT_character */);
void tao_set_ran_state_cmd(std::string state_string);
extern "C" void fortran_tao_set_real_value(
    c_Real& var /* 0D_NOT_real */,
    c_Char var_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Bool& error /* 0D_NOT_logical */,
    c_Real* min_val /* 0D_NOT_real */,
    c_Real* max_val /* 0D_NOT_real */,
    c_Int* dflt_uni /* 0D_NOT_integer */);
struct TaoSetRealValue {
  double var;
  bool error;
};
TaoSetRealValue tao_set_real_value(
    std::string var_str,
    std::string value_str,
    std::optional<double> min_val = std::nullopt,
    std::optional<double> max_val = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt);
extern "C" void fortran_tao_set_region_cmd(
    c_Char region_name /* 0D_NOT_character */,
    c_Char component /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_region_cmd(
    std::string region_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_space_charge_com_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_space_charge_com_cmd(std::string who, std::string value_str);

// Skipped unusable routine tao_set_switch_value:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_set_symbolic_number_cmd(
    c_Char sym_str /* 0D_NOT_character */,
    c_Char num_str /* 0D_NOT_character */,
    c_Real* val /* 0D_NOT_real */);
void tao_set_symbolic_number_cmd(
    std::string sym_str,
    std::optional<std::string> num_str = std::nullopt,
    std::optional<double> val = std::nullopt);
extern "C" void fortran_tao_set_tune_cmd(
    c_Char branch_str /* 0D_NOT_character */,
    c_Char mask_str /* 0D_NOT_character */,
    c_Bool& print_list /* 0D_NOT_logical */,
    c_Char qa_str /* 0D_NOT_character */,
    c_Char qb_str /* 0D_NOT_character */,
    c_Bool& delta_input /* 0D_NOT_logical */);
void tao_set_tune_cmd(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string qa_str,
    std::string qb_str,
    bool delta_input);
extern "C" void fortran_tao_set_universe_cmd(
    c_Char uni /* 0D_NOT_character */,
    c_Char who /* 0D_NOT_character */,
    c_Char what /* 0D_NOT_character */);
void tao_set_universe_cmd(std::string uni, std::string who, std::string what);
extern "C" void fortran_tao_set_var_cmd(
    c_Char var_str /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */);
void tao_set_var_cmd(std::string var_str, std::string value_str);

// Skipped unusable routine tao_set_var_model_value:
// Untranslated type: TaoVarProxy (0D_NOT_type)
extern "C" void fortran_tao_set_var_useit_opt();
void tao_set_var_useit_opt();
extern "C" void fortran_tao_set_wave_cmd(
    c_Char who /* 0D_NOT_character */,
    c_Char value_str /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */);
bool tao_set_wave_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_z_tune_cmd(
    c_Char branch_str /* 0D_NOT_character */,
    c_Char q_str /* 0D_NOT_character */,
    c_Bool& delta_input /* 0D_NOT_logical */);
void tao_set_z_tune_cmd(
    std::string branch_str,
    std::string q_str,
    bool delta_input);
extern "C" void fortran_tao_setup_key_table();
void tao_setup_key_table();

// Skipped unusable routine tao_shape_init:
// Untranslated type: TaoEleShapeProxy (0D_NOT_type)

// Skipped unusable routine tao_shape_pattern_point_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_shape_pattern_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_show_cmd(c_Char what /* 0D_NOT_character */);
void tao_show_cmd(std::string what);
extern "C" void fortran_tao_show_constraints(
    c_Int& iunit /* 0D_NOT_integer */,
    c_Char form /* 0D_NOT_character */);
void tao_show_constraints(int iunit, std::string form);

// Skipped unusable routine tao_show_this:
// Variable-sized out character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_single_mode(c_Char char_ /* 0D_NOT_character */);
void tao_single_mode(std::string char_);
extern "C" void fortran_tao_single_track(
    void* tao_lat /* 0D_NOT_type */,
    c_Bool& calc_ok /* 0D_NOT_logical */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Bool* print_err /* 0D_NOT_logical */);
bool tao_single_track(
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine tao_spin_dn_dpz_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_spin_ele_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_spin_map_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_spin_matrices_calc_needed(
    c_Char data_type /* 0D_NOT_character */,
    c_Char data_source /* 0D_NOT_character */,
    c_Bool& do_calc /* 0D_NOT_logical */);
void tao_spin_matrices_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_calc);

// Skipped unusable routine tao_spin_matrix_calc:
// Routine in configuration skip list

// Skipped unusable routine tao_spin_polarization_calc:
// Routine in configuration skip list

// Skipped unusable routine tao_spin_polarization_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_spin_tracking_turn_on();
void tao_spin_tracking_turn_on();

// Skipped unusable routine tao_split_component:
// Untranslated type: TaoDataVarComponentProxy (1D_ALLOC_type)
extern "C" bool fortran_tao_srdt_calc_needed(
    c_Char data_type /* 0D_NOT_character */,
    c_Char data_source /* 0D_NOT_character */,
    c_Int& do_srdt /* 0D_NOT_integer */);
void tao_srdt_calc_needed(
    std::string data_type,
    std::string data_source,
    int do_srdt);

// Skipped unusable routine tao_string_array_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_subin_uni_number(
    c_Char name_in /* 0D_NOT_character */,
    c_Int& ix_uni /* 0D_NOT_integer */,
    c_Char name_out /* 0D_NOT_character */,
    c_Bool& ok /* 0D_NOT_logical */);
std::string tao_subin_uni_number(std::string name_in, int ix_uni, bool ok);

// Skipped unusable routine tao_super_universe_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_svd_func:
// Variable inout sized array: dy_da(:,:) 2D_NOT_real
extern "C" void fortran_tao_svd_optimizer(c_Bool& abort /* 0D_NOT_logical */);
bool tao_svd_optimizer();
extern "C" void fortran_tao_symbol_import_from_lat(void* lat /* 0D_NOT_type */);
void tao_symbol_import_from_lat(LatProxy& lat);
extern "C" void fortran_tao_taper_cmd(
    c_Char except /* 0D_NOT_character */,
    c_Char uni_names /* 0D_NOT_character */);
void tao_taper_cmd(std::string except, std::string uni_names);

// Skipped unusable routine tao_timer:
// Module name unset? Internal error

// Skipped unusable routine tao_title_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_to_change_number(
    c_Char num_str /* 0D_NOT_character */,
    c_Int& n_size /* 0D_NOT_integer */,
    void* change_number /* 1D_ALLOC_real */,
    c_Char abs_or_rel /* 0D_NOT_character */,
    c_Bool& err /* 0D_NOT_logical */);
void tao_to_change_number(
    std::string num_str,
    int n_size,
    RealAlloc1D& change_number,
    std::string abs_or_rel,
    bool err);
extern "C" void fortran_tao_to_int(
    c_Char str /* 0D_NOT_character */,
    c_Int& i_int /* 0D_NOT_integer */,
    c_Bool& err /* 0D_NOT_logical */);
void tao_to_int(std::string str, int i_int, bool err);
extern "C" void fortran_tao_to_phase_and_coupling_reading(
    void* ele /* 0D_NOT_type */,
    void* bpm_data /* 0D_NOT_type */,
    c_Bool& valid_value /* 0D_NOT_logical */,
    c_Char why_invalid /* 0D_NOT_character */,
    void* datum /* 0D_NOT_type */);
struct TaoToPhaseAndCouplingReading {
  BpmPhaseCouplingProxy bpm_data;
  bool valid_value;
};
TaoToPhaseAndCouplingReading tao_to_phase_and_coupling_reading(
    EleProxy& ele,
    std::string why_invalid,
    TaoDataProxy& datum);
extern "C" void fortran_tao_to_real(
    c_Char expression /* 0D_NOT_character */,
    c_Real& value /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */);
struct TaoToReal {
  double value;
  bool err_flag;
};
TaoToReal tao_to_real(std::string expression);

// Skipped unusable routine tao_to_top10:
// Untranslated type: TaoTop10Proxy (1D_ALLOC_type)
extern "C" bool fortran_tao_too_many_particles_lost(
    void* beam /* 0D_NOT_type */,
    c_Bool& no_beam /* 0D_NOT_logical */);
void tao_too_many_particles_lost(BeamProxy& beam, bool no_beam);
extern "C" void fortran_tao_top10_derivative_print();
void tao_top10_derivative_print();
extern "C" void fortran_tao_top10_merit_categories_print(
    c_Int& iunit /* 0D_NOT_integer */);
void tao_top10_merit_categories_print(int iunit);

// Skipped unusable routine tao_top10_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_top_level(
    c_Char command /* 0D_NOT_character */,
    c_Int& errcode /* 0D_NOT_integer */);
int tao_top_level(std::optional<std::string> command = std::nullopt);
extern "C" bool fortran_tao_tracking_ele_index(
    void* ele /* 0D_PTR_type */,
    void* datum /* 0D_NOT_type */,
    c_Int& ix_branch /* 0D_NOT_integer */,
    c_Int& ix_ele /* 0D_NOT_integer */);
struct TaoTrackingEleIndex {
  int ix_branch;
  int ix_ele;
};
TaoTrackingEleIndex tao_tracking_ele_index(EleProxy& ele, TaoDataProxy& datum);
extern "C" void fortran_tao_turn_on_special_calcs_if_needed_for_plotting();
void tao_turn_on_special_calcs_if_needed_for_plotting();

// Skipped unusable routine tao_type_expression_tree:
// Untranslated type: TaoEvalNodeProxy (0D_NOT_type)
extern "C" bool fortran_tao_uni_atsign_index(
    c_Char string /* 0D_NOT_character */,
    c_Int& ix_amp /* 0D_NOT_integer */);
int tao_uni_atsign_index(std::string string);

// Skipped unusable routine tao_universe_calc_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" bool fortran_tao_universe_index(
    c_Int& i_uni /* 0D_NOT_integer */,
    c_Bool* neg2_to_default /* 0D_NOT_logical */,
    c_Int& i_this_uni /* 0D_NOT_integer */);
void tao_universe_index(
    int i_uni,
    std::optional<bool> neg2_to_default,
    int i_this_uni);

// Skipped unusable routine tao_universe_pointer_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_universe_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_use_data(
    c_Char action /* 0D_NOT_character */,
    c_Char data_name /* 0D_NOT_character */);
void tao_use_data(std::string action, std::string data_name);
extern "C" void fortran_tao_use_var(
    c_Char action /* 0D_NOT_character */,
    c_Char var_name /* 0D_NOT_character */);
void tao_use_var(std::string action, std::string var_name);
extern "C" bool fortran_tao_user_is_terminating_optimization(
    c_Bool& is_terminating /* 0D_NOT_logical */);
bool tao_user_is_terminating_optimization();

// Skipped unusable routine tao_v1_var_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_v1_var_input_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_v1_var_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_var1_name:
// Untranslated type: TaoVarProxy (0D_NOT_type)

// Skipped unusable routine tao_var_array_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_var_attrib_name:
// Untranslated type: TaoVarProxy (0D_NOT_type)

// Skipped unusable routine tao_var_check:
// Untranslated type: ElePointerProxy (1D_ALLOC_type)

// Skipped unusable routine tao_var_input_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_var_repoint();
void tao_var_repoint();

// Skipped unusable routine tao_var_show_use:
// Untranslated type: TaoV1VarProxy (0D_NOT_type)
// Variable-sized inout character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_var_slave_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_var_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_var_stuffit1:
// Untranslated type: TaoVarInputProxy (1D_NOT_type)
// Untranslated type: TaoV1VarProxy (0D_PTR_type)
// Untranslated type: TaoV1VarInputProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_var_stuffit2:
// Untranslated type: TaoVarProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_var_target_calc();
void tao_var_target_calc();

// Skipped unusable routine tao_var_useit_plot_calc:
// Untranslated type: TaoGraphProxy (0D_NOT_type)
// Untranslated type: TaoVarProxy (1D_ALLOC_type)
extern "C" void fortran_tao_var_write(
    c_Char out_file /* 0D_NOT_character */,
    c_Bool* show_good_opt_only /* 0D_NOT_logical */,
    c_Bool* tao_format /* 0D_NOT_logical */);
void tao_var_write(
    std::string out_file,
    std::optional<bool> show_good_opt_only = std::nullopt,
    std::optional<bool> tao_format = std::nullopt);
extern "C" void fortran_tao_veto_vars_with_zero_dmodel();
void tao_veto_vars_with_zero_dmodel();

// Skipped unusable routine tao_wave_analysis:
// Untranslated type: TaoPlotProxy (0D_NOT_type)
extern "C" void fortran_tao_wave_cmd(
    c_Char curve_name /* 0D_NOT_character */,
    c_Char plot_place /* 0D_NOT_character */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void tao_wave_cmd(
    std::string curve_name,
    std::string plot_place,
    bool err_flag);

// Skipped unusable routine tao_wave_fit:
// Untranslated type: TaoCurveProxy (0D_NOT_type)

// Skipped unusable routine tao_wave_kick_pt_struct_to_json:
// Routine module (tao_json) in configuration skip list

// Skipped unusable routine tao_wave_struct_to_json:
// Routine module (tao_json) in configuration skip list
extern "C" void fortran_tao_write_cmd(c_Char what /* 0D_NOT_character */);
void tao_write_cmd(std::string what);

// Skipped unusable routine tao_write_lines:
// Variable-sized in character array: line(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_x_axis_cmd(
    c_Char where /* 0D_NOT_character */,
    c_Char what /* 0D_NOT_character */);
void tao_x_axis_cmd(std::string where, std::string what);
extern "C" void fortran_tao_x_scale_cmd(
    c_Char where /* 0D_NOT_character */,
    c_Real& x_min_in /* 0D_NOT_real */,
    c_Real& x_max_in /* 0D_NOT_real */,
    c_Bool& err /* 0D_NOT_logical */,
    c_Bool* include_wall /* 0D_NOT_logical */,
    c_Char gang /* 0D_NOT_character */,
    c_Bool* exact /* 0D_NOT_logical */,
    c_Bool* turn_autoscale_off /* 0D_NOT_logical */);
bool tao_x_scale_cmd(
    std::string where,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> exact = std::nullopt,
    std::optional<bool> turn_autoscale_off = std::nullopt);

// Skipped unusable routine tao_x_scale_graph:
// Untranslated type: TaoGraphProxy (0D_NOT_type)

// Skipped unusable routine tao_x_scale_plot:
// Untranslated type: TaoPlotProxy (0D_NOT_type)

// Skipped unusable routine user_signal:
// Module name unset? Internal error
} // namespace Tao
