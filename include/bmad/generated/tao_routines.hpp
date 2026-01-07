#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"

using namespace Bmad;

namespace Tao {

// Skipped unusable routine avv:
// - Routine in configuration skip list

// Skipped unusable routine callback:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_integrate_max(
    int& ix_start /* 0D_NOT_integer inout */,
    int& ix_ele /* 0D_NOT_integer inout */,
    double& datum_value /* 0D_NOT_real inout */,
    int& ix_m /* 0D_NOT_integer inout */,
    void* branch /* 0D_NOT_type inout */,
    void* vec /* 1D_ALLOC_real inout */,
    void* datum /* 0D_NOT_type inout */);
void integrate_max(
    int& ix_start,
    int& ix_ele,
    double& datum_value,
    int& ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum);
extern "C" void fortran_integrate_min(
    int& ix_start /* 0D_NOT_integer inout */,
    int& ix_ele /* 0D_NOT_integer inout */,
    double& datum_value /* 0D_NOT_real inout */,
    int& ix_m /* 0D_NOT_integer inout */,
    void* branch /* 0D_NOT_type inout */,
    void* vec /* 1D_ALLOC_real inout */,
    void* datum /* 0D_NOT_type inout */);
void integrate_min(
    int& ix_start,
    int& ix_ele,
    double& datum_value,
    int& ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum);

// Skipped unusable routine jacobian:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_re_allocate_c_double(
    void* re /* 1D_ALLOC_real inout */,
    int& n /* 0D_NOT_integer in */,
    bool* exact /* 0D_NOT_logical in */,
    double* init_val /* 0D_NOT_real inout */);
void re_allocate_c_double(
    RealAlloc1D& re,
    int n,
    std::optional<bool> exact = std::nullopt,
    optional_ref<double> init_val = std::nullopt);
extern "C" void fortran_tao_abort_command_file(
    bool* force_abort /* 0D_NOT_logical in */);
void tao_abort_command_file(std::optional<bool> force_abort = std::nullopt);
extern "C" void fortran_tao_add_to_normal_mode_h_array(
    const char* h_str /* 0D_NOT_character in */,
    void* h_array /* 1D_ALLOC_type out */);
ResonanceHProxyAlloc1D tao_add_to_normal_mode_h_array(std::string h_str);
extern "C" void fortran_tao_alias_cmd(
    const char* alias /* 0D_NOT_character in */,
    const char* string /* 0D_NOT_character in */);
void tao_alias_cmd(std::string alias, std::string string);
extern "C" void fortran_tao_allocate_data_array(
    void* u /* 0D_NOT_type inout */,
    int& n_data /* 0D_NOT_integer inout */,
    bool* exact /* 0D_NOT_logical inout */);
void tao_allocate_data_array(
    TaoUniverseProxy& u,
    int& n_data,
    optional_ref<bool> exact = std::nullopt);
extern "C" void fortran_tao_allocate_v1_var(
    int& n_v1 /* 0D_NOT_integer inout */,
    bool& save_old /* 0D_NOT_logical inout */);
void tao_allocate_v1_var(int& n_v1, bool& save_old);
extern "C" void fortran_tao_allocate_var_array(
    int& n_var /* 0D_NOT_integer in */,
    bool& default_good_user /* 0D_NOT_logical inout */);
void tao_allocate_var_array(int n_var, bool& default_good_user);
extern "C" bool fortran_tao_beam_emit_calc(
    int& plane /* 0D_NOT_integer in */,
    int& emit_type /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* bunch_params /* 0D_NOT_type in */,
    double& emit /* 0D_NOT_real inout */);
void tao_beam_emit_calc(
    int plane,
    int emit_type,
    EleProxy& ele,
    BunchParamsProxy& bunch_params,
    double& emit);
extern "C" void fortran_tao_beam_track(
    void* u /* 0D_NOT_type in */,
    void* tao_lat /* 0D_NOT_type in */,
    int& ix_branch /* 0D_NOT_integer in */,
    void* beam /* 0D_NOT_type inout */,
    bool& calc_ok /* 0D_NOT_logical out */);
bool tao_beam_track(
    TaoUniverseProxy& u,
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    BeamProxy& beam);
extern "C" bool fortran_tao_beam_track_endpoint(
    const char* ele_id /* 0D_NOT_character in */,
    void* lat /* 0D_NOT_type in */,
    const char* branch_str /* 0D_NOT_character in */,
    const char* where /* 0D_NOT_character in */,
    void* u /* 0D_NOT_type in */,
    void* ele /* 0D_PTR_type inout */);
void tao_beam_track_endpoint(
    std::string ele_id,
    LatProxy& lat,
    std::string branch_str,
    std::string where,
    TaoUniverseProxy& u,
    EleProxy& ele);
extern "C" bool fortran_tao_branch_index(
    int& ix_branch /* 0D_NOT_integer in */,
    int& ix_this /* 0D_NOT_integer inout */);
void tao_branch_index(int ix_branch, int& ix_this);

// Skipped unusable routine tao_c_command:
// - Argument not defined: c_str (have: [])
// - Argument not defined: tao_c_command (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_beam_track_element:
// - Argument not defined: tao_c_get_beam_track_element (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_integer_array:
// - Argument not defined: tao_c_get_integer_array (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_get_real_array:
// - Argument not defined: tao_c_get_real_array (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_init_tao:
// - Argument not defined: c_str (have: [])
// - Argument not defined: tao_c_init_tao (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_integer_array_size:
// - Argument not defined: tao_c_integer_array_size (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_out_io_buffer_get_line:
// - Argument not defined: n (have: [])
// - Argument not defined: tao_c_out_io_buffer_get_line (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_c_out_io_buffer_num_lines:
// - Argument not defined: tao_c_out_io_buffer_num_lines (have: [])
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_c_out_io_buffer_reset();
void tao_c_out_io_buffer_reset();

// Skipped unusable routine tao_c_real_array_size:
// - Argument not defined: tao_c_real_array_size (have: [])
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_calc_data_at_s_pts(
    void* tao_lat /* 0D_NOT_type inout */,
    void* curve /* 0D_NOT_type inout */,
    double& comp_sign /* 0D_NOT_real inout */,
    void* good /* 1D_ALLOC_logical inout */);
void tao_calc_data_at_s_pts(
    TaoLatticeProxy& tao_lat,
    TaoCurveProxy& curve,
    double& comp_sign,
    BoolAlloc1D& good);

// Skipped unusable routine tao_call_cmd:
// - Variable-sized in character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_cbar_wave_anal(void* plot /* 0D_NOT_type inout */);
void tao_cbar_wave_anal(TaoPlotProxy& plot);
extern "C" void fortran_tao_change_ele(
    const char* ele_name /* 0D_NOT_character in */,
    const char* attrib_name /* 0D_NOT_character in */,
    const char* num_str /* 0D_NOT_character in */,
    bool& update /* 0D_NOT_logical inout */,
    bool& err_flag /* 0D_NOT_logical out */);
bool tao_change_ele(
    std::string ele_name,
    std::string attrib_name,
    std::string num_str,
    bool& update);
extern "C" void fortran_tao_change_tune(
    const char* branch_str /* 0D_NOT_character in */,
    const char* mask_str /* 0D_NOT_character in */,
    bool& print_list /* 0D_NOT_logical in */,
    const char* dqa_str /* 0D_NOT_character in */,
    const char* dqb_str /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool tao_change_tune(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string dqa_str,
    std::string dqb_str);
extern "C" void fortran_tao_change_var(
    const char* name /* 0D_NOT_character in */,
    const char* num_str /* 0D_NOT_character in */,
    bool& silent /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool tao_change_var(std::string name, std::string num_str, bool silent);
extern "C" void fortran_tao_change_z_tune(
    const char* branch_str /* 0D_NOT_character in */,
    const char* dq_str /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool tao_change_z_tune(std::string branch_str, std::string dq_str);
extern "C" bool fortran_tao_chrom_calc_needed(
    const char* data_type /* 0D_NOT_character inout */,
    const char* data_source /* 0D_NOT_character inout */,
    bool& do_chrom /* 0D_NOT_logical inout */);
void tao_chrom_calc_needed(
    std::string& data_type,
    std::string& data_source,
    bool& do_chrom);
extern "C" void fortran_tao_clear_cmd(
    const char* cmd_line /* 0D_NOT_character in */);
void tao_clear_cmd(std::string cmd_line);
extern "C" void fortran_tao_clip_cmd(
    bool& gang /* 0D_NOT_logical in */,
    const char* where /* 0D_NOT_character in */,
    double& value1 /* 0D_NOT_real inout */,
    double& value2 /* 0D_NOT_real inout */);
void tao_clip_cmd(bool gang, std::string where, double& value1, double& value2);
extern "C" void fortran_tao_close_command_file();
void tao_close_command_file();

// Skipped unusable routine tao_cmd_end_calc:
// - Module name unset
extern "C" void fortran_tao_cmd_history_record(
    const char* cmd /* 0D_NOT_character inout */);
void tao_cmd_history_record(std::string& cmd);

// Skipped unusable routine tao_cmd_split:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_command(
    const char* command_line /* 0D_NOT_character in */,
    bool& err /* 0D_NOT_logical inout */,
    bool& err_is_fatal /* 0D_NOT_logical out */);
bool tao_command(std::string command_line, bool& err);
extern "C" bool fortran_tao_constraint_type_name(
    void* datum /* 0D_NOT_type in */,
    const char* datum_name /* 0D_NOT_character inout */);
void tao_constraint_type_name(TaoDataProxy& datum, std::string& datum_name);
extern "C" void fortran_tao_control_tree_list(
    void* ele /* 0D_NOT_type in */,
    void* tree /* 1D_ALLOC_type in */);
void tao_control_tree_list(EleProxy& ele, ElePointerProxyAlloc1D& tree);
extern "C" void fortran_tao_count_strings(
    const char* string /* 0D_NOT_character in */,
    const char* pattern /* 0D_NOT_character in */,
    int& num /* 0D_NOT_integer out */);
int tao_count_strings(std::string string, std::string pattern);
extern "C" void fortran_tao_create_plot_window();
void tao_create_plot_window();
extern "C" void fortran_tao_curve_beam_ellipse_setup(
    void* curve /* 0D_NOT_type inout */);
void tao_curve_beam_ellipse_setup(TaoCurveProxy& curve);
extern "C" bool fortran_tao_curve_check_universe(
    void* curve /* 0D_NOT_type inout */,
    void* uni /* 0D_PTR_type in */,
    bool& is_ok /* 0D_NOT_logical out */);
bool tao_curve_check_universe(TaoCurveProxy& curve, TaoUniverseProxy& uni);
extern "C" void fortran_tao_curve_data_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */,
    void* curve /* 0D_NOT_type inout */);
void tao_curve_data_setup(
    TaoPlotProxy& plot,
    TaoGraphProxy& graph,
    TaoCurveProxy& curve);
extern "C" void fortran_tao_curve_datum_calc(
    void* eles /* 1D_ALLOC_type in */,
    void* plot /* 0D_NOT_type in */,
    void* curve /* 0D_NOT_type inout */,
    const char* who /* 0D_NOT_character in */);
void tao_curve_datum_calc(
    ElePointerProxyAlloc1D& eles,
    TaoPlotProxy& plot,
    TaoCurveProxy& curve,
    std::string who);
extern "C" bool fortran_tao_curve_ele_ref(
    void* curve /* 0D_NOT_type in */,
    bool& point_to_ele_ref /* 0D_NOT_logical inout */,
    void* ele_track /* 0D_PTR_type inout */);
void tao_curve_ele_ref(
    TaoCurveProxy& curve,
    bool& point_to_ele_ref,
    EleProxy& ele_track);
extern "C" bool fortran_tao_curve_ix_uni(
    void* curve /* 0D_NOT_type in */,
    int& ix_uni /* 0D_NOT_integer inout */);
void tao_curve_ix_uni(TaoCurveProxy& curve, int& ix_uni);
extern "C" bool fortran_tao_curve_name(
    void* curve /* 0D_NOT_type in */,
    bool* use_region /* 0D_NOT_logical in */,
    const char* curve_name /* 0D_NOT_character inout */);
void tao_curve_name(
    TaoCurveProxy& curve,
    std::optional<bool> use_region,
    std::string& curve_name);
extern "C" void fortran_tao_curve_rms_calc(
    void* curve /* 0D_NOT_type in */,
    const char* who /* 0D_NOT_character in */,
    double& rms /* 0D_NOT_real out */,
    double& mean /* 0D_NOT_real out */);
struct TaoCurveRmsCalc {
  double rms;
  double mean;
};
Tao::TaoCurveRmsCalc tao_curve_rms_calc(TaoCurveProxy& curve, std::string who);
extern "C" bool fortran_tao_d2_d1_name(
    void* d1 /* 0D_NOT_type in */,
    bool* show_universe /* 0D_NOT_logical in */,
    const char* d2_d1_name /* 0D_NOT_character inout */);
void tao_d2_d1_name(
    TaoD1DataProxy& d1,
    std::optional<bool> show_universe,
    std::string& d2_d1_name);
extern "C" void fortran_tao_d2_data_stuffit(
    void* u /* 0D_NOT_type inout */,
    const char* d2_name /* 0D_NOT_character inout */,
    int& n_d1_data /* 0D_NOT_integer inout */);
void tao_d2_data_stuffit(
    TaoUniverseProxy& u,
    std::string& d2_name,
    int& n_d1_data);
extern "C" void fortran_tao_data_check(bool& err /* 0D_NOT_logical inout */);
void tao_data_check(bool& err);
extern "C" void fortran_tao_data_coupling_init(
    void* branch /* 0D_NOT_type in */);
void tao_data_coupling_init(BranchProxy& branch);
extern "C" bool fortran_tao_data_sanity_check(
    void* datum /* 0D_NOT_type in */,
    bool& print_err /* 0D_NOT_logical in */,
    const char* default_data_type /* 0D_NOT_character in */,
    void* uni /* 0D_NOT_type in */,
    bool& is_valid /* 0D_NOT_logical inout */);
void tao_data_sanity_check(
    TaoDataProxy& datum,
    bool print_err,
    std::string default_data_type,
    optional_ref<TaoUniverseProxy> uni,
    bool& is_valid);

// Skipped unusable routine tao_data_show_use:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_data_type_substitute(
    const char* template_ /* 0D_NOT_character in */,
    const char* str_out /* 0D_NOT_character out */,
    void* curve /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
std::string tao_data_type_substitute(
    std::string template_,
    TaoCurveProxy& curve,
    TaoGraphProxy& graph);
extern "C" void fortran_tao_data_useit_plot_calc(
    void* curve /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */,
    void* data /* 1D_ALLOC_type out */,
    bool& check_s_position /* 0D_NOT_logical in */,
    const char* most_invalid /* 0D_NOT_character out */);
struct TaoDataUseitPlotCalc {
  TaoDataProxyAlloc1D data;
  std::string most_invalid;
};
Tao::TaoDataUseitPlotCalc tao_data_useit_plot_calc(
    TaoCurveProxy& curve,
    TaoGraphProxy& graph,
    bool check_s_position);
extern "C" bool fortran_tao_datum_has_associated_ele(
    const char* data_type /* 0D_NOT_character in */,
    int* branch_geometry /* 0D_NOT_integer in */,
    int& has_associated_ele /* 0D_NOT_integer inout */);
void tao_datum_has_associated_ele(
    std::string data_type,
    std::optional<int> branch_geometry,
    int& has_associated_ele);
extern "C" bool fortran_tao_datum_integrate(
    void* datum /* 0D_NOT_type in */,
    void* branch /* 0D_NOT_type in */,
    void* s_pos /* 1D_ALLOC_real in */,
    void* values /* 1D_ALLOC_real in */,
    bool& valid_value /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character out */,
    double& result /* 0D_NOT_real out */);
struct TaoDatumIntegrate {
  bool valid_value;
  std::string why_invalid;
  double result;
};
Tao::TaoDatumIntegrate tao_datum_integrate(
    TaoDataProxy& datum,
    BranchProxy& branch,
    RealAlloc1D& s_pos,
    RealAlloc1D& values);
extern "C" bool fortran_tao_datum_name(
    void* datum /* 0D_NOT_type in */,
    bool* show_universe /* 0D_NOT_logical in */,
    const char* datum_name /* 0D_NOT_character inout */);
void tao_datum_name(
    TaoDataProxy& datum,
    std::optional<bool> show_universe,
    std::string& datum_name);
extern "C" bool fortran_tao_datum_s_position(
    void* datum /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    double& s_pos /* 0D_NOT_real out */);
double tao_datum_s_position(TaoDataProxy& datum, EleProxy& ele);
extern "C" void fortran_tao_de_optimizer(bool& abort /* 0D_NOT_logical out */);
bool tao_de_optimizer();
extern "C" void fortran_tao_deallocate_plot_cache(
    void* plot_cache /* 1D_ALLOC_type inout */);
void tao_deallocate_plot_cache(TaoPlotCacheProxyAlloc1D& plot_cache);

// Skipped unusable routine tao_deallocate_tree:
// - Untranslated type: tao_eval_node_struct (0D)
extern "C" void fortran_tao_destroy_plot_window();
void tao_destroy_plot_window();
extern "C" void fortran_tao_dmerit_calc();
void tao_dmerit_calc();
extern "C" void fortran_tao_dmodel_dvar_calc(
    bool& force_calc /* 0D_NOT_logical in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool tao_dmodel_dvar_calc(bool force_calc);
extern "C" bool fortran_tao_do_wire_scan(
    void* ele /* 0D_NOT_type in */,
    double& theta /* 0D_NOT_real in */,
    void* beam /* 0D_NOT_type in */,
    double& moment /* 0D_NOT_real out */);
double tao_do_wire_scan(EleProxy& ele, double theta, BeamProxy& beam);
extern "C" void fortran_tao_draw_beam_chamber_wall(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_draw_beam_chamber_wall(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_draw_curve_data(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */,
    void* curve /* 0D_NOT_type in */,
    bool& have_data /* 0D_NOT_logical inout */);
void tao_draw_curve_data(
    TaoPlotProxy& plot,
    TaoGraphProxy& graph,
    TaoCurveProxy& curve,
    bool& have_data);
extern "C" void fortran_tao_draw_ele_for_floor_plan(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */,
    void* tao_lat /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    void* ele_shape /* 0D_PTR_type in */,
    const char* label_name /* 0D_NOT_character in */,
    double& offset1 /* 0D_NOT_real inout */,
    double& offset2 /* 0D_NOT_real inout */);
void tao_draw_ele_for_floor_plan(
    TaoPlotProxy& plot,
    TaoGraphProxy& graph,
    TaoLatticeProxy& tao_lat,
    EleProxy& ele,
    TaoEleShapeProxy& ele_shape,
    std::string label_name,
    double& offset1,
    double& offset2);
extern "C" void fortran_tao_draw_floor_plan(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_draw_floor_plan(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_draw_graph_axes(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_draw_graph_axes(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_draw_histogram_data(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */,
    void* curve /* 0D_NOT_type in */,
    bool& have_data /* 0D_NOT_logical inout */);
void tao_draw_histogram_data(
    TaoPlotProxy& plot,
    TaoGraphProxy& graph,
    TaoCurveProxy& curve,
    bool& have_data);
extern "C" void fortran_tao_draw_lat_layout(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_draw_lat_layout(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_draw_plots(bool* do_clear /* 0D_NOT_logical in */);
void tao_draw_plots(std::optional<bool> do_clear = std::nullopt);
extern "C" bool fortran_tao_ele_geometry_with_misalignments(
    void* datum /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    bool& valid_value /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character out */,
    double& value /* 0D_NOT_real out */);
struct TaoEleGeometryWithMisalignments {
  bool valid_value;
  std::string why_invalid;
  double value;
};
Tao::TaoEleGeometryWithMisalignments tao_ele_geometry_with_misalignments(
    TaoDataProxy& datum,
    EleProxy& ele);
extern "C" void fortran_tao_ele_shape_info(
    int& ix_uni /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* ele_shapes /* 1D_ALLOC_type in */,
    void* e_shape /* 0D_PTR_type out */,
    const char* label_name /* 0D_NOT_character out */,
    double& y1 /* 0D_NOT_real inout */,
    double& y2 /* 0D_NOT_real inout */,
    int* ix_shape_min /* 0D_NOT_integer inout */);
struct TaoEleShapeInfo {
  TaoEleShapeProxy e_shape;
  std::string label_name;
};
Tao::TaoEleShapeInfo tao_ele_shape_info(
    int ix_uni,
    EleProxy& ele,
    TaoEleShapeProxyAlloc1D& ele_shapes,
    double& y1,
    double& y2,
    optional_ref<int> ix_shape_min = std::nullopt);

// Skipped unusable routine tao_ele_shape_input_to_struct:
// - Untranslated type: tao_ele_shape_input (0D)

// Skipped unusable routine tao_ele_shape_struct_to_input:
// - Untranslated type: tao_ele_shape_input (0D)
extern "C" bool fortran_tao_eval_floor_orbit(
    void* datum /* 0D_NOT_type in */,
    void* ele /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    void* bunch_params /* 0D_NOT_type in */,
    bool& valid_value /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character out */,
    double& value /* 0D_NOT_real out */);
struct TaoEvalFloorOrbit {
  bool valid_value;
  std::string why_invalid;
  double value;
};
Tao::TaoEvalFloorOrbit tao_eval_floor_orbit(
    TaoDataProxy& datum,
    EleProxy& ele,
    CoordProxy& orbit,
    BunchParamsProxy& bunch_params);
extern "C" void fortran_tao_evaluate_a_datum(
    void* datum /* 0D_NOT_type inout */,
    void* u /* 0D_NOT_type in */,
    void* tao_lat /* 0D_NOT_type in */,
    double& datum_value /* 0D_NOT_real out */,
    bool& valid_value /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character out */,
    bool* called_from_lat_calc /* 0D_NOT_logical in */,
    bool* print_err /* 0D_NOT_logical in */);
struct TaoEvaluateADatum {
  double datum_value;
  bool valid_value;
  std::string why_invalid;
};
Tao::TaoEvaluateADatum tao_evaluate_a_datum(
    TaoDataProxy& datum,
    TaoUniverseProxy& u,
    TaoLatticeProxy& tao_lat,
    std::optional<bool> called_from_lat_calc = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" bool fortran_tao_evaluate_datum_at_s(
    void* datum /* 0D_NOT_type in */,
    void* tao_lat /* 0D_NOT_type in */,
    void* ele /* 0D_PTR_type in */,
    void* ele_ref /* 0D_PTR_type in */,
    bool& valid_value /* 0D_NOT_logical in */,
    const char* err_str /* 0D_NOT_character out */,
    bool& bad_datum /* 0D_NOT_logical out */,
    double& value /* 0D_NOT_real out */);
struct TaoEvaluateDatumAtS {
  std::string err_str;
  bool bad_datum;
  double value;
};
Tao::TaoEvaluateDatumAtS tao_evaluate_datum_at_s(
    TaoDataProxy& datum,
    TaoLatticeProxy& tao_lat,
    EleProxy& ele,
    EleProxy& ele_ref,
    bool valid_value);

// Skipped unusable routine tao_evaluate_element_parameters:
// - Untranslated type: tao_expression_info_struct (1D)

// Skipped unusable routine tao_evaluate_expression:
// - Untranslated type: tao_expression_info_struct (1D)
// - Untranslated type: tao_eval_node_struct (1D)

// Skipped unusable routine tao_evaluate_expression_new:
// - Untranslated type: tao_expression_info_struct (1D)
// - Untranslated type: tao_eval_node_struct (1D)

// Skipped unusable routine tao_evaluate_expression_old:
// - Untranslated type: tao_expression_info_struct (1D)
// - Untranslated type: tao_eval_node_struct (1D)
extern "C" void fortran_tao_evaluate_lat_or_beam_data(
    bool& err /* 0D_NOT_logical out */,
    const char* data_name /* 0D_NOT_character in */,
    void* values /* 1D_ALLOC_real out */,
    bool& print_err /* 0D_NOT_logical in */,
    const char* default_source /* 0D_NOT_character inout */,
    void* dflt_ele_ref /* 0D_PTR_type in */,
    void* dflt_ele_start /* 0D_PTR_type in */,
    void* dflt_ele /* 0D_PTR_type in */,
    const char* dflt_component /* 0D_NOT_character in */,
    int* dflt_uni /* 0D_NOT_integer in */,
    int* dflt_eval_point /* 0D_NOT_integer in */,
    double* dflt_s_offset /* 0D_NOT_real in */);
struct TaoEvaluateLatOrBeamData {
  bool err;
  RealAlloc1D values;
};
Tao::TaoEvaluateLatOrBeamData tao_evaluate_lat_or_beam_data(
    std::string data_name,
    bool print_err,
    std::string& default_source,
    optional_ref<EleProxy> dflt_ele_ref = std::nullopt,
    optional_ref<EleProxy> dflt_ele_start = std::nullopt,
    optional_ref<EleProxy> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt);

// Skipped unusable routine tao_evaluate_stack_old:
// - Untranslated type: tao_eval_node_struct (1D)
// - Untranslated type: tao_expression_info_struct (1D)

// Skipped unusable routine tao_evaluate_tree:
// - Untranslated type: tao_eval_node_struct (0D)
// - Untranslated type: tao_expression_info_struct (1D)
extern "C" bool fortran_tao_evaluate_tune(
    const char* q_str /* 0D_NOT_character in */,
    double& q0 /* 0D_NOT_real in */,
    bool& delta_input /* 0D_NOT_logical in */,
    double& q_val /* 0D_NOT_real inout */);
void tao_evaluate_tune(
    std::string q_str,
    double q0,
    bool delta_input,
    double& q_val);
extern "C" void fortran_tao_expression_hash_substitute(
    const char* expression_in /* 0D_NOT_character in */,
    const char* expression_out /* 0D_NOT_character out */,
    void* eval_ele /* 0D_PTR_type in */);
std::string tao_expression_hash_substitute(
    std::string expression_in,
    optional_ref<EleProxy> eval_ele = std::nullopt);

// Skipped unusable routine tao_expression_tree_to_string:
// - Untranslated type: tao_eval_node_struct (0D)
// - Untranslated type: tao_eval_node_struct (0D)

// Skipped unusable routine tao_find_data:
// - Untranslated type: tao_d2_data_array_struct (1D)
// - Untranslated type: tao_d1_data_array_struct (1D)
// - Untranslated type: tao_data_array_struct (1D)
// - Untranslated type: tao_real_pointer_struct (1D)
// - Untranslated type: tao_logical_array_struct (1D)
// - Untranslated type: tao_string_array_struct (1D)
// - Untranslated type: tao_integer_array_struct (1D)
extern "C" void fortran_tao_find_plot_region(
    bool& err /* 0D_NOT_logical out */,
    const char* where /* 0D_NOT_character in */,
    void* region /* 0D_PTR_type out */,
    bool* print_flag /* 0D_NOT_logical in */);
struct TaoFindPlotRegion {
  bool err;
  TaoPlotRegionProxy region;
};
Tao::TaoFindPlotRegion tao_find_plot_region(
    std::string where,
    std::optional<bool> print_flag = std::nullopt);

// Skipped unusable routine tao_find_plots:
// - Untranslated type: tao_plot_array_struct (1D)
// - Untranslated type: tao_graph_array_struct (1D)
// - Untranslated type: tao_curve_array_struct (1D)

// Skipped unusable routine tao_find_var:
// - Untranslated type: tao_v1_var_array_struct (1D)
// - Untranslated type: tao_var_array_struct (1D)
// - Untranslated type: tao_real_pointer_struct (1D)
// - Untranslated type: tao_logical_array_struct (1D)
// - Untranslated type: tao_string_array_struct (1D)
extern "C" void fortran_tao_fixer(
    const char* switch_ /* 0D_NOT_character in */,
    const char* word1 /* 0D_NOT_character in */,
    const char* word2 /* 0D_NOT_character in */);
void tao_fixer(std::string switch_, std::string word1, std::string word2);
extern "C" void fortran_tao_floor_to_screen(
    void* graph /* 0D_NOT_type in */,
    double* r_floor /* 1D_NOT_real inout */,
    double& x_screen /* 0D_NOT_real out */,
    double& y_screen /* 0D_NOT_real out */);
struct TaoFloorToScreen {
  double x_screen;
  double y_screen;
};
Tao::TaoFloorToScreen tao_floor_to_screen(
    TaoGraphProxy& graph,
    FixedArray1D<Real, 3> r_floor);
extern "C" void fortran_tao_floor_to_screen_coords(
    void* graph /* 0D_NOT_type in */,
    void* floor /* 0D_NOT_type in */,
    void* screen /* 0D_NOT_type out */);
FloorPositionProxy tao_floor_to_screen_coords(
    TaoGraphProxy& graph,
    FloorPositionProxy& floor);

// Skipped unusable routine tao_geo_lm_func:
// - Routine in configuration skip list
extern "C" void fortran_tao_geodesic_lm_optimizer(
    bool& abort /* 0D_NOT_logical out */);
bool tao_geodesic_lm_optimizer();
extern "C" void fortran_tao_get_data(
    void* data_value /* 1D_ALLOC_real out */,
    void* data_weight /* 1D_ALLOC_real out */,
    void* data_meas_value /* 1D_ALLOC_real out */,
    void* data_ix_dModel /* 1D_ALLOC_integer out */);
struct TaoGetData {
  RealAlloc1D data_value;
  RealAlloc1D data_weight;
  RealAlloc1D data_meas_value;
  IntAlloc1D data_ix_dModel;
};
Tao::TaoGetData tao_get_data();
extern "C" void fortran_tao_get_opt_vars(
    void* var_value /* 1D_ALLOC_real out */,
    void* var_step /* 1D_ALLOC_real out */,
    void* var_delta /* 1D_ALLOC_real out */,
    void* var_weight /* 1D_ALLOC_real out */,
    void* var_ix /* 1D_ALLOC_integer out */,
    bool& ignore_if_weight_is_zero /* 0D_NOT_logical out */,
    bool& ignore_if_not_limited /* 0D_NOT_logical out */);
struct TaoGetOptVars {
  RealAlloc1D var_value;
  RealAlloc1D var_step;
  RealAlloc1D var_delta;
  RealAlloc1D var_weight;
  IntAlloc1D var_ix;
  bool ignore_if_weight_is_zero;
  bool ignore_if_not_limited;
};
Tao::TaoGetOptVars tao_get_opt_vars();
extern "C" void fortran_tao_get_user_input(
    const char* cmd_out /* 0D_NOT_character out */,
    const char* prompt_str /* 0D_NOT_character in */,
    bool* wait_flag /* 0D_NOT_logical in */,
    const char* cmd_in /* 0D_NOT_character in */);
std::string tao_get_user_input(
    std::optional<std::string> prompt_str = std::nullopt,
    std::optional<bool> wait_flag = std::nullopt,
    std::optional<std::string> cmd_in = std::nullopt);
extern "C" void fortran_tao_graph_controller_setup(
    void* graph /* 0D_NOT_type inout */);
void tao_graph_controller_setup(TaoGraphProxy& graph);
extern "C" void fortran_tao_graph_data_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */);
void tao_graph_data_setup(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_graph_data_slice_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */);
void tao_graph_data_slice_setup(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_graph_dynamic_aperture_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */);
void tao_graph_dynamic_aperture_setup(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_graph_histogram_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */);
void tao_graph_histogram_setup(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" bool fortran_tao_graph_name(
    void* graph /* 0D_NOT_type in */,
    bool* use_region /* 0D_NOT_logical in */,
    const char* graph_name /* 0D_NOT_character inout */);
void tao_graph_name(
    TaoGraphProxy& graph,
    std::optional<bool> use_region,
    std::string& graph_name);
extern "C" void fortran_tao_graph_phase_space_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */);
void tao_graph_phase_space_setup(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_graph_s_min_max_calc(
    void* graph /* 0D_NOT_type in */,
    void* branch /* 0D_NOT_type in */,
    double& s_min /* 0D_NOT_real out */,
    double& s_max /* 0D_NOT_real out */);
struct TaoGraphSMinMaxCalc {
  double s_min;
  double s_max;
};
Tao::TaoGraphSMinMaxCalc tao_graph_s_min_max_calc(
    TaoGraphProxy& graph,
    BranchProxy& branch);
extern "C" void fortran_tao_graph_setup(
    void* plot /* 0D_NOT_type inout */,
    void* graph /* 0D_NOT_type inout */);
void tao_graph_setup(TaoPlotProxy& plot, TaoGraphProxy& graph);

// Skipped unusable routine tao_help:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_hook_branch_calc_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_command_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_curve_s_pt_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_data_sanity_check_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_draw_floor_plan_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_draw_graph_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_evaluate_a_datum_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_graph_postsetup_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_graph_setup_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init1_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init2_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_beam_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_data_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_global_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_lattice_post_parse_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_plotting_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_read_lattice_info_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_init_var_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_lattice_calc_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_merit_data_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_merit_var_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_optimizer_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_parse_command_args_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_plot_setup_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_post_process_data_def:
// - Routine in configuration skip list

// Skipped unusable routine tao_hook_show_cmd_def:
// - Routine in configuration skip list
extern "C" void fortran_tao_init(bool& err_flag /* 0D_NOT_logical out */);
bool tao_init();
extern "C" void fortran_tao_init_beam_in_universe(
    void* u /* 0D_NOT_type inout */,
    void* beam_init /* 0D_NOT_type inout */,
    const char* track_start /* 0D_NOT_character inout */,
    const char* track_end /* 0D_NOT_character inout */,
    double& comb_ds_save /* 0D_NOT_real inout */);
void tao_init_beam_in_universe(
    TaoUniverseProxy& u,
    BeamInitProxy& beam_init,
    std::string& track_start,
    std::string& track_end,
    double& comb_ds_save);
extern "C" void fortran_tao_init_beams(
    const char* init_file /* 0D_NOT_character in */);
void tao_init_beams(std::string init_file);

// Skipped unusable routine tao_init_building_wall:
// - Module name unset
extern "C" void fortran_tao_init_data(
    const char* data_file /* 0D_NOT_character in */);
void tao_init_data(std::string data_file);
extern "C" void fortran_tao_init_data_end_stuff();
void tao_init_data_end_stuff();
extern "C" void fortran_tao_init_data_in_universe(
    void* u /* 0D_NOT_type inout */,
    int& n_d2_add /* 0D_NOT_integer inout */,
    bool* keep_existing_data /* 0D_NOT_logical inout */);
void tao_init_data_in_universe(
    TaoUniverseProxy& u,
    int& n_d2_add,
    optional_ref<bool> keep_existing_data = std::nullopt);
extern "C" void fortran_tao_init_dynamic_aperture(
    const char* init_file /* 0D_NOT_character in */);
void tao_init_dynamic_aperture(std::string init_file);
extern "C" void fortran_tao_init_find_elements(
    void* u /* 0D_NOT_type in */,
    const char* search_string /* 0D_NOT_character in */,
    void* eles /* 1D_ALLOC_type out */,
    const char* attribute /* 0D_NOT_character in */,
    bool& found_one /* 0D_NOT_logical out */);
struct TaoInitFindElements {
  ElePointerProxyAlloc1D eles;
  bool found_one;
};
Tao::TaoInitFindElements tao_init_find_elements(
    TaoUniverseProxy& u,
    std::string search_string,
    std::optional<std::string> attribute = std::nullopt);
extern "C" void fortran_tao_init_global(
    const char* init_file /* 0D_NOT_character in */);
void tao_init_global(std::string init_file);
extern "C" void fortran_tao_init_lattice(
    const char* lat_file /* 0D_NOT_character inout */,
    bool& err_flag /* 0D_NOT_logical inout */);
void tao_init_lattice(std::string& lat_file, bool& err_flag);
extern "C" void fortran_tao_init_plotting(
    const char* plot_file /* 0D_NOT_character inout */);
void tao_init_plotting(std::string& plot_file);
extern "C" void fortran_tao_init_variables(
    const char* var_file /* 0D_NOT_character in */);
void tao_init_variables(std::string var_file);
extern "C" void fortran_tao_inject_beam(
    void* u /* 0D_NOT_type in */,
    void* model /* 0D_NOT_type in */,
    int& ix_branch /* 0D_NOT_integer in */,
    void* beam /* 0D_NOT_type out */,
    bool& init_ok /* 0D_NOT_logical out */);
struct TaoInjectBeam {
  BeamProxy beam;
  bool init_ok;
};
Tao::TaoInjectBeam tao_inject_beam(
    TaoUniverseProxy& u,
    TaoLatticeProxy& model,
    int ix_branch);
extern "C" void fortran_tao_inject_particle(
    void* u /* 0D_NOT_type inout */,
    void* model /* 0D_NOT_type inout */,
    int& ix_branch /* 0D_NOT_integer inout */);
void tao_inject_particle(
    TaoUniverseProxy& u,
    TaoLatticeProxy& model,
    int& ix_branch);
extern "C" bool fortran_tao_is_valid_name(
    const char* name /* 0D_NOT_character in */,
    const char* why_invalid /* 0D_NOT_character out */,
    bool& is_valid /* 0D_NOT_logical inout */);
std::string tao_is_valid_name(std::string name, bool& is_valid);
extern "C" void fortran_tao_json_cmd(
    const char* input_str /* 0D_NOT_character in */);
void tao_json_cmd(std::string input_str);
extern "C" void fortran_tao_key_info_to_str(
    int& ix_key /* 0D_NOT_integer inout */,
    int& ix_min_key /* 0D_NOT_integer inout */,
    int& ix_max_key /* 0D_NOT_integer inout */,
    const char* key_str /* 0D_NOT_character inout */,
    const char* header_str /* 0D_NOT_character inout */);
void tao_key_info_to_str(
    int& ix_key,
    int& ix_min_key,
    int& ix_max_key,
    std::string& key_str,
    std::string& header_str);
extern "C" void fortran_tao_lat_bookkeeper(
    void* u /* 0D_NOT_type in */,
    void* tao_lat /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */);
bool tao_lat_bookkeeper(TaoUniverseProxy& u, TaoLatticeProxy& tao_lat);
extern "C" bool fortran_tao_lat_emit_calc(
    int& plane /* 0D_NOT_integer in */,
    int& emit_type /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* modes /* 0D_NOT_type in */,
    double& emit /* 0D_NOT_real inout */);
void tao_lat_emit_calc(
    int plane,
    int emit_type,
    EleProxy& ele,
    NormalModesProxy& modes,
    double& emit);
extern "C" bool fortran_tao_lat_sigma_calc_needed(
    const char* data_type /* 0D_NOT_character inout */,
    const char* data_source /* 0D_NOT_character inout */,
    bool& do_lat_sigma /* 0D_NOT_logical inout */);
void tao_lat_sigma_calc_needed(
    std::string& data_type,
    std::string& data_source,
    bool& do_lat_sigma);
extern "C" void fortran_tao_lat_sigma_track(
    void* tao_lat /* 0D_NOT_type in */,
    bool& calc_ok /* 0D_NOT_logical out */,
    int& ix_branch /* 0D_NOT_integer in */,
    bool* print_err /* 0D_NOT_logical in */,
    bool* force_calc /* 0D_NOT_logical in */);
bool tao_lat_sigma_track(
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> force_calc = std::nullopt);
extern "C" void fortran_tao_lattice_branches_equal_tao_lattice_branches(
    void* tlb1 /* 1D_ALLOC_type inout */,
    void* tlb2 /* 1D_ALLOC_type in */);
void tao_lattice_branches_equal_tao_lattice_branches(
    TaoLatticeBranchProxyAlloc1D& tlb1,
    TaoLatticeBranchProxyAlloc1D& tlb2);
extern "C" void fortran_tao_lattice_calc(
    bool& calc_ok /* 0D_NOT_logical out */,
    bool& print_err /* 0D_NOT_logical out */);
struct TaoLatticeCalc {
  bool calc_ok;
  bool print_err;
};
Tao::TaoLatticeCalc tao_lattice_calc();
extern "C" void fortran_tao_lattice_equal_tao_lattice(
    void* lat1 /* 0D_NOT_type inout */,
    void* lat2 /* 0D_NOT_type in */);
void tao_lattice_equal_tao_lattice(
    TaoLatticeProxy& lat1,
    TaoLatticeProxy& lat2);
extern "C" void fortran_tao_limit_calc(bool& limited /* 0D_NOT_logical out */);
bool tao_limit_calc();
extern "C" void fortran_tao_lm_optimizer(bool& abort /* 0D_NOT_logical out */);
bool tao_lm_optimizer();
extern "C" void fortran_tao_lmdif_optimizer(
    bool& abort /* 0D_NOT_logical out */);
bool tao_lmdif_optimizer();
extern "C" void fortran_tao_load_this_datum(
    void* vec /* 1D_ALLOC_real inout */,
    void* ele_ref /* 0D_PTR_type inout */,
    void* ele_start /* 0D_PTR_type inout */,
    void* ele /* 0D_PTR_type inout */,
    double& datum_value /* 0D_NOT_real inout */,
    bool& valid_value /* 0D_NOT_logical inout */,
    void* datum /* 0D_NOT_type inout */,
    void* branch /* 0D_NOT_type inout */,
    const char* why_invalid /* 0D_NOT_character inout */,
    void* good /* 1D_ALLOC_logical inout */);
void tao_load_this_datum(
    RealAlloc1D& vec,
    EleProxy& ele_ref,
    EleProxy& ele_start,
    EleProxy& ele,
    double& datum_value,
    bool& valid_value,
    TaoDataProxy& datum,
    BranchProxy& branch,
    optional_ref<std::string> why_invalid = std::nullopt,
    optional_ref<BoolAlloc1D> good = std::nullopt);
extern "C" void fortran_tao_locate_all_elements(
    const char* ele_list /* 0D_NOT_character in */,
    void* eles /* 1D_ALLOC_type out */,
    bool& err /* 0D_NOT_logical out */,
    bool* ignore_blank /* 0D_NOT_logical in */);
struct TaoLocateAllElements {
  ElePointerProxyAlloc1D eles;
  bool err;
};
Tao::TaoLocateAllElements tao_locate_all_elements(
    std::string ele_list,
    std::optional<bool> ignore_blank = std::nullopt);
extern "C" void fortran_tao_locate_elements(
    const char* ele_list /* 0D_NOT_character in */,
    int& ix_universe /* 0D_NOT_integer in */,
    void* eles /* 1D_ALLOC_type out */,
    bool& err /* 0D_NOT_logical out */,
    int* lat_type /* 0D_NOT_integer in */,
    bool* ignore_blank /* 0D_NOT_logical in */,
    int* err_stat_level /* 0D_NOT_integer in */,
    bool* above_ubound_is_err /* 0D_NOT_logical inout */,
    int* ix_branch /* 0D_NOT_integer in */,
    bool* multiple_eles_is_err /* 0D_NOT_logical in */);
struct TaoLocateElements {
  ElePointerProxyAlloc1D eles;
  bool err;
};
Tao::TaoLocateElements tao_locate_elements(
    std::string ele_list,
    int ix_universe,
    std::optional<int> lat_type = std::nullopt,
    std::optional<bool> ignore_blank = std::nullopt,
    std::optional<int> err_stat_level = std::nullopt,
    optional_ref<bool> above_ubound_is_err = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> multiple_eles_is_err = std::nullopt);
extern "C" void fortran_tao_mark_lattice_ele(void* lat /* 0D_NOT_type inout */);
void tao_mark_lattice_ele(LatProxy& lat);
extern "C" bool fortran_tao_merit(
    bool& calc_ok /* 0D_NOT_logical out */,
    double& this_merit /* 0D_NOT_real inout */);
bool tao_merit(double& this_merit);

// Skipped unusable routine tao_mrq_func:
// - Variable out sized array: 2D_NOT_real

// Skipped unusable routine tao_next_switch:
// - Variable-sized in character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_next_word(
    const char* line /* 0D_NOT_character inout */,
    const char* word /* 0D_NOT_character out */);
std::string tao_next_word(std::string& line);
extern "C" bool fortran_tao_one_turn_map_calc_needed(
    const char* data_type /* 0D_NOT_character inout */,
    const char* data_source /* 0D_NOT_character inout */,
    bool& do_one_turn_map /* 0D_NOT_logical inout */);
void tao_one_turn_map_calc_needed(
    std::string& data_type,
    std::string& data_source,
    bool& do_one_turn_map);
extern "C" void fortran_tao_open_file(
    const char* file /* 0D_NOT_character inout */,
    int& iunit /* 0D_NOT_integer out */,
    const char* file_name /* 0D_NOT_character in */,
    int& error_severity /* 0D_NOT_integer in */,
    bool* binary /* 0D_NOT_logical in */);
int tao_open_file(
    std::string& file,
    std::string file_name,
    int error_severity,
    std::optional<bool> binary = std::nullopt);
extern "C" bool fortran_tao_open_scratch_file(
    bool& err /* 0D_NOT_logical out */,
    int& iu /* 0D_NOT_integer inout */);
bool tao_open_scratch_file(int& iu);
extern "C" bool fortran_tao_optimization_status(
    void* datum /* 0D_NOT_type in */,
    const char* why_str /* 0D_NOT_character inout */);
void tao_optimization_status(TaoDataProxy& datum, std::string& why_str);
extern "C" void fortran_tao_orbit_beta_wave_anal(
    void* plot /* 0D_NOT_type inout */);
void tao_orbit_beta_wave_anal(TaoPlotProxy& plot);
extern "C" bool fortran_tao_oreint_building_wall_pt(
    void* pt_in /* 0D_NOT_type in */,
    void* pt_out /* 0D_NOT_type inout */);
void tao_oreint_building_wall_pt(
    TaoBuildingWallPointProxy& pt_in,
    TaoBuildingWallPointProxy& pt_out);
extern "C" bool fortran_tao_param_value_at_s(
    const char* dat_name /* 0D_NOT_character inout */,
    void* ele_to_s /* 0D_NOT_type in */,
    void* ele_here /* 0D_NOT_type in */,
    void* orbit /* 0D_NOT_type in */,
    bool& err_flag /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character out */,
    bool& print_err /* 0D_NOT_logical out */,
    bool& bad_datum /* 0D_NOT_logical out */,
    double& value /* 0D_NOT_real inout */);
struct TaoParamValueAtS {
  bool err_flag;
  std::string why_invalid;
  bool print_err;
  bool bad_datum;
};
Tao::TaoParamValueAtS tao_param_value_at_s(
    std::string& dat_name,
    EleProxy& ele_to_s,
    EleProxy& ele_here,
    CoordProxy& orbit,
    double& value);

// Skipped unusable routine tao_param_value_routine:
// - Untranslated type: tao_eval_node_struct (0D)
extern "C" void fortran_tao_parse_command_args(
    bool& error /* 0D_NOT_logical out */,
    const char* cmd_line /* 0D_NOT_character inout */);
bool tao_parse_command_args(optional_ref<std::string> cmd_line = std::nullopt);
extern "C" void fortran_tao_parse_element_param_str(
    bool& err /* 0D_NOT_logical out */,
    const char* in_str /* 0D_NOT_character in */,
    const char* uni /* 0D_NOT_character out */,
    const char* element /* 0D_NOT_character out */,
    const char* parameter /* 0D_NOT_character out */,
    int& where /* 0D_NOT_integer out */,
    const char* component /* 0D_NOT_character out */);
struct TaoParseElementParamStr {
  bool err;
  std::string uni;
  std::string element;
  std::string parameter;
  int where;
  std::string component;
};
Tao::TaoParseElementParamStr tao_parse_element_param_str(std::string in_str);
extern "C" void fortran_tao_particle_data_value(
    const char* data_type /* 0D_NOT_character in */,
    void* p /* 1D_ALLOC_type in */,
    void* value /* 1D_ALLOC_real out */,
    bool& err /* 0D_NOT_logical out */,
    void* ele /* 0D_NOT_type in */,
    int& ix_bunch /* 0D_NOT_integer in */);
struct TaoParticleDataValue {
  RealAlloc1D value;
  bool err;
};
Tao::TaoParticleDataValue tao_particle_data_value(
    std::string data_type,
    CoordProxyAlloc1D& p,
    EleProxy& ele,
    int ix_bunch);
extern "C" void fortran_tao_pause_cmd(double& time /* 0D_NOT_real in */);
void tao_pause_cmd(double time);
extern "C" bool fortran_tao_phase_space_axis_index(
    const char* data_type /* 0D_NOT_character in */,
    bool& err /* 0D_NOT_logical in */,
    int& ix_axis /* 0D_NOT_integer out */);
int tao_phase_space_axis_index(std::string data_type, bool err);
extern "C" void fortran_tao_phase_wave_anal(void* plot /* 0D_NOT_type inout */);
void tao_phase_wave_anal(TaoPlotProxy& plot);
extern "C" void fortran_tao_pick_universe(
    const char* name_in /* 0D_NOT_character in */,
    const char* name_out /* 0D_NOT_character out */,
    void* picked /* 1D_ALLOC_logical out */,
    bool& err /* 0D_NOT_logical out */,
    int& ix_uni /* 0D_NOT_integer out */,
    bool& explicit_uni /* 0D_NOT_logical out */,
    int* dflt_uni /* 0D_NOT_integer in */,
    bool* pure_uni /* 0D_NOT_logical in */);
struct TaoPickUniverse {
  std::string name_out;
  BoolAlloc1D picked;
  bool err;
  int ix_uni;
  bool explicit_uni;
};
Tao::TaoPickUniverse tao_pick_universe(
    std::string name_in,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<bool> pure_uni = std::nullopt);
extern "C" void fortran_tao_pipe_cmd(
    const char* input_str /* 0D_NOT_character in */);
void tao_pipe_cmd(std::string input_str);
extern "C" void fortran_tao_place_cmd(
    const char* where /* 0D_NOT_character in */,
    const char* who /* 0D_NOT_character in */,
    bool* no_buffer /* 0D_NOT_logical in */);
void tao_place_cmd(
    std::string where,
    std::string who,
    std::optional<bool> no_buffer = std::nullopt);
extern "C" void fortran_tao_plot_cmd(
    const char* where /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */);
void tao_plot_cmd(std::string where, std::string component);
extern "C" void fortran_tao_plot_data(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_plot_data(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_plot_histogram(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_plot_histogram(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_plot_key_table(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_plot_key_table(TaoPlotProxy& plot, TaoGraphProxy& graph);
extern "C" void fortran_tao_plot_setup();
void tao_plot_setup();
extern "C" void fortran_tao_plot_struct_transfer(
    void* plot_in /* 0D_NOT_type in */,
    void* plot_out /* 0D_NOT_type out */);
TaoPlotProxy tao_plot_struct_transfer(TaoPlotProxy& plot_in);
extern "C" void fortran_tao_plot_wave(
    void* plot /* 0D_NOT_type in */,
    void* graph /* 0D_NOT_type in */);
void tao_plot_wave(TaoPlotProxy& plot, TaoGraphProxy& graph);

// Skipped unusable routine tao_point_d1_to_data:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_point_v1_to_var:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_pointer_to_branches:
// - Untranslated type: branch_pointer_struct (1D)
extern "C" bool fortran_tao_pointer_to_building_wall_shape(
    const char* wall_name /* 0D_NOT_character in */,
    void* e_shape /* 0D_PTR_type inout */);
void tao_pointer_to_building_wall_shape(
    std::string wall_name,
    TaoEleShapeProxy& e_shape);
extern "C" bool fortran_tao_pointer_to_datum(
    void* d1 /* 0D_NOT_type in */,
    const char* ele_name /* 0D_NOT_character in */,
    void* datum_ptr /* 0D_PTR_type inout */);
void tao_pointer_to_datum(
    TaoD1DataProxy& d1,
    std::string ele_name,
    TaoDataProxy& datum_ptr);
extern "C" bool fortran_tao_pointer_to_datum_ele(
    void* lat /* 0D_NOT_type in */,
    const char* ele_name /* 0D_NOT_character inout */,
    int& ix_ele /* 0D_NOT_integer in */,
    void* datum /* 0D_NOT_type in */,
    bool& valid /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character out */,
    bool* print_err /* 0D_NOT_logical in */,
    void* ele /* 0D_PTR_type out */);
struct TaoPointerToDatumEle {
  bool valid;
  std::string why_invalid;
  EleProxy ele;
};
Tao::TaoPointerToDatumEle tao_pointer_to_datum_ele(
    LatProxy& lat,
    std::string& ele_name,
    int ix_ele,
    TaoDataProxy& datum,
    std::optional<bool> print_err = std::nullopt);
extern "C" bool fortran_tao_pointer_to_ele_shape(
    int& ix_uni /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type in */,
    void* ele_shape /* 1D_ALLOC_type in */,
    const char* dat_var_name /* 0D_NOT_character out */,
    double& dat_var_value /* 0D_NOT_real out */,
    int* ix_shape_min /* 0D_NOT_integer inout */,
    void* e_shape /* 0D_PTR_type inout */);
struct TaoPointerToEleShape {
  std::string dat_var_name;
  double dat_var_value;
};
Tao::TaoPointerToEleShape tao_pointer_to_ele_shape(
    int ix_uni,
    EleProxy& ele,
    TaoEleShapeProxyAlloc1D& ele_shape,
    optional_ref<int> ix_shape_min,
    TaoEleShapeProxy& e_shape);
extern "C" bool fortran_tao_pointer_to_tao_lat(
    void* u /* 0D_NOT_type in */,
    int* lat_type /* 0D_NOT_integer in */,
    void* tao_lat /* 0D_PTR_type inout */);
void tao_pointer_to_tao_lat(
    TaoUniverseProxy& u,
    std::optional<int> lat_type,
    TaoLatticeProxy& tao_lat);
extern "C" void fortran_tao_pointer_to_universes(
    const char* name_in /* 0D_NOT_character in */,
    void* unis /* 1D_ALLOC_type out */,
    bool& err /* 0D_NOT_logical out */,
    const char* name_out /* 0D_NOT_character out */,
    bool& explicit_uni /* 0D_NOT_logical out */,
    int* dflt_uni /* 0D_NOT_integer in */);
struct TaoPointerToUniverses {
  TaoUniversePointerProxyAlloc1D unis;
  bool err;
  std::string name_out;
  bool explicit_uni;
};
Tao::TaoPointerToUniverses tao_pointer_to_universes(
    std::string name_in,
    std::optional<int> dflt_uni = std::nullopt);
extern "C" void fortran_tao_pointer_to_var_in_lattice(
    void* var /* 0D_NOT_type in */,
    int& ix_uni /* 0D_NOT_integer in */,
    void* ele /* 0D_NOT_type inout */,
    bool& err /* 0D_NOT_logical out */);
bool tao_pointer_to_var_in_lattice(TaoVarProxy& var, int ix_uni, EleProxy& ele);
extern "C" void fortran_tao_pointer_to_var_in_lattice2(
    void* var /* 0D_NOT_type in */,
    int& ix_uni /* 0D_NOT_integer in */,
    bool& err /* 0D_NOT_logical out */);
bool tao_pointer_to_var_in_lattice2(TaoVarProxy& var, int ix_uni);
extern "C" void fortran_tao_print_command_line_info();
void tao_print_command_line_info();

// Skipped unusable routine tao_print_vars:
// - Untranslated type: tao_var_array_struct (1D)

// Skipped unusable routine tao_ptc_cmd:
// - Module name unset
extern "C" void fortran_tao_ptc_normal_form(
    bool& do_calc /* 0D_NOT_logical in */,
    void* tao_lat /* 0D_NOT_type in */,
    int& ix_branch /* 0D_NOT_integer in */,
    int* rf_on /* 0D_NOT_integer in */);
void tao_ptc_normal_form(
    bool do_calc,
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<int> rf_on = std::nullopt);
extern "C" void fortran_tao_python_cmd(
    const char* input_str /* 0D_NOT_character in */);
void tao_python_cmd(std::string input_str);
extern "C" void fortran_tao_quiet_set(
    const char* set /* 0D_NOT_character in */);
void tao_quiet_set(std::string set);
extern "C" bool fortran_tao_rad_int_calc_needed(
    const char* data_type /* 0D_NOT_character inout */,
    const char* data_source /* 0D_NOT_character inout */,
    bool& do_rad_int /* 0D_NOT_logical inout */);
void tao_rad_int_calc_needed(
    std::string& data_type,
    std::string& data_source,
    bool& do_rad_int);

// Skipped unusable routine tao_re_allocate_expression_info:
// - Untranslated type: tao_expression_info_struct (1D)

// Skipped unusable routine tao_re_associate_node_array:
// - Untranslated type: tao_eval_node_struct (0D)
extern "C" void fortran_tao_re_execute(
    const char* string /* 0D_NOT_character inout */,
    bool& err /* 0D_NOT_logical inout */);
void tao_re_execute(std::string& string, bool& err);
extern "C" void fortran_tao_read_cmd(
    const char* which /* 0D_NOT_character inout */,
    const char* unis /* 0D_NOT_character in */,
    const char* file /* 0D_NOT_character inout */,
    bool& silent /* 0D_NOT_logical in */);
void tao_read_cmd(
    std::string& which,
    std::string unis,
    std::string& file,
    bool silent);

// Skipped unusable routine tao_read_in_patterns:
// - Module name unset
extern "C" bool fortran_tao_read_phase_space_index(
    const char* name /* 0D_NOT_character in */,
    int& ixc /* 0D_NOT_integer in */,
    bool* print_err /* 0D_NOT_logical in */,
    int& ix_ps /* 0D_NOT_integer inout */);
void tao_read_phase_space_index(
    std::string name,
    int ixc,
    std::optional<bool> print_err,
    int& ix_ps);
extern "C" void fortran_tao_regression_test();
void tao_regression_test();
extern "C" void fortran_tao_remove_blank_characters(
    const char* str /* 0D_NOT_character inout */);
void tao_remove_blank_characters(std::string& str);
extern "C" void fortran_tao_run_cmd(
    const char* which /* 0D_NOT_character in */,
    bool& abort /* 0D_NOT_logical out */);
bool tao_run_cmd(std::string which);
extern "C" void fortran_tao_scale_cmd(
    const char* where /* 0D_NOT_character in */,
    double& y_min_in /* 0D_NOT_real in */,
    double& y_max_in /* 0D_NOT_real in */,
    const char* axis /* 0D_NOT_character in */,
    bool* include_wall /* 0D_NOT_logical in */,
    const char* gang /* 0D_NOT_character in */,
    bool* exact /* 0D_NOT_logical in */,
    bool* turn_autoscale_off /* 0D_NOT_logical in */);
void tao_scale_cmd(
    std::string where,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> exact = std::nullopt,
    std::optional<bool> turn_autoscale_off = std::nullopt);
extern "C" void fortran_tao_scale_graph(
    void* graph /* 0D_NOT_type inout */,
    double& y_min /* 0D_NOT_real in */,
    double& y_max /* 0D_NOT_real in */,
    const char* axis /* 0D_NOT_character in */,
    bool* include_wall /* 0D_NOT_logical in */,
    double* y_range /* 1D_NOT_real out */,
    double* y2_range /* 1D_NOT_real out */);
struct TaoScaleGraph {
  FixedArray1D<Real, 2> y_range;
  FixedArray1D<Real, 2> y2_range;
};
Tao::TaoScaleGraph tao_scale_graph(
    TaoGraphProxy& graph,
    double y_min,
    double y_max,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt);
extern "C" void fortran_tao_scale_ping_data(void* u /* 0D_NOT_type inout */);
void tao_scale_ping_data(TaoUniverseProxy& u);
extern "C" void fortran_tao_scale_plot(
    void* plot /* 0D_NOT_type inout */,
    double& y_min_in /* 0D_NOT_real in */,
    double& y_max_in /* 0D_NOT_real in */,
    const char* axis /* 0D_NOT_character in */,
    bool* include_wall /* 0D_NOT_logical in */,
    const char* gang /* 0D_NOT_character in */,
    bool* skip_lat_layout /* 0D_NOT_logical in */);
void tao_scale_plot(
    TaoPlotProxy& plot,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> skip_lat_layout = std::nullopt);
extern "C" void fortran_tao_scratch_values_calc(
    void* ele_ref /* 0D_PTR_type inout */,
    void* ele_start /* 0D_PTR_type inout */,
    void* ele /* 0D_PTR_type inout */,
    void* datum /* 0D_NOT_type inout */,
    void* branch /* 0D_NOT_type inout */,
    void* orbit /* 1D_ALLOC_type inout */);
void tao_scratch_values_calc(
    EleProxy& ele_ref,
    EleProxy& ele_start,
    EleProxy& ele,
    TaoDataProxy& datum,
    BranchProxy& branch,
    CoordProxyAlloc1D& orbit);
extern "C" void fortran_tao_set_beam_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    const char* branch_str /* 0D_NOT_character in */);
void tao_set_beam_cmd(
    std::string who,
    std::string value_str,
    std::string branch_str);
extern "C" void fortran_tao_set_beam_init_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    const char* branch_str /* 0D_NOT_character in */);
void tao_set_beam_init_cmd(
    std::string who,
    std::string value_str,
    std::string branch_str);
extern "C" void fortran_tao_set_bmad_com_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_bmad_com_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_branch_cmd(
    const char* branch_str /* 0D_NOT_character in */,
    const char* component_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_branch_cmd(
    std::string branch_str,
    std::string component_str,
    std::string value_str);
extern "C" void fortran_tao_set_calculate_cmd(
    const char* switch_ /* 0D_NOT_character inout */);
void tao_set_calculate_cmd(optional_ref<std::string> switch_ = std::nullopt);
extern "C" void fortran_tao_set_curve_cmd(
    const char* curve_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_curve_cmd(
    std::string curve_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_curve_invalid(
    void* curve /* 0D_NOT_type inout */,
    const char* why_invalid /* 0D_NOT_character in */,
    bool* print_err /* 0D_NOT_logical in */);
void tao_set_curve_invalid(
    TaoCurveProxy& curve,
    std::string why_invalid,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_set_data_cmd(
    const char* who_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    bool* silent /* 0D_NOT_logical inout */);
void tao_set_data_cmd(
    std::string who_str,
    std::string value_str,
    optional_ref<bool> silent = std::nullopt);
extern "C" void fortran_tao_set_data_useit_opt(
    void* data /* 1D_ALLOC_type in */);
void tao_set_data_useit_opt(
    optional_ref<TaoDataProxyAlloc1D> data = std::nullopt);
extern "C" void fortran_tao_set_default_cmd(
    const char* who_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_default_cmd(std::string who_str, std::string value_str);
extern "C" void fortran_tao_set_drawing_cmd(
    void* drawing /* 0D_NOT_type in */,
    const char* component /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_drawing_cmd(
    TaoDrawingProxy& drawing,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_dynamic_aperture_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_dynamic_aperture_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_elements_cmd(
    const char* ele_list /* 0D_NOT_character in */,
    const char* attribute /* 0D_NOT_character in */,
    const char* value /* 0D_NOT_character in */,
    bool& update /* 0D_NOT_logical inout */);
void tao_set_elements_cmd(
    std::string ele_list,
    std::string attribute,
    std::string value,
    bool& update);

// Skipped unusable routine tao_set_flags_for_changed_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" void fortran_tao_set_floor_plan_axis_label(
    void* graph /* 0D_NOT_type inout */,
    void* axis_in /* 0D_NOT_type inout */,
    void* axis_out /* 0D_NOT_type inout */,
    const char* which /* 0D_NOT_character inout */);
void tao_set_floor_plan_axis_label(
    TaoGraphProxy& graph,
    QpAxisProxy& axis_in,
    QpAxisProxy& axis_out,
    std::string& which);
extern "C" void fortran_tao_set_geodesic_lm_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_geodesic_lm_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_global_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_global_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_graph_cmd(
    const char* graph_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_graph_cmd(
    std::string graph_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_integer_value(
    int& var /* 0D_NOT_integer out */,
    const char* var_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */,
    int* min_val /* 0D_NOT_integer in */,
    int* max_val /* 0D_NOT_integer in */,
    bool* print_err /* 0D_NOT_logical in */);
struct TaoSetIntegerValue {
  int var;
  bool error;
};
Tao::TaoSetIntegerValue tao_set_integer_value(
    std::string var_str,
    std::string value_str,
    std::optional<int> min_val = std::nullopt,
    std::optional<int> max_val = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_set_invalid(
    void* datum /* 0D_NOT_type in */,
    const char* message /* 0D_NOT_character in */,
    const char* why_invalid /* 0D_NOT_character out */,
    bool* exterminate /* 0D_NOT_logical in */,
    int* err_level /* 0D_NOT_integer in */,
    bool* print_err /* 0D_NOT_logical in */);
std::string tao_set_invalid(
    TaoDataProxy& datum,
    std::string message,
    std::optional<bool> exterminate = std::nullopt,
    std::optional<int> err_level = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_set_key_cmd(
    const char* key_str /* 0D_NOT_character in */,
    const char* cmd_str /* 0D_NOT_character in */);
void tao_set_key_cmd(std::string key_str, std::string cmd_str);
extern "C" void fortran_tao_set_lattice_cmd(
    const char* dest_lat /* 0D_NOT_character in */,
    const char* source_lat /* 0D_NOT_character in */);
void tao_set_lattice_cmd(std::string dest_lat, std::string source_lat);
extern "C" void fortran_tao_set_logical_value(
    bool& var /* 0D_NOT_logical out */,
    const char* var_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */);
struct TaoSetLogicalValue {
  bool var;
  bool error;
};
Tao::TaoSetLogicalValue tao_set_logical_value(
    std::string var_str,
    std::string value_str);
extern "C" void fortran_tao_set_openmp_n_threads(
    int& n_threads /* 0D_NOT_integer in */);
void tao_set_openmp_n_threads(int n_threads);
extern "C" void fortran_tao_set_opt_vars(
    void* var_vec /* 1D_ALLOC_real in */,
    bool* print_limit_warning /* 0D_NOT_logical in */);
void tao_set_opt_vars(
    RealAlloc1D& var_vec,
    std::optional<bool> print_limit_warning = std::nullopt);
extern "C" void fortran_tao_set_opti_de_param_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_opti_de_param_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_particle_start_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_particle_start_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_plot_cmd(
    const char* plot_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_plot_cmd(
    std::string plot_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_plot_page_cmd(
    const char* component /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    const char* value_str2 /* 0D_NOT_character in */);
void tao_set_plot_page_cmd(
    std::string component,
    std::string value_str,
    std::optional<std::string> value_str2 = std::nullopt);

// Skipped unusable routine tao_set_plotting:
// - Untranslated type: tao_plot_page_input (0D)
extern "C" void fortran_tao_set_ptc_com_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_ptc_com_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_qp_axis_struct(
    const char* qp_axis_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    void* qp_axis /* 0D_NOT_type inout */,
    const char* value /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */,
    int& ix_uni /* 0D_NOT_integer out */);
struct TaoSetQpAxisStruct {
  bool error;
  int ix_uni;
};
Tao::TaoSetQpAxisStruct tao_set_qp_axis_struct(
    std::string qp_axis_name,
    std::string component,
    QpAxisProxy& qp_axis,
    std::string value);
extern "C" void fortran_tao_set_qp_point_struct(
    const char* qp_point_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    void* qp_point /* 0D_NOT_type inout */,
    const char* value /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */,
    int& ix_uni /* 0D_NOT_integer out */);
struct TaoSetQpPointStruct {
  bool error;
  int ix_uni;
};
Tao::TaoSetQpPointStruct tao_set_qp_point_struct(
    std::string qp_point_name,
    std::string component,
    QpPointProxy& qp_point,
    std::string value);
extern "C" void fortran_tao_set_qp_rect_struct(
    const char* qp_rect_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    void* qp_rect /* 0D_NOT_type inout */,
    const char* value /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */,
    int& ix_uni /* 0D_NOT_integer out */);
struct TaoSetQpRectStruct {
  bool error;
  int ix_uni;
};
Tao::TaoSetQpRectStruct tao_set_qp_rect_struct(
    std::string qp_rect_name,
    std::string component,
    QpRectProxy& qp_rect,
    std::string value);
extern "C" void fortran_tao_set_ran_state_cmd(
    const char* state_string /* 0D_NOT_character in */);
void tao_set_ran_state_cmd(std::string state_string);
extern "C" void fortran_tao_set_real_value(
    double& var /* 0D_NOT_real out */,
    const char* var_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    bool& error /* 0D_NOT_logical out */,
    double* min_val /* 0D_NOT_real in */,
    double* max_val /* 0D_NOT_real in */,
    int* dflt_uni /* 0D_NOT_integer in */);
struct TaoSetRealValue {
  double var;
  bool error;
};
Tao::TaoSetRealValue tao_set_real_value(
    std::string var_str,
    std::string value_str,
    std::optional<double> min_val = std::nullopt,
    std::optional<double> max_val = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt);
extern "C" void fortran_tao_set_region_cmd(
    const char* region_name /* 0D_NOT_character in */,
    const char* component /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_region_cmd(
    std::string region_name,
    std::string component,
    std::string value_str);
extern "C" void fortran_tao_set_space_charge_com_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_space_charge_com_cmd(std::string who, std::string value_str);

// Skipped unusable routine tao_set_switch_value:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_set_symbolic_number_cmd(
    const char* sym_str /* 0D_NOT_character in */,
    const char* num_str /* 0D_NOT_character in */,
    double* val /* 0D_NOT_real in */);
void tao_set_symbolic_number_cmd(
    std::string sym_str,
    std::optional<std::string> num_str = std::nullopt,
    std::optional<double> val = std::nullopt);
extern "C" void fortran_tao_set_tune_cmd(
    const char* branch_str /* 0D_NOT_character in */,
    const char* mask_str /* 0D_NOT_character in */,
    bool& print_list /* 0D_NOT_logical in */,
    const char* qa_str /* 0D_NOT_character in */,
    const char* qb_str /* 0D_NOT_character in */,
    bool& delta_input /* 0D_NOT_logical in */);
void tao_set_tune_cmd(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string qa_str,
    std::string qb_str,
    bool delta_input);
extern "C" void fortran_tao_set_universe_cmd(
    const char* uni /* 0D_NOT_character in */,
    const char* who /* 0D_NOT_character in */,
    const char* what /* 0D_NOT_character in */);
void tao_set_universe_cmd(std::string uni, std::string who, std::string what);
extern "C" void fortran_tao_set_var_cmd(
    const char* var_str /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */);
void tao_set_var_cmd(std::string var_str, std::string value_str);
extern "C" void fortran_tao_set_var_model_value(
    void* var /* 0D_NOT_type in */,
    double& value /* 0D_NOT_real in */,
    bool* print_limit_warning /* 0D_NOT_logical in */);
void tao_set_var_model_value(
    TaoVarProxy& var,
    double value,
    std::optional<bool> print_limit_warning = std::nullopt);
extern "C" void fortran_tao_set_var_useit_opt();
void tao_set_var_useit_opt();
extern "C" void fortran_tao_set_wave_cmd(
    const char* who /* 0D_NOT_character in */,
    const char* value_str /* 0D_NOT_character in */,
    bool& err /* 0D_NOT_logical out */);
bool tao_set_wave_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_z_tune_cmd(
    const char* branch_str /* 0D_NOT_character in */,
    const char* q_str /* 0D_NOT_character in */,
    bool& delta_input /* 0D_NOT_logical in */);
void tao_set_z_tune_cmd(
    std::string branch_str,
    std::string q_str,
    bool delta_input);
extern "C" void fortran_tao_setup_key_table();
void tao_setup_key_table();
extern "C" void fortran_tao_shape_init(
    void* shape /* 0D_NOT_type inout */,
    bool& err /* 0D_NOT_logical out */,
    bool* print_err /* 0D_NOT_logical in */);
bool tao_shape_init(
    TaoEleShapeProxy& shape,
    std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_show_cmd(
    const char* what /* 0D_NOT_character in */);
void tao_show_cmd(std::string what);
extern "C" void fortran_tao_show_constraints(
    int& iunit /* 0D_NOT_integer in */,
    const char* form /* 0D_NOT_character in */);
void tao_show_constraints(int iunit, std::string form);

// Skipped unusable routine tao_show_this:
// - Variable-sized out character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_single_mode(
    const char* char_ /* 0D_NOT_character in */);
void tao_single_mode(std::string char_);
extern "C" void fortran_tao_single_track(
    void* tao_lat /* 0D_NOT_type in */,
    bool& calc_ok /* 0D_NOT_logical out */,
    int& ix_branch /* 0D_NOT_integer in */,
    bool* print_err /* 0D_NOT_logical in */);
bool tao_single_track(
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<bool> print_err = std::nullopt);
extern "C" bool fortran_tao_spin_matrices_calc_needed(
    const char* data_type /* 0D_NOT_character inout */,
    const char* data_source /* 0D_NOT_character inout */,
    bool& do_calc /* 0D_NOT_logical inout */);
void tao_spin_matrices_calc_needed(
    std::string& data_type,
    std::string& data_source,
    bool& do_calc);

// Skipped unusable routine tao_spin_matrix_calc:
// - Routine in configuration skip list

// Skipped unusable routine tao_spin_polarization_calc:
// - Routine in configuration skip list
extern "C" void fortran_tao_spin_tracking_turn_on();
void tao_spin_tracking_turn_on();
extern "C" void fortran_tao_split_component(
    const char* comp_str /* 0D_NOT_character in */,
    void* comp /* 1D_ALLOC_type out */,
    bool& err /* 0D_NOT_logical out */);
struct TaoSplitComponent {
  TaoDataVarComponentProxyAlloc1D comp;
  bool err;
};
Tao::TaoSplitComponent tao_split_component(std::string comp_str);
extern "C" bool fortran_tao_srdt_calc_needed(
    const char* data_type /* 0D_NOT_character inout */,
    const char* data_source /* 0D_NOT_character inout */,
    int& do_srdt /* 0D_NOT_integer inout */);
void tao_srdt_calc_needed(
    std::string& data_type,
    std::string& data_source,
    int& do_srdt);
extern "C" bool fortran_tao_subin_uni_number(
    const char* name_in /* 0D_NOT_character in */,
    int& ix_uni /* 0D_NOT_integer in */,
    const char* name_out /* 0D_NOT_character out */,
    bool& ok /* 0D_NOT_logical inout */);
std::string tao_subin_uni_number(std::string name_in, int ix_uni, bool& ok);

// Skipped unusable routine tao_svd_func:
// - Variable out sized array: 2D_NOT_real
extern "C" void fortran_tao_svd_optimizer(bool& abort /* 0D_NOT_logical out */);
bool tao_svd_optimizer();
extern "C" void fortran_tao_symbol_import_from_lat(
    void* lat /* 0D_NOT_type inout */);
void tao_symbol_import_from_lat(LatProxy& lat);
extern "C" void fortran_tao_taper_cmd(
    const char* except /* 0D_NOT_character in */,
    const char* uni_names /* 0D_NOT_character in */);
void tao_taper_cmd(std::string except, std::string uni_names);

// Skipped unusable routine tao_timer:
// - Module name unset
extern "C" void fortran_tao_to_change_number(
    const char* num_str /* 0D_NOT_character inout */,
    int& n_size /* 0D_NOT_integer inout */,
    void* change_number /* 1D_ALLOC_real inout */,
    const char* abs_or_rel /* 0D_NOT_character inout */,
    bool& err /* 0D_NOT_logical inout */);
void tao_to_change_number(
    std::string& num_str,
    int& n_size,
    RealAlloc1D& change_number,
    std::string& abs_or_rel,
    bool& err);
extern "C" void fortran_tao_to_int(
    const char* str /* 0D_NOT_character inout */,
    int& i_int /* 0D_NOT_integer inout */,
    bool& err /* 0D_NOT_logical inout */);
void tao_to_int(std::string& str, int& i_int, bool& err);
extern "C" void fortran_tao_to_phase_and_coupling_reading(
    void* ele /* 0D_NOT_type in */,
    void* bpm_data /* 0D_NOT_type out */,
    bool& valid_value /* 0D_NOT_logical out */,
    const char* why_invalid /* 0D_NOT_character inout */,
    void* datum /* 0D_NOT_type inout */);
struct TaoToPhaseAndCouplingReading {
  BpmPhaseCouplingProxy bpm_data;
  bool valid_value;
};
Tao::TaoToPhaseAndCouplingReading tao_to_phase_and_coupling_reading(
    EleProxy& ele,
    std::string& why_invalid,
    TaoDataProxy& datum);
extern "C" void fortran_tao_to_real(
    const char* expression /* 0D_NOT_character in */,
    double& value /* 0D_NOT_real out */,
    bool& err_flag /* 0D_NOT_logical out */);
struct TaoToReal {
  double value;
  bool err_flag;
};
Tao::TaoToReal tao_to_real(std::string expression);

// Skipped unusable routine tao_to_top10:
// - Untranslated type: tao_top10_struct (1D)
extern "C" bool fortran_tao_too_many_particles_lost(
    void* beam /* 0D_NOT_type inout */,
    bool& no_beam /* 0D_NOT_logical inout */);
void tao_too_many_particles_lost(BeamProxy& beam, bool& no_beam);
extern "C" void fortran_tao_top10_derivative_print();
void tao_top10_derivative_print();
extern "C" void fortran_tao_top10_merit_categories_print(
    int& iunit /* 0D_NOT_integer in */);
void tao_top10_merit_categories_print(int iunit);
extern "C" void fortran_tao_top_level(
    const char* command /* 0D_NOT_character in */,
    int& errcode /* 0D_NOT_integer out */);
int tao_top_level(std::optional<std::string> command = std::nullopt);
extern "C" bool fortran_tao_tracking_ele_index(
    void* ele /* 0D_PTR_type in */,
    void* datum /* 0D_NOT_type in */,
    int& ix_branch /* 0D_NOT_integer out */,
    int& ix_ele /* 0D_NOT_integer out */);
struct TaoTrackingEleIndex {
  int ix_branch;
  int ix_ele;
};
Tao::TaoTrackingEleIndex tao_tracking_ele_index(
    EleProxy& ele,
    TaoDataProxy& datum);
extern "C" void fortran_tao_turn_on_special_calcs_if_needed_for_plotting();
void tao_turn_on_special_calcs_if_needed_for_plotting();

// Skipped unusable routine tao_type_expression_tree:
// - Untranslated type: tao_eval_node_struct (0D)
extern "C" bool fortran_tao_uni_atsign_index(
    const char* string /* 0D_NOT_character in */,
    int& ix_amp /* 0D_NOT_integer out */);
int tao_uni_atsign_index(std::string string);
extern "C" bool fortran_tao_universe_index(
    int& i_uni /* 0D_NOT_integer in */,
    bool* neg2_to_default /* 0D_NOT_logical in */,
    int& i_this_uni /* 0D_NOT_integer inout */);
void tao_universe_index(
    int i_uni,
    std::optional<bool> neg2_to_default,
    int& i_this_uni);
extern "C" void fortran_tao_use_data(
    const char* action /* 0D_NOT_character in */,
    const char* data_name /* 0D_NOT_character in */);
void tao_use_data(std::string action, std::string data_name);
extern "C" void fortran_tao_use_var(
    const char* action /* 0D_NOT_character in */,
    const char* var_name /* 0D_NOT_character in */);
void tao_use_var(std::string action, std::string var_name);
extern "C" bool fortran_tao_user_is_terminating_optimization(
    bool& is_terminating /* 0D_NOT_logical out */);
bool tao_user_is_terminating_optimization();
extern "C" bool fortran_tao_var1_name(
    void* var /* 0D_NOT_type in */,
    const char* var1_name /* 0D_NOT_character inout */);
void tao_var1_name(TaoVarProxy& var, std::string& var1_name);
extern "C" bool fortran_tao_var_attrib_name(
    void* var /* 0D_NOT_type in */,
    const char* var_attrib_name /* 0D_NOT_character inout */);
void tao_var_attrib_name(TaoVarProxy& var, std::string& var_attrib_name);
extern "C" void fortran_tao_var_check(
    void* eles /* 1D_ALLOC_type in */,
    const char* attribute /* 0D_NOT_character in */,
    bool& silent /* 0D_NOT_logical in */);
void tao_var_check(
    ElePointerProxyAlloc1D& eles,
    std::string attribute,
    bool silent);
extern "C" void fortran_tao_var_repoint();
void tao_var_repoint();

// Skipped unusable routine tao_var_show_use:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_var_stuffit1:
// - Untranslated type: tao_var_input (1D)
// - Untranslated type: tao_v1_var_input (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_var_stuffit2:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_var_target_calc();
void tao_var_target_calc();
extern "C" void fortran_tao_var_useit_plot_calc(
    void* graph /* 0D_NOT_type inout */,
    void* var /* 1D_ALLOC_type out */);
TaoVarProxyAlloc1D tao_var_useit_plot_calc(TaoGraphProxy& graph);
extern "C" void fortran_tao_var_write(
    const char* out_file /* 0D_NOT_character in */,
    bool* show_good_opt_only /* 0D_NOT_logical in */,
    bool* tao_format /* 0D_NOT_logical in */);
void tao_var_write(
    std::string out_file,
    std::optional<bool> show_good_opt_only = std::nullopt,
    std::optional<bool> tao_format = std::nullopt);
extern "C" void fortran_tao_veto_vars_with_zero_dmodel();
void tao_veto_vars_with_zero_dmodel();
extern "C" void fortran_tao_wave_analysis(void* plot /* 0D_NOT_type inout */);
void tao_wave_analysis(TaoPlotProxy& plot);
extern "C" void fortran_tao_wave_cmd(
    const char* curve_name /* 0D_NOT_character in */,
    const char* plot_place /* 0D_NOT_character in */,
    bool& err_flag /* 0D_NOT_logical inout */);
void tao_wave_cmd(
    std::string curve_name,
    std::string plot_place,
    bool& err_flag);
extern "C" void fortran_tao_wave_fit(
    void* curve /* 0D_NOT_type in */,
    int& ix1 /* 0D_NOT_integer in */,
    int& n_dat /* 0D_NOT_integer in */,
    void* coef /* 1D_ALLOC_real out */,
    void* rms /* 1D_ALLOC_real out */,
    void* f1 /* 1D_ALLOC_real in */,
    void* f2 /* 1D_ALLOC_real in */,
    void* f3 /* 1D_ALLOC_real in */,
    void* f4 /* 1D_ALLOC_real in */);
struct TaoWaveFit {
  RealAlloc1D coef;
  RealAlloc1D rms;
};
Tao::TaoWaveFit tao_wave_fit(
    TaoCurveProxy& curve,
    int ix1,
    int n_dat,
    RealAlloc1D& f1,
    optional_ref<RealAlloc1D> f2 = std::nullopt,
    optional_ref<RealAlloc1D> f3 = std::nullopt,
    optional_ref<RealAlloc1D> f4 = std::nullopt);
extern "C" void fortran_tao_write_cmd(
    const char* what /* 0D_NOT_character in */);
void tao_write_cmd(std::string what);

// Skipped unusable routine tao_write_lines:
// - Variable-sized in character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_x_axis_cmd(
    const char* where /* 0D_NOT_character in */,
    const char* what /* 0D_NOT_character in */);
void tao_x_axis_cmd(std::string where, std::string what);
extern "C" void fortran_tao_x_scale_cmd(
    const char* where /* 0D_NOT_character in */,
    double& x_min_in /* 0D_NOT_real in */,
    double& x_max_in /* 0D_NOT_real in */,
    bool& err /* 0D_NOT_logical out */,
    bool* include_wall /* 0D_NOT_logical in */,
    const char* gang /* 0D_NOT_character in */,
    bool* exact /* 0D_NOT_logical in */,
    bool* turn_autoscale_off /* 0D_NOT_logical in */);
bool tao_x_scale_cmd(
    std::string where,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> exact = std::nullopt,
    std::optional<bool> turn_autoscale_off = std::nullopt);
extern "C" void fortran_tao_x_scale_graph(
    void* graph /* 0D_NOT_type inout */,
    double& x_min /* 0D_NOT_real inout */,
    double& x_max /* 0D_NOT_real inout */,
    bool* include_wall /* 0D_NOT_logical inout */,
    bool* have_scaled /* 0D_NOT_logical inout */);
void tao_x_scale_graph(
    TaoGraphProxy& graph,
    double& x_min,
    double& x_max,
    optional_ref<bool> include_wall = std::nullopt,
    optional_ref<bool> have_scaled = std::nullopt);
extern "C" void fortran_tao_x_scale_plot(
    void* plot /* 0D_NOT_type in */,
    double& x_min_in /* 0D_NOT_real in */,
    double& x_max_in /* 0D_NOT_real in */,
    bool* include_wall /* 0D_NOT_logical in */,
    const char* gang /* 0D_NOT_character in */,
    bool& have_scaled /* 0D_NOT_logical out */);
bool tao_x_scale_plot(
    TaoPlotProxy& plot,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt);

// Skipped unusable routine user_signal:
// - Module name unset
} // namespace Tao
