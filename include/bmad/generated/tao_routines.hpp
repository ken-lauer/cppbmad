#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

namespace Tao {

// Skipped unusable routine avv:
// - Routine in configuration skip list

// Skipped unusable routine callback:
// - Array bounds handling: "Enum 'N' found in bounds 'n' but not in provided map."
// - Array bounds handling: "Enum 'M' found in bounds 'm' but not in provided map."
// - Array bounds handling: "Enum 'M' found in bounds 'm' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_integrate_max(
    int &ix_start /* 0D_NOT_integer in */,
    int &ix_ele /* 0D_NOT_integer in */,
    double &datum_value /* 0D_NOT_real in */,
    int &ix_m /* 0D_NOT_integer in */,
    void *branch /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &vec /* 1D_NOT_real inout */,
    void *datum /* 0D_NOT_type inout */
);
void integrate_max(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchStruct &branch,
    FArray1D<Real> &vec,
    TaoDataStruct &datum
);
extern "C" void fortran_integrate_min(
    int &ix_start /* 0D_NOT_integer in */,
    int &ix_ele /* 0D_NOT_integer in */,
    double &datum_value /* 0D_NOT_real in */,
    int &ix_m /* 0D_NOT_integer in */,
    void *branch /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &vec /* 1D_NOT_real inout */,
    void *datum /* 0D_NOT_type inout */
);
void integrate_min(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchStruct &branch,
    FArray1D<Real> &vec,
    TaoDataStruct &datum
);

// Skipped unusable routine jacobian:
// - Array bounds handling: "Enum 'NN' found in bounds 'nn' but not in provided map."
// - Array bounds handling: "Enum 'MM' found in bounds 'mm' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_abort_command_file(bool *force_abort /* 0D_NOT_logical in */);
void tao_abort_command_file(std::optional<bool> force_abort = std::nullopt);
extern "C" void fortran_tao_add_to_normal_mode_h_array(
    const char *h_str /* 0D_NOT_character in */,
    void *h_array /* 1D_ALLOC_type inout */
);
void tao_add_to_normal_mode_h_array(std::string h_str, ResonanceHStructAlloc1D h_array);
extern "C" void fortran_tao_alias_cmd(
    const char *alias /* 0D_NOT_character in */,
    const char *string /* 0D_NOT_character in */
);
void tao_alias_cmd(std::string alias, std::string string);
extern "C" void fortran_tao_allocate_data_array(
    void *u /* 0D_NOT_type inout */,
    int &n_data /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void tao_allocate_data_array(
    TaoUniverseStruct &u,
    int n_data,
    std::optional<bool> exact = std::nullopt
);
extern "C" void fortran_tao_allocate_v1_var(
    int &n_v1 /* 0D_NOT_integer in */,
    bool &save_old /* 0D_NOT_logical in */
);
void tao_allocate_v1_var(int n_v1, bool save_old);
extern "C" void fortran_tao_allocate_var_array(
    int &n_var /* 0D_NOT_integer in */,
    bool &default_good_user /* 0D_NOT_logical in */
);
void tao_allocate_var_array(int n_var, bool default_good_user);
extern "C" bool fortran_tao_beam_emit_calc(
    int &plane /* 0D_NOT_integer in */,
    int &emit_type /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    void *bunch_params /* 0D_NOT_type in */,
    double &emit /* 0D_NOT_real out */
);
double
tao_beam_emit_calc(int plane, int emit_type, EleStruct &ele, BunchParamsStruct &bunch_params);
extern "C" void fortran_tao_beam_track(
    void *u /* 0D_NOT_type in */,
    void *tao_lat /* 0D_NOT_type in */,
    int &ix_branch /* 0D_NOT_integer in */,
    void *beam /* 0D_NOT_type inout */,
    bool &calc_ok /* 0D_NOT_logical out */
);
bool tao_beam_track(
    TaoUniverseStruct &u,
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    BeamStruct &beam
);
extern "C" bool fortran_tao_beam_track_endpoint(
    const char *ele_id /* 0D_NOT_character in */,
    void *lat /* 0D_NOT_type in */,
    const char *branch_str /* 0D_NOT_character in */,
    const char *where /* 0D_NOT_character in */,
    void *u /* 0D_NOT_type in */,
    void *ele /* 0D_PTR_type out */
);
std::optional<EleStruct> tao_beam_track_endpoint(
    std::string ele_id,
    LatStruct &lat,
    std::string branch_str,
    std::string where,
    TaoUniverseStruct &u
);
extern "C" bool fortran_tao_branch_index(
    int &ix_branch /* 0D_NOT_integer in */,
    int &ix_this /* 0D_NOT_integer out */
);
int tao_branch_index(int ix_branch);
extern "C" void fortran_tao_calc_data_at_s_pts(
    void *tao_lat /* 0D_NOT_type inout */,
    void *curve /* 0D_NOT_type inout */,
    double &comp_sign /* 0D_NOT_real in */,
    void *good /* 1D_ALLOC_logical inout */
);
void tao_calc_data_at_s_pts(
    TaoLatticeStruct &tao_lat,
    TaoCurveStruct &curve,
    double comp_sign,
    BoolAlloc1D &good
);

// Skipped unusable routine tao_call_cmd:
// - Variable-sized in character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_cbar_wave_anal(void *plot /* 0D_NOT_type inout */);
void tao_cbar_wave_anal(TaoPlotStruct &plot);
extern "C" void fortran_tao_change_ele(
    const char *ele_name /* 0D_NOT_character in */,
    const char *attrib_name /* 0D_NOT_character in */,
    const char *num_str /* 0D_NOT_character in */,
    bool &update /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool tao_change_ele(
    std::string ele_name,
    std::string attrib_name,
    std::string num_str,
    bool update
);
extern "C" void fortran_tao_change_tune(
    const char *branch_str /* 0D_NOT_character in */,
    const char *mask_str /* 0D_NOT_character in */,
    bool &print_list /* 0D_NOT_logical in */,
    const char *dqa_str /* 0D_NOT_character in */,
    const char *dqb_str /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool tao_change_tune(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string dqa_str,
    std::string dqb_str
);
extern "C" void fortran_tao_change_var(
    const char *name /* 0D_NOT_character in */,
    const char *num_str /* 0D_NOT_character in */,
    bool &silent /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool tao_change_var(std::string name, std::string num_str, bool silent);
extern "C" void fortran_tao_change_z_tune(
    const char *branch_str /* 0D_NOT_character in */,
    const char *dq_str /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool tao_change_z_tune(std::string branch_str, std::string dq_str);
extern "C" bool fortran_tao_chrom_calc_needed(
    const char *data_type /* 0D_NOT_character in */,
    const char *data_source /* 0D_NOT_character in */,
    bool &do_chrom /* 0D_NOT_logical in */
);
void tao_chrom_calc_needed(std::string data_type, std::string data_source, bool do_chrom);
extern "C" void fortran_tao_clear_cmd(const char *cmd_line /* 0D_NOT_character in */);
void tao_clear_cmd(std::string cmd_line);
extern "C" void fortran_tao_clip_cmd(
    bool &gang /* 0D_NOT_logical in */,
    const char *where /* 0D_NOT_character in */,
    double &value1 /* 0D_NOT_real in */,
    double &value2 /* 0D_NOT_real in */
);
void tao_clip_cmd(bool gang, std::string where, double value1, double value2);
extern "C" void fortran_tao_close_command_file();
void tao_close_command_file();

// Skipped unusable routine tao_cmd_end_calc:
// - Module name unset
extern "C" void fortran_tao_cmd_history_record(const char *cmd /* 0D_NOT_character in */);
void tao_cmd_history_record(std::string cmd);

// Skipped unusable routine tao_cmd_split:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_command(
    const char *command_line /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical in */,
    bool &err_is_fatal /* 0D_NOT_logical out */
);
bool tao_command(std::string command_line, bool err);
extern "C" bool fortran_tao_constraint_type_name(
    void *datum /* 0D_NOT_type in */,
    const char *datum_name /* 0D_NOT_character out */
);
std::string tao_constraint_type_name(TaoDataStruct &datum);
extern "C" void
fortran_tao_control_tree_list(void *ele /* 0D_NOT_type in */, void *tree /* 1D_ALLOC_type inout */);
void tao_control_tree_list(EleStruct &ele, ElePointerStructAlloc1D tree);
extern "C" void fortran_tao_count_strings(
    const char *string /* 0D_NOT_character in */,
    const char *pattern /* 0D_NOT_character in */,
    int &num /* 0D_NOT_integer out */
);
int tao_count_strings(std::string string, std::string pattern);
extern "C" void fortran_tao_create_plot_window();
void tao_create_plot_window();
extern "C" void fortran_tao_curve_beam_ellipse_setup(void *curve /* 0D_NOT_type inout */);
void tao_curve_beam_ellipse_setup(TaoCurveStruct &curve);
extern "C" bool fortran_tao_curve_check_universe(
    void *curve /* 0D_NOT_type inout */,
    void *uni /* 0D_PTR_type in */,
    bool &is_ok /* 0D_NOT_logical out */
);
bool tao_curve_check_universe(TaoCurveStruct &curve, TaoUniverseStruct &uni);
extern "C" void fortran_tao_curve_data_setup(
    void *plot /* 0D_NOT_type inout */,
    void *graph /* 0D_NOT_type inout */,
    void *curve /* 0D_NOT_type inout */
);
void tao_curve_data_setup(TaoPlotStruct &plot, TaoGraphStruct &graph, TaoCurveStruct &curve);
extern "C" void fortran_tao_curve_datum_calc(
    void *eles /* 1D_ALLOC_type in */,
    void *plot /* 0D_NOT_type in */,
    void *curve /* 0D_NOT_type inout */,
    const char *who /* 0D_NOT_character in */
);
void tao_curve_datum_calc(
    ElePointerStructAlloc1D eles,
    TaoPlotStruct &plot,
    TaoCurveStruct &curve,
    std::string who
);
extern "C" bool fortran_tao_curve_ele_ref(
    void *curve /* 0D_NOT_type in */,
    bool &point_to_ele_ref /* 0D_NOT_logical in */,
    void *ele_track /* 0D_PTR_type inout */
);
void tao_curve_ele_ref(TaoCurveStruct &curve, bool point_to_ele_ref, EleStruct &ele_track);
extern "C" bool
fortran_tao_curve_ix_uni(void *curve /* 0D_NOT_type in */, int &ix_uni /* 0D_NOT_integer out */);
int tao_curve_ix_uni(TaoCurveStruct &curve);
extern "C" bool fortran_tao_curve_name(
    void *curve /* 0D_NOT_type in */,
    bool *use_region /* 0D_NOT_logical in */,
    const char *curve_name /* 0D_NOT_character out */
);
std::string tao_curve_name(TaoCurveStruct &curve, std::optional<bool> use_region = std::nullopt);
extern "C" void fortran_tao_curve_rms_calc(
    void *curve /* 0D_NOT_type in */,
    const char *who /* 0D_NOT_character in */,
    double &rms /* 0D_NOT_real out */,
    double &mean /* 0D_NOT_real out */
);
struct TaoCurveRmsCalc {
  double rms;
  double mean;
};
Tao::TaoCurveRmsCalc tao_curve_rms_calc(TaoCurveStruct &curve, std::string who);
extern "C" bool fortran_tao_d2_d1_name(
    void *d1 /* 0D_NOT_type in */,
    bool *show_universe /* 0D_NOT_logical in */,
    const char *d2_d1_name /* 0D_NOT_character out */
);
std::string tao_d2_d1_name(TaoD1DataStruct &d1, std::optional<bool> show_universe = std::nullopt);
extern "C" void fortran_tao_d2_data_stuffit(
    void *u /* 0D_NOT_type inout */,
    const char *d2_name /* 0D_NOT_character in */,
    int &n_d1_data /* 0D_NOT_integer in */
);
void tao_d2_data_stuffit(TaoUniverseStruct &u, std::string d2_name, int n_d1_data);
extern "C" void fortran_tao_data_check(bool &err /* 0D_NOT_logical in */);
void tao_data_check(bool err);
extern "C" void fortran_tao_data_coupling_init(void *branch /* 0D_NOT_type in */);
void tao_data_coupling_init(BranchStruct &branch);
extern "C" bool fortran_tao_data_sanity_check(
    void *datum /* 0D_NOT_type in */,
    bool &print_err /* 0D_NOT_logical in */,
    const char *default_data_type /* 0D_NOT_character in */,
    void *uni /* 0D_NOT_type in */,
    bool &is_valid /* 0D_NOT_logical out */
);
bool tao_data_sanity_check(
    TaoDataStruct &datum,
    bool print_err,
    std::string default_data_type,
    optional_ref<TaoUniverseStruct> uni = std::nullopt
);

// Skipped unusable routine tao_data_show_use:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_data_type_substitute(
    const char *template_ /* 0D_NOT_character in */,
    const char *str_out /* 0D_NOT_character out */,
    void *curve /* 0D_NOT_type in */,
    void *graph /* 0D_NOT_type in */
);
std::string
tao_data_type_substitute(std::string template_, TaoCurveStruct &curve, TaoGraphStruct &graph);
extern "C" void fortran_tao_data_useit_plot_calc(
    void *curve /* 0D_NOT_type in */,
    void *graph /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &data /* 1D_NOT_type inout */,
    bool &check_s_position /* 0D_NOT_logical in */,
    const char *most_invalid /* 0D_NOT_character out */
);
std::string tao_data_useit_plot_calc(
    TaoCurveStruct &curve,
    TaoGraphStruct &graph,
    TaoDataStructArray1D data,
    bool check_s_position
);
extern "C" bool fortran_tao_datum_has_associated_ele(
    const char *data_type /* 0D_NOT_character in */,
    int *branch_geometry /* 0D_NOT_integer in */,
    int &has_associated_ele /* 0D_NOT_integer out */
);
int tao_datum_has_associated_ele(
    std::string data_type,
    std::optional<int> branch_geometry = std::nullopt
);
extern "C" bool fortran_tao_datum_integrate(
    void *datum /* 0D_NOT_type in */,
    void *branch /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &s_pos /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &values /* 1D_NOT_real in */,
    bool &valid_value /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character out */,
    double &result /* 0D_NOT_real out */
);
struct TaoDatumIntegrate {
  bool valid_value;
  std::string why_invalid;
  double result;
};
Tao::TaoDatumIntegrate tao_datum_integrate(
    TaoDataStruct &datum,
    BranchStruct &branch,
    FArray1D<Real> &s_pos,
    FArray1D<Real> &values
);
extern "C" bool fortran_tao_datum_name(
    void *datum /* 0D_NOT_type in */,
    bool *show_universe /* 0D_NOT_logical in */,
    const char *datum_name /* 0D_NOT_character out */
);
std::string tao_datum_name(TaoDataStruct &datum, std::optional<bool> show_universe = std::nullopt);
extern "C" bool fortran_tao_datum_s_position(
    void *datum /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    double &s_pos /* 0D_NOT_real out */
);
double tao_datum_s_position(TaoDataStruct &datum, EleStruct &ele);
extern "C" void fortran_tao_de_optimizer(bool &abort /* 0D_NOT_logical out */);
bool tao_de_optimizer();
extern "C" void fortran_tao_deallocate_plot_cache(void *plot_cache /* 1D_ALLOC_type inout */);
void tao_deallocate_plot_cache(TaoPlotCacheStructAlloc1D plot_cache);
extern "C" void fortran_tao_deallocate_tree(void *tree /* 0D_NOT_type inout */);
void tao_deallocate_tree(TaoEvalNodeStruct &tree);
extern "C" void fortran_tao_destroy_plot_window();
void tao_destroy_plot_window();
extern "C" void fortran_tao_dmerit_calc();
void tao_dmerit_calc();
extern "C" void fortran_tao_dmodel_dvar_calc(
    bool &force_calc /* 0D_NOT_logical in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool tao_dmodel_dvar_calc(bool force_calc);
extern "C" bool fortran_tao_do_wire_scan(
    void *ele /* 0D_NOT_type in */,
    double &theta /* 0D_NOT_real in */,
    void *beam /* 0D_NOT_type in */,
    double &moment /* 0D_NOT_real out */
);
double tao_do_wire_scan(EleStruct &ele, double theta, BeamStruct &beam);
extern "C" void fortran_tao_draw_beam_chamber_wall(
    void *plot /* 0D_NOT_type in */,
    void *graph /* 0D_NOT_type in */
);
void tao_draw_beam_chamber_wall(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_draw_curve_data(
    void *plot /* 0D_NOT_type in */,
    void *graph /* 0D_NOT_type in */,
    void *curve /* 0D_NOT_type in */,
    bool &have_data /* 0D_NOT_logical inout */
);
void tao_draw_curve_data(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoCurveStruct &curve,
    bool &have_data
);
extern "C" void fortran_tao_draw_ele_for_floor_plan(
    void *plot /* 0D_NOT_type in */,
    void *graph /* 0D_NOT_type in */,
    void *tao_lat /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *ele_shape /* 0D_PTR_type in */,
    const char *label_name /* 0D_NOT_character in */,
    double &offset1 /* 0D_NOT_real in */,
    double &offset2 /* 0D_NOT_real in */
);
void tao_draw_ele_for_floor_plan(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoLatticeStruct &tao_lat,
    EleStruct &ele,
    TaoEleShapeStruct &ele_shape,
    std::string label_name,
    double offset1,
    double offset2
);
extern "C" void
fortran_tao_draw_floor_plan(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_draw_floor_plan(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void
fortran_tao_draw_graph_axes(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_draw_graph_axes(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_draw_histogram_data(
    void *plot /* 0D_NOT_type in */,
    void *graph /* 0D_NOT_type in */,
    void *curve /* 0D_NOT_type in */,
    bool &have_data /* 0D_NOT_logical inout */
);
void tao_draw_histogram_data(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoCurveStruct &curve,
    bool &have_data
);
extern "C" void
fortran_tao_draw_lat_layout(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_draw_lat_layout(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_draw_plots(bool *do_clear /* 0D_NOT_logical in */);
void tao_draw_plots(std::optional<bool> do_clear = std::nullopt);
extern "C" bool fortran_tao_ele_geometry_with_misalignments(
    void *datum /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    bool &valid_value /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character out */,
    double &value /* 0D_NOT_real out */
);
struct TaoEleGeometryWithMisalignments {
  bool valid_value;
  std::string why_invalid;
  double value;
};
Tao::TaoEleGeometryWithMisalignments
tao_ele_geometry_with_misalignments(TaoDataStruct &datum, EleStruct &ele);
extern "C" void fortran_tao_ele_shape_info(
    int &ix_uni /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &ele_shapes /* 1D_NOT_type in */,
    void *e_shape /* 0D_PTR_type out */,
    const char *label_name /* 0D_NOT_character out */,
    double &y1 /* 0D_NOT_real out */,
    double &y2 /* 0D_NOT_real out */,
    int *ix_shape_min /* 0D_NOT_integer inout */
);
struct TaoEleShapeInfo {
  std::optional<TaoEleShapeStruct> e_shape;
  std::string label_name;
  double y1;
  double y2;
};
Tao::TaoEleShapeInfo tao_ele_shape_info(
    int ix_uni,
    EleStruct &ele,
    TaoEleShapeStructArray1D ele_shapes,
    optional_ref<int> ix_shape_min = std::nullopt
);

// Skipped unusable routine tao_ele_shape_input_to_struct:
// - Untranslated type: tao_ele_shape_input (0D)

// Skipped unusable routine tao_ele_shape_struct_to_input:
// - Untranslated type: tao_ele_shape_input (0D)
extern "C" bool fortran_tao_eval_floor_orbit(
    void *datum /* 0D_NOT_type in */,
    void *ele /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    void *bunch_params /* 0D_NOT_type in */,
    bool &valid_value /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character out */,
    double &value /* 0D_NOT_real out */
);
struct TaoEvalFloorOrbit {
  bool valid_value;
  std::string why_invalid;
  double value;
};
Tao::TaoEvalFloorOrbit tao_eval_floor_orbit(
    TaoDataStruct &datum,
    EleStruct &ele,
    CoordStruct &orbit,
    BunchParamsStruct &bunch_params
);
extern "C" void fortran_tao_evaluate_a_datum(
    void *datum /* 0D_NOT_type inout */,
    void *u /* 0D_NOT_type in */,
    void *tao_lat /* 0D_NOT_type in */,
    double &datum_value /* 0D_NOT_real out */,
    bool &valid_value /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character out */,
    bool *called_from_lat_calc /* 0D_NOT_logical in */,
    bool *print_err /* 0D_NOT_logical in */
);
struct TaoEvaluateADatum {
  double datum_value;
  bool valid_value;
  std::string why_invalid;
};
Tao::TaoEvaluateADatum tao_evaluate_a_datum(
    TaoDataStruct &datum,
    TaoUniverseStruct &u,
    TaoLatticeStruct &tao_lat,
    std::optional<bool> called_from_lat_calc = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" bool fortran_tao_evaluate_datum_at_s(
    void *datum /* 0D_NOT_type in */,
    void *tao_lat /* 0D_NOT_type in */,
    void *ele /* 0D_PTR_type in */,
    void *ele_ref /* 0D_PTR_type in */,
    bool &valid_value /* 0D_NOT_logical in */,
    const char *err_str /* 0D_NOT_character out */,
    bool &bad_datum /* 0D_NOT_logical out */,
    double &value /* 0D_NOT_real out */
);
struct TaoEvaluateDatumAtS {
  std::string err_str;
  bool bad_datum;
  double value;
};
Tao::TaoEvaluateDatumAtS tao_evaluate_datum_at_s(
    TaoDataStruct &datum,
    TaoLatticeStruct &tao_lat,
    EleStruct &ele,
    EleStruct &ele_ref,
    bool valid_value
);
extern "C" void fortran_tao_evaluate_element_parameters(
    bool &err /* 0D_NOT_logical out */,
    const char *param_name /* 0D_NOT_character in */,
    void *values /* 1D_ALLOC_real inout */,
    bool &print_err /* 0D_NOT_logical in */,
    void *dflt_ele /* 0D_PTR_type in */,
    const char *dflt_source /* 0D_NOT_character in */,
    const char *dflt_component /* 0D_NOT_character in */,
    int *dflt_uni /* 0D_NOT_integer in */,
    int *eval_point /* 0D_NOT_integer in */,
    void *info /* 1D_ALLOC_type inout */
);
bool tao_evaluate_element_parameters(
    std::string param_name,
    RealAlloc1D &values,
    bool print_err,
    optional_ref<EleStruct> dflt_ele,
    std::string dflt_source,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> eval_point = std::nullopt,
    std::optional<TaoExpressionInfoStructAlloc1D> info = std::nullopt
);
extern "C" void fortran_tao_evaluate_expression(
    const char *expression /* 0D_NOT_character in */,
    int &n_size /* 0D_NOT_integer in */,
    bool &use_good_user /* 0D_NOT_logical in */,
    void *value /* 1D_ALLOC_real inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    void *info /* 1D_ALLOC_type inout */,
    void *stack /* 1D_ALLOC_type inout */,
    const char *dflt_component /* 0D_NOT_character in */,
    const char *dflt_source /* 0D_NOT_character in */,
    void *dflt_ele_ref /* 0D_PTR_type in */,
    void *dflt_ele_start /* 0D_PTR_type in */,
    void *dflt_ele /* 0D_PTR_type in */,
    const char *dflt_dat_or_var_index /* 0D_NOT_character in */,
    int *dflt_uni /* 0D_NOT_integer in */,
    int *dflt_eval_point /* 0D_NOT_integer in */,
    double *dflt_s_offset /* 0D_NOT_real in */,
    void *dflt_orbit /* 0D_NOT_type in */,
    void *datum /* 0D_NOT_type in */
);
bool tao_evaluate_expression(
    std::string expression,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    std::optional<bool> print_err = std::nullopt,
    std::optional<TaoExpressionInfoStructAlloc1D> info = std::nullopt,
    std::optional<TaoEvalNodeStructAlloc1D> stack = std::nullopt,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<std::string> dflt_source = std::nullopt,
    optional_ref<EleStruct> dflt_ele_ref = std::nullopt,
    optional_ref<EleStruct> dflt_ele_start = std::nullopt,
    optional_ref<EleStruct> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_dat_or_var_index = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt,
    optional_ref<CoordStruct> dflt_orbit = std::nullopt,
    optional_ref<TaoDataStruct> datum = std::nullopt
);
extern "C" void fortran_tao_evaluate_expression_new(
    const char *expression /* 0D_NOT_character in */,
    int &n_size /* 0D_NOT_integer in */,
    bool &use_good_user /* 0D_NOT_logical in */,
    void *value /* 1D_ALLOC_real inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    void *info /* 1D_ALLOC_type inout */,
    void *stack /* 1D_ALLOC_type inout */,
    const char *dflt_component /* 0D_NOT_character in */,
    const char *dflt_source /* 0D_NOT_character in */,
    void *dflt_ele_ref /* 0D_PTR_type in */,
    void *dflt_ele_start /* 0D_PTR_type in */,
    void *dflt_ele /* 0D_PTR_type in */,
    const char *dflt_dat_or_var_index /* 0D_NOT_character in */,
    int *dflt_uni /* 0D_NOT_integer in */,
    int *dflt_eval_point /* 0D_NOT_integer in */,
    double *dflt_s_offset /* 0D_NOT_real in */,
    void *dflt_orbit /* 0D_NOT_type in */,
    void *datum /* 0D_NOT_type in */
);
bool tao_evaluate_expression_new(
    std::string expression,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    std::optional<bool> print_err = std::nullopt,
    std::optional<TaoExpressionInfoStructAlloc1D> info = std::nullopt,
    std::optional<TaoEvalNodeStructAlloc1D> stack = std::nullopt,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<std::string> dflt_source = std::nullopt,
    optional_ref<EleStruct> dflt_ele_ref = std::nullopt,
    optional_ref<EleStruct> dflt_ele_start = std::nullopt,
    optional_ref<EleStruct> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_dat_or_var_index = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt,
    optional_ref<CoordStruct> dflt_orbit = std::nullopt,
    optional_ref<TaoDataStruct> datum = std::nullopt
);
extern "C" void fortran_tao_evaluate_expression_old(
    const char *expression /* 0D_NOT_character in */,
    int &n_size /* 0D_NOT_integer in */,
    bool &use_good_user /* 0D_NOT_logical in */,
    void *value /* 1D_ALLOC_real inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */,
    void *info /* 1D_ALLOC_type inout */,
    void *stack /* 1D_ALLOC_type inout */,
    const char *dflt_component /* 0D_NOT_character in */,
    const char *dflt_source /* 0D_NOT_character in */,
    void *dflt_ele_ref /* 0D_PTR_type in */,
    void *dflt_ele_start /* 0D_PTR_type in */,
    void *dflt_ele /* 0D_PTR_type in */,
    const char *dflt_dat_or_var_index /* 0D_NOT_character in */,
    int *dflt_uni /* 0D_NOT_integer in */,
    int *dflt_eval_point /* 0D_NOT_integer in */,
    double *dflt_s_offset /* 0D_NOT_real in */,
    void *dflt_orbit /* 0D_NOT_type in */,
    void *datum /* 0D_NOT_type in */
);
bool tao_evaluate_expression_old(
    std::string expression,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    std::optional<bool> print_err = std::nullopt,
    std::optional<TaoExpressionInfoStructAlloc1D> info = std::nullopt,
    std::optional<TaoEvalNodeStructAlloc1D> stack = std::nullopt,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<std::string> dflt_source = std::nullopt,
    optional_ref<EleStruct> dflt_ele_ref = std::nullopt,
    optional_ref<EleStruct> dflt_ele_start = std::nullopt,
    optional_ref<EleStruct> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_dat_or_var_index = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt,
    optional_ref<CoordStruct> dflt_orbit = std::nullopt,
    optional_ref<TaoDataStruct> datum = std::nullopt
);
extern "C" void fortran_tao_evaluate_lat_or_beam_data(
    bool &err /* 0D_NOT_logical out */,
    const char *data_name /* 0D_NOT_character in */,
    void *values /* 1D_ALLOC_real inout */,
    bool &print_err /* 0D_NOT_logical in */,
    const char *default_source /* 0D_NOT_character in */,
    void *dflt_ele_ref /* 0D_PTR_type in */,
    void *dflt_ele_start /* 0D_PTR_type in */,
    void *dflt_ele /* 0D_PTR_type in */,
    const char *dflt_component /* 0D_NOT_character in */,
    int *dflt_uni /* 0D_NOT_integer in */,
    int *dflt_eval_point /* 0D_NOT_integer in */,
    double *dflt_s_offset /* 0D_NOT_real in */
);
bool tao_evaluate_lat_or_beam_data(
    std::string data_name,
    RealAlloc1D &values,
    bool print_err,
    std::string default_source,
    optional_ref<EleStruct> dflt_ele_ref = std::nullopt,
    optional_ref<EleStruct> dflt_ele_start = std::nullopt,
    optional_ref<EleStruct> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt
);
extern "C" void fortran_tao_evaluate_stack_old(
    Bmad::array_descriptor_t &stack /* 1D_NOT_type in */,
    int &n_size_in /* 0D_NOT_integer in */,
    bool &use_good_user /* 0D_NOT_logical in */,
    void *value /* 1D_ALLOC_real inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool &print_err /* 0D_NOT_logical in */,
    const char *expression /* 0D_NOT_character in */,
    void *info_in /* 1D_ALLOC_type inout */
);
bool tao_evaluate_stack_old(
    TaoEvalNodeStructArray1D stack,
    int n_size_in,
    bool use_good_user,
    RealAlloc1D &value,
    bool print_err,
    std::string expression,
    std::optional<TaoExpressionInfoStructAlloc1D> info_in = std::nullopt
);
extern "C" void fortran_tao_evaluate_tree(
    void *tao_tree /* 0D_NOT_type in */,
    int &n_size /* 0D_NOT_integer in */,
    bool &use_good_user /* 0D_NOT_logical in */,
    void *value /* 1D_ALLOC_real inout */,
    bool &err_flag /* 0D_NOT_logical out */,
    bool &print_err /* 0D_NOT_logical in */,
    const char *expression /* 0D_NOT_character in */,
    void *info_in /* 1D_ALLOC_type inout */
);
bool tao_evaluate_tree(
    TaoEvalNodeStruct &tao_tree,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    bool print_err,
    std::string expression,
    std::optional<TaoExpressionInfoStructAlloc1D> info_in = std::nullopt
);
extern "C" bool fortran_tao_evaluate_tune(
    const char *q_str /* 0D_NOT_character in */,
    double &q0 /* 0D_NOT_real in */,
    bool &delta_input /* 0D_NOT_logical in */,
    double &q_val /* 0D_NOT_real out */
);
double tao_evaluate_tune(std::string q_str, double q0, bool delta_input);
extern "C" void fortran_tao_expression_hash_substitute(
    const char *expression_in /* 0D_NOT_character in */,
    const char *expression_out /* 0D_NOT_character out */,
    void *eval_ele /* 0D_PTR_type in */
);
std::string tao_expression_hash_substitute(
    std::string expression_in,
    optional_ref<EleStruct> eval_ele = std::nullopt
);
extern "C" bool fortran_tao_expression_tree_to_string(
    void *tree /* 0D_NOT_type in */,
    bool *include_root /* 0D_NOT_logical in */,
    int *n_node /* 0D_NOT_integer in */,
    void *parent /* 0D_NOT_type in */,
    const char *str_out /* 0D_ALLOC_character out */
);
std::string tao_expression_tree_to_string(
    TaoEvalNodeStruct &tree,
    std::optional<bool> include_root = std::nullopt,
    std::optional<int> n_node = std::nullopt,
    optional_ref<TaoEvalNodeStruct> parent = std::nullopt
);

// Skipped unusable routine tao_find_data:
// - Untranslated type: tao_d2_data_array_struct (1D)
// - Untranslated type: tao_d1_data_array_struct (1D)
// - Untranslated type: tao_data_array_struct (1D)
// - Untranslated type: tao_real_pointer_struct (1D)
// - Untranslated type: tao_logical_array_struct (1D)
// - Untranslated type: tao_string_array_struct (1D)
// - Untranslated type: tao_integer_array_struct (1D)
extern "C" void fortran_tao_find_plot_region(
    bool &err /* 0D_NOT_logical out */,
    const char *where /* 0D_NOT_character in */,
    void *region /* 0D_PTR_type out */,
    bool *print_flag /* 0D_NOT_logical in */
);
struct TaoFindPlotRegion {
  bool err;
  std::optional<TaoPlotRegionStruct> region;
};
Tao::TaoFindPlotRegion
tao_find_plot_region(std::string where, std::optional<bool> print_flag = std::nullopt);

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
    const char *switch_ /* 0D_NOT_character in */,
    const char *word1 /* 0D_NOT_character in */,
    const char *word2 /* 0D_NOT_character in */
);
void tao_fixer(std::string switch_, std::string word1, std::string word2);
extern "C" void fortran_tao_floor_to_screen(
    void *graph /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &r_floor /* 1D_NOT_real inout */,
    double &x_screen /* 0D_NOT_real out */,
    double &y_screen /* 0D_NOT_real out */
);
struct TaoFloorToScreen {
  double x_screen;
  double y_screen;
};
Tao::TaoFloorToScreen tao_floor_to_screen(TaoGraphStruct &graph, FixedArray1D<Real, 3> r_floor);
extern "C" void fortran_tao_floor_to_screen_coords(
    void *graph /* 0D_NOT_type in */,
    void *floor /* 0D_NOT_type in */,
    void *screen /* 0D_NOT_type out */
);
FloorPositionStruct tao_floor_to_screen_coords(TaoGraphStruct &graph, FloorPositionStruct &floor);

// Skipped unusable routine tao_geo_lm_func:
// - Routine in configuration skip list
extern "C" void fortran_tao_geodesic_lm_optimizer(bool &abort /* 0D_NOT_logical out */);
bool tao_geodesic_lm_optimizer();
extern "C" void fortran_tao_get_data(
    void *data_value /* 1D_ALLOC_real inout */,
    void *data_weight /* 1D_ALLOC_real inout */,
    void *data_meas_value /* 1D_ALLOC_real inout */,
    void *data_ix_dModel /* 1D_ALLOC_integer inout */
);
void tao_get_data(
    optional_ref<RealAlloc1D> data_value = std::nullopt,
    optional_ref<RealAlloc1D> data_weight = std::nullopt,
    optional_ref<RealAlloc1D> data_meas_value = std::nullopt,
    optional_ref<IntAlloc1D> data_ix_dModel = std::nullopt
);
extern "C" void fortran_tao_get_opt_vars(
    void *var_value /* 1D_ALLOC_real inout */,
    void *var_step /* 1D_ALLOC_real inout */,
    void *var_delta /* 1D_ALLOC_real inout */,
    void *var_weight /* 1D_ALLOC_real inout */,
    void *var_ix /* 1D_ALLOC_integer inout */,
    bool &ignore_if_weight_is_zero /* 0D_NOT_logical out */,
    bool &ignore_if_not_limited /* 0D_NOT_logical out */
);
struct TaoGetOptVars {
  bool ignore_if_weight_is_zero;
  bool ignore_if_not_limited;
};
Tao::TaoGetOptVars tao_get_opt_vars(
    optional_ref<RealAlloc1D> var_value = std::nullopt,
    optional_ref<RealAlloc1D> var_step = std::nullopt,
    optional_ref<RealAlloc1D> var_delta = std::nullopt,
    optional_ref<RealAlloc1D> var_weight = std::nullopt,
    optional_ref<IntAlloc1D> var_ix = std::nullopt
);
extern "C" void fortran_tao_get_user_input(
    const char *cmd_out /* 0D_NOT_character out */,
    const char *prompt_str /* 0D_NOT_character in */,
    bool *wait_flag /* 0D_NOT_logical in */,
    const char *cmd_in /* 0D_NOT_character in */
);
std::string tao_get_user_input(
    std::optional<std::string> prompt_str = std::nullopt,
    std::optional<bool> wait_flag = std::nullopt,
    std::optional<std::string> cmd_in = std::nullopt
);
extern "C" void fortran_tao_graph_controller_setup(void *graph /* 0D_NOT_type inout */);
void tao_graph_controller_setup(TaoGraphStruct &graph);
extern "C" void fortran_tao_graph_data_setup(
    void *plot /* 0D_NOT_type inout */,
    void *graph /* 0D_NOT_type inout */
);
void tao_graph_data_setup(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_graph_data_slice_setup(
    void *plot /* 0D_NOT_type inout */,
    void *graph /* 0D_NOT_type inout */
);
void tao_graph_data_slice_setup(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_graph_dynamic_aperture_setup(
    void *plot /* 0D_NOT_type inout */,
    void *graph /* 0D_NOT_type inout */
);
void tao_graph_dynamic_aperture_setup(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_graph_histogram_setup(
    void *plot /* 0D_NOT_type inout */,
    void *graph /* 0D_NOT_type inout */
);
void tao_graph_histogram_setup(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" bool fortran_tao_graph_name(
    void *graph /* 0D_NOT_type in */,
    bool *use_region /* 0D_NOT_logical in */,
    const char *graph_name /* 0D_NOT_character out */
);
std::string tao_graph_name(TaoGraphStruct &graph, std::optional<bool> use_region = std::nullopt);
extern "C" void fortran_tao_graph_phase_space_setup(
    void *plot /* 0D_NOT_type inout */,
    void *graph /* 0D_NOT_type inout */
);
void tao_graph_phase_space_setup(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_graph_s_min_max_calc(
    void *graph /* 0D_NOT_type in */,
    void *branch /* 0D_NOT_type in */,
    double &s_min /* 0D_NOT_real out */,
    double &s_max /* 0D_NOT_real out */
);
struct TaoGraphSMinMaxCalc {
  double s_min;
  double s_max;
};
Tao::TaoGraphSMinMaxCalc tao_graph_s_min_max_calc(TaoGraphStruct &graph, BranchStruct &branch);
extern "C" void
fortran_tao_graph_setup(void *plot /* 0D_NOT_type inout */, void *graph /* 0D_NOT_type inout */);
void tao_graph_setup(TaoPlotStruct &plot, TaoGraphStruct &graph);

// Skipped unusable routine tao_help:
// - Variable-sized inout character array: 1D_ALLOC_character
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

// Skipped unusable routine tao_hook_universe_calc_post_process_def:
// - Routine in configuration skip list
extern "C" void fortran_tao_init(bool &err_flag /* 0D_NOT_logical out */);
bool tao_init();
extern "C" void fortran_tao_init_beam_in_universe(
    void *u /* 0D_NOT_type inout */,
    void *beam_init /* 0D_NOT_type inout */,
    const char *track_start /* 0D_NOT_character in */,
    const char *track_end /* 0D_NOT_character in */,
    double &comb_ds_save /* 0D_NOT_real in */
);
void tao_init_beam_in_universe(
    TaoUniverseStruct &u,
    BeamInitStruct &beam_init,
    std::string track_start,
    std::string track_end,
    double comb_ds_save
);
extern "C" void fortran_tao_init_beams(const char *init_file /* 0D_NOT_character in */);
void tao_init_beams(std::string init_file);

// Skipped unusable routine tao_init_building_wall:
// - Module name unset
extern "C" void fortran_tao_init_data(const char *data_file /* 0D_NOT_character in */);
void tao_init_data(std::string data_file);
extern "C" void fortran_tao_init_data_end_stuff();
void tao_init_data_end_stuff();
extern "C" void fortran_tao_init_data_in_universe(
    void *u /* 0D_NOT_type inout */,
    int &n_d2_add /* 0D_NOT_integer in */,
    bool *keep_existing_data /* 0D_NOT_logical in */
);
void tao_init_data_in_universe(
    TaoUniverseStruct &u,
    int n_d2_add,
    std::optional<bool> keep_existing_data = std::nullopt
);
extern "C" void fortran_tao_init_dynamic_aperture(const char *init_file /* 0D_NOT_character in */);
void tao_init_dynamic_aperture(std::string init_file);
extern "C" void fortran_tao_init_find_elements(
    void *u /* 0D_NOT_type in */,
    const char *search_string /* 0D_NOT_character in */,
    void *eles /* 1D_ALLOC_type inout */,
    const char *attribute /* 0D_NOT_character in */,
    bool &found_one /* 0D_NOT_logical out */
);
bool tao_init_find_elements(
    TaoUniverseStruct &u,
    std::string search_string,
    ElePointerStructAlloc1D eles,
    std::optional<std::string> attribute = std::nullopt
);
extern "C" void fortran_tao_init_global(const char *init_file /* 0D_NOT_character in */);
void tao_init_global(std::string init_file);
extern "C" void fortran_tao_init_lattice(
    const char *lat_file /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void tao_init_lattice(std::string lat_file, bool err_flag);
extern "C" void fortran_tao_init_plotting(const char *plot_file /* 0D_NOT_character in */);
void tao_init_plotting(std::string plot_file);
extern "C" void fortran_tao_init_variables(const char *var_file /* 0D_NOT_character in */);
void tao_init_variables(std::string var_file);
extern "C" void fortran_tao_inject_beam(
    void *u /* 0D_NOT_type in */,
    void *model /* 0D_NOT_type in */,
    int &ix_branch /* 0D_NOT_integer in */,
    void *beam /* 0D_NOT_type out */,
    bool &init_ok /* 0D_NOT_logical out */
);
struct TaoInjectBeam {
  BeamStruct beam;
  bool init_ok;
};
Tao::TaoInjectBeam tao_inject_beam(TaoUniverseStruct &u, TaoLatticeStruct &model, int ix_branch);
extern "C" void fortran_tao_inject_particle(
    void *u /* 0D_NOT_type inout */,
    void *model /* 0D_NOT_type inout */,
    int &ix_branch /* 0D_NOT_integer in */
);
void tao_inject_particle(TaoUniverseStruct &u, TaoLatticeStruct &model, int ix_branch);
extern "C" bool fortran_tao_is_valid_name(
    const char *name /* 0D_NOT_character in */,
    const char *why_invalid /* 0D_NOT_character out */,
    bool &is_valid /* 0D_NOT_logical out */
);
struct TaoIsValidName {
  std::string why_invalid;
  bool is_valid;
};
Tao::TaoIsValidName tao_is_valid_name(std::string name);
extern "C" void fortran_tao_json_cmd(const char *input_str /* 0D_NOT_character in */);
void tao_json_cmd(std::string input_str);
extern "C" void fortran_tao_key_info_to_str(
    int &ix_key /* 0D_NOT_integer in */,
    int &ix_min_key /* 0D_NOT_integer in */,
    int &ix_max_key /* 0D_NOT_integer in */,
    const char *key_str /* 0D_NOT_character in */,
    const char *header_str /* 0D_NOT_character in */
);
void tao_key_info_to_str(
    int ix_key,
    int ix_min_key,
    int ix_max_key,
    std::string key_str,
    std::string header_str
);
extern "C" void fortran_tao_lat_bookkeeper(
    void *u /* 0D_NOT_type in */,
    void *tao_lat /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */
);
bool tao_lat_bookkeeper(TaoUniverseStruct &u, TaoLatticeStruct &tao_lat);
extern "C" bool fortran_tao_lat_emit_calc(
    int &plane /* 0D_NOT_integer in */,
    int &emit_type /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    void *modes /* 0D_NOT_type in */,
    double &emit /* 0D_NOT_real out */
);
double tao_lat_emit_calc(int plane, int emit_type, EleStruct &ele, NormalModesStruct &modes);
extern "C" bool fortran_tao_lat_sigma_calc_needed(
    const char *data_type /* 0D_NOT_character in */,
    const char *data_source /* 0D_NOT_character in */,
    bool &do_lat_sigma /* 0D_NOT_logical in */
);
void tao_lat_sigma_calc_needed(std::string data_type, std::string data_source, bool do_lat_sigma);
extern "C" void fortran_tao_lat_sigma_track(
    void *tao_lat /* 0D_NOT_type in */,
    bool &calc_ok /* 0D_NOT_logical out */,
    int &ix_branch /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */,
    bool *force_calc /* 0D_NOT_logical in */
);
bool tao_lat_sigma_track(
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    std::optional<bool> print_err = std::nullopt,
    std::optional<bool> force_calc = std::nullopt
);
extern "C" void fortran_tao_lattice_branches_equal_tao_lattice_branches(
    Bmad::array_descriptor_t &tlb1 /* 1D_NOT_type inout */,
    Bmad::array_descriptor_t &tlb2 /* 1D_NOT_type in */
);
void tao_lattice_branches_equal_tao_lattice_branches(
    TaoLatticeBranchStructArray1D tlb1,
    TaoLatticeBranchStructArray1D tlb2
);
extern "C" void fortran_tao_lattice_calc(
    bool &calc_ok /* 0D_NOT_logical out */,
    bool &print_err /* 0D_NOT_logical out */
);
struct TaoLatticeCalc {
  bool calc_ok;
  bool print_err;
};
Tao::TaoLatticeCalc tao_lattice_calc();
extern "C" void fortran_tao_lattice_equal_tao_lattice(
    void *lat1 /* 0D_NOT_type inout */,
    void *lat2 /* 0D_NOT_type in */
);
void tao_lattice_equal_tao_lattice(TaoLatticeStruct &lat1, TaoLatticeStruct &lat2);
extern "C" void fortran_tao_limit_calc(bool &limited /* 0D_NOT_logical out */);
bool tao_limit_calc();
extern "C" void fortran_tao_lm_optimizer(bool &abort /* 0D_NOT_logical out */);
bool tao_lm_optimizer();
extern "C" void fortran_tao_lmdif_optimizer(bool &abort /* 0D_NOT_logical out */);
bool tao_lmdif_optimizer();
extern "C" void fortran_tao_load_this_datum(
    Bmad::array_descriptor_t &vec /* 1D_NOT_real inout */,
    void *ele_ref /* 0D_PTR_type inout */,
    void *ele_start /* 0D_PTR_type inout */,
    void *ele /* 0D_PTR_type inout */,
    double &datum_value /* 0D_NOT_real in */,
    bool &valid_value /* 0D_NOT_logical in */,
    void *datum /* 0D_NOT_type inout */,
    void *branch /* 0D_NOT_type inout */,
    const char *why_invalid /* 0D_NOT_character in */,
    void *good /* 1D_ALLOC_logical inout */
);
void tao_load_this_datum(
    FArray1D<Real> &vec,
    EleStruct &ele_ref,
    EleStruct &ele_start,
    EleStruct &ele,
    double datum_value,
    bool valid_value,
    TaoDataStruct &datum,
    BranchStruct &branch,
    std::optional<std::string> why_invalid = std::nullopt,
    optional_ref<BoolAlloc1D> good = std::nullopt
);
extern "C" void fortran_tao_locate_all_elements(
    const char *ele_list /* 0D_NOT_character in */,
    void *eles /* 1D_ALLOC_type inout */,
    bool &err /* 0D_NOT_logical out */,
    bool *ignore_blank /* 0D_NOT_logical in */
);
bool tao_locate_all_elements(
    std::string ele_list,
    ElePointerStructAlloc1D eles,
    std::optional<bool> ignore_blank = std::nullopt
);
extern "C" void fortran_tao_locate_elements(
    const char *ele_list /* 0D_NOT_character in */,
    int &ix_universe /* 0D_NOT_integer in */,
    void *eles /* 1D_ALLOC_type inout */,
    bool &err /* 0D_NOT_logical out */,
    int *lat_type /* 0D_NOT_integer in */,
    bool *ignore_blank /* 0D_NOT_logical in */,
    int *err_stat_level /* 0D_NOT_integer in */,
    bool *above_ubound_is_err /* 0D_NOT_logical in */,
    int *ix_branch /* 0D_NOT_integer in */,
    bool *multiple_eles_is_err /* 0D_NOT_logical in */
);
bool tao_locate_elements(
    std::string ele_list,
    int ix_universe,
    ElePointerStructAlloc1D eles,
    std::optional<int> lat_type = std::nullopt,
    std::optional<bool> ignore_blank = std::nullopt,
    std::optional<int> err_stat_level = std::nullopt,
    std::optional<bool> above_ubound_is_err = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    std::optional<bool> multiple_eles_is_err = std::nullopt
);
extern "C" void fortran_tao_mark_lattice_ele(void *lat /* 0D_NOT_type inout */);
void tao_mark_lattice_ele(LatStruct &lat);
extern "C" bool
fortran_tao_merit(bool &calc_ok /* 0D_NOT_logical out */, double &this_merit /* 0D_NOT_real out */);
struct TaoMerit {
  bool calc_ok;
  double this_merit;
};
Tao::TaoMerit tao_merit();

// Skipped unusable routine tao_mrq_func:
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine tao_next_switch:
// - Variable-sized in character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_next_word(
    const char *line /* 0D_NOT_character inout */,
    const char *word /* 0D_NOT_character out */
);
std::string tao_next_word(std::string &line);
extern "C" bool fortran_tao_one_turn_map_calc_needed(
    const char *data_type /* 0D_NOT_character in */,
    const char *data_source /* 0D_NOT_character in */,
    bool &do_one_turn_map /* 0D_NOT_logical in */
);
void tao_one_turn_map_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_one_turn_map
);
extern "C" void fortran_tao_open_file(
    const char *file /* 0D_NOT_character in */,
    int &iunit /* 0D_NOT_integer out */,
    const char *file_name /* 0D_NOT_character in */,
    int &error_severity /* 0D_NOT_integer in */,
    bool *binary /* 0D_NOT_logical in */
);
int tao_open_file(
    std::string file,
    std::string file_name,
    int error_severity,
    std::optional<bool> binary = std::nullopt
);
extern "C" bool
fortran_tao_open_scratch_file(bool &err /* 0D_NOT_logical out */, int &iu /* 0D_NOT_integer out */);
struct TaoOpenScratchFile {
  bool err;
  int iu;
};
Tao::TaoOpenScratchFile tao_open_scratch_file();
extern "C" bool fortran_tao_optimization_status(
    void *datum /* 0D_NOT_type in */,
    const char *why_str /* 0D_NOT_character out */
);
std::string tao_optimization_status(TaoDataStruct &datum);
extern "C" void fortran_tao_orbit_beta_wave_anal(void *plot /* 0D_NOT_type inout */);
void tao_orbit_beta_wave_anal(TaoPlotStruct &plot);
extern "C" bool fortran_tao_oreint_building_wall_pt(
    void *pt_in /* 0D_NOT_type in */,
    void *pt_out /* 0D_NOT_type out */
);
TaoBuildingWallPointStruct tao_oreint_building_wall_pt(TaoBuildingWallPointStruct &pt_in);
extern "C" bool fortran_tao_param_value_at_s(
    const char *dat_name /* 0D_NOT_character in */,
    void *ele_to_s /* 0D_NOT_type in */,
    void *ele_here /* 0D_NOT_type in */,
    void *orbit /* 0D_NOT_type in */,
    bool &err_flag /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character out */,
    bool &print_err /* 0D_NOT_logical out */,
    bool &bad_datum /* 0D_NOT_logical out */,
    double &value /* 0D_NOT_real out */
);
struct TaoParamValueAtS {
  bool err_flag;
  std::string why_invalid;
  bool print_err;
  bool bad_datum;
  double value;
};
Tao::TaoParamValueAtS tao_param_value_at_s(
    std::string dat_name,
    EleStruct &ele_to_s,
    EleStruct &ele_here,
    CoordStruct &orbit
);
extern "C" void fortran_tao_param_value_routine(
    const char *str /* 0D_NOT_character in */,
    bool &use_good_user /* 0D_NOT_logical in */,
    const char *saved_prefix /* 0D_NOT_character in */,
    void *stack /* 0D_NOT_type inout */,
    bool &err_flag /* 0D_NOT_logical in */,
    bool &print_err /* 0D_NOT_logical in */,
    const char *dflt_component /* 0D_NOT_character in */,
    const char *dflt_source /* 0D_NOT_character in */,
    void *dflt_ele_ref /* 0D_PTR_type inout */,
    void *dflt_ele_start /* 0D_PTR_type inout */,
    void *dflt_ele /* 0D_PTR_type inout */,
    const char *dflt_dat_or_var_index /* 0D_NOT_character in */,
    int *dflt_uni /* 0D_NOT_integer in */,
    int *dflt_eval_point /* 0D_NOT_integer in */,
    double *dflt_s_offset /* 0D_NOT_real in */,
    void *dflt_orbit /* 0D_NOT_type inout */,
    void *datum /* 0D_NOT_type inout */
);
void tao_param_value_routine(
    std::string str,
    bool use_good_user,
    std::string saved_prefix,
    TaoEvalNodeStruct &stack,
    bool err_flag,
    bool print_err,
    std::optional<std::string> dflt_component = std::nullopt,
    std::optional<std::string> dflt_source = std::nullopt,
    optional_ref<EleStruct> dflt_ele_ref = std::nullopt,
    optional_ref<EleStruct> dflt_ele_start = std::nullopt,
    optional_ref<EleStruct> dflt_ele = std::nullopt,
    std::optional<std::string> dflt_dat_or_var_index = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<int> dflt_eval_point = std::nullopt,
    std::optional<double> dflt_s_offset = std::nullopt,
    optional_ref<CoordStruct> dflt_orbit = std::nullopt,
    optional_ref<TaoDataStruct> datum = std::nullopt
);
extern "C" void fortran_tao_parse_command_args(
    bool &error /* 0D_NOT_logical out */,
    const char *cmd_line /* 0D_NOT_character in */
);
bool tao_parse_command_args(std::optional<std::string> cmd_line = std::nullopt);
extern "C" void fortran_tao_parse_element_param_str(
    bool &err /* 0D_NOT_logical out */,
    const char *in_str /* 0D_NOT_character in */,
    const char *uni /* 0D_NOT_character out */,
    const char *element /* 0D_NOT_character out */,
    const char *parameter /* 0D_NOT_character out */,
    int &where /* 0D_NOT_integer out */,
    const char *component /* 0D_NOT_character out */
);
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
    const char *data_type /* 0D_NOT_character in */,
    Bmad::array_descriptor_t &p /* 1D_NOT_type in */,
    void *value /* 1D_ALLOC_real inout */,
    bool &err /* 0D_NOT_logical out */,
    void *ele /* 0D_NOT_type in */,
    int &ix_bunch /* 0D_NOT_integer in */
);
bool tao_particle_data_value(
    std::string data_type,
    CoordStructArray1D p,
    RealAlloc1D &value,
    EleStruct &ele,
    int ix_bunch
);
extern "C" void fortran_tao_pause_cmd(double &time /* 0D_NOT_real in */);
void tao_pause_cmd(double time);
extern "C" bool fortran_tao_phase_space_axis_index(
    const char *data_type /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical in */,
    int &ix_axis /* 0D_NOT_integer out */
);
int tao_phase_space_axis_index(std::string data_type, bool err);
extern "C" void fortran_tao_phase_wave_anal(void *plot /* 0D_NOT_type inout */);
void tao_phase_wave_anal(TaoPlotStruct &plot);
extern "C" void fortran_tao_pick_universe(
    const char *name_in /* 0D_NOT_character in */,
    const char *name_out /* 0D_NOT_character out */,
    void *picked /* 1D_ALLOC_logical inout */,
    bool &err /* 0D_NOT_logical out */,
    int &ix_uni /* 0D_NOT_integer out */,
    bool &explicit_uni /* 0D_NOT_logical out */,
    int *dflt_uni /* 0D_NOT_integer in */,
    bool *pure_uni /* 0D_NOT_logical in */
);
struct TaoPickUniverse {
  std::string name_out;
  bool err;
  int ix_uni;
  bool explicit_uni;
};
Tao::TaoPickUniverse tao_pick_universe(
    std::string name_in,
    BoolAlloc1D &picked,
    std::optional<int> dflt_uni = std::nullopt,
    std::optional<bool> pure_uni = std::nullopt
);
extern "C" void fortran_tao_pipe_cmd(const char *input_str /* 0D_NOT_character in */);
void tao_pipe_cmd(std::string input_str);
extern "C" void fortran_tao_place_cmd(
    const char *where /* 0D_NOT_character in */,
    const char *who /* 0D_NOT_character in */,
    bool *no_buffer /* 0D_NOT_logical in */
);
void tao_place_cmd(
    std::string where,
    std::string who,
    std::optional<bool> no_buffer = std::nullopt
);
extern "C" void fortran_tao_plot_cmd(
    const char *where /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */
);
void tao_plot_cmd(std::string where, std::string component);
extern "C" void
fortran_tao_plot_data(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_plot_data(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void
fortran_tao_plot_histogram(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_plot_histogram(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void
fortran_tao_plot_key_table(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_plot_key_table(TaoPlotStruct &plot, TaoGraphStruct &graph);
extern "C" void fortran_tao_plot_setup();
void tao_plot_setup();
extern "C" void fortran_tao_plot_struct_transfer(
    void *plot_in /* 0D_NOT_type in */,
    void *plot_out /* 0D_NOT_type out */
);
TaoPlotStruct tao_plot_struct_transfer(TaoPlotStruct &plot_in);
extern "C" void
fortran_tao_plot_wave(void *plot /* 0D_NOT_type in */, void *graph /* 0D_NOT_type in */);
void tao_plot_wave(TaoPlotStruct &plot, TaoGraphStruct &graph);

// Skipped unusable routine tao_point_d1_to_data:
// - Array bounds handling: "Enum 'N_MIN' found in bounds 'n_min' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_point_v1_to_var:
// - Array bounds handling: "Enum 'N' found in bounds 'n' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_pointer_to_branches:
// - Untranslated type: branch_pointer_struct (1D)
extern "C" bool fortran_tao_pointer_to_building_wall_shape(
    const char *wall_name /* 0D_NOT_character in */,
    void *e_shape /* 0D_PTR_type out */
);
std::optional<TaoEleShapeStruct> tao_pointer_to_building_wall_shape(std::string wall_name);
extern "C" bool fortran_tao_pointer_to_datum(
    void *d1 /* 0D_NOT_type in */,
    const char *ele_name /* 0D_NOT_character in */,
    void *datum_ptr /* 0D_PTR_type out */
);
std::optional<TaoDataStruct> tao_pointer_to_datum(TaoD1DataStruct &d1, std::string ele_name);
extern "C" bool fortran_tao_pointer_to_datum_ele(
    void *lat /* 0D_NOT_type in */,
    const char *ele_name /* 0D_NOT_character in */,
    int &ix_ele /* 0D_NOT_integer in */,
    void *datum /* 0D_NOT_type in */,
    bool &valid /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character out */,
    bool *print_err /* 0D_NOT_logical in */,
    void *ele /* 0D_PTR_type out */
);
struct TaoPointerToDatumEle {
  bool valid;
  std::string why_invalid;
  std::optional<EleStruct> ele;
};
Tao::TaoPointerToDatumEle tao_pointer_to_datum_ele(
    LatStruct &lat,
    std::string ele_name,
    int ix_ele,
    TaoDataStruct &datum,
    std::optional<bool> print_err = std::nullopt
);
extern "C" bool fortran_tao_pointer_to_ele_shape(
    int &ix_uni /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type in */,
    Bmad::array_descriptor_t &ele_shape /* 1D_NOT_type in */,
    const char *dat_var_name /* 0D_NOT_character out */,
    double &dat_var_value /* 0D_NOT_real out */,
    int *ix_shape_min /* 0D_NOT_integer inout */,
    void *e_shape /* 0D_PTR_type out */
);
struct TaoPointerToEleShape {
  std::string dat_var_name;
  double dat_var_value;
  std::optional<TaoEleShapeStruct> e_shape;
};
Tao::TaoPointerToEleShape tao_pointer_to_ele_shape(
    int ix_uni,
    EleStruct &ele,
    TaoEleShapeStructArray1D ele_shape,
    optional_ref<int> ix_shape_min = std::nullopt
);
extern "C" bool fortran_tao_pointer_to_tao_lat(
    void *u /* 0D_NOT_type in */,
    int *lat_type /* 0D_NOT_integer in */,
    void *tao_lat /* 0D_PTR_type out */
);
std::optional<TaoLatticeStruct>
tao_pointer_to_tao_lat(TaoUniverseStruct &u, std::optional<int> lat_type = std::nullopt);
extern "C" bool fortran_tao_pointer_to_universe_int(
    int &ix_uni /* 0D_NOT_integer in */,
    bool *neg2_to_default /* 0D_NOT_logical in */,
    void *u /* 0D_PTR_type out */
);
std::optional<TaoUniverseStruct>
tao_pointer_to_universe(int ix_uni, std::optional<bool> neg2_to_default = std::nullopt);
extern "C" bool fortran_tao_pointer_to_universe_str(
    const char *string /* 0D_NOT_character inout */,
    bool *neg2_to_default /* 0D_NOT_logical in */,
    void *u /* 0D_PTR_type out */
);
std::optional<TaoUniverseStruct>
tao_pointer_to_universe(std::string &string, std::optional<bool> neg2_to_default = std::nullopt);
extern "C" void fortran_tao_pointer_to_universes(
    const char *name_in /* 0D_NOT_character in */,
    void *unis /* 1D_ALLOC_type inout */,
    bool &err /* 0D_NOT_logical out */,
    const char *name_out /* 0D_NOT_character out */,
    bool &explicit_uni /* 0D_NOT_logical out */,
    int *dflt_uni /* 0D_NOT_integer in */
);
struct TaoPointerToUniverses {
  bool err;
  std::string name_out;
  bool explicit_uni;
};
Tao::TaoPointerToUniverses tao_pointer_to_universes(
    std::string name_in,
    TaoUniversePointerStructAlloc1D unis,
    std::optional<int> dflt_uni = std::nullopt
);
extern "C" void fortran_tao_pointer_to_var_in_lattice(
    void *var /* 0D_NOT_type in */,
    int &ix_uni /* 0D_NOT_integer in */,
    void *ele /* 0D_NOT_type inout */,
    bool &err /* 0D_NOT_logical out */
);
bool tao_pointer_to_var_in_lattice(TaoVarStruct &var, int ix_uni, EleStruct &ele);
extern "C" void fortran_tao_pointer_to_var_in_lattice2(
    void *var /* 0D_NOT_type in */,
    int &ix_uni /* 0D_NOT_integer in */,
    bool &err /* 0D_NOT_logical out */
);
bool tao_pointer_to_var_in_lattice2(TaoVarStruct &var, int ix_uni);
extern "C" void fortran_tao_print_command_line_info();
void tao_print_command_line_info();

// Skipped unusable routine tao_print_vars:
// - Untranslated type: tao_var_array_struct (1D)

// Skipped unusable routine tao_ptc_cmd:
// - Module name unset
extern "C" void fortran_tao_ptc_normal_form(
    bool &do_calc /* 0D_NOT_logical in */,
    void *tao_lat /* 0D_NOT_type in */,
    int &ix_branch /* 0D_NOT_integer in */,
    int *rf_on /* 0D_NOT_integer in */
);
void tao_ptc_normal_form(
    bool do_calc,
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    std::optional<int> rf_on = std::nullopt
);
extern "C" void fortran_tao_python_cmd(const char *input_str /* 0D_NOT_character in */);
void tao_python_cmd(std::string input_str);
extern "C" void fortran_tao_quiet_set(const char *set /* 0D_NOT_character in */);
void tao_quiet_set(std::string set);
extern "C" bool fortran_tao_rad_int_calc_needed(
    const char *data_type /* 0D_NOT_character in */,
    const char *data_source /* 0D_NOT_character in */,
    bool &do_rad_int /* 0D_NOT_logical in */
);
void tao_rad_int_calc_needed(std::string data_type, std::string data_source, bool do_rad_int);
extern "C" void fortran_tao_re_allocate_expression_info(
    void *info /* 1D_ALLOC_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void tao_re_allocate_expression_info(
    TaoExpressionInfoStructAlloc1D info,
    int n,
    std::optional<bool> exact = std::nullopt
);
extern "C" void fortran_tao_re_associate_node_array(
    void *tree /* 0D_NOT_type inout */,
    int &n /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void tao_re_associate_node_array(
    TaoEvalNodeStruct &tree,
    int n,
    std::optional<bool> exact = std::nullopt
);
extern "C" void fortran_tao_re_execute(
    const char *string /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical in */
);
void tao_re_execute(std::string string, bool err);
extern "C" void fortran_tao_read_cmd(
    const char *which /* 0D_NOT_character in */,
    const char *unis /* 0D_NOT_character in */,
    const char *file /* 0D_NOT_character in */,
    bool &silent /* 0D_NOT_logical in */
);
void tao_read_cmd(std::string which, std::string unis, std::string file, bool silent);

// Skipped unusable routine tao_read_in_patterns:
// - Module name unset
extern "C" bool fortran_tao_read_phase_space_index(
    const char *name /* 0D_NOT_character in */,
    int &ixc /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */,
    int &ix_ps /* 0D_NOT_integer out */
);
int tao_read_phase_space_index(
    std::string name,
    int ixc,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_tao_regression_test(const char *cmd_str /* 0D_NOT_character in */);
void tao_regression_test(std::string cmd_str);
extern "C" void fortran_tao_remove_blank_characters(const char *str /* 0D_NOT_character inout */);
void tao_remove_blank_characters(std::string &str);
extern "C" void fortran_tao_run_cmd(
    const char *which /* 0D_NOT_character in */,
    bool &abort /* 0D_NOT_logical out */
);
bool tao_run_cmd(std::string which);
extern "C" void fortran_tao_scale_cmd(
    const char *where /* 0D_NOT_character in */,
    double &y_min_in /* 0D_NOT_real in */,
    double &y_max_in /* 0D_NOT_real in */,
    const char *axis /* 0D_NOT_character in */,
    bool *include_wall /* 0D_NOT_logical in */,
    const char *gang /* 0D_NOT_character in */,
    bool *exact /* 0D_NOT_logical in */,
    bool *turn_autoscale_off /* 0D_NOT_logical in */
);
void tao_scale_cmd(
    std::string where,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> exact = std::nullopt,
    std::optional<bool> turn_autoscale_off = std::nullopt
);
extern "C" void fortran_tao_scale_graph(
    void *graph /* 0D_NOT_type inout */,
    double &y_min /* 0D_NOT_real in */,
    double &y_max /* 0D_NOT_real in */,
    const char *axis /* 0D_NOT_character in */,
    bool *include_wall /* 0D_NOT_logical in */,
    Bmad::array_descriptor_t &y_range /* 1D_NOT_real out */,
    Bmad::array_descriptor_t &y2_range /* 1D_NOT_real out */
);
struct TaoScaleGraph {
  FixedArray1D<Real, 2> y_range;
  FixedArray1D<Real, 2> y2_range;
};
Tao::TaoScaleGraph tao_scale_graph(
    TaoGraphStruct &graph,
    double y_min,
    double y_max,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt
);
extern "C" void fortran_tao_scale_ping_data(void *u /* 0D_NOT_type inout */);
void tao_scale_ping_data(TaoUniverseStruct &u);
extern "C" void fortran_tao_scale_plot(
    void *plot /* 0D_NOT_type inout */,
    double &y_min_in /* 0D_NOT_real in */,
    double &y_max_in /* 0D_NOT_real in */,
    const char *axis /* 0D_NOT_character in */,
    bool *include_wall /* 0D_NOT_logical in */,
    const char *gang /* 0D_NOT_character in */,
    bool *skip_lat_layout /* 0D_NOT_logical in */
);
void tao_scale_plot(
    TaoPlotStruct &plot,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis = std::nullopt,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> skip_lat_layout = std::nullopt
);
extern "C" void fortran_tao_scratch_values_calc(
    void *ele_ref /* 0D_PTR_type inout */,
    void *ele_start /* 0D_PTR_type inout */,
    void *ele /* 0D_PTR_type inout */,
    void *datum /* 0D_NOT_type inout */,
    void *branch /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &orbit /* 1D_NOT_type inout */
);
void tao_scratch_values_calc(
    EleStruct &ele_ref,
    EleStruct &ele_start,
    EleStruct &ele,
    TaoDataStruct &datum,
    BranchStruct &branch,
    CoordStructArray1D orbit
);
extern "C" void fortran_tao_set_beam_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    const char *branch_str /* 0D_NOT_character in */
);
void tao_set_beam_cmd(std::string who, std::string value_str, std::string branch_str);
extern "C" void fortran_tao_set_beam_init_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    const char *branch_str /* 0D_NOT_character in */
);
void tao_set_beam_init_cmd(std::string who, std::string value_str, std::string branch_str);
extern "C" void fortran_tao_set_bmad_com_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_bmad_com_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_branch_cmd(
    const char *branch_str /* 0D_NOT_character in */,
    const char *component_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_branch_cmd(std::string branch_str, std::string component_str, std::string value_str);
extern "C" void fortran_tao_set_calculate_cmd(const char *switch_ /* 0D_NOT_character in */);
void tao_set_calculate_cmd(std::optional<std::string> switch_ = std::nullopt);
extern "C" void fortran_tao_set_curve_cmd(
    const char *curve_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_curve_cmd(std::string curve_name, std::string component, std::string value_str);
extern "C" void fortran_tao_set_curve_invalid(
    void *curve /* 0D_NOT_type inout */,
    const char *why_invalid /* 0D_NOT_character in */,
    bool *print_err /* 0D_NOT_logical in */
);
void tao_set_curve_invalid(
    TaoCurveStruct &curve,
    std::string why_invalid,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_tao_set_data_cmd(
    const char *who_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    bool *silent /* 0D_NOT_logical in */
);
void tao_set_data_cmd(
    std::string who_str,
    std::string value_str,
    std::optional<bool> silent = std::nullopt
);
extern "C" void fortran_tao_set_data_useit_opt(Bmad::array_descriptor_t &data /* 1D_NOT_type in */);
void tao_set_data_useit_opt(std::optional<TaoDataStructArray1D> data = std::nullopt);
extern "C" void fortran_tao_set_default_cmd(
    const char *who_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_default_cmd(std::string who_str, std::string value_str);
extern "C" void fortran_tao_set_drawing_cmd(
    void *drawing /* 0D_NOT_type in */,
    const char *component /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_drawing_cmd(TaoDrawingStruct &drawing, std::string component, std::string value_str);
extern "C" void fortran_tao_set_dynamic_aperture_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_dynamic_aperture_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_elements_cmd(
    const char *ele_list /* 0D_NOT_character in */,
    const char *attribute /* 0D_NOT_character in */,
    const char *value /* 0D_NOT_character in */,
    bool &update /* 0D_NOT_logical in */
);
void tao_set_elements_cmd(
    std::string ele_list,
    std::string attribute,
    std::string value,
    bool update
);

// Skipped unusable routine tao_set_flags_for_changed_attribute:
// - Untranslated type: all_pointer_struct (0D)
extern "C" void fortran_tao_set_floor_plan_axis_label(
    void *graph /* 0D_NOT_type inout */,
    void *axis_in /* 0D_NOT_type inout */,
    void *axis_out /* 0D_NOT_type inout */,
    const char *which /* 0D_NOT_character in */
);
void tao_set_floor_plan_axis_label(
    TaoGraphStruct &graph,
    QpAxisStruct &axis_in,
    QpAxisStruct &axis_out,
    std::string which
);
extern "C" void fortran_tao_set_geodesic_lm_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_geodesic_lm_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_global_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_global_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_graph_cmd(
    const char *graph_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_graph_cmd(std::string graph_name, std::string component, std::string value_str);
extern "C" void fortran_tao_set_integer_value(
    int &var /* 0D_NOT_integer out */,
    const char *var_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */,
    int *min_val /* 0D_NOT_integer in */,
    int *max_val /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */
);
struct TaoSetIntegerValue {
  int var;
  bool error;
};
Tao::TaoSetIntegerValue tao_set_integer_value(
    std::string var_str,
    std::string value_str,
    std::optional<int> min_val = std::nullopt,
    std::optional<int> max_val = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_tao_set_invalid(
    void *datum /* 0D_NOT_type in */,
    const char *message /* 0D_NOT_character in */,
    const char *why_invalid /* 0D_NOT_character out */,
    bool *exterminate /* 0D_NOT_logical in */,
    int *err_level /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */
);
std::string tao_set_invalid(
    TaoDataStruct &datum,
    std::string message,
    std::optional<bool> exterminate = std::nullopt,
    std::optional<int> err_level = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" void fortran_tao_set_key_cmd(
    const char *key_str /* 0D_NOT_character in */,
    const char *cmd_str /* 0D_NOT_character in */
);
void tao_set_key_cmd(std::string key_str, std::string cmd_str);
extern "C" void fortran_tao_set_lattice_cmd(
    const char *dest_lat /* 0D_NOT_character in */,
    const char *source_lat /* 0D_NOT_character in */
);
void tao_set_lattice_cmd(std::string dest_lat, std::string source_lat);
extern "C" void fortran_tao_set_logical_value(
    bool &var /* 0D_NOT_logical out */,
    const char *var_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */
);
struct TaoSetLogicalValue {
  bool var;
  bool error;
};
Tao::TaoSetLogicalValue tao_set_logical_value(std::string var_str, std::string value_str);
extern "C" void fortran_tao_set_openmp_n_threads(int &n_threads /* 0D_NOT_integer in */);
void tao_set_openmp_n_threads(int n_threads);
extern "C" void fortran_tao_set_opt_vars(
    Bmad::array_descriptor_t &var_vec /* 1D_NOT_real in */,
    bool *print_limit_warning /* 0D_NOT_logical in */
);
void tao_set_opt_vars(
    FArray1D<Real> &var_vec,
    std::optional<bool> print_limit_warning = std::nullopt
);
extern "C" void fortran_tao_set_opti_de_param_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_opti_de_param_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_particle_start_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_particle_start_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_plot_cmd(
    const char *plot_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_plot_cmd(std::string plot_name, std::string component, std::string value_str);
extern "C" void fortran_tao_set_plot_page_cmd(
    const char *component /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    const char *value_str2 /* 0D_NOT_character in */
);
void tao_set_plot_page_cmd(
    std::string component,
    std::string value_str,
    std::optional<std::string> value_str2 = std::nullopt
);

// Skipped unusable routine tao_set_plotting:
// - Untranslated type: tao_plot_page_input (0D)
extern "C" void fortran_tao_set_ptc_com_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_ptc_com_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_qp_axis_struct(
    const char *qp_axis_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    void *qp_axis /* 0D_NOT_type inout */,
    const char *value /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */,
    int &ix_uni /* 0D_NOT_integer out */
);
struct TaoSetQpAxisStruct {
  bool error;
  int ix_uni;
};
Tao::TaoSetQpAxisStruct tao_set_qp_axis_struct(
    std::string qp_axis_name,
    std::string component,
    QpAxisStruct &qp_axis,
    std::string value
);
extern "C" void fortran_tao_set_qp_point_struct(
    const char *qp_point_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    void *qp_point /* 0D_NOT_type inout */,
    const char *value /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */,
    int &ix_uni /* 0D_NOT_integer out */
);
struct TaoSetQpPointStruct {
  bool error;
  int ix_uni;
};
Tao::TaoSetQpPointStruct tao_set_qp_point_struct(
    std::string qp_point_name,
    std::string component,
    QpPointStruct &qp_point,
    std::string value
);
extern "C" void fortran_tao_set_qp_rect_struct(
    const char *qp_rect_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    void *qp_rect /* 0D_NOT_type inout */,
    const char *value /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */,
    int &ix_uni /* 0D_NOT_integer out */
);
struct TaoSetQpRectStruct {
  bool error;
  int ix_uni;
};
Tao::TaoSetQpRectStruct tao_set_qp_rect_struct(
    std::string qp_rect_name,
    std::string component,
    QpRectStruct &qp_rect,
    std::string value
);
extern "C" void fortran_tao_set_ran_state_cmd(const char *state_string /* 0D_NOT_character in */);
void tao_set_ran_state_cmd(std::string state_string);
extern "C" void fortran_tao_set_real_value(
    double &var /* 0D_NOT_real out */,
    const char *var_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    bool &error /* 0D_NOT_logical out */,
    double *min_val /* 0D_NOT_real in */,
    double *max_val /* 0D_NOT_real in */,
    int *dflt_uni /* 0D_NOT_integer in */
);
struct TaoSetRealValue {
  double var;
  bool error;
};
Tao::TaoSetRealValue tao_set_real_value(
    std::string var_str,
    std::string value_str,
    std::optional<double> min_val = std::nullopt,
    std::optional<double> max_val = std::nullopt,
    std::optional<int> dflt_uni = std::nullopt
);
extern "C" void fortran_tao_set_region_cmd(
    const char *region_name /* 0D_NOT_character in */,
    const char *component /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_region_cmd(std::string region_name, std::string component, std::string value_str);
extern "C" void fortran_tao_set_space_charge_com_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_space_charge_com_cmd(std::string who, std::string value_str);

// Skipped unusable routine tao_set_switch_value:
// - Array bounds handling: "Enum 'L_BOUND' found in bounds 'l_bound' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_set_symbolic_number_cmd(
    const char *sym_str /* 0D_NOT_character in */,
    const char *num_str /* 0D_NOT_character in */,
    double *val /* 0D_NOT_real in */
);
void tao_set_symbolic_number_cmd(
    std::string sym_str,
    std::optional<std::string> num_str = std::nullopt,
    std::optional<double> val = std::nullopt
);
extern "C" void fortran_tao_set_tune_cmd(
    const char *branch_str /* 0D_NOT_character in */,
    const char *mask_str /* 0D_NOT_character in */,
    bool &print_list /* 0D_NOT_logical in */,
    const char *qa_str /* 0D_NOT_character in */,
    const char *qb_str /* 0D_NOT_character in */,
    bool &delta_input /* 0D_NOT_logical in */
);
void tao_set_tune_cmd(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string qa_str,
    std::string qb_str,
    bool delta_input
);
extern "C" void fortran_tao_set_universe_cmd(
    const char *uni /* 0D_NOT_character in */,
    const char *who /* 0D_NOT_character in */,
    const char *what /* 0D_NOT_character in */
);
void tao_set_universe_cmd(std::string uni, std::string who, std::string what);
extern "C" void fortran_tao_set_var_cmd(
    const char *var_str /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */
);
void tao_set_var_cmd(std::string var_str, std::string value_str);
extern "C" void fortran_tao_set_var_model_value(
    void *var /* 0D_NOT_type in */,
    double &value /* 0D_NOT_real in */,
    bool *print_limit_warning /* 0D_NOT_logical in */
);
void tao_set_var_model_value(
    TaoVarStruct &var,
    double value,
    std::optional<bool> print_limit_warning = std::nullopt
);
extern "C" void fortran_tao_set_var_useit_opt();
void tao_set_var_useit_opt();
extern "C" void fortran_tao_set_wave_cmd(
    const char *who /* 0D_NOT_character in */,
    const char *value_str /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical out */
);
bool tao_set_wave_cmd(std::string who, std::string value_str);
extern "C" void fortran_tao_set_z_tune_cmd(
    const char *branch_str /* 0D_NOT_character in */,
    const char *q_str /* 0D_NOT_character in */,
    bool &delta_input /* 0D_NOT_logical in */
);
void tao_set_z_tune_cmd(std::string branch_str, std::string q_str, bool delta_input);
extern "C" void fortran_tao_setup_key_table();
void tao_setup_key_table();
extern "C" void fortran_tao_shape_init(
    void *shape /* 0D_NOT_type inout */,
    bool &err /* 0D_NOT_logical out */,
    bool *print_err /* 0D_NOT_logical in */
);
bool tao_shape_init(TaoEleShapeStruct &shape, std::optional<bool> print_err = std::nullopt);
extern "C" void fortran_tao_show_cmd(const char *what /* 0D_NOT_character in */);
void tao_show_cmd(std::string what);
extern "C" void fortran_tao_show_constraints(
    int &iunit /* 0D_NOT_integer in */,
    const char *form /* 0D_NOT_character in */
);
void tao_show_constraints(int iunit, std::string form);

// Skipped unusable routine tao_show_this:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_single_mode(const char *char_ /* 0D_NOT_character in */);
void tao_single_mode(std::string char_);
extern "C" void fortran_tao_single_track(
    void *tao_lat /* 0D_NOT_type in */,
    bool &calc_ok /* 0D_NOT_logical out */,
    int &ix_branch /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */
);
bool tao_single_track(
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    std::optional<bool> print_err = std::nullopt
);
extern "C" bool fortran_tao_spin_matrices_calc_needed(
    const char *data_type /* 0D_NOT_character in */,
    const char *data_source /* 0D_NOT_character in */,
    bool &do_calc /* 0D_NOT_logical in */
);
void tao_spin_matrices_calc_needed(std::string data_type, std::string data_source, bool do_calc);

// Skipped unusable routine tao_spin_matrix_calc:
// - Routine in configuration skip list

// Skipped unusable routine tao_spin_polarization_calc:
// - Routine in configuration skip list
extern "C" void fortran_tao_spin_tracking_turn_on();
void tao_spin_tracking_turn_on();
extern "C" void fortran_tao_split_component(
    const char *comp_str /* 0D_NOT_character in */,
    void *comp /* 1D_ALLOC_type inout */,
    bool &err /* 0D_NOT_logical out */
);
bool tao_split_component(std::string comp_str, TaoDataVarComponentStructAlloc1D comp);
extern "C" bool fortran_tao_srdt_calc_needed(
    const char *data_type /* 0D_NOT_character in */,
    const char *data_source /* 0D_NOT_character in */,
    int &do_srdt /* 0D_NOT_integer in */
);
void tao_srdt_calc_needed(std::string data_type, std::string data_source, int do_srdt);
extern "C" bool fortran_tao_subin_uni_number(
    const char *name_in /* 0D_NOT_character in */,
    int &ix_uni /* 0D_NOT_integer in */,
    const char *name_out /* 0D_NOT_character out */,
    bool &ok /* 0D_NOT_logical out */
);
struct TaoSubinUniNumber {
  std::string name_out;
  bool ok;
};
Tao::TaoSubinUniNumber tao_subin_uni_number(std::string name_in, int ix_uni);

// Skipped unusable routine tao_svd_func:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_tao_svd_optimizer(bool &abort /* 0D_NOT_logical out */);
bool tao_svd_optimizer();
extern "C" void fortran_tao_symbol_import_from_lat(void *lat /* 0D_NOT_type inout */);
void tao_symbol_import_from_lat(LatStruct &lat);
extern "C" void fortran_tao_taper_cmd(
    const char *except /* 0D_NOT_character in */,
    const char *uni_names /* 0D_NOT_character in */
);
void tao_taper_cmd(std::string except, std::string uni_names);

// Skipped unusable routine tao_timer:
// - Module name unset
extern "C" void fortran_tao_to_change_number(
    const char *num_str /* 0D_NOT_character in */,
    int &n_size /* 0D_NOT_integer in */,
    void *change_number /* 1D_ALLOC_real inout */,
    const char *abs_or_rel /* 0D_NOT_character in */,
    bool &err /* 0D_NOT_logical in */
);
void tao_to_change_number(
    std::string num_str,
    int n_size,
    RealAlloc1D &change_number,
    std::string abs_or_rel,
    bool err
);
extern "C" void fortran_tao_to_int(
    const char *str /* 0D_NOT_character in */,
    int &i_int /* 0D_NOT_integer in */,
    bool &err /* 0D_NOT_logical in */
);
void tao_to_int(std::string str, int i_int, bool err);
extern "C" void fortran_tao_to_phase_and_coupling_reading(
    void *ele /* 0D_NOT_type in */,
    void *bpm_data /* 0D_NOT_type out */,
    bool &valid_value /* 0D_NOT_logical out */,
    const char *why_invalid /* 0D_NOT_character in */,
    void *datum /* 0D_NOT_type inout */
);
struct TaoToPhaseAndCouplingReading {
  BpmPhaseCouplingStruct bpm_data;
  bool valid_value;
};
Tao::TaoToPhaseAndCouplingReading
tao_to_phase_and_coupling_reading(EleStruct &ele, std::string why_invalid, TaoDataStruct &datum);
extern "C" void fortran_tao_to_real(
    const char *expression /* 0D_NOT_character in */,
    double &value /* 0D_NOT_real out */,
    bool &err_flag /* 0D_NOT_logical out */
);
struct TaoToReal {
  double value;
  bool err_flag;
};
Tao::TaoToReal tao_to_real(std::string expression);

// Skipped unusable routine tao_to_top10:
// - Untranslated type: tao_top10_struct (1D)
extern "C" bool fortran_tao_too_many_particles_lost(
    void *beam /* 0D_NOT_type inout */,
    bool &no_beam /* 0D_NOT_logical in */
);
void tao_too_many_particles_lost(BeamStruct &beam, bool no_beam);
extern "C" void fortran_tao_top10_derivative_print();
void tao_top10_derivative_print();
extern "C" void fortran_tao_top10_merit_categories_print(int &iunit /* 0D_NOT_integer in */);
void tao_top10_merit_categories_print(int iunit);
extern "C" void fortran_tao_top_level(
    const char *command /* 0D_NOT_character in */,
    int &errcode /* 0D_NOT_integer out */
);
int tao_top_level(std::optional<std::string> command = std::nullopt);
extern "C" bool fortran_tao_tracking_ele_index(
    void *ele /* 0D_PTR_type in */,
    void *datum /* 0D_NOT_type in */,
    int &ix_branch /* 0D_NOT_integer out */,
    int &ix_ele /* 0D_NOT_integer out */
);
struct TaoTrackingEleIndex {
  int ix_branch;
  int ix_ele;
};
Tao::TaoTrackingEleIndex tao_tracking_ele_index(EleStruct &ele, TaoDataStruct &datum);
extern "C" void fortran_tao_turn_on_special_calcs_if_needed_for_plotting();
void tao_turn_on_special_calcs_if_needed_for_plotting();
extern "C" void fortran_tao_type_expression_tree(
    void *tree /* 0D_NOT_type in */,
    int *indent /* 0D_NOT_integer in */
);
void tao_type_expression_tree(TaoEvalNodeStruct &tree, std::optional<int> indent = std::nullopt);
extern "C" bool fortran_tao_uni_atsign_index(
    const char *string /* 0D_NOT_character in */,
    int &ix_amp /* 0D_NOT_integer out */
);
int tao_uni_atsign_index(std::string string);
extern "C" bool fortran_tao_universe_index(
    int &i_uni /* 0D_NOT_integer in */,
    bool *neg2_to_default /* 0D_NOT_logical in */,
    int &i_this_uni /* 0D_NOT_integer out */
);
int tao_universe_index(int i_uni, std::optional<bool> neg2_to_default = std::nullopt);
extern "C" void fortran_tao_use_data(
    const char *action /* 0D_NOT_character in */,
    const char *data_name /* 0D_NOT_character in */
);
void tao_use_data(std::string action, std::string data_name);
extern "C" void fortran_tao_use_var(
    const char *action /* 0D_NOT_character in */,
    const char *var_name /* 0D_NOT_character in */
);
void tao_use_var(std::string action, std::string var_name);
extern "C" bool
fortran_tao_user_is_terminating_optimization(bool &is_terminating /* 0D_NOT_logical out */);
bool tao_user_is_terminating_optimization();
extern "C" bool fortran_tao_var1_name(
    void *var /* 0D_NOT_type in */,
    const char *var1_name /* 0D_NOT_character out */
);
std::string tao_var1_name(TaoVarStruct &var);
extern "C" bool fortran_tao_var_attrib_name(
    void *var /* 0D_NOT_type in */,
    const char *var_attrib_name /* 0D_NOT_character out */
);
std::string tao_var_attrib_name(TaoVarStruct &var);
extern "C" void fortran_tao_var_check(
    void *eles /* 1D_ALLOC_type in */,
    const char *attribute /* 0D_NOT_character in */,
    bool &silent /* 0D_NOT_logical in */
);
void tao_var_check(ElePointerStructAlloc1D eles, std::string attribute, bool silent);
extern "C" void fortran_tao_var_repoint();
void tao_var_repoint();

// Skipped unusable routine tao_var_show_use:
// - Variable-sized inout character array: 1D_ALLOC_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_var_stuffit1:
// - Untranslated type: tao_var_input (1D)
// - Array bounds handling: "Enum 'N_VAR_MINN' found in bounds 'n_var_minn' but not in provided
// map."
// - Untranslated type: tao_v1_var_input (0D)
// - 2D array handling not supported for logical: dflt_good_unis(lbound(s%u,1):) 2D_ALLOC_logical
// - Array bounds handling: "Enum '1)' found in bounds '1)' but not in provided map."
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tao_var_stuffit2:
// - 2D array handling not supported for logical: good_unis(lbound(s%u,1):) 2D_ALLOC_logical
// - Array bounds handling: "Enum '1)' found in bounds '1)' but not in provided map."
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_var_target_calc();
void tao_var_target_calc();
extern "C" void fortran_tao_var_useit_plot_calc(
    void *graph /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &var /* 1D_NOT_type inout */
);
void tao_var_useit_plot_calc(TaoGraphStruct &graph, TaoVarStructArray1D var);
extern "C" void fortran_tao_var_write(
    const char *out_file /* 0D_NOT_character in */,
    bool *show_good_opt_only /* 0D_NOT_logical in */,
    bool *tao_format /* 0D_NOT_logical in */
);
void tao_var_write(
    std::string out_file,
    std::optional<bool> show_good_opt_only = std::nullopt,
    std::optional<bool> tao_format = std::nullopt
);
extern "C" void fortran_tao_veto_vars_with_zero_dmodel();
void tao_veto_vars_with_zero_dmodel();
extern "C" void fortran_tao_wave_analysis(void *plot /* 0D_NOT_type inout */);
void tao_wave_analysis(TaoPlotStruct &plot);
extern "C" void fortran_tao_wave_cmd(
    const char *curve_name /* 0D_NOT_character in */,
    const char *plot_place /* 0D_NOT_character in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void tao_wave_cmd(std::string curve_name, std::string plot_place, bool err_flag);
extern "C" void fortran_tao_wave_fit(
    void *curve /* 0D_NOT_type in */,
    int &ix1 /* 0D_NOT_integer in */,
    int &n_dat /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &coef /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &rms /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &f1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &f2 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &f3 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &f4 /* 1D_NOT_real in */
);
void tao_wave_fit(
    TaoCurveStruct &curve,
    int ix1,
    int n_dat,
    FArray1D<Real> &coef,
    FArray1D<Real> &rms,
    FArray1D<Real> &f1,
    std::optional<FArray1D<Real>> f2 = std::nullopt,
    std::optional<FArray1D<Real>> f3 = std::nullopt,
    std::optional<FArray1D<Real>> f4 = std::nullopt
);
extern "C" void fortran_tao_write_cmd(const char *what /* 0D_NOT_character in */);
void tao_write_cmd(std::string what);

// Skipped unusable routine tao_write_lines:
// - Variable-sized in character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_tao_x_axis_cmd(
    const char *where /* 0D_NOT_character in */,
    const char *what /* 0D_NOT_character in */
);
void tao_x_axis_cmd(std::string where, std::string what);
extern "C" void fortran_tao_x_scale_cmd(
    const char *where /* 0D_NOT_character in */,
    double &x_min_in /* 0D_NOT_real in */,
    double &x_max_in /* 0D_NOT_real in */,
    bool &err /* 0D_NOT_logical out */,
    bool *include_wall /* 0D_NOT_logical in */,
    const char *gang /* 0D_NOT_character in */,
    bool *exact /* 0D_NOT_logical in */,
    bool *turn_autoscale_off /* 0D_NOT_logical in */
);
bool tao_x_scale_cmd(
    std::string where,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt,
    std::optional<bool> exact = std::nullopt,
    std::optional<bool> turn_autoscale_off = std::nullopt
);
extern "C" void fortran_tao_x_scale_graph(
    void *graph /* 0D_NOT_type inout */,
    double &x_min /* 0D_NOT_real in */,
    double &x_max /* 0D_NOT_real in */,
    bool *include_wall /* 0D_NOT_logical in */,
    bool *have_scaled /* 0D_NOT_logical in */
);
void tao_x_scale_graph(
    TaoGraphStruct &graph,
    double x_min,
    double x_max,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<bool> have_scaled = std::nullopt
);
extern "C" void fortran_tao_x_scale_plot(
    void *plot /* 0D_NOT_type in */,
    double &x_min_in /* 0D_NOT_real in */,
    double &x_max_in /* 0D_NOT_real in */,
    bool *include_wall /* 0D_NOT_logical in */,
    const char *gang /* 0D_NOT_character in */,
    bool &have_scaled /* 0D_NOT_logical out */
);
bool tao_x_scale_plot(
    TaoPlotStruct &plot,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall = std::nullopt,
    std::optional<std::string> gang = std::nullopt
);

// Skipped unusable routine user_signal:
// - Module name unset
} // namespace Tao
