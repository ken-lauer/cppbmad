#include "bmad/generated/tao_routines.hpp"

#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/proxy.hpp"
#include "bmad/json.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

using json = nlohmann::json;
void Tao::integrate_max(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchStruct &branch,
    FArray1D<Real> &vec,
    TaoDataStruct &datum
) {
  // vec: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  fortran_integrate_max(
      /* int& */ ix_start,
      /* int& */ ix_ele,
      /* double& */ datum_value,
      /* int& */ ix_m,
      /* void* */ branch.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* void* */ datum.get_fortran_ptr()
  );
}
void Tao::integrate_min(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchStruct &branch,
    FArray1D<Real> &vec,
    TaoDataStruct &datum
) {
  // vec: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  fortran_integrate_min(
      /* int& */ ix_start,
      /* int& */ ix_ele,
      /* double& */ datum_value,
      /* int& */ ix_m,
      /* void* */ branch.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* void* */ datum.get_fortran_ptr()
  );
}
void Tao::tao_abort_command_file(std::optional<bool> force_abort) {
  bool force_abort_lvalue;
  auto *_force_abort{&force_abort_lvalue};
  if (force_abort.has_value()) {
    force_abort_lvalue = force_abort.value();
  } else {
    _force_abort = nullptr;
  }
  fortran_tao_abort_command_file(/* bool* */ _force_abort);
}
void Tao::tao_add_to_normal_mode_h_array(std::string h_str, ResonanceHStructAlloc1D h_array) {
  auto _h_str = h_str.c_str();
  // intent=inout allocatable type array
  fortran_tao_add_to_normal_mode_h_array(
      /* const char* */ _h_str,
      /* void* */ h_array.get_fortran_ptr()
  );
}
void Tao::tao_alias_cmd(std::string alias, std::string string) {
  auto _alias = alias.c_str();
  auto _string = string.c_str();
  fortran_tao_alias_cmd(/* const char* */ _alias, /* const char* */ _string);
}
void Tao::tao_allocate_data_array(TaoUniverseStruct &u, int n_data, std::optional<bool> exact) {
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_tao_allocate_data_array(
      /* void* */ u.get_fortran_ptr(),
      /* int& */ n_data,
      /* bool* */ _exact
  );
}
void Tao::tao_allocate_v1_var(int n_v1, bool save_old) {
  fortran_tao_allocate_v1_var(/* int& */ n_v1, /* bool& */ save_old);
}
void Tao::tao_allocate_var_array(int n_var, bool default_good_user) {
  fortran_tao_allocate_var_array(/* int& */ n_var, /* bool& */ default_good_user);
}
double
Tao::tao_beam_emit_calc(int plane, int emit_type, EleStruct &ele, BunchParamsStruct &bunch_params) {
  double _emit{};
  fortran_tao_beam_emit_calc(
      /* int& */ plane,
      /* int& */ emit_type,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* double& */ _emit
  );
  return _emit;
}
bool Tao::tao_beam_track(
    TaoUniverseStruct &u,
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    BeamStruct &beam
) {
  bool _calc_ok{};
  fortran_tao_beam_track(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* int& */ ix_branch,
      /* void* */ beam.get_fortran_ptr(),
      /* bool& */ _calc_ok
  );
  return _calc_ok;
}
std::optional<EleStruct> Tao::tao_beam_track_endpoint(
    std::string ele_id,
    LatStruct &lat,
    std::string branch_str,
    std::string where,
    TaoUniverseStruct &u
) {
  auto _ele_id = ele_id.c_str();
  auto _branch_str = branch_str.c_str();
  auto _where = where.c_str();
  void *_ele;
  fortran_tao_beam_track_endpoint(
      /* const char* */ _ele_id,
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _branch_str,
      /* const char* */ _where,
      /* void* */ u.get_fortran_ptr(),
      /* void* */ &_ele
  );
  return std::move((_ele ? std::make_optional<EleStruct>(_ele) : std::nullopt));
}
int Tao::tao_branch_index(int ix_branch) {
  int _ix_this{};
  fortran_tao_branch_index(/* int& */ ix_branch, /* int& */ _ix_this);
  return _ix_this;
}
void Tao::tao_calc_data_at_s_pts(
    TaoLatticeStruct &tao_lat,
    TaoCurveStruct &curve,
    double comp_sign,
    BoolAlloc1D &good
) {
  // intent=inout allocatable general array
  fortran_tao_calc_data_at_s_pts(
      /* void* */ tao_lat.get_fortran_ptr(),
      /* void* */ curve.get_fortran_ptr(),
      /* double& */ comp_sign,
      /* void* */ good.get_fortran_ptr()
  );
}
void Tao::tao_cbar_wave_anal(TaoPlotStruct &plot) {
  fortran_tao_cbar_wave_anal(/* void* */ plot.get_fortran_ptr());
}
bool Tao::tao_change_ele(
    std::string ele_name,
    std::string attrib_name,
    std::string num_str,
    bool update
) {
  auto _ele_name = ele_name.c_str();
  auto _attrib_name = attrib_name.c_str();
  auto _num_str = num_str.c_str();
  bool _err_flag{};
  fortran_tao_change_ele(
      /* const char* */ _ele_name,
      /* const char* */ _attrib_name,
      /* const char* */ _num_str,
      /* bool& */ update,
      /* bool& */ _err_flag
  );
  return _err_flag;
}
bool Tao::tao_change_tune(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string dqa_str,
    std::string dqb_str
) {
  auto _branch_str = branch_str.c_str();
  auto _mask_str = mask_str.c_str();
  auto _dqa_str = dqa_str.c_str();
  auto _dqb_str = dqb_str.c_str();
  bool _err_flag{};
  fortran_tao_change_tune(
      /* const char* */ _branch_str,
      /* const char* */ _mask_str,
      /* bool& */ print_list,
      /* const char* */ _dqa_str,
      /* const char* */ _dqb_str,
      /* bool& */ _err_flag
  );
  return _err_flag;
}
bool Tao::tao_change_var(std::string name, std::string num_str, bool silent) {
  auto _name = name.c_str();
  auto _num_str = num_str.c_str();
  bool _err_flag{};
  fortran_tao_change_var(
      /* const char* */ _name,
      /* const char* */ _num_str,
      /* bool& */ silent,
      /* bool& */ _err_flag
  );
  return _err_flag;
}
bool Tao::tao_change_z_tune(std::string branch_str, std::string dq_str) {
  auto _branch_str = branch_str.c_str();
  auto _dq_str = dq_str.c_str();
  bool _err_flag{};
  fortran_tao_change_z_tune(
      /* const char* */ _branch_str,
      /* const char* */ _dq_str,
      /* bool& */ _err_flag
  );
  return _err_flag;
}
void Tao::tao_chrom_calc_needed(std::string data_type, std::string data_source, bool do_chrom) {
  auto _data_type = data_type.c_str();
  auto _data_source = data_source.c_str();
  fortran_tao_chrom_calc_needed(
      /* const char* */ _data_type,
      /* const char* */ _data_source,
      /* bool& */ do_chrom
  );
}
void Tao::tao_clear_cmd(std::string cmd_line) {
  auto _cmd_line = cmd_line.c_str();
  fortran_tao_clear_cmd(/* const char* */ _cmd_line);
}
void Tao::tao_clip_cmd(bool gang, std::string where, double value1, double value2) {
  auto _where = where.c_str();
  fortran_tao_clip_cmd(
      /* bool& */ gang,
      /* const char* */ _where,
      /* double& */ value1,
      /* double& */ value2
  );
}
void Tao::tao_close_command_file() { fortran_tao_close_command_file(); }
void Tao::tao_cmd_history_record(std::string cmd) {
  auto _cmd = cmd.c_str();
  fortran_tao_cmd_history_record(/* const char* */ _cmd);
}
bool Tao::tao_command(std::string command_line, bool err) {
  auto _command_line = command_line.c_str();
  bool _err_is_fatal{};
  fortran_tao_command(/* const char* */ _command_line, /* bool& */ err, /* bool& */ _err_is_fatal);
  return _err_is_fatal;
}
std::string Tao::tao_constraint_type_name(TaoDataStruct &datum) {
  char _datum_name[4096];
  fortran_tao_constraint_type_name(
      /* void* */ datum.get_fortran_ptr(),
      /* const char* */ _datum_name
  );
  return _datum_name;
}
void Tao::tao_control_tree_list(EleStruct &ele, ElePointerStructAlloc1D tree) {
  // intent=in allocatable type array
  fortran_tao_control_tree_list(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ tree.get_fortran_ptr()
  );
}
int Tao::tao_count_strings(std::string string, std::string pattern) {
  auto _string = string.c_str();
  auto _pattern = pattern.c_str();
  int _num{};
  fortran_tao_count_strings(/* const char* */ _string, /* const char* */ _pattern, /* int& */ _num);
  return _num;
}
void Tao::tao_create_plot_window() { fortran_tao_create_plot_window(); }
void Tao::tao_curve_beam_ellipse_setup(TaoCurveStruct &curve) {
  fortran_tao_curve_beam_ellipse_setup(/* void* */ curve.get_fortran_ptr());
}
bool Tao::tao_curve_check_universe(TaoCurveStruct &curve, TaoUniverseStruct &uni) {
  auto _uni = &uni; // input, required, pointer
  bool _is_ok{};
  fortran_tao_curve_check_universe(
      /* void* */ curve.get_fortran_ptr(),
      /* void* */ &uni,
      /* bool& */ _is_ok
  );
  return _is_ok;
}
void Tao::tao_curve_data_setup(TaoPlotStruct &plot, TaoGraphStruct &graph, TaoCurveStruct &curve) {
  fortran_tao_curve_data_setup(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ curve.get_fortran_ptr()
  );
}
void Tao::tao_curve_datum_calc(
    ElePointerStructAlloc1D eles,
    TaoPlotStruct &plot,
    TaoCurveStruct &curve,
    std::string who
) {
  // intent=in allocatable type array
  auto _who = who.c_str();
  fortran_tao_curve_datum_calc(
      /* void* */ eles.get_fortran_ptr(),
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ curve.get_fortran_ptr(),
      /* const char* */ _who
  );
}
void Tao::tao_curve_ele_ref(TaoCurveStruct &curve, bool point_to_ele_ref, EleStruct &ele_track) {
  auto _ele_track = &ele_track; // input, required, pointer
  fortran_tao_curve_ele_ref(
      /* void* */ curve.get_fortran_ptr(),
      /* bool& */ point_to_ele_ref,
      /* void* */ &ele_track
  );
}
int Tao::tao_curve_ix_uni(TaoCurveStruct &curve) {
  int _ix_uni{};
  fortran_tao_curve_ix_uni(/* void* */ curve.get_fortran_ptr(), /* int& */ _ix_uni);
  return _ix_uni;
}
std::string Tao::tao_curve_name(TaoCurveStruct &curve, std::optional<bool> use_region) {
  bool use_region_lvalue;
  auto *_use_region{&use_region_lvalue};
  if (use_region.has_value()) {
    use_region_lvalue = use_region.value();
  } else {
    _use_region = nullptr;
  }
  char _curve_name[4096];
  fortran_tao_curve_name(
      /* void* */ curve.get_fortran_ptr(),
      /* bool* */ _use_region,
      /* const char* */ _curve_name
  );
  return _curve_name;
}
Tao::TaoCurveRmsCalc Tao::tao_curve_rms_calc(TaoCurveStruct &curve, std::string who) {
  auto _who = who.c_str();
  double _rms{};
  double _mean{};
  fortran_tao_curve_rms_calc(
      /* void* */ curve.get_fortran_ptr(),
      /* const char* */ _who,
      /* double& */ _rms,
      /* double& */ _mean
  );
  return TaoCurveRmsCalc{_rms, _mean};
}
std::string Tao::tao_d2_d1_name(TaoD1DataStruct &d1, std::optional<bool> show_universe) {
  bool show_universe_lvalue;
  auto *_show_universe{&show_universe_lvalue};
  if (show_universe.has_value()) {
    show_universe_lvalue = show_universe.value();
  } else {
    _show_universe = nullptr;
  }
  char _d2_d1_name[4096];
  fortran_tao_d2_d1_name(
      /* void* */ d1.get_fortran_ptr(),
      /* bool* */ _show_universe,
      /* const char* */ _d2_d1_name
  );
  return _d2_d1_name;
}
void Tao::tao_d2_data_stuffit(TaoUniverseStruct &u, std::string d2_name, int n_d1_data) {
  auto _d2_name = d2_name.c_str();
  fortran_tao_d2_data_stuffit(
      /* void* */ u.get_fortran_ptr(),
      /* const char* */ _d2_name,
      /* int& */ n_d1_data
  );
}
void Tao::tao_data_check(bool err) { fortran_tao_data_check(/* bool& */ err); }
void Tao::tao_data_coupling_init(BranchStruct &branch) {
  fortran_tao_data_coupling_init(/* void* */ branch.get_fortran_ptr());
}
bool Tao::tao_data_sanity_check(
    TaoDataStruct &datum,
    bool print_err,
    std::string default_data_type,
    optional_ref<TaoUniverseStruct> uni
) {
  auto _default_data_type = default_data_type.c_str();
  auto *_uni = uni.has_value() ? uni->get().get_fortran_ptr() : nullptr; // input, optional
  bool _is_valid{};
  fortran_tao_data_sanity_check(
      /* void* */ datum.get_fortran_ptr(),
      /* bool& */ print_err,
      /* const char* */ _default_data_type,
      /* void* */ _uni,
      /* bool& */ _is_valid
  );
  return _is_valid;
}
std::string
Tao::tao_data_type_substitute(std::string template_, TaoCurveStruct &curve, TaoGraphStruct &graph) {
  auto _template_ = template_.c_str();
  char _str_out[4096];
  fortran_tao_data_type_substitute(
      /* const char* */ _template_,
      /* const char* */ _str_out,
      /* void* */ curve.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
  return _str_out;
}
std::string Tao::tao_data_useit_plot_calc(
    TaoCurveStruct &curve,
    TaoGraphStruct &graph,
    TaoDataStructArray1D data,
    bool check_s_position
) {
  // data: TaoDataStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _data_desc;
  _data_desc.rank = 1;
  _data_desc.data_ptr = data.data();
  _data_desc.dims[0] = data.size();
  _data_desc.strides[0] = 1;
  char _most_invalid[4096];
  fortran_tao_data_useit_plot_calc(
      /* void* */ curve.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _data_desc,
      /* bool& */ check_s_position,
      /* const char* */ _most_invalid
  );
  return _most_invalid;
}
int Tao::tao_datum_has_associated_ele(std::string data_type, std::optional<int> branch_geometry) {
  auto _data_type = data_type.c_str();
  int branch_geometry_lvalue;
  auto *_branch_geometry{&branch_geometry_lvalue};
  if (branch_geometry.has_value()) {
    branch_geometry_lvalue = branch_geometry.value();
  } else {
    _branch_geometry = nullptr;
  }
  int _has_associated_ele{};
  fortran_tao_datum_has_associated_ele(
      /* const char* */ _data_type,
      /* int* */ _branch_geometry,
      /* int& */ _has_associated_ele
  );
  return _has_associated_ele;
}
Tao::TaoDatumIntegrate Tao::tao_datum_integrate(
    TaoDataStruct &datum,
    BranchStruct &branch,
    FArray1D<Real> &s_pos,
    FArray1D<Real> &values
) {
  // s_pos: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _s_pos_desc;
  _s_pos_desc.rank = 1;
  _s_pos_desc.data_ptr = s_pos.data();
  _s_pos_desc.dims[0] = s_pos.size();
  // values: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _values_desc;
  _values_desc.rank = 1;
  _values_desc.data_ptr = values.data();
  _values_desc.dims[0] = values.size();
  bool _valid_value{};
  char _why_invalid[4096];
  double _result{};
  fortran_tao_datum_integrate(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _s_pos_desc,
      /* Bmad::array_descriptor_t& */ _values_desc,
      /* bool& */ _valid_value,
      /* const char* */ _why_invalid,
      /* double& */ _result
  );
  return TaoDatumIntegrate{_valid_value, _why_invalid, _result};
}
std::string Tao::tao_datum_name(TaoDataStruct &datum, std::optional<bool> show_universe) {
  bool show_universe_lvalue;
  auto *_show_universe{&show_universe_lvalue};
  if (show_universe.has_value()) {
    show_universe_lvalue = show_universe.value();
  } else {
    _show_universe = nullptr;
  }
  char _datum_name[4096];
  fortran_tao_datum_name(
      /* void* */ datum.get_fortran_ptr(),
      /* bool* */ _show_universe,
      /* const char* */ _datum_name
  );
  return _datum_name;
}
double Tao::tao_datum_s_position(TaoDataStruct &datum, EleStruct &ele) {
  double _s_pos{};
  fortran_tao_datum_s_position(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ _s_pos
  );
  return _s_pos;
}
bool Tao::tao_de_optimizer() {
  bool _abort{};
  fortran_tao_de_optimizer(/* bool& */ _abort);
  return _abort;
}
void Tao::tao_deallocate_plot_cache(TaoPlotCacheStructAlloc1D plot_cache) {
  // intent=inout allocatable type array
  fortran_tao_deallocate_plot_cache(/* void* */ plot_cache.get_fortran_ptr());
}
void Tao::tao_deallocate_tree(TaoEvalNodeStruct &tree) {
  fortran_tao_deallocate_tree(/* void* */ tree.get_fortran_ptr());
}
void Tao::tao_destroy_plot_window() { fortran_tao_destroy_plot_window(); }
void Tao::tao_dmerit_calc() { fortran_tao_dmerit_calc(); }
bool Tao::tao_dmodel_dvar_calc(bool force_calc) {
  bool _err_flag{};
  fortran_tao_dmodel_dvar_calc(/* bool& */ force_calc, /* bool& */ _err_flag);
  return _err_flag;
}
double Tao::tao_do_wire_scan(EleStruct &ele, double theta, BeamStruct &beam) {
  double _moment{};
  fortran_tao_do_wire_scan(
      /* void* */ ele.get_fortran_ptr(),
      /* double& */ theta,
      /* void* */ beam.get_fortran_ptr(),
      /* double& */ _moment
  );
  return _moment;
}
void Tao::tao_draw_beam_chamber_wall(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_draw_beam_chamber_wall(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_draw_curve_data(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoCurveStruct &curve,
    bool have_data
) {
  fortran_tao_draw_curve_data(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ curve.get_fortran_ptr(),
      /* bool& */ have_data
  );
}
void Tao::tao_draw_ele_for_floor_plan(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoLatticeStruct &tao_lat,
    EleStruct &ele,
    TaoEleShapeStruct &ele_shape,
    std::string label_name,
    double offset1,
    double offset2
) {
  auto _ele_shape = &ele_shape; // input, required, pointer
  auto _label_name = label_name.c_str();
  fortran_tao_draw_ele_for_floor_plan(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ &ele_shape,
      /* const char* */ _label_name,
      /* double& */ offset1,
      /* double& */ offset2
  );
}
void Tao::tao_draw_floor_plan(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_draw_floor_plan(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_draw_graph_axes(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_draw_graph_axes(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_draw_histogram_data(
    TaoPlotStruct &plot,
    TaoGraphStruct &graph,
    TaoCurveStruct &curve,
    bool have_data
) {
  fortran_tao_draw_histogram_data(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ curve.get_fortran_ptr(),
      /* bool& */ have_data
  );
}
void Tao::tao_draw_lat_layout(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_draw_lat_layout(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_draw_plots(std::optional<bool> do_clear) {
  bool do_clear_lvalue;
  auto *_do_clear{&do_clear_lvalue};
  if (do_clear.has_value()) {
    do_clear_lvalue = do_clear.value();
  } else {
    _do_clear = nullptr;
  }
  fortran_tao_draw_plots(/* bool* */ _do_clear);
}
Tao::TaoEleGeometryWithMisalignments
Tao::tao_ele_geometry_with_misalignments(TaoDataStruct &datum, EleStruct &ele) {
  bool _valid_value{};
  char _why_invalid[4096];
  double _value{};
  fortran_tao_ele_geometry_with_misalignments(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _valid_value,
      /* const char* */ _why_invalid,
      /* double& */ _value
  );
  return TaoEleGeometryWithMisalignments{_valid_value, _why_invalid, _value};
}
Tao::TaoEleShapeInfo Tao::tao_ele_shape_info(
    int ix_uni,
    EleStruct &ele,
    TaoEleShapeStructArray1D ele_shapes,
    double y1,
    double y2,
    std::optional<int> ix_shape_min
) {
  // ele_shapes: TaoEleShapeStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ele_shapes_desc;
  _ele_shapes_desc.rank = 1;
  _ele_shapes_desc.data_ptr = ele_shapes.data();
  _ele_shapes_desc.dims[0] = ele_shapes.size();
  _ele_shapes_desc.strides[0] = 1;
  void *_e_shape;
  char _label_name[4096];
  int ix_shape_min_lvalue;
  auto *_ix_shape_min{&ix_shape_min_lvalue};
  if (ix_shape_min.has_value()) {
    ix_shape_min_lvalue = ix_shape_min.value();
  } else {
    _ix_shape_min = nullptr;
  }
  fortran_tao_ele_shape_info(
      /* int& */ ix_uni,
      /* void* */ ele.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _ele_shapes_desc,
      /* void* */ &_e_shape,
      /* const char* */ _label_name,
      /* double& */ y1,
      /* double& */ y2,
      /* int* */ _ix_shape_min
  );
  return TaoEleShapeInfo{
      std::move((_e_shape ? std::make_optional<TaoEleShapeStruct>(_e_shape) : std::nullopt)),
      _label_name
  };
}
Tao::TaoEvalFloorOrbit Tao::tao_eval_floor_orbit(
    TaoDataStruct &datum,
    EleStruct &ele,
    CoordStruct &orbit,
    BunchParamsStruct &bunch_params
) {
  bool _valid_value{};
  char _why_invalid[4096];
  double _value{};
  fortran_tao_eval_floor_orbit(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* bool& */ _valid_value,
      /* const char* */ _why_invalid,
      /* double& */ _value
  );
  return TaoEvalFloorOrbit{_valid_value, _why_invalid, _value};
}
Tao::TaoEvaluateADatum Tao::tao_evaluate_a_datum(
    TaoDataStruct &datum,
    TaoUniverseStruct &u,
    TaoLatticeStruct &tao_lat,
    std::optional<bool> called_from_lat_calc,
    std::optional<bool> print_err
) {
  double _datum_value{};
  bool _valid_value{};
  char _why_invalid[4096];
  bool called_from_lat_calc_lvalue;
  auto *_called_from_lat_calc{&called_from_lat_calc_lvalue};
  if (called_from_lat_calc.has_value()) {
    called_from_lat_calc_lvalue = called_from_lat_calc.value();
  } else {
    _called_from_lat_calc = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_evaluate_a_datum(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ u.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* double& */ _datum_value,
      /* bool& */ _valid_value,
      /* const char* */ _why_invalid,
      /* bool* */ _called_from_lat_calc,
      /* bool* */ _print_err
  );
  return TaoEvaluateADatum{_datum_value, _valid_value, _why_invalid};
}
Tao::TaoEvaluateDatumAtS Tao::tao_evaluate_datum_at_s(
    TaoDataStruct &datum,
    TaoLatticeStruct &tao_lat,
    EleStruct &ele,
    EleStruct &ele_ref,
    bool valid_value
) {
  auto _ele = &ele; // input, required, pointer
  auto _ele_ref = &ele_ref; // input, required, pointer
  char _err_str[4096];
  bool _bad_datum{};
  double _value{};
  fortran_tao_evaluate_datum_at_s(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* void* */ &ele,
      /* void* */ &ele_ref,
      /* bool& */ valid_value,
      /* const char* */ _err_str,
      /* bool& */ _bad_datum,
      /* double& */ _value
  );
  return TaoEvaluateDatumAtS{_err_str, _bad_datum, _value};
}
bool Tao::tao_evaluate_element_parameters(
    std::string param_name,
    RealAlloc1D &values,
    bool print_err,
    optional_ref<EleStruct> dflt_ele,
    std::string dflt_source,
    std::optional<std::string> dflt_component,
    std::optional<int> dflt_uni,
    std::optional<int> eval_point,
    std::optional<TaoExpressionInfoStructAlloc1D> info
) {
  bool _err{};
  auto _param_name = param_name.c_str();
  // intent=inout allocatable general array
  auto *_dflt_ele =
      dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr() : nullptr; // input, optional
  auto _dflt_source = dflt_source.c_str();
  const char *_dflt_component = dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int eval_point_lvalue;
  auto *_eval_point{&eval_point_lvalue};
  if (eval_point.has_value()) {
    eval_point_lvalue = eval_point.value();
  } else {
    _eval_point = nullptr;
  }
  // intent=inout allocatable type array
  auto *_info = info.has_value() ? info->get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_evaluate_element_parameters(
      /* bool& */ _err,
      /* const char* */ _param_name,
      /* void* */ values.get_fortran_ptr(),
      /* bool& */ print_err,
      /* void* */ _dflt_ele,
      /* const char* */ _dflt_source,
      /* const char* */ _dflt_component,
      /* int* */ _dflt_uni,
      /* int* */ _eval_point,
      /* void* */ _info
  );
  return _err;
}
bool Tao::tao_evaluate_expression(
    std::string expression,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    std::optional<bool> print_err,
    std::optional<TaoExpressionInfoStructAlloc1D> info,
    std::optional<TaoEvalNodeStructAlloc1D> stack,
    std::optional<std::string> dflt_component,
    std::optional<std::string> dflt_source,
    optional_ref<EleStruct> dflt_ele_ref,
    optional_ref<EleStruct> dflt_ele_start,
    optional_ref<EleStruct> dflt_ele,
    std::optional<std::string> dflt_dat_or_var_index,
    std::optional<int> dflt_uni,
    std::optional<int> dflt_eval_point,
    std::optional<double> dflt_s_offset,
    optional_ref<CoordStruct> dflt_orbit,
    optional_ref<TaoDataStruct> datum
) {
  auto _expression = expression.c_str();
  // intent=inout allocatable general array
  bool _err_flag{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  // intent=inout allocatable type array
  auto *_info = info.has_value() ? info->get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable type array
  auto *_stack = stack.has_value() ? stack->get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_component = dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  const char *_dflt_source = dflt_source.has_value() ? dflt_source->c_str() : nullptr;
  auto *_dflt_ele_ref =
      dflt_ele_ref.has_value() ? dflt_ele_ref->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_dflt_ele_start = dflt_ele_start.has_value() ? dflt_ele_start->get().get_fortran_ptr()
                                                     : nullptr; // input, optional
  auto *_dflt_ele =
      dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_dat_or_var_index =
      dflt_dat_or_var_index.has_value() ? dflt_dat_or_var_index->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int dflt_eval_point_lvalue;
  auto *_dflt_eval_point{&dflt_eval_point_lvalue};
  if (dflt_eval_point.has_value()) {
    dflt_eval_point_lvalue = dflt_eval_point.value();
  } else {
    _dflt_eval_point = nullptr;
  }
  double dflt_s_offset_lvalue;
  auto *_dflt_s_offset{&dflt_s_offset_lvalue};
  if (dflt_s_offset.has_value()) {
    dflt_s_offset_lvalue = dflt_s_offset.value();
  } else {
    _dflt_s_offset = nullptr;
  }
  auto *_dflt_orbit =
      dflt_orbit.has_value() ? dflt_orbit->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_datum = datum.has_value() ? datum->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_evaluate_expression(
      /* const char* */ _expression,
      /* int& */ n_size,
      /* bool& */ use_good_user,
      /* void* */ value.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _print_err,
      /* void* */ _info,
      /* void* */ _stack,
      /* const char* */ _dflt_component,
      /* const char* */ _dflt_source,
      /* void* */ _dflt_ele_ref,
      /* void* */ _dflt_ele_start,
      /* void* */ _dflt_ele,
      /* const char* */ _dflt_dat_or_var_index,
      /* int* */ _dflt_uni,
      /* int* */ _dflt_eval_point,
      /* double* */ _dflt_s_offset,
      /* void* */ _dflt_orbit,
      /* void* */ _datum
  );
  return _err_flag;
}
bool Tao::tao_evaluate_expression_new(
    std::string expression,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    std::optional<bool> print_err,
    std::optional<TaoExpressionInfoStructAlloc1D> info,
    std::optional<TaoEvalNodeStructAlloc1D> stack,
    std::optional<std::string> dflt_component,
    std::optional<std::string> dflt_source,
    optional_ref<EleStruct> dflt_ele_ref,
    optional_ref<EleStruct> dflt_ele_start,
    optional_ref<EleStruct> dflt_ele,
    std::optional<std::string> dflt_dat_or_var_index,
    std::optional<int> dflt_uni,
    std::optional<int> dflt_eval_point,
    std::optional<double> dflt_s_offset,
    optional_ref<CoordStruct> dflt_orbit,
    optional_ref<TaoDataStruct> datum
) {
  auto _expression = expression.c_str();
  // intent=inout allocatable general array
  bool _err_flag{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  // intent=inout allocatable type array
  auto *_info = info.has_value() ? info->get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable type array
  auto *_stack = stack.has_value() ? stack->get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_component = dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  const char *_dflt_source = dflt_source.has_value() ? dflt_source->c_str() : nullptr;
  auto *_dflt_ele_ref =
      dflt_ele_ref.has_value() ? dflt_ele_ref->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_dflt_ele_start = dflt_ele_start.has_value() ? dflt_ele_start->get().get_fortran_ptr()
                                                     : nullptr; // input, optional
  auto *_dflt_ele =
      dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_dat_or_var_index =
      dflt_dat_or_var_index.has_value() ? dflt_dat_or_var_index->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int dflt_eval_point_lvalue;
  auto *_dflt_eval_point{&dflt_eval_point_lvalue};
  if (dflt_eval_point.has_value()) {
    dflt_eval_point_lvalue = dflt_eval_point.value();
  } else {
    _dflt_eval_point = nullptr;
  }
  double dflt_s_offset_lvalue;
  auto *_dflt_s_offset{&dflt_s_offset_lvalue};
  if (dflt_s_offset.has_value()) {
    dflt_s_offset_lvalue = dflt_s_offset.value();
  } else {
    _dflt_s_offset = nullptr;
  }
  auto *_dflt_orbit =
      dflt_orbit.has_value() ? dflt_orbit->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_datum = datum.has_value() ? datum->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_evaluate_expression_new(
      /* const char* */ _expression,
      /* int& */ n_size,
      /* bool& */ use_good_user,
      /* void* */ value.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _print_err,
      /* void* */ _info,
      /* void* */ _stack,
      /* const char* */ _dflt_component,
      /* const char* */ _dflt_source,
      /* void* */ _dflt_ele_ref,
      /* void* */ _dflt_ele_start,
      /* void* */ _dflt_ele,
      /* const char* */ _dflt_dat_or_var_index,
      /* int* */ _dflt_uni,
      /* int* */ _dflt_eval_point,
      /* double* */ _dflt_s_offset,
      /* void* */ _dflt_orbit,
      /* void* */ _datum
  );
  return _err_flag;
}
bool Tao::tao_evaluate_expression_old(
    std::string expression,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    std::optional<bool> print_err,
    std::optional<TaoExpressionInfoStructAlloc1D> info,
    std::optional<TaoEvalNodeStructAlloc1D> stack,
    std::optional<std::string> dflt_component,
    std::optional<std::string> dflt_source,
    optional_ref<EleStruct> dflt_ele_ref,
    optional_ref<EleStruct> dflt_ele_start,
    optional_ref<EleStruct> dflt_ele,
    std::optional<std::string> dflt_dat_or_var_index,
    std::optional<int> dflt_uni,
    std::optional<int> dflt_eval_point,
    std::optional<double> dflt_s_offset,
    optional_ref<CoordStruct> dflt_orbit,
    optional_ref<TaoDataStruct> datum
) {
  auto _expression = expression.c_str();
  // intent=inout allocatable general array
  bool _err_flag{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  // intent=inout allocatable type array
  auto *_info = info.has_value() ? info->get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable type array
  auto *_stack = stack.has_value() ? stack->get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_component = dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  const char *_dflt_source = dflt_source.has_value() ? dflt_source->c_str() : nullptr;
  auto *_dflt_ele_ref =
      dflt_ele_ref.has_value() ? dflt_ele_ref->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_dflt_ele_start = dflt_ele_start.has_value() ? dflt_ele_start->get().get_fortran_ptr()
                                                     : nullptr; // input, optional
  auto *_dflt_ele =
      dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_dat_or_var_index =
      dflt_dat_or_var_index.has_value() ? dflt_dat_or_var_index->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int dflt_eval_point_lvalue;
  auto *_dflt_eval_point{&dflt_eval_point_lvalue};
  if (dflt_eval_point.has_value()) {
    dflt_eval_point_lvalue = dflt_eval_point.value();
  } else {
    _dflt_eval_point = nullptr;
  }
  double dflt_s_offset_lvalue;
  auto *_dflt_s_offset{&dflt_s_offset_lvalue};
  if (dflt_s_offset.has_value()) {
    dflt_s_offset_lvalue = dflt_s_offset.value();
  } else {
    _dflt_s_offset = nullptr;
  }
  auto *_dflt_orbit =
      dflt_orbit.has_value() ? dflt_orbit->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_datum = datum.has_value() ? datum->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_evaluate_expression_old(
      /* const char* */ _expression,
      /* int& */ n_size,
      /* bool& */ use_good_user,
      /* void* */ value.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool* */ _print_err,
      /* void* */ _info,
      /* void* */ _stack,
      /* const char* */ _dflt_component,
      /* const char* */ _dflt_source,
      /* void* */ _dflt_ele_ref,
      /* void* */ _dflt_ele_start,
      /* void* */ _dflt_ele,
      /* const char* */ _dflt_dat_or_var_index,
      /* int* */ _dflt_uni,
      /* int* */ _dflt_eval_point,
      /* double* */ _dflt_s_offset,
      /* void* */ _dflt_orbit,
      /* void* */ _datum
  );
  return _err_flag;
}
bool Tao::tao_evaluate_lat_or_beam_data(
    std::string data_name,
    RealAlloc1D &values,
    bool print_err,
    std::string default_source,
    optional_ref<EleStruct> dflt_ele_ref,
    optional_ref<EleStruct> dflt_ele_start,
    optional_ref<EleStruct> dflt_ele,
    std::optional<std::string> dflt_component,
    std::optional<int> dflt_uni,
    std::optional<int> dflt_eval_point,
    std::optional<double> dflt_s_offset
) {
  bool _err{};
  auto _data_name = data_name.c_str();
  // intent=inout allocatable general array
  auto _default_source = default_source.c_str();
  auto *_dflt_ele_ref =
      dflt_ele_ref.has_value() ? dflt_ele_ref->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_dflt_ele_start = dflt_ele_start.has_value() ? dflt_ele_start->get().get_fortran_ptr()
                                                     : nullptr; // input, optional
  auto *_dflt_ele =
      dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_component = dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int dflt_eval_point_lvalue;
  auto *_dflt_eval_point{&dflt_eval_point_lvalue};
  if (dflt_eval_point.has_value()) {
    dflt_eval_point_lvalue = dflt_eval_point.value();
  } else {
    _dflt_eval_point = nullptr;
  }
  double dflt_s_offset_lvalue;
  auto *_dflt_s_offset{&dflt_s_offset_lvalue};
  if (dflt_s_offset.has_value()) {
    dflt_s_offset_lvalue = dflt_s_offset.value();
  } else {
    _dflt_s_offset = nullptr;
  }
  fortran_tao_evaluate_lat_or_beam_data(
      /* bool& */ _err,
      /* const char* */ _data_name,
      /* void* */ values.get_fortran_ptr(),
      /* bool& */ print_err,
      /* const char* */ _default_source,
      /* void* */ _dflt_ele_ref,
      /* void* */ _dflt_ele_start,
      /* void* */ _dflt_ele,
      /* const char* */ _dflt_component,
      /* int* */ _dflt_uni,
      /* int* */ _dflt_eval_point,
      /* double* */ _dflt_s_offset
  );
  return _err;
}
bool Tao::tao_evaluate_stack_old(
    TaoEvalNodeStructArray1D stack,
    int n_size_in,
    bool use_good_user,
    RealAlloc1D &value,
    bool print_err,
    std::string expression,
    std::optional<TaoExpressionInfoStructAlloc1D> info_in
) {
  // stack: TaoEvalNodeStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _stack_desc;
  _stack_desc.rank = 1;
  _stack_desc.data_ptr = stack.data();
  _stack_desc.dims[0] = stack.size();
  _stack_desc.strides[0] = 1;
  // intent=inout allocatable general array
  bool _err_flag{};
  auto _expression = expression.c_str();
  // intent=inout allocatable type array
  auto *_info_in = info_in.has_value() ? info_in->get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_evaluate_stack_old(
      /* Bmad::array_descriptor_t& */ _stack_desc,
      /* int& */ n_size_in,
      /* bool& */ use_good_user,
      /* void* */ value.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool& */ print_err,
      /* const char* */ _expression,
      /* void* */ _info_in
  );
  return _err_flag;
}
bool Tao::tao_evaluate_tree(
    TaoEvalNodeStruct &tao_tree,
    int n_size,
    bool use_good_user,
    RealAlloc1D &value,
    bool print_err,
    std::string expression,
    std::optional<TaoExpressionInfoStructAlloc1D> info_in
) {
  // intent=inout allocatable general array
  bool _err_flag{};
  auto _expression = expression.c_str();
  // intent=inout allocatable type array
  auto *_info_in = info_in.has_value() ? info_in->get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_evaluate_tree(
      /* void* */ tao_tree.get_fortran_ptr(),
      /* int& */ n_size,
      /* bool& */ use_good_user,
      /* void* */ value.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* bool& */ print_err,
      /* const char* */ _expression,
      /* void* */ _info_in
  );
  return _err_flag;
}
void Tao::tao_evaluate_tune(std::string q_str, double q0, bool delta_input, double q_val) {
  auto _q_str = q_str.c_str();
  fortran_tao_evaluate_tune(
      /* const char* */ _q_str,
      /* double& */ q0,
      /* bool& */ delta_input,
      /* double& */ q_val
  );
}
std::string
Tao::tao_expression_hash_substitute(std::string expression_in, optional_ref<EleStruct> eval_ele) {
  auto _expression_in = expression_in.c_str();
  char _expression_out[4096];
  auto *_eval_ele =
      eval_ele.has_value() ? eval_ele->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_expression_hash_substitute(
      /* const char* */ _expression_in,
      /* const char* */ _expression_out,
      /* void* */ _eval_ele
  );
  return _expression_out;
}
std::string Tao::tao_expression_tree_to_string(
    TaoEvalNodeStruct &tree,
    std::optional<bool> include_root,
    std::optional<int> n_node,
    optional_ref<TaoEvalNodeStruct> parent
) {
  bool include_root_lvalue;
  auto *_include_root{&include_root_lvalue};
  if (include_root.has_value()) {
    include_root_lvalue = include_root.value();
  } else {
    _include_root = nullptr;
  }
  int n_node_lvalue;
  auto *_n_node{&n_node_lvalue};
  if (n_node.has_value()) {
    n_node_lvalue = n_node.value();
  } else {
    _n_node = nullptr;
  }
  auto *_parent = parent.has_value() ? parent->get().get_fortran_ptr() : nullptr; // input, optional
  char _str_out[4096];
  fortran_tao_expression_tree_to_string(
      /* void* */ tree.get_fortran_ptr(),
      /* bool* */ _include_root,
      /* int* */ _n_node,
      /* void* */ _parent,
      /* const char* */ _str_out
  );
  return _str_out;
}
Tao::TaoFindPlotRegion
Tao::tao_find_plot_region(std::string where, std::optional<bool> print_flag) {
  bool _err{};
  auto _where = where.c_str();
  void *_region;
  bool print_flag_lvalue;
  auto *_print_flag{&print_flag_lvalue};
  if (print_flag.has_value()) {
    print_flag_lvalue = print_flag.value();
  } else {
    _print_flag = nullptr;
  }
  fortran_tao_find_plot_region(
      /* bool& */ _err,
      /* const char* */ _where,
      /* void* */ &_region,
      /* bool* */ _print_flag
  );
  return TaoFindPlotRegion{
      _err,
      std::move((_region ? std::make_optional<TaoPlotRegionStruct>(_region) : std::nullopt))
  };
}
void Tao::tao_fixer(std::string switch_, std::string word1, std::string word2) {
  auto _switch_ = switch_.c_str();
  auto _word1 = word1.c_str();
  auto _word2 = word2.c_str();
  fortran_tao_fixer(/* const char* */ _switch_, /* const char* */ _word1, /* const char* */ _word2);
}
Tao::TaoFloorToScreen
Tao::tao_floor_to_screen(TaoGraphStruct &graph, FixedArray1D<Real, 3> r_floor) {
  // r_floor: inout NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _r_floor_desc;
  _r_floor_desc.rank = 1;
  _r_floor_desc.data_ptr = r_floor.data();
  _r_floor_desc.dims[0] = r_floor.size();
  double _x_screen{};
  double _y_screen{};
  fortran_tao_floor_to_screen(
      /* void* */ graph.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _r_floor_desc,
      /* double& */ _x_screen,
      /* double& */ _y_screen
  );
  return TaoFloorToScreen{_x_screen, _y_screen};
}
FloorPositionStruct
Tao::tao_floor_to_screen_coords(TaoGraphStruct &graph, FloorPositionStruct &floor) {
  FloorPositionStruct _screen;
  fortran_tao_floor_to_screen_coords(
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ floor.get_fortran_ptr(),
      /* void* */ _screen.get_fortran_ptr()
  );
  return std::move(_screen);
}
bool Tao::tao_geodesic_lm_optimizer() {
  bool _abort{};
  fortran_tao_geodesic_lm_optimizer(/* bool& */ _abort);
  return _abort;
}
void Tao::tao_get_data(
    optional_ref<RealAlloc1D> data_value,
    optional_ref<RealAlloc1D> data_weight,
    optional_ref<RealAlloc1D> data_meas_value,
    optional_ref<IntAlloc1D> data_ix_dModel
) {
  // intent=inout allocatable general array
  auto *_data_value =
      data_value.has_value() ? data_value->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_data_weight =
      data_weight.has_value() ? data_weight->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_data_meas_value = data_meas_value.has_value() ? data_meas_value->get().get_fortran_ptr()
                                                       : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_data_ix_dModel = data_ix_dModel.has_value() ? data_ix_dModel->get().get_fortran_ptr()
                                                     : nullptr; // input, optional
  fortran_tao_get_data(
      /* void* */ _data_value,
      /* void* */ _data_weight,
      /* void* */ _data_meas_value,
      /* void* */ _data_ix_dModel
  );
}
Tao::TaoGetOptVars Tao::tao_get_opt_vars(
    optional_ref<RealAlloc1D> var_value,
    optional_ref<RealAlloc1D> var_step,
    optional_ref<RealAlloc1D> var_delta,
    optional_ref<RealAlloc1D> var_weight,
    optional_ref<IntAlloc1D> var_ix
) {
  // intent=inout allocatable general array
  auto *_var_value =
      var_value.has_value() ? var_value->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_var_step =
      var_step.has_value() ? var_step->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_var_delta =
      var_delta.has_value() ? var_delta->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_var_weight =
      var_weight.has_value() ? var_weight->get().get_fortran_ptr() : nullptr; // input, optional
  // intent=inout allocatable general array
  auto *_var_ix = var_ix.has_value() ? var_ix->get().get_fortran_ptr() : nullptr; // input, optional
  bool _ignore_if_weight_is_zero{};
  bool _ignore_if_not_limited{};
  fortran_tao_get_opt_vars(
      /* void* */ _var_value,
      /* void* */ _var_step,
      /* void* */ _var_delta,
      /* void* */ _var_weight,
      /* void* */ _var_ix,
      /* bool& */ _ignore_if_weight_is_zero,
      /* bool& */ _ignore_if_not_limited
  );
  return TaoGetOptVars{_ignore_if_weight_is_zero, _ignore_if_not_limited};
}
std::string Tao::tao_get_user_input(
    std::optional<std::string> prompt_str,
    std::optional<bool> wait_flag,
    std::optional<std::string> cmd_in
) {
  char _cmd_out[4096];
  const char *_prompt_str = prompt_str.has_value() ? prompt_str->c_str() : nullptr;
  bool wait_flag_lvalue;
  auto *_wait_flag{&wait_flag_lvalue};
  if (wait_flag.has_value()) {
    wait_flag_lvalue = wait_flag.value();
  } else {
    _wait_flag = nullptr;
  }
  const char *_cmd_in = cmd_in.has_value() ? cmd_in->c_str() : nullptr;
  fortran_tao_get_user_input(
      /* const char* */ _cmd_out,
      /* const char* */ _prompt_str,
      /* bool* */ _wait_flag,
      /* const char* */ _cmd_in
  );
  return _cmd_out;
}
void Tao::tao_graph_controller_setup(TaoGraphStruct &graph) {
  fortran_tao_graph_controller_setup(/* void* */ graph.get_fortran_ptr());
}
void Tao::tao_graph_data_setup(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_graph_data_setup(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_graph_data_slice_setup(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_graph_data_slice_setup(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_graph_dynamic_aperture_setup(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_graph_dynamic_aperture_setup(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_graph_histogram_setup(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_graph_histogram_setup(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
std::string Tao::tao_graph_name(TaoGraphStruct &graph, std::optional<bool> use_region) {
  bool use_region_lvalue;
  auto *_use_region{&use_region_lvalue};
  if (use_region.has_value()) {
    use_region_lvalue = use_region.value();
  } else {
    _use_region = nullptr;
  }
  char _graph_name[4096];
  fortran_tao_graph_name(
      /* void* */ graph.get_fortran_ptr(),
      /* bool* */ _use_region,
      /* const char* */ _graph_name
  );
  return _graph_name;
}
void Tao::tao_graph_phase_space_setup(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_graph_phase_space_setup(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
Tao::TaoGraphSMinMaxCalc
Tao::tao_graph_s_min_max_calc(TaoGraphStruct &graph, BranchStruct &branch) {
  double _s_min{};
  double _s_max{};
  fortran_tao_graph_s_min_max_calc(
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* double& */ _s_min,
      /* double& */ _s_max
  );
  return TaoGraphSMinMaxCalc{_s_min, _s_max};
}
void Tao::tao_graph_setup(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_graph_setup(/* void* */ plot.get_fortran_ptr(), /* void* */ graph.get_fortran_ptr());
}
bool Tao::tao_init() {
  bool _err_flag{};
  fortran_tao_init(/* bool& */ _err_flag);
  return _err_flag;
}
void Tao::tao_init_beam_in_universe(
    TaoUniverseStruct &u,
    BeamInitStruct &beam_init,
    std::string track_start,
    std::string track_end,
    double comb_ds_save
) {
  auto _track_start = track_start.c_str();
  auto _track_end = track_end.c_str();
  fortran_tao_init_beam_in_universe(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* const char* */ _track_start,
      /* const char* */ _track_end,
      /* double& */ comb_ds_save
  );
}
void Tao::tao_init_beams(std::string init_file) {
  auto _init_file = init_file.c_str();
  fortran_tao_init_beams(/* const char* */ _init_file);
}
void Tao::tao_init_data(std::string data_file) {
  auto _data_file = data_file.c_str();
  fortran_tao_init_data(/* const char* */ _data_file);
}
void Tao::tao_init_data_end_stuff() { fortran_tao_init_data_end_stuff(); }
void Tao::tao_init_data_in_universe(
    TaoUniverseStruct &u,
    int n_d2_add,
    std::optional<bool> keep_existing_data
) {
  bool keep_existing_data_lvalue;
  auto *_keep_existing_data{&keep_existing_data_lvalue};
  if (keep_existing_data.has_value()) {
    keep_existing_data_lvalue = keep_existing_data.value();
  } else {
    _keep_existing_data = nullptr;
  }
  fortran_tao_init_data_in_universe(
      /* void* */ u.get_fortran_ptr(),
      /* int& */ n_d2_add,
      /* bool* */ _keep_existing_data
  );
}
void Tao::tao_init_dynamic_aperture(std::string init_file) {
  auto _init_file = init_file.c_str();
  fortran_tao_init_dynamic_aperture(/* const char* */ _init_file);
}
bool Tao::tao_init_find_elements(
    TaoUniverseStruct &u,
    std::string search_string,
    ElePointerStructAlloc1D eles,
    std::optional<std::string> attribute
) {
  auto _search_string = search_string.c_str();
  // intent=inout allocatable type array
  const char *_attribute = attribute.has_value() ? attribute->c_str() : nullptr;
  bool _found_one{};
  fortran_tao_init_find_elements(
      /* void* */ u.get_fortran_ptr(),
      /* const char* */ _search_string,
      /* void* */ eles.get_fortran_ptr(),
      /* const char* */ _attribute,
      /* bool& */ _found_one
  );
  return _found_one;
}
void Tao::tao_init_global(std::string init_file) {
  auto _init_file = init_file.c_str();
  fortran_tao_init_global(/* const char* */ _init_file);
}
void Tao::tao_init_lattice(std::string lat_file, bool err_flag) {
  auto _lat_file = lat_file.c_str();
  fortran_tao_init_lattice(/* const char* */ _lat_file, /* bool& */ err_flag);
}
void Tao::tao_init_plotting(std::string plot_file) {
  auto _plot_file = plot_file.c_str();
  fortran_tao_init_plotting(/* const char* */ _plot_file);
}
void Tao::tao_init_variables(std::string var_file) {
  auto _var_file = var_file.c_str();
  fortran_tao_init_variables(/* const char* */ _var_file);
}
Tao::TaoInjectBeam
Tao::tao_inject_beam(TaoUniverseStruct &u, TaoLatticeStruct &model, int ix_branch) {
  BeamStruct _beam;
  bool _init_ok{};
  fortran_tao_inject_beam(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ model.get_fortran_ptr(),
      /* int& */ ix_branch,
      /* void* */ _beam.get_fortran_ptr(),
      /* bool& */ _init_ok
  );
  return TaoInjectBeam{std::move(_beam), _init_ok};
}
void Tao::tao_inject_particle(TaoUniverseStruct &u, TaoLatticeStruct &model, int ix_branch) {
  fortran_tao_inject_particle(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ model.get_fortran_ptr(),
      /* int& */ ix_branch
  );
}
Tao::TaoIsValidName Tao::tao_is_valid_name(std::string name) {
  auto _name = name.c_str();
  char _why_invalid[4096];
  bool _is_valid{};
  fortran_tao_is_valid_name(
      /* const char* */ _name,
      /* const char* */ _why_invalid,
      /* bool& */ _is_valid
  );
  return TaoIsValidName{_why_invalid, _is_valid};
}
void Tao::tao_json_cmd(std::string input_str) {
  auto _input_str = input_str.c_str();
  fortran_tao_json_cmd(/* const char* */ _input_str);
}
void Tao::tao_key_info_to_str(
    int ix_key,
    int ix_min_key,
    int ix_max_key,
    std::string key_str,
    std::string header_str
) {
  auto _key_str = key_str.c_str();
  auto _header_str = header_str.c_str();
  fortran_tao_key_info_to_str(
      /* int& */ ix_key,
      /* int& */ ix_min_key,
      /* int& */ ix_max_key,
      /* const char* */ _key_str,
      /* const char* */ _header_str
  );
}
bool Tao::tao_lat_bookkeeper(TaoUniverseStruct &u, TaoLatticeStruct &tao_lat) {
  bool _err_flag{};
  fortran_tao_lat_bookkeeper(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* bool& */ _err_flag
  );
  return _err_flag;
}
double Tao::tao_lat_emit_calc(int plane, int emit_type, EleStruct &ele, NormalModesStruct &modes) {
  double _emit{};
  fortran_tao_lat_emit_calc(
      /* int& */ plane,
      /* int& */ emit_type,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ modes.get_fortran_ptr(),
      /* double& */ _emit
  );
  return _emit;
}
void Tao::tao_lat_sigma_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_lat_sigma
) {
  auto _data_type = data_type.c_str();
  auto _data_source = data_source.c_str();
  fortran_tao_lat_sigma_calc_needed(
      /* const char* */ _data_type,
      /* const char* */ _data_source,
      /* bool& */ do_lat_sigma
  );
}
bool Tao::tao_lat_sigma_track(
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    std::optional<bool> print_err,
    std::optional<bool> force_calc
) {
  bool _calc_ok{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool force_calc_lvalue;
  auto *_force_calc{&force_calc_lvalue};
  if (force_calc.has_value()) {
    force_calc_lvalue = force_calc.value();
  } else {
    _force_calc = nullptr;
  }
  fortran_tao_lat_sigma_track(
      /* void* */ tao_lat.get_fortran_ptr(),
      /* bool& */ _calc_ok,
      /* int& */ ix_branch,
      /* bool* */ _print_err,
      /* bool* */ _force_calc
  );
  return _calc_ok;
}
void Tao::tao_lattice_branches_equal_tao_lattice_branches(
    TaoLatticeBranchStructArray1D tlb1,
    TaoLatticeBranchStructArray1D tlb2
) {
  // tlb1: TaoLatticeBranchStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _tlb1_desc;
  _tlb1_desc.rank = 1;
  _tlb1_desc.data_ptr = tlb1.data();
  _tlb1_desc.dims[0] = tlb1.size();
  _tlb1_desc.strides[0] = 1;
  // tlb2: TaoLatticeBranchStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _tlb2_desc;
  _tlb2_desc.rank = 1;
  _tlb2_desc.data_ptr = tlb2.data();
  _tlb2_desc.dims[0] = tlb2.size();
  _tlb2_desc.strides[0] = 1;
  fortran_tao_lattice_branches_equal_tao_lattice_branches(
      /* Bmad::array_descriptor_t& */ _tlb1_desc,
      /* Bmad::array_descriptor_t& */ _tlb2_desc
  );
}
Tao::TaoLatticeCalc Tao::tao_lattice_calc() {
  bool _calc_ok{};
  bool _print_err{};
  fortran_tao_lattice_calc(/* bool& */ _calc_ok, /* bool& */ _print_err);
  return TaoLatticeCalc{_calc_ok, _print_err};
}
void Tao::tao_lattice_equal_tao_lattice(TaoLatticeStruct &lat1, TaoLatticeStruct &lat2) {
  fortran_tao_lattice_equal_tao_lattice(
      /* void* */ lat1.get_fortran_ptr(),
      /* void* */ lat2.get_fortran_ptr()
  );
}
bool Tao::tao_limit_calc() {
  bool _limited{};
  fortran_tao_limit_calc(/* bool& */ _limited);
  return _limited;
}
bool Tao::tao_lm_optimizer() {
  bool _abort{};
  fortran_tao_lm_optimizer(/* bool& */ _abort);
  return _abort;
}
bool Tao::tao_lmdif_optimizer() {
  bool _abort{};
  fortran_tao_lmdif_optimizer(/* bool& */ _abort);
  return _abort;
}
void Tao::tao_load_this_datum(
    FArray1D<Real> &vec,
    EleStruct &ele_ref,
    EleStruct &ele_start,
    EleStruct &ele,
    double datum_value,
    bool valid_value,
    TaoDataStruct &datum,
    BranchStruct &branch,
    std::optional<std::string> why_invalid,
    optional_ref<BoolAlloc1D> good
) {
  // vec: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  auto _ele_ref = &ele_ref; // input, required, pointer
  auto _ele_start = &ele_start; // input, required, pointer
  auto _ele = &ele; // input, required, pointer
  const char *_why_invalid = why_invalid.has_value() ? why_invalid->c_str() : nullptr;
  // intent=inout allocatable general array
  auto *_good = good.has_value() ? good->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_load_this_datum(
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* void* */ &ele_ref,
      /* void* */ &ele_start,
      /* void* */ &ele,
      /* double& */ datum_value,
      /* bool& */ valid_value,
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* const char* */ _why_invalid,
      /* void* */ _good
  );
}
bool Tao::tao_locate_all_elements(
    std::string ele_list,
    ElePointerStructAlloc1D eles,
    std::optional<bool> ignore_blank
) {
  auto _ele_list = ele_list.c_str();
  // intent=inout allocatable type array
  bool _err{};
  bool ignore_blank_lvalue;
  auto *_ignore_blank{&ignore_blank_lvalue};
  if (ignore_blank.has_value()) {
    ignore_blank_lvalue = ignore_blank.value();
  } else {
    _ignore_blank = nullptr;
  }
  fortran_tao_locate_all_elements(
      /* const char* */ _ele_list,
      /* void* */ eles.get_fortran_ptr(),
      /* bool& */ _err,
      /* bool* */ _ignore_blank
  );
  return _err;
}
bool Tao::tao_locate_elements(
    std::string ele_list,
    int ix_universe,
    ElePointerStructAlloc1D eles,
    std::optional<int> lat_type,
    std::optional<bool> ignore_blank,
    std::optional<int> err_stat_level,
    std::optional<bool> above_ubound_is_err,
    std::optional<int> ix_branch,
    std::optional<bool> multiple_eles_is_err
) {
  auto _ele_list = ele_list.c_str();
  // intent=inout allocatable type array
  bool _err{};
  int lat_type_lvalue;
  auto *_lat_type{&lat_type_lvalue};
  if (lat_type.has_value()) {
    lat_type_lvalue = lat_type.value();
  } else {
    _lat_type = nullptr;
  }
  bool ignore_blank_lvalue;
  auto *_ignore_blank{&ignore_blank_lvalue};
  if (ignore_blank.has_value()) {
    ignore_blank_lvalue = ignore_blank.value();
  } else {
    _ignore_blank = nullptr;
  }
  int err_stat_level_lvalue;
  auto *_err_stat_level{&err_stat_level_lvalue};
  if (err_stat_level.has_value()) {
    err_stat_level_lvalue = err_stat_level.value();
  } else {
    _err_stat_level = nullptr;
  }
  bool above_ubound_is_err_lvalue;
  auto *_above_ubound_is_err{&above_ubound_is_err_lvalue};
  if (above_ubound_is_err.has_value()) {
    above_ubound_is_err_lvalue = above_ubound_is_err.value();
  } else {
    _above_ubound_is_err = nullptr;
  }
  int ix_branch_lvalue;
  auto *_ix_branch{&ix_branch_lvalue};
  if (ix_branch.has_value()) {
    ix_branch_lvalue = ix_branch.value();
  } else {
    _ix_branch = nullptr;
  }
  bool multiple_eles_is_err_lvalue;
  auto *_multiple_eles_is_err{&multiple_eles_is_err_lvalue};
  if (multiple_eles_is_err.has_value()) {
    multiple_eles_is_err_lvalue = multiple_eles_is_err.value();
  } else {
    _multiple_eles_is_err = nullptr;
  }
  fortran_tao_locate_elements(
      /* const char* */ _ele_list,
      /* int& */ ix_universe,
      /* void* */ eles.get_fortran_ptr(),
      /* bool& */ _err,
      /* int* */ _lat_type,
      /* bool* */ _ignore_blank,
      /* int* */ _err_stat_level,
      /* bool* */ _above_ubound_is_err,
      /* int* */ _ix_branch,
      /* bool* */ _multiple_eles_is_err
  );
  return _err;
}
void Tao::tao_mark_lattice_ele(LatStruct &lat) {
  fortran_tao_mark_lattice_ele(/* void* */ lat.get_fortran_ptr());
}
Tao::TaoMerit Tao::tao_merit() {
  bool _calc_ok{};
  double _this_merit{};
  fortran_tao_merit(/* bool& */ _calc_ok, /* double& */ _this_merit);
  return TaoMerit{_calc_ok, _this_merit};
}
std::string Tao::tao_next_word(std::string line) {
  auto _line = line.c_str();
  char _word[4096];
  fortran_tao_next_word(/* const char* */ _line, /* const char* */ _word);
  return _word;
}
void Tao::tao_one_turn_map_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_one_turn_map
) {
  auto _data_type = data_type.c_str();
  auto _data_source = data_source.c_str();
  fortran_tao_one_turn_map_calc_needed(
      /* const char* */ _data_type,
      /* const char* */ _data_source,
      /* bool& */ do_one_turn_map
  );
}
int Tao::tao_open_file(
    std::string file,
    std::string file_name,
    int error_severity,
    std::optional<bool> binary
) {
  auto _file = file.c_str();
  int _iunit{};
  auto _file_name = file_name.c_str();
  bool binary_lvalue;
  auto *_binary{&binary_lvalue};
  if (binary.has_value()) {
    binary_lvalue = binary.value();
  } else {
    _binary = nullptr;
  }
  fortran_tao_open_file(
      /* const char* */ _file,
      /* int& */ _iunit,
      /* const char* */ _file_name,
      /* int& */ error_severity,
      /* bool* */ _binary
  );
  return _iunit;
}
Tao::TaoOpenScratchFile Tao::tao_open_scratch_file() {
  bool _err{};
  int _iu{};
  fortran_tao_open_scratch_file(/* bool& */ _err, /* int& */ _iu);
  return TaoOpenScratchFile{_err, _iu};
}
std::string Tao::tao_optimization_status(TaoDataStruct &datum) {
  char _why_str[4096];
  fortran_tao_optimization_status(/* void* */ datum.get_fortran_ptr(), /* const char* */ _why_str);
  return _why_str;
}
void Tao::tao_orbit_beta_wave_anal(TaoPlotStruct &plot) {
  fortran_tao_orbit_beta_wave_anal(/* void* */ plot.get_fortran_ptr());
}
TaoBuildingWallPointStruct Tao::tao_oreint_building_wall_pt(TaoBuildingWallPointStruct &pt_in) {
  TaoBuildingWallPointStruct _pt_out;
  fortran_tao_oreint_building_wall_pt(
      /* void* */ pt_in.get_fortran_ptr(),
      /* void* */ _pt_out.get_fortran_ptr()
  );
  return std::move(_pt_out);
}
Tao::TaoParamValueAtS Tao::tao_param_value_at_s(
    std::string dat_name,
    EleStruct &ele_to_s,
    EleStruct &ele_here,
    CoordStruct &orbit
) {
  auto _dat_name = dat_name.c_str();
  bool _err_flag{};
  char _why_invalid[4096];
  bool _print_err{};
  bool _bad_datum{};
  double _value{};
  fortran_tao_param_value_at_s(
      /* const char* */ _dat_name,
      /* void* */ ele_to_s.get_fortran_ptr(),
      /* void* */ ele_here.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* bool& */ _err_flag,
      /* const char* */ _why_invalid,
      /* bool& */ _print_err,
      /* bool& */ _bad_datum,
      /* double& */ _value
  );
  return TaoParamValueAtS{_err_flag, _why_invalid, _print_err, _bad_datum, _value};
}
void Tao::tao_param_value_routine(
    std::string str,
    bool use_good_user,
    std::string saved_prefix,
    TaoEvalNodeStruct &stack,
    bool err_flag,
    bool print_err,
    std::optional<std::string> dflt_component,
    std::optional<std::string> dflt_source,
    optional_ref<EleStruct> dflt_ele_ref,
    optional_ref<EleStruct> dflt_ele_start,
    optional_ref<EleStruct> dflt_ele,
    std::optional<std::string> dflt_dat_or_var_index,
    std::optional<int> dflt_uni,
    std::optional<int> dflt_eval_point,
    std::optional<double> dflt_s_offset,
    optional_ref<CoordStruct> dflt_orbit,
    optional_ref<TaoDataStruct> datum
) {
  auto _str = str.c_str();
  auto _saved_prefix = saved_prefix.c_str();
  const char *_dflt_component = dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  const char *_dflt_source = dflt_source.has_value() ? dflt_source->c_str() : nullptr;
  auto *_dflt_ele_ref =
      dflt_ele_ref.has_value() ? dflt_ele_ref->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_dflt_ele_start = dflt_ele_start.has_value() ? dflt_ele_start->get().get_fortran_ptr()
                                                     : nullptr; // input, optional
  auto *_dflt_ele =
      dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr() : nullptr; // input, optional
  const char *_dflt_dat_or_var_index =
      dflt_dat_or_var_index.has_value() ? dflt_dat_or_var_index->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int dflt_eval_point_lvalue;
  auto *_dflt_eval_point{&dflt_eval_point_lvalue};
  if (dflt_eval_point.has_value()) {
    dflt_eval_point_lvalue = dflt_eval_point.value();
  } else {
    _dflt_eval_point = nullptr;
  }
  double dflt_s_offset_lvalue;
  auto *_dflt_s_offset{&dflt_s_offset_lvalue};
  if (dflt_s_offset.has_value()) {
    dflt_s_offset_lvalue = dflt_s_offset.value();
  } else {
    _dflt_s_offset = nullptr;
  }
  auto *_dflt_orbit =
      dflt_orbit.has_value() ? dflt_orbit->get().get_fortran_ptr() : nullptr; // input, optional
  auto *_datum = datum.has_value() ? datum->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_tao_param_value_routine(
      /* const char* */ _str,
      /* bool& */ use_good_user,
      /* const char* */ _saved_prefix,
      /* void* */ stack.get_fortran_ptr(),
      /* bool& */ err_flag,
      /* bool& */ print_err,
      /* const char* */ _dflt_component,
      /* const char* */ _dflt_source,
      /* void* */ _dflt_ele_ref,
      /* void* */ _dflt_ele_start,
      /* void* */ _dflt_ele,
      /* const char* */ _dflt_dat_or_var_index,
      /* int* */ _dflt_uni,
      /* int* */ _dflt_eval_point,
      /* double* */ _dflt_s_offset,
      /* void* */ _dflt_orbit,
      /* void* */ _datum
  );
}
bool Tao::tao_parse_command_args(std::optional<std::string> cmd_line) {
  bool _error{};
  const char *_cmd_line = cmd_line.has_value() ? cmd_line->c_str() : nullptr;
  fortran_tao_parse_command_args(/* bool& */ _error, /* const char* */ _cmd_line);
  return _error;
}
Tao::TaoParseElementParamStr Tao::tao_parse_element_param_str(std::string in_str) {
  bool _err{};
  auto _in_str = in_str.c_str();
  char _uni[4096];
  char _element[4096];
  char _parameter[4096];
  int _where{};
  char _component[4096];
  fortran_tao_parse_element_param_str(
      /* bool& */ _err,
      /* const char* */ _in_str,
      /* const char* */ _uni,
      /* const char* */ _element,
      /* const char* */ _parameter,
      /* int& */ _where,
      /* const char* */ _component
  );
  return TaoParseElementParamStr{_err, _uni, _element, _parameter, _where, _component};
}
bool Tao::tao_particle_data_value(
    std::string data_type,
    CoordStructArray1D p,
    RealAlloc1D &value,
    EleStruct &ele,
    int ix_bunch
) {
  auto _data_type = data_type.c_str();
  // p: CoordStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _p_desc;
  _p_desc.rank = 1;
  _p_desc.data_ptr = p.data();
  _p_desc.dims[0] = p.size();
  _p_desc.strides[0] = 1;
  // intent=inout allocatable general array
  bool _err{};
  fortran_tao_particle_data_value(
      /* const char* */ _data_type,
      /* Bmad::array_descriptor_t& */ _p_desc,
      /* void* */ value.get_fortran_ptr(),
      /* bool& */ _err,
      /* void* */ ele.get_fortran_ptr(),
      /* int& */ ix_bunch
  );
  return _err;
}
void Tao::tao_pause_cmd(double time) { fortran_tao_pause_cmd(/* double& */ time); }
int Tao::tao_phase_space_axis_index(std::string data_type, bool err) {
  auto _data_type = data_type.c_str();
  int _ix_axis{};
  fortran_tao_phase_space_axis_index(
      /* const char* */ _data_type,
      /* bool& */ err,
      /* int& */ _ix_axis
  );
  return _ix_axis;
}
void Tao::tao_phase_wave_anal(TaoPlotStruct &plot) {
  fortran_tao_phase_wave_anal(/* void* */ plot.get_fortran_ptr());
}
Tao::TaoPickUniverse Tao::tao_pick_universe(
    std::string name_in,
    BoolAlloc1D &picked,
    std::optional<int> dflt_uni,
    std::optional<bool> pure_uni
) {
  auto _name_in = name_in.c_str();
  char _name_out[4096];
  // intent=inout allocatable general array
  bool _err{};
  int _ix_uni{};
  bool _explicit_uni{};
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  bool pure_uni_lvalue;
  auto *_pure_uni{&pure_uni_lvalue};
  if (pure_uni.has_value()) {
    pure_uni_lvalue = pure_uni.value();
  } else {
    _pure_uni = nullptr;
  }
  fortran_tao_pick_universe(
      /* const char* */ _name_in,
      /* const char* */ _name_out,
      /* void* */ picked.get_fortran_ptr(),
      /* bool& */ _err,
      /* int& */ _ix_uni,
      /* bool& */ _explicit_uni,
      /* int* */ _dflt_uni,
      /* bool* */ _pure_uni
  );
  return TaoPickUniverse{_name_out, _err, _ix_uni, _explicit_uni};
}
void Tao::tao_pipe_cmd(std::string input_str) {
  auto _input_str = input_str.c_str();
  fortran_tao_pipe_cmd(/* const char* */ _input_str);
}
void Tao::tao_place_cmd(std::string where, std::string who, std::optional<bool> no_buffer) {
  auto _where = where.c_str();
  auto _who = who.c_str();
  bool no_buffer_lvalue;
  auto *_no_buffer{&no_buffer_lvalue};
  if (no_buffer.has_value()) {
    no_buffer_lvalue = no_buffer.value();
  } else {
    _no_buffer = nullptr;
  }
  fortran_tao_place_cmd(/* const char* */ _where, /* const char* */ _who, /* bool* */ _no_buffer);
}
void Tao::tao_plot_cmd(std::string where, std::string component) {
  auto _where = where.c_str();
  auto _component = component.c_str();
  fortran_tao_plot_cmd(/* const char* */ _where, /* const char* */ _component);
}
void Tao::tao_plot_data(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_plot_data(/* void* */ plot.get_fortran_ptr(), /* void* */ graph.get_fortran_ptr());
}
void Tao::tao_plot_histogram(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_plot_histogram(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_plot_key_table(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_plot_key_table(
      /* void* */ plot.get_fortran_ptr(),
      /* void* */ graph.get_fortran_ptr()
  );
}
void Tao::tao_plot_setup() { fortran_tao_plot_setup(); }
TaoPlotStruct Tao::tao_plot_struct_transfer(TaoPlotStruct &plot_in) {
  TaoPlotStruct _plot_out;
  fortran_tao_plot_struct_transfer(
      /* void* */ plot_in.get_fortran_ptr(),
      /* void* */ _plot_out.get_fortran_ptr()
  );
  return std::move(_plot_out);
}
void Tao::tao_plot_wave(TaoPlotStruct &plot, TaoGraphStruct &graph) {
  fortran_tao_plot_wave(/* void* */ plot.get_fortran_ptr(), /* void* */ graph.get_fortran_ptr());
}
std::optional<TaoEleShapeStruct> Tao::tao_pointer_to_building_wall_shape(std::string wall_name) {
  auto _wall_name = wall_name.c_str();
  void *_e_shape;
  fortran_tao_pointer_to_building_wall_shape(/* const char* */ _wall_name, /* void* */ &_e_shape);
  return std::move((_e_shape ? std::make_optional<TaoEleShapeStruct>(_e_shape) : std::nullopt));
}
void Tao::tao_pointer_to_datum(
    TaoD1DataStruct &d1,
    std::string ele_name,
    TaoDataStruct &datum_ptr
) {
  auto _ele_name = ele_name.c_str();
  auto _datum_ptr = &datum_ptr; // input, required, pointer
  fortran_tao_pointer_to_datum(
      /* void* */ d1.get_fortran_ptr(),
      /* const char* */ _ele_name,
      /* void* */ &datum_ptr
  );
}
Tao::TaoPointerToDatumEle Tao::tao_pointer_to_datum_ele(
    LatStruct &lat,
    std::string ele_name,
    int ix_ele,
    TaoDataStruct &datum,
    std::optional<bool> print_err
) {
  auto _ele_name = ele_name.c_str();
  bool _valid{};
  char _why_invalid[4096];
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  void *_ele;
  fortran_tao_pointer_to_datum_ele(
      /* void* */ lat.get_fortran_ptr(),
      /* const char* */ _ele_name,
      /* int& */ ix_ele,
      /* void* */ datum.get_fortran_ptr(),
      /* bool& */ _valid,
      /* const char* */ _why_invalid,
      /* bool* */ _print_err,
      /* void* */ &_ele
  );
  return TaoPointerToDatumEle{
      _valid,
      _why_invalid,
      std::move((_ele ? std::make_optional<EleStruct>(_ele) : std::nullopt))
  };
}
Tao::TaoPointerToEleShape Tao::tao_pointer_to_ele_shape(
    int ix_uni,
    EleStruct &ele,
    TaoEleShapeStructArray1D ele_shape,
    std::optional<int> ix_shape_min
) {
  // ele_shape: TaoEleShapeStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _ele_shape_desc;
  _ele_shape_desc.rank = 1;
  _ele_shape_desc.data_ptr = ele_shape.data();
  _ele_shape_desc.dims[0] = ele_shape.size();
  _ele_shape_desc.strides[0] = 1;
  char _dat_var_name[4096];
  double _dat_var_value{};
  int ix_shape_min_lvalue;
  auto *_ix_shape_min{&ix_shape_min_lvalue};
  if (ix_shape_min.has_value()) {
    ix_shape_min_lvalue = ix_shape_min.value();
  } else {
    _ix_shape_min = nullptr;
  }
  void *_e_shape;
  fortran_tao_pointer_to_ele_shape(
      /* int& */ ix_uni,
      /* void* */ ele.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _ele_shape_desc,
      /* const char* */ _dat_var_name,
      /* double& */ _dat_var_value,
      /* int* */ _ix_shape_min,
      /* void* */ &_e_shape
  );
  return TaoPointerToEleShape{
      _dat_var_name,
      _dat_var_value,
      std::move((_e_shape ? std::make_optional<TaoEleShapeStruct>(_e_shape) : std::nullopt))
  };
}
std::optional<TaoLatticeStruct>
Tao::tao_pointer_to_tao_lat(TaoUniverseStruct &u, std::optional<int> lat_type) {
  int lat_type_lvalue;
  auto *_lat_type{&lat_type_lvalue};
  if (lat_type.has_value()) {
    lat_type_lvalue = lat_type.value();
  } else {
    _lat_type = nullptr;
  }
  void *_tao_lat;
  fortran_tao_pointer_to_tao_lat(
      /* void* */ u.get_fortran_ptr(),
      /* int* */ _lat_type,
      /* void* */ &_tao_lat
  );
  return std::move((_tao_lat ? std::make_optional<TaoLatticeStruct>(_tao_lat) : std::nullopt));
}
std::optional<TaoUniverseStruct>
Tao::tao_pointer_to_universe(int ix_uni, std::optional<bool> neg2_to_default) {
  bool neg2_to_default_lvalue;
  auto *_neg2_to_default{&neg2_to_default_lvalue};
  if (neg2_to_default.has_value()) {
    neg2_to_default_lvalue = neg2_to_default.value();
  } else {
    _neg2_to_default = nullptr;
  }
  void *_u;
  fortran_tao_pointer_to_universe_int(
      /* int& */ ix_uni,
      /* bool* */ _neg2_to_default,
      /* void* */ &_u
  );
  return std::move((_u ? std::make_optional<TaoUniverseStruct>(_u) : std::nullopt));
}
std::optional<TaoUniverseStruct>
Tao::tao_pointer_to_universe(std::string string, std::optional<bool> neg2_to_default) {
  auto _string = string.c_str();
  bool neg2_to_default_lvalue;
  auto *_neg2_to_default{&neg2_to_default_lvalue};
  if (neg2_to_default.has_value()) {
    neg2_to_default_lvalue = neg2_to_default.value();
  } else {
    _neg2_to_default = nullptr;
  }
  void *_u;
  fortran_tao_pointer_to_universe_str(
      /* const char* */ _string,
      /* bool* */ _neg2_to_default,
      /* void* */ &_u
  );
  return std::move((_u ? std::make_optional<TaoUniverseStruct>(_u) : std::nullopt));
}
Tao::TaoPointerToUniverses Tao::tao_pointer_to_universes(
    std::string name_in,
    TaoUniversePointerStructAlloc1D unis,
    std::optional<int> dflt_uni
) {
  auto _name_in = name_in.c_str();
  // intent=inout allocatable type array
  bool _err{};
  char _name_out[4096];
  bool _explicit_uni{};
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  fortran_tao_pointer_to_universes(
      /* const char* */ _name_in,
      /* void* */ unis.get_fortran_ptr(),
      /* bool& */ _err,
      /* const char* */ _name_out,
      /* bool& */ _explicit_uni,
      /* int* */ _dflt_uni
  );
  return TaoPointerToUniverses{_err, _name_out, _explicit_uni};
}
bool Tao::tao_pointer_to_var_in_lattice(TaoVarStruct &var, int ix_uni, EleStruct &ele) {
  bool _err{};
  fortran_tao_pointer_to_var_in_lattice(
      /* void* */ var.get_fortran_ptr(),
      /* int& */ ix_uni,
      /* void* */ ele.get_fortran_ptr(),
      /* bool& */ _err
  );
  return _err;
}
bool Tao::tao_pointer_to_var_in_lattice2(TaoVarStruct &var, int ix_uni) {
  bool _err{};
  fortran_tao_pointer_to_var_in_lattice2(
      /* void* */ var.get_fortran_ptr(),
      /* int& */ ix_uni,
      /* bool& */ _err
  );
  return _err;
}
void Tao::tao_print_command_line_info() { fortran_tao_print_command_line_info(); }
void Tao::tao_ptc_normal_form(
    bool do_calc,
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    std::optional<int> rf_on
) {
  int rf_on_lvalue;
  auto *_rf_on{&rf_on_lvalue};
  if (rf_on.has_value()) {
    rf_on_lvalue = rf_on.value();
  } else {
    _rf_on = nullptr;
  }
  fortran_tao_ptc_normal_form(
      /* bool& */ do_calc,
      /* void* */ tao_lat.get_fortran_ptr(),
      /* int& */ ix_branch,
      /* int* */ _rf_on
  );
}
void Tao::tao_python_cmd(std::string input_str) {
  auto _input_str = input_str.c_str();
  fortran_tao_python_cmd(/* const char* */ _input_str);
}
void Tao::tao_quiet_set(std::string set) {
  auto _set = set.c_str();
  fortran_tao_quiet_set(/* const char* */ _set);
}
void Tao::tao_rad_int_calc_needed(std::string data_type, std::string data_source, bool do_rad_int) {
  auto _data_type = data_type.c_str();
  auto _data_source = data_source.c_str();
  fortran_tao_rad_int_calc_needed(
      /* const char* */ _data_type,
      /* const char* */ _data_source,
      /* bool& */ do_rad_int
  );
}
void Tao::tao_re_allocate_expression_info(
    TaoExpressionInfoStructAlloc1D info,
    int n,
    std::optional<bool> exact
) {
  // intent=inout allocatable type array
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_tao_re_allocate_expression_info(
      /* void* */ info.get_fortran_ptr(),
      /* int& */ n,
      /* bool* */ _exact
  );
}
void Tao::tao_re_associate_node_array(TaoEvalNodeStruct &tree, int n, std::optional<bool> exact) {
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_tao_re_associate_node_array(
      /* void* */ tree.get_fortran_ptr(),
      /* int& */ n,
      /* bool* */ _exact
  );
}
void Tao::tao_re_execute(std::string string, bool err) {
  auto _string = string.c_str();
  fortran_tao_re_execute(/* const char* */ _string, /* bool& */ err);
}
void Tao::tao_read_cmd(std::string which, std::string unis, std::string file, bool silent) {
  auto _which = which.c_str();
  auto _unis = unis.c_str();
  auto _file = file.c_str();
  fortran_tao_read_cmd(
      /* const char* */ _which,
      /* const char* */ _unis,
      /* const char* */ _file,
      /* bool& */ silent
  );
}
int Tao::tao_read_phase_space_index(std::string name, int ixc, std::optional<bool> print_err) {
  auto _name = name.c_str();
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _ix_ps{};
  fortran_tao_read_phase_space_index(
      /* const char* */ _name,
      /* int& */ ixc,
      /* bool* */ _print_err,
      /* int& */ _ix_ps
  );
  return _ix_ps;
}
void Tao::tao_regression_test() { fortran_tao_regression_test(); }
void Tao::tao_remove_blank_characters(std::string str) {
  auto _str = str.c_str();
  fortran_tao_remove_blank_characters(/* const char* */ _str);
}
bool Tao::tao_run_cmd(std::string which) {
  auto _which = which.c_str();
  bool _abort{};
  fortran_tao_run_cmd(/* const char* */ _which, /* bool& */ _abort);
  return _abort;
}
void Tao::tao_scale_cmd(
    std::string where,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis,
    std::optional<bool> include_wall,
    std::optional<std::string> gang,
    std::optional<bool> exact,
    std::optional<bool> turn_autoscale_off
) {
  auto _where = where.c_str();
  const char *_axis = axis.has_value() ? axis->c_str() : nullptr;
  bool include_wall_lvalue;
  auto *_include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  const char *_gang = gang.has_value() ? gang->c_str() : nullptr;
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  bool turn_autoscale_off_lvalue;
  auto *_turn_autoscale_off{&turn_autoscale_off_lvalue};
  if (turn_autoscale_off.has_value()) {
    turn_autoscale_off_lvalue = turn_autoscale_off.value();
  } else {
    _turn_autoscale_off = nullptr;
  }
  fortran_tao_scale_cmd(
      /* const char* */ _where,
      /* double& */ y_min_in,
      /* double& */ y_max_in,
      /* const char* */ _axis,
      /* bool* */ _include_wall,
      /* const char* */ _gang,
      /* bool* */ _exact,
      /* bool* */ _turn_autoscale_off
  );
}
Tao::TaoScaleGraph Tao::tao_scale_graph(
    TaoGraphStruct &graph,
    double y_min,
    double y_max,
    std::optional<std::string> axis,
    std::optional<bool> include_wall
) {
  const char *_axis = axis.has_value() ? axis->c_str() : nullptr;
  bool include_wall_lvalue;
  auto *_include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  // y_range: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _y_range_desc;
  _y_range_desc.rank = 1;
  FixedArray1D<Real, 2> _y_range;
  _y_range_desc.data_ptr = _y_range.data();
  _y_range_desc.dims[0] = _y_range.size();
  // y2_range: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _y2_range_desc;
  _y2_range_desc.rank = 1;
  FixedArray1D<Real, 2> _y2_range;
  _y2_range_desc.data_ptr = _y2_range.data();
  _y2_range_desc.dims[0] = _y2_range.size();
  fortran_tao_scale_graph(
      /* void* */ graph.get_fortran_ptr(),
      /* double& */ y_min,
      /* double& */ y_max,
      /* const char* */ _axis,
      /* bool* */ _include_wall,
      /* Bmad::array_descriptor_t& */ _y_range_desc,
      /* Bmad::array_descriptor_t& */ _y2_range_desc
  );
  return TaoScaleGraph{_y_range, _y2_range};
}
void Tao::tao_scale_ping_data(TaoUniverseStruct &u) {
  fortran_tao_scale_ping_data(/* void* */ u.get_fortran_ptr());
}
void Tao::tao_scale_plot(
    TaoPlotStruct &plot,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis,
    std::optional<bool> include_wall,
    std::optional<std::string> gang,
    std::optional<bool> skip_lat_layout
) {
  const char *_axis = axis.has_value() ? axis->c_str() : nullptr;
  bool include_wall_lvalue;
  auto *_include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  const char *_gang = gang.has_value() ? gang->c_str() : nullptr;
  bool skip_lat_layout_lvalue;
  auto *_skip_lat_layout{&skip_lat_layout_lvalue};
  if (skip_lat_layout.has_value()) {
    skip_lat_layout_lvalue = skip_lat_layout.value();
  } else {
    _skip_lat_layout = nullptr;
  }
  fortran_tao_scale_plot(
      /* void* */ plot.get_fortran_ptr(),
      /* double& */ y_min_in,
      /* double& */ y_max_in,
      /* const char* */ _axis,
      /* bool* */ _include_wall,
      /* const char* */ _gang,
      /* bool* */ _skip_lat_layout
  );
}
void Tao::tao_scratch_values_calc(
    EleStruct &ele_ref,
    EleStruct &ele_start,
    EleStruct &ele,
    TaoDataStruct &datum,
    BranchStruct &branch,
    CoordStructArray1D orbit
) {
  auto _ele_ref = &ele_ref; // input, required, pointer
  auto _ele_start = &ele_start; // input, required, pointer
  auto _ele = &ele; // input, required, pointer
  // orbit: CoordStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _orbit_desc;
  _orbit_desc.rank = 1;
  _orbit_desc.data_ptr = orbit.data();
  _orbit_desc.dims[0] = orbit.size();
  _orbit_desc.strides[0] = 1;
  fortran_tao_scratch_values_calc(
      /* void* */ &ele_ref,
      /* void* */ &ele_start,
      /* void* */ &ele,
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _orbit_desc
  );
}
void Tao::tao_set_beam_cmd(std::string who, std::string value_str, std::string branch_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  auto _branch_str = branch_str.c_str();
  fortran_tao_set_beam_cmd(
      /* const char* */ _who,
      /* const char* */ _value_str,
      /* const char* */ _branch_str
  );
}
void Tao::tao_set_beam_init_cmd(std::string who, std::string value_str, std::string branch_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  auto _branch_str = branch_str.c_str();
  fortran_tao_set_beam_init_cmd(
      /* const char* */ _who,
      /* const char* */ _value_str,
      /* const char* */ _branch_str
  );
}
void Tao::tao_set_bmad_com_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_bmad_com_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_branch_cmd(
    std::string branch_str,
    std::string component_str,
    std::string value_str
) {
  auto _branch_str = branch_str.c_str();
  auto _component_str = component_str.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_branch_cmd(
      /* const char* */ _branch_str,
      /* const char* */ _component_str,
      /* const char* */ _value_str
  );
}
void Tao::tao_set_calculate_cmd(std::optional<std::string> switch_) {
  const char *_switch_ = switch_.has_value() ? switch_->c_str() : nullptr;
  fortran_tao_set_calculate_cmd(/* const char* */ _switch_);
}
void Tao::tao_set_curve_cmd(std::string curve_name, std::string component, std::string value_str) {
  auto _curve_name = curve_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_curve_cmd(
      /* const char* */ _curve_name,
      /* const char* */ _component,
      /* const char* */ _value_str
  );
}
void Tao::tao_set_curve_invalid(
    TaoCurveStruct &curve,
    std::string why_invalid,
    std::optional<bool> print_err
) {
  auto _why_invalid = why_invalid.c_str();
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_set_curve_invalid(
      /* void* */ curve.get_fortran_ptr(),
      /* const char* */ _why_invalid,
      /* bool* */ _print_err
  );
}
void Tao::tao_set_data_cmd(std::string who_str, std::string value_str, std::optional<bool> silent) {
  auto _who_str = who_str.c_str();
  auto _value_str = value_str.c_str();
  bool silent_lvalue;
  auto *_silent{&silent_lvalue};
  if (silent.has_value()) {
    silent_lvalue = silent.value();
  } else {
    _silent = nullptr;
  }
  fortran_tao_set_data_cmd(
      /* const char* */ _who_str,
      /* const char* */ _value_str,
      /* bool* */ _silent
  );
}
void Tao::tao_set_data_useit_opt(std::optional<TaoDataStructArray1D> data) {
  // data: TaoDataStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _data_desc;
  _data_desc.rank = 1;
  if (data) {
    _data_desc.data_ptr = data->data();
    _data_desc.dims[0] = data->size();
  } else {
    _data_desc.data_ptr = nullptr;
    _data_desc.dims[0] = 0;
  }
  _data_desc.strides[0] = 1;
  fortran_tao_set_data_useit_opt(/* Bmad::array_descriptor_t& */ _data_desc);
}
void Tao::tao_set_default_cmd(std::string who_str, std::string value_str) {
  auto _who_str = who_str.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_default_cmd(/* const char* */ _who_str, /* const char* */ _value_str);
}
void Tao::tao_set_drawing_cmd(
    TaoDrawingStruct &drawing,
    std::string component,
    std::string value_str
) {
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_drawing_cmd(
      /* void* */ drawing.get_fortran_ptr(),
      /* const char* */ _component,
      /* const char* */ _value_str
  );
}
void Tao::tao_set_dynamic_aperture_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_dynamic_aperture_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_elements_cmd(
    std::string ele_list,
    std::string attribute,
    std::string value,
    bool update
) {
  auto _ele_list = ele_list.c_str();
  auto _attribute = attribute.c_str();
  auto _value = value.c_str();
  fortran_tao_set_elements_cmd(
      /* const char* */ _ele_list,
      /* const char* */ _attribute,
      /* const char* */ _value,
      /* bool& */ update
  );
}
void Tao::tao_set_floor_plan_axis_label(
    TaoGraphStruct &graph,
    QpAxisStruct &axis_in,
    QpAxisStruct &axis_out,
    std::string which
) {
  auto _which = which.c_str();
  fortran_tao_set_floor_plan_axis_label(
      /* void* */ graph.get_fortran_ptr(),
      /* void* */ axis_in.get_fortran_ptr(),
      /* void* */ axis_out.get_fortran_ptr(),
      /* const char* */ _which
  );
}
void Tao::tao_set_geodesic_lm_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_geodesic_lm_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_global_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_global_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_graph_cmd(std::string graph_name, std::string component, std::string value_str) {
  auto _graph_name = graph_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_graph_cmd(
      /* const char* */ _graph_name,
      /* const char* */ _component,
      /* const char* */ _value_str
  );
}
Tao::TaoSetIntegerValue Tao::tao_set_integer_value(
    std::string var_str,
    std::string value_str,
    std::optional<int> min_val,
    std::optional<int> max_val,
    std::optional<bool> print_err
) {
  int _var{};
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  bool _error{};
  int min_val_lvalue;
  auto *_min_val{&min_val_lvalue};
  if (min_val.has_value()) {
    min_val_lvalue = min_val.value();
  } else {
    _min_val = nullptr;
  }
  int max_val_lvalue;
  auto *_max_val{&max_val_lvalue};
  if (max_val.has_value()) {
    max_val_lvalue = max_val.value();
  } else {
    _max_val = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_set_integer_value(
      /* int& */ _var,
      /* const char* */ _var_str,
      /* const char* */ _value_str,
      /* bool& */ _error,
      /* int* */ _min_val,
      /* int* */ _max_val,
      /* bool* */ _print_err
  );
  return TaoSetIntegerValue{_var, _error};
}
std::string Tao::tao_set_invalid(
    TaoDataStruct &datum,
    std::string message,
    std::optional<bool> exterminate,
    std::optional<int> err_level,
    std::optional<bool> print_err
) {
  auto _message = message.c_str();
  char _why_invalid[4096];
  bool exterminate_lvalue;
  auto *_exterminate{&exterminate_lvalue};
  if (exterminate.has_value()) {
    exterminate_lvalue = exterminate.value();
  } else {
    _exterminate = nullptr;
  }
  int err_level_lvalue;
  auto *_err_level{&err_level_lvalue};
  if (err_level.has_value()) {
    err_level_lvalue = err_level.value();
  } else {
    _err_level = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_set_invalid(
      /* void* */ datum.get_fortran_ptr(),
      /* const char* */ _message,
      /* const char* */ _why_invalid,
      /* bool* */ _exterminate,
      /* int* */ _err_level,
      /* bool* */ _print_err
  );
  return _why_invalid;
}
void Tao::tao_set_key_cmd(std::string key_str, std::string cmd_str) {
  auto _key_str = key_str.c_str();
  auto _cmd_str = cmd_str.c_str();
  fortran_tao_set_key_cmd(/* const char* */ _key_str, /* const char* */ _cmd_str);
}
void Tao::tao_set_lattice_cmd(std::string dest_lat, std::string source_lat) {
  auto _dest_lat = dest_lat.c_str();
  auto _source_lat = source_lat.c_str();
  fortran_tao_set_lattice_cmd(/* const char* */ _dest_lat, /* const char* */ _source_lat);
}
Tao::TaoSetLogicalValue Tao::tao_set_logical_value(std::string var_str, std::string value_str) {
  bool _var{};
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  bool _error{};
  fortran_tao_set_logical_value(
      /* bool& */ _var,
      /* const char* */ _var_str,
      /* const char* */ _value_str,
      /* bool& */ _error
  );
  return TaoSetLogicalValue{_var, _error};
}
void Tao::tao_set_openmp_n_threads(int n_threads) {
  fortran_tao_set_openmp_n_threads(/* int& */ n_threads);
}
void Tao::tao_set_opt_vars(FArray1D<Real> &var_vec, std::optional<bool> print_limit_warning) {
  // var_vec: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _var_vec_desc;
  _var_vec_desc.rank = 1;
  _var_vec_desc.data_ptr = var_vec.data();
  _var_vec_desc.dims[0] = var_vec.size();
  bool print_limit_warning_lvalue;
  auto *_print_limit_warning{&print_limit_warning_lvalue};
  if (print_limit_warning.has_value()) {
    print_limit_warning_lvalue = print_limit_warning.value();
  } else {
    _print_limit_warning = nullptr;
  }
  fortran_tao_set_opt_vars(
      /* Bmad::array_descriptor_t& */ _var_vec_desc,
      /* bool* */ _print_limit_warning
  );
}
void Tao::tao_set_opti_de_param_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_opti_de_param_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_particle_start_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_particle_start_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_plot_cmd(std::string plot_name, std::string component, std::string value_str) {
  auto _plot_name = plot_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_plot_cmd(
      /* const char* */ _plot_name,
      /* const char* */ _component,
      /* const char* */ _value_str
  );
}
void Tao::tao_set_plot_page_cmd(
    std::string component,
    std::string value_str,
    std::optional<std::string> value_str2
) {
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  const char *_value_str2 = value_str2.has_value() ? value_str2->c_str() : nullptr;
  fortran_tao_set_plot_page_cmd(
      /* const char* */ _component,
      /* const char* */ _value_str,
      /* const char* */ _value_str2
  );
}
void Tao::tao_set_ptc_com_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_ptc_com_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
Tao::TaoSetQpAxisStruct Tao::tao_set_qp_axis_struct(
    std::string qp_axis_name,
    std::string component,
    QpAxisStruct &qp_axis,
    std::string value
) {
  auto _qp_axis_name = qp_axis_name.c_str();
  auto _component = component.c_str();
  auto _value = value.c_str();
  bool _error{};
  int _ix_uni{};
  fortran_tao_set_qp_axis_struct(
      /* const char* */ _qp_axis_name,
      /* const char* */ _component,
      /* void* */ qp_axis.get_fortran_ptr(),
      /* const char* */ _value,
      /* bool& */ _error,
      /* int& */ _ix_uni
  );
  return TaoSetQpAxisStruct{_error, _ix_uni};
}
Tao::TaoSetQpPointStruct Tao::tao_set_qp_point_struct(
    std::string qp_point_name,
    std::string component,
    QpPointStruct &qp_point,
    std::string value
) {
  auto _qp_point_name = qp_point_name.c_str();
  auto _component = component.c_str();
  auto _value = value.c_str();
  bool _error{};
  int _ix_uni{};
  fortran_tao_set_qp_point_struct(
      /* const char* */ _qp_point_name,
      /* const char* */ _component,
      /* void* */ qp_point.get_fortran_ptr(),
      /* const char* */ _value,
      /* bool& */ _error,
      /* int& */ _ix_uni
  );
  return TaoSetQpPointStruct{_error, _ix_uni};
}
Tao::TaoSetQpRectStruct Tao::tao_set_qp_rect_struct(
    std::string qp_rect_name,
    std::string component,
    QpRectStruct &qp_rect,
    std::string value
) {
  auto _qp_rect_name = qp_rect_name.c_str();
  auto _component = component.c_str();
  auto _value = value.c_str();
  bool _error{};
  int _ix_uni{};
  fortran_tao_set_qp_rect_struct(
      /* const char* */ _qp_rect_name,
      /* const char* */ _component,
      /* void* */ qp_rect.get_fortran_ptr(),
      /* const char* */ _value,
      /* bool& */ _error,
      /* int& */ _ix_uni
  );
  return TaoSetQpRectStruct{_error, _ix_uni};
}
void Tao::tao_set_ran_state_cmd(std::string state_string) {
  auto _state_string = state_string.c_str();
  fortran_tao_set_ran_state_cmd(/* const char* */ _state_string);
}
Tao::TaoSetRealValue Tao::tao_set_real_value(
    std::string var_str,
    std::string value_str,
    std::optional<double> min_val,
    std::optional<double> max_val,
    std::optional<int> dflt_uni
) {
  double _var{};
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  bool _error{};
  double min_val_lvalue;
  auto *_min_val{&min_val_lvalue};
  if (min_val.has_value()) {
    min_val_lvalue = min_val.value();
  } else {
    _min_val = nullptr;
  }
  double max_val_lvalue;
  auto *_max_val{&max_val_lvalue};
  if (max_val.has_value()) {
    max_val_lvalue = max_val.value();
  } else {
    _max_val = nullptr;
  }
  int dflt_uni_lvalue;
  auto *_dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  fortran_tao_set_real_value(
      /* double& */ _var,
      /* const char* */ _var_str,
      /* const char* */ _value_str,
      /* bool& */ _error,
      /* double* */ _min_val,
      /* double* */ _max_val,
      /* int* */ _dflt_uni
  );
  return TaoSetRealValue{_var, _error};
}
void Tao::tao_set_region_cmd(
    std::string region_name,
    std::string component,
    std::string value_str
) {
  auto _region_name = region_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_region_cmd(
      /* const char* */ _region_name,
      /* const char* */ _component,
      /* const char* */ _value_str
  );
}
void Tao::tao_set_space_charge_com_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_space_charge_com_cmd(/* const char* */ _who, /* const char* */ _value_str);
}
void Tao::tao_set_symbolic_number_cmd(
    std::string sym_str,
    std::optional<std::string> num_str,
    std::optional<double> val
) {
  auto _sym_str = sym_str.c_str();
  const char *_num_str = num_str.has_value() ? num_str->c_str() : nullptr;
  double val_lvalue;
  auto *_val{&val_lvalue};
  if (val.has_value()) {
    val_lvalue = val.value();
  } else {
    _val = nullptr;
  }
  fortran_tao_set_symbolic_number_cmd(
      /* const char* */ _sym_str,
      /* const char* */ _num_str,
      /* double* */ _val
  );
}
void Tao::tao_set_tune_cmd(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string qa_str,
    std::string qb_str,
    bool delta_input
) {
  auto _branch_str = branch_str.c_str();
  auto _mask_str = mask_str.c_str();
  auto _qa_str = qa_str.c_str();
  auto _qb_str = qb_str.c_str();
  fortran_tao_set_tune_cmd(
      /* const char* */ _branch_str,
      /* const char* */ _mask_str,
      /* bool& */ print_list,
      /* const char* */ _qa_str,
      /* const char* */ _qb_str,
      /* bool& */ delta_input
  );
}
void Tao::tao_set_universe_cmd(std::string uni, std::string who, std::string what) {
  auto _uni = uni.c_str();
  auto _who = who.c_str();
  auto _what = what.c_str();
  fortran_tao_set_universe_cmd(
      /* const char* */ _uni,
      /* const char* */ _who,
      /* const char* */ _what
  );
}
void Tao::tao_set_var_cmd(std::string var_str, std::string value_str) {
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_var_cmd(/* const char* */ _var_str, /* const char* */ _value_str);
}
void Tao::tao_set_var_model_value(
    TaoVarStruct &var,
    double value,
    std::optional<bool> print_limit_warning
) {
  bool print_limit_warning_lvalue;
  auto *_print_limit_warning{&print_limit_warning_lvalue};
  if (print_limit_warning.has_value()) {
    print_limit_warning_lvalue = print_limit_warning.value();
  } else {
    _print_limit_warning = nullptr;
  }
  fortran_tao_set_var_model_value(
      /* void* */ var.get_fortran_ptr(),
      /* double& */ value,
      /* bool* */ _print_limit_warning
  );
}
void Tao::tao_set_var_useit_opt() { fortran_tao_set_var_useit_opt(); }
bool Tao::tao_set_wave_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  bool _err{};
  fortran_tao_set_wave_cmd(/* const char* */ _who, /* const char* */ _value_str, /* bool& */ _err);
  return _err;
}
void Tao::tao_set_z_tune_cmd(std::string branch_str, std::string q_str, bool delta_input) {
  auto _branch_str = branch_str.c_str();
  auto _q_str = q_str.c_str();
  fortran_tao_set_z_tune_cmd(
      /* const char* */ _branch_str,
      /* const char* */ _q_str,
      /* bool& */ delta_input
  );
}
void Tao::tao_setup_key_table() { fortran_tao_setup_key_table(); }
bool Tao::tao_shape_init(TaoEleShapeStruct &shape, std::optional<bool> print_err) {
  bool _err{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_shape_init(
      /* void* */ shape.get_fortran_ptr(),
      /* bool& */ _err,
      /* bool* */ _print_err
  );
  return _err;
}
void Tao::tao_show_cmd(std::string what) {
  auto _what = what.c_str();
  fortran_tao_show_cmd(/* const char* */ _what);
}
void Tao::tao_show_constraints(int iunit, std::string form) {
  auto _form = form.c_str();
  fortran_tao_show_constraints(/* int& */ iunit, /* const char* */ _form);
}
void Tao::tao_single_mode(std::string char_) {
  auto _char_ = char_.c_str();
  fortran_tao_single_mode(/* const char* */ _char_);
}
bool Tao::tao_single_track(
    TaoLatticeStruct &tao_lat,
    int ix_branch,
    std::optional<bool> print_err
) {
  bool _calc_ok{};
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_single_track(
      /* void* */ tao_lat.get_fortran_ptr(),
      /* bool& */ _calc_ok,
      /* int& */ ix_branch,
      /* bool* */ _print_err
  );
  return _calc_ok;
}
void Tao::tao_spin_matrices_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_calc
) {
  auto _data_type = data_type.c_str();
  auto _data_source = data_source.c_str();
  fortran_tao_spin_matrices_calc_needed(
      /* const char* */ _data_type,
      /* const char* */ _data_source,
      /* bool& */ do_calc
  );
}
void Tao::tao_spin_tracking_turn_on() { fortran_tao_spin_tracking_turn_on(); }
bool Tao::tao_split_component(std::string comp_str, TaoDataVarComponentStructAlloc1D comp) {
  auto _comp_str = comp_str.c_str();
  // intent=inout allocatable type array
  bool _err{};
  fortran_tao_split_component(
      /* const char* */ _comp_str,
      /* void* */ comp.get_fortran_ptr(),
      /* bool& */ _err
  );
  return _err;
}
void Tao::tao_srdt_calc_needed(std::string data_type, std::string data_source, int do_srdt) {
  auto _data_type = data_type.c_str();
  auto _data_source = data_source.c_str();
  fortran_tao_srdt_calc_needed(
      /* const char* */ _data_type,
      /* const char* */ _data_source,
      /* int& */ do_srdt
  );
}
Tao::TaoSubinUniNumber Tao::tao_subin_uni_number(std::string name_in, int ix_uni) {
  auto _name_in = name_in.c_str();
  char _name_out[4096];
  bool _ok{};
  fortran_tao_subin_uni_number(
      /* const char* */ _name_in,
      /* int& */ ix_uni,
      /* const char* */ _name_out,
      /* bool& */ _ok
  );
  return TaoSubinUniNumber{_name_out, _ok};
}
bool Tao::tao_svd_optimizer() {
  bool _abort{};
  fortran_tao_svd_optimizer(/* bool& */ _abort);
  return _abort;
}
void Tao::tao_symbol_import_from_lat(LatStruct &lat) {
  fortran_tao_symbol_import_from_lat(/* void* */ lat.get_fortran_ptr());
}
void Tao::tao_taper_cmd(std::string except, std::string uni_names) {
  auto _except = except.c_str();
  auto _uni_names = uni_names.c_str();
  fortran_tao_taper_cmd(/* const char* */ _except, /* const char* */ _uni_names);
}
void Tao::tao_to_change_number(
    std::string num_str,
    int n_size,
    RealAlloc1D &change_number,
    std::string abs_or_rel,
    bool err
) {
  auto _num_str = num_str.c_str();
  // intent=inout allocatable general array
  auto _abs_or_rel = abs_or_rel.c_str();
  fortran_tao_to_change_number(
      /* const char* */ _num_str,
      /* int& */ n_size,
      /* void* */ change_number.get_fortran_ptr(),
      /* const char* */ _abs_or_rel,
      /* bool& */ err
  );
}
void Tao::tao_to_int(std::string str, int i_int, bool err) {
  auto _str = str.c_str();
  fortran_tao_to_int(/* const char* */ _str, /* int& */ i_int, /* bool& */ err);
}
Tao::TaoToPhaseAndCouplingReading Tao::tao_to_phase_and_coupling_reading(
    EleStruct &ele,
    std::string why_invalid,
    TaoDataStruct &datum
) {
  BpmPhaseCouplingStruct _bpm_data;
  bool _valid_value{};
  auto _why_invalid = why_invalid.c_str();
  fortran_tao_to_phase_and_coupling_reading(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _bpm_data.get_fortran_ptr(),
      /* bool& */ _valid_value,
      /* const char* */ _why_invalid,
      /* void* */ datum.get_fortran_ptr()
  );
  return TaoToPhaseAndCouplingReading{std::move(_bpm_data), _valid_value};
}
Tao::TaoToReal Tao::tao_to_real(std::string expression) {
  auto _expression = expression.c_str();
  double _value{};
  bool _err_flag{};
  fortran_tao_to_real(/* const char* */ _expression, /* double& */ _value, /* bool& */ _err_flag);
  return TaoToReal{_value, _err_flag};
}
void Tao::tao_too_many_particles_lost(BeamStruct &beam, bool no_beam) {
  fortran_tao_too_many_particles_lost(/* void* */ beam.get_fortran_ptr(), /* bool& */ no_beam);
}
void Tao::tao_top10_derivative_print() { fortran_tao_top10_derivative_print(); }
void Tao::tao_top10_merit_categories_print(int iunit) {
  fortran_tao_top10_merit_categories_print(/* int& */ iunit);
}
int Tao::tao_top_level(std::optional<std::string> command) {
  const char *_command = command.has_value() ? command->c_str() : nullptr;
  int _errcode{};
  fortran_tao_top_level(/* const char* */ _command, /* int& */ _errcode);
  return _errcode;
}
Tao::TaoTrackingEleIndex Tao::tao_tracking_ele_index(EleStruct &ele, TaoDataStruct &datum) {
  auto _ele = &ele; // input, required, pointer
  int _ix_branch{};
  int _ix_ele{};
  fortran_tao_tracking_ele_index(
      /* void* */ &ele,
      /* void* */ datum.get_fortran_ptr(),
      /* int& */ _ix_branch,
      /* int& */ _ix_ele
  );
  return TaoTrackingEleIndex{_ix_branch, _ix_ele};
}
void Tao::tao_turn_on_special_calcs_if_needed_for_plotting() {
  fortran_tao_turn_on_special_calcs_if_needed_for_plotting();
}
void Tao::tao_type_expression_tree(TaoEvalNodeStruct &tree, std::optional<int> indent) {
  int indent_lvalue;
  auto *_indent{&indent_lvalue};
  if (indent.has_value()) {
    indent_lvalue = indent.value();
  } else {
    _indent = nullptr;
  }
  fortran_tao_type_expression_tree(/* void* */ tree.get_fortran_ptr(), /* int* */ _indent);
}
int Tao::tao_uni_atsign_index(std::string string) {
  auto _string = string.c_str();
  int _ix_amp{};
  fortran_tao_uni_atsign_index(/* const char* */ _string, /* int& */ _ix_amp);
  return _ix_amp;
}
int Tao::tao_universe_index(int i_uni, std::optional<bool> neg2_to_default) {
  bool neg2_to_default_lvalue;
  auto *_neg2_to_default{&neg2_to_default_lvalue};
  if (neg2_to_default.has_value()) {
    neg2_to_default_lvalue = neg2_to_default.value();
  } else {
    _neg2_to_default = nullptr;
  }
  int _i_this_uni{};
  fortran_tao_universe_index(
      /* int& */ i_uni,
      /* bool* */ _neg2_to_default,
      /* int& */ _i_this_uni
  );
  return _i_this_uni;
}
void Tao::tao_use_data(std::string action, std::string data_name) {
  auto _action = action.c_str();
  auto _data_name = data_name.c_str();
  fortran_tao_use_data(/* const char* */ _action, /* const char* */ _data_name);
}
void Tao::tao_use_var(std::string action, std::string var_name) {
  auto _action = action.c_str();
  auto _var_name = var_name.c_str();
  fortran_tao_use_var(/* const char* */ _action, /* const char* */ _var_name);
}
bool Tao::tao_user_is_terminating_optimization() {
  bool _is_terminating{};
  fortran_tao_user_is_terminating_optimization(/* bool& */ _is_terminating);
  return _is_terminating;
}
std::string Tao::tao_var1_name(TaoVarStruct &var) {
  char _var1_name[4096];
  fortran_tao_var1_name(/* void* */ var.get_fortran_ptr(), /* const char* */ _var1_name);
  return _var1_name;
}
std::string Tao::tao_var_attrib_name(TaoVarStruct &var) {
  char _var_attrib_name[4096];
  fortran_tao_var_attrib_name(
      /* void* */ var.get_fortran_ptr(),
      /* const char* */ _var_attrib_name
  );
  return _var_attrib_name;
}
void Tao::tao_var_check(ElePointerStructAlloc1D eles, std::string attribute, bool silent) {
  // intent=in allocatable type array
  auto _attribute = attribute.c_str();
  fortran_tao_var_check(
      /* void* */ eles.get_fortran_ptr(),
      /* const char* */ _attribute,
      /* bool& */ silent
  );
}
void Tao::tao_var_repoint() { fortran_tao_var_repoint(); }
void Tao::tao_var_target_calc() { fortran_tao_var_target_calc(); }
void Tao::tao_var_useit_plot_calc(TaoGraphStruct &graph, TaoVarStructArray1D var) {
  // var: TaoVarStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _var_desc;
  _var_desc.rank = 1;
  _var_desc.data_ptr = var.data();
  _var_desc.dims[0] = var.size();
  _var_desc.strides[0] = 1;
  fortran_tao_var_useit_plot_calc(
      /* void* */ graph.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _var_desc
  );
}
void Tao::tao_var_write(
    std::string out_file,
    std::optional<bool> show_good_opt_only,
    std::optional<bool> tao_format
) {
  auto _out_file = out_file.c_str();
  bool show_good_opt_only_lvalue;
  auto *_show_good_opt_only{&show_good_opt_only_lvalue};
  if (show_good_opt_only.has_value()) {
    show_good_opt_only_lvalue = show_good_opt_only.value();
  } else {
    _show_good_opt_only = nullptr;
  }
  bool tao_format_lvalue;
  auto *_tao_format{&tao_format_lvalue};
  if (tao_format.has_value()) {
    tao_format_lvalue = tao_format.value();
  } else {
    _tao_format = nullptr;
  }
  fortran_tao_var_write(
      /* const char* */ _out_file,
      /* bool* */ _show_good_opt_only,
      /* bool* */ _tao_format
  );
}
void Tao::tao_veto_vars_with_zero_dmodel() { fortran_tao_veto_vars_with_zero_dmodel(); }
void Tao::tao_wave_analysis(TaoPlotStruct &plot) {
  fortran_tao_wave_analysis(/* void* */ plot.get_fortran_ptr());
}
void Tao::tao_wave_cmd(std::string curve_name, std::string plot_place, bool err_flag) {
  auto _curve_name = curve_name.c_str();
  auto _plot_place = plot_place.c_str();
  fortran_tao_wave_cmd(
      /* const char* */ _curve_name,
      /* const char* */ _plot_place,
      /* bool& */ err_flag
  );
}
void Tao::tao_wave_fit(
    TaoCurveStruct &curve,
    int ix1,
    int n_dat,
    FArray1D<Real> &coef,
    FArray1D<Real> &rms,
    FArray1D<Real> &f1,
    std::optional<FArray1D<Real>> f2,
    std::optional<FArray1D<Real>> f3,
    std::optional<FArray1D<Real>> f4
) {
  // coef: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _coef_desc;
  _coef_desc.rank = 1;
  _coef_desc.data_ptr = coef.data();
  _coef_desc.dims[0] = coef.size();
  // rms: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _rms_desc;
  _rms_desc.rank = 1;
  _rms_desc.data_ptr = rms.data();
  _rms_desc.dims[0] = rms.size();
  // f1: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _f1_desc;
  _f1_desc.rank = 1;
  _f1_desc.data_ptr = f1.data();
  _f1_desc.dims[0] = f1.size();
  // f2: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _f2_desc;
  _f2_desc.rank = 1;
  if (f2.has_value()) {
    _f2_desc.data_ptr = f2->data();
    _f2_desc.dims[0] = f2->size();
  } else {
    _f2_desc.data_ptr = nullptr;
    _f2_desc.dims[0] = 0;
  }
  // f3: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _f3_desc;
  _f3_desc.rank = 1;
  if (f3.has_value()) {
    _f3_desc.data_ptr = f3->data();
    _f3_desc.dims[0] = f3->size();
  } else {
    _f3_desc.data_ptr = nullptr;
    _f3_desc.dims[0] = 0;
  }
  // f4: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _f4_desc;
  _f4_desc.rank = 1;
  if (f4.has_value()) {
    _f4_desc.data_ptr = f4->data();
    _f4_desc.dims[0] = f4->size();
  } else {
    _f4_desc.data_ptr = nullptr;
    _f4_desc.dims[0] = 0;
  }
  fortran_tao_wave_fit(
      /* void* */ curve.get_fortran_ptr(),
      /* int& */ ix1,
      /* int& */ n_dat,
      /* Bmad::array_descriptor_t& */ _coef_desc,
      /* Bmad::array_descriptor_t& */ _rms_desc,
      /* Bmad::array_descriptor_t& */ _f1_desc,
      /* Bmad::array_descriptor_t& */ _f2_desc,
      /* Bmad::array_descriptor_t& */ _f3_desc,
      /* Bmad::array_descriptor_t& */ _f4_desc
  );
}
void Tao::tao_write_cmd(std::string what) {
  auto _what = what.c_str();
  fortran_tao_write_cmd(/* const char* */ _what);
}
void Tao::tao_x_axis_cmd(std::string where, std::string what) {
  auto _where = where.c_str();
  auto _what = what.c_str();
  fortran_tao_x_axis_cmd(/* const char* */ _where, /* const char* */ _what);
}
bool Tao::tao_x_scale_cmd(
    std::string where,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall,
    std::optional<std::string> gang,
    std::optional<bool> exact,
    std::optional<bool> turn_autoscale_off
) {
  auto _where = where.c_str();
  bool _err{};
  bool include_wall_lvalue;
  auto *_include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  const char *_gang = gang.has_value() ? gang->c_str() : nullptr;
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  bool turn_autoscale_off_lvalue;
  auto *_turn_autoscale_off{&turn_autoscale_off_lvalue};
  if (turn_autoscale_off.has_value()) {
    turn_autoscale_off_lvalue = turn_autoscale_off.value();
  } else {
    _turn_autoscale_off = nullptr;
  }
  fortran_tao_x_scale_cmd(
      /* const char* */ _where,
      /* double& */ x_min_in,
      /* double& */ x_max_in,
      /* bool& */ _err,
      /* bool* */ _include_wall,
      /* const char* */ _gang,
      /* bool* */ _exact,
      /* bool* */ _turn_autoscale_off
  );
  return _err;
}
void Tao::tao_x_scale_graph(
    TaoGraphStruct &graph,
    double x_min,
    double x_max,
    std::optional<bool> include_wall,
    std::optional<bool> have_scaled
) {
  bool include_wall_lvalue;
  auto *_include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  bool have_scaled_lvalue;
  auto *_have_scaled{&have_scaled_lvalue};
  if (have_scaled.has_value()) {
    have_scaled_lvalue = have_scaled.value();
  } else {
    _have_scaled = nullptr;
  }
  fortran_tao_x_scale_graph(
      /* void* */ graph.get_fortran_ptr(),
      /* double& */ x_min,
      /* double& */ x_max,
      /* bool* */ _include_wall,
      /* bool* */ _have_scaled
  );
}
bool Tao::tao_x_scale_plot(
    TaoPlotStruct &plot,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall,
    std::optional<std::string> gang
) {
  bool include_wall_lvalue;
  auto *_include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  const char *_gang = gang.has_value() ? gang->c_str() : nullptr;
  bool _have_scaled{};
  fortran_tao_x_scale_plot(
      /* void* */ plot.get_fortran_ptr(),
      /* double& */ x_min_in,
      /* double& */ x_max_in,
      /* bool* */ _include_wall,
      /* const char* */ _gang,
      /* bool& */ _have_scaled
  );
  return _have_scaled;
}