#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/tao_routines.hpp"
#include "bmad/types.h"
#include "json.hpp"

using namespace Bmad;

using json = nlohmann::json;
void Tao::tao_deallocate_plot_cache(TaoPlotCacheProxyAlloc1D& plot_cache) {
  // intent=inout allocatable type array
  fortran_tao_deallocate_plot_cache(/* void* */ plot_cache.get_fortran_ptr());
}
void Tao::tao_lattice_branches_equal_tao_lattice_branches(
    TaoLatticeBranchProxyAlloc1D& tlb1,
    TaoLatticeBranchProxyAlloc1D& tlb2) {
  // intent=inout allocatable type array
  // intent=inout allocatable type array
  fortran_tao_lattice_branches_equal_tao_lattice_branches(
      /* void* */ tlb1.get_fortran_ptr(), /* void* */ tlb2.get_fortran_ptr());
}
void Tao::tao_lattice_equal_tao_lattice(
    TaoLatticeProxy& lat1,
    TaoLatticeProxy& lat2) {
  fortran_tao_lattice_equal_tao_lattice(
      /* void* */ lat1.get_fortran_ptr(), /* void* */ lat2.get_fortran_ptr());
}
void Tao::tao_turn_on_special_calcs_if_needed_for_plotting() {
  fortran_tao_turn_on_special_calcs_if_needed_for_plotting();
}
void Tao::tao_pointer_to_tao_lat(
    TaoUniverseProxy& u,
    std::optional<int> lat_type,
    TaoLatticeProxy& tao_lat) {
  int lat_type_lvalue;
  auto* _lat_type{&lat_type_lvalue};
  if (lat_type.has_value()) {
    lat_type_lvalue = lat_type.value();
  } else {
    _lat_type = nullptr;
  }
  auto _tao_lat = &tao_lat; // input, required, pointer
  fortran_tao_pointer_to_tao_lat(
      /* void* */ u.get_fortran_ptr(),
      /* c_Int* */ _lat_type,
      /* void* */ &tao_lat);
}
bool Tao::tao_de_optimizer() {
  bool _abort{};
  fortran_tao_de_optimizer(/* c_Bool& */ _abort);
  return _abort;
}
bool Tao::tao_geodesic_lm_optimizer() {
  bool _abort{};
  fortran_tao_geodesic_lm_optimizer(/* c_Bool& */ _abort);
  return _abort;
}
void Tao::tao_branch_index(int ix_branch, int ix_this) {
  fortran_tao_branch_index(/* c_Int& */ ix_branch, /* c_Int& */ ix_this);
}
bool Tao::tao_limit_calc() {
  bool _limited{};
  fortran_tao_limit_calc(/* c_Bool& */ _limited);
  return _limited;
}
void Tao::tao_alias_cmd(std::string alias, std::string string) {
  auto _alias = alias.c_str();
  auto _string = string.c_str();
  fortran_tao_alias_cmd(/* c_Char */ _alias, /* c_Char */ _string);
}
int Tao::tao_top_level(std::optional<std::string> command) {
  const char* _command = command.has_value() ? command->c_str() : nullptr;
  int _errcode{};
  fortran_tao_top_level(/* c_Char */ _command, /* c_Int& */ _errcode);
  return _errcode;
}
void Tao::tao_universe_index(
    int i_uni,
    std::optional<bool> neg2_to_default,
    int i_this_uni) {
  bool neg2_to_default_lvalue;
  auto* _neg2_to_default{&neg2_to_default_lvalue};
  if (neg2_to_default.has_value()) {
    neg2_to_default_lvalue = neg2_to_default.value();
  } else {
    _neg2_to_default = nullptr;
  }
  fortran_tao_universe_index(
      /* c_Int& */ i_uni,
      /* c_Bool* */ _neg2_to_default,
      /* c_Int& */ i_this_uni);
}
Tao::TaoToReal Tao::tao_to_real(std::string expression) {
  auto _expression = expression.c_str();
  double _value{};
  bool _err_flag{};
  fortran_tao_to_real(
      /* c_Char */ _expression, /* c_Real& */ _value, /* c_Bool& */ _err_flag);
  return TaoToReal{_value, _err_flag};
}
void Tao::tao_data_sanity_check(
    TaoDataProxy& datum,
    bool print_err,
    std::string default_data_type,
    optional_ref<TaoUniverseProxy> uni,
    bool is_valid) {
  auto _default_data_type = default_data_type.c_str();
  auto* _uni = uni.has_value() ? uni->get().get_fortran_ptr()
                               : nullptr; // input, optional
  fortran_tao_data_sanity_check(
      /* void* */ datum.get_fortran_ptr(),
      /* c_Bool& */ print_err,
      /* c_Char */ _default_data_type,
      /* void* */ _uni,
      /* c_Bool& */ is_valid);
}
bool Tao::tao_command(std::string command_line, bool err) {
  auto _command_line = command_line.c_str();
  bool _err_is_fatal{};
  fortran_tao_command(
      /* c_Char */ _command_line,
      /* c_Bool& */ err,
      /* c_Bool& */ _err_is_fatal);
  return _err_is_fatal;
}
void Tao::tao_scale_cmd(
    std::string where,
    double y_min_in,
    double y_max_in,
    std::optional<std::string> axis,
    std::optional<bool> include_wall,
    std::optional<std::string> gang,
    std::optional<bool> exact,
    std::optional<bool> turn_autoscale_off) {
  auto _where = where.c_str();
  const char* _axis = axis.has_value() ? axis->c_str() : nullptr;
  bool include_wall_lvalue;
  auto* _include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  const char* _gang = gang.has_value() ? gang->c_str() : nullptr;
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  bool turn_autoscale_off_lvalue;
  auto* _turn_autoscale_off{&turn_autoscale_off_lvalue};
  if (turn_autoscale_off.has_value()) {
    turn_autoscale_off_lvalue = turn_autoscale_off.value();
  } else {
    _turn_autoscale_off = nullptr;
  }
  fortran_tao_scale_cmd(
      /* c_Char */ _where,
      /* c_Real& */ y_min_in,
      /* c_Real& */ y_max_in,
      /* c_Char */ _axis,
      /* c_Bool* */ _include_wall,
      /* c_Char */ _gang,
      /* c_Bool* */ _exact,
      /* c_Bool* */ _turn_autoscale_off);
}
void Tao::tao_pointer_to_datum(
    TaoD1DataProxy& d1,
    std::string ele_name,
    TaoDataProxy& datum_ptr) {
  auto _ele_name = ele_name.c_str();
  auto _datum_ptr = &datum_ptr; // input, required, pointer
  fortran_tao_pointer_to_datum(
      /* void* */ d1.get_fortran_ptr(),
      /* c_Char */ _ele_name,
      /* void* */ &datum_ptr);
}
void Tao::tao_lat_sigma_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_lat_sigma) {
  auto _data_type = data_type.c_str(); // ptr, inout, required
  auto _data_source = data_source.c_str(); // ptr, inout, required
  fortran_tao_lat_sigma_calc_needed(
      /* c_Char */ _data_type,
      /* c_Char */ _data_source,
      /* c_Bool& */ do_lat_sigma);
}
bool Tao::tao_single_track(
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<bool> print_err) {
  bool _calc_ok{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_single_track(
      /* void* */ tao_lat.get_fortran_ptr(),
      /* c_Bool& */ _calc_ok,
      /* c_Int& */ ix_branch,
      /* c_Bool* */ _print_err);
  return _calc_ok;
}
bool Tao::tao_lat_sigma_track(
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<bool> print_err,
    std::optional<bool> force_calc) {
  bool _calc_ok{};
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool force_calc_lvalue;
  auto* _force_calc{&force_calc_lvalue};
  if (force_calc.has_value()) {
    force_calc_lvalue = force_calc.value();
  } else {
    _force_calc = nullptr;
  }
  fortran_tao_lat_sigma_track(
      /* void* */ tao_lat.get_fortran_ptr(),
      /* c_Bool& */ _calc_ok,
      /* c_Int& */ ix_branch,
      /* c_Bool* */ _print_err,
      /* c_Bool* */ _force_calc);
  return _calc_ok;
}
bool Tao::tao_beam_track(
    TaoUniverseProxy& u,
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    BeamProxy& beam) {
  bool _calc_ok{};
  fortran_tao_beam_track(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* c_Int& */ ix_branch,
      /* void* */ beam.get_fortran_ptr(),
      /* c_Bool& */ _calc_ok);
  return _calc_ok;
}
void Tao::tao_too_many_particles_lost(BeamProxy& beam, bool no_beam) {
  fortran_tao_too_many_particles_lost(
      /* void* */ beam.get_fortran_ptr(), /* c_Bool& */ no_beam);
}
void Tao::tao_inject_particle(
    TaoUniverseProxy& u,
    TaoLatticeProxy& model,
    int ix_branch) {
  fortran_tao_inject_particle(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ model.get_fortran_ptr(),
      /* c_Int& */ ix_branch);
}
Tao::TaoInjectBeam Tao::tao_inject_beam(
    TaoUniverseProxy& u,
    TaoLatticeProxy& model,
    int ix_branch) {
  BeamProxy _beam;
  bool _init_ok{};
  fortran_tao_inject_beam(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ model.get_fortran_ptr(),
      /* c_Int& */ ix_branch,
      /* void* */ _beam.get_fortran_ptr(),
      /* c_Bool& */ _init_ok);
  return TaoInjectBeam{std::move(_beam), _init_ok};
}
void Tao::tao_beam_track_endpoint(
    std::string ele_id,
    LatProxy& lat,
    std::string branch_str,
    std::string where,
    TaoUniverseProxy& u,
    EleProxy& ele) {
  auto _ele_id = ele_id.c_str();
  auto _branch_str = branch_str.c_str();
  auto _where = where.c_str();
  auto _ele = &ele; // input, required, pointer
  fortran_tao_beam_track_endpoint(
      /* c_Char */ _ele_id,
      /* void* */ lat.get_fortran_ptr(),
      /* c_Char */ _branch_str,
      /* c_Char */ _where,
      /* void* */ u.get_fortran_ptr(),
      /* void* */ &ele);
}
int Tao::tao_open_file(
    std::string file,
    std::string file_name,
    int error_severity,
    std::optional<bool> binary) {
  auto _file = file.c_str(); // ptr, inout, required
  int _iunit{};
  auto _file_name = file_name.c_str();
  bool binary_lvalue;
  auto* _binary{&binary_lvalue};
  if (binary.has_value()) {
    binary_lvalue = binary.value();
  } else {
    _binary = nullptr;
  }
  fortran_tao_open_file(
      /* c_Char */ _file,
      /* c_Int& */ _iunit,
      /* c_Char */ _file_name,
      /* c_Int& */ error_severity,
      /* c_Bool* */ _binary);
  return _iunit;
}
void Tao::tao_json_cmd(std::string input_str) {
  auto _input_str = input_str.c_str(); // ptr, inout, required
  fortran_tao_json_cmd(/* c_Char */ _input_str);
}
void Tao::tao_datum_name(
    TaoDataProxy& datum,
    std::optional<bool> show_universe,
    std::string datum_name) {
  bool show_universe_lvalue;
  auto* _show_universe{&show_universe_lvalue};
  if (show_universe.has_value()) {
    show_universe_lvalue = show_universe.value();
  } else {
    _show_universe = nullptr;
  }
  auto _datum_name = datum_name.c_str(); // ptr, inout, required
  fortran_tao_datum_name(
      /* void* */ datum.get_fortran_ptr(),
      /* c_Bool* */ _show_universe,
      /* c_Char */ _datum_name);
}
void Tao::tao_rad_int_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_rad_int) {
  auto _data_type = data_type.c_str(); // ptr, inout, required
  auto _data_source = data_source.c_str(); // ptr, inout, required
  fortran_tao_rad_int_calc_needed(
      /* c_Char */ _data_type,
      /* c_Char */ _data_source,
      /* c_Bool& */ do_rad_int);
}
void Tao::tao_var_target_calc() {
  fortran_tao_var_target_calc();
}
void Tao::tao_cmd_history_record(std::string cmd) {
  auto _cmd = cmd.c_str(); // ptr, inout, required
  fortran_tao_cmd_history_record(/* c_Char */ _cmd);
}
void Tao::tao_re_execute(std::string string, bool err) {
  auto _string = string.c_str(); // ptr, inout, required
  fortran_tao_re_execute(/* c_Char */ _string, /* c_Bool& */ err);
}
std::string Tao::tao_next_word(std::string line) {
  auto _line = line.c_str(); // ptr, inout, required
  char _word[4096];
  fortran_tao_next_word(/* c_Char */ _line, /* c_Char */ _word);
  return _word;
}
void Tao::tao_spin_matrices_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_calc) {
  auto _data_type = data_type.c_str(); // ptr, inout, required
  auto _data_source = data_source.c_str(); // ptr, inout, required
  fortran_tao_spin_matrices_calc_needed(
      /* c_Char */ _data_type,
      /* c_Char */ _data_source,
      /* c_Bool& */ do_calc);
}
void Tao::tao_init_lattice(std::string lat_file, bool err_flag) {
  auto _lat_file = lat_file.c_str(); // ptr, inout, required
  fortran_tao_init_lattice(/* c_Char */ _lat_file, /* c_Bool& */ err_flag);
}
void Tao::tao_pause_cmd(double time) {
  fortran_tao_pause_cmd(/* c_Real& */ time);
}
void Tao::tao_draw_plots(std::optional<bool> do_clear) {
  bool do_clear_lvalue;
  auto* _do_clear{&do_clear_lvalue};
  if (do_clear.has_value()) {
    do_clear_lvalue = do_clear.value();
  } else {
    _do_clear = nullptr;
  }
  fortran_tao_draw_plots(/* c_Bool* */ _do_clear);
}
bool Tao::tao_open_scratch_file(int iu) {
  bool _err{};
  fortran_tao_open_scratch_file(/* c_Bool& */ _err, /* c_Int& */ iu);
  return _err;
}
void Tao::tao_create_plot_window() {
  fortran_tao_create_plot_window();
}
void Tao::tao_destroy_plot_window() {
  fortran_tao_destroy_plot_window();
}
void Tao::tao_var_repoint() {
  fortran_tao_var_repoint();
}
void Tao::tao_quiet_set(std::string set) {
  auto _set = set.c_str();
  fortran_tao_quiet_set(/* c_Char */ _set);
}
void Tao::tao_mark_lattice_ele(LatProxy& lat) {
  fortran_tao_mark_lattice_ele(/* void* */ lat.get_fortran_ptr());
}
void Tao::tao_use_data(std::string action, std::string data_name) {
  auto _action = action.c_str();
  auto _data_name = data_name.c_str();
  fortran_tao_use_data(/* c_Char */ _action, /* c_Char */ _data_name);
}
void Tao::tao_set_opt_vars(
    RealAlloc1D& var_vec,
    std::optional<bool> print_limit_warning) {
  // intent=in allocatable general array
  bool print_limit_warning_lvalue;
  auto* _print_limit_warning{&print_limit_warning_lvalue};
  if (print_limit_warning.has_value()) {
    print_limit_warning_lvalue = print_limit_warning.value();
  } else {
    _print_limit_warning = nullptr;
  }
  fortran_tao_set_opt_vars(
      /* void* */ var_vec.get_fortran_ptr(),
      /* c_Bool* */ _print_limit_warning);
}
void Tao::tao_close_command_file() {
  fortran_tao_close_command_file();
}
Tao::TaoLatticeCalc Tao::tao_lattice_calc() {
  bool _calc_ok{};
  bool _print_err{};
  fortran_tao_lattice_calc(/* c_Bool& */ _calc_ok, /* c_Bool& */ _print_err);
  return TaoLatticeCalc{_calc_ok, _print_err};
}
std::string Tao::tao_set_invalid(
    TaoDataProxy& datum,
    std::string message,
    std::optional<bool> exterminate,
    std::optional<int> err_level,
    std::optional<bool> print_err) {
  auto _message = message.c_str();
  char _why_invalid[4096];
  bool exterminate_lvalue;
  auto* _exterminate{&exterminate_lvalue};
  if (exterminate.has_value()) {
    exterminate_lvalue = exterminate.value();
  } else {
    _exterminate = nullptr;
  }
  int err_level_lvalue;
  auto* _err_level{&err_level_lvalue};
  if (err_level.has_value()) {
    err_level_lvalue = err_level.value();
  } else {
    _err_level = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_set_invalid(
      /* void* */ datum.get_fortran_ptr(),
      /* c_Char */ _message,
      /* c_Char */ _why_invalid,
      /* c_Bool* */ _exterminate,
      /* c_Int* */ _err_level,
      /* c_Bool* */ _print_err);
  return _why_invalid;
}
void Tao::tao_lat_emit_calc(
    int plane,
    int emit_type,
    EleProxy& ele,
    NormalModesProxy& modes,
    double emit) {
  fortran_tao_lat_emit_calc(
      /* c_Int& */ plane,
      /* c_Int& */ emit_type,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ modes.get_fortran_ptr(),
      /* c_Real& */ emit);
}
std::string Tao::tao_is_valid_name(std::string name, bool is_valid) {
  auto _name = name.c_str();
  char _why_invalid[4096];
  fortran_tao_is_valid_name(
      /* c_Char */ _name, /* c_Char */ _why_invalid, /* c_Bool& */ is_valid);
  return _why_invalid;
}
void Tao::tao_python_cmd(std::string input_str) {
  auto _input_str = input_str.c_str();
  fortran_tao_python_cmd(/* c_Char */ _input_str);
}
void Tao::tao_scale_ping_data(TaoUniverseProxy& u) {
  fortran_tao_scale_ping_data(/* void* */ u.get_fortran_ptr());
}
void Tao::tao_plot_cmd(std::string where, std::string component) {
  auto _where = where.c_str();
  auto _component = component.c_str();
  fortran_tao_plot_cmd(/* c_Char */ _where, /* c_Char */ _component);
}
Tao::TaoEvaluateADatum Tao::tao_evaluate_a_datum(
    TaoDataProxy& datum,
    TaoUniverseProxy& u,
    TaoLatticeProxy& tao_lat,
    std::optional<bool> called_from_lat_calc,
    std::optional<bool> print_err) {
  double _datum_value{};
  bool _valid_value{};
  char _why_invalid[4096];
  bool called_from_lat_calc_lvalue;
  auto* _called_from_lat_calc{&called_from_lat_calc_lvalue};
  if (called_from_lat_calc.has_value()) {
    called_from_lat_calc_lvalue = called_from_lat_calc.value();
  } else {
    _called_from_lat_calc = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_evaluate_a_datum(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ u.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* c_Real& */ _datum_value,
      /* c_Bool& */ _valid_value,
      /* c_Char */ _why_invalid,
      /* c_Bool* */ _called_from_lat_calc,
      /* c_Bool* */ _print_err);
  return TaoEvaluateADatum{_datum_value, _valid_value, _why_invalid};
}
void Tao::tao_srdt_calc_needed(
    std::string data_type,
    std::string data_source,
    int do_srdt) {
  auto _data_type = data_type.c_str(); // ptr, inout, required
  auto _data_source = data_source.c_str(); // ptr, inout, required
  fortran_tao_srdt_calc_needed(
      /* c_Char */ _data_type, /* c_Char */ _data_source, /* c_Int& */ do_srdt);
}
void Tao::tao_regression_test() {
  fortran_tao_regression_test();
}
void Tao::tao_pipe_cmd(std::string input_str) {
  auto _input_str = input_str.c_str();
  fortran_tao_pipe_cmd(/* c_Char */ _input_str);
}
void Tao::tao_optimization_status(TaoDataProxy& datum, std::string why_str) {
  auto _why_str = why_str.c_str(); // ptr, inout, required
  fortran_tao_optimization_status(
      /* void* */ datum.get_fortran_ptr(), /* c_Char */ _why_str);
}
Tao::TaoEvaluateLatOrBeamData Tao::tao_evaluate_lat_or_beam_data(
    std::string data_name,
    bool print_err,
    std::string default_source,
    optional_ref<EleProxy> dflt_ele_ref,
    optional_ref<EleProxy> dflt_ele_start,
    optional_ref<EleProxy> dflt_ele,
    std::optional<std::string> dflt_component,
    std::optional<int> dflt_uni,
    std::optional<int> dflt_eval_point,
    std::optional<double> dflt_s_offset) {
  bool _err{};
  auto _data_name = data_name.c_str();
  // intent=out allocatable general array
  auto values{RealAlloc1D()};
  auto _default_source = default_source.c_str(); // ptr, inout, required
  auto* _dflt_ele_ref = dflt_ele_ref.has_value()
      ? dflt_ele_ref->get().get_fortran_ptr()
      : nullptr; // input, optional
  auto* _dflt_ele_start = dflt_ele_start.has_value()
      ? dflt_ele_start->get().get_fortran_ptr()
      : nullptr; // input, optional
  auto* _dflt_ele = dflt_ele.has_value() ? dflt_ele->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  const char* _dflt_component =
      dflt_component.has_value() ? dflt_component->c_str() : nullptr;
  int dflt_uni_lvalue;
  auto* _dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  int dflt_eval_point_lvalue;
  auto* _dflt_eval_point{&dflt_eval_point_lvalue};
  if (dflt_eval_point.has_value()) {
    dflt_eval_point_lvalue = dflt_eval_point.value();
  } else {
    _dflt_eval_point = nullptr;
  }
  double dflt_s_offset_lvalue;
  auto* _dflt_s_offset{&dflt_s_offset_lvalue};
  if (dflt_s_offset.has_value()) {
    dflt_s_offset_lvalue = dflt_s_offset.value();
  } else {
    _dflt_s_offset = nullptr;
  }
  fortran_tao_evaluate_lat_or_beam_data(
      /* c_Bool& */ _err,
      /* c_Char */ _data_name,
      /* void* */ values.get_fortran_ptr(),
      /* c_Bool& */ print_err,
      /* c_Char */ _default_source,
      /* void* */ _dflt_ele_ref,
      /* void* */ _dflt_ele_start,
      /* void* */ _dflt_ele,
      /* c_Char */ _dflt_component,
      /* c_Int* */ _dflt_uni,
      /* c_Int* */ _dflt_eval_point,
      /* c_Real* */ _dflt_s_offset);
  return TaoEvaluateLatOrBeamData{_err, std::move(values)};
}
Tao::TaoToPhaseAndCouplingReading Tao::tao_to_phase_and_coupling_reading(
    EleProxy& ele,
    std::string why_invalid,
    TaoDataProxy& datum) {
  BpmPhaseCouplingProxy _bpm_data;
  bool _valid_value{};
  auto _why_invalid = why_invalid.c_str(); // ptr, inout, required
  fortran_tao_to_phase_and_coupling_reading(
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ _bpm_data.get_fortran_ptr(),
      /* c_Bool& */ _valid_value,
      /* c_Char */ _why_invalid,
      /* void* */ datum.get_fortran_ptr());
  return TaoToPhaseAndCouplingReading{std::move(_bpm_data), _valid_value};
}
Tao::TaoGetData Tao::tao_get_data() {
  // intent=out allocatable general array
  auto data_value{RealAlloc1D()};
  // intent=out allocatable general array
  auto data_weight{RealAlloc1D()};
  // intent=out allocatable general array
  auto data_meas_value{RealAlloc1D()};
  // intent=out allocatable general array
  auto data_ix_dModel{IntAlloc1D()};
  fortran_tao_get_data(
      /* void* */ data_value.get_fortran_ptr(),
      /* void* */ data_weight.get_fortran_ptr(),
      /* void* */ data_meas_value.get_fortran_ptr(),
      /* void* */ data_ix_dModel.get_fortran_ptr());
  return TaoGetData{
      std::move(data_value),
      std::move(data_weight),
      std::move(data_meas_value),
      std::move(data_ix_dModel)};
}
std::string Tao::tao_expression_hash_substitute(
    std::string expression_in,
    optional_ref<EleProxy> eval_ele) {
  auto _expression_in = expression_in.c_str();
  char _expression_out[4096];
  auto* _eval_ele = eval_ele.has_value() ? eval_ele->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  fortran_tao_expression_hash_substitute(
      /* c_Char */ _expression_in,
      /* c_Char */ _expression_out,
      /* void* */ _eval_ele);
  return _expression_out;
}
void Tao::tao_load_this_datum(
    RealAlloc1D& vec,
    EleProxy& ele_ref,
    EleProxy& ele_start,
    EleProxy& ele,
    double datum_value,
    bool valid_value,
    TaoDataProxy& datum,
    BranchProxy& branch,
    optional_ref<std::string> why_invalid,
    optional_ref<BoolAlloc1D> good) {
  // intent=inout allocatable general array
  auto _ele_ref = &ele_ref; // input, required, pointer
  auto _ele_start = &ele_start; // input, required, pointer
  auto _ele = &ele; // input, required, pointer
  const char* _why_invalid =
      why_invalid.has_value() ? why_invalid->get().c_str() : nullptr;
  // intent=inout allocatable general array
  auto* _good = good.has_value() ? good->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  fortran_tao_load_this_datum(
      /* void* */ vec.get_fortran_ptr(),
      /* void* */ &ele_ref,
      /* void* */ &ele_start,
      /* void* */ &ele,
      /* c_Real& */ datum_value,
      /* c_Bool& */ valid_value,
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* c_Char */ _why_invalid,
      /* void* */ _good);
}
double Tao::tao_datum_s_position(TaoDataProxy& datum, EleProxy& ele) {
  double _s_pos{};
  fortran_tao_datum_s_position(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* c_Real& */ _s_pos);
  return _s_pos;
}
Tao::TaoDatumIntegrate Tao::tao_datum_integrate(
    TaoDataProxy& datum,
    BranchProxy& branch,
    RealAlloc1D& s_pos,
    RealAlloc1D& values) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  bool _valid_value{};
  char _why_invalid[4096];
  double _result{};
  fortran_tao_datum_integrate(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ s_pos.get_fortran_ptr(),
      /* void* */ values.get_fortran_ptr(),
      /* c_Bool& */ _valid_value,
      /* c_Char */ _why_invalid,
      /* c_Real& */ _result);
  return TaoDatumIntegrate{_valid_value, _why_invalid, _result};
}
Tao::TaoTrackingEleIndex Tao::tao_tracking_ele_index(
    EleProxy& ele,
    TaoDataProxy& datum) {
  auto _ele = &ele; // input, required, pointer
  int _ix_branch{};
  int _ix_ele{};
  fortran_tao_tracking_ele_index(
      /* void* */ &ele,
      /* void* */ datum.get_fortran_ptr(),
      /* c_Int& */ _ix_branch,
      /* c_Int& */ _ix_ele);
  return TaoTrackingEleIndex{_ix_branch, _ix_ele};
}
void Tao::integrate_min(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum) {
  // intent=inout allocatable general array
  fortran_integrate_min(
      /* c_Int& */ ix_start,
      /* c_Int& */ ix_ele,
      /* c_Real& */ datum_value,
      /* c_Int& */ ix_m,
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ vec.get_fortran_ptr(),
      /* void* */ datum.get_fortran_ptr());
}
void Tao::integrate_max(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum) {
  // intent=inout allocatable general array
  fortran_integrate_max(
      /* c_Int& */ ix_start,
      /* c_Int& */ ix_ele,
      /* c_Real& */ datum_value,
      /* c_Int& */ ix_m,
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ vec.get_fortran_ptr(),
      /* void* */ datum.get_fortran_ptr());
}
void Tao::tao_scratch_values_calc(
    EleProxy& ele_ref,
    EleProxy& ele_start,
    EleProxy& ele,
    TaoDataProxy& datum,
    BranchProxy& branch,
    CoordProxyAlloc1D& orbit) {
  auto _ele_ref = &ele_ref; // input, required, pointer
  auto _ele_start = &ele_start; // input, required, pointer
  auto _ele = &ele; // input, required, pointer
  // intent=inout allocatable type array
  fortran_tao_scratch_values_calc(
      /* void* */ &ele_ref,
      /* void* */ &ele_start,
      /* void* */ &ele,
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ branch.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr());
}
double Tao::tao_do_wire_scan(EleProxy& ele, double theta, BeamProxy& beam) {
  double _moment{};
  fortran_tao_do_wire_scan(
      /* void* */ ele.get_fortran_ptr(),
      /* c_Real& */ theta,
      /* void* */ beam.get_fortran_ptr(),
      /* c_Real& */ _moment);
  return _moment;
}
Tao::TaoPointerToDatumEle Tao::tao_pointer_to_datum_ele(
    LatProxy& lat,
    std::string ele_name,
    int ix_ele,
    TaoDataProxy& datum,
    std::optional<bool> print_err) {
  auto _ele_name = ele_name.c_str(); // ptr, inout, required
  bool _valid{};
  char _why_invalid[4096];
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  EleProxy _ele;
  fortran_tao_pointer_to_datum_ele(
      /* void* */ lat.get_fortran_ptr(),
      /* c_Char */ _ele_name,
      /* c_Int& */ ix_ele,
      /* void* */ datum.get_fortran_ptr(),
      /* c_Bool& */ _valid,
      /* c_Char */ _why_invalid,
      /* c_Bool* */ _print_err,
      /* void* */ _ele.get_fortran_ptr());
  return TaoPointerToDatumEle{_valid, _why_invalid, std::move(_ele)};
}
Tao::TaoEvaluateDatumAtS Tao::tao_evaluate_datum_at_s(
    TaoDataProxy& datum,
    TaoLatticeProxy& tao_lat,
    EleProxy& ele,
    EleProxy& ele_ref,
    bool valid_value) {
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
      /* c_Bool& */ valid_value,
      /* c_Char */ _err_str,
      /* c_Bool& */ _bad_datum,
      /* c_Real& */ _value);
  return TaoEvaluateDatumAtS{_err_str, _bad_datum, _value};
}
void Tao::tao_to_int(std::string str, int i_int, bool err) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_tao_to_int(/* c_Char */ _str, /* c_Int& */ i_int, /* c_Bool& */ err);
}
Tao::TaoEleGeometryWithMisalignments Tao::tao_ele_geometry_with_misalignments(
    TaoDataProxy& datum,
    EleProxy& ele) {
  bool _valid_value{};
  char _why_invalid[4096];
  double _value{};
  fortran_tao_ele_geometry_with_misalignments(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* c_Bool& */ _valid_value,
      /* c_Char */ _why_invalid,
      /* c_Real& */ _value);
  return TaoEleGeometryWithMisalignments{_valid_value, _why_invalid, _value};
}
Tao::TaoEvalFloorOrbit Tao::tao_eval_floor_orbit(
    TaoDataProxy& datum,
    EleProxy& ele,
    CoordProxy& orbit,
    BunchParamsProxy& bunch_params) {
  bool _valid_value{};
  char _why_invalid[4096];
  double _value{};
  fortran_tao_eval_floor_orbit(
      /* void* */ datum.get_fortran_ptr(),
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* c_Bool& */ _valid_value,
      /* c_Char */ _why_invalid,
      /* c_Real& */ _value);
  return TaoEvalFloorOrbit{_valid_value, _why_invalid, _value};
}
void Tao::tao_use_var(std::string action, std::string var_name) {
  auto _action = action.c_str();
  auto _var_name = var_name.c_str();
  fortran_tao_use_var(/* c_Char */ _action, /* c_Char */ _var_name);
}
Tao::TaoParseElementParamStr Tao::tao_parse_element_param_str(
    std::string in_str) {
  bool _err{};
  auto _in_str = in_str.c_str();
  char _uni[4096];
  char _element[4096];
  char _parameter[4096];
  int _where{};
  char _component[4096];
  fortran_tao_parse_element_param_str(
      /* c_Bool& */ _err,
      /* c_Char */ _in_str,
      /* c_Char */ _uni,
      /* c_Char */ _element,
      /* c_Char */ _parameter,
      /* c_Int& */ _where,
      /* c_Char */ _component);
  return TaoParseElementParamStr{
      _err, _uni, _element, _parameter, _where, _component};
}
void Tao::tao_data_coupling_init(BranchProxy& branch) {
  fortran_tao_data_coupling_init(/* void* */ branch.get_fortran_ptr());
}
void Tao::tao_write_cmd(std::string what) {
  auto _what = what.c_str();
  fortran_tao_write_cmd(/* c_Char */ _what);
}
void Tao::tao_x_axis_cmd(std::string where, std::string what) {
  auto _where = where.c_str();
  auto _what = what.c_str();
  fortran_tao_x_axis_cmd(/* c_Char */ _where, /* c_Char */ _what);
}
bool Tao::tao_merit(double this_merit) {
  bool _calc_ok{};
  fortran_tao_merit(/* c_Bool& */ _calc_ok, /* c_Real& */ this_merit);
  return _calc_ok;
}
void Tao::tao_d2_d1_name(
    TaoD1DataProxy& d1,
    std::optional<bool> show_universe,
    std::string d2_d1_name) {
  bool show_universe_lvalue;
  auto* _show_universe{&show_universe_lvalue};
  if (show_universe.has_value()) {
    show_universe_lvalue = show_universe.value();
  } else {
    _show_universe = nullptr;
  }
  auto _d2_d1_name = d2_d1_name.c_str(); // ptr, inout, required
  fortran_tao_d2_d1_name(
      /* void* */ d1.get_fortran_ptr(),
      /* c_Bool* */ _show_universe,
      /* c_Char */ _d2_d1_name);
}
Tao::TaoParamValueAtS Tao::tao_param_value_at_s(
    std::string dat_name,
    EleProxy& ele_to_s,
    EleProxy& ele_here,
    CoordProxy& orbit,
    double value) {
  auto _dat_name = dat_name.c_str(); // ptr, inout, required
  bool _err_flag{};
  char _why_invalid[4096];
  bool _print_err{};
  bool _bad_datum{};
  fortran_tao_param_value_at_s(
      /* c_Char */ _dat_name,
      /* void* */ ele_to_s.get_fortran_ptr(),
      /* void* */ ele_here.get_fortran_ptr(),
      /* void* */ orbit.get_fortran_ptr(),
      /* c_Bool& */ _err_flag,
      /* c_Char */ _why_invalid,
      /* c_Bool& */ _print_err,
      /* c_Bool& */ _bad_datum,
      /* c_Real& */ value);
  return TaoParamValueAtS{_err_flag, _why_invalid, _print_err, _bad_datum};
}
void Tao::tao_show_cmd(std::string what) {
  auto _what = what.c_str();
  fortran_tao_show_cmd(/* c_Char */ _what);
}
void Tao::tao_evaluate_tune(
    std::string q_str,
    double q0,
    bool delta_input,
    double q_val) {
  auto _q_str = q_str.c_str();
  fortran_tao_evaluate_tune(
      /* c_Char */ _q_str,
      /* c_Real& */ q0,
      /* c_Bool& */ delta_input,
      /* c_Real& */ q_val);
}
void Tao::re_allocate_c_double(
    RealAlloc1D& re,
    int n,
    std::optional<bool> exact,
    optional_ref<double> init_val) {
  // intent=inout allocatable general array
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  auto* _init_val =
      init_val.has_value() ? &init_val->get() : nullptr; // inout, optional
  fortran_re_allocate_c_double(
      /* void* */ re.get_fortran_ptr(),
      /* c_Int& */ n,
      /* c_Bool* */ _exact,
      /* c_Real* */ _init_val);
}
void Tao::tao_c_out_io_buffer_reset() {
  fortran_tao_c_out_io_buffer_reset();
}
std::string Tao::tao_subin_uni_number(
    std::string name_in,
    int ix_uni,
    bool ok) {
  auto _name_in = name_in.c_str();
  char _name_out[4096];
  fortran_tao_subin_uni_number(
      /* c_Char */ _name_in,
      /* c_Int& */ ix_uni,
      /* c_Char */ _name_out,
      /* c_Bool& */ ok);
  return _name_out;
}
bool Tao::tao_parse_command_args(optional_ref<std::string> cmd_line) {
  bool _error{};
  const char* _cmd_line =
      cmd_line.has_value() ? cmd_line->get().c_str() : nullptr;
  fortran_tao_parse_command_args(/* c_Bool& */ _error, /* c_Char */ _cmd_line);
  return _error;
}
void Tao::tao_chrom_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_chrom) {
  auto _data_type = data_type.c_str(); // ptr, inout, required
  auto _data_source = data_source.c_str(); // ptr, inout, required
  fortran_tao_chrom_calc_needed(
      /* c_Char */ _data_type,
      /* c_Char */ _data_source,
      /* c_Bool& */ do_chrom);
}
bool Tao::tao_lm_optimizer() {
  bool _abort{};
  fortran_tao_lm_optimizer(/* c_Bool& */ _abort);
  return _abort;
}
std::string Tao::tao_get_user_input(
    std::optional<std::string> prompt_str,
    std::optional<bool> wait_flag,
    std::optional<std::string> cmd_in) {
  char _cmd_out[4096];
  const char* _prompt_str =
      prompt_str.has_value() ? prompt_str->c_str() : nullptr;
  bool wait_flag_lvalue;
  auto* _wait_flag{&wait_flag_lvalue};
  if (wait_flag.has_value()) {
    wait_flag_lvalue = wait_flag.value();
  } else {
    _wait_flag = nullptr;
  }
  const char* _cmd_in = cmd_in.has_value() ? cmd_in->c_str() : nullptr;
  fortran_tao_get_user_input(
      /* c_Char */ _cmd_out,
      /* c_Char */ _prompt_str,
      /* c_Bool* */ _wait_flag,
      /* c_Char */ _cmd_in);
  return _cmd_out;
}
void Tao::tao_top10_merit_categories_print(int iunit) {
  fortran_tao_top10_merit_categories_print(/* c_Int& */ iunit);
}
void Tao::tao_top10_derivative_print() {
  fortran_tao_top10_derivative_print();
}
void Tao::tao_show_constraints(int iunit, std::string form) {
  auto _form = form.c_str();
  fortran_tao_show_constraints(/* c_Int& */ iunit, /* c_Char */ _form);
}
void Tao::tao_var_write(
    std::string out_file,
    std::optional<bool> show_good_opt_only,
    std::optional<bool> tao_format) {
  auto _out_file = out_file.c_str();
  bool show_good_opt_only_lvalue;
  auto* _show_good_opt_only{&show_good_opt_only_lvalue};
  if (show_good_opt_only.has_value()) {
    show_good_opt_only_lvalue = show_good_opt_only.value();
  } else {
    _show_good_opt_only = nullptr;
  }
  bool tao_format_lvalue;
  auto* _tao_format{&tao_format_lvalue};
  if (tao_format.has_value()) {
    tao_format_lvalue = tao_format.value();
  } else {
    _tao_format = nullptr;
  }
  fortran_tao_var_write(
      /* c_Char */ _out_file,
      /* c_Bool* */ _show_good_opt_only,
      /* c_Bool* */ _tao_format);
}
void Tao::tao_spin_tracking_turn_on() {
  fortran_tao_spin_tracking_turn_on();
}
void Tao::tao_init_data(std::string data_file) {
  auto _data_file = data_file.c_str();
  fortran_tao_init_data(/* c_Char */ _data_file);
}
void Tao::tao_init_data_end_stuff() {
  fortran_tao_init_data_end_stuff();
}
void Tao::tao_allocate_data_array(
    TaoUniverseProxy& u,
    int n_data,
    optional_ref<bool> exact) {
  auto* _exact = exact.has_value() ? &exact->get() : nullptr; // inout, optional
  fortran_tao_allocate_data_array(
      /* void* */ u.get_fortran_ptr(),
      /* c_Int& */ n_data,
      /* c_Bool* */ _exact);
}
void Tao::tao_d2_data_stuffit(
    TaoUniverseProxy& u,
    std::string d2_name,
    int n_d1_data) {
  auto _d2_name = d2_name.c_str(); // ptr, inout, required
  fortran_tao_d2_data_stuffit(
      /* void* */ u.get_fortran_ptr(),
      /* c_Char */ _d2_name,
      /* c_Int& */ n_d1_data);
}
void Tao::tao_init_data_in_universe(
    TaoUniverseProxy& u,
    int n_d2_add,
    optional_ref<bool> keep_existing_data) {
  auto* _keep_existing_data = keep_existing_data.has_value()
      ? &keep_existing_data->get()
      : nullptr; // inout, optional
  fortran_tao_init_data_in_universe(
      /* void* */ u.get_fortran_ptr(),
      /* c_Int& */ n_d2_add,
      /* c_Bool* */ _keep_existing_data);
}
ResonanceHProxyAlloc1D Tao::tao_add_to_normal_mode_h_array(std::string h_str) {
  auto _h_str = h_str.c_str();
  // intent=out allocatable type array
  auto h_array{ResonanceHProxyAlloc1D()};
  fortran_tao_add_to_normal_mode_h_array(
      /* c_Char */ _h_str, /* void* */ h_array.get_fortran_ptr());
  return std::move(h_array);
}
void Tao::tao_plot_setup() {
  fortran_tao_plot_setup();
}
void Tao::tao_set_var_useit_opt() {
  fortran_tao_set_var_useit_opt();
}
int Tao::tao_phase_space_axis_index(std::string data_type, bool err) {
  auto _data_type = data_type.c_str();
  int _ix_axis{};
  fortran_tao_phase_space_axis_index(
      /* c_Char */ _data_type, /* c_Bool& */ err, /* c_Int& */ _ix_axis);
  return _ix_axis;
}
Tao::TaoParticleDataValue Tao::tao_particle_data_value(
    std::string data_type,
    CoordProxyAlloc1D& p,
    EleProxy& ele,
    int ix_bunch) {
  auto _data_type = data_type.c_str();
  // intent=in allocatable type array
  // intent=out allocatable general array
  auto value{RealAlloc1D()};
  bool _err{};
  fortran_tao_particle_data_value(
      /* c_Char */ _data_type,
      /* void* */ p.get_fortran_ptr(),
      /* void* */ value.get_fortran_ptr(),
      /* c_Bool& */ _err,
      /* void* */ ele.get_fortran_ptr(),
      /* c_Int& */ ix_bunch);
  return TaoParticleDataValue{std::move(value), _err};
}
void Tao::tao_constraint_type_name(
    TaoDataProxy& datum,
    std::string datum_name) {
  auto _datum_name = datum_name.c_str(); // ptr, inout, required
  fortran_tao_constraint_type_name(
      /* void* */ datum.get_fortran_ptr(), /* c_Char */ _datum_name);
}
void Tao::tao_read_cmd(
    std::string which,
    std::string unis,
    std::string file,
    bool silent) {
  auto _which = which.c_str(); // ptr, inout, required
  auto _unis = unis.c_str();
  auto _file = file.c_str(); // ptr, inout, required
  fortran_tao_read_cmd(
      /* c_Char */ _which,
      /* c_Char */ _unis,
      /* c_Char */ _file,
      /* c_Bool& */ silent);
}
Tao::TaoGetOptVars Tao::tao_get_opt_vars() {
  // intent=out allocatable general array
  auto var_value{RealAlloc1D()};
  // intent=out allocatable general array
  auto var_step{RealAlloc1D()};
  // intent=out allocatable general array
  auto var_delta{RealAlloc1D()};
  // intent=out allocatable general array
  auto var_weight{RealAlloc1D()};
  // intent=out allocatable general array
  auto var_ix{IntAlloc1D()};
  bool _ignore_if_weight_is_zero{};
  bool _ignore_if_not_limited{};
  fortran_tao_get_opt_vars(
      /* void* */ var_value.get_fortran_ptr(),
      /* void* */ var_step.get_fortran_ptr(),
      /* void* */ var_delta.get_fortran_ptr(),
      /* void* */ var_weight.get_fortran_ptr(),
      /* void* */ var_ix.get_fortran_ptr(),
      /* c_Bool& */ _ignore_if_weight_is_zero,
      /* c_Bool& */ _ignore_if_not_limited);
  return TaoGetOptVars{
      std::move(var_value),
      std::move(var_step),
      std::move(var_delta),
      std::move(var_weight),
      std::move(var_ix),
      _ignore_if_weight_is_zero,
      _ignore_if_not_limited};
}
void Tao::tao_init_variables(std::string var_file) {
  auto _var_file = var_file.c_str();
  fortran_tao_init_variables(/* c_Char */ _var_file);
}
void Tao::tao_allocate_v1_var(int n_v1, bool save_old) {
  fortran_tao_allocate_v1_var(/* c_Int& */ n_v1, /* c_Bool& */ save_old);
}
void Tao::tao_allocate_var_array(int n_var, bool default_good_user) {
  fortran_tao_allocate_var_array(
      /* c_Int& */ n_var, /* c_Bool& */ default_good_user);
}
void Tao::tao_init_plotting(std::string plot_file) {
  auto _plot_file = plot_file.c_str(); // ptr, inout, required
  fortran_tao_init_plotting(/* c_Char */ _plot_file);
}
bool Tao::tao_x_scale_cmd(
    std::string where,
    double x_min_in,
    double x_max_in,
    std::optional<bool> include_wall,
    std::optional<std::string> gang,
    std::optional<bool> exact,
    std::optional<bool> turn_autoscale_off) {
  auto _where = where.c_str();
  bool _err{};
  bool include_wall_lvalue;
  auto* _include_wall{&include_wall_lvalue};
  if (include_wall.has_value()) {
    include_wall_lvalue = include_wall.value();
  } else {
    _include_wall = nullptr;
  }
  const char* _gang = gang.has_value() ? gang->c_str() : nullptr;
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  bool turn_autoscale_off_lvalue;
  auto* _turn_autoscale_off{&turn_autoscale_off_lvalue};
  if (turn_autoscale_off.has_value()) {
    turn_autoscale_off_lvalue = turn_autoscale_off.value();
  } else {
    _turn_autoscale_off = nullptr;
  }
  fortran_tao_x_scale_cmd(
      /* c_Char */ _where,
      /* c_Real& */ x_min_in,
      /* c_Real& */ x_max_in,
      /* c_Bool& */ _err,
      /* c_Bool* */ _include_wall,
      /* c_Char */ _gang,
      /* c_Bool* */ _exact,
      /* c_Bool* */ _turn_autoscale_off);
  return _err;
}
Tao::TaoPickUniverse Tao::tao_pick_universe(
    std::string name_in,
    std::optional<int> dflt_uni,
    std::optional<bool> pure_uni) {
  auto _name_in = name_in.c_str();
  char _name_out[4096];
  // intent=out allocatable general array
  auto picked{BoolAlloc1D()};
  bool _err{};
  int _ix_uni{};
  bool _explicit_uni{};
  int dflt_uni_lvalue;
  auto* _dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  bool pure_uni_lvalue;
  auto* _pure_uni{&pure_uni_lvalue};
  if (pure_uni.has_value()) {
    pure_uni_lvalue = pure_uni.value();
  } else {
    _pure_uni = nullptr;
  }
  fortran_tao_pick_universe(
      /* c_Char */ _name_in,
      /* c_Char */ _name_out,
      /* void* */ picked.get_fortran_ptr(),
      /* c_Bool& */ _err,
      /* c_Int& */ _ix_uni,
      /* c_Bool& */ _explicit_uni,
      /* c_Int* */ _dflt_uni,
      /* c_Bool* */ _pure_uni);
  return TaoPickUniverse{
      _name_out, std::move(picked), _err, _ix_uni, _explicit_uni};
}
void Tao::tao_wave_cmd(
    std::string curve_name,
    std::string plot_place,
    bool err_flag) {
  auto _curve_name = curve_name.c_str();
  auto _plot_place = plot_place.c_str();
  fortran_tao_wave_cmd(
      /* c_Char */ _curve_name,
      /* c_Char */ _plot_place,
      /* c_Bool& */ err_flag);
}
void Tao::tao_clear_cmd(std::string cmd_line) {
  auto _cmd_line = cmd_line.c_str();
  fortran_tao_clear_cmd(/* c_Char */ _cmd_line);
}
void Tao::tao_symbol_import_from_lat(LatProxy& lat) {
  fortran_tao_symbol_import_from_lat(/* void* */ lat.get_fortran_ptr());
}
void Tao::tao_datum_has_associated_ele(
    std::string data_type,
    std::optional<int> branch_geometry,
    int has_associated_ele) {
  auto _data_type = data_type.c_str();
  int branch_geometry_lvalue;
  auto* _branch_geometry{&branch_geometry_lvalue};
  if (branch_geometry.has_value()) {
    branch_geometry_lvalue = branch_geometry.value();
  } else {
    _branch_geometry = nullptr;
  }
  fortran_tao_datum_has_associated_ele(
      /* c_Char */ _data_type,
      /* c_Int* */ _branch_geometry,
      /* c_Int& */ has_associated_ele);
}
int Tao::tao_count_strings(std::string string, std::string pattern) {
  auto _string = string.c_str();
  auto _pattern = pattern.c_str();
  int _num{};
  fortran_tao_count_strings(
      /* c_Char */ _string, /* c_Char */ _pattern, /* c_Int& */ _num);
  return _num;
}
void Tao::tao_ptc_normal_form(
    bool do_calc,
    TaoLatticeProxy& tao_lat,
    int ix_branch,
    std::optional<int> rf_on) {
  int rf_on_lvalue;
  auto* _rf_on{&rf_on_lvalue};
  if (rf_on.has_value()) {
    rf_on_lvalue = rf_on.value();
  } else {
    _rf_on = nullptr;
  }
  fortran_tao_ptc_normal_form(
      /* c_Bool& */ do_calc,
      /* void* */ tao_lat.get_fortran_ptr(),
      /* c_Int& */ ix_branch,
      /* c_Int* */ _rf_on);
}
void Tao::tao_one_turn_map_calc_needed(
    std::string data_type,
    std::string data_source,
    bool do_one_turn_map) {
  auto _data_type = data_type.c_str(); // ptr, inout, required
  auto _data_source = data_source.c_str(); // ptr, inout, required
  fortran_tao_one_turn_map_calc_needed(
      /* c_Char */ _data_type,
      /* c_Char */ _data_source,
      /* c_Bool& */ do_one_turn_map);
}
bool Tao::tao_lat_bookkeeper(TaoUniverseProxy& u, TaoLatticeProxy& tao_lat) {
  bool _err_flag{};
  fortran_tao_lat_bookkeeper(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ tao_lat.get_fortran_ptr(),
      /* c_Bool& */ _err_flag);
  return _err_flag;
}
void Tao::tao_setup_key_table() {
  fortran_tao_setup_key_table();
}
void Tao::tao_abort_command_file(std::optional<bool> force_abort) {
  bool force_abort_lvalue;
  auto* _force_abort{&force_abort_lvalue};
  if (force_abort.has_value()) {
    force_abort_lvalue = force_abort.value();
  } else {
    _force_abort = nullptr;
  }
  fortran_tao_abort_command_file(/* c_Bool* */ _force_abort);
}
void Tao::tao_data_check(bool err) {
  fortran_tao_data_check(/* c_Bool& */ err);
}
void Tao::tao_beam_emit_calc(
    int plane,
    int emit_type,
    EleProxy& ele,
    BunchParamsProxy& bunch_params,
    double emit) {
  fortran_tao_beam_emit_calc(
      /* c_Int& */ plane,
      /* c_Int& */ emit_type,
      /* void* */ ele.get_fortran_ptr(),
      /* void* */ bunch_params.get_fortran_ptr(),
      /* c_Real& */ emit);
}
void Tao::tao_clip_cmd(
    bool gang,
    std::string where,
    double value1,
    double value2) {
  auto _where = where.c_str();
  fortran_tao_clip_cmd(
      /* c_Bool& */ gang,
      /* c_Char */ _where,
      /* c_Real& */ value1,
      /* c_Real& */ value2);
}
void Tao::tao_fixer(std::string switch_, std::string word1, std::string word2) {
  auto _switch_ = switch_.c_str();
  auto _word1 = word1.c_str();
  auto _word2 = word2.c_str();
  fortran_tao_fixer(
      /* c_Char */ _switch_, /* c_Char */ _word1, /* c_Char */ _word2);
}
bool Tao::tao_init() {
  bool _err_flag{};
  fortran_tao_init(/* c_Bool& */ _err_flag);
  return _err_flag;
}
void Tao::tao_key_info_to_str(
    int ix_key,
    int ix_min_key,
    int ix_max_key,
    std::string key_str,
    std::string header_str) {
  auto _key_str = key_str.c_str(); // ptr, inout, required
  auto _header_str = header_str.c_str(); // ptr, inout, required
  fortran_tao_key_info_to_str(
      /* c_Int& */ ix_key,
      /* c_Int& */ ix_min_key,
      /* c_Int& */ ix_max_key,
      /* c_Char */ _key_str,
      /* c_Char */ _header_str);
}
bool Tao::tao_lmdif_optimizer() {
  bool _abort{};
  fortran_tao_lmdif_optimizer(/* c_Bool& */ _abort);
  return _abort;
}
void Tao::tao_place_cmd(
    std::string where,
    std::string who,
    std::optional<bool> no_buffer) {
  auto _where = where.c_str();
  auto _who = who.c_str();
  bool no_buffer_lvalue;
  auto* _no_buffer{&no_buffer_lvalue};
  if (no_buffer.has_value()) {
    no_buffer_lvalue = no_buffer.value();
  } else {
    _no_buffer = nullptr;
  }
  fortran_tao_place_cmd(
      /* c_Char */ _where, /* c_Char */ _who, /* c_Bool* */ _no_buffer);
}
void Tao::tao_print_command_line_info() {
  fortran_tao_print_command_line_info();
}
void Tao::tao_remove_blank_characters(std::string str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_tao_remove_blank_characters(/* c_Char */ _str);
}
void Tao::tao_read_phase_space_index(
    std::string name,
    int ixc,
    std::optional<bool> print_err,
    int ix_ps) {
  auto _name = name.c_str();
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_read_phase_space_index(
      /* c_Char */ _name,
      /* c_Int& */ ixc,
      /* c_Bool* */ _print_err,
      /* c_Int& */ ix_ps);
}
bool Tao::tao_run_cmd(std::string which) {
  auto _which = which.c_str();
  bool _abort{};
  fortran_tao_run_cmd(/* c_Char */ _which, /* c_Bool& */ _abort);
  return _abort;
}
void Tao::tao_set_data_useit_opt(optional_ref<TaoDataProxyAlloc1D> data) {
  // intent=in allocatable type array
  auto* _data = data.has_value() ? data->get().get_fortran_ptr()
                                 : nullptr; // input, optional
  fortran_tao_set_data_useit_opt(/* void* */ _data);
}
void Tao::tao_single_mode(std::string char_) {
  auto _char_ = char_.c_str();
  fortran_tao_single_mode(/* c_Char */ _char_);
}
void Tao::tao_taper_cmd(std::string except, std::string uni_names) {
  auto _except = except.c_str();
  auto _uni_names = uni_names.c_str();
  fortran_tao_taper_cmd(/* c_Char */ _except, /* c_Char */ _uni_names);
}
bool Tao::tao_user_is_terminating_optimization() {
  bool _is_terminating{};
  fortran_tao_user_is_terminating_optimization(/* c_Bool& */ _is_terminating);
  return _is_terminating;
}
int Tao::tao_uni_atsign_index(std::string string) {
  auto _string = string.c_str();
  int _ix_amp{};
  fortran_tao_uni_atsign_index(/* c_Char */ _string, /* c_Int& */ _ix_amp);
  return _ix_amp;
}
bool Tao::tao_dmodel_dvar_calc(bool force_calc) {
  bool _err_flag{};
  fortran_tao_dmodel_dvar_calc(
      /* c_Bool& */ force_calc, /* c_Bool& */ _err_flag);
  return _err_flag;
}
void Tao::tao_veto_vars_with_zero_dmodel() {
  fortran_tao_veto_vars_with_zero_dmodel();
}
void Tao::tao_dmerit_calc() {
  fortran_tao_dmerit_calc();
}
void Tao::tao_init_global(std::string init_file) {
  auto _init_file = init_file.c_str();
  fortran_tao_init_global(/* c_Char */ _init_file);
}
void Tao::tao_init_beams(std::string init_file) {
  auto _init_file = init_file.c_str();
  fortran_tao_init_beams(/* c_Char */ _init_file);
}
void Tao::tao_init_beam_in_universe(
    TaoUniverseProxy& u,
    BeamInitProxy& beam_init,
    std::string track_start,
    std::string track_end,
    double comb_ds_save) {
  auto _track_start = track_start.c_str(); // ptr, inout, required
  auto _track_end = track_end.c_str(); // ptr, inout, required
  fortran_tao_init_beam_in_universe(
      /* void* */ u.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* c_Char */ _track_start,
      /* c_Char */ _track_end,
      /* c_Real& */ comb_ds_save);
}
void Tao::tao_init_dynamic_aperture(std::string init_file) {
  auto _init_file = init_file.c_str();
  fortran_tao_init_dynamic_aperture(/* c_Char */ _init_file);
}
bool Tao::tao_svd_optimizer() {
  bool _abort{};
  fortran_tao_svd_optimizer(/* c_Bool& */ _abort);
  return _abort;
}
bool Tao::tao_change_tune(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string dqa_str,
    std::string dqb_str) {
  auto _branch_str = branch_str.c_str();
  auto _mask_str = mask_str.c_str();
  auto _dqa_str = dqa_str.c_str();
  auto _dqb_str = dqb_str.c_str();
  bool _err_flag{};
  fortran_tao_change_tune(
      /* c_Char */ _branch_str,
      /* c_Char */ _mask_str,
      /* c_Bool& */ print_list,
      /* c_Char */ _dqa_str,
      /* c_Char */ _dqb_str,
      /* c_Bool& */ _err_flag);
  return _err_flag;
}
bool Tao::tao_change_z_tune(std::string branch_str, std::string dq_str) {
  auto _branch_str = branch_str.c_str();
  auto _dq_str = dq_str.c_str();
  bool _err_flag{};
  fortran_tao_change_z_tune(
      /* c_Char */ _branch_str, /* c_Char */ _dq_str, /* c_Bool& */ _err_flag);
  return _err_flag;
}
bool Tao::tao_change_var(std::string name, std::string num_str, bool silent) {
  auto _name = name.c_str();
  auto _num_str = num_str.c_str();
  bool _err_flag{};
  fortran_tao_change_var(
      /* c_Char */ _name,
      /* c_Char */ _num_str,
      /* c_Bool& */ silent,
      /* c_Bool& */ _err_flag);
  return _err_flag;
}
bool Tao::tao_change_ele(
    std::string ele_name,
    std::string attrib_name,
    std::string num_str,
    bool update) {
  auto _ele_name = ele_name.c_str();
  auto _attrib_name = attrib_name.c_str();
  auto _num_str = num_str.c_str();
  bool _err_flag{};
  fortran_tao_change_ele(
      /* c_Char */ _ele_name,
      /* c_Char */ _attrib_name,
      /* c_Char */ _num_str,
      /* c_Bool& */ update,
      /* c_Bool& */ _err_flag);
  return _err_flag;
}
void Tao::tao_to_change_number(
    std::string num_str,
    int n_size,
    RealAlloc1D& change_number,
    std::string abs_or_rel,
    bool err) {
  auto _num_str = num_str.c_str(); // ptr, inout, required
  // intent=inout allocatable general array
  auto _abs_or_rel = abs_or_rel.c_str(); // ptr, inout, required
  fortran_tao_to_change_number(
      /* c_Char */ _num_str,
      /* c_Int& */ n_size,
      /* void* */ change_number.get_fortran_ptr(),
      /* c_Char */ _abs_or_rel,
      /* c_Bool& */ err);
}
void Tao::tao_set_tune_cmd(
    std::string branch_str,
    std::string mask_str,
    bool print_list,
    std::string qa_str,
    std::string qb_str,
    bool delta_input) {
  auto _branch_str = branch_str.c_str();
  auto _mask_str = mask_str.c_str();
  auto _qa_str = qa_str.c_str();
  auto _qb_str = qb_str.c_str();
  fortran_tao_set_tune_cmd(
      /* c_Char */ _branch_str,
      /* c_Char */ _mask_str,
      /* c_Bool& */ print_list,
      /* c_Char */ _qa_str,
      /* c_Char */ _qb_str,
      /* c_Bool& */ delta_input);
}
void Tao::tao_set_z_tune_cmd(
    std::string branch_str,
    std::string q_str,
    bool delta_input) {
  auto _branch_str = branch_str.c_str();
  auto _q_str = q_str.c_str();
  fortran_tao_set_z_tune_cmd(
      /* c_Char */ _branch_str, /* c_Char */ _q_str, /* c_Bool& */ delta_input);
}
void Tao::tao_set_openmp_n_threads(int n_threads) {
  fortran_tao_set_openmp_n_threads(/* c_Int& */ n_threads);
}
void Tao::tao_set_calculate_cmd(optional_ref<std::string> switch_) {
  const char* _switch_ = switch_.has_value() ? switch_->get().c_str() : nullptr;
  fortran_tao_set_calculate_cmd(/* c_Char */ _switch_);
}
void Tao::tao_set_key_cmd(std::string key_str, std::string cmd_str) {
  auto _key_str = key_str.c_str();
  auto _cmd_str = cmd_str.c_str();
  fortran_tao_set_key_cmd(/* c_Char */ _key_str, /* c_Char */ _cmd_str);
}
void Tao::tao_set_ran_state_cmd(std::string state_string) {
  auto _state_string = state_string.c_str();
  fortran_tao_set_ran_state_cmd(/* c_Char */ _state_string);
}
void Tao::tao_set_lattice_cmd(std::string dest_lat, std::string source_lat) {
  auto _dest_lat = dest_lat.c_str();
  auto _source_lat = source_lat.c_str();
  fortran_tao_set_lattice_cmd(/* c_Char */ _dest_lat, /* c_Char */ _source_lat);
}
void Tao::tao_set_global_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_global_cmd(/* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_space_charge_com_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_space_charge_com_cmd(
      /* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_bmad_com_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_bmad_com_cmd(/* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_ptc_com_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_ptc_com_cmd(/* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_geodesic_lm_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_geodesic_lm_cmd(/* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_opti_de_param_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_opti_de_param_cmd(/* c_Char */ _who, /* c_Char */ _value_str);
}
bool Tao::tao_set_wave_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  bool _err{};
  fortran_tao_set_wave_cmd(
      /* c_Char */ _who, /* c_Char */ _value_str, /* c_Bool& */ _err);
  return _err;
}
void Tao::tao_set_beam_cmd(
    std::string who,
    std::string value_str,
    std::string branch_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  auto _branch_str = branch_str.c_str();
  fortran_tao_set_beam_cmd(
      /* c_Char */ _who, /* c_Char */ _value_str, /* c_Char */ _branch_str);
}
void Tao::tao_set_beam_init_cmd(
    std::string who,
    std::string value_str,
    std::string branch_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  auto _branch_str = branch_str.c_str();
  fortran_tao_set_beam_init_cmd(
      /* c_Char */ _who, /* c_Char */ _value_str, /* c_Char */ _branch_str);
}
void Tao::tao_set_particle_start_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_particle_start_cmd(
      /* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_plot_page_cmd(
    std::string component,
    std::string value_str,
    std::optional<std::string> value_str2) {
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  const char* _value_str2 =
      value_str2.has_value() ? value_str2->c_str() : nullptr;
  fortran_tao_set_plot_page_cmd(
      /* c_Char */ _component,
      /* c_Char */ _value_str,
      /* c_Char */ _value_str2);
}
void Tao::tao_set_curve_cmd(
    std::string curve_name,
    std::string component,
    std::string value_str) {
  auto _curve_name = curve_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_curve_cmd(
      /* c_Char */ _curve_name,
      /* c_Char */ _component,
      /* c_Char */ _value_str);
}
void Tao::tao_set_plot_cmd(
    std::string plot_name,
    std::string component,
    std::string value_str) {
  auto _plot_name = plot_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_plot_cmd(
      /* c_Char */ _plot_name,
      /* c_Char */ _component,
      /* c_Char */ _value_str);
}
void Tao::tao_set_region_cmd(
    std::string region_name,
    std::string component,
    std::string value_str) {
  auto _region_name = region_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_region_cmd(
      /* c_Char */ _region_name,
      /* c_Char */ _component,
      /* c_Char */ _value_str);
}
void Tao::tao_set_graph_cmd(
    std::string graph_name,
    std::string component,
    std::string value_str) {
  auto _graph_name = graph_name.c_str();
  auto _component = component.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_graph_cmd(
      /* c_Char */ _graph_name,
      /* c_Char */ _component,
      /* c_Char */ _value_str);
}
void Tao::tao_set_var_cmd(std::string var_str, std::string value_str) {
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_var_cmd(/* c_Char */ _var_str, /* c_Char */ _value_str);
}
void Tao::tao_set_branch_cmd(
    std::string branch_str,
    std::string component_str,
    std::string value_str) {
  auto _branch_str = branch_str.c_str();
  auto _component_str = component_str.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_branch_cmd(
      /* c_Char */ _branch_str,
      /* c_Char */ _component_str,
      /* c_Char */ _value_str);
}
void Tao::tao_set_data_cmd(
    std::string who_str,
    std::string value_str,
    optional_ref<bool> silent) {
  auto _who_str = who_str.c_str();
  auto _value_str = value_str.c_str();
  auto* _silent =
      silent.has_value() ? &silent->get() : nullptr; // inout, optional
  fortran_tao_set_data_cmd(
      /* c_Char */ _who_str, /* c_Char */ _value_str, /* c_Bool* */ _silent);
}
void Tao::tao_set_default_cmd(std::string who_str, std::string value_str) {
  auto _who_str = who_str.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_default_cmd(/* c_Char */ _who_str, /* c_Char */ _value_str);
}
void Tao::tao_set_dynamic_aperture_cmd(std::string who, std::string value_str) {
  auto _who = who.c_str();
  auto _value_str = value_str.c_str();
  fortran_tao_set_dynamic_aperture_cmd(
      /* c_Char */ _who, /* c_Char */ _value_str);
}
void Tao::tao_set_universe_cmd(
    std::string uni,
    std::string who,
    std::string what) {
  auto _uni = uni.c_str();
  auto _who = who.c_str();
  auto _what = what.c_str();
  fortran_tao_set_universe_cmd(
      /* c_Char */ _uni, /* c_Char */ _who, /* c_Char */ _what);
}
void Tao::tao_set_elements_cmd(
    std::string ele_list,
    std::string attribute,
    std::string value,
    bool update) {
  auto _ele_list = ele_list.c_str();
  auto _attribute = attribute.c_str();
  auto _value = value.c_str();
  fortran_tao_set_elements_cmd(
      /* c_Char */ _ele_list,
      /* c_Char */ _attribute,
      /* c_Char */ _value,
      /* c_Bool& */ update);
}
Tao::TaoSetLogicalValue Tao::tao_set_logical_value(
    std::string var_str,
    std::string value_str) {
  bool _var{};
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  bool _error{};
  fortran_tao_set_logical_value(
      /* c_Bool& */ _var,
      /* c_Char */ _var_str,
      /* c_Char */ _value_str,
      /* c_Bool& */ _error);
  return TaoSetLogicalValue{_var, _error};
}
Tao::TaoSetIntegerValue Tao::tao_set_integer_value(
    std::string var_str,
    std::string value_str,
    std::optional<int> min_val,
    std::optional<int> max_val,
    std::optional<bool> print_err) {
  int _var{};
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  bool _error{};
  int min_val_lvalue;
  auto* _min_val{&min_val_lvalue};
  if (min_val.has_value()) {
    min_val_lvalue = min_val.value();
  } else {
    _min_val = nullptr;
  }
  int max_val_lvalue;
  auto* _max_val{&max_val_lvalue};
  if (max_val.has_value()) {
    max_val_lvalue = max_val.value();
  } else {
    _max_val = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_tao_set_integer_value(
      /* c_Int& */ _var,
      /* c_Char */ _var_str,
      /* c_Char */ _value_str,
      /* c_Bool& */ _error,
      /* c_Int* */ _min_val,
      /* c_Int* */ _max_val,
      /* c_Bool* */ _print_err);
  return TaoSetIntegerValue{_var, _error};
}
Tao::TaoSetRealValue Tao::tao_set_real_value(
    std::string var_str,
    std::string value_str,
    std::optional<double> min_val,
    std::optional<double> max_val,
    std::optional<int> dflt_uni) {
  double _var{};
  auto _var_str = var_str.c_str();
  auto _value_str = value_str.c_str();
  bool _error{};
  double min_val_lvalue;
  auto* _min_val{&min_val_lvalue};
  if (min_val.has_value()) {
    min_val_lvalue = min_val.value();
  } else {
    _min_val = nullptr;
  }
  double max_val_lvalue;
  auto* _max_val{&max_val_lvalue};
  if (max_val.has_value()) {
    max_val_lvalue = max_val.value();
  } else {
    _max_val = nullptr;
  }
  int dflt_uni_lvalue;
  auto* _dflt_uni{&dflt_uni_lvalue};
  if (dflt_uni.has_value()) {
    dflt_uni_lvalue = dflt_uni.value();
  } else {
    _dflt_uni = nullptr;
  }
  fortran_tao_set_real_value(
      /* c_Real& */ _var,
      /* c_Char */ _var_str,
      /* c_Char */ _value_str,
      /* c_Bool& */ _error,
      /* c_Real* */ _min_val,
      /* c_Real* */ _max_val,
      /* c_Int* */ _dflt_uni);
  return TaoSetRealValue{_var, _error};
}
void Tao::tao_set_symbolic_number_cmd(
    std::string sym_str,
    std::optional<std::string> num_str,
    std::optional<double> val) {
  auto _sym_str = sym_str.c_str();
  const char* _num_str = num_str.has_value() ? num_str->c_str() : nullptr;
  double val_lvalue;
  auto* _val{&val_lvalue};
  if (val.has_value()) {
    val_lvalue = val.value();
  } else {
    _val = nullptr;
  }
  fortran_tao_set_symbolic_number_cmd(
      /* c_Char */ _sym_str, /* c_Char */ _num_str, /* c_Real* */ _val);
}