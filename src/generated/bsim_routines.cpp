#include "bmad/generated/bsim_routines.hpp"

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
void bsim::bbu_add_a_bunch(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    BeamInitStruct &beam_init
) {
  fortran_bbu_add_a_bunch(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ bbu_beam.get_fortran_ptr(),
      /* void* */ bbu_param.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr()
  );
}
void bsim::bbu_hom_voltage_calc(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    int n_period,
    int ix_stage_last_tracked
) {
  fortran_bbu_hom_voltage_calc(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ bbu_beam.get_fortran_ptr(),
      /* int& */ n_period,
      /* int& */ ix_stage_last_tracked
  );
}
void bsim::bbu_remove_head_bunch(BbuBeamStruct &bbu_beam) {
  fortran_bbu_remove_head_bunch(/* void* */ bbu_beam.get_fortran_ptr());
}
void bsim::bbu_setup(
    LatStruct &lat,
    double dt_bunch,
    BbuParamStruct &bbu_param,
    BbuBeamStruct &bbu_beam
) {
  fortran_bbu_setup(
      /* void* */ lat.get_fortran_ptr(),
      /* double& */ dt_bunch,
      /* void* */ bbu_param.get_fortran_ptr(),
      /* void* */ bbu_beam.get_fortran_ptr()
  );
}
void bsim::bbu_track_a_stage(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    bool lost,
    int ix_stage_tracked
) {
  fortran_bbu_track_a_stage(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ bbu_beam.get_fortran_ptr(),
      /* void* */ bbu_param.get_fortran_ptr(),
      /* bool& */ lost,
      /* int& */ ix_stage_tracked
  );
}
void bsim::bbu_track_all(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    BeamInitStruct &beam_init,
    double hom_voltage_normalized,
    double growth_rate,
    bool lost,
    int irep
) {
  fortran_bbu_track_all(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ bbu_beam.get_fortran_ptr(),
      /* void* */ bbu_param.get_fortran_ptr(),
      /* void* */ beam_init.get_fortran_ptr(),
      /* double& */ hom_voltage_normalized,
      /* double& */ growth_rate,
      /* bool& */ lost,
      /* int& */ irep
  );
}
void bsim::check_rf_freq(LatStruct &lat, double fb) {
  fortran_check_rf_freq(/* void* */ lat.get_fortran_ptr(), /* double& */ fb);
}
int bsim::count_lines_in_file(std::string file_name) {
  auto _file_name = file_name.c_str();
  int _lines{};
  fortran_count_lines_in_file(/* const char* */ _file_name, /* int& */ _lines);
  return _lines;
}
void bsim::hom_voltage(WakeLrModeStruct &lr_wake, double voltage) {
  fortran_hom_voltage(/* void* */ lr_wake.get_fortran_ptr(), /* double& */ voltage);
}
void bsim::insert_phase_trombone(BranchStruct &branch) {
  fortran_insert_phase_trombone(/* void* */ branch.get_fortran_ptr());
}
void bsim::logical_to_python(bool logic, std::string string) {
  auto _string = string.c_str();
  fortran_logical_to_python(/* bool& */ logic, /* const char* */ _string);
}
void bsim::rf_cav_names(LatStruct &lat) { fortran_rf_cav_names(/* void* */ lat.get_fortran_ptr()); }
void bsim::set_tune_3d(
    BranchStruct &branch,
    FixedArray1D<Real, 3> target_tunes,
    std::optional<std::string> mask,
    std::optional<bool> use_phase_trombone,
    std::optional<bool> z_tune_set,
    std::optional<FixedArray1D<string, 2>> group_knobs,
    std::optional<bool> print_err,
    bool everything_ok
) {
  // target_tunes: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _target_tunes_desc;
  _target_tunes_desc.rank = 1;
  _target_tunes_desc.data_ptr = target_tunes.data();
  _target_tunes_desc.dims[0] = target_tunes.size();
  const char *_mask = mask.has_value() ? mask->c_str() : nullptr;
  bool use_phase_trombone_lvalue;
  auto *_use_phase_trombone{&use_phase_trombone_lvalue};
  if (use_phase_trombone.has_value()) {
    use_phase_trombone_lvalue = use_phase_trombone.value();
  } else {
    _use_phase_trombone = nullptr;
  }
  bool z_tune_set_lvalue;
  auto *_z_tune_set{&z_tune_set_lvalue};
  if (z_tune_set.has_value()) {
    z_tune_set_lvalue = z_tune_set.value();
  } else {
    _z_tune_set = nullptr;
  }
  // Optional string array
  std::vector<const char *> _group_knobs{group_knobs.has_value() ? group_knobs->size() : 0};
  for (size_t i{0}; i < _group_knobs.size(); i++) {
    _group_knobs.push_back(group_knobs->data()[i].data());
  }

  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_set_tune_3d(
      /* void* */ branch.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _target_tunes_desc,
      /* const char* */ _mask,
      /* bool* */ _use_phase_trombone,
      /* bool* */ _z_tune_set,
      /* const char** */ _group_knobs.size() ? _group_knobs.data() : nullptr,
      /* bool* */ _print_err,
      /* bool& */ everything_ok
  );
}
void bsim::write_bunch_by_bunch_info(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    BbuStageStruct &this_stage
) {
  auto _this_stage = &this_stage; // input, required, pointer
  fortran_write_bunch_by_bunch_info(
      /* void* */ lat.get_fortran_ptr(),
      /* void* */ bbu_beam.get_fortran_ptr(),
      /* void* */ bbu_param.get_fortran_ptr(),
      /* void* */ &this_stage
  );
}