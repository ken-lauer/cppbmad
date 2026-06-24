#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

namespace bsim {
extern "C" void fortran_bbu_add_a_bunch(
    void *lat /* 0D_NOT_type inout */,
    void *bbu_beam /* 0D_NOT_type inout */,
    void *bbu_param /* 0D_NOT_type inout */,
    void *beam_init /* 0D_NOT_type inout */
);
void bbu_add_a_bunch(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    BeamInitStruct &beam_init
);
extern "C" void fortran_bbu_hom_voltage_calc(
    void *lat /* 0D_NOT_type inout */,
    void *bbu_beam /* 0D_NOT_type inout */,
    int &n_period /* 0D_NOT_integer in */,
    int &ix_stage_last_tracked /* 0D_NOT_integer in */
);
void bbu_hom_voltage_calc(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    int n_period,
    int ix_stage_last_tracked
);
extern "C" void fortran_bbu_remove_head_bunch(void *bbu_beam /* 0D_NOT_type inout */);
void bbu_remove_head_bunch(BbuBeamStruct &bbu_beam);
extern "C" void fortran_bbu_setup(
    void *lat /* 0D_NOT_type inout */,
    double &dt_bunch /* 0D_NOT_real in */,
    void *bbu_param /* 0D_NOT_type inout */,
    void *bbu_beam /* 0D_NOT_type inout */
);
void bbu_setup(LatStruct &lat, double dt_bunch, BbuParamStruct &bbu_param, BbuBeamStruct &bbu_beam);
extern "C" void fortran_bbu_track_a_stage(
    void *lat /* 0D_NOT_type inout */,
    void *bbu_beam /* 0D_NOT_type inout */,
    void *bbu_param /* 0D_NOT_type inout */,
    bool &lost /* 0D_NOT_logical in */,
    int &ix_stage_tracked /* 0D_NOT_integer in */
);
void bbu_track_a_stage(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    bool lost,
    int ix_stage_tracked
);
extern "C" void fortran_bbu_track_all(
    void *lat /* 0D_NOT_type inout */,
    void *bbu_beam /* 0D_NOT_type inout */,
    void *bbu_param /* 0D_NOT_type inout */,
    void *beam_init /* 0D_NOT_type inout */,
    double &hom_voltage_normalized /* 0D_NOT_real out */,
    double &growth_rate /* 0D_NOT_real out */,
    bool &lost /* 0D_NOT_logical out */,
    int &irep /* 0D_NOT_integer out */
);
struct BbuTrackAll {
  double hom_voltage_normalized;
  double growth_rate;
  bool lost;
  int irep;
};
bsim::BbuTrackAll bbu_track_all(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    BeamInitStruct &beam_init
);
extern "C" void
fortran_check_rf_freq(void *lat /* 0D_NOT_type inout */, double &fb /* 0D_NOT_real in */);
void check_rf_freq(LatStruct &lat, double fb);
extern "C" void fortran_count_lines_in_file(
    const char *file_name /* 0D_NOT_character in */,
    int &lines /* 0D_NOT_integer out */
);
int count_lines_in_file(std::string file_name);
extern "C" bool
fortran_hom_voltage(void *lr_wake /* 0D_NOT_type inout */, double &voltage /* 0D_NOT_real out */);
double hom_voltage(WakeLrModeStruct &lr_wake);
extern "C" void fortran_insert_phase_trombone(void *branch /* 0D_NOT_type inout */);
void insert_phase_trombone(BranchStruct &branch);
extern "C" bool fortran_logical_to_python(
    bool &logic /* 0D_NOT_logical in */,
    const char *string /* 0D_NOT_character out */
);
std::string logical_to_python(bool logic);

// Skipped unusable routine longitudinal_beta:
// - Module name unset
extern "C" void fortran_rf_cav_names(void *lat /* 0D_NOT_type inout */);
void rf_cav_names(LatStruct &lat);
extern "C" bool fortran_set_tune_3d(
    void *branch /* 0D_NOT_type inout */,
    Bmad::array_descriptor_t &target_tunes /* 1D_NOT_real in */,
    const char *mask /* 0D_NOT_character in */,
    bool *use_phase_trombone /* 0D_NOT_logical in */,
    bool *z_tune_set /* 0D_NOT_logical in */,
    const char **group_knobs /* 1D_NOT_character in */,
    bool *print_err /* 0D_NOT_logical in */,
    bool &everything_ok /* 0D_NOT_logical out */
);
bool set_tune_3d(
    BranchStruct &branch,
    FixedArray1D<Real, 3> target_tunes,
    std::optional<std::string> mask = std::nullopt,
    std::optional<bool> use_phase_trombone = std::nullopt,
    std::optional<bool> z_tune_set = std::nullopt,
    std::optional<FixedArray1D<string, 2>> group_knobs = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);

// Skipped unusable routine track_s_to_s:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_write_bunch_by_bunch_info(
    void *lat /* 0D_NOT_type inout */,
    void *bbu_beam /* 0D_NOT_type inout */,
    void *bbu_param /* 0D_NOT_type inout */,
    void *this_stage /* 0D_PTR_type inout */
);
void write_bunch_by_bunch_info(
    LatStruct &lat,
    BbuBeamStruct &bbu_beam,
    BbuParamStruct &bbu_param,
    BbuStageStruct &this_stage
);
} // namespace bsim
