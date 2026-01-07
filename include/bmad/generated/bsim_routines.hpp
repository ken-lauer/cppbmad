#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"

using namespace Bmad;

namespace bsim {
extern "C" void fortran_bbu_add_a_bunch(
    void* lat /* 0D_NOT_type inout */,
    void* bbu_beam /* 0D_NOT_type inout */,
    void* bbu_param /* 0D_NOT_type inout */,
    void* beam_init /* 0D_NOT_type inout */);
void bbu_add_a_bunch(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    BbuParamProxy& bbu_param,
    BeamInitProxy& beam_init);
extern "C" void fortran_bbu_hom_voltage_calc(
    void* lat /* 0D_NOT_type inout */,
    void* bbu_beam /* 0D_NOT_type inout */,
    int& n_period /* 0D_NOT_integer inout */,
    int& ix_stage_last_tracked /* 0D_NOT_integer inout */);
void bbu_hom_voltage_calc(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    int& n_period,
    int& ix_stage_last_tracked);
extern "C" void fortran_bbu_remove_head_bunch(
    void* bbu_beam /* 0D_NOT_type inout */);
void bbu_remove_head_bunch(BbuBeamProxy& bbu_beam);
extern "C" void fortran_bbu_setup(
    void* lat /* 0D_NOT_type inout */,
    double& dt_bunch /* 0D_NOT_real inout */,
    void* bbu_param /* 0D_NOT_type inout */,
    void* bbu_beam /* 0D_NOT_type inout */);
void bbu_setup(
    LatProxy& lat,
    double& dt_bunch,
    BbuParamProxy& bbu_param,
    BbuBeamProxy& bbu_beam);
extern "C" void fortran_bbu_track_a_stage(
    void* lat /* 0D_NOT_type inout */,
    void* bbu_beam /* 0D_NOT_type inout */,
    void* bbu_param /* 0D_NOT_type inout */,
    bool& lost /* 0D_NOT_logical inout */,
    int& ix_stage_tracked /* 0D_NOT_integer inout */);
void bbu_track_a_stage(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    BbuParamProxy& bbu_param,
    bool& lost,
    int& ix_stage_tracked);
extern "C" void fortran_bbu_track_all(
    void* lat /* 0D_NOT_type inout */,
    void* bbu_beam /* 0D_NOT_type inout */,
    void* bbu_param /* 0D_NOT_type inout */,
    void* beam_init /* 0D_NOT_type inout */,
    double& hom_voltage_normalized /* 0D_NOT_real inout */,
    double& growth_rate /* 0D_NOT_real inout */,
    bool& lost /* 0D_NOT_logical inout */,
    int& irep /* 0D_NOT_integer inout */);
void bbu_track_all(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    BbuParamProxy& bbu_param,
    BeamInitProxy& beam_init,
    double& hom_voltage_normalized,
    double& growth_rate,
    bool& lost,
    int& irep);
extern "C" void fortran_check_rf_freq(
    void* lat /* 0D_NOT_type inout */,
    double& fb /* 0D_NOT_real inout */);
void check_rf_freq(LatProxy& lat, double& fb);
extern "C" void fortran_count_lines_in_file(
    const char* file_name /* 0D_NOT_character in */,
    int& lines /* 0D_NOT_integer out */);
int count_lines_in_file(std::string file_name);
extern "C" bool fortran_hom_voltage(
    void* lr_wake /* 0D_NOT_type inout */,
    double& voltage /* 0D_NOT_real inout */);
void hom_voltage(WakeLrModeProxy& lr_wake, double& voltage);
extern "C" void fortran_insert_phase_trombone(
    void* branch /* 0D_NOT_type inout */);
void insert_phase_trombone(BranchProxy& branch);
extern "C" bool fortran_logical_to_python(
    bool& logic /* 0D_NOT_logical inout */,
    const char* string /* 0D_NOT_character inout */);
void logical_to_python(bool& logic, std::string& string);

// Skipped unusable routine longitudinal_beta:
// - Module name unset
extern "C" void fortran_rf_cav_names(void* lat /* 0D_NOT_type inout */);
void rf_cav_names(LatProxy& lat);
extern "C" bool fortran_set_tune_3d(
    void* branch /* 0D_NOT_type inout */,
    double* target_tunes /* 1D_NOT_real in */,
    const char* mask /* 0D_NOT_character inout */,
    bool* use_phase_trombone /* 0D_NOT_logical in */,
    bool* z_tune_set /* 0D_NOT_logical in */,
    const char** group_knobs /* 1D_NOT_character in */,
    bool* print_err /* 0D_NOT_logical in */,
    bool& everything_ok /* 0D_NOT_logical inout */);
void set_tune_3d(
    BranchProxy& branch,
    FixedArray1D<Real, 3> target_tunes,
    optional_ref<std::string> mask,
    std::optional<bool> use_phase_trombone,
    std::optional<bool> z_tune_set,
    std::optional<FixedArray1D<string, 2>> group_knobs,
    std::optional<bool> print_err,
    bool& everything_ok);

// Skipped unusable routine track_s_to_s:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_write_bunch_by_bunch_info(
    void* lat /* 0D_NOT_type inout */,
    void* bbu_beam /* 0D_NOT_type inout */,
    void* bbu_param /* 0D_NOT_type inout */,
    void* this_stage /* 0D_PTR_type inout */);
void write_bunch_by_bunch_info(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    BbuParamProxy& bbu_param,
    BbuStageProxy& this_stage);
} // namespace bsim
