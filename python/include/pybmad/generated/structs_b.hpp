#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_bbu_beam_struct(nb::module_ &m, nb::class_<BbuBeamStruct> &class_);
void init_bbu_param_struct(nb::module_ &m, nb::class_<BbuParamStruct> &class_);
void init_bbu_stage_struct(nb::module_ &m, nb::class_<BbuStageStruct> &class_);
void init_beam_init_struct(nb::module_ &m, nb::class_<BeamInitStruct> &class_);
void init_beam_struct(nb::module_ &m, nb::class_<BeamStruct> &class_);
void init_bmad_common_struct(nb::module_ &m, nb::class_<BmadCommonStruct> &class_);
void init_bmad_normal_form_struct(nb::module_ &m, nb::class_<BmadNormalFormStruct> &class_);
void init_bookkeeping_state_struct(nb::module_ &m, nb::class_<BookkeepingStateStruct> &class_);
void init_bpm_phase_coupling_struct(nb::module_ &m, nb::class_<BpmPhaseCouplingStruct> &class_);
void init_branch_struct(nb::module_ &m, nb::class_<BranchStruct> &class_);
void init_bunch_params_struct(nb::module_ &m, nb::class_<BunchParamsStruct> &class_);
void init_bunch_struct(nb::module_ &m, nb::class_<BunchStruct> &class_);
void init_bunch_track_struct(nb::module_ &m, nb::class_<BunchTrackStruct> &class_);
void init_bicubic_cmplx_coef_struct(nb::module_ &m, nb::class_<BicubicCmplxCoefStruct> &class_);
