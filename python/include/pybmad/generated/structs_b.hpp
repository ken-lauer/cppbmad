#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_bbu_beam_struct(py::module &m, py::class_<BbuBeamStruct> &class_);
void init_bbu_param_struct(py::module &m, py::class_<BbuParamStruct> &class_);
void init_bbu_stage_struct(py::module &m, py::class_<BbuStageStruct> &class_);
void init_beam_init_struct(py::module &m, py::class_<BeamInitStruct> &class_);
void init_beam_struct(py::module &m, py::class_<BeamStruct> &class_);
void init_bmad_common_struct(py::module &m, py::class_<BmadCommonStruct> &class_);
void init_bmad_normal_form_struct(py::module &m, py::class_<BmadNormalFormStruct> &class_);
void init_bookkeeping_state_struct(py::module &m, py::class_<BookkeepingStateStruct> &class_);
void init_bpm_phase_coupling_struct(py::module &m, py::class_<BpmPhaseCouplingStruct> &class_);
void init_branch_struct(py::module &m, py::class_<BranchStruct> &class_);
void init_bunch_params_struct(py::module &m, py::class_<BunchParamsStruct> &class_);
void init_bunch_struct(py::module &m, py::class_<BunchStruct> &class_);
void init_bunch_track_struct(py::module &m, py::class_<BunchTrackStruct> &class_);
void init_bicubic_cmplx_coef_struct(py::module &m, py::class_<BicubicCmplxCoefStruct> &class_);
