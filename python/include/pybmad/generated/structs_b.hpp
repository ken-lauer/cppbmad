#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_bbu_beam_struct(py::module& m, py::class_<BbuBeamProxy>& class_);
void init_bbu_param_struct(py::module& m, py::class_<BbuParamProxy>& class_);
void init_bbu_stage_struct(py::module& m, py::class_<BbuStageProxy>& class_);
void init_beam_init_struct(py::module& m, py::class_<BeamInitProxy>& class_);
void init_beam_struct(py::module& m, py::class_<BeamProxy>& class_);
void init_bmad_common_struct(
    py::module& m,
    py::class_<BmadCommonProxy>& class_);
void init_bmad_normal_form_struct(
    py::module& m,
    py::class_<BmadNormalFormProxy>& class_);
void init_bookkeeping_state_struct(
    py::module& m,
    py::class_<BookkeepingStateProxy>& class_);
void init_bpm_phase_coupling_struct(
    py::module& m,
    py::class_<BpmPhaseCouplingProxy>& class_);
void init_branch_struct(py::module& m, py::class_<BranchProxy>& class_);
void init_bunch_params_struct(
    py::module& m,
    py::class_<BunchParamsProxy>& class_);
void init_bunch_struct(py::module& m, py::class_<BunchProxy>& class_);
void init_bunch_track_struct(
    py::module& m,
    py::class_<BunchTrackProxy>& class_);
void init_bicubic_cmplx_coef_struct(
    py::module& m,
    py::class_<BicubicCmplxCoefProxy>& class_);
