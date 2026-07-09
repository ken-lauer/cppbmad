#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_general_bin_struct(nb::module_ &m, nb::class_<GeneralBinStruct> &class_);
void init_gen_grad_curve_struct(nb::module_ &m, nb::class_<GenGradCurveStruct> &class_);
void init_gen_gradients_struct(nb::module_ &m, nb::class_<GenGradientsStruct> &class_);
void init_gg_taylor_struct(nb::module_ &m, nb::class_<GgTaylorStruct> &class_);
void init_gg_taylor_term_struct(nb::module_ &m, nb::class_<GgTaylorTermStruct> &class_);
void init_grid_beam_init_struct(nb::module_ &m, nb::class_<GridBeamInitStruct> &class_);
void init_grid_field_pt1_struct(nb::module_ &m, nb::class_<GridFieldPt1Struct> &class_);
void init_grid_field_pt_struct(nb::module_ &m, nb::class_<GridFieldPtStruct> &class_);
void init_grid_field_struct(nb::module_ &m, nb::class_<GridFieldStruct> &class_);
void init_gpt_lat_param_struct(nb::module_ &m, nb::class_<GptLatParamStruct> &class_);
