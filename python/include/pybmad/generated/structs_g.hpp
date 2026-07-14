#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_general_bin_struct(nb::module_ &m, nb::class_<GeneralBinStruct> &class_);
void init_gen_grad1_struct(nb::module_ &m, nb::class_<GenGrad1Struct> &class_);
void init_gen_grad_map_struct(nb::module_ &m, nb::class_<GenGradMapStruct> &class_);
void init_grid_beam_init_struct(nb::module_ &m, nb::class_<GridBeamInitStruct> &class_);
void init_grid_field_pt1_struct(nb::module_ &m, nb::class_<GridFieldPt1Struct> &class_);
void init_grid_field_pt_struct(nb::module_ &m, nb::class_<GridFieldPtStruct> &class_);
void init_grid_field_struct(nb::module_ &m, nb::class_<GridFieldStruct> &class_);
void init_gpt_lat_param_struct(nb::module_ &m, nb::class_<GptLatParamStruct> &class_);
