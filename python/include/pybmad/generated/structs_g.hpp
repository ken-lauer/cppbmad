#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_gen_grad1_struct(nb::module_ &m, nb::class_<GenGrad1Struct> &class_);
void init_gen_grad_map_struct(nb::module_ &m, nb::class_<GenGradMapStruct> &class_);
void init_gg_taylor_struct(nb::module_ &m, nb::class_<GgTaylorStruct> &class_);
void init_gg_taylor_term_struct(nb::module_ &m, nb::class_<GgTaylorTermStruct> &class_);
void init_grid_beam_init_struct(nb::module_ &m, nb::class_<GridBeamInitStruct> &class_);
void init_grid_field_pt1_struct(nb::module_ &m, nb::class_<GridFieldPt1Struct> &class_);
void init_grid_field_pt_struct(nb::module_ &m, nb::class_<GridFieldPtStruct> &class_);
void init_grid_field_struct(nb::module_ &m, nb::class_<GridFieldStruct> &class_);
