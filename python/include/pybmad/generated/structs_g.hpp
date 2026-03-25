#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_gen_grad1_struct(py::module &m, py::class_<GenGrad1Struct> &class_);
void init_gen_grad_map_struct(py::module &m, py::class_<GenGradMapStruct> &class_);
void init_gg_taylor_struct(py::module &m, py::class_<GgTaylorStruct> &class_);
void init_gg_taylor_term_struct(py::module &m, py::class_<GgTaylorTermStruct> &class_);
void init_grid_beam_init_struct(py::module &m, py::class_<GridBeamInitStruct> &class_);
void init_grid_field_pt1_struct(py::module &m, py::class_<GridFieldPt1Struct> &class_);
void init_grid_field_pt_struct(py::module &m, py::class_<GridFieldPtStruct> &class_);
void init_grid_field_struct(py::module &m, py::class_<GridFieldStruct> &class_);
