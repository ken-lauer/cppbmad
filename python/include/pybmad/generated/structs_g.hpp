#pragma once
#include <pybind11/pybind11.h>
#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_gen_grad1_struct(py::module& m, py::class_<GenGrad1Proxy>& class_);
void init_gen_grad_map_struct(
    py::module& m,
    py::class_<GenGradMapProxy>& class_);
void init_grid_beam_init_struct(
    py::module& m,
    py::class_<GridBeamInitProxy>& class_);
void init_grid_field_pt1_struct(
    py::module& m,
    py::class_<GridFieldPt1Proxy>& class_);
void init_grid_field_pt_struct(
    py::module& m,
    py::class_<GridFieldPtProxy>& class_);
void init_grid_field_struct(py::module& m, py::class_<GridFieldProxy>& class_);
