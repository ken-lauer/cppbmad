#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_wake_lr_mode_struct(py::module &m, py::class_<WakeLrModeStruct> &class_);
void init_wake_lr_struct(py::module &m, py::class_<WakeLrStruct> &class_);
void init_wake_sr_mode_struct(py::module &m, py::class_<WakeSrModeStruct> &class_);
void init_wake_sr_struct(py::module &m, py::class_<WakeSrStruct> &class_);
void init_wake_sr_z_long_struct(py::module &m, py::class_<WakeSrZLongStruct> &class_);
void init_wake_struct(py::module &m, py::class_<WakeStruct> &class_);
void init_wall3d_section_struct(py::module &m, py::class_<Wall3dSectionStruct> &class_);
void init_wall3d_struct(py::module &m, py::class_<Wall3dStruct> &class_);
void init_wall3d_vertex_struct(py::module &m, py::class_<Wall3dVertexStruct> &class_);
