#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_wake_lr_mode_struct(
    py::module& m,
    py::class_<WakeLrModeProxy>& class_);
void init_wake_lr_struct(py::module& m, py::class_<WakeLrProxy>& class_);
void init_wake_sr_mode_struct(
    py::module& m,
    py::class_<WakeSrModeProxy>& class_);
void init_wake_sr_struct(py::module& m, py::class_<WakeSrProxy>& class_);
void init_wake_sr_z_long_struct(
    py::module& m,
    py::class_<WakeSrZLongProxy>& class_);
void init_wake_struct(py::module& m, py::class_<WakeProxy>& class_);
void init_wall3d_section_struct(
    py::module& m,
    py::class_<Wall3dSectionProxy>& class_);
void init_wall3d_struct(py::module& m, py::class_<Wall3dProxy>& class_);
void init_wall3d_vertex_struct(
    py::module& m,
    py::class_<Wall3dVertexProxy>& class_);
