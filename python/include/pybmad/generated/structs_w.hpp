#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_wake_lr_mode_struct(nb::module_ &m, nb::class_<WakeLrModeStruct> &class_);
void init_wake_lr_struct(nb::module_ &m, nb::class_<WakeLrStruct> &class_);
void init_wake_sr_mode_struct(nb::module_ &m, nb::class_<WakeSrModeStruct> &class_);
void init_wake_sr_struct(nb::module_ &m, nb::class_<WakeSrStruct> &class_);
void init_wake_sr_z_long_struct(nb::module_ &m, nb::class_<WakeSrZLongStruct> &class_);
void init_wake_struct(nb::module_ &m, nb::class_<WakeStruct> &class_);
void init_wall3d_section_struct(nb::module_ &m, nb::class_<Wall3dSectionStruct> &class_);
void init_wall3d_struct(nb::module_ &m, nb::class_<Wall3dStruct> &class_);
void init_wall3d_vertex_struct(nb::module_ &m, nb::class_<Wall3dVertexStruct> &class_);
