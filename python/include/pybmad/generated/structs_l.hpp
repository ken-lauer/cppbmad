#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_lat_ele_loc_struct(nb::module_ &m, nb::class_<LatEleLocStruct> &class_);
void init_lat_ele_order1_struct(nb::module_ &m, nb::class_<LatEleOrder1Struct> &class_);
void init_lat_ele_order_array_struct(nb::module_ &m, nb::class_<LatEleOrderArrayStruct> &class_);
void init_lat_ele_order_struct(nb::module_ &m, nb::class_<LatEleOrderStruct> &class_);
void init_lat_param_struct(nb::module_ &m, nb::class_<LatParamStruct> &class_);
void init_lat_pointer_struct(nb::module_ &m, nb::class_<LatPointerStruct> &class_);
void init_lat_struct(nb::module_ &m, nb::class_<LatStruct> &class_);
void init_linac_normal_mode_struct(nb::module_ &m, nb::class_<LinacNormalModeStruct> &class_);
void init_linear_ele_isf_struct(nb::module_ &m, nb::class_<LinearEleIsfStruct> &class_);
void init_linear_isf1_struct(nb::module_ &m, nb::class_<LinearIsf1Struct> &class_);
void init_layout(nb::module_ &m, nb::class_<Layout> &class_);
