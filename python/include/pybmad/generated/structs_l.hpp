#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_lat_ele_loc_struct(py::module &m, py::class_<LatEleLocStruct> &class_);
void init_lat_ele_order1_struct(py::module &m, py::class_<LatEleOrder1Struct> &class_);
void init_lat_ele_order_array_struct(py::module &m, py::class_<LatEleOrderArrayStruct> &class_);
void init_lat_ele_order_struct(py::module &m, py::class_<LatEleOrderStruct> &class_);
void init_lat_param_struct(py::module &m, py::class_<LatParamStruct> &class_);
void init_lat_struct(py::module &m, py::class_<LatStruct> &class_);
void init_linac_normal_mode_struct(py::module &m, py::class_<LinacNormalModeStruct> &class_);
void init_layout(py::module &m, py::class_<Layout> &class_);
