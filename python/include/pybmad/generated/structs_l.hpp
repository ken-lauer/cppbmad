#pragma once
#include <pybind11/pybind11.h>
#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_lat_ele_loc_struct(py::module& m, py::class_<LatEleLocProxy>& class_);
void init_lat_ele_order1_struct(
    py::module& m,
    py::class_<LatEleOrder1Proxy>& class_);
void init_lat_ele_order_array_struct(
    py::module& m,
    py::class_<LatEleOrderArrayProxy>& class_);
void init_lat_ele_order_struct(
    py::module& m,
    py::class_<LatEleOrderProxy>& class_);
void init_lat_param_struct(py::module& m, py::class_<LatParamProxy>& class_);
void init_lat_struct(py::module& m, py::class_<LatProxy>& class_);
void init_linac_normal_mode_struct(
    py::module& m,
    py::class_<LinacNormalModeProxy>& class_);
void init_layout(py::module& m, py::class_<LayoutRawStruct>& class_);
