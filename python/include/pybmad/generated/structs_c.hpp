#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_cartesian_map_struct(py::module &m, py::class_<CartesianMapStruct> &class_);
void init_cartesian_map_term1_struct(py::module &m, py::class_<CartesianMapTerm1Struct> &class_);
void init_cartesian_map_term_struct(py::module &m, py::class_<CartesianMapTermStruct> &class_);
void init_complex_taylor_struct(py::module &m, py::class_<ComplexTaylorStruct> &class_);
void init_complex_taylor_term_struct(py::module &m, py::class_<ComplexTaylorTermStruct> &class_);
void init_control_ramp1_struct(py::module &m, py::class_<ControlRamp1Struct> &class_);
void init_control_struct(py::module &m, py::class_<ControlStruct> &class_);
void init_control_var1_struct(py::module &m, py::class_<ControlVar1Struct> &class_);
void init_controller_struct(py::module &m, py::class_<ControllerStruct> &class_);
void init_coord_array_struct(py::module &m, py::class_<CoordArrayStruct> &class_);
void init_coord_struct(py::module &m, py::class_<CoordStruct> &class_);
void init_cylindrical_map_struct(py::module &m, py::class_<CylindricalMapStruct> &class_);
void init_cylindrical_map_term1_struct(
    py::module &m,
    py::class_<CylindricalMapTerm1Struct> &class_
);
void init_cylindrical_map_term_struct(py::module &m, py::class_<CylindricalMapTermStruct> &class_);
