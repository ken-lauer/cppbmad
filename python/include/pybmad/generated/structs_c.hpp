#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_cartesian_map_struct(
    py::module& m,
    py::class_<CartesianMapProxy>& class_);
void init_cartesian_map_term1_struct(
    py::module& m,
    py::class_<CartesianMapTerm1Proxy>& class_);
void init_cartesian_map_term_struct(
    py::module& m,
    py::class_<CartesianMapTermProxy>& class_);
void init_complex_taylor_struct(
    py::module& m,
    py::class_<ComplexTaylorProxy>& class_);
void init_complex_taylor_term_struct(
    py::module& m,
    py::class_<ComplexTaylorTermProxy>& class_);
void init_control_ramp1_struct(
    py::module& m,
    py::class_<ControlRamp1Proxy>& class_);
void init_control_struct(py::module& m, py::class_<ControlProxy>& class_);
void init_control_var1_struct(
    py::module& m,
    py::class_<ControlVar1Proxy>& class_);
void init_controller_struct(py::module& m, py::class_<ControllerProxy>& class_);
void init_coord_array_struct(
    py::module& m,
    py::class_<CoordArrayProxy>& class_);
void init_coord_struct(py::module& m, py::class_<CoordProxy>& class_);
void init_cylindrical_map_struct(
    py::module& m,
    py::class_<CylindricalMapProxy>& class_);
void init_cylindrical_map_term1_struct(
    py::module& m,
    py::class_<CylindricalMapTerm1Proxy>& class_);
void init_cylindrical_map_term_struct(
    py::module& m,
    py::class_<CylindricalMapTermProxy>& class_);
