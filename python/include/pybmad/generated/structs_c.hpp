#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_cartesian_map_struct(nb::module_ &m, nb::class_<CartesianMapStruct> &class_);
void init_cartesian_map_term1_struct(nb::module_ &m, nb::class_<CartesianMapTerm1Struct> &class_);
void init_cartesian_map_term_struct(nb::module_ &m, nb::class_<CartesianMapTermStruct> &class_);
void init_complex_taylor_struct(nb::module_ &m, nb::class_<ComplexTaylorStruct> &class_);
void init_complex_taylor_term_struct(nb::module_ &m, nb::class_<ComplexTaylorTermStruct> &class_);
void init_control_ramp1_struct(nb::module_ &m, nb::class_<ControlRamp1Struct> &class_);
void init_control_struct(nb::module_ &m, nb::class_<ControlStruct> &class_);
void init_control_var1_struct(nb::module_ &m, nb::class_<ControlVar1Struct> &class_);
void init_controller_struct(nb::module_ &m, nb::class_<ControllerStruct> &class_);
void init_coord_array_struct(nb::module_ &m, nb::class_<CoordArrayStruct> &class_);
void init_coord_struct(nb::module_ &m, nb::class_<CoordStruct> &class_);
void init_cylindrical_map_struct(nb::module_ &m, nb::class_<CylindricalMapStruct> &class_);
void init_cylindrical_map_term1_struct(
    nb::module_ &m,
    nb::class_<CylindricalMapTerm1Struct> &class_
);
void init_cylindrical_map_term_struct(nb::module_ &m, nb::class_<CylindricalMapTermStruct> &class_);
