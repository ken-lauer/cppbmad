#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_rad_int1_struct(py::module& m, py::class_<RadInt1Proxy>& class_);
void init_rad_int_all_ele_struct(
    py::module& m,
    py::class_<RadIntAllEleProxy>& class_);
void init_rad_int_branch_struct(
    py::module& m,
    py::class_<RadIntBranchProxy>& class_);
void init_rad_map_ele_struct(py::module& m, py::class_<RadMapEleProxy>& class_);
void init_rad_map_struct(py::module& m, py::class_<RadMapProxy>& class_);
void init_ramper_lord_struct(
    py::module& m,
    py::class_<RamperLordProxy>& class_);
void init_resonance_h_struct(
    py::module& m,
    py::class_<ResonanceHProxy>& class_);
void init_rf_ele_struct(py::module& m, py::class_<RfEleProxy>& class_);
void init_rf_stair_step_struct(
    py::module& m,
    py::class_<RfStairStepProxy>& class_);
void init_random_state_struct(
    py::module& m,
    py::class_<RandomStateProxy>& class_);
