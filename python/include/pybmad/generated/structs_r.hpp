#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_rad_int1_struct(py::module &m, py::class_<RadInt1Struct> &class_);
void init_rad_int_all_ele_struct(py::module &m, py::class_<RadIntAllEleStruct> &class_);
void init_rad_int_branch_struct(py::module &m, py::class_<RadIntBranchStruct> &class_);
void init_rad_map_ele_struct(py::module &m, py::class_<RadMapEleStruct> &class_);
void init_rad_map_struct(py::module &m, py::class_<RadMapStruct> &class_);
void init_ramper_lord_struct(py::module &m, py::class_<RamperLordStruct> &class_);
void init_resonance_h_struct(py::module &m, py::class_<ResonanceHStruct> &class_);
void init_rf_ele_struct(py::module &m, py::class_<RfEleStruct> &class_);
void init_rf_stair_step_struct(py::module &m, py::class_<RfStairStepStruct> &class_);
void init_random_state_struct(py::module &m, py::class_<RandomStateStruct> &class_);
