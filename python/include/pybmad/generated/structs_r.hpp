#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_rad_int1_struct(nb::module_ &m, nb::class_<RadInt1Struct> &class_);
void init_rad_int_all_ele_struct(nb::module_ &m, nb::class_<RadIntAllEleStruct> &class_);
void init_rad_int_branch_struct(nb::module_ &m, nb::class_<RadIntBranchStruct> &class_);
void init_rad_map_ele_struct(nb::module_ &m, nb::class_<RadMapEleStruct> &class_);
void init_rad_map_struct(nb::module_ &m, nb::class_<RadMapStruct> &class_);
void init_ramper_lord_struct(nb::module_ &m, nb::class_<RamperLordStruct> &class_);
void init_resonance_h_struct(nb::module_ &m, nb::class_<ResonanceHStruct> &class_);
void init_rf_ele_struct(nb::module_ &m, nb::class_<RfEleStruct> &class_);
void init_rf_stair_step_struct(nb::module_ &m, nb::class_<RfStairStepStruct> &class_);
void init_random_state_struct(nb::module_ &m, nb::class_<RandomStateStruct> &class_);
