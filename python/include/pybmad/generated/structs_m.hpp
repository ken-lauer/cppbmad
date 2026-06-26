#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_mode3_struct(nb::module_ &m, nb::class_<Mode3Struct> &class_);
void init_mode_info_struct(nb::module_ &m, nb::class_<ModeInfoStruct> &class_);
void init_mad_energy_struct(nb::module_ &m, nb::class_<MadEnergyStruct> &class_);
void init_mad_map_struct(nb::module_ &m, nb::class_<MadMapStruct> &class_);
