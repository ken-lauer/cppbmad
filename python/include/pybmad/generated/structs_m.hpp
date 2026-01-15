#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_mode3_struct(py::module &m, py::class_<Mode3Struct> &class_);
void init_mode_info_struct(py::module &m, py::class_<ModeInfoStruct> &class_);
void init_mad_energy_struct(py::module &m, py::class_<MadEnergyStruct> &class_);
void init_mad_map_struct(py::module &m, py::class_<MadMapStruct> &class_);
