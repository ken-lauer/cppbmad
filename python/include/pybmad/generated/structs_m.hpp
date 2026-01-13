#pragma once
#include <pybind11/pybind11.h>
#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_mode3_struct(py::module& m, py::class_<Mode3Proxy>& class_);
void init_mode_info_struct(py::module& m, py::class_<ModeInfoProxy>& class_);
void init_mad_energy_struct(py::module& m, py::class_<MadEnergyProxy>& class_);
void init_mad_map_struct(py::module& m, py::class_<MadMapProxy>& class_);
