#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_high_energy_space_charge_struct(
    py::module &m,
    py::class_<HighEnergySpaceChargeStruct> &class_
);
