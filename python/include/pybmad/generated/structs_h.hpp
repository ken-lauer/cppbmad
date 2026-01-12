#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_high_energy_space_charge_struct(
    py::module& m,
    py::class_<HighEnergySpaceChargeProxy>& class_);
