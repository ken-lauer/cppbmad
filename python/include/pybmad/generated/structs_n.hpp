#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_normal_modes_struct(py::module &m, py::class_<NormalModesStruct> &class_);
void init_nametable_struct(py::module &m, py::class_<NametableStruct> &class_);
