#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_floor_position_struct(py::module &m, py::class_<FloorPositionStruct> &class_);
void init_fibre(py::module &m, py::class_<Fibre> &class_);
