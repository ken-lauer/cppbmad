#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_structs_f(py::module& m);

// Per-struct init functions
void init_floor_position_struct(
    py::module& m,
    py::class_<FloorPositionProxy>& class_);
