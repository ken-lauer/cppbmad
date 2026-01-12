#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_structs_x(py::module& m);

// Per-struct init functions
void init_xy_disp_struct(py::module& m, py::class_<XyDispProxy>& class_);
