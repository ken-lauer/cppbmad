#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_structs_i(py::module& m);

// Per-struct init functions
void init_interval1_coef_struct(
    py::module& m,
    py::class_<Interval1CoefProxy>& class_);
