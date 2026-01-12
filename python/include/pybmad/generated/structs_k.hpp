#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_kv_beam_init_struct(
    py::module& m,
    py::class_<KvBeamInitProxy>& class_);
