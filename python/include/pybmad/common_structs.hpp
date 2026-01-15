#pragma once

#include <bmad/common_structs.hpp>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace Pybmad {

void init_common_structs(py::module &m);
}; // namespace Pybmad
