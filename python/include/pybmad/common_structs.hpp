#pragma once

#include <bmad/common_structs.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;

namespace Pybmad {

void init_common_structs(nb::module_ &m);
}; // namespace Pybmad
