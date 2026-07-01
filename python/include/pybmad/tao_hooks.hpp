#pragma once

#include <bmad/tao_hooks.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;

namespace Pybmad {

// Registers `tao.set_<x>_hook` / `tao.clear_<x>_hook` on the given (tao) module.
void init_tao_hooks(nb::module_ &m);

}; // namespace Pybmad
