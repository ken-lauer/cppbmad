#pragma once

#include <bmad/bmad_hooks.hpp>
#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace Pybmad {

// Registers `bmad.set_<x>_hook` / `bmad.clear_<x>_hook` on the given (bmad) module.
void init_bmad_hooks(nb::module_ &m);

}; // namespace Pybmad
