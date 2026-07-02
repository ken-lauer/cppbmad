#pragma once

#include <bmad/bmad_hooks.hpp>
#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace Pybmad {

// Binds the `bmad.hooks` registry (a BmadHooks instance whose properties install
// or clear each tracking hook) onto the given (bmad) module.
void init_bmad_hooks(nb::module_ &m);

}; // namespace Pybmad
