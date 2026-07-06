#pragma once

#include <bmad/tao_hooks.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;

namespace Pybmad {

// Binds the `tao.hooks` registry (a TaoHooks instance whose properties install
// or clear each Tao hook) onto the given (tao) module.
void init_tao_hooks(nb::module_ &m);

}; // namespace Pybmad
