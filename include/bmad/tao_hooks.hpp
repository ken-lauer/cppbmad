#pragma once

// Hand-written C++ interface to Tao's hook procedure pointers.
//
// Each Tao hook `tao_hook_<x>` can be given a C++ callable via
// `Tao::set_<x>_hook(...)`; `Tao::clear_<x>_hook()` restores Tao's default
// behavior. The callable receives non-owning proxy views of the Fortran
// derived-type arguments, valid only for the duration of the call. Callables
// must not throw across the boundary; exceptions are caught and swallowed at the
// trampoline (see src/tao_hooks/tao_hooks.cpp). See plans/hook-plan.md.
//
// This layer has no nanobind dependency and is part of libcppbmad. The Python
// bindings (python/src/tao_hooks.cpp) build on top of it.

#include <functional>

#include "bmad/generated/proxy.hpp"

namespace Tao {

// Receives the current value; returns the new `calc_ok` / `abort` value.
using LatticeCalcHook = std::function<bool(bool)>;
void set_lattice_calc_hook(LatticeCalcHook fn);
void clear_lattice_calc_hook();

using OptimizerHook = std::function<bool(bool)>;
void set_optimizer_hook(OptimizerHook fn);
void clear_optimizer_hook();

using MeritVarHook = std::function<void(int, int, Bmad::TaoVarStruct &)>;
void set_merit_var_hook(MeritVarHook fn);
void clear_merit_var_hook();

// Receives (i_uni, j_data, datum, current valid_value_set); returns new value.
using MeritDataHook = std::function<bool(int, int, Bmad::TaoDataStruct &, bool)>;
void set_merit_data_hook(MeritDataHook fn);
void clear_merit_data_hook();

} // namespace Tao
