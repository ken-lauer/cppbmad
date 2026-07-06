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
//
// The `std::function` typedefs (e.g. Tao::MeritDataHook) and the set_<x>_hook /
// clear_<x>_hook declarations are generated from the [[hooks]] spec
// (codegen/hooks.py) into bmad/generated/hook_abi.hpp, together with the
// matching Fortran ci_<name> interfaces, so the C ABI cannot drift. The Tao
// hooks receive the current value and return the new one (calc_ok / abort /
// valid_value_set).

#include "bmad/generated/hook_abi.hpp"
