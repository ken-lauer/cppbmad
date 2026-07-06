#pragma once

// Hand-written C++ interface to Bmad's tracking hook / custom procedure
// pointers (bmad_routine_interface).
//
// Each hook `<x>` can be given a C++ callable via `Bmad::set_<x>_hook(...)`;
// `Bmad::clear_<x>_hook()` restores Bmad's default behavior. The callable
// receives non-owning proxy views of the Fortran derived-type arguments, valid
// only for the duration of the call, and out-parameters by reference (write into
// them to affect Bmad). Optional derived-type / integer / logical arguments are
// passed as pointers that are null when the argument is absent; optional and
// required `coord_struct` arrays are passed as live `CoordStructArray1D` views.
//
// Callables must not throw across the boundary; exceptions are caught and
// swallowed at the trampoline (see src/hooks/bmad_hooks.cpp). This layer has no
// nanobind dependency. See plans/hook-plan.md.
//
// The `std::function` typedefs (e.g. Bmad::Track1BunchHook) and the
// set_<x>_hook / clear_<x>_hook declarations are generated from the [[hooks]]
// spec (codegen/hooks.py) into bmad/generated/hook_abi.hpp, together with the
// matching Fortran ci_<name> interfaces, so the C ABI cannot drift.

#include "bmad/generated/hook_abi.hpp"
