// Hand-written C++ glue between Tao's Fortran hook trampolines
// (src/hooks/tao_hooks.f90) and C++ callables declared in
// include/bmad/tao_hooks.hpp.
//
// For each hook we keep a file-local std::function slot and an extern "C"
// trampoline matching the Fortran-side C ABI. The trampoline wraps each void*
// argument in a non-owning proxy, invokes the std::function inside a
// try/catch(...) so no exception unwinds into Fortran, and writes back any
// out-parameters. Registration installs/removes the trampoline via the
// `tao_hook_register_*` bind(c) entry points.

#include "bmad/tao_hooks.hpp"

#include <utility>

#include "bmad/hook_util.hpp"

using Bmad::make_fortran_proxy;
using Bmad::hooks::log_hook_error;

// --- Fortran registration entry points -------------------------------------
extern "C" {
void tao_hook_register_lattice_calc(void (*)(int *));
void tao_hook_register_optimizer(void (*)(int *));
void tao_hook_register_merit_var(void (*)(int, int, void *));
void tao_hook_register_merit_data(void (*)(int, int, void *, int *));
}

namespace {

Tao::LatticeCalcHook g_lattice_calc;
extern "C" void tramp_lattice_calc(int *calc_ok) {
  if (!g_lattice_calc)
    return;
  try {
    *calc_ok = g_lattice_calc(*calc_ok != 0) ? 1 : 0;
  } catch (...) {
    log_hook_error("lattice_calc");
  }
}

Tao::OptimizerHook g_optimizer;
extern "C" void tramp_optimizer(int *abort) {
  if (!g_optimizer)
    return;
  try {
    *abort = g_optimizer(*abort != 0) ? 1 : 0;
  } catch (...) {
    log_hook_error("optimizer");
  }
}

Tao::MeritVarHook g_merit_var;
extern "C" void tramp_merit_var(int i_uni, int j_var, void *var) {
  if (!g_merit_var)
    return;
  try {
    auto pv = make_fortran_proxy<Bmad::TaoVarStruct>(var);
    g_merit_var(i_uni, j_var, pv);
  } catch (...) {
    log_hook_error("merit_var");
  }
}

Tao::MeritDataHook g_merit_data;
extern "C" void tramp_merit_data(int i_uni, int j_data, void *data, int *valid_value_set) {
  if (!g_merit_data)
    return;
  try {
    auto pd = make_fortran_proxy<Bmad::TaoDataStruct>(data);
    *valid_value_set = g_merit_data(i_uni, j_data, pd, *valid_value_set != 0) ? 1 : 0;
  } catch (...) {
    log_hook_error("merit_data");
  }
}

} // namespace

// --- set_*/clear_* -----------------------------------------------------------

#define TAO_HOOK_SETCLEAR(name, HookType, slot, regfn, tramp) \
  void Tao::set_##name##_hook(HookType fn) {                  \
    slot = std::move(fn);                                     \
    regfn(slot ? &tramp : nullptr);                           \
  }                                                           \
  void Tao::clear_##name##_hook() {                           \
    slot = nullptr;                                           \
    regfn(nullptr);                                           \
  }

TAO_HOOK_SETCLEAR(
    lattice_calc,
    LatticeCalcHook,
    g_lattice_calc,
    tao_hook_register_lattice_calc,
    tramp_lattice_calc
)
TAO_HOOK_SETCLEAR(
    optimizer,
    OptimizerHook,
    g_optimizer,
    tao_hook_register_optimizer,
    tramp_optimizer
)
TAO_HOOK_SETCLEAR(
    merit_var,
    MeritVarHook,
    g_merit_var,
    tao_hook_register_merit_var,
    tramp_merit_var
)
TAO_HOOK_SETCLEAR(
    merit_data,
    MeritDataHook,
    g_merit_data,
    tao_hook_register_merit_data,
    tramp_merit_data
)

#undef TAO_HOOK_SETCLEAR
