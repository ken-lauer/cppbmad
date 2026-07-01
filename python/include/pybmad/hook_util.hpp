#pragma once

// Shared nanobind helpers for the hand-written hook bindings
// (python/src/tao_hooks.cpp, python/src/bmad_hooks.cpp).

#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace Pybmad {
namespace hooks {

// Wrap a live proxy as a non-owning Python object (valid only for the call).
template <class T>
inline nb::object ref(T &x) {
  return nb::cast(&x, nb::rv_policy::reference);
}

// Python truthiness of a hook return value (None -> false).
inline bool obj_to_bool(const nb::object &r) { return !r.is_none() && nb::cast<bool>(r); }

// Report a Python exception so it never unwinds into the C trampoline / Fortran.
inline void report_error(nb::python_error &e) {
  e.restore();
  PyErr_Print();
}

} // namespace hooks
} // namespace Pybmad
