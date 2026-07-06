#pragma once

// Shared nanobind helpers for the hand-written hook bindings
// (python/src/tao_hooks.cpp, python/src/bmad_hooks.cpp).

#include <nanobind/nanobind.h>

#include <cstdio>
#include <exception>

namespace nb = nanobind;

namespace Pybmad {
namespace hooks {

// Wrap a live proxy as a non-owning Python object (valid only for the call).
template <class T>
inline nb::object ref(T &x) {
  return nb::cast(&x, nb::rv_policy::reference);
}

// Report a Python exception raised inside a hook callback. Uses the
// unraisable-hook path (not PyErr_Print) so that a SystemExit / KeyboardInterrupt
// does not terminate the process, and the error never unwinds into the C
// trampoline / Fortran.
inline void report_error(nb::python_error &e) {
  e.restore();
  PyErr_WriteUnraisable(nullptr);
}

// Report a C++ exception from a hook callback -- most commonly a return value of
// the wrong shape (a nanobind cast_error). Names the hook and the reason so the
// user sees more than the generic trampoline "exception was swallowed".
inline void report_error(const char *hook, const std::exception &e) {
  std::fprintf(stderr, "[cppbmad] hook '%s' returned an unexpected value: %s\n", hook, e.what());
}

} // namespace hooks
} // namespace Pybmad
