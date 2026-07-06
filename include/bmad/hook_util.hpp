#pragma once

// Shared helpers for the hand-written hook glue (src/hooks/*.cpp).

#include <cstddef>
#include <cstdio>

#include "bmad/generated/proxy.hpp"

namespace Bmad {
namespace hooks {

inline void log_hook_error(const char *hook) {
  std::fprintf(stderr, "[cppbmad] exception in hook '%s' was swallowed\n", hook);
}

// Build a non-owning live view of a Fortran coord_struct array from its base
// address, inclusive bounds, and element size. Returns an empty (invalid) view
// when the array is absent or empty. Mutations through the view's element
// proxies write straight back into the Fortran storage.
inline CoordStructArray1D make_coord_view(void *data, int lb, int ub, std::size_t esize) {
  if (!data || ub < lb)
    return CoordStructArray1D();
  int n = ub - lb + 1;
  return CoordStructArray1D(data, n, lb, ub, true, esize);
}

} // namespace hooks
} // namespace Bmad
