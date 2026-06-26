#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_Bmad_routines_r(nb::module_ &m);

struct PyRadiationIntegrals : public Bmad::RadiationIntegrals {
  std::optional<int> ix_cache;
  PyRadiationIntegrals(Bmad::RadiationIntegrals _base, std::optional<int> ix_cache)
      : Bmad::RadiationIntegrals(std::move(_base))
      , ix_cache(ix_cache) {}
};
struct PyReleaseRadIntCache {
  int ix_cache;
};
