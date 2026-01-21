#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_r(py::module &m);

struct PyRadiationIntegrals : public Bmad::RadiationIntegrals {
  std::optional<int> ix_cache;
  PyRadiationIntegrals(Bmad::RadiationIntegrals _base, std::optional<int> ix_cache)
      : Bmad::RadiationIntegrals(std::move(_base))
      , ix_cache(ix_cache) {}
};
struct PyReleaseRadIntCache {
  int ix_cache;
};
