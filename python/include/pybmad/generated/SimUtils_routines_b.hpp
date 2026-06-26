#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_SimUtils_routines_b(nb::module_ &m);

struct PyBinXCenter {
  double x_center;
  int ix_bin;
};
struct PyBitSet {
  int word;
};
