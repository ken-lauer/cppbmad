#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_Bmad_routines_s(nb::module_ &m);

struct PyScAdaptiveStep {
  double dt_next;
  bool include_image;
  double dt_step;
};

struct PyScStep {
  int n_emit;
  bool include_image;
};
struct PySetFringeOnOff {
  double fringe_at;
};
struct PySetPtcQuiet {
  int old_val;
};
