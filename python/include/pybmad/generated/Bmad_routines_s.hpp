#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_s(py::module &m);

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
