#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_bsim_routines_b(py::module& m);

struct PyBbuHomVoltageCalc {
  int n_period;
  int ix_stage_last_tracked;
};
struct PyBbuSetup {
  double dt_bunch;
};
struct PyBbuTrackAStage {
  bool lost;
  int ix_stage_tracked;
};
struct PyBbuTrackAll {
  double hom_voltage_normalized;
  double growth_rate;
  bool lost;
  int irep;
};
