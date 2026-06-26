#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_Bmad_routines_t(nb::module_ &m);

struct PyTargetMinMaxCalc {
  double y_min;
  double y_max;
  double phi_min;
  double phi_max;
};
struct PyTrack1TimeRungeKutta : public Bmad::Track1TimeRungeKutta {
  std::optional<double> dt_step;
  PyTrack1TimeRungeKutta(Bmad::Track1TimeRungeKutta _base, std::optional<double> dt_step)
      : Bmad::Track1TimeRungeKutta(std::move(_base))
      , dt_step(dt_step) {}
};
struct PyTrackADrift {
  std::optional<double> time;
};
struct PyTypeTaylors {
  std::optional<int> n_lines;
};
