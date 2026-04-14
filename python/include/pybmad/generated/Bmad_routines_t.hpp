#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_t(py::module &m);

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
