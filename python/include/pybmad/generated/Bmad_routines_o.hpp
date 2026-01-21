#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_o(py::module &m);

struct PyOdeintBmadTime : public Bmad::OdeintBmadTime {
  double rf_time;
  PyOdeintBmadTime(Bmad::OdeintBmadTime _base, double rf_time)
      : Bmad::OdeintBmadTime(std::move(_base))
      , rf_time(rf_time) {}
};
struct PyOffsetParticle : public Bmad::OffsetParticle {
  std::optional<double> time;
  PyOffsetParticle(Bmad::OffsetParticle _base, std::optional<double> time)
      : Bmad::OffsetParticle(std::move(_base))
      , time(time) {}
};
