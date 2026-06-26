#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_Bmad_routines_o(nb::module_ &m);

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
