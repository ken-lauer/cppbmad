#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_Bmad_routines_c(nb::module_ &m);

struct PyChromCalc : public Bmad::ChromCalc {
  double delta_e;
  PyChromCalc(Bmad::ChromCalc _base, double delta_e)
      : Bmad::ChromCalc(std::move(_base))
      , delta_e(delta_e) {}
};

struct PyChromTune {
  bool err_flag;
  double delta_e;
};
