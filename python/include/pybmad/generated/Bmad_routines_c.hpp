#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_c(py::module &m);

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
