#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_p(py::module& m);

struct PyParseFortranFormat {
  std::string format_str;
  int n_repeat;
  int power;
  std::string descrip;
  int width;
  int digits;
};
struct PyPolyEval {
  double y;
};
struct PyProbabilityFunct {
  double prob;
};
struct PyProjdd {
  std::complex<double> func_retval__;
};
