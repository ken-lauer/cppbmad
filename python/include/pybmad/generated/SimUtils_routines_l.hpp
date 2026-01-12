#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_l(py::module& m);

struct PyLinearFit {
  int n_data;
  double a;
  double b;
  double sig_a;
  double sig_b;
};
struct PyLogicStr {
  bool logic;
  std::string str;
};
