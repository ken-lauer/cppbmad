#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_n(py::module& m);

struct PyNormalFormComplexTaylors {
  bool rf_on;
  std::optional<int> order;
};
struct PyNumFieldEles {
  int n_field_ele;
};
struct PyNumLords {
  int num;
};
