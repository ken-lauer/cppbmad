#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_t(py::module& m);

struct PyToStr {
  double num;
  std::optional<int> max_signif;
  std::string string;
};
struct PyTypeThisFile {
  std::string filename;
};
