#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_q(py::module& m);

struct PyQueryString {
  std::string query_str;
  bool upcase;
  std::string return_str;
  int ix;
  int ios;
};
struct PyQuote {
  std::string str;
  std::string q_str;
};
