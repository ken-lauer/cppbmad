#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_m(py::module& m);

struct PyMakeLegalComment {
  std::string comment_in;
  std::string comment_out;
};
struct PyMatchReg {
  std::string str;
  std::string pat;
  bool is_match;
};
struct PyMatchWild {
  std::string string;
  std::string template_;
  bool is_match;
};
struct PyMaximizeProjection {
  double seed;
  double func_retval__;
};
struct PyMilliSleep {
  int milli_sec;
};
