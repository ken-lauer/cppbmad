#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_a(py::module& m);

struct PyApfft {
  std::string window;
  double phase;
  std::optional<int> diag;
};
struct PyApfftExt {
  std::string window;
  double phase;
  double amp;
  double freq;
  std::optional<int> diag;
};
struct PyAsinc {
  double y;
};
struct PyAssertEqual {
  std::string err_str;
  int ival;
};
