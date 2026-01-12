#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_d(py::module& m);

struct PyDateAndTimeStamp {
  std::string string;
  std::optional<bool> numeric_month;
  std::optional<bool> include_zone;
};
struct PyDetab {
  std::string str;
};
struct PyDisplaySizeAndResolution {
  int ix_screen;
  double x_size;
  double y_size;
  double x_res;
  double y_res;
};
struct PyDjBessel {
  int m;
  double arg;
  double dj_bes;
};
struct PyDjbHash {
  std::string str;
  std::optional<int> old_hash;
  int hash;
};
struct PyDjbStrHash {
  std::string in_str;
  std::string hash_str;
};
struct PyDowncaseString {
  std::string string;
};
