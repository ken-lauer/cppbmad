#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_r(py::module& m);

struct PyRanGaussScalar {
  double harvest;
  std::optional<int> index_quasi;
};

struct PyRanUniformScalar {
  double harvest;
  std::optional<int> index_quasi;
};
struct PyRealNumFortranFormat {
  double number;
  int width;
  std::optional<int> n_blanks;
  std::string fmt_str;
};
struct PyRealPath {
  std::string path_in;
  std::string path_out;
  bool is_ok;
};
struct PyRealStr {
  double r_num;
  std::optional<int> n_signif;
  std::optional<int> n_decimal;
  std::string str;
};
struct PyRealToString {
  double real_num;
  int width;
  std::optional<int> n_signif;
  std::optional<int> n_decimal;
  std::string str;
};

struct PyRmsValue {
  double ave_val;
  double rms_val;
};
struct PyRunTimer {
  std::string command;
  std::optional<double> time;
  std::optional<double> time0;
};
