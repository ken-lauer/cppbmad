#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_g(py::module& m);

struct PyGenCompleteElliptic {
  double kc;
  double p;
  double c;
  double s;
  std::optional<double> err_tol;
  double value;
};
struct PyGetFileNumber {
  std::string file_name;
  std::string cnum_in;
  int num_out;
  bool err_flag;
};
struct PyGetFileTimeStamp {
  std::string file;
  std::string time_stamp;
};
