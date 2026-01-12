#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_c(py::module& m);

struct PyCalcFileNumber {
  std::string file_name;
  int num_in;
  int num_out;
  bool err_flag;
};
struct PyChangeFileNumber {
  std::string file_name;
  int change;
};

struct PyCoarseFrequencyEstimate {
  double frequency;
  std::optional<bool> error;
};
struct PyComplexErrorFunction {
  double wr;
  double wi;
  double zr;
  double zi;
};
struct PyCosOne {
  double cos1;
};
struct PyCosc {
  double y;
};
