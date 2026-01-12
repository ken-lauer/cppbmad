#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_v(py::module& m);

struct PyValidFieldCalc {
  bool is_valid;
};
struct PyValidFringeType {
  bool is_valid;
};
struct PyValidMat6CalcMethod {
  bool is_valid;
};
struct PyValidSpinTrackingMethod {
  bool is_valid;
};
struct PyValidTrackingMethod {
  bool is_valid;
};

struct PyValueOfAttribute {
  bool err_flag;
  double value;
};
struct PyValueToLine {
  std::string line;
  double value;
  std::string str;
  std::string typ;
  std::optional<bool> ignore_if_zero;
  std::optional<bool> use_comma;
};
