#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_a(py::module& m);

struct PyAbsoluteTimeTracking {
  bool is_abs_time;
};
struct PyAcKickerAmp {
  double ac_amp;
};
struct PyAddThisTaylorTerm {
  int i_out;
  double coef;
};
struct PyAdjustSuperSlaveNames {
  int ix1_lord;
  int ix2_lord;
  std::optional<bool> first_time;
};
struct PyAngleBetweenPolars {
  double angle;
};
struct PyArrayReStr {
  std::optional<std::string> parens_in;
  std::string str_out;
};
struct PyAstraMaxFieldReference {
  double field_value;
};
struct PyAtThisEleEnd {
  bool is_at_this_end;
};
