#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_g(py::module& m);

struct PyGammaRef {
  double gamma;
};
struct PyGenGradField {
  double rho;
  double theta;
};
struct PyGetCalledFile {
  std::string delim;
  std::string call_file;
  bool err;
};
struct PyGptFieldGridScaling {
  int dimensions;
  double field_scale;
  double ref_time;
};
struct PyGptMaxFieldReference {
  double field_value;
};
struct PyGradientShiftSrWake {
  double grad_shift;
};
