#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_l(py::module& m);

struct PyLafun {
  double x;
  double y;
  double z;
  double res;
};

struct PyLatEleLocator {
  bool err;
  int n_loc;
};
struct PyLordEdgeAligned {
  bool is_aligned;
};
struct PyLowEnergyZCorrection {
  double dz;
};
