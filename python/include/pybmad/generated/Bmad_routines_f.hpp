#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_f(py::module& m);

struct PyFibreToEle {
  bool err_flag;
  int ix_ele;
};
struct PyFringeHere {
  bool is_here;
};
