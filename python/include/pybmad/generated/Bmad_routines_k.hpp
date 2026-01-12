#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_k(py::module& m);

struct PyKeyNameToKeyIndex {
  int key_index;
};
struct PyKickVectorCalc : public Bmad::KickVectorCalc {
  std::optional<bool> print_err;
  PyKickVectorCalc(Bmad::KickVectorCalc _base, std::optional<bool> print_err)
      : Bmad::KickVectorCalc(std::move(_base)), print_err(print_err) {}
};

struct PyKnotInterpolate {
  bool err_flag;
  double y_pt;
};
struct PyKnotsToString {
  std::string str;
};
