#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_m(py::module& m);

struct PyMakeMat6Taylor {
  CoordProxy end_orb;
  std::optional<bool> err_flag;
};
struct PyMakeupControlSlave {
  bool err_flag;
};
struct PyMakeupGroupLord {
  bool err_flag;
};
struct PyMakeupMultipassSlave {
  bool err_flag;
};
struct PyMakeupSuperSlave {
  bool err_flag;
};
struct PyMasterParameterValue {
  double value;
};

struct PyMat4Multipole {
  FixedArray2D<Real, 4, 4> kick_mat;
  int n;
};
struct PyMexp {
  double this_exp;
};
struct PyMomentumCompaction {
  double mom_comp;
};
struct PyMytan {
  double y;
  double x;
  double arg;
};
