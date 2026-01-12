#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_n(py::module& m);

struct PyNBinsAutomatic {
  int n_data;
  int n;
};
struct PyNChooseK {
  int n;
  int k;
  double nck;
};
struct PyNaff {
  std::optional<int> opt_dump_spectra;
  std::optional<bool> opt_zero_first;
};
struct PyNametableAdd {
  std::string name;
  int ix_name;
};
struct PyNametableBracketIndexx {
  std::string name;
  std::optional<int> n_match;
  int ix_max;
};
struct PyNametableChange1 {
  std::string name;
  int ix_name;
};
struct PyNametableInit {
  std::optional<int> n_min;
  std::optional<int> n_max;
};
struct PyNametableRemove {
  int ix_name;
};
