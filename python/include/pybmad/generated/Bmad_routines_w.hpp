#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_w(py::module& m);

struct PyWordToValue {
  std::string word;
  double value;
  bool err_flag;
};
struct PyWriteAstraBend {
  int iu;
  double strength;
  int id;
};
struct PyWriteBlenderEle {
  int iu;
  std::optional<bool> old_format;
};
struct PyWriteBlenderLatLayout {
  std::string file_name;
};
struct PyWriteLatLine {
  std::string line;
};
struct PyWriteLatticeInSadFormat {
  std::string out_file_name;
  std::optional<bool> include_apertures;
  std::optional<int> ix_branch;
  std::optional<bool> err;
};
struct PyWriteLineElement {
  std::string line;
  int iu;
};
