#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_f(py::module& m);

struct PyFactorial {
  int n;
  double fact;
};
struct PyFileDirectorizer {
  std::string in_file;
  std::string out_file;
  std::string directory;
  bool add_switch;
};
struct PyFileGet {
  std::string string;
  std::string dflt_file_name;
  std::string file_name;
};
struct PyFileGetOpen {
  std::string string;
  std::string dflt_file_name;
  std::string file_name;
  int file_unit;
  bool readonly;
};
struct PyFileSuffixer {
  std::string in_file_name;
  std::string out_file_name;
  std::string suffix;
  bool add_switch;
};
struct PyFindLocationInt {
  int value;
  int ix_match;
};
struct PyFindLocationLogic {
  bool value;
  int ix_match;
};
struct PyFindLocationReal {
  int ix_match;
};
struct PyFixedwindowls {
  double z;
};
