#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_w(py::module& m);

struct PyWordLen {
  std::string wording;
  int wlen;
};
struct PyWordRead {
  std::string in_str;
  std::string delim_list;
  std::string word;
  int ix_word;
  std::string delim;
  bool delim_found;
  std::string out_str;
  std::optional<bool> ignore_interior;
};
