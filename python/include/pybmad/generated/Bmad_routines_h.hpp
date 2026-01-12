#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_h(py::module& m);

struct PyHasAttribute {
  std::string attrib;
  bool has_it;
};
struct PyHdf5WriteBeam {
  std::string file_name;
  bool append;
  bool error;
  std::optional<bool> alive_only;
};
struct PyHdf5WriteGridField {
  std::string file_name;
  bool err_flag;
};
