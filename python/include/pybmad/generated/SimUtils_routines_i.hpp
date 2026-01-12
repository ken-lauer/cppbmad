#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_i(py::module& m);

struct PyIBessel {
  int m;
  double arg;
  double i_bes;
};
struct PyIBesselExtended {
  int m;
  double arg;
  std::complex<double> i_bes;
};
struct PyIncrementFileNumber {
  std::string file_name;
  int digits;
  int number;
  std::string cnumber;
};
struct PyIndexNocase {
  std::string string1;
  std::string string2;
  int indx;
};
struct PyIntStr {
  int int_;
  std::optional<int> width;
  std::string str;
};
struct PyInterpolatedFft {
  bool calc_ok;
  std::optional<int> opt_dump_spectrum;
  std::optional<int> opt_dump_index;
  double this_fft;
};
struct PyInterpolatedFftGsl {
  bool calc_ok;
  std::optional<int> opt_dump_spectrum;
  std::optional<int> opt_dump_index;
  double this_fft;
};
struct PyIsAlphabetic {
  std::string string;
  std::optional<std::string> valid_chars;
  bool is_alpha;
};
struct PyIsDecreasingSequence {
  bool is_decreasing;
};
struct PyIsIncreasingSequence {
  bool is_increasing;
};
struct PyIsInteger {
  std::string string;
  std::optional<int> int_;
  std::optional<std::string> delims;
  std::optional<int> ix_word;
  bool valid;
};
struct PyIsLogical {
  std::string string;
  std::optional<bool> ignore;
  bool valid;
};
struct PyIsReal {
  std::string string;
  std::optional<bool> ignore;
  std::optional<double> real_num;
  bool valid;
};
