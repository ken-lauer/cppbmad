#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_SimUtils_routines_s(py::module& m);

struct PySetParameterInt {
  int param_val;
  int set_val;
  int save_val;
};
struct PySetParameterLogic {
  bool param_val;
  bool set_val;
  bool save_val;
};
struct PySetParameterReal {
  double param_val;
  double set_val;
  double save_val;
};
struct PySinc {
  double y;
};
struct PySincc {
  double y;
};
struct PySinhxX {
  double y;
};
struct PySkipHeader {
  int ix_unit;
  bool error_flag;
};
struct PySqrtAlpha {
  double y;
};
struct PySqrtOne {
  double ds1;
};
struct PyStrCount {
  std::string str;
  std::string match;
  int num;
};
struct PyStrFirstInSet {
  std::string line;
  std::string set;
  std::optional<bool> ignore_clauses;
  int ix_match;
};
struct PyStrFirstNotInSet {
  std::string line;
  std::string set;
  int ix_match;
};
struct PyStrLastInSet {
  std::string line;
  std::string set;
  int ix_match;
};
struct PyStrLastNotInSet {
  std::string line;
  std::string set;
  int ix_match;
};
struct PyStrMatchWild {
  std::string str;
  std::string pat;
  bool a_match;
};
struct PyStrSubstitute {
  std::string string;
  std::optional<std::string> str_match;
  std::optional<std::string> str_replace;
  std::optional<bool> do_trim;
  std::optional<bool> ignore_escaped;
};
struct PyStringToInt {
  std::string line;
  int default_;
  bool err_flag;
  std::optional<bool> err_print_flag;
  int value;
};
struct PyStringToReal {
  std::string line;
  double default_;
  bool err_flag;
  std::optional<bool> err_print_flag;
  double value;
};
struct PyStringTrim {
  std::string in_string;
  std::string out_string;
  int word_len;
};
struct PyStringTrim2 {
  std::string in_str;
  std::string delimitors;
  std::string out_str;
  int ix_word;
  std::string delim;
  int ix_next;
};
struct PySystemCommand {
  std::string line;
  std::optional<bool> err_flag;
};
