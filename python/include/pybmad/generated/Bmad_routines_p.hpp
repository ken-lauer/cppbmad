#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_p(py::module& m);

struct PyParseCartesianMap {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
struct PyParseCylindricalMap {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
struct PyParseGenGradMap {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
struct PyParseGridField {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
struct PyParseIntegerList {
  std::string err_str;
  bool exact_size;
  std::string delim;
  bool delim_found;
  std::optional<std::string> open_delim;
  std::optional<std::string> separator;
  std::optional<std::string> close_delim;
  std::optional<int> default_value;
  bool is_ok;
};
struct PyParseIntegerList2 : public Bmad::ParseIntegerList2 {
  std::optional<int> num_expected;
  std::optional<std::string> open_delim;
  std::optional<std::string> separator;
  std::optional<std::string> close_delim;
  std::optional<int> default_value;
  PyParseIntegerList2(
      Bmad::ParseIntegerList2 _base,
      std::optional<int> num_expected,
      std::optional<std::string> open_delim,
      std::optional<std::string> separator,
      std::optional<std::string> close_delim,
      std::optional<int> default_value)
      : Bmad::ParseIntegerList2(std::move(_base)),
        num_expected(num_expected),
        open_delim(open_delim),
        separator(separator),
        close_delim(close_delim),
        default_value(default_value) {}
};
struct PyParseRealList2 : public Bmad::ParseRealList2 {
  std::optional<int> num_expected;
  std::optional<std::string> open_brace;
  std::optional<std::string> separator;
  std::optional<std::string> close_brace;
  std::optional<double> default_value;
  std::optional<bool> single_value;
  PyParseRealList2(
      Bmad::ParseRealList2 _base,
      std::optional<int> num_expected,
      std::optional<std::string> open_brace,
      std::optional<std::string> separator,
      std::optional<std::string> close_brace,
      std::optional<double> default_value,
      std::optional<bool> single_value)
      : Bmad::ParseRealList2(std::move(_base)),
        num_expected(num_expected),
        open_brace(open_brace),
        separator(separator),
        close_brace(close_brace),
        default_value(default_value),
        single_value(single_value) {}
};
struct PyParserAddConstant {
  std::string word;
  bool redef_is_error;
};
struct PyParserCallCheck {
  std::string word;
  int ix_word;
  std::string delim;
  bool delim_found;
  bool call_found;
  std::optional<bool> err_flag;
};
struct PyParserFastIntegerRead {
  std::string delim_wanted;
  std::string err_str;
  bool is_ok;
};
struct PyParserFileStack {
  std::string how;
  std::optional<std::string> file_name_in;
  std::optional<bool> finished;
  std::optional<bool> err;
  std::optional<bool> open_file;
  std::optional<bool> abort_on_open_error;
};
struct PyParserGetInteger {
  int int_val;
  std::string word;
  int ix_word;
  std::string delim;
  bool delim_found;
  bool err;
  std::optional<std::string> str1;
  std::optional<std::string> str2;
};
struct PyParserGetLogical {
  std::string attrib_name;
  bool this_logic;
  std::string ele_name;
  std::string delim;
  bool delim_found;
  bool err;
};
struct PyParserPrintLine {
  bool end_of_file;
};
struct PyParserReadLrWake {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
struct PyParserReadSrWake {
  std::string delim;
  bool delim_found;
  bool err_flag;
};
struct PyParticleIsMovingBackwards {
  bool is_moving_backwards;
};
struct PyParticleIsMovingForward {
  bool is_moving_forward;
};
struct PyParticleRfTime {
  long double time;
};
struct PyPatchFlipsPropagationDirection {
  bool is_flip;
};
struct PyPatchLength {
  double length;
};
struct PyPhotonAddToDetectorStatistics {
  std::optional<int> ix_pt;
  std::optional<int> iy_pt;
};

struct PyPhotonTargetCornerCalc {
  TargetPointProxy corner;
  double x_lim;
  double y_lim;
  double z_lim;
};
struct PyPhysicalEleEnd {
  int physical_end;
};

struct PyPointerToSurfaceDisplacementPt {
  SurfaceDisplacementPtProxy pt;
  double x;
  double y;
  std::optional<int> ix;
  std::optional<int> iy;
  std::optional<double> xx;
  std::optional<double> yy;
};

struct PyPointerToSurfaceSegmentedPt {
  SurfaceSegmentedPtProxy pt;
  double x;
  double y;
  std::optional<int> ix;
  std::optional<int> iy;
  std::optional<double> xx;
  std::optional<double> yy;
};
