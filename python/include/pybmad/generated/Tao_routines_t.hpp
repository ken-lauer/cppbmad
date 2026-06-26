#pragma once
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace nb = nanobind;

void init_Tao_routines_t(nb::module_ &m);

struct PyTaoDrawCurveData {
  bool have_data;
};
struct PyTaoDrawHistogramData {
  bool have_data;
};
struct PyTaoEleShapeInfo : public Tao::TaoEleShapeInfo {
  std::optional<int> ix_shape_min;
  PyTaoEleShapeInfo(Tao::TaoEleShapeInfo _base, std::optional<int> ix_shape_min)
      : Tao::TaoEleShapeInfo(std::move(_base))
      , ix_shape_min(ix_shape_min) {}
};
struct PyTaoNextSwitch : public Tao::TaoNextSwitch {
  std::string line;
  PyTaoNextSwitch(Tao::TaoNextSwitch _base, std::string line)
      : Tao::TaoNextSwitch(std::move(_base))
      , line(line) {}
};

struct PyTaoNextWord {
  std::string word;
  std::string line;
};
struct PyTaoPointerToEleShape : public Tao::TaoPointerToEleShape {
  std::optional<int> ix_shape_min;
  PyTaoPointerToEleShape(Tao::TaoPointerToEleShape _base, std::optional<int> ix_shape_min)
      : Tao::TaoPointerToEleShape(std::move(_base))
      , ix_shape_min(ix_shape_min) {}
};

struct PyTaoPointerToUniverseStr {
  std::optional<TaoUniverseStruct> u;
  std::string string;
};
struct PyTaoRemoveBlankCharacters {
  std::string str;
};
