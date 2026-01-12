#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_t(py::module& m);

struct PyTaperMagStrengths {
  std::optional<bool> err_flag;
};
struct PyTargetMinMaxCalc {
  double y_min;
  double y_max;
  double phi_min;
  double phi_max;
};
struct PyToFieldmapCoords {
  double x;
  double y;
  double z;
  double cos_ang;
  double sin_ang;
};
struct PyTouschekRate1Zap {
  double rate;
  std::optional<int> ix;
  std::optional<double> s;
};
struct PyTrack1Spin : public Bmad::Track1Spin {
  std::optional<bool> make_quaternion;
  PyTrack1Spin(Bmad::Track1Spin _base, std::optional<bool> make_quaternion)
      : Bmad::Track1Spin(std::move(_base)), make_quaternion(make_quaternion) {}
};
struct PyTrack1TimeRungeKutta : public Bmad::Track1TimeRungeKutta {
  std::optional<double> dt_step;
  PyTrack1TimeRungeKutta(
      Bmad::Track1TimeRungeKutta _base,
      std::optional<double> dt_step)
      : Bmad::Track1TimeRungeKutta(std::move(_base)), dt_step(dt_step) {}
};
struct PyTrackADrift {
  std::optional<double> time;
};

struct PyTrackAMatch {
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
  std::optional<bool> err_flag;
};

struct PyTrackAPickup {
  std::optional<FixedArray2D<Real, 6, 6>> mat6;
  std::optional<bool> err_flag;
};

struct PyTwiss3AtStart {
  FixedArray1D<Real, 3> tune3;
  bool err_flag;
};
struct PyTwiss3Propagate1 {
  bool err_flag;
};
