#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_c(py::module& m);

struct PyCalcBunchSigmaMatrixEtc {
  BunchParamsProxy bunch_params;
  std::optional<bool> is_time_coords;
};
struct PyCheckForSuperimposeProblem {
  bool err_flag;
  bool wrap;
};
struct PyChromCalc : public Bmad::ChromCalc {
  double delta_e;
  PyChromCalc(Bmad::ChromCalc _base, double delta_e)
      : Bmad::ChromCalc(std::move(_base)), delta_e(delta_e) {}
};

struct PyChromTune {
  bool err_flag;
  double delta_e;
};
struct PyClassicalRadius {
  double radius;
};
struct PyCmplxReStr {
  std::complex<double> cmp;
  std::string str_out;
};
struct PyComplexTaylorCoef1 {
  std::complex<double> coef;
};
struct PyComplexTaylorCoef2 {
  std::complex<double> coef;
};
struct PyConvertLocalCartesianToLocalCurvilinear {
  double x;
  double z;
  double g;
  double xout;
  double sout;
};
struct PyConvertLocalCurvilinearToLocalCartesian {
  double x;
  double s;
  double g;
  double xout;
  double zout;
};

struct PyCoordStateName {
  std::string state_str;
  std::optional<bool> one_word;
};
struct PyCoordsRelativeToFloor {
  std::optional<double> theta;
  std::optional<double> phi;
  std::optional<double> psi;
};
struct PyCoulombfun {
  double u;
  double v;
  double w;
  double gam;
  double res;
};
struct PyCreateConcatenatedWall3d {
  bool err;
};
struct PyCreateGirder {
  bool err_flag;
};
