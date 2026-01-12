#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_r(py::module& m);

struct PyRadGIntegrals {
  FixedArray1D<Real, 2> int_g;
  double int_g2;
  double int_g3;
};
struct PyRadiationIntegrals : public Bmad::RadiationIntegrals {
  std::optional<int> ix_cache;
  PyRadiationIntegrals(
      Bmad::RadiationIntegrals _base,
      std::optional<int> ix_cache)
      : Bmad::RadiationIntegrals(std::move(_base)), ix_cache(ix_cache) {}
};

struct PyRamperValue {
  bool err_flag;
  double value;
};
struct PyRchomp {
  double rel;
  int plc;
  std::string out;
};
struct PyReAllocateWall3dSectionArray {
  std::optional<bool> exact;
};
struct PyReAllocateWall3dVertexArray {
  std::optional<bool> exact;
};
struct PyReStrQp {
  long double rel;
  std::string str_out;
};
struct PyReStrRp {
  double rel;
  std::string str_out;
};
struct PyReadBeamFile : public Bmad::ReadBeamFile {
  std::optional<bool> conserve_momentum;
  PyReadBeamFile(
      Bmad::ReadBeamFile _base,
      std::optional<bool> conserve_momentum)
      : Bmad::ReadBeamFile(std::move(_base)),
        conserve_momentum(conserve_momentum) {}
};
struct PyReallocateBeam {
  std::optional<bool> extend;
};
struct PyRelTrackingChargeToMass {
  double rel_charge;
};
struct PyRelativeModeFlip {
  bool func_retval__;
};
struct PyReleaseRadIntCache {
  int ix_cache;
};
struct PyRfIsOn {
  bool is_on;
};
struct PyRfRefTimeOffset {
  double time;
};
struct PyRfun {
  double u;
  double v;
  double w;
  double gam;
  double a;
  double b;
  double hz;
  int i;
  int j;
  double res;
};
struct PyRkAdaptiveTimeStep {
  int t_dir;
  double rf_time;
  double dt_try;
  double dt_did;
  double dt_next;
  bool err_flag;
};

struct PyRkTimeStep1 {
  FixedArray1D<Real, 10> r_err;
  bool err_flag;
  std::optional<bool> print_err;
};
struct PyRotate3 {
  double angle;
};
struct PyRotateFieldZx {
  double theta;
};
