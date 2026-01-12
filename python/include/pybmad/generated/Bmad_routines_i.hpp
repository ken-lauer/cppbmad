#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_i(py::module& m);

struct PyIbsMatrixC {
  bool tail_cut;
  double tau;
  double energy;
  double n_part;
  int species;
};
struct PyIgfcoulombfun {
  double u;
  double v;
  double w;
  double gam;
  double dx;
  double dy;
  double dz;
  double res;
};
struct PyIgfexfun {
  double u;
  double v;
  double w;
  double gam;
  double dx;
  double dy;
  double dz;
  double res;
};
struct PyIgfeyfun {
  double u;
  double v;
  double w;
  double gam;
  double dx;
  double dy;
  double dz;
  double res;
};
struct PyIgfezfun {
  double u;
  double v;
  double w;
  double gam;
  double dx;
  double dy;
  double dz;
  double res;
};
struct PyInitAttributeName1 {
  bool is_ok;
};
struct PyInitBeamDistribution : public Bmad::InitBeamDistribution {
  std::optional<bool> conserve_momentum;
  PyInitBeamDistribution(
      Bmad::InitBeamDistribution _base,
      std::optional<bool> conserve_momentum)
      : Bmad::InitBeamDistribution(std::move(_base)),
        conserve_momentum(conserve_momentum) {}
};
struct PyInitBunchDistribution : public Bmad::InitBunchDistribution {
  std::optional<bool> conserve_momentum;
  PyInitBunchDistribution(
      Bmad::InitBunchDistribution _base,
      std::optional<bool> conserve_momentum)
      : Bmad::InitBunchDistribution(std::move(_base)),
        conserve_momentum(conserve_momentum) {}
};
struct PyInitSurfaceSegment {
  int ix;
  int iy;
};
struct PyIntegrandBase {
  double func_retval__;
};
struct PyIntegrationTimerEle {
  double tol;
};
