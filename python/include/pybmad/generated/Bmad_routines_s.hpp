#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_s(py::module& m);

struct PySBodyCalc {
  double s_body;
};

struct PyScAdaptiveStep {
  double dt_next;
  bool include_image;
  double dt_step;
};

struct PyScStep {
  int n_emit;
  bool include_image;
};
struct PySetFlagsForChangedIntegerAttribute {
  int attrib;
};
struct PySetFlagsForChangedLogicalAttribute {
  bool attrib;
};
struct PySetFlagsForChangedRealAttribute {
  std::optional<double> attrib;
};
struct PySetFringeOnOff {
  double fringe_at;
};
struct PySetPtcQuiet {
  int old_val;
};
struct PySetPtcVerbose {
  bool on;
};
struct PySetTune {
  bool ok;
};
struct PySignificantDifference {
  bool is_different;
};
struct PySkipEleBlender {
  bool skip;
};
struct PySolQuadMat6Calc {
  double ks_in;
  double k1_in;
};
struct PySpinOmega {
  int sign_z_vel;
  std::optional<bool> phase_space_coords;
};
struct PyStreamEleEnd {
  int stream_end;
};
struct PyStrongBeamStrength {
  double strength;
};
struct PySurfaceGridDisplacement {
  double x;
  double y;
};
