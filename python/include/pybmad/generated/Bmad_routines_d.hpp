#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_d(py::module& m);

struct PyDampingMatrixD {
  double gamma;
  double g_tot;
  double B0;
  double B1;
  double delta;
  int species;
};
struct PyDefaultTrackingSpecies {
  int species;
};
struct PyDiffractionPlateOrMaskHitSpot {
  int ix_section;
};
struct PyDiffusionMatrixB {
  double gamma;
  double g_tot;
  int species;
};

struct PyDistanceToAperture {
  bool no_aperture_here;
  double dist;
};
struct PyDpcGivenDe {
  double pc_old;
  double mass;
  double dE;
  double dpc;
};
