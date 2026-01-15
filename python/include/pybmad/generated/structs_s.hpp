#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_space_charge_common_struct(py::module &m, py::class_<SpaceChargeCommonStruct> &class_);
void init_spin_axis_struct(py::module &m, py::class_<SpinAxisStruct> &class_);
void init_spin_orbit_map1_struct(py::module &m, py::class_<SpinOrbitMap1Struct> &class_);
void init_spin_polar_struct(py::module &m, py::class_<SpinPolarStruct> &class_);
void init_strong_beam_struct(py::module &m, py::class_<StrongBeamStruct> &class_);
void init_surface_curvature_struct(py::module &m, py::class_<SurfaceCurvatureStruct> &class_);
void init_surface_displacement_pt_struct(
    py::module &m,
    py::class_<SurfaceDisplacementPtStruct> &class_
);
void init_surface_displacement_struct(py::module &m, py::class_<SurfaceDisplacementStruct> &class_);
void init_surface_h_misalign_pt_struct(py::module &m, py::class_<SurfaceHMisalignPtStruct> &class_);
void init_surface_h_misalign_struct(py::module &m, py::class_<SurfaceHMisalignStruct> &class_);
void init_surface_segmented_pt_struct(py::module &m, py::class_<SurfaceSegmentedPtStruct> &class_);
void init_surface_segmented_struct(py::module &m, py::class_<SurfaceSegmentedStruct> &class_);
void init_spline_struct(py::module &m, py::class_<SplineStruct> &class_);
void init_summation_rdt_struct(py::module &m, py::class_<SummationRdtStruct> &class_);
