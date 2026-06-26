#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_space_charge_common_struct(nb::module_ &m, nb::class_<SpaceChargeCommonStruct> &class_);
void init_spin_axis_struct(nb::module_ &m, nb::class_<SpinAxisStruct> &class_);
void init_spin_orbit_map1_struct(nb::module_ &m, nb::class_<SpinOrbitMap1Struct> &class_);
void init_spin_polar_struct(nb::module_ &m, nb::class_<SpinPolarStruct> &class_);
void init_strong_beam_struct(nb::module_ &m, nb::class_<StrongBeamStruct> &class_);
void init_surface_curvature_struct(nb::module_ &m, nb::class_<SurfaceCurvatureStruct> &class_);
void init_surface_displacement_pt_struct(
    nb::module_ &m,
    nb::class_<SurfaceDisplacementPtStruct> &class_
);
void init_surface_displacement_struct(
    nb::module_ &m,
    nb::class_<SurfaceDisplacementStruct> &class_
);
void init_surface_h_misalign_pt_struct(
    nb::module_ &m,
    nb::class_<SurfaceHMisalignPtStruct> &class_
);
void init_surface_h_misalign_struct(nb::module_ &m, nb::class_<SurfaceHMisalignStruct> &class_);
void init_surface_segmented_pt_struct(nb::module_ &m, nb::class_<SurfaceSegmentedPtStruct> &class_);
void init_surface_segmented_struct(nb::module_ &m, nb::class_<SurfaceSegmentedStruct> &class_);
void init_spline_struct(nb::module_ &m, nb::class_<SplineStruct> &class_);
void init_summation_rdt_struct(nb::module_ &m, nb::class_<SummationRdtStruct> &class_);
