#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_space_charge_common_struct(
    py::module& m,
    py::class_<SpaceChargeCommonProxy>& class_);
void init_spin_axis_struct(py::module& m, py::class_<SpinAxisProxy>& class_);
void init_spin_orbit_map1_struct(
    py::module& m,
    py::class_<SpinOrbitMap1Proxy>& class_);
void init_spin_polar_struct(py::module& m, py::class_<SpinPolarProxy>& class_);
void init_strong_beam_struct(
    py::module& m,
    py::class_<StrongBeamProxy>& class_);
void init_surface_curvature_struct(
    py::module& m,
    py::class_<SurfaceCurvatureProxy>& class_);
void init_surface_displacement_pt_struct(
    py::module& m,
    py::class_<SurfaceDisplacementPtProxy>& class_);
void init_surface_displacement_struct(
    py::module& m,
    py::class_<SurfaceDisplacementProxy>& class_);
void init_surface_h_misalign_pt_struct(
    py::module& m,
    py::class_<SurfaceHMisalignPtProxy>& class_);
void init_surface_h_misalign_struct(
    py::module& m,
    py::class_<SurfaceHMisalignProxy>& class_);
void init_surface_segmented_pt_struct(
    py::module& m,
    py::class_<SurfaceSegmentedPtProxy>& class_);
void init_surface_segmented_struct(
    py::module& m,
    py::class_<SurfaceSegmentedProxy>& class_);
void init_spline_struct(py::module& m, py::class_<SplineProxy>& class_);
void init_summation_rdt_struct(
    py::module& m,
    py::class_<SummationRdtProxy>& class_);
