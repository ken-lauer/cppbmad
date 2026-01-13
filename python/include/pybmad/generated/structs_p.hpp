#pragma once
#include <pybind11/pybind11.h>
#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_photon_element_struct(
    py::module& m,
    py::class_<PhotonElementProxy>& class_);
void init_photon_material_struct(
    py::module& m,
    py::class_<PhotonMaterialProxy>& class_);
void init_photon_reflect_surface_struct(
    py::module& m,
    py::class_<PhotonReflectSurfaceProxy>& class_);
void init_photon_reflect_table_struct(
    py::module& m,
    py::class_<PhotonReflectTableProxy>& class_);
void init_photon_target_struct(
    py::module& m,
    py::class_<PhotonTargetProxy>& class_);
void init_pixel_detec_struct(
    py::module& m,
    py::class_<PixelDetecProxy>& class_);
void init_pixel_pt_struct(py::module& m, py::class_<PixelPtProxy>& class_);
void init_pre_tracker_struct(
    py::module& m,
    py::class_<PreTrackerProxy>& class_);
void init_ptc_normal_form_struct(
    py::module& m,
    py::class_<PtcNormalFormProxy>& class_);
