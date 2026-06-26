#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_parser_controller_struct(nb::module_ &m, nb::class_<ParserControllerStruct> &class_);
void init_parser_ele_struct(nb::module_ &m, nb::class_<ParserEleStruct> &class_);
void init_photon_element_struct(nb::module_ &m, nb::class_<PhotonElementStruct> &class_);
void init_photon_material_struct(nb::module_ &m, nb::class_<PhotonMaterialStruct> &class_);
void init_photon_reflect_surface_struct(
    nb::module_ &m,
    nb::class_<PhotonReflectSurfaceStruct> &class_
);
void init_photon_reflect_table_struct(nb::module_ &m, nb::class_<PhotonReflectTableStruct> &class_);
void init_photon_target_struct(nb::module_ &m, nb::class_<PhotonTargetStruct> &class_);
void init_pixel_detec_struct(nb::module_ &m, nb::class_<PixelDetecStruct> &class_);
void init_pixel_pt_struct(nb::module_ &m, nb::class_<PixelPtStruct> &class_);
void init_pre_tracker_struct(nb::module_ &m, nb::class_<PreTrackerStruct> &class_);
void init_ptc_normal_form_struct(nb::module_ &m, nb::class_<PtcNormalFormStruct> &class_);
