#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_parser_controller_struct(nb::module_ &m, nb::class_<ParserControllerStruct> &class_);
void init_parser_ele_struct(nb::module_ &m, nb::class_<ParserEleStruct> &class_);
void init_parser_lat_struct(nb::module_ &m, nb::class_<ParserLatStruct> &class_);
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
void init_pmd_header_struct(nb::module_ &m, nb::class_<PmdHeaderStruct> &class_);
void init_pre_tracker_struct(nb::module_ &m, nb::class_<PreTrackerStruct> &class_);
void init_ptc_branch1_struct(nb::module_ &m, nb::class_<PtcBranch1Struct> &class_);
void init_ptc_layout_pointer_struct(nb::module_ &m, nb::class_<PtcLayoutPointerStruct> &class_);
void init_ptc_normal_form_struct(nb::module_ &m, nb::class_<PtcNormalFormStruct> &class_);
void init_photon_coord_struct(nb::module_ &m, nb::class_<PhotonCoordStruct> &class_);
void init_photon_track_struct(nb::module_ &m, nb::class_<PhotonTrackStruct> &class_);
void init_photon_init_splines_struct(nb::module_ &m, nb::class_<PhotonInitSplinesStruct> &class_);
void init_photon_init_x_angle_spline_struct(
    nb::module_ &m,
    nb::class_<PhotonInitXAngleSplineStruct> &class_
);
void init_photon_init_y_angle_spline_struct(
    nb::module_ &m,
    nb::class_<PhotonInitYAngleSplineStruct> &class_
);
void init_ptc_rad_map_struct(nb::module_ &m, nb::class_<PtcRadMapStruct> &class_);
