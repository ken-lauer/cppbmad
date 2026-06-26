#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_material_struct(nb::module_ &m, nb::class_<MaterialStruct> &class_);
void init_mode3_struct(nb::module_ &m, nb::class_<Mode3Struct> &class_);
void init_mode_info_struct(nb::module_ &m, nb::class_<ModeInfoStruct> &class_);
void init_multipass_all_info_struct(nb::module_ &m, nb::class_<MultipassAllInfoStruct> &class_);
void init_multipass_branch_info_struct(
    nb::module_ &m,
    nb::class_<MultipassBranchInfoStruct> &class_
);
void init_multipass_ele_info_struct(nb::module_ &m, nb::class_<MultipassEleInfoStruct> &class_);
void init_multipass_lord_info_struct(nb::module_ &m, nb::class_<MultipassLordInfoStruct> &class_);
void init_multipole_cache_struct(nb::module_ &m, nb::class_<MultipoleCacheStruct> &class_);
void init_mad_energy_struct(nb::module_ &m, nb::class_<MadEnergyStruct> &class_);
void init_mad_map_struct(nb::module_ &m, nb::class_<MadMapStruct> &class_);
void init_mesh3d_struct(nb::module_ &m, nb::class_<Mesh3dStruct> &class_);
void init_molecular_component_struct(nb::module_ &m, nb::class_<MolecularComponentStruct> &class_);
void init_momentum_aperture_struct(nb::module_ &m, nb::class_<MomentumApertureStruct> &class_);
void init_multipass_region_branch_struct(
    nb::module_ &m,
    nb::class_<MultipassRegionBranchStruct> &class_
);
void init_multipass_region_ele_struct(nb::module_ &m, nb::class_<MultipassRegionEleStruct> &class_);
void init_multipass_region_lat_struct(nb::module_ &m, nb::class_<MultipassRegionLatStruct> &class_);
