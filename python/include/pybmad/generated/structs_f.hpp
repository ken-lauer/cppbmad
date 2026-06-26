#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_floor_position_struct(nb::module_ &m, nb::class_<FloorPositionStruct> &class_);
void init_foil_struct(nb::module_ &m, nb::class_<FoilStruct> &class_);
void init_fringe_field_info_struct(nb::module_ &m, nb::class_<FringeFieldInfoStruct> &class_);
void init_field1_at_2D_pt_struct(nb::module_ &m, nb::class_<Field1At2dPtStruct> &class_);
void init_field1_at_3D_pt_struct(nb::module_ &m, nb::class_<Field1At3dPtStruct> &class_);
void init_field_at_2D_box_struct(nb::module_ &m, nb::class_<FieldAt2dBoxStruct> &class_);
void init_field_at_3D_box_struct(nb::module_ &m, nb::class_<FieldAt3dBoxStruct> &class_);
void init_fibre(nb::module_ &m, nb::class_<Fibre> &class_);
