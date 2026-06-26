#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_cartesian_map_struct(nb::module_ &m, nb::class_<CartesianMapStruct> &class_);
void init_cartesian_map_term1_struct(nb::module_ &m, nb::class_<CartesianMapTerm1Struct> &class_);
void init_cartesian_map_term_struct(nb::module_ &m, nb::class_<CartesianMapTermStruct> &class_);
void init_complex_taylor_struct(nb::module_ &m, nb::class_<ComplexTaylorStruct> &class_);
void init_complex_taylor_term_struct(nb::module_ &m, nb::class_<ComplexTaylorTermStruct> &class_);
void init_control_ramp1_struct(nb::module_ &m, nb::class_<ControlRamp1Struct> &class_);
void init_control_struct(nb::module_ &m, nb::class_<ControlStruct> &class_);
void init_control_var1_struct(nb::module_ &m, nb::class_<ControlVar1Struct> &class_);
void init_controller_struct(nb::module_ &m, nb::class_<ControllerStruct> &class_);
void init_converter_dir_1D_struct(nb::module_ &m, nb::class_<ConverterDir1dStruct> &class_);
void init_converter_dir_2D_struct(nb::module_ &m, nb::class_<ConverterDir2dStruct> &class_);
void init_converter_dir_coef_struct(nb::module_ &m, nb::class_<ConverterDirCoefStruct> &class_);
void init_converter_direction_out_struct(
    nb::module_ &m,
    nb::class_<ConverterDirectionOutStruct> &class_
);
void init_converter_distribution_struct(
    nb::module_ &m,
    nb::class_<ConverterDistributionStruct> &class_
);
void init_converter_prob_pc_r_struct(nb::module_ &m, nb::class_<ConverterProbPcRStruct> &class_);
void init_converter_struct(nb::module_ &m, nb::class_<ConverterStruct> &class_);
void init_converter_sub_distribution_struct(
    nb::module_ &m,
    nb::class_<ConverterSubDistributionStruct> &class_
);
void init_coord_array_struct(nb::module_ &m, nb::class_<CoordArrayStruct> &class_);
void init_coord_struct(nb::module_ &m, nb::class_<CoordStruct> &class_);
void init_cylindrical_map_struct(nb::module_ &m, nb::class_<CylindricalMapStruct> &class_);
void init_cylindrical_map_term1_struct(
    nb::module_ &m,
    nb::class_<CylindricalMapTerm1Struct> &class_
);
void init_cylindrical_map_term_struct(nb::module_ &m, nb::class_<CylindricalMapTermStruct> &class_);
void init_csr_bunch_slice_struct(nb::module_ &m, nb::class_<CsrBunchSliceStruct> &class_);
void init_csr_ele_info_struct(nb::module_ &m, nb::class_<CsrEleInfoStruct> &class_);
void init_csr_kick1_struct(nb::module_ &m, nb::class_<CsrKick1Struct> &class_);
void init_csr_particle_position_struct(
    nb::module_ &m,
    nb::class_<CsrParticlePositionStruct> &class_
);
void init_csr_struct(nb::module_ &m, nb::class_<CsrStruct> &class_);
void init_cmplx_field1_at_2D_pt_struct(nb::module_ &m, nb::class_<CmplxField1At2dPtStruct> &class_);
void init_cmplx_field1_at_3D_pt_struct(nb::module_ &m, nb::class_<CmplxField1At3dPtStruct> &class_);
void init_cmplx_field_at_2D_box_struct(nb::module_ &m, nb::class_<CmplxFieldAt2dBoxStruct> &class_);
void init_cmplx_field_at_3D_box_struct(nb::module_ &m, nb::class_<CmplxFieldAt3dBoxStruct> &class_);
void init_crystal_param_struct(nb::module_ &m, nb::class_<CrystalParamStruct> &class_);
