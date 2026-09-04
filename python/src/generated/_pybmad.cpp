#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <string>

#include "bmad/json.hpp"
#include "pybmad/bmad_hooks.hpp"
#include "pybmad/common_structs.hpp"
#include "pybmad/generated/init.hpp"
#include "pybmad/tao_hooks.hpp"
#include "pybmad/util.hpp"

using namespace Bmad;
using namespace Pybmad;

namespace nb = nanobind;
using namespace nanobind::literals;

#include "pybmad/generated/routines.hpp"
#include "pybmad/generated/structs.hpp"

NB_MODULE(_pybmad, m) {
  // Generated definitions:
  m.doc() = "pybmad";

  bind_standard_arrays(m);

  // Structures
  auto py_AcKickerFreqStruct = nb::class_<AcKickerFreqStruct>(
      m,
      "AcKickerFreqStruct",
      "Fortran struct: ac_kicker_freq_struct"
  );
  auto py_AcKickerStruct =
      nb::class_<AcKickerStruct>(m, "AcKickerStruct", "Fortran struct: ac_kicker_struct");
  auto py_AcKickerTimeStruct = nb::class_<AcKickerTimeStruct>(
      m,
      "AcKickerTimeStruct",
      "Fortran struct: ac_kicker_time_struct"
  );
  auto py_AllEncompassingStruct = nb::class_<AllEncompassingStruct>(
      m,
      "AllEncompassingStruct",
      "Fortran struct: all_encompassing_struct"
  );
  auto py_AllPointerStruct =
      nb::class_<AllPointerStruct>(m, "AllPointerStruct", "Fortran struct: all_pointer_struct");
  auto py_AnormalModeStruct =
      nb::class_<AnormalModeStruct>(m, "AnormalModeStruct", "Fortran struct: anormal_mode_struct");
  auto py_ApertureParamStruct = nb::class_<ApertureParamStruct>(
      m,
      "ApertureParamStruct",
      "Fortran struct: aperture_param_struct"
  );
  auto py_AperturePointStruct = nb::class_<AperturePointStruct>(
      m,
      "AperturePointStruct",
      "Fortran struct: aperture_point_struct"
  );
  auto py_ApertureScanStruct = nb::class_<ApertureScanStruct>(
      m,
      "ApertureScanStruct",
      "Fortran struct: aperture_scan_struct"
  );
  auto py_AstraLatticeParamStruct = nb::class_<AstraLatticeParamStruct>(
      m,
      "AstraLatticeParamStruct",
      "Fortran struct: astra_lattice_param_struct"
  );
  auto py_BaseLineEleStruct =
      nb::class_<BaseLineEleStruct>(m, "BaseLineEleStruct", "Fortran struct: base_line_ele_struct");
  auto py_BbuBeamStruct =
      nb::class_<BbuBeamStruct>(m, "BbuBeamStruct", "Fortran struct: bbu_beam_struct");
  auto py_BbuParamStruct =
      nb::class_<BbuParamStruct>(m, "BbuParamStruct", "Fortran struct: bbu_param_struct");
  auto py_BbuStageStruct =
      nb::class_<BbuStageStruct>(m, "BbuStageStruct", "Fortran struct: bbu_stage_struct");
  auto py_BeamInitStruct =
      nb::class_<BeamInitStruct>(m, "BeamInitStruct", "Fortran struct: beam_init_struct");
  auto py_BeamStruct = nb::class_<BeamStruct>(m, "BeamStruct", "Fortran struct: beam_struct");
  auto py_BicubicCmplxCoefStruct = nb::class_<BicubicCmplxCoefStruct>(
      m,
      "BicubicCmplxCoefStruct",
      "Fortran struct: bicubic_cmplx_coef_struct"
  );
  auto py_BicubicCoefStruct =
      nb::class_<BicubicCoefStruct>(m, "BicubicCoefStruct", "Fortran struct: bicubic_coef_struct");
  auto py_BinStruct = nb::class_<BinStruct>(m, "BinStruct", "Fortran struct: bin_struct");
  auto py_BmadCommonStruct =
      nb::class_<BmadCommonStruct>(m, "BmadCommonStruct", "Fortran struct: bmad_common_struct");
  auto py_BmadNormalFormStruct = nb::class_<BmadNormalFormStruct>(
      m,
      "BmadNormalFormStruct",
      "Fortran struct: bmad_normal_form_struct"
  );
  auto py_BookkeepingStateStruct = nb::class_<BookkeepingStateStruct>(
      m,
      "BookkeepingStateStruct",
      "Fortran struct: bookkeeping_state_struct"
  );
  auto py_BpmPhaseCouplingStruct = nb::class_<BpmPhaseCouplingStruct>(
      m,
      "BpmPhaseCouplingStruct",
      "Fortran struct: bpm_phase_coupling_struct"
  );
  auto py_BranchPointerStruct = nb::class_<BranchPointerStruct>(
      m,
      "BranchPointerStruct",
      "Fortran struct: branch_pointer_struct"
  );
  auto py_BranchStruct =
      nb::class_<BranchStruct>(m, "BranchStruct", "Fortran struct: branch_struct");
  auto py_BunchParamsStruct =
      nb::class_<BunchParamsStruct>(m, "BunchParamsStruct", "Fortran struct: bunch_params_struct");
  auto py_BunchStruct = nb::class_<BunchStruct>(m, "BunchStruct", "Fortran struct: bunch_struct");
  auto py_BunchTrackStruct =
      nb::class_<BunchTrackStruct>(m, "BunchTrackStruct", "Fortran struct: bunch_track_struct");
  auto py_CartesianMapStruct = nb::class_<CartesianMapStruct>(
      m,
      "CartesianMapStruct",
      "Fortran struct: cartesian_map_struct"
  );
  auto py_CartesianMapTerm1Struct = nb::class_<CartesianMapTerm1Struct>(
      m,
      "CartesianMapTerm1Struct",
      "Fortran struct: cartesian_map_term1_struct"
  );
  auto py_CartesianMapTermStruct = nb::class_<CartesianMapTermStruct>(
      m,
      "CartesianMapTermStruct",
      "Fortran struct: cartesian_map_term_struct"
  );
  auto py_CmplxField1At2DPtStruct = nb::class_<CmplxField1At2dPtStruct>(
      m,
      "CmplxField1At2DPtStruct",
      "Fortran struct: cmplx_field1_at_2D_pt_struct"
  );
  auto py_CmplxField1At3DPtStruct = nb::class_<CmplxField1At3dPtStruct>(
      m,
      "CmplxField1At3DPtStruct",
      "Fortran struct: cmplx_field1_at_3D_pt_struct"
  );
  auto py_CmplxFieldAt2DBoxStruct = nb::class_<CmplxFieldAt2dBoxStruct>(
      m,
      "CmplxFieldAt2DBoxStruct",
      "Fortran struct: cmplx_field_at_2D_box_struct"
  );
  auto py_CmplxFieldAt3DBoxStruct = nb::class_<CmplxFieldAt3dBoxStruct>(
      m,
      "CmplxFieldAt3DBoxStruct",
      "Fortran struct: cmplx_field_at_3D_box_struct"
  );
  auto py_ComplexTaylorStruct = nb::class_<ComplexTaylorStruct>(
      m,
      "ComplexTaylorStruct",
      "Fortran struct: complex_taylor_struct"
  );
  auto py_ComplexTaylorTermStruct = nb::class_<ComplexTaylorTermStruct>(
      m,
      "ComplexTaylorTermStruct",
      "Fortran struct: complex_taylor_term_struct"
  );
  auto py_ControlRamp1Struct = nb::class_<ControlRamp1Struct>(
      m,
      "ControlRamp1Struct",
      "Fortran struct: control_ramp1_struct"
  );
  auto py_ControlStruct =
      nb::class_<ControlStruct>(m, "ControlStruct", "Fortran struct: control_struct");
  auto py_ControlVar1Struct =
      nb::class_<ControlVar1Struct>(m, "ControlVar1Struct", "Fortran struct: control_var1_struct");
  auto py_ControllerStruct =
      nb::class_<ControllerStruct>(m, "ControllerStruct", "Fortran struct: controller_struct");
  auto py_ConverterDir1DStruct = nb::class_<ConverterDir1dStruct>(
      m,
      "ConverterDir1DStruct",
      "Fortran struct: converter_dir_1D_struct"
  );
  auto py_ConverterDir2DStruct = nb::class_<ConverterDir2dStruct>(
      m,
      "ConverterDir2DStruct",
      "Fortran struct: converter_dir_2D_struct"
  );
  auto py_ConverterDirCoefStruct = nb::class_<ConverterDirCoefStruct>(
      m,
      "ConverterDirCoefStruct",
      "Fortran struct: converter_dir_coef_struct"
  );
  auto py_ConverterDirectionOutStruct = nb::class_<ConverterDirectionOutStruct>(
      m,
      "ConverterDirectionOutStruct",
      "Fortran struct: converter_direction_out_struct"
  );
  auto py_ConverterDistributionStruct = nb::class_<ConverterDistributionStruct>(
      m,
      "ConverterDistributionStruct",
      "Fortran struct: converter_distribution_struct"
  );
  auto py_ConverterProbPcRStruct = nb::class_<ConverterProbPcRStruct>(
      m,
      "ConverterProbPcRStruct",
      "Fortran struct: converter_prob_pc_r_struct"
  );
  auto py_ConverterStruct =
      nb::class_<ConverterStruct>(m, "ConverterStruct", "Fortran struct: converter_struct");
  auto py_ConverterSubDistributionStruct = nb::class_<ConverterSubDistributionStruct>(
      m,
      "ConverterSubDistributionStruct",
      "Fortran struct: converter_sub_distribution_struct"
  );
  auto py_CoordArrayStruct =
      nb::class_<CoordArrayStruct>(m, "CoordArrayStruct", "Fortran struct: coord_array_struct");
  auto py_CoordStruct = nb::class_<CoordStruct>(m, "CoordStruct", "Fortran struct: coord_struct");
  auto py_CrystalParamStruct = nb::class_<CrystalParamStruct>(
      m,
      "CrystalParamStruct",
      "Fortran struct: crystal_param_struct"
  );
  auto py_CsrBunchSliceStruct = nb::class_<CsrBunchSliceStruct>(
      m,
      "CsrBunchSliceStruct",
      "Fortran struct: csr_bunch_slice_struct"
  );
  auto py_CsrEleInfoStruct =
      nb::class_<CsrEleInfoStruct>(m, "CsrEleInfoStruct", "Fortran struct: csr_ele_info_struct");
  auto py_CsrKick1Struct =
      nb::class_<CsrKick1Struct>(m, "CsrKick1Struct", "Fortran struct: csr_kick1_struct");
  auto py_CsrParticlePositionStruct = nb::class_<CsrParticlePositionStruct>(
      m,
      "CsrParticlePositionStruct",
      "Fortran struct: csr_particle_position_struct"
  );
  auto py_CsrStruct = nb::class_<CsrStruct>(m, "CsrStruct", "Fortran struct: csr_struct");
  auto py_CylindricalMapStruct = nb::class_<CylindricalMapStruct>(
      m,
      "CylindricalMapStruct",
      "Fortran struct: cylindrical_map_struct"
  );
  auto py_CylindricalMapTerm1Struct = nb::class_<CylindricalMapTerm1Struct>(
      m,
      "CylindricalMapTerm1Struct",
      "Fortran struct: cylindrical_map_term1_struct"
  );
  auto py_CylindricalMapTermStruct = nb::class_<CylindricalMapTermStruct>(
      m,
      "CylindricalMapTermStruct",
      "Fortran struct: cylindrical_map_term_struct"
  );
  auto py_DiffuseParamStruct = nb::class_<DiffuseParamStruct>(
      m,
      "DiffuseParamStruct",
      "Fortran struct: diffuse_param_struct"
  );
  auto py_DoLoopStruct =
      nb::class_<DoLoopStruct>(m, "DoLoopStruct", "Fortran struct: do_loop_struct");
  auto py_EleAttributeStruct = nb::class_<EleAttributeStruct>(
      m,
      "EleAttributeStruct",
      "Fortran struct: ele_attribute_struct"
  );
  auto py_ElePointerStruct =
      nb::class_<ElePointerStruct>(m, "ElePointerStruct", "Fortran struct: ele_pointer_struct");
  auto py_EleStruct = nb::class_<EleStruct>(m, "EleStruct", "Fortran struct: ele_struct");
  auto py_EllipseBeamInitStruct = nb::class_<EllipseBeamInitStruct>(
      m,
      "EllipseBeamInitStruct",
      "Fortran struct: ellipse_beam_init_struct"
  );
  auto py_EmFieldStruct =
      nb::class_<EmFieldStruct>(m, "EmFieldStruct", "Fortran struct: em_field_struct");
  auto py_ExpressionAtomStruct = nb::class_<ExpressionAtomStruct>(
      m,
      "ExpressionAtomStruct",
      "Fortran struct: expression_atom_struct"
  );
  auto py_ExpressionTreeStruct = nb::class_<ExpressionTreeStruct>(
      m,
      "ExpressionTreeStruct",
      "Fortran struct: expression_tree_struct"
  );
  auto py_ExtraParsingInfoStruct = nb::class_<ExtraParsingInfoStruct>(
      m,
      "ExtraParsingInfoStruct",
      "Fortran struct: extra_parsing_info_struct"
  );
  auto py_Fibre = nb::class_<Fibre>(m, "Fibre", "Fortran struct: fibre");
  auto py_Field1At2DPtStruct = nb::class_<Field1At2dPtStruct>(
      m,
      "Field1At2DPtStruct",
      "Fortran struct: field1_at_2D_pt_struct"
  );
  auto py_Field1At3DPtStruct = nb::class_<Field1At3dPtStruct>(
      m,
      "Field1At3DPtStruct",
      "Fortran struct: field1_at_3D_pt_struct"
  );
  auto py_FieldAt2DBoxStruct = nb::class_<FieldAt2dBoxStruct>(
      m,
      "FieldAt2DBoxStruct",
      "Fortran struct: field_at_2D_box_struct"
  );
  auto py_FieldAt3DBoxStruct = nb::class_<FieldAt3dBoxStruct>(
      m,
      "FieldAt3DBoxStruct",
      "Fortran struct: field_at_3D_box_struct"
  );
  auto py_FloorPositionStruct = nb::class_<FloorPositionStruct>(
      m,
      "FloorPositionStruct",
      "Fortran struct: floor_position_struct"
  );
  auto py_FoilStruct = nb::class_<FoilStruct>(m, "FoilStruct", "Fortran struct: foil_struct");
  auto py_FringeFieldInfoStruct = nb::class_<FringeFieldInfoStruct>(
      m,
      "FringeFieldInfoStruct",
      "Fortran struct: fringe_field_info_struct"
  );
  auto py_GenGradCurveStruct = nb::class_<GenGradCurveStruct>(
      m,
      "GenGradCurveStruct",
      "Fortran struct: gen_grad_curve_struct"
  );
  auto py_GenGradientsStruct = nb::class_<GenGradientsStruct>(
      m,
      "GenGradientsStruct",
      "Fortran struct: gen_gradients_struct"
  );
  auto py_GeneralBinStruct =
      nb::class_<GeneralBinStruct>(m, "GeneralBinStruct", "Fortran struct: general_bin_struct");
  auto py_GgTaylorStruct =
      nb::class_<GgTaylorStruct>(m, "GgTaylorStruct", "Fortran struct: gg_taylor_struct");
  auto py_GgTaylorTermStruct = nb::class_<GgTaylorTermStruct>(
      m,
      "GgTaylorTermStruct",
      "Fortran struct: gg_taylor_term_struct"
  );
  auto py_GptLatParamStruct =
      nb::class_<GptLatParamStruct>(m, "GptLatParamStruct", "Fortran struct: gpt_lat_param_struct");
  auto py_GridBeamInitStruct = nb::class_<GridBeamInitStruct>(
      m,
      "GridBeamInitStruct",
      "Fortran struct: grid_beam_init_struct"
  );
  auto py_GridFieldPt1Struct = nb::class_<GridFieldPt1Struct>(
      m,
      "GridFieldPt1Struct",
      "Fortran struct: grid_field_pt1_struct"
  );
  auto py_GridFieldPtStruct =
      nb::class_<GridFieldPtStruct>(m, "GridFieldPtStruct", "Fortran struct: grid_field_pt_struct");
  auto py_GridFieldStruct =
      nb::class_<GridFieldStruct>(m, "GridFieldStruct", "Fortran struct: grid_field_struct");
  auto py_HighEnergySpaceChargeStruct = nb::class_<HighEnergySpaceChargeStruct>(
      m,
      "HighEnergySpaceChargeStruct",
      "Fortran struct: high_energy_space_charge_struct"
  );
  auto py_IbsLifetimeStruct =
      nb::class_<IbsLifetimeStruct>(m, "IbsLifetimeStruct", "Fortran struct: ibs_lifetime_struct");
  auto py_IbsMaxratioStruct =
      nb::class_<IbsMaxratioStruct>(m, "IbsMaxratioStruct", "Fortran struct: ibs_maxratio_struct");
  auto py_IbsSimParamStruct =
      nb::class_<IbsSimParamStruct>(m, "IbsSimParamStruct", "Fortran struct: ibs_sim_param_struct");
  auto py_IbsStruct = nb::class_<IbsStruct>(m, "IbsStruct", "Fortran struct: ibs_struct");
  auto py_Interval1CoefStruct = nb::class_<Interval1CoefStruct>(
      m,
      "Interval1CoefStruct",
      "Fortran struct: interval1_coef_struct"
  );
  auto py_KvBeamInitStruct =
      nb::class_<KvBeamInitStruct>(m, "KvBeamInitStruct", "Fortran struct: kv_beam_init_struct");
  auto py_LatEleLocStruct =
      nb::class_<LatEleLocStruct>(m, "LatEleLocStruct", "Fortran struct: lat_ele_loc_struct");
  auto py_LatEleOrder1Struct = nb::class_<LatEleOrder1Struct>(
      m,
      "LatEleOrder1Struct",
      "Fortran struct: lat_ele_order1_struct"
  );
  auto py_LatEleOrderArrayStruct = nb::class_<LatEleOrderArrayStruct>(
      m,
      "LatEleOrderArrayStruct",
      "Fortran struct: lat_ele_order_array_struct"
  );
  auto py_LatEleOrderStruct =
      nb::class_<LatEleOrderStruct>(m, "LatEleOrderStruct", "Fortran struct: lat_ele_order_struct");
  auto py_LatParamStruct =
      nb::class_<LatParamStruct>(m, "LatParamStruct", "Fortran struct: lat_param_struct");
  auto py_LatPointerStruct =
      nb::class_<LatPointerStruct>(m, "LatPointerStruct", "Fortran struct: lat_pointer_struct");
  auto py_LatStruct = nb::class_<LatStruct>(m, "LatStruct", "Fortran struct: lat_struct");
  auto py_Layout = nb::class_<Layout>(m, "Layout", "Fortran struct: layout");
  auto py_LinacNormalModeStruct = nb::class_<LinacNormalModeStruct>(
      m,
      "LinacNormalModeStruct",
      "Fortran struct: linac_normal_mode_struct"
  );
  auto py_LinearEleIsfStruct = nb::class_<LinearEleIsfStruct>(
      m,
      "LinearEleIsfStruct",
      "Fortran struct: linear_ele_isf_struct"
  );
  auto py_LinearIsf1Struct =
      nb::class_<LinearIsf1Struct>(m, "LinearIsf1Struct", "Fortran struct: linear_isf1_struct");
  auto py_MadEnergyStruct =
      nb::class_<MadEnergyStruct>(m, "MadEnergyStruct", "Fortran struct: mad_energy_struct");
  auto py_MadMapStruct =
      nb::class_<MadMapStruct>(m, "MadMapStruct", "Fortran struct: mad_map_struct");
  auto py_MaterialStruct =
      nb::class_<MaterialStruct>(m, "MaterialStruct", "Fortran struct: material_struct");
  auto py_Mesh3DStruct =
      nb::class_<Mesh3dStruct>(m, "Mesh3DStruct", "Fortran struct: mesh3d_struct");
  auto py_Mode3Struct = nb::class_<Mode3Struct>(m, "Mode3Struct", "Fortran struct: mode3_struct");
  auto py_ModeInfoStruct =
      nb::class_<ModeInfoStruct>(m, "ModeInfoStruct", "Fortran struct: mode_info_struct");
  auto py_MolecularComponentStruct = nb::class_<MolecularComponentStruct>(
      m,
      "MolecularComponentStruct",
      "Fortran struct: molecular_component_struct"
  );
  auto py_MomentumApertureStruct = nb::class_<MomentumApertureStruct>(
      m,
      "MomentumApertureStruct",
      "Fortran struct: momentum_aperture_struct"
  );
  auto py_MultipassAllInfoStruct = nb::class_<MultipassAllInfoStruct>(
      m,
      "MultipassAllInfoStruct",
      "Fortran struct: multipass_all_info_struct"
  );
  auto py_MultipassBranchInfoStruct = nb::class_<MultipassBranchInfoStruct>(
      m,
      "MultipassBranchInfoStruct",
      "Fortran struct: multipass_branch_info_struct"
  );
  auto py_MultipassEleInfoStruct = nb::class_<MultipassEleInfoStruct>(
      m,
      "MultipassEleInfoStruct",
      "Fortran struct: multipass_ele_info_struct"
  );
  auto py_MultipassLordInfoStruct = nb::class_<MultipassLordInfoStruct>(
      m,
      "MultipassLordInfoStruct",
      "Fortran struct: multipass_lord_info_struct"
  );
  auto py_MultipassRegionBranchStruct = nb::class_<MultipassRegionBranchStruct>(
      m,
      "MultipassRegionBranchStruct",
      "Fortran struct: multipass_region_branch_struct"
  );
  auto py_MultipassRegionEleStruct = nb::class_<MultipassRegionEleStruct>(
      m,
      "MultipassRegionEleStruct",
      "Fortran struct: multipass_region_ele_struct"
  );
  auto py_MultipassRegionLatStruct = nb::class_<MultipassRegionLatStruct>(
      m,
      "MultipassRegionLatStruct",
      "Fortran struct: multipass_region_lat_struct"
  );
  auto py_MultipoleCacheStruct = nb::class_<MultipoleCacheStruct>(
      m,
      "MultipoleCacheStruct",
      "Fortran struct: multipole_cache_struct"
  );
  auto py_NamedNumberStruct =
      nb::class_<NamedNumberStruct>(m, "NamedNumberStruct", "Fortran struct: named_number_struct");
  auto py_NametableStruct =
      nb::class_<NametableStruct>(m, "NametableStruct", "Fortran struct: nametable_struct");
  auto py_NormalModesStruct =
      nb::class_<NormalModesStruct>(m, "NormalModesStruct", "Fortran struct: normal_modes_struct");
  auto py_OutIoOutputDirectStruct = nb::class_<OutIoOutputDirectStruct>(
      m,
      "OutIoOutputDirectStruct",
      "Fortran struct: out_io_output_direct_struct"
  );
  auto py_ParserControllerStruct = nb::class_<ParserControllerStruct>(
      m,
      "ParserControllerStruct",
      "Fortran struct: parser_controller_struct"
  );
  auto py_ParserEleStruct =
      nb::class_<ParserEleStruct>(m, "ParserEleStruct", "Fortran struct: parser_ele_struct");
  auto py_ParserLatStruct =
      nb::class_<ParserLatStruct>(m, "ParserLatStruct", "Fortran struct: parser_lat_struct");
  auto py_PhotonCoordStruct =
      nb::class_<PhotonCoordStruct>(m, "PhotonCoordStruct", "Fortran struct: photon_coord_struct");
  auto py_PhotonElementStruct = nb::class_<PhotonElementStruct>(
      m,
      "PhotonElementStruct",
      "Fortran struct: photon_element_struct"
  );
  auto py_PhotonInitSplinesStruct = nb::class_<PhotonInitSplinesStruct>(
      m,
      "PhotonInitSplinesStruct",
      "Fortran struct: photon_init_splines_struct"
  );
  auto py_PhotonInitXAngleSplineStruct = nb::class_<PhotonInitXAngleSplineStruct>(
      m,
      "PhotonInitXAngleSplineStruct",
      "Fortran struct: photon_init_x_angle_spline_struct"
  );
  auto py_PhotonInitYAngleSplineStruct = nb::class_<PhotonInitYAngleSplineStruct>(
      m,
      "PhotonInitYAngleSplineStruct",
      "Fortran struct: photon_init_y_angle_spline_struct"
  );
  auto py_PhotonMaterialStruct = nb::class_<PhotonMaterialStruct>(
      m,
      "PhotonMaterialStruct",
      "Fortran struct: photon_material_struct"
  );
  auto py_PhotonReflectSurfaceStruct = nb::class_<PhotonReflectSurfaceStruct>(
      m,
      "PhotonReflectSurfaceStruct",
      "Fortran struct: photon_reflect_surface_struct"
  );
  auto py_PhotonReflectTableStruct = nb::class_<PhotonReflectTableStruct>(
      m,
      "PhotonReflectTableStruct",
      "Fortran struct: photon_reflect_table_struct"
  );
  auto py_PhotonTargetStruct = nb::class_<PhotonTargetStruct>(
      m,
      "PhotonTargetStruct",
      "Fortran struct: photon_target_struct"
  );
  auto py_PhotonTrackStruct =
      nb::class_<PhotonTrackStruct>(m, "PhotonTrackStruct", "Fortran struct: photon_track_struct");
  auto py_PixelDetecStruct =
      nb::class_<PixelDetecStruct>(m, "PixelDetecStruct", "Fortran struct: pixel_detec_struct");
  auto py_PixelPtStruct =
      nb::class_<PixelPtStruct>(m, "PixelPtStruct", "Fortran struct: pixel_pt_struct");
  auto py_PmdHeaderStruct =
      nb::class_<PmdHeaderStruct>(m, "PmdHeaderStruct", "Fortran struct: pmd_header_struct");
  auto py_PreTrackerStruct =
      nb::class_<PreTrackerStruct>(m, "PreTrackerStruct", "Fortran struct: pre_tracker_struct");
  auto py_PtcBranch1Struct =
      nb::class_<PtcBranch1Struct>(m, "PtcBranch1Struct", "Fortran struct: ptc_branch1_struct");
  auto py_PtcLayoutPointerStruct = nb::class_<PtcLayoutPointerStruct>(
      m,
      "PtcLayoutPointerStruct",
      "Fortran struct: ptc_layout_pointer_struct"
  );
  auto py_PtcNormalFormStruct = nb::class_<PtcNormalFormStruct>(
      m,
      "PtcNormalFormStruct",
      "Fortran struct: ptc_normal_form_struct"
  );
  auto py_PtcRadMapStruct =
      nb::class_<PtcRadMapStruct>(m, "PtcRadMapStruct", "Fortran struct: ptc_rad_map_struct");
  auto py_QpAxisStruct =
      nb::class_<QpAxisStruct>(m, "QpAxisStruct", "Fortran struct: qp_axis_struct");
  auto py_QpLegendStruct =
      nb::class_<QpLegendStruct>(m, "QpLegendStruct", "Fortran struct: qp_legend_struct");
  auto py_QpLineStruct =
      nb::class_<QpLineStruct>(m, "QpLineStruct", "Fortran struct: qp_line_struct");
  auto py_QpPointStruct =
      nb::class_<QpPointStruct>(m, "QpPointStruct", "Fortran struct: qp_point_struct");
  auto py_QpRectStruct =
      nb::class_<QpRectStruct>(m, "QpRectStruct", "Fortran struct: qp_rect_struct");
  auto py_QpSymbolStruct =
      nb::class_<QpSymbolStruct>(m, "QpSymbolStruct", "Fortran struct: qp_symbol_struct");
  auto py_RadInt1Struct =
      nb::class_<RadInt1Struct>(m, "RadInt1Struct", "Fortran struct: rad_int1_struct");
  auto py_RadIntAllEleStruct = nb::class_<RadIntAllEleStruct>(
      m,
      "RadIntAllEleStruct",
      "Fortran struct: rad_int_all_ele_struct"
  );
  auto py_RadIntBranchStruct = nb::class_<RadIntBranchStruct>(
      m,
      "RadIntBranchStruct",
      "Fortran struct: rad_int_branch_struct"
  );
  auto py_RadIntCache1Struct = nb::class_<RadIntCache1Struct>(
      m,
      "RadIntCache1Struct",
      "Fortran struct: rad_int_cache1_struct"
  );
  auto py_RadIntInfoStruct =
      nb::class_<RadIntInfoStruct>(m, "RadIntInfoStruct", "Fortran struct: rad_int_info_struct");
  auto py_RadIntTrackPointStruct = nb::class_<RadIntTrackPointStruct>(
      m,
      "RadIntTrackPointStruct",
      "Fortran struct: rad_int_track_point_struct"
  );
  auto py_RadMapEleStruct =
      nb::class_<RadMapEleStruct>(m, "RadMapEleStruct", "Fortran struct: rad_map_ele_struct");
  auto py_RadMapStruct =
      nb::class_<RadMapStruct>(m, "RadMapStruct", "Fortran struct: rad_map_struct");
  auto py_RamperLordStruct =
      nb::class_<RamperLordStruct>(m, "RamperLordStruct", "Fortran struct: ramper_lord_struct");
  auto py_RandomStateStruct =
      nb::class_<RandomStateStruct>(m, "RandomStateStruct", "Fortran struct: random_state_struct");
  auto py_ResonanceHStruct =
      nb::class_<ResonanceHStruct>(m, "ResonanceHStruct", "Fortran struct: resonance_h_struct");
  auto py_RfEleStruct = nb::class_<RfEleStruct>(m, "RfEleStruct", "Fortran struct: rf_ele_struct");
  auto py_RfStairStepStruct =
      nb::class_<RfStairStepStruct>(m, "RfStairStepStruct", "Fortran struct: rf_stair_step_struct");
  auto py_SeqEleStruct =
      nb::class_<SeqEleStruct>(m, "SeqEleStruct", "Fortran struct: seq_ele_struct");
  auto py_SeqStruct = nb::class_<SeqStruct>(m, "SeqStruct", "Fortran struct: seq_struct");
  auto py_SpaceChargeCommonStruct = nb::class_<SpaceChargeCommonStruct>(
      m,
      "SpaceChargeCommonStruct",
      "Fortran struct: space_charge_common_struct"
  );
  auto py_SpinAxisStruct =
      nb::class_<SpinAxisStruct>(m, "SpinAxisStruct", "Fortran struct: spin_axis_struct");
  auto py_SpinEigenStruct =
      nb::class_<SpinEigenStruct>(m, "SpinEigenStruct", "Fortran struct: spin_eigen_struct");
  auto py_SpinMatchingStruct = nb::class_<SpinMatchingStruct>(
      m,
      "SpinMatchingStruct",
      "Fortran struct: spin_matching_struct"
  );
  auto py_SpinOrbitMap1Struct = nb::class_<SpinOrbitMap1Struct>(
      m,
      "SpinOrbitMap1Struct",
      "Fortran struct: spin_orbit_map1_struct"
  );
  auto py_SpinPolarStruct =
      nb::class_<SpinPolarStruct>(m, "SpinPolarStruct", "Fortran struct: spin_polar_struct");
  auto py_SplineStruct =
      nb::class_<SplineStruct>(m, "SplineStruct", "Fortran struct: spline_struct");
  auto py_StrIndexStruct =
      nb::class_<StrIndexStruct>(m, "StrIndexStruct", "Fortran struct: str_index_struct");
  auto py_StrongBeamStruct =
      nb::class_<StrongBeamStruct>(m, "StrongBeamStruct", "Fortran struct: strong_beam_struct");
  auto py_SummationRdtStruct = nb::class_<SummationRdtStruct>(
      m,
      "SummationRdtStruct",
      "Fortran struct: summation_rdt_struct"
  );
  auto py_SurfaceCurvatureStruct = nb::class_<SurfaceCurvatureStruct>(
      m,
      "SurfaceCurvatureStruct",
      "Fortran struct: surface_curvature_struct"
  );
  auto py_SurfaceDisplacementPtStruct = nb::class_<SurfaceDisplacementPtStruct>(
      m,
      "SurfaceDisplacementPtStruct",
      "Fortran struct: surface_displacement_pt_struct"
  );
  auto py_SurfaceDisplacementStruct = nb::class_<SurfaceDisplacementStruct>(
      m,
      "SurfaceDisplacementStruct",
      "Fortran struct: surface_displacement_struct"
  );
  auto py_SurfaceHMisalignPtStruct = nb::class_<SurfaceHMisalignPtStruct>(
      m,
      "SurfaceHMisalignPtStruct",
      "Fortran struct: surface_h_misalign_pt_struct"
  );
  auto py_SurfaceHMisalignStruct = nb::class_<SurfaceHMisalignStruct>(
      m,
      "SurfaceHMisalignStruct",
      "Fortran struct: surface_h_misalign_struct"
  );
  auto py_SurfaceSegmentedPtStruct = nb::class_<SurfaceSegmentedPtStruct>(
      m,
      "SurfaceSegmentedPtStruct",
      "Fortran struct: surface_segmented_pt_struct"
  );
  auto py_SurfaceSegmentedStruct = nb::class_<SurfaceSegmentedStruct>(
      m,
      "SurfaceSegmentedStruct",
      "Fortran struct: surface_segmented_struct"
  );
  auto py_TaoAliasStruct =
      nb::class_<TaoAliasStruct>(m, "TaoAliasStruct", "Fortran struct: tao_alias_struct");
  auto py_TaoBeamBranchStruct = nb::class_<TaoBeamBranchStruct>(
      m,
      "TaoBeamBranchStruct",
      "Fortran struct: tao_beam_branch_struct"
  );
  auto py_TaoBeamUniStruct =
      nb::class_<TaoBeamUniStruct>(m, "TaoBeamUniStruct", "Fortran struct: tao_beam_uni_struct");
  auto py_TaoBuildingWallOrientationStruct = nb::class_<TaoBuildingWallOrientationStruct>(
      m,
      "TaoBuildingWallOrientationStruct",
      "Fortran struct: tao_building_wall_orientation_struct"
  );
  auto py_TaoBuildingWallPointStruct = nb::class_<TaoBuildingWallPointStruct>(
      m,
      "TaoBuildingWallPointStruct",
      "Fortran struct: tao_building_wall_point_struct"
  );
  auto py_TaoBuildingWallSectionStruct = nb::class_<TaoBuildingWallSectionStruct>(
      m,
      "TaoBuildingWallSectionStruct",
      "Fortran struct: tao_building_wall_section_struct"
  );
  auto py_TaoBuildingWallStruct = nb::class_<TaoBuildingWallStruct>(
      m,
      "TaoBuildingWallStruct",
      "Fortran struct: tao_building_wall_struct"
  );
  auto py_TaoCmdHistoryStruct = nb::class_<TaoCmdHistoryStruct>(
      m,
      "TaoCmdHistoryStruct",
      "Fortran struct: tao_cmd_history_struct"
  );
  auto py_TaoCommandFileStruct = nb::class_<TaoCommandFileStruct>(
      m,
      "TaoCommandFileStruct",
      "Fortran struct: tao_command_file_struct"
  );
  auto py_TaoCommonStruct =
      nb::class_<TaoCommonStruct>(m, "TaoCommonStruct", "Fortran struct: tao_common_struct");
  auto py_TaoCurveArrayStruct = nb::class_<TaoCurveArrayStruct>(
      m,
      "TaoCurveArrayStruct",
      "Fortran struct: tao_curve_array_struct"
  );
  auto py_TaoCurveColorStruct = nb::class_<TaoCurveColorStruct>(
      m,
      "TaoCurveColorStruct",
      "Fortran struct: tao_curve_color_struct"
  );
  auto py_TaoCurveOrbitStruct = nb::class_<TaoCurveOrbitStruct>(
      m,
      "TaoCurveOrbitStruct",
      "Fortran struct: tao_curve_orbit_struct"
  );
  auto py_TaoCurveStruct =
      nb::class_<TaoCurveStruct>(m, "TaoCurveStruct", "Fortran struct: tao_curve_struct");
  auto py_TaoD1DataArrayStruct = nb::class_<TaoD1DataArrayStruct>(
      m,
      "TaoD1DataArrayStruct",
      "Fortran struct: tao_d1_data_array_struct"
  );
  auto py_TaoD1DataStruct =
      nb::class_<TaoD1DataStruct>(m, "TaoD1DataStruct", "Fortran struct: tao_d1_data_struct");
  auto py_TaoD2DataArrayStruct = nb::class_<TaoD2DataArrayStruct>(
      m,
      "TaoD2DataArrayStruct",
      "Fortran struct: tao_d2_data_array_struct"
  );
  auto py_TaoD2DataStruct =
      nb::class_<TaoD2DataStruct>(m, "TaoD2DataStruct", "Fortran struct: tao_d2_data_struct");
  auto py_TaoDataArrayStruct = nb::class_<TaoDataArrayStruct>(
      m,
      "TaoDataArrayStruct",
      "Fortran struct: tao_data_array_struct"
  );
  auto py_TaoDataStruct =
      nb::class_<TaoDataStruct>(m, "TaoDataStruct", "Fortran struct: tao_data_struct");
  auto py_TaoDataVarComponentStruct = nb::class_<TaoDataVarComponentStruct>(
      m,
      "TaoDataVarComponentStruct",
      "Fortran struct: tao_data_var_component_struct"
  );
  auto py_TaoDrawingStruct =
      nb::class_<TaoDrawingStruct>(m, "TaoDrawingStruct", "Fortran struct: tao_drawing_struct");
  auto py_TaoDynamicApertureStruct = nb::class_<TaoDynamicApertureStruct>(
      m,
      "TaoDynamicApertureStruct",
      "Fortran struct: tao_dynamic_aperture_struct"
  );
  auto py_TaoElePointerStruct = nb::class_<TaoElePointerStruct>(
      m,
      "TaoElePointerStruct",
      "Fortran struct: tao_ele_pointer_struct"
  );
  auto py_TaoEleShapeInput =
      nb::class_<TaoEleShapeInput>(m, "TaoEleShapeInput", "Fortran struct: tao_ele_shape_input");
  auto py_TaoEleShapeStruct =
      nb::class_<TaoEleShapeStruct>(m, "TaoEleShapeStruct", "Fortran struct: tao_ele_shape_struct");
  auto py_TaoEvalNodeStruct =
      nb::class_<TaoEvalNodeStruct>(m, "TaoEvalNodeStruct", "Fortran struct: tao_eval_node_struct");
  auto py_TaoExpressionInfoStruct = nb::class_<TaoExpressionInfoStruct>(
      m,
      "TaoExpressionInfoStruct",
      "Fortran struct: tao_expression_info_struct"
  );
  auto py_TaoFloorPlanStruct = nb::class_<TaoFloorPlanStruct>(
      m,
      "TaoFloorPlanStruct",
      "Fortran struct: tao_floor_plan_struct"
  );
  auto py_TaoGlobalStruct =
      nb::class_<TaoGlobalStruct>(m, "TaoGlobalStruct", "Fortran struct: tao_global_struct");
  auto py_TaoGraphArrayStruct = nb::class_<TaoGraphArrayStruct>(
      m,
      "TaoGraphArrayStruct",
      "Fortran struct: tao_graph_array_struct"
  );
  auto py_TaoGraphStruct =
      nb::class_<TaoGraphStruct>(m, "TaoGraphStruct", "Fortran struct: tao_graph_struct");
  auto py_TaoHistogramStruct = nb::class_<TaoHistogramStruct>(
      m,
      "TaoHistogramStruct",
      "Fortran struct: tao_histogram_struct"
  );
  auto py_TaoInitStruct =
      nb::class_<TaoInitStruct>(m, "TaoInitStruct", "Fortran struct: tao_init_struct");
  auto py_TaoIntegerArrayStruct = nb::class_<TaoIntegerArrayStruct>(
      m,
      "TaoIntegerArrayStruct",
      "Fortran struct: tao_integer_array_struct"
  );
  auto py_TaoLatSigmaStruct =
      nb::class_<TaoLatSigmaStruct>(m, "TaoLatSigmaStruct", "Fortran struct: tao_lat_sigma_struct");
  auto py_TaoLatticeBranchStruct = nb::class_<TaoLatticeBranchStruct>(
      m,
      "TaoLatticeBranchStruct",
      "Fortran struct: tao_lattice_branch_struct"
  );
  auto py_TaoLatticeStruct =
      nb::class_<TaoLatticeStruct>(m, "TaoLatticeStruct", "Fortran struct: tao_lattice_struct");
  auto py_TaoLogicalArrayStruct = nb::class_<TaoLogicalArrayStruct>(
      m,
      "TaoLogicalArrayStruct",
      "Fortran struct: tao_logical_array_struct"
  );
  auto py_TaoModelBranchStruct = nb::class_<TaoModelBranchStruct>(
      m,
      "TaoModelBranchStruct",
      "Fortran struct: tao_model_branch_struct"
  );
  auto py_TaoModelElementStruct = nb::class_<TaoModelElementStruct>(
      m,
      "TaoModelElementStruct",
      "Fortran struct: tao_model_element_struct"
  );
  auto py_TaoPingScaleStruct = nb::class_<TaoPingScaleStruct>(
      m,
      "TaoPingScaleStruct",
      "Fortran struct: tao_ping_scale_struct"
  );
  auto py_TaoPlotArrayStruct = nb::class_<TaoPlotArrayStruct>(
      m,
      "TaoPlotArrayStruct",
      "Fortran struct: tao_plot_array_struct"
  );
  auto py_TaoPlotCacheStruct = nb::class_<TaoPlotCacheStruct>(
      m,
      "TaoPlotCacheStruct",
      "Fortran struct: tao_plot_cache_struct"
  );
  auto py_TaoPlotPageInput =
      nb::class_<TaoPlotPageInput>(m, "TaoPlotPageInput", "Fortran struct: tao_plot_page_input");
  auto py_TaoPlotPageStruct =
      nb::class_<TaoPlotPageStruct>(m, "TaoPlotPageStruct", "Fortran struct: tao_plot_page_struct");
  auto py_TaoPlotRegionStruct = nb::class_<TaoPlotRegionStruct>(
      m,
      "TaoPlotRegionStruct",
      "Fortran struct: tao_plot_region_struct"
  );
  auto py_TaoPlotStruct =
      nb::class_<TaoPlotStruct>(m, "TaoPlotStruct", "Fortran struct: tao_plot_struct");
  auto py_TaoRealPointerStruct = nb::class_<TaoRealPointerStruct>(
      m,
      "TaoRealPointerStruct",
      "Fortran struct: tao_real_pointer_struct"
  );
  auto py_TaoShapePatternPointStruct = nb::class_<TaoShapePatternPointStruct>(
      m,
      "TaoShapePatternPointStruct",
      "Fortran struct: tao_shape_pattern_point_struct"
  );
  auto py_TaoShapePatternStruct = nb::class_<TaoShapePatternStruct>(
      m,
      "TaoShapePatternStruct",
      "Fortran struct: tao_shape_pattern_struct"
  );
  auto py_TaoSpinDnDpzStruct = nb::class_<TaoSpinDnDpzStruct>(
      m,
      "TaoSpinDnDpzStruct",
      "Fortran struct: tao_spin_dn_dpz_struct"
  );
  auto py_TaoSpinEleStruct =
      nb::class_<TaoSpinEleStruct>(m, "TaoSpinEleStruct", "Fortran struct: tao_spin_ele_struct");
  auto py_TaoSpinMapStruct =
      nb::class_<TaoSpinMapStruct>(m, "TaoSpinMapStruct", "Fortran struct: tao_spin_map_struct");
  auto py_TaoSpinPolarizationStruct = nb::class_<TaoSpinPolarizationStruct>(
      m,
      "TaoSpinPolarizationStruct",
      "Fortran struct: tao_spin_polarization_struct"
  );
  auto py_TaoStringArrayStruct = nb::class_<TaoStringArrayStruct>(
      m,
      "TaoStringArrayStruct",
      "Fortran struct: tao_string_array_struct"
  );
  auto py_TaoSuperUniverseStruct = nb::class_<TaoSuperUniverseStruct>(
      m,
      "TaoSuperUniverseStruct",
      "Fortran struct: tao_super_universe_struct"
  );
  auto py_TaoTitleStruct =
      nb::class_<TaoTitleStruct>(m, "TaoTitleStruct", "Fortran struct: tao_title_struct");
  auto py_TaoTop10Struct =
      nb::class_<TaoTop10Struct>(m, "TaoTop10Struct", "Fortran struct: tao_top10_struct");
  auto py_TaoUniverseCalcStruct = nb::class_<TaoUniverseCalcStruct>(
      m,
      "TaoUniverseCalcStruct",
      "Fortran struct: tao_universe_calc_struct"
  );
  auto py_TaoUniversePointerStruct = nb::class_<TaoUniversePointerStruct>(
      m,
      "TaoUniversePointerStruct",
      "Fortran struct: tao_universe_pointer_struct"
  );
  auto py_TaoUniverseStruct =
      nb::class_<TaoUniverseStruct>(m, "TaoUniverseStruct", "Fortran struct: tao_universe_struct");
  auto py_TaoV1VarArrayStruct = nb::class_<TaoV1VarArrayStruct>(
      m,
      "TaoV1VarArrayStruct",
      "Fortran struct: tao_v1_var_array_struct"
  );
  auto py_TaoV1VarStruct =
      nb::class_<TaoV1VarStruct>(m, "TaoV1VarStruct", "Fortran struct: tao_v1_var_struct");
  auto py_TaoVarArrayStruct =
      nb::class_<TaoVarArrayStruct>(m, "TaoVarArrayStruct", "Fortran struct: tao_var_array_struct");
  auto py_TaoVarSlaveStruct =
      nb::class_<TaoVarSlaveStruct>(m, "TaoVarSlaveStruct", "Fortran struct: tao_var_slave_struct");
  auto py_TaoVarStruct =
      nb::class_<TaoVarStruct>(m, "TaoVarStruct", "Fortran struct: tao_var_struct");
  auto py_TaoWaveKickPtStruct = nb::class_<TaoWaveKickPtStruct>(
      m,
      "TaoWaveKickPtStruct",
      "Fortran struct: tao_wave_kick_pt_struct"
  );
  auto py_TaoWaveStruct =
      nb::class_<TaoWaveStruct>(m, "TaoWaveStruct", "Fortran struct: tao_wave_struct");
  auto py_TargetPointStruct =
      nb::class_<TargetPointStruct>(m, "TargetPointStruct", "Fortran struct: target_point_struct");
  auto py_TaylorStruct =
      nb::class_<TaylorStruct>(m, "TaylorStruct", "Fortran struct: taylor_struct");
  auto py_TaylorTermStruct =
      nb::class_<TaylorTermStruct>(m, "TaylorTermStruct", "Fortran struct: taylor_term_struct");
  auto py_TestSubStruct =
      nb::class_<TestSubStruct>(m, "TestSubStruct", "Fortran struct: test_sub_struct");
  auto py_TestSubSubStruct =
      nb::class_<TestSubSubStruct>(m, "TestSubSubStruct", "Fortran struct: test_sub_sub_struct");
  auto py_TrackPointStruct =
      nb::class_<TrackPointStruct>(m, "TrackPointStruct", "Fortran struct: track_point_struct");
  auto py_TrackStruct = nb::class_<TrackStruct>(m, "TrackStruct", "Fortran struct: track_struct");
  auto py_TricubicCmplxCoefStruct = nb::class_<TricubicCmplxCoefStruct>(
      m,
      "TricubicCmplxCoefStruct",
      "Fortran struct: tricubic_cmplx_coef_struct"
  );
  auto py_TricubicCoefStruct = nb::class_<TricubicCoefStruct>(
      m,
      "TricubicCoefStruct",
      "Fortran struct: tricubic_coef_struct"
  );
  auto py_TwissStruct = nb::class_<TwissStruct>(m, "TwissStruct", "Fortran struct: twiss_struct");
  auto py_VarLengthStringStruct = nb::class_<VarLengthStringStruct>(
      m,
      "VarLengthStringStruct",
      "Fortran struct: var_length_string_struct"
  );
  auto py_WakeLrModeStruct =
      nb::class_<WakeLrModeStruct>(m, "WakeLrModeStruct", "Fortran struct: wake_lr_mode_struct");
  auto py_WakeLrStruct =
      nb::class_<WakeLrStruct>(m, "WakeLrStruct", "Fortran struct: wake_lr_struct");
  auto py_WakeSrModeStruct =
      nb::class_<WakeSrModeStruct>(m, "WakeSrModeStruct", "Fortran struct: wake_sr_mode_struct");
  auto py_WakeSrStruct =
      nb::class_<WakeSrStruct>(m, "WakeSrStruct", "Fortran struct: wake_sr_struct");
  auto py_WakeSrZLongStruct = nb::class_<WakeSrZLongStruct>(
      m,
      "WakeSrZLongStruct",
      "Fortran struct: wake_sr_z_long_struct"
  );
  auto py_WakeStruct = nb::class_<WakeStruct>(m, "WakeStruct", "Fortran struct: wake_struct");
  auto py_Wall3DSectionStruct = nb::class_<Wall3dSectionStruct>(
      m,
      "Wall3DSectionStruct",
      "Fortran struct: wall3d_section_struct"
  );
  auto py_Wall3DStruct =
      nb::class_<Wall3dStruct>(m, "Wall3DStruct", "Fortran struct: wall3d_struct");
  auto py_Wall3DVertexStruct = nb::class_<Wall3dVertexStruct>(
      m,
      "Wall3DVertexStruct",
      "Fortran struct: wall3d_vertex_struct"
  );
  auto py_XyDispStruct =
      nb::class_<XyDispStruct>(m, "XyDispStruct", "Fortran struct: xy_disp_struct");
  init_ac_kicker_freq_struct(m, py_AcKickerFreqStruct);
  init_ac_kicker_struct(m, py_AcKickerStruct);
  init_ac_kicker_time_struct(m, py_AcKickerTimeStruct);
  init_all_encompassing_struct(m, py_AllEncompassingStruct);
  init_all_pointer_struct(m, py_AllPointerStruct);
  init_anormal_mode_struct(m, py_AnormalModeStruct);
  init_aperture_param_struct(m, py_ApertureParamStruct);
  init_aperture_point_struct(m, py_AperturePointStruct);
  init_aperture_scan_struct(m, py_ApertureScanStruct);
  init_astra_lattice_param_struct(m, py_AstraLatticeParamStruct);
  init_base_line_ele_struct(m, py_BaseLineEleStruct);
  init_bbu_beam_struct(m, py_BbuBeamStruct);
  init_bbu_param_struct(m, py_BbuParamStruct);
  init_bbu_stage_struct(m, py_BbuStageStruct);
  init_beam_init_struct(m, py_BeamInitStruct);
  init_beam_struct(m, py_BeamStruct);
  init_bicubic_cmplx_coef_struct(m, py_BicubicCmplxCoefStruct);
  init_bicubic_coef_struct(m, py_BicubicCoefStruct);
  init_bin_struct(m, py_BinStruct);
  init_bmad_common_struct(m, py_BmadCommonStruct);
  init_bmad_normal_form_struct(m, py_BmadNormalFormStruct);
  init_bookkeeping_state_struct(m, py_BookkeepingStateStruct);
  init_bpm_phase_coupling_struct(m, py_BpmPhaseCouplingStruct);
  init_branch_pointer_struct(m, py_BranchPointerStruct);
  init_branch_struct(m, py_BranchStruct);
  init_bunch_params_struct(m, py_BunchParamsStruct);
  init_bunch_struct(m, py_BunchStruct);
  init_bunch_track_struct(m, py_BunchTrackStruct);
  init_cartesian_map_struct(m, py_CartesianMapStruct);
  init_cartesian_map_term1_struct(m, py_CartesianMapTerm1Struct);
  init_cartesian_map_term_struct(m, py_CartesianMapTermStruct);
  init_cmplx_field1_at_2D_pt_struct(m, py_CmplxField1At2DPtStruct);
  init_cmplx_field1_at_3D_pt_struct(m, py_CmplxField1At3DPtStruct);
  init_cmplx_field_at_2D_box_struct(m, py_CmplxFieldAt2DBoxStruct);
  init_cmplx_field_at_3D_box_struct(m, py_CmplxFieldAt3DBoxStruct);
  init_complex_taylor_struct(m, py_ComplexTaylorStruct);
  init_complex_taylor_term_struct(m, py_ComplexTaylorTermStruct);
  init_control_ramp1_struct(m, py_ControlRamp1Struct);
  init_control_struct(m, py_ControlStruct);
  init_control_var1_struct(m, py_ControlVar1Struct);
  init_controller_struct(m, py_ControllerStruct);
  init_converter_dir_1D_struct(m, py_ConverterDir1DStruct);
  init_converter_dir_2D_struct(m, py_ConverterDir2DStruct);
  init_converter_dir_coef_struct(m, py_ConverterDirCoefStruct);
  init_converter_direction_out_struct(m, py_ConverterDirectionOutStruct);
  init_converter_distribution_struct(m, py_ConverterDistributionStruct);
  init_converter_prob_pc_r_struct(m, py_ConverterProbPcRStruct);
  init_converter_struct(m, py_ConverterStruct);
  init_converter_sub_distribution_struct(m, py_ConverterSubDistributionStruct);
  init_coord_array_struct(m, py_CoordArrayStruct);
  init_coord_struct(m, py_CoordStruct);
  init_crystal_param_struct(m, py_CrystalParamStruct);
  init_csr_bunch_slice_struct(m, py_CsrBunchSliceStruct);
  init_csr_ele_info_struct(m, py_CsrEleInfoStruct);
  init_csr_kick1_struct(m, py_CsrKick1Struct);
  init_csr_particle_position_struct(m, py_CsrParticlePositionStruct);
  init_csr_struct(m, py_CsrStruct);
  init_cylindrical_map_struct(m, py_CylindricalMapStruct);
  init_cylindrical_map_term1_struct(m, py_CylindricalMapTerm1Struct);
  init_cylindrical_map_term_struct(m, py_CylindricalMapTermStruct);
  init_diffuse_param_struct(m, py_DiffuseParamStruct);
  init_do_loop_struct(m, py_DoLoopStruct);
  init_ele_attribute_struct(m, py_EleAttributeStruct);
  init_ele_pointer_struct(m, py_ElePointerStruct);
  init_ele_struct(m, py_EleStruct);
  init_ellipse_beam_init_struct(m, py_EllipseBeamInitStruct);
  init_em_field_struct(m, py_EmFieldStruct);
  init_expression_atom_struct(m, py_ExpressionAtomStruct);
  init_expression_tree_struct(m, py_ExpressionTreeStruct);
  init_extra_parsing_info_struct(m, py_ExtraParsingInfoStruct);
  init_fibre(m, py_Fibre);
  init_field1_at_2D_pt_struct(m, py_Field1At2DPtStruct);
  init_field1_at_3D_pt_struct(m, py_Field1At3DPtStruct);
  init_field_at_2D_box_struct(m, py_FieldAt2DBoxStruct);
  init_field_at_3D_box_struct(m, py_FieldAt3DBoxStruct);
  init_floor_position_struct(m, py_FloorPositionStruct);
  init_foil_struct(m, py_FoilStruct);
  init_fringe_field_info_struct(m, py_FringeFieldInfoStruct);
  init_gen_grad_curve_struct(m, py_GenGradCurveStruct);
  init_gen_gradients_struct(m, py_GenGradientsStruct);
  init_general_bin_struct(m, py_GeneralBinStruct);
  init_gg_taylor_struct(m, py_GgTaylorStruct);
  init_gg_taylor_term_struct(m, py_GgTaylorTermStruct);
  init_gpt_lat_param_struct(m, py_GptLatParamStruct);
  init_grid_beam_init_struct(m, py_GridBeamInitStruct);
  init_grid_field_pt1_struct(m, py_GridFieldPt1Struct);
  init_grid_field_pt_struct(m, py_GridFieldPtStruct);
  init_grid_field_struct(m, py_GridFieldStruct);
  init_high_energy_space_charge_struct(m, py_HighEnergySpaceChargeStruct);
  init_ibs_lifetime_struct(m, py_IbsLifetimeStruct);
  init_ibs_maxratio_struct(m, py_IbsMaxratioStruct);
  init_ibs_sim_param_struct(m, py_IbsSimParamStruct);
  init_ibs_struct(m, py_IbsStruct);
  init_interval1_coef_struct(m, py_Interval1CoefStruct);
  init_kv_beam_init_struct(m, py_KvBeamInitStruct);
  init_lat_ele_loc_struct(m, py_LatEleLocStruct);
  init_lat_ele_order1_struct(m, py_LatEleOrder1Struct);
  init_lat_ele_order_array_struct(m, py_LatEleOrderArrayStruct);
  init_lat_ele_order_struct(m, py_LatEleOrderStruct);
  init_lat_param_struct(m, py_LatParamStruct);
  init_lat_pointer_struct(m, py_LatPointerStruct);
  init_lat_struct(m, py_LatStruct);
  init_layout(m, py_Layout);
  init_linac_normal_mode_struct(m, py_LinacNormalModeStruct);
  init_linear_ele_isf_struct(m, py_LinearEleIsfStruct);
  init_linear_isf1_struct(m, py_LinearIsf1Struct);
  init_mad_energy_struct(m, py_MadEnergyStruct);
  init_mad_map_struct(m, py_MadMapStruct);
  init_material_struct(m, py_MaterialStruct);
  init_mesh3d_struct(m, py_Mesh3DStruct);
  init_mode3_struct(m, py_Mode3Struct);
  init_mode_info_struct(m, py_ModeInfoStruct);
  init_molecular_component_struct(m, py_MolecularComponentStruct);
  init_momentum_aperture_struct(m, py_MomentumApertureStruct);
  init_multipass_all_info_struct(m, py_MultipassAllInfoStruct);
  init_multipass_branch_info_struct(m, py_MultipassBranchInfoStruct);
  init_multipass_ele_info_struct(m, py_MultipassEleInfoStruct);
  init_multipass_lord_info_struct(m, py_MultipassLordInfoStruct);
  init_multipass_region_branch_struct(m, py_MultipassRegionBranchStruct);
  init_multipass_region_ele_struct(m, py_MultipassRegionEleStruct);
  init_multipass_region_lat_struct(m, py_MultipassRegionLatStruct);
  init_multipole_cache_struct(m, py_MultipoleCacheStruct);
  init_named_number_struct(m, py_NamedNumberStruct);
  init_nametable_struct(m, py_NametableStruct);
  init_normal_modes_struct(m, py_NormalModesStruct);
  init_out_io_output_direct_struct(m, py_OutIoOutputDirectStruct);
  init_parser_controller_struct(m, py_ParserControllerStruct);
  init_parser_ele_struct(m, py_ParserEleStruct);
  init_parser_lat_struct(m, py_ParserLatStruct);
  init_photon_coord_struct(m, py_PhotonCoordStruct);
  init_photon_element_struct(m, py_PhotonElementStruct);
  init_photon_init_splines_struct(m, py_PhotonInitSplinesStruct);
  init_photon_init_x_angle_spline_struct(m, py_PhotonInitXAngleSplineStruct);
  init_photon_init_y_angle_spline_struct(m, py_PhotonInitYAngleSplineStruct);
  init_photon_material_struct(m, py_PhotonMaterialStruct);
  init_photon_reflect_surface_struct(m, py_PhotonReflectSurfaceStruct);
  init_photon_reflect_table_struct(m, py_PhotonReflectTableStruct);
  init_photon_target_struct(m, py_PhotonTargetStruct);
  init_photon_track_struct(m, py_PhotonTrackStruct);
  init_pixel_detec_struct(m, py_PixelDetecStruct);
  init_pixel_pt_struct(m, py_PixelPtStruct);
  init_pmd_header_struct(m, py_PmdHeaderStruct);
  init_pre_tracker_struct(m, py_PreTrackerStruct);
  init_ptc_branch1_struct(m, py_PtcBranch1Struct);
  init_ptc_layout_pointer_struct(m, py_PtcLayoutPointerStruct);
  init_ptc_normal_form_struct(m, py_PtcNormalFormStruct);
  init_ptc_rad_map_struct(m, py_PtcRadMapStruct);
  init_qp_axis_struct(m, py_QpAxisStruct);
  init_qp_legend_struct(m, py_QpLegendStruct);
  init_qp_line_struct(m, py_QpLineStruct);
  init_qp_point_struct(m, py_QpPointStruct);
  init_qp_rect_struct(m, py_QpRectStruct);
  init_qp_symbol_struct(m, py_QpSymbolStruct);
  init_rad_int1_struct(m, py_RadInt1Struct);
  init_rad_int_all_ele_struct(m, py_RadIntAllEleStruct);
  init_rad_int_branch_struct(m, py_RadIntBranchStruct);
  init_rad_int_cache1_struct(m, py_RadIntCache1Struct);
  init_rad_int_info_struct(m, py_RadIntInfoStruct);
  init_rad_int_track_point_struct(m, py_RadIntTrackPointStruct);
  init_rad_map_ele_struct(m, py_RadMapEleStruct);
  init_rad_map_struct(m, py_RadMapStruct);
  init_ramper_lord_struct(m, py_RamperLordStruct);
  init_random_state_struct(m, py_RandomStateStruct);
  init_resonance_h_struct(m, py_ResonanceHStruct);
  init_rf_ele_struct(m, py_RfEleStruct);
  init_rf_stair_step_struct(m, py_RfStairStepStruct);
  init_seq_ele_struct(m, py_SeqEleStruct);
  init_seq_struct(m, py_SeqStruct);
  init_space_charge_common_struct(m, py_SpaceChargeCommonStruct);
  init_spin_axis_struct(m, py_SpinAxisStruct);
  init_spin_eigen_struct(m, py_SpinEigenStruct);
  init_spin_matching_struct(m, py_SpinMatchingStruct);
  init_spin_orbit_map1_struct(m, py_SpinOrbitMap1Struct);
  init_spin_polar_struct(m, py_SpinPolarStruct);
  init_spline_struct(m, py_SplineStruct);
  init_str_index_struct(m, py_StrIndexStruct);
  init_strong_beam_struct(m, py_StrongBeamStruct);
  init_summation_rdt_struct(m, py_SummationRdtStruct);
  init_surface_curvature_struct(m, py_SurfaceCurvatureStruct);
  init_surface_displacement_pt_struct(m, py_SurfaceDisplacementPtStruct);
  init_surface_displacement_struct(m, py_SurfaceDisplacementStruct);
  init_surface_h_misalign_pt_struct(m, py_SurfaceHMisalignPtStruct);
  init_surface_h_misalign_struct(m, py_SurfaceHMisalignStruct);
  init_surface_segmented_pt_struct(m, py_SurfaceSegmentedPtStruct);
  init_surface_segmented_struct(m, py_SurfaceSegmentedStruct);
  init_tao_alias_struct(m, py_TaoAliasStruct);
  init_tao_beam_branch_struct(m, py_TaoBeamBranchStruct);
  init_tao_beam_uni_struct(m, py_TaoBeamUniStruct);
  init_tao_building_wall_orientation_struct(m, py_TaoBuildingWallOrientationStruct);
  init_tao_building_wall_point_struct(m, py_TaoBuildingWallPointStruct);
  init_tao_building_wall_section_struct(m, py_TaoBuildingWallSectionStruct);
  init_tao_building_wall_struct(m, py_TaoBuildingWallStruct);
  init_tao_cmd_history_struct(m, py_TaoCmdHistoryStruct);
  init_tao_command_file_struct(m, py_TaoCommandFileStruct);
  init_tao_common_struct(m, py_TaoCommonStruct);
  init_tao_curve_array_struct(m, py_TaoCurveArrayStruct);
  init_tao_curve_color_struct(m, py_TaoCurveColorStruct);
  init_tao_curve_orbit_struct(m, py_TaoCurveOrbitStruct);
  init_tao_curve_struct(m, py_TaoCurveStruct);
  init_tao_d1_data_array_struct(m, py_TaoD1DataArrayStruct);
  init_tao_d1_data_struct(m, py_TaoD1DataStruct);
  init_tao_d2_data_array_struct(m, py_TaoD2DataArrayStruct);
  init_tao_d2_data_struct(m, py_TaoD2DataStruct);
  init_tao_data_array_struct(m, py_TaoDataArrayStruct);
  init_tao_data_struct(m, py_TaoDataStruct);
  init_tao_data_var_component_struct(m, py_TaoDataVarComponentStruct);
  init_tao_drawing_struct(m, py_TaoDrawingStruct);
  init_tao_dynamic_aperture_struct(m, py_TaoDynamicApertureStruct);
  init_tao_ele_pointer_struct(m, py_TaoElePointerStruct);
  init_tao_ele_shape_input(m, py_TaoEleShapeInput);
  init_tao_ele_shape_struct(m, py_TaoEleShapeStruct);
  init_tao_eval_node_struct(m, py_TaoEvalNodeStruct);
  init_tao_expression_info_struct(m, py_TaoExpressionInfoStruct);
  init_tao_floor_plan_struct(m, py_TaoFloorPlanStruct);
  init_tao_global_struct(m, py_TaoGlobalStruct);
  init_tao_graph_array_struct(m, py_TaoGraphArrayStruct);
  init_tao_graph_struct(m, py_TaoGraphStruct);
  init_tao_histogram_struct(m, py_TaoHistogramStruct);
  init_tao_init_struct(m, py_TaoInitStruct);
  init_tao_integer_array_struct(m, py_TaoIntegerArrayStruct);
  init_tao_lat_sigma_struct(m, py_TaoLatSigmaStruct);
  init_tao_lattice_branch_struct(m, py_TaoLatticeBranchStruct);
  init_tao_lattice_struct(m, py_TaoLatticeStruct);
  init_tao_logical_array_struct(m, py_TaoLogicalArrayStruct);
  init_tao_model_branch_struct(m, py_TaoModelBranchStruct);
  init_tao_model_element_struct(m, py_TaoModelElementStruct);
  init_tao_ping_scale_struct(m, py_TaoPingScaleStruct);
  init_tao_plot_array_struct(m, py_TaoPlotArrayStruct);
  init_tao_plot_cache_struct(m, py_TaoPlotCacheStruct);
  init_tao_plot_page_input(m, py_TaoPlotPageInput);
  init_tao_plot_page_struct(m, py_TaoPlotPageStruct);
  init_tao_plot_region_struct(m, py_TaoPlotRegionStruct);
  init_tao_plot_struct(m, py_TaoPlotStruct);
  init_tao_real_pointer_struct(m, py_TaoRealPointerStruct);
  init_tao_shape_pattern_point_struct(m, py_TaoShapePatternPointStruct);
  init_tao_shape_pattern_struct(m, py_TaoShapePatternStruct);
  init_tao_spin_dn_dpz_struct(m, py_TaoSpinDnDpzStruct);
  init_tao_spin_ele_struct(m, py_TaoSpinEleStruct);
  init_tao_spin_map_struct(m, py_TaoSpinMapStruct);
  init_tao_spin_polarization_struct(m, py_TaoSpinPolarizationStruct);
  init_tao_string_array_struct(m, py_TaoStringArrayStruct);
  init_tao_super_universe_struct(m, py_TaoSuperUniverseStruct);
  init_tao_title_struct(m, py_TaoTitleStruct);
  init_tao_top10_struct(m, py_TaoTop10Struct);
  init_tao_universe_calc_struct(m, py_TaoUniverseCalcStruct);
  init_tao_universe_pointer_struct(m, py_TaoUniversePointerStruct);
  init_tao_universe_struct(m, py_TaoUniverseStruct);
  init_tao_v1_var_array_struct(m, py_TaoV1VarArrayStruct);
  init_tao_v1_var_struct(m, py_TaoV1VarStruct);
  init_tao_var_array_struct(m, py_TaoVarArrayStruct);
  init_tao_var_slave_struct(m, py_TaoVarSlaveStruct);
  init_tao_var_struct(m, py_TaoVarStruct);
  init_tao_wave_kick_pt_struct(m, py_TaoWaveKickPtStruct);
  init_tao_wave_struct(m, py_TaoWaveStruct);
  init_target_point_struct(m, py_TargetPointStruct);
  init_taylor_struct(m, py_TaylorStruct);
  init_taylor_term_struct(m, py_TaylorTermStruct);
  init_test_sub_struct(m, py_TestSubStruct);
  init_test_sub_sub_struct(m, py_TestSubSubStruct);
  init_track_point_struct(m, py_TrackPointStruct);
  init_track_struct(m, py_TrackStruct);
  init_tricubic_cmplx_coef_struct(m, py_TricubicCmplxCoefStruct);
  init_tricubic_coef_struct(m, py_TricubicCoefStruct);
  init_twiss_struct(m, py_TwissStruct);
  init_var_length_string_struct(m, py_VarLengthStringStruct);
  init_wake_lr_mode_struct(m, py_WakeLrModeStruct);
  init_wake_lr_struct(m, py_WakeLrStruct);
  init_wake_sr_mode_struct(m, py_WakeSrModeStruct);
  init_wake_sr_struct(m, py_WakeSrStruct);
  init_wake_sr_z_long_struct(m, py_WakeSrZLongStruct);
  init_wake_struct(m, py_WakeStruct);
  init_wall3d_section_struct(m, py_Wall3DSectionStruct);
  init_wall3d_struct(m, py_Wall3DStruct);
  init_wall3d_vertex_struct(m, py_Wall3DVertexStruct);
  init_xy_disp_struct(m, py_XyDispStruct);

  // Hand-written bindings
  init_common_structs(m);

  // Routine submodules (one per C++ namespace)
  auto m_bmad = m.def_submodule("bmad", "Bmad routines");
  auto m_extra = m.def_submodule("extra", "CppBmadExtra routines");
  auto m_test = m.def_submodule("test", "CppBmadTest routines");
  auto m_simutils = m.def_submodule("simutils", "SimUtils routines");
  auto m_tao = m.def_submodule("tao", "Tao routines");
  auto m_bsim = m.def_submodule("bsim", "bsim routines");

  // Hand-written submodule bindings
  init_tao_hooks(m_tao);
  init_bmad_hooks(m_bmad);

  // Routine initializers
  init_Bmad_routines_a(m_bmad);
  init_Bmad_routines_b(m_bmad);
  init_Bmad_routines_c(m_bmad);
  init_Bmad_routines_d(m_bmad);
  init_Bmad_routines_e(m_bmad);
  init_Bmad_routines_f(m_bmad);
  init_Bmad_routines_g(m_bmad);
  init_Bmad_routines_h(m_bmad);
  init_Bmad_routines_i(m_bmad);
  init_Bmad_routines_k(m_bmad);
  init_Bmad_routines_l(m_bmad);
  init_Bmad_routines_m(m_bmad);
  init_Bmad_routines_n(m_bmad);
  init_Bmad_routines_o(m_bmad);
  init_Bmad_routines_p(m_bmad);
  init_Bmad_routines_r(m_bmad);
  init_Bmad_routines_s(m_bmad);
  init_Bmad_routines_t(m_bmad);
  init_Bmad_routines_u(m_bmad);
  init_Bmad_routines_v(m_bmad);
  init_Bmad_routines_w(m_bmad);
  init_Bmad_routines_x(m_bmad);
  init_Bmad_routines_y(m_bmad);
  init_Bmad_routines_z(m_bmad);
  init_CppBmadExtra_routines_s(m_extra);
  init_CppBmadTest_routines_t(m_test);
  init_SimUtils_routines_a(m_simutils);
  init_SimUtils_routines_b(m_simutils);
  init_SimUtils_routines_c(m_simutils);
  init_SimUtils_routines_d(m_simutils);
  init_SimUtils_routines_e(m_simutils);
  init_SimUtils_routines_f(m_simutils);
  init_SimUtils_routines_g(m_simutils);
  init_SimUtils_routines_h(m_simutils);
  init_SimUtils_routines_i(m_simutils);
  init_SimUtils_routines_j(m_simutils);
  init_SimUtils_routines_l(m_simutils);
  init_SimUtils_routines_m(m_simutils);
  init_SimUtils_routines_n(m_simutils);
  init_SimUtils_routines_o(m_simutils);
  init_SimUtils_routines_p(m_simutils);
  init_SimUtils_routines_q(m_simutils);
  init_SimUtils_routines_r(m_simutils);
  init_SimUtils_routines_s(m_simutils);
  init_SimUtils_routines_t(m_simutils);
  init_SimUtils_routines_u(m_simutils);
  init_SimUtils_routines_v(m_simutils);
  init_SimUtils_routines_w(m_simutils);
  init_SimUtils_routines_x(m_simutils);
  init_SimUtils_routines_z(m_simutils);
  init_Tao_routines_d(m_tao);
  init_Tao_routines_i(m_tao);
  init_Tao_routines_t(m_tao);
  init_bsim_routines_b(m_bsim);
  init_bsim_routines_c(m_bsim);
  init_bsim_routines_h(m_bsim);
  init_bsim_routines_i(m_bsim);
  init_bsim_routines_l(m_bsim);
  init_bsim_routines_r(m_bsim);
  init_bsim_routines_s(m_bsim);
  init_bsim_routines_w(m_bsim);
}
