#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <string>

#include "bmad/json.hpp"
#include "pybmad/common_structs.hpp"
#include "pybmad/generated/init.hpp"
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
  auto py_SplineStruct =
      nb::class_<SplineStruct>(m, "SplineStruct", "Fortran struct: spline_struct");
  auto py_SpinPolarStruct =
      nb::class_<SpinPolarStruct>(m, "SpinPolarStruct", "Fortran struct: spin_polar_struct");
  auto py_AcKickerTimeStruct = nb::class_<AcKickerTimeStruct>(
      m,
      "AcKickerTimeStruct",
      "Fortran struct: ac_kicker_time_struct"
  );
  auto py_AcKickerFreqStruct = nb::class_<AcKickerFreqStruct>(
      m,
      "AcKickerFreqStruct",
      "Fortran struct: ac_kicker_freq_struct"
  );
  auto py_AcKickerStruct =
      nb::class_<AcKickerStruct>(m, "AcKickerStruct", "Fortran struct: ac_kicker_struct");
  auto py_Interval1CoefStruct = nb::class_<Interval1CoefStruct>(
      m,
      "Interval1CoefStruct",
      "Fortran struct: interval1_coef_struct"
  );
  auto py_PhotonReflectTableStruct = nb::class_<PhotonReflectTableStruct>(
      m,
      "PhotonReflectTableStruct",
      "Fortran struct: photon_reflect_table_struct"
  );
  auto py_PhotonReflectSurfaceStruct = nb::class_<PhotonReflectSurfaceStruct>(
      m,
      "PhotonReflectSurfaceStruct",
      "Fortran struct: photon_reflect_surface_struct"
  );
  auto py_CoordStruct = nb::class_<CoordStruct>(m, "CoordStruct", "Fortran struct: coord_struct");
  auto py_CoordArrayStruct =
      nb::class_<CoordArrayStruct>(m, "CoordArrayStruct", "Fortran struct: coord_array_struct");
  auto py_BpmPhaseCouplingStruct = nb::class_<BpmPhaseCouplingStruct>(
      m,
      "BpmPhaseCouplingStruct",
      "Fortran struct: bpm_phase_coupling_struct"
  );
  auto py_ExpressionAtomStruct = nb::class_<ExpressionAtomStruct>(
      m,
      "ExpressionAtomStruct",
      "Fortran struct: expression_atom_struct"
  );
  auto py_WakeSrZLongStruct = nb::class_<WakeSrZLongStruct>(
      m,
      "WakeSrZLongStruct",
      "Fortran struct: wake_sr_z_long_struct"
  );
  auto py_WakeSrModeStruct =
      nb::class_<WakeSrModeStruct>(m, "WakeSrModeStruct", "Fortran struct: wake_sr_mode_struct");
  auto py_WakeSrStruct =
      nb::class_<WakeSrStruct>(m, "WakeSrStruct", "Fortran struct: wake_sr_struct");
  auto py_WakeLrModeStruct =
      nb::class_<WakeLrModeStruct>(m, "WakeLrModeStruct", "Fortran struct: wake_lr_mode_struct");
  auto py_WakeLrStruct =
      nb::class_<WakeLrStruct>(m, "WakeLrStruct", "Fortran struct: wake_lr_struct");
  auto py_LatEleLocStruct =
      nb::class_<LatEleLocStruct>(m, "LatEleLocStruct", "Fortran struct: lat_ele_loc_struct");
  auto py_WakeStruct = nb::class_<WakeStruct>(m, "WakeStruct", "Fortran struct: wake_struct");
  auto py_TaylorTermStruct =
      nb::class_<TaylorTermStruct>(m, "TaylorTermStruct", "Fortran struct: taylor_term_struct");
  auto py_TaylorStruct =
      nb::class_<TaylorStruct>(m, "TaylorStruct", "Fortran struct: taylor_struct");
  auto py_GgTaylorTermStruct = nb::class_<GgTaylorTermStruct>(
      m,
      "GgTaylorTermStruct",
      "Fortran struct: gg_taylor_term_struct"
  );
  auto py_GgTaylorStruct =
      nb::class_<GgTaylorStruct>(m, "GgTaylorStruct", "Fortran struct: gg_taylor_struct");
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
  auto py_CartesianMapStruct = nb::class_<CartesianMapStruct>(
      m,
      "CartesianMapStruct",
      "Fortran struct: cartesian_map_struct"
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
  auto py_CylindricalMapStruct = nb::class_<CylindricalMapStruct>(
      m,
      "CylindricalMapStruct",
      "Fortran struct: cylindrical_map_struct"
  );
  auto py_BicubicCmplxCoefStruct = nb::class_<BicubicCmplxCoefStruct>(
      m,
      "BicubicCmplxCoefStruct",
      "Fortran struct: bicubic_cmplx_coef_struct"
  );
  auto py_TricubicCmplxCoefStruct = nb::class_<TricubicCmplxCoefStruct>(
      m,
      "TricubicCmplxCoefStruct",
      "Fortran struct: tricubic_cmplx_coef_struct"
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
  auto py_FloorPositionStruct = nb::class_<FloorPositionStruct>(
      m,
      "FloorPositionStruct",
      "Fortran struct: floor_position_struct"
  );
  auto py_HighEnergySpaceChargeStruct = nb::class_<HighEnergySpaceChargeStruct>(
      m,
      "HighEnergySpaceChargeStruct",
      "Fortran struct: high_energy_space_charge_struct"
  );
  auto py_XyDispStruct =
      nb::class_<XyDispStruct>(m, "XyDispStruct", "Fortran struct: xy_disp_struct");
  auto py_TwissStruct = nb::class_<TwissStruct>(m, "TwissStruct", "Fortran struct: twiss_struct");
  auto py_Mode3Struct = nb::class_<Mode3Struct>(m, "Mode3Struct", "Fortran struct: mode3_struct");
  auto py_BookkeepingStateStruct = nb::class_<BookkeepingStateStruct>(
      m,
      "BookkeepingStateStruct",
      "Fortran struct: bookkeeping_state_struct"
  );
  auto py_RadMapStruct =
      nb::class_<RadMapStruct>(m, "RadMapStruct", "Fortran struct: rad_map_struct");
  auto py_RadMapEleStruct =
      nb::class_<RadMapEleStruct>(m, "RadMapEleStruct", "Fortran struct: rad_map_ele_struct");
  auto py_GenGrad1Struct =
      nb::class_<GenGrad1Struct>(m, "GenGrad1Struct", "Fortran struct: gen_grad1_struct");
  auto py_GenGradMapStruct =
      nb::class_<GenGradMapStruct>(m, "GenGradMapStruct", "Fortran struct: gen_grad_map_struct");
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
  auto py_TargetPointStruct =
      nb::class_<TargetPointStruct>(m, "TargetPointStruct", "Fortran struct: target_point_struct");
  auto py_SurfaceCurvatureStruct = nb::class_<SurfaceCurvatureStruct>(
      m,
      "SurfaceCurvatureStruct",
      "Fortran struct: surface_curvature_struct"
  );
  auto py_PhotonTargetStruct = nb::class_<PhotonTargetStruct>(
      m,
      "PhotonTargetStruct",
      "Fortran struct: photon_target_struct"
  );
  auto py_PhotonMaterialStruct = nb::class_<PhotonMaterialStruct>(
      m,
      "PhotonMaterialStruct",
      "Fortran struct: photon_material_struct"
  );
  auto py_PixelPtStruct =
      nb::class_<PixelPtStruct>(m, "PixelPtStruct", "Fortran struct: pixel_pt_struct");
  auto py_PixelDetecStruct =
      nb::class_<PixelDetecStruct>(m, "PixelDetecStruct", "Fortran struct: pixel_detec_struct");
  auto py_PhotonElementStruct = nb::class_<PhotonElementStruct>(
      m,
      "PhotonElementStruct",
      "Fortran struct: photon_element_struct"
  );
  auto py_Wall3DVertexStruct = nb::class_<Wall3dVertexStruct>(
      m,
      "Wall3DVertexStruct",
      "Fortran struct: wall3d_vertex_struct"
  );
  auto py_Wall3DSectionStruct = nb::class_<Wall3dSectionStruct>(
      m,
      "Wall3DSectionStruct",
      "Fortran struct: wall3d_section_struct"
  );
  auto py_Wall3DStruct =
      nb::class_<Wall3dStruct>(m, "Wall3DStruct", "Fortran struct: wall3d_struct");
  auto py_RamperLordStruct =
      nb::class_<RamperLordStruct>(m, "RamperLordStruct", "Fortran struct: ramper_lord_struct");
  auto py_ControlStruct =
      nb::class_<ControlStruct>(m, "ControlStruct", "Fortran struct: control_struct");
  auto py_ControlVar1Struct =
      nb::class_<ControlVar1Struct>(m, "ControlVar1Struct", "Fortran struct: control_var1_struct");
  auto py_ControlRamp1Struct = nb::class_<ControlRamp1Struct>(
      m,
      "ControlRamp1Struct",
      "Fortran struct: control_ramp1_struct"
  );
  auto py_ControllerStruct =
      nb::class_<ControllerStruct>(m, "ControllerStruct", "Fortran struct: controller_struct");
  auto py_EllipseBeamInitStruct = nb::class_<EllipseBeamInitStruct>(
      m,
      "EllipseBeamInitStruct",
      "Fortran struct: ellipse_beam_init_struct"
  );
  auto py_KvBeamInitStruct =
      nb::class_<KvBeamInitStruct>(m, "KvBeamInitStruct", "Fortran struct: kv_beam_init_struct");
  auto py_GridBeamInitStruct = nb::class_<GridBeamInitStruct>(
      m,
      "GridBeamInitStruct",
      "Fortran struct: grid_beam_init_struct"
  );
  auto py_BeamInitStruct =
      nb::class_<BeamInitStruct>(m, "BeamInitStruct", "Fortran struct: beam_init_struct");
  auto py_LatParamStruct =
      nb::class_<LatParamStruct>(m, "LatParamStruct", "Fortran struct: lat_param_struct");
  auto py_ModeInfoStruct =
      nb::class_<ModeInfoStruct>(m, "ModeInfoStruct", "Fortran struct: mode_info_struct");
  auto py_PreTrackerStruct =
      nb::class_<PreTrackerStruct>(m, "PreTrackerStruct", "Fortran struct: pre_tracker_struct");
  auto py_AnormalModeStruct =
      nb::class_<AnormalModeStruct>(m, "AnormalModeStruct", "Fortran struct: anormal_mode_struct");
  auto py_LinacNormalModeStruct = nb::class_<LinacNormalModeStruct>(
      m,
      "LinacNormalModeStruct",
      "Fortran struct: linac_normal_mode_struct"
  );
  auto py_NormalModesStruct =
      nb::class_<NormalModesStruct>(m, "NormalModesStruct", "Fortran struct: normal_modes_struct");
  auto py_EmFieldStruct =
      nb::class_<EmFieldStruct>(m, "EmFieldStruct", "Fortran struct: em_field_struct");
  auto py_StrongBeamStruct =
      nb::class_<StrongBeamStruct>(m, "StrongBeamStruct", "Fortran struct: strong_beam_struct");
  auto py_TrackPointStruct =
      nb::class_<TrackPointStruct>(m, "TrackPointStruct", "Fortran struct: track_point_struct");
  auto py_TrackStruct = nb::class_<TrackStruct>(m, "TrackStruct", "Fortran struct: track_struct");
  auto py_SpaceChargeCommonStruct = nb::class_<SpaceChargeCommonStruct>(
      m,
      "SpaceChargeCommonStruct",
      "Fortran struct: space_charge_common_struct"
  );
  auto py_BmadCommonStruct =
      nb::class_<BmadCommonStruct>(m, "BmadCommonStruct", "Fortran struct: bmad_common_struct");
  auto py_RadInt1Struct =
      nb::class_<RadInt1Struct>(m, "RadInt1Struct", "Fortran struct: rad_int1_struct");
  auto py_RadIntBranchStruct = nb::class_<RadIntBranchStruct>(
      m,
      "RadIntBranchStruct",
      "Fortran struct: rad_int_branch_struct"
  );
  auto py_RadIntAllEleStruct = nb::class_<RadIntAllEleStruct>(
      m,
      "RadIntAllEleStruct",
      "Fortran struct: rad_int_all_ele_struct"
  );
  auto py_RfStairStepStruct =
      nb::class_<RfStairStepStruct>(m, "RfStairStepStruct", "Fortran struct: rf_stair_step_struct");
  auto py_RfEleStruct = nb::class_<RfEleStruct>(m, "RfEleStruct", "Fortran struct: rf_ele_struct");
  auto py_EleStruct = nb::class_<EleStruct>(m, "EleStruct", "Fortran struct: ele_struct");
  auto py_ComplexTaylorTermStruct = nb::class_<ComplexTaylorTermStruct>(
      m,
      "ComplexTaylorTermStruct",
      "Fortran struct: complex_taylor_term_struct"
  );
  auto py_ComplexTaylorStruct = nb::class_<ComplexTaylorStruct>(
      m,
      "ComplexTaylorStruct",
      "Fortran struct: complex_taylor_struct"
  );
  auto py_BranchStruct =
      nb::class_<BranchStruct>(m, "BranchStruct", "Fortran struct: branch_struct");
  auto py_LatStruct = nb::class_<LatStruct>(m, "LatStruct", "Fortran struct: lat_struct");
  auto py_BunchStruct = nb::class_<BunchStruct>(m, "BunchStruct", "Fortran struct: bunch_struct");
  auto py_BunchParamsStruct =
      nb::class_<BunchParamsStruct>(m, "BunchParamsStruct", "Fortran struct: bunch_params_struct");
  auto py_BeamStruct = nb::class_<BeamStruct>(m, "BeamStruct", "Fortran struct: beam_struct");
  auto py_AperturePointStruct = nb::class_<AperturePointStruct>(
      m,
      "AperturePointStruct",
      "Fortran struct: aperture_point_struct"
  );
  auto py_ApertureParamStruct = nb::class_<ApertureParamStruct>(
      m,
      "ApertureParamStruct",
      "Fortran struct: aperture_param_struct"
  );
  auto py_ApertureScanStruct = nb::class_<ApertureScanStruct>(
      m,
      "ApertureScanStruct",
      "Fortran struct: aperture_scan_struct"
  );
  auto py_ElePointerStruct =
      nb::class_<ElePointerStruct>(m, "ElePointerStruct", "Fortran struct: ele_pointer_struct");
  auto py_ExpressionTreeStruct = nb::class_<ExpressionTreeStruct>(
      m,
      "ExpressionTreeStruct",
      "Fortran struct: expression_tree_struct"
  );
  auto py_NametableStruct =
      nb::class_<NametableStruct>(m, "NametableStruct", "Fortran struct: nametable_struct");
  auto py_AllPointerStruct =
      nb::class_<AllPointerStruct>(m, "AllPointerStruct", "Fortran struct: all_pointer_struct");
  auto py_ParserEleStruct =
      nb::class_<ParserEleStruct>(m, "ParserEleStruct", "Fortran struct: parser_ele_struct");
  auto py_ParserControllerStruct = nb::class_<ParserControllerStruct>(
      m,
      "ParserControllerStruct",
      "Fortran struct: parser_controller_struct"
  );
  auto py_TaoSpinDnDpzStruct = nb::class_<TaoSpinDnDpzStruct>(
      m,
      "TaoSpinDnDpzStruct",
      "Fortran struct: tao_spin_dn_dpz_struct"
  );
  auto py_ResonanceHStruct =
      nb::class_<ResonanceHStruct>(m, "ResonanceHStruct", "Fortran struct: resonance_h_struct");
  auto py_SpinOrbitMap1Struct = nb::class_<SpinOrbitMap1Struct>(
      m,
      "SpinOrbitMap1Struct",
      "Fortran struct: spin_orbit_map1_struct"
  );
  auto py_SpinAxisStruct =
      nb::class_<SpinAxisStruct>(m, "SpinAxisStruct", "Fortran struct: spin_axis_struct");
  auto py_PtcNormalFormStruct = nb::class_<PtcNormalFormStruct>(
      m,
      "PtcNormalFormStruct",
      "Fortran struct: ptc_normal_form_struct"
  );
  auto py_BmadNormalFormStruct = nb::class_<BmadNormalFormStruct>(
      m,
      "BmadNormalFormStruct",
      "Fortran struct: bmad_normal_form_struct"
  );
  auto py_BunchTrackStruct =
      nb::class_<BunchTrackStruct>(m, "BunchTrackStruct", "Fortran struct: bunch_track_struct");
  auto py_SummationRdtStruct = nb::class_<SummationRdtStruct>(
      m,
      "SummationRdtStruct",
      "Fortran struct: summation_rdt_struct"
  );
  auto py_TaoEleShapeStruct =
      nb::class_<TaoEleShapeStruct>(m, "TaoEleShapeStruct", "Fortran struct: tao_ele_shape_struct");
  auto py_TaoElePointerStruct = nb::class_<TaoElePointerStruct>(
      m,
      "TaoElePointerStruct",
      "Fortran struct: tao_ele_pointer_struct"
  );
  auto py_TaoCurveStruct =
      nb::class_<TaoCurveStruct>(m, "TaoCurveStruct", "Fortran struct: tao_curve_struct");
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
  auto py_TaoHistogramStruct = nb::class_<TaoHistogramStruct>(
      m,
      "TaoHistogramStruct",
      "Fortran struct: tao_histogram_struct"
  );
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
  auto py_TaoLatSigmaStruct =
      nb::class_<TaoLatSigmaStruct>(m, "TaoLatSigmaStruct", "Fortran struct: tao_lat_sigma_struct");
  auto py_TaoSpinEleStruct =
      nb::class_<TaoSpinEleStruct>(m, "TaoSpinEleStruct", "Fortran struct: tao_spin_ele_struct");
  auto py_TaoPlotCacheStruct = nb::class_<TaoPlotCacheStruct>(
      m,
      "TaoPlotCacheStruct",
      "Fortran struct: tao_plot_cache_struct"
  );
  auto py_TaoSpinPolarizationStruct = nb::class_<TaoSpinPolarizationStruct>(
      m,
      "TaoSpinPolarizationStruct",
      "Fortran struct: tao_spin_polarization_struct"
  );
  auto py_TaoLatticeBranchStruct = nb::class_<TaoLatticeBranchStruct>(
      m,
      "TaoLatticeBranchStruct",
      "Fortran struct: tao_lattice_branch_struct"
  );
  auto py_TaoModelElementStruct = nb::class_<TaoModelElementStruct>(
      m,
      "TaoModelElementStruct",
      "Fortran struct: tao_model_element_struct"
  );
  auto py_TaoBeamBranchStruct = nb::class_<TaoBeamBranchStruct>(
      m,
      "TaoBeamBranchStruct",
      "Fortran struct: tao_beam_branch_struct"
  );
  auto py_TaoD1DataStruct =
      nb::class_<TaoD1DataStruct>(m, "TaoD1DataStruct", "Fortran struct: tao_d1_data_struct");
  auto py_TaoD2DataStruct =
      nb::class_<TaoD2DataStruct>(m, "TaoD2DataStruct", "Fortran struct: tao_d2_data_struct");
  auto py_TaoDataVarComponentStruct = nb::class_<TaoDataVarComponentStruct>(
      m,
      "TaoDataVarComponentStruct",
      "Fortran struct: tao_data_var_component_struct"
  );
  auto py_TaoGraphStruct =
      nb::class_<TaoGraphStruct>(m, "TaoGraphStruct", "Fortran struct: tao_graph_struct");
  auto py_TaoPlotStruct =
      nb::class_<TaoPlotStruct>(m, "TaoPlotStruct", "Fortran struct: tao_plot_struct");
  auto py_TaoPlotRegionStruct = nb::class_<TaoPlotRegionStruct>(
      m,
      "TaoPlotRegionStruct",
      "Fortran struct: tao_plot_region_struct"
  );
  auto py_TaoUniversePointerStruct = nb::class_<TaoUniversePointerStruct>(
      m,
      "TaoUniversePointerStruct",
      "Fortran struct: tao_universe_pointer_struct"
  );
  auto py_TaoSuperUniverseStruct = nb::class_<TaoSuperUniverseStruct>(
      m,
      "TaoSuperUniverseStruct",
      "Fortran struct: tao_super_universe_struct"
  );
  auto py_TaoVarStruct =
      nb::class_<TaoVarStruct>(m, "TaoVarStruct", "Fortran struct: tao_var_struct");
  auto py_TaoVarSlaveStruct =
      nb::class_<TaoVarSlaveStruct>(m, "TaoVarSlaveStruct", "Fortran struct: tao_var_slave_struct");
  auto py_TaoLatticeStruct =
      nb::class_<TaoLatticeStruct>(m, "TaoLatticeStruct", "Fortran struct: tao_lattice_struct");
  auto py_TaoBeamUniStruct =
      nb::class_<TaoBeamUniStruct>(m, "TaoBeamUniStruct", "Fortran struct: tao_beam_uni_struct");
  auto py_TaoDynamicApertureStruct = nb::class_<TaoDynamicApertureStruct>(
      m,
      "TaoDynamicApertureStruct",
      "Fortran struct: tao_dynamic_aperture_struct"
  );
  auto py_TaoModelBranchStruct = nb::class_<TaoModelBranchStruct>(
      m,
      "TaoModelBranchStruct",
      "Fortran struct: tao_model_branch_struct"
  );
  auto py_TaoSpinMapStruct =
      nb::class_<TaoSpinMapStruct>(m, "TaoSpinMapStruct", "Fortran struct: tao_spin_map_struct");
  auto py_TaoDataStruct =
      nb::class_<TaoDataStruct>(m, "TaoDataStruct", "Fortran struct: tao_data_struct");
  auto py_TaoPingScaleStruct = nb::class_<TaoPingScaleStruct>(
      m,
      "TaoPingScaleStruct",
      "Fortran struct: tao_ping_scale_struct"
  );
  auto py_TaoUniverseCalcStruct = nb::class_<TaoUniverseCalcStruct>(
      m,
      "TaoUniverseCalcStruct",
      "Fortran struct: tao_universe_calc_struct"
  );
  auto py_LatEleOrderStruct =
      nb::class_<LatEleOrderStruct>(m, "LatEleOrderStruct", "Fortran struct: lat_ele_order_struct");
  auto py_TaoExpressionInfoStruct = nb::class_<TaoExpressionInfoStruct>(
      m,
      "TaoExpressionInfoStruct",
      "Fortran struct: tao_expression_info_struct"
  );
  auto py_TaoEvalNodeStruct =
      nb::class_<TaoEvalNodeStruct>(m, "TaoEvalNodeStruct", "Fortran struct: tao_eval_node_struct");
  auto py_TaoTitleStruct =
      nb::class_<TaoTitleStruct>(m, "TaoTitleStruct", "Fortran struct: tao_title_struct");
  auto py_QpRectStruct =
      nb::class_<QpRectStruct>(m, "QpRectStruct", "Fortran struct: qp_rect_struct");
  auto py_TaoDrawingStruct =
      nb::class_<TaoDrawingStruct>(m, "TaoDrawingStruct", "Fortran struct: tao_drawing_struct");
  auto py_TaoShapePatternStruct = nb::class_<TaoShapePatternStruct>(
      m,
      "TaoShapePatternStruct",
      "Fortran struct: tao_shape_pattern_struct"
  );
  auto py_TaoShapePatternPointStruct = nb::class_<TaoShapePatternPointStruct>(
      m,
      "TaoShapePatternPointStruct",
      "Fortran struct: tao_shape_pattern_point_struct"
  );
  auto py_QpAxisStruct =
      nb::class_<QpAxisStruct>(m, "QpAxisStruct", "Fortran struct: qp_axis_struct");
  auto py_QpLegendStruct =
      nb::class_<QpLegendStruct>(m, "QpLegendStruct", "Fortran struct: qp_legend_struct");
  auto py_QpPointStruct =
      nb::class_<QpPointStruct>(m, "QpPointStruct", "Fortran struct: qp_point_struct");
  auto py_QpLineStruct =
      nb::class_<QpLineStruct>(m, "QpLineStruct", "Fortran struct: qp_line_struct");
  auto py_QpSymbolStruct =
      nb::class_<QpSymbolStruct>(m, "QpSymbolStruct", "Fortran struct: qp_symbol_struct");
  auto py_TaoFloorPlanStruct = nb::class_<TaoFloorPlanStruct>(
      m,
      "TaoFloorPlanStruct",
      "Fortran struct: tao_floor_plan_struct"
  );
  auto py_TaoV1VarStruct =
      nb::class_<TaoV1VarStruct>(m, "TaoV1VarStruct", "Fortran struct: tao_v1_var_struct");
  auto py_TaoGlobalStruct =
      nb::class_<TaoGlobalStruct>(m, "TaoGlobalStruct", "Fortran struct: tao_global_struct");
  auto py_TaoInitStruct =
      nb::class_<TaoInitStruct>(m, "TaoInitStruct", "Fortran struct: tao_init_struct");
  auto py_TaoCommonStruct =
      nb::class_<TaoCommonStruct>(m, "TaoCommonStruct", "Fortran struct: tao_common_struct");
  auto py_TaoPlotPageStruct =
      nb::class_<TaoPlotPageStruct>(m, "TaoPlotPageStruct", "Fortran struct: tao_plot_page_struct");
  auto py_TaoBuildingWallStruct = nb::class_<TaoBuildingWallStruct>(
      m,
      "TaoBuildingWallStruct",
      "Fortran struct: tao_building_wall_struct"
  );
  auto py_TaoBuildingWallOrientationStruct = nb::class_<TaoBuildingWallOrientationStruct>(
      m,
      "TaoBuildingWallOrientationStruct",
      "Fortran struct: tao_building_wall_orientation_struct"
  );
  auto py_TaoBuildingWallSectionStruct = nb::class_<TaoBuildingWallSectionStruct>(
      m,
      "TaoBuildingWallSectionStruct",
      "Fortran struct: tao_building_wall_section_struct"
  );
  auto py_TaoBuildingWallPointStruct = nb::class_<TaoBuildingWallPointStruct>(
      m,
      "TaoBuildingWallPointStruct",
      "Fortran struct: tao_building_wall_point_struct"
  );
  auto py_TaoWaveStruct =
      nb::class_<TaoWaveStruct>(m, "TaoWaveStruct", "Fortran struct: tao_wave_struct");
  auto py_TaoWaveKickPtStruct = nb::class_<TaoWaveKickPtStruct>(
      m,
      "TaoWaveKickPtStruct",
      "Fortran struct: tao_wave_kick_pt_struct"
  );
  auto py_TaoCmdHistoryStruct = nb::class_<TaoCmdHistoryStruct>(
      m,
      "TaoCmdHistoryStruct",
      "Fortran struct: tao_cmd_history_struct"
  );
  auto py_TaoUniverseStruct =
      nb::class_<TaoUniverseStruct>(m, "TaoUniverseStruct", "Fortran struct: tao_universe_struct");
  auto py_MadEnergyStruct =
      nb::class_<MadEnergyStruct>(m, "MadEnergyStruct", "Fortran struct: mad_energy_struct");
  auto py_MadMapStruct =
      nb::class_<MadMapStruct>(m, "MadMapStruct", "Fortran struct: mad_map_struct");
  auto py_RandomStateStruct =
      nb::class_<RandomStateStruct>(m, "RandomStateStruct", "Fortran struct: random_state_struct");
  auto py_BbuStageStruct =
      nb::class_<BbuStageStruct>(m, "BbuStageStruct", "Fortran struct: bbu_stage_struct");
  auto py_BbuBeamStruct =
      nb::class_<BbuBeamStruct>(m, "BbuBeamStruct", "Fortran struct: bbu_beam_struct");
  auto py_BbuParamStruct =
      nb::class_<BbuParamStruct>(m, "BbuParamStruct", "Fortran struct: bbu_param_struct");
  auto py_Fibre = nb::class_<Fibre>(m, "Fibre", "Fortran struct: fibre");
  auto py_Layout = nb::class_<Layout>(m, "Layout", "Fortran struct: layout");
  auto py_AllEncompassingStruct = nb::class_<AllEncompassingStruct>(
      m,
      "AllEncompassingStruct",
      "Fortran struct: all_encompassing_struct"
  );
  auto py_TestSubStruct =
      nb::class_<TestSubStruct>(m, "TestSubStruct", "Fortran struct: test_sub_struct");
  auto py_TestSubSubStruct =
      nb::class_<TestSubSubStruct>(m, "TestSubSubStruct", "Fortran struct: test_sub_sub_struct");
  init_spline_struct(m, py_SplineStruct);
  init_spin_polar_struct(m, py_SpinPolarStruct);
  init_ac_kicker_time_struct(m, py_AcKickerTimeStruct);
  init_ac_kicker_freq_struct(m, py_AcKickerFreqStruct);
  init_ac_kicker_struct(m, py_AcKickerStruct);
  init_interval1_coef_struct(m, py_Interval1CoefStruct);
  init_photon_reflect_table_struct(m, py_PhotonReflectTableStruct);
  init_photon_reflect_surface_struct(m, py_PhotonReflectSurfaceStruct);
  init_coord_struct(m, py_CoordStruct);
  init_coord_array_struct(m, py_CoordArrayStruct);
  init_bpm_phase_coupling_struct(m, py_BpmPhaseCouplingStruct);
  init_expression_atom_struct(m, py_ExpressionAtomStruct);
  init_wake_sr_z_long_struct(m, py_WakeSrZLongStruct);
  init_wake_sr_mode_struct(m, py_WakeSrModeStruct);
  init_wake_sr_struct(m, py_WakeSrStruct);
  init_wake_lr_mode_struct(m, py_WakeLrModeStruct);
  init_wake_lr_struct(m, py_WakeLrStruct);
  init_lat_ele_loc_struct(m, py_LatEleLocStruct);
  init_wake_struct(m, py_WakeStruct);
  init_taylor_term_struct(m, py_TaylorTermStruct);
  init_taylor_struct(m, py_TaylorStruct);
  init_gg_taylor_term_struct(m, py_GgTaylorTermStruct);
  init_gg_taylor_struct(m, py_GgTaylorStruct);
  init_cartesian_map_term1_struct(m, py_CartesianMapTerm1Struct);
  init_cartesian_map_term_struct(m, py_CartesianMapTermStruct);
  init_cartesian_map_struct(m, py_CartesianMapStruct);
  init_cylindrical_map_term1_struct(m, py_CylindricalMapTerm1Struct);
  init_cylindrical_map_term_struct(m, py_CylindricalMapTermStruct);
  init_cylindrical_map_struct(m, py_CylindricalMapStruct);
  init_bicubic_cmplx_coef_struct(m, py_BicubicCmplxCoefStruct);
  init_tricubic_cmplx_coef_struct(m, py_TricubicCmplxCoefStruct);
  init_grid_field_pt1_struct(m, py_GridFieldPt1Struct);
  init_grid_field_pt_struct(m, py_GridFieldPtStruct);
  init_grid_field_struct(m, py_GridFieldStruct);
  init_floor_position_struct(m, py_FloorPositionStruct);
  init_high_energy_space_charge_struct(m, py_HighEnergySpaceChargeStruct);
  init_xy_disp_struct(m, py_XyDispStruct);
  init_twiss_struct(m, py_TwissStruct);
  init_mode3_struct(m, py_Mode3Struct);
  init_bookkeeping_state_struct(m, py_BookkeepingStateStruct);
  init_rad_map_struct(m, py_RadMapStruct);
  init_rad_map_ele_struct(m, py_RadMapEleStruct);
  init_gen_grad1_struct(m, py_GenGrad1Struct);
  init_gen_grad_map_struct(m, py_GenGradMapStruct);
  init_surface_segmented_pt_struct(m, py_SurfaceSegmentedPtStruct);
  init_surface_segmented_struct(m, py_SurfaceSegmentedStruct);
  init_surface_h_misalign_pt_struct(m, py_SurfaceHMisalignPtStruct);
  init_surface_h_misalign_struct(m, py_SurfaceHMisalignStruct);
  init_surface_displacement_pt_struct(m, py_SurfaceDisplacementPtStruct);
  init_surface_displacement_struct(m, py_SurfaceDisplacementStruct);
  init_target_point_struct(m, py_TargetPointStruct);
  init_surface_curvature_struct(m, py_SurfaceCurvatureStruct);
  init_photon_target_struct(m, py_PhotonTargetStruct);
  init_photon_material_struct(m, py_PhotonMaterialStruct);
  init_pixel_pt_struct(m, py_PixelPtStruct);
  init_pixel_detec_struct(m, py_PixelDetecStruct);
  init_photon_element_struct(m, py_PhotonElementStruct);
  init_wall3d_vertex_struct(m, py_Wall3DVertexStruct);
  init_wall3d_section_struct(m, py_Wall3DSectionStruct);
  init_wall3d_struct(m, py_Wall3DStruct);
  init_ramper_lord_struct(m, py_RamperLordStruct);
  init_control_struct(m, py_ControlStruct);
  init_control_var1_struct(m, py_ControlVar1Struct);
  init_control_ramp1_struct(m, py_ControlRamp1Struct);
  init_controller_struct(m, py_ControllerStruct);
  init_ellipse_beam_init_struct(m, py_EllipseBeamInitStruct);
  init_kv_beam_init_struct(m, py_KvBeamInitStruct);
  init_grid_beam_init_struct(m, py_GridBeamInitStruct);
  init_beam_init_struct(m, py_BeamInitStruct);
  init_lat_param_struct(m, py_LatParamStruct);
  init_mode_info_struct(m, py_ModeInfoStruct);
  init_pre_tracker_struct(m, py_PreTrackerStruct);
  init_anormal_mode_struct(m, py_AnormalModeStruct);
  init_linac_normal_mode_struct(m, py_LinacNormalModeStruct);
  init_normal_modes_struct(m, py_NormalModesStruct);
  init_em_field_struct(m, py_EmFieldStruct);
  init_strong_beam_struct(m, py_StrongBeamStruct);
  init_track_point_struct(m, py_TrackPointStruct);
  init_track_struct(m, py_TrackStruct);
  init_space_charge_common_struct(m, py_SpaceChargeCommonStruct);
  init_bmad_common_struct(m, py_BmadCommonStruct);
  init_rad_int1_struct(m, py_RadInt1Struct);
  init_rad_int_branch_struct(m, py_RadIntBranchStruct);
  init_rad_int_all_ele_struct(m, py_RadIntAllEleStruct);
  init_rf_stair_step_struct(m, py_RfStairStepStruct);
  init_rf_ele_struct(m, py_RfEleStruct);
  init_ele_struct(m, py_EleStruct);
  init_complex_taylor_term_struct(m, py_ComplexTaylorTermStruct);
  init_complex_taylor_struct(m, py_ComplexTaylorStruct);
  init_branch_struct(m, py_BranchStruct);
  init_lat_struct(m, py_LatStruct);
  init_bunch_struct(m, py_BunchStruct);
  init_bunch_params_struct(m, py_BunchParamsStruct);
  init_beam_struct(m, py_BeamStruct);
  init_aperture_point_struct(m, py_AperturePointStruct);
  init_aperture_param_struct(m, py_ApertureParamStruct);
  init_aperture_scan_struct(m, py_ApertureScanStruct);
  init_ele_pointer_struct(m, py_ElePointerStruct);
  init_expression_tree_struct(m, py_ExpressionTreeStruct);
  init_nametable_struct(m, py_NametableStruct);
  init_all_pointer_struct(m, py_AllPointerStruct);
  init_parser_ele_struct(m, py_ParserEleStruct);
  init_parser_controller_struct(m, py_ParserControllerStruct);
  init_tao_spin_dn_dpz_struct(m, py_TaoSpinDnDpzStruct);
  init_resonance_h_struct(m, py_ResonanceHStruct);
  init_spin_orbit_map1_struct(m, py_SpinOrbitMap1Struct);
  init_spin_axis_struct(m, py_SpinAxisStruct);
  init_ptc_normal_form_struct(m, py_PtcNormalFormStruct);
  init_bmad_normal_form_struct(m, py_BmadNormalFormStruct);
  init_bunch_track_struct(m, py_BunchTrackStruct);
  init_summation_rdt_struct(m, py_SummationRdtStruct);
  init_tao_ele_shape_struct(m, py_TaoEleShapeStruct);
  init_tao_ele_pointer_struct(m, py_TaoElePointerStruct);
  init_tao_curve_struct(m, py_TaoCurveStruct);
  init_tao_curve_color_struct(m, py_TaoCurveColorStruct);
  init_tao_curve_orbit_struct(m, py_TaoCurveOrbitStruct);
  init_tao_histogram_struct(m, py_TaoHistogramStruct);
  init_lat_ele_order1_struct(m, py_LatEleOrder1Struct);
  init_lat_ele_order_array_struct(m, py_LatEleOrderArrayStruct);
  init_tao_lat_sigma_struct(m, py_TaoLatSigmaStruct);
  init_tao_spin_ele_struct(m, py_TaoSpinEleStruct);
  init_tao_plot_cache_struct(m, py_TaoPlotCacheStruct);
  init_tao_spin_polarization_struct(m, py_TaoSpinPolarizationStruct);
  init_tao_lattice_branch_struct(m, py_TaoLatticeBranchStruct);
  init_tao_model_element_struct(m, py_TaoModelElementStruct);
  init_tao_beam_branch_struct(m, py_TaoBeamBranchStruct);
  init_tao_d1_data_struct(m, py_TaoD1DataStruct);
  init_tao_d2_data_struct(m, py_TaoD2DataStruct);
  init_tao_data_var_component_struct(m, py_TaoDataVarComponentStruct);
  init_tao_graph_struct(m, py_TaoGraphStruct);
  init_tao_plot_struct(m, py_TaoPlotStruct);
  init_tao_plot_region_struct(m, py_TaoPlotRegionStruct);
  init_tao_universe_pointer_struct(m, py_TaoUniversePointerStruct);
  init_tao_super_universe_struct(m, py_TaoSuperUniverseStruct);
  init_tao_var_struct(m, py_TaoVarStruct);
  init_tao_var_slave_struct(m, py_TaoVarSlaveStruct);
  init_tao_lattice_struct(m, py_TaoLatticeStruct);
  init_tao_beam_uni_struct(m, py_TaoBeamUniStruct);
  init_tao_dynamic_aperture_struct(m, py_TaoDynamicApertureStruct);
  init_tao_model_branch_struct(m, py_TaoModelBranchStruct);
  init_tao_spin_map_struct(m, py_TaoSpinMapStruct);
  init_tao_data_struct(m, py_TaoDataStruct);
  init_tao_ping_scale_struct(m, py_TaoPingScaleStruct);
  init_tao_universe_calc_struct(m, py_TaoUniverseCalcStruct);
  init_lat_ele_order_struct(m, py_LatEleOrderStruct);
  init_tao_expression_info_struct(m, py_TaoExpressionInfoStruct);
  init_tao_eval_node_struct(m, py_TaoEvalNodeStruct);
  init_tao_title_struct(m, py_TaoTitleStruct);
  init_qp_rect_struct(m, py_QpRectStruct);
  init_tao_drawing_struct(m, py_TaoDrawingStruct);
  init_tao_shape_pattern_struct(m, py_TaoShapePatternStruct);
  init_tao_shape_pattern_point_struct(m, py_TaoShapePatternPointStruct);
  init_qp_axis_struct(m, py_QpAxisStruct);
  init_qp_legend_struct(m, py_QpLegendStruct);
  init_qp_point_struct(m, py_QpPointStruct);
  init_qp_line_struct(m, py_QpLineStruct);
  init_qp_symbol_struct(m, py_QpSymbolStruct);
  init_tao_floor_plan_struct(m, py_TaoFloorPlanStruct);
  init_tao_v1_var_struct(m, py_TaoV1VarStruct);
  init_tao_global_struct(m, py_TaoGlobalStruct);
  init_tao_init_struct(m, py_TaoInitStruct);
  init_tao_common_struct(m, py_TaoCommonStruct);
  init_tao_plot_page_struct(m, py_TaoPlotPageStruct);
  init_tao_building_wall_struct(m, py_TaoBuildingWallStruct);
  init_tao_building_wall_orientation_struct(m, py_TaoBuildingWallOrientationStruct);
  init_tao_building_wall_section_struct(m, py_TaoBuildingWallSectionStruct);
  init_tao_building_wall_point_struct(m, py_TaoBuildingWallPointStruct);
  init_tao_wave_struct(m, py_TaoWaveStruct);
  init_tao_wave_kick_pt_struct(m, py_TaoWaveKickPtStruct);
  init_tao_cmd_history_struct(m, py_TaoCmdHistoryStruct);
  init_tao_universe_struct(m, py_TaoUniverseStruct);
  init_mad_energy_struct(m, py_MadEnergyStruct);
  init_mad_map_struct(m, py_MadMapStruct);
  init_random_state_struct(m, py_RandomStateStruct);
  init_bbu_stage_struct(m, py_BbuStageStruct);
  init_bbu_beam_struct(m, py_BbuBeamStruct);
  init_bbu_param_struct(m, py_BbuParamStruct);
  init_fibre(m, py_Fibre);
  init_layout(m, py_Layout);
  init_all_encompassing_struct(m, py_AllEncompassingStruct);
  init_test_sub_struct(m, py_TestSubStruct);
  init_test_sub_sub_struct(m, py_TestSubSubStruct);

  // Hand-written bindings
  init_common_structs(m);

  // Routine submodules (one per C++ namespace)
  auto m_bmad = m.def_submodule("bmad", "Bmad routines");
  auto m_extra = m.def_submodule("extra", "CppBmadExtra routines");
  auto m_test = m.def_submodule("test", "CppBmadTest routines");
  auto m_simutils = m.def_submodule("simutils", "SimUtils routines");
  auto m_tao = m.def_submodule("tao", "Tao routines");
  auto m_bsim = m.def_submodule("bsim", "bsim routines");

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
