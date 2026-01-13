#include <pybind11/complex.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <string>

#include "bmad/json.hpp"
#include "pybmad/common_structs.hpp"
#include "pybmad/generated/init.hpp"
#include "pybmad/util.hpp"

using namespace Bmad;
using namespace Pybmad;

namespace py = pybind11;
using namespace pybind11::literals;

#include "pybmad/generated/routines.hpp"
#include "pybmad/generated/structs.hpp"

PYBIND11_MODULE(_pybmad, m) {
  // Generated definitions:
  m.doc() = "pybmad";

  bind_standard_arrays(m);

  // Structures
  auto py_SplineStruct = py::class_<SplineProxy>(
      m, "SplineStruct", "Fortran struct: spline_struct");
  auto py_SpinPolarStruct = py::class_<SpinPolarProxy>(
      m, "SpinPolarStruct", "Fortran struct: spin_polar_struct");
  auto py_AcKickerTimeStruct = py::class_<AcKickerTimeProxy>(
      m, "AcKickerTimeStruct", "Fortran struct: ac_kicker_time_struct");
  auto py_AcKickerFreqStruct = py::class_<AcKickerFreqProxy>(
      m, "AcKickerFreqStruct", "Fortran struct: ac_kicker_freq_struct");
  auto py_AcKickerStruct = py::class_<AcKickerProxy>(
      m, "AcKickerStruct", "Fortran struct: ac_kicker_struct");
  auto py_Interval1CoefStruct = py::class_<Interval1CoefProxy>(
      m, "Interval1CoefStruct", "Fortran struct: interval1_coef_struct");
  auto py_PhotonReflectTableStruct = py::class_<PhotonReflectTableProxy>(
      m,
      "PhotonReflectTableStruct",
      "Fortran struct: photon_reflect_table_struct");
  auto py_PhotonReflectSurfaceStruct = py::class_<PhotonReflectSurfaceProxy>(
      m,
      "PhotonReflectSurfaceStruct",
      "Fortran struct: photon_reflect_surface_struct");
  auto py_CoordStruct =
      py::class_<CoordProxy>(m, "CoordStruct", "Fortran struct: coord_struct");
  auto py_CoordArrayStruct = py::class_<CoordArrayProxy>(
      m, "CoordArrayStruct", "Fortran struct: coord_array_struct");
  auto py_BpmPhaseCouplingStruct = py::class_<BpmPhaseCouplingProxy>(
      m, "BpmPhaseCouplingStruct", "Fortran struct: bpm_phase_coupling_struct");
  auto py_ExpressionAtomStruct = py::class_<ExpressionAtomProxy>(
      m, "ExpressionAtomStruct", "Fortran struct: expression_atom_struct");
  auto py_WakeSrZLongStruct = py::class_<WakeSrZLongProxy>(
      m, "WakeSrZLongStruct", "Fortran struct: wake_sr_z_long_struct");
  auto py_WakeSrModeStruct = py::class_<WakeSrModeProxy>(
      m, "WakeSrModeStruct", "Fortran struct: wake_sr_mode_struct");
  auto py_WakeSrStruct = py::class_<WakeSrProxy>(
      m, "WakeSrStruct", "Fortran struct: wake_sr_struct");
  auto py_WakeLrModeStruct = py::class_<WakeLrModeProxy>(
      m, "WakeLrModeStruct", "Fortran struct: wake_lr_mode_struct");
  auto py_WakeLrStruct = py::class_<WakeLrProxy>(
      m, "WakeLrStruct", "Fortran struct: wake_lr_struct");
  auto py_LatEleLocStruct = py::class_<LatEleLocProxy>(
      m, "LatEleLocStruct", "Fortran struct: lat_ele_loc_struct");
  auto py_WakeStruct =
      py::class_<WakeProxy>(m, "WakeStruct", "Fortran struct: wake_struct");
  auto py_TaylorTermStruct = py::class_<TaylorTermProxy>(
      m, "TaylorTermStruct", "Fortran struct: taylor_term_struct");
  auto py_TaylorStruct = py::class_<TaylorProxy>(
      m, "TaylorStruct", "Fortran struct: taylor_struct");
  auto py_EmTaylorTermStruct = py::class_<EmTaylorTermProxy>(
      m, "EmTaylorTermStruct", "Fortran struct: em_taylor_term_struct");
  auto py_EmTaylorStruct = py::class_<EmTaylorProxy>(
      m, "EmTaylorStruct", "Fortran struct: em_taylor_struct");
  auto py_CartesianMapTerm1Struct = py::class_<CartesianMapTerm1Proxy>(
      m,
      "CartesianMapTerm1Struct",
      "Fortran struct: cartesian_map_term1_struct");
  auto py_CartesianMapTermStruct = py::class_<CartesianMapTermProxy>(
      m, "CartesianMapTermStruct", "Fortran struct: cartesian_map_term_struct");
  auto py_CartesianMapStruct = py::class_<CartesianMapProxy>(
      m, "CartesianMapStruct", "Fortran struct: cartesian_map_struct");
  auto py_CylindricalMapTerm1Struct = py::class_<CylindricalMapTerm1Proxy>(
      m,
      "CylindricalMapTerm1Struct",
      "Fortran struct: cylindrical_map_term1_struct");
  auto py_CylindricalMapTermStruct = py::class_<CylindricalMapTermProxy>(
      m,
      "CylindricalMapTermStruct",
      "Fortran struct: cylindrical_map_term_struct");
  auto py_CylindricalMapStruct = py::class_<CylindricalMapProxy>(
      m, "CylindricalMapStruct", "Fortran struct: cylindrical_map_struct");
  auto py_BicubicCmplxCoefStruct = py::class_<BicubicCmplxCoefProxy>(
      m, "BicubicCmplxCoefStruct", "Fortran struct: bicubic_cmplx_coef_struct");
  auto py_TricubicCmplxCoefStruct = py::class_<TricubicCmplxCoefProxy>(
      m,
      "TricubicCmplxCoefStruct",
      "Fortran struct: tricubic_cmplx_coef_struct");
  auto py_GridFieldPt1Struct = py::class_<GridFieldPt1Proxy>(
      m, "GridFieldPt1Struct", "Fortran struct: grid_field_pt1_struct");
  auto py_GridFieldPtStruct = py::class_<GridFieldPtProxy>(
      m, "GridFieldPtStruct", "Fortran struct: grid_field_pt_struct");
  auto py_GridFieldStruct = py::class_<GridFieldProxy>(
      m, "GridFieldStruct", "Fortran struct: grid_field_struct");
  auto py_FloorPositionStruct = py::class_<FloorPositionProxy>(
      m, "FloorPositionStruct", "Fortran struct: floor_position_struct");
  auto py_HighEnergySpaceChargeStruct = py::class_<HighEnergySpaceChargeProxy>(
      m,
      "HighEnergySpaceChargeStruct",
      "Fortran struct: high_energy_space_charge_struct");
  auto py_XyDispStruct = py::class_<XyDispProxy>(
      m, "XyDispStruct", "Fortran struct: xy_disp_struct");
  auto py_TwissStruct =
      py::class_<TwissProxy>(m, "TwissStruct", "Fortran struct: twiss_struct");
  auto py_Mode3Struct =
      py::class_<Mode3Proxy>(m, "Mode3Struct", "Fortran struct: mode3_struct");
  auto py_BookkeepingStateStruct = py::class_<BookkeepingStateProxy>(
      m, "BookkeepingStateStruct", "Fortran struct: bookkeeping_state_struct");
  auto py_RadMapStruct = py::class_<RadMapProxy>(
      m, "RadMapStruct", "Fortran struct: rad_map_struct");
  auto py_RadMapEleStruct = py::class_<RadMapEleProxy>(
      m, "RadMapEleStruct", "Fortran struct: rad_map_ele_struct");
  auto py_GenGrad1Struct = py::class_<GenGrad1Proxy>(
      m, "GenGrad1Struct", "Fortran struct: gen_grad1_struct");
  auto py_GenGradMapStruct = py::class_<GenGradMapProxy>(
      m, "GenGradMapStruct", "Fortran struct: gen_grad_map_struct");
  auto py_SurfaceSegmentedPtStruct = py::class_<SurfaceSegmentedPtProxy>(
      m,
      "SurfaceSegmentedPtStruct",
      "Fortran struct: surface_segmented_pt_struct");
  auto py_SurfaceSegmentedStruct = py::class_<SurfaceSegmentedProxy>(
      m, "SurfaceSegmentedStruct", "Fortran struct: surface_segmented_struct");
  auto py_SurfaceHMisalignPtStruct = py::class_<SurfaceHMisalignPtProxy>(
      m,
      "SurfaceHMisalignPtStruct",
      "Fortran struct: surface_h_misalign_pt_struct");
  auto py_SurfaceHMisalignStruct = py::class_<SurfaceHMisalignProxy>(
      m, "SurfaceHMisalignStruct", "Fortran struct: surface_h_misalign_struct");
  auto py_SurfaceDisplacementPtStruct = py::class_<SurfaceDisplacementPtProxy>(
      m,
      "SurfaceDisplacementPtStruct",
      "Fortran struct: surface_displacement_pt_struct");
  auto py_SurfaceDisplacementStruct = py::class_<SurfaceDisplacementProxy>(
      m,
      "SurfaceDisplacementStruct",
      "Fortran struct: surface_displacement_struct");
  auto py_TargetPointStruct = py::class_<TargetPointProxy>(
      m, "TargetPointStruct", "Fortran struct: target_point_struct");
  auto py_SurfaceCurvatureStruct = py::class_<SurfaceCurvatureProxy>(
      m, "SurfaceCurvatureStruct", "Fortran struct: surface_curvature_struct");
  auto py_PhotonTargetStruct = py::class_<PhotonTargetProxy>(
      m, "PhotonTargetStruct", "Fortran struct: photon_target_struct");
  auto py_PhotonMaterialStruct = py::class_<PhotonMaterialProxy>(
      m, "PhotonMaterialStruct", "Fortran struct: photon_material_struct");
  auto py_PixelPtStruct = py::class_<PixelPtProxy>(
      m, "PixelPtStruct", "Fortran struct: pixel_pt_struct");
  auto py_PixelDetecStruct = py::class_<PixelDetecProxy>(
      m, "PixelDetecStruct", "Fortran struct: pixel_detec_struct");
  auto py_PhotonElementStruct = py::class_<PhotonElementProxy>(
      m, "PhotonElementStruct", "Fortran struct: photon_element_struct");
  auto py_Wall3DVertexStruct = py::class_<Wall3dVertexProxy>(
      m, "Wall3DVertexStruct", "Fortran struct: wall3d_vertex_struct");
  auto py_Wall3DSectionStruct = py::class_<Wall3dSectionProxy>(
      m, "Wall3DSectionStruct", "Fortran struct: wall3d_section_struct");
  auto py_Wall3DStruct = py::class_<Wall3dProxy>(
      m, "Wall3DStruct", "Fortran struct: wall3d_struct");
  auto py_RamperLordStruct = py::class_<RamperLordProxy>(
      m, "RamperLordStruct", "Fortran struct: ramper_lord_struct");
  auto py_ControlStruct = py::class_<ControlProxy>(
      m, "ControlStruct", "Fortran struct: control_struct");
  auto py_ControlVar1Struct = py::class_<ControlVar1Proxy>(
      m, "ControlVar1Struct", "Fortran struct: control_var1_struct");
  auto py_ControlRamp1Struct = py::class_<ControlRamp1Proxy>(
      m, "ControlRamp1Struct", "Fortran struct: control_ramp1_struct");
  auto py_ControllerStruct = py::class_<ControllerProxy>(
      m, "ControllerStruct", "Fortran struct: controller_struct");
  auto py_EllipseBeamInitStruct = py::class_<EllipseBeamInitProxy>(
      m, "EllipseBeamInitStruct", "Fortran struct: ellipse_beam_init_struct");
  auto py_KvBeamInitStruct = py::class_<KvBeamInitProxy>(
      m, "KvBeamInitStruct", "Fortran struct: kv_beam_init_struct");
  auto py_GridBeamInitStruct = py::class_<GridBeamInitProxy>(
      m, "GridBeamInitStruct", "Fortran struct: grid_beam_init_struct");
  auto py_BeamInitStruct = py::class_<BeamInitProxy>(
      m, "BeamInitStruct", "Fortran struct: beam_init_struct");
  auto py_LatParamStruct = py::class_<LatParamProxy>(
      m, "LatParamStruct", "Fortran struct: lat_param_struct");
  auto py_ModeInfoStruct = py::class_<ModeInfoProxy>(
      m, "ModeInfoStruct", "Fortran struct: mode_info_struct");
  auto py_PreTrackerStruct = py::class_<PreTrackerProxy>(
      m, "PreTrackerStruct", "Fortran struct: pre_tracker_struct");
  auto py_AnormalModeStruct = py::class_<AnormalModeProxy>(
      m, "AnormalModeStruct", "Fortran struct: anormal_mode_struct");
  auto py_LinacNormalModeStruct = py::class_<LinacNormalModeProxy>(
      m, "LinacNormalModeStruct", "Fortran struct: linac_normal_mode_struct");
  auto py_NormalModesStruct = py::class_<NormalModesProxy>(
      m, "NormalModesStruct", "Fortran struct: normal_modes_struct");
  auto py_EmFieldStruct = py::class_<EmFieldProxy>(
      m, "EmFieldStruct", "Fortran struct: em_field_struct");
  auto py_StrongBeamStruct = py::class_<StrongBeamProxy>(
      m, "StrongBeamStruct", "Fortran struct: strong_beam_struct");
  auto py_TrackPointStruct = py::class_<TrackPointProxy>(
      m, "TrackPointStruct", "Fortran struct: track_point_struct");
  auto py_TrackStruct =
      py::class_<TrackProxy>(m, "TrackStruct", "Fortran struct: track_struct");
  auto py_SpaceChargeCommonStruct = py::class_<SpaceChargeCommonProxy>(
      m,
      "SpaceChargeCommonStruct",
      "Fortran struct: space_charge_common_struct");
  auto py_BmadCommonStruct = py::class_<BmadCommonProxy>(
      m, "BmadCommonStruct", "Fortran struct: bmad_common_struct");
  auto py_RadInt1Struct = py::class_<RadInt1Proxy>(
      m, "RadInt1Struct", "Fortran struct: rad_int1_struct");
  auto py_RadIntBranchStruct = py::class_<RadIntBranchProxy>(
      m, "RadIntBranchStruct", "Fortran struct: rad_int_branch_struct");
  auto py_RadIntAllEleStruct = py::class_<RadIntAllEleProxy>(
      m, "RadIntAllEleStruct", "Fortran struct: rad_int_all_ele_struct");
  auto py_RfStairStepStruct = py::class_<RfStairStepProxy>(
      m, "RfStairStepStruct", "Fortran struct: rf_stair_step_struct");
  auto py_RfEleStruct =
      py::class_<RfEleProxy>(m, "RfEleStruct", "Fortran struct: rf_ele_struct");
  auto py_EleStruct =
      py::class_<EleProxy>(m, "EleStruct", "Fortran struct: ele_struct");
  auto py_ComplexTaylorTermStruct = py::class_<ComplexTaylorTermProxy>(
      m,
      "ComplexTaylorTermStruct",
      "Fortran struct: complex_taylor_term_struct");
  auto py_ComplexTaylorStruct = py::class_<ComplexTaylorProxy>(
      m, "ComplexTaylorStruct", "Fortran struct: complex_taylor_struct");
  auto py_BranchStruct = py::class_<BranchProxy>(
      m, "BranchStruct", "Fortran struct: branch_struct");
  auto py_LatStruct =
      py::class_<LatProxy>(m, "LatStruct", "Fortran struct: lat_struct");
  auto py_BunchStruct =
      py::class_<BunchProxy>(m, "BunchStruct", "Fortran struct: bunch_struct");
  auto py_BunchParamsStruct = py::class_<BunchParamsProxy>(
      m, "BunchParamsStruct", "Fortran struct: bunch_params_struct");
  auto py_BeamStruct =
      py::class_<BeamProxy>(m, "BeamStruct", "Fortran struct: beam_struct");
  auto py_AperturePointStruct = py::class_<AperturePointProxy>(
      m, "AperturePointStruct", "Fortran struct: aperture_point_struct");
  auto py_ApertureParamStruct = py::class_<ApertureParamProxy>(
      m, "ApertureParamStruct", "Fortran struct: aperture_param_struct");
  auto py_ApertureScanStruct = py::class_<ApertureScanProxy>(
      m, "ApertureScanStruct", "Fortran struct: aperture_scan_struct");
  auto py_ElePointerStruct = py::class_<ElePointerProxy>(
      m, "ElePointerStruct", "Fortran struct: ele_pointer_struct");
  auto py_ExpressionTreeStruct = py::class_<ExpressionTreeProxy>(
      m, "ExpressionTreeStruct", "Fortran struct: expression_tree_struct");
  auto py_NametableStruct = py::class_<NametableProxy>(
      m, "NametableStruct", "Fortran struct: nametable_struct");
  auto py_TaoSpinDnDpzStruct = py::class_<TaoSpinDnDpzProxy>(
      m, "TaoSpinDnDpzStruct", "Fortran struct: tao_spin_dn_dpz_struct");
  auto py_ResonanceHStruct = py::class_<ResonanceHProxy>(
      m, "ResonanceHStruct", "Fortran struct: resonance_h_struct");
  auto py_SpinOrbitMap1Struct = py::class_<SpinOrbitMap1Proxy>(
      m, "SpinOrbitMap1Struct", "Fortran struct: spin_orbit_map1_struct");
  auto py_SpinAxisStruct = py::class_<SpinAxisProxy>(
      m, "SpinAxisStruct", "Fortran struct: spin_axis_struct");
  auto py_PtcNormalFormStruct = py::class_<PtcNormalFormProxy>(
      m, "PtcNormalFormStruct", "Fortran struct: ptc_normal_form_struct");
  auto py_BmadNormalFormStruct = py::class_<BmadNormalFormProxy>(
      m, "BmadNormalFormStruct", "Fortran struct: bmad_normal_form_struct");
  auto py_BunchTrackStruct = py::class_<BunchTrackProxy>(
      m, "BunchTrackStruct", "Fortran struct: bunch_track_struct");
  auto py_SummationRdtStruct = py::class_<SummationRdtProxy>(
      m, "SummationRdtStruct", "Fortran struct: summation_rdt_struct");
  auto py_TaoEleShapeStruct = py::class_<TaoEleShapeProxy>(
      m, "TaoEleShapeStruct", "Fortran struct: tao_ele_shape_struct");
  auto py_TaoElePointerStruct = py::class_<TaoElePointerProxy>(
      m, "TaoElePointerStruct", "Fortran struct: tao_ele_pointer_struct");
  auto py_TaoCurveStruct = py::class_<TaoCurveProxy>(
      m, "TaoCurveStruct", "Fortran struct: tao_curve_struct");
  auto py_TaoCurveColorStruct = py::class_<TaoCurveColorProxy>(
      m, "TaoCurveColorStruct", "Fortran struct: tao_curve_color_struct");
  auto py_TaoCurveOrbitStruct = py::class_<TaoCurveOrbitProxy>(
      m, "TaoCurveOrbitStruct", "Fortran struct: tao_curve_orbit_struct");
  auto py_TaoHistogramStruct = py::class_<TaoHistogramProxy>(
      m, "TaoHistogramStruct", "Fortran struct: tao_histogram_struct");
  auto py_LatEleOrder1Struct = py::class_<LatEleOrder1Proxy>(
      m, "LatEleOrder1Struct", "Fortran struct: lat_ele_order1_struct");
  auto py_LatEleOrderArrayStruct = py::class_<LatEleOrderArrayProxy>(
      m,
      "LatEleOrderArrayStruct",
      "Fortran struct: lat_ele_order_array_struct");
  auto py_TaoLatSigmaStruct = py::class_<TaoLatSigmaProxy>(
      m, "TaoLatSigmaStruct", "Fortran struct: tao_lat_sigma_struct");
  auto py_TaoSpinEleStruct = py::class_<TaoSpinEleProxy>(
      m, "TaoSpinEleStruct", "Fortran struct: tao_spin_ele_struct");
  auto py_TaoPlotCacheStruct = py::class_<TaoPlotCacheProxy>(
      m, "TaoPlotCacheStruct", "Fortran struct: tao_plot_cache_struct");
  auto py_TaoSpinPolarizationStruct = py::class_<TaoSpinPolarizationProxy>(
      m,
      "TaoSpinPolarizationStruct",
      "Fortran struct: tao_spin_polarization_struct");
  auto py_TaoLatticeBranchStruct = py::class_<TaoLatticeBranchProxy>(
      m, "TaoLatticeBranchStruct", "Fortran struct: tao_lattice_branch_struct");
  auto py_TaoModelElementStruct = py::class_<TaoModelElementProxy>(
      m, "TaoModelElementStruct", "Fortran struct: tao_model_element_struct");
  auto py_TaoBeamBranchStruct = py::class_<TaoBeamBranchProxy>(
      m, "TaoBeamBranchStruct", "Fortran struct: tao_beam_branch_struct");
  auto py_TaoD1DataStruct = py::class_<TaoD1DataProxy>(
      m, "TaoD1DataStruct", "Fortran struct: tao_d1_data_struct");
  auto py_TaoD2DataStruct = py::class_<TaoD2DataProxy>(
      m, "TaoD2DataStruct", "Fortran struct: tao_d2_data_struct");
  auto py_TaoDataVarComponentStruct = py::class_<TaoDataVarComponentProxy>(
      m,
      "TaoDataVarComponentStruct",
      "Fortran struct: tao_data_var_component_struct");
  auto py_TaoGraphStruct = py::class_<TaoGraphProxy>(
      m, "TaoGraphStruct", "Fortran struct: tao_graph_struct");
  auto py_TaoPlotStruct = py::class_<TaoPlotProxy>(
      m, "TaoPlotStruct", "Fortran struct: tao_plot_struct");
  auto py_TaoPlotRegionStruct = py::class_<TaoPlotRegionProxy>(
      m, "TaoPlotRegionStruct", "Fortran struct: tao_plot_region_struct");
  auto py_TaoUniversePointerStruct = py::class_<TaoUniversePointerProxy>(
      m,
      "TaoUniversePointerStruct",
      "Fortran struct: tao_universe_pointer_struct");
  auto py_TaoSuperUniverseStruct = py::class_<TaoSuperUniverseProxy>(
      m, "TaoSuperUniverseStruct", "Fortran struct: tao_super_universe_struct");
  auto py_TaoVarStruct = py::class_<TaoVarProxy>(
      m, "TaoVarStruct", "Fortran struct: tao_var_struct");
  auto py_TaoVarSlaveStruct = py::class_<TaoVarSlaveProxy>(
      m, "TaoVarSlaveStruct", "Fortran struct: tao_var_slave_struct");
  auto py_TaoLatticeStruct = py::class_<TaoLatticeProxy>(
      m, "TaoLatticeStruct", "Fortran struct: tao_lattice_struct");
  auto py_TaoBeamUniStruct = py::class_<TaoBeamUniProxy>(
      m, "TaoBeamUniStruct", "Fortran struct: tao_beam_uni_struct");
  auto py_TaoDynamicApertureStruct = py::class_<TaoDynamicApertureProxy>(
      m,
      "TaoDynamicApertureStruct",
      "Fortran struct: tao_dynamic_aperture_struct");
  auto py_TaoModelBranchStruct = py::class_<TaoModelBranchProxy>(
      m, "TaoModelBranchStruct", "Fortran struct: tao_model_branch_struct");
  auto py_TaoSpinMapStruct = py::class_<TaoSpinMapProxy>(
      m, "TaoSpinMapStruct", "Fortran struct: tao_spin_map_struct");
  auto py_TaoDataStruct = py::class_<TaoDataProxy>(
      m, "TaoDataStruct", "Fortran struct: tao_data_struct");
  auto py_TaoPingScaleStruct = py::class_<TaoPingScaleProxy>(
      m, "TaoPingScaleStruct", "Fortran struct: tao_ping_scale_struct");
  auto py_TaoUniverseCalcStruct = py::class_<TaoUniverseCalcProxy>(
      m, "TaoUniverseCalcStruct", "Fortran struct: tao_universe_calc_struct");
  auto py_LatEleOrderStruct = py::class_<LatEleOrderProxy>(
      m, "LatEleOrderStruct", "Fortran struct: lat_ele_order_struct");
  auto py_TaoTitleStruct = py::class_<TaoTitleProxy>(
      m, "TaoTitleStruct", "Fortran struct: tao_title_struct");
  auto py_QpRectStruct = py::class_<QpRectProxy>(
      m, "QpRectStruct", "Fortran struct: qp_rect_struct");
  auto py_TaoDrawingStruct = py::class_<TaoDrawingProxy>(
      m, "TaoDrawingStruct", "Fortran struct: tao_drawing_struct");
  auto py_TaoShapePatternStruct = py::class_<TaoShapePatternProxy>(
      m, "TaoShapePatternStruct", "Fortran struct: tao_shape_pattern_struct");
  auto py_TaoShapePatternPointStruct = py::class_<TaoShapePatternPointProxy>(
      m,
      "TaoShapePatternPointStruct",
      "Fortran struct: tao_shape_pattern_point_struct");
  auto py_QpAxisStruct = py::class_<QpAxisProxy>(
      m, "QpAxisStruct", "Fortran struct: qp_axis_struct");
  auto py_QpLegendStruct = py::class_<QpLegendProxy>(
      m, "QpLegendStruct", "Fortran struct: qp_legend_struct");
  auto py_QpPointStruct = py::class_<QpPointProxy>(
      m, "QpPointStruct", "Fortran struct: qp_point_struct");
  auto py_QpLineStruct = py::class_<QpLineProxy>(
      m, "QpLineStruct", "Fortran struct: qp_line_struct");
  auto py_QpSymbolStruct = py::class_<QpSymbolProxy>(
      m, "QpSymbolStruct", "Fortran struct: qp_symbol_struct");
  auto py_TaoFloorPlanStruct = py::class_<TaoFloorPlanProxy>(
      m, "TaoFloorPlanStruct", "Fortran struct: tao_floor_plan_struct");
  auto py_TaoV1VarStruct = py::class_<TaoV1VarProxy>(
      m, "TaoV1VarStruct", "Fortran struct: tao_v1_var_struct");
  auto py_TaoGlobalStruct = py::class_<TaoGlobalProxy>(
      m, "TaoGlobalStruct", "Fortran struct: tao_global_struct");
  auto py_TaoInitStruct = py::class_<TaoInitProxy>(
      m, "TaoInitStruct", "Fortran struct: tao_init_struct");
  auto py_TaoCommonStruct = py::class_<TaoCommonProxy>(
      m, "TaoCommonStruct", "Fortran struct: tao_common_struct");
  auto py_TaoPlotPageStruct = py::class_<TaoPlotPageProxy>(
      m, "TaoPlotPageStruct", "Fortran struct: tao_plot_page_struct");
  auto py_TaoBuildingWallStruct = py::class_<TaoBuildingWallProxy>(
      m, "TaoBuildingWallStruct", "Fortran struct: tao_building_wall_struct");
  auto py_TaoBuildingWallOrientationStruct =
      py::class_<TaoBuildingWallOrientationProxy>(
          m,
          "TaoBuildingWallOrientationStruct",
          "Fortran struct: tao_building_wall_orientation_struct");
  auto py_TaoBuildingWallSectionStruct =
      py::class_<TaoBuildingWallSectionProxy>(
          m,
          "TaoBuildingWallSectionStruct",
          "Fortran struct: tao_building_wall_section_struct");
  auto py_TaoBuildingWallPointStruct = py::class_<TaoBuildingWallPointProxy>(
      m,
      "TaoBuildingWallPointStruct",
      "Fortran struct: tao_building_wall_point_struct");
  auto py_TaoWaveStruct = py::class_<TaoWaveProxy>(
      m, "TaoWaveStruct", "Fortran struct: tao_wave_struct");
  auto py_TaoWaveKickPtStruct = py::class_<TaoWaveKickPtProxy>(
      m, "TaoWaveKickPtStruct", "Fortran struct: tao_wave_kick_pt_struct");
  auto py_TaoCmdHistoryStruct = py::class_<TaoCmdHistoryProxy>(
      m, "TaoCmdHistoryStruct", "Fortran struct: tao_cmd_history_struct");
  auto py_TaoUniverseStruct = py::class_<TaoUniverseProxy>(
      m, "TaoUniverseStruct", "Fortran struct: tao_universe_struct");
  auto py_MadEnergyStruct = py::class_<MadEnergyProxy>(
      m, "MadEnergyStruct", "Fortran struct: mad_energy_struct");
  auto py_MadMapStruct = py::class_<MadMapProxy>(
      m, "MadMapStruct", "Fortran struct: mad_map_struct");
  auto py_RandomStateStruct = py::class_<RandomStateProxy>(
      m, "RandomStateStruct", "Fortran struct: random_state_struct");
  auto py_BbuStageStruct = py::class_<BbuStageProxy>(
      m, "BbuStageStruct", "Fortran struct: bbu_stage_struct");
  auto py_BbuBeamStruct = py::class_<BbuBeamProxy>(
      m, "BbuBeamStruct", "Fortran struct: bbu_beam_struct");
  auto py_BbuParamStruct = py::class_<BbuParamProxy>(
      m, "BbuParamStruct", "Fortran struct: bbu_param_struct");
  auto py_AllEncompassingStruct = py::class_<AllEncompassingProxy>(
      m, "AllEncompassingStruct", "Fortran struct: all_encompassing_struct");
  auto py_TestSubStruct = py::class_<TestSubProxy>(
      m, "TestSubStruct", "Fortran struct: test_sub_struct");
  auto py_TestSubSubStruct = py::class_<TestSubSubProxy>(
      m, "TestSubSubStruct", "Fortran struct: test_sub_sub_struct");
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
  init_em_taylor_term_struct(m, py_EmTaylorTermStruct);
  init_em_taylor_struct(m, py_EmTaylorStruct);
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
  init_tao_building_wall_orientation_struct(
      m, py_TaoBuildingWallOrientationStruct);
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
  init_all_encompassing_struct(m, py_AllEncompassingStruct);
  init_test_sub_struct(m, py_TestSubStruct);
  init_test_sub_sub_struct(m, py_TestSubSubStruct);

  // Hand-written bindings
  init_common_structs(m);

  // Routine initializers
  init_Bmad_routines_a(m);
  init_Bmad_routines_b(m);
  init_Bmad_routines_c(m);
  init_Bmad_routines_d(m);
  init_Bmad_routines_e(m);
  init_Bmad_routines_f(m);
  init_Bmad_routines_g(m);
  init_Bmad_routines_h(m);
  init_Bmad_routines_i(m);
  init_Bmad_routines_j(m);
  init_Bmad_routines_k(m);
  init_Bmad_routines_l(m);
  init_Bmad_routines_m(m);
  init_Bmad_routines_n(m);
  init_Bmad_routines_o(m);
  init_Bmad_routines_p(m);
  init_Bmad_routines_q(m);
  init_Bmad_routines_r(m);
  init_Bmad_routines_s(m);
  init_Bmad_routines_t(m);
  init_Bmad_routines_u(m);
  init_Bmad_routines_v(m);
  init_Bmad_routines_w(m);
  init_Bmad_routines_x(m);
  init_Bmad_routines_y(m);
  init_Bmad_routines_z(m);
  init_CppBmadTest_routines_t(m);
  init_SimUtils_routines_a(m);
  init_SimUtils_routines_b(m);
  init_SimUtils_routines_c(m);
  init_SimUtils_routines_d(m);
  init_SimUtils_routines_e(m);
  init_SimUtils_routines_f(m);
  init_SimUtils_routines_g(m);
  init_SimUtils_routines_h(m);
  init_SimUtils_routines_i(m);
  init_SimUtils_routines_j(m);
  init_SimUtils_routines_l(m);
  init_SimUtils_routines_m(m);
  init_SimUtils_routines_n(m);
  init_SimUtils_routines_o(m);
  init_SimUtils_routines_p(m);
  init_SimUtils_routines_q(m);
  init_SimUtils_routines_r(m);
  init_SimUtils_routines_s(m);
  init_SimUtils_routines_t(m);
  init_SimUtils_routines_u(m);
  init_SimUtils_routines_v(m);
  init_SimUtils_routines_w(m);
  init_SimUtils_routines_x(m);
  init_Tao_routines_a(m);
  init_Tao_routines_c(m);
  init_Tao_routines_i(m);
  init_Tao_routines_j(m);
  init_Tao_routines_r(m);
  init_Tao_routines_t(m);
  init_Tao_routines_u(m);
  init_bsim_routines_b(m);
  init_bsim_routines_c(m);
  init_bsim_routines_h(m);
  init_bsim_routines_i(m);
  init_bsim_routines_l(m);
  init_bsim_routines_r(m);
  init_bsim_routines_s(m);
  init_bsim_routines_t(m);
  init_bsim_routines_w(m);

  // Enums
  py::native_enum<EleAttribute>(m, "EleAttribute", "enum.IntEnum")
      .value(
          "L",
          EleAttribute::L,
          "Assumed unique. Do not assign 1 to another attribute.")
      .value("TILT", EleAttribute::TILT, "Important: tilt$ = roll$")
      .value("ROLL", EleAttribute::ROLL)
      .value("N_PART", EleAttribute::N_PART)
      .value("INHERIT_FROM_FORK", EleAttribute::INHERIT_FROM_FORK)
      .value("REF_TILT", EleAttribute::REF_TILT)
      .value("DIRECTION", EleAttribute::DIRECTION)
      .value("REPETITION_FREQUENCY", EleAttribute::REPETITION_FREQUENCY)
      .value("DETA_DS_MASTER", EleAttribute::DETA_DS_MASTER)
      .value("KICK", EleAttribute::KICK)
      .value("X_GAIN_ERR", EleAttribute::X_GAIN_ERR)
      .value("TAYLOR_ORDER", EleAttribute::TAYLOR_ORDER)
      .value("R_SOLENOID", EleAttribute::R_SOLENOID)
      .value("FINAL_CHARGE", EleAttribute::FINAL_CHARGE)
      .value("K1", EleAttribute::K1)
      .value("KX", EleAttribute::KX)
      .value("HARMON", EleAttribute::HARMON)
      .value("H_DISPLACE", EleAttribute::H_DISPLACE)
      .value("Y_GAIN_ERR", EleAttribute::Y_GAIN_ERR)
      .value("S_TWISS_REF", EleAttribute::S_TWISS_REF)
      .value("CRITICAL_ANGLE_FACTOR", EleAttribute::CRITICAL_ANGLE_FACTOR)
      .value("TILT_CORR", EleAttribute::TILT_CORR)
      .value("REF_COORDS", EleAttribute::REF_COORDS)
      .value("DT_MAX", EleAttribute::DT_MAX)
      .value("IX_FIXER", EleAttribute::IX_FIXER)
      .value("GRAZE_ANGLE", EleAttribute::GRAZE_ANGLE)
      .value("K2", EleAttribute::K2)
      .value("B_MAX", EleAttribute::B_MAX)
      .value("V_DISPLACE", EleAttribute::V_DISPLACE)
      .value("GRADIENT_TOT", EleAttribute::GRADIENT_TOT)
      .value("HARMON_MASTER", EleAttribute::HARMON_MASTER)
      .value("FLEXIBLE", EleAttribute::FLEXIBLE)
      .value("CRUNCH", EleAttribute::CRUNCH)
      .value("REF_ORBIT_FOLLOWS", EleAttribute::REF_ORBIT_FOLLOWS)
      .value("PC_OUT_MIN", EleAttribute::PC_OUT_MIN)
      .value("GRADIENT", EleAttribute::GRADIENT)
      .value("K3", EleAttribute::K3)
      .value("NOISE", EleAttribute::NOISE)
      .value("NEW_BRANCH", EleAttribute::NEW_BRANCH)
      .value("IX_BRANCH", EleAttribute::IX_BRANCH)
      .value("G_MAX", EleAttribute::G_MAX)
      .value("G", EleAttribute::G)
      .value("SYMMETRY", EleAttribute::SYMMETRY)
      .value("FIELD_SCALE_FACTOR", EleAttribute::FIELD_SCALE_FACTOR)
      .value("PC_OUT_MAX", EleAttribute::PC_OUT_MAX)
      .value("DG", EleAttribute::DG)
      .value("BBI_CONST", EleAttribute::BBI_CONST)
      .value("OSC_AMPLITUDE", EleAttribute::OSC_AMPLITUDE)
      .value("IX_TO_BRANCH", EleAttribute::IX_TO_BRANCH)
      .value("ANGLE_OUT_MAX", EleAttribute::ANGLE_OUT_MAX)
      .value("GRADIENT_ERR", EleAttribute::GRADIENT_ERR)
      .value("CRITICAL_ANGLE", EleAttribute::CRITICAL_ANGLE)
      .value("BRAGG_ANGLE_IN", EleAttribute::BRAGG_ANGLE_IN)
      .value("SPIN_DN_DPZ_X", EleAttribute::SPIN_DN_DPZ_X)
      .value("DELTA_E_REF", EleAttribute::DELTA_E_REF)
      .value("INTERPOLATION", EleAttribute::INTERPOLATION)
      .value("BRAGG_ANGLE_OUT", EleAttribute::BRAGG_ANGLE_OUT)
      .value("K1X", EleAttribute::K1X)
      .value("SPIN_DN_DPZ_Y", EleAttribute::SPIN_DN_DPZ_Y)
      .value("CHARGE", EleAttribute::CHARGE)
      .value("X_GAIN_CALIB", EleAttribute::X_GAIN_CALIB)
      .value("IX_TO_ELEMENT", EleAttribute::IX_TO_ELEMENT)
      .value("VOLTAGE", EleAttribute::VOLTAGE)
      .value("G_TOT", EleAttribute::G_TOT)
      .value("RHO", EleAttribute::RHO)
      .value("VOLTAGE_ERR", EleAttribute::VOLTAGE_ERR)
      .value("BRAGG_ANGLE", EleAttribute::BRAGG_ANGLE)
      .value("K1Y", EleAttribute::K1Y)
      .value("N_PARTICLE", EleAttribute::N_PARTICLE)
      .value("SPIN_DN_DPZ_Z", EleAttribute::SPIN_DN_DPZ_Z)
      .value("FRINGE_TYPE", EleAttribute::FRINGE_TYPE)
      .value("DBRAGG_ANGLE_DE", EleAttribute::DBRAGG_ANGLE_DE)
      .value("FRINGE_AT", EleAttribute::FRINGE_AT)
      .value("GANG", EleAttribute::GANG)
      .value("DARWIN_WIDTH_SIGMA", EleAttribute::DARWIN_WIDTH_SIGMA)
      .value("DARWIN_WIDTH_PI", EleAttribute::DARWIN_WIDTH_PI)
      .value("SPIN_FRINGE_ON", EleAttribute::SPIN_FRINGE_ON)
      .value(
          "PENDELLOSUNG_PERIOD_SIGMA", EleAttribute::PENDELLOSUNG_PERIOD_SIGMA)
      .value("SIG_X", EleAttribute::SIG_X)
      .value("EXACT_MULTIPOLES", EleAttribute::EXACT_MULTIPOLES)
      .value("PENDELLOSUNG_PERIOD_PI", EleAttribute::PENDELLOSUNG_PERIOD_PI)
      .value("SIG_Y", EleAttribute::SIG_Y)
      .value("GRAZE_ANGLE_IN", EleAttribute::GRAZE_ANGLE_IN)
      .value("R0_ELEC", EleAttribute::R0_ELEC)
      .value("RF_FREQUENCY", EleAttribute::RF_FREQUENCY)
      .value("SIG_Z", EleAttribute::SIG_Z)
      .value("GRAZE_ANGLE_OUT", EleAttribute::GRAZE_ANGLE_OUT)
      .value("R0_MAG", EleAttribute::R0_MAG)
      .value("RF_WAVELENGTH", EleAttribute::RF_WAVELENGTH)
      .value("SIG_VX", EleAttribute::SIG_VX)
      .value("SIG_VY", EleAttribute::SIG_VY)
      .value("CONSTANT_REF_ENERGY", EleAttribute::CONSTANT_REF_ENERGY)
      .value("KS", EleAttribute::KS)
      .value("SIG_E", EleAttribute::SIG_E)
      .value("SIG_PZ", EleAttribute::SIG_PZ)
      .value("AUTOSCALE_AMPLITUDE", EleAttribute::AUTOSCALE_AMPLITUDE)
      .value("D1_THICKNESS", EleAttribute::D1_THICKNESS)
      .value("DEFAULT_TRACKING_SPECIES", EleAttribute::DEFAULT_TRACKING_SPECIES)
      .value("AUTOSCALE_PHASE", EleAttribute::AUTOSCALE_PHASE)
      .value("N_SLICE", EleAttribute::N_SLICE)
      .value("Y_GAIN_CALIB", EleAttribute::Y_GAIN_CALIB)
      .value("SIG_E2", EleAttribute::SIG_E2)
      .value("FB1", EleAttribute::FB1)
      .value("POLARITY", EleAttribute::POLARITY)
      .value("CRUNCH_CALIB", EleAttribute::CRUNCH_CALIB)
      .value("ALPHA_ANGLE", EleAttribute::ALPHA_ANGLE)
      .value("D2_THICKNESS", EleAttribute::D2_THICKNESS)
      .value("BETA_A_STRONG", EleAttribute::BETA_A_STRONG)
      .value("BETA_A_OUT", EleAttribute::BETA_A_OUT)
      .value("E_LOSS", EleAttribute::E_LOSS)
      .value("GAP", EleAttribute::GAP)
      .value("SPIN_X", EleAttribute::SPIN_X)
      .value("E_CENTER", EleAttribute::E_CENTER)
      .value("SCATTER_TEST", EleAttribute::SCATTER_TEST)
      .value("FB2", EleAttribute::FB2)
      .value("X_OFFSET_CALIB", EleAttribute::X_OFFSET_CALIB)
      .value("V1_UNITCELL", EleAttribute::V1_UNITCELL)
      .value("PSI_ANGLE", EleAttribute::PSI_ANGLE)
      .value("CAVITY_TYPE", EleAttribute::CAVITY_TYPE)
      .value("BETA_B_STRONG", EleAttribute::BETA_B_STRONG)
      .value("BETA_B_OUT", EleAttribute::BETA_B_OUT)
      .value("SPIN_Y", EleAttribute::SPIN_Y)
      .value("E2_CENTER", EleAttribute::E2_CENTER)
      .value("N_PERIOD", EleAttribute::N_PERIOD)
      .value("EMIT_FRACTION", EleAttribute::EMIT_FRACTION)
      .value("X1_EDGE", EleAttribute::X1_EDGE)
      .value("Y_OFFSET_CALIB", EleAttribute::Y_OFFSET_CALIB)
      .value("V_UNITCELL", EleAttribute::V_UNITCELL)
      .value("V2_UNITCELL", EleAttribute::V2_UNITCELL)
      .value("SPIN_Z", EleAttribute::SPIN_Z)
      .value("L_PERIOD", EleAttribute::L_PERIOD)
      .value("FQ1", EleAttribute::FQ1)
      .value("ALPHA_A_STRONG", EleAttribute::ALPHA_A_STRONG)
      .value("ALPHA_A_OUT", EleAttribute::ALPHA_A_OUT)
      .value("E2_PROBABILITY", EleAttribute::E2_PROBABILITY)
      .value("PHI0_MAX", EleAttribute::PHI0_MAX)
      .value("X2_EDGE", EleAttribute::X2_EDGE)
      .value("FQ2", EleAttribute::FQ2)
      .value("PHI0", EleAttribute::PHI0)
      .value("TILT_CALIB", EleAttribute::TILT_CALIB)
      .value("E_CENTER_RELATIVE_TO_REF", EleAttribute::E_CENTER_RELATIVE_TO_REF)
      .value("Y1_EDGE", EleAttribute::Y1_EDGE)
      .value("ALPHA_B_STRONG", EleAttribute::ALPHA_B_STRONG)
      .value("ALPHA_B_OUT", EleAttribute::ALPHA_B_OUT)
      .value("IS_MOSAIC", EleAttribute::IS_MOSAIC)
      .value("PX_APERTURE_WIDTH2", EleAttribute::PX_APERTURE_WIDTH2)
      .value("PHI0_ERR", EleAttribute::PHI0_ERR)
      .value("CURRENT", EleAttribute::CURRENT)
      .value("MOSAIC_THICKNESS", EleAttribute::MOSAIC_THICKNESS)
      .value("PX_APERTURE_CENTER", EleAttribute::PX_APERTURE_CENTER)
      .value("ETA_X_OUT", EleAttribute::ETA_X_OUT)
      .value("QUAD_TILT", EleAttribute::QUAD_TILT)
      .value("DE_ETA_MEAS", EleAttribute::DE_ETA_MEAS)
      .value("SPATIAL_DISTRIBUTION", EleAttribute::SPATIAL_DISTRIBUTION)
      .value("Y2_EDGE", EleAttribute::Y2_EDGE)
      .value("SPECIES_STRONG", EleAttribute::SPECIES_STRONG)
      .value("ETA_Y_OUT", EleAttribute::ETA_Y_OUT)
      .value("MODE", EleAttribute::MODE)
      .value("VELOCITY_DISTRIBUTION", EleAttribute::VELOCITY_DISTRIBUTION)
      .value("PY_APERTURE_WIDTH2", EleAttribute::PY_APERTURE_WIDTH2)
      .value("PHI0_MULTIPASS", EleAttribute::PHI0_MULTIPASS)
      .value("N_SAMPLE", EleAttribute::N_SAMPLE)
      .value("ORIGIN_ELE_REF_PT", EleAttribute::ORIGIN_ELE_REF_PT)
      .value(
          "MOSAIC_ANGLE_RMS_IN_PLANE", EleAttribute::MOSAIC_ANGLE_RMS_IN_PLANE)
      .value("EPS_STEP_SCALE", EleAttribute::EPS_STEP_SCALE)
      .value("E_TOT_STRONG", EleAttribute::E_TOT_STRONG)
      .value("DTHICKNESS_DX", EleAttribute::DTHICKNESS_DX)
      .value("BEND_TILT", EleAttribute::BEND_TILT)
      .value("ETAP_X_OUT", EleAttribute::ETAP_X_OUT)
      .value("PHI0_AUTOSCALE", EleAttribute::PHI0_AUTOSCALE)
      .value("DX_ORIGIN", EleAttribute::DX_ORIGIN)
      .value("ENERGY_DISTRIBUTION", EleAttribute::ENERGY_DISTRIBUTION)
      .value("X_QUAD", EleAttribute::X_QUAD)
      .value("DS_PHOTON_SLICE", EleAttribute::DS_PHOTON_SLICE)
      .value(
          "MOSAIC_ANGLE_RMS_OUT_PLANE",
          EleAttribute::MOSAIC_ANGLE_RMS_OUT_PLANE)
      .value("PY_APERTURE_CENTER", EleAttribute::PY_APERTURE_CENTER)
      .value("X_DISPERSION_ERR", EleAttribute::X_DISPERSION_ERR)
      .value("L_RECTANGLE", EleAttribute::L_RECTANGLE)
      .value("PC_STRONG", EleAttribute::PC_STRONG)
      .value("ETAP_Y_OUT", EleAttribute::ETAP_Y_OUT)
      .value("DY_ORIGIN", EleAttribute::DY_ORIGIN)
      .value("Y_QUAD", EleAttribute::Y_QUAD)
      .value("E_FIELD_X", EleAttribute::E_FIELD_X)
      .value("Y_DISPERSION_ERR", EleAttribute::Y_DISPERSION_ERR)
      .value("Z_APERTURE_WIDTH2", EleAttribute::Z_APERTURE_WIDTH2)
      .value("USER_SETS_LENGTH", EleAttribute::USER_SETS_LENGTH)
      .value("B_FIELD_TOT", EleAttribute::B_FIELD_TOT)
      .value("UPSTREAM_COORD_DIR", EleAttribute::UPSTREAM_COORD_DIR)
      .value("DZ_ORIGIN", EleAttribute::DZ_ORIGIN)
      .value("MOSAIC_DIFFRACTION_NUM", EleAttribute::MOSAIC_DIFFRACTION_NUM)
      .value("CMAT_11", EleAttribute::CMAT_11)
      .value("FIELD_AUTOSCALE", EleAttribute::FIELD_AUTOSCALE)
      .value("L_SAGITTA", EleAttribute::L_SAGITTA)
      .value("E_FIELD_Y", EleAttribute::E_FIELD_Y)
      .value("X_DISPERSION_CALIB", EleAttribute::X_DISPERSION_CALIB)
      .value("Z_APERTURE_CENTER", EleAttribute::Z_APERTURE_CENTER)
      .value("F_FACTOR", EleAttribute::F_FACTOR)
      .value("CMAT_12", EleAttribute::CMAT_12)
      .value("DTHETA_ORIGIN", EleAttribute::DTHETA_ORIGIN)
      .value("B_PARAM", EleAttribute::B_PARAM)
      .value("L_CHORD", EleAttribute::L_CHORD)
      .value("DOWNSTREAM_COORD_DIR", EleAttribute::DOWNSTREAM_COORD_DIR)
      .value("PZ_APERTURE_WIDTH2", EleAttribute::PZ_APERTURE_WIDTH2)
      .value("Y_DISPERSION_CALIB", EleAttribute::Y_DISPERSION_CALIB)
      .value("SCALE_FIELD_TO_ONE", EleAttribute::SCALE_FIELD_TO_ONE)
      .value("VOLTAGE_TOT", EleAttribute::VOLTAGE_TOT)
      .value("SCATTER_METHOD", EleAttribute::SCATTER_METHOD)
      .value("CMAT_21", EleAttribute::CMAT_21)
      .value("L_ACTIVE", EleAttribute::L_ACTIVE)
      .value("DPHI_ORIGIN", EleAttribute::DPHI_ORIGIN)
      .value("SPLIT_ID", EleAttribute::SPLIT_ID)
      .value("REF_CAP_GAMMA", EleAttribute::REF_CAP_GAMMA)
      .value("L_SOFT_EDGE", EleAttribute::L_SOFT_EDGE)
      .value("TRANSVERSE_SIGMA_CUT", EleAttribute::TRANSVERSE_SIGMA_CUT)
      .value("PZ_APERTURE_CENTER", EleAttribute::PZ_APERTURE_CENTER)
      .value("MEAN_EXCITATION_ENERGY", EleAttribute::MEAN_EXCITATION_ENERGY)
      .value("FIDUCIAL_PT", EleAttribute::FIDUCIAL_PT)
      .value("CMAT_22", EleAttribute::CMAT_22)
      .value("DPSI_ORIGIN", EleAttribute::DPSI_ORIGIN)
      .value("T_OFFSET", EleAttribute::T_OFFSET)
      .value("DS_SLICE", EleAttribute::DS_SLICE)
      .value("USE_REFLECTIVITY_TABLE", EleAttribute::USE_REFLECTIVITY_TABLE)
      .value("INIT_NEEDED", EleAttribute::INIT_NEEDED)
      .value("LONGITUDINAL_MODE", EleAttribute::LONGITUDINAL_MODE)
      .value("ANGLE", EleAttribute::ANGLE)
      .value("N_CELL", EleAttribute::N_CELL)
      .value("MODE_FLIP", EleAttribute::MODE_FLIP)
      .value("CROSSING_TIME", EleAttribute::CROSSING_TIME)
      .value("X_KICK", EleAttribute::X_KICK)
      .value(
          "X_PITCH",
          EleAttribute::X_PITCH,
          "Note: [x_kick$, px_kick$, ..., pz_kick$] must be in order.")
      .value("PX_KICK", EleAttribute::PX_KICK)
      .value("Y_PITCH", EleAttribute::Y_PITCH)
      .value("Y_KICK", EleAttribute::Y_KICK)
      .value("X_OFFSET", EleAttribute::X_OFFSET)
      .value("PY_KICK", EleAttribute::PY_KICK)
      .value("Y_OFFSET", EleAttribute::Y_OFFSET)
      .value("Z_KICK", EleAttribute::Z_KICK)
      .value("Z_OFFSET", EleAttribute::Z_OFFSET)
      .value("PZ_KICK", EleAttribute::PZ_KICK)
      .value("HKICK", EleAttribute::HKICK)
      .value("D_SPACING", EleAttribute::D_SPACING)
      .value("X_OFFSET_MULT", EleAttribute::X_OFFSET_MULT)
      .value("EMITTANCE_A", EleAttribute::EMITTANCE_A)
      .value("CRAB_X1", EleAttribute::CRAB_X1)
      .value("VKICK", EleAttribute::VKICK)
      .value("Y_OFFSET_MULT", EleAttribute::Y_OFFSET_MULT)
      .value("P0C_REF_INIT", EleAttribute::P0C_REF_INIT)
      .value("EMITTANCE_B", EleAttribute::EMITTANCE_B)
      .value("CRAB_X2", EleAttribute::CRAB_X2)
      .value("BL_HKICK", EleAttribute::BL_HKICK)
      .value("E_TOT_REF_INIT", EleAttribute::E_TOT_REF_INIT)
      .value("EMITTANCE_Z", EleAttribute::EMITTANCE_Z)
      .value("CRAB_X3", EleAttribute::CRAB_X3)
      .value("BL_VKICK", EleAttribute::BL_VKICK)
      .value("CRAB_TILT", EleAttribute::CRAB_TILT)
      .value("BL_KICK", EleAttribute::BL_KICK)
      .value("B_FIELD", EleAttribute::B_FIELD)
      .value("E_FIELD", EleAttribute::E_FIELD)
      .value(
          "HIGH_ENERGY_SPACE_CHARGE_ON",
          EleAttribute::HIGH_ENERGY_SPACE_CHARGE_ON)
      .value("CRAB_X4", EleAttribute::CRAB_X4)
      .value("N_RF_STEPS", EleAttribute::N_RF_STEPS)
      .value("PHOTON_TYPE", EleAttribute::PHOTON_TYPE)
      .value("COUPLER_PHASE", EleAttribute::COUPLER_PHASE)
      .value("DB_FIELD", EleAttribute::DB_FIELD)
      .value("CRAB_X5", EleAttribute::CRAB_X5)
      .value("LATTICE_TYPE", EleAttribute::LATTICE_TYPE)
      .value("B1_GRADIENT", EleAttribute::B1_GRADIENT)
      .value("E1_GRADIENT", EleAttribute::E1_GRADIENT)
      .value("COUPLER_ANGLE", EleAttribute::COUPLER_ANGLE)
      .value("LIVE_BRANCH", EleAttribute::LIVE_BRANCH)
      .value("B2_GRADIENT", EleAttribute::B2_GRADIENT)
      .value("E2_GRADIENT", EleAttribute::E2_GRADIENT)
      .value("COUPLER_STRENGTH", EleAttribute::COUPLER_STRENGTH)
      .value("GEOMETRY", EleAttribute::GEOMETRY)
      .value("COUPLER_AT", EleAttribute::COUPLER_AT)
      .value("E_TOT_OFFSET", EleAttribute::E_TOT_OFFSET)
      .value("PTC_CANONICAL_COORDS", EleAttribute::PTC_CANONICAL_COORDS)
      .value("B3_GRADIENT", EleAttribute::B3_GRADIENT)
      .value("E3_GRADIENT", EleAttribute::E3_GRADIENT)
      .value("PTC_FRINGE_GEOMETRY", EleAttribute::PTC_FRINGE_GEOMETRY)
      .value("E_TOT_SET", EleAttribute::E_TOT_SET)
      .value("BS_FIELD", EleAttribute::BS_FIELD)
      .value("P0C_SET", EleAttribute::P0C_SET)
      .value("PTC_FIELD_GEOMETRY", EleAttribute::PTC_FIELD_GEOMETRY)
      .value("DELTA_REF_TIME_USER_SET", EleAttribute::DELTA_REF_TIME_USER_SET)
      .value("DELTA_REF_TIME", EleAttribute::DELTA_REF_TIME)
      .value("P0C_START", EleAttribute::P0C_START)
      .value("E_TOT_START", EleAttribute::E_TOT_START)
      .value("P0C", EleAttribute::P0C)
      .value("E_TOT", EleAttribute::E_TOT)
      .value("X_PITCH_TOT", EleAttribute::X_PITCH_TOT)
      .value("NO_END_MARKER", EleAttribute::NO_END_MARKER)
      .value("Y_PITCH_TOT", EleAttribute::Y_PITCH_TOT)
      .value("X_OFFSET_TOT", EleAttribute::X_OFFSET_TOT)
      .value("Y_OFFSET_TOT", EleAttribute::Y_OFFSET_TOT)
      .value("Z_OFFSET_TOT", EleAttribute::Z_OFFSET_TOT)
      .value(
          "TILT_TOT",
          EleAttribute::TILT_TOT,
          "Important: tilt_tot$ = roll_tot$")
      .value("ROLL_TOT", EleAttribute::ROLL_TOT)
      .value("REF_TILT_TOT", EleAttribute::REF_TILT_TOT)
      .value("MULTIPASS_REF_ENERGY", EleAttribute::MULTIPASS_REF_ENERGY)
      .value("DISPATCH", EleAttribute::DISPATCH)
      .value("REF_TIME_START", EleAttribute::REF_TIME_START)
      .value(
          "THICKNESS",
          EleAttribute::THICKNESS,
          "For Etiennes' PTC: 2, 4, 6, or 8.")
      .value("INTEGRATOR_ORDER", EleAttribute::INTEGRATOR_ORDER)
      .value(
          "NUM_STEPS",
          EleAttribute::NUM_STEPS,
          "Assumed unique by set_flags_for_changed_real_attribute")
      .value(
          "DS_STEP",
          EleAttribute::DS_STEP,
          "Assumed unique by set_flags_for_changed_real_attribute")
      .value("CSR_DS_STEP", EleAttribute::CSR_DS_STEP)
      .value("LORD_PAD1", EleAttribute::LORD_PAD1)
      .value("LORD_PAD2", EleAttribute::LORD_PAD2)
      .value("REF_WAVELENGTH", EleAttribute::REF_WAVELENGTH)
      .value("X1_LIMIT", EleAttribute::X1_LIMIT)
      .value("X2_LIMIT", EleAttribute::X2_LIMIT)
      .value("Y1_LIMIT", EleAttribute::Y1_LIMIT)
      .value("Y2_LIMIT", EleAttribute::Y2_LIMIT)
      .value("CHECK_SUM", EleAttribute::CHECK_SUM)
      .export_values()
      .finalize();
  py::native_enum<EleKey>(m, "EleKey", "enum.IntEnum")
      .value("DRIFT", EleKey::DRIFT)
      .value("SBEND", EleKey::SBEND)
      .value("QUADRUPOLE", EleKey::QUADRUPOLE)
      .value("GROUP", EleKey::GROUP)
      .value("SEXTUPOLE", EleKey::SEXTUPOLE)
      .value("OVERLAY", EleKey::OVERLAY)
      .value("CUSTOM", EleKey::CUSTOM)
      .value("TAYLOR", EleKey::TAYLOR)
      .value("RFCAVITY", EleKey::RFCAVITY)
      .value("ELSEPARATOR", EleKey::ELSEPARATOR)
      .value("BEAMBEAM", EleKey::BEAMBEAM)
      .value("WIGGLER", EleKey::WIGGLER)
      .value("SOL_QUAD", EleKey::SOL_QUAD)
      .value("MARKER", EleKey::MARKER)
      .value("KICKER", EleKey::KICKER)
      .value("HYBRID", EleKey::HYBRID)
      .value("OCTUPOLE", EleKey::OCTUPOLE)
      .value("RBEND", EleKey::RBEND)
      .value("MULTIPOLE", EleKey::MULTIPOLE)
      .value("DEF_BMAD_COM", EleKey::DEF_BMAD_COM)
      .value("DEF_MAD_BEAM", EleKey::DEF_MAD_BEAM)
      .value("AB_MULTIPOLE", EleKey::AB_MULTIPOLE)
      .value("SOLENOID", EleKey::SOLENOID)
      .value("PATCH", EleKey::PATCH)
      .value("LCAVITY", EleKey::LCAVITY)
      .value("DEF_PARAMETER", EleKey::DEF_PARAMETER)
      .value("NULL_ELE", EleKey::NULL_ELE)
      .value("BEGINNING_ELE", EleKey::BEGINNING_ELE)
      .value("DEF_LINE", EleKey::DEF_LINE)
      .value("MATCH", EleKey::MATCH)
      .value("MONITOR", EleKey::MONITOR)
      .value("INSTRUMENT", EleKey::INSTRUMENT)
      .value("HKICKER", EleKey::HKICKER)
      .value("VKICKER", EleKey::VKICKER)
      .value("RCOLLIMATOR", EleKey::RCOLLIMATOR)
      .value("ECOLLIMATOR", EleKey::ECOLLIMATOR)
      .value("GIRDER", EleKey::GIRDER)
      .value("CONVERTER", EleKey::CONVERTER)
      .value("DEF_PARTICLE_START", EleKey::DEF_PARTICLE_START)
      .value("PHOTON_FORK", EleKey::PHOTON_FORK)
      .value("FORK", EleKey::FORK)
      .value("MIRROR", EleKey::MIRROR)
      .value("CRYSTAL", EleKey::CRYSTAL)
      .value("PIPE", EleKey::PIPE)
      .value("CAPILLARY", EleKey::CAPILLARY)
      .value("MULTILAYER_MIRROR", EleKey::MULTILAYER_MIRROR)
      .value("E_GUN", EleKey::E_GUN)
      .value("EM_FIELD", EleKey::EM_FIELD)
      .value("FLOOR_SHIFT", EleKey::FLOOR_SHIFT)
      .value("FIDUCIAL", EleKey::FIDUCIAL)
      .value("UNDULATOR", EleKey::UNDULATOR)
      .value("DIFFRACTION_PLATE", EleKey::DIFFRACTION_PLATE)
      .value("PHOTON_INIT", EleKey::PHOTON_INIT)
      .value("SAMPLE", EleKey::SAMPLE)
      .value("DETECTOR", EleKey::DETECTOR)
      .value("SAD_MULT", EleKey::SAD_MULT)
      .value("MASK", EleKey::MASK)
      .value("AC_KICKER", EleKey::AC_KICKER)
      .value("LENS", EleKey::LENS)
      .value("DEF_SPACE_CHARGE_COM", EleKey::DEF_SPACE_CHARGE_COM)
      .value("CRAB_CAVITY", EleKey::CRAB_CAVITY)
      .value("RAMPER", EleKey::RAMPER)
      .value("DEF_PTC_COM", EleKey::DEF_PTC_COM)
      .value("RF_BEND", EleKey::RF_BEND)
      .value("GKICKER", EleKey::GKICKER)
      .value("FOIL", EleKey::FOIL)
      .value("THICK_MULTIPOLE", EleKey::THICK_MULTIPOLE)
      .value("PICKUP", EleKey::PICKUP)
      .value("FEEDBACK", EleKey::FEEDBACK)
      .value("FIXER", EleKey::FIXER)
      .value("N_KEY", EleKey::N_KEY)
      .export_values()
      .finalize();

  // Enums from bmad_struct.f90
  m.attr("BMAD_INC_VERSION") = py::int_(Bmad::BMAD_INC_VERSION);
  m.attr("NONE") = py::int_(Bmad::NONE);
  // maximum multipole order
  m.attr("N_POLE_MAXX") = py::int_(Bmad::N_POLE_MAXX);
  // For indexing into ele%control%var(:) array
  m.attr("OLD_CONTROL_VAR_OFFSET") = py::int_(Bmad::OLD_CONTROL_VAR_OFFSET);
  // Important: var_offset$ > old_control_var_offset$
  m.attr("VAR_OFFSET") = py::int_(Bmad::VAR_OFFSET);
  // Maximum number of variables per controller.
  m.attr("N_VAR_MAX") = py::int_(Bmad::N_VAR_MAX);
  // Taylor term index offset.
  m.attr("TAYLOR_OFFSET") = py::int_(Bmad::TAYLOR_OFFSET);
  m.attr("BMAD_STANDARD") = py::int_(Bmad::BMAD_STANDARD);
  m.attr("SYMP_LIE_PTC") = py::int_(Bmad::SYMP_LIE_PTC);
  m.attr("RUNGE_KUTTA") = py::int_(Bmad::RUNGE_KUTTA);
  m.attr("LINEAR") = py::int_(Bmad::LINEAR);
  m.attr("TRACKING") = py::int_(Bmad::TRACKING);
  m.attr("TIME_RUNGE_KUTTA") = py::int_(Bmad::TIME_RUNGE_KUTTA);
  m.attr("FIXED_STEP_RUNGE_KUTTA") = py::int_(Bmad::FIXED_STEP_RUNGE_KUTTA);
  m.attr("SYMP_LIE_BMAD") = py::int_(Bmad::SYMP_LIE_BMAD);
  m.attr("MAGNUS") = py::int_(Bmad::MAGNUS);
  m.attr("AUTO") = py::int_(Bmad::AUTO);
  m.attr("SPRINT") = py::int_(Bmad::SPRINT);
  m.attr("FIXED_STEP_TIME_RUNGE_KUTTA") =
      py::int_(Bmad::FIXED_STEP_TIME_RUNGE_KUTTA);
  m.attr("MAD") = py::int_(Bmad::MAD);
  m.attr("TRANSVERSE_KICK") = py::int_(Bmad::TRANSVERSE_KICK);
  m.attr("SPIN_INTEGRATION") = py::int_(Bmad::SPIN_INTEGRATION);
  m.attr("DRIFT_KICK") = py::int_(Bmad::DRIFT_KICK);
  m.attr("MATRIX_KICK") = py::int_(Bmad::MATRIX_KICK);
  m.attr("RIPKEN_KICK") = py::int_(Bmad::RIPKEN_KICK);
  m.attr("SECTOR") = py::int_(Bmad::SECTOR);
  m.attr("STRAIGHT") = py::int_(Bmad::STRAIGHT);
  m.attr("FIELDMAP") = py::int_(Bmad::FIELDMAP);
  m.attr("PLANAR_MODEL") = py::int_(Bmad::PLANAR_MODEL);
  m.attr("REFER_TO_LORDS") = py::int_(Bmad::REFER_TO_LORDS);
  m.attr("NO_FIELD") = py::int_(Bmad::NO_FIELD);
  m.attr("HELICAL_MODEL") = py::int_(Bmad::HELICAL_MODEL);
  m.attr("SOFT_EDGE") = py::int_(Bmad::SOFT_EDGE);
  m.attr("UNIFORM") = py::int_(Bmad::UNIFORM);
  m.attr("GAUSSIAN") = py::int_(Bmad::GAUSSIAN);
  m.attr("SPHERICAL") = py::int_(Bmad::SPHERICAL);
  m.attr("CURVE") = py::int_(Bmad::CURVE);
  // Index to set slice_slave%ix_ele to.
  m.attr("IX_SLICE_SLAVE") = py::int_(Bmad::IX_SLICE_SLAVE);
  m.attr("MINOR_SLAVE") = py::int_(Bmad::MINOR_SLAVE);
  m.attr("SUPER_SLAVE") = py::int_(Bmad::SUPER_SLAVE);
  m.attr("FREE") = py::int_(Bmad::FREE);
  m.attr("GROUP_LORD") = py::int_(Bmad::GROUP_LORD);
  m.attr("SUPER_LORD") = py::int_(Bmad::SUPER_LORD);
  m.attr("OVERLAY_LORD") = py::int_(Bmad::OVERLAY_LORD);
  m.attr("GIRDER_LORD") = py::int_(Bmad::GIRDER_LORD);
  m.attr("MULTIPASS_LORD") = py::int_(Bmad::MULTIPASS_LORD);
  m.attr("MULTIPASS_SLAVE") = py::int_(Bmad::MULTIPASS_SLAVE);
  m.attr("NOT_A_LORD") = py::int_(Bmad::NOT_A_LORD);
  m.attr("SLICE_SLAVE") = py::int_(Bmad::SLICE_SLAVE);
  m.attr("CONTROL_LORD") = py::int_(Bmad::CONTROL_LORD);
  m.attr("RAMPER_LORD") = py::int_(Bmad::RAMPER_LORD);
  // governor$ = Union of overlay and group lords.
  m.attr("GOVERNOR") = py::int_(Bmad::GOVERNOR);
  m.attr("FIELD_LORD") = py::int_(Bmad::FIELD_LORD);
  // Used with pointer_to_lord(...)
  m.attr("MULTIPOLE_SOURCE") = py::int_(Bmad::MULTIPOLE_SOURCE);
  m.attr("AUTO_APERTURE") = py::int_(Bmad::AUTO_APERTURE);
  m.attr("RECTANGULAR") = py::int_(Bmad::RECTANGULAR);
  m.attr("ELLIPTICAL") = py::int_(Bmad::ELLIPTICAL);
  m.attr("WALL3D") = py::int_(Bmad::WALL3D);
  m.attr("CUSTOM_APERTURE") = py::int_(Bmad::CUSTOM_APERTURE);
  m.attr("LORD_DEFINED") = py::int_(Bmad::LORD_DEFINED);
  m.attr("SOFT_EDGE_ONLY") = py::int_(Bmad::SOFT_EDGE_ONLY);
  m.attr("HARD_EDGE_ONLY") = py::int_(Bmad::HARD_EDGE_ONLY);
  m.attr("FULL") = py::int_(Bmad::FULL);
  m.attr("SAD_FULL") = py::int_(Bmad::SAD_FULL);
  m.attr("LINEAR_EDGE") = py::int_(Bmad::LINEAR_EDGE);
  m.attr("BASIC_BEND") = py::int_(Bmad::BASIC_BEND);
  m.attr("STANDING_WAVE") = py::int_(Bmad::STANDING_WAVE);
  m.attr("TRAVELING_WAVE") = py::int_(Bmad::TRAVELING_WAVE);
  m.attr("PTC_STANDARD") = py::int_(Bmad::PTC_STANDARD);
  m.attr("X_INVARIANT") = py::int_(Bmad::X_INVARIANT);
  m.attr("MULTIPOLE_SYMMETRY") = py::int_(Bmad::MULTIPOLE_SYMMETRY);
  m.attr("CONTROL_VAR") = py::int_(Bmad::CONTROL_VAR);
  m.attr("OLD_CONTROL_VAR") = py::int_(Bmad::OLD_CONTROL_VAR);
  m.attr("ALL_CONTROL_VAR") = py::int_(Bmad::ALL_CONTROL_VAR);
  m.attr("ELEC_MULTIPOLE") = py::int_(Bmad::ELEC_MULTIPOLE);
  m.attr("OK") = py::int_(Bmad::OK);
  m.attr("IN_STOP_BAND") = py::int_(Bmad::IN_STOP_BAND);
  m.attr("NON_SYMPLECTIC") = py::int_(Bmad::NON_SYMPLECTIC);
  m.attr("UNSTABLE") = py::int_(Bmad::UNSTABLE);
  m.attr("UNSTABLE_A") = py::int_(Bmad::UNSTABLE_A);
  m.attr("UNSTABLE_B") = py::int_(Bmad::UNSTABLE_B);
  m.attr("XFER_MAT_CALC_FAILURE") = py::int_(Bmad::XFER_MAT_CALC_FAILURE);
  m.attr("TWISS_PROPAGATE_FAILURE") = py::int_(Bmad::TWISS_PROPAGATE_FAILURE);
  m.attr("NO_CLOSED_ORBIT") = py::int_(Bmad::NO_CLOSED_ORBIT);
  m.attr("NO_COMPLETE_ORBIT") = py::int_(Bmad::NO_COMPLETE_ORBIT);
  m.attr("INCLUDE_KICKS") = py::int_(Bmad::INCLUDE_KICKS);
  m.attr("SHORT") = py::int_(Bmad::SHORT);
  m.attr("USER_SET") = py::int_(Bmad::USER_SET);
  m.attr("FIRST_PASS") = py::int_(Bmad::FIRST_PASS);
  m.attr("HIGHLAND") = py::int_(Bmad::HIGHLAND);
  m.attr("LYNCH_DAHL") = py::int_(Bmad::LYNCH_DAHL);
  m.attr("INCOHERENT") = py::int_(Bmad::INCOHERENT);
  m.attr("COHERENT") = py::int_(Bmad::COHERENT);
  m.attr("ASCII") = py::int_(Bmad::ASCII);
  m.attr("BINARY") = py::int_(Bmad::BINARY);
  m.attr("HDF5") = py::int_(Bmad::HDF5);
  m.attr("ONE_FILE") = py::int_(Bmad::ONE_FILE);
  // For testing purposes.
  m.attr("OLD_ASCII") = py::int_(Bmad::OLD_ASCII);
  m.attr("NUM_ELE_ATTRIB") = py::int_(Bmad::NUM_ELE_ATTRIB);
  m.attr("OFF") = py::int_(Bmad::OFF);
  m.attr("ON") = py::int_(Bmad::ON);
  m.attr("SAVE_STATE") = py::int_(Bmad::SAVE_STATE);
  m.attr("RESTORE_STATE") = py::int_(Bmad::RESTORE_STATE);
  m.attr("OFF_AND_SAVE") = py::int_(Bmad::OFF_AND_SAVE);
  m.attr("HORIZONTALLY_PURE") = py::int_(Bmad::HORIZONTALLY_PURE);
  m.attr("VERTICALLY_PURE") = py::int_(Bmad::VERTICALLY_PURE);
  m.attr("ONE_DIM") = py::int_(Bmad::ONE_DIM);
  m.attr("STEADY_STATE_3D") = py::int_(Bmad::STEADY_STATE_3D);
  m.attr("SLICE") = py::int_(Bmad::SLICE);
  m.attr("FFT_3D") = py::int_(Bmad::FFT_3D);
  m.attr("CATHODE_FFT_3D") = py::int_(Bmad::CATHODE_FFT_3D);
  m.attr("MAGNETIC") = py::int_(Bmad::MAGNETIC);
  m.attr("ELECTRIC") = py::int_(Bmad::ELECTRIC);
  m.attr("MIXED") = py::int_(Bmad::MIXED);
  m.attr("BRAGG_DIFFRACTED") = py::int_(Bmad::BRAGG_DIFFRACTED);
  m.attr("FORWARD_DIFFRACTED") = py::int_(Bmad::FORWARD_DIFFRACTED);
  m.attr("UNDIFFRACTED") = py::int_(Bmad::UNDIFFRACTED);
  m.attr("REFLECTION") = py::int_(Bmad::REFLECTION);
  m.attr("TRANSMISSION") = py::int_(Bmad::TRANSMISSION);
  m.attr("ANCHOR_BEGINNING") = py::int_(Bmad::ANCHOR_BEGINNING);
  m.attr("ANCHOR_CENTER") = py::int_(Bmad::ANCHOR_CENTER);
  m.attr("ANCHOR_END") = py::int_(Bmad::ANCHOR_END);
  m.attr("NONE_PT") = py::int_(Bmad::NONE_PT);
  m.attr("ENTRANCE_END") = py::int_(Bmad::ENTRANCE_END);
  m.attr("EXIT_END") = py::int_(Bmad::EXIT_END);
  m.attr("BOTH_ENDS") = py::int_(Bmad::BOTH_ENDS);
  m.attr("NO_END") = py::int_(Bmad::NO_END);
  m.attr("NO_APERTURE") = py::int_(Bmad::NO_APERTURE);
  m.attr("NOWHERE") = py::int_(Bmad::NOWHERE);
  m.attr("CONTINUOUS") = py::int_(Bmad::CONTINUOUS);
  m.attr("SURFACE") = py::int_(Bmad::SURFACE);
  m.attr("WALL_TRANSITION") = py::int_(Bmad::WALL_TRANSITION);
  m.attr("UPSTREAM_END") = py::int_(Bmad::UPSTREAM_END);
  m.attr("DOWNSTREAM_END") = py::int_(Bmad::DOWNSTREAM_END);
  m.attr("INSIDE") = py::int_(Bmad::INSIDE);
  m.attr("CENTER_PT") = py::int_(Bmad::CENTER_PT);
  m.attr("START_END") = py::int_(Bmad::START_END);
  // Must be different from upstream_end$, downstream_end$
  m.attr("FIRST_TRACK_EDGE") = py::int_(Bmad::FIRST_TRACK_EDGE);
  m.attr("SECOND_TRACK_EDGE") = py::int_(Bmad::SECOND_TRACK_EDGE);
  m.attr("IN_BETWEEN") = py::int_(Bmad::IN_BETWEEN);
  m.attr("NORMAL") = py::int_(Bmad::NORMAL);
  m.attr("CLEAR") = py::int_(Bmad::CLEAR);
  m.attr("OPAQUE") = py::int_(Bmad::OPAQUE);
  m.attr("WALL_START") = py::int_(Bmad::WALL_START);
  m.attr("WALL_END") = py::int_(Bmad::WALL_END);
  m.attr("ABSOLUTE") = py::int_(Bmad::ABSOLUTE);
  m.attr("RELATIVE") = py::int_(Bmad::RELATIVE);
  m.attr("SHIFTED_TO_RELATIVE") = py::int_(Bmad::SHIFTED_TO_RELATIVE);
  m.attr("CHAMBER_WALL") = py::int_(Bmad::CHAMBER_WALL);
  m.attr("MASK_PLATE") = py::int_(Bmad::MASK_PLATE);
  m.attr("X_PLANE") = py::int_(Bmad::X_PLANE);
  m.attr("Y_PLANE") = py::int_(Bmad::Y_PLANE);
  m.attr("Z_PLANE") = py::int_(Bmad::Z_PLANE);
  m.attr("N_PLANE") = py::int_(Bmad::N_PLANE);
  m.attr("S_PLANE") = py::int_(Bmad::S_PLANE);
  m.attr("MOVING_FORWARD") = py::int_(Bmad::MOVING_FORWARD);
  // EG: before cathode emission. Conforms to OpenPMD standard.
  m.attr("PRE_BORN") = py::int_(Bmad::PRE_BORN);
  // Conforms to OpenPMD standard.
  m.attr("ALIVE") = py::int_(Bmad::ALIVE);
  m.attr("LOST") = py::int_(Bmad::LOST);
  m.attr("LOST_NEG_X") = py::int_(Bmad::LOST_NEG_X);
  m.attr("LOST_POS_X") = py::int_(Bmad::LOST_POS_X);
  m.attr("LOST_NEG_Y") = py::int_(Bmad::LOST_NEG_Y);
  m.attr("LOST_POS_Y") = py::int_(Bmad::LOST_POS_Y);
  m.attr("LOST_Z") = py::int_(Bmad::LOST_Z);
  // Particle "turned around" when not tracking with time_runge_kutta.
  m.attr("LOST_PZ") = py::int_(Bmad::LOST_PZ);
  // old names.
  m.attr("LOST_NEG_X_APERTURE") = py::int_(Bmad::LOST_NEG_X_APERTURE);
  m.attr("LOST_POS_X_APERTURE") = py::int_(Bmad::LOST_POS_X_APERTURE);
  m.attr("LOST_NEG_Y_APERTURE") = py::int_(Bmad::LOST_NEG_Y_APERTURE);
  m.attr("LOST_POS_Y_APERTURE") = py::int_(Bmad::LOST_POS_Y_APERTURE);
  m.attr("LOST_Z_APERTURE") = py::int_(Bmad::LOST_Z_APERTURE);
  // Particle "turned around" when not tracking with time_runge_kutta.
  m.attr("LOST_PZ_APERTURE") = py::int_(Bmad::LOST_PZ_APERTURE);
  m.attr("NO_MISALIGNMENT") = py::float_(Bmad::NO_MISALIGNMENT);
  m.attr("X_POLARIZATION") = py::int_(Bmad::X_POLARIZATION);
  m.attr("Y_POLARIZATION") = py::int_(Bmad::Y_POLARIZATION);
  m.attr("XY") = py::int_(Bmad::XY);
  m.attr("LEADING") = py::int_(Bmad::LEADING);
  m.attr("TRAILING") = py::int_(Bmad::TRAILING);
  m.attr("X_LEADING") = py::int_(Bmad::X_LEADING);
  m.attr("Y_LEADING") = py::int_(Bmad::Y_LEADING);
  m.attr("X_TRAILING") = py::int_(Bmad::X_TRAILING);
  m.attr("Y_TRAILING") = py::int_(Bmad::Y_TRAILING);
  m.attr("FAMILY_Y") = py::int_(Bmad::FAMILY_Y);
  m.attr("FAMILY_X") = py::int_(Bmad::FAMILY_X);
  m.attr("FAMILY_QU") = py::int_(Bmad::FAMILY_QU);
  m.attr("FAMILY_SQ") = py::int_(Bmad::FAMILY_SQ);
  m.attr("HYPER_Y") = py::int_(Bmad::HYPER_Y);
  m.attr("HYPER_XY") = py::int_(Bmad::HYPER_XY);
  m.attr("HYPER_X") = py::int_(Bmad::HYPER_X);
  m.attr("SUPER_OK") = py::int_(Bmad::SUPER_OK);
  m.attr("STALE") = py::int_(Bmad::STALE);
  m.attr("ATTRIBUTE_GROUP") = py::int_(Bmad::ATTRIBUTE_GROUP);
  m.attr("CONTROL_GROUP") = py::int_(Bmad::CONTROL_GROUP);
  m.attr("FLOOR_POSITION_GROUP") = py::int_(Bmad::FLOOR_POSITION_GROUP);
  m.attr("S_POSITION_GROUP") = py::int_(Bmad::S_POSITION_GROUP);
  m.attr("REF_ENERGY_GROUP") = py::int_(Bmad::REF_ENERGY_GROUP);
  m.attr("MAT6_GROUP") = py::int_(Bmad::MAT6_GROUP);
  m.attr("RAD_INT_GROUP") = py::int_(Bmad::RAD_INT_GROUP);
  m.attr("ALL_GROUPS") = py::int_(Bmad::ALL_GROUPS);
  m.attr("S_AND_FLOOR_POSITION_GROUP") =
      py::int_(Bmad::S_AND_FLOOR_POSITION_GROUP);
  m.attr("POLARIZED") = py::int_(Bmad::POLARIZED);
  m.attr("UNPOLARIZED") = py::int_(Bmad::UNPOLARIZED);
  m.attr("CUBIC") = py::int_(Bmad::CUBIC);
  m.attr("OPAL") = py::int_(Bmad::OPAL);
  m.attr("IMPACTT") = py::int_(Bmad::IMPACTT);
  m.attr("DRIFT") = py::int_(Bmad::DRIFT);
  m.attr("SBEND") = py::int_(Bmad::SBEND);
  m.attr("QUADRUPOLE") = py::int_(Bmad::QUADRUPOLE);
  m.attr("GROUP") = py::int_(Bmad::GROUP);
  m.attr("SEXTUPOLE") = py::int_(Bmad::SEXTUPOLE);
  m.attr("OVERLAY") = py::int_(Bmad::OVERLAY);
  m.attr("CUSTOM") = py::int_(Bmad::CUSTOM);
  m.attr("TAYLOR") = py::int_(Bmad::TAYLOR);
  m.attr("RFCAVITY") = py::int_(Bmad::RFCAVITY);
  m.attr("ELSEPARATOR") = py::int_(Bmad::ELSEPARATOR);
  m.attr("BEAMBEAM") = py::int_(Bmad::BEAMBEAM);
  m.attr("WIGGLER") = py::int_(Bmad::WIGGLER);
  m.attr("SOL_QUAD") = py::int_(Bmad::SOL_QUAD);
  m.attr("MARKER") = py::int_(Bmad::MARKER);
  m.attr("KICKER") = py::int_(Bmad::KICKER);
  m.attr("HYBRID") = py::int_(Bmad::HYBRID);
  m.attr("OCTUPOLE") = py::int_(Bmad::OCTUPOLE);
  m.attr("RBEND") = py::int_(Bmad::RBEND);
  m.attr("MULTIPOLE") = py::int_(Bmad::MULTIPOLE);
  m.attr("DEF_BMAD_COM") = py::int_(Bmad::DEF_BMAD_COM);
  m.attr("DEF_MAD_BEAM") = py::int_(Bmad::DEF_MAD_BEAM);
  m.attr("AB_MULTIPOLE") = py::int_(Bmad::AB_MULTIPOLE);
  m.attr("SOLENOID") = py::int_(Bmad::SOLENOID);
  m.attr("PATCH") = py::int_(Bmad::PATCH);
  m.attr("LCAVITY") = py::int_(Bmad::LCAVITY);
  m.attr("DEF_PARAMETER") = py::int_(Bmad::DEF_PARAMETER);
  m.attr("NULL_ELE") = py::int_(Bmad::NULL_ELE);
  m.attr("BEGINNING_ELE") = py::int_(Bmad::BEGINNING_ELE);
  m.attr("DEF_LINE") = py::int_(Bmad::DEF_LINE);
  m.attr("MATCH") = py::int_(Bmad::MATCH);
  m.attr("MONITOR") = py::int_(Bmad::MONITOR);
  m.attr("INSTRUMENT") = py::int_(Bmad::INSTRUMENT);
  m.attr("HKICKER") = py::int_(Bmad::HKICKER);
  m.attr("VKICKER") = py::int_(Bmad::VKICKER);
  m.attr("RCOLLIMATOR") = py::int_(Bmad::RCOLLIMATOR);
  m.attr("ECOLLIMATOR") = py::int_(Bmad::ECOLLIMATOR);
  m.attr("GIRDER") = py::int_(Bmad::GIRDER);
  m.attr("CONVERTER") = py::int_(Bmad::CONVERTER);
  m.attr("DEF_PARTICLE_START") = py::int_(Bmad::DEF_PARTICLE_START);
  m.attr("PHOTON_FORK") = py::int_(Bmad::PHOTON_FORK);
  m.attr("FORK") = py::int_(Bmad::FORK);
  m.attr("MIRROR") = py::int_(Bmad::MIRROR);
  m.attr("CRYSTAL") = py::int_(Bmad::CRYSTAL);
  m.attr("PIPE") = py::int_(Bmad::PIPE);
  m.attr("CAPILLARY") = py::int_(Bmad::CAPILLARY);
  m.attr("MULTILAYER_MIRROR") = py::int_(Bmad::MULTILAYER_MIRROR);
  m.attr("E_GUN") = py::int_(Bmad::E_GUN);
  m.attr("EM_FIELD") = py::int_(Bmad::EM_FIELD);
  m.attr("FLOOR_SHIFT") = py::int_(Bmad::FLOOR_SHIFT);
  m.attr("FIDUCIAL") = py::int_(Bmad::FIDUCIAL);
  m.attr("UNDULATOR") = py::int_(Bmad::UNDULATOR);
  m.attr("DIFFRACTION_PLATE") = py::int_(Bmad::DIFFRACTION_PLATE);
  m.attr("PHOTON_INIT") = py::int_(Bmad::PHOTON_INIT);
  m.attr("SAMPLE") = py::int_(Bmad::SAMPLE);
  m.attr("DETECTOR") = py::int_(Bmad::DETECTOR);
  m.attr("SAD_MULT") = py::int_(Bmad::SAD_MULT);
  m.attr("MASK") = py::int_(Bmad::MASK);
  m.attr("AC_KICKER") = py::int_(Bmad::AC_KICKER);
  m.attr("LENS") = py::int_(Bmad::LENS);
  m.attr("DEF_SPACE_CHARGE_COM") = py::int_(Bmad::DEF_SPACE_CHARGE_COM);
  m.attr("CRAB_CAVITY") = py::int_(Bmad::CRAB_CAVITY);
  m.attr("RAMPER") = py::int_(Bmad::RAMPER);
  m.attr("DEF_PTC_COM") = py::int_(Bmad::DEF_PTC_COM);
  m.attr("RF_BEND") = py::int_(Bmad::RF_BEND);
  m.attr("GKICKER") = py::int_(Bmad::GKICKER);
  m.attr("FOIL") = py::int_(Bmad::FOIL);
  m.attr("THICK_MULTIPOLE") = py::int_(Bmad::THICK_MULTIPOLE);
  m.attr("PICKUP") = py::int_(Bmad::PICKUP);
  m.attr("FEEDBACK") = py::int_(Bmad::FEEDBACK);
  m.attr("FIXER") = py::int_(Bmad::FIXER);
  m.attr("N_KEY") = py::int_(Bmad::N_KEY);
  m.attr("STANDARD") = py::int_(Bmad::STANDARD);
  m.attr("MATCH_TWISS") = py::int_(Bmad::MATCH_TWISS);
  m.attr("IDENTITY") = py::int_(Bmad::IDENTITY);
  m.attr("PHASE_TROMBONE") = py::int_(Bmad::PHASE_TROMBONE);
  m.attr("MATCH_ORBIT") = py::int_(Bmad::MATCH_ORBIT);
  m.attr("ZERO") = py::int_(Bmad::ZERO);
  m.attr("VAL1") = py::int_(Bmad::VAL1);
  m.attr("VAL2") = py::int_(Bmad::VAL2);
  m.attr("VAL3") = py::int_(Bmad::VAL3);
  m.attr("VAL4") = py::int_(Bmad::VAL4);
  m.attr("VAL5") = py::int_(Bmad::VAL5);
  m.attr("VAL6") = py::int_(Bmad::VAL6);
  m.attr("VAL7") = py::int_(Bmad::VAL7);
  m.attr("VAL8") = py::int_(Bmad::VAL8);
  m.attr("VAL9") = py::int_(Bmad::VAL9);
  m.attr("VAL10") = py::int_(Bmad::VAL10);
  m.attr("VAL11") = py::int_(Bmad::VAL11);
  m.attr("VAL12") = py::int_(Bmad::VAL12);
  m.attr("BETA_A0") = py::int_(Bmad::BETA_A0);
  m.attr("ALPHA_A0") = py::int_(Bmad::ALPHA_A0);
  m.attr("BETA_B0") = py::int_(Bmad::BETA_B0);
  m.attr("ALPHA_B0") = py::int_(Bmad::ALPHA_B0);
  m.attr("BETA_A1") = py::int_(Bmad::BETA_A1);
  m.attr("ALPHA_A1") = py::int_(Bmad::ALPHA_A1);
  m.attr("BETA_B1") = py::int_(Bmad::BETA_B1);
  m.attr("ALPHA_B1") = py::int_(Bmad::ALPHA_B1);
  m.attr("DPHI_A") = py::int_(Bmad::DPHI_A);
  m.attr("DPHI_B") = py::int_(Bmad::DPHI_B);
  m.attr("ETA_X0") = py::int_(Bmad::ETA_X0);
  m.attr("ETAP_X0") = py::int_(Bmad::ETAP_X0);
  m.attr("ETA_Y0") = py::int_(Bmad::ETA_Y0);
  m.attr("ETAP_Y0") = py::int_(Bmad::ETAP_Y0);
  m.attr("ETA_X1") = py::int_(Bmad::ETA_X1);
  m.attr("ETAP_X1") = py::int_(Bmad::ETAP_X1);
  m.attr("ETA_Y1") = py::int_(Bmad::ETA_Y1);
  m.attr("ETAP_Y1") = py::int_(Bmad::ETAP_Y1);
  m.attr("C11_MAT0") = py::int_(Bmad::C11_MAT0);
  m.attr("C12_MAT0") = py::int_(Bmad::C12_MAT0);
  m.attr("C21_MAT0") = py::int_(Bmad::C21_MAT0);
  m.attr("C22_MAT0") = py::int_(Bmad::C22_MAT0);
  m.attr("MODE_FLIP0") = py::int_(Bmad::MODE_FLIP0);
  m.attr("C11_MAT1") = py::int_(Bmad::C11_MAT1);
  m.attr("C12_MAT1") = py::int_(Bmad::C12_MAT1);
  m.attr("C21_MAT1") = py::int_(Bmad::C21_MAT1);
  m.attr("C22_MAT1") = py::int_(Bmad::C22_MAT1);
  m.attr("MODE_FLIP1") = py::int_(Bmad::MODE_FLIP1);
  m.attr("X0") = py::int_(Bmad::X0);
  m.attr("PX0") = py::int_(Bmad::PX0);
  m.attr("Y0") = py::int_(Bmad::Y0);
  m.attr("PY0") = py::int_(Bmad::PY0);
  m.attr("Z0") = py::int_(Bmad::Z0);
  m.attr("PZ0") = py::int_(Bmad::PZ0);
  m.attr("X1") = py::int_(Bmad::X1);
  m.attr("PX1") = py::int_(Bmad::PX1);
  m.attr("Y1") = py::int_(Bmad::Y1);
  m.attr("PY1") = py::int_(Bmad::PY1);
  m.attr("Z1") = py::int_(Bmad::Z1);
  m.attr("PZ1") = py::int_(Bmad::PZ1);
  m.attr("MATRIX") = py::int_(Bmad::MATRIX);
  m.attr("KICK0") = py::int_(Bmad::KICK0);
  m.attr("RECALC") = py::int_(Bmad::RECALC);
  m.attr("DELTA_TIME") = py::int_(Bmad::DELTA_TIME);
  m.attr("X") = py::int_(Bmad::X);
  m.attr("PX") = py::int_(Bmad::PX);
  m.attr("Y") = py::int_(Bmad::Y);
  m.attr("PY") = py::int_(Bmad::PY);
  m.attr("Z") = py::int_(Bmad::Z);
  m.attr("PZ") = py::int_(Bmad::PZ);
  m.attr("T") = py::int_(Bmad::T);
  m.attr("FIELD_X") = py::int_(Bmad::FIELD_X);
  m.attr("FIELD_Y") = py::int_(Bmad::FIELD_Y);
  m.attr("PHASE_X") = py::int_(Bmad::PHASE_X);
  m.attr("PHASE_Y") = py::int_(Bmad::PHASE_Y);
  m.attr("E_PHOTON") = py::int_(Bmad::E_PHOTON);
  m.attr("E1") = py::int_(Bmad::E1);
  m.attr("E2") = py::int_(Bmad::E2);
  m.attr("FINT") = py::int_(Bmad::FINT);
  m.attr("FINTX") = py::int_(Bmad::FINTX);
  m.attr("HGAP") = py::int_(Bmad::HGAP);
  m.attr("HGAPX") = py::int_(Bmad::HGAPX);
  m.attr("H1") = py::int_(Bmad::H1);
  m.attr("H2") = py::int_(Bmad::H2);
  m.attr("SPIN_X_STORED") = py::int_(Bmad::SPIN_X_STORED);
  m.attr("SPIN_Y_STORED") = py::int_(Bmad::SPIN_Y_STORED);
  m.attr("SPIN_Z_STORED") = py::int_(Bmad::SPIN_Z_STORED);
  m.attr("X_STORED") = py::int_(Bmad::X_STORED);
  m.attr("PX_STORED") = py::int_(Bmad::PX_STORED);
  m.attr("Y_STORED") = py::int_(Bmad::Y_STORED);
  m.attr("PY_STORED") = py::int_(Bmad::PY_STORED);
  m.attr("Z_STORED") = py::int_(Bmad::Z_STORED);
  m.attr("PZ_STORED") = py::int_(Bmad::PZ_STORED);
  m.attr("BETA_A_STORED") = py::int_(Bmad::BETA_A_STORED);
  m.attr("ALPHA_A_STORED") = py::int_(Bmad::ALPHA_A_STORED);
  m.attr("BETA_B_STORED") = py::int_(Bmad::BETA_B_STORED);
  m.attr("ALPHA_B_STORED") = py::int_(Bmad::ALPHA_B_STORED);
  m.attr("PHI_A_STORED") = py::int_(Bmad::PHI_A_STORED);
  m.attr("PHI_B_STORED") = py::int_(Bmad::PHI_B_STORED);
  m.attr("MODE_FLIP_STORED") = py::int_(Bmad::MODE_FLIP_STORED);
  m.attr("ETA_X_STORED") = py::int_(Bmad::ETA_X_STORED);
  m.attr("ETAP_X_STORED") = py::int_(Bmad::ETAP_X_STORED);
  m.attr("ETA_Y_STORED") = py::int_(Bmad::ETA_Y_STORED);
  m.attr("ETAP_Y_STORED") = py::int_(Bmad::ETAP_Y_STORED);
  m.attr("CMAT_11_STORED") = py::int_(Bmad::CMAT_11_STORED);
  m.attr("CMAT_12_STORED") = py::int_(Bmad::CMAT_12_STORED);
  m.attr("CMAT_21_STORED") = py::int_(Bmad::CMAT_21_STORED);
  m.attr("CMAT_22_STORED") = py::int_(Bmad::CMAT_22_STORED);
  m.attr("DBETA_DPZ_A_STORED") = py::int_(Bmad::DBETA_DPZ_A_STORED);
  m.attr("DBETA_DPZ_B_STORED") = py::int_(Bmad::DBETA_DPZ_B_STORED);
  m.attr("DALPHA_DPZ_A_STORED") = py::int_(Bmad::DALPHA_DPZ_A_STORED);
  m.attr("DALPHA_DPZ_B_STORED") = py::int_(Bmad::DALPHA_DPZ_B_STORED);
  m.attr("DETA_DPZ_X_STORED") = py::int_(Bmad::DETA_DPZ_X_STORED);
  m.attr("DETA_DPZ_Y_STORED") = py::int_(Bmad::DETA_DPZ_Y_STORED);
  m.attr("DETAP_DPZ_X_STORED") = py::int_(Bmad::DETAP_DPZ_X_STORED);
  m.attr("DETAP_DPZ_Y_STORED") = py::int_(Bmad::DETAP_DPZ_Y_STORED);
  m.attr("DCMAT_DPZ_11_STORED") = py::int_(Bmad::DCMAT_DPZ_11_STORED);
  m.attr("DCMAT_DPZ_12_STORED") = py::int_(Bmad::DCMAT_DPZ_12_STORED);
  m.attr("DCMAT_DPZ_21_STORED") = py::int_(Bmad::DCMAT_DPZ_21_STORED);
  m.attr("DCMAT_DPZ_22_STORED") = py::int_(Bmad::DCMAT_DPZ_22_STORED);
  m.attr("RADIUS") = py::int_(Bmad::RADIUS);
  m.attr("FOCAL_STRENGTH") = py::int_(Bmad::FOCAL_STRENGTH);
  // Assumed unique. Do not assign 1 to another attribute.
  m.attr("L") = py::int_(Bmad::L);
  // Important: tilt$ = roll$
  m.attr("TILT") = py::int_(Bmad::TILT);
  m.attr("ROLL") = py::int_(Bmad::ROLL);
  m.attr("N_PART") = py::int_(Bmad::N_PART);
  m.attr("INHERIT_FROM_FORK") = py::int_(Bmad::INHERIT_FROM_FORK);
  m.attr("REF_TILT") = py::int_(Bmad::REF_TILT);
  m.attr("DIRECTION") = py::int_(Bmad::DIRECTION);
  m.attr("REPETITION_FREQUENCY") = py::int_(Bmad::REPETITION_FREQUENCY);
  m.attr("DETA_DS_MASTER") = py::int_(Bmad::DETA_DS_MASTER);
  m.attr("KICK") = py::int_(Bmad::KICK);
  m.attr("X_GAIN_ERR") = py::int_(Bmad::X_GAIN_ERR);
  m.attr("TAYLOR_ORDER") = py::int_(Bmad::TAYLOR_ORDER);
  m.attr("R_SOLENOID") = py::int_(Bmad::R_SOLENOID);
  m.attr("FINAL_CHARGE") = py::int_(Bmad::FINAL_CHARGE);
  m.attr("K1") = py::int_(Bmad::K1);
  m.attr("KX") = py::int_(Bmad::KX);
  m.attr("HARMON") = py::int_(Bmad::HARMON);
  m.attr("H_DISPLACE") = py::int_(Bmad::H_DISPLACE);
  m.attr("Y_GAIN_ERR") = py::int_(Bmad::Y_GAIN_ERR);
  m.attr("S_TWISS_REF") = py::int_(Bmad::S_TWISS_REF);
  m.attr("CRITICAL_ANGLE_FACTOR") = py::int_(Bmad::CRITICAL_ANGLE_FACTOR);
  m.attr("TILT_CORR") = py::int_(Bmad::TILT_CORR);
  m.attr("REF_COORDS") = py::int_(Bmad::REF_COORDS);
  m.attr("DT_MAX") = py::int_(Bmad::DT_MAX);
  m.attr("IX_FIXER") = py::int_(Bmad::IX_FIXER);
  m.attr("GRAZE_ANGLE") = py::int_(Bmad::GRAZE_ANGLE);
  m.attr("K2") = py::int_(Bmad::K2);
  m.attr("B_MAX") = py::int_(Bmad::B_MAX);
  m.attr("V_DISPLACE") = py::int_(Bmad::V_DISPLACE);
  m.attr("GRADIENT_TOT") = py::int_(Bmad::GRADIENT_TOT);
  m.attr("HARMON_MASTER") = py::int_(Bmad::HARMON_MASTER);
  m.attr("FLEXIBLE") = py::int_(Bmad::FLEXIBLE);
  m.attr("CRUNCH") = py::int_(Bmad::CRUNCH);
  m.attr("REF_ORBIT_FOLLOWS") = py::int_(Bmad::REF_ORBIT_FOLLOWS);
  m.attr("PC_OUT_MIN") = py::int_(Bmad::PC_OUT_MIN);
  m.attr("GRADIENT") = py::int_(Bmad::GRADIENT);
  m.attr("K3") = py::int_(Bmad::K3);
  m.attr("NOISE") = py::int_(Bmad::NOISE);
  m.attr("NEW_BRANCH") = py::int_(Bmad::NEW_BRANCH);
  m.attr("IX_BRANCH") = py::int_(Bmad::IX_BRANCH);
  m.attr("G_MAX") = py::int_(Bmad::G_MAX);
  m.attr("G") = py::int_(Bmad::G);
  m.attr("SYMMETRY") = py::int_(Bmad::SYMMETRY);
  m.attr("FIELD_SCALE_FACTOR") = py::int_(Bmad::FIELD_SCALE_FACTOR);
  m.attr("PC_OUT_MAX") = py::int_(Bmad::PC_OUT_MAX);
  m.attr("DG") = py::int_(Bmad::DG);
  m.attr("BBI_CONST") = py::int_(Bmad::BBI_CONST);
  m.attr("OSC_AMPLITUDE") = py::int_(Bmad::OSC_AMPLITUDE);
  m.attr("IX_TO_BRANCH") = py::int_(Bmad::IX_TO_BRANCH);
  m.attr("ANGLE_OUT_MAX") = py::int_(Bmad::ANGLE_OUT_MAX);
  m.attr("GRADIENT_ERR") = py::int_(Bmad::GRADIENT_ERR);
  m.attr("CRITICAL_ANGLE") = py::int_(Bmad::CRITICAL_ANGLE);
  m.attr("BRAGG_ANGLE_IN") = py::int_(Bmad::BRAGG_ANGLE_IN);
  m.attr("SPIN_DN_DPZ_X") = py::int_(Bmad::SPIN_DN_DPZ_X);
  m.attr("DELTA_E_REF") = py::int_(Bmad::DELTA_E_REF);
  m.attr("INTERPOLATION") = py::int_(Bmad::INTERPOLATION);
  m.attr("BRAGG_ANGLE_OUT") = py::int_(Bmad::BRAGG_ANGLE_OUT);
  m.attr("K1X") = py::int_(Bmad::K1X);
  m.attr("SPIN_DN_DPZ_Y") = py::int_(Bmad::SPIN_DN_DPZ_Y);
  m.attr("CHARGE") = py::int_(Bmad::CHARGE);
  m.attr("X_GAIN_CALIB") = py::int_(Bmad::X_GAIN_CALIB);
  m.attr("IX_TO_ELEMENT") = py::int_(Bmad::IX_TO_ELEMENT);
  m.attr("VOLTAGE") = py::int_(Bmad::VOLTAGE);
  m.attr("G_TOT") = py::int_(Bmad::G_TOT);
  m.attr("RHO") = py::int_(Bmad::RHO);
  m.attr("VOLTAGE_ERR") = py::int_(Bmad::VOLTAGE_ERR);
  m.attr("BRAGG_ANGLE") = py::int_(Bmad::BRAGG_ANGLE);
  m.attr("K1Y") = py::int_(Bmad::K1Y);
  m.attr("N_PARTICLE") = py::int_(Bmad::N_PARTICLE);
  m.attr("SPIN_DN_DPZ_Z") = py::int_(Bmad::SPIN_DN_DPZ_Z);
  m.attr("FRINGE_TYPE") = py::int_(Bmad::FRINGE_TYPE);
  m.attr("DBRAGG_ANGLE_DE") = py::int_(Bmad::DBRAGG_ANGLE_DE);
  m.attr("FRINGE_AT") = py::int_(Bmad::FRINGE_AT);
  m.attr("GANG") = py::int_(Bmad::GANG);
  m.attr("DARWIN_WIDTH_SIGMA") = py::int_(Bmad::DARWIN_WIDTH_SIGMA);
  m.attr("DARWIN_WIDTH_PI") = py::int_(Bmad::DARWIN_WIDTH_PI);
  m.attr("SPIN_FRINGE_ON") = py::int_(Bmad::SPIN_FRINGE_ON);
  m.attr("PENDELLOSUNG_PERIOD_SIGMA") =
      py::int_(Bmad::PENDELLOSUNG_PERIOD_SIGMA);
  m.attr("SIG_X") = py::int_(Bmad::SIG_X);
  m.attr("EXACT_MULTIPOLES") = py::int_(Bmad::EXACT_MULTIPOLES);
  m.attr("PENDELLOSUNG_PERIOD_PI") = py::int_(Bmad::PENDELLOSUNG_PERIOD_PI);
  m.attr("SIG_Y") = py::int_(Bmad::SIG_Y);
  m.attr("GRAZE_ANGLE_IN") = py::int_(Bmad::GRAZE_ANGLE_IN);
  m.attr("R0_ELEC") = py::int_(Bmad::R0_ELEC);
  m.attr("RF_FREQUENCY") = py::int_(Bmad::RF_FREQUENCY);
  m.attr("SIG_Z") = py::int_(Bmad::SIG_Z);
  m.attr("GRAZE_ANGLE_OUT") = py::int_(Bmad::GRAZE_ANGLE_OUT);
  m.attr("R0_MAG") = py::int_(Bmad::R0_MAG);
  m.attr("RF_WAVELENGTH") = py::int_(Bmad::RF_WAVELENGTH);
  m.attr("SIG_VX") = py::int_(Bmad::SIG_VX);
  m.attr("SIG_VY") = py::int_(Bmad::SIG_VY);
  m.attr("CONSTANT_REF_ENERGY") = py::int_(Bmad::CONSTANT_REF_ENERGY);
  m.attr("KS") = py::int_(Bmad::KS);
  m.attr("SIG_E") = py::int_(Bmad::SIG_E);
  m.attr("SIG_PZ") = py::int_(Bmad::SIG_PZ);
  m.attr("AUTOSCALE_AMPLITUDE") = py::int_(Bmad::AUTOSCALE_AMPLITUDE);
  m.attr("D1_THICKNESS") = py::int_(Bmad::D1_THICKNESS);
  m.attr("DEFAULT_TRACKING_SPECIES") = py::int_(Bmad::DEFAULT_TRACKING_SPECIES);
  m.attr("AUTOSCALE_PHASE") = py::int_(Bmad::AUTOSCALE_PHASE);
  m.attr("N_SLICE") = py::int_(Bmad::N_SLICE);
  m.attr("Y_GAIN_CALIB") = py::int_(Bmad::Y_GAIN_CALIB);
  m.attr("SIG_E2") = py::int_(Bmad::SIG_E2);
  m.attr("FB1") = py::int_(Bmad::FB1);
  m.attr("POLARITY") = py::int_(Bmad::POLARITY);
  m.attr("CRUNCH_CALIB") = py::int_(Bmad::CRUNCH_CALIB);
  m.attr("ALPHA_ANGLE") = py::int_(Bmad::ALPHA_ANGLE);
  m.attr("D2_THICKNESS") = py::int_(Bmad::D2_THICKNESS);
  m.attr("BETA_A_STRONG") = py::int_(Bmad::BETA_A_STRONG);
  m.attr("BETA_A_OUT") = py::int_(Bmad::BETA_A_OUT);
  m.attr("E_LOSS") = py::int_(Bmad::E_LOSS);
  m.attr("GAP") = py::int_(Bmad::GAP);
  m.attr("SPIN_X") = py::int_(Bmad::SPIN_X);
  m.attr("E_CENTER") = py::int_(Bmad::E_CENTER);
  m.attr("SCATTER_TEST") = py::int_(Bmad::SCATTER_TEST);
  m.attr("FB2") = py::int_(Bmad::FB2);
  m.attr("X_OFFSET_CALIB") = py::int_(Bmad::X_OFFSET_CALIB);
  m.attr("V1_UNITCELL") = py::int_(Bmad::V1_UNITCELL);
  m.attr("PSI_ANGLE") = py::int_(Bmad::PSI_ANGLE);
  m.attr("CAVITY_TYPE") = py::int_(Bmad::CAVITY_TYPE);
  m.attr("BETA_B_STRONG") = py::int_(Bmad::BETA_B_STRONG);
  m.attr("BETA_B_OUT") = py::int_(Bmad::BETA_B_OUT);
  m.attr("SPIN_Y") = py::int_(Bmad::SPIN_Y);
  m.attr("E2_CENTER") = py::int_(Bmad::E2_CENTER);
  m.attr("N_PERIOD") = py::int_(Bmad::N_PERIOD);
  m.attr("EMIT_FRACTION") = py::int_(Bmad::EMIT_FRACTION);
  m.attr("X1_EDGE") = py::int_(Bmad::X1_EDGE);
  m.attr("Y_OFFSET_CALIB") = py::int_(Bmad::Y_OFFSET_CALIB);
  m.attr("V_UNITCELL") = py::int_(Bmad::V_UNITCELL);
  m.attr("V2_UNITCELL") = py::int_(Bmad::V2_UNITCELL);
  m.attr("SPIN_Z") = py::int_(Bmad::SPIN_Z);
  m.attr("L_PERIOD") = py::int_(Bmad::L_PERIOD);
  m.attr("FQ1") = py::int_(Bmad::FQ1);
  m.attr("ALPHA_A_STRONG") = py::int_(Bmad::ALPHA_A_STRONG);
  m.attr("ALPHA_A_OUT") = py::int_(Bmad::ALPHA_A_OUT);
  m.attr("E2_PROBABILITY") = py::int_(Bmad::E2_PROBABILITY);
  m.attr("PHI0_MAX") = py::int_(Bmad::PHI0_MAX);
  m.attr("X2_EDGE") = py::int_(Bmad::X2_EDGE);
  m.attr("FQ2") = py::int_(Bmad::FQ2);
  m.attr("PHI0") = py::int_(Bmad::PHI0);
  m.attr("TILT_CALIB") = py::int_(Bmad::TILT_CALIB);
  m.attr("E_CENTER_RELATIVE_TO_REF") = py::int_(Bmad::E_CENTER_RELATIVE_TO_REF);
  m.attr("Y1_EDGE") = py::int_(Bmad::Y1_EDGE);
  m.attr("ALPHA_B_STRONG") = py::int_(Bmad::ALPHA_B_STRONG);
  m.attr("ALPHA_B_OUT") = py::int_(Bmad::ALPHA_B_OUT);
  m.attr("IS_MOSAIC") = py::int_(Bmad::IS_MOSAIC);
  m.attr("PX_APERTURE_WIDTH2") = py::int_(Bmad::PX_APERTURE_WIDTH2);
  m.attr("PHI0_ERR") = py::int_(Bmad::PHI0_ERR);
  m.attr("CURRENT") = py::int_(Bmad::CURRENT);
  m.attr("MOSAIC_THICKNESS") = py::int_(Bmad::MOSAIC_THICKNESS);
  m.attr("PX_APERTURE_CENTER") = py::int_(Bmad::PX_APERTURE_CENTER);
  m.attr("ETA_X_OUT") = py::int_(Bmad::ETA_X_OUT);
  m.attr("QUAD_TILT") = py::int_(Bmad::QUAD_TILT);
  m.attr("DE_ETA_MEAS") = py::int_(Bmad::DE_ETA_MEAS);
  m.attr("SPATIAL_DISTRIBUTION") = py::int_(Bmad::SPATIAL_DISTRIBUTION);
  m.attr("Y2_EDGE") = py::int_(Bmad::Y2_EDGE);
  m.attr("SPECIES_STRONG") = py::int_(Bmad::SPECIES_STRONG);
  m.attr("ETA_Y_OUT") = py::int_(Bmad::ETA_Y_OUT);
  m.attr("MODE") = py::int_(Bmad::MODE);
  m.attr("VELOCITY_DISTRIBUTION") = py::int_(Bmad::VELOCITY_DISTRIBUTION);
  m.attr("PY_APERTURE_WIDTH2") = py::int_(Bmad::PY_APERTURE_WIDTH2);
  m.attr("PHI0_MULTIPASS") = py::int_(Bmad::PHI0_MULTIPASS);
  m.attr("N_SAMPLE") = py::int_(Bmad::N_SAMPLE);
  m.attr("ORIGIN_ELE_REF_PT") = py::int_(Bmad::ORIGIN_ELE_REF_PT);
  m.attr("MOSAIC_ANGLE_RMS_IN_PLANE") =
      py::int_(Bmad::MOSAIC_ANGLE_RMS_IN_PLANE);
  m.attr("EPS_STEP_SCALE") = py::int_(Bmad::EPS_STEP_SCALE);
  m.attr("E_TOT_STRONG") = py::int_(Bmad::E_TOT_STRONG);
  m.attr("DTHICKNESS_DX") = py::int_(Bmad::DTHICKNESS_DX);
  m.attr("BEND_TILT") = py::int_(Bmad::BEND_TILT);
  m.attr("ETAP_X_OUT") = py::int_(Bmad::ETAP_X_OUT);
  m.attr("PHI0_AUTOSCALE") = py::int_(Bmad::PHI0_AUTOSCALE);
  m.attr("DX_ORIGIN") = py::int_(Bmad::DX_ORIGIN);
  m.attr("ENERGY_DISTRIBUTION") = py::int_(Bmad::ENERGY_DISTRIBUTION);
  m.attr("X_QUAD") = py::int_(Bmad::X_QUAD);
  m.attr("DS_PHOTON_SLICE") = py::int_(Bmad::DS_PHOTON_SLICE);
  m.attr("MOSAIC_ANGLE_RMS_OUT_PLANE") =
      py::int_(Bmad::MOSAIC_ANGLE_RMS_OUT_PLANE);
  m.attr("PY_APERTURE_CENTER") = py::int_(Bmad::PY_APERTURE_CENTER);
  m.attr("X_DISPERSION_ERR") = py::int_(Bmad::X_DISPERSION_ERR);
  m.attr("L_RECTANGLE") = py::int_(Bmad::L_RECTANGLE);
  m.attr("PC_STRONG") = py::int_(Bmad::PC_STRONG);
  m.attr("ETAP_Y_OUT") = py::int_(Bmad::ETAP_Y_OUT);
  m.attr("DY_ORIGIN") = py::int_(Bmad::DY_ORIGIN);
  m.attr("Y_QUAD") = py::int_(Bmad::Y_QUAD);
  m.attr("E_FIELD_X") = py::int_(Bmad::E_FIELD_X);
  m.attr("Y_DISPERSION_ERR") = py::int_(Bmad::Y_DISPERSION_ERR);
  m.attr("Z_APERTURE_WIDTH2") = py::int_(Bmad::Z_APERTURE_WIDTH2);
  m.attr("USER_SETS_LENGTH") = py::int_(Bmad::USER_SETS_LENGTH);
  m.attr("B_FIELD_TOT") = py::int_(Bmad::B_FIELD_TOT);
  m.attr("UPSTREAM_COORD_DIR") = py::int_(Bmad::UPSTREAM_COORD_DIR);
  m.attr("DZ_ORIGIN") = py::int_(Bmad::DZ_ORIGIN);
  m.attr("MOSAIC_DIFFRACTION_NUM") = py::int_(Bmad::MOSAIC_DIFFRACTION_NUM);
  m.attr("CMAT_11") = py::int_(Bmad::CMAT_11);
  m.attr("FIELD_AUTOSCALE") = py::int_(Bmad::FIELD_AUTOSCALE);
  m.attr("L_SAGITTA") = py::int_(Bmad::L_SAGITTA);
  m.attr("E_FIELD_Y") = py::int_(Bmad::E_FIELD_Y);
  m.attr("X_DISPERSION_CALIB") = py::int_(Bmad::X_DISPERSION_CALIB);
  m.attr("Z_APERTURE_CENTER") = py::int_(Bmad::Z_APERTURE_CENTER);
  m.attr("F_FACTOR") = py::int_(Bmad::F_FACTOR);
  m.attr("CMAT_12") = py::int_(Bmad::CMAT_12);
  m.attr("DTHETA_ORIGIN") = py::int_(Bmad::DTHETA_ORIGIN);
  m.attr("B_PARAM") = py::int_(Bmad::B_PARAM);
  m.attr("L_CHORD") = py::int_(Bmad::L_CHORD);
  m.attr("DOWNSTREAM_COORD_DIR") = py::int_(Bmad::DOWNSTREAM_COORD_DIR);
  m.attr("PZ_APERTURE_WIDTH2") = py::int_(Bmad::PZ_APERTURE_WIDTH2);
  m.attr("Y_DISPERSION_CALIB") = py::int_(Bmad::Y_DISPERSION_CALIB);
  m.attr("SCALE_FIELD_TO_ONE") = py::int_(Bmad::SCALE_FIELD_TO_ONE);
  m.attr("VOLTAGE_TOT") = py::int_(Bmad::VOLTAGE_TOT);
  m.attr("SCATTER_METHOD") = py::int_(Bmad::SCATTER_METHOD);
  m.attr("CMAT_21") = py::int_(Bmad::CMAT_21);
  m.attr("L_ACTIVE") = py::int_(Bmad::L_ACTIVE);
  m.attr("DPHI_ORIGIN") = py::int_(Bmad::DPHI_ORIGIN);
  m.attr("SPLIT_ID") = py::int_(Bmad::SPLIT_ID);
  m.attr("REF_CAP_GAMMA") = py::int_(Bmad::REF_CAP_GAMMA);
  m.attr("L_SOFT_EDGE") = py::int_(Bmad::L_SOFT_EDGE);
  m.attr("TRANSVERSE_SIGMA_CUT") = py::int_(Bmad::TRANSVERSE_SIGMA_CUT);
  m.attr("PZ_APERTURE_CENTER") = py::int_(Bmad::PZ_APERTURE_CENTER);
  m.attr("MEAN_EXCITATION_ENERGY") = py::int_(Bmad::MEAN_EXCITATION_ENERGY);
  m.attr("FIDUCIAL_PT") = py::int_(Bmad::FIDUCIAL_PT);
  m.attr("CMAT_22") = py::int_(Bmad::CMAT_22);
  m.attr("DPSI_ORIGIN") = py::int_(Bmad::DPSI_ORIGIN);
  m.attr("T_OFFSET") = py::int_(Bmad::T_OFFSET);
  m.attr("DS_SLICE") = py::int_(Bmad::DS_SLICE);
  m.attr("USE_REFLECTIVITY_TABLE") = py::int_(Bmad::USE_REFLECTIVITY_TABLE);
  m.attr("INIT_NEEDED") = py::int_(Bmad::INIT_NEEDED);
  m.attr("LONGITUDINAL_MODE") = py::int_(Bmad::LONGITUDINAL_MODE);
  m.attr("ANGLE") = py::int_(Bmad::ANGLE);
  m.attr("N_CELL") = py::int_(Bmad::N_CELL);
  m.attr("MODE_FLIP") = py::int_(Bmad::MODE_FLIP);
  m.attr("CROSSING_TIME") = py::int_(Bmad::CROSSING_TIME);
  m.attr("X_KICK") = py::int_(Bmad::X_KICK);
  // Note: [x_kick$, px_kick$, ..., pz_kick$] must be in order.
  m.attr("X_PITCH") = py::int_(Bmad::X_PITCH);
  m.attr("PX_KICK") = py::int_(Bmad::PX_KICK);
  m.attr("Y_PITCH") = py::int_(Bmad::Y_PITCH);
  m.attr("Y_KICK") = py::int_(Bmad::Y_KICK);
  m.attr("X_OFFSET") = py::int_(Bmad::X_OFFSET);
  m.attr("PY_KICK") = py::int_(Bmad::PY_KICK);
  m.attr("Y_OFFSET") = py::int_(Bmad::Y_OFFSET);
  m.attr("Z_KICK") = py::int_(Bmad::Z_KICK);
  m.attr("Z_OFFSET") = py::int_(Bmad::Z_OFFSET);
  m.attr("PZ_KICK") = py::int_(Bmad::PZ_KICK);
  m.attr("HKICK") = py::int_(Bmad::HKICK);
  m.attr("D_SPACING") = py::int_(Bmad::D_SPACING);
  m.attr("X_OFFSET_MULT") = py::int_(Bmad::X_OFFSET_MULT);
  m.attr("EMITTANCE_A") = py::int_(Bmad::EMITTANCE_A);
  m.attr("CRAB_X1") = py::int_(Bmad::CRAB_X1);
  m.attr("VKICK") = py::int_(Bmad::VKICK);
  m.attr("Y_OFFSET_MULT") = py::int_(Bmad::Y_OFFSET_MULT);
  m.attr("P0C_REF_INIT") = py::int_(Bmad::P0C_REF_INIT);
  m.attr("EMITTANCE_B") = py::int_(Bmad::EMITTANCE_B);
  m.attr("CRAB_X2") = py::int_(Bmad::CRAB_X2);
  m.attr("BL_HKICK") = py::int_(Bmad::BL_HKICK);
  m.attr("E_TOT_REF_INIT") = py::int_(Bmad::E_TOT_REF_INIT);
  m.attr("EMITTANCE_Z") = py::int_(Bmad::EMITTANCE_Z);
  m.attr("CRAB_X3") = py::int_(Bmad::CRAB_X3);
  m.attr("BL_VKICK") = py::int_(Bmad::BL_VKICK);
  m.attr("CRAB_TILT") = py::int_(Bmad::CRAB_TILT);
  m.attr("BL_KICK") = py::int_(Bmad::BL_KICK);
  m.attr("B_FIELD") = py::int_(Bmad::B_FIELD);
  m.attr("E_FIELD") = py::int_(Bmad::E_FIELD);
  m.attr("HIGH_ENERGY_SPACE_CHARGE_ON") =
      py::int_(Bmad::HIGH_ENERGY_SPACE_CHARGE_ON);
  m.attr("CRAB_X4") = py::int_(Bmad::CRAB_X4);
  m.attr("N_RF_STEPS") = py::int_(Bmad::N_RF_STEPS);
  m.attr("PHOTON_TYPE") = py::int_(Bmad::PHOTON_TYPE);
  m.attr("COUPLER_PHASE") = py::int_(Bmad::COUPLER_PHASE);
  m.attr("DB_FIELD") = py::int_(Bmad::DB_FIELD);
  m.attr("CRAB_X5") = py::int_(Bmad::CRAB_X5);
  m.attr("LATTICE_TYPE") = py::int_(Bmad::LATTICE_TYPE);
  m.attr("B1_GRADIENT") = py::int_(Bmad::B1_GRADIENT);
  m.attr("E1_GRADIENT") = py::int_(Bmad::E1_GRADIENT);
  m.attr("COUPLER_ANGLE") = py::int_(Bmad::COUPLER_ANGLE);
  m.attr("LIVE_BRANCH") = py::int_(Bmad::LIVE_BRANCH);
  m.attr("B2_GRADIENT") = py::int_(Bmad::B2_GRADIENT);
  m.attr("E2_GRADIENT") = py::int_(Bmad::E2_GRADIENT);
  m.attr("COUPLER_STRENGTH") = py::int_(Bmad::COUPLER_STRENGTH);
  m.attr("GEOMETRY") = py::int_(Bmad::GEOMETRY);
  m.attr("COUPLER_AT") = py::int_(Bmad::COUPLER_AT);
  m.attr("E_TOT_OFFSET") = py::int_(Bmad::E_TOT_OFFSET);
  m.attr("PTC_CANONICAL_COORDS") = py::int_(Bmad::PTC_CANONICAL_COORDS);
  m.attr("B3_GRADIENT") = py::int_(Bmad::B3_GRADIENT);
  m.attr("E3_GRADIENT") = py::int_(Bmad::E3_GRADIENT);
  m.attr("PTC_FRINGE_GEOMETRY") = py::int_(Bmad::PTC_FRINGE_GEOMETRY);
  m.attr("E_TOT_SET") = py::int_(Bmad::E_TOT_SET);
  m.attr("BS_FIELD") = py::int_(Bmad::BS_FIELD);
  m.attr("P0C_SET") = py::int_(Bmad::P0C_SET);
  m.attr("PTC_FIELD_GEOMETRY") = py::int_(Bmad::PTC_FIELD_GEOMETRY);
  m.attr("DELTA_REF_TIME_USER_SET") = py::int_(Bmad::DELTA_REF_TIME_USER_SET);
  m.attr("DELTA_REF_TIME") = py::int_(Bmad::DELTA_REF_TIME);
  m.attr("P0C_START") = py::int_(Bmad::P0C_START);
  m.attr("E_TOT_START") = py::int_(Bmad::E_TOT_START);
  m.attr("P0C") = py::int_(Bmad::P0C);
  m.attr("E_TOT") = py::int_(Bmad::E_TOT);
  m.attr("X_PITCH_TOT") = py::int_(Bmad::X_PITCH_TOT);
  m.attr("NO_END_MARKER") = py::int_(Bmad::NO_END_MARKER);
  m.attr("Y_PITCH_TOT") = py::int_(Bmad::Y_PITCH_TOT);
  m.attr("X_OFFSET_TOT") = py::int_(Bmad::X_OFFSET_TOT);
  m.attr("Y_OFFSET_TOT") = py::int_(Bmad::Y_OFFSET_TOT);
  m.attr("Z_OFFSET_TOT") = py::int_(Bmad::Z_OFFSET_TOT);
  // Important: tilt_tot$ = roll_tot$
  m.attr("TILT_TOT") = py::int_(Bmad::TILT_TOT);
  m.attr("ROLL_TOT") = py::int_(Bmad::ROLL_TOT);
  m.attr("REF_TILT_TOT") = py::int_(Bmad::REF_TILT_TOT);
  m.attr("MULTIPASS_REF_ENERGY") = py::int_(Bmad::MULTIPASS_REF_ENERGY);
  m.attr("DISPATCH") = py::int_(Bmad::DISPATCH);
  m.attr("REF_TIME_START") = py::int_(Bmad::REF_TIME_START);
  // For Etiennes' PTC: 2, 4, 6, or 8.
  m.attr("THICKNESS") = py::int_(Bmad::THICKNESS);
  m.attr("INTEGRATOR_ORDER") = py::int_(Bmad::INTEGRATOR_ORDER);
  // Assumed unique by set_flags_for_changed_real_attribute
  m.attr("NUM_STEPS") = py::int_(Bmad::NUM_STEPS);
  // Assumed unique by set_flags_for_changed_real_attribute
  m.attr("DS_STEP") = py::int_(Bmad::DS_STEP);
  m.attr("CSR_DS_STEP") = py::int_(Bmad::CSR_DS_STEP);
  m.attr("LORD_PAD1") = py::int_(Bmad::LORD_PAD1);
  m.attr("LORD_PAD2") = py::int_(Bmad::LORD_PAD2);
  m.attr("REF_WAVELENGTH") = py::int_(Bmad::REF_WAVELENGTH);
  m.attr("X1_LIMIT") = py::int_(Bmad::X1_LIMIT);
  m.attr("X2_LIMIT") = py::int_(Bmad::X2_LIMIT);
  m.attr("Y1_LIMIT") = py::int_(Bmad::Y1_LIMIT);
  m.attr("Y2_LIMIT") = py::int_(Bmad::Y2_LIMIT);
  m.attr("CHECK_SUM") = py::int_(Bmad::CHECK_SUM);
  m.attr("IS_ON") = py::int_(Bmad::IS_ON);
  m.attr("ALIAS") = py::int_(Bmad::ALIAS);
  m.attr("DISTRIBUTION") = py::int_(Bmad::DISTRIBUTION);
  m.attr("TT") = py::int_(Bmad::TT);
  m.attr("X_KNOT") = py::int_(Bmad::X_KNOT);
  m.attr("MAX_FRINGE_ORDER") = py::int_(Bmad::MAX_FRINGE_ORDER);
  m.attr("ETA_X") = py::int_(Bmad::ETA_X);
  m.attr("ELECTRIC_DIPOLE_MOMENT") = py::int_(Bmad::ELECTRIC_DIPOLE_MOMENT);
  m.attr("LR_SELF_WAKE_ON") = py::int_(Bmad::LR_SELF_WAKE_ON);
  m.attr("X_REF") = py::int_(Bmad::X_REF);
  m.attr("SPECIES_OUT") = py::int_(Bmad::SPECIES_OUT);
  m.attr("Y_KNOT") = py::int_(Bmad::Y_KNOT);
  m.attr("ETA_Y") = py::int_(Bmad::ETA_Y);
  m.attr("DENSITY") = py::int_(Bmad::DENSITY);
  m.attr("LR_WAKE_FILE") = py::int_(Bmad::LR_WAKE_FILE);
  m.attr("PX_REF") = py::int_(Bmad::PX_REF);
  m.attr("ETAP_X") = py::int_(Bmad::ETAP_X);
  m.attr("SLAVE") = py::int_(Bmad::SLAVE);
  m.attr("DENSITY_USED") = py::int_(Bmad::DENSITY_USED);
  m.attr("LR_FREQ_SPREAD") = py::int_(Bmad::LR_FREQ_SPREAD);
  m.attr("Y_REF") = py::int_(Bmad::Y_REF);
  m.attr("ETAP_Y") = py::int_(Bmad::ETAP_Y);
  m.attr("AREA_DENSITY") = py::int_(Bmad::AREA_DENSITY);
  m.attr("INPUT_ELE") = py::int_(Bmad::INPUT_ELE);
  m.attr("LATTICE") = py::int_(Bmad::LATTICE);
  m.attr("PHI_A") = py::int_(Bmad::PHI_A);
  m.attr("MULTIPOLES_ON") = py::int_(Bmad::MULTIPOLES_ON);
  m.attr("PY_REF") = py::int_(Bmad::PY_REF);
  m.attr("AREA_DENSITY_USED") = py::int_(Bmad::AREA_DENSITY_USED);
  m.attr("OUTPUT_ELE") = py::int_(Bmad::OUTPUT_ELE);
  m.attr("APERTURE_TYPE") = py::int_(Bmad::APERTURE_TYPE);
  m.attr("ETA_Z") = py::int_(Bmad::ETA_Z);
  m.attr("MACHINE") = py::int_(Bmad::MACHINE);
  m.attr("TAYLOR_MAP_INCLUDES_OFFSETS") =
      py::int_(Bmad::TAYLOR_MAP_INCLUDES_OFFSETS);
  m.attr("PIXEL") = py::int_(Bmad::PIXEL);
  m.attr("P88") = py::int_(Bmad::P88);
  m.attr("RADIATION_LENGTH") = py::int_(Bmad::RADIATION_LENGTH);
  m.attr("DETA_DPZ_X") = py::int_(Bmad::DETA_DPZ_X);
  m.attr("CSR_METHOD") = py::int_(Bmad::CSR_METHOD);
  m.attr("VAR") = py::int_(Bmad::VAR);
  m.attr("Z_REF") = py::int_(Bmad::Z_REF);
  m.attr("P89") = py::int_(Bmad::P89);
  m.attr("RADIATION_LENGTH_USED") = py::int_(Bmad::RADIATION_LENGTH_USED);
  m.attr("PZ_REF") = py::int_(Bmad::PZ_REF);
  m.attr("SPACE_CHARGE_METHOD") = py::int_(Bmad::SPACE_CHARGE_METHOD);
  m.attr("P90") = py::int_(Bmad::P90);
  m.attr("DETAP_DPZ_X") = py::int_(Bmad::DETAP_DPZ_X);
  m.attr("MAT6_CALC_METHOD") = py::int_(Bmad::MAT6_CALC_METHOD);
  m.attr("TRACKING_METHOD") = py::int_(Bmad::TRACKING_METHOD);
  m.attr("REF_TIME") = py::int_(Bmad::REF_TIME);
  m.attr("PTC_INTEGRATION_TYPE") = py::int_(Bmad::PTC_INTEGRATION_TYPE);
  m.attr("SPIN_TRACKING_METHOD") = py::int_(Bmad::SPIN_TRACKING_METHOD);
  m.attr("ETA_A") = py::int_(Bmad::ETA_A);
  m.attr("APERTURE") = py::int_(Bmad::APERTURE);
  m.attr("ETAP_A") = py::int_(Bmad::ETAP_A);
  m.attr("DETA_DPZ_Y") = py::int_(Bmad::DETA_DPZ_Y);
  m.attr("X_LIMIT") = py::int_(Bmad::X_LIMIT);
  m.attr("ABSOLUTE_TIME_TRACKING") = py::int_(Bmad::ABSOLUTE_TIME_TRACKING);
  m.attr("ETA_B") = py::int_(Bmad::ETA_B);
  m.attr("DETAP_DPZ_Y") = py::int_(Bmad::DETAP_DPZ_Y);
  m.attr("Y_LIMIT") = py::int_(Bmad::Y_LIMIT);
  m.attr("ETAP_B") = py::int_(Bmad::ETAP_B);
  m.attr("OFFSET_MOVES_APERTURE") = py::int_(Bmad::OFFSET_MOVES_APERTURE);
  m.attr("ALPHA_A") = py::int_(Bmad::ALPHA_A);
  m.attr("REFLECTIVITY_TABLE") = py::int_(Bmad::REFLECTIVITY_TABLE);
  m.attr("ENERGY_PROBABILITY_CURVE") = py::int_(Bmad::ENERGY_PROBABILITY_CURVE);
  m.attr("EXACT_MISALIGN") = py::int_(Bmad::EXACT_MISALIGN);
  m.attr("PHYSICAL_SOURCE") = py::int_(Bmad::PHYSICAL_SOURCE);
  m.attr("SR_WAKE_FILE") = py::int_(Bmad::SR_WAKE_FILE);
  m.attr("ALPHA_B") = py::int_(Bmad::ALPHA_B);
  m.attr("TERM") = py::int_(Bmad::TERM);
  m.attr("FREQUENCIES") = py::int_(Bmad::FREQUENCIES);
  m.attr("OLD_INTEGRATOR") = py::int_(Bmad::OLD_INTEGRATOR);
  m.attr("CURVATURE") = py::int_(Bmad::CURVATURE);
  m.attr("S_LONG") = py::int_(Bmad::S_LONG);
  m.attr("X_POSITION") = py::int_(Bmad::X_POSITION);
  m.attr("EXACT_MODEL") = py::int_(Bmad::EXACT_MODEL);
  m.attr("SYMPLECTIFY") = py::int_(Bmad::SYMPLECTIFY);
  m.attr("Y_POSITION") = py::int_(Bmad::Y_POSITION);
  m.attr("N_SLICE_SPLINE") = py::int_(Bmad::N_SLICE_SPLINE);
  m.attr("Z_POSITION") = py::int_(Bmad::Z_POSITION);
  m.attr("AMP_VS_TIME") = py::int_(Bmad::AMP_VS_TIME);
  m.attr("THETA_POSITION") = py::int_(Bmad::THETA_POSITION);
  m.attr("VERTICAL_KICK") = py::int_(Bmad::VERTICAL_KICK);
  m.attr("FIELD_CALC") = py::int_(Bmad::FIELD_CALC);
  m.attr("PHI_POSITION") = py::int_(Bmad::PHI_POSITION);
  m.attr("PSI_POSITION") = py::int_(Bmad::PSI_POSITION);
  m.attr("WALL") = py::int_(Bmad::WALL);
  m.attr("APERTURE_AT") = py::int_(Bmad::APERTURE_AT);
  m.attr("BETA_A") = py::int_(Bmad::BETA_A);
  m.attr("RAN_SEED") = py::int_(Bmad::RAN_SEED);
  m.attr("ORIGIN_ELE") = py::int_(Bmad::ORIGIN_ELE);
  m.attr("BETA_B") = py::int_(Bmad::BETA_B);
  m.attr("TO_LINE") = py::int_(Bmad::TO_LINE);
  m.attr("FIELD_OVERLAPS") = py::int_(Bmad::FIELD_OVERLAPS);
  m.attr("DBETA_DPZ_A") = py::int_(Bmad::DBETA_DPZ_A);
  m.attr("FIELD_MASTER") = py::int_(Bmad::FIELD_MASTER);
  m.attr("TO_ELEMENT") = py::int_(Bmad::TO_ELEMENT);
  m.attr("DBETA_DPZ_B") = py::int_(Bmad::DBETA_DPZ_B);
  m.attr("DESCRIP") = py::int_(Bmad::DESCRIP);
  m.attr("SCALE_MULTIPOLES") = py::int_(Bmad::SCALE_MULTIPOLES);
  m.attr("DALPHA_DPZ_A") = py::int_(Bmad::DALPHA_DPZ_A);
  m.attr("SR_WAKE") = py::int_(Bmad::SR_WAKE);
  m.attr("DALPHA_DPZ_B") = py::int_(Bmad::DALPHA_DPZ_B);
  m.attr("REF_ORBIT") = py::int_(Bmad::REF_ORBIT);
  m.attr("LR_WAKE") = py::int_(Bmad::LR_WAKE);
  m.attr("PHI_B") = py::int_(Bmad::PHI_B);
  m.attr("CRYSTAL_TYPE") = py::int_(Bmad::CRYSTAL_TYPE);
  m.attr("MATERIAL_TYPE") = py::int_(Bmad::MATERIAL_TYPE);
  m.attr("TYPE") = py::int_(Bmad::TYPE);
  m.attr("REF_ORIGIN") = py::int_(Bmad::REF_ORIGIN);
  m.attr("ELE_ORIGIN") = py::int_(Bmad::ELE_ORIGIN);
  m.attr("SUPERIMPOSE") = py::int_(Bmad::SUPERIMPOSE);
  m.attr("SUPER_OFFSET") = py::int_(Bmad::SUPER_OFFSET);
  m.attr("REFERENCE") = py::int_(Bmad::REFERENCE);
  m.attr("CARTESIAN_MAP") = py::int_(Bmad::CARTESIAN_MAP);
  m.attr("CYLINDRICAL_MAP") = py::int_(Bmad::CYLINDRICAL_MAP);
  m.attr("GRID_FIELD") = py::int_(Bmad::GRID_FIELD);
  m.attr("GEN_GRAD_MAP") = py::int_(Bmad::GEN_GRAD_MAP);
  m.attr("CREATE_JUMBO_SLAVE") = py::int_(Bmad::CREATE_JUMBO_SLAVE);
  m.attr("ACCORDION_EDGE") = py::int_(Bmad::ACCORDION_EDGE);
  m.attr("START_EDGE") = py::int_(Bmad::START_EDGE);
  m.attr("END_EDGE") = py::int_(Bmad::END_EDGE);
  m.attr("S_POSITION") = py::int_(Bmad::S_POSITION);
  m.attr("REF_SPECIES") = py::int_(Bmad::REF_SPECIES);
  m.attr("PARTICLE") = py::int_(Bmad::PARTICLE);
  m.attr("WRAP_SUPERIMPOSE") = py::int_(Bmad::WRAP_SUPERIMPOSE);
  m.attr("A0") = py::int_(Bmad::A0);
  m.attr("A21") = py::int_(Bmad::A21);
  m.attr("B0") = py::int_(Bmad::B0);
  m.attr("B21") = py::int_(Bmad::B21);
  m.attr("K0L") = py::int_(Bmad::K0L);
  m.attr("K21L") = py::int_(Bmad::K21L);
  m.attr("T0") = py::int_(Bmad::T0);
  m.attr("T21") = py::int_(Bmad::T21);
  m.attr("K0SL") = py::int_(Bmad::K0SL);
  m.attr("K21SL") = py::int_(Bmad::K21SL);
  m.attr("A0_ELEC") = py::int_(Bmad::A0_ELEC);
  m.attr("A21_ELEC") = py::int_(Bmad::A21_ELEC);
  m.attr("B0_ELEC") = py::int_(Bmad::B0_ELEC);
  m.attr("B21_ELEC") = py::int_(Bmad::B21_ELEC);
  m.attr("CUSTOM_ATTRIBUTE0") = py::int_(Bmad::CUSTOM_ATTRIBUTE0);
  m.attr("CUSTOM_ATTRIBUTE_NUM") = py::int_(Bmad::CUSTOM_ATTRIBUTE_NUM);
  m.attr("NUM_ELE_ATTRIB_EXTENDED") = py::int_(Bmad::NUM_ELE_ATTRIB_EXTENDED);
  // For backwards compatibility.
  m.attr("G_ERR") = py::int_(Bmad::G_ERR);
  // For backwards compatibility
  m.attr("B_FIELD_ERR") = py::int_(Bmad::B_FIELD_ERR);
  m.attr("OPEN") = py::int_(Bmad::OPEN);
  m.attr("CLOSED") = py::int_(Bmad::CLOSED);
  m.attr("BENDS") = py::int_(Bmad::BENDS);
  m.attr("WIGGLERS") = py::int_(Bmad::WIGGLERS);
  m.attr("ALL") = py::int_(Bmad::ALL);
  m.attr("UPSTREAM") = py::int_(Bmad::UPSTREAM);
  m.attr("DOWNSTREAM") = py::int_(Bmad::DOWNSTREAM);
  m.attr("RADIANS") = py::int_(Bmad::RADIANS);
  m.attr("DEGREES") = py::int_(Bmad::DEGREES);
  m.attr("CYCLES") = py::int_(Bmad::CYCLES);
  m.attr("RADIANS_OVER_2PI") = py::int_(Bmad::RADIANS_OVER_2PI);
  m.attr("ROTATIONALLY_SYMMETRIC_RZ") =
      py::int_(Bmad::ROTATIONALLY_SYMMETRIC_RZ);
  m.attr("XYZ") = py::int_(Bmad::XYZ);
  m.attr("INVALID_NAME") = py::int_(Bmad::INVALID_NAME);
  m.attr("IS_LOGICAL") = py::int_(Bmad::IS_LOGICAL);
  m.attr("IS_INTEGER") = py::int_(Bmad::IS_INTEGER);
  m.attr("IS_REAL") = py::int_(Bmad::IS_REAL);
  m.attr("IS_SWITCH") = py::int_(Bmad::IS_SWITCH);
  m.attr("IS_STRING") = py::int_(Bmad::IS_STRING);
  m.attr("IS_STRUCT") = py::int_(Bmad::IS_STRUCT);
  m.attr("UNKNOWN") = py::int_(Bmad::UNKNOWN);
  m.attr("PATCH_PROBLEM") = py::int_(Bmad::PATCH_PROBLEM);
  m.attr("CANNOT_FIND") = py::int_(Bmad::CANNOT_FIND);
  m.attr("OUTSIDE") = py::int_(Bmad::OUTSIDE);
  m.attr("SMALL_REL_CHANGE") = py::float_(Bmad::SMALL_REL_CHANGE);
  m.attr("END_STACK") = py::int_(Bmad::END_STACK);
  m.attr("PLUS") = py::int_(Bmad::PLUS);
  m.attr("MINUS") = py::int_(Bmad::MINUS);
  m.attr("TIMES") = py::int_(Bmad::TIMES);
  m.attr("DIVIDE") = py::int_(Bmad::DIVIDE);
  m.attr("L_PARENS") = py::int_(Bmad::L_PARENS);
  m.attr("R_PARENS") = py::int_(Bmad::R_PARENS);
  m.attr("POWER") = py::int_(Bmad::POWER);
  m.attr("UNARY_MINUS") = py::int_(Bmad::UNARY_MINUS);
  m.attr("UNARY_PLUS") = py::int_(Bmad::UNARY_PLUS);
  m.attr("NO_DELIM") = py::int_(Bmad::NO_DELIM);
  m.attr("SIN") = py::int_(Bmad::SIN);
  m.attr("COS") = py::int_(Bmad::COS);
  m.attr("TAN") = py::int_(Bmad::TAN);
  m.attr("ASIN") = py::int_(Bmad::ASIN);
  m.attr("ACOS") = py::int_(Bmad::ACOS);
  m.attr("ATAN") = py::int_(Bmad::ATAN);
  m.attr("ABS") = py::int_(Bmad::ABS);
  m.attr("SQRT") = py::int_(Bmad::SQRT);
  m.attr("LOG") = py::int_(Bmad::LOG);
  m.attr("EXP") = py::int_(Bmad::EXP);
  m.attr("RAN") = py::int_(Bmad::RAN);
  m.attr("RAN_GAUSS") = py::int_(Bmad::RAN_GAUSS);
  m.attr("ATAN2") = py::int_(Bmad::ATAN2);
  m.attr("FACTORIAL") = py::int_(Bmad::FACTORIAL);
  m.attr("INT") = py::int_(Bmad::INT);
  m.attr("NINT") = py::int_(Bmad::NINT);
  m.attr("FLOOR") = py::int_(Bmad::FLOOR);
  m.attr("CEILING") = py::int_(Bmad::CEILING);
  m.attr("NUMERIC") = py::int_(Bmad::NUMERIC);
  m.attr("VARIABLE") = py::int_(Bmad::VARIABLE);
  m.attr("MASS_OF") = py::int_(Bmad::MASS_OF);
  m.attr("CHARGE_OF") = py::int_(Bmad::CHARGE_OF);
  m.attr("ANOMALOUS_MOMENT_OF") = py::int_(Bmad::ANOMALOUS_MOMENT_OF);
  m.attr("SPECIES") = py::int_(Bmad::SPECIES);
  m.attr("SPECIES_CONST") = py::int_(Bmad::SPECIES_CONST);
  m.attr("SINC") = py::int_(Bmad::SINC);
  m.attr("CONSTANT") = py::int_(Bmad::CONSTANT);
  m.attr("COMMA") = py::int_(Bmad::COMMA);
  m.attr("RMS") = py::int_(Bmad::RMS);
  m.attr("AVERAGE") = py::int_(Bmad::AVERAGE);
  m.attr("SUM") = py::int_(Bmad::SUM);
  m.attr("ARG_COUNT") = py::int_(Bmad::ARG_COUNT);
  m.attr("ANTIPARTICLE") = py::int_(Bmad::ANTIPARTICLE);
  m.attr("COT") = py::int_(Bmad::COT);
  m.attr("SEC") = py::int_(Bmad::SEC);
  m.attr("CSC") = py::int_(Bmad::CSC);
  m.attr("SIGN") = py::int_(Bmad::SIGN);
  m.attr("L_FUNC_PARENS") = py::int_(Bmad::L_FUNC_PARENS);
  m.attr("SINH") = py::int_(Bmad::SINH);
  m.attr("COSH") = py::int_(Bmad::COSH);
  m.attr("TANH") = py::int_(Bmad::TANH);
  m.attr("COTH") = py::int_(Bmad::COTH);
  m.attr("ASINH") = py::int_(Bmad::ASINH);
  m.attr("ACOSH") = py::int_(Bmad::ACOSH);
  m.attr("ATANH") = py::int_(Bmad::ATANH);
  m.attr("ACOTH") = py::int_(Bmad::ACOTH);
  m.attr("MIN") = py::int_(Bmad::MIN);
  m.attr("MAX") = py::int_(Bmad::MAX);
  m.attr("MODULO") = py::int_(Bmad::MODULO);
  m.attr("ROOT") = py::int_(Bmad::ROOT);
  m.attr("PARENS") = py::int_(Bmad::PARENS);
  m.attr("SQUARE_BRACKETS") = py::int_(Bmad::SQUARE_BRACKETS);
  m.attr("CURLY_BRACKETS") = py::int_(Bmad::CURLY_BRACKETS);
  m.attr("FUNC_PARENS") = py::int_(Bmad::FUNC_PARENS);
  m.attr("ARROW") = py::int_(Bmad::ARROW);
  m.attr("EQUAL") = py::int_(Bmad::EQUAL);
  m.attr("COLON") = py::int_(Bmad::COLON);
  m.attr("DOUBLE_COLON") = py::int_(Bmad::DOUBLE_COLON);
  m.attr("COMPOUND") = py::int_(Bmad::COMPOUND);
  m.attr("FUNCTION") = py::int_(Bmad::FUNCTION);
  m.attr("VERTICAL_BAR") = py::int_(Bmad::VERTICAL_BAR);
  m.attr("BLANK") = py::int_(Bmad::BLANK);
  m.attr("AMPERSAND") = py::int_(Bmad::AMPERSAND);

  // Enums from output_mod.f90
  // No message printed. Used to override a status level variable.
  m.attr("S_NOOUTPUT") = py::int_(Bmad::S_NOOUTPUT);
  // Information message. The routine name is not printed.
  m.attr("S_BLANK") = py::int_(Bmad::S_BLANK);
  // Informational message.
  m.attr("S_INFO") = py::int_(Bmad::S_INFO);
  // Info message (w/timestamp).
  m.attr("S_DINFO") = py::int_(Bmad::S_DINFO);
  // Successful completion.
  m.attr("S_SUCCESS") = py::int_(Bmad::S_SUCCESS);
  // Warning of a possible problem.
  m.attr("S_WARN") = py::int_(Bmad::S_WARN);
  // Warning of a possible problem (w/timestamp).
  m.attr("S_DWARN") = py::int_(Bmad::S_DWARN);
  // An error as occurred [EG: bad user input] (w/ timestamp).
  m.attr("S_ERROR") = py::int_(Bmad::S_ERROR);
  // A fatal error has occurred so that computations
  m.attr("S_FATAL") = py::int_(Bmad::S_FATAL);
  // A severe error has occurred and
  m.attr("S_ABORT") = py::int_(Bmad::S_ABORT);
  // An important message.
  m.attr("S_IMPORTANT") = py::int_(Bmad::S_IMPORTANT);

  // Enums from physical_constants.f90
  m.attr("PI") = py::float_(Bmad::PI);
  m.attr("TWOPI") = py::float_(Bmad::TWOPI);
  m.attr("FOURPI") = py::float_(Bmad::FOURPI);
  m.attr("SQRT_2") = py::float_(Bmad::SQRT_2);
  m.attr("SQRT_3") = py::float_(Bmad::SQRT_3);
  // Mass [eV]
  m.attr("M_ELECTRON") = py::float_(Bmad::M_ELECTRON);
  // Mass [eV]
  m.attr("M_PROTON") = py::float_(Bmad::M_PROTON);
  // Mass [eV]
  m.attr("M_NEUTRON") = py::float_(Bmad::M_NEUTRON);
  // Mass [eV]
  m.attr("M_MUON") = py::float_(Bmad::M_MUON);
  // Mass He3 nucleus
  m.attr("M_HELION") = py::float_(Bmad::M_HELION);
  // [GeV] FOR MAD COMPATIBILITY USE ONLY. USE M_ELECTRON INSTEAD.
  m.attr("E_MASS") = py::float_(Bmad::E_MASS);
  // [GeV] FOR MAD COMPATIBILITY USE ONLY. USE M_PROTON INSTEAD.
  m.attr("P_MASS") = py::float_(Bmad::P_MASS);
  // Mass [eV]
  m.attr("M_PION_0") = py::float_(Bmad::M_PION_0);
  // Mass [eV]
  m.attr("M_PION_CHARGED") = py::float_(Bmad::M_PION_CHARGED);
  // Mass [eV]
  m.attr("M_DEUTERON") = py::float_(Bmad::M_DEUTERON);
  // unified atomic mass unit u (or dalton) in [eV]
  m.attr("ATOMIC_MASS_UNIT") = py::float_(Bmad::ATOMIC_MASS_UNIT);
  // speed of light
  m.attr("C_LIGHT") = py::float_(Bmad::C_LIGHT);
  // classical electron radius
  m.attr("R_E") = py::float_(Bmad::R_E);
  // proton radius
  m.attr("R_P") = py::float_(Bmad::R_P);
  // electron charge [Coul]
  m.attr("E_CHARGE") = py::float_(Bmad::E_CHARGE);
  // Planck's constant [eV*sec]
  m.attr("H_PLANCK") = py::float_(Bmad::H_PLANCK);
  // h_planck/twopi [eV*sec]
  m.attr("H_BAR_PLANCK") = py::float_(Bmad::H_BAR_PLANCK);
  // Vacuum permeability 2018 CODATA.
  m.attr("MU_0_VAC") = py::float_(Bmad::MU_0_VAC);
  // e^2 / (4 pi eps_0) [m*eV]
  m.attr("CLASSICAL_RADIUS_FACTOR") = py::float_(Bmad::CLASSICAL_RADIUS_FACTOR);
  // Number / mole  (exact)
  m.attr("N_AVOGADRO") = py::float_(Bmad::N_AVOGADRO);
  m.attr("FINE_STRUCTURE_CONSTANT") = py::float_(Bmad::FINE_STRUCTURE_CONSTANT);
  m.attr("ANOMALOUS_MAG_MOMENT_ELECTRON") =
      py::float_(Bmad::ANOMALOUS_MAG_MOMENT_ELECTRON);
  m.attr("ANOMALOUS_MAG_MOMENT_PROTON") =
      py::float_(Bmad::ANOMALOUS_MAG_MOMENT_PROTON);
  // ~fine_structure_constant / twopi
  m.attr("ANOMALOUS_MAG_MOMENT_MUON") =
      py::float_(Bmad::ANOMALOUS_MAG_MOMENT_MUON);
  m.attr("ANOMALOUS_MAG_MOMENT_DEUTERON") =
      py::float_(Bmad::ANOMALOUS_MAG_MOMENT_DEUTERON);
  m.attr("ANOMALOUS_MAG_MOMENT_NEUTRON") =
      py::float_(Bmad::ANOMALOUS_MAG_MOMENT_NEUTRON);
  m.attr("ANOMALOUS_MAG_MOMENT_HE3") =
      py::float_(Bmad::ANOMALOUS_MAG_MOMENT_HE3);

  // Enums from particle_species_mod.f90
  m.attr("PION_0") = py::int_(Bmad::PION_0);
  m.attr("HELION") = py::int_(Bmad::HELION);
  m.attr("REF_PARTICLE") = py::int_(Bmad::REF_PARTICLE);
  m.attr("NEUTRON") = py::int_(Bmad::NEUTRON);
  m.attr("DEUTERON") = py::int_(Bmad::DEUTERON);
  m.attr("PION_PLUS") = py::int_(Bmad::PION_PLUS);
  m.attr("ANTIMUON") = py::int_(Bmad::ANTIMUON);
  m.attr("PROTON") = py::int_(Bmad::PROTON);
  m.attr("POSITRON") = py::int_(Bmad::POSITRON);
  m.attr("PHOTON") = py::int_(Bmad::PHOTON);
  m.attr("ELECTRON") = py::int_(Bmad::ELECTRON);
  m.attr("ANTIPROTON") = py::int_(Bmad::ANTIPROTON);
  m.attr("MUON") = py::int_(Bmad::MUON);
  m.attr("PION_MINUS") = py::int_(Bmad::PION_MINUS);
  m.attr("ANTI_DEUTERON") = py::int_(Bmad::ANTI_DEUTERON);
  m.attr("ANTI_NEUTRON") = py::int_(Bmad::ANTI_NEUTRON);
  m.attr("ANTI_REF_PARTICLE") = py::int_(Bmad::ANTI_REF_PARTICLE);
  m.attr("ANTI_HELION") = py::int_(Bmad::ANTI_HELION);
  m.attr("LB_SUBATOMIC") = py::int_(Bmad::LB_SUBATOMIC);
  m.attr("UB_SUBATOMIC") = py::int_(Bmad::UB_SUBATOMIC);
  m.attr("ANTI_ATOM") = py::int_(Bmad::ANTI_ATOM);

  // Enums from sim_utils_struct.f90
  m.attr("INT_GARBAGE") = py::int_(Bmad::INT_GARBAGE);
  m.attr("REAL_GARBAGE") = py::float_(Bmad::REAL_GARBAGE);
  m.attr("INVALID") = py::int_(Bmad::INVALID);
  m.attr("NOT_SET") = py::int_(Bmad::NOT_SET);
  m.attr("X_AXIS") = py::int_(Bmad::X_AXIS);
  m.attr("Y_AXIS") = py::int_(Bmad::Y_AXIS);
  m.attr("Z_AXIS") = py::int_(Bmad::Z_AXIS);
  m.attr("XY_AXIS") = py::int_(Bmad::XY_AXIS);
  m.attr("TRUE_") = py::float_(Bmad::TRUE_);
  m.attr("FALSE_") = py::float_(Bmad::FALSE_);
  m.attr("TRUE_INT") = py::int_(Bmad::TRUE_INT);
  m.attr("FALSE_INT") = py::int_(Bmad::FALSE_INT);
  m.attr("YES") = py::int_(Bmad::YES);
  m.attr("NO") = py::int_(Bmad::NO);
  m.attr("MAYBE") = py::int_(Bmad::MAYBE);
  m.attr("PROVISIONAL") = py::int_(Bmad::PROVISIONAL);

  // Enums from quick_plot_struct.f90
  m.attr("WHITE") = py::int_(Bmad::WHITE);
  m.attr("BLACK") = py::int_(Bmad::BLACK);
  m.attr("RED") = py::int_(Bmad::RED);
  m.attr("GREEN") = py::int_(Bmad::GREEN);
  m.attr("BLUE") = py::int_(Bmad::BLUE);
  m.attr("CYAN") = py::int_(Bmad::CYAN);
  m.attr("MAGENTA") = py::int_(Bmad::MAGENTA);
  m.attr("YELLOW") = py::int_(Bmad::YELLOW);
  m.attr("ORANGE") = py::int_(Bmad::ORANGE);
  m.attr("YELLOW_GREEN") = py::int_(Bmad::YELLOW_GREEN);
  m.attr("LIGHT_GREEN") = py::int_(Bmad::LIGHT_GREEN);
  m.attr("NAVY_BLUE") = py::int_(Bmad::NAVY_BLUE);
  m.attr("PURPLE") = py::int_(Bmad::PURPLE);
  m.attr("REDDISH_PURPLE") = py::int_(Bmad::REDDISH_PURPLE);
  m.attr("DARK_GREY") = py::int_(Bmad::DARK_GREY);
  m.attr("LIGHT_GREY") = py::int_(Bmad::LIGHT_GREY);
  m.attr("TRANSPARENT") = py::int_(Bmad::TRANSPARENT);
  m.attr("SOLID") = py::int_(Bmad::SOLID);
  m.attr("DASHED") = py::int_(Bmad::DASHED);
  m.attr("DASH_DOT") = py::int_(Bmad::DASH_DOT);
  m.attr("DOTTED") = py::int_(Bmad::DOTTED);
  m.attr("DASH_DOT3") = py::int_(Bmad::DASH_DOT3);
  m.attr("SOLID_FILL") = py::int_(Bmad::SOLID_FILL);
  m.attr("NO_FILL") = py::int_(Bmad::NO_FILL);
  m.attr("HATCHED") = py::int_(Bmad::HATCHED);
  m.attr("CROSS_HATCHED") = py::int_(Bmad::CROSS_HATCHED);
  m.attr("SQUARE_SYM") = py::int_(Bmad::SQUARE_SYM);
  m.attr("DOT_SYM") = py::int_(Bmad::DOT_SYM);
  m.attr("PLUS_SYM") = py::int_(Bmad::PLUS_SYM);
  m.attr("TIMES_SYM") = py::int_(Bmad::TIMES_SYM);
  m.attr("CIRCLE_SYM") = py::int_(Bmad::CIRCLE_SYM);
  m.attr("X_SYMBOL_SYM") = py::int_(Bmad::X_SYMBOL_SYM);
  m.attr("TRIANGLE_SYM") = py::int_(Bmad::TRIANGLE_SYM);
  m.attr("CIRCLE_PLUS_SYM") = py::int_(Bmad::CIRCLE_PLUS_SYM);
  m.attr("CIRCLE_DOT_SYM") = py::int_(Bmad::CIRCLE_DOT_SYM);
  m.attr("SQUARE_CONCAVE_SYM") = py::int_(Bmad::SQUARE_CONCAVE_SYM);
  m.attr("DIAMOND_SYM") = py::int_(Bmad::DIAMOND_SYM);
  m.attr("STAR5_SYM") = py::int_(Bmad::STAR5_SYM);
  m.attr("TRIANGLE_FILLED_SYM") = py::int_(Bmad::TRIANGLE_FILLED_SYM);
  m.attr("RED_CROSS_SYM") = py::int_(Bmad::RED_CROSS_SYM);
  m.attr("STAR_OF_DAVID_SYM") = py::int_(Bmad::STAR_OF_DAVID_SYM);
  m.attr("SQUARE_FILLED_SYM") = py::int_(Bmad::SQUARE_FILLED_SYM);
  m.attr("CIRCLE_FILLED_SYM") = py::int_(Bmad::CIRCLE_FILLED_SYM);
  m.attr("STAR5_FILLED_SYM") = py::int_(Bmad::STAR5_FILLED_SYM);
  m.attr("DFLT_DRAW") = py::int_(Bmad::DFLT_DRAW);
  m.attr("DFLT_SET") = py::int_(Bmad::DFLT_SET);
  m.attr("PRINT_PAGE_LONG_LEN") = py::float_(Bmad::PRINT_PAGE_LONG_LEN);
  m.attr("PRINT_PAGE_SHORT_LEN") = py::float_(Bmad::PRINT_PAGE_SHORT_LEN);
  m.attr("FILLED_ARROW_HEAD") = py::int_(Bmad::FILLED_ARROW_HEAD);
  m.attr("OUTLINE_ARROW_HEAD") = py::int_(Bmad::OUTLINE_ARROW_HEAD);
}
