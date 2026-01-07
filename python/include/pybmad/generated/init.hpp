#pragma once

#include <functional>

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// #include "bmad_enums.h"
// #include "bmad_proxy_routines.hpp"
// #include "bmad_std_typedef.h"
// #include "converter_templates.h"
#include "bmad/routines.hpp"
#include "pybmad/arrays.hpp"
// #include "tao_proxies.hpp"

// using namespace Bmad;
// using namespace Tao;
namespace py = pybind11;
using namespace pybind11::literals;

namespace Pybmad {
void init_spline_struct(py::module&, py::class_<SplineProxy>&);
void init_spin_polar_struct(py::module&, py::class_<SpinPolarProxy>&);
void init_ac_kicker_time_struct(py::module&, py::class_<AcKickerTimeProxy>&);
void init_ac_kicker_freq_struct(py::module&, py::class_<AcKickerFreqProxy>&);
void init_ac_kicker_struct(py::module&, py::class_<AcKickerProxy>&);
void init_interval1_coef_struct(py::module&, py::class_<Interval1CoefProxy>&);
void init_photon_reflect_table_struct(
    py::module&,
    py::class_<PhotonReflectTableProxy>&);
void init_photon_reflect_surface_struct(
    py::module&,
    py::class_<PhotonReflectSurfaceProxy>&);
void init_coord_struct(py::module&, py::class_<CoordProxy>&);
void init_coord_array_struct(py::module&, py::class_<CoordArrayProxy>&);
void init_bpm_phase_coupling_struct(
    py::module&,
    py::class_<BpmPhaseCouplingProxy>&);
void init_expression_atom_struct(py::module&, py::class_<ExpressionAtomProxy>&);
void init_wake_sr_z_long_struct(py::module&, py::class_<WakeSrZLongProxy>&);
void init_wake_sr_mode_struct(py::module&, py::class_<WakeSrModeProxy>&);
void init_wake_sr_struct(py::module&, py::class_<WakeSrProxy>&);
void init_wake_lr_mode_struct(py::module&, py::class_<WakeLrModeProxy>&);
void init_wake_lr_struct(py::module&, py::class_<WakeLrProxy>&);
void init_lat_ele_loc_struct(py::module&, py::class_<LatEleLocProxy>&);
void init_wake_struct(py::module&, py::class_<WakeProxy>&);
void init_taylor_term_struct(py::module&, py::class_<TaylorTermProxy>&);
void init_taylor_struct(py::module&, py::class_<TaylorProxy>&);
void init_em_taylor_term_struct(py::module&, py::class_<EmTaylorTermProxy>&);
void init_em_taylor_struct(py::module&, py::class_<EmTaylorProxy>&);
void init_cartesian_map_term1_struct(
    py::module&,
    py::class_<CartesianMapTerm1Proxy>&);
void init_cartesian_map_term_struct(
    py::module&,
    py::class_<CartesianMapTermProxy>&);
void init_cartesian_map_struct(py::module&, py::class_<CartesianMapProxy>&);
void init_cylindrical_map_term1_struct(
    py::module&,
    py::class_<CylindricalMapTerm1Proxy>&);
void init_cylindrical_map_term_struct(
    py::module&,
    py::class_<CylindricalMapTermProxy>&);
void init_cylindrical_map_struct(py::module&, py::class_<CylindricalMapProxy>&);
void init_bicubic_cmplx_coef_struct(
    py::module&,
    py::class_<BicubicCmplxCoefProxy>&);
void init_tricubic_cmplx_coef_struct(
    py::module&,
    py::class_<TricubicCmplxCoefProxy>&);
void init_grid_field_pt1_struct(py::module&, py::class_<GridFieldPt1Proxy>&);
void init_grid_field_pt_struct(py::module&, py::class_<GridFieldPtProxy>&);
void init_grid_field_struct(py::module&, py::class_<GridFieldProxy>&);
void init_floor_position_struct(py::module&, py::class_<FloorPositionProxy>&);
void init_high_energy_space_charge_struct(
    py::module&,
    py::class_<HighEnergySpaceChargeProxy>&);
void init_xy_disp_struct(py::module&, py::class_<XyDispProxy>&);
void init_twiss_struct(py::module&, py::class_<TwissProxy>&);
void init_mode3_struct(py::module&, py::class_<Mode3Proxy>&);
void init_bookkeeping_state_struct(
    py::module&,
    py::class_<BookkeepingStateProxy>&);
void init_rad_map_struct(py::module&, py::class_<RadMapProxy>&);
void init_rad_map_ele_struct(py::module&, py::class_<RadMapEleProxy>&);
void init_gen_grad1_struct(py::module&, py::class_<GenGrad1Proxy>&);
void init_gen_grad_map_struct(py::module&, py::class_<GenGradMapProxy>&);
void init_surface_segmented_pt_struct(
    py::module&,
    py::class_<SurfaceSegmentedPtProxy>&);
void init_surface_segmented_struct(
    py::module&,
    py::class_<SurfaceSegmentedProxy>&);
void init_surface_h_misalign_pt_struct(
    py::module&,
    py::class_<SurfaceHMisalignPtProxy>&);
void init_surface_h_misalign_struct(
    py::module&,
    py::class_<SurfaceHMisalignProxy>&);
void init_surface_displacement_pt_struct(
    py::module&,
    py::class_<SurfaceDisplacementPtProxy>&);
void init_surface_displacement_struct(
    py::module&,
    py::class_<SurfaceDisplacementProxy>&);
void init_target_point_struct(py::module&, py::class_<TargetPointProxy>&);
void init_surface_curvature_struct(
    py::module&,
    py::class_<SurfaceCurvatureProxy>&);
void init_photon_target_struct(py::module&, py::class_<PhotonTargetProxy>&);
void init_photon_material_struct(py::module&, py::class_<PhotonMaterialProxy>&);
void init_pixel_pt_struct(py::module&, py::class_<PixelPtProxy>&);
void init_pixel_detec_struct(py::module&, py::class_<PixelDetecProxy>&);
void init_photon_element_struct(py::module&, py::class_<PhotonElementProxy>&);
void init_wall3d_vertex_struct(py::module&, py::class_<Wall3dVertexProxy>&);
void init_wall3d_section_struct(py::module&, py::class_<Wall3dSectionProxy>&);
void init_wall3d_struct(py::module&, py::class_<Wall3dProxy>&);
void init_ramper_lord_struct(py::module&, py::class_<RamperLordProxy>&);
void init_control_struct(py::module&, py::class_<ControlProxy>&);
void init_control_var1_struct(py::module&, py::class_<ControlVar1Proxy>&);
void init_control_ramp1_struct(py::module&, py::class_<ControlRamp1Proxy>&);
void init_controller_struct(py::module&, py::class_<ControllerProxy>&);
void init_ellipse_beam_init_struct(
    py::module&,
    py::class_<EllipseBeamInitProxy>&);
void init_kv_beam_init_struct(py::module&, py::class_<KvBeamInitProxy>&);
void init_grid_beam_init_struct(py::module&, py::class_<GridBeamInitProxy>&);
void init_beam_init_struct(py::module&, py::class_<BeamInitProxy>&);
void init_lat_param_struct(py::module&, py::class_<LatParamProxy>&);
void init_mode_info_struct(py::module&, py::class_<ModeInfoProxy>&);
void init_pre_tracker_struct(py::module&, py::class_<PreTrackerProxy>&);
void init_anormal_mode_struct(py::module&, py::class_<AnormalModeProxy>&);
void init_linac_normal_mode_struct(
    py::module&,
    py::class_<LinacNormalModeProxy>&);
void init_normal_modes_struct(py::module&, py::class_<NormalModesProxy>&);
void init_em_field_struct(py::module&, py::class_<EmFieldProxy>&);
void init_strong_beam_struct(py::module&, py::class_<StrongBeamProxy>&);
void init_track_point_struct(py::module&, py::class_<TrackPointProxy>&);
void init_track_struct(py::module&, py::class_<TrackProxy>&);
void init_space_charge_common_struct(
    py::module&,
    py::class_<SpaceChargeCommonProxy>&);
void init_bmad_common_struct(py::module&, py::class_<BmadCommonProxy>&);
void init_rad_int1_struct(py::module&, py::class_<RadInt1Proxy>&);
void init_rad_int_branch_struct(py::module&, py::class_<RadIntBranchProxy>&);
void init_rad_int_all_ele_struct(py::module&, py::class_<RadIntAllEleProxy>&);
void init_rf_stair_step_struct(py::module&, py::class_<RfStairStepProxy>&);
void init_rf_ele_struct(py::module&, py::class_<RfEleProxy>&);
void init_ele_struct(py::module&, py::class_<EleProxy>&);
void init_complex_taylor_term_struct(
    py::module&,
    py::class_<ComplexTaylorTermProxy>&);
void init_complex_taylor_struct(py::module&, py::class_<ComplexTaylorProxy>&);
void init_branch_struct(py::module&, py::class_<BranchProxy>&);
void init_lat_struct(py::module&, py::class_<LatProxy>&);
void init_bunch_struct(py::module&, py::class_<BunchProxy>&);
void init_bunch_params_struct(py::module&, py::class_<BunchParamsProxy>&);
void init_beam_struct(py::module&, py::class_<BeamProxy>&);
void init_aperture_point_struct(py::module&, py::class_<AperturePointProxy>&);
void init_aperture_param_struct(py::module&, py::class_<ApertureParamProxy>&);
void init_aperture_scan_struct(py::module&, py::class_<ApertureScanProxy>&);
void init_ele_pointer_struct(py::module&, py::class_<ElePointerProxy>&);
void init_expression_tree_struct(py::module&, py::class_<ExpressionTreeProxy>&);
void init_nametable_struct(py::module&, py::class_<NametableProxy>&);
void init_tao_spin_dn_dpz_struct(py::module&, py::class_<TaoSpinDnDpzProxy>&);
void init_resonance_h_struct(py::module&, py::class_<ResonanceHProxy>&);
void init_spin_orbit_map1_struct(py::module&, py::class_<SpinOrbitMap1Proxy>&);
void init_spin_axis_struct(py::module&, py::class_<SpinAxisProxy>&);
void init_ptc_normal_form_struct(py::module&, py::class_<PtcNormalFormProxy>&);
void init_bmad_normal_form_struct(
    py::module&,
    py::class_<BmadNormalFormProxy>&);
void init_bunch_track_struct(py::module&, py::class_<BunchTrackProxy>&);
void init_summation_rdt_struct(py::module&, py::class_<SummationRdtProxy>&);
void init_tao_ele_shape_struct(py::module&, py::class_<TaoEleShapeProxy>&);
void init_tao_ele_pointer_struct(py::module&, py::class_<TaoElePointerProxy>&);
void init_tao_curve_struct(py::module&, py::class_<TaoCurveProxy>&);
void init_tao_curve_color_struct(py::module&, py::class_<TaoCurveColorProxy>&);
void init_tao_curve_orbit_struct(py::module&, py::class_<TaoCurveOrbitProxy>&);
void init_tao_histogram_struct(py::module&, py::class_<TaoHistogramProxy>&);
void init_lat_ele_order1_struct(py::module&, py::class_<LatEleOrder1Proxy>&);
void init_lat_ele_order_array_struct(
    py::module&,
    py::class_<LatEleOrderArrayProxy>&);
void init_tao_lat_sigma_struct(py::module&, py::class_<TaoLatSigmaProxy>&);
void init_tao_spin_ele_struct(py::module&, py::class_<TaoSpinEleProxy>&);
void init_tao_plot_cache_struct(py::module&, py::class_<TaoPlotCacheProxy>&);
void init_tao_spin_polarization_struct(
    py::module&,
    py::class_<TaoSpinPolarizationProxy>&);
void init_tao_lattice_branch_struct(
    py::module&,
    py::class_<TaoLatticeBranchProxy>&);
void init_tao_model_element_struct(
    py::module&,
    py::class_<TaoModelElementProxy>&);
void init_tao_beam_branch_struct(py::module&, py::class_<TaoBeamBranchProxy>&);
void init_tao_d1_data_struct(py::module&, py::class_<TaoD1DataProxy>&);
void init_tao_d2_data_struct(py::module&, py::class_<TaoD2DataProxy>&);
void init_tao_data_var_component_struct(
    py::module&,
    py::class_<TaoDataVarComponentProxy>&);
void init_tao_graph_struct(py::module&, py::class_<TaoGraphProxy>&);
void init_tao_plot_struct(py::module&, py::class_<TaoPlotProxy>&);
void init_tao_plot_region_struct(py::module&, py::class_<TaoPlotRegionProxy>&);
void init_tao_universe_pointer_struct(
    py::module&,
    py::class_<TaoUniversePointerProxy>&);
void init_tao_super_universe_struct(
    py::module&,
    py::class_<TaoSuperUniverseProxy>&);
void init_tao_var_struct(py::module&, py::class_<TaoVarProxy>&);
void init_tao_var_slave_struct(py::module&, py::class_<TaoVarSlaveProxy>&);
void init_tao_lattice_struct(py::module&, py::class_<TaoLatticeProxy>&);
void init_tao_beam_uni_struct(py::module&, py::class_<TaoBeamUniProxy>&);
void init_tao_dynamic_aperture_struct(
    py::module&,
    py::class_<TaoDynamicApertureProxy>&);
void init_tao_model_branch_struct(
    py::module&,
    py::class_<TaoModelBranchProxy>&);
void init_tao_spin_map_struct(py::module&, py::class_<TaoSpinMapProxy>&);
void init_tao_data_struct(py::module&, py::class_<TaoDataProxy>&);
void init_tao_ping_scale_struct(py::module&, py::class_<TaoPingScaleProxy>&);
void init_tao_universe_calc_struct(
    py::module&,
    py::class_<TaoUniverseCalcProxy>&);
void init_lat_ele_order_struct(py::module&, py::class_<LatEleOrderProxy>&);
void init_tao_title_struct(py::module&, py::class_<TaoTitleProxy>&);
void init_qp_rect_struct(py::module&, py::class_<QpRectProxy>&);
void init_tao_drawing_struct(py::module&, py::class_<TaoDrawingProxy>&);
void init_tao_shape_pattern_struct(
    py::module&,
    py::class_<TaoShapePatternProxy>&);
void init_tao_shape_pattern_point_struct(
    py::module&,
    py::class_<TaoShapePatternPointProxy>&);
void init_qp_axis_struct(py::module&, py::class_<QpAxisProxy>&);
void init_qp_legend_struct(py::module&, py::class_<QpLegendProxy>&);
void init_qp_point_struct(py::module&, py::class_<QpPointProxy>&);
void init_qp_line_struct(py::module&, py::class_<QpLineProxy>&);
void init_qp_symbol_struct(py::module&, py::class_<QpSymbolProxy>&);
void init_tao_floor_plan_struct(py::module&, py::class_<TaoFloorPlanProxy>&);
void init_tao_v1_var_struct(py::module&, py::class_<TaoV1VarProxy>&);
void init_tao_global_struct(py::module&, py::class_<TaoGlobalProxy>&);
void init_tao_init_struct(py::module&, py::class_<TaoInitProxy>&);
void init_tao_common_struct(py::module&, py::class_<TaoCommonProxy>&);
void init_tao_plot_page_struct(py::module&, py::class_<TaoPlotPageProxy>&);
void init_tao_building_wall_struct(
    py::module&,
    py::class_<TaoBuildingWallProxy>&);
void init_tao_building_wall_orientation_struct(
    py::module&,
    py::class_<TaoBuildingWallOrientationProxy>&);
void init_tao_building_wall_section_struct(
    py::module&,
    py::class_<TaoBuildingWallSectionProxy>&);
void init_tao_building_wall_point_struct(
    py::module&,
    py::class_<TaoBuildingWallPointProxy>&);
void init_tao_wave_struct(py::module&, py::class_<TaoWaveProxy>&);
void init_tao_wave_kick_pt_struct(py::module&, py::class_<TaoWaveKickPtProxy>&);
void init_tao_cmd_history_struct(py::module&, py::class_<TaoCmdHistoryProxy>&);
void init_tao_universe_struct(py::module&, py::class_<TaoUniverseProxy>&);
void init_mad_energy_struct(py::module&, py::class_<MadEnergyProxy>&);
void init_mad_map_struct(py::module&, py::class_<MadMapProxy>&);
void init_bbu_stage_struct(py::module&, py::class_<BbuStageProxy>&);
void init_bbu_beam_struct(py::module&, py::class_<BbuBeamProxy>&);
void init_bbu_param_struct(py::module&, py::class_<BbuParamProxy>&);
void init_all_encompassing_struct(
    py::module&,
    py::class_<AllEncompassingProxy>&);
void init_test_sub_struct(py::module&, py::class_<TestSubProxy>&);
void init_test_sub_sub_struct(py::module&, py::class_<TestSubSubProxy>&);
} // namespace Pybmad
