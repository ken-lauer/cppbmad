#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_target_point_struct(
    py::module& m,
    py::class_<TargetPointProxy>& class_);
void init_taylor_struct(py::module& m, py::class_<TaylorProxy>& class_);
void init_taylor_term_struct(
    py::module& m,
    py::class_<TaylorTermProxy>& class_);
void init_track_point_struct(
    py::module& m,
    py::class_<TrackPointProxy>& class_);
void init_track_struct(py::module& m, py::class_<TrackProxy>& class_);
void init_twiss_struct(py::module& m, py::class_<TwissProxy>& class_);
void init_tricubic_cmplx_coef_struct(
    py::module& m,
    py::class_<TricubicCmplxCoefProxy>& class_);
void init_tao_beam_branch_struct(
    py::module& m,
    py::class_<TaoBeamBranchProxy>& class_);
void init_tao_beam_uni_struct(
    py::module& m,
    py::class_<TaoBeamUniProxy>& class_);
void init_tao_building_wall_orientation_struct(
    py::module& m,
    py::class_<TaoBuildingWallOrientationProxy>& class_);
void init_tao_building_wall_point_struct(
    py::module& m,
    py::class_<TaoBuildingWallPointProxy>& class_);
void init_tao_building_wall_section_struct(
    py::module& m,
    py::class_<TaoBuildingWallSectionProxy>& class_);
void init_tao_building_wall_struct(
    py::module& m,
    py::class_<TaoBuildingWallProxy>& class_);
void init_tao_cmd_history_struct(
    py::module& m,
    py::class_<TaoCmdHistoryProxy>& class_);
void init_tao_common_struct(py::module& m, py::class_<TaoCommonProxy>& class_);
void init_tao_curve_color_struct(
    py::module& m,
    py::class_<TaoCurveColorProxy>& class_);
void init_tao_curve_orbit_struct(
    py::module& m,
    py::class_<TaoCurveOrbitProxy>& class_);
void init_tao_curve_struct(py::module& m, py::class_<TaoCurveProxy>& class_);
void init_tao_d1_data_struct(py::module& m, py::class_<TaoD1DataProxy>& class_);
void init_tao_d2_data_struct(py::module& m, py::class_<TaoD2DataProxy>& class_);
void init_tao_data_struct(py::module& m, py::class_<TaoDataProxy>& class_);
void init_tao_data_var_component_struct(
    py::module& m,
    py::class_<TaoDataVarComponentProxy>& class_);
void init_tao_drawing_struct(
    py::module& m,
    py::class_<TaoDrawingProxy>& class_);
void init_tao_dynamic_aperture_struct(
    py::module& m,
    py::class_<TaoDynamicApertureProxy>& class_);
void init_tao_ele_pointer_struct(
    py::module& m,
    py::class_<TaoElePointerProxy>& class_);
void init_tao_ele_shape_struct(
    py::module& m,
    py::class_<TaoEleShapeProxy>& class_);
void init_tao_floor_plan_struct(
    py::module& m,
    py::class_<TaoFloorPlanProxy>& class_);
void init_tao_global_struct(py::module& m, py::class_<TaoGlobalProxy>& class_);
void init_tao_graph_struct(py::module& m, py::class_<TaoGraphProxy>& class_);
void init_tao_histogram_struct(
    py::module& m,
    py::class_<TaoHistogramProxy>& class_);
void init_tao_init_struct(py::module& m, py::class_<TaoInitProxy>& class_);
void init_tao_lat_sigma_struct(
    py::module& m,
    py::class_<TaoLatSigmaProxy>& class_);
void init_tao_lattice_branch_struct(
    py::module& m,
    py::class_<TaoLatticeBranchProxy>& class_);
void init_tao_lattice_struct(
    py::module& m,
    py::class_<TaoLatticeProxy>& class_);
void init_tao_model_branch_struct(
    py::module& m,
    py::class_<TaoModelBranchProxy>& class_);
void init_tao_model_element_struct(
    py::module& m,
    py::class_<TaoModelElementProxy>& class_);
void init_tao_ping_scale_struct(
    py::module& m,
    py::class_<TaoPingScaleProxy>& class_);
void init_tao_plot_cache_struct(
    py::module& m,
    py::class_<TaoPlotCacheProxy>& class_);
void init_tao_plot_page_struct(
    py::module& m,
    py::class_<TaoPlotPageProxy>& class_);
void init_tao_plot_region_struct(
    py::module& m,
    py::class_<TaoPlotRegionProxy>& class_);
void init_tao_plot_struct(py::module& m, py::class_<TaoPlotProxy>& class_);
void init_tao_shape_pattern_point_struct(
    py::module& m,
    py::class_<TaoShapePatternPointProxy>& class_);
void init_tao_shape_pattern_struct(
    py::module& m,
    py::class_<TaoShapePatternProxy>& class_);
void init_tao_spin_dn_dpz_struct(
    py::module& m,
    py::class_<TaoSpinDnDpzProxy>& class_);
void init_tao_spin_ele_struct(
    py::module& m,
    py::class_<TaoSpinEleProxy>& class_);
void init_tao_spin_map_struct(
    py::module& m,
    py::class_<TaoSpinMapProxy>& class_);
void init_tao_spin_polarization_struct(
    py::module& m,
    py::class_<TaoSpinPolarizationProxy>& class_);
void init_tao_super_universe_struct(
    py::module& m,
    py::class_<TaoSuperUniverseProxy>& class_);
void init_tao_title_struct(py::module& m, py::class_<TaoTitleProxy>& class_);
void init_tao_universe_calc_struct(
    py::module& m,
    py::class_<TaoUniverseCalcProxy>& class_);
void init_tao_universe_pointer_struct(
    py::module& m,
    py::class_<TaoUniversePointerProxy>& class_);
void init_tao_universe_struct(
    py::module& m,
    py::class_<TaoUniverseProxy>& class_);
void init_tao_v1_var_struct(py::module& m, py::class_<TaoV1VarProxy>& class_);
void init_tao_var_slave_struct(
    py::module& m,
    py::class_<TaoVarSlaveProxy>& class_);
void init_tao_var_struct(py::module& m, py::class_<TaoVarProxy>& class_);
void init_tao_wave_kick_pt_struct(
    py::module& m,
    py::class_<TaoWaveKickPtProxy>& class_);
void init_tao_wave_struct(py::module& m, py::class_<TaoWaveProxy>& class_);
void init_test_sub_struct(py::module& m, py::class_<TestSubProxy>& class_);
void init_test_sub_sub_struct(
    py::module& m,
    py::class_<TestSubSubProxy>& class_);
