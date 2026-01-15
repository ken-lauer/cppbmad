#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_target_point_struct(py::module &m, py::class_<TargetPointStruct> &class_);
void init_taylor_struct(py::module &m, py::class_<TaylorStruct> &class_);
void init_taylor_term_struct(py::module &m, py::class_<TaylorTermStruct> &class_);
void init_track_point_struct(py::module &m, py::class_<TrackPointStruct> &class_);
void init_track_struct(py::module &m, py::class_<TrackStruct> &class_);
void init_twiss_struct(py::module &m, py::class_<TwissStruct> &class_);
void init_tricubic_cmplx_coef_struct(py::module &m, py::class_<TricubicCmplxCoefStruct> &class_);
void init_tao_beam_branch_struct(py::module &m, py::class_<TaoBeamBranchStruct> &class_);
void init_tao_beam_uni_struct(py::module &m, py::class_<TaoBeamUniStruct> &class_);
void init_tao_building_wall_orientation_struct(
    py::module &m,
    py::class_<TaoBuildingWallOrientationStruct> &class_
);
void init_tao_building_wall_point_struct(
    py::module &m,
    py::class_<TaoBuildingWallPointStruct> &class_
);
void init_tao_building_wall_section_struct(
    py::module &m,
    py::class_<TaoBuildingWallSectionStruct> &class_
);
void init_tao_building_wall_struct(py::module &m, py::class_<TaoBuildingWallStruct> &class_);
void init_tao_cmd_history_struct(py::module &m, py::class_<TaoCmdHistoryStruct> &class_);
void init_tao_common_struct(py::module &m, py::class_<TaoCommonStruct> &class_);
void init_tao_curve_color_struct(py::module &m, py::class_<TaoCurveColorStruct> &class_);
void init_tao_curve_orbit_struct(py::module &m, py::class_<TaoCurveOrbitStruct> &class_);
void init_tao_curve_struct(py::module &m, py::class_<TaoCurveStruct> &class_);
void init_tao_d1_data_struct(py::module &m, py::class_<TaoD1DataStruct> &class_);
void init_tao_d2_data_struct(py::module &m, py::class_<TaoD2DataStruct> &class_);
void init_tao_data_struct(py::module &m, py::class_<TaoDataStruct> &class_);
void init_tao_data_var_component_struct(
    py::module &m,
    py::class_<TaoDataVarComponentStruct> &class_
);
void init_tao_drawing_struct(py::module &m, py::class_<TaoDrawingStruct> &class_);
void init_tao_dynamic_aperture_struct(py::module &m, py::class_<TaoDynamicApertureStruct> &class_);
void init_tao_ele_pointer_struct(py::module &m, py::class_<TaoElePointerStruct> &class_);
void init_tao_ele_shape_struct(py::module &m, py::class_<TaoEleShapeStruct> &class_);
void init_tao_eval_node_struct(py::module &m, py::class_<TaoEvalNodeStruct> &class_);
void init_tao_expression_info_struct(py::module &m, py::class_<TaoExpressionInfoStruct> &class_);
void init_tao_floor_plan_struct(py::module &m, py::class_<TaoFloorPlanStruct> &class_);
void init_tao_global_struct(py::module &m, py::class_<TaoGlobalStruct> &class_);
void init_tao_graph_struct(py::module &m, py::class_<TaoGraphStruct> &class_);
void init_tao_histogram_struct(py::module &m, py::class_<TaoHistogramStruct> &class_);
void init_tao_init_struct(py::module &m, py::class_<TaoInitStruct> &class_);
void init_tao_lat_sigma_struct(py::module &m, py::class_<TaoLatSigmaStruct> &class_);
void init_tao_lattice_branch_struct(py::module &m, py::class_<TaoLatticeBranchStruct> &class_);
void init_tao_lattice_struct(py::module &m, py::class_<TaoLatticeStruct> &class_);
void init_tao_model_branch_struct(py::module &m, py::class_<TaoModelBranchStruct> &class_);
void init_tao_model_element_struct(py::module &m, py::class_<TaoModelElementStruct> &class_);
void init_tao_ping_scale_struct(py::module &m, py::class_<TaoPingScaleStruct> &class_);
void init_tao_plot_cache_struct(py::module &m, py::class_<TaoPlotCacheStruct> &class_);
void init_tao_plot_page_struct(py::module &m, py::class_<TaoPlotPageStruct> &class_);
void init_tao_plot_region_struct(py::module &m, py::class_<TaoPlotRegionStruct> &class_);
void init_tao_plot_struct(py::module &m, py::class_<TaoPlotStruct> &class_);
void init_tao_shape_pattern_point_struct(
    py::module &m,
    py::class_<TaoShapePatternPointStruct> &class_
);
void init_tao_shape_pattern_struct(py::module &m, py::class_<TaoShapePatternStruct> &class_);
void init_tao_spin_dn_dpz_struct(py::module &m, py::class_<TaoSpinDnDpzStruct> &class_);
void init_tao_spin_ele_struct(py::module &m, py::class_<TaoSpinEleStruct> &class_);
void init_tao_spin_map_struct(py::module &m, py::class_<TaoSpinMapStruct> &class_);
void init_tao_spin_polarization_struct(
    py::module &m,
    py::class_<TaoSpinPolarizationStruct> &class_
);
void init_tao_super_universe_struct(py::module &m, py::class_<TaoSuperUniverseStruct> &class_);
void init_tao_title_struct(py::module &m, py::class_<TaoTitleStruct> &class_);
void init_tao_universe_calc_struct(py::module &m, py::class_<TaoUniverseCalcStruct> &class_);
void init_tao_universe_pointer_struct(py::module &m, py::class_<TaoUniversePointerStruct> &class_);
void init_tao_universe_struct(py::module &m, py::class_<TaoUniverseStruct> &class_);
void init_tao_v1_var_struct(py::module &m, py::class_<TaoV1VarStruct> &class_);
void init_tao_var_slave_struct(py::module &m, py::class_<TaoVarSlaveStruct> &class_);
void init_tao_var_struct(py::module &m, py::class_<TaoVarStruct> &class_);
void init_tao_wave_kick_pt_struct(py::module &m, py::class_<TaoWaveKickPtStruct> &class_);
void init_tao_wave_struct(py::module &m, py::class_<TaoWaveStruct> &class_);
void init_test_sub_struct(py::module &m, py::class_<TestSubStruct> &class_);
void init_test_sub_sub_struct(py::module &m, py::class_<TestSubSubStruct> &class_);
