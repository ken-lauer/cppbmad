#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_target_point_struct(nb::module_ &m, nb::class_<TargetPointStruct> &class_);
void init_taylor_struct(nb::module_ &m, nb::class_<TaylorStruct> &class_);
void init_taylor_term_struct(nb::module_ &m, nb::class_<TaylorTermStruct> &class_);
void init_track_point_struct(nb::module_ &m, nb::class_<TrackPointStruct> &class_);
void init_track_struct(nb::module_ &m, nb::class_<TrackStruct> &class_);
void init_twiss_struct(nb::module_ &m, nb::class_<TwissStruct> &class_);
void init_tricubic_cmplx_coef_struct(nb::module_ &m, nb::class_<TricubicCmplxCoefStruct> &class_);
void init_tao_beam_branch_struct(nb::module_ &m, nb::class_<TaoBeamBranchStruct> &class_);
void init_tao_beam_uni_struct(nb::module_ &m, nb::class_<TaoBeamUniStruct> &class_);
void init_tao_building_wall_orientation_struct(
    nb::module_ &m,
    nb::class_<TaoBuildingWallOrientationStruct> &class_
);
void init_tao_building_wall_point_struct(
    nb::module_ &m,
    nb::class_<TaoBuildingWallPointStruct> &class_
);
void init_tao_building_wall_section_struct(
    nb::module_ &m,
    nb::class_<TaoBuildingWallSectionStruct> &class_
);
void init_tao_building_wall_struct(nb::module_ &m, nb::class_<TaoBuildingWallStruct> &class_);
void init_tao_cmd_history_struct(nb::module_ &m, nb::class_<TaoCmdHistoryStruct> &class_);
void init_tao_common_struct(nb::module_ &m, nb::class_<TaoCommonStruct> &class_);
void init_tao_curve_color_struct(nb::module_ &m, nb::class_<TaoCurveColorStruct> &class_);
void init_tao_curve_orbit_struct(nb::module_ &m, nb::class_<TaoCurveOrbitStruct> &class_);
void init_tao_curve_struct(nb::module_ &m, nb::class_<TaoCurveStruct> &class_);
void init_tao_d1_data_struct(nb::module_ &m, nb::class_<TaoD1DataStruct> &class_);
void init_tao_d2_data_struct(nb::module_ &m, nb::class_<TaoD2DataStruct> &class_);
void init_tao_data_struct(nb::module_ &m, nb::class_<TaoDataStruct> &class_);
void init_tao_data_var_component_struct(
    nb::module_ &m,
    nb::class_<TaoDataVarComponentStruct> &class_
);
void init_tao_drawing_struct(nb::module_ &m, nb::class_<TaoDrawingStruct> &class_);
void init_tao_dynamic_aperture_struct(nb::module_ &m, nb::class_<TaoDynamicApertureStruct> &class_);
void init_tao_ele_pointer_struct(nb::module_ &m, nb::class_<TaoElePointerStruct> &class_);
void init_tao_ele_shape_struct(nb::module_ &m, nb::class_<TaoEleShapeStruct> &class_);
void init_tao_eval_node_struct(nb::module_ &m, nb::class_<TaoEvalNodeStruct> &class_);
void init_tao_expression_info_struct(nb::module_ &m, nb::class_<TaoExpressionInfoStruct> &class_);
void init_tao_floor_plan_struct(nb::module_ &m, nb::class_<TaoFloorPlanStruct> &class_);
void init_tao_global_struct(nb::module_ &m, nb::class_<TaoGlobalStruct> &class_);
void init_tao_graph_struct(nb::module_ &m, nb::class_<TaoGraphStruct> &class_);
void init_tao_histogram_struct(nb::module_ &m, nb::class_<TaoHistogramStruct> &class_);
void init_tao_init_struct(nb::module_ &m, nb::class_<TaoInitStruct> &class_);
void init_tao_lat_sigma_struct(nb::module_ &m, nb::class_<TaoLatSigmaStruct> &class_);
void init_tao_lattice_branch_struct(nb::module_ &m, nb::class_<TaoLatticeBranchStruct> &class_);
void init_tao_lattice_struct(nb::module_ &m, nb::class_<TaoLatticeStruct> &class_);
void init_tao_model_branch_struct(nb::module_ &m, nb::class_<TaoModelBranchStruct> &class_);
void init_tao_model_element_struct(nb::module_ &m, nb::class_<TaoModelElementStruct> &class_);
void init_tao_ping_scale_struct(nb::module_ &m, nb::class_<TaoPingScaleStruct> &class_);
void init_tao_plot_cache_struct(nb::module_ &m, nb::class_<TaoPlotCacheStruct> &class_);
void init_tao_plot_page_struct(nb::module_ &m, nb::class_<TaoPlotPageStruct> &class_);
void init_tao_plot_region_struct(nb::module_ &m, nb::class_<TaoPlotRegionStruct> &class_);
void init_tao_plot_struct(nb::module_ &m, nb::class_<TaoPlotStruct> &class_);
void init_tao_shape_pattern_point_struct(
    nb::module_ &m,
    nb::class_<TaoShapePatternPointStruct> &class_
);
void init_tao_shape_pattern_struct(nb::module_ &m, nb::class_<TaoShapePatternStruct> &class_);
void init_tao_spin_dn_dpz_struct(nb::module_ &m, nb::class_<TaoSpinDnDpzStruct> &class_);
void init_tao_spin_ele_struct(nb::module_ &m, nb::class_<TaoSpinEleStruct> &class_);
void init_tao_spin_map_struct(nb::module_ &m, nb::class_<TaoSpinMapStruct> &class_);
void init_tao_spin_polarization_struct(
    nb::module_ &m,
    nb::class_<TaoSpinPolarizationStruct> &class_
);
void init_tao_super_universe_struct(nb::module_ &m, nb::class_<TaoSuperUniverseStruct> &class_);
void init_tao_title_struct(nb::module_ &m, nb::class_<TaoTitleStruct> &class_);
void init_tao_universe_calc_struct(nb::module_ &m, nb::class_<TaoUniverseCalcStruct> &class_);
void init_tao_universe_pointer_struct(nb::module_ &m, nb::class_<TaoUniversePointerStruct> &class_);
void init_tao_universe_struct(nb::module_ &m, nb::class_<TaoUniverseStruct> &class_);
void init_tao_v1_var_struct(nb::module_ &m, nb::class_<TaoV1VarStruct> &class_);
void init_tao_var_slave_struct(nb::module_ &m, nb::class_<TaoVarSlaveStruct> &class_);
void init_tao_var_struct(nb::module_ &m, nb::class_<TaoVarStruct> &class_);
void init_tao_wave_kick_pt_struct(nb::module_ &m, nb::class_<TaoWaveKickPtStruct> &class_);
void init_tao_wave_struct(nb::module_ &m, nb::class_<TaoWaveStruct> &class_);
void init_test_sub_struct(nb::module_ &m, nb::class_<TestSubStruct> &class_);
void init_test_sub_sub_struct(nb::module_ &m, nb::class_<TestSubSubStruct> &class_);
