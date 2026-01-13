#include "pybmad/generated/Bmad_routines_e.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyEAccelField python_e_accel_field(
    EleProxy& ele,
    int voltage_or_gradient,
    std::optional<bool> bmad_standard_tracking,
    double field) {
  Bmad::e_accel_field(ele, voltage_or_gradient, bmad_standard_tracking, field);
  auto py_result{PyEAccelField{field}};
  return py_result;
}
PyEleFullName python_ele_full_name(
    EleProxy& ele,
    std::optional<std::string> template_,
    std::string str) {
  Bmad::ele_full_name(ele, template_, str);
  auto py_result{PyEleFullName{str}};
  return py_result;
}
PyEleHasConstantDsDtRef python_ele_has_constant_ds_dt_ref(
    EleProxy& ele,
    bool is_const) {
  Bmad::ele_has_constant_ds_dt_ref(ele, is_const);
  auto py_result{PyEleHasConstantDsDtRef{is_const}};
  return py_result;
}
PyEleHasNonzeroKick python_ele_has_nonzero_kick(bool has_kick) {
  auto _result = Bmad::ele_has_nonzero_kick(has_kick);
  auto py_result{PyEleHasNonzeroKick{_result, has_kick}};
  return py_result;
}
PyEleHasNonzeroOffset python_ele_has_nonzero_offset(
    EleProxy& ele,
    bool has_offset) {
  Bmad::ele_has_nonzero_offset(ele, has_offset);
  auto py_result{PyEleHasNonzeroOffset{has_offset}};
  return py_result;
}
PyEleLocName python_ele_loc_name(
    EleProxy& ele,
    std::optional<bool> show_branch0,
    std::optional<std::string> parens,
    std::string str) {
  Bmad::ele_loc_name(ele, show_branch0, parens, str);
  auto py_result{PyEleLocName{str}};
  return py_result;
}
PyEleNametableIndex python_ele_nametable_index(EleProxy& ele, int ix_nt) {
  Bmad::ele_nametable_index(ele, ix_nt);
  auto py_result{PyEleNametableIndex{ix_nt}};
  return py_result;
}
PyEleRfStepIndex python_ele_rf_step_index(
    double E_ref,
    double s_rel,
    EleProxy& ele,
    int ix_step) {
  Bmad::ele_rf_step_index(E_ref, s_rel, ele, ix_step);
  auto py_result{PyEleRfStepIndex{ix_step}};
  return py_result;
}
PyEleUniqueName python_ele_unique_name(
    EleProxy& ele,
    LatEleOrderProxy& order,
    std::string unique_name) {
  Bmad::ele_unique_name(ele, order, unique_name);
  auto py_result{PyEleUniqueName{unique_name}};
  return py_result;
}
PyEleValueHasChanged python_ele_value_has_changed(
    EleProxy& ele,
    IntAlloc1D& list,
    RealAlloc1D& abs_tol,
    bool set_old,
    bool has_changed) {
  Bmad::ele_value_has_changed(ele, list, abs_tol, set_old, has_changed);
  auto py_result{PyEleValueHasChanged{has_changed}};
  return py_result;
}
PyEmFieldDerivatives python_em_field_derivatives(
    EleProxy& ele,
    LatParamProxy& param,
    double s_pos,
    CoordProxy& orbit,
    bool local_ref_frame,
    std::optional<bool> grid_allow_s_out_of_bounds = std::nullopt,
    std::optional<double> rf_time = std::nullopt) {
  auto _result = Bmad::em_field_derivatives(
      ele,
      param,
      s_pos,
      orbit,
      local_ref_frame,
      make_opt_ref(grid_allow_s_out_of_bounds),
      make_opt_ref(rf_time));
  auto py_result{PyEmFieldDerivatives{
      _result, s_pos, local_ref_frame, grid_allow_s_out_of_bounds, rf_time}};
  return py_result;
}
PyEnteringElement python_entering_element(
    CoordProxy& orbit,
    int particle_at,
    bool is_entering) {
  Bmad::entering_element(orbit, particle_at, is_entering);
  auto py_result{PyEnteringElement{is_entering}};
  return py_result;
}
PyEqAcKicker python_eq_ac_kicker(
    AcKickerProxy& f1,
    AcKickerProxy& f2,
    bool is_eq) {
  Bmad::eq_ac_kicker(f1, f2, is_eq);
  auto py_result{PyEqAcKicker{is_eq}};
  return py_result;
}
PyEqAcKickerFreq python_eq_ac_kicker_freq(
    AcKickerFreqProxy& f1,
    AcKickerFreqProxy& f2,
    bool is_eq) {
  Bmad::eq_ac_kicker_freq(f1, f2, is_eq);
  auto py_result{PyEqAcKickerFreq{is_eq}};
  return py_result;
}
PyEqAcKickerTime python_eq_ac_kicker_time(
    AcKickerTimeProxy& f1,
    AcKickerTimeProxy& f2,
    bool is_eq) {
  Bmad::eq_ac_kicker_time(f1, f2, is_eq);
  auto py_result{PyEqAcKickerTime{is_eq}};
  return py_result;
}
PyEqAnormalMode python_eq_anormal_mode(
    AnormalModeProxy& f1,
    AnormalModeProxy& f2,
    bool is_eq) {
  Bmad::eq_anormal_mode(f1, f2, is_eq);
  auto py_result{PyEqAnormalMode{is_eq}};
  return py_result;
}
PyEqApertureParam python_eq_aperture_param(
    ApertureParamProxy& f1,
    ApertureParamProxy& f2,
    bool is_eq) {
  Bmad::eq_aperture_param(f1, f2, is_eq);
  auto py_result{PyEqApertureParam{is_eq}};
  return py_result;
}
PyEqAperturePoint python_eq_aperture_point(
    AperturePointProxy& f1,
    AperturePointProxy& f2,
    bool is_eq) {
  Bmad::eq_aperture_point(f1, f2, is_eq);
  auto py_result{PyEqAperturePoint{is_eq}};
  return py_result;
}
PyEqApertureScan python_eq_aperture_scan(
    ApertureScanProxy& f1,
    ApertureScanProxy& f2,
    bool is_eq) {
  Bmad::eq_aperture_scan(f1, f2, is_eq);
  auto py_result{PyEqApertureScan{is_eq}};
  return py_result;
}
PyEqBeam python_eq_beam(BeamProxy& f1, BeamProxy& f2, bool is_eq) {
  Bmad::eq_beam(f1, f2, is_eq);
  auto py_result{PyEqBeam{is_eq}};
  return py_result;
}
PyEqBeamInit python_eq_beam_init(
    BeamInitProxy& f1,
    BeamInitProxy& f2,
    bool is_eq) {
  Bmad::eq_beam_init(f1, f2, is_eq);
  auto py_result{PyEqBeamInit{is_eq}};
  return py_result;
}
PyEqBmadCommon python_eq_bmad_common(
    BmadCommonProxy& f1,
    BmadCommonProxy& f2,
    bool is_eq) {
  Bmad::eq_bmad_common(f1, f2, is_eq);
  auto py_result{PyEqBmadCommon{is_eq}};
  return py_result;
}
PyEqBookkeepingState python_eq_bookkeeping_state(
    BookkeepingStateProxy& f1,
    BookkeepingStateProxy& f2,
    bool is_eq) {
  Bmad::eq_bookkeeping_state(f1, f2, is_eq);
  auto py_result{PyEqBookkeepingState{is_eq}};
  return py_result;
}
PyEqBpmPhaseCoupling python_eq_bpm_phase_coupling(
    BpmPhaseCouplingProxy& f1,
    BpmPhaseCouplingProxy& f2,
    bool is_eq) {
  Bmad::eq_bpm_phase_coupling(f1, f2, is_eq);
  auto py_result{PyEqBpmPhaseCoupling{is_eq}};
  return py_result;
}
PyEqBranch python_eq_branch(BranchProxy& f1, BranchProxy& f2, bool is_eq) {
  Bmad::eq_branch(f1, f2, is_eq);
  auto py_result{PyEqBranch{is_eq}};
  return py_result;
}
PyEqBunch python_eq_bunch(BunchProxy& f1, BunchProxy& f2, bool is_eq) {
  Bmad::eq_bunch(f1, f2, is_eq);
  auto py_result{PyEqBunch{is_eq}};
  return py_result;
}
PyEqBunchParams python_eq_bunch_params(
    BunchParamsProxy& f1,
    BunchParamsProxy& f2,
    bool is_eq) {
  Bmad::eq_bunch_params(f1, f2, is_eq);
  auto py_result{PyEqBunchParams{is_eq}};
  return py_result;
}
PyEqCartesianMap python_eq_cartesian_map(
    CartesianMapProxy& f1,
    CartesianMapProxy& f2,
    bool is_eq) {
  Bmad::eq_cartesian_map(f1, f2, is_eq);
  auto py_result{PyEqCartesianMap{is_eq}};
  return py_result;
}
PyEqCartesianMapTerm python_eq_cartesian_map_term(
    CartesianMapTermProxy& f1,
    CartesianMapTermProxy& f2,
    bool is_eq) {
  Bmad::eq_cartesian_map_term(f1, f2, is_eq);
  auto py_result{PyEqCartesianMapTerm{is_eq}};
  return py_result;
}
PyEqCartesianMapTerm1 python_eq_cartesian_map_term1(
    CartesianMapTerm1Proxy& f1,
    CartesianMapTerm1Proxy& f2,
    bool is_eq) {
  Bmad::eq_cartesian_map_term1(f1, f2, is_eq);
  auto py_result{PyEqCartesianMapTerm1{is_eq}};
  return py_result;
}
PyEqComplexTaylor python_eq_complex_taylor(
    ComplexTaylorProxy& f1,
    ComplexTaylorProxy& f2,
    bool is_eq) {
  Bmad::eq_complex_taylor(f1, f2, is_eq);
  auto py_result{PyEqComplexTaylor{is_eq}};
  return py_result;
}
PyEqComplexTaylorTerm python_eq_complex_taylor_term(
    ComplexTaylorTermProxy& f1,
    ComplexTaylorTermProxy& f2,
    bool is_eq) {
  Bmad::eq_complex_taylor_term(f1, f2, is_eq);
  auto py_result{PyEqComplexTaylorTerm{is_eq}};
  return py_result;
}
PyEqControl python_eq_control(ControlProxy& f1, ControlProxy& f2, bool is_eq) {
  Bmad::eq_control(f1, f2, is_eq);
  auto py_result{PyEqControl{is_eq}};
  return py_result;
}
PyEqControlRamp1 python_eq_control_ramp1(
    ControlRamp1Proxy& f1,
    ControlRamp1Proxy& f2,
    bool is_eq) {
  Bmad::eq_control_ramp1(f1, f2, is_eq);
  auto py_result{PyEqControlRamp1{is_eq}};
  return py_result;
}
PyEqControlVar1 python_eq_control_var1(
    ControlVar1Proxy& f1,
    ControlVar1Proxy& f2,
    bool is_eq) {
  Bmad::eq_control_var1(f1, f2, is_eq);
  auto py_result{PyEqControlVar1{is_eq}};
  return py_result;
}
PyEqController python_eq_controller(
    ControllerProxy& f1,
    ControllerProxy& f2,
    bool is_eq) {
  Bmad::eq_controller(f1, f2, is_eq);
  auto py_result{PyEqController{is_eq}};
  return py_result;
}
PyEqCoord python_eq_coord(CoordProxy& f1, CoordProxy& f2, bool is_eq) {
  Bmad::eq_coord(f1, f2, is_eq);
  auto py_result{PyEqCoord{is_eq}};
  return py_result;
}
PyEqCoordArray python_eq_coord_array(
    CoordArrayProxy& f1,
    CoordArrayProxy& f2,
    bool is_eq) {
  Bmad::eq_coord_array(f1, f2, is_eq);
  auto py_result{PyEqCoordArray{is_eq}};
  return py_result;
}
PyEqCylindricalMap python_eq_cylindrical_map(
    CylindricalMapProxy& f1,
    CylindricalMapProxy& f2,
    bool is_eq) {
  Bmad::eq_cylindrical_map(f1, f2, is_eq);
  auto py_result{PyEqCylindricalMap{is_eq}};
  return py_result;
}
PyEqCylindricalMapTerm python_eq_cylindrical_map_term(
    CylindricalMapTermProxy& f1,
    CylindricalMapTermProxy& f2,
    bool is_eq) {
  Bmad::eq_cylindrical_map_term(f1, f2, is_eq);
  auto py_result{PyEqCylindricalMapTerm{is_eq}};
  return py_result;
}
PyEqCylindricalMapTerm1 python_eq_cylindrical_map_term1(
    CylindricalMapTerm1Proxy& f1,
    CylindricalMapTerm1Proxy& f2,
    bool is_eq) {
  Bmad::eq_cylindrical_map_term1(f1, f2, is_eq);
  auto py_result{PyEqCylindricalMapTerm1{is_eq}};
  return py_result;
}
PyEqEle python_eq_ele(EleProxy& f1, EleProxy& f2, bool is_eq) {
  Bmad::eq_ele(f1, f2, is_eq);
  auto py_result{PyEqEle{is_eq}};
  return py_result;
}
PyEqEllipseBeamInit python_eq_ellipse_beam_init(
    EllipseBeamInitProxy& f1,
    EllipseBeamInitProxy& f2,
    bool is_eq) {
  Bmad::eq_ellipse_beam_init(f1, f2, is_eq);
  auto py_result{PyEqEllipseBeamInit{is_eq}};
  return py_result;
}
PyEqEmField python_eq_em_field(EmFieldProxy& f1, EmFieldProxy& f2, bool is_eq) {
  Bmad::eq_em_field(f1, f2, is_eq);
  auto py_result{PyEqEmField{is_eq}};
  return py_result;
}
PyEqEmTaylor python_eq_em_taylor(
    EmTaylorProxy& f1,
    EmTaylorProxy& f2,
    bool is_eq) {
  Bmad::eq_em_taylor(f1, f2, is_eq);
  auto py_result{PyEqEmTaylor{is_eq}};
  return py_result;
}
PyEqEmTaylorTerm python_eq_em_taylor_term(
    EmTaylorTermProxy& f1,
    EmTaylorTermProxy& f2,
    bool is_eq) {
  Bmad::eq_em_taylor_term(f1, f2, is_eq);
  auto py_result{PyEqEmTaylorTerm{is_eq}};
  return py_result;
}
PyEqExpressionAtom python_eq_expression_atom(
    ExpressionAtomProxy& f1,
    ExpressionAtomProxy& f2,
    bool is_eq) {
  Bmad::eq_expression_atom(f1, f2, is_eq);
  auto py_result{PyEqExpressionAtom{is_eq}};
  return py_result;
}
PyEqFloorPosition python_eq_floor_position(
    FloorPositionProxy& f1,
    FloorPositionProxy& f2,
    bool is_eq) {
  Bmad::eq_floor_position(f1, f2, is_eq);
  auto py_result{PyEqFloorPosition{is_eq}};
  return py_result;
}
PyEqGenGrad1 python_eq_gen_grad1(
    GenGrad1Proxy& f1,
    GenGrad1Proxy& f2,
    bool is_eq) {
  Bmad::eq_gen_grad1(f1, f2, is_eq);
  auto py_result{PyEqGenGrad1{is_eq}};
  return py_result;
}
PyEqGenGradMap python_eq_gen_grad_map(
    GenGradMapProxy& f1,
    GenGradMapProxy& f2,
    bool is_eq) {
  Bmad::eq_gen_grad_map(f1, f2, is_eq);
  auto py_result{PyEqGenGradMap{is_eq}};
  return py_result;
}
PyEqGridBeamInit python_eq_grid_beam_init(
    GridBeamInitProxy& f1,
    GridBeamInitProxy& f2,
    bool is_eq) {
  Bmad::eq_grid_beam_init(f1, f2, is_eq);
  auto py_result{PyEqGridBeamInit{is_eq}};
  return py_result;
}
PyEqGridField python_eq_grid_field(
    GridFieldProxy& f1,
    GridFieldProxy& f2,
    bool is_eq) {
  Bmad::eq_grid_field(f1, f2, is_eq);
  auto py_result{PyEqGridField{is_eq}};
  return py_result;
}
PyEqGridFieldPt python_eq_grid_field_pt(
    GridFieldPtProxy& f1,
    GridFieldPtProxy& f2,
    bool is_eq) {
  Bmad::eq_grid_field_pt(f1, f2, is_eq);
  auto py_result{PyEqGridFieldPt{is_eq}};
  return py_result;
}
PyEqGridFieldPt1 python_eq_grid_field_pt1(
    GridFieldPt1Proxy& f1,
    GridFieldPt1Proxy& f2,
    bool is_eq) {
  Bmad::eq_grid_field_pt1(f1, f2, is_eq);
  auto py_result{PyEqGridFieldPt1{is_eq}};
  return py_result;
}
PyEqHighEnergySpaceCharge python_eq_high_energy_space_charge(
    HighEnergySpaceChargeProxy& f1,
    HighEnergySpaceChargeProxy& f2,
    bool is_eq) {
  Bmad::eq_high_energy_space_charge(f1, f2, is_eq);
  auto py_result{PyEqHighEnergySpaceCharge{is_eq}};
  return py_result;
}
PyEqInterval1Coef python_eq_interval1_coef(
    Interval1CoefProxy& f1,
    Interval1CoefProxy& f2,
    bool is_eq) {
  Bmad::eq_interval1_coef(f1, f2, is_eq);
  auto py_result{PyEqInterval1Coef{is_eq}};
  return py_result;
}
PyEqKvBeamInit python_eq_kv_beam_init(
    KvBeamInitProxy& f1,
    KvBeamInitProxy& f2,
    bool is_eq) {
  Bmad::eq_kv_beam_init(f1, f2, is_eq);
  auto py_result{PyEqKvBeamInit{is_eq}};
  return py_result;
}
PyEqLat python_eq_lat(LatProxy& f1, LatProxy& f2, bool is_eq) {
  Bmad::eq_lat(f1, f2, is_eq);
  auto py_result{PyEqLat{is_eq}};
  return py_result;
}
PyEqLatEleLoc python_eq_lat_ele_loc(
    LatEleLocProxy& f1,
    LatEleLocProxy& f2,
    bool is_eq) {
  Bmad::eq_lat_ele_loc(f1, f2, is_eq);
  auto py_result{PyEqLatEleLoc{is_eq}};
  return py_result;
}
PyEqLatParam python_eq_lat_param(
    LatParamProxy& f1,
    LatParamProxy& f2,
    bool is_eq) {
  Bmad::eq_lat_param(f1, f2, is_eq);
  auto py_result{PyEqLatParam{is_eq}};
  return py_result;
}
PyEqLinacNormalMode python_eq_linac_normal_mode(
    LinacNormalModeProxy& f1,
    LinacNormalModeProxy& f2,
    bool is_eq) {
  Bmad::eq_linac_normal_mode(f1, f2, is_eq);
  auto py_result{PyEqLinacNormalMode{is_eq}};
  return py_result;
}
PyEqMode3 python_eq_mode3(Mode3Proxy& f1, Mode3Proxy& f2, bool is_eq) {
  Bmad::eq_mode3(f1, f2, is_eq);
  auto py_result{PyEqMode3{is_eq}};
  return py_result;
}
PyEqModeInfo python_eq_mode_info(
    ModeInfoProxy& f1,
    ModeInfoProxy& f2,
    bool is_eq) {
  Bmad::eq_mode_info(f1, f2, is_eq);
  auto py_result{PyEqModeInfo{is_eq}};
  return py_result;
}
PyEqNormalModes python_eq_normal_modes(
    NormalModesProxy& f1,
    NormalModesProxy& f2,
    bool is_eq) {
  Bmad::eq_normal_modes(f1, f2, is_eq);
  auto py_result{PyEqNormalModes{is_eq}};
  return py_result;
}
PyEqPhotonElement python_eq_photon_element(
    PhotonElementProxy& f1,
    PhotonElementProxy& f2,
    bool is_eq) {
  Bmad::eq_photon_element(f1, f2, is_eq);
  auto py_result{PyEqPhotonElement{is_eq}};
  return py_result;
}
PyEqPhotonMaterial python_eq_photon_material(
    PhotonMaterialProxy& f1,
    PhotonMaterialProxy& f2,
    bool is_eq) {
  Bmad::eq_photon_material(f1, f2, is_eq);
  auto py_result{PyEqPhotonMaterial{is_eq}};
  return py_result;
}
PyEqPhotonReflectSurface python_eq_photon_reflect_surface(
    PhotonReflectSurfaceProxy& f1,
    PhotonReflectSurfaceProxy& f2,
    bool is_eq) {
  Bmad::eq_photon_reflect_surface(f1, f2, is_eq);
  auto py_result{PyEqPhotonReflectSurface{is_eq}};
  return py_result;
}
PyEqPhotonReflectTable python_eq_photon_reflect_table(
    PhotonReflectTableProxy& f1,
    PhotonReflectTableProxy& f2,
    bool is_eq) {
  Bmad::eq_photon_reflect_table(f1, f2, is_eq);
  auto py_result{PyEqPhotonReflectTable{is_eq}};
  return py_result;
}
PyEqPhotonTarget python_eq_photon_target(
    PhotonTargetProxy& f1,
    PhotonTargetProxy& f2,
    bool is_eq) {
  Bmad::eq_photon_target(f1, f2, is_eq);
  auto py_result{PyEqPhotonTarget{is_eq}};
  return py_result;
}
PyEqPixelDetec python_eq_pixel_detec(
    PixelDetecProxy& f1,
    PixelDetecProxy& f2,
    bool is_eq) {
  Bmad::eq_pixel_detec(f1, f2, is_eq);
  auto py_result{PyEqPixelDetec{is_eq}};
  return py_result;
}
PyEqPixelPt python_eq_pixel_pt(PixelPtProxy& f1, PixelPtProxy& f2, bool is_eq) {
  Bmad::eq_pixel_pt(f1, f2, is_eq);
  auto py_result{PyEqPixelPt{is_eq}};
  return py_result;
}
PyEqPreTracker python_eq_pre_tracker(
    PreTrackerProxy& f1,
    PreTrackerProxy& f2,
    bool is_eq) {
  Bmad::eq_pre_tracker(f1, f2, is_eq);
  auto py_result{PyEqPreTracker{is_eq}};
  return py_result;
}
PyEqRadInt1 python_eq_rad_int1(RadInt1Proxy& f1, RadInt1Proxy& f2, bool is_eq) {
  Bmad::eq_rad_int1(f1, f2, is_eq);
  auto py_result{PyEqRadInt1{is_eq}};
  return py_result;
}
PyEqRadIntAllEle python_eq_rad_int_all_ele(
    RadIntAllEleProxy& f1,
    RadIntAllEleProxy& f2,
    bool is_eq) {
  Bmad::eq_rad_int_all_ele(f1, f2, is_eq);
  auto py_result{PyEqRadIntAllEle{is_eq}};
  return py_result;
}
PyEqRadIntBranch python_eq_rad_int_branch(
    RadIntBranchProxy& f1,
    RadIntBranchProxy& f2,
    bool is_eq) {
  Bmad::eq_rad_int_branch(f1, f2, is_eq);
  auto py_result{PyEqRadIntBranch{is_eq}};
  return py_result;
}
PyEqRadMap python_eq_rad_map(RadMapProxy& f1, RadMapProxy& f2, bool is_eq) {
  Bmad::eq_rad_map(f1, f2, is_eq);
  auto py_result{PyEqRadMap{is_eq}};
  return py_result;
}
PyEqRadMapEle python_eq_rad_map_ele(
    RadMapEleProxy& f1,
    RadMapEleProxy& f2,
    bool is_eq) {
  Bmad::eq_rad_map_ele(f1, f2, is_eq);
  auto py_result{PyEqRadMapEle{is_eq}};
  return py_result;
}
PyEqRamperLord python_eq_ramper_lord(
    RamperLordProxy& f1,
    RamperLordProxy& f2,
    bool is_eq) {
  Bmad::eq_ramper_lord(f1, f2, is_eq);
  auto py_result{PyEqRamperLord{is_eq}};
  return py_result;
}
PyEqSpaceChargeCommon python_eq_space_charge_common(
    SpaceChargeCommonProxy& f1,
    SpaceChargeCommonProxy& f2,
    bool is_eq) {
  Bmad::eq_space_charge_common(f1, f2, is_eq);
  auto py_result{PyEqSpaceChargeCommon{is_eq}};
  return py_result;
}
PyEqSpinPolar python_eq_spin_polar(
    SpinPolarProxy& f1,
    SpinPolarProxy& f2,
    bool is_eq) {
  Bmad::eq_spin_polar(f1, f2, is_eq);
  auto py_result{PyEqSpinPolar{is_eq}};
  return py_result;
}
PyEqSpline python_eq_spline(SplineProxy& f1, SplineProxy& f2, bool is_eq) {
  Bmad::eq_spline(f1, f2, is_eq);
  auto py_result{PyEqSpline{is_eq}};
  return py_result;
}
PyEqStrongBeam python_eq_strong_beam(
    StrongBeamProxy& f1,
    StrongBeamProxy& f2,
    bool is_eq) {
  Bmad::eq_strong_beam(f1, f2, is_eq);
  auto py_result{PyEqStrongBeam{is_eq}};
  return py_result;
}
PyEqSurfaceCurvature python_eq_surface_curvature(
    SurfaceCurvatureProxy& f1,
    SurfaceCurvatureProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_curvature(f1, f2, is_eq);
  auto py_result{PyEqSurfaceCurvature{is_eq}};
  return py_result;
}
PyEqSurfaceDisplacement python_eq_surface_displacement(
    SurfaceDisplacementProxy& f1,
    SurfaceDisplacementProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_displacement(f1, f2, is_eq);
  auto py_result{PyEqSurfaceDisplacement{is_eq}};
  return py_result;
}
PyEqSurfaceDisplacementPt python_eq_surface_displacement_pt(
    SurfaceDisplacementPtProxy& f1,
    SurfaceDisplacementPtProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_displacement_pt(f1, f2, is_eq);
  auto py_result{PyEqSurfaceDisplacementPt{is_eq}};
  return py_result;
}
PyEqSurfaceHMisalign python_eq_surface_h_misalign(
    SurfaceHMisalignProxy& f1,
    SurfaceHMisalignProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_h_misalign(f1, f2, is_eq);
  auto py_result{PyEqSurfaceHMisalign{is_eq}};
  return py_result;
}
PyEqSurfaceHMisalignPt python_eq_surface_h_misalign_pt(
    SurfaceHMisalignPtProxy& f1,
    SurfaceHMisalignPtProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_h_misalign_pt(f1, f2, is_eq);
  auto py_result{PyEqSurfaceHMisalignPt{is_eq}};
  return py_result;
}
PyEqSurfaceSegmented python_eq_surface_segmented(
    SurfaceSegmentedProxy& f1,
    SurfaceSegmentedProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_segmented(f1, f2, is_eq);
  auto py_result{PyEqSurfaceSegmented{is_eq}};
  return py_result;
}
PyEqSurfaceSegmentedPt python_eq_surface_segmented_pt(
    SurfaceSegmentedPtProxy& f1,
    SurfaceSegmentedPtProxy& f2,
    bool is_eq) {
  Bmad::eq_surface_segmented_pt(f1, f2, is_eq);
  auto py_result{PyEqSurfaceSegmentedPt{is_eq}};
  return py_result;
}
PyEqTargetPoint python_eq_target_point(
    TargetPointProxy& f1,
    TargetPointProxy& f2,
    bool is_eq) {
  Bmad::eq_target_point(f1, f2, is_eq);
  auto py_result{PyEqTargetPoint{is_eq}};
  return py_result;
}
PyEqTaylor python_eq_taylor(TaylorProxy& f1, TaylorProxy& f2, bool is_eq) {
  Bmad::eq_taylor(f1, f2, is_eq);
  auto py_result{PyEqTaylor{is_eq}};
  return py_result;
}
PyEqTaylorTerm python_eq_taylor_term(
    TaylorTermProxy& f1,
    TaylorTermProxy& f2,
    bool is_eq) {
  Bmad::eq_taylor_term(f1, f2, is_eq);
  auto py_result{PyEqTaylorTerm{is_eq}};
  return py_result;
}
PyEqTrack python_eq_track(TrackProxy& f1, TrackProxy& f2, bool is_eq) {
  Bmad::eq_track(f1, f2, is_eq);
  auto py_result{PyEqTrack{is_eq}};
  return py_result;
}
PyEqTrackPoint python_eq_track_point(
    TrackPointProxy& f1,
    TrackPointProxy& f2,
    bool is_eq) {
  Bmad::eq_track_point(f1, f2, is_eq);
  auto py_result{PyEqTrackPoint{is_eq}};
  return py_result;
}
PyEqTwiss python_eq_twiss(TwissProxy& f1, TwissProxy& f2, bool is_eq) {
  Bmad::eq_twiss(f1, f2, is_eq);
  auto py_result{PyEqTwiss{is_eq}};
  return py_result;
}
PyEqWake python_eq_wake(WakeProxy& f1, WakeProxy& f2, bool is_eq) {
  Bmad::eq_wake(f1, f2, is_eq);
  auto py_result{PyEqWake{is_eq}};
  return py_result;
}
PyEqWakeLr python_eq_wake_lr(WakeLrProxy& f1, WakeLrProxy& f2, bool is_eq) {
  Bmad::eq_wake_lr(f1, f2, is_eq);
  auto py_result{PyEqWakeLr{is_eq}};
  return py_result;
}
PyEqWakeLrMode python_eq_wake_lr_mode(
    WakeLrModeProxy& f1,
    WakeLrModeProxy& f2,
    bool is_eq) {
  Bmad::eq_wake_lr_mode(f1, f2, is_eq);
  auto py_result{PyEqWakeLrMode{is_eq}};
  return py_result;
}
PyEqWakeSr python_eq_wake_sr(WakeSrProxy& f1, WakeSrProxy& f2, bool is_eq) {
  Bmad::eq_wake_sr(f1, f2, is_eq);
  auto py_result{PyEqWakeSr{is_eq}};
  return py_result;
}
PyEqWakeSrMode python_eq_wake_sr_mode(
    WakeSrModeProxy& f1,
    WakeSrModeProxy& f2,
    bool is_eq) {
  Bmad::eq_wake_sr_mode(f1, f2, is_eq);
  auto py_result{PyEqWakeSrMode{is_eq}};
  return py_result;
}
PyEqWakeSrZLong python_eq_wake_sr_z_long(
    WakeSrZLongProxy& f1,
    WakeSrZLongProxy& f2,
    bool is_eq) {
  Bmad::eq_wake_sr_z_long(f1, f2, is_eq);
  auto py_result{PyEqWakeSrZLong{is_eq}};
  return py_result;
}
PyEqWall3d python_eq_wall3d(Wall3dProxy& f1, Wall3dProxy& f2, bool is_eq) {
  Bmad::eq_wall3d(f1, f2, is_eq);
  auto py_result{PyEqWall3d{is_eq}};
  return py_result;
}
PyEqWall3dSection python_eq_wall3d_section(
    Wall3dSectionProxy& f1,
    Wall3dSectionProxy& f2,
    bool is_eq) {
  Bmad::eq_wall3d_section(f1, f2, is_eq);
  auto py_result{PyEqWall3dSection{is_eq}};
  return py_result;
}
PyEqWall3dVertex python_eq_wall3d_vertex(
    Wall3dVertexProxy& f1,
    Wall3dVertexProxy& f2,
    bool is_eq) {
  Bmad::eq_wall3d_vertex(f1, f2, is_eq);
  auto py_result{PyEqWall3dVertex{is_eq}};
  return py_result;
}
PyEqXyDisp python_eq_xy_disp(XyDispProxy& f1, XyDispProxy& f2, bool is_eq) {
  Bmad::eq_xy_disp(f1, f2, is_eq);
  auto py_result{PyEqXyDisp{is_eq}};
  return py_result;
}
PyEqualSignHere python_equal_sign_here(
    EleProxy& ele,
    std::string delim,
    bool is_here) {
  Bmad::equal_sign_here(ele, delim, is_here);
  auto py_result{PyEqualSignHere{delim, is_here}};
  return py_result;
}
PyEquivalentTaylorAttributes python_equivalent_taylor_attributes(
    EleProxy& ele_taylor,
    EleProxy& ele2,
    bool equiv) {
  Bmad::equivalent_taylor_attributes(ele_taylor, ele2, equiv);
  auto py_result{PyEquivalentTaylorAttributes{equiv}};
  return py_result;
}
PyEtdiv python_etdiv(
    double A,
    double B,
    double C,
    double D,
    double E,
    double F) {
  Bmad::etdiv(A, B, C, D, E, F);
  auto py_result{PyEtdiv{A, B, C, D, E, F}};
  return py_result;
}
PyExpectOneOf python_expect_one_of(
    std::string delim_list,
    bool check_input_delim,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool is_ok) {
  Bmad::expect_one_of(
      delim_list, check_input_delim, ele_name, delim, delim_found, is_ok);
  auto py_result{PyExpectOneOf{delim, delim_found, is_ok}};
  return py_result;
}

void init_Bmad_routines_e(py::module& m) {
  py::class_<PyEAccelField, std::unique_ptr<PyEAccelField>>(
      m, "EAccelField", "e_accel_field return type")
      .def_readonly("field", &PyEAccelField::field)
      .def("__len__", [](const PyEAccelField&) { return 1; })
      .def("__getitem__", [](const PyEAccelField& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.field);
        throw py::index_error();
      });
  m.def(
      "e_accel_field",
      &python_e_accel_field,
      py::arg("ele"),
      py::arg("voltage_or_gradient"),
      py::arg("bmad_standard_tracking") = py::none(),
      py::arg("field"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lcavity or rfcavity element.
  voltage_or_gradient : int
      voltage$ or gradient$
  bmad_standard_tracking : bool, optional
      Using bmad_standard tracking? Default is False.
  field : 
  )""");
  m.def(
      "e_crit_photon",
      &Bmad::e_crit_photon,
      py::arg("gamma"),
      py::arg("g_bend"),
      R"""(Function E_crit_photon (gamma, g_bend) result (E_crit)

  Routine to calculate the photon critical energy in a bend.

  Parameters
  ----------
  gamma : float
      Gamma factor of charged particle emitting photon.
  g_bend : float
      1/radius bending strength.

  Returns
  -------
  E_crit : float
      Critical photon energy.
  )""");
  py::class_<Bmad::EigenDecomp6mat, std::unique_ptr<Bmad::EigenDecomp6mat>>(
      m, "EigenDecomp6mat", "eigen_decomp_6mat return type")
      .def_readonly("eval", &Bmad::EigenDecomp6mat::eval)
      .def_readonly("evec", &Bmad::EigenDecomp6mat::evec)
      .def_readonly("err_flag", &Bmad::EigenDecomp6mat::err_flag)
      .def_readonly("tunes", &Bmad::EigenDecomp6mat::tunes)
      .def("__len__", [](const Bmad::EigenDecomp6mat&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::EigenDecomp6mat& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.eval);
            if (i == 1)
              return py::cast(s.evec);
            if (i == 2)
              return py::cast(s.err_flag);
            if (i == 3)
              return py::cast(s.tunes);
            throw py::index_error();
          });
  m.def(
      "eigen_decomp_6mat",
      &Bmad::eigen_decomp_6mat,
      py::arg("mat"),
      R"""(Subroutine eigen_decomp_6mat(mat, eval, evec, tunes, err_flag)

  Compute eigenvalues and eigenvectors of a real 6x6 matrix.
  The evals and evecs are in general complex.

  Parameters
  ----------
  mat : float
      6x6 real matrix.  Usually a transfer matrix or sigma matrix.

  Returns
  -------
  eval : complex
      complex eigenvalues.
  evec : complex
      complex eigenvectors arranged down columns.
  err_flag : bool
      set to true if an error has occured.
  tunes : float
      Mode tunes, in radians.
  )""");
  m.def(
      "ele_compute_ref_energy_and_time",
      &Bmad::ele_compute_ref_energy_and_time,
      py::arg("ele0"),
      py::arg("ele"),
      py::arg("param"),
      py::arg("err_flag"),
      R"""(Parameters
  ----------
  ele0 : EleStruct
      Previous element in lattice with starting energy and time values.
  ele : EleStruct
      Lattice element .time_ref_orb_in  -- Starting orbit for ref time calc.
      This parameter is an input/output and is modified in-place. As an output: Lattice element with reference
      energy and time.
  param : LatParamStruct
      Lattice parameters.
  err_flag : bool
      Set true if there is an error. False otherwise.
  )""");
  m.def(
      "ele_equal_ele",
      &Bmad::ele_equal_ele,
      py::arg("ele_out"),
      py::arg("ele_in"),
      R"""(Parameters
  ----------
  ele_out : 
  ele_in : 
  )""");
  m.def(
      "ele_equals_ele",
      &Bmad::ele_equals_ele,
      py::arg("ele_in"),
      py::arg("update_nametable"),
      R"""(Subroutine ele_equals_ele (ele_out, ele_in, update_nametable)

  Subroutine that is used to set an element equal to another.
  Note: Use ele_equal_ele instead unless you know what you are doing.


  Parameters
  ----------
  ele_in : EleStruct
      Input element.
  update_nametable : bool
      If true, update the nametable. If false, do not. Note: nametable updates can take time if this routine is
      called a many times. See remove_eles_from_lat as an example.

  Returns
  -------
  ele_out : EleStruct
      Output element.
  )""");
  m.def(
      "ele_finalizer",
      &Bmad::ele_finalizer,
      py::arg("ele"),
      R"""(Subroutine ele_finalizer(ele)

  Finalizer routine for ele_struct instances.
  NOTE: Not currently used.

  Parameters
  ----------
  ele : EleStruct
      Element to cleanup.
      This parameter is an input/output and is modified in-place. As an output: Element with pointers
      deallocated as needed.
  )""");
  py::class_<PyEleFullName, std::unique_ptr<PyEleFullName>>(
      m, "EleFullName", "ele_full_name return type")
      .def_readonly("str", &PyEleFullName::str)
      .def("__len__", [](const PyEleFullName&) { return 1; })
      .def("__getitem__", [](const PyEleFullName& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "ele_full_name",
      &python_ele_full_name,
      py::arg("ele"),
      py::arg("template_") = py::none(),
      py::arg("str"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element in a lattice
  template : unknown, optional
      Encoding template. Default is "@N (&#)".
  str : 
  )""");
  m.def(
      "ele_geometry",
      &Bmad::ele_geometry,
      py::arg("floor_start"),
      py::arg("ele"),
      py::arg("len_scale") = py::none(),
      py::arg("ignore_patch_err") = py::none(),
      R"""(Parameters
  ----------
  floor_start : 
      Starting floor coordinates at upstream end. Not used for fiducial and girder elements.
  ele : EleStruct
      Element to propagate the geometry through.
  floor_end : FloorPositionStruct
      Output floor position. If not present then ele.floor will be used and ele.bookkeeping_state.floor_position
      will be set to ok$. .r(3)              -- X, Y, Z Floor position at end of element .w(3,3)            -- W
      matrix corresponding to orientation angles .theta, phi, .psi  -- Orientation angles
  len_scale : float, optional
      factor to scale the length of the element. 1.0_rp => Output is geometry at end of element (default).
      0.5_rp => Output is geometry at center of element. -1.0_rp => Used to propagate geometry in reverse.
  ignore_patch_err : bool, optional
      If present and True, ignore flexible patch errors. This is used by ele_compute_ref_energy_and_time to
      suppress unnecessary messages.
  )""");
  m.def(
      "ele_geometry_with_misalignments",
      &Bmad::ele_geometry_with_misalignments,
      py::arg("ele"),
      py::arg("len_scale") = py::none(),
      py::arg("floor"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element under consideration.
  len_scale : float, optional
      factor to scale the length of the element. 1.0_rp => Output is geometry at end of element (default).
      0.5_rp => Output is geometry at center of element. -1.0_rp => Used to propagate geometry in reverse.
  floor : 
  )""");
  py::class_<PyEleHasConstantDsDtRef, std::unique_ptr<PyEleHasConstantDsDtRef>>(
      m, "EleHasConstantDsDtRef", "ele_has_constant_ds_dt_ref return type")
      .def_readonly("is_const", &PyEleHasConstantDsDtRef::is_const)
      .def("__len__", [](const PyEleHasConstantDsDtRef&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEleHasConstantDsDtRef& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_const);
            throw py::index_error();
          });
  m.def(
      "ele_has_constant_ds_dt_ref",
      &python_ele_has_constant_ds_dt_ref,
      py::arg("ele"),
      py::arg("is_const"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element.
  is_const : 
  )""");
  py::class_<PyEleHasNonzeroKick, std::unique_ptr<PyEleHasNonzeroKick>>(
      m, "EleHasNonzeroKick", "ele_has_nonzero_kick return type")
      .def_readonly("ele", &PyEleHasNonzeroKick::ele)
      .def_readonly("has_kick", &PyEleHasNonzeroKick::has_kick)
      .def("__len__", [](const PyEleHasNonzeroKick&) { return 2; })
      .def(
          "__getitem__", [](const PyEleHasNonzeroKick& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.ele);
            if (i == 1)
              return py::cast(s.has_kick);
            throw py::index_error();
          });
  m.def(
      "ele_has_nonzero_kick",
      &python_ele_has_nonzero_kick,
      py::arg("has_kick"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element with no kicks.
  has_kick : 
  )""");
  py::class_<PyEleHasNonzeroOffset, std::unique_ptr<PyEleHasNonzeroOffset>>(
      m, "EleHasNonzeroOffset", "ele_has_nonzero_offset return type")
      .def_readonly("has_offset", &PyEleHasNonzeroOffset::has_offset)
      .def("__len__", [](const PyEleHasNonzeroOffset&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEleHasNonzeroOffset& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.has_offset);
            throw py::index_error();
          });
  m.def(
      "ele_has_nonzero_offset",
      &python_ele_has_nonzero_offset,
      py::arg("ele"),
      py::arg("has_offset"),
      R"""(Parameters
  ----------
  ele : 
  has_offset : 
  )""");
  m.def(
      "ele_is_monitor",
      &Bmad::ele_is_monitor,
      py::arg("ele"),
      py::arg("print_warning") = py::none(),
      R"""(Function ele_is_monitor (ele, print_warning) result (is_monitor)

  Routine to check that an element is either a detector, instrument, monitor, or marker.
  These are the elements where measurement errors can be defined.

  Parameters
  ----------
  ele : EleStruct
      Lattice element.
  print_warning : bool, optional
      If True print a warning message if the element not a monitor like element. Default is True.

  Returns
  -------
  is_monitor : bool
      Set True if the element is a monitor like element.
  )""");
  m.def(
      "ele_loc",
      &Bmad::ele_loc,
      py::arg("ele"),
      py::arg("loc"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to be identified
  loc : 
  )""");
  py::class_<PyEleLocName, std::unique_ptr<PyEleLocName>>(
      m, "EleLocName", "ele_loc_name return type")
      .def_readonly("str", &PyEleLocName::str)
      .def("__len__", [](const PyEleLocName&) { return 1; })
      .def("__getitem__", [](const PyEleLocName& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "ele_loc_name",
      &python_ele_loc_name,
      py::arg("ele"),
      py::arg("show_branch0") = py::none(),
      py::arg("parens") = py::none(),
      py::arg("str"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element in a lattice
  show_branch0 : bool, optional
      Explicitly show branch for main lattice elements? Default is False.
  parens : unknown, optional
      If present, enclose location string using the two characters supplied. Typically parens will be set to
      "()" or "[]".
  str : 
  )""");
  py::class_<
      Bmad::EleMisalignmentLSCalc,
      std::unique_ptr<Bmad::EleMisalignmentLSCalc>>(
      m, "EleMisalignmentLSCalc", "ele_misalignment_l_s_calc return type")
      .def_readonly("L_mis", &Bmad::EleMisalignmentLSCalc::L_mis)
      .def_readonly("S_mis", &Bmad::EleMisalignmentLSCalc::S_mis)
      .def("__len__", [](const Bmad::EleMisalignmentLSCalc&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::EleMisalignmentLSCalc& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.L_mis);
            if (i == 1)
              return py::cast(s.S_mis);
            throw py::index_error();
          });
  m.def(
      "ele_misalignment_l_s_calc",
      &Bmad::ele_misalignment_l_s_calc,
      py::arg("ele"),
      R"""(Parameters
  ----------
  ele : float
      Element
  L_mis : float
      Misalignment vector relative to center of element
  S_mis : float
      Misalignment matrix relative to center of element
  )""");
  py::class_<PyEleNametableIndex, std::unique_ptr<PyEleNametableIndex>>(
      m, "EleNametableIndex", "ele_nametable_index return type")
      .def_readonly("ix_nt", &PyEleNametableIndex::ix_nt)
      .def("__len__", [](const PyEleNametableIndex&) { return 1; })
      .def(
          "__getitem__", [](const PyEleNametableIndex& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.ix_nt);
            throw py::index_error();
          });
  m.def(
      "ele_nametable_index",
      &python_ele_nametable_index,
      py::arg("ele"),
      py::arg("ix_nt"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element in a lattice.
  ix_nt : 
  )""");
  m.def(
      "ele_order_calc",
      &Bmad::ele_order_calc,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Lattice to analyze.
  order : LatEleOrderStruct
      Structure holding the element order information.
  )""");
  m.def(
      "ele_reference_energy_correction",
      &Bmad::ele_reference_energy_correction,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("particle_at"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element being tracked through.
  orbit : CoordStruct
      Coordinates to correct.
  particle_at : int
      first_track_edge$ (that is, entering the element), or second_track_edge$ (that is, leaving the element),
      or upstream_end$ (inherit ele.value(p0c_start$) ref), or downstream_end$ (inherit ele.value(p0c$)).
      inside$ (or anything else) -> Do nothing.
  mat6 : float, optional
      Transfer matrix before correction.
      This parameter is an input/output and is modified in-place. As an output: Transfer matrix transfer matrix
      including correction.
  make_matrix : bool, optional
      Propagate the transfer matrix? Default is false.
  )""");
  py::class_<PyEleRfStepIndex, std::unique_ptr<PyEleRfStepIndex>>(
      m, "EleRfStepIndex", "ele_rf_step_index return type")
      .def_readonly("ix_step", &PyEleRfStepIndex::ix_step)
      .def("__len__", [](const PyEleRfStepIndex&) { return 1; })
      .def("__getitem__", [](const PyEleRfStepIndex& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.ix_step);
        throw py::index_error();
      });
  m.def(
      "ele_rf_step_index",
      &python_ele_rf_step_index,
      py::arg("E_ref"),
      py::arg("s_rel"),
      py::arg("ele"),
      py::arg("ix_step"),
      R"""(Parameters
  ----------
  E_ref : float
      Reference energy of step. If negative, ignore and use s_rel.
  s_rel : float
      S-position relative to the beginning of the element
  ele : float
      RF cavity.
  ix_step : 
  )""");
  py::class_<Bmad::EleToFibre, std::unique_ptr<Bmad::EleToFibre>>(
      m, "EleToFibre", "ele_to_fibre return type")
      .def_readonly("ptc_fibre", &Bmad::EleToFibre::ptc_fibre)
      .def_readonly("err_flag", &Bmad::EleToFibre::err_flag)
      .def("__len__", [](const Bmad::EleToFibre&) { return 2; })
      .def("__getitem__", [](const Bmad::EleToFibre& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ptc_fibre);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "ele_to_fibre",
      &Bmad::ele_to_fibre,
      py::arg("ele"),
      py::arg("use_offsets"),
      py::arg("integ_order") = py::none(),
      py::arg("steps") = py::none(),
      py::arg("for_layout") = py::none(),
      py::arg("ref_in") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Bmad element.
  ptc_fibre : unknown
      PTC fibre element.
  use_offsets : bool
      Does ptc_fibre include element offsets, pitches and tilt?
  err_flag : bool
      Set True if setup OK. False otherwise.
  integ_order : int, optional
      Order for the sympletic integrator. Possibilities are: 2, 4, or 6 Overrides ele.value(integrator_order$).
      default = 2 (if not set with set_ptc).
  steps : int, optional
      Number of integration steps. Overrides ele.value(ds_step$).
  for_layout : bool, optional
      If True then fibre will be put in the PTC layout. Default is False.
  ref_in : CoordStruct, optional
      Particle to be tracked. ref_particle$, electron$, etc. This argument should only be present when the fibre
      is not to be put in a layout.
  )""");
  py::class_<
      Bmad::EleToPtcMagneticBnAn,
      std::unique_ptr<Bmad::EleToPtcMagneticBnAn>>(
      m, "EleToPtcMagneticBnAn", "ele_to_ptc_magnetic_bn_an return type")
      .def_readonly("bn", &Bmad::EleToPtcMagneticBnAn::bn)
      .def_readonly("an", &Bmad::EleToPtcMagneticBnAn::an)
      .def_readonly("n_max", &Bmad::EleToPtcMagneticBnAn::n_max)
      .def("__len__", [](const Bmad::EleToPtcMagneticBnAn&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::EleToPtcMagneticBnAn& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.bn);
            if (i == 1)
              return py::cast(s.an);
            if (i == 2)
              return py::cast(s.n_max);
            throw py::index_error();
          });
  m.def(
      "ele_to_ptc_magnetic_bn_an",
      &Bmad::ele_to_ptc_magnetic_bn_an,
      py::arg("ele"),
      R"""(Subroutine ele_to_ptc_magnetic_bn_an (ele, bn, an, n_max)

  Routine to compute the a(n) and b(n) magnetic multipole components of a magnet.
  This is used to interface between eles and PTC fibres

  Note: The multipole index uses the PTC convention of starting from 1 instead of zero.

  Note: On the PTC side bn(1) is error field when creating a fibre but
  is the total field when the fibre is being modified. This routine returns the error field.

  Parameters
  ----------
  ele : EleStruct
      Bmad Element.

  Returns
  -------
  bn : float
      Normal multipole component.
  an : float
      Skew multipole component.
  n_max : int
      Maximum non-zero multipole component. Set to zero if there are no multipoles.
  )""");
  m.def(
      "ele_to_spin_taylor",
      &Bmad::ele_to_spin_taylor,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orb0"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
      This parameter is an input/output and is modified in-place. As an output: Element with spin map.
  param : unknown
      Branch parameters.
  orb0 : CoordStruct
      Starting ref coords.
  )""");
  py::class_<Bmad::EleToTaylor, std::unique_ptr<Bmad::EleToTaylor>>(
      m, "EleToTaylor", "ele_to_taylor return type")
      .def_readonly("orbital_taylor", &Bmad::EleToTaylor::orbital_taylor)
      .def_readonly("spin_taylor", &Bmad::EleToTaylor::spin_taylor)
      .def("__len__", [](const Bmad::EleToTaylor&) { return 2; })
      .def("__getitem__", [](const Bmad::EleToTaylor& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.orbital_taylor);
        if (i == 1)
          return py::cast(s.spin_taylor);
        throw py::index_error();
      });
  m.def(
      "ele_to_taylor",
      &Bmad::ele_to_taylor,
      py::arg("ele"),
      py::arg("orb0") = py::none(),
      py::arg("taylor_map_includes_offsets") = py::none(),
      py::arg("include_damping") = py::none(),
      R"""(Parameters
  ----------
  ele : ElementStruct
      Element to construct map for.
  orb0 : CoordStruct, optional
      Starting coords around which the Taylor map is evaluated. Default is the zero orbit.
  taylor_map_includes_offsets : bool, optional
      If present then value overrides ele.taylor_map_includes_offsets. -- Logical, optional: If present then
      value overrides ele.taylor_map_includes_offsets.
  include_damping : bool, optional
      Sets if radiation damping is included. Default is what is set in ptc_private.base_state.
  orbital_taylor : TaylorStruct
      Orbital taylor map. If not present then the map is put in ele.taylor.
  spin_taylor : TaylorStruct
      Spin taylor map. If not present then the map is put in ele.spin_taylor.
  )""");
  py::class_<PyEleUniqueName, std::unique_ptr<PyEleUniqueName>>(
      m, "EleUniqueName", "ele_unique_name return type")
      .def_readonly("unique_name", &PyEleUniqueName::unique_name)
      .def("__len__", [](const PyEleUniqueName&) { return 1; })
      .def("__getitem__", [](const PyEleUniqueName& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.unique_name);
        throw py::index_error();
      });
  m.def(
      "ele_unique_name",
      &python_ele_unique_name,
      py::arg("ele"),
      py::arg("order"),
      py::arg("unique_name"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to construct a unique name for.
  order : LatEleOrderStruct
      Information on element ordering. Before calling this routine, use the routine ele_order_calc to compute
      this argument.
  unique_name : 
  )""");
  py::class_<PyEleValueHasChanged, std::unique_ptr<PyEleValueHasChanged>>(
      m, "EleValueHasChanged", "ele_value_has_changed return type")
      .def_readonly("has_changed", &PyEleValueHasChanged::has_changed)
      .def("__len__", [](const PyEleValueHasChanged&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEleValueHasChanged& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.has_changed);
            throw py::index_error();
          });
  m.def(
      "ele_value_has_changed",
      &python_ele_value_has_changed,
      py::arg("ele"),
      py::arg("list"),
      py::arg("abs_tol"),
      py::arg("set_old"),
      py::arg("has_changed"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element under consideration.
      This parameter is an input/output and is modified in-place. As an output: ele.old_value may be set
      depending upon setting of set_old
  list : int
      List of indexes of ele.value(:) array to check.
  abs_tol : float
      List of values such that if the change in parameter value is less than this it is not considered to have
      changed significantly.
  set_old : bool
      If True then set ele.old_value(j) = ele.value(j) for j in list
  has_changed : 
  )""");
  m.def(
      "ele_vec_equal_ele_vec",
      &Bmad::ele_vec_equal_ele_vec,
      py::arg("ele1"),
      py::arg("ele2"),
      R"""(Parameters
  ----------
  ele1 : 
  ele2 : 
  )""");
  py::class_<
      Bmad::ElecMultipoleField,
      std::unique_ptr<Bmad::ElecMultipoleField>>(
      m, "ElecMultipoleField", "elec_multipole_field return type")
      .def_readonly("Ex", &Bmad::ElecMultipoleField::Ex)
      .def_readonly("Ey", &Bmad::ElecMultipoleField::Ey)
      .def_readonly("dE", &Bmad::ElecMultipoleField::dE)
      .def_readonly("compute_dE", &Bmad::ElecMultipoleField::compute_dE)
      .def("__len__", [](const Bmad::ElecMultipoleField&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::ElecMultipoleField& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.Ex);
            if (i == 1)
              return py::cast(s.Ey);
            if (i == 2)
              return py::cast(s.dE);
            if (i == 3)
              return py::cast(s.compute_dE);
            throw py::index_error();
          });
  m.def(
      "elec_multipole_field",
      &Bmad::elec_multipole_field,
      py::arg("a"),
      py::arg("b"),
      py::arg("n"),
      py::arg("coord"),
      R"""(Parameters
  ----------
  a : float
      Multipole skew component.
  b : float
      Multipole normal component.
  n : float
      Multipole order.
  coord : CoordStruct
  Ex : float
      X field component
  Ey : float
      Y field component.
  dE : float
      Field derivatives: dfield(x,y)/d(x,y).
  compute_dE : bool
      If False, do not compute the field derivatives even if dE is present. Default is True.
  )""");
  py::class_<Bmad::ElementAtSBranch, std::unique_ptr<Bmad::ElementAtSBranch>>(
      m, "ElementAtSBranch", "element_at_s_branch return type")
      .def_readonly("err_flag", &Bmad::ElementAtSBranch::err_flag)
      .def_readonly("s_eff", &Bmad::ElementAtSBranch::s_eff)
      .def_readonly("position", &Bmad::ElementAtSBranch::position)
      .def_readonly("ix_ele", &Bmad::ElementAtSBranch::ix_ele)
      .def("__len__", [](const Bmad::ElementAtSBranch&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::ElementAtSBranch& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.s_eff);
            if (i == 2)
              return py::cast(s.position);
            if (i == 3)
              return py::cast(s.ix_ele);
            throw py::index_error();
          });
  m.def(
      "element_at_s",
      py::overload_cast<BranchProxy&, double, bool, std::optional<bool>>(
          &Bmad::element_at_s),
      py::arg("branch"),
      py::arg("s"),
      py::arg("choose_max"),
      py::arg("print_err") = py::none(),
      R"""(Function element_at_s (...) result (ix_ele)

  Function to return the index of the element at position s.

  element_at_s is an overloaded name for:
    function element_at_s_lat (lat, s, choose_max, ix_branch, err_flag, s_eff, position, print_err) result (ix_ele)
    function element_at_s_branch (branch, s, choose_max, err_flag, s_eff, position, print_err) result (ix_ele)

  The differnce between these two routine is that with element_at_s_lat, the branch is given by the lat
    and ix_ele arguments: branch = lat%branch(ix_ele). With element_at_s_branch, the branch is an argument.

  Also see: pointer_to_element_at_s

  ix_ele is choisen such that:
  If choose_max = True:
      If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
      Else: branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
  If choose_max = False:
      If s = branch%ele(0)%s: ix_ele = 0
      Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
  That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
      choose_max = True  => ix_ele = ix2
      choose_max = False => ix_ele = ix1

  The setting of choose_max only makes a difference when s corresponds to an element boundary.

  Note: For a circular lattice, s is evaluated at the effective s which
  is modulo the branch length:
      s_eff = s - branch_length * floor(s/branch_length)

  Note: If there are multiple elements that are at the given s position due to the presence of
  an element with a negative length, which of the possible elements is actually chosen is ill-defined.

  Parameters
  ----------
  lat : LatStruct
      Lattice of elements.
  branch : BranchStruct
      Branch to use
  s : float
      Longitudinal position.
  choose_max : bool
      See above
  ix_branch : int, optional
      Branch index. Default is 0.
  print_err : bool, optional
      Print error message if there is an error? Default is True.

  Returns
  -------
  ix_ele : int
      Index of element at s.
  err_flag : bool
      Set True if s is out of bounds. False otherwise.
  s_eff : float
      Effective s. Equal to s with a open lattice. See above.
  position : CoordStruct
      Positional information. .s         -- Same as input s. .ix_ele    -- Same as output ix_ele .location  --
      Location relative to element. Upstream_end$, downstream_end$, or inside$

  Notes
  -----
  Related routines:
  pointer_to_element_at_s ix_ele = ix_end_of_branch branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
  ix_ele = 0 branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s choose_max = True => ix_ele = ix2 choose_max
  = False => ix_ele = ix1 The setting of choose_max only makes a difference when s corresponds to an element
  boundary. For a circular lattice s is evaluated at the effective s which s_eff = s - branch_length *
  floor(s/branch_length) If there are multiple elements that are at the given s position due to the presence of
  an element with a negative length which of the possible elements is actually chosen is ill-defined.
  Overloaded versions:
  )""");
  py::class_<Bmad::ElementAtSLat, std::unique_ptr<Bmad::ElementAtSLat>>(
      m, "ElementAtSLat", "element_at_s_lat return type")
      .def_readonly("err_flag", &Bmad::ElementAtSLat::err_flag)
      .def_readonly("s_eff", &Bmad::ElementAtSLat::s_eff)
      .def_readonly("position", &Bmad::ElementAtSLat::position)
      .def_readonly("ix_ele", &Bmad::ElementAtSLat::ix_ele)
      .def("__len__", [](const Bmad::ElementAtSLat&) { return 4; })
      .def(
          "__getitem__", [](const Bmad::ElementAtSLat& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.s_eff);
            if (i == 2)
              return py::cast(s.position);
            if (i == 3)
              return py::cast(s.ix_ele);
            throw py::index_error();
          });
  m.def(
      "element_at_s",
      py::overload_cast<
          LatProxy&,
          double,
          bool,
          std::optional<int>,
          std::optional<bool>>(&Bmad::element_at_s),
      py::arg("lat"),
      py::arg("s"),
      py::arg("choose_max"),
      py::arg("ix_branch") = py::none(),
      py::arg("print_err") = py::none(),
      R"""(Function element_at_s (...) result (ix_ele)

  Function to return the index of the element at position s.

  element_at_s is an overloaded name for:
    function element_at_s_lat (lat, s, choose_max, ix_branch, err_flag, s_eff, position, print_err) result (ix_ele)
    function element_at_s_branch (branch, s, choose_max, err_flag, s_eff, position, print_err) result (ix_ele)

  The differnce between these two routine is that with element_at_s_lat, the branch is given by the lat
    and ix_ele arguments: branch = lat%branch(ix_ele). With element_at_s_branch, the branch is an argument.

  Also see: pointer_to_element_at_s

  ix_ele is choisen such that:
  If choose_max = True:
      If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
      Else: branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
  If choose_max = False:
      If s = branch%ele(0)%s: ix_ele = 0
      Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
  That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
      choose_max = True  => ix_ele = ix2
      choose_max = False => ix_ele = ix1

  The setting of choose_max only makes a difference when s corresponds to an element boundary.

  Note: For a circular lattice, s is evaluated at the effective s which
  is modulo the branch length:
      s_eff = s - branch_length * floor(s/branch_length)

  Note: If there are multiple elements that are at the given s position due to the presence of
  an element with a negative length, which of the possible elements is actually chosen is ill-defined.

  Parameters
  ----------
  lat : LatStruct
      Lattice of elements.
  branch : BranchStruct
      Branch to use
  s : float
      Longitudinal position.
  choose_max : bool
      See above
  ix_branch : int, optional
      Branch index. Default is 0.
  print_err : bool, optional
      Print error message if there is an error? Default is True.

  Returns
  -------
  ix_ele : int
      Index of element at s.
  err_flag : bool
      Set True if s is out of bounds. False otherwise.
  s_eff : float
      Effective s. Equal to s with a open lattice. See above.
  position : CoordStruct
      Positional information. .s         -- Same as input s. .ix_ele    -- Same as output ix_ele .location  --
      Location relative to element. Upstream_end$, downstream_end$, or inside$

  Notes
  -----
  Related routines:
  pointer_to_element_at_s ix_ele = ix_end_of_branch branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
  ix_ele = 0 branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s choose_max = True => ix_ele = ix2 choose_max
  = False => ix_ele = ix1 The setting of choose_max only makes a difference when s corresponds to an element
  boundary. For a circular lattice s is evaluated at the effective s which s_eff = s - branch_length *
  floor(s/branch_length) If there are multiple elements that are at the given s position due to the presence of
  an element with a negative length which of the possible elements is actually chosen is ill-defined.
  Overloaded versions:
  )""");
  m.def(
      "element_slice_iterator",
      &Bmad::element_slice_iterator,
      py::arg("ele"),
      py::arg("param"),
      py::arg("i_slice"),
      py::arg("n_slice_tot"),
      py::arg("sliced_ele"),
      py::arg("s_start") = py::none(),
      py::arg("s_end") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to slice and dice.
  param : LatParamStruct
      Lattice parameters
  i_slice : int
      Slice index
  n_slice_tot : int
      Total number of slices.
  sliced_ele : 
  s_start : float, optional
      Starting edge of slice relative to beginning of element.
  s_end : float, optional
      Ending edge of slice relative to beginning of element.
  )""");
  m.def("ellipinc_test", &Bmad::ellipinc_test, R"""()""");
  py::class_<Bmad::EmFieldCalc, std::unique_ptr<Bmad::EmFieldCalc>>(
      m, "EmFieldCalc", "em_field_calc return type")
      .def_readonly("field", &Bmad::EmFieldCalc::field)
      .def_readonly("err_flag", &Bmad::EmFieldCalc::err_flag)
      .def("__len__", [](const Bmad::EmFieldCalc&) { return 2; })
      .def("__getitem__", [](const Bmad::EmFieldCalc& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.field);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "em_field_calc",
      &Bmad::em_field_calc,
      py::arg("ele"),
      py::arg("param"),
      py::arg("s_pos"),
      py::arg("orbit"),
      py::arg("local_ref_frame"),
      py::arg("calc_dfield") = py::none(),
      py::arg("calc_potential") = py::none(),
      py::arg("use_overlap") = py::none(),
      py::arg("grid_allow_s_out_of_bounds") = py::none(),
      py::arg("rf_time") = py::none(),
      py::arg("used_eles") = py::none(),
      py::arg("print_err") = py::none(),
      py::arg("original_ele") = py::none(),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  param : LatParamStruct
      Lattice parameters.
  s_pos : float
      Longitudinal position. If local_ref_frame = T: In Body coords relative to the entrance edge of the
      element. If local_ref_frame = F: In Lab coords relative to the upstream edge of the element.
  orbit : CoordStruct
      Transverse coordinates. .vec(1), .vec(3)    -- Transverse coords. .t                  -- Used with
      absolute time tracking. .vec(5)             -- Used with relative time tracking (except with time Runge-
      Kutta).
  local_ref_frame : 
      Logical, If True then take the input coordinates and output fields as being with respect to the frame of
      referene of the element (ignore misalignments).
  field : EmFieldStruct
      E and B fields and derivatives.
  calc_dfield : bool, optional
      If present and True then calculate the field derivatives.
  err_flag : bool
      Set True if there is an error. False otherwise.
  calc_potential : bool, optional
      Calc electric and magnetic potentials? Default is false. This is experimental and only implemented for
      wigglers at present.
  use_overlap : bool, optional
      Add in overlap fields from other elements? Default is True.
  grid_allow_s_out_of_bounds : bool, optional
      For grids, allow s-coordinate to be grossly out of bounds -- logical, optional: For grids, allow
      s-coordinate to be grossly out of bounds and return zero instead of an error? Default: False. Used
      internally for overlapping fields.
  rf_time : float, optional
      Set the time relative to the RF clock. Normally this time is calculated using orbit.t or orbit.vec(5) but
      sometimes it is convenient to be able to override this. For example, time_runge_kutta uses this.
  used_eles : ElePointerStruct, optional
      For internal use only when this routine is called recursively. Used to prevent double counting when there
      is field overlap.
  print_err : bool, optional
      Print an error message? Default is True. For example, if the particle is out of bounds when the field is
      defined on a grid.
  original_ele : EleStruct, optional
      Used with recursive calls that pass the lord as the ele argument. In this case original_ele is the
      original ele argument.
  )""");
  py::class_<PyEmFieldDerivatives, std::unique_ptr<PyEmFieldDerivatives>>(
      m, "EmFieldDerivatives", "em_field_derivatives return type")
      .def_readonly("dfield", &PyEmFieldDerivatives::dfield)
      .def_readonly("s_pos", &PyEmFieldDerivatives::s_pos)
      .def_readonly("local_ref_frame", &PyEmFieldDerivatives::local_ref_frame)
      .def_readonly(
          "grid_allow_s_out_of_bounds",
          &PyEmFieldDerivatives::grid_allow_s_out_of_bounds)
      .def_readonly("rf_time", &PyEmFieldDerivatives::rf_time)
      .def("__len__", [](const PyEmFieldDerivatives&) { return 5; })
      .def(
          "__getitem__",
          [](const PyEmFieldDerivatives& s, int i) -> py::object {
            if (i < 0)
              i += 5;
            if (i == 0)
              return py::cast(s.dfield);
            if (i == 1)
              return py::cast(s.s_pos);
            if (i == 2)
              return py::cast(s.local_ref_frame);
            if (i == 3)
              return py::cast(s.grid_allow_s_out_of_bounds);
            if (i == 4)
              return py::cast(s.rf_time);
            throw py::index_error();
          });
  m.def(
      "em_field_derivatives",
      &python_em_field_derivatives,
      py::arg("ele"),
      py::arg("param"),
      py::arg("s_pos"),
      py::arg("orbit"),
      py::arg("local_ref_frame"),
      py::arg("grid_allow_s_out_of_bounds") = py::none(),
      py::arg("rf_time") = py::none(),
      R"""(Subroutine em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, dfield, grid_allow_s_out_of_bounds, rf_time)

  Routine to calculate field derivatives.
  In theory this should be handled by em_filed_calc. In practice, em_field_calc is currently incomplete.

  Input
    ele             -- Ele_struct: Element
    param           -- lat_param_struct: Lattice parameters.
    s_pos           -- Real(rp): Longitudinal position relative to the upstream edge of the element.
    time            -- Real(rp): Particle time.
                        For absolute time tracking this is the absolute time.
                        For relative time tracking this is relative to the reference particle entering the element.
    orbit           -- Coord_struct: Transverse coordinates.
      %vec(1), %vec(3)  -- Transverse coords. These are the only components used in the calculation.
    local_ref_frame     -- Logical, If True then take the input coordinates and output fields
                                    as being with respect to the frame of referene of the element (ignore misalignments).
    grid_allow_s_out_of_bounds
                     -- logical, optional: For grids, allow s-coordinate to be grossly out of bounds
                          and return zero instead of an error? Default: False. Used internally for overlapping fields.
    rf_time          -- real(rp), optional: RF clock time. If not present then the time will be calculated using the standard algorithm.


  Returns
  -------
  dfield : EmFieldStruct
      E and B field derivatives. dfield.E and dfield.B are not touched.
  )""");
  m.def(
      "em_field_kick_vector_time",
      &Bmad::em_field_kick_vector_time,
      py::arg("ele"),
      py::arg("param"),
      py::arg("rf_time"),
      py::arg("orbit"),
      py::arg("err_flag"),
      py::arg("print_err") = py::none(),
      py::arg("extra_field") = py::none(),
      R"""(Subroutine em_field_kick_vector_time (ele, param, rf_time, orbit, dvec_dt, err_flag, print_err, extra_field))

  Subroutine to convert particle coordinates from t-based to s-based system.

  Parameters
  ----------
  ele : CoordStruct
      input particle
  param : float
      Reference momentum. The sign indicates direction of p_s.
  rf_time : float
      RF time.
  orbit : CoordStruct
      in t-based system
  err_flag : bool
      Set True if there is an error. False otherwise.
  print_err : bool, optional
      Passed to em_field_calc
  extra_field : EmFieldStruct, optional
      Static field to be added to the element field. Eg used with space charge.

  Returns
  -------
  dvec_dt : float
      Derivatives.
  )""");
  m.def(
      "em_field_plus_em_field",
      &Bmad::em_field_plus_em_field,
      py::arg("field1"),
      py::arg("field2"),
      py::arg("field_tot"),
      R"""(Parameters
  ----------
  field1 : 
  field2 : 
  field_tot : 
  )""");
  m.def(
      "em_taylor_equal_em_taylor",
      &Bmad::em_taylor_equal_em_taylor,
      py::arg("em_taylor1"),
      py::arg("em_taylor2"),
      R"""(Parameters
  ----------
  em_taylor1 : 
  em_taylor2 : 
  )""");
  m.def(
      "em_taylors_equal_em_taylors",
      &Bmad::em_taylors_equal_em_taylors,
      py::arg("em_taylor1"),
      py::arg("em_taylor2"),
      R"""(Parameters
  ----------
  em_taylor1 : 
  em_taylor2 : 
  )""");
  py::class_<Bmad::Emit6d, std::unique_ptr<Bmad::Emit6d>>(
      m, "Emit6d", "emit_6d return type")
      .def_readonly("mode", &Bmad::Emit6d::mode)
      .def_readonly("sigma_mat", &Bmad::Emit6d::sigma_mat)
      .def_readonly("rad_int_by_ele", &Bmad::Emit6d::rad_int_by_ele)
      .def("__len__", [](const Bmad::Emit6d&) { return 3; })
      .def("__getitem__", [](const Bmad::Emit6d& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.mode);
        if (i == 1)
          return py::cast(s.sigma_mat);
        if (i == 2)
          return py::cast(s.rad_int_by_ele);
        throw py::index_error();
      });
  m.def(
      "emit_6d",
      &Bmad::emit_6d,
      py::arg("ele_ref"),
      py::arg("include_opening_angle"),
      py::arg("closed_orbit") = py::none(),
      R"""(Subroutine emit_6d (ele_ref, include_opening_angle, mode, sigma_mat, closed_orbit, rad_int_by_ele)

  Routine to calculate the three normal mode emittances, damping partition numbers, radiation integrals, etc.
  Since the emattances, etc. are only an invariant in the limit of zero damping, the calculated
  values will vary depending upon the reference element.

  If the lattice geometry is open, only the radiation integrals is computed.

  Parameters
  ----------
  ele_ref : EleStruct
      Origin of the 1-turn maps used to evaluate the emittances.
  include_opening_angle : bool
      If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
      comparing against other codes.
  closed_orbit : CoordStruct, optional
      Closed orbit. If not present this routine will calculate it.

  Returns
  -------
  mode : NormalModesStruct
      Emittance and other info.
  sigma_mat : float
      Sigma matrix.
  rad_int_by_ele : RadIntAllEleStruct
      Radiation integrals element-by-element.
  )""");
  py::class_<PyEnteringElement, std::unique_ptr<PyEnteringElement>>(
      m, "EnteringElement", "entering_element return type")
      .def_readonly("is_entering", &PyEnteringElement::is_entering)
      .def("__len__", [](const PyEnteringElement&) { return 1; })
      .def("__getitem__", [](const PyEnteringElement& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_entering);
        throw py::index_error();
      });
  m.def(
      "entering_element",
      &python_entering_element,
      py::arg("orbit"),
      py::arg("particle_at"),
      py::arg("is_entering"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle orbit.
  particle_at : int
      First_track_edge$ or second_track_edge$
  is_entering : 
  )""");
  m.def(
      "envelope_radints",
      &Bmad::envelope_radints,
      py::arg("Lambda"),
      py::arg("Theta"),
      py::arg("Iota"),
      py::arg("alpha"),
      py::arg("emit"),
      R"""(subroutine envelope_radints(Lambda,Theta,Iota,alpha,emit)

  Calculates damping decrement and emittance of the three
  normal modes from the integrate diffusion, damping, and vertical
  excitation matrices names Lambda, Theta, and Iota, respectively.
  These three matrices are obtained from the subroutine integrated_mats.

  The damping times can obtained from alpha using:
     tau = lattice_length/c_light/alpha

  )""");
  py::class_<
      Bmad::EnvelopeRadintsIbs,
      std::unique_ptr<Bmad::EnvelopeRadintsIbs>>(
      m, "EnvelopeRadintsIbs", "envelope_radints_ibs return type")
      .def_readonly("alpha", &Bmad::EnvelopeRadintsIbs::alpha)
      .def_readonly("emit", &Bmad::EnvelopeRadintsIbs::emit)
      .def("__len__", [](const Bmad::EnvelopeRadintsIbs&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::EnvelopeRadintsIbs& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.alpha);
            if (i == 1)
              return py::cast(s.emit);
            throw py::index_error();
          });
  m.def(
      "envelope_radints_ibs",
      &Bmad::envelope_radints_ibs,
      py::arg("Lambda"),
      py::arg("Theta"),
      py::arg("Iota"),
      py::arg("eles"),
      py::arg("mode"),
      py::arg("tail_cut"),
      py::arg("npart"),
      py::arg("species"),
      R"""(subroutine envelope_radints_ibs(Lambda, Theta, Iota, eles, alpha, emit, mode, tail_cut, npart, species)

  Calculates damping decrement and emittance of the three
  normal modes by integrating the IBS, SR diffusion, and SR damping matrices.

  The IBS depends on the envelope, and so this routine iterates to
  locate the equilibrium beam envelope. This iterative process can fail to converge.

  The damping times can obtained from alpha using:
     tau = lattice_length/c_light/alpha

  alpha and emit are quantities for the three normal modes.
  alpha and emit are ordered by plane dominance.

  Only radiation from sbends and rbends is taken into account.
  The one-turn transfer matrix at each element (slice) is obtained
  by concatenating the individual element transfer matrices.

  Parameters
  ----------
  Lambda : float
      Integrated damping matrix.
  Theta : float
      Integrated diffusion matrix.
  Iota : float
      Integrated vertical excitation matrix.
  eles : EleStruct
      array of element structures representing ring. .mat6(6,6)            -- real(rp): element transfer matrix.
      .value(l$)            -- real(rp): element (slice) length. .value(E_TOT$)        -- real(rp): Beam energy
      in element.
  mode : unknown
      tune of a-mode. .b.tune                  -- real(rp): tune of b-mode. .z.tune                  --
      real(rp): tune of z-mode.
  tail_cut : bool
      apply tail cut.
  npart : float
      number of particles in bunch.
  species : int
      Particle species.

  Returns
  -------
  alpha : float
      Normal mode damping decrements.
  emit : float
      Normal mode emittances.
  )""");
  py::class_<PyEqAcKicker, std::unique_ptr<PyEqAcKicker>>(
      m, "EqAcKicker", "eq_ac_kicker return type")
      .def_readonly("is_eq", &PyEqAcKicker::is_eq)
      .def("__len__", [](const PyEqAcKicker&) { return 1; })
      .def("__getitem__", [](const PyEqAcKicker& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_ac_kicker",
      &python_eq_ac_kicker,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqAcKickerFreq, std::unique_ptr<PyEqAcKickerFreq>>(
      m, "EqAcKickerFreq", "eq_ac_kicker_freq return type")
      .def_readonly("is_eq", &PyEqAcKickerFreq::is_eq)
      .def("__len__", [](const PyEqAcKickerFreq&) { return 1; })
      .def("__getitem__", [](const PyEqAcKickerFreq& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_ac_kicker_freq",
      &python_eq_ac_kicker_freq,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqAcKickerTime, std::unique_ptr<PyEqAcKickerTime>>(
      m, "EqAcKickerTime", "eq_ac_kicker_time return type")
      .def_readonly("is_eq", &PyEqAcKickerTime::is_eq)
      .def("__len__", [](const PyEqAcKickerTime&) { return 1; })
      .def("__getitem__", [](const PyEqAcKickerTime& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_ac_kicker_time",
      &python_eq_ac_kicker_time,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqAnormalMode, std::unique_ptr<PyEqAnormalMode>>(
      m, "EqAnormalMode", "eq_anormal_mode return type")
      .def_readonly("is_eq", &PyEqAnormalMode::is_eq)
      .def("__len__", [](const PyEqAnormalMode&) { return 1; })
      .def("__getitem__", [](const PyEqAnormalMode& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_anormal_mode",
      &python_eq_anormal_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqApertureParam, std::unique_ptr<PyEqApertureParam>>(
      m, "EqApertureParam", "eq_aperture_param return type")
      .def_readonly("is_eq", &PyEqApertureParam::is_eq)
      .def("__len__", [](const PyEqApertureParam&) { return 1; })
      .def("__getitem__", [](const PyEqApertureParam& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_aperture_param",
      &python_eq_aperture_param,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqAperturePoint, std::unique_ptr<PyEqAperturePoint>>(
      m, "EqAperturePoint", "eq_aperture_point return type")
      .def_readonly("is_eq", &PyEqAperturePoint::is_eq)
      .def("__len__", [](const PyEqAperturePoint&) { return 1; })
      .def("__getitem__", [](const PyEqAperturePoint& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_aperture_point",
      &python_eq_aperture_point,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqApertureScan, std::unique_ptr<PyEqApertureScan>>(
      m, "EqApertureScan", "eq_aperture_scan return type")
      .def_readonly("is_eq", &PyEqApertureScan::is_eq)
      .def("__len__", [](const PyEqApertureScan&) { return 1; })
      .def("__getitem__", [](const PyEqApertureScan& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_aperture_scan",
      &python_eq_aperture_scan,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBeam, std::unique_ptr<PyEqBeam>>(
      m, "EqBeam", "eq_beam return type")
      .def_readonly("is_eq", &PyEqBeam::is_eq)
      .def("__len__", [](const PyEqBeam&) { return 1; })
      .def("__getitem__", [](const PyEqBeam& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_beam",
      &python_eq_beam,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBeamInit, std::unique_ptr<PyEqBeamInit>>(
      m, "EqBeamInit", "eq_beam_init return type")
      .def_readonly("is_eq", &PyEqBeamInit::is_eq)
      .def("__len__", [](const PyEqBeamInit&) { return 1; })
      .def("__getitem__", [](const PyEqBeamInit& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_beam_init",
      &python_eq_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBmadCommon, std::unique_ptr<PyEqBmadCommon>>(
      m, "EqBmadCommon", "eq_bmad_common return type")
      .def_readonly("is_eq", &PyEqBmadCommon::is_eq)
      .def("__len__", [](const PyEqBmadCommon&) { return 1; })
      .def("__getitem__", [](const PyEqBmadCommon& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_bmad_common",
      &python_eq_bmad_common,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBookkeepingState, std::unique_ptr<PyEqBookkeepingState>>(
      m, "EqBookkeepingState", "eq_bookkeeping_state return type")
      .def_readonly("is_eq", &PyEqBookkeepingState::is_eq)
      .def("__len__", [](const PyEqBookkeepingState&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqBookkeepingState& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_bookkeeping_state",
      &python_eq_bookkeeping_state,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBpmPhaseCoupling, std::unique_ptr<PyEqBpmPhaseCoupling>>(
      m, "EqBpmPhaseCoupling", "eq_bpm_phase_coupling return type")
      .def_readonly("is_eq", &PyEqBpmPhaseCoupling::is_eq)
      .def("__len__", [](const PyEqBpmPhaseCoupling&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqBpmPhaseCoupling& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_bpm_phase_coupling",
      &python_eq_bpm_phase_coupling,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBranch, std::unique_ptr<PyEqBranch>>(
      m, "EqBranch", "eq_branch return type")
      .def_readonly("is_eq", &PyEqBranch::is_eq)
      .def("__len__", [](const PyEqBranch&) { return 1; })
      .def("__getitem__", [](const PyEqBranch& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_branch",
      &python_eq_branch,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBunch, std::unique_ptr<PyEqBunch>>(
      m, "EqBunch", "eq_bunch return type")
      .def_readonly("is_eq", &PyEqBunch::is_eq)
      .def("__len__", [](const PyEqBunch&) { return 1; })
      .def("__getitem__", [](const PyEqBunch& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_bunch",
      &python_eq_bunch,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqBunchParams, std::unique_ptr<PyEqBunchParams>>(
      m, "EqBunchParams", "eq_bunch_params return type")
      .def_readonly("is_eq", &PyEqBunchParams::is_eq)
      .def("__len__", [](const PyEqBunchParams&) { return 1; })
      .def("__getitem__", [](const PyEqBunchParams& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_bunch_params",
      &python_eq_bunch_params,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCartesianMap, std::unique_ptr<PyEqCartesianMap>>(
      m, "EqCartesianMap", "eq_cartesian_map return type")
      .def_readonly("is_eq", &PyEqCartesianMap::is_eq)
      .def("__len__", [](const PyEqCartesianMap&) { return 1; })
      .def("__getitem__", [](const PyEqCartesianMap& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_cartesian_map",
      &python_eq_cartesian_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCartesianMapTerm, std::unique_ptr<PyEqCartesianMapTerm>>(
      m, "EqCartesianMapTerm", "eq_cartesian_map_term return type")
      .def_readonly("is_eq", &PyEqCartesianMapTerm::is_eq)
      .def("__len__", [](const PyEqCartesianMapTerm&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqCartesianMapTerm& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_cartesian_map_term",
      &python_eq_cartesian_map_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCartesianMapTerm1, std::unique_ptr<PyEqCartesianMapTerm1>>(
      m, "EqCartesianMapTerm1", "eq_cartesian_map_term1 return type")
      .def_readonly("is_eq", &PyEqCartesianMapTerm1::is_eq)
      .def("__len__", [](const PyEqCartesianMapTerm1&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqCartesianMapTerm1& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_cartesian_map_term1",
      &python_eq_cartesian_map_term1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqComplexTaylor, std::unique_ptr<PyEqComplexTaylor>>(
      m, "EqComplexTaylor", "eq_complex_taylor return type")
      .def_readonly("is_eq", &PyEqComplexTaylor::is_eq)
      .def("__len__", [](const PyEqComplexTaylor&) { return 1; })
      .def("__getitem__", [](const PyEqComplexTaylor& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_complex_taylor",
      &python_eq_complex_taylor,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqComplexTaylorTerm, std::unique_ptr<PyEqComplexTaylorTerm>>(
      m, "EqComplexTaylorTerm", "eq_complex_taylor_term return type")
      .def_readonly("is_eq", &PyEqComplexTaylorTerm::is_eq)
      .def("__len__", [](const PyEqComplexTaylorTerm&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqComplexTaylorTerm& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_complex_taylor_term",
      &python_eq_complex_taylor_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqControl, std::unique_ptr<PyEqControl>>(
      m, "EqControl", "eq_control return type")
      .def_readonly("is_eq", &PyEqControl::is_eq)
      .def("__len__", [](const PyEqControl&) { return 1; })
      .def("__getitem__", [](const PyEqControl& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_control",
      &python_eq_control,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqControlRamp1, std::unique_ptr<PyEqControlRamp1>>(
      m, "EqControlRamp1", "eq_control_ramp1 return type")
      .def_readonly("is_eq", &PyEqControlRamp1::is_eq)
      .def("__len__", [](const PyEqControlRamp1&) { return 1; })
      .def("__getitem__", [](const PyEqControlRamp1& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_control_ramp1",
      &python_eq_control_ramp1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqControlVar1, std::unique_ptr<PyEqControlVar1>>(
      m, "EqControlVar1", "eq_control_var1 return type")
      .def_readonly("is_eq", &PyEqControlVar1::is_eq)
      .def("__len__", [](const PyEqControlVar1&) { return 1; })
      .def("__getitem__", [](const PyEqControlVar1& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_control_var1",
      &python_eq_control_var1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqController, std::unique_ptr<PyEqController>>(
      m, "EqController", "eq_controller return type")
      .def_readonly("is_eq", &PyEqController::is_eq)
      .def("__len__", [](const PyEqController&) { return 1; })
      .def("__getitem__", [](const PyEqController& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_controller",
      &python_eq_controller,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCoord, std::unique_ptr<PyEqCoord>>(
      m, "EqCoord", "eq_coord return type")
      .def_readonly("is_eq", &PyEqCoord::is_eq)
      .def("__len__", [](const PyEqCoord&) { return 1; })
      .def("__getitem__", [](const PyEqCoord& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_coord",
      &python_eq_coord,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCoordArray, std::unique_ptr<PyEqCoordArray>>(
      m, "EqCoordArray", "eq_coord_array return type")
      .def_readonly("is_eq", &PyEqCoordArray::is_eq)
      .def("__len__", [](const PyEqCoordArray&) { return 1; })
      .def("__getitem__", [](const PyEqCoordArray& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_coord_array",
      &python_eq_coord_array,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCylindricalMap, std::unique_ptr<PyEqCylindricalMap>>(
      m, "EqCylindricalMap", "eq_cylindrical_map return type")
      .def_readonly("is_eq", &PyEqCylindricalMap::is_eq)
      .def("__len__", [](const PyEqCylindricalMap&) { return 1; })
      .def("__getitem__", [](const PyEqCylindricalMap& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_cylindrical_map",
      &python_eq_cylindrical_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCylindricalMapTerm, std::unique_ptr<PyEqCylindricalMapTerm>>(
      m, "EqCylindricalMapTerm", "eq_cylindrical_map_term return type")
      .def_readonly("is_eq", &PyEqCylindricalMapTerm::is_eq)
      .def("__len__", [](const PyEqCylindricalMapTerm&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqCylindricalMapTerm& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_cylindrical_map_term",
      &python_eq_cylindrical_map_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqCylindricalMapTerm1, std::unique_ptr<PyEqCylindricalMapTerm1>>(
      m, "EqCylindricalMapTerm1", "eq_cylindrical_map_term1 return type")
      .def_readonly("is_eq", &PyEqCylindricalMapTerm1::is_eq)
      .def("__len__", [](const PyEqCylindricalMapTerm1&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqCylindricalMapTerm1& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_cylindrical_map_term1",
      &python_eq_cylindrical_map_term1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqEle, std::unique_ptr<PyEqEle>>(
      m, "EqEle", "eq_ele return type")
      .def_readonly("is_eq", &PyEqEle::is_eq)
      .def("__len__", [](const PyEqEle&) { return 1; })
      .def("__getitem__", [](const PyEqEle& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_ele",
      &python_eq_ele,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqEllipseBeamInit, std::unique_ptr<PyEqEllipseBeamInit>>(
      m, "EqEllipseBeamInit", "eq_ellipse_beam_init return type")
      .def_readonly("is_eq", &PyEqEllipseBeamInit::is_eq)
      .def("__len__", [](const PyEqEllipseBeamInit&) { return 1; })
      .def(
          "__getitem__", [](const PyEqEllipseBeamInit& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_ellipse_beam_init",
      &python_eq_ellipse_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqEmField, std::unique_ptr<PyEqEmField>>(
      m, "EqEmField", "eq_em_field return type")
      .def_readonly("is_eq", &PyEqEmField::is_eq)
      .def("__len__", [](const PyEqEmField&) { return 1; })
      .def("__getitem__", [](const PyEqEmField& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_em_field",
      &python_eq_em_field,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqEmTaylor, std::unique_ptr<PyEqEmTaylor>>(
      m, "EqEmTaylor", "eq_em_taylor return type")
      .def_readonly("is_eq", &PyEqEmTaylor::is_eq)
      .def("__len__", [](const PyEqEmTaylor&) { return 1; })
      .def("__getitem__", [](const PyEqEmTaylor& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_em_taylor",
      &python_eq_em_taylor,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqEmTaylorTerm, std::unique_ptr<PyEqEmTaylorTerm>>(
      m, "EqEmTaylorTerm", "eq_em_taylor_term return type")
      .def_readonly("is_eq", &PyEqEmTaylorTerm::is_eq)
      .def("__len__", [](const PyEqEmTaylorTerm&) { return 1; })
      .def("__getitem__", [](const PyEqEmTaylorTerm& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_em_taylor_term",
      &python_eq_em_taylor_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqExpressionAtom, std::unique_ptr<PyEqExpressionAtom>>(
      m, "EqExpressionAtom", "eq_expression_atom return type")
      .def_readonly("is_eq", &PyEqExpressionAtom::is_eq)
      .def("__len__", [](const PyEqExpressionAtom&) { return 1; })
      .def("__getitem__", [](const PyEqExpressionAtom& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_expression_atom",
      &python_eq_expression_atom,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqFloorPosition, std::unique_ptr<PyEqFloorPosition>>(
      m, "EqFloorPosition", "eq_floor_position return type")
      .def_readonly("is_eq", &PyEqFloorPosition::is_eq)
      .def("__len__", [](const PyEqFloorPosition&) { return 1; })
      .def("__getitem__", [](const PyEqFloorPosition& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_floor_position",
      &python_eq_floor_position,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqGenGrad1, std::unique_ptr<PyEqGenGrad1>>(
      m, "EqGenGrad1", "eq_gen_grad1 return type")
      .def_readonly("is_eq", &PyEqGenGrad1::is_eq)
      .def("__len__", [](const PyEqGenGrad1&) { return 1; })
      .def("__getitem__", [](const PyEqGenGrad1& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_gen_grad1",
      &python_eq_gen_grad1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqGenGradMap, std::unique_ptr<PyEqGenGradMap>>(
      m, "EqGenGradMap", "eq_gen_grad_map return type")
      .def_readonly("is_eq", &PyEqGenGradMap::is_eq)
      .def("__len__", [](const PyEqGenGradMap&) { return 1; })
      .def("__getitem__", [](const PyEqGenGradMap& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_gen_grad_map",
      &python_eq_gen_grad_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqGridBeamInit, std::unique_ptr<PyEqGridBeamInit>>(
      m, "EqGridBeamInit", "eq_grid_beam_init return type")
      .def_readonly("is_eq", &PyEqGridBeamInit::is_eq)
      .def("__len__", [](const PyEqGridBeamInit&) { return 1; })
      .def("__getitem__", [](const PyEqGridBeamInit& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_grid_beam_init",
      &python_eq_grid_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqGridField, std::unique_ptr<PyEqGridField>>(
      m, "EqGridField", "eq_grid_field return type")
      .def_readonly("is_eq", &PyEqGridField::is_eq)
      .def("__len__", [](const PyEqGridField&) { return 1; })
      .def("__getitem__", [](const PyEqGridField& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_grid_field",
      &python_eq_grid_field,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqGridFieldPt, std::unique_ptr<PyEqGridFieldPt>>(
      m, "EqGridFieldPt", "eq_grid_field_pt return type")
      .def_readonly("is_eq", &PyEqGridFieldPt::is_eq)
      .def("__len__", [](const PyEqGridFieldPt&) { return 1; })
      .def("__getitem__", [](const PyEqGridFieldPt& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_grid_field_pt",
      &python_eq_grid_field_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqGridFieldPt1, std::unique_ptr<PyEqGridFieldPt1>>(
      m, "EqGridFieldPt1", "eq_grid_field_pt1 return type")
      .def_readonly("is_eq", &PyEqGridFieldPt1::is_eq)
      .def("__len__", [](const PyEqGridFieldPt1&) { return 1; })
      .def("__getitem__", [](const PyEqGridFieldPt1& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_grid_field_pt1",
      &python_eq_grid_field_pt1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<
      PyEqHighEnergySpaceCharge,
      std::unique_ptr<PyEqHighEnergySpaceCharge>>(
      m, "EqHighEnergySpaceCharge", "eq_high_energy_space_charge return type")
      .def_readonly("is_eq", &PyEqHighEnergySpaceCharge::is_eq)
      .def("__len__", [](const PyEqHighEnergySpaceCharge&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqHighEnergySpaceCharge& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_high_energy_space_charge",
      &python_eq_high_energy_space_charge,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqInterval1Coef, std::unique_ptr<PyEqInterval1Coef>>(
      m, "EqInterval1Coef", "eq_interval1_coef return type")
      .def_readonly("is_eq", &PyEqInterval1Coef::is_eq)
      .def("__len__", [](const PyEqInterval1Coef&) { return 1; })
      .def("__getitem__", [](const PyEqInterval1Coef& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_interval1_coef",
      &python_eq_interval1_coef,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqKvBeamInit, std::unique_ptr<PyEqKvBeamInit>>(
      m, "EqKvBeamInit", "eq_kv_beam_init return type")
      .def_readonly("is_eq", &PyEqKvBeamInit::is_eq)
      .def("__len__", [](const PyEqKvBeamInit&) { return 1; })
      .def("__getitem__", [](const PyEqKvBeamInit& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_kv_beam_init",
      &python_eq_kv_beam_init,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqLat, std::unique_ptr<PyEqLat>>(
      m, "EqLat", "eq_lat return type")
      .def_readonly("is_eq", &PyEqLat::is_eq)
      .def("__len__", [](const PyEqLat&) { return 1; })
      .def("__getitem__", [](const PyEqLat& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_lat",
      &python_eq_lat,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqLatEleLoc, std::unique_ptr<PyEqLatEleLoc>>(
      m, "EqLatEleLoc", "eq_lat_ele_loc return type")
      .def_readonly("is_eq", &PyEqLatEleLoc::is_eq)
      .def("__len__", [](const PyEqLatEleLoc&) { return 1; })
      .def("__getitem__", [](const PyEqLatEleLoc& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_lat_ele_loc",
      &python_eq_lat_ele_loc,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqLatParam, std::unique_ptr<PyEqLatParam>>(
      m, "EqLatParam", "eq_lat_param return type")
      .def_readonly("is_eq", &PyEqLatParam::is_eq)
      .def("__len__", [](const PyEqLatParam&) { return 1; })
      .def("__getitem__", [](const PyEqLatParam& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_lat_param",
      &python_eq_lat_param,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqLinacNormalMode, std::unique_ptr<PyEqLinacNormalMode>>(
      m, "EqLinacNormalMode", "eq_linac_normal_mode return type")
      .def_readonly("is_eq", &PyEqLinacNormalMode::is_eq)
      .def("__len__", [](const PyEqLinacNormalMode&) { return 1; })
      .def(
          "__getitem__", [](const PyEqLinacNormalMode& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_linac_normal_mode",
      &python_eq_linac_normal_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqMode3, std::unique_ptr<PyEqMode3>>(
      m, "EqMode3", "eq_mode3 return type")
      .def_readonly("is_eq", &PyEqMode3::is_eq)
      .def("__len__", [](const PyEqMode3&) { return 1; })
      .def("__getitem__", [](const PyEqMode3& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_mode3",
      &python_eq_mode3,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqModeInfo, std::unique_ptr<PyEqModeInfo>>(
      m, "EqModeInfo", "eq_mode_info return type")
      .def_readonly("is_eq", &PyEqModeInfo::is_eq)
      .def("__len__", [](const PyEqModeInfo&) { return 1; })
      .def("__getitem__", [](const PyEqModeInfo& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_mode_info",
      &python_eq_mode_info,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqNormalModes, std::unique_ptr<PyEqNormalModes>>(
      m, "EqNormalModes", "eq_normal_modes return type")
      .def_readonly("is_eq", &PyEqNormalModes::is_eq)
      .def("__len__", [](const PyEqNormalModes&) { return 1; })
      .def("__getitem__", [](const PyEqNormalModes& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_normal_modes",
      &python_eq_normal_modes,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPhotonElement, std::unique_ptr<PyEqPhotonElement>>(
      m, "EqPhotonElement", "eq_photon_element return type")
      .def_readonly("is_eq", &PyEqPhotonElement::is_eq)
      .def("__len__", [](const PyEqPhotonElement&) { return 1; })
      .def("__getitem__", [](const PyEqPhotonElement& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_photon_element",
      &python_eq_photon_element,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPhotonMaterial, std::unique_ptr<PyEqPhotonMaterial>>(
      m, "EqPhotonMaterial", "eq_photon_material return type")
      .def_readonly("is_eq", &PyEqPhotonMaterial::is_eq)
      .def("__len__", [](const PyEqPhotonMaterial&) { return 1; })
      .def("__getitem__", [](const PyEqPhotonMaterial& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_photon_material",
      &python_eq_photon_material,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<
      PyEqPhotonReflectSurface,
      std::unique_ptr<PyEqPhotonReflectSurface>>(
      m, "EqPhotonReflectSurface", "eq_photon_reflect_surface return type")
      .def_readonly("is_eq", &PyEqPhotonReflectSurface::is_eq)
      .def("__len__", [](const PyEqPhotonReflectSurface&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqPhotonReflectSurface& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_photon_reflect_surface",
      &python_eq_photon_reflect_surface,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPhotonReflectTable, std::unique_ptr<PyEqPhotonReflectTable>>(
      m, "EqPhotonReflectTable", "eq_photon_reflect_table return type")
      .def_readonly("is_eq", &PyEqPhotonReflectTable::is_eq)
      .def("__len__", [](const PyEqPhotonReflectTable&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqPhotonReflectTable& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_photon_reflect_table",
      &python_eq_photon_reflect_table,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPhotonTarget, std::unique_ptr<PyEqPhotonTarget>>(
      m, "EqPhotonTarget", "eq_photon_target return type")
      .def_readonly("is_eq", &PyEqPhotonTarget::is_eq)
      .def("__len__", [](const PyEqPhotonTarget&) { return 1; })
      .def("__getitem__", [](const PyEqPhotonTarget& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_photon_target",
      &python_eq_photon_target,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPixelDetec, std::unique_ptr<PyEqPixelDetec>>(
      m, "EqPixelDetec", "eq_pixel_detec return type")
      .def_readonly("is_eq", &PyEqPixelDetec::is_eq)
      .def("__len__", [](const PyEqPixelDetec&) { return 1; })
      .def("__getitem__", [](const PyEqPixelDetec& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_pixel_detec",
      &python_eq_pixel_detec,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPixelPt, std::unique_ptr<PyEqPixelPt>>(
      m, "EqPixelPt", "eq_pixel_pt return type")
      .def_readonly("is_eq", &PyEqPixelPt::is_eq)
      .def("__len__", [](const PyEqPixelPt&) { return 1; })
      .def("__getitem__", [](const PyEqPixelPt& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_pixel_pt",
      &python_eq_pixel_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqPreTracker, std::unique_ptr<PyEqPreTracker>>(
      m, "EqPreTracker", "eq_pre_tracker return type")
      .def_readonly("is_eq", &PyEqPreTracker::is_eq)
      .def("__len__", [](const PyEqPreTracker&) { return 1; })
      .def("__getitem__", [](const PyEqPreTracker& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_pre_tracker",
      &python_eq_pre_tracker,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqRadInt1, std::unique_ptr<PyEqRadInt1>>(
      m, "EqRadInt1", "eq_rad_int1 return type")
      .def_readonly("is_eq", &PyEqRadInt1::is_eq)
      .def("__len__", [](const PyEqRadInt1&) { return 1; })
      .def("__getitem__", [](const PyEqRadInt1& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_rad_int1",
      &python_eq_rad_int1,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqRadIntAllEle, std::unique_ptr<PyEqRadIntAllEle>>(
      m, "EqRadIntAllEle", "eq_rad_int_all_ele return type")
      .def_readonly("is_eq", &PyEqRadIntAllEle::is_eq)
      .def("__len__", [](const PyEqRadIntAllEle&) { return 1; })
      .def("__getitem__", [](const PyEqRadIntAllEle& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_rad_int_all_ele",
      &python_eq_rad_int_all_ele,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqRadIntBranch, std::unique_ptr<PyEqRadIntBranch>>(
      m, "EqRadIntBranch", "eq_rad_int_branch return type")
      .def_readonly("is_eq", &PyEqRadIntBranch::is_eq)
      .def("__len__", [](const PyEqRadIntBranch&) { return 1; })
      .def("__getitem__", [](const PyEqRadIntBranch& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_rad_int_branch",
      &python_eq_rad_int_branch,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqRadMap, std::unique_ptr<PyEqRadMap>>(
      m, "EqRadMap", "eq_rad_map return type")
      .def_readonly("is_eq", &PyEqRadMap::is_eq)
      .def("__len__", [](const PyEqRadMap&) { return 1; })
      .def("__getitem__", [](const PyEqRadMap& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_rad_map",
      &python_eq_rad_map,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqRadMapEle, std::unique_ptr<PyEqRadMapEle>>(
      m, "EqRadMapEle", "eq_rad_map_ele return type")
      .def_readonly("is_eq", &PyEqRadMapEle::is_eq)
      .def("__len__", [](const PyEqRadMapEle&) { return 1; })
      .def("__getitem__", [](const PyEqRadMapEle& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_rad_map_ele",
      &python_eq_rad_map_ele,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqRamperLord, std::unique_ptr<PyEqRamperLord>>(
      m, "EqRamperLord", "eq_ramper_lord return type")
      .def_readonly("is_eq", &PyEqRamperLord::is_eq)
      .def("__len__", [](const PyEqRamperLord&) { return 1; })
      .def("__getitem__", [](const PyEqRamperLord& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_ramper_lord",
      &python_eq_ramper_lord,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSpaceChargeCommon, std::unique_ptr<PyEqSpaceChargeCommon>>(
      m, "EqSpaceChargeCommon", "eq_space_charge_common return type")
      .def_readonly("is_eq", &PyEqSpaceChargeCommon::is_eq)
      .def("__len__", [](const PyEqSpaceChargeCommon&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSpaceChargeCommon& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_space_charge_common",
      &python_eq_space_charge_common,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSpinPolar, std::unique_ptr<PyEqSpinPolar>>(
      m, "EqSpinPolar", "eq_spin_polar return type")
      .def_readonly("is_eq", &PyEqSpinPolar::is_eq)
      .def("__len__", [](const PyEqSpinPolar&) { return 1; })
      .def("__getitem__", [](const PyEqSpinPolar& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_spin_polar",
      &python_eq_spin_polar,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSpline, std::unique_ptr<PyEqSpline>>(
      m, "EqSpline", "eq_spline return type")
      .def_readonly("is_eq", &PyEqSpline::is_eq)
      .def("__len__", [](const PyEqSpline&) { return 1; })
      .def("__getitem__", [](const PyEqSpline& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_spline",
      &python_eq_spline,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqStrongBeam, std::unique_ptr<PyEqStrongBeam>>(
      m, "EqStrongBeam", "eq_strong_beam return type")
      .def_readonly("is_eq", &PyEqStrongBeam::is_eq)
      .def("__len__", [](const PyEqStrongBeam&) { return 1; })
      .def("__getitem__", [](const PyEqStrongBeam& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_strong_beam",
      &python_eq_strong_beam,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSurfaceCurvature, std::unique_ptr<PyEqSurfaceCurvature>>(
      m, "EqSurfaceCurvature", "eq_surface_curvature return type")
      .def_readonly("is_eq", &PyEqSurfaceCurvature::is_eq)
      .def("__len__", [](const PyEqSurfaceCurvature&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceCurvature& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_curvature",
      &python_eq_surface_curvature,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSurfaceDisplacement, std::unique_ptr<PyEqSurfaceDisplacement>>(
      m, "EqSurfaceDisplacement", "eq_surface_displacement return type")
      .def_readonly("is_eq", &PyEqSurfaceDisplacement::is_eq)
      .def("__len__", [](const PyEqSurfaceDisplacement&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceDisplacement& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_displacement",
      &python_eq_surface_displacement,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<
      PyEqSurfaceDisplacementPt,
      std::unique_ptr<PyEqSurfaceDisplacementPt>>(
      m, "EqSurfaceDisplacementPt", "eq_surface_displacement_pt return type")
      .def_readonly("is_eq", &PyEqSurfaceDisplacementPt::is_eq)
      .def("__len__", [](const PyEqSurfaceDisplacementPt&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceDisplacementPt& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_displacement_pt",
      &python_eq_surface_displacement_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSurfaceHMisalign, std::unique_ptr<PyEqSurfaceHMisalign>>(
      m, "EqSurfaceHMisalign", "eq_surface_h_misalign return type")
      .def_readonly("is_eq", &PyEqSurfaceHMisalign::is_eq)
      .def("__len__", [](const PyEqSurfaceHMisalign&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceHMisalign& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_h_misalign",
      &python_eq_surface_h_misalign,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSurfaceHMisalignPt, std::unique_ptr<PyEqSurfaceHMisalignPt>>(
      m, "EqSurfaceHMisalignPt", "eq_surface_h_misalign_pt return type")
      .def_readonly("is_eq", &PyEqSurfaceHMisalignPt::is_eq)
      .def("__len__", [](const PyEqSurfaceHMisalignPt&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceHMisalignPt& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_h_misalign_pt",
      &python_eq_surface_h_misalign_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSurfaceSegmented, std::unique_ptr<PyEqSurfaceSegmented>>(
      m, "EqSurfaceSegmented", "eq_surface_segmented return type")
      .def_readonly("is_eq", &PyEqSurfaceSegmented::is_eq)
      .def("__len__", [](const PyEqSurfaceSegmented&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceSegmented& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_segmented",
      &python_eq_surface_segmented,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqSurfaceSegmentedPt, std::unique_ptr<PyEqSurfaceSegmentedPt>>(
      m, "EqSurfaceSegmentedPt", "eq_surface_segmented_pt return type")
      .def_readonly("is_eq", &PyEqSurfaceSegmentedPt::is_eq)
      .def("__len__", [](const PyEqSurfaceSegmentedPt&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEqSurfaceSegmentedPt& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_eq);
            throw py::index_error();
          });
  m.def(
      "eq_surface_segmented_pt",
      &python_eq_surface_segmented_pt,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqTargetPoint, std::unique_ptr<PyEqTargetPoint>>(
      m, "EqTargetPoint", "eq_target_point return type")
      .def_readonly("is_eq", &PyEqTargetPoint::is_eq)
      .def("__len__", [](const PyEqTargetPoint&) { return 1; })
      .def("__getitem__", [](const PyEqTargetPoint& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_target_point",
      &python_eq_target_point,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqTaylor, std::unique_ptr<PyEqTaylor>>(
      m, "EqTaylor", "eq_taylor return type")
      .def_readonly("is_eq", &PyEqTaylor::is_eq)
      .def("__len__", [](const PyEqTaylor&) { return 1; })
      .def("__getitem__", [](const PyEqTaylor& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_taylor",
      &python_eq_taylor,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqTaylorTerm, std::unique_ptr<PyEqTaylorTerm>>(
      m, "EqTaylorTerm", "eq_taylor_term return type")
      .def_readonly("is_eq", &PyEqTaylorTerm::is_eq)
      .def("__len__", [](const PyEqTaylorTerm&) { return 1; })
      .def("__getitem__", [](const PyEqTaylorTerm& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_taylor_term",
      &python_eq_taylor_term,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqTrack, std::unique_ptr<PyEqTrack>>(
      m, "EqTrack", "eq_track return type")
      .def_readonly("is_eq", &PyEqTrack::is_eq)
      .def("__len__", [](const PyEqTrack&) { return 1; })
      .def("__getitem__", [](const PyEqTrack& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_track",
      &python_eq_track,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqTrackPoint, std::unique_ptr<PyEqTrackPoint>>(
      m, "EqTrackPoint", "eq_track_point return type")
      .def_readonly("is_eq", &PyEqTrackPoint::is_eq)
      .def("__len__", [](const PyEqTrackPoint&) { return 1; })
      .def("__getitem__", [](const PyEqTrackPoint& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_track_point",
      &python_eq_track_point,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqTwiss, std::unique_ptr<PyEqTwiss>>(
      m, "EqTwiss", "eq_twiss return type")
      .def_readonly("is_eq", &PyEqTwiss::is_eq)
      .def("__len__", [](const PyEqTwiss&) { return 1; })
      .def("__getitem__", [](const PyEqTwiss& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_twiss",
      &python_eq_twiss,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWake, std::unique_ptr<PyEqWake>>(
      m, "EqWake", "eq_wake return type")
      .def_readonly("is_eq", &PyEqWake::is_eq)
      .def("__len__", [](const PyEqWake&) { return 1; })
      .def("__getitem__", [](const PyEqWake& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wake",
      &python_eq_wake,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWakeLr, std::unique_ptr<PyEqWakeLr>>(
      m, "EqWakeLr", "eq_wake_lr return type")
      .def_readonly("is_eq", &PyEqWakeLr::is_eq)
      .def("__len__", [](const PyEqWakeLr&) { return 1; })
      .def("__getitem__", [](const PyEqWakeLr& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wake_lr",
      &python_eq_wake_lr,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWakeLrMode, std::unique_ptr<PyEqWakeLrMode>>(
      m, "EqWakeLrMode", "eq_wake_lr_mode return type")
      .def_readonly("is_eq", &PyEqWakeLrMode::is_eq)
      .def("__len__", [](const PyEqWakeLrMode&) { return 1; })
      .def("__getitem__", [](const PyEqWakeLrMode& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wake_lr_mode",
      &python_eq_wake_lr_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWakeSr, std::unique_ptr<PyEqWakeSr>>(
      m, "EqWakeSr", "eq_wake_sr return type")
      .def_readonly("is_eq", &PyEqWakeSr::is_eq)
      .def("__len__", [](const PyEqWakeSr&) { return 1; })
      .def("__getitem__", [](const PyEqWakeSr& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wake_sr",
      &python_eq_wake_sr,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWakeSrMode, std::unique_ptr<PyEqWakeSrMode>>(
      m, "EqWakeSrMode", "eq_wake_sr_mode return type")
      .def_readonly("is_eq", &PyEqWakeSrMode::is_eq)
      .def("__len__", [](const PyEqWakeSrMode&) { return 1; })
      .def("__getitem__", [](const PyEqWakeSrMode& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wake_sr_mode",
      &python_eq_wake_sr_mode,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWakeSrZLong, std::unique_ptr<PyEqWakeSrZLong>>(
      m, "EqWakeSrZLong", "eq_wake_sr_z_long return type")
      .def_readonly("is_eq", &PyEqWakeSrZLong::is_eq)
      .def("__len__", [](const PyEqWakeSrZLong&) { return 1; })
      .def("__getitem__", [](const PyEqWakeSrZLong& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wake_sr_z_long",
      &python_eq_wake_sr_z_long,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWall3d, std::unique_ptr<PyEqWall3d>>(
      m, "EqWall3d", "eq_wall3d return type")
      .def_readonly("is_eq", &PyEqWall3d::is_eq)
      .def("__len__", [](const PyEqWall3d&) { return 1; })
      .def("__getitem__", [](const PyEqWall3d& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wall3d",
      &python_eq_wall3d,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWall3dSection, std::unique_ptr<PyEqWall3dSection>>(
      m, "EqWall3dSection", "eq_wall3d_section return type")
      .def_readonly("is_eq", &PyEqWall3dSection::is_eq)
      .def("__len__", [](const PyEqWall3dSection&) { return 1; })
      .def("__getitem__", [](const PyEqWall3dSection& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wall3d_section",
      &python_eq_wall3d_section,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqWall3dVertex, std::unique_ptr<PyEqWall3dVertex>>(
      m, "EqWall3dVertex", "eq_wall3d_vertex return type")
      .def_readonly("is_eq", &PyEqWall3dVertex::is_eq)
      .def("__len__", [](const PyEqWall3dVertex&) { return 1; })
      .def("__getitem__", [](const PyEqWall3dVertex& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_wall3d_vertex",
      &python_eq_wall3d_vertex,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqXyDisp, std::unique_ptr<PyEqXyDisp>>(
      m, "EqXyDisp", "eq_xy_disp return type")
      .def_readonly("is_eq", &PyEqXyDisp::is_eq)
      .def("__len__", [](const PyEqXyDisp&) { return 1; })
      .def("__getitem__", [](const PyEqXyDisp& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.is_eq);
        throw py::index_error();
      });
  m.def(
      "eq_xy_disp",
      &python_eq_xy_disp,
      py::arg("f1"),
      py::arg("f2"),
      py::arg("is_eq"),
      R"""(Parameters
  ----------
  f1 : 
  f2 : 
  is_eq : 
  )""");
  py::class_<PyEqualSignHere, std::unique_ptr<PyEqualSignHere>>(
      m, "EqualSignHere", "equal_sign_here return type")
      .def_readonly("delim", &PyEqualSignHere::delim)
      .def_readonly("is_here", &PyEqualSignHere::is_here)
      .def("__len__", [](const PyEqualSignHere&) { return 2; })
      .def("__getitem__", [](const PyEqualSignHere& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.is_here);
        throw py::index_error();
      });
  m.def(
      "equal_sign_here",
      &python_equal_sign_here,
      py::arg("ele"),
      py::arg("delim"),
      py::arg("is_here"),
      R"""(Parameters
  ----------
  ele : 
  delim : 
  is_here : 
  )""");
  py::class_<
      PyEquivalentTaylorAttributes,
      std::unique_ptr<PyEquivalentTaylorAttributes>>(
      m,
      "EquivalentTaylorAttributes",
      "equivalent_taylor_attributes return type")
      .def_readonly("equiv", &PyEquivalentTaylorAttributes::equiv)
      .def("__len__", [](const PyEquivalentTaylorAttributes&) { return 1; })
      .def(
          "__getitem__",
          [](const PyEquivalentTaylorAttributes& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.equiv);
            throw py::index_error();
          });
  m.def(
      "equivalent_taylor_attributes",
      &python_equivalent_taylor_attributes,
      py::arg("ele_taylor"),
      py::arg("ele2"),
      py::arg("equiv"),
      R"""(Parameters
  ----------
  ele_taylor : EleStruct
      Element with a Taylor map
  ele2 : EleStruct
      Element that might receive the Taylor map from ele_taylor.
  equiv : 
  )""");
  py::class_<PyEtdiv, std::unique_ptr<PyEtdiv>>(m, "Etdiv", "etdiv return type")
      .def_readonly("A", &PyEtdiv::A)
      .def_readonly("B", &PyEtdiv::B)
      .def_readonly("C", &PyEtdiv::C)
      .def_readonly("D", &PyEtdiv::D)
      .def_readonly("E", &PyEtdiv::E)
      .def_readonly("F", &PyEtdiv::F)
      .def("__len__", [](const PyEtdiv&) { return 6; })
      .def("__getitem__", [](const PyEtdiv& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.A);
        if (i == 1)
          return py::cast(s.B);
        if (i == 2)
          return py::cast(s.C);
        if (i == 3)
          return py::cast(s.D);
        if (i == 4)
          return py::cast(s.E);
        if (i == 5)
          return py::cast(s.F);
        throw py::index_error();
      });
  m.def(
      "etdiv",
      &python_etdiv,
      py::arg("A"),
      py::arg("B"),
      py::arg("C"),
      py::arg("D"),
      py::arg("E"),
      py::arg("F"),
      R"""(Parameters
  ----------
  A : 
  B : 
  C : 
  D : 
  E : 
  F : 
  )""");
  py::class_<
      Bmad::EvaluateArrayIndex,
      std::unique_ptr<Bmad::EvaluateArrayIndex>>(
      m, "EvaluateArrayIndex", "evaluate_array_index return type")
      .def_readonly("err_flag", &Bmad::EvaluateArrayIndex::err_flag)
      .def_readonly("word2", &Bmad::EvaluateArrayIndex::word2)
      .def_readonly("delim2", &Bmad::EvaluateArrayIndex::delim2)
      .def_readonly("this_index", &Bmad::EvaluateArrayIndex::this_index)
      .def("__len__", [](const Bmad::EvaluateArrayIndex&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::EvaluateArrayIndex& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.word2);
            if (i == 2)
              return py::cast(s.delim2);
            if (i == 3)
              return py::cast(s.this_index);
            throw py::index_error();
          });
  m.def(
      "evaluate_array_index",
      &Bmad::evaluate_array_index,
      py::arg("delim_list1"),
      py::arg("delim_list2"),
      R"""(Function evaluate_array_index (err_flag, delim_list1, word2, delim_list2, delim2) result (this_index)

  Function of evaluate the index of an array. Typically the text being parsed looks like:
       "5) = ..."         or
       "6).COMP = ..."

  Parameters
  ----------
  delim_list1 : unknown
      Delimitor after the integer. Normally ')'.
  delim_list2 : unknown
      Delimitor list to mark the end of word2. Normally '='.

  Returns
  -------
  err_flag : bool
      Set True if there is an error. False otherwise.
  word2 : unknown
      Word found after delim1. Normally this should be blank.
  delim2 : unknown
      Actual delimitor found after word2.
  this_index : int
      Integer value
  )""");
  py::class_<Bmad::EvaluateLogical, std::unique_ptr<Bmad::EvaluateLogical>>(
      m, "EvaluateLogical", "evaluate_logical return type")
      .def_readonly("iostat", &Bmad::EvaluateLogical::iostat)
      .def_readonly("this_logic", &Bmad::EvaluateLogical::this_logic)
      .def("__len__", [](const Bmad::EvaluateLogical&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::EvaluateLogical& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.iostat);
            if (i == 1)
              return py::cast(s.this_logic);
            throw py::index_error();
          });
  m.def(
      "evaluate_logical",
      &Bmad::evaluate_logical,
      py::arg("word"),
      R"""(Function evaluate_logical (word, iostat) result (this_logic)

  Function of convert a string into a logical value.
  Accepted possibilities are:
    .TRUE.  .FALSE.
     TRUE    FALSE
     T       F

  Parameters
  ----------
  word : unknown
      Input string.

  Returns
  -------
  this_logic : bool
      Result.
  iostat : int
      Status: Returns 0 if conversion successful.
  )""");
  m.def(
      "exact_bend_edge_kick",
      &Bmad::exact_bend_edge_kick,
      py::arg("ele"),
      py::arg("param"),
      py::arg("particle_at"),
      py::arg("orb"),
      py::arg("mat6") = py::none(),
      py::arg("make_matrix") = py::none(),
      R"""(Subroutine exact_bend_edge_kick (ele, param, particle_at, orb, mat6, make_matrix)

  Subroutine to track through the edge field of an sbend.
  Uses routines adapted from PTC

  Parameters
  ----------
  ele : EleStruct
      SBend element.
  param : LatParamStruct
  particle_at : int
      first_track_edge$, or second_track_edge$.
  orb : CoordStruct
      Coords after tracking.
  mat6 : float, optional
      Transfer matrix up to the edge.
      This parameter is an input/output and is modified in-place. As an output: Transfer matrix through the
      edge.
  make_matrix : float, optional
      Propagate the transfer matrix? Default is False.
  )""");
  m.def(
      "exp_bessi0",
      &Bmad::exp_bessi0,
      py::arg("t"),
      py::arg("B1"),
      py::arg("B2"),
      R"""(Function exp_bessi0(t, B1, B2)

  This is essentially the Numercal Recipes bessi0 function multiplied by exp(-B1*t).

  This overcomes an issue where exp(B2*t) may be huge and exp(-B1*t) may be small.
  Evaluating exp(B2*t) may result in overflow, but exp((B2-B1)*t) has a moderate value.
  Simplifying the algebra of B2-B1 suggests that is should always have a moderate magnitude.

  Parameters
  ----------
  t : float
      Scalar agrument to evaluate function at.
  B1 : float
      Scalar value.  Eq. 33 from Piwinski's paper.
  B2 : float
      Scalar value.  Eq. 34 from Piwinski's paper. <return value> -- Real(rp): Scalar return value.
  )""");
  py::class_<PyExpectOneOf, std::unique_ptr<PyExpectOneOf>>(
      m, "ExpectOneOf", "expect_one_of return type")
      .def_readonly("delim", &PyExpectOneOf::delim)
      .def_readonly("delim_found", &PyExpectOneOf::delim_found)
      .def_readonly("is_ok", &PyExpectOneOf::is_ok)
      .def("__len__", [](const PyExpectOneOf&) { return 3; })
      .def("__getitem__", [](const PyExpectOneOf& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "expect_one_of",
      &python_expect_one_of,
      py::arg("delim_list"),
      py::arg("check_input_delim"),
      py::arg("ele_name"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("is_ok"),
      R"""(Function expect_one_of (delim_list, check_input_delim, ele_name, delim, delim_found) result (is_ok)

  Routine to check either that the current delimitor or the next character in the parse stream is the
  expected delimitor.
  This routine is used for Bmad lattice file parsing and is not meant for general use.

  Also see: expect_this

  Parameters
  ----------
  delim_list : unknown
      List of expected (valid) delimitors. If list contains a space character then no delimitor (indicating the
      end of the command) is a valid possibility.
  check_input_delim : bool
      If True, then check if delim argument is in the delim_list. -- logical: If True, then check if delim
      argument is in the delim_list. If False, check that the next character in the parse stream is an expected
      delimitor.
  ele_name : unknown
      Lattice element under construction. Used for error messages.
  delim : unknown
      Current delimitor that will be checked if check_input_delim = .true.
      This parameter is an input/output and is modified in-place. As an output: Next delim if check_input_delim
      = False.

  Returns
  -------
  is_ok

  Notes
  -----
  Related routines:
  expect_this
  )""");
  py::class_<Bmad::ExpectThis, std::unique_ptr<Bmad::ExpectThis>>(
      m, "ExpectThis", "expect_this return type")
      .def_readonly("delim", &Bmad::ExpectThis::delim)
      .def_readonly("delim_found", &Bmad::ExpectThis::delim_found)
      .def_readonly("is_ok", &Bmad::ExpectThis::is_ok)
      .def("__len__", [](const Bmad::ExpectThis&) { return 3; })
      .def("__getitem__", [](const Bmad::ExpectThis& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "expect_this",
      &Bmad::expect_this,
      py::arg("expecting"),
      py::arg("check_delim"),
      py::arg("call_check"),
      py::arg("err_str"),
      py::arg("ele"),
      R"""(Function expect_this (expecting, check_delim, call_check, err_str, ele, delim, delim_found) result (is_ok)

  Checks that the next character or characters in the parse stream corresponds to the
  characters in the expecting argument. For example, if expecting is ')={' these three characters
  should be the next non-blank characters in the parse stream.

  Also see: expect_one_of

  Parameters
  ----------
  expecting : unknown
      list of characters that are expected to be next in the parse stream.
  check_delim : bool
      If True then use delim argument as first token to check. A blank character indicates end of command is
      expected.
  call_check : bool
      If True then check for 'call::<filename>' construct.
  err_str : unknown
      String used for error messages.
  ele : EleStruct
      Element parameters being parsed.

  Returns
  -------
  delim : unknown
      Final delim
  delim_found : bool
      Is there a final delim (as opposed to end of command).

  Notes
  -----
  Related routines:
  expect_one_of
  )""");
  m.def(
      "expression_stack_to_string",
      &Bmad::expression_stack_to_string,
      py::arg("stack"),
      py::arg("polish") = py::none(),
      R"""(Function expression_stack_to_string (stack, polish) result (str)

  Routine to convert an expression stack to a string

  Parameters
  ----------
  stack : ExpressionAtomStruct
      arithmetic expression
  polish : , optional
      logical, optional, Construct expression in reverse polish? Default is False.

  Returns
  -------
  str : unknown
      : Expression in string form.
  )""");
  py::class_<
      Bmad::ExpressionStackValue,
      std::unique_ptr<Bmad::ExpressionStackValue>>(
      m, "ExpressionStackValue", "expression_stack_value return type")
      .def_readonly("err_flag", &Bmad::ExpressionStackValue::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionStackValue::err_str)
      .def_readonly("value", &Bmad::ExpressionStackValue::value)
      .def("__len__", [](const Bmad::ExpressionStackValue&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::ExpressionStackValue& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.err_str);
            if (i == 2)
              return py::cast(s.value);
            throw py::index_error();
          });
  m.def(
      "expression_stack_value",
      &Bmad::expression_stack_value,
      py::arg("stack"),
      py::arg("var") = py::none(),
      py::arg("use_old") = py::none(),
      R"""(Function expression_stack_value (stack, err_flag, err_str, var, use_old) result (value)

  Routine to evaluate a mathematical expression represented by an "expression stack".
  Expression stacks are created by expression_string_to_stack.

  Note: Stack elements with stack(i)%type == variable$ need to be evalauated before
  calling this routine and the value placed in stack(i)%value.

  Parameters
  ----------
  stack : ExpressionAtomStruct
      Expression to evaluate.
  var : ControlVar1Struct, optional
      Array of control variables. Used with Bmad controller elements.
  use_old : bool, optional
      Use var.old_value? Must be present if var(:) is present.

  Returns
  -------
  value : float
      Value of the expression.
  err_flag : bool
      True if there is an evaluation problem. False otherwise.
  err_str : unknown
      Error string explaining error if there is one.

  Notes
  -----
  Related routines:
  expression_value expression_string_to_stack
  )""");
  py::class_<
      Bmad::ExpressionStringToStack,
      std::unique_ptr<Bmad::ExpressionStringToStack>>(
      m, "ExpressionStringToStack", "expression_string_to_stack return type")
      .def_readonly("stack", &Bmad::ExpressionStringToStack::stack)
      .def_readonly("n_stack", &Bmad::ExpressionStringToStack::n_stack)
      .def_readonly("err_flag", &Bmad::ExpressionStringToStack::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionStringToStack::err_str)
      .def("__len__", [](const Bmad::ExpressionStringToStack&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::ExpressionStringToStack& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.stack);
            if (i == 1)
              return py::cast(s.n_stack);
            if (i == 2)
              return py::cast(s.err_flag);
            if (i == 3)
              return py::cast(s.err_str);
            throw py::index_error();
          });
  m.def(
      "expression_string_to_stack",
      &Bmad::expression_string_to_stack,
      py::arg("string"),
      R"""(Subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)

  This routine creates an expression stack array which can be used
  to evaluate an arithmethic expression.

  Stack end elements not used are marked stack(i)%type = end_stack$

  Stack elements with stack(i)%type = variable$ are elements that need
  to be evaluated before calling expression_stack_value.

  Parameters
  ----------
  string : unknown
      Expression to be converted.

  Returns
  -------
  stack : ExpressionAtomStruct
      Expression evaluation stack.
  n_stack : int
      number of "atoms" used by the expression
  err_flag : bool
      Set True if there is an error (EG divide by 0).
  err_str : unknown
      String describing the error.

  Notes
  -----
  Related routines:
  expression_value expression_stack_value
  )""");
  py::class_<
      Bmad::ExpressionStringToTree,
      std::unique_ptr<Bmad::ExpressionStringToTree>>(
      m, "ExpressionStringToTree", "expression_string_to_tree return type")
      .def_readonly("err_flag", &Bmad::ExpressionStringToTree::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionStringToTree::err_str)
      .def("__len__", [](const Bmad::ExpressionStringToTree&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::ExpressionStringToTree& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.err_str);
            throw py::index_error();
          });
  m.def(
      "expression_string_to_tree",
      &Bmad::expression_string_to_tree,
      py::arg("string"),
      py::arg("root_tree"),
      R"""(Subroutine expression_string_to_tree (string, root_tree, err_flag, err_str)

  Routine to create an expression tree array which can be used
  to evaluate an arithmethic expression.

  Parameters
  ----------
  root_tree : ExpressionTreeStruct
      Only used when recursively called.
  string : unknown
      Expression to be converted.

  Returns
  -------
  tree : ExpressionTreeStruct
      Expression evaluation tree.
  err_flag : bool
      Set True if there is an error (EG divide by 0).
  err_str : unknown
      String describing the error. Make length large to hold the expression.

  Notes
  -----
  Related routines:
  expression_value expression_tree_value deallocate_expression_tree Important! trees use pointers as opposed to
  allocatable arrays due to the ifort compiler not being able to ) being an allocatable array. Thus
  deallocate_expression_tree must be called before any tree instance goes out of scope. plus$ minus$ times$
  divide$ power$ unary_minus$ unary_plus$ constant$ numeric$ variable$ function$ root$ parens$ func_parens$
  square_brackets$ curly_brackets$ arrow$ equal$ colon$ double_colon$ vertical_bar$ compound$ "->" "::" + - * /
  ^ = : & [] () {} Root node name is "root" is of type root$ Brackets in the expression string must be matched.
  "[]" / square_brackets$ "()" / parens$ func_parens$ "{}" / curley_brackets$ The root node equal nodes all
  bracket nodes will have an array of child nodes all of which will be comma nodes. "[A B]" will translate to a
  "[]" node with two comma children the first comma child will have a single child "A" the second comma child
  will have a single child "B". "(A)" will translate to a "()" node with one comma child this comma child will
  have a single child "A". If the string is an equation. For example "A B = C D Z". In this case the root node
  will have two equal node children (and not comma children) The first equal node represents the left hand side
  of the equation this node will have two comma children. The second equal node child will have three comma
  children. :orxit.x" (this is a Tao construct) which get :" "orbit.x" functions line "atan()" are considered
  compound vars with children "atan" "()" The funciton argument of a species related function like "He++" in the
  construct "mass_of(He++)" will not get split will get marked as a species_const$.
  )""");
  m.def(
      "expression_tree_to_string",
      &Bmad::expression_tree_to_string,
      py::arg("tree"),
      py::arg("include_root") = py::none(),
      py::arg("n_node") = py::none(),
      py::arg("parent") = py::none(),
      R"""(Function expression_tree_to_string (tree, include_root, n_node, parent) result(str_out)

  Routine to convert an expression tree to a expression string.

  Parameters
  ----------
  tree : ExpressionTreeStruct
      Root of tree to print.
  include_root : bool, optional
      Default is True. If True, do not inculde in the output string the root node. Note: If the root node is of
      type root$, this node is always ignored.
  n_node : int, optional
      Node index. parent.node(n_node) === tree. Internal use only. Used with recursive calls.
  parent : ExpressionTreeStruct, optional
      Internal use only. Used with recusive calls.

  Returns
  -------
  str_out : unknown
      Expression string.
  )""");
  py::class_<Bmad::ExpressionValue, std::unique_ptr<Bmad::ExpressionValue>>(
      m, "ExpressionValue", "expression_value return type")
      .def_readonly("err_flag", &Bmad::ExpressionValue::err_flag)
      .def_readonly("err_str", &Bmad::ExpressionValue::err_str)
      .def_readonly("value", &Bmad::ExpressionValue::value)
      .def("__len__", [](const Bmad::ExpressionValue&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::ExpressionValue& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.err_str);
            if (i == 2)
              return py::cast(s.value);
            throw py::index_error();
          });
  m.def(
      "expression_value",
      &Bmad::expression_value,
      py::arg("expression"),
      py::arg("var") = py::none(),
      py::arg("use_old") = py::none(),
      R"""(Function expression_value (expression, err_flag, err_str, var, use_old) result (value)

  Routine to evaluate a mathematical expression encoded in a string.

  Parameters
  ----------
  expression : unknown
      Expression string.
  var : ControlVar1Struct, optional
      Array of control variables. Used with Bmad controller elements.
  use_old : bool, optional
      Use var.old_value? Must be present if var(:) is present.

  Returns
  -------
  value : float
      Value of the expression.
  err_flag : bool
      True if there is an evaluation problem. False otherwise.
  err_str : unknown
      Error string explaining error if there is one.

  Notes
  -----
  Related routines:
  expression_string_to_stack expression_stack_value
  )""");
}
