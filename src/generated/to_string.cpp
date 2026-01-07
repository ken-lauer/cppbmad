// clang-format off
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
// clang-format on
// TODO: some issue where clang-format is rearranging these

namespace Bmad {
std::string to_string(const SplineProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SplineProxy",
      {std::pair{"x0", to_string(self.x0())},
       std::pair{"y0", to_string(self.y0())},
       std::pair{"x1", to_string(self.x1())},
       std::pair{"coef", to_string(self.coef())}});
}
std::string to_string(const SpinPolarProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SpinPolarProxy",
      {std::pair{"polarization", to_string(self.polarization())},
       std::pair{"theta", to_string(self.theta())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"xi", to_string(self.xi())}});
}
std::string to_string(const AcKickerTimeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "AcKickerTimeProxy",
      {std::pair{"amp", to_string(self.amp())},
       std::pair{"time", to_string(self.time())},
       std::pair{"spline", to_string(self.spline())}});
}
std::string to_string(const AcKickerFreqProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "AcKickerFreqProxy",
      {std::pair{"f", to_string(self.f())},
       std::pair{"amp", to_string(self.amp())},
       std::pair{"phi", to_string(self.phi())}});
}
std::string to_string(const AcKickerProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "AcKickerProxy",
      {std::pair{"amp_vs_time", "[...]"}, std::pair{"frequency", "[...]"}});
}
std::string to_string(const Interval1CoefProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "Interval1CoefProxy",
      {std::pair{"c0", to_string(self.c0())},
       std::pair{"c1", to_string(self.c1())},
       std::pair{"n_exp", to_string(self.n_exp())}});
}
std::string to_string(const PhotonReflectTableProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PhotonReflectTableProxy",
      {std::pair{"angle", to_string(self.angle())},
       std::pair{"energy", to_string(self.energy())},
       std::pair{"int1", "[...]"},
       std::pair{"p_reflect", to_string(self.p_reflect())},
       std::pair{"max_energy", to_string(self.max_energy())},
       std::pair{"p_reflect_scratch", to_string(self.p_reflect_scratch())},
       std::pair{"bragg_angle", to_string(self.bragg_angle())}});
}
std::string to_string(const PhotonReflectSurfaceProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PhotonReflectSurfaceProxy",
      {std::pair{"name", self.name()},
       std::pair{"description", self.description()},
       std::pair{"reflectivity_file", self.reflectivity_file()},
       std::pair{"table", "[...]"},
       std::pair{
           "surface_roughness_rms", to_string(self.surface_roughness_rms())},
       std::pair{
           "roughness_correlation_len",
           to_string(self.roughness_correlation_len())},
       std::pair{"ix_surface", to_string(self.ix_surface())}});
}
std::string to_string(const CoordProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CoordProxy",
      {std::pair{"vec", to_string(self.vec())},
       std::pair{"s", to_string(self.s())},
       std::pair{"t", to_string(self.t())},
       std::pair{"spin", to_string(self.spin())},
       std::pair{"field", to_string(self.field())},
       std::pair{"phase", to_string(self.phase())},
       std::pair{"charge", to_string(self.charge())},
       std::pair{"dt_ref", to_string(self.dt_ref())},
       std::pair{"r", to_string(self.r())},
       std::pair{"p0c", to_string(self.p0c())},
       std::pair{"E_potential", to_string(self.E_potential())},
       std::pair{"beta", to_string(self.beta())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_turn", to_string(self.ix_turn())},
       std::pair{"ix_user", to_string(self.ix_user())},
       std::pair{"state", to_string(self.state())},
       std::pair{"direction", to_string(self.direction())},
       std::pair{"time_dir", to_string(self.time_dir())},
       std::pair{"species", to_string(self.species())},
       std::pair{"location", to_string(self.location())}});
}
std::string to_string(const CoordArrayProxy& self) {
  return repr(
      self.get_fortran_ptr(), "CoordArrayProxy", {std::pair{"orbit", "[...]"}});
}
std::string to_string(const BpmPhaseCouplingProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BpmPhaseCouplingProxy",
      {std::pair{"K_22a", to_string(self.K_22a())},
       std::pair{"K_12a", to_string(self.K_12a())},
       std::pair{"K_11b", to_string(self.K_11b())},
       std::pair{"K_12b", to_string(self.K_12b())},
       std::pair{"Cbar22_a", to_string(self.Cbar22_a())},
       std::pair{"Cbar12_a", to_string(self.Cbar12_a())},
       std::pair{"Cbar11_b", to_string(self.Cbar11_b())},
       std::pair{"Cbar12_b", to_string(self.Cbar12_b())},
       std::pair{"phi_a", to_string(self.phi_a())},
       std::pair{"phi_b", to_string(self.phi_b())}});
}
std::string to_string(const ExpressionAtomProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ExpressionAtomProxy",
      {std::pair{"name", self.name()},
       std::pair{"type", to_string(self.type())},
       std::pair{"value", to_string(self.value())}});
}
std::string to_string(const WakeSrZLongProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "WakeSrZLongProxy",
      {std::pair{"w", to_string(self.w())},
       std::pair{"fw", to_string(self.fw())},
       std::pair{"fbunch", to_string(self.fbunch())},
       std::pair{"w_out", to_string(self.w_out())},
       std::pair{"dz", to_string(self.dz())},
       std::pair{"z0", to_string(self.z0())},
       std::pair{"smoothing_sigma", to_string(self.smoothing_sigma())},
       std::pair{"position_dependence", to_string(self.position_dependence())},
       std::pair{"time_based", to_string(self.time_based())}});
}
std::string to_string(const WakeSrModeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "WakeSrModeProxy",
      {std::pair{"amp", to_string(self.amp())},
       std::pair{"damp", to_string(self.damp())},
       std::pair{"k", to_string(self.k())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"b_sin", to_string(self.b_sin())},
       std::pair{"b_cos", to_string(self.b_cos())},
       std::pair{"a_sin", to_string(self.a_sin())},
       std::pair{"a_cos", to_string(self.a_cos())},
       std::pair{"polarization", to_string(self.polarization())},
       std::pair{
           "position_dependence", to_string(self.position_dependence())}});
}
std::string to_string(const WakeSrProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "WakeSrProxy",
      {std::pair{"file", self.file()},
       std::pair{"z_long", to_string(self.z_long())},
       std::pair{"long_wake", "[...]"},
       std::pair{"trans_wake", "[...]"},
       std::pair{"z_ref_long", to_string(self.z_ref_long())},
       std::pair{"z_ref_trans", to_string(self.z_ref_trans())},
       std::pair{"z_max", to_string(self.z_max())},
       std::pair{"amp_scale", to_string(self.amp_scale())},
       std::pair{"z_scale", to_string(self.z_scale())},
       std::pair{"scale_with_length", to_string(self.scale_with_length())}});
}
std::string to_string(const WakeLrModeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "WakeLrModeProxy",
      {std::pair{"freq", to_string(self.freq())},
       std::pair{"freq_in", to_string(self.freq_in())},
       std::pair{"R_over_Q", to_string(self.R_over_Q())},
       std::pair{"Q", to_string(self.Q())},
       std::pair{"damp", to_string(self.damp())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"angle", to_string(self.angle())},
       std::pair{"b_sin", to_string(self.b_sin())},
       std::pair{"b_cos", to_string(self.b_cos())},
       std::pair{"a_sin", to_string(self.a_sin())},
       std::pair{"a_cos", to_string(self.a_cos())},
       std::pair{"m", to_string(self.m())},
       std::pair{"polarized", to_string(self.polarized())}});
}
std::string to_string(const WakeLrProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "WakeLrProxy",
      {std::pair{"file", self.file()},
       std::pair{"mode", "[...]"},
       std::pair{"t_ref", to_string(self.t_ref())},
       std::pair{"freq_spread", to_string(self.freq_spread())},
       std::pair{"amp_scale", to_string(self.amp_scale())},
       std::pair{"time_scale", to_string(self.time_scale())},
       std::pair{"self_wake_on", to_string(self.self_wake_on())}});
}
std::string to_string(const LatEleLocProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LatEleLocProxy",
      {std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_branch", to_string(self.ix_branch())}});
}
std::string to_string(const WakeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "WakeProxy",
      {std::pair{"sr", to_string(self.sr())},
       std::pair{"lr", to_string(self.lr())}});
}
std::string to_string(const TaylorTermProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaylorTermProxy",
      {std::pair{"coef", to_string(self.coef())},
       std::pair{"expn", to_string(self.expn())}});
}
std::string to_string(const TaylorProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaylorProxy",
      {std::pair{"ref", to_string(self.ref())}, std::pair{"term", "[...]"}});
}
std::string to_string(const EmTaylorTermProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "EmTaylorTermProxy",
      {std::pair{"coef", to_string(self.coef())},
       std::pair{"expn", to_string(self.expn())}});
}
std::string to_string(const EmTaylorProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "EmTaylorProxy",
      {std::pair{"ref", to_string(self.ref())}, std::pair{"term", "[...]"}});
}
std::string to_string(const CartesianMapTerm1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CartesianMapTerm1Proxy",
      {std::pair{"coef", to_string(self.coef())},
       std::pair{"kx", to_string(self.kx())},
       std::pair{"ky", to_string(self.ky())},
       std::pair{"kz", to_string(self.kz())},
       std::pair{"x0", to_string(self.x0())},
       std::pair{"y0", to_string(self.y0())},
       std::pair{"phi_z", to_string(self.phi_z())},
       std::pair{"family", to_string(self.family())},
       std::pair{"form", to_string(self.form())}});
}
std::string to_string(const CartesianMapTermProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CartesianMapTermProxy",
      {std::pair{"file", self.file()},
       std::pair{"n_link", to_string(self.n_link())},
       std::pair{"term", "[...]"}});
}
std::string to_string(const CartesianMapProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CartesianMapProxy",
      {std::pair{"field_scale", to_string(self.field_scale())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"master_parameter", to_string(self.master_parameter())},
       std::pair{"ele_anchor_pt", to_string(self.ele_anchor_pt())},
       std::pair{"field_type", to_string(self.field_type())},
       std::pair{"ptr", to_string(self.ptr())}});
}
std::string to_string(const CylindricalMapTerm1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CylindricalMapTerm1Proxy",
      {std::pair{"e_coef", to_string(self.e_coef())},
       std::pair{"b_coef", to_string(self.b_coef())}});
}
std::string to_string(const CylindricalMapTermProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CylindricalMapTermProxy",
      {std::pair{"file", self.file()},
       std::pair{"n_link", to_string(self.n_link())},
       std::pair{"term", "[...]"}});
}
std::string to_string(const CylindricalMapProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "CylindricalMapProxy",
      {std::pair{"m", to_string(self.m())},
       std::pair{"harmonic", to_string(self.harmonic())},
       std::pair{"phi0_fieldmap", to_string(self.phi0_fieldmap())},
       std::pair{"theta0_azimuth", to_string(self.theta0_azimuth())},
       std::pair{"field_scale", to_string(self.field_scale())},
       std::pair{"master_parameter", to_string(self.master_parameter())},
       std::pair{"ele_anchor_pt", to_string(self.ele_anchor_pt())},
       std::pair{"dz", to_string(self.dz())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"ptr", to_string(self.ptr())}});
}
std::string to_string(const BicubicCmplxCoefProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BicubicCmplxCoefProxy",
      {std::pair{"coef", to_string(self.coef())},
       std::pair{"i_box", to_string(self.i_box())}});
}
std::string to_string(const TricubicCmplxCoefProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TricubicCmplxCoefProxy",
      {std::pair{"coef", to_string(self.coef())},
       std::pair{"i_box", to_string(self.i_box())}});
}
std::string to_string(const GridFieldPt1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "GridFieldPt1Proxy",
      {std::pair{"E", to_string(self.E())},
       std::pair{"B", to_string(self.B())}});
}
std::string to_string(const GridFieldPtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "GridFieldPtProxy",
      {std::pair{"file", self.file()},
       std::pair{"n_link", to_string(self.n_link())},
       std::pair{"pt", "[...]"}});
}
std::string to_string(const GridFieldProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "GridFieldProxy",
      {std::pair{"geometry", to_string(self.geometry())},
       std::pair{"harmonic", to_string(self.harmonic())},
       std::pair{"phi0_fieldmap", to_string(self.phi0_fieldmap())},
       std::pair{"field_scale", to_string(self.field_scale())},
       std::pair{"field_type", to_string(self.field_type())},
       std::pair{"master_parameter", to_string(self.master_parameter())},
       std::pair{"ele_anchor_pt", to_string(self.ele_anchor_pt())},
       std::pair{"interpolation_order", to_string(self.interpolation_order())},
       std::pair{"dr", to_string(self.dr())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"curved_ref_frame", to_string(self.curved_ref_frame())},
       std::pair{"ptr", to_string(self.ptr())},
       std::pair{"bi_coef", "[...]"},
       std::pair{"tri_coef", "[...]"}});
}
std::string to_string(const FloorPositionProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "FloorPositionProxy",
      {std::pair{"r", to_string(self.r())},
       std::pair{"w", to_string(self.w())},
       std::pair{"theta", to_string(self.theta())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"psi", to_string(self.psi())}});
}
std::string to_string(const HighEnergySpaceChargeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "HighEnergySpaceChargeProxy",
      {std::pair{"closed_orb", to_string(self.closed_orb())},
       std::pair{"kick_const", to_string(self.kick_const())},
       std::pair{"sig_x", to_string(self.sig_x())},
       std::pair{"sig_y", to_string(self.sig_y())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"sin_phi", to_string(self.sin_phi())},
       std::pair{"cos_phi", to_string(self.cos_phi())},
       std::pair{"sig_z", to_string(self.sig_z())}});
}
std::string to_string(const XyDispProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "XyDispProxy",
      {std::pair{"eta", to_string(self.eta())},
       std::pair{"etap", to_string(self.etap())},
       std::pair{"deta_ds", to_string(self.deta_ds())},
       std::pair{"sigma", to_string(self.sigma())},
       std::pair{"deta_dpz", to_string(self.deta_dpz())},
       std::pair{"detap_dpz", to_string(self.detap_dpz())}});
}
std::string to_string(const TwissProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TwissProxy",
      {std::pair{"beta", to_string(self.beta())},
       std::pair{"alpha", to_string(self.alpha())},
       std::pair{"gamma", to_string(self.gamma())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"eta", to_string(self.eta())},
       std::pair{"etap", to_string(self.etap())},
       std::pair{"deta_ds", to_string(self.deta_ds())},
       std::pair{"sigma", to_string(self.sigma())},
       std::pair{"sigma_p", to_string(self.sigma_p())},
       std::pair{"emit", to_string(self.emit())},
       std::pair{"norm_emit", to_string(self.norm_emit())},
       std::pair{"chrom", to_string(self.chrom())},
       std::pair{"dbeta_dpz", to_string(self.dbeta_dpz())},
       std::pair{"dalpha_dpz", to_string(self.dalpha_dpz())},
       std::pair{"deta_dpz", to_string(self.deta_dpz())},
       std::pair{"detap_dpz", to_string(self.detap_dpz())}});
}
std::string to_string(const Mode3Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "Mode3Proxy",
      {std::pair{"v", to_string(self.v())},
       std::pair{"a", to_string(self.a())},
       std::pair{"b", to_string(self.b())},
       std::pair{"c", to_string(self.c())},
       std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())}});
}
std::string to_string(const BookkeepingStateProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BookkeepingStateProxy",
      {std::pair{"attributes", to_string(self.attributes())},
       std::pair{"control", to_string(self.control())},
       std::pair{"floor_position", to_string(self.floor_position())},
       std::pair{"s_position", to_string(self.s_position())},
       std::pair{"ref_energy", to_string(self.ref_energy())},
       std::pair{"mat6", to_string(self.mat6())},
       std::pair{"rad_int", to_string(self.rad_int())},
       std::pair{"ptc", to_string(self.ptc())},
       std::pair{"has_misalign", to_string(self.has_misalign())}});
}
std::string to_string(const RadMapProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RadMapProxy",
      {std::pair{"ref_orb", to_string(self.ref_orb())},
       std::pair{"damp_dmat", to_string(self.damp_dmat())},
       std::pair{"xfer_damp_vec", to_string(self.xfer_damp_vec())},
       std::pair{"xfer_damp_mat", to_string(self.xfer_damp_mat())},
       std::pair{"stoc_mat", to_string(self.stoc_mat())}});
}
std::string to_string(const RadMapEleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RadMapEleProxy",
      {std::pair{"rm0", to_string(self.rm0())},
       std::pair{"rm1", to_string(self.rm1())},
       std::pair{"stale", to_string(self.stale())}});
}
std::string to_string(const GenGrad1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "GenGrad1Proxy",
      {std::pair{"m", to_string(self.m())},
       std::pair{"sincos", to_string(self.sincos())},
       std::pair{"n_deriv_max", to_string(self.n_deriv_max())},
       std::pair{"deriv", to_string(self.deriv())}});
}
std::string to_string(const GenGradMapProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "GenGradMapProxy",
      {std::pair{"file", self.file()},
       std::pair{"gg", "[...]"},
       std::pair{"ele_anchor_pt", to_string(self.ele_anchor_pt())},
       std::pair{"field_type", to_string(self.field_type())},
       std::pair{"iz0", to_string(self.iz0())},
       std::pair{"iz1", to_string(self.iz1())},
       std::pair{"dz", to_string(self.dz())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"field_scale", to_string(self.field_scale())},
       std::pair{"master_parameter", to_string(self.master_parameter())},
       std::pair{"curved_ref_frame", to_string(self.curved_ref_frame())}});
}
std::string to_string(const SurfaceSegmentedPtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceSegmentedPtProxy",
      {std::pair{"x0", to_string(self.x0())},
       std::pair{"y0", to_string(self.y0())},
       std::pair{"z0", to_string(self.z0())},
       std::pair{"dz_dx", to_string(self.dz_dx())},
       std::pair{"dz_dy", to_string(self.dz_dy())}});
}
std::string to_string(const SurfaceSegmentedProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceSegmentedProxy",
      {std::pair{"active", to_string(self.active())},
       std::pair{"dr", to_string(self.dr())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"pt", "[...]"}});
}
std::string to_string(const SurfaceHMisalignPtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceHMisalignPtProxy",
      {std::pair{"x0", to_string(self.x0())},
       std::pair{"y0", to_string(self.y0())},
       std::pair{"rot_y", to_string(self.rot_y())},
       std::pair{"rot_t", to_string(self.rot_t())},
       std::pair{"rot_y_rms", to_string(self.rot_y_rms())},
       std::pair{"rot_t_rms", to_string(self.rot_t_rms())}});
}
std::string to_string(const SurfaceHMisalignProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceHMisalignProxy",
      {std::pair{"active", to_string(self.active())},
       std::pair{"dr", to_string(self.dr())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"pt", "[...]"}});
}
std::string to_string(const SurfaceDisplacementPtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceDisplacementPtProxy",
      {std::pair{"x0", to_string(self.x0())},
       std::pair{"y0", to_string(self.y0())},
       std::pair{"z0", to_string(self.z0())},
       std::pair{"dz_dx", to_string(self.dz_dx())},
       std::pair{"dz_dy", to_string(self.dz_dy())},
       std::pair{"d2z_dxdy", to_string(self.d2z_dxdy())}});
}
std::string to_string(const SurfaceDisplacementProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceDisplacementProxy",
      {std::pair{"active", to_string(self.active())},
       std::pair{"dr", to_string(self.dr())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"pt", "[...]"}});
}
std::string to_string(const TargetPointProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TargetPointProxy",
      {std::pair{"r", to_string(self.r())}});
}
std::string to_string(const SurfaceCurvatureProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SurfaceCurvatureProxy",
      {std::pair{"xy", to_string(self.xy())},
       std::pair{"spherical", to_string(self.spherical())},
       std::pair{"elliptical", to_string(self.elliptical())},
       std::pair{"has_curvature", to_string(self.has_curvature())}});
}
std::string to_string(const PhotonTargetProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PhotonTargetProxy",
      {std::pair{"type", to_string(self.type())},
       std::pair{"n_corner", to_string(self.n_corner())},
       std::pair{"ele_loc", to_string(self.ele_loc())},
       std::pair{"corner", "[...]"},
       std::pair{"center", to_string(self.center())}});
}
std::string to_string(const PhotonMaterialProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PhotonMaterialProxy",
      {std::pair{"f0_m1", to_string(self.f0_m1())},
       std::pair{"f0_m2", to_string(self.f0_m2())},
       std::pair{"f_0", to_string(self.f_0())},
       std::pair{"f_h", to_string(self.f_h())},
       std::pair{"f_hbar", to_string(self.f_hbar())},
       std::pair{"f_hkl", to_string(self.f_hkl())},
       std::pair{"h_norm", to_string(self.h_norm())},
       std::pair{"l_ref", to_string(self.l_ref())}});
}
std::string to_string(const PixelPtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PixelPtProxy",
      {std::pair{"n_photon", to_string(self.n_photon())},
       std::pair{"E_x", to_string(self.E_x())},
       std::pair{"E_y", to_string(self.E_y())},
       std::pair{"intensity_x", to_string(self.intensity_x())},
       std::pair{"intensity_y", to_string(self.intensity_y())},
       std::pair{"intensity", to_string(self.intensity())},
       std::pair{"orbit", to_string(self.orbit())},
       std::pair{"orbit_rms", to_string(self.orbit_rms())},
       std::pair{"init_orbit", to_string(self.init_orbit())},
       std::pair{"init_orbit_rms", to_string(self.init_orbit_rms())}});
}
std::string to_string(const PixelDetecProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PixelDetecProxy",
      {std::pair{"dr", to_string(self.dr())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"n_track_tot", to_string(self.n_track_tot())},
       std::pair{"n_hit_detec", to_string(self.n_hit_detec())},
       std::pair{"n_hit_pixel", to_string(self.n_hit_pixel())},
       std::pair{"pt", "[...]"}});
}
std::string to_string(const PhotonElementProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PhotonElementProxy",
      {std::pair{"curvature", to_string(self.curvature())},
       std::pair{"target", to_string(self.target())},
       std::pair{"material", to_string(self.material())},
       std::pair{"segmented", to_string(self.segmented())},
       std::pair{"h_misalign", to_string(self.h_misalign())},
       std::pair{"displacement", to_string(self.displacement())},
       std::pair{"pixel", to_string(self.pixel())},
       std::pair{
           "reflectivity_table_type",
           to_string(self.reflectivity_table_type())},
       std::pair{
           "reflectivity_table_sigma",
           to_string(self.reflectivity_table_sigma())},
       std::pair{
           "reflectivity_table_pi", to_string(self.reflectivity_table_pi())},
       std::pair{"init_energy_prob", "[...]"},
       std::pair{
           "integrated_init_energy_prob",
           to_string(self.integrated_init_energy_prob())}});
}
std::string to_string(const Wall3dVertexProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "Wall3dVertexProxy",
      {std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"radius_x", to_string(self.radius_x())},
       std::pair{"radius_y", to_string(self.radius_y())},
       std::pair{"tilt", to_string(self.tilt())},
       std::pair{"angle", to_string(self.angle())},
       std::pair{"x0", to_string(self.x0())},
       std::pair{"y0", to_string(self.y0())},
       std::pair{"type", to_string(self.type())}});
}
std::string to_string(const Wall3dSectionProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "Wall3dSectionProxy",
      {std::pair{"name", self.name()},
       std::pair{"material", self.material()},
       std::pair{"v", "[...]"},
       std::pair{"surface", to_string(self.surface())},
       std::pair{"type", to_string(self.type())},
       std::pair{"n_vertex_input", to_string(self.n_vertex_input())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"vertices_state", to_string(self.vertices_state())},
       std::pair{"patch_in_region", to_string(self.patch_in_region())},
       std::pair{"thickness", to_string(self.thickness())},
       std::pair{"s", to_string(self.s())},
       std::pair{"r0", to_string(self.r0())},
       std::pair{"dx0_ds", to_string(self.dx0_ds())},
       std::pair{"dy0_ds", to_string(self.dy0_ds())},
       std::pair{"x0_coef", to_string(self.x0_coef())},
       std::pair{"y0_coef", to_string(self.y0_coef())},
       std::pair{"dr_ds", to_string(self.dr_ds())},
       std::pair{"p1_coef", to_string(self.p1_coef())},
       std::pair{"p2_coef", to_string(self.p2_coef())}});
}
std::string to_string(const Wall3dProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "Wall3dProxy",
      {std::pair{"name", self.name()},
       std::pair{"type", to_string(self.type())},
       std::pair{"ix_wall3d", to_string(self.ix_wall3d())},
       std::pair{"n_link", to_string(self.n_link())},
       std::pair{"thickness", to_string(self.thickness())},
       std::pair{"clear_material", self.clear_material()},
       std::pair{"opaque_material", self.opaque_material()},
       std::pair{"superimpose", to_string(self.superimpose())},
       std::pair{"ele_anchor_pt", to_string(self.ele_anchor_pt())},
       std::pair{"section", "[...]"}});
}
std::string to_string(const RamperLordProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RamperLordProxy",
      {std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_con", to_string(self.ix_con())},
       std::pair{"attrib_ptr", to_string(self.attrib_ptr())}});
}
std::string to_string(const ControlProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ControlProxy",
      {std::pair{"value", to_string(self.value())},
       std::pair{"y_knot", to_string(self.y_knot())},
       std::pair{"stack", "[...]"},
       std::pair{"slave", to_string(self.slave())},
       std::pair{"lord", to_string(self.lord())},
       std::pair{"slave_name", self.slave_name()},
       std::pair{"attribute", self.attribute()},
       std::pair{"ix_attrib", to_string(self.ix_attrib())}});
}
std::string to_string(const ControlVar1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ControlVar1Proxy",
      {std::pair{"name", self.name()},
       std::pair{"value", to_string(self.value())},
       std::pair{"old_value", to_string(self.old_value())}});
}
std::string to_string(const ControlRamp1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ControlRamp1Proxy",
      {std::pair{"y_knot", to_string(self.y_knot())},
       std::pair{"stack", "[...]"},
       std::pair{"attribute", self.attribute()},
       std::pair{"slave_name", self.slave_name()},
       std::pair{"is_controller", to_string(self.is_controller())}});
}
std::string to_string(const ControllerProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ControllerProxy",
      {std::pair{"var", "[...]"},
       std::pair{"ramp", "[...]"},
       std::pair{"ramper_lord", "[...]"},
       std::pair{"x_knot", to_string(self.x_knot())}});
}
std::string to_string(const EllipseBeamInitProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "EllipseBeamInitProxy",
      {std::pair{"part_per_ellipse", to_string(self.part_per_ellipse())},
       std::pair{"n_ellipse", to_string(self.n_ellipse())},
       std::pair{"sigma_cutoff", to_string(self.sigma_cutoff())}});
}
std::string to_string(const KvBeamInitProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "KvBeamInitProxy",
      {std::pair{"part_per_phi", to_string(self.part_per_phi())},
       std::pair{"n_I2", to_string(self.n_I2())},
       std::pair{"A", to_string(self.A())}});
}
std::string to_string(const GridBeamInitProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "GridBeamInitProxy",
      {std::pair{"n_x", to_string(self.n_x())},
       std::pair{"n_px", to_string(self.n_px())},
       std::pair{"x_min", to_string(self.x_min())},
       std::pair{"x_max", to_string(self.x_max())},
       std::pair{"px_min", to_string(self.px_min())},
       std::pair{"px_max", to_string(self.px_max())}});
}
std::string to_string(const BeamInitProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BeamInitProxy",
      {std::pair{"position_file", self.position_file()},
       std::pair{"distribution_type", to_string(self.distribution_type())},
       std::pair{"spin", to_string(self.spin())},
       std::pair{"ellipse", "[...]"},
       std::pair{"KV", to_string(self.KV())},
       std::pair{"grid", "[...]"},
       std::pair{"center_jitter", to_string(self.center_jitter())},
       std::pair{"emit_jitter", to_string(self.emit_jitter())},
       std::pair{"sig_z_jitter", to_string(self.sig_z_jitter())},
       std::pair{"sig_pz_jitter", to_string(self.sig_pz_jitter())},
       std::pair{"n_particle", to_string(self.n_particle())},
       std::pair{"renorm_center", to_string(self.renorm_center())},
       std::pair{"renorm_sigma", to_string(self.renorm_sigma())},
       std::pair{"random_engine", self.random_engine()},
       std::pair{"random_gauss_converter", self.random_gauss_converter()},
       std::pair{"random_sigma_cutoff", to_string(self.random_sigma_cutoff())},
       std::pair{"a_norm_emit", to_string(self.a_norm_emit())},
       std::pair{"b_norm_emit", to_string(self.b_norm_emit())},
       std::pair{"a_emit", to_string(self.a_emit())},
       std::pair{"b_emit", to_string(self.b_emit())},
       std::pair{"dPz_dz", to_string(self.dPz_dz())},
       std::pair{"center", to_string(self.center())},
       std::pair{"t_offset", to_string(self.t_offset())},
       std::pair{"dt_bunch", to_string(self.dt_bunch())},
       std::pair{"sig_z", to_string(self.sig_z())},
       std::pair{"sig_pz", to_string(self.sig_pz())},
       std::pair{"bunch_charge", to_string(self.bunch_charge())},
       std::pair{"n_bunch", to_string(self.n_bunch())},
       std::pair{"ix_turn", to_string(self.ix_turn())},
       std::pair{"species", self.species()},
       std::pair{
           "full_6D_coupling_calc", to_string(self.full_6D_coupling_calc())},
       std::pair{"use_particle_start", to_string(self.use_particle_start())},
       std::pair{"use_t_coords", to_string(self.use_t_coords())},
       std::pair{"use_z_as_t", to_string(self.use_z_as_t())},
       std::pair{"file_name", self.file_name()}});
}
std::string to_string(const LatParamProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LatParamProxy",
      {std::pair{"n_part", to_string(self.n_part())},
       std::pair{"total_length", to_string(self.total_length())},
       std::pair{"unstable_factor", to_string(self.unstable_factor())},
       std::pair{"t1_with_RF", to_string(self.t1_with_RF())},
       std::pair{"t1_no_RF", to_string(self.t1_no_RF())},
       std::pair{"spin_tune", to_string(self.spin_tune())},
       std::pair{"particle", to_string(self.particle())},
       std::pair{
           "default_tracking_species",
           to_string(self.default_tracking_species())},
       std::pair{"geometry", to_string(self.geometry())},
       std::pair{"ixx", to_string(self.ixx())},
       std::pair{"stable", to_string(self.stable())},
       std::pair{"live_branch", to_string(self.live_branch())},
       std::pair{"g1_integral", to_string(self.g1_integral())},
       std::pair{"g2_integral", to_string(self.g2_integral())},
       std::pair{"g3_integral", to_string(self.g3_integral())},
       std::pair{"bookkeeping_state", to_string(self.bookkeeping_state())},
       std::pair{"beam_init", to_string(self.beam_init())}});
}
std::string to_string(const ModeInfoProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ModeInfoProxy",
      {std::pair{"stable", to_string(self.stable())},
       std::pair{"tune", to_string(self.tune())},
       std::pair{"emit", to_string(self.emit())},
       std::pair{"chrom", to_string(self.chrom())},
       std::pair{"sigma", to_string(self.sigma())},
       std::pair{"sigmap", to_string(self.sigmap())}});
}
std::string to_string(const PreTrackerProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PreTrackerProxy",
      {std::pair{"who", to_string(self.who())},
       std::pair{"ix_ele_start", to_string(self.ix_ele_start())},
       std::pair{"ix_ele_end", to_string(self.ix_ele_end())},
       std::pair{"input_file", self.input_file()}});
}
std::string to_string(const AnormalModeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "AnormalModeProxy",
      {std::pair{"emittance", to_string(self.emittance())},
       std::pair{"emittance_no_vert", to_string(self.emittance_no_vert())},
       std::pair{"synch_int", to_string(self.synch_int())},
       std::pair{"j_damp", to_string(self.j_damp())},
       std::pair{"alpha_damp", to_string(self.alpha_damp())},
       std::pair{"chrom", to_string(self.chrom())},
       std::pair{"tune", to_string(self.tune())}});
}
std::string to_string(const LinacNormalModeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LinacNormalModeProxy",
      {std::pair{"i2_E4", to_string(self.i2_E4())},
       std::pair{"i3_E7", to_string(self.i3_E7())},
       std::pair{"i5a_E6", to_string(self.i5a_E6())},
       std::pair{"i5b_E6", to_string(self.i5b_E6())},
       std::pair{"sig_E1", to_string(self.sig_E1())},
       std::pair{"a_emittance_end", to_string(self.a_emittance_end())},
       std::pair{"b_emittance_end", to_string(self.b_emittance_end())}});
}
std::string to_string(const NormalModesProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "NormalModesProxy",
      {std::pair{"synch_int", to_string(self.synch_int())},
       std::pair{"sigE_E", to_string(self.sigE_E())},
       std::pair{"sig_z", to_string(self.sig_z())},
       std::pair{"e_loss", to_string(self.e_loss())},
       std::pair{"rf_voltage", to_string(self.rf_voltage())},
       std::pair{"pz_aperture", to_string(self.pz_aperture())},
       std::pair{"pz_average", to_string(self.pz_average())},
       std::pair{"momentum_compaction", to_string(self.momentum_compaction())},
       std::pair{"dpz_damp", to_string(self.dpz_damp())},
       std::pair{"a", to_string(self.a())},
       std::pair{"b", to_string(self.b())},
       std::pair{"z", to_string(self.z())},
       std::pair{"lin", to_string(self.lin())}});
}
std::string to_string(const EmFieldProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "EmFieldProxy",
      {std::pair{"E", to_string(self.E())},
       std::pair{"B", to_string(self.B())},
       std::pair{"dE", to_string(self.dE())},
       std::pair{"dB", to_string(self.dB())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"phi_B", to_string(self.phi_B())},
       std::pair{"A", to_string(self.A())}});
}
std::string to_string(const StrongBeamProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "StrongBeamProxy",
      {std::pair{"ix_slice", to_string(self.ix_slice())},
       std::pair{"x_center", to_string(self.x_center())},
       std::pair{"y_center", to_string(self.y_center())},
       std::pair{"x_sigma", to_string(self.x_sigma())},
       std::pair{"y_sigma", to_string(self.y_sigma())},
       std::pair{"dx", to_string(self.dx())},
       std::pair{"dy", to_string(self.dy())}});
}
std::string to_string(const TrackPointProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TrackPointProxy",
      {std::pair{"s_lab", to_string(self.s_lab())},
       std::pair{"s_body", to_string(self.s_body())},
       std::pair{"orb", to_string(self.orb())},
       std::pair{"field", to_string(self.field())},
       std::pair{"strong_beam", to_string(self.strong_beam())},
       std::pair{"vec0", to_string(self.vec0())},
       std::pair{"mat6", to_string(self.mat6())}});
}
std::string to_string(const TrackProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TrackProxy",
      {std::pair{"pt", "[...]"},
       std::pair{"ds_save", to_string(self.ds_save())},
       std::pair{"n_pt", to_string(self.n_pt())},
       std::pair{"n_bad", to_string(self.n_bad())},
       std::pair{"n_ok", to_string(self.n_ok())}});
}
std::string to_string(const SpaceChargeCommonProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SpaceChargeCommonProxy",
      {std::pair{"ds_track_step", to_string(self.ds_track_step())},
       std::pair{"dt_track_step", to_string(self.dt_track_step())},
       std::pair{
           "cathode_strength_cutoff",
           to_string(self.cathode_strength_cutoff())},
       std::pair{"rel_tol_tracking", to_string(self.rel_tol_tracking())},
       std::pair{"abs_tol_tracking", to_string(self.abs_tol_tracking())},
       std::pair{"beam_chamber_height", to_string(self.beam_chamber_height())},
       std::pair{"lsc_sigma_cutoff", to_string(self.lsc_sigma_cutoff())},
       std::pair{
           "particle_sigma_cutoff", to_string(self.particle_sigma_cutoff())},
       std::pair{
           "space_charge_mesh_size", to_string(self.space_charge_mesh_size())},
       std::pair{"csr3d_mesh_size", to_string(self.csr3d_mesh_size())},
       std::pair{"n_bin", to_string(self.n_bin())},
       std::pair{"particle_bin_span", to_string(self.particle_bin_span())},
       std::pair{"n_shield_images", to_string(self.n_shield_images())},
       std::pair{"sc_min_in_bin", to_string(self.sc_min_in_bin())},
       std::pair{
           "lsc_kick_transverse_dependence",
           to_string(self.lsc_kick_transverse_dependence())},
       std::pair{"debug", to_string(self.debug())},
       std::pair{"diagnostic_output_file", self.diagnostic_output_file()}});
}
std::string to_string(const BmadCommonProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BmadCommonProxy",
      {std::pair{"max_aperture_limit", to_string(self.max_aperture_limit())},
       std::pair{"d_orb", to_string(self.d_orb())},
       std::pair{"default_ds_step", to_string(self.default_ds_step())},
       std::pair{"significant_length", to_string(self.significant_length())},
       std::pair{"rel_tol_tracking", to_string(self.rel_tol_tracking())},
       std::pair{"abs_tol_tracking", to_string(self.abs_tol_tracking())},
       std::pair{
           "rel_tol_adaptive_tracking",
           to_string(self.rel_tol_adaptive_tracking())},
       std::pair{
           "abs_tol_adaptive_tracking",
           to_string(self.abs_tol_adaptive_tracking())},
       std::pair{
           "init_ds_adaptive_tracking",
           to_string(self.init_ds_adaptive_tracking())},
       std::pair{
           "min_ds_adaptive_tracking",
           to_string(self.min_ds_adaptive_tracking())},
       std::pair{
           "fatal_ds_adaptive_tracking",
           to_string(self.fatal_ds_adaptive_tracking())},
       std::pair{
           "autoscale_amp_abs_tol", to_string(self.autoscale_amp_abs_tol())},
       std::pair{
           "autoscale_amp_rel_tol", to_string(self.autoscale_amp_rel_tol())},
       std::pair{"autoscale_phase_tol", to_string(self.autoscale_phase_tol())},
       std::pair{
           "electric_dipole_moment", to_string(self.electric_dipole_moment())},
       std::pair{"synch_rad_scale", to_string(self.synch_rad_scale())},
       std::pair{"sad_eps_scale", to_string(self.sad_eps_scale())},
       std::pair{"sad_amp_max", to_string(self.sad_amp_max())},
       std::pair{"sad_n_div_max", to_string(self.sad_n_div_max())},
       std::pair{"taylor_order", to_string(self.taylor_order())},
       std::pair{"runge_kutta_order", to_string(self.runge_kutta_order())},
       std::pair{"default_integ_order", to_string(self.default_integ_order())},
       std::pair{
           "max_num_runge_kutta_step",
           to_string(self.max_num_runge_kutta_step())},
       std::pair{
           "rf_phase_below_transition_ref",
           to_string(self.rf_phase_below_transition_ref())},
       std::pair{"sr_wakes_on", to_string(self.sr_wakes_on())},
       std::pair{"lr_wakes_on", to_string(self.lr_wakes_on())},
       std::pair{"auto_bookkeeper", to_string(self.auto_bookkeeper())},
       std::pair{
           "high_energy_space_charge_on",
           to_string(self.high_energy_space_charge_on())},
       std::pair{
           "csr_and_space_charge_on",
           to_string(self.csr_and_space_charge_on())},
       std::pair{"spin_tracking_on", to_string(self.spin_tracking_on())},
       std::pair{
           "spin_sokolov_ternov_flipping_on",
           to_string(self.spin_sokolov_ternov_flipping_on())},
       std::pair{
           "radiation_damping_on", to_string(self.radiation_damping_on())},
       std::pair{
           "radiation_zero_average", to_string(self.radiation_zero_average())},
       std::pair{
           "radiation_fluctuations_on",
           to_string(self.radiation_fluctuations_on())},
       std::pair{
           "conserve_taylor_maps", to_string(self.conserve_taylor_maps())},
       std::pair{
           "absolute_time_tracking", to_string(self.absolute_time_tracking())},
       std::pair{
           "absolute_time_ref_shift",
           to_string(self.absolute_time_ref_shift())},
       std::pair{
           "convert_to_kinetic_momentum",
           to_string(self.convert_to_kinetic_momentum())},
       std::pair{"normalize_twiss", to_string(self.normalize_twiss())},
       std::pair{"aperture_limit_on", to_string(self.aperture_limit_on())},
       std::pair{
           "spin_n0_direction_user_set",
           to_string(self.spin_n0_direction_user_set())},
       std::pair{"debug", to_string(self.debug())}});
}
std::string to_string(const RadInt1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RadInt1Proxy",
      {std::pair{"i0", to_string(self.i0())},
       std::pair{"i1", to_string(self.i1())},
       std::pair{"i2", to_string(self.i2())},
       std::pair{"i3", to_string(self.i3())},
       std::pair{"i4a", to_string(self.i4a())},
       std::pair{"i4b", to_string(self.i4b())},
       std::pair{"i4z", to_string(self.i4z())},
       std::pair{"i5a", to_string(self.i5a())},
       std::pair{"i5b", to_string(self.i5b())},
       std::pair{"i6b", to_string(self.i6b())},
       std::pair{"lin_i2_E4", to_string(self.lin_i2_E4())},
       std::pair{"lin_i3_E7", to_string(self.lin_i3_E7())},
       std::pair{"lin_i5a_E6", to_string(self.lin_i5a_E6())},
       std::pair{"lin_i5b_E6", to_string(self.lin_i5b_E6())},
       std::pair{"lin_norm_emit_a", to_string(self.lin_norm_emit_a())},
       std::pair{"lin_norm_emit_b", to_string(self.lin_norm_emit_b())},
       std::pair{"lin_sig_E", to_string(self.lin_sig_E())},
       std::pair{"n_steps", to_string(self.n_steps())}});
}
std::string to_string(const RadIntBranchProxy& self) {
  return repr(
      self.get_fortran_ptr(), "RadIntBranchProxy", {std::pair{"ele", "[...]"}});
}
std::string to_string(const RadIntAllEleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RadIntAllEleProxy",
      {std::pair{"branch", "[...]"}});
}
std::string to_string(const RfStairStepProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RfStairStepProxy",
      {std::pair{"E_tot0", to_string(self.E_tot0())},
       std::pair{"E_tot1", to_string(self.E_tot1())},
       std::pair{"p0c", to_string(self.p0c())},
       std::pair{"p1c", to_string(self.p1c())},
       std::pair{"scale", to_string(self.scale())},
       std::pair{"time", to_string(self.time())},
       std::pair{"s0", to_string(self.s0())},
       std::pair{"s", to_string(self.s())},
       std::pair{"ix_step", to_string(self.ix_step())}});
}
std::string to_string(const RfEleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "RfEleProxy",
      {std::pair{"steps", "[...]"},
       std::pair{"ds_step", to_string(self.ds_step())}});
}
std::string to_string(const EleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "EleProxy",
      {std::pair{"name", self.name()},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_ele", to_string(self.ix_ele())}});
}
std::string to_string(const ComplexTaylorTermProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ComplexTaylorTermProxy",
      {std::pair{"coef", to_string(self.coef())},
       std::pair{"expn", to_string(self.expn())}});
}
std::string to_string(const ComplexTaylorProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ComplexTaylorProxy",
      {std::pair{"ref", to_string(self.ref())}, std::pair{"term", "[...]"}});
}
std::string to_string(const BranchProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BranchProxy",
      {std::pair{"name", self.name()},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_from_branch", to_string(self.ix_from_branch())},
       std::pair{"ix_from_ele", to_string(self.ix_from_ele())},
       std::pair{"ix_to_ele", to_string(self.ix_to_ele())},
       std::pair{"ix_fixer", to_string(self.ix_fixer())},
       std::pair{"n_ele_track", to_string(self.n_ele_track())},
       std::pair{"n_ele_max", to_string(self.n_ele_max())},
       std::pair{"lat", to_string(self.lat())},
       std::pair{"a", to_string(self.a())},
       std::pair{"b", to_string(self.b())},
       std::pair{"z", to_string(self.z())},
       std::pair{"ele", "[...]"},
       std::pair{"param", to_string(self.param())},
       std::pair{"particle_start", to_string(self.particle_start())},
       std::pair{"wall3d", "[...]"}});
}
std::string to_string(const LatProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LatProxy",
      {std::pair{"use_name", self.use_name()},
       std::pair{"#branch", to_string(self.branch().size())}});
}
std::string to_string(const BunchProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BunchProxy",
      {std::pair{"particle", "[...]"},
       std::pair{"ix_z", to_string(self.ix_z())},
       std::pair{"charge_tot", to_string(self.charge_tot())},
       std::pair{"charge_live", to_string(self.charge_live())},
       std::pair{"z_center", to_string(self.z_center())},
       std::pair{"t_center", to_string(self.t_center())},
       std::pair{"t0", to_string(self.t0())},
       std::pair{
           "drift_between_t_and_s", to_string(self.drift_between_t_and_s())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_bunch", to_string(self.ix_bunch())},
       std::pair{"ix_turn", to_string(self.ix_turn())},
       std::pair{"n_live", to_string(self.n_live())},
       std::pair{"n_good", to_string(self.n_good())},
       std::pair{"n_bad", to_string(self.n_bad())}});
}
std::string to_string(const BunchParamsProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BunchParamsProxy",
      {std::pair{"centroid", to_string(self.centroid())},
       std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"z", to_string(self.z())},
       std::pair{"a", to_string(self.a())},
       std::pair{"b", to_string(self.b())},
       std::pair{"c", to_string(self.c())},
       std::pair{"sigma", to_string(self.sigma())},
       std::pair{"rel_max", to_string(self.rel_max())},
       std::pair{"rel_min", to_string(self.rel_min())},
       std::pair{"s", to_string(self.s())},
       std::pair{"t", to_string(self.t())},
       std::pair{"sigma_t", to_string(self.sigma_t())},
       std::pair{"charge_live", to_string(self.charge_live())},
       std::pair{"charge_tot", to_string(self.charge_tot())},
       std::pair{"n_particle_tot", to_string(self.n_particle_tot())},
       std::pair{"n_particle_live", to_string(self.n_particle_live())},
       std::pair{
           "n_particle_lost_in_ele", to_string(self.n_particle_lost_in_ele())},
       std::pair{"n_good_steps", to_string(self.n_good_steps())},
       std::pair{"n_bad_steps", to_string(self.n_bad_steps())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"location", to_string(self.location())},
       std::pair{"twiss_valid", to_string(self.twiss_valid())}});
}
std::string to_string(const BeamProxy& self) {
  return repr(
      self.get_fortran_ptr(), "BeamProxy", {std::pair{"bunch", "[...]"}});
}
std::string to_string(const AperturePointProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "AperturePointProxy",
      {std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"plane", to_string(self.plane())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"i_turn", to_string(self.i_turn())}});
}
std::string to_string(const ApertureParamProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ApertureParamProxy",
      {std::pair{"min_angle", to_string(self.min_angle())},
       std::pair{"max_angle", to_string(self.max_angle())},
       std::pair{"n_angle", to_string(self.n_angle())},
       std::pair{"n_turn", to_string(self.n_turn())},
       std::pair{"x_init", to_string(self.x_init())},
       std::pair{"y_init", to_string(self.y_init())},
       std::pair{"rel_accuracy", to_string(self.rel_accuracy())},
       std::pair{"abs_accuracy", to_string(self.abs_accuracy())},
       std::pair{"start_ele", self.start_ele()}});
}
std::string to_string(const ApertureScanProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ApertureScanProxy",
      {std::pair{"point", "[...]"},
       std::pair{"ref_orb", to_string(self.ref_orb())},
       std::pair{"pz_start", to_string(self.pz_start())}});
}
std::string to_string(const ElePointerProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ElePointerProxy",
      {std::pair{"ele", to_string(self.ele())},
       std::pair{"loc", to_string(self.loc())},
       std::pair{"id", to_string(self.id())}});
}
std::string to_string(const ExpressionTreeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ExpressionTreeProxy",
      {std::pair{"name", self.name()},
       std::pair{"type", to_string(self.type())},
       std::pair{"value", to_string(self.value())},
       std::pair{"node", "[...]"}});
}
std::string to_string(const NametableProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "NametableProxy",
      {std::pair{"name", to_string(self.name())},
       std::pair{"index", to_string(self.index())},
       std::pair{"n_min", to_string(self.n_min())},
       std::pair{"n_max", to_string(self.n_max())}});
}
std::string to_string(const TaoSpinDnDpzProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoSpinDnDpzProxy",
      {std::pair{"vec", to_string(self.vec())},
       std::pair{"partial", to_string(self.partial())},
       std::pair{"partial2", to_string(self.partial2())}});
}
std::string to_string(const ResonanceHProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "ResonanceHProxy",
      {std::pair{"id", self.id()},
       std::pair{"c_val", to_string(self.c_val())}});
}
std::string to_string(const SpinOrbitMap1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SpinOrbitMap1Proxy",
      {std::pair{"orb_mat", to_string(self.orb_mat())},
       std::pair{"vec0", to_string(self.vec0())},
       std::pair{"spin_q", to_string(self.spin_q())}});
}
std::string to_string(const SpinAxisProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SpinAxisProxy",
      {std::pair{"l", to_string(self.l())},
       std::pair{"n0", to_string(self.n0())},
       std::pair{"m", to_string(self.m())}});
}
std::string to_string(const PtcNormalFormProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "PtcNormalFormProxy",
      {std::pair{"ele_origin", to_string(self.ele_origin())},
       std::pair{"orb0", to_string(self.orb0())},
       std::pair{"valid_map", to_string(self.valid_map())}});
}
std::string to_string(const BmadNormalFormProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BmadNormalFormProxy",
      {std::pair{"ele_origin", to_string(self.ele_origin())},
       std::pair{"M", "[...]"},
       std::pair{"A", "[...]"},
       std::pair{"A_inv", "[...]"},
       std::pair{"dhdj", "[...]"},
       std::pair{"F", "[...]"},
       std::pair{"L", "[...]"},
       std::pair{"h", "[...]"}});
}
std::string to_string(const BunchTrackProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BunchTrackProxy",
      {std::pair{"pt", "[...]"},
       std::pair{"ds_save", to_string(self.ds_save())},
       std::pair{"n_pt", to_string(self.n_pt())}});
}
std::string to_string(const SummationRdtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "SummationRdtProxy",
      {std::pair{"h11001", to_string(self.h11001())},
       std::pair{"h00111", to_string(self.h00111())},
       std::pair{"h20001", to_string(self.h20001())},
       std::pair{"h00201", to_string(self.h00201())},
       std::pair{"h10002", to_string(self.h10002())},
       std::pair{"h21000", to_string(self.h21000())},
       std::pair{"h30000", to_string(self.h30000())},
       std::pair{"h10110", to_string(self.h10110())},
       std::pair{"h10020", to_string(self.h10020())},
       std::pair{"h10200", to_string(self.h10200())},
       std::pair{"h31000", to_string(self.h31000())},
       std::pair{"h40000", to_string(self.h40000())},
       std::pair{"h20110", to_string(self.h20110())},
       std::pair{"h11200", to_string(self.h11200())},
       std::pair{"h20020", to_string(self.h20020())},
       std::pair{"h20200", to_string(self.h20200())},
       std::pair{"h00310", to_string(self.h00310())},
       std::pair{"h00400", to_string(self.h00400())},
       std::pair{"h22000", to_string(self.h22000())},
       std::pair{"h00220", to_string(self.h00220())},
       std::pair{"h11110", to_string(self.h11110())}});
}
std::string to_string(const TaoEleShapeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoEleShapeProxy",
      {std::pair{"ele_id", self.ele_id()},
       std::pair{"shape", self.shape()},
       std::pair{"color", self.color()},
       std::pair{"size", to_string(self.size())},
       std::pair{"label", self.label()},
       std::pair{"draw", to_string(self.draw())},
       std::pair{"multi", to_string(self.multi())},
       std::pair{"line_width", to_string(self.line_width())},
       std::pair{"offset", to_string(self.offset())},
       std::pair{"ix_key", to_string(self.ix_key())},
       std::pair{"name_ele", self.name_ele()},
       std::pair{"uni", "[...]"}});
}
std::string to_string(const TaoElePointerProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoElePointerProxy",
      {std::pair{"eles", "[...]"},
       std::pair{"n_loc", to_string(self.n_loc())}});
}
std::string to_string(const TaoCurveProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoCurveProxy",
      {std::pair{"name", self.name()},
       std::pair{"data_source", self.data_source()},
       std::pair{"data_index", self.data_index()},
       std::pair{"data_type_x", self.data_type_x()},
       std::pair{"data_type", self.data_type()},
       std::pair{"ele_ref_name", self.ele_ref_name()},
       std::pair{"legend_text", self.legend_text()},
       std::pair{"message_text", self.message_text()},
       std::pair{"component", self.component()},
       std::pair{"why_invalid", self.why_invalid()},
       std::pair{"g", "..."},
       std::pair{"hist", to_string(self.hist())},
       std::pair{"z_color", to_string(self.z_color())},
       std::pair{"x_line", to_string(self.x_line())},
       std::pair{"y_line", to_string(self.y_line())},
       std::pair{"y2_line", to_string(self.y2_line())},
       std::pair{"ix_line", to_string(self.ix_line())},
       std::pair{"x_symb", to_string(self.x_symb())},
       std::pair{"y_symb", to_string(self.y_symb())},
       std::pair{"z_symb", to_string(self.z_symb())},
       std::pair{"err_symb", to_string(self.err_symb())},
       std::pair{"symb_size", to_string(self.symb_size())},
       std::pair{"ix_symb", to_string(self.ix_symb())},
       std::pair{"y_axis_scale_factor", to_string(self.y_axis_scale_factor())},
       std::pair{"line", to_string(self.line())},
       std::pair{"symbol", to_string(self.symbol())},
       std::pair{"orbit", to_string(self.orbit())},
       std::pair{"ix_universe", to_string(self.ix_universe())},
       std::pair{"symbol_every", to_string(self.symbol_every())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_bunch", to_string(self.ix_bunch())},
       std::pair{"n_turn", to_string(self.n_turn())},
       std::pair{"use_y2", to_string(self.use_y2())},
       std::pair{"draw_line", to_string(self.draw_line())},
       std::pair{"draw_symbols", to_string(self.draw_symbols())},
       std::pair{"draw_symbol_index", to_string(self.draw_symbol_index())},
       std::pair{"draw_error_bars", to_string(self.draw_error_bars())},
       std::pair{"smooth_line_calc", to_string(self.smooth_line_calc())},
       std::pair{"valid", to_string(self.valid())}});
}
std::string to_string(const TaoCurveColorProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoCurveColorProxy",
      {std::pair{"data_type", self.data_type()},
       std::pair{"is_on", to_string(self.is_on())},
       std::pair{"min", to_string(self.min())},
       std::pair{"max", to_string(self.max())},
       std::pair{"autoscale", to_string(self.autoscale())}});
}
std::string to_string(const TaoCurveOrbitProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoCurveOrbitProxy",
      {std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"t", to_string(self.t())}});
}
std::string to_string(const TaoHistogramProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoHistogramProxy",
      {std::pair{"density_normalized", to_string(self.density_normalized())},
       std::pair{"weight_by_charge", to_string(self.weight_by_charge())},
       std::pair{"minimum", to_string(self.minimum())},
       std::pair{"maximum", to_string(self.maximum())},
       std::pair{"width", to_string(self.width())},
       std::pair{"center", to_string(self.center())},
       std::pair{"number", to_string(self.number())}});
}
std::string to_string(const LatEleOrder1Proxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LatEleOrder1Proxy",
      {std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_order", to_string(self.ix_order())}});
}
std::string to_string(const LatEleOrderArrayProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LatEleOrderArrayProxy",
      {std::pair{"ele", "[...]"}});
}
std::string to_string(const TaoLatSigmaProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoLatSigmaProxy",
      {std::pair{"mat", to_string(self.mat())}});
}
std::string to_string(const TaoSpinEleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoSpinEleProxy",
      {std::pair{"dn_dpz", to_string(self.dn_dpz())},
       std::pair{"orb_eigen_val", to_string(self.orb_eigen_val())},
       std::pair{"orb_eigen_vec", to_string(self.orb_eigen_vec())},
       std::pair{"spin_eigen_vec", to_string(self.spin_eigen_vec())},
       std::pair{"valid", to_string(self.valid())}});
}
std::string to_string(const TaoPlotCacheProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoPlotCacheProxy",
      {std::pair{"ele_to_s", to_string(self.ele_to_s())},
       std::pair{"orbit", to_string(self.orbit())},
       std::pair{"err", to_string(self.err())}});
}
std::string to_string(const TaoSpinPolarizationProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoSpinPolarizationProxy",
      {std::pair{"tune", to_string(self.tune())},
       std::pair{"pol_limit_st", to_string(self.pol_limit_st())},
       std::pair{"pol_limit_dk", to_string(self.pol_limit_dk())},
       std::pair{
           "pol_limit_dk_partial", to_string(self.pol_limit_dk_partial())},
       std::pair{
           "pol_limit_dk_partial2", to_string(self.pol_limit_dk_partial2())},
       std::pair{"pol_rate_bks", to_string(self.pol_rate_bks())},
       std::pair{"depol_rate", to_string(self.depol_rate())},
       std::pair{"depol_rate_partial", to_string(self.depol_rate_partial())},
       std::pair{"depol_rate_partial2", to_string(self.depol_rate_partial2())},
       std::pair{"integral_bn", to_string(self.integral_bn())},
       std::pair{"integral_bdn", to_string(self.integral_bdn())},
       std::pair{"integral_1ns", to_string(self.integral_1ns())},
       std::pair{"integral_dn2", to_string(self.integral_dn2())},
       std::pair{"valid", to_string(self.valid())},
       std::pair{"q_1turn", to_string(self.q_1turn())},
       std::pair{"q_ele", "[...]"}});
}
std::string to_string(const TaoLatticeBranchProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoLatticeBranchProxy",
      {std::pair{"tao_lat", "..."},
       std::pair{"lat_sigma", "[...]"},
       std::pair{"spin_ele", "[...]"},
       std::pair{"bunch_params", "[...]"},
       std::pair{"bunch_params_comb", "[...]"},
       std::pair{"orbit", "[...]"},
       std::pair{"plot_cache", "[...]"},
       std::pair{"spin", to_string(self.spin())},
       std::pair{"srdt", to_string(self.srdt())},
       std::pair{"orb0", to_string(self.orb0())},
       std::pair{"modes_ri", to_string(self.modes_ri())},
       std::pair{"modes_6d", to_string(self.modes_6d())},
       std::pair{"ptc_normal_form", to_string(self.ptc_normal_form())},
       std::pair{"bmad_normal_form", to_string(self.bmad_normal_form())},
       std::pair{"high_E_orb", "[...]"},
       std::pair{"low_E_orb", "[...]"},
       std::pair{"taylor_save", "[...]"},
       std::pair{"cache_x_min", to_string(self.cache_x_min())},
       std::pair{"cache_x_max", to_string(self.cache_x_max())},
       std::pair{"comb_ds_save", to_string(self.comb_ds_save())},
       std::pair{"ix_ref_taylor", to_string(self.ix_ref_taylor())},
       std::pair{"ix_ele_taylor", to_string(self.ix_ele_taylor())},
       std::pair{"track_state", to_string(self.track_state())},
       std::pair{"cache_n_pts", to_string(self.cache_n_pts())},
       std::pair{"ix_rad_int_cache", to_string(self.ix_rad_int_cache())},
       std::pair{
           "has_open_match_element", to_string(self.has_open_match_element())},
       std::pair{"plot_cache_valid", to_string(self.plot_cache_valid())},
       std::pair{"spin_map_valid", to_string(self.spin_map_valid())},
       std::pair{"twiss_valid", to_string(self.twiss_valid())},
       std::pair{"mode_flip_here", to_string(self.mode_flip_here())},
       std::pair{"chrom_calc_ok", to_string(self.chrom_calc_ok())},
       std::pair{"rad_int_calc_ok", to_string(self.rad_int_calc_ok())},
       std::pair{"emit_6d_calc_ok", to_string(self.emit_6d_calc_ok())},
       std::pair{"sigma_track_ok", to_string(self.sigma_track_ok())}});
}
std::string to_string(const TaoModelElementProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoModelElementProxy",
      {std::pair{"beam", to_string(self.beam())},
       std::pair{
           "save_beam_internally", to_string(self.save_beam_internally())},
       std::pair{"save_beam_to_file", to_string(self.save_beam_to_file())}});
}
std::string to_string(const TaoBeamBranchProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoBeamBranchProxy",
      {std::pair{"beam_at_start", to_string(self.beam_at_start())},
       std::pair{"beam_init", to_string(self.beam_init())},
       std::pair{"beam_init_used", to_string(self.beam_init_used())},
       std::pair{
           "init_starting_distribution",
           to_string(self.init_starting_distribution())},
       std::pair{"track_start", self.track_start()},
       std::pair{"track_end", self.track_end()},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_track_start", to_string(self.ix_track_start())},
       std::pair{"ix_track_end", to_string(self.ix_track_end())}});
}
std::string to_string(const TaoD1DataProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoD1DataProxy",
      {std::pair{"name", self.name()},
       std::pair{"d2", "..."},
       std::pair{"d", "[...]"}});
}
std::string to_string(const TaoD2DataProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoD2DataProxy",
      {std::pair{"name", self.name()},
       std::pair{"data_file_name", self.data_file_name()},
       std::pair{"ref_file_name", self.ref_file_name()},
       std::pair{"data_date", self.data_date()},
       std::pair{"ref_date", self.ref_date()},
       std::pair{"descrip", to_string(self.descrip())},
       std::pair{"d1", "[...]"},
       std::pair{"ix_universe", to_string(self.ix_universe())},
       std::pair{"ix_d2_data", to_string(self.ix_d2_data())},
       std::pair{"ix_ref", to_string(self.ix_ref())},
       std::pair{"data_read_in", to_string(self.data_read_in())},
       std::pair{"ref_read_in", to_string(self.ref_read_in())}});
}
std::string to_string(const TaoDataVarComponentProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoDataVarComponentProxy",
      {std::pair{"name", self.name()},
       std::pair{"sign", to_string(self.sign())}});
}
std::string to_string(const TaoGraphProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoGraphProxy",
      {std::pair{"name", self.name()},
       std::pair{"type", self.type()},
       std::pair{"title", self.title()},
       std::pair{"title_suffix", self.title_suffix()},
       std::pair{"text_legend", to_string(self.text_legend())},
       std::pair{"text_legend_out", to_string(self.text_legend_out())},
       std::pair{"why_invalid", self.why_invalid()},
       std::pair{"curve", "[...]"},
       std::pair{"p", "..."},
       std::pair{"floor_plan", to_string(self.floor_plan())},
       std::pair{"text_legend_origin", to_string(self.text_legend_origin())},
       std::pair{"curve_legend_origin", to_string(self.curve_legend_origin())},
       std::pair{"curve_legend", to_string(self.curve_legend())},
       std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"x2", to_string(self.x2())},
       std::pair{"y2", to_string(self.y2())},
       std::pair{"margin", to_string(self.margin())},
       std::pair{"scale_margin", to_string(self.scale_margin())},
       std::pair{"x_axis_scale_factor", to_string(self.x_axis_scale_factor())},
       std::pair{"symbol_size_scale", to_string(self.symbol_size_scale())},
       std::pair{"box", to_string(self.box())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_universe", to_string(self.ix_universe())},
       std::pair{"clip", to_string(self.clip())},
       std::pair{"y2_mirrors_y", to_string(self.y2_mirrors_y())},
       std::pair{"limited", to_string(self.limited())},
       std::pair{"draw_axes", to_string(self.draw_axes())},
       std::pair{"draw_curve_legend", to_string(self.draw_curve_legend())},
       std::pair{"draw_grid", to_string(self.draw_grid())},
       std::pair{"draw_title", to_string(self.draw_title())},
       std::pair{
           "draw_only_good_user_data_or_vars",
           to_string(self.draw_only_good_user_data_or_vars())},
       std::pair{"allow_wrap_around", to_string(self.allow_wrap_around())},
       std::pair{"is_valid", to_string(self.is_valid())}});
}
std::string to_string(const TaoPlotProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoPlotProxy",
      {std::pair{"name", self.name()},
       std::pair{"description", self.description()},
       std::pair{"graph", "[...]"},
       std::pair{"r", "..."},
       std::pair{"ix_plot", to_string(self.ix_plot())},
       std::pair{"n_curve_pts", to_string(self.n_curve_pts())},
       std::pair{"type", self.type()},
       std::pair{"x_axis_type", self.x_axis_type()},
       std::pair{"autoscale_x", to_string(self.autoscale_x())},
       std::pair{"autoscale_y", to_string(self.autoscale_y())},
       std::pair{"autoscale_gang_x", to_string(self.autoscale_gang_x())},
       std::pair{"autoscale_gang_y", to_string(self.autoscale_gang_y())},
       std::pair{
           "list_with_show_plot_command",
           to_string(self.list_with_show_plot_command())},
       std::pair{"phantom", to_string(self.phantom())},
       std::pair{"default_plot", to_string(self.default_plot())}});
}
std::string to_string(const TaoPlotRegionProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoPlotRegionProxy",
      {std::pair{"name", self.name()},
       std::pair{"plot", to_string(self.plot())},
       std::pair{"location", to_string(self.location())},
       std::pair{"visible", to_string(self.visible())},
       std::pair{
           "list_with_show_plot_command",
           to_string(self.list_with_show_plot_command())},
       std::pair{"setup_done", to_string(self.setup_done())}});
}
std::string to_string(const TaoUniversePointerProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoUniversePointerProxy",
      {std::pair{"u", to_string(self.u())}});
}
std::string to_string(const TaoSuperUniverseProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoSuperUniverseProxy",
      {std::pair{"global", to_string(self.global())},
       std::pair{"init", to_string(self.init())},
       std::pair{"com", to_string(self.com())},
       std::pair{"plot_page", to_string(self.plot_page())},
       std::pair{"v1_var", "[...]"},
       std::pair{"var", "[...]"},
       std::pair{"u", "[...]"},
       std::pair{"key", to_string(self.key())},
       std::pair{"building_wall", to_string(self.building_wall())},
       std::pair{"wave", to_string(self.wave())},
       std::pair{"n_var_used", to_string(self.n_var_used())},
       std::pair{"n_v1_var_used", to_string(self.n_v1_var_used())},
       std::pair{"history", "[...]"},
       std::pair{"initialized", to_string(self.initialized())}});
}
std::string to_string(const TaoVarProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoVarProxy",
      {std::pair{"ele_name", self.ele_name()},
       std::pair{"attrib_name", self.attrib_name()},
       std::pair{"id", self.id()},
       std::pair{"slave", "[...]"},
       std::pair{"ix_v1", to_string(self.ix_v1())},
       std::pair{"ix_var", to_string(self.ix_var())},
       std::pair{"ix_dvar", to_string(self.ix_dvar())},
       std::pair{"ix_attrib", to_string(self.ix_attrib())},
       std::pair{"ix_key_table", to_string(self.ix_key_table())},
       std::pair{"model_value", to_string(self.model_value())},
       std::pair{"base_value", to_string(self.base_value())},
       std::pair{"design_value", to_string(self.design_value())},
       std::pair{"scratch_value", to_string(self.scratch_value())},
       std::pair{"old_value", to_string(self.old_value())},
       std::pair{"meas_value", to_string(self.meas_value())},
       std::pair{"ref_value", to_string(self.ref_value())},
       std::pair{"correction_value", to_string(self.correction_value())},
       std::pair{"high_lim", to_string(self.high_lim())},
       std::pair{"low_lim", to_string(self.low_lim())},
       std::pair{"step", to_string(self.step())},
       std::pair{"weight", to_string(self.weight())},
       std::pair{"delta_merit", to_string(self.delta_merit())},
       std::pair{"merit", to_string(self.merit())},
       std::pair{"dMerit_dVar", to_string(self.dMerit_dVar())},
       std::pair{"key_val0", to_string(self.key_val0())},
       std::pair{"key_delta", to_string(self.key_delta())},
       std::pair{"s", to_string(self.s())},
       std::pair{"extend_val", to_string(self.extend_val())},
       std::pair{"merit_type", self.merit_type()},
       std::pair{"exists", to_string(self.exists())},
       std::pair{"good_var", to_string(self.good_var())},
       std::pair{"good_user", to_string(self.good_user())},
       std::pair{"good_opt", to_string(self.good_opt())},
       std::pair{"good_plot", to_string(self.good_plot())},
       std::pair{"useit_opt", to_string(self.useit_opt())},
       std::pair{"useit_plot", to_string(self.useit_plot())},
       std::pair{"key_bound", to_string(self.key_bound())},
       std::pair{"v1", "..."}});
}
std::string to_string(const TaoVarSlaveProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoVarSlaveProxy",
      {std::pair{"ix_uni", to_string(self.ix_uni())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"model_value", to_string(self.model_value())},
       std::pair{"base_value", to_string(self.base_value())}});
}
std::string to_string(const TaoLatticeProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoLatticeProxy",
      {std::pair{"name", self.name()},
       std::pair{"lat", to_string(self.lat())},
       std::pair{"high_E_lat", to_string(self.high_E_lat())},
       std::pair{"low_E_lat", to_string(self.low_E_lat())},
       std::pair{"rad_int_by_ele_ri", to_string(self.rad_int_by_ele_ri())},
       std::pair{"rad_int_by_ele_6d", to_string(self.rad_int_by_ele_6d())},
       std::pair{"tao_branch", "[...]"}});
}
std::string to_string(const TaoBeamUniProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoBeamUniProxy",
      {std::pair{"saved_at", self.saved_at()},
       std::pair{"dump_file", self.dump_file()},
       std::pair{"dump_at", self.dump_at()},
       std::pair{
           "track_beam_in_universe", to_string(self.track_beam_in_universe())},
       std::pair{"always_reinit", to_string(self.always_reinit())}});
}
std::string to_string(const TaoDynamicApertureProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoDynamicApertureProxy",
      {std::pair{"param", to_string(self.param())},
       std::pair{"scan", "[...]"},
       std::pair{"pz", to_string(self.pz())},
       std::pair{"ellipse_scale", to_string(self.ellipse_scale())},
       std::pair{"a_emit", to_string(self.a_emit())},
       std::pair{"b_emit", to_string(self.b_emit())}});
}
std::string to_string(const TaoModelBranchProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoModelBranchProxy",
      {std::pair{"ele", "[...]"}, std::pair{"beam", to_string(self.beam())}});
}
std::string to_string(const TaoSpinMapProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoSpinMapProxy",
      {std::pair{"valid", to_string(self.valid())},
       std::pair{"map1", to_string(self.map1())},
       std::pair{"axis_input", to_string(self.axis_input())},
       std::pair{"axis0", to_string(self.axis0())},
       std::pair{"axis1", to_string(self.axis1())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_ref", to_string(self.ix_ref())},
       std::pair{"ix_uni", to_string(self.ix_uni())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"mat8", to_string(self.mat8())}});
}
std::string to_string(const TaoDataProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoDataProxy",
      {std::pair{"ele_name", self.ele_name()},
       std::pair{"ele_start_name", self.ele_start_name()},
       std::pair{"ele_ref_name", self.ele_ref_name()},
       std::pair{"data_type", self.data_type()},
       std::pair{"merit_type", self.merit_type()},
       std::pair{"id", self.id()},
       std::pair{"data_source", self.data_source()},
       std::pair{"why_invalid", self.why_invalid()},
       std::pair{"ix_uni", to_string(self.ix_uni())},
       std::pair{"ix_bunch", to_string(self.ix_bunch())},
       std::pair{"ix_branch", to_string(self.ix_branch())},
       std::pair{"ix_ele", to_string(self.ix_ele())},
       std::pair{"ix_ele_start", to_string(self.ix_ele_start())},
       std::pair{"ix_ele_ref", to_string(self.ix_ele_ref())},
       std::pair{"ix_ele_merit", to_string(self.ix_ele_merit())},
       std::pair{"ix_d1", to_string(self.ix_d1())},
       std::pair{"ix_data", to_string(self.ix_data())},
       std::pair{"ix_dModel", to_string(self.ix_dModel())},
       std::pair{"eval_point", to_string(self.eval_point())},
       std::pair{"meas_value", to_string(self.meas_value())},
       std::pair{"ref_value", to_string(self.ref_value())},
       std::pair{"model_value", to_string(self.model_value())},
       std::pair{"design_value", to_string(self.design_value())},
       std::pair{"old_value", to_string(self.old_value())},
       std::pair{"base_value", to_string(self.base_value())},
       std::pair{"error_rms", to_string(self.error_rms())},
       std::pair{"delta_merit", to_string(self.delta_merit())},
       std::pair{"weight", to_string(self.weight())},
       std::pair{"invalid_value", to_string(self.invalid_value())},
       std::pair{"merit", to_string(self.merit())},
       std::pair{"s", to_string(self.s())},
       std::pair{"s_offset", to_string(self.s_offset())},
       std::pair{"ref_s_offset", to_string(self.ref_s_offset())},
       std::pair{"err_message_printed", to_string(self.err_message_printed())},
       std::pair{"exists", to_string(self.exists())},
       std::pair{"good_model", to_string(self.good_model())},
       std::pair{"good_base", to_string(self.good_base())},
       std::pair{"good_design", to_string(self.good_design())},
       std::pair{"good_meas", to_string(self.good_meas())},
       std::pair{"good_ref", to_string(self.good_ref())},
       std::pair{"good_user", to_string(self.good_user())},
       std::pair{"good_opt", to_string(self.good_opt())},
       std::pair{"good_plot", to_string(self.good_plot())},
       std::pair{"useit_plot", to_string(self.useit_plot())},
       std::pair{"useit_opt", to_string(self.useit_opt())},
       std::pair{"spin_map", to_string(self.spin_map())},
       std::pair{"d1", "..."}});
}
std::string to_string(const TaoPingScaleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoPingScaleProxy",
      {std::pair{"a_mode_meas", to_string(self.a_mode_meas())},
       std::pair{"a_mode_ref", to_string(self.a_mode_ref())},
       std::pair{"b_mode_meas", to_string(self.b_mode_meas())},
       std::pair{"b_mode_ref", to_string(self.b_mode_ref())}});
}
std::string to_string(const TaoUniverseCalcProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoUniverseCalcProxy",
      {std::pair{"srdt_for_data", to_string(self.srdt_for_data())},
       std::pair{"rad_int_for_data", to_string(self.rad_int_for_data())},
       std::pair{
           "rad_int_for_plotting", to_string(self.rad_int_for_plotting())},
       std::pair{"chrom_for_data", to_string(self.chrom_for_data())},
       std::pair{"chrom_for_plotting", to_string(self.chrom_for_plotting())},
       std::pair{"lat_sigma_for_data", to_string(self.lat_sigma_for_data())},
       std::pair{
           "lat_sigma_for_plotting", to_string(self.lat_sigma_for_plotting())},
       std::pair{"dynamic_aperture", to_string(self.dynamic_aperture())},
       std::pair{"one_turn_map", to_string(self.one_turn_map())},
       std::pair{"lattice", to_string(self.lattice())},
       std::pair{"twiss", to_string(self.twiss())},
       std::pair{"track", to_string(self.track())},
       std::pair{"spin_matrices", to_string(self.spin_matrices())}});
}
std::string to_string(const LatEleOrderProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "LatEleOrderProxy",
      {std::pair{"branch", "[...]"}});
}
std::string to_string(const TaoTitleProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoTitleProxy",
      {std::pair{"string", self.string()},
       std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"units", self.units()},
       std::pair{"justify", self.justify()},
       std::pair{"draw_it", to_string(self.draw_it())}});
}
std::string to_string(const QpRectProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "QpRectProxy",
      {std::pair{"x1", to_string(self.x1())},
       std::pair{"x2", to_string(self.x2())},
       std::pair{"y1", to_string(self.y1())},
       std::pair{"y2", to_string(self.y2())},
       std::pair{"units", self.units()}});
}
std::string to_string(const TaoDrawingProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoDrawingProxy",
      {std::pair{"ele_shape", "[...]"}});
}
std::string to_string(const TaoShapePatternProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoShapePatternProxy",
      {std::pair{"name", self.name()},
       std::pair{"line", to_string(self.line())},
       std::pair{"pt", "[...]"}});
}
std::string to_string(const TaoShapePatternPointProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoShapePatternPointProxy",
      {std::pair{"s", to_string(self.s())},
       std::pair{"y", to_string(self.y())},
       std::pair{"radius", to_string(self.radius())}});
}
std::string to_string(const QpAxisProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "QpAxisProxy",
      {std::pair{"label", self.label()},
       std::pair{"min", to_string(self.min())},
       std::pair{"max", to_string(self.max())},
       std::pair{"tick_min", to_string(self.tick_min())},
       std::pair{"tick_max", to_string(self.tick_max())},
       std::pair{"eval_min", to_string(self.eval_min())},
       std::pair{"eval_max", to_string(self.eval_max())},
       std::pair{"dtick", to_string(self.dtick())},
       std::pair{"number_offset", to_string(self.number_offset())},
       std::pair{"label_offset", to_string(self.label_offset())},
       std::pair{"major_tick_len", to_string(self.major_tick_len())},
       std::pair{"minor_tick_len", to_string(self.minor_tick_len())},
       std::pair{"label_color", self.label_color()},
       std::pair{"major_div", to_string(self.major_div())},
       std::pair{"major_div_nominal", to_string(self.major_div_nominal())},
       std::pair{"minor_div", to_string(self.minor_div())},
       std::pair{"minor_div_max", to_string(self.minor_div_max())},
       std::pair{"places", to_string(self.places())},
       std::pair{"type", self.type()},
       std::pair{"bounds", self.bounds()},
       std::pair{"tick_side", to_string(self.tick_side())},
       std::pair{"number_side", to_string(self.number_side())},
       std::pair{"draw_label", to_string(self.draw_label())},
       std::pair{"draw_numbers", to_string(self.draw_numbers())}});
}
std::string to_string(const QpLegendProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "QpLegendProxy",
      {std::pair{"row_spacing", to_string(self.row_spacing())},
       std::pair{"line_length", to_string(self.line_length())},
       std::pair{"text_offset", to_string(self.text_offset())},
       std::pair{"draw_line", to_string(self.draw_line())},
       std::pair{"draw_symbol", to_string(self.draw_symbol())},
       std::pair{"draw_text", to_string(self.draw_text())}});
}
std::string to_string(const QpPointProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "QpPointProxy",
      {std::pair{"x", to_string(self.x())},
       std::pair{"y", to_string(self.y())},
       std::pair{"units", self.units()}});
}
std::string to_string(const QpLineProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "QpLineProxy",
      {std::pair{"width", to_string(self.width())},
       std::pair{"color", self.color()},
       std::pair{"pattern", self.pattern()}});
}
std::string to_string(const QpSymbolProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "QpSymbolProxy",
      {std::pair{"type", self.type()},
       std::pair{"height", to_string(self.height())},
       std::pair{"color", self.color()},
       std::pair{"fill_pattern", self.fill_pattern()},
       std::pair{"line_width", to_string(self.line_width())}});
}
std::string to_string(const TaoFloorPlanProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoFloorPlanProxy",
      {std::pair{"view", self.view()},
       std::pair{"rotation", to_string(self.rotation())},
       std::pair{"correct_distortion", to_string(self.correct_distortion())},
       std::pair{"flip_label_side", to_string(self.flip_label_side())},
       std::pair{"size_is_absolute", to_string(self.size_is_absolute())},
       std::pair{
           "draw_only_first_pass", to_string(self.draw_only_first_pass())},
       std::pair{"draw_building_wall", to_string(self.draw_building_wall())},
       std::pair{"orbit_scale", to_string(self.orbit_scale())},
       std::pair{"orbit_color", self.orbit_color()},
       std::pair{"orbit_pattern", self.orbit_pattern()},
       std::pair{"orbit_lattice", self.orbit_lattice()},
       std::pair{"orbit_width", to_string(self.orbit_width())}});
}
std::string to_string(const TaoV1VarProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoV1VarProxy",
      {std::pair{"name", self.name()},
       std::pair{"ix_v1_var", to_string(self.ix_v1_var())},
       std::pair{"v", "[...]"}});
}
std::string to_string(const TaoGlobalProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoGlobalProxy",
      {std::pair{"beam_dead_cutoff", to_string(self.beam_dead_cutoff())},
       std::pair{"lm_opt_deriv_reinit", to_string(self.lm_opt_deriv_reinit())},
       std::pair{"de_lm_step_ratio", to_string(self.de_lm_step_ratio())},
       std::pair{
           "de_var_to_population_factor",
           to_string(self.de_var_to_population_factor())},
       std::pair{"lmdif_eps", to_string(self.lmdif_eps())},
       std::pair{
           "lmdif_negligible_merit", to_string(self.lmdif_negligible_merit())},
       std::pair{"svd_cutoff", to_string(self.svd_cutoff())},
       std::pair{"unstable_penalty", to_string(self.unstable_penalty())},
       std::pair{"merit_stop_value", to_string(self.merit_stop_value())},
       std::pair{"dmerit_stop_value", to_string(self.dmerit_stop_value())},
       std::pair{"random_sigma_cutoff", to_string(self.random_sigma_cutoff())},
       std::pair{"delta_e_chrom", to_string(self.delta_e_chrom())},
       std::pair{"max_plot_time", to_string(self.max_plot_time())},
       std::pair{"default_universe", to_string(self.default_universe())},
       std::pair{"default_branch", to_string(self.default_branch())},
       std::pair{"n_opti_cycles", to_string(self.n_opti_cycles())},
       std::pair{"n_opti_loops", to_string(self.n_opti_loops())},
       std::pair{"n_threads", to_string(self.n_threads())},
       std::pair{"phase_units", to_string(self.phase_units())},
       std::pair{"bunch_to_plot", to_string(self.bunch_to_plot())},
       std::pair{"random_seed", to_string(self.random_seed())},
       std::pair{"n_top10_merit", to_string(self.n_top10_merit())},
       std::pair{"srdt_gen_n_slices", to_string(self.srdt_gen_n_slices())},
       std::pair{
           "datum_err_messages_max", to_string(self.datum_err_messages_max())},
       std::pair{"srdt_sxt_n_slices", to_string(self.srdt_sxt_n_slices())},
       std::pair{"srdt_use_cache", to_string(self.srdt_use_cache())},
       std::pair{"quiet", self.quiet()},
       std::pair{"random_engine", self.random_engine()},
       std::pair{"random_gauss_converter", self.random_gauss_converter()},
       std::pair{"track_type", self.track_type()},
       std::pair{
           "lat_sigma_calc_uses_emit_from",
           self.lat_sigma_calc_uses_emit_from()},
       std::pair{"prompt_string", self.prompt_string()},
       std::pair{"prompt_color", self.prompt_color()},
       std::pair{"optimizer", self.optimizer()},
       std::pair{"print_command", self.print_command()},
       std::pair{"var_out_file", self.var_out_file()},
       std::pair{"history_file", self.history_file()},
       std::pair{"beam_timer_on", to_string(self.beam_timer_on())},
       std::pair{"box_plots", to_string(self.box_plots())},
       std::pair{
           "blank_line_between_commands",
           to_string(self.blank_line_between_commands())},
       std::pair{
           "cmd_file_abort_on_error",
           to_string(self.cmd_file_abort_on_error())},
       std::pair{"concatenate_maps", to_string(self.concatenate_maps())},
       std::pair{"derivative_recalc", to_string(self.derivative_recalc())},
       std::pair{
           "derivative_uses_design", to_string(self.derivative_uses_design())},
       std::pair{
           "disable_smooth_line_calc",
           to_string(self.disable_smooth_line_calc())},
       std::pair{
           "draw_curve_off_scale_warn",
           to_string(self.draw_curve_off_scale_warn())},
       std::pair{"external_plotting", to_string(self.external_plotting())},
       std::pair{
           "label_lattice_elements", to_string(self.label_lattice_elements())},
       std::pair{"label_keys", to_string(self.label_keys())},
       std::pair{"lattice_calc_on", to_string(self.lattice_calc_on())},
       std::pair{"only_limit_opt_vars", to_string(self.only_limit_opt_vars())},
       std::pair{"opt_with_ref", to_string(self.opt_with_ref())},
       std::pair{"opt_with_base", to_string(self.opt_with_base())},
       std::pair{
           "opt_match_auto_recalc", to_string(self.opt_match_auto_recalc())},
       std::pair{"opti_write_var_file", to_string(self.opti_write_var_file())},
       std::pair{
           "optimizer_allow_user_abort",
           to_string(self.optimizer_allow_user_abort())},
       std::pair{
           "optimizer_var_limit_warn",
           to_string(self.optimizer_var_limit_warn())},
       std::pair{"plot_on", to_string(self.plot_on())},
       std::pair{
           "rad_int_user_calc_on", to_string(self.rad_int_user_calc_on())},
       std::pair{"rf_on", to_string(self.rf_on())},
       std::pair{"single_step", to_string(self.single_step())},
       std::pair{"stop_on_error", to_string(self.stop_on_error())},
       std::pair{
           "svd_retreat_on_merit_increase",
           to_string(self.svd_retreat_on_merit_increase())},
       std::pair{"var_limits_on", to_string(self.var_limits_on())},
       std::pair{
           "wait_for_CR_in_single_mode",
           to_string(self.wait_for_CR_in_single_mode())},
       std::pair{"symbol_import", to_string(self.symbol_import())},
       std::pair{"debug_on", to_string(self.debug_on())},
       std::pair{"expression_tree_on", to_string(self.expression_tree_on())},
       std::pair{"verbose_on", to_string(self.verbose_on())}});
}
std::string to_string(const TaoInitProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoInitProxy",
      {std::pair{"parse_cmd_args", to_string(self.parse_cmd_args())},
       std::pair{"debug_switch", to_string(self.debug_switch())},
       std::pair{
           "external_plotting_switch",
           to_string(self.external_plotting_switch())},
       std::pair{"init_name", self.init_name()},
       std::pair{"hook_init_file", self.hook_init_file()},
       std::pair{"hook_lat_file", self.hook_lat_file()},
       std::pair{"hook_beam_file", self.hook_beam_file()},
       std::pair{"hook_data_file", self.hook_data_file()},
       std::pair{"hook_plot_file", self.hook_plot_file()},
       std::pair{"hook_startup_file", self.hook_startup_file()},
       std::pair{"hook_var_file", self.hook_var_file()},
       std::pair{"hook_building_wall_file", self.hook_building_wall_file()},
       std::pair{"init_file_arg_path", self.init_file_arg_path()},
       std::pair{"lattice_file_arg", self.lattice_file_arg()},
       std::pair{"hook_init_file_arg", self.hook_init_file_arg()},
       std::pair{"init_file_arg", self.init_file_arg()},
       std::pair{"beam_file_arg", self.beam_file_arg()},
       std::pair{
           "beam_init_position_file_arg", self.beam_init_position_file_arg()},
       std::pair{"command_arg", self.command_arg()},
       std::pair{"data_file_arg", self.data_file_arg()},
       std::pair{"plot_file_arg", self.plot_file_arg()},
       std::pair{"startup_file_arg", self.startup_file_arg()},
       std::pair{"var_file_arg", self.var_file_arg()},
       std::pair{"building_wall_file_arg", self.building_wall_file_arg()},
       std::pair{"geometry_arg", self.geometry_arg()},
       std::pair{"slice_lattice_arg", self.slice_lattice_arg()},
       std::pair{"start_branch_at_arg", self.start_branch_at_arg()},
       std::pair{"log_startup_arg", self.log_startup_arg()},
       std::pair{"no_stopping_arg", self.no_stopping_arg()},
       std::pair{"noplot_arg", self.noplot_arg()},
       std::pair{"no_rad_int_arg", self.no_rad_int_arg()},
       std::pair{"reverse_arg", self.reverse_arg()},
       std::pair{"debug_arg", self.debug_arg()},
       std::pair{
           "disable_smooth_line_calc_arg", self.disable_smooth_line_calc_arg()},
       std::pair{"rf_on_arg", self.rf_on_arg()},
       std::pair{"prompt_color_arg", self.prompt_color_arg()},
       std::pair{"quiet_arg", self.quiet_arg()},
       std::pair{"noinit_arg", self.noinit_arg()},
       std::pair{"nostartup_arg", self.nostartup_arg()},
       std::pair{"symbol_import_arg", self.symbol_import_arg()},
       std::pair{"unique_name_suffix", self.unique_name_suffix()}});
}
std::string to_string(const TaoCommonProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoCommonProxy",
      {std::pair{"plot_place_buffer", "[...]"},
       std::pair{"covar", to_string(self.covar())},
       std::pair{"alpha", to_string(self.alpha())},
       std::pair{"dummy_target", to_string(self.dummy_target())},
       std::pair{"n_alias", to_string(self.n_alias())},
       std::pair{"cmd_file_level", to_string(self.cmd_file_level())},
       std::pair{"ix_key_bank", to_string(self.ix_key_bank())},
       std::pair{"ix_history", to_string(self.ix_history())},
       std::pair{"n_history", to_string(self.n_history())},
       std::pair{"lev_loop", to_string(self.lev_loop())},
       std::pair{
           "n_err_messages_printed", to_string(self.n_err_messages_printed())},
       std::pair{"n_universes", to_string(self.n_universes())},
       std::pair{
           "ix_beam_track_active_element",
           to_string(self.ix_beam_track_active_element())},
       std::pair{"cmd_file_paused", to_string(self.cmd_file_paused())},
       std::pair{"use_cmd_here", to_string(self.use_cmd_here())},
       std::pair{"cmd_from_cmd_file", to_string(self.cmd_from_cmd_file())},
       std::pair{
           "use_saved_beam_in_tracking",
           to_string(self.use_saved_beam_in_tracking())},
       std::pair{"single_mode", to_string(self.single_mode())},
       std::pair{
           "combine_consecutive_elements_of_like_name",
           to_string(self.combine_consecutive_elements_of_like_name())},
       std::pair{"have_tracked_beam", to_string(self.have_tracked_beam())},
       std::pair{"init_plot_needed", to_string(self.init_plot_needed())},
       std::pair{"init_beam", to_string(self.init_beam())},
       std::pair{"init_var", to_string(self.init_var())},
       std::pair{"init_read_lat_info", to_string(self.init_read_lat_info())},
       std::pair{"optimizer_running", to_string(self.optimizer_running())},
       std::pair{
           "have_datums_using_expressions",
           to_string(self.have_datums_using_expressions())},
       std::pair{"print_to_terminal", to_string(self.print_to_terminal())},
       std::pair{"lattice_calc_done", to_string(self.lattice_calc_done())},
       std::pair{
           "add_measurement_noise", to_string(self.add_measurement_noise())},
       std::pair{
           "command_arg_has_been_executed",
           to_string(self.command_arg_has_been_executed())},
       std::pair{
           "all_merit_weights_positive",
           to_string(self.all_merit_weights_positive())},
       std::pair{
           "multi_turn_orbit_is_plotted",
           to_string(self.multi_turn_orbit_is_plotted())},
       std::pair{"force_chrom_calc", to_string(self.force_chrom_calc())},
       std::pair{"force_rad_int_calc", to_string(self.force_rad_int_calc())},
       std::pair{"rad_int_ri_calc_on", to_string(self.rad_int_ri_calc_on())},
       std::pair{"rad_int_6d_calc_on", to_string(self.rad_int_6d_calc_on())},
       std::pair{"valid_plot_who", to_string(self.valid_plot_who())},
       std::pair{"single_mode_buffer", self.single_mode_buffer()},
       std::pair{"cmd", self.cmd()},
       std::pair{"saved_cmd_line", self.saved_cmd_line()}});
}
std::string to_string(const TaoPlotPageProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoPlotPageProxy",
      {std::pair{"title", to_string(self.title())},
       std::pair{"subtitle", to_string(self.subtitle())},
       std::pair{"border", to_string(self.border())},
       std::pair{"floor_plan", to_string(self.floor_plan())},
       std::pair{"lat_layout", to_string(self.lat_layout())},
       std::pair{"pattern", "[...]"},
       std::pair{"template_", "[...]"},
       std::pair{"region", "[...]"},
       std::pair{"plot_display_type", self.plot_display_type()},
       std::pair{"size", to_string(self.size())},
       std::pair{"text_height", to_string(self.text_height())},
       std::pair{
           "main_title_text_scale", to_string(self.main_title_text_scale())},
       std::pair{
           "graph_title_text_scale", to_string(self.graph_title_text_scale())},
       std::pair{
           "axis_number_text_scale", to_string(self.axis_number_text_scale())},
       std::pair{
           "axis_label_text_scale", to_string(self.axis_label_text_scale())},
       std::pair{"legend_text_scale", to_string(self.legend_text_scale())},
       std::pair{
           "key_table_text_scale", to_string(self.key_table_text_scale())},
       std::pair{
           "floor_plan_shape_scale", to_string(self.floor_plan_shape_scale())},
       std::pair{
           "floor_plan_text_scale", to_string(self.floor_plan_text_scale())},
       std::pair{
           "lat_layout_shape_scale", to_string(self.lat_layout_shape_scale())},
       std::pair{
           "lat_layout_text_scale", to_string(self.lat_layout_text_scale())},
       std::pair{"n_curve_pts", to_string(self.n_curve_pts())},
       std::pair{"id_window", to_string(self.id_window())},
       std::pair{
           "delete_overlapping_plots",
           to_string(self.delete_overlapping_plots())},
       std::pair{
           "draw_graph_title_suffix",
           to_string(self.draw_graph_title_suffix())}});
}
std::string to_string(const TaoBuildingWallProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoBuildingWallProxy",
      {std::pair{"orientation", to_string(self.orientation())},
       std::pair{"section", "[...]"}});
}
std::string to_string(const TaoBuildingWallOrientationProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoBuildingWallOrientationProxy",
      {std::pair{"theta", to_string(self.theta())},
       std::pair{"x_offset", to_string(self.x_offset())},
       std::pair{"z_offset", to_string(self.z_offset())}});
}
std::string to_string(const TaoBuildingWallSectionProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoBuildingWallSectionProxy",
      {std::pair{"name", self.name()},
       std::pair{"constraint", self.constraint()},
       std::pair{"point", "[...]"}});
}
std::string to_string(const TaoBuildingWallPointProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoBuildingWallPointProxy",
      {std::pair{"z", to_string(self.z())},
       std::pair{"x", to_string(self.x())},
       std::pair{"radius", to_string(self.radius())},
       std::pair{"z_center", to_string(self.z_center())},
       std::pair{"x_center", to_string(self.x_center())}});
}
std::string to_string(const TaoWaveProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoWaveProxy",
      {std::pair{"data_type", self.data_type()},
       std::pair{"rms_rel_a", to_string(self.rms_rel_a())},
       std::pair{"rms_rel_b", to_string(self.rms_rel_b())},
       std::pair{"rms_rel_as", to_string(self.rms_rel_as())},
       std::pair{"rms_rel_bs", to_string(self.rms_rel_bs())},
       std::pair{"rms_rel_ar", to_string(self.rms_rel_ar())},
       std::pair{"rms_rel_br", to_string(self.rms_rel_br())},
       std::pair{"rms_rel_k", to_string(self.rms_rel_k())},
       std::pair{"rms_rel_ks", to_string(self.rms_rel_ks())},
       std::pair{"rms_rel_kr", to_string(self.rms_rel_kr())},
       std::pair{"rms_phi", to_string(self.rms_phi())},
       std::pair{"rms_phi_s", to_string(self.rms_phi_s())},
       std::pair{"rms_phi_r", to_string(self.rms_phi_r())},
       std::pair{"amp_ba_s", to_string(self.amp_ba_s())},
       std::pair{"amp_ba_r", to_string(self.amp_ba_r())},
       std::pair{"chi_a", to_string(self.chi_a())},
       std::pair{"chi_c", to_string(self.chi_c())},
       std::pair{"chi_ba", to_string(self.chi_ba())},
       std::pair{"amp_a", to_string(self.amp_a())},
       std::pair{"amp_b", to_string(self.amp_b())},
       std::pair{"amp_ba", to_string(self.amp_ba())},
       std::pair{"coef_a", to_string(self.coef_a())},
       std::pair{"coef_b", to_string(self.coef_b())},
       std::pair{"coef_ba", to_string(self.coef_ba())},
       std::pair{"n_func", to_string(self.n_func())},
       std::pair{"ix_a1", to_string(self.ix_a1())},
       std::pair{"ix_a2", to_string(self.ix_a2())},
       std::pair{"ix_b1", to_string(self.ix_b1())},
       std::pair{"ix_b2", to_string(self.ix_b2())},
       std::pair{"i_a1", to_string(self.i_a1())},
       std::pair{"i_a2", to_string(self.i_a2())},
       std::pair{"i_b1", to_string(self.i_b1())},
       std::pair{"i_b2", to_string(self.i_b2())},
       std::pair{"n_a", to_string(self.n_a())},
       std::pair{"n_b", to_string(self.n_b())},
       std::pair{"i_curve_wrap_pt", to_string(self.i_curve_wrap_pt())},
       std::pair{"ix_data", to_string(self.ix_data())},
       std::pair{"n_kick", to_string(self.n_kick())},
       std::pair{"kick", "[...]"},
       std::pair{"base_graph", to_string(self.base_graph())},
       std::pair{"region", to_string(self.region())},
       std::pair{"d1_dat", to_string(self.d1_dat())}});
}
std::string to_string(const TaoWaveKickPtProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoWaveKickPtProxy",
      {std::pair{"phi_s", to_string(self.phi_s())},
       std::pair{"phi_r", to_string(self.phi_r())},
       std::pair{"phi", to_string(self.phi())},
       std::pair{"amp", to_string(self.amp())},
       std::pair{"s", to_string(self.s())},
       std::pair{"ix_dat_before_kick", to_string(self.ix_dat_before_kick())},
       std::pair{"ele", to_string(self.ele())}});
}
std::string to_string(const TaoCmdHistoryProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoCmdHistoryProxy",
      {std::pair{"cmd", self.cmd()}, std::pair{"ix", to_string(self.ix())}});
}
std::string to_string(const TaoUniverseProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TaoUniverseProxy",
      {std::pair{"model", to_string(self.model())},
       std::pair{"design", to_string(self.design())},
       std::pair{"base", to_string(self.base())},
       std::pair{"beam", to_string(self.beam())},
       std::pair{"dynamic_aperture", to_string(self.dynamic_aperture())},
       std::pair{"model_branch", "[...]"},
       std::pair{"d2_data", "[...]"},
       std::pair{"data", "[...]"},
       std::pair{"ping_scale", to_string(self.ping_scale())},
       std::pair{"scratch_lat", to_string(self.scratch_lat())},
       std::pair{"calc", to_string(self.calc())},
       std::pair{"ele_order", to_string(self.ele_order())},
       std::pair{"spin_map", to_string(self.spin_map())},
       std::pair{"dModel_dVar", to_string(self.dModel_dVar())},
       std::pair{"ix_uni", to_string(self.ix_uni())},
       std::pair{"n_d2_data_used", to_string(self.n_d2_data_used())},
       std::pair{"n_data_used", to_string(self.n_data_used())},
       std::pair{"is_on", to_string(self.is_on())},
       std::pair{
           "design_same_as_previous",
           to_string(self.design_same_as_previous())},
       std::pair{"picked_uni", to_string(self.picked_uni())}});
}
std::string to_string(const MadEnergyProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "MadEnergyProxy",
      {std::pair{"total", to_string(self.total())},
       std::pair{"beta", to_string(self.beta())},
       std::pair{"gamma", to_string(self.gamma())},
       std::pair{"kinetic", to_string(self.kinetic())},
       std::pair{"p0c", to_string(self.p0c())},
       std::pair{"particle", to_string(self.particle())}});
}
std::string to_string(const MadMapProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "MadMapProxy",
      {std::pair{"k", to_string(self.k())},
       std::pair{"r", to_string(self.r())},
       std::pair{"t", to_string(self.t())}});
}
std::string to_string(const BbuStageProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BbuStageProxy",
      {std::pair{"ix_ele_lr_wake", to_string(self.ix_ele_lr_wake())},
       std::pair{"ix_ele_stage_end", to_string(self.ix_ele_stage_end())},
       std::pair{"ix_pass", to_string(self.ix_pass())},
       std::pair{"ix_stage_pass1", to_string(self.ix_stage_pass1())},
       std::pair{"ix_head_bunch", to_string(self.ix_head_bunch())},
       std::pair{"ix_hom_max", to_string(self.ix_hom_max())},
       std::pair{"hom_voltage_max", to_string(self.hom_voltage_max())},
       std::pair{"time_at_wake_ele", to_string(self.time_at_wake_ele())},
       std::pair{"ave_orb", to_string(self.ave_orb())},
       std::pair{"rms_orb", to_string(self.rms_orb())},
       std::pair{"min_orb", to_string(self.min_orb())},
       std::pair{"max_orb", to_string(self.max_orb())},
       std::pair{"n_orb", to_string(self.n_orb())}});
}
std::string to_string(const BbuBeamProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BbuBeamProxy",
      {std::pair{"bunch", "[...]"},
       std::pair{"stage", "[...]"},
       std::pair{"ix_ele_bunch", to_string(self.ix_ele_bunch())},
       std::pair{"ix_bunch_head", to_string(self.ix_bunch_head())},
       std::pair{"ix_bunch_end", to_string(self.ix_bunch_end())},
       std::pair{"n_bunch_in_lat", to_string(self.n_bunch_in_lat())},
       std::pair{
           "ix_stage_voltage_max", to_string(self.ix_stage_voltage_max())},
       std::pair{"hom_voltage_max", to_string(self.hom_voltage_max())},
       std::pair{"time_now", to_string(self.time_now())},
       std::pair{"one_turn_time", to_string(self.one_turn_time())},
       std::pair{"rf_wavelength_max", to_string(self.rf_wavelength_max())}});
}
std::string to_string(const BbuParamProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "BbuParamProxy",
      {std::pair{"lat_filename", self.lat_filename()},
       std::pair{"lat2_filename", self.lat2_filename()},
       std::pair{"bunch_by_bunch_info_file", self.bunch_by_bunch_info_file()},
       std::pair{"hybridize", to_string(self.hybridize())},
       std::pair{
           "write_digested_hybrid_lat",
           to_string(self.write_digested_hybrid_lat())},
       std::pair{
           "write_voltage_vs_time_dat",
           to_string(self.write_voltage_vs_time_dat())},
       std::pair{
           "keep_overlays_and_groups",
           to_string(self.keep_overlays_and_groups())},
       std::pair{"keep_all_lcavities", to_string(self.keep_all_lcavities())},
       std::pair{
           "use_taylor_for_hybrids", to_string(self.use_taylor_for_hybrids())},
       std::pair{"stable_orbit_anal", to_string(self.stable_orbit_anal())},
       std::pair{"limit_factor", to_string(self.limit_factor())},
       std::pair{
           "simulation_turns_max", to_string(self.simulation_turns_max())},
       std::pair{"bunch_freq", to_string(self.bunch_freq())},
       std::pair{
           "init_particle_offset", to_string(self.init_particle_offset())},
       std::pair{"current", to_string(self.current())},
       std::pair{"rel_tol", to_string(self.rel_tol())},
       std::pair{"drscan", to_string(self.drscan())},
       std::pair{
           "use_interpolated_threshold",
           to_string(self.use_interpolated_threshold())},
       std::pair{"write_hom_info", to_string(self.write_hom_info())},
       std::pair{"elindex", to_string(self.elindex())},
       std::pair{"elname", self.elname()},
       std::pair{"nstep", to_string(self.nstep())},
       std::pair{"begdr", to_string(self.begdr())},
       std::pair{"enddr", to_string(self.enddr())},
       std::pair{"nrep", to_string(self.nrep())},
       std::pair{"ran_seed", to_string(self.ran_seed())},
       std::pair{"hom_order_cutoff", to_string(self.hom_order_cutoff())},
       std::pair{"ran_gauss_sigma_cut", to_string(self.ran_gauss_sigma_cut())},
       std::pair{"ele_track_end", self.ele_track_end()},
       std::pair{"ix_ele_track_end", to_string(self.ix_ele_track_end())},
       std::pair{"regression", to_string(self.regression())},
       std::pair{"normalize_z_to_rf", to_string(self.normalize_z_to_rf())},
       std::pair{"ramp_on", to_string(self.ramp_on())},
       std::pair{"ramp_pattern", to_string(self.ramp_pattern())},
       std::pair{"ramp_n_start", to_string(self.ramp_n_start())},
       std::pair{"n_ramp_pattern", to_string(self.n_ramp_pattern())}});
}
std::string to_string(const AllEncompassingProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "AllEncompassingProxy",
      {std::pair{"real_rp_0d", to_string(self.real_rp_0d())},
       std::pair{"real_rp_1d", to_string(self.real_rp_1d())},
       std::pair{"real_rp_2d", to_string(self.real_rp_2d())},
       std::pair{"real_rp_3d", to_string(self.real_rp_3d())},
       std::pair{"real_rp_0d_ptr", to_string(self.real_rp_0d_ptr())},
       std::pair{"real_rp_1d_ptr", to_string(self.real_rp_1d_ptr())},
       std::pair{"real_rp_2d_ptr", to_string(self.real_rp_2d_ptr())},
       std::pair{"real_rp_3d_ptr", to_string(self.real_rp_3d_ptr())},
       std::pair{"real_rp_1d_alloc", to_string(self.real_rp_1d_alloc())},
       std::pair{"real_rp_2d_alloc", to_string(self.real_rp_2d_alloc())},
       std::pair{"real_rp_3d_alloc", to_string(self.real_rp_3d_alloc())},
       std::pair{"real_dp_0d", to_string(self.real_dp_0d())},
       std::pair{"real_dp_1d", to_string(self.real_dp_1d())},
       std::pair{"real_dp_2d", to_string(self.real_dp_2d())},
       std::pair{"real_dp_3d", to_string(self.real_dp_3d())},
       std::pair{"real_dp_0d_ptr", to_string(self.real_dp_0d_ptr())},
       std::pair{"real_dp_1d_ptr", to_string(self.real_dp_1d_ptr())},
       std::pair{"real_dp_2d_ptr", to_string(self.real_dp_2d_ptr())},
       std::pair{"real_dp_3d_ptr", to_string(self.real_dp_3d_ptr())},
       std::pair{"real_dp_1d_alloc", to_string(self.real_dp_1d_alloc())},
       std::pair{"real_dp_2d_alloc", to_string(self.real_dp_2d_alloc())},
       std::pair{"real_dp_3d_alloc", to_string(self.real_dp_3d_alloc())},
       std::pair{"complex_dp_0d", to_string(self.complex_dp_0d())},
       std::pair{"complex_dp_1d", to_string(self.complex_dp_1d())},
       std::pair{"complex_dp_2d", to_string(self.complex_dp_2d())},
       std::pair{"complex_dp_3d", to_string(self.complex_dp_3d())},
       std::pair{"complex_dp_1d_ptr", to_string(self.complex_dp_1d_ptr())},
       std::pair{"complex_dp_2d_ptr", to_string(self.complex_dp_2d_ptr())},
       std::pair{"complex_dp_3d_ptr", to_string(self.complex_dp_3d_ptr())},
       std::pair{"complex_dp_1d_alloc", to_string(self.complex_dp_1d_alloc())},
       std::pair{"complex_dp_2d_alloc", to_string(self.complex_dp_2d_alloc())},
       std::pair{"complex_dp_3d_alloc", to_string(self.complex_dp_3d_alloc())},
       std::pair{"int_0d", to_string(self.int_0d())},
       std::pair{"int_1d", to_string(self.int_1d())},
       std::pair{"int_2d", to_string(self.int_2d())},
       std::pair{"int_3d", to_string(self.int_3d())},
       std::pair{"int_0d_ptr", to_string(self.int_0d_ptr())},
       std::pair{"int_1d_ptr", to_string(self.int_1d_ptr())},
       std::pair{"int_2d_ptr", to_string(self.int_2d_ptr())},
       std::pair{"int_3d_ptr", to_string(self.int_3d_ptr())},
       std::pair{"int_1d_alloc", to_string(self.int_1d_alloc())},
       std::pair{"int_2d_alloc", to_string(self.int_2d_alloc())},
       std::pair{"int_3d_alloc", to_string(self.int_3d_alloc())},
       std::pair{"int8_0d", to_string(self.int8_0d())},
       std::pair{"int8_0d_ptr", to_string(self.int8_0d_ptr())},
       std::pair{"logical_0d", to_string(self.logical_0d())},
       std::pair{"logical_0d_ptr", to_string(self.logical_0d_ptr())},
       std::pair{"type_0d", to_string(self.type_0d())},
       std::pair{"type_1d", "[...]"},
       std::pair{"type_2d", "[...]"},
       std::pair{"type_3d", "[...]"},
       std::pair{"type_0d_ptr", to_string(self.type_0d_ptr())},
       std::pair{"type_1d_ptr", "[...]"},
       std::pair{"type_2d_ptr", "[...]"},
       std::pair{"type_3d_ptr", "[...]"},
       std::pair{"type_1d_alloc", "[...]"},
       std::pair{"type_2d_alloc", "[...]"},
       std::pair{"type_3d_alloc", "[...]"}});
}
std::string to_string(const TestSubProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TestSubProxy",
      {std::pair{"sr", to_string(self.sr())}});
}
std::string to_string(const TestSubSubProxy& self) {
  return repr(
      self.get_fortran_ptr(),
      "TestSubSubProxy",
      {std::pair{"aaa", to_string(self.aaa())},
       std::pair{"bbb", to_string(self.bbb())},
       std::pair{"file", self.file()},
       std::pair{"t_ref", to_string(self.t_ref())},
       std::pair{"freq_spread", to_string(self.freq_spread())}});
}
std::string to_string(const Bmad::AbMultipoleKick& self) {
  return repr(
      &self,
      "Bmad::AbMultipoleKick",
      {std::pair{"kx", to_string(self.kx)},
       std::pair{"ky", to_string(self.ky)},
       std::pair{"dk", to_string(self.dk)}});
}
std::string to_string(const Bmad::ActionToXyz& self) {
  return repr(
      &self,
      "Bmad::ActionToXyz",
      {std::pair{"X", to_string(self.X)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::AddSuperimpose& self) {
  return repr(
      &self,
      "Bmad::AddSuperimpose",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"super_ele_out", to_string(self.super_ele_out)}});
}
std::string to_string(const SimUtils::ApfftCorr& self) {
  return repr(
      &self,
      "SimUtils::ApfftCorr",
      {std::pair{"phase", to_string(self.phase)},
       std::pair{"amp", to_string(self.amp)},
       std::pair{"freq", to_string(self.freq)}});
}
std::string to_string(const Bmad::BbiKick& self) {
  return repr(
      &self,
      "Bmad::BbiKick",
      {std::pair{"nk", to_string(self.nk)},
       std::pair{"dnk", to_string(self.dnk)}});
}
std::string to_string(const Bmad::BeamTilts& self) {
  return repr(
      &self,
      "Bmad::BeamTilts",
      {std::pair{"angle_xy", to_string(self.angle_xy)},
       std::pair{"angle_xz", to_string(self.angle_xz)},
       std::pair{"angle_yz", to_string(self.angle_yz)},
       std::pair{"angle_xpz", to_string(self.angle_xpz)},
       std::pair{"angle_ypz", to_string(self.angle_ypz)}});
}
std::string to_string(const SimUtils::BicubicCmplxEval& self) {
  return repr(
      &self,
      "SimUtils::BicubicCmplxEval",
      {std::pair{"df_dx", to_string(self.df_dx)},
       std::pair{"df_dy", to_string(self.df_dy)},
       std::pair{"f_val", to_string(self.f_val)}});
}
std::string to_string(const Bmad::BmadParser& self) {
  return repr(
      &self,
      "Bmad::BmadParser",
      {std::pair{"lat", to_string(self.lat)},
       std::pair{"digested_read_ok", to_string(self.digested_read_ok)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"parse_lat", to_string(self.parse_lat)}});
}
std::string to_string(const SimUtils::BracketIndexForSpline& self) {
  return repr(
      &self,
      "SimUtils::BracketIndexForSpline",
      {std::pair{"ix0", to_string(self.ix0)},
       std::pair{"ok", to_string(self.ok)}});
}
std::string to_string(const Bmad::CalcEmittancesAndTwissFromSigmaMatrix& self) {
  return repr(
      &self,
      "Bmad::CalcEmittancesAndTwissFromSigmaMatrix",
      {std::pair{"bunch_params", to_string(self.bunch_params)},
       std::pair{"error", to_string(self.error)},
       std::pair{"n_mat", to_string(self.n_mat)}});
}
std::string to_string(const Bmad::CalcWallRadius& self) {
  return repr(
      &self,
      "Bmad::CalcWallRadius",
      {std::pair{"r_wall", to_string(self.r_wall)},
       std::pair{"dr_dtheta", to_string(self.dr_dtheta)},
       std::pair{"ix_vertex", to_string(self.ix_vertex)}});
}
std::string to_string(const Bmad::CheckIfSInBounds& self) {
  return repr(
      &self,
      "Bmad::CheckIfSInBounds",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"translated_s", to_string(self.translated_s)}});
}
std::string to_string(const Bmad::ChooseQuadsForSetTune& self) {
  return repr(
      &self,
      "Bmad::ChooseQuadsForSetTune",
      {std::pair{"dk1", to_string(self.dk1)},
       std::pair{"eles", "[...]"},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ChromCalc& self) {
  return repr(
      &self,
      "Bmad::ChromCalc",
      {std::pair{"chrom_a", to_string(self.chrom_a)},
       std::pair{"chrom_b", to_string(self.chrom_b)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"low_E_lat", to_string(self.low_E_lat)},
       std::pair{"high_E_lat", to_string(self.high_E_lat)},
       std::pair{"low_E_orb", "[...]"},
       std::pair{"high_E_orb", "[...]"}});
}
std::string to_string(const Bmad::ClosedOrbitFromTracking& self) {
  return repr(
      &self,
      "Bmad::ClosedOrbitFromTracking",
      {std::pair{"closed_orb", "[...]"},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ComplexTaylorToMat6& self) {
  return repr(
      &self,
      "Bmad::ComplexTaylorToMat6",
      {std::pair{"vec0", to_string(self.vec0)},
       std::pair{"mat6", to_string(self.mat6)},
       std::pair{"r_out", to_string(self.r_out)}});
}
std::string to_string(const Bmad::ConvertCoords& self) {
  return repr(
      &self,
      "Bmad::ConvertCoords",
      {std::pair{"out_type_str", self.out_type_str},
       std::pair{"coord_out", to_string(self.coord_out)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ConvertPcTo& self) {
  return repr(
      &self,
      "Bmad::ConvertPcTo",
      {std::pair{"E_tot", to_string(self.E_tot)},
       std::pair{"gamma", to_string(self.gamma)},
       std::pair{"kinetic", to_string(self.kinetic)},
       std::pair{"beta", to_string(self.beta)},
       std::pair{"brho", to_string(self.brho)},
       std::pair{"beta1", to_string(self.beta1)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ConvertTotalEnergyTo& self) {
  return repr(
      &self,
      "Bmad::ConvertTotalEnergyTo",
      {std::pair{"gamma", to_string(self.gamma)},
       std::pair{"kinetic", to_string(self.kinetic)},
       std::pair{"beta", to_string(self.beta)},
       std::pair{"pc", to_string(self.pc)},
       std::pair{"brho", to_string(self.brho)},
       std::pair{"beta1", to_string(self.beta1)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ConverterDistributionParser& self) {
  return repr(
      &self,
      "Bmad::ConverterDistributionParser",
      {std::pair{"delim", self.delim},
       std::pair{"delim_found", to_string(self.delim_found)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::CoordsFloorToCurvilinear& self) {
  return repr(
      &self,
      "Bmad::CoordsFloorToCurvilinear",
      {std::pair{"ele1", to_string(self.ele1)},
       std::pair{"status", to_string(self.status)},
       std::pair{"w_mat", to_string(self.w_mat)}});
}
std::string to_string(const Bmad::CoordsFloorToLocalCurvilinear& self) {
  return repr(
      &self,
      "Bmad::CoordsFloorToLocalCurvilinear",
      {std::pair{"status", to_string(self.status)},
       std::pair{"w_mat", to_string(self.w_mat)}});
}
std::string to_string(const Bmad::CreateElementSlice& self) {
  return repr(
      &self,
      "Bmad::CreateElementSlice",
      {std::pair{"sliced_ele", to_string(self.sliced_ele)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::CreatePlanarWigglerModel& self) {
  return repr(
      &self,
      "Bmad::CreatePlanarWigglerModel",
      {std::pair{"lat", to_string(self.lat)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::EigenDecomp6mat& self) {
  return repr(
      &self,
      "Bmad::EigenDecomp6mat",
      {std::pair{"eval", to_string(self.eval)},
       std::pair{"evec", to_string(self.evec)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"tunes", to_string(self.tunes)}});
}
std::string to_string(const Bmad::EleMisalignmentLSCalc& self) {
  return repr(
      &self,
      "Bmad::EleMisalignmentLSCalc",
      {std::pair{"L_mis", to_string(self.L_mis)},
       std::pair{"S_mis", to_string(self.S_mis)}});
}
std::string to_string(const Bmad::EleToPtcMagneticBnAn& self) {
  return repr(
      &self,
      "Bmad::EleToPtcMagneticBnAn",
      {std::pair{"bn", to_string(self.bn)},
       std::pair{"an", to_string(self.an)},
       std::pair{"n_max", to_string(self.n_max)}});
}
std::string to_string(const Bmad::EleToTaylor& self) {
  return repr(
      &self,
      "Bmad::EleToTaylor",
      {std::pair{"orbital_taylor", "[...]"},
       std::pair{"spin_taylor", "[...]"}});
}
std::string to_string(const Bmad::ElecMultipoleField& self) {
  return repr(
      &self,
      "Bmad::ElecMultipoleField",
      {std::pair{"Ex", to_string(self.Ex)},
       std::pair{"Ey", to_string(self.Ey)},
       std::pair{"dE", to_string(self.dE)},
       std::pair{"compute_dE", to_string(self.compute_dE)}});
}
std::string to_string(const Bmad::EmFieldCalc& self) {
  return repr(
      &self,
      "Bmad::EmFieldCalc",
      {std::pair{"field", to_string(self.field)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::Emit6d& self) {
  return repr(
      &self,
      "Bmad::Emit6d",
      {std::pair{"mode", to_string(self.mode)},
       std::pair{"sigma_mat", to_string(self.sigma_mat)},
       std::pair{"rad_int_by_ele", to_string(self.rad_int_by_ele)}});
}
std::string to_string(const Bmad::EnvelopeRadintsIbs& self) {
  return repr(
      &self,
      "Bmad::EnvelopeRadintsIbs",
      {std::pair{"alpha", to_string(self.alpha)},
       std::pair{"emit", to_string(self.emit)}});
}
std::string to_string(const Bmad::EvaluateArrayIndex& self) {
  return repr(
      &self,
      "Bmad::EvaluateArrayIndex",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"word2", self.word2},
       std::pair{"delim2", self.delim2},
       std::pair{"this_index", to_string(self.this_index)}});
}
std::string to_string(const Bmad::EvaluateLogical& self) {
  return repr(
      &self,
      "Bmad::EvaluateLogical",
      {std::pair{"iostat", to_string(self.iostat)},
       std::pair{"this_logic", to_string(self.this_logic)}});
}
std::string to_string(const Bmad::ExpectThis& self) {
  return repr(
      &self,
      "Bmad::ExpectThis",
      {std::pair{"delim", self.delim},
       std::pair{"delim_found", to_string(self.delim_found)},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::ExpressionStackValue& self) {
  return repr(
      &self,
      "Bmad::ExpressionStackValue",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"err_str", self.err_str},
       std::pair{"value", to_string(self.value)}});
}
std::string to_string(const Bmad::ExpressionStringToStack& self) {
  return repr(
      &self,
      "Bmad::ExpressionStringToStack",
      {std::pair{"stack", "[...]"},
       std::pair{"n_stack", to_string(self.n_stack)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"err_str", self.err_str}});
}
std::string to_string(const Bmad::ExpressionStringToTree& self) {
  return repr(
      &self,
      "Bmad::ExpressionStringToTree",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"err_str", self.err_str}});
}
std::string to_string(const Bmad::ExpressionValue& self) {
  return repr(
      &self,
      "Bmad::ExpressionValue",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"err_str", self.err_str},
       std::pair{"value", to_string(self.value)}});
}
std::string to_string(const Bmad::FindElementEnds& self) {
  return repr(
      &self,
      "Bmad::FindElementEnds",
      {std::pair{"ele1", to_string(self.ele1)},
       std::pair{"ele2", to_string(self.ele2)}});
}
std::string to_string(const Bmad::FindMatchingFieldmap& self) {
  return repr(
      &self,
      "Bmad::FindMatchingFieldmap",
      {std::pair{"match_ele", to_string(self.match_ele)},
       std::pair{"ix_field", to_string(self.ix_field)}});
}
std::string to_string(const Bmad::FloorAnglesToWMat& self) {
  return repr(
      &self,
      "Bmad::FloorAnglesToWMat",
      {std::pair{"w_mat", to_string(self.w_mat)},
       std::pair{"w_mat_inv", to_string(self.w_mat_inv)}});
}
std::string to_string(const Bmad::FloorWMatToAngles& self) {
  return repr(
      &self,
      "Bmad::FloorWMatToAngles",
      {std::pair{"theta", to_string(self.theta)},
       std::pair{"phi", to_string(self.phi)},
       std::pair{"psi", to_string(self.psi)}});
}
std::string to_string(const Bmad::FormDigestedBmadFileName& self) {
  return repr(
      &self,
      "Bmad::FormDigestedBmadFileName",
      {std::pair{"digested_file", self.digested_file},
       std::pair{"full_lat_file", self.full_lat_file}});
}
std::string to_string(const SimUtils::FourierAmplitude& self) {
  return repr(
      &self,
      "SimUtils::FourierAmplitude",
      {std::pair{"cos_amp", to_string(self.cos_amp)},
       std::pair{"sin_amp", to_string(self.sin_amp)},
       std::pair{"dcos_amp", to_string(self.dcos_amp)},
       std::pair{"dsin_amp", to_string(self.dsin_amp)}});
}
std::string to_string(const Bmad::GBendingStrengthFromEmField& self) {
  return repr(
      &self,
      "Bmad::GBendingStrengthFromEmField",
      {std::pair{"g", to_string(self.g)}, std::pair{"dg", to_string(self.dg)}});
}
std::string to_string(const Bmad::GetEmitFromSigmaMat& self) {
  return repr(
      &self,
      "Bmad::GetEmitFromSigmaMat",
      {std::pair{"normal", to_string(self.normal)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::GetSlaveList& self) {
  return repr(
      &self,
      "Bmad::GetSlaveList",
      {std::pair{"slaves", "[...]"},
       std::pair{"n_slave", to_string(self.n_slave)}});
}
std::string to_string(const Bmad::GptToParticleBunch& self) {
  return repr(
      &self,
      "Bmad::GptToParticleBunch",
      {std::pair{"bunch", to_string(self.bunch)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::InitBeamDistribution& self) {
  return repr(
      &self,
      "Bmad::InitBeamDistribution",
      {std::pair{"beam", to_string(self.beam)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"beam_init_set", to_string(self.beam_init_set)}});
}
std::string to_string(const Bmad::InitBunchDistribution& self) {
  return repr(
      &self,
      "Bmad::InitBunchDistribution",
      {std::pair{"bunch", to_string(self.bunch)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"beam_init_used", to_string(self.beam_init_used)}});
}
std::string to_string(const Bmad::InitPhotonIntegProb& self) {
  return repr(
      &self,
      "Bmad::InitPhotonIntegProb",
      {std::pair{"E_photon", to_string(self.E_photon)},
       std::pair{"integ_prob", to_string(self.integ_prob)}});
}
std::string to_string(const Bmad::KickVectorCalc& self) {
  return repr(
      &self,
      "Bmad::KickVectorCalc",
      {std::pair{"dr_ds", to_string(self.dr_ds)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::LinearCoef& self) {
  return repr(
      &self,
      "Bmad::LinearCoef",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"coef", to_string(self.coef)}});
}
std::string to_string(const Bmad::LoadParseLine& self) {
  return repr(
      &self,
      "Bmad::LoadParseLine",
      {std::pair{"end_of_file", to_string(self.end_of_file)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::MadTmfoc& self) {
  return repr(
      &self,
      "Bmad::MadTmfoc",
      {std::pair{"c", to_string(self.c)},
       std::pair{"s", to_string(self.s)},
       std::pair{"d", to_string(self.d)},
       std::pair{"f", to_string(self.f)}});
}
std::string to_string(const Bmad::MakeGMats& self) {
  return repr(
      &self,
      "Bmad::MakeGMats",
      {std::pair{"g_mat", to_string(self.g_mat)},
       std::pair{"g_inv_mat", to_string(self.g_inv_mat)}});
}
std::string to_string(const Bmad::MakeHvbp& self) {
  return repr(
      &self,
      "Bmad::MakeHvbp",
      {std::pair{"B", to_string(self.B)},
       std::pair{"V", to_string(self.V)},
       std::pair{"H", to_string(self.H)},
       std::pair{"Vbar", to_string(self.Vbar)},
       std::pair{"Hbar", to_string(self.Hbar)}});
}
std::string to_string(const Bmad::MakeMadMap& self) {
  return repr(
      &self,
      "Bmad::MakeMadMap",
      {std::pair{"energy", to_string(self.energy)},
       std::pair{"map", to_string(self.map)}});
}
std::string to_string(const Bmad::MakeMat6& self) {
  return repr(
      &self,
      "Bmad::MakeMat6",
      {std::pair{"end_orb", to_string(self.end_orb)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::MakeMat6Bmad& self) {
  return repr(
      &self,
      "Bmad::MakeMat6Bmad",
      {std::pair{"end_orb", to_string(self.end_orb)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::MakeMat6BmadPhoton& self) {
  return repr(
      &self,
      "Bmad::MakeMat6BmadPhoton",
      {std::pair{"end_orb", to_string(self.end_orb)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::MakeMat6Tracking& self) {
  return repr(
      &self,
      "Bmad::MakeMat6Tracking",
      {std::pair{"end_orb", to_string(self.end_orb)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::MakeN& self) {
  return repr(
      &self,
      "Bmad::MakeN",
      {std::pair{"N", to_string(self.N)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"tunes_out", to_string(self.tunes_out)},
       std::pair{"U", to_string(self.U)}});
}
std::string to_string(const Bmad::MakePbrh& self) {
  return repr(
      &self,
      "Bmad::MakePbrh",
      {std::pair{"P", to_string(self.P)},
       std::pair{"Bp", to_string(self.Bp)},
       std::pair{"R", to_string(self.R)},
       std::pair{"H", to_string(self.H)}});
}
std::string to_string(const Bmad::MakeSmatFromAbc& self) {
  return repr(
      &self,
      "Bmad::MakeSmatFromAbc",
      {std::pair{"sigma_mat", to_string(self.sigma_mat)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"Nout", to_string(self.Nout)}});
}
std::string to_string(const Bmad::MakeVMats& self) {
  return repr(
      &self,
      "Bmad::MakeVMats",
      {std::pair{"v_mat", to_string(self.v_mat)},
       std::pair{"v_inv_mat", to_string(self.v_inv_mat)}});
}
std::string to_string(const Bmad::MatSympDecouple& self) {
  return repr(
      &self,
      "Bmad::MatSympDecouple",
      {std::pair{"stat", to_string(self.stat)},
       std::pair{"twiss1", to_string(self.twiss1)},
       std::pair{"twiss2", to_string(self.twiss2)},
       std::pair{"gamma", to_string(self.gamma)}});
}
std::string to_string(const Bmad::MatchEleToMat6& self) {
  return repr(
      &self,
      "Bmad::MatchEleToMat6",
      {std::pair{"mat6", to_string(self.mat6)},
       std::pair{"vec0", to_string(self.vec0)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::MultiTurnTrackingAnalysis& self) {
  return repr(
      &self,
      "Bmad::MultiTurnTrackingAnalysis",
      {std::pair{"track0", to_string(self.track0)},
       std::pair{"ele", to_string(self.ele)},
       std::pair{"stable", to_string(self.stable)},
       std::pair{"growth_rate", to_string(self.growth_rate)},
       std::pair{"chi", to_string(self.chi)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::Multipole1AbToKt& self) {
  return repr(
      &self,
      "Bmad::Multipole1AbToKt",
      {std::pair{"knl", to_string(self.knl)},
       std::pair{"tn", to_string(self.tn)}});
}
std::string to_string(const Bmad::Multipole1KtToAb& self) {
  return repr(
      &self,
      "Bmad::Multipole1KtToAb",
      {std::pair{"an", to_string(self.an)},
       std::pair{"bn", to_string(self.bn)}});
}
std::string to_string(const Bmad::MultipoleAbToKt& self) {
  return repr(
      &self,
      "Bmad::MultipoleAbToKt",
      {std::pair{"knl", to_string(self.knl)},
       std::pair{"tn", to_string(self.tn)}});
}
std::string to_string(const Bmad::MultipoleEleToAb& self) {
  return repr(
      &self,
      "Bmad::MultipoleEleToAb",
      {std::pair{"ix_pole_max", to_string(self.ix_pole_max)},
       std::pair{"a", to_string(self.a)},
       std::pair{"b", to_string(self.b)},
       std::pair{"b1", to_string(self.b1)}});
}
std::string to_string(const Bmad::MultipoleEleToKt& self) {
  return repr(
      &self,
      "Bmad::MultipoleEleToKt",
      {std::pair{"ix_pole_max", to_string(self.ix_pole_max)},
       std::pair{"knl", to_string(self.knl)},
       std::pair{"tilt", to_string(self.tilt)}});
}
std::string to_string(const Bmad::MultipoleKtToAb& self) {
  return repr(
      &self,
      "Bmad::MultipoleKtToAb",
      {std::pair{"an", to_string(self.an)},
       std::pair{"bn", to_string(self.bn)}});
}
std::string to_string(const Bmad::NormalFormTaylors& self) {
  return repr(
      &self,
      "Bmad::NormalFormTaylors",
      {std::pair{"dhdj", "[...]"},
       std::pair{"A", "[...]"},
       std::pair{"A_inverse", "[...]"}});
}
std::string to_string(const Bmad::NormalMode3Calc& self) {
  return repr(
      &self,
      "Bmad::NormalMode3Calc",
      {std::pair{"tune", to_string(self.tune)},
       std::pair{"B", to_string(self.B)},
       std::pair{"HV", to_string(self.HV)}});
}
std::string to_string(const Bmad::OdeintBmad& self) {
  return repr(
      &self,
      "Bmad::OdeintBmad",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"track", to_string(self.track)}});
}
std::string to_string(const Bmad::OdeintBmadTime& self) {
  return repr(
      &self,
      "Bmad::OdeintBmadTime",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"dt_step", to_string(self.dt_step)}});
}
std::string to_string(const Bmad::OffsetParticle& self) {
  return repr(
      &self,
      "Bmad::OffsetParticle",
      {std::pair{"s_out", to_string(self.s_out)},
       std::pair{"spin_qrot", to_string(self.spin_qrot)}});
}
std::string to_string(const Bmad::OpenBinaryFile& self) {
  return repr(
      &self,
      "Bmad::OpenBinaryFile",
      {std::pair{"iu", to_string(self.iu)},
       std::pair{"iver", to_string(self.iver)},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::OrbitAmplitudeCalc& self) {
  return repr(
      &self,
      "Bmad::OrbitAmplitudeCalc",
      {std::pair{"amp_a", to_string(self.amp_a)},
       std::pair{"amp_b", to_string(self.amp_b)},
       std::pair{"amp_na", to_string(self.amp_na)},
       std::pair{"amp_nb", to_string(self.amp_nb)}});
}
std::string to_string(const Bmad::OrderEvecsByNSimilarity& self) {
  return repr(
      &self,
      "Bmad::OrderEvecsByNSimilarity",
      {std::pair{"evec", to_string(self.evec)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ParseIntegerList2& self) {
  return repr(
      &self,
      "Bmad::ParseIntegerList2",
      {std::pair{"num_found", to_string(self.num_found)},
       std::pair{"delim", self.delim},
       std::pair{"delim_found", to_string(self.delim_found)},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::ParseRealList& self) {
  return repr(
      &self,
      "Bmad::ParseRealList",
      {std::pair{"real_array", to_string(self.real_array)},
       std::pair{"delim", self.delim},
       std::pair{"delim_found", to_string(self.delim_found)},
       std::pair{"num_found", to_string(self.num_found)},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::ParseRealList2& self) {
  return repr(
      &self,
      "Bmad::ParseRealList2",
      {std::pair{"num_found", to_string(self.num_found)},
       std::pair{"delim", self.delim},
       std::pair{"delim_found", to_string(self.delim_found)},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::ParserFastComplexRead& self) {
  return repr(
      &self,
      "Bmad::ParserFastComplexRead",
      {std::pair{"cmplx_vec", to_string(self.cmplx_vec)},
       std::pair{"delim", self.delim},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::ParserFastRealRead& self) {
  return repr(
      &self,
      "Bmad::ParserFastRealRead",
      {std::pair{"real_vec", to_string(self.real_vec)},
       std::pair{"delim", self.delim},
       std::pair{"n_real", to_string(self.n_real)},
       std::pair{"is_ok", to_string(self.is_ok)}});
}
std::string to_string(const Bmad::PhotonAbsorptionAndPhaseShift& self) {
  return repr(
      &self,
      "Bmad::PhotonAbsorptionAndPhaseShift",
      {std::pair{"absorption", to_string(self.absorption)},
       std::pair{"phase_shift", to_string(self.phase_shift)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::PhotonReflection& self) {
  return repr(
      &self,
      "Bmad::PhotonReflection",
      {std::pair{"graze_angle_out", to_string(self.graze_angle_out)},
       std::pair{"phi_out", to_string(self.phi_out)}});
}
std::string to_string(const Bmad::PhotonReflectivity& self) {
  return repr(
      &self,
      "Bmad::PhotonReflectivity",
      {std::pair{"p_reflect", to_string(self.p_reflect)},
       std::pair{"rel_p_specular", to_string(self.rel_p_specular)}});
}
std::string to_string(const Bmad::PointerToElementAtS& self) {
  return repr(
      &self,
      "Bmad::PointerToElementAtS",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"s_eff", to_string(self.s_eff)},
       std::pair{"position", to_string(self.position)},
       std::pair{"ele", to_string(self.ele)}});
}
std::string to_string(const Bmad::PointerToLord& self) {
  return repr(
      &self,
      "Bmad::PointerToLord",
      {std::pair{"control", to_string(self.control)},
       std::pair{"ix_slave_back", to_string(self.ix_slave_back)},
       std::pair{"ix_control", to_string(self.ix_control)},
       std::pair{"ix_ic", to_string(self.ix_ic)}});
}
std::string to_string(const Bmad::PointerToMultipassLord& self) {
  return repr(
      &self,
      "Bmad::PointerToMultipassLord",
      {std::pair{"ix_pass", to_string(self.ix_pass)},
       std::pair{"super_lord", to_string(self.super_lord)}});
}
std::string to_string(const Bmad::PointerToSlave& self) {
  return repr(
      &self,
      "Bmad::PointerToSlave",
      {std::pair{"control", to_string(self.control)},
       std::pair{"ix_lord_back", to_string(self.ix_lord_back)},
       std::pair{"ix_control", to_string(self.ix_control)},
       std::pair{"ix_ic", to_string(self.ix_ic)},
       std::pair{"slave_ptr", to_string(self.slave_ptr)}});
}
std::string to_string(const Bmad::PointerToSuperLord& self) {
  return repr(
      &self,
      "Bmad::PointerToSuperLord",
      {std::pair{"control", to_string(self.control)},
       std::pair{"ix_slave_back", to_string(self.ix_slave_back)},
       std::pair{"ix_control", to_string(self.ix_control)},
       std::pair{"ix_ic", to_string(self.ix_ic)}});
}
std::string to_string(const Bmad::PointerToWall3d& self) {
  return repr(
      &self,
      "Bmad::PointerToWall3d",
      {std::pair{"ds_offset", to_string(self.ds_offset)},
       std::pair{"is_branch_wall", to_string(self.is_branch_wall)},
       std::pair{"wall3d", to_string(self.wall3d)}});
}
std::string to_string(const Bmad::ProjectEmitToXyz& self) {
  return repr(
      &self,
      "Bmad::ProjectEmitToXyz",
      {std::pair{"sigma_x", to_string(self.sigma_x)},
       std::pair{"sigma_y", to_string(self.sigma_y)},
       std::pair{"sigma_z", to_string(self.sigma_z)}});
}
std::string to_string(const Bmad::PtcEmitCalc& self) {
  return repr(
      &self,
      "Bmad::PtcEmitCalc",
      {std::pair{"norm_mode", to_string(self.norm_mode)},
       std::pair{"closed_orb", to_string(self.closed_orb)}});
}
std::string to_string(const Bmad::PtcSpinCalc& self) {
  return repr(
      &self,
      "Bmad::PtcSpinCalc",
      {std::pair{"norm_mode", to_string(self.norm_mode)},
       std::pair{"closed_orb", to_string(self.closed_orb)}});
}
std::string to_string(const Bmad::PtcTrackAll& self) {
  return repr(
      &self,
      "Bmad::PtcTrackAll",
      {std::pair{"track_state", to_string(self.track_state)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const SimUtils::QuatToAxisAngle& self) {
  return repr(
      &self,
      "SimUtils::QuatToAxisAngle",
      {std::pair{"axis", to_string(self.axis)},
       std::pair{"angle", to_string(self.angle)}});
}
std::string to_string(const Bmad::Rad1DampAndStocMats& self) {
  return repr(
      &self,
      "Bmad::Rad1DampAndStocMats",
      {std::pair{"rad_map", to_string(self.rad_map)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"rad_int1", to_string(self.rad_int1)}});
}
std::string to_string(const Bmad::RadDampAndStocMats& self) {
  return repr(
      &self,
      "Bmad::RadDampAndStocMats",
      {std::pair{"rmap", to_string(self.rmap)},
       std::pair{"mode", to_string(self.mode)},
       std::pair{"xfer_nodamp_mat", to_string(self.xfer_nodamp_mat)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"rad_int_branch", to_string(self.rad_int_branch)}});
}
std::string to_string(const Bmad::RadiationIntegrals& self) {
  return repr(
      &self,
      "Bmad::RadiationIntegrals",
      {std::pair{"mode", to_string(self.mode)},
       std::pair{"rad_int_by_ele", to_string(self.rad_int_by_ele)}});
}
std::string to_string(const Bmad::ReadBeamAscii& self) {
  return repr(
      &self,
      "Bmad::ReadBeamAscii",
      {std::pair{"beam", to_string(self.beam)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::ReadBeamFile& self) {
  return repr(
      &self,
      "Bmad::ReadBeamFile",
      {std::pair{"beam", to_string(self.beam)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::RemoveConstantTaylor& self) {
  return repr(
      &self,
      "Bmad::RemoveConstantTaylor",
      {std::pair{"taylor_out", "[...]"}, std::pair{"c0", to_string(self.c0)}});
}
std::string to_string(const Bmad::SetEleAttribute& self) {
  return repr(
      &self,
      "Bmad::SetEleAttribute",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"err_id", to_string(self.err_id)}});
}
std::string to_string(const Bmad::SetEleStatusStale& self) {
  return repr(
      &self,
      "Bmad::SetEleStatusStale",
      {std::pair{"ele", to_string(self.ele)},
       std::pair{"status_group", to_string(self.status_group)},
       std::pair{"set_slaves", to_string(self.set_slaves)}});
}
std::string to_string(const Bmad::SolvePsiFixedSteps& self) {
  return repr(
      &self,
      "Bmad::SolvePsiFixedSteps",
      {std::pair{"t", to_string(self.t)}, std::pair{"p", to_string(self.p)}});
}
std::string to_string(const Bmad::SpinMat8ResonanceStrengths& self) {
  return repr(
      &self,
      "Bmad::SpinMat8ResonanceStrengths",
      {std::pair{"xi_sum", to_string(self.xi_sum)},
       std::pair{"xi_diff", to_string(self.xi_diff)}});
}
std::string to_string(const Bmad::SpinMatToEigen& self) {
  return repr(
      &self,
      "Bmad::SpinMatToEigen",
      {std::pair{"orb_eval", to_string(self.orb_eval)},
       std::pair{"orb_evec", to_string(self.orb_evec)},
       std::pair{"n0", to_string(self.n0)},
       std::pair{"spin_evec", to_string(self.spin_evec)},
       std::pair{"error", to_string(self.error)}});
}
std::string to_string(const Bmad::SpinQuatResonanceStrengths& self) {
  return repr(
      &self,
      "Bmad::SpinQuatResonanceStrengths",
      {std::pair{"xi_sum", to_string(self.xi_sum)},
       std::pair{"xi_diff", to_string(self.xi_diff)}});
}
std::string to_string(const SimUtils::SplineAkimaInterpolate& self) {
  return repr(
      &self,
      "SimUtils::SplineAkimaInterpolate",
      {std::pair{"ok", to_string(self.ok)},
       std::pair{"y", to_string(self.y)},
       std::pair{"dy", to_string(self.dy)}});
}
std::string to_string(const SimUtils::SplineEvaluate& self) {
  return repr(
      &self,
      "SimUtils::SplineEvaluate",
      {std::pair{"ok", to_string(self.ok)},
       std::pair{"y", to_string(self.y)},
       std::pair{"dy", to_string(self.dy)}});
}
std::string to_string(const Bmad::SplitLat& self) {
  return repr(
      &self,
      "Bmad::SplitLat",
      {std::pair{"ix_split", to_string(self.ix_split)},
       std::pair{"split_done", to_string(self.split_done)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::StrongBeamSigmaCalc& self) {
  return repr(
      &self,
      "Bmad::StrongBeamSigmaCalc",
      {std::pair{"sigma", to_string(self.sigma)},
       std::pair{"bbi_const", to_string(self.bbi_const)},
       std::pair{"dsigma_ds", to_string(self.dsigma_ds)}});
}
std::string to_string(const SimUtils::SuperBicubicInterpolation& self) {
  return repr(
      &self,
      "SimUtils::SuperBicubicInterpolation",
      {std::pair{"ansy", to_string(self.ansy)},
       std::pair{"ansy1", to_string(self.ansy1)},
       std::pair{"ansy2", to_string(self.ansy2)}});
}
std::string to_string(const SimUtils::SuperPolint& self) {
  return repr(
      &self,
      "SimUtils::SuperPolint",
      {std::pair{"y", to_string(self.y)}, std::pair{"dy", to_string(self.dy)}});
}
std::string to_string(const Bmad::T6ToB123& self) {
  return repr(
      &self,
      "Bmad::T6ToB123",
      {std::pair{"B1", to_string(self.B1)},
       std::pair{"B2", to_string(self.B2)},
       std::pair{"B3", to_string(self.B3)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Tao::TaoCurveRmsCalc& self) {
  return repr(
      &self,
      "Tao::TaoCurveRmsCalc",
      {std::pair{"rms", to_string(self.rms)},
       std::pair{"mean", to_string(self.mean)}});
}
std::string to_string(const Tao::TaoDataUseitPlotCalc& self) {
  return repr(
      &self,
      "Tao::TaoDataUseitPlotCalc",
      {std::pair{"data", "[...]"},
       std::pair{"most_invalid", self.most_invalid}});
}
std::string to_string(const Tao::TaoDatumIntegrate& self) {
  return repr(
      &self,
      "Tao::TaoDatumIntegrate",
      {std::pair{"valid_value", to_string(self.valid_value)},
       std::pair{"why_invalid", self.why_invalid},
       std::pair{"result", to_string(self.result)}});
}
std::string to_string(const Tao::TaoEleGeometryWithMisalignments& self) {
  return repr(
      &self,
      "Tao::TaoEleGeometryWithMisalignments",
      {std::pair{"valid_value", to_string(self.valid_value)},
       std::pair{"why_invalid", self.why_invalid},
       std::pair{"value", to_string(self.value)}});
}
std::string to_string(const Tao::TaoEleShapeInfo& self) {
  return repr(
      &self,
      "Tao::TaoEleShapeInfo",
      {std::pair{"e_shape", to_string(self.e_shape)},
       std::pair{"label_name", self.label_name}});
}
std::string to_string(const Tao::TaoEvalFloorOrbit& self) {
  return repr(
      &self,
      "Tao::TaoEvalFloorOrbit",
      {std::pair{"valid_value", to_string(self.valid_value)},
       std::pair{"why_invalid", self.why_invalid},
       std::pair{"value", to_string(self.value)}});
}
std::string to_string(const Tao::TaoEvaluateADatum& self) {
  return repr(
      &self,
      "Tao::TaoEvaluateADatum",
      {std::pair{"datum_value", to_string(self.datum_value)},
       std::pair{"valid_value", to_string(self.valid_value)},
       std::pair{"why_invalid", self.why_invalid}});
}
std::string to_string(const Tao::TaoEvaluateDatumAtS& self) {
  return repr(
      &self,
      "Tao::TaoEvaluateDatumAtS",
      {std::pair{"err_str", self.err_str},
       std::pair{"bad_datum", to_string(self.bad_datum)},
       std::pair{"value", to_string(self.value)}});
}
std::string to_string(const Tao::TaoEvaluateLatOrBeamData& self) {
  return repr(
      &self,
      "Tao::TaoEvaluateLatOrBeamData",
      {std::pair{"err", to_string(self.err)},
       std::pair{"values", to_string(self.values)}});
}
std::string to_string(const Tao::TaoFindPlotRegion& self) {
  return repr(
      &self,
      "Tao::TaoFindPlotRegion",
      {std::pair{"err", to_string(self.err)},
       std::pair{"region", to_string(self.region)}});
}
std::string to_string(const Tao::TaoFloorToScreen& self) {
  return repr(
      &self,
      "Tao::TaoFloorToScreen",
      {std::pair{"x_screen", to_string(self.x_screen)},
       std::pair{"y_screen", to_string(self.y_screen)}});
}
std::string to_string(const Tao::TaoGetData& self) {
  return repr(
      &self,
      "Tao::TaoGetData",
      {std::pair{"data_value", to_string(self.data_value)},
       std::pair{"data_weight", to_string(self.data_weight)},
       std::pair{"data_meas_value", to_string(self.data_meas_value)},
       std::pair{"data_ix_dModel", to_string(self.data_ix_dModel)}});
}
std::string to_string(const Tao::TaoGetOptVars& self) {
  return repr(
      &self,
      "Tao::TaoGetOptVars",
      {std::pair{"var_value", to_string(self.var_value)},
       std::pair{"var_step", to_string(self.var_step)},
       std::pair{"var_delta", to_string(self.var_delta)},
       std::pair{"var_weight", to_string(self.var_weight)},
       std::pair{"var_ix", to_string(self.var_ix)},
       std::pair{
           "ignore_if_weight_is_zero",
           to_string(self.ignore_if_weight_is_zero)},
       std::pair{
           "ignore_if_not_limited", to_string(self.ignore_if_not_limited)}});
}
std::string to_string(const Tao::TaoGraphSMinMaxCalc& self) {
  return repr(
      &self,
      "Tao::TaoGraphSMinMaxCalc",
      {std::pair{"s_min", to_string(self.s_min)},
       std::pair{"s_max", to_string(self.s_max)}});
}
std::string to_string(const Tao::TaoInitFindElements& self) {
  return repr(
      &self,
      "Tao::TaoInitFindElements",
      {std::pair{"eles", "[...]"},
       std::pair{"found_one", to_string(self.found_one)}});
}
std::string to_string(const Tao::TaoInjectBeam& self) {
  return repr(
      &self,
      "Tao::TaoInjectBeam",
      {std::pair{"beam", to_string(self.beam)},
       std::pair{"init_ok", to_string(self.init_ok)}});
}
std::string to_string(const Tao::TaoLatticeCalc& self) {
  return repr(
      &self,
      "Tao::TaoLatticeCalc",
      {std::pair{"calc_ok", to_string(self.calc_ok)},
       std::pair{"print_err", to_string(self.print_err)}});
}
std::string to_string(const Tao::TaoLocateAllElements& self) {
  return repr(
      &self,
      "Tao::TaoLocateAllElements",
      {std::pair{"eles", "[...]"}, std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Tao::TaoLocateElements& self) {
  return repr(
      &self,
      "Tao::TaoLocateElements",
      {std::pair{"eles", "[...]"}, std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Tao::TaoParamValueAtS& self) {
  return repr(
      &self,
      "Tao::TaoParamValueAtS",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"why_invalid", self.why_invalid},
       std::pair{"print_err", to_string(self.print_err)},
       std::pair{"bad_datum", to_string(self.bad_datum)}});
}
std::string to_string(const Tao::TaoParseElementParamStr& self) {
  return repr(
      &self,
      "Tao::TaoParseElementParamStr",
      {std::pair{"err", to_string(self.err)},
       std::pair{"uni", self.uni},
       std::pair{"element", self.element},
       std::pair{"parameter", self.parameter},
       std::pair{"where", to_string(self.where)},
       std::pair{"component", self.component}});
}
std::string to_string(const Tao::TaoParticleDataValue& self) {
  return repr(
      &self,
      "Tao::TaoParticleDataValue",
      {std::pair{"value", to_string(self.value)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Tao::TaoPickUniverse& self) {
  return repr(
      &self,
      "Tao::TaoPickUniverse",
      {std::pair{"name_out", self.name_out},
       std::pair{"err", to_string(self.err)},
       std::pair{"ix_uni", to_string(self.ix_uni)},
       std::pair{"explicit_uni", to_string(self.explicit_uni)}});
}
std::string to_string(const Tao::TaoPointerToDatumEle& self) {
  return repr(
      &self,
      "Tao::TaoPointerToDatumEle",
      {std::pair{"valid", to_string(self.valid)},
       std::pair{"why_invalid", self.why_invalid},
       std::pair{"ele", to_string(self.ele)}});
}
std::string to_string(const Tao::TaoPointerToEleShape& self) {
  return repr(
      &self,
      "Tao::TaoPointerToEleShape",
      {std::pair{"dat_var_name", self.dat_var_name},
       std::pair{"dat_var_value", to_string(self.dat_var_value)}});
}
std::string to_string(const Tao::TaoPointerToUniverses& self) {
  return repr(
      &self,
      "Tao::TaoPointerToUniverses",
      {std::pair{"unis", "[...]"},
       std::pair{"err", to_string(self.err)},
       std::pair{"name_out", self.name_out},
       std::pair{"explicit_uni", to_string(self.explicit_uni)}});
}
std::string to_string(const Tao::TaoScaleGraph& self) {
  return repr(
      &self,
      "Tao::TaoScaleGraph",
      {std::pair{"y_range", to_string(self.y_range)},
       std::pair{"y2_range", to_string(self.y2_range)}});
}
std::string to_string(const Tao::TaoSetIntegerValue& self) {
  return repr(
      &self,
      "Tao::TaoSetIntegerValue",
      {std::pair{"var", to_string(self.var)},
       std::pair{"error", to_string(self.error)}});
}
std::string to_string(const Tao::TaoSetLogicalValue& self) {
  return repr(
      &self,
      "Tao::TaoSetLogicalValue",
      {std::pair{"var", to_string(self.var)},
       std::pair{"error", to_string(self.error)}});
}
std::string to_string(const Tao::TaoSetQpAxisStruct& self) {
  return repr(
      &self,
      "Tao::TaoSetQpAxisStruct",
      {std::pair{"error", to_string(self.error)},
       std::pair{"ix_uni", to_string(self.ix_uni)}});
}
std::string to_string(const Tao::TaoSetQpPointStruct& self) {
  return repr(
      &self,
      "Tao::TaoSetQpPointStruct",
      {std::pair{"error", to_string(self.error)},
       std::pair{"ix_uni", to_string(self.ix_uni)}});
}
std::string to_string(const Tao::TaoSetQpRectStruct& self) {
  return repr(
      &self,
      "Tao::TaoSetQpRectStruct",
      {std::pair{"error", to_string(self.error)},
       std::pair{"ix_uni", to_string(self.ix_uni)}});
}
std::string to_string(const Tao::TaoSetRealValue& self) {
  return repr(
      &self,
      "Tao::TaoSetRealValue",
      {std::pair{"var", to_string(self.var)},
       std::pair{"error", to_string(self.error)}});
}
std::string to_string(const Tao::TaoSplitComponent& self) {
  return repr(
      &self,
      "Tao::TaoSplitComponent",
      {std::pair{"comp", "[...]"}, std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Tao::TaoToPhaseAndCouplingReading& self) {
  return repr(
      &self,
      "Tao::TaoToPhaseAndCouplingReading",
      {std::pair{"bpm_data", to_string(self.bpm_data)},
       std::pair{"valid_value", to_string(self.valid_value)}});
}
std::string to_string(const Tao::TaoToReal& self) {
  return repr(
      &self,
      "Tao::TaoToReal",
      {std::pair{"value", to_string(self.value)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Tao::TaoTrackingEleIndex& self) {
  return repr(
      &self,
      "Tao::TaoTrackingEleIndex",
      {std::pair{"ix_branch", to_string(self.ix_branch)},
       std::pair{"ix_ele", to_string(self.ix_ele)}});
}
std::string to_string(const Tao::TaoWaveFit& self) {
  return repr(
      &self,
      "Tao::TaoWaveFit",
      {std::pair{"coef", to_string(self.coef)},
       std::pair{"rms", to_string(self.rms)}});
}
std::string to_string(const Bmad::TargetRotMats& self) {
  return repr(
      &self,
      "Bmad::TargetRotMats",
      {std::pair{"w_to_target", to_string(self.w_to_target)},
       std::pair{"w_to_ele", to_string(self.w_to_ele)}});
}
std::string to_string(const Bmad::TaylorInverse& self) {
  return repr(
      &self,
      "Bmad::TaylorInverse",
      {std::pair{"taylor_inv", "[...]"},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const CppBmadTest::TestBunchStructArray& self) {
  return repr(
      &self,
      "CppBmadTest::TestBunchStructArray",
      {std::pair{"arr_out", "[...]"},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestBunchStructScalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestBunchStructScalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestCharacterScalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestCharacterScalar",
      {std::pair{"val_out", self.val_out},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestComplexArray& self) {
  return repr(
      &self,
      "CppBmadTest::TestComplexArray",
      {std::pair{"arr_out", to_string(self.arr_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestComplexScalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestComplexScalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestInteger8Array& self) {
  return repr(
      &self,
      "CppBmadTest::TestInteger8Array",
      {std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestInteger8Scalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestInteger8Scalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestIntegerArray& self) {
  return repr(
      &self,
      "CppBmadTest::TestIntegerArray",
      {std::pair{"arr_out", to_string(self.arr_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestIntegerScalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestIntegerScalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestLogicalArray& self) {
  return repr(
      &self,
      "CppBmadTest::TestLogicalArray",
      {std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestLogicalScalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestLogicalScalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestReal16Array& self) {
  return repr(
      &self,
      "CppBmadTest::TestReal16Array",
      {std::pair{"arr_out", to_string(self.arr_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestReal16Scalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestReal16Scalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestRealArray& self) {
  return repr(
      &self,
      "CppBmadTest::TestRealArray",
      {std::pair{"arr_out", to_string(self.arr_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const CppBmadTest::TestRealScalar& self) {
  return repr(
      &self,
      "CppBmadTest::TestRealScalar",
      {std::pair{"val_out", to_string(self.val_out)},
       std::pair{"opt_status", to_string(self.opt_status)}});
}
std::string to_string(const Bmad::ToEtaReading& self) {
  return repr(
      &self,
      "Bmad::ToEtaReading",
      {std::pair{"reading", to_string(self.reading)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::ToOrbitReading& self) {
  return repr(
      &self,
      "Bmad::ToOrbitReading",
      {std::pair{"reading", to_string(self.reading)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::ToPhaseAndCouplingReading& self) {
  return repr(
      &self,
      "Bmad::ToPhaseAndCouplingReading",
      {std::pair{"reading", to_string(self.reading)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::Track1& self) {
  return repr(
      &self,
      "Bmad::Track1",
      {std::pair{"end_orb", to_string(self.end_orb)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::Track1Bmad& self) {
  return repr(
      &self,
      "Bmad::Track1Bmad",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"track", to_string(self.track)}});
}
std::string to_string(const Bmad::Track1RungeKutta& self) {
  return repr(
      &self,
      "Bmad::Track1RungeKutta",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"track", to_string(self.track)}});
}
std::string to_string(const Bmad::Track1Spin& self) {
  return repr(
      &self,
      "Bmad::Track1Spin",
      {std::pair{"ele", to_string(self.ele)},
       std::pair{"end_orb", to_string(self.end_orb)}});
}
std::string to_string(const Bmad::Track1TimeRungeKutta& self) {
  return repr(
      &self,
      "Bmad::Track1TimeRungeKutta",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"track", to_string(self.track)}});
}
std::string to_string(const Bmad::TrackABeambeam& self) {
  return repr(
      &self,
      "Bmad::TrackABeambeam",
      {std::pair{"track", to_string(self.track)},
       std::pair{"mat6", to_string(self.mat6)}});
}
std::string to_string(const Bmad::TrackAPatch& self) {
  return repr(
      &self,
      "Bmad::TrackAPatch",
      {std::pair{"s_ent", to_string(self.s_ent)},
       std::pair{"ds_ref", to_string(self.ds_ref)},
       std::pair{"mat6", to_string(self.mat6)}});
}
std::string to_string(const Bmad::TrackAZeroLengthElement& self) {
  return repr(
      &self,
      "Bmad::TrackAZeroLengthElement",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"track", to_string(self.track)}});
}
std::string to_string(const Bmad::TrackAll& self) {
  return repr(
      &self,
      "Bmad::TrackAll",
      {std::pair{"track_state", to_string(self.track_state)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"orbit0", "[...]"}});
}
std::string to_string(const Bmad::TrackFromSToS& self) {
  return repr(
      &self,
      "Bmad::TrackFromSToS",
      {std::pair{"orbit_end", to_string(self.orbit_end)},
       std::pair{"all_orb", "[...]"},
       std::pair{"track_state", to_string(self.track_state)}});
}
std::string to_string(const Bmad::TrackUntilDead& self) {
  return repr(
      &self,
      "Bmad::TrackUntilDead",
      {std::pair{"end_orb", to_string(self.end_orb)},
       std::pair{"track", to_string(self.track)}});
}
std::string to_string(const Bmad::TrackingRadMapSetup& self) {
  return repr(
      &self,
      "Bmad::TrackingRadMapSetup",
      {std::pair{"rad_map", to_string(self.rad_map)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::TransferMapFromSToS& self) {
  return repr(
      &self,
      "Bmad::TransferMapFromSToS",
      {std::pair{"ref_orb_out", to_string(self.ref_orb_out)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const SimUtils::TricubicCmplxEval& self) {
  return repr(
      &self,
      "SimUtils::TricubicCmplxEval",
      {std::pair{"df_dx", to_string(self.df_dx)},
       std::pair{"df_dy", to_string(self.df_dy)},
       std::pair{"df_dz", to_string(self.df_dz)},
       std::pair{"f_val", to_string(self.f_val)}});
}
std::string to_string(const Bmad::Twiss1Propagate& self) {
  return repr(
      &self,
      "Bmad::Twiss1Propagate",
      {std::pair{"twiss2", to_string(self.twiss2)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::TwissAndTrackFromSToS& self) {
  return repr(
      &self,
      "Bmad::TwissAndTrackFromSToS",
      {std::pair{"orbit_end", to_string(self.orbit_end)},
       std::pair{"ele_end", to_string(self.ele_end)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::TwissAndTrackIntraEle& self) {
  return repr(
      &self,
      "Bmad::TwissAndTrackIntraEle",
      {std::pair{"orbit_end", to_string(self.orbit_end)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::TwissAtElement& self) {
  return repr(
      &self,
      "Bmad::TwissAtElement",
      {std::pair{"start", to_string(self.start)},
       std::pair{"end", to_string(self.end)},
       std::pair{"average", to_string(self.average)}});
}
std::string to_string(const Bmad::TwissFromTracking& self) {
  return repr(
      &self,
      "Bmad::TwissFromTracking",
      {std::pair{"symp_err", to_string(self.symp_err)},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const SimUtils::WMatToAxisAngle& self) {
  return repr(
      &self,
      "SimUtils::WMatToAxisAngle",
      {std::pair{"axis", to_string(self.axis)},
       std::pair{"angle", to_string(self.angle)}});
}
std::string to_string(const Bmad::Wall3dDRadius& self) {
  return repr(
      &self,
      "Bmad::Wall3dDRadius",
      {std::pair{"perp", to_string(self.perp)},
       std::pair{"ix_section", to_string(self.ix_section)},
       std::pair{"no_wall_here", to_string(self.no_wall_here)},
       std::pair{"origin", to_string(self.origin)},
       std::pair{"radius_wall", to_string(self.radius_wall)},
       std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"d_radius", to_string(self.d_radius)}});
}
std::string to_string(const Bmad::WriteAstraFieldGridFile& self) {
  return repr(
      &self,
      "Bmad::WriteAstraFieldGridFile",
      {std::pair{"maxfield", to_string(self.maxfield)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::WriteAstraFieldGridFile3d& self) {
  return repr(
      &self,
      "Bmad::WriteAstraFieldGridFile3d",
      {std::pair{"maxfield", to_string(self.maxfield)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::WriteGptFieldGridFile1d& self) {
  return repr(
      &self,
      "Bmad::WriteGptFieldGridFile1d",
      {std::pair{"maxfield", to_string(self.maxfield)},
       std::pair{"ref_time", to_string(self.ref_time)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::WriteGptFieldGridFile2d& self) {
  return repr(
      &self,
      "Bmad::WriteGptFieldGridFile2d",
      {std::pair{"maxfield", to_string(self.maxfield)},
       std::pair{"ref_time", to_string(self.ref_time)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::WriteGptFieldGridFile3d& self) {
  return repr(
      &self,
      "Bmad::WriteGptFieldGridFile3d",
      {std::pair{"maxfield", to_string(self.maxfield)},
       std::pair{"ref_time", to_string(self.ref_time)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::WriteLatticeInScibmad& self) {
  return repr(
      &self,
      "Bmad::WriteLatticeInScibmad",
      {std::pair{"scibmad_file", self.scibmad_file},
       std::pair{"err_flag", to_string(self.err_flag)}});
}
std::string to_string(const Bmad::WriteOpalFieldGridFile& self) {
  return repr(
      &self,
      "Bmad::WriteOpalFieldGridFile",
      {std::pair{"maxfield", to_string(self.maxfield)},
       std::pair{"err", to_string(self.err)}});
}
std::string to_string(const Bmad::ZAtSurface& self) {
  return repr(
      &self,
      "Bmad::ZAtSurface",
      {std::pair{"err_flag", to_string(self.err_flag)},
       std::pair{"dz_dxy", to_string(self.dz_dxy)},
       std::pair{"z", to_string(self.z)}});
}
} // namespace Bmad
