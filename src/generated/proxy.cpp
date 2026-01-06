#include "bmad/generated/proxy.hpp"

using namespace Bmad;
double SplineProxy::x0() const {
  double value;
  spline_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SplineProxy::set_x0(double value) {
  spline_struct_set_x0(fortran_ptr_, value);
}
double SplineProxy::y0() const {
  double value;
  spline_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SplineProxy::set_y0(double value) {
  spline_struct_set_y0(fortran_ptr_, value);
}
double SplineProxy::x1() const {
  double value;
  spline_struct_get_x1(fortran_ptr_, &value);
  return value;
}
void SplineProxy::set_x1(double value) {
  spline_struct_set_x1(fortran_ptr_, value);
}
FArray1D<double> SplineProxy::coef() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, spline_struct_get_coef_info);
}
double SpinPolarProxy::polarization() const {
  double value;
  spin_polar_struct_get_polarization(fortran_ptr_, &value);
  return value;
}
void SpinPolarProxy::set_polarization(double value) {
  spin_polar_struct_set_polarization(fortran_ptr_, value);
}
double SpinPolarProxy::theta() const {
  double value;
  spin_polar_struct_get_theta(fortran_ptr_, &value);
  return value;
}
void SpinPolarProxy::set_theta(double value) {
  spin_polar_struct_set_theta(fortran_ptr_, value);
}
double SpinPolarProxy::phi() const {
  double value;
  spin_polar_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void SpinPolarProxy::set_phi(double value) {
  spin_polar_struct_set_phi(fortran_ptr_, value);
}
double SpinPolarProxy::xi() const {
  double value;
  spin_polar_struct_get_xi(fortran_ptr_, &value);
  return value;
}
void SpinPolarProxy::set_xi(double value) {
  spin_polar_struct_set_xi(fortran_ptr_, value);
}
double AcKickerTimeProxy::amp() const {
  double value;
  ac_kicker_time_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void AcKickerTimeProxy::set_amp(double value) {
  ac_kicker_time_struct_set_amp(fortran_ptr_, value);
}
double AcKickerTimeProxy::time() const {
  double value;
  ac_kicker_time_struct_get_time(fortran_ptr_, &value);
  return value;
}
void AcKickerTimeProxy::set_time(double value) {
  ac_kicker_time_struct_set_time(fortran_ptr_, value);
}
SplineProxy AcKickerTimeProxy::spline() const {
  void* ptr;
  ac_kicker_time_struct_get_spline(fortran_ptr_, &ptr);
  return SplineProxy(ptr);
}
void AcKickerTimeProxy::set_spline(const SplineProxy& src) {
  ac_kicker_time_struct_set_spline(fortran_ptr_, src.get_fortran_ptr());
}
double AcKickerFreqProxy::f() const {
  double value;
  ac_kicker_freq_struct_get_f(fortran_ptr_, &value);
  return value;
}
void AcKickerFreqProxy::set_f(double value) {
  ac_kicker_freq_struct_set_f(fortran_ptr_, value);
}
double AcKickerFreqProxy::amp() const {
  double value;
  ac_kicker_freq_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void AcKickerFreqProxy::set_amp(double value) {
  ac_kicker_freq_struct_set_amp(fortran_ptr_, value);
}
double AcKickerFreqProxy::phi() const {
  double value;
  ac_kicker_freq_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void AcKickerFreqProxy::set_phi(double value) {
  ac_kicker_freq_struct_set_phi(fortran_ptr_, value);
}
AcKickerTimeProxyArray1D AcKickerProxy::amp_vs_time() const {
  return ProxyHelpers::get_type_array_1d<AcKickerTimeProxyArray1D>(
      fortran_ptr_, ac_kicker_struct_get_amp_vs_time_info);
}
AcKickerFreqProxyArray1D AcKickerProxy::frequency() const {
  return ProxyHelpers::get_type_array_1d<AcKickerFreqProxyArray1D>(
      fortran_ptr_, ac_kicker_struct_get_frequency_info);
}
double Interval1CoefProxy::c0() const {
  double value;
  interval1_coef_struct_get_c0(fortran_ptr_, &value);
  return value;
}
void Interval1CoefProxy::set_c0(double value) {
  interval1_coef_struct_set_c0(fortran_ptr_, value);
}
double Interval1CoefProxy::c1() const {
  double value;
  interval1_coef_struct_get_c1(fortran_ptr_, &value);
  return value;
}
void Interval1CoefProxy::set_c1(double value) {
  interval1_coef_struct_set_c1(fortran_ptr_, value);
}
double Interval1CoefProxy::n_exp() const {
  double value;
  interval1_coef_struct_get_n_exp(fortran_ptr_, &value);
  return value;
}
void Interval1CoefProxy::set_n_exp(double value) {
  interval1_coef_struct_set_n_exp(fortran_ptr_, value);
}
FArray1D<double> PhotonReflectTableProxy::angle() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_reflect_table_struct_get_angle_info);
}
FArray1D<double> PhotonReflectTableProxy::energy() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_reflect_table_struct_get_energy_info);
}
Interval1CoefProxyArray1D PhotonReflectTableProxy::int1() const {
  return ProxyHelpers::get_type_array_1d<Interval1CoefProxyArray1D>(
      fortran_ptr_, photon_reflect_table_struct_get_int1_info);
}
FArray2D<double> PhotonReflectTableProxy::p_reflect() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, photon_reflect_table_struct_get_p_reflect_info);
}
double PhotonReflectTableProxy::max_energy() const {
  double value;
  photon_reflect_table_struct_get_max_energy(fortran_ptr_, &value);
  return value;
}
void PhotonReflectTableProxy::set_max_energy(double value) {
  photon_reflect_table_struct_set_max_energy(fortran_ptr_, value);
}
FArray1D<double> PhotonReflectTableProxy::p_reflect_scratch() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_reflect_table_struct_get_p_reflect_scratch_info);
}
FArray1D<double> PhotonReflectTableProxy::bragg_angle() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_reflect_table_struct_get_bragg_angle_info);
}
std::string PhotonReflectSurfaceProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, photon_reflect_surface_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void PhotonReflectSurfaceProxy::set_name(const std::string& value) {
  photon_reflect_surface_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string PhotonReflectSurfaceProxy::description() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, photon_reflect_surface_struct_get_description_info);
  return std::string(arr.data(), arr.size());
}
void PhotonReflectSurfaceProxy::set_description(const std::string& value) {
  photon_reflect_surface_struct_set_description(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string PhotonReflectSurfaceProxy::reflectivity_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, photon_reflect_surface_struct_get_reflectivity_file_info);
  return std::string(arr.data(), arr.size());
}
void PhotonReflectSurfaceProxy::set_reflectivity_file(
    const std::string& value) {
  photon_reflect_surface_struct_set_reflectivity_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
PhotonReflectTableProxyArray1D PhotonReflectSurfaceProxy::table() const {
  return ProxyHelpers::get_type_array_1d<PhotonReflectTableProxyArray1D>(
      fortran_ptr_, photon_reflect_surface_struct_get_table_info);
}
double PhotonReflectSurfaceProxy::surface_roughness_rms() const {
  double value;
  photon_reflect_surface_struct_get_surface_roughness_rms(fortran_ptr_, &value);
  return value;
}
void PhotonReflectSurfaceProxy::set_surface_roughness_rms(double value) {
  photon_reflect_surface_struct_set_surface_roughness_rms(fortran_ptr_, value);
}
double PhotonReflectSurfaceProxy::roughness_correlation_len() const {
  double value;
  photon_reflect_surface_struct_get_roughness_correlation_len(
      fortran_ptr_, &value);
  return value;
}
void PhotonReflectSurfaceProxy::set_roughness_correlation_len(double value) {
  photon_reflect_surface_struct_set_roughness_correlation_len(
      fortran_ptr_, value);
}
int PhotonReflectSurfaceProxy::ix_surface() const {
  int value;
  photon_reflect_surface_struct_get_ix_surface(fortran_ptr_, &value);
  return value;
}
void PhotonReflectSurfaceProxy::set_ix_surface(int value) {
  photon_reflect_surface_struct_set_ix_surface(fortran_ptr_, value);
}
FArray1D<double> CoordProxy::vec() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, coord_struct_get_vec_info);
}
double CoordProxy::s() const {
  double value;
  coord_struct_get_s(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_s(double value) {
  coord_struct_set_s(fortran_ptr_, value);
}
long double CoordProxy::t() const {
  long double value;
  coord_struct_get_t(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_t(long double value) {
  coord_struct_set_t(fortran_ptr_, value);
}
FArray1D<double> CoordProxy::spin() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, coord_struct_get_spin_info);
}
FArray1D<double> CoordProxy::field() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, coord_struct_get_field_info);
}
FArray1D<double> CoordProxy::phase() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, coord_struct_get_phase_info);
}
double CoordProxy::charge() const {
  double value;
  coord_struct_get_charge(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_charge(double value) {
  coord_struct_set_charge(fortran_ptr_, value);
}
double CoordProxy::dt_ref() const {
  double value;
  coord_struct_get_dt_ref(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_dt_ref(double value) {
  coord_struct_set_dt_ref(fortran_ptr_, value);
}
double CoordProxy::r() const {
  double value;
  coord_struct_get_r(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_r(double value) {
  coord_struct_set_r(fortran_ptr_, value);
}
double CoordProxy::p0c() const {
  double value;
  coord_struct_get_p0c(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_p0c(double value) {
  coord_struct_set_p0c(fortran_ptr_, value);
}
double CoordProxy::E_potential() const {
  double value;
  coord_struct_get_E_potential(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_E_potential(double value) {
  coord_struct_set_E_potential(fortran_ptr_, value);
}
double CoordProxy::beta() const {
  double value;
  coord_struct_get_beta(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_beta(double value) {
  coord_struct_set_beta(fortran_ptr_, value);
}
int CoordProxy::ix_ele() const {
  int value;
  coord_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_ix_ele(int value) {
  coord_struct_set_ix_ele(fortran_ptr_, value);
}
int CoordProxy::ix_branch() const {
  int value;
  coord_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_ix_branch(int value) {
  coord_struct_set_ix_branch(fortran_ptr_, value);
}
int CoordProxy::ix_turn() const {
  int value;
  coord_struct_get_ix_turn(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_ix_turn(int value) {
  coord_struct_set_ix_turn(fortran_ptr_, value);
}
int CoordProxy::ix_user() const {
  int value;
  coord_struct_get_ix_user(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_ix_user(int value) {
  coord_struct_set_ix_user(fortran_ptr_, value);
}
int CoordProxy::state() const {
  int value;
  coord_struct_get_state(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_state(int value) {
  coord_struct_set_state(fortran_ptr_, value);
}
int CoordProxy::direction() const {
  int value;
  coord_struct_get_direction(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_direction(int value) {
  coord_struct_set_direction(fortran_ptr_, value);
}
int CoordProxy::time_dir() const {
  int value;
  coord_struct_get_time_dir(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_time_dir(int value) {
  coord_struct_set_time_dir(fortran_ptr_, value);
}
int CoordProxy::species() const {
  int value;
  coord_struct_get_species(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_species(int value) {
  coord_struct_set_species(fortran_ptr_, value);
}
int CoordProxy::location() const {
  int value;
  coord_struct_get_location(fortran_ptr_, &value);
  return value;
}
void CoordProxy::set_location(int value) {
  coord_struct_set_location(fortran_ptr_, value);
}
CoordProxyArray1D CoordArrayProxy::orbit() const {
  return ProxyHelpers::get_type_array_1d<CoordProxyArray1D>(
      fortran_ptr_, coord_array_struct_get_orbit_info);
}
double BpmPhaseCouplingProxy::K_22a() const {
  double value;
  bpm_phase_coupling_struct_get_K_22a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_K_22a(double value) {
  bpm_phase_coupling_struct_set_K_22a(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::K_12a() const {
  double value;
  bpm_phase_coupling_struct_get_K_12a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_K_12a(double value) {
  bpm_phase_coupling_struct_set_K_12a(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::K_11b() const {
  double value;
  bpm_phase_coupling_struct_get_K_11b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_K_11b(double value) {
  bpm_phase_coupling_struct_set_K_11b(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::K_12b() const {
  double value;
  bpm_phase_coupling_struct_get_K_12b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_K_12b(double value) {
  bpm_phase_coupling_struct_set_K_12b(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::Cbar22_a() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar22_a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_Cbar22_a(double value) {
  bpm_phase_coupling_struct_set_Cbar22_a(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::Cbar12_a() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar12_a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_Cbar12_a(double value) {
  bpm_phase_coupling_struct_set_Cbar12_a(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::Cbar11_b() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar11_b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_Cbar11_b(double value) {
  bpm_phase_coupling_struct_set_Cbar11_b(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::Cbar12_b() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar12_b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_Cbar12_b(double value) {
  bpm_phase_coupling_struct_set_Cbar12_b(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::phi_a() const {
  double value;
  bpm_phase_coupling_struct_get_phi_a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_phi_a(double value) {
  bpm_phase_coupling_struct_set_phi_a(fortran_ptr_, value);
}
double BpmPhaseCouplingProxy::phi_b() const {
  double value;
  bpm_phase_coupling_struct_get_phi_b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingProxy::set_phi_b(double value) {
  bpm_phase_coupling_struct_set_phi_b(fortran_ptr_, value);
}
std::string ExpressionAtomProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, expression_atom_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void ExpressionAtomProxy::set_name(const std::string& value) {
  expression_atom_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int ExpressionAtomProxy::type() const {
  int value;
  expression_atom_struct_get_type(fortran_ptr_, &value);
  return value;
}
void ExpressionAtomProxy::set_type(int value) {
  expression_atom_struct_set_type(fortran_ptr_, value);
}
double ExpressionAtomProxy::value() const {
  double value;
  expression_atom_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ExpressionAtomProxy::set_value(double value) {
  expression_atom_struct_set_value(fortran_ptr_, value);
}
FArray1D<double> WakeSrZLongProxy::w() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, wake_sr_z_long_struct_get_w_info);
}
FArray1D<std::complex<double>> WakeSrZLongProxy::fw() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, wake_sr_z_long_struct_get_fw_info);
}
FArray1D<std::complex<double>> WakeSrZLongProxy::fbunch() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, wake_sr_z_long_struct_get_fbunch_info);
}
FArray1D<std::complex<double>> WakeSrZLongProxy::w_out() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, wake_sr_z_long_struct_get_w_out_info);
}
double WakeSrZLongProxy::dz() const {
  double value;
  wake_sr_z_long_struct_get_dz(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongProxy::set_dz(double value) {
  wake_sr_z_long_struct_set_dz(fortran_ptr_, value);
}
double WakeSrZLongProxy::z0() const {
  double value;
  wake_sr_z_long_struct_get_z0(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongProxy::set_z0(double value) {
  wake_sr_z_long_struct_set_z0(fortran_ptr_, value);
}
double WakeSrZLongProxy::smoothing_sigma() const {
  double value;
  wake_sr_z_long_struct_get_smoothing_sigma(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongProxy::set_smoothing_sigma(double value) {
  wake_sr_z_long_struct_set_smoothing_sigma(fortran_ptr_, value);
}
int WakeSrZLongProxy::position_dependence() const {
  int value;
  wake_sr_z_long_struct_get_position_dependence(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongProxy::set_position_dependence(int value) {
  wake_sr_z_long_struct_set_position_dependence(fortran_ptr_, value);
}
bool WakeSrZLongProxy::time_based() const {
  bool value;
  wake_sr_z_long_struct_get_time_based(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongProxy::set_time_based(bool value) {
  wake_sr_z_long_struct_set_time_based(fortran_ptr_, value);
}
double WakeSrModeProxy::amp() const {
  double value;
  wake_sr_mode_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_amp(double value) {
  wake_sr_mode_struct_set_amp(fortran_ptr_, value);
}
double WakeSrModeProxy::damp() const {
  double value;
  wake_sr_mode_struct_get_damp(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_damp(double value) {
  wake_sr_mode_struct_set_damp(fortran_ptr_, value);
}
double WakeSrModeProxy::k() const {
  double value;
  wake_sr_mode_struct_get_k(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_k(double value) {
  wake_sr_mode_struct_set_k(fortran_ptr_, value);
}
double WakeSrModeProxy::phi() const {
  double value;
  wake_sr_mode_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_phi(double value) {
  wake_sr_mode_struct_set_phi(fortran_ptr_, value);
}
double WakeSrModeProxy::b_sin() const {
  double value;
  wake_sr_mode_struct_get_b_sin(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_b_sin(double value) {
  wake_sr_mode_struct_set_b_sin(fortran_ptr_, value);
}
double WakeSrModeProxy::b_cos() const {
  double value;
  wake_sr_mode_struct_get_b_cos(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_b_cos(double value) {
  wake_sr_mode_struct_set_b_cos(fortran_ptr_, value);
}
double WakeSrModeProxy::a_sin() const {
  double value;
  wake_sr_mode_struct_get_a_sin(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_a_sin(double value) {
  wake_sr_mode_struct_set_a_sin(fortran_ptr_, value);
}
double WakeSrModeProxy::a_cos() const {
  double value;
  wake_sr_mode_struct_get_a_cos(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_a_cos(double value) {
  wake_sr_mode_struct_set_a_cos(fortran_ptr_, value);
}
int WakeSrModeProxy::polarization() const {
  int value;
  wake_sr_mode_struct_get_polarization(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_polarization(int value) {
  wake_sr_mode_struct_set_polarization(fortran_ptr_, value);
}
int WakeSrModeProxy::position_dependence() const {
  int value;
  wake_sr_mode_struct_get_position_dependence(fortran_ptr_, &value);
  return value;
}
void WakeSrModeProxy::set_position_dependence(int value) {
  wake_sr_mode_struct_set_position_dependence(fortran_ptr_, value);
}
std::string WakeSrProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wake_sr_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void WakeSrProxy::set_file(const std::string& value) {
  wake_sr_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
WakeSrZLongProxy WakeSrProxy::z_long() const {
  void* ptr;
  wake_sr_struct_get_z_long(fortran_ptr_, &ptr);
  return WakeSrZLongProxy(ptr);
}
void WakeSrProxy::set_z_long(const WakeSrZLongProxy& src) {
  wake_sr_struct_set_z_long(fortran_ptr_, src.get_fortran_ptr());
}
WakeSrModeProxyArray1D WakeSrProxy::long_wake() const {
  return ProxyHelpers::get_type_array_1d<WakeSrModeProxyArray1D>(
      fortran_ptr_, wake_sr_struct_get_long_info);
}
WakeSrModeProxyArray1D WakeSrProxy::trans_wake() const {
  return ProxyHelpers::get_type_array_1d<WakeSrModeProxyArray1D>(
      fortran_ptr_, wake_sr_struct_get_trans_info);
}
double WakeSrProxy::z_ref_long() const {
  double value;
  wake_sr_struct_get_z_ref_long(fortran_ptr_, &value);
  return value;
}
void WakeSrProxy::set_z_ref_long(double value) {
  wake_sr_struct_set_z_ref_long(fortran_ptr_, value);
}
double WakeSrProxy::z_ref_trans() const {
  double value;
  wake_sr_struct_get_z_ref_trans(fortran_ptr_, &value);
  return value;
}
void WakeSrProxy::set_z_ref_trans(double value) {
  wake_sr_struct_set_z_ref_trans(fortran_ptr_, value);
}
double WakeSrProxy::z_max() const {
  double value;
  wake_sr_struct_get_z_max(fortran_ptr_, &value);
  return value;
}
void WakeSrProxy::set_z_max(double value) {
  wake_sr_struct_set_z_max(fortran_ptr_, value);
}
double WakeSrProxy::amp_scale() const {
  double value;
  wake_sr_struct_get_amp_scale(fortran_ptr_, &value);
  return value;
}
void WakeSrProxy::set_amp_scale(double value) {
  wake_sr_struct_set_amp_scale(fortran_ptr_, value);
}
double WakeSrProxy::z_scale() const {
  double value;
  wake_sr_struct_get_z_scale(fortran_ptr_, &value);
  return value;
}
void WakeSrProxy::set_z_scale(double value) {
  wake_sr_struct_set_z_scale(fortran_ptr_, value);
}
bool WakeSrProxy::scale_with_length() const {
  bool value;
  wake_sr_struct_get_scale_with_length(fortran_ptr_, &value);
  return value;
}
void WakeSrProxy::set_scale_with_length(bool value) {
  wake_sr_struct_set_scale_with_length(fortran_ptr_, value);
}
double WakeLrModeProxy::freq() const {
  double value;
  wake_lr_mode_struct_get_freq(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_freq(double value) {
  wake_lr_mode_struct_set_freq(fortran_ptr_, value);
}
double WakeLrModeProxy::freq_in() const {
  double value;
  wake_lr_mode_struct_get_freq_in(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_freq_in(double value) {
  wake_lr_mode_struct_set_freq_in(fortran_ptr_, value);
}
double WakeLrModeProxy::R_over_Q() const {
  double value;
  wake_lr_mode_struct_get_R_over_Q(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_R_over_Q(double value) {
  wake_lr_mode_struct_set_R_over_Q(fortran_ptr_, value);
}
double WakeLrModeProxy::Q() const {
  double value;
  wake_lr_mode_struct_get_Q(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_Q(double value) {
  wake_lr_mode_struct_set_Q(fortran_ptr_, value);
}
double WakeLrModeProxy::damp() const {
  double value;
  wake_lr_mode_struct_get_damp(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_damp(double value) {
  wake_lr_mode_struct_set_damp(fortran_ptr_, value);
}
double WakeLrModeProxy::phi() const {
  double value;
  wake_lr_mode_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_phi(double value) {
  wake_lr_mode_struct_set_phi(fortran_ptr_, value);
}
double WakeLrModeProxy::angle() const {
  double value;
  wake_lr_mode_struct_get_angle(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_angle(double value) {
  wake_lr_mode_struct_set_angle(fortran_ptr_, value);
}
double WakeLrModeProxy::b_sin() const {
  double value;
  wake_lr_mode_struct_get_b_sin(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_b_sin(double value) {
  wake_lr_mode_struct_set_b_sin(fortran_ptr_, value);
}
double WakeLrModeProxy::b_cos() const {
  double value;
  wake_lr_mode_struct_get_b_cos(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_b_cos(double value) {
  wake_lr_mode_struct_set_b_cos(fortran_ptr_, value);
}
double WakeLrModeProxy::a_sin() const {
  double value;
  wake_lr_mode_struct_get_a_sin(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_a_sin(double value) {
  wake_lr_mode_struct_set_a_sin(fortran_ptr_, value);
}
double WakeLrModeProxy::a_cos() const {
  double value;
  wake_lr_mode_struct_get_a_cos(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_a_cos(double value) {
  wake_lr_mode_struct_set_a_cos(fortran_ptr_, value);
}
int WakeLrModeProxy::m() const {
  int value;
  wake_lr_mode_struct_get_m(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_m(int value) {
  wake_lr_mode_struct_set_m(fortran_ptr_, value);
}
bool WakeLrModeProxy::polarized() const {
  bool value;
  wake_lr_mode_struct_get_polarized(fortran_ptr_, &value);
  return value;
}
void WakeLrModeProxy::set_polarized(bool value) {
  wake_lr_mode_struct_set_polarized(fortran_ptr_, value);
}
std::string WakeLrProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wake_lr_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void WakeLrProxy::set_file(const std::string& value) {
  wake_lr_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
WakeLrModeProxyArray1D WakeLrProxy::mode() const {
  return ProxyHelpers::get_type_array_1d<WakeLrModeProxyArray1D>(
      fortran_ptr_, wake_lr_struct_get_mode_info);
}
double WakeLrProxy::t_ref() const {
  double value;
  wake_lr_struct_get_t_ref(fortran_ptr_, &value);
  return value;
}
void WakeLrProxy::set_t_ref(double value) {
  wake_lr_struct_set_t_ref(fortran_ptr_, value);
}
double WakeLrProxy::freq_spread() const {
  double value;
  wake_lr_struct_get_freq_spread(fortran_ptr_, &value);
  return value;
}
void WakeLrProxy::set_freq_spread(double value) {
  wake_lr_struct_set_freq_spread(fortran_ptr_, value);
}
double WakeLrProxy::amp_scale() const {
  double value;
  wake_lr_struct_get_amp_scale(fortran_ptr_, &value);
  return value;
}
void WakeLrProxy::set_amp_scale(double value) {
  wake_lr_struct_set_amp_scale(fortran_ptr_, value);
}
double WakeLrProxy::time_scale() const {
  double value;
  wake_lr_struct_get_time_scale(fortran_ptr_, &value);
  return value;
}
void WakeLrProxy::set_time_scale(double value) {
  wake_lr_struct_set_time_scale(fortran_ptr_, value);
}
bool WakeLrProxy::self_wake_on() const {
  bool value;
  wake_lr_struct_get_self_wake_on(fortran_ptr_, &value);
  return value;
}
void WakeLrProxy::set_self_wake_on(bool value) {
  wake_lr_struct_set_self_wake_on(fortran_ptr_, value);
}
int LatEleLocProxy::ix_ele() const {
  int value;
  lat_ele_loc_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void LatEleLocProxy::set_ix_ele(int value) {
  lat_ele_loc_struct_set_ix_ele(fortran_ptr_, value);
}
int LatEleLocProxy::ix_branch() const {
  int value;
  lat_ele_loc_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void LatEleLocProxy::set_ix_branch(int value) {
  lat_ele_loc_struct_set_ix_branch(fortran_ptr_, value);
}
WakeSrProxy WakeProxy::sr() const {
  void* ptr;
  wake_struct_get_sr(fortran_ptr_, &ptr);
  return WakeSrProxy(ptr);
}
void WakeProxy::set_sr(const WakeSrProxy& src) {
  wake_struct_set_sr(fortran_ptr_, src.get_fortran_ptr());
}
WakeLrProxy WakeProxy::lr() const {
  void* ptr;
  wake_struct_get_lr(fortran_ptr_, &ptr);
  return WakeLrProxy(ptr);
}
void WakeProxy::set_lr(const WakeLrProxy& src) {
  wake_struct_set_lr(fortran_ptr_, src.get_fortran_ptr());
}
double TaylorTermProxy::coef() const {
  double value;
  taylor_term_struct_get_coef(fortran_ptr_, &value);
  return value;
}
void TaylorTermProxy::set_coef(double value) {
  taylor_term_struct_set_coef(fortran_ptr_, value);
}
FArray1D<int> TaylorTermProxy::expn() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, taylor_term_struct_get_expn_info);
}
double TaylorProxy::ref() const {
  double value;
  taylor_struct_get_ref(fortran_ptr_, &value);
  return value;
}
void TaylorProxy::set_ref(double value) {
  taylor_struct_set_ref(fortran_ptr_, value);
}
TaylorTermProxyArray1D TaylorProxy::term() const {
  return ProxyHelpers::get_type_array_1d<TaylorTermProxyArray1D>(
      fortran_ptr_, taylor_struct_get_term_info);
}
double EmTaylorTermProxy::coef() const {
  double value;
  em_taylor_term_struct_get_coef(fortran_ptr_, &value);
  return value;
}
void EmTaylorTermProxy::set_coef(double value) {
  em_taylor_term_struct_set_coef(fortran_ptr_, value);
}
FArray1D<int> EmTaylorTermProxy::expn() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, em_taylor_term_struct_get_expn_info);
}
double EmTaylorProxy::ref() const {
  double value;
  em_taylor_struct_get_ref(fortran_ptr_, &value);
  return value;
}
void EmTaylorProxy::set_ref(double value) {
  em_taylor_struct_set_ref(fortran_ptr_, value);
}
EmTaylorTermProxyArray1D EmTaylorProxy::term() const {
  return ProxyHelpers::get_type_array_1d<EmTaylorTermProxyArray1D>(
      fortran_ptr_, em_taylor_struct_get_term_info);
}
double CartesianMapTerm1Proxy::coef() const {
  double value;
  cartesian_map_term1_struct_get_coef(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_coef(double value) {
  cartesian_map_term1_struct_set_coef(fortran_ptr_, value);
}
double CartesianMapTerm1Proxy::kx() const {
  double value;
  cartesian_map_term1_struct_get_kx(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_kx(double value) {
  cartesian_map_term1_struct_set_kx(fortran_ptr_, value);
}
double CartesianMapTerm1Proxy::ky() const {
  double value;
  cartesian_map_term1_struct_get_ky(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_ky(double value) {
  cartesian_map_term1_struct_set_ky(fortran_ptr_, value);
}
double CartesianMapTerm1Proxy::kz() const {
  double value;
  cartesian_map_term1_struct_get_kz(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_kz(double value) {
  cartesian_map_term1_struct_set_kz(fortran_ptr_, value);
}
double CartesianMapTerm1Proxy::x0() const {
  double value;
  cartesian_map_term1_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_x0(double value) {
  cartesian_map_term1_struct_set_x0(fortran_ptr_, value);
}
double CartesianMapTerm1Proxy::y0() const {
  double value;
  cartesian_map_term1_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_y0(double value) {
  cartesian_map_term1_struct_set_y0(fortran_ptr_, value);
}
double CartesianMapTerm1Proxy::phi_z() const {
  double value;
  cartesian_map_term1_struct_get_phi_z(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_phi_z(double value) {
  cartesian_map_term1_struct_set_phi_z(fortran_ptr_, value);
}
int CartesianMapTerm1Proxy::family() const {
  int value;
  cartesian_map_term1_struct_get_family(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_family(int value) {
  cartesian_map_term1_struct_set_family(fortran_ptr_, value);
}
int CartesianMapTerm1Proxy::form() const {
  int value;
  cartesian_map_term1_struct_get_form(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Proxy::set_form(int value) {
  cartesian_map_term1_struct_set_form(fortran_ptr_, value);
}
std::string CartesianMapTermProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, cartesian_map_term_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void CartesianMapTermProxy::set_file(const std::string& value) {
  cartesian_map_term_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int CartesianMapTermProxy::n_link() const {
  int value;
  cartesian_map_term_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void CartesianMapTermProxy::set_n_link(int value) {
  cartesian_map_term_struct_set_n_link(fortran_ptr_, value);
}
CartesianMapTerm1ProxyArray1D CartesianMapTermProxy::term() const {
  return ProxyHelpers::get_type_array_1d<CartesianMapTerm1ProxyArray1D>(
      fortran_ptr_, cartesian_map_term_struct_get_term_info);
}
double CartesianMapProxy::field_scale() const {
  double value;
  cartesian_map_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void CartesianMapProxy::set_field_scale(double value) {
  cartesian_map_struct_set_field_scale(fortran_ptr_, value);
}
FArray1D<double> CartesianMapProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, cartesian_map_struct_get_r0_info);
}
int CartesianMapProxy::master_parameter() const {
  int value;
  cartesian_map_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void CartesianMapProxy::set_master_parameter(int value) {
  cartesian_map_struct_set_master_parameter(fortran_ptr_, value);
}
int CartesianMapProxy::ele_anchor_pt() const {
  int value;
  cartesian_map_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void CartesianMapProxy::set_ele_anchor_pt(int value) {
  cartesian_map_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
int CartesianMapProxy::field_type() const {
  int value;
  cartesian_map_struct_get_field_type(fortran_ptr_, &value);
  return value;
}
void CartesianMapProxy::set_field_type(int value) {
  cartesian_map_struct_set_field_type(fortran_ptr_, value);
}
std::optional<CartesianMapTermProxy> CartesianMapProxy::ptr() const {
  void* ptr;
  cartesian_map_struct_get_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return CartesianMapTermProxy(ptr);
}
void CartesianMapProxy::set_ptr(const CartesianMapTermProxy& src) {
  cartesian_map_struct_set_ptr(fortran_ptr_, src.get_fortran_ptr());
}
std::complex<double> CylindricalMapTerm1Proxy::e_coef() const {
  std::complex<double> c_value;
  cylindrical_map_term1_struct_get_e_coef(fortran_ptr_, &c_value);
  return c_value;
}
void CylindricalMapTerm1Proxy::set_e_coef(std::complex<double> value) {
  cylindrical_map_term1_struct_set_e_coef(fortran_ptr_, value);
}
std::complex<double> CylindricalMapTerm1Proxy::b_coef() const {
  std::complex<double> c_value;
  cylindrical_map_term1_struct_get_b_coef(fortran_ptr_, &c_value);
  return c_value;
}
void CylindricalMapTerm1Proxy::set_b_coef(std::complex<double> value) {
  cylindrical_map_term1_struct_set_b_coef(fortran_ptr_, value);
}
std::string CylindricalMapTermProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, cylindrical_map_term_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void CylindricalMapTermProxy::set_file(const std::string& value) {
  cylindrical_map_term_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int CylindricalMapTermProxy::n_link() const {
  int value;
  cylindrical_map_term_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void CylindricalMapTermProxy::set_n_link(int value) {
  cylindrical_map_term_struct_set_n_link(fortran_ptr_, value);
}
CylindricalMapTerm1ProxyArray1D CylindricalMapTermProxy::term() const {
  return ProxyHelpers::get_type_array_1d<CylindricalMapTerm1ProxyArray1D>(
      fortran_ptr_, cylindrical_map_term_struct_get_term_info);
}
int CylindricalMapProxy::m() const {
  int value;
  cylindrical_map_struct_get_m(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_m(int value) {
  cylindrical_map_struct_set_m(fortran_ptr_, value);
}
int CylindricalMapProxy::harmonic() const {
  int value;
  cylindrical_map_struct_get_harmonic(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_harmonic(int value) {
  cylindrical_map_struct_set_harmonic(fortran_ptr_, value);
}
double CylindricalMapProxy::phi0_fieldmap() const {
  double value;
  cylindrical_map_struct_get_phi0_fieldmap(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_phi0_fieldmap(double value) {
  cylindrical_map_struct_set_phi0_fieldmap(fortran_ptr_, value);
}
double CylindricalMapProxy::theta0_azimuth() const {
  double value;
  cylindrical_map_struct_get_theta0_azimuth(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_theta0_azimuth(double value) {
  cylindrical_map_struct_set_theta0_azimuth(fortran_ptr_, value);
}
double CylindricalMapProxy::field_scale() const {
  double value;
  cylindrical_map_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_field_scale(double value) {
  cylindrical_map_struct_set_field_scale(fortran_ptr_, value);
}
int CylindricalMapProxy::master_parameter() const {
  int value;
  cylindrical_map_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_master_parameter(int value) {
  cylindrical_map_struct_set_master_parameter(fortran_ptr_, value);
}
int CylindricalMapProxy::ele_anchor_pt() const {
  int value;
  cylindrical_map_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_ele_anchor_pt(int value) {
  cylindrical_map_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
double CylindricalMapProxy::dz() const {
  double value;
  cylindrical_map_struct_get_dz(fortran_ptr_, &value);
  return value;
}
void CylindricalMapProxy::set_dz(double value) {
  cylindrical_map_struct_set_dz(fortran_ptr_, value);
}
FArray1D<double> CylindricalMapProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, cylindrical_map_struct_get_r0_info);
}
std::optional<CylindricalMapTermProxy> CylindricalMapProxy::ptr() const {
  void* ptr;
  cylindrical_map_struct_get_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return CylindricalMapTermProxy(ptr);
}
void CylindricalMapProxy::set_ptr(const CylindricalMapTermProxy& src) {
  cylindrical_map_struct_set_ptr(fortran_ptr_, src.get_fortran_ptr());
}
FArray2D<std::complex<double>> BicubicCmplxCoefProxy::coef() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_, bicubic_cmplx_coef_struct_get_coef_info);
}
FArray1D<int> BicubicCmplxCoefProxy::i_box() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, bicubic_cmplx_coef_struct_get_i_box_info);
}
FArray3D<std::complex<double>> TricubicCmplxCoefProxy::coef() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_, tricubic_cmplx_coef_struct_get_coef_info);
}
FArray1D<int> TricubicCmplxCoefProxy::i_box() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, tricubic_cmplx_coef_struct_get_i_box_info);
}
FArray1D<std::complex<double>> GridFieldPt1Proxy::E() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, grid_field_pt1_struct_get_E_info);
}
FArray1D<std::complex<double>> GridFieldPt1Proxy::B() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, grid_field_pt1_struct_get_B_info);
}
std::string GridFieldPtProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, grid_field_pt_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void GridFieldPtProxy::set_file(const std::string& value) {
  grid_field_pt_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int GridFieldPtProxy::n_link() const {
  int value;
  grid_field_pt_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void GridFieldPtProxy::set_n_link(int value) {
  grid_field_pt_struct_set_n_link(fortran_ptr_, value);
}
GridFieldPt1ProxyArray3D GridFieldPtProxy::pt() const {
  return ProxyHelpers::get_type_array_3d<GridFieldPt1ProxyArray3D>(
      fortran_ptr_, grid_field_pt_struct_get_pt_info);
}
int GridFieldProxy::geometry() const {
  int value;
  grid_field_struct_get_geometry(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_geometry(int value) {
  grid_field_struct_set_geometry(fortran_ptr_, value);
}
int GridFieldProxy::harmonic() const {
  int value;
  grid_field_struct_get_harmonic(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_harmonic(int value) {
  grid_field_struct_set_harmonic(fortran_ptr_, value);
}
double GridFieldProxy::phi0_fieldmap() const {
  double value;
  grid_field_struct_get_phi0_fieldmap(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_phi0_fieldmap(double value) {
  grid_field_struct_set_phi0_fieldmap(fortran_ptr_, value);
}
double GridFieldProxy::field_scale() const {
  double value;
  grid_field_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_field_scale(double value) {
  grid_field_struct_set_field_scale(fortran_ptr_, value);
}
int GridFieldProxy::field_type() const {
  int value;
  grid_field_struct_get_field_type(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_field_type(int value) {
  grid_field_struct_set_field_type(fortran_ptr_, value);
}
int GridFieldProxy::master_parameter() const {
  int value;
  grid_field_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_master_parameter(int value) {
  grid_field_struct_set_master_parameter(fortran_ptr_, value);
}
int GridFieldProxy::ele_anchor_pt() const {
  int value;
  grid_field_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_ele_anchor_pt(int value) {
  grid_field_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
int GridFieldProxy::interpolation_order() const {
  int value;
  grid_field_struct_get_interpolation_order(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_interpolation_order(int value) {
  grid_field_struct_set_interpolation_order(fortran_ptr_, value);
}
FArray1D<double> GridFieldProxy::dr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, grid_field_struct_get_dr_info);
}
FArray1D<double> GridFieldProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, grid_field_struct_get_r0_info);
}
bool GridFieldProxy::curved_ref_frame() const {
  bool value;
  grid_field_struct_get_curved_ref_frame(fortran_ptr_, &value);
  return value;
}
void GridFieldProxy::set_curved_ref_frame(bool value) {
  grid_field_struct_set_curved_ref_frame(fortran_ptr_, value);
}
std::optional<GridFieldPtProxy> GridFieldProxy::ptr() const {
  void* ptr;
  grid_field_struct_get_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return GridFieldPtProxy(ptr);
}
void GridFieldProxy::set_ptr(const GridFieldPtProxy& src) {
  grid_field_struct_set_ptr(fortran_ptr_, src.get_fortran_ptr());
}
BicubicCmplxCoefProxyArray3D GridFieldProxy::bi_coef() const {
  return ProxyHelpers::get_type_array_3d<BicubicCmplxCoefProxyArray3D>(
      fortran_ptr_, grid_field_struct_get_bi_coef_info);
}
TricubicCmplxCoefProxyArray3D GridFieldProxy::tri_coef() const {
  return ProxyHelpers::get_type_array_3d<TricubicCmplxCoefProxyArray3D>(
      fortran_ptr_, grid_field_struct_get_tri_coef_info);
}
FArray1D<double> FloorPositionProxy::r() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, floor_position_struct_get_r_info);
}
FArray2D<double> FloorPositionProxy::w() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, floor_position_struct_get_w_info);
}
double FloorPositionProxy::theta() const {
  double value;
  floor_position_struct_get_theta(fortran_ptr_, &value);
  return value;
}
void FloorPositionProxy::set_theta(double value) {
  floor_position_struct_set_theta(fortran_ptr_, value);
}
double FloorPositionProxy::phi() const {
  double value;
  floor_position_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void FloorPositionProxy::set_phi(double value) {
  floor_position_struct_set_phi(fortran_ptr_, value);
}
double FloorPositionProxy::psi() const {
  double value;
  floor_position_struct_get_psi(fortran_ptr_, &value);
  return value;
}
void FloorPositionProxy::set_psi(double value) {
  floor_position_struct_set_psi(fortran_ptr_, value);
}
CoordProxy HighEnergySpaceChargeProxy::closed_orb() const {
  void* ptr;
  high_energy_space_charge_struct_get_closed_orb(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void HighEnergySpaceChargeProxy::set_closed_orb(const CoordProxy& src) {
  high_energy_space_charge_struct_set_closed_orb(
      fortran_ptr_, src.get_fortran_ptr());
}
double HighEnergySpaceChargeProxy::kick_const() const {
  double value;
  high_energy_space_charge_struct_get_kick_const(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_kick_const(double value) {
  high_energy_space_charge_struct_set_kick_const(fortran_ptr_, value);
}
double HighEnergySpaceChargeProxy::sig_x() const {
  double value;
  high_energy_space_charge_struct_get_sig_x(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_sig_x(double value) {
  high_energy_space_charge_struct_set_sig_x(fortran_ptr_, value);
}
double HighEnergySpaceChargeProxy::sig_y() const {
  double value;
  high_energy_space_charge_struct_get_sig_y(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_sig_y(double value) {
  high_energy_space_charge_struct_set_sig_y(fortran_ptr_, value);
}
double HighEnergySpaceChargeProxy::phi() const {
  double value;
  high_energy_space_charge_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_phi(double value) {
  high_energy_space_charge_struct_set_phi(fortran_ptr_, value);
}
double HighEnergySpaceChargeProxy::sin_phi() const {
  double value;
  high_energy_space_charge_struct_get_sin_phi(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_sin_phi(double value) {
  high_energy_space_charge_struct_set_sin_phi(fortran_ptr_, value);
}
double HighEnergySpaceChargeProxy::cos_phi() const {
  double value;
  high_energy_space_charge_struct_get_cos_phi(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_cos_phi(double value) {
  high_energy_space_charge_struct_set_cos_phi(fortran_ptr_, value);
}
double HighEnergySpaceChargeProxy::sig_z() const {
  double value;
  high_energy_space_charge_struct_get_sig_z(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeProxy::set_sig_z(double value) {
  high_energy_space_charge_struct_set_sig_z(fortran_ptr_, value);
}
double XyDispProxy::eta() const {
  double value;
  xy_disp_struct_get_eta(fortran_ptr_, &value);
  return value;
}
void XyDispProxy::set_eta(double value) {
  xy_disp_struct_set_eta(fortran_ptr_, value);
}
double XyDispProxy::etap() const {
  double value;
  xy_disp_struct_get_etap(fortran_ptr_, &value);
  return value;
}
void XyDispProxy::set_etap(double value) {
  xy_disp_struct_set_etap(fortran_ptr_, value);
}
double XyDispProxy::deta_ds() const {
  double value;
  xy_disp_struct_get_deta_ds(fortran_ptr_, &value);
  return value;
}
void XyDispProxy::set_deta_ds(double value) {
  xy_disp_struct_set_deta_ds(fortran_ptr_, value);
}
double XyDispProxy::sigma() const {
  double value;
  xy_disp_struct_get_sigma(fortran_ptr_, &value);
  return value;
}
void XyDispProxy::set_sigma(double value) {
  xy_disp_struct_set_sigma(fortran_ptr_, value);
}
double XyDispProxy::deta_dpz() const {
  double value;
  xy_disp_struct_get_deta_dpz(fortran_ptr_, &value);
  return value;
}
void XyDispProxy::set_deta_dpz(double value) {
  xy_disp_struct_set_deta_dpz(fortran_ptr_, value);
}
double XyDispProxy::detap_dpz() const {
  double value;
  xy_disp_struct_get_detap_dpz(fortran_ptr_, &value);
  return value;
}
void XyDispProxy::set_detap_dpz(double value) {
  xy_disp_struct_set_detap_dpz(fortran_ptr_, value);
}
double TwissProxy::beta() const {
  double value;
  twiss_struct_get_beta(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_beta(double value) {
  twiss_struct_set_beta(fortran_ptr_, value);
}
double TwissProxy::alpha() const {
  double value;
  twiss_struct_get_alpha(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_alpha(double value) {
  twiss_struct_set_alpha(fortran_ptr_, value);
}
double TwissProxy::gamma() const {
  double value;
  twiss_struct_get_gamma(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_gamma(double value) {
  twiss_struct_set_gamma(fortran_ptr_, value);
}
double TwissProxy::phi() const {
  double value;
  twiss_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_phi(double value) {
  twiss_struct_set_phi(fortran_ptr_, value);
}
double TwissProxy::eta() const {
  double value;
  twiss_struct_get_eta(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_eta(double value) {
  twiss_struct_set_eta(fortran_ptr_, value);
}
double TwissProxy::etap() const {
  double value;
  twiss_struct_get_etap(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_etap(double value) {
  twiss_struct_set_etap(fortran_ptr_, value);
}
double TwissProxy::deta_ds() const {
  double value;
  twiss_struct_get_deta_ds(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_deta_ds(double value) {
  twiss_struct_set_deta_ds(fortran_ptr_, value);
}
double TwissProxy::sigma() const {
  double value;
  twiss_struct_get_sigma(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_sigma(double value) {
  twiss_struct_set_sigma(fortran_ptr_, value);
}
double TwissProxy::sigma_p() const {
  double value;
  twiss_struct_get_sigma_p(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_sigma_p(double value) {
  twiss_struct_set_sigma_p(fortran_ptr_, value);
}
double TwissProxy::emit() const {
  double value;
  twiss_struct_get_emit(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_emit(double value) {
  twiss_struct_set_emit(fortran_ptr_, value);
}
double TwissProxy::norm_emit() const {
  double value;
  twiss_struct_get_norm_emit(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_norm_emit(double value) {
  twiss_struct_set_norm_emit(fortran_ptr_, value);
}
double TwissProxy::chrom() const {
  double value;
  twiss_struct_get_chrom(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_chrom(double value) {
  twiss_struct_set_chrom(fortran_ptr_, value);
}
double TwissProxy::dbeta_dpz() const {
  double value;
  twiss_struct_get_dbeta_dpz(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_dbeta_dpz(double value) {
  twiss_struct_set_dbeta_dpz(fortran_ptr_, value);
}
double TwissProxy::dalpha_dpz() const {
  double value;
  twiss_struct_get_dalpha_dpz(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_dalpha_dpz(double value) {
  twiss_struct_set_dalpha_dpz(fortran_ptr_, value);
}
double TwissProxy::deta_dpz() const {
  double value;
  twiss_struct_get_deta_dpz(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_deta_dpz(double value) {
  twiss_struct_set_deta_dpz(fortran_ptr_, value);
}
double TwissProxy::detap_dpz() const {
  double value;
  twiss_struct_get_detap_dpz(fortran_ptr_, &value);
  return value;
}
void TwissProxy::set_detap_dpz(double value) {
  twiss_struct_set_detap_dpz(fortran_ptr_, value);
}
FArray2D<double> Mode3Proxy::v() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, mode3_struct_get_v_info);
}
TwissProxy Mode3Proxy::a() const {
  void* ptr;
  mode3_struct_get_a(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void Mode3Proxy::set_a(const TwissProxy& src) {
  mode3_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy Mode3Proxy::b() const {
  void* ptr;
  mode3_struct_get_b(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void Mode3Proxy::set_b(const TwissProxy& src) {
  mode3_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy Mode3Proxy::c() const {
  void* ptr;
  mode3_struct_get_c(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void Mode3Proxy::set_c(const TwissProxy& src) {
  mode3_struct_set_c(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy Mode3Proxy::x() const {
  void* ptr;
  mode3_struct_get_x(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void Mode3Proxy::set_x(const TwissProxy& src) {
  mode3_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy Mode3Proxy::y() const {
  void* ptr;
  mode3_struct_get_y(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void Mode3Proxy::set_y(const TwissProxy& src) {
  mode3_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
int BookkeepingStateProxy::attributes() const {
  int value;
  bookkeeping_state_struct_get_attributes(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_attributes(int value) {
  bookkeeping_state_struct_set_attributes(fortran_ptr_, value);
}
int BookkeepingStateProxy::control() const {
  int value;
  bookkeeping_state_struct_get_control(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_control(int value) {
  bookkeeping_state_struct_set_control(fortran_ptr_, value);
}
int BookkeepingStateProxy::floor_position() const {
  int value;
  bookkeeping_state_struct_get_floor_position(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_floor_position(int value) {
  bookkeeping_state_struct_set_floor_position(fortran_ptr_, value);
}
int BookkeepingStateProxy::s_position() const {
  int value;
  bookkeeping_state_struct_get_s_position(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_s_position(int value) {
  bookkeeping_state_struct_set_s_position(fortran_ptr_, value);
}
int BookkeepingStateProxy::ref_energy() const {
  int value;
  bookkeeping_state_struct_get_ref_energy(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_ref_energy(int value) {
  bookkeeping_state_struct_set_ref_energy(fortran_ptr_, value);
}
int BookkeepingStateProxy::mat6() const {
  int value;
  bookkeeping_state_struct_get_mat6(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_mat6(int value) {
  bookkeeping_state_struct_set_mat6(fortran_ptr_, value);
}
int BookkeepingStateProxy::rad_int() const {
  int value;
  bookkeeping_state_struct_get_rad_int(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_rad_int(int value) {
  bookkeeping_state_struct_set_rad_int(fortran_ptr_, value);
}
int BookkeepingStateProxy::ptc() const {
  int value;
  bookkeeping_state_struct_get_ptc(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_ptc(int value) {
  bookkeeping_state_struct_set_ptc(fortran_ptr_, value);
}
bool BookkeepingStateProxy::has_misalign() const {
  bool value;
  bookkeeping_state_struct_get_has_misalign(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateProxy::set_has_misalign(bool value) {
  bookkeeping_state_struct_set_has_misalign(fortran_ptr_, value);
}
FArray1D<double> RadMapProxy::ref_orb() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, rad_map_struct_get_ref_orb_info);
}
FArray2D<double> RadMapProxy::damp_dmat() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, rad_map_struct_get_damp_dmat_info);
}
FArray1D<double> RadMapProxy::xfer_damp_vec() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, rad_map_struct_get_xfer_damp_vec_info);
}
FArray2D<double> RadMapProxy::xfer_damp_mat() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, rad_map_struct_get_xfer_damp_mat_info);
}
FArray2D<double> RadMapProxy::stoc_mat() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, rad_map_struct_get_stoc_mat_info);
}
RadMapProxy RadMapEleProxy::rm0() const {
  void* ptr;
  rad_map_ele_struct_get_rm0(fortran_ptr_, &ptr);
  return RadMapProxy(ptr);
}
void RadMapEleProxy::set_rm0(const RadMapProxy& src) {
  rad_map_ele_struct_set_rm0(fortran_ptr_, src.get_fortran_ptr());
}
RadMapProxy RadMapEleProxy::rm1() const {
  void* ptr;
  rad_map_ele_struct_get_rm1(fortran_ptr_, &ptr);
  return RadMapProxy(ptr);
}
void RadMapEleProxy::set_rm1(const RadMapProxy& src) {
  rad_map_ele_struct_set_rm1(fortran_ptr_, src.get_fortran_ptr());
}
bool RadMapEleProxy::stale() const {
  bool value;
  rad_map_ele_struct_get_stale(fortran_ptr_, &value);
  return value;
}
void RadMapEleProxy::set_stale(bool value) {
  rad_map_ele_struct_set_stale(fortran_ptr_, value);
}
int GenGrad1Proxy::m() const {
  int value;
  gen_grad1_struct_get_m(fortran_ptr_, &value);
  return value;
}
void GenGrad1Proxy::set_m(int value) {
  gen_grad1_struct_set_m(fortran_ptr_, value);
}
int GenGrad1Proxy::sincos() const {
  int value;
  gen_grad1_struct_get_sincos(fortran_ptr_, &value);
  return value;
}
void GenGrad1Proxy::set_sincos(int value) {
  gen_grad1_struct_set_sincos(fortran_ptr_, value);
}
int GenGrad1Proxy::n_deriv_max() const {
  int value;
  gen_grad1_struct_get_n_deriv_max(fortran_ptr_, &value);
  return value;
}
void GenGrad1Proxy::set_n_deriv_max(int value) {
  gen_grad1_struct_set_n_deriv_max(fortran_ptr_, value);
}
FArray2D<double> GenGrad1Proxy::deriv() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, gen_grad1_struct_get_deriv_info);
}
std::string GenGradMapProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, gen_grad_map_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void GenGradMapProxy::set_file(const std::string& value) {
  gen_grad_map_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
GenGrad1ProxyArray1D GenGradMapProxy::gg() const {
  return ProxyHelpers::get_type_array_1d<GenGrad1ProxyArray1D>(
      fortran_ptr_, gen_grad_map_struct_get_gg_info);
}
int GenGradMapProxy::ele_anchor_pt() const {
  int value;
  gen_grad_map_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_ele_anchor_pt(int value) {
  gen_grad_map_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
int GenGradMapProxy::field_type() const {
  int value;
  gen_grad_map_struct_get_field_type(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_field_type(int value) {
  gen_grad_map_struct_set_field_type(fortran_ptr_, value);
}
int GenGradMapProxy::iz0() const {
  int value;
  gen_grad_map_struct_get_iz0(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_iz0(int value) {
  gen_grad_map_struct_set_iz0(fortran_ptr_, value);
}
int GenGradMapProxy::iz1() const {
  int value;
  gen_grad_map_struct_get_iz1(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_iz1(int value) {
  gen_grad_map_struct_set_iz1(fortran_ptr_, value);
}
double GenGradMapProxy::dz() const {
  double value;
  gen_grad_map_struct_get_dz(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_dz(double value) {
  gen_grad_map_struct_set_dz(fortran_ptr_, value);
}
FArray1D<double> GenGradMapProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, gen_grad_map_struct_get_r0_info);
}
double GenGradMapProxy::field_scale() const {
  double value;
  gen_grad_map_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_field_scale(double value) {
  gen_grad_map_struct_set_field_scale(fortran_ptr_, value);
}
int GenGradMapProxy::master_parameter() const {
  int value;
  gen_grad_map_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_master_parameter(int value) {
  gen_grad_map_struct_set_master_parameter(fortran_ptr_, value);
}
bool GenGradMapProxy::curved_ref_frame() const {
  bool value;
  gen_grad_map_struct_get_curved_ref_frame(fortran_ptr_, &value);
  return value;
}
void GenGradMapProxy::set_curved_ref_frame(bool value) {
  gen_grad_map_struct_set_curved_ref_frame(fortran_ptr_, value);
}
double SurfaceSegmentedPtProxy::x0() const {
  double value;
  surface_segmented_pt_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtProxy::set_x0(double value) {
  surface_segmented_pt_struct_set_x0(fortran_ptr_, value);
}
double SurfaceSegmentedPtProxy::y0() const {
  double value;
  surface_segmented_pt_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtProxy::set_y0(double value) {
  surface_segmented_pt_struct_set_y0(fortran_ptr_, value);
}
double SurfaceSegmentedPtProxy::z0() const {
  double value;
  surface_segmented_pt_struct_get_z0(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtProxy::set_z0(double value) {
  surface_segmented_pt_struct_set_z0(fortran_ptr_, value);
}
double SurfaceSegmentedPtProxy::dz_dx() const {
  double value;
  surface_segmented_pt_struct_get_dz_dx(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtProxy::set_dz_dx(double value) {
  surface_segmented_pt_struct_set_dz_dx(fortran_ptr_, value);
}
double SurfaceSegmentedPtProxy::dz_dy() const {
  double value;
  surface_segmented_pt_struct_get_dz_dy(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtProxy::set_dz_dy(double value) {
  surface_segmented_pt_struct_set_dz_dy(fortran_ptr_, value);
}
bool SurfaceSegmentedProxy::active() const {
  bool value;
  surface_segmented_struct_get_active(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedProxy::set_active(bool value) {
  surface_segmented_struct_set_active(fortran_ptr_, value);
}
FArray1D<double> SurfaceSegmentedProxy::dr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_segmented_struct_get_dr_info);
}
FArray1D<double> SurfaceSegmentedProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_segmented_struct_get_r0_info);
}
SurfaceSegmentedPtProxyArray2D SurfaceSegmentedProxy::pt() const {
  return ProxyHelpers::get_type_array_2d<SurfaceSegmentedPtProxyArray2D>(
      fortran_ptr_, surface_segmented_struct_get_pt_info);
}
double SurfaceHMisalignPtProxy::x0() const {
  double value;
  surface_h_misalign_pt_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtProxy::set_x0(double value) {
  surface_h_misalign_pt_struct_set_x0(fortran_ptr_, value);
}
double SurfaceHMisalignPtProxy::y0() const {
  double value;
  surface_h_misalign_pt_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtProxy::set_y0(double value) {
  surface_h_misalign_pt_struct_set_y0(fortran_ptr_, value);
}
double SurfaceHMisalignPtProxy::rot_y() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_y(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtProxy::set_rot_y(double value) {
  surface_h_misalign_pt_struct_set_rot_y(fortran_ptr_, value);
}
double SurfaceHMisalignPtProxy::rot_t() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_t(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtProxy::set_rot_t(double value) {
  surface_h_misalign_pt_struct_set_rot_t(fortran_ptr_, value);
}
double SurfaceHMisalignPtProxy::rot_y_rms() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_y_rms(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtProxy::set_rot_y_rms(double value) {
  surface_h_misalign_pt_struct_set_rot_y_rms(fortran_ptr_, value);
}
double SurfaceHMisalignPtProxy::rot_t_rms() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_t_rms(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtProxy::set_rot_t_rms(double value) {
  surface_h_misalign_pt_struct_set_rot_t_rms(fortran_ptr_, value);
}
bool SurfaceHMisalignProxy::active() const {
  bool value;
  surface_h_misalign_struct_get_active(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignProxy::set_active(bool value) {
  surface_h_misalign_struct_set_active(fortran_ptr_, value);
}
FArray1D<double> SurfaceHMisalignProxy::dr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_h_misalign_struct_get_dr_info);
}
FArray1D<double> SurfaceHMisalignProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_h_misalign_struct_get_r0_info);
}
SurfaceHMisalignPtProxyArray2D SurfaceHMisalignProxy::pt() const {
  return ProxyHelpers::get_type_array_2d<SurfaceHMisalignPtProxyArray2D>(
      fortran_ptr_, surface_h_misalign_struct_get_pt_info);
}
double SurfaceDisplacementPtProxy::x0() const {
  double value;
  surface_displacement_pt_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtProxy::set_x0(double value) {
  surface_displacement_pt_struct_set_x0(fortran_ptr_, value);
}
double SurfaceDisplacementPtProxy::y0() const {
  double value;
  surface_displacement_pt_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtProxy::set_y0(double value) {
  surface_displacement_pt_struct_set_y0(fortran_ptr_, value);
}
double SurfaceDisplacementPtProxy::z0() const {
  double value;
  surface_displacement_pt_struct_get_z0(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtProxy::set_z0(double value) {
  surface_displacement_pt_struct_set_z0(fortran_ptr_, value);
}
double SurfaceDisplacementPtProxy::dz_dx() const {
  double value;
  surface_displacement_pt_struct_get_dz_dx(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtProxy::set_dz_dx(double value) {
  surface_displacement_pt_struct_set_dz_dx(fortran_ptr_, value);
}
double SurfaceDisplacementPtProxy::dz_dy() const {
  double value;
  surface_displacement_pt_struct_get_dz_dy(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtProxy::set_dz_dy(double value) {
  surface_displacement_pt_struct_set_dz_dy(fortran_ptr_, value);
}
double SurfaceDisplacementPtProxy::d2z_dxdy() const {
  double value;
  surface_displacement_pt_struct_get_d2z_dxdy(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtProxy::set_d2z_dxdy(double value) {
  surface_displacement_pt_struct_set_d2z_dxdy(fortran_ptr_, value);
}
bool SurfaceDisplacementProxy::active() const {
  bool value;
  surface_displacement_struct_get_active(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementProxy::set_active(bool value) {
  surface_displacement_struct_set_active(fortran_ptr_, value);
}
FArray1D<double> SurfaceDisplacementProxy::dr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_displacement_struct_get_dr_info);
}
FArray1D<double> SurfaceDisplacementProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_displacement_struct_get_r0_info);
}
SurfaceDisplacementPtProxyArray2D SurfaceDisplacementProxy::pt() const {
  return ProxyHelpers::get_type_array_2d<SurfaceDisplacementPtProxyArray2D>(
      fortran_ptr_, surface_displacement_struct_get_pt_info);
}
FArray1D<double> TargetPointProxy::r() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, target_point_struct_get_r_info);
}
FArray2D<double> SurfaceCurvatureProxy::xy() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, surface_curvature_struct_get_xy_info);
}
double SurfaceCurvatureProxy::spherical() const {
  double value;
  surface_curvature_struct_get_spherical(fortran_ptr_, &value);
  return value;
}
void SurfaceCurvatureProxy::set_spherical(double value) {
  surface_curvature_struct_set_spherical(fortran_ptr_, value);
}
FArray1D<double> SurfaceCurvatureProxy::elliptical() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, surface_curvature_struct_get_elliptical_info);
}
bool SurfaceCurvatureProxy::has_curvature() const {
  bool value;
  surface_curvature_struct_get_has_curvature(fortran_ptr_, &value);
  return value;
}
void SurfaceCurvatureProxy::set_has_curvature(bool value) {
  surface_curvature_struct_set_has_curvature(fortran_ptr_, value);
}
int PhotonTargetProxy::type() const {
  int value;
  photon_target_struct_get_type(fortran_ptr_, &value);
  return value;
}
void PhotonTargetProxy::set_type(int value) {
  photon_target_struct_set_type(fortran_ptr_, value);
}
int PhotonTargetProxy::n_corner() const {
  int value;
  photon_target_struct_get_n_corner(fortran_ptr_, &value);
  return value;
}
void PhotonTargetProxy::set_n_corner(int value) {
  photon_target_struct_set_n_corner(fortran_ptr_, value);
}
LatEleLocProxy PhotonTargetProxy::ele_loc() const {
  void* ptr;
  photon_target_struct_get_ele_loc(fortran_ptr_, &ptr);
  return LatEleLocProxy(ptr);
}
void PhotonTargetProxy::set_ele_loc(const LatEleLocProxy& src) {
  photon_target_struct_set_ele_loc(fortran_ptr_, src.get_fortran_ptr());
}
TargetPointProxyArray1D PhotonTargetProxy::corner() const {
  return ProxyHelpers::get_type_array_1d<TargetPointProxyArray1D>(
      fortran_ptr_, photon_target_struct_get_corner_info);
}
TargetPointProxy PhotonTargetProxy::center() const {
  void* ptr;
  photon_target_struct_get_center(fortran_ptr_, &ptr);
  return TargetPointProxy(ptr);
}
void PhotonTargetProxy::set_center(const TargetPointProxy& src) {
  photon_target_struct_set_center(fortran_ptr_, src.get_fortran_ptr());
}
std::complex<double> PhotonMaterialProxy::f0_m1() const {
  std::complex<double> c_value;
  photon_material_struct_get_f0_m1(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialProxy::set_f0_m1(std::complex<double> value) {
  photon_material_struct_set_f0_m1(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialProxy::f0_m2() const {
  std::complex<double> c_value;
  photon_material_struct_get_f0_m2(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialProxy::set_f0_m2(std::complex<double> value) {
  photon_material_struct_set_f0_m2(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialProxy::f_0() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_0(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialProxy::set_f_0(std::complex<double> value) {
  photon_material_struct_set_f_0(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialProxy::f_h() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_h(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialProxy::set_f_h(std::complex<double> value) {
  photon_material_struct_set_f_h(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialProxy::f_hbar() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_hbar(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialProxy::set_f_hbar(std::complex<double> value) {
  photon_material_struct_set_f_hbar(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialProxy::f_hkl() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_hkl(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialProxy::set_f_hkl(std::complex<double> value) {
  photon_material_struct_set_f_hkl(fortran_ptr_, value);
}
FArray1D<double> PhotonMaterialProxy::h_norm() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_material_struct_get_h_norm_info);
}
FArray1D<double> PhotonMaterialProxy::l_ref() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_material_struct_get_l_ref_info);
}
int64_t PixelPtProxy::n_photon() const {
  int64_t value;
  pixel_pt_struct_get_n_photon(fortran_ptr_, &value);
  return value;
}
void PixelPtProxy::set_n_photon(int64_t value) {
  pixel_pt_struct_set_n_photon(fortran_ptr_, value);
}
std::complex<double> PixelPtProxy::E_x() const {
  std::complex<double> c_value;
  pixel_pt_struct_get_E_x(fortran_ptr_, &c_value);
  return c_value;
}
void PixelPtProxy::set_E_x(std::complex<double> value) {
  pixel_pt_struct_set_E_x(fortran_ptr_, value);
}
std::complex<double> PixelPtProxy::E_y() const {
  std::complex<double> c_value;
  pixel_pt_struct_get_E_y(fortran_ptr_, &c_value);
  return c_value;
}
void PixelPtProxy::set_E_y(std::complex<double> value) {
  pixel_pt_struct_set_E_y(fortran_ptr_, value);
}
double PixelPtProxy::intensity_x() const {
  double value;
  pixel_pt_struct_get_intensity_x(fortran_ptr_, &value);
  return value;
}
void PixelPtProxy::set_intensity_x(double value) {
  pixel_pt_struct_set_intensity_x(fortran_ptr_, value);
}
double PixelPtProxy::intensity_y() const {
  double value;
  pixel_pt_struct_get_intensity_y(fortran_ptr_, &value);
  return value;
}
void PixelPtProxy::set_intensity_y(double value) {
  pixel_pt_struct_set_intensity_y(fortran_ptr_, value);
}
double PixelPtProxy::intensity() const {
  double value;
  pixel_pt_struct_get_intensity(fortran_ptr_, &value);
  return value;
}
void PixelPtProxy::set_intensity(double value) {
  pixel_pt_struct_set_intensity(fortran_ptr_, value);
}
FArray1D<double> PixelPtProxy::orbit() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, pixel_pt_struct_get_orbit_info);
}
FArray1D<double> PixelPtProxy::orbit_rms() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, pixel_pt_struct_get_orbit_rms_info);
}
FArray1D<double> PixelPtProxy::init_orbit() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, pixel_pt_struct_get_init_orbit_info);
}
FArray1D<double> PixelPtProxy::init_orbit_rms() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, pixel_pt_struct_get_init_orbit_rms_info);
}
FArray1D<double> PixelDetecProxy::dr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, pixel_detec_struct_get_dr_info);
}
FArray1D<double> PixelDetecProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, pixel_detec_struct_get_r0_info);
}
int64_t PixelDetecProxy::n_track_tot() const {
  int64_t value;
  pixel_detec_struct_get_n_track_tot(fortran_ptr_, &value);
  return value;
}
void PixelDetecProxy::set_n_track_tot(int64_t value) {
  pixel_detec_struct_set_n_track_tot(fortran_ptr_, value);
}
int64_t PixelDetecProxy::n_hit_detec() const {
  int64_t value;
  pixel_detec_struct_get_n_hit_detec(fortran_ptr_, &value);
  return value;
}
void PixelDetecProxy::set_n_hit_detec(int64_t value) {
  pixel_detec_struct_set_n_hit_detec(fortran_ptr_, value);
}
int64_t PixelDetecProxy::n_hit_pixel() const {
  int64_t value;
  pixel_detec_struct_get_n_hit_pixel(fortran_ptr_, &value);
  return value;
}
void PixelDetecProxy::set_n_hit_pixel(int64_t value) {
  pixel_detec_struct_set_n_hit_pixel(fortran_ptr_, value);
}
PixelPtProxyArray2D PixelDetecProxy::pt() const {
  return ProxyHelpers::get_type_array_2d<PixelPtProxyArray2D>(
      fortran_ptr_, pixel_detec_struct_get_pt_info);
}
SurfaceCurvatureProxy PhotonElementProxy::curvature() const {
  void* ptr;
  photon_element_struct_get_curvature(fortran_ptr_, &ptr);
  return SurfaceCurvatureProxy(ptr);
}
void PhotonElementProxy::set_curvature(const SurfaceCurvatureProxy& src) {
  photon_element_struct_set_curvature(fortran_ptr_, src.get_fortran_ptr());
}
PhotonTargetProxy PhotonElementProxy::target() const {
  void* ptr;
  photon_element_struct_get_target(fortran_ptr_, &ptr);
  return PhotonTargetProxy(ptr);
}
void PhotonElementProxy::set_target(const PhotonTargetProxy& src) {
  photon_element_struct_set_target(fortran_ptr_, src.get_fortran_ptr());
}
PhotonMaterialProxy PhotonElementProxy::material() const {
  void* ptr;
  photon_element_struct_get_material(fortran_ptr_, &ptr);
  return PhotonMaterialProxy(ptr);
}
void PhotonElementProxy::set_material(const PhotonMaterialProxy& src) {
  photon_element_struct_set_material(fortran_ptr_, src.get_fortran_ptr());
}
SurfaceSegmentedProxy PhotonElementProxy::segmented() const {
  void* ptr;
  photon_element_struct_get_segmented(fortran_ptr_, &ptr);
  return SurfaceSegmentedProxy(ptr);
}
void PhotonElementProxy::set_segmented(const SurfaceSegmentedProxy& src) {
  photon_element_struct_set_segmented(fortran_ptr_, src.get_fortran_ptr());
}
SurfaceHMisalignProxy PhotonElementProxy::h_misalign() const {
  void* ptr;
  photon_element_struct_get_h_misalign(fortran_ptr_, &ptr);
  return SurfaceHMisalignProxy(ptr);
}
void PhotonElementProxy::set_h_misalign(const SurfaceHMisalignProxy& src) {
  photon_element_struct_set_h_misalign(fortran_ptr_, src.get_fortran_ptr());
}
SurfaceDisplacementProxy PhotonElementProxy::displacement() const {
  void* ptr;
  photon_element_struct_get_displacement(fortran_ptr_, &ptr);
  return SurfaceDisplacementProxy(ptr);
}
void PhotonElementProxy::set_displacement(const SurfaceDisplacementProxy& src) {
  photon_element_struct_set_displacement(fortran_ptr_, src.get_fortran_ptr());
}
PixelDetecProxy PhotonElementProxy::pixel() const {
  void* ptr;
  photon_element_struct_get_pixel(fortran_ptr_, &ptr);
  return PixelDetecProxy(ptr);
}
void PhotonElementProxy::set_pixel(const PixelDetecProxy& src) {
  photon_element_struct_set_pixel(fortran_ptr_, src.get_fortran_ptr());
}
int PhotonElementProxy::reflectivity_table_type() const {
  int value;
  photon_element_struct_get_reflectivity_table_type(fortran_ptr_, &value);
  return value;
}
void PhotonElementProxy::set_reflectivity_table_type(int value) {
  photon_element_struct_set_reflectivity_table_type(fortran_ptr_, value);
}
PhotonReflectTableProxy PhotonElementProxy::reflectivity_table_sigma() const {
  void* ptr;
  photon_element_struct_get_reflectivity_table_sigma(fortran_ptr_, &ptr);
  return PhotonReflectTableProxy(ptr);
}
void PhotonElementProxy::set_reflectivity_table_sigma(
    const PhotonReflectTableProxy& src) {
  photon_element_struct_set_reflectivity_table_sigma(
      fortran_ptr_, src.get_fortran_ptr());
}
PhotonReflectTableProxy PhotonElementProxy::reflectivity_table_pi() const {
  void* ptr;
  photon_element_struct_get_reflectivity_table_pi(fortran_ptr_, &ptr);
  return PhotonReflectTableProxy(ptr);
}
void PhotonElementProxy::set_reflectivity_table_pi(
    const PhotonReflectTableProxy& src) {
  photon_element_struct_set_reflectivity_table_pi(
      fortran_ptr_, src.get_fortran_ptr());
}
SplineProxyArray1D PhotonElementProxy::init_energy_prob() const {
  return ProxyHelpers::get_type_array_1d<SplineProxyArray1D>(
      fortran_ptr_, photon_element_struct_get_init_energy_prob_info);
}
FArray1D<double> PhotonElementProxy::integrated_init_energy_prob() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, photon_element_struct_get_integrated_init_energy_prob_info);
}
double Wall3dVertexProxy::x() const {
  double value;
  wall3d_vertex_struct_get_x(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_x(double value) {
  wall3d_vertex_struct_set_x(fortran_ptr_, value);
}
double Wall3dVertexProxy::y() const {
  double value;
  wall3d_vertex_struct_get_y(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_y(double value) {
  wall3d_vertex_struct_set_y(fortran_ptr_, value);
}
double Wall3dVertexProxy::radius_x() const {
  double value;
  wall3d_vertex_struct_get_radius_x(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_radius_x(double value) {
  wall3d_vertex_struct_set_radius_x(fortran_ptr_, value);
}
double Wall3dVertexProxy::radius_y() const {
  double value;
  wall3d_vertex_struct_get_radius_y(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_radius_y(double value) {
  wall3d_vertex_struct_set_radius_y(fortran_ptr_, value);
}
double Wall3dVertexProxy::tilt() const {
  double value;
  wall3d_vertex_struct_get_tilt(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_tilt(double value) {
  wall3d_vertex_struct_set_tilt(fortran_ptr_, value);
}
double Wall3dVertexProxy::angle() const {
  double value;
  wall3d_vertex_struct_get_angle(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_angle(double value) {
  wall3d_vertex_struct_set_angle(fortran_ptr_, value);
}
double Wall3dVertexProxy::x0() const {
  double value;
  wall3d_vertex_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_x0(double value) {
  wall3d_vertex_struct_set_x0(fortran_ptr_, value);
}
double Wall3dVertexProxy::y0() const {
  double value;
  wall3d_vertex_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_y0(double value) {
  wall3d_vertex_struct_set_y0(fortran_ptr_, value);
}
int Wall3dVertexProxy::type() const {
  int value;
  wall3d_vertex_struct_get_type(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexProxy::set_type(int value) {
  wall3d_vertex_struct_set_type(fortran_ptr_, value);
}
std::string Wall3dSectionProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wall3d_section_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dSectionProxy::set_name(const std::string& value) {
  wall3d_section_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string Wall3dSectionProxy::material() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wall3d_section_struct_get_material_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dSectionProxy::set_material(const std::string& value) {
  wall3d_section_struct_set_material(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
Wall3dVertexProxyArray1D Wall3dSectionProxy::v() const {
  return ProxyHelpers::get_type_array_1d<Wall3dVertexProxyArray1D>(
      fortran_ptr_, wall3d_section_struct_get_v_info);
}
std::optional<PhotonReflectSurfaceProxy> Wall3dSectionProxy::surface() const {
  void* ptr;
  wall3d_section_struct_get_surface(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return PhotonReflectSurfaceProxy(ptr);
}
void Wall3dSectionProxy::set_surface(const PhotonReflectSurfaceProxy& src) {
  wall3d_section_struct_set_surface(fortran_ptr_, src.get_fortran_ptr());
}
int Wall3dSectionProxy::type() const {
  int value;
  wall3d_section_struct_get_type(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_type(int value) {
  wall3d_section_struct_set_type(fortran_ptr_, value);
}
int Wall3dSectionProxy::n_vertex_input() const {
  int value;
  wall3d_section_struct_get_n_vertex_input(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_n_vertex_input(int value) {
  wall3d_section_struct_set_n_vertex_input(fortran_ptr_, value);
}
int Wall3dSectionProxy::ix_ele() const {
  int value;
  wall3d_section_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_ix_ele(int value) {
  wall3d_section_struct_set_ix_ele(fortran_ptr_, value);
}
int Wall3dSectionProxy::ix_branch() const {
  int value;
  wall3d_section_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_ix_branch(int value) {
  wall3d_section_struct_set_ix_branch(fortran_ptr_, value);
}
int Wall3dSectionProxy::vertices_state() const {
  int value;
  wall3d_section_struct_get_vertices_state(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_vertices_state(int value) {
  wall3d_section_struct_set_vertices_state(fortran_ptr_, value);
}
bool Wall3dSectionProxy::patch_in_region() const {
  bool value;
  wall3d_section_struct_get_patch_in_region(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_patch_in_region(bool value) {
  wall3d_section_struct_set_patch_in_region(fortran_ptr_, value);
}
double Wall3dSectionProxy::thickness() const {
  double value;
  wall3d_section_struct_get_thickness(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_thickness(double value) {
  wall3d_section_struct_set_thickness(fortran_ptr_, value);
}
double Wall3dSectionProxy::s() const {
  double value;
  wall3d_section_struct_get_s(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_s(double value) {
  wall3d_section_struct_set_s(fortran_ptr_, value);
}
FArray1D<double> Wall3dSectionProxy::r0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, wall3d_section_struct_get_r0_info);
}
double Wall3dSectionProxy::dx0_ds() const {
  double value;
  wall3d_section_struct_get_dx0_ds(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_dx0_ds(double value) {
  wall3d_section_struct_set_dx0_ds(fortran_ptr_, value);
}
double Wall3dSectionProxy::dy0_ds() const {
  double value;
  wall3d_section_struct_get_dy0_ds(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_dy0_ds(double value) {
  wall3d_section_struct_set_dy0_ds(fortran_ptr_, value);
}
FArray1D<double> Wall3dSectionProxy::x0_coef() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, wall3d_section_struct_get_x0_coef_info);
}
FArray1D<double> Wall3dSectionProxy::y0_coef() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, wall3d_section_struct_get_y0_coef_info);
}
double Wall3dSectionProxy::dr_ds() const {
  double value;
  wall3d_section_struct_get_dr_ds(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionProxy::set_dr_ds(double value) {
  wall3d_section_struct_set_dr_ds(fortran_ptr_, value);
}
FArray1D<double> Wall3dSectionProxy::p1_coef() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, wall3d_section_struct_get_p1_coef_info);
}
FArray1D<double> Wall3dSectionProxy::p2_coef() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, wall3d_section_struct_get_p2_coef_info);
}
std::string Wall3dProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wall3d_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dProxy::set_name(const std::string& value) {
  wall3d_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int Wall3dProxy::type() const {
  int value;
  wall3d_struct_get_type(fortran_ptr_, &value);
  return value;
}
void Wall3dProxy::set_type(int value) {
  wall3d_struct_set_type(fortran_ptr_, value);
}
int Wall3dProxy::ix_wall3d() const {
  int value;
  wall3d_struct_get_ix_wall3d(fortran_ptr_, &value);
  return value;
}
void Wall3dProxy::set_ix_wall3d(int value) {
  wall3d_struct_set_ix_wall3d(fortran_ptr_, value);
}
int Wall3dProxy::n_link() const {
  int value;
  wall3d_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void Wall3dProxy::set_n_link(int value) {
  wall3d_struct_set_n_link(fortran_ptr_, value);
}
double Wall3dProxy::thickness() const {
  double value;
  wall3d_struct_get_thickness(fortran_ptr_, &value);
  return value;
}
void Wall3dProxy::set_thickness(double value) {
  wall3d_struct_set_thickness(fortran_ptr_, value);
}
std::string Wall3dProxy::clear_material() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wall3d_struct_get_clear_material_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dProxy::set_clear_material(const std::string& value) {
  wall3d_struct_set_clear_material(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string Wall3dProxy::opaque_material() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, wall3d_struct_get_opaque_material_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dProxy::set_opaque_material(const std::string& value) {
  wall3d_struct_set_opaque_material(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool Wall3dProxy::superimpose() const {
  bool value;
  wall3d_struct_get_superimpose(fortran_ptr_, &value);
  return value;
}
void Wall3dProxy::set_superimpose(bool value) {
  wall3d_struct_set_superimpose(fortran_ptr_, value);
}
int Wall3dProxy::ele_anchor_pt() const {
  int value;
  wall3d_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void Wall3dProxy::set_ele_anchor_pt(int value) {
  wall3d_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
Wall3dSectionProxyArray1D Wall3dProxy::section() const {
  return ProxyHelpers::get_type_array_1d<Wall3dSectionProxyArray1D>(
      fortran_ptr_, wall3d_struct_get_section_info);
}
int RamperLordProxy::ix_ele() const {
  int value;
  ramper_lord_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void RamperLordProxy::set_ix_ele(int value) {
  ramper_lord_struct_set_ix_ele(fortran_ptr_, value);
}
int RamperLordProxy::ix_con() const {
  int value;
  ramper_lord_struct_get_ix_con(fortran_ptr_, &value);
  return value;
}
void RamperLordProxy::set_ix_con(int value) {
  ramper_lord_struct_set_ix_con(fortran_ptr_, value);
}
double* RamperLordProxy::attrib_ptr() const {
  double* ptr;
  ramper_lord_struct_get_attrib_ptr(fortran_ptr_, &ptr);
  return ptr;
}
void RamperLordProxy::set_attrib_ptr(double value) {
  ramper_lord_struct_set_attrib_ptr(fortran_ptr_, value);
}
double ControlProxy::value() const {
  double value;
  control_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ControlProxy::set_value(double value) {
  control_struct_set_value(fortran_ptr_, value);
}
FArray1D<double> ControlProxy::y_knot() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, control_struct_get_y_knot_info);
}
ExpressionAtomProxyArray1D ControlProxy::stack() const {
  return ProxyHelpers::get_type_array_1d<ExpressionAtomProxyArray1D>(
      fortran_ptr_, control_struct_get_stack_info);
}
LatEleLocProxy ControlProxy::slave() const {
  void* ptr;
  control_struct_get_slave(fortran_ptr_, &ptr);
  return LatEleLocProxy(ptr);
}
void ControlProxy::set_slave(const LatEleLocProxy& src) {
  control_struct_set_slave(fortran_ptr_, src.get_fortran_ptr());
}
LatEleLocProxy ControlProxy::lord() const {
  void* ptr;
  control_struct_get_lord(fortran_ptr_, &ptr);
  return LatEleLocProxy(ptr);
}
void ControlProxy::set_lord(const LatEleLocProxy& src) {
  control_struct_set_lord(fortran_ptr_, src.get_fortran_ptr());
}
std::string ControlProxy::slave_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, control_struct_get_slave_name_info);
  return std::string(arr.data(), arr.size());
}
void ControlProxy::set_slave_name(const std::string& value) {
  control_struct_set_slave_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string ControlProxy::attribute() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, control_struct_get_attribute_info);
  return std::string(arr.data(), arr.size());
}
void ControlProxy::set_attribute(const std::string& value) {
  control_struct_set_attribute(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int ControlProxy::ix_attrib() const {
  int value;
  control_struct_get_ix_attrib(fortran_ptr_, &value);
  return value;
}
void ControlProxy::set_ix_attrib(int value) {
  control_struct_set_ix_attrib(fortran_ptr_, value);
}
std::string ControlVar1Proxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, control_var1_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void ControlVar1Proxy::set_name(const std::string& value) {
  control_var1_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double ControlVar1Proxy::value() const {
  double value;
  control_var1_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ControlVar1Proxy::set_value(double value) {
  control_var1_struct_set_value(fortran_ptr_, value);
}
double ControlVar1Proxy::old_value() const {
  double value;
  control_var1_struct_get_old_value(fortran_ptr_, &value);
  return value;
}
void ControlVar1Proxy::set_old_value(double value) {
  control_var1_struct_set_old_value(fortran_ptr_, value);
}
FArray1D<double> ControlRamp1Proxy::y_knot() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, control_ramp1_struct_get_y_knot_info);
}
ExpressionAtomProxyArray1D ControlRamp1Proxy::stack() const {
  return ProxyHelpers::get_type_array_1d<ExpressionAtomProxyArray1D>(
      fortran_ptr_, control_ramp1_struct_get_stack_info);
}
std::string ControlRamp1Proxy::attribute() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, control_ramp1_struct_get_attribute_info);
  return std::string(arr.data(), arr.size());
}
void ControlRamp1Proxy::set_attribute(const std::string& value) {
  control_ramp1_struct_set_attribute(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string ControlRamp1Proxy::slave_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, control_ramp1_struct_get_slave_name_info);
  return std::string(arr.data(), arr.size());
}
void ControlRamp1Proxy::set_slave_name(const std::string& value) {
  control_ramp1_struct_set_slave_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool ControlRamp1Proxy::is_controller() const {
  bool value;
  control_ramp1_struct_get_is_controller(fortran_ptr_, &value);
  return value;
}
void ControlRamp1Proxy::set_is_controller(bool value) {
  control_ramp1_struct_set_is_controller(fortran_ptr_, value);
}
ControlVar1ProxyArray1D ControllerProxy::var() const {
  return ProxyHelpers::get_type_array_1d<ControlVar1ProxyArray1D>(
      fortran_ptr_, controller_struct_get_var_info);
}
ControlRamp1ProxyArray1D ControllerProxy::ramp() const {
  return ProxyHelpers::get_type_array_1d<ControlRamp1ProxyArray1D>(
      fortran_ptr_, controller_struct_get_ramp_info);
}
RamperLordProxyArray1D ControllerProxy::ramper_lord() const {
  return ProxyHelpers::get_type_array_1d<RamperLordProxyArray1D>(
      fortran_ptr_, controller_struct_get_ramper_lord_info);
}
FArray1D<double> ControllerProxy::x_knot() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, controller_struct_get_x_knot_info);
}
int EllipseBeamInitProxy::part_per_ellipse() const {
  int value;
  ellipse_beam_init_struct_get_part_per_ellipse(fortran_ptr_, &value);
  return value;
}
void EllipseBeamInitProxy::set_part_per_ellipse(int value) {
  ellipse_beam_init_struct_set_part_per_ellipse(fortran_ptr_, value);
}
int EllipseBeamInitProxy::n_ellipse() const {
  int value;
  ellipse_beam_init_struct_get_n_ellipse(fortran_ptr_, &value);
  return value;
}
void EllipseBeamInitProxy::set_n_ellipse(int value) {
  ellipse_beam_init_struct_set_n_ellipse(fortran_ptr_, value);
}
double EllipseBeamInitProxy::sigma_cutoff() const {
  double value;
  ellipse_beam_init_struct_get_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void EllipseBeamInitProxy::set_sigma_cutoff(double value) {
  ellipse_beam_init_struct_set_sigma_cutoff(fortran_ptr_, value);
}
FArray1D<int> KvBeamInitProxy::part_per_phi() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, kv_beam_init_struct_get_part_per_phi_info);
}
int KvBeamInitProxy::n_I2() const {
  int value;
  kv_beam_init_struct_get_n_I2(fortran_ptr_, &value);
  return value;
}
void KvBeamInitProxy::set_n_I2(int value) {
  kv_beam_init_struct_set_n_I2(fortran_ptr_, value);
}
double KvBeamInitProxy::A() const {
  double value;
  kv_beam_init_struct_get_A(fortran_ptr_, &value);
  return value;
}
void KvBeamInitProxy::set_A(double value) {
  kv_beam_init_struct_set_A(fortran_ptr_, value);
}
int GridBeamInitProxy::n_x() const {
  int value;
  grid_beam_init_struct_get_n_x(fortran_ptr_, &value);
  return value;
}
void GridBeamInitProxy::set_n_x(int value) {
  grid_beam_init_struct_set_n_x(fortran_ptr_, value);
}
int GridBeamInitProxy::n_px() const {
  int value;
  grid_beam_init_struct_get_n_px(fortran_ptr_, &value);
  return value;
}
void GridBeamInitProxy::set_n_px(int value) {
  grid_beam_init_struct_set_n_px(fortran_ptr_, value);
}
double GridBeamInitProxy::x_min() const {
  double value;
  grid_beam_init_struct_get_x_min(fortran_ptr_, &value);
  return value;
}
void GridBeamInitProxy::set_x_min(double value) {
  grid_beam_init_struct_set_x_min(fortran_ptr_, value);
}
double GridBeamInitProxy::x_max() const {
  double value;
  grid_beam_init_struct_get_x_max(fortran_ptr_, &value);
  return value;
}
void GridBeamInitProxy::set_x_max(double value) {
  grid_beam_init_struct_set_x_max(fortran_ptr_, value);
}
double GridBeamInitProxy::px_min() const {
  double value;
  grid_beam_init_struct_get_px_min(fortran_ptr_, &value);
  return value;
}
void GridBeamInitProxy::set_px_min(double value) {
  grid_beam_init_struct_set_px_min(fortran_ptr_, value);
}
double GridBeamInitProxy::px_max() const {
  double value;
  grid_beam_init_struct_get_px_max(fortran_ptr_, &value);
  return value;
}
void GridBeamInitProxy::set_px_max(double value) {
  grid_beam_init_struct_set_px_max(fortran_ptr_, value);
}
std::string BeamInitProxy::position_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, beam_init_struct_get_position_file_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitProxy::set_position_file(const std::string& value) {
  beam_init_struct_set_position_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D BeamInitProxy::distribution_type() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, beam_init_struct_get_distribution_type_info);
}
FArray1D<double> BeamInitProxy::spin() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, beam_init_struct_get_spin_info);
}
EllipseBeamInitProxyArray1D BeamInitProxy::ellipse() const {
  return ProxyHelpers::get_type_array_1d<EllipseBeamInitProxyArray1D>(
      fortran_ptr_, beam_init_struct_get_ellipse_info);
}
KvBeamInitProxy BeamInitProxy::KV() const {
  void* ptr;
  beam_init_struct_get_KV(fortran_ptr_, &ptr);
  return KvBeamInitProxy(ptr);
}
void BeamInitProxy::set_KV(const KvBeamInitProxy& src) {
  beam_init_struct_set_KV(fortran_ptr_, src.get_fortran_ptr());
}
GridBeamInitProxyArray1D BeamInitProxy::grid() const {
  return ProxyHelpers::get_type_array_1d<GridBeamInitProxyArray1D>(
      fortran_ptr_, beam_init_struct_get_grid_info);
}
FArray1D<double> BeamInitProxy::center_jitter() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, beam_init_struct_get_center_jitter_info);
}
FArray1D<double> BeamInitProxy::emit_jitter() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, beam_init_struct_get_emit_jitter_info);
}
double BeamInitProxy::sig_z_jitter() const {
  double value;
  beam_init_struct_get_sig_z_jitter(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_sig_z_jitter(double value) {
  beam_init_struct_set_sig_z_jitter(fortran_ptr_, value);
}
double BeamInitProxy::sig_pz_jitter() const {
  double value;
  beam_init_struct_get_sig_pz_jitter(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_sig_pz_jitter(double value) {
  beam_init_struct_set_sig_pz_jitter(fortran_ptr_, value);
}
int BeamInitProxy::n_particle() const {
  int value;
  beam_init_struct_get_n_particle(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_n_particle(int value) {
  beam_init_struct_set_n_particle(fortran_ptr_, value);
}
bool BeamInitProxy::renorm_center() const {
  bool value;
  beam_init_struct_get_renorm_center(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_renorm_center(bool value) {
  beam_init_struct_set_renorm_center(fortran_ptr_, value);
}
bool BeamInitProxy::renorm_sigma() const {
  bool value;
  beam_init_struct_get_renorm_sigma(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_renorm_sigma(bool value) {
  beam_init_struct_set_renorm_sigma(fortran_ptr_, value);
}
std::string BeamInitProxy::random_engine() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, beam_init_struct_get_random_engine_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitProxy::set_random_engine(const std::string& value) {
  beam_init_struct_set_random_engine(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string BeamInitProxy::random_gauss_converter() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, beam_init_struct_get_random_gauss_converter_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitProxy::set_random_gauss_converter(const std::string& value) {
  beam_init_struct_set_random_gauss_converter(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double BeamInitProxy::random_sigma_cutoff() const {
  double value;
  beam_init_struct_get_random_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_random_sigma_cutoff(double value) {
  beam_init_struct_set_random_sigma_cutoff(fortran_ptr_, value);
}
double BeamInitProxy::a_norm_emit() const {
  double value;
  beam_init_struct_get_a_norm_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_a_norm_emit(double value) {
  beam_init_struct_set_a_norm_emit(fortran_ptr_, value);
}
double BeamInitProxy::b_norm_emit() const {
  double value;
  beam_init_struct_get_b_norm_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_b_norm_emit(double value) {
  beam_init_struct_set_b_norm_emit(fortran_ptr_, value);
}
double BeamInitProxy::a_emit() const {
  double value;
  beam_init_struct_get_a_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_a_emit(double value) {
  beam_init_struct_set_a_emit(fortran_ptr_, value);
}
double BeamInitProxy::b_emit() const {
  double value;
  beam_init_struct_get_b_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_b_emit(double value) {
  beam_init_struct_set_b_emit(fortran_ptr_, value);
}
double BeamInitProxy::dPz_dz() const {
  double value;
  beam_init_struct_get_dPz_dz(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_dPz_dz(double value) {
  beam_init_struct_set_dPz_dz(fortran_ptr_, value);
}
FArray1D<double> BeamInitProxy::center() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, beam_init_struct_get_center_info);
}
double BeamInitProxy::t_offset() const {
  double value;
  beam_init_struct_get_t_offset(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_t_offset(double value) {
  beam_init_struct_set_t_offset(fortran_ptr_, value);
}
double BeamInitProxy::dt_bunch() const {
  double value;
  beam_init_struct_get_dt_bunch(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_dt_bunch(double value) {
  beam_init_struct_set_dt_bunch(fortran_ptr_, value);
}
double BeamInitProxy::sig_z() const {
  double value;
  beam_init_struct_get_sig_z(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_sig_z(double value) {
  beam_init_struct_set_sig_z(fortran_ptr_, value);
}
double BeamInitProxy::sig_pz() const {
  double value;
  beam_init_struct_get_sig_pz(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_sig_pz(double value) {
  beam_init_struct_set_sig_pz(fortran_ptr_, value);
}
double BeamInitProxy::bunch_charge() const {
  double value;
  beam_init_struct_get_bunch_charge(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_bunch_charge(double value) {
  beam_init_struct_set_bunch_charge(fortran_ptr_, value);
}
int BeamInitProxy::n_bunch() const {
  int value;
  beam_init_struct_get_n_bunch(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_n_bunch(int value) {
  beam_init_struct_set_n_bunch(fortran_ptr_, value);
}
int BeamInitProxy::ix_turn() const {
  int value;
  beam_init_struct_get_ix_turn(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_ix_turn(int value) {
  beam_init_struct_set_ix_turn(fortran_ptr_, value);
}
std::string BeamInitProxy::species() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, beam_init_struct_get_species_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitProxy::set_species(const std::string& value) {
  beam_init_struct_set_species(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool BeamInitProxy::full_6D_coupling_calc() const {
  bool value;
  beam_init_struct_get_full_6D_coupling_calc(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_full_6D_coupling_calc(bool value) {
  beam_init_struct_set_full_6D_coupling_calc(fortran_ptr_, value);
}
bool BeamInitProxy::use_particle_start() const {
  bool value;
  beam_init_struct_get_use_particle_start(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_use_particle_start(bool value) {
  beam_init_struct_set_use_particle_start(fortran_ptr_, value);
}
bool BeamInitProxy::use_t_coords() const {
  bool value;
  beam_init_struct_get_use_t_coords(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_use_t_coords(bool value) {
  beam_init_struct_set_use_t_coords(fortran_ptr_, value);
}
bool BeamInitProxy::use_z_as_t() const {
  bool value;
  beam_init_struct_get_use_z_as_t(fortran_ptr_, &value);
  return value;
}
void BeamInitProxy::set_use_z_as_t(bool value) {
  beam_init_struct_set_use_z_as_t(fortran_ptr_, value);
}
std::string BeamInitProxy::file_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, beam_init_struct_get_file_name_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitProxy::set_file_name(const std::string& value) {
  beam_init_struct_set_file_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double LatParamProxy::n_part() const {
  double value;
  lat_param_struct_get_n_part(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_n_part(double value) {
  lat_param_struct_set_n_part(fortran_ptr_, value);
}
double LatParamProxy::total_length() const {
  double value;
  lat_param_struct_get_total_length(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_total_length(double value) {
  lat_param_struct_set_total_length(fortran_ptr_, value);
}
double LatParamProxy::unstable_factor() const {
  double value;
  lat_param_struct_get_unstable_factor(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_unstable_factor(double value) {
  lat_param_struct_set_unstable_factor(fortran_ptr_, value);
}
FArray2D<double> LatParamProxy::t1_with_RF() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, lat_param_struct_get_t1_with_RF_info);
}
FArray2D<double> LatParamProxy::t1_no_RF() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, lat_param_struct_get_t1_no_RF_info);
}
double LatParamProxy::spin_tune() const {
  double value;
  lat_param_struct_get_spin_tune(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_spin_tune(double value) {
  lat_param_struct_set_spin_tune(fortran_ptr_, value);
}
int LatParamProxy::particle() const {
  int value;
  lat_param_struct_get_particle(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_particle(int value) {
  lat_param_struct_set_particle(fortran_ptr_, value);
}
int LatParamProxy::default_tracking_species() const {
  int value;
  lat_param_struct_get_default_tracking_species(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_default_tracking_species(int value) {
  lat_param_struct_set_default_tracking_species(fortran_ptr_, value);
}
int LatParamProxy::geometry() const {
  int value;
  lat_param_struct_get_geometry(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_geometry(int value) {
  lat_param_struct_set_geometry(fortran_ptr_, value);
}
int LatParamProxy::ixx() const {
  int value;
  lat_param_struct_get_ixx(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_ixx(int value) {
  lat_param_struct_set_ixx(fortran_ptr_, value);
}
bool LatParamProxy::stable() const {
  bool value;
  lat_param_struct_get_stable(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_stable(bool value) {
  lat_param_struct_set_stable(fortran_ptr_, value);
}
bool LatParamProxy::live_branch() const {
  bool value;
  lat_param_struct_get_live_branch(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_live_branch(bool value) {
  lat_param_struct_set_live_branch(fortran_ptr_, value);
}
double LatParamProxy::g1_integral() const {
  double value;
  lat_param_struct_get_g1_integral(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_g1_integral(double value) {
  lat_param_struct_set_g1_integral(fortran_ptr_, value);
}
double LatParamProxy::g2_integral() const {
  double value;
  lat_param_struct_get_g2_integral(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_g2_integral(double value) {
  lat_param_struct_set_g2_integral(fortran_ptr_, value);
}
double LatParamProxy::g3_integral() const {
  double value;
  lat_param_struct_get_g3_integral(fortran_ptr_, &value);
  return value;
}
void LatParamProxy::set_g3_integral(double value) {
  lat_param_struct_set_g3_integral(fortran_ptr_, value);
}
BookkeepingStateProxy LatParamProxy::bookkeeping_state() const {
  void* ptr;
  lat_param_struct_get_bookkeeping_state(fortran_ptr_, &ptr);
  return BookkeepingStateProxy(ptr);
}
void LatParamProxy::set_bookkeeping_state(const BookkeepingStateProxy& src) {
  lat_param_struct_set_bookkeeping_state(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitProxy LatParamProxy::beam_init() const {
  void* ptr;
  lat_param_struct_get_beam_init(fortran_ptr_, &ptr);
  return BeamInitProxy(ptr);
}
void LatParamProxy::set_beam_init(const BeamInitProxy& src) {
  lat_param_struct_set_beam_init(fortran_ptr_, src.get_fortran_ptr());
}
bool ModeInfoProxy::stable() const {
  bool value;
  mode_info_struct_get_stable(fortran_ptr_, &value);
  return value;
}
void ModeInfoProxy::set_stable(bool value) {
  mode_info_struct_set_stable(fortran_ptr_, value);
}
double ModeInfoProxy::tune() const {
  double value;
  mode_info_struct_get_tune(fortran_ptr_, &value);
  return value;
}
void ModeInfoProxy::set_tune(double value) {
  mode_info_struct_set_tune(fortran_ptr_, value);
}
double ModeInfoProxy::emit() const {
  double value;
  mode_info_struct_get_emit(fortran_ptr_, &value);
  return value;
}
void ModeInfoProxy::set_emit(double value) {
  mode_info_struct_set_emit(fortran_ptr_, value);
}
double ModeInfoProxy::chrom() const {
  double value;
  mode_info_struct_get_chrom(fortran_ptr_, &value);
  return value;
}
void ModeInfoProxy::set_chrom(double value) {
  mode_info_struct_set_chrom(fortran_ptr_, value);
}
double ModeInfoProxy::sigma() const {
  double value;
  mode_info_struct_get_sigma(fortran_ptr_, &value);
  return value;
}
void ModeInfoProxy::set_sigma(double value) {
  mode_info_struct_set_sigma(fortran_ptr_, value);
}
double ModeInfoProxy::sigmap() const {
  double value;
  mode_info_struct_get_sigmap(fortran_ptr_, &value);
  return value;
}
void ModeInfoProxy::set_sigmap(double value) {
  mode_info_struct_set_sigmap(fortran_ptr_, value);
}
int PreTrackerProxy::who() const {
  int value;
  pre_tracker_struct_get_who(fortran_ptr_, &value);
  return value;
}
void PreTrackerProxy::set_who(int value) {
  pre_tracker_struct_set_who(fortran_ptr_, value);
}
int PreTrackerProxy::ix_ele_start() const {
  int value;
  pre_tracker_struct_get_ix_ele_start(fortran_ptr_, &value);
  return value;
}
void PreTrackerProxy::set_ix_ele_start(int value) {
  pre_tracker_struct_set_ix_ele_start(fortran_ptr_, value);
}
int PreTrackerProxy::ix_ele_end() const {
  int value;
  pre_tracker_struct_get_ix_ele_end(fortran_ptr_, &value);
  return value;
}
void PreTrackerProxy::set_ix_ele_end(int value) {
  pre_tracker_struct_set_ix_ele_end(fortran_ptr_, value);
}
std::string PreTrackerProxy::input_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, pre_tracker_struct_get_input_file_info);
  return std::string(arr.data(), arr.size());
}
void PreTrackerProxy::set_input_file(const std::string& value) {
  pre_tracker_struct_set_input_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double AnormalModeProxy::emittance() const {
  double value;
  anormal_mode_struct_get_emittance(fortran_ptr_, &value);
  return value;
}
void AnormalModeProxy::set_emittance(double value) {
  anormal_mode_struct_set_emittance(fortran_ptr_, value);
}
double AnormalModeProxy::emittance_no_vert() const {
  double value;
  anormal_mode_struct_get_emittance_no_vert(fortran_ptr_, &value);
  return value;
}
void AnormalModeProxy::set_emittance_no_vert(double value) {
  anormal_mode_struct_set_emittance_no_vert(fortran_ptr_, value);
}
FArray1D<double> AnormalModeProxy::synch_int() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, anormal_mode_struct_get_synch_int_info);
}
double AnormalModeProxy::j_damp() const {
  double value;
  anormal_mode_struct_get_j_damp(fortran_ptr_, &value);
  return value;
}
void AnormalModeProxy::set_j_damp(double value) {
  anormal_mode_struct_set_j_damp(fortran_ptr_, value);
}
double AnormalModeProxy::alpha_damp() const {
  double value;
  anormal_mode_struct_get_alpha_damp(fortran_ptr_, &value);
  return value;
}
void AnormalModeProxy::set_alpha_damp(double value) {
  anormal_mode_struct_set_alpha_damp(fortran_ptr_, value);
}
double AnormalModeProxy::chrom() const {
  double value;
  anormal_mode_struct_get_chrom(fortran_ptr_, &value);
  return value;
}
void AnormalModeProxy::set_chrom(double value) {
  anormal_mode_struct_set_chrom(fortran_ptr_, value);
}
double AnormalModeProxy::tune() const {
  double value;
  anormal_mode_struct_get_tune(fortran_ptr_, &value);
  return value;
}
void AnormalModeProxy::set_tune(double value) {
  anormal_mode_struct_set_tune(fortran_ptr_, value);
}
double LinacNormalModeProxy::i2_E4() const {
  double value;
  linac_normal_mode_struct_get_i2_E4(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_i2_E4(double value) {
  linac_normal_mode_struct_set_i2_E4(fortran_ptr_, value);
}
double LinacNormalModeProxy::i3_E7() const {
  double value;
  linac_normal_mode_struct_get_i3_E7(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_i3_E7(double value) {
  linac_normal_mode_struct_set_i3_E7(fortran_ptr_, value);
}
double LinacNormalModeProxy::i5a_E6() const {
  double value;
  linac_normal_mode_struct_get_i5a_E6(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_i5a_E6(double value) {
  linac_normal_mode_struct_set_i5a_E6(fortran_ptr_, value);
}
double LinacNormalModeProxy::i5b_E6() const {
  double value;
  linac_normal_mode_struct_get_i5b_E6(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_i5b_E6(double value) {
  linac_normal_mode_struct_set_i5b_E6(fortran_ptr_, value);
}
double LinacNormalModeProxy::sig_E1() const {
  double value;
  linac_normal_mode_struct_get_sig_E1(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_sig_E1(double value) {
  linac_normal_mode_struct_set_sig_E1(fortran_ptr_, value);
}
double LinacNormalModeProxy::a_emittance_end() const {
  double value;
  linac_normal_mode_struct_get_a_emittance_end(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_a_emittance_end(double value) {
  linac_normal_mode_struct_set_a_emittance_end(fortran_ptr_, value);
}
double LinacNormalModeProxy::b_emittance_end() const {
  double value;
  linac_normal_mode_struct_get_b_emittance_end(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeProxy::set_b_emittance_end(double value) {
  linac_normal_mode_struct_set_b_emittance_end(fortran_ptr_, value);
}
FArray1D<double> NormalModesProxy::synch_int() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, normal_modes_struct_get_synch_int_info);
}
double NormalModesProxy::sigE_E() const {
  double value;
  normal_modes_struct_get_sigE_E(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_sigE_E(double value) {
  normal_modes_struct_set_sigE_E(fortran_ptr_, value);
}
double NormalModesProxy::sig_z() const {
  double value;
  normal_modes_struct_get_sig_z(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_sig_z(double value) {
  normal_modes_struct_set_sig_z(fortran_ptr_, value);
}
double NormalModesProxy::e_loss() const {
  double value;
  normal_modes_struct_get_e_loss(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_e_loss(double value) {
  normal_modes_struct_set_e_loss(fortran_ptr_, value);
}
double NormalModesProxy::rf_voltage() const {
  double value;
  normal_modes_struct_get_rf_voltage(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_rf_voltage(double value) {
  normal_modes_struct_set_rf_voltage(fortran_ptr_, value);
}
double NormalModesProxy::pz_aperture() const {
  double value;
  normal_modes_struct_get_pz_aperture(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_pz_aperture(double value) {
  normal_modes_struct_set_pz_aperture(fortran_ptr_, value);
}
double NormalModesProxy::pz_average() const {
  double value;
  normal_modes_struct_get_pz_average(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_pz_average(double value) {
  normal_modes_struct_set_pz_average(fortran_ptr_, value);
}
double NormalModesProxy::momentum_compaction() const {
  double value;
  normal_modes_struct_get_momentum_compaction(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_momentum_compaction(double value) {
  normal_modes_struct_set_momentum_compaction(fortran_ptr_, value);
}
double NormalModesProxy::dpz_damp() const {
  double value;
  normal_modes_struct_get_dpz_damp(fortran_ptr_, &value);
  return value;
}
void NormalModesProxy::set_dpz_damp(double value) {
  normal_modes_struct_set_dpz_damp(fortran_ptr_, value);
}
AnormalModeProxy NormalModesProxy::a() const {
  void* ptr;
  normal_modes_struct_get_a(fortran_ptr_, &ptr);
  return AnormalModeProxy(ptr);
}
void NormalModesProxy::set_a(const AnormalModeProxy& src) {
  normal_modes_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
AnormalModeProxy NormalModesProxy::b() const {
  void* ptr;
  normal_modes_struct_get_b(fortran_ptr_, &ptr);
  return AnormalModeProxy(ptr);
}
void NormalModesProxy::set_b(const AnormalModeProxy& src) {
  normal_modes_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
AnormalModeProxy NormalModesProxy::z() const {
  void* ptr;
  normal_modes_struct_get_z(fortran_ptr_, &ptr);
  return AnormalModeProxy(ptr);
}
void NormalModesProxy::set_z(const AnormalModeProxy& src) {
  normal_modes_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
LinacNormalModeProxy NormalModesProxy::lin() const {
  void* ptr;
  normal_modes_struct_get_lin(fortran_ptr_, &ptr);
  return LinacNormalModeProxy(ptr);
}
void NormalModesProxy::set_lin(const LinacNormalModeProxy& src) {
  normal_modes_struct_set_lin(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> EmFieldProxy::E() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, em_field_struct_get_E_info);
}
FArray1D<double> EmFieldProxy::B() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, em_field_struct_get_B_info);
}
FArray2D<double> EmFieldProxy::dE() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, em_field_struct_get_dE_info);
}
FArray2D<double> EmFieldProxy::dB() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, em_field_struct_get_dB_info);
}
double EmFieldProxy::phi() const {
  double value;
  em_field_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void EmFieldProxy::set_phi(double value) {
  em_field_struct_set_phi(fortran_ptr_, value);
}
double EmFieldProxy::phi_B() const {
  double value;
  em_field_struct_get_phi_B(fortran_ptr_, &value);
  return value;
}
void EmFieldProxy::set_phi_B(double value) {
  em_field_struct_set_phi_B(fortran_ptr_, value);
}
FArray1D<double> EmFieldProxy::A() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, em_field_struct_get_A_info);
}
int StrongBeamProxy::ix_slice() const {
  int value;
  strong_beam_struct_get_ix_slice(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_ix_slice(int value) {
  strong_beam_struct_set_ix_slice(fortran_ptr_, value);
}
double StrongBeamProxy::x_center() const {
  double value;
  strong_beam_struct_get_x_center(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_x_center(double value) {
  strong_beam_struct_set_x_center(fortran_ptr_, value);
}
double StrongBeamProxy::y_center() const {
  double value;
  strong_beam_struct_get_y_center(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_y_center(double value) {
  strong_beam_struct_set_y_center(fortran_ptr_, value);
}
double StrongBeamProxy::x_sigma() const {
  double value;
  strong_beam_struct_get_x_sigma(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_x_sigma(double value) {
  strong_beam_struct_set_x_sigma(fortran_ptr_, value);
}
double StrongBeamProxy::y_sigma() const {
  double value;
  strong_beam_struct_get_y_sigma(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_y_sigma(double value) {
  strong_beam_struct_set_y_sigma(fortran_ptr_, value);
}
double StrongBeamProxy::dx() const {
  double value;
  strong_beam_struct_get_dx(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_dx(double value) {
  strong_beam_struct_set_dx(fortran_ptr_, value);
}
double StrongBeamProxy::dy() const {
  double value;
  strong_beam_struct_get_dy(fortran_ptr_, &value);
  return value;
}
void StrongBeamProxy::set_dy(double value) {
  strong_beam_struct_set_dy(fortran_ptr_, value);
}
double TrackPointProxy::s_lab() const {
  double value;
  track_point_struct_get_s_lab(fortran_ptr_, &value);
  return value;
}
void TrackPointProxy::set_s_lab(double value) {
  track_point_struct_set_s_lab(fortran_ptr_, value);
}
double TrackPointProxy::s_body() const {
  double value;
  track_point_struct_get_s_body(fortran_ptr_, &value);
  return value;
}
void TrackPointProxy::set_s_body(double value) {
  track_point_struct_set_s_body(fortran_ptr_, value);
}
CoordProxy TrackPointProxy::orb() const {
  void* ptr;
  track_point_struct_get_orb(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void TrackPointProxy::set_orb(const CoordProxy& src) {
  track_point_struct_set_orb(fortran_ptr_, src.get_fortran_ptr());
}
EmFieldProxy TrackPointProxy::field() const {
  void* ptr;
  track_point_struct_get_field(fortran_ptr_, &ptr);
  return EmFieldProxy(ptr);
}
void TrackPointProxy::set_field(const EmFieldProxy& src) {
  track_point_struct_set_field(fortran_ptr_, src.get_fortran_ptr());
}
StrongBeamProxy TrackPointProxy::strong_beam() const {
  void* ptr;
  track_point_struct_get_strong_beam(fortran_ptr_, &ptr);
  return StrongBeamProxy(ptr);
}
void TrackPointProxy::set_strong_beam(const StrongBeamProxy& src) {
  track_point_struct_set_strong_beam(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TrackPointProxy::vec0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, track_point_struct_get_vec0_info);
}
FArray2D<double> TrackPointProxy::mat6() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, track_point_struct_get_mat6_info);
}
TrackPointProxyArray1D TrackProxy::pt() const {
  return ProxyHelpers::get_type_array_1d<TrackPointProxyArray1D>(
      fortran_ptr_, track_struct_get_pt_info);
}
double TrackProxy::ds_save() const {
  double value;
  track_struct_get_ds_save(fortran_ptr_, &value);
  return value;
}
void TrackProxy::set_ds_save(double value) {
  track_struct_set_ds_save(fortran_ptr_, value);
}
int TrackProxy::n_pt() const {
  int value;
  track_struct_get_n_pt(fortran_ptr_, &value);
  return value;
}
void TrackProxy::set_n_pt(int value) {
  track_struct_set_n_pt(fortran_ptr_, value);
}
int TrackProxy::n_bad() const {
  int value;
  track_struct_get_n_bad(fortran_ptr_, &value);
  return value;
}
void TrackProxy::set_n_bad(int value) {
  track_struct_set_n_bad(fortran_ptr_, value);
}
int TrackProxy::n_ok() const {
  int value;
  track_struct_get_n_ok(fortran_ptr_, &value);
  return value;
}
void TrackProxy::set_n_ok(int value) {
  track_struct_set_n_ok(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::ds_track_step() const {
  double value;
  space_charge_common_struct_get_ds_track_step(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_ds_track_step(double value) {
  space_charge_common_struct_set_ds_track_step(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::dt_track_step() const {
  double value;
  space_charge_common_struct_get_dt_track_step(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_dt_track_step(double value) {
  space_charge_common_struct_set_dt_track_step(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::cathode_strength_cutoff() const {
  double value;
  space_charge_common_struct_get_cathode_strength_cutoff(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_cathode_strength_cutoff(double value) {
  space_charge_common_struct_set_cathode_strength_cutoff(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::rel_tol_tracking() const {
  double value;
  space_charge_common_struct_get_rel_tol_tracking(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_rel_tol_tracking(double value) {
  space_charge_common_struct_set_rel_tol_tracking(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::abs_tol_tracking() const {
  double value;
  space_charge_common_struct_get_abs_tol_tracking(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_abs_tol_tracking(double value) {
  space_charge_common_struct_set_abs_tol_tracking(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::beam_chamber_height() const {
  double value;
  space_charge_common_struct_get_beam_chamber_height(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_beam_chamber_height(double value) {
  space_charge_common_struct_set_beam_chamber_height(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::lsc_sigma_cutoff() const {
  double value;
  space_charge_common_struct_get_lsc_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_lsc_sigma_cutoff(double value) {
  space_charge_common_struct_set_lsc_sigma_cutoff(fortran_ptr_, value);
}
double SpaceChargeCommonProxy::particle_sigma_cutoff() const {
  double value;
  space_charge_common_struct_get_particle_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_particle_sigma_cutoff(double value) {
  space_charge_common_struct_set_particle_sigma_cutoff(fortran_ptr_, value);
}
FArray1D<int> SpaceChargeCommonProxy::space_charge_mesh_size() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, space_charge_common_struct_get_space_charge_mesh_size_info);
}
FArray1D<int> SpaceChargeCommonProxy::csr3d_mesh_size() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, space_charge_common_struct_get_csr3d_mesh_size_info);
}
int SpaceChargeCommonProxy::n_bin() const {
  int value;
  space_charge_common_struct_get_n_bin(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_n_bin(int value) {
  space_charge_common_struct_set_n_bin(fortran_ptr_, value);
}
int SpaceChargeCommonProxy::particle_bin_span() const {
  int value;
  space_charge_common_struct_get_particle_bin_span(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_particle_bin_span(int value) {
  space_charge_common_struct_set_particle_bin_span(fortran_ptr_, value);
}
int SpaceChargeCommonProxy::n_shield_images() const {
  int value;
  space_charge_common_struct_get_n_shield_images(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_n_shield_images(int value) {
  space_charge_common_struct_set_n_shield_images(fortran_ptr_, value);
}
int SpaceChargeCommonProxy::sc_min_in_bin() const {
  int value;
  space_charge_common_struct_get_sc_min_in_bin(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_sc_min_in_bin(int value) {
  space_charge_common_struct_set_sc_min_in_bin(fortran_ptr_, value);
}
bool SpaceChargeCommonProxy::lsc_kick_transverse_dependence() const {
  bool value;
  space_charge_common_struct_get_lsc_kick_transverse_dependence(
      fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_lsc_kick_transverse_dependence(bool value) {
  space_charge_common_struct_set_lsc_kick_transverse_dependence(
      fortran_ptr_, value);
}
bool SpaceChargeCommonProxy::debug() const {
  bool value;
  space_charge_common_struct_get_debug(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonProxy::set_debug(bool value) {
  space_charge_common_struct_set_debug(fortran_ptr_, value);
}
std::string SpaceChargeCommonProxy::diagnostic_output_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, space_charge_common_struct_get_diagnostic_output_file_info);
  return std::string(arr.data(), arr.size());
}
void SpaceChargeCommonProxy::set_diagnostic_output_file(
    const std::string& value) {
  space_charge_common_struct_set_diagnostic_output_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double BmadCommonProxy::max_aperture_limit() const {
  double value;
  bmad_common_struct_get_max_aperture_limit(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_max_aperture_limit(double value) {
  bmad_common_struct_set_max_aperture_limit(fortran_ptr_, value);
}
FArray1D<double> BmadCommonProxy::d_orb() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, bmad_common_struct_get_d_orb_info);
}
double BmadCommonProxy::default_ds_step() const {
  double value;
  bmad_common_struct_get_default_ds_step(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_default_ds_step(double value) {
  bmad_common_struct_set_default_ds_step(fortran_ptr_, value);
}
double BmadCommonProxy::significant_length() const {
  double value;
  bmad_common_struct_get_significant_length(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_significant_length(double value) {
  bmad_common_struct_set_significant_length(fortran_ptr_, value);
}
double BmadCommonProxy::rel_tol_tracking() const {
  double value;
  bmad_common_struct_get_rel_tol_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_rel_tol_tracking(double value) {
  bmad_common_struct_set_rel_tol_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::abs_tol_tracking() const {
  double value;
  bmad_common_struct_get_abs_tol_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_abs_tol_tracking(double value) {
  bmad_common_struct_set_abs_tol_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::rel_tol_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_rel_tol_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_rel_tol_adaptive_tracking(double value) {
  bmad_common_struct_set_rel_tol_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::abs_tol_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_abs_tol_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_abs_tol_adaptive_tracking(double value) {
  bmad_common_struct_set_abs_tol_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::init_ds_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_init_ds_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_init_ds_adaptive_tracking(double value) {
  bmad_common_struct_set_init_ds_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::min_ds_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_min_ds_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_min_ds_adaptive_tracking(double value) {
  bmad_common_struct_set_min_ds_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::fatal_ds_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_fatal_ds_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_fatal_ds_adaptive_tracking(double value) {
  bmad_common_struct_set_fatal_ds_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonProxy::autoscale_amp_abs_tol() const {
  double value;
  bmad_common_struct_get_autoscale_amp_abs_tol(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_autoscale_amp_abs_tol(double value) {
  bmad_common_struct_set_autoscale_amp_abs_tol(fortran_ptr_, value);
}
double BmadCommonProxy::autoscale_amp_rel_tol() const {
  double value;
  bmad_common_struct_get_autoscale_amp_rel_tol(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_autoscale_amp_rel_tol(double value) {
  bmad_common_struct_set_autoscale_amp_rel_tol(fortran_ptr_, value);
}
double BmadCommonProxy::autoscale_phase_tol() const {
  double value;
  bmad_common_struct_get_autoscale_phase_tol(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_autoscale_phase_tol(double value) {
  bmad_common_struct_set_autoscale_phase_tol(fortran_ptr_, value);
}
double BmadCommonProxy::electric_dipole_moment() const {
  double value;
  bmad_common_struct_get_electric_dipole_moment(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_electric_dipole_moment(double value) {
  bmad_common_struct_set_electric_dipole_moment(fortran_ptr_, value);
}
double BmadCommonProxy::synch_rad_scale() const {
  double value;
  bmad_common_struct_get_synch_rad_scale(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_synch_rad_scale(double value) {
  bmad_common_struct_set_synch_rad_scale(fortran_ptr_, value);
}
double BmadCommonProxy::sad_eps_scale() const {
  double value;
  bmad_common_struct_get_sad_eps_scale(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_sad_eps_scale(double value) {
  bmad_common_struct_set_sad_eps_scale(fortran_ptr_, value);
}
double BmadCommonProxy::sad_amp_max() const {
  double value;
  bmad_common_struct_get_sad_amp_max(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_sad_amp_max(double value) {
  bmad_common_struct_set_sad_amp_max(fortran_ptr_, value);
}
int BmadCommonProxy::sad_n_div_max() const {
  int value;
  bmad_common_struct_get_sad_n_div_max(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_sad_n_div_max(int value) {
  bmad_common_struct_set_sad_n_div_max(fortran_ptr_, value);
}
int BmadCommonProxy::taylor_order() const {
  int value;
  bmad_common_struct_get_taylor_order(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_taylor_order(int value) {
  bmad_common_struct_set_taylor_order(fortran_ptr_, value);
}
int BmadCommonProxy::runge_kutta_order() const {
  int value;
  bmad_common_struct_get_runge_kutta_order(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_runge_kutta_order(int value) {
  bmad_common_struct_set_runge_kutta_order(fortran_ptr_, value);
}
int BmadCommonProxy::default_integ_order() const {
  int value;
  bmad_common_struct_get_default_integ_order(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_default_integ_order(int value) {
  bmad_common_struct_set_default_integ_order(fortran_ptr_, value);
}
int BmadCommonProxy::max_num_runge_kutta_step() const {
  int value;
  bmad_common_struct_get_max_num_runge_kutta_step(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_max_num_runge_kutta_step(int value) {
  bmad_common_struct_set_max_num_runge_kutta_step(fortran_ptr_, value);
}
bool BmadCommonProxy::rf_phase_below_transition_ref() const {
  bool value;
  bmad_common_struct_get_rf_phase_below_transition_ref(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_rf_phase_below_transition_ref(bool value) {
  bmad_common_struct_set_rf_phase_below_transition_ref(fortran_ptr_, value);
}
bool BmadCommonProxy::sr_wakes_on() const {
  bool value;
  bmad_common_struct_get_sr_wakes_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_sr_wakes_on(bool value) {
  bmad_common_struct_set_sr_wakes_on(fortran_ptr_, value);
}
bool BmadCommonProxy::lr_wakes_on() const {
  bool value;
  bmad_common_struct_get_lr_wakes_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_lr_wakes_on(bool value) {
  bmad_common_struct_set_lr_wakes_on(fortran_ptr_, value);
}
bool BmadCommonProxy::auto_bookkeeper() const {
  bool value;
  bmad_common_struct_get_auto_bookkeeper(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_auto_bookkeeper(bool value) {
  bmad_common_struct_set_auto_bookkeeper(fortran_ptr_, value);
}
bool BmadCommonProxy::high_energy_space_charge_on() const {
  bool value;
  bmad_common_struct_get_high_energy_space_charge_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_high_energy_space_charge_on(bool value) {
  bmad_common_struct_set_high_energy_space_charge_on(fortran_ptr_, value);
}
bool BmadCommonProxy::csr_and_space_charge_on() const {
  bool value;
  bmad_common_struct_get_csr_and_space_charge_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_csr_and_space_charge_on(bool value) {
  bmad_common_struct_set_csr_and_space_charge_on(fortran_ptr_, value);
}
bool BmadCommonProxy::spin_tracking_on() const {
  bool value;
  bmad_common_struct_get_spin_tracking_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_spin_tracking_on(bool value) {
  bmad_common_struct_set_spin_tracking_on(fortran_ptr_, value);
}
bool BmadCommonProxy::spin_sokolov_ternov_flipping_on() const {
  bool value;
  bmad_common_struct_get_spin_sokolov_ternov_flipping_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_spin_sokolov_ternov_flipping_on(bool value) {
  bmad_common_struct_set_spin_sokolov_ternov_flipping_on(fortran_ptr_, value);
}
bool BmadCommonProxy::radiation_damping_on() const {
  bool value;
  bmad_common_struct_get_radiation_damping_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_radiation_damping_on(bool value) {
  bmad_common_struct_set_radiation_damping_on(fortran_ptr_, value);
}
bool BmadCommonProxy::radiation_zero_average() const {
  bool value;
  bmad_common_struct_get_radiation_zero_average(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_radiation_zero_average(bool value) {
  bmad_common_struct_set_radiation_zero_average(fortran_ptr_, value);
}
bool BmadCommonProxy::radiation_fluctuations_on() const {
  bool value;
  bmad_common_struct_get_radiation_fluctuations_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_radiation_fluctuations_on(bool value) {
  bmad_common_struct_set_radiation_fluctuations_on(fortran_ptr_, value);
}
bool BmadCommonProxy::conserve_taylor_maps() const {
  bool value;
  bmad_common_struct_get_conserve_taylor_maps(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_conserve_taylor_maps(bool value) {
  bmad_common_struct_set_conserve_taylor_maps(fortran_ptr_, value);
}
bool BmadCommonProxy::absolute_time_tracking() const {
  bool value;
  bmad_common_struct_get_absolute_time_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_absolute_time_tracking(bool value) {
  bmad_common_struct_set_absolute_time_tracking(fortran_ptr_, value);
}
bool BmadCommonProxy::absolute_time_ref_shift() const {
  bool value;
  bmad_common_struct_get_absolute_time_ref_shift(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_absolute_time_ref_shift(bool value) {
  bmad_common_struct_set_absolute_time_ref_shift(fortran_ptr_, value);
}
bool BmadCommonProxy::convert_to_kinetic_momentum() const {
  bool value;
  bmad_common_struct_get_convert_to_kinetic_momentum(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_convert_to_kinetic_momentum(bool value) {
  bmad_common_struct_set_convert_to_kinetic_momentum(fortran_ptr_, value);
}
bool BmadCommonProxy::normalize_twiss() const {
  bool value;
  bmad_common_struct_get_normalize_twiss(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_normalize_twiss(bool value) {
  bmad_common_struct_set_normalize_twiss(fortran_ptr_, value);
}
bool BmadCommonProxy::aperture_limit_on() const {
  bool value;
  bmad_common_struct_get_aperture_limit_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_aperture_limit_on(bool value) {
  bmad_common_struct_set_aperture_limit_on(fortran_ptr_, value);
}
bool BmadCommonProxy::spin_n0_direction_user_set() const {
  bool value;
  bmad_common_struct_get_spin_n0_direction_user_set(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_spin_n0_direction_user_set(bool value) {
  bmad_common_struct_set_spin_n0_direction_user_set(fortran_ptr_, value);
}
bool BmadCommonProxy::debug() const {
  bool value;
  bmad_common_struct_get_debug(fortran_ptr_, &value);
  return value;
}
void BmadCommonProxy::set_debug(bool value) {
  bmad_common_struct_set_debug(fortran_ptr_, value);
}
double RadInt1Proxy::i0() const {
  double value;
  rad_int1_struct_get_i0(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i0(double value) {
  rad_int1_struct_set_i0(fortran_ptr_, value);
}
double RadInt1Proxy::i1() const {
  double value;
  rad_int1_struct_get_i1(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i1(double value) {
  rad_int1_struct_set_i1(fortran_ptr_, value);
}
double RadInt1Proxy::i2() const {
  double value;
  rad_int1_struct_get_i2(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i2(double value) {
  rad_int1_struct_set_i2(fortran_ptr_, value);
}
double RadInt1Proxy::i3() const {
  double value;
  rad_int1_struct_get_i3(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i3(double value) {
  rad_int1_struct_set_i3(fortran_ptr_, value);
}
double RadInt1Proxy::i4a() const {
  double value;
  rad_int1_struct_get_i4a(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i4a(double value) {
  rad_int1_struct_set_i4a(fortran_ptr_, value);
}
double RadInt1Proxy::i4b() const {
  double value;
  rad_int1_struct_get_i4b(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i4b(double value) {
  rad_int1_struct_set_i4b(fortran_ptr_, value);
}
double RadInt1Proxy::i4z() const {
  double value;
  rad_int1_struct_get_i4z(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i4z(double value) {
  rad_int1_struct_set_i4z(fortran_ptr_, value);
}
double RadInt1Proxy::i5a() const {
  double value;
  rad_int1_struct_get_i5a(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i5a(double value) {
  rad_int1_struct_set_i5a(fortran_ptr_, value);
}
double RadInt1Proxy::i5b() const {
  double value;
  rad_int1_struct_get_i5b(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i5b(double value) {
  rad_int1_struct_set_i5b(fortran_ptr_, value);
}
double RadInt1Proxy::i6b() const {
  double value;
  rad_int1_struct_get_i6b(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_i6b(double value) {
  rad_int1_struct_set_i6b(fortran_ptr_, value);
}
double RadInt1Proxy::lin_i2_E4() const {
  double value;
  rad_int1_struct_get_lin_i2_E4(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_i2_E4(double value) {
  rad_int1_struct_set_lin_i2_E4(fortran_ptr_, value);
}
double RadInt1Proxy::lin_i3_E7() const {
  double value;
  rad_int1_struct_get_lin_i3_E7(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_i3_E7(double value) {
  rad_int1_struct_set_lin_i3_E7(fortran_ptr_, value);
}
double RadInt1Proxy::lin_i5a_E6() const {
  double value;
  rad_int1_struct_get_lin_i5a_E6(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_i5a_E6(double value) {
  rad_int1_struct_set_lin_i5a_E6(fortran_ptr_, value);
}
double RadInt1Proxy::lin_i5b_E6() const {
  double value;
  rad_int1_struct_get_lin_i5b_E6(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_i5b_E6(double value) {
  rad_int1_struct_set_lin_i5b_E6(fortran_ptr_, value);
}
double RadInt1Proxy::lin_norm_emit_a() const {
  double value;
  rad_int1_struct_get_lin_norm_emit_a(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_norm_emit_a(double value) {
  rad_int1_struct_set_lin_norm_emit_a(fortran_ptr_, value);
}
double RadInt1Proxy::lin_norm_emit_b() const {
  double value;
  rad_int1_struct_get_lin_norm_emit_b(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_norm_emit_b(double value) {
  rad_int1_struct_set_lin_norm_emit_b(fortran_ptr_, value);
}
double RadInt1Proxy::lin_sig_E() const {
  double value;
  rad_int1_struct_get_lin_sig_E(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_lin_sig_E(double value) {
  rad_int1_struct_set_lin_sig_E(fortran_ptr_, value);
}
double RadInt1Proxy::n_steps() const {
  double value;
  rad_int1_struct_get_n_steps(fortran_ptr_, &value);
  return value;
}
void RadInt1Proxy::set_n_steps(double value) {
  rad_int1_struct_set_n_steps(fortran_ptr_, value);
}
RadInt1ProxyArray1D RadIntBranchProxy::ele() const {
  return ProxyHelpers::get_type_array_1d<RadInt1ProxyArray1D>(
      fortran_ptr_, rad_int_branch_struct_get_ele_info);
}
RadIntBranchProxyArray1D RadIntAllEleProxy::branch() const {
  return ProxyHelpers::get_type_array_1d<RadIntBranchProxyArray1D>(
      fortran_ptr_, rad_int_all_ele_struct_get_branch_info);
}
double RfStairStepProxy::E_tot0() const {
  double value;
  rf_stair_step_struct_get_E_tot0(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_E_tot0(double value) {
  rf_stair_step_struct_set_E_tot0(fortran_ptr_, value);
}
double RfStairStepProxy::E_tot1() const {
  double value;
  rf_stair_step_struct_get_E_tot1(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_E_tot1(double value) {
  rf_stair_step_struct_set_E_tot1(fortran_ptr_, value);
}
double RfStairStepProxy::p0c() const {
  double value;
  rf_stair_step_struct_get_p0c(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_p0c(double value) {
  rf_stair_step_struct_set_p0c(fortran_ptr_, value);
}
double RfStairStepProxy::p1c() const {
  double value;
  rf_stair_step_struct_get_p1c(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_p1c(double value) {
  rf_stair_step_struct_set_p1c(fortran_ptr_, value);
}
double RfStairStepProxy::scale() const {
  double value;
  rf_stair_step_struct_get_scale(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_scale(double value) {
  rf_stair_step_struct_set_scale(fortran_ptr_, value);
}
double RfStairStepProxy::time() const {
  double value;
  rf_stair_step_struct_get_time(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_time(double value) {
  rf_stair_step_struct_set_time(fortran_ptr_, value);
}
double RfStairStepProxy::s0() const {
  double value;
  rf_stair_step_struct_get_s0(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_s0(double value) {
  rf_stair_step_struct_set_s0(fortran_ptr_, value);
}
double RfStairStepProxy::s() const {
  double value;
  rf_stair_step_struct_get_s(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_s(double value) {
  rf_stair_step_struct_set_s(fortran_ptr_, value);
}
int RfStairStepProxy::ix_step() const {
  int value;
  rf_stair_step_struct_get_ix_step(fortran_ptr_, &value);
  return value;
}
void RfStairStepProxy::set_ix_step(int value) {
  rf_stair_step_struct_set_ix_step(fortran_ptr_, value);
}
RfStairStepProxyArray1D RfEleProxy::steps() const {
  return ProxyHelpers::get_type_array_1d<RfStairStepProxyArray1D>(
      fortran_ptr_, rf_ele_struct_get_steps_info);
}
double RfEleProxy::ds_step() const {
  double value;
  rf_ele_struct_get_ds_step(fortran_ptr_, &value);
  return value;
}
void RfEleProxy::set_ds_step(double value) {
  rf_ele_struct_set_ds_step(fortran_ptr_, value);
}
std::string EleProxy::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void EleProxy::set_name(const std::string& value) {
  ele_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleProxy::type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void EleProxy::set_type(const std::string& value) {
  ele_struct_set_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleProxy::alias() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_alias_info);
  return std::string(arr.data(), arr.size());
}
void EleProxy::set_alias(const std::string& value) {
  ele_struct_set_alias(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleProxy::component_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, ele_struct_get_component_name_info);
  return std::string(arr.data(), arr.size());
}
void EleProxy::set_component_name(const std::string& value) {
  ele_struct_set_component_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleProxy::descrip() const {
  return ProxyHelpers::get_string(fortran_ptr_, ele_struct_get_descrip_info);
}
void EleProxy::set_descrip(const std::string& value) {
  ele_struct_set_descrip(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TwissProxy EleProxy::a() const {
  void* ptr;
  ele_struct_get_a(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void EleProxy::set_a(const TwissProxy& src) {
  ele_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy EleProxy::b() const {
  void* ptr;
  ele_struct_get_b(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void EleProxy::set_b(const TwissProxy& src) {
  ele_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy EleProxy::z() const {
  void* ptr;
  ele_struct_get_z(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void EleProxy::set_z(const TwissProxy& src) {
  ele_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
XyDispProxy EleProxy::x() const {
  void* ptr;
  ele_struct_get_x(fortran_ptr_, &ptr);
  return XyDispProxy(ptr);
}
void EleProxy::set_x(const XyDispProxy& src) {
  ele_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
XyDispProxy EleProxy::y() const {
  void* ptr;
  ele_struct_get_y(fortran_ptr_, &ptr);
  return XyDispProxy(ptr);
}
void EleProxy::set_y(const XyDispProxy& src) {
  ele_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<AcKickerProxy> EleProxy::ac_kick() const {
  void* ptr;
  ele_struct_get_ac_kick(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return AcKickerProxy(ptr);
}
void EleProxy::set_ac_kick(const AcKickerProxy& src) {
  ele_struct_set_ac_kick(fortran_ptr_, src.get_fortran_ptr());
}
BookkeepingStateProxy EleProxy::bookkeeping_state() const {
  void* ptr;
  ele_struct_get_bookkeeping_state(fortran_ptr_, &ptr);
  return BookkeepingStateProxy(ptr);
}
void EleProxy::set_bookkeeping_state(const BookkeepingStateProxy& src) {
  ele_struct_set_bookkeeping_state(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<BranchProxy> EleProxy::branch() const {
  void* ptr;
  ele_struct_get_branch(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return BranchProxy(ptr);
}
void EleProxy::set_branch(const BranchProxy& src) {
  ele_struct_set_branch(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<ControllerProxy> EleProxy::control() const {
  void* ptr;
  ele_struct_get_control(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ControllerProxy(ptr);
}
void EleProxy::set_control(const ControllerProxy& src) {
  ele_struct_set_control(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<RfEleProxy> EleProxy::rf() const {
  void* ptr;
  ele_struct_get_rf(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return RfEleProxy(ptr);
}
void EleProxy::set_rf(const RfEleProxy& src) {
  ele_struct_set_rf(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<EleProxy> EleProxy::lord() const {
  void* ptr;
  ele_struct_get_lord(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleProxy(ptr);
}
void EleProxy::set_lord(const EleProxy& src) {
  ele_struct_set_lord(fortran_ptr_, src.get_fortran_ptr());
}
FloorPositionProxy EleProxy::floor() const {
  void* ptr;
  ele_struct_get_floor(fortran_ptr_, &ptr);
  return FloorPositionProxy(ptr);
}
void EleProxy::set_floor(const FloorPositionProxy& src) {
  ele_struct_set_floor(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<HighEnergySpaceChargeProxy> EleProxy::high_energy_space_charge()
    const {
  void* ptr;
  ele_struct_get_high_energy_space_charge(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return HighEnergySpaceChargeProxy(ptr);
}
void EleProxy::set_high_energy_space_charge(
    const HighEnergySpaceChargeProxy& src) {
  ele_struct_set_high_energy_space_charge(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<Mode3Proxy> EleProxy::mode3() const {
  void* ptr;
  ele_struct_get_mode3(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Mode3Proxy(ptr);
}
void EleProxy::set_mode3(const Mode3Proxy& src) {
  ele_struct_set_mode3(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<PhotonElementProxy> EleProxy::photon() const {
  void* ptr;
  ele_struct_get_photon(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return PhotonElementProxy(ptr);
}
void EleProxy::set_photon(const PhotonElementProxy& src) {
  ele_struct_set_photon(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<RadMapEleProxy> EleProxy::rad_map() const {
  void* ptr;
  ele_struct_get_rad_map(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return RadMapEleProxy(ptr);
}
void EleProxy::set_rad_map(const RadMapEleProxy& src) {
  ele_struct_set_rad_map(fortran_ptr_, src.get_fortran_ptr());
}
TaylorProxyArray1D EleProxy::taylor() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, ele_struct_get_taylor_info);
}
FArray1D<double> EleProxy::spin_taylor_ref_orb_in() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_spin_taylor_ref_orb_in_info);
}
TaylorProxyArray1D EleProxy::spin_taylor() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, ele_struct_get_spin_taylor_info);
}
std::optional<WakeProxy> EleProxy::wake() const {
  void* ptr;
  ele_struct_get_wake(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return WakeProxy(ptr);
}
void EleProxy::set_wake(const WakeProxy& src) {
  ele_struct_set_wake(fortran_ptr_, src.get_fortran_ptr());
}
Wall3dProxyArray1D EleProxy::wall3d() const {
  return ProxyHelpers::get_type_array_1d<Wall3dProxyArray1D>(
      fortran_ptr_, ele_struct_get_wall3d_info);
}
CartesianMapProxyArray1D EleProxy::cartesian_map() const {
  return ProxyHelpers::get_type_array_1d<CartesianMapProxyArray1D>(
      fortran_ptr_, ele_struct_get_cartesian_map_info);
}
CylindricalMapProxyArray1D EleProxy::cylindrical_map() const {
  return ProxyHelpers::get_type_array_1d<CylindricalMapProxyArray1D>(
      fortran_ptr_, ele_struct_get_cylindrical_map_info);
}
GenGradMapProxyArray1D EleProxy::gen_grad_map() const {
  return ProxyHelpers::get_type_array_1d<GenGradMapProxyArray1D>(
      fortran_ptr_, ele_struct_get_gen_grad_map_info);
}
GridFieldProxyArray1D EleProxy::grid_field() const {
  return ProxyHelpers::get_type_array_1d<GridFieldProxyArray1D>(
      fortran_ptr_, ele_struct_get_grid_field_info);
}
CoordProxy EleProxy::map_ref_orb_in() const {
  void* ptr;
  ele_struct_get_map_ref_orb_in(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void EleProxy::set_map_ref_orb_in(const CoordProxy& src) {
  ele_struct_set_map_ref_orb_in(fortran_ptr_, src.get_fortran_ptr());
}
CoordProxy EleProxy::map_ref_orb_out() const {
  void* ptr;
  ele_struct_get_map_ref_orb_out(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void EleProxy::set_map_ref_orb_out(const CoordProxy& src) {
  ele_struct_set_map_ref_orb_out(fortran_ptr_, src.get_fortran_ptr());
}
CoordProxy EleProxy::time_ref_orb_in() const {
  void* ptr;
  ele_struct_get_time_ref_orb_in(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void EleProxy::set_time_ref_orb_in(const CoordProxy& src) {
  ele_struct_set_time_ref_orb_in(fortran_ptr_, src.get_fortran_ptr());
}
CoordProxy EleProxy::time_ref_orb_out() const {
  void* ptr;
  ele_struct_get_time_ref_orb_out(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void EleProxy::set_time_ref_orb_out(const CoordProxy& src) {
  ele_struct_set_time_ref_orb_out(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> EleProxy::value() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_value_info);
}
FArray1D<double> EleProxy::old_value() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_old_value_info);
}
FArray2D<double> EleProxy::spin_q() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, ele_struct_get_spin_q_info);
}
FArray1D<double> EleProxy::vec0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_vec0_info);
}
FArray2D<double> EleProxy::mat6() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, ele_struct_get_mat6_info);
}
FArray2D<double> EleProxy::c_mat() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, ele_struct_get_c_mat_info);
}
FArray2D<double> EleProxy::dc_mat_dpz() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, ele_struct_get_dc_mat_dpz_info);
}
double EleProxy::gamma_c() const {
  double value;
  ele_struct_get_gamma_c(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_gamma_c(double value) {
  ele_struct_set_gamma_c(fortran_ptr_, value);
}
double EleProxy::s_start() const {
  double value;
  ele_struct_get_s_start(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_s_start(double value) {
  ele_struct_set_s_start(fortran_ptr_, value);
}
double EleProxy::s() const {
  double value;
  ele_struct_get_s(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_s(double value) {
  ele_struct_set_s(fortran_ptr_, value);
}
double EleProxy::ref_time() const {
  double value;
  ele_struct_get_ref_time(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ref_time(double value) {
  ele_struct_set_ref_time(fortran_ptr_, value);
}
FArray1D<double> EleProxy::a_pole() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_a_pole_info);
}
FArray1D<double> EleProxy::b_pole() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_b_pole_info);
}
FArray1D<double> EleProxy::a_pole_elec() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_a_pole_elec_info);
}
FArray1D<double> EleProxy::b_pole_elec() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_b_pole_elec_info);
}
FArray1D<double> EleProxy::custom() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ele_struct_get_custom_info);
}
FArray3D<double> EleProxy::r() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, ele_struct_get_r_info);
}
int EleProxy::key() const {
  int value;
  ele_struct_get_key(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_key(int value) {
  ele_struct_set_key(fortran_ptr_, value);
}
int EleProxy::sub_key() const {
  int value;
  ele_struct_get_sub_key(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_sub_key(int value) {
  ele_struct_set_sub_key(fortran_ptr_, value);
}
int EleProxy::ix_ele() const {
  int value;
  ele_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ix_ele(int value) {
  ele_struct_set_ix_ele(fortran_ptr_, value);
}
int EleProxy::ix_branch() const {
  int value;
  ele_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ix_branch(int value) {
  ele_struct_set_ix_branch(fortran_ptr_, value);
}
int EleProxy::lord_status() const {
  int value;
  ele_struct_get_lord_status(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_lord_status(int value) {
  ele_struct_set_lord_status(fortran_ptr_, value);
}
int EleProxy::n_slave() const {
  int value;
  ele_struct_get_n_slave(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_n_slave(int value) {
  ele_struct_set_n_slave(fortran_ptr_, value);
}
int EleProxy::n_slave_field() const {
  int value;
  ele_struct_get_n_slave_field(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_n_slave_field(int value) {
  ele_struct_set_n_slave_field(fortran_ptr_, value);
}
int EleProxy::ix1_slave() const {
  int value;
  ele_struct_get_ix1_slave(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ix1_slave(int value) {
  ele_struct_set_ix1_slave(fortran_ptr_, value);
}
int EleProxy::slave_status() const {
  int value;
  ele_struct_get_slave_status(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_slave_status(int value) {
  ele_struct_set_slave_status(fortran_ptr_, value);
}
int EleProxy::n_lord() const {
  int value;
  ele_struct_get_n_lord(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_n_lord(int value) {
  ele_struct_set_n_lord(fortran_ptr_, value);
}
int EleProxy::n_lord_field() const {
  int value;
  ele_struct_get_n_lord_field(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_n_lord_field(int value) {
  ele_struct_set_n_lord_field(fortran_ptr_, value);
}
int EleProxy::n_lord_ramper() const {
  int value;
  ele_struct_get_n_lord_ramper(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_n_lord_ramper(int value) {
  ele_struct_set_n_lord_ramper(fortran_ptr_, value);
}
int EleProxy::ic1_lord() const {
  int value;
  ele_struct_get_ic1_lord(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ic1_lord(int value) {
  ele_struct_set_ic1_lord(fortran_ptr_, value);
}
int EleProxy::ix_pointer() const {
  int value;
  ele_struct_get_ix_pointer(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ix_pointer(int value) {
  ele_struct_set_ix_pointer(fortran_ptr_, value);
}
int EleProxy::ixx() const {
  int value;
  ele_struct_get_ixx(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ixx(int value) {
  ele_struct_set_ixx(fortran_ptr_, value);
}
int EleProxy::iyy() const {
  int value;
  ele_struct_get_iyy(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_iyy(int value) {
  ele_struct_set_iyy(fortran_ptr_, value);
}
int EleProxy::izz() const {
  int value;
  ele_struct_get_izz(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_izz(int value) {
  ele_struct_set_izz(fortran_ptr_, value);
}
int EleProxy::mat6_calc_method() const {
  int value;
  ele_struct_get_mat6_calc_method(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_mat6_calc_method(int value) {
  ele_struct_set_mat6_calc_method(fortran_ptr_, value);
}
int EleProxy::tracking_method() const {
  int value;
  ele_struct_get_tracking_method(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_tracking_method(int value) {
  ele_struct_set_tracking_method(fortran_ptr_, value);
}
int EleProxy::spin_tracking_method() const {
  int value;
  ele_struct_get_spin_tracking_method(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_spin_tracking_method(int value) {
  ele_struct_set_spin_tracking_method(fortran_ptr_, value);
}
int EleProxy::csr_method() const {
  int value;
  ele_struct_get_csr_method(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_csr_method(int value) {
  ele_struct_set_csr_method(fortran_ptr_, value);
}
int EleProxy::space_charge_method() const {
  int value;
  ele_struct_get_space_charge_method(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_space_charge_method(int value) {
  ele_struct_set_space_charge_method(fortran_ptr_, value);
}
int EleProxy::ptc_integration_type() const {
  int value;
  ele_struct_get_ptc_integration_type(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ptc_integration_type(int value) {
  ele_struct_set_ptc_integration_type(fortran_ptr_, value);
}
int EleProxy::field_calc() const {
  int value;
  ele_struct_get_field_calc(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_field_calc(int value) {
  ele_struct_set_field_calc(fortran_ptr_, value);
}
int EleProxy::aperture_at() const {
  int value;
  ele_struct_get_aperture_at(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_aperture_at(int value) {
  ele_struct_set_aperture_at(fortran_ptr_, value);
}
int EleProxy::aperture_type() const {
  int value;
  ele_struct_get_aperture_type(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_aperture_type(int value) {
  ele_struct_set_aperture_type(fortran_ptr_, value);
}
int EleProxy::ref_species() const {
  int value;
  ele_struct_get_ref_species(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_ref_species(int value) {
  ele_struct_set_ref_species(fortran_ptr_, value);
}
int EleProxy::orientation() const {
  int value;
  ele_struct_get_orientation(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_orientation(int value) {
  ele_struct_set_orientation(fortran_ptr_, value);
}
bool EleProxy::symplectify() const {
  bool value;
  ele_struct_get_symplectify(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_symplectify(bool value) {
  ele_struct_set_symplectify(fortran_ptr_, value);
}
bool EleProxy::mode_flip() const {
  bool value;
  ele_struct_get_mode_flip(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_mode_flip(bool value) {
  ele_struct_set_mode_flip(fortran_ptr_, value);
}
bool EleProxy::multipoles_on() const {
  bool value;
  ele_struct_get_multipoles_on(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_multipoles_on(bool value) {
  ele_struct_set_multipoles_on(fortran_ptr_, value);
}
bool EleProxy::scale_multipoles() const {
  bool value;
  ele_struct_get_scale_multipoles(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_scale_multipoles(bool value) {
  ele_struct_set_scale_multipoles(fortran_ptr_, value);
}
bool EleProxy::taylor_map_includes_offsets() const {
  bool value;
  ele_struct_get_taylor_map_includes_offsets(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_taylor_map_includes_offsets(bool value) {
  ele_struct_set_taylor_map_includes_offsets(fortran_ptr_, value);
}
bool EleProxy::field_master() const {
  bool value;
  ele_struct_get_field_master(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_field_master(bool value) {
  ele_struct_set_field_master(fortran_ptr_, value);
}
bool EleProxy::is_on() const {
  bool value;
  ele_struct_get_is_on(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_is_on(bool value) {
  ele_struct_set_is_on(fortran_ptr_, value);
}
bool EleProxy::logic() const {
  bool value;
  ele_struct_get_logic(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_logic(bool value) {
  ele_struct_set_logic(fortran_ptr_, value);
}
bool EleProxy::bmad_logic() const {
  bool value;
  ele_struct_get_bmad_logic(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_bmad_logic(bool value) {
  ele_struct_set_bmad_logic(fortran_ptr_, value);
}
bool EleProxy::select() const {
  bool value;
  ele_struct_get_select(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_select(bool value) {
  ele_struct_set_select(fortran_ptr_, value);
}
bool EleProxy::offset_moves_aperture() const {
  bool value;
  ele_struct_get_offset_moves_aperture(fortran_ptr_, &value);
  return value;
}
void EleProxy::set_offset_moves_aperture(bool value) {
  ele_struct_set_offset_moves_aperture(fortran_ptr_, value);
}
std::complex<double> ComplexTaylorTermProxy::coef() const {
  std::complex<double> c_value;
  complex_taylor_term_struct_get_coef(fortran_ptr_, &c_value);
  return c_value;
}
void ComplexTaylorTermProxy::set_coef(std::complex<double> value) {
  complex_taylor_term_struct_set_coef(fortran_ptr_, value);
}
FArray1D<int> ComplexTaylorTermProxy::expn() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, complex_taylor_term_struct_get_expn_info);
}
std::complex<double> ComplexTaylorProxy::ref() const {
  std::complex<double> c_value;
  complex_taylor_struct_get_ref(fortran_ptr_, &c_value);
  return c_value;
}
void ComplexTaylorProxy::set_ref(std::complex<double> value) {
  complex_taylor_struct_set_ref(fortran_ptr_, value);
}
ComplexTaylorTermProxyArray1D ComplexTaylorProxy::term() const {
  return ProxyHelpers::get_type_array_1d<ComplexTaylorTermProxyArray1D>(
      fortran_ptr_, complex_taylor_struct_get_term_info);
}
std::string BranchProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, branch_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void BranchProxy::set_name(const std::string& value) {
  branch_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int BranchProxy::ix_branch() const {
  int value;
  branch_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_ix_branch(int value) {
  branch_struct_set_ix_branch(fortran_ptr_, value);
}
int BranchProxy::ix_from_branch() const {
  int value;
  branch_struct_get_ix_from_branch(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_ix_from_branch(int value) {
  branch_struct_set_ix_from_branch(fortran_ptr_, value);
}
int BranchProxy::ix_from_ele() const {
  int value;
  branch_struct_get_ix_from_ele(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_ix_from_ele(int value) {
  branch_struct_set_ix_from_ele(fortran_ptr_, value);
}
int BranchProxy::ix_to_ele() const {
  int value;
  branch_struct_get_ix_to_ele(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_ix_to_ele(int value) {
  branch_struct_set_ix_to_ele(fortran_ptr_, value);
}
int BranchProxy::ix_fixer() const {
  int value;
  branch_struct_get_ix_fixer(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_ix_fixer(int value) {
  branch_struct_set_ix_fixer(fortran_ptr_, value);
}
int BranchProxy::n_ele_track() const {
  int value;
  branch_struct_get_n_ele_track(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_n_ele_track(int value) {
  branch_struct_set_n_ele_track(fortran_ptr_, value);
}
int BranchProxy::n_ele_max() const {
  int value;
  branch_struct_get_n_ele_max(fortran_ptr_, &value);
  return value;
}
void BranchProxy::set_n_ele_max(int value) {
  branch_struct_set_n_ele_max(fortran_ptr_, value);
}
std::optional<LatProxy> BranchProxy::lat() const {
  void* ptr;
  branch_struct_get_lat(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return LatProxy(ptr);
}
void BranchProxy::set_lat(const LatProxy& src) {
  branch_struct_set_lat(fortran_ptr_, src.get_fortran_ptr());
}
ModeInfoProxy BranchProxy::a() const {
  void* ptr;
  branch_struct_get_a(fortran_ptr_, &ptr);
  return ModeInfoProxy(ptr);
}
void BranchProxy::set_a(const ModeInfoProxy& src) {
  branch_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
ModeInfoProxy BranchProxy::b() const {
  void* ptr;
  branch_struct_get_b(fortran_ptr_, &ptr);
  return ModeInfoProxy(ptr);
}
void BranchProxy::set_b(const ModeInfoProxy& src) {
  branch_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
ModeInfoProxy BranchProxy::z() const {
  void* ptr;
  branch_struct_get_z(fortran_ptr_, &ptr);
  return ModeInfoProxy(ptr);
}
void BranchProxy::set_z(const ModeInfoProxy& src) {
  branch_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
EleProxyArray1D BranchProxy::ele() const {
  return ProxyHelpers::get_type_array_1d<EleProxyArray1D>(
      fortran_ptr_, branch_struct_get_ele_info);
}
LatParamProxy BranchProxy::param() const {
  void* ptr;
  branch_struct_get_param(fortran_ptr_, &ptr);
  return LatParamProxy(ptr);
}
void BranchProxy::set_param(const LatParamProxy& src) {
  branch_struct_set_param(fortran_ptr_, src.get_fortran_ptr());
}
CoordProxy BranchProxy::particle_start() const {
  void* ptr;
  branch_struct_get_particle_start(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void BranchProxy::set_particle_start(const CoordProxy& src) {
  branch_struct_set_particle_start(fortran_ptr_, src.get_fortran_ptr());
}
Wall3dProxyArray1D BranchProxy::wall3d() const {
  return ProxyHelpers::get_type_array_1d<Wall3dProxyArray1D>(
      fortran_ptr_, branch_struct_get_wall3d_info);
}
std::string LatProxy::use_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, lat_struct_get_use_name_info);
  return std::string(arr.data(), arr.size());
}
void LatProxy::set_use_name(const std::string& value) {
  lat_struct_set_use_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatProxy::lattice() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, lat_struct_get_lattice_info);
  return std::string(arr.data(), arr.size());
}
void LatProxy::set_lattice(const std::string& value) {
  lat_struct_set_lattice(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatProxy::machine() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, lat_struct_get_machine_info);
  return std::string(arr.data(), arr.size());
}
void LatProxy::set_machine(const std::string& value) {
  lat_struct_set_machine(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatProxy::input_file_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, lat_struct_get_input_file_name_info);
  return std::string(arr.data(), arr.size());
}
void LatProxy::set_input_file_name(const std::string& value) {
  lat_struct_set_input_file_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatProxy::title() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, lat_struct_get_title_info);
  return std::string(arr.data(), arr.size());
}
void LatProxy::set_title(const std::string& value) {
  lat_struct_set_title(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D LatProxy::print_str() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, lat_struct_get_print_str_info);
}
ExpressionAtomProxyArray1D LatProxy::constant() const {
  return ProxyHelpers::get_type_array_1d<ExpressionAtomProxyArray1D>(
      fortran_ptr_, lat_struct_get_constant_info);
}
std::optional<ModeInfoProxy> LatProxy::a() const {
  void* ptr;
  lat_struct_get_a(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ModeInfoProxy(ptr);
}
void LatProxy::set_a(const ModeInfoProxy& src) {
  lat_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<ModeInfoProxy> LatProxy::b() const {
  void* ptr;
  lat_struct_get_b(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ModeInfoProxy(ptr);
}
void LatProxy::set_b(const ModeInfoProxy& src) {
  lat_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<ModeInfoProxy> LatProxy::z() const {
  void* ptr;
  lat_struct_get_z(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ModeInfoProxy(ptr);
}
void LatProxy::set_z(const ModeInfoProxy& src) {
  lat_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<LatParamProxy> LatProxy::param() const {
  void* ptr;
  lat_struct_get_param(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return LatParamProxy(ptr);
}
void LatProxy::set_param(const LatParamProxy& src) {
  lat_struct_set_param(fortran_ptr_, src.get_fortran_ptr());
}
BookkeepingStateProxy LatProxy::lord_state() const {
  void* ptr;
  lat_struct_get_lord_state(fortran_ptr_, &ptr);
  return BookkeepingStateProxy(ptr);
}
void LatProxy::set_lord_state(const BookkeepingStateProxy& src) {
  lat_struct_set_lord_state(fortran_ptr_, src.get_fortran_ptr());
}
EleProxy LatProxy::ele_init() const {
  void* ptr;
  lat_struct_get_ele_init(fortran_ptr_, &ptr);
  return EleProxy(ptr);
}
void LatProxy::set_ele_init(const EleProxy& src) {
  lat_struct_set_ele_init(fortran_ptr_, src.get_fortran_ptr());
}
EleProxyArray1D LatProxy::ele() const {
  return ProxyHelpers::get_type_array_1d<EleProxyArray1D>(
      fortran_ptr_, lat_struct_get_ele_info);
}
BranchProxyArray1D LatProxy::branch() const {
  return ProxyHelpers::get_type_array_1d<BranchProxyArray1D>(
      fortran_ptr_, lat_struct_get_branch_info);
}
ControlProxyArray1D LatProxy::control() const {
  return ProxyHelpers::get_type_array_1d<ControlProxyArray1D>(
      fortran_ptr_, lat_struct_get_control_info);
}
std::optional<CoordProxy> LatProxy::particle_start() const {
  void* ptr;
  lat_struct_get_particle_start(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return CoordProxy(ptr);
}
void LatProxy::set_particle_start(const CoordProxy& src) {
  lat_struct_set_particle_start(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitProxy LatProxy::beam_init() const {
  void* ptr;
  lat_struct_get_beam_init(fortran_ptr_, &ptr);
  return BeamInitProxy(ptr);
}
void LatProxy::set_beam_init(const BeamInitProxy& src) {
  lat_struct_set_beam_init(fortran_ptr_, src.get_fortran_ptr());
}
PreTrackerProxy LatProxy::pre_tracker() const {
  void* ptr;
  lat_struct_get_pre_tracker(fortran_ptr_, &ptr);
  return PreTrackerProxy(ptr);
}
void LatProxy::set_pre_tracker(const PreTrackerProxy& src) {
  lat_struct_set_pre_tracker(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> LatProxy::custom() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, lat_struct_get_custom_info);
}
int LatProxy::version() const {
  int value;
  lat_struct_get_version(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_version(int value) {
  lat_struct_set_version(fortran_ptr_, value);
}
int* LatProxy::n_ele_track() const {
  int* ptr;
  lat_struct_get_n_ele_track(fortran_ptr_, &ptr);
  return ptr;
}
void LatProxy::set_n_ele_track(int value) {
  lat_struct_set_n_ele_track(fortran_ptr_, value);
}
int* LatProxy::n_ele_max() const {
  int* ptr;
  lat_struct_get_n_ele_max(fortran_ptr_, &ptr);
  return ptr;
}
void LatProxy::set_n_ele_max(int value) {
  lat_struct_set_n_ele_max(fortran_ptr_, value);
}
int LatProxy::n_control_max() const {
  int value;
  lat_struct_get_n_control_max(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_n_control_max(int value) {
  lat_struct_set_n_control_max(fortran_ptr_, value);
}
int LatProxy::n_ic_max() const {
  int value;
  lat_struct_get_n_ic_max(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_n_ic_max(int value) {
  lat_struct_set_n_ic_max(fortran_ptr_, value);
}
int LatProxy::input_taylor_order() const {
  int value;
  lat_struct_get_input_taylor_order(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_input_taylor_order(int value) {
  lat_struct_set_input_taylor_order(fortran_ptr_, value);
}
FArray1D<int> LatProxy::ic() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, lat_struct_get_ic_info);
}
int LatProxy::photon_type() const {
  int value;
  lat_struct_get_photon_type(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_photon_type(int value) {
  lat_struct_set_photon_type(fortran_ptr_, value);
}
int LatProxy::creation_hash() const {
  int value;
  lat_struct_get_creation_hash(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_creation_hash(int value) {
  lat_struct_set_creation_hash(fortran_ptr_, value);
}
int LatProxy::ramper_slave_bookkeeping() const {
  int value;
  lat_struct_get_ramper_slave_bookkeeping(fortran_ptr_, &value);
  return value;
}
void LatProxy::set_ramper_slave_bookkeeping(int value) {
  lat_struct_set_ramper_slave_bookkeeping(fortran_ptr_, value);
}
CoordProxyArray1D BunchProxy::particle() const {
  return ProxyHelpers::get_type_array_1d<CoordProxyArray1D>(
      fortran_ptr_, bunch_struct_get_particle_info);
}
FArray1D<int> BunchProxy::ix_z() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, bunch_struct_get_ix_z_info);
}
double BunchProxy::charge_tot() const {
  double value;
  bunch_struct_get_charge_tot(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_charge_tot(double value) {
  bunch_struct_set_charge_tot(fortran_ptr_, value);
}
double BunchProxy::charge_live() const {
  double value;
  bunch_struct_get_charge_live(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_charge_live(double value) {
  bunch_struct_set_charge_live(fortran_ptr_, value);
}
double BunchProxy::z_center() const {
  double value;
  bunch_struct_get_z_center(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_z_center(double value) {
  bunch_struct_set_z_center(fortran_ptr_, value);
}
double BunchProxy::t_center() const {
  double value;
  bunch_struct_get_t_center(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_t_center(double value) {
  bunch_struct_set_t_center(fortran_ptr_, value);
}
double BunchProxy::t0() const {
  double value;
  bunch_struct_get_t0(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_t0(double value) {
  bunch_struct_set_t0(fortran_ptr_, value);
}
bool BunchProxy::drift_between_t_and_s() const {
  bool value;
  bunch_struct_get_drift_between_t_and_s(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_drift_between_t_and_s(bool value) {
  bunch_struct_set_drift_between_t_and_s(fortran_ptr_, value);
}
int BunchProxy::ix_ele() const {
  int value;
  bunch_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_ix_ele(int value) {
  bunch_struct_set_ix_ele(fortran_ptr_, value);
}
int BunchProxy::ix_bunch() const {
  int value;
  bunch_struct_get_ix_bunch(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_ix_bunch(int value) {
  bunch_struct_set_ix_bunch(fortran_ptr_, value);
}
int BunchProxy::ix_turn() const {
  int value;
  bunch_struct_get_ix_turn(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_ix_turn(int value) {
  bunch_struct_set_ix_turn(fortran_ptr_, value);
}
int BunchProxy::n_live() const {
  int value;
  bunch_struct_get_n_live(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_n_live(int value) {
  bunch_struct_set_n_live(fortran_ptr_, value);
}
int BunchProxy::n_good() const {
  int value;
  bunch_struct_get_n_good(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_n_good(int value) {
  bunch_struct_set_n_good(fortran_ptr_, value);
}
int BunchProxy::n_bad() const {
  int value;
  bunch_struct_get_n_bad(fortran_ptr_, &value);
  return value;
}
void BunchProxy::set_n_bad(int value) {
  bunch_struct_set_n_bad(fortran_ptr_, value);
}
CoordProxy BunchParamsProxy::centroid() const {
  void* ptr;
  bunch_params_struct_get_centroid(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void BunchParamsProxy::set_centroid(const CoordProxy& src) {
  bunch_params_struct_set_centroid(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy BunchParamsProxy::x() const {
  void* ptr;
  bunch_params_struct_get_x(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void BunchParamsProxy::set_x(const TwissProxy& src) {
  bunch_params_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy BunchParamsProxy::y() const {
  void* ptr;
  bunch_params_struct_get_y(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void BunchParamsProxy::set_y(const TwissProxy& src) {
  bunch_params_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy BunchParamsProxy::z() const {
  void* ptr;
  bunch_params_struct_get_z(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void BunchParamsProxy::set_z(const TwissProxy& src) {
  bunch_params_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy BunchParamsProxy::a() const {
  void* ptr;
  bunch_params_struct_get_a(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void BunchParamsProxy::set_a(const TwissProxy& src) {
  bunch_params_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy BunchParamsProxy::b() const {
  void* ptr;
  bunch_params_struct_get_b(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void BunchParamsProxy::set_b(const TwissProxy& src) {
  bunch_params_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
TwissProxy BunchParamsProxy::c() const {
  void* ptr;
  bunch_params_struct_get_c(fortran_ptr_, &ptr);
  return TwissProxy(ptr);
}
void BunchParamsProxy::set_c(const TwissProxy& src) {
  bunch_params_struct_set_c(fortran_ptr_, src.get_fortran_ptr());
}
FArray2D<double> BunchParamsProxy::sigma() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, bunch_params_struct_get_sigma_info);
}
FArray1D<double> BunchParamsProxy::rel_max() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, bunch_params_struct_get_rel_max_info);
}
FArray1D<double> BunchParamsProxy::rel_min() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, bunch_params_struct_get_rel_min_info);
}
double BunchParamsProxy::s() const {
  double value;
  bunch_params_struct_get_s(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_s(double value) {
  bunch_params_struct_set_s(fortran_ptr_, value);
}
double BunchParamsProxy::t() const {
  double value;
  bunch_params_struct_get_t(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_t(double value) {
  bunch_params_struct_set_t(fortran_ptr_, value);
}
double BunchParamsProxy::sigma_t() const {
  double value;
  bunch_params_struct_get_sigma_t(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_sigma_t(double value) {
  bunch_params_struct_set_sigma_t(fortran_ptr_, value);
}
double BunchParamsProxy::charge_live() const {
  double value;
  bunch_params_struct_get_charge_live(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_charge_live(double value) {
  bunch_params_struct_set_charge_live(fortran_ptr_, value);
}
double BunchParamsProxy::charge_tot() const {
  double value;
  bunch_params_struct_get_charge_tot(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_charge_tot(double value) {
  bunch_params_struct_set_charge_tot(fortran_ptr_, value);
}
int BunchParamsProxy::n_particle_tot() const {
  int value;
  bunch_params_struct_get_n_particle_tot(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_n_particle_tot(int value) {
  bunch_params_struct_set_n_particle_tot(fortran_ptr_, value);
}
int BunchParamsProxy::n_particle_live() const {
  int value;
  bunch_params_struct_get_n_particle_live(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_n_particle_live(int value) {
  bunch_params_struct_set_n_particle_live(fortran_ptr_, value);
}
int BunchParamsProxy::n_particle_lost_in_ele() const {
  int value;
  bunch_params_struct_get_n_particle_lost_in_ele(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_n_particle_lost_in_ele(int value) {
  bunch_params_struct_set_n_particle_lost_in_ele(fortran_ptr_, value);
}
int BunchParamsProxy::n_good_steps() const {
  int value;
  bunch_params_struct_get_n_good_steps(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_n_good_steps(int value) {
  bunch_params_struct_set_n_good_steps(fortran_ptr_, value);
}
int BunchParamsProxy::n_bad_steps() const {
  int value;
  bunch_params_struct_get_n_bad_steps(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_n_bad_steps(int value) {
  bunch_params_struct_set_n_bad_steps(fortran_ptr_, value);
}
int BunchParamsProxy::ix_ele() const {
  int value;
  bunch_params_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_ix_ele(int value) {
  bunch_params_struct_set_ix_ele(fortran_ptr_, value);
}
int BunchParamsProxy::location() const {
  int value;
  bunch_params_struct_get_location(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_location(int value) {
  bunch_params_struct_set_location(fortran_ptr_, value);
}
bool BunchParamsProxy::twiss_valid() const {
  bool value;
  bunch_params_struct_get_twiss_valid(fortran_ptr_, &value);
  return value;
}
void BunchParamsProxy::set_twiss_valid(bool value) {
  bunch_params_struct_set_twiss_valid(fortran_ptr_, value);
}
BunchProxyArray1D BeamProxy::bunch() const {
  return ProxyHelpers::get_type_array_1d<BunchProxyArray1D>(
      fortran_ptr_, beam_struct_get_bunch_info);
}
double AperturePointProxy::x() const {
  double value;
  aperture_point_struct_get_x(fortran_ptr_, &value);
  return value;
}
void AperturePointProxy::set_x(double value) {
  aperture_point_struct_set_x(fortran_ptr_, value);
}
double AperturePointProxy::y() const {
  double value;
  aperture_point_struct_get_y(fortran_ptr_, &value);
  return value;
}
void AperturePointProxy::set_y(double value) {
  aperture_point_struct_set_y(fortran_ptr_, value);
}
int AperturePointProxy::plane() const {
  int value;
  aperture_point_struct_get_plane(fortran_ptr_, &value);
  return value;
}
void AperturePointProxy::set_plane(int value) {
  aperture_point_struct_set_plane(fortran_ptr_, value);
}
int AperturePointProxy::ix_ele() const {
  int value;
  aperture_point_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void AperturePointProxy::set_ix_ele(int value) {
  aperture_point_struct_set_ix_ele(fortran_ptr_, value);
}
int AperturePointProxy::i_turn() const {
  int value;
  aperture_point_struct_get_i_turn(fortran_ptr_, &value);
  return value;
}
void AperturePointProxy::set_i_turn(int value) {
  aperture_point_struct_set_i_turn(fortran_ptr_, value);
}
double ApertureParamProxy::min_angle() const {
  double value;
  aperture_param_struct_get_min_angle(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_min_angle(double value) {
  aperture_param_struct_set_min_angle(fortran_ptr_, value);
}
double ApertureParamProxy::max_angle() const {
  double value;
  aperture_param_struct_get_max_angle(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_max_angle(double value) {
  aperture_param_struct_set_max_angle(fortran_ptr_, value);
}
int ApertureParamProxy::n_angle() const {
  int value;
  aperture_param_struct_get_n_angle(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_n_angle(int value) {
  aperture_param_struct_set_n_angle(fortran_ptr_, value);
}
int ApertureParamProxy::n_turn() const {
  int value;
  aperture_param_struct_get_n_turn(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_n_turn(int value) {
  aperture_param_struct_set_n_turn(fortran_ptr_, value);
}
double ApertureParamProxy::x_init() const {
  double value;
  aperture_param_struct_get_x_init(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_x_init(double value) {
  aperture_param_struct_set_x_init(fortran_ptr_, value);
}
double ApertureParamProxy::y_init() const {
  double value;
  aperture_param_struct_get_y_init(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_y_init(double value) {
  aperture_param_struct_set_y_init(fortran_ptr_, value);
}
double ApertureParamProxy::rel_accuracy() const {
  double value;
  aperture_param_struct_get_rel_accuracy(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_rel_accuracy(double value) {
  aperture_param_struct_set_rel_accuracy(fortran_ptr_, value);
}
double ApertureParamProxy::abs_accuracy() const {
  double value;
  aperture_param_struct_get_abs_accuracy(fortran_ptr_, &value);
  return value;
}
void ApertureParamProxy::set_abs_accuracy(double value) {
  aperture_param_struct_set_abs_accuracy(fortran_ptr_, value);
}
std::string ApertureParamProxy::start_ele() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, aperture_param_struct_get_start_ele_info);
  return std::string(arr.data(), arr.size());
}
void ApertureParamProxy::set_start_ele(const std::string& value) {
  aperture_param_struct_set_start_ele(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
AperturePointProxyArray1D ApertureScanProxy::point() const {
  return ProxyHelpers::get_type_array_1d<AperturePointProxyArray1D>(
      fortran_ptr_, aperture_scan_struct_get_point_info);
}
CoordProxy ApertureScanProxy::ref_orb() const {
  void* ptr;
  aperture_scan_struct_get_ref_orb(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void ApertureScanProxy::set_ref_orb(const CoordProxy& src) {
  aperture_scan_struct_set_ref_orb(fortran_ptr_, src.get_fortran_ptr());
}
double ApertureScanProxy::pz_start() const {
  double value;
  aperture_scan_struct_get_pz_start(fortran_ptr_, &value);
  return value;
}
void ApertureScanProxy::set_pz_start(double value) {
  aperture_scan_struct_set_pz_start(fortran_ptr_, value);
}
std::optional<EleProxy> ElePointerProxy::ele() const {
  void* ptr;
  ele_pointer_struct_get_ele(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleProxy(ptr);
}
void ElePointerProxy::set_ele(const EleProxy& src) {
  ele_pointer_struct_set_ele(fortran_ptr_, src.get_fortran_ptr());
}
LatEleLocProxy ElePointerProxy::loc() const {
  void* ptr;
  ele_pointer_struct_get_loc(fortran_ptr_, &ptr);
  return LatEleLocProxy(ptr);
}
void ElePointerProxy::set_loc(const LatEleLocProxy& src) {
  ele_pointer_struct_set_loc(fortran_ptr_, src.get_fortran_ptr());
}
int ElePointerProxy::id() const {
  int value;
  ele_pointer_struct_get_id(fortran_ptr_, &value);
  return value;
}
void ElePointerProxy::set_id(int value) {
  ele_pointer_struct_set_id(fortran_ptr_, value);
}
std::string ExpressionTreeProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, expression_tree_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void ExpressionTreeProxy::set_name(const std::string& value) {
  expression_tree_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int ExpressionTreeProxy::type() const {
  int value;
  expression_tree_struct_get_type(fortran_ptr_, &value);
  return value;
}
void ExpressionTreeProxy::set_type(int value) {
  expression_tree_struct_set_type(fortran_ptr_, value);
}
double ExpressionTreeProxy::value() const {
  double value;
  expression_tree_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ExpressionTreeProxy::set_value(double value) {
  expression_tree_struct_set_value(fortran_ptr_, value);
}
ExpressionTreeProxyArray1D ExpressionTreeProxy::node() const {
  return ProxyHelpers::get_type_array_1d<ExpressionTreeProxyArray1D>(
      fortran_ptr_, expression_tree_struct_get_node_info);
}
FCharArray1D NametableProxy::name() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, nametable_struct_get_name_info);
}
FArray1D<int> NametableProxy::index() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, nametable_struct_get_index_info);
}
int NametableProxy::n_min() const {
  int value;
  nametable_struct_get_n_min(fortran_ptr_, &value);
  return value;
}
void NametableProxy::set_n_min(int value) {
  nametable_struct_set_n_min(fortran_ptr_, value);
}
int NametableProxy::n_max() const {
  int value;
  nametable_struct_get_n_max(fortran_ptr_, &value);
  return value;
}
void NametableProxy::set_n_max(int value) {
  nametable_struct_set_n_max(fortran_ptr_, value);
}
FArray1D<double> TaoSpinDnDpzProxy::vec() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_spin_dn_dpz_struct_get_vec_info);
}
FArray2D<double> TaoSpinDnDpzProxy::partial() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_spin_dn_dpz_struct_get_partial_info);
}
FArray2D<double> TaoSpinDnDpzProxy::partial2() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_spin_dn_dpz_struct_get_partial2_info);
}
std::string ResonanceHProxy::id() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, resonance_h_struct_get_id_info);
  return std::string(arr.data(), arr.size());
}
void ResonanceHProxy::set_id(const std::string& value) {
  resonance_h_struct_set_id(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::complex<double> ResonanceHProxy::c_val() const {
  std::complex<double> c_value;
  resonance_h_struct_get_c_val(fortran_ptr_, &c_value);
  return c_value;
}
void ResonanceHProxy::set_c_val(std::complex<double> value) {
  resonance_h_struct_set_c_val(fortran_ptr_, value);
}
FArray2D<double> SpinOrbitMap1Proxy::orb_mat() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, spin_orbit_map1_struct_get_orb_mat_info);
}
FArray1D<double> SpinOrbitMap1Proxy::vec0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, spin_orbit_map1_struct_get_vec0_info);
}
FArray2D<double> SpinOrbitMap1Proxy::spin_q() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, spin_orbit_map1_struct_get_spin_q_info);
}
FArray1D<double> SpinAxisProxy::l() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, spin_axis_struct_get_l_info);
}
FArray1D<double> SpinAxisProxy::n0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, spin_axis_struct_get_n0_info);
}
FArray1D<double> SpinAxisProxy::m() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, spin_axis_struct_get_m_info);
}
std::optional<EleProxy> PtcNormalFormProxy::ele_origin() const {
  void* ptr;
  ptc_normal_form_struct_get_ele_origin(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleProxy(ptr);
}
void PtcNormalFormProxy::set_ele_origin(const EleProxy& src) {
  ptc_normal_form_struct_set_ele_origin(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> PtcNormalFormProxy::orb0() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, ptc_normal_form_struct_get_orb0_info);
}
bool PtcNormalFormProxy::valid_map() const {
  bool value;
  ptc_normal_form_struct_get_valid_map(fortran_ptr_, &value);
  return value;
}
void PtcNormalFormProxy::set_valid_map(bool value) {
  ptc_normal_form_struct_set_valid_map(fortran_ptr_, value);
}
std::optional<EleProxy> BmadNormalFormProxy::ele_origin() const {
  void* ptr;
  bmad_normal_form_struct_get_ele_origin(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleProxy(ptr);
}
void BmadNormalFormProxy::set_ele_origin(const EleProxy& src) {
  bmad_normal_form_struct_set_ele_origin(fortran_ptr_, src.get_fortran_ptr());
}
TaylorProxyArray1D BmadNormalFormProxy::M() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_M_info);
}
TaylorProxyArray1D BmadNormalFormProxy::A() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_A_info);
}
TaylorProxyArray1D BmadNormalFormProxy::A_inv() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_A_inv_info);
}
TaylorProxyArray1D BmadNormalFormProxy::dhdj() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_dhdj_info);
}
ComplexTaylorProxyArray1D BmadNormalFormProxy::F() const {
  return ProxyHelpers::get_type_array_1d<ComplexTaylorProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_F_info);
}
ComplexTaylorProxyArray1D BmadNormalFormProxy::L() const {
  return ProxyHelpers::get_type_array_1d<ComplexTaylorProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_L_info);
}
ResonanceHProxyArray1D BmadNormalFormProxy::h() const {
  return ProxyHelpers::get_type_array_1d<ResonanceHProxyArray1D>(
      fortran_ptr_, bmad_normal_form_struct_get_h_info);
}
BunchParamsProxyArray1D BunchTrackProxy::pt() const {
  return ProxyHelpers::get_type_array_1d<BunchParamsProxyArray1D>(
      fortran_ptr_, bunch_track_struct_get_pt_info);
}
double BunchTrackProxy::ds_save() const {
  double value;
  bunch_track_struct_get_ds_save(fortran_ptr_, &value);
  return value;
}
void BunchTrackProxy::set_ds_save(double value) {
  bunch_track_struct_set_ds_save(fortran_ptr_, value);
}
int BunchTrackProxy::n_pt() const {
  int value;
  bunch_track_struct_get_n_pt(fortran_ptr_, &value);
  return value;
}
void BunchTrackProxy::set_n_pt(int value) {
  bunch_track_struct_set_n_pt(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h11001() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h11001(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h11001(std::complex<double> value) {
  summation_rdt_struct_set_h11001(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h00111() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00111(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h00111(std::complex<double> value) {
  summation_rdt_struct_set_h00111(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h20001() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20001(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h20001(std::complex<double> value) {
  summation_rdt_struct_set_h20001(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h00201() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00201(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h00201(std::complex<double> value) {
  summation_rdt_struct_set_h00201(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h10002() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10002(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h10002(std::complex<double> value) {
  summation_rdt_struct_set_h10002(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h21000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h21000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h21000(std::complex<double> value) {
  summation_rdt_struct_set_h21000(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h30000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h30000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h30000(std::complex<double> value) {
  summation_rdt_struct_set_h30000(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h10110() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10110(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h10110(std::complex<double> value) {
  summation_rdt_struct_set_h10110(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h10020() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10020(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h10020(std::complex<double> value) {
  summation_rdt_struct_set_h10020(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h10200() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10200(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h10200(std::complex<double> value) {
  summation_rdt_struct_set_h10200(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h31000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h31000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h31000(std::complex<double> value) {
  summation_rdt_struct_set_h31000(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h40000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h40000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h40000(std::complex<double> value) {
  summation_rdt_struct_set_h40000(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h20110() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20110(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h20110(std::complex<double> value) {
  summation_rdt_struct_set_h20110(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h11200() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h11200(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h11200(std::complex<double> value) {
  summation_rdt_struct_set_h11200(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h20020() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20020(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h20020(std::complex<double> value) {
  summation_rdt_struct_set_h20020(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h20200() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20200(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h20200(std::complex<double> value) {
  summation_rdt_struct_set_h20200(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h00310() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00310(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h00310(std::complex<double> value) {
  summation_rdt_struct_set_h00310(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h00400() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00400(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h00400(std::complex<double> value) {
  summation_rdt_struct_set_h00400(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h22000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h22000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h22000(std::complex<double> value) {
  summation_rdt_struct_set_h22000(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h00220() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00220(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h00220(std::complex<double> value) {
  summation_rdt_struct_set_h00220(fortran_ptr_, value);
}
std::complex<double> SummationRdtProxy::h11110() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h11110(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtProxy::set_h11110(std::complex<double> value) {
  summation_rdt_struct_set_h11110(fortran_ptr_, value);
}
std::string TaoEleShapeProxy::ele_id() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_ele_shape_struct_get_ele_id_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeProxy::set_ele_id(const std::string& value) {
  tao_ele_shape_struct_set_ele_id(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoEleShapeProxy::shape() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_ele_shape_struct_get_shape_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeProxy::set_shape(const std::string& value) {
  tao_ele_shape_struct_set_shape(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoEleShapeProxy::color() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_ele_shape_struct_get_color_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeProxy::set_color(const std::string& value) {
  tao_ele_shape_struct_set_color(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoEleShapeProxy::size() const {
  double value;
  tao_ele_shape_struct_get_size(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeProxy::set_size(double value) {
  tao_ele_shape_struct_set_size(fortran_ptr_, value);
}
std::string TaoEleShapeProxy::label() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_ele_shape_struct_get_label_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeProxy::set_label(const std::string& value) {
  tao_ele_shape_struct_set_label(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoEleShapeProxy::draw() const {
  bool value;
  tao_ele_shape_struct_get_draw(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeProxy::set_draw(bool value) {
  tao_ele_shape_struct_set_draw(fortran_ptr_, value);
}
bool TaoEleShapeProxy::multi() const {
  bool value;
  tao_ele_shape_struct_get_multi(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeProxy::set_multi(bool value) {
  tao_ele_shape_struct_set_multi(fortran_ptr_, value);
}
int TaoEleShapeProxy::line_width() const {
  int value;
  tao_ele_shape_struct_get_line_width(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeProxy::set_line_width(int value) {
  tao_ele_shape_struct_set_line_width(fortran_ptr_, value);
}
double TaoEleShapeProxy::offset() const {
  double value;
  tao_ele_shape_struct_get_offset(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeProxy::set_offset(double value) {
  tao_ele_shape_struct_set_offset(fortran_ptr_, value);
}
int TaoEleShapeProxy::ix_key() const {
  int value;
  tao_ele_shape_struct_get_ix_key(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeProxy::set_ix_key(int value) {
  tao_ele_shape_struct_set_ix_key(fortran_ptr_, value);
}
std::string TaoEleShapeProxy::name_ele() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_ele_shape_struct_get_name_ele_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeProxy::set_name_ele(const std::string& value) {
  tao_ele_shape_struct_set_name_ele(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoElePointerProxyArray1D TaoEleShapeProxy::uni() const {
  return ProxyHelpers::get_type_array_1d<TaoElePointerProxyArray1D>(
      fortran_ptr_, tao_ele_shape_struct_get_uni_info);
}
ElePointerProxyArray1D TaoElePointerProxy::eles() const {
  return ProxyHelpers::get_type_array_1d<ElePointerProxyArray1D>(
      fortran_ptr_, tao_ele_pointer_struct_get_eles_info);
}
int TaoElePointerProxy::n_loc() const {
  int value;
  tao_ele_pointer_struct_get_n_loc(fortran_ptr_, &value);
  return value;
}
void TaoElePointerProxy::set_n_loc(int value) {
  tao_ele_pointer_struct_set_n_loc(fortran_ptr_, value);
}
std::string TaoCurveProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_name(const std::string& value) {
  tao_curve_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::data_source() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_data_source_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_data_source(const std::string& value) {
  tao_curve_struct_set_data_source(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::data_index() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_data_index_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_data_index(const std::string& value) {
  tao_curve_struct_set_data_index(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::data_type_x() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_data_type_x_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_data_type_x(const std::string& value) {
  tao_curve_struct_set_data_type_x(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::data_type() const {
  return ProxyHelpers::get_string(
      fortran_ptr_, tao_curve_struct_get_data_type_info);
}
void TaoCurveProxy::set_data_type(const std::string& value) {
  tao_curve_struct_set_data_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::ele_ref_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_ele_ref_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_ele_ref_name(const std::string& value) {
  tao_curve_struct_set_ele_ref_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::legend_text() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_legend_text_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_legend_text(const std::string& value) {
  tao_curve_struct_set_legend_text(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::message_text() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_message_text_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_message_text(const std::string& value) {
  tao_curve_struct_set_message_text(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::component() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_component_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_component(const std::string& value) {
  tao_curve_struct_set_component(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveProxy::why_invalid() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_struct_get_why_invalid_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveProxy::set_why_invalid(const std::string& value) {
  tao_curve_struct_set_why_invalid(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::optional<TaoGraphProxy> TaoCurveProxy::g() const {
  void* ptr;
  tao_curve_struct_get_g(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoGraphProxy(ptr);
}
void TaoCurveProxy::set_g(const TaoGraphProxy& src) {
  tao_curve_struct_set_g(fortran_ptr_, src.get_fortran_ptr());
}
TaoHistogramProxy TaoCurveProxy::hist() const {
  void* ptr;
  tao_curve_struct_get_hist(fortran_ptr_, &ptr);
  return TaoHistogramProxy(ptr);
}
void TaoCurveProxy::set_hist(const TaoHistogramProxy& src) {
  tao_curve_struct_set_hist(fortran_ptr_, src.get_fortran_ptr());
}
TaoCurveColorProxy TaoCurveProxy::z_color() const {
  void* ptr;
  tao_curve_struct_get_z_color(fortran_ptr_, &ptr);
  return TaoCurveColorProxy(ptr);
}
void TaoCurveProxy::set_z_color(const TaoCurveColorProxy& src) {
  tao_curve_struct_set_z_color(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TaoCurveProxy::x_line() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_x_line_info);
}
FArray1D<double> TaoCurveProxy::y_line() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_y_line_info);
}
FArray1D<double> TaoCurveProxy::y2_line() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_y2_line_info);
}
FArray1D<int> TaoCurveProxy::ix_line() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, tao_curve_struct_get_ix_line_info);
}
FArray1D<double> TaoCurveProxy::x_symb() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_x_symb_info);
}
FArray1D<double> TaoCurveProxy::y_symb() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_y_symb_info);
}
FArray1D<double> TaoCurveProxy::z_symb() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_z_symb_info);
}
FArray1D<double> TaoCurveProxy::err_symb() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_err_symb_info);
}
FArray1D<double> TaoCurveProxy::symb_size() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_curve_struct_get_symb_size_info);
}
FArray1D<int> TaoCurveProxy::ix_symb() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, tao_curve_struct_get_ix_symb_info);
}
double TaoCurveProxy::y_axis_scale_factor() const {
  double value;
  tao_curve_struct_get_y_axis_scale_factor(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_y_axis_scale_factor(double value) {
  tao_curve_struct_set_y_axis_scale_factor(fortran_ptr_, value);
}
QpLineProxy TaoCurveProxy::line() const {
  void* ptr;
  tao_curve_struct_get_line(fortran_ptr_, &ptr);
  return QpLineProxy(ptr);
}
void TaoCurveProxy::set_line(const QpLineProxy& src) {
  tao_curve_struct_set_line(fortran_ptr_, src.get_fortran_ptr());
}
QpSymbolProxy TaoCurveProxy::symbol() const {
  void* ptr;
  tao_curve_struct_get_symbol(fortran_ptr_, &ptr);
  return QpSymbolProxy(ptr);
}
void TaoCurveProxy::set_symbol(const QpSymbolProxy& src) {
  tao_curve_struct_set_symbol(fortran_ptr_, src.get_fortran_ptr());
}
TaoCurveOrbitProxy TaoCurveProxy::orbit() const {
  void* ptr;
  tao_curve_struct_get_orbit(fortran_ptr_, &ptr);
  return TaoCurveOrbitProxy(ptr);
}
void TaoCurveProxy::set_orbit(const TaoCurveOrbitProxy& src) {
  tao_curve_struct_set_orbit(fortran_ptr_, src.get_fortran_ptr());
}
int TaoCurveProxy::ix_universe() const {
  int value;
  tao_curve_struct_get_ix_universe(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_ix_universe(int value) {
  tao_curve_struct_set_ix_universe(fortran_ptr_, value);
}
int TaoCurveProxy::symbol_every() const {
  int value;
  tao_curve_struct_get_symbol_every(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_symbol_every(int value) {
  tao_curve_struct_set_symbol_every(fortran_ptr_, value);
}
int TaoCurveProxy::ix_branch() const {
  int value;
  tao_curve_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_ix_branch(int value) {
  tao_curve_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoCurveProxy::ix_bunch() const {
  int value;
  tao_curve_struct_get_ix_bunch(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_ix_bunch(int value) {
  tao_curve_struct_set_ix_bunch(fortran_ptr_, value);
}
int TaoCurveProxy::n_turn() const {
  int value;
  tao_curve_struct_get_n_turn(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_n_turn(int value) {
  tao_curve_struct_set_n_turn(fortran_ptr_, value);
}
bool TaoCurveProxy::use_y2() const {
  bool value;
  tao_curve_struct_get_use_y2(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_use_y2(bool value) {
  tao_curve_struct_set_use_y2(fortran_ptr_, value);
}
bool TaoCurveProxy::draw_line() const {
  bool value;
  tao_curve_struct_get_draw_line(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_draw_line(bool value) {
  tao_curve_struct_set_draw_line(fortran_ptr_, value);
}
bool TaoCurveProxy::draw_symbols() const {
  bool value;
  tao_curve_struct_get_draw_symbols(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_draw_symbols(bool value) {
  tao_curve_struct_set_draw_symbols(fortran_ptr_, value);
}
bool TaoCurveProxy::draw_symbol_index() const {
  bool value;
  tao_curve_struct_get_draw_symbol_index(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_draw_symbol_index(bool value) {
  tao_curve_struct_set_draw_symbol_index(fortran_ptr_, value);
}
bool TaoCurveProxy::draw_error_bars() const {
  bool value;
  tao_curve_struct_get_draw_error_bars(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_draw_error_bars(bool value) {
  tao_curve_struct_set_draw_error_bars(fortran_ptr_, value);
}
bool TaoCurveProxy::smooth_line_calc() const {
  bool value;
  tao_curve_struct_get_smooth_line_calc(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_smooth_line_calc(bool value) {
  tao_curve_struct_set_smooth_line_calc(fortran_ptr_, value);
}
bool TaoCurveProxy::valid() const {
  bool value;
  tao_curve_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoCurveProxy::set_valid(bool value) {
  tao_curve_struct_set_valid(fortran_ptr_, value);
}
std::string TaoCurveColorProxy::data_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_curve_color_struct_get_data_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveColorProxy::set_data_type(const std::string& value) {
  tao_curve_color_struct_set_data_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoCurveColorProxy::is_on() const {
  bool value;
  tao_curve_color_struct_get_is_on(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorProxy::set_is_on(bool value) {
  tao_curve_color_struct_set_is_on(fortran_ptr_, value);
}
double TaoCurveColorProxy::min() const {
  double value;
  tao_curve_color_struct_get_min(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorProxy::set_min(double value) {
  tao_curve_color_struct_set_min(fortran_ptr_, value);
}
double TaoCurveColorProxy::max() const {
  double value;
  tao_curve_color_struct_get_max(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorProxy::set_max(double value) {
  tao_curve_color_struct_set_max(fortran_ptr_, value);
}
bool TaoCurveColorProxy::autoscale() const {
  bool value;
  tao_curve_color_struct_get_autoscale(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorProxy::set_autoscale(bool value) {
  tao_curve_color_struct_set_autoscale(fortran_ptr_, value);
}
double TaoCurveOrbitProxy::x() const {
  double value;
  tao_curve_orbit_struct_get_x(fortran_ptr_, &value);
  return value;
}
void TaoCurveOrbitProxy::set_x(double value) {
  tao_curve_orbit_struct_set_x(fortran_ptr_, value);
}
double TaoCurveOrbitProxy::y() const {
  double value;
  tao_curve_orbit_struct_get_y(fortran_ptr_, &value);
  return value;
}
void TaoCurveOrbitProxy::set_y(double value) {
  tao_curve_orbit_struct_set_y(fortran_ptr_, value);
}
double TaoCurveOrbitProxy::t() const {
  double value;
  tao_curve_orbit_struct_get_t(fortran_ptr_, &value);
  return value;
}
void TaoCurveOrbitProxy::set_t(double value) {
  tao_curve_orbit_struct_set_t(fortran_ptr_, value);
}
bool TaoHistogramProxy::density_normalized() const {
  bool value;
  tao_histogram_struct_get_density_normalized(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_density_normalized(bool value) {
  tao_histogram_struct_set_density_normalized(fortran_ptr_, value);
}
bool TaoHistogramProxy::weight_by_charge() const {
  bool value;
  tao_histogram_struct_get_weight_by_charge(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_weight_by_charge(bool value) {
  tao_histogram_struct_set_weight_by_charge(fortran_ptr_, value);
}
double TaoHistogramProxy::minimum() const {
  double value;
  tao_histogram_struct_get_minimum(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_minimum(double value) {
  tao_histogram_struct_set_minimum(fortran_ptr_, value);
}
double TaoHistogramProxy::maximum() const {
  double value;
  tao_histogram_struct_get_maximum(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_maximum(double value) {
  tao_histogram_struct_set_maximum(fortran_ptr_, value);
}
double TaoHistogramProxy::width() const {
  double value;
  tao_histogram_struct_get_width(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_width(double value) {
  tao_histogram_struct_set_width(fortran_ptr_, value);
}
double TaoHistogramProxy::center() const {
  double value;
  tao_histogram_struct_get_center(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_center(double value) {
  tao_histogram_struct_set_center(fortran_ptr_, value);
}
int TaoHistogramProxy::number() const {
  int value;
  tao_histogram_struct_get_number(fortran_ptr_, &value);
  return value;
}
void TaoHistogramProxy::set_number(int value) {
  tao_histogram_struct_set_number(fortran_ptr_, value);
}
int LatEleOrder1Proxy::ix_branch() const {
  int value;
  lat_ele_order1_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void LatEleOrder1Proxy::set_ix_branch(int value) {
  lat_ele_order1_struct_set_ix_branch(fortran_ptr_, value);
}
int LatEleOrder1Proxy::ix_order() const {
  int value;
  lat_ele_order1_struct_get_ix_order(fortran_ptr_, &value);
  return value;
}
void LatEleOrder1Proxy::set_ix_order(int value) {
  lat_ele_order1_struct_set_ix_order(fortran_ptr_, value);
}
LatEleOrder1ProxyArray1D LatEleOrderArrayProxy::ele() const {
  return ProxyHelpers::get_type_array_1d<LatEleOrder1ProxyArray1D>(
      fortran_ptr_, lat_ele_order_array_struct_get_ele_info);
}
FArray2D<double> TaoLatSigmaProxy::mat() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_lat_sigma_struct_get_mat_info);
}
TaoSpinDnDpzProxy TaoSpinEleProxy::dn_dpz() const {
  void* ptr;
  tao_spin_ele_struct_get_dn_dpz(fortran_ptr_, &ptr);
  return TaoSpinDnDpzProxy(ptr);
}
void TaoSpinEleProxy::set_dn_dpz(const TaoSpinDnDpzProxy& src) {
  tao_spin_ele_struct_set_dn_dpz(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TaoSpinEleProxy::orb_eigen_val() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_spin_ele_struct_get_orb_eigen_val_info);
}
FArray2D<double> TaoSpinEleProxy::orb_eigen_vec() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_spin_ele_struct_get_orb_eigen_vec_info);
}
FArray2D<double> TaoSpinEleProxy::spin_eigen_vec() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_spin_ele_struct_get_spin_eigen_vec_info);
}
bool TaoSpinEleProxy::valid() const {
  bool value;
  tao_spin_ele_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoSpinEleProxy::set_valid(bool value) {
  tao_spin_ele_struct_set_valid(fortran_ptr_, value);
}
EleProxy TaoPlotCacheProxy::ele_to_s() const {
  void* ptr;
  tao_plot_cache_struct_get_ele_to_s(fortran_ptr_, &ptr);
  return EleProxy(ptr);
}
void TaoPlotCacheProxy::set_ele_to_s(const EleProxy& src) {
  tao_plot_cache_struct_set_ele_to_s(fortran_ptr_, src.get_fortran_ptr());
}
CoordProxy TaoPlotCacheProxy::orbit() const {
  void* ptr;
  tao_plot_cache_struct_get_orbit(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void TaoPlotCacheProxy::set_orbit(const CoordProxy& src) {
  tao_plot_cache_struct_set_orbit(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoPlotCacheProxy::err() const {
  bool value;
  tao_plot_cache_struct_get_err(fortran_ptr_, &value);
  return value;
}
void TaoPlotCacheProxy::set_err(bool value) {
  tao_plot_cache_struct_set_err(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::tune() const {
  double value;
  tao_spin_polarization_struct_get_tune(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_tune(double value) {
  tao_spin_polarization_struct_set_tune(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::pol_limit_st() const {
  double value;
  tao_spin_polarization_struct_get_pol_limit_st(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_pol_limit_st(double value) {
  tao_spin_polarization_struct_set_pol_limit_st(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::pol_limit_dk() const {
  double value;
  tao_spin_polarization_struct_get_pol_limit_dk(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_pol_limit_dk(double value) {
  tao_spin_polarization_struct_set_pol_limit_dk(fortran_ptr_, value);
}
FArray1D<double> TaoSpinPolarizationProxy::pol_limit_dk_partial() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_spin_polarization_struct_get_pol_limit_dk_partial_info);
}
FArray1D<double> TaoSpinPolarizationProxy::pol_limit_dk_partial2() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      tao_spin_polarization_struct_get_pol_limit_dk_partial2_info);
}
double TaoSpinPolarizationProxy::pol_rate_bks() const {
  double value;
  tao_spin_polarization_struct_get_pol_rate_bks(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_pol_rate_bks(double value) {
  tao_spin_polarization_struct_set_pol_rate_bks(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::depol_rate() const {
  double value;
  tao_spin_polarization_struct_get_depol_rate(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_depol_rate(double value) {
  tao_spin_polarization_struct_set_depol_rate(fortran_ptr_, value);
}
FArray1D<double> TaoSpinPolarizationProxy::depol_rate_partial() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_spin_polarization_struct_get_depol_rate_partial_info);
}
FArray1D<double> TaoSpinPolarizationProxy::depol_rate_partial2() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_spin_polarization_struct_get_depol_rate_partial2_info);
}
double TaoSpinPolarizationProxy::integral_bn() const {
  double value;
  tao_spin_polarization_struct_get_integral_bn(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_integral_bn(double value) {
  tao_spin_polarization_struct_set_integral_bn(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::integral_bdn() const {
  double value;
  tao_spin_polarization_struct_get_integral_bdn(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_integral_bdn(double value) {
  tao_spin_polarization_struct_set_integral_bdn(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::integral_1ns() const {
  double value;
  tao_spin_polarization_struct_get_integral_1ns(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_integral_1ns(double value) {
  tao_spin_polarization_struct_set_integral_1ns(fortran_ptr_, value);
}
double TaoSpinPolarizationProxy::integral_dn2() const {
  double value;
  tao_spin_polarization_struct_get_integral_dn2(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_integral_dn2(double value) {
  tao_spin_polarization_struct_set_integral_dn2(fortran_ptr_, value);
}
bool TaoSpinPolarizationProxy::valid() const {
  bool value;
  tao_spin_polarization_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationProxy::set_valid(bool value) {
  tao_spin_polarization_struct_set_valid(fortran_ptr_, value);
}
SpinOrbitMap1Proxy TaoSpinPolarizationProxy::q_1turn() const {
  void* ptr;
  tao_spin_polarization_struct_get_q_1turn(fortran_ptr_, &ptr);
  return SpinOrbitMap1Proxy(ptr);
}
void TaoSpinPolarizationProxy::set_q_1turn(const SpinOrbitMap1Proxy& src) {
  tao_spin_polarization_struct_set_q_1turn(fortran_ptr_, src.get_fortran_ptr());
}
SpinOrbitMap1ProxyArray1D TaoSpinPolarizationProxy::q_ele() const {
  return ProxyHelpers::get_type_array_1d<SpinOrbitMap1ProxyArray1D>(
      fortran_ptr_, tao_spin_polarization_struct_get_q_ele_info);
}
std::optional<TaoLatticeProxy> TaoLatticeBranchProxy::tao_lat() const {
  void* ptr;
  tao_lattice_branch_struct_get_tao_lat(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeProxy(ptr);
}
void TaoLatticeBranchProxy::set_tao_lat(const TaoLatticeProxy& src) {
  tao_lattice_branch_struct_set_tao_lat(fortran_ptr_, src.get_fortran_ptr());
}
TaoLatSigmaProxyArray1D TaoLatticeBranchProxy::lat_sigma() const {
  return ProxyHelpers::get_type_array_1d<TaoLatSigmaProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_lat_sigma_info);
}
TaoSpinEleProxyArray1D TaoLatticeBranchProxy::spin_ele() const {
  return ProxyHelpers::get_type_array_1d<TaoSpinEleProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_spin_ele_info);
}
BunchParamsProxyArray1D TaoLatticeBranchProxy::bunch_params() const {
  return ProxyHelpers::get_type_array_1d<BunchParamsProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_bunch_params_info);
}
BunchTrackProxyArray1D TaoLatticeBranchProxy::bunch_params_comb() const {
  return ProxyHelpers::get_type_array_1d<BunchTrackProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_bunch_params_comb_info);
}
CoordProxyArray1D TaoLatticeBranchProxy::orbit() const {
  return ProxyHelpers::get_type_array_1d<CoordProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_orbit_info);
}
TaoPlotCacheProxyArray1D TaoLatticeBranchProxy::plot_cache() const {
  return ProxyHelpers::get_type_array_1d<TaoPlotCacheProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_plot_cache_info);
}
TaoSpinPolarizationProxy TaoLatticeBranchProxy::spin() const {
  void* ptr;
  tao_lattice_branch_struct_get_spin(fortran_ptr_, &ptr);
  return TaoSpinPolarizationProxy(ptr);
}
void TaoLatticeBranchProxy::set_spin(const TaoSpinPolarizationProxy& src) {
  tao_lattice_branch_struct_set_spin(fortran_ptr_, src.get_fortran_ptr());
}
SummationRdtProxy TaoLatticeBranchProxy::srdt() const {
  void* ptr;
  tao_lattice_branch_struct_get_srdt(fortran_ptr_, &ptr);
  return SummationRdtProxy(ptr);
}
void TaoLatticeBranchProxy::set_srdt(const SummationRdtProxy& src) {
  tao_lattice_branch_struct_set_srdt(fortran_ptr_, src.get_fortran_ptr());
}
CoordProxy TaoLatticeBranchProxy::orb0() const {
  void* ptr;
  tao_lattice_branch_struct_get_orb0(fortran_ptr_, &ptr);
  return CoordProxy(ptr);
}
void TaoLatticeBranchProxy::set_orb0(const CoordProxy& src) {
  tao_lattice_branch_struct_set_orb0(fortran_ptr_, src.get_fortran_ptr());
}
NormalModesProxy TaoLatticeBranchProxy::modes_ri() const {
  void* ptr;
  tao_lattice_branch_struct_get_modes_ri(fortran_ptr_, &ptr);
  return NormalModesProxy(ptr);
}
void TaoLatticeBranchProxy::set_modes_ri(const NormalModesProxy& src) {
  tao_lattice_branch_struct_set_modes_ri(fortran_ptr_, src.get_fortran_ptr());
}
NormalModesProxy TaoLatticeBranchProxy::modes_6d() const {
  void* ptr;
  tao_lattice_branch_struct_get_modes_6d(fortran_ptr_, &ptr);
  return NormalModesProxy(ptr);
}
void TaoLatticeBranchProxy::set_modes_6d(const NormalModesProxy& src) {
  tao_lattice_branch_struct_set_modes_6d(fortran_ptr_, src.get_fortran_ptr());
}
PtcNormalFormProxy TaoLatticeBranchProxy::ptc_normal_form() const {
  void* ptr;
  tao_lattice_branch_struct_get_ptc_normal_form(fortran_ptr_, &ptr);
  return PtcNormalFormProxy(ptr);
}
void TaoLatticeBranchProxy::set_ptc_normal_form(const PtcNormalFormProxy& src) {
  tao_lattice_branch_struct_set_ptc_normal_form(
      fortran_ptr_, src.get_fortran_ptr());
}
BmadNormalFormProxy TaoLatticeBranchProxy::bmad_normal_form() const {
  void* ptr;
  tao_lattice_branch_struct_get_bmad_normal_form(fortran_ptr_, &ptr);
  return BmadNormalFormProxy(ptr);
}
void TaoLatticeBranchProxy::set_bmad_normal_form(
    const BmadNormalFormProxy& src) {
  tao_lattice_branch_struct_set_bmad_normal_form(
      fortran_ptr_, src.get_fortran_ptr());
}
CoordProxyArray1D TaoLatticeBranchProxy::high_E_orb() const {
  return ProxyHelpers::get_type_array_1d<CoordProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_high_E_orb_info);
}
CoordProxyArray1D TaoLatticeBranchProxy::low_E_orb() const {
  return ProxyHelpers::get_type_array_1d<CoordProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_low_E_orb_info);
}
TaylorProxyArray1D TaoLatticeBranchProxy::taylor_save() const {
  return ProxyHelpers::get_type_array_1d<TaylorProxyArray1D>(
      fortran_ptr_, tao_lattice_branch_struct_get_taylor_save_info);
}
double TaoLatticeBranchProxy::cache_x_min() const {
  double value;
  tao_lattice_branch_struct_get_cache_x_min(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_cache_x_min(double value) {
  tao_lattice_branch_struct_set_cache_x_min(fortran_ptr_, value);
}
double TaoLatticeBranchProxy::cache_x_max() const {
  double value;
  tao_lattice_branch_struct_get_cache_x_max(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_cache_x_max(double value) {
  tao_lattice_branch_struct_set_cache_x_max(fortran_ptr_, value);
}
double TaoLatticeBranchProxy::comb_ds_save() const {
  double value;
  tao_lattice_branch_struct_get_comb_ds_save(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_comb_ds_save(double value) {
  tao_lattice_branch_struct_set_comb_ds_save(fortran_ptr_, value);
}
int TaoLatticeBranchProxy::ix_ref_taylor() const {
  int value;
  tao_lattice_branch_struct_get_ix_ref_taylor(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_ix_ref_taylor(int value) {
  tao_lattice_branch_struct_set_ix_ref_taylor(fortran_ptr_, value);
}
int TaoLatticeBranchProxy::ix_ele_taylor() const {
  int value;
  tao_lattice_branch_struct_get_ix_ele_taylor(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_ix_ele_taylor(int value) {
  tao_lattice_branch_struct_set_ix_ele_taylor(fortran_ptr_, value);
}
int TaoLatticeBranchProxy::track_state() const {
  int value;
  tao_lattice_branch_struct_get_track_state(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_track_state(int value) {
  tao_lattice_branch_struct_set_track_state(fortran_ptr_, value);
}
int TaoLatticeBranchProxy::cache_n_pts() const {
  int value;
  tao_lattice_branch_struct_get_cache_n_pts(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_cache_n_pts(int value) {
  tao_lattice_branch_struct_set_cache_n_pts(fortran_ptr_, value);
}
int TaoLatticeBranchProxy::ix_rad_int_cache() const {
  int value;
  tao_lattice_branch_struct_get_ix_rad_int_cache(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_ix_rad_int_cache(int value) {
  tao_lattice_branch_struct_set_ix_rad_int_cache(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::has_open_match_element() const {
  bool value;
  tao_lattice_branch_struct_get_has_open_match_element(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_has_open_match_element(bool value) {
  tao_lattice_branch_struct_set_has_open_match_element(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::plot_cache_valid() const {
  bool value;
  tao_lattice_branch_struct_get_plot_cache_valid(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_plot_cache_valid(bool value) {
  tao_lattice_branch_struct_set_plot_cache_valid(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::spin_map_valid() const {
  bool value;
  tao_lattice_branch_struct_get_spin_map_valid(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_spin_map_valid(bool value) {
  tao_lattice_branch_struct_set_spin_map_valid(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::twiss_valid() const {
  bool value;
  tao_lattice_branch_struct_get_twiss_valid(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_twiss_valid(bool value) {
  tao_lattice_branch_struct_set_twiss_valid(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::mode_flip_here() const {
  bool value;
  tao_lattice_branch_struct_get_mode_flip_here(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_mode_flip_here(bool value) {
  tao_lattice_branch_struct_set_mode_flip_here(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::chrom_calc_ok() const {
  bool value;
  tao_lattice_branch_struct_get_chrom_calc_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_chrom_calc_ok(bool value) {
  tao_lattice_branch_struct_set_chrom_calc_ok(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::rad_int_calc_ok() const {
  bool value;
  tao_lattice_branch_struct_get_rad_int_calc_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_rad_int_calc_ok(bool value) {
  tao_lattice_branch_struct_set_rad_int_calc_ok(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::emit_6d_calc_ok() const {
  bool value;
  tao_lattice_branch_struct_get_emit_6d_calc_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_emit_6d_calc_ok(bool value) {
  tao_lattice_branch_struct_set_emit_6d_calc_ok(fortran_ptr_, value);
}
bool TaoLatticeBranchProxy::sigma_track_ok() const {
  bool value;
  tao_lattice_branch_struct_get_sigma_track_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchProxy::set_sigma_track_ok(bool value) {
  tao_lattice_branch_struct_set_sigma_track_ok(fortran_ptr_, value);
}
BeamProxy TaoModelElementProxy::beam() const {
  void* ptr;
  tao_model_element_struct_get_beam(fortran_ptr_, &ptr);
  return BeamProxy(ptr);
}
void TaoModelElementProxy::set_beam(const BeamProxy& src) {
  tao_model_element_struct_set_beam(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoModelElementProxy::save_beam_internally() const {
  bool value;
  tao_model_element_struct_get_save_beam_internally(fortran_ptr_, &value);
  return value;
}
void TaoModelElementProxy::set_save_beam_internally(bool value) {
  tao_model_element_struct_set_save_beam_internally(fortran_ptr_, value);
}
bool TaoModelElementProxy::save_beam_to_file() const {
  bool value;
  tao_model_element_struct_get_save_beam_to_file(fortran_ptr_, &value);
  return value;
}
void TaoModelElementProxy::set_save_beam_to_file(bool value) {
  tao_model_element_struct_set_save_beam_to_file(fortran_ptr_, value);
}
BeamProxy TaoBeamBranchProxy::beam_at_start() const {
  void* ptr;
  tao_beam_branch_struct_get_beam_at_start(fortran_ptr_, &ptr);
  return BeamProxy(ptr);
}
void TaoBeamBranchProxy::set_beam_at_start(const BeamProxy& src) {
  tao_beam_branch_struct_set_beam_at_start(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitProxy TaoBeamBranchProxy::beam_init() const {
  void* ptr;
  tao_beam_branch_struct_get_beam_init(fortran_ptr_, &ptr);
  return BeamInitProxy(ptr);
}
void TaoBeamBranchProxy::set_beam_init(const BeamInitProxy& src) {
  tao_beam_branch_struct_set_beam_init(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitProxy TaoBeamBranchProxy::beam_init_used() const {
  void* ptr;
  tao_beam_branch_struct_get_beam_init_used(fortran_ptr_, &ptr);
  return BeamInitProxy(ptr);
}
void TaoBeamBranchProxy::set_beam_init_used(const BeamInitProxy& src) {
  tao_beam_branch_struct_set_beam_init_used(
      fortran_ptr_, src.get_fortran_ptr());
}
bool TaoBeamBranchProxy::init_starting_distribution() const {
  bool value;
  tao_beam_branch_struct_get_init_starting_distribution(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchProxy::set_init_starting_distribution(bool value) {
  tao_beam_branch_struct_set_init_starting_distribution(fortran_ptr_, value);
}
std::string TaoBeamBranchProxy::track_start() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_beam_branch_struct_get_track_start_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamBranchProxy::set_track_start(const std::string& value) {
  tao_beam_branch_struct_set_track_start(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoBeamBranchProxy::track_end() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_beam_branch_struct_get_track_end_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamBranchProxy::set_track_end(const std::string& value) {
  tao_beam_branch_struct_set_track_end(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoBeamBranchProxy::ix_branch() const {
  int value;
  tao_beam_branch_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchProxy::set_ix_branch(int value) {
  tao_beam_branch_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoBeamBranchProxy::ix_track_start() const {
  int value;
  tao_beam_branch_struct_get_ix_track_start(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchProxy::set_ix_track_start(int value) {
  tao_beam_branch_struct_set_ix_track_start(fortran_ptr_, value);
}
int TaoBeamBranchProxy::ix_track_end() const {
  int value;
  tao_beam_branch_struct_get_ix_track_end(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchProxy::set_ix_track_end(int value) {
  tao_beam_branch_struct_set_ix_track_end(fortran_ptr_, value);
}
std::string TaoD1DataProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_d1_data_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD1DataProxy::set_name(const std::string& value) {
  tao_d1_data_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::optional<TaoD2DataProxy> TaoD1DataProxy::d2() const {
  void* ptr;
  tao_d1_data_struct_get_d2(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoD2DataProxy(ptr);
}
void TaoD1DataProxy::set_d2(const TaoD2DataProxy& src) {
  tao_d1_data_struct_set_d2(fortran_ptr_, src.get_fortran_ptr());
}
TaoDataProxyArray1D TaoD1DataProxy::d() const {
  return ProxyHelpers::get_type_array_1d<TaoDataProxyArray1D>(
      fortran_ptr_, tao_d1_data_struct_get_d_info);
}
std::string TaoD2DataProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_d2_data_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataProxy::set_name(const std::string& value) {
  tao_d2_data_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoD2DataProxy::data_file_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_d2_data_struct_get_data_file_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataProxy::set_data_file_name(const std::string& value) {
  tao_d2_data_struct_set_data_file_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoD2DataProxy::ref_file_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_d2_data_struct_get_ref_file_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataProxy::set_ref_file_name(const std::string& value) {
  tao_d2_data_struct_set_ref_file_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoD2DataProxy::data_date() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_d2_data_struct_get_data_date_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataProxy::set_data_date(const std::string& value) {
  tao_d2_data_struct_set_data_date(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoD2DataProxy::ref_date() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_d2_data_struct_get_ref_date_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataProxy::set_ref_date(const std::string& value) {
  tao_d2_data_struct_set_ref_date(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D TaoD2DataProxy::descrip() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, tao_d2_data_struct_get_descrip_info);
}
TaoD1DataProxyArray1D TaoD2DataProxy::d1() const {
  return ProxyHelpers::get_type_array_1d<TaoD1DataProxyArray1D>(
      fortran_ptr_, tao_d2_data_struct_get_d1_info);
}
int TaoD2DataProxy::ix_universe() const {
  int value;
  tao_d2_data_struct_get_ix_universe(fortran_ptr_, &value);
  return value;
}
void TaoD2DataProxy::set_ix_universe(int value) {
  tao_d2_data_struct_set_ix_universe(fortran_ptr_, value);
}
int TaoD2DataProxy::ix_d2_data() const {
  int value;
  tao_d2_data_struct_get_ix_d2_data(fortran_ptr_, &value);
  return value;
}
void TaoD2DataProxy::set_ix_d2_data(int value) {
  tao_d2_data_struct_set_ix_d2_data(fortran_ptr_, value);
}
int TaoD2DataProxy::ix_ref() const {
  int value;
  tao_d2_data_struct_get_ix_ref(fortran_ptr_, &value);
  return value;
}
void TaoD2DataProxy::set_ix_ref(int value) {
  tao_d2_data_struct_set_ix_ref(fortran_ptr_, value);
}
bool TaoD2DataProxy::data_read_in() const {
  bool value;
  tao_d2_data_struct_get_data_read_in(fortran_ptr_, &value);
  return value;
}
void TaoD2DataProxy::set_data_read_in(bool value) {
  tao_d2_data_struct_set_data_read_in(fortran_ptr_, value);
}
bool TaoD2DataProxy::ref_read_in() const {
  bool value;
  tao_d2_data_struct_get_ref_read_in(fortran_ptr_, &value);
  return value;
}
void TaoD2DataProxy::set_ref_read_in(bool value) {
  tao_d2_data_struct_set_ref_read_in(fortran_ptr_, value);
}
std::string TaoDataVarComponentProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_var_component_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataVarComponentProxy::set_name(const std::string& value) {
  tao_data_var_component_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoDataVarComponentProxy::sign() const {
  double value;
  tao_data_var_component_struct_get_sign(fortran_ptr_, &value);
  return value;
}
void TaoDataVarComponentProxy::set_sign(double value) {
  tao_data_var_component_struct_set_sign(fortran_ptr_, value);
}
std::string TaoGraphProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_graph_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphProxy::set_name(const std::string& value) {
  tao_graph_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGraphProxy::type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_graph_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphProxy::set_type(const std::string& value) {
  tao_graph_struct_set_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGraphProxy::title() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_graph_struct_get_title_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphProxy::set_title(const std::string& value) {
  tao_graph_struct_set_title(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGraphProxy::title_suffix() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_graph_struct_get_title_suffix_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphProxy::set_title_suffix(const std::string& value) {
  tao_graph_struct_set_title_suffix(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D TaoGraphProxy::text_legend() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, tao_graph_struct_get_text_legend_info);
}
FCharArray1D TaoGraphProxy::text_legend_out() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, tao_graph_struct_get_text_legend_out_info);
}
std::string TaoGraphProxy::why_invalid() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_graph_struct_get_why_invalid_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphProxy::set_why_invalid(const std::string& value) {
  tao_graph_struct_set_why_invalid(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoCurveProxyArray1D TaoGraphProxy::curve() const {
  return ProxyHelpers::get_type_array_1d<TaoCurveProxyArray1D>(
      fortran_ptr_, tao_graph_struct_get_curve_info);
}
std::optional<TaoPlotProxy> TaoGraphProxy::p() const {
  void* ptr;
  tao_graph_struct_get_p(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoPlotProxy(ptr);
}
void TaoGraphProxy::set_p(const TaoPlotProxy& src) {
  tao_graph_struct_set_p(fortran_ptr_, src.get_fortran_ptr());
}
TaoFloorPlanProxy TaoGraphProxy::floor_plan() const {
  void* ptr;
  tao_graph_struct_get_floor_plan(fortran_ptr_, &ptr);
  return TaoFloorPlanProxy(ptr);
}
void TaoGraphProxy::set_floor_plan(const TaoFloorPlanProxy& src) {
  tao_graph_struct_set_floor_plan(fortran_ptr_, src.get_fortran_ptr());
}
QpPointProxy TaoGraphProxy::text_legend_origin() const {
  void* ptr;
  tao_graph_struct_get_text_legend_origin(fortran_ptr_, &ptr);
  return QpPointProxy(ptr);
}
void TaoGraphProxy::set_text_legend_origin(const QpPointProxy& src) {
  tao_graph_struct_set_text_legend_origin(fortran_ptr_, src.get_fortran_ptr());
}
QpPointProxy TaoGraphProxy::curve_legend_origin() const {
  void* ptr;
  tao_graph_struct_get_curve_legend_origin(fortran_ptr_, &ptr);
  return QpPointProxy(ptr);
}
void TaoGraphProxy::set_curve_legend_origin(const QpPointProxy& src) {
  tao_graph_struct_set_curve_legend_origin(fortran_ptr_, src.get_fortran_ptr());
}
QpLegendProxy TaoGraphProxy::curve_legend() const {
  void* ptr;
  tao_graph_struct_get_curve_legend(fortran_ptr_, &ptr);
  return QpLegendProxy(ptr);
}
void TaoGraphProxy::set_curve_legend(const QpLegendProxy& src) {
  tao_graph_struct_set_curve_legend(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisProxy TaoGraphProxy::x() const {
  void* ptr;
  tao_graph_struct_get_x(fortran_ptr_, &ptr);
  return QpAxisProxy(ptr);
}
void TaoGraphProxy::set_x(const QpAxisProxy& src) {
  tao_graph_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisProxy TaoGraphProxy::y() const {
  void* ptr;
  tao_graph_struct_get_y(fortran_ptr_, &ptr);
  return QpAxisProxy(ptr);
}
void TaoGraphProxy::set_y(const QpAxisProxy& src) {
  tao_graph_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisProxy TaoGraphProxy::x2() const {
  void* ptr;
  tao_graph_struct_get_x2(fortran_ptr_, &ptr);
  return QpAxisProxy(ptr);
}
void TaoGraphProxy::set_x2(const QpAxisProxy& src) {
  tao_graph_struct_set_x2(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisProxy TaoGraphProxy::y2() const {
  void* ptr;
  tao_graph_struct_get_y2(fortran_ptr_, &ptr);
  return QpAxisProxy(ptr);
}
void TaoGraphProxy::set_y2(const QpAxisProxy& src) {
  tao_graph_struct_set_y2(fortran_ptr_, src.get_fortran_ptr());
}
QpRectProxy TaoGraphProxy::margin() const {
  void* ptr;
  tao_graph_struct_get_margin(fortran_ptr_, &ptr);
  return QpRectProxy(ptr);
}
void TaoGraphProxy::set_margin(const QpRectProxy& src) {
  tao_graph_struct_set_margin(fortran_ptr_, src.get_fortran_ptr());
}
QpRectProxy TaoGraphProxy::scale_margin() const {
  void* ptr;
  tao_graph_struct_get_scale_margin(fortran_ptr_, &ptr);
  return QpRectProxy(ptr);
}
void TaoGraphProxy::set_scale_margin(const QpRectProxy& src) {
  tao_graph_struct_set_scale_margin(fortran_ptr_, src.get_fortran_ptr());
}
double TaoGraphProxy::x_axis_scale_factor() const {
  double value;
  tao_graph_struct_get_x_axis_scale_factor(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_x_axis_scale_factor(double value) {
  tao_graph_struct_set_x_axis_scale_factor(fortran_ptr_, value);
}
double TaoGraphProxy::symbol_size_scale() const {
  double value;
  tao_graph_struct_get_symbol_size_scale(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_symbol_size_scale(double value) {
  tao_graph_struct_set_symbol_size_scale(fortran_ptr_, value);
}
FArray1D<int> TaoGraphProxy::box() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, tao_graph_struct_get_box_info);
}
int TaoGraphProxy::ix_branch() const {
  int value;
  tao_graph_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_ix_branch(int value) {
  tao_graph_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoGraphProxy::ix_universe() const {
  int value;
  tao_graph_struct_get_ix_universe(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_ix_universe(int value) {
  tao_graph_struct_set_ix_universe(fortran_ptr_, value);
}
bool TaoGraphProxy::clip() const {
  bool value;
  tao_graph_struct_get_clip(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_clip(bool value) {
  tao_graph_struct_set_clip(fortran_ptr_, value);
}
bool TaoGraphProxy::y2_mirrors_y() const {
  bool value;
  tao_graph_struct_get_y2_mirrors_y(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_y2_mirrors_y(bool value) {
  tao_graph_struct_set_y2_mirrors_y(fortran_ptr_, value);
}
bool TaoGraphProxy::limited() const {
  bool value;
  tao_graph_struct_get_limited(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_limited(bool value) {
  tao_graph_struct_set_limited(fortran_ptr_, value);
}
bool TaoGraphProxy::draw_axes() const {
  bool value;
  tao_graph_struct_get_draw_axes(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_draw_axes(bool value) {
  tao_graph_struct_set_draw_axes(fortran_ptr_, value);
}
bool TaoGraphProxy::draw_curve_legend() const {
  bool value;
  tao_graph_struct_get_draw_curve_legend(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_draw_curve_legend(bool value) {
  tao_graph_struct_set_draw_curve_legend(fortran_ptr_, value);
}
bool TaoGraphProxy::draw_grid() const {
  bool value;
  tao_graph_struct_get_draw_grid(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_draw_grid(bool value) {
  tao_graph_struct_set_draw_grid(fortran_ptr_, value);
}
bool TaoGraphProxy::draw_title() const {
  bool value;
  tao_graph_struct_get_draw_title(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_draw_title(bool value) {
  tao_graph_struct_set_draw_title(fortran_ptr_, value);
}
bool TaoGraphProxy::draw_only_good_user_data_or_vars() const {
  bool value;
  tao_graph_struct_get_draw_only_good_user_data_or_vars(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_draw_only_good_user_data_or_vars(bool value) {
  tao_graph_struct_set_draw_only_good_user_data_or_vars(fortran_ptr_, value);
}
bool TaoGraphProxy::allow_wrap_around() const {
  bool value;
  tao_graph_struct_get_allow_wrap_around(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_allow_wrap_around(bool value) {
  tao_graph_struct_set_allow_wrap_around(fortran_ptr_, value);
}
bool TaoGraphProxy::is_valid() const {
  bool value;
  tao_graph_struct_get_is_valid(fortran_ptr_, &value);
  return value;
}
void TaoGraphProxy::set_is_valid(bool value) {
  tao_graph_struct_set_is_valid(fortran_ptr_, value);
}
std::string TaoPlotProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_plot_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotProxy::set_name(const std::string& value) {
  tao_plot_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoPlotProxy::description() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_plot_struct_get_description_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotProxy::set_description(const std::string& value) {
  tao_plot_struct_set_description(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoGraphProxyArray1D TaoPlotProxy::graph() const {
  return ProxyHelpers::get_type_array_1d<TaoGraphProxyArray1D>(
      fortran_ptr_, tao_plot_struct_get_graph_info);
}
std::optional<TaoPlotRegionProxy> TaoPlotProxy::r() const {
  void* ptr;
  tao_plot_struct_get_r(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoPlotRegionProxy(ptr);
}
void TaoPlotProxy::set_r(const TaoPlotRegionProxy& src) {
  tao_plot_struct_set_r(fortran_ptr_, src.get_fortran_ptr());
}
int TaoPlotProxy::ix_plot() const {
  int value;
  tao_plot_struct_get_ix_plot(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_ix_plot(int value) {
  tao_plot_struct_set_ix_plot(fortran_ptr_, value);
}
int TaoPlotProxy::n_curve_pts() const {
  int value;
  tao_plot_struct_get_n_curve_pts(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_n_curve_pts(int value) {
  tao_plot_struct_set_n_curve_pts(fortran_ptr_, value);
}
std::string TaoPlotProxy::type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_plot_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotProxy::set_type(const std::string& value) {
  tao_plot_struct_set_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoPlotProxy::x_axis_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_plot_struct_get_x_axis_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotProxy::set_x_axis_type(const std::string& value) {
  tao_plot_struct_set_x_axis_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoPlotProxy::autoscale_x() const {
  bool value;
  tao_plot_struct_get_autoscale_x(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_autoscale_x(bool value) {
  tao_plot_struct_set_autoscale_x(fortran_ptr_, value);
}
bool TaoPlotProxy::autoscale_y() const {
  bool value;
  tao_plot_struct_get_autoscale_y(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_autoscale_y(bool value) {
  tao_plot_struct_set_autoscale_y(fortran_ptr_, value);
}
bool TaoPlotProxy::autoscale_gang_x() const {
  bool value;
  tao_plot_struct_get_autoscale_gang_x(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_autoscale_gang_x(bool value) {
  tao_plot_struct_set_autoscale_gang_x(fortran_ptr_, value);
}
bool TaoPlotProxy::autoscale_gang_y() const {
  bool value;
  tao_plot_struct_get_autoscale_gang_y(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_autoscale_gang_y(bool value) {
  tao_plot_struct_set_autoscale_gang_y(fortran_ptr_, value);
}
bool TaoPlotProxy::list_with_show_plot_command() const {
  bool value;
  tao_plot_struct_get_list_with_show_plot_command(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_list_with_show_plot_command(bool value) {
  tao_plot_struct_set_list_with_show_plot_command(fortran_ptr_, value);
}
bool TaoPlotProxy::phantom() const {
  bool value;
  tao_plot_struct_get_phantom(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_phantom(bool value) {
  tao_plot_struct_set_phantom(fortran_ptr_, value);
}
bool TaoPlotProxy::default_plot() const {
  bool value;
  tao_plot_struct_get_default_plot(fortran_ptr_, &value);
  return value;
}
void TaoPlotProxy::set_default_plot(bool value) {
  tao_plot_struct_set_default_plot(fortran_ptr_, value);
}
std::string TaoPlotRegionProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_plot_region_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotRegionProxy::set_name(const std::string& value) {
  tao_plot_region_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoPlotProxy TaoPlotRegionProxy::plot() const {
  void* ptr;
  tao_plot_region_struct_get_plot(fortran_ptr_, &ptr);
  return TaoPlotProxy(ptr);
}
void TaoPlotRegionProxy::set_plot(const TaoPlotProxy& src) {
  tao_plot_region_struct_set_plot(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TaoPlotRegionProxy::location() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_plot_region_struct_get_location_info);
}
bool TaoPlotRegionProxy::visible() const {
  bool value;
  tao_plot_region_struct_get_visible(fortran_ptr_, &value);
  return value;
}
void TaoPlotRegionProxy::set_visible(bool value) {
  tao_plot_region_struct_set_visible(fortran_ptr_, value);
}
bool TaoPlotRegionProxy::list_with_show_plot_command() const {
  bool value;
  tao_plot_region_struct_get_list_with_show_plot_command(fortran_ptr_, &value);
  return value;
}
void TaoPlotRegionProxy::set_list_with_show_plot_command(bool value) {
  tao_plot_region_struct_set_list_with_show_plot_command(fortran_ptr_, value);
}
bool TaoPlotRegionProxy::setup_done() const {
  bool value;
  tao_plot_region_struct_get_setup_done(fortran_ptr_, &value);
  return value;
}
void TaoPlotRegionProxy::set_setup_done(bool value) {
  tao_plot_region_struct_set_setup_done(fortran_ptr_, value);
}
std::optional<TaoUniverseProxy> TaoUniversePointerProxy::u() const {
  void* ptr;
  tao_universe_pointer_struct_get_u(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoUniverseProxy(ptr);
}
void TaoUniversePointerProxy::set_u(const TaoUniverseProxy& src) {
  tao_universe_pointer_struct_set_u(fortran_ptr_, src.get_fortran_ptr());
}
TaoGlobalProxy TaoSuperUniverseProxy::global() const {
  void* ptr;
  tao_super_universe_struct_get_global(fortran_ptr_, &ptr);
  return TaoGlobalProxy(ptr);
}
void TaoSuperUniverseProxy::set_global(const TaoGlobalProxy& src) {
  tao_super_universe_struct_set_global(fortran_ptr_, src.get_fortran_ptr());
}
TaoInitProxy TaoSuperUniverseProxy::init() const {
  void* ptr;
  tao_super_universe_struct_get_init(fortran_ptr_, &ptr);
  return TaoInitProxy(ptr);
}
void TaoSuperUniverseProxy::set_init(const TaoInitProxy& src) {
  tao_super_universe_struct_set_init(fortran_ptr_, src.get_fortran_ptr());
}
TaoCommonProxy TaoSuperUniverseProxy::com() const {
  void* ptr;
  tao_super_universe_struct_get_com(fortran_ptr_, &ptr);
  return TaoCommonProxy(ptr);
}
void TaoSuperUniverseProxy::set_com(const TaoCommonProxy& src) {
  tao_super_universe_struct_set_com(fortran_ptr_, src.get_fortran_ptr());
}
TaoPlotPageProxy TaoSuperUniverseProxy::plot_page() const {
  void* ptr;
  tao_super_universe_struct_get_plot_page(fortran_ptr_, &ptr);
  return TaoPlotPageProxy(ptr);
}
void TaoSuperUniverseProxy::set_plot_page(const TaoPlotPageProxy& src) {
  tao_super_universe_struct_set_plot_page(fortran_ptr_, src.get_fortran_ptr());
}
TaoV1VarProxyArray1D TaoSuperUniverseProxy::v1_var() const {
  return ProxyHelpers::get_type_array_1d<TaoV1VarProxyArray1D>(
      fortran_ptr_, tao_super_universe_struct_get_v1_var_info);
}
TaoVarProxyArray1D TaoSuperUniverseProxy::var() const {
  return ProxyHelpers::get_type_array_1d<TaoVarProxyArray1D>(
      fortran_ptr_, tao_super_universe_struct_get_var_info);
}
TaoUniverseProxyArray1D TaoSuperUniverseProxy::u() const {
  return ProxyHelpers::get_type_array_1d<TaoUniverseProxyArray1D>(
      fortran_ptr_, tao_super_universe_struct_get_u_info);
}
FArray1D<int> TaoSuperUniverseProxy::key() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, tao_super_universe_struct_get_key_info);
}
TaoBuildingWallProxy TaoSuperUniverseProxy::building_wall() const {
  void* ptr;
  tao_super_universe_struct_get_building_wall(fortran_ptr_, &ptr);
  return TaoBuildingWallProxy(ptr);
}
void TaoSuperUniverseProxy::set_building_wall(const TaoBuildingWallProxy& src) {
  tao_super_universe_struct_set_building_wall(
      fortran_ptr_, src.get_fortran_ptr());
}
TaoWaveProxy TaoSuperUniverseProxy::wave() const {
  void* ptr;
  tao_super_universe_struct_get_wave(fortran_ptr_, &ptr);
  return TaoWaveProxy(ptr);
}
void TaoSuperUniverseProxy::set_wave(const TaoWaveProxy& src) {
  tao_super_universe_struct_set_wave(fortran_ptr_, src.get_fortran_ptr());
}
int TaoSuperUniverseProxy::n_var_used() const {
  int value;
  tao_super_universe_struct_get_n_var_used(fortran_ptr_, &value);
  return value;
}
void TaoSuperUniverseProxy::set_n_var_used(int value) {
  tao_super_universe_struct_set_n_var_used(fortran_ptr_, value);
}
int TaoSuperUniverseProxy::n_v1_var_used() const {
  int value;
  tao_super_universe_struct_get_n_v1_var_used(fortran_ptr_, &value);
  return value;
}
void TaoSuperUniverseProxy::set_n_v1_var_used(int value) {
  tao_super_universe_struct_set_n_v1_var_used(fortran_ptr_, value);
}
TaoCmdHistoryProxyArray1D TaoSuperUniverseProxy::history() const {
  return ProxyHelpers::get_type_array_1d<TaoCmdHistoryProxyArray1D>(
      fortran_ptr_, tao_super_universe_struct_get_history_info);
}
bool TaoSuperUniverseProxy::initialized() const {
  bool value;
  tao_super_universe_struct_get_initialized(fortran_ptr_, &value);
  return value;
}
void TaoSuperUniverseProxy::set_initialized(bool value) {
  tao_super_universe_struct_set_initialized(fortran_ptr_, value);
}
std::string TaoVarProxy::ele_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_var_struct_get_ele_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarProxy::set_ele_name(const std::string& value) {
  tao_var_struct_set_ele_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoVarProxy::attrib_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_var_struct_get_attrib_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarProxy::set_attrib_name(const std::string& value) {
  tao_var_struct_set_attrib_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoVarProxy::id() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_var_struct_get_id_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarProxy::set_id(const std::string& value) {
  tao_var_struct_set_id(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoVarSlaveProxyArray1D TaoVarProxy::slave() const {
  return ProxyHelpers::get_type_array_1d<TaoVarSlaveProxyArray1D>(
      fortran_ptr_, tao_var_struct_get_slave_info);
}
int TaoVarProxy::ix_v1() const {
  int value;
  tao_var_struct_get_ix_v1(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_ix_v1(int value) {
  tao_var_struct_set_ix_v1(fortran_ptr_, value);
}
int TaoVarProxy::ix_var() const {
  int value;
  tao_var_struct_get_ix_var(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_ix_var(int value) {
  tao_var_struct_set_ix_var(fortran_ptr_, value);
}
int TaoVarProxy::ix_dvar() const {
  int value;
  tao_var_struct_get_ix_dvar(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_ix_dvar(int value) {
  tao_var_struct_set_ix_dvar(fortran_ptr_, value);
}
int TaoVarProxy::ix_attrib() const {
  int value;
  tao_var_struct_get_ix_attrib(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_ix_attrib(int value) {
  tao_var_struct_set_ix_attrib(fortran_ptr_, value);
}
int TaoVarProxy::ix_key_table() const {
  int value;
  tao_var_struct_get_ix_key_table(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_ix_key_table(int value) {
  tao_var_struct_set_ix_key_table(fortran_ptr_, value);
}
double* TaoVarProxy::model_value() const {
  double* ptr;
  tao_var_struct_get_model_value(fortran_ptr_, &ptr);
  return ptr;
}
void TaoVarProxy::set_model_value(double value) {
  tao_var_struct_set_model_value(fortran_ptr_, value);
}
double* TaoVarProxy::base_value() const {
  double* ptr;
  tao_var_struct_get_base_value(fortran_ptr_, &ptr);
  return ptr;
}
void TaoVarProxy::set_base_value(double value) {
  tao_var_struct_set_base_value(fortran_ptr_, value);
}
double TaoVarProxy::design_value() const {
  double value;
  tao_var_struct_get_design_value(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_design_value(double value) {
  tao_var_struct_set_design_value(fortran_ptr_, value);
}
double TaoVarProxy::scratch_value() const {
  double value;
  tao_var_struct_get_scratch_value(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_scratch_value(double value) {
  tao_var_struct_set_scratch_value(fortran_ptr_, value);
}
double TaoVarProxy::old_value() const {
  double value;
  tao_var_struct_get_old_value(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_old_value(double value) {
  tao_var_struct_set_old_value(fortran_ptr_, value);
}
double TaoVarProxy::meas_value() const {
  double value;
  tao_var_struct_get_meas_value(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_meas_value(double value) {
  tao_var_struct_set_meas_value(fortran_ptr_, value);
}
double TaoVarProxy::ref_value() const {
  double value;
  tao_var_struct_get_ref_value(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_ref_value(double value) {
  tao_var_struct_set_ref_value(fortran_ptr_, value);
}
double TaoVarProxy::correction_value() const {
  double value;
  tao_var_struct_get_correction_value(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_correction_value(double value) {
  tao_var_struct_set_correction_value(fortran_ptr_, value);
}
double TaoVarProxy::high_lim() const {
  double value;
  tao_var_struct_get_high_lim(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_high_lim(double value) {
  tao_var_struct_set_high_lim(fortran_ptr_, value);
}
double TaoVarProxy::low_lim() const {
  double value;
  tao_var_struct_get_low_lim(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_low_lim(double value) {
  tao_var_struct_set_low_lim(fortran_ptr_, value);
}
double TaoVarProxy::step() const {
  double value;
  tao_var_struct_get_step(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_step(double value) {
  tao_var_struct_set_step(fortran_ptr_, value);
}
double TaoVarProxy::weight() const {
  double value;
  tao_var_struct_get_weight(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_weight(double value) {
  tao_var_struct_set_weight(fortran_ptr_, value);
}
double TaoVarProxy::delta_merit() const {
  double value;
  tao_var_struct_get_delta_merit(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_delta_merit(double value) {
  tao_var_struct_set_delta_merit(fortran_ptr_, value);
}
double TaoVarProxy::merit() const {
  double value;
  tao_var_struct_get_merit(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_merit(double value) {
  tao_var_struct_set_merit(fortran_ptr_, value);
}
double TaoVarProxy::dMerit_dVar() const {
  double value;
  tao_var_struct_get_dMerit_dVar(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_dMerit_dVar(double value) {
  tao_var_struct_set_dMerit_dVar(fortran_ptr_, value);
}
double TaoVarProxy::key_val0() const {
  double value;
  tao_var_struct_get_key_val0(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_key_val0(double value) {
  tao_var_struct_set_key_val0(fortran_ptr_, value);
}
double TaoVarProxy::key_delta() const {
  double value;
  tao_var_struct_get_key_delta(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_key_delta(double value) {
  tao_var_struct_set_key_delta(fortran_ptr_, value);
}
double TaoVarProxy::s() const {
  double value;
  tao_var_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_s(double value) {
  tao_var_struct_set_s(fortran_ptr_, value);
}
double TaoVarProxy::extend_val() const {
  double value;
  tao_var_struct_get_extend_val(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_extend_val(double value) {
  tao_var_struct_set_extend_val(fortran_ptr_, value);
}
std::string TaoVarProxy::merit_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_var_struct_get_merit_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarProxy::set_merit_type(const std::string& value) {
  tao_var_struct_set_merit_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoVarProxy::exists() const {
  bool value;
  tao_var_struct_get_exists(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_exists(bool value) {
  tao_var_struct_set_exists(fortran_ptr_, value);
}
bool TaoVarProxy::good_var() const {
  bool value;
  tao_var_struct_get_good_var(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_good_var(bool value) {
  tao_var_struct_set_good_var(fortran_ptr_, value);
}
bool TaoVarProxy::good_user() const {
  bool value;
  tao_var_struct_get_good_user(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_good_user(bool value) {
  tao_var_struct_set_good_user(fortran_ptr_, value);
}
bool TaoVarProxy::good_opt() const {
  bool value;
  tao_var_struct_get_good_opt(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_good_opt(bool value) {
  tao_var_struct_set_good_opt(fortran_ptr_, value);
}
bool TaoVarProxy::good_plot() const {
  bool value;
  tao_var_struct_get_good_plot(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_good_plot(bool value) {
  tao_var_struct_set_good_plot(fortran_ptr_, value);
}
bool TaoVarProxy::useit_opt() const {
  bool value;
  tao_var_struct_get_useit_opt(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_useit_opt(bool value) {
  tao_var_struct_set_useit_opt(fortran_ptr_, value);
}
bool TaoVarProxy::useit_plot() const {
  bool value;
  tao_var_struct_get_useit_plot(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_useit_plot(bool value) {
  tao_var_struct_set_useit_plot(fortran_ptr_, value);
}
bool TaoVarProxy::key_bound() const {
  bool value;
  tao_var_struct_get_key_bound(fortran_ptr_, &value);
  return value;
}
void TaoVarProxy::set_key_bound(bool value) {
  tao_var_struct_set_key_bound(fortran_ptr_, value);
}
std::optional<TaoV1VarProxy> TaoVarProxy::v1() const {
  void* ptr;
  tao_var_struct_get_v1(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoV1VarProxy(ptr);
}
void TaoVarProxy::set_v1(const TaoV1VarProxy& src) {
  tao_var_struct_set_v1(fortran_ptr_, src.get_fortran_ptr());
}
int TaoVarSlaveProxy::ix_uni() const {
  int value;
  tao_var_slave_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoVarSlaveProxy::set_ix_uni(int value) {
  tao_var_slave_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoVarSlaveProxy::ix_branch() const {
  int value;
  tao_var_slave_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoVarSlaveProxy::set_ix_branch(int value) {
  tao_var_slave_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoVarSlaveProxy::ix_ele() const {
  int value;
  tao_var_slave_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void TaoVarSlaveProxy::set_ix_ele(int value) {
  tao_var_slave_struct_set_ix_ele(fortran_ptr_, value);
}
double* TaoVarSlaveProxy::model_value() const {
  double* ptr;
  tao_var_slave_struct_get_model_value(fortran_ptr_, &ptr);
  return ptr;
}
void TaoVarSlaveProxy::set_model_value(double value) {
  tao_var_slave_struct_set_model_value(fortran_ptr_, value);
}
double* TaoVarSlaveProxy::base_value() const {
  double* ptr;
  tao_var_slave_struct_get_base_value(fortran_ptr_, &ptr);
  return ptr;
}
void TaoVarSlaveProxy::set_base_value(double value) {
  tao_var_slave_struct_set_base_value(fortran_ptr_, value);
}
std::string TaoLatticeProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_lattice_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoLatticeProxy::set_name(const std::string& value) {
  tao_lattice_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
LatProxy TaoLatticeProxy::lat() const {
  void* ptr;
  tao_lattice_struct_get_lat(fortran_ptr_, &ptr);
  return LatProxy(ptr);
}
void TaoLatticeProxy::set_lat(const LatProxy& src) {
  tao_lattice_struct_set_lat(fortran_ptr_, src.get_fortran_ptr());
}
LatProxy TaoLatticeProxy::high_E_lat() const {
  void* ptr;
  tao_lattice_struct_get_high_E_lat(fortran_ptr_, &ptr);
  return LatProxy(ptr);
}
void TaoLatticeProxy::set_high_E_lat(const LatProxy& src) {
  tao_lattice_struct_set_high_E_lat(fortran_ptr_, src.get_fortran_ptr());
}
LatProxy TaoLatticeProxy::low_E_lat() const {
  void* ptr;
  tao_lattice_struct_get_low_E_lat(fortran_ptr_, &ptr);
  return LatProxy(ptr);
}
void TaoLatticeProxy::set_low_E_lat(const LatProxy& src) {
  tao_lattice_struct_set_low_E_lat(fortran_ptr_, src.get_fortran_ptr());
}
RadIntAllEleProxy TaoLatticeProxy::rad_int_by_ele_ri() const {
  void* ptr;
  tao_lattice_struct_get_rad_int_by_ele_ri(fortran_ptr_, &ptr);
  return RadIntAllEleProxy(ptr);
}
void TaoLatticeProxy::set_rad_int_by_ele_ri(const RadIntAllEleProxy& src) {
  tao_lattice_struct_set_rad_int_by_ele_ri(fortran_ptr_, src.get_fortran_ptr());
}
RadIntAllEleProxy TaoLatticeProxy::rad_int_by_ele_6d() const {
  void* ptr;
  tao_lattice_struct_get_rad_int_by_ele_6d(fortran_ptr_, &ptr);
  return RadIntAllEleProxy(ptr);
}
void TaoLatticeProxy::set_rad_int_by_ele_6d(const RadIntAllEleProxy& src) {
  tao_lattice_struct_set_rad_int_by_ele_6d(fortran_ptr_, src.get_fortran_ptr());
}
TaoLatticeBranchProxyArray1D TaoLatticeProxy::tao_branch() const {
  return ProxyHelpers::get_type_array_1d<TaoLatticeBranchProxyArray1D>(
      fortran_ptr_, tao_lattice_struct_get_tao_branch_info);
}
std::string TaoBeamUniProxy::saved_at() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_beam_uni_struct_get_saved_at_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamUniProxy::set_saved_at(const std::string& value) {
  tao_beam_uni_struct_set_saved_at(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoBeamUniProxy::dump_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_beam_uni_struct_get_dump_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamUniProxy::set_dump_file(const std::string& value) {
  tao_beam_uni_struct_set_dump_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoBeamUniProxy::dump_at() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_beam_uni_struct_get_dump_at_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamUniProxy::set_dump_at(const std::string& value) {
  tao_beam_uni_struct_set_dump_at(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoBeamUniProxy::track_beam_in_universe() const {
  bool value;
  tao_beam_uni_struct_get_track_beam_in_universe(fortran_ptr_, &value);
  return value;
}
void TaoBeamUniProxy::set_track_beam_in_universe(bool value) {
  tao_beam_uni_struct_set_track_beam_in_universe(fortran_ptr_, value);
}
bool TaoBeamUniProxy::always_reinit() const {
  bool value;
  tao_beam_uni_struct_get_always_reinit(fortran_ptr_, &value);
  return value;
}
void TaoBeamUniProxy::set_always_reinit(bool value) {
  tao_beam_uni_struct_set_always_reinit(fortran_ptr_, value);
}
ApertureParamProxy TaoDynamicApertureProxy::param() const {
  void* ptr;
  tao_dynamic_aperture_struct_get_param(fortran_ptr_, &ptr);
  return ApertureParamProxy(ptr);
}
void TaoDynamicApertureProxy::set_param(const ApertureParamProxy& src) {
  tao_dynamic_aperture_struct_set_param(fortran_ptr_, src.get_fortran_ptr());
}
ApertureScanProxyArray1D TaoDynamicApertureProxy::scan() const {
  return ProxyHelpers::get_type_array_1d<ApertureScanProxyArray1D>(
      fortran_ptr_, tao_dynamic_aperture_struct_get_scan_info);
}
FArray1D<double> TaoDynamicApertureProxy::pz() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_dynamic_aperture_struct_get_pz_info);
}
double TaoDynamicApertureProxy::ellipse_scale() const {
  double value;
  tao_dynamic_aperture_struct_get_ellipse_scale(fortran_ptr_, &value);
  return value;
}
void TaoDynamicApertureProxy::set_ellipse_scale(double value) {
  tao_dynamic_aperture_struct_set_ellipse_scale(fortran_ptr_, value);
}
double TaoDynamicApertureProxy::a_emit() const {
  double value;
  tao_dynamic_aperture_struct_get_a_emit(fortran_ptr_, &value);
  return value;
}
void TaoDynamicApertureProxy::set_a_emit(double value) {
  tao_dynamic_aperture_struct_set_a_emit(fortran_ptr_, value);
}
double TaoDynamicApertureProxy::b_emit() const {
  double value;
  tao_dynamic_aperture_struct_get_b_emit(fortran_ptr_, &value);
  return value;
}
void TaoDynamicApertureProxy::set_b_emit(double value) {
  tao_dynamic_aperture_struct_set_b_emit(fortran_ptr_, value);
}
TaoModelElementProxyArray1D TaoModelBranchProxy::ele() const {
  return ProxyHelpers::get_type_array_1d<TaoModelElementProxyArray1D>(
      fortran_ptr_, tao_model_branch_struct_get_ele_info);
}
TaoBeamBranchProxy TaoModelBranchProxy::beam() const {
  void* ptr;
  tao_model_branch_struct_get_beam(fortran_ptr_, &ptr);
  return TaoBeamBranchProxy(ptr);
}
void TaoModelBranchProxy::set_beam(const TaoBeamBranchProxy& src) {
  tao_model_branch_struct_set_beam(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoSpinMapProxy::valid() const {
  bool value;
  tao_spin_map_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapProxy::set_valid(bool value) {
  tao_spin_map_struct_set_valid(fortran_ptr_, value);
}
SpinOrbitMap1Proxy TaoSpinMapProxy::map1() const {
  void* ptr;
  tao_spin_map_struct_get_map1(fortran_ptr_, &ptr);
  return SpinOrbitMap1Proxy(ptr);
}
void TaoSpinMapProxy::set_map1(const SpinOrbitMap1Proxy& src) {
  tao_spin_map_struct_set_map1(fortran_ptr_, src.get_fortran_ptr());
}
SpinAxisProxy TaoSpinMapProxy::axis_input() const {
  void* ptr;
  tao_spin_map_struct_get_axis_input(fortran_ptr_, &ptr);
  return SpinAxisProxy(ptr);
}
void TaoSpinMapProxy::set_axis_input(const SpinAxisProxy& src) {
  tao_spin_map_struct_set_axis_input(fortran_ptr_, src.get_fortran_ptr());
}
SpinAxisProxy TaoSpinMapProxy::axis0() const {
  void* ptr;
  tao_spin_map_struct_get_axis0(fortran_ptr_, &ptr);
  return SpinAxisProxy(ptr);
}
void TaoSpinMapProxy::set_axis0(const SpinAxisProxy& src) {
  tao_spin_map_struct_set_axis0(fortran_ptr_, src.get_fortran_ptr());
}
SpinAxisProxy TaoSpinMapProxy::axis1() const {
  void* ptr;
  tao_spin_map_struct_get_axis1(fortran_ptr_, &ptr);
  return SpinAxisProxy(ptr);
}
void TaoSpinMapProxy::set_axis1(const SpinAxisProxy& src) {
  tao_spin_map_struct_set_axis1(fortran_ptr_, src.get_fortran_ptr());
}
int TaoSpinMapProxy::ix_ele() const {
  int value;
  tao_spin_map_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapProxy::set_ix_ele(int value) {
  tao_spin_map_struct_set_ix_ele(fortran_ptr_, value);
}
int TaoSpinMapProxy::ix_ref() const {
  int value;
  tao_spin_map_struct_get_ix_ref(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapProxy::set_ix_ref(int value) {
  tao_spin_map_struct_set_ix_ref(fortran_ptr_, value);
}
int TaoSpinMapProxy::ix_uni() const {
  int value;
  tao_spin_map_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapProxy::set_ix_uni(int value) {
  tao_spin_map_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoSpinMapProxy::ix_branch() const {
  int value;
  tao_spin_map_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapProxy::set_ix_branch(int value) {
  tao_spin_map_struct_set_ix_branch(fortran_ptr_, value);
}
FArray2D<double> TaoSpinMapProxy::mat8() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_spin_map_struct_get_mat8_info);
}
std::string TaoDataProxy::ele_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_ele_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_ele_name(const std::string& value) {
  tao_data_struct_set_ele_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::ele_start_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_ele_start_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_ele_start_name(const std::string& value) {
  tao_data_struct_set_ele_start_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::ele_ref_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_ele_ref_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_ele_ref_name(const std::string& value) {
  tao_data_struct_set_ele_ref_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::data_type() const {
  return ProxyHelpers::get_string(
      fortran_ptr_, tao_data_struct_get_data_type_info);
}
void TaoDataProxy::set_data_type(const std::string& value) {
  tao_data_struct_set_data_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::merit_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_merit_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_merit_type(const std::string& value) {
  tao_data_struct_set_merit_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::id() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_id_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_id(const std::string& value) {
  tao_data_struct_set_id(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::data_source() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_data_source_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_data_source(const std::string& value) {
  tao_data_struct_set_data_source(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataProxy::why_invalid() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_data_struct_get_why_invalid_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataProxy::set_why_invalid(const std::string& value) {
  tao_data_struct_set_why_invalid(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoDataProxy::ix_uni() const {
  int value;
  tao_data_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_uni(int value) {
  tao_data_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoDataProxy::ix_bunch() const {
  int value;
  tao_data_struct_get_ix_bunch(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_bunch(int value) {
  tao_data_struct_set_ix_bunch(fortran_ptr_, value);
}
int TaoDataProxy::ix_branch() const {
  int value;
  tao_data_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_branch(int value) {
  tao_data_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoDataProxy::ix_ele() const {
  int value;
  tao_data_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_ele(int value) {
  tao_data_struct_set_ix_ele(fortran_ptr_, value);
}
int TaoDataProxy::ix_ele_start() const {
  int value;
  tao_data_struct_get_ix_ele_start(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_ele_start(int value) {
  tao_data_struct_set_ix_ele_start(fortran_ptr_, value);
}
int TaoDataProxy::ix_ele_ref() const {
  int value;
  tao_data_struct_get_ix_ele_ref(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_ele_ref(int value) {
  tao_data_struct_set_ix_ele_ref(fortran_ptr_, value);
}
int TaoDataProxy::ix_ele_merit() const {
  int value;
  tao_data_struct_get_ix_ele_merit(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_ele_merit(int value) {
  tao_data_struct_set_ix_ele_merit(fortran_ptr_, value);
}
int TaoDataProxy::ix_d1() const {
  int value;
  tao_data_struct_get_ix_d1(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_d1(int value) {
  tao_data_struct_set_ix_d1(fortran_ptr_, value);
}
int TaoDataProxy::ix_data() const {
  int value;
  tao_data_struct_get_ix_data(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_data(int value) {
  tao_data_struct_set_ix_data(fortran_ptr_, value);
}
int TaoDataProxy::ix_dModel() const {
  int value;
  tao_data_struct_get_ix_dModel(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ix_dModel(int value) {
  tao_data_struct_set_ix_dModel(fortran_ptr_, value);
}
int TaoDataProxy::eval_point() const {
  int value;
  tao_data_struct_get_eval_point(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_eval_point(int value) {
  tao_data_struct_set_eval_point(fortran_ptr_, value);
}
double TaoDataProxy::meas_value() const {
  double value;
  tao_data_struct_get_meas_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_meas_value(double value) {
  tao_data_struct_set_meas_value(fortran_ptr_, value);
}
double TaoDataProxy::ref_value() const {
  double value;
  tao_data_struct_get_ref_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ref_value(double value) {
  tao_data_struct_set_ref_value(fortran_ptr_, value);
}
double TaoDataProxy::model_value() const {
  double value;
  tao_data_struct_get_model_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_model_value(double value) {
  tao_data_struct_set_model_value(fortran_ptr_, value);
}
double TaoDataProxy::design_value() const {
  double value;
  tao_data_struct_get_design_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_design_value(double value) {
  tao_data_struct_set_design_value(fortran_ptr_, value);
}
double TaoDataProxy::old_value() const {
  double value;
  tao_data_struct_get_old_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_old_value(double value) {
  tao_data_struct_set_old_value(fortran_ptr_, value);
}
double TaoDataProxy::base_value() const {
  double value;
  tao_data_struct_get_base_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_base_value(double value) {
  tao_data_struct_set_base_value(fortran_ptr_, value);
}
double TaoDataProxy::error_rms() const {
  double value;
  tao_data_struct_get_error_rms(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_error_rms(double value) {
  tao_data_struct_set_error_rms(fortran_ptr_, value);
}
double TaoDataProxy::delta_merit() const {
  double value;
  tao_data_struct_get_delta_merit(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_delta_merit(double value) {
  tao_data_struct_set_delta_merit(fortran_ptr_, value);
}
double TaoDataProxy::weight() const {
  double value;
  tao_data_struct_get_weight(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_weight(double value) {
  tao_data_struct_set_weight(fortran_ptr_, value);
}
double TaoDataProxy::invalid_value() const {
  double value;
  tao_data_struct_get_invalid_value(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_invalid_value(double value) {
  tao_data_struct_set_invalid_value(fortran_ptr_, value);
}
double TaoDataProxy::merit() const {
  double value;
  tao_data_struct_get_merit(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_merit(double value) {
  tao_data_struct_set_merit(fortran_ptr_, value);
}
double TaoDataProxy::s() const {
  double value;
  tao_data_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_s(double value) {
  tao_data_struct_set_s(fortran_ptr_, value);
}
double TaoDataProxy::s_offset() const {
  double value;
  tao_data_struct_get_s_offset(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_s_offset(double value) {
  tao_data_struct_set_s_offset(fortran_ptr_, value);
}
double TaoDataProxy::ref_s_offset() const {
  double value;
  tao_data_struct_get_ref_s_offset(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_ref_s_offset(double value) {
  tao_data_struct_set_ref_s_offset(fortran_ptr_, value);
}
bool TaoDataProxy::err_message_printed() const {
  bool value;
  tao_data_struct_get_err_message_printed(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_err_message_printed(bool value) {
  tao_data_struct_set_err_message_printed(fortran_ptr_, value);
}
bool TaoDataProxy::exists() const {
  bool value;
  tao_data_struct_get_exists(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_exists(bool value) {
  tao_data_struct_set_exists(fortran_ptr_, value);
}
bool TaoDataProxy::good_model() const {
  bool value;
  tao_data_struct_get_good_model(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_model(bool value) {
  tao_data_struct_set_good_model(fortran_ptr_, value);
}
bool TaoDataProxy::good_base() const {
  bool value;
  tao_data_struct_get_good_base(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_base(bool value) {
  tao_data_struct_set_good_base(fortran_ptr_, value);
}
bool TaoDataProxy::good_design() const {
  bool value;
  tao_data_struct_get_good_design(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_design(bool value) {
  tao_data_struct_set_good_design(fortran_ptr_, value);
}
bool TaoDataProxy::good_meas() const {
  bool value;
  tao_data_struct_get_good_meas(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_meas(bool value) {
  tao_data_struct_set_good_meas(fortran_ptr_, value);
}
bool TaoDataProxy::good_ref() const {
  bool value;
  tao_data_struct_get_good_ref(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_ref(bool value) {
  tao_data_struct_set_good_ref(fortran_ptr_, value);
}
bool TaoDataProxy::good_user() const {
  bool value;
  tao_data_struct_get_good_user(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_user(bool value) {
  tao_data_struct_set_good_user(fortran_ptr_, value);
}
bool TaoDataProxy::good_opt() const {
  bool value;
  tao_data_struct_get_good_opt(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_opt(bool value) {
  tao_data_struct_set_good_opt(fortran_ptr_, value);
}
bool TaoDataProxy::good_plot() const {
  bool value;
  tao_data_struct_get_good_plot(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_good_plot(bool value) {
  tao_data_struct_set_good_plot(fortran_ptr_, value);
}
bool TaoDataProxy::useit_plot() const {
  bool value;
  tao_data_struct_get_useit_plot(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_useit_plot(bool value) {
  tao_data_struct_set_useit_plot(fortran_ptr_, value);
}
bool TaoDataProxy::useit_opt() const {
  bool value;
  tao_data_struct_get_useit_opt(fortran_ptr_, &value);
  return value;
}
void TaoDataProxy::set_useit_opt(bool value) {
  tao_data_struct_set_useit_opt(fortran_ptr_, value);
}
TaoSpinMapProxy TaoDataProxy::spin_map() const {
  void* ptr;
  tao_data_struct_get_spin_map(fortran_ptr_, &ptr);
  return TaoSpinMapProxy(ptr);
}
void TaoDataProxy::set_spin_map(const TaoSpinMapProxy& src) {
  tao_data_struct_set_spin_map(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoD1DataProxy> TaoDataProxy::d1() const {
  void* ptr;
  tao_data_struct_get_d1(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoD1DataProxy(ptr);
}
void TaoDataProxy::set_d1(const TaoD1DataProxy& src) {
  tao_data_struct_set_d1(fortran_ptr_, src.get_fortran_ptr());
}
double TaoPingScaleProxy::a_mode_meas() const {
  double value;
  tao_ping_scale_struct_get_a_mode_meas(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleProxy::set_a_mode_meas(double value) {
  tao_ping_scale_struct_set_a_mode_meas(fortran_ptr_, value);
}
double TaoPingScaleProxy::a_mode_ref() const {
  double value;
  tao_ping_scale_struct_get_a_mode_ref(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleProxy::set_a_mode_ref(double value) {
  tao_ping_scale_struct_set_a_mode_ref(fortran_ptr_, value);
}
double TaoPingScaleProxy::b_mode_meas() const {
  double value;
  tao_ping_scale_struct_get_b_mode_meas(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleProxy::set_b_mode_meas(double value) {
  tao_ping_scale_struct_set_b_mode_meas(fortran_ptr_, value);
}
double TaoPingScaleProxy::b_mode_ref() const {
  double value;
  tao_ping_scale_struct_get_b_mode_ref(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleProxy::set_b_mode_ref(double value) {
  tao_ping_scale_struct_set_b_mode_ref(fortran_ptr_, value);
}
int TaoUniverseCalcProxy::srdt_for_data() const {
  int value;
  tao_universe_calc_struct_get_srdt_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_srdt_for_data(int value) {
  tao_universe_calc_struct_set_srdt_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::rad_int_for_data() const {
  bool value;
  tao_universe_calc_struct_get_rad_int_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_rad_int_for_data(bool value) {
  tao_universe_calc_struct_set_rad_int_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::rad_int_for_plotting() const {
  bool value;
  tao_universe_calc_struct_get_rad_int_for_plotting(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_rad_int_for_plotting(bool value) {
  tao_universe_calc_struct_set_rad_int_for_plotting(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::chrom_for_data() const {
  bool value;
  tao_universe_calc_struct_get_chrom_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_chrom_for_data(bool value) {
  tao_universe_calc_struct_set_chrom_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::chrom_for_plotting() const {
  bool value;
  tao_universe_calc_struct_get_chrom_for_plotting(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_chrom_for_plotting(bool value) {
  tao_universe_calc_struct_set_chrom_for_plotting(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::lat_sigma_for_data() const {
  bool value;
  tao_universe_calc_struct_get_lat_sigma_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_lat_sigma_for_data(bool value) {
  tao_universe_calc_struct_set_lat_sigma_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::lat_sigma_for_plotting() const {
  bool value;
  tao_universe_calc_struct_get_lat_sigma_for_plotting(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_lat_sigma_for_plotting(bool value) {
  tao_universe_calc_struct_set_lat_sigma_for_plotting(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::dynamic_aperture() const {
  bool value;
  tao_universe_calc_struct_get_dynamic_aperture(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_dynamic_aperture(bool value) {
  tao_universe_calc_struct_set_dynamic_aperture(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::one_turn_map() const {
  bool value;
  tao_universe_calc_struct_get_one_turn_map(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_one_turn_map(bool value) {
  tao_universe_calc_struct_set_one_turn_map(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::lattice() const {
  bool value;
  tao_universe_calc_struct_get_lattice(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_lattice(bool value) {
  tao_universe_calc_struct_set_lattice(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::twiss() const {
  bool value;
  tao_universe_calc_struct_get_twiss(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_twiss(bool value) {
  tao_universe_calc_struct_set_twiss(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::track() const {
  bool value;
  tao_universe_calc_struct_get_track(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_track(bool value) {
  tao_universe_calc_struct_set_track(fortran_ptr_, value);
}
bool TaoUniverseCalcProxy::spin_matrices() const {
  bool value;
  tao_universe_calc_struct_get_spin_matrices(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcProxy::set_spin_matrices(bool value) {
  tao_universe_calc_struct_set_spin_matrices(fortran_ptr_, value);
}
LatEleOrderArrayProxyArray1D LatEleOrderProxy::branch() const {
  return ProxyHelpers::get_type_array_1d<LatEleOrderArrayProxyArray1D>(
      fortran_ptr_, lat_ele_order_struct_get_branch_info);
}
std::string TaoTitleProxy::string() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_title_struct_get_string_info);
  return std::string(arr.data(), arr.size());
}
void TaoTitleProxy::set_string(const std::string& value) {
  tao_title_struct_set_string(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoTitleProxy::x() const {
  double value;
  tao_title_struct_get_x(fortran_ptr_, &value);
  return value;
}
void TaoTitleProxy::set_x(double value) {
  tao_title_struct_set_x(fortran_ptr_, value);
}
double TaoTitleProxy::y() const {
  double value;
  tao_title_struct_get_y(fortran_ptr_, &value);
  return value;
}
void TaoTitleProxy::set_y(double value) {
  tao_title_struct_set_y(fortran_ptr_, value);
}
std::string TaoTitleProxy::units() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_title_struct_get_units_info);
  return std::string(arr.data(), arr.size());
}
void TaoTitleProxy::set_units(const std::string& value) {
  tao_title_struct_set_units(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoTitleProxy::justify() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_title_struct_get_justify_info);
  return std::string(arr.data(), arr.size());
}
void TaoTitleProxy::set_justify(const std::string& value) {
  tao_title_struct_set_justify(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoTitleProxy::draw_it() const {
  bool value;
  tao_title_struct_get_draw_it(fortran_ptr_, &value);
  return value;
}
void TaoTitleProxy::set_draw_it(bool value) {
  tao_title_struct_set_draw_it(fortran_ptr_, value);
}
double QpRectProxy::x1() const {
  double value;
  qp_rect_struct_get_x1(fortran_ptr_, &value);
  return value;
}
void QpRectProxy::set_x1(double value) {
  qp_rect_struct_set_x1(fortran_ptr_, value);
}
double QpRectProxy::x2() const {
  double value;
  qp_rect_struct_get_x2(fortran_ptr_, &value);
  return value;
}
void QpRectProxy::set_x2(double value) {
  qp_rect_struct_set_x2(fortran_ptr_, value);
}
double QpRectProxy::y1() const {
  double value;
  qp_rect_struct_get_y1(fortran_ptr_, &value);
  return value;
}
void QpRectProxy::set_y1(double value) {
  qp_rect_struct_set_y1(fortran_ptr_, value);
}
double QpRectProxy::y2() const {
  double value;
  qp_rect_struct_get_y2(fortran_ptr_, &value);
  return value;
}
void QpRectProxy::set_y2(double value) {
  qp_rect_struct_set_y2(fortran_ptr_, value);
}
std::string QpRectProxy::units() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_rect_struct_get_units_info);
  return std::string(arr.data(), arr.size());
}
void QpRectProxy::set_units(const std::string& value) {
  qp_rect_struct_set_units(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoEleShapeProxyArray1D TaoDrawingProxy::ele_shape() const {
  return ProxyHelpers::get_type_array_1d<TaoEleShapeProxyArray1D>(
      fortran_ptr_, tao_drawing_struct_get_ele_shape_info);
}
std::string TaoShapePatternProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_shape_pattern_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoShapePatternProxy::set_name(const std::string& value) {
  tao_shape_pattern_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
QpLineProxy TaoShapePatternProxy::line() const {
  void* ptr;
  tao_shape_pattern_struct_get_line(fortran_ptr_, &ptr);
  return QpLineProxy(ptr);
}
void TaoShapePatternProxy::set_line(const QpLineProxy& src) {
  tao_shape_pattern_struct_set_line(fortran_ptr_, src.get_fortran_ptr());
}
TaoShapePatternPointProxyArray1D TaoShapePatternProxy::pt() const {
  return ProxyHelpers::get_type_array_1d<TaoShapePatternPointProxyArray1D>(
      fortran_ptr_, tao_shape_pattern_struct_get_pt_info);
}
double TaoShapePatternPointProxy::s() const {
  double value;
  tao_shape_pattern_point_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoShapePatternPointProxy::set_s(double value) {
  tao_shape_pattern_point_struct_set_s(fortran_ptr_, value);
}
double TaoShapePatternPointProxy::y() const {
  double value;
  tao_shape_pattern_point_struct_get_y(fortran_ptr_, &value);
  return value;
}
void TaoShapePatternPointProxy::set_y(double value) {
  tao_shape_pattern_point_struct_set_y(fortran_ptr_, value);
}
double TaoShapePatternPointProxy::radius() const {
  double value;
  tao_shape_pattern_point_struct_get_radius(fortran_ptr_, &value);
  return value;
}
void TaoShapePatternPointProxy::set_radius(double value) {
  tao_shape_pattern_point_struct_set_radius(fortran_ptr_, value);
}
std::string QpAxisProxy::label() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_axis_struct_get_label_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisProxy::set_label(const std::string& value) {
  qp_axis_struct_set_label(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double QpAxisProxy::min() const {
  double value;
  qp_axis_struct_get_min(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_min(double value) {
  qp_axis_struct_set_min(fortran_ptr_, value);
}
double QpAxisProxy::max() const {
  double value;
  qp_axis_struct_get_max(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_max(double value) {
  qp_axis_struct_set_max(fortran_ptr_, value);
}
double QpAxisProxy::tick_min() const {
  double value;
  qp_axis_struct_get_tick_min(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_tick_min(double value) {
  qp_axis_struct_set_tick_min(fortran_ptr_, value);
}
double QpAxisProxy::tick_max() const {
  double value;
  qp_axis_struct_get_tick_max(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_tick_max(double value) {
  qp_axis_struct_set_tick_max(fortran_ptr_, value);
}
double QpAxisProxy::eval_min() const {
  double value;
  qp_axis_struct_get_eval_min(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_eval_min(double value) {
  qp_axis_struct_set_eval_min(fortran_ptr_, value);
}
double QpAxisProxy::eval_max() const {
  double value;
  qp_axis_struct_get_eval_max(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_eval_max(double value) {
  qp_axis_struct_set_eval_max(fortran_ptr_, value);
}
double QpAxisProxy::dtick() const {
  double value;
  qp_axis_struct_get_dtick(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_dtick(double value) {
  qp_axis_struct_set_dtick(fortran_ptr_, value);
}
double QpAxisProxy::number_offset() const {
  double value;
  qp_axis_struct_get_number_offset(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_number_offset(double value) {
  qp_axis_struct_set_number_offset(fortran_ptr_, value);
}
double QpAxisProxy::label_offset() const {
  double value;
  qp_axis_struct_get_label_offset(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_label_offset(double value) {
  qp_axis_struct_set_label_offset(fortran_ptr_, value);
}
double QpAxisProxy::major_tick_len() const {
  double value;
  qp_axis_struct_get_major_tick_len(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_major_tick_len(double value) {
  qp_axis_struct_set_major_tick_len(fortran_ptr_, value);
}
double QpAxisProxy::minor_tick_len() const {
  double value;
  qp_axis_struct_get_minor_tick_len(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_minor_tick_len(double value) {
  qp_axis_struct_set_minor_tick_len(fortran_ptr_, value);
}
std::string QpAxisProxy::label_color() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_axis_struct_get_label_color_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisProxy::set_label_color(const std::string& value) {
  qp_axis_struct_set_label_color(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpAxisProxy::major_div() const {
  int value;
  qp_axis_struct_get_major_div(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_major_div(int value) {
  qp_axis_struct_set_major_div(fortran_ptr_, value);
}
int QpAxisProxy::major_div_nominal() const {
  int value;
  qp_axis_struct_get_major_div_nominal(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_major_div_nominal(int value) {
  qp_axis_struct_set_major_div_nominal(fortran_ptr_, value);
}
int QpAxisProxy::minor_div() const {
  int value;
  qp_axis_struct_get_minor_div(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_minor_div(int value) {
  qp_axis_struct_set_minor_div(fortran_ptr_, value);
}
int QpAxisProxy::minor_div_max() const {
  int value;
  qp_axis_struct_get_minor_div_max(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_minor_div_max(int value) {
  qp_axis_struct_set_minor_div_max(fortran_ptr_, value);
}
int QpAxisProxy::places() const {
  int value;
  qp_axis_struct_get_places(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_places(int value) {
  qp_axis_struct_set_places(fortran_ptr_, value);
}
std::string QpAxisProxy::type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_axis_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisProxy::set_type(const std::string& value) {
  qp_axis_struct_set_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpAxisProxy::bounds() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_axis_struct_get_bounds_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisProxy::set_bounds(const std::string& value) {
  qp_axis_struct_set_bounds(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpAxisProxy::tick_side() const {
  int value;
  qp_axis_struct_get_tick_side(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_tick_side(int value) {
  qp_axis_struct_set_tick_side(fortran_ptr_, value);
}
int QpAxisProxy::number_side() const {
  int value;
  qp_axis_struct_get_number_side(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_number_side(int value) {
  qp_axis_struct_set_number_side(fortran_ptr_, value);
}
bool QpAxisProxy::draw_label() const {
  bool value;
  qp_axis_struct_get_draw_label(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_draw_label(bool value) {
  qp_axis_struct_set_draw_label(fortran_ptr_, value);
}
bool QpAxisProxy::draw_numbers() const {
  bool value;
  qp_axis_struct_get_draw_numbers(fortran_ptr_, &value);
  return value;
}
void QpAxisProxy::set_draw_numbers(bool value) {
  qp_axis_struct_set_draw_numbers(fortran_ptr_, value);
}
double QpLegendProxy::row_spacing() const {
  double value;
  qp_legend_struct_get_row_spacing(fortran_ptr_, &value);
  return value;
}
void QpLegendProxy::set_row_spacing(double value) {
  qp_legend_struct_set_row_spacing(fortran_ptr_, value);
}
double QpLegendProxy::line_length() const {
  double value;
  qp_legend_struct_get_line_length(fortran_ptr_, &value);
  return value;
}
void QpLegendProxy::set_line_length(double value) {
  qp_legend_struct_set_line_length(fortran_ptr_, value);
}
double QpLegendProxy::text_offset() const {
  double value;
  qp_legend_struct_get_text_offset(fortran_ptr_, &value);
  return value;
}
void QpLegendProxy::set_text_offset(double value) {
  qp_legend_struct_set_text_offset(fortran_ptr_, value);
}
bool QpLegendProxy::draw_line() const {
  bool value;
  qp_legend_struct_get_draw_line(fortran_ptr_, &value);
  return value;
}
void QpLegendProxy::set_draw_line(bool value) {
  qp_legend_struct_set_draw_line(fortran_ptr_, value);
}
bool QpLegendProxy::draw_symbol() const {
  bool value;
  qp_legend_struct_get_draw_symbol(fortran_ptr_, &value);
  return value;
}
void QpLegendProxy::set_draw_symbol(bool value) {
  qp_legend_struct_set_draw_symbol(fortran_ptr_, value);
}
bool QpLegendProxy::draw_text() const {
  bool value;
  qp_legend_struct_get_draw_text(fortran_ptr_, &value);
  return value;
}
void QpLegendProxy::set_draw_text(bool value) {
  qp_legend_struct_set_draw_text(fortran_ptr_, value);
}
double QpPointProxy::x() const {
  double value;
  qp_point_struct_get_x(fortran_ptr_, &value);
  return value;
}
void QpPointProxy::set_x(double value) {
  qp_point_struct_set_x(fortran_ptr_, value);
}
double QpPointProxy::y() const {
  double value;
  qp_point_struct_get_y(fortran_ptr_, &value);
  return value;
}
void QpPointProxy::set_y(double value) {
  qp_point_struct_set_y(fortran_ptr_, value);
}
std::string QpPointProxy::units() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_point_struct_get_units_info);
  return std::string(arr.data(), arr.size());
}
void QpPointProxy::set_units(const std::string& value) {
  qp_point_struct_set_units(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpLineProxy::width() const {
  int value;
  qp_line_struct_get_width(fortran_ptr_, &value);
  return value;
}
void QpLineProxy::set_width(int value) {
  qp_line_struct_set_width(fortran_ptr_, value);
}
std::string QpLineProxy::color() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_line_struct_get_color_info);
  return std::string(arr.data(), arr.size());
}
void QpLineProxy::set_color(const std::string& value) {
  qp_line_struct_set_color(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpLineProxy::pattern() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_line_struct_get_pattern_info);
  return std::string(arr.data(), arr.size());
}
void QpLineProxy::set_pattern(const std::string& value) {
  qp_line_struct_set_pattern(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpSymbolProxy::type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_symbol_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void QpSymbolProxy::set_type(const std::string& value) {
  qp_symbol_struct_set_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double QpSymbolProxy::height() const {
  double value;
  qp_symbol_struct_get_height(fortran_ptr_, &value);
  return value;
}
void QpSymbolProxy::set_height(double value) {
  qp_symbol_struct_set_height(fortran_ptr_, value);
}
std::string QpSymbolProxy::color() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_symbol_struct_get_color_info);
  return std::string(arr.data(), arr.size());
}
void QpSymbolProxy::set_color(const std::string& value) {
  qp_symbol_struct_set_color(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpSymbolProxy::fill_pattern() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, qp_symbol_struct_get_fill_pattern_info);
  return std::string(arr.data(), arr.size());
}
void QpSymbolProxy::set_fill_pattern(const std::string& value) {
  qp_symbol_struct_set_fill_pattern(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpSymbolProxy::line_width() const {
  int value;
  qp_symbol_struct_get_line_width(fortran_ptr_, &value);
  return value;
}
void QpSymbolProxy::set_line_width(int value) {
  qp_symbol_struct_set_line_width(fortran_ptr_, value);
}
std::string TaoFloorPlanProxy::view() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_floor_plan_struct_get_view_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanProxy::set_view(const std::string& value) {
  tao_floor_plan_struct_set_view(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoFloorPlanProxy::rotation() const {
  double value;
  tao_floor_plan_struct_get_rotation(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_rotation(double value) {
  tao_floor_plan_struct_set_rotation(fortran_ptr_, value);
}
bool TaoFloorPlanProxy::correct_distortion() const {
  bool value;
  tao_floor_plan_struct_get_correct_distortion(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_correct_distortion(bool value) {
  tao_floor_plan_struct_set_correct_distortion(fortran_ptr_, value);
}
bool TaoFloorPlanProxy::flip_label_side() const {
  bool value;
  tao_floor_plan_struct_get_flip_label_side(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_flip_label_side(bool value) {
  tao_floor_plan_struct_set_flip_label_side(fortran_ptr_, value);
}
bool TaoFloorPlanProxy::size_is_absolute() const {
  bool value;
  tao_floor_plan_struct_get_size_is_absolute(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_size_is_absolute(bool value) {
  tao_floor_plan_struct_set_size_is_absolute(fortran_ptr_, value);
}
bool TaoFloorPlanProxy::draw_only_first_pass() const {
  bool value;
  tao_floor_plan_struct_get_draw_only_first_pass(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_draw_only_first_pass(bool value) {
  tao_floor_plan_struct_set_draw_only_first_pass(fortran_ptr_, value);
}
bool TaoFloorPlanProxy::draw_building_wall() const {
  bool value;
  tao_floor_plan_struct_get_draw_building_wall(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_draw_building_wall(bool value) {
  tao_floor_plan_struct_set_draw_building_wall(fortran_ptr_, value);
}
double TaoFloorPlanProxy::orbit_scale() const {
  double value;
  tao_floor_plan_struct_get_orbit_scale(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_orbit_scale(double value) {
  tao_floor_plan_struct_set_orbit_scale(fortran_ptr_, value);
}
std::string TaoFloorPlanProxy::orbit_color() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_floor_plan_struct_get_orbit_color_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanProxy::set_orbit_color(const std::string& value) {
  tao_floor_plan_struct_set_orbit_color(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoFloorPlanProxy::orbit_pattern() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_floor_plan_struct_get_orbit_pattern_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanProxy::set_orbit_pattern(const std::string& value) {
  tao_floor_plan_struct_set_orbit_pattern(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoFloorPlanProxy::orbit_lattice() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_floor_plan_struct_get_orbit_lattice_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanProxy::set_orbit_lattice(const std::string& value) {
  tao_floor_plan_struct_set_orbit_lattice(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoFloorPlanProxy::orbit_width() const {
  int value;
  tao_floor_plan_struct_get_orbit_width(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanProxy::set_orbit_width(int value) {
  tao_floor_plan_struct_set_orbit_width(fortran_ptr_, value);
}
std::string TaoV1VarProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_v1_var_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoV1VarProxy::set_name(const std::string& value) {
  tao_v1_var_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoV1VarProxy::ix_v1_var() const {
  int value;
  tao_v1_var_struct_get_ix_v1_var(fortran_ptr_, &value);
  return value;
}
void TaoV1VarProxy::set_ix_v1_var(int value) {
  tao_v1_var_struct_set_ix_v1_var(fortran_ptr_, value);
}
TaoVarProxyArray1D TaoV1VarProxy::v() const {
  return ProxyHelpers::get_type_array_1d<TaoVarProxyArray1D>(
      fortran_ptr_, tao_v1_var_struct_get_v_info);
}
double TaoGlobalProxy::beam_dead_cutoff() const {
  double value;
  tao_global_struct_get_beam_dead_cutoff(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_beam_dead_cutoff(double value) {
  tao_global_struct_set_beam_dead_cutoff(fortran_ptr_, value);
}
double TaoGlobalProxy::lm_opt_deriv_reinit() const {
  double value;
  tao_global_struct_get_lm_opt_deriv_reinit(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_lm_opt_deriv_reinit(double value) {
  tao_global_struct_set_lm_opt_deriv_reinit(fortran_ptr_, value);
}
double TaoGlobalProxy::de_lm_step_ratio() const {
  double value;
  tao_global_struct_get_de_lm_step_ratio(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_de_lm_step_ratio(double value) {
  tao_global_struct_set_de_lm_step_ratio(fortran_ptr_, value);
}
double TaoGlobalProxy::de_var_to_population_factor() const {
  double value;
  tao_global_struct_get_de_var_to_population_factor(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_de_var_to_population_factor(double value) {
  tao_global_struct_set_de_var_to_population_factor(fortran_ptr_, value);
}
double TaoGlobalProxy::lmdif_eps() const {
  double value;
  tao_global_struct_get_lmdif_eps(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_lmdif_eps(double value) {
  tao_global_struct_set_lmdif_eps(fortran_ptr_, value);
}
double TaoGlobalProxy::lmdif_negligible_merit() const {
  double value;
  tao_global_struct_get_lmdif_negligible_merit(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_lmdif_negligible_merit(double value) {
  tao_global_struct_set_lmdif_negligible_merit(fortran_ptr_, value);
}
double TaoGlobalProxy::svd_cutoff() const {
  double value;
  tao_global_struct_get_svd_cutoff(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_svd_cutoff(double value) {
  tao_global_struct_set_svd_cutoff(fortran_ptr_, value);
}
double TaoGlobalProxy::unstable_penalty() const {
  double value;
  tao_global_struct_get_unstable_penalty(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_unstable_penalty(double value) {
  tao_global_struct_set_unstable_penalty(fortran_ptr_, value);
}
double TaoGlobalProxy::merit_stop_value() const {
  double value;
  tao_global_struct_get_merit_stop_value(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_merit_stop_value(double value) {
  tao_global_struct_set_merit_stop_value(fortran_ptr_, value);
}
double TaoGlobalProxy::dmerit_stop_value() const {
  double value;
  tao_global_struct_get_dmerit_stop_value(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_dmerit_stop_value(double value) {
  tao_global_struct_set_dmerit_stop_value(fortran_ptr_, value);
}
double TaoGlobalProxy::random_sigma_cutoff() const {
  double value;
  tao_global_struct_get_random_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_random_sigma_cutoff(double value) {
  tao_global_struct_set_random_sigma_cutoff(fortran_ptr_, value);
}
double TaoGlobalProxy::delta_e_chrom() const {
  double value;
  tao_global_struct_get_delta_e_chrom(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_delta_e_chrom(double value) {
  tao_global_struct_set_delta_e_chrom(fortran_ptr_, value);
}
double TaoGlobalProxy::max_plot_time() const {
  double value;
  tao_global_struct_get_max_plot_time(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_max_plot_time(double value) {
  tao_global_struct_set_max_plot_time(fortran_ptr_, value);
}
int TaoGlobalProxy::default_universe() const {
  int value;
  tao_global_struct_get_default_universe(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_default_universe(int value) {
  tao_global_struct_set_default_universe(fortran_ptr_, value);
}
int TaoGlobalProxy::default_branch() const {
  int value;
  tao_global_struct_get_default_branch(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_default_branch(int value) {
  tao_global_struct_set_default_branch(fortran_ptr_, value);
}
int TaoGlobalProxy::n_opti_cycles() const {
  int value;
  tao_global_struct_get_n_opti_cycles(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_n_opti_cycles(int value) {
  tao_global_struct_set_n_opti_cycles(fortran_ptr_, value);
}
int TaoGlobalProxy::n_opti_loops() const {
  int value;
  tao_global_struct_get_n_opti_loops(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_n_opti_loops(int value) {
  tao_global_struct_set_n_opti_loops(fortran_ptr_, value);
}
int TaoGlobalProxy::n_threads() const {
  int value;
  tao_global_struct_get_n_threads(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_n_threads(int value) {
  tao_global_struct_set_n_threads(fortran_ptr_, value);
}
int TaoGlobalProxy::phase_units() const {
  int value;
  tao_global_struct_get_phase_units(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_phase_units(int value) {
  tao_global_struct_set_phase_units(fortran_ptr_, value);
}
int TaoGlobalProxy::bunch_to_plot() const {
  int value;
  tao_global_struct_get_bunch_to_plot(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_bunch_to_plot(int value) {
  tao_global_struct_set_bunch_to_plot(fortran_ptr_, value);
}
int TaoGlobalProxy::random_seed() const {
  int value;
  tao_global_struct_get_random_seed(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_random_seed(int value) {
  tao_global_struct_set_random_seed(fortran_ptr_, value);
}
int TaoGlobalProxy::n_top10_merit() const {
  int value;
  tao_global_struct_get_n_top10_merit(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_n_top10_merit(int value) {
  tao_global_struct_set_n_top10_merit(fortran_ptr_, value);
}
int TaoGlobalProxy::srdt_gen_n_slices() const {
  int value;
  tao_global_struct_get_srdt_gen_n_slices(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_srdt_gen_n_slices(int value) {
  tao_global_struct_set_srdt_gen_n_slices(fortran_ptr_, value);
}
int TaoGlobalProxy::datum_err_messages_max() const {
  int value;
  tao_global_struct_get_datum_err_messages_max(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_datum_err_messages_max(int value) {
  tao_global_struct_set_datum_err_messages_max(fortran_ptr_, value);
}
int TaoGlobalProxy::srdt_sxt_n_slices() const {
  int value;
  tao_global_struct_get_srdt_sxt_n_slices(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_srdt_sxt_n_slices(int value) {
  tao_global_struct_set_srdt_sxt_n_slices(fortran_ptr_, value);
}
bool TaoGlobalProxy::srdt_use_cache() const {
  bool value;
  tao_global_struct_get_srdt_use_cache(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_srdt_use_cache(bool value) {
  tao_global_struct_set_srdt_use_cache(fortran_ptr_, value);
}
std::string TaoGlobalProxy::quiet() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_quiet_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_quiet(const std::string& value) {
  tao_global_struct_set_quiet(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::random_engine() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_random_engine_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_random_engine(const std::string& value) {
  tao_global_struct_set_random_engine(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::random_gauss_converter() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_random_gauss_converter_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_random_gauss_converter(const std::string& value) {
  tao_global_struct_set_random_gauss_converter(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::track_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_track_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_track_type(const std::string& value) {
  tao_global_struct_set_track_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::lat_sigma_calc_uses_emit_from() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_lat_sigma_calc_uses_emit_from_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_lat_sigma_calc_uses_emit_from(
    const std::string& value) {
  tao_global_struct_set_lat_sigma_calc_uses_emit_from(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::prompt_string() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_prompt_string_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_prompt_string(const std::string& value) {
  tao_global_struct_set_prompt_string(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::prompt_color() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_prompt_color_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_prompt_color(const std::string& value) {
  tao_global_struct_set_prompt_color(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::optimizer() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_optimizer_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_optimizer(const std::string& value) {
  tao_global_struct_set_optimizer(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::print_command() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_print_command_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_print_command(const std::string& value) {
  tao_global_struct_set_print_command(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::var_out_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_var_out_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_var_out_file(const std::string& value) {
  tao_global_struct_set_var_out_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalProxy::history_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_global_struct_get_history_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalProxy::set_history_file(const std::string& value) {
  tao_global_struct_set_history_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoGlobalProxy::beam_timer_on() const {
  bool value;
  tao_global_struct_get_beam_timer_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_beam_timer_on(bool value) {
  tao_global_struct_set_beam_timer_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::box_plots() const {
  bool value;
  tao_global_struct_get_box_plots(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_box_plots(bool value) {
  tao_global_struct_set_box_plots(fortran_ptr_, value);
}
bool TaoGlobalProxy::blank_line_between_commands() const {
  bool value;
  tao_global_struct_get_blank_line_between_commands(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_blank_line_between_commands(bool value) {
  tao_global_struct_set_blank_line_between_commands(fortran_ptr_, value);
}
bool TaoGlobalProxy::cmd_file_abort_on_error() const {
  bool value;
  tao_global_struct_get_cmd_file_abort_on_error(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_cmd_file_abort_on_error(bool value) {
  tao_global_struct_set_cmd_file_abort_on_error(fortran_ptr_, value);
}
bool TaoGlobalProxy::concatenate_maps() const {
  bool value;
  tao_global_struct_get_concatenate_maps(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_concatenate_maps(bool value) {
  tao_global_struct_set_concatenate_maps(fortran_ptr_, value);
}
bool TaoGlobalProxy::derivative_recalc() const {
  bool value;
  tao_global_struct_get_derivative_recalc(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_derivative_recalc(bool value) {
  tao_global_struct_set_derivative_recalc(fortran_ptr_, value);
}
bool TaoGlobalProxy::derivative_uses_design() const {
  bool value;
  tao_global_struct_get_derivative_uses_design(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_derivative_uses_design(bool value) {
  tao_global_struct_set_derivative_uses_design(fortran_ptr_, value);
}
bool TaoGlobalProxy::disable_smooth_line_calc() const {
  bool value;
  tao_global_struct_get_disable_smooth_line_calc(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_disable_smooth_line_calc(bool value) {
  tao_global_struct_set_disable_smooth_line_calc(fortran_ptr_, value);
}
bool TaoGlobalProxy::draw_curve_off_scale_warn() const {
  bool value;
  tao_global_struct_get_draw_curve_off_scale_warn(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_draw_curve_off_scale_warn(bool value) {
  tao_global_struct_set_draw_curve_off_scale_warn(fortran_ptr_, value);
}
bool TaoGlobalProxy::external_plotting() const {
  bool value;
  tao_global_struct_get_external_plotting(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_external_plotting(bool value) {
  tao_global_struct_set_external_plotting(fortran_ptr_, value);
}
bool TaoGlobalProxy::label_lattice_elements() const {
  bool value;
  tao_global_struct_get_label_lattice_elements(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_label_lattice_elements(bool value) {
  tao_global_struct_set_label_lattice_elements(fortran_ptr_, value);
}
bool TaoGlobalProxy::label_keys() const {
  bool value;
  tao_global_struct_get_label_keys(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_label_keys(bool value) {
  tao_global_struct_set_label_keys(fortran_ptr_, value);
}
bool TaoGlobalProxy::lattice_calc_on() const {
  bool value;
  tao_global_struct_get_lattice_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_lattice_calc_on(bool value) {
  tao_global_struct_set_lattice_calc_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::only_limit_opt_vars() const {
  bool value;
  tao_global_struct_get_only_limit_opt_vars(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_only_limit_opt_vars(bool value) {
  tao_global_struct_set_only_limit_opt_vars(fortran_ptr_, value);
}
bool TaoGlobalProxy::opt_with_ref() const {
  bool value;
  tao_global_struct_get_opt_with_ref(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_opt_with_ref(bool value) {
  tao_global_struct_set_opt_with_ref(fortran_ptr_, value);
}
bool TaoGlobalProxy::opt_with_base() const {
  bool value;
  tao_global_struct_get_opt_with_base(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_opt_with_base(bool value) {
  tao_global_struct_set_opt_with_base(fortran_ptr_, value);
}
bool TaoGlobalProxy::opt_match_auto_recalc() const {
  bool value;
  tao_global_struct_get_opt_match_auto_recalc(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_opt_match_auto_recalc(bool value) {
  tao_global_struct_set_opt_match_auto_recalc(fortran_ptr_, value);
}
bool TaoGlobalProxy::opti_write_var_file() const {
  bool value;
  tao_global_struct_get_opti_write_var_file(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_opti_write_var_file(bool value) {
  tao_global_struct_set_opti_write_var_file(fortran_ptr_, value);
}
bool TaoGlobalProxy::optimizer_allow_user_abort() const {
  bool value;
  tao_global_struct_get_optimizer_allow_user_abort(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_optimizer_allow_user_abort(bool value) {
  tao_global_struct_set_optimizer_allow_user_abort(fortran_ptr_, value);
}
bool TaoGlobalProxy::optimizer_var_limit_warn() const {
  bool value;
  tao_global_struct_get_optimizer_var_limit_warn(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_optimizer_var_limit_warn(bool value) {
  tao_global_struct_set_optimizer_var_limit_warn(fortran_ptr_, value);
}
bool TaoGlobalProxy::plot_on() const {
  bool value;
  tao_global_struct_get_plot_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_plot_on(bool value) {
  tao_global_struct_set_plot_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::rad_int_user_calc_on() const {
  bool value;
  tao_global_struct_get_rad_int_user_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_rad_int_user_calc_on(bool value) {
  tao_global_struct_set_rad_int_user_calc_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::rf_on() const {
  bool value;
  tao_global_struct_get_rf_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_rf_on(bool value) {
  tao_global_struct_set_rf_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::single_step() const {
  bool value;
  tao_global_struct_get_single_step(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_single_step(bool value) {
  tao_global_struct_set_single_step(fortran_ptr_, value);
}
bool TaoGlobalProxy::stop_on_error() const {
  bool value;
  tao_global_struct_get_stop_on_error(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_stop_on_error(bool value) {
  tao_global_struct_set_stop_on_error(fortran_ptr_, value);
}
bool TaoGlobalProxy::svd_retreat_on_merit_increase() const {
  bool value;
  tao_global_struct_get_svd_retreat_on_merit_increase(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_svd_retreat_on_merit_increase(bool value) {
  tao_global_struct_set_svd_retreat_on_merit_increase(fortran_ptr_, value);
}
bool TaoGlobalProxy::var_limits_on() const {
  bool value;
  tao_global_struct_get_var_limits_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_var_limits_on(bool value) {
  tao_global_struct_set_var_limits_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::wait_for_CR_in_single_mode() const {
  bool value;
  tao_global_struct_get_wait_for_CR_in_single_mode(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_wait_for_CR_in_single_mode(bool value) {
  tao_global_struct_set_wait_for_CR_in_single_mode(fortran_ptr_, value);
}
bool TaoGlobalProxy::symbol_import() const {
  bool value;
  tao_global_struct_get_symbol_import(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_symbol_import(bool value) {
  tao_global_struct_set_symbol_import(fortran_ptr_, value);
}
bool TaoGlobalProxy::debug_on() const {
  bool value;
  tao_global_struct_get_debug_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_debug_on(bool value) {
  tao_global_struct_set_debug_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::expression_tree_on() const {
  bool value;
  tao_global_struct_get_expression_tree_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_expression_tree_on(bool value) {
  tao_global_struct_set_expression_tree_on(fortran_ptr_, value);
}
bool TaoGlobalProxy::verbose_on() const {
  bool value;
  tao_global_struct_get_verbose_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalProxy::set_verbose_on(bool value) {
  tao_global_struct_set_verbose_on(fortran_ptr_, value);
}
bool TaoInitProxy::parse_cmd_args() const {
  bool value;
  tao_init_struct_get_parse_cmd_args(fortran_ptr_, &value);
  return value;
}
void TaoInitProxy::set_parse_cmd_args(bool value) {
  tao_init_struct_set_parse_cmd_args(fortran_ptr_, value);
}
bool TaoInitProxy::debug_switch() const {
  bool value;
  tao_init_struct_get_debug_switch(fortran_ptr_, &value);
  return value;
}
void TaoInitProxy::set_debug_switch(bool value) {
  tao_init_struct_set_debug_switch(fortran_ptr_, value);
}
bool TaoInitProxy::external_plotting_switch() const {
  bool value;
  tao_init_struct_get_external_plotting_switch(fortran_ptr_, &value);
  return value;
}
void TaoInitProxy::set_external_plotting_switch(bool value) {
  tao_init_struct_set_external_plotting_switch(fortran_ptr_, value);
}
std::string TaoInitProxy::init_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_init_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_init_name(const std::string& value) {
  tao_init_struct_set_init_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_init_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_init_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_init_file(const std::string& value) {
  tao_init_struct_set_hook_init_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_lat_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_lat_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_lat_file(const std::string& value) {
  tao_init_struct_set_hook_lat_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_beam_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_beam_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_beam_file(const std::string& value) {
  tao_init_struct_set_hook_beam_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_data_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_data_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_data_file(const std::string& value) {
  tao_init_struct_set_hook_data_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_plot_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_plot_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_plot_file(const std::string& value) {
  tao_init_struct_set_hook_plot_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_startup_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_startup_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_startup_file(const std::string& value) {
  tao_init_struct_set_hook_startup_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_var_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_var_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_var_file(const std::string& value) {
  tao_init_struct_set_hook_var_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_building_wall_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_building_wall_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_building_wall_file(const std::string& value) {
  tao_init_struct_set_hook_building_wall_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::init_file_arg_path() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_init_file_arg_path_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_init_file_arg_path(const std::string& value) {
  tao_init_struct_set_init_file_arg_path(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::lattice_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_lattice_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_lattice_file_arg(const std::string& value) {
  tao_init_struct_set_lattice_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::hook_init_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_hook_init_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_hook_init_file_arg(const std::string& value) {
  tao_init_struct_set_hook_init_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::init_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_init_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_init_file_arg(const std::string& value) {
  tao_init_struct_set_init_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::beam_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_beam_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_beam_file_arg(const std::string& value) {
  tao_init_struct_set_beam_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::beam_init_position_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_beam_init_position_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_beam_init_position_file_arg(const std::string& value) {
  tao_init_struct_set_beam_init_position_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::command_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_command_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_command_arg(const std::string& value) {
  tao_init_struct_set_command_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::data_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_data_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_data_file_arg(const std::string& value) {
  tao_init_struct_set_data_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::plot_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_plot_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_plot_file_arg(const std::string& value) {
  tao_init_struct_set_plot_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::startup_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_startup_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_startup_file_arg(const std::string& value) {
  tao_init_struct_set_startup_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::var_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_var_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_var_file_arg(const std::string& value) {
  tao_init_struct_set_var_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::building_wall_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_building_wall_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_building_wall_file_arg(const std::string& value) {
  tao_init_struct_set_building_wall_file_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::geometry_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_geometry_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_geometry_arg(const std::string& value) {
  tao_init_struct_set_geometry_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::slice_lattice_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_slice_lattice_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_slice_lattice_arg(const std::string& value) {
  tao_init_struct_set_slice_lattice_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::start_branch_at_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_start_branch_at_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_start_branch_at_arg(const std::string& value) {
  tao_init_struct_set_start_branch_at_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::log_startup_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_log_startup_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_log_startup_arg(const std::string& value) {
  tao_init_struct_set_log_startup_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::no_stopping_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_no_stopping_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_no_stopping_arg(const std::string& value) {
  tao_init_struct_set_no_stopping_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::noplot_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_noplot_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_noplot_arg(const std::string& value) {
  tao_init_struct_set_noplot_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::no_rad_int_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_no_rad_int_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_no_rad_int_arg(const std::string& value) {
  tao_init_struct_set_no_rad_int_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::reverse_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_reverse_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_reverse_arg(const std::string& value) {
  tao_init_struct_set_reverse_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::debug_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_debug_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_debug_arg(const std::string& value) {
  tao_init_struct_set_debug_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::disable_smooth_line_calc_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_disable_smooth_line_calc_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_disable_smooth_line_calc_arg(const std::string& value) {
  tao_init_struct_set_disable_smooth_line_calc_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::rf_on_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_rf_on_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_rf_on_arg(const std::string& value) {
  tao_init_struct_set_rf_on_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::prompt_color_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_prompt_color_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_prompt_color_arg(const std::string& value) {
  tao_init_struct_set_prompt_color_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::quiet_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_quiet_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_quiet_arg(const std::string& value) {
  tao_init_struct_set_quiet_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::noinit_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_noinit_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_noinit_arg(const std::string& value) {
  tao_init_struct_set_noinit_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::nostartup_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_nostartup_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_nostartup_arg(const std::string& value) {
  tao_init_struct_set_nostartup_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::symbol_import_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_symbol_import_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_symbol_import_arg(const std::string& value) {
  tao_init_struct_set_symbol_import_arg(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitProxy::unique_name_suffix() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_init_struct_get_unique_name_suffix_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitProxy::set_unique_name_suffix(const std::string& value) {
  tao_init_struct_set_unique_name_suffix(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoPlotRegionProxyArray1D TaoCommonProxy::plot_place_buffer() const {
  return ProxyHelpers::get_type_array_1d<TaoPlotRegionProxyArray1D>(
      fortran_ptr_, tao_common_struct_get_plot_place_buffer_info);
}
FArray2D<double> TaoCommonProxy::covar() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_common_struct_get_covar_info);
}
FArray2D<double> TaoCommonProxy::alpha() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_common_struct_get_alpha_info);
}
double TaoCommonProxy::dummy_target() const {
  double value;
  tao_common_struct_get_dummy_target(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_dummy_target(double value) {
  tao_common_struct_set_dummy_target(fortran_ptr_, value);
}
int TaoCommonProxy::n_alias() const {
  int value;
  tao_common_struct_get_n_alias(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_n_alias(int value) {
  tao_common_struct_set_n_alias(fortran_ptr_, value);
}
int TaoCommonProxy::cmd_file_level() const {
  int value;
  tao_common_struct_get_cmd_file_level(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_cmd_file_level(int value) {
  tao_common_struct_set_cmd_file_level(fortran_ptr_, value);
}
int TaoCommonProxy::ix_key_bank() const {
  int value;
  tao_common_struct_get_ix_key_bank(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_ix_key_bank(int value) {
  tao_common_struct_set_ix_key_bank(fortran_ptr_, value);
}
int TaoCommonProxy::ix_history() const {
  int value;
  tao_common_struct_get_ix_history(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_ix_history(int value) {
  tao_common_struct_set_ix_history(fortran_ptr_, value);
}
int TaoCommonProxy::n_history() const {
  int value;
  tao_common_struct_get_n_history(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_n_history(int value) {
  tao_common_struct_set_n_history(fortran_ptr_, value);
}
int TaoCommonProxy::lev_loop() const {
  int value;
  tao_common_struct_get_lev_loop(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_lev_loop(int value) {
  tao_common_struct_set_lev_loop(fortran_ptr_, value);
}
int TaoCommonProxy::n_err_messages_printed() const {
  int value;
  tao_common_struct_get_n_err_messages_printed(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_n_err_messages_printed(int value) {
  tao_common_struct_set_n_err_messages_printed(fortran_ptr_, value);
}
int TaoCommonProxy::n_universes() const {
  int value;
  tao_common_struct_get_n_universes(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_n_universes(int value) {
  tao_common_struct_set_n_universes(fortran_ptr_, value);
}
int TaoCommonProxy::ix_beam_track_active_element() const {
  int value;
  tao_common_struct_get_ix_beam_track_active_element(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_ix_beam_track_active_element(int value) {
  tao_common_struct_set_ix_beam_track_active_element(fortran_ptr_, value);
}
bool TaoCommonProxy::cmd_file_paused() const {
  bool value;
  tao_common_struct_get_cmd_file_paused(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_cmd_file_paused(bool value) {
  tao_common_struct_set_cmd_file_paused(fortran_ptr_, value);
}
bool TaoCommonProxy::use_cmd_here() const {
  bool value;
  tao_common_struct_get_use_cmd_here(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_use_cmd_here(bool value) {
  tao_common_struct_set_use_cmd_here(fortran_ptr_, value);
}
bool TaoCommonProxy::cmd_from_cmd_file() const {
  bool value;
  tao_common_struct_get_cmd_from_cmd_file(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_cmd_from_cmd_file(bool value) {
  tao_common_struct_set_cmd_from_cmd_file(fortran_ptr_, value);
}
bool TaoCommonProxy::use_saved_beam_in_tracking() const {
  bool value;
  tao_common_struct_get_use_saved_beam_in_tracking(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_use_saved_beam_in_tracking(bool value) {
  tao_common_struct_set_use_saved_beam_in_tracking(fortran_ptr_, value);
}
bool TaoCommonProxy::single_mode() const {
  bool value;
  tao_common_struct_get_single_mode(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_single_mode(bool value) {
  tao_common_struct_set_single_mode(fortran_ptr_, value);
}
bool TaoCommonProxy::combine_consecutive_elements_of_like_name() const {
  bool value;
  tao_common_struct_get_combine_consecutive_elements_of_like_name(
      fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_combine_consecutive_elements_of_like_name(bool value) {
  tao_common_struct_set_combine_consecutive_elements_of_like_name(
      fortran_ptr_, value);
}
bool TaoCommonProxy::have_tracked_beam() const {
  bool value;
  tao_common_struct_get_have_tracked_beam(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_have_tracked_beam(bool value) {
  tao_common_struct_set_have_tracked_beam(fortran_ptr_, value);
}
bool TaoCommonProxy::init_plot_needed() const {
  bool value;
  tao_common_struct_get_init_plot_needed(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_init_plot_needed(bool value) {
  tao_common_struct_set_init_plot_needed(fortran_ptr_, value);
}
bool TaoCommonProxy::init_beam() const {
  bool value;
  tao_common_struct_get_init_beam(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_init_beam(bool value) {
  tao_common_struct_set_init_beam(fortran_ptr_, value);
}
bool TaoCommonProxy::init_var() const {
  bool value;
  tao_common_struct_get_init_var(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_init_var(bool value) {
  tao_common_struct_set_init_var(fortran_ptr_, value);
}
bool TaoCommonProxy::init_read_lat_info() const {
  bool value;
  tao_common_struct_get_init_read_lat_info(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_init_read_lat_info(bool value) {
  tao_common_struct_set_init_read_lat_info(fortran_ptr_, value);
}
bool TaoCommonProxy::optimizer_running() const {
  bool value;
  tao_common_struct_get_optimizer_running(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_optimizer_running(bool value) {
  tao_common_struct_set_optimizer_running(fortran_ptr_, value);
}
bool TaoCommonProxy::have_datums_using_expressions() const {
  bool value;
  tao_common_struct_get_have_datums_using_expressions(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_have_datums_using_expressions(bool value) {
  tao_common_struct_set_have_datums_using_expressions(fortran_ptr_, value);
}
bool TaoCommonProxy::print_to_terminal() const {
  bool value;
  tao_common_struct_get_print_to_terminal(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_print_to_terminal(bool value) {
  tao_common_struct_set_print_to_terminal(fortran_ptr_, value);
}
bool TaoCommonProxy::lattice_calc_done() const {
  bool value;
  tao_common_struct_get_lattice_calc_done(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_lattice_calc_done(bool value) {
  tao_common_struct_set_lattice_calc_done(fortran_ptr_, value);
}
bool TaoCommonProxy::add_measurement_noise() const {
  bool value;
  tao_common_struct_get_add_measurement_noise(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_add_measurement_noise(bool value) {
  tao_common_struct_set_add_measurement_noise(fortran_ptr_, value);
}
bool TaoCommonProxy::command_arg_has_been_executed() const {
  bool value;
  tao_common_struct_get_command_arg_has_been_executed(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_command_arg_has_been_executed(bool value) {
  tao_common_struct_set_command_arg_has_been_executed(fortran_ptr_, value);
}
bool TaoCommonProxy::all_merit_weights_positive() const {
  bool value;
  tao_common_struct_get_all_merit_weights_positive(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_all_merit_weights_positive(bool value) {
  tao_common_struct_set_all_merit_weights_positive(fortran_ptr_, value);
}
bool TaoCommonProxy::multi_turn_orbit_is_plotted() const {
  bool value;
  tao_common_struct_get_multi_turn_orbit_is_plotted(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_multi_turn_orbit_is_plotted(bool value) {
  tao_common_struct_set_multi_turn_orbit_is_plotted(fortran_ptr_, value);
}
bool TaoCommonProxy::force_chrom_calc() const {
  bool value;
  tao_common_struct_get_force_chrom_calc(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_force_chrom_calc(bool value) {
  tao_common_struct_set_force_chrom_calc(fortran_ptr_, value);
}
bool TaoCommonProxy::force_rad_int_calc() const {
  bool value;
  tao_common_struct_get_force_rad_int_calc(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_force_rad_int_calc(bool value) {
  tao_common_struct_set_force_rad_int_calc(fortran_ptr_, value);
}
bool TaoCommonProxy::rad_int_ri_calc_on() const {
  bool value;
  tao_common_struct_get_rad_int_ri_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_rad_int_ri_calc_on(bool value) {
  tao_common_struct_set_rad_int_ri_calc_on(fortran_ptr_, value);
}
bool TaoCommonProxy::rad_int_6d_calc_on() const {
  bool value;
  tao_common_struct_get_rad_int_6d_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoCommonProxy::set_rad_int_6d_calc_on(bool value) {
  tao_common_struct_set_rad_int_6d_calc_on(fortran_ptr_, value);
}
FCharArray1D TaoCommonProxy::valid_plot_who() const {
  return ProxyHelpers::get_char_array_1d(
      fortran_ptr_, tao_common_struct_get_valid_plot_who_info);
}
std::string TaoCommonProxy::single_mode_buffer() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_common_struct_get_single_mode_buffer_info);
  return std::string(arr.data(), arr.size());
}
void TaoCommonProxy::set_single_mode_buffer(const std::string& value) {
  tao_common_struct_set_single_mode_buffer(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCommonProxy::cmd() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_common_struct_get_cmd_info);
  return std::string(arr.data(), arr.size());
}
void TaoCommonProxy::set_cmd(const std::string& value) {
  tao_common_struct_set_cmd(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCommonProxy::saved_cmd_line() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_common_struct_get_saved_cmd_line_info);
  return std::string(arr.data(), arr.size());
}
void TaoCommonProxy::set_saved_cmd_line(const std::string& value) {
  tao_common_struct_set_saved_cmd_line(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoTitleProxy TaoPlotPageProxy::title() const {
  void* ptr;
  tao_plot_page_struct_get_title(fortran_ptr_, &ptr);
  return TaoTitleProxy(ptr);
}
void TaoPlotPageProxy::set_title(const TaoTitleProxy& src) {
  tao_plot_page_struct_set_title(fortran_ptr_, src.get_fortran_ptr());
}
TaoTitleProxy TaoPlotPageProxy::subtitle() const {
  void* ptr;
  tao_plot_page_struct_get_subtitle(fortran_ptr_, &ptr);
  return TaoTitleProxy(ptr);
}
void TaoPlotPageProxy::set_subtitle(const TaoTitleProxy& src) {
  tao_plot_page_struct_set_subtitle(fortran_ptr_, src.get_fortran_ptr());
}
QpRectProxy TaoPlotPageProxy::border() const {
  void* ptr;
  tao_plot_page_struct_get_border(fortran_ptr_, &ptr);
  return QpRectProxy(ptr);
}
void TaoPlotPageProxy::set_border(const QpRectProxy& src) {
  tao_plot_page_struct_set_border(fortran_ptr_, src.get_fortran_ptr());
}
TaoDrawingProxy TaoPlotPageProxy::floor_plan() const {
  void* ptr;
  tao_plot_page_struct_get_floor_plan(fortran_ptr_, &ptr);
  return TaoDrawingProxy(ptr);
}
void TaoPlotPageProxy::set_floor_plan(const TaoDrawingProxy& src) {
  tao_plot_page_struct_set_floor_plan(fortran_ptr_, src.get_fortran_ptr());
}
TaoDrawingProxy TaoPlotPageProxy::lat_layout() const {
  void* ptr;
  tao_plot_page_struct_get_lat_layout(fortran_ptr_, &ptr);
  return TaoDrawingProxy(ptr);
}
void TaoPlotPageProxy::set_lat_layout(const TaoDrawingProxy& src) {
  tao_plot_page_struct_set_lat_layout(fortran_ptr_, src.get_fortran_ptr());
}
TaoShapePatternProxyArray1D TaoPlotPageProxy::pattern() const {
  return ProxyHelpers::get_type_array_1d<TaoShapePatternProxyArray1D>(
      fortran_ptr_, tao_plot_page_struct_get_pattern_info);
}
TaoPlotProxyArray1D TaoPlotPageProxy::template_() const {
  return ProxyHelpers::get_type_array_1d<TaoPlotProxyArray1D>(
      fortran_ptr_, tao_plot_page_struct_get_template_info);
}
TaoPlotRegionProxyArray1D TaoPlotPageProxy::region() const {
  return ProxyHelpers::get_type_array_1d<TaoPlotRegionProxyArray1D>(
      fortran_ptr_, tao_plot_page_struct_get_region_info);
}
std::string TaoPlotPageProxy::plot_display_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_plot_page_struct_get_plot_display_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotPageProxy::set_plot_display_type(const std::string& value) {
  tao_plot_page_struct_set_plot_display_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FArray1D<double> TaoPlotPageProxy::size() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_plot_page_struct_get_size_info);
}
double TaoPlotPageProxy::text_height() const {
  double value;
  tao_plot_page_struct_get_text_height(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_text_height(double value) {
  tao_plot_page_struct_set_text_height(fortran_ptr_, value);
}
double TaoPlotPageProxy::main_title_text_scale() const {
  double value;
  tao_plot_page_struct_get_main_title_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_main_title_text_scale(double value) {
  tao_plot_page_struct_set_main_title_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::graph_title_text_scale() const {
  double value;
  tao_plot_page_struct_get_graph_title_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_graph_title_text_scale(double value) {
  tao_plot_page_struct_set_graph_title_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::axis_number_text_scale() const {
  double value;
  tao_plot_page_struct_get_axis_number_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_axis_number_text_scale(double value) {
  tao_plot_page_struct_set_axis_number_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::axis_label_text_scale() const {
  double value;
  tao_plot_page_struct_get_axis_label_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_axis_label_text_scale(double value) {
  tao_plot_page_struct_set_axis_label_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::legend_text_scale() const {
  double value;
  tao_plot_page_struct_get_legend_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_legend_text_scale(double value) {
  tao_plot_page_struct_set_legend_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::key_table_text_scale() const {
  double value;
  tao_plot_page_struct_get_key_table_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_key_table_text_scale(double value) {
  tao_plot_page_struct_set_key_table_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::floor_plan_shape_scale() const {
  double value;
  tao_plot_page_struct_get_floor_plan_shape_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_floor_plan_shape_scale(double value) {
  tao_plot_page_struct_set_floor_plan_shape_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::floor_plan_text_scale() const {
  double value;
  tao_plot_page_struct_get_floor_plan_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_floor_plan_text_scale(double value) {
  tao_plot_page_struct_set_floor_plan_text_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::lat_layout_shape_scale() const {
  double value;
  tao_plot_page_struct_get_lat_layout_shape_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_lat_layout_shape_scale(double value) {
  tao_plot_page_struct_set_lat_layout_shape_scale(fortran_ptr_, value);
}
double TaoPlotPageProxy::lat_layout_text_scale() const {
  double value;
  tao_plot_page_struct_get_lat_layout_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_lat_layout_text_scale(double value) {
  tao_plot_page_struct_set_lat_layout_text_scale(fortran_ptr_, value);
}
int TaoPlotPageProxy::n_curve_pts() const {
  int value;
  tao_plot_page_struct_get_n_curve_pts(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_n_curve_pts(int value) {
  tao_plot_page_struct_set_n_curve_pts(fortran_ptr_, value);
}
int TaoPlotPageProxy::id_window() const {
  int value;
  tao_plot_page_struct_get_id_window(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_id_window(int value) {
  tao_plot_page_struct_set_id_window(fortran_ptr_, value);
}
bool TaoPlotPageProxy::delete_overlapping_plots() const {
  bool value;
  tao_plot_page_struct_get_delete_overlapping_plots(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_delete_overlapping_plots(bool value) {
  tao_plot_page_struct_set_delete_overlapping_plots(fortran_ptr_, value);
}
bool TaoPlotPageProxy::draw_graph_title_suffix() const {
  bool value;
  tao_plot_page_struct_get_draw_graph_title_suffix(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageProxy::set_draw_graph_title_suffix(bool value) {
  tao_plot_page_struct_set_draw_graph_title_suffix(fortran_ptr_, value);
}
TaoBuildingWallOrientationProxy TaoBuildingWallProxy::orientation() const {
  void* ptr;
  tao_building_wall_struct_get_orientation(fortran_ptr_, &ptr);
  return TaoBuildingWallOrientationProxy(ptr);
}
void TaoBuildingWallProxy::set_orientation(
    const TaoBuildingWallOrientationProxy& src) {
  tao_building_wall_struct_set_orientation(fortran_ptr_, src.get_fortran_ptr());
}
TaoBuildingWallSectionProxyArray1D TaoBuildingWallProxy::section() const {
  return ProxyHelpers::get_type_array_1d<TaoBuildingWallSectionProxyArray1D>(
      fortran_ptr_, tao_building_wall_struct_get_section_info);
}
double TaoBuildingWallOrientationProxy::theta() const {
  double value;
  tao_building_wall_orientation_struct_get_theta(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallOrientationProxy::set_theta(double value) {
  tao_building_wall_orientation_struct_set_theta(fortran_ptr_, value);
}
double TaoBuildingWallOrientationProxy::x_offset() const {
  double value;
  tao_building_wall_orientation_struct_get_x_offset(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallOrientationProxy::set_x_offset(double value) {
  tao_building_wall_orientation_struct_set_x_offset(fortran_ptr_, value);
}
double TaoBuildingWallOrientationProxy::z_offset() const {
  double value;
  tao_building_wall_orientation_struct_get_z_offset(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallOrientationProxy::set_z_offset(double value) {
  tao_building_wall_orientation_struct_set_z_offset(fortran_ptr_, value);
}
std::string TaoBuildingWallSectionProxy::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_building_wall_section_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoBuildingWallSectionProxy::set_name(const std::string& value) {
  tao_building_wall_section_struct_set_name(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoBuildingWallSectionProxy::constraint() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_building_wall_section_struct_get_constraint_info);
  return std::string(arr.data(), arr.size());
}
void TaoBuildingWallSectionProxy::set_constraint(const std::string& value) {
  tao_building_wall_section_struct_set_constraint(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoBuildingWallPointProxyArray1D TaoBuildingWallSectionProxy::point() const {
  return ProxyHelpers::get_type_array_1d<TaoBuildingWallPointProxyArray1D>(
      fortran_ptr_, tao_building_wall_section_struct_get_point_info);
}
double TaoBuildingWallPointProxy::z() const {
  double value;
  tao_building_wall_point_struct_get_z(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointProxy::set_z(double value) {
  tao_building_wall_point_struct_set_z(fortran_ptr_, value);
}
double TaoBuildingWallPointProxy::x() const {
  double value;
  tao_building_wall_point_struct_get_x(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointProxy::set_x(double value) {
  tao_building_wall_point_struct_set_x(fortran_ptr_, value);
}
double TaoBuildingWallPointProxy::radius() const {
  double value;
  tao_building_wall_point_struct_get_radius(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointProxy::set_radius(double value) {
  tao_building_wall_point_struct_set_radius(fortran_ptr_, value);
}
double TaoBuildingWallPointProxy::z_center() const {
  double value;
  tao_building_wall_point_struct_get_z_center(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointProxy::set_z_center(double value) {
  tao_building_wall_point_struct_set_z_center(fortran_ptr_, value);
}
double TaoBuildingWallPointProxy::x_center() const {
  double value;
  tao_building_wall_point_struct_get_x_center(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointProxy::set_x_center(double value) {
  tao_building_wall_point_struct_set_x_center(fortran_ptr_, value);
}
std::string TaoWaveProxy::data_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, tao_wave_struct_get_data_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoWaveProxy::set_data_type(const std::string& value) {
  tao_wave_struct_set_data_type(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoWaveProxy::rms_rel_a() const {
  double value;
  tao_wave_struct_get_rms_rel_a(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_a(double value) {
  tao_wave_struct_set_rms_rel_a(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_b() const {
  double value;
  tao_wave_struct_get_rms_rel_b(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_b(double value) {
  tao_wave_struct_set_rms_rel_b(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_as() const {
  double value;
  tao_wave_struct_get_rms_rel_as(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_as(double value) {
  tao_wave_struct_set_rms_rel_as(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_bs() const {
  double value;
  tao_wave_struct_get_rms_rel_bs(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_bs(double value) {
  tao_wave_struct_set_rms_rel_bs(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_ar() const {
  double value;
  tao_wave_struct_get_rms_rel_ar(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_ar(double value) {
  tao_wave_struct_set_rms_rel_ar(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_br() const {
  double value;
  tao_wave_struct_get_rms_rel_br(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_br(double value) {
  tao_wave_struct_set_rms_rel_br(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_k() const {
  double value;
  tao_wave_struct_get_rms_rel_k(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_k(double value) {
  tao_wave_struct_set_rms_rel_k(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_ks() const {
  double value;
  tao_wave_struct_get_rms_rel_ks(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_ks(double value) {
  tao_wave_struct_set_rms_rel_ks(fortran_ptr_, value);
}
double TaoWaveProxy::rms_rel_kr() const {
  double value;
  tao_wave_struct_get_rms_rel_kr(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_rel_kr(double value) {
  tao_wave_struct_set_rms_rel_kr(fortran_ptr_, value);
}
double TaoWaveProxy::rms_phi() const {
  double value;
  tao_wave_struct_get_rms_phi(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_phi(double value) {
  tao_wave_struct_set_rms_phi(fortran_ptr_, value);
}
double TaoWaveProxy::rms_phi_s() const {
  double value;
  tao_wave_struct_get_rms_phi_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_phi_s(double value) {
  tao_wave_struct_set_rms_phi_s(fortran_ptr_, value);
}
double TaoWaveProxy::rms_phi_r() const {
  double value;
  tao_wave_struct_get_rms_phi_r(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_rms_phi_r(double value) {
  tao_wave_struct_set_rms_phi_r(fortran_ptr_, value);
}
double TaoWaveProxy::amp_ba_s() const {
  double value;
  tao_wave_struct_get_amp_ba_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_amp_ba_s(double value) {
  tao_wave_struct_set_amp_ba_s(fortran_ptr_, value);
}
double TaoWaveProxy::amp_ba_r() const {
  double value;
  tao_wave_struct_get_amp_ba_r(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_amp_ba_r(double value) {
  tao_wave_struct_set_amp_ba_r(fortran_ptr_, value);
}
double TaoWaveProxy::chi_a() const {
  double value;
  tao_wave_struct_get_chi_a(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_chi_a(double value) {
  tao_wave_struct_set_chi_a(fortran_ptr_, value);
}
double TaoWaveProxy::chi_c() const {
  double value;
  tao_wave_struct_get_chi_c(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_chi_c(double value) {
  tao_wave_struct_set_chi_c(fortran_ptr_, value);
}
double TaoWaveProxy::chi_ba() const {
  double value;
  tao_wave_struct_get_chi_ba(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_chi_ba(double value) {
  tao_wave_struct_set_chi_ba(fortran_ptr_, value);
}
FArray1D<double> TaoWaveProxy::amp_a() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_wave_struct_get_amp_a_info);
}
FArray1D<double> TaoWaveProxy::amp_b() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_wave_struct_get_amp_b_info);
}
FArray1D<double> TaoWaveProxy::amp_ba() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_wave_struct_get_amp_ba_info);
}
FArray1D<double> TaoWaveProxy::coef_a() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_wave_struct_get_coef_a_info);
}
FArray1D<double> TaoWaveProxy::coef_b() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_wave_struct_get_coef_b_info);
}
FArray1D<double> TaoWaveProxy::coef_ba() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, tao_wave_struct_get_coef_ba_info);
}
int TaoWaveProxy::n_func() const {
  int value;
  tao_wave_struct_get_n_func(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_n_func(int value) {
  tao_wave_struct_set_n_func(fortran_ptr_, value);
}
int TaoWaveProxy::ix_a1() const {
  int value;
  tao_wave_struct_get_ix_a1(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_ix_a1(int value) {
  tao_wave_struct_set_ix_a1(fortran_ptr_, value);
}
int TaoWaveProxy::ix_a2() const {
  int value;
  tao_wave_struct_get_ix_a2(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_ix_a2(int value) {
  tao_wave_struct_set_ix_a2(fortran_ptr_, value);
}
int TaoWaveProxy::ix_b1() const {
  int value;
  tao_wave_struct_get_ix_b1(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_ix_b1(int value) {
  tao_wave_struct_set_ix_b1(fortran_ptr_, value);
}
int TaoWaveProxy::ix_b2() const {
  int value;
  tao_wave_struct_get_ix_b2(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_ix_b2(int value) {
  tao_wave_struct_set_ix_b2(fortran_ptr_, value);
}
int TaoWaveProxy::i_a1() const {
  int value;
  tao_wave_struct_get_i_a1(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_i_a1(int value) {
  tao_wave_struct_set_i_a1(fortran_ptr_, value);
}
int TaoWaveProxy::i_a2() const {
  int value;
  tao_wave_struct_get_i_a2(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_i_a2(int value) {
  tao_wave_struct_set_i_a2(fortran_ptr_, value);
}
int TaoWaveProxy::i_b1() const {
  int value;
  tao_wave_struct_get_i_b1(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_i_b1(int value) {
  tao_wave_struct_set_i_b1(fortran_ptr_, value);
}
int TaoWaveProxy::i_b2() const {
  int value;
  tao_wave_struct_get_i_b2(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_i_b2(int value) {
  tao_wave_struct_set_i_b2(fortran_ptr_, value);
}
int TaoWaveProxy::n_a() const {
  int value;
  tao_wave_struct_get_n_a(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_n_a(int value) {
  tao_wave_struct_set_n_a(fortran_ptr_, value);
}
int TaoWaveProxy::n_b() const {
  int value;
  tao_wave_struct_get_n_b(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_n_b(int value) {
  tao_wave_struct_set_n_b(fortran_ptr_, value);
}
int TaoWaveProxy::i_curve_wrap_pt() const {
  int value;
  tao_wave_struct_get_i_curve_wrap_pt(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_i_curve_wrap_pt(int value) {
  tao_wave_struct_set_i_curve_wrap_pt(fortran_ptr_, value);
}
FArray1D<int> TaoWaveProxy::ix_data() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, tao_wave_struct_get_ix_data_info);
}
int TaoWaveProxy::n_kick() const {
  int value;
  tao_wave_struct_get_n_kick(fortran_ptr_, &value);
  return value;
}
void TaoWaveProxy::set_n_kick(int value) {
  tao_wave_struct_set_n_kick(fortran_ptr_, value);
}
TaoWaveKickPtProxyArray1D TaoWaveProxy::kick() const {
  return ProxyHelpers::get_type_array_1d<TaoWaveKickPtProxyArray1D>(
      fortran_ptr_, tao_wave_struct_get_kick_info);
}
TaoGraphProxy TaoWaveProxy::base_graph() const {
  void* ptr;
  tao_wave_struct_get_base_graph(fortran_ptr_, &ptr);
  return TaoGraphProxy(ptr);
}
void TaoWaveProxy::set_base_graph(const TaoGraphProxy& src) {
  tao_wave_struct_set_base_graph(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoPlotRegionProxy> TaoWaveProxy::region() const {
  void* ptr;
  tao_wave_struct_get_region(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoPlotRegionProxy(ptr);
}
void TaoWaveProxy::set_region(const TaoPlotRegionProxy& src) {
  tao_wave_struct_set_region(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoD1DataProxy> TaoWaveProxy::d1_dat() const {
  void* ptr;
  tao_wave_struct_get_d1_dat(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoD1DataProxy(ptr);
}
void TaoWaveProxy::set_d1_dat(const TaoD1DataProxy& src) {
  tao_wave_struct_set_d1_dat(fortran_ptr_, src.get_fortran_ptr());
}
double TaoWaveKickPtProxy::phi_s() const {
  double value;
  tao_wave_kick_pt_struct_get_phi_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtProxy::set_phi_s(double value) {
  tao_wave_kick_pt_struct_set_phi_s(fortran_ptr_, value);
}
double TaoWaveKickPtProxy::phi_r() const {
  double value;
  tao_wave_kick_pt_struct_get_phi_r(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtProxy::set_phi_r(double value) {
  tao_wave_kick_pt_struct_set_phi_r(fortran_ptr_, value);
}
double TaoWaveKickPtProxy::phi() const {
  double value;
  tao_wave_kick_pt_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtProxy::set_phi(double value) {
  tao_wave_kick_pt_struct_set_phi(fortran_ptr_, value);
}
double TaoWaveKickPtProxy::amp() const {
  double value;
  tao_wave_kick_pt_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtProxy::set_amp(double value) {
  tao_wave_kick_pt_struct_set_amp(fortran_ptr_, value);
}
double TaoWaveKickPtProxy::s() const {
  double value;
  tao_wave_kick_pt_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtProxy::set_s(double value) {
  tao_wave_kick_pt_struct_set_s(fortran_ptr_, value);
}
int TaoWaveKickPtProxy::ix_dat_before_kick() const {
  int value;
  tao_wave_kick_pt_struct_get_ix_dat_before_kick(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtProxy::set_ix_dat_before_kick(int value) {
  tao_wave_kick_pt_struct_set_ix_dat_before_kick(fortran_ptr_, value);
}
std::optional<EleProxy> TaoWaveKickPtProxy::ele() const {
  void* ptr;
  tao_wave_kick_pt_struct_get_ele(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleProxy(ptr);
}
void TaoWaveKickPtProxy::set_ele(const EleProxy& src) {
  tao_wave_kick_pt_struct_set_ele(fortran_ptr_, src.get_fortran_ptr());
}
std::string TaoCmdHistoryProxy::cmd() const {
  return ProxyHelpers::get_string(
      fortran_ptr_, tao_cmd_history_struct_get_cmd_info);
}
void TaoCmdHistoryProxy::set_cmd(const std::string& value) {
  tao_cmd_history_struct_set_cmd(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoCmdHistoryProxy::ix() const {
  int value;
  tao_cmd_history_struct_get_ix(fortran_ptr_, &value);
  return value;
}
void TaoCmdHistoryProxy::set_ix(int value) {
  tao_cmd_history_struct_set_ix(fortran_ptr_, value);
}
std::optional<TaoLatticeProxy> TaoUniverseProxy::model() const {
  void* ptr;
  tao_universe_struct_get_model(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeProxy(ptr);
}
void TaoUniverseProxy::set_model(const TaoLatticeProxy& src) {
  tao_universe_struct_set_model(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoLatticeProxy> TaoUniverseProxy::design() const {
  void* ptr;
  tao_universe_struct_get_design(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeProxy(ptr);
}
void TaoUniverseProxy::set_design(const TaoLatticeProxy& src) {
  tao_universe_struct_set_design(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoLatticeProxy> TaoUniverseProxy::base() const {
  void* ptr;
  tao_universe_struct_get_base(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeProxy(ptr);
}
void TaoUniverseProxy::set_base(const TaoLatticeProxy& src) {
  tao_universe_struct_set_base(fortran_ptr_, src.get_fortran_ptr());
}
TaoBeamUniProxy TaoUniverseProxy::beam() const {
  void* ptr;
  tao_universe_struct_get_beam(fortran_ptr_, &ptr);
  return TaoBeamUniProxy(ptr);
}
void TaoUniverseProxy::set_beam(const TaoBeamUniProxy& src) {
  tao_universe_struct_set_beam(fortran_ptr_, src.get_fortran_ptr());
}
TaoDynamicApertureProxy TaoUniverseProxy::dynamic_aperture() const {
  void* ptr;
  tao_universe_struct_get_dynamic_aperture(fortran_ptr_, &ptr);
  return TaoDynamicApertureProxy(ptr);
}
void TaoUniverseProxy::set_dynamic_aperture(
    const TaoDynamicApertureProxy& src) {
  tao_universe_struct_set_dynamic_aperture(fortran_ptr_, src.get_fortran_ptr());
}
TaoModelBranchProxyArray1D TaoUniverseProxy::model_branch() const {
  return ProxyHelpers::get_type_array_1d<TaoModelBranchProxyArray1D>(
      fortran_ptr_, tao_universe_struct_get_model_branch_info);
}
TaoD2DataProxyArray1D TaoUniverseProxy::d2_data() const {
  return ProxyHelpers::get_type_array_1d<TaoD2DataProxyArray1D>(
      fortran_ptr_, tao_universe_struct_get_d2_data_info);
}
TaoDataProxyArray1D TaoUniverseProxy::data() const {
  return ProxyHelpers::get_type_array_1d<TaoDataProxyArray1D>(
      fortran_ptr_, tao_universe_struct_get_data_info);
}
TaoPingScaleProxy TaoUniverseProxy::ping_scale() const {
  void* ptr;
  tao_universe_struct_get_ping_scale(fortran_ptr_, &ptr);
  return TaoPingScaleProxy(ptr);
}
void TaoUniverseProxy::set_ping_scale(const TaoPingScaleProxy& src) {
  tao_universe_struct_set_ping_scale(fortran_ptr_, src.get_fortran_ptr());
}
LatProxy TaoUniverseProxy::scratch_lat() const {
  void* ptr;
  tao_universe_struct_get_scratch_lat(fortran_ptr_, &ptr);
  return LatProxy(ptr);
}
void TaoUniverseProxy::set_scratch_lat(const LatProxy& src) {
  tao_universe_struct_set_scratch_lat(fortran_ptr_, src.get_fortran_ptr());
}
TaoUniverseCalcProxy TaoUniverseProxy::calc() const {
  void* ptr;
  tao_universe_struct_get_calc(fortran_ptr_, &ptr);
  return TaoUniverseCalcProxy(ptr);
}
void TaoUniverseProxy::set_calc(const TaoUniverseCalcProxy& src) {
  tao_universe_struct_set_calc(fortran_ptr_, src.get_fortran_ptr());
}
LatEleOrderProxy TaoUniverseProxy::ele_order() const {
  void* ptr;
  tao_universe_struct_get_ele_order(fortran_ptr_, &ptr);
  return LatEleOrderProxy(ptr);
}
void TaoUniverseProxy::set_ele_order(const LatEleOrderProxy& src) {
  tao_universe_struct_set_ele_order(fortran_ptr_, src.get_fortran_ptr());
}
TaoSpinMapProxy TaoUniverseProxy::spin_map() const {
  void* ptr;
  tao_universe_struct_get_spin_map(fortran_ptr_, &ptr);
  return TaoSpinMapProxy(ptr);
}
void TaoUniverseProxy::set_spin_map(const TaoSpinMapProxy& src) {
  tao_universe_struct_set_spin_map(fortran_ptr_, src.get_fortran_ptr());
}
FArray2D<double> TaoUniverseProxy::dModel_dVar() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, tao_universe_struct_get_dModel_dVar_info);
}
int TaoUniverseProxy::ix_uni() const {
  int value;
  tao_universe_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoUniverseProxy::set_ix_uni(int value) {
  tao_universe_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoUniverseProxy::n_d2_data_used() const {
  int value;
  tao_universe_struct_get_n_d2_data_used(fortran_ptr_, &value);
  return value;
}
void TaoUniverseProxy::set_n_d2_data_used(int value) {
  tao_universe_struct_set_n_d2_data_used(fortran_ptr_, value);
}
int TaoUniverseProxy::n_data_used() const {
  int value;
  tao_universe_struct_get_n_data_used(fortran_ptr_, &value);
  return value;
}
void TaoUniverseProxy::set_n_data_used(int value) {
  tao_universe_struct_set_n_data_used(fortran_ptr_, value);
}
bool TaoUniverseProxy::is_on() const {
  bool value;
  tao_universe_struct_get_is_on(fortran_ptr_, &value);
  return value;
}
void TaoUniverseProxy::set_is_on(bool value) {
  tao_universe_struct_set_is_on(fortran_ptr_, value);
}
bool TaoUniverseProxy::design_same_as_previous() const {
  bool value;
  tao_universe_struct_get_design_same_as_previous(fortran_ptr_, &value);
  return value;
}
void TaoUniverseProxy::set_design_same_as_previous(bool value) {
  tao_universe_struct_set_design_same_as_previous(fortran_ptr_, value);
}
bool TaoUniverseProxy::picked_uni() const {
  bool value;
  tao_universe_struct_get_picked_uni(fortran_ptr_, &value);
  return value;
}
void TaoUniverseProxy::set_picked_uni(bool value) {
  tao_universe_struct_set_picked_uni(fortran_ptr_, value);
}
double MadEnergyProxy::total() const {
  double value;
  mad_energy_struct_get_total(fortran_ptr_, &value);
  return value;
}
void MadEnergyProxy::set_total(double value) {
  mad_energy_struct_set_total(fortran_ptr_, value);
}
double MadEnergyProxy::beta() const {
  double value;
  mad_energy_struct_get_beta(fortran_ptr_, &value);
  return value;
}
void MadEnergyProxy::set_beta(double value) {
  mad_energy_struct_set_beta(fortran_ptr_, value);
}
double MadEnergyProxy::gamma() const {
  double value;
  mad_energy_struct_get_gamma(fortran_ptr_, &value);
  return value;
}
void MadEnergyProxy::set_gamma(double value) {
  mad_energy_struct_set_gamma(fortran_ptr_, value);
}
double MadEnergyProxy::kinetic() const {
  double value;
  mad_energy_struct_get_kinetic(fortran_ptr_, &value);
  return value;
}
void MadEnergyProxy::set_kinetic(double value) {
  mad_energy_struct_set_kinetic(fortran_ptr_, value);
}
double MadEnergyProxy::p0c() const {
  double value;
  mad_energy_struct_get_p0c(fortran_ptr_, &value);
  return value;
}
void MadEnergyProxy::set_p0c(double value) {
  mad_energy_struct_set_p0c(fortran_ptr_, value);
}
int MadEnergyProxy::particle() const {
  int value;
  mad_energy_struct_get_particle(fortran_ptr_, &value);
  return value;
}
void MadEnergyProxy::set_particle(int value) {
  mad_energy_struct_set_particle(fortran_ptr_, value);
}
FArray1D<double> MadMapProxy::k() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, mad_map_struct_get_k_info);
}
FArray2D<double> MadMapProxy::r() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, mad_map_struct_get_r_info);
}
FArray3D<double> MadMapProxy::t() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, mad_map_struct_get_t_info);
}
double AllEncompassingProxy::real_rp_0d() const {
  double value;
  all_encompassing_struct_get_real_rp_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingProxy::set_real_rp_0d(double value) {
  all_encompassing_struct_set_real_rp_0d(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingProxy::real_rp_1d() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_1d_info);
}
FArray2D<double> AllEncompassingProxy::real_rp_2d() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_2d_info);
}
FArray3D<double> AllEncompassingProxy::real_rp_3d() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_3d_info);
}
double* AllEncompassingProxy::real_rp_0d_ptr() const {
  double* ptr;
  all_encompassing_struct_get_real_rp_0d_ptr(fortran_ptr_, &ptr);
  return ptr;
}
void AllEncompassingProxy::set_real_rp_0d_ptr(double value) {
  all_encompassing_struct_set_real_rp_0d_ptr(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingProxy::real_rp_1d_ptr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_1d_ptr_info);
}
FArray2D<double> AllEncompassingProxy::real_rp_2d_ptr() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_2d_ptr_info);
}
FArray3D<double> AllEncompassingProxy::real_rp_3d_ptr() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_3d_ptr_info);
}
FArray1D<double> AllEncompassingProxy::real_rp_1d_alloc() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_1d_alloc_info);
}
FArray2D<double> AllEncompassingProxy::real_rp_2d_alloc() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_2d_alloc_info);
}
FArray3D<double> AllEncompassingProxy::real_rp_3d_alloc() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_rp_3d_alloc_info);
}
double AllEncompassingProxy::real_dp_0d() const {
  double value;
  all_encompassing_struct_get_real_dp_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingProxy::set_real_dp_0d(double value) {
  all_encompassing_struct_set_real_dp_0d(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingProxy::real_dp_1d() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_1d_info);
}
FArray2D<double> AllEncompassingProxy::real_dp_2d() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_2d_info);
}
FArray3D<double> AllEncompassingProxy::real_dp_3d() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_3d_info);
}
double* AllEncompassingProxy::real_dp_0d_ptr() const {
  double* ptr;
  all_encompassing_struct_get_real_dp_0d_ptr(fortran_ptr_, &ptr);
  return ptr;
}
void AllEncompassingProxy::set_real_dp_0d_ptr(double value) {
  all_encompassing_struct_set_real_dp_0d_ptr(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingProxy::real_dp_1d_ptr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_1d_ptr_info);
}
FArray2D<double> AllEncompassingProxy::real_dp_2d_ptr() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_2d_ptr_info);
}
FArray3D<double> AllEncompassingProxy::real_dp_3d_ptr() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_3d_ptr_info);
}
FArray1D<double> AllEncompassingProxy::real_dp_1d_alloc() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_1d_alloc_info);
}
FArray2D<double> AllEncompassingProxy::real_dp_2d_alloc() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_2d_alloc_info);
}
FArray3D<double> AllEncompassingProxy::real_dp_3d_alloc() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_, all_encompassing_struct_get_real_dp_3d_alloc_info);
}
std::complex<double> AllEncompassingProxy::complex_dp_0d() const {
  std::complex<double> c_value;
  all_encompassing_struct_get_complex_dp_0d(fortran_ptr_, &c_value);
  return c_value;
}
void AllEncompassingProxy::set_complex_dp_0d(std::complex<double> value) {
  all_encompassing_struct_set_complex_dp_0d(fortran_ptr_, value);
}
FArray1D<std::complex<double>> AllEncompassingProxy::complex_dp_1d() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_1d_info);
}
FArray2D<std::complex<double>> AllEncompassingProxy::complex_dp_2d() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_2d_info);
}
FArray3D<std::complex<double>> AllEncompassingProxy::complex_dp_3d() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_3d_info);
}
FArray1D<std::complex<double>> AllEncompassingProxy::complex_dp_1d_ptr() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_1d_ptr_info);
}
FArray2D<std::complex<double>> AllEncompassingProxy::complex_dp_2d_ptr() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_2d_ptr_info);
}
FArray3D<std::complex<double>> AllEncompassingProxy::complex_dp_3d_ptr() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_3d_ptr_info);
}
FArray1D<std::complex<double>> AllEncompassingProxy::complex_dp_1d_alloc()
    const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_1d_alloc_info);
}
FArray2D<std::complex<double>> AllEncompassingProxy::complex_dp_2d_alloc()
    const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_2d_alloc_info);
}
FArray3D<std::complex<double>> AllEncompassingProxy::complex_dp_3d_alloc()
    const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_, all_encompassing_struct_get_complex_dp_3d_alloc_info);
}
int AllEncompassingProxy::int_0d() const {
  int value;
  all_encompassing_struct_get_int_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingProxy::set_int_0d(int value) {
  all_encompassing_struct_set_int_0d(fortran_ptr_, value);
}
FArray1D<int> AllEncompassingProxy::int_1d() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_1d_info);
}
FArray2D<int> AllEncompassingProxy::int_2d() const {
  return ProxyHelpers::get_array_2d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_2d_info);
}
FArray3D<int> AllEncompassingProxy::int_3d() const {
  return ProxyHelpers::get_array_3d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_3d_info);
}
int* AllEncompassingProxy::int_0d_ptr() const {
  int* ptr;
  all_encompassing_struct_get_int_0d_ptr(fortran_ptr_, &ptr);
  return ptr;
}
void AllEncompassingProxy::set_int_0d_ptr(int value) {
  all_encompassing_struct_set_int_0d_ptr(fortran_ptr_, value);
}
FArray1D<int> AllEncompassingProxy::int_1d_ptr() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_1d_ptr_info);
}
FArray2D<int> AllEncompassingProxy::int_2d_ptr() const {
  return ProxyHelpers::get_array_2d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_2d_ptr_info);
}
FArray3D<int> AllEncompassingProxy::int_3d_ptr() const {
  return ProxyHelpers::get_array_3d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_3d_ptr_info);
}
FArray1D<int> AllEncompassingProxy::int_1d_alloc() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_1d_alloc_info);
}
FArray2D<int> AllEncompassingProxy::int_2d_alloc() const {
  return ProxyHelpers::get_array_2d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_2d_alloc_info);
}
FArray3D<int> AllEncompassingProxy::int_3d_alloc() const {
  return ProxyHelpers::get_array_3d<int>(
      fortran_ptr_, all_encompassing_struct_get_int_3d_alloc_info);
}
int64_t AllEncompassingProxy::int8_0d() const {
  int64_t value;
  all_encompassing_struct_get_int8_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingProxy::set_int8_0d(int64_t value) {
  all_encompassing_struct_set_int8_0d(fortran_ptr_, value);
}
int64_t* AllEncompassingProxy::int8_0d_ptr() const {
  int64_t* ptr;
  all_encompassing_struct_get_int8_0d_ptr(fortran_ptr_, &ptr);
  return ptr;
}
void AllEncompassingProxy::set_int8_0d_ptr(int64_t value) {
  all_encompassing_struct_set_int8_0d_ptr(fortran_ptr_, value);
}
bool AllEncompassingProxy::logical_0d() const {
  bool value;
  all_encompassing_struct_get_logical_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingProxy::set_logical_0d(bool value) {
  all_encompassing_struct_set_logical_0d(fortran_ptr_, value);
}
bool* AllEncompassingProxy::logical_0d_ptr() const {
  bool* ptr;
  all_encompassing_struct_get_logical_0d_ptr(fortran_ptr_, &ptr);
  return ptr;
}
void AllEncompassingProxy::set_logical_0d_ptr(bool value) {
  all_encompassing_struct_set_logical_0d_ptr(fortran_ptr_, value);
}
TestSubProxy AllEncompassingProxy::type_0d() const {
  void* ptr;
  all_encompassing_struct_get_type_0d(fortran_ptr_, &ptr);
  return TestSubProxy(ptr);
}
void AllEncompassingProxy::set_type_0d(const TestSubProxy& src) {
  all_encompassing_struct_set_type_0d(fortran_ptr_, src.get_fortran_ptr());
}
TestSubProxyArray1D AllEncompassingProxy::type_1d() const {
  return ProxyHelpers::get_type_array_1d<TestSubProxyArray1D>(
      fortran_ptr_, all_encompassing_struct_get_type_1d_info);
}
TestSubProxyArray2D AllEncompassingProxy::type_2d() const {
  return ProxyHelpers::get_type_array_2d<TestSubProxyArray2D>(
      fortran_ptr_, all_encompassing_struct_get_type_2d_info);
}
TestSubProxyArray3D AllEncompassingProxy::type_3d() const {
  return ProxyHelpers::get_type_array_3d<TestSubProxyArray3D>(
      fortran_ptr_, all_encompassing_struct_get_type_3d_info);
}
std::optional<TestSubProxy> AllEncompassingProxy::type_0d_ptr() const {
  void* ptr;
  all_encompassing_struct_get_type_0d_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TestSubProxy(ptr);
}
void AllEncompassingProxy::set_type_0d_ptr(const TestSubProxy& src) {
  all_encompassing_struct_set_type_0d_ptr(fortran_ptr_, src.get_fortran_ptr());
}
TestSubProxyArray1D AllEncompassingProxy::type_1d_ptr() const {
  return ProxyHelpers::get_type_array_1d<TestSubProxyArray1D>(
      fortran_ptr_, all_encompassing_struct_get_type_1d_ptr_info);
}
TestSubProxyArray2D AllEncompassingProxy::type_2d_ptr() const {
  return ProxyHelpers::get_type_array_2d<TestSubProxyArray2D>(
      fortran_ptr_, all_encompassing_struct_get_type_2d_ptr_info);
}
TestSubProxyArray3D AllEncompassingProxy::type_3d_ptr() const {
  return ProxyHelpers::get_type_array_3d<TestSubProxyArray3D>(
      fortran_ptr_, all_encompassing_struct_get_type_3d_ptr_info);
}
TestSubProxyArray1D AllEncompassingProxy::type_1d_alloc() const {
  return ProxyHelpers::get_type_array_1d<TestSubProxyArray1D>(
      fortran_ptr_, all_encompassing_struct_get_type_1d_alloc_info);
}
TestSubProxyArray2D AllEncompassingProxy::type_2d_alloc() const {
  return ProxyHelpers::get_type_array_2d<TestSubProxyArray2D>(
      fortran_ptr_, all_encompassing_struct_get_type_2d_alloc_info);
}
TestSubProxyArray3D AllEncompassingProxy::type_3d_alloc() const {
  return ProxyHelpers::get_type_array_3d<TestSubProxyArray3D>(
      fortran_ptr_, all_encompassing_struct_get_type_3d_alloc_info);
}
TestSubSubProxy TestSubProxy::sr() const {
  void* ptr;
  test_sub_struct_get_sr(fortran_ptr_, &ptr);
  return TestSubSubProxy(ptr);
}
void TestSubProxy::set_sr(const TestSubSubProxy& src) {
  test_sub_struct_set_sr(fortran_ptr_, src.get_fortran_ptr());
}
int64_t TestSubSubProxy::aaa() const {
  int64_t value;
  test_sub_sub_struct_get_aaa(fortran_ptr_, &value);
  return value;
}
void TestSubSubProxy::set_aaa(int64_t value) {
  test_sub_sub_struct_set_aaa(fortran_ptr_, value);
}
int TestSubSubProxy::bbb() const {
  int value;
  test_sub_sub_struct_get_bbb(fortran_ptr_, &value);
  return value;
}
void TestSubSubProxy::set_bbb(int value) {
  test_sub_sub_struct_set_bbb(fortran_ptr_, value);
}
std::string TestSubSubProxy::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_, test_sub_sub_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void TestSubSubProxy::set_file(const std::string& value) {
  test_sub_sub_struct_set_file(
      fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TestSubSubProxy::t_ref() const {
  double value;
  test_sub_sub_struct_get_t_ref(fortran_ptr_, &value);
  return value;
}
void TestSubSubProxy::set_t_ref(double value) {
  test_sub_sub_struct_set_t_ref(fortran_ptr_, value);
}
double TestSubSubProxy::freq_spread() const {
  double value;
  test_sub_sub_struct_get_freq_spread(fortran_ptr_, &value);
  return value;
}
void TestSubSubProxy::set_freq_spread(double value) {
  test_sub_sub_struct_set_freq_spread(fortran_ptr_, value);
}