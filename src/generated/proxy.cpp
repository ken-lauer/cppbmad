#include "bmad/generated/proxy.hpp"

using namespace Bmad;
double SplineStruct::x0() const {
  double value;
  spline_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SplineStruct::set_x0(double value) { spline_struct_set_x0(fortran_ptr_, value); }
double SplineStruct::y0() const {
  double value;
  spline_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SplineStruct::set_y0(double value) { spline_struct_set_y0(fortran_ptr_, value); }
double SplineStruct::x1() const {
  double value;
  spline_struct_get_x1(fortran_ptr_, &value);
  return value;
}
void SplineStruct::set_x1(double value) { spline_struct_set_x1(fortran_ptr_, value); }
FArray1D<double> SplineStruct::coef() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, spline_struct_get_coef_info);
}
void SplineStruct::set_coef(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  spline_struct_set_coef(fortran_ptr_, v.data(), shape);
}
double SpinPolarStruct::polarization() const {
  double value;
  spin_polar_struct_get_polarization(fortran_ptr_, &value);
  return value;
}
void SpinPolarStruct::set_polarization(double value) {
  spin_polar_struct_set_polarization(fortran_ptr_, value);
}
double SpinPolarStruct::theta() const {
  double value;
  spin_polar_struct_get_theta(fortran_ptr_, &value);
  return value;
}
void SpinPolarStruct::set_theta(double value) { spin_polar_struct_set_theta(fortran_ptr_, value); }
double SpinPolarStruct::phi() const {
  double value;
  spin_polar_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void SpinPolarStruct::set_phi(double value) { spin_polar_struct_set_phi(fortran_ptr_, value); }
double SpinPolarStruct::xi() const {
  double value;
  spin_polar_struct_get_xi(fortran_ptr_, &value);
  return value;
}
void SpinPolarStruct::set_xi(double value) { spin_polar_struct_set_xi(fortran_ptr_, value); }
double AcKickerTimeStruct::amp() const {
  double value;
  ac_kicker_time_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void AcKickerTimeStruct::set_amp(double value) {
  ac_kicker_time_struct_set_amp(fortran_ptr_, value);
}
double AcKickerTimeStruct::time() const {
  double value;
  ac_kicker_time_struct_get_time(fortran_ptr_, &value);
  return value;
}
void AcKickerTimeStruct::set_time(double value) {
  ac_kicker_time_struct_set_time(fortran_ptr_, value);
}
SplineStruct AcKickerTimeStruct::spline() const {
  void *ptr;
  ac_kicker_time_struct_get_spline(fortran_ptr_, &ptr);
  return SplineStruct(ptr);
}
void AcKickerTimeStruct::set_spline(const SplineStruct &src) {
  ac_kicker_time_struct_set_spline(fortran_ptr_, src.get_fortran_ptr());
}
double AcKickerFreqStruct::f() const {
  double value;
  ac_kicker_freq_struct_get_f(fortran_ptr_, &value);
  return value;
}
void AcKickerFreqStruct::set_f(double value) { ac_kicker_freq_struct_set_f(fortran_ptr_, value); }
double AcKickerFreqStruct::amp() const {
  double value;
  ac_kicker_freq_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void AcKickerFreqStruct::set_amp(double value) {
  ac_kicker_freq_struct_set_amp(fortran_ptr_, value);
}
double AcKickerFreqStruct::phi() const {
  double value;
  ac_kicker_freq_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void AcKickerFreqStruct::set_phi(double value) {
  ac_kicker_freq_struct_set_phi(fortran_ptr_, value);
}
AcKickerTimeStructAlloc1D AcKickerStruct::amp_vs_time() const {
  return AcKickerTimeStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      ac_kicker_struct_reallocate_amp_vs_time,
      ac_kicker_struct_get_amp_vs_time_info
  );
}
AcKickerFreqStructAlloc1D AcKickerStruct::frequency() const {
  return AcKickerFreqStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      ac_kicker_struct_reallocate_frequency,
      ac_kicker_struct_get_frequency_info
  );
}
double Interval1CoefStruct::c0() const {
  double value;
  interval1_coef_struct_get_c0(fortran_ptr_, &value);
  return value;
}
void Interval1CoefStruct::set_c0(double value) {
  interval1_coef_struct_set_c0(fortran_ptr_, value);
}
double Interval1CoefStruct::c1() const {
  double value;
  interval1_coef_struct_get_c1(fortran_ptr_, &value);
  return value;
}
void Interval1CoefStruct::set_c1(double value) {
  interval1_coef_struct_set_c1(fortran_ptr_, value);
}
double Interval1CoefStruct::n_exp() const {
  double value;
  interval1_coef_struct_get_n_exp(fortran_ptr_, &value);
  return value;
}
void Interval1CoefStruct::set_n_exp(double value) {
  interval1_coef_struct_set_n_exp(fortran_ptr_, value);
}
RealAlloc1D PhotonReflectTableStruct::angle() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_reflect_table_struct_reallocate_angle,
      photon_reflect_table_struct_get_angle_info
  );
}
void PhotonReflectTableStruct::set_angle(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_reflect_table_struct_set_angle(fortran_ptr_, v.data(), shape);
}
RealAlloc1D PhotonReflectTableStruct::energy() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_reflect_table_struct_reallocate_energy,
      photon_reflect_table_struct_get_energy_info
  );
}
void PhotonReflectTableStruct::set_energy(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_reflect_table_struct_set_energy(fortran_ptr_, v.data(), shape);
}
Interval1CoefStructAlloc1D PhotonReflectTableStruct::int1() const {
  return Interval1CoefStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_reflect_table_struct_reallocate_int1,
      photon_reflect_table_struct_get_int1_info
  );
}
FArray2D<double> PhotonReflectTableStruct::p_reflect() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      photon_reflect_table_struct_get_p_reflect_info
  );
}
void PhotonReflectTableStruct::set_p_reflect(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, photon_reflect_table_struct_set_p_reflect, v);
}
double PhotonReflectTableStruct::max_energy() const {
  double value;
  photon_reflect_table_struct_get_max_energy(fortran_ptr_, &value);
  return value;
}
void PhotonReflectTableStruct::set_max_energy(double value) {
  photon_reflect_table_struct_set_max_energy(fortran_ptr_, value);
}
RealAlloc1D PhotonReflectTableStruct::p_reflect_scratch() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_reflect_table_struct_reallocate_p_reflect_scratch,
      photon_reflect_table_struct_get_p_reflect_scratch_info
  );
}
void PhotonReflectTableStruct::set_p_reflect_scratch(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_reflect_table_struct_set_p_reflect_scratch(fortran_ptr_, v.data(), shape);
}
RealAlloc1D PhotonReflectTableStruct::bragg_angle() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_reflect_table_struct_reallocate_bragg_angle,
      photon_reflect_table_struct_get_bragg_angle_info
  );
}
void PhotonReflectTableStruct::set_bragg_angle(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_reflect_table_struct_set_bragg_angle(fortran_ptr_, v.data(), shape);
}
std::string PhotonReflectSurfaceStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, photon_reflect_surface_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void PhotonReflectSurfaceStruct::set_name(const std::string &value) {
  photon_reflect_surface_struct_set_name(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string PhotonReflectSurfaceStruct::description() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      photon_reflect_surface_struct_get_description_info
  );
  return std::string(arr.data(), arr.size());
}
void PhotonReflectSurfaceStruct::set_description(const std::string &value) {
  photon_reflect_surface_struct_set_description(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string PhotonReflectSurfaceStruct::reflectivity_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      photon_reflect_surface_struct_get_reflectivity_file_info
  );
  return std::string(arr.data(), arr.size());
}
void PhotonReflectSurfaceStruct::set_reflectivity_file(const std::string &value) {
  photon_reflect_surface_struct_set_reflectivity_file(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
PhotonReflectTableStructAlloc1D PhotonReflectSurfaceStruct::table() const {
  return PhotonReflectTableStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_reflect_surface_struct_reallocate_table,
      photon_reflect_surface_struct_get_table_info
  );
}
double PhotonReflectSurfaceStruct::surface_roughness_rms() const {
  double value;
  photon_reflect_surface_struct_get_surface_roughness_rms(fortran_ptr_, &value);
  return value;
}
void PhotonReflectSurfaceStruct::set_surface_roughness_rms(double value) {
  photon_reflect_surface_struct_set_surface_roughness_rms(fortran_ptr_, value);
}
double PhotonReflectSurfaceStruct::roughness_correlation_len() const {
  double value;
  photon_reflect_surface_struct_get_roughness_correlation_len(fortran_ptr_, &value);
  return value;
}
void PhotonReflectSurfaceStruct::set_roughness_correlation_len(double value) {
  photon_reflect_surface_struct_set_roughness_correlation_len(fortran_ptr_, value);
}
int PhotonReflectSurfaceStruct::ix_surface() const {
  int value;
  photon_reflect_surface_struct_get_ix_surface(fortran_ptr_, &value);
  return value;
}
void PhotonReflectSurfaceStruct::set_ix_surface(int value) {
  photon_reflect_surface_struct_set_ix_surface(fortran_ptr_, value);
}
FArray1D<double> CoordStruct::vec() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, coord_struct_get_vec_info);
}
void CoordStruct::set_vec(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  coord_struct_set_vec(fortran_ptr_, v.data(), shape);
}
double CoordStruct::s() const {
  double value;
  coord_struct_get_s(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_s(double value) { coord_struct_set_s(fortran_ptr_, value); }
long double CoordStruct::t() const {
  long double value;
  coord_struct_get_t(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_t(long double value) { coord_struct_set_t(fortran_ptr_, value); }
FArray1D<double> CoordStruct::spin() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, coord_struct_get_spin_info);
}
void CoordStruct::set_spin(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  coord_struct_set_spin(fortran_ptr_, v.data(), shape);
}
FArray1D<double> CoordStruct::field() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, coord_struct_get_field_info);
}
void CoordStruct::set_field(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  coord_struct_set_field(fortran_ptr_, v.data(), shape);
}
FArray1D<double> CoordStruct::phase() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, coord_struct_get_phase_info);
}
void CoordStruct::set_phase(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  coord_struct_set_phase(fortran_ptr_, v.data(), shape);
}
double CoordStruct::charge() const {
  double value;
  coord_struct_get_charge(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_charge(double value) { coord_struct_set_charge(fortran_ptr_, value); }
double CoordStruct::dt_ref() const {
  double value;
  coord_struct_get_dt_ref(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_dt_ref(double value) { coord_struct_set_dt_ref(fortran_ptr_, value); }
double CoordStruct::r() const {
  double value;
  coord_struct_get_r(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_r(double value) { coord_struct_set_r(fortran_ptr_, value); }
double CoordStruct::p0c() const {
  double value;
  coord_struct_get_p0c(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_p0c(double value) { coord_struct_set_p0c(fortran_ptr_, value); }
double CoordStruct::E_potential() const {
  double value;
  coord_struct_get_E_potential(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_E_potential(double value) {
  coord_struct_set_E_potential(fortran_ptr_, value);
}
double CoordStruct::beta() const {
  double value;
  coord_struct_get_beta(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_beta(double value) { coord_struct_set_beta(fortran_ptr_, value); }
int CoordStruct::ix_ele() const {
  int value;
  coord_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_ix_ele(int value) { coord_struct_set_ix_ele(fortran_ptr_, value); }
int CoordStruct::ix_branch() const {
  int value;
  coord_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_ix_branch(int value) { coord_struct_set_ix_branch(fortran_ptr_, value); }
int CoordStruct::ix_turn() const {
  int value;
  coord_struct_get_ix_turn(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_ix_turn(int value) { coord_struct_set_ix_turn(fortran_ptr_, value); }
int CoordStruct::ix_user() const {
  int value;
  coord_struct_get_ix_user(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_ix_user(int value) { coord_struct_set_ix_user(fortran_ptr_, value); }
int CoordStruct::state() const {
  int value;
  coord_struct_get_state(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_state(int value) { coord_struct_set_state(fortran_ptr_, value); }
int CoordStruct::direction() const {
  int value;
  coord_struct_get_direction(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_direction(int value) { coord_struct_set_direction(fortran_ptr_, value); }
int CoordStruct::time_dir() const {
  int value;
  coord_struct_get_time_dir(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_time_dir(int value) { coord_struct_set_time_dir(fortran_ptr_, value); }
int CoordStruct::species() const {
  int value;
  coord_struct_get_species(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_species(int value) { coord_struct_set_species(fortran_ptr_, value); }
int CoordStruct::location() const {
  int value;
  coord_struct_get_location(fortran_ptr_, &value);
  return value;
}
void CoordStruct::set_location(int value) { coord_struct_set_location(fortran_ptr_, value); }
CoordStructAlloc1D CoordArrayStruct::orbit() const {
  return CoordStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      coord_array_struct_reallocate_orbit,
      coord_array_struct_get_orbit_info
  );
}
double BpmPhaseCouplingStruct::K_22a() const {
  double value;
  bpm_phase_coupling_struct_get_K_22a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_K_22a(double value) {
  bpm_phase_coupling_struct_set_K_22a(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::K_12a() const {
  double value;
  bpm_phase_coupling_struct_get_K_12a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_K_12a(double value) {
  bpm_phase_coupling_struct_set_K_12a(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::K_11b() const {
  double value;
  bpm_phase_coupling_struct_get_K_11b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_K_11b(double value) {
  bpm_phase_coupling_struct_set_K_11b(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::K_12b() const {
  double value;
  bpm_phase_coupling_struct_get_K_12b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_K_12b(double value) {
  bpm_phase_coupling_struct_set_K_12b(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::Cbar22_a() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar22_a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_Cbar22_a(double value) {
  bpm_phase_coupling_struct_set_Cbar22_a(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::Cbar12_a() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar12_a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_Cbar12_a(double value) {
  bpm_phase_coupling_struct_set_Cbar12_a(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::Cbar11_b() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar11_b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_Cbar11_b(double value) {
  bpm_phase_coupling_struct_set_Cbar11_b(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::Cbar12_b() const {
  double value;
  bpm_phase_coupling_struct_get_Cbar12_b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_Cbar12_b(double value) {
  bpm_phase_coupling_struct_set_Cbar12_b(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::phi_a() const {
  double value;
  bpm_phase_coupling_struct_get_phi_a(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_phi_a(double value) {
  bpm_phase_coupling_struct_set_phi_a(fortran_ptr_, value);
}
double BpmPhaseCouplingStruct::phi_b() const {
  double value;
  bpm_phase_coupling_struct_get_phi_b(fortran_ptr_, &value);
  return value;
}
void BpmPhaseCouplingStruct::set_phi_b(double value) {
  bpm_phase_coupling_struct_set_phi_b(fortran_ptr_, value);
}
std::string ExpressionAtomStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, expression_atom_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void ExpressionAtomStruct::set_name(const std::string &value) {
  expression_atom_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int ExpressionAtomStruct::type() const {
  int value;
  expression_atom_struct_get_type(fortran_ptr_, &value);
  return value;
}
void ExpressionAtomStruct::set_type(int value) {
  expression_atom_struct_set_type(fortran_ptr_, value);
}
double ExpressionAtomStruct::value() const {
  double value;
  expression_atom_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ExpressionAtomStruct::set_value(double value) {
  expression_atom_struct_set_value(fortran_ptr_, value);
}
RealAlloc1D WakeSrZLongStruct::w() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_sr_z_long_struct_reallocate_w,
      wake_sr_z_long_struct_get_w_info
  );
}
void WakeSrZLongStruct::set_w(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wake_sr_z_long_struct_set_w(fortran_ptr_, v.data(), shape);
}
ComplexAlloc1D WakeSrZLongStruct::fw() const {
  return ComplexAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_sr_z_long_struct_reallocate_fw,
      wake_sr_z_long_struct_get_fw_info
  );
}
void WakeSrZLongStruct::set_fw(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wake_sr_z_long_struct_set_fw(fortran_ptr_, v.data(), shape);
}
ComplexAlloc1D WakeSrZLongStruct::fbunch() const {
  return ComplexAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_sr_z_long_struct_reallocate_fbunch,
      wake_sr_z_long_struct_get_fbunch_info
  );
}
void WakeSrZLongStruct::set_fbunch(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wake_sr_z_long_struct_set_fbunch(fortran_ptr_, v.data(), shape);
}
ComplexAlloc1D WakeSrZLongStruct::w_out() const {
  return ComplexAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_sr_z_long_struct_reallocate_w_out,
      wake_sr_z_long_struct_get_w_out_info
  );
}
void WakeSrZLongStruct::set_w_out(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wake_sr_z_long_struct_set_w_out(fortran_ptr_, v.data(), shape);
}
double WakeSrZLongStruct::dz() const {
  double value;
  wake_sr_z_long_struct_get_dz(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongStruct::set_dz(double value) { wake_sr_z_long_struct_set_dz(fortran_ptr_, value); }
double WakeSrZLongStruct::z0() const {
  double value;
  wake_sr_z_long_struct_get_z0(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongStruct::set_z0(double value) { wake_sr_z_long_struct_set_z0(fortran_ptr_, value); }
double WakeSrZLongStruct::smoothing_sigma() const {
  double value;
  wake_sr_z_long_struct_get_smoothing_sigma(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongStruct::set_smoothing_sigma(double value) {
  wake_sr_z_long_struct_set_smoothing_sigma(fortran_ptr_, value);
}
int WakeSrZLongStruct::position_dependence() const {
  int value;
  wake_sr_z_long_struct_get_position_dependence(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongStruct::set_position_dependence(int value) {
  wake_sr_z_long_struct_set_position_dependence(fortran_ptr_, value);
}
bool WakeSrZLongStruct::time_based() const {
  bool value;
  wake_sr_z_long_struct_get_time_based(fortran_ptr_, &value);
  return value;
}
void WakeSrZLongStruct::set_time_based(bool value) {
  wake_sr_z_long_struct_set_time_based(fortran_ptr_, value);
}
double WakeSrModeStruct::amp() const {
  double value;
  wake_sr_mode_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_amp(double value) { wake_sr_mode_struct_set_amp(fortran_ptr_, value); }
double WakeSrModeStruct::damp() const {
  double value;
  wake_sr_mode_struct_get_damp(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_damp(double value) { wake_sr_mode_struct_set_damp(fortran_ptr_, value); }
double WakeSrModeStruct::k() const {
  double value;
  wake_sr_mode_struct_get_k(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_k(double value) { wake_sr_mode_struct_set_k(fortran_ptr_, value); }
double WakeSrModeStruct::phi() const {
  double value;
  wake_sr_mode_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_phi(double value) { wake_sr_mode_struct_set_phi(fortran_ptr_, value); }
double WakeSrModeStruct::b_sin() const {
  double value;
  wake_sr_mode_struct_get_b_sin(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_b_sin(double value) {
  wake_sr_mode_struct_set_b_sin(fortran_ptr_, value);
}
double WakeSrModeStruct::b_cos() const {
  double value;
  wake_sr_mode_struct_get_b_cos(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_b_cos(double value) {
  wake_sr_mode_struct_set_b_cos(fortran_ptr_, value);
}
double WakeSrModeStruct::a_sin() const {
  double value;
  wake_sr_mode_struct_get_a_sin(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_a_sin(double value) {
  wake_sr_mode_struct_set_a_sin(fortran_ptr_, value);
}
double WakeSrModeStruct::a_cos() const {
  double value;
  wake_sr_mode_struct_get_a_cos(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_a_cos(double value) {
  wake_sr_mode_struct_set_a_cos(fortran_ptr_, value);
}
int WakeSrModeStruct::polarization() const {
  int value;
  wake_sr_mode_struct_get_polarization(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_polarization(int value) {
  wake_sr_mode_struct_set_polarization(fortran_ptr_, value);
}
int WakeSrModeStruct::position_dependence() const {
  int value;
  wake_sr_mode_struct_get_position_dependence(fortran_ptr_, &value);
  return value;
}
void WakeSrModeStruct::set_position_dependence(int value) {
  wake_sr_mode_struct_set_position_dependence(fortran_ptr_, value);
}
std::string WakeSrStruct::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, wake_sr_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void WakeSrStruct::set_file(const std::string &value) {
  wake_sr_struct_set_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
WakeSrZLongStruct WakeSrStruct::z_long() const {
  void *ptr;
  wake_sr_struct_get_z_long(fortran_ptr_, &ptr);
  return WakeSrZLongStruct(ptr);
}
void WakeSrStruct::set_z_long(const WakeSrZLongStruct &src) {
  wake_sr_struct_set_z_long(fortran_ptr_, src.get_fortran_ptr());
}
WakeSrModeStructAlloc1D WakeSrStruct::long_wake() const {
  return WakeSrModeStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_sr_struct_reallocate_long,
      wake_sr_struct_get_long_info
  );
}
WakeSrModeStructAlloc1D WakeSrStruct::trans_wake() const {
  return WakeSrModeStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_sr_struct_reallocate_trans,
      wake_sr_struct_get_trans_info
  );
}
double WakeSrStruct::z_ref_long() const {
  double value;
  wake_sr_struct_get_z_ref_long(fortran_ptr_, &value);
  return value;
}
void WakeSrStruct::set_z_ref_long(double value) {
  wake_sr_struct_set_z_ref_long(fortran_ptr_, value);
}
double WakeSrStruct::z_ref_trans() const {
  double value;
  wake_sr_struct_get_z_ref_trans(fortran_ptr_, &value);
  return value;
}
void WakeSrStruct::set_z_ref_trans(double value) {
  wake_sr_struct_set_z_ref_trans(fortran_ptr_, value);
}
double WakeSrStruct::z_max() const {
  double value;
  wake_sr_struct_get_z_max(fortran_ptr_, &value);
  return value;
}
void WakeSrStruct::set_z_max(double value) { wake_sr_struct_set_z_max(fortran_ptr_, value); }
double WakeSrStruct::amp_scale() const {
  double value;
  wake_sr_struct_get_amp_scale(fortran_ptr_, &value);
  return value;
}
void WakeSrStruct::set_amp_scale(double value) {
  wake_sr_struct_set_amp_scale(fortran_ptr_, value);
}
double WakeSrStruct::z_scale() const {
  double value;
  wake_sr_struct_get_z_scale(fortran_ptr_, &value);
  return value;
}
void WakeSrStruct::set_z_scale(double value) { wake_sr_struct_set_z_scale(fortran_ptr_, value); }
bool WakeSrStruct::scale_with_length() const {
  bool value;
  wake_sr_struct_get_scale_with_length(fortran_ptr_, &value);
  return value;
}
void WakeSrStruct::set_scale_with_length(bool value) {
  wake_sr_struct_set_scale_with_length(fortran_ptr_, value);
}
double WakeLrModeStruct::freq() const {
  double value;
  wake_lr_mode_struct_get_freq(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_freq(double value) { wake_lr_mode_struct_set_freq(fortran_ptr_, value); }
double WakeLrModeStruct::freq_in() const {
  double value;
  wake_lr_mode_struct_get_freq_in(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_freq_in(double value) {
  wake_lr_mode_struct_set_freq_in(fortran_ptr_, value);
}
double WakeLrModeStruct::R_over_Q() const {
  double value;
  wake_lr_mode_struct_get_R_over_Q(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_R_over_Q(double value) {
  wake_lr_mode_struct_set_R_over_Q(fortran_ptr_, value);
}
double WakeLrModeStruct::Q() const {
  double value;
  wake_lr_mode_struct_get_Q(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_Q(double value) { wake_lr_mode_struct_set_Q(fortran_ptr_, value); }
double WakeLrModeStruct::damp() const {
  double value;
  wake_lr_mode_struct_get_damp(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_damp(double value) { wake_lr_mode_struct_set_damp(fortran_ptr_, value); }
double WakeLrModeStruct::phi() const {
  double value;
  wake_lr_mode_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_phi(double value) { wake_lr_mode_struct_set_phi(fortran_ptr_, value); }
double WakeLrModeStruct::angle() const {
  double value;
  wake_lr_mode_struct_get_angle(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_angle(double value) {
  wake_lr_mode_struct_set_angle(fortran_ptr_, value);
}
double WakeLrModeStruct::b_sin() const {
  double value;
  wake_lr_mode_struct_get_b_sin(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_b_sin(double value) {
  wake_lr_mode_struct_set_b_sin(fortran_ptr_, value);
}
double WakeLrModeStruct::b_cos() const {
  double value;
  wake_lr_mode_struct_get_b_cos(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_b_cos(double value) {
  wake_lr_mode_struct_set_b_cos(fortran_ptr_, value);
}
double WakeLrModeStruct::a_sin() const {
  double value;
  wake_lr_mode_struct_get_a_sin(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_a_sin(double value) {
  wake_lr_mode_struct_set_a_sin(fortran_ptr_, value);
}
double WakeLrModeStruct::a_cos() const {
  double value;
  wake_lr_mode_struct_get_a_cos(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_a_cos(double value) {
  wake_lr_mode_struct_set_a_cos(fortran_ptr_, value);
}
int WakeLrModeStruct::m() const {
  int value;
  wake_lr_mode_struct_get_m(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_m(int value) { wake_lr_mode_struct_set_m(fortran_ptr_, value); }
bool WakeLrModeStruct::polarized() const {
  bool value;
  wake_lr_mode_struct_get_polarized(fortran_ptr_, &value);
  return value;
}
void WakeLrModeStruct::set_polarized(bool value) {
  wake_lr_mode_struct_set_polarized(fortran_ptr_, value);
}
std::string WakeLrStruct::file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, wake_lr_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void WakeLrStruct::set_file(const std::string &value) {
  wake_lr_struct_set_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
WakeLrModeStructAlloc1D WakeLrStruct::mode() const {
  return WakeLrModeStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wake_lr_struct_reallocate_mode,
      wake_lr_struct_get_mode_info
  );
}
double WakeLrStruct::t_ref() const {
  double value;
  wake_lr_struct_get_t_ref(fortran_ptr_, &value);
  return value;
}
void WakeLrStruct::set_t_ref(double value) { wake_lr_struct_set_t_ref(fortran_ptr_, value); }
double WakeLrStruct::freq_spread() const {
  double value;
  wake_lr_struct_get_freq_spread(fortran_ptr_, &value);
  return value;
}
void WakeLrStruct::set_freq_spread(double value) {
  wake_lr_struct_set_freq_spread(fortran_ptr_, value);
}
double WakeLrStruct::amp_scale() const {
  double value;
  wake_lr_struct_get_amp_scale(fortran_ptr_, &value);
  return value;
}
void WakeLrStruct::set_amp_scale(double value) {
  wake_lr_struct_set_amp_scale(fortran_ptr_, value);
}
double WakeLrStruct::time_scale() const {
  double value;
  wake_lr_struct_get_time_scale(fortran_ptr_, &value);
  return value;
}
void WakeLrStruct::set_time_scale(double value) {
  wake_lr_struct_set_time_scale(fortran_ptr_, value);
}
bool WakeLrStruct::self_wake_on() const {
  bool value;
  wake_lr_struct_get_self_wake_on(fortran_ptr_, &value);
  return value;
}
void WakeLrStruct::set_self_wake_on(bool value) {
  wake_lr_struct_set_self_wake_on(fortran_ptr_, value);
}
int LatEleLocStruct::ix_ele() const {
  int value;
  lat_ele_loc_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void LatEleLocStruct::set_ix_ele(int value) { lat_ele_loc_struct_set_ix_ele(fortran_ptr_, value); }
int LatEleLocStruct::ix_branch() const {
  int value;
  lat_ele_loc_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void LatEleLocStruct::set_ix_branch(int value) {
  lat_ele_loc_struct_set_ix_branch(fortran_ptr_, value);
}
WakeSrStruct WakeStruct::sr() const {
  void *ptr;
  wake_struct_get_sr(fortran_ptr_, &ptr);
  return WakeSrStruct(ptr);
}
void WakeStruct::set_sr(const WakeSrStruct &src) {
  wake_struct_set_sr(fortran_ptr_, src.get_fortran_ptr());
}
WakeLrStruct WakeStruct::lr() const {
  void *ptr;
  wake_struct_get_lr(fortran_ptr_, &ptr);
  return WakeLrStruct(ptr);
}
void WakeStruct::set_lr(const WakeLrStruct &src) {
  wake_struct_set_lr(fortran_ptr_, src.get_fortran_ptr());
}
double TaylorTermStruct::coef() const {
  double value;
  taylor_term_struct_get_coef(fortran_ptr_, &value);
  return value;
}
void TaylorTermStruct::set_coef(double value) { taylor_term_struct_set_coef(fortran_ptr_, value); }
FArray1D<int> TaylorTermStruct::expn() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, taylor_term_struct_get_expn_info);
}
void TaylorTermStruct::set_expn(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  taylor_term_struct_set_expn(fortran_ptr_, v.data(), shape);
}
double TaylorStruct::ref() const {
  double value;
  taylor_struct_get_ref(fortran_ptr_, &value);
  return value;
}
void TaylorStruct::set_ref(double value) { taylor_struct_set_ref(fortran_ptr_, value); }
TaylorTermStructArray1D TaylorStruct::term() const {
  return ProxyHelpers::get_type_array_1d<TaylorTermStructArray1D>(
      fortran_ptr_,
      taylor_struct_get_term_info
  );
}
double EmTaylorTermStruct::coef() const {
  double value;
  em_taylor_term_struct_get_coef(fortran_ptr_, &value);
  return value;
}
void EmTaylorTermStruct::set_coef(double value) {
  em_taylor_term_struct_set_coef(fortran_ptr_, value);
}
FArray1D<int> EmTaylorTermStruct::expn() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, em_taylor_term_struct_get_expn_info);
}
void EmTaylorTermStruct::set_expn(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  em_taylor_term_struct_set_expn(fortran_ptr_, v.data(), shape);
}
double EmTaylorStruct::ref() const {
  double value;
  em_taylor_struct_get_ref(fortran_ptr_, &value);
  return value;
}
void EmTaylorStruct::set_ref(double value) { em_taylor_struct_set_ref(fortran_ptr_, value); }
EmTaylorTermStructAlloc1D EmTaylorStruct::term() const {
  return EmTaylorTermStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      em_taylor_struct_reallocate_term,
      em_taylor_struct_get_term_info
  );
}
double CartesianMapTerm1Struct::coef() const {
  double value;
  cartesian_map_term1_struct_get_coef(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_coef(double value) {
  cartesian_map_term1_struct_set_coef(fortran_ptr_, value);
}
double CartesianMapTerm1Struct::kx() const {
  double value;
  cartesian_map_term1_struct_get_kx(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_kx(double value) {
  cartesian_map_term1_struct_set_kx(fortran_ptr_, value);
}
double CartesianMapTerm1Struct::ky() const {
  double value;
  cartesian_map_term1_struct_get_ky(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_ky(double value) {
  cartesian_map_term1_struct_set_ky(fortran_ptr_, value);
}
double CartesianMapTerm1Struct::kz() const {
  double value;
  cartesian_map_term1_struct_get_kz(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_kz(double value) {
  cartesian_map_term1_struct_set_kz(fortran_ptr_, value);
}
double CartesianMapTerm1Struct::x0() const {
  double value;
  cartesian_map_term1_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_x0(double value) {
  cartesian_map_term1_struct_set_x0(fortran_ptr_, value);
}
double CartesianMapTerm1Struct::y0() const {
  double value;
  cartesian_map_term1_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_y0(double value) {
  cartesian_map_term1_struct_set_y0(fortran_ptr_, value);
}
double CartesianMapTerm1Struct::phi_z() const {
  double value;
  cartesian_map_term1_struct_get_phi_z(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_phi_z(double value) {
  cartesian_map_term1_struct_set_phi_z(fortran_ptr_, value);
}
int CartesianMapTerm1Struct::family() const {
  int value;
  cartesian_map_term1_struct_get_family(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_family(int value) {
  cartesian_map_term1_struct_set_family(fortran_ptr_, value);
}
int CartesianMapTerm1Struct::form() const {
  int value;
  cartesian_map_term1_struct_get_form(fortran_ptr_, &value);
  return value;
}
void CartesianMapTerm1Struct::set_form(int value) {
  cartesian_map_term1_struct_set_form(fortran_ptr_, value);
}
std::string CartesianMapTermStruct::file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, cartesian_map_term_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void CartesianMapTermStruct::set_file(const std::string &value) {
  cartesian_map_term_struct_set_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int CartesianMapTermStruct::n_link() const {
  int value;
  cartesian_map_term_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void CartesianMapTermStruct::set_n_link(int value) {
  cartesian_map_term_struct_set_n_link(fortran_ptr_, value);
}
CartesianMapTerm1StructAlloc1D CartesianMapTermStruct::term() const {
  return CartesianMapTerm1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      cartesian_map_term_struct_reallocate_term,
      cartesian_map_term_struct_get_term_info
  );
}
double CartesianMapStruct::field_scale() const {
  double value;
  cartesian_map_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void CartesianMapStruct::set_field_scale(double value) {
  cartesian_map_struct_set_field_scale(fortran_ptr_, value);
}
FArray1D<double> CartesianMapStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, cartesian_map_struct_get_r0_info);
}
void CartesianMapStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  cartesian_map_struct_set_r0(fortran_ptr_, v.data(), shape);
}
int CartesianMapStruct::master_parameter() const {
  int value;
  cartesian_map_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void CartesianMapStruct::set_master_parameter(int value) {
  cartesian_map_struct_set_master_parameter(fortran_ptr_, value);
}
int CartesianMapStruct::ele_anchor_pt() const {
  int value;
  cartesian_map_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void CartesianMapStruct::set_ele_anchor_pt(int value) {
  cartesian_map_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
int CartesianMapStruct::field_type() const {
  int value;
  cartesian_map_struct_get_field_type(fortran_ptr_, &value);
  return value;
}
void CartesianMapStruct::set_field_type(int value) {
  cartesian_map_struct_set_field_type(fortran_ptr_, value);
}
std::optional<CartesianMapTermStruct> CartesianMapStruct::ptr() const {
  void *ptr;
  cartesian_map_struct_get_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return CartesianMapTermStruct(ptr);
}
void CartesianMapStruct::set_ptr(const CartesianMapTermStruct &src) {
  cartesian_map_struct_set_ptr(fortran_ptr_, src.get_fortran_ptr());
}
std::complex<double> CylindricalMapTerm1Struct::e_coef() const {
  std::complex<double> c_value;
  cylindrical_map_term1_struct_get_e_coef(fortran_ptr_, &c_value);
  return c_value;
}
void CylindricalMapTerm1Struct::set_e_coef(std::complex<double> value) {
  cylindrical_map_term1_struct_set_e_coef(fortran_ptr_, value);
}
std::complex<double> CylindricalMapTerm1Struct::b_coef() const {
  std::complex<double> c_value;
  cylindrical_map_term1_struct_get_b_coef(fortran_ptr_, &c_value);
  return c_value;
}
void CylindricalMapTerm1Struct::set_b_coef(std::complex<double> value) {
  cylindrical_map_term1_struct_set_b_coef(fortran_ptr_, value);
}
std::string CylindricalMapTermStruct::file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, cylindrical_map_term_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void CylindricalMapTermStruct::set_file(const std::string &value) {
  cylindrical_map_term_struct_set_file(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
int CylindricalMapTermStruct::n_link() const {
  int value;
  cylindrical_map_term_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void CylindricalMapTermStruct::set_n_link(int value) {
  cylindrical_map_term_struct_set_n_link(fortran_ptr_, value);
}
CylindricalMapTerm1StructAlloc1D CylindricalMapTermStruct::term() const {
  return CylindricalMapTerm1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      cylindrical_map_term_struct_reallocate_term,
      cylindrical_map_term_struct_get_term_info
  );
}
int CylindricalMapStruct::m() const {
  int value;
  cylindrical_map_struct_get_m(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_m(int value) { cylindrical_map_struct_set_m(fortran_ptr_, value); }
int CylindricalMapStruct::harmonic() const {
  int value;
  cylindrical_map_struct_get_harmonic(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_harmonic(int value) {
  cylindrical_map_struct_set_harmonic(fortran_ptr_, value);
}
double CylindricalMapStruct::phi0_fieldmap() const {
  double value;
  cylindrical_map_struct_get_phi0_fieldmap(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_phi0_fieldmap(double value) {
  cylindrical_map_struct_set_phi0_fieldmap(fortran_ptr_, value);
}
double CylindricalMapStruct::theta0_azimuth() const {
  double value;
  cylindrical_map_struct_get_theta0_azimuth(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_theta0_azimuth(double value) {
  cylindrical_map_struct_set_theta0_azimuth(fortran_ptr_, value);
}
double CylindricalMapStruct::field_scale() const {
  double value;
  cylindrical_map_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_field_scale(double value) {
  cylindrical_map_struct_set_field_scale(fortran_ptr_, value);
}
int CylindricalMapStruct::master_parameter() const {
  int value;
  cylindrical_map_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_master_parameter(int value) {
  cylindrical_map_struct_set_master_parameter(fortran_ptr_, value);
}
int CylindricalMapStruct::ele_anchor_pt() const {
  int value;
  cylindrical_map_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_ele_anchor_pt(int value) {
  cylindrical_map_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
double CylindricalMapStruct::dz() const {
  double value;
  cylindrical_map_struct_get_dz(fortran_ptr_, &value);
  return value;
}
void CylindricalMapStruct::set_dz(double value) {
  cylindrical_map_struct_set_dz(fortran_ptr_, value);
}
FArray1D<double> CylindricalMapStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, cylindrical_map_struct_get_r0_info);
}
void CylindricalMapStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  cylindrical_map_struct_set_r0(fortran_ptr_, v.data(), shape);
}
std::optional<CylindricalMapTermStruct> CylindricalMapStruct::ptr() const {
  void *ptr;
  cylindrical_map_struct_get_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return CylindricalMapTermStruct(ptr);
}
void CylindricalMapStruct::set_ptr(const CylindricalMapTermStruct &src) {
  cylindrical_map_struct_set_ptr(fortran_ptr_, src.get_fortran_ptr());
}
FArray2D<std::complex<double>> BicubicCmplxCoefStruct::coef() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_,
      bicubic_cmplx_coef_struct_get_coef_info
  );
}
void BicubicCmplxCoefStruct::set_coef(const std::vector<std::vector<std::complex<double>>> &v) {
  ProxyHelpers::set_array_2d<std::complex<double>>(
      fortran_ptr_,
      bicubic_cmplx_coef_struct_set_coef,
      v
  );
}
FArray1D<int> BicubicCmplxCoefStruct::i_box() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, bicubic_cmplx_coef_struct_get_i_box_info);
}
void BicubicCmplxCoefStruct::set_i_box(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bicubic_cmplx_coef_struct_set_i_box(fortran_ptr_, v.data(), shape);
}
FArray3D<std::complex<double>> TricubicCmplxCoefStruct::coef() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_,
      tricubic_cmplx_coef_struct_get_coef_info
  );
}
void TricubicCmplxCoefStruct::set_coef(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &v
) {
  ProxyHelpers::set_array_3d<std::complex<double>>(
      fortran_ptr_,
      tricubic_cmplx_coef_struct_set_coef,
      v
  );
}
FArray1D<int> TricubicCmplxCoefStruct::i_box() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, tricubic_cmplx_coef_struct_get_i_box_info);
}
void TricubicCmplxCoefStruct::set_i_box(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tricubic_cmplx_coef_struct_set_i_box(fortran_ptr_, v.data(), shape);
}
FArray1D<std::complex<double>> GridFieldPt1Struct::E() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_,
      grid_field_pt1_struct_get_E_info
  );
}
void GridFieldPt1Struct::set_E(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  grid_field_pt1_struct_set_E(fortran_ptr_, v.data(), shape);
}
FArray1D<std::complex<double>> GridFieldPt1Struct::B() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_,
      grid_field_pt1_struct_get_B_info
  );
}
void GridFieldPt1Struct::set_B(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  grid_field_pt1_struct_set_B(fortran_ptr_, v.data(), shape);
}
std::string GridFieldPtStruct::file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, grid_field_pt_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void GridFieldPtStruct::set_file(const std::string &value) {
  grid_field_pt_struct_set_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int GridFieldPtStruct::n_link() const {
  int value;
  grid_field_pt_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void GridFieldPtStruct::set_n_link(int value) {
  grid_field_pt_struct_set_n_link(fortran_ptr_, value);
}
GridFieldPt1StructArray3D GridFieldPtStruct::pt() const {
  return ProxyHelpers::get_type_array_3d<GridFieldPt1StructArray3D>(
      fortran_ptr_,
      grid_field_pt_struct_get_pt_info
  );
}
int GridFieldStruct::geometry() const {
  int value;
  grid_field_struct_get_geometry(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_geometry(int value) {
  grid_field_struct_set_geometry(fortran_ptr_, value);
}
int GridFieldStruct::harmonic() const {
  int value;
  grid_field_struct_get_harmonic(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_harmonic(int value) {
  grid_field_struct_set_harmonic(fortran_ptr_, value);
}
double GridFieldStruct::phi0_fieldmap() const {
  double value;
  grid_field_struct_get_phi0_fieldmap(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_phi0_fieldmap(double value) {
  grid_field_struct_set_phi0_fieldmap(fortran_ptr_, value);
}
double GridFieldStruct::field_scale() const {
  double value;
  grid_field_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_field_scale(double value) {
  grid_field_struct_set_field_scale(fortran_ptr_, value);
}
int GridFieldStruct::field_type() const {
  int value;
  grid_field_struct_get_field_type(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_field_type(int value) {
  grid_field_struct_set_field_type(fortran_ptr_, value);
}
int GridFieldStruct::master_parameter() const {
  int value;
  grid_field_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_master_parameter(int value) {
  grid_field_struct_set_master_parameter(fortran_ptr_, value);
}
int GridFieldStruct::ele_anchor_pt() const {
  int value;
  grid_field_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_ele_anchor_pt(int value) {
  grid_field_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
int GridFieldStruct::interpolation_order() const {
  int value;
  grid_field_struct_get_interpolation_order(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_interpolation_order(int value) {
  grid_field_struct_set_interpolation_order(fortran_ptr_, value);
}
FArray1D<double> GridFieldStruct::dr() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, grid_field_struct_get_dr_info);
}
void GridFieldStruct::set_dr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  grid_field_struct_set_dr(fortran_ptr_, v.data(), shape);
}
FArray1D<double> GridFieldStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, grid_field_struct_get_r0_info);
}
void GridFieldStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  grid_field_struct_set_r0(fortran_ptr_, v.data(), shape);
}
bool GridFieldStruct::curved_ref_frame() const {
  bool value;
  grid_field_struct_get_curved_ref_frame(fortran_ptr_, &value);
  return value;
}
void GridFieldStruct::set_curved_ref_frame(bool value) {
  grid_field_struct_set_curved_ref_frame(fortran_ptr_, value);
}
std::optional<GridFieldPtStruct> GridFieldStruct::ptr() const {
  void *ptr;
  grid_field_struct_get_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return GridFieldPtStruct(ptr);
}
void GridFieldStruct::set_ptr(const GridFieldPtStruct &src) {
  grid_field_struct_set_ptr(fortran_ptr_, src.get_fortran_ptr());
}
BicubicCmplxCoefStructArray3D GridFieldStruct::bi_coef() const {
  return ProxyHelpers::get_type_array_3d<BicubicCmplxCoefStructArray3D>(
      fortran_ptr_,
      grid_field_struct_get_bi_coef_info
  );
}
TricubicCmplxCoefStructArray3D GridFieldStruct::tri_coef() const {
  return ProxyHelpers::get_type_array_3d<TricubicCmplxCoefStructArray3D>(
      fortran_ptr_,
      grid_field_struct_get_tri_coef_info
  );
}
FArray1D<double> FloorPositionStruct::r() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, floor_position_struct_get_r_info);
}
void FloorPositionStruct::set_r(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  floor_position_struct_set_r(fortran_ptr_, v.data(), shape);
}
FArray2D<double> FloorPositionStruct::w() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, floor_position_struct_get_w_info);
}
void FloorPositionStruct::set_w(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, floor_position_struct_set_w, v);
}
double FloorPositionStruct::theta() const {
  double value;
  floor_position_struct_get_theta(fortran_ptr_, &value);
  return value;
}
void FloorPositionStruct::set_theta(double value) {
  floor_position_struct_set_theta(fortran_ptr_, value);
}
double FloorPositionStruct::phi() const {
  double value;
  floor_position_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void FloorPositionStruct::set_phi(double value) {
  floor_position_struct_set_phi(fortran_ptr_, value);
}
double FloorPositionStruct::psi() const {
  double value;
  floor_position_struct_get_psi(fortran_ptr_, &value);
  return value;
}
void FloorPositionStruct::set_psi(double value) {
  floor_position_struct_set_psi(fortran_ptr_, value);
}
CoordStruct HighEnergySpaceChargeStruct::closed_orb() const {
  void *ptr;
  high_energy_space_charge_struct_get_closed_orb(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void HighEnergySpaceChargeStruct::set_closed_orb(const CoordStruct &src) {
  high_energy_space_charge_struct_set_closed_orb(fortran_ptr_, src.get_fortran_ptr());
}
double HighEnergySpaceChargeStruct::kick_const() const {
  double value;
  high_energy_space_charge_struct_get_kick_const(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_kick_const(double value) {
  high_energy_space_charge_struct_set_kick_const(fortran_ptr_, value);
}
double HighEnergySpaceChargeStruct::sig_x() const {
  double value;
  high_energy_space_charge_struct_get_sig_x(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_sig_x(double value) {
  high_energy_space_charge_struct_set_sig_x(fortran_ptr_, value);
}
double HighEnergySpaceChargeStruct::sig_y() const {
  double value;
  high_energy_space_charge_struct_get_sig_y(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_sig_y(double value) {
  high_energy_space_charge_struct_set_sig_y(fortran_ptr_, value);
}
double HighEnergySpaceChargeStruct::phi() const {
  double value;
  high_energy_space_charge_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_phi(double value) {
  high_energy_space_charge_struct_set_phi(fortran_ptr_, value);
}
double HighEnergySpaceChargeStruct::sin_phi() const {
  double value;
  high_energy_space_charge_struct_get_sin_phi(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_sin_phi(double value) {
  high_energy_space_charge_struct_set_sin_phi(fortran_ptr_, value);
}
double HighEnergySpaceChargeStruct::cos_phi() const {
  double value;
  high_energy_space_charge_struct_get_cos_phi(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_cos_phi(double value) {
  high_energy_space_charge_struct_set_cos_phi(fortran_ptr_, value);
}
double HighEnergySpaceChargeStruct::sig_z() const {
  double value;
  high_energy_space_charge_struct_get_sig_z(fortran_ptr_, &value);
  return value;
}
void HighEnergySpaceChargeStruct::set_sig_z(double value) {
  high_energy_space_charge_struct_set_sig_z(fortran_ptr_, value);
}
double XyDispStruct::eta() const {
  double value;
  xy_disp_struct_get_eta(fortran_ptr_, &value);
  return value;
}
void XyDispStruct::set_eta(double value) { xy_disp_struct_set_eta(fortran_ptr_, value); }
double XyDispStruct::etap() const {
  double value;
  xy_disp_struct_get_etap(fortran_ptr_, &value);
  return value;
}
void XyDispStruct::set_etap(double value) { xy_disp_struct_set_etap(fortran_ptr_, value); }
double XyDispStruct::deta_ds() const {
  double value;
  xy_disp_struct_get_deta_ds(fortran_ptr_, &value);
  return value;
}
void XyDispStruct::set_deta_ds(double value) { xy_disp_struct_set_deta_ds(fortran_ptr_, value); }
double XyDispStruct::sigma() const {
  double value;
  xy_disp_struct_get_sigma(fortran_ptr_, &value);
  return value;
}
void XyDispStruct::set_sigma(double value) { xy_disp_struct_set_sigma(fortran_ptr_, value); }
double XyDispStruct::deta_dpz() const {
  double value;
  xy_disp_struct_get_deta_dpz(fortran_ptr_, &value);
  return value;
}
void XyDispStruct::set_deta_dpz(double value) { xy_disp_struct_set_deta_dpz(fortran_ptr_, value); }
double XyDispStruct::detap_dpz() const {
  double value;
  xy_disp_struct_get_detap_dpz(fortran_ptr_, &value);
  return value;
}
void XyDispStruct::set_detap_dpz(double value) {
  xy_disp_struct_set_detap_dpz(fortran_ptr_, value);
}
double TwissStruct::beta() const {
  double value;
  twiss_struct_get_beta(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_beta(double value) { twiss_struct_set_beta(fortran_ptr_, value); }
double TwissStruct::alpha() const {
  double value;
  twiss_struct_get_alpha(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_alpha(double value) { twiss_struct_set_alpha(fortran_ptr_, value); }
double TwissStruct::gamma() const {
  double value;
  twiss_struct_get_gamma(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_gamma(double value) { twiss_struct_set_gamma(fortran_ptr_, value); }
double TwissStruct::phi() const {
  double value;
  twiss_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_phi(double value) { twiss_struct_set_phi(fortran_ptr_, value); }
double TwissStruct::eta() const {
  double value;
  twiss_struct_get_eta(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_eta(double value) { twiss_struct_set_eta(fortran_ptr_, value); }
double TwissStruct::etap() const {
  double value;
  twiss_struct_get_etap(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_etap(double value) { twiss_struct_set_etap(fortran_ptr_, value); }
double TwissStruct::deta_ds() const {
  double value;
  twiss_struct_get_deta_ds(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_deta_ds(double value) { twiss_struct_set_deta_ds(fortran_ptr_, value); }
double TwissStruct::sigma() const {
  double value;
  twiss_struct_get_sigma(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_sigma(double value) { twiss_struct_set_sigma(fortran_ptr_, value); }
double TwissStruct::sigma_p() const {
  double value;
  twiss_struct_get_sigma_p(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_sigma_p(double value) { twiss_struct_set_sigma_p(fortran_ptr_, value); }
double TwissStruct::emit() const {
  double value;
  twiss_struct_get_emit(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_emit(double value) { twiss_struct_set_emit(fortran_ptr_, value); }
double TwissStruct::norm_emit() const {
  double value;
  twiss_struct_get_norm_emit(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_norm_emit(double value) { twiss_struct_set_norm_emit(fortran_ptr_, value); }
double TwissStruct::chrom() const {
  double value;
  twiss_struct_get_chrom(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_chrom(double value) { twiss_struct_set_chrom(fortran_ptr_, value); }
double TwissStruct::dbeta_dpz() const {
  double value;
  twiss_struct_get_dbeta_dpz(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_dbeta_dpz(double value) { twiss_struct_set_dbeta_dpz(fortran_ptr_, value); }
double TwissStruct::dalpha_dpz() const {
  double value;
  twiss_struct_get_dalpha_dpz(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_dalpha_dpz(double value) { twiss_struct_set_dalpha_dpz(fortran_ptr_, value); }
double TwissStruct::deta_dpz() const {
  double value;
  twiss_struct_get_deta_dpz(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_deta_dpz(double value) { twiss_struct_set_deta_dpz(fortran_ptr_, value); }
double TwissStruct::detap_dpz() const {
  double value;
  twiss_struct_get_detap_dpz(fortran_ptr_, &value);
  return value;
}
void TwissStruct::set_detap_dpz(double value) { twiss_struct_set_detap_dpz(fortran_ptr_, value); }
FArray2D<double> Mode3Struct::v() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, mode3_struct_get_v_info);
}
void Mode3Struct::set_v(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, mode3_struct_set_v, v);
}
TwissStruct Mode3Struct::a() const {
  void *ptr;
  mode3_struct_get_a(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void Mode3Struct::set_a(const TwissStruct &src) {
  mode3_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct Mode3Struct::b() const {
  void *ptr;
  mode3_struct_get_b(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void Mode3Struct::set_b(const TwissStruct &src) {
  mode3_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct Mode3Struct::c() const {
  void *ptr;
  mode3_struct_get_c(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void Mode3Struct::set_c(const TwissStruct &src) {
  mode3_struct_set_c(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct Mode3Struct::x() const {
  void *ptr;
  mode3_struct_get_x(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void Mode3Struct::set_x(const TwissStruct &src) {
  mode3_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct Mode3Struct::y() const {
  void *ptr;
  mode3_struct_get_y(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void Mode3Struct::set_y(const TwissStruct &src) {
  mode3_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
int BookkeepingStateStruct::attributes() const {
  int value;
  bookkeeping_state_struct_get_attributes(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_attributes(int value) {
  bookkeeping_state_struct_set_attributes(fortran_ptr_, value);
}
int BookkeepingStateStruct::control() const {
  int value;
  bookkeeping_state_struct_get_control(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_control(int value) {
  bookkeeping_state_struct_set_control(fortran_ptr_, value);
}
int BookkeepingStateStruct::floor_position() const {
  int value;
  bookkeeping_state_struct_get_floor_position(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_floor_position(int value) {
  bookkeeping_state_struct_set_floor_position(fortran_ptr_, value);
}
int BookkeepingStateStruct::s_position() const {
  int value;
  bookkeeping_state_struct_get_s_position(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_s_position(int value) {
  bookkeeping_state_struct_set_s_position(fortran_ptr_, value);
}
int BookkeepingStateStruct::ref_energy() const {
  int value;
  bookkeeping_state_struct_get_ref_energy(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_ref_energy(int value) {
  bookkeeping_state_struct_set_ref_energy(fortran_ptr_, value);
}
int BookkeepingStateStruct::mat6() const {
  int value;
  bookkeeping_state_struct_get_mat6(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_mat6(int value) {
  bookkeeping_state_struct_set_mat6(fortran_ptr_, value);
}
int BookkeepingStateStruct::rad_int() const {
  int value;
  bookkeeping_state_struct_get_rad_int(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_rad_int(int value) {
  bookkeeping_state_struct_set_rad_int(fortran_ptr_, value);
}
int BookkeepingStateStruct::ptc() const {
  int value;
  bookkeeping_state_struct_get_ptc(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_ptc(int value) {
  bookkeeping_state_struct_set_ptc(fortran_ptr_, value);
}
bool BookkeepingStateStruct::has_misalign() const {
  bool value;
  bookkeeping_state_struct_get_has_misalign(fortran_ptr_, &value);
  return value;
}
void BookkeepingStateStruct::set_has_misalign(bool value) {
  bookkeeping_state_struct_set_has_misalign(fortran_ptr_, value);
}
FArray1D<double> RadMapStruct::ref_orb() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, rad_map_struct_get_ref_orb_info);
}
void RadMapStruct::set_ref_orb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  rad_map_struct_set_ref_orb(fortran_ptr_, v.data(), shape);
}
FArray2D<double> RadMapStruct::damp_dmat() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, rad_map_struct_get_damp_dmat_info);
}
void RadMapStruct::set_damp_dmat(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, rad_map_struct_set_damp_dmat, v);
}
FArray1D<double> RadMapStruct::xfer_damp_vec() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, rad_map_struct_get_xfer_damp_vec_info);
}
void RadMapStruct::set_xfer_damp_vec(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  rad_map_struct_set_xfer_damp_vec(fortran_ptr_, v.data(), shape);
}
FArray2D<double> RadMapStruct::xfer_damp_mat() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, rad_map_struct_get_xfer_damp_mat_info);
}
void RadMapStruct::set_xfer_damp_mat(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, rad_map_struct_set_xfer_damp_mat, v);
}
FArray2D<double> RadMapStruct::stoc_mat() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, rad_map_struct_get_stoc_mat_info);
}
void RadMapStruct::set_stoc_mat(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, rad_map_struct_set_stoc_mat, v);
}
RadMapStruct RadMapEleStruct::rm0() const {
  void *ptr;
  rad_map_ele_struct_get_rm0(fortran_ptr_, &ptr);
  return RadMapStruct(ptr);
}
void RadMapEleStruct::set_rm0(const RadMapStruct &src) {
  rad_map_ele_struct_set_rm0(fortran_ptr_, src.get_fortran_ptr());
}
RadMapStruct RadMapEleStruct::rm1() const {
  void *ptr;
  rad_map_ele_struct_get_rm1(fortran_ptr_, &ptr);
  return RadMapStruct(ptr);
}
void RadMapEleStruct::set_rm1(const RadMapStruct &src) {
  rad_map_ele_struct_set_rm1(fortran_ptr_, src.get_fortran_ptr());
}
bool RadMapEleStruct::stale() const {
  bool value;
  rad_map_ele_struct_get_stale(fortran_ptr_, &value);
  return value;
}
void RadMapEleStruct::set_stale(bool value) { rad_map_ele_struct_set_stale(fortran_ptr_, value); }
int GenGrad1Struct::m() const {
  int value;
  gen_grad1_struct_get_m(fortran_ptr_, &value);
  return value;
}
void GenGrad1Struct::set_m(int value) { gen_grad1_struct_set_m(fortran_ptr_, value); }
int GenGrad1Struct::sincos() const {
  int value;
  gen_grad1_struct_get_sincos(fortran_ptr_, &value);
  return value;
}
void GenGrad1Struct::set_sincos(int value) { gen_grad1_struct_set_sincos(fortran_ptr_, value); }
int GenGrad1Struct::n_deriv_max() const {
  int value;
  gen_grad1_struct_get_n_deriv_max(fortran_ptr_, &value);
  return value;
}
void GenGrad1Struct::set_n_deriv_max(int value) {
  gen_grad1_struct_set_n_deriv_max(fortran_ptr_, value);
}
FArray2D<double> GenGrad1Struct::deriv() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, gen_grad1_struct_get_deriv_info);
}
void GenGrad1Struct::set_deriv(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, gen_grad1_struct_set_deriv, v);
}
std::string GenGradMapStruct::file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, gen_grad_map_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void GenGradMapStruct::set_file(const std::string &value) {
  gen_grad_map_struct_set_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
GenGrad1StructAlloc1D GenGradMapStruct::gg() const {
  return GenGrad1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      gen_grad_map_struct_reallocate_gg,
      gen_grad_map_struct_get_gg_info
  );
}
int GenGradMapStruct::ele_anchor_pt() const {
  int value;
  gen_grad_map_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_ele_anchor_pt(int value) {
  gen_grad_map_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
int GenGradMapStruct::field_type() const {
  int value;
  gen_grad_map_struct_get_field_type(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_field_type(int value) {
  gen_grad_map_struct_set_field_type(fortran_ptr_, value);
}
int GenGradMapStruct::iz0() const {
  int value;
  gen_grad_map_struct_get_iz0(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_iz0(int value) { gen_grad_map_struct_set_iz0(fortran_ptr_, value); }
int GenGradMapStruct::iz1() const {
  int value;
  gen_grad_map_struct_get_iz1(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_iz1(int value) { gen_grad_map_struct_set_iz1(fortran_ptr_, value); }
double GenGradMapStruct::dz() const {
  double value;
  gen_grad_map_struct_get_dz(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_dz(double value) { gen_grad_map_struct_set_dz(fortran_ptr_, value); }
FArray1D<double> GenGradMapStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, gen_grad_map_struct_get_r0_info);
}
void GenGradMapStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  gen_grad_map_struct_set_r0(fortran_ptr_, v.data(), shape);
}
double GenGradMapStruct::field_scale() const {
  double value;
  gen_grad_map_struct_get_field_scale(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_field_scale(double value) {
  gen_grad_map_struct_set_field_scale(fortran_ptr_, value);
}
int GenGradMapStruct::master_parameter() const {
  int value;
  gen_grad_map_struct_get_master_parameter(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_master_parameter(int value) {
  gen_grad_map_struct_set_master_parameter(fortran_ptr_, value);
}
bool GenGradMapStruct::curved_ref_frame() const {
  bool value;
  gen_grad_map_struct_get_curved_ref_frame(fortran_ptr_, &value);
  return value;
}
void GenGradMapStruct::set_curved_ref_frame(bool value) {
  gen_grad_map_struct_set_curved_ref_frame(fortran_ptr_, value);
}
double SurfaceSegmentedPtStruct::x0() const {
  double value;
  surface_segmented_pt_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtStruct::set_x0(double value) {
  surface_segmented_pt_struct_set_x0(fortran_ptr_, value);
}
double SurfaceSegmentedPtStruct::y0() const {
  double value;
  surface_segmented_pt_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtStruct::set_y0(double value) {
  surface_segmented_pt_struct_set_y0(fortran_ptr_, value);
}
double SurfaceSegmentedPtStruct::z0() const {
  double value;
  surface_segmented_pt_struct_get_z0(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtStruct::set_z0(double value) {
  surface_segmented_pt_struct_set_z0(fortran_ptr_, value);
}
double SurfaceSegmentedPtStruct::dz_dx() const {
  double value;
  surface_segmented_pt_struct_get_dz_dx(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtStruct::set_dz_dx(double value) {
  surface_segmented_pt_struct_set_dz_dx(fortran_ptr_, value);
}
double SurfaceSegmentedPtStruct::dz_dy() const {
  double value;
  surface_segmented_pt_struct_get_dz_dy(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedPtStruct::set_dz_dy(double value) {
  surface_segmented_pt_struct_set_dz_dy(fortran_ptr_, value);
}
bool SurfaceSegmentedStruct::active() const {
  bool value;
  surface_segmented_struct_get_active(fortran_ptr_, &value);
  return value;
}
void SurfaceSegmentedStruct::set_active(bool value) {
  surface_segmented_struct_set_active(fortran_ptr_, value);
}
FArray1D<double> SurfaceSegmentedStruct::dr() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, surface_segmented_struct_get_dr_info);
}
void SurfaceSegmentedStruct::set_dr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_segmented_struct_set_dr(fortran_ptr_, v.data(), shape);
}
FArray1D<double> SurfaceSegmentedStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, surface_segmented_struct_get_r0_info);
}
void SurfaceSegmentedStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_segmented_struct_set_r0(fortran_ptr_, v.data(), shape);
}
SurfaceSegmentedPtStructArray2D SurfaceSegmentedStruct::pt() const {
  return ProxyHelpers::get_type_array_2d<SurfaceSegmentedPtStructArray2D>(
      fortran_ptr_,
      surface_segmented_struct_get_pt_info
  );
}
double SurfaceHMisalignPtStruct::x0() const {
  double value;
  surface_h_misalign_pt_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtStruct::set_x0(double value) {
  surface_h_misalign_pt_struct_set_x0(fortran_ptr_, value);
}
double SurfaceHMisalignPtStruct::y0() const {
  double value;
  surface_h_misalign_pt_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtStruct::set_y0(double value) {
  surface_h_misalign_pt_struct_set_y0(fortran_ptr_, value);
}
double SurfaceHMisalignPtStruct::rot_y() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_y(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtStruct::set_rot_y(double value) {
  surface_h_misalign_pt_struct_set_rot_y(fortran_ptr_, value);
}
double SurfaceHMisalignPtStruct::rot_t() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_t(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtStruct::set_rot_t(double value) {
  surface_h_misalign_pt_struct_set_rot_t(fortran_ptr_, value);
}
double SurfaceHMisalignPtStruct::rot_y_rms() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_y_rms(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtStruct::set_rot_y_rms(double value) {
  surface_h_misalign_pt_struct_set_rot_y_rms(fortran_ptr_, value);
}
double SurfaceHMisalignPtStruct::rot_t_rms() const {
  double value;
  surface_h_misalign_pt_struct_get_rot_t_rms(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignPtStruct::set_rot_t_rms(double value) {
  surface_h_misalign_pt_struct_set_rot_t_rms(fortran_ptr_, value);
}
bool SurfaceHMisalignStruct::active() const {
  bool value;
  surface_h_misalign_struct_get_active(fortran_ptr_, &value);
  return value;
}
void SurfaceHMisalignStruct::set_active(bool value) {
  surface_h_misalign_struct_set_active(fortran_ptr_, value);
}
FArray1D<double> SurfaceHMisalignStruct::dr() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, surface_h_misalign_struct_get_dr_info);
}
void SurfaceHMisalignStruct::set_dr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_h_misalign_struct_set_dr(fortran_ptr_, v.data(), shape);
}
FArray1D<double> SurfaceHMisalignStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, surface_h_misalign_struct_get_r0_info);
}
void SurfaceHMisalignStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_h_misalign_struct_set_r0(fortran_ptr_, v.data(), shape);
}
SurfaceHMisalignPtStructArray2D SurfaceHMisalignStruct::pt() const {
  return ProxyHelpers::get_type_array_2d<SurfaceHMisalignPtStructArray2D>(
      fortran_ptr_,
      surface_h_misalign_struct_get_pt_info
  );
}
double SurfaceDisplacementPtStruct::x0() const {
  double value;
  surface_displacement_pt_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtStruct::set_x0(double value) {
  surface_displacement_pt_struct_set_x0(fortran_ptr_, value);
}
double SurfaceDisplacementPtStruct::y0() const {
  double value;
  surface_displacement_pt_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtStruct::set_y0(double value) {
  surface_displacement_pt_struct_set_y0(fortran_ptr_, value);
}
double SurfaceDisplacementPtStruct::z0() const {
  double value;
  surface_displacement_pt_struct_get_z0(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtStruct::set_z0(double value) {
  surface_displacement_pt_struct_set_z0(fortran_ptr_, value);
}
double SurfaceDisplacementPtStruct::dz_dx() const {
  double value;
  surface_displacement_pt_struct_get_dz_dx(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtStruct::set_dz_dx(double value) {
  surface_displacement_pt_struct_set_dz_dx(fortran_ptr_, value);
}
double SurfaceDisplacementPtStruct::dz_dy() const {
  double value;
  surface_displacement_pt_struct_get_dz_dy(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtStruct::set_dz_dy(double value) {
  surface_displacement_pt_struct_set_dz_dy(fortran_ptr_, value);
}
double SurfaceDisplacementPtStruct::d2z_dxdy() const {
  double value;
  surface_displacement_pt_struct_get_d2z_dxdy(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementPtStruct::set_d2z_dxdy(double value) {
  surface_displacement_pt_struct_set_d2z_dxdy(fortran_ptr_, value);
}
bool SurfaceDisplacementStruct::active() const {
  bool value;
  surface_displacement_struct_get_active(fortran_ptr_, &value);
  return value;
}
void SurfaceDisplacementStruct::set_active(bool value) {
  surface_displacement_struct_set_active(fortran_ptr_, value);
}
FArray1D<double> SurfaceDisplacementStruct::dr() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, surface_displacement_struct_get_dr_info);
}
void SurfaceDisplacementStruct::set_dr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_displacement_struct_set_dr(fortran_ptr_, v.data(), shape);
}
FArray1D<double> SurfaceDisplacementStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, surface_displacement_struct_get_r0_info);
}
void SurfaceDisplacementStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_displacement_struct_set_r0(fortran_ptr_, v.data(), shape);
}
SurfaceDisplacementPtStructArray2D SurfaceDisplacementStruct::pt() const {
  return ProxyHelpers::get_type_array_2d<SurfaceDisplacementPtStructArray2D>(
      fortran_ptr_,
      surface_displacement_struct_get_pt_info
  );
}
FArray1D<double> TargetPointStruct::r() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, target_point_struct_get_r_info);
}
void TargetPointStruct::set_r(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  target_point_struct_set_r(fortran_ptr_, v.data(), shape);
}
FArray2D<double> SurfaceCurvatureStruct::xy() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, surface_curvature_struct_get_xy_info);
}
void SurfaceCurvatureStruct::set_xy(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, surface_curvature_struct_set_xy, v);
}
double SurfaceCurvatureStruct::spherical() const {
  double value;
  surface_curvature_struct_get_spherical(fortran_ptr_, &value);
  return value;
}
void SurfaceCurvatureStruct::set_spherical(double value) {
  surface_curvature_struct_set_spherical(fortran_ptr_, value);
}
FArray1D<double> SurfaceCurvatureStruct::elliptical() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      surface_curvature_struct_get_elliptical_info
  );
}
void SurfaceCurvatureStruct::set_elliptical(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  surface_curvature_struct_set_elliptical(fortran_ptr_, v.data(), shape);
}
bool SurfaceCurvatureStruct::has_curvature() const {
  bool value;
  surface_curvature_struct_get_has_curvature(fortran_ptr_, &value);
  return value;
}
void SurfaceCurvatureStruct::set_has_curvature(bool value) {
  surface_curvature_struct_set_has_curvature(fortran_ptr_, value);
}
int PhotonTargetStruct::type() const {
  int value;
  photon_target_struct_get_type(fortran_ptr_, &value);
  return value;
}
void PhotonTargetStruct::set_type(int value) { photon_target_struct_set_type(fortran_ptr_, value); }
int PhotonTargetStruct::n_corner() const {
  int value;
  photon_target_struct_get_n_corner(fortran_ptr_, &value);
  return value;
}
void PhotonTargetStruct::set_n_corner(int value) {
  photon_target_struct_set_n_corner(fortran_ptr_, value);
}
LatEleLocStruct PhotonTargetStruct::ele_loc() const {
  void *ptr;
  photon_target_struct_get_ele_loc(fortran_ptr_, &ptr);
  return LatEleLocStruct(ptr);
}
void PhotonTargetStruct::set_ele_loc(const LatEleLocStruct &src) {
  photon_target_struct_set_ele_loc(fortran_ptr_, src.get_fortran_ptr());
}
TargetPointStructArray1D PhotonTargetStruct::corner() const {
  return ProxyHelpers::get_type_array_1d<TargetPointStructArray1D>(
      fortran_ptr_,
      photon_target_struct_get_corner_info
  );
}
TargetPointStruct PhotonTargetStruct::center() const {
  void *ptr;
  photon_target_struct_get_center(fortran_ptr_, &ptr);
  return TargetPointStruct(ptr);
}
void PhotonTargetStruct::set_center(const TargetPointStruct &src) {
  photon_target_struct_set_center(fortran_ptr_, src.get_fortran_ptr());
}
std::complex<double> PhotonMaterialStruct::f0_m1() const {
  std::complex<double> c_value;
  photon_material_struct_get_f0_m1(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialStruct::set_f0_m1(std::complex<double> value) {
  photon_material_struct_set_f0_m1(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialStruct::f0_m2() const {
  std::complex<double> c_value;
  photon_material_struct_get_f0_m2(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialStruct::set_f0_m2(std::complex<double> value) {
  photon_material_struct_set_f0_m2(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialStruct::f_0() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_0(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialStruct::set_f_0(std::complex<double> value) {
  photon_material_struct_set_f_0(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialStruct::f_h() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_h(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialStruct::set_f_h(std::complex<double> value) {
  photon_material_struct_set_f_h(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialStruct::f_hbar() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_hbar(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialStruct::set_f_hbar(std::complex<double> value) {
  photon_material_struct_set_f_hbar(fortran_ptr_, value);
}
std::complex<double> PhotonMaterialStruct::f_hkl() const {
  std::complex<double> c_value;
  photon_material_struct_get_f_hkl(fortran_ptr_, &c_value);
  return c_value;
}
void PhotonMaterialStruct::set_f_hkl(std::complex<double> value) {
  photon_material_struct_set_f_hkl(fortran_ptr_, value);
}
FArray1D<double> PhotonMaterialStruct::h_norm() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, photon_material_struct_get_h_norm_info);
}
void PhotonMaterialStruct::set_h_norm(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_material_struct_set_h_norm(fortran_ptr_, v.data(), shape);
}
FArray1D<double> PhotonMaterialStruct::l_ref() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, photon_material_struct_get_l_ref_info);
}
void PhotonMaterialStruct::set_l_ref(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_material_struct_set_l_ref(fortran_ptr_, v.data(), shape);
}
int64_t PixelPtStruct::n_photon() const {
  int64_t value;
  pixel_pt_struct_get_n_photon(fortran_ptr_, &value);
  return value;
}
void PixelPtStruct::set_n_photon(int64_t value) {
  pixel_pt_struct_set_n_photon(fortran_ptr_, value);
}
std::complex<double> PixelPtStruct::E_x() const {
  std::complex<double> c_value;
  pixel_pt_struct_get_E_x(fortran_ptr_, &c_value);
  return c_value;
}
void PixelPtStruct::set_E_x(std::complex<double> value) {
  pixel_pt_struct_set_E_x(fortran_ptr_, value);
}
std::complex<double> PixelPtStruct::E_y() const {
  std::complex<double> c_value;
  pixel_pt_struct_get_E_y(fortran_ptr_, &c_value);
  return c_value;
}
void PixelPtStruct::set_E_y(std::complex<double> value) {
  pixel_pt_struct_set_E_y(fortran_ptr_, value);
}
double PixelPtStruct::intensity_x() const {
  double value;
  pixel_pt_struct_get_intensity_x(fortran_ptr_, &value);
  return value;
}
void PixelPtStruct::set_intensity_x(double value) {
  pixel_pt_struct_set_intensity_x(fortran_ptr_, value);
}
double PixelPtStruct::intensity_y() const {
  double value;
  pixel_pt_struct_get_intensity_y(fortran_ptr_, &value);
  return value;
}
void PixelPtStruct::set_intensity_y(double value) {
  pixel_pt_struct_set_intensity_y(fortran_ptr_, value);
}
double PixelPtStruct::intensity() const {
  double value;
  pixel_pt_struct_get_intensity(fortran_ptr_, &value);
  return value;
}
void PixelPtStruct::set_intensity(double value) {
  pixel_pt_struct_set_intensity(fortran_ptr_, value);
}
FArray1D<double> PixelPtStruct::orbit() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, pixel_pt_struct_get_orbit_info);
}
void PixelPtStruct::set_orbit(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  pixel_pt_struct_set_orbit(fortran_ptr_, v.data(), shape);
}
FArray1D<double> PixelPtStruct::orbit_rms() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, pixel_pt_struct_get_orbit_rms_info);
}
void PixelPtStruct::set_orbit_rms(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  pixel_pt_struct_set_orbit_rms(fortran_ptr_, v.data(), shape);
}
FArray1D<double> PixelPtStruct::init_orbit() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, pixel_pt_struct_get_init_orbit_info);
}
void PixelPtStruct::set_init_orbit(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  pixel_pt_struct_set_init_orbit(fortran_ptr_, v.data(), shape);
}
FArray1D<double> PixelPtStruct::init_orbit_rms() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, pixel_pt_struct_get_init_orbit_rms_info);
}
void PixelPtStruct::set_init_orbit_rms(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  pixel_pt_struct_set_init_orbit_rms(fortran_ptr_, v.data(), shape);
}
FArray1D<double> PixelDetecStruct::dr() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, pixel_detec_struct_get_dr_info);
}
void PixelDetecStruct::set_dr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  pixel_detec_struct_set_dr(fortran_ptr_, v.data(), shape);
}
FArray1D<double> PixelDetecStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, pixel_detec_struct_get_r0_info);
}
void PixelDetecStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  pixel_detec_struct_set_r0(fortran_ptr_, v.data(), shape);
}
int64_t PixelDetecStruct::n_track_tot() const {
  int64_t value;
  pixel_detec_struct_get_n_track_tot(fortran_ptr_, &value);
  return value;
}
void PixelDetecStruct::set_n_track_tot(int64_t value) {
  pixel_detec_struct_set_n_track_tot(fortran_ptr_, value);
}
int64_t PixelDetecStruct::n_hit_detec() const {
  int64_t value;
  pixel_detec_struct_get_n_hit_detec(fortran_ptr_, &value);
  return value;
}
void PixelDetecStruct::set_n_hit_detec(int64_t value) {
  pixel_detec_struct_set_n_hit_detec(fortran_ptr_, value);
}
int64_t PixelDetecStruct::n_hit_pixel() const {
  int64_t value;
  pixel_detec_struct_get_n_hit_pixel(fortran_ptr_, &value);
  return value;
}
void PixelDetecStruct::set_n_hit_pixel(int64_t value) {
  pixel_detec_struct_set_n_hit_pixel(fortran_ptr_, value);
}
PixelPtStructArray2D PixelDetecStruct::pt() const {
  return ProxyHelpers::get_type_array_2d<PixelPtStructArray2D>(
      fortran_ptr_,
      pixel_detec_struct_get_pt_info
  );
}
SurfaceCurvatureStruct PhotonElementStruct::curvature() const {
  void *ptr;
  photon_element_struct_get_curvature(fortran_ptr_, &ptr);
  return SurfaceCurvatureStruct(ptr);
}
void PhotonElementStruct::set_curvature(const SurfaceCurvatureStruct &src) {
  photon_element_struct_set_curvature(fortran_ptr_, src.get_fortran_ptr());
}
PhotonTargetStruct PhotonElementStruct::target() const {
  void *ptr;
  photon_element_struct_get_target(fortran_ptr_, &ptr);
  return PhotonTargetStruct(ptr);
}
void PhotonElementStruct::set_target(const PhotonTargetStruct &src) {
  photon_element_struct_set_target(fortran_ptr_, src.get_fortran_ptr());
}
PhotonMaterialStruct PhotonElementStruct::material() const {
  void *ptr;
  photon_element_struct_get_material(fortran_ptr_, &ptr);
  return PhotonMaterialStruct(ptr);
}
void PhotonElementStruct::set_material(const PhotonMaterialStruct &src) {
  photon_element_struct_set_material(fortran_ptr_, src.get_fortran_ptr());
}
SurfaceSegmentedStruct PhotonElementStruct::segmented() const {
  void *ptr;
  photon_element_struct_get_segmented(fortran_ptr_, &ptr);
  return SurfaceSegmentedStruct(ptr);
}
void PhotonElementStruct::set_segmented(const SurfaceSegmentedStruct &src) {
  photon_element_struct_set_segmented(fortran_ptr_, src.get_fortran_ptr());
}
SurfaceHMisalignStruct PhotonElementStruct::h_misalign() const {
  void *ptr;
  photon_element_struct_get_h_misalign(fortran_ptr_, &ptr);
  return SurfaceHMisalignStruct(ptr);
}
void PhotonElementStruct::set_h_misalign(const SurfaceHMisalignStruct &src) {
  photon_element_struct_set_h_misalign(fortran_ptr_, src.get_fortran_ptr());
}
SurfaceDisplacementStruct PhotonElementStruct::displacement() const {
  void *ptr;
  photon_element_struct_get_displacement(fortran_ptr_, &ptr);
  return SurfaceDisplacementStruct(ptr);
}
void PhotonElementStruct::set_displacement(const SurfaceDisplacementStruct &src) {
  photon_element_struct_set_displacement(fortran_ptr_, src.get_fortran_ptr());
}
PixelDetecStruct PhotonElementStruct::pixel() const {
  void *ptr;
  photon_element_struct_get_pixel(fortran_ptr_, &ptr);
  return PixelDetecStruct(ptr);
}
void PhotonElementStruct::set_pixel(const PixelDetecStruct &src) {
  photon_element_struct_set_pixel(fortran_ptr_, src.get_fortran_ptr());
}
int PhotonElementStruct::reflectivity_table_type() const {
  int value;
  photon_element_struct_get_reflectivity_table_type(fortran_ptr_, &value);
  return value;
}
void PhotonElementStruct::set_reflectivity_table_type(int value) {
  photon_element_struct_set_reflectivity_table_type(fortran_ptr_, value);
}
PhotonReflectTableStruct PhotonElementStruct::reflectivity_table_sigma() const {
  void *ptr;
  photon_element_struct_get_reflectivity_table_sigma(fortran_ptr_, &ptr);
  return PhotonReflectTableStruct(ptr);
}
void PhotonElementStruct::set_reflectivity_table_sigma(const PhotonReflectTableStruct &src) {
  photon_element_struct_set_reflectivity_table_sigma(fortran_ptr_, src.get_fortran_ptr());
}
PhotonReflectTableStruct PhotonElementStruct::reflectivity_table_pi() const {
  void *ptr;
  photon_element_struct_get_reflectivity_table_pi(fortran_ptr_, &ptr);
  return PhotonReflectTableStruct(ptr);
}
void PhotonElementStruct::set_reflectivity_table_pi(const PhotonReflectTableStruct &src) {
  photon_element_struct_set_reflectivity_table_pi(fortran_ptr_, src.get_fortran_ptr());
}
SplineStructAlloc1D PhotonElementStruct::init_energy_prob() const {
  return SplineStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_element_struct_reallocate_init_energy_prob,
      photon_element_struct_get_init_energy_prob_info
  );
}
RealAlloc1D PhotonElementStruct::integrated_init_energy_prob() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      photon_element_struct_reallocate_integrated_init_energy_prob,
      photon_element_struct_get_integrated_init_energy_prob_info
  );
}
void PhotonElementStruct::set_integrated_init_energy_prob(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  photon_element_struct_set_integrated_init_energy_prob(fortran_ptr_, v.data(), shape);
}
double Wall3dVertexStruct::x() const {
  double value;
  wall3d_vertex_struct_get_x(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_x(double value) { wall3d_vertex_struct_set_x(fortran_ptr_, value); }
double Wall3dVertexStruct::y() const {
  double value;
  wall3d_vertex_struct_get_y(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_y(double value) { wall3d_vertex_struct_set_y(fortran_ptr_, value); }
double Wall3dVertexStruct::radius_x() const {
  double value;
  wall3d_vertex_struct_get_radius_x(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_radius_x(double value) {
  wall3d_vertex_struct_set_radius_x(fortran_ptr_, value);
}
double Wall3dVertexStruct::radius_y() const {
  double value;
  wall3d_vertex_struct_get_radius_y(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_radius_y(double value) {
  wall3d_vertex_struct_set_radius_y(fortran_ptr_, value);
}
double Wall3dVertexStruct::tilt() const {
  double value;
  wall3d_vertex_struct_get_tilt(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_tilt(double value) {
  wall3d_vertex_struct_set_tilt(fortran_ptr_, value);
}
double Wall3dVertexStruct::angle() const {
  double value;
  wall3d_vertex_struct_get_angle(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_angle(double value) {
  wall3d_vertex_struct_set_angle(fortran_ptr_, value);
}
double Wall3dVertexStruct::x0() const {
  double value;
  wall3d_vertex_struct_get_x0(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_x0(double value) { wall3d_vertex_struct_set_x0(fortran_ptr_, value); }
double Wall3dVertexStruct::y0() const {
  double value;
  wall3d_vertex_struct_get_y0(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_y0(double value) { wall3d_vertex_struct_set_y0(fortran_ptr_, value); }
int Wall3dVertexStruct::type() const {
  int value;
  wall3d_vertex_struct_get_type(fortran_ptr_, &value);
  return value;
}
void Wall3dVertexStruct::set_type(int value) { wall3d_vertex_struct_set_type(fortran_ptr_, value); }
std::string Wall3dSectionStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, wall3d_section_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dSectionStruct::set_name(const std::string &value) {
  wall3d_section_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string Wall3dSectionStruct::material() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, wall3d_section_struct_get_material_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dSectionStruct::set_material(const std::string &value) {
  wall3d_section_struct_set_material(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
Wall3dVertexStructAlloc1D Wall3dSectionStruct::v() const {
  return Wall3dVertexStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wall3d_section_struct_reallocate_v,
      wall3d_section_struct_get_v_info
  );
}
std::optional<PhotonReflectSurfaceStruct> Wall3dSectionStruct::surface() const {
  void *ptr;
  wall3d_section_struct_get_surface(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return PhotonReflectSurfaceStruct(ptr);
}
void Wall3dSectionStruct::set_surface(const PhotonReflectSurfaceStruct &src) {
  wall3d_section_struct_set_surface(fortran_ptr_, src.get_fortran_ptr());
}
int Wall3dSectionStruct::type() const {
  int value;
  wall3d_section_struct_get_type(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_type(int value) {
  wall3d_section_struct_set_type(fortran_ptr_, value);
}
int Wall3dSectionStruct::n_vertex_input() const {
  int value;
  wall3d_section_struct_get_n_vertex_input(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_n_vertex_input(int value) {
  wall3d_section_struct_set_n_vertex_input(fortran_ptr_, value);
}
int Wall3dSectionStruct::ix_ele() const {
  int value;
  wall3d_section_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_ix_ele(int value) {
  wall3d_section_struct_set_ix_ele(fortran_ptr_, value);
}
int Wall3dSectionStruct::ix_branch() const {
  int value;
  wall3d_section_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_ix_branch(int value) {
  wall3d_section_struct_set_ix_branch(fortran_ptr_, value);
}
int Wall3dSectionStruct::vertices_state() const {
  int value;
  wall3d_section_struct_get_vertices_state(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_vertices_state(int value) {
  wall3d_section_struct_set_vertices_state(fortran_ptr_, value);
}
bool Wall3dSectionStruct::patch_in_region() const {
  bool value;
  wall3d_section_struct_get_patch_in_region(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_patch_in_region(bool value) {
  wall3d_section_struct_set_patch_in_region(fortran_ptr_, value);
}
double Wall3dSectionStruct::thickness() const {
  double value;
  wall3d_section_struct_get_thickness(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_thickness(double value) {
  wall3d_section_struct_set_thickness(fortran_ptr_, value);
}
double Wall3dSectionStruct::s() const {
  double value;
  wall3d_section_struct_get_s(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_s(double value) { wall3d_section_struct_set_s(fortran_ptr_, value); }
FArray1D<double> Wall3dSectionStruct::r0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, wall3d_section_struct_get_r0_info);
}
void Wall3dSectionStruct::set_r0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wall3d_section_struct_set_r0(fortran_ptr_, v.data(), shape);
}
double Wall3dSectionStruct::dx0_ds() const {
  double value;
  wall3d_section_struct_get_dx0_ds(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_dx0_ds(double value) {
  wall3d_section_struct_set_dx0_ds(fortran_ptr_, value);
}
double Wall3dSectionStruct::dy0_ds() const {
  double value;
  wall3d_section_struct_get_dy0_ds(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_dy0_ds(double value) {
  wall3d_section_struct_set_dy0_ds(fortran_ptr_, value);
}
FArray1D<double> Wall3dSectionStruct::x0_coef() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, wall3d_section_struct_get_x0_coef_info);
}
void Wall3dSectionStruct::set_x0_coef(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wall3d_section_struct_set_x0_coef(fortran_ptr_, v.data(), shape);
}
FArray1D<double> Wall3dSectionStruct::y0_coef() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, wall3d_section_struct_get_y0_coef_info);
}
void Wall3dSectionStruct::set_y0_coef(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wall3d_section_struct_set_y0_coef(fortran_ptr_, v.data(), shape);
}
double Wall3dSectionStruct::dr_ds() const {
  double value;
  wall3d_section_struct_get_dr_ds(fortran_ptr_, &value);
  return value;
}
void Wall3dSectionStruct::set_dr_ds(double value) {
  wall3d_section_struct_set_dr_ds(fortran_ptr_, value);
}
FArray1D<double> Wall3dSectionStruct::p1_coef() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, wall3d_section_struct_get_p1_coef_info);
}
void Wall3dSectionStruct::set_p1_coef(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wall3d_section_struct_set_p1_coef(fortran_ptr_, v.data(), shape);
}
FArray1D<double> Wall3dSectionStruct::p2_coef() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, wall3d_section_struct_get_p2_coef_info);
}
void Wall3dSectionStruct::set_p2_coef(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  wall3d_section_struct_set_p2_coef(fortran_ptr_, v.data(), shape);
}
std::string Wall3dStruct::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, wall3d_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dStruct::set_name(const std::string &value) {
  wall3d_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int Wall3dStruct::type() const {
  int value;
  wall3d_struct_get_type(fortran_ptr_, &value);
  return value;
}
void Wall3dStruct::set_type(int value) { wall3d_struct_set_type(fortran_ptr_, value); }
int Wall3dStruct::ix_wall3d() const {
  int value;
  wall3d_struct_get_ix_wall3d(fortran_ptr_, &value);
  return value;
}
void Wall3dStruct::set_ix_wall3d(int value) { wall3d_struct_set_ix_wall3d(fortran_ptr_, value); }
int Wall3dStruct::n_link() const {
  int value;
  wall3d_struct_get_n_link(fortran_ptr_, &value);
  return value;
}
void Wall3dStruct::set_n_link(int value) { wall3d_struct_set_n_link(fortran_ptr_, value); }
double Wall3dStruct::thickness() const {
  double value;
  wall3d_struct_get_thickness(fortran_ptr_, &value);
  return value;
}
void Wall3dStruct::set_thickness(double value) { wall3d_struct_set_thickness(fortran_ptr_, value); }
std::string Wall3dStruct::clear_material() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, wall3d_struct_get_clear_material_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dStruct::set_clear_material(const std::string &value) {
  wall3d_struct_set_clear_material(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string Wall3dStruct::opaque_material() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, wall3d_struct_get_opaque_material_info);
  return std::string(arr.data(), arr.size());
}
void Wall3dStruct::set_opaque_material(const std::string &value) {
  wall3d_struct_set_opaque_material(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool Wall3dStruct::superimpose() const {
  bool value;
  wall3d_struct_get_superimpose(fortran_ptr_, &value);
  return value;
}
void Wall3dStruct::set_superimpose(bool value) {
  wall3d_struct_set_superimpose(fortran_ptr_, value);
}
int Wall3dStruct::ele_anchor_pt() const {
  int value;
  wall3d_struct_get_ele_anchor_pt(fortran_ptr_, &value);
  return value;
}
void Wall3dStruct::set_ele_anchor_pt(int value) {
  wall3d_struct_set_ele_anchor_pt(fortran_ptr_, value);
}
Wall3dSectionStructAlloc1D Wall3dStruct::section() const {
  return Wall3dSectionStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      wall3d_struct_reallocate_section,
      wall3d_struct_get_section_info
  );
}
int RamperLordStruct::ix_ele() const {
  int value;
  ramper_lord_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void RamperLordStruct::set_ix_ele(int value) { ramper_lord_struct_set_ix_ele(fortran_ptr_, value); }
int RamperLordStruct::ix_con() const {
  int value;
  ramper_lord_struct_get_ix_con(fortran_ptr_, &value);
  return value;
}
void RamperLordStruct::set_ix_con(int value) { ramper_lord_struct_set_ix_con(fortran_ptr_, value); }
std::optional<double> RamperLordStruct::attrib_ptr() const {
  double value;
  bool is_valid;
  ramper_lord_struct_get_attrib_ptr(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void RamperLordStruct::set_attrib_ptr(double value) {
  ramper_lord_struct_set_attrib_ptr(fortran_ptr_, value);
}
double ControlStruct::value() const {
  double value;
  control_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ControlStruct::set_value(double value) { control_struct_set_value(fortran_ptr_, value); }
RealAlloc1D ControlStruct::y_knot() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      control_struct_reallocate_y_knot,
      control_struct_get_y_knot_info
  );
}
void ControlStruct::set_y_knot(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  control_struct_set_y_knot(fortran_ptr_, v.data(), shape);
}
ExpressionAtomStructAlloc1D ControlStruct::stack() const {
  return ExpressionAtomStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      control_struct_reallocate_stack,
      control_struct_get_stack_info
  );
}
LatEleLocStruct ControlStruct::slave() const {
  void *ptr;
  control_struct_get_slave(fortran_ptr_, &ptr);
  return LatEleLocStruct(ptr);
}
void ControlStruct::set_slave(const LatEleLocStruct &src) {
  control_struct_set_slave(fortran_ptr_, src.get_fortran_ptr());
}
LatEleLocStruct ControlStruct::lord() const {
  void *ptr;
  control_struct_get_lord(fortran_ptr_, &ptr);
  return LatEleLocStruct(ptr);
}
void ControlStruct::set_lord(const LatEleLocStruct &src) {
  control_struct_set_lord(fortran_ptr_, src.get_fortran_ptr());
}
std::string ControlStruct::slave_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, control_struct_get_slave_name_info);
  return std::string(arr.data(), arr.size());
}
void ControlStruct::set_slave_name(const std::string &value) {
  control_struct_set_slave_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string ControlStruct::attribute() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, control_struct_get_attribute_info);
  return std::string(arr.data(), arr.size());
}
void ControlStruct::set_attribute(const std::string &value) {
  control_struct_set_attribute(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int ControlStruct::ix_attrib() const {
  int value;
  control_struct_get_ix_attrib(fortran_ptr_, &value);
  return value;
}
void ControlStruct::set_ix_attrib(int value) { control_struct_set_ix_attrib(fortran_ptr_, value); }
std::string ControlVar1Struct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, control_var1_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void ControlVar1Struct::set_name(const std::string &value) {
  control_var1_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double ControlVar1Struct::value() const {
  double value;
  control_var1_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ControlVar1Struct::set_value(double value) {
  control_var1_struct_set_value(fortran_ptr_, value);
}
double ControlVar1Struct::old_value() const {
  double value;
  control_var1_struct_get_old_value(fortran_ptr_, &value);
  return value;
}
void ControlVar1Struct::set_old_value(double value) {
  control_var1_struct_set_old_value(fortran_ptr_, value);
}
RealAlloc1D ControlRamp1Struct::y_knot() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      control_ramp1_struct_reallocate_y_knot,
      control_ramp1_struct_get_y_knot_info
  );
}
void ControlRamp1Struct::set_y_knot(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  control_ramp1_struct_set_y_knot(fortran_ptr_, v.data(), shape);
}
ExpressionAtomStructAlloc1D ControlRamp1Struct::stack() const {
  return ExpressionAtomStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      control_ramp1_struct_reallocate_stack,
      control_ramp1_struct_get_stack_info
  );
}
std::string ControlRamp1Struct::attribute() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, control_ramp1_struct_get_attribute_info);
  return std::string(arr.data(), arr.size());
}
void ControlRamp1Struct::set_attribute(const std::string &value) {
  control_ramp1_struct_set_attribute(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string ControlRamp1Struct::slave_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, control_ramp1_struct_get_slave_name_info);
  return std::string(arr.data(), arr.size());
}
void ControlRamp1Struct::set_slave_name(const std::string &value) {
  control_ramp1_struct_set_slave_name(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
bool ControlRamp1Struct::is_controller() const {
  bool value;
  control_ramp1_struct_get_is_controller(fortran_ptr_, &value);
  return value;
}
void ControlRamp1Struct::set_is_controller(bool value) {
  control_ramp1_struct_set_is_controller(fortran_ptr_, value);
}
ControlVar1StructAlloc1D ControllerStruct::var() const {
  return ControlVar1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      controller_struct_reallocate_var,
      controller_struct_get_var_info
  );
}
ControlRamp1StructAlloc1D ControllerStruct::ramp() const {
  return ControlRamp1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      controller_struct_reallocate_ramp,
      controller_struct_get_ramp_info
  );
}
RamperLordStructAlloc1D ControllerStruct::ramper_lord() const {
  return RamperLordStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      controller_struct_reallocate_ramper_lord,
      controller_struct_get_ramper_lord_info
  );
}
RealAlloc1D ControllerStruct::x_knot() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      controller_struct_reallocate_x_knot,
      controller_struct_get_x_knot_info
  );
}
void ControllerStruct::set_x_knot(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  controller_struct_set_x_knot(fortran_ptr_, v.data(), shape);
}
int EllipseBeamInitStruct::part_per_ellipse() const {
  int value;
  ellipse_beam_init_struct_get_part_per_ellipse(fortran_ptr_, &value);
  return value;
}
void EllipseBeamInitStruct::set_part_per_ellipse(int value) {
  ellipse_beam_init_struct_set_part_per_ellipse(fortran_ptr_, value);
}
int EllipseBeamInitStruct::n_ellipse() const {
  int value;
  ellipse_beam_init_struct_get_n_ellipse(fortran_ptr_, &value);
  return value;
}
void EllipseBeamInitStruct::set_n_ellipse(int value) {
  ellipse_beam_init_struct_set_n_ellipse(fortran_ptr_, value);
}
double EllipseBeamInitStruct::sigma_cutoff() const {
  double value;
  ellipse_beam_init_struct_get_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void EllipseBeamInitStruct::set_sigma_cutoff(double value) {
  ellipse_beam_init_struct_set_sigma_cutoff(fortran_ptr_, value);
}
FArray1D<int> KvBeamInitStruct::part_per_phi() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, kv_beam_init_struct_get_part_per_phi_info);
}
void KvBeamInitStruct::set_part_per_phi(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  kv_beam_init_struct_set_part_per_phi(fortran_ptr_, v.data(), shape);
}
int KvBeamInitStruct::n_I2() const {
  int value;
  kv_beam_init_struct_get_n_I2(fortran_ptr_, &value);
  return value;
}
void KvBeamInitStruct::set_n_I2(int value) { kv_beam_init_struct_set_n_I2(fortran_ptr_, value); }
double KvBeamInitStruct::A() const {
  double value;
  kv_beam_init_struct_get_A(fortran_ptr_, &value);
  return value;
}
void KvBeamInitStruct::set_A(double value) { kv_beam_init_struct_set_A(fortran_ptr_, value); }
int GridBeamInitStruct::n_x() const {
  int value;
  grid_beam_init_struct_get_n_x(fortran_ptr_, &value);
  return value;
}
void GridBeamInitStruct::set_n_x(int value) { grid_beam_init_struct_set_n_x(fortran_ptr_, value); }
int GridBeamInitStruct::n_px() const {
  int value;
  grid_beam_init_struct_get_n_px(fortran_ptr_, &value);
  return value;
}
void GridBeamInitStruct::set_n_px(int value) {
  grid_beam_init_struct_set_n_px(fortran_ptr_, value);
}
double GridBeamInitStruct::x_min() const {
  double value;
  grid_beam_init_struct_get_x_min(fortran_ptr_, &value);
  return value;
}
void GridBeamInitStruct::set_x_min(double value) {
  grid_beam_init_struct_set_x_min(fortran_ptr_, value);
}
double GridBeamInitStruct::x_max() const {
  double value;
  grid_beam_init_struct_get_x_max(fortran_ptr_, &value);
  return value;
}
void GridBeamInitStruct::set_x_max(double value) {
  grid_beam_init_struct_set_x_max(fortran_ptr_, value);
}
double GridBeamInitStruct::px_min() const {
  double value;
  grid_beam_init_struct_get_px_min(fortran_ptr_, &value);
  return value;
}
void GridBeamInitStruct::set_px_min(double value) {
  grid_beam_init_struct_set_px_min(fortran_ptr_, value);
}
double GridBeamInitStruct::px_max() const {
  double value;
  grid_beam_init_struct_get_px_max(fortran_ptr_, &value);
  return value;
}
void GridBeamInitStruct::set_px_max(double value) {
  grid_beam_init_struct_set_px_max(fortran_ptr_, value);
}
std::string BeamInitStruct::position_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, beam_init_struct_get_position_file_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitStruct::set_position_file(const std::string &value) {
  beam_init_struct_set_position_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D BeamInitStruct::distribution_type() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, beam_init_struct_get_distribution_type_info);
}
FArray1D<double> BeamInitStruct::spin() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, beam_init_struct_get_spin_info);
}
void BeamInitStruct::set_spin(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  beam_init_struct_set_spin(fortran_ptr_, v.data(), shape);
}
EllipseBeamInitStructArray1D BeamInitStruct::ellipse() const {
  return ProxyHelpers::get_type_array_1d<EllipseBeamInitStructArray1D>(
      fortran_ptr_,
      beam_init_struct_get_ellipse_info
  );
}
KvBeamInitStruct BeamInitStruct::KV() const {
  void *ptr;
  beam_init_struct_get_KV(fortran_ptr_, &ptr);
  return KvBeamInitStruct(ptr);
}
void BeamInitStruct::set_KV(const KvBeamInitStruct &src) {
  beam_init_struct_set_KV(fortran_ptr_, src.get_fortran_ptr());
}
GridBeamInitStructArray1D BeamInitStruct::grid() const {
  return ProxyHelpers::get_type_array_1d<GridBeamInitStructArray1D>(
      fortran_ptr_,
      beam_init_struct_get_grid_info
  );
}
FArray1D<double> BeamInitStruct::center_jitter() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, beam_init_struct_get_center_jitter_info);
}
void BeamInitStruct::set_center_jitter(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  beam_init_struct_set_center_jitter(fortran_ptr_, v.data(), shape);
}
FArray1D<double> BeamInitStruct::emit_jitter() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, beam_init_struct_get_emit_jitter_info);
}
void BeamInitStruct::set_emit_jitter(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  beam_init_struct_set_emit_jitter(fortran_ptr_, v.data(), shape);
}
double BeamInitStruct::sig_z_jitter() const {
  double value;
  beam_init_struct_get_sig_z_jitter(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_sig_z_jitter(double value) {
  beam_init_struct_set_sig_z_jitter(fortran_ptr_, value);
}
double BeamInitStruct::sig_pz_jitter() const {
  double value;
  beam_init_struct_get_sig_pz_jitter(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_sig_pz_jitter(double value) {
  beam_init_struct_set_sig_pz_jitter(fortran_ptr_, value);
}
int BeamInitStruct::n_particle() const {
  int value;
  beam_init_struct_get_n_particle(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_n_particle(int value) {
  beam_init_struct_set_n_particle(fortran_ptr_, value);
}
bool BeamInitStruct::renorm_center() const {
  bool value;
  beam_init_struct_get_renorm_center(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_renorm_center(bool value) {
  beam_init_struct_set_renorm_center(fortran_ptr_, value);
}
bool BeamInitStruct::renorm_sigma() const {
  bool value;
  beam_init_struct_get_renorm_sigma(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_renorm_sigma(bool value) {
  beam_init_struct_set_renorm_sigma(fortran_ptr_, value);
}
std::string BeamInitStruct::random_engine() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, beam_init_struct_get_random_engine_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitStruct::set_random_engine(const std::string &value) {
  beam_init_struct_set_random_engine(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string BeamInitStruct::random_gauss_converter() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      beam_init_struct_get_random_gauss_converter_info
  );
  return std::string(arr.data(), arr.size());
}
void BeamInitStruct::set_random_gauss_converter(const std::string &value) {
  beam_init_struct_set_random_gauss_converter(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
double BeamInitStruct::random_sigma_cutoff() const {
  double value;
  beam_init_struct_get_random_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_random_sigma_cutoff(double value) {
  beam_init_struct_set_random_sigma_cutoff(fortran_ptr_, value);
}
double BeamInitStruct::a_norm_emit() const {
  double value;
  beam_init_struct_get_a_norm_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_a_norm_emit(double value) {
  beam_init_struct_set_a_norm_emit(fortran_ptr_, value);
}
double BeamInitStruct::b_norm_emit() const {
  double value;
  beam_init_struct_get_b_norm_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_b_norm_emit(double value) {
  beam_init_struct_set_b_norm_emit(fortran_ptr_, value);
}
double BeamInitStruct::a_emit() const {
  double value;
  beam_init_struct_get_a_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_a_emit(double value) { beam_init_struct_set_a_emit(fortran_ptr_, value); }
double BeamInitStruct::b_emit() const {
  double value;
  beam_init_struct_get_b_emit(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_b_emit(double value) { beam_init_struct_set_b_emit(fortran_ptr_, value); }
double BeamInitStruct::dPz_dz() const {
  double value;
  beam_init_struct_get_dPz_dz(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_dPz_dz(double value) { beam_init_struct_set_dPz_dz(fortran_ptr_, value); }
FArray1D<double> BeamInitStruct::center() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, beam_init_struct_get_center_info);
}
void BeamInitStruct::set_center(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  beam_init_struct_set_center(fortran_ptr_, v.data(), shape);
}
double BeamInitStruct::t_offset() const {
  double value;
  beam_init_struct_get_t_offset(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_t_offset(double value) {
  beam_init_struct_set_t_offset(fortran_ptr_, value);
}
double BeamInitStruct::dt_bunch() const {
  double value;
  beam_init_struct_get_dt_bunch(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_dt_bunch(double value) {
  beam_init_struct_set_dt_bunch(fortran_ptr_, value);
}
double BeamInitStruct::sig_z() const {
  double value;
  beam_init_struct_get_sig_z(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_sig_z(double value) { beam_init_struct_set_sig_z(fortran_ptr_, value); }
double BeamInitStruct::sig_pz() const {
  double value;
  beam_init_struct_get_sig_pz(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_sig_pz(double value) { beam_init_struct_set_sig_pz(fortran_ptr_, value); }
double BeamInitStruct::bunch_charge() const {
  double value;
  beam_init_struct_get_bunch_charge(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_bunch_charge(double value) {
  beam_init_struct_set_bunch_charge(fortran_ptr_, value);
}
int BeamInitStruct::n_bunch() const {
  int value;
  beam_init_struct_get_n_bunch(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_n_bunch(int value) { beam_init_struct_set_n_bunch(fortran_ptr_, value); }
int BeamInitStruct::ix_turn() const {
  int value;
  beam_init_struct_get_ix_turn(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_ix_turn(int value) { beam_init_struct_set_ix_turn(fortran_ptr_, value); }
std::string BeamInitStruct::species() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, beam_init_struct_get_species_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitStruct::set_species(const std::string &value) {
  beam_init_struct_set_species(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool BeamInitStruct::full_6D_coupling_calc() const {
  bool value;
  beam_init_struct_get_full_6D_coupling_calc(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_full_6D_coupling_calc(bool value) {
  beam_init_struct_set_full_6D_coupling_calc(fortran_ptr_, value);
}
bool BeamInitStruct::use_particle_start() const {
  bool value;
  beam_init_struct_get_use_particle_start(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_use_particle_start(bool value) {
  beam_init_struct_set_use_particle_start(fortran_ptr_, value);
}
bool BeamInitStruct::use_t_coords() const {
  bool value;
  beam_init_struct_get_use_t_coords(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_use_t_coords(bool value) {
  beam_init_struct_set_use_t_coords(fortran_ptr_, value);
}
bool BeamInitStruct::use_z_as_t() const {
  bool value;
  beam_init_struct_get_use_z_as_t(fortran_ptr_, &value);
  return value;
}
void BeamInitStruct::set_use_z_as_t(bool value) {
  beam_init_struct_set_use_z_as_t(fortran_ptr_, value);
}
std::string BeamInitStruct::file_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, beam_init_struct_get_file_name_info);
  return std::string(arr.data(), arr.size());
}
void BeamInitStruct::set_file_name(const std::string &value) {
  beam_init_struct_set_file_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double LatParamStruct::n_part() const {
  double value;
  lat_param_struct_get_n_part(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_n_part(double value) { lat_param_struct_set_n_part(fortran_ptr_, value); }
double LatParamStruct::total_length() const {
  double value;
  lat_param_struct_get_total_length(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_total_length(double value) {
  lat_param_struct_set_total_length(fortran_ptr_, value);
}
double LatParamStruct::unstable_factor() const {
  double value;
  lat_param_struct_get_unstable_factor(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_unstable_factor(double value) {
  lat_param_struct_set_unstable_factor(fortran_ptr_, value);
}
FArray2D<double> LatParamStruct::t1_with_RF() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, lat_param_struct_get_t1_with_RF_info);
}
void LatParamStruct::set_t1_with_RF(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, lat_param_struct_set_t1_with_RF, v);
}
FArray2D<double> LatParamStruct::t1_no_RF() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, lat_param_struct_get_t1_no_RF_info);
}
void LatParamStruct::set_t1_no_RF(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, lat_param_struct_set_t1_no_RF, v);
}
double LatParamStruct::spin_tune() const {
  double value;
  lat_param_struct_get_spin_tune(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_spin_tune(double value) {
  lat_param_struct_set_spin_tune(fortran_ptr_, value);
}
int LatParamStruct::particle() const {
  int value;
  lat_param_struct_get_particle(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_particle(int value) { lat_param_struct_set_particle(fortran_ptr_, value); }
int LatParamStruct::default_tracking_species() const {
  int value;
  lat_param_struct_get_default_tracking_species(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_default_tracking_species(int value) {
  lat_param_struct_set_default_tracking_species(fortran_ptr_, value);
}
int LatParamStruct::geometry() const {
  int value;
  lat_param_struct_get_geometry(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_geometry(int value) { lat_param_struct_set_geometry(fortran_ptr_, value); }
int LatParamStruct::ixx() const {
  int value;
  lat_param_struct_get_ixx(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_ixx(int value) { lat_param_struct_set_ixx(fortran_ptr_, value); }
bool LatParamStruct::stable() const {
  bool value;
  lat_param_struct_get_stable(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_stable(bool value) { lat_param_struct_set_stable(fortran_ptr_, value); }
bool LatParamStruct::live_branch() const {
  bool value;
  lat_param_struct_get_live_branch(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_live_branch(bool value) {
  lat_param_struct_set_live_branch(fortran_ptr_, value);
}
double LatParamStruct::g1_integral() const {
  double value;
  lat_param_struct_get_g1_integral(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_g1_integral(double value) {
  lat_param_struct_set_g1_integral(fortran_ptr_, value);
}
double LatParamStruct::g2_integral() const {
  double value;
  lat_param_struct_get_g2_integral(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_g2_integral(double value) {
  lat_param_struct_set_g2_integral(fortran_ptr_, value);
}
double LatParamStruct::g3_integral() const {
  double value;
  lat_param_struct_get_g3_integral(fortran_ptr_, &value);
  return value;
}
void LatParamStruct::set_g3_integral(double value) {
  lat_param_struct_set_g3_integral(fortran_ptr_, value);
}
BookkeepingStateStruct LatParamStruct::bookkeeping_state() const {
  void *ptr;
  lat_param_struct_get_bookkeeping_state(fortran_ptr_, &ptr);
  return BookkeepingStateStruct(ptr);
}
void LatParamStruct::set_bookkeeping_state(const BookkeepingStateStruct &src) {
  lat_param_struct_set_bookkeeping_state(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitStruct LatParamStruct::beam_init() const {
  void *ptr;
  lat_param_struct_get_beam_init(fortran_ptr_, &ptr);
  return BeamInitStruct(ptr);
}
void LatParamStruct::set_beam_init(const BeamInitStruct &src) {
  lat_param_struct_set_beam_init(fortran_ptr_, src.get_fortran_ptr());
}
bool ModeInfoStruct::stable() const {
  bool value;
  mode_info_struct_get_stable(fortran_ptr_, &value);
  return value;
}
void ModeInfoStruct::set_stable(bool value) { mode_info_struct_set_stable(fortran_ptr_, value); }
double ModeInfoStruct::tune() const {
  double value;
  mode_info_struct_get_tune(fortran_ptr_, &value);
  return value;
}
void ModeInfoStruct::set_tune(double value) { mode_info_struct_set_tune(fortran_ptr_, value); }
double ModeInfoStruct::emit() const {
  double value;
  mode_info_struct_get_emit(fortran_ptr_, &value);
  return value;
}
void ModeInfoStruct::set_emit(double value) { mode_info_struct_set_emit(fortran_ptr_, value); }
double ModeInfoStruct::chrom() const {
  double value;
  mode_info_struct_get_chrom(fortran_ptr_, &value);
  return value;
}
void ModeInfoStruct::set_chrom(double value) { mode_info_struct_set_chrom(fortran_ptr_, value); }
double ModeInfoStruct::sigma() const {
  double value;
  mode_info_struct_get_sigma(fortran_ptr_, &value);
  return value;
}
void ModeInfoStruct::set_sigma(double value) { mode_info_struct_set_sigma(fortran_ptr_, value); }
double ModeInfoStruct::sigmap() const {
  double value;
  mode_info_struct_get_sigmap(fortran_ptr_, &value);
  return value;
}
void ModeInfoStruct::set_sigmap(double value) { mode_info_struct_set_sigmap(fortran_ptr_, value); }
int PreTrackerStruct::who() const {
  int value;
  pre_tracker_struct_get_who(fortran_ptr_, &value);
  return value;
}
void PreTrackerStruct::set_who(int value) { pre_tracker_struct_set_who(fortran_ptr_, value); }
int PreTrackerStruct::ix_ele_start() const {
  int value;
  pre_tracker_struct_get_ix_ele_start(fortran_ptr_, &value);
  return value;
}
void PreTrackerStruct::set_ix_ele_start(int value) {
  pre_tracker_struct_set_ix_ele_start(fortran_ptr_, value);
}
int PreTrackerStruct::ix_ele_end() const {
  int value;
  pre_tracker_struct_get_ix_ele_end(fortran_ptr_, &value);
  return value;
}
void PreTrackerStruct::set_ix_ele_end(int value) {
  pre_tracker_struct_set_ix_ele_end(fortran_ptr_, value);
}
std::string PreTrackerStruct::input_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, pre_tracker_struct_get_input_file_info);
  return std::string(arr.data(), arr.size());
}
void PreTrackerStruct::set_input_file(const std::string &value) {
  pre_tracker_struct_set_input_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double AnormalModeStruct::emittance() const {
  double value;
  anormal_mode_struct_get_emittance(fortran_ptr_, &value);
  return value;
}
void AnormalModeStruct::set_emittance(double value) {
  anormal_mode_struct_set_emittance(fortran_ptr_, value);
}
double AnormalModeStruct::emittance_no_vert() const {
  double value;
  anormal_mode_struct_get_emittance_no_vert(fortran_ptr_, &value);
  return value;
}
void AnormalModeStruct::set_emittance_no_vert(double value) {
  anormal_mode_struct_set_emittance_no_vert(fortran_ptr_, value);
}
FArray1D<double> AnormalModeStruct::synch_int() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, anormal_mode_struct_get_synch_int_info);
}
void AnormalModeStruct::set_synch_int(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  anormal_mode_struct_set_synch_int(fortran_ptr_, v.data(), shape);
}
double AnormalModeStruct::j_damp() const {
  double value;
  anormal_mode_struct_get_j_damp(fortran_ptr_, &value);
  return value;
}
void AnormalModeStruct::set_j_damp(double value) {
  anormal_mode_struct_set_j_damp(fortran_ptr_, value);
}
double AnormalModeStruct::alpha_damp() const {
  double value;
  anormal_mode_struct_get_alpha_damp(fortran_ptr_, &value);
  return value;
}
void AnormalModeStruct::set_alpha_damp(double value) {
  anormal_mode_struct_set_alpha_damp(fortran_ptr_, value);
}
double AnormalModeStruct::chrom() const {
  double value;
  anormal_mode_struct_get_chrom(fortran_ptr_, &value);
  return value;
}
void AnormalModeStruct::set_chrom(double value) {
  anormal_mode_struct_set_chrom(fortran_ptr_, value);
}
double AnormalModeStruct::tune() const {
  double value;
  anormal_mode_struct_get_tune(fortran_ptr_, &value);
  return value;
}
void AnormalModeStruct::set_tune(double value) {
  anormal_mode_struct_set_tune(fortran_ptr_, value);
}
double LinacNormalModeStruct::i2_E4() const {
  double value;
  linac_normal_mode_struct_get_i2_E4(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_i2_E4(double value) {
  linac_normal_mode_struct_set_i2_E4(fortran_ptr_, value);
}
double LinacNormalModeStruct::i3_E7() const {
  double value;
  linac_normal_mode_struct_get_i3_E7(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_i3_E7(double value) {
  linac_normal_mode_struct_set_i3_E7(fortran_ptr_, value);
}
double LinacNormalModeStruct::i5a_E6() const {
  double value;
  linac_normal_mode_struct_get_i5a_E6(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_i5a_E6(double value) {
  linac_normal_mode_struct_set_i5a_E6(fortran_ptr_, value);
}
double LinacNormalModeStruct::i5b_E6() const {
  double value;
  linac_normal_mode_struct_get_i5b_E6(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_i5b_E6(double value) {
  linac_normal_mode_struct_set_i5b_E6(fortran_ptr_, value);
}
double LinacNormalModeStruct::sig_E1() const {
  double value;
  linac_normal_mode_struct_get_sig_E1(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_sig_E1(double value) {
  linac_normal_mode_struct_set_sig_E1(fortran_ptr_, value);
}
double LinacNormalModeStruct::a_emittance_end() const {
  double value;
  linac_normal_mode_struct_get_a_emittance_end(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_a_emittance_end(double value) {
  linac_normal_mode_struct_set_a_emittance_end(fortran_ptr_, value);
}
double LinacNormalModeStruct::b_emittance_end() const {
  double value;
  linac_normal_mode_struct_get_b_emittance_end(fortran_ptr_, &value);
  return value;
}
void LinacNormalModeStruct::set_b_emittance_end(double value) {
  linac_normal_mode_struct_set_b_emittance_end(fortran_ptr_, value);
}
FArray1D<double> NormalModesStruct::synch_int() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, normal_modes_struct_get_synch_int_info);
}
void NormalModesStruct::set_synch_int(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  normal_modes_struct_set_synch_int(fortran_ptr_, v.data(), shape);
}
double NormalModesStruct::sigE_E() const {
  double value;
  normal_modes_struct_get_sigE_E(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_sigE_E(double value) {
  normal_modes_struct_set_sigE_E(fortran_ptr_, value);
}
double NormalModesStruct::sig_z() const {
  double value;
  normal_modes_struct_get_sig_z(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_sig_z(double value) {
  normal_modes_struct_set_sig_z(fortran_ptr_, value);
}
double NormalModesStruct::e_loss() const {
  double value;
  normal_modes_struct_get_e_loss(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_e_loss(double value) {
  normal_modes_struct_set_e_loss(fortran_ptr_, value);
}
double NormalModesStruct::rf_voltage() const {
  double value;
  normal_modes_struct_get_rf_voltage(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_rf_voltage(double value) {
  normal_modes_struct_set_rf_voltage(fortran_ptr_, value);
}
double NormalModesStruct::pz_aperture() const {
  double value;
  normal_modes_struct_get_pz_aperture(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_pz_aperture(double value) {
  normal_modes_struct_set_pz_aperture(fortran_ptr_, value);
}
double NormalModesStruct::pz_average() const {
  double value;
  normal_modes_struct_get_pz_average(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_pz_average(double value) {
  normal_modes_struct_set_pz_average(fortran_ptr_, value);
}
double NormalModesStruct::momentum_compaction() const {
  double value;
  normal_modes_struct_get_momentum_compaction(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_momentum_compaction(double value) {
  normal_modes_struct_set_momentum_compaction(fortran_ptr_, value);
}
double NormalModesStruct::dpz_damp() const {
  double value;
  normal_modes_struct_get_dpz_damp(fortran_ptr_, &value);
  return value;
}
void NormalModesStruct::set_dpz_damp(double value) {
  normal_modes_struct_set_dpz_damp(fortran_ptr_, value);
}
AnormalModeStruct NormalModesStruct::a() const {
  void *ptr;
  normal_modes_struct_get_a(fortran_ptr_, &ptr);
  return AnormalModeStruct(ptr);
}
void NormalModesStruct::set_a(const AnormalModeStruct &src) {
  normal_modes_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
AnormalModeStruct NormalModesStruct::b() const {
  void *ptr;
  normal_modes_struct_get_b(fortran_ptr_, &ptr);
  return AnormalModeStruct(ptr);
}
void NormalModesStruct::set_b(const AnormalModeStruct &src) {
  normal_modes_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
AnormalModeStruct NormalModesStruct::z() const {
  void *ptr;
  normal_modes_struct_get_z(fortran_ptr_, &ptr);
  return AnormalModeStruct(ptr);
}
void NormalModesStruct::set_z(const AnormalModeStruct &src) {
  normal_modes_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
LinacNormalModeStruct NormalModesStruct::lin() const {
  void *ptr;
  normal_modes_struct_get_lin(fortran_ptr_, &ptr);
  return LinacNormalModeStruct(ptr);
}
void NormalModesStruct::set_lin(const LinacNormalModeStruct &src) {
  normal_modes_struct_set_lin(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> EmFieldStruct::E() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, em_field_struct_get_E_info);
}
void EmFieldStruct::set_E(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  em_field_struct_set_E(fortran_ptr_, v.data(), shape);
}
FArray1D<double> EmFieldStruct::B() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, em_field_struct_get_B_info);
}
void EmFieldStruct::set_B(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  em_field_struct_set_B(fortran_ptr_, v.data(), shape);
}
FArray2D<double> EmFieldStruct::dE() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, em_field_struct_get_dE_info);
}
void EmFieldStruct::set_dE(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, em_field_struct_set_dE, v);
}
FArray2D<double> EmFieldStruct::dB() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, em_field_struct_get_dB_info);
}
void EmFieldStruct::set_dB(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, em_field_struct_set_dB, v);
}
double EmFieldStruct::phi() const {
  double value;
  em_field_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void EmFieldStruct::set_phi(double value) { em_field_struct_set_phi(fortran_ptr_, value); }
double EmFieldStruct::phi_B() const {
  double value;
  em_field_struct_get_phi_B(fortran_ptr_, &value);
  return value;
}
void EmFieldStruct::set_phi_B(double value) { em_field_struct_set_phi_B(fortran_ptr_, value); }
FArray1D<double> EmFieldStruct::A() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, em_field_struct_get_A_info);
}
void EmFieldStruct::set_A(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  em_field_struct_set_A(fortran_ptr_, v.data(), shape);
}
int StrongBeamStruct::ix_slice() const {
  int value;
  strong_beam_struct_get_ix_slice(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_ix_slice(int value) {
  strong_beam_struct_set_ix_slice(fortran_ptr_, value);
}
double StrongBeamStruct::x_center() const {
  double value;
  strong_beam_struct_get_x_center(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_x_center(double value) {
  strong_beam_struct_set_x_center(fortran_ptr_, value);
}
double StrongBeamStruct::y_center() const {
  double value;
  strong_beam_struct_get_y_center(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_y_center(double value) {
  strong_beam_struct_set_y_center(fortran_ptr_, value);
}
double StrongBeamStruct::x_sigma() const {
  double value;
  strong_beam_struct_get_x_sigma(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_x_sigma(double value) {
  strong_beam_struct_set_x_sigma(fortran_ptr_, value);
}
double StrongBeamStruct::y_sigma() const {
  double value;
  strong_beam_struct_get_y_sigma(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_y_sigma(double value) {
  strong_beam_struct_set_y_sigma(fortran_ptr_, value);
}
double StrongBeamStruct::dx() const {
  double value;
  strong_beam_struct_get_dx(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_dx(double value) { strong_beam_struct_set_dx(fortran_ptr_, value); }
double StrongBeamStruct::dy() const {
  double value;
  strong_beam_struct_get_dy(fortran_ptr_, &value);
  return value;
}
void StrongBeamStruct::set_dy(double value) { strong_beam_struct_set_dy(fortran_ptr_, value); }
double TrackPointStruct::s_lab() const {
  double value;
  track_point_struct_get_s_lab(fortran_ptr_, &value);
  return value;
}
void TrackPointStruct::set_s_lab(double value) {
  track_point_struct_set_s_lab(fortran_ptr_, value);
}
double TrackPointStruct::s_body() const {
  double value;
  track_point_struct_get_s_body(fortran_ptr_, &value);
  return value;
}
void TrackPointStruct::set_s_body(double value) {
  track_point_struct_set_s_body(fortran_ptr_, value);
}
CoordStruct TrackPointStruct::orb() const {
  void *ptr;
  track_point_struct_get_orb(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void TrackPointStruct::set_orb(const CoordStruct &src) {
  track_point_struct_set_orb(fortran_ptr_, src.get_fortran_ptr());
}
EmFieldStruct TrackPointStruct::field() const {
  void *ptr;
  track_point_struct_get_field(fortran_ptr_, &ptr);
  return EmFieldStruct(ptr);
}
void TrackPointStruct::set_field(const EmFieldStruct &src) {
  track_point_struct_set_field(fortran_ptr_, src.get_fortran_ptr());
}
StrongBeamStruct TrackPointStruct::strong_beam() const {
  void *ptr;
  track_point_struct_get_strong_beam(fortran_ptr_, &ptr);
  return StrongBeamStruct(ptr);
}
void TrackPointStruct::set_strong_beam(const StrongBeamStruct &src) {
  track_point_struct_set_strong_beam(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TrackPointStruct::vec0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, track_point_struct_get_vec0_info);
}
void TrackPointStruct::set_vec0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  track_point_struct_set_vec0(fortran_ptr_, v.data(), shape);
}
FArray2D<double> TrackPointStruct::mat6() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, track_point_struct_get_mat6_info);
}
void TrackPointStruct::set_mat6(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, track_point_struct_set_mat6, v);
}
TrackPointStructAlloc1D TrackStruct::pt() const {
  return TrackPointStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      track_struct_reallocate_pt,
      track_struct_get_pt_info
  );
}
double TrackStruct::ds_save() const {
  double value;
  track_struct_get_ds_save(fortran_ptr_, &value);
  return value;
}
void TrackStruct::set_ds_save(double value) { track_struct_set_ds_save(fortran_ptr_, value); }
int TrackStruct::n_pt() const {
  int value;
  track_struct_get_n_pt(fortran_ptr_, &value);
  return value;
}
void TrackStruct::set_n_pt(int value) { track_struct_set_n_pt(fortran_ptr_, value); }
int TrackStruct::n_bad() const {
  int value;
  track_struct_get_n_bad(fortran_ptr_, &value);
  return value;
}
void TrackStruct::set_n_bad(int value) { track_struct_set_n_bad(fortran_ptr_, value); }
int TrackStruct::n_ok() const {
  int value;
  track_struct_get_n_ok(fortran_ptr_, &value);
  return value;
}
void TrackStruct::set_n_ok(int value) { track_struct_set_n_ok(fortran_ptr_, value); }
double SpaceChargeCommonStruct::ds_track_step() const {
  double value;
  space_charge_common_struct_get_ds_track_step(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_ds_track_step(double value) {
  space_charge_common_struct_set_ds_track_step(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::dt_track_step() const {
  double value;
  space_charge_common_struct_get_dt_track_step(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_dt_track_step(double value) {
  space_charge_common_struct_set_dt_track_step(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::cathode_strength_cutoff() const {
  double value;
  space_charge_common_struct_get_cathode_strength_cutoff(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_cathode_strength_cutoff(double value) {
  space_charge_common_struct_set_cathode_strength_cutoff(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::rel_tol_tracking() const {
  double value;
  space_charge_common_struct_get_rel_tol_tracking(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_rel_tol_tracking(double value) {
  space_charge_common_struct_set_rel_tol_tracking(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::abs_tol_tracking() const {
  double value;
  space_charge_common_struct_get_abs_tol_tracking(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_abs_tol_tracking(double value) {
  space_charge_common_struct_set_abs_tol_tracking(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::beam_chamber_height() const {
  double value;
  space_charge_common_struct_get_beam_chamber_height(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_beam_chamber_height(double value) {
  space_charge_common_struct_set_beam_chamber_height(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::lsc_sigma_cutoff() const {
  double value;
  space_charge_common_struct_get_lsc_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_lsc_sigma_cutoff(double value) {
  space_charge_common_struct_set_lsc_sigma_cutoff(fortran_ptr_, value);
}
double SpaceChargeCommonStruct::particle_sigma_cutoff() const {
  double value;
  space_charge_common_struct_get_particle_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_particle_sigma_cutoff(double value) {
  space_charge_common_struct_set_particle_sigma_cutoff(fortran_ptr_, value);
}
FArray1D<int> SpaceChargeCommonStruct::space_charge_mesh_size() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_,
      space_charge_common_struct_get_space_charge_mesh_size_info
  );
}
void SpaceChargeCommonStruct::set_space_charge_mesh_size(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  space_charge_common_struct_set_space_charge_mesh_size(fortran_ptr_, v.data(), shape);
}
FArray1D<int> SpaceChargeCommonStruct::csr3d_mesh_size() const {
  return ProxyHelpers::get_array_1d<int>(
      fortran_ptr_,
      space_charge_common_struct_get_csr3d_mesh_size_info
  );
}
void SpaceChargeCommonStruct::set_csr3d_mesh_size(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  space_charge_common_struct_set_csr3d_mesh_size(fortran_ptr_, v.data(), shape);
}
int SpaceChargeCommonStruct::n_bin() const {
  int value;
  space_charge_common_struct_get_n_bin(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_n_bin(int value) {
  space_charge_common_struct_set_n_bin(fortran_ptr_, value);
}
int SpaceChargeCommonStruct::particle_bin_span() const {
  int value;
  space_charge_common_struct_get_particle_bin_span(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_particle_bin_span(int value) {
  space_charge_common_struct_set_particle_bin_span(fortran_ptr_, value);
}
int SpaceChargeCommonStruct::n_shield_images() const {
  int value;
  space_charge_common_struct_get_n_shield_images(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_n_shield_images(int value) {
  space_charge_common_struct_set_n_shield_images(fortran_ptr_, value);
}
int SpaceChargeCommonStruct::sc_min_in_bin() const {
  int value;
  space_charge_common_struct_get_sc_min_in_bin(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_sc_min_in_bin(int value) {
  space_charge_common_struct_set_sc_min_in_bin(fortran_ptr_, value);
}
bool SpaceChargeCommonStruct::lsc_kick_transverse_dependence() const {
  bool value;
  space_charge_common_struct_get_lsc_kick_transverse_dependence(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_lsc_kick_transverse_dependence(bool value) {
  space_charge_common_struct_set_lsc_kick_transverse_dependence(fortran_ptr_, value);
}
bool SpaceChargeCommonStruct::debug() const {
  bool value;
  space_charge_common_struct_get_debug(fortran_ptr_, &value);
  return value;
}
void SpaceChargeCommonStruct::set_debug(bool value) {
  space_charge_common_struct_set_debug(fortran_ptr_, value);
}
std::string SpaceChargeCommonStruct::diagnostic_output_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      space_charge_common_struct_get_diagnostic_output_file_info
  );
  return std::string(arr.data(), arr.size());
}
void SpaceChargeCommonStruct::set_diagnostic_output_file(const std::string &value) {
  space_charge_common_struct_set_diagnostic_output_file(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
double BmadCommonStruct::max_aperture_limit() const {
  double value;
  bmad_common_struct_get_max_aperture_limit(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_max_aperture_limit(double value) {
  bmad_common_struct_set_max_aperture_limit(fortran_ptr_, value);
}
FArray1D<double> BmadCommonStruct::d_orb() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bmad_common_struct_get_d_orb_info);
}
void BmadCommonStruct::set_d_orb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bmad_common_struct_set_d_orb(fortran_ptr_, v.data(), shape);
}
double BmadCommonStruct::default_ds_step() const {
  double value;
  bmad_common_struct_get_default_ds_step(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_default_ds_step(double value) {
  bmad_common_struct_set_default_ds_step(fortran_ptr_, value);
}
double BmadCommonStruct::significant_length() const {
  double value;
  bmad_common_struct_get_significant_length(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_significant_length(double value) {
  bmad_common_struct_set_significant_length(fortran_ptr_, value);
}
double BmadCommonStruct::rel_tol_tracking() const {
  double value;
  bmad_common_struct_get_rel_tol_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_rel_tol_tracking(double value) {
  bmad_common_struct_set_rel_tol_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::abs_tol_tracking() const {
  double value;
  bmad_common_struct_get_abs_tol_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_abs_tol_tracking(double value) {
  bmad_common_struct_set_abs_tol_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::rel_tol_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_rel_tol_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_rel_tol_adaptive_tracking(double value) {
  bmad_common_struct_set_rel_tol_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::abs_tol_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_abs_tol_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_abs_tol_adaptive_tracking(double value) {
  bmad_common_struct_set_abs_tol_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::init_ds_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_init_ds_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_init_ds_adaptive_tracking(double value) {
  bmad_common_struct_set_init_ds_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::min_ds_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_min_ds_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_min_ds_adaptive_tracking(double value) {
  bmad_common_struct_set_min_ds_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::fatal_ds_adaptive_tracking() const {
  double value;
  bmad_common_struct_get_fatal_ds_adaptive_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_fatal_ds_adaptive_tracking(double value) {
  bmad_common_struct_set_fatal_ds_adaptive_tracking(fortran_ptr_, value);
}
double BmadCommonStruct::autoscale_amp_abs_tol() const {
  double value;
  bmad_common_struct_get_autoscale_amp_abs_tol(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_autoscale_amp_abs_tol(double value) {
  bmad_common_struct_set_autoscale_amp_abs_tol(fortran_ptr_, value);
}
double BmadCommonStruct::autoscale_amp_rel_tol() const {
  double value;
  bmad_common_struct_get_autoscale_amp_rel_tol(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_autoscale_amp_rel_tol(double value) {
  bmad_common_struct_set_autoscale_amp_rel_tol(fortran_ptr_, value);
}
double BmadCommonStruct::autoscale_phase_tol() const {
  double value;
  bmad_common_struct_get_autoscale_phase_tol(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_autoscale_phase_tol(double value) {
  bmad_common_struct_set_autoscale_phase_tol(fortran_ptr_, value);
}
double BmadCommonStruct::electric_dipole_moment() const {
  double value;
  bmad_common_struct_get_electric_dipole_moment(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_electric_dipole_moment(double value) {
  bmad_common_struct_set_electric_dipole_moment(fortran_ptr_, value);
}
double BmadCommonStruct::synch_rad_scale() const {
  double value;
  bmad_common_struct_get_synch_rad_scale(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_synch_rad_scale(double value) {
  bmad_common_struct_set_synch_rad_scale(fortran_ptr_, value);
}
double BmadCommonStruct::sad_eps_scale() const {
  double value;
  bmad_common_struct_get_sad_eps_scale(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_sad_eps_scale(double value) {
  bmad_common_struct_set_sad_eps_scale(fortran_ptr_, value);
}
double BmadCommonStruct::sad_amp_max() const {
  double value;
  bmad_common_struct_get_sad_amp_max(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_sad_amp_max(double value) {
  bmad_common_struct_set_sad_amp_max(fortran_ptr_, value);
}
int BmadCommonStruct::sad_n_div_max() const {
  int value;
  bmad_common_struct_get_sad_n_div_max(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_sad_n_div_max(int value) {
  bmad_common_struct_set_sad_n_div_max(fortran_ptr_, value);
}
int BmadCommonStruct::taylor_order() const {
  int value;
  bmad_common_struct_get_taylor_order(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_taylor_order(int value) {
  bmad_common_struct_set_taylor_order(fortran_ptr_, value);
}
int BmadCommonStruct::runge_kutta_order() const {
  int value;
  bmad_common_struct_get_runge_kutta_order(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_runge_kutta_order(int value) {
  bmad_common_struct_set_runge_kutta_order(fortran_ptr_, value);
}
int BmadCommonStruct::default_integ_order() const {
  int value;
  bmad_common_struct_get_default_integ_order(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_default_integ_order(int value) {
  bmad_common_struct_set_default_integ_order(fortran_ptr_, value);
}
int BmadCommonStruct::max_num_runge_kutta_step() const {
  int value;
  bmad_common_struct_get_max_num_runge_kutta_step(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_max_num_runge_kutta_step(int value) {
  bmad_common_struct_set_max_num_runge_kutta_step(fortran_ptr_, value);
}
bool BmadCommonStruct::rf_phase_below_transition_ref() const {
  bool value;
  bmad_common_struct_get_rf_phase_below_transition_ref(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_rf_phase_below_transition_ref(bool value) {
  bmad_common_struct_set_rf_phase_below_transition_ref(fortran_ptr_, value);
}
bool BmadCommonStruct::sr_wakes_on() const {
  bool value;
  bmad_common_struct_get_sr_wakes_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_sr_wakes_on(bool value) {
  bmad_common_struct_set_sr_wakes_on(fortran_ptr_, value);
}
bool BmadCommonStruct::lr_wakes_on() const {
  bool value;
  bmad_common_struct_get_lr_wakes_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_lr_wakes_on(bool value) {
  bmad_common_struct_set_lr_wakes_on(fortran_ptr_, value);
}
bool BmadCommonStruct::auto_bookkeeper() const {
  bool value;
  bmad_common_struct_get_auto_bookkeeper(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_auto_bookkeeper(bool value) {
  bmad_common_struct_set_auto_bookkeeper(fortran_ptr_, value);
}
bool BmadCommonStruct::high_energy_space_charge_on() const {
  bool value;
  bmad_common_struct_get_high_energy_space_charge_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_high_energy_space_charge_on(bool value) {
  bmad_common_struct_set_high_energy_space_charge_on(fortran_ptr_, value);
}
bool BmadCommonStruct::high_energy_space_charge_linear() const {
  bool value;
  bmad_common_struct_get_high_energy_space_charge_linear(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_high_energy_space_charge_linear(bool value) {
  bmad_common_struct_set_high_energy_space_charge_linear(fortran_ptr_, value);
}
bool BmadCommonStruct::csr_and_space_charge_on() const {
  bool value;
  bmad_common_struct_get_csr_and_space_charge_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_csr_and_space_charge_on(bool value) {
  bmad_common_struct_set_csr_and_space_charge_on(fortran_ptr_, value);
}
bool BmadCommonStruct::spin_tracking_on() const {
  bool value;
  bmad_common_struct_get_spin_tracking_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_spin_tracking_on(bool value) {
  bmad_common_struct_set_spin_tracking_on(fortran_ptr_, value);
}
bool BmadCommonStruct::spin_sokolov_ternov_flipping_on() const {
  bool value;
  bmad_common_struct_get_spin_sokolov_ternov_flipping_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_spin_sokolov_ternov_flipping_on(bool value) {
  bmad_common_struct_set_spin_sokolov_ternov_flipping_on(fortran_ptr_, value);
}
bool BmadCommonStruct::radiation_damping_on() const {
  bool value;
  bmad_common_struct_get_radiation_damping_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_radiation_damping_on(bool value) {
  bmad_common_struct_set_radiation_damping_on(fortran_ptr_, value);
}
bool BmadCommonStruct::radiation_zero_average() const {
  bool value;
  bmad_common_struct_get_radiation_zero_average(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_radiation_zero_average(bool value) {
  bmad_common_struct_set_radiation_zero_average(fortran_ptr_, value);
}
bool BmadCommonStruct::radiation_fluctuations_on() const {
  bool value;
  bmad_common_struct_get_radiation_fluctuations_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_radiation_fluctuations_on(bool value) {
  bmad_common_struct_set_radiation_fluctuations_on(fortran_ptr_, value);
}
bool BmadCommonStruct::conserve_taylor_maps() const {
  bool value;
  bmad_common_struct_get_conserve_taylor_maps(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_conserve_taylor_maps(bool value) {
  bmad_common_struct_set_conserve_taylor_maps(fortran_ptr_, value);
}
bool BmadCommonStruct::absolute_time_tracking() const {
  bool value;
  bmad_common_struct_get_absolute_time_tracking(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_absolute_time_tracking(bool value) {
  bmad_common_struct_set_absolute_time_tracking(fortran_ptr_, value);
}
bool BmadCommonStruct::absolute_time_ref_shift() const {
  bool value;
  bmad_common_struct_get_absolute_time_ref_shift(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_absolute_time_ref_shift(bool value) {
  bmad_common_struct_set_absolute_time_ref_shift(fortran_ptr_, value);
}
bool BmadCommonStruct::convert_to_kinetic_momentum() const {
  bool value;
  bmad_common_struct_get_convert_to_kinetic_momentum(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_convert_to_kinetic_momentum(bool value) {
  bmad_common_struct_set_convert_to_kinetic_momentum(fortran_ptr_, value);
}
bool BmadCommonStruct::normalize_twiss() const {
  bool value;
  bmad_common_struct_get_normalize_twiss(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_normalize_twiss(bool value) {
  bmad_common_struct_set_normalize_twiss(fortran_ptr_, value);
}
bool BmadCommonStruct::aperture_limit_on() const {
  bool value;
  bmad_common_struct_get_aperture_limit_on(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_aperture_limit_on(bool value) {
  bmad_common_struct_set_aperture_limit_on(fortran_ptr_, value);
}
bool BmadCommonStruct::spin_n0_direction_user_set() const {
  bool value;
  bmad_common_struct_get_spin_n0_direction_user_set(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_spin_n0_direction_user_set(bool value) {
  bmad_common_struct_set_spin_n0_direction_user_set(fortran_ptr_, value);
}
bool BmadCommonStruct::debug() const {
  bool value;
  bmad_common_struct_get_debug(fortran_ptr_, &value);
  return value;
}
void BmadCommonStruct::set_debug(bool value) { bmad_common_struct_set_debug(fortran_ptr_, value); }
double RadInt1Struct::i0() const {
  double value;
  rad_int1_struct_get_i0(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i0(double value) { rad_int1_struct_set_i0(fortran_ptr_, value); }
double RadInt1Struct::i1() const {
  double value;
  rad_int1_struct_get_i1(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i1(double value) { rad_int1_struct_set_i1(fortran_ptr_, value); }
double RadInt1Struct::i2() const {
  double value;
  rad_int1_struct_get_i2(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i2(double value) { rad_int1_struct_set_i2(fortran_ptr_, value); }
double RadInt1Struct::i3() const {
  double value;
  rad_int1_struct_get_i3(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i3(double value) { rad_int1_struct_set_i3(fortran_ptr_, value); }
double RadInt1Struct::i4a() const {
  double value;
  rad_int1_struct_get_i4a(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i4a(double value) { rad_int1_struct_set_i4a(fortran_ptr_, value); }
double RadInt1Struct::i4b() const {
  double value;
  rad_int1_struct_get_i4b(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i4b(double value) { rad_int1_struct_set_i4b(fortran_ptr_, value); }
double RadInt1Struct::i4z() const {
  double value;
  rad_int1_struct_get_i4z(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i4z(double value) { rad_int1_struct_set_i4z(fortran_ptr_, value); }
double RadInt1Struct::i5a() const {
  double value;
  rad_int1_struct_get_i5a(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i5a(double value) { rad_int1_struct_set_i5a(fortran_ptr_, value); }
double RadInt1Struct::i5b() const {
  double value;
  rad_int1_struct_get_i5b(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i5b(double value) { rad_int1_struct_set_i5b(fortran_ptr_, value); }
double RadInt1Struct::i6b() const {
  double value;
  rad_int1_struct_get_i6b(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_i6b(double value) { rad_int1_struct_set_i6b(fortran_ptr_, value); }
double RadInt1Struct::lin_i2_E4() const {
  double value;
  rad_int1_struct_get_lin_i2_E4(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_i2_E4(double value) {
  rad_int1_struct_set_lin_i2_E4(fortran_ptr_, value);
}
double RadInt1Struct::lin_i3_E7() const {
  double value;
  rad_int1_struct_get_lin_i3_E7(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_i3_E7(double value) {
  rad_int1_struct_set_lin_i3_E7(fortran_ptr_, value);
}
double RadInt1Struct::lin_i5a_E6() const {
  double value;
  rad_int1_struct_get_lin_i5a_E6(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_i5a_E6(double value) {
  rad_int1_struct_set_lin_i5a_E6(fortran_ptr_, value);
}
double RadInt1Struct::lin_i5b_E6() const {
  double value;
  rad_int1_struct_get_lin_i5b_E6(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_i5b_E6(double value) {
  rad_int1_struct_set_lin_i5b_E6(fortran_ptr_, value);
}
double RadInt1Struct::lin_norm_emit_a() const {
  double value;
  rad_int1_struct_get_lin_norm_emit_a(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_norm_emit_a(double value) {
  rad_int1_struct_set_lin_norm_emit_a(fortran_ptr_, value);
}
double RadInt1Struct::lin_norm_emit_b() const {
  double value;
  rad_int1_struct_get_lin_norm_emit_b(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_norm_emit_b(double value) {
  rad_int1_struct_set_lin_norm_emit_b(fortran_ptr_, value);
}
double RadInt1Struct::lin_sig_E() const {
  double value;
  rad_int1_struct_get_lin_sig_E(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_lin_sig_E(double value) {
  rad_int1_struct_set_lin_sig_E(fortran_ptr_, value);
}
double RadInt1Struct::n_steps() const {
  double value;
  rad_int1_struct_get_n_steps(fortran_ptr_, &value);
  return value;
}
void RadInt1Struct::set_n_steps(double value) { rad_int1_struct_set_n_steps(fortran_ptr_, value); }
RadInt1StructAlloc1D RadIntBranchStruct::ele() const {
  return RadInt1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      rad_int_branch_struct_reallocate_ele,
      rad_int_branch_struct_get_ele_info
  );
}
RadIntBranchStructAlloc1D RadIntAllEleStruct::branch() const {
  return RadIntBranchStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      rad_int_all_ele_struct_reallocate_branch,
      rad_int_all_ele_struct_get_branch_info
  );
}
double RfStairStepStruct::E_tot0() const {
  double value;
  rf_stair_step_struct_get_E_tot0(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_E_tot0(double value) {
  rf_stair_step_struct_set_E_tot0(fortran_ptr_, value);
}
double RfStairStepStruct::E_tot1() const {
  double value;
  rf_stair_step_struct_get_E_tot1(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_E_tot1(double value) {
  rf_stair_step_struct_set_E_tot1(fortran_ptr_, value);
}
double RfStairStepStruct::p0c() const {
  double value;
  rf_stair_step_struct_get_p0c(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_p0c(double value) { rf_stair_step_struct_set_p0c(fortran_ptr_, value); }
double RfStairStepStruct::p1c() const {
  double value;
  rf_stair_step_struct_get_p1c(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_p1c(double value) { rf_stair_step_struct_set_p1c(fortran_ptr_, value); }
double RfStairStepStruct::scale() const {
  double value;
  rf_stair_step_struct_get_scale(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_scale(double value) {
  rf_stair_step_struct_set_scale(fortran_ptr_, value);
}
double RfStairStepStruct::time() const {
  double value;
  rf_stair_step_struct_get_time(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_time(double value) {
  rf_stair_step_struct_set_time(fortran_ptr_, value);
}
double RfStairStepStruct::s0() const {
  double value;
  rf_stair_step_struct_get_s0(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_s0(double value) { rf_stair_step_struct_set_s0(fortran_ptr_, value); }
double RfStairStepStruct::s() const {
  double value;
  rf_stair_step_struct_get_s(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_s(double value) { rf_stair_step_struct_set_s(fortran_ptr_, value); }
int RfStairStepStruct::ix_step() const {
  int value;
  rf_stair_step_struct_get_ix_step(fortran_ptr_, &value);
  return value;
}
void RfStairStepStruct::set_ix_step(int value) {
  rf_stair_step_struct_set_ix_step(fortran_ptr_, value);
}
RfStairStepStructAlloc1D RfEleStruct::steps() const {
  return RfStairStepStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      rf_ele_struct_reallocate_steps,
      rf_ele_struct_get_steps_info
  );
}
double RfEleStruct::ds_step() const {
  double value;
  rf_ele_struct_get_ds_step(fortran_ptr_, &value);
  return value;
}
void RfEleStruct::set_ds_step(double value) { rf_ele_struct_set_ds_step(fortran_ptr_, value); }
std::string EleStruct::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void EleStruct::set_name(const std::string &value) {
  ele_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleStruct::type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void EleStruct::set_type(const std::string &value) {
  ele_struct_set_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleStruct::alias() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_alias_info);
  return std::string(arr.data(), arr.size());
}
void EleStruct::set_alias(const std::string &value) {
  ele_struct_set_alias(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleStruct::component_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, ele_struct_get_component_name_info);
  return std::string(arr.data(), arr.size());
}
void EleStruct::set_component_name(const std::string &value) {
  ele_struct_set_component_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string EleStruct::descrip() const {
  return ProxyHelpers::get_string(fortran_ptr_, ele_struct_get_descrip_info);
}
void EleStruct::set_descrip(const std::string &value) {
  ele_struct_set_descrip(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TwissStruct EleStruct::a() const {
  void *ptr;
  ele_struct_get_a(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void EleStruct::set_a(const TwissStruct &src) {
  ele_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct EleStruct::b() const {
  void *ptr;
  ele_struct_get_b(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void EleStruct::set_b(const TwissStruct &src) {
  ele_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct EleStruct::z() const {
  void *ptr;
  ele_struct_get_z(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void EleStruct::set_z(const TwissStruct &src) {
  ele_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
XyDispStruct EleStruct::x() const {
  void *ptr;
  ele_struct_get_x(fortran_ptr_, &ptr);
  return XyDispStruct(ptr);
}
void EleStruct::set_x(const XyDispStruct &src) {
  ele_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
XyDispStruct EleStruct::y() const {
  void *ptr;
  ele_struct_get_y(fortran_ptr_, &ptr);
  return XyDispStruct(ptr);
}
void EleStruct::set_y(const XyDispStruct &src) {
  ele_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<AcKickerStruct> EleStruct::ac_kick() const {
  void *ptr;
  ele_struct_get_ac_kick(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return AcKickerStruct(ptr);
}
void EleStruct::set_ac_kick(const AcKickerStruct &src) {
  ele_struct_set_ac_kick(fortran_ptr_, src.get_fortran_ptr());
}
BookkeepingStateStruct EleStruct::bookkeeping_state() const {
  void *ptr;
  ele_struct_get_bookkeeping_state(fortran_ptr_, &ptr);
  return BookkeepingStateStruct(ptr);
}
void EleStruct::set_bookkeeping_state(const BookkeepingStateStruct &src) {
  ele_struct_set_bookkeeping_state(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<BranchStruct> EleStruct::branch() const {
  void *ptr;
  ele_struct_get_branch(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return BranchStruct(ptr);
}
void EleStruct::set_branch(const BranchStruct &src) {
  ele_struct_set_branch(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<ControllerStruct> EleStruct::control() const {
  void *ptr;
  ele_struct_get_control(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ControllerStruct(ptr);
}
void EleStruct::set_control(const ControllerStruct &src) {
  ele_struct_set_control(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<RfEleStruct> EleStruct::rf() const {
  void *ptr;
  ele_struct_get_rf(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return RfEleStruct(ptr);
}
void EleStruct::set_rf(const RfEleStruct &src) {
  ele_struct_set_rf(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<EleStruct> EleStruct::lord() const {
  void *ptr;
  ele_struct_get_lord(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleStruct(ptr);
}
void EleStruct::set_lord(const EleStruct &src) {
  ele_struct_set_lord(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<Fibre> EleStruct::ptc_fibre() const {
  void *ptr;
  ele_struct_get_ptc_fibre(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void EleStruct::set_ptc_fibre(const Fibre &src) {
  ele_struct_set_ptc_fibre(fortran_ptr_, src.get_fortran_ptr());
}
FloorPositionStruct EleStruct::floor() const {
  void *ptr;
  ele_struct_get_floor(fortran_ptr_, &ptr);
  return FloorPositionStruct(ptr);
}
void EleStruct::set_floor(const FloorPositionStruct &src) {
  ele_struct_set_floor(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<HighEnergySpaceChargeStruct> EleStruct::high_energy_space_charge() const {
  void *ptr;
  ele_struct_get_high_energy_space_charge(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return HighEnergySpaceChargeStruct(ptr);
}
void EleStruct::set_high_energy_space_charge(const HighEnergySpaceChargeStruct &src) {
  ele_struct_set_high_energy_space_charge(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<Mode3Struct> EleStruct::mode3() const {
  void *ptr;
  ele_struct_get_mode3(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Mode3Struct(ptr);
}
void EleStruct::set_mode3(const Mode3Struct &src) {
  ele_struct_set_mode3(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<PhotonElementStruct> EleStruct::photon() const {
  void *ptr;
  ele_struct_get_photon(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return PhotonElementStruct(ptr);
}
void EleStruct::set_photon(const PhotonElementStruct &src) {
  ele_struct_set_photon(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<RadMapEleStruct> EleStruct::rad_map() const {
  void *ptr;
  ele_struct_get_rad_map(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return RadMapEleStruct(ptr);
}
void EleStruct::set_rad_map(const RadMapEleStruct &src) {
  ele_struct_set_rad_map(fortran_ptr_, src.get_fortran_ptr());
}
TaylorStructArray1D EleStruct::taylor() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      ele_struct_get_taylor_info
  );
}
FArray1D<double> EleStruct::spin_taylor_ref_orb_in() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      ele_struct_get_spin_taylor_ref_orb_in_info
  );
}
void EleStruct::set_spin_taylor_ref_orb_in(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_spin_taylor_ref_orb_in(fortran_ptr_, v.data(), shape);
}
TaylorStructArray1D EleStruct::spin_taylor() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      ele_struct_get_spin_taylor_info
  );
}
std::optional<WakeStruct> EleStruct::wake() const {
  void *ptr;
  ele_struct_get_wake(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return WakeStruct(ptr);
}
void EleStruct::set_wake(const WakeStruct &src) {
  ele_struct_set_wake(fortran_ptr_, src.get_fortran_ptr());
}
Wall3dStructArray1D EleStruct::wall3d() const {
  return ProxyHelpers::get_type_array_1d<Wall3dStructArray1D>(
      fortran_ptr_,
      ele_struct_get_wall3d_info
  );
}
CartesianMapStructArray1D EleStruct::cartesian_map() const {
  return ProxyHelpers::get_type_array_1d<CartesianMapStructArray1D>(
      fortran_ptr_,
      ele_struct_get_cartesian_map_info
  );
}
CylindricalMapStructArray1D EleStruct::cylindrical_map() const {
  return ProxyHelpers::get_type_array_1d<CylindricalMapStructArray1D>(
      fortran_ptr_,
      ele_struct_get_cylindrical_map_info
  );
}
GenGradMapStructArray1D EleStruct::gen_grad_map() const {
  return ProxyHelpers::get_type_array_1d<GenGradMapStructArray1D>(
      fortran_ptr_,
      ele_struct_get_gen_grad_map_info
  );
}
GridFieldStructArray1D EleStruct::grid_field() const {
  return ProxyHelpers::get_type_array_1d<GridFieldStructArray1D>(
      fortran_ptr_,
      ele_struct_get_grid_field_info
  );
}
CoordStruct EleStruct::map_ref_orb_in() const {
  void *ptr;
  ele_struct_get_map_ref_orb_in(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void EleStruct::set_map_ref_orb_in(const CoordStruct &src) {
  ele_struct_set_map_ref_orb_in(fortran_ptr_, src.get_fortran_ptr());
}
CoordStruct EleStruct::map_ref_orb_out() const {
  void *ptr;
  ele_struct_get_map_ref_orb_out(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void EleStruct::set_map_ref_orb_out(const CoordStruct &src) {
  ele_struct_set_map_ref_orb_out(fortran_ptr_, src.get_fortran_ptr());
}
CoordStruct EleStruct::time_ref_orb_in() const {
  void *ptr;
  ele_struct_get_time_ref_orb_in(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void EleStruct::set_time_ref_orb_in(const CoordStruct &src) {
  ele_struct_set_time_ref_orb_in(fortran_ptr_, src.get_fortran_ptr());
}
CoordStruct EleStruct::time_ref_orb_out() const {
  void *ptr;
  ele_struct_get_time_ref_orb_out(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void EleStruct::set_time_ref_orb_out(const CoordStruct &src) {
  ele_struct_set_time_ref_orb_out(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> EleStruct::value() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_value_info);
}
void EleStruct::set_value(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_value(fortran_ptr_, v.data(), shape);
}
FArray1D<double> EleStruct::old_value() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_old_value_info);
}
void EleStruct::set_old_value(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_old_value(fortran_ptr_, v.data(), shape);
}
FArray2D<double> EleStruct::spin_q() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, ele_struct_get_spin_q_info);
}
void EleStruct::set_spin_q(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, ele_struct_set_spin_q, v);
}
FArray1D<double> EleStruct::vec0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_vec0_info);
}
void EleStruct::set_vec0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_vec0(fortran_ptr_, v.data(), shape);
}
FArray2D<double> EleStruct::mat6() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, ele_struct_get_mat6_info);
}
void EleStruct::set_mat6(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, ele_struct_set_mat6, v);
}
FArray2D<double> EleStruct::c_mat() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, ele_struct_get_c_mat_info);
}
void EleStruct::set_c_mat(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, ele_struct_set_c_mat, v);
}
FArray2D<double> EleStruct::dc_mat_dpz() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, ele_struct_get_dc_mat_dpz_info);
}
void EleStruct::set_dc_mat_dpz(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, ele_struct_set_dc_mat_dpz, v);
}
double EleStruct::gamma_c() const {
  double value;
  ele_struct_get_gamma_c(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_gamma_c(double value) { ele_struct_set_gamma_c(fortran_ptr_, value); }
double EleStruct::s_start() const {
  double value;
  ele_struct_get_s_start(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_s_start(double value) { ele_struct_set_s_start(fortran_ptr_, value); }
double EleStruct::s() const {
  double value;
  ele_struct_get_s(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_s(double value) { ele_struct_set_s(fortran_ptr_, value); }
double EleStruct::ref_time() const {
  double value;
  ele_struct_get_ref_time(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ref_time(double value) { ele_struct_set_ref_time(fortran_ptr_, value); }
FArray1D<double> EleStruct::a_pole() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_a_pole_info);
}
void EleStruct::set_a_pole(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_a_pole(fortran_ptr_, v.data(), shape);
}
FArray1D<double> EleStruct::b_pole() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_b_pole_info);
}
void EleStruct::set_b_pole(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_b_pole(fortran_ptr_, v.data(), shape);
}
FArray1D<double> EleStruct::a_pole_elec() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_a_pole_elec_info);
}
void EleStruct::set_a_pole_elec(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_a_pole_elec(fortran_ptr_, v.data(), shape);
}
FArray1D<double> EleStruct::b_pole_elec() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_b_pole_elec_info);
}
void EleStruct::set_b_pole_elec(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_b_pole_elec(fortran_ptr_, v.data(), shape);
}
FArray1D<double> EleStruct::custom() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ele_struct_get_custom_info);
}
void EleStruct::set_custom(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ele_struct_set_custom(fortran_ptr_, v.data(), shape);
}
FArray3D<double> EleStruct::r() const {
  return ProxyHelpers::get_array_3d<double>(fortran_ptr_, ele_struct_get_r_info);
}
void EleStruct::set_r(const std::vector<std::vector<std::vector<double>>> &v) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, ele_struct_set_r, v);
}
int EleStruct::key() const {
  int value;
  ele_struct_get_key(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_key(int value) { ele_struct_set_key(fortran_ptr_, value); }
int EleStruct::sub_key() const {
  int value;
  ele_struct_get_sub_key(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_sub_key(int value) { ele_struct_set_sub_key(fortran_ptr_, value); }
int EleStruct::ix_ele() const {
  int value;
  ele_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ix_ele(int value) { ele_struct_set_ix_ele(fortran_ptr_, value); }
int EleStruct::ix_branch() const {
  int value;
  ele_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ix_branch(int value) { ele_struct_set_ix_branch(fortran_ptr_, value); }
int EleStruct::lord_status() const {
  int value;
  ele_struct_get_lord_status(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_lord_status(int value) { ele_struct_set_lord_status(fortran_ptr_, value); }
int EleStruct::n_slave() const {
  int value;
  ele_struct_get_n_slave(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_n_slave(int value) { ele_struct_set_n_slave(fortran_ptr_, value); }
int EleStruct::n_slave_field() const {
  int value;
  ele_struct_get_n_slave_field(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_n_slave_field(int value) { ele_struct_set_n_slave_field(fortran_ptr_, value); }
int EleStruct::ix1_slave() const {
  int value;
  ele_struct_get_ix1_slave(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ix1_slave(int value) { ele_struct_set_ix1_slave(fortran_ptr_, value); }
int EleStruct::slave_status() const {
  int value;
  ele_struct_get_slave_status(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_slave_status(int value) { ele_struct_set_slave_status(fortran_ptr_, value); }
int EleStruct::n_lord() const {
  int value;
  ele_struct_get_n_lord(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_n_lord(int value) { ele_struct_set_n_lord(fortran_ptr_, value); }
int EleStruct::n_lord_field() const {
  int value;
  ele_struct_get_n_lord_field(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_n_lord_field(int value) { ele_struct_set_n_lord_field(fortran_ptr_, value); }
int EleStruct::n_lord_ramper() const {
  int value;
  ele_struct_get_n_lord_ramper(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_n_lord_ramper(int value) { ele_struct_set_n_lord_ramper(fortran_ptr_, value); }
int EleStruct::ic1_lord() const {
  int value;
  ele_struct_get_ic1_lord(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ic1_lord(int value) { ele_struct_set_ic1_lord(fortran_ptr_, value); }
int EleStruct::ix_pointer() const {
  int value;
  ele_struct_get_ix_pointer(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ix_pointer(int value) { ele_struct_set_ix_pointer(fortran_ptr_, value); }
int EleStruct::ixx() const {
  int value;
  ele_struct_get_ixx(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ixx(int value) { ele_struct_set_ixx(fortran_ptr_, value); }
int EleStruct::iyy() const {
  int value;
  ele_struct_get_iyy(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_iyy(int value) { ele_struct_set_iyy(fortran_ptr_, value); }
int EleStruct::izz() const {
  int value;
  ele_struct_get_izz(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_izz(int value) { ele_struct_set_izz(fortran_ptr_, value); }
int EleStruct::mat6_calc_method() const {
  int value;
  ele_struct_get_mat6_calc_method(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_mat6_calc_method(int value) {
  ele_struct_set_mat6_calc_method(fortran_ptr_, value);
}
int EleStruct::tracking_method() const {
  int value;
  ele_struct_get_tracking_method(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_tracking_method(int value) {
  ele_struct_set_tracking_method(fortran_ptr_, value);
}
int EleStruct::spin_tracking_method() const {
  int value;
  ele_struct_get_spin_tracking_method(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_spin_tracking_method(int value) {
  ele_struct_set_spin_tracking_method(fortran_ptr_, value);
}
int EleStruct::csr_method() const {
  int value;
  ele_struct_get_csr_method(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_csr_method(int value) { ele_struct_set_csr_method(fortran_ptr_, value); }
int EleStruct::space_charge_method() const {
  int value;
  ele_struct_get_space_charge_method(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_space_charge_method(int value) {
  ele_struct_set_space_charge_method(fortran_ptr_, value);
}
int EleStruct::ptc_integration_type() const {
  int value;
  ele_struct_get_ptc_integration_type(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ptc_integration_type(int value) {
  ele_struct_set_ptc_integration_type(fortran_ptr_, value);
}
int EleStruct::field_calc() const {
  int value;
  ele_struct_get_field_calc(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_field_calc(int value) { ele_struct_set_field_calc(fortran_ptr_, value); }
int EleStruct::aperture_at() const {
  int value;
  ele_struct_get_aperture_at(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_aperture_at(int value) { ele_struct_set_aperture_at(fortran_ptr_, value); }
int EleStruct::aperture_type() const {
  int value;
  ele_struct_get_aperture_type(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_aperture_type(int value) { ele_struct_set_aperture_type(fortran_ptr_, value); }
int EleStruct::ref_species() const {
  int value;
  ele_struct_get_ref_species(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_ref_species(int value) { ele_struct_set_ref_species(fortran_ptr_, value); }
int EleStruct::orientation() const {
  int value;
  ele_struct_get_orientation(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_orientation(int value) { ele_struct_set_orientation(fortran_ptr_, value); }
bool EleStruct::symplectify() const {
  bool value;
  ele_struct_get_symplectify(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_symplectify(bool value) { ele_struct_set_symplectify(fortran_ptr_, value); }
bool EleStruct::mode_flip() const {
  bool value;
  ele_struct_get_mode_flip(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_mode_flip(bool value) { ele_struct_set_mode_flip(fortran_ptr_, value); }
bool EleStruct::multipoles_on() const {
  bool value;
  ele_struct_get_multipoles_on(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_multipoles_on(bool value) { ele_struct_set_multipoles_on(fortran_ptr_, value); }
bool EleStruct::scale_multipoles() const {
  bool value;
  ele_struct_get_scale_multipoles(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_scale_multipoles(bool value) {
  ele_struct_set_scale_multipoles(fortran_ptr_, value);
}
bool EleStruct::taylor_map_includes_offsets() const {
  bool value;
  ele_struct_get_taylor_map_includes_offsets(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_taylor_map_includes_offsets(bool value) {
  ele_struct_set_taylor_map_includes_offsets(fortran_ptr_, value);
}
bool EleStruct::field_master() const {
  bool value;
  ele_struct_get_field_master(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_field_master(bool value) { ele_struct_set_field_master(fortran_ptr_, value); }
bool EleStruct::is_on() const {
  bool value;
  ele_struct_get_is_on(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_is_on(bool value) { ele_struct_set_is_on(fortran_ptr_, value); }
bool EleStruct::logic() const {
  bool value;
  ele_struct_get_logic(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_logic(bool value) { ele_struct_set_logic(fortran_ptr_, value); }
bool EleStruct::bmad_logic() const {
  bool value;
  ele_struct_get_bmad_logic(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_bmad_logic(bool value) { ele_struct_set_bmad_logic(fortran_ptr_, value); }
bool EleStruct::select() const {
  bool value;
  ele_struct_get_select(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_select(bool value) { ele_struct_set_select(fortran_ptr_, value); }
bool EleStruct::offset_moves_aperture() const {
  bool value;
  ele_struct_get_offset_moves_aperture(fortran_ptr_, &value);
  return value;
}
void EleStruct::set_offset_moves_aperture(bool value) {
  ele_struct_set_offset_moves_aperture(fortran_ptr_, value);
}
std::complex<double> ComplexTaylorTermStruct::coef() const {
  std::complex<double> c_value;
  complex_taylor_term_struct_get_coef(fortran_ptr_, &c_value);
  return c_value;
}
void ComplexTaylorTermStruct::set_coef(std::complex<double> value) {
  complex_taylor_term_struct_set_coef(fortran_ptr_, value);
}
FArray1D<int> ComplexTaylorTermStruct::expn() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, complex_taylor_term_struct_get_expn_info);
}
void ComplexTaylorTermStruct::set_expn(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  complex_taylor_term_struct_set_expn(fortran_ptr_, v.data(), shape);
}
std::complex<double> ComplexTaylorStruct::ref() const {
  std::complex<double> c_value;
  complex_taylor_struct_get_ref(fortran_ptr_, &c_value);
  return c_value;
}
void ComplexTaylorStruct::set_ref(std::complex<double> value) {
  complex_taylor_struct_set_ref(fortran_ptr_, value);
}
ComplexTaylorTermStructArray1D ComplexTaylorStruct::term() const {
  return ProxyHelpers::get_type_array_1d<ComplexTaylorTermStructArray1D>(
      fortran_ptr_,
      complex_taylor_struct_get_term_info
  );
}
std::string BranchStruct::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, branch_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void BranchStruct::set_name(const std::string &value) {
  branch_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int BranchStruct::ix_branch() const {
  int value;
  branch_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_ix_branch(int value) { branch_struct_set_ix_branch(fortran_ptr_, value); }
int BranchStruct::ix_from_branch() const {
  int value;
  branch_struct_get_ix_from_branch(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_ix_from_branch(int value) {
  branch_struct_set_ix_from_branch(fortran_ptr_, value);
}
int BranchStruct::ix_from_ele() const {
  int value;
  branch_struct_get_ix_from_ele(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_ix_from_ele(int value) {
  branch_struct_set_ix_from_ele(fortran_ptr_, value);
}
int BranchStruct::ix_to_ele() const {
  int value;
  branch_struct_get_ix_to_ele(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_ix_to_ele(int value) { branch_struct_set_ix_to_ele(fortran_ptr_, value); }
int BranchStruct::ix_fixer() const {
  int value;
  branch_struct_get_ix_fixer(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_ix_fixer(int value) { branch_struct_set_ix_fixer(fortran_ptr_, value); }
int BranchStruct::n_ele_track() const {
  int value;
  branch_struct_get_n_ele_track(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_n_ele_track(int value) {
  branch_struct_set_n_ele_track(fortran_ptr_, value);
}
int BranchStruct::n_ele_max() const {
  int value;
  branch_struct_get_n_ele_max(fortran_ptr_, &value);
  return value;
}
void BranchStruct::set_n_ele_max(int value) { branch_struct_set_n_ele_max(fortran_ptr_, value); }
std::optional<LatStruct> BranchStruct::lat() const {
  void *ptr;
  branch_struct_get_lat(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return LatStruct(ptr);
}
void BranchStruct::set_lat(const LatStruct &src) {
  branch_struct_set_lat(fortran_ptr_, src.get_fortran_ptr());
}
ModeInfoStruct BranchStruct::a() const {
  void *ptr;
  branch_struct_get_a(fortran_ptr_, &ptr);
  return ModeInfoStruct(ptr);
}
void BranchStruct::set_a(const ModeInfoStruct &src) {
  branch_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
ModeInfoStruct BranchStruct::b() const {
  void *ptr;
  branch_struct_get_b(fortran_ptr_, &ptr);
  return ModeInfoStruct(ptr);
}
void BranchStruct::set_b(const ModeInfoStruct &src) {
  branch_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
ModeInfoStruct BranchStruct::z() const {
  void *ptr;
  branch_struct_get_z(fortran_ptr_, &ptr);
  return ModeInfoStruct(ptr);
}
void BranchStruct::set_z(const ModeInfoStruct &src) {
  branch_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
EleStructArray1D BranchStruct::ele() const {
  return ProxyHelpers::get_type_array_1d<EleStructArray1D>(
      fortran_ptr_,
      branch_struct_get_ele_info
  );
}
LatParamStruct BranchStruct::param() const {
  void *ptr;
  branch_struct_get_param(fortran_ptr_, &ptr);
  return LatParamStruct(ptr);
}
void BranchStruct::set_param(const LatParamStruct &src) {
  branch_struct_set_param(fortran_ptr_, src.get_fortran_ptr());
}
CoordStruct BranchStruct::particle_start() const {
  void *ptr;
  branch_struct_get_particle_start(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void BranchStruct::set_particle_start(const CoordStruct &src) {
  branch_struct_set_particle_start(fortran_ptr_, src.get_fortran_ptr());
}
Wall3dStructArray1D BranchStruct::wall3d() const {
  return ProxyHelpers::get_type_array_1d<Wall3dStructArray1D>(
      fortran_ptr_,
      branch_struct_get_wall3d_info
  );
}
std::string LatStruct::use_name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, lat_struct_get_use_name_info);
  return std::string(arr.data(), arr.size());
}
void LatStruct::set_use_name(const std::string &value) {
  lat_struct_set_use_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatStruct::lattice() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, lat_struct_get_lattice_info);
  return std::string(arr.data(), arr.size());
}
void LatStruct::set_lattice(const std::string &value) {
  lat_struct_set_lattice(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatStruct::machine() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, lat_struct_get_machine_info);
  return std::string(arr.data(), arr.size());
}
void LatStruct::set_machine(const std::string &value) {
  lat_struct_set_machine(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatStruct::input_file_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, lat_struct_get_input_file_name_info);
  return std::string(arr.data(), arr.size());
}
void LatStruct::set_input_file_name(const std::string &value) {
  lat_struct_set_input_file_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string LatStruct::title() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, lat_struct_get_title_info);
  return std::string(arr.data(), arr.size());
}
void LatStruct::set_title(const std::string &value) {
  lat_struct_set_title(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D LatStruct::print_str() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, lat_struct_get_print_str_info);
}
ExpressionAtomStructAlloc1D LatStruct::constant() const {
  return ExpressionAtomStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_struct_reallocate_constant,
      lat_struct_get_constant_info
  );
}
std::optional<ModeInfoStruct> LatStruct::a() const {
  void *ptr;
  lat_struct_get_a(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ModeInfoStruct(ptr);
}
void LatStruct::set_a(const ModeInfoStruct &src) {
  lat_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<ModeInfoStruct> LatStruct::b() const {
  void *ptr;
  lat_struct_get_b(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ModeInfoStruct(ptr);
}
void LatStruct::set_b(const ModeInfoStruct &src) {
  lat_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<ModeInfoStruct> LatStruct::z() const {
  void *ptr;
  lat_struct_get_z(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return ModeInfoStruct(ptr);
}
void LatStruct::set_z(const ModeInfoStruct &src) {
  lat_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<LatParamStruct> LatStruct::param() const {
  void *ptr;
  lat_struct_get_param(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return LatParamStruct(ptr);
}
void LatStruct::set_param(const LatParamStruct &src) {
  lat_struct_set_param(fortran_ptr_, src.get_fortran_ptr());
}
BookkeepingStateStruct LatStruct::lord_state() const {
  void *ptr;
  lat_struct_get_lord_state(fortran_ptr_, &ptr);
  return BookkeepingStateStruct(ptr);
}
void LatStruct::set_lord_state(const BookkeepingStateStruct &src) {
  lat_struct_set_lord_state(fortran_ptr_, src.get_fortran_ptr());
}
EleStruct LatStruct::ele_init() const {
  void *ptr;
  lat_struct_get_ele_init(fortran_ptr_, &ptr);
  return EleStruct(ptr);
}
void LatStruct::set_ele_init(const EleStruct &src) {
  lat_struct_set_ele_init(fortran_ptr_, src.get_fortran_ptr());
}
EleStructArray1D LatStruct::ele() const {
  return ProxyHelpers::get_type_array_1d<EleStructArray1D>(fortran_ptr_, lat_struct_get_ele_info);
}
BranchStructAlloc1D LatStruct::branch() const {
  return BranchStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_struct_reallocate_branch,
      lat_struct_get_branch_info
  );
}
ControlStructAlloc1D LatStruct::control() const {
  return ControlStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_struct_reallocate_control,
      lat_struct_get_control_info
  );
}
std::optional<CoordStruct> LatStruct::particle_start() const {
  void *ptr;
  lat_struct_get_particle_start(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return CoordStruct(ptr);
}
void LatStruct::set_particle_start(const CoordStruct &src) {
  lat_struct_set_particle_start(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitStruct LatStruct::beam_init() const {
  void *ptr;
  lat_struct_get_beam_init(fortran_ptr_, &ptr);
  return BeamInitStruct(ptr);
}
void LatStruct::set_beam_init(const BeamInitStruct &src) {
  lat_struct_set_beam_init(fortran_ptr_, src.get_fortran_ptr());
}
PreTrackerStruct LatStruct::pre_tracker() const {
  void *ptr;
  lat_struct_get_pre_tracker(fortran_ptr_, &ptr);
  return PreTrackerStruct(ptr);
}
void LatStruct::set_pre_tracker(const PreTrackerStruct &src) {
  lat_struct_set_pre_tracker(fortran_ptr_, src.get_fortran_ptr());
}
RealAlloc1D LatStruct::custom() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_struct_reallocate_custom,
      lat_struct_get_custom_info
  );
}
void LatStruct::set_custom(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  lat_struct_set_custom(fortran_ptr_, v.data(), shape);
}
int LatStruct::version() const {
  int value;
  lat_struct_get_version(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_version(int value) { lat_struct_set_version(fortran_ptr_, value); }
std::optional<int> LatStruct::n_ele_track() const {
  int value;
  bool is_valid;
  lat_struct_get_n_ele_track(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void LatStruct::set_n_ele_track(int value) { lat_struct_set_n_ele_track(fortran_ptr_, value); }
std::optional<int> LatStruct::n_ele_max() const {
  int value;
  bool is_valid;
  lat_struct_get_n_ele_max(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void LatStruct::set_n_ele_max(int value) { lat_struct_set_n_ele_max(fortran_ptr_, value); }
int LatStruct::n_control_max() const {
  int value;
  lat_struct_get_n_control_max(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_n_control_max(int value) { lat_struct_set_n_control_max(fortran_ptr_, value); }
int LatStruct::n_ic_max() const {
  int value;
  lat_struct_get_n_ic_max(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_n_ic_max(int value) { lat_struct_set_n_ic_max(fortran_ptr_, value); }
int LatStruct::input_taylor_order() const {
  int value;
  lat_struct_get_input_taylor_order(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_input_taylor_order(int value) {
  lat_struct_set_input_taylor_order(fortran_ptr_, value);
}
IntAlloc1D LatStruct::ic() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_struct_reallocate_ic,
      lat_struct_get_ic_info
  );
}
void LatStruct::set_ic(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  lat_struct_set_ic(fortran_ptr_, v.data(), shape);
}
int LatStruct::photon_type() const {
  int value;
  lat_struct_get_photon_type(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_photon_type(int value) { lat_struct_set_photon_type(fortran_ptr_, value); }
int LatStruct::creation_hash() const {
  int value;
  lat_struct_get_creation_hash(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_creation_hash(int value) { lat_struct_set_creation_hash(fortran_ptr_, value); }
int LatStruct::ramper_slave_bookkeeping() const {
  int value;
  lat_struct_get_ramper_slave_bookkeeping(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_ramper_slave_bookkeeping(int value) {
  lat_struct_set_ramper_slave_bookkeeping(fortran_ptr_, value);
}
bool LatStruct::parser_make_xfer_mats() const {
  bool value;
  lat_struct_get_parser_make_xfer_mats(fortran_ptr_, &value);
  return value;
}
void LatStruct::set_parser_make_xfer_mats(bool value) {
  lat_struct_set_parser_make_xfer_mats(fortran_ptr_, value);
}
CoordStructAlloc1D BunchStruct::particle() const {
  return CoordStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bunch_struct_reallocate_particle,
      bunch_struct_get_particle_info
  );
}
IntAlloc1D BunchStruct::ix_z() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bunch_struct_reallocate_ix_z,
      bunch_struct_get_ix_z_info
  );
}
void BunchStruct::set_ix_z(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bunch_struct_set_ix_z(fortran_ptr_, v.data(), shape);
}
double BunchStruct::charge_tot() const {
  double value;
  bunch_struct_get_charge_tot(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_charge_tot(double value) { bunch_struct_set_charge_tot(fortran_ptr_, value); }
double BunchStruct::charge_live() const {
  double value;
  bunch_struct_get_charge_live(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_charge_live(double value) {
  bunch_struct_set_charge_live(fortran_ptr_, value);
}
double BunchStruct::z_center() const {
  double value;
  bunch_struct_get_z_center(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_z_center(double value) { bunch_struct_set_z_center(fortran_ptr_, value); }
double BunchStruct::t_center() const {
  double value;
  bunch_struct_get_t_center(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_t_center(double value) { bunch_struct_set_t_center(fortran_ptr_, value); }
double BunchStruct::t0() const {
  double value;
  bunch_struct_get_t0(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_t0(double value) { bunch_struct_set_t0(fortran_ptr_, value); }
bool BunchStruct::drift_between_t_and_s() const {
  bool value;
  bunch_struct_get_drift_between_t_and_s(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_drift_between_t_and_s(bool value) {
  bunch_struct_set_drift_between_t_and_s(fortran_ptr_, value);
}
int BunchStruct::ix_ele() const {
  int value;
  bunch_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_ix_ele(int value) { bunch_struct_set_ix_ele(fortran_ptr_, value); }
int BunchStruct::ix_bunch() const {
  int value;
  bunch_struct_get_ix_bunch(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_ix_bunch(int value) { bunch_struct_set_ix_bunch(fortran_ptr_, value); }
int BunchStruct::ix_turn() const {
  int value;
  bunch_struct_get_ix_turn(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_ix_turn(int value) { bunch_struct_set_ix_turn(fortran_ptr_, value); }
int BunchStruct::n_live() const {
  int value;
  bunch_struct_get_n_live(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_n_live(int value) { bunch_struct_set_n_live(fortran_ptr_, value); }
int BunchStruct::n_good() const {
  int value;
  bunch_struct_get_n_good(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_n_good(int value) { bunch_struct_set_n_good(fortran_ptr_, value); }
int BunchStruct::n_bad() const {
  int value;
  bunch_struct_get_n_bad(fortran_ptr_, &value);
  return value;
}
void BunchStruct::set_n_bad(int value) { bunch_struct_set_n_bad(fortran_ptr_, value); }
CoordStruct BunchParamsStruct::centroid() const {
  void *ptr;
  bunch_params_struct_get_centroid(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void BunchParamsStruct::set_centroid(const CoordStruct &src) {
  bunch_params_struct_set_centroid(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct BunchParamsStruct::x() const {
  void *ptr;
  bunch_params_struct_get_x(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void BunchParamsStruct::set_x(const TwissStruct &src) {
  bunch_params_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct BunchParamsStruct::y() const {
  void *ptr;
  bunch_params_struct_get_y(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void BunchParamsStruct::set_y(const TwissStruct &src) {
  bunch_params_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct BunchParamsStruct::z() const {
  void *ptr;
  bunch_params_struct_get_z(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void BunchParamsStruct::set_z(const TwissStruct &src) {
  bunch_params_struct_set_z(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct BunchParamsStruct::a() const {
  void *ptr;
  bunch_params_struct_get_a(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void BunchParamsStruct::set_a(const TwissStruct &src) {
  bunch_params_struct_set_a(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct BunchParamsStruct::b() const {
  void *ptr;
  bunch_params_struct_get_b(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void BunchParamsStruct::set_b(const TwissStruct &src) {
  bunch_params_struct_set_b(fortran_ptr_, src.get_fortran_ptr());
}
TwissStruct BunchParamsStruct::c() const {
  void *ptr;
  bunch_params_struct_get_c(fortran_ptr_, &ptr);
  return TwissStruct(ptr);
}
void BunchParamsStruct::set_c(const TwissStruct &src) {
  bunch_params_struct_set_c(fortran_ptr_, src.get_fortran_ptr());
}
FArray2D<double> BunchParamsStruct::sigma() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, bunch_params_struct_get_sigma_info);
}
void BunchParamsStruct::set_sigma(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, bunch_params_struct_set_sigma, v);
}
FArray1D<double> BunchParamsStruct::rel_max() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bunch_params_struct_get_rel_max_info);
}
void BunchParamsStruct::set_rel_max(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bunch_params_struct_set_rel_max(fortran_ptr_, v.data(), shape);
}
FArray1D<double> BunchParamsStruct::rel_min() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bunch_params_struct_get_rel_min_info);
}
void BunchParamsStruct::set_rel_min(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bunch_params_struct_set_rel_min(fortran_ptr_, v.data(), shape);
}
double BunchParamsStruct::s() const {
  double value;
  bunch_params_struct_get_s(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_s(double value) { bunch_params_struct_set_s(fortran_ptr_, value); }
double BunchParamsStruct::t() const {
  double value;
  bunch_params_struct_get_t(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_t(double value) { bunch_params_struct_set_t(fortran_ptr_, value); }
double BunchParamsStruct::sigma_t() const {
  double value;
  bunch_params_struct_get_sigma_t(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_sigma_t(double value) {
  bunch_params_struct_set_sigma_t(fortran_ptr_, value);
}
double BunchParamsStruct::charge_live() const {
  double value;
  bunch_params_struct_get_charge_live(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_charge_live(double value) {
  bunch_params_struct_set_charge_live(fortran_ptr_, value);
}
double BunchParamsStruct::charge_tot() const {
  double value;
  bunch_params_struct_get_charge_tot(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_charge_tot(double value) {
  bunch_params_struct_set_charge_tot(fortran_ptr_, value);
}
int BunchParamsStruct::n_particle_tot() const {
  int value;
  bunch_params_struct_get_n_particle_tot(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_n_particle_tot(int value) {
  bunch_params_struct_set_n_particle_tot(fortran_ptr_, value);
}
int BunchParamsStruct::n_particle_live() const {
  int value;
  bunch_params_struct_get_n_particle_live(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_n_particle_live(int value) {
  bunch_params_struct_set_n_particle_live(fortran_ptr_, value);
}
int BunchParamsStruct::n_particle_lost_in_ele() const {
  int value;
  bunch_params_struct_get_n_particle_lost_in_ele(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_n_particle_lost_in_ele(int value) {
  bunch_params_struct_set_n_particle_lost_in_ele(fortran_ptr_, value);
}
int BunchParamsStruct::n_good_steps() const {
  int value;
  bunch_params_struct_get_n_good_steps(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_n_good_steps(int value) {
  bunch_params_struct_set_n_good_steps(fortran_ptr_, value);
}
int BunchParamsStruct::n_bad_steps() const {
  int value;
  bunch_params_struct_get_n_bad_steps(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_n_bad_steps(int value) {
  bunch_params_struct_set_n_bad_steps(fortran_ptr_, value);
}
int BunchParamsStruct::ix_ele() const {
  int value;
  bunch_params_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_ix_ele(int value) {
  bunch_params_struct_set_ix_ele(fortran_ptr_, value);
}
int BunchParamsStruct::location() const {
  int value;
  bunch_params_struct_get_location(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_location(int value) {
  bunch_params_struct_set_location(fortran_ptr_, value);
}
bool BunchParamsStruct::twiss_valid() const {
  bool value;
  bunch_params_struct_get_twiss_valid(fortran_ptr_, &value);
  return value;
}
void BunchParamsStruct::set_twiss_valid(bool value) {
  bunch_params_struct_set_twiss_valid(fortran_ptr_, value);
}
BunchStructAlloc1D BeamStruct::bunch() const {
  return BunchStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      beam_struct_reallocate_bunch,
      beam_struct_get_bunch_info
  );
}
double AperturePointStruct::x() const {
  double value;
  aperture_point_struct_get_x(fortran_ptr_, &value);
  return value;
}
void AperturePointStruct::set_x(double value) { aperture_point_struct_set_x(fortran_ptr_, value); }
double AperturePointStruct::y() const {
  double value;
  aperture_point_struct_get_y(fortran_ptr_, &value);
  return value;
}
void AperturePointStruct::set_y(double value) { aperture_point_struct_set_y(fortran_ptr_, value); }
int AperturePointStruct::plane() const {
  int value;
  aperture_point_struct_get_plane(fortran_ptr_, &value);
  return value;
}
void AperturePointStruct::set_plane(int value) {
  aperture_point_struct_set_plane(fortran_ptr_, value);
}
int AperturePointStruct::ix_ele() const {
  int value;
  aperture_point_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void AperturePointStruct::set_ix_ele(int value) {
  aperture_point_struct_set_ix_ele(fortran_ptr_, value);
}
int AperturePointStruct::i_turn() const {
  int value;
  aperture_point_struct_get_i_turn(fortran_ptr_, &value);
  return value;
}
void AperturePointStruct::set_i_turn(int value) {
  aperture_point_struct_set_i_turn(fortran_ptr_, value);
}
double ApertureParamStruct::min_angle() const {
  double value;
  aperture_param_struct_get_min_angle(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_min_angle(double value) {
  aperture_param_struct_set_min_angle(fortran_ptr_, value);
}
double ApertureParamStruct::max_angle() const {
  double value;
  aperture_param_struct_get_max_angle(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_max_angle(double value) {
  aperture_param_struct_set_max_angle(fortran_ptr_, value);
}
int ApertureParamStruct::n_angle() const {
  int value;
  aperture_param_struct_get_n_angle(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_n_angle(int value) {
  aperture_param_struct_set_n_angle(fortran_ptr_, value);
}
int ApertureParamStruct::n_turn() const {
  int value;
  aperture_param_struct_get_n_turn(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_n_turn(int value) {
  aperture_param_struct_set_n_turn(fortran_ptr_, value);
}
double ApertureParamStruct::x_init() const {
  double value;
  aperture_param_struct_get_x_init(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_x_init(double value) {
  aperture_param_struct_set_x_init(fortran_ptr_, value);
}
double ApertureParamStruct::y_init() const {
  double value;
  aperture_param_struct_get_y_init(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_y_init(double value) {
  aperture_param_struct_set_y_init(fortran_ptr_, value);
}
double ApertureParamStruct::rel_accuracy() const {
  double value;
  aperture_param_struct_get_rel_accuracy(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_rel_accuracy(double value) {
  aperture_param_struct_set_rel_accuracy(fortran_ptr_, value);
}
double ApertureParamStruct::abs_accuracy() const {
  double value;
  aperture_param_struct_get_abs_accuracy(fortran_ptr_, &value);
  return value;
}
void ApertureParamStruct::set_abs_accuracy(double value) {
  aperture_param_struct_set_abs_accuracy(fortran_ptr_, value);
}
std::string ApertureParamStruct::start_ele() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, aperture_param_struct_get_start_ele_info);
  return std::string(arr.data(), arr.size());
}
void ApertureParamStruct::set_start_ele(const std::string &value) {
  aperture_param_struct_set_start_ele(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
AperturePointStructAlloc1D ApertureScanStruct::point() const {
  return AperturePointStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      aperture_scan_struct_reallocate_point,
      aperture_scan_struct_get_point_info
  );
}
CoordStruct ApertureScanStruct::ref_orb() const {
  void *ptr;
  aperture_scan_struct_get_ref_orb(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void ApertureScanStruct::set_ref_orb(const CoordStruct &src) {
  aperture_scan_struct_set_ref_orb(fortran_ptr_, src.get_fortran_ptr());
}
double ApertureScanStruct::pz_start() const {
  double value;
  aperture_scan_struct_get_pz_start(fortran_ptr_, &value);
  return value;
}
void ApertureScanStruct::set_pz_start(double value) {
  aperture_scan_struct_set_pz_start(fortran_ptr_, value);
}
std::optional<EleStruct> ElePointerStruct::ele() const {
  void *ptr;
  ele_pointer_struct_get_ele(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleStruct(ptr);
}
void ElePointerStruct::set_ele(const EleStruct &src) {
  ele_pointer_struct_set_ele(fortran_ptr_, src.get_fortran_ptr());
}
LatEleLocStruct ElePointerStruct::loc() const {
  void *ptr;
  ele_pointer_struct_get_loc(fortran_ptr_, &ptr);
  return LatEleLocStruct(ptr);
}
void ElePointerStruct::set_loc(const LatEleLocStruct &src) {
  ele_pointer_struct_set_loc(fortran_ptr_, src.get_fortran_ptr());
}
int ElePointerStruct::id() const {
  int value;
  ele_pointer_struct_get_id(fortran_ptr_, &value);
  return value;
}
void ElePointerStruct::set_id(int value) { ele_pointer_struct_set_id(fortran_ptr_, value); }
std::string ExpressionTreeStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, expression_tree_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void ExpressionTreeStruct::set_name(const std::string &value) {
  expression_tree_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int ExpressionTreeStruct::type() const {
  int value;
  expression_tree_struct_get_type(fortran_ptr_, &value);
  return value;
}
void ExpressionTreeStruct::set_type(int value) {
  expression_tree_struct_set_type(fortran_ptr_, value);
}
double ExpressionTreeStruct::value() const {
  double value;
  expression_tree_struct_get_value(fortran_ptr_, &value);
  return value;
}
void ExpressionTreeStruct::set_value(double value) {
  expression_tree_struct_set_value(fortran_ptr_, value);
}
ExpressionTreeStructArray1D ExpressionTreeStruct::node() const {
  return ProxyHelpers::get_type_array_1d<ExpressionTreeStructArray1D>(
      fortran_ptr_,
      expression_tree_struct_get_node_info
  );
}
FCharArray1D NametableStruct::name() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, nametable_struct_get_name_info);
}
IntAlloc1D NametableStruct::index() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      nametable_struct_reallocate_index,
      nametable_struct_get_index_info
  );
}
void NametableStruct::set_index(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  nametable_struct_set_index(fortran_ptr_, v.data(), shape);
}
int NametableStruct::n_min() const {
  int value;
  nametable_struct_get_n_min(fortran_ptr_, &value);
  return value;
}
void NametableStruct::set_n_min(int value) { nametable_struct_set_n_min(fortran_ptr_, value); }
int NametableStruct::n_max() const {
  int value;
  nametable_struct_get_n_max(fortran_ptr_, &value);
  return value;
}
void NametableStruct::set_n_max(int value) { nametable_struct_set_n_max(fortran_ptr_, value); }
FArray1D<double> TaoSpinDnDpzStruct::vec() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_spin_dn_dpz_struct_get_vec_info);
}
void TaoSpinDnDpzStruct::set_vec(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_spin_dn_dpz_struct_set_vec(fortran_ptr_, v.data(), shape);
}
FArray2D<double> TaoSpinDnDpzStruct::partial() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_spin_dn_dpz_struct_get_partial_info);
}
void TaoSpinDnDpzStruct::set_partial(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_spin_dn_dpz_struct_set_partial, v);
}
FArray2D<double> TaoSpinDnDpzStruct::partial2() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_spin_dn_dpz_struct_get_partial2_info);
}
void TaoSpinDnDpzStruct::set_partial2(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_spin_dn_dpz_struct_set_partial2, v);
}
std::string ResonanceHStruct::id() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, resonance_h_struct_get_id_info);
  return std::string(arr.data(), arr.size());
}
void ResonanceHStruct::set_id(const std::string &value) {
  resonance_h_struct_set_id(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::complex<double> ResonanceHStruct::c_val() const {
  std::complex<double> c_value;
  resonance_h_struct_get_c_val(fortran_ptr_, &c_value);
  return c_value;
}
void ResonanceHStruct::set_c_val(std::complex<double> value) {
  resonance_h_struct_set_c_val(fortran_ptr_, value);
}
FArray2D<double> SpinOrbitMap1Struct::orb_mat() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, spin_orbit_map1_struct_get_orb_mat_info);
}
void SpinOrbitMap1Struct::set_orb_mat(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, spin_orbit_map1_struct_set_orb_mat, v);
}
FArray1D<double> SpinOrbitMap1Struct::vec0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, spin_orbit_map1_struct_get_vec0_info);
}
void SpinOrbitMap1Struct::set_vec0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  spin_orbit_map1_struct_set_vec0(fortran_ptr_, v.data(), shape);
}
FArray2D<double> SpinOrbitMap1Struct::spin_q() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, spin_orbit_map1_struct_get_spin_q_info);
}
void SpinOrbitMap1Struct::set_spin_q(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, spin_orbit_map1_struct_set_spin_q, v);
}
FArray1D<double> SpinAxisStruct::l() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, spin_axis_struct_get_l_info);
}
void SpinAxisStruct::set_l(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  spin_axis_struct_set_l(fortran_ptr_, v.data(), shape);
}
FArray1D<double> SpinAxisStruct::n0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, spin_axis_struct_get_n0_info);
}
void SpinAxisStruct::set_n0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  spin_axis_struct_set_n0(fortran_ptr_, v.data(), shape);
}
FArray1D<double> SpinAxisStruct::m() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, spin_axis_struct_get_m_info);
}
void SpinAxisStruct::set_m(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  spin_axis_struct_set_m(fortran_ptr_, v.data(), shape);
}
std::optional<EleStruct> PtcNormalFormStruct::ele_origin() const {
  void *ptr;
  ptc_normal_form_struct_get_ele_origin(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleStruct(ptr);
}
void PtcNormalFormStruct::set_ele_origin(const EleStruct &src) {
  ptc_normal_form_struct_set_ele_origin(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> PtcNormalFormStruct::orb0() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, ptc_normal_form_struct_get_orb0_info);
}
void PtcNormalFormStruct::set_orb0(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  ptc_normal_form_struct_set_orb0(fortran_ptr_, v.data(), shape);
}
bool PtcNormalFormStruct::valid_map() const {
  bool value;
  ptc_normal_form_struct_get_valid_map(fortran_ptr_, &value);
  return value;
}
void PtcNormalFormStruct::set_valid_map(bool value) {
  ptc_normal_form_struct_set_valid_map(fortran_ptr_, value);
}
std::optional<EleStruct> BmadNormalFormStruct::ele_origin() const {
  void *ptr;
  bmad_normal_form_struct_get_ele_origin(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleStruct(ptr);
}
void BmadNormalFormStruct::set_ele_origin(const EleStruct &src) {
  bmad_normal_form_struct_set_ele_origin(fortran_ptr_, src.get_fortran_ptr());
}
TaylorStructArray1D BmadNormalFormStruct::M() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      bmad_normal_form_struct_get_M_info
  );
}
TaylorStructArray1D BmadNormalFormStruct::A() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      bmad_normal_form_struct_get_A_info
  );
}
TaylorStructArray1D BmadNormalFormStruct::A_inv() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      bmad_normal_form_struct_get_A_inv_info
  );
}
TaylorStructArray1D BmadNormalFormStruct::dhdj() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      bmad_normal_form_struct_get_dhdj_info
  );
}
ComplexTaylorStructArray1D BmadNormalFormStruct::F() const {
  return ProxyHelpers::get_type_array_1d<ComplexTaylorStructArray1D>(
      fortran_ptr_,
      bmad_normal_form_struct_get_F_info
  );
}
ComplexTaylorStructArray1D BmadNormalFormStruct::L() const {
  return ProxyHelpers::get_type_array_1d<ComplexTaylorStructArray1D>(
      fortran_ptr_,
      bmad_normal_form_struct_get_L_info
  );
}
ResonanceHStructAlloc1D BmadNormalFormStruct::h() const {
  return ResonanceHStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bmad_normal_form_struct_reallocate_h,
      bmad_normal_form_struct_get_h_info
  );
}
BunchParamsStructAlloc1D BunchTrackStruct::pt() const {
  return BunchParamsStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bunch_track_struct_reallocate_pt,
      bunch_track_struct_get_pt_info
  );
}
double BunchTrackStruct::ds_save() const {
  double value;
  bunch_track_struct_get_ds_save(fortran_ptr_, &value);
  return value;
}
void BunchTrackStruct::set_ds_save(double value) {
  bunch_track_struct_set_ds_save(fortran_ptr_, value);
}
int BunchTrackStruct::n_pt() const {
  int value;
  bunch_track_struct_get_n_pt(fortran_ptr_, &value);
  return value;
}
void BunchTrackStruct::set_n_pt(int value) { bunch_track_struct_set_n_pt(fortran_ptr_, value); }
std::complex<double> SummationRdtStruct::h11001() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h11001(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h11001(std::complex<double> value) {
  summation_rdt_struct_set_h11001(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h00111() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00111(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h00111(std::complex<double> value) {
  summation_rdt_struct_set_h00111(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h20001() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20001(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h20001(std::complex<double> value) {
  summation_rdt_struct_set_h20001(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h00201() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00201(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h00201(std::complex<double> value) {
  summation_rdt_struct_set_h00201(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h10002() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10002(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h10002(std::complex<double> value) {
  summation_rdt_struct_set_h10002(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h21000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h21000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h21000(std::complex<double> value) {
  summation_rdt_struct_set_h21000(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h30000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h30000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h30000(std::complex<double> value) {
  summation_rdt_struct_set_h30000(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h10110() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10110(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h10110(std::complex<double> value) {
  summation_rdt_struct_set_h10110(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h10020() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10020(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h10020(std::complex<double> value) {
  summation_rdt_struct_set_h10020(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h10200() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h10200(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h10200(std::complex<double> value) {
  summation_rdt_struct_set_h10200(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h31000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h31000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h31000(std::complex<double> value) {
  summation_rdt_struct_set_h31000(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h40000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h40000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h40000(std::complex<double> value) {
  summation_rdt_struct_set_h40000(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h20110() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20110(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h20110(std::complex<double> value) {
  summation_rdt_struct_set_h20110(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h11200() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h11200(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h11200(std::complex<double> value) {
  summation_rdt_struct_set_h11200(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h20020() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20020(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h20020(std::complex<double> value) {
  summation_rdt_struct_set_h20020(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h20200() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h20200(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h20200(std::complex<double> value) {
  summation_rdt_struct_set_h20200(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h00310() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00310(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h00310(std::complex<double> value) {
  summation_rdt_struct_set_h00310(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h00400() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00400(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h00400(std::complex<double> value) {
  summation_rdt_struct_set_h00400(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h22000() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h22000(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h22000(std::complex<double> value) {
  summation_rdt_struct_set_h22000(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h00220() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h00220(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h00220(std::complex<double> value) {
  summation_rdt_struct_set_h00220(fortran_ptr_, value);
}
std::complex<double> SummationRdtStruct::h11110() const {
  std::complex<double> c_value;
  summation_rdt_struct_get_h11110(fortran_ptr_, &c_value);
  return c_value;
}
void SummationRdtStruct::set_h11110(std::complex<double> value) {
  summation_rdt_struct_set_h11110(fortran_ptr_, value);
}
std::string TaoEleShapeStruct::ele_id() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_ele_shape_struct_get_ele_id_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeStruct::set_ele_id(const std::string &value) {
  tao_ele_shape_struct_set_ele_id(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoEleShapeStruct::shape() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_ele_shape_struct_get_shape_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeStruct::set_shape(const std::string &value) {
  tao_ele_shape_struct_set_shape(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoEleShapeStruct::color() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_ele_shape_struct_get_color_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeStruct::set_color(const std::string &value) {
  tao_ele_shape_struct_set_color(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoEleShapeStruct::size() const {
  double value;
  tao_ele_shape_struct_get_size(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeStruct::set_size(double value) {
  tao_ele_shape_struct_set_size(fortran_ptr_, value);
}
std::string TaoEleShapeStruct::label() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_ele_shape_struct_get_label_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeStruct::set_label(const std::string &value) {
  tao_ele_shape_struct_set_label(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoEleShapeStruct::draw() const {
  bool value;
  tao_ele_shape_struct_get_draw(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeStruct::set_draw(bool value) { tao_ele_shape_struct_set_draw(fortran_ptr_, value); }
bool TaoEleShapeStruct::multi() const {
  bool value;
  tao_ele_shape_struct_get_multi(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeStruct::set_multi(bool value) {
  tao_ele_shape_struct_set_multi(fortran_ptr_, value);
}
int TaoEleShapeStruct::line_width() const {
  int value;
  tao_ele_shape_struct_get_line_width(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeStruct::set_line_width(int value) {
  tao_ele_shape_struct_set_line_width(fortran_ptr_, value);
}
double TaoEleShapeStruct::offset() const {
  double value;
  tao_ele_shape_struct_get_offset(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeStruct::set_offset(double value) {
  tao_ele_shape_struct_set_offset(fortran_ptr_, value);
}
int TaoEleShapeStruct::ix_key() const {
  int value;
  tao_ele_shape_struct_get_ix_key(fortran_ptr_, &value);
  return value;
}
void TaoEleShapeStruct::set_ix_key(int value) {
  tao_ele_shape_struct_set_ix_key(fortran_ptr_, value);
}
std::string TaoEleShapeStruct::name_ele() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_ele_shape_struct_get_name_ele_info);
  return std::string(arr.data(), arr.size());
}
void TaoEleShapeStruct::set_name_ele(const std::string &value) {
  tao_ele_shape_struct_set_name_ele(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoElePointerStructAlloc1D TaoEleShapeStruct::uni() const {
  return TaoElePointerStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_ele_shape_struct_reallocate_uni,
      tao_ele_shape_struct_get_uni_info
  );
}
ElePointerStructAlloc1D TaoElePointerStruct::eles() const {
  return ElePointerStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_ele_pointer_struct_reallocate_eles,
      tao_ele_pointer_struct_get_eles_info
  );
}
int TaoElePointerStruct::n_loc() const {
  int value;
  tao_ele_pointer_struct_get_n_loc(fortran_ptr_, &value);
  return value;
}
void TaoElePointerStruct::set_n_loc(int value) {
  tao_ele_pointer_struct_set_n_loc(fortran_ptr_, value);
}
std::string TaoCurveStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_name(const std::string &value) {
  tao_curve_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::data_source() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_data_source_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_data_source(const std::string &value) {
  tao_curve_struct_set_data_source(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::data_index() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_data_index_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_data_index(const std::string &value) {
  tao_curve_struct_set_data_index(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::data_type_x() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_data_type_x_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_data_type_x(const std::string &value) {
  tao_curve_struct_set_data_type_x(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::data_type() const {
  return ProxyHelpers::get_string(fortran_ptr_, tao_curve_struct_get_data_type_info);
}
void TaoCurveStruct::set_data_type(const std::string &value) {
  tao_curve_struct_set_data_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::ele_ref_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_ele_ref_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_ele_ref_name(const std::string &value) {
  tao_curve_struct_set_ele_ref_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::legend_text() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_legend_text_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_legend_text(const std::string &value) {
  tao_curve_struct_set_legend_text(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::message_text() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_message_text_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_message_text(const std::string &value) {
  tao_curve_struct_set_message_text(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::component() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_component_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_component(const std::string &value) {
  tao_curve_struct_set_component(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoCurveStruct::why_invalid() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_struct_get_why_invalid_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveStruct::set_why_invalid(const std::string &value) {
  tao_curve_struct_set_why_invalid(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::optional<TaoGraphStruct> TaoCurveStruct::g() const {
  void *ptr;
  tao_curve_struct_get_g(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoGraphStruct(ptr);
}
void TaoCurveStruct::set_g(const TaoGraphStruct &src) {
  tao_curve_struct_set_g(fortran_ptr_, src.get_fortran_ptr());
}
TaoHistogramStruct TaoCurveStruct::hist() const {
  void *ptr;
  tao_curve_struct_get_hist(fortran_ptr_, &ptr);
  return TaoHistogramStruct(ptr);
}
void TaoCurveStruct::set_hist(const TaoHistogramStruct &src) {
  tao_curve_struct_set_hist(fortran_ptr_, src.get_fortran_ptr());
}
TaoCurveColorStruct TaoCurveStruct::z_color() const {
  void *ptr;
  tao_curve_struct_get_z_color(fortran_ptr_, &ptr);
  return TaoCurveColorStruct(ptr);
}
void TaoCurveStruct::set_z_color(const TaoCurveColorStruct &src) {
  tao_curve_struct_set_z_color(fortran_ptr_, src.get_fortran_ptr());
}
RealAlloc1D TaoCurveStruct::x_line() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_x_line,
      tao_curve_struct_get_x_line_info
  );
}
void TaoCurveStruct::set_x_line(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_x_line(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::y_line() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_y_line,
      tao_curve_struct_get_y_line_info
  );
}
void TaoCurveStruct::set_y_line(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_y_line(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::y2_line() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_y2_line,
      tao_curve_struct_get_y2_line_info
  );
}
void TaoCurveStruct::set_y2_line(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_y2_line(fortran_ptr_, v.data(), shape);
}
IntAlloc1D TaoCurveStruct::ix_line() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_ix_line,
      tao_curve_struct_get_ix_line_info
  );
}
void TaoCurveStruct::set_ix_line(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_ix_line(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::x_symb() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_x_symb,
      tao_curve_struct_get_x_symb_info
  );
}
void TaoCurveStruct::set_x_symb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_x_symb(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::y_symb() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_y_symb,
      tao_curve_struct_get_y_symb_info
  );
}
void TaoCurveStruct::set_y_symb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_y_symb(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::z_symb() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_z_symb,
      tao_curve_struct_get_z_symb_info
  );
}
void TaoCurveStruct::set_z_symb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_z_symb(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::err_symb() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_err_symb,
      tao_curve_struct_get_err_symb_info
  );
}
void TaoCurveStruct::set_err_symb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_err_symb(fortran_ptr_, v.data(), shape);
}
RealAlloc1D TaoCurveStruct::symb_size() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_symb_size,
      tao_curve_struct_get_symb_size_info
  );
}
void TaoCurveStruct::set_symb_size(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_symb_size(fortran_ptr_, v.data(), shape);
}
IntAlloc1D TaoCurveStruct::ix_symb() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_curve_struct_reallocate_ix_symb,
      tao_curve_struct_get_ix_symb_info
  );
}
void TaoCurveStruct::set_ix_symb(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_curve_struct_set_ix_symb(fortran_ptr_, v.data(), shape);
}
double TaoCurveStruct::y_axis_scale_factor() const {
  double value;
  tao_curve_struct_get_y_axis_scale_factor(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_y_axis_scale_factor(double value) {
  tao_curve_struct_set_y_axis_scale_factor(fortran_ptr_, value);
}
QpLineStruct TaoCurveStruct::line() const {
  void *ptr;
  tao_curve_struct_get_line(fortran_ptr_, &ptr);
  return QpLineStruct(ptr);
}
void TaoCurveStruct::set_line(const QpLineStruct &src) {
  tao_curve_struct_set_line(fortran_ptr_, src.get_fortran_ptr());
}
QpSymbolStruct TaoCurveStruct::symbol() const {
  void *ptr;
  tao_curve_struct_get_symbol(fortran_ptr_, &ptr);
  return QpSymbolStruct(ptr);
}
void TaoCurveStruct::set_symbol(const QpSymbolStruct &src) {
  tao_curve_struct_set_symbol(fortran_ptr_, src.get_fortran_ptr());
}
TaoCurveOrbitStruct TaoCurveStruct::orbit() const {
  void *ptr;
  tao_curve_struct_get_orbit(fortran_ptr_, &ptr);
  return TaoCurveOrbitStruct(ptr);
}
void TaoCurveStruct::set_orbit(const TaoCurveOrbitStruct &src) {
  tao_curve_struct_set_orbit(fortran_ptr_, src.get_fortran_ptr());
}
int TaoCurveStruct::ix_universe() const {
  int value;
  tao_curve_struct_get_ix_universe(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_ix_universe(int value) {
  tao_curve_struct_set_ix_universe(fortran_ptr_, value);
}
int TaoCurveStruct::symbol_every() const {
  int value;
  tao_curve_struct_get_symbol_every(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_symbol_every(int value) {
  tao_curve_struct_set_symbol_every(fortran_ptr_, value);
}
int TaoCurveStruct::ix_branch() const {
  int value;
  tao_curve_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_ix_branch(int value) {
  tao_curve_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoCurveStruct::ix_bunch() const {
  int value;
  tao_curve_struct_get_ix_bunch(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_ix_bunch(int value) { tao_curve_struct_set_ix_bunch(fortran_ptr_, value); }
int TaoCurveStruct::n_turn() const {
  int value;
  tao_curve_struct_get_n_turn(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_n_turn(int value) { tao_curve_struct_set_n_turn(fortran_ptr_, value); }
bool TaoCurveStruct::use_y2() const {
  bool value;
  tao_curve_struct_get_use_y2(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_use_y2(bool value) { tao_curve_struct_set_use_y2(fortran_ptr_, value); }
bool TaoCurveStruct::draw_line() const {
  bool value;
  tao_curve_struct_get_draw_line(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_draw_line(bool value) {
  tao_curve_struct_set_draw_line(fortran_ptr_, value);
}
bool TaoCurveStruct::draw_symbols() const {
  bool value;
  tao_curve_struct_get_draw_symbols(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_draw_symbols(bool value) {
  tao_curve_struct_set_draw_symbols(fortran_ptr_, value);
}
bool TaoCurveStruct::draw_symbol_index() const {
  bool value;
  tao_curve_struct_get_draw_symbol_index(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_draw_symbol_index(bool value) {
  tao_curve_struct_set_draw_symbol_index(fortran_ptr_, value);
}
bool TaoCurveStruct::draw_error_bars() const {
  bool value;
  tao_curve_struct_get_draw_error_bars(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_draw_error_bars(bool value) {
  tao_curve_struct_set_draw_error_bars(fortran_ptr_, value);
}
bool TaoCurveStruct::smooth_line_calc() const {
  bool value;
  tao_curve_struct_get_smooth_line_calc(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_smooth_line_calc(bool value) {
  tao_curve_struct_set_smooth_line_calc(fortran_ptr_, value);
}
bool TaoCurveStruct::valid() const {
  bool value;
  tao_curve_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoCurveStruct::set_valid(bool value) { tao_curve_struct_set_valid(fortran_ptr_, value); }
std::string TaoCurveColorStruct::data_type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_curve_color_struct_get_data_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoCurveColorStruct::set_data_type(const std::string &value) {
  tao_curve_color_struct_set_data_type(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
bool TaoCurveColorStruct::is_on() const {
  bool value;
  tao_curve_color_struct_get_is_on(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorStruct::set_is_on(bool value) {
  tao_curve_color_struct_set_is_on(fortran_ptr_, value);
}
double TaoCurveColorStruct::min() const {
  double value;
  tao_curve_color_struct_get_min(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorStruct::set_min(double value) {
  tao_curve_color_struct_set_min(fortran_ptr_, value);
}
double TaoCurveColorStruct::max() const {
  double value;
  tao_curve_color_struct_get_max(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorStruct::set_max(double value) {
  tao_curve_color_struct_set_max(fortran_ptr_, value);
}
bool TaoCurveColorStruct::autoscale() const {
  bool value;
  tao_curve_color_struct_get_autoscale(fortran_ptr_, &value);
  return value;
}
void TaoCurveColorStruct::set_autoscale(bool value) {
  tao_curve_color_struct_set_autoscale(fortran_ptr_, value);
}
double TaoCurveOrbitStruct::x() const {
  double value;
  tao_curve_orbit_struct_get_x(fortran_ptr_, &value);
  return value;
}
void TaoCurveOrbitStruct::set_x(double value) { tao_curve_orbit_struct_set_x(fortran_ptr_, value); }
double TaoCurveOrbitStruct::y() const {
  double value;
  tao_curve_orbit_struct_get_y(fortran_ptr_, &value);
  return value;
}
void TaoCurveOrbitStruct::set_y(double value) { tao_curve_orbit_struct_set_y(fortran_ptr_, value); }
double TaoCurveOrbitStruct::t() const {
  double value;
  tao_curve_orbit_struct_get_t(fortran_ptr_, &value);
  return value;
}
void TaoCurveOrbitStruct::set_t(double value) { tao_curve_orbit_struct_set_t(fortran_ptr_, value); }
bool TaoHistogramStruct::density_normalized() const {
  bool value;
  tao_histogram_struct_get_density_normalized(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_density_normalized(bool value) {
  tao_histogram_struct_set_density_normalized(fortran_ptr_, value);
}
bool TaoHistogramStruct::weight_by_charge() const {
  bool value;
  tao_histogram_struct_get_weight_by_charge(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_weight_by_charge(bool value) {
  tao_histogram_struct_set_weight_by_charge(fortran_ptr_, value);
}
double TaoHistogramStruct::minimum() const {
  double value;
  tao_histogram_struct_get_minimum(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_minimum(double value) {
  tao_histogram_struct_set_minimum(fortran_ptr_, value);
}
double TaoHistogramStruct::maximum() const {
  double value;
  tao_histogram_struct_get_maximum(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_maximum(double value) {
  tao_histogram_struct_set_maximum(fortran_ptr_, value);
}
double TaoHistogramStruct::width() const {
  double value;
  tao_histogram_struct_get_width(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_width(double value) {
  tao_histogram_struct_set_width(fortran_ptr_, value);
}
double TaoHistogramStruct::center() const {
  double value;
  tao_histogram_struct_get_center(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_center(double value) {
  tao_histogram_struct_set_center(fortran_ptr_, value);
}
int TaoHistogramStruct::number() const {
  int value;
  tao_histogram_struct_get_number(fortran_ptr_, &value);
  return value;
}
void TaoHistogramStruct::set_number(int value) {
  tao_histogram_struct_set_number(fortran_ptr_, value);
}
int LatEleOrder1Struct::ix_branch() const {
  int value;
  lat_ele_order1_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void LatEleOrder1Struct::set_ix_branch(int value) {
  lat_ele_order1_struct_set_ix_branch(fortran_ptr_, value);
}
int LatEleOrder1Struct::ix_order() const {
  int value;
  lat_ele_order1_struct_get_ix_order(fortran_ptr_, &value);
  return value;
}
void LatEleOrder1Struct::set_ix_order(int value) {
  lat_ele_order1_struct_set_ix_order(fortran_ptr_, value);
}
LatEleOrder1StructAlloc1D LatEleOrderArrayStruct::ele() const {
  return LatEleOrder1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_ele_order_array_struct_reallocate_ele,
      lat_ele_order_array_struct_get_ele_info
  );
}
FArray2D<double> TaoLatSigmaStruct::mat() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_lat_sigma_struct_get_mat_info);
}
void TaoLatSigmaStruct::set_mat(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_lat_sigma_struct_set_mat, v);
}
TaoSpinDnDpzStruct TaoSpinEleStruct::dn_dpz() const {
  void *ptr;
  tao_spin_ele_struct_get_dn_dpz(fortran_ptr_, &ptr);
  return TaoSpinDnDpzStruct(ptr);
}
void TaoSpinEleStruct::set_dn_dpz(const TaoSpinDnDpzStruct &src) {
  tao_spin_ele_struct_set_dn_dpz(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TaoSpinEleStruct::orb_eigen_val() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      tao_spin_ele_struct_get_orb_eigen_val_info
  );
}
void TaoSpinEleStruct::set_orb_eigen_val(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_spin_ele_struct_set_orb_eigen_val(fortran_ptr_, v.data(), shape);
}
FArray2D<double> TaoSpinEleStruct::orb_eigen_vec() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      tao_spin_ele_struct_get_orb_eigen_vec_info
  );
}
void TaoSpinEleStruct::set_orb_eigen_vec(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_spin_ele_struct_set_orb_eigen_vec, v);
}
FArray2D<double> TaoSpinEleStruct::spin_eigen_vec() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      tao_spin_ele_struct_get_spin_eigen_vec_info
  );
}
void TaoSpinEleStruct::set_spin_eigen_vec(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_spin_ele_struct_set_spin_eigen_vec, v);
}
bool TaoSpinEleStruct::valid() const {
  bool value;
  tao_spin_ele_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoSpinEleStruct::set_valid(bool value) { tao_spin_ele_struct_set_valid(fortran_ptr_, value); }
EleStruct TaoPlotCacheStruct::ele_to_s() const {
  void *ptr;
  tao_plot_cache_struct_get_ele_to_s(fortran_ptr_, &ptr);
  return EleStruct(ptr);
}
void TaoPlotCacheStruct::set_ele_to_s(const EleStruct &src) {
  tao_plot_cache_struct_set_ele_to_s(fortran_ptr_, src.get_fortran_ptr());
}
CoordStruct TaoPlotCacheStruct::orbit() const {
  void *ptr;
  tao_plot_cache_struct_get_orbit(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void TaoPlotCacheStruct::set_orbit(const CoordStruct &src) {
  tao_plot_cache_struct_set_orbit(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoPlotCacheStruct::err() const {
  bool value;
  tao_plot_cache_struct_get_err(fortran_ptr_, &value);
  return value;
}
void TaoPlotCacheStruct::set_err(bool value) { tao_plot_cache_struct_set_err(fortran_ptr_, value); }
double TaoSpinPolarizationStruct::tune() const {
  double value;
  tao_spin_polarization_struct_get_tune(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_tune(double value) {
  tao_spin_polarization_struct_set_tune(fortran_ptr_, value);
}
double TaoSpinPolarizationStruct::pol_limit_st() const {
  double value;
  tao_spin_polarization_struct_get_pol_limit_st(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_pol_limit_st(double value) {
  tao_spin_polarization_struct_set_pol_limit_st(fortran_ptr_, value);
}
double TaoSpinPolarizationStruct::pol_limit_dk() const {
  double value;
  tao_spin_polarization_struct_get_pol_limit_dk(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_pol_limit_dk(double value) {
  tao_spin_polarization_struct_set_pol_limit_dk(fortran_ptr_, value);
}
FArray1D<double> TaoSpinPolarizationStruct::pol_limit_dk_partial() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      tao_spin_polarization_struct_get_pol_limit_dk_partial_info
  );
}
void TaoSpinPolarizationStruct::set_pol_limit_dk_partial(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_spin_polarization_struct_set_pol_limit_dk_partial(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoSpinPolarizationStruct::pol_limit_dk_partial2() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      tao_spin_polarization_struct_get_pol_limit_dk_partial2_info
  );
}
void TaoSpinPolarizationStruct::set_pol_limit_dk_partial2(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_spin_polarization_struct_set_pol_limit_dk_partial2(fortran_ptr_, v.data(), shape);
}
double TaoSpinPolarizationStruct::pol_rate_bks() const {
  double value;
  tao_spin_polarization_struct_get_pol_rate_bks(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_pol_rate_bks(double value) {
  tao_spin_polarization_struct_set_pol_rate_bks(fortran_ptr_, value);
}
double TaoSpinPolarizationStruct::depol_rate() const {
  double value;
  tao_spin_polarization_struct_get_depol_rate(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_depol_rate(double value) {
  tao_spin_polarization_struct_set_depol_rate(fortran_ptr_, value);
}
FArray1D<double> TaoSpinPolarizationStruct::depol_rate_partial() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      tao_spin_polarization_struct_get_depol_rate_partial_info
  );
}
void TaoSpinPolarizationStruct::set_depol_rate_partial(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_spin_polarization_struct_set_depol_rate_partial(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoSpinPolarizationStruct::depol_rate_partial2() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      tao_spin_polarization_struct_get_depol_rate_partial2_info
  );
}
void TaoSpinPolarizationStruct::set_depol_rate_partial2(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_spin_polarization_struct_set_depol_rate_partial2(fortran_ptr_, v.data(), shape);
}
double TaoSpinPolarizationStruct::integral_bn() const {
  double value;
  tao_spin_polarization_struct_get_integral_bn(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_integral_bn(double value) {
  tao_spin_polarization_struct_set_integral_bn(fortran_ptr_, value);
}
double TaoSpinPolarizationStruct::integral_bdn() const {
  double value;
  tao_spin_polarization_struct_get_integral_bdn(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_integral_bdn(double value) {
  tao_spin_polarization_struct_set_integral_bdn(fortran_ptr_, value);
}
double TaoSpinPolarizationStruct::integral_1ns() const {
  double value;
  tao_spin_polarization_struct_get_integral_1ns(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_integral_1ns(double value) {
  tao_spin_polarization_struct_set_integral_1ns(fortran_ptr_, value);
}
double TaoSpinPolarizationStruct::integral_dn2() const {
  double value;
  tao_spin_polarization_struct_get_integral_dn2(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_integral_dn2(double value) {
  tao_spin_polarization_struct_set_integral_dn2(fortran_ptr_, value);
}
bool TaoSpinPolarizationStruct::valid() const {
  bool value;
  tao_spin_polarization_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoSpinPolarizationStruct::set_valid(bool value) {
  tao_spin_polarization_struct_set_valid(fortran_ptr_, value);
}
SpinOrbitMap1Struct TaoSpinPolarizationStruct::q_1turn() const {
  void *ptr;
  tao_spin_polarization_struct_get_q_1turn(fortran_ptr_, &ptr);
  return SpinOrbitMap1Struct(ptr);
}
void TaoSpinPolarizationStruct::set_q_1turn(const SpinOrbitMap1Struct &src) {
  tao_spin_polarization_struct_set_q_1turn(fortran_ptr_, src.get_fortran_ptr());
}
SpinOrbitMap1StructAlloc1D TaoSpinPolarizationStruct::q_ele() const {
  return SpinOrbitMap1StructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_spin_polarization_struct_reallocate_q_ele,
      tao_spin_polarization_struct_get_q_ele_info
  );
}
std::optional<TaoLatticeStruct> TaoLatticeBranchStruct::tao_lat() const {
  void *ptr;
  tao_lattice_branch_struct_get_tao_lat(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeStruct(ptr);
}
void TaoLatticeBranchStruct::set_tao_lat(const TaoLatticeStruct &src) {
  tao_lattice_branch_struct_set_tao_lat(fortran_ptr_, src.get_fortran_ptr());
}
TaoLatSigmaStructAlloc1D TaoLatticeBranchStruct::lat_sigma() const {
  return TaoLatSigmaStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_lat_sigma,
      tao_lattice_branch_struct_get_lat_sigma_info
  );
}
TaoSpinEleStructAlloc1D TaoLatticeBranchStruct::spin_ele() const {
  return TaoSpinEleStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_spin_ele,
      tao_lattice_branch_struct_get_spin_ele_info
  );
}
BunchParamsStructAlloc1D TaoLatticeBranchStruct::bunch_params() const {
  return BunchParamsStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_bunch_params,
      tao_lattice_branch_struct_get_bunch_params_info
  );
}
BunchTrackStructAlloc1D TaoLatticeBranchStruct::bunch_params_comb() const {
  return BunchTrackStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_bunch_params_comb,
      tao_lattice_branch_struct_get_bunch_params_comb_info
  );
}
CoordStructAlloc1D TaoLatticeBranchStruct::orbit() const {
  return CoordStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_orbit,
      tao_lattice_branch_struct_get_orbit_info
  );
}
TaoPlotCacheStructAlloc1D TaoLatticeBranchStruct::plot_cache() const {
  return TaoPlotCacheStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_plot_cache,
      tao_lattice_branch_struct_get_plot_cache_info
  );
}
TaoSpinPolarizationStruct TaoLatticeBranchStruct::spin() const {
  void *ptr;
  tao_lattice_branch_struct_get_spin(fortran_ptr_, &ptr);
  return TaoSpinPolarizationStruct(ptr);
}
void TaoLatticeBranchStruct::set_spin(const TaoSpinPolarizationStruct &src) {
  tao_lattice_branch_struct_set_spin(fortran_ptr_, src.get_fortran_ptr());
}
SummationRdtStruct TaoLatticeBranchStruct::srdt() const {
  void *ptr;
  tao_lattice_branch_struct_get_srdt(fortran_ptr_, &ptr);
  return SummationRdtStruct(ptr);
}
void TaoLatticeBranchStruct::set_srdt(const SummationRdtStruct &src) {
  tao_lattice_branch_struct_set_srdt(fortran_ptr_, src.get_fortran_ptr());
}
CoordStruct TaoLatticeBranchStruct::orb0() const {
  void *ptr;
  tao_lattice_branch_struct_get_orb0(fortran_ptr_, &ptr);
  return CoordStruct(ptr);
}
void TaoLatticeBranchStruct::set_orb0(const CoordStruct &src) {
  tao_lattice_branch_struct_set_orb0(fortran_ptr_, src.get_fortran_ptr());
}
NormalModesStruct TaoLatticeBranchStruct::modes_ri() const {
  void *ptr;
  tao_lattice_branch_struct_get_modes_ri(fortran_ptr_, &ptr);
  return NormalModesStruct(ptr);
}
void TaoLatticeBranchStruct::set_modes_ri(const NormalModesStruct &src) {
  tao_lattice_branch_struct_set_modes_ri(fortran_ptr_, src.get_fortran_ptr());
}
NormalModesStruct TaoLatticeBranchStruct::modes_6d() const {
  void *ptr;
  tao_lattice_branch_struct_get_modes_6d(fortran_ptr_, &ptr);
  return NormalModesStruct(ptr);
}
void TaoLatticeBranchStruct::set_modes_6d(const NormalModesStruct &src) {
  tao_lattice_branch_struct_set_modes_6d(fortran_ptr_, src.get_fortran_ptr());
}
PtcNormalFormStruct TaoLatticeBranchStruct::ptc_normal_form() const {
  void *ptr;
  tao_lattice_branch_struct_get_ptc_normal_form(fortran_ptr_, &ptr);
  return PtcNormalFormStruct(ptr);
}
void TaoLatticeBranchStruct::set_ptc_normal_form(const PtcNormalFormStruct &src) {
  tao_lattice_branch_struct_set_ptc_normal_form(fortran_ptr_, src.get_fortran_ptr());
}
BmadNormalFormStruct TaoLatticeBranchStruct::bmad_normal_form() const {
  void *ptr;
  tao_lattice_branch_struct_get_bmad_normal_form(fortran_ptr_, &ptr);
  return BmadNormalFormStruct(ptr);
}
void TaoLatticeBranchStruct::set_bmad_normal_form(const BmadNormalFormStruct &src) {
  tao_lattice_branch_struct_set_bmad_normal_form(fortran_ptr_, src.get_fortran_ptr());
}
CoordStructAlloc1D TaoLatticeBranchStruct::high_E_orb() const {
  return CoordStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_high_E_orb,
      tao_lattice_branch_struct_get_high_E_orb_info
  );
}
CoordStructAlloc1D TaoLatticeBranchStruct::low_E_orb() const {
  return CoordStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_branch_struct_reallocate_low_E_orb,
      tao_lattice_branch_struct_get_low_E_orb_info
  );
}
TaylorStructArray1D TaoLatticeBranchStruct::taylor_save() const {
  return ProxyHelpers::get_type_array_1d<TaylorStructArray1D>(
      fortran_ptr_,
      tao_lattice_branch_struct_get_taylor_save_info
  );
}
double TaoLatticeBranchStruct::cache_x_min() const {
  double value;
  tao_lattice_branch_struct_get_cache_x_min(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_cache_x_min(double value) {
  tao_lattice_branch_struct_set_cache_x_min(fortran_ptr_, value);
}
double TaoLatticeBranchStruct::cache_x_max() const {
  double value;
  tao_lattice_branch_struct_get_cache_x_max(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_cache_x_max(double value) {
  tao_lattice_branch_struct_set_cache_x_max(fortran_ptr_, value);
}
double TaoLatticeBranchStruct::comb_ds_save() const {
  double value;
  tao_lattice_branch_struct_get_comb_ds_save(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_comb_ds_save(double value) {
  tao_lattice_branch_struct_set_comb_ds_save(fortran_ptr_, value);
}
int TaoLatticeBranchStruct::ix_ref_taylor() const {
  int value;
  tao_lattice_branch_struct_get_ix_ref_taylor(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_ix_ref_taylor(int value) {
  tao_lattice_branch_struct_set_ix_ref_taylor(fortran_ptr_, value);
}
int TaoLatticeBranchStruct::ix_ele_taylor() const {
  int value;
  tao_lattice_branch_struct_get_ix_ele_taylor(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_ix_ele_taylor(int value) {
  tao_lattice_branch_struct_set_ix_ele_taylor(fortran_ptr_, value);
}
int TaoLatticeBranchStruct::track_state() const {
  int value;
  tao_lattice_branch_struct_get_track_state(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_track_state(int value) {
  tao_lattice_branch_struct_set_track_state(fortran_ptr_, value);
}
int TaoLatticeBranchStruct::cache_n_pts() const {
  int value;
  tao_lattice_branch_struct_get_cache_n_pts(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_cache_n_pts(int value) {
  tao_lattice_branch_struct_set_cache_n_pts(fortran_ptr_, value);
}
int TaoLatticeBranchStruct::ix_rad_int_cache() const {
  int value;
  tao_lattice_branch_struct_get_ix_rad_int_cache(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_ix_rad_int_cache(int value) {
  tao_lattice_branch_struct_set_ix_rad_int_cache(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::has_open_match_element() const {
  bool value;
  tao_lattice_branch_struct_get_has_open_match_element(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_has_open_match_element(bool value) {
  tao_lattice_branch_struct_set_has_open_match_element(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::plot_cache_valid() const {
  bool value;
  tao_lattice_branch_struct_get_plot_cache_valid(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_plot_cache_valid(bool value) {
  tao_lattice_branch_struct_set_plot_cache_valid(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::spin_map_valid() const {
  bool value;
  tao_lattice_branch_struct_get_spin_map_valid(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_spin_map_valid(bool value) {
  tao_lattice_branch_struct_set_spin_map_valid(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::twiss_valid() const {
  bool value;
  tao_lattice_branch_struct_get_twiss_valid(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_twiss_valid(bool value) {
  tao_lattice_branch_struct_set_twiss_valid(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::mode_flip_here() const {
  bool value;
  tao_lattice_branch_struct_get_mode_flip_here(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_mode_flip_here(bool value) {
  tao_lattice_branch_struct_set_mode_flip_here(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::chrom_calc_ok() const {
  bool value;
  tao_lattice_branch_struct_get_chrom_calc_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_chrom_calc_ok(bool value) {
  tao_lattice_branch_struct_set_chrom_calc_ok(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::rad_int_calc_ok() const {
  bool value;
  tao_lattice_branch_struct_get_rad_int_calc_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_rad_int_calc_ok(bool value) {
  tao_lattice_branch_struct_set_rad_int_calc_ok(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::emit_6d_calc_ok() const {
  bool value;
  tao_lattice_branch_struct_get_emit_6d_calc_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_emit_6d_calc_ok(bool value) {
  tao_lattice_branch_struct_set_emit_6d_calc_ok(fortran_ptr_, value);
}
bool TaoLatticeBranchStruct::sigma_track_ok() const {
  bool value;
  tao_lattice_branch_struct_get_sigma_track_ok(fortran_ptr_, &value);
  return value;
}
void TaoLatticeBranchStruct::set_sigma_track_ok(bool value) {
  tao_lattice_branch_struct_set_sigma_track_ok(fortran_ptr_, value);
}
BeamStruct TaoModelElementStruct::beam() const {
  void *ptr;
  tao_model_element_struct_get_beam(fortran_ptr_, &ptr);
  return BeamStruct(ptr);
}
void TaoModelElementStruct::set_beam(const BeamStruct &src) {
  tao_model_element_struct_set_beam(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoModelElementStruct::save_beam_internally() const {
  bool value;
  tao_model_element_struct_get_save_beam_internally(fortran_ptr_, &value);
  return value;
}
void TaoModelElementStruct::set_save_beam_internally(bool value) {
  tao_model_element_struct_set_save_beam_internally(fortran_ptr_, value);
}
bool TaoModelElementStruct::save_beam_to_file() const {
  bool value;
  tao_model_element_struct_get_save_beam_to_file(fortran_ptr_, &value);
  return value;
}
void TaoModelElementStruct::set_save_beam_to_file(bool value) {
  tao_model_element_struct_set_save_beam_to_file(fortran_ptr_, value);
}
BeamStruct TaoBeamBranchStruct::beam_at_start() const {
  void *ptr;
  tao_beam_branch_struct_get_beam_at_start(fortran_ptr_, &ptr);
  return BeamStruct(ptr);
}
void TaoBeamBranchStruct::set_beam_at_start(const BeamStruct &src) {
  tao_beam_branch_struct_set_beam_at_start(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitStruct TaoBeamBranchStruct::beam_init() const {
  void *ptr;
  tao_beam_branch_struct_get_beam_init(fortran_ptr_, &ptr);
  return BeamInitStruct(ptr);
}
void TaoBeamBranchStruct::set_beam_init(const BeamInitStruct &src) {
  tao_beam_branch_struct_set_beam_init(fortran_ptr_, src.get_fortran_ptr());
}
BeamInitStruct TaoBeamBranchStruct::beam_init_used() const {
  void *ptr;
  tao_beam_branch_struct_get_beam_init_used(fortran_ptr_, &ptr);
  return BeamInitStruct(ptr);
}
void TaoBeamBranchStruct::set_beam_init_used(const BeamInitStruct &src) {
  tao_beam_branch_struct_set_beam_init_used(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoBeamBranchStruct::init_starting_distribution() const {
  bool value;
  tao_beam_branch_struct_get_init_starting_distribution(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchStruct::set_init_starting_distribution(bool value) {
  tao_beam_branch_struct_set_init_starting_distribution(fortran_ptr_, value);
}
std::string TaoBeamBranchStruct::track_start() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_beam_branch_struct_get_track_start_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamBranchStruct::set_track_start(const std::string &value) {
  tao_beam_branch_struct_set_track_start(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoBeamBranchStruct::track_end() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_beam_branch_struct_get_track_end_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamBranchStruct::set_track_end(const std::string &value) {
  tao_beam_branch_struct_set_track_end(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
int TaoBeamBranchStruct::ix_branch() const {
  int value;
  tao_beam_branch_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchStruct::set_ix_branch(int value) {
  tao_beam_branch_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoBeamBranchStruct::ix_track_start() const {
  int value;
  tao_beam_branch_struct_get_ix_track_start(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchStruct::set_ix_track_start(int value) {
  tao_beam_branch_struct_set_ix_track_start(fortran_ptr_, value);
}
int TaoBeamBranchStruct::ix_track_end() const {
  int value;
  tao_beam_branch_struct_get_ix_track_end(fortran_ptr_, &value);
  return value;
}
void TaoBeamBranchStruct::set_ix_track_end(int value) {
  tao_beam_branch_struct_set_ix_track_end(fortran_ptr_, value);
}
std::string TaoD1DataStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_d1_data_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD1DataStruct::set_name(const std::string &value) {
  tao_d1_data_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::optional<TaoD2DataStruct> TaoD1DataStruct::d2() const {
  void *ptr;
  tao_d1_data_struct_get_d2(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoD2DataStruct(ptr);
}
void TaoD1DataStruct::set_d2(const TaoD2DataStruct &src) {
  tao_d1_data_struct_set_d2(fortran_ptr_, src.get_fortran_ptr());
}
TaoDataStructArray1D TaoD1DataStruct::d() const {
  return ProxyHelpers::get_type_array_1d<TaoDataStructArray1D>(
      fortran_ptr_,
      tao_d1_data_struct_get_d_info
  );
}
std::string TaoD2DataStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_d2_data_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataStruct::set_name(const std::string &value) {
  tao_d2_data_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoD2DataStruct::data_file_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_d2_data_struct_get_data_file_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataStruct::set_data_file_name(const std::string &value) {
  tao_d2_data_struct_set_data_file_name(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoD2DataStruct::ref_file_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_d2_data_struct_get_ref_file_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataStruct::set_ref_file_name(const std::string &value) {
  tao_d2_data_struct_set_ref_file_name(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoD2DataStruct::data_date() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_d2_data_struct_get_data_date_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataStruct::set_data_date(const std::string &value) {
  tao_d2_data_struct_set_data_date(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoD2DataStruct::ref_date() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_d2_data_struct_get_ref_date_info);
  return std::string(arr.data(), arr.size());
}
void TaoD2DataStruct::set_ref_date(const std::string &value) {
  tao_d2_data_struct_set_ref_date(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D TaoD2DataStruct::descrip() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, tao_d2_data_struct_get_descrip_info);
}
TaoD1DataStructAlloc1D TaoD2DataStruct::d1() const {
  return TaoD1DataStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_d2_data_struct_reallocate_d1,
      tao_d2_data_struct_get_d1_info
  );
}
int TaoD2DataStruct::ix_universe() const {
  int value;
  tao_d2_data_struct_get_ix_universe(fortran_ptr_, &value);
  return value;
}
void TaoD2DataStruct::set_ix_universe(int value) {
  tao_d2_data_struct_set_ix_universe(fortran_ptr_, value);
}
int TaoD2DataStruct::ix_d2_data() const {
  int value;
  tao_d2_data_struct_get_ix_d2_data(fortran_ptr_, &value);
  return value;
}
void TaoD2DataStruct::set_ix_d2_data(int value) {
  tao_d2_data_struct_set_ix_d2_data(fortran_ptr_, value);
}
int TaoD2DataStruct::ix_ref() const {
  int value;
  tao_d2_data_struct_get_ix_ref(fortran_ptr_, &value);
  return value;
}
void TaoD2DataStruct::set_ix_ref(int value) { tao_d2_data_struct_set_ix_ref(fortran_ptr_, value); }
bool TaoD2DataStruct::data_read_in() const {
  bool value;
  tao_d2_data_struct_get_data_read_in(fortran_ptr_, &value);
  return value;
}
void TaoD2DataStruct::set_data_read_in(bool value) {
  tao_d2_data_struct_set_data_read_in(fortran_ptr_, value);
}
bool TaoD2DataStruct::ref_read_in() const {
  bool value;
  tao_d2_data_struct_get_ref_read_in(fortran_ptr_, &value);
  return value;
}
void TaoD2DataStruct::set_ref_read_in(bool value) {
  tao_d2_data_struct_set_ref_read_in(fortran_ptr_, value);
}
std::string TaoDataVarComponentStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_var_component_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataVarComponentStruct::set_name(const std::string &value) {
  tao_data_var_component_struct_set_name(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
double TaoDataVarComponentStruct::sign() const {
  double value;
  tao_data_var_component_struct_get_sign(fortran_ptr_, &value);
  return value;
}
void TaoDataVarComponentStruct::set_sign(double value) {
  tao_data_var_component_struct_set_sign(fortran_ptr_, value);
}
std::string TaoGraphStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_graph_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphStruct::set_name(const std::string &value) {
  tao_graph_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGraphStruct::type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_graph_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphStruct::set_type(const std::string &value) {
  tao_graph_struct_set_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGraphStruct::title() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_graph_struct_get_title_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphStruct::set_title(const std::string &value) {
  tao_graph_struct_set_title(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGraphStruct::title_suffix() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_graph_struct_get_title_suffix_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphStruct::set_title_suffix(const std::string &value) {
  tao_graph_struct_set_title_suffix(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
FCharArray1D TaoGraphStruct::text_legend() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, tao_graph_struct_get_text_legend_info);
}
FCharArray1D TaoGraphStruct::text_legend_out() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, tao_graph_struct_get_text_legend_out_info);
}
std::string TaoGraphStruct::why_invalid() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_graph_struct_get_why_invalid_info);
  return std::string(arr.data(), arr.size());
}
void TaoGraphStruct::set_why_invalid(const std::string &value) {
  tao_graph_struct_set_why_invalid(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoCurveStructAlloc1D TaoGraphStruct::curve() const {
  return TaoCurveStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_graph_struct_reallocate_curve,
      tao_graph_struct_get_curve_info
  );
}
std::optional<TaoPlotStruct> TaoGraphStruct::p() const {
  void *ptr;
  tao_graph_struct_get_p(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoPlotStruct(ptr);
}
void TaoGraphStruct::set_p(const TaoPlotStruct &src) {
  tao_graph_struct_set_p(fortran_ptr_, src.get_fortran_ptr());
}
TaoFloorPlanStruct TaoGraphStruct::floor_plan() const {
  void *ptr;
  tao_graph_struct_get_floor_plan(fortran_ptr_, &ptr);
  return TaoFloorPlanStruct(ptr);
}
void TaoGraphStruct::set_floor_plan(const TaoFloorPlanStruct &src) {
  tao_graph_struct_set_floor_plan(fortran_ptr_, src.get_fortran_ptr());
}
QpPointStruct TaoGraphStruct::text_legend_origin() const {
  void *ptr;
  tao_graph_struct_get_text_legend_origin(fortran_ptr_, &ptr);
  return QpPointStruct(ptr);
}
void TaoGraphStruct::set_text_legend_origin(const QpPointStruct &src) {
  tao_graph_struct_set_text_legend_origin(fortran_ptr_, src.get_fortran_ptr());
}
QpPointStruct TaoGraphStruct::curve_legend_origin() const {
  void *ptr;
  tao_graph_struct_get_curve_legend_origin(fortran_ptr_, &ptr);
  return QpPointStruct(ptr);
}
void TaoGraphStruct::set_curve_legend_origin(const QpPointStruct &src) {
  tao_graph_struct_set_curve_legend_origin(fortran_ptr_, src.get_fortran_ptr());
}
QpLegendStruct TaoGraphStruct::curve_legend() const {
  void *ptr;
  tao_graph_struct_get_curve_legend(fortran_ptr_, &ptr);
  return QpLegendStruct(ptr);
}
void TaoGraphStruct::set_curve_legend(const QpLegendStruct &src) {
  tao_graph_struct_set_curve_legend(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisStruct TaoGraphStruct::x() const {
  void *ptr;
  tao_graph_struct_get_x(fortran_ptr_, &ptr);
  return QpAxisStruct(ptr);
}
void TaoGraphStruct::set_x(const QpAxisStruct &src) {
  tao_graph_struct_set_x(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisStruct TaoGraphStruct::y() const {
  void *ptr;
  tao_graph_struct_get_y(fortran_ptr_, &ptr);
  return QpAxisStruct(ptr);
}
void TaoGraphStruct::set_y(const QpAxisStruct &src) {
  tao_graph_struct_set_y(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisStruct TaoGraphStruct::x2() const {
  void *ptr;
  tao_graph_struct_get_x2(fortran_ptr_, &ptr);
  return QpAxisStruct(ptr);
}
void TaoGraphStruct::set_x2(const QpAxisStruct &src) {
  tao_graph_struct_set_x2(fortran_ptr_, src.get_fortran_ptr());
}
QpAxisStruct TaoGraphStruct::y2() const {
  void *ptr;
  tao_graph_struct_get_y2(fortran_ptr_, &ptr);
  return QpAxisStruct(ptr);
}
void TaoGraphStruct::set_y2(const QpAxisStruct &src) {
  tao_graph_struct_set_y2(fortran_ptr_, src.get_fortran_ptr());
}
QpRectStruct TaoGraphStruct::margin() const {
  void *ptr;
  tao_graph_struct_get_margin(fortran_ptr_, &ptr);
  return QpRectStruct(ptr);
}
void TaoGraphStruct::set_margin(const QpRectStruct &src) {
  tao_graph_struct_set_margin(fortran_ptr_, src.get_fortran_ptr());
}
QpRectStruct TaoGraphStruct::scale_margin() const {
  void *ptr;
  tao_graph_struct_get_scale_margin(fortran_ptr_, &ptr);
  return QpRectStruct(ptr);
}
void TaoGraphStruct::set_scale_margin(const QpRectStruct &src) {
  tao_graph_struct_set_scale_margin(fortran_ptr_, src.get_fortran_ptr());
}
double TaoGraphStruct::x_axis_scale_factor() const {
  double value;
  tao_graph_struct_get_x_axis_scale_factor(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_x_axis_scale_factor(double value) {
  tao_graph_struct_set_x_axis_scale_factor(fortran_ptr_, value);
}
double TaoGraphStruct::symbol_size_scale() const {
  double value;
  tao_graph_struct_get_symbol_size_scale(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_symbol_size_scale(double value) {
  tao_graph_struct_set_symbol_size_scale(fortran_ptr_, value);
}
FArray1D<int> TaoGraphStruct::box() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, tao_graph_struct_get_box_info);
}
void TaoGraphStruct::set_box(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_graph_struct_set_box(fortran_ptr_, v.data(), shape);
}
int TaoGraphStruct::ix_branch() const {
  int value;
  tao_graph_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_ix_branch(int value) {
  tao_graph_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoGraphStruct::ix_universe() const {
  int value;
  tao_graph_struct_get_ix_universe(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_ix_universe(int value) {
  tao_graph_struct_set_ix_universe(fortran_ptr_, value);
}
bool TaoGraphStruct::clip() const {
  bool value;
  tao_graph_struct_get_clip(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_clip(bool value) { tao_graph_struct_set_clip(fortran_ptr_, value); }
bool TaoGraphStruct::y2_mirrors_y() const {
  bool value;
  tao_graph_struct_get_y2_mirrors_y(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_y2_mirrors_y(bool value) {
  tao_graph_struct_set_y2_mirrors_y(fortran_ptr_, value);
}
bool TaoGraphStruct::limited() const {
  bool value;
  tao_graph_struct_get_limited(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_limited(bool value) { tao_graph_struct_set_limited(fortran_ptr_, value); }
bool TaoGraphStruct::draw_axes() const {
  bool value;
  tao_graph_struct_get_draw_axes(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_draw_axes(bool value) {
  tao_graph_struct_set_draw_axes(fortran_ptr_, value);
}
bool TaoGraphStruct::draw_curve_legend() const {
  bool value;
  tao_graph_struct_get_draw_curve_legend(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_draw_curve_legend(bool value) {
  tao_graph_struct_set_draw_curve_legend(fortran_ptr_, value);
}
bool TaoGraphStruct::draw_grid() const {
  bool value;
  tao_graph_struct_get_draw_grid(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_draw_grid(bool value) {
  tao_graph_struct_set_draw_grid(fortran_ptr_, value);
}
bool TaoGraphStruct::draw_title() const {
  bool value;
  tao_graph_struct_get_draw_title(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_draw_title(bool value) {
  tao_graph_struct_set_draw_title(fortran_ptr_, value);
}
bool TaoGraphStruct::draw_only_good_user_data_or_vars() const {
  bool value;
  tao_graph_struct_get_draw_only_good_user_data_or_vars(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_draw_only_good_user_data_or_vars(bool value) {
  tao_graph_struct_set_draw_only_good_user_data_or_vars(fortran_ptr_, value);
}
bool TaoGraphStruct::allow_wrap_around() const {
  bool value;
  tao_graph_struct_get_allow_wrap_around(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_allow_wrap_around(bool value) {
  tao_graph_struct_set_allow_wrap_around(fortran_ptr_, value);
}
bool TaoGraphStruct::is_valid() const {
  bool value;
  tao_graph_struct_get_is_valid(fortran_ptr_, &value);
  return value;
}
void TaoGraphStruct::set_is_valid(bool value) {
  tao_graph_struct_set_is_valid(fortran_ptr_, value);
}
std::string TaoPlotStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_plot_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotStruct::set_name(const std::string &value) {
  tao_plot_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoPlotStruct::description() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_plot_struct_get_description_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotStruct::set_description(const std::string &value) {
  tao_plot_struct_set_description(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoGraphStructAlloc1D TaoPlotStruct::graph() const {
  return TaoGraphStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_plot_struct_reallocate_graph,
      tao_plot_struct_get_graph_info
  );
}
std::optional<TaoPlotRegionStruct> TaoPlotStruct::r() const {
  void *ptr;
  tao_plot_struct_get_r(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoPlotRegionStruct(ptr);
}
void TaoPlotStruct::set_r(const TaoPlotRegionStruct &src) {
  tao_plot_struct_set_r(fortran_ptr_, src.get_fortran_ptr());
}
int TaoPlotStruct::ix_plot() const {
  int value;
  tao_plot_struct_get_ix_plot(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_ix_plot(int value) { tao_plot_struct_set_ix_plot(fortran_ptr_, value); }
int TaoPlotStruct::n_curve_pts() const {
  int value;
  tao_plot_struct_get_n_curve_pts(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_n_curve_pts(int value) {
  tao_plot_struct_set_n_curve_pts(fortran_ptr_, value);
}
std::string TaoPlotStruct::type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_plot_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotStruct::set_type(const std::string &value) {
  tao_plot_struct_set_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoPlotStruct::x_axis_type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_plot_struct_get_x_axis_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotStruct::set_x_axis_type(const std::string &value) {
  tao_plot_struct_set_x_axis_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoPlotStruct::autoscale_x() const {
  bool value;
  tao_plot_struct_get_autoscale_x(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_autoscale_x(bool value) {
  tao_plot_struct_set_autoscale_x(fortran_ptr_, value);
}
bool TaoPlotStruct::autoscale_y() const {
  bool value;
  tao_plot_struct_get_autoscale_y(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_autoscale_y(bool value) {
  tao_plot_struct_set_autoscale_y(fortran_ptr_, value);
}
bool TaoPlotStruct::autoscale_gang_x() const {
  bool value;
  tao_plot_struct_get_autoscale_gang_x(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_autoscale_gang_x(bool value) {
  tao_plot_struct_set_autoscale_gang_x(fortran_ptr_, value);
}
bool TaoPlotStruct::autoscale_gang_y() const {
  bool value;
  tao_plot_struct_get_autoscale_gang_y(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_autoscale_gang_y(bool value) {
  tao_plot_struct_set_autoscale_gang_y(fortran_ptr_, value);
}
bool TaoPlotStruct::list_with_show_plot_command() const {
  bool value;
  tao_plot_struct_get_list_with_show_plot_command(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_list_with_show_plot_command(bool value) {
  tao_plot_struct_set_list_with_show_plot_command(fortran_ptr_, value);
}
bool TaoPlotStruct::phantom() const {
  bool value;
  tao_plot_struct_get_phantom(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_phantom(bool value) { tao_plot_struct_set_phantom(fortran_ptr_, value); }
bool TaoPlotStruct::default_plot() const {
  bool value;
  tao_plot_struct_get_default_plot(fortran_ptr_, &value);
  return value;
}
void TaoPlotStruct::set_default_plot(bool value) {
  tao_plot_struct_set_default_plot(fortran_ptr_, value);
}
std::string TaoPlotRegionStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_plot_region_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoPlotRegionStruct::set_name(const std::string &value) {
  tao_plot_region_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoPlotStruct TaoPlotRegionStruct::plot() const {
  void *ptr;
  tao_plot_region_struct_get_plot(fortran_ptr_, &ptr);
  return TaoPlotStruct(ptr);
}
void TaoPlotRegionStruct::set_plot(const TaoPlotStruct &src) {
  tao_plot_region_struct_set_plot(fortran_ptr_, src.get_fortran_ptr());
}
FArray1D<double> TaoPlotRegionStruct::location() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_plot_region_struct_get_location_info);
}
void TaoPlotRegionStruct::set_location(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_plot_region_struct_set_location(fortran_ptr_, v.data(), shape);
}
bool TaoPlotRegionStruct::visible() const {
  bool value;
  tao_plot_region_struct_get_visible(fortran_ptr_, &value);
  return value;
}
void TaoPlotRegionStruct::set_visible(bool value) {
  tao_plot_region_struct_set_visible(fortran_ptr_, value);
}
bool TaoPlotRegionStruct::list_with_show_plot_command() const {
  bool value;
  tao_plot_region_struct_get_list_with_show_plot_command(fortran_ptr_, &value);
  return value;
}
void TaoPlotRegionStruct::set_list_with_show_plot_command(bool value) {
  tao_plot_region_struct_set_list_with_show_plot_command(fortran_ptr_, value);
}
bool TaoPlotRegionStruct::setup_done() const {
  bool value;
  tao_plot_region_struct_get_setup_done(fortran_ptr_, &value);
  return value;
}
void TaoPlotRegionStruct::set_setup_done(bool value) {
  tao_plot_region_struct_set_setup_done(fortran_ptr_, value);
}
std::optional<TaoUniverseStruct> TaoUniversePointerStruct::u() const {
  void *ptr;
  tao_universe_pointer_struct_get_u(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoUniverseStruct(ptr);
}
void TaoUniversePointerStruct::set_u(const TaoUniverseStruct &src) {
  tao_universe_pointer_struct_set_u(fortran_ptr_, src.get_fortran_ptr());
}
TaoGlobalStruct TaoSuperUniverseStruct::global() const {
  void *ptr;
  tao_super_universe_struct_get_global(fortran_ptr_, &ptr);
  return TaoGlobalStruct(ptr);
}
void TaoSuperUniverseStruct::set_global(const TaoGlobalStruct &src) {
  tao_super_universe_struct_set_global(fortran_ptr_, src.get_fortran_ptr());
}
TaoInitStruct TaoSuperUniverseStruct::init() const {
  void *ptr;
  tao_super_universe_struct_get_init(fortran_ptr_, &ptr);
  return TaoInitStruct(ptr);
}
void TaoSuperUniverseStruct::set_init(const TaoInitStruct &src) {
  tao_super_universe_struct_set_init(fortran_ptr_, src.get_fortran_ptr());
}
TaoCommonStruct TaoSuperUniverseStruct::com() const {
  void *ptr;
  tao_super_universe_struct_get_com(fortran_ptr_, &ptr);
  return TaoCommonStruct(ptr);
}
void TaoSuperUniverseStruct::set_com(const TaoCommonStruct &src) {
  tao_super_universe_struct_set_com(fortran_ptr_, src.get_fortran_ptr());
}
TaoPlotPageStruct TaoSuperUniverseStruct::plot_page() const {
  void *ptr;
  tao_super_universe_struct_get_plot_page(fortran_ptr_, &ptr);
  return TaoPlotPageStruct(ptr);
}
void TaoSuperUniverseStruct::set_plot_page(const TaoPlotPageStruct &src) {
  tao_super_universe_struct_set_plot_page(fortran_ptr_, src.get_fortran_ptr());
}
TaoV1VarStructAlloc1D TaoSuperUniverseStruct::v1_var() const {
  return TaoV1VarStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_super_universe_struct_reallocate_v1_var,
      tao_super_universe_struct_get_v1_var_info
  );
}
TaoVarStructAlloc1D TaoSuperUniverseStruct::var() const {
  return TaoVarStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_super_universe_struct_reallocate_var,
      tao_super_universe_struct_get_var_info
  );
}
TaoUniverseStructAlloc1D TaoSuperUniverseStruct::u() const {
  return TaoUniverseStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_super_universe_struct_reallocate_u,
      tao_super_universe_struct_get_u_info
  );
}
IntAlloc1D TaoSuperUniverseStruct::key() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_super_universe_struct_reallocate_key,
      tao_super_universe_struct_get_key_info
  );
}
void TaoSuperUniverseStruct::set_key(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_super_universe_struct_set_key(fortran_ptr_, v.data(), shape);
}
TaoBuildingWallStruct TaoSuperUniverseStruct::building_wall() const {
  void *ptr;
  tao_super_universe_struct_get_building_wall(fortran_ptr_, &ptr);
  return TaoBuildingWallStruct(ptr);
}
void TaoSuperUniverseStruct::set_building_wall(const TaoBuildingWallStruct &src) {
  tao_super_universe_struct_set_building_wall(fortran_ptr_, src.get_fortran_ptr());
}
TaoWaveStruct TaoSuperUniverseStruct::wave() const {
  void *ptr;
  tao_super_universe_struct_get_wave(fortran_ptr_, &ptr);
  return TaoWaveStruct(ptr);
}
void TaoSuperUniverseStruct::set_wave(const TaoWaveStruct &src) {
  tao_super_universe_struct_set_wave(fortran_ptr_, src.get_fortran_ptr());
}
int TaoSuperUniverseStruct::n_var_used() const {
  int value;
  tao_super_universe_struct_get_n_var_used(fortran_ptr_, &value);
  return value;
}
void TaoSuperUniverseStruct::set_n_var_used(int value) {
  tao_super_universe_struct_set_n_var_used(fortran_ptr_, value);
}
int TaoSuperUniverseStruct::n_v1_var_used() const {
  int value;
  tao_super_universe_struct_get_n_v1_var_used(fortran_ptr_, &value);
  return value;
}
void TaoSuperUniverseStruct::set_n_v1_var_used(int value) {
  tao_super_universe_struct_set_n_v1_var_used(fortran_ptr_, value);
}
TaoCmdHistoryStructArray1D TaoSuperUniverseStruct::history() const {
  return ProxyHelpers::get_type_array_1d<TaoCmdHistoryStructArray1D>(
      fortran_ptr_,
      tao_super_universe_struct_get_history_info
  );
}
bool TaoSuperUniverseStruct::initialized() const {
  bool value;
  tao_super_universe_struct_get_initialized(fortran_ptr_, &value);
  return value;
}
void TaoSuperUniverseStruct::set_initialized(bool value) {
  tao_super_universe_struct_set_initialized(fortran_ptr_, value);
}
std::string TaoVarStruct::ele_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_var_struct_get_ele_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarStruct::set_ele_name(const std::string &value) {
  tao_var_struct_set_ele_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoVarStruct::attrib_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_var_struct_get_attrib_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarStruct::set_attrib_name(const std::string &value) {
  tao_var_struct_set_attrib_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoVarStruct::id() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_var_struct_get_id_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarStruct::set_id(const std::string &value) {
  tao_var_struct_set_id(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoVarSlaveStructAlloc1D TaoVarStruct::slave() const {
  return TaoVarSlaveStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_var_struct_reallocate_slave,
      tao_var_struct_get_slave_info
  );
}
int TaoVarStruct::ix_v1() const {
  int value;
  tao_var_struct_get_ix_v1(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_ix_v1(int value) { tao_var_struct_set_ix_v1(fortran_ptr_, value); }
int TaoVarStruct::ix_var() const {
  int value;
  tao_var_struct_get_ix_var(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_ix_var(int value) { tao_var_struct_set_ix_var(fortran_ptr_, value); }
int TaoVarStruct::ix_dvar() const {
  int value;
  tao_var_struct_get_ix_dvar(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_ix_dvar(int value) { tao_var_struct_set_ix_dvar(fortran_ptr_, value); }
int TaoVarStruct::ix_attrib() const {
  int value;
  tao_var_struct_get_ix_attrib(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_ix_attrib(int value) { tao_var_struct_set_ix_attrib(fortran_ptr_, value); }
int TaoVarStruct::ix_key_table() const {
  int value;
  tao_var_struct_get_ix_key_table(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_ix_key_table(int value) {
  tao_var_struct_set_ix_key_table(fortran_ptr_, value);
}
std::optional<double> TaoVarStruct::model_value() const {
  double value;
  bool is_valid;
  tao_var_struct_get_model_value(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void TaoVarStruct::set_model_value(double value) {
  tao_var_struct_set_model_value(fortran_ptr_, value);
}
std::optional<double> TaoVarStruct::base_value() const {
  double value;
  bool is_valid;
  tao_var_struct_get_base_value(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void TaoVarStruct::set_base_value(double value) {
  tao_var_struct_set_base_value(fortran_ptr_, value);
}
double TaoVarStruct::design_value() const {
  double value;
  tao_var_struct_get_design_value(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_design_value(double value) {
  tao_var_struct_set_design_value(fortran_ptr_, value);
}
double TaoVarStruct::scratch_value() const {
  double value;
  tao_var_struct_get_scratch_value(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_scratch_value(double value) {
  tao_var_struct_set_scratch_value(fortran_ptr_, value);
}
double TaoVarStruct::old_value() const {
  double value;
  tao_var_struct_get_old_value(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_old_value(double value) {
  tao_var_struct_set_old_value(fortran_ptr_, value);
}
double TaoVarStruct::meas_value() const {
  double value;
  tao_var_struct_get_meas_value(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_meas_value(double value) {
  tao_var_struct_set_meas_value(fortran_ptr_, value);
}
double TaoVarStruct::ref_value() const {
  double value;
  tao_var_struct_get_ref_value(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_ref_value(double value) {
  tao_var_struct_set_ref_value(fortran_ptr_, value);
}
double TaoVarStruct::correction_value() const {
  double value;
  tao_var_struct_get_correction_value(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_correction_value(double value) {
  tao_var_struct_set_correction_value(fortran_ptr_, value);
}
double TaoVarStruct::high_lim() const {
  double value;
  tao_var_struct_get_high_lim(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_high_lim(double value) { tao_var_struct_set_high_lim(fortran_ptr_, value); }
double TaoVarStruct::low_lim() const {
  double value;
  tao_var_struct_get_low_lim(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_low_lim(double value) { tao_var_struct_set_low_lim(fortran_ptr_, value); }
double TaoVarStruct::step() const {
  double value;
  tao_var_struct_get_step(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_step(double value) { tao_var_struct_set_step(fortran_ptr_, value); }
double TaoVarStruct::weight() const {
  double value;
  tao_var_struct_get_weight(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_weight(double value) { tao_var_struct_set_weight(fortran_ptr_, value); }
double TaoVarStruct::delta_merit() const {
  double value;
  tao_var_struct_get_delta_merit(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_delta_merit(double value) {
  tao_var_struct_set_delta_merit(fortran_ptr_, value);
}
double TaoVarStruct::merit() const {
  double value;
  tao_var_struct_get_merit(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_merit(double value) { tao_var_struct_set_merit(fortran_ptr_, value); }
double TaoVarStruct::dMerit_dVar() const {
  double value;
  tao_var_struct_get_dMerit_dVar(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_dMerit_dVar(double value) {
  tao_var_struct_set_dMerit_dVar(fortran_ptr_, value);
}
double TaoVarStruct::key_val0() const {
  double value;
  tao_var_struct_get_key_val0(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_key_val0(double value) { tao_var_struct_set_key_val0(fortran_ptr_, value); }
double TaoVarStruct::key_delta() const {
  double value;
  tao_var_struct_get_key_delta(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_key_delta(double value) {
  tao_var_struct_set_key_delta(fortran_ptr_, value);
}
double TaoVarStruct::s() const {
  double value;
  tao_var_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_s(double value) { tao_var_struct_set_s(fortran_ptr_, value); }
double TaoVarStruct::extend_val() const {
  double value;
  tao_var_struct_get_extend_val(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_extend_val(double value) {
  tao_var_struct_set_extend_val(fortran_ptr_, value);
}
std::string TaoVarStruct::merit_type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_var_struct_get_merit_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoVarStruct::set_merit_type(const std::string &value) {
  tao_var_struct_set_merit_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoVarStruct::exists() const {
  bool value;
  tao_var_struct_get_exists(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_exists(bool value) { tao_var_struct_set_exists(fortran_ptr_, value); }
bool TaoVarStruct::good_var() const {
  bool value;
  tao_var_struct_get_good_var(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_good_var(bool value) { tao_var_struct_set_good_var(fortran_ptr_, value); }
bool TaoVarStruct::good_user() const {
  bool value;
  tao_var_struct_get_good_user(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_good_user(bool value) { tao_var_struct_set_good_user(fortran_ptr_, value); }
bool TaoVarStruct::good_opt() const {
  bool value;
  tao_var_struct_get_good_opt(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_good_opt(bool value) { tao_var_struct_set_good_opt(fortran_ptr_, value); }
bool TaoVarStruct::good_plot() const {
  bool value;
  tao_var_struct_get_good_plot(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_good_plot(bool value) { tao_var_struct_set_good_plot(fortran_ptr_, value); }
bool TaoVarStruct::useit_opt() const {
  bool value;
  tao_var_struct_get_useit_opt(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_useit_opt(bool value) { tao_var_struct_set_useit_opt(fortran_ptr_, value); }
bool TaoVarStruct::useit_plot() const {
  bool value;
  tao_var_struct_get_useit_plot(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_useit_plot(bool value) {
  tao_var_struct_set_useit_plot(fortran_ptr_, value);
}
bool TaoVarStruct::key_bound() const {
  bool value;
  tao_var_struct_get_key_bound(fortran_ptr_, &value);
  return value;
}
void TaoVarStruct::set_key_bound(bool value) { tao_var_struct_set_key_bound(fortran_ptr_, value); }
std::optional<TaoV1VarStruct> TaoVarStruct::v1() const {
  void *ptr;
  tao_var_struct_get_v1(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoV1VarStruct(ptr);
}
void TaoVarStruct::set_v1(const TaoV1VarStruct &src) {
  tao_var_struct_set_v1(fortran_ptr_, src.get_fortran_ptr());
}
int TaoVarSlaveStruct::ix_uni() const {
  int value;
  tao_var_slave_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoVarSlaveStruct::set_ix_uni(int value) {
  tao_var_slave_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoVarSlaveStruct::ix_branch() const {
  int value;
  tao_var_slave_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoVarSlaveStruct::set_ix_branch(int value) {
  tao_var_slave_struct_set_ix_branch(fortran_ptr_, value);
}
int TaoVarSlaveStruct::ix_ele() const {
  int value;
  tao_var_slave_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void TaoVarSlaveStruct::set_ix_ele(int value) {
  tao_var_slave_struct_set_ix_ele(fortran_ptr_, value);
}
std::optional<double> TaoVarSlaveStruct::model_value() const {
  double value;
  bool is_valid;
  tao_var_slave_struct_get_model_value(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void TaoVarSlaveStruct::set_model_value(double value) {
  tao_var_slave_struct_set_model_value(fortran_ptr_, value);
}
std::optional<double> TaoVarSlaveStruct::base_value() const {
  double value;
  bool is_valid;
  tao_var_slave_struct_get_base_value(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void TaoVarSlaveStruct::set_base_value(double value) {
  tao_var_slave_struct_set_base_value(fortran_ptr_, value);
}
std::string TaoLatticeStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_lattice_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoLatticeStruct::set_name(const std::string &value) {
  tao_lattice_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
LatStruct TaoLatticeStruct::lat() const {
  void *ptr;
  tao_lattice_struct_get_lat(fortran_ptr_, &ptr);
  return LatStruct(ptr);
}
void TaoLatticeStruct::set_lat(const LatStruct &src) {
  tao_lattice_struct_set_lat(fortran_ptr_, src.get_fortran_ptr());
}
LatStruct TaoLatticeStruct::high_E_lat() const {
  void *ptr;
  tao_lattice_struct_get_high_E_lat(fortran_ptr_, &ptr);
  return LatStruct(ptr);
}
void TaoLatticeStruct::set_high_E_lat(const LatStruct &src) {
  tao_lattice_struct_set_high_E_lat(fortran_ptr_, src.get_fortran_ptr());
}
LatStruct TaoLatticeStruct::low_E_lat() const {
  void *ptr;
  tao_lattice_struct_get_low_E_lat(fortran_ptr_, &ptr);
  return LatStruct(ptr);
}
void TaoLatticeStruct::set_low_E_lat(const LatStruct &src) {
  tao_lattice_struct_set_low_E_lat(fortran_ptr_, src.get_fortran_ptr());
}
RadIntAllEleStruct TaoLatticeStruct::rad_int_by_ele_ri() const {
  void *ptr;
  tao_lattice_struct_get_rad_int_by_ele_ri(fortran_ptr_, &ptr);
  return RadIntAllEleStruct(ptr);
}
void TaoLatticeStruct::set_rad_int_by_ele_ri(const RadIntAllEleStruct &src) {
  tao_lattice_struct_set_rad_int_by_ele_ri(fortran_ptr_, src.get_fortran_ptr());
}
RadIntAllEleStruct TaoLatticeStruct::rad_int_by_ele_6d() const {
  void *ptr;
  tao_lattice_struct_get_rad_int_by_ele_6d(fortran_ptr_, &ptr);
  return RadIntAllEleStruct(ptr);
}
void TaoLatticeStruct::set_rad_int_by_ele_6d(const RadIntAllEleStruct &src) {
  tao_lattice_struct_set_rad_int_by_ele_6d(fortran_ptr_, src.get_fortran_ptr());
}
TaoLatticeBranchStructAlloc1D TaoLatticeStruct::tao_branch() const {
  return TaoLatticeBranchStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_lattice_struct_reallocate_tao_branch,
      tao_lattice_struct_get_tao_branch_info
  );
}
std::string TaoBeamUniStruct::saved_at() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_beam_uni_struct_get_saved_at_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamUniStruct::set_saved_at(const std::string &value) {
  tao_beam_uni_struct_set_saved_at(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoBeamUniStruct::dump_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_beam_uni_struct_get_dump_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamUniStruct::set_dump_file(const std::string &value) {
  tao_beam_uni_struct_set_dump_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoBeamUniStruct::dump_at() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_beam_uni_struct_get_dump_at_info);
  return std::string(arr.data(), arr.size());
}
void TaoBeamUniStruct::set_dump_at(const std::string &value) {
  tao_beam_uni_struct_set_dump_at(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoBeamUniStruct::track_beam_in_universe() const {
  bool value;
  tao_beam_uni_struct_get_track_beam_in_universe(fortran_ptr_, &value);
  return value;
}
void TaoBeamUniStruct::set_track_beam_in_universe(bool value) {
  tao_beam_uni_struct_set_track_beam_in_universe(fortran_ptr_, value);
}
bool TaoBeamUniStruct::always_reinit() const {
  bool value;
  tao_beam_uni_struct_get_always_reinit(fortran_ptr_, &value);
  return value;
}
void TaoBeamUniStruct::set_always_reinit(bool value) {
  tao_beam_uni_struct_set_always_reinit(fortran_ptr_, value);
}
ApertureParamStruct TaoDynamicApertureStruct::param() const {
  void *ptr;
  tao_dynamic_aperture_struct_get_param(fortran_ptr_, &ptr);
  return ApertureParamStruct(ptr);
}
void TaoDynamicApertureStruct::set_param(const ApertureParamStruct &src) {
  tao_dynamic_aperture_struct_set_param(fortran_ptr_, src.get_fortran_ptr());
}
ApertureScanStructAlloc1D TaoDynamicApertureStruct::scan() const {
  return ApertureScanStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_dynamic_aperture_struct_reallocate_scan,
      tao_dynamic_aperture_struct_get_scan_info
  );
}
RealAlloc1D TaoDynamicApertureStruct::pz() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_dynamic_aperture_struct_reallocate_pz,
      tao_dynamic_aperture_struct_get_pz_info
  );
}
void TaoDynamicApertureStruct::set_pz(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_dynamic_aperture_struct_set_pz(fortran_ptr_, v.data(), shape);
}
double TaoDynamicApertureStruct::ellipse_scale() const {
  double value;
  tao_dynamic_aperture_struct_get_ellipse_scale(fortran_ptr_, &value);
  return value;
}
void TaoDynamicApertureStruct::set_ellipse_scale(double value) {
  tao_dynamic_aperture_struct_set_ellipse_scale(fortran_ptr_, value);
}
double TaoDynamicApertureStruct::a_emit() const {
  double value;
  tao_dynamic_aperture_struct_get_a_emit(fortran_ptr_, &value);
  return value;
}
void TaoDynamicApertureStruct::set_a_emit(double value) {
  tao_dynamic_aperture_struct_set_a_emit(fortran_ptr_, value);
}
double TaoDynamicApertureStruct::b_emit() const {
  double value;
  tao_dynamic_aperture_struct_get_b_emit(fortran_ptr_, &value);
  return value;
}
void TaoDynamicApertureStruct::set_b_emit(double value) {
  tao_dynamic_aperture_struct_set_b_emit(fortran_ptr_, value);
}
TaoModelElementStructAlloc1D TaoModelBranchStruct::ele() const {
  return TaoModelElementStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_model_branch_struct_reallocate_ele,
      tao_model_branch_struct_get_ele_info
  );
}
TaoBeamBranchStruct TaoModelBranchStruct::beam() const {
  void *ptr;
  tao_model_branch_struct_get_beam(fortran_ptr_, &ptr);
  return TaoBeamBranchStruct(ptr);
}
void TaoModelBranchStruct::set_beam(const TaoBeamBranchStruct &src) {
  tao_model_branch_struct_set_beam(fortran_ptr_, src.get_fortran_ptr());
}
bool TaoSpinMapStruct::valid() const {
  bool value;
  tao_spin_map_struct_get_valid(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapStruct::set_valid(bool value) { tao_spin_map_struct_set_valid(fortran_ptr_, value); }
SpinOrbitMap1Struct TaoSpinMapStruct::map1() const {
  void *ptr;
  tao_spin_map_struct_get_map1(fortran_ptr_, &ptr);
  return SpinOrbitMap1Struct(ptr);
}
void TaoSpinMapStruct::set_map1(const SpinOrbitMap1Struct &src) {
  tao_spin_map_struct_set_map1(fortran_ptr_, src.get_fortran_ptr());
}
SpinAxisStruct TaoSpinMapStruct::axis_input() const {
  void *ptr;
  tao_spin_map_struct_get_axis_input(fortran_ptr_, &ptr);
  return SpinAxisStruct(ptr);
}
void TaoSpinMapStruct::set_axis_input(const SpinAxisStruct &src) {
  tao_spin_map_struct_set_axis_input(fortran_ptr_, src.get_fortran_ptr());
}
SpinAxisStruct TaoSpinMapStruct::axis0() const {
  void *ptr;
  tao_spin_map_struct_get_axis0(fortran_ptr_, &ptr);
  return SpinAxisStruct(ptr);
}
void TaoSpinMapStruct::set_axis0(const SpinAxisStruct &src) {
  tao_spin_map_struct_set_axis0(fortran_ptr_, src.get_fortran_ptr());
}
SpinAxisStruct TaoSpinMapStruct::axis1() const {
  void *ptr;
  tao_spin_map_struct_get_axis1(fortran_ptr_, &ptr);
  return SpinAxisStruct(ptr);
}
void TaoSpinMapStruct::set_axis1(const SpinAxisStruct &src) {
  tao_spin_map_struct_set_axis1(fortran_ptr_, src.get_fortran_ptr());
}
int TaoSpinMapStruct::ix_ele() const {
  int value;
  tao_spin_map_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapStruct::set_ix_ele(int value) {
  tao_spin_map_struct_set_ix_ele(fortran_ptr_, value);
}
int TaoSpinMapStruct::ix_ref() const {
  int value;
  tao_spin_map_struct_get_ix_ref(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapStruct::set_ix_ref(int value) {
  tao_spin_map_struct_set_ix_ref(fortran_ptr_, value);
}
int TaoSpinMapStruct::ix_uni() const {
  int value;
  tao_spin_map_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapStruct::set_ix_uni(int value) {
  tao_spin_map_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoSpinMapStruct::ix_branch() const {
  int value;
  tao_spin_map_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoSpinMapStruct::set_ix_branch(int value) {
  tao_spin_map_struct_set_ix_branch(fortran_ptr_, value);
}
FArray2D<double> TaoSpinMapStruct::mat8() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_spin_map_struct_get_mat8_info);
}
void TaoSpinMapStruct::set_mat8(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_spin_map_struct_set_mat8, v);
}
std::string TaoDataStruct::ele_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_ele_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_ele_name(const std::string &value) {
  tao_data_struct_set_ele_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::ele_start_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_ele_start_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_ele_start_name(const std::string &value) {
  tao_data_struct_set_ele_start_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::ele_ref_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_ele_ref_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_ele_ref_name(const std::string &value) {
  tao_data_struct_set_ele_ref_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::data_type() const {
  return ProxyHelpers::get_string(fortran_ptr_, tao_data_struct_get_data_type_info);
}
void TaoDataStruct::set_data_type(const std::string &value) {
  tao_data_struct_set_data_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::merit_type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_merit_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_merit_type(const std::string &value) {
  tao_data_struct_set_merit_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::id() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_id_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_id(const std::string &value) {
  tao_data_struct_set_id(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::data_source() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_data_source_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_data_source(const std::string &value) {
  tao_data_struct_set_data_source(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoDataStruct::why_invalid() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_data_struct_get_why_invalid_info);
  return std::string(arr.data(), arr.size());
}
void TaoDataStruct::set_why_invalid(const std::string &value) {
  tao_data_struct_set_why_invalid(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoDataStruct::ix_uni() const {
  int value;
  tao_data_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_uni(int value) { tao_data_struct_set_ix_uni(fortran_ptr_, value); }
int TaoDataStruct::ix_bunch() const {
  int value;
  tao_data_struct_get_ix_bunch(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_bunch(int value) { tao_data_struct_set_ix_bunch(fortran_ptr_, value); }
int TaoDataStruct::ix_branch() const {
  int value;
  tao_data_struct_get_ix_branch(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_branch(int value) { tao_data_struct_set_ix_branch(fortran_ptr_, value); }
int TaoDataStruct::ix_ele() const {
  int value;
  tao_data_struct_get_ix_ele(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_ele(int value) { tao_data_struct_set_ix_ele(fortran_ptr_, value); }
int TaoDataStruct::ix_ele_start() const {
  int value;
  tao_data_struct_get_ix_ele_start(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_ele_start(int value) {
  tao_data_struct_set_ix_ele_start(fortran_ptr_, value);
}
int TaoDataStruct::ix_ele_ref() const {
  int value;
  tao_data_struct_get_ix_ele_ref(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_ele_ref(int value) {
  tao_data_struct_set_ix_ele_ref(fortran_ptr_, value);
}
int TaoDataStruct::ix_ele_merit() const {
  int value;
  tao_data_struct_get_ix_ele_merit(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_ele_merit(int value) {
  tao_data_struct_set_ix_ele_merit(fortran_ptr_, value);
}
int TaoDataStruct::ix_d1() const {
  int value;
  tao_data_struct_get_ix_d1(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_d1(int value) { tao_data_struct_set_ix_d1(fortran_ptr_, value); }
int TaoDataStruct::ix_data() const {
  int value;
  tao_data_struct_get_ix_data(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_data(int value) { tao_data_struct_set_ix_data(fortran_ptr_, value); }
int TaoDataStruct::ix_dModel() const {
  int value;
  tao_data_struct_get_ix_dModel(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ix_dModel(int value) { tao_data_struct_set_ix_dModel(fortran_ptr_, value); }
int TaoDataStruct::eval_point() const {
  int value;
  tao_data_struct_get_eval_point(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_eval_point(int value) {
  tao_data_struct_set_eval_point(fortran_ptr_, value);
}
double TaoDataStruct::meas_value() const {
  double value;
  tao_data_struct_get_meas_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_meas_value(double value) {
  tao_data_struct_set_meas_value(fortran_ptr_, value);
}
double TaoDataStruct::ref_value() const {
  double value;
  tao_data_struct_get_ref_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ref_value(double value) {
  tao_data_struct_set_ref_value(fortran_ptr_, value);
}
double TaoDataStruct::model_value() const {
  double value;
  tao_data_struct_get_model_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_model_value(double value) {
  tao_data_struct_set_model_value(fortran_ptr_, value);
}
double TaoDataStruct::design_value() const {
  double value;
  tao_data_struct_get_design_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_design_value(double value) {
  tao_data_struct_set_design_value(fortran_ptr_, value);
}
double TaoDataStruct::old_value() const {
  double value;
  tao_data_struct_get_old_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_old_value(double value) {
  tao_data_struct_set_old_value(fortran_ptr_, value);
}
double TaoDataStruct::base_value() const {
  double value;
  tao_data_struct_get_base_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_base_value(double value) {
  tao_data_struct_set_base_value(fortran_ptr_, value);
}
double TaoDataStruct::error_rms() const {
  double value;
  tao_data_struct_get_error_rms(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_error_rms(double value) {
  tao_data_struct_set_error_rms(fortran_ptr_, value);
}
double TaoDataStruct::delta_merit() const {
  double value;
  tao_data_struct_get_delta_merit(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_delta_merit(double value) {
  tao_data_struct_set_delta_merit(fortran_ptr_, value);
}
double TaoDataStruct::weight() const {
  double value;
  tao_data_struct_get_weight(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_weight(double value) { tao_data_struct_set_weight(fortran_ptr_, value); }
double TaoDataStruct::invalid_value() const {
  double value;
  tao_data_struct_get_invalid_value(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_invalid_value(double value) {
  tao_data_struct_set_invalid_value(fortran_ptr_, value);
}
double TaoDataStruct::merit() const {
  double value;
  tao_data_struct_get_merit(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_merit(double value) { tao_data_struct_set_merit(fortran_ptr_, value); }
double TaoDataStruct::s() const {
  double value;
  tao_data_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_s(double value) { tao_data_struct_set_s(fortran_ptr_, value); }
double TaoDataStruct::s_offset() const {
  double value;
  tao_data_struct_get_s_offset(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_s_offset(double value) {
  tao_data_struct_set_s_offset(fortran_ptr_, value);
}
double TaoDataStruct::ref_s_offset() const {
  double value;
  tao_data_struct_get_ref_s_offset(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_ref_s_offset(double value) {
  tao_data_struct_set_ref_s_offset(fortran_ptr_, value);
}
bool TaoDataStruct::err_message_printed() const {
  bool value;
  tao_data_struct_get_err_message_printed(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_err_message_printed(bool value) {
  tao_data_struct_set_err_message_printed(fortran_ptr_, value);
}
bool TaoDataStruct::exists() const {
  bool value;
  tao_data_struct_get_exists(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_exists(bool value) { tao_data_struct_set_exists(fortran_ptr_, value); }
bool TaoDataStruct::good_model() const {
  bool value;
  tao_data_struct_get_good_model(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_model(bool value) {
  tao_data_struct_set_good_model(fortran_ptr_, value);
}
bool TaoDataStruct::good_base() const {
  bool value;
  tao_data_struct_get_good_base(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_base(bool value) {
  tao_data_struct_set_good_base(fortran_ptr_, value);
}
bool TaoDataStruct::good_design() const {
  bool value;
  tao_data_struct_get_good_design(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_design(bool value) {
  tao_data_struct_set_good_design(fortran_ptr_, value);
}
bool TaoDataStruct::good_meas() const {
  bool value;
  tao_data_struct_get_good_meas(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_meas(bool value) {
  tao_data_struct_set_good_meas(fortran_ptr_, value);
}
bool TaoDataStruct::good_ref() const {
  bool value;
  tao_data_struct_get_good_ref(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_ref(bool value) { tao_data_struct_set_good_ref(fortran_ptr_, value); }
bool TaoDataStruct::good_user() const {
  bool value;
  tao_data_struct_get_good_user(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_user(bool value) {
  tao_data_struct_set_good_user(fortran_ptr_, value);
}
bool TaoDataStruct::good_opt() const {
  bool value;
  tao_data_struct_get_good_opt(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_opt(bool value) { tao_data_struct_set_good_opt(fortran_ptr_, value); }
bool TaoDataStruct::good_plot() const {
  bool value;
  tao_data_struct_get_good_plot(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_good_plot(bool value) {
  tao_data_struct_set_good_plot(fortran_ptr_, value);
}
bool TaoDataStruct::useit_plot() const {
  bool value;
  tao_data_struct_get_useit_plot(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_useit_plot(bool value) {
  tao_data_struct_set_useit_plot(fortran_ptr_, value);
}
bool TaoDataStruct::useit_opt() const {
  bool value;
  tao_data_struct_get_useit_opt(fortran_ptr_, &value);
  return value;
}
void TaoDataStruct::set_useit_opt(bool value) {
  tao_data_struct_set_useit_opt(fortran_ptr_, value);
}
TaoSpinMapStruct TaoDataStruct::spin_map() const {
  void *ptr;
  tao_data_struct_get_spin_map(fortran_ptr_, &ptr);
  return TaoSpinMapStruct(ptr);
}
void TaoDataStruct::set_spin_map(const TaoSpinMapStruct &src) {
  tao_data_struct_set_spin_map(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoD1DataStruct> TaoDataStruct::d1() const {
  void *ptr;
  tao_data_struct_get_d1(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoD1DataStruct(ptr);
}
void TaoDataStruct::set_d1(const TaoD1DataStruct &src) {
  tao_data_struct_set_d1(fortran_ptr_, src.get_fortran_ptr());
}
double TaoPingScaleStruct::a_mode_meas() const {
  double value;
  tao_ping_scale_struct_get_a_mode_meas(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleStruct::set_a_mode_meas(double value) {
  tao_ping_scale_struct_set_a_mode_meas(fortran_ptr_, value);
}
double TaoPingScaleStruct::a_mode_ref() const {
  double value;
  tao_ping_scale_struct_get_a_mode_ref(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleStruct::set_a_mode_ref(double value) {
  tao_ping_scale_struct_set_a_mode_ref(fortran_ptr_, value);
}
double TaoPingScaleStruct::b_mode_meas() const {
  double value;
  tao_ping_scale_struct_get_b_mode_meas(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleStruct::set_b_mode_meas(double value) {
  tao_ping_scale_struct_set_b_mode_meas(fortran_ptr_, value);
}
double TaoPingScaleStruct::b_mode_ref() const {
  double value;
  tao_ping_scale_struct_get_b_mode_ref(fortran_ptr_, &value);
  return value;
}
void TaoPingScaleStruct::set_b_mode_ref(double value) {
  tao_ping_scale_struct_set_b_mode_ref(fortran_ptr_, value);
}
int TaoUniverseCalcStruct::srdt_for_data() const {
  int value;
  tao_universe_calc_struct_get_srdt_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_srdt_for_data(int value) {
  tao_universe_calc_struct_set_srdt_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::rad_int_for_data() const {
  bool value;
  tao_universe_calc_struct_get_rad_int_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_rad_int_for_data(bool value) {
  tao_universe_calc_struct_set_rad_int_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::rad_int_for_plotting() const {
  bool value;
  tao_universe_calc_struct_get_rad_int_for_plotting(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_rad_int_for_plotting(bool value) {
  tao_universe_calc_struct_set_rad_int_for_plotting(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::chrom_for_data() const {
  bool value;
  tao_universe_calc_struct_get_chrom_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_chrom_for_data(bool value) {
  tao_universe_calc_struct_set_chrom_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::chrom_for_plotting() const {
  bool value;
  tao_universe_calc_struct_get_chrom_for_plotting(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_chrom_for_plotting(bool value) {
  tao_universe_calc_struct_set_chrom_for_plotting(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::lat_sigma_for_data() const {
  bool value;
  tao_universe_calc_struct_get_lat_sigma_for_data(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_lat_sigma_for_data(bool value) {
  tao_universe_calc_struct_set_lat_sigma_for_data(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::lat_sigma_for_plotting() const {
  bool value;
  tao_universe_calc_struct_get_lat_sigma_for_plotting(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_lat_sigma_for_plotting(bool value) {
  tao_universe_calc_struct_set_lat_sigma_for_plotting(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::dynamic_aperture() const {
  bool value;
  tao_universe_calc_struct_get_dynamic_aperture(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_dynamic_aperture(bool value) {
  tao_universe_calc_struct_set_dynamic_aperture(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::one_turn_map() const {
  bool value;
  tao_universe_calc_struct_get_one_turn_map(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_one_turn_map(bool value) {
  tao_universe_calc_struct_set_one_turn_map(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::lattice() const {
  bool value;
  tao_universe_calc_struct_get_lattice(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_lattice(bool value) {
  tao_universe_calc_struct_set_lattice(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::twiss() const {
  bool value;
  tao_universe_calc_struct_get_twiss(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_twiss(bool value) {
  tao_universe_calc_struct_set_twiss(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::track() const {
  bool value;
  tao_universe_calc_struct_get_track(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_track(bool value) {
  tao_universe_calc_struct_set_track(fortran_ptr_, value);
}
bool TaoUniverseCalcStruct::spin_matrices() const {
  bool value;
  tao_universe_calc_struct_get_spin_matrices(fortran_ptr_, &value);
  return value;
}
void TaoUniverseCalcStruct::set_spin_matrices(bool value) {
  tao_universe_calc_struct_set_spin_matrices(fortran_ptr_, value);
}
LatEleOrderArrayStructAlloc1D LatEleOrderStruct::branch() const {
  return LatEleOrderArrayStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      lat_ele_order_struct_reallocate_branch,
      lat_ele_order_struct_get_branch_info
  );
}
bool TaoExpressionInfoStruct::good() const {
  bool value;
  tao_expression_info_struct_get_good(fortran_ptr_, &value);
  return value;
}
void TaoExpressionInfoStruct::set_good(bool value) {
  tao_expression_info_struct_set_good(fortran_ptr_, value);
}
std::optional<EleStruct> TaoExpressionInfoStruct::ele() const {
  void *ptr;
  tao_expression_info_struct_get_ele(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleStruct(ptr);
}
void TaoExpressionInfoStruct::set_ele(const EleStruct &src) {
  tao_expression_info_struct_set_ele(fortran_ptr_, src.get_fortran_ptr());
}
double TaoExpressionInfoStruct::s() const {
  double value;
  tao_expression_info_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoExpressionInfoStruct::set_s(double value) {
  tao_expression_info_struct_set_s(fortran_ptr_, value);
}
int TaoEvalNodeStruct::type() const {
  int value;
  tao_eval_node_struct_get_type(fortran_ptr_, &value);
  return value;
}
void TaoEvalNodeStruct::set_type(int value) { tao_eval_node_struct_set_type(fortran_ptr_, value); }
std::string TaoEvalNodeStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_eval_node_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoEvalNodeStruct::set_name(const std::string &value) {
  tao_eval_node_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoEvalNodeStruct::scale() const {
  double value;
  tao_eval_node_struct_get_scale(fortran_ptr_, &value);
  return value;
}
void TaoEvalNodeStruct::set_scale(double value) {
  tao_eval_node_struct_set_scale(fortran_ptr_, value);
}
RealAlloc1D TaoEvalNodeStruct::value() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_eval_node_struct_reallocate_value,
      tao_eval_node_struct_get_value_info
  );
}
void TaoEvalNodeStruct::set_value(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_eval_node_struct_set_value(fortran_ptr_, v.data(), shape);
}
TaoExpressionInfoStructAlloc1D TaoEvalNodeStruct::info() const {
  return TaoExpressionInfoStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_eval_node_struct_reallocate_info,
      tao_eval_node_struct_get_info_info
  );
}
TaoEvalNodeStructArray1D TaoEvalNodeStruct::node() const {
  return ProxyHelpers::get_type_array_1d<TaoEvalNodeStructArray1D>(
      fortran_ptr_,
      tao_eval_node_struct_get_node_info
  );
}
std::string TaoTitleStruct::string() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_title_struct_get_string_info);
  return std::string(arr.data(), arr.size());
}
void TaoTitleStruct::set_string(const std::string &value) {
  tao_title_struct_set_string(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoTitleStruct::x() const {
  double value;
  tao_title_struct_get_x(fortran_ptr_, &value);
  return value;
}
void TaoTitleStruct::set_x(double value) { tao_title_struct_set_x(fortran_ptr_, value); }
double TaoTitleStruct::y() const {
  double value;
  tao_title_struct_get_y(fortran_ptr_, &value);
  return value;
}
void TaoTitleStruct::set_y(double value) { tao_title_struct_set_y(fortran_ptr_, value); }
std::string TaoTitleStruct::units() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_title_struct_get_units_info);
  return std::string(arr.data(), arr.size());
}
void TaoTitleStruct::set_units(const std::string &value) {
  tao_title_struct_set_units(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoTitleStruct::justify() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_title_struct_get_justify_info);
  return std::string(arr.data(), arr.size());
}
void TaoTitleStruct::set_justify(const std::string &value) {
  tao_title_struct_set_justify(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoTitleStruct::draw_it() const {
  bool value;
  tao_title_struct_get_draw_it(fortran_ptr_, &value);
  return value;
}
void TaoTitleStruct::set_draw_it(bool value) { tao_title_struct_set_draw_it(fortran_ptr_, value); }
double QpRectStruct::x1() const {
  double value;
  qp_rect_struct_get_x1(fortran_ptr_, &value);
  return value;
}
void QpRectStruct::set_x1(double value) { qp_rect_struct_set_x1(fortran_ptr_, value); }
double QpRectStruct::x2() const {
  double value;
  qp_rect_struct_get_x2(fortran_ptr_, &value);
  return value;
}
void QpRectStruct::set_x2(double value) { qp_rect_struct_set_x2(fortran_ptr_, value); }
double QpRectStruct::y1() const {
  double value;
  qp_rect_struct_get_y1(fortran_ptr_, &value);
  return value;
}
void QpRectStruct::set_y1(double value) { qp_rect_struct_set_y1(fortran_ptr_, value); }
double QpRectStruct::y2() const {
  double value;
  qp_rect_struct_get_y2(fortran_ptr_, &value);
  return value;
}
void QpRectStruct::set_y2(double value) { qp_rect_struct_set_y2(fortran_ptr_, value); }
std::string QpRectStruct::units() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_rect_struct_get_units_info);
  return std::string(arr.data(), arr.size());
}
void QpRectStruct::set_units(const std::string &value) {
  qp_rect_struct_set_units(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoEleShapeStructAlloc1D TaoDrawingStruct::ele_shape() const {
  return TaoEleShapeStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_drawing_struct_reallocate_ele_shape,
      tao_drawing_struct_get_ele_shape_info
  );
}
std::string TaoShapePatternStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_shape_pattern_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoShapePatternStruct::set_name(const std::string &value) {
  tao_shape_pattern_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
QpLineStruct TaoShapePatternStruct::line() const {
  void *ptr;
  tao_shape_pattern_struct_get_line(fortran_ptr_, &ptr);
  return QpLineStruct(ptr);
}
void TaoShapePatternStruct::set_line(const QpLineStruct &src) {
  tao_shape_pattern_struct_set_line(fortran_ptr_, src.get_fortran_ptr());
}
TaoShapePatternPointStructAlloc1D TaoShapePatternStruct::pt() const {
  return TaoShapePatternPointStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_shape_pattern_struct_reallocate_pt,
      tao_shape_pattern_struct_get_pt_info
  );
}
double TaoShapePatternPointStruct::s() const {
  double value;
  tao_shape_pattern_point_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoShapePatternPointStruct::set_s(double value) {
  tao_shape_pattern_point_struct_set_s(fortran_ptr_, value);
}
double TaoShapePatternPointStruct::y() const {
  double value;
  tao_shape_pattern_point_struct_get_y(fortran_ptr_, &value);
  return value;
}
void TaoShapePatternPointStruct::set_y(double value) {
  tao_shape_pattern_point_struct_set_y(fortran_ptr_, value);
}
double TaoShapePatternPointStruct::radius() const {
  double value;
  tao_shape_pattern_point_struct_get_radius(fortran_ptr_, &value);
  return value;
}
void TaoShapePatternPointStruct::set_radius(double value) {
  tao_shape_pattern_point_struct_set_radius(fortran_ptr_, value);
}
std::string QpAxisStruct::label() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_axis_struct_get_label_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisStruct::set_label(const std::string &value) {
  qp_axis_struct_set_label(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double QpAxisStruct::min() const {
  double value;
  qp_axis_struct_get_min(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_min(double value) { qp_axis_struct_set_min(fortran_ptr_, value); }
double QpAxisStruct::max() const {
  double value;
  qp_axis_struct_get_max(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_max(double value) { qp_axis_struct_set_max(fortran_ptr_, value); }
double QpAxisStruct::tick_min() const {
  double value;
  qp_axis_struct_get_tick_min(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_tick_min(double value) { qp_axis_struct_set_tick_min(fortran_ptr_, value); }
double QpAxisStruct::tick_max() const {
  double value;
  qp_axis_struct_get_tick_max(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_tick_max(double value) { qp_axis_struct_set_tick_max(fortran_ptr_, value); }
double QpAxisStruct::eval_min() const {
  double value;
  qp_axis_struct_get_eval_min(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_eval_min(double value) { qp_axis_struct_set_eval_min(fortran_ptr_, value); }
double QpAxisStruct::eval_max() const {
  double value;
  qp_axis_struct_get_eval_max(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_eval_max(double value) { qp_axis_struct_set_eval_max(fortran_ptr_, value); }
double QpAxisStruct::dtick() const {
  double value;
  qp_axis_struct_get_dtick(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_dtick(double value) { qp_axis_struct_set_dtick(fortran_ptr_, value); }
double QpAxisStruct::number_offset() const {
  double value;
  qp_axis_struct_get_number_offset(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_number_offset(double value) {
  qp_axis_struct_set_number_offset(fortran_ptr_, value);
}
double QpAxisStruct::label_offset() const {
  double value;
  qp_axis_struct_get_label_offset(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_label_offset(double value) {
  qp_axis_struct_set_label_offset(fortran_ptr_, value);
}
double QpAxisStruct::major_tick_len() const {
  double value;
  qp_axis_struct_get_major_tick_len(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_major_tick_len(double value) {
  qp_axis_struct_set_major_tick_len(fortran_ptr_, value);
}
double QpAxisStruct::minor_tick_len() const {
  double value;
  qp_axis_struct_get_minor_tick_len(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_minor_tick_len(double value) {
  qp_axis_struct_set_minor_tick_len(fortran_ptr_, value);
}
std::string QpAxisStruct::label_color() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_axis_struct_get_label_color_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisStruct::set_label_color(const std::string &value) {
  qp_axis_struct_set_label_color(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpAxisStruct::major_div() const {
  int value;
  qp_axis_struct_get_major_div(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_major_div(int value) { qp_axis_struct_set_major_div(fortran_ptr_, value); }
int QpAxisStruct::major_div_nominal() const {
  int value;
  qp_axis_struct_get_major_div_nominal(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_major_div_nominal(int value) {
  qp_axis_struct_set_major_div_nominal(fortran_ptr_, value);
}
int QpAxisStruct::minor_div() const {
  int value;
  qp_axis_struct_get_minor_div(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_minor_div(int value) { qp_axis_struct_set_minor_div(fortran_ptr_, value); }
int QpAxisStruct::minor_div_max() const {
  int value;
  qp_axis_struct_get_minor_div_max(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_minor_div_max(int value) {
  qp_axis_struct_set_minor_div_max(fortran_ptr_, value);
}
int QpAxisStruct::places() const {
  int value;
  qp_axis_struct_get_places(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_places(int value) { qp_axis_struct_set_places(fortran_ptr_, value); }
std::string QpAxisStruct::type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_axis_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisStruct::set_type(const std::string &value) {
  qp_axis_struct_set_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpAxisStruct::bounds() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_axis_struct_get_bounds_info);
  return std::string(arr.data(), arr.size());
}
void QpAxisStruct::set_bounds(const std::string &value) {
  qp_axis_struct_set_bounds(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpAxisStruct::tick_side() const {
  int value;
  qp_axis_struct_get_tick_side(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_tick_side(int value) { qp_axis_struct_set_tick_side(fortran_ptr_, value); }
int QpAxisStruct::number_side() const {
  int value;
  qp_axis_struct_get_number_side(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_number_side(int value) {
  qp_axis_struct_set_number_side(fortran_ptr_, value);
}
bool QpAxisStruct::draw_label() const {
  bool value;
  qp_axis_struct_get_draw_label(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_draw_label(bool value) {
  qp_axis_struct_set_draw_label(fortran_ptr_, value);
}
bool QpAxisStruct::draw_numbers() const {
  bool value;
  qp_axis_struct_get_draw_numbers(fortran_ptr_, &value);
  return value;
}
void QpAxisStruct::set_draw_numbers(bool value) {
  qp_axis_struct_set_draw_numbers(fortran_ptr_, value);
}
double QpLegendStruct::row_spacing() const {
  double value;
  qp_legend_struct_get_row_spacing(fortran_ptr_, &value);
  return value;
}
void QpLegendStruct::set_row_spacing(double value) {
  qp_legend_struct_set_row_spacing(fortran_ptr_, value);
}
double QpLegendStruct::line_length() const {
  double value;
  qp_legend_struct_get_line_length(fortran_ptr_, &value);
  return value;
}
void QpLegendStruct::set_line_length(double value) {
  qp_legend_struct_set_line_length(fortran_ptr_, value);
}
double QpLegendStruct::text_offset() const {
  double value;
  qp_legend_struct_get_text_offset(fortran_ptr_, &value);
  return value;
}
void QpLegendStruct::set_text_offset(double value) {
  qp_legend_struct_set_text_offset(fortran_ptr_, value);
}
bool QpLegendStruct::draw_line() const {
  bool value;
  qp_legend_struct_get_draw_line(fortran_ptr_, &value);
  return value;
}
void QpLegendStruct::set_draw_line(bool value) {
  qp_legend_struct_set_draw_line(fortran_ptr_, value);
}
bool QpLegendStruct::draw_symbol() const {
  bool value;
  qp_legend_struct_get_draw_symbol(fortran_ptr_, &value);
  return value;
}
void QpLegendStruct::set_draw_symbol(bool value) {
  qp_legend_struct_set_draw_symbol(fortran_ptr_, value);
}
bool QpLegendStruct::draw_text() const {
  bool value;
  qp_legend_struct_get_draw_text(fortran_ptr_, &value);
  return value;
}
void QpLegendStruct::set_draw_text(bool value) {
  qp_legend_struct_set_draw_text(fortran_ptr_, value);
}
double QpPointStruct::x() const {
  double value;
  qp_point_struct_get_x(fortran_ptr_, &value);
  return value;
}
void QpPointStruct::set_x(double value) { qp_point_struct_set_x(fortran_ptr_, value); }
double QpPointStruct::y() const {
  double value;
  qp_point_struct_get_y(fortran_ptr_, &value);
  return value;
}
void QpPointStruct::set_y(double value) { qp_point_struct_set_y(fortran_ptr_, value); }
std::string QpPointStruct::units() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_point_struct_get_units_info);
  return std::string(arr.data(), arr.size());
}
void QpPointStruct::set_units(const std::string &value) {
  qp_point_struct_set_units(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpLineStruct::width() const {
  int value;
  qp_line_struct_get_width(fortran_ptr_, &value);
  return value;
}
void QpLineStruct::set_width(int value) { qp_line_struct_set_width(fortran_ptr_, value); }
std::string QpLineStruct::color() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_line_struct_get_color_info);
  return std::string(arr.data(), arr.size());
}
void QpLineStruct::set_color(const std::string &value) {
  qp_line_struct_set_color(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpLineStruct::pattern() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_line_struct_get_pattern_info);
  return std::string(arr.data(), arr.size());
}
void QpLineStruct::set_pattern(const std::string &value) {
  qp_line_struct_set_pattern(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpSymbolStruct::type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_symbol_struct_get_type_info);
  return std::string(arr.data(), arr.size());
}
void QpSymbolStruct::set_type(const std::string &value) {
  qp_symbol_struct_set_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double QpSymbolStruct::height() const {
  double value;
  qp_symbol_struct_get_height(fortran_ptr_, &value);
  return value;
}
void QpSymbolStruct::set_height(double value) { qp_symbol_struct_set_height(fortran_ptr_, value); }
std::string QpSymbolStruct::color() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_symbol_struct_get_color_info);
  return std::string(arr.data(), arr.size());
}
void QpSymbolStruct::set_color(const std::string &value) {
  qp_symbol_struct_set_color(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string QpSymbolStruct::fill_pattern() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, qp_symbol_struct_get_fill_pattern_info);
  return std::string(arr.data(), arr.size());
}
void QpSymbolStruct::set_fill_pattern(const std::string &value) {
  qp_symbol_struct_set_fill_pattern(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int QpSymbolStruct::line_width() const {
  int value;
  qp_symbol_struct_get_line_width(fortran_ptr_, &value);
  return value;
}
void QpSymbolStruct::set_line_width(int value) {
  qp_symbol_struct_set_line_width(fortran_ptr_, value);
}
std::string TaoFloorPlanStruct::view() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_floor_plan_struct_get_view_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanStruct::set_view(const std::string &value) {
  tao_floor_plan_struct_set_view(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoFloorPlanStruct::rotation() const {
  double value;
  tao_floor_plan_struct_get_rotation(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_rotation(double value) {
  tao_floor_plan_struct_set_rotation(fortran_ptr_, value);
}
bool TaoFloorPlanStruct::correct_distortion() const {
  bool value;
  tao_floor_plan_struct_get_correct_distortion(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_correct_distortion(bool value) {
  tao_floor_plan_struct_set_correct_distortion(fortran_ptr_, value);
}
bool TaoFloorPlanStruct::flip_label_side() const {
  bool value;
  tao_floor_plan_struct_get_flip_label_side(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_flip_label_side(bool value) {
  tao_floor_plan_struct_set_flip_label_side(fortran_ptr_, value);
}
bool TaoFloorPlanStruct::size_is_absolute() const {
  bool value;
  tao_floor_plan_struct_get_size_is_absolute(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_size_is_absolute(bool value) {
  tao_floor_plan_struct_set_size_is_absolute(fortran_ptr_, value);
}
bool TaoFloorPlanStruct::draw_only_first_pass() const {
  bool value;
  tao_floor_plan_struct_get_draw_only_first_pass(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_draw_only_first_pass(bool value) {
  tao_floor_plan_struct_set_draw_only_first_pass(fortran_ptr_, value);
}
bool TaoFloorPlanStruct::draw_building_wall() const {
  bool value;
  tao_floor_plan_struct_get_draw_building_wall(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_draw_building_wall(bool value) {
  tao_floor_plan_struct_set_draw_building_wall(fortran_ptr_, value);
}
double TaoFloorPlanStruct::orbit_scale() const {
  double value;
  tao_floor_plan_struct_get_orbit_scale(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_orbit_scale(double value) {
  tao_floor_plan_struct_set_orbit_scale(fortran_ptr_, value);
}
std::string TaoFloorPlanStruct::orbit_color() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_floor_plan_struct_get_orbit_color_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanStruct::set_orbit_color(const std::string &value) {
  tao_floor_plan_struct_set_orbit_color(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoFloorPlanStruct::orbit_pattern() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_floor_plan_struct_get_orbit_pattern_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanStruct::set_orbit_pattern(const std::string &value) {
  tao_floor_plan_struct_set_orbit_pattern(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoFloorPlanStruct::orbit_lattice() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_floor_plan_struct_get_orbit_lattice_info);
  return std::string(arr.data(), arr.size());
}
void TaoFloorPlanStruct::set_orbit_lattice(const std::string &value) {
  tao_floor_plan_struct_set_orbit_lattice(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
int TaoFloorPlanStruct::orbit_width() const {
  int value;
  tao_floor_plan_struct_get_orbit_width(fortran_ptr_, &value);
  return value;
}
void TaoFloorPlanStruct::set_orbit_width(int value) {
  tao_floor_plan_struct_set_orbit_width(fortran_ptr_, value);
}
std::string TaoV1VarStruct::name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_v1_var_struct_get_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoV1VarStruct::set_name(const std::string &value) {
  tao_v1_var_struct_set_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoV1VarStruct::ix_v1_var() const {
  int value;
  tao_v1_var_struct_get_ix_v1_var(fortran_ptr_, &value);
  return value;
}
void TaoV1VarStruct::set_ix_v1_var(int value) {
  tao_v1_var_struct_set_ix_v1_var(fortran_ptr_, value);
}
TaoVarStructArray1D TaoV1VarStruct::v() const {
  return ProxyHelpers::get_type_array_1d<TaoVarStructArray1D>(
      fortran_ptr_,
      tao_v1_var_struct_get_v_info
  );
}
double TaoGlobalStruct::beam_dead_cutoff() const {
  double value;
  tao_global_struct_get_beam_dead_cutoff(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_beam_dead_cutoff(double value) {
  tao_global_struct_set_beam_dead_cutoff(fortran_ptr_, value);
}
double TaoGlobalStruct::lm_opt_deriv_reinit() const {
  double value;
  tao_global_struct_get_lm_opt_deriv_reinit(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_lm_opt_deriv_reinit(double value) {
  tao_global_struct_set_lm_opt_deriv_reinit(fortran_ptr_, value);
}
double TaoGlobalStruct::de_lm_step_ratio() const {
  double value;
  tao_global_struct_get_de_lm_step_ratio(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_de_lm_step_ratio(double value) {
  tao_global_struct_set_de_lm_step_ratio(fortran_ptr_, value);
}
double TaoGlobalStruct::de_var_to_population_factor() const {
  double value;
  tao_global_struct_get_de_var_to_population_factor(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_de_var_to_population_factor(double value) {
  tao_global_struct_set_de_var_to_population_factor(fortran_ptr_, value);
}
double TaoGlobalStruct::lmdif_eps() const {
  double value;
  tao_global_struct_get_lmdif_eps(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_lmdif_eps(double value) {
  tao_global_struct_set_lmdif_eps(fortran_ptr_, value);
}
double TaoGlobalStruct::lmdif_negligible_merit() const {
  double value;
  tao_global_struct_get_lmdif_negligible_merit(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_lmdif_negligible_merit(double value) {
  tao_global_struct_set_lmdif_negligible_merit(fortran_ptr_, value);
}
double TaoGlobalStruct::svd_cutoff() const {
  double value;
  tao_global_struct_get_svd_cutoff(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_svd_cutoff(double value) {
  tao_global_struct_set_svd_cutoff(fortran_ptr_, value);
}
double TaoGlobalStruct::unstable_penalty() const {
  double value;
  tao_global_struct_get_unstable_penalty(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_unstable_penalty(double value) {
  tao_global_struct_set_unstable_penalty(fortran_ptr_, value);
}
double TaoGlobalStruct::merit_stop_value() const {
  double value;
  tao_global_struct_get_merit_stop_value(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_merit_stop_value(double value) {
  tao_global_struct_set_merit_stop_value(fortran_ptr_, value);
}
double TaoGlobalStruct::dmerit_stop_value() const {
  double value;
  tao_global_struct_get_dmerit_stop_value(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_dmerit_stop_value(double value) {
  tao_global_struct_set_dmerit_stop_value(fortran_ptr_, value);
}
double TaoGlobalStruct::random_sigma_cutoff() const {
  double value;
  tao_global_struct_get_random_sigma_cutoff(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_random_sigma_cutoff(double value) {
  tao_global_struct_set_random_sigma_cutoff(fortran_ptr_, value);
}
double TaoGlobalStruct::delta_e_chrom() const {
  double value;
  tao_global_struct_get_delta_e_chrom(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_delta_e_chrom(double value) {
  tao_global_struct_set_delta_e_chrom(fortran_ptr_, value);
}
double TaoGlobalStruct::max_plot_time() const {
  double value;
  tao_global_struct_get_max_plot_time(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_max_plot_time(double value) {
  tao_global_struct_set_max_plot_time(fortran_ptr_, value);
}
int TaoGlobalStruct::default_universe() const {
  int value;
  tao_global_struct_get_default_universe(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_default_universe(int value) {
  tao_global_struct_set_default_universe(fortran_ptr_, value);
}
int TaoGlobalStruct::default_branch() const {
  int value;
  tao_global_struct_get_default_branch(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_default_branch(int value) {
  tao_global_struct_set_default_branch(fortran_ptr_, value);
}
int TaoGlobalStruct::n_opti_cycles() const {
  int value;
  tao_global_struct_get_n_opti_cycles(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_n_opti_cycles(int value) {
  tao_global_struct_set_n_opti_cycles(fortran_ptr_, value);
}
int TaoGlobalStruct::n_opti_loops() const {
  int value;
  tao_global_struct_get_n_opti_loops(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_n_opti_loops(int value) {
  tao_global_struct_set_n_opti_loops(fortran_ptr_, value);
}
int TaoGlobalStruct::n_threads() const {
  int value;
  tao_global_struct_get_n_threads(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_n_threads(int value) {
  tao_global_struct_set_n_threads(fortran_ptr_, value);
}
int TaoGlobalStruct::phase_units() const {
  int value;
  tao_global_struct_get_phase_units(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_phase_units(int value) {
  tao_global_struct_set_phase_units(fortran_ptr_, value);
}
int TaoGlobalStruct::bunch_to_plot() const {
  int value;
  tao_global_struct_get_bunch_to_plot(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_bunch_to_plot(int value) {
  tao_global_struct_set_bunch_to_plot(fortran_ptr_, value);
}
int TaoGlobalStruct::random_seed() const {
  int value;
  tao_global_struct_get_random_seed(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_random_seed(int value) {
  tao_global_struct_set_random_seed(fortran_ptr_, value);
}
int TaoGlobalStruct::n_top10_merit() const {
  int value;
  tao_global_struct_get_n_top10_merit(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_n_top10_merit(int value) {
  tao_global_struct_set_n_top10_merit(fortran_ptr_, value);
}
int TaoGlobalStruct::srdt_gen_n_slices() const {
  int value;
  tao_global_struct_get_srdt_gen_n_slices(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_srdt_gen_n_slices(int value) {
  tao_global_struct_set_srdt_gen_n_slices(fortran_ptr_, value);
}
int TaoGlobalStruct::datum_err_messages_max() const {
  int value;
  tao_global_struct_get_datum_err_messages_max(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_datum_err_messages_max(int value) {
  tao_global_struct_set_datum_err_messages_max(fortran_ptr_, value);
}
int TaoGlobalStruct::srdt_sxt_n_slices() const {
  int value;
  tao_global_struct_get_srdt_sxt_n_slices(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_srdt_sxt_n_slices(int value) {
  tao_global_struct_set_srdt_sxt_n_slices(fortran_ptr_, value);
}
bool TaoGlobalStruct::srdt_use_cache() const {
  bool value;
  tao_global_struct_get_srdt_use_cache(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_srdt_use_cache(bool value) {
  tao_global_struct_set_srdt_use_cache(fortran_ptr_, value);
}
std::string TaoGlobalStruct::quiet() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_quiet_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_quiet(const std::string &value) {
  tao_global_struct_set_quiet(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalStruct::random_engine() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_random_engine_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_random_engine(const std::string &value) {
  tao_global_struct_set_random_engine(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoGlobalStruct::random_gauss_converter() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_global_struct_get_random_gauss_converter_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_random_gauss_converter(const std::string &value) {
  tao_global_struct_set_random_gauss_converter(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoGlobalStruct::track_type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_track_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_track_type(const std::string &value) {
  tao_global_struct_set_track_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalStruct::lat_sigma_calc_uses_emit_from() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_global_struct_get_lat_sigma_calc_uses_emit_from_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_lat_sigma_calc_uses_emit_from(const std::string &value) {
  tao_global_struct_set_lat_sigma_calc_uses_emit_from(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoGlobalStruct::prompt_string() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_prompt_string_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_prompt_string(const std::string &value) {
  tao_global_struct_set_prompt_string(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoGlobalStruct::prompt_color() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_prompt_color_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_prompt_color(const std::string &value) {
  tao_global_struct_set_prompt_color(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalStruct::optimizer() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_optimizer_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_optimizer(const std::string &value) {
  tao_global_struct_set_optimizer(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalStruct::print_command() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_print_command_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_print_command(const std::string &value) {
  tao_global_struct_set_print_command(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoGlobalStruct::var_out_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_var_out_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_var_out_file(const std::string &value) {
  tao_global_struct_set_var_out_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoGlobalStruct::history_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_global_struct_get_history_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoGlobalStruct::set_history_file(const std::string &value) {
  tao_global_struct_set_history_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
bool TaoGlobalStruct::beam_timer_on() const {
  bool value;
  tao_global_struct_get_beam_timer_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_beam_timer_on(bool value) {
  tao_global_struct_set_beam_timer_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::box_plots() const {
  bool value;
  tao_global_struct_get_box_plots(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_box_plots(bool value) {
  tao_global_struct_set_box_plots(fortran_ptr_, value);
}
bool TaoGlobalStruct::blank_line_between_commands() const {
  bool value;
  tao_global_struct_get_blank_line_between_commands(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_blank_line_between_commands(bool value) {
  tao_global_struct_set_blank_line_between_commands(fortran_ptr_, value);
}
bool TaoGlobalStruct::cmd_file_abort_on_error() const {
  bool value;
  tao_global_struct_get_cmd_file_abort_on_error(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_cmd_file_abort_on_error(bool value) {
  tao_global_struct_set_cmd_file_abort_on_error(fortran_ptr_, value);
}
bool TaoGlobalStruct::concatenate_maps() const {
  bool value;
  tao_global_struct_get_concatenate_maps(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_concatenate_maps(bool value) {
  tao_global_struct_set_concatenate_maps(fortran_ptr_, value);
}
bool TaoGlobalStruct::derivative_recalc() const {
  bool value;
  tao_global_struct_get_derivative_recalc(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_derivative_recalc(bool value) {
  tao_global_struct_set_derivative_recalc(fortran_ptr_, value);
}
bool TaoGlobalStruct::derivative_uses_design() const {
  bool value;
  tao_global_struct_get_derivative_uses_design(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_derivative_uses_design(bool value) {
  tao_global_struct_set_derivative_uses_design(fortran_ptr_, value);
}
bool TaoGlobalStruct::disable_smooth_line_calc() const {
  bool value;
  tao_global_struct_get_disable_smooth_line_calc(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_disable_smooth_line_calc(bool value) {
  tao_global_struct_set_disable_smooth_line_calc(fortran_ptr_, value);
}
bool TaoGlobalStruct::draw_curve_off_scale_warn() const {
  bool value;
  tao_global_struct_get_draw_curve_off_scale_warn(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_draw_curve_off_scale_warn(bool value) {
  tao_global_struct_set_draw_curve_off_scale_warn(fortran_ptr_, value);
}
bool TaoGlobalStruct::external_plotting() const {
  bool value;
  tao_global_struct_get_external_plotting(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_external_plotting(bool value) {
  tao_global_struct_set_external_plotting(fortran_ptr_, value);
}
bool TaoGlobalStruct::label_lattice_elements() const {
  bool value;
  tao_global_struct_get_label_lattice_elements(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_label_lattice_elements(bool value) {
  tao_global_struct_set_label_lattice_elements(fortran_ptr_, value);
}
bool TaoGlobalStruct::label_keys() const {
  bool value;
  tao_global_struct_get_label_keys(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_label_keys(bool value) {
  tao_global_struct_set_label_keys(fortran_ptr_, value);
}
bool TaoGlobalStruct::lattice_calc_on() const {
  bool value;
  tao_global_struct_get_lattice_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_lattice_calc_on(bool value) {
  tao_global_struct_set_lattice_calc_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::only_limit_opt_vars() const {
  bool value;
  tao_global_struct_get_only_limit_opt_vars(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_only_limit_opt_vars(bool value) {
  tao_global_struct_set_only_limit_opt_vars(fortran_ptr_, value);
}
bool TaoGlobalStruct::opt_with_ref() const {
  bool value;
  tao_global_struct_get_opt_with_ref(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_opt_with_ref(bool value) {
  tao_global_struct_set_opt_with_ref(fortran_ptr_, value);
}
bool TaoGlobalStruct::opt_with_base() const {
  bool value;
  tao_global_struct_get_opt_with_base(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_opt_with_base(bool value) {
  tao_global_struct_set_opt_with_base(fortran_ptr_, value);
}
bool TaoGlobalStruct::opt_match_auto_recalc() const {
  bool value;
  tao_global_struct_get_opt_match_auto_recalc(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_opt_match_auto_recalc(bool value) {
  tao_global_struct_set_opt_match_auto_recalc(fortran_ptr_, value);
}
bool TaoGlobalStruct::opti_write_var_file() const {
  bool value;
  tao_global_struct_get_opti_write_var_file(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_opti_write_var_file(bool value) {
  tao_global_struct_set_opti_write_var_file(fortran_ptr_, value);
}
bool TaoGlobalStruct::optimizer_allow_user_abort() const {
  bool value;
  tao_global_struct_get_optimizer_allow_user_abort(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_optimizer_allow_user_abort(bool value) {
  tao_global_struct_set_optimizer_allow_user_abort(fortran_ptr_, value);
}
bool TaoGlobalStruct::optimizer_var_limit_warn() const {
  bool value;
  tao_global_struct_get_optimizer_var_limit_warn(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_optimizer_var_limit_warn(bool value) {
  tao_global_struct_set_optimizer_var_limit_warn(fortran_ptr_, value);
}
bool TaoGlobalStruct::plot_on() const {
  bool value;
  tao_global_struct_get_plot_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_plot_on(bool value) {
  tao_global_struct_set_plot_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::rad_int_user_calc_on() const {
  bool value;
  tao_global_struct_get_rad_int_user_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_rad_int_user_calc_on(bool value) {
  tao_global_struct_set_rad_int_user_calc_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::rf_on() const {
  bool value;
  tao_global_struct_get_rf_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_rf_on(bool value) { tao_global_struct_set_rf_on(fortran_ptr_, value); }
bool TaoGlobalStruct::single_step() const {
  bool value;
  tao_global_struct_get_single_step(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_single_step(bool value) {
  tao_global_struct_set_single_step(fortran_ptr_, value);
}
bool TaoGlobalStruct::stop_on_error() const {
  bool value;
  tao_global_struct_get_stop_on_error(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_stop_on_error(bool value) {
  tao_global_struct_set_stop_on_error(fortran_ptr_, value);
}
bool TaoGlobalStruct::svd_retreat_on_merit_increase() const {
  bool value;
  tao_global_struct_get_svd_retreat_on_merit_increase(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_svd_retreat_on_merit_increase(bool value) {
  tao_global_struct_set_svd_retreat_on_merit_increase(fortran_ptr_, value);
}
bool TaoGlobalStruct::var_limits_on() const {
  bool value;
  tao_global_struct_get_var_limits_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_var_limits_on(bool value) {
  tao_global_struct_set_var_limits_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::wait_for_CR_in_single_mode() const {
  bool value;
  tao_global_struct_get_wait_for_CR_in_single_mode(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_wait_for_CR_in_single_mode(bool value) {
  tao_global_struct_set_wait_for_CR_in_single_mode(fortran_ptr_, value);
}
bool TaoGlobalStruct::symbol_import() const {
  bool value;
  tao_global_struct_get_symbol_import(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_symbol_import(bool value) {
  tao_global_struct_set_symbol_import(fortran_ptr_, value);
}
bool TaoGlobalStruct::debug_on() const {
  bool value;
  tao_global_struct_get_debug_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_debug_on(bool value) {
  tao_global_struct_set_debug_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::expression_tree_on() const {
  bool value;
  tao_global_struct_get_expression_tree_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_expression_tree_on(bool value) {
  tao_global_struct_set_expression_tree_on(fortran_ptr_, value);
}
bool TaoGlobalStruct::verbose_on() const {
  bool value;
  tao_global_struct_get_verbose_on(fortran_ptr_, &value);
  return value;
}
void TaoGlobalStruct::set_verbose_on(bool value) {
  tao_global_struct_set_verbose_on(fortran_ptr_, value);
}
bool TaoInitStruct::parse_cmd_args() const {
  bool value;
  tao_init_struct_get_parse_cmd_args(fortran_ptr_, &value);
  return value;
}
void TaoInitStruct::set_parse_cmd_args(bool value) {
  tao_init_struct_set_parse_cmd_args(fortran_ptr_, value);
}
bool TaoInitStruct::debug_switch() const {
  bool value;
  tao_init_struct_get_debug_switch(fortran_ptr_, &value);
  return value;
}
void TaoInitStruct::set_debug_switch(bool value) {
  tao_init_struct_set_debug_switch(fortran_ptr_, value);
}
bool TaoInitStruct::external_plotting_switch() const {
  bool value;
  tao_init_struct_get_external_plotting_switch(fortran_ptr_, &value);
  return value;
}
void TaoInitStruct::set_external_plotting_switch(bool value) {
  tao_init_struct_set_external_plotting_switch(fortran_ptr_, value);
}
std::string TaoInitStruct::init_name() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_init_name_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_init_name(const std::string &value) {
  tao_init_struct_set_init_name(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_init_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_init_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_init_file(const std::string &value) {
  tao_init_struct_set_hook_init_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_lat_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_lat_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_lat_file(const std::string &value) {
  tao_init_struct_set_hook_lat_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_beam_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_beam_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_beam_file(const std::string &value) {
  tao_init_struct_set_hook_beam_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_data_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_data_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_data_file(const std::string &value) {
  tao_init_struct_set_hook_data_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_plot_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_plot_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_plot_file(const std::string &value) {
  tao_init_struct_set_hook_plot_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_startup_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_startup_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_startup_file(const std::string &value) {
  tao_init_struct_set_hook_startup_file(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::hook_var_file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_var_file_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_var_file(const std::string &value) {
  tao_init_struct_set_hook_var_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::hook_building_wall_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_init_struct_get_hook_building_wall_file_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_building_wall_file(const std::string &value) {
  tao_init_struct_set_hook_building_wall_file(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::init_file_arg_path() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_init_file_arg_path_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_init_file_arg_path(const std::string &value) {
  tao_init_struct_set_init_file_arg_path(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::lattice_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_lattice_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_lattice_file_arg(const std::string &value) {
  tao_init_struct_set_lattice_file_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::hook_init_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_hook_init_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_hook_init_file_arg(const std::string &value) {
  tao_init_struct_set_hook_init_file_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::init_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_init_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_init_file_arg(const std::string &value) {
  tao_init_struct_set_init_file_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::beam_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_beam_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_beam_file_arg(const std::string &value) {
  tao_init_struct_set_beam_file_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::beam_init_position_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_init_struct_get_beam_init_position_file_arg_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_beam_init_position_file_arg(const std::string &value) {
  tao_init_struct_set_beam_init_position_file_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::command_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_command_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_command_arg(const std::string &value) {
  tao_init_struct_set_command_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::data_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_data_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_data_file_arg(const std::string &value) {
  tao_init_struct_set_data_file_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::plot_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_plot_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_plot_file_arg(const std::string &value) {
  tao_init_struct_set_plot_file_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::startup_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_startup_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_startup_file_arg(const std::string &value) {
  tao_init_struct_set_startup_file_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::var_file_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_var_file_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_var_file_arg(const std::string &value) {
  tao_init_struct_set_var_file_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::building_wall_file_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_init_struct_get_building_wall_file_arg_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_building_wall_file_arg(const std::string &value) {
  tao_init_struct_set_building_wall_file_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::geometry_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_geometry_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_geometry_arg(const std::string &value) {
  tao_init_struct_set_geometry_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::slice_lattice_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_slice_lattice_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_slice_lattice_arg(const std::string &value) {
  tao_init_struct_set_slice_lattice_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::start_branch_at_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_start_branch_at_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_start_branch_at_arg(const std::string &value) {
  tao_init_struct_set_start_branch_at_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::log_startup_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_log_startup_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_log_startup_arg(const std::string &value) {
  tao_init_struct_set_log_startup_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::no_stopping_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_no_stopping_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_no_stopping_arg(const std::string &value) {
  tao_init_struct_set_no_stopping_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::noplot_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_noplot_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_noplot_arg(const std::string &value) {
  tao_init_struct_set_noplot_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::no_rad_int_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_no_rad_int_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_no_rad_int_arg(const std::string &value) {
  tao_init_struct_set_no_rad_int_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::reverse_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_reverse_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_reverse_arg(const std::string &value) {
  tao_init_struct_set_reverse_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::debug_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_debug_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_debug_arg(const std::string &value) {
  tao_init_struct_set_debug_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::disable_smooth_line_calc_arg() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_init_struct_get_disable_smooth_line_calc_arg_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_disable_smooth_line_calc_arg(const std::string &value) {
  tao_init_struct_set_disable_smooth_line_calc_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::rf_on_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_rf_on_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_rf_on_arg(const std::string &value) {
  tao_init_struct_set_rf_on_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::prompt_color_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_prompt_color_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_prompt_color_arg(const std::string &value) {
  tao_init_struct_set_prompt_color_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::quiet_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_quiet_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_quiet_arg(const std::string &value) {
  tao_init_struct_set_quiet_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::noinit_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_noinit_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_noinit_arg(const std::string &value) {
  tao_init_struct_set_noinit_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::nostartup_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_nostartup_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_nostartup_arg(const std::string &value) {
  tao_init_struct_set_nostartup_arg(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string TaoInitStruct::symbol_import_arg() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_symbol_import_arg_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_symbol_import_arg(const std::string &value) {
  tao_init_struct_set_symbol_import_arg(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoInitStruct::unique_name_suffix() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_init_struct_get_unique_name_suffix_info);
  return std::string(arr.data(), arr.size());
}
void TaoInitStruct::set_unique_name_suffix(const std::string &value) {
  tao_init_struct_set_unique_name_suffix(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
TaoPlotRegionStructAlloc1D TaoCommonStruct::plot_place_buffer() const {
  return TaoPlotRegionStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_common_struct_reallocate_plot_place_buffer,
      tao_common_struct_get_plot_place_buffer_info
  );
}
FArray2D<double> TaoCommonStruct::covar() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_common_struct_get_covar_info);
}
void TaoCommonStruct::set_covar(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_common_struct_set_covar, v);
}
FArray2D<double> TaoCommonStruct::alpha() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_common_struct_get_alpha_info);
}
void TaoCommonStruct::set_alpha(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_common_struct_set_alpha, v);
}
double TaoCommonStruct::dummy_target() const {
  double value;
  tao_common_struct_get_dummy_target(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_dummy_target(double value) {
  tao_common_struct_set_dummy_target(fortran_ptr_, value);
}
int TaoCommonStruct::n_alias() const {
  int value;
  tao_common_struct_get_n_alias(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_n_alias(int value) { tao_common_struct_set_n_alias(fortran_ptr_, value); }
int TaoCommonStruct::cmd_file_level() const {
  int value;
  tao_common_struct_get_cmd_file_level(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_cmd_file_level(int value) {
  tao_common_struct_set_cmd_file_level(fortran_ptr_, value);
}
int TaoCommonStruct::ix_key_bank() const {
  int value;
  tao_common_struct_get_ix_key_bank(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_ix_key_bank(int value) {
  tao_common_struct_set_ix_key_bank(fortran_ptr_, value);
}
int TaoCommonStruct::ix_history() const {
  int value;
  tao_common_struct_get_ix_history(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_ix_history(int value) {
  tao_common_struct_set_ix_history(fortran_ptr_, value);
}
int TaoCommonStruct::n_history() const {
  int value;
  tao_common_struct_get_n_history(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_n_history(int value) {
  tao_common_struct_set_n_history(fortran_ptr_, value);
}
int TaoCommonStruct::lev_loop() const {
  int value;
  tao_common_struct_get_lev_loop(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_lev_loop(int value) {
  tao_common_struct_set_lev_loop(fortran_ptr_, value);
}
int TaoCommonStruct::n_err_messages_printed() const {
  int value;
  tao_common_struct_get_n_err_messages_printed(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_n_err_messages_printed(int value) {
  tao_common_struct_set_n_err_messages_printed(fortran_ptr_, value);
}
int TaoCommonStruct::n_universes() const {
  int value;
  tao_common_struct_get_n_universes(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_n_universes(int value) {
  tao_common_struct_set_n_universes(fortran_ptr_, value);
}
int TaoCommonStruct::ix_beam_track_active_element() const {
  int value;
  tao_common_struct_get_ix_beam_track_active_element(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_ix_beam_track_active_element(int value) {
  tao_common_struct_set_ix_beam_track_active_element(fortran_ptr_, value);
}
bool TaoCommonStruct::cmd_file_paused() const {
  bool value;
  tao_common_struct_get_cmd_file_paused(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_cmd_file_paused(bool value) {
  tao_common_struct_set_cmd_file_paused(fortran_ptr_, value);
}
bool TaoCommonStruct::use_cmd_here() const {
  bool value;
  tao_common_struct_get_use_cmd_here(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_use_cmd_here(bool value) {
  tao_common_struct_set_use_cmd_here(fortran_ptr_, value);
}
bool TaoCommonStruct::cmd_from_cmd_file() const {
  bool value;
  tao_common_struct_get_cmd_from_cmd_file(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_cmd_from_cmd_file(bool value) {
  tao_common_struct_set_cmd_from_cmd_file(fortran_ptr_, value);
}
bool TaoCommonStruct::use_saved_beam_in_tracking() const {
  bool value;
  tao_common_struct_get_use_saved_beam_in_tracking(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_use_saved_beam_in_tracking(bool value) {
  tao_common_struct_set_use_saved_beam_in_tracking(fortran_ptr_, value);
}
bool TaoCommonStruct::single_mode() const {
  bool value;
  tao_common_struct_get_single_mode(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_single_mode(bool value) {
  tao_common_struct_set_single_mode(fortran_ptr_, value);
}
bool TaoCommonStruct::combine_consecutive_elements_of_like_name() const {
  bool value;
  tao_common_struct_get_combine_consecutive_elements_of_like_name(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_combine_consecutive_elements_of_like_name(bool value) {
  tao_common_struct_set_combine_consecutive_elements_of_like_name(fortran_ptr_, value);
}
bool TaoCommonStruct::have_tracked_beam() const {
  bool value;
  tao_common_struct_get_have_tracked_beam(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_have_tracked_beam(bool value) {
  tao_common_struct_set_have_tracked_beam(fortran_ptr_, value);
}
bool TaoCommonStruct::init_plot_needed() const {
  bool value;
  tao_common_struct_get_init_plot_needed(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_init_plot_needed(bool value) {
  tao_common_struct_set_init_plot_needed(fortran_ptr_, value);
}
bool TaoCommonStruct::init_beam() const {
  bool value;
  tao_common_struct_get_init_beam(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_init_beam(bool value) {
  tao_common_struct_set_init_beam(fortran_ptr_, value);
}
bool TaoCommonStruct::init_var() const {
  bool value;
  tao_common_struct_get_init_var(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_init_var(bool value) {
  tao_common_struct_set_init_var(fortran_ptr_, value);
}
bool TaoCommonStruct::init_read_lat_info() const {
  bool value;
  tao_common_struct_get_init_read_lat_info(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_init_read_lat_info(bool value) {
  tao_common_struct_set_init_read_lat_info(fortran_ptr_, value);
}
bool TaoCommonStruct::optimizer_running() const {
  bool value;
  tao_common_struct_get_optimizer_running(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_optimizer_running(bool value) {
  tao_common_struct_set_optimizer_running(fortran_ptr_, value);
}
bool TaoCommonStruct::have_datums_using_expressions() const {
  bool value;
  tao_common_struct_get_have_datums_using_expressions(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_have_datums_using_expressions(bool value) {
  tao_common_struct_set_have_datums_using_expressions(fortran_ptr_, value);
}
bool TaoCommonStruct::print_to_terminal() const {
  bool value;
  tao_common_struct_get_print_to_terminal(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_print_to_terminal(bool value) {
  tao_common_struct_set_print_to_terminal(fortran_ptr_, value);
}
bool TaoCommonStruct::lattice_calc_done() const {
  bool value;
  tao_common_struct_get_lattice_calc_done(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_lattice_calc_done(bool value) {
  tao_common_struct_set_lattice_calc_done(fortran_ptr_, value);
}
bool TaoCommonStruct::add_measurement_noise() const {
  bool value;
  tao_common_struct_get_add_measurement_noise(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_add_measurement_noise(bool value) {
  tao_common_struct_set_add_measurement_noise(fortran_ptr_, value);
}
FArray1D<bool> TaoCommonStruct::is_err_message_printed() const {
  return ProxyHelpers::get_array_1d<bool>(
      fortran_ptr_,
      tao_common_struct_get_is_err_message_printed_info
  );
}
void TaoCommonStruct::set_is_err_message_printed(const std::vector<bool> &v) {
  int shape[] = {static_cast<int>(v.size())};
  std::vector<int> bv(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    bv[i] = v[i] ? 1 : 0;
  tao_common_struct_set_is_err_message_printed(fortran_ptr_, bv.data(), shape);
}
bool TaoCommonStruct::command_arg_has_been_executed() const {
  bool value;
  tao_common_struct_get_command_arg_has_been_executed(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_command_arg_has_been_executed(bool value) {
  tao_common_struct_set_command_arg_has_been_executed(fortran_ptr_, value);
}
bool TaoCommonStruct::all_merit_weights_positive() const {
  bool value;
  tao_common_struct_get_all_merit_weights_positive(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_all_merit_weights_positive(bool value) {
  tao_common_struct_set_all_merit_weights_positive(fortran_ptr_, value);
}
bool TaoCommonStruct::multi_turn_orbit_is_plotted() const {
  bool value;
  tao_common_struct_get_multi_turn_orbit_is_plotted(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_multi_turn_orbit_is_plotted(bool value) {
  tao_common_struct_set_multi_turn_orbit_is_plotted(fortran_ptr_, value);
}
bool TaoCommonStruct::force_chrom_calc() const {
  bool value;
  tao_common_struct_get_force_chrom_calc(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_force_chrom_calc(bool value) {
  tao_common_struct_set_force_chrom_calc(fortran_ptr_, value);
}
bool TaoCommonStruct::force_rad_int_calc() const {
  bool value;
  tao_common_struct_get_force_rad_int_calc(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_force_rad_int_calc(bool value) {
  tao_common_struct_set_force_rad_int_calc(fortran_ptr_, value);
}
bool TaoCommonStruct::rad_int_ri_calc_on() const {
  bool value;
  tao_common_struct_get_rad_int_ri_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_rad_int_ri_calc_on(bool value) {
  tao_common_struct_set_rad_int_ri_calc_on(fortran_ptr_, value);
}
bool TaoCommonStruct::rad_int_6d_calc_on() const {
  bool value;
  tao_common_struct_get_rad_int_6d_calc_on(fortran_ptr_, &value);
  return value;
}
void TaoCommonStruct::set_rad_int_6d_calc_on(bool value) {
  tao_common_struct_set_rad_int_6d_calc_on(fortran_ptr_, value);
}
FCharArray1D TaoCommonStruct::valid_plot_who() const {
  return ProxyHelpers::get_char_array_1d(fortran_ptr_, tao_common_struct_get_valid_plot_who_info);
}
std::string TaoCommonStruct::single_mode_buffer() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_common_struct_get_single_mode_buffer_info);
  return std::string(arr.data(), arr.size());
}
void TaoCommonStruct::set_single_mode_buffer(const std::string &value) {
  tao_common_struct_set_single_mode_buffer(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoCommonStruct::cmd() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_common_struct_get_cmd_info);
  return std::string(arr.data(), arr.size());
}
void TaoCommonStruct::set_cmd(const std::string &value) {
  tao_common_struct_set_cmd(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
TaoTitleStruct TaoPlotPageStruct::title() const {
  void *ptr;
  tao_plot_page_struct_get_title(fortran_ptr_, &ptr);
  return TaoTitleStruct(ptr);
}
void TaoPlotPageStruct::set_title(const TaoTitleStruct &src) {
  tao_plot_page_struct_set_title(fortran_ptr_, src.get_fortran_ptr());
}
TaoTitleStruct TaoPlotPageStruct::subtitle() const {
  void *ptr;
  tao_plot_page_struct_get_subtitle(fortran_ptr_, &ptr);
  return TaoTitleStruct(ptr);
}
void TaoPlotPageStruct::set_subtitle(const TaoTitleStruct &src) {
  tao_plot_page_struct_set_subtitle(fortran_ptr_, src.get_fortran_ptr());
}
QpRectStruct TaoPlotPageStruct::border() const {
  void *ptr;
  tao_plot_page_struct_get_border(fortran_ptr_, &ptr);
  return QpRectStruct(ptr);
}
void TaoPlotPageStruct::set_border(const QpRectStruct &src) {
  tao_plot_page_struct_set_border(fortran_ptr_, src.get_fortran_ptr());
}
TaoDrawingStruct TaoPlotPageStruct::floor_plan() const {
  void *ptr;
  tao_plot_page_struct_get_floor_plan(fortran_ptr_, &ptr);
  return TaoDrawingStruct(ptr);
}
void TaoPlotPageStruct::set_floor_plan(const TaoDrawingStruct &src) {
  tao_plot_page_struct_set_floor_plan(fortran_ptr_, src.get_fortran_ptr());
}
TaoDrawingStruct TaoPlotPageStruct::lat_layout() const {
  void *ptr;
  tao_plot_page_struct_get_lat_layout(fortran_ptr_, &ptr);
  return TaoDrawingStruct(ptr);
}
void TaoPlotPageStruct::set_lat_layout(const TaoDrawingStruct &src) {
  tao_plot_page_struct_set_lat_layout(fortran_ptr_, src.get_fortran_ptr());
}
TaoShapePatternStructAlloc1D TaoPlotPageStruct::pattern() const {
  return TaoShapePatternStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_plot_page_struct_reallocate_pattern,
      tao_plot_page_struct_get_pattern_info
  );
}
TaoPlotStructAlloc1D TaoPlotPageStruct::template_() const {
  return TaoPlotStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_plot_page_struct_reallocate_template,
      tao_plot_page_struct_get_template_info
  );
}
TaoPlotRegionStructAlloc1D TaoPlotPageStruct::region() const {
  return TaoPlotRegionStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_plot_page_struct_reallocate_region,
      tao_plot_page_struct_get_region_info
  );
}
std::string TaoPlotPageStruct::plot_display_type() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_plot_page_struct_get_plot_display_type_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoPlotPageStruct::set_plot_display_type(const std::string &value) {
  tao_plot_page_struct_set_plot_display_type(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
FArray1D<double> TaoPlotPageStruct::size() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_plot_page_struct_get_size_info);
}
void TaoPlotPageStruct::set_size(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_plot_page_struct_set_size(fortran_ptr_, v.data(), shape);
}
double TaoPlotPageStruct::text_height() const {
  double value;
  tao_plot_page_struct_get_text_height(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_text_height(double value) {
  tao_plot_page_struct_set_text_height(fortran_ptr_, value);
}
double TaoPlotPageStruct::main_title_text_scale() const {
  double value;
  tao_plot_page_struct_get_main_title_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_main_title_text_scale(double value) {
  tao_plot_page_struct_set_main_title_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::graph_title_text_scale() const {
  double value;
  tao_plot_page_struct_get_graph_title_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_graph_title_text_scale(double value) {
  tao_plot_page_struct_set_graph_title_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::axis_number_text_scale() const {
  double value;
  tao_plot_page_struct_get_axis_number_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_axis_number_text_scale(double value) {
  tao_plot_page_struct_set_axis_number_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::axis_label_text_scale() const {
  double value;
  tao_plot_page_struct_get_axis_label_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_axis_label_text_scale(double value) {
  tao_plot_page_struct_set_axis_label_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::legend_text_scale() const {
  double value;
  tao_plot_page_struct_get_legend_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_legend_text_scale(double value) {
  tao_plot_page_struct_set_legend_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::key_table_text_scale() const {
  double value;
  tao_plot_page_struct_get_key_table_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_key_table_text_scale(double value) {
  tao_plot_page_struct_set_key_table_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::floor_plan_shape_scale() const {
  double value;
  tao_plot_page_struct_get_floor_plan_shape_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_floor_plan_shape_scale(double value) {
  tao_plot_page_struct_set_floor_plan_shape_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::floor_plan_text_scale() const {
  double value;
  tao_plot_page_struct_get_floor_plan_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_floor_plan_text_scale(double value) {
  tao_plot_page_struct_set_floor_plan_text_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::lat_layout_shape_scale() const {
  double value;
  tao_plot_page_struct_get_lat_layout_shape_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_lat_layout_shape_scale(double value) {
  tao_plot_page_struct_set_lat_layout_shape_scale(fortran_ptr_, value);
}
double TaoPlotPageStruct::lat_layout_text_scale() const {
  double value;
  tao_plot_page_struct_get_lat_layout_text_scale(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_lat_layout_text_scale(double value) {
  tao_plot_page_struct_set_lat_layout_text_scale(fortran_ptr_, value);
}
int TaoPlotPageStruct::n_curve_pts() const {
  int value;
  tao_plot_page_struct_get_n_curve_pts(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_n_curve_pts(int value) {
  tao_plot_page_struct_set_n_curve_pts(fortran_ptr_, value);
}
int TaoPlotPageStruct::id_window() const {
  int value;
  tao_plot_page_struct_get_id_window(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_id_window(int value) {
  tao_plot_page_struct_set_id_window(fortran_ptr_, value);
}
bool TaoPlotPageStruct::delete_overlapping_plots() const {
  bool value;
  tao_plot_page_struct_get_delete_overlapping_plots(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_delete_overlapping_plots(bool value) {
  tao_plot_page_struct_set_delete_overlapping_plots(fortran_ptr_, value);
}
bool TaoPlotPageStruct::draw_graph_title_suffix() const {
  bool value;
  tao_plot_page_struct_get_draw_graph_title_suffix(fortran_ptr_, &value);
  return value;
}
void TaoPlotPageStruct::set_draw_graph_title_suffix(bool value) {
  tao_plot_page_struct_set_draw_graph_title_suffix(fortran_ptr_, value);
}
TaoBuildingWallOrientationStruct TaoBuildingWallStruct::orientation() const {
  void *ptr;
  tao_building_wall_struct_get_orientation(fortran_ptr_, &ptr);
  return TaoBuildingWallOrientationStruct(ptr);
}
void TaoBuildingWallStruct::set_orientation(const TaoBuildingWallOrientationStruct &src) {
  tao_building_wall_struct_set_orientation(fortran_ptr_, src.get_fortran_ptr());
}
TaoBuildingWallSectionStructAlloc1D TaoBuildingWallStruct::section() const {
  return TaoBuildingWallSectionStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_building_wall_struct_reallocate_section,
      tao_building_wall_struct_get_section_info
  );
}
double TaoBuildingWallOrientationStruct::theta() const {
  double value;
  tao_building_wall_orientation_struct_get_theta(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallOrientationStruct::set_theta(double value) {
  tao_building_wall_orientation_struct_set_theta(fortran_ptr_, value);
}
double TaoBuildingWallOrientationStruct::x_offset() const {
  double value;
  tao_building_wall_orientation_struct_get_x_offset(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallOrientationStruct::set_x_offset(double value) {
  tao_building_wall_orientation_struct_set_x_offset(fortran_ptr_, value);
}
double TaoBuildingWallOrientationStruct::z_offset() const {
  double value;
  tao_building_wall_orientation_struct_get_z_offset(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallOrientationStruct::set_z_offset(double value) {
  tao_building_wall_orientation_struct_set_z_offset(fortran_ptr_, value);
}
std::string TaoBuildingWallSectionStruct::name() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_building_wall_section_struct_get_name_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoBuildingWallSectionStruct::set_name(const std::string &value) {
  tao_building_wall_section_struct_set_name(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
std::string TaoBuildingWallSectionStruct::constraint() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      tao_building_wall_section_struct_get_constraint_info
  );
  return std::string(arr.data(), arr.size());
}
void TaoBuildingWallSectionStruct::set_constraint(const std::string &value) {
  tao_building_wall_section_struct_set_constraint(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
TaoBuildingWallPointStructAlloc1D TaoBuildingWallSectionStruct::point() const {
  return TaoBuildingWallPointStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_building_wall_section_struct_reallocate_point,
      tao_building_wall_section_struct_get_point_info
  );
}
double TaoBuildingWallPointStruct::z() const {
  double value;
  tao_building_wall_point_struct_get_z(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointStruct::set_z(double value) {
  tao_building_wall_point_struct_set_z(fortran_ptr_, value);
}
double TaoBuildingWallPointStruct::x() const {
  double value;
  tao_building_wall_point_struct_get_x(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointStruct::set_x(double value) {
  tao_building_wall_point_struct_set_x(fortran_ptr_, value);
}
double TaoBuildingWallPointStruct::radius() const {
  double value;
  tao_building_wall_point_struct_get_radius(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointStruct::set_radius(double value) {
  tao_building_wall_point_struct_set_radius(fortran_ptr_, value);
}
double TaoBuildingWallPointStruct::z_center() const {
  double value;
  tao_building_wall_point_struct_get_z_center(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointStruct::set_z_center(double value) {
  tao_building_wall_point_struct_set_z_center(fortran_ptr_, value);
}
double TaoBuildingWallPointStruct::x_center() const {
  double value;
  tao_building_wall_point_struct_get_x_center(fortran_ptr_, &value);
  return value;
}
void TaoBuildingWallPointStruct::set_x_center(double value) {
  tao_building_wall_point_struct_set_x_center(fortran_ptr_, value);
}
std::string TaoWaveStruct::data_type() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, tao_wave_struct_get_data_type_info);
  return std::string(arr.data(), arr.size());
}
void TaoWaveStruct::set_data_type(const std::string &value) {
  tao_wave_struct_set_data_type(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TaoWaveStruct::rms_rel_a() const {
  double value;
  tao_wave_struct_get_rms_rel_a(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_a(double value) {
  tao_wave_struct_set_rms_rel_a(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_b() const {
  double value;
  tao_wave_struct_get_rms_rel_b(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_b(double value) {
  tao_wave_struct_set_rms_rel_b(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_as() const {
  double value;
  tao_wave_struct_get_rms_rel_as(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_as(double value) {
  tao_wave_struct_set_rms_rel_as(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_bs() const {
  double value;
  tao_wave_struct_get_rms_rel_bs(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_bs(double value) {
  tao_wave_struct_set_rms_rel_bs(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_ar() const {
  double value;
  tao_wave_struct_get_rms_rel_ar(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_ar(double value) {
  tao_wave_struct_set_rms_rel_ar(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_br() const {
  double value;
  tao_wave_struct_get_rms_rel_br(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_br(double value) {
  tao_wave_struct_set_rms_rel_br(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_k() const {
  double value;
  tao_wave_struct_get_rms_rel_k(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_k(double value) {
  tao_wave_struct_set_rms_rel_k(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_ks() const {
  double value;
  tao_wave_struct_get_rms_rel_ks(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_ks(double value) {
  tao_wave_struct_set_rms_rel_ks(fortran_ptr_, value);
}
double TaoWaveStruct::rms_rel_kr() const {
  double value;
  tao_wave_struct_get_rms_rel_kr(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_rel_kr(double value) {
  tao_wave_struct_set_rms_rel_kr(fortran_ptr_, value);
}
double TaoWaveStruct::rms_phi() const {
  double value;
  tao_wave_struct_get_rms_phi(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_phi(double value) { tao_wave_struct_set_rms_phi(fortran_ptr_, value); }
double TaoWaveStruct::rms_phi_s() const {
  double value;
  tao_wave_struct_get_rms_phi_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_phi_s(double value) {
  tao_wave_struct_set_rms_phi_s(fortran_ptr_, value);
}
double TaoWaveStruct::rms_phi_r() const {
  double value;
  tao_wave_struct_get_rms_phi_r(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_rms_phi_r(double value) {
  tao_wave_struct_set_rms_phi_r(fortran_ptr_, value);
}
double TaoWaveStruct::amp_ba_s() const {
  double value;
  tao_wave_struct_get_amp_ba_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_amp_ba_s(double value) {
  tao_wave_struct_set_amp_ba_s(fortran_ptr_, value);
}
double TaoWaveStruct::amp_ba_r() const {
  double value;
  tao_wave_struct_get_amp_ba_r(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_amp_ba_r(double value) {
  tao_wave_struct_set_amp_ba_r(fortran_ptr_, value);
}
double TaoWaveStruct::chi_a() const {
  double value;
  tao_wave_struct_get_chi_a(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_chi_a(double value) { tao_wave_struct_set_chi_a(fortran_ptr_, value); }
double TaoWaveStruct::chi_c() const {
  double value;
  tao_wave_struct_get_chi_c(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_chi_c(double value) { tao_wave_struct_set_chi_c(fortran_ptr_, value); }
double TaoWaveStruct::chi_ba() const {
  double value;
  tao_wave_struct_get_chi_ba(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_chi_ba(double value) { tao_wave_struct_set_chi_ba(fortran_ptr_, value); }
FArray1D<double> TaoWaveStruct::amp_a() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_wave_struct_get_amp_a_info);
}
void TaoWaveStruct::set_amp_a(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_amp_a(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoWaveStruct::amp_b() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_wave_struct_get_amp_b_info);
}
void TaoWaveStruct::set_amp_b(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_amp_b(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoWaveStruct::amp_ba() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_wave_struct_get_amp_ba_info);
}
void TaoWaveStruct::set_amp_ba(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_amp_ba(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoWaveStruct::coef_a() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_wave_struct_get_coef_a_info);
}
void TaoWaveStruct::set_coef_a(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_coef_a(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoWaveStruct::coef_b() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_wave_struct_get_coef_b_info);
}
void TaoWaveStruct::set_coef_b(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_coef_b(fortran_ptr_, v.data(), shape);
}
FArray1D<double> TaoWaveStruct::coef_ba() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, tao_wave_struct_get_coef_ba_info);
}
void TaoWaveStruct::set_coef_ba(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_coef_ba(fortran_ptr_, v.data(), shape);
}
int TaoWaveStruct::n_func() const {
  int value;
  tao_wave_struct_get_n_func(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_n_func(int value) { tao_wave_struct_set_n_func(fortran_ptr_, value); }
int TaoWaveStruct::ix_a1() const {
  int value;
  tao_wave_struct_get_ix_a1(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_ix_a1(int value) { tao_wave_struct_set_ix_a1(fortran_ptr_, value); }
int TaoWaveStruct::ix_a2() const {
  int value;
  tao_wave_struct_get_ix_a2(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_ix_a2(int value) { tao_wave_struct_set_ix_a2(fortran_ptr_, value); }
int TaoWaveStruct::ix_b1() const {
  int value;
  tao_wave_struct_get_ix_b1(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_ix_b1(int value) { tao_wave_struct_set_ix_b1(fortran_ptr_, value); }
int TaoWaveStruct::ix_b2() const {
  int value;
  tao_wave_struct_get_ix_b2(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_ix_b2(int value) { tao_wave_struct_set_ix_b2(fortran_ptr_, value); }
int TaoWaveStruct::i_a1() const {
  int value;
  tao_wave_struct_get_i_a1(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_i_a1(int value) { tao_wave_struct_set_i_a1(fortran_ptr_, value); }
int TaoWaveStruct::i_a2() const {
  int value;
  tao_wave_struct_get_i_a2(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_i_a2(int value) { tao_wave_struct_set_i_a2(fortran_ptr_, value); }
int TaoWaveStruct::i_b1() const {
  int value;
  tao_wave_struct_get_i_b1(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_i_b1(int value) { tao_wave_struct_set_i_b1(fortran_ptr_, value); }
int TaoWaveStruct::i_b2() const {
  int value;
  tao_wave_struct_get_i_b2(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_i_b2(int value) { tao_wave_struct_set_i_b2(fortran_ptr_, value); }
int TaoWaveStruct::n_a() const {
  int value;
  tao_wave_struct_get_n_a(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_n_a(int value) { tao_wave_struct_set_n_a(fortran_ptr_, value); }
int TaoWaveStruct::n_b() const {
  int value;
  tao_wave_struct_get_n_b(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_n_b(int value) { tao_wave_struct_set_n_b(fortran_ptr_, value); }
int TaoWaveStruct::i_curve_wrap_pt() const {
  int value;
  tao_wave_struct_get_i_curve_wrap_pt(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_i_curve_wrap_pt(int value) {
  tao_wave_struct_set_i_curve_wrap_pt(fortran_ptr_, value);
}
IntAlloc1D TaoWaveStruct::ix_data() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_wave_struct_reallocate_ix_data,
      tao_wave_struct_get_ix_data_info
  );
}
void TaoWaveStruct::set_ix_data(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  tao_wave_struct_set_ix_data(fortran_ptr_, v.data(), shape);
}
int TaoWaveStruct::n_kick() const {
  int value;
  tao_wave_struct_get_n_kick(fortran_ptr_, &value);
  return value;
}
void TaoWaveStruct::set_n_kick(int value) { tao_wave_struct_set_n_kick(fortran_ptr_, value); }
TaoWaveKickPtStructAlloc1D TaoWaveStruct::kick() const {
  return TaoWaveKickPtStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_wave_struct_reallocate_kick,
      tao_wave_struct_get_kick_info
  );
}
TaoGraphStruct TaoWaveStruct::base_graph() const {
  void *ptr;
  tao_wave_struct_get_base_graph(fortran_ptr_, &ptr);
  return TaoGraphStruct(ptr);
}
void TaoWaveStruct::set_base_graph(const TaoGraphStruct &src) {
  tao_wave_struct_set_base_graph(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoPlotRegionStruct> TaoWaveStruct::region() const {
  void *ptr;
  tao_wave_struct_get_region(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoPlotRegionStruct(ptr);
}
void TaoWaveStruct::set_region(const TaoPlotRegionStruct &src) {
  tao_wave_struct_set_region(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoD1DataStruct> TaoWaveStruct::d1_dat() const {
  void *ptr;
  tao_wave_struct_get_d1_dat(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoD1DataStruct(ptr);
}
void TaoWaveStruct::set_d1_dat(const TaoD1DataStruct &src) {
  tao_wave_struct_set_d1_dat(fortran_ptr_, src.get_fortran_ptr());
}
double TaoWaveKickPtStruct::phi_s() const {
  double value;
  tao_wave_kick_pt_struct_get_phi_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtStruct::set_phi_s(double value) {
  tao_wave_kick_pt_struct_set_phi_s(fortran_ptr_, value);
}
double TaoWaveKickPtStruct::phi_r() const {
  double value;
  tao_wave_kick_pt_struct_get_phi_r(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtStruct::set_phi_r(double value) {
  tao_wave_kick_pt_struct_set_phi_r(fortran_ptr_, value);
}
double TaoWaveKickPtStruct::phi() const {
  double value;
  tao_wave_kick_pt_struct_get_phi(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtStruct::set_phi(double value) {
  tao_wave_kick_pt_struct_set_phi(fortran_ptr_, value);
}
double TaoWaveKickPtStruct::amp() const {
  double value;
  tao_wave_kick_pt_struct_get_amp(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtStruct::set_amp(double value) {
  tao_wave_kick_pt_struct_set_amp(fortran_ptr_, value);
}
double TaoWaveKickPtStruct::s() const {
  double value;
  tao_wave_kick_pt_struct_get_s(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtStruct::set_s(double value) {
  tao_wave_kick_pt_struct_set_s(fortran_ptr_, value);
}
int TaoWaveKickPtStruct::ix_dat_before_kick() const {
  int value;
  tao_wave_kick_pt_struct_get_ix_dat_before_kick(fortran_ptr_, &value);
  return value;
}
void TaoWaveKickPtStruct::set_ix_dat_before_kick(int value) {
  tao_wave_kick_pt_struct_set_ix_dat_before_kick(fortran_ptr_, value);
}
std::optional<EleStruct> TaoWaveKickPtStruct::ele() const {
  void *ptr;
  tao_wave_kick_pt_struct_get_ele(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return EleStruct(ptr);
}
void TaoWaveKickPtStruct::set_ele(const EleStruct &src) {
  tao_wave_kick_pt_struct_set_ele(fortran_ptr_, src.get_fortran_ptr());
}
std::string TaoCmdHistoryStruct::cmd() const {
  return ProxyHelpers::get_string(fortran_ptr_, tao_cmd_history_struct_get_cmd_info);
}
void TaoCmdHistoryStruct::set_cmd(const std::string &value) {
  tao_cmd_history_struct_set_cmd(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int TaoCmdHistoryStruct::ix() const {
  int value;
  tao_cmd_history_struct_get_ix(fortran_ptr_, &value);
  return value;
}
void TaoCmdHistoryStruct::set_ix(int value) { tao_cmd_history_struct_set_ix(fortran_ptr_, value); }
std::optional<TaoLatticeStruct> TaoUniverseStruct::model() const {
  void *ptr;
  tao_universe_struct_get_model(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeStruct(ptr);
}
void TaoUniverseStruct::set_model(const TaoLatticeStruct &src) {
  tao_universe_struct_set_model(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoLatticeStruct> TaoUniverseStruct::design() const {
  void *ptr;
  tao_universe_struct_get_design(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeStruct(ptr);
}
void TaoUniverseStruct::set_design(const TaoLatticeStruct &src) {
  tao_universe_struct_set_design(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<TaoLatticeStruct> TaoUniverseStruct::base() const {
  void *ptr;
  tao_universe_struct_get_base(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TaoLatticeStruct(ptr);
}
void TaoUniverseStruct::set_base(const TaoLatticeStruct &src) {
  tao_universe_struct_set_base(fortran_ptr_, src.get_fortran_ptr());
}
TaoBeamUniStruct TaoUniverseStruct::beam() const {
  void *ptr;
  tao_universe_struct_get_beam(fortran_ptr_, &ptr);
  return TaoBeamUniStruct(ptr);
}
void TaoUniverseStruct::set_beam(const TaoBeamUniStruct &src) {
  tao_universe_struct_set_beam(fortran_ptr_, src.get_fortran_ptr());
}
TaoDynamicApertureStruct TaoUniverseStruct::dynamic_aperture() const {
  void *ptr;
  tao_universe_struct_get_dynamic_aperture(fortran_ptr_, &ptr);
  return TaoDynamicApertureStruct(ptr);
}
void TaoUniverseStruct::set_dynamic_aperture(const TaoDynamicApertureStruct &src) {
  tao_universe_struct_set_dynamic_aperture(fortran_ptr_, src.get_fortran_ptr());
}
TaoModelBranchStructArray1D TaoUniverseStruct::model_branch() const {
  return ProxyHelpers::get_type_array_1d<TaoModelBranchStructArray1D>(
      fortran_ptr_,
      tao_universe_struct_get_model_branch_info
  );
}
TaoD2DataStructAlloc1D TaoUniverseStruct::d2_data() const {
  return TaoD2DataStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_universe_struct_reallocate_d2_data,
      tao_universe_struct_get_d2_data_info
  );
}
TaoDataStructAlloc1D TaoUniverseStruct::data() const {
  return TaoDataStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      tao_universe_struct_reallocate_data,
      tao_universe_struct_get_data_info
  );
}
TaoPingScaleStruct TaoUniverseStruct::ping_scale() const {
  void *ptr;
  tao_universe_struct_get_ping_scale(fortran_ptr_, &ptr);
  return TaoPingScaleStruct(ptr);
}
void TaoUniverseStruct::set_ping_scale(const TaoPingScaleStruct &src) {
  tao_universe_struct_set_ping_scale(fortran_ptr_, src.get_fortran_ptr());
}
LatStruct TaoUniverseStruct::scratch_lat() const {
  void *ptr;
  tao_universe_struct_get_scratch_lat(fortran_ptr_, &ptr);
  return LatStruct(ptr);
}
void TaoUniverseStruct::set_scratch_lat(const LatStruct &src) {
  tao_universe_struct_set_scratch_lat(fortran_ptr_, src.get_fortran_ptr());
}
TaoUniverseCalcStruct TaoUniverseStruct::calc() const {
  void *ptr;
  tao_universe_struct_get_calc(fortran_ptr_, &ptr);
  return TaoUniverseCalcStruct(ptr);
}
void TaoUniverseStruct::set_calc(const TaoUniverseCalcStruct &src) {
  tao_universe_struct_set_calc(fortran_ptr_, src.get_fortran_ptr());
}
LatEleOrderStruct TaoUniverseStruct::ele_order() const {
  void *ptr;
  tao_universe_struct_get_ele_order(fortran_ptr_, &ptr);
  return LatEleOrderStruct(ptr);
}
void TaoUniverseStruct::set_ele_order(const LatEleOrderStruct &src) {
  tao_universe_struct_set_ele_order(fortran_ptr_, src.get_fortran_ptr());
}
TaoSpinMapStruct TaoUniverseStruct::spin_map() const {
  void *ptr;
  tao_universe_struct_get_spin_map(fortran_ptr_, &ptr);
  return TaoSpinMapStruct(ptr);
}
void TaoUniverseStruct::set_spin_map(const TaoSpinMapStruct &src) {
  tao_universe_struct_set_spin_map(fortran_ptr_, src.get_fortran_ptr());
}
FArray2D<double> TaoUniverseStruct::dModel_dVar() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, tao_universe_struct_get_dModel_dVar_info);
}
void TaoUniverseStruct::set_dModel_dVar(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, tao_universe_struct_set_dModel_dVar, v);
}
int TaoUniverseStruct::ix_uni() const {
  int value;
  tao_universe_struct_get_ix_uni(fortran_ptr_, &value);
  return value;
}
void TaoUniverseStruct::set_ix_uni(int value) {
  tao_universe_struct_set_ix_uni(fortran_ptr_, value);
}
int TaoUniverseStruct::n_d2_data_used() const {
  int value;
  tao_universe_struct_get_n_d2_data_used(fortran_ptr_, &value);
  return value;
}
void TaoUniverseStruct::set_n_d2_data_used(int value) {
  tao_universe_struct_set_n_d2_data_used(fortran_ptr_, value);
}
int TaoUniverseStruct::n_data_used() const {
  int value;
  tao_universe_struct_get_n_data_used(fortran_ptr_, &value);
  return value;
}
void TaoUniverseStruct::set_n_data_used(int value) {
  tao_universe_struct_set_n_data_used(fortran_ptr_, value);
}
bool TaoUniverseStruct::is_on() const {
  bool value;
  tao_universe_struct_get_is_on(fortran_ptr_, &value);
  return value;
}
void TaoUniverseStruct::set_is_on(bool value) {
  tao_universe_struct_set_is_on(fortran_ptr_, value);
}
bool TaoUniverseStruct::design_same_as_previous() const {
  bool value;
  tao_universe_struct_get_design_same_as_previous(fortran_ptr_, &value);
  return value;
}
void TaoUniverseStruct::set_design_same_as_previous(bool value) {
  tao_universe_struct_set_design_same_as_previous(fortran_ptr_, value);
}
bool TaoUniverseStruct::picked_uni() const {
  bool value;
  tao_universe_struct_get_picked_uni(fortran_ptr_, &value);
  return value;
}
void TaoUniverseStruct::set_picked_uni(bool value) {
  tao_universe_struct_set_picked_uni(fortran_ptr_, value);
}
double MadEnergyStruct::total() const {
  double value;
  mad_energy_struct_get_total(fortran_ptr_, &value);
  return value;
}
void MadEnergyStruct::set_total(double value) { mad_energy_struct_set_total(fortran_ptr_, value); }
double MadEnergyStruct::beta() const {
  double value;
  mad_energy_struct_get_beta(fortran_ptr_, &value);
  return value;
}
void MadEnergyStruct::set_beta(double value) { mad_energy_struct_set_beta(fortran_ptr_, value); }
double MadEnergyStruct::gamma() const {
  double value;
  mad_energy_struct_get_gamma(fortran_ptr_, &value);
  return value;
}
void MadEnergyStruct::set_gamma(double value) { mad_energy_struct_set_gamma(fortran_ptr_, value); }
double MadEnergyStruct::kinetic() const {
  double value;
  mad_energy_struct_get_kinetic(fortran_ptr_, &value);
  return value;
}
void MadEnergyStruct::set_kinetic(double value) {
  mad_energy_struct_set_kinetic(fortran_ptr_, value);
}
double MadEnergyStruct::p0c() const {
  double value;
  mad_energy_struct_get_p0c(fortran_ptr_, &value);
  return value;
}
void MadEnergyStruct::set_p0c(double value) { mad_energy_struct_set_p0c(fortran_ptr_, value); }
int MadEnergyStruct::particle() const {
  int value;
  mad_energy_struct_get_particle(fortran_ptr_, &value);
  return value;
}
void MadEnergyStruct::set_particle(int value) {
  mad_energy_struct_set_particle(fortran_ptr_, value);
}
FArray1D<double> MadMapStruct::k() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, mad_map_struct_get_k_info);
}
void MadMapStruct::set_k(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  mad_map_struct_set_k(fortran_ptr_, v.data(), shape);
}
FArray2D<double> MadMapStruct::r() const {
  return ProxyHelpers::get_array_2d<double>(fortran_ptr_, mad_map_struct_get_r_info);
}
void MadMapStruct::set_r(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, mad_map_struct_set_r, v);
}
FArray3D<double> MadMapStruct::t() const {
  return ProxyHelpers::get_array_3d<double>(fortran_ptr_, mad_map_struct_get_t_info);
}
void MadMapStruct::set_t(const std::vector<std::vector<std::vector<double>>> &v) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, mad_map_struct_set_t, v);
}
int64_t RandomStateStruct::ix() const {
  int64_t value;
  random_state_struct_get_ix(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_ix(int64_t value) { random_state_struct_set_ix(fortran_ptr_, value); }
int64_t RandomStateStruct::iy() const {
  int64_t value;
  random_state_struct_get_iy(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_iy(int64_t value) { random_state_struct_set_iy(fortran_ptr_, value); }
bool RandomStateStruct::number_stored() const {
  bool value;
  random_state_struct_get_number_stored(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_number_stored(bool value) {
  random_state_struct_set_number_stored(fortran_ptr_, value);
}
double RandomStateStruct::h_saved() const {
  double value;
  random_state_struct_get_h_saved(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_h_saved(double value) {
  random_state_struct_set_h_saved(fortran_ptr_, value);
}
int RandomStateStruct::engine() const {
  int value;
  random_state_struct_get_engine(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_engine(int value) {
  random_state_struct_set_engine(fortran_ptr_, value);
}
int RandomStateStruct::seed() const {
  int value;
  random_state_struct_get_seed(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_seed(int value) { random_state_struct_set_seed(fortran_ptr_, value); }
double RandomStateStruct::am() const {
  double value;
  random_state_struct_get_am(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_am(double value) { random_state_struct_set_am(fortran_ptr_, value); }
int RandomStateStruct::gauss_converter() const {
  int value;
  random_state_struct_get_gauss_converter(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_gauss_converter(int value) {
  random_state_struct_set_gauss_converter(fortran_ptr_, value);
}
double RandomStateStruct::gauss_sigma_cut() const {
  double value;
  random_state_struct_get_gauss_sigma_cut(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_gauss_sigma_cut(double value) {
  random_state_struct_set_gauss_sigma_cut(fortran_ptr_, value);
}
int64_t RandomStateStruct::in_sobseq() const {
  int64_t value;
  random_state_struct_get_in_sobseq(fortran_ptr_, &value);
  return value;
}
void RandomStateStruct::set_in_sobseq(int64_t value) {
  random_state_struct_set_in_sobseq(fortran_ptr_, value);
}
FArray1D<int64_t> RandomStateStruct::ix_sobseq() const {
  return ProxyHelpers::get_array_1d<int64_t>(fortran_ptr_, random_state_struct_get_ix_sobseq_info);
}
void RandomStateStruct::set_ix_sobseq(const std::vector<int64_t> &v) {
  int shape[] = {static_cast<int>(v.size())};
  random_state_struct_set_ix_sobseq(fortran_ptr_, v.data(), shape);
}
FArray1D<double> RandomStateStruct::x_sobseq() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, random_state_struct_get_x_sobseq_info);
}
void RandomStateStruct::set_x_sobseq(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  random_state_struct_set_x_sobseq(fortran_ptr_, v.data(), shape);
}
int BbuStageStruct::ix_ele_lr_wake() const {
  int value;
  bbu_stage_struct_get_ix_ele_lr_wake(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_ix_ele_lr_wake(int value) {
  bbu_stage_struct_set_ix_ele_lr_wake(fortran_ptr_, value);
}
int BbuStageStruct::ix_ele_stage_end() const {
  int value;
  bbu_stage_struct_get_ix_ele_stage_end(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_ix_ele_stage_end(int value) {
  bbu_stage_struct_set_ix_ele_stage_end(fortran_ptr_, value);
}
int BbuStageStruct::ix_pass() const {
  int value;
  bbu_stage_struct_get_ix_pass(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_ix_pass(int value) { bbu_stage_struct_set_ix_pass(fortran_ptr_, value); }
int BbuStageStruct::ix_stage_pass1() const {
  int value;
  bbu_stage_struct_get_ix_stage_pass1(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_ix_stage_pass1(int value) {
  bbu_stage_struct_set_ix_stage_pass1(fortran_ptr_, value);
}
int BbuStageStruct::ix_head_bunch() const {
  int value;
  bbu_stage_struct_get_ix_head_bunch(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_ix_head_bunch(int value) {
  bbu_stage_struct_set_ix_head_bunch(fortran_ptr_, value);
}
int BbuStageStruct::ix_hom_max() const {
  int value;
  bbu_stage_struct_get_ix_hom_max(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_ix_hom_max(int value) {
  bbu_stage_struct_set_ix_hom_max(fortran_ptr_, value);
}
double BbuStageStruct::hom_voltage_max() const {
  double value;
  bbu_stage_struct_get_hom_voltage_max(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_hom_voltage_max(double value) {
  bbu_stage_struct_set_hom_voltage_max(fortran_ptr_, value);
}
double BbuStageStruct::time_at_wake_ele() const {
  double value;
  bbu_stage_struct_get_time_at_wake_ele(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_time_at_wake_ele(double value) {
  bbu_stage_struct_set_time_at_wake_ele(fortran_ptr_, value);
}
FArray1D<double> BbuStageStruct::ave_orb() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bbu_stage_struct_get_ave_orb_info);
}
void BbuStageStruct::set_ave_orb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bbu_stage_struct_set_ave_orb(fortran_ptr_, v.data(), shape);
}
FArray1D<double> BbuStageStruct::rms_orb() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bbu_stage_struct_get_rms_orb_info);
}
void BbuStageStruct::set_rms_orb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bbu_stage_struct_set_rms_orb(fortran_ptr_, v.data(), shape);
}
FArray1D<double> BbuStageStruct::min_orb() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bbu_stage_struct_get_min_orb_info);
}
void BbuStageStruct::set_min_orb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bbu_stage_struct_set_min_orb(fortran_ptr_, v.data(), shape);
}
FArray1D<double> BbuStageStruct::max_orb() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bbu_stage_struct_get_max_orb_info);
}
void BbuStageStruct::set_max_orb(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bbu_stage_struct_set_max_orb(fortran_ptr_, v.data(), shape);
}
int BbuStageStruct::n_orb() const {
  int value;
  bbu_stage_struct_get_n_orb(fortran_ptr_, &value);
  return value;
}
void BbuStageStruct::set_n_orb(int value) { bbu_stage_struct_set_n_orb(fortran_ptr_, value); }
BunchStructAlloc1D BbuBeamStruct::bunch() const {
  return BunchStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bbu_beam_struct_reallocate_bunch,
      bbu_beam_struct_get_bunch_info
  );
}
BbuStageStructAlloc1D BbuBeamStruct::stage() const {
  return BbuStageStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bbu_beam_struct_reallocate_stage,
      bbu_beam_struct_get_stage_info
  );
}
IntAlloc1D BbuBeamStruct::ix_ele_bunch() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      bbu_beam_struct_reallocate_ix_ele_bunch,
      bbu_beam_struct_get_ix_ele_bunch_info
  );
}
void BbuBeamStruct::set_ix_ele_bunch(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bbu_beam_struct_set_ix_ele_bunch(fortran_ptr_, v.data(), shape);
}
int BbuBeamStruct::ix_bunch_head() const {
  int value;
  bbu_beam_struct_get_ix_bunch_head(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_ix_bunch_head(int value) {
  bbu_beam_struct_set_ix_bunch_head(fortran_ptr_, value);
}
int BbuBeamStruct::ix_bunch_end() const {
  int value;
  bbu_beam_struct_get_ix_bunch_end(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_ix_bunch_end(int value) {
  bbu_beam_struct_set_ix_bunch_end(fortran_ptr_, value);
}
int BbuBeamStruct::n_bunch_in_lat() const {
  int value;
  bbu_beam_struct_get_n_bunch_in_lat(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_n_bunch_in_lat(int value) {
  bbu_beam_struct_set_n_bunch_in_lat(fortran_ptr_, value);
}
int BbuBeamStruct::ix_stage_voltage_max() const {
  int value;
  bbu_beam_struct_get_ix_stage_voltage_max(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_ix_stage_voltage_max(int value) {
  bbu_beam_struct_set_ix_stage_voltage_max(fortran_ptr_, value);
}
double BbuBeamStruct::hom_voltage_max() const {
  double value;
  bbu_beam_struct_get_hom_voltage_max(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_hom_voltage_max(double value) {
  bbu_beam_struct_set_hom_voltage_max(fortran_ptr_, value);
}
double BbuBeamStruct::time_now() const {
  double value;
  bbu_beam_struct_get_time_now(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_time_now(double value) {
  bbu_beam_struct_set_time_now(fortran_ptr_, value);
}
double BbuBeamStruct::one_turn_time() const {
  double value;
  bbu_beam_struct_get_one_turn_time(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_one_turn_time(double value) {
  bbu_beam_struct_set_one_turn_time(fortran_ptr_, value);
}
double BbuBeamStruct::rf_wavelength_max() const {
  double value;
  bbu_beam_struct_get_rf_wavelength_max(fortran_ptr_, &value);
  return value;
}
void BbuBeamStruct::set_rf_wavelength_max(double value) {
  bbu_beam_struct_set_rf_wavelength_max(fortran_ptr_, value);
}
std::string BbuParamStruct::lat_filename() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, bbu_param_struct_get_lat_filename_info);
  return std::string(arr.data(), arr.size());
}
void BbuParamStruct::set_lat_filename(const std::string &value) {
  bbu_param_struct_set_lat_filename(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string BbuParamStruct::lat2_filename() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, bbu_param_struct_get_lat2_filename_info);
  return std::string(arr.data(), arr.size());
}
void BbuParamStruct::set_lat2_filename(const std::string &value) {
  bbu_param_struct_set_lat2_filename(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::string BbuParamStruct::bunch_by_bunch_info_file() const {
  FArray1D<char> arr = ProxyHelpers::get_array_1d<char>(
      fortran_ptr_,
      bbu_param_struct_get_bunch_by_bunch_info_file_info
  );
  return std::string(arr.data(), arr.size());
}
void BbuParamStruct::set_bunch_by_bunch_info_file(const std::string &value) {
  bbu_param_struct_set_bunch_by_bunch_info_file(
      fortran_ptr_,
      value.c_str(),
      static_cast<int>(value.length())
  );
}
bool BbuParamStruct::hybridize() const {
  bool value;
  bbu_param_struct_get_hybridize(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_hybridize(bool value) {
  bbu_param_struct_set_hybridize(fortran_ptr_, value);
}
bool BbuParamStruct::write_digested_hybrid_lat() const {
  bool value;
  bbu_param_struct_get_write_digested_hybrid_lat(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_write_digested_hybrid_lat(bool value) {
  bbu_param_struct_set_write_digested_hybrid_lat(fortran_ptr_, value);
}
bool BbuParamStruct::write_voltage_vs_time_dat() const {
  bool value;
  bbu_param_struct_get_write_voltage_vs_time_dat(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_write_voltage_vs_time_dat(bool value) {
  bbu_param_struct_set_write_voltage_vs_time_dat(fortran_ptr_, value);
}
bool BbuParamStruct::keep_overlays_and_groups() const {
  bool value;
  bbu_param_struct_get_keep_overlays_and_groups(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_keep_overlays_and_groups(bool value) {
  bbu_param_struct_set_keep_overlays_and_groups(fortran_ptr_, value);
}
bool BbuParamStruct::keep_all_lcavities() const {
  bool value;
  bbu_param_struct_get_keep_all_lcavities(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_keep_all_lcavities(bool value) {
  bbu_param_struct_set_keep_all_lcavities(fortran_ptr_, value);
}
bool BbuParamStruct::use_taylor_for_hybrids() const {
  bool value;
  bbu_param_struct_get_use_taylor_for_hybrids(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_use_taylor_for_hybrids(bool value) {
  bbu_param_struct_set_use_taylor_for_hybrids(fortran_ptr_, value);
}
bool BbuParamStruct::stable_orbit_anal() const {
  bool value;
  bbu_param_struct_get_stable_orbit_anal(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_stable_orbit_anal(bool value) {
  bbu_param_struct_set_stable_orbit_anal(fortran_ptr_, value);
}
double BbuParamStruct::limit_factor() const {
  double value;
  bbu_param_struct_get_limit_factor(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_limit_factor(double value) {
  bbu_param_struct_set_limit_factor(fortran_ptr_, value);
}
double BbuParamStruct::simulation_turns_max() const {
  double value;
  bbu_param_struct_get_simulation_turns_max(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_simulation_turns_max(double value) {
  bbu_param_struct_set_simulation_turns_max(fortran_ptr_, value);
}
double BbuParamStruct::bunch_freq() const {
  double value;
  bbu_param_struct_get_bunch_freq(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_bunch_freq(double value) {
  bbu_param_struct_set_bunch_freq(fortran_ptr_, value);
}
double BbuParamStruct::init_particle_offset() const {
  double value;
  bbu_param_struct_get_init_particle_offset(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_init_particle_offset(double value) {
  bbu_param_struct_set_init_particle_offset(fortran_ptr_, value);
}
double BbuParamStruct::current() const {
  double value;
  bbu_param_struct_get_current(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_current(double value) {
  bbu_param_struct_set_current(fortran_ptr_, value);
}
double BbuParamStruct::rel_tol() const {
  double value;
  bbu_param_struct_get_rel_tol(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_rel_tol(double value) {
  bbu_param_struct_set_rel_tol(fortran_ptr_, value);
}
bool BbuParamStruct::drscan() const {
  bool value;
  bbu_param_struct_get_drscan(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_drscan(bool value) { bbu_param_struct_set_drscan(fortran_ptr_, value); }
bool BbuParamStruct::use_interpolated_threshold() const {
  bool value;
  bbu_param_struct_get_use_interpolated_threshold(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_use_interpolated_threshold(bool value) {
  bbu_param_struct_set_use_interpolated_threshold(fortran_ptr_, value);
}
bool BbuParamStruct::write_hom_info() const {
  bool value;
  bbu_param_struct_get_write_hom_info(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_write_hom_info(bool value) {
  bbu_param_struct_set_write_hom_info(fortran_ptr_, value);
}
int BbuParamStruct::elindex() const {
  int value;
  bbu_param_struct_get_elindex(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_elindex(int value) { bbu_param_struct_set_elindex(fortran_ptr_, value); }
std::string BbuParamStruct::elname() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, bbu_param_struct_get_elname_info);
  return std::string(arr.data(), arr.size());
}
void BbuParamStruct::set_elname(const std::string &value) {
  bbu_param_struct_set_elname(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int BbuParamStruct::nstep() const {
  int value;
  bbu_param_struct_get_nstep(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_nstep(int value) { bbu_param_struct_set_nstep(fortran_ptr_, value); }
double BbuParamStruct::begdr() const {
  double value;
  bbu_param_struct_get_begdr(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_begdr(double value) { bbu_param_struct_set_begdr(fortran_ptr_, value); }
double BbuParamStruct::enddr() const {
  double value;
  bbu_param_struct_get_enddr(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_enddr(double value) { bbu_param_struct_set_enddr(fortran_ptr_, value); }
int BbuParamStruct::nrep() const {
  int value;
  bbu_param_struct_get_nrep(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_nrep(int value) { bbu_param_struct_set_nrep(fortran_ptr_, value); }
int BbuParamStruct::ran_seed() const {
  int value;
  bbu_param_struct_get_ran_seed(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_ran_seed(int value) { bbu_param_struct_set_ran_seed(fortran_ptr_, value); }
int BbuParamStruct::hom_order_cutoff() const {
  int value;
  bbu_param_struct_get_hom_order_cutoff(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_hom_order_cutoff(int value) {
  bbu_param_struct_set_hom_order_cutoff(fortran_ptr_, value);
}
double BbuParamStruct::ran_gauss_sigma_cut() const {
  double value;
  bbu_param_struct_get_ran_gauss_sigma_cut(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_ran_gauss_sigma_cut(double value) {
  bbu_param_struct_set_ran_gauss_sigma_cut(fortran_ptr_, value);
}
std::string BbuParamStruct::ele_track_end() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, bbu_param_struct_get_ele_track_end_info);
  return std::string(arr.data(), arr.size());
}
void BbuParamStruct::set_ele_track_end(const std::string &value) {
  bbu_param_struct_set_ele_track_end(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
int BbuParamStruct::ix_ele_track_end() const {
  int value;
  bbu_param_struct_get_ix_ele_track_end(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_ix_ele_track_end(int value) {
  bbu_param_struct_set_ix_ele_track_end(fortran_ptr_, value);
}
bool BbuParamStruct::regression() const {
  bool value;
  bbu_param_struct_get_regression(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_regression(bool value) {
  bbu_param_struct_set_regression(fortran_ptr_, value);
}
bool BbuParamStruct::normalize_z_to_rf() const {
  bool value;
  bbu_param_struct_get_normalize_z_to_rf(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_normalize_z_to_rf(bool value) {
  bbu_param_struct_set_normalize_z_to_rf(fortran_ptr_, value);
}
bool BbuParamStruct::ramp_on() const {
  bool value;
  bbu_param_struct_get_ramp_on(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_ramp_on(bool value) { bbu_param_struct_set_ramp_on(fortran_ptr_, value); }
FArray1D<double> BbuParamStruct::ramp_pattern() const {
  return ProxyHelpers::get_array_1d<double>(fortran_ptr_, bbu_param_struct_get_ramp_pattern_info);
}
void BbuParamStruct::set_ramp_pattern(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  bbu_param_struct_set_ramp_pattern(fortran_ptr_, v.data(), shape);
}
int BbuParamStruct::ramp_n_start() const {
  int value;
  bbu_param_struct_get_ramp_n_start(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_ramp_n_start(int value) {
  bbu_param_struct_set_ramp_n_start(fortran_ptr_, value);
}
int BbuParamStruct::n_ramp_pattern() const {
  int value;
  bbu_param_struct_get_n_ramp_pattern(fortran_ptr_, &value);
  return value;
}
void BbuParamStruct::set_n_ramp_pattern(int value) {
  bbu_param_struct_set_n_ramp_pattern(fortran_ptr_, value);
}
std::optional<int> Fibre::DIR() const {
  int value;
  bool is_valid;
  fibre_get_DIR(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_DIR(int value) { fibre_set_DIR(fortran_ptr_, value); }
std::optional<Fibre> Fibre::PREVIOUS() const {
  void *ptr;
  fibre_get_PREVIOUS(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Fibre::set_PREVIOUS(const Fibre &src) {
  fibre_set_PREVIOUS(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<Fibre> Fibre::NEXT() const {
  void *ptr;
  fibre_get_NEXT(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Fibre::set_NEXT(const Fibre &src) { fibre_set_NEXT(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<Layout> Fibre::PARENT_LAYOUT() const {
  void *ptr;
  fibre_get_PARENT_LAYOUT(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Layout(ptr);
}
void Fibre::set_PARENT_LAYOUT(const Layout &src) {
  fibre_set_PARENT_LAYOUT(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<int> Fibre::pos() const {
  int value;
  bool is_valid;
  fibre_get_pos(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_pos(int value) { fibre_set_pos(fortran_ptr_, value); }
std::optional<double> Fibre::BETA0() const {
  double value;
  bool is_valid;
  fibre_get_BETA0(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_BETA0(double value) { fibre_set_BETA0(fortran_ptr_, value); }
std::optional<double> Fibre::GAMMA0I() const {
  double value;
  bool is_valid;
  fibre_get_GAMMA0I(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_GAMMA0I(double value) { fibre_set_GAMMA0I(fortran_ptr_, value); }
std::optional<double> Fibre::GAMBET() const {
  double value;
  bool is_valid;
  fibre_get_GAMBET(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_GAMBET(double value) { fibre_set_GAMBET(fortran_ptr_, value); }
std::optional<double> Fibre::MASS() const {
  double value;
  bool is_valid;
  fibre_get_MASS(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_MASS(double value) { fibre_set_MASS(fortran_ptr_, value); }
std::optional<double> Fibre::CHARGE() const {
  double value;
  bool is_valid;
  fibre_get_CHARGE(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_CHARGE(double value) { fibre_set_CHARGE(fortran_ptr_, value); }
std::optional<double> Fibre::AG() const {
  double value;
  bool is_valid;
  fibre_get_AG(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_AG(double value) { fibre_set_AG(fortran_ptr_, value); }
std::optional<Fibre> Fibre::P() const {
  void *ptr;
  fibre_get_P(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Fibre::set_P(const Fibre &src) { fibre_set_P(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<Fibre> Fibre::N() const {
  void *ptr;
  fibre_get_N(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Fibre::set_N(const Fibre &src) { fibre_set_N(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<int> Fibre::loc() const {
  int value;
  bool is_valid;
  fibre_get_loc(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Fibre::set_loc(int value) { fibre_set_loc(fortran_ptr_, value); }
std::string Layout::NAME() const {
  return ProxyHelpers::get_string(fortran_ptr_, layout_get_NAME_info);
}
void Layout::set_NAME(const std::string &value) {
  layout_set_NAME(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
std::optional<int> Layout::INDEX() const {
  int value;
  bool is_valid;
  layout_get_INDEX(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_INDEX(int value) { layout_set_INDEX(fortran_ptr_, value); }
std::optional<double> Layout::HARMONIC_NUMBER() const {
  double value;
  bool is_valid;
  layout_get_HARMONIC_NUMBER(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_HARMONIC_NUMBER(double value) { layout_set_HARMONIC_NUMBER(fortran_ptr_, value); }
std::optional<bool> Layout::CLOSED() const {
  bool value;
  bool is_valid;
  layout_get_CLOSED(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_CLOSED(bool value) { layout_set_CLOSED(fortran_ptr_, value); }
std::optional<int> Layout::N() const {
  int value;
  bool is_valid;
  layout_get_N(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_N(int value) { layout_set_N(fortran_ptr_, value); }
std::optional<int> Layout::NTHIN() const {
  int value;
  bool is_valid;
  layout_get_NTHIN(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_NTHIN(int value) { layout_set_NTHIN(fortran_ptr_, value); }
std::optional<double> Layout::THIN() const {
  double value;
  bool is_valid;
  layout_get_THIN(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_THIN(double value) { layout_set_THIN(fortran_ptr_, value); }
std::optional<int> Layout::LASTPOS() const {
  int value;
  bool is_valid;
  layout_get_LASTPOS(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void Layout::set_LASTPOS(int value) { layout_set_LASTPOS(fortran_ptr_, value); }
std::optional<Fibre> Layout::LAST() const {
  void *ptr;
  layout_get_LAST(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Layout::set_LAST(const Fibre &src) { layout_set_LAST(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<Fibre> Layout::END() const {
  void *ptr;
  layout_get_END(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Layout::set_END(const Fibre &src) { layout_set_END(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<Fibre> Layout::START() const {
  void *ptr;
  layout_get_START(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Layout::set_START(const Fibre &src) { layout_set_START(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<Fibre> Layout::START_GROUND() const {
  void *ptr;
  layout_get_START_GROUND(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Layout::set_START_GROUND(const Fibre &src) {
  layout_set_START_GROUND(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<Fibre> Layout::END_GROUND() const {
  void *ptr;
  layout_get_END_GROUND(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Fibre(ptr);
}
void Layout::set_END_GROUND(const Fibre &src) {
  layout_set_END_GROUND(fortran_ptr_, src.get_fortran_ptr());
}
std::optional<Layout> Layout::NEXT() const {
  void *ptr;
  layout_get_NEXT(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Layout(ptr);
}
void Layout::set_NEXT(const Layout &src) { layout_set_NEXT(fortran_ptr_, src.get_fortran_ptr()); }
std::optional<Layout> Layout::PREVIOUS() const {
  void *ptr;
  layout_get_PREVIOUS(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return Layout(ptr);
}
void Layout::set_PREVIOUS(const Layout &src) {
  layout_set_PREVIOUS(fortran_ptr_, src.get_fortran_ptr());
}
double AllEncompassingStruct::real_rp_0d() const {
  double value;
  all_encompassing_struct_get_real_rp_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingStruct::set_real_rp_0d(double value) {
  all_encompassing_struct_set_real_rp_0d(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingStruct::real_rp_1d() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_1d_info
  );
}
void AllEncompassingStruct::set_real_rp_1d(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_real_rp_1d(fortran_ptr_, v.data(), shape);
}
FArray2D<double> AllEncompassingStruct::real_rp_2d() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_2d_info
  );
}
void AllEncompassingStruct::set_real_rp_2d(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, all_encompassing_struct_set_real_rp_2d, v);
}
FArray3D<double> AllEncompassingStruct::real_rp_3d() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_3d_info
  );
}
void AllEncompassingStruct::set_real_rp_3d(const std::vector<std::vector<std::vector<double>>> &v) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, all_encompassing_struct_set_real_rp_3d, v);
}
std::optional<double> AllEncompassingStruct::real_rp_0d_ptr() const {
  double value;
  bool is_valid;
  all_encompassing_struct_get_real_rp_0d_ptr(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void AllEncompassingStruct::set_real_rp_0d_ptr(double value) {
  all_encompassing_struct_set_real_rp_0d_ptr(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingStruct::real_rp_1d_ptr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_1d_ptr_info
  );
}
void AllEncompassingStruct::set_real_rp_1d_ptr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_real_rp_1d_ptr(fortran_ptr_, v.data(), shape);
}
FArray2D<double> AllEncompassingStruct::real_rp_2d_ptr() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_2d_ptr_info
  );
}
void AllEncompassingStruct::set_real_rp_2d_ptr(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, all_encompassing_struct_set_real_rp_2d_ptr, v);
}
FArray3D<double> AllEncompassingStruct::real_rp_3d_ptr() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_3d_ptr_info
  );
}
void AllEncompassingStruct::set_real_rp_3d_ptr(
    const std::vector<std::vector<std::vector<double>>> &v
) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, all_encompassing_struct_set_real_rp_3d_ptr, v);
}
RealAlloc1D AllEncompassingStruct::real_rp_1d_alloc() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      all_encompassing_struct_reallocate_real_rp_1d_alloc,
      all_encompassing_struct_get_real_rp_1d_alloc_info
  );
}
void AllEncompassingStruct::set_real_rp_1d_alloc(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_real_rp_1d_alloc(fortran_ptr_, v.data(), shape);
}
FArray2D<double> AllEncompassingStruct::real_rp_2d_alloc() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_2d_alloc_info
  );
}
void AllEncompassingStruct::set_real_rp_2d_alloc(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, all_encompassing_struct_set_real_rp_2d_alloc, v);
}
FArray3D<double> AllEncompassingStruct::real_rp_3d_alloc() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_rp_3d_alloc_info
  );
}
void AllEncompassingStruct::set_real_rp_3d_alloc(
    const std::vector<std::vector<std::vector<double>>> &v
) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, all_encompassing_struct_set_real_rp_3d_alloc, v);
}
double AllEncompassingStruct::real_dp_0d() const {
  double value;
  all_encompassing_struct_get_real_dp_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingStruct::set_real_dp_0d(double value) {
  all_encompassing_struct_set_real_dp_0d(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingStruct::real_dp_1d() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_1d_info
  );
}
void AllEncompassingStruct::set_real_dp_1d(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_real_dp_1d(fortran_ptr_, v.data(), shape);
}
FArray2D<double> AllEncompassingStruct::real_dp_2d() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_2d_info
  );
}
void AllEncompassingStruct::set_real_dp_2d(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, all_encompassing_struct_set_real_dp_2d, v);
}
FArray3D<double> AllEncompassingStruct::real_dp_3d() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_3d_info
  );
}
void AllEncompassingStruct::set_real_dp_3d(const std::vector<std::vector<std::vector<double>>> &v) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, all_encompassing_struct_set_real_dp_3d, v);
}
std::optional<double> AllEncompassingStruct::real_dp_0d_ptr() const {
  double value;
  bool is_valid;
  all_encompassing_struct_get_real_dp_0d_ptr(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void AllEncompassingStruct::set_real_dp_0d_ptr(double value) {
  all_encompassing_struct_set_real_dp_0d_ptr(fortran_ptr_, value);
}
FArray1D<double> AllEncompassingStruct::real_dp_1d_ptr() const {
  return ProxyHelpers::get_array_1d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_1d_ptr_info
  );
}
void AllEncompassingStruct::set_real_dp_1d_ptr(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_real_dp_1d_ptr(fortran_ptr_, v.data(), shape);
}
FArray2D<double> AllEncompassingStruct::real_dp_2d_ptr() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_2d_ptr_info
  );
}
void AllEncompassingStruct::set_real_dp_2d_ptr(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, all_encompassing_struct_set_real_dp_2d_ptr, v);
}
FArray3D<double> AllEncompassingStruct::real_dp_3d_ptr() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_3d_ptr_info
  );
}
void AllEncompassingStruct::set_real_dp_3d_ptr(
    const std::vector<std::vector<std::vector<double>>> &v
) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, all_encompassing_struct_set_real_dp_3d_ptr, v);
}
RealAlloc1D AllEncompassingStruct::real_dp_1d_alloc() const {
  return RealAlloc1D(
      const_cast<void *>(fortran_ptr_),
      all_encompassing_struct_reallocate_real_dp_1d_alloc,
      all_encompassing_struct_get_real_dp_1d_alloc_info
  );
}
void AllEncompassingStruct::set_real_dp_1d_alloc(const std::vector<double> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_real_dp_1d_alloc(fortran_ptr_, v.data(), shape);
}
FArray2D<double> AllEncompassingStruct::real_dp_2d_alloc() const {
  return ProxyHelpers::get_array_2d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_2d_alloc_info
  );
}
void AllEncompassingStruct::set_real_dp_2d_alloc(const std::vector<std::vector<double>> &v) {
  ProxyHelpers::set_array_2d<double>(fortran_ptr_, all_encompassing_struct_set_real_dp_2d_alloc, v);
}
FArray3D<double> AllEncompassingStruct::real_dp_3d_alloc() const {
  return ProxyHelpers::get_array_3d<double>(
      fortran_ptr_,
      all_encompassing_struct_get_real_dp_3d_alloc_info
  );
}
void AllEncompassingStruct::set_real_dp_3d_alloc(
    const std::vector<std::vector<std::vector<double>>> &v
) {
  ProxyHelpers::set_array_3d<double>(fortran_ptr_, all_encompassing_struct_set_real_dp_3d_alloc, v);
}
std::complex<double> AllEncompassingStruct::complex_dp_0d() const {
  std::complex<double> c_value;
  all_encompassing_struct_get_complex_dp_0d(fortran_ptr_, &c_value);
  return c_value;
}
void AllEncompassingStruct::set_complex_dp_0d(std::complex<double> value) {
  all_encompassing_struct_set_complex_dp_0d(fortran_ptr_, value);
}
FArray1D<std::complex<double>> AllEncompassingStruct::complex_dp_1d() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_1d_info
  );
}
void AllEncompassingStruct::set_complex_dp_1d(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_complex_dp_1d(fortran_ptr_, v.data(), shape);
}
FArray2D<std::complex<double>> AllEncompassingStruct::complex_dp_2d() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_2d_info
  );
}
void AllEncompassingStruct::set_complex_dp_2d(
    const std::vector<std::vector<std::complex<double>>> &v
) {
  ProxyHelpers::set_array_2d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_set_complex_dp_2d,
      v
  );
}
FArray3D<std::complex<double>> AllEncompassingStruct::complex_dp_3d() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_3d_info
  );
}
void AllEncompassingStruct::set_complex_dp_3d(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &v
) {
  ProxyHelpers::set_array_3d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_set_complex_dp_3d,
      v
  );
}
std::optional<std::complex<double>> AllEncompassingStruct::complex_dp_0d_ptr() const {
  std::complex<double> val;
  bool is_valid;
  all_encompassing_struct_get_complex_dp_0d_ptr(
      fortran_ptr_,
      reinterpret_cast<double _Complex *>(&val),
      &is_valid
  );
  if (is_valid)
    return val;
  return std::nullopt;
}
void AllEncompassingStruct::set_complex_dp_0d_ptr(std::complex<double> value) {
  all_encompassing_struct_set_complex_dp_0d_ptr(fortran_ptr_, value);
}
FArray1D<std::complex<double>> AllEncompassingStruct::complex_dp_1d_ptr() const {
  return ProxyHelpers::get_array_1d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_1d_ptr_info
  );
}
void AllEncompassingStruct::set_complex_dp_1d_ptr(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_complex_dp_1d_ptr(fortran_ptr_, v.data(), shape);
}
FArray2D<std::complex<double>> AllEncompassingStruct::complex_dp_2d_ptr() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_2d_ptr_info
  );
}
void AllEncompassingStruct::set_complex_dp_2d_ptr(
    const std::vector<std::vector<std::complex<double>>> &v
) {
  ProxyHelpers::set_array_2d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_set_complex_dp_2d_ptr,
      v
  );
}
FArray3D<std::complex<double>> AllEncompassingStruct::complex_dp_3d_ptr() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_3d_ptr_info
  );
}
void AllEncompassingStruct::set_complex_dp_3d_ptr(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &v
) {
  ProxyHelpers::set_array_3d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_set_complex_dp_3d_ptr,
      v
  );
}
ComplexAlloc1D AllEncompassingStruct::complex_dp_1d_alloc() const {
  return ComplexAlloc1D(
      const_cast<void *>(fortran_ptr_),
      all_encompassing_struct_reallocate_complex_dp_1d_alloc,
      all_encompassing_struct_get_complex_dp_1d_alloc_info
  );
}
void AllEncompassingStruct::set_complex_dp_1d_alloc(const std::vector<std::complex<double>> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_complex_dp_1d_alloc(fortran_ptr_, v.data(), shape);
}
FArray2D<std::complex<double>> AllEncompassingStruct::complex_dp_2d_alloc() const {
  return ProxyHelpers::get_array_2d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_2d_alloc_info
  );
}
void AllEncompassingStruct::set_complex_dp_2d_alloc(
    const std::vector<std::vector<std::complex<double>>> &v
) {
  ProxyHelpers::set_array_2d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_set_complex_dp_2d_alloc,
      v
  );
}
FArray3D<std::complex<double>> AllEncompassingStruct::complex_dp_3d_alloc() const {
  return ProxyHelpers::get_array_3d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_get_complex_dp_3d_alloc_info
  );
}
void AllEncompassingStruct::set_complex_dp_3d_alloc(
    const std::vector<std::vector<std::vector<std::complex<double>>>> &v
) {
  ProxyHelpers::set_array_3d<std::complex<double>>(
      fortran_ptr_,
      all_encompassing_struct_set_complex_dp_3d_alloc,
      v
  );
}
int AllEncompassingStruct::int_0d() const {
  int value;
  all_encompassing_struct_get_int_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingStruct::set_int_0d(int value) {
  all_encompassing_struct_set_int_0d(fortran_ptr_, value);
}
FArray1D<int> AllEncompassingStruct::int_1d() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, all_encompassing_struct_get_int_1d_info);
}
void AllEncompassingStruct::set_int_1d(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_int_1d(fortran_ptr_, v.data(), shape);
}
FArray2D<int> AllEncompassingStruct::int_2d() const {
  return ProxyHelpers::get_array_2d<int>(fortran_ptr_, all_encompassing_struct_get_int_2d_info);
}
void AllEncompassingStruct::set_int_2d(const std::vector<std::vector<int>> &v) {
  ProxyHelpers::set_array_2d<int>(fortran_ptr_, all_encompassing_struct_set_int_2d, v);
}
FArray3D<int> AllEncompassingStruct::int_3d() const {
  return ProxyHelpers::get_array_3d<int>(fortran_ptr_, all_encompassing_struct_get_int_3d_info);
}
void AllEncompassingStruct::set_int_3d(const std::vector<std::vector<std::vector<int>>> &v) {
  ProxyHelpers::set_array_3d<int>(fortran_ptr_, all_encompassing_struct_set_int_3d, v);
}
std::optional<int> AllEncompassingStruct::int_0d_ptr() const {
  int value;
  bool is_valid;
  all_encompassing_struct_get_int_0d_ptr(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void AllEncompassingStruct::set_int_0d_ptr(int value) {
  all_encompassing_struct_set_int_0d_ptr(fortran_ptr_, value);
}
FArray1D<int> AllEncompassingStruct::int_1d_ptr() const {
  return ProxyHelpers::get_array_1d<int>(fortran_ptr_, all_encompassing_struct_get_int_1d_ptr_info);
}
void AllEncompassingStruct::set_int_1d_ptr(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_int_1d_ptr(fortran_ptr_, v.data(), shape);
}
FArray2D<int> AllEncompassingStruct::int_2d_ptr() const {
  return ProxyHelpers::get_array_2d<int>(fortran_ptr_, all_encompassing_struct_get_int_2d_ptr_info);
}
void AllEncompassingStruct::set_int_2d_ptr(const std::vector<std::vector<int>> &v) {
  ProxyHelpers::set_array_2d<int>(fortran_ptr_, all_encompassing_struct_set_int_2d_ptr, v);
}
FArray3D<int> AllEncompassingStruct::int_3d_ptr() const {
  return ProxyHelpers::get_array_3d<int>(fortran_ptr_, all_encompassing_struct_get_int_3d_ptr_info);
}
void AllEncompassingStruct::set_int_3d_ptr(const std::vector<std::vector<std::vector<int>>> &v) {
  ProxyHelpers::set_array_3d<int>(fortran_ptr_, all_encompassing_struct_set_int_3d_ptr, v);
}
IntAlloc1D AllEncompassingStruct::int_1d_alloc() const {
  return IntAlloc1D(
      const_cast<void *>(fortran_ptr_),
      all_encompassing_struct_reallocate_int_1d_alloc,
      all_encompassing_struct_get_int_1d_alloc_info
  );
}
void AllEncompassingStruct::set_int_1d_alloc(const std::vector<int> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_int_1d_alloc(fortran_ptr_, v.data(), shape);
}
FArray2D<int> AllEncompassingStruct::int_2d_alloc() const {
  return ProxyHelpers::get_array_2d<int>(
      fortran_ptr_,
      all_encompassing_struct_get_int_2d_alloc_info
  );
}
void AllEncompassingStruct::set_int_2d_alloc(const std::vector<std::vector<int>> &v) {
  ProxyHelpers::set_array_2d<int>(fortran_ptr_, all_encompassing_struct_set_int_2d_alloc, v);
}
FArray3D<int> AllEncompassingStruct::int_3d_alloc() const {
  return ProxyHelpers::get_array_3d<int>(
      fortran_ptr_,
      all_encompassing_struct_get_int_3d_alloc_info
  );
}
void AllEncompassingStruct::set_int_3d_alloc(const std::vector<std::vector<std::vector<int>>> &v) {
  ProxyHelpers::set_array_3d<int>(fortran_ptr_, all_encompassing_struct_set_int_3d_alloc, v);
}
int64_t AllEncompassingStruct::int8_0d() const {
  int64_t value;
  all_encompassing_struct_get_int8_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingStruct::set_int8_0d(int64_t value) {
  all_encompassing_struct_set_int8_0d(fortran_ptr_, value);
}
FArray1D<int64_t> AllEncompassingStruct::int8_1d() const {
  return ProxyHelpers::get_array_1d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_1d_info
  );
}
void AllEncompassingStruct::set_int8_1d(const std::vector<int64_t> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_int8_1d(fortran_ptr_, v.data(), shape);
}
FArray2D<int64_t> AllEncompassingStruct::int8_2d() const {
  return ProxyHelpers::get_array_2d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_2d_info
  );
}
void AllEncompassingStruct::set_int8_2d(const std::vector<std::vector<int64_t>> &v) {
  ProxyHelpers::set_array_2d<int64_t>(fortran_ptr_, all_encompassing_struct_set_int8_2d, v);
}
FArray3D<int64_t> AllEncompassingStruct::int8_3d() const {
  return ProxyHelpers::get_array_3d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_3d_info
  );
}
void AllEncompassingStruct::set_int8_3d(const std::vector<std::vector<std::vector<int64_t>>> &v) {
  ProxyHelpers::set_array_3d<int64_t>(fortran_ptr_, all_encompassing_struct_set_int8_3d, v);
}
std::optional<int64_t> AllEncompassingStruct::int8_0d_ptr() const {
  int64_t value;
  bool is_valid;
  all_encompassing_struct_get_int8_0d_ptr(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void AllEncompassingStruct::set_int8_0d_ptr(int64_t value) {
  all_encompassing_struct_set_int8_0d_ptr(fortran_ptr_, value);
}
FArray1D<int64_t> AllEncompassingStruct::int8_1d_ptr() const {
  return ProxyHelpers::get_array_1d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_1d_ptr_info
  );
}
void AllEncompassingStruct::set_int8_1d_ptr(const std::vector<int64_t> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_int8_1d_ptr(fortran_ptr_, v.data(), shape);
}
FArray2D<int64_t> AllEncompassingStruct::int8_2d_ptr() const {
  return ProxyHelpers::get_array_2d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_2d_ptr_info
  );
}
void AllEncompassingStruct::set_int8_2d_ptr(const std::vector<std::vector<int64_t>> &v) {
  ProxyHelpers::set_array_2d<int64_t>(fortran_ptr_, all_encompassing_struct_set_int8_2d_ptr, v);
}
FArray3D<int64_t> AllEncompassingStruct::int8_3d_ptr() const {
  return ProxyHelpers::get_array_3d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_3d_ptr_info
  );
}
void AllEncompassingStruct::set_int8_3d_ptr(const std::vector<std::vector<std::vector<int64_t>>> &v
) {
  ProxyHelpers::set_array_3d<int64_t>(fortran_ptr_, all_encompassing_struct_set_int8_3d_ptr, v);
}
Int8Alloc1D AllEncompassingStruct::int8_1d_alloc() const {
  return Int8Alloc1D(
      const_cast<void *>(fortran_ptr_),
      all_encompassing_struct_reallocate_int8_1d_alloc,
      all_encompassing_struct_get_int8_1d_alloc_info
  );
}
void AllEncompassingStruct::set_int8_1d_alloc(const std::vector<int64_t> &v) {
  int shape[] = {static_cast<int>(v.size())};
  all_encompassing_struct_set_int8_1d_alloc(fortran_ptr_, v.data(), shape);
}
FArray2D<int64_t> AllEncompassingStruct::int8_2d_alloc() const {
  return ProxyHelpers::get_array_2d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_2d_alloc_info
  );
}
void AllEncompassingStruct::set_int8_2d_alloc(const std::vector<std::vector<int64_t>> &v) {
  ProxyHelpers::set_array_2d<int64_t>(fortran_ptr_, all_encompassing_struct_set_int8_2d_alloc, v);
}
FArray3D<int64_t> AllEncompassingStruct::int8_3d_alloc() const {
  return ProxyHelpers::get_array_3d<int64_t>(
      fortran_ptr_,
      all_encompassing_struct_get_int8_3d_alloc_info
  );
}
void AllEncompassingStruct::set_int8_3d_alloc(
    const std::vector<std::vector<std::vector<int64_t>>> &v
) {
  ProxyHelpers::set_array_3d<int64_t>(fortran_ptr_, all_encompassing_struct_set_int8_3d_alloc, v);
}
bool AllEncompassingStruct::logical_0d() const {
  bool value;
  all_encompassing_struct_get_logical_0d(fortran_ptr_, &value);
  return value;
}
void AllEncompassingStruct::set_logical_0d(bool value) {
  all_encompassing_struct_set_logical_0d(fortran_ptr_, value);
}
FArray1D<bool> AllEncompassingStruct::logical_1d() const {
  return ProxyHelpers::get_array_1d<bool>(
      fortran_ptr_,
      all_encompassing_struct_get_logical_1d_info
  );
}
void AllEncompassingStruct::set_logical_1d(const std::vector<bool> &v) {
  int shape[] = {static_cast<int>(v.size())};
  std::vector<int> bv(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    bv[i] = v[i] ? 1 : 0;
  all_encompassing_struct_set_logical_1d(fortran_ptr_, bv.data(), shape);
}
FArray2D<bool> AllEncompassingStruct::logical_2d() const {
  return ProxyHelpers::get_array_2d<bool>(
      fortran_ptr_,
      all_encompassing_struct_get_logical_2d_info
  );
}
void AllEncompassingStruct::set_logical_2d(const std::vector<std::vector<bool>> &v) {
  int rows = static_cast<int>(v.size());
  int cols = rows > 0 ? static_cast<int>(v[0].size()) : 0;
  int shape[] = {cols, rows};

  std::vector<int> flat;
  flat.reserve(rows * cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      flat.push_back(v[i][j] ? 1 : 0);
    }
  }
  all_encompassing_struct_set_logical_2d(fortran_ptr_, flat.data(), shape);
}
FArray3D<bool> AllEncompassingStruct::logical_3d() const {
  return ProxyHelpers::get_array_3d<bool>(
      fortran_ptr_,
      all_encompassing_struct_get_logical_3d_info
  );
}
void AllEncompassingStruct::set_logical_3d(const std::vector<std::vector<std::vector<bool>>> &v) {
  int n3 = static_cast<int>(v.size());
  int n2 = n3 > 0 ? static_cast<int>(v[0].size()) : 0;
  int n1 = n2 > 0 ? static_cast<int>(v[0][0].size()) : 0;
  int shape[] = {n1, n2, n3};

  std::vector<int> flat;
  flat.reserve(n1 * n2 * n3);
  for (int k = 0; k < n3; ++k) {
    for (int j = 0; j < n2; ++j) {
      for (int i = 0; i < n1; ++i) {
        flat.push_back(v[k][j][i] ? 1 : 0);
      }
    }
  }
  all_encompassing_struct_set_logical_3d(fortran_ptr_, flat.data(), shape);
}
std::optional<bool> AllEncompassingStruct::logical_0d_ptr() const {
  bool value;
  bool is_valid;
  all_encompassing_struct_get_logical_0d_ptr(fortran_ptr_, &value, &is_valid);
  if (is_valid)
    return value;
  return std::nullopt;
}
void AllEncompassingStruct::set_logical_0d_ptr(bool value) {
  all_encompassing_struct_set_logical_0d_ptr(fortran_ptr_, value);
}
FArray1D<bool> AllEncompassingStruct::logical_1d_ptr() const {
  return ProxyHelpers::get_array_1d<bool>(
      fortran_ptr_,
      all_encompassing_struct_get_logical_1d_ptr_info
  );
}
void AllEncompassingStruct::set_logical_1d_ptr(const std::vector<bool> &v) {
  int shape[] = {static_cast<int>(v.size())};
  std::vector<int> bv(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    bv[i] = v[i] ? 1 : 0;
  all_encompassing_struct_set_logical_1d_ptr(fortran_ptr_, bv.data(), shape);
}
TestSubStruct AllEncompassingStruct::type_0d() const {
  void *ptr;
  all_encompassing_struct_get_type_0d(fortran_ptr_, &ptr);
  return TestSubStruct(ptr);
}
void AllEncompassingStruct::set_type_0d(const TestSubStruct &src) {
  all_encompassing_struct_set_type_0d(fortran_ptr_, src.get_fortran_ptr());
}
TestSubStructArray1D AllEncompassingStruct::type_1d() const {
  return ProxyHelpers::get_type_array_1d<TestSubStructArray1D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_1d_info
  );
}
TestSubStructArray2D AllEncompassingStruct::type_2d() const {
  return ProxyHelpers::get_type_array_2d<TestSubStructArray2D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_2d_info
  );
}
TestSubStructArray3D AllEncompassingStruct::type_3d() const {
  return ProxyHelpers::get_type_array_3d<TestSubStructArray3D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_3d_info
  );
}
std::optional<TestSubStruct> AllEncompassingStruct::type_0d_ptr() const {
  void *ptr;
  all_encompassing_struct_get_type_0d_ptr(fortran_ptr_, &ptr);
  if (!ptr)
    return std::nullopt;
  return TestSubStruct(ptr);
}
void AllEncompassingStruct::set_type_0d_ptr(const TestSubStruct &src) {
  all_encompassing_struct_set_type_0d_ptr(fortran_ptr_, src.get_fortran_ptr());
}
TestSubStructArray1D AllEncompassingStruct::type_1d_ptr() const {
  return ProxyHelpers::get_type_array_1d<TestSubStructArray1D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_1d_ptr_info
  );
}
TestSubStructArray2D AllEncompassingStruct::type_2d_ptr() const {
  return ProxyHelpers::get_type_array_2d<TestSubStructArray2D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_2d_ptr_info
  );
}
TestSubStructArray3D AllEncompassingStruct::type_3d_ptr() const {
  return ProxyHelpers::get_type_array_3d<TestSubStructArray3D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_3d_ptr_info
  );
}
TestSubStructAlloc1D AllEncompassingStruct::type_1d_alloc() const {
  return TestSubStructAlloc1D(
      const_cast<void *>(fortran_ptr_),
      all_encompassing_struct_reallocate_type_1d_alloc,
      all_encompassing_struct_get_type_1d_alloc_info
  );
}
TestSubStructArray2D AllEncompassingStruct::type_2d_alloc() const {
  return ProxyHelpers::get_type_array_2d<TestSubStructArray2D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_2d_alloc_info
  );
}
TestSubStructArray3D AllEncompassingStruct::type_3d_alloc() const {
  return ProxyHelpers::get_type_array_3d<TestSubStructArray3D>(
      fortran_ptr_,
      all_encompassing_struct_get_type_3d_alloc_info
  );
}
TestSubSubStruct TestSubStruct::sr() const {
  void *ptr;
  test_sub_struct_get_sr(fortran_ptr_, &ptr);
  return TestSubSubStruct(ptr);
}
void TestSubStruct::set_sr(const TestSubSubStruct &src) {
  test_sub_struct_set_sr(fortran_ptr_, src.get_fortran_ptr());
}
int64_t TestSubSubStruct::aaa() const {
  int64_t value;
  test_sub_sub_struct_get_aaa(fortran_ptr_, &value);
  return value;
}
void TestSubSubStruct::set_aaa(int64_t value) { test_sub_sub_struct_set_aaa(fortran_ptr_, value); }
int TestSubSubStruct::bbb() const {
  int value;
  test_sub_sub_struct_get_bbb(fortran_ptr_, &value);
  return value;
}
void TestSubSubStruct::set_bbb(int value) { test_sub_sub_struct_set_bbb(fortran_ptr_, value); }
std::string TestSubSubStruct::file() const {
  FArray1D<char> arr =
      ProxyHelpers::get_array_1d<char>(fortran_ptr_, test_sub_sub_struct_get_file_info);
  return std::string(arr.data(), arr.size());
}
void TestSubSubStruct::set_file(const std::string &value) {
  test_sub_sub_struct_set_file(fortran_ptr_, value.c_str(), static_cast<int>(value.length()));
}
double TestSubSubStruct::t_ref() const {
  double value;
  test_sub_sub_struct_get_t_ref(fortran_ptr_, &value);
  return value;
}
void TestSubSubStruct::set_t_ref(double value) {
  test_sub_sub_struct_set_t_ref(fortran_ptr_, value);
}
double TestSubSubStruct::freq_spread() const {
  double value;
  test_sub_sub_struct_get_freq_spread(fortran_ptr_, &value);
  return value;
}
void TestSubSubStruct::set_freq_spread(double value) {
  test_sub_sub_struct_set_freq_spread(fortran_ptr_, value);
}