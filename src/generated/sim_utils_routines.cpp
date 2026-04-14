#include "bmad/generated/sim_utils_routines.hpp"

#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/proxy.hpp"
#include "bmad/json.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

using json = nlohmann::json;
void SimUtils::allocate_thread_states() { fortran_allocate_thread_states(); }
double SimUtils::anomalous_moment_of(int species) {
  double _moment{};
  fortran_anomalous_moment_of(/* int& */ species, /* double& */ _moment);
  return _moment;
}
int SimUtils::antiparticle(int species) {
  int _anti_species{};
  fortran_antiparticle(/* int& */ species, /* int& */ _anti_species);
  return _anti_species;
}
void SimUtils::apfft(
    FArray1D<Real> &rdata_in,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    std::optional<int> diag
) {
  // rdata_in: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _rdata_in_desc;
  _rdata_in_desc.rank = 1;
  _rdata_in_desc.data_ptr = rdata_in.data();
  _rdata_in_desc.dims[0] = rdata_in.size();
  // bounds: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _bounds_desc;
  _bounds_desc.rank = 1;
  _bounds_desc.data_ptr = bounds.data();
  _bounds_desc.dims[0] = bounds.size();
  auto _window = window.c_str();
  int diag_lvalue;
  auto *_diag{&diag_lvalue};
  if (diag.has_value()) {
    diag_lvalue = diag.value();
  } else {
    _diag = nullptr;
  }
  fortran_apfft(
      /* Bmad::array_descriptor_t& */ _rdata_in_desc,
      /* Bmad::array_descriptor_t& */ _bounds_desc,
      /* const char* */ _window,
      /* double& */ phase,
      /* int* */ _diag
  );
}
SimUtils::ApfftCorr SimUtils::apfft_corr(
    FArray1D<Real> &rdata_in,
    std::string window,
    std::optional<FixedArray1D<Real, 2>> bounds,
    std::optional<int> diag
) {
  // rdata_in: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _rdata_in_desc;
  _rdata_in_desc.rank = 1;
  _rdata_in_desc.data_ptr = rdata_in.data();
  _rdata_in_desc.dims[0] = rdata_in.size();
  // bounds: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _bounds_desc;
  _bounds_desc.rank = 1;
  if (bounds.has_value()) {
    _bounds_desc.data_ptr = bounds->data();
    _bounds_desc.dims[0] = bounds->size();
  } else {
    _bounds_desc.data_ptr = nullptr;
    _bounds_desc.dims[0] = 0;
  }
  auto _window = window.c_str();
  double _phase{};
  double _amp{};
  double _freq{};
  int diag_lvalue;
  auto *_diag{&diag_lvalue};
  if (diag.has_value()) {
    diag_lvalue = diag.value();
  } else {
    _diag = nullptr;
  }
  fortran_apfft_corr(
      /* Bmad::array_descriptor_t& */ _rdata_in_desc,
      /* Bmad::array_descriptor_t& */ _bounds_desc,
      /* const char* */ _window,
      /* double& */ _phase,
      /* double& */ _amp,
      /* double& */ _freq,
      /* int* */ _diag
  );
  return ApfftCorr{_phase, _amp, _freq};
}
void SimUtils::apfft_ext(
    FArray1D<Real> &rdata,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    double amp,
    double freq,
    std::optional<int> diag
) {
  // rdata: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _rdata_desc;
  _rdata_desc.rank = 1;
  _rdata_desc.data_ptr = rdata.data();
  _rdata_desc.dims[0] = rdata.size();
  // bounds: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _bounds_desc;
  _bounds_desc.rank = 1;
  _bounds_desc.data_ptr = bounds.data();
  _bounds_desc.dims[0] = bounds.size();
  auto _window = window.c_str();
  int diag_lvalue;
  auto *_diag{&diag_lvalue};
  if (diag.has_value()) {
    diag_lvalue = diag.value();
  } else {
    _diag = nullptr;
  }
  fortran_apfft_ext(
      /* Bmad::array_descriptor_t& */ _rdata_desc,
      /* Bmad::array_descriptor_t& */ _bounds_desc,
      /* const char* */ _window,
      /* double& */ phase,
      /* double& */ amp,
      /* double& */ freq,
      /* int* */ _diag
  );
}
double SimUtils::asinc(double x, std::optional<int> nd) {
  int nd_lvalue;
  auto *_nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  double _y{};
  fortran_asinc(/* double& */ x, /* int* */ _nd, /* double& */ _y);
  return _y;
}
void SimUtils::assert_equal(FArray1D<Int> &int_arr, std::string err_str, int ival) {
  // int_arr: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _int_arr_desc;
  _int_arr_desc.rank = 1;
  _int_arr_desc.data_ptr = int_arr.data();
  _int_arr_desc.dims[0] = int_arr.size();
  auto _err_str = err_str.c_str();
  fortran_assert_equal(
      /* Bmad::array_descriptor_t& */ _int_arr_desc,
      /* const char* */ _err_str,
      /* int& */ ival
  );
}
int SimUtils::atomic_number(int species) {
  int _atomic_num{};
  fortran_atomic_number(/* int& */ species, /* int& */ _atomic_num);
  return _atomic_num;
}
int SimUtils::atomic_species_id(int charge, bool is_anti, int atomic_num, int n_nuc) {
  int _species_id{};
  fortran_atomic_species_id(
      /* int& */ charge,
      /* bool& */ is_anti,
      /* int& */ atomic_num,
      /* int& */ n_nuc,
      /* int& */ _species_id
  );
  return _species_id;
}
FixedArray1D<Real, 4> SimUtils::axis_angle_to_quat(FixedArray1D<Real, 3> axis, double angle) {
  // axis: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _axis_desc;
  _axis_desc.rank = 1;
  _axis_desc.data_ptr = axis.data();
  _axis_desc.dims[0] = axis.size();
  // quat: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  FixedArray1D<Real, 4> _quat;
  _quat_desc.data_ptr = _quat.data();
  _quat_desc.dims[0] = _quat.size();
  fortran_axis_angle_to_quat(
      /* Bmad::array_descriptor_t& */ _axis_desc,
      /* double& */ angle,
      /* Bmad::array_descriptor_t& */ _quat_desc
  );
  return _quat;
}
FixedArray2D<Real, 3, 3> SimUtils::axis_angle_to_w_mat(FixedArray1D<Real, 3> axis, double angle) {
  // axis: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _axis_desc;
  _axis_desc.rank = 1;
  _axis_desc.data_ptr = axis.data();
  _axis_desc.dims[0] = axis.size();
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  fortran_axis_angle_to_w_mat(
      /* Bmad::array_descriptor_t& */ _axis_desc,
      /* double& */ angle,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
SimUtils::BicubicCmplxEval
SimUtils::bicubic_cmplx_eval(double x_norm, double y_norm, BicubicCmplxCoefStruct &bi_coef) {
  std::complex<double> _df_dx{};
  std::complex<double> _df_dy{};
  std::complex<double> _f_val{};
  fortran_bicubic_cmplx_eval(
      /* double& */ x_norm,
      /* double& */ y_norm,
      /* void* */ bi_coef.get_fortran_ptr(),
      /* std::complex<double>& */ _df_dx,
      /* std::complex<double>& */ _df_dy,
      /* std::complex<double>& */ _f_val
  );
  return BicubicCmplxEval{_df_dx, _df_dy, _f_val};
}
int SimUtils::bin_index(double x, double bin1_x_min, double bin_delta) {
  int _ix_bin{};
  fortran_bin_index(
      /* double& */ x,
      /* double& */ bin1_x_min,
      /* double& */ bin_delta,
      /* int& */ _ix_bin
  );
  return _ix_bin;
}
void SimUtils::bin_x_center(int &ix_bin, double bin1_x_min, double bin_delta, double x_center) {
  fortran_bin_x_center(
      /* int& */ ix_bin,
      /* double& */ bin1_x_min,
      /* double& */ bin_delta,
      /* double& */ x_center
  );
}
void SimUtils::bit_set(int &word, int pos, bool set_to_1) {
  fortran_bit_set(/* int& */ word, /* int& */ pos, /* bool& */ set_to_1);
}
SimUtils::BracketIndexForSpline SimUtils::bracket_index_for_spline(
    FArray1D<Real> &x_knot,
    double x,
    std::optional<bool> strict,
    std::optional<bool> print_err
) {
  // x_knot: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_knot_desc;
  _x_knot_desc.rank = 1;
  _x_knot_desc.data_ptr = x_knot.data();
  _x_knot_desc.dims[0] = x_knot.size();
  int _ix0{};
  bool strict_lvalue;
  auto *_strict{&strict_lvalue};
  if (strict.has_value()) {
    strict_lvalue = strict.value();
  } else {
    _strict = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool _ok{};
  fortran_bracket_index_for_spline(
      /* Bmad::array_descriptor_t& */ _x_knot_desc,
      /* double& */ x,
      /* int& */ _ix0,
      /* bool* */ _strict,
      /* bool* */ _print_err,
      /* bool& */ _ok
  );
  return BracketIndexForSpline{_ix0, _ok};
}
void SimUtils::calc_file_number(std::string file_name, int num_in, int num_out, bool err_flag) {
  auto _file_name = file_name.c_str();
  fortran_calc_file_number(
      /* const char* */ _file_name,
      /* int& */ num_in,
      /* int& */ num_out,
      /* bool& */ err_flag
  );
}
SimUtils::Celbd SimUtils::celbd(double mc) {
  double _elb{};
  double _eld{};
  fortran_celbd(/* double& */ mc, /* double& */ _elb, /* double& */ _eld);
  return Celbd{_elb, _eld};
}
std::string SimUtils::cesr_getarg(int i_arg) {
  char _arg[4096];
  fortran_cesr_getarg(/* int& */ i_arg, /* const char* */ _arg);
  return _arg;
}
void SimUtils::cesr_iargc(int func_retval__) { fortran_cesr_iargc(/* int& */ func_retval__); }
void SimUtils::change_file_number(std::string file_name, int change) {
  auto _file_name = file_name.c_str();
  fortran_change_file_number(/* const char* */ _file_name, /* int& */ change);
}
int SimUtils::charge_of(int species, std::optional<int> default_) {
  int default__lvalue;
  auto *_default_{&default__lvalue};
  if (default_.has_value()) {
    default__lvalue = default_.value();
  } else {
    _default_ = nullptr;
  }
  int _charge{};
  fortran_charge_of(/* int& */ species, /* int* */ _default_, /* int& */ _charge);
  return _charge;
}
double SimUtils::charge_to_mass_of(int species) {
  double _charge_mass_ratio{};
  fortran_charge_to_mass_of(/* int& */ species, /* double& */ _charge_mass_ratio);
  return _charge_mass_ratio;
}
double SimUtils::coarse_frequency_estimate(FArray1D<Real> &data, std::optional<bool> error) {
  // data: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _data_desc;
  _data_desc.rank = 1;
  _data_desc.data_ptr = data.data();
  _data_desc.dims[0] = data.size();
  bool error_lvalue;
  auto *_error{&error_lvalue};
  if (error.has_value()) {
    error_lvalue = error.value();
  } else {
    _error = nullptr;
  }
  double _frequency{};
  fortran_coarse_frequency_estimate(
      /* Bmad::array_descriptor_t& */ _data_desc,
      /* bool* */ _error,
      /* double& */ _frequency
  );
  return _frequency;
}
void SimUtils::complex_error_function(double wr, double wi, double zr, double zi) {
  fortran_complex_error_function(
      /* double& */ wr,
      /* double& */ wi,
      /* double& */ zr,
      /* double& */ zi
  );
}
double SimUtils::cos_one(double angle) {
  double _cos1{};
  fortran_cos_one(/* double& */ angle, /* double& */ _cos1);
  return _cos1;
}
double SimUtils::cosc(double x, std::optional<int> nd) {
  int nd_lvalue;
  auto *_nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  double _y{};
  fortran_cosc(/* double& */ x, /* int* */ _nd, /* double& */ _y);
  return _y;
}
SplineStruct
SimUtils::create_a_spline(FArray1D<Real> &r0, FArray1D<Real> &r1, double slope0, double slope1) {
  // r0: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _r0_desc;
  _r0_desc.rank = 1;
  _r0_desc.data_ptr = r0.data();
  _r0_desc.dims[0] = r0.size();
  // r1: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _r1_desc;
  _r1_desc.rank = 1;
  _r1_desc.data_ptr = r1.data();
  _r1_desc.dims[0] = r1.size();
  SplineStruct _spline;
  fortran_create_a_spline(
      /* Bmad::array_descriptor_t& */ _r0_desc,
      /* Bmad::array_descriptor_t& */ _r1_desc,
      /* double& */ slope0,
      /* double& */ slope1,
      /* void* */ _spline.get_fortran_ptr()
  );
  return std::move(_spline);
}
FixedArray1D<Real, 3> SimUtils::cross_product(FArray1D<Real> &a, FArray1D<Real> &b) {
  // a: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _a_desc;
  _a_desc.rank = 1;
  _a_desc.data_ptr = a.data();
  _a_desc.dims[0] = a.size();
  // b: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _b_desc;
  _b_desc.rank = 1;
  _b_desc.data_ptr = b.data();
  _b_desc.dims[0] = b.size();
  // c: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _c_desc;
  _c_desc.rank = 1;
  FixedArray1D<Real, 3> _c;
  _c_desc.data_ptr = _c.data();
  _c_desc.dims[0] = _c.size();
  fortran_cross_product(
      /* Bmad::array_descriptor_t& */ _a_desc,
      /* Bmad::array_descriptor_t& */ _b_desc,
      /* Bmad::array_descriptor_t& */ _c_desc
  );
  return _c;
}
void SimUtils::date_and_time_stamp(
    std::string string,
    std::optional<bool> numeric_month,
    std::optional<bool> include_zone
) {
  auto _string = string.c_str();
  bool numeric_month_lvalue;
  auto *_numeric_month{&numeric_month_lvalue};
  if (numeric_month.has_value()) {
    numeric_month_lvalue = numeric_month.value();
  } else {
    _numeric_month = nullptr;
  }
  bool include_zone_lvalue;
  auto *_include_zone{&include_zone_lvalue};
  if (include_zone.has_value()) {
    include_zone_lvalue = include_zone.value();
  } else {
    _include_zone = nullptr;
  }
  fortran_date_and_time_stamp(
      /* const char* */ _string,
      /* bool* */ _numeric_month,
      /* bool* */ _include_zone
  );
}
void SimUtils::destfixedwindowls(int id) { fortran_destfixedwindowls(/* int& */ id); }
void SimUtils::detab(std::string str) {
  auto _str = str.c_str();
  fortran_detab(/* const char* */ _str);
}
void SimUtils::display_size_and_resolution(
    int ix_screen,
    double x_size,
    double y_size,
    double x_res,
    double y_res
) {
  fortran_display_size_and_resolution(
      /* int& */ ix_screen,
      /* double& */ x_size,
      /* double& */ y_size,
      /* double& */ x_res,
      /* double& */ y_res
  );
}
double SimUtils::dj_bessel(int m, double arg) {
  double _dj_bes{};
  fortran_dj_bessel(/* int& */ m, /* double& */ arg, /* double& */ _dj_bes);
  return _dj_bes;
}
void SimUtils::djb_hash(std::string str, int hash, std::optional<int> old_hash) {
  auto _str = str.c_str();
  int old_hash_lvalue;
  auto *_old_hash{&old_hash_lvalue};
  if (old_hash.has_value()) {
    old_hash_lvalue = old_hash.value();
  } else {
    _old_hash = nullptr;
  }
  fortran_djb_hash(/* const char* */ _str, /* int* */ _old_hash, /* int& */ hash);
}
void SimUtils::djb_str_hash(std::string in_str, std::string hash_str) {
  auto _in_str = in_str.c_str();
  auto _hash_str = hash_str.c_str();
  fortran_djb_str_hash(/* const char* */ _in_str, /* const char* */ _hash_str);
}
void SimUtils::downcase_string(std::string string) {
  auto _string = string.c_str();
  fortran_downcase_string(/* const char* */ _string);
}
SimUtils::Elbd SimUtils::elbd(double phi, double phic, double mc) {
  double _b{};
  double _d{};
  fortran_elbd(
      /* double& */ phi,
      /* double& */ phic,
      /* double& */ mc,
      /* double& */ _b,
      /* double& */ _d
  );
  return Elbd{_b, _d};
}
SimUtils::Elcbd SimUtils::elcbd(double c0, double mc) {
  double _b{};
  double _dx{};
  fortran_elcbd(/* double& */ c0, /* double& */ mc, /* double& */ _b, /* double& */ _dx);
  return Elcbd{_b, _dx};
}
SimUtils::Ellipinc SimUtils::ellipinc(double phi, double m) {
  double _ellipkinc{};
  double _ellipeinc{};
  fortran_ellipinc(
      /* double& */ phi,
      /* double& */ m,
      /* double& */ _ellipkinc,
      /* double& */ _ellipeinc
  );
  return Ellipinc{_ellipkinc, _ellipeinc};
}
SimUtils::Elsbd SimUtils::elsbd(double s0, double mc) {
  double _b{};
  double _d{};
  fortran_elsbd(/* double& */ s0, /* double& */ mc, /* double& */ _b, /* double& */ _d);
  return Elsbd{_b, _d};
}
void SimUtils::end_akima_spline_calc(SplineStructArray1D spline, int which_end) {
  // spline: SplineStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spline_desc;
  _spline_desc.rank = 1;
  _spline_desc.data_ptr = spline.data();
  _spline_desc.dims[0] = spline.size();
  _spline_desc.strides[0] = 1;
  fortran_end_akima_spline_calc(/* Bmad::array_descriptor_t& */ _spline_desc, /* int& */ which_end);
}
void SimUtils::err_exit(std::optional<std::string> err_str) {
  const char *_err_str = err_str.has_value() ? err_str->c_str() : nullptr;
  fortran_err_exit(/* const char* */ _err_str);
}
double SimUtils::factorial(int n) {
  double _fact{};
  fortran_factorial(/* int& */ n, /* double& */ _fact);
  return _fact;
}
void SimUtils::faddeeva_function(
    FixedArray1D<Real, 2> z,
    FixedArray1D<Real, 2> w,
    FixedArray2D<Real, 2, 2> dw
) {
  // z: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _z_desc;
  _z_desc.rank = 1;
  _z_desc.data_ptr = z.data();
  _z_desc.dims[0] = z.size();
  // w: inout NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _w_desc;
  _w_desc.rank = 1;
  _w_desc.data_ptr = w.data();
  _w_desc.dims[0] = w.size();
  // dw: inout NOT (CppWrapperGeneralArgumentArray) (['2', '2'])
  Bmad::array_descriptor_t _dw_desc;
  _dw_desc.rank = 2;
  double _dw_vec[2 * 2];
  _dw_desc.data_ptr = _dw_vec;
  _dw_desc.dims[0] = 2;
  _dw_desc.dims[1] = 2;
  matrix_to_vec(dw, _dw_vec);
  fortran_faddeeva_function(
      /* Bmad::array_descriptor_t& */ _z_desc,
      /* Bmad::array_descriptor_t& */ _w_desc,
      /* Bmad::array_descriptor_t& */ _dw_desc
  );
  vec_to_matrix(_dw_vec, dw);
}
void SimUtils::fft_1d(FArray1D<Complex> &arr, int isign) {
  // arr: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_desc;
  _arr_desc.rank = 1;
  _arr_desc.data_ptr = arr.data();
  _arr_desc.dims[0] = arr.size();
  fortran_fft_1d(/* Bmad::array_descriptor_t& */ _arr_desc, /* int& */ isign);
}
void SimUtils::file_directorizer(
    std::string in_file,
    std::string out_file,
    std::string directory,
    bool add_switch
) {
  auto _in_file = in_file.c_str();
  auto _out_file = out_file.c_str();
  auto _directory = directory.c_str();
  fortran_file_directorizer(
      /* const char* */ _in_file,
      /* const char* */ _out_file,
      /* const char* */ _directory,
      /* bool& */ add_switch
  );
}
void SimUtils::file_get(std::string string, std::string dflt_file_name, std::string file_name) {
  auto _string = string.c_str();
  auto _dflt_file_name = dflt_file_name.c_str();
  auto _file_name = file_name.c_str();
  fortran_file_get(
      /* const char* */ _string,
      /* const char* */ _dflt_file_name,
      /* const char* */ _file_name
  );
}
void SimUtils::file_get_open(
    std::string string,
    std::string dflt_file_name,
    std::string file_name,
    int file_unit,
    bool readonly
) {
  auto _string = string.c_str();
  auto _dflt_file_name = dflt_file_name.c_str();
  auto _file_name = file_name.c_str();
  fortran_file_get_open(
      /* const char* */ _string,
      /* const char* */ _dflt_file_name,
      /* const char* */ _file_name,
      /* int& */ file_unit,
      /* bool& */ readonly
  );
}
void SimUtils::file_suffixer(
    std::string in_file_name,
    std::string out_file_name,
    std::string suffix,
    bool add_switch
) {
  auto _in_file_name = in_file_name.c_str();
  auto _out_file_name = out_file_name.c_str();
  auto _suffix = suffix.c_str();
  fortran_file_suffixer(
      /* const char* */ _in_file_name,
      /* const char* */ _out_file_name,
      /* const char* */ _suffix,
      /* bool& */ add_switch
  );
}
void SimUtils::find_location(FArray1D<Int> &arr, int value, int ix_match) {
  // arr: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_desc;
  _arr_desc.rank = 1;
  _arr_desc.data_ptr = arr.data();
  _arr_desc.dims[0] = arr.size();
  fortran_find_location_int(
      /* Bmad::array_descriptor_t& */ _arr_desc,
      /* int& */ value,
      /* int& */ ix_match
  );
}
void SimUtils::find_location(BoolAlloc1D &arr, bool value, int ix_match) {
  // intent=inout allocatable general array
  fortran_find_location_logic(
      /* void* */ arr.get_fortran_ptr(),
      /* bool& */ value,
      /* int& */ ix_match
  );
}
int SimUtils::find_location(FArray1D<Real> &arr, double value) {
  // arr: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_desc;
  _arr_desc.rank = 1;
  _arr_desc.data_ptr = arr.data();
  _arr_desc.dims[0] = arr.size();
  int _ix_match{};
  fortran_find_location_real(
      /* Bmad::array_descriptor_t& */ _arr_desc,
      /* double& */ value,
      /* int& */ _ix_match
  );
  return _ix_match;
}
void SimUtils::find_location(CharacterAlloc1D &arr, std::string value, int ix_match) {
  // intent=inout character array container
  auto _value = value.c_str();
  fortran_find_location_str(
      /* void* */ arr.get_fortran_ptr(),
      /* const char* */ _value,
      /* int& */ ix_match
  );
}
double SimUtils::fine_frequency_estimate(FArray1D<Real> &data) {
  // data: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _data_desc;
  _data_desc.rank = 1;
  _data_desc.data_ptr = data.data();
  _data_desc.dims[0] = data.size();
  double _frequency{};
  fortran_fine_frequency_estimate(
      /* Bmad::array_descriptor_t& */ _data_desc,
      /* double& */ _frequency
  );
  return _frequency;
}
void SimUtils::fixedwindowls(double ynew, int id, double z) {
  fortran_fixedwindowls(/* double& */ ynew, /* int& */ id, /* double& */ z);
}
SimUtils::FourierAmplitude SimUtils::fourier_amplitude(FArray1D<Real> &data, double frequency) {
  // data: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _data_desc;
  _data_desc.rank = 1;
  _data_desc.data_ptr = data.data();
  _data_desc.dims[0] = data.size();
  double _cos_amp{};
  double _sin_amp{};
  double _dcos_amp{};
  double _dsin_amp{};
  fortran_fourier_amplitude(
      /* Bmad::array_descriptor_t& */ _data_desc,
      /* double& */ frequency,
      /* double& */ _cos_amp,
      /* double& */ _sin_amp,
      /* double& */ _dcos_amp,
      /* double& */ _dsin_amp
  );
  return FourierAmplitude{_cos_amp, _sin_amp, _dcos_amp, _dsin_amp};
}
SimUtils::Gelbd SimUtils::gelbd(double phi, double mc) {
  double _elb{};
  double _eld{};
  fortran_gelbd(/* double& */ phi, /* double& */ mc, /* double& */ _elb, /* double& */ _eld);
  return Gelbd{_elb, _eld};
}
double SimUtils::gen_complete_elliptic(
    double kc,
    double p,
    double c,
    double s,
    std::optional<double> err_tol
) {
  double err_tol_lvalue;
  auto *_err_tol{&err_tol_lvalue};
  if (err_tol.has_value()) {
    err_tol_lvalue = err_tol.value();
  } else {
    _err_tol = nullptr;
  }
  double _value{};
  fortran_gen_complete_elliptic(
      /* double& */ kc,
      /* double& */ p,
      /* double& */ c,
      /* double& */ s,
      /* double* */ _err_tol,
      /* double& */ _value
  );
  return _value;
}
std::string SimUtils::get_a_char(bool wait, optional_ref<CharacterAlloc1D> ignore_this) {
  char _this_char[4096];
  // intent=in character array container
  auto *_ignore_this =
      ignore_this.has_value() ? ignore_this->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_get_a_char(/* const char* */ _this_char, /* bool& */ wait, /* void* */ _ignore_this);
  return _this_char;
}
void SimUtils::get_file_number(
    std::string file_name,
    std::string cnum_in,
    int num_out,
    bool err_flag
) {
  auto _file_name = file_name.c_str();
  auto _cnum_in = cnum_in.c_str();
  fortran_get_file_number(
      /* const char* */ _file_name,
      /* const char* */ _cnum_in,
      /* int& */ num_out,
      /* bool& */ err_flag
  );
}
void SimUtils::get_file_time_stamp(std::string file, std::string time_stamp) {
  auto _file = file.c_str();
  auto _time_stamp = time_stamp.c_str();
  fortran_get_file_time_stamp(/* const char* */ _file, /* const char* */ _time_stamp);
}
std::string SimUtils::get_tty_char(bool wait, bool flush) {
  char _this_char[4096];
  fortran_get_tty_char(/* const char* */ _this_char, /* bool& */ wait, /* bool& */ flush);
  return _this_char;
}
void SimUtils::hanhan(int N, FArray1D<Real> &hh) {
  // hh: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _hh_desc;
  _hh_desc.rank = 1;
  _hh_desc.data_ptr = hh.data();
  _hh_desc.dims[0] = hh.size();
  fortran_hanhan(/* int& */ N, /* Bmad::array_descriptor_t& */ _hh_desc);
}
double SimUtils::i_bessel(int m, double arg) {
  double _i_bes{};
  fortran_i_bessel(/* int& */ m, /* double& */ arg, /* double& */ _i_bes);
  return _i_bes;
}
std::complex<double> SimUtils::i_bessel_extended(int m, double arg) {
  std::complex<double> _i_bes{};
  fortran_i_bessel_extended(/* int& */ m, /* double& */ arg, /* std::complex<double>& */ _i_bes);
  return _i_bes;
}
void SimUtils::increment_file_number(
    std::string file_name,
    int digits,
    int number,
    std::string cnumber
) {
  auto _file_name = file_name.c_str();
  auto _cnumber = cnumber.c_str();
  fortran_increment_file_number(
      /* const char* */ _file_name,
      /* int& */ digits,
      /* int& */ number,
      /* const char* */ _cnumber
  );
}
void SimUtils::index_nocase(std::string string1, std::string string2, int indx) {
  auto _string1 = string1.c_str();
  auto _string2 = string2.c_str();
  fortran_index_nocase(/* const char* */ _string1, /* const char* */ _string2, /* int& */ indx);
}
void SimUtils::initfixedwindowls(int N, double dt, int order, int der, int id) {
  fortran_initfixedwindowls(
      /* int& */ N,
      /* double& */ dt,
      /* int& */ order,
      /* int& */ der,
      /* int& */ id
  );
}
void SimUtils::initial_lmdif() { fortran_initial_lmdif(); }
void SimUtils::int_str(int int_, std::string str, std::optional<int> width) {
  int width_lvalue;
  auto *_width{&width_lvalue};
  if (width.has_value()) {
    width_lvalue = width.value();
  } else {
    _width = nullptr;
  }
  auto _str = str.c_str();
  fortran_int_str(/* int& */ int_, /* int* */ _width, /* const char* */ _str);
}
void SimUtils::interpolated_fft(
    FArray1D<Complex> &cdata,
    bool calc_ok,
    double this_fft,
    std::optional<int> opt_dump_spectrum,
    std::optional<int> opt_dump_index
) {
  // cdata: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _cdata_desc;
  _cdata_desc.rank = 1;
  _cdata_desc.data_ptr = cdata.data();
  _cdata_desc.dims[0] = cdata.size();
  int opt_dump_spectrum_lvalue;
  auto *_opt_dump_spectrum{&opt_dump_spectrum_lvalue};
  if (opt_dump_spectrum.has_value()) {
    opt_dump_spectrum_lvalue = opt_dump_spectrum.value();
  } else {
    _opt_dump_spectrum = nullptr;
  }
  int opt_dump_index_lvalue;
  auto *_opt_dump_index{&opt_dump_index_lvalue};
  if (opt_dump_index.has_value()) {
    opt_dump_index_lvalue = opt_dump_index.value();
  } else {
    _opt_dump_index = nullptr;
  }
  fortran_interpolated_fft(
      /* Bmad::array_descriptor_t& */ _cdata_desc,
      /* bool& */ calc_ok,
      /* int* */ _opt_dump_spectrum,
      /* int* */ _opt_dump_index,
      /* double& */ this_fft
  );
}
void SimUtils::interpolated_fft_gsl(
    FArray1D<Complex> &cdata,
    bool calc_ok,
    double this_fft,
    std::optional<int> opt_dump_spectrum,
    std::optional<int> opt_dump_index
) {
  // cdata: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _cdata_desc;
  _cdata_desc.rank = 1;
  _cdata_desc.data_ptr = cdata.data();
  _cdata_desc.dims[0] = cdata.size();
  int opt_dump_spectrum_lvalue;
  auto *_opt_dump_spectrum{&opt_dump_spectrum_lvalue};
  if (opt_dump_spectrum.has_value()) {
    opt_dump_spectrum_lvalue = opt_dump_spectrum.value();
  } else {
    _opt_dump_spectrum = nullptr;
  }
  int opt_dump_index_lvalue;
  auto *_opt_dump_index{&opt_dump_index_lvalue};
  if (opt_dump_index.has_value()) {
    opt_dump_index_lvalue = opt_dump_index.value();
  } else {
    _opt_dump_index = nullptr;
  }
  fortran_interpolated_fft_gsl(
      /* Bmad::array_descriptor_t& */ _cdata_desc,
      /* bool& */ calc_ok,
      /* int* */ _opt_dump_spectrum,
      /* int* */ _opt_dump_index,
      /* double& */ this_fft
  );
}
void SimUtils::is_alphabetic(
    std::string string,
    bool is_alpha,
    std::optional<std::string> valid_chars
) {
  auto _string = string.c_str();
  const char *_valid_chars = valid_chars.has_value() ? valid_chars->c_str() : nullptr;
  fortran_is_alphabetic(
      /* const char* */ _string,
      /* const char* */ _valid_chars,
      /* bool& */ is_alpha
  );
}
bool SimUtils::is_decreasing_sequence(FArray1D<Real> &array, std::optional<bool> strict) {
  // array: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _array_desc;
  _array_desc.rank = 1;
  _array_desc.data_ptr = array.data();
  _array_desc.dims[0] = array.size();
  bool strict_lvalue;
  auto *_strict{&strict_lvalue};
  if (strict.has_value()) {
    strict_lvalue = strict.value();
  } else {
    _strict = nullptr;
  }
  bool _is_decreasing{};
  fortran_is_decreasing_sequence(
      /* Bmad::array_descriptor_t& */ _array_desc,
      /* bool* */ _strict,
      /* bool& */ _is_decreasing
  );
  return _is_decreasing;
}
bool SimUtils::is_false(double param) {
  bool _this_false{};
  fortran_is_false(/* double& */ param, /* bool& */ _this_false);
  return _this_false;
}
bool SimUtils::is_increasing_sequence(FArray1D<Real> &array, std::optional<bool> strict) {
  // array: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _array_desc;
  _array_desc.rank = 1;
  _array_desc.data_ptr = array.data();
  _array_desc.dims[0] = array.size();
  bool strict_lvalue;
  auto *_strict{&strict_lvalue};
  if (strict.has_value()) {
    strict_lvalue = strict.value();
  } else {
    _strict = nullptr;
  }
  bool _is_increasing{};
  fortran_is_increasing_sequence(
      /* Bmad::array_descriptor_t& */ _array_desc,
      /* bool* */ _strict,
      /* bool& */ _is_increasing
  );
  return _is_increasing;
}
void SimUtils::is_integer(
    std::string string,
    bool valid,
    std::optional<int> int_,
    std::optional<std::string> delims,
    std::optional<int> ix_word
) {
  auto _string = string.c_str();
  int int__lvalue;
  auto *_int_{&int__lvalue};
  if (int_.has_value()) {
    int__lvalue = int_.value();
  } else {
    _int_ = nullptr;
  }
  const char *_delims = delims.has_value() ? delims->c_str() : nullptr;
  int ix_word_lvalue;
  auto *_ix_word{&ix_word_lvalue};
  if (ix_word.has_value()) {
    ix_word_lvalue = ix_word.value();
  } else {
    _ix_word = nullptr;
  }
  fortran_is_integer(
      /* const char* */ _string,
      /* int* */ _int_,
      /* const char* */ _delims,
      /* int* */ _ix_word,
      /* bool& */ valid
  );
}
void SimUtils::is_logical(std::string string, bool valid, std::optional<bool> ignore) {
  auto _string = string.c_str();
  bool ignore_lvalue;
  auto *_ignore{&ignore_lvalue};
  if (ignore.has_value()) {
    ignore_lvalue = ignore.value();
  } else {
    _ignore = nullptr;
  }
  fortran_is_logical(/* const char* */ _string, /* bool* */ _ignore, /* bool& */ valid);
}
void SimUtils::is_real(
    std::string string,
    bool valid,
    std::optional<bool> ignore,
    std::optional<double> real_num
) {
  auto _string = string.c_str();
  bool ignore_lvalue;
  auto *_ignore{&ignore_lvalue};
  if (ignore.has_value()) {
    ignore_lvalue = ignore.value();
  } else {
    _ignore = nullptr;
  }
  double real_num_lvalue;
  auto *_real_num{&real_num_lvalue};
  if (real_num.has_value()) {
    real_num_lvalue = real_num.value();
  } else {
    _real_num = nullptr;
  }
  fortran_is_real(
      /* const char* */ _string,
      /* bool* */ _ignore,
      /* double* */ _real_num,
      /* bool& */ valid
  );
}
bool SimUtils::is_subatomic_species(int species) {
  bool _is_subatomic{};
  fortran_is_subatomic_species(/* int& */ species, /* bool& */ _is_subatomic);
  return _is_subatomic;
}
bool SimUtils::is_true(double param) {
  bool _this_true{};
  fortran_is_true(/* double& */ param, /* bool& */ _this_true);
  return _this_true;
}
double SimUtils::j_bessel(int m, double arg) {
  double _j_bes{};
  fortran_j_bessel(/* int& */ m, /* double& */ arg, /* double& */ _j_bes);
  return _j_bes;
}
void SimUtils::linear_fit(
    FArray1D<Real> &x,
    FArray1D<Real> &y,
    int n_data,
    double a,
    double b,
    double sig_a,
    double sig_b
) {
  // x: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_desc;
  _x_desc.rank = 1;
  _x_desc.data_ptr = x.data();
  _x_desc.dims[0] = x.size();
  // y: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _y_desc;
  _y_desc.rank = 1;
  _y_desc.data_ptr = y.data();
  _y_desc.dims[0] = y.size();
  fortran_linear_fit(
      /* Bmad::array_descriptor_t& */ _x_desc,
      /* Bmad::array_descriptor_t& */ _y_desc,
      /* int& */ n_data,
      /* double& */ a,
      /* double& */ b,
      /* double& */ sig_a,
      /* double& */ sig_b
  );
}
FixedArray1D<Real, 3>
SimUtils::linear_fit_2d(FArray1D<Real> &x, FArray1D<Real> &y, FArray1D<Real> &z) {
  // x: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_desc;
  _x_desc.rank = 1;
  _x_desc.data_ptr = x.data();
  _x_desc.dims[0] = x.size();
  // y: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _y_desc;
  _y_desc.rank = 1;
  _y_desc.data_ptr = y.data();
  _y_desc.dims[0] = y.size();
  // z: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _z_desc;
  _z_desc.rank = 1;
  _z_desc.data_ptr = z.data();
  _z_desc.dims[0] = z.size();
  // coef: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _coef_desc;
  _coef_desc.rank = 1;
  FixedArray1D<Real, 3> _coef;
  _coef_desc.data_ptr = _coef.data();
  _coef_desc.dims[0] = _coef.size();
  fortran_linear_fit_2d(
      /* Bmad::array_descriptor_t& */ _x_desc,
      /* Bmad::array_descriptor_t& */ _y_desc,
      /* Bmad::array_descriptor_t& */ _z_desc,
      /* Bmad::array_descriptor_t& */ _coef_desc
  );
  return _coef;
}
void SimUtils::logic_str(bool logic, std::string str) {
  auto _str = str.c_str();
  fortran_logic_str(/* bool& */ logic, /* const char* */ _str);
}
void SimUtils::lunget(int func_retval__) { fortran_lunget(/* int& */ func_retval__); }
void SimUtils::make_legal_comment(std::string comment_in, std::string comment_out) {
  auto _comment_in = comment_in.c_str();
  auto _comment_out = comment_out.c_str();
  fortran_make_legal_comment(/* const char* */ _comment_in, /* const char* */ _comment_out);
}
double SimUtils::mass_of(int species) {
  double _mass{};
  fortran_mass_of(/* int& */ species, /* double& */ _mass);
  return _mass;
}
void SimUtils::match_reg(std::string str, std::string pat, bool is_match) {
  auto _str = str.c_str();
  auto _pat = pat.c_str();
  fortran_match_reg(/* const char* */ _str, /* const char* */ _pat, /* bool& */ is_match);
}
void SimUtils::match_wild(std::string string, std::string template_, bool is_match) {
  auto _string = string.c_str();
  auto _template_ = template_.c_str();
  fortran_match_wild(/* const char* */ _string, /* const char* */ _template_, /* bool& */ is_match);
}
void SimUtils::match_word(
    std::string string,
    CharacterAlloc1D &names,
    int ix,
    std::optional<bool> exact_case,
    std::optional<bool> can_abbreviate,
    std::optional<std::string> matched_name
) {
  auto _string = string.c_str();
  // intent=inout character array container
  bool exact_case_lvalue;
  auto *_exact_case{&exact_case_lvalue};
  if (exact_case.has_value()) {
    exact_case_lvalue = exact_case.value();
  } else {
    _exact_case = nullptr;
  }
  bool can_abbreviate_lvalue;
  auto *_can_abbreviate{&can_abbreviate_lvalue};
  if (can_abbreviate.has_value()) {
    can_abbreviate_lvalue = can_abbreviate.value();
  } else {
    _can_abbreviate = nullptr;
  }
  const char *_matched_name = matched_name.has_value() ? matched_name->c_str() : nullptr;
  fortran_match_word(
      /* const char* */ _string,
      /* void* */ names.get_fortran_ptr(),
      /* int& */ ix,
      /* bool* */ _exact_case,
      /* bool* */ _can_abbreviate,
      /* const char* */ _matched_name
  );
}
void SimUtils::maximize_projection(double seed, FArray1D<Complex> &cdata, double func_retval__) {
  // cdata: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _cdata_desc;
  _cdata_desc.rank = 1;
  _cdata_desc.data_ptr = cdata.data();
  _cdata_desc.dims[0] = cdata.size();
  fortran_maximize_projection(
      /* double& */ seed,
      /* Bmad::array_descriptor_t& */ _cdata_desc,
      /* double& */ func_retval__
  );
}
void SimUtils::milli_sleep(int milli_sec) { fortran_milli_sleep(/* int& */ milli_sec); }
double SimUtils::modulo2_dp(double x, double amp) {
  double _mod2{};
  fortran_modulo2_dp(/* double& */ x, /* double& */ amp, /* double& */ _mod2);
  return _mod2;
}
int SimUtils::modulo2_int(int x, int amp) {
  int _mod2{};
  fortran_modulo2_int(/* int& */ x, /* int& */ amp, /* int& */ _mod2);
  return _mod2;
}
double SimUtils::modulo2_qp(double x, double amp) {
  double _mod2{};
  fortran_modulo2_qp(/* double& */ x, /* double& */ amp, /* double& */ _mod2);
  return _mod2;
}
double SimUtils::modulo2_sp(double x, double amp) {
  double _mod2{};
  fortran_modulo2_sp(/* double& */ x, /* double& */ amp, /* double& */ _mod2);
  return _mod2;
}
void SimUtils::n_bins_automatic(int n_data, int n) {
  fortran_n_bins_automatic(/* int& */ n_data, /* int& */ n);
}
double SimUtils::n_choose_k(int n, int k) {
  double _nck{};
  fortran_n_choose_k(/* int& */ n, /* int& */ k, /* double& */ _nck);
  return _nck;
}
void SimUtils::n_spline_create(
    FArray1D<Real> &deriv0,
    FArray1D<Real> &deriv1,
    double x1,
    FArray1D<Real> &n_spline
) {
  // deriv0: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _deriv0_desc;
  _deriv0_desc.rank = 1;
  _deriv0_desc.data_ptr = deriv0.data();
  _deriv0_desc.dims[0] = deriv0.size();
  // deriv1: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _deriv1_desc;
  _deriv1_desc.rank = 1;
  _deriv1_desc.data_ptr = deriv1.data();
  _deriv1_desc.dims[0] = deriv1.size();
  // n_spline: inout NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _n_spline_desc;
  _n_spline_desc.rank = 1;
  _n_spline_desc.data_ptr = n_spline.data();
  _n_spline_desc.dims[0] = n_spline.size();
  fortran_n_spline_create(
      /* Bmad::array_descriptor_t& */ _deriv0_desc,
      /* Bmad::array_descriptor_t& */ _deriv1_desc,
      /* double& */ x1,
      /* Bmad::array_descriptor_t& */ _n_spline_desc
  );
}
void SimUtils::naff(
    FArray1D<Complex> &cdata,
    FArray1D<Real> &freqs,
    FArray1D<Complex> &amps,
    std::optional<int> opt_dump_spectra,
    std::optional<bool> opt_zero_first
) {
  // cdata: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _cdata_desc;
  _cdata_desc.rank = 1;
  _cdata_desc.data_ptr = cdata.data();
  _cdata_desc.dims[0] = cdata.size();
  // freqs: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _freqs_desc;
  _freqs_desc.rank = 1;
  _freqs_desc.data_ptr = freqs.data();
  _freqs_desc.dims[0] = freqs.size();
  // amps: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _amps_desc;
  _amps_desc.rank = 1;
  _amps_desc.data_ptr = amps.data();
  _amps_desc.dims[0] = amps.size();
  int opt_dump_spectra_lvalue;
  auto *_opt_dump_spectra{&opt_dump_spectra_lvalue};
  if (opt_dump_spectra.has_value()) {
    opt_dump_spectra_lvalue = opt_dump_spectra.value();
  } else {
    _opt_dump_spectra = nullptr;
  }
  bool opt_zero_first_lvalue;
  auto *_opt_zero_first{&opt_zero_first_lvalue};
  if (opt_zero_first.has_value()) {
    opt_zero_first_lvalue = opt_zero_first.value();
  } else {
    _opt_zero_first = nullptr;
  }
  fortran_naff(
      /* Bmad::array_descriptor_t& */ _cdata_desc,
      /* Bmad::array_descriptor_t& */ _freqs_desc,
      /* Bmad::array_descriptor_t& */ _amps_desc,
      /* int* */ _opt_dump_spectra,
      /* bool* */ _opt_zero_first
  );
}
void SimUtils::nametable_add(NametableStruct &nametable, std::string name, int ix_name) {
  auto _name = name.c_str();
  fortran_nametable_add(
      /* void* */ nametable.get_fortran_ptr(),
      /* const char* */ _name,
      /* int& */ ix_name
  );
}
void SimUtils::nametable_bracket_indexx(
    NametableStruct &nametable,
    std::string name,
    int ix_max,
    std::optional<int> n_match
) {
  auto _name = name.c_str();
  int n_match_lvalue;
  auto *_n_match{&n_match_lvalue};
  if (n_match.has_value()) {
    n_match_lvalue = n_match.value();
  } else {
    _n_match = nullptr;
  }
  fortran_nametable_bracket_indexx(
      /* void* */ nametable.get_fortran_ptr(),
      /* const char* */ _name,
      /* int* */ _n_match,
      /* int& */ ix_max
  );
}
void SimUtils::nametable_change1(NametableStruct &nametable, std::string name, int ix_name) {
  auto _name = name.c_str();
  fortran_nametable_change1(
      /* void* */ nametable.get_fortran_ptr(),
      /* const char* */ _name,
      /* int& */ ix_name
  );
}
void SimUtils::nametable_init(
    NametableStruct &nametable,
    std::optional<int> n_min,
    std::optional<int> n_max
) {
  int n_min_lvalue;
  auto *_n_min{&n_min_lvalue};
  if (n_min.has_value()) {
    n_min_lvalue = n_min.value();
  } else {
    _n_min = nullptr;
  }
  int n_max_lvalue;
  auto *_n_max{&n_max_lvalue};
  if (n_max.has_value()) {
    n_max_lvalue = n_max.value();
  } else {
    _n_max = nullptr;
  }
  fortran_nametable_init(
      /* void* */ nametable.get_fortran_ptr(),
      /* int* */ _n_min,
      /* int* */ _n_max
  );
}
void SimUtils::nametable_remove(NametableStruct &nametable, int ix_name) {
  fortran_nametable_remove(/* void* */ nametable.get_fortran_ptr(), /* int& */ ix_name);
}
void SimUtils::negative_ampsquared(double frequency, double amp, std::optional<int> status) {
  int status_lvalue;
  auto *_status{&status_lvalue};
  if (status.has_value()) {
    status_lvalue = status.value();
  } else {
    _status = nullptr;
  }
  fortran_negative_ampsquared(/* double& */ frequency, /* int* */ _status, /* double& */ amp);
}
void SimUtils::negative_dampsquared(double frequency, double damp, std::optional<int> status) {
  int status_lvalue;
  auto *_status{&status_lvalue};
  if (status.has_value()) {
    status_lvalue = status.value();
  } else {
    _status = nullptr;
  }
  fortran_negative_dampsquared(/* double& */ frequency, /* int* */ _status, /* double& */ damp);
}
FixedArray1D<Real, 4> SimUtils::omega_to_quat(FixedArray1D<Real, 3> omega) {
  // omega: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _omega_desc;
  _omega_desc.rank = 1;
  _omega_desc.data_ptr = omega.data();
  _omega_desc.dims[0] = omega.size();
  // quat: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  FixedArray1D<Real, 4> _quat;
  _quat_desc.data_ptr = _quat.data();
  _quat_desc.dims[0] = _quat.size();
  fortran_omega_to_quat(
      /* Bmad::array_descriptor_t& */ _omega_desc,
      /* Bmad::array_descriptor_t& */ _quat_desc
  );
  return _quat;
}
std::string SimUtils::openpmd_species_name(int species) {
  char _pmd_name[4096];
  fortran_openpmd_species_name(/* int& */ species, /* const char* */ _pmd_name);
  return _pmd_name;
}
void SimUtils::ordinal_str(int n, std::string str) {
  auto _str = str.c_str();
  fortran_ordinal_str(/* int& */ n, /* const char* */ _str);
}
void SimUtils::out_io_buffer_get_line(int ix_line, std::string line) {
  auto _line = line.c_str();
  fortran_out_io_buffer_get_line(/* int& */ ix_line, /* const char* */ _line);
}
void SimUtils::out_io_buffer_num_lines(int n_lines) {
  fortran_out_io_buffer_num_lines(/* int& */ n_lines);
}
void SimUtils::out_io_buffer_reset() { fortran_out_io_buffer_reset(); }
void SimUtils::out_io(
    int level,
    std::string routine_name,
    std::string line,
    int i_num,
    std::optional<bool> insert_tag_line
) {
  auto _routine_name = routine_name.c_str();
  auto _line = line.c_str();
  bool insert_tag_line_lvalue;
  auto *_insert_tag_line{&insert_tag_line_lvalue};
  if (insert_tag_line.has_value()) {
    insert_tag_line_lvalue = insert_tag_line.value();
  } else {
    _insert_tag_line = nullptr;
  }
  fortran_out_io_int(
      /* int& */ level,
      /* const char* */ _routine_name,
      /* const char* */ _line,
      /* int& */ i_num,
      /* bool* */ _insert_tag_line
  );
}
void SimUtils::out_io(
    int level,
    std::string routine_name,
    std::string line1,
    std::optional<std::string> line2,
    std::optional<std::string> line3,
    std::optional<std::string> line4,
    std::optional<std::string> line5,
    std::optional<std::string> line6,
    std::optional<std::string> line7,
    std::optional<std::string> line8,
    std::optional<std::string> line9,
    std::optional<std::string> line10,
    std::optional<std::string> line11,
    std::optional<std::string> line12,
    std::optional<FArray1D<Real>> r_array,
    std::optional<FArray1D<Int>> i_array,
    optional_ref<BoolAlloc1D> l_array,
    std::optional<bool> insert_tag_line
) {
  auto _routine_name = routine_name.c_str();
  auto _line1 = line1.c_str();
  const char *_line2 = line2.has_value() ? line2->c_str() : nullptr;
  const char *_line3 = line3.has_value() ? line3->c_str() : nullptr;
  const char *_line4 = line4.has_value() ? line4->c_str() : nullptr;
  const char *_line5 = line5.has_value() ? line5->c_str() : nullptr;
  const char *_line6 = line6.has_value() ? line6->c_str() : nullptr;
  const char *_line7 = line7.has_value() ? line7->c_str() : nullptr;
  const char *_line8 = line8.has_value() ? line8->c_str() : nullptr;
  const char *_line9 = line9.has_value() ? line9->c_str() : nullptr;
  const char *_line10 = line10.has_value() ? line10->c_str() : nullptr;
  const char *_line11 = line11.has_value() ? line11->c_str() : nullptr;
  const char *_line12 = line12.has_value() ? line12->c_str() : nullptr;
  // r_array: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _r_array_desc;
  _r_array_desc.rank = 1;
  if (r_array.has_value()) {
    _r_array_desc.data_ptr = r_array->data();
    _r_array_desc.dims[0] = r_array->size();
  } else {
    _r_array_desc.data_ptr = nullptr;
    _r_array_desc.dims[0] = 0;
  }
  // i_array: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _i_array_desc;
  _i_array_desc.rank = 1;
  if (i_array.has_value()) {
    _i_array_desc.data_ptr = i_array->data();
    _i_array_desc.dims[0] = i_array->size();
  } else {
    _i_array_desc.data_ptr = nullptr;
    _i_array_desc.dims[0] = 0;
  }
  // intent=inout allocatable general array
  auto *_l_array =
      l_array.has_value() ? l_array->get().get_fortran_ptr() : nullptr; // input, optional
  bool insert_tag_line_lvalue;
  auto *_insert_tag_line{&insert_tag_line_lvalue};
  if (insert_tag_line.has_value()) {
    insert_tag_line_lvalue = insert_tag_line.value();
  } else {
    _insert_tag_line = nullptr;
  }
  fortran_out_io_line12(
      /* int& */ level,
      /* const char* */ _routine_name,
      /* const char* */ _line1,
      /* const char* */ _line2,
      /* const char* */ _line3,
      /* const char* */ _line4,
      /* const char* */ _line5,
      /* const char* */ _line6,
      /* const char* */ _line7,
      /* const char* */ _line8,
      /* const char* */ _line9,
      /* const char* */ _line10,
      /* const char* */ _line11,
      /* const char* */ _line12,
      /* Bmad::array_descriptor_t& */ _r_array_desc,
      /* Bmad::array_descriptor_t& */ _i_array_desc,
      /* void* */ _l_array,
      /* bool* */ _insert_tag_line
  );
}
void SimUtils::out_io(
    int level,
    std::string routine_name,
    CharacterAlloc1D &lines,
    std::optional<FArray1D<Real>> r_array,
    std::optional<FArray1D<Int>> i_array,
    optional_ref<BoolAlloc1D> l_array,
    std::optional<bool> insert_tag_line
) {
  auto _routine_name = routine_name.c_str();
  // intent=inout character array container
  // r_array: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _r_array_desc;
  _r_array_desc.rank = 1;
  if (r_array.has_value()) {
    _r_array_desc.data_ptr = r_array->data();
    _r_array_desc.dims[0] = r_array->size();
  } else {
    _r_array_desc.data_ptr = nullptr;
    _r_array_desc.dims[0] = 0;
  }
  // i_array: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _i_array_desc;
  _i_array_desc.rank = 1;
  if (i_array.has_value()) {
    _i_array_desc.data_ptr = i_array->data();
    _i_array_desc.dims[0] = i_array->size();
  } else {
    _i_array_desc.data_ptr = nullptr;
    _i_array_desc.dims[0] = 0;
  }
  // intent=inout allocatable general array
  auto *_l_array =
      l_array.has_value() ? l_array->get().get_fortran_ptr() : nullptr; // input, optional
  bool insert_tag_line_lvalue;
  auto *_insert_tag_line{&insert_tag_line_lvalue};
  if (insert_tag_line.has_value()) {
    insert_tag_line_lvalue = insert_tag_line.value();
  } else {
    _insert_tag_line = nullptr;
  }
  fortran_out_io_lines(
      /* int& */ level,
      /* const char* */ _routine_name,
      /* void* */ lines.get_fortran_ptr(),
      /* Bmad::array_descriptor_t& */ _r_array_desc,
      /* Bmad::array_descriptor_t& */ _i_array_desc,
      /* void* */ _l_array,
      /* bool* */ _insert_tag_line
  );
}
void SimUtils::out_io(
    int level,
    std::string routine_name,
    std::string line,
    bool l_num,
    std::optional<bool> insert_tag_line
) {
  auto _routine_name = routine_name.c_str();
  auto _line = line.c_str();
  bool insert_tag_line_lvalue;
  auto *_insert_tag_line{&insert_tag_line_lvalue};
  if (insert_tag_line.has_value()) {
    insert_tag_line_lvalue = insert_tag_line.value();
  } else {
    _insert_tag_line = nullptr;
  }
  fortran_out_io_logical(
      /* int& */ level,
      /* const char* */ _routine_name,
      /* const char* */ _line,
      /* bool& */ l_num,
      /* bool* */ _insert_tag_line
  );
}
void SimUtils::out_io_print_and_capture_setup(
    std::optional<bool> print_on,
    std::optional<std::string> capture_state,
    std::optional<bool> capture_add_null
) {
  bool print_on_lvalue;
  auto *_print_on{&print_on_lvalue};
  if (print_on.has_value()) {
    print_on_lvalue = print_on.value();
  } else {
    _print_on = nullptr;
  }
  const char *_capture_state = capture_state.has_value() ? capture_state->c_str() : nullptr;
  bool capture_add_null_lvalue;
  auto *_capture_add_null{&capture_add_null_lvalue};
  if (capture_add_null.has_value()) {
    capture_add_null_lvalue = capture_add_null.value();
  } else {
    _capture_add_null = nullptr;
  }
  fortran_out_io_print_and_capture_setup(
      /* bool* */ _print_on,
      /* const char* */ _capture_state,
      /* bool* */ _capture_add_null
  );
}
void SimUtils::out_io(
    int level,
    std::string routine_name,
    std::string line,
    double r_num,
    std::optional<bool> insert_tag_line
) {
  auto _routine_name = routine_name.c_str();
  auto _line = line.c_str();
  bool insert_tag_line_lvalue;
  auto *_insert_tag_line{&insert_tag_line_lvalue};
  if (insert_tag_line.has_value()) {
    insert_tag_line_lvalue = insert_tag_line.value();
  } else {
    _insert_tag_line = nullptr;
  }
  fortran_out_io_real(
      /* int& */ level,
      /* const char* */ _routine_name,
      /* const char* */ _line,
      /* double& */ r_num,
      /* bool* */ _insert_tag_line
  );
}
void SimUtils::parse_fortran_format(
    std::string format_str,
    int n_repeat,
    int power,
    std::string descrip,
    int width,
    int digits
) {
  auto _format_str = format_str.c_str();
  auto _descrip = descrip.c_str();
  fortran_parse_fortran_format(
      /* const char* */ _format_str,
      /* int& */ n_repeat,
      /* int& */ power,
      /* const char* */ _descrip,
      /* int& */ width,
      /* int& */ digits
  );
}
void SimUtils::pointer_to_locations(
    std::string string,
    IntAlloc1D &array,
    int num,
    int ix_min,
    int ix_max,
    optional_ref<CharacterAlloc1D> names,
    std::optional<bool> exact_case,
    std::optional<bool> print_err
) {
  auto _string = string.c_str();
  // intent=inout allocatable general array
  // intent=inout character array container
  auto *_names = names.has_value() ? names->get().get_fortran_ptr() : nullptr; // input, optional
  bool exact_case_lvalue;
  auto *_exact_case{&exact_case_lvalue};
  if (exact_case.has_value()) {
    exact_case_lvalue = exact_case.value();
  } else {
    _exact_case = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  fortran_pointer_to_locations(
      /* const char* */ _string,
      /* void* */ array.get_fortran_ptr(),
      /* int& */ num,
      /* int& */ ix_min,
      /* int& */ ix_max,
      /* void* */ _names,
      /* bool* */ _exact_case,
      /* bool* */ _print_err
  );
}
std::optional<RandomStateStruct> SimUtils::pointer_to_ran_state(
    optional_ref<RandomStateStruct> ran_state,
    std::optional<int> ix_thread
) {
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  int ix_thread_lvalue;
  auto *_ix_thread{&ix_thread_lvalue};
  if (ix_thread.has_value()) {
    ix_thread_lvalue = ix_thread.value();
  } else {
    _ix_thread = nullptr;
  }
  void *_ran_state_ptr;
  fortran_pointer_to_ran_state(
      /* void* */ _ran_state,
      /* int* */ _ix_thread,
      /* void* */ &_ran_state_ptr
  );
  return std::move(
      (_ran_state_ptr ? std::make_optional<RandomStateStruct>(_ran_state_ptr) : std::nullopt)
  );
}
double SimUtils::poly_eval(FArray1D<Real> &poly, double x, std::optional<bool> diff_coef) {
  // poly: in NOT (CppWrapperGeneralArgumentArray) (['0:'])
  Bmad::array_descriptor_t _poly_desc;
  _poly_desc.rank = 1;
  _poly_desc.data_ptr = poly.data();
  _poly_desc.dims[0] = poly.size();
  bool diff_coef_lvalue;
  auto *_diff_coef{&diff_coef_lvalue};
  if (diff_coef.has_value()) {
    diff_coef_lvalue = diff_coef.value();
  } else {
    _diff_coef = nullptr;
  }
  double _y{};
  fortran_poly_eval(
      /* Bmad::array_descriptor_t& */ _poly_desc,
      /* double& */ x,
      /* bool* */ _diff_coef,
      /* double& */ _y
  );
  return _y;
}
void SimUtils::probability_funct(double x, double prob) {
  fortran_probability_funct(/* double& */ x, /* double& */ prob);
}
void SimUtils::projdd(
    FArray1D<Complex> &a,
    FArray1D<Complex> &b,
    std::complex<double> func_retval__
) {
  // a: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _a_desc;
  _a_desc.rank = 1;
  _a_desc.data_ptr = a.data();
  _a_desc.dims[0] = a.size();
  // b: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _b_desc;
  _b_desc.rank = 1;
  _b_desc.data_ptr = b.data();
  _b_desc.dims[0] = b.size();
  fortran_projdd(
      /* Bmad::array_descriptor_t& */ _a_desc,
      /* Bmad::array_descriptor_t& */ _b_desc,
      /* std::complex<double>& */ func_retval__
  );
}
FixedArray1D<Complex, 2> SimUtils::quadratic_roots(FixedArray1D<Real, 3> coefs) {
  // coefs: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _coefs_desc;
  _coefs_desc.rank = 1;
  _coefs_desc.data_ptr = coefs.data();
  _coefs_desc.dims[0] = coefs.size();
  // root: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _root_desc;
  _root_desc.rank = 1;
  FixedArray1D<Complex, 2> _root;
  _root_desc.data_ptr = _root.data();
  _root_desc.dims[0] = _root.size();
  fortran_quadratic_roots(
      /* Bmad::array_descriptor_t& */ _coefs_desc,
      /* Bmad::array_descriptor_t& */ _root_desc
  );
  return _root;
}
FixedArray1D<Complex, 4> SimUtils::quat_conj(FixedArray1D<Complex, 4> q_in) {
  // q_in: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_in_desc;
  _q_in_desc.rank = 1;
  _q_in_desc.data_ptr = q_in.data();
  _q_in_desc.dims[0] = q_in.size();
  // q_out: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_out_desc;
  _q_out_desc.rank = 1;
  FixedArray1D<Complex, 4> _q_out;
  _q_out_desc.data_ptr = _q_out.data();
  _q_out_desc.dims[0] = _q_out.size();
  fortran_quat_conj_complex(
      /* Bmad::array_descriptor_t& */ _q_in_desc,
      /* Bmad::array_descriptor_t& */ _q_out_desc
  );
  return _q_out;
}
FixedArray1D<Real, 4> SimUtils::quat_conj(FixedArray1D<Real, 4> q_in) {
  // q_in: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_in_desc;
  _q_in_desc.rank = 1;
  _q_in_desc.data_ptr = q_in.data();
  _q_in_desc.dims[0] = q_in.size();
  // q_out: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_out_desc;
  _q_out_desc.rank = 1;
  FixedArray1D<Real, 4> _q_out;
  _q_out_desc.data_ptr = _q_out.data();
  _q_out_desc.dims[0] = _q_out.size();
  fortran_quat_conj_real(
      /* Bmad::array_descriptor_t& */ _q_in_desc,
      /* Bmad::array_descriptor_t& */ _q_out_desc
  );
  return _q_out;
}
FixedArray1D<Real, 4> SimUtils::quat_inverse(FixedArray1D<Real, 4> q_in) {
  // q_in: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_in_desc;
  _q_in_desc.rank = 1;
  _q_in_desc.data_ptr = q_in.data();
  _q_in_desc.dims[0] = q_in.size();
  // q_out: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_out_desc;
  _q_out_desc.rank = 1;
  FixedArray1D<Real, 4> _q_out;
  _q_out_desc.data_ptr = _q_out.data();
  _q_out_desc.dims[0] = _q_out.size();
  fortran_quat_inverse(
      /* Bmad::array_descriptor_t& */ _q_in_desc,
      /* Bmad::array_descriptor_t& */ _q_out_desc
  );
  return _q_out;
}
FixedArray1D<Complex, 4> SimUtils::quat_mul(
    FixedArray1D<Complex, 4> q1,
    FixedArray1D<Complex, 4> q2,
    std::optional<FixedArray1D<Complex, 4>> q3,
    std::optional<FixedArray1D<Complex, 4>> q4,
    std::optional<FixedArray1D<Complex, 4>> q5,
    std::optional<FixedArray1D<Complex, 4>> q6,
    std::optional<FixedArray1D<Complex, 4>> q7,
    std::optional<FixedArray1D<Complex, 4>> q8,
    std::optional<FixedArray1D<Complex, 4>> q9
) {
  // q1: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q1_desc;
  _q1_desc.rank = 1;
  _q1_desc.data_ptr = q1.data();
  _q1_desc.dims[0] = q1.size();
  // q2: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q2_desc;
  _q2_desc.rank = 1;
  _q2_desc.data_ptr = q2.data();
  _q2_desc.dims[0] = q2.size();
  // q3: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q3_desc;
  _q3_desc.rank = 1;
  if (q3.has_value()) {
    _q3_desc.data_ptr = q3->data();
    _q3_desc.dims[0] = q3->size();
  } else {
    _q3_desc.data_ptr = nullptr;
    _q3_desc.dims[0] = 0;
  }
  // q4: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q4_desc;
  _q4_desc.rank = 1;
  if (q4.has_value()) {
    _q4_desc.data_ptr = q4->data();
    _q4_desc.dims[0] = q4->size();
  } else {
    _q4_desc.data_ptr = nullptr;
    _q4_desc.dims[0] = 0;
  }
  // q5: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q5_desc;
  _q5_desc.rank = 1;
  if (q5.has_value()) {
    _q5_desc.data_ptr = q5->data();
    _q5_desc.dims[0] = q5->size();
  } else {
    _q5_desc.data_ptr = nullptr;
    _q5_desc.dims[0] = 0;
  }
  // q6: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q6_desc;
  _q6_desc.rank = 1;
  if (q6.has_value()) {
    _q6_desc.data_ptr = q6->data();
    _q6_desc.dims[0] = q6->size();
  } else {
    _q6_desc.data_ptr = nullptr;
    _q6_desc.dims[0] = 0;
  }
  // q7: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q7_desc;
  _q7_desc.rank = 1;
  if (q7.has_value()) {
    _q7_desc.data_ptr = q7->data();
    _q7_desc.dims[0] = q7->size();
  } else {
    _q7_desc.data_ptr = nullptr;
    _q7_desc.dims[0] = 0;
  }
  // q8: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q8_desc;
  _q8_desc.rank = 1;
  if (q8.has_value()) {
    _q8_desc.data_ptr = q8->data();
    _q8_desc.dims[0] = q8->size();
  } else {
    _q8_desc.data_ptr = nullptr;
    _q8_desc.dims[0] = 0;
  }
  // q9: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q9_desc;
  _q9_desc.rank = 1;
  if (q9.has_value()) {
    _q9_desc.data_ptr = q9->data();
    _q9_desc.dims[0] = q9->size();
  } else {
    _q9_desc.data_ptr = nullptr;
    _q9_desc.dims[0] = 0;
  }
  // q_out: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_out_desc;
  _q_out_desc.rank = 1;
  FixedArray1D<Complex, 4> _q_out;
  _q_out_desc.data_ptr = _q_out.data();
  _q_out_desc.dims[0] = _q_out.size();
  fortran_quat_mul_complex(
      /* Bmad::array_descriptor_t& */ _q1_desc,
      /* Bmad::array_descriptor_t& */ _q2_desc,
      /* Bmad::array_descriptor_t& */ _q3_desc,
      /* Bmad::array_descriptor_t& */ _q4_desc,
      /* Bmad::array_descriptor_t& */ _q5_desc,
      /* Bmad::array_descriptor_t& */ _q6_desc,
      /* Bmad::array_descriptor_t& */ _q7_desc,
      /* Bmad::array_descriptor_t& */ _q8_desc,
      /* Bmad::array_descriptor_t& */ _q9_desc,
      /* Bmad::array_descriptor_t& */ _q_out_desc
  );
  return _q_out;
}
FixedArray1D<Real, 4> SimUtils::quat_mul(
    FixedArray1D<Real, 4> q1,
    FixedArray1D<Real, 4> q2,
    std::optional<FixedArray1D<Real, 4>> q3,
    std::optional<FixedArray1D<Real, 4>> q4,
    std::optional<FixedArray1D<Real, 4>> q5,
    std::optional<FixedArray1D<Real, 4>> q6,
    std::optional<FixedArray1D<Real, 4>> q7,
    std::optional<FixedArray1D<Real, 4>> q8,
    std::optional<FixedArray1D<Real, 4>> q9
) {
  // q1: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q1_desc;
  _q1_desc.rank = 1;
  _q1_desc.data_ptr = q1.data();
  _q1_desc.dims[0] = q1.size();
  // q2: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q2_desc;
  _q2_desc.rank = 1;
  _q2_desc.data_ptr = q2.data();
  _q2_desc.dims[0] = q2.size();
  // q3: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q3_desc;
  _q3_desc.rank = 1;
  if (q3.has_value()) {
    _q3_desc.data_ptr = q3->data();
    _q3_desc.dims[0] = q3->size();
  } else {
    _q3_desc.data_ptr = nullptr;
    _q3_desc.dims[0] = 0;
  }
  // q4: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q4_desc;
  _q4_desc.rank = 1;
  if (q4.has_value()) {
    _q4_desc.data_ptr = q4->data();
    _q4_desc.dims[0] = q4->size();
  } else {
    _q4_desc.data_ptr = nullptr;
    _q4_desc.dims[0] = 0;
  }
  // q5: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q5_desc;
  _q5_desc.rank = 1;
  if (q5.has_value()) {
    _q5_desc.data_ptr = q5->data();
    _q5_desc.dims[0] = q5->size();
  } else {
    _q5_desc.data_ptr = nullptr;
    _q5_desc.dims[0] = 0;
  }
  // q6: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q6_desc;
  _q6_desc.rank = 1;
  if (q6.has_value()) {
    _q6_desc.data_ptr = q6->data();
    _q6_desc.dims[0] = q6->size();
  } else {
    _q6_desc.data_ptr = nullptr;
    _q6_desc.dims[0] = 0;
  }
  // q7: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q7_desc;
  _q7_desc.rank = 1;
  if (q7.has_value()) {
    _q7_desc.data_ptr = q7->data();
    _q7_desc.dims[0] = q7->size();
  } else {
    _q7_desc.data_ptr = nullptr;
    _q7_desc.dims[0] = 0;
  }
  // q8: inout NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q8_desc;
  _q8_desc.rank = 1;
  if (q8.has_value()) {
    _q8_desc.data_ptr = q8->data();
    _q8_desc.dims[0] = q8->size();
  } else {
    _q8_desc.data_ptr = nullptr;
    _q8_desc.dims[0] = 0;
  }
  // q9: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q9_desc;
  _q9_desc.rank = 1;
  if (q9.has_value()) {
    _q9_desc.data_ptr = q9->data();
    _q9_desc.dims[0] = q9->size();
  } else {
    _q9_desc.data_ptr = nullptr;
    _q9_desc.dims[0] = 0;
  }
  // q_out: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _q_out_desc;
  _q_out_desc.rank = 1;
  FixedArray1D<Real, 4> _q_out;
  _q_out_desc.data_ptr = _q_out.data();
  _q_out_desc.dims[0] = _q_out.size();
  fortran_quat_mul_real(
      /* Bmad::array_descriptor_t& */ _q1_desc,
      /* Bmad::array_descriptor_t& */ _q2_desc,
      /* Bmad::array_descriptor_t& */ _q3_desc,
      /* Bmad::array_descriptor_t& */ _q4_desc,
      /* Bmad::array_descriptor_t& */ _q5_desc,
      /* Bmad::array_descriptor_t& */ _q6_desc,
      /* Bmad::array_descriptor_t& */ _q7_desc,
      /* Bmad::array_descriptor_t& */ _q8_desc,
      /* Bmad::array_descriptor_t& */ _q9_desc,
      /* Bmad::array_descriptor_t& */ _q_out_desc
  );
  return _q_out;
}
FixedArray1D<Complex, 3>
SimUtils::quat_rotate(FixedArray1D<Complex, 4> quat, FixedArray1D<Complex, 3> vec_in) {
  // quat: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  _quat_desc.data_ptr = quat.data();
  _quat_desc.dims[0] = quat.size();
  // vec_in: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_in_desc;
  _vec_in_desc.rank = 1;
  _vec_in_desc.data_ptr = vec_in.data();
  _vec_in_desc.dims[0] = vec_in.size();
  // vec_out: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_out_desc;
  _vec_out_desc.rank = 1;
  FixedArray1D<Complex, 3> _vec_out;
  _vec_out_desc.data_ptr = _vec_out.data();
  _vec_out_desc.dims[0] = _vec_out.size();
  fortran_quat_rotate_complex(
      /* Bmad::array_descriptor_t& */ _quat_desc,
      /* Bmad::array_descriptor_t& */ _vec_in_desc,
      /* Bmad::array_descriptor_t& */ _vec_out_desc
  );
  return _vec_out;
}
FixedArray1D<Real, 3>
SimUtils::quat_rotate(FixedArray1D<Real, 4> quat, FixedArray1D<Real, 3> vec_in) {
  // quat: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  _quat_desc.data_ptr = quat.data();
  _quat_desc.dims[0] = quat.size();
  // vec_in: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_in_desc;
  _vec_in_desc.rank = 1;
  _vec_in_desc.data_ptr = vec_in.data();
  _vec_in_desc.dims[0] = vec_in.size();
  // vec_out: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_out_desc;
  _vec_out_desc.rank = 1;
  FixedArray1D<Real, 3> _vec_out;
  _vec_out_desc.data_ptr = _vec_out.data();
  _vec_out_desc.dims[0] = _vec_out.size();
  fortran_quat_rotate_real(
      /* Bmad::array_descriptor_t& */ _quat_desc,
      /* Bmad::array_descriptor_t& */ _vec_in_desc,
      /* Bmad::array_descriptor_t& */ _vec_out_desc
  );
  return _vec_out;
}
SimUtils::QuatToAxisAngle SimUtils::quat_to_axis_angle(FixedArray1D<Real, 4> quat) {
  // quat: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  _quat_desc.data_ptr = quat.data();
  _quat_desc.dims[0] = quat.size();
  // axis: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _axis_desc;
  _axis_desc.rank = 1;
  FixedArray1D<Real, 3> _axis;
  _axis_desc.data_ptr = _axis.data();
  _axis_desc.dims[0] = _axis.size();
  double _angle{};
  fortran_quat_to_axis_angle(
      /* Bmad::array_descriptor_t& */ _quat_desc,
      /* Bmad::array_descriptor_t& */ _axis_desc,
      /* double& */ _angle
  );
  return QuatToAxisAngle{_axis, _angle};
}
FixedArray1D<Real, 3> SimUtils::quat_to_omega(FixedArray1D<Real, 4> quat) {
  // quat: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  _quat_desc.data_ptr = quat.data();
  _quat_desc.dims[0] = quat.size();
  // omega: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _omega_desc;
  _omega_desc.rank = 1;
  FixedArray1D<Real, 3> _omega;
  _omega_desc.data_ptr = _omega.data();
  _omega_desc.dims[0] = _omega.size();
  fortran_quat_to_omega(
      /* Bmad::array_descriptor_t& */ _quat_desc,
      /* Bmad::array_descriptor_t& */ _omega_desc
  );
  return _omega;
}
FixedArray2D<Real, 3, 3> SimUtils::quat_to_w_mat(FixedArray1D<Real, 4> quat) {
  // quat: in NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  _quat_desc.data_ptr = quat.data();
  _quat_desc.dims[0] = quat.size();
  // w_mat: out NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  fortran_quat_to_w_mat(
      /* Bmad::array_descriptor_t& */ _quat_desc,
      /* Bmad::array_descriptor_t& */ _w_mat_desc
  );
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
void SimUtils::query_string(
    std::string query_str,
    bool upcase,
    std::string return_str,
    int ix,
    int ios
) {
  auto _query_str = query_str.c_str();
  auto _return_str = return_str.c_str();
  fortran_query_string(
      /* const char* */ _query_str,
      /* bool& */ upcase,
      /* const char* */ _return_str,
      /* int& */ ix,
      /* int& */ ios
  );
}
void SimUtils::quote(std::string str, std::string q_str) {
  auto _str = str.c_str();
  auto _q_str = q_str.c_str();
  fortran_quote(/* const char* */ _str, /* const char* */ _q_str);
}
void SimUtils::quoten(CharacterAlloc1D &str, std::string q_str, std::optional<std::string> delim) {
  // intent=inout character array container
  const char *_delim = delim.has_value() ? delim->c_str() : nullptr;
  auto _q_str = q_str.c_str();
  fortran_quoten(
      /* void* */ str.get_fortran_ptr(),
      /* const char* */ _delim,
      /* const char* */ _q_str
  );
}
RandomStateStruct SimUtils::ran_default_state(optional_ref<RandomStateStruct> set_state) {
  auto *_set_state =
      set_state.has_value() ? set_state->get().get_fortran_ptr() : nullptr; // input, optional
  RandomStateStruct _get_state;
  fortran_ran_default_state(/* void* */ _set_state, /* void* */ _get_state.get_fortran_ptr());
  return std::move(_get_state);
}
void SimUtils::ran_engine(
    std::optional<std::string> set,
    std::optional<std::string> get,
    optional_ref<RandomStateStruct> ran_state
) {
  const char *_set = set.has_value() ? set->c_str() : nullptr;
  const char *_get = get.has_value() ? get->c_str() : nullptr;
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_ran_engine(/* const char* */ _set, /* const char* */ _get, /* void* */ _ran_state);
}
SimUtils::RanGaussConverter SimUtils::ran_gauss_converter(
    std::optional<std::string> set,
    std::optional<double> set_sigma_cut,
    optional_ref<RandomStateStruct> ran_state
) {
  const char *_set = set.has_value() ? set->c_str() : nullptr;
  double set_sigma_cut_lvalue;
  auto *_set_sigma_cut{&set_sigma_cut_lvalue};
  if (set_sigma_cut.has_value()) {
    set_sigma_cut_lvalue = set_sigma_cut.value();
  } else {
    _set_sigma_cut = nullptr;
  }
  char _get[4096];
  double _get_sigma_cut{};
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_ran_gauss_converter(
      /* const char* */ _set,
      /* double* */ _set_sigma_cut,
      /* const char* */ _get,
      /* double& */ _get_sigma_cut,
      /* void* */ _ran_state
  );
  return RanGaussConverter{_get, _get_sigma_cut};
}
double SimUtils::ran_gauss_scalar(
    optional_ref<RandomStateStruct> ran_state,
    std::optional<double> sigma_cut,
    std::optional<int> index_quasi
) {
  double _harvest{};
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  double sigma_cut_lvalue;
  auto *_sigma_cut{&sigma_cut_lvalue};
  if (sigma_cut.has_value()) {
    sigma_cut_lvalue = sigma_cut.value();
  } else {
    _sigma_cut = nullptr;
  }
  int index_quasi_lvalue;
  auto *_index_quasi{&index_quasi_lvalue};
  if (index_quasi.has_value()) {
    index_quasi_lvalue = index_quasi.value();
  } else {
    _index_quasi = nullptr;
  }
  fortran_ran_gauss_scalar(
      /* double& */ _harvest,
      /* void* */ _ran_state,
      /* double* */ _sigma_cut,
      /* int* */ _index_quasi
  );
  return _harvest;
}
void SimUtils::ran_gauss_vector(
    FArray1D<Real> &harvest,
    optional_ref<RandomStateStruct> ran_state,
    std::optional<double> sigma_cut
) {
  // harvest: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _harvest_desc;
  _harvest_desc.rank = 1;
  _harvest_desc.data_ptr = harvest.data();
  _harvest_desc.dims[0] = harvest.size();
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  double sigma_cut_lvalue;
  auto *_sigma_cut{&sigma_cut_lvalue};
  if (sigma_cut.has_value()) {
    sigma_cut_lvalue = sigma_cut.value();
  } else {
    _sigma_cut = nullptr;
  }
  fortran_ran_gauss_vector(
      /* Bmad::array_descriptor_t& */ _harvest_desc,
      /* void* */ _ran_state,
      /* double* */ _sigma_cut
  );
}
int SimUtils::ran_seed_get() {
  int _seed{};
  fortran_ran_seed_get(/* int& */ _seed);
  return _seed;
}
void SimUtils::ran_seed_put(int seed, std::optional<int> mpi_offset) {
  int mpi_offset_lvalue;
  auto *_mpi_offset{&mpi_offset_lvalue};
  if (mpi_offset.has_value()) {
    mpi_offset_lvalue = mpi_offset.value();
  } else {
    _mpi_offset = nullptr;
  }
  fortran_ran_seed_put(/* int& */ seed, /* int* */ _mpi_offset);
}
double
SimUtils::ran_uniform(optional_ref<RandomStateStruct> ran_state, std::optional<int> index_quasi) {
  double _harvest{};
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  int index_quasi_lvalue;
  auto *_index_quasi{&index_quasi_lvalue};
  if (index_quasi.has_value()) {
    index_quasi_lvalue = index_quasi.value();
  } else {
    _index_quasi = nullptr;
  }
  fortran_ran_uniform_scalar(
      /* double& */ _harvest,
      /* void* */ _ran_state,
      /* int* */ _index_quasi
  );
  return _harvest;
}
void SimUtils::ran_uniform(FArray1D<Real> &harvest, optional_ref<RandomStateStruct> ran_state) {
  // harvest: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _harvest_desc;
  _harvest_desc.rank = 1;
  _harvest_desc.data_ptr = harvest.data();
  _harvest_desc.dims[0] = harvest.size();
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_ran_uniform_vector(/* Bmad::array_descriptor_t& */ _harvest_desc, /* void* */ _ran_state);
}
void SimUtils::rcelbd(double mc, double elb, double eld) {
  fortran_rcelbd(/* double& */ mc, /* double& */ elb, /* double& */ eld);
}
std::string SimUtils::read_a_line(
    std::string prompt,
    std::optional<bool> trim_prompt,
    std::optional<std::string> prompt_color,
    std::optional<bool> prompt_bold,
    std::optional<std::string> history_file
) {
  auto _prompt = prompt.c_str();
  char _line_out[4096];
  bool trim_prompt_lvalue;
  auto *_trim_prompt{&trim_prompt_lvalue};
  if (trim_prompt.has_value()) {
    trim_prompt_lvalue = trim_prompt.value();
  } else {
    _trim_prompt = nullptr;
  }
  const char *_prompt_color = prompt_color.has_value() ? prompt_color->c_str() : nullptr;
  bool prompt_bold_lvalue;
  auto *_prompt_bold{&prompt_bold_lvalue};
  if (prompt_bold.has_value()) {
    prompt_bold_lvalue = prompt_bold.value();
  } else {
    _prompt_bold = nullptr;
  }
  const char *_history_file = history_file.has_value() ? history_file->c_str() : nullptr;
  fortran_read_a_line(
      /* const char* */ _prompt,
      /* const char* */ _line_out,
      /* bool* */ _trim_prompt,
      /* const char* */ _prompt_color,
      /* bool* */ _prompt_bold,
      /* const char* */ _history_file
  );
  return _line_out;
}
int SimUtils::readline_read_history(std::string history_file) {
  auto _history_file = history_file.c_str();
  int _status{};
  fortran_readline_read_history(/* const char* */ _history_file, /* int& */ _status);
  return _status;
}
int SimUtils::readline_write_history(std::string history_file) {
  auto _history_file = history_file.c_str();
  int _status{};
  fortran_readline_write_history(/* const char* */ _history_file, /* int& */ _status);
  return _status;
}
void SimUtils::real_num_fortran_format(
    double number,
    int width,
    std::string fmt_str,
    std::optional<int> n_blanks
) {
  int n_blanks_lvalue;
  auto *_n_blanks{&n_blanks_lvalue};
  if (n_blanks.has_value()) {
    n_blanks_lvalue = n_blanks.value();
  } else {
    _n_blanks = nullptr;
  }
  auto _fmt_str = fmt_str.c_str();
  fortran_real_num_fortran_format(
      /* double& */ number,
      /* int& */ width,
      /* int* */ _n_blanks,
      /* const char* */ _fmt_str
  );
}
void SimUtils::real_path(std::string path_in, std::string path_out, bool is_ok) {
  auto _path_in = path_in.c_str();
  auto _path_out = path_out.c_str();
  fortran_real_path(/* const char* */ _path_in, /* const char* */ _path_out, /* bool& */ is_ok);
}
void SimUtils::real_str(
    double r_num,
    std::string str,
    std::optional<int> n_signif,
    std::optional<int> n_decimal
) {
  int n_signif_lvalue;
  auto *_n_signif{&n_signif_lvalue};
  if (n_signif.has_value()) {
    n_signif_lvalue = n_signif.value();
  } else {
    _n_signif = nullptr;
  }
  int n_decimal_lvalue;
  auto *_n_decimal{&n_decimal_lvalue};
  if (n_decimal.has_value()) {
    n_decimal_lvalue = n_decimal.value();
  } else {
    _n_decimal = nullptr;
  }
  auto _str = str.c_str();
  fortran_real_str(
      /* double& */ r_num,
      /* int* */ _n_signif,
      /* int* */ _n_decimal,
      /* const char* */ _str
  );
}
void SimUtils::real_to_string(
    double real_num,
    int width,
    std::string str,
    std::optional<int> n_signif,
    std::optional<int> n_decimal
) {
  int n_signif_lvalue;
  auto *_n_signif{&n_signif_lvalue};
  if (n_signif.has_value()) {
    n_signif_lvalue = n_signif.value();
  } else {
    _n_signif = nullptr;
  }
  int n_decimal_lvalue;
  auto *_n_decimal{&n_decimal_lvalue};
  if (n_decimal.has_value()) {
    n_decimal_lvalue = n_decimal.value();
  } else {
    _n_decimal = nullptr;
  }
  auto _str = str.c_str();
  fortran_real_to_string(
      /* double& */ real_num,
      /* int& */ width,
      /* int* */ _n_signif,
      /* int* */ _n_decimal,
      /* const char* */ _str
  );
}
void SimUtils::reallocate_spline(
    SplineStructAlloc1D spline,
    int n,
    std::optional<int> n_min,
    std::optional<bool> exact
) {
  // intent=inout allocatable type array
  int n_min_lvalue;
  auto *_n_min{&n_min_lvalue};
  if (n_min.has_value()) {
    n_min_lvalue = n_min.value();
  } else {
    _n_min = nullptr;
  }
  bool exact_lvalue;
  auto *_exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_reallocate_spline(
      /* void* */ spline.get_fortran_ptr(),
      /* int& */ n,
      /* int* */ _n_min,
      /* bool* */ _exact
  );
}
void SimUtils::relbd(double phi, double phic, double mc, double b, double d) {
  fortran_relbd(
      /* double& */ phi,
      /* double& */ phic,
      /* double& */ mc,
      /* double& */ b,
      /* double& */ d
  );
}
void SimUtils::relcbd(double c0, double mc, double b, double dx) {
  fortran_relcbd(/* double& */ c0, /* double& */ mc, /* double& */ b, /* double& */ dx);
}
void SimUtils::relsbd(double s0, double mc, double b, double d) {
  fortran_relsbd(/* double& */ s0, /* double& */ mc, /* double& */ b, /* double& */ d);
}
void SimUtils::rgelbd(double phi, double mc, double elb, double eld) {
  fortran_rgelbd(/* double& */ phi, /* double& */ mc, /* double& */ elb, /* double& */ eld);
}
SimUtils::RmsValue
SimUtils::rms_value(FArray1D<Real> &val_arr, optional_ref<BoolAlloc1D> good_val) {
  // val_arr: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _val_arr_desc;
  _val_arr_desc.rank = 1;
  _val_arr_desc.data_ptr = val_arr.data();
  _val_arr_desc.dims[0] = val_arr.size();
  // intent=in allocatable general array
  auto *_good_val =
      good_val.has_value() ? good_val->get().get_fortran_ptr() : nullptr; // input, optional
  double _ave_val{};
  double _rms_val{};
  fortran_rms_value(
      /* Bmad::array_descriptor_t& */ _val_arr_desc,
      /* void* */ _good_val,
      /* double& */ _ave_val,
      /* double& */ _rms_val
  );
  return RmsValue{_ave_val, _rms_val};
}
FixedArray1D<Real, 2> SimUtils::rot_2d(FixedArray1D<Real, 2> vec_in, double angle) {
  // vec_in: in NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _vec_in_desc;
  _vec_in_desc.rank = 1;
  _vec_in_desc.data_ptr = vec_in.data();
  _vec_in_desc.dims[0] = vec_in.size();
  // vec_out: out NOT (CppWrapperGeneralArgumentArray) (['2'])
  Bmad::array_descriptor_t _vec_out_desc;
  _vec_out_desc.rank = 1;
  FixedArray1D<Real, 2> _vec_out;
  _vec_out_desc.data_ptr = _vec_out.data();
  _vec_out_desc.dims[0] = _vec_out.size();
  fortran_rot_2d(
      /* Bmad::array_descriptor_t& */ _vec_in_desc,
      /* double& */ angle,
      /* Bmad::array_descriptor_t& */ _vec_out_desc
  );
  return _vec_out;
}
void SimUtils::rotate_vec(FArray1D<Real> &vec, int axis, double angle) {
  // vec: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _vec_desc;
  _vec_desc.rank = 1;
  _vec_desc.data_ptr = vec.data();
  _vec_desc.dims[0] = vec.size();
  fortran_rotate_vec(
      /* Bmad::array_descriptor_t& */ _vec_desc,
      /* int& */ axis,
      /* double& */ angle
  );
}
FixedArray1D<Real, 3> SimUtils::rotate_vec_given_axis_angle(
    FixedArray1D<Real, 3> vec_in,
    FArray1D<Real> &axis,
    double angle
) {
  // vec_in: in NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_in_desc;
  _vec_in_desc.rank = 1;
  _vec_in_desc.data_ptr = vec_in.data();
  _vec_in_desc.dims[0] = vec_in.size();
  // axis: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _axis_desc;
  _axis_desc.rank = 1;
  _axis_desc.data_ptr = axis.data();
  _axis_desc.dims[0] = axis.size();
  // vec_out: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _vec_out_desc;
  _vec_out_desc.rank = 1;
  FixedArray1D<Real, 3> _vec_out;
  _vec_out_desc.data_ptr = _vec_out.data();
  _vec_out_desc.dims[0] = _vec_out.size();
  fortran_rotate_vec_given_axis_angle(
      /* Bmad::array_descriptor_t& */ _vec_in_desc,
      /* Bmad::array_descriptor_t& */ _axis_desc,
      /* double& */ angle,
      /* Bmad::array_descriptor_t& */ _vec_out_desc
  );
  return _vec_out;
}
double SimUtils::rp8(int int_in) {
  double _re_out{};
  fortran_rp8(/* int& */ int_in, /* double& */ _re_out);
  return _re_out;
}
void SimUtils::rserbd(double y, double m, double b, double d) {
  fortran_rserbd(/* double& */ y, /* double& */ m, /* double& */ b, /* double& */ d);
}
void SimUtils::run_timer(
    std::string command,
    std::optional<double> time,
    std::optional<double> time0
) {
  auto _command = command.c_str();
  double time_lvalue;
  auto *_time{&time_lvalue};
  if (time.has_value()) {
    time_lvalue = time.value();
  } else {
    _time = nullptr;
  }
  double time0_lvalue;
  auto *_time0{&time0_lvalue};
  if (time0.has_value()) {
    time0_lvalue = time0.value();
  } else {
    _time0 = nullptr;
  }
  fortran_run_timer(/* const char* */ _command, /* double* */ _time, /* double* */ _time0);
}
SimUtils::Serbd SimUtils::serbd(double y, double m) {
  double _b{};
  double _d{};
  fortran_serbd(/* double& */ y, /* double& */ m, /* double& */ _b, /* double& */ _d);
  return Serbd{_b, _d};
}
void SimUtils::set_env(std::string env_name, std::string env_value, bool err_flag) {
  auto _env_name = env_name.c_str();
  auto _env_value = env_value.c_str();
  fortran_set_env(/* const char* */ _env_name, /* const char* */ _env_value, /* bool& */ err_flag);
}
void SimUtils::set_parameter(int param_val, int set_val, int save_val) {
  fortran_set_parameter_int(/* int& */ param_val, /* int& */ set_val, /* int& */ save_val);
}
void SimUtils::set_parameter(bool param_val, bool set_val, bool save_val) {
  fortran_set_parameter_logic(/* bool& */ param_val, /* bool& */ set_val, /* bool& */ save_val);
}
void SimUtils::set_parameter(double param_val, double set_val, double save_val) {
  fortran_set_parameter_real(
      /* double& */ param_val,
      /* double& */ set_val,
      /* double& */ save_val
  );
}
int SimUtils::set_species_charge(int species_in, int charge) {
  int _species_charged{};
  fortran_set_species_charge(/* int& */ species_in, /* int& */ charge, /* int& */ _species_charged);
  return _species_charged;
}
int SimUtils::sign_of(int num, std::optional<bool> zero_is_zero) {
  bool zero_is_zero_lvalue;
  auto *_zero_is_zero{&zero_is_zero_lvalue};
  if (zero_is_zero.has_value()) {
    zero_is_zero_lvalue = zero_is_zero.value();
  } else {
    _zero_is_zero = nullptr;
  }
  int _num_sign{};
  fortran_sign_of_int(/* int& */ num, /* bool* */ _zero_is_zero, /* int& */ _num_sign);
  return _num_sign;
}
int SimUtils::sign_of(double num, std::optional<bool> zero_is_zero) {
  bool zero_is_zero_lvalue;
  auto *_zero_is_zero{&zero_is_zero_lvalue};
  if (zero_is_zero.has_value()) {
    zero_is_zero_lvalue = zero_is_zero.value();
  } else {
    _zero_is_zero = nullptr;
  }
  int _num_sign{};
  fortran_sign_of_real(/* double& */ num, /* bool* */ _zero_is_zero, /* int& */ _num_sign);
  return _num_sign;
}
double SimUtils::sinc(double x, std::optional<int> nd) {
  int nd_lvalue;
  auto *_nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  double _y{};
  fortran_sinc(/* double& */ x, /* int* */ _nd, /* double& */ _y);
  return _y;
}
double SimUtils::sincc(double x, std::optional<int> nd) {
  int nd_lvalue;
  auto *_nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  double _y{};
  fortran_sincc(/* double& */ x, /* int* */ _nd, /* double& */ _y);
  return _y;
}
double SimUtils::sinhx_x(double x, std::optional<int> nd) {
  int nd_lvalue;
  auto *_nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  double _y{};
  fortran_sinhx_x(/* double& */ x, /* int* */ _nd, /* double& */ _y);
  return _y;
}
void SimUtils::skip_header(int ix_unit, bool error_flag) {
  fortran_skip_header(/* int& */ ix_unit, /* bool& */ error_flag);
}
void SimUtils::special_projection(double f, double func_retval__, std::optional<int> status) {
  int status_lvalue;
  auto *_status{&status_lvalue};
  if (status.has_value()) {
    status_lvalue = status.value();
  } else {
    _status = nullptr;
  }
  fortran_special_projection(/* double& */ f, /* int* */ _status, /* double& */ func_retval__);
}
int SimUtils::species_id(
    std::string name,
    std::optional<int> default_,
    std::optional<bool> print_err
) {
  auto _name = name.c_str();
  int default__lvalue;
  auto *_default_{&default__lvalue};
  if (default_.has_value()) {
    default__lvalue = default_.value();
  } else {
    _default_ = nullptr;
  }
  bool print_err_lvalue;
  auto *_print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _species{};
  fortran_species_id(
      /* const char* */ _name,
      /* int* */ _default_,
      /* bool* */ _print_err,
      /* int& */ _species
  );
  return _species;
}
int SimUtils::species_id_from_openpmd(std::string pmd_name, int charge) {
  auto _pmd_name = pmd_name.c_str();
  int _species{};
  fortran_species_id_from_openpmd(
      /* const char* */ _pmd_name,
      /* int& */ charge,
      /* int& */ _species
  );
  return _species;
}
std::string SimUtils::species_name(int species) {
  char _name[4096];
  fortran_species_name(/* int& */ species, /* const char* */ _name);
  return _name;
}
int SimUtils::species_of(double mass, int charge) {
  int _species{};
  fortran_species_of(/* double& */ mass, /* int& */ charge, /* int& */ _species);
  return _species;
}
double SimUtils::spin_of(int species, std::optional<double> non_subatomic_default) {
  double non_subatomic_default_lvalue;
  auto *_non_subatomic_default{&non_subatomic_default_lvalue};
  if (non_subatomic_default.has_value()) {
    non_subatomic_default_lvalue = non_subatomic_default.value();
  } else {
    _non_subatomic_default = nullptr;
  }
  double _spin{};
  fortran_spin_of(/* int& */ species, /* double* */ _non_subatomic_default, /* double& */ _spin);
  return _spin;
}
double SimUtils::spline1(SplineStruct &a_spline, double x, std::optional<int> n) {
  int n_lvalue;
  auto *_n{&n_lvalue};
  if (n.has_value()) {
    n_lvalue = n.value();
  } else {
    _n = nullptr;
  }
  double _y{};
  fortran_spline1(
      /* void* */ a_spline.get_fortran_ptr(),
      /* double& */ x,
      /* int* */ _n,
      /* double& */ _y
  );
  return _y;
}
bool SimUtils::spline_akima(SplineStructArray1D spline) {
  // spline: SplineStruct inout (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spline_desc;
  _spline_desc.rank = 1;
  _spline_desc.data_ptr = spline.data();
  _spline_desc.dims[0] = spline.size();
  _spline_desc.strides[0] = 1;
  bool _ok{};
  fortran_spline_akima(/* Bmad::array_descriptor_t& */ _spline_desc, /* bool& */ _ok);
  return _ok;
}
SimUtils::SplineAkimaInterpolate
SimUtils::spline_akima_interpolate(FArray1D<Real> &x_knot, FArray1D<Real> &y_knot, double x) {
  // x_knot: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_knot_desc;
  _x_knot_desc.rank = 1;
  _x_knot_desc.data_ptr = x_knot.data();
  _x_knot_desc.dims[0] = x_knot.size();
  // y_knot: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _y_knot_desc;
  _y_knot_desc.rank = 1;
  _y_knot_desc.data_ptr = y_knot.data();
  _y_knot_desc.dims[0] = y_knot.size();
  bool _ok{};
  double _y{};
  double _dy{};
  fortran_spline_akima_interpolate(
      /* Bmad::array_descriptor_t& */ _x_knot_desc,
      /* Bmad::array_descriptor_t& */ _y_knot_desc,
      /* double& */ x,
      /* bool& */ _ok,
      /* double& */ _y,
      /* double& */ _dy
  );
  return SplineAkimaInterpolate{_ok, _y, _dy};
}
SimUtils::SplineEvaluate SimUtils::spline_evaluate(SplineStructArray1D spline, double x) {
  // spline: SplineStruct in (CppWrapperTypeArgumentArray)
  Bmad::array_descriptor_t _spline_desc;
  _spline_desc.rank = 1;
  _spline_desc.data_ptr = spline.data();
  _spline_desc.dims[0] = spline.size();
  _spline_desc.strides[0] = 1;
  bool _ok{};
  double _y{};
  double _dy{};
  fortran_spline_evaluate(
      /* Bmad::array_descriptor_t& */ _spline_desc,
      /* double& */ x,
      /* bool& */ _ok,
      /* double& */ _y,
      /* double& */ _dy
  );
  return SplineEvaluate{_ok, _y, _dy};
}
double SimUtils::sqrt_alpha(double alpha, double x) {
  double _y{};
  fortran_sqrt_alpha(/* double& */ alpha, /* double& */ x, /* double& */ _y);
  return _y;
}
void SimUtils::sqrt_one(double x, double ds1, std::optional<int> nd) {
  int nd_lvalue;
  auto *_nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_sqrt_one(/* double& */ x, /* int* */ _nd, /* double& */ ds1);
}
void SimUtils::str_count(std::string str, std::string match, int num) {
  auto _str = str.c_str();
  auto _match = match.c_str();
  fortran_str_count(/* const char* */ _str, /* const char* */ _match, /* int& */ num);
}
std::string SimUtils::str_downcase(std::string src) {
  char _dst[4096];
  auto _src = src.c_str();
  fortran_str_downcase(/* const char* */ _dst, /* const char* */ _src);
  return _dst;
}
void SimUtils::str_first_in_set(
    std::string line,
    std::string set,
    int ix_match,
    std::optional<bool> ignore_clauses
) {
  auto _line = line.c_str();
  auto _set = set.c_str();
  bool ignore_clauses_lvalue;
  auto *_ignore_clauses{&ignore_clauses_lvalue};
  if (ignore_clauses.has_value()) {
    ignore_clauses_lvalue = ignore_clauses.value();
  } else {
    _ignore_clauses = nullptr;
  }
  fortran_str_first_in_set(
      /* const char* */ _line,
      /* const char* */ _set,
      /* bool* */ _ignore_clauses,
      /* int& */ ix_match
  );
}
void SimUtils::str_first_not_in_set(std::string line, std::string set, int ix_match) {
  auto _line = line.c_str();
  auto _set = set.c_str();
  fortran_str_first_not_in_set(
      /* const char* */ _line,
      /* const char* */ _set,
      /* int& */ ix_match
  );
}
void SimUtils::str_last_in_set(std::string line, std::string set, int ix_match) {
  auto _line = line.c_str();
  auto _set = set.c_str();
  fortran_str_last_in_set(/* const char* */ _line, /* const char* */ _set, /* int& */ ix_match);
}
void SimUtils::str_last_not_in_set(std::string line, std::string set, int ix_match) {
  auto _line = line.c_str();
  auto _set = set.c_str();
  fortran_str_last_not_in_set(/* const char* */ _line, /* const char* */ _set, /* int& */ ix_match);
}
void SimUtils::str_match_wild(std::string str, std::string pat, bool a_match) {
  auto _str = str.c_str();
  auto _pat = pat.c_str();
  fortran_str_match_wild(/* const char* */ _str, /* const char* */ _pat, /* bool& */ a_match);
}
void SimUtils::str_substitute(
    std::string string,
    std::optional<std::string> str_match,
    std::optional<std::string> str_replace,
    std::optional<bool> do_trim,
    std::optional<bool> ignore_escaped
) {
  auto _string = string.c_str();
  const char *_str_match = str_match.has_value() ? str_match->c_str() : nullptr;
  const char *_str_replace = str_replace.has_value() ? str_replace->c_str() : nullptr;
  bool do_trim_lvalue;
  auto *_do_trim{&do_trim_lvalue};
  if (do_trim.has_value()) {
    do_trim_lvalue = do_trim.value();
  } else {
    _do_trim = nullptr;
  }
  bool ignore_escaped_lvalue;
  auto *_ignore_escaped{&ignore_escaped_lvalue};
  if (ignore_escaped.has_value()) {
    ignore_escaped_lvalue = ignore_escaped.value();
  } else {
    _ignore_escaped = nullptr;
  }
  fortran_str_substitute(
      /* const char* */ _string,
      /* const char* */ _str_match,
      /* const char* */ _str_replace,
      /* bool* */ _do_trim,
      /* bool* */ _ignore_escaped
  );
}
std::string SimUtils::str_upcase(std::string src) {
  char _dst[4096];
  auto _src = src.c_str();
  fortran_str_upcase(/* const char* */ _dst, /* const char* */ _src);
  return _dst;
}
void SimUtils::string_to_int(
    std::string line,
    int default_,
    bool err_flag,
    int value,
    std::optional<bool> err_print_flag
) {
  auto _line = line.c_str();
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  fortran_string_to_int(
      /* const char* */ _line,
      /* int& */ default_,
      /* bool& */ err_flag,
      /* bool* */ _err_print_flag,
      /* int& */ value
  );
}
void SimUtils::string_to_real(
    std::string line,
    double default_,
    bool err_flag,
    double value,
    std::optional<bool> err_print_flag
) {
  auto _line = line.c_str();
  bool err_print_flag_lvalue;
  auto *_err_print_flag{&err_print_flag_lvalue};
  if (err_print_flag.has_value()) {
    err_print_flag_lvalue = err_print_flag.value();
  } else {
    _err_print_flag = nullptr;
  }
  fortran_string_to_real(
      /* const char* */ _line,
      /* double& */ default_,
      /* bool& */ err_flag,
      /* bool* */ _err_print_flag,
      /* double& */ value
  );
}
void SimUtils::string_trim(std::string in_string, std::string out_string, int word_len) {
  auto _in_string = in_string.c_str();
  auto _out_string = out_string.c_str();
  fortran_string_trim(
      /* const char* */ _in_string,
      /* const char* */ _out_string,
      /* int& */ word_len
  );
}
void SimUtils::string_trim2(
    std::string in_str,
    std::string delimitors,
    std::string out_str,
    int ix_word,
    std::string delim,
    int ix_next
) {
  auto _in_str = in_str.c_str();
  auto _delimitors = delimitors.c_str();
  auto _out_str = out_str.c_str();
  auto _delim = delim.c_str();
  fortran_string_trim2(
      /* const char* */ _in_str,
      /* const char* */ _delimitors,
      /* const char* */ _out_str,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* int& */ ix_next
  );
}
bool SimUtils::suggest_lmdif(
    FArray1D<Real> &XV,
    FArray1D<Real> &FV,
    double eps,
    int itermx,
    std::optional<bool> reset_flag
) {
  // XV: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _XV_desc;
  _XV_desc.rank = 1;
  _XV_desc.data_ptr = XV.data();
  _XV_desc.dims[0] = XV.size();
  // FV: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _FV_desc;
  _FV_desc.rank = 1;
  _FV_desc.data_ptr = FV.data();
  _FV_desc.dims[0] = FV.size();
  bool _at_end{};
  bool reset_flag_lvalue;
  auto *_reset_flag{&reset_flag_lvalue};
  if (reset_flag.has_value()) {
    reset_flag_lvalue = reset_flag.value();
  } else {
    _reset_flag = nullptr;
  }
  fortran_suggest_lmdif(
      /* Bmad::array_descriptor_t& */ _XV_desc,
      /* Bmad::array_descriptor_t& */ _FV_desc,
      /* double& */ eps,
      /* int& */ itermx,
      /* bool& */ _at_end,
      /* bool* */ _reset_flag
  );
  return _at_end;
}
FixedArray2D<Real, 4, 4> SimUtils::super_bicubic_coef(
    FixedArray1D<Real, 4> y,
    FixedArray1D<Real, 4> y1,
    FixedArray1D<Real, 4> y2,
    FixedArray1D<Real, 4> y12,
    double d1,
    double d2
) {
  // y: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y_desc;
  _y_desc.rank = 1;
  _y_desc.data_ptr = y.data();
  _y_desc.dims[0] = y.size();
  // y1: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y1_desc;
  _y1_desc.rank = 1;
  _y1_desc.data_ptr = y1.data();
  _y1_desc.dims[0] = y1.size();
  // y2: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y2_desc;
  _y2_desc.rank = 1;
  _y2_desc.data_ptr = y2.data();
  _y2_desc.dims[0] = y2.size();
  // y12: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y12_desc;
  _y12_desc.rank = 1;
  _y12_desc.data_ptr = y12.data();
  _y12_desc.dims[0] = y12.size();
  // c: out NOT (CppWrapperGeneralArgumentArray) (['4', '4'])
  Bmad::array_descriptor_t _c_desc;
  _c_desc.rank = 2;
  FixedArray2D<Real, 4, 4> c;
  double _c_vec[4 * 4];
  _c_desc.data_ptr = _c_vec;
  _c_desc.dims[0] = 4;
  _c_desc.dims[1] = 4;
  fortran_super_bicubic_coef(
      /* Bmad::array_descriptor_t& */ _y_desc,
      /* Bmad::array_descriptor_t& */ _y1_desc,
      /* Bmad::array_descriptor_t& */ _y2_desc,
      /* Bmad::array_descriptor_t& */ _y12_desc,
      /* double& */ d1,
      /* double& */ d2,
      /* Bmad::array_descriptor_t& */ _c_desc
  );
  vec_to_matrix(_c_vec, c);
  return c;
}
SimUtils::SuperBicubicInterpolation SimUtils::super_bicubic_interpolation(
    FixedArray1D<Real, 4> y,
    FixedArray1D<Real, 4> y1,
    FixedArray1D<Real, 4> y2,
    FixedArray1D<Real, 4> y12,
    double x1l,
    double x1u,
    double x2l,
    double x2u,
    double x1,
    double x2
) {
  // y: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y_desc;
  _y_desc.rank = 1;
  _y_desc.data_ptr = y.data();
  _y_desc.dims[0] = y.size();
  // y1: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y1_desc;
  _y1_desc.rank = 1;
  _y1_desc.data_ptr = y1.data();
  _y1_desc.dims[0] = y1.size();
  // y2: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y2_desc;
  _y2_desc.rank = 1;
  _y2_desc.data_ptr = y2.data();
  _y2_desc.dims[0] = y2.size();
  // y12: in NOT (CppWrapperGeneralArgumentArray) (['4'])
  Bmad::array_descriptor_t _y12_desc;
  _y12_desc.rank = 1;
  _y12_desc.data_ptr = y12.data();
  _y12_desc.dims[0] = y12.size();
  double _ansy{};
  double _ansy1{};
  double _ansy2{};
  fortran_super_bicubic_interpolation(
      /* Bmad::array_descriptor_t& */ _y_desc,
      /* Bmad::array_descriptor_t& */ _y1_desc,
      /* Bmad::array_descriptor_t& */ _y2_desc,
      /* Bmad::array_descriptor_t& */ _y12_desc,
      /* double& */ x1l,
      /* double& */ x1u,
      /* double& */ x2l,
      /* double& */ x2u,
      /* double& */ x1,
      /* double& */ x2,
      /* double& */ _ansy,
      /* double& */ _ansy1,
      /* double& */ _ansy2
  );
  return SuperBicubicInterpolation{_ansy, _ansy1, _ansy2};
}
SimUtils::SuperPolint SimUtils::super_polint(FArray1D<Real> &xa, FArray1D<Real> &ya, double x) {
  // xa: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _xa_desc;
  _xa_desc.rank = 1;
  _xa_desc.data_ptr = xa.data();
  _xa_desc.dims[0] = xa.size();
  // ya: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _ya_desc;
  _ya_desc.rank = 1;
  _ya_desc.data_ptr = ya.data();
  _ya_desc.dims[0] = ya.size();
  double _y{};
  double _dy{};
  fortran_super_polint(
      /* Bmad::array_descriptor_t& */ _xa_desc,
      /* Bmad::array_descriptor_t& */ _ya_desc,
      /* double& */ x,
      /* double& */ _y,
      /* double& */ _dy
  );
  return SuperPolint{_y, _dy};
}
double SimUtils::super_poly(double x, FArray1D<Real> &coeffs) {
  // coeffs: in NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _coeffs_desc;
  _coeffs_desc.rank = 1;
  _coeffs_desc.data_ptr = coeffs.data();
  _coeffs_desc.dims[0] = coeffs.size();
  double _value{};
  fortran_super_poly(
      /* double& */ x,
      /* Bmad::array_descriptor_t& */ _coeffs_desc,
      /* double& */ _value
  );
  return _value;
}
void SimUtils::super_sobseq(FArray1D<Real> &x, optional_ref<RandomStateStruct> ran_state) {
  // x: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _x_desc;
  _x_desc.rank = 1;
  _x_desc.data_ptr = x.data();
  _x_desc.dims[0] = x.size();
  auto *_ran_state =
      ran_state.has_value() ? ran_state->get().get_fortran_ptr() : nullptr; // input, optional
  fortran_super_sobseq(/* Bmad::array_descriptor_t& */ _x_desc, /* void* */ _ran_state);
}
void SimUtils::super_sort(FArray1D<Int> &arr) {
  // arr: inout NOT (CppWrapperGeneralArgumentArray) ([':'])
  Bmad::array_descriptor_t _arr_desc;
  _arr_desc.rank = 1;
  _arr_desc.data_ptr = arr.data();
  _arr_desc.dims[0] = arr.size();
  fortran_super_sort(/* Bmad::array_descriptor_t& */ _arr_desc);
}
bool SimUtils::system_command(std::string line) {
  auto _line = line.c_str();
  bool _err_flag{};
  fortran_system_command(/* const char* */ _line, /* bool& */ _err_flag);
  return _err_flag;
}
void SimUtils::test_xgelbd() { fortran_test_xgelbd(); }
void SimUtils::to_str(double num, std::string string, std::optional<int> max_signif) {
  int max_signif_lvalue;
  auto *_max_signif{&max_signif_lvalue};
  if (max_signif.has_value()) {
    max_signif_lvalue = max_signif.value();
  } else {
    _max_signif = nullptr;
  }
  auto _string = string.c_str();
  fortran_to_str(/* double& */ num, /* int* */ _max_signif, /* const char* */ _string);
}
SimUtils::TricubicCmplxEval SimUtils::tricubic_cmplx_eval(
    double x_norm,
    double y_norm,
    double z_norm,
    TricubicCmplxCoefStruct &tri_coef
) {
  std::complex<double> _df_dx{};
  std::complex<double> _df_dy{};
  std::complex<double> _df_dz{};
  std::complex<double> _f_val{};
  fortran_tricubic_cmplx_eval(
      /* double& */ x_norm,
      /* double& */ y_norm,
      /* double& */ z_norm,
      /* void* */ tri_coef.get_fortran_ptr(),
      /* std::complex<double>& */ _df_dx,
      /* std::complex<double>& */ _df_dy,
      /* std::complex<double>& */ _df_dz,
      /* std::complex<double>& */ _f_val
  );
  return TricubicCmplxEval{_df_dx, _df_dy, _df_dz, _f_val};
}
void SimUtils::type_this_file(std::string filename) {
  auto _filename = filename.c_str();
  fortran_type_this_file(/* const char* */ _filename);
}
void SimUtils::upcase_string(std::string string) {
  auto _string = string.c_str();
  fortran_upcase_string(/* const char* */ _string);
}
void SimUtils::virtual_memory_usage(int usage) { fortran_virtual_memory_usage(/* int& */ usage); }
SimUtils::WMatToAxisAngle SimUtils::w_mat_to_axis_angle(FixedArray2D<Real, 3, 3> w_mat) {
  // w_mat: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  matrix_to_vec(w_mat, _w_mat_vec);
  // axis: out NOT (CppWrapperGeneralArgumentArray) (['3'])
  Bmad::array_descriptor_t _axis_desc;
  _axis_desc.rank = 1;
  FixedArray1D<Real, 3> _axis;
  _axis_desc.data_ptr = _axis.data();
  _axis_desc.dims[0] = _axis.size();
  double _angle{};
  fortran_w_mat_to_axis_angle(
      /* Bmad::array_descriptor_t& */ _w_mat_desc,
      /* Bmad::array_descriptor_t& */ _axis_desc,
      /* double& */ _angle
  );
  return WMatToAxisAngle{_axis, _angle};
}
FixedArray1D<Real, 4> SimUtils::w_mat_to_quat(FixedArray2D<Real, 3, 3> w_mat) {
  // w_mat: in NOT (CppWrapperGeneralArgumentArray) (['3', '3'])
  Bmad::array_descriptor_t _w_mat_desc;
  _w_mat_desc.rank = 2;
  double _w_mat_vec[3 * 3];
  _w_mat_desc.data_ptr = _w_mat_vec;
  _w_mat_desc.dims[0] = 3;
  _w_mat_desc.dims[1] = 3;
  matrix_to_vec(w_mat, _w_mat_vec);
  // quat: out NOT (CppWrapperGeneralArgumentArray) (['0:3'])
  Bmad::array_descriptor_t _quat_desc;
  _quat_desc.rank = 1;
  FixedArray1D<Real, 4> _quat;
  _quat_desc.data_ptr = _quat.data();
  _quat_desc.dims[0] = _quat.size();
  fortran_w_mat_to_quat(
      /* Bmad::array_descriptor_t& */ _w_mat_desc,
      /* Bmad::array_descriptor_t& */ _quat_desc
  );
  return _quat;
}
void SimUtils::word_len(std::string wording, int wlen) {
  auto _wording = wording.c_str();
  fortran_word_len(/* const char* */ _wording, /* int& */ wlen);
}
void SimUtils::word_read(
    std::string in_str,
    std::string delim_list,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    std::string out_str,
    std::optional<bool> ignore_interior
) {
  auto _in_str = in_str.c_str();
  auto _delim_list = delim_list.c_str();
  auto _word = word.c_str();
  auto _delim = delim.c_str();
  auto _out_str = out_str.c_str();
  bool ignore_interior_lvalue;
  auto *_ignore_interior{&ignore_interior_lvalue};
  if (ignore_interior.has_value()) {
    ignore_interior_lvalue = ignore_interior.value();
  } else {
    _ignore_interior = nullptr;
  }
  fortran_word_read(
      /* const char* */ _in_str,
      /* const char* */ _delim_list,
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* const char* */ _out_str,
      /* bool* */ _ignore_interior
  );
}
double SimUtils::x0_radiation_length(int species) {
  double _x0{};
  fortran_x0_radiation_length(/* int& */ species, /* double& */ _x0);
  return _x0;
}
void SimUtils::zig_table_init() { fortran_zig_table_init(); }