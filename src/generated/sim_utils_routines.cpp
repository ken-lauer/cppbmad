#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/sim_utils_routines.hpp"
#include "bmad/types.h"
#include "json.hpp"

using namespace Bmad;

using json = nlohmann::json;
int SimUtils::antiparticle(int species) {
  int _anti_species{};
  fortran_antiparticle(/* c_Int& */ species, /* c_Int& */ _anti_species);
  return _anti_species;
}
int SimUtils::species_of(double mass, int charge) {
  int _species{};
  fortran_species_of(
      /* c_Real& */ mass, /* c_Int& */ charge, /* c_Int& */ _species);
  return _species;
}
int SimUtils::species_id(
    std::string name,
    std::optional<int> default_,
    std::optional<bool> print_err) {
  auto _name = name.c_str();
  int default__lvalue;
  auto* _default_{&default__lvalue};
  if (default_.has_value()) {
    default__lvalue = default_.value();
  } else {
    _default_ = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  int _species{};
  fortran_species_id(
      /* c_Char */ _name,
      /* c_Int* */ _default_,
      /* c_Bool* */ _print_err,
      /* c_Int& */ _species);
  return _species;
}
int SimUtils::atomic_species_id(
    int charge,
    bool is_anti,
    int atomic_num,
    int n_nuc) {
  int _species_id{};
  fortran_atomic_species_id(
      /* c_Int& */ charge,
      /* c_Bool& */ is_anti,
      /* c_Int& */ atomic_num,
      /* c_Int& */ n_nuc,
      /* c_Int& */ _species_id);
  return _species_id;
}
std::string SimUtils::species_name(int species) {
  char _name[4096];
  fortran_species_name(/* c_Int& */ species, /* c_Char */ _name);
  return _name;
}
int SimUtils::species_id_from_openpmd(std::string pmd_name, int charge) {
  auto _pmd_name = pmd_name.c_str();
  int _species{};
  fortran_species_id_from_openpmd(
      /* c_Char */ _pmd_name, /* c_Int& */ charge, /* c_Int& */ _species);
  return _species;
}
std::string SimUtils::openpmd_species_name(int species) {
  char _pmd_name[4096];
  fortran_openpmd_species_name(/* c_Int& */ species, /* c_Char */ _pmd_name);
  return _pmd_name;
}
double SimUtils::anomalous_moment_of(int species) {
  double _moment{};
  fortran_anomalous_moment_of(/* c_Int& */ species, /* c_Real& */ _moment);
  return _moment;
}
double SimUtils::spin_of(
    int species,
    std::optional<double> non_subatomic_default) {
  double non_subatomic_default_lvalue;
  auto* _non_subatomic_default{&non_subatomic_default_lvalue};
  if (non_subatomic_default.has_value()) {
    non_subatomic_default_lvalue = non_subatomic_default.value();
  } else {
    _non_subatomic_default = nullptr;
  }
  double _spin{};
  fortran_spin_of(
      /* c_Int& */ species,
      /* c_Real* */ _non_subatomic_default,
      /* c_Real& */ _spin);
  return _spin;
}
int SimUtils::charge_of(int species, std::optional<int> default_) {
  int default__lvalue;
  auto* _default_{&default__lvalue};
  if (default_.has_value()) {
    default__lvalue = default_.value();
  } else {
    _default_ = nullptr;
  }
  int _charge{};
  fortran_charge_of(
      /* c_Int& */ species, /* c_Int* */ _default_, /* c_Int& */ _charge);
  return _charge;
}
double SimUtils::mass_of(int species) {
  double _mass{};
  fortran_mass_of(/* c_Int& */ species, /* c_Real& */ _mass);
  return _mass;
}
double SimUtils::charge_to_mass_of(int species) {
  double _charge_mass_ratio{};
  fortran_charge_to_mass_of(
      /* c_Int& */ species, /* c_Real& */ _charge_mass_ratio);
  return _charge_mass_ratio;
}
int SimUtils::set_species_charge(int species_in, int charge) {
  int _species_charged{};
  fortran_set_species_charge(
      /* c_Int& */ species_in,
      /* c_Int& */ charge,
      /* c_Int& */ _species_charged);
  return _species_charged;
}
double SimUtils::x0_radiation_length(int species) {
  double _x0{};
  fortran_x0_radiation_length(/* c_Int& */ species, /* c_Real& */ _x0);
  return _x0;
}
int SimUtils::atomic_number(int species) {
  int _atomic_num{};
  fortran_atomic_number(/* c_Int& */ species, /* c_Int& */ _atomic_num);
  return _atomic_num;
}
bool SimUtils::is_subatomic_species(int species) {
  bool _is_subatomic{};
  fortran_is_subatomic_species(
      /* c_Int& */ species, /* c_Bool& */ _is_subatomic);
  return _is_subatomic;
}
bool SimUtils::is_true(double param) {
  bool _this_true{};
  fortran_is_true(/* c_Real& */ param, /* c_Bool& */ _this_true);
  return _this_true;
}
bool SimUtils::is_false(double param) {
  bool _this_false{};
  fortran_is_false(/* c_Real& */ param, /* c_Bool& */ _this_false);
  return _this_false;
}
double SimUtils::rp8(int int_in) {
  double _re_out{};
  fortran_rp8(/* c_Int& */ int_in, /* c_Real& */ _re_out);
  return _re_out;
}
void SimUtils::set_parameter_real(
    double param_val,
    double set_val,
    double save_val) {
  fortran_set_parameter_real(
      /* c_Real& */ param_val, /* c_Real& */ set_val, /* c_Real& */ save_val);
}
void SimUtils::set_parameter_int(int param_val, int set_val, int save_val) {
  fortran_set_parameter_int(
      /* c_Int& */ param_val, /* c_Int& */ set_val, /* c_Int& */ save_val);
}
void SimUtils::set_parameter_logic(
    bool param_val,
    bool set_val,
    bool save_val) {
  fortran_set_parameter_logic(
      /* c_Bool& */ param_val, /* c_Bool& */ set_val, /* c_Bool& */ save_val);
}
void SimUtils::cosc(double x, optional_ref<int> nd, double y) {
  auto* _nd = nd.has_value() ? &nd->get() : nullptr; // inout, optional
  fortran_cosc(/* c_Real& */ x, /* c_Int* */ _nd, /* c_Real& */ y);
}
void SimUtils::asinc(double x, optional_ref<int> nd, double y) {
  auto* _nd = nd.has_value() ? &nd->get() : nullptr; // inout, optional
  fortran_asinc(/* c_Real& */ x, /* c_Int* */ _nd, /* c_Real& */ y);
}
void SimUtils::assert_equal(
    IntAllocatable1D& int_arr,
    std::string err_str,
    int ival) {
  // intent=inout allocatable general array
  auto _err_str = err_str.c_str(); // ptr, inout, required
  fortran_assert_equal(
      /* void* */ int_arr.get_fortran_ptr(),
      /* c_Char */ _err_str,
      /* c_Int& */ ival);
}
void SimUtils::calc_file_number(
    std::string file_name,
    int num_in,
    int num_out,
    bool err_flag) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_calc_file_number(
      /* c_Char */ _file_name,
      /* c_Int& */ num_in,
      /* c_Int& */ num_out,
      /* c_Bool& */ err_flag);
}
void SimUtils::change_file_number(std::string file_name, int change) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_change_file_number(/* c_Char */ _file_name, /* c_Int& */ change);
}
void SimUtils::cos_one(double angle, double cos1) {
  fortran_cos_one(/* c_Real& */ angle, /* c_Real& */ cos1);
}
void SimUtils::complex_error_function(
    double wr,
    double wi,
    double zr,
    double zi) {
  fortran_complex_error_function(
      /* c_Real& */ wr, /* c_Real& */ wi, /* c_Real& */ zr, /* c_Real& */ zi);
}
void SimUtils::cross_product(
    RealAllocatable1D& a,
    RealAllocatable1D& b,
    FixedArray1D<Real, 3> c) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  auto* _c = c.data(); // CppWrapperGeneralArgument
  fortran_cross_product(
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* c_RealArr */ _c);
}
void SimUtils::date_and_time_stamp(
    std::string string,
    optional_ref<bool> numeric_month,
    optional_ref<bool> include_zone) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _numeric_month = numeric_month.has_value() ? &numeric_month->get()
                                                   : nullptr; // inout, optional
  auto* _include_zone = include_zone.has_value() ? &include_zone->get()
                                                 : nullptr; // inout, optional
  fortran_date_and_time_stamp(
      /* c_Char */ _string,
      /* c_Bool* */ _numeric_month,
      /* c_Bool* */ _include_zone);
}
void SimUtils::detab(std::string str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_detab(/* c_Char */ _str);
}
void SimUtils::display_size_and_resolution(
    int ix_screen,
    double x_size,
    double y_size,
    double x_res,
    double y_res) {
  fortran_display_size_and_resolution(
      /* c_Int& */ ix_screen,
      /* c_Real& */ x_size,
      /* c_Real& */ y_size,
      /* c_Real& */ x_res,
      /* c_Real& */ y_res);
}
void SimUtils::dj_bessel(int m, double arg, double dj_bes) {
  fortran_dj_bessel(/* c_Int& */ m, /* c_Real& */ arg, /* c_Real& */ dj_bes);
}
void SimUtils::djb_hash(std::string str, optional_ref<int> old_hash, int hash) {
  auto _str = str.c_str(); // ptr, inout, required
  auto* _old_hash =
      old_hash.has_value() ? &old_hash->get() : nullptr; // inout, optional
  fortran_djb_hash(
      /* c_Char */ _str, /* c_Int* */ _old_hash, /* c_Int& */ hash);
}
void SimUtils::djb_str_hash(std::string in_str, std::string hash_str) {
  auto _in_str = in_str.c_str(); // ptr, inout, required
  auto _hash_str = hash_str.c_str(); // ptr, inout, required
  fortran_djb_str_hash(/* c_Char */ _in_str, /* c_Char */ _hash_str);
}
void SimUtils::downcase_string(std::string string) {
  auto _string = string.c_str(); // ptr, inout, required
  fortran_downcase_string(/* c_Char */ _string);
}
void SimUtils::err_exit(optional_ref<std::string> err_str) {
  const char* _err_str = err_str.has_value() ? err_str->get().c_str() : nullptr;
  fortran_err_exit(/* c_Char */ _err_str);
}
void SimUtils::factorial(int n, double fact) {
  fortran_factorial(/* c_Int& */ n, /* c_Real& */ fact);
}
void SimUtils::faddeeva_function(
    FixedArray1D<Real, 2> z,
    FixedArray1D<Real, 2> w,
    FixedArray2D<Real, 2, 2> dw) {
  auto* _z = z.data(); // CppWrapperGeneralArgument
  auto* _w = w.data(); // CppWrapperGeneralArgument
  double _dw_vec[2 * 2];
  matrix_to_vec(dw, _dw_vec);
  fortran_faddeeva_function(
      /* c_RealArr */ _z, /* c_RealArr */ _w, /* c_RealArr */ _dw_vec);
  vec_to_matrix(_dw_vec, dw);
}
void SimUtils::fft_1d(ComplexAllocatable1D& arr, int isign) {
  // intent=inout allocatable general array
  fortran_fft_1d(/* void* */ arr.get_fortran_ptr(), /* c_Int& */ isign);
}
void SimUtils::file_directorizer(
    std::string in_file,
    std::string out_file,
    std::string directory,
    bool add_switch) {
  auto _in_file = in_file.c_str(); // ptr, inout, required
  auto _out_file = out_file.c_str(); // ptr, inout, required
  auto _directory = directory.c_str(); // ptr, inout, required
  fortran_file_directorizer(
      /* c_Char */ _in_file,
      /* c_Char */ _out_file,
      /* c_Char */ _directory,
      /* c_Bool& */ add_switch);
}
void SimUtils::file_get(
    std::string string,
    std::string dflt_file_name,
    std::string file_name) {
  auto _string = string.c_str(); // ptr, inout, required
  auto _dflt_file_name = dflt_file_name.c_str(); // ptr, inout, required
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_file_get(
      /* c_Char */ _string,
      /* c_Char */ _dflt_file_name,
      /* c_Char */ _file_name);
}
void SimUtils::file_get_open(
    std::string string,
    std::string dflt_file_name,
    std::string file_name,
    int file_unit,
    bool readonly) {
  auto _string = string.c_str(); // ptr, inout, required
  auto _dflt_file_name = dflt_file_name.c_str(); // ptr, inout, required
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_file_get_open(
      /* c_Char */ _string,
      /* c_Char */ _dflt_file_name,
      /* c_Char */ _file_name,
      /* c_Int& */ file_unit,
      /* c_Bool& */ readonly);
}
void SimUtils::file_suffixer(
    std::string in_file_name,
    std::string out_file_name,
    std::string suffix,
    bool add_switch) {
  auto _in_file_name = in_file_name.c_str(); // ptr, inout, required
  auto _out_file_name = out_file_name.c_str(); // ptr, inout, required
  auto _suffix = suffix.c_str(); // ptr, inout, required
  fortran_file_suffixer(
      /* c_Char */ _in_file_name,
      /* c_Char */ _out_file_name,
      /* c_Char */ _suffix,
      /* c_Bool& */ add_switch);
}
void SimUtils::gen_complete_elliptic(
    double kc,
    double p,
    double c,
    double s,
    optional_ref<double> err_tol,
    double value) {
  auto* _err_tol =
      err_tol.has_value() ? &err_tol->get() : nullptr; // inout, optional
  fortran_gen_complete_elliptic(
      /* c_Real& */ kc,
      /* c_Real& */ p,
      /* c_Real& */ c,
      /* c_Real& */ s,
      /* c_Real* */ _err_tol,
      /* c_Real& */ value);
}
void SimUtils::get_file_number(
    std::string file_name,
    std::string cnum_in,
    int num_out,
    bool err_flag) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  auto _cnum_in = cnum_in.c_str(); // ptr, inout, required
  fortran_get_file_number(
      /* c_Char */ _file_name,
      /* c_Char */ _cnum_in,
      /* c_Int& */ num_out,
      /* c_Bool& */ err_flag);
}
void SimUtils::get_file_time_stamp(std::string file, std::string time_stamp) {
  auto _file = file.c_str(); // ptr, inout, required
  auto _time_stamp = time_stamp.c_str(); // ptr, inout, required
  fortran_get_file_time_stamp(/* c_Char */ _file, /* c_Char */ _time_stamp);
}
void SimUtils::i_bessel(int m, double arg, double i_bes) {
  fortran_i_bessel(/* c_Int& */ m, /* c_Real& */ arg, /* c_Real& */ i_bes);
}
void SimUtils::i_bessel_extended(
    int m,
    double arg,
    std::complex<double> i_bes) {
  fortran_i_bessel_extended(
      /* c_Int& */ m, /* c_Real& */ arg, /* c_Complex& */ i_bes);
}
void SimUtils::increment_file_number(
    std::string file_name,
    int digits,
    int number,
    std::string cnumber) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  auto _cnumber = cnumber.c_str(); // ptr, inout, required
  fortran_increment_file_number(
      /* c_Char */ _file_name,
      /* c_Int& */ digits,
      /* c_Int& */ number,
      /* c_Char */ _cnumber);
}
void SimUtils::index_nocase(
    std::string string1,
    std::string string2,
    int indx) {
  auto _string1 = string1.c_str(); // ptr, inout, required
  auto _string2 = string2.c_str(); // ptr, inout, required
  fortran_index_nocase(
      /* c_Char */ _string1, /* c_Char */ _string2, /* c_Int& */ indx);
}
void SimUtils::int_str(int int_, optional_ref<int> width, std::string str) {
  auto* _width = width.has_value() ? &width->get() : nullptr; // inout, optional
  auto _str = str.c_str(); // ptr, inout, required
  fortran_int_str(/* c_Int& */ int_, /* c_Int* */ _width, /* c_Char */ _str);
}
void SimUtils::is_alphabetic(
    std::string string,
    optional_ref<std::string> valid_chars,
    bool is_alpha) {
  auto _string = string.c_str(); // ptr, inout, required
  const char* _valid_chars =
      valid_chars.has_value() ? valid_chars->get().c_str() : nullptr;
  fortran_is_alphabetic(
      /* c_Char */ _string, /* c_Char */ _valid_chars, /* c_Bool& */ is_alpha);
}
void SimUtils::is_decreasing_sequence(
    RealAllocatable1D& array,
    std::optional<bool> strict,
    bool is_decreasing) {
  // intent=in allocatable general array
  bool strict_lvalue;
  auto* _strict{&strict_lvalue};
  if (strict.has_value()) {
    strict_lvalue = strict.value();
  } else {
    _strict = nullptr;
  }
  fortran_is_decreasing_sequence(
      /* void* */ array.get_fortran_ptr(),
      /* c_Bool* */ _strict,
      /* c_Bool& */ is_decreasing);
}
void SimUtils::is_increasing_sequence(
    RealAllocatable1D& array,
    std::optional<bool> strict,
    bool is_increasing) {
  // intent=in allocatable general array
  bool strict_lvalue;
  auto* _strict{&strict_lvalue};
  if (strict.has_value()) {
    strict_lvalue = strict.value();
  } else {
    _strict = nullptr;
  }
  fortran_is_increasing_sequence(
      /* void* */ array.get_fortran_ptr(),
      /* c_Bool* */ _strict,
      /* c_Bool& */ is_increasing);
}
void SimUtils::is_integer(
    std::string string,
    optional_ref<int> int_,
    optional_ref<std::string> delims,
    optional_ref<int> ix_word,
    bool valid) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _int_ = int_.has_value() ? &int_->get() : nullptr; // inout, optional
  const char* _delims = delims.has_value() ? delims->get().c_str() : nullptr;
  auto* _ix_word =
      ix_word.has_value() ? &ix_word->get() : nullptr; // inout, optional
  fortran_is_integer(
      /* c_Char */ _string,
      /* c_Int* */ _int_,
      /* c_Char */ _delims,
      /* c_Int* */ _ix_word,
      /* c_Bool& */ valid);
}
void SimUtils::is_logical(
    std::string string,
    optional_ref<bool> ignore,
    bool valid) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _ignore =
      ignore.has_value() ? &ignore->get() : nullptr; // inout, optional
  fortran_is_logical(
      /* c_Char */ _string, /* c_Bool* */ _ignore, /* c_Bool& */ valid);
}
void SimUtils::is_real(
    std::string string,
    optional_ref<bool> ignore,
    optional_ref<double> real_num,
    bool valid) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _ignore =
      ignore.has_value() ? &ignore->get() : nullptr; // inout, optional
  auto* _real_num =
      real_num.has_value() ? &real_num->get() : nullptr; // inout, optional
  fortran_is_real(
      /* c_Char */ _string,
      /* c_Bool* */ _ignore,
      /* c_Real* */ _real_num,
      /* c_Bool& */ valid);
}
void SimUtils::j_bessel(int m, double arg, double j_bes) {
  fortran_j_bessel(/* c_Int& */ m, /* c_Real& */ arg, /* c_Real& */ j_bes);
}
void SimUtils::linear_fit(
    RealAllocatable1D& x,
    RealAllocatable1D& y,
    int n_data,
    double a,
    double b,
    double sig_a,
    double sig_b) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  fortran_linear_fit(
      /* void* */ x.get_fortran_ptr(),
      /* void* */ y.get_fortran_ptr(),
      /* c_Int& */ n_data,
      /* c_Real& */ a,
      /* c_Real& */ b,
      /* c_Real& */ sig_a,
      /* c_Real& */ sig_b);
}
FixedArray1D<Real, 3> SimUtils::linear_fit_2d(
    RealAllocatable1D& x,
    RealAllocatable1D& y,
    RealAllocatable1D& z) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=in allocatable general array
  FixedArray1D<Real, 3> _coef;
  fortran_linear_fit_2d(
      /* void* */ x.get_fortran_ptr(),
      /* void* */ y.get_fortran_ptr(),
      /* void* */ z.get_fortran_ptr(),
      /* c_RealArr */ _coef.data());
  return _coef;
}
void SimUtils::logic_str(bool logic, std::string str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_logic_str(/* c_Bool& */ logic, /* c_Char */ _str);
}
int SimUtils::lunget() {
  int _func_retval__{};
  fortran_lunget(/* c_Int& */ _func_retval__);
  return _func_retval__;
}
void SimUtils::match_reg(std::string str, std::string pat, bool is_match) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _pat = pat.c_str(); // ptr, inout, required
  fortran_match_reg(
      /* c_Char */ _str, /* c_Char */ _pat, /* c_Bool& */ is_match);
}
void SimUtils::milli_sleep(int milli_sec) {
  fortran_milli_sleep(/* c_Int& */ milli_sec);
}
void SimUtils::make_legal_comment(
    std::string comment_in,
    std::string comment_out) {
  auto _comment_in = comment_in.c_str(); // ptr, inout, required
  auto _comment_out = comment_out.c_str(); // ptr, inout, required
  fortran_make_legal_comment(
      /* c_Char */ _comment_in, /* c_Char */ _comment_out);
}
void SimUtils::match_wild(
    std::string string,
    std::string template_,
    bool is_match) {
  auto _string = string.c_str(); // ptr, inout, required
  auto _template_ = template_.c_str(); // ptr, inout, required
  fortran_match_wild(
      /* c_Char */ _string, /* c_Char */ _template_, /* c_Bool& */ is_match);
}
void SimUtils::n_choose_k(int n, int k, double nck) {
  fortran_n_choose_k(/* c_Int& */ n, /* c_Int& */ k, /* c_Real& */ nck);
}
RealAllocatable1D SimUtils::n_spline_create(
    RealAllocatable1D& deriv0,
    RealAllocatable1D& deriv1,
    double x1) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=out allocatable general array
  auto n_spline{RealAllocatable1D()};
  fortran_n_spline_create(
      /* void* */ deriv0.get_fortran_ptr(),
      /* void* */ deriv1.get_fortran_ptr(),
      /* c_Real& */ x1,
      /* void* */ n_spline.get_fortran_ptr());
  return std::move(n_spline);
}
void SimUtils::ordinal_str(int n, std::string str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_ordinal_str(/* c_Int& */ n, /* c_Char */ _str);
}
void SimUtils::parse_fortran_format(
    std::string format_str,
    int n_repeat,
    int power,
    std::string descrip,
    int width,
    int digits) {
  auto _format_str = format_str.c_str(); // ptr, inout, required
  auto _descrip = descrip.c_str(); // ptr, inout, required
  fortran_parse_fortran_format(
      /* c_Char */ _format_str,
      /* c_Int& */ n_repeat,
      /* c_Int& */ power,
      /* c_Char */ _descrip,
      /* c_Int& */ width,
      /* c_Int& */ digits);
}
void SimUtils::poly_eval(
    RealAllocatable1D& poly,
    double x,
    std::optional<bool> diff_coef,
    double y) {
  // intent=in allocatable general array
  bool diff_coef_lvalue;
  auto* _diff_coef{&diff_coef_lvalue};
  if (diff_coef.has_value()) {
    diff_coef_lvalue = diff_coef.value();
  } else {
    _diff_coef = nullptr;
  }
  fortran_poly_eval(
      /* void* */ poly.get_fortran_ptr(),
      /* c_Real& */ x,
      /* c_Bool* */ _diff_coef,
      /* c_Real& */ y);
}
void SimUtils::probability_funct(double x, double prob) {
  fortran_probability_funct(/* c_Real& */ x, /* c_Real& */ prob);
}
void SimUtils::quadratic_roots(
    FixedArray1D<Real, 3> coefs,
    FixedArray1D<Complex, 2> root) {
  auto* _coefs = coefs.data(); // CppWrapperGeneralArgument
  auto* _root = root.data(); // CppWrapperGeneralArgument
  fortran_quadratic_roots(/* c_RealArr */ _coefs, /* c_ComplexArr */ _root);
}
void SimUtils::query_string(
    std::string query_str,
    bool upcase,
    std::string return_str,
    int ix,
    int ios) {
  auto _query_str = query_str.c_str(); // ptr, inout, required
  auto _return_str = return_str.c_str(); // ptr, inout, required
  fortran_query_string(
      /* c_Char */ _query_str,
      /* c_Bool& */ upcase,
      /* c_Char */ _return_str,
      /* c_Int& */ ix,
      /* c_Int& */ ios);
}
void SimUtils::quote(std::string str, std::string q_str) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _q_str = q_str.c_str(); // ptr, inout, required
  fortran_quote(/* c_Char */ _str, /* c_Char */ _q_str);
}
void SimUtils::real_to_string(
    double real_num,
    int width,
    optional_ref<int> n_signif,
    optional_ref<int> n_decimal,
    std::string str) {
  auto* _n_signif =
      n_signif.has_value() ? &n_signif->get() : nullptr; // inout, optional
  auto* _n_decimal =
      n_decimal.has_value() ? &n_decimal->get() : nullptr; // inout, optional
  auto _str = str.c_str(); // ptr, inout, required
  fortran_real_to_string(
      /* c_Real& */ real_num,
      /* c_Int& */ width,
      /* c_Int* */ _n_signif,
      /* c_Int* */ _n_decimal,
      /* c_Char */ _str);
}
void SimUtils::real_num_fortran_format(
    double number,
    int width,
    optional_ref<int> n_blanks,
    std::string fmt_str) {
  auto* _n_blanks =
      n_blanks.has_value() ? &n_blanks->get() : nullptr; // inout, optional
  auto _fmt_str = fmt_str.c_str(); // ptr, inout, required
  fortran_real_num_fortran_format(
      /* c_Real& */ number,
      /* c_Int& */ width,
      /* c_Int* */ _n_blanks,
      /* c_Char */ _fmt_str);
}
void SimUtils::str_count(std::string str, std::string match, int num) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _match = match.c_str(); // ptr, inout, required
  fortran_str_count(/* c_Char */ _str, /* c_Char */ _match, /* c_Int& */ num);
}
void SimUtils::real_path(
    std::string path_in,
    std::string path_out,
    bool is_ok) {
  auto _path_in = path_in.c_str(); // ptr, inout, required
  auto _path_out = path_out.c_str(); // ptr, inout, required
  fortran_real_path(
      /* c_Char */ _path_in, /* c_Char */ _path_out, /* c_Bool& */ is_ok);
}
void SimUtils::real_str(
    double r_num,
    optional_ref<int> n_signif,
    optional_ref<int> n_decimal,
    std::string str) {
  auto* _n_signif =
      n_signif.has_value() ? &n_signif->get() : nullptr; // inout, optional
  auto* _n_decimal =
      n_decimal.has_value() ? &n_decimal->get() : nullptr; // inout, optional
  auto _str = str.c_str(); // ptr, inout, required
  fortran_real_str(
      /* c_Real& */ r_num,
      /* c_Int* */ _n_signif,
      /* c_Int* */ _n_decimal,
      /* c_Char */ _str);
}
double SimUtils::rms_value(
    RealAllocatable1D& val_arr,
    optional_ref<BoolAllocatable1D> good_val,
    double rms_val) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  auto* _good_val = good_val.has_value() ? good_val->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  double _ave_val{};
  fortran_rms_value(
      /* void* */ val_arr.get_fortran_ptr(),
      /* void* */ _good_val,
      /* c_Real& */ _ave_val,
      /* c_Real& */ rms_val);
  return _ave_val;
}
void SimUtils::rot_2d(
    FixedArray1D<Real, 2> vec_in,
    double angle,
    FixedArray1D<Real, 2> vec_out) {
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  auto* _vec_out = vec_out.data(); // CppWrapperGeneralArgument
  fortran_rot_2d(
      /* c_RealArr */ _vec_in, /* c_Real& */ angle, /* c_RealArr */ _vec_out);
}
void SimUtils::run_timer(
    std::string command,
    optional_ref<double> time,
    optional_ref<double> time0) {
  auto _command = command.c_str(); // ptr, inout, required
  auto* _time = time.has_value() ? &time->get() : nullptr; // inout, optional
  auto* _time0 = time0.has_value() ? &time0->get() : nullptr; // inout, optional
  fortran_run_timer(
      /* c_Char */ _command, /* c_Real* */ _time, /* c_Real* */ _time0);
}
void SimUtils::sinc(double x, optional_ref<int> nd, double y) {
  auto* _nd = nd.has_value() ? &nd->get() : nullptr; // inout, optional
  fortran_sinc(/* c_Real& */ x, /* c_Int* */ _nd, /* c_Real& */ y);
}
void SimUtils::sincc(double x, optional_ref<int> nd, double y) {
  auto* _nd = nd.has_value() ? &nd->get() : nullptr; // inout, optional
  fortran_sincc(/* c_Real& */ x, /* c_Int* */ _nd, /* c_Real& */ y);
}
void SimUtils::sinhx_x(double x, optional_ref<int> nd, double y) {
  auto* _nd = nd.has_value() ? &nd->get() : nullptr; // inout, optional
  fortran_sinhx_x(/* c_Real& */ x, /* c_Int* */ _nd, /* c_Real& */ y);
}
void SimUtils::skip_header(int ix_unit, bool error_flag) {
  fortran_skip_header(/* c_Int& */ ix_unit, /* c_Bool& */ error_flag);
}
void SimUtils::sqrt_one(double x, optional_ref<int> nd, double ds1) {
  auto* _nd = nd.has_value() ? &nd->get() : nullptr; // inout, optional
  fortran_sqrt_one(/* c_Real& */ x, /* c_Int* */ _nd, /* c_Real& */ ds1);
}
void SimUtils::sqrt_alpha(double alpha, double x, double y) {
  fortran_sqrt_alpha(/* c_Real& */ alpha, /* c_Real& */ x, /* c_Real& */ y);
}
void SimUtils::str_first_in_set(
    std::string line,
    std::string set,
    optional_ref<bool> ignore_clauses,
    int ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  auto* _ignore_clauses = ignore_clauses.has_value()
      ? &ignore_clauses->get()
      : nullptr; // inout, optional
  fortran_str_first_in_set(
      /* c_Char */ _line,
      /* c_Char */ _set,
      /* c_Bool* */ _ignore_clauses,
      /* c_Int& */ ix_match);
}
void SimUtils::str_first_not_in_set(
    std::string line,
    std::string set,
    int ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  fortran_str_first_not_in_set(
      /* c_Char */ _line, /* c_Char */ _set, /* c_Int& */ ix_match);
}
void SimUtils::str_last_in_set(
    std::string line,
    std::string set,
    int ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  fortran_str_last_in_set(
      /* c_Char */ _line, /* c_Char */ _set, /* c_Int& */ ix_match);
}
void SimUtils::str_last_not_in_set(
    std::string line,
    std::string set,
    int ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  fortran_str_last_not_in_set(
      /* c_Char */ _line, /* c_Char */ _set, /* c_Int& */ ix_match);
}
void SimUtils::string_to_int(
    std::string line,
    int default_,
    bool err_flag,
    optional_ref<bool> err_print_flag,
    int value) {
  auto _line = line.c_str(); // ptr, inout, required
  auto* _err_print_flag = err_print_flag.has_value()
      ? &err_print_flag->get()
      : nullptr; // inout, optional
  fortran_string_to_int(
      /* c_Char */ _line,
      /* c_Int& */ default_,
      /* c_Bool& */ err_flag,
      /* c_Bool* */ _err_print_flag,
      /* c_Int& */ value);
}
void SimUtils::string_to_real(
    std::string line,
    double default_,
    bool err_flag,
    optional_ref<bool> err_print_flag,
    double value) {
  auto _line = line.c_str(); // ptr, inout, required
  auto* _err_print_flag = err_print_flag.has_value()
      ? &err_print_flag->get()
      : nullptr; // inout, optional
  fortran_string_to_real(
      /* c_Char */ _line,
      /* c_Real& */ default_,
      /* c_Bool& */ err_flag,
      /* c_Bool* */ _err_print_flag,
      /* c_Real& */ value);
}
void SimUtils::string_trim2(
    std::string in_str,
    std::string delimitors,
    std::string out_str,
    int ix_word,
    std::string delim,
    int ix_next) {
  auto _in_str = in_str.c_str(); // ptr, inout, required
  auto _delimitors = delimitors.c_str(); // ptr, inout, required
  auto _out_str = out_str.c_str(); // ptr, inout, required
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_string_trim2(
      /* c_Char */ _in_str,
      /* c_Char */ _delimitors,
      /* c_Char */ _out_str,
      /* c_Int& */ ix_word,
      /* c_Char */ _delim,
      /* c_Int& */ ix_next);
}
void SimUtils::to_str(
    double num,
    optional_ref<int> max_signif,
    std::string string) {
  auto* _max_signif =
      max_signif.has_value() ? &max_signif->get() : nullptr; // inout, optional
  auto _string = string.c_str(); // ptr, inout, required
  fortran_to_str(
      /* c_Real& */ num, /* c_Int* */ _max_signif, /* c_Char */ _string);
}
void SimUtils::type_this_file(std::string filename) {
  auto _filename = filename.c_str(); // ptr, inout, required
  fortran_type_this_file(/* c_Char */ _filename);
}
void SimUtils::upcase_string(std::string string) {
  auto _string = string.c_str(); // ptr, inout, required
  fortran_upcase_string(/* c_Char */ _string);
}
void SimUtils::word_len(std::string wording, int wlen) {
  auto _wording = wording.c_str(); // ptr, inout, required
  fortran_word_len(/* c_Char */ _wording, /* c_Int& */ wlen);
}
void SimUtils::word_read(
    std::string in_str,
    std::string delim_list,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    std::string out_str,
    optional_ref<bool> ignore_interior) {
  auto _in_str = in_str.c_str(); // ptr, inout, required
  auto _delim_list = delim_list.c_str(); // ptr, inout, required
  auto _word = word.c_str(); // ptr, inout, required
  auto _delim = delim.c_str(); // ptr, inout, required
  auto _out_str = out_str.c_str(); // ptr, inout, required
  auto* _ignore_interior = ignore_interior.has_value()
      ? &ignore_interior->get()
      : nullptr; // inout, optional
  fortran_word_read(
      /* c_Char */ _in_str,
      /* c_Char */ _delim_list,
      /* c_Char */ _word,
      /* c_Int& */ ix_word,
      /* c_Char */ _delim,
      /* c_Bool& */ delim_found,
      /* c_Char */ _out_str,
      /* c_Bool* */ _ignore_interior);
}
void SimUtils::str_substitute(
    std::string string,
    optional_ref<std::string> str_match,
    optional_ref<std::string> str_replace,
    optional_ref<bool> do_trim,
    optional_ref<bool> ignore_escaped) {
  auto _string = string.c_str(); // ptr, inout, required
  const char* _str_match =
      str_match.has_value() ? str_match->get().c_str() : nullptr;
  const char* _str_replace =
      str_replace.has_value() ? str_replace->get().c_str() : nullptr;
  auto* _do_trim =
      do_trim.has_value() ? &do_trim->get() : nullptr; // inout, optional
  auto* _ignore_escaped = ignore_escaped.has_value()
      ? &ignore_escaped->get()
      : nullptr; // inout, optional
  fortran_str_substitute(
      /* c_Char */ _string,
      /* c_Char */ _str_match,
      /* c_Char */ _str_replace,
      /* c_Bool* */ _do_trim,
      /* c_Bool* */ _ignore_escaped);
}
void SimUtils::str_match_wild(std::string str, std::string pat, bool a_match) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _pat = pat.c_str(); // ptr, inout, required
  fortran_str_match_wild(
      /* c_Char */ _str, /* c_Char */ _pat, /* c_Bool& */ a_match);
}
void SimUtils::str_upcase(std::string dst, std::string src) {
  auto _dst = dst.c_str(); // ptr, inout, required
  auto _src = src.c_str(); // ptr, inout, required
  fortran_str_upcase(/* c_Char */ _dst, /* c_Char */ _src);
}
void SimUtils::str_downcase(std::string dst, std::string src) {
  auto _dst = dst.c_str(); // ptr, inout, required
  auto _src = src.c_str(); // ptr, inout, required
  fortran_str_downcase(/* c_Char */ _dst, /* c_Char */ _src);
}
void SimUtils::system_command(std::string line, optional_ref<bool> err_flag) {
  auto _line = line.c_str(); // ptr, inout, required
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  fortran_system_command(/* c_Char */ _line, /* c_Bool* */ _err_flag);
}
void SimUtils::string_trim(
    std::string in_string,
    std::string out_string,
    int word_len) {
  auto _in_string = in_string.c_str(); // ptr, inout, required
  auto _out_string = out_string.c_str(); // ptr, inout, required
  fortran_string_trim(
      /* c_Char */ _in_string, /* c_Char */ _out_string, /* c_Int& */ word_len);
}
int SimUtils::virtual_memory_usage() {
  int _usage{};
  fortran_virtual_memory_usage(/* c_Int& */ _usage);
  return _usage;
}
void SimUtils::find_location_real(
    RealAllocatable1D& arr,
    double value,
    int ix_match) {
  // intent=in allocatable general array
  fortran_find_location_real(
      /* void* */ arr.get_fortran_ptr(),
      /* c_Real& */ value,
      /* c_Int& */ ix_match);
}
void SimUtils::find_location_int(
    IntAllocatable1D& arr,
    int value,
    int ix_match) {
  // intent=inout allocatable general array
  fortran_find_location_int(
      /* void* */ arr.get_fortran_ptr(),
      /* c_Int& */ value,
      /* c_Int& */ ix_match);
}
void SimUtils::find_location_logic(
    BoolAllocatable1D& arr,
    bool value,
    int ix_match) {
  // intent=inout allocatable general array
  fortran_find_location_logic(
      /* void* */ arr.get_fortran_ptr(),
      /* c_Bool& */ value,
      /* c_Int& */ ix_match);
}
double SimUtils::coarse_frequency_estimate(
    RealAllocatable1D& data,
    optional_ref<bool> error) {
  // intent=in allocatable general array
  auto* _error = error.has_value() ? &error->get() : nullptr; // inout, optional
  double _frequency{};
  fortran_coarse_frequency_estimate(
      /* void* */ data.get_fortran_ptr(),
      /* c_Bool* */ _error,
      /* c_Real& */ _frequency);
  return _frequency;
}
double SimUtils::fine_frequency_estimate(RealAllocatable1D& data) {
  // intent=in allocatable general array
  double _frequency{};
  fortran_fine_frequency_estimate(
      /* void* */ data.get_fortran_ptr(), /* c_Real& */ _frequency);
  return _frequency;
}
SimUtils::FourierAmplitude SimUtils::fourier_amplitude(
    RealAllocatable1D& data,
    double frequency) {
  // intent=in allocatable general array
  double _cos_amp{};
  double _sin_amp{};
  double _dcos_amp{};
  double _dsin_amp{};
  fortran_fourier_amplitude(
      /* void* */ data.get_fortran_ptr(),
      /* c_Real& */ frequency,
      /* c_Real& */ _cos_amp,
      /* c_Real& */ _sin_amp,
      /* c_Real& */ _dcos_amp,
      /* c_Real& */ _dsin_amp);
  return FourierAmplitude{_cos_amp, _sin_amp, _dcos_amp, _dsin_amp};
}
SimUtils::WMatToAxisAngle SimUtils::w_mat_to_axis_angle(
    FixedArray2D<Real, 3, 3> w_mat) {
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  FixedArray1D<Real, 3> _axis;
  double _angle{};
  fortran_w_mat_to_axis_angle(
      /* c_RealArr */ _w_mat_vec,
      /* c_RealArr */ _axis.data(),
      /* c_Real& */ _angle);
  return WMatToAxisAngle{_axis, _angle};
}
FixedArray1D<Real, 4> SimUtils::w_mat_to_quat(FixedArray2D<Real, 3, 3> w_mat) {
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  FixedArray1D<Real, 4> _quat;
  fortran_w_mat_to_quat(
      /* c_RealArr */ _w_mat_vec, /* c_RealArr */ _quat.data());
  return _quat;
}
FixedArray2D<Real, 3, 3> SimUtils::quat_to_w_mat(FixedArray1D<Real, 4> quat) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  fortran_quat_to_w_mat(/* c_RealArr */ _quat, /* c_RealArr */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
FixedArray2D<Real, 3, 3> SimUtils::axis_angle_to_w_mat(
    FixedArray1D<Real, 3> axis,
    double angle) {
  auto* _axis = axis.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  fortran_axis_angle_to_w_mat(
      /* c_RealArr */ _axis, /* c_Real& */ angle, /* c_RealArr */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
FixedArray1D<Real, 3> SimUtils::quat_to_omega(FixedArray1D<Real, 4> quat) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _omega;
  fortran_quat_to_omega(/* c_RealArr */ _quat, /* c_RealArr */ _omega.data());
  return _omega;
}
FixedArray1D<Real, 4> SimUtils::omega_to_quat(FixedArray1D<Real, 3> omega) {
  auto* _omega = omega.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _quat;
  fortran_omega_to_quat(/* c_RealArr */ _omega, /* c_RealArr */ _quat.data());
  return _quat;
}
SimUtils::QuatToAxisAngle SimUtils::quat_to_axis_angle(
    FixedArray1D<Real, 4> quat) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _axis;
  double _angle{};
  fortran_quat_to_axis_angle(
      /* c_RealArr */ _quat,
      /* c_RealArr */ _axis.data(),
      /* c_Real& */ _angle);
  return QuatToAxisAngle{_axis, _angle};
}
FixedArray1D<Real, 4> SimUtils::axis_angle_to_quat(
    FixedArray1D<Real, 3> axis,
    double angle) {
  auto* _axis = axis.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _quat;
  fortran_axis_angle_to_quat(
      /* c_RealArr */ _axis, /* c_Real& */ angle, /* c_RealArr */ _quat.data());
  return _quat;
}
FixedArray1D<Real, 4> SimUtils::quat_conj_real(FixedArray1D<Real, 4> q_in) {
  auto* _q_in = q_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _q_out;
  fortran_quat_conj_real(/* c_RealArr */ _q_in, /* c_RealArr */ _q_out.data());
  return _q_out;
}
FixedArray1D<Complex, 4> SimUtils::quat_conj_complex(
    FixedArray1D<Complex, 4> q_in) {
  auto* _q_in = q_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Complex, 4> _q_out;
  fortran_quat_conj_complex(
      /* c_ComplexArr */ _q_in, /* c_ComplexArr */ _q_out.data());
  return _q_out;
}
FixedArray1D<Real, 4> SimUtils::quat_inverse(FixedArray1D<Real, 4> q_in) {
  auto* _q_in = q_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _q_out;
  fortran_quat_inverse(/* c_RealArr */ _q_in, /* c_RealArr */ _q_out.data());
  return _q_out;
}
FixedArray1D<Real, 4> SimUtils::quat_mul_real(
    FixedArray1D<Real, 4> q1,
    FixedArray1D<Real, 4> q2,
    std::optional<FixedArray1D<Real, 4>> q3,
    std::optional<FixedArray1D<Real, 4>> q4,
    std::optional<FixedArray1D<Real, 4>> q5,
    std::optional<FixedArray1D<Real, 4>> q6,
    std::optional<FixedArray1D<Real, 4>> q7,
    std::optional<FixedArray1D<Real, 4>> q8,
    std::optional<FixedArray1D<Real, 4>> q9) {
  auto* _q1 = q1.data(); // CppWrapperGeneralArgument
  auto* _q2 = q2.data(); // CppWrapperGeneralArgument
  c_RealArr _q3 = q3.has_value() ? q3.value().data() : nullptr;
  c_RealArr _q4 = q4.has_value() ? q4.value().data() : nullptr;
  c_RealArr _q5 = q5.has_value() ? q5.value().data() : nullptr;
  c_RealArr _q6 = q6.has_value() ? q6.value().data() : nullptr;
  c_RealArr _q7 = q7.has_value() ? q7.value().data() : nullptr;
  c_RealArr _q8 = q8.has_value() ? q8.value().data() : nullptr;
  c_RealArr _q9 = q9.has_value() ? q9.value().data() : nullptr;
  FixedArray1D<Real, 4> _q_out;
  fortran_quat_mul_real(
      /* c_RealArr */ _q1,
      /* c_RealArr */ _q2,
      /* c_RealArr */ _q3,
      /* c_RealArr */ _q4,
      /* c_RealArr */ _q5,
      /* c_RealArr */ _q6,
      /* c_RealArr */ _q7,
      /* c_RealArr */ _q8,
      /* c_RealArr */ _q9,
      /* c_RealArr */ _q_out.data());
  return _q_out;
}
FixedArray1D<Complex, 4> SimUtils::quat_mul_complex(
    FixedArray1D<Complex, 4> q1,
    FixedArray1D<Complex, 4> q2,
    std::optional<FixedArray1D<Complex, 4>> q3,
    std::optional<FixedArray1D<Complex, 4>> q4,
    std::optional<FixedArray1D<Complex, 4>> q5,
    std::optional<FixedArray1D<Complex, 4>> q6,
    std::optional<FixedArray1D<Complex, 4>> q7,
    std::optional<FixedArray1D<Complex, 4>> q8,
    std::optional<FixedArray1D<Complex, 4>> q9) {
  auto* _q1 = q1.data(); // CppWrapperGeneralArgument
  auto* _q2 = q2.data(); // CppWrapperGeneralArgument
  c_ComplexArr _q3 = q3.has_value() ? q3.value().data() : nullptr;
  c_ComplexArr _q4 = q4.has_value() ? q4.value().data() : nullptr;
  c_ComplexArr _q5 = q5.has_value() ? q5.value().data() : nullptr;
  c_ComplexArr _q6 = q6.has_value() ? q6.value().data() : nullptr;
  c_ComplexArr _q7 = q7.has_value() ? q7.value().data() : nullptr;
  c_ComplexArr _q8 = q8.has_value() ? q8.value().data() : nullptr;
  c_ComplexArr _q9 = q9.has_value() ? q9.value().data() : nullptr;
  FixedArray1D<Complex, 4> _q_out;
  fortran_quat_mul_complex(
      /* c_ComplexArr */ _q1,
      /* c_ComplexArr */ _q2,
      /* c_ComplexArr */ _q3,
      /* c_ComplexArr */ _q4,
      /* c_ComplexArr */ _q5,
      /* c_ComplexArr */ _q6,
      /* c_ComplexArr */ _q7,
      /* c_ComplexArr */ _q8,
      /* c_ComplexArr */ _q9,
      /* c_ComplexArr */ _q_out.data());
  return _q_out;
}
FixedArray1D<Real, 3> SimUtils::quat_rotate_real(
    FixedArray1D<Real, 4> quat,
    FixedArray1D<Real, 3> vec_in) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _vec_out;
  fortran_quat_rotate_real(
      /* c_RealArr */ _quat,
      /* c_RealArr */ _vec_in,
      /* c_RealArr */ _vec_out.data());
  return _vec_out;
}
FixedArray1D<Complex, 3> SimUtils::quat_rotate_complex(
    FixedArray1D<Complex, 4> quat,
    FixedArray1D<Complex, 3> vec_in) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Complex, 3> _vec_out;
  fortran_quat_rotate_complex(
      /* c_ComplexArr */ _quat,
      /* c_ComplexArr */ _vec_in,
      /* c_ComplexArr */ _vec_out.data());
  return _vec_out;
}
FixedArray1D<Real, 3> SimUtils::rotate_vec_given_axis_angle(
    FixedArray1D<Real, 3> vec_in,
    RealAllocatable1D& axis,
    double angle) {
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  // intent=in allocatable general array
  FixedArray1D<Real, 3> _vec_out;
  fortran_rotate_vec_given_axis_angle(
      /* c_RealArr */ _vec_in,
      /* void* */ axis.get_fortran_ptr(),
      /* c_Real& */ angle,
      /* c_RealArr */ _vec_out.data());
  return _vec_out;
}
void SimUtils::rotate_vec(RealAllocatable1D& vec, int axis, double angle) {
  // intent=inout allocatable general array
  fortran_rotate_vec(
      /* void* */ vec.get_fortran_ptr(),
      /* c_Int& */ axis,
      /* c_Real& */ angle);
}
void SimUtils::naff(
    ComplexAllocatable1D& cdata,
    RealAllocatable1D& freqs,
    ComplexAllocatable1D& amps,
    optional_ref<int> opt_dump_spectra,
    optional_ref<bool> opt_zero_first) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  auto* _opt_dump_spectra = opt_dump_spectra.has_value()
      ? &opt_dump_spectra->get()
      : nullptr; // inout, optional
  auto* _opt_zero_first = opt_zero_first.has_value()
      ? &opt_zero_first->get()
      : nullptr; // inout, optional
  fortran_naff(
      /* void* */ cdata.get_fortran_ptr(),
      /* void* */ freqs.get_fortran_ptr(),
      /* void* */ amps.get_fortran_ptr(),
      /* c_Int* */ _opt_dump_spectra,
      /* c_Bool* */ _opt_zero_first);
}
void SimUtils::projdd(
    ComplexAllocatable1D& a,
    ComplexAllocatable1D& b,
    std::complex<double> func_retval__) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  fortran_projdd(
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* c_Complex& */ func_retval__);
}
void SimUtils::maximize_projection(
    double seed,
    ComplexAllocatable1D& cdata,
    double func_retval__) {
  // intent=inout allocatable general array
  fortran_maximize_projection(
      /* c_Real& */ seed,
      /* void* */ cdata.get_fortran_ptr(),
      /* c_Real& */ func_retval__);
}
void SimUtils::interpolated_fft_gsl(
    ComplexAllocatable1D& cdata,
    bool calc_ok,
    optional_ref<int> opt_dump_spectrum,
    optional_ref<int> opt_dump_index,
    double this_fft) {
  // intent=inout allocatable general array
  auto* _opt_dump_spectrum = opt_dump_spectrum.has_value()
      ? &opt_dump_spectrum->get()
      : nullptr; // inout, optional
  auto* _opt_dump_index = opt_dump_index.has_value()
      ? &opt_dump_index->get()
      : nullptr; // inout, optional
  fortran_interpolated_fft_gsl(
      /* void* */ cdata.get_fortran_ptr(),
      /* c_Bool& */ calc_ok,
      /* c_Int* */ _opt_dump_spectrum,
      /* c_Int* */ _opt_dump_index,
      /* c_Real& */ this_fft);
}
void SimUtils::interpolated_fft(
    ComplexAllocatable1D& cdata,
    bool calc_ok,
    optional_ref<int> opt_dump_spectrum,
    optional_ref<int> opt_dump_index,
    double this_fft) {
  // intent=inout allocatable general array
  auto* _opt_dump_spectrum = opt_dump_spectrum.has_value()
      ? &opt_dump_spectrum->get()
      : nullptr; // inout, optional
  auto* _opt_dump_index = opt_dump_index.has_value()
      ? &opt_dump_index->get()
      : nullptr; // inout, optional
  fortran_interpolated_fft(
      /* void* */ cdata.get_fortran_ptr(),
      /* c_Bool& */ calc_ok,
      /* c_Int* */ _opt_dump_spectrum,
      /* c_Int* */ _opt_dump_index,
      /* c_Real& */ this_fft);
}
int SimUtils::initfixedwindowls(int N, double dt, int order, int der) {
  int _id{};
  fortran_initfixedwindowls(
      /* c_Int& */ N,
      /* c_Real& */ dt,
      /* c_Int& */ order,
      /* c_Int& */ der,
      /* c_Int& */ _id);
  return _id;
}
void SimUtils::destfixedwindowls(int id) {
  fortran_destfixedwindowls(/* c_Int& */ id);
}
void SimUtils::fixedwindowls(double ynew, int id, double z) {
  fortran_fixedwindowls(/* c_Real& */ ynew, /* c_Int& */ id, /* c_Real& */ z);
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
    double x2) {
  auto* _y = y.data(); // CppWrapperGeneralArgument
  auto* _y1 = y1.data(); // CppWrapperGeneralArgument
  auto* _y2 = y2.data(); // CppWrapperGeneralArgument
  auto* _y12 = y12.data(); // CppWrapperGeneralArgument
  double _ansy{};
  double _ansy1{};
  double _ansy2{};
  fortran_super_bicubic_interpolation(
      /* c_RealArr */ _y,
      /* c_RealArr */ _y1,
      /* c_RealArr */ _y2,
      /* c_RealArr */ _y12,
      /* c_Real& */ x1l,
      /* c_Real& */ x1u,
      /* c_Real& */ x2l,
      /* c_Real& */ x2u,
      /* c_Real& */ x1,
      /* c_Real& */ x2,
      /* c_Real& */ _ansy,
      /* c_Real& */ _ansy1,
      /* c_Real& */ _ansy2);
  return SuperBicubicInterpolation{_ansy, _ansy1, _ansy2};
}
FixedArray2D<Real, 4, 4> SimUtils::super_bicubic_coef(
    FixedArray1D<Real, 4> y,
    FixedArray1D<Real, 4> y1,
    FixedArray1D<Real, 4> y2,
    FixedArray1D<Real, 4> y12,
    double d1,
    double d2) {
  auto* _y = y.data(); // CppWrapperGeneralArgument
  auto* _y1 = y1.data(); // CppWrapperGeneralArgument
  auto* _y2 = y2.data(); // CppWrapperGeneralArgument
  auto* _y12 = y12.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 4, 4> c;
  double _c_vec[4 * 4];
  fortran_super_bicubic_coef(
      /* c_RealArr */ _y,
      /* c_RealArr */ _y1,
      /* c_RealArr */ _y2,
      /* c_RealArr */ _y12,
      /* c_Real& */ d1,
      /* c_Real& */ d2,
      /* c_RealArr */ _c_vec);
  vec_to_matrix(_c_vec, c);
  return c;
}
void SimUtils::super_sort(IntAllocatable1D& arr) {
  // intent=inout allocatable general array
  fortran_super_sort(/* void* */ arr.get_fortran_ptr());
}
SimUtils::SuperPolint SimUtils::super_polint(
    RealAllocatable1D& xa,
    RealAllocatable1D& ya,
    double x) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  double _y{};
  double _dy{};
  fortran_super_polint(
      /* void* */ xa.get_fortran_ptr(),
      /* void* */ ya.get_fortran_ptr(),
      /* c_Real& */ x,
      /* c_Real& */ _y,
      /* c_Real& */ _dy);
  return SuperPolint{_y, _dy};
}
double SimUtils::super_poly(double x, RealAllocatable1D& coeffs) {
  // intent=inout allocatable general array
  double _value{};
  fortran_super_poly(
      /* c_Real& */ x,
      /* void* */ coeffs.get_fortran_ptr(),
      /* c_Real& */ _value);
  return _value;
}
void SimUtils::ran_seed_put(int seed, std::optional<int> mpi_offset) {
  int mpi_offset_lvalue;
  auto* _mpi_offset{&mpi_offset_lvalue};
  if (mpi_offset.has_value()) {
    mpi_offset_lvalue = mpi_offset.value();
  } else {
    _mpi_offset = nullptr;
  }
  fortran_ran_seed_put(/* c_Int& */ seed, /* c_Int* */ _mpi_offset);
}
int SimUtils::ran_seed_get() {
  int _seed{};
  fortran_ran_seed_get(/* c_Int& */ _seed);
  return _seed;
}
void SimUtils::allocate_thread_states() {
  fortran_allocate_thread_states();
}
SimUtils::BicubicCmplxEval SimUtils::bicubic_cmplx_eval(
    double x_norm,
    double y_norm,
    BicubicCmplxCoefProxy& bi_coef) {
  std::complex<double> _df_dx{};
  std::complex<double> _df_dy{};
  std::complex<double> _f_val{};
  fortran_bicubic_cmplx_eval(
      /* c_Real& */ x_norm,
      /* c_Real& */ y_norm,
      /* void* */ bi_coef.get_fortran_ptr(),
      /* c_Complex& */ _df_dx,
      /* c_Complex& */ _df_dy,
      /* c_Complex& */ _f_val);
  return BicubicCmplxEval{_df_dx, _df_dy, _f_val};
}
SimUtils::TricubicCmplxEval SimUtils::tricubic_cmplx_eval(
    double x_norm,
    double y_norm,
    double z_norm,
    TricubicCmplxCoefProxy& tri_coef) {
  std::complex<double> _df_dx{};
  std::complex<double> _df_dy{};
  std::complex<double> _df_dz{};
  std::complex<double> _f_val{};
  fortran_tricubic_cmplx_eval(
      /* c_Real& */ x_norm,
      /* c_Real& */ y_norm,
      /* c_Real& */ z_norm,
      /* void* */ tri_coef.get_fortran_ptr(),
      /* c_Complex& */ _df_dx,
      /* c_Complex& */ _df_dy,
      /* c_Complex& */ _df_dz,
      /* c_Complex& */ _f_val);
  return TricubicCmplxEval{_df_dx, _df_dy, _df_dz, _f_val};
}
int SimUtils::bin_index(double x, double bin1_x_min, double bin_delta) {
  int _ix_bin{};
  fortran_bin_index(
      /* c_Real& */ x,
      /* c_Real& */ bin1_x_min,
      /* c_Real& */ bin_delta,
      /* c_Int& */ _ix_bin);
  return _ix_bin;
}
double SimUtils::bin_x_center(int ix_bin, double bin1_x_min, double bin_delta) {
  double _x_center{};
  fortran_bin_x_center(
      /* c_Int& */ ix_bin,
      /* c_Real& */ bin1_x_min,
      /* c_Real& */ bin_delta,
      /* c_Real& */ _x_center);
  return _x_center;
}
void SimUtils::n_bins_automatic(int n_data, int n) {
  fortran_n_bins_automatic(/* c_Int& */ n_data, /* c_Int& */ n);
}
void SimUtils::reallocate_spline(
    SplineProxyAllocatable1D& spline,
    int n,
    std::optional<int> n_min,
    std::optional<bool> exact) {
  // intent=inout allocatable type array
  int n_min_lvalue;
  auto* _n_min{&n_min_lvalue};
  if (n_min.has_value()) {
    n_min_lvalue = n_min.value();
  } else {
    _n_min = nullptr;
  }
  bool exact_lvalue;
  auto* _exact{&exact_lvalue};
  if (exact.has_value()) {
    exact_lvalue = exact.value();
  } else {
    _exact = nullptr;
  }
  fortran_reallocate_spline(
      /* void* */ spline.get_fortran_ptr(),
      /* c_Int& */ n,
      /* c_Int* */ _n_min,
      /* c_Bool* */ _exact);
}
SplineProxy SimUtils::create_a_spline(
    RealAllocatable1D& r0,
    RealAllocatable1D& r1,
    double slope0,
    double slope1) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  SplineProxy _spline;
  fortran_create_a_spline(
      /* void* */ r0.get_fortran_ptr(),
      /* void* */ r1.get_fortran_ptr(),
      /* c_Real& */ slope0,
      /* c_Real& */ slope1,
      /* void* */ _spline.get_fortran_ptr());
  return std::move(_spline);
}
SimUtils::SplineAkimaInterpolate SimUtils::spline_akima_interpolate(
    RealAllocatable1D& x_knot,
    RealAllocatable1D& y_knot,
    double x) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  bool _ok{};
  double _y{};
  double _dy{};
  fortran_spline_akima_interpolate(
      /* void* */ x_knot.get_fortran_ptr(),
      /* void* */ y_knot.get_fortran_ptr(),
      /* c_Real& */ x,
      /* c_Bool& */ _ok,
      /* c_Real& */ _y,
      /* c_Real& */ _dy);
  return SplineAkimaInterpolate{_ok, _y, _dy};
}
SimUtils::SplineEvaluate SimUtils::spline_evaluate(
    SplineProxyAllocatable1D& spline,
    double x) {
  // intent=in allocatable type array
  bool _ok{};
  double _y{};
  double _dy{};
  fortran_spline_evaluate(
      /* void* */ spline.get_fortran_ptr(),
      /* c_Real& */ x,
      /* c_Bool& */ _ok,
      /* c_Real& */ _y,
      /* c_Real& */ _dy);
  return SplineEvaluate{_ok, _y, _dy};
}
SimUtils::BracketIndexForSpline SimUtils::bracket_index_for_spline(
    RealAllocatable1D& x_knot,
    double x,
    std::optional<bool> strict,
    std::optional<bool> print_err) {
  // intent=in allocatable general array
  int _ix0{};
  bool strict_lvalue;
  auto* _strict{&strict_lvalue};
  if (strict.has_value()) {
    strict_lvalue = strict.value();
  } else {
    _strict = nullptr;
  }
  bool print_err_lvalue;
  auto* _print_err{&print_err_lvalue};
  if (print_err.has_value()) {
    print_err_lvalue = print_err.value();
  } else {
    _print_err = nullptr;
  }
  bool _ok{};
  fortran_bracket_index_for_spline(
      /* void* */ x_knot.get_fortran_ptr(),
      /* c_Real& */ x,
      /* c_Int& */ _ix0,
      /* c_Bool* */ _strict,
      /* c_Bool* */ _print_err,
      /* c_Bool& */ _ok);
  return BracketIndexForSpline{_ix0, _ok};
}
double SimUtils::spline1(
    SplineProxy& a_spline,
    double x,
    std::optional<int> n) {
  int n_lvalue;
  auto* _n{&n_lvalue};
  if (n.has_value()) {
    n_lvalue = n.value();
  } else {
    _n = nullptr;
  }
  double _y{};
  fortran_spline1(
      /* void* */ a_spline.get_fortran_ptr(),
      /* c_Real& */ x,
      /* c_Int* */ _n,
      /* c_Real& */ _y);
  return _y;
}
bool SimUtils::spline_akima(SplineProxyAllocatable1D& spline) {
  // intent=inout allocatable type array
  bool _ok{};
  fortran_spline_akima(/* void* */ spline.get_fortran_ptr(), /* c_Bool& */ _ok);
  return _ok;
}
void SimUtils::end_akima_spline_calc(
    SplineProxyAllocatable1D& spline,
    int which_end) {
  // intent=inout allocatable type array
  fortran_end_akima_spline_calc(
      /* void* */ spline.get_fortran_ptr(), /* c_Int& */ which_end);
}
SimUtils::ApfftCorr SimUtils::apfft_corr(
    RealAllocatable1D& rdata_in,
    std::optional<FixedArray1D<Real, 2>> bounds,
    std::string window,
    std::optional<int> diag) {
  // intent=in allocatable general array
  c_RealArr _bounds = bounds.has_value() ? bounds.value().data() : nullptr;
  auto _window = window.c_str();
  double _phase{};
  double _amp{};
  double _freq{};
  int diag_lvalue;
  auto* _diag{&diag_lvalue};
  if (diag.has_value()) {
    diag_lvalue = diag.value();
  } else {
    _diag = nullptr;
  }
  fortran_apfft_corr(
      /* void* */ rdata_in.get_fortran_ptr(),
      /* c_RealArr */ _bounds,
      /* c_Char */ _window,
      /* c_Real& */ _phase,
      /* c_Real& */ _amp,
      /* c_Real& */ _freq,
      /* c_Int* */ _diag);
  return ApfftCorr{_phase, _amp, _freq};
}
void SimUtils::apfft(
    RealAllocatable1D& rdata_in,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    optional_ref<int> diag) {
  // intent=inout allocatable general array
  auto* _bounds = bounds.data(); // CppWrapperGeneralArgument
  auto _window = window.c_str(); // ptr, inout, required
  auto* _diag = diag.has_value() ? &diag->get() : nullptr; // inout, optional
  fortran_apfft(
      /* void* */ rdata_in.get_fortran_ptr(),
      /* c_RealArr */ _bounds,
      /* c_Char */ _window,
      /* c_Real& */ phase,
      /* c_Int* */ _diag);
}
void SimUtils::apfft_ext(
    RealAllocatable1D& rdata,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    double amp,
    double freq,
    optional_ref<int> diag) {
  // intent=inout allocatable general array
  auto* _bounds = bounds.data(); // CppWrapperGeneralArgument
  auto _window = window.c_str(); // ptr, inout, required
  auto* _diag = diag.has_value() ? &diag->get() : nullptr; // inout, optional
  fortran_apfft_ext(
      /* void* */ rdata.get_fortran_ptr(),
      /* c_RealArr */ _bounds,
      /* c_Char */ _window,
      /* c_Real& */ phase,
      /* c_Real& */ amp,
      /* c_Real& */ freq,
      /* c_Int* */ _diag);
}
void SimUtils::hanhan(int N, RealAllocatable1D& hh) {
  // intent=inout allocatable general array
  fortran_hanhan(/* c_Int& */ N, /* void* */ hh.get_fortran_ptr());
}
void SimUtils::bit_set(int word, int pos, bool set_to_1) {
  fortran_bit_set(/* c_Int& */ word, /* c_Int& */ pos, /* c_Bool& */ set_to_1);
}