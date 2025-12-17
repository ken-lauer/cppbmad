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
  fortran_antiparticle(/* int& */ species, /* int& */ _anti_species);
  return _anti_species;
}
int SimUtils::species_of(double mass, int charge) {
  int _species{};
  fortran_species_of(
      /* double& */ mass, /* int& */ charge, /* int& */ _species);
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
      /* const char* */ _name,
      /* int* */ _default_,
      /* bool* */ _print_err,
      /* int& */ _species);
  return _species;
}
int SimUtils::atomic_species_id(
    int charge,
    bool is_anti,
    int atomic_num,
    int n_nuc) {
  int _species_id{};
  fortran_atomic_species_id(
      /* int& */ charge,
      /* bool& */ is_anti,
      /* int& */ atomic_num,
      /* int& */ n_nuc,
      /* int& */ _species_id);
  return _species_id;
}
std::string SimUtils::species_name(int species) {
  char _name[4096];
  fortran_species_name(/* int& */ species, /* const char* */ _name);
  return _name;
}
int SimUtils::species_id_from_openpmd(std::string pmd_name, int charge) {
  auto _pmd_name = pmd_name.c_str();
  int _species{};
  fortran_species_id_from_openpmd(
      /* const char* */ _pmd_name, /* int& */ charge, /* int& */ _species);
  return _species;
}
std::string SimUtils::openpmd_species_name(int species) {
  char _pmd_name[4096];
  fortran_openpmd_species_name(/* int& */ species, /* const char* */ _pmd_name);
  return _pmd_name;
}
double SimUtils::anomalous_moment_of(int species) {
  double _moment{};
  fortran_anomalous_moment_of(/* int& */ species, /* double& */ _moment);
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
      /* int& */ species,
      /* double* */ _non_subatomic_default,
      /* double& */ _spin);
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
      /* int& */ species, /* int* */ _default_, /* int& */ _charge);
  return _charge;
}
double SimUtils::mass_of(int species) {
  double _mass{};
  fortran_mass_of(/* int& */ species, /* double& */ _mass);
  return _mass;
}
double SimUtils::charge_to_mass_of(int species) {
  double _charge_mass_ratio{};
  fortran_charge_to_mass_of(
      /* int& */ species, /* double& */ _charge_mass_ratio);
  return _charge_mass_ratio;
}
int SimUtils::set_species_charge(int species_in, int charge) {
  int _species_charged{};
  fortran_set_species_charge(
      /* int& */ species_in, /* int& */ charge, /* int& */ _species_charged);
  return _species_charged;
}
double SimUtils::x0_radiation_length(int species) {
  double _x0{};
  fortran_x0_radiation_length(/* int& */ species, /* double& */ _x0);
  return _x0;
}
int SimUtils::atomic_number(int species) {
  int _atomic_num{};
  fortran_atomic_number(/* int& */ species, /* int& */ _atomic_num);
  return _atomic_num;
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
bool SimUtils::is_false(double param) {
  bool _this_false{};
  fortran_is_false(/* double& */ param, /* bool& */ _this_false);
  return _this_false;
}
double SimUtils::rp8(int int_in) {
  double _re_out{};
  fortran_rp8(/* int& */ int_in, /* double& */ _re_out);
  return _re_out;
}
void SimUtils::set_parameter_real(
    double& param_val,
    double& set_val,
    double& save_val) {
  fortran_set_parameter_real(
      /* double& */ param_val, /* double& */ set_val, /* double& */ save_val);
}
void SimUtils::set_parameter_int(int& param_val, int& set_val, int& save_val) {
  fortran_set_parameter_int(
      /* int& */ param_val, /* int& */ set_val, /* int& */ save_val);
}
void SimUtils::set_parameter_logic(
    bool& param_val,
    bool& set_val,
    bool& save_val) {
  fortran_set_parameter_logic(
      /* bool& */ param_val, /* bool& */ set_val, /* bool& */ save_val);
}
void SimUtils::cosc(double x, std::optional<int> nd, double& y) {
  int nd_lvalue;
  auto* _nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_cosc(/* double& */ x, /* int* */ _nd, /* double& */ y);
}
void SimUtils::asinc(double x, std::optional<int> nd, double& y) {
  int nd_lvalue;
  auto* _nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_asinc(/* double& */ x, /* int* */ _nd, /* double& */ y);
}
void SimUtils::assert_equal(
    IntAlloc1D& int_arr,
    std::string& err_str,
    int& ival) {
  // intent=in allocatable general array
  auto _err_str = err_str.c_str(); // ptr, inout, required
  fortran_assert_equal(
      /* void* */ int_arr.get_fortran_ptr(),
      /* const char* */ _err_str,
      /* int& */ ival);
}
void SimUtils::calc_file_number(
    std::string& file_name,
    int& num_in,
    int& num_out,
    bool& err_flag) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_calc_file_number(
      /* const char* */ _file_name,
      /* int& */ num_in,
      /* int& */ num_out,
      /* bool& */ err_flag);
}
void SimUtils::change_file_number(std::string& file_name, int& change) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_change_file_number(/* const char* */ _file_name, /* int& */ change);
}
void SimUtils::cos_one(double angle, double& cos1) {
  fortran_cos_one(/* double& */ angle, /* double& */ cos1);
}
void SimUtils::complex_error_function(
    double& wr,
    double& wi,
    double& zr,
    double& zi) {
  fortran_complex_error_function(
      /* double& */ wr, /* double& */ wi, /* double& */ zr, /* double& */ zi);
}
void SimUtils::cross_product(
    RealAlloc1D& a,
    RealAlloc1D& b,
    FixedArray1D<Real, 3> c) {
  // intent=in allocatable general array
  // intent=inout allocatable general array
  auto* _c = c.data(); // CppWrapperGeneralArgument
  fortran_cross_product(
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* double* */ _c);
}
void SimUtils::date_and_time_stamp(
    std::string& string,
    optional_ref<bool> numeric_month,
    optional_ref<bool> include_zone) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _numeric_month = numeric_month.has_value() ? &numeric_month->get()
                                                   : nullptr; // inout, optional
  auto* _include_zone = include_zone.has_value() ? &include_zone->get()
                                                 : nullptr; // inout, optional
  fortran_date_and_time_stamp(
      /* const char* */ _string,
      /* bool* */ _numeric_month,
      /* bool* */ _include_zone);
}
void SimUtils::detab(std::string& str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_detab(/* const char* */ _str);
}
void SimUtils::display_size_and_resolution(
    int& ix_screen,
    double& x_size,
    double& y_size,
    double& x_res,
    double& y_res) {
  fortran_display_size_and_resolution(
      /* int& */ ix_screen,
      /* double& */ x_size,
      /* double& */ y_size,
      /* double& */ x_res,
      /* double& */ y_res);
}
void SimUtils::dj_bessel(int& m, double& arg, double& dj_bes) {
  fortran_dj_bessel(/* int& */ m, /* double& */ arg, /* double& */ dj_bes);
}
void SimUtils::djb_hash(
    std::string& str,
    optional_ref<int> old_hash,
    int& hash) {
  auto _str = str.c_str(); // ptr, inout, required
  auto* _old_hash =
      old_hash.has_value() ? &old_hash->get() : nullptr; // inout, optional
  fortran_djb_hash(
      /* const char* */ _str, /* int* */ _old_hash, /* int& */ hash);
}
void SimUtils::djb_str_hash(std::string& in_str, std::string& hash_str) {
  auto _in_str = in_str.c_str(); // ptr, inout, required
  auto _hash_str = hash_str.c_str(); // ptr, inout, required
  fortran_djb_str_hash(/* const char* */ _in_str, /* const char* */ _hash_str);
}
void SimUtils::downcase_string(std::string& string) {
  auto _string = string.c_str(); // ptr, inout, required
  fortran_downcase_string(/* const char* */ _string);
}
void SimUtils::err_exit(optional_ref<std::string> err_str) {
  const char* _err_str = err_str.has_value() ? err_str->get().c_str() : nullptr;
  fortran_err_exit(/* const char* */ _err_str);
}
void SimUtils::factorial(int& n, double& fact) {
  fortran_factorial(/* int& */ n, /* double& */ fact);
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
      /* double* */ _z, /* double* */ _w, /* double* */ _dw_vec);
  vec_to_matrix(_dw_vec, dw);
}
void SimUtils::fft_1d(ComplexAlloc1D& arr, int isign) {
  // intent=inout allocatable general array
  fortran_fft_1d(/* void* */ arr.get_fortran_ptr(), /* int& */ isign);
}
void SimUtils::file_directorizer(
    std::string& in_file,
    std::string& out_file,
    std::string& directory,
    bool& add_switch) {
  auto _in_file = in_file.c_str(); // ptr, inout, required
  auto _out_file = out_file.c_str(); // ptr, inout, required
  auto _directory = directory.c_str(); // ptr, inout, required
  fortran_file_directorizer(
      /* const char* */ _in_file,
      /* const char* */ _out_file,
      /* const char* */ _directory,
      /* bool& */ add_switch);
}
void SimUtils::file_get(
    std::string& string,
    std::string& dflt_file_name,
    std::string& file_name) {
  auto _string = string.c_str(); // ptr, inout, required
  auto _dflt_file_name = dflt_file_name.c_str(); // ptr, inout, required
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_file_get(
      /* const char* */ _string,
      /* const char* */ _dflt_file_name,
      /* const char* */ _file_name);
}
void SimUtils::file_get_open(
    std::string& string,
    std::string& dflt_file_name,
    std::string& file_name,
    int& file_unit,
    bool& readonly) {
  auto _string = string.c_str(); // ptr, inout, required
  auto _dflt_file_name = dflt_file_name.c_str(); // ptr, inout, required
  auto _file_name = file_name.c_str(); // ptr, inout, required
  fortran_file_get_open(
      /* const char* */ _string,
      /* const char* */ _dflt_file_name,
      /* const char* */ _file_name,
      /* int& */ file_unit,
      /* bool& */ readonly);
}
void SimUtils::file_suffixer(
    std::string& in_file_name,
    std::string& out_file_name,
    std::string& suffix,
    bool& add_switch) {
  auto _in_file_name = in_file_name.c_str(); // ptr, inout, required
  auto _out_file_name = out_file_name.c_str(); // ptr, inout, required
  auto _suffix = suffix.c_str(); // ptr, inout, required
  fortran_file_suffixer(
      /* const char* */ _in_file_name,
      /* const char* */ _out_file_name,
      /* const char* */ _suffix,
      /* bool& */ add_switch);
}
void SimUtils::gen_complete_elliptic(
    double& kc,
    double& p,
    double& c,
    double& s,
    optional_ref<double> err_tol,
    double& value) {
  auto* _err_tol =
      err_tol.has_value() ? &err_tol->get() : nullptr; // inout, optional
  fortran_gen_complete_elliptic(
      /* double& */ kc,
      /* double& */ p,
      /* double& */ c,
      /* double& */ s,
      /* double* */ _err_tol,
      /* double& */ value);
}
void SimUtils::get_file_number(
    std::string& file_name,
    std::string& cnum_in,
    int& num_out,
    bool& err_flag) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  auto _cnum_in = cnum_in.c_str(); // ptr, inout, required
  fortran_get_file_number(
      /* const char* */ _file_name,
      /* const char* */ _cnum_in,
      /* int& */ num_out,
      /* bool& */ err_flag);
}
void SimUtils::get_file_time_stamp(std::string& file, std::string& time_stamp) {
  auto _file = file.c_str(); // ptr, inout, required
  auto _time_stamp = time_stamp.c_str(); // ptr, inout, required
  fortran_get_file_time_stamp(
      /* const char* */ _file, /* const char* */ _time_stamp);
}
void SimUtils::i_bessel(int& m, double& arg, double& i_bes) {
  fortran_i_bessel(/* int& */ m, /* double& */ arg, /* double& */ i_bes);
}
void SimUtils::i_bessel_extended(
    int& m,
    double& arg,
    std::complex<double>& i_bes) {
  fortran_i_bessel_extended(
      /* int& */ m, /* double& */ arg, /* std::complex<double>& */ i_bes);
}
void SimUtils::increment_file_number(
    std::string& file_name,
    int& digits,
    int& number,
    std::string& cnumber) {
  auto _file_name = file_name.c_str(); // ptr, inout, required
  auto _cnumber = cnumber.c_str(); // ptr, inout, required
  fortran_increment_file_number(
      /* const char* */ _file_name,
      /* int& */ digits,
      /* int& */ number,
      /* const char* */ _cnumber);
}
void SimUtils::index_nocase(
    std::string& string1,
    std::string& string2,
    int& indx) {
  auto _string1 = string1.c_str(); // ptr, inout, required
  auto _string2 = string2.c_str(); // ptr, inout, required
  fortran_index_nocase(
      /* const char* */ _string1, /* const char* */ _string2, /* int& */ indx);
}
void SimUtils::int_str(int& int_, optional_ref<int> width, std::string& str) {
  auto* _width = width.has_value() ? &width->get() : nullptr; // inout, optional
  auto _str = str.c_str(); // ptr, inout, required
  fortran_int_str(/* int& */ int_, /* int* */ _width, /* const char* */ _str);
}
void SimUtils::is_alphabetic(
    std::string& string,
    optional_ref<std::string> valid_chars,
    bool& is_alpha) {
  auto _string = string.c_str(); // ptr, inout, required
  const char* _valid_chars =
      valid_chars.has_value() ? valid_chars->get().c_str() : nullptr;
  fortran_is_alphabetic(
      /* const char* */ _string,
      /* const char* */ _valid_chars,
      /* bool& */ is_alpha);
}
void SimUtils::is_decreasing_sequence(
    RealAlloc1D& array,
    std::optional<bool> strict,
    bool& is_decreasing) {
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
      /* bool* */ _strict,
      /* bool& */ is_decreasing);
}
void SimUtils::is_increasing_sequence(
    RealAlloc1D& array,
    std::optional<bool> strict,
    bool& is_increasing) {
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
      /* bool* */ _strict,
      /* bool& */ is_increasing);
}
void SimUtils::is_integer(
    std::string& string,
    optional_ref<int> int_,
    optional_ref<std::string> delims,
    optional_ref<int> ix_word,
    bool& valid) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _int_ = int_.has_value() ? &int_->get() : nullptr; // inout, optional
  const char* _delims = delims.has_value() ? delims->get().c_str() : nullptr;
  auto* _ix_word =
      ix_word.has_value() ? &ix_word->get() : nullptr; // inout, optional
  fortran_is_integer(
      /* const char* */ _string,
      /* int* */ _int_,
      /* const char* */ _delims,
      /* int* */ _ix_word,
      /* bool& */ valid);
}
void SimUtils::is_logical(
    std::string& string,
    optional_ref<bool> ignore,
    bool& valid) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _ignore =
      ignore.has_value() ? &ignore->get() : nullptr; // inout, optional
  fortran_is_logical(
      /* const char* */ _string, /* bool* */ _ignore, /* bool& */ valid);
}
void SimUtils::is_real(
    std::string& string,
    optional_ref<bool> ignore,
    optional_ref<double> real_num,
    bool& valid) {
  auto _string = string.c_str(); // ptr, inout, required
  auto* _ignore =
      ignore.has_value() ? &ignore->get() : nullptr; // inout, optional
  auto* _real_num =
      real_num.has_value() ? &real_num->get() : nullptr; // inout, optional
  fortran_is_real(
      /* const char* */ _string,
      /* bool* */ _ignore,
      /* double* */ _real_num,
      /* bool& */ valid);
}
void SimUtils::j_bessel(int& m, double& arg, double& j_bes) {
  fortran_j_bessel(/* int& */ m, /* double& */ arg, /* double& */ j_bes);
}
void SimUtils::linear_fit(
    RealAlloc1D& x,
    RealAlloc1D& y,
    int& n_data,
    double& a,
    double& b,
    double& sig_a,
    double& sig_b) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  fortran_linear_fit(
      /* void* */ x.get_fortran_ptr(),
      /* void* */ y.get_fortran_ptr(),
      /* int& */ n_data,
      /* double& */ a,
      /* double& */ b,
      /* double& */ sig_a,
      /* double& */ sig_b);
}
FixedArray1D<Real, 3> SimUtils::linear_fit_2d(
    RealAlloc1D& x,
    RealAlloc1D& y,
    RealAlloc1D& z) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=in allocatable general array
  FixedArray1D<Real, 3> _coef;
  fortran_linear_fit_2d(
      /* void* */ x.get_fortran_ptr(),
      /* void* */ y.get_fortran_ptr(),
      /* void* */ z.get_fortran_ptr(),
      /* double* */ _coef.data());
  return _coef;
}
void SimUtils::logic_str(bool& logic, std::string& str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_logic_str(/* bool& */ logic, /* const char* */ _str);
}
int SimUtils::lunget() {
  int _func_retval__{};
  fortran_lunget(/* int& */ _func_retval__);
  return _func_retval__;
}
void SimUtils::match_reg(std::string& str, std::string& pat, bool& is_match) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _pat = pat.c_str(); // ptr, inout, required
  fortran_match_reg(
      /* const char* */ _str, /* const char* */ _pat, /* bool& */ is_match);
}
void SimUtils::milli_sleep(int& milli_sec) {
  fortran_milli_sleep(/* int& */ milli_sec);
}
void SimUtils::make_legal_comment(
    std::string& comment_in,
    std::string& comment_out) {
  auto _comment_in = comment_in.c_str(); // ptr, inout, required
  auto _comment_out = comment_out.c_str(); // ptr, inout, required
  fortran_make_legal_comment(
      /* const char* */ _comment_in, /* const char* */ _comment_out);
}
void SimUtils::match_wild(
    std::string& string,
    std::string& template_,
    bool& is_match) {
  auto _string = string.c_str(); // ptr, inout, required
  auto _template_ = template_.c_str(); // ptr, inout, required
  fortran_match_wild(
      /* const char* */ _string,
      /* const char* */ _template_,
      /* bool& */ is_match);
}
void SimUtils::n_choose_k(int& n, int& k, double& nck) {
  fortran_n_choose_k(/* int& */ n, /* int& */ k, /* double& */ nck);
}
RealAlloc1D SimUtils::n_spline_create(
    RealAlloc1D& deriv0,
    RealAlloc1D& deriv1,
    double x1) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  // intent=out allocatable general array
  auto n_spline{RealAlloc1D()};
  fortran_n_spline_create(
      /* void* */ deriv0.get_fortran_ptr(),
      /* void* */ deriv1.get_fortran_ptr(),
      /* double& */ x1,
      /* void* */ n_spline.get_fortran_ptr());
  return std::move(n_spline);
}
void SimUtils::nametable_add(
    NametableProxy& nametable,
    std::string& name,
    int& ix_name) {
  auto _name = name.c_str(); // ptr, inout, required
  fortran_nametable_add(
      /* void* */ nametable.get_fortran_ptr(),
      /* const char* */ _name,
      /* int& */ ix_name);
}
void SimUtils::nametable_bracket_indexx(
    NametableProxy& nametable,
    std::string& name,
    optional_ref<int> n_match,
    int& ix_max) {
  auto _name = name.c_str(); // ptr, inout, required
  auto* _n_match =
      n_match.has_value() ? &n_match->get() : nullptr; // inout, optional
  fortran_nametable_bracket_indexx(
      /* void* */ nametable.get_fortran_ptr(),
      /* const char* */ _name,
      /* int* */ _n_match,
      /* int& */ ix_max);
}
void SimUtils::nametable_change1(
    NametableProxy& nametable,
    std::string& name,
    int& ix_name) {
  auto _name = name.c_str(); // ptr, inout, required
  fortran_nametable_change1(
      /* void* */ nametable.get_fortran_ptr(),
      /* const char* */ _name,
      /* int& */ ix_name);
}
void SimUtils::nametable_init(
    NametableProxy& nametable,
    optional_ref<int> n_min,
    optional_ref<int> n_max) {
  auto* _n_min = n_min.has_value() ? &n_min->get() : nullptr; // inout, optional
  auto* _n_max = n_max.has_value() ? &n_max->get() : nullptr; // inout, optional
  fortran_nametable_init(
      /* void* */ nametable.get_fortran_ptr(),
      /* int* */ _n_min,
      /* int* */ _n_max);
}
void SimUtils::nametable_remove(NametableProxy& nametable, int& ix_name) {
  fortran_nametable_remove(
      /* void* */ nametable.get_fortran_ptr(), /* int& */ ix_name);
}
void SimUtils::ordinal_str(int& n, std::string& str) {
  auto _str = str.c_str(); // ptr, inout, required
  fortran_ordinal_str(/* int& */ n, /* const char* */ _str);
}
void SimUtils::parse_fortran_format(
    std::string& format_str,
    int& n_repeat,
    int& power,
    std::string& descrip,
    int& width,
    int& digits) {
  auto _format_str = format_str.c_str(); // ptr, inout, required
  auto _descrip = descrip.c_str(); // ptr, inout, required
  fortran_parse_fortran_format(
      /* const char* */ _format_str,
      /* int& */ n_repeat,
      /* int& */ power,
      /* const char* */ _descrip,
      /* int& */ width,
      /* int& */ digits);
}
void SimUtils::poly_eval(
    RealAlloc1D& poly,
    double x,
    std::optional<bool> diff_coef,
    double& y) {
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
      /* double& */ x,
      /* bool* */ _diff_coef,
      /* double& */ y);
}
void SimUtils::probability_funct(double x, double& prob) {
  fortran_probability_funct(/* double& */ x, /* double& */ prob);
}
void SimUtils::quadratic_roots(
    FixedArray1D<Real, 3> coefs,
    FixedArray1D<Complex, 2> root) {
  auto* _coefs = coefs.data(); // CppWrapperGeneralArgument
  auto* _root = root.data(); // CppWrapperGeneralArgument
  fortran_quadratic_roots(
      /* double* */ _coefs, /* std::complex<double>* */ _root);
}
void SimUtils::query_string(
    std::string& query_str,
    bool& upcase,
    std::string& return_str,
    int& ix,
    int& ios) {
  auto _query_str = query_str.c_str(); // ptr, inout, required
  auto _return_str = return_str.c_str(); // ptr, inout, required
  fortran_query_string(
      /* const char* */ _query_str,
      /* bool& */ upcase,
      /* const char* */ _return_str,
      /* int& */ ix,
      /* int& */ ios);
}
void SimUtils::quote(std::string& str, std::string& q_str) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _q_str = q_str.c_str(); // ptr, inout, required
  fortran_quote(/* const char* */ _str, /* const char* */ _q_str);
}
void SimUtils::real_to_string(
    double& real_num,
    int& width,
    optional_ref<int> n_signif,
    optional_ref<int> n_decimal,
    std::string& str) {
  auto* _n_signif =
      n_signif.has_value() ? &n_signif->get() : nullptr; // inout, optional
  auto* _n_decimal =
      n_decimal.has_value() ? &n_decimal->get() : nullptr; // inout, optional
  auto _str = str.c_str(); // ptr, inout, required
  fortran_real_to_string(
      /* double& */ real_num,
      /* int& */ width,
      /* int* */ _n_signif,
      /* int* */ _n_decimal,
      /* const char* */ _str);
}
void SimUtils::real_num_fortran_format(
    double& number,
    int& width,
    optional_ref<int> n_blanks,
    std::string& fmt_str) {
  auto* _n_blanks =
      n_blanks.has_value() ? &n_blanks->get() : nullptr; // inout, optional
  auto _fmt_str = fmt_str.c_str(); // ptr, inout, required
  fortran_real_num_fortran_format(
      /* double& */ number,
      /* int& */ width,
      /* int* */ _n_blanks,
      /* const char* */ _fmt_str);
}
void SimUtils::str_count(std::string& str, std::string& match, int& num) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _match = match.c_str(); // ptr, inout, required
  fortran_str_count(
      /* const char* */ _str, /* const char* */ _match, /* int& */ num);
}
void SimUtils::real_path(
    std::string& path_in,
    std::string& path_out,
    bool& is_ok) {
  auto _path_in = path_in.c_str(); // ptr, inout, required
  auto _path_out = path_out.c_str(); // ptr, inout, required
  fortran_real_path(
      /* const char* */ _path_in,
      /* const char* */ _path_out,
      /* bool& */ is_ok);
}
void SimUtils::real_str(
    double& r_num,
    optional_ref<int> n_signif,
    optional_ref<int> n_decimal,
    std::string& str) {
  auto* _n_signif =
      n_signif.has_value() ? &n_signif->get() : nullptr; // inout, optional
  auto* _n_decimal =
      n_decimal.has_value() ? &n_decimal->get() : nullptr; // inout, optional
  auto _str = str.c_str(); // ptr, inout, required
  fortran_real_str(
      /* double& */ r_num,
      /* int* */ _n_signif,
      /* int* */ _n_decimal,
      /* const char* */ _str);
}
double SimUtils::rms_value(
    RealAlloc1D& val_arr,
    optional_ref<BoolAlloc1D> good_val,
    double& rms_val) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  auto* _good_val = good_val.has_value() ? good_val->get().get_fortran_ptr()
                                         : nullptr; // input, optional
  double _ave_val{};
  fortran_rms_value(
      /* void* */ val_arr.get_fortran_ptr(),
      /* void* */ _good_val,
      /* double& */ _ave_val,
      /* double& */ rms_val);
  return _ave_val;
}
void SimUtils::rot_2d(
    FixedArray1D<Real, 2> vec_in,
    double angle,
    FixedArray1D<Real, 2> vec_out) {
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  auto* _vec_out = vec_out.data(); // CppWrapperGeneralArgument
  fortran_rot_2d(
      /* double* */ _vec_in, /* double& */ angle, /* double* */ _vec_out);
}
void SimUtils::run_timer(
    std::string& command,
    optional_ref<double> time,
    optional_ref<double> time0) {
  auto _command = command.c_str(); // ptr, inout, required
  auto* _time = time.has_value() ? &time->get() : nullptr; // inout, optional
  auto* _time0 = time0.has_value() ? &time0->get() : nullptr; // inout, optional
  fortran_run_timer(
      /* const char* */ _command, /* double* */ _time, /* double* */ _time0);
}
void SimUtils::sinc(double x, std::optional<int> nd, double& y) {
  int nd_lvalue;
  auto* _nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_sinc(/* double& */ x, /* int* */ _nd, /* double& */ y);
}
void SimUtils::sincc(double x, std::optional<int> nd, double& y) {
  int nd_lvalue;
  auto* _nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_sincc(/* double& */ x, /* int* */ _nd, /* double& */ y);
}
void SimUtils::sinhx_x(double x, std::optional<int> nd, double& y) {
  int nd_lvalue;
  auto* _nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_sinhx_x(/* double& */ x, /* int* */ _nd, /* double& */ y);
}
void SimUtils::skip_header(int& ix_unit, bool& error_flag) {
  fortran_skip_header(/* int& */ ix_unit, /* bool& */ error_flag);
}
void SimUtils::sqrt_one(double x, std::optional<int> nd, double& ds1) {
  int nd_lvalue;
  auto* _nd{&nd_lvalue};
  if (nd.has_value()) {
    nd_lvalue = nd.value();
  } else {
    _nd = nullptr;
  }
  fortran_sqrt_one(/* double& */ x, /* int* */ _nd, /* double& */ ds1);
}
void SimUtils::sqrt_alpha(double alpha, double x, double& y) {
  fortran_sqrt_alpha(/* double& */ alpha, /* double& */ x, /* double& */ y);
}
void SimUtils::str_first_in_set(
    std::string& line,
    std::string& set,
    optional_ref<bool> ignore_clauses,
    int& ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  auto* _ignore_clauses = ignore_clauses.has_value()
      ? &ignore_clauses->get()
      : nullptr; // inout, optional
  fortran_str_first_in_set(
      /* const char* */ _line,
      /* const char* */ _set,
      /* bool* */ _ignore_clauses,
      /* int& */ ix_match);
}
void SimUtils::str_first_not_in_set(
    std::string& line,
    std::string& set,
    int& ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  fortran_str_first_not_in_set(
      /* const char* */ _line, /* const char* */ _set, /* int& */ ix_match);
}
void SimUtils::str_last_in_set(
    std::string& line,
    std::string& set,
    int& ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  fortran_str_last_in_set(
      /* const char* */ _line, /* const char* */ _set, /* int& */ ix_match);
}
void SimUtils::str_last_not_in_set(
    std::string& line,
    std::string& set,
    int& ix_match) {
  auto _line = line.c_str(); // ptr, inout, required
  auto _set = set.c_str(); // ptr, inout, required
  fortran_str_last_not_in_set(
      /* const char* */ _line, /* const char* */ _set, /* int& */ ix_match);
}
void SimUtils::string_to_int(
    std::string& line,
    int& default_,
    bool& err_flag,
    optional_ref<bool> err_print_flag,
    int& value) {
  auto _line = line.c_str(); // ptr, inout, required
  auto* _err_print_flag = err_print_flag.has_value()
      ? &err_print_flag->get()
      : nullptr; // inout, optional
  fortran_string_to_int(
      /* const char* */ _line,
      /* int& */ default_,
      /* bool& */ err_flag,
      /* bool* */ _err_print_flag,
      /* int& */ value);
}
void SimUtils::string_to_real(
    std::string& line,
    double& default_,
    bool& err_flag,
    optional_ref<bool> err_print_flag,
    double& value) {
  auto _line = line.c_str(); // ptr, inout, required
  auto* _err_print_flag = err_print_flag.has_value()
      ? &err_print_flag->get()
      : nullptr; // inout, optional
  fortran_string_to_real(
      /* const char* */ _line,
      /* double& */ default_,
      /* bool& */ err_flag,
      /* bool* */ _err_print_flag,
      /* double& */ value);
}
void SimUtils::string_trim2(
    std::string& in_str,
    std::string& delimitors,
    std::string& out_str,
    int& ix_word,
    std::string& delim,
    int& ix_next) {
  auto _in_str = in_str.c_str(); // ptr, inout, required
  auto _delimitors = delimitors.c_str(); // ptr, inout, required
  auto _out_str = out_str.c_str(); // ptr, inout, required
  auto _delim = delim.c_str(); // ptr, inout, required
  fortran_string_trim2(
      /* const char* */ _in_str,
      /* const char* */ _delimitors,
      /* const char* */ _out_str,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* int& */ ix_next);
}
void SimUtils::to_str(
    double& num,
    optional_ref<int> max_signif,
    std::string& string) {
  auto* _max_signif =
      max_signif.has_value() ? &max_signif->get() : nullptr; // inout, optional
  auto _string = string.c_str(); // ptr, inout, required
  fortran_to_str(
      /* double& */ num, /* int* */ _max_signif, /* const char* */ _string);
}
void SimUtils::type_this_file(std::string& filename) {
  auto _filename = filename.c_str(); // ptr, inout, required
  fortran_type_this_file(/* const char* */ _filename);
}
void SimUtils::upcase_string(std::string& string) {
  auto _string = string.c_str(); // ptr, inout, required
  fortran_upcase_string(/* const char* */ _string);
}
void SimUtils::word_len(std::string& wording, int& wlen) {
  auto _wording = wording.c_str(); // ptr, inout, required
  fortran_word_len(/* const char* */ _wording, /* int& */ wlen);
}
void SimUtils::word_read(
    std::string& in_str,
    std::string& delim_list,
    std::string& word,
    int& ix_word,
    std::string& delim,
    bool& delim_found,
    std::string& out_str,
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
      /* const char* */ _in_str,
      /* const char* */ _delim_list,
      /* const char* */ _word,
      /* int& */ ix_word,
      /* const char* */ _delim,
      /* bool& */ delim_found,
      /* const char* */ _out_str,
      /* bool* */ _ignore_interior);
}
void SimUtils::str_substitute(
    std::string& string,
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
      /* const char* */ _string,
      /* const char* */ _str_match,
      /* const char* */ _str_replace,
      /* bool* */ _do_trim,
      /* bool* */ _ignore_escaped);
}
void SimUtils::str_match_wild(
    std::string& str,
    std::string& pat,
    bool& a_match) {
  auto _str = str.c_str(); // ptr, inout, required
  auto _pat = pat.c_str(); // ptr, inout, required
  fortran_str_match_wild(
      /* const char* */ _str, /* const char* */ _pat, /* bool& */ a_match);
}
std::string SimUtils::str_upcase(std::string src) {
  char _dst[4096];
  auto _src = src.c_str();
  fortran_str_upcase(/* const char* */ _dst, /* const char* */ _src);
  return _dst;
}
std::string SimUtils::str_downcase(std::string src) {
  char _dst[4096];
  auto _src = src.c_str();
  fortran_str_downcase(/* const char* */ _dst, /* const char* */ _src);
  return _dst;
}
void SimUtils::system_command(std::string& line, optional_ref<bool> err_flag) {
  auto _line = line.c_str(); // ptr, inout, required
  auto* _err_flag =
      err_flag.has_value() ? &err_flag->get() : nullptr; // inout, optional
  fortran_system_command(/* const char* */ _line, /* bool* */ _err_flag);
}
void SimUtils::string_trim(
    std::string& in_string,
    std::string& out_string,
    int& word_len) {
  auto _in_string = in_string.c_str(); // ptr, inout, required
  auto _out_string = out_string.c_str(); // ptr, inout, required
  fortran_string_trim(
      /* const char* */ _in_string,
      /* const char* */ _out_string,
      /* int& */ word_len);
}
int SimUtils::virtual_memory_usage() {
  int _usage{};
  fortran_virtual_memory_usage(/* int& */ _usage);
  return _usage;
}
void SimUtils::find_location_real(
    RealAlloc1D& arr,
    double value,
    int& ix_match) {
  // intent=in allocatable general array
  fortran_find_location_real(
      /* void* */ arr.get_fortran_ptr(),
      /* double& */ value,
      /* int& */ ix_match);
}
void SimUtils::find_location_int(IntAlloc1D& arr, int& value, int& ix_match) {
  // intent=inout allocatable general array
  fortran_find_location_int(
      /* void* */ arr.get_fortran_ptr(), /* int& */ value, /* int& */ ix_match);
}
void SimUtils::find_location_logic(
    BoolAlloc1D& arr,
    bool& value,
    int& ix_match) {
  // intent=inout allocatable general array
  fortran_find_location_logic(
      /* void* */ arr.get_fortran_ptr(),
      /* bool& */ value,
      /* int& */ ix_match);
}
double SimUtils::coarse_frequency_estimate(
    RealAlloc1D& data,
    optional_ref<bool> error) {
  // intent=in allocatable general array
  auto* _error = error.has_value() ? &error->get() : nullptr; // inout, optional
  double _frequency{};
  fortran_coarse_frequency_estimate(
      /* void* */ data.get_fortran_ptr(),
      /* bool* */ _error,
      /* double& */ _frequency);
  return _frequency;
}
double SimUtils::fine_frequency_estimate(RealAlloc1D& data) {
  // intent=in allocatable general array
  double _frequency{};
  fortran_fine_frequency_estimate(
      /* void* */ data.get_fortran_ptr(), /* double& */ _frequency);
  return _frequency;
}
SimUtils::FourierAmplitude SimUtils::fourier_amplitude(
    RealAlloc1D& data,
    double frequency) {
  // intent=in allocatable general array
  double _cos_amp{};
  double _sin_amp{};
  double _dcos_amp{};
  double _dsin_amp{};
  fortran_fourier_amplitude(
      /* void* */ data.get_fortran_ptr(),
      /* double& */ frequency,
      /* double& */ _cos_amp,
      /* double& */ _sin_amp,
      /* double& */ _dcos_amp,
      /* double& */ _dsin_amp);
  return FourierAmplitude{_cos_amp, _sin_amp, _dcos_amp, _dsin_amp};
}
SimUtils::WMatToAxisAngle SimUtils::w_mat_to_axis_angle(
    FixedArray2D<Real, 3, 3> w_mat) {
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  FixedArray1D<Real, 3> _axis;
  double _angle{};
  fortran_w_mat_to_axis_angle(
      /* double* */ _w_mat_vec,
      /* double* */ _axis.data(),
      /* double& */ _angle);
  return WMatToAxisAngle{_axis, _angle};
}
FixedArray1D<Real, 4> SimUtils::w_mat_to_quat(FixedArray2D<Real, 3, 3> w_mat) {
  double _w_mat_vec[3 * 3];
  matrix_to_vec(w_mat, _w_mat_vec);
  FixedArray1D<Real, 4> _quat;
  fortran_w_mat_to_quat(/* double* */ _w_mat_vec, /* double* */ _quat.data());
  return _quat;
}
FixedArray2D<Real, 3, 3> SimUtils::quat_to_w_mat(FixedArray1D<Real, 4> quat) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  FixedArray2D<Real, 3, 3> w_mat;
  double _w_mat_vec[3 * 3];
  fortran_quat_to_w_mat(/* double* */ _quat, /* double* */ _w_mat_vec);
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
      /* double* */ _axis, /* double& */ angle, /* double* */ _w_mat_vec);
  vec_to_matrix(_w_mat_vec, w_mat);
  return w_mat;
}
FixedArray1D<Real, 3> SimUtils::quat_to_omega(FixedArray1D<Real, 4> quat) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _omega;
  fortran_quat_to_omega(/* double* */ _quat, /* double* */ _omega.data());
  return _omega;
}
FixedArray1D<Real, 4> SimUtils::omega_to_quat(FixedArray1D<Real, 3> omega) {
  auto* _omega = omega.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _quat;
  fortran_omega_to_quat(/* double* */ _omega, /* double* */ _quat.data());
  return _quat;
}
SimUtils::QuatToAxisAngle SimUtils::quat_to_axis_angle(
    FixedArray1D<Real, 4> quat) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _axis;
  double _angle{};
  fortran_quat_to_axis_angle(
      /* double* */ _quat, /* double* */ _axis.data(), /* double& */ _angle);
  return QuatToAxisAngle{_axis, _angle};
}
FixedArray1D<Real, 4> SimUtils::axis_angle_to_quat(
    FixedArray1D<Real, 3> axis,
    double angle) {
  auto* _axis = axis.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _quat;
  fortran_axis_angle_to_quat(
      /* double* */ _axis, /* double& */ angle, /* double* */ _quat.data());
  return _quat;
}
FixedArray1D<Real, 4> SimUtils::quat_conj_real(FixedArray1D<Real, 4> q_in) {
  auto* _q_in = q_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _q_out;
  fortran_quat_conj_real(/* double* */ _q_in, /* double* */ _q_out.data());
  return _q_out;
}
FixedArray1D<Complex, 4> SimUtils::quat_conj_complex(
    FixedArray1D<Complex, 4> q_in) {
  auto* _q_in = q_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Complex, 4> _q_out;
  fortran_quat_conj_complex(
      /* std::complex<double>* */ _q_in,
      /* std::complex<double>* */ _q_out.data());
  return _q_out;
}
FixedArray1D<Real, 4> SimUtils::quat_inverse(FixedArray1D<Real, 4> q_in) {
  auto* _q_in = q_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 4> _q_out;
  fortran_quat_inverse(/* double* */ _q_in, /* double* */ _q_out.data());
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
  double* _q3 = q3.has_value() ? q3.value().data() : nullptr;
  double* _q4 = q4.has_value() ? q4.value().data() : nullptr;
  double* _q5 = q5.has_value() ? q5.value().data() : nullptr;
  double* _q6 = q6.has_value() ? q6.value().data() : nullptr;
  double* _q7 = q7.has_value() ? q7.value().data() : nullptr;
  double* _q8 = q8.has_value() ? q8.value().data() : nullptr;
  double* _q9 = q9.has_value() ? q9.value().data() : nullptr;
  FixedArray1D<Real, 4> _q_out;
  fortran_quat_mul_real(
      /* double* */ _q1,
      /* double* */ _q2,
      /* double* */ _q3,
      /* double* */ _q4,
      /* double* */ _q5,
      /* double* */ _q6,
      /* double* */ _q7,
      /* double* */ _q8,
      /* double* */ _q9,
      /* double* */ _q_out.data());
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
  std::complex<double>* _q3 = q3.has_value() ? q3.value().data() : nullptr;
  std::complex<double>* _q4 = q4.has_value() ? q4.value().data() : nullptr;
  std::complex<double>* _q5 = q5.has_value() ? q5.value().data() : nullptr;
  std::complex<double>* _q6 = q6.has_value() ? q6.value().data() : nullptr;
  std::complex<double>* _q7 = q7.has_value() ? q7.value().data() : nullptr;
  std::complex<double>* _q8 = q8.has_value() ? q8.value().data() : nullptr;
  std::complex<double>* _q9 = q9.has_value() ? q9.value().data() : nullptr;
  FixedArray1D<Complex, 4> _q_out;
  fortran_quat_mul_complex(
      /* std::complex<double>* */ _q1,
      /* std::complex<double>* */ _q2,
      /* std::complex<double>* */ _q3,
      /* std::complex<double>* */ _q4,
      /* std::complex<double>* */ _q5,
      /* std::complex<double>* */ _q6,
      /* std::complex<double>* */ _q7,
      /* std::complex<double>* */ _q8,
      /* std::complex<double>* */ _q9,
      /* std::complex<double>* */ _q_out.data());
  return _q_out;
}
FixedArray1D<Real, 3> SimUtils::quat_rotate_real(
    FixedArray1D<Real, 4> quat,
    FixedArray1D<Real, 3> vec_in) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Real, 3> _vec_out;
  fortran_quat_rotate_real(
      /* double* */ _quat,
      /* double* */ _vec_in,
      /* double* */ _vec_out.data());
  return _vec_out;
}
FixedArray1D<Complex, 3> SimUtils::quat_rotate_complex(
    FixedArray1D<Complex, 4> quat,
    FixedArray1D<Complex, 3> vec_in) {
  auto* _quat = quat.data(); // CppWrapperGeneralArgument
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  FixedArray1D<Complex, 3> _vec_out;
  fortran_quat_rotate_complex(
      /* std::complex<double>* */ _quat,
      /* std::complex<double>* */ _vec_in,
      /* std::complex<double>* */ _vec_out.data());
  return _vec_out;
}
FixedArray1D<Real, 3> SimUtils::rotate_vec_given_axis_angle(
    FixedArray1D<Real, 3> vec_in,
    RealAlloc1D& axis,
    double angle) {
  auto* _vec_in = vec_in.data(); // CppWrapperGeneralArgument
  // intent=in allocatable general array
  FixedArray1D<Real, 3> _vec_out;
  fortran_rotate_vec_given_axis_angle(
      /* double* */ _vec_in,
      /* void* */ axis.get_fortran_ptr(),
      /* double& */ angle,
      /* double* */ _vec_out.data());
  return _vec_out;
}
void SimUtils::rotate_vec(RealAlloc1D& vec, int axis, double angle) {
  // intent=inout allocatable general array
  fortran_rotate_vec(
      /* void* */ vec.get_fortran_ptr(), /* int& */ axis, /* double& */ angle);
}
void SimUtils::naff(
    ComplexAlloc1D& cdata,
    RealAlloc1D& freqs,
    ComplexAlloc1D& amps,
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
      /* int* */ _opt_dump_spectra,
      /* bool* */ _opt_zero_first);
}
void SimUtils::projdd(
    ComplexAlloc1D& a,
    ComplexAlloc1D& b,
    std::complex<double>& func_retval__) {
  // intent=inout allocatable general array
  // intent=inout allocatable general array
  fortran_projdd(
      /* void* */ a.get_fortran_ptr(),
      /* void* */ b.get_fortran_ptr(),
      /* std::complex<double>& */ func_retval__);
}
void SimUtils::maximize_projection(
    double& seed,
    ComplexAlloc1D& cdata,
    double& func_retval__) {
  // intent=inout allocatable general array
  fortran_maximize_projection(
      /* double& */ seed,
      /* void* */ cdata.get_fortran_ptr(),
      /* double& */ func_retval__);
}
void SimUtils::interpolated_fft_gsl(
    ComplexAlloc1D& cdata,
    bool& calc_ok,
    optional_ref<int> opt_dump_spectrum,
    optional_ref<int> opt_dump_index,
    double& this_fft) {
  // intent=inout allocatable general array
  auto* _opt_dump_spectrum = opt_dump_spectrum.has_value()
      ? &opt_dump_spectrum->get()
      : nullptr; // inout, optional
  auto* _opt_dump_index = opt_dump_index.has_value()
      ? &opt_dump_index->get()
      : nullptr; // inout, optional
  fortran_interpolated_fft_gsl(
      /* void* */ cdata.get_fortran_ptr(),
      /* bool& */ calc_ok,
      /* int* */ _opt_dump_spectrum,
      /* int* */ _opt_dump_index,
      /* double& */ this_fft);
}
void SimUtils::interpolated_fft(
    ComplexAlloc1D& cdata,
    bool& calc_ok,
    optional_ref<int> opt_dump_spectrum,
    optional_ref<int> opt_dump_index,
    double& this_fft) {
  // intent=inout allocatable general array
  auto* _opt_dump_spectrum = opt_dump_spectrum.has_value()
      ? &opt_dump_spectrum->get()
      : nullptr; // inout, optional
  auto* _opt_dump_index = opt_dump_index.has_value()
      ? &opt_dump_index->get()
      : nullptr; // inout, optional
  fortran_interpolated_fft(
      /* void* */ cdata.get_fortran_ptr(),
      /* bool& */ calc_ok,
      /* int* */ _opt_dump_spectrum,
      /* int* */ _opt_dump_index,
      /* double& */ this_fft);
}
int SimUtils::initfixedwindowls(int N, double dt, int order, int der) {
  int _id{};
  fortran_initfixedwindowls(
      /* int& */ N,
      /* double& */ dt,
      /* int& */ order,
      /* int& */ der,
      /* int& */ _id);
  return _id;
}
void SimUtils::destfixedwindowls(int id) {
  fortran_destfixedwindowls(/* int& */ id);
}
void SimUtils::fixedwindowls(double ynew, int id, double& z) {
  fortran_fixedwindowls(/* double& */ ynew, /* int& */ id, /* double& */ z);
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
      /* double* */ _y,
      /* double* */ _y1,
      /* double* */ _y2,
      /* double* */ _y12,
      /* double& */ x1l,
      /* double& */ x1u,
      /* double& */ x2l,
      /* double& */ x2u,
      /* double& */ x1,
      /* double& */ x2,
      /* double& */ _ansy,
      /* double& */ _ansy1,
      /* double& */ _ansy2);
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
      /* double* */ _y,
      /* double* */ _y1,
      /* double* */ _y2,
      /* double* */ _y12,
      /* double& */ d1,
      /* double& */ d2,
      /* double* */ _c_vec);
  vec_to_matrix(_c_vec, c);
  return c;
}
void SimUtils::super_sort(IntAlloc1D& arr) {
  // intent=inout allocatable general array
  fortran_super_sort(/* void* */ arr.get_fortran_ptr());
}
SimUtils::SuperPolint SimUtils::super_polint(
    RealAlloc1D& xa,
    RealAlloc1D& ya,
    double x) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  double _y{};
  double _dy{};
  fortran_super_polint(
      /* void* */ xa.get_fortran_ptr(),
      /* void* */ ya.get_fortran_ptr(),
      /* double& */ x,
      /* double& */ _y,
      /* double& */ _dy);
  return SuperPolint{_y, _dy};
}
double SimUtils::super_poly(double x, RealAlloc1D& coeffs) {
  // intent=in allocatable general array
  double _value{};
  fortran_super_poly(
      /* double& */ x,
      /* void* */ coeffs.get_fortran_ptr(),
      /* double& */ _value);
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
  fortran_ran_seed_put(/* int& */ seed, /* int* */ _mpi_offset);
}
int SimUtils::ran_seed_get() {
  int _seed{};
  fortran_ran_seed_get(/* int& */ _seed);
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
      /* double& */ x_norm,
      /* double& */ y_norm,
      /* void* */ bi_coef.get_fortran_ptr(),
      /* std::complex<double>& */ _df_dx,
      /* std::complex<double>& */ _df_dy,
      /* std::complex<double>& */ _f_val);
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
      /* double& */ x_norm,
      /* double& */ y_norm,
      /* double& */ z_norm,
      /* void* */ tri_coef.get_fortran_ptr(),
      /* std::complex<double>& */ _df_dx,
      /* std::complex<double>& */ _df_dy,
      /* std::complex<double>& */ _df_dz,
      /* std::complex<double>& */ _f_val);
  return TricubicCmplxEval{_df_dx, _df_dy, _df_dz, _f_val};
}
int SimUtils::bin_index(double x, double bin1_x_min, double bin_delta) {
  int _ix_bin{};
  fortran_bin_index(
      /* double& */ x,
      /* double& */ bin1_x_min,
      /* double& */ bin_delta,
      /* int& */ _ix_bin);
  return _ix_bin;
}
double SimUtils::bin_x_center(
    int& ix_bin,
    double bin1_x_min,
    double bin_delta) {
  double _x_center{};
  fortran_bin_x_center(
      /* int& */ ix_bin,
      /* double& */ bin1_x_min,
      /* double& */ bin_delta,
      /* double& */ _x_center);
  return _x_center;
}
void SimUtils::n_bins_automatic(int& n_data, int& n) {
  fortran_n_bins_automatic(/* int& */ n_data, /* int& */ n);
}
void SimUtils::reallocate_spline(
    SplineProxyAlloc1D& spline,
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
      /* int& */ n,
      /* int* */ _n_min,
      /* bool* */ _exact);
}
SplineProxy SimUtils::create_a_spline(
    RealAlloc1D& r0,
    RealAlloc1D& r1,
    double slope0,
    double slope1) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  SplineProxy _spline;
  fortran_create_a_spline(
      /* void* */ r0.get_fortran_ptr(),
      /* void* */ r1.get_fortran_ptr(),
      /* double& */ slope0,
      /* double& */ slope1,
      /* void* */ _spline.get_fortran_ptr());
  return std::move(_spline);
}
SimUtils::SplineAkimaInterpolate SimUtils::spline_akima_interpolate(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    double x) {
  // intent=in allocatable general array
  // intent=in allocatable general array
  bool _ok{};
  double _y{};
  double _dy{};
  fortran_spline_akima_interpolate(
      /* void* */ x_knot.get_fortran_ptr(),
      /* void* */ y_knot.get_fortran_ptr(),
      /* double& */ x,
      /* bool& */ _ok,
      /* double& */ _y,
      /* double& */ _dy);
  return SplineAkimaInterpolate{_ok, _y, _dy};
}
SimUtils::SplineEvaluate SimUtils::spline_evaluate(
    SplineProxyAlloc1D& spline,
    double x) {
  // intent=in allocatable type array
  bool _ok{};
  double _y{};
  double _dy{};
  fortran_spline_evaluate(
      /* void* */ spline.get_fortran_ptr(),
      /* double& */ x,
      /* bool& */ _ok,
      /* double& */ _y,
      /* double& */ _dy);
  return SplineEvaluate{_ok, _y, _dy};
}
SimUtils::BracketIndexForSpline SimUtils::bracket_index_for_spline(
    RealAlloc1D& x_knot,
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
      /* double& */ x,
      /* int& */ _ix0,
      /* bool* */ _strict,
      /* bool* */ _print_err,
      /* bool& */ _ok);
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
      /* double& */ x,
      /* int* */ _n,
      /* double& */ _y);
  return _y;
}
bool SimUtils::spline_akima(SplineProxyAlloc1D& spline) {
  // intent=inout allocatable type array
  bool _ok{};
  fortran_spline_akima(/* void* */ spline.get_fortran_ptr(), /* bool& */ _ok);
  return _ok;
}
void SimUtils::end_akima_spline_calc(
    SplineProxyAlloc1D& spline,
    int which_end) {
  // intent=inout allocatable type array
  fortran_end_akima_spline_calc(
      /* void* */ spline.get_fortran_ptr(), /* int& */ which_end);
}
SimUtils::ApfftCorr SimUtils::apfft_corr(
    RealAlloc1D& rdata_in,
    std::optional<FixedArray1D<Real, 2>> bounds,
    std::string window,
    std::optional<int> diag) {
  // intent=in allocatable general array
  double* _bounds = bounds.has_value() ? bounds.value().data() : nullptr;
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
      /* double* */ _bounds,
      /* const char* */ _window,
      /* double& */ _phase,
      /* double& */ _amp,
      /* double& */ _freq,
      /* int* */ _diag);
  return ApfftCorr{_phase, _amp, _freq};
}
void SimUtils::apfft(
    RealAlloc1D& rdata_in,
    FixedArray1D<Real, 2> bounds,
    std::string& window,
    double& phase,
    optional_ref<int> diag) {
  // intent=inout allocatable general array
  auto* _bounds = bounds.data(); // CppWrapperGeneralArgument
  auto _window = window.c_str(); // ptr, inout, required
  auto* _diag = diag.has_value() ? &diag->get() : nullptr; // inout, optional
  fortran_apfft(
      /* void* */ rdata_in.get_fortran_ptr(),
      /* double* */ _bounds,
      /* const char* */ _window,
      /* double& */ phase,
      /* int* */ _diag);
}
void SimUtils::apfft_ext(
    RealAlloc1D& rdata,
    FixedArray1D<Real, 2> bounds,
    std::string& window,
    double& phase,
    double& amp,
    double& freq,
    optional_ref<int> diag) {
  // intent=inout allocatable general array
  auto* _bounds = bounds.data(); // CppWrapperGeneralArgument
  auto _window = window.c_str(); // ptr, inout, required
  auto* _diag = diag.has_value() ? &diag->get() : nullptr; // inout, optional
  fortran_apfft_ext(
      /* void* */ rdata.get_fortran_ptr(),
      /* double* */ _bounds,
      /* const char* */ _window,
      /* double& */ phase,
      /* double& */ amp,
      /* double& */ freq,
      /* int* */ _diag);
}
void SimUtils::hanhan(int& N, RealAlloc1D& hh) {
  // intent=inout allocatable general array
  fortran_hanhan(/* int& */ N, /* void* */ hh.get_fortran_ptr());
}
void SimUtils::bit_set(int& word, int pos, bool set_to_1) {
  fortran_bit_set(/* int& */ word, /* int& */ pos, /* bool& */ set_to_1);
}