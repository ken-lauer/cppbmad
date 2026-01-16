#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.hpp"

using namespace Bmad;

namespace SimUtils {

// Skipped unusable routine all_pointer_to_string:
// - Untranslated type: all_pointer_struct (0D)
extern "C" void fortran_allocate_thread_states();
void allocate_thread_states();
extern "C" bool fortran_anomalous_moment_of(
    int &species /* 0D_NOT_integer in */,
    double &moment /* 0D_NOT_real out */
);
double anomalous_moment_of(int species);
extern "C" bool fortran_antiparticle(
    int &species /* 0D_NOT_integer in */,
    int &anti_species /* 0D_NOT_integer out */
);
int antiparticle(int species);
extern "C" void fortran_apfft(
    Bmad::array_descriptor_t &rdata_in /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &bounds /* 1D_NOT_real inout */,
    const char *window /* 0D_NOT_character in */,
    double &phase /* 0D_NOT_real in */,
    int *diag /* 0D_NOT_integer in */
);
void apfft(
    FArray1D<Real> &rdata_in,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    std::optional<int> diag = std::nullopt
);
extern "C" void fortran_apfft_corr(
    Bmad::array_descriptor_t &rdata_in /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &bounds /* 1D_NOT_real in */,
    const char *window /* 0D_NOT_character in */,
    double &phase /* 0D_NOT_real out */,
    double &amp /* 0D_NOT_real out */,
    double &freq /* 0D_NOT_real out */,
    int *diag /* 0D_NOT_integer in */
);
struct ApfftCorr {
  double phase;
  double amp;
  double freq;
};
SimUtils::ApfftCorr apfft_corr(
    FArray1D<Real> &rdata_in,
    std::optional<FixedArray1D<Real, 2>> bounds,
    std::string window,
    std::optional<int> diag = std::nullopt
);
extern "C" void fortran_apfft_ext(
    Bmad::array_descriptor_t &rdata /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &bounds /* 1D_NOT_real inout */,
    const char *window /* 0D_NOT_character in */,
    double &phase /* 0D_NOT_real in */,
    double &amp /* 0D_NOT_real in */,
    double &freq /* 0D_NOT_real in */,
    int *diag /* 0D_NOT_integer in */
);
void apfft_ext(
    FArray1D<Real> &rdata,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    double amp,
    double freq,
    std::optional<int> diag = std::nullopt
);
extern "C" bool fortran_asinc(
    double &x /* 0D_NOT_real in */,
    int *nd /* 0D_NOT_integer in */,
    double &y /* 0D_NOT_real in */
);
void asinc(double x, std::optional<int> nd, double y);
extern "C" bool fortran_assert_equal(
    Bmad::array_descriptor_t &int_arr /* 1D_NOT_integer in */,
    const char *err_str /* 0D_NOT_character in */,
    int &ival /* 0D_NOT_integer in */
);
void assert_equal(FArray1D<Int> &int_arr, std::string err_str, int ival);
extern "C" bool fortran_atomic_number(
    int &species /* 0D_NOT_integer in */,
    int &atomic_num /* 0D_NOT_integer out */
);
int atomic_number(int species);
extern "C" bool fortran_atomic_species_id(
    int &charge /* 0D_NOT_integer in */,
    bool &is_anti /* 0D_NOT_logical in */,
    int &atomic_num /* 0D_NOT_integer in */,
    int &n_nuc /* 0D_NOT_integer in */,
    int &species_id /* 0D_NOT_integer out */
);
int atomic_species_id(int charge, bool is_anti, int atomic_num, int n_nuc);
extern "C" bool fortran_axis_angle_to_quat(
    Bmad::array_descriptor_t &axis /* 1D_NOT_real in */,
    double &angle /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &quat /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> axis_angle_to_quat(FixedArray1D<Real, 3> axis, double angle);
extern "C" void fortran_axis_angle_to_w_mat(
    Bmad::array_descriptor_t &axis /* 1D_NOT_real in */,
    double &angle /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3> axis_angle_to_w_mat(FixedArray1D<Real, 3> axis, double angle);
extern "C" bool fortran_bicubic_cmplx_eval(
    double &x_norm /* 0D_NOT_real in */,
    double &y_norm /* 0D_NOT_real in */,
    void *bi_coef /* 0D_NOT_type in */,
    std::complex<double> &df_dx /* 0D_NOT_complex out */,
    std::complex<double> &df_dy /* 0D_NOT_complex out */,
    std::complex<double> &f_val /* 0D_NOT_complex out */
);
struct BicubicCmplxEval {
  std::complex<double> df_dx;
  std::complex<double> df_dy;
  std::complex<double> f_val;
};
SimUtils::BicubicCmplxEval
bicubic_cmplx_eval(double x_norm, double y_norm, BicubicCmplxCoefStruct &bi_coef);

// Skipped unusable routine bicubic_compute_cmplx_field_at_2d_box:
// - Untranslated type: cmplx_field_at_2d_box_struct (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine bicubic_compute_field_at_2d_box:
// - Untranslated type: field_at_2d_box_struct (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine bicubic_eval:
// - Untranslated type: bicubic_coef_struct (0D)

// Skipped unusable routine bicubic_interpolation_cmplx_coefs:
// - Untranslated type: cmplx_field_at_2d_box_struct (0D)

// Skipped unusable routine bicubic_interpolation_coefs:
// - Untranslated type: field_at_2d_box_struct (0D)
// - Untranslated type: bicubic_coef_struct (0D)

// Skipped unusable routine bin_2d:
// - Untranslated type: general_bin_struct (0D)

// Skipped unusable routine bin_data:
// - Untranslated type: bin_struct (0D)

// Skipped unusable routine bin_data_density:
// - Untranslated type: bin_struct (0D)

// Skipped unusable routine bin_data_density_2d:
// - Untranslated type: general_bin_struct (0D)
extern "C" bool fortran_bin_index(
    double &x /* 0D_NOT_real in */,
    double &bin1_x_min /* 0D_NOT_real in */,
    double &bin_delta /* 0D_NOT_real in */,
    int &ix_bin /* 0D_NOT_integer out */
);
int bin_index(double x, double bin1_x_min, double bin_delta);
extern "C" bool fortran_bin_x_center(
    int &ix_bin /* 0D_NOT_integer in */,
    double &bin1_x_min /* 0D_NOT_real in */,
    double &bin_delta /* 0D_NOT_real in */,
    double &x_center /* 0D_NOT_real out */
);
double bin_x_center(int ix_bin, double bin1_x_min, double bin_delta);
extern "C" void fortran_bit_set(
    int &word /* 0D_NOT_integer in */,
    int &pos /* 0D_NOT_integer in */,
    bool &set_to_1 /* 0D_NOT_logical in */
);
void bit_set(int word, int pos, bool set_to_1);

// Skipped unusable routine bracket_index:
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine bracket_index2:
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_bracket_index_for_spline(
    Bmad::array_descriptor_t &x_knot /* 1D_NOT_real in */,
    double &x /* 0D_NOT_real in */,
    int &ix0 /* 0D_NOT_integer out */,
    bool *strict /* 0D_NOT_logical in */,
    bool *print_err /* 0D_NOT_logical in */,
    bool &ok /* 0D_NOT_logical out */
);
struct BracketIndexForSpline {
  int ix0;
  bool ok;
};
SimUtils::BracketIndexForSpline bracket_index_for_spline(
    FArray1D<Real> &x_knot,
    double x,
    std::optional<bool> strict = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);

// Skipped unusable routine bracket_index_int:
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_calc_file_number(
    const char *file_name /* 0D_NOT_character in */,
    int &num_in /* 0D_NOT_integer in */,
    int &num_out /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void calc_file_number(std::string file_name, int num_in, int num_out, bool err_flag);
extern "C" void fortran_change_file_number(
    const char *file_name /* 0D_NOT_character in */,
    int &change /* 0D_NOT_integer in */
);
void change_file_number(std::string file_name, int change);
extern "C" bool fortran_charge_of(
    int &species /* 0D_NOT_integer in */,
    int *default_ /* 0D_NOT_integer in */,
    int &charge /* 0D_NOT_integer out */
);
int charge_of(int species, std::optional<int> default_ = std::nullopt);
extern "C" bool fortran_charge_to_mass_of(
    int &species /* 0D_NOT_integer in */,
    double &charge_mass_ratio /* 0D_NOT_real out */
);
double charge_to_mass_of(int species);
extern "C" bool fortran_coarse_frequency_estimate(
    Bmad::array_descriptor_t &data /* 1D_NOT_real in */,
    bool *error /* 0D_NOT_logical in */,
    double &frequency /* 0D_NOT_real out */
);
double coarse_frequency_estimate(FArray1D<Real> &data, std::optional<bool> error = std::nullopt);
extern "C" void fortran_complex_error_function(
    double &wr /* 0D_NOT_real in */,
    double &wi /* 0D_NOT_real in */,
    double &zr /* 0D_NOT_real in */,
    double &zi /* 0D_NOT_real in */
);
void complex_error_function(double wr, double wi, double zr, double zi);
extern "C" bool
fortran_cos_one(double &angle /* 0D_NOT_real in */, double &cos1 /* 0D_NOT_real in */);
void cos_one(double angle, double cos1);
extern "C" bool fortran_cosc(
    double &x /* 0D_NOT_real in */,
    int *nd /* 0D_NOT_integer in */,
    double &y /* 0D_NOT_real in */
);
void cosc(double x, std::optional<int> nd, double y);

// Skipped unusable routine count_at_index:
// - Untranslated type: bin_struct (0D)

// Skipped unusable routine covar_expand:
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine cplx_mat_inverse:
// - Variable inout sized array: 2D_NOT_complex
// - Variable inout sized array: 2D_NOT_complex

// Skipped unusable routine cplx_mat_make_unit:
// - Variable inout sized array: 2D_NOT_complex
extern "C" bool fortran_create_a_spline(
    Bmad::array_descriptor_t &r0 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &r1 /* 1D_NOT_real in */,
    double &slope0 /* 0D_NOT_real in */,
    double &slope1 /* 0D_NOT_real in */,
    void *spline /* 0D_NOT_type out */
);
SplineStruct create_a_spline(FArray1D<Real> &r0, FArray1D<Real> &r1, double slope0, double slope1);
extern "C" bool fortran_cross_product(
    Bmad::array_descriptor_t &a /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &b /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &c /* 1D_NOT_real inout */
);
void cross_product(FArray1D<Real> &a, FArray1D<Real> &b, FixedArray1D<Real, 3> c);

// Skipped unusable routine da2_div:
// - Variable in sized array: 2D_NOT_real
// - Variable in sized array: 2D_NOT_real
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine da2_evaluate:
// - Variable in sized array: 2D_NOT_real

// Skipped unusable routine da2_inverse:
// - Variable in sized array: 2D_NOT_real
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine da2_mult:
// - Variable in sized array: 2D_NOT_real
// - Variable in sized array: 2D_NOT_real
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_date_and_time_stamp(
    const char *string /* 0D_NOT_character in */,
    bool *numeric_month /* 0D_NOT_logical in */,
    bool *include_zone /* 0D_NOT_logical in */
);
void date_and_time_stamp(
    std::string string,
    std::optional<bool> numeric_month = std::nullopt,
    std::optional<bool> include_zone = std::nullopt
);
extern "C" void fortran_destfixedwindowls(int &id /* 0D_NOT_integer in */);
void destfixedwindowls(int id);
extern "C" void fortran_detab(const char *str /* 0D_NOT_character in */);
void detab(std::string str);

// Skipped unusable routine determinant:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_display_size_and_resolution(
    int &ix_screen /* 0D_NOT_integer in */,
    double &x_size /* 0D_NOT_real in */,
    double &y_size /* 0D_NOT_real in */,
    double &x_res /* 0D_NOT_real in */,
    double &y_res /* 0D_NOT_real in */
);
void display_size_and_resolution(
    int ix_screen,
    double x_size,
    double y_size,
    double x_res,
    double y_res
);
extern "C" bool fortran_dj_bessel(
    int &m /* 0D_NOT_integer in */,
    double &arg /* 0D_NOT_real in */,
    double &dj_bes /* 0D_NOT_real in */
);
void dj_bessel(int m, double arg, double dj_bes);
extern "C" bool fortran_djb_hash(
    const char *str /* 0D_NOT_character in */,
    int *old_hash /* 0D_NOT_integer in */,
    int &hash /* 0D_NOT_integer in */
);
void djb_hash(std::string str, std::optional<int> old_hash, int hash);
extern "C" bool fortran_djb_str_hash(
    const char *in_str /* 0D_NOT_character in */,
    const char *hash_str /* 0D_NOT_character in */
);
void djb_str_hash(std::string in_str, std::string hash_str);

// Skipped unusable routine doubleup_quotes:
// - No matching docstring

// Skipped unusable routine downcase:
// - No matching docstring
extern "C" void fortran_downcase_string(const char *string /* 0D_NOT_character in */);
void downcase_string(std::string string);

// Skipped unusable routine ed:
// - Routine in configuration skip list
extern "C" void fortran_end_akima_spline_calc(
    Bmad::array_descriptor_t &spline /* 1D_NOT_type inout */,
    int &which_end /* 0D_NOT_integer in */
);
void end_akima_spline_calc(SplineStructArray1D spline, int which_end);
extern "C" void fortran_err_exit(const char *err_str /* 0D_NOT_character in */);
void err_exit(std::optional<std::string> err_str = std::nullopt);
extern "C" bool
fortran_factorial(int &n /* 0D_NOT_integer in */, double &fact /* 0D_NOT_real in */);
void factorial(int n, double fact);
extern "C" void fortran_faddeeva_function(
    Bmad::array_descriptor_t &z /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &w /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &dw /* 2D_NOT_real inout */
);
void faddeeva_function(
    FixedArray1D<Real, 2> z,
    FixedArray1D<Real, 2> w,
    FixedArray2D<Real, 2, 2> dw
);
extern "C" void fortran_fft_1d(
    Bmad::array_descriptor_t &arr /* 1D_NOT_complex inout */,
    int &isign /* 0D_NOT_integer in */
);
void fft_1d(FArray1D<Complex> &arr, int isign);
extern "C" void fortran_file_directorizer(
    const char *in_file /* 0D_NOT_character in */,
    const char *out_file /* 0D_NOT_character in */,
    const char *directory /* 0D_NOT_character in */,
    bool &add_switch /* 0D_NOT_logical in */
);
void file_directorizer(
    std::string in_file,
    std::string out_file,
    std::string directory,
    bool add_switch
);
extern "C" void fortran_file_get(
    const char *string /* 0D_NOT_character in */,
    const char *dflt_file_name /* 0D_NOT_character in */,
    const char *file_name /* 0D_NOT_character in */
);
void file_get(std::string string, std::string dflt_file_name, std::string file_name);
extern "C" void fortran_file_get_open(
    const char *string /* 0D_NOT_character in */,
    const char *dflt_file_name /* 0D_NOT_character in */,
    const char *file_name /* 0D_NOT_character in */,
    int &file_unit /* 0D_NOT_integer in */,
    bool &readonly /* 0D_NOT_logical in */
);
void file_get_open(
    std::string string,
    std::string dflt_file_name,
    std::string file_name,
    int file_unit,
    bool readonly
);
extern "C" void fortran_file_suffixer(
    const char *in_file_name /* 0D_NOT_character in */,
    const char *out_file_name /* 0D_NOT_character in */,
    const char *suffix /* 0D_NOT_character in */,
    bool &add_switch /* 0D_NOT_logical in */
);
void file_suffixer(
    std::string in_file_name,
    std::string out_file_name,
    std::string suffix,
    bool add_switch
);
extern "C" bool fortran_find_location_int(
    Bmad::array_descriptor_t &arr /* 1D_NOT_integer inout */,
    int &value /* 0D_NOT_integer in */,
    int &ix_match /* 0D_NOT_integer in */
);
void find_location(FArray1D<Int> &arr, int value, int ix_match);
extern "C" bool fortran_find_location_logic(
    void *arr /* 1D_ALLOC_logical inout */,
    bool &value /* 0D_NOT_logical in */,
    int &ix_match /* 0D_NOT_integer in */
);
void find_location(BoolAlloc1D &arr, bool value, int ix_match);
extern "C" bool fortran_find_location_real(
    Bmad::array_descriptor_t &arr /* 1D_NOT_real in */,
    double &value /* 0D_NOT_real in */,
    int &ix_match /* 0D_NOT_integer in */
);
void find_location(FArray1D<Real> &arr, double value, int ix_match);

// Skipped unusable routine find_location_str:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_fine_frequency_estimate(
    Bmad::array_descriptor_t &data /* 1D_NOT_real in */,
    double &frequency /* 0D_NOT_real out */
);
double fine_frequency_estimate(FArray1D<Real> &data);
extern "C" bool fortran_fixedwindowls(
    double &ynew /* 0D_NOT_real in */,
    int &id /* 0D_NOT_integer in */,
    double &z /* 0D_NOT_real in */
);
void fixedwindowls(double ynew, int id, double z);
extern "C" void fortran_fourier_amplitude(
    Bmad::array_descriptor_t &data /* 1D_NOT_real in */,
    double &frequency /* 0D_NOT_real in */,
    double &cos_amp /* 0D_NOT_real out */,
    double &sin_amp /* 0D_NOT_real out */,
    double &dcos_amp /* 0D_NOT_real out */,
    double &dsin_amp /* 0D_NOT_real out */
);
struct FourierAmplitude {
  double cos_amp;
  double sin_amp;
  double dcos_amp;
  double dsin_amp;
};
SimUtils::FourierAmplitude fourier_amplitude(FArray1D<Real> &data, double frequency);
extern "C" bool fortran_gen_complete_elliptic(
    double &kc /* 0D_NOT_real in */,
    double &p /* 0D_NOT_real in */,
    double &c /* 0D_NOT_real in */,
    double &s /* 0D_NOT_real in */,
    double *err_tol /* 0D_NOT_real in */,
    double &value /* 0D_NOT_real in */
);
void gen_complete_elliptic(
    double kc,
    double p,
    double c,
    double s,
    std::optional<double> err_tol,
    double value
);

// Skipped unusable routine general_bin_count:
// - Untranslated type: general_bin_struct (0D)

// Skipped unusable routine general_bin_index:
// - Untranslated type: general_bin_struct (0D)

// Skipped unusable routine general_bin_index_in_bounds:
// - Untranslated type: general_bin_struct (0D)
extern "C" void fortran_get_file_number(
    const char *file_name /* 0D_NOT_character in */,
    const char *cnum_in /* 0D_NOT_character in */,
    int &num_out /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical in */
);
void get_file_number(std::string file_name, std::string cnum_in, int num_out, bool err_flag);
extern "C" void fortran_get_file_time_stamp(
    const char *file /* 0D_NOT_character in */,
    const char *time_stamp /* 0D_NOT_character in */
);
void get_file_time_stamp(std::string file, std::string time_stamp);

// Skipped unusable routine han:
// - Routine in configuration skip list
extern "C" void fortran_hanhan(
    int &N /* 0D_NOT_integer in */,
    Bmad::array_descriptor_t &hh /* 1D_NOT_real inout */
);
void hanhan(int N, FArray1D<Real> &hh);
extern "C" bool fortran_i_bessel(
    int &m /* 0D_NOT_integer in */,
    double &arg /* 0D_NOT_real in */,
    double &i_bes /* 0D_NOT_real in */
);
void i_bessel(int m, double arg, double i_bes);
extern "C" bool fortran_i_bessel_extended(
    int &m /* 0D_NOT_integer in */,
    double &arg /* 0D_NOT_real in */,
    std::complex<double> &i_bes /* 0D_NOT_complex in */
);
void i_bessel_extended(int m, double arg, std::complex<double> i_bes);
extern "C" void fortran_increment_file_number(
    const char *file_name /* 0D_NOT_character in */,
    int &digits /* 0D_NOT_integer in */,
    int &number /* 0D_NOT_integer in */,
    const char *cnumber /* 0D_NOT_character in */
);
void increment_file_number(std::string file_name, int digits, int number, std::string cnumber);
extern "C" bool fortran_index_nocase(
    const char *string1 /* 0D_NOT_character in */,
    const char *string2 /* 0D_NOT_character in */,
    int &indx /* 0D_NOT_integer in */
);
void index_nocase(std::string string1, std::string string2, int indx);
extern "C" bool fortran_initfixedwindowls(
    int &N /* 0D_NOT_integer in */,
    double &dt /* 0D_NOT_real in */,
    int &order /* 0D_NOT_integer in */,
    int &der /* 0D_NOT_integer in */,
    int &id /* 0D_NOT_integer out */
);
int initfixedwindowls(int N, double dt, int order, int der);

// Skipped unusable routine int_logic:
// - Routine in configuration skip list
extern "C" bool fortran_int_str(
    int &int_ /* 0D_NOT_integer in */,
    int *width /* 0D_NOT_integer in */,
    const char *str /* 0D_ALLOC_character in */
);
void int_str(int int_, std::optional<int> width, std::string str);
extern "C" bool fortran_interpolated_fft(
    Bmad::array_descriptor_t &cdata /* 1D_NOT_complex inout */,
    bool &calc_ok /* 0D_NOT_logical in */,
    int *opt_dump_spectrum /* 0D_NOT_integer in */,
    int *opt_dump_index /* 0D_NOT_integer in */,
    double &this_fft /* 0D_NOT_real in */
);
void interpolated_fft(
    FArray1D<Complex> &cdata,
    bool calc_ok,
    std::optional<int> opt_dump_spectrum,
    std::optional<int> opt_dump_index,
    double this_fft
);
extern "C" bool fortran_interpolated_fft_gsl(
    Bmad::array_descriptor_t &cdata /* 1D_NOT_complex inout */,
    bool &calc_ok /* 0D_NOT_logical in */,
    int *opt_dump_spectrum /* 0D_NOT_integer in */,
    int *opt_dump_index /* 0D_NOT_integer in */,
    double &this_fft /* 0D_NOT_real in */
);
void interpolated_fft_gsl(
    FArray1D<Complex> &cdata,
    bool calc_ok,
    std::optional<int> opt_dump_spectrum,
    std::optional<int> opt_dump_index,
    double this_fft
);

// Skipped unusable routine inverse:
// - Argument not defined: funct (have: [])
// - Argument not defined: y (have: [])
// - Argument not defined: x1 (have: [])
// - Argument not defined: x2 (have: [])
// - Argument not defined: tol (have: [])
// - Argument not defined: x (have: [])
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_is_alphabetic(
    const char *string /* 0D_NOT_character in */,
    const char *valid_chars /* 0D_NOT_character in */,
    bool &is_alpha /* 0D_NOT_logical in */
);
void is_alphabetic(std::string string, std::optional<std::string> valid_chars, bool is_alpha);
extern "C" bool fortran_is_decreasing_sequence(
    Bmad::array_descriptor_t &array /* 1D_NOT_real in */,
    bool *strict /* 0D_NOT_logical in */,
    bool &is_decreasing /* 0D_NOT_logical in */
);
void is_decreasing_sequence(FArray1D<Real> &array, std::optional<bool> strict, bool is_decreasing);
extern "C" bool
fortran_is_false(double &param /* 0D_NOT_real in */, bool &this_false /* 0D_NOT_logical out */);
bool is_false(double param);
extern "C" bool fortran_is_increasing_sequence(
    Bmad::array_descriptor_t &array /* 1D_NOT_real in */,
    bool *strict /* 0D_NOT_logical in */,
    bool &is_increasing /* 0D_NOT_logical in */
);
void is_increasing_sequence(FArray1D<Real> &array, std::optional<bool> strict, bool is_increasing);
extern "C" bool fortran_is_integer(
    const char *string /* 0D_NOT_character in */,
    int *int_ /* 0D_NOT_integer in */,
    const char *delims /* 0D_NOT_character in */,
    int *ix_word /* 0D_NOT_integer in */,
    bool &valid /* 0D_NOT_logical in */
);
void is_integer(
    std::string string,
    std::optional<int> int_,
    std::optional<std::string> delims,
    std::optional<int> ix_word,
    bool valid
);
extern "C" bool fortran_is_logical(
    const char *string /* 0D_NOT_character in */,
    bool *ignore /* 0D_NOT_logical in */,
    bool &valid /* 0D_NOT_logical in */
);
void is_logical(std::string string, std::optional<bool> ignore, bool valid);
extern "C" bool fortran_is_real(
    const char *string /* 0D_NOT_character in */,
    bool *ignore /* 0D_NOT_logical in */,
    double *real_num /* 0D_NOT_real in */,
    bool &valid /* 0D_NOT_logical in */
);
void is_real(
    std::string string,
    std::optional<bool> ignore,
    std::optional<double> real_num,
    bool valid
);
extern "C" bool fortran_is_subatomic_species(
    int &species /* 0D_NOT_integer in */,
    bool &is_subatomic /* 0D_NOT_logical out */
);
bool is_subatomic_species(int species);
extern "C" bool
fortran_is_true(double &param /* 0D_NOT_real in */, bool &this_true /* 0D_NOT_logical out */);
bool is_true(double param);

// Skipped unusable routine isatty:
// - Routine in configuration skip list
extern "C" bool fortran_j_bessel(
    int &m /* 0D_NOT_integer in */,
    double &arg /* 0D_NOT_real in */,
    double &j_bes /* 0D_NOT_real in */
);
void j_bessel(int m, double arg, double j_bes);
extern "C" void fortran_linear_fit(
    Bmad::array_descriptor_t &x /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &y /* 1D_NOT_real inout */,
    int &n_data /* 0D_NOT_integer in */,
    double &a /* 0D_NOT_real in */,
    double &b /* 0D_NOT_real in */,
    double &sig_a /* 0D_NOT_real in */,
    double &sig_b /* 0D_NOT_real in */
);
void linear_fit(
    FArray1D<Real> &x,
    FArray1D<Real> &y,
    int n_data,
    double a,
    double b,
    double sig_a,
    double sig_b
);
extern "C" void fortran_linear_fit_2d(
    Bmad::array_descriptor_t &x /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &z /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &coef /* 1D_NOT_real out */
);
FixedArray1D<Real, 3> linear_fit_2d(FArray1D<Real> &x, FArray1D<Real> &y, FArray1D<Real> &z);

// Skipped unusable routine location_decode:
// - Translated arg count mismatch (unsupported?)
extern "C" bool
fortran_logic_str(bool &logic /* 0D_NOT_logical in */, const char *str /* 0D_NOT_character in */);
void logic_str(bool logic, std::string str);
extern "C" bool fortran_lunget(int &func_retval__ /* 0D_NOT_integer out */);
int lunget();
extern "C" void fortran_make_legal_comment(
    const char *comment_in /* 0D_NOT_character in */,
    const char *comment_out /* 0D_NOT_character in */
);
void make_legal_comment(std::string comment_in, std::string comment_out);
extern "C" bool
fortran_mass_of(int &species /* 0D_NOT_integer in */, double &mass /* 0D_NOT_real out */);
double mass_of(int species);

// Skipped unusable routine mat_eigen:
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_complex

// Skipped unusable routine mat_inverse:
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine mat_make_unit:
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine mat_pseudoinverse:
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine mat_rotation:
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine mat_scale_p0:
// - Variable in sized array: 2D_NOT_real
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine mat_symp_conj:
// - Variable inout sized array: 2D_NOT_real
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine mat_symp_conj_i:
// - Variable inout sized array: 2D_NOT_complex
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine mat_symp_error:
// - Variable in sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine mat_symplectify:
// - Variable in sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine mat_type:
// - Variable inout sized array: 2D_NOT_real
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_match_reg(
    const char *str /* 0D_NOT_character in */,
    const char *pat /* 0D_NOT_character in */,
    bool &is_match /* 0D_NOT_logical in */
);
void match_reg(std::string str, std::string pat, bool is_match);
extern "C" bool fortran_match_wild(
    const char *string /* 0D_NOT_character in */,
    const char *template_ /* 0D_NOT_character in */,
    bool &is_match /* 0D_NOT_logical in */
);
void match_wild(std::string string, std::string template_, bool is_match);

// Skipped unusable routine match_word:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine max_nonzero:
// - Module name unset
extern "C" bool fortran_maximize_projection(
    double &seed /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &cdata /* 1D_NOT_complex inout */,
    double &func_retval__ /* 0D_NOT_real in */
);
void maximize_projection(double seed, FArray1D<Complex> &cdata, double func_retval__);
extern "C" void fortran_milli_sleep(int &milli_sec /* 0D_NOT_integer in */);
void milli_sleep(int milli_sec);

// Skipped unusable routine molecular_components:
// - Untranslated type: molecular_component_struct (1D)
extern "C" bool
fortran_n_bins_automatic(int &n_data /* 0D_NOT_integer in */, int &n /* 0D_NOT_integer in */);
void n_bins_automatic(int n_data, int n);
extern "C" bool fortran_n_choose_k(
    int &n /* 0D_NOT_integer in */,
    int &k /* 0D_NOT_integer in */,
    double &nck /* 0D_NOT_real in */
);
void n_choose_k(int n, int k, double nck);
extern "C" void fortran_n_spline_create(
    Bmad::array_descriptor_t &deriv0 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &deriv1 /* 1D_NOT_real in */,
    double &x1 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &n_spline /* 1D_NOT_real inout */
);
void n_spline_create(
    FArray1D<Real> &deriv0,
    FArray1D<Real> &deriv1,
    double x1,
    FArray1D<Real> &n_spline
);
extern "C" void fortran_naff(
    Bmad::array_descriptor_t &cdata /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &freqs /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &amps /* 1D_NOT_complex inout */,
    int *opt_dump_spectra /* 0D_NOT_integer in */,
    bool *opt_zero_first /* 0D_NOT_logical in */
);
void naff(
    FArray1D<Complex> &cdata,
    FArray1D<Real> &freqs,
    FArray1D<Complex> &amps,
    std::optional<int> opt_dump_spectra = std::nullopt,
    std::optional<bool> opt_zero_first = std::nullopt
);
extern "C" void fortran_nametable_add(
    void *nametable /* 0D_NOT_type inout */,
    const char *name /* 0D_NOT_character in */,
    int &ix_name /* 0D_NOT_integer in */
);
void nametable_add(NametableStruct &nametable, std::string name, int ix_name);
extern "C" bool fortran_nametable_bracket_indexx(
    void *nametable /* 0D_NOT_type inout */,
    const char *name /* 0D_NOT_character in */,
    int *n_match /* 0D_NOT_integer in */,
    int &ix_max /* 0D_NOT_integer in */
);
void nametable_bracket_indexx(
    NametableStruct &nametable,
    std::string name,
    std::optional<int> n_match,
    int ix_max
);
extern "C" void fortran_nametable_change1(
    void *nametable /* 0D_NOT_type inout */,
    const char *name /* 0D_NOT_character in */,
    int &ix_name /* 0D_NOT_integer in */
);
void nametable_change1(NametableStruct &nametable, std::string name, int ix_name);
extern "C" void fortran_nametable_init(
    void *nametable /* 0D_NOT_type inout */,
    int *n_min /* 0D_NOT_integer in */,
    int *n_max /* 0D_NOT_integer in */
);
void nametable_init(
    NametableStruct &nametable,
    std::optional<int> n_min = std::nullopt,
    std::optional<int> n_max = std::nullopt
);
extern "C" void fortran_nametable_remove(
    void *nametable /* 0D_NOT_type inout */,
    int &ix_name /* 0D_NOT_integer in */
);
void nametable_remove(NametableStruct &nametable, int ix_name);

// Skipped unusable routine node_put:
// - Routine in configuration skip list
extern "C" bool fortran_omega_to_quat(
    Bmad::array_descriptor_t &omega /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &quat /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> omega_to_quat(FixedArray1D<Real, 3> omega);
extern "C" bool fortran_openpmd_species_name(
    int &species /* 0D_NOT_integer in */,
    const char *pmd_name /* 0D_NOT_character out */
);
std::string openpmd_species_name(int species);
extern "C" bool
fortran_ordinal_str(int &n /* 0D_NOT_integer in */, const char *str /* 0D_ALLOC_character in */);
void ordinal_str(int n, std::string str);

// Skipped unusable routine outer_product:
// - Module name unset
extern "C" void fortran_parse_fortran_format(
    const char *format_str /* 0D_NOT_character in */,
    int &n_repeat /* 0D_NOT_integer in */,
    int &power /* 0D_NOT_integer in */,
    const char *descrip /* 0D_NOT_character in */,
    int &width /* 0D_NOT_integer in */,
    int &digits /* 0D_NOT_integer in */
);
void parse_fortran_format(
    std::string format_str,
    int n_repeat,
    int power,
    std::string descrip,
    int width,
    int digits
);

// Skipped unusable routine pointer_to_locations:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" bool fortran_pointer_to_ran_state(
    void *ran_state /* 0D_NOT_type in */,
    int *ix_thread /* 0D_NOT_integer in */,
    void *ran_state_ptr /* 0D_PTR_type out */
);
RandomStateStruct pointer_to_ran_state(
    optional_ref<RandomStateStruct> ran_state = std::nullopt,
    std::optional<int> ix_thread = std::nullopt
);
extern "C" bool fortran_poly_eval(
    Bmad::array_descriptor_t &poly /* 1D_NOT_real in */,
    double &x /* 0D_NOT_real in */,
    bool *diff_coef /* 0D_NOT_logical in */,
    double &y /* 0D_NOT_real in */
);
void poly_eval(FArray1D<Real> &poly, double x, std::optional<bool> diff_coef, double y);
extern "C" bool
fortran_probability_funct(double &x /* 0D_NOT_real in */, double &prob /* 0D_NOT_real in */);
void probability_funct(double x, double prob);
extern "C" bool fortran_projdd(
    Bmad::array_descriptor_t &a /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &b /* 1D_NOT_complex inout */,
    std::complex<double> &func_retval__ /* 0D_NOT_complex in */
);
void projdd(FArray1D<Complex> &a, FArray1D<Complex> &b, std::complex<double> func_retval__);
extern "C" bool fortran_quadratic_roots(
    Bmad::array_descriptor_t &coefs /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &root /* 1D_NOT_complex inout */
);
void quadratic_roots(FixedArray1D<Real, 3> coefs, FixedArray1D<Complex, 2> root);
extern "C" bool fortran_quat_conj_complex(
    Bmad::array_descriptor_t &q_in /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &q_out /* 1D_NOT_complex out */
);
FixedArray1D<Complex, 4> quat_conj(FixedArray1D<Complex, 4> q_in);
extern "C" bool fortran_quat_conj_real(
    Bmad::array_descriptor_t &q_in /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &q_out /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> quat_conj(FixedArray1D<Real, 4> q_in);
extern "C" bool fortran_quat_inverse(
    Bmad::array_descriptor_t &q_in /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &q_out /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> quat_inverse(FixedArray1D<Real, 4> q_in);
extern "C" bool fortran_quat_mul_complex(
    Bmad::array_descriptor_t &q1 /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &q2 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q3 /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &q4 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q5 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q6 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q7 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q8 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q9 /* 1D_NOT_complex inout */,
    Bmad::array_descriptor_t &q_out /* 1D_NOT_complex out */
);
FixedArray1D<Complex, 4> quat_mul(
    FixedArray1D<Complex, 4> q1,
    FixedArray1D<Complex, 4> q2,
    std::optional<FixedArray1D<Complex, 4>> q3 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q4 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q5 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q6 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q7 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q8 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q9 = std::nullopt
);
extern "C" bool fortran_quat_mul_real(
    Bmad::array_descriptor_t &q1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &q2 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q3 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &q4 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q5 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q6 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q7 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q8 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q9 /* 1D_NOT_real inout */,
    Bmad::array_descriptor_t &q_out /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> quat_mul(
    FixedArray1D<Real, 4> q1,
    FixedArray1D<Real, 4> q2,
    std::optional<FixedArray1D<Real, 4>> q3 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q4 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q5 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q6 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q7 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q8 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q9 = std::nullopt
);
extern "C" bool fortran_quat_rotate_complex(
    Bmad::array_descriptor_t &quat /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &vec_in /* 1D_NOT_complex in */,
    Bmad::array_descriptor_t &vec_out /* 1D_NOT_complex out */
);
FixedArray1D<Complex, 3>
quat_rotate(FixedArray1D<Complex, 4> quat, FixedArray1D<Complex, 3> vec_in);
extern "C" bool fortran_quat_rotate_real(
    Bmad::array_descriptor_t &quat /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &vec_in /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &vec_out /* 1D_NOT_real out */
);
FixedArray1D<Real, 3> quat_rotate(FixedArray1D<Real, 4> quat, FixedArray1D<Real, 3> vec_in);
extern "C" void fortran_quat_to_axis_angle(
    Bmad::array_descriptor_t &quat /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &axis /* 1D_NOT_real out */,
    double &angle /* 0D_NOT_real out */
);
struct QuatToAxisAngle {
  FixedArray1D<Real, 3> axis;
  double angle;
};
SimUtils::QuatToAxisAngle quat_to_axis_angle(FixedArray1D<Real, 4> quat);
extern "C" bool fortran_quat_to_omega(
    Bmad::array_descriptor_t &quat /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &omega /* 1D_NOT_real out */
);
FixedArray1D<Real, 3> quat_to_omega(FixedArray1D<Real, 4> quat);
extern "C" bool fortran_quat_to_w_mat(
    Bmad::array_descriptor_t &quat /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real out */
);
FixedArray2D<Real, 3, 3> quat_to_w_mat(FixedArray1D<Real, 4> quat);
extern "C" void fortran_query_string(
    const char *query_str /* 0D_NOT_character in */,
    bool &upcase /* 0D_NOT_logical in */,
    const char *return_str /* 0D_NOT_character in */,
    int &ix /* 0D_NOT_integer in */,
    int &ios /* 0D_NOT_integer in */
);
void query_string(std::string query_str, bool upcase, std::string return_str, int ix, int ios);
extern "C" bool fortran_quote(
    const char *str /* 0D_NOT_character in */,
    const char *q_str /* 0D_ALLOC_character in */
);
void quote(std::string str, std::string q_str);

// Skipped unusable routine quoten:
// - Variable-sized inout character array: 1D_NOT_character
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_ran_default_state(
    void *set_state /* 0D_NOT_type in */,
    void *get_state /* 0D_NOT_type out */
);
RandomStateStruct ran_default_state(optional_ref<RandomStateStruct> set_state = std::nullopt);
extern "C" void fortran_ran_engine(
    const char *set /* 0D_NOT_character in */,
    const char *get /* 0D_NOT_character in */,
    void *ran_state /* 0D_NOT_type in */
);
void ran_engine(
    std::optional<std::string> set = std::nullopt,
    std::optional<std::string> get = std::nullopt,
    optional_ref<RandomStateStruct> ran_state = std::nullopt
);
extern "C" void fortran_ran_gauss_converter(
    const char *set /* 0D_NOT_character in */,
    double *set_sigma_cut /* 0D_NOT_real in */,
    const char *get /* 0D_NOT_character out */,
    double &get_sigma_cut /* 0D_NOT_real out */,
    void *ran_state /* 0D_NOT_type in */
);
struct RanGaussConverter {
  std::string get;
  double get_sigma_cut;
};
SimUtils::RanGaussConverter ran_gauss_converter(
    std::optional<std::string> set = std::nullopt,
    std::optional<double> set_sigma_cut = std::nullopt,
    optional_ref<RandomStateStruct> ran_state = std::nullopt
);
extern "C" void fortran_ran_gauss_scalar(
    double &harvest /* 0D_NOT_real out */,
    void *ran_state /* 0D_NOT_type in */,
    double *sigma_cut /* 0D_NOT_real in */,
    int *index_quasi /* 0D_NOT_integer in */
);
double ran_gauss_scalar(
    optional_ref<RandomStateStruct> ran_state = std::nullopt,
    std::optional<double> sigma_cut = std::nullopt,
    std::optional<int> index_quasi = std::nullopt
);
extern "C" void fortran_ran_gauss_vector(
    Bmad::array_descriptor_t &harvest /* 1D_NOT_real inout */,
    void *ran_state /* 0D_NOT_type in */,
    double *sigma_cut /* 0D_NOT_real in */
);
void ran_gauss_vector(
    FArray1D<Real> &harvest,
    optional_ref<RandomStateStruct> ran_state = std::nullopt,
    std::optional<double> sigma_cut = std::nullopt
);
extern "C" void fortran_ran_seed_get(int &seed /* 0D_NOT_integer out */);
int ran_seed_get();
extern "C" void
fortran_ran_seed_put(int &seed /* 0D_NOT_integer in */, int *mpi_offset /* 0D_NOT_integer in */);
void ran_seed_put(int seed, std::optional<int> mpi_offset = std::nullopt);
extern "C" void fortran_ran_uniform_scalar(
    double &harvest /* 0D_NOT_real out */,
    void *ran_state /* 0D_NOT_type in */,
    int *index_quasi /* 0D_NOT_integer in */
);
double ran_uniform(
    optional_ref<RandomStateStruct> ran_state = std::nullopt,
    std::optional<int> index_quasi = std::nullopt
);
extern "C" void fortran_ran_uniform_vector(
    Bmad::array_descriptor_t &harvest /* 1D_NOT_real inout */,
    void *ran_state /* 0D_NOT_type in */
);
void ran_uniform(FArray1D<Real> &harvest, optional_ref<RandomStateStruct> ran_state = std::nullopt);
extern "C" bool fortran_real_num_fortran_format(
    double &number /* 0D_NOT_real in */,
    int &width /* 0D_NOT_integer in */,
    int *n_blanks /* 0D_NOT_integer in */,
    const char *fmt_str /* 0D_NOT_character in */
);
void real_num_fortran_format(
    double number,
    int width,
    std::optional<int> n_blanks,
    std::string fmt_str
);
extern "C" bool fortran_real_path(
    const char *path_in /* 0D_NOT_character in */,
    const char *path_out /* 0D_NOT_character in */,
    bool &is_ok /* 0D_NOT_logical in */
);
void real_path(std::string path_in, std::string path_out, bool is_ok);
extern "C" bool fortran_real_str(
    double &r_num /* 0D_NOT_real in */,
    int *n_signif /* 0D_NOT_integer in */,
    int *n_decimal /* 0D_NOT_integer in */,
    const char *str /* 0D_ALLOC_character in */
);
void real_str(
    double r_num,
    std::optional<int> n_signif,
    std::optional<int> n_decimal,
    std::string str
);
extern "C" bool fortran_real_to_string(
    double &real_num /* 0D_NOT_real in */,
    int &width /* 0D_NOT_integer in */,
    int *n_signif /* 0D_NOT_integer in */,
    int *n_decimal /* 0D_NOT_integer in */,
    const char *str /* 0D_NOT_character in */
);
void real_to_string(
    double real_num,
    int width,
    std::optional<int> n_signif,
    std::optional<int> n_decimal,
    std::string str
);
extern "C" void fortran_reallocate_spline(
    void *spline /* 1D_ALLOC_type inout */,
    int &n /* 0D_NOT_integer in */,
    int *n_min /* 0D_NOT_integer in */,
    bool *exact /* 0D_NOT_logical in */
);
void reallocate_spline(
    SplineStructAlloc1D spline,
    int n,
    std::optional<int> n_min = std::nullopt,
    std::optional<bool> exact = std::nullopt
);

// Skipped unusable routine reals_to_string:
// - No matching docstring

// Skipped unusable routine reals_to_table_row:
// - No matching docstring
extern "C" bool fortran_rms_value(
    Bmad::array_descriptor_t &val_arr /* 1D_NOT_real in */,
    void *good_val /* 1D_ALLOC_logical in */,
    double &ave_val /* 0D_NOT_real out */,
    double &rms_val /* 0D_NOT_real in */
);
double rms_value(FArray1D<Real> &val_arr, optional_ref<BoolAlloc1D> good_val, double rms_val);
extern "C" bool fortran_rot_2d(
    Bmad::array_descriptor_t &vec_in /* 1D_NOT_real in */,
    double &angle /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &vec_out /* 1D_NOT_real inout */
);
void rot_2d(FixedArray1D<Real, 2> vec_in, double angle, FixedArray1D<Real, 2> vec_out);

// Skipped unusable routine rotate_mat:
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_rotate_vec(
    Bmad::array_descriptor_t &vec /* 1D_NOT_real inout */,
    int &axis /* 0D_NOT_integer in */,
    double &angle /* 0D_NOT_real in */
);
void rotate_vec(FArray1D<Real> &vec, int axis, double angle);
extern "C" bool fortran_rotate_vec_given_axis_angle(
    Bmad::array_descriptor_t &vec_in /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &axis /* 1D_NOT_real in */,
    double &angle /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &vec_out /* 1D_NOT_real out */
);
FixedArray1D<Real, 3>
rotate_vec_given_axis_angle(FixedArray1D<Real, 3> vec_in, FArray1D<Real> &axis, double angle);
extern "C" bool
fortran_rp8(int &int_in /* 0D_NOT_integer in */, double &re_out /* 0D_NOT_real out */);
double rp8(int int_in);
extern "C" void fortran_run_timer(
    const char *command /* 0D_NOT_character in */,
    double *time /* 0D_NOT_real in */,
    double *time0 /* 0D_NOT_real in */
);
void run_timer(
    std::string command,
    std::optional<double> time = std::nullopt,
    std::optional<double> time0 = std::nullopt
);

// Skipped unusable routine set_all_ptr:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_set_parameter_int(
    int &param_val /* 0D_NOT_integer in */,
    int &set_val /* 0D_NOT_integer in */,
    int &save_val /* 0D_NOT_integer in */
);
void set_parameter(int param_val, int set_val, int save_val);
extern "C" bool fortran_set_parameter_logic(
    bool &param_val /* 0D_NOT_logical in */,
    bool &set_val /* 0D_NOT_logical in */,
    bool &save_val /* 0D_NOT_logical in */
);
void set_parameter(bool param_val, bool set_val, bool save_val);
extern "C" bool fortran_set_parameter_real(
    double &param_val /* 0D_NOT_real in */,
    double &set_val /* 0D_NOT_real in */,
    double &save_val /* 0D_NOT_real in */
);
void set_parameter(double param_val, double set_val, double save_val);
extern "C" bool fortran_set_species_charge(
    int &species_in /* 0D_NOT_integer in */,
    int &charge /* 0D_NOT_integer in */,
    int &species_charged /* 0D_NOT_integer out */
);
int set_species_charge(int species_in, int charge);
extern "C" bool fortran_sinc(
    double &x /* 0D_NOT_real in */,
    int *nd /* 0D_NOT_integer in */,
    double &y /* 0D_NOT_real in */
);
void sinc(double x, std::optional<int> nd, double y);
extern "C" bool fortran_sincc(
    double &x /* 0D_NOT_real in */,
    int *nd /* 0D_NOT_integer in */,
    double &y /* 0D_NOT_real in */
);
void sincc(double x, std::optional<int> nd, double y);
extern "C" bool fortran_sinhx_x(
    double &x /* 0D_NOT_real in */,
    int *nd /* 0D_NOT_integer in */,
    double &y /* 0D_NOT_real in */
);
void sinhx_x(double x, std::optional<int> nd, double y);
extern "C" void
fortran_skip_header(int &ix_unit /* 0D_NOT_integer in */, bool &error_flag /* 0D_NOT_logical in */);
void skip_header(int ix_unit, bool error_flag);
extern "C" bool fortran_species_id(
    const char *name /* 0D_NOT_character in */,
    int *default_ /* 0D_NOT_integer in */,
    bool *print_err /* 0D_NOT_logical in */,
    int &species /* 0D_NOT_integer out */
);
int species_id(
    std::string name,
    std::optional<int> default_ = std::nullopt,
    std::optional<bool> print_err = std::nullopt
);
extern "C" bool fortran_species_id_from_openpmd(
    const char *pmd_name /* 0D_NOT_character in */,
    int &charge /* 0D_NOT_integer in */,
    int &species /* 0D_NOT_integer out */
);
int species_id_from_openpmd(std::string pmd_name, int charge);
extern "C" bool fortran_species_name(
    int &species /* 0D_NOT_integer in */,
    const char *name /* 0D_NOT_character out */
);
std::string species_name(int species);
extern "C" bool fortran_species_of(
    double &mass /* 0D_NOT_real in */,
    int &charge /* 0D_NOT_integer in */,
    int &species /* 0D_NOT_integer out */
);
int species_of(double mass, int charge);
extern "C" bool fortran_spin_of(
    int &species /* 0D_NOT_integer in */,
    double *non_subatomic_default /* 0D_NOT_real in */,
    double &spin /* 0D_NOT_real out */
);
double spin_of(int species, std::optional<double> non_subatomic_default = std::nullopt);
extern "C" bool fortran_spline1(
    void *a_spline /* 0D_NOT_type in */,
    double &x /* 0D_NOT_real in */,
    int *n /* 0D_NOT_integer in */,
    double &y /* 0D_NOT_real out */
);
double spline1(SplineStruct &a_spline, double x, std::optional<int> n = std::nullopt);
extern "C" void fortran_spline_akima(
    Bmad::array_descriptor_t &spline /* 1D_NOT_type inout */,
    bool &ok /* 0D_NOT_logical out */
);
bool spline_akima(SplineStructArray1D spline);
extern "C" void fortran_spline_akima_interpolate(
    Bmad::array_descriptor_t &x_knot /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y_knot /* 1D_NOT_real in */,
    double &x /* 0D_NOT_real in */,
    bool &ok /* 0D_NOT_logical out */,
    double &y /* 0D_NOT_real out */,
    double &dy /* 0D_NOT_real out */
);
struct SplineAkimaInterpolate {
  bool ok;
  double y;
  double dy;
};
SimUtils::SplineAkimaInterpolate
spline_akima_interpolate(FArray1D<Real> &x_knot, FArray1D<Real> &y_knot, double x);
extern "C" void fortran_spline_evaluate(
    Bmad::array_descriptor_t &spline /* 1D_NOT_type in */,
    double &x /* 0D_NOT_real in */,
    bool &ok /* 0D_NOT_logical out */,
    double &y /* 0D_NOT_real out */,
    double &dy /* 0D_NOT_real out */
);
struct SplineEvaluate {
  bool ok;
  double y;
  double dy;
};
SimUtils::SplineEvaluate spline_evaluate(SplineStructArray1D spline, double x);
extern "C" bool fortran_sqrt_alpha(
    double &alpha /* 0D_NOT_real in */,
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real in */
);
void sqrt_alpha(double alpha, double x, double y);
extern "C" bool fortran_sqrt_one(
    double &x /* 0D_NOT_real in */,
    int *nd /* 0D_NOT_integer in */,
    double &ds1 /* 0D_NOT_real in */
);
void sqrt_one(double x, std::optional<int> nd, double ds1);
extern "C" bool fortran_str_count(
    const char *str /* 0D_NOT_character in */,
    const char *match /* 0D_NOT_character in */,
    int &num /* 0D_NOT_integer in */
);
void str_count(std::string str, std::string match, int num);
extern "C" void fortran_str_downcase(
    const char *dst /* 0D_NOT_character out */,
    const char *src /* 0D_NOT_character in */
);
std::string str_downcase(std::string src);
extern "C" bool fortran_str_first_in_set(
    const char *line /* 0D_NOT_character in */,
    const char *set /* 0D_NOT_character in */,
    bool *ignore_clauses /* 0D_NOT_logical in */,
    int &ix_match /* 0D_NOT_integer in */
);
void str_first_in_set(
    std::string line,
    std::string set,
    std::optional<bool> ignore_clauses,
    int ix_match
);
extern "C" bool fortran_str_first_not_in_set(
    const char *line /* 0D_NOT_character in */,
    const char *set /* 0D_NOT_character in */,
    int &ix_match /* 0D_NOT_integer in */
);
void str_first_not_in_set(std::string line, std::string set, int ix_match);
extern "C" bool fortran_str_last_in_set(
    const char *line /* 0D_NOT_character in */,
    const char *set /* 0D_NOT_character in */,
    int &ix_match /* 0D_NOT_integer in */
);
void str_last_in_set(std::string line, std::string set, int ix_match);
extern "C" bool fortran_str_last_not_in_set(
    const char *line /* 0D_NOT_character in */,
    const char *set /* 0D_NOT_character in */,
    int &ix_match /* 0D_NOT_integer in */
);
void str_last_not_in_set(std::string line, std::string set, int ix_match);
extern "C" bool fortran_str_match_wild(
    const char *str /* 0D_NOT_character in */,
    const char *pat /* 0D_NOT_character in */,
    bool &a_match /* 0D_NOT_logical in */
);
void str_match_wild(std::string str, std::string pat, bool a_match);

// Skipped unusable routine str_set:
// - Routine in configuration skip list
extern "C" void fortran_str_substitute(
    const char *string /* 0D_NOT_character in */,
    const char *str_match /* 0D_NOT_character in */,
    const char *str_replace /* 0D_NOT_character in */,
    bool *do_trim /* 0D_NOT_logical in */,
    bool *ignore_escaped /* 0D_NOT_logical in */
);
void str_substitute(
    std::string string,
    std::optional<std::string> str_match = std::nullopt,
    std::optional<std::string> str_replace = std::nullopt,
    std::optional<bool> do_trim = std::nullopt,
    std::optional<bool> ignore_escaped = std::nullopt
);
extern "C" void fortran_str_upcase(
    const char *dst /* 0D_NOT_character out */,
    const char *src /* 0D_NOT_character in */
);
std::string str_upcase(std::string src);
extern "C" bool fortran_string_to_int(
    const char *line /* 0D_NOT_character in */,
    int &default_ /* 0D_NOT_integer in */,
    bool &err_flag /* 0D_NOT_logical in */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    int &value /* 0D_NOT_integer in */
);
void string_to_int(
    std::string line,
    int default_,
    bool err_flag,
    std::optional<bool> err_print_flag,
    int value
);
extern "C" bool fortran_string_to_real(
    const char *line /* 0D_NOT_character in */,
    double &default_ /* 0D_NOT_real in */,
    bool &err_flag /* 0D_NOT_logical in */,
    bool *err_print_flag /* 0D_NOT_logical in */,
    double &value /* 0D_NOT_real in */
);
void string_to_real(
    std::string line,
    double default_,
    bool err_flag,
    std::optional<bool> err_print_flag,
    double value
);
extern "C" void fortran_string_trim(
    const char *in_string /* 0D_NOT_character in */,
    const char *out_string /* 0D_NOT_character in */,
    int &word_len /* 0D_NOT_integer in */
);
void string_trim(std::string in_string, std::string out_string, int word_len);
extern "C" void fortran_string_trim2(
    const char *in_str /* 0D_NOT_character in */,
    const char *delimitors /* 0D_NOT_character in */,
    const char *out_str /* 0D_NOT_character in */,
    int &ix_word /* 0D_NOT_integer in */,
    const char *delim /* 0D_NOT_character in */,
    int &ix_next /* 0D_NOT_integer in */
);
void string_trim2(
    std::string in_str,
    std::string delimitors,
    std::string out_str,
    int ix_word,
    std::string delim,
    int ix_next
);

// Skipped unusable routine substr:
// - Routine in configuration skip list
extern "C" void fortran_super_bicubic_coef(
    Bmad::array_descriptor_t &y /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y2 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y12 /* 1D_NOT_real in */,
    double &d1 /* 0D_NOT_real in */,
    double &d2 /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &c /* 2D_NOT_real out */
);
FixedArray2D<Real, 4, 4> super_bicubic_coef(
    FixedArray1D<Real, 4> y,
    FixedArray1D<Real, 4> y1,
    FixedArray1D<Real, 4> y2,
    FixedArray1D<Real, 4> y12,
    double d1,
    double d2
);
extern "C" void fortran_super_bicubic_interpolation(
    Bmad::array_descriptor_t &y /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y1 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y2 /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &y12 /* 1D_NOT_real in */,
    double &x1l /* 0D_NOT_real in */,
    double &x1u /* 0D_NOT_real in */,
    double &x2l /* 0D_NOT_real in */,
    double &x2u /* 0D_NOT_real in */,
    double &x1 /* 0D_NOT_real in */,
    double &x2 /* 0D_NOT_real in */,
    double &ansy /* 0D_NOT_real out */,
    double &ansy1 /* 0D_NOT_real out */,
    double &ansy2 /* 0D_NOT_real out */
);
struct SuperBicubicInterpolation {
  double ansy;
  double ansy1;
  double ansy2;
};
SimUtils::SuperBicubicInterpolation super_bicubic_interpolation(
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
);

// Skipped unusable routine super_bracket_root:
// - Argument not defined: func (have: [])
// - Argument not defined: x1 (have: [])
// - Argument not defined: x2 (have: [])
// - Argument not defined: status (have: [])
// - Argument not defined: x_range (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_brent:
// - Argument not defined: ax (have: [])
// - Argument not defined: bx (have: [])
// - Argument not defined: cx (have: [])
// - Argument not defined: func (have: [])
// - Argument not defined: rel_tol (have: [])
// - Argument not defined: abs_tol (have: [])
// - Argument not defined: xmin (have: [])
// - Argument not defined: status (have: [])
// - Argument not defined: ymin (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_dbrent:
// - Argument not defined: ax (have: [])
// - Argument not defined: bx (have: [])
// - Argument not defined: cx (have: [])
// - Argument not defined: func (have: [])
// - Argument not defined: dfunc (have: [])
// - Argument not defined: rel_tol (have: [])
// - Argument not defined: abs_tol (have: [])
// - Argument not defined: xmin (have: [])
// - Argument not defined: status (have: [])
// - Argument not defined: func_min (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_gaussj:
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real

// Skipped unusable routine super_ludcmp:
// - Variable in sized array: 2D_NOT_real

// Skipped unusable routine super_mnbrak:
// - Argument not defined: ax (have: [])
// - Argument not defined: bx (have: [])
// - Argument not defined: cx (have: [])
// - Argument not defined: fa (have: [])
// - Argument not defined: fb (have: [])
// - Argument not defined: fc (have: [])
// - Argument not defined: func (have: [])
// - Argument not defined: status (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_mrqcof:
// - Argument not defined: a (have: [])
// - Argument not defined: y (have: [])
// - Argument not defined: co_alpha (have: [])
// - Argument not defined: da_beta (have: [])
// - Argument not defined: weight (have: [])
// - Argument not defined: chisq (have: [])
// - Argument not defined: funcs (have: [])
// - Argument not defined: storage (have: [])
// - Argument not defined: status (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_mrqmin:
// - Argument not defined: y (have: [])
// - Argument not defined: weight (have: [])
// - Argument not defined: a (have: [])
// - Argument not defined: chisq (have: [])
// - Argument not defined: funcs (have: [])
// - Argument not defined: storage (have: [])
// - Argument not defined: alamda (have: [])
// - Argument not defined: status (have: [])
// - Argument not defined: maska (have: [])
// - Argument not defined: print_err (have: [])
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_super_polint(
    Bmad::array_descriptor_t &xa /* 1D_NOT_real in */,
    Bmad::array_descriptor_t &ya /* 1D_NOT_real in */,
    double &x /* 0D_NOT_real in */,
    double &y /* 0D_NOT_real out */,
    double &dy /* 0D_NOT_real out */
);
struct SuperPolint {
  double y;
  double dy;
};
SimUtils::SuperPolint super_polint(FArray1D<Real> &xa, FArray1D<Real> &ya, double x);
extern "C" bool fortran_super_poly(
    double &x /* 0D_NOT_real in */,
    Bmad::array_descriptor_t &coeffs /* 1D_NOT_real in */,
    double &value /* 0D_NOT_real out */
);
double super_poly(double x, FArray1D<Real> &coeffs);

// Skipped unusable routine super_qromb:
// - Argument not defined: func (have: [])
// - Argument not defined: a (have: [])
// - Argument not defined: b (have: [])
// - Argument not defined: rel_tol (have: [])
// - Argument not defined: abs_tol (have: [])
// - Argument not defined: k_order (have: [])
// - Argument not defined: err_flag (have: [])
// - Argument not defined: integral (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_qromb_2d:
// - Argument not defined: func (have: [])
// - Argument not defined: ax (have: [])
// - Argument not defined: bx (have: [])
// - Argument not defined: ay (have: [])
// - Argument not defined: by (have: [])
// - Argument not defined: rel_tol (have: [])
// - Argument not defined: abs_tol (have: [])
// - Argument not defined: k_order (have: [])
// - Argument not defined: err_flag (have: [])
// - Argument not defined: integral (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_rtsafe:
// - Argument not defined: funcs (have: [])
// - Argument not defined: x1 (have: [])
// - Argument not defined: x2 (have: [])
// - Argument not defined: rel_tol (have: [])
// - Argument not defined: abs_tol (have: [])
// - Argument not defined: status (have: [])
// - Argument not defined: x_zero (have: [])
// - Translated arg count mismatch (unsupported?)
extern "C" void fortran_super_sobseq(
    Bmad::array_descriptor_t &x /* 1D_NOT_real inout */,
    void *ran_state /* 0D_NOT_type inout */
);
void super_sobseq(FArray1D<Real> &x, optional_ref<RandomStateStruct> ran_state = std::nullopt);
extern "C" void fortran_super_sort(Bmad::array_descriptor_t &arr /* 1D_NOT_integer inout */);
void super_sort(FArray1D<Int> &arr);

// Skipped unusable routine super_trapzd:
// - Argument not defined: func (have: [])
// - Argument not defined: a (have: [])
// - Argument not defined: b (have: [])
// - Argument not defined: s (have: [])
// - Argument not defined: n (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_zbrent:
// - Argument not defined: func (have: [])
// - Argument not defined: x1 (have: [])
// - Argument not defined: x2 (have: [])
// - Argument not defined: rel_tol (have: [])
// - Argument not defined: abs_tol (have: [])
// - Argument not defined: status (have: [])
// - Argument not defined: func_val (have: [])
// - Argument not defined: x_zero (have: [])
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine svd_fit:
// - Variable inout sized array: 2D_NOT_real
// - Variable inout sized array: 2D_NOT_real
extern "C" void fortran_system_command(
    const char *line /* 0D_NOT_character in */,
    bool *err_flag /* 0D_NOT_logical in */
);
void system_command(std::string line, std::optional<bool> err_flag = std::nullopt);
extern "C" bool fortran_to_str(
    double &num /* 0D_NOT_real in */,
    int *max_signif /* 0D_NOT_integer in */,
    const char *string /* 0D_ALLOC_character in */
);
void to_str(double num, std::optional<int> max_signif, std::string string);
extern "C" bool fortran_tricubic_cmplx_eval(
    double &x_norm /* 0D_NOT_real in */,
    double &y_norm /* 0D_NOT_real in */,
    double &z_norm /* 0D_NOT_real in */,
    void *tri_coef /* 0D_NOT_type in */,
    std::complex<double> &df_dx /* 0D_NOT_complex out */,
    std::complex<double> &df_dy /* 0D_NOT_complex out */,
    std::complex<double> &df_dz /* 0D_NOT_complex out */,
    std::complex<double> &f_val /* 0D_NOT_complex out */
);
struct TricubicCmplxEval {
  std::complex<double> df_dx;
  std::complex<double> df_dy;
  std::complex<double> df_dz;
  std::complex<double> f_val;
};
SimUtils::TricubicCmplxEval
tricubic_cmplx_eval(double x_norm, double y_norm, double z_norm, TricubicCmplxCoefStruct &tri_coef);

// Skipped unusable routine tricubic_compute_cmplx_field_at_3d_box:
// - Untranslated type: cmplx_field_at_3d_box_struct (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tricubic_compute_field_at_3d_box:
// - Untranslated type: field_at_3d_box_struct (0D)
// - Translated arg count mismatch (unsupported?)

// Skipped unusable routine tricubic_eval:
// - Untranslated type: tricubic_coef_struct (0D)

// Skipped unusable routine tricubic_interpolation_cmplx_coefs:
// - Untranslated type: cmplx_field_at_3d_box_struct (0D)

// Skipped unusable routine tricubic_interpolation_coefs:
// - Untranslated type: field_at_3d_box_struct (0D)
// - Untranslated type: tricubic_coef_struct (0D)
extern "C" void fortran_type_this_file(const char *filename /* 0D_NOT_character in */);
void type_this_file(std::string filename);

// Skipped unusable routine unquote:
// - No matching docstring

// Skipped unusable routine upcase:
// - No matching docstring
extern "C" void fortran_upcase_string(const char *string /* 0D_NOT_character in */);
void upcase_string(std::string string);

// Skipped unusable routine value_of_all_ptr:
// - Untranslated type: all_pointer_struct (0D)
extern "C" bool fortran_virtual_memory_usage(int &usage /* 0D_NOT_integer out */);
int virtual_memory_usage();
extern "C" void fortran_w_mat_to_axis_angle(
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &axis /* 1D_NOT_real out */,
    double &angle /* 0D_NOT_real out */
);
struct WMatToAxisAngle {
  FixedArray1D<Real, 3> axis;
  double angle;
};
SimUtils::WMatToAxisAngle w_mat_to_axis_angle(FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_to_quat(
    Bmad::array_descriptor_t &w_mat /* 2D_NOT_real in */,
    Bmad::array_descriptor_t &quat /* 1D_NOT_real out */
);
FixedArray1D<Real, 4> w_mat_to_quat(FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool
fortran_word_len(const char *wording /* 0D_NOT_character in */, int &wlen /* 0D_NOT_integer in */);
void word_len(std::string wording, int wlen);
extern "C" void fortran_word_read(
    const char *in_str /* 0D_NOT_character in */,
    const char *delim_list /* 0D_NOT_character in */,
    const char *word /* 0D_NOT_character in */,
    int &ix_word /* 0D_NOT_integer in */,
    const char *delim /* 0D_NOT_character in */,
    bool &delim_found /* 0D_NOT_logical in */,
    const char *out_str /* 0D_NOT_character in */,
    bool *ignore_interior /* 0D_NOT_logical in */
);
void word_read(
    std::string in_str,
    std::string delim_list,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    std::string out_str,
    std::optional<bool> ignore_interior = std::nullopt
);
extern "C" bool
fortran_x0_radiation_length(int &species /* 0D_NOT_integer in */, double &x0 /* 0D_NOT_real out */);
double x0_radiation_length(int species);
} // namespace SimUtils
