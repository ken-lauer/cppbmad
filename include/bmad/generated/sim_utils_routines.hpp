#pragma once

#include <functional>

#include "bmad/convert.h"
#include "bmad/generated/enums.h"
#include "bmad/generated/proxy.hpp"
#include "bmad/types.h"

using namespace Bmad;

namespace SimUtils {

// Skipped unusable routine all_pointer_to_string:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" void fortran_allocate_thread_states();
void allocate_thread_states();
extern "C" bool fortran_anomalous_moment_of(
    c_Int& species /* 0D_NOT_integer */,
    c_Real& moment /* 0D_NOT_real */);
double anomalous_moment_of(int species);
extern "C" bool fortran_antiparticle(
    c_Int& species /* 0D_NOT_integer */,
    c_Int& anti_species /* 0D_NOT_integer */);
int antiparticle(int species);
extern "C" void fortran_apfft(
    void* rdata_in /* 1D_ALLOC_real */,
    c_RealArr bounds /* 1D_NOT_real */,
    c_Char window /* 0D_NOT_character */,
    c_Real& phase /* 0D_NOT_real */,
    c_Int* diag /* 0D_NOT_integer */);
void apfft(
    RealAlloc1D& rdata_in,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    optional_ref<int> diag = std::nullopt);
extern "C" void fortran_apfft_corr(
    void* rdata_in /* 1D_ALLOC_real */,
    c_RealArr bounds /* 1D_NOT_real */,
    c_Char window /* 0D_NOT_character */,
    c_Real& phase /* 0D_NOT_real */,
    c_Real& amp /* 0D_NOT_real */,
    c_Real& freq /* 0D_NOT_real */,
    c_Int* diag /* 0D_NOT_integer */);
struct ApfftCorr {
  double phase;
  double amp;
  double freq;
};
ApfftCorr apfft_corr(
    RealAlloc1D& rdata_in,
    std::optional<FixedArray1D<Real, 2>> bounds,
    std::string window,
    std::optional<int> diag = std::nullopt);
extern "C" void fortran_apfft_ext(
    void* rdata /* 1D_ALLOC_real */,
    c_RealArr bounds /* 1D_NOT_real */,
    c_Char window /* 0D_NOT_character */,
    c_Real& phase /* 0D_NOT_real */,
    c_Real& amp /* 0D_NOT_real */,
    c_Real& freq /* 0D_NOT_real */,
    c_Int* diag /* 0D_NOT_integer */);
void apfft_ext(
    RealAlloc1D& rdata,
    FixedArray1D<Real, 2> bounds,
    std::string window,
    double phase,
    double amp,
    double freq,
    optional_ref<int> diag = std::nullopt);
extern "C" bool fortran_asinc(
    c_Real& x /* 0D_NOT_real */,
    c_Int* nd /* 0D_NOT_integer */,
    c_Real& y /* 0D_NOT_real */);
void asinc(double x, optional_ref<int> nd, double y);
extern "C" bool fortran_assert_equal(
    void* int_arr /* 1D_ALLOC_integer */,
    c_Char err_str /* 0D_NOT_character */,
    c_Int& ival /* 0D_NOT_integer */);
void assert_equal(IntAlloc1D& int_arr, std::string err_str, int ival);
extern "C" bool fortran_atomic_number(
    c_Int& species /* 0D_NOT_integer */,
    c_Int& atomic_num /* 0D_NOT_integer */);
int atomic_number(int species);
extern "C" bool fortran_atomic_species_id(
    c_Int& charge /* 0D_NOT_integer */,
    c_Bool& is_anti /* 0D_NOT_logical */,
    c_Int& atomic_num /* 0D_NOT_integer */,
    c_Int& n_nuc /* 0D_NOT_integer */,
    c_Int& species_id /* 0D_NOT_integer */);
int atomic_species_id(int charge, bool is_anti, int atomic_num, int n_nuc);
extern "C" bool fortran_axis_angle_to_quat(
    c_RealArr axis /* 1D_NOT_real */,
    c_Real& angle /* 0D_NOT_real */,
    c_RealArr quat /* 1D_NOT_real */);
FixedArray1D<Real, 4> axis_angle_to_quat(
    FixedArray1D<Real, 3> axis,
    double angle);
extern "C" void fortran_axis_angle_to_w_mat(
    c_RealArr axis /* 1D_NOT_real */,
    c_Real& angle /* 0D_NOT_real */,
    c_RealArr w_mat /* 2D_NOT_real */);
FixedArray2D<Real, 3, 3> axis_angle_to_w_mat(
    FixedArray1D<Real, 3> axis,
    double angle);
extern "C" bool fortran_bicubic_cmplx_eval(
    c_Real& x_norm /* 0D_NOT_real */,
    c_Real& y_norm /* 0D_NOT_real */,
    void* bi_coef /* 0D_NOT_type */,
    c_Complex& df_dx /* 0D_NOT_complex */,
    c_Complex& df_dy /* 0D_NOT_complex */,
    c_Complex& f_val /* 0D_NOT_complex */);
struct BicubicCmplxEval {
  std::complex<double> df_dx;
  std::complex<double> df_dy;
  std::complex<double> f_val;
};
BicubicCmplxEval bicubic_cmplx_eval(
    double x_norm,
    double y_norm,
    BicubicCmplxCoefProxy& bi_coef);

// Skipped unusable routine bicubic_compute_cmplx_field_at_2d_box:
// Untranslated type: CmplxFieldAt2dBoxProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine bicubic_compute_field_at_2d_box:
// Untranslated type: FieldAt2dBoxProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine bicubic_eval:
// Untranslated type: BicubicCoefProxy (0D_NOT_type)

// Skipped unusable routine bicubic_interpolation_cmplx_coefs:
// Untranslated type: CmplxFieldAt2dBoxProxy (0D_NOT_type)

// Skipped unusable routine bicubic_interpolation_coefs:
// Untranslated type: FieldAt2dBoxProxy (0D_NOT_type)
// Untranslated type: BicubicCoefProxy (0D_NOT_type)

// Skipped unusable routine bin_2d:
// Untranslated type: GeneralBinProxy (0D_NOT_type)

// Skipped unusable routine bin_data:
// Untranslated type: BinProxy (0D_NOT_type)

// Skipped unusable routine bin_data_density:
// Untranslated type: BinProxy (0D_NOT_type)

// Skipped unusable routine bin_data_density_2d:
// Untranslated type: GeneralBinProxy (0D_NOT_type)
extern "C" bool fortran_bin_index(
    c_Real& x /* 0D_NOT_real */,
    c_Real& bin1_x_min /* 0D_NOT_real */,
    c_Real& bin_delta /* 0D_NOT_real */,
    c_Int& ix_bin /* 0D_NOT_integer */);
int bin_index(double x, double bin1_x_min, double bin_delta);
extern "C" bool fortran_bin_x_center(
    c_Int& ix_bin /* 0D_NOT_integer */,
    c_Real& bin1_x_min /* 0D_NOT_real */,
    c_Real& bin_delta /* 0D_NOT_real */,
    c_Real& x_center /* 0D_NOT_real */);
double bin_x_center(int ix_bin, double bin1_x_min, double bin_delta);
extern "C" void fortran_bit_set(
    c_Int& word /* 0D_NOT_integer */,
    c_Int& pos /* 0D_NOT_integer */,
    c_Bool& set_to_1 /* 0D_NOT_logical */);
void bit_set(int word, int pos, bool set_to_1);

// Skipped unusable routine bracket_index:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine bracket_index2:
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_bracket_index_for_spline(
    void* x_knot /* 1D_ALLOC_real */,
    c_Real& x /* 0D_NOT_real */,
    c_Int& ix0 /* 0D_NOT_integer */,
    c_Bool* strict /* 0D_NOT_logical */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_Bool& ok /* 0D_NOT_logical */);
struct BracketIndexForSpline {
  int ix0;
  bool ok;
};
BracketIndexForSpline bracket_index_for_spline(
    RealAlloc1D& x_knot,
    double x,
    std::optional<bool> strict = std::nullopt,
    std::optional<bool> print_err = std::nullopt);

// Skipped unusable routine bracket_index_int:
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_calc_file_number(
    c_Char file_name /* 0D_NOT_character */,
    c_Int& num_in /* 0D_NOT_integer */,
    c_Int& num_out /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void calc_file_number(
    std::string file_name,
    int num_in,
    int num_out,
    bool err_flag);
extern "C" void fortran_change_file_number(
    c_Char file_name /* 0D_NOT_character */,
    c_Int& change /* 0D_NOT_integer */);
void change_file_number(std::string file_name, int change);
extern "C" bool fortran_charge_of(
    c_Int& species /* 0D_NOT_integer */,
    c_Int* default_ /* 0D_NOT_integer */,
    c_Int& charge /* 0D_NOT_integer */);
int charge_of(int species, std::optional<int> default_ = std::nullopt);
extern "C" bool fortran_charge_to_mass_of(
    c_Int& species /* 0D_NOT_integer */,
    c_Real& charge_mass_ratio /* 0D_NOT_real */);
double charge_to_mass_of(int species);
extern "C" bool fortran_coarse_frequency_estimate(
    void* data /* 1D_ALLOC_real */,
    c_Bool* error /* 0D_NOT_logical */,
    c_Real& frequency /* 0D_NOT_real */);
double coarse_frequency_estimate(
    RealAlloc1D& data,
    optional_ref<bool> error = std::nullopt);
extern "C" void fortran_complex_error_function(
    c_Real& wr /* 0D_NOT_real */,
    c_Real& wi /* 0D_NOT_real */,
    c_Real& zr /* 0D_NOT_real */,
    c_Real& zi /* 0D_NOT_real */);
void complex_error_function(double wr, double wi, double zr, double zi);
extern "C" bool fortran_cos_one(
    c_Real& angle /* 0D_NOT_real */,
    c_Real& cos1 /* 0D_NOT_real */);
void cos_one(double angle, double cos1);
extern "C" bool fortran_cosc(
    c_Real& x /* 0D_NOT_real */,
    c_Int* nd /* 0D_NOT_integer */,
    c_Real& y /* 0D_NOT_real */);
void cosc(double x, optional_ref<int> nd, double y);

// Skipped unusable routine count_at_index:
// Untranslated type: BinProxy (0D_NOT_type)

// Skipped unusable routine covar_expand:
// Variable inout sized array: covar(:,:) 2D_NOT_real

// Skipped unusable routine cplx_mat_inverse:
// Variable inout sized array: mat(:,:) 2D_NOT_complex
// Variable inout sized array: mat_inv(:,:) 2D_NOT_complex

// Skipped unusable routine cplx_mat_make_unit:
// Variable inout sized array: mat(:,:) 2D_NOT_complex
extern "C" bool fortran_create_a_spline(
    void* r0 /* 1D_ALLOC_real */,
    void* r1 /* 1D_ALLOC_real */,
    c_Real& slope0 /* 0D_NOT_real */,
    c_Real& slope1 /* 0D_NOT_real */,
    void* spline /* 0D_NOT_type */);
SplineProxy create_a_spline(
    RealAlloc1D& r0,
    RealAlloc1D& r1,
    double slope0,
    double slope1);
extern "C" bool fortran_cross_product(
    void* a /* 1D_ALLOC_real */,
    void* b /* 1D_ALLOC_real */,
    c_RealArr c /* 1D_NOT_real */);
void cross_product(RealAlloc1D& a, RealAlloc1D& b, FixedArray1D<Real, 3> c);

// Skipped unusable routine da2_div:
// Variable in sized array: ta(0:,0:) 2D_NOT_real
// Variable in sized array: tb(0:,0:) 2D_NOT_real
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine da2_evaluate:
// Variable in sized array: ta(0:,0:) 2D_NOT_real

// Skipped unusable routine da2_inverse:
// Variable in sized array: ta(0:,0:) 2D_NOT_real
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine da2_mult:
// Variable in sized array: ta(0:,0:) 2D_NOT_real
// Variable in sized array: tb(0:,0:) 2D_NOT_real
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_date_and_time_stamp(
    c_Char string /* 0D_NOT_character */,
    c_Bool* numeric_month /* 0D_NOT_logical */,
    c_Bool* include_zone /* 0D_NOT_logical */);
void date_and_time_stamp(
    std::string string,
    optional_ref<bool> numeric_month = std::nullopt,
    optional_ref<bool> include_zone = std::nullopt);
extern "C" void fortran_destfixedwindowls(c_Int& id /* 0D_NOT_integer */);
void destfixedwindowls(int id);
extern "C" void fortran_detab(c_Char str /* 0D_NOT_character */);
void detab(std::string str);

// Skipped unusable routine determinant:
// Variable inout sized array: mat(:,:) 2D_NOT_real
extern "C" void fortran_display_size_and_resolution(
    c_Int& ix_screen /* 0D_NOT_integer */,
    c_Real& x_size /* 0D_NOT_real */,
    c_Real& y_size /* 0D_NOT_real */,
    c_Real& x_res /* 0D_NOT_real */,
    c_Real& y_res /* 0D_NOT_real */);
void display_size_and_resolution(
    int ix_screen,
    double x_size,
    double y_size,
    double x_res,
    double y_res);
extern "C" bool fortran_dj_bessel(
    c_Int& m /* 0D_NOT_integer */,
    c_Real& arg /* 0D_NOT_real */,
    c_Real& dj_bes /* 0D_NOT_real */);
void dj_bessel(int m, double arg, double dj_bes);
extern "C" bool fortran_djb_hash(
    c_Char str /* 0D_NOT_character */,
    c_Int* old_hash /* 0D_NOT_integer */,
    c_Int& hash /* 0D_NOT_integer */);
void djb_hash(std::string str, optional_ref<int> old_hash, int hash);
extern "C" bool fortran_djb_str_hash(
    c_Char in_str /* 0D_NOT_character */,
    c_Char hash_str /* 0D_NOT_character */);
void djb_str_hash(std::string in_str, std::string hash_str);

// Skipped unusable routine doubleup_quotes:
// No matching docstring

// Skipped unusable routine downcase:
// No matching docstring
extern "C" void fortran_downcase_string(c_Char string /* 0D_NOT_character */);
void downcase_string(std::string string);

// Skipped unusable routine ed:
// Routine in configuration skip list
extern "C" void fortran_end_akima_spline_calc(
    void* spline /* 1D_ALLOC_type */,
    c_Int& which_end /* 0D_NOT_integer */);
void end_akima_spline_calc(SplineProxyAlloc1D& spline, int which_end);
extern "C" void fortran_err_exit(c_Char err_str /* 0D_NOT_character */);
void err_exit(optional_ref<std::string> err_str = std::nullopt);
extern "C" bool fortran_factorial(
    c_Int& n /* 0D_NOT_integer */,
    c_Real& fact /* 0D_NOT_real */);
void factorial(int n, double fact);
extern "C" void fortran_faddeeva_function(
    c_RealArr z /* 1D_NOT_real */,
    c_RealArr w /* 1D_NOT_real */,
    c_RealArr dw /* 2D_NOT_real */);
void faddeeva_function(
    FixedArray1D<Real, 2> z,
    FixedArray1D<Real, 2> w,
    FixedArray2D<Real, 2, 2> dw);
extern "C" void fortran_fft_1d(
    void* arr /* 1D_ALLOC_complex */,
    c_Int& isign /* 0D_NOT_integer */);
void fft_1d(ComplexAlloc1D& arr, int isign);
extern "C" void fortran_file_directorizer(
    c_Char in_file /* 0D_NOT_character */,
    c_Char out_file /* 0D_NOT_character */,
    c_Char directory /* 0D_NOT_character */,
    c_Bool& add_switch /* 0D_NOT_logical */);
void file_directorizer(
    std::string in_file,
    std::string out_file,
    std::string directory,
    bool add_switch);
extern "C" void fortran_file_get(
    c_Char string /* 0D_NOT_character */,
    c_Char dflt_file_name /* 0D_NOT_character */,
    c_Char file_name /* 0D_NOT_character */);
void file_get(
    std::string string,
    std::string dflt_file_name,
    std::string file_name);
extern "C" void fortran_file_get_open(
    c_Char string /* 0D_NOT_character */,
    c_Char dflt_file_name /* 0D_NOT_character */,
    c_Char file_name /* 0D_NOT_character */,
    c_Int& file_unit /* 0D_NOT_integer */,
    c_Bool& readonly /* 0D_NOT_logical */);
void file_get_open(
    std::string string,
    std::string dflt_file_name,
    std::string file_name,
    int file_unit,
    bool readonly);
extern "C" void fortran_file_suffixer(
    c_Char in_file_name /* 0D_NOT_character */,
    c_Char out_file_name /* 0D_NOT_character */,
    c_Char suffix /* 0D_NOT_character */,
    c_Bool& add_switch /* 0D_NOT_logical */);
void file_suffixer(
    std::string in_file_name,
    std::string out_file_name,
    std::string suffix,
    bool add_switch);
extern "C" bool fortran_find_location_int(
    void* arr /* 1D_ALLOC_integer */,
    c_Int& value /* 0D_NOT_integer */,
    c_Int& ix_match /* 0D_NOT_integer */);
void find_location_int(IntAlloc1D& arr, int value, int ix_match);
extern "C" bool fortran_find_location_logic(
    void* arr /* 1D_ALLOC_logical */,
    c_Bool& value /* 0D_NOT_logical */,
    c_Int& ix_match /* 0D_NOT_integer */);
void find_location_logic(BoolAlloc1D& arr, bool value, int ix_match);
extern "C" bool fortran_find_location_real(
    void* arr /* 1D_ALLOC_real */,
    c_Real& value /* 0D_NOT_real */,
    c_Int& ix_match /* 0D_NOT_integer */);
void find_location_real(RealAlloc1D& arr, double value, int ix_match);

// Skipped unusable routine find_location_str:
// Variable-sized inout character array: arr(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_fine_frequency_estimate(
    void* data /* 1D_ALLOC_real */,
    c_Real& frequency /* 0D_NOT_real */);
double fine_frequency_estimate(RealAlloc1D& data);
extern "C" bool fortran_fixedwindowls(
    c_Real& ynew /* 0D_NOT_real */,
    c_Int& id /* 0D_NOT_integer */,
    c_Real& z /* 0D_NOT_real */);
void fixedwindowls(double ynew, int id, double z);
extern "C" void fortran_fourier_amplitude(
    void* data /* 1D_ALLOC_real */,
    c_Real& frequency /* 0D_NOT_real */,
    c_Real& cos_amp /* 0D_NOT_real */,
    c_Real& sin_amp /* 0D_NOT_real */,
    c_Real& dcos_amp /* 0D_NOT_real */,
    c_Real& dsin_amp /* 0D_NOT_real */);
struct FourierAmplitude {
  double cos_amp;
  double sin_amp;
  double dcos_amp;
  double dsin_amp;
};
FourierAmplitude fourier_amplitude(RealAlloc1D& data, double frequency);
extern "C" bool fortran_gen_complete_elliptic(
    c_Real& kc /* 0D_NOT_real */,
    c_Real& p /* 0D_NOT_real */,
    c_Real& c /* 0D_NOT_real */,
    c_Real& s /* 0D_NOT_real */,
    c_Real* err_tol /* 0D_NOT_real */,
    c_Real& value /* 0D_NOT_real */);
void gen_complete_elliptic(
    double kc,
    double p,
    double c,
    double s,
    optional_ref<double> err_tol,
    double value);

// Skipped unusable routine general_bin_count:
// Untranslated type: GeneralBinProxy (0D_NOT_type)

// Skipped unusable routine general_bin_index:
// Untranslated type: GeneralBinProxy (0D_NOT_type)

// Skipped unusable routine general_bin_index_in_bounds:
// Untranslated type: GeneralBinProxy (0D_NOT_type)
extern "C" void fortran_get_file_number(
    c_Char file_name /* 0D_NOT_character */,
    c_Char cnum_in /* 0D_NOT_character */,
    c_Int& num_out /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */);
void get_file_number(
    std::string file_name,
    std::string cnum_in,
    int num_out,
    bool err_flag);
extern "C" void fortran_get_file_time_stamp(
    c_Char file /* 0D_NOT_character */,
    c_Char time_stamp /* 0D_NOT_character */);
void get_file_time_stamp(std::string file, std::string time_stamp);

// Skipped unusable routine han:
// Routine in configuration skip list
extern "C" void fortran_hanhan(
    c_Int& N /* 0D_NOT_integer */,
    void* hh /* 1D_ALLOC_real */);
void hanhan(int N, RealAlloc1D& hh);
extern "C" bool fortran_i_bessel(
    c_Int& m /* 0D_NOT_integer */,
    c_Real& arg /* 0D_NOT_real */,
    c_Real& i_bes /* 0D_NOT_real */);
void i_bessel(int m, double arg, double i_bes);
extern "C" bool fortran_i_bessel_extended(
    c_Int& m /* 0D_NOT_integer */,
    c_Real& arg /* 0D_NOT_real */,
    c_Complex& i_bes /* 0D_NOT_complex */);
void i_bessel_extended(int m, double arg, std::complex<double> i_bes);
extern "C" void fortran_increment_file_number(
    c_Char file_name /* 0D_NOT_character */,
    c_Int& digits /* 0D_NOT_integer */,
    c_Int& number /* 0D_NOT_integer */,
    c_Char cnumber /* 0D_NOT_character */);
void increment_file_number(
    std::string file_name,
    int digits,
    int number,
    std::string cnumber);
extern "C" bool fortran_index_nocase(
    c_Char string1 /* 0D_NOT_character */,
    c_Char string2 /* 0D_NOT_character */,
    c_Int& indx /* 0D_NOT_integer */);
void index_nocase(std::string string1, std::string string2, int indx);
extern "C" bool fortran_initfixedwindowls(
    c_Int& N /* 0D_NOT_integer */,
    c_Real& dt /* 0D_NOT_real */,
    c_Int& order /* 0D_NOT_integer */,
    c_Int& der /* 0D_NOT_integer */,
    c_Int& id /* 0D_NOT_integer */);
int initfixedwindowls(int N, double dt, int order, int der);

// Skipped unusable routine int_logic:
// Routine in configuration skip list
extern "C" bool fortran_int_str(
    c_Int& int_ /* 0D_NOT_integer */,
    c_Int* width /* 0D_NOT_integer */,
    c_Char str /* 0D_ALLOC_character */);
void int_str(int int_, optional_ref<int> width, std::string str);
extern "C" bool fortran_interpolated_fft(
    void* cdata /* 1D_ALLOC_complex */,
    c_Bool& calc_ok /* 0D_NOT_logical */,
    c_Int* opt_dump_spectrum /* 0D_NOT_integer */,
    c_Int* opt_dump_index /* 0D_NOT_integer */,
    c_Real& this_fft /* 0D_NOT_real */);
void interpolated_fft(
    ComplexAlloc1D& cdata,
    bool calc_ok,
    optional_ref<int> opt_dump_spectrum,
    optional_ref<int> opt_dump_index,
    double this_fft);
extern "C" bool fortran_interpolated_fft_gsl(
    void* cdata /* 1D_ALLOC_complex */,
    c_Bool& calc_ok /* 0D_NOT_logical */,
    c_Int* opt_dump_spectrum /* 0D_NOT_integer */,
    c_Int* opt_dump_index /* 0D_NOT_integer */,
    c_Real& this_fft /* 0D_NOT_real */);
void interpolated_fft_gsl(
    ComplexAlloc1D& cdata,
    bool calc_ok,
    optional_ref<int> opt_dump_spectrum,
    optional_ref<int> opt_dump_index,
    double this_fft);

// Skipped unusable routine inverse:
// Argument not defined: funct (have: [])
// Argument not defined: y (have: [])
// Argument not defined: x1 (have: [])
// Argument not defined: x2 (have: [])
// Argument not defined: tol (have: [])
// Argument not defined: x (have: [])
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_is_alphabetic(
    c_Char string /* 0D_NOT_character */,
    c_Char valid_chars /* 0D_NOT_character */,
    c_Bool& is_alpha /* 0D_NOT_logical */);
void is_alphabetic(
    std::string string,
    optional_ref<std::string> valid_chars,
    bool is_alpha);
extern "C" bool fortran_is_decreasing_sequence(
    void* array /* 1D_ALLOC_real */,
    c_Bool* strict /* 0D_NOT_logical */,
    c_Bool& is_decreasing /* 0D_NOT_logical */);
void is_decreasing_sequence(
    RealAlloc1D& array,
    std::optional<bool> strict,
    bool is_decreasing);
extern "C" bool fortran_is_false(
    c_Real& param /* 0D_NOT_real */,
    c_Bool& this_false /* 0D_NOT_logical */);
bool is_false(double param);
extern "C" bool fortran_is_increasing_sequence(
    void* array /* 1D_ALLOC_real */,
    c_Bool* strict /* 0D_NOT_logical */,
    c_Bool& is_increasing /* 0D_NOT_logical */);
void is_increasing_sequence(
    RealAlloc1D& array,
    std::optional<bool> strict,
    bool is_increasing);
extern "C" bool fortran_is_integer(
    c_Char string /* 0D_NOT_character */,
    c_Int* int_ /* 0D_NOT_integer */,
    c_Char delims /* 0D_NOT_character */,
    c_Int* ix_word /* 0D_NOT_integer */,
    c_Bool& valid /* 0D_NOT_logical */);
void is_integer(
    std::string string,
    optional_ref<int> int_,
    optional_ref<std::string> delims,
    optional_ref<int> ix_word,
    bool valid);
extern "C" bool fortran_is_logical(
    c_Char string /* 0D_NOT_character */,
    c_Bool* ignore /* 0D_NOT_logical */,
    c_Bool& valid /* 0D_NOT_logical */);
void is_logical(std::string string, optional_ref<bool> ignore, bool valid);
extern "C" bool fortran_is_real(
    c_Char string /* 0D_NOT_character */,
    c_Bool* ignore /* 0D_NOT_logical */,
    c_Real* real_num /* 0D_NOT_real */,
    c_Bool& valid /* 0D_NOT_logical */);
void is_real(
    std::string string,
    optional_ref<bool> ignore,
    optional_ref<double> real_num,
    bool valid);
extern "C" bool fortran_is_subatomic_species(
    c_Int& species /* 0D_NOT_integer */,
    c_Bool& is_subatomic /* 0D_NOT_logical */);
bool is_subatomic_species(int species);
extern "C" bool fortran_is_true(
    c_Real& param /* 0D_NOT_real */,
    c_Bool& this_true /* 0D_NOT_logical */);
bool is_true(double param);

// Skipped unusable routine isatty:
// Routine in configuration skip list
extern "C" bool fortran_j_bessel(
    c_Int& m /* 0D_NOT_integer */,
    c_Real& arg /* 0D_NOT_real */,
    c_Real& j_bes /* 0D_NOT_real */);
void j_bessel(int m, double arg, double j_bes);
extern "C" void fortran_linear_fit(
    void* x /* 1D_ALLOC_real */,
    void* y /* 1D_ALLOC_real */,
    c_Int& n_data /* 0D_NOT_integer */,
    c_Real& a /* 0D_NOT_real */,
    c_Real& b /* 0D_NOT_real */,
    c_Real& sig_a /* 0D_NOT_real */,
    c_Real& sig_b /* 0D_NOT_real */);
void linear_fit(
    RealAlloc1D& x,
    RealAlloc1D& y,
    int n_data,
    double a,
    double b,
    double sig_a,
    double sig_b);
extern "C" void fortran_linear_fit_2d(
    void* x /* 1D_ALLOC_real */,
    void* y /* 1D_ALLOC_real */,
    void* z /* 1D_ALLOC_real */,
    c_RealArr coef /* 1D_NOT_real */);
FixedArray1D<Real, 3> linear_fit_2d(
    RealAlloc1D& x,
    RealAlloc1D& y,
    RealAlloc1D& z);

// Skipped unusable routine location_decode:
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_logic_str(
    c_Bool& logic /* 0D_NOT_logical */,
    c_Char str /* 0D_NOT_character */);
void logic_str(bool logic, std::string str);
extern "C" bool fortran_lunget(c_Int& func_retval__ /* 0D_NOT_integer */);
int lunget();
extern "C" void fortran_make_legal_comment(
    c_Char comment_in /* 0D_NOT_character */,
    c_Char comment_out /* 0D_NOT_character */);
void make_legal_comment(std::string comment_in, std::string comment_out);
extern "C" bool fortran_mass_of(
    c_Int& species /* 0D_NOT_integer */,
    c_Real& mass /* 0D_NOT_real */);
double mass_of(int species);

// Skipped unusable routine mat_eigen:
// Variable inout sized array: mat(:,:) 2D_NOT_real
// Variable inout sized array: eigen_vec(:,:) 2D_NOT_complex

// Skipped unusable routine mat_inverse:
// Variable inout sized array: mat(:,:) 2D_NOT_real
// Variable inout sized array: mat_inv(:,:) 2D_NOT_real

// Skipped unusable routine mat_make_unit:
// Variable inout sized array: mat(:,:) 2D_NOT_real

// Skipped unusable routine mat_pseudoinverse:
// Variable inout sized array: A(:,:) 2D_NOT_real
// Variable inout sized array: Ap(:,:) 2D_NOT_real

// Skipped unusable routine mat_rotation:
// Variable inout sized array: mat(:,:) 2D_NOT_real

// Skipped unusable routine mat_scale_p0:
// Variable inout sized array: mat_in(:,:) 2D_NOT_real
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine mat_symp_conj:
// Variable inout sized array: mat(:,:) 2D_NOT_real
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine mat_symp_conj_i:
// Variable inout sized array: mat(:,:) 2D_NOT_complex
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine mat_symp_error:
// Variable inout sized array: mat(:,:) 2D_NOT_real
// Variable inout sized array: err_mat(:,:) 2D_NOT_real

// Skipped unusable routine mat_symplectify:
// Variable inout sized array: mat_in(:,:) 2D_NOT_real
// Variable inout sized array: mat_symp(:,:) 2D_NOT_real

// Skipped unusable routine mat_type:
// Variable inout sized array: mat(:,:) 2D_NOT_real
// Variable-sized inout character array: lines(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)
extern "C" bool fortran_match_reg(
    c_Char str /* 0D_NOT_character */,
    c_Char pat /* 0D_NOT_character */,
    c_Bool& is_match /* 0D_NOT_logical */);
void match_reg(std::string str, std::string pat, bool is_match);
extern "C" bool fortran_match_wild(
    c_Char string /* 0D_NOT_character */,
    c_Char template_ /* 0D_NOT_character */,
    c_Bool& is_match /* 0D_NOT_logical */);
void match_wild(std::string string, std::string template_, bool is_match);

// Skipped unusable routine match_word:
// Variable-sized inout character array: names(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine max_nonzero:
// Module name unset? Internal error
extern "C" bool fortran_maximize_projection(
    c_Real& seed /* 0D_NOT_real */,
    void* cdata /* 1D_ALLOC_complex */,
    c_Real& func_retval__ /* 0D_NOT_real */);
void maximize_projection(
    double seed,
    ComplexAlloc1D& cdata,
    double func_retval__);
extern "C" void fortran_milli_sleep(c_Int& milli_sec /* 0D_NOT_integer */);
void milli_sleep(int milli_sec);

// Skipped unusable routine molecular_components:
// Untranslated type: MolecularComponentProxy (1D_ALLOC_type)
extern "C" bool fortran_n_bins_automatic(
    c_Int& n_data /* 0D_NOT_integer */,
    c_Int& n /* 0D_NOT_integer */);
void n_bins_automatic(int n_data, int n);
extern "C" bool fortran_n_choose_k(
    c_Int& n /* 0D_NOT_integer */,
    c_Int& k /* 0D_NOT_integer */,
    c_Real& nck /* 0D_NOT_real */);
void n_choose_k(int n, int k, double nck);
extern "C" void fortran_n_spline_create(
    void* deriv0 /* 1D_ALLOC_real */,
    void* deriv1 /* 1D_ALLOC_real */,
    c_Real& x1 /* 0D_NOT_real */,
    void* n_spline /* 1D_ALLOC_real */);
RealAlloc1D n_spline_create(
    RealAlloc1D& deriv0,
    RealAlloc1D& deriv1,
    double x1);
extern "C" void fortran_naff(
    void* cdata /* 1D_ALLOC_complex */,
    void* freqs /* 1D_ALLOC_real */,
    void* amps /* 1D_ALLOC_complex */,
    c_Int* opt_dump_spectra /* 0D_NOT_integer */,
    c_Bool* opt_zero_first /* 0D_NOT_logical */);
void naff(
    ComplexAlloc1D& cdata,
    RealAlloc1D& freqs,
    ComplexAlloc1D& amps,
    optional_ref<int> opt_dump_spectra = std::nullopt,
    optional_ref<bool> opt_zero_first = std::nullopt);

// Skipped unusable routine nametable_add:
// Untranslated type: NametableProxy (0D_NOT_type)

// Skipped unusable routine nametable_bracket_indexx:
// Untranslated type: NametableProxy (0D_NOT_type)

// Skipped unusable routine nametable_change1:
// Untranslated type: NametableProxy (0D_NOT_type)

// Skipped unusable routine nametable_init:
// Untranslated type: NametableProxy (0D_NOT_type)

// Skipped unusable routine nametable_remove:
// Untranslated type: NametableProxy (0D_NOT_type)

// Skipped unusable routine node_put:
// Routine in configuration skip list
extern "C" bool fortran_omega_to_quat(
    c_RealArr omega /* 1D_NOT_real */,
    c_RealArr quat /* 1D_NOT_real */);
FixedArray1D<Real, 4> omega_to_quat(FixedArray1D<Real, 3> omega);
extern "C" bool fortran_openpmd_species_name(
    c_Int& species /* 0D_NOT_integer */,
    c_Char pmd_name /* 0D_NOT_character */);
std::string openpmd_species_name(int species);
extern "C" bool fortran_ordinal_str(
    c_Int& n /* 0D_NOT_integer */,
    c_Char str /* 0D_ALLOC_character */);
void ordinal_str(int n, std::string str);

// Skipped unusable routine outer_product:
// Module name unset? Internal error
extern "C" void fortran_parse_fortran_format(
    c_Char format_str /* 0D_NOT_character */,
    c_Int& n_repeat /* 0D_NOT_integer */,
    c_Int& power /* 0D_NOT_integer */,
    c_Char descrip /* 0D_NOT_character */,
    c_Int& width /* 0D_NOT_integer */,
    c_Int& digits /* 0D_NOT_integer */);
void parse_fortran_format(
    std::string format_str,
    int n_repeat,
    int power,
    std::string descrip,
    int width,
    int digits);

// Skipped unusable routine pointer_to_locations:
// Variable-sized inout character array: names(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine pointer_to_ran_state:
// Untranslated type: RandomStateProxy (0D_NOT_type)
// Untranslated type: RandomStateProxy (0D_PTR_type)
extern "C" bool fortran_poly_eval(
    void* poly /* 1D_ALLOC_real */,
    c_Real& x /* 0D_NOT_real */,
    c_Bool* diff_coef /* 0D_NOT_logical */,
    c_Real& y /* 0D_NOT_real */);
void poly_eval(
    RealAlloc1D& poly,
    double x,
    std::optional<bool> diff_coef,
    double y);
extern "C" bool fortran_probability_funct(
    c_Real& x /* 0D_NOT_real */,
    c_Real& prob /* 0D_NOT_real */);
void probability_funct(double x, double prob);
extern "C" bool fortran_projdd(
    void* a /* 1D_ALLOC_complex */,
    void* b /* 1D_ALLOC_complex */,
    c_Complex& func_retval__ /* 0D_NOT_complex */);
void projdd(
    ComplexAlloc1D& a,
    ComplexAlloc1D& b,
    std::complex<double> func_retval__);
extern "C" bool fortran_quadratic_roots(
    c_RealArr coefs /* 1D_NOT_real */,
    c_ComplexArr root /* 1D_NOT_complex */);
void quadratic_roots(
    FixedArray1D<Real, 3> coefs,
    FixedArray1D<Complex, 2> root);
extern "C" bool fortran_quat_conj_complex(
    c_ComplexArr q_in /* 1D_NOT_complex */,
    c_ComplexArr q_out /* 1D_NOT_complex */);
FixedArray1D<Complex, 4> quat_conj_complex(FixedArray1D<Complex, 4> q_in);
extern "C" bool fortran_quat_conj_real(
    c_RealArr q_in /* 1D_NOT_real */,
    c_RealArr q_out /* 1D_NOT_real */);
FixedArray1D<Real, 4> quat_conj_real(FixedArray1D<Real, 4> q_in);
extern "C" bool fortran_quat_inverse(
    c_RealArr q_in /* 1D_NOT_real */,
    c_RealArr q_out /* 1D_NOT_real */);
FixedArray1D<Real, 4> quat_inverse(FixedArray1D<Real, 4> q_in);
extern "C" bool fortran_quat_mul_complex(
    c_ComplexArr q1 /* 1D_NOT_complex */,
    c_ComplexArr q2 /* 1D_NOT_complex */,
    c_ComplexArr q3 /* 1D_NOT_complex */,
    c_ComplexArr q4 /* 1D_NOT_complex */,
    c_ComplexArr q5 /* 1D_NOT_complex */,
    c_ComplexArr q6 /* 1D_NOT_complex */,
    c_ComplexArr q7 /* 1D_NOT_complex */,
    c_ComplexArr q8 /* 1D_NOT_complex */,
    c_ComplexArr q9 /* 1D_NOT_complex */,
    c_ComplexArr q_out /* 1D_NOT_complex */);
FixedArray1D<Complex, 4> quat_mul_complex(
    FixedArray1D<Complex, 4> q1,
    FixedArray1D<Complex, 4> q2,
    std::optional<FixedArray1D<Complex, 4>> q3 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q4 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q5 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q6 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q7 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q8 = std::nullopt,
    std::optional<FixedArray1D<Complex, 4>> q9 = std::nullopt);
extern "C" bool fortran_quat_mul_real(
    c_RealArr q1 /* 1D_NOT_real */,
    c_RealArr q2 /* 1D_NOT_real */,
    c_RealArr q3 /* 1D_NOT_real */,
    c_RealArr q4 /* 1D_NOT_real */,
    c_RealArr q5 /* 1D_NOT_real */,
    c_RealArr q6 /* 1D_NOT_real */,
    c_RealArr q7 /* 1D_NOT_real */,
    c_RealArr q8 /* 1D_NOT_real */,
    c_RealArr q9 /* 1D_NOT_real */,
    c_RealArr q_out /* 1D_NOT_real */);
FixedArray1D<Real, 4> quat_mul_real(
    FixedArray1D<Real, 4> q1,
    FixedArray1D<Real, 4> q2,
    std::optional<FixedArray1D<Real, 4>> q3 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q4 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q5 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q6 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q7 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q8 = std::nullopt,
    std::optional<FixedArray1D<Real, 4>> q9 = std::nullopt);
extern "C" bool fortran_quat_rotate_complex(
    c_ComplexArr quat /* 1D_NOT_complex */,
    c_ComplexArr vec_in /* 1D_NOT_complex */,
    c_ComplexArr vec_out /* 1D_NOT_complex */);
FixedArray1D<Complex, 3> quat_rotate_complex(
    FixedArray1D<Complex, 4> quat,
    FixedArray1D<Complex, 3> vec_in);
extern "C" bool fortran_quat_rotate_real(
    c_RealArr quat /* 1D_NOT_real */,
    c_RealArr vec_in /* 1D_NOT_real */,
    c_RealArr vec_out /* 1D_NOT_real */);
FixedArray1D<Real, 3> quat_rotate_real(
    FixedArray1D<Real, 4> quat,
    FixedArray1D<Real, 3> vec_in);
extern "C" void fortran_quat_to_axis_angle(
    c_RealArr quat /* 1D_NOT_real */,
    c_RealArr axis /* 1D_NOT_real */,
    c_Real& angle /* 0D_NOT_real */);
struct QuatToAxisAngle {
  FixedArray1D<Real, 3> axis;
  double angle;
};
QuatToAxisAngle quat_to_axis_angle(FixedArray1D<Real, 4> quat);
extern "C" bool fortran_quat_to_omega(
    c_RealArr quat /* 1D_NOT_real */,
    c_RealArr omega /* 1D_NOT_real */);
FixedArray1D<Real, 3> quat_to_omega(FixedArray1D<Real, 4> quat);
extern "C" bool fortran_quat_to_w_mat(
    c_RealArr quat /* 1D_NOT_real */,
    c_RealArr w_mat /* 2D_NOT_real */);
FixedArray2D<Real, 3, 3> quat_to_w_mat(FixedArray1D<Real, 4> quat);
extern "C" void fortran_query_string(
    c_Char query_str /* 0D_NOT_character */,
    c_Bool& upcase /* 0D_NOT_logical */,
    c_Char return_str /* 0D_NOT_character */,
    c_Int& ix /* 0D_NOT_integer */,
    c_Int& ios /* 0D_NOT_integer */);
void query_string(
    std::string query_str,
    bool upcase,
    std::string return_str,
    int ix,
    int ios);
extern "C" bool fortran_quote(
    c_Char str /* 0D_NOT_character */,
    c_Char q_str /* 0D_ALLOC_character */);
void quote(std::string str, std::string q_str);

// Skipped unusable routine quoten:
// Variable-sized inout character array: str(:) 1D_ALLOC_character
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine ran_default_state:
// Untranslated type: RandomStateProxy (0D_NOT_type)
// Untranslated type: RandomStateProxy (0D_NOT_type)

// Skipped unusable routine ran_engine:
// Untranslated type: RandomStateProxy (0D_NOT_type)

// Skipped unusable routine ran_gauss_converter:
// Untranslated type: RandomStateProxy (0D_NOT_type)

// Skipped unusable routine ran_gauss_scalar:
// Untranslated type: RandomStateProxy (0D_NOT_type)

// Skipped unusable routine ran_gauss_vector:
// Untranslated type: RandomStateProxy (0D_NOT_type)
extern "C" void fortran_ran_seed_get(c_Int& seed /* 0D_NOT_integer */);
int ran_seed_get();
extern "C" void fortran_ran_seed_put(
    c_Int& seed /* 0D_NOT_integer */,
    c_Int* mpi_offset /* 0D_NOT_integer */);
void ran_seed_put(int seed, std::optional<int> mpi_offset = std::nullopt);

// Skipped unusable routine ran_uniform_scalar:
// Untranslated type: RandomStateProxy (0D_NOT_type)

// Skipped unusable routine ran_uniform_vector:
// Untranslated type: RandomStateProxy (0D_NOT_type)
extern "C" bool fortran_real_num_fortran_format(
    c_Real& number /* 0D_NOT_real */,
    c_Int& width /* 0D_NOT_integer */,
    c_Int* n_blanks /* 0D_NOT_integer */,
    c_Char fmt_str /* 0D_NOT_character */);
void real_num_fortran_format(
    double number,
    int width,
    optional_ref<int> n_blanks,
    std::string fmt_str);
extern "C" bool fortran_real_path(
    c_Char path_in /* 0D_NOT_character */,
    c_Char path_out /* 0D_NOT_character */,
    c_Bool& is_ok /* 0D_NOT_logical */);
void real_path(std::string path_in, std::string path_out, bool is_ok);
extern "C" bool fortran_real_str(
    c_Real& r_num /* 0D_NOT_real */,
    c_Int* n_signif /* 0D_NOT_integer */,
    c_Int* n_decimal /* 0D_NOT_integer */,
    c_Char str /* 0D_ALLOC_character */);
void real_str(
    double r_num,
    optional_ref<int> n_signif,
    optional_ref<int> n_decimal,
    std::string str);
extern "C" bool fortran_real_to_string(
    c_Real& real_num /* 0D_NOT_real */,
    c_Int& width /* 0D_NOT_integer */,
    c_Int* n_signif /* 0D_NOT_integer */,
    c_Int* n_decimal /* 0D_NOT_integer */,
    c_Char str /* 0D_NOT_character */);
void real_to_string(
    double real_num,
    int width,
    optional_ref<int> n_signif,
    optional_ref<int> n_decimal,
    std::string str);
extern "C" void fortran_reallocate_spline(
    void* spline /* 1D_ALLOC_type */,
    c_Int& n /* 0D_NOT_integer */,
    c_Int* n_min /* 0D_NOT_integer */,
    c_Bool* exact /* 0D_NOT_logical */);
void reallocate_spline(
    SplineProxyAlloc1D& spline,
    int n,
    std::optional<int> n_min = std::nullopt,
    std::optional<bool> exact = std::nullopt);

// Skipped unusable routine reals_to_string:
// No matching docstring

// Skipped unusable routine reals_to_table_row:
// No matching docstring
extern "C" bool fortran_rms_value(
    void* val_arr /* 1D_ALLOC_real */,
    void* good_val /* 1D_ALLOC_logical */,
    c_Real& ave_val /* 0D_NOT_real */,
    c_Real& rms_val /* 0D_NOT_real */);
double rms_value(
    RealAlloc1D& val_arr,
    optional_ref<BoolAlloc1D> good_val,
    double rms_val);
extern "C" bool fortran_rot_2d(
    c_RealArr vec_in /* 1D_NOT_real */,
    c_Real& angle /* 0D_NOT_real */,
    c_RealArr vec_out /* 1D_NOT_real */);
void rot_2d(
    FixedArray1D<Real, 2> vec_in,
    double angle,
    FixedArray1D<Real, 2> vec_out);

// Skipped unusable routine rotate_mat:
// Variable inout sized array: mat(:,:) 2D_NOT_real
extern "C" void fortran_rotate_vec(
    void* vec /* 1D_ALLOC_real */,
    c_Int& axis /* 0D_NOT_integer */,
    c_Real& angle /* 0D_NOT_real */);
void rotate_vec(RealAlloc1D& vec, int axis, double angle);
extern "C" bool fortran_rotate_vec_given_axis_angle(
    c_RealArr vec_in /* 1D_NOT_real */,
    void* axis /* 1D_ALLOC_real */,
    c_Real& angle /* 0D_NOT_real */,
    c_RealArr vec_out /* 1D_NOT_real */);
FixedArray1D<Real, 3> rotate_vec_given_axis_angle(
    FixedArray1D<Real, 3> vec_in,
    RealAlloc1D& axis,
    double angle);
extern "C" bool fortran_rp8(
    c_Int& int_in /* 0D_NOT_integer */,
    c_Real& re_out /* 0D_NOT_real */);
double rp8(int int_in);
extern "C" void fortran_run_timer(
    c_Char command /* 0D_NOT_character */,
    c_Real* time /* 0D_NOT_real */,
    c_Real* time0 /* 0D_NOT_real */);
void run_timer(
    std::string command,
    optional_ref<double> time = std::nullopt,
    optional_ref<double> time0 = std::nullopt);

// Skipped unusable routine set_all_ptr:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" bool fortran_set_parameter_int(
    c_Int& param_val /* 0D_NOT_integer */,
    c_Int& set_val /* 0D_NOT_integer */,
    c_Int& save_val /* 0D_NOT_integer */);
void set_parameter_int(int param_val, int set_val, int save_val);
extern "C" bool fortran_set_parameter_logic(
    c_Bool& param_val /* 0D_NOT_logical */,
    c_Bool& set_val /* 0D_NOT_logical */,
    c_Bool& save_val /* 0D_NOT_logical */);
void set_parameter_logic(bool param_val, bool set_val, bool save_val);
extern "C" bool fortran_set_parameter_real(
    c_Real& param_val /* 0D_NOT_real */,
    c_Real& set_val /* 0D_NOT_real */,
    c_Real& save_val /* 0D_NOT_real */);
void set_parameter_real(double param_val, double set_val, double save_val);
extern "C" bool fortran_set_species_charge(
    c_Int& species_in /* 0D_NOT_integer */,
    c_Int& charge /* 0D_NOT_integer */,
    c_Int& species_charged /* 0D_NOT_integer */);
int set_species_charge(int species_in, int charge);
extern "C" bool fortran_sinc(
    c_Real& x /* 0D_NOT_real */,
    c_Int* nd /* 0D_NOT_integer */,
    c_Real& y /* 0D_NOT_real */);
void sinc(double x, optional_ref<int> nd, double y);
extern "C" bool fortran_sincc(
    c_Real& x /* 0D_NOT_real */,
    c_Int* nd /* 0D_NOT_integer */,
    c_Real& y /* 0D_NOT_real */);
void sincc(double x, optional_ref<int> nd, double y);
extern "C" bool fortran_sinhx_x(
    c_Real& x /* 0D_NOT_real */,
    c_Int* nd /* 0D_NOT_integer */,
    c_Real& y /* 0D_NOT_real */);
void sinhx_x(double x, optional_ref<int> nd, double y);
extern "C" void fortran_skip_header(
    c_Int& ix_unit /* 0D_NOT_integer */,
    c_Bool& error_flag /* 0D_NOT_logical */);
void skip_header(int ix_unit, bool error_flag);
extern "C" bool fortran_species_id(
    c_Char name /* 0D_NOT_character */,
    c_Int* default_ /* 0D_NOT_integer */,
    c_Bool* print_err /* 0D_NOT_logical */,
    c_Int& species /* 0D_NOT_integer */);
int species_id(
    std::string name,
    std::optional<int> default_ = std::nullopt,
    std::optional<bool> print_err = std::nullopt);
extern "C" bool fortran_species_id_from_openpmd(
    c_Char pmd_name /* 0D_NOT_character */,
    c_Int& charge /* 0D_NOT_integer */,
    c_Int& species /* 0D_NOT_integer */);
int species_id_from_openpmd(std::string pmd_name, int charge);
extern "C" bool fortran_species_name(
    c_Int& species /* 0D_NOT_integer */,
    c_Char name /* 0D_NOT_character */);
std::string species_name(int species);
extern "C" bool fortran_species_of(
    c_Real& mass /* 0D_NOT_real */,
    c_Int& charge /* 0D_NOT_integer */,
    c_Int& species /* 0D_NOT_integer */);
int species_of(double mass, int charge);
extern "C" bool fortran_spin_of(
    c_Int& species /* 0D_NOT_integer */,
    c_Real* non_subatomic_default /* 0D_NOT_real */,
    c_Real& spin /* 0D_NOT_real */);
double spin_of(
    int species,
    std::optional<double> non_subatomic_default = std::nullopt);
extern "C" bool fortran_spline1(
    void* a_spline /* 0D_NOT_type */,
    c_Real& x /* 0D_NOT_real */,
    c_Int* n /* 0D_NOT_integer */,
    c_Real& y /* 0D_NOT_real */);
double spline1(
    SplineProxy& a_spline,
    double x,
    std::optional<int> n = std::nullopt);
extern "C" void fortran_spline_akima(
    void* spline /* 1D_ALLOC_type */,
    c_Bool& ok /* 0D_NOT_logical */);
bool spline_akima(SplineProxyAlloc1D& spline);
extern "C" void fortran_spline_akima_interpolate(
    void* x_knot /* 1D_ALLOC_real */,
    void* y_knot /* 1D_ALLOC_real */,
    c_Real& x /* 0D_NOT_real */,
    c_Bool& ok /* 0D_NOT_logical */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */);
struct SplineAkimaInterpolate {
  bool ok;
  double y;
  double dy;
};
SplineAkimaInterpolate spline_akima_interpolate(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    double x);
extern "C" void fortran_spline_evaluate(
    void* spline /* 1D_ALLOC_type */,
    c_Real& x /* 0D_NOT_real */,
    c_Bool& ok /* 0D_NOT_logical */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */);
struct SplineEvaluate {
  bool ok;
  double y;
  double dy;
};
SplineEvaluate spline_evaluate(SplineProxyAlloc1D& spline, double x);
extern "C" bool fortran_sqrt_alpha(
    c_Real& alpha /* 0D_NOT_real */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */);
void sqrt_alpha(double alpha, double x, double y);
extern "C" bool fortran_sqrt_one(
    c_Real& x /* 0D_NOT_real */,
    c_Int* nd /* 0D_NOT_integer */,
    c_Real& ds1 /* 0D_NOT_real */);
void sqrt_one(double x, optional_ref<int> nd, double ds1);
extern "C" bool fortran_str_count(
    c_Char str /* 0D_NOT_character */,
    c_Char match /* 0D_NOT_character */,
    c_Int& num /* 0D_NOT_integer */);
void str_count(std::string str, std::string match, int num);
extern "C" void fortran_str_downcase(
    c_Char dst /* 0D_NOT_character */,
    c_Char src /* 0D_NOT_character */);
void str_downcase(std::string dst, std::string src);
extern "C" bool fortran_str_first_in_set(
    c_Char line /* 0D_NOT_character */,
    c_Char set /* 0D_NOT_character */,
    c_Bool* ignore_clauses /* 0D_NOT_logical */,
    c_Int& ix_match /* 0D_NOT_integer */);
void str_first_in_set(
    std::string line,
    std::string set,
    optional_ref<bool> ignore_clauses,
    int ix_match);
extern "C" bool fortran_str_first_not_in_set(
    c_Char line /* 0D_NOT_character */,
    c_Char set /* 0D_NOT_character */,
    c_Int& ix_match /* 0D_NOT_integer */);
void str_first_not_in_set(std::string line, std::string set, int ix_match);
extern "C" bool fortran_str_last_in_set(
    c_Char line /* 0D_NOT_character */,
    c_Char set /* 0D_NOT_character */,
    c_Int& ix_match /* 0D_NOT_integer */);
void str_last_in_set(std::string line, std::string set, int ix_match);
extern "C" bool fortran_str_last_not_in_set(
    c_Char line /* 0D_NOT_character */,
    c_Char set /* 0D_NOT_character */,
    c_Int& ix_match /* 0D_NOT_integer */);
void str_last_not_in_set(std::string line, std::string set, int ix_match);
extern "C" bool fortran_str_match_wild(
    c_Char str /* 0D_NOT_character */,
    c_Char pat /* 0D_NOT_character */,
    c_Bool& a_match /* 0D_NOT_logical */);
void str_match_wild(std::string str, std::string pat, bool a_match);

// Skipped unusable routine str_set:
// Routine in configuration skip list
extern "C" void fortran_str_substitute(
    c_Char string /* 0D_NOT_character */,
    c_Char str_match /* 0D_NOT_character */,
    c_Char str_replace /* 0D_NOT_character */,
    c_Bool* do_trim /* 0D_NOT_logical */,
    c_Bool* ignore_escaped /* 0D_NOT_logical */);
void str_substitute(
    std::string string,
    optional_ref<std::string> str_match = std::nullopt,
    optional_ref<std::string> str_replace = std::nullopt,
    optional_ref<bool> do_trim = std::nullopt,
    optional_ref<bool> ignore_escaped = std::nullopt);
extern "C" void fortran_str_upcase(
    c_Char dst /* 0D_NOT_character */,
    c_Char src /* 0D_NOT_character */);
void str_upcase(std::string dst, std::string src);
extern "C" bool fortran_string_to_int(
    c_Char line /* 0D_NOT_character */,
    c_Int& default_ /* 0D_NOT_integer */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Int& value /* 0D_NOT_integer */);
void string_to_int(
    std::string line,
    int default_,
    bool err_flag,
    optional_ref<bool> err_print_flag,
    int value);
extern "C" bool fortran_string_to_real(
    c_Char line /* 0D_NOT_character */,
    c_Real& default_ /* 0D_NOT_real */,
    c_Bool& err_flag /* 0D_NOT_logical */,
    c_Bool* err_print_flag /* 0D_NOT_logical */,
    c_Real& value /* 0D_NOT_real */);
void string_to_real(
    std::string line,
    double default_,
    bool err_flag,
    optional_ref<bool> err_print_flag,
    double value);
extern "C" void fortran_string_trim(
    c_Char in_string /* 0D_NOT_character */,
    c_Char out_string /* 0D_NOT_character */,
    c_Int& word_len /* 0D_NOT_integer */);
void string_trim(std::string in_string, std::string out_string, int word_len);
extern "C" void fortran_string_trim2(
    c_Char in_str /* 0D_NOT_character */,
    c_Char delimitors /* 0D_NOT_character */,
    c_Char out_str /* 0D_NOT_character */,
    c_Int& ix_word /* 0D_NOT_integer */,
    c_Char delim /* 0D_NOT_character */,
    c_Int& ix_next /* 0D_NOT_integer */);
void string_trim2(
    std::string in_str,
    std::string delimitors,
    std::string out_str,
    int ix_word,
    std::string delim,
    int ix_next);

// Skipped unusable routine substr:
// Routine in configuration skip list
extern "C" void fortran_super_bicubic_coef(
    c_RealArr y /* 1D_NOT_real */,
    c_RealArr y1 /* 1D_NOT_real */,
    c_RealArr y2 /* 1D_NOT_real */,
    c_RealArr y12 /* 1D_NOT_real */,
    c_Real& d1 /* 0D_NOT_real */,
    c_Real& d2 /* 0D_NOT_real */,
    c_RealArr c /* 2D_NOT_real */);
FixedArray2D<Real, 4, 4> super_bicubic_coef(
    FixedArray1D<Real, 4> y,
    FixedArray1D<Real, 4> y1,
    FixedArray1D<Real, 4> y2,
    FixedArray1D<Real, 4> y12,
    double d1,
    double d2);
extern "C" void fortran_super_bicubic_interpolation(
    c_RealArr y /* 1D_NOT_real */,
    c_RealArr y1 /* 1D_NOT_real */,
    c_RealArr y2 /* 1D_NOT_real */,
    c_RealArr y12 /* 1D_NOT_real */,
    c_Real& x1l /* 0D_NOT_real */,
    c_Real& x1u /* 0D_NOT_real */,
    c_Real& x2l /* 0D_NOT_real */,
    c_Real& x2u /* 0D_NOT_real */,
    c_Real& x1 /* 0D_NOT_real */,
    c_Real& x2 /* 0D_NOT_real */,
    c_Real& ansy /* 0D_NOT_real */,
    c_Real& ansy1 /* 0D_NOT_real */,
    c_Real& ansy2 /* 0D_NOT_real */);
struct SuperBicubicInterpolation {
  double ansy;
  double ansy1;
  double ansy2;
};
SuperBicubicInterpolation super_bicubic_interpolation(
    FixedArray1D<Real, 4> y,
    FixedArray1D<Real, 4> y1,
    FixedArray1D<Real, 4> y2,
    FixedArray1D<Real, 4> y12,
    double x1l,
    double x1u,
    double x2l,
    double x2u,
    double x1,
    double x2);

// Skipped unusable routine super_bracket_root:
// Argument not defined: func (have: [])
// Argument not defined: x1 (have: [])
// Argument not defined: x2 (have: [])
// Argument not defined: status (have: [])
// Argument not defined: x_range (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_brent:
// Argument not defined: ax (have: [])
// Argument not defined: bx (have: [])
// Argument not defined: cx (have: [])
// Argument not defined: func (have: [])
// Argument not defined: rel_tol (have: [])
// Argument not defined: abs_tol (have: [])
// Argument not defined: xmin (have: [])
// Argument not defined: status (have: [])
// Argument not defined: ymin (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_dbrent:
// Argument not defined: ax (have: [])
// Argument not defined: bx (have: [])
// Argument not defined: cx (have: [])
// Argument not defined: func (have: [])
// Argument not defined: dfunc (have: [])
// Argument not defined: rel_tol (have: [])
// Argument not defined: abs_tol (have: [])
// Argument not defined: xmin (have: [])
// Argument not defined: status (have: [])
// Argument not defined: func_min (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_gaussj:
// Variable inout sized array: a(:,:) 2D_NOT_real
// Variable inout sized array: b(:,:) 2D_NOT_real

// Skipped unusable routine super_ludcmp:
// Variable in sized array: a(:,:) 2D_NOT_real

// Skipped unusable routine super_mnbrak:
// Argument not defined: ax (have: [])
// Argument not defined: bx (have: [])
// Argument not defined: cx (have: [])
// Argument not defined: fa (have: [])
// Argument not defined: fb (have: [])
// Argument not defined: fc (have: [])
// Argument not defined: func (have: [])
// Argument not defined: status (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_mrqcof:
// Variable inout sized array: co_alpha(:,:) 2D_NOT_real
// Untranslated type: SuperMrqminStorageProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_mrqmin:
// Argument not defined: y (have: [])
// Argument not defined: weight (have: [])
// Argument not defined: a (have: [])
// Argument not defined: chisq (have: [])
// Argument not defined: funcs (have: [])
// Argument not defined: storage (have: [])
// Argument not defined: alamda (have: [])
// Argument not defined: status (have: [])
// Argument not defined: maska (have: [])
// Argument not defined: print_err (have: [])
// Translated arg count mismatch (unsupported?)
extern "C" void fortran_super_polint(
    void* xa /* 1D_ALLOC_real */,
    void* ya /* 1D_ALLOC_real */,
    c_Real& x /* 0D_NOT_real */,
    c_Real& y /* 0D_NOT_real */,
    c_Real& dy /* 0D_NOT_real */);
struct SuperPolint {
  double y;
  double dy;
};
SuperPolint super_polint(RealAlloc1D& xa, RealAlloc1D& ya, double x);
extern "C" bool fortran_super_poly(
    c_Real& x /* 0D_NOT_real */,
    void* coeffs /* 1D_ALLOC_real */,
    c_Real& value /* 0D_NOT_real */);
double super_poly(double x, RealAlloc1D& coeffs);

// Skipped unusable routine super_qromb:
// Argument not defined: func (have: [])
// Argument not defined: a (have: [])
// Argument not defined: b (have: [])
// Argument not defined: rel_tol (have: [])
// Argument not defined: abs_tol (have: [])
// Argument not defined: k_order (have: [])
// Argument not defined: err_flag (have: [])
// Argument not defined: integral (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_qromb_2d:
// Argument not defined: func (have: [])
// Argument not defined: ax (have: [])
// Argument not defined: bx (have: [])
// Argument not defined: ay (have: [])
// Argument not defined: by (have: [])
// Argument not defined: rel_tol (have: [])
// Argument not defined: abs_tol (have: [])
// Argument not defined: k_order (have: [])
// Argument not defined: err_flag (have: [])
// Argument not defined: integral (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_rtsafe:
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_sobseq:
// Untranslated type: RandomStateProxy (0D_NOT_type)
extern "C" void fortran_super_sort(void* arr /* 1D_ALLOC_integer */);
void super_sort(IntAlloc1D& arr);

// Skipped unusable routine super_trapzd:
// Argument not defined: func (have: [])
// Argument not defined: a (have: [])
// Argument not defined: b (have: [])
// Argument not defined: s (have: [])
// Argument not defined: n (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine super_zbrent:
// Argument not defined: func (have: [])
// Argument not defined: x1 (have: [])
// Argument not defined: x2 (have: [])
// Argument not defined: rel_tol (have: [])
// Argument not defined: abs_tol (have: [])
// Argument not defined: status (have: [])
// Argument not defined: func_val (have: [])
// Argument not defined: x_zero (have: [])
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine svd_fit:
// Variable inout sized array: A(:,:) 2D_NOT_real
// Variable inout sized array: v_mat(:,:) 2D_NOT_real
extern "C" void fortran_system_command(
    c_Char line /* 0D_NOT_character */,
    c_Bool* err_flag /* 0D_NOT_logical */);
void system_command(
    std::string line,
    optional_ref<bool> err_flag = std::nullopt);
extern "C" bool fortran_to_str(
    c_Real& num /* 0D_NOT_real */,
    c_Int* max_signif /* 0D_NOT_integer */,
    c_Char string /* 0D_ALLOC_character */);
void to_str(double num, optional_ref<int> max_signif, std::string string);
extern "C" bool fortran_tricubic_cmplx_eval(
    c_Real& x_norm /* 0D_NOT_real */,
    c_Real& y_norm /* 0D_NOT_real */,
    c_Real& z_norm /* 0D_NOT_real */,
    void* tri_coef /* 0D_NOT_type */,
    c_Complex& df_dx /* 0D_NOT_complex */,
    c_Complex& df_dy /* 0D_NOT_complex */,
    c_Complex& df_dz /* 0D_NOT_complex */,
    c_Complex& f_val /* 0D_NOT_complex */);
struct TricubicCmplxEval {
  std::complex<double> df_dx;
  std::complex<double> df_dy;
  std::complex<double> df_dz;
  std::complex<double> f_val;
};
TricubicCmplxEval tricubic_cmplx_eval(
    double x_norm,
    double y_norm,
    double z_norm,
    TricubicCmplxCoefProxy& tri_coef);

// Skipped unusable routine tricubic_compute_cmplx_field_at_3d_box:
// Untranslated type: CmplxFieldAt3dBoxProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tricubic_compute_field_at_3d_box:
// Untranslated type: FieldAt3dBoxProxy (0D_NOT_type)
// Translated arg count mismatch (unsupported?)

// Skipped unusable routine tricubic_eval:
// Untranslated type: TricubicCoefProxy (0D_NOT_type)

// Skipped unusable routine tricubic_interpolation_cmplx_coefs:
// Untranslated type: CmplxFieldAt3dBoxProxy (0D_NOT_type)

// Skipped unusable routine tricubic_interpolation_coefs:
// Untranslated type: FieldAt3dBoxProxy (0D_NOT_type)
// Untranslated type: TricubicCoefProxy (0D_NOT_type)
extern "C" void fortran_type_this_file(c_Char filename /* 0D_NOT_character */);
void type_this_file(std::string filename);

// Skipped unusable routine unquote:
// No matching docstring

// Skipped unusable routine upcase:
// No matching docstring
extern "C" void fortran_upcase_string(c_Char string /* 0D_NOT_character */);
void upcase_string(std::string string);

// Skipped unusable routine value_of_all_ptr:
// Untranslated type: AllPointerProxy (0D_NOT_type)
extern "C" bool fortran_virtual_memory_usage(c_Int& usage /* 0D_NOT_integer */);
int virtual_memory_usage();
extern "C" void fortran_w_mat_to_axis_angle(
    c_RealArr w_mat /* 2D_NOT_real */,
    c_RealArr axis /* 1D_NOT_real */,
    c_Real& angle /* 0D_NOT_real */);
struct WMatToAxisAngle {
  FixedArray1D<Real, 3> axis;
  double angle;
};
WMatToAxisAngle w_mat_to_axis_angle(FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_w_mat_to_quat(
    c_RealArr w_mat /* 2D_NOT_real */,
    c_RealArr quat /* 1D_NOT_real */);
FixedArray1D<Real, 4> w_mat_to_quat(FixedArray2D<Real, 3, 3> w_mat);
extern "C" bool fortran_word_len(
    c_Char wording /* 0D_NOT_character */,
    c_Int& wlen /* 0D_NOT_integer */);
void word_len(std::string wording, int wlen);
extern "C" void fortran_word_read(
    c_Char in_str /* 0D_NOT_character */,
    c_Char delim_list /* 0D_NOT_character */,
    c_Char word /* 0D_NOT_character */,
    c_Int& ix_word /* 0D_NOT_integer */,
    c_Char delim /* 0D_NOT_character */,
    c_Bool& delim_found /* 0D_NOT_logical */,
    c_Char out_str /* 0D_NOT_character */,
    c_Bool* ignore_interior /* 0D_NOT_logical */);
void word_read(
    std::string in_str,
    std::string delim_list,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    std::string out_str,
    optional_ref<bool> ignore_interior = std::nullopt);
extern "C" bool fortran_x0_radiation_length(
    c_Int& species /* 0D_NOT_integer */,
    c_Real& x0 /* 0D_NOT_real */);
double x0_radiation_length(int species);
} // namespace SimUtils
