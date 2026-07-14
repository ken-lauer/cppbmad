module cppbmad_sim_utils_routines

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

use sim_utils_interface, only: all_pointer_to_string, asinc, assert_equal, calc_file_number, &
    change_file_number, complex_error_function, cos_one, cosc, cross_product, &
    date_and_time_stamp, detab, display_size_and_resolution, dj_bessel, djb_hash, djb_str_hash, &
    downcase_string, err_exit, factorial, faddeeva_function, fff_sub, fft_1d, &
    file_directorizer, file_get, file_get_open, file_suffixer, find_location, &
    gen_complete_elliptic, get_file_number, get_file_time_stamp, get_next_number, i_bessel, &
    i_bessel_extended, increment_file_number, index_nocase, int_str, inverse_prob, &
    is_alphabetic, is_decreasing_sequence, is_increasing_sequence, is_integer, is_logical, &
    is_real, j_bessel, linear_fit, linear_fit_2d, logic_str, lunget, make_legal_comment, &
    match_reg, match_wild, match_word, milli_sleep, n_choose_k, n_spline_create, nametable_add, &
    nametable_bracket_indexx, nametable_change1, nametable_init, nametable_remove, ordinal_str, &
    parse_fortran_format, pointer_to_locations, poly_eval, probability_funct, quadratic_roots, &
    query_string, quote, quoten, real_num_fortran_format, real_path, real_str, real_to_string, &
    rms_value, rot_2d, run_timer, set_all_ptr, set_parameter, sinc, sincc, sinhx_x, &
    skip_header, sqrt_alpha, sqrt_one, str_count, str_downcase, str_first_in_set, &
    str_first_not_in_set, str_last_in_set, str_last_not_in_set, str_match_wild, str_substitute, &
    str_upcase, string_to_int, string_to_real, string_trim, string_trim2, system_command, &
    test_tune_tracker_lock, to_str, type_this_file, upcase_string, value_of_all_ptr, &
    virtual_memory_usage, word_len, word_read

use random_mod, only: allocate_thread_states, pointer_to_ran_state, ran_default_state, &
    ran_engine, ran_gauss_converter, ran_gauss_scalar, ran_gauss_vector, ran_seed_get, &
    ran_seed_put, ran_uniform, super_sobseq

use particle_species_mod, only: anomalous_moment_of, antiparticle, atomic_number, &
    atomic_species_id, charge_of, charge_to_mass_of, is_subatomic_species, mass_of, &
    molecular_components, openpmd_species_name, set_species_charge, species_id, &
    species_id_from_openpmd, species_name, species_of, spin_of, x0_radiation_length

use all_phase_fft, only: apfft, apfft_corr, apfft_ext, hanhan

use rotation_3d_mod, only: axis_angle_to_quat, axis_angle_to_w_mat, omega_to_quat, quat_conj, &
    quat_inverse, quat_mul, quat_rotate, quat_to_axis_angle, quat_to_omega, quat_to_w_mat, &
    rotate_vec, rotate_vec_given_axis_angle, w_mat_to_axis_angle, w_mat_to_quat

use cubic_interpolation_mod, only: bicubic_cmplx_eval, bicubic_eval, &
    bicubic_interpolation_cmplx_coefs, bicubic_interpolation_coefs, tricubic_cmplx_eval, &
    tricubic_eval, tricubic_interpolation_cmplx_coefs, tricubic_interpolation_coefs

use bin_mod, only: bin_2d, bin_data, bin_data_density, bin_data_density_2d, bin_index, &
    bin_x_center, count_at_index, general_bin_count, general_bin_index, &
    general_bin_index_in_bounds, n_bins_automatic

use bit_mod, only: bit_set

use spline_mod, only: bracket_index_for_spline, create_a_spline, end_akima_spline_calc, &
    reallocate_spline, spline1, spline_akima, spline_akima_interpolate, spline_evaluate

use elliptic_integral_mod, only: celbd, elbd, elcbd, ellipinc, elsbd, gelbd, rcelbd, relbd, &
    relcbd, relsbd, rgelbd, rserbd, serbd, test_xgelbd

use command_line_mod, only: cesr_getarg, cesr_iargc

use fourier_mod, only: coarse_frequency_estimate, fine_frequency_estimate, fourier_amplitude

use windowls_mod, only: destfixedwindowls, fixedwindowls, initfixedwindowls

use input_mod, only: get_a_char, get_tty_char, read_a_line, readline_read_history, &
    readline_write_history

use lmdif_mod, only: initial_lmdif, suggest_lmdif

use naff_mod, only: interpolated_fft, interpolated_fft_gsl, maximize_projection, naff, projdd

use sim_utils_struct, only: is_false, is_true

use modulo2_mod, only: modulo2_dp, modulo2_int, modulo2_qp, modulo2_sp

use output_mod, only: out_io, out_io_buffer_get_line, out_io_buffer_num_lines, &
    out_io_buffer_reset, out_io_print_and_capture_setup, output_direct

use precision_def, only: rp8

use sign_of_mod, only: sign_of

use super_recipes_mod, only: super_bicubic_coef, super_bicubic_interpolation, super_polint, &
    super_poly, super_sort


use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

subroutine fortran_all_pointer_to_string (a_ptr, err, str) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: all_pointer_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err
  logical, target :: f_err_native
  logical, pointer :: f_err_native_ptr
  logical(c_bool), pointer :: f_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), value :: a_ptr  ! 0D_NOT_type
  type(all_pointer_struct), pointer :: f_a_ptr
  ! ** End of parameters **
  ! inout: f_a_ptr 0D_NOT_type
  if (.not. c_associated(a_ptr)) return
  call c_f_pointer(a_ptr, f_a_ptr)
  ! in: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_native = f_err_ptr
    f_err_native_ptr => f_err_native
  else
    f_err_native_ptr => null()
  endif
  f_str = all_pointer_to_string(f_a_ptr, f_err_native_ptr)

  ! out: f_str 0D_NOT_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1])
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_allocate_thread_states () bind(c)

  use array_desc_mod
  implicit none
  ! ** End of parameters **
  call allocate_thread_states()

end subroutine
subroutine fortran_anomalous_moment_of (species, moment) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: moment  ! 0D_NOT_real
  real(rp) :: f_moment
  real(c_double), pointer :: f_moment_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_moment = anomalous_moment_of(f_species)

  ! out: f_moment 0D_NOT_real
  call c_f_pointer(moment, f_moment_ptr)
  f_moment_ptr = f_moment
end subroutine
subroutine fortran_antiparticle (species, anti_species) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: anti_species  ! 0D_NOT_integer
  integer :: f_anti_species
  integer(c_int), pointer :: f_anti_species_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_anti_species = antiparticle(f_species)

  ! out: f_anti_species 0D_NOT_integer
  call c_f_pointer(anti_species, f_anti_species_ptr)
  f_anti_species_ptr = f_anti_species
end subroutine
subroutine fortran_apfft (rdata_in, bounds, window, phase, diag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: window
  character(len=4096), target :: f_window
  character(kind=c_char), pointer :: f_window_ptr(:)
  real(c_double) :: phase  ! 0D_NOT_real
  real(rp) :: f_phase
  type(c_ptr), intent(in), value :: diag  ! 0D_NOT_integer
  integer(c_int) :: f_diag
  integer(c_int), pointer :: f_diag_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: rdata_in
  real(rp), pointer :: f_rdata_in(:)
  real(c_double), pointer :: f_rdata_in_ptr(:)
  type(array_descriptor_t), intent(in) :: bounds
  real(rp) :: f_bounds(2)
  real(c_double), pointer :: f_bounds_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(rdata_in%data_ptr)) then
    call c_f_pointer(rdata_in%data_ptr, f_rdata_in_ptr, [rdata_in%dims(1)])
    f_rdata_in => f_rdata_in_ptr
  else
    f_rdata_in => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(bounds%data_ptr)) then
    call c_f_pointer(bounds%data_ptr, f_bounds_ptr, [bounds%dims(1)])
    f_bounds = f_bounds_ptr(:)
  else
    f_bounds_ptr => null()
  endif
  ! in: f_window 0D_NOT_character
  if (.not. c_associated(window)) return
  call c_f_pointer(window, f_window_ptr, [huge(0)])
  call to_f_str(f_window_ptr, f_window)
  ! in: f_phase 0D_NOT_real
  f_phase = phase
  ! in: f_diag 0D_NOT_integer
  if (c_associated(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
  else
    f_diag_ptr => null()
  endif
  call apfft(f_rdata_in, f_bounds, f_window, f_phase, f_diag_ptr)

end subroutine
subroutine fortran_apfft_corr (rdata_in, bounds, window, phase, amp, freq, diag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: rdata_in
  real(rp), pointer :: f_rdata_in(:)
  real(c_double), pointer :: f_rdata_in_ptr(:)
  type(array_descriptor_t), intent(in) :: bounds
  real(rp) :: f_bounds(2)
  real(c_double), pointer :: f_bounds_ptr(:)
  type(c_ptr), intent(in), value :: window
  character(len=4096), target :: f_window
  character(kind=c_char), pointer :: f_window_ptr(:)
  type(c_ptr), intent(in), value :: diag  ! 0D_NOT_integer
  integer(c_int) :: f_diag
  integer(c_int), pointer :: f_diag_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: phase  ! 0D_NOT_real
  real(rp) :: f_phase
  real(c_double), pointer :: f_phase_ptr
  type(c_ptr), intent(in), value :: amp  ! 0D_NOT_real
  real(rp) :: f_amp
  real(c_double), pointer :: f_amp_ptr
  type(c_ptr), intent(in), value :: freq  ! 0D_NOT_real
  real(rp) :: f_freq
  real(c_double), pointer :: f_freq_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(rdata_in%data_ptr)) then
    call c_f_pointer(rdata_in%data_ptr, f_rdata_in_ptr, [rdata_in%dims(1)])
    f_rdata_in => f_rdata_in_ptr
  else
    f_rdata_in => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(bounds%data_ptr)) then
    call c_f_pointer(bounds%data_ptr, f_bounds_ptr, [bounds%dims(1)])
    f_bounds = f_bounds_ptr(:)
  else
    f_bounds_ptr => null()
  endif
  ! in: f_window 0D_NOT_character
  if (.not. c_associated(window)) return
  call c_f_pointer(window, f_window_ptr, [huge(0)])
  call to_f_str(f_window_ptr, f_window)
  ! in: f_diag 0D_NOT_integer
  if (c_associated(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
  else
    f_diag_ptr => null()
  endif
  call apfft_corr(f_rdata_in, f_bounds, f_window, f_phase, f_amp, f_freq, f_diag_ptr)

  ! out: f_phase 0D_NOT_real
  call c_f_pointer(phase, f_phase_ptr)
  f_phase_ptr = f_phase
  ! out: f_amp 0D_NOT_real
  call c_f_pointer(amp, f_amp_ptr)
  f_amp_ptr = f_amp
  ! out: f_freq 0D_NOT_real
  call c_f_pointer(freq, f_freq_ptr)
  f_freq_ptr = f_freq
end subroutine
subroutine fortran_apfft_ext (rdata, bounds, window, phase, amp, freq, diag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: window
  character(len=4096), target :: f_window
  character(kind=c_char), pointer :: f_window_ptr(:)
  real(c_double) :: phase  ! 0D_NOT_real
  real(rp) :: f_phase
  real(c_double) :: amp  ! 0D_NOT_real
  real(rp) :: f_amp
  real(c_double) :: freq  ! 0D_NOT_real
  real(rp) :: f_freq
  type(c_ptr), intent(in), value :: diag  ! 0D_NOT_integer
  integer(c_int) :: f_diag
  integer(c_int), pointer :: f_diag_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: rdata
  real(rp), pointer :: f_rdata(:)
  real(c_double), pointer :: f_rdata_ptr(:)
  type(array_descriptor_t), intent(in) :: bounds
  real(rp) :: f_bounds(2)
  real(c_double), pointer :: f_bounds_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(rdata%data_ptr)) then
    call c_f_pointer(rdata%data_ptr, f_rdata_ptr, [rdata%dims(1)])
    f_rdata => f_rdata_ptr
  else
    f_rdata => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(bounds%data_ptr)) then
    call c_f_pointer(bounds%data_ptr, f_bounds_ptr, [bounds%dims(1)])
    f_bounds = f_bounds_ptr(:)
  else
    f_bounds_ptr => null()
  endif
  ! in: f_window 0D_NOT_character
  if (.not. c_associated(window)) return
  call c_f_pointer(window, f_window_ptr, [huge(0)])
  call to_f_str(f_window_ptr, f_window)
  ! in: f_phase 0D_NOT_real
  f_phase = phase
  ! in: f_amp 0D_NOT_real
  f_amp = amp
  ! in: f_freq 0D_NOT_real
  f_freq = freq
  ! in: f_diag 0D_NOT_integer
  if (c_associated(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
  else
    f_diag_ptr => null()
  endif
  call apfft_ext(f_rdata, f_bounds, f_window, f_phase, f_amp, f_freq, f_diag_ptr)

end subroutine
subroutine fortran_asinc (x, nd, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_nd 0D_NOT_integer
  if (c_associated(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
  else
    f_nd_ptr => null()
  endif
  f_y = asinc(f_x, f_nd_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_assert_equal (int_arr, err_str, ival) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: int_arr
  integer, pointer :: f_int_arr(:)
  integer(c_int), pointer :: f_int_arr_ptr(:)
  type(c_ptr), intent(in), value :: err_str
  character(len=4096), target :: f_err_str
  character(kind=c_char), pointer :: f_err_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ival  ! 0D_NOT_integer
  integer :: f_ival
  integer(c_int), pointer :: f_ival_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_integer) in
  if (c_associated(int_arr%data_ptr)) then
    call c_f_pointer(int_arr%data_ptr, f_int_arr_ptr, [int_arr%dims(1)])
    f_int_arr => f_int_arr_ptr
  else
    f_int_arr => null()
  endif
  ! in: f_err_str 0D_NOT_character
  if (.not. c_associated(err_str)) return
  call c_f_pointer(err_str, f_err_str_ptr, [huge(0)])
  call to_f_str(f_err_str_ptr, f_err_str)
  f_ival = assert_equal(f_int_arr, f_err_str)

  ! out: f_ival 0D_NOT_integer
  call c_f_pointer(ival, f_ival_ptr)
  f_ival_ptr = f_ival
end subroutine
subroutine fortran_atomic_number (species, atomic_num) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: atomic_num  ! 0D_NOT_integer
  integer :: f_atomic_num
  integer(c_int), pointer :: f_atomic_num_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_atomic_num = atomic_number(f_species)

  ! out: f_atomic_num 0D_NOT_integer
  call c_f_pointer(atomic_num, f_atomic_num_ptr)
  f_atomic_num_ptr = f_atomic_num
end subroutine
subroutine fortran_atomic_species_id (charge, is_anti, atomic_num, n_nuc, species_id) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: charge  ! 0D_NOT_integer
  integer :: f_charge
  logical(c_bool) :: is_anti  ! 0D_NOT_logical
  logical :: f_is_anti
  integer(c_int) :: atomic_num  ! 0D_NOT_integer
  integer :: f_atomic_num
  integer(c_int) :: n_nuc  ! 0D_NOT_integer
  integer :: f_n_nuc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species_id  ! 0D_NOT_integer
  integer :: f_species_id
  integer(c_int), pointer :: f_species_id_ptr
  ! ** End of parameters **
  ! in: f_charge 0D_NOT_integer
  f_charge = charge
  ! in: f_is_anti 0D_NOT_logical
  f_is_anti = is_anti
  ! in: f_atomic_num 0D_NOT_integer
  f_atomic_num = atomic_num
  ! in: f_n_nuc 0D_NOT_integer
  f_n_nuc = n_nuc
  f_species_id = atomic_species_id(f_charge, f_is_anti, f_atomic_num, f_n_nuc)

  ! out: f_species_id 0D_NOT_integer
  call c_f_pointer(species_id, f_species_id_ptr)
  f_species_id_ptr = f_species_id
end subroutine
subroutine fortran_axis_angle_to_quat (axis, angle, quat) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    f_axis = f_axis_ptr(:)
  else
    f_axis_ptr => null()
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  !! general array (1D_NOT_real) out
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    ! output-only
  else
    f_quat_ptr => null()
  endif
  f_quat = axis_angle_to_quat(f_axis, f_angle)

  ! out: f_quat 1D_NOT_real
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat_ptr = f_quat(:)
  endif
end subroutine
subroutine fortran_axis_angle_to_w_mat (axis, angle, w_mat) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    f_axis = f_axis_ptr(:)
  else
    f_axis_ptr => null()
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  !! general array (2D_NOT_real) out
  if (c_associated(w_mat%data_ptr)) then
    call c_f_pointer(w_mat%data_ptr, f_w_mat_ptr, [product(w_mat%dims(1:w_mat%rank))])
    ! output-only
  else
    f_w_mat_ptr => null()
  endif
  call axis_angle_to_w_mat(f_axis, f_angle, f_w_mat)

  ! out: f_w_mat 2D_NOT_real
  if (c_associated(w_mat%data_ptr)) f_w_mat_ptr = mat2vec(f_w_mat, product(w_mat%dims(1:w_mat%rank)))
end subroutine
subroutine fortran_bicubic_cmplx_eval (x_norm, y_norm, bi_coef, df_dx, df_dy, f_val) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: bicubic_cmplx_coef_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: x_norm  ! 0D_NOT_real
  real(rp) :: f_x_norm
  real(c_double) :: y_norm  ! 0D_NOT_real
  real(rp) :: f_y_norm
  type(c_ptr), value :: bi_coef  ! 0D_NOT_type
  type(bicubic_cmplx_coef_struct), pointer :: f_bi_coef
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: df_dx  ! 0D_NOT_complex
  complex(rp) :: f_df_dx
  complex(c_double_complex), pointer :: f_df_dx_ptr
  type(c_ptr), intent(in), value :: df_dy  ! 0D_NOT_complex
  complex(rp) :: f_df_dy
  complex(c_double_complex), pointer :: f_df_dy_ptr
  type(c_ptr), intent(in), value :: f_val  ! 0D_NOT_complex
  complex(rp) :: f_f_val
  complex(c_double_complex), pointer :: f_f_val_ptr
  ! ** End of parameters **
  ! in: f_x_norm 0D_NOT_real
  f_x_norm = x_norm
  ! in: f_y_norm 0D_NOT_real
  f_y_norm = y_norm
  ! in: f_bi_coef 0D_NOT_type
  if (.not. c_associated(bi_coef)) return
  call c_f_pointer(bi_coef, f_bi_coef)
  ! out: f_df_dx 0D_NOT_complex
  if (c_associated(df_dx)) then
    call c_f_pointer(df_dx, f_df_dx_ptr)
  else
    f_df_dx_ptr => null()
  endif
  ! out: f_df_dy 0D_NOT_complex
  if (c_associated(df_dy)) then
    call c_f_pointer(df_dy, f_df_dy_ptr)
  else
    f_df_dy_ptr => null()
  endif
  f_f_val = bicubic_cmplx_eval(f_x_norm, f_y_norm, f_bi_coef, f_df_dx, f_df_dy)

  ! out: f_df_dx 0D_NOT_complex
  call c_f_pointer(df_dx, f_df_dx_ptr)
  f_df_dx_ptr = f_df_dx
  ! out: f_df_dy 0D_NOT_complex
  call c_f_pointer(df_dy, f_df_dy_ptr)
  f_df_dy_ptr = f_df_dy
  ! out: f_f_val 0D_NOT_complex
  call c_f_pointer(f_val, f_f_val_ptr)
  f_f_val_ptr = f_f_val
end subroutine
subroutine fortran_bicubic_eval (x_norm, y_norm, bi_coef, df_dx, df_dy, f_val) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: bicubic_coef_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: x_norm  ! 0D_NOT_real
  real(rp) :: f_x_norm
  real(c_double) :: y_norm  ! 0D_NOT_real
  real(rp) :: f_y_norm
  type(c_ptr), value :: bi_coef  ! 0D_NOT_type
  type(bicubic_coef_struct), pointer :: f_bi_coef
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: df_dx  ! 0D_NOT_real
  real(rp) :: f_df_dx
  real(c_double), pointer :: f_df_dx_ptr
  type(c_ptr), intent(in), value :: df_dy  ! 0D_NOT_real
  real(rp) :: f_df_dy
  real(c_double), pointer :: f_df_dy_ptr
  type(c_ptr), intent(in), value :: f_val  ! 0D_NOT_real
  real(rp) :: f_f_val
  real(c_double), pointer :: f_f_val_ptr
  ! ** End of parameters **
  ! in: f_x_norm 0D_NOT_real
  f_x_norm = x_norm
  ! in: f_y_norm 0D_NOT_real
  f_y_norm = y_norm
  ! in: f_bi_coef 0D_NOT_type
  if (.not. c_associated(bi_coef)) return
  call c_f_pointer(bi_coef, f_bi_coef)
  ! out: f_df_dx 0D_NOT_real
  if (c_associated(df_dx)) then
    call c_f_pointer(df_dx, f_df_dx_ptr)
  else
    f_df_dx_ptr => null()
  endif
  ! out: f_df_dy 0D_NOT_real
  if (c_associated(df_dy)) then
    call c_f_pointer(df_dy, f_df_dy_ptr)
  else
    f_df_dy_ptr => null()
  endif
  f_f_val = bicubic_eval(f_x_norm, f_y_norm, f_bi_coef, f_df_dx, f_df_dy)

  ! out: f_df_dx 0D_NOT_real
  call c_f_pointer(df_dx, f_df_dx_ptr)
  f_df_dx_ptr = f_df_dx
  ! out: f_df_dy 0D_NOT_real
  call c_f_pointer(df_dy, f_df_dy_ptr)
  f_df_dy_ptr = f_df_dy
  ! out: f_f_val 0D_NOT_real
  call c_f_pointer(f_val, f_f_val_ptr)
  f_f_val_ptr = f_f_val
end subroutine
subroutine fortran_bicubic_interpolation_cmplx_coefs (field_at_box, bi_coef) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: bicubic_cmplx_coef_struct, cmplx_field_at_2D_box_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: field_at_box  ! 0D_NOT_type
  type(cmplx_field_at_2D_box_struct), pointer :: f_field_at_box
  ! ** Out parameters **
  type(c_ptr), value :: bi_coef  ! 0D_NOT_type
  type(bicubic_cmplx_coef_struct), pointer :: f_bi_coef
  ! ** End of parameters **
  ! in: f_field_at_box 0D_NOT_type
  if (.not. c_associated(field_at_box)) return
  call c_f_pointer(field_at_box, f_field_at_box)
  ! out: f_bi_coef 0D_NOT_type
  if (.not. c_associated(bi_coef)) return
  call c_f_pointer(bi_coef, f_bi_coef)
  call bicubic_interpolation_cmplx_coefs(f_field_at_box, f_bi_coef)

  ! out: f_bi_coef 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_bicubic_interpolation_coefs (field_at_box, bi_coef) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: bicubic_coef_struct, field_at_2D_box_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: field_at_box  ! 0D_NOT_type
  type(field_at_2D_box_struct), pointer :: f_field_at_box
  ! ** Out parameters **
  type(c_ptr), value :: bi_coef  ! 0D_NOT_type
  type(bicubic_coef_struct), pointer :: f_bi_coef
  ! ** End of parameters **
  ! in: f_field_at_box 0D_NOT_type
  if (.not. c_associated(field_at_box)) return
  call c_f_pointer(field_at_box, f_field_at_box)
  ! out: f_bi_coef 0D_NOT_type
  if (.not. c_associated(bi_coef)) return
  call c_f_pointer(bi_coef, f_bi_coef)
  call bicubic_interpolation_coefs(f_field_at_box, f_bi_coef)

  ! out: f_bi_coef 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_bin_2d (data1, data2, weight, min1, max1, min2, max2, n_bins1, n_bins2, &
    bin_data) bind(c)

  use array_desc_mod
  use bin_mod, only: general_bin_struct
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: data1
  real(rp), pointer :: f_data1(:)
  real(c_double), pointer :: f_data1_ptr(:)
  type(array_descriptor_t), intent(in) :: data2
  real(rp), pointer :: f_data2(:)
  real(c_double), pointer :: f_data2_ptr(:)
  type(array_descriptor_t), intent(in) :: weight
  real(rp), pointer :: f_weight(:)
  real(c_double), pointer :: f_weight_ptr(:)
  type(c_ptr), intent(in), value :: min1  ! 0D_NOT_real
  real(c_double) :: f_min1
  real(c_double), pointer :: f_min1_ptr
  type(c_ptr), intent(in), value :: max1  ! 0D_NOT_real
  real(c_double) :: f_max1
  real(c_double), pointer :: f_max1_ptr
  type(c_ptr), intent(in), value :: min2  ! 0D_NOT_real
  real(c_double) :: f_min2
  real(c_double), pointer :: f_min2_ptr
  type(c_ptr), intent(in), value :: max2  ! 0D_NOT_real
  real(c_double) :: f_max2
  real(c_double), pointer :: f_max2_ptr
  type(c_ptr), intent(in), value :: n_bins1  ! 0D_NOT_integer
  integer(c_int) :: f_n_bins1
  integer(c_int), pointer :: f_n_bins1_ptr
  type(c_ptr), intent(in), value :: n_bins2  ! 0D_NOT_integer
  integer(c_int) :: f_n_bins2
  integer(c_int), pointer :: f_n_bins2_ptr
  ! ** Out parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(general_bin_struct), pointer :: f_bin_data
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(data1%data_ptr)) then
    call c_f_pointer(data1%data_ptr, f_data1_ptr, [data1%dims(1)])
    f_data1 => f_data1_ptr
  else
    f_data1 => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(data2%data_ptr)) then
    call c_f_pointer(data2%data_ptr, f_data2_ptr, [data2%dims(1)])
    f_data2 => f_data2_ptr
  else
    f_data2 => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(weight%data_ptr)) then
    call c_f_pointer(weight%data_ptr, f_weight_ptr, [weight%dims(1)])
    f_weight => f_weight_ptr
  else
    f_weight => null()
  endif
  ! in: f_min1 0D_NOT_real
  if (c_associated(min1)) then
    call c_f_pointer(min1, f_min1_ptr)
  else
    f_min1_ptr => null()
  endif
  ! in: f_max1 0D_NOT_real
  if (c_associated(max1)) then
    call c_f_pointer(max1, f_max1_ptr)
  else
    f_max1_ptr => null()
  endif
  ! in: f_min2 0D_NOT_real
  if (c_associated(min2)) then
    call c_f_pointer(min2, f_min2_ptr)
  else
    f_min2_ptr => null()
  endif
  ! in: f_max2 0D_NOT_real
  if (c_associated(max2)) then
    call c_f_pointer(max2, f_max2_ptr)
  else
    f_max2_ptr => null()
  endif
  ! in: f_n_bins1 0D_NOT_integer
  if (c_associated(n_bins1)) then
    call c_f_pointer(n_bins1, f_n_bins1_ptr)
  else
    f_n_bins1_ptr => null()
  endif
  ! in: f_n_bins2 0D_NOT_integer
  if (c_associated(n_bins2)) then
    call c_f_pointer(n_bins2, f_n_bins2_ptr)
  else
    f_n_bins2_ptr => null()
  endif
  ! out: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  f_bin_data = bin_2d(f_data1, f_data2, f_weight, f_min1_ptr, f_max1_ptr, f_min2_ptr, &
      f_max2_ptr, f_n_bins1_ptr, f_n_bins2_ptr)

  ! out: f_bin_data 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_bin_data (data, weight, min, max, n_bins, binned_data) bind(c)

  use array_desc_mod
  use bin_mod, only: bin_struct
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: data
  real(rp), pointer :: f_data(:)
  real(c_double), pointer :: f_data_ptr(:)
  type(array_descriptor_t), intent(in) :: weight
  real(rp), pointer :: f_weight(:)
  real(c_double), pointer :: f_weight_ptr(:)
  type(c_ptr), intent(in), value :: min  ! 0D_NOT_real
  real(c_double) :: f_min
  real(c_double), pointer :: f_min_ptr
  type(c_ptr), intent(in), value :: max  ! 0D_NOT_real
  real(c_double) :: f_max
  real(c_double), pointer :: f_max_ptr
  type(c_ptr), intent(in), value :: n_bins  ! 0D_NOT_integer
  integer(c_int) :: f_n_bins
  integer(c_int), pointer :: f_n_bins_ptr
  ! ** Out parameters **
  type(c_ptr), value :: binned_data  ! 0D_NOT_type
  type(bin_struct), pointer :: f_binned_data
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(data%data_ptr)) then
    call c_f_pointer(data%data_ptr, f_data_ptr, [data%dims(1)])
    f_data => f_data_ptr
  else
    f_data => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(weight%data_ptr)) then
    call c_f_pointer(weight%data_ptr, f_weight_ptr, [weight%dims(1)])
    f_weight => f_weight_ptr
  else
    f_weight => null()
  endif
  ! in: f_min 0D_NOT_real
  if (c_associated(min)) then
    call c_f_pointer(min, f_min_ptr)
  else
    f_min_ptr => null()
  endif
  ! in: f_max 0D_NOT_real
  if (c_associated(max)) then
    call c_f_pointer(max, f_max_ptr)
  else
    f_max_ptr => null()
  endif
  ! in: f_n_bins 0D_NOT_integer
  if (c_associated(n_bins)) then
    call c_f_pointer(n_bins, f_n_bins_ptr)
  else
    f_n_bins_ptr => null()
  endif
  ! out: f_binned_data 0D_NOT_type
  if (.not. c_associated(binned_data)) return
  call c_f_pointer(binned_data, f_binned_data)
  f_binned_data = bin_data(f_data, f_weight, f_min_ptr, f_max_ptr, f_n_bins_ptr)

  ! out: f_binned_data 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_bin_data_density (bin_data, x, order, r) bind(c)

  use array_desc_mod
  use bin_mod, only: bin_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(bin_struct), pointer :: f_bin_data
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: order  ! 0D_NOT_integer
  integer(c_int) :: f_order
  integer(c_int), pointer :: f_order_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: r  ! 0D_NOT_real
  real(rp) :: f_r
  real(c_double), pointer :: f_r_ptr
  ! ** End of parameters **
  ! in: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_order 0D_NOT_integer
  if (c_associated(order)) then
    call c_f_pointer(order, f_order_ptr)
  else
    f_order_ptr => null()
  endif
  f_r = bin_data_density(f_bin_data, f_x, f_order_ptr)

  ! out: f_r 0D_NOT_real
  call c_f_pointer(r, f_r_ptr)
  f_r_ptr = f_r
end subroutine
subroutine fortran_bin_data_density_2d (bin_data, x, y, order, r0) bind(c)

  use array_desc_mod
  use bin_mod, only: general_bin_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(general_bin_struct), pointer :: f_bin_data
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  real(c_double) :: y  ! 0D_NOT_real
  real(rp) :: f_y
  type(c_ptr), intent(in), value :: order  ! 0D_NOT_integer
  integer(c_int) :: f_order
  integer(c_int), pointer :: f_order_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: r0  ! 0D_NOT_real
  real(rp) :: f_r0
  real(c_double), pointer :: f_r0_ptr
  ! ** End of parameters **
  ! in: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_y 0D_NOT_real
  f_y = y
  ! in: f_order 0D_NOT_integer
  if (c_associated(order)) then
    call c_f_pointer(order, f_order_ptr)
  else
    f_order_ptr => null()
  endif
  f_r0 = bin_data_density_2d(f_bin_data, f_x, f_y, f_order_ptr)

  ! out: f_r0 0D_NOT_real
  call c_f_pointer(r0, f_r0_ptr)
  f_r0_ptr = f_r0
end subroutine
subroutine fortran_bin_index (x, bin1_x_min, bin_delta, ix_bin) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  real(c_double) :: bin1_x_min  ! 0D_NOT_real
  real(rp) :: f_bin1_x_min
  real(c_double) :: bin_delta  ! 0D_NOT_real
  real(rp) :: f_bin_delta
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_bin  ! 0D_NOT_integer
  integer :: f_ix_bin
  integer(c_int), pointer :: f_ix_bin_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_bin1_x_min 0D_NOT_real
  f_bin1_x_min = bin1_x_min
  ! in: f_bin_delta 0D_NOT_real
  f_bin_delta = bin_delta
  f_ix_bin = bin_index(f_x, f_bin1_x_min, f_bin_delta)

  ! out: f_ix_bin 0D_NOT_integer
  call c_f_pointer(ix_bin, f_ix_bin_ptr)
  f_ix_bin_ptr = f_ix_bin
end subroutine
subroutine fortran_bin_x_center (ix_bin, bin1_x_min, bin_delta, x_center) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: bin1_x_min  ! 0D_NOT_real
  real(rp) :: f_bin1_x_min
  real(c_double) :: bin_delta  ! 0D_NOT_real
  real(rp) :: f_bin_delta
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: x_center  ! 0D_NOT_real
  real(rp) :: f_x_center
  real(c_double), pointer :: f_x_center_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_bin  ! 0D_NOT_integer
  integer(c_int) :: f_ix_bin
  integer(c_int), pointer :: f_ix_bin_ptr
  ! ** End of parameters **
  ! inout: f_ix_bin 0D_NOT_integer
  if (c_associated(ix_bin)) then
    call c_f_pointer(ix_bin, f_ix_bin_ptr)
  else
    f_ix_bin_ptr => null()
  endif
  ! in: f_bin1_x_min 0D_NOT_real
  f_bin1_x_min = bin1_x_min
  ! in: f_bin_delta 0D_NOT_real
  f_bin_delta = bin_delta
  f_x_center = bin_x_center(f_ix_bin_ptr, f_bin1_x_min, f_bin_delta)

  ! inout: f_ix_bin 0D_NOT_integer
  ! no output conversion for f_ix_bin
  ! out: f_x_center 0D_NOT_real
  call c_f_pointer(x_center, f_x_center_ptr)
  f_x_center_ptr = f_x_center
end subroutine
subroutine fortran_bit_set (word, pos, set_to_1) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: pos  ! 0D_NOT_integer
  integer :: f_pos
  logical(c_bool) :: set_to_1  ! 0D_NOT_logical
  logical :: f_set_to_1
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: word  ! 0D_NOT_integer
  integer(c_int) :: f_word
  integer(c_int), pointer :: f_word_ptr
  ! ** End of parameters **
  ! inout: f_word 0D_NOT_integer
  if (c_associated(word)) then
    call c_f_pointer(word, f_word_ptr)
  else
    f_word_ptr => null()
  endif
  ! in: f_pos 0D_NOT_integer
  f_pos = pos
  ! in: f_set_to_1 0D_NOT_logical
  f_set_to_1 = set_to_1
  call bit_set(f_word_ptr, f_pos, f_set_to_1)

  ! inout: f_word 0D_NOT_integer
  ! no output conversion for f_word
end subroutine
subroutine fortran_bracket_index_for_spline (x_knot, x, ix0, strict, print_err, ok) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: x_knot
  real(rp), pointer :: f_x_knot(:)
  real(c_double), pointer :: f_x_knot_ptr(:)
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: strict  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_strict
  logical, target :: f_strict_native
  logical, pointer :: f_strict_native_ptr
  logical(c_bool), pointer :: f_strict_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix0  ! 0D_NOT_integer
  integer :: f_ix0
  integer(c_int), pointer :: f_ix0_ptr
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(x_knot%data_ptr)) then
    call c_f_pointer(x_knot%data_ptr, f_x_knot_ptr, [x_knot%dims(1)])
    f_x_knot => f_x_knot_ptr
  else
    f_x_knot => null()
  endif
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_strict 0D_NOT_logical
  if (c_associated(strict)) then
    call c_f_pointer(strict, f_strict_ptr)
    f_strict_native = f_strict_ptr
    f_strict_native_ptr => f_strict_native
  else
    f_strict_native_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  f_ok = bracket_index_for_spline(f_x_knot, f_x, f_ix0, f_strict_native_ptr, &
      f_print_err_native_ptr)

  ! out: f_ix0 0D_NOT_integer
  call c_f_pointer(ix0, f_ix0_ptr)
  f_ix0_ptr = f_ix0
  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
end subroutine
subroutine fortran_calc_file_number (file_name, num_in, num_out, err_flag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  integer(c_int) :: num_in  ! 0D_NOT_integer
  integer :: f_num_in
  integer(c_int) :: num_out  ! 0D_NOT_integer
  integer :: f_num_out
  logical(c_bool) :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  ! ** End of parameters **
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! in: f_num_in 0D_NOT_integer
  f_num_in = num_in
  ! in: f_num_out 0D_NOT_integer
  f_num_out = num_out
  ! in: f_err_flag 0D_NOT_logical
  f_err_flag = err_flag
  call calc_file_number(f_file_name, f_num_in, f_num_out, f_err_flag)

end subroutine
subroutine fortran_celbd (mc, elb, eld) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: mc  ! 0D_NOT_real
  real(dp) :: f_mc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: elb  ! 0D_NOT_real
  real(dp) :: f_elb
  real(c_double), pointer :: f_elb_ptr
  type(c_ptr), intent(in), value :: eld  ! 0D_NOT_real
  real(dp) :: f_eld
  real(c_double), pointer :: f_eld_ptr
  ! ** End of parameters **
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  call celbd(f_mc, f_elb, f_eld)

  ! out: f_elb 0D_NOT_real
  call c_f_pointer(elb, f_elb_ptr)
  f_elb_ptr = f_elb
  ! out: f_eld 0D_NOT_real
  call c_f_pointer(eld, f_eld_ptr)
  f_eld_ptr = f_eld
end subroutine
subroutine fortran_cesr_getarg (i_arg, arg) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: i_arg  ! 0D_NOT_integer
  integer :: f_i_arg
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arg
  character(len=4096), target :: f_arg
  character(kind=c_char), pointer :: f_arg_ptr(:)
  ! ** End of parameters **
  ! in: f_i_arg 0D_NOT_integer
  f_i_arg = i_arg
  call cesr_getarg(f_i_arg, f_arg)

  ! out: f_arg 0D_NOT_character
  call c_f_pointer(arg, f_arg_ptr, [len_trim(f_arg) + 1])
  call to_c_str(f_arg, f_arg_ptr)
end subroutine
subroutine fortran_cesr_iargc (func_retval__) bind(c)

  use array_desc_mod
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: func_retval__  ! 0D_NOT_integer
  integer :: f_func_retval__
  integer(c_int), pointer :: f_func_retval___ptr
  ! ** End of parameters **
  f_func_retval__ = cesr_iargc()

  ! out: f_func_retval__ 0D_NOT_integer
  call c_f_pointer(func_retval__, f_func_retval___ptr)
  f_func_retval___ptr = f_func_retval__
end subroutine
subroutine fortran_change_file_number (file_name, change) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  integer(c_int) :: change  ! 0D_NOT_integer
  integer :: f_change
  ! ** End of parameters **
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! in: f_change 0D_NOT_integer
  f_change = change
  call change_file_number(f_file_name, f_change)

end subroutine
subroutine fortran_charge_of (species, default_, charge) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  type(c_ptr), intent(in), value :: default_  ! 0D_NOT_integer
  integer(c_int) :: f_default
  integer(c_int), pointer :: f_default_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: charge  ! 0D_NOT_integer
  integer :: f_charge
  integer(c_int), pointer :: f_charge_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  ! in: f_default 0D_NOT_integer
  if (c_associated(default_)) then
    call c_f_pointer(default_, f_default_ptr)
  else
    f_default_ptr => null()
  endif
  f_charge = charge_of(f_species, f_default_ptr)

  ! out: f_charge 0D_NOT_integer
  call c_f_pointer(charge, f_charge_ptr)
  f_charge_ptr = f_charge
end subroutine
subroutine fortran_charge_to_mass_of (species, charge_mass_ratio) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: charge_mass_ratio  ! 0D_NOT_real
  real(rp) :: f_charge_mass_ratio
  real(c_double), pointer :: f_charge_mass_ratio_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_charge_mass_ratio = charge_to_mass_of(f_species)

  ! out: f_charge_mass_ratio 0D_NOT_real
  call c_f_pointer(charge_mass_ratio, f_charge_mass_ratio_ptr)
  f_charge_mass_ratio_ptr = f_charge_mass_ratio
end subroutine
subroutine fortran_coarse_frequency_estimate (data, error, frequency) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: data
  real(rp), pointer :: f_data(:)
  real(c_double), pointer :: f_data_ptr(:)
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_error
  logical, target :: f_error_native
  logical, pointer :: f_error_native_ptr
  logical(c_bool), pointer :: f_error_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: frequency  ! 0D_NOT_real
  real(rp) :: f_frequency
  real(c_double), pointer :: f_frequency_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(data%data_ptr)) then
    call c_f_pointer(data%data_ptr, f_data_ptr, [data%dims(1)])
    f_data => f_data_ptr
  else
    f_data => null()
  endif
  ! in: f_error 0D_NOT_logical
  if (c_associated(error)) then
    call c_f_pointer(error, f_error_ptr)
    f_error_native = f_error_ptr
    f_error_native_ptr => f_error_native
  else
    f_error_native_ptr => null()
  endif
  f_frequency = coarse_frequency_estimate(f_data, f_error_native_ptr)

  ! out: f_frequency 0D_NOT_real
  call c_f_pointer(frequency, f_frequency_ptr)
  f_frequency_ptr = f_frequency
end subroutine
subroutine fortran_complex_error_function (wr, wi, zr, zi) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: wr  ! 0D_NOT_real
  real(rp) :: f_wr
  real(c_double) :: wi  ! 0D_NOT_real
  real(rp) :: f_wi
  real(c_double) :: zr  ! 0D_NOT_real
  real(rp) :: f_zr
  real(c_double) :: zi  ! 0D_NOT_real
  real(rp) :: f_zi
  ! ** End of parameters **
  ! in: f_wr 0D_NOT_real
  f_wr = wr
  ! in: f_wi 0D_NOT_real
  f_wi = wi
  ! in: f_zr 0D_NOT_real
  f_zr = zr
  ! in: f_zi 0D_NOT_real
  f_zi = zi
  call complex_error_function(f_wr, f_wi, f_zr, f_zi)

end subroutine
subroutine fortran_cos_one (angle, cos1) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: cos1  ! 0D_NOT_real
  real(rp) :: f_cos1
  real(c_double), pointer :: f_cos1_ptr
  ! ** End of parameters **
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  f_cos1 = cos_one(f_angle)

  ! out: f_cos1 0D_NOT_real
  call c_f_pointer(cos1, f_cos1_ptr)
  f_cos1_ptr = f_cos1
end subroutine
subroutine fortran_cosc (x, nd, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_nd 0D_NOT_integer
  if (c_associated(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
  else
    f_nd_ptr => null()
  endif
  f_y = cosc(f_x, f_nd_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_count_at_index (bin_data, index, c) bind(c)

  use array_desc_mod
  use bin_mod, only: bin_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: index  ! 0D_NOT_integer
  integer :: f_index
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: c  ! 0D_NOT_real
  real(rp) :: f_c
  real(c_double), pointer :: f_c_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(bin_struct), pointer :: f_bin_data
  ! ** End of parameters **
  ! inout: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  ! in: f_index 0D_NOT_integer
  f_index = index
  f_c = count_at_index(f_bin_data, f_index)

  ! out: f_c 0D_NOT_real
  call c_f_pointer(c, f_c_ptr)
  f_c_ptr = f_c
end subroutine
subroutine fortran_create_a_spline (r0, r1, slope0, slope1, spline) bind(c)

  use array_desc_mod
  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: r0
  real(rp), pointer :: f_r0(:)
  real(c_double), pointer :: f_r0_ptr(:)
  type(array_descriptor_t), intent(in) :: r1
  real(rp), pointer :: f_r1(:)
  real(c_double), pointer :: f_r1_ptr(:)
  real(c_double) :: slope0  ! 0D_NOT_real
  real(rp) :: f_slope0
  real(c_double) :: slope1  ! 0D_NOT_real
  real(rp) :: f_slope1
  ! ** Out parameters **
  type(c_ptr), value :: spline  ! 0D_NOT_type
  type(spline_struct), pointer :: f_spline
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(r0%data_ptr)) then
    call c_f_pointer(r0%data_ptr, f_r0_ptr, [r0%dims(1)])
    f_r0 => f_r0_ptr
  else
    f_r0 => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(r1%data_ptr)) then
    call c_f_pointer(r1%data_ptr, f_r1_ptr, [r1%dims(1)])
    f_r1 => f_r1_ptr
  else
    f_r1 => null()
  endif
  ! in: f_slope0 0D_NOT_real
  f_slope0 = slope0
  ! in: f_slope1 0D_NOT_real
  f_slope1 = slope1
  ! out: f_spline 0D_NOT_type
  if (.not. c_associated(spline)) return
  call c_f_pointer(spline, f_spline)
  f_spline = create_a_spline(f_r0, f_r1, f_slope0, f_slope1)

  ! out: f_spline 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_cross_product (a, b, c) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: a
  real(rp), pointer :: f_a(:)
  real(c_double), pointer :: f_a_ptr(:)
  type(array_descriptor_t), intent(in) :: b
  real(rp), pointer :: f_b(:)
  real(c_double), pointer :: f_b_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: c
  real(rp) :: f_c(3)
  real(c_double), pointer :: f_c_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(a%data_ptr)) then
    call c_f_pointer(a%data_ptr, f_a_ptr, [a%dims(1)])
    f_a => f_a_ptr
  else
    f_a => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(b%data_ptr)) then
    call c_f_pointer(b%data_ptr, f_b_ptr, [b%dims(1)])
    f_b => f_b_ptr
  else
    f_b => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(c%data_ptr)) then
    call c_f_pointer(c%data_ptr, f_c_ptr, [c%dims(1)])
    ! output-only
  else
    f_c_ptr => null()
  endif
  f_c = cross_product(f_a, f_b)

  ! out: f_c 1D_NOT_real
  if (c_associated(c%data_ptr)) then
    call c_f_pointer(c%data_ptr, f_c_ptr, [c%dims(1)])
    f_c_ptr = f_c(:)
  endif
end subroutine
subroutine fortran_date_and_time_stamp (string, numeric_month, include_zone) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: numeric_month  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_numeric_month
  logical, target :: f_numeric_month_native
  logical, pointer :: f_numeric_month_native_ptr
  logical(c_bool), pointer :: f_numeric_month_ptr
  type(c_ptr), intent(in), value :: include_zone  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_zone
  logical, target :: f_include_zone_native
  logical, pointer :: f_include_zone_native_ptr
  logical(c_bool), pointer :: f_include_zone_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_numeric_month 0D_NOT_logical
  if (c_associated(numeric_month)) then
    call c_f_pointer(numeric_month, f_numeric_month_ptr)
    f_numeric_month_native = f_numeric_month_ptr
    f_numeric_month_native_ptr => f_numeric_month_native
  else
    f_numeric_month_native_ptr => null()
  endif
  ! in: f_include_zone 0D_NOT_logical
  if (c_associated(include_zone)) then
    call c_f_pointer(include_zone, f_include_zone_ptr)
    f_include_zone_native = f_include_zone_ptr
    f_include_zone_native_ptr => f_include_zone_native
  else
    f_include_zone_native_ptr => null()
  endif
  call date_and_time_stamp(f_string, f_numeric_month_native_ptr, f_include_zone_native_ptr)

end subroutine
subroutine fortran_destfixedwindowls (id) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: id  ! 0D_NOT_integer
  INTEGER :: f_id
  ! ** End of parameters **
  ! in: f_id 0D_NOT_integer
  f_id = id
  call destfixedwindowls(f_id)

end subroutine
subroutine fortran_detab (str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! in: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  call detab(f_str)

end subroutine
subroutine fortran_display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_screen  ! 0D_NOT_integer
  integer :: f_ix_screen
  real(c_double) :: x_size  ! 0D_NOT_real
  real(rp) :: f_x_size
  real(c_double) :: y_size  ! 0D_NOT_real
  real(rp) :: f_y_size
  real(c_double) :: x_res  ! 0D_NOT_real
  real(rp) :: f_x_res
  real(c_double) :: y_res  ! 0D_NOT_real
  real(rp) :: f_y_res
  ! ** End of parameters **
  ! in: f_ix_screen 0D_NOT_integer
  f_ix_screen = ix_screen
  ! in: f_x_size 0D_NOT_real
  f_x_size = x_size
  ! in: f_y_size 0D_NOT_real
  f_y_size = y_size
  ! in: f_x_res 0D_NOT_real
  f_x_res = x_res
  ! in: f_y_res 0D_NOT_real
  f_y_res = y_res
  call display_size_and_resolution(f_ix_screen, f_x_size, f_y_size, f_x_res, f_y_res)

end subroutine
subroutine fortran_dj_bessel (m, arg, dj_bes) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: m  ! 0D_NOT_integer
  integer :: f_m
  real(c_double) :: arg  ! 0D_NOT_real
  real(rp) :: f_arg
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: dj_bes  ! 0D_NOT_real
  real(rp) :: f_dj_bes
  real(c_double), pointer :: f_dj_bes_ptr
  ! ** End of parameters **
  ! in: f_m 0D_NOT_integer
  f_m = m
  ! in: f_arg 0D_NOT_real
  f_arg = arg
  f_dj_bes = dj_bessel(f_m, f_arg)

  ! out: f_dj_bes 0D_NOT_real
  call c_f_pointer(dj_bes, f_dj_bes_ptr)
  f_dj_bes_ptr = f_dj_bes
end subroutine
subroutine fortran_djb_hash (str, old_hash, hash) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: old_hash  ! 0D_NOT_integer
  integer(c_int) :: f_old_hash
  integer(c_int), pointer :: f_old_hash_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: hash  ! 0D_NOT_integer
  integer :: f_hash
  integer(c_int), pointer :: f_hash_ptr
  ! ** End of parameters **
  ! in: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! in: f_old_hash 0D_NOT_integer
  if (c_associated(old_hash)) then
    call c_f_pointer(old_hash, f_old_hash_ptr)
  else
    f_old_hash_ptr => null()
  endif
  f_hash = djb_hash(f_str, f_old_hash_ptr)

  ! out: f_hash 0D_NOT_integer
  call c_f_pointer(hash, f_hash_ptr)
  f_hash_ptr = f_hash
end subroutine
subroutine fortran_djb_str_hash (in_str, hash_str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096), target :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: hash_str
  character(len=4096), target :: f_hash_str
  character(kind=c_char), pointer :: f_hash_str_ptr(:)
  ! ** End of parameters **
  ! in: f_in_str 0D_NOT_character
  if (.not. c_associated(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  f_hash_str = djb_str_hash(f_in_str)

  ! out: f_hash_str 0D_NOT_character
  call c_f_pointer(hash_str, f_hash_str_ptr, [len_trim(f_hash_str) + 1])
  call to_c_str(f_hash_str, f_hash_str_ptr)
end subroutine
subroutine fortran_downcase_string (string) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  call downcase_string(f_string)

end subroutine
subroutine fortran_elbd (phi, phic, mc, b, d) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: phi  ! 0D_NOT_real
  real(dp) :: f_phi
  real(c_double) :: phic  ! 0D_NOT_real
  real(dp) :: f_phic
  real(c_double) :: mc  ! 0D_NOT_real
  real(dp) :: f_mc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: b  ! 0D_NOT_real
  real(dp) :: f_b
  real(c_double), pointer :: f_b_ptr
  type(c_ptr), intent(in), value :: d  ! 0D_NOT_real
  real(dp) :: f_d
  real(c_double), pointer :: f_d_ptr
  ! ** End of parameters **
  ! in: f_phi 0D_NOT_real
  f_phi = phi
  ! in: f_phic 0D_NOT_real
  f_phic = phic
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  call elbd(f_phi, f_phic, f_mc, f_b, f_d)

  ! out: f_b 0D_NOT_real
  call c_f_pointer(b, f_b_ptr)
  f_b_ptr = f_b
  ! out: f_d 0D_NOT_real
  call c_f_pointer(d, f_d_ptr)
  f_d_ptr = f_d
end subroutine
subroutine fortran_elcbd (c0, mc, b, dx) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: c0  ! 0D_NOT_real
  real(dp) :: f_c0
  real(c_double) :: mc  ! 0D_NOT_real
  real(dp) :: f_mc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: b  ! 0D_NOT_real
  real(dp) :: f_b
  real(c_double), pointer :: f_b_ptr
  type(c_ptr), intent(in), value :: dx  ! 0D_NOT_real
  real(dp) :: f_dx
  real(c_double), pointer :: f_dx_ptr
  ! ** End of parameters **
  ! in: f_c0 0D_NOT_real
  f_c0 = c0
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  call elcbd(f_c0, f_mc, f_b, f_dx)

  ! out: f_b 0D_NOT_real
  call c_f_pointer(b, f_b_ptr)
  f_b_ptr = f_b
  ! out: f_dx 0D_NOT_real
  call c_f_pointer(dx, f_dx_ptr)
  f_dx_ptr = f_dx
end subroutine
subroutine fortran_ellipinc (phi, m, ellipkinc, ellipeinc) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: phi  ! 0D_NOT_real
  real(dp) :: f_phi
  real(c_double) :: m  ! 0D_NOT_real
  real(dp) :: f_m
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ellipkinc  ! 0D_NOT_real
  real(dp) :: f_ellipkinc
  real(c_double), pointer :: f_ellipkinc_ptr
  type(c_ptr), intent(in), value :: ellipeinc  ! 0D_NOT_real
  real(dp) :: f_ellipeinc
  real(c_double), pointer :: f_ellipeinc_ptr
  ! ** End of parameters **
  ! in: f_phi 0D_NOT_real
  f_phi = phi
  ! in: f_m 0D_NOT_real
  f_m = m
  call ellipinc(f_phi, f_m, f_ellipkinc, f_ellipeinc)

  ! out: f_ellipkinc 0D_NOT_real
  call c_f_pointer(ellipkinc, f_ellipkinc_ptr)
  f_ellipkinc_ptr = f_ellipkinc
  ! out: f_ellipeinc 0D_NOT_real
  call c_f_pointer(ellipeinc, f_ellipeinc_ptr)
  f_ellipeinc_ptr = f_ellipeinc
end subroutine
subroutine fortran_elsbd (s0, mc, b, d) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: s0  ! 0D_NOT_real
  real(dp) :: f_s0
  real(c_double) :: mc  ! 0D_NOT_real
  real(dp) :: f_mc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: b  ! 0D_NOT_real
  real(dp) :: f_b
  real(c_double), pointer :: f_b_ptr
  type(c_ptr), intent(in), value :: d  ! 0D_NOT_real
  real(dp) :: f_d
  real(c_double), pointer :: f_d_ptr
  ! ** End of parameters **
  ! in: f_s0 0D_NOT_real
  f_s0 = s0
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  call elsbd(f_s0, f_mc, f_b, f_d)

  ! out: f_b 0D_NOT_real
  call c_f_pointer(b, f_b_ptr)
  f_b_ptr = f_b
  ! out: f_d 0D_NOT_real
  call c_f_pointer(d, f_d_ptr)
  f_d_ptr = f_d
end subroutine
subroutine fortran_end_akima_spline_calc (spline, which_end) bind(c)

  use array_desc_mod
  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: which_end  ! 0D_NOT_integer
  integer :: f_which_end
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: spline
  type(spline_struct), pointer :: f_spline(:)
  type(spline_struct), pointer :: f_spline_ptr(:)
  ! ** End of parameters **
  !! type array (1D_NOT_type)
  if (c_associated(spline%data_ptr)) then
    call c_f_pointer(spline%data_ptr, f_spline_ptr, [spline%dims(1)])
    f_spline => f_spline_ptr
  else
    f_spline => null()
  endif
  ! in: f_which_end 0D_NOT_integer
  f_which_end = which_end
  call end_akima_spline_calc(f_spline, f_which_end)

end subroutine
subroutine fortran_err_exit (err_str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: err_str
  character(len=4096), target :: f_err_str
  character(kind=c_char), pointer :: f_err_str_ptr(:)
  character(len=4096), pointer :: f_err_str_call_ptr
  ! ** End of parameters **
  ! in: f_err_str 0D_NOT_character
  if (c_associated(err_str)) then
    call c_f_pointer(err_str, f_err_str_ptr, [huge(0)])
    call to_f_str(f_err_str_ptr, f_err_str)
    f_err_str_call_ptr => f_err_str
  else
    f_err_str_call_ptr => null()
  endif
  call err_exit(f_err_str_call_ptr)

end subroutine
subroutine fortran_factorial (n, fact) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: n  ! 0D_NOT_integer
  integer :: f_n
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: fact  ! 0D_NOT_real
  real(rp) :: f_fact
  real(c_double), pointer :: f_fact_ptr
  ! ** End of parameters **
  ! in: f_n 0D_NOT_integer
  f_n = n
  f_fact = factorial(f_n)

  ! out: f_fact 0D_NOT_real
  call c_f_pointer(fact, f_fact_ptr)
  f_fact_ptr = f_fact
end subroutine
subroutine fortran_faddeeva_function (z, w, dw) bind(c)

  use array_desc_mod
  implicit none
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: z
  real(rp) :: f_z(2)
  real(c_double), pointer :: f_z_ptr(:)
  type(array_descriptor_t), intent(in) :: w
  real(rp) :: f_w(2)
  real(c_double), pointer :: f_w_ptr(:)
  type(array_descriptor_t), intent(in) :: dw
  real(rp) :: f_dw(2,2)
  real(c_double), pointer :: f_dw_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(z%data_ptr)) then
    call c_f_pointer(z%data_ptr, f_z_ptr, [z%dims(1)])
    f_z = f_z_ptr(:)
  else
    f_z_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(w%data_ptr)) then
    call c_f_pointer(w%data_ptr, f_w_ptr, [w%dims(1)])
    f_w = f_w_ptr(:)
  else
    f_w_ptr => null()
  endif
  !! general array (2D_NOT_real) inout
  if (c_associated(dw%data_ptr)) then
    call c_f_pointer(dw%data_ptr, f_dw_ptr, [product(dw%dims(1:dw%rank))])
    call vec2mat(f_dw_ptr, f_dw)
  else
    f_dw_ptr => null()
  endif
  call faddeeva_function(f_z, f_w, f_dw)

end subroutine
subroutine fortran_fff_sub (line, error) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  logical(c_bool) :: error  ! 0D_NOT_logical
  logical :: f_error
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_error 0D_NOT_logical
  f_error = error
  call fff_sub(f_line, f_error)

end subroutine
subroutine fortran_fft_1d (arr, isign) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: isign  ! 0D_NOT_integer
  integer :: f_isign
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: arr
  complex(rp), pointer :: f_arr(:)
  complex(c_double_complex), pointer :: f_arr_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) inout
  if (c_associated(arr%data_ptr)) then
    call c_f_pointer(arr%data_ptr, f_arr_ptr, [arr%dims(1)])
    f_arr => f_arr_ptr
  else
    f_arr => null()
  endif
  ! in: f_isign 0D_NOT_integer
  f_isign = isign
  call fft_1d(f_arr, f_isign)

end subroutine
subroutine fortran_file_directorizer (in_file, out_file, directory, add_switch) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_file
  character(len=4096), target :: f_in_file
  character(kind=c_char), pointer :: f_in_file_ptr(:)
  type(c_ptr), intent(in), value :: out_file
  character(len=4096), target :: f_out_file
  character(kind=c_char), pointer :: f_out_file_ptr(:)
  type(c_ptr), intent(in), value :: directory
  character(len=4096), target :: f_directory
  character(kind=c_char), pointer :: f_directory_ptr(:)
  logical(c_bool) :: add_switch  ! 0D_NOT_logical
  logical :: f_add_switch
  ! ** End of parameters **
  ! in: f_in_file 0D_NOT_character
  if (.not. c_associated(in_file)) return
  call c_f_pointer(in_file, f_in_file_ptr, [huge(0)])
  call to_f_str(f_in_file_ptr, f_in_file)
  ! in: f_out_file 0D_NOT_character
  if (.not. c_associated(out_file)) return
  call c_f_pointer(out_file, f_out_file_ptr, [huge(0)])
  call to_f_str(f_out_file_ptr, f_out_file)
  ! in: f_directory 0D_NOT_character
  if (.not. c_associated(directory)) return
  call c_f_pointer(directory, f_directory_ptr, [huge(0)])
  call to_f_str(f_directory_ptr, f_directory)
  ! in: f_add_switch 0D_NOT_logical
  f_add_switch = add_switch
  call file_directorizer(f_in_file, f_out_file, f_directory, f_add_switch)

end subroutine
subroutine fortran_file_get (string, dflt_file_name, file_name) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: dflt_file_name
  character(len=4096), target :: f_dflt_file_name
  character(kind=c_char), pointer :: f_dflt_file_name_ptr(:)
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_dflt_file_name 0D_NOT_character
  if (.not. c_associated(dflt_file_name)) return
  call c_f_pointer(dflt_file_name, f_dflt_file_name_ptr, [huge(0)])
  call to_f_str(f_dflt_file_name_ptr, f_dflt_file_name)
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  call file_get(f_string, f_dflt_file_name, f_file_name)

end subroutine
subroutine fortran_file_get_open (string, dflt_file_name, file_name, file_unit, readonly) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: dflt_file_name
  character(len=4096), target :: f_dflt_file_name
  character(kind=c_char), pointer :: f_dflt_file_name_ptr(:)
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  integer(c_int) :: file_unit  ! 0D_NOT_integer
  integer :: f_file_unit
  logical(c_bool) :: readonly  ! 0D_NOT_logical
  logical :: f_readonly
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_dflt_file_name 0D_NOT_character
  if (.not. c_associated(dflt_file_name)) return
  call c_f_pointer(dflt_file_name, f_dflt_file_name_ptr, [huge(0)])
  call to_f_str(f_dflt_file_name_ptr, f_dflt_file_name)
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! in: f_file_unit 0D_NOT_integer
  f_file_unit = file_unit
  ! in: f_readonly 0D_NOT_logical
  f_readonly = readonly
  call file_get_open(f_string, f_dflt_file_name, f_file_name, f_file_unit, f_readonly)

end subroutine
subroutine fortran_file_suffixer (in_file_name, out_file_name, suffix, add_switch) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_file_name
  character(len=4096), target :: f_in_file_name
  character(kind=c_char), pointer :: f_in_file_name_ptr(:)
  type(c_ptr), intent(in), value :: out_file_name
  character(len=4096), target :: f_out_file_name
  character(kind=c_char), pointer :: f_out_file_name_ptr(:)
  type(c_ptr), intent(in), value :: suffix
  character(len=4096), target :: f_suffix
  character(kind=c_char), pointer :: f_suffix_ptr(:)
  logical(c_bool) :: add_switch  ! 0D_NOT_logical
  logical :: f_add_switch
  ! ** End of parameters **
  ! in: f_in_file_name 0D_NOT_character
  if (.not. c_associated(in_file_name)) return
  call c_f_pointer(in_file_name, f_in_file_name_ptr, [huge(0)])
  call to_f_str(f_in_file_name_ptr, f_in_file_name)
  ! in: f_out_file_name 0D_NOT_character
  if (.not. c_associated(out_file_name)) return
  call c_f_pointer(out_file_name, f_out_file_name_ptr, [huge(0)])
  call to_f_str(f_out_file_name_ptr, f_out_file_name)
  ! in: f_suffix 0D_NOT_character
  if (.not. c_associated(suffix)) return
  call c_f_pointer(suffix, f_suffix_ptr, [huge(0)])
  call to_f_str(f_suffix_ptr, f_suffix)
  ! in: f_add_switch 0D_NOT_logical
  f_add_switch = add_switch
  call file_suffixer(f_in_file_name, f_out_file_name, f_suffix, f_add_switch)

end subroutine
subroutine fortran_find_location_int (arr, value, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: value  ! 0D_NOT_integer
  integer :: f_value
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: arr
  integer, pointer :: f_arr(:)
  integer(c_int), pointer :: f_arr_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_integer) inout
  if (c_associated(arr%data_ptr)) then
    call c_f_pointer(arr%data_ptr, f_arr_ptr, [arr%dims(1)])
    f_arr => f_arr_ptr
  else
    f_arr => null()
  endif
  ! in: f_value 0D_NOT_integer
  f_value = value
  f_ix_match = find_location(f_arr, f_value)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_find_location_logic (arr, value, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: value  ! 0D_NOT_logical
  logical :: f_value
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr
  type(logical_container_alloc), pointer :: f_arr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_logical)
  if (c_associated(arr))   call c_f_pointer(arr, f_arr)
  ! in: f_value 0D_NOT_logical
  f_value = value
  f_ix_match = find_location(f_arr%data, f_value)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_find_location_real (arr, value, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: arr
  real(rp), pointer :: f_arr(:)
  real(c_double), pointer :: f_arr_ptr(:)
  real(c_double) :: value  ! 0D_NOT_real
  real(rp) :: f_value
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(arr%data_ptr)) then
    call c_f_pointer(arr%data_ptr, f_arr_ptr, [arr%dims(1)])
    f_arr => f_arr_ptr
  else
    f_arr => null()
  endif
  ! in: f_value 0D_NOT_real
  f_value = value
  f_ix_match = find_location(f_arr, f_value)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_find_location_str (arr, value, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: value
  character(len=4096), target :: f_value
  character(kind=c_char), pointer :: f_value_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr
  type(character_container_alloc), pointer :: f_arr
  character(200), allocatable :: f_arr_local(:)
  ! ** End of parameters **
  !! container character array (1D_NOT_character)
  if (c_associated(arr))   call c_f_pointer(arr, f_arr)
  if (c_associated(arr) .and. allocated(f_arr%data)) then
    allocate(f_arr_local, mold=f_arr%data)
    f_arr_local = ''
  endif
  ! in: f_value 0D_NOT_character
  if (.not. c_associated(value)) return
  call c_f_pointer(value, f_value_ptr, [huge(0)])
  call to_f_str(f_value_ptr, f_value)
  f_ix_match = find_location(f_arr_local, f_value)

  !! copy allocatable character result into container
  if (c_associated(arr) .and. allocated(f_arr_local)) then
    if (allocated(f_arr%data)) deallocate(f_arr%data)
    allocate(f_arr%data, source=f_arr_local)
  endif
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_fine_frequency_estimate (data, frequency) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: data
  real(rp), pointer :: f_data(:)
  real(c_double), pointer :: f_data_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: frequency  ! 0D_NOT_real
  real(rp) :: f_frequency
  real(c_double), pointer :: f_frequency_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(data%data_ptr)) then
    call c_f_pointer(data%data_ptr, f_data_ptr, [data%dims(1)])
    f_data => f_data_ptr
  else
    f_data => null()
  endif
  f_frequency = fine_frequency_estimate(f_data)

  ! out: f_frequency 0D_NOT_real
  call c_f_pointer(frequency, f_frequency_ptr)
  f_frequency_ptr = f_frequency
end subroutine
subroutine fortran_fixedwindowls (ynew, id, z) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: ynew  ! 0D_NOT_real
  REAL(rp) :: f_ynew
  integer(c_int) :: id  ! 0D_NOT_integer
  INTEGER :: f_id
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: z  ! 0D_NOT_real
  REAL(rp) :: f_z
  real(c_double), pointer :: f_z_ptr
  ! ** End of parameters **
  ! in: f_ynew 0D_NOT_real
  f_ynew = ynew
  ! in: f_id 0D_NOT_integer
  f_id = id
  f_z = fixedwindowls(f_ynew, f_id)

  ! out: f_z 0D_NOT_real
  call c_f_pointer(z, f_z_ptr)
  f_z_ptr = f_z
end subroutine
subroutine fortran_fourier_amplitude (data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: data
  real(rp), pointer :: f_data(:)
  real(c_double), pointer :: f_data_ptr(:)
  real(c_double) :: frequency  ! 0D_NOT_real
  real(rp) :: f_frequency
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: cos_amp  ! 0D_NOT_real
  real(rp) :: f_cos_amp
  real(c_double), pointer :: f_cos_amp_ptr
  type(c_ptr), intent(in), value :: sin_amp  ! 0D_NOT_real
  real(rp) :: f_sin_amp
  real(c_double), pointer :: f_sin_amp_ptr
  type(c_ptr), intent(in), value :: dcos_amp  ! 0D_NOT_real
  real(rp) :: f_dcos_amp
  real(c_double), pointer :: f_dcos_amp_ptr
  type(c_ptr), intent(in), value :: dsin_amp  ! 0D_NOT_real
  real(rp) :: f_dsin_amp
  real(c_double), pointer :: f_dsin_amp_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(data%data_ptr)) then
    call c_f_pointer(data%data_ptr, f_data_ptr, [data%dims(1)])
    f_data => f_data_ptr
  else
    f_data => null()
  endif
  ! in: f_frequency 0D_NOT_real
  f_frequency = frequency
  ! out: f_dcos_amp 0D_NOT_real
  if (c_associated(dcos_amp)) then
    call c_f_pointer(dcos_amp, f_dcos_amp_ptr)
  else
    f_dcos_amp_ptr => null()
  endif
  ! out: f_dsin_amp 0D_NOT_real
  if (c_associated(dsin_amp)) then
    call c_f_pointer(dsin_amp, f_dsin_amp_ptr)
  else
    f_dsin_amp_ptr => null()
  endif
  call fourier_amplitude(f_data, f_frequency, f_cos_amp, f_sin_amp, f_dcos_amp, f_dsin_amp)

  ! out: f_cos_amp 0D_NOT_real
  call c_f_pointer(cos_amp, f_cos_amp_ptr)
  f_cos_amp_ptr = f_cos_amp
  ! out: f_sin_amp 0D_NOT_real
  call c_f_pointer(sin_amp, f_sin_amp_ptr)
  f_sin_amp_ptr = f_sin_amp
  ! out: f_dcos_amp 0D_NOT_real
  call c_f_pointer(dcos_amp, f_dcos_amp_ptr)
  f_dcos_amp_ptr = f_dcos_amp
  ! out: f_dsin_amp 0D_NOT_real
  call c_f_pointer(dsin_amp, f_dsin_amp_ptr)
  f_dsin_amp_ptr = f_dsin_amp
end subroutine
subroutine fortran_gelbd (phi, mc, elb, eld) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: phi  ! 0D_NOT_real
  real(dp) :: f_phi
  real(c_double) :: mc  ! 0D_NOT_real
  real(dp) :: f_mc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: elb  ! 0D_NOT_real
  real(dp) :: f_elb
  real(c_double), pointer :: f_elb_ptr
  type(c_ptr), intent(in), value :: eld  ! 0D_NOT_real
  real(dp) :: f_eld
  real(c_double), pointer :: f_eld_ptr
  ! ** End of parameters **
  ! in: f_phi 0D_NOT_real
  f_phi = phi
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  call gelbd(f_phi, f_mc, f_elb, f_eld)

  ! out: f_elb 0D_NOT_real
  call c_f_pointer(elb, f_elb_ptr)
  f_elb_ptr = f_elb
  ! out: f_eld 0D_NOT_real
  call c_f_pointer(eld, f_eld_ptr)
  f_eld_ptr = f_eld
end subroutine
subroutine fortran_gen_complete_elliptic (kc, p, c, s, err_tol, value) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: kc  ! 0D_NOT_real
  real(rp) :: f_kc
  real(c_double) :: p  ! 0D_NOT_real
  real(rp) :: f_p
  real(c_double) :: c  ! 0D_NOT_real
  real(rp) :: f_c
  real(c_double) :: s  ! 0D_NOT_real
  real(rp) :: f_s
  type(c_ptr), intent(in), value :: err_tol  ! 0D_NOT_real
  real(c_double) :: f_err_tol
  real(c_double), pointer :: f_err_tol_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_kc 0D_NOT_real
  f_kc = kc
  ! in: f_p 0D_NOT_real
  f_p = p
  ! in: f_c 0D_NOT_real
  f_c = c
  ! in: f_s 0D_NOT_real
  f_s = s
  ! in: f_err_tol 0D_NOT_real
  if (c_associated(err_tol)) then
    call c_f_pointer(err_tol, f_err_tol_ptr)
  else
    f_err_tol_ptr => null()
  endif
  f_value = gen_complete_elliptic(f_kc, f_p, f_c, f_s, f_err_tol_ptr)

  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_general_bin_count (bin_data, ix1, ix2, ix3, count) bind(c)

  use array_desc_mod
  use bin_mod, only: general_bin_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix1  ! 0D_NOT_integer
  integer :: f_ix1
  type(c_ptr), intent(in), value :: ix2  ! 0D_NOT_integer
  integer(c_int) :: f_ix2
  integer(c_int), pointer :: f_ix2_ptr
  type(c_ptr), intent(in), value :: ix3  ! 0D_NOT_integer
  integer(c_int) :: f_ix3
  integer(c_int), pointer :: f_ix3_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: count  ! 0D_NOT_real
  real(rp) :: f_count
  real(c_double), pointer :: f_count_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(general_bin_struct), pointer :: f_bin_data
  ! ** End of parameters **
  ! inout: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  ! in: f_ix1 0D_NOT_integer
  f_ix1 = ix1
  ! in: f_ix2 0D_NOT_integer
  if (c_associated(ix2)) then
    call c_f_pointer(ix2, f_ix2_ptr)
  else
    f_ix2_ptr => null()
  endif
  ! in: f_ix3 0D_NOT_integer
  if (c_associated(ix3)) then
    call c_f_pointer(ix3, f_ix3_ptr)
  else
    f_ix3_ptr => null()
  endif
  f_count = general_bin_count(f_bin_data, f_ix1, f_ix2_ptr, f_ix3_ptr)

  ! out: f_count 0D_NOT_real
  call c_f_pointer(count, f_count_ptr)
  f_count_ptr = f_count
end subroutine
subroutine fortran_general_bin_index (bin_data, ix1, ix2, ix3, index) bind(c)

  use array_desc_mod
  use bin_mod, only: general_bin_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix1  ! 0D_NOT_integer
  integer :: f_ix1
  type(c_ptr), intent(in), value :: ix2  ! 0D_NOT_integer
  integer(c_int) :: f_ix2
  integer(c_int), pointer :: f_ix2_ptr
  type(c_ptr), intent(in), value :: ix3  ! 0D_NOT_integer
  integer(c_int) :: f_ix3
  integer(c_int), pointer :: f_ix3_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: index  ! 0D_NOT_integer
  integer :: f_index
  integer(c_int), pointer :: f_index_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(general_bin_struct), pointer :: f_bin_data
  ! ** End of parameters **
  ! inout: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  ! in: f_ix1 0D_NOT_integer
  f_ix1 = ix1
  ! in: f_ix2 0D_NOT_integer
  if (c_associated(ix2)) then
    call c_f_pointer(ix2, f_ix2_ptr)
  else
    f_ix2_ptr => null()
  endif
  ! in: f_ix3 0D_NOT_integer
  if (c_associated(ix3)) then
    call c_f_pointer(ix3, f_ix3_ptr)
  else
    f_ix3_ptr => null()
  endif
  f_index = general_bin_index(f_bin_data, f_ix1, f_ix2_ptr, f_ix3_ptr)

  ! out: f_index 0D_NOT_integer
  call c_f_pointer(index, f_index_ptr)
  f_index_ptr = f_index
end subroutine
subroutine fortran_general_bin_index_in_bounds (bin_data, ix1, ix2, ix3, in_bounds) bind(c)

  use array_desc_mod
  use bin_mod, only: general_bin_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix1  ! 0D_NOT_integer
  integer :: f_ix1
  type(c_ptr), intent(in), value :: ix2  ! 0D_NOT_integer
  integer(c_int) :: f_ix2
  integer(c_int), pointer :: f_ix2_ptr
  type(c_ptr), intent(in), value :: ix3  ! 0D_NOT_integer
  integer(c_int) :: f_ix3
  integer(c_int), pointer :: f_ix3_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: in_bounds  ! 0D_NOT_logical
  logical :: f_in_bounds
  logical(c_bool), pointer :: f_in_bounds_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: bin_data  ! 0D_NOT_type
  type(general_bin_struct), pointer :: f_bin_data
  ! ** End of parameters **
  ! inout: f_bin_data 0D_NOT_type
  if (.not. c_associated(bin_data)) return
  call c_f_pointer(bin_data, f_bin_data)
  ! in: f_ix1 0D_NOT_integer
  f_ix1 = ix1
  ! in: f_ix2 0D_NOT_integer
  if (c_associated(ix2)) then
    call c_f_pointer(ix2, f_ix2_ptr)
  else
    f_ix2_ptr => null()
  endif
  ! in: f_ix3 0D_NOT_integer
  if (c_associated(ix3)) then
    call c_f_pointer(ix3, f_ix3_ptr)
  else
    f_ix3_ptr => null()
  endif
  f_in_bounds = general_bin_index_in_bounds(f_bin_data, f_ix1, f_ix2_ptr, f_ix3_ptr)

  ! out: f_in_bounds 0D_NOT_logical
  call c_f_pointer(in_bounds, f_in_bounds_ptr)
  f_in_bounds_ptr = f_in_bounds
end subroutine
subroutine fortran_get_a_char (this_char, wait, ignore_this) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: wait  ! 0D_NOT_logical
  logical :: f_wait
  type(c_ptr), intent(in), value :: ignore_this
  type(character_container_alloc), pointer :: f_ignore_this
  character(200), allocatable :: f_ignore_this_local(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_char
  character(len=4096), target :: f_this_char
  character(kind=c_char), pointer :: f_this_char_ptr(:)
  ! ** End of parameters **
  ! in: f_wait 0D_NOT_logical
  f_wait = wait
  !! container character array (1D_NOT_character)
  if (c_associated(ignore_this))   call c_f_pointer(ignore_this, f_ignore_this)
  if (c_associated(ignore_this) .and. allocated(f_ignore_this%data)) then
    allocate(f_ignore_this_local, mold=f_ignore_this%data)
    f_ignore_this_local = ''
  endif
  call get_a_char(f_this_char, f_wait, f_ignore_this_local)

  ! out: f_this_char 0D_NOT_character
  call c_f_pointer(this_char, f_this_char_ptr, [len_trim(f_this_char) + 1])
  call to_c_str(f_this_char, f_this_char_ptr)
  !! copy allocatable character result into container
  if (c_associated(ignore_this) .and. allocated(f_ignore_this_local)) then
    if (allocated(f_ignore_this%data)) deallocate(f_ignore_this%data)
    allocate(f_ignore_this%data, source=f_ignore_this_local)
  endif
end subroutine
subroutine fortran_get_file_number (file_name, cnum_in, num_out, err_flag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  type(c_ptr), intent(in), value :: cnum_in
  character(len=4096), target :: f_cnum_in
  character(kind=c_char), pointer :: f_cnum_in_ptr(:)
  integer(c_int) :: num_out  ! 0D_NOT_integer
  integer :: f_num_out
  logical(c_bool) :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  ! ** End of parameters **
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! in: f_cnum_in 0D_NOT_character
  if (.not. c_associated(cnum_in)) return
  call c_f_pointer(cnum_in, f_cnum_in_ptr, [huge(0)])
  call to_f_str(f_cnum_in_ptr, f_cnum_in)
  ! in: f_num_out 0D_NOT_integer
  f_num_out = num_out
  ! in: f_err_flag 0D_NOT_logical
  f_err_flag = err_flag
  call get_file_number(f_file_name, f_cnum_in, f_num_out, f_err_flag)

end subroutine
subroutine fortran_get_file_time_stamp (file, time_stamp) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file
  character(len=4096), target :: f_file
  character(kind=c_char), pointer :: f_file_ptr(:)
  type(c_ptr), intent(in), value :: time_stamp
  character(len=4096), target :: f_time_stamp
  character(kind=c_char), pointer :: f_time_stamp_ptr(:)
  ! ** End of parameters **
  ! in: f_file 0D_NOT_character
  if (.not. c_associated(file)) return
  call c_f_pointer(file, f_file_ptr, [huge(0)])
  call to_f_str(f_file_ptr, f_file)
  ! in: f_time_stamp 0D_NOT_character
  if (.not. c_associated(time_stamp)) return
  call c_f_pointer(time_stamp, f_time_stamp_ptr, [huge(0)])
  call to_f_str(f_time_stamp_ptr, f_time_stamp)
  call get_file_time_stamp(f_file, f_time_stamp)

end subroutine
subroutine fortran_get_next_number (filein, cnum, digits) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: filein
  character(len=4096), target :: f_filein
  character(kind=c_char), pointer :: f_filein_ptr(:)
  type(c_ptr), intent(in), value :: cnum
  character(len=4096), target :: f_cnum
  character(kind=c_char), pointer :: f_cnum_ptr(:)
  integer(c_int) :: digits  ! 0D_NOT_integer
  integer :: f_digits
  ! ** End of parameters **
  ! in: f_filein 0D_NOT_character
  if (.not. c_associated(filein)) return
  call c_f_pointer(filein, f_filein_ptr, [huge(0)])
  call to_f_str(f_filein_ptr, f_filein)
  ! in: f_cnum 0D_NOT_character
  if (.not. c_associated(cnum)) return
  call c_f_pointer(cnum, f_cnum_ptr, [huge(0)])
  call to_f_str(f_cnum_ptr, f_cnum)
  ! in: f_digits 0D_NOT_integer
  f_digits = digits
  call get_next_number(f_filein, f_cnum, f_digits)

end subroutine
subroutine fortran_get_tty_char (this_char, wait, flush) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: wait  ! 0D_NOT_logical
  logical :: f_wait
  logical(c_bool) :: flush  ! 0D_NOT_logical
  logical :: f_flush
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_char
  character(len=4096), target :: f_this_char
  character(kind=c_char), pointer :: f_this_char_ptr(:)
  ! ** End of parameters **
  ! in: f_wait 0D_NOT_logical
  f_wait = wait
  ! in: f_flush 0D_NOT_logical
  f_flush = flush
  call get_tty_char(f_this_char, f_wait, f_flush)

  ! out: f_this_char 0D_NOT_character
  call c_f_pointer(this_char, f_this_char_ptr, [len_trim(f_this_char) + 1])
  call to_c_str(f_this_char, f_this_char_ptr)
end subroutine
subroutine fortran_hanhan (N, hh) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: N  ! 0D_NOT_integer
  integer :: f_N
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: hh
  real(rp), pointer :: f_hh(:)
  real(c_double), pointer :: f_hh_ptr(:)
  ! ** End of parameters **
  ! in: f_N 0D_NOT_integer
  f_N = N
  !! general array (1D_NOT_real) inout
  if (c_associated(hh%data_ptr)) then
    call c_f_pointer(hh%data_ptr, f_hh_ptr, [hh%dims(1)])
    f_hh => f_hh_ptr
  else
    f_hh => null()
  endif
  call hanhan(f_N, f_hh)

end subroutine
subroutine fortran_i_bessel (m, arg, i_bes) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: m  ! 0D_NOT_integer
  integer :: f_m
  real(c_double) :: arg  ! 0D_NOT_real
  real(rp) :: f_arg
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: i_bes  ! 0D_NOT_real
  real(rp) :: f_i_bes
  real(c_double), pointer :: f_i_bes_ptr
  ! ** End of parameters **
  ! in: f_m 0D_NOT_integer
  f_m = m
  ! in: f_arg 0D_NOT_real
  f_arg = arg
  f_i_bes = i_bessel(f_m, f_arg)

  ! out: f_i_bes 0D_NOT_real
  call c_f_pointer(i_bes, f_i_bes_ptr)
  f_i_bes_ptr = f_i_bes
end subroutine
subroutine fortran_i_bessel_extended (m, arg, i_bes) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: m  ! 0D_NOT_integer
  integer :: f_m
  real(c_double) :: arg  ! 0D_NOT_real
  real(rp) :: f_arg
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: i_bes  ! 0D_NOT_complex
  complex(rp) :: f_i_bes
  complex(c_double_complex), pointer :: f_i_bes_ptr
  ! ** End of parameters **
  ! in: f_m 0D_NOT_integer
  f_m = m
  ! in: f_arg 0D_NOT_real
  f_arg = arg
  f_i_bes = i_bessel_extended(f_m, f_arg)

  ! out: f_i_bes 0D_NOT_complex
  call c_f_pointer(i_bes, f_i_bes_ptr)
  f_i_bes_ptr = f_i_bes
end subroutine
subroutine fortran_increment_file_number (file_name, digits, number, cnumber) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  integer(c_int) :: digits  ! 0D_NOT_integer
  integer :: f_digits
  integer(c_int) :: number  ! 0D_NOT_integer
  integer :: f_number
  type(c_ptr), intent(in), value :: cnumber
  character(len=4096), target :: f_cnumber
  character(kind=c_char), pointer :: f_cnumber_ptr(:)
  ! ** End of parameters **
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! in: f_digits 0D_NOT_integer
  f_digits = digits
  ! in: f_number 0D_NOT_integer
  f_number = number
  ! in: f_cnumber 0D_NOT_character
  if (.not. c_associated(cnumber)) return
  call c_f_pointer(cnumber, f_cnumber_ptr, [huge(0)])
  call to_f_str(f_cnumber_ptr, f_cnumber)
  call increment_file_number(f_file_name, f_digits, f_number, f_cnumber)

end subroutine
subroutine fortran_index_nocase (string1, string2, indx) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string1
  character(len=4096), target :: f_string1
  character(kind=c_char), pointer :: f_string1_ptr(:)
  type(c_ptr), intent(in), value :: string2
  character(len=4096), target :: f_string2
  character(kind=c_char), pointer :: f_string2_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: indx  ! 0D_NOT_integer
  integer :: f_indx
  integer(c_int), pointer :: f_indx_ptr
  ! ** End of parameters **
  ! in: f_string1 0D_NOT_character
  if (.not. c_associated(string1)) return
  call c_f_pointer(string1, f_string1_ptr, [huge(0)])
  call to_f_str(f_string1_ptr, f_string1)
  ! in: f_string2 0D_NOT_character
  if (.not. c_associated(string2)) return
  call c_f_pointer(string2, f_string2_ptr, [huge(0)])
  call to_f_str(f_string2_ptr, f_string2)
  f_indx = index_nocase(f_string1, f_string2)

  ! out: f_indx 0D_NOT_integer
  call c_f_pointer(indx, f_indx_ptr)
  f_indx_ptr = f_indx
end subroutine
subroutine fortran_initfixedwindowls (N, dt, order, der, id) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: N  ! 0D_NOT_integer
  INTEGER :: f_N
  real(c_double) :: dt  ! 0D_NOT_real
  REAL(rp) :: f_dt
  integer(c_int) :: order  ! 0D_NOT_integer
  INTEGER :: f_order
  integer(c_int) :: der  ! 0D_NOT_integer
  INTEGER :: f_der
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: id  ! 0D_NOT_integer
  INTEGER :: f_id
  integer(c_int), pointer :: f_id_ptr
  ! ** End of parameters **
  ! in: f_N 0D_NOT_integer
  f_N = N
  ! in: f_dt 0D_NOT_real
  f_dt = dt
  ! in: f_order 0D_NOT_integer
  f_order = order
  ! in: f_der 0D_NOT_integer
  f_der = der
  f_id = initfixedwindowls(f_N, f_dt, f_order, f_der)

  ! out: f_id 0D_NOT_integer
  call c_f_pointer(id, f_id_ptr)
  f_id_ptr = f_id
end subroutine
subroutine fortran_initial_lmdif () bind(c)

  use array_desc_mod
  implicit none
  ! ** End of parameters **
  call initial_lmdif()

end subroutine
subroutine fortran_int_str (int_, width, str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: int_  ! 0D_NOT_integer
  integer :: f_int
  type(c_ptr), intent(in), value :: width  ! 0D_NOT_integer
  integer(c_int) :: f_width
  integer(c_int), pointer :: f_width_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! in: f_int 0D_NOT_integer
  f_int = int_
  ! in: f_width 0D_NOT_integer
  if (c_associated(width)) then
    call c_f_pointer(width, f_width_ptr)
  else
    f_width_ptr => null()
  endif
  f_str = int_str(f_int, f_width_ptr)

  ! out: f_str 0D_ALLOC_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1])
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_interpolated_fft (cdata, calc_ok, opt_dump_spectrum, opt_dump_index, &
    this_fft) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  type(c_ptr), intent(in), value :: opt_dump_spectrum  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_spectrum
  integer(c_int), pointer :: f_opt_dump_spectrum_ptr
  type(c_ptr), intent(in), value :: opt_dump_index  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_index
  integer(c_int), pointer :: f_opt_dump_index_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_fft  ! 0D_NOT_real
  real(rp) :: f_this_fft
  real(c_double), pointer :: f_this_fft_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: cdata
  complex(rp), pointer :: f_cdata(:)
  complex(c_double_complex), pointer :: f_cdata_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) inout
  if (c_associated(cdata%data_ptr)) then
    call c_f_pointer(cdata%data_ptr, f_cdata_ptr, [cdata%dims(1)])
    f_cdata => f_cdata_ptr
  else
    f_cdata => null()
  endif
  ! in: f_calc_ok 0D_NOT_logical
  f_calc_ok = calc_ok
  ! in: f_opt_dump_spectrum 0D_NOT_integer
  if (c_associated(opt_dump_spectrum)) then
    call c_f_pointer(opt_dump_spectrum, f_opt_dump_spectrum_ptr)
  else
    f_opt_dump_spectrum_ptr => null()
  endif
  ! in: f_opt_dump_index 0D_NOT_integer
  if (c_associated(opt_dump_index)) then
    call c_f_pointer(opt_dump_index, f_opt_dump_index_ptr)
  else
    f_opt_dump_index_ptr => null()
  endif
  f_this_fft = interpolated_fft(f_cdata, f_calc_ok, f_opt_dump_spectrum_ptr, &
      f_opt_dump_index_ptr)

  ! out: f_this_fft 0D_NOT_real
  call c_f_pointer(this_fft, f_this_fft_ptr)
  f_this_fft_ptr = f_this_fft
end subroutine
subroutine fortran_interpolated_fft_gsl (cdata, calc_ok, opt_dump_spectrum, opt_dump_index, &
    this_fft) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  type(c_ptr), intent(in), value :: opt_dump_spectrum  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_spectrum
  integer(c_int), pointer :: f_opt_dump_spectrum_ptr
  type(c_ptr), intent(in), value :: opt_dump_index  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_index
  integer(c_int), pointer :: f_opt_dump_index_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_fft  ! 0D_NOT_real
  real(rp) :: f_this_fft
  real(c_double), pointer :: f_this_fft_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: cdata
  complex(rp), pointer :: f_cdata(:)
  complex(c_double_complex), pointer :: f_cdata_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) inout
  if (c_associated(cdata%data_ptr)) then
    call c_f_pointer(cdata%data_ptr, f_cdata_ptr, [cdata%dims(1)])
    f_cdata => f_cdata_ptr
  else
    f_cdata => null()
  endif
  ! in: f_calc_ok 0D_NOT_logical
  f_calc_ok = calc_ok
  ! in: f_opt_dump_spectrum 0D_NOT_integer
  if (c_associated(opt_dump_spectrum)) then
    call c_f_pointer(opt_dump_spectrum, f_opt_dump_spectrum_ptr)
  else
    f_opt_dump_spectrum_ptr => null()
  endif
  ! in: f_opt_dump_index 0D_NOT_integer
  if (c_associated(opt_dump_index)) then
    call c_f_pointer(opt_dump_index, f_opt_dump_index_ptr)
  else
    f_opt_dump_index_ptr => null()
  endif
  f_this_fft = interpolated_fft_gsl(f_cdata, f_calc_ok, f_opt_dump_spectrum_ptr, &
      f_opt_dump_index_ptr)

  ! out: f_this_fft 0D_NOT_real
  call c_f_pointer(this_fft, f_this_fft_ptr)
  f_this_fft_ptr = f_this_fft
end subroutine
subroutine fortran_inverse_prob (val, prob) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: val  ! 0D_NOT_real
  real(rp) :: f_val
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: prob  ! 0D_NOT_real
  real(rp) :: f_prob
  real(c_double), pointer :: f_prob_ptr
  ! ** End of parameters **
  ! in: f_val 0D_NOT_real
  f_val = val
  f_prob = inverse_prob(f_val)

  ! out: f_prob 0D_NOT_real
  call c_f_pointer(prob, f_prob_ptr)
  f_prob_ptr = f_prob
end subroutine
subroutine fortran_is_alphabetic (string, valid_chars, is_alpha) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: valid_chars
  character(len=4096), target :: f_valid_chars
  character(kind=c_char), pointer :: f_valid_chars_ptr(:)
  character(len=4096), pointer :: f_valid_chars_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_alpha  ! 0D_NOT_logical
  logical :: f_is_alpha
  logical(c_bool), pointer :: f_is_alpha_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_valid_chars 0D_NOT_character
  if (c_associated(valid_chars)) then
    call c_f_pointer(valid_chars, f_valid_chars_ptr, [huge(0)])
    call to_f_str(f_valid_chars_ptr, f_valid_chars)
    f_valid_chars_call_ptr => f_valid_chars
  else
    f_valid_chars_call_ptr => null()
  endif
  f_is_alpha = is_alphabetic(f_string, f_valid_chars_call_ptr)

  ! out: f_is_alpha 0D_NOT_logical
  call c_f_pointer(is_alpha, f_is_alpha_ptr)
  f_is_alpha_ptr = f_is_alpha
end subroutine
subroutine fortran_is_decreasing_sequence (array, strict, is_decreasing) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: array
  real(rp), pointer :: f_array(:)
  real(c_double), pointer :: f_array_ptr(:)
  type(c_ptr), intent(in), value :: strict  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_strict
  logical, target :: f_strict_native
  logical, pointer :: f_strict_native_ptr
  logical(c_bool), pointer :: f_strict_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_decreasing  ! 0D_NOT_logical
  logical :: f_is_decreasing
  logical(c_bool), pointer :: f_is_decreasing_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(array%data_ptr)) then
    call c_f_pointer(array%data_ptr, f_array_ptr, [array%dims(1)])
    f_array => f_array_ptr
  else
    f_array => null()
  endif
  ! in: f_strict 0D_NOT_logical
  if (c_associated(strict)) then
    call c_f_pointer(strict, f_strict_ptr)
    f_strict_native = f_strict_ptr
    f_strict_native_ptr => f_strict_native
  else
    f_strict_native_ptr => null()
  endif
  f_is_decreasing = is_decreasing_sequence(f_array, f_strict_native_ptr)

  ! out: f_is_decreasing 0D_NOT_logical
  call c_f_pointer(is_decreasing, f_is_decreasing_ptr)
  f_is_decreasing_ptr = f_is_decreasing
end subroutine
subroutine fortran_is_false (param, this_false) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: param  ! 0D_NOT_real
  real(rp) :: f_param
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_false  ! 0D_NOT_logical
  logical :: f_this_false
  logical(c_bool), pointer :: f_this_false_ptr
  ! ** End of parameters **
  ! in: f_param 0D_NOT_real
  f_param = param
  f_this_false = is_false(f_param)

  ! out: f_this_false 0D_NOT_logical
  call c_f_pointer(this_false, f_this_false_ptr)
  f_this_false_ptr = f_this_false
end subroutine
subroutine fortran_is_increasing_sequence (array, strict, is_increasing) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: array
  real(rp), pointer :: f_array(:)
  real(c_double), pointer :: f_array_ptr(:)
  type(c_ptr), intent(in), value :: strict  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_strict
  logical, target :: f_strict_native
  logical, pointer :: f_strict_native_ptr
  logical(c_bool), pointer :: f_strict_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_increasing  ! 0D_NOT_logical
  logical :: f_is_increasing
  logical(c_bool), pointer :: f_is_increasing_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(array%data_ptr)) then
    call c_f_pointer(array%data_ptr, f_array_ptr, [array%dims(1)])
    f_array => f_array_ptr
  else
    f_array => null()
  endif
  ! in: f_strict 0D_NOT_logical
  if (c_associated(strict)) then
    call c_f_pointer(strict, f_strict_ptr)
    f_strict_native = f_strict_ptr
    f_strict_native_ptr => f_strict_native
  else
    f_strict_native_ptr => null()
  endif
  f_is_increasing = is_increasing_sequence(f_array, f_strict_native_ptr)

  ! out: f_is_increasing 0D_NOT_logical
  call c_f_pointer(is_increasing, f_is_increasing_ptr)
  f_is_increasing_ptr = f_is_increasing
end subroutine
subroutine fortran_is_integer (string, int_, delims, ix_word, valid) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: int_  ! 0D_NOT_integer
  integer(c_int) :: f_int
  integer(c_int), pointer :: f_int_ptr
  type(c_ptr), intent(in), value :: delims
  character(len=4096), target :: f_delims
  character(kind=c_char), pointer :: f_delims_ptr(:)
  character(len=4096), pointer :: f_delims_call_ptr
  type(c_ptr), intent(in), value :: ix_word  ! 0D_NOT_integer
  integer(c_int) :: f_ix_word
  integer(c_int), pointer :: f_ix_word_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_int 0D_NOT_integer
  if (c_associated(int_)) then
    call c_f_pointer(int_, f_int_ptr)
  else
    f_int_ptr => null()
  endif
  ! in: f_delims 0D_NOT_character
  if (c_associated(delims)) then
    call c_f_pointer(delims, f_delims_ptr, [huge(0)])
    call to_f_str(f_delims_ptr, f_delims)
    f_delims_call_ptr => f_delims
  else
    f_delims_call_ptr => null()
  endif
  ! in: f_ix_word 0D_NOT_integer
  if (c_associated(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
  else
    f_ix_word_ptr => null()
  endif
  f_valid = is_integer(f_string, f_int_ptr, f_delims_call_ptr, f_ix_word_ptr)

  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
end subroutine
subroutine fortran_is_logical (string, ignore, valid) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: ignore  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore
  logical, target :: f_ignore_native
  logical, pointer :: f_ignore_native_ptr
  logical(c_bool), pointer :: f_ignore_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_ignore 0D_NOT_logical
  if (c_associated(ignore)) then
    call c_f_pointer(ignore, f_ignore_ptr)
    f_ignore_native = f_ignore_ptr
    f_ignore_native_ptr => f_ignore_native
  else
    f_ignore_native_ptr => null()
  endif
  f_valid = is_logical(f_string, f_ignore_native_ptr)

  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
end subroutine
subroutine fortran_is_real (string, ignore, real_num, valid) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: ignore  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore
  logical, target :: f_ignore_native
  logical, pointer :: f_ignore_native_ptr
  logical(c_bool), pointer :: f_ignore_ptr
  type(c_ptr), intent(in), value :: real_num  ! 0D_NOT_real
  real(c_double) :: f_real_num
  real(c_double), pointer :: f_real_num_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_ignore 0D_NOT_logical
  if (c_associated(ignore)) then
    call c_f_pointer(ignore, f_ignore_ptr)
    f_ignore_native = f_ignore_ptr
    f_ignore_native_ptr => f_ignore_native
  else
    f_ignore_native_ptr => null()
  endif
  ! in: f_real_num 0D_NOT_real
  if (c_associated(real_num)) then
    call c_f_pointer(real_num, f_real_num_ptr)
  else
    f_real_num_ptr => null()
  endif
  f_valid = is_real(f_string, f_ignore_native_ptr, f_real_num_ptr)

  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
end subroutine
subroutine fortran_is_subatomic_species (species, is_subatomic) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_subatomic  ! 0D_NOT_logical
  logical :: f_is_subatomic
  logical(c_bool), pointer :: f_is_subatomic_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_is_subatomic = is_subatomic_species(f_species)

  ! out: f_is_subatomic 0D_NOT_logical
  call c_f_pointer(is_subatomic, f_is_subatomic_ptr)
  f_is_subatomic_ptr = f_is_subatomic
end subroutine
subroutine fortran_is_true (param, this_true) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: param  ! 0D_NOT_real
  real(rp) :: f_param
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_true  ! 0D_NOT_logical
  logical :: f_this_true
  logical(c_bool), pointer :: f_this_true_ptr
  ! ** End of parameters **
  ! in: f_param 0D_NOT_real
  f_param = param
  f_this_true = is_true(f_param)

  ! out: f_this_true 0D_NOT_logical
  call c_f_pointer(this_true, f_this_true_ptr)
  f_this_true_ptr = f_this_true
end subroutine
subroutine fortran_j_bessel (m, arg, j_bes) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: m  ! 0D_NOT_integer
  integer :: f_m
  real(c_double) :: arg  ! 0D_NOT_real
  real(rp) :: f_arg
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: j_bes  ! 0D_NOT_real
  real(rp) :: f_j_bes
  real(c_double), pointer :: f_j_bes_ptr
  ! ** End of parameters **
  ! in: f_m 0D_NOT_integer
  f_m = m
  ! in: f_arg 0D_NOT_real
  f_arg = arg
  f_j_bes = j_bessel(f_m, f_arg)

  ! out: f_j_bes 0D_NOT_real
  call c_f_pointer(j_bes, f_j_bes_ptr)
  f_j_bes_ptr = f_j_bes
end subroutine
subroutine fortran_linear_fit (x, y, n_data, a, b, sig_a, sig_b) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: n_data  ! 0D_NOT_integer
  integer :: f_n_data
  real(c_double) :: a  ! 0D_NOT_real
  real(rp) :: f_a
  real(c_double) :: b  ! 0D_NOT_real
  real(rp) :: f_b
  real(c_double) :: sig_a  ! 0D_NOT_real
  real(rp) :: f_sig_a
  real(c_double) :: sig_b  ! 0D_NOT_real
  real(rp) :: f_sig_b
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: x
  real(rp), pointer :: f_x(:)
  real(c_double), pointer :: f_x_ptr(:)
  type(array_descriptor_t), intent(in) :: y
  real(rp), pointer :: f_y(:)
  real(c_double), pointer :: f_y_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(x%data_ptr)) then
    call c_f_pointer(x%data_ptr, f_x_ptr, [x%dims(1)])
    f_x => f_x_ptr
  else
    f_x => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(y%data_ptr)) then
    call c_f_pointer(y%data_ptr, f_y_ptr, [y%dims(1)])
    f_y => f_y_ptr
  else
    f_y => null()
  endif
  ! in: f_n_data 0D_NOT_integer
  f_n_data = n_data
  ! in: f_a 0D_NOT_real
  f_a = a
  ! in: f_b 0D_NOT_real
  f_b = b
  ! in: f_sig_a 0D_NOT_real
  f_sig_a = sig_a
  ! in: f_sig_b 0D_NOT_real
  f_sig_b = sig_b
  call linear_fit(f_x, f_y, f_n_data, f_a, f_b, f_sig_a, f_sig_b)

end subroutine
subroutine fortran_linear_fit_2d (x, y, z, coef) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: x
  real(rp), pointer :: f_x(:)
  real(c_double), pointer :: f_x_ptr(:)
  type(array_descriptor_t), intent(in) :: y
  real(rp), pointer :: f_y(:)
  real(c_double), pointer :: f_y_ptr(:)
  type(array_descriptor_t), intent(in) :: z
  real(rp), pointer :: f_z(:)
  real(c_double), pointer :: f_z_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: coef
  real(rp) :: f_coef(3)
  real(c_double), pointer :: f_coef_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(x%data_ptr)) then
    call c_f_pointer(x%data_ptr, f_x_ptr, [x%dims(1)])
    f_x => f_x_ptr
  else
    f_x => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y%data_ptr)) then
    call c_f_pointer(y%data_ptr, f_y_ptr, [y%dims(1)])
    f_y => f_y_ptr
  else
    f_y => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(z%data_ptr)) then
    call c_f_pointer(z%data_ptr, f_z_ptr, [z%dims(1)])
    f_z => f_z_ptr
  else
    f_z => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(coef%data_ptr)) then
    call c_f_pointer(coef%data_ptr, f_coef_ptr, [coef%dims(1)])
    ! output-only
  else
    f_coef_ptr => null()
  endif
  call linear_fit_2d(f_x, f_y, f_z, f_coef)

  ! out: f_coef 1D_NOT_real
  if (c_associated(coef%data_ptr)) then
    call c_f_pointer(coef%data_ptr, f_coef_ptr, [coef%dims(1)])
    f_coef_ptr = f_coef(:)
  endif
end subroutine
subroutine fortran_logic_str (logic, str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: logic  ! 0D_NOT_logical
  logical :: f_logic
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! in: f_logic 0D_NOT_logical
  f_logic = logic
  f_str = logic_str(f_logic)

  ! out: f_str 0D_NOT_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1])
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_lunget (func_retval__) bind(c)

  use array_desc_mod
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: func_retval__  ! 0D_NOT_integer
  integer :: f_func_retval__
  integer(c_int), pointer :: f_func_retval___ptr
  ! ** End of parameters **
  f_func_retval__ = lunget()

  ! out: f_func_retval__ 0D_NOT_integer
  call c_f_pointer(func_retval__, f_func_retval___ptr)
  f_func_retval___ptr = f_func_retval__
end subroutine
subroutine fortran_make_legal_comment (comment_in, comment_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: comment_in
  character(len=4096), target :: f_comment_in
  character(kind=c_char), pointer :: f_comment_in_ptr(:)
  type(c_ptr), intent(in), value :: comment_out
  character(len=4096), target :: f_comment_out
  character(kind=c_char), pointer :: f_comment_out_ptr(:)
  ! ** End of parameters **
  ! in: f_comment_in 0D_NOT_character
  if (.not. c_associated(comment_in)) return
  call c_f_pointer(comment_in, f_comment_in_ptr, [huge(0)])
  call to_f_str(f_comment_in_ptr, f_comment_in)
  ! in: f_comment_out 0D_NOT_character
  if (.not. c_associated(comment_out)) return
  call c_f_pointer(comment_out, f_comment_out_ptr, [huge(0)])
  call to_f_str(f_comment_out_ptr, f_comment_out)
  call make_legal_comment(f_comment_in, f_comment_out)

end subroutine
subroutine fortran_mass_of (species, mass) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: mass  ! 0D_NOT_real
  real(rp) :: f_mass
  real(c_double), pointer :: f_mass_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_mass = mass_of(f_species)

  ! out: f_mass 0D_NOT_real
  call c_f_pointer(mass, f_mass_ptr)
  f_mass_ptr = f_mass
end subroutine
subroutine fortran_match_reg (str, pat, is_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: pat
  character(len=4096), target :: f_pat
  character(kind=c_char), pointer :: f_pat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_match  ! 0D_NOT_logical
  logical :: f_is_match
  logical(c_bool), pointer :: f_is_match_ptr
  ! ** End of parameters **
  ! in: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! in: f_pat 0D_NOT_character
  if (.not. c_associated(pat)) return
  call c_f_pointer(pat, f_pat_ptr, [huge(0)])
  call to_f_str(f_pat_ptr, f_pat)
  f_is_match = match_reg(f_str, f_pat)

  ! out: f_is_match 0D_NOT_logical
  call c_f_pointer(is_match, f_is_match_ptr)
  f_is_match_ptr = f_is_match
end subroutine
subroutine fortran_match_wild (string, template_, is_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: template_
  character(len=4096), target :: f_template
  character(kind=c_char), pointer :: f_template_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_match  ! 0D_NOT_logical
  logical :: f_is_match
  logical(c_bool), pointer :: f_is_match_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_template 0D_NOT_character
  if (.not. c_associated(template_)) return
  call c_f_pointer(template_, f_template_ptr, [huge(0)])
  call to_f_str(f_template_ptr, f_template)
  f_is_match = match_wild(f_string, f_template)

  ! out: f_is_match 0D_NOT_logical
  call c_f_pointer(is_match, f_is_match_ptr)
  f_is_match_ptr = f_is_match
end subroutine
subroutine fortran_match_word (string, names, ix, exact_case, can_abbreviate, matched_name) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  integer(c_int) :: ix  ! 0D_NOT_integer
  integer :: f_ix
  type(c_ptr), intent(in), value :: exact_case  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact_case
  logical, target :: f_exact_case_native
  logical, pointer :: f_exact_case_native_ptr
  logical(c_bool), pointer :: f_exact_case_ptr
  type(c_ptr), intent(in), value :: can_abbreviate  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_can_abbreviate
  logical, target :: f_can_abbreviate_native
  logical, pointer :: f_can_abbreviate_native_ptr
  logical(c_bool), pointer :: f_can_abbreviate_ptr
  type(c_ptr), intent(in), value :: matched_name
  character(len=4096), target :: f_matched_name
  character(kind=c_char), pointer :: f_matched_name_ptr(:)
  character(len=4096), pointer :: f_matched_name_call_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: names
  type(character_container_alloc), pointer :: f_names
  character(200), allocatable :: f_names_local(:)
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  !! container character array (1D_NOT_character)
  if (c_associated(names))   call c_f_pointer(names, f_names)
  if (c_associated(names) .and. allocated(f_names%data)) then
    allocate(f_names_local, mold=f_names%data)
    f_names_local = ''
  endif
  ! in: f_ix 0D_NOT_integer
  f_ix = ix
  ! in: f_exact_case 0D_NOT_logical
  if (c_associated(exact_case)) then
    call c_f_pointer(exact_case, f_exact_case_ptr)
    f_exact_case_native = f_exact_case_ptr
    f_exact_case_native_ptr => f_exact_case_native
  else
    f_exact_case_native_ptr => null()
  endif
  ! in: f_can_abbreviate 0D_NOT_logical
  if (c_associated(can_abbreviate)) then
    call c_f_pointer(can_abbreviate, f_can_abbreviate_ptr)
    f_can_abbreviate_native = f_can_abbreviate_ptr
    f_can_abbreviate_native_ptr => f_can_abbreviate_native
  else
    f_can_abbreviate_native_ptr => null()
  endif
  ! in: f_matched_name 0D_NOT_character
  if (c_associated(matched_name)) then
    call c_f_pointer(matched_name, f_matched_name_ptr, [huge(0)])
    call to_f_str(f_matched_name_ptr, f_matched_name)
    f_matched_name_call_ptr => f_matched_name
  else
    f_matched_name_call_ptr => null()
  endif
  call match_word(f_string, f_names_local, f_ix, f_exact_case_native_ptr, &
      f_can_abbreviate_native_ptr, f_matched_name_call_ptr)

  !! copy allocatable character result into container
  if (c_associated(names) .and. allocated(f_names_local)) then
    if (allocated(f_names%data)) deallocate(f_names%data)
    allocate(f_names%data, source=f_names_local)
  endif
end subroutine
subroutine fortran_maximize_projection (seed, cdata, func_retval__) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: seed  ! 0D_NOT_real
  real(rp) :: f_seed
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: func_retval__  ! 0D_NOT_real
  real(rp) :: f_func_retval__
  real(c_double), pointer :: f_func_retval___ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: cdata
  complex(rp), pointer :: f_cdata(:)
  complex(c_double_complex), pointer :: f_cdata_ptr(:)
  ! ** End of parameters **
  ! in: f_seed 0D_NOT_real
  f_seed = seed
  !! general array (1D_NOT_complex) inout
  if (c_associated(cdata%data_ptr)) then
    call c_f_pointer(cdata%data_ptr, f_cdata_ptr, [cdata%dims(1)])
    f_cdata => f_cdata_ptr
  else
    f_cdata => null()
  endif
  f_func_retval__ = maximize_projection(f_seed, f_cdata)

  ! out: f_func_retval__ 0D_NOT_real
  call c_f_pointer(func_retval__, f_func_retval___ptr)
  f_func_retval___ptr = f_func_retval__
end subroutine
subroutine fortran_milli_sleep (milli_sec) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: milli_sec  ! 0D_NOT_integer
  integer :: f_milli_sec
  ! ** End of parameters **
  ! in: f_milli_sec 0D_NOT_integer
  f_milli_sec = milli_sec
  call milli_sleep(f_milli_sec)

end subroutine
subroutine fortran_modulo2_dp (x, amp, mod2) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  real(c_double) :: amp  ! 0D_NOT_real
  real(rp) :: f_amp
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: mod2  ! 0D_NOT_real
  real(rp) :: f_mod2
  real(c_double), pointer :: f_mod2_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_amp 0D_NOT_real
  f_amp = amp
  f_mod2 = modulo2_dp(f_x, f_amp)

  ! out: f_mod2 0D_NOT_real
  call c_f_pointer(mod2, f_mod2_ptr)
  f_mod2_ptr = f_mod2
end subroutine
subroutine fortran_modulo2_int (x, amp, mod2) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: x  ! 0D_NOT_integer
  integer :: f_x
  integer(c_int) :: amp  ! 0D_NOT_integer
  integer :: f_amp
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: mod2  ! 0D_NOT_integer
  integer :: f_mod2
  integer(c_int), pointer :: f_mod2_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_integer
  f_x = x
  ! in: f_amp 0D_NOT_integer
  f_amp = amp
  f_mod2 = modulo2_int(f_x, f_amp)

  ! out: f_mod2 0D_NOT_integer
  call c_f_pointer(mod2, f_mod2_ptr)
  f_mod2_ptr = f_mod2
end subroutine
subroutine fortran_modulo2_qp (x, amp, mod2) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(quadp) :: f_x
  real(c_double) :: amp  ! 0D_NOT_real
  real(quadp) :: f_amp
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: mod2  ! 0D_NOT_real
  real(quadp) :: f_mod2
  real(c_double), pointer :: f_mod2_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_amp 0D_NOT_real
  f_amp = amp
  f_mod2 = modulo2_qp(f_x, f_amp)

  ! out: f_mod2 0D_NOT_real
  call c_f_pointer(mod2, f_mod2_ptr)
  f_mod2_ptr = f_mod2
end subroutine
subroutine fortran_modulo2_sp (x, amp, mod2) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(sp) :: f_x
  real(c_double) :: amp  ! 0D_NOT_real
  real(sp) :: f_amp
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: mod2  ! 0D_NOT_real
  real(sp) :: f_mod2
  real(c_double), pointer :: f_mod2_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_amp 0D_NOT_real
  f_amp = amp
  f_mod2 = modulo2_sp(f_x, f_amp)

  ! out: f_mod2 0D_NOT_real
  call c_f_pointer(mod2, f_mod2_ptr)
  f_mod2_ptr = f_mod2
end subroutine
subroutine fortran_molecular_components (molecule, component) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: molecular_component_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: molecule
  character(len=4096), target :: f_molecule
  character(kind=c_char), pointer :: f_molecule_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: component
  type(molecular_component_struct_container_alloc), pointer :: f_component
  ! ** End of parameters **
  ! in: f_molecule 0D_NOT_character
  if (.not. c_associated(molecule)) return
  call c_f_pointer(molecule, f_molecule_ptr, [huge(0)])
  call to_f_str(f_molecule_ptr, f_molecule)
  !! container type array (1D_ALLOC_type)
  if (c_associated(component))   call c_f_pointer(component, f_component)
  call molecular_components(f_molecule, f_component%data)

end subroutine
subroutine fortran_n_bins_automatic (n_data, n) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: n_data  ! 0D_NOT_integer
  integer :: f_n_data
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: n  ! 0D_NOT_integer
  integer :: f_n
  integer(c_int), pointer :: f_n_ptr
  ! ** End of parameters **
  ! in: f_n_data 0D_NOT_integer
  f_n_data = n_data
  f_n = n_bins_automatic(f_n_data)

  ! out: f_n 0D_NOT_integer
  call c_f_pointer(n, f_n_ptr)
  f_n_ptr = f_n
end subroutine
subroutine fortran_n_choose_k (n, k, nck) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: n  ! 0D_NOT_integer
  integer :: f_n
  integer(c_int) :: k  ! 0D_NOT_integer
  integer :: f_k
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: nck  ! 0D_NOT_real
  real(rp) :: f_nck
  real(c_double), pointer :: f_nck_ptr
  ! ** End of parameters **
  ! in: f_n 0D_NOT_integer
  f_n = n
  ! in: f_k 0D_NOT_integer
  f_k = k
  f_nck = n_choose_k(f_n, f_k)

  ! out: f_nck 0D_NOT_real
  call c_f_pointer(nck, f_nck_ptr)
  f_nck_ptr = f_nck
end subroutine
subroutine fortran_n_spline_create (deriv0, deriv1, x1, n_spline) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: deriv0
  real(rp), pointer :: f_deriv0(:)
  real(c_double), pointer :: f_deriv0_ptr(:)
  type(array_descriptor_t), intent(in) :: deriv1
  real(rp), pointer :: f_deriv1(:)
  real(c_double), pointer :: f_deriv1_ptr(:)
  real(c_double) :: x1  ! 0D_NOT_real
  real(rp) :: f_x1
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: n_spline
  real(rp), pointer :: f_n_spline(:)
  real(c_double), pointer :: f_n_spline_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(deriv0%data_ptr)) then
    call c_f_pointer(deriv0%data_ptr, f_deriv0_ptr, [deriv0%dims(1)])
    f_deriv0(0:) => f_deriv0_ptr
  else
    f_deriv0 => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(deriv1%data_ptr)) then
    call c_f_pointer(deriv1%data_ptr, f_deriv1_ptr, [deriv1%dims(1)])
    f_deriv1(0:) => f_deriv1_ptr
  else
    f_deriv1 => null()
  endif
  ! in: f_x1 0D_NOT_real
  f_x1 = x1
  !! general array (1D_NOT_real) inout
  if (c_associated(n_spline%data_ptr)) then
    call c_f_pointer(n_spline%data_ptr, f_n_spline_ptr, [n_spline%dims(1)])
    f_n_spline(0:) => f_n_spline_ptr
  else
    f_n_spline => null()
  endif
  call n_spline_create(f_deriv0(0:), f_deriv1(0:), f_x1, f_n_spline(0:))

end subroutine
subroutine fortran_naff (cdata, freqs, amps, opt_dump_spectra, opt_zero_first) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: opt_dump_spectra  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_spectra
  integer(c_int), pointer :: f_opt_dump_spectra_ptr
  type(c_ptr), intent(in), value :: opt_zero_first  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_opt_zero_first
  logical, target :: f_opt_zero_first_native
  logical, pointer :: f_opt_zero_first_native_ptr
  logical(c_bool), pointer :: f_opt_zero_first_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: cdata
  complex(rp), pointer :: f_cdata(:)
  complex(c_double_complex), pointer :: f_cdata_ptr(:)
  type(array_descriptor_t), intent(in) :: freqs
  real(rp), pointer :: f_freqs(:)
  real(c_double), pointer :: f_freqs_ptr(:)
  type(array_descriptor_t), intent(in) :: amps
  complex(rp), pointer :: f_amps(:)
  complex(c_double_complex), pointer :: f_amps_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) inout
  if (c_associated(cdata%data_ptr)) then
    call c_f_pointer(cdata%data_ptr, f_cdata_ptr, [cdata%dims(1)])
    f_cdata => f_cdata_ptr
  else
    f_cdata => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(freqs%data_ptr)) then
    call c_f_pointer(freqs%data_ptr, f_freqs_ptr, [freqs%dims(1)])
    f_freqs => f_freqs_ptr
  else
    f_freqs => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(amps%data_ptr)) then
    call c_f_pointer(amps%data_ptr, f_amps_ptr, [amps%dims(1)])
    f_amps => f_amps_ptr
  else
    f_amps => null()
  endif
  ! in: f_opt_dump_spectra 0D_NOT_integer
  if (c_associated(opt_dump_spectra)) then
    call c_f_pointer(opt_dump_spectra, f_opt_dump_spectra_ptr)
  else
    f_opt_dump_spectra_ptr => null()
  endif
  ! in: f_opt_zero_first 0D_NOT_logical
  if (c_associated(opt_zero_first)) then
    call c_f_pointer(opt_zero_first, f_opt_zero_first_ptr)
    f_opt_zero_first_native = f_opt_zero_first_ptr
    f_opt_zero_first_native_ptr => f_opt_zero_first_native
  else
    f_opt_zero_first_native_ptr => null()
  endif
  call naff(f_cdata, f_freqs, f_amps, f_opt_dump_spectra_ptr, f_opt_zero_first_native_ptr)

end subroutine
subroutine fortran_nametable_add (nametable, name, ix_name) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: nametable_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  integer(c_int) :: ix_name  ! 0D_NOT_integer
  integer :: f_ix_name
  ! ** Inout parameters **
  type(c_ptr), value :: nametable  ! 0D_NOT_type
  type(nametable_struct), pointer :: f_nametable
  ! ** End of parameters **
  ! inout: f_nametable 0D_NOT_type
  if (.not. c_associated(nametable)) return
  call c_f_pointer(nametable, f_nametable)
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_ix_name 0D_NOT_integer
  f_ix_name = ix_name
  call nametable_add(f_nametable, f_name, f_ix_name)

end subroutine
subroutine fortran_nametable_bracket_indexx (nametable, name, n_match, ix_max) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: nametable_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  type(c_ptr), intent(in), value :: n_match  ! 0D_NOT_integer
  integer(c_int) :: f_n_match
  integer(c_int), pointer :: f_n_match_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_max  ! 0D_NOT_integer
  integer :: f_ix_max
  integer(c_int), pointer :: f_ix_max_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: nametable  ! 0D_NOT_type
  type(nametable_struct), pointer :: f_nametable
  ! ** End of parameters **
  ! inout: f_nametable 0D_NOT_type
  if (.not. c_associated(nametable)) return
  call c_f_pointer(nametable, f_nametable)
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_n_match 0D_NOT_integer
  if (c_associated(n_match)) then
    call c_f_pointer(n_match, f_n_match_ptr)
  else
    f_n_match_ptr => null()
  endif
  f_ix_max = nametable_bracket_indexx(f_nametable, f_name, f_n_match_ptr)

  ! out: f_ix_max 0D_NOT_integer
  call c_f_pointer(ix_max, f_ix_max_ptr)
  f_ix_max_ptr = f_ix_max
end subroutine
subroutine fortran_nametable_change1 (nametable, name, ix_name) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: nametable_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  integer(c_int) :: ix_name  ! 0D_NOT_integer
  integer :: f_ix_name
  ! ** Inout parameters **
  type(c_ptr), value :: nametable  ! 0D_NOT_type
  type(nametable_struct), pointer :: f_nametable
  ! ** End of parameters **
  ! inout: f_nametable 0D_NOT_type
  if (.not. c_associated(nametable)) return
  call c_f_pointer(nametable, f_nametable)
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_ix_name 0D_NOT_integer
  f_ix_name = ix_name
  call nametable_change1(f_nametable, f_name, f_ix_name)

end subroutine
subroutine fortran_nametable_init (nametable, n_min, n_max) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: nametable_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: n_min  ! 0D_NOT_integer
  integer(c_int) :: f_n_min
  integer(c_int), pointer :: f_n_min_ptr
  type(c_ptr), intent(in), value :: n_max  ! 0D_NOT_integer
  integer(c_int) :: f_n_max
  integer(c_int), pointer :: f_n_max_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: nametable  ! 0D_NOT_type
  type(nametable_struct), pointer :: f_nametable
  ! ** End of parameters **
  ! inout: f_nametable 0D_NOT_type
  if (.not. c_associated(nametable)) return
  call c_f_pointer(nametable, f_nametable)
  ! in: f_n_min 0D_NOT_integer
  if (c_associated(n_min)) then
    call c_f_pointer(n_min, f_n_min_ptr)
  else
    f_n_min_ptr => null()
  endif
  ! in: f_n_max 0D_NOT_integer
  if (c_associated(n_max)) then
    call c_f_pointer(n_max, f_n_max_ptr)
  else
    f_n_max_ptr => null()
  endif
  call nametable_init(f_nametable, f_n_min_ptr, f_n_max_ptr)

end subroutine
subroutine fortran_nametable_remove (nametable, ix_name) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: nametable_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_name  ! 0D_NOT_integer
  integer :: f_ix_name
  ! ** Inout parameters **
  type(c_ptr), value :: nametable  ! 0D_NOT_type
  type(nametable_struct), pointer :: f_nametable
  ! ** End of parameters **
  ! inout: f_nametable 0D_NOT_type
  if (.not. c_associated(nametable)) return
  call c_f_pointer(nametable, f_nametable)
  ! in: f_ix_name 0D_NOT_integer
  f_ix_name = ix_name
  call nametable_remove(f_nametable, f_ix_name)

end subroutine
subroutine fortran_omega_to_quat (omega, quat) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: omega
  real(rp) :: f_omega(3)
  real(c_double), pointer :: f_omega_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(omega%data_ptr)) then
    call c_f_pointer(omega%data_ptr, f_omega_ptr, [omega%dims(1)])
    f_omega = f_omega_ptr(:)
  else
    f_omega_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    ! output-only
  else
    f_quat_ptr => null()
  endif
  f_quat = omega_to_quat(f_omega)

  ! out: f_quat 1D_NOT_real
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat_ptr = f_quat(:)
  endif
end subroutine
subroutine fortran_openpmd_species_name (species, pmd_name) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: pmd_name
  character(len=4096), target :: f_pmd_name
  character(kind=c_char), pointer :: f_pmd_name_ptr(:)
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_pmd_name = openpmd_species_name(f_species)

  ! out: f_pmd_name 0D_NOT_character
  call c_f_pointer(pmd_name, f_pmd_name_ptr, [len_trim(f_pmd_name) + 1])
  call to_c_str(f_pmd_name, f_pmd_name_ptr)
end subroutine
subroutine fortran_ordinal_str (n, str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: n  ! 0D_NOT_integer
  integer :: f_n
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! in: f_n 0D_NOT_integer
  f_n = n
  f_str = ordinal_str(f_n)

  ! out: f_str 0D_ALLOC_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1])
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_out_io_buffer_get_line (ix_line, line) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_line  ! 0D_NOT_integer
  integer :: f_ix_line
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  ! ** End of parameters **
  ! in: f_ix_line 0D_NOT_integer
  f_ix_line = ix_line
  f_line = out_io_buffer_get_line(f_ix_line)

  ! out: f_line 0D_NOT_character
  call c_f_pointer(line, f_line_ptr, [len_trim(f_line) + 1])
  call to_c_str(f_line, f_line_ptr)
end subroutine
subroutine fortran_out_io_buffer_num_lines (n_lines) bind(c)

  use array_desc_mod
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: n_lines  ! 0D_NOT_integer
  integer :: f_n_lines
  integer(c_int), pointer :: f_n_lines_ptr
  ! ** End of parameters **
  f_n_lines = out_io_buffer_num_lines()

  ! out: f_n_lines 0D_NOT_integer
  call c_f_pointer(n_lines, f_n_lines_ptr)
  f_n_lines_ptr = f_n_lines
end subroutine
subroutine fortran_out_io_buffer_reset () bind(c)

  use array_desc_mod
  implicit none
  ! ** End of parameters **
  call out_io_buffer_reset()

end subroutine
subroutine fortran_out_io_int (level, routine_name, line, i_num, insert_tag_line) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: level  ! 0D_NOT_integer
  integer :: f_level
  type(c_ptr), intent(in), value :: routine_name
  character(len=4096), target :: f_routine_name
  character(kind=c_char), pointer :: f_routine_name_ptr(:)
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  integer(c_int) :: i_num  ! 0D_NOT_integer
  integer :: f_i_num
  type(c_ptr), intent(in), value :: insert_tag_line  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_insert_tag_line
  logical, target :: f_insert_tag_line_native
  logical, pointer :: f_insert_tag_line_native_ptr
  logical(c_bool), pointer :: f_insert_tag_line_ptr
  ! ** End of parameters **
  ! in: f_level 0D_NOT_integer
  f_level = level
  ! in: f_routine_name 0D_NOT_character
  if (.not. c_associated(routine_name)) return
  call c_f_pointer(routine_name, f_routine_name_ptr, [huge(0)])
  call to_f_str(f_routine_name_ptr, f_routine_name)
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_i_num 0D_NOT_integer
  f_i_num = i_num
  ! in: f_insert_tag_line 0D_NOT_logical
  if (c_associated(insert_tag_line)) then
    call c_f_pointer(insert_tag_line, f_insert_tag_line_ptr)
    f_insert_tag_line_native = f_insert_tag_line_ptr
    f_insert_tag_line_native_ptr => f_insert_tag_line_native
  else
    f_insert_tag_line_native_ptr => null()
  endif
  call out_io(f_level, f_routine_name, f_line, f_i_num, f_insert_tag_line_native_ptr)

end subroutine
subroutine fortran_out_io_line12 (level, routine_name, line1, line2, line3, line4, line5, &
    line6, line7, line8, line9, line10, line11, line12, r_array, i_array, l_array, &
    insert_tag_line) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: level  ! 0D_NOT_integer
  integer :: f_level
  type(c_ptr), intent(in), value :: routine_name
  character(len=4096), target :: f_routine_name
  character(kind=c_char), pointer :: f_routine_name_ptr(:)
  type(c_ptr), intent(in), value :: line1
  character(len=4096), target :: f_line1
  character(kind=c_char), pointer :: f_line1_ptr(:)
  type(c_ptr), intent(in), value :: line2
  character(len=4096), target :: f_line2
  character(kind=c_char), pointer :: f_line2_ptr(:)
  character(len=4096), pointer :: f_line2_call_ptr
  type(c_ptr), intent(in), value :: line3
  character(len=4096), target :: f_line3
  character(kind=c_char), pointer :: f_line3_ptr(:)
  character(len=4096), pointer :: f_line3_call_ptr
  type(c_ptr), intent(in), value :: line4
  character(len=4096), target :: f_line4
  character(kind=c_char), pointer :: f_line4_ptr(:)
  character(len=4096), pointer :: f_line4_call_ptr
  type(c_ptr), intent(in), value :: line5
  character(len=4096), target :: f_line5
  character(kind=c_char), pointer :: f_line5_ptr(:)
  character(len=4096), pointer :: f_line5_call_ptr
  type(c_ptr), intent(in), value :: line6
  character(len=4096), target :: f_line6
  character(kind=c_char), pointer :: f_line6_ptr(:)
  character(len=4096), pointer :: f_line6_call_ptr
  type(c_ptr), intent(in), value :: line7
  character(len=4096), target :: f_line7
  character(kind=c_char), pointer :: f_line7_ptr(:)
  character(len=4096), pointer :: f_line7_call_ptr
  type(c_ptr), intent(in), value :: line8
  character(len=4096), target :: f_line8
  character(kind=c_char), pointer :: f_line8_ptr(:)
  character(len=4096), pointer :: f_line8_call_ptr
  type(c_ptr), intent(in), value :: line9
  character(len=4096), target :: f_line9
  character(kind=c_char), pointer :: f_line9_ptr(:)
  character(len=4096), pointer :: f_line9_call_ptr
  type(c_ptr), intent(in), value :: line10
  character(len=4096), target :: f_line10
  character(kind=c_char), pointer :: f_line10_ptr(:)
  character(len=4096), pointer :: f_line10_call_ptr
  type(c_ptr), intent(in), value :: line11
  character(len=4096), target :: f_line11
  character(kind=c_char), pointer :: f_line11_ptr(:)
  character(len=4096), pointer :: f_line11_call_ptr
  type(c_ptr), intent(in), value :: line12
  character(len=4096), target :: f_line12
  character(kind=c_char), pointer :: f_line12_ptr(:)
  character(len=4096), pointer :: f_line12_call_ptr
  type(c_ptr), intent(in), value :: insert_tag_line  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_insert_tag_line
  logical, target :: f_insert_tag_line_native
  logical, pointer :: f_insert_tag_line_native_ptr
  logical(c_bool), pointer :: f_insert_tag_line_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: r_array
  real(rp), pointer :: f_r_array(:)
  real(c_double), pointer :: f_r_array_ptr(:)
  type(array_descriptor_t), intent(in) :: i_array
  integer, pointer :: f_i_array(:)
  integer(c_int), pointer :: f_i_array_ptr(:)
  type(c_ptr), intent(in), value :: l_array
  type(logical_container_alloc), pointer :: f_l_array
  ! ** End of parameters **
  ! in: f_level 0D_NOT_integer
  f_level = level
  ! in: f_routine_name 0D_NOT_character
  if (.not. c_associated(routine_name)) return
  call c_f_pointer(routine_name, f_routine_name_ptr, [huge(0)])
  call to_f_str(f_routine_name_ptr, f_routine_name)
  ! in: f_line1 0D_NOT_character
  if (.not. c_associated(line1)) return
  call c_f_pointer(line1, f_line1_ptr, [huge(0)])
  call to_f_str(f_line1_ptr, f_line1)
  ! in: f_line2 0D_NOT_character
  if (c_associated(line2)) then
    call c_f_pointer(line2, f_line2_ptr, [huge(0)])
    call to_f_str(f_line2_ptr, f_line2)
    f_line2_call_ptr => f_line2
  else
    f_line2_call_ptr => null()
  endif
  ! in: f_line3 0D_NOT_character
  if (c_associated(line3)) then
    call c_f_pointer(line3, f_line3_ptr, [huge(0)])
    call to_f_str(f_line3_ptr, f_line3)
    f_line3_call_ptr => f_line3
  else
    f_line3_call_ptr => null()
  endif
  ! in: f_line4 0D_NOT_character
  if (c_associated(line4)) then
    call c_f_pointer(line4, f_line4_ptr, [huge(0)])
    call to_f_str(f_line4_ptr, f_line4)
    f_line4_call_ptr => f_line4
  else
    f_line4_call_ptr => null()
  endif
  ! in: f_line5 0D_NOT_character
  if (c_associated(line5)) then
    call c_f_pointer(line5, f_line5_ptr, [huge(0)])
    call to_f_str(f_line5_ptr, f_line5)
    f_line5_call_ptr => f_line5
  else
    f_line5_call_ptr => null()
  endif
  ! in: f_line6 0D_NOT_character
  if (c_associated(line6)) then
    call c_f_pointer(line6, f_line6_ptr, [huge(0)])
    call to_f_str(f_line6_ptr, f_line6)
    f_line6_call_ptr => f_line6
  else
    f_line6_call_ptr => null()
  endif
  ! in: f_line7 0D_NOT_character
  if (c_associated(line7)) then
    call c_f_pointer(line7, f_line7_ptr, [huge(0)])
    call to_f_str(f_line7_ptr, f_line7)
    f_line7_call_ptr => f_line7
  else
    f_line7_call_ptr => null()
  endif
  ! in: f_line8 0D_NOT_character
  if (c_associated(line8)) then
    call c_f_pointer(line8, f_line8_ptr, [huge(0)])
    call to_f_str(f_line8_ptr, f_line8)
    f_line8_call_ptr => f_line8
  else
    f_line8_call_ptr => null()
  endif
  ! in: f_line9 0D_NOT_character
  if (c_associated(line9)) then
    call c_f_pointer(line9, f_line9_ptr, [huge(0)])
    call to_f_str(f_line9_ptr, f_line9)
    f_line9_call_ptr => f_line9
  else
    f_line9_call_ptr => null()
  endif
  ! in: f_line10 0D_NOT_character
  if (c_associated(line10)) then
    call c_f_pointer(line10, f_line10_ptr, [huge(0)])
    call to_f_str(f_line10_ptr, f_line10)
    f_line10_call_ptr => f_line10
  else
    f_line10_call_ptr => null()
  endif
  ! in: f_line11 0D_NOT_character
  if (c_associated(line11)) then
    call c_f_pointer(line11, f_line11_ptr, [huge(0)])
    call to_f_str(f_line11_ptr, f_line11)
    f_line11_call_ptr => f_line11
  else
    f_line11_call_ptr => null()
  endif
  ! in: f_line12 0D_NOT_character
  if (c_associated(line12)) then
    call c_f_pointer(line12, f_line12_ptr, [huge(0)])
    call to_f_str(f_line12_ptr, f_line12)
    f_line12_call_ptr => f_line12
  else
    f_line12_call_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(r_array%data_ptr)) then
    call c_f_pointer(r_array%data_ptr, f_r_array_ptr, [r_array%dims(1)])
    f_r_array => f_r_array_ptr
  else
    f_r_array => null()
  endif
  !! general array (1D_NOT_integer) inout
  if (c_associated(i_array%data_ptr)) then
    call c_f_pointer(i_array%data_ptr, f_i_array_ptr, [i_array%dims(1)])
    f_i_array => f_i_array_ptr
  else
    f_i_array => null()
  endif
  !! container general array (1D_ALLOC_logical)
  if (c_associated(l_array))   call c_f_pointer(l_array, f_l_array)
  ! in: f_insert_tag_line 0D_NOT_logical
  if (c_associated(insert_tag_line)) then
    call c_f_pointer(insert_tag_line, f_insert_tag_line_ptr)
    f_insert_tag_line_native = f_insert_tag_line_ptr
    f_insert_tag_line_native_ptr => f_insert_tag_line_native
  else
    f_insert_tag_line_native_ptr => null()
  endif
  call out_io(f_level, f_routine_name, f_line1, f_line2_call_ptr, f_line3_call_ptr, &
      f_line4_call_ptr, f_line5_call_ptr, f_line6_call_ptr, f_line7_call_ptr, f_line8_call_ptr, &
      f_line9_call_ptr, f_line10_call_ptr, f_line11_call_ptr, f_line12_call_ptr, f_r_array, &
      f_i_array, f_l_array%data, f_insert_tag_line_native_ptr)

end subroutine
subroutine fortran_out_io_lines (level, routine_name, lines, r_array, i_array, l_array, &
    insert_tag_line) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: level  ! 0D_NOT_integer
  integer :: f_level
  type(c_ptr), intent(in), value :: routine_name
  character(len=4096), target :: f_routine_name
  character(kind=c_char), pointer :: f_routine_name_ptr(:)
  type(c_ptr), intent(in), value :: insert_tag_line  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_insert_tag_line
  logical, target :: f_insert_tag_line_native
  logical, pointer :: f_insert_tag_line_native_ptr
  logical(c_bool), pointer :: f_insert_tag_line_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: lines
  type(character_container_alloc), pointer :: f_lines
  character(200), allocatable :: f_lines_local(:)
  type(array_descriptor_t), intent(in) :: r_array
  real(rp), pointer :: f_r_array(:)
  real(c_double), pointer :: f_r_array_ptr(:)
  type(array_descriptor_t), intent(in) :: i_array
  integer, pointer :: f_i_array(:)
  integer(c_int), pointer :: f_i_array_ptr(:)
  type(c_ptr), intent(in), value :: l_array
  type(logical_container_alloc), pointer :: f_l_array
  ! ** End of parameters **
  ! in: f_level 0D_NOT_integer
  f_level = level
  ! in: f_routine_name 0D_NOT_character
  if (.not. c_associated(routine_name)) return
  call c_f_pointer(routine_name, f_routine_name_ptr, [huge(0)])
  call to_f_str(f_routine_name_ptr, f_routine_name)
  !! container character array (1D_NOT_character)
  if (c_associated(lines))   call c_f_pointer(lines, f_lines)
  if (c_associated(lines) .and. allocated(f_lines%data)) then
    allocate(f_lines_local, mold=f_lines%data)
    f_lines_local = ''
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(r_array%data_ptr)) then
    call c_f_pointer(r_array%data_ptr, f_r_array_ptr, [r_array%dims(1)])
    f_r_array => f_r_array_ptr
  else
    f_r_array => null()
  endif
  !! general array (1D_NOT_integer) inout
  if (c_associated(i_array%data_ptr)) then
    call c_f_pointer(i_array%data_ptr, f_i_array_ptr, [i_array%dims(1)])
    f_i_array => f_i_array_ptr
  else
    f_i_array => null()
  endif
  !! container general array (1D_ALLOC_logical)
  if (c_associated(l_array))   call c_f_pointer(l_array, f_l_array)
  ! in: f_insert_tag_line 0D_NOT_logical
  if (c_associated(insert_tag_line)) then
    call c_f_pointer(insert_tag_line, f_insert_tag_line_ptr)
    f_insert_tag_line_native = f_insert_tag_line_ptr
    f_insert_tag_line_native_ptr => f_insert_tag_line_native
  else
    f_insert_tag_line_native_ptr => null()
  endif
  call out_io(f_level, f_routine_name, f_lines_local, f_r_array, f_i_array, f_l_array%data, &
      f_insert_tag_line_native_ptr)

  !! copy allocatable character result into container
  if (c_associated(lines) .and. allocated(f_lines_local)) then
    if (allocated(f_lines%data)) deallocate(f_lines%data)
    allocate(f_lines%data, source=f_lines_local)
  endif
end subroutine
subroutine fortran_out_io_logical (level, routine_name, line, l_num, insert_tag_line) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: level  ! 0D_NOT_integer
  integer :: f_level
  type(c_ptr), intent(in), value :: routine_name
  character(len=4096), target :: f_routine_name
  character(kind=c_char), pointer :: f_routine_name_ptr(:)
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  logical(c_bool) :: l_num  ! 0D_NOT_logical
  logical :: f_l_num
  type(c_ptr), intent(in), value :: insert_tag_line  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_insert_tag_line
  logical, target :: f_insert_tag_line_native
  logical, pointer :: f_insert_tag_line_native_ptr
  logical(c_bool), pointer :: f_insert_tag_line_ptr
  ! ** End of parameters **
  ! in: f_level 0D_NOT_integer
  f_level = level
  ! in: f_routine_name 0D_NOT_character
  if (.not. c_associated(routine_name)) return
  call c_f_pointer(routine_name, f_routine_name_ptr, [huge(0)])
  call to_f_str(f_routine_name_ptr, f_routine_name)
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_l_num 0D_NOT_logical
  f_l_num = l_num
  ! in: f_insert_tag_line 0D_NOT_logical
  if (c_associated(insert_tag_line)) then
    call c_f_pointer(insert_tag_line, f_insert_tag_line_ptr)
    f_insert_tag_line_native = f_insert_tag_line_ptr
    f_insert_tag_line_native_ptr => f_insert_tag_line_native
  else
    f_insert_tag_line_native_ptr => null()
  endif
  call out_io(f_level, f_routine_name, f_line, f_l_num, f_insert_tag_line_native_ptr)

end subroutine
subroutine fortran_out_io_print_and_capture_setup (print_on, capture_state, capture_add_null) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: print_on  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_on
  logical, target :: f_print_on_native
  logical, pointer :: f_print_on_native_ptr
  logical(c_bool), pointer :: f_print_on_ptr
  type(c_ptr), intent(in), value :: capture_state
  character(len=4096), target :: f_capture_state
  character(kind=c_char), pointer :: f_capture_state_ptr(:)
  character(len=4096), pointer :: f_capture_state_call_ptr
  type(c_ptr), intent(in), value :: capture_add_null  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_capture_add_null
  logical, target :: f_capture_add_null_native
  logical, pointer :: f_capture_add_null_native_ptr
  logical(c_bool), pointer :: f_capture_add_null_ptr
  ! ** End of parameters **
  ! in: f_print_on 0D_NOT_logical
  if (c_associated(print_on)) then
    call c_f_pointer(print_on, f_print_on_ptr)
    f_print_on_native = f_print_on_ptr
    f_print_on_native_ptr => f_print_on_native
  else
    f_print_on_native_ptr => null()
  endif
  ! in: f_capture_state 0D_NOT_character
  if (c_associated(capture_state)) then
    call c_f_pointer(capture_state, f_capture_state_ptr, [huge(0)])
    call to_f_str(f_capture_state_ptr, f_capture_state)
    f_capture_state_call_ptr => f_capture_state
  else
    f_capture_state_call_ptr => null()
  endif
  ! in: f_capture_add_null 0D_NOT_logical
  if (c_associated(capture_add_null)) then
    call c_f_pointer(capture_add_null, f_capture_add_null_ptr)
    f_capture_add_null_native = f_capture_add_null_ptr
    f_capture_add_null_native_ptr => f_capture_add_null_native
  else
    f_capture_add_null_native_ptr => null()
  endif
  call out_io_print_and_capture_setup(f_print_on_native_ptr, f_capture_state_call_ptr, &
      f_capture_add_null_native_ptr)

end subroutine
subroutine fortran_out_io_real (level, routine_name, line, r_num, insert_tag_line) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: level  ! 0D_NOT_integer
  integer :: f_level
  type(c_ptr), intent(in), value :: routine_name
  character(len=4096), target :: f_routine_name
  character(kind=c_char), pointer :: f_routine_name_ptr(:)
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  real(c_double) :: r_num  ! 0D_NOT_real
  real(rp) :: f_r_num
  type(c_ptr), intent(in), value :: insert_tag_line  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_insert_tag_line
  logical, target :: f_insert_tag_line_native
  logical, pointer :: f_insert_tag_line_native_ptr
  logical(c_bool), pointer :: f_insert_tag_line_ptr
  ! ** End of parameters **
  ! in: f_level 0D_NOT_integer
  f_level = level
  ! in: f_routine_name 0D_NOT_character
  if (.not. c_associated(routine_name)) return
  call c_f_pointer(routine_name, f_routine_name_ptr, [huge(0)])
  call to_f_str(f_routine_name_ptr, f_routine_name)
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_r_num 0D_NOT_real
  f_r_num = r_num
  ! in: f_insert_tag_line 0D_NOT_logical
  if (c_associated(insert_tag_line)) then
    call c_f_pointer(insert_tag_line, f_insert_tag_line_ptr)
    f_insert_tag_line_native = f_insert_tag_line_ptr
    f_insert_tag_line_native_ptr => f_insert_tag_line_native
  else
    f_insert_tag_line_native_ptr => null()
  endif
  call out_io(f_level, f_routine_name, f_line, f_r_num, f_insert_tag_line_native_ptr)

end subroutine
subroutine fortran_output_direct (file_unit, print_and_capture, min_level, max_level, set, get) &
    bind(c)

  use array_desc_mod
  use output_mod, only: out_io_output_direct_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_unit  ! 0D_NOT_integer
  integer(c_int) :: f_file_unit
  integer(c_int), pointer :: f_file_unit_ptr
  type(c_ptr), intent(in), value :: print_and_capture  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_and_capture
  logical, target :: f_print_and_capture_native
  logical, pointer :: f_print_and_capture_native_ptr
  logical(c_bool), pointer :: f_print_and_capture_ptr
  type(c_ptr), intent(in), value :: min_level  ! 0D_NOT_integer
  integer(c_int) :: f_min_level
  integer(c_int), pointer :: f_min_level_ptr
  type(c_ptr), intent(in), value :: max_level  ! 0D_NOT_integer
  integer(c_int) :: f_max_level
  integer(c_int), pointer :: f_max_level_ptr
  type(c_ptr), value :: set  ! 0D_NOT_type
  type(out_io_output_direct_struct), pointer :: f_set
  ! ** Out parameters **
  type(c_ptr), value :: get  ! 0D_NOT_type
  type(out_io_output_direct_struct), pointer :: f_get
  ! ** End of parameters **
  ! in: f_file_unit 0D_NOT_integer
  if (c_associated(file_unit)) then
    call c_f_pointer(file_unit, f_file_unit_ptr)
  else
    f_file_unit_ptr => null()
  endif
  ! in: f_print_and_capture 0D_NOT_logical
  if (c_associated(print_and_capture)) then
    call c_f_pointer(print_and_capture, f_print_and_capture_ptr)
    f_print_and_capture_native = f_print_and_capture_ptr
    f_print_and_capture_native_ptr => f_print_and_capture_native
  else
    f_print_and_capture_native_ptr => null()
  endif
  ! in: f_min_level 0D_NOT_integer
  if (c_associated(min_level)) then
    call c_f_pointer(min_level, f_min_level_ptr)
  else
    f_min_level_ptr => null()
  endif
  ! in: f_max_level 0D_NOT_integer
  if (c_associated(max_level)) then
    call c_f_pointer(max_level, f_max_level_ptr)
  else
    f_max_level_ptr => null()
  endif
  ! in: f_set 0D_NOT_type
  if (c_associated(set))   call c_f_pointer(set, f_set)
  ! out: f_get 0D_NOT_type
  if (c_associated(get))   call c_f_pointer(get, f_get)
  call output_direct(f_file_unit_ptr, f_print_and_capture_native_ptr, f_min_level_ptr, &
      f_max_level_ptr, f_set, f_get)

  ! out: f_get 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_parse_fortran_format (format_str, n_repeat, power, descrip, width, digits) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: format_str
  character(len=4096), target :: f_format_str
  character(kind=c_char), pointer :: f_format_str_ptr(:)
  integer(c_int) :: n_repeat  ! 0D_NOT_integer
  integer :: f_n_repeat
  integer(c_int) :: power  ! 0D_NOT_integer
  integer :: f_power
  type(c_ptr), intent(in), value :: descrip
  character(len=4096), target :: f_descrip
  character(kind=c_char), pointer :: f_descrip_ptr(:)
  integer(c_int) :: width  ! 0D_NOT_integer
  integer :: f_width
  integer(c_int) :: digits  ! 0D_NOT_integer
  integer :: f_digits
  ! ** End of parameters **
  ! in: f_format_str 0D_NOT_character
  if (.not. c_associated(format_str)) return
  call c_f_pointer(format_str, f_format_str_ptr, [huge(0)])
  call to_f_str(f_format_str_ptr, f_format_str)
  ! in: f_n_repeat 0D_NOT_integer
  f_n_repeat = n_repeat
  ! in: f_power 0D_NOT_integer
  f_power = power
  ! in: f_descrip 0D_NOT_character
  if (.not. c_associated(descrip)) return
  call c_f_pointer(descrip, f_descrip_ptr, [huge(0)])
  call to_f_str(f_descrip_ptr, f_descrip)
  ! in: f_width 0D_NOT_integer
  f_width = width
  ! in: f_digits 0D_NOT_integer
  f_digits = digits
  call parse_fortran_format(f_format_str, f_n_repeat, f_power, f_descrip, f_width, f_digits)

end subroutine
subroutine fortran_pointer_to_locations (string, array, num, ix_min, ix_max, names, exact_case, &
    print_err) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  integer(c_int) :: num  ! 0D_NOT_integer
  integer :: f_num
  integer(c_int) :: ix_min  ! 0D_NOT_integer
  integer :: f_ix_min
  integer(c_int) :: ix_max  ! 0D_NOT_integer
  integer :: f_ix_max
  type(c_ptr), intent(in), value :: exact_case  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact_case
  logical, target :: f_exact_case_native
  logical, pointer :: f_exact_case_native_ptr
  logical(c_bool), pointer :: f_exact_case_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: array
  type(integer_container_alloc), pointer :: f_array
  type(c_ptr), intent(in), value :: names
  type(character_container_alloc), pointer :: f_names
  character(200), allocatable :: f_names_local(:)
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(array))   call c_f_pointer(array, f_array)
  ! in: f_num 0D_NOT_integer
  f_num = num
  ! in: f_ix_min 0D_NOT_integer
  f_ix_min = ix_min
  ! in: f_ix_max 0D_NOT_integer
  f_ix_max = ix_max
  !! container character array (1D_NOT_character)
  if (c_associated(names))   call c_f_pointer(names, f_names)
  if (c_associated(names) .and. allocated(f_names%data)) then
    allocate(f_names_local, mold=f_names%data)
    f_names_local = ''
  endif
  ! in: f_exact_case 0D_NOT_logical
  if (c_associated(exact_case)) then
    call c_f_pointer(exact_case, f_exact_case_ptr)
    f_exact_case_native = f_exact_case_ptr
    f_exact_case_native_ptr => f_exact_case_native
  else
    f_exact_case_native_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call pointer_to_locations(f_string, f_array%data, f_num, f_ix_min, f_ix_max, f_names_local, &
      f_exact_case_native_ptr, f_print_err_native_ptr)

  !! copy allocatable character result into container
  if (c_associated(names) .and. allocated(f_names_local)) then
    if (allocated(f_names%data)) deallocate(f_names%data)
    allocate(f_names%data, source=f_names_local)
  endif
end subroutine
subroutine fortran_pointer_to_ran_state (ran_state, ix_thread, ran_state_ptr) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  type(c_ptr), intent(in), value :: ix_thread  ! 0D_NOT_integer
  integer(c_int) :: f_ix_thread
  integer(c_int), pointer :: f_ix_thread_ptr
  ! ** Out parameters **
  type(c_ptr) :: ran_state_ptr  ! 0D_PTR_type
  type(random_state_struct), pointer :: f_ran_state_ptr
  ! ** End of parameters **
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  ! in: f_ix_thread 0D_NOT_integer
  if (c_associated(ix_thread)) then
    call c_f_pointer(ix_thread, f_ix_thread_ptr)
  else
    f_ix_thread_ptr => null()
  endif
  ! out: f_ran_state_ptr 0D_PTR_type
  if (c_associated(ran_state_ptr))   call c_f_pointer(ran_state_ptr, f_ran_state_ptr)
  f_ran_state_ptr => pointer_to_ran_state(f_ran_state, f_ix_thread_ptr)

  ! out: f_ran_state_ptr 0D_PTR_type
  ran_state_ptr = c_loc(f_ran_state_ptr)
end subroutine
subroutine fortran_poly_eval (poly, x, diff_coef, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: poly
  real(rp), pointer :: f_poly(:)
  real(c_double), pointer :: f_poly_ptr(:)
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: diff_coef  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_diff_coef
  logical, target :: f_diff_coef_native
  logical, pointer :: f_diff_coef_native_ptr
  logical(c_bool), pointer :: f_diff_coef_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(poly%data_ptr)) then
    call c_f_pointer(poly%data_ptr, f_poly_ptr, [poly%dims(1)])
    f_poly(0:) => f_poly_ptr
  else
    f_poly => null()
  endif
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_diff_coef 0D_NOT_logical
  if (c_associated(diff_coef)) then
    call c_f_pointer(diff_coef, f_diff_coef_ptr)
    f_diff_coef_native = f_diff_coef_ptr
    f_diff_coef_native_ptr => f_diff_coef_native
  else
    f_diff_coef_native_ptr => null()
  endif
  f_y = poly_eval(f_poly(0:), f_x, f_diff_coef_native_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_probability_funct (x, prob) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: prob  ! 0D_NOT_real
  real(rp) :: f_prob
  real(c_double), pointer :: f_prob_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  f_prob = probability_funct(f_x)

  ! out: f_prob 0D_NOT_real
  call c_f_pointer(prob, f_prob_ptr)
  f_prob_ptr = f_prob
end subroutine
subroutine fortran_projdd (a, b, func_retval__) bind(c)

  use array_desc_mod
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: func_retval__  ! 0D_NOT_complex
  complex(rp) :: f_func_retval__
  complex(c_double_complex), pointer :: f_func_retval___ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: a
  complex(rp), pointer :: f_a(:)
  complex(c_double_complex), pointer :: f_a_ptr(:)
  type(array_descriptor_t), intent(in) :: b
  complex(rp), pointer :: f_b(:)
  complex(c_double_complex), pointer :: f_b_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) inout
  if (c_associated(a%data_ptr)) then
    call c_f_pointer(a%data_ptr, f_a_ptr, [a%dims(1)])
    f_a => f_a_ptr
  else
    f_a => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(b%data_ptr)) then
    call c_f_pointer(b%data_ptr, f_b_ptr, [b%dims(1)])
    f_b => f_b_ptr
  else
    f_b => null()
  endif
  f_func_retval__ = projdd(f_a, f_b)

  ! out: f_func_retval__ 0D_NOT_complex
  call c_f_pointer(func_retval__, f_func_retval___ptr)
  f_func_retval___ptr = f_func_retval__
end subroutine
subroutine fortran_quadratic_roots (coefs, root) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: coefs
  real(rp) :: f_coefs(3)
  real(c_double), pointer :: f_coefs_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: root
  complex(rp) :: f_root(2)
  complex(c_double_complex), pointer :: f_root_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(coefs%data_ptr)) then
    call c_f_pointer(coefs%data_ptr, f_coefs_ptr, [coefs%dims(1)])
    f_coefs = f_coefs_ptr(:)
  else
    f_coefs_ptr => null()
  endif
  !! general array (1D_NOT_complex) out
  if (c_associated(root%data_ptr)) then
    call c_f_pointer(root%data_ptr, f_root_ptr, [root%dims(1)])
    ! output-only
  else
    f_root_ptr => null()
  endif
  f_root = quadratic_roots(f_coefs)

  ! out: f_root 1D_NOT_complex
  if (c_associated(root%data_ptr)) then
    call c_f_pointer(root%data_ptr, f_root_ptr, [root%dims(1)])
    f_root_ptr = f_root(:)
  endif
end subroutine
subroutine fortran_quat_conj_complex (q_in, q_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: q_in
  complex(rp) :: f_q_in(0:3)
  complex(c_double_complex), pointer :: f_q_in_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: q_out
  complex(rp) :: f_q_out(0:3)
  complex(c_double_complex), pointer :: f_q_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) in
  if (c_associated(q_in%data_ptr)) then
    call c_f_pointer(q_in%data_ptr, f_q_in_ptr, [q_in%dims(1)])
    f_q_in = f_q_in_ptr(:)
  else
    f_q_in_ptr => null()
  endif
  !! general array (1D_NOT_complex) out
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    ! output-only
  else
    f_q_out_ptr => null()
  endif
  f_q_out = quat_conj(f_q_in)

  ! out: f_q_out 1D_NOT_complex
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    f_q_out_ptr = f_q_out(:)
  endif
end subroutine
subroutine fortran_quat_conj_real (q_in, q_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: q_in
  real(rp) :: f_q_in(0:3)
  real(c_double), pointer :: f_q_in_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: q_out
  real(rp) :: f_q_out(0:3)
  real(c_double), pointer :: f_q_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(q_in%data_ptr)) then
    call c_f_pointer(q_in%data_ptr, f_q_in_ptr, [q_in%dims(1)])
    f_q_in = f_q_in_ptr(:)
  else
    f_q_in_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    ! output-only
  else
    f_q_out_ptr => null()
  endif
  f_q_out = quat_conj(f_q_in)

  ! out: f_q_out 1D_NOT_real
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    f_q_out_ptr = f_q_out(:)
  endif
end subroutine
subroutine fortran_quat_inverse (q_in, q_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: q_in
  real(rp) :: f_q_in(0:3)
  real(c_double), pointer :: f_q_in_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: q_out
  real(rp) :: f_q_out(0:3)
  real(c_double), pointer :: f_q_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(q_in%data_ptr)) then
    call c_f_pointer(q_in%data_ptr, f_q_in_ptr, [q_in%dims(1)])
    f_q_in = f_q_in_ptr(:)
  else
    f_q_in_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    ! output-only
  else
    f_q_out_ptr => null()
  endif
  f_q_out = quat_inverse(f_q_in)

  ! out: f_q_out 1D_NOT_real
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    f_q_out_ptr = f_q_out(:)
  endif
end subroutine
subroutine fortran_quat_mul_complex (q1, q2, q3, q4, q5, q6, q7, q8, q9, q_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: q1
  complex(rp) :: f_q1(0:3)
  complex(c_double_complex), pointer :: f_q1_ptr(:)
  type(array_descriptor_t), intent(in) :: q2
  complex(rp) :: f_q2(0:3)
  complex(c_double_complex), pointer :: f_q2_ptr(:)
  type(array_descriptor_t), intent(in) :: q3
  complex(rp) :: f_q3(0:3)
  complex(c_double_complex), pointer :: f_q3_ptr(:)
  type(array_descriptor_t), intent(in) :: q9
  complex(rp) :: f_q9(0:3)
  complex(c_double_complex), pointer :: f_q9_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: q_out
  complex(rp) :: f_q_out(0:3)
  complex(c_double_complex), pointer :: f_q_out_ptr(:)
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: q4
  complex(rp) :: f_q4(0:3)
  complex(c_double_complex), pointer :: f_q4_ptr(:)
  type(array_descriptor_t), intent(in) :: q5
  complex(rp) :: f_q5(0:3)
  complex(c_double_complex), pointer :: f_q5_ptr(:)
  type(array_descriptor_t), intent(in) :: q6
  complex(rp) :: f_q6(0:3)
  complex(c_double_complex), pointer :: f_q6_ptr(:)
  type(array_descriptor_t), intent(in) :: q7
  complex(rp) :: f_q7(0:3)
  complex(c_double_complex), pointer :: f_q7_ptr(:)
  type(array_descriptor_t), intent(in) :: q8
  complex(rp) :: f_q8(0:3)
  complex(c_double_complex), pointer :: f_q8_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) in
  if (c_associated(q1%data_ptr)) then
    call c_f_pointer(q1%data_ptr, f_q1_ptr, [q1%dims(1)])
    f_q1 = f_q1_ptr(:)
  else
    f_q1_ptr => null()
  endif
  !! general array (1D_NOT_complex) in
  if (c_associated(q2%data_ptr)) then
    call c_f_pointer(q2%data_ptr, f_q2_ptr, [q2%dims(1)])
    f_q2 = f_q2_ptr(:)
  else
    f_q2_ptr => null()
  endif
  !! general array (1D_NOT_complex) in
  if (c_associated(q3%data_ptr)) then
    call c_f_pointer(q3%data_ptr, f_q3_ptr, [q3%dims(1)])
    f_q3 = f_q3_ptr(:)
  else
    f_q3_ptr => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(q4%data_ptr)) then
    call c_f_pointer(q4%data_ptr, f_q4_ptr, [q4%dims(1)])
    f_q4 = f_q4_ptr(:)
  else
    f_q4_ptr => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(q5%data_ptr)) then
    call c_f_pointer(q5%data_ptr, f_q5_ptr, [q5%dims(1)])
    f_q5 = f_q5_ptr(:)
  else
    f_q5_ptr => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(q6%data_ptr)) then
    call c_f_pointer(q6%data_ptr, f_q6_ptr, [q6%dims(1)])
    f_q6 = f_q6_ptr(:)
  else
    f_q6_ptr => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(q7%data_ptr)) then
    call c_f_pointer(q7%data_ptr, f_q7_ptr, [q7%dims(1)])
    f_q7 = f_q7_ptr(:)
  else
    f_q7_ptr => null()
  endif
  !! general array (1D_NOT_complex) inout
  if (c_associated(q8%data_ptr)) then
    call c_f_pointer(q8%data_ptr, f_q8_ptr, [q8%dims(1)])
    f_q8 = f_q8_ptr(:)
  else
    f_q8_ptr => null()
  endif
  !! general array (1D_NOT_complex) in
  if (c_associated(q9%data_ptr)) then
    call c_f_pointer(q9%data_ptr, f_q9_ptr, [q9%dims(1)])
    f_q9 = f_q9_ptr(:)
  else
    f_q9_ptr => null()
  endif
  !! general array (1D_NOT_complex) out
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    ! output-only
  else
    f_q_out_ptr => null()
  endif
  f_q_out = quat_mul(f_q1, f_q2, f_q3, f_q4, f_q5, f_q6, f_q7, f_q8, f_q9)

  ! out: f_q_out 1D_NOT_complex
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    f_q_out_ptr = f_q_out(:)
  endif
end subroutine
subroutine fortran_quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9, q_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: q1
  real(rp) :: f_q1(0:3)
  real(c_double), pointer :: f_q1_ptr(:)
  type(array_descriptor_t), intent(in) :: q2
  real(rp) :: f_q2(0:3)
  real(c_double), pointer :: f_q2_ptr(:)
  type(array_descriptor_t), intent(in) :: q3
  real(rp) :: f_q3(0:3)
  real(c_double), pointer :: f_q3_ptr(:)
  type(array_descriptor_t), intent(in) :: q9
  real(rp) :: f_q9(0:3)
  real(c_double), pointer :: f_q9_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: q_out
  real(rp) :: f_q_out(0:3)
  real(c_double), pointer :: f_q_out_ptr(:)
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: q4
  real(rp) :: f_q4(0:3)
  real(c_double), pointer :: f_q4_ptr(:)
  type(array_descriptor_t), intent(in) :: q5
  real(rp) :: f_q5(0:3)
  real(c_double), pointer :: f_q5_ptr(:)
  type(array_descriptor_t), intent(in) :: q6
  real(rp) :: f_q6(0:3)
  real(c_double), pointer :: f_q6_ptr(:)
  type(array_descriptor_t), intent(in) :: q7
  real(rp) :: f_q7(0:3)
  real(c_double), pointer :: f_q7_ptr(:)
  type(array_descriptor_t), intent(in) :: q8
  real(rp) :: f_q8(0:3)
  real(c_double), pointer :: f_q8_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(q1%data_ptr)) then
    call c_f_pointer(q1%data_ptr, f_q1_ptr, [q1%dims(1)])
    f_q1 = f_q1_ptr(:)
  else
    f_q1_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(q2%data_ptr)) then
    call c_f_pointer(q2%data_ptr, f_q2_ptr, [q2%dims(1)])
    f_q2 = f_q2_ptr(:)
  else
    f_q2_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(q3%data_ptr)) then
    call c_f_pointer(q3%data_ptr, f_q3_ptr, [q3%dims(1)])
    f_q3 = f_q3_ptr(:)
  else
    f_q3_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(q4%data_ptr)) then
    call c_f_pointer(q4%data_ptr, f_q4_ptr, [q4%dims(1)])
    f_q4 = f_q4_ptr(:)
  else
    f_q4_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(q5%data_ptr)) then
    call c_f_pointer(q5%data_ptr, f_q5_ptr, [q5%dims(1)])
    f_q5 = f_q5_ptr(:)
  else
    f_q5_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(q6%data_ptr)) then
    call c_f_pointer(q6%data_ptr, f_q6_ptr, [q6%dims(1)])
    f_q6 = f_q6_ptr(:)
  else
    f_q6_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(q7%data_ptr)) then
    call c_f_pointer(q7%data_ptr, f_q7_ptr, [q7%dims(1)])
    f_q7 = f_q7_ptr(:)
  else
    f_q7_ptr => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(q8%data_ptr)) then
    call c_f_pointer(q8%data_ptr, f_q8_ptr, [q8%dims(1)])
    f_q8 = f_q8_ptr(:)
  else
    f_q8_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(q9%data_ptr)) then
    call c_f_pointer(q9%data_ptr, f_q9_ptr, [q9%dims(1)])
    f_q9 = f_q9_ptr(:)
  else
    f_q9_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    ! output-only
  else
    f_q_out_ptr => null()
  endif
  f_q_out = quat_mul(f_q1, f_q2, f_q3, f_q4, f_q5, f_q6, f_q7, f_q8, f_q9)

  ! out: f_q_out 1D_NOT_real
  if (c_associated(q_out%data_ptr)) then
    call c_f_pointer(q_out%data_ptr, f_q_out_ptr, [q_out%dims(1)])
    f_q_out_ptr = f_q_out(:)
  endif
end subroutine
subroutine fortran_quat_rotate_complex (quat, vec_in, vec_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: quat
  complex(rp) :: f_quat(0:3)
  complex(c_double_complex), pointer :: f_quat_ptr(:)
  type(array_descriptor_t), intent(in) :: vec_in
  complex(rp) :: f_vec_in(3)
  complex(c_double_complex), pointer :: f_vec_in_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: vec_out
  complex(rp) :: f_vec_out(3)
  complex(c_double_complex), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex) in
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat = f_quat_ptr(:)
  else
    f_quat_ptr => null()
  endif
  !! general array (1D_NOT_complex) in
  if (c_associated(vec_in%data_ptr)) then
    call c_f_pointer(vec_in%data_ptr, f_vec_in_ptr, [vec_in%dims(1)])
    f_vec_in = f_vec_in_ptr(:)
  else
    f_vec_in_ptr => null()
  endif
  !! general array (1D_NOT_complex) out
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    ! output-only
  else
    f_vec_out_ptr => null()
  endif
  f_vec_out = quat_rotate(f_quat, f_vec_in)

  ! out: f_vec_out 1D_NOT_complex
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    f_vec_out_ptr = f_vec_out(:)
  endif
end subroutine
subroutine fortran_quat_rotate_real (quat, vec_in, vec_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  type(array_descriptor_t), intent(in) :: vec_in
  real(rp) :: f_vec_in(3)
  real(c_double), pointer :: f_vec_in_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: vec_out
  real(rp) :: f_vec_out(3)
  real(c_double), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat = f_quat_ptr(:)
  else
    f_quat_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(vec_in%data_ptr)) then
    call c_f_pointer(vec_in%data_ptr, f_vec_in_ptr, [vec_in%dims(1)])
    f_vec_in = f_vec_in_ptr(:)
  else
    f_vec_in_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    ! output-only
  else
    f_vec_out_ptr => null()
  endif
  f_vec_out = quat_rotate(f_quat, f_vec_in)

  ! out: f_vec_out 1D_NOT_real
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    f_vec_out_ptr = f_vec_out(:)
  endif
end subroutine
subroutine fortran_quat_to_axis_angle (quat, axis, angle) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  type(c_ptr), intent(in), value :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  real(c_double), pointer :: f_angle_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat = f_quat_ptr(:)
  else
    f_quat_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    ! output-only
  else
    f_axis_ptr => null()
  endif
  call quat_to_axis_angle(f_quat, f_axis, f_angle)

  ! out: f_axis 1D_NOT_real
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    f_axis_ptr = f_axis(:)
  endif
  ! out: f_angle 0D_NOT_real
  call c_f_pointer(angle, f_angle_ptr)
  f_angle_ptr = f_angle
end subroutine
subroutine fortran_quat_to_omega (quat, omega) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: omega
  real(rp) :: f_omega(3)
  real(c_double), pointer :: f_omega_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat = f_quat_ptr(:)
  else
    f_quat_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(omega%data_ptr)) then
    call c_f_pointer(omega%data_ptr, f_omega_ptr, [omega%dims(1)])
    ! output-only
  else
    f_omega_ptr => null()
  endif
  f_omega = quat_to_omega(f_quat)

  ! out: f_omega 1D_NOT_real
  if (c_associated(omega%data_ptr)) then
    call c_f_pointer(omega%data_ptr, f_omega_ptr, [omega%dims(1)])
    f_omega_ptr = f_omega(:)
  endif
end subroutine
subroutine fortran_quat_to_w_mat (quat, w_mat) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat = f_quat_ptr(:)
  else
    f_quat_ptr => null()
  endif
  !! general array (2D_NOT_real) out
  if (c_associated(w_mat%data_ptr)) then
    call c_f_pointer(w_mat%data_ptr, f_w_mat_ptr, [product(w_mat%dims(1:w_mat%rank))])
    ! output-only
  else
    f_w_mat_ptr => null()
  endif
  f_w_mat = quat_to_w_mat(f_quat)

  ! out: f_w_mat 2D_NOT_real
  if (c_associated(w_mat%data_ptr)) f_w_mat_ptr = mat2vec(f_w_mat, product(w_mat%dims(1:w_mat%rank)))
end subroutine
subroutine fortran_query_string (query_str, upcase, return_str, ix, ios) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: query_str
  character(len=4096), target :: f_query_str
  character(kind=c_char), pointer :: f_query_str_ptr(:)
  logical(c_bool) :: upcase  ! 0D_NOT_logical
  logical :: f_upcase
  type(c_ptr), intent(in), value :: return_str
  character(len=4096), target :: f_return_str
  character(kind=c_char), pointer :: f_return_str_ptr(:)
  integer(c_int) :: ix  ! 0D_NOT_integer
  integer :: f_ix
  integer(c_int) :: ios  ! 0D_NOT_integer
  integer :: f_ios
  ! ** End of parameters **
  ! in: f_query_str 0D_NOT_character
  if (.not. c_associated(query_str)) return
  call c_f_pointer(query_str, f_query_str_ptr, [huge(0)])
  call to_f_str(f_query_str_ptr, f_query_str)
  ! in: f_upcase 0D_NOT_logical
  f_upcase = upcase
  ! in: f_return_str 0D_NOT_character
  if (.not. c_associated(return_str)) return
  call c_f_pointer(return_str, f_return_str_ptr, [huge(0)])
  call to_f_str(f_return_str_ptr, f_return_str)
  ! in: f_ix 0D_NOT_integer
  f_ix = ix
  ! in: f_ios 0D_NOT_integer
  f_ios = ios
  call query_string(f_query_str, f_upcase, f_return_str, f_ix, f_ios)

end subroutine
subroutine fortran_quote (str, q_str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_str
  character(len=4096), target :: f_q_str
  character(kind=c_char), pointer :: f_q_str_ptr(:)
  ! ** End of parameters **
  ! in: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  f_q_str = quote(f_str)

  ! out: f_q_str 0D_ALLOC_character
  call c_f_pointer(q_str, f_q_str_ptr, [len_trim(f_q_str) + 1])
  call to_c_str(f_q_str, f_q_str_ptr)
end subroutine
subroutine fortran_quoten (str, delim, q_str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: delim
  character(len=4096), target :: f_delim
  character(kind=c_char), pointer :: f_delim_ptr(:)
  character(len=4096), pointer :: f_delim_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_str
  character(len=4096), target :: f_q_str
  character(kind=c_char), pointer :: f_q_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  type(character_container_alloc), pointer :: f_str
  character(200), allocatable :: f_str_local(:)
  ! ** End of parameters **
  !! container character array (1D_NOT_character)
  if (c_associated(str))   call c_f_pointer(str, f_str)
  if (c_associated(str) .and. allocated(f_str%data)) then
    allocate(f_str_local, mold=f_str%data)
    f_str_local = ''
  endif
  ! in: f_delim 0D_NOT_character
  if (c_associated(delim)) then
    call c_f_pointer(delim, f_delim_ptr, [huge(0)])
    call to_f_str(f_delim_ptr, f_delim)
    f_delim_call_ptr => f_delim
  else
    f_delim_call_ptr => null()
  endif
  f_q_str = quoten(f_str_local, f_delim_call_ptr)

  !! copy allocatable character result into container
  if (c_associated(str) .and. allocated(f_str_local)) then
    if (allocated(f_str%data)) deallocate(f_str%data)
    allocate(f_str%data, source=f_str_local)
  endif
  ! out: f_q_str 0D_ALLOC_character
  call c_f_pointer(q_str, f_q_str_ptr, [len_trim(f_q_str) + 1])
  call to_c_str(f_q_str, f_q_str_ptr)
end subroutine
subroutine fortran_ran_default_state (set_state, get_state) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: set_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_set_state
  ! ** Out parameters **
  type(c_ptr), value :: get_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_get_state
  ! ** End of parameters **
  ! in: f_set_state 0D_NOT_type
  if (c_associated(set_state))   call c_f_pointer(set_state, f_set_state)
  ! out: f_get_state 0D_NOT_type
  if (c_associated(get_state))   call c_f_pointer(get_state, f_get_state)
  call ran_default_state(f_set_state, f_get_state)

  ! out: f_get_state 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_ran_engine (set, get, ran_state) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  character(len=4096), pointer :: f_set_call_ptr
  type(c_ptr), intent(in), value :: get
  character(len=4096), target :: f_get
  character(kind=c_char), pointer :: f_get_ptr(:)
  character(len=4096), pointer :: f_get_call_ptr
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  ! ** End of parameters **
  ! in: f_set 0D_NOT_character
  if (c_associated(set)) then
    call c_f_pointer(set, f_set_ptr, [huge(0)])
    call to_f_str(f_set_ptr, f_set)
    f_set_call_ptr => f_set
  else
    f_set_call_ptr => null()
  endif
  ! in: f_get 0D_NOT_character
  if (c_associated(get)) then
    call c_f_pointer(get, f_get_ptr, [huge(0)])
    call to_f_str(f_get_ptr, f_get)
    f_get_call_ptr => f_get
  else
    f_get_call_ptr => null()
  endif
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  call ran_engine(f_set_call_ptr, f_get_call_ptr, f_ran_state)

end subroutine
subroutine fortran_ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut, ran_state) &
    bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  character(len=4096), pointer :: f_set_call_ptr
  type(c_ptr), intent(in), value :: set_sigma_cut  ! 0D_NOT_real
  real(c_double) :: f_set_sigma_cut
  real(c_double), pointer :: f_set_sigma_cut_ptr
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: get
  character(len=4096), target :: f_get
  character(kind=c_char), pointer :: f_get_ptr(:)
  character(len=4096), pointer :: f_get_call_ptr
  type(c_ptr), intent(in), value :: get_sigma_cut  ! 0D_NOT_real
  real(rp) :: f_get_sigma_cut
  real(c_double), pointer :: f_get_sigma_cut_ptr
  ! ** End of parameters **
  ! in: f_set 0D_NOT_character
  if (c_associated(set)) then
    call c_f_pointer(set, f_set_ptr, [huge(0)])
    call to_f_str(f_set_ptr, f_set)
    f_set_call_ptr => f_set
  else
    f_set_call_ptr => null()
  endif
  ! in: f_set_sigma_cut 0D_NOT_real
  if (c_associated(set_sigma_cut)) then
    call c_f_pointer(set_sigma_cut, f_set_sigma_cut_ptr)
  else
    f_set_sigma_cut_ptr => null()
  endif
  ! out: f_get 0D_NOT_character
  if (c_associated(get)) then
    call c_f_pointer(get, f_get_ptr, [huge(0)])
    f_get_call_ptr => f_get
  else
    f_get_call_ptr => null()
  endif
  ! out: f_get_sigma_cut 0D_NOT_real
  if (c_associated(get_sigma_cut)) then
    call c_f_pointer(get_sigma_cut, f_get_sigma_cut_ptr)
  else
    f_get_sigma_cut_ptr => null()
  endif
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  call ran_gauss_converter(f_set_call_ptr, f_set_sigma_cut_ptr, f_get_call_ptr, &
      f_get_sigma_cut, f_ran_state)

  ! out: f_get 0D_NOT_character
  if (c_associated(get)) then
    call c_f_pointer(get, f_get_ptr, [len_trim(f_get) + 1])
    call to_c_str(f_get, f_get_ptr)
  endif
  ! out: f_get_sigma_cut 0D_NOT_real
  call c_f_pointer(get_sigma_cut, f_get_sigma_cut_ptr)
  f_get_sigma_cut_ptr = f_get_sigma_cut
end subroutine
subroutine fortran_ran_gauss_scalar (harvest, ran_state, sigma_cut, index_quasi) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  type(c_ptr), intent(in), value :: sigma_cut  ! 0D_NOT_real
  real(c_double) :: f_sigma_cut
  real(c_double), pointer :: f_sigma_cut_ptr
  type(c_ptr), intent(in), value :: index_quasi  ! 0D_NOT_integer
  integer(c_int) :: f_index_quasi
  integer(c_int), pointer :: f_index_quasi_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: harvest  ! 0D_NOT_real
  real(rp) :: f_harvest
  real(c_double), pointer :: f_harvest_ptr
  ! ** End of parameters **
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  ! in: f_sigma_cut 0D_NOT_real
  if (c_associated(sigma_cut)) then
    call c_f_pointer(sigma_cut, f_sigma_cut_ptr)
  else
    f_sigma_cut_ptr => null()
  endif
  ! in: f_index_quasi 0D_NOT_integer
  if (c_associated(index_quasi)) then
    call c_f_pointer(index_quasi, f_index_quasi_ptr)
  else
    f_index_quasi_ptr => null()
  endif
  call ran_gauss_scalar(f_harvest, f_ran_state, f_sigma_cut_ptr, f_index_quasi_ptr)

  ! out: f_harvest 0D_NOT_real
  call c_f_pointer(harvest, f_harvest_ptr)
  f_harvest_ptr = f_harvest
end subroutine
subroutine fortran_ran_gauss_vector (harvest, ran_state, sigma_cut) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  type(c_ptr), intent(in), value :: sigma_cut  ! 0D_NOT_real
  real(c_double) :: f_sigma_cut
  real(c_double), pointer :: f_sigma_cut_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: harvest
  real(rp), pointer :: f_harvest(:)
  real(c_double), pointer :: f_harvest_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(harvest%data_ptr)) then
    call c_f_pointer(harvest%data_ptr, f_harvest_ptr, [harvest%dims(1)])
    f_harvest => f_harvest_ptr
  else
    f_harvest => null()
  endif
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  ! in: f_sigma_cut 0D_NOT_real
  if (c_associated(sigma_cut)) then
    call c_f_pointer(sigma_cut, f_sigma_cut_ptr)
  else
    f_sigma_cut_ptr => null()
  endif
  call ran_gauss_vector(f_harvest, f_ran_state, f_sigma_cut_ptr)

end subroutine
subroutine fortran_ran_seed_get (seed) bind(c)

  use array_desc_mod
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: seed  ! 0D_NOT_integer
  integer :: f_seed
  integer(c_int), pointer :: f_seed_ptr
  ! ** End of parameters **
  call ran_seed_get(f_seed)

  ! out: f_seed 0D_NOT_integer
  call c_f_pointer(seed, f_seed_ptr)
  f_seed_ptr = f_seed
end subroutine
subroutine fortran_ran_seed_put (seed, mpi_offset) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: seed  ! 0D_NOT_integer
  integer :: f_seed
  type(c_ptr), intent(in), value :: mpi_offset  ! 0D_NOT_integer
  integer(c_int) :: f_mpi_offset
  integer(c_int), pointer :: f_mpi_offset_ptr
  ! ** End of parameters **
  ! in: f_seed 0D_NOT_integer
  f_seed = seed
  ! in: f_mpi_offset 0D_NOT_integer
  if (c_associated(mpi_offset)) then
    call c_f_pointer(mpi_offset, f_mpi_offset_ptr)
  else
    f_mpi_offset_ptr => null()
  endif
  call ran_seed_put(f_seed, f_mpi_offset_ptr)

end subroutine
subroutine fortran_ran_uniform_scalar (harvest, ran_state, index_quasi) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  type(c_ptr), intent(in), value :: index_quasi  ! 0D_NOT_integer
  integer(c_int) :: f_index_quasi
  integer(c_int), pointer :: f_index_quasi_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: harvest  ! 0D_NOT_real
  real(rp) :: f_harvest
  real(c_double), pointer :: f_harvest_ptr
  ! ** End of parameters **
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  ! in: f_index_quasi 0D_NOT_integer
  if (c_associated(index_quasi)) then
    call c_f_pointer(index_quasi, f_index_quasi_ptr)
  else
    f_index_quasi_ptr => null()
  endif
  call ran_uniform(f_harvest, f_ran_state, f_index_quasi_ptr)

  ! out: f_harvest 0D_NOT_real
  call c_f_pointer(harvest, f_harvest_ptr)
  f_harvest_ptr = f_harvest
end subroutine
subroutine fortran_ran_uniform_vector (harvest, ran_state) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: harvest
  real(rp), pointer :: f_harvest(:)
  real(c_double), pointer :: f_harvest_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(harvest%data_ptr)) then
    call c_f_pointer(harvest%data_ptr, f_harvest_ptr, [harvest%dims(1)])
    f_harvest => f_harvest_ptr
  else
    f_harvest => null()
  endif
  ! in: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  call ran_uniform(f_harvest, f_ran_state)

end subroutine
subroutine fortran_rcelbd (mc, elb, eld) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: mc  ! 0D_NOT_real
  real(sp) :: f_mc
  real(c_double) :: elb  ! 0D_NOT_real
  real(sp) :: f_elb
  real(c_double) :: eld  ! 0D_NOT_real
  real(sp) :: f_eld
  ! ** End of parameters **
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  ! in: f_elb 0D_NOT_real
  f_elb = elb
  ! in: f_eld 0D_NOT_real
  f_eld = eld
  call rcelbd(f_mc, f_elb, f_eld)

end subroutine
subroutine fortran_read_a_line (prompt, line_out, trim_prompt, prompt_color, prompt_bold, &
    history_file) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: prompt
  character(len=4096), target :: f_prompt
  character(kind=c_char), pointer :: f_prompt_ptr(:)
  type(c_ptr), intent(in), value :: trim_prompt  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_trim_prompt
  logical, target :: f_trim_prompt_native
  logical, pointer :: f_trim_prompt_native_ptr
  logical(c_bool), pointer :: f_trim_prompt_ptr
  type(c_ptr), intent(in), value :: prompt_color
  character(len=4096), target :: f_prompt_color
  character(kind=c_char), pointer :: f_prompt_color_ptr(:)
  character(len=4096), pointer :: f_prompt_color_call_ptr
  type(c_ptr), intent(in), value :: prompt_bold  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_prompt_bold
  logical, target :: f_prompt_bold_native
  logical, pointer :: f_prompt_bold_native_ptr
  logical(c_bool), pointer :: f_prompt_bold_ptr
  type(c_ptr), intent(in), value :: history_file
  character(len=4096), target :: f_history_file
  character(kind=c_char), pointer :: f_history_file_ptr(:)
  character(len=4096), pointer :: f_history_file_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: line_out
  character(len=4096), target :: f_line_out
  character(kind=c_char), pointer :: f_line_out_ptr(:)
  ! ** End of parameters **
  ! in: f_prompt 0D_NOT_character
  if (.not. c_associated(prompt)) return
  call c_f_pointer(prompt, f_prompt_ptr, [huge(0)])
  call to_f_str(f_prompt_ptr, f_prompt)
  ! in: f_trim_prompt 0D_NOT_logical
  if (c_associated(trim_prompt)) then
    call c_f_pointer(trim_prompt, f_trim_prompt_ptr)
    f_trim_prompt_native = f_trim_prompt_ptr
    f_trim_prompt_native_ptr => f_trim_prompt_native
  else
    f_trim_prompt_native_ptr => null()
  endif
  ! in: f_prompt_color 0D_NOT_character
  if (c_associated(prompt_color)) then
    call c_f_pointer(prompt_color, f_prompt_color_ptr, [huge(0)])
    call to_f_str(f_prompt_color_ptr, f_prompt_color)
    f_prompt_color_call_ptr => f_prompt_color
  else
    f_prompt_color_call_ptr => null()
  endif
  ! in: f_prompt_bold 0D_NOT_logical
  if (c_associated(prompt_bold)) then
    call c_f_pointer(prompt_bold, f_prompt_bold_ptr)
    f_prompt_bold_native = f_prompt_bold_ptr
    f_prompt_bold_native_ptr => f_prompt_bold_native
  else
    f_prompt_bold_native_ptr => null()
  endif
  ! in: f_history_file 0D_NOT_character
  if (c_associated(history_file)) then
    call c_f_pointer(history_file, f_history_file_ptr, [huge(0)])
    call to_f_str(f_history_file_ptr, f_history_file)
    f_history_file_call_ptr => f_history_file
  else
    f_history_file_call_ptr => null()
  endif
  call read_a_line(f_prompt, f_line_out, f_trim_prompt_native_ptr, f_prompt_color_call_ptr, &
      f_prompt_bold_native_ptr, f_history_file_call_ptr)

  ! out: f_line_out 0D_NOT_character
  call c_f_pointer(line_out, f_line_out_ptr, [len_trim(f_line_out) + 1])
  call to_c_str(f_line_out, f_line_out_ptr)
end subroutine
subroutine fortran_readline_read_history (history_file, status) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: history_file
  character(len=4096), target :: f_history_file
  character(kind=c_char), pointer :: f_history_file_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: status  ! 0D_NOT_integer
  integer :: f_status
  integer(c_int), pointer :: f_status_ptr
  ! ** End of parameters **
  ! in: f_history_file 0D_NOT_character
  if (.not. c_associated(history_file)) return
  call c_f_pointer(history_file, f_history_file_ptr, [huge(0)])
  call to_f_str(f_history_file_ptr, f_history_file)
  call readline_read_history(f_history_file, f_status)

  ! out: f_status 0D_NOT_integer
  call c_f_pointer(status, f_status_ptr)
  f_status_ptr = f_status
end subroutine
subroutine fortran_readline_write_history (history_file, status) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: history_file
  character(len=4096), target :: f_history_file
  character(kind=c_char), pointer :: f_history_file_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: status  ! 0D_NOT_integer
  integer :: f_status
  integer(c_int), pointer :: f_status_ptr
  ! ** End of parameters **
  ! in: f_history_file 0D_NOT_character
  if (.not. c_associated(history_file)) return
  call c_f_pointer(history_file, f_history_file_ptr, [huge(0)])
  call to_f_str(f_history_file_ptr, f_history_file)
  call readline_write_history(f_history_file, f_status)

  ! out: f_status 0D_NOT_integer
  call c_f_pointer(status, f_status_ptr)
  f_status_ptr = f_status
end subroutine
subroutine fortran_real_num_fortran_format (number, width, n_blanks, fmt_str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: number  ! 0D_NOT_real
  real(rp) :: f_number
  integer(c_int) :: width  ! 0D_NOT_integer
  integer :: f_width
  type(c_ptr), intent(in), value :: n_blanks  ! 0D_NOT_integer
  integer(c_int) :: f_n_blanks
  integer(c_int), pointer :: f_n_blanks_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: fmt_str
  character(len=4096), target :: f_fmt_str
  character(kind=c_char), pointer :: f_fmt_str_ptr(:)
  ! ** End of parameters **
  ! in: f_number 0D_NOT_real
  f_number = number
  ! in: f_width 0D_NOT_integer
  f_width = width
  ! in: f_n_blanks 0D_NOT_integer
  if (c_associated(n_blanks)) then
    call c_f_pointer(n_blanks, f_n_blanks_ptr)
  else
    f_n_blanks_ptr => null()
  endif
  f_fmt_str = real_num_fortran_format(f_number, f_width, f_n_blanks_ptr)

  ! out: f_fmt_str 0D_NOT_character
  call c_f_pointer(fmt_str, f_fmt_str_ptr, [len_trim(f_fmt_str) + 1])
  call to_c_str(f_fmt_str, f_fmt_str_ptr)
end subroutine
subroutine fortran_real_path (path_in, path_out, is_ok) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: path_in
  character(len=4096), target :: f_path_in
  character(kind=c_char), pointer :: f_path_in_ptr(:)
  type(c_ptr), intent(in), value :: path_out
  character(len=4096), target :: f_path_out
  character(kind=c_char), pointer :: f_path_out_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_ok  ! 0D_NOT_logical
  logical :: f_is_ok
  logical(c_bool), pointer :: f_is_ok_ptr
  ! ** End of parameters **
  ! in: f_path_in 0D_NOT_character
  if (.not. c_associated(path_in)) return
  call c_f_pointer(path_in, f_path_in_ptr, [huge(0)])
  call to_f_str(f_path_in_ptr, f_path_in)
  ! in: f_path_out 0D_NOT_character
  if (.not. c_associated(path_out)) return
  call c_f_pointer(path_out, f_path_out_ptr, [huge(0)])
  call to_f_str(f_path_out_ptr, f_path_out)
  f_is_ok = real_path(f_path_in, f_path_out)

  ! out: f_is_ok 0D_NOT_logical
  call c_f_pointer(is_ok, f_is_ok_ptr)
  f_is_ok_ptr = f_is_ok
end subroutine
subroutine fortran_real_str (r_num, n_signif, n_decimal, str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: r_num  ! 0D_NOT_real
  real(rp) :: f_r_num
  type(c_ptr), intent(in), value :: n_signif  ! 0D_NOT_integer
  integer(c_int) :: f_n_signif
  integer(c_int), pointer :: f_n_signif_ptr
  type(c_ptr), intent(in), value :: n_decimal  ! 0D_NOT_integer
  integer(c_int) :: f_n_decimal
  integer(c_int), pointer :: f_n_decimal_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! in: f_r_num 0D_NOT_real
  f_r_num = r_num
  ! in: f_n_signif 0D_NOT_integer
  if (c_associated(n_signif)) then
    call c_f_pointer(n_signif, f_n_signif_ptr)
  else
    f_n_signif_ptr => null()
  endif
  ! in: f_n_decimal 0D_NOT_integer
  if (c_associated(n_decimal)) then
    call c_f_pointer(n_decimal, f_n_decimal_ptr)
  else
    f_n_decimal_ptr => null()
  endif
  f_str = real_str(f_r_num, f_n_signif_ptr, f_n_decimal_ptr)

  ! out: f_str 0D_ALLOC_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1])
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_real_to_string (real_num, width, n_signif, n_decimal, str) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: real_num  ! 0D_NOT_real
  real(rp) :: f_real_num
  integer(c_int) :: width  ! 0D_NOT_integer
  integer :: f_width
  type(c_ptr), intent(in), value :: n_signif  ! 0D_NOT_integer
  integer(c_int) :: f_n_signif
  integer(c_int), pointer :: f_n_signif_ptr
  type(c_ptr), intent(in), value :: n_decimal  ! 0D_NOT_integer
  integer(c_int) :: f_n_decimal
  integer(c_int), pointer :: f_n_decimal_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! in: f_real_num 0D_NOT_real
  f_real_num = real_num
  ! in: f_width 0D_NOT_integer
  f_width = width
  ! in: f_n_signif 0D_NOT_integer
  if (c_associated(n_signif)) then
    call c_f_pointer(n_signif, f_n_signif_ptr)
  else
    f_n_signif_ptr => null()
  endif
  ! in: f_n_decimal 0D_NOT_integer
  if (c_associated(n_decimal)) then
    call c_f_pointer(n_decimal, f_n_decimal_ptr)
  else
    f_n_decimal_ptr => null()
  endif
  f_str = real_to_string(f_real_num, f_width, f_n_signif_ptr, f_n_decimal_ptr)

  ! out: f_str 0D_NOT_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1])
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_reallocate_spline (spline, n, n_min, exact) bind(c)

  use array_desc_mod
  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: n  ! 0D_NOT_integer
  integer :: f_n
  type(c_ptr), intent(in), value :: n_min  ! 0D_NOT_integer
  integer(c_int) :: f_n_min
  integer(c_int), pointer :: f_n_min_ptr
  type(c_ptr), intent(in), value :: exact  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact
  logical, target :: f_exact_native
  logical, pointer :: f_exact_native_ptr
  logical(c_bool), pointer :: f_exact_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: spline
  type(spline_struct_container_alloc), pointer :: f_spline
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(spline))   call c_f_pointer(spline, f_spline)
  ! in: f_n 0D_NOT_integer
  f_n = n
  ! in: f_n_min 0D_NOT_integer
  if (c_associated(n_min)) then
    call c_f_pointer(n_min, f_n_min_ptr)
  else
    f_n_min_ptr => null()
  endif
  ! in: f_exact 0D_NOT_logical
  if (c_associated(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_native = f_exact_ptr
    f_exact_native_ptr => f_exact_native
  else
    f_exact_native_ptr => null()
  endif
  call reallocate_spline(f_spline%data, f_n, f_n_min_ptr, f_exact_native_ptr)

end subroutine
subroutine fortran_relbd (phi, phic, mc, b, d) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: phi  ! 0D_NOT_real
  real(sp) :: f_phi
  real(c_double) :: phic  ! 0D_NOT_real
  real(sp) :: f_phic
  real(c_double) :: mc  ! 0D_NOT_real
  real(sp) :: f_mc
  real(c_double) :: b  ! 0D_NOT_real
  real(sp) :: f_b
  real(c_double) :: d  ! 0D_NOT_real
  real(sp) :: f_d
  ! ** End of parameters **
  ! in: f_phi 0D_NOT_real
  f_phi = phi
  ! in: f_phic 0D_NOT_real
  f_phic = phic
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  ! in: f_b 0D_NOT_real
  f_b = b
  ! in: f_d 0D_NOT_real
  f_d = d
  call relbd(f_phi, f_phic, f_mc, f_b, f_d)

end subroutine
subroutine fortran_relcbd (c0, mc, b, dx) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: c0  ! 0D_NOT_real
  real(sp) :: f_c0
  real(c_double) :: mc  ! 0D_NOT_real
  real(sp) :: f_mc
  real(c_double) :: b  ! 0D_NOT_real
  real(sp) :: f_b
  real(c_double) :: dx  ! 0D_NOT_real
  real(sp) :: f_dx
  ! ** End of parameters **
  ! in: f_c0 0D_NOT_real
  f_c0 = c0
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  ! in: f_b 0D_NOT_real
  f_b = b
  ! in: f_dx 0D_NOT_real
  f_dx = dx
  call relcbd(f_c0, f_mc, f_b, f_dx)

end subroutine
subroutine fortran_relsbd (s0, mc, b, d) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: s0  ! 0D_NOT_real
  real(sp) :: f_s0
  real(c_double) :: mc  ! 0D_NOT_real
  real(sp) :: f_mc
  real(c_double) :: b  ! 0D_NOT_real
  real(sp) :: f_b
  real(c_double) :: d  ! 0D_NOT_real
  real(sp) :: f_d
  ! ** End of parameters **
  ! in: f_s0 0D_NOT_real
  f_s0 = s0
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  ! in: f_b 0D_NOT_real
  f_b = b
  ! in: f_d 0D_NOT_real
  f_d = d
  call relsbd(f_s0, f_mc, f_b, f_d)

end subroutine
subroutine fortran_rgelbd (phi, mc, elb, eld) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: phi  ! 0D_NOT_real
  real(sp) :: f_phi
  real(c_double) :: mc  ! 0D_NOT_real
  real(sp) :: f_mc
  real(c_double) :: elb  ! 0D_NOT_real
  real(sp) :: f_elb
  real(c_double) :: eld  ! 0D_NOT_real
  real(sp) :: f_eld
  ! ** End of parameters **
  ! in: f_phi 0D_NOT_real
  f_phi = phi
  ! in: f_mc 0D_NOT_real
  f_mc = mc
  ! in: f_elb 0D_NOT_real
  f_elb = elb
  ! in: f_eld 0D_NOT_real
  f_eld = eld
  call rgelbd(f_phi, f_mc, f_elb, f_eld)

end subroutine
subroutine fortran_rms_value (val_arr, good_val, ave_val, rms_val) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: val_arr
  real(rp), pointer :: f_val_arr(:)
  real(c_double), pointer :: f_val_arr_ptr(:)
  type(c_ptr), intent(in), value :: good_val
  type(logical_container_alloc), pointer :: f_good_val
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ave_val  ! 0D_NOT_real
  real(rp) :: f_ave_val
  real(c_double), pointer :: f_ave_val_ptr
  type(c_ptr), intent(in), value :: rms_val  ! 0D_NOT_real
  real(rp) :: f_rms_val
  real(c_double), pointer :: f_rms_val_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(val_arr%data_ptr)) then
    call c_f_pointer(val_arr%data_ptr, f_val_arr_ptr, [val_arr%dims(1)])
    f_val_arr => f_val_arr_ptr
  else
    f_val_arr => null()
  endif
  !! container general array (1D_ALLOC_logical)
  if (c_associated(good_val))   call c_f_pointer(good_val, f_good_val)
  ! out: f_ave_val 0D_NOT_real
  if (c_associated(ave_val)) then
    call c_f_pointer(ave_val, f_ave_val_ptr)
  else
    f_ave_val_ptr => null()
  endif
  f_rms_val = rms_value(f_val_arr, f_good_val%data, f_ave_val)

  ! out: f_ave_val 0D_NOT_real
  call c_f_pointer(ave_val, f_ave_val_ptr)
  f_ave_val_ptr = f_ave_val
  ! out: f_rms_val 0D_NOT_real
  call c_f_pointer(rms_val, f_rms_val_ptr)
  f_rms_val_ptr = f_rms_val
end subroutine
subroutine fortran_rot_2d (vec_in, angle, vec_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: vec_in
  real(rp) :: f_vec_in(2)
  real(c_double), pointer :: f_vec_in_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: vec_out
  real(rp) :: f_vec_out(2)
  real(c_double), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(vec_in%data_ptr)) then
    call c_f_pointer(vec_in%data_ptr, f_vec_in_ptr, [vec_in%dims(1)])
    f_vec_in = f_vec_in_ptr(:)
  else
    f_vec_in_ptr => null()
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  !! general array (1D_NOT_real) out
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    ! output-only
  else
    f_vec_out_ptr => null()
  endif
  f_vec_out = rot_2d(f_vec_in, f_angle)

  ! out: f_vec_out 1D_NOT_real
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    f_vec_out_ptr = f_vec_out(:)
  endif
end subroutine
subroutine fortran_rotate_vec (vec, axis, angle) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: axis  ! 0D_NOT_integer
  integer :: f_axis
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: vec
  real(rp), pointer :: f_vec(:)
  real(c_double), pointer :: f_vec_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(vec%data_ptr)) then
    call c_f_pointer(vec%data_ptr, f_vec_ptr, [vec%dims(1)])
    f_vec => f_vec_ptr
  else
    f_vec => null()
  endif
  ! in: f_axis 0D_NOT_integer
  f_axis = axis
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  call rotate_vec(f_vec, f_axis, f_angle)

end subroutine
subroutine fortran_rotate_vec_given_axis_angle (vec_in, axis, angle, vec_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: vec_in
  real(rp) :: f_vec_in(3)
  real(c_double), pointer :: f_vec_in_ptr(:)
  type(array_descriptor_t), intent(in) :: axis
  real(rp), pointer :: f_axis(:)
  real(c_double), pointer :: f_axis_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: vec_out
  real(rp) :: f_vec_out(3)
  real(c_double), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(vec_in%data_ptr)) then
    call c_f_pointer(vec_in%data_ptr, f_vec_in_ptr, [vec_in%dims(1)])
    f_vec_in = f_vec_in_ptr(:)
  else
    f_vec_in_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    f_axis => f_axis_ptr
  else
    f_axis => null()
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  !! general array (1D_NOT_real) out
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    ! output-only
  else
    f_vec_out_ptr => null()
  endif
  f_vec_out = rotate_vec_given_axis_angle(f_vec_in, f_axis, f_angle)

  ! out: f_vec_out 1D_NOT_real
  if (c_associated(vec_out%data_ptr)) then
    call c_f_pointer(vec_out%data_ptr, f_vec_out_ptr, [vec_out%dims(1)])
    f_vec_out_ptr = f_vec_out(:)
  endif
end subroutine
subroutine fortran_rp8 (int_in, re_out) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: int_in  ! 0D_NOT_integer
  integer :: f_int_in
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: re_out  ! 0D_NOT_real
  real(rp) :: f_re_out
  real(c_double), pointer :: f_re_out_ptr
  ! ** End of parameters **
  ! in: f_int_in 0D_NOT_integer
  f_int_in = int_in
  f_re_out = rp8(f_int_in)

  ! out: f_re_out 0D_NOT_real
  call c_f_pointer(re_out, f_re_out_ptr)
  f_re_out_ptr = f_re_out
end subroutine
subroutine fortran_rserbd (y, m, b, d) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: y  ! 0D_NOT_real
  real(sp) :: f_y
  real(c_double) :: m  ! 0D_NOT_real
  real(sp) :: f_m
  real(c_double) :: b  ! 0D_NOT_real
  real(sp) :: f_b
  real(c_double) :: d  ! 0D_NOT_real
  real(sp) :: f_d
  ! ** End of parameters **
  ! in: f_y 0D_NOT_real
  f_y = y
  ! in: f_m 0D_NOT_real
  f_m = m
  ! in: f_b 0D_NOT_real
  f_b = b
  ! in: f_d 0D_NOT_real
  f_d = d
  call rserbd(f_y, f_m, f_b, f_d)

end subroutine
subroutine fortran_run_timer (command, time, time0) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: command
  character(len=4096), target :: f_command
  character(kind=c_char), pointer :: f_command_ptr(:)
  type(c_ptr), intent(in), value :: time  ! 0D_NOT_real
  real(c_double) :: f_time
  real(c_double), pointer :: f_time_ptr
  type(c_ptr), intent(in), value :: time0  ! 0D_NOT_real
  real(c_double) :: f_time0
  real(c_double), pointer :: f_time0_ptr
  ! ** End of parameters **
  ! in: f_command 0D_NOT_character
  if (.not. c_associated(command)) return
  call c_f_pointer(command, f_command_ptr, [huge(0)])
  call to_f_str(f_command_ptr, f_command)
  ! in: f_time 0D_NOT_real
  if (c_associated(time)) then
    call c_f_pointer(time, f_time_ptr)
  else
    f_time_ptr => null()
  endif
  ! in: f_time0 0D_NOT_real
  if (c_associated(time0)) then
    call c_f_pointer(time0, f_time0_ptr)
  else
    f_time0_ptr => null()
  endif
  call run_timer(f_command, f_time_ptr, f_time0_ptr)

end subroutine
subroutine fortran_serbd (y, m, b, d) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: y  ! 0D_NOT_real
  real(dp) :: f_y
  real(c_double) :: m  ! 0D_NOT_real
  real(dp) :: f_m
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: b  ! 0D_NOT_real
  real(dp) :: f_b
  real(c_double), pointer :: f_b_ptr
  type(c_ptr), intent(in), value :: d  ! 0D_NOT_real
  real(dp) :: f_d
  real(c_double), pointer :: f_d_ptr
  ! ** End of parameters **
  ! in: f_y 0D_NOT_real
  f_y = y
  ! in: f_m 0D_NOT_real
  f_m = m
  call serbd(f_y, f_m, f_b, f_d)

  ! out: f_b 0D_NOT_real
  call c_f_pointer(b, f_b_ptr)
  f_b_ptr = f_b
  ! out: f_d 0D_NOT_real
  call c_f_pointer(d, f_d_ptr)
  f_d_ptr = f_d
end subroutine
subroutine fortran_set_all_ptr (a_ptr, value, delta, value_set) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: all_pointer_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: value  ! 0D_NOT_real
  real(rp) :: f_value
  type(c_ptr), intent(in), value :: delta  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_delta
  logical, target :: f_delta_native
  logical, pointer :: f_delta_native_ptr
  logical(c_bool), pointer :: f_delta_ptr
  type(c_ptr), intent(in), value :: value_set  ! 0D_NOT_real
  real(c_double) :: f_value_set
  real(c_double), pointer :: f_value_set_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: a_ptr  ! 0D_NOT_type
  type(all_pointer_struct), pointer :: f_a_ptr
  ! ** End of parameters **
  ! inout: f_a_ptr 0D_NOT_type
  if (.not. c_associated(a_ptr)) return
  call c_f_pointer(a_ptr, f_a_ptr)
  ! in: f_value 0D_NOT_real
  f_value = value
  ! in: f_delta 0D_NOT_logical
  if (c_associated(delta)) then
    call c_f_pointer(delta, f_delta_ptr)
    f_delta_native = f_delta_ptr
    f_delta_native_ptr => f_delta_native
  else
    f_delta_native_ptr => null()
  endif
  ! in: f_value_set 0D_NOT_real
  if (c_associated(value_set)) then
    call c_f_pointer(value_set, f_value_set_ptr)
  else
    f_value_set_ptr => null()
  endif
  call set_all_ptr(f_a_ptr, f_value, f_delta_native_ptr, f_value_set_ptr)

end subroutine
subroutine fortran_set_parameter_int (param_val, set_val, save_val) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: param_val  ! 0D_NOT_integer
  integer :: f_param_val
  integer(c_int) :: set_val  ! 0D_NOT_integer
  integer :: f_set_val
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: save_val  ! 0D_NOT_integer
  integer :: f_save_val
  integer(c_int), pointer :: f_save_val_ptr
  ! ** End of parameters **
  ! in: f_param_val 0D_NOT_integer
  f_param_val = param_val
  ! in: f_set_val 0D_NOT_integer
  f_set_val = set_val
  f_save_val = set_parameter(f_param_val, f_set_val)

  ! out: f_save_val 0D_NOT_integer
  call c_f_pointer(save_val, f_save_val_ptr)
  f_save_val_ptr = f_save_val
end subroutine
subroutine fortran_set_parameter_logic (param_val, set_val, save_val) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  logical(c_bool) :: param_val  ! 0D_NOT_logical
  logical :: f_param_val
  logical(c_bool) :: set_val  ! 0D_NOT_logical
  logical :: f_set_val
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: save_val  ! 0D_NOT_logical
  logical :: f_save_val
  logical(c_bool), pointer :: f_save_val_ptr
  ! ** End of parameters **
  ! in: f_param_val 0D_NOT_logical
  f_param_val = param_val
  ! in: f_set_val 0D_NOT_logical
  f_set_val = set_val
  f_save_val = set_parameter(f_param_val, f_set_val)

  ! out: f_save_val 0D_NOT_logical
  call c_f_pointer(save_val, f_save_val_ptr)
  f_save_val_ptr = f_save_val
end subroutine
subroutine fortran_set_parameter_real (param_val, set_val, save_val) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: param_val  ! 0D_NOT_real
  real(rp) :: f_param_val
  real(c_double) :: set_val  ! 0D_NOT_real
  real(rp) :: f_set_val
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: save_val  ! 0D_NOT_real
  real(rp) :: f_save_val
  real(c_double), pointer :: f_save_val_ptr
  ! ** End of parameters **
  ! in: f_param_val 0D_NOT_real
  f_param_val = param_val
  ! in: f_set_val 0D_NOT_real
  f_set_val = set_val
  f_save_val = set_parameter(f_param_val, f_set_val)

  ! out: f_save_val 0D_NOT_real
  call c_f_pointer(save_val, f_save_val_ptr)
  f_save_val_ptr = f_save_val
end subroutine
subroutine fortran_set_species_charge (species_in, charge, species_charged) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species_in  ! 0D_NOT_integer
  integer :: f_species_in
  integer(c_int) :: charge  ! 0D_NOT_integer
  integer :: f_charge
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species_charged  ! 0D_NOT_integer
  integer :: f_species_charged
  integer(c_int), pointer :: f_species_charged_ptr
  ! ** End of parameters **
  ! in: f_species_in 0D_NOT_integer
  f_species_in = species_in
  ! in: f_charge 0D_NOT_integer
  f_charge = charge
  f_species_charged = set_species_charge(f_species_in, f_charge)

  ! out: f_species_charged 0D_NOT_integer
  call c_f_pointer(species_charged, f_species_charged_ptr)
  f_species_charged_ptr = f_species_charged
end subroutine
subroutine fortran_sign_of_int (num, zero_is_zero, num_sign) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: num  ! 0D_NOT_integer
  integer :: f_num
  type(c_ptr), intent(in), value :: zero_is_zero  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_zero_is_zero
  logical, target :: f_zero_is_zero_native
  logical, pointer :: f_zero_is_zero_native_ptr
  logical(c_bool), pointer :: f_zero_is_zero_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: num_sign  ! 0D_NOT_integer
  integer :: f_num_sign
  integer(c_int), pointer :: f_num_sign_ptr
  ! ** End of parameters **
  ! in: f_num 0D_NOT_integer
  f_num = num
  ! in: f_zero_is_zero 0D_NOT_logical
  if (c_associated(zero_is_zero)) then
    call c_f_pointer(zero_is_zero, f_zero_is_zero_ptr)
    f_zero_is_zero_native = f_zero_is_zero_ptr
    f_zero_is_zero_native_ptr => f_zero_is_zero_native
  else
    f_zero_is_zero_native_ptr => null()
  endif
  f_num_sign = sign_of(f_num, f_zero_is_zero_native_ptr)

  ! out: f_num_sign 0D_NOT_integer
  call c_f_pointer(num_sign, f_num_sign_ptr)
  f_num_sign_ptr = f_num_sign
end subroutine
subroutine fortran_sign_of_real (num, zero_is_zero, num_sign) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: num  ! 0D_NOT_real
  real(rp) :: f_num
  type(c_ptr), intent(in), value :: zero_is_zero  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_zero_is_zero
  logical, target :: f_zero_is_zero_native
  logical, pointer :: f_zero_is_zero_native_ptr
  logical(c_bool), pointer :: f_zero_is_zero_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: num_sign  ! 0D_NOT_integer
  integer :: f_num_sign
  integer(c_int), pointer :: f_num_sign_ptr
  ! ** End of parameters **
  ! in: f_num 0D_NOT_real
  f_num = num
  ! in: f_zero_is_zero 0D_NOT_logical
  if (c_associated(zero_is_zero)) then
    call c_f_pointer(zero_is_zero, f_zero_is_zero_ptr)
    f_zero_is_zero_native = f_zero_is_zero_ptr
    f_zero_is_zero_native_ptr => f_zero_is_zero_native
  else
    f_zero_is_zero_native_ptr => null()
  endif
  f_num_sign = sign_of(f_num, f_zero_is_zero_native_ptr)

  ! out: f_num_sign 0D_NOT_integer
  call c_f_pointer(num_sign, f_num_sign_ptr)
  f_num_sign_ptr = f_num_sign
end subroutine
subroutine fortran_sinc (x, nd, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_nd 0D_NOT_integer
  if (c_associated(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
  else
    f_nd_ptr => null()
  endif
  f_y = sinc(f_x, f_nd_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_sincc (x, nd, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_nd 0D_NOT_integer
  if (c_associated(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
  else
    f_nd_ptr => null()
  endif
  f_y = sincc(f_x, f_nd_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_sinhx_x (x, nd, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_nd 0D_NOT_integer
  if (c_associated(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
  else
    f_nd_ptr => null()
  endif
  f_y = sinhx_x(f_x, f_nd_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_skip_header (ix_unit, error_flag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_unit  ! 0D_NOT_integer
  integer :: f_ix_unit
  logical(c_bool) :: error_flag  ! 0D_NOT_logical
  logical :: f_error_flag
  ! ** End of parameters **
  ! in: f_ix_unit 0D_NOT_integer
  f_ix_unit = ix_unit
  ! in: f_error_flag 0D_NOT_logical
  f_error_flag = error_flag
  call skip_header(f_ix_unit, f_error_flag)

end subroutine
subroutine fortran_species_id (name, default_, print_err, species) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  type(c_ptr), intent(in), value :: default_  ! 0D_NOT_integer
  integer(c_int) :: f_default
  integer(c_int), pointer :: f_default_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species  ! 0D_NOT_integer
  integer :: f_species
  integer(c_int), pointer :: f_species_ptr
  ! ** End of parameters **
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_default 0D_NOT_integer
  if (c_associated(default_)) then
    call c_f_pointer(default_, f_default_ptr)
  else
    f_default_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  f_species = species_id(f_name, f_default_ptr, f_print_err_native_ptr)

  ! out: f_species 0D_NOT_integer
  call c_f_pointer(species, f_species_ptr)
  f_species_ptr = f_species
end subroutine
subroutine fortran_species_id_from_openpmd (pmd_name, charge, species) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: pmd_name
  character(len=4096), target :: f_pmd_name
  character(kind=c_char), pointer :: f_pmd_name_ptr(:)
  integer(c_int) :: charge  ! 0D_NOT_integer
  integer :: f_charge
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species  ! 0D_NOT_integer
  integer :: f_species
  integer(c_int), pointer :: f_species_ptr
  ! ** End of parameters **
  ! in: f_pmd_name 0D_NOT_character
  if (.not. c_associated(pmd_name)) return
  call c_f_pointer(pmd_name, f_pmd_name_ptr, [huge(0)])
  call to_f_str(f_pmd_name_ptr, f_pmd_name)
  ! in: f_charge 0D_NOT_integer
  f_charge = charge
  f_species = species_id_from_openpmd(f_pmd_name, f_charge)

  ! out: f_species 0D_NOT_integer
  call c_f_pointer(species, f_species_ptr)
  f_species_ptr = f_species
end subroutine
subroutine fortran_species_name (species, name) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_name = species_name(f_species)

  ! out: f_name 0D_NOT_character
  call c_f_pointer(name, f_name_ptr, [len_trim(f_name) + 1])
  call to_c_str(f_name, f_name_ptr)
end subroutine
subroutine fortran_species_of (mass, charge, species) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: mass  ! 0D_NOT_real
  real(rp) :: f_mass
  integer(c_int) :: charge  ! 0D_NOT_integer
  integer :: f_charge
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species  ! 0D_NOT_integer
  integer :: f_species
  integer(c_int), pointer :: f_species_ptr
  ! ** End of parameters **
  ! in: f_mass 0D_NOT_real
  f_mass = mass
  ! in: f_charge 0D_NOT_integer
  f_charge = charge
  f_species = species_of(f_mass, f_charge)

  ! out: f_species 0D_NOT_integer
  call c_f_pointer(species, f_species_ptr)
  f_species_ptr = f_species
end subroutine
subroutine fortran_spin_of (species, non_subatomic_default, spin) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  type(c_ptr), intent(in), value :: non_subatomic_default  ! 0D_NOT_real
  real(c_double) :: f_non_subatomic_default
  real(c_double), pointer :: f_non_subatomic_default_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: spin  ! 0D_NOT_real
  real(rp) :: f_spin
  real(c_double), pointer :: f_spin_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  ! in: f_non_subatomic_default 0D_NOT_real
  if (c_associated(non_subatomic_default)) then
    call c_f_pointer(non_subatomic_default, f_non_subatomic_default_ptr)
  else
    f_non_subatomic_default_ptr => null()
  endif
  f_spin = spin_of(f_species, f_non_subatomic_default_ptr)

  ! out: f_spin 0D_NOT_real
  call c_f_pointer(spin, f_spin_ptr)
  f_spin_ptr = f_spin
end subroutine
subroutine fortran_spline1 (a_spline, x, n, y) bind(c)

  use array_desc_mod
  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: a_spline  ! 0D_NOT_type
  type(spline_struct), pointer :: f_a_spline
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: n  ! 0D_NOT_integer
  integer(c_int) :: f_n
  integer(c_int), pointer :: f_n_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_a_spline 0D_NOT_type
  if (.not. c_associated(a_spline)) return
  call c_f_pointer(a_spline, f_a_spline)
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_n 0D_NOT_integer
  if (c_associated(n)) then
    call c_f_pointer(n, f_n_ptr)
  else
    f_n_ptr => null()
  endif
  f_y = spline1(f_a_spline, f_x, f_n_ptr)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_spline_akima (spline, ok) bind(c)

  use array_desc_mod
  use spline_mod, only: spline_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: spline
  type(spline_struct), pointer :: f_spline(:)
  type(spline_struct), pointer :: f_spline_ptr(:)
  ! ** End of parameters **
  !! type array (1D_NOT_type)
  if (c_associated(spline%data_ptr)) then
    call c_f_pointer(spline%data_ptr, f_spline_ptr, [spline%dims(1)])
    f_spline => f_spline_ptr
  else
    f_spline => null()
  endif
  call spline_akima(f_spline, f_ok)

  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
end subroutine
subroutine fortran_spline_akima_interpolate (x_knot, y_knot, x, ok, y, dy) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: x_knot
  real(rp), pointer :: f_x_knot(:)
  real(c_double), pointer :: f_x_knot_ptr(:)
  type(array_descriptor_t), intent(in) :: y_knot
  real(rp), pointer :: f_y_knot(:)
  real(c_double), pointer :: f_y_knot_ptr(:)
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  type(c_ptr), intent(in), value :: dy  ! 0D_NOT_real
  real(rp) :: f_dy
  real(c_double), pointer :: f_dy_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(x_knot%data_ptr)) then
    call c_f_pointer(x_knot%data_ptr, f_x_knot_ptr, [x_knot%dims(1)])
    f_x_knot => f_x_knot_ptr
  else
    f_x_knot => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y_knot%data_ptr)) then
    call c_f_pointer(y_knot%data_ptr, f_y_knot_ptr, [y_knot%dims(1)])
    f_y_knot => f_y_knot_ptr
  else
    f_y_knot => null()
  endif
  ! in: f_x 0D_NOT_real
  f_x = x
  ! out: f_y 0D_NOT_real
  if (c_associated(y)) then
    call c_f_pointer(y, f_y_ptr)
  else
    f_y_ptr => null()
  endif
  ! out: f_dy 0D_NOT_real
  if (c_associated(dy)) then
    call c_f_pointer(dy, f_dy_ptr)
  else
    f_dy_ptr => null()
  endif
  call spline_akima_interpolate(f_x_knot, f_y_knot, f_x, f_ok, f_y, f_dy)

  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
  ! out: f_dy 0D_NOT_real
  call c_f_pointer(dy, f_dy_ptr)
  f_dy_ptr = f_dy
end subroutine
subroutine fortran_spline_evaluate (spline, x, ok, y, dy) bind(c)

  use array_desc_mod
  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: spline
  type(spline_struct), pointer :: f_spline(:)
  type(spline_struct), pointer :: f_spline_ptr(:)
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  type(c_ptr), intent(in), value :: dy  ! 0D_NOT_real
  real(rp) :: f_dy
  real(c_double), pointer :: f_dy_ptr
  ! ** End of parameters **
  !! type array (1D_NOT_type)
  if (c_associated(spline%data_ptr)) then
    call c_f_pointer(spline%data_ptr, f_spline_ptr, [spline%dims(1)])
    f_spline => f_spline_ptr
  else
    f_spline => null()
  endif
  ! in: f_x 0D_NOT_real
  f_x = x
  ! out: f_y 0D_NOT_real
  if (c_associated(y)) then
    call c_f_pointer(y, f_y_ptr)
  else
    f_y_ptr => null()
  endif
  ! out: f_dy 0D_NOT_real
  if (c_associated(dy)) then
    call c_f_pointer(dy, f_dy_ptr)
  else
    f_dy_ptr => null()
  endif
  call spline_evaluate(f_spline, f_x, f_ok, f_y, f_dy)

  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
  ! out: f_dy 0D_NOT_real
  call c_f_pointer(dy, f_dy_ptr)
  f_dy_ptr = f_dy
end subroutine
subroutine fortran_sqrt_alpha (alpha, x, y) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: alpha  ! 0D_NOT_real
  real(rp) :: f_alpha
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  ! in: f_alpha 0D_NOT_real
  f_alpha = alpha
  ! in: f_x 0D_NOT_real
  f_x = x
  f_y = sqrt_alpha(f_alpha, f_x)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_sqrt_one (x, nd, ds1) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ds1  ! 0D_NOT_real
  real(rp) :: f_ds1
  real(c_double), pointer :: f_ds1_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_nd 0D_NOT_integer
  if (c_associated(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
  else
    f_nd_ptr => null()
  endif
  f_ds1 = sqrt_one(f_x, f_nd_ptr)

  ! out: f_ds1 0D_NOT_real
  call c_f_pointer(ds1, f_ds1_ptr)
  f_ds1_ptr = f_ds1
end subroutine
subroutine fortran_str_count (str, match, num) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: match
  character(len=4096), target :: f_match
  character(kind=c_char), pointer :: f_match_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: num  ! 0D_NOT_integer
  integer :: f_num
  integer(c_int), pointer :: f_num_ptr
  ! ** End of parameters **
  ! in: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! in: f_match 0D_NOT_character
  if (.not. c_associated(match)) return
  call c_f_pointer(match, f_match_ptr, [huge(0)])
  call to_f_str(f_match_ptr, f_match)
  f_num = str_count(f_str, f_match)

  ! out: f_num 0D_NOT_integer
  call c_f_pointer(num, f_num_ptr)
  f_num_ptr = f_num
end subroutine
subroutine fortran_str_downcase (dst, src) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: src
  character(len=4096), target :: f_src
  character(kind=c_char), pointer :: f_src_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: dst
  character(len=4096), target :: f_dst
  character(kind=c_char), pointer :: f_dst_ptr(:)
  ! ** End of parameters **
  ! in: f_src 0D_NOT_character
  if (.not. c_associated(src)) return
  call c_f_pointer(src, f_src_ptr, [huge(0)])
  call to_f_str(f_src_ptr, f_src)
  call str_downcase(f_dst, f_src)

  ! out: f_dst 0D_NOT_character
  call c_f_pointer(dst, f_dst_ptr, [len_trim(f_dst) + 1])
  call to_c_str(f_dst, f_dst_ptr)
end subroutine
subroutine fortran_str_first_in_set (line, set, ignore_clauses, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  type(c_ptr), intent(in), value :: ignore_clauses  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_clauses
  logical, target :: f_ignore_clauses_native
  logical, pointer :: f_ignore_clauses_native_ptr
  logical(c_bool), pointer :: f_ignore_clauses_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_set 0D_NOT_character
  if (.not. c_associated(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  ! in: f_ignore_clauses 0D_NOT_logical
  if (c_associated(ignore_clauses)) then
    call c_f_pointer(ignore_clauses, f_ignore_clauses_ptr)
    f_ignore_clauses_native = f_ignore_clauses_ptr
    f_ignore_clauses_native_ptr => f_ignore_clauses_native
  else
    f_ignore_clauses_native_ptr => null()
  endif
  f_ix_match = str_first_in_set(f_line, f_set, f_ignore_clauses_native_ptr)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_first_not_in_set (line, set, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_set 0D_NOT_character
  if (.not. c_associated(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  f_ix_match = str_first_not_in_set(f_line, f_set)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_last_in_set (line, set, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_set 0D_NOT_character
  if (.not. c_associated(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  f_ix_match = str_last_in_set(f_line, f_set)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_last_not_in_set (line, set, ix_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_set 0D_NOT_character
  if (.not. c_associated(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  f_ix_match = str_last_not_in_set(f_line, f_set)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_match_wild (str, pat, a_match) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: pat
  character(len=4096), target :: f_pat
  character(kind=c_char), pointer :: f_pat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: a_match  ! 0D_NOT_logical
  logical :: f_a_match
  logical(c_bool), pointer :: f_a_match_ptr
  ! ** End of parameters **
  ! in: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! in: f_pat 0D_NOT_character
  if (.not. c_associated(pat)) return
  call c_f_pointer(pat, f_pat_ptr, [huge(0)])
  call to_f_str(f_pat_ptr, f_pat)
  f_a_match = str_match_wild(f_str, f_pat)

  ! out: f_a_match 0D_NOT_logical
  call c_f_pointer(a_match, f_a_match_ptr)
  f_a_match_ptr = f_a_match
end subroutine
subroutine fortran_str_substitute (string, str_match, str_replace, do_trim, ignore_escaped) &
    bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: str_match
  character(len=4096), target :: f_str_match
  character(kind=c_char), pointer :: f_str_match_ptr(:)
  character(len=4096), pointer :: f_str_match_call_ptr
  type(c_ptr), intent(in), value :: str_replace
  character(len=4096), target :: f_str_replace
  character(kind=c_char), pointer :: f_str_replace_ptr(:)
  character(len=4096), pointer :: f_str_replace_call_ptr
  type(c_ptr), intent(in), value :: do_trim  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_do_trim
  logical, target :: f_do_trim_native
  logical, pointer :: f_do_trim_native_ptr
  logical(c_bool), pointer :: f_do_trim_ptr
  type(c_ptr), intent(in), value :: ignore_escaped  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_escaped
  logical, target :: f_ignore_escaped_native
  logical, pointer :: f_ignore_escaped_native_ptr
  logical(c_bool), pointer :: f_ignore_escaped_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_str_match 0D_NOT_character
  if (c_associated(str_match)) then
    call c_f_pointer(str_match, f_str_match_ptr, [huge(0)])
    call to_f_str(f_str_match_ptr, f_str_match)
    f_str_match_call_ptr => f_str_match
  else
    f_str_match_call_ptr => null()
  endif
  ! in: f_str_replace 0D_NOT_character
  if (c_associated(str_replace)) then
    call c_f_pointer(str_replace, f_str_replace_ptr, [huge(0)])
    call to_f_str(f_str_replace_ptr, f_str_replace)
    f_str_replace_call_ptr => f_str_replace
  else
    f_str_replace_call_ptr => null()
  endif
  ! in: f_do_trim 0D_NOT_logical
  if (c_associated(do_trim)) then
    call c_f_pointer(do_trim, f_do_trim_ptr)
    f_do_trim_native = f_do_trim_ptr
    f_do_trim_native_ptr => f_do_trim_native
  else
    f_do_trim_native_ptr => null()
  endif
  ! in: f_ignore_escaped 0D_NOT_logical
  if (c_associated(ignore_escaped)) then
    call c_f_pointer(ignore_escaped, f_ignore_escaped_ptr)
    f_ignore_escaped_native = f_ignore_escaped_ptr
    f_ignore_escaped_native_ptr => f_ignore_escaped_native
  else
    f_ignore_escaped_native_ptr => null()
  endif
  call str_substitute(f_string, f_str_match_call_ptr, f_str_replace_call_ptr, &
      f_do_trim_native_ptr, f_ignore_escaped_native_ptr)

end subroutine
subroutine fortran_str_upcase (dst, src) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: src
  character(len=4096), target :: f_src
  character(kind=c_char), pointer :: f_src_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: dst
  character(len=4096), target :: f_dst
  character(kind=c_char), pointer :: f_dst_ptr(:)
  ! ** End of parameters **
  ! in: f_src 0D_NOT_character
  if (.not. c_associated(src)) return
  call c_f_pointer(src, f_src_ptr, [huge(0)])
  call to_f_str(f_src_ptr, f_src)
  call str_upcase(f_dst, f_src)

  ! out: f_dst 0D_NOT_character
  call c_f_pointer(dst, f_dst_ptr, [len_trim(f_dst) + 1])
  call to_c_str(f_dst, f_dst_ptr)
end subroutine
subroutine fortran_string_to_int (line, default_, err_flag, err_print_flag, value) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  integer(c_int) :: default_  ! 0D_NOT_integer
  integer :: f_default
  logical(c_bool) :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  type(c_ptr), intent(in), value :: err_print_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_print_flag
  logical, target :: f_err_print_flag_native
  logical, pointer :: f_err_print_flag_native_ptr
  logical(c_bool), pointer :: f_err_print_flag_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_integer
  integer :: f_value
  integer(c_int), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_default 0D_NOT_integer
  f_default = default_
  ! in: f_err_flag 0D_NOT_logical
  f_err_flag = err_flag
  ! in: f_err_print_flag 0D_NOT_logical
  if (c_associated(err_print_flag)) then
    call c_f_pointer(err_print_flag, f_err_print_flag_ptr)
    f_err_print_flag_native = f_err_print_flag_ptr
    f_err_print_flag_native_ptr => f_err_print_flag_native
  else
    f_err_print_flag_native_ptr => null()
  endif
  f_value = string_to_int(f_line, f_default, f_err_flag, f_err_print_flag_native_ptr)

  ! out: f_value 0D_NOT_integer
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_string_to_real (line, default_, err_flag, err_print_flag, value) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  real(c_double) :: default_  ! 0D_NOT_real
  real(rp) :: f_default
  logical(c_bool) :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  type(c_ptr), intent(in), value :: err_print_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_print_flag
  logical, target :: f_err_print_flag_native
  logical, pointer :: f_err_print_flag_native_ptr
  logical(c_bool), pointer :: f_err_print_flag_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! in: f_default 0D_NOT_real
  f_default = default_
  ! in: f_err_flag 0D_NOT_logical
  f_err_flag = err_flag
  ! in: f_err_print_flag 0D_NOT_logical
  if (c_associated(err_print_flag)) then
    call c_f_pointer(err_print_flag, f_err_print_flag_ptr)
    f_err_print_flag_native = f_err_print_flag_ptr
    f_err_print_flag_native_ptr => f_err_print_flag_native
  else
    f_err_print_flag_native_ptr => null()
  endif
  f_value = string_to_real(f_line, f_default, f_err_flag, f_err_print_flag_native_ptr)

  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_string_trim (in_string, out_string, word_len) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_string
  character(len=4096), target :: f_in_string
  character(kind=c_char), pointer :: f_in_string_ptr(:)
  type(c_ptr), intent(in), value :: out_string
  character(len=4096), target :: f_out_string
  character(kind=c_char), pointer :: f_out_string_ptr(:)
  integer(c_int) :: word_len  ! 0D_NOT_integer
  integer :: f_word_len
  ! ** End of parameters **
  ! in: f_in_string 0D_NOT_character
  if (.not. c_associated(in_string)) return
  call c_f_pointer(in_string, f_in_string_ptr, [huge(0)])
  call to_f_str(f_in_string_ptr, f_in_string)
  ! in: f_out_string 0D_NOT_character
  if (.not. c_associated(out_string)) return
  call c_f_pointer(out_string, f_out_string_ptr, [huge(0)])
  call to_f_str(f_out_string_ptr, f_out_string)
  ! in: f_word_len 0D_NOT_integer
  f_word_len = word_len
  call string_trim(f_in_string, f_out_string, f_word_len)

end subroutine
subroutine fortran_string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096), target :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  type(c_ptr), intent(in), value :: delimitors
  character(len=4096), target :: f_delimitors
  character(kind=c_char), pointer :: f_delimitors_ptr(:)
  type(c_ptr), intent(in), value :: out_str
  character(len=4096), target :: f_out_str
  character(kind=c_char), pointer :: f_out_str_ptr(:)
  integer(c_int) :: ix_word  ! 0D_NOT_integer
  integer :: f_ix_word
  type(c_ptr), intent(in), value :: delim
  character(len=4096), target :: f_delim
  character(kind=c_char), pointer :: f_delim_ptr(:)
  integer(c_int) :: ix_next  ! 0D_NOT_integer
  integer :: f_ix_next
  ! ** End of parameters **
  ! in: f_in_str 0D_NOT_character
  if (.not. c_associated(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  ! in: f_delimitors 0D_NOT_character
  if (.not. c_associated(delimitors)) return
  call c_f_pointer(delimitors, f_delimitors_ptr, [huge(0)])
  call to_f_str(f_delimitors_ptr, f_delimitors)
  ! in: f_out_str 0D_NOT_character
  if (.not. c_associated(out_str)) return
  call c_f_pointer(out_str, f_out_str_ptr, [huge(0)])
  call to_f_str(f_out_str_ptr, f_out_str)
  ! in: f_ix_word 0D_NOT_integer
  f_ix_word = ix_word
  ! in: f_delim 0D_NOT_character
  if (.not. c_associated(delim)) return
  call c_f_pointer(delim, f_delim_ptr, [huge(0)])
  call to_f_str(f_delim_ptr, f_delim)
  ! in: f_ix_next 0D_NOT_integer
  f_ix_next = ix_next
  call string_trim2(f_in_str, f_delimitors, f_out_str, f_ix_word, f_delim, f_ix_next)

end subroutine
subroutine fortran_suggest_lmdif (XV, FV, eps, itermx, at_end, reset_flag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: eps  ! 0D_NOT_real
  real(dp) :: f_eps
  integer(c_int) :: itermx  ! 0D_NOT_integer
  integer :: f_itermx
  type(c_ptr), intent(in), value :: reset_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_reset_flag
  logical, target :: f_reset_flag_native
  logical, pointer :: f_reset_flag_native_ptr
  logical(c_bool), pointer :: f_reset_flag_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: at_end  ! 0D_NOT_logical
  logical :: f_at_end
  logical(c_bool), pointer :: f_at_end_ptr
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: XV
  real(dp), pointer :: f_XV(:)
  real(c_double), pointer :: f_XV_ptr(:)
  type(array_descriptor_t), intent(in) :: FV
  real(dp), pointer :: f_FV(:)
  real(c_double), pointer :: f_FV_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(XV%data_ptr)) then
    call c_f_pointer(XV%data_ptr, f_XV_ptr, [XV%dims(1)])
    f_XV => f_XV_ptr
  else
    f_XV => null()
  endif
  !! general array (1D_NOT_real) inout
  if (c_associated(FV%data_ptr)) then
    call c_f_pointer(FV%data_ptr, f_FV_ptr, [FV%dims(1)])
    f_FV => f_FV_ptr
  else
    f_FV => null()
  endif
  ! in: f_eps 0D_NOT_real
  f_eps = eps
  ! in: f_itermx 0D_NOT_integer
  f_itermx = itermx
  ! in: f_reset_flag 0D_NOT_logical
  if (c_associated(reset_flag)) then
    call c_f_pointer(reset_flag, f_reset_flag_ptr)
    f_reset_flag_native = f_reset_flag_ptr
    f_reset_flag_native_ptr => f_reset_flag_native
  else
    f_reset_flag_native_ptr => null()
  endif
  call suggest_lmdif(f_XV, f_FV, f_eps, f_itermx, f_at_end, f_reset_flag_native_ptr)

  ! out: f_at_end 0D_NOT_logical
  call c_f_pointer(at_end, f_at_end_ptr)
  f_at_end_ptr = f_at_end
end subroutine
subroutine fortran_super_bicubic_coef (y, y1, y2, y12, d1, d2, c) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: y
  real(dp) :: f_y(4)
  real(c_double), pointer :: f_y_ptr(:)
  type(array_descriptor_t), intent(in) :: y1
  real(dp) :: f_y1(4)
  real(c_double), pointer :: f_y1_ptr(:)
  type(array_descriptor_t), intent(in) :: y2
  real(dp) :: f_y2(4)
  real(c_double), pointer :: f_y2_ptr(:)
  type(array_descriptor_t), intent(in) :: y12
  real(dp) :: f_y12(4)
  real(c_double), pointer :: f_y12_ptr(:)
  real(c_double) :: d1  ! 0D_NOT_real
  real(dp) :: f_d1
  real(c_double) :: d2  ! 0D_NOT_real
  real(dp) :: f_d2
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: c
  real(dp) :: f_c(4,4)
  real(c_double), pointer :: f_c_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(y%data_ptr)) then
    call c_f_pointer(y%data_ptr, f_y_ptr, [y%dims(1)])
    f_y = f_y_ptr(:)
  else
    f_y_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y1%data_ptr)) then
    call c_f_pointer(y1%data_ptr, f_y1_ptr, [y1%dims(1)])
    f_y1 = f_y1_ptr(:)
  else
    f_y1_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y2%data_ptr)) then
    call c_f_pointer(y2%data_ptr, f_y2_ptr, [y2%dims(1)])
    f_y2 = f_y2_ptr(:)
  else
    f_y2_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y12%data_ptr)) then
    call c_f_pointer(y12%data_ptr, f_y12_ptr, [y12%dims(1)])
    f_y12 = f_y12_ptr(:)
  else
    f_y12_ptr => null()
  endif
  ! in: f_d1 0D_NOT_real
  f_d1 = d1
  ! in: f_d2 0D_NOT_real
  f_d2 = d2
  !! general array (2D_NOT_real) out
  if (c_associated(c%data_ptr)) then
    call c_f_pointer(c%data_ptr, f_c_ptr, [product(c%dims(1:c%rank))])
    ! output-only
  else
    f_c_ptr => null()
  endif
  call super_bicubic_coef(f_y, f_y1, f_y2, f_y12, f_d1, f_d2, f_c)

  ! out: f_c 2D_NOT_real
  if (c_associated(c%data_ptr)) f_c_ptr = mat2vec(f_c, product(c%dims(1:c%rank)))
end subroutine
subroutine fortran_super_bicubic_interpolation (y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, &
    ansy, ansy1, ansy2) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: y
  real(rp) :: f_y(4)
  real(c_double), pointer :: f_y_ptr(:)
  type(array_descriptor_t), intent(in) :: y1
  real(rp) :: f_y1(4)
  real(c_double), pointer :: f_y1_ptr(:)
  type(array_descriptor_t), intent(in) :: y2
  real(rp) :: f_y2(4)
  real(c_double), pointer :: f_y2_ptr(:)
  type(array_descriptor_t), intent(in) :: y12
  real(rp) :: f_y12(4)
  real(c_double), pointer :: f_y12_ptr(:)
  real(c_double) :: x1l  ! 0D_NOT_real
  real(rp) :: f_x1l
  real(c_double) :: x1u  ! 0D_NOT_real
  real(rp) :: f_x1u
  real(c_double) :: x2l  ! 0D_NOT_real
  real(rp) :: f_x2l
  real(c_double) :: x2u  ! 0D_NOT_real
  real(rp) :: f_x2u
  real(c_double) :: x1  ! 0D_NOT_real
  real(rp) :: f_x1
  real(c_double) :: x2  ! 0D_NOT_real
  real(rp) :: f_x2
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ansy  ! 0D_NOT_real
  real(rp) :: f_ansy
  real(c_double), pointer :: f_ansy_ptr
  type(c_ptr), intent(in), value :: ansy1  ! 0D_NOT_real
  real(rp) :: f_ansy1
  real(c_double), pointer :: f_ansy1_ptr
  type(c_ptr), intent(in), value :: ansy2  ! 0D_NOT_real
  real(rp) :: f_ansy2
  real(c_double), pointer :: f_ansy2_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(y%data_ptr)) then
    call c_f_pointer(y%data_ptr, f_y_ptr, [y%dims(1)])
    f_y = f_y_ptr(:)
  else
    f_y_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y1%data_ptr)) then
    call c_f_pointer(y1%data_ptr, f_y1_ptr, [y1%dims(1)])
    f_y1 = f_y1_ptr(:)
  else
    f_y1_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y2%data_ptr)) then
    call c_f_pointer(y2%data_ptr, f_y2_ptr, [y2%dims(1)])
    f_y2 = f_y2_ptr(:)
  else
    f_y2_ptr => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(y12%data_ptr)) then
    call c_f_pointer(y12%data_ptr, f_y12_ptr, [y12%dims(1)])
    f_y12 = f_y12_ptr(:)
  else
    f_y12_ptr => null()
  endif
  ! in: f_x1l 0D_NOT_real
  f_x1l = x1l
  ! in: f_x1u 0D_NOT_real
  f_x1u = x1u
  ! in: f_x2l 0D_NOT_real
  f_x2l = x2l
  ! in: f_x2u 0D_NOT_real
  f_x2u = x2u
  ! in: f_x1 0D_NOT_real
  f_x1 = x1
  ! in: f_x2 0D_NOT_real
  f_x2 = x2
  call super_bicubic_interpolation(f_y, f_y1, f_y2, f_y12, f_x1l, f_x1u, f_x2l, f_x2u, f_x1, &
      f_x2, f_ansy, f_ansy1, f_ansy2)

  ! out: f_ansy 0D_NOT_real
  call c_f_pointer(ansy, f_ansy_ptr)
  f_ansy_ptr = f_ansy
  ! out: f_ansy1 0D_NOT_real
  call c_f_pointer(ansy1, f_ansy1_ptr)
  f_ansy1_ptr = f_ansy1
  ! out: f_ansy2 0D_NOT_real
  call c_f_pointer(ansy2, f_ansy2_ptr)
  f_ansy2_ptr = f_ansy2
end subroutine
subroutine fortran_super_polint (xa, ya, x, y, dy) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: xa
  real(rp), pointer :: f_xa(:)
  real(c_double), pointer :: f_xa_ptr(:)
  type(array_descriptor_t), intent(in) :: ya
  real(rp), pointer :: f_ya(:)
  real(c_double), pointer :: f_ya_ptr(:)
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  type(c_ptr), intent(in), value :: dy  ! 0D_NOT_real
  real(rp) :: f_dy
  real(c_double), pointer :: f_dy_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real) in
  if (c_associated(xa%data_ptr)) then
    call c_f_pointer(xa%data_ptr, f_xa_ptr, [xa%dims(1)])
    f_xa => f_xa_ptr
  else
    f_xa => null()
  endif
  !! general array (1D_NOT_real) in
  if (c_associated(ya%data_ptr)) then
    call c_f_pointer(ya%data_ptr, f_ya_ptr, [ya%dims(1)])
    f_ya => f_ya_ptr
  else
    f_ya => null()
  endif
  ! in: f_x 0D_NOT_real
  f_x = x
  call super_polint(f_xa, f_ya, f_x, f_y, f_dy)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
  ! out: f_dy 0D_NOT_real
  call c_f_pointer(dy, f_dy_ptr)
  f_dy_ptr = f_dy
end subroutine
subroutine fortran_super_poly (x, coeffs, value) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(array_descriptor_t), intent(in) :: coeffs
  real(rp), pointer :: f_coeffs(:)
  real(c_double), pointer :: f_coeffs_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  !! general array (1D_NOT_real) in
  if (c_associated(coeffs%data_ptr)) then
    call c_f_pointer(coeffs%data_ptr, f_coeffs_ptr, [coeffs%dims(1)])
    f_coeffs => f_coeffs_ptr
  else
    f_coeffs => null()
  endif
  f_value = super_poly(f_x, f_coeffs)

  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_super_sobseq (x, ran_state) bind(c)

  use array_desc_mod
  use random_mod, only: random_state_struct
  implicit none
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: x
  real(rp), pointer :: f_x(:)
  real(c_double), pointer :: f_x_ptr(:)
  type(c_ptr), value :: ran_state  ! 0D_NOT_type
  type(random_state_struct), pointer :: f_ran_state
  ! ** End of parameters **
  !! general array (1D_NOT_real) inout
  if (c_associated(x%data_ptr)) then
    call c_f_pointer(x%data_ptr, f_x_ptr, [x%dims(1)])
    f_x => f_x_ptr
  else
    f_x => null()
  endif
  ! inout: f_ran_state 0D_NOT_type
  if (c_associated(ran_state))   call c_f_pointer(ran_state, f_ran_state)
  call super_sobseq(f_x, f_ran_state)

end subroutine
subroutine fortran_super_sort (arr) bind(c)

  use array_desc_mod
  implicit none
  ! ** Inout parameters **
  type(array_descriptor_t), intent(in) :: arr
  integer, pointer :: f_arr(:)
  integer(c_int), pointer :: f_arr_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_integer) inout
  if (c_associated(arr%data_ptr)) then
    call c_f_pointer(arr%data_ptr, f_arr_ptr, [arr%dims(1)])
    f_arr => f_arr_ptr
  else
    f_arr => null()
  endif
  call super_sort(f_arr)

end subroutine
subroutine fortran_system_command (line, err_flag) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_line 0D_NOT_character
  if (.not. c_associated(line)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! out: f_err_flag 0D_NOT_logical
  if (c_associated(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
  else
    f_err_flag_ptr => null()
  endif
  call system_command(f_line, f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_test_tune_tracker_lock (tracker_locked) bind(c)

  use array_desc_mod
  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: tracker_locked
  type(logical_container_alloc), pointer :: f_tracker_locked
  ! ** End of parameters **
  !! container general array (1D_ALLOC_logical)
  if (c_associated(tracker_locked))   call c_f_pointer(tracker_locked, f_tracker_locked)
  call test_tune_tracker_lock(f_tracker_locked%data)

end subroutine
subroutine fortran_test_xgelbd () bind(c)

  use array_desc_mod
  implicit none
  ! ** End of parameters **
  call test_xgelbd()

end subroutine
subroutine fortran_to_str (num, max_signif, string) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  real(c_double) :: num  ! 0D_NOT_real
  real(rp) :: f_num
  type(c_ptr), intent(in), value :: max_signif  ! 0D_NOT_integer
  integer(c_int) :: f_max_signif
  integer(c_int), pointer :: f_max_signif_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** End of parameters **
  ! in: f_num 0D_NOT_real
  f_num = num
  ! in: f_max_signif 0D_NOT_integer
  if (c_associated(max_signif)) then
    call c_f_pointer(max_signif, f_max_signif_ptr)
  else
    f_max_signif_ptr => null()
  endif
  f_string = to_str(f_num, f_max_signif_ptr)

  ! out: f_string 0D_ALLOC_character
  call c_f_pointer(string, f_string_ptr, [len_trim(f_string) + 1])
  call to_c_str(f_string, f_string_ptr)
end subroutine
subroutine fortran_tricubic_cmplx_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz, &
    f_val) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: tricubic_cmplx_coef_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: x_norm  ! 0D_NOT_real
  real(rp) :: f_x_norm
  real(c_double) :: y_norm  ! 0D_NOT_real
  real(rp) :: f_y_norm
  real(c_double) :: z_norm  ! 0D_NOT_real
  real(rp) :: f_z_norm
  type(c_ptr), value :: tri_coef  ! 0D_NOT_type
  type(tricubic_cmplx_coef_struct), pointer :: f_tri_coef
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: df_dx  ! 0D_NOT_complex
  complex(rp) :: f_df_dx
  complex(c_double_complex), pointer :: f_df_dx_ptr
  type(c_ptr), intent(in), value :: df_dy  ! 0D_NOT_complex
  complex(rp) :: f_df_dy
  complex(c_double_complex), pointer :: f_df_dy_ptr
  type(c_ptr), intent(in), value :: df_dz  ! 0D_NOT_complex
  complex(rp) :: f_df_dz
  complex(c_double_complex), pointer :: f_df_dz_ptr
  type(c_ptr), intent(in), value :: f_val  ! 0D_NOT_complex
  complex(rp) :: f_f_val
  complex(c_double_complex), pointer :: f_f_val_ptr
  ! ** End of parameters **
  ! in: f_x_norm 0D_NOT_real
  f_x_norm = x_norm
  ! in: f_y_norm 0D_NOT_real
  f_y_norm = y_norm
  ! in: f_z_norm 0D_NOT_real
  f_z_norm = z_norm
  ! in: f_tri_coef 0D_NOT_type
  if (.not. c_associated(tri_coef)) return
  call c_f_pointer(tri_coef, f_tri_coef)
  ! out: f_df_dx 0D_NOT_complex
  if (c_associated(df_dx)) then
    call c_f_pointer(df_dx, f_df_dx_ptr)
  else
    f_df_dx_ptr => null()
  endif
  ! out: f_df_dy 0D_NOT_complex
  if (c_associated(df_dy)) then
    call c_f_pointer(df_dy, f_df_dy_ptr)
  else
    f_df_dy_ptr => null()
  endif
  ! out: f_df_dz 0D_NOT_complex
  if (c_associated(df_dz)) then
    call c_f_pointer(df_dz, f_df_dz_ptr)
  else
    f_df_dz_ptr => null()
  endif
  f_f_val = tricubic_cmplx_eval(f_x_norm, f_y_norm, f_z_norm, f_tri_coef, f_df_dx, f_df_dy, &
      f_df_dz)

  ! out: f_df_dx 0D_NOT_complex
  call c_f_pointer(df_dx, f_df_dx_ptr)
  f_df_dx_ptr = f_df_dx
  ! out: f_df_dy 0D_NOT_complex
  call c_f_pointer(df_dy, f_df_dy_ptr)
  f_df_dy_ptr = f_df_dy
  ! out: f_df_dz 0D_NOT_complex
  call c_f_pointer(df_dz, f_df_dz_ptr)
  f_df_dz_ptr = f_df_dz
  ! out: f_f_val 0D_NOT_complex
  call c_f_pointer(f_val, f_f_val_ptr)
  f_f_val_ptr = f_f_val
end subroutine
subroutine fortran_tricubic_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz, f_val) &
    bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: tricubic_coef_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: x_norm  ! 0D_NOT_real
  real(rp) :: f_x_norm
  real(c_double) :: y_norm  ! 0D_NOT_real
  real(rp) :: f_y_norm
  real(c_double) :: z_norm  ! 0D_NOT_real
  real(rp) :: f_z_norm
  type(c_ptr), value :: tri_coef  ! 0D_NOT_type
  type(tricubic_coef_struct), pointer :: f_tri_coef
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: df_dx  ! 0D_NOT_real
  real(rp) :: f_df_dx
  real(c_double), pointer :: f_df_dx_ptr
  type(c_ptr), intent(in), value :: df_dy  ! 0D_NOT_real
  real(rp) :: f_df_dy
  real(c_double), pointer :: f_df_dy_ptr
  type(c_ptr), intent(in), value :: df_dz  ! 0D_NOT_real
  real(rp) :: f_df_dz
  real(c_double), pointer :: f_df_dz_ptr
  type(c_ptr), intent(in), value :: f_val  ! 0D_NOT_real
  real(rp) :: f_f_val
  real(c_double), pointer :: f_f_val_ptr
  ! ** End of parameters **
  ! in: f_x_norm 0D_NOT_real
  f_x_norm = x_norm
  ! in: f_y_norm 0D_NOT_real
  f_y_norm = y_norm
  ! in: f_z_norm 0D_NOT_real
  f_z_norm = z_norm
  ! in: f_tri_coef 0D_NOT_type
  if (.not. c_associated(tri_coef)) return
  call c_f_pointer(tri_coef, f_tri_coef)
  ! out: f_df_dx 0D_NOT_real
  if (c_associated(df_dx)) then
    call c_f_pointer(df_dx, f_df_dx_ptr)
  else
    f_df_dx_ptr => null()
  endif
  ! out: f_df_dy 0D_NOT_real
  if (c_associated(df_dy)) then
    call c_f_pointer(df_dy, f_df_dy_ptr)
  else
    f_df_dy_ptr => null()
  endif
  ! out: f_df_dz 0D_NOT_real
  if (c_associated(df_dz)) then
    call c_f_pointer(df_dz, f_df_dz_ptr)
  else
    f_df_dz_ptr => null()
  endif
  f_f_val = tricubic_eval(f_x_norm, f_y_norm, f_z_norm, f_tri_coef, f_df_dx, f_df_dy, f_df_dz)

  ! out: f_df_dx 0D_NOT_real
  call c_f_pointer(df_dx, f_df_dx_ptr)
  f_df_dx_ptr = f_df_dx
  ! out: f_df_dy 0D_NOT_real
  call c_f_pointer(df_dy, f_df_dy_ptr)
  f_df_dy_ptr = f_df_dy
  ! out: f_df_dz 0D_NOT_real
  call c_f_pointer(df_dz, f_df_dz_ptr)
  f_df_dz_ptr = f_df_dz
  ! out: f_f_val 0D_NOT_real
  call c_f_pointer(f_val, f_f_val_ptr)
  f_f_val_ptr = f_f_val
end subroutine
subroutine fortran_tricubic_interpolation_cmplx_coefs (field_at_box, tri_coef) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: cmplx_field_at_3D_box_struct, tricubic_cmplx_coef_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: field_at_box  ! 0D_NOT_type
  type(cmplx_field_at_3D_box_struct), pointer :: f_field_at_box
  ! ** Out parameters **
  type(c_ptr), value :: tri_coef  ! 0D_NOT_type
  type(tricubic_cmplx_coef_struct), pointer :: f_tri_coef
  ! ** End of parameters **
  ! in: f_field_at_box 0D_NOT_type
  if (.not. c_associated(field_at_box)) return
  call c_f_pointer(field_at_box, f_field_at_box)
  ! out: f_tri_coef 0D_NOT_type
  if (.not. c_associated(tri_coef)) return
  call c_f_pointer(tri_coef, f_tri_coef)
  call tricubic_interpolation_cmplx_coefs(f_field_at_box, f_tri_coef)

  ! out: f_tri_coef 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_tricubic_interpolation_coefs (field_at_box, tri_coef) bind(c)

  use array_desc_mod
  use cubic_interpolation_mod, only: field_at_3D_box_struct, tricubic_coef_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: field_at_box  ! 0D_NOT_type
  type(field_at_3D_box_struct), pointer :: f_field_at_box
  ! ** Out parameters **
  type(c_ptr), value :: tri_coef  ! 0D_NOT_type
  type(tricubic_coef_struct), pointer :: f_tri_coef
  ! ** End of parameters **
  ! in: f_field_at_box 0D_NOT_type
  if (.not. c_associated(field_at_box)) return
  call c_f_pointer(field_at_box, f_field_at_box)
  ! out: f_tri_coef 0D_NOT_type
  if (.not. c_associated(tri_coef)) return
  call c_f_pointer(tri_coef, f_tri_coef)
  call tricubic_interpolation_coefs(f_field_at_box, f_tri_coef)

  ! out: f_tri_coef 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_type_this_file (filename) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: filename
  character(len=4096), target :: f_filename
  character(kind=c_char), pointer :: f_filename_ptr(:)
  ! ** End of parameters **
  ! in: f_filename 0D_NOT_character
  if (.not. c_associated(filename)) return
  call c_f_pointer(filename, f_filename_ptr, [huge(0)])
  call to_f_str(f_filename_ptr, f_filename)
  call type_this_file(f_filename)

end subroutine
subroutine fortran_upcase_string (string) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  call upcase_string(f_string)

end subroutine
subroutine fortran_value_of_all_ptr (a_ptr, value) bind(c)

  use array_desc_mod
  use sim_utils_struct, only: all_pointer_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: a_ptr  ! 0D_NOT_type
  type(all_pointer_struct), pointer :: f_a_ptr
  ! ** End of parameters **
  ! inout: f_a_ptr 0D_NOT_type
  if (.not. c_associated(a_ptr)) return
  call c_f_pointer(a_ptr, f_a_ptr)
  f_value = value_of_all_ptr(f_a_ptr)

  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_virtual_memory_usage (usage) bind(c)

  use array_desc_mod
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: usage  ! 0D_NOT_integer
  integer :: f_usage
  integer(c_int), pointer :: f_usage_ptr
  ! ** End of parameters **
  f_usage = virtual_memory_usage()

  ! out: f_usage 0D_NOT_integer
  call c_f_pointer(usage, f_usage_ptr)
  f_usage_ptr = f_usage
end subroutine
subroutine fortran_w_mat_to_axis_angle (w_mat, axis, angle) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  type(c_ptr), intent(in), value :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  real(c_double), pointer :: f_angle_ptr
  ! ** End of parameters **
  !! general array (2D_NOT_real) in
  if (c_associated(w_mat%data_ptr)) then
    call c_f_pointer(w_mat%data_ptr, f_w_mat_ptr, [product(w_mat%dims(1:w_mat%rank))])
    call vec2mat(f_w_mat_ptr, f_w_mat)
  else
    f_w_mat_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    ! output-only
  else
    f_axis_ptr => null()
  endif
  call w_mat_to_axis_angle(f_w_mat, f_axis, f_angle)

  ! out: f_axis 1D_NOT_real
  if (c_associated(axis%data_ptr)) then
    call c_f_pointer(axis%data_ptr, f_axis_ptr, [axis%dims(1)])
    f_axis_ptr = f_axis(:)
  endif
  ! out: f_angle 0D_NOT_real
  call c_f_pointer(angle, f_angle_ptr)
  f_angle_ptr = f_angle
end subroutine
subroutine fortran_w_mat_to_quat (w_mat, quat) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(array_descriptor_t), intent(in) :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** Out parameters **
  type(array_descriptor_t), intent(in) :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** End of parameters **
  !! general array (2D_NOT_real) in
  if (c_associated(w_mat%data_ptr)) then
    call c_f_pointer(w_mat%data_ptr, f_w_mat_ptr, [product(w_mat%dims(1:w_mat%rank))])
    call vec2mat(f_w_mat_ptr, f_w_mat)
  else
    f_w_mat_ptr => null()
  endif
  !! general array (1D_NOT_real) out
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    ! output-only
  else
    f_quat_ptr => null()
  endif
  f_quat = w_mat_to_quat(f_w_mat)

  ! out: f_quat 1D_NOT_real
  if (c_associated(quat%data_ptr)) then
    call c_f_pointer(quat%data_ptr, f_quat_ptr, [quat%dims(1)])
    f_quat_ptr = f_quat(:)
  endif
end subroutine
subroutine fortran_word_len (wording, wlen) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: wording
  character(len=4096), target :: f_wording
  character(kind=c_char), pointer :: f_wording_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: wlen  ! 0D_NOT_integer
  integer :: f_wlen
  integer(c_int), pointer :: f_wlen_ptr
  ! ** End of parameters **
  ! in: f_wording 0D_NOT_character
  if (.not. c_associated(wording)) return
  call c_f_pointer(wording, f_wording_ptr, [huge(0)])
  call to_f_str(f_wording_ptr, f_wording)
  f_wlen = word_len(f_wording)

  ! out: f_wlen 0D_NOT_integer
  call c_f_pointer(wlen, f_wlen_ptr)
  f_wlen_ptr = f_wlen
end subroutine
subroutine fortran_word_read (in_str, delim_list, word, ix_word, delim, delim_found, out_str, &
    ignore_interior) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096), target :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  type(c_ptr), intent(in), value :: delim_list
  character(len=4096), target :: f_delim_list
  character(kind=c_char), pointer :: f_delim_list_ptr(:)
  type(c_ptr), intent(in), value :: word
  character(len=4096), target :: f_word
  character(kind=c_char), pointer :: f_word_ptr(:)
  integer(c_int) :: ix_word  ! 0D_NOT_integer
  integer :: f_ix_word
  type(c_ptr), intent(in), value :: delim
  character(len=4096), target :: f_delim
  character(kind=c_char), pointer :: f_delim_ptr(:)
  logical(c_bool) :: delim_found  ! 0D_NOT_logical
  logical :: f_delim_found
  type(c_ptr), intent(in), value :: out_str
  character(len=4096), target :: f_out_str
  character(kind=c_char), pointer :: f_out_str_ptr(:)
  type(c_ptr), intent(in), value :: ignore_interior  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_interior
  logical, target :: f_ignore_interior_native
  logical, pointer :: f_ignore_interior_native_ptr
  logical(c_bool), pointer :: f_ignore_interior_ptr
  ! ** End of parameters **
  ! in: f_in_str 0D_NOT_character
  if (.not. c_associated(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  ! in: f_delim_list 0D_NOT_character
  if (.not. c_associated(delim_list)) return
  call c_f_pointer(delim_list, f_delim_list_ptr, [huge(0)])
  call to_f_str(f_delim_list_ptr, f_delim_list)
  ! in: f_word 0D_NOT_character
  if (.not. c_associated(word)) return
  call c_f_pointer(word, f_word_ptr, [huge(0)])
  call to_f_str(f_word_ptr, f_word)
  ! in: f_ix_word 0D_NOT_integer
  f_ix_word = ix_word
  ! in: f_delim 0D_NOT_character
  if (.not. c_associated(delim)) return
  call c_f_pointer(delim, f_delim_ptr, [huge(0)])
  call to_f_str(f_delim_ptr, f_delim)
  ! in: f_delim_found 0D_NOT_logical
  f_delim_found = delim_found
  ! in: f_out_str 0D_NOT_character
  if (.not. c_associated(out_str)) return
  call c_f_pointer(out_str, f_out_str_ptr, [huge(0)])
  call to_f_str(f_out_str_ptr, f_out_str)
  ! in: f_ignore_interior 0D_NOT_logical
  if (c_associated(ignore_interior)) then
    call c_f_pointer(ignore_interior, f_ignore_interior_ptr)
    f_ignore_interior_native = f_ignore_interior_ptr
    f_ignore_interior_native_ptr => f_ignore_interior_native
  else
    f_ignore_interior_native_ptr => null()
  endif
  call word_read(f_in_str, f_delim_list, f_word, f_ix_word, f_delim, f_delim_found, f_out_str, &
      f_ignore_interior_native_ptr)

end subroutine
subroutine fortran_x0_radiation_length (species, x0) bind(c)

  use array_desc_mod
  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: x0  ! 0D_NOT_real
  real(rp) :: f_x0
  real(c_double), pointer :: f_x0_ptr
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_x0 = x0_radiation_length(f_species)

  ! out: f_x0 0D_NOT_real
  call c_f_pointer(x0, f_x0_ptr)
  f_x0_ptr = f_x0
end subroutine

end module cppbmad_sim_utils_routines
