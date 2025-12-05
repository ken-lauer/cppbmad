module cppbmad_sim_utils_routines

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

use particle_species_mod, only: anomalous_moment_of, antiparticle, atomic_number, &
    atomic_species_id, charge_of, charge_to_mass_of, is_subatomic_species, mass_of, &
    openpmd_species_name, set_species_charge, species_id, species_id_from_openpmd, &
    species_name, species_of, spin_of, x0_radiation_length

use sim_utils_struct, only: is_false, is_true

use precision_def, only: rp8

use sim_utils_interface, only: asinc, assert_equal, calc_file_number, change_file_number, &
    complex_error_function, cos_one, cosc, cross_product, date_and_time_stamp, detab, &
    display_size_and_resolution, dj_bessel, djb_hash, djb_str_hash, downcase_string, err_exit, &
    factorial, faddeeva_function, fft_1d, file_directorizer, file_get, file_get_open, &
    file_suffixer, find_location_int, find_location_logic, find_location_real, &
    gen_complete_elliptic, get_file_number, get_file_time_stamp, i_bessel, i_bessel_extended, &
    increment_file_number, index_nocase, int_str, is_alphabetic, is_decreasing_sequence, &
    is_increasing_sequence, is_integer, is_logical, is_real, j_bessel, linear_fit, &
    linear_fit_2d, logic_str, lunget, make_legal_comment, match_reg, match_wild, milli_sleep, &
    n_choose_k, n_spline_create, ordinal_str, parse_fortran_format, poly_eval, &
    probability_funct, quadratic_roots, query_string, quote, real_num_fortran_format, &
    real_path, real_str, real_to_string, rms_value, rot_2d, run_timer, set_parameter_int, &
    set_parameter_logic, set_parameter_real, sinc, sincc, sinhx_x, skip_header, sqrt_alpha, &
    sqrt_one, str_count, str_downcase, str_first_in_set, str_first_not_in_set, str_last_in_set, &
    str_last_not_in_set, str_match_wild, str_substitute, str_upcase, string_to_int, &
    string_to_real, string_trim, string_trim2, system_command, to_str, type_this_file, &
    upcase_string, virtual_memory_usage, word_len, word_read

use fourier_mod, only: coarse_frequency_estimate, fine_frequency_estimate, fourier_amplitude

use rotation_3d_mod, only: axis_angle_to_quat, axis_angle_to_w_mat, omega_to_quat, &
    quat_conj_complex, quat_conj_real, quat_inverse, quat_mul_complex, quat_mul_real, &
    quat_rotate_complex, quat_rotate_real, quat_to_axis_angle, quat_to_omega, quat_to_w_mat, &
    rotate_vec, rotate_vec_given_axis_angle, w_mat_to_axis_angle, w_mat_to_quat

use naff_mod, only: interpolated_fft, interpolated_fft_gsl, maximize_projection, naff, projdd

use windowls_mod, only: destfixedwindowls, fixedwindowls, initfixedwindowls

use super_recipes_mod, only: super_bicubic_coef, super_bicubic_interpolation, super_polint, &
    super_poly, super_sort

use random_mod, only: allocate_thread_states, ran_seed_get, ran_seed_put

use cubic_interpolation_mod, only: bicubic_cmplx_eval, tricubic_cmplx_eval

use bin_mod, only: bin_index, bin_x_center, n_bins_automatic

use spline_mod, only: bracket_index_for_spline, create_a_spline, end_akima_spline_calc, &
    reallocate_spline, spline1, spline_akima, spline_akima_interpolate, spline_evaluate

use all_phase_fft, only: apfft, apfft_corr, apfft_ext, hanhan

use bit_mod, only: bit_set


use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

subroutine fortran_antiparticle (species, anti_species) bind(c)

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
  f_anti_species = antiparticle(species=f_species)

  ! out: f_anti_species 0D_NOT_integer
  call c_f_pointer(anti_species, f_anti_species_ptr)
  f_anti_species_ptr = f_anti_species
end subroutine
subroutine fortran_species_of (mass, charge, species) bind(c)

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
  f_species = species_of(mass=f_mass, charge=f_charge)

  ! out: f_species 0D_NOT_integer
  call c_f_pointer(species, f_species_ptr)
  f_species_ptr = f_species
end subroutine
subroutine fortran_species_id (name, default_, print_err, species) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096) :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  type(c_ptr), intent(in), value :: default_  ! 0D_NOT_integer
  integer(c_int) :: f_default
  integer(c_int), pointer :: f_default_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical :: f_print_err_native
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species  ! 0D_NOT_integer
  integer :: f_species
  integer(c_int), pointer :: f_species_ptr
  ! ** End of parameters **
  ! in: f_name 0D_NOT_character
  if (.not. assc(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_default 0D_NOT_integer
  if (assc(default_)) then
    call c_f_pointer(default_, f_default_ptr)
    f_default = f_default_ptr
  else
    ! default_ unset
  endif
  ! in: f_print_err 0D_NOT_logical
  if (assc(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
  else
    ! print_err unset
  endif
  if (assc(default_) .and. assc(print_err)) then

    f_species = species_id(name=f_name, default=f_default, print_err=f_print_err_native)

  elseif (assc(default_)) then

    f_species = species_id(name=f_name, default=f_default)

  elseif (assc(print_err)) then

    f_species = species_id(name=f_name, print_err=f_print_err_native)

  else
    f_species = species_id(name=f_name)

  endif
  ! out: f_species 0D_NOT_integer
  call c_f_pointer(species, f_species_ptr)
  f_species_ptr = f_species
end subroutine
subroutine fortran_atomic_species_id (charge, is_anti, atomic_num, n_nuc, species_id) bind(c)

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
  f_species_id = atomic_species_id(charge=f_charge, is_anti=f_is_anti, atomic_num=f_atomic_num, &
      n_nuc=f_n_nuc)

  ! out: f_species_id 0D_NOT_integer
  call c_f_pointer(species_id, f_species_id_ptr)
  f_species_id_ptr = f_species_id
end subroutine
subroutine fortran_species_name (species, name) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096) :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_name = species_name(species=f_species)

  ! out: f_name 0D_NOT_character
  call c_f_pointer(name, f_name_ptr, [len_trim(f_name) + 1]) ! output-only string
  call to_c_str(f_name, f_name_ptr)
end subroutine
subroutine fortran_species_id_from_openpmd (pmd_name, charge, species) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: pmd_name
  character(len=4096) :: f_pmd_name
  character(kind=c_char), pointer :: f_pmd_name_ptr(:)
  integer(c_int) :: charge  ! 0D_NOT_integer
  integer :: f_charge
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: species  ! 0D_NOT_integer
  integer :: f_species
  integer(c_int), pointer :: f_species_ptr
  ! ** End of parameters **
  ! in: f_pmd_name 0D_NOT_character
  if (.not. assc(pmd_name)) return
  call c_f_pointer(pmd_name, f_pmd_name_ptr, [huge(0)])
  call to_f_str(f_pmd_name_ptr, f_pmd_name)
  ! in: f_charge 0D_NOT_integer
  f_charge = charge
  f_species = species_id_from_openpmd(pmd_name=f_pmd_name, charge=f_charge)

  ! out: f_species 0D_NOT_integer
  call c_f_pointer(species, f_species_ptr)
  f_species_ptr = f_species
end subroutine
subroutine fortran_openpmd_species_name (species, pmd_name) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: species  ! 0D_NOT_integer
  integer :: f_species
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: pmd_name
  character(len=4096) :: f_pmd_name
  character(kind=c_char), pointer :: f_pmd_name_ptr(:)
  ! ** End of parameters **
  ! in: f_species 0D_NOT_integer
  f_species = species
  f_pmd_name = openpmd_species_name(species=f_species)

  ! out: f_pmd_name 0D_NOT_character
  call c_f_pointer(pmd_name, f_pmd_name_ptr, [len_trim(f_pmd_name) + 1]) ! output-only string
  call to_c_str(f_pmd_name, f_pmd_name_ptr)
end subroutine
subroutine fortran_anomalous_moment_of (species, moment) bind(c)

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
  f_moment = anomalous_moment_of(species=f_species)

  ! out: f_moment 0D_NOT_real
  call c_f_pointer(moment, f_moment_ptr)
  f_moment_ptr = f_moment
end subroutine
subroutine fortran_spin_of (species, non_subatomic_default, spin) bind(c)

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
  if (assc(non_subatomic_default)) then
    call c_f_pointer(non_subatomic_default, f_non_subatomic_default_ptr)
    f_non_subatomic_default = f_non_subatomic_default_ptr
  else
    ! non_subatomic_default unset
  endif
  if (assc(non_subatomic_default)) then

    f_spin = spin_of(species=f_species, non_subatomic_default=f_non_subatomic_default)

  else
    f_spin = spin_of(species=f_species)

  endif
  ! out: f_spin 0D_NOT_real
  call c_f_pointer(spin, f_spin_ptr)
  f_spin_ptr = f_spin
end subroutine
subroutine fortran_charge_of (species, default_, charge) bind(c)

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
  if (assc(default_)) then
    call c_f_pointer(default_, f_default_ptr)
    f_default = f_default_ptr
  else
    ! default_ unset
  endif
  if (assc(default_)) then

    f_charge = charge_of(species=f_species, default=f_default)

  else
    f_charge = charge_of(species=f_species)

  endif
  ! out: f_charge 0D_NOT_integer
  call c_f_pointer(charge, f_charge_ptr)
  f_charge_ptr = f_charge
end subroutine
subroutine fortran_mass_of (species, mass) bind(c)

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
  f_mass = mass_of(species=f_species)

  ! out: f_mass 0D_NOT_real
  call c_f_pointer(mass, f_mass_ptr)
  f_mass_ptr = f_mass
end subroutine
subroutine fortran_charge_to_mass_of (species, charge_mass_ratio) bind(c)

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
  f_charge_mass_ratio = charge_to_mass_of(species=f_species)

  ! out: f_charge_mass_ratio 0D_NOT_real
  call c_f_pointer(charge_mass_ratio, f_charge_mass_ratio_ptr)
  f_charge_mass_ratio_ptr = f_charge_mass_ratio
end subroutine
subroutine fortran_set_species_charge (species_in, charge, species_charged) bind(c)

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
  f_species_charged = set_species_charge(species_in=f_species_in, charge=f_charge)

  ! out: f_species_charged 0D_NOT_integer
  call c_f_pointer(species_charged, f_species_charged_ptr)
  f_species_charged_ptr = f_species_charged
end subroutine
subroutine fortran_x0_radiation_length (species, x0) bind(c)

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
  f_x0 = x0_radiation_length(species=f_species)

  ! out: f_x0 0D_NOT_real
  call c_f_pointer(x0, f_x0_ptr)
  f_x0_ptr = f_x0
end subroutine
subroutine fortran_atomic_number (species, atomic_num) bind(c)

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
  f_atomic_num = atomic_number(species=f_species)

  ! out: f_atomic_num 0D_NOT_integer
  call c_f_pointer(atomic_num, f_atomic_num_ptr)
  f_atomic_num_ptr = f_atomic_num
end subroutine
subroutine fortran_is_subatomic_species (species, is_subatomic) bind(c)

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
  f_is_subatomic = is_subatomic_species(species=f_species)

  ! out: f_is_subatomic 0D_NOT_logical
  call c_f_pointer(is_subatomic, f_is_subatomic_ptr)
  f_is_subatomic_ptr = f_is_subatomic
end subroutine
subroutine fortran_is_true (param, this_true) bind(c)

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
  f_this_true = is_true(param=f_param)

  ! out: f_this_true 0D_NOT_logical
  call c_f_pointer(this_true, f_this_true_ptr)
  f_this_true_ptr = f_this_true
end subroutine
subroutine fortran_is_false (param, this_false) bind(c)

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
  f_this_false = is_false(param=f_param)

  ! out: f_this_false 0D_NOT_logical
  call c_f_pointer(this_false, f_this_false_ptr)
  f_this_false_ptr = f_this_false
end subroutine
subroutine fortran_rp8 (int_in, re_out) bind(c)

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
  f_re_out = rp8(int_in=f_int_in)

  ! out: f_re_out 0D_NOT_real
  call c_f_pointer(re_out, f_re_out_ptr)
  f_re_out_ptr = f_re_out
end subroutine
subroutine fortran_set_parameter_real (param_val, set_val, save_val) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: save_val  ! 0D_NOT_real
  real(rp) :: f_save_val
  real(c_double), pointer :: f_save_val_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: param_val  ! 0D_NOT_real
  real(c_double) :: f_param_val
  real(c_double), pointer :: f_param_val_ptr
  type(c_ptr), intent(in), value :: set_val  ! 0D_NOT_real
  real(c_double) :: f_set_val
  real(c_double), pointer :: f_set_val_ptr
  ! ** End of parameters **
  ! inout: f_param_val 0D_NOT_real
  if (assc(param_val)) then
    call c_f_pointer(param_val, f_param_val_ptr)
    f_param_val = f_param_val_ptr
  else
    ! param_val unset
  endif
  ! inout: f_set_val 0D_NOT_real
  if (assc(set_val)) then
    call c_f_pointer(set_val, f_set_val_ptr)
    f_set_val = f_set_val_ptr
  else
    ! set_val unset
  endif
  f_save_val = set_parameter_real(param_val=f_param_val, set_val=f_set_val)

  ! inout: f_param_val 0D_NOT_real
  if (assc(param_val)) then
    call c_f_pointer(param_val, f_param_val_ptr)
    f_param_val_ptr = f_param_val
  endif
  ! inout: f_set_val 0D_NOT_real
  if (assc(set_val)) then
    call c_f_pointer(set_val, f_set_val_ptr)
    f_set_val_ptr = f_set_val
  endif
  ! out: f_save_val 0D_NOT_real
  call c_f_pointer(save_val, f_save_val_ptr)
  f_save_val_ptr = f_save_val
end subroutine
subroutine fortran_set_parameter_int (param_val, set_val, save_val) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: save_val  ! 0D_NOT_integer
  integer :: f_save_val
  integer(c_int), pointer :: f_save_val_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: param_val  ! 0D_NOT_integer
  integer(c_int) :: f_param_val
  integer(c_int), pointer :: f_param_val_ptr
  type(c_ptr), intent(in), value :: set_val  ! 0D_NOT_integer
  integer(c_int) :: f_set_val
  integer(c_int), pointer :: f_set_val_ptr
  ! ** End of parameters **
  ! inout: f_param_val 0D_NOT_integer
  if (assc(param_val)) then
    call c_f_pointer(param_val, f_param_val_ptr)
    f_param_val = f_param_val_ptr
  else
    ! param_val unset
  endif
  ! inout: f_set_val 0D_NOT_integer
  if (assc(set_val)) then
    call c_f_pointer(set_val, f_set_val_ptr)
    f_set_val = f_set_val_ptr
  else
    ! set_val unset
  endif
  f_save_val = set_parameter_int(param_val=f_param_val, set_val=f_set_val)

  ! inout: f_param_val 0D_NOT_integer
  if (assc(param_val)) then
    call c_f_pointer(param_val, f_param_val_ptr)
    f_param_val_ptr = f_param_val
  endif
  ! inout: f_set_val 0D_NOT_integer
  if (assc(set_val)) then
    call c_f_pointer(set_val, f_set_val_ptr)
    f_set_val_ptr = f_set_val
  endif
  ! out: f_save_val 0D_NOT_integer
  call c_f_pointer(save_val, f_save_val_ptr)
  f_save_val_ptr = f_save_val
end subroutine
subroutine fortran_set_parameter_logic (param_val, set_val, save_val) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: save_val  ! 0D_NOT_logical
  logical :: f_save_val
  logical(c_bool), pointer :: f_save_val_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: param_val  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_param_val
  logical :: f_param_val_native
  logical(c_bool), pointer :: f_param_val_ptr
  type(c_ptr), intent(in), value :: set_val  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_set_val
  logical :: f_set_val_native
  logical(c_bool), pointer :: f_set_val_ptr
  ! ** End of parameters **
  ! inout: f_param_val 0D_NOT_logical
  if (assc(param_val)) then
    call c_f_pointer(param_val, f_param_val_ptr)
    f_param_val_native = f_param_val_ptr
  else
    ! param_val unset
  endif
  ! inout: f_set_val 0D_NOT_logical
  if (assc(set_val)) then
    call c_f_pointer(set_val, f_set_val_ptr)
    f_set_val_native = f_set_val_ptr
  else
    ! set_val unset
  endif
  f_save_val = set_parameter_logic(param_val=f_param_val_native, set_val=f_set_val_native)

  ! inout: f_param_val 0D_NOT_logical
  if (assc(param_val)) then
    call c_f_pointer(param_val, f_param_val_ptr)
    f_param_val_ptr = f_param_val_native
  else
    ! f_param_val unset
  endif
  ! inout: f_set_val 0D_NOT_logical
  if (assc(set_val)) then
    call c_f_pointer(set_val, f_set_val_ptr)
    f_set_val_ptr = f_set_val_native
  else
    ! f_set_val unset
  endif
  ! out: f_save_val 0D_NOT_logical
  call c_f_pointer(save_val, f_save_val_ptr)
  f_save_val_ptr = f_save_val
end subroutine
subroutine fortran_cosc (x, nd, y) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** End of parameters **
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd = f_nd_ptr
  else
    ! nd unset
  endif
  if (assc(nd)) then

    f_y = cosc(x=f_x, nd=f_nd)

  else
    f_y = cosc(x=f_x)

  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd_ptr = f_nd
  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_asinc (x, nd, y) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** End of parameters **
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd = f_nd_ptr
  else
    ! nd unset
  endif
  if (assc(nd)) then

    f_y = asinc(x=f_x, nd=f_nd)

  else
    f_y = asinc(x=f_x)

  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd_ptr = f_nd
  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_assert_equal (int_arr, err_str, ival) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ival  ! 0D_NOT_integer
  integer :: f_ival
  integer(c_int), pointer :: f_ival_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: int_arr
  type(integer_container_alloc), pointer :: f_int_arr
  type(c_ptr), intent(in), value :: err_str
  character(len=4096) :: f_err_str
  character(kind=c_char), pointer :: f_err_str_ptr(:)
  ! ** End of parameters **
  !! container general array (1D_ALLOC_integer)
  if (assc(int_arr))   call c_f_pointer(int_arr, f_int_arr)
  ! inout: f_err_str 0D_NOT_character
  if (.not. assc(err_str)) return
  call c_f_pointer(err_str, f_err_str_ptr, [huge(0)])
  call to_f_str(f_err_str_ptr, f_err_str)
  f_ival = assert_equal(int_arr=f_int_arr%data, err_str=f_err_str)

  ! inout: f_err_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_ival 0D_NOT_integer
  call c_f_pointer(ival, f_ival_ptr)
  f_ival_ptr = f_ival
end subroutine
subroutine fortran_calc_file_number (file_name, num_in, num_out, err_flag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096) :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  type(c_ptr), intent(in), value :: num_in  ! 0D_NOT_integer
  integer(c_int) :: f_num_in
  integer(c_int), pointer :: f_num_in_ptr
  type(c_ptr), intent(in), value :: num_out  ! 0D_NOT_integer
  integer(c_int) :: f_num_out
  integer(c_int), pointer :: f_num_out_ptr
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical :: f_err_flag_native
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! inout: f_file_name 0D_NOT_character
  if (.not. assc(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! inout: f_num_in 0D_NOT_integer
  if (assc(num_in)) then
    call c_f_pointer(num_in, f_num_in_ptr)
    f_num_in = f_num_in_ptr
  else
    ! num_in unset
  endif
  ! inout: f_num_out 0D_NOT_integer
  if (assc(num_out)) then
    call c_f_pointer(num_out, f_num_out_ptr)
    f_num_out = f_num_out_ptr
  else
    ! num_out unset
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
  else
    ! err_flag unset
  endif
  call calc_file_number(file_name=f_file_name, num_in=f_num_in, num_out=f_num_out, &
      err_flag=f_err_flag_native)

  ! inout: f_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_num_in 0D_NOT_integer
  if (assc(num_in)) then
    call c_f_pointer(num_in, f_num_in_ptr)
    f_num_in_ptr = f_num_in
  endif
  ! inout: f_num_out 0D_NOT_integer
  if (assc(num_out)) then
    call c_f_pointer(num_out, f_num_out_ptr)
    f_num_out_ptr = f_num_out
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
end subroutine
subroutine fortran_change_file_number (file_name, change) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096) :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  type(c_ptr), intent(in), value :: change  ! 0D_NOT_integer
  integer(c_int) :: f_change
  integer(c_int), pointer :: f_change_ptr
  ! ** End of parameters **
  ! inout: f_file_name 0D_NOT_character
  if (.not. assc(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! inout: f_change 0D_NOT_integer
  if (assc(change)) then
    call c_f_pointer(change, f_change_ptr)
    f_change = f_change_ptr
  else
    ! change unset
  endif
  call change_file_number(file_name=f_file_name, change=f_change)

  ! inout: f_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_change 0D_NOT_integer
  if (assc(change)) then
    call c_f_pointer(change, f_change_ptr)
    f_change_ptr = f_change
  endif
end subroutine
subroutine fortran_cos_one (angle, cos1) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: cos1  ! 0D_NOT_real
  real(rp) :: f_cos1
  real(c_double), pointer :: f_cos1_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: angle  ! 0D_NOT_real
  real(c_double) :: f_angle
  real(c_double), pointer :: f_angle_ptr
  ! ** End of parameters **
  ! inout: f_angle 0D_NOT_real
  if (assc(angle)) then
    call c_f_pointer(angle, f_angle_ptr)
    f_angle = f_angle_ptr
  else
    ! angle unset
  endif
  f_cos1 = cos_one(angle=f_angle)

  ! inout: f_angle 0D_NOT_real
  if (assc(angle)) then
    call c_f_pointer(angle, f_angle_ptr)
    f_angle_ptr = f_angle
  endif
  ! out: f_cos1 0D_NOT_real
  call c_f_pointer(cos1, f_cos1_ptr)
  f_cos1_ptr = f_cos1
end subroutine
subroutine fortran_complex_error_function (wr, wi, zr, zi) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: wr  ! 0D_NOT_real
  real(c_double) :: f_wr
  real(c_double), pointer :: f_wr_ptr
  type(c_ptr), intent(in), value :: wi  ! 0D_NOT_real
  real(c_double) :: f_wi
  real(c_double), pointer :: f_wi_ptr
  type(c_ptr), intent(in), value :: zr  ! 0D_NOT_real
  real(c_double) :: f_zr
  real(c_double), pointer :: f_zr_ptr
  type(c_ptr), intent(in), value :: zi  ! 0D_NOT_real
  real(c_double) :: f_zi
  real(c_double), pointer :: f_zi_ptr
  ! ** End of parameters **
  ! inout: f_wr 0D_NOT_real
  if (assc(wr)) then
    call c_f_pointer(wr, f_wr_ptr)
    f_wr = f_wr_ptr
  else
    ! wr unset
  endif
  ! inout: f_wi 0D_NOT_real
  if (assc(wi)) then
    call c_f_pointer(wi, f_wi_ptr)
    f_wi = f_wi_ptr
  else
    ! wi unset
  endif
  ! inout: f_zr 0D_NOT_real
  if (assc(zr)) then
    call c_f_pointer(zr, f_zr_ptr)
    f_zr = f_zr_ptr
  else
    ! zr unset
  endif
  ! inout: f_zi 0D_NOT_real
  if (assc(zi)) then
    call c_f_pointer(zi, f_zi_ptr)
    f_zi = f_zi_ptr
  else
    ! zi unset
  endif
  call complex_error_function(wr=f_wr, wi=f_wi, zr=f_zr, zi=f_zi)

  ! inout: f_wr 0D_NOT_real
  if (assc(wr)) then
    call c_f_pointer(wr, f_wr_ptr)
    f_wr_ptr = f_wr
  endif
  ! inout: f_wi 0D_NOT_real
  if (assc(wi)) then
    call c_f_pointer(wi, f_wi_ptr)
    f_wi_ptr = f_wi
  endif
  ! inout: f_zr 0D_NOT_real
  if (assc(zr)) then
    call c_f_pointer(zr, f_zr_ptr)
    f_zr_ptr = f_zr
  endif
  ! inout: f_zi 0D_NOT_real
  if (assc(zi)) then
    call c_f_pointer(zi, f_zi_ptr)
    f_zi_ptr = f_zi
  endif
end subroutine
subroutine fortran_cross_product (a, b, c) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: a
  type(real_container_alloc), pointer :: f_a
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: c
  real(rp) :: f_c(3)
  real(c_double), pointer :: f_c_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: b
  type(real_container_alloc), pointer :: f_b
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(a))   call c_f_pointer(a, f_a)
  !! container general array (1D_ALLOC_real)
  if (assc(b))   call c_f_pointer(b, f_b)
  f_c = cross_product(a=f_a%data, b=f_b%data)

  ! out: f_c 1D_NOT_real
  if (assc(c)) f_c_ptr = f_c
end subroutine
subroutine fortran_date_and_time_stamp (string, numeric_month, include_zone) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: numeric_month  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_numeric_month
  logical :: f_numeric_month_native
  logical(c_bool), pointer :: f_numeric_month_ptr
  type(c_ptr), intent(in), value :: include_zone  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_zone
  logical :: f_include_zone_native
  logical(c_bool), pointer :: f_include_zone_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_numeric_month 0D_NOT_logical
  if (assc(numeric_month)) then
    call c_f_pointer(numeric_month, f_numeric_month_ptr)
    f_numeric_month_native = f_numeric_month_ptr
  else
    ! numeric_month unset
  endif
  ! inout: f_include_zone 0D_NOT_logical
  if (assc(include_zone)) then
    call c_f_pointer(include_zone, f_include_zone_ptr)
    f_include_zone_native = f_include_zone_ptr
  else
    ! include_zone unset
  endif
  if (assc(numeric_month) .and. assc(include_zone)) then

    call date_and_time_stamp(string=f_string, numeric_month=f_numeric_month_native, &
        include_zone=f_include_zone_native)

  elseif (assc(numeric_month)) then

    call date_and_time_stamp(string=f_string, numeric_month=f_numeric_month_native)

  elseif (assc(include_zone)) then

    call date_and_time_stamp(string=f_string, include_zone=f_include_zone_native)

  else
    call date_and_time_stamp(string=f_string)

  endif
  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_numeric_month 0D_NOT_logical
  if (assc(numeric_month)) then
    call c_f_pointer(numeric_month, f_numeric_month_ptr)
    f_numeric_month_ptr = f_numeric_month_native
  else
    ! f_numeric_month unset
  endif
  ! inout: f_include_zone 0D_NOT_logical
  if (assc(include_zone)) then
    call c_f_pointer(include_zone, f_include_zone_ptr)
    f_include_zone_ptr = f_include_zone_native
  else
    ! f_include_zone unset
  endif
end subroutine
subroutine fortran_detab (str) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. assc(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  call detab(str=f_str)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res) &
    bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_screen  ! 0D_NOT_integer
  integer(c_int) :: f_ix_screen
  integer(c_int), pointer :: f_ix_screen_ptr
  type(c_ptr), intent(in), value :: x_size  ! 0D_NOT_real
  real(c_double) :: f_x_size
  real(c_double), pointer :: f_x_size_ptr
  type(c_ptr), intent(in), value :: y_size  ! 0D_NOT_real
  real(c_double) :: f_y_size
  real(c_double), pointer :: f_y_size_ptr
  type(c_ptr), intent(in), value :: x_res  ! 0D_NOT_real
  real(c_double) :: f_x_res
  real(c_double), pointer :: f_x_res_ptr
  type(c_ptr), intent(in), value :: y_res  ! 0D_NOT_real
  real(c_double) :: f_y_res
  real(c_double), pointer :: f_y_res_ptr
  ! ** End of parameters **
  ! inout: f_ix_screen 0D_NOT_integer
  if (assc(ix_screen)) then
    call c_f_pointer(ix_screen, f_ix_screen_ptr)
    f_ix_screen = f_ix_screen_ptr
  else
    ! ix_screen unset
  endif
  ! inout: f_x_size 0D_NOT_real
  if (assc(x_size)) then
    call c_f_pointer(x_size, f_x_size_ptr)
    f_x_size = f_x_size_ptr
  else
    ! x_size unset
  endif
  ! inout: f_y_size 0D_NOT_real
  if (assc(y_size)) then
    call c_f_pointer(y_size, f_y_size_ptr)
    f_y_size = f_y_size_ptr
  else
    ! y_size unset
  endif
  ! inout: f_x_res 0D_NOT_real
  if (assc(x_res)) then
    call c_f_pointer(x_res, f_x_res_ptr)
    f_x_res = f_x_res_ptr
  else
    ! x_res unset
  endif
  ! inout: f_y_res 0D_NOT_real
  if (assc(y_res)) then
    call c_f_pointer(y_res, f_y_res_ptr)
    f_y_res = f_y_res_ptr
  else
    ! y_res unset
  endif
  call display_size_and_resolution(ix_screen=f_ix_screen, x_size=f_x_size, y_size=f_y_size, &
      x_res=f_x_res, y_res=f_y_res)

  ! inout: f_ix_screen 0D_NOT_integer
  if (assc(ix_screen)) then
    call c_f_pointer(ix_screen, f_ix_screen_ptr)
    f_ix_screen_ptr = f_ix_screen
  endif
  ! inout: f_x_size 0D_NOT_real
  if (assc(x_size)) then
    call c_f_pointer(x_size, f_x_size_ptr)
    f_x_size_ptr = f_x_size
  endif
  ! inout: f_y_size 0D_NOT_real
  if (assc(y_size)) then
    call c_f_pointer(y_size, f_y_size_ptr)
    f_y_size_ptr = f_y_size
  endif
  ! inout: f_x_res 0D_NOT_real
  if (assc(x_res)) then
    call c_f_pointer(x_res, f_x_res_ptr)
    f_x_res_ptr = f_x_res
  endif
  ! inout: f_y_res 0D_NOT_real
  if (assc(y_res)) then
    call c_f_pointer(y_res, f_y_res_ptr)
    f_y_res_ptr = f_y_res
  endif
end subroutine
subroutine fortran_dj_bessel (m, arg, dj_bes) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: dj_bes  ! 0D_NOT_real
  real(rp) :: f_dj_bes
  real(c_double), pointer :: f_dj_bes_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: m  ! 0D_NOT_integer
  integer(c_int) :: f_m
  integer(c_int), pointer :: f_m_ptr
  type(c_ptr), intent(in), value :: arg  ! 0D_NOT_real
  real(c_double) :: f_arg
  real(c_double), pointer :: f_arg_ptr
  ! ** End of parameters **
  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m = f_m_ptr
  else
    ! m unset
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg = f_arg_ptr
  else
    ! arg unset
  endif
  f_dj_bes = dj_bessel(m=f_m, arg=f_arg)

  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m_ptr = f_m
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg_ptr = f_arg
  endif
  ! out: f_dj_bes 0D_NOT_real
  call c_f_pointer(dj_bes, f_dj_bes_ptr)
  f_dj_bes_ptr = f_dj_bes
end subroutine
subroutine fortran_djb_hash (str, old_hash, hash) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: hash  ! 0D_NOT_integer
  integer :: f_hash
  integer(c_int), pointer :: f_hash_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: old_hash  ! 0D_NOT_integer
  integer(c_int) :: f_old_hash
  integer(c_int), pointer :: f_old_hash_ptr
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. assc(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! inout: f_old_hash 0D_NOT_integer
  if (assc(old_hash)) then
    call c_f_pointer(old_hash, f_old_hash_ptr)
    f_old_hash = f_old_hash_ptr
  else
    ! old_hash unset
  endif
  if (assc(old_hash)) then

    f_hash = djb_hash(str=f_str, old_hash=f_old_hash)

  else
    f_hash = djb_hash(str=f_str)

  endif
  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_old_hash 0D_NOT_integer
  if (assc(old_hash)) then
    call c_f_pointer(old_hash, f_old_hash_ptr)
    f_old_hash_ptr = f_old_hash
  endif
  ! out: f_hash 0D_NOT_integer
  call c_f_pointer(hash, f_hash_ptr)
  f_hash_ptr = f_hash
end subroutine
subroutine fortran_djb_str_hash (in_str, hash_str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: hash_str
  character(len=4096) :: f_hash_str
  character(kind=c_char), pointer :: f_hash_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096) :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  ! ** End of parameters **
  ! inout: f_in_str 0D_NOT_character
  if (.not. assc(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  f_hash_str = djb_str_hash(in_str=f_in_str)

  ! inout: f_in_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_hash_str 0D_NOT_character
  call c_f_pointer(hash_str, f_hash_str_ptr, [len_trim(f_hash_str) + 1]) ! output-only string
  call to_c_str(f_hash_str, f_hash_str_ptr)
end subroutine
subroutine fortran_downcase_string (string) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  call downcase_string(string=f_string)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_err_exit (err_str) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: err_str
  character(len=4096) :: f_err_str
  character(kind=c_char), pointer :: f_err_str_ptr(:)
  ! ** End of parameters **
  ! inout: f_err_str 0D_NOT_character
  if (assc(err_str)) then
    call c_f_pointer(err_str, f_err_str_ptr, [huge(0)])
    call to_f_str(f_err_str_ptr, f_err_str)
  endif
  if (assc(err_str)) then

    call err_exit(err_str=f_err_str)

  else
    call err_exit()

  endif
  ! inout: f_err_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_factorial (n, fact) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: fact  ! 0D_NOT_real
  real(rp) :: f_fact
  real(c_double), pointer :: f_fact_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: n  ! 0D_NOT_integer
  integer(c_int) :: f_n
  integer(c_int), pointer :: f_n_ptr
  ! ** End of parameters **
  ! inout: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n = f_n_ptr
  else
    ! n unset
  endif
  f_fact = factorial(n=f_n)

  ! inout: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n_ptr = f_n
  endif
  ! out: f_fact 0D_NOT_real
  call c_f_pointer(fact, f_fact_ptr)
  f_fact_ptr = f_fact
end subroutine
subroutine fortran_faddeeva_function (z, w, dw) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: z
  real(rp) :: f_z(2)
  real(c_double), pointer :: f_z_ptr(:)
  type(c_ptr), intent(in), value :: w
  real(rp) :: f_w(2)
  real(c_double), pointer :: f_w_ptr(:)
  type(c_ptr), intent(in), value :: dw
  real(rp) :: f_dw(2,2)
  real(c_double), pointer :: f_dw_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(z)) then
    call c_f_pointer(z, f_z_ptr, [2])
    f_z = f_z_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(w)) then
    call c_f_pointer(w, f_w_ptr, [2])
    f_w = f_w_ptr(:)
  endif
  !! general array (2D_NOT_real)
  if (assc(dw)) then
    call c_f_pointer(dw, f_dw_ptr, [2*2])
    call vec2mat(f_dw_ptr, f_dw)
  endif
  call faddeeva_function(z=f_z, w=f_w, dw=f_dw)

end subroutine
subroutine fortran_fft_1d (arr, isign) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: isign  ! 0D_NOT_integer
  integer :: f_isign
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr
  type(complex_container_alloc), pointer :: f_arr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_complex)
  if (assc(arr))   call c_f_pointer(arr, f_arr)
  ! in: f_isign 0D_NOT_integer
  f_isign = isign
  call fft_1d(arr=f_arr%data, isign=f_isign)

end subroutine
subroutine fortran_file_directorizer (in_file, out_file, directory, add_switch) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: in_file
  character(len=4096) :: f_in_file
  character(kind=c_char), pointer :: f_in_file_ptr(:)
  type(c_ptr), intent(in), value :: out_file
  character(len=4096) :: f_out_file
  character(kind=c_char), pointer :: f_out_file_ptr(:)
  type(c_ptr), intent(in), value :: directory
  character(len=4096) :: f_directory
  character(kind=c_char), pointer :: f_directory_ptr(:)
  type(c_ptr), intent(in), value :: add_switch  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_add_switch
  logical :: f_add_switch_native
  logical(c_bool), pointer :: f_add_switch_ptr
  ! ** End of parameters **
  ! inout: f_in_file 0D_NOT_character
  if (.not. assc(in_file)) return
  call c_f_pointer(in_file, f_in_file_ptr, [huge(0)])
  call to_f_str(f_in_file_ptr, f_in_file)
  ! inout: f_out_file 0D_NOT_character
  if (.not. assc(out_file)) return
  call c_f_pointer(out_file, f_out_file_ptr, [huge(0)])
  call to_f_str(f_out_file_ptr, f_out_file)
  ! inout: f_directory 0D_NOT_character
  if (.not. assc(directory)) return
  call c_f_pointer(directory, f_directory_ptr, [huge(0)])
  call to_f_str(f_directory_ptr, f_directory)
  ! inout: f_add_switch 0D_NOT_logical
  if (assc(add_switch)) then
    call c_f_pointer(add_switch, f_add_switch_ptr)
    f_add_switch_native = f_add_switch_ptr
  else
    ! add_switch unset
  endif
  call file_directorizer(in_file=f_in_file, out_file=f_out_file, directory=f_directory, &
      add_switch=f_add_switch_native)

  ! inout: f_in_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_out_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_directory 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_add_switch 0D_NOT_logical
  if (assc(add_switch)) then
    call c_f_pointer(add_switch, f_add_switch_ptr)
    f_add_switch_ptr = f_add_switch_native
  else
    ! f_add_switch unset
  endif
end subroutine
subroutine fortran_file_get (string, dflt_file_name, file_name) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: dflt_file_name
  character(len=4096) :: f_dflt_file_name
  character(kind=c_char), pointer :: f_dflt_file_name_ptr(:)
  type(c_ptr), intent(in), value :: file_name
  character(len=4096) :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_dflt_file_name 0D_NOT_character
  if (.not. assc(dflt_file_name)) return
  call c_f_pointer(dflt_file_name, f_dflt_file_name_ptr, [huge(0)])
  call to_f_str(f_dflt_file_name_ptr, f_dflt_file_name)
  ! inout: f_file_name 0D_NOT_character
  if (.not. assc(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  call file_get(string=f_string, dflt_file_name=f_dflt_file_name, file_name=f_file_name)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_dflt_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_file_get_open (string, dflt_file_name, file_name, file_unit, readonly) &
    bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: dflt_file_name
  character(len=4096) :: f_dflt_file_name
  character(kind=c_char), pointer :: f_dflt_file_name_ptr(:)
  type(c_ptr), intent(in), value :: file_name
  character(len=4096) :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  type(c_ptr), intent(in), value :: file_unit  ! 0D_NOT_integer
  integer(c_int) :: f_file_unit
  integer(c_int), pointer :: f_file_unit_ptr
  type(c_ptr), intent(in), value :: readonly  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_readonly
  logical :: f_readonly_native
  logical(c_bool), pointer :: f_readonly_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_dflt_file_name 0D_NOT_character
  if (.not. assc(dflt_file_name)) return
  call c_f_pointer(dflt_file_name, f_dflt_file_name_ptr, [huge(0)])
  call to_f_str(f_dflt_file_name_ptr, f_dflt_file_name)
  ! inout: f_file_name 0D_NOT_character
  if (.not. assc(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! inout: f_file_unit 0D_NOT_integer
  if (assc(file_unit)) then
    call c_f_pointer(file_unit, f_file_unit_ptr)
    f_file_unit = f_file_unit_ptr
  else
    ! file_unit unset
  endif
  ! inout: f_readonly 0D_NOT_logical
  if (assc(readonly)) then
    call c_f_pointer(readonly, f_readonly_ptr)
    f_readonly_native = f_readonly_ptr
  else
    ! readonly unset
  endif
  call file_get_open(string=f_string, dflt_file_name=f_dflt_file_name, file_name=f_file_name, &
      file_unit=f_file_unit, readonly=f_readonly_native)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_dflt_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_file_unit 0D_NOT_integer
  if (assc(file_unit)) then
    call c_f_pointer(file_unit, f_file_unit_ptr)
    f_file_unit_ptr = f_file_unit
  endif
  ! inout: f_readonly 0D_NOT_logical
  if (assc(readonly)) then
    call c_f_pointer(readonly, f_readonly_ptr)
    f_readonly_ptr = f_readonly_native
  else
    ! f_readonly unset
  endif
end subroutine
subroutine fortran_file_suffixer (in_file_name, out_file_name, suffix, add_switch) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: in_file_name
  character(len=4096) :: f_in_file_name
  character(kind=c_char), pointer :: f_in_file_name_ptr(:)
  type(c_ptr), intent(in), value :: out_file_name
  character(len=4096) :: f_out_file_name
  character(kind=c_char), pointer :: f_out_file_name_ptr(:)
  type(c_ptr), intent(in), value :: suffix
  character(len=4096) :: f_suffix
  character(kind=c_char), pointer :: f_suffix_ptr(:)
  type(c_ptr), intent(in), value :: add_switch  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_add_switch
  logical :: f_add_switch_native
  logical(c_bool), pointer :: f_add_switch_ptr
  ! ** End of parameters **
  ! inout: f_in_file_name 0D_NOT_character
  if (.not. assc(in_file_name)) return
  call c_f_pointer(in_file_name, f_in_file_name_ptr, [huge(0)])
  call to_f_str(f_in_file_name_ptr, f_in_file_name)
  ! inout: f_out_file_name 0D_NOT_character
  if (.not. assc(out_file_name)) return
  call c_f_pointer(out_file_name, f_out_file_name_ptr, [huge(0)])
  call to_f_str(f_out_file_name_ptr, f_out_file_name)
  ! inout: f_suffix 0D_NOT_character
  if (.not. assc(suffix)) return
  call c_f_pointer(suffix, f_suffix_ptr, [huge(0)])
  call to_f_str(f_suffix_ptr, f_suffix)
  ! inout: f_add_switch 0D_NOT_logical
  if (assc(add_switch)) then
    call c_f_pointer(add_switch, f_add_switch_ptr)
    f_add_switch_native = f_add_switch_ptr
  else
    ! add_switch unset
  endif
  call file_suffixer(in_file_name=f_in_file_name, out_file_name=f_out_file_name, &
      suffix=f_suffix, add_switch=f_add_switch_native)

  ! inout: f_in_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_out_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_suffix 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_add_switch 0D_NOT_logical
  if (assc(add_switch)) then
    call c_f_pointer(add_switch, f_add_switch_ptr)
    f_add_switch_ptr = f_add_switch_native
  else
    ! f_add_switch unset
  endif
end subroutine
subroutine fortran_gen_complete_elliptic (kc, p, c, s, err_tol, value) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: kc  ! 0D_NOT_real
  real(c_double) :: f_kc
  real(c_double), pointer :: f_kc_ptr
  type(c_ptr), intent(in), value :: p  ! 0D_NOT_real
  real(c_double) :: f_p
  real(c_double), pointer :: f_p_ptr
  type(c_ptr), intent(in), value :: c  ! 0D_NOT_real
  real(c_double) :: f_c
  real(c_double), pointer :: f_c_ptr
  type(c_ptr), intent(in), value :: s  ! 0D_NOT_real
  real(c_double) :: f_s
  real(c_double), pointer :: f_s_ptr
  type(c_ptr), intent(in), value :: err_tol  ! 0D_NOT_real
  real(c_double) :: f_err_tol
  real(c_double), pointer :: f_err_tol_ptr
  ! ** End of parameters **
  ! inout: f_kc 0D_NOT_real
  if (assc(kc)) then
    call c_f_pointer(kc, f_kc_ptr)
    f_kc = f_kc_ptr
  else
    ! kc unset
  endif
  ! inout: f_p 0D_NOT_real
  if (assc(p)) then
    call c_f_pointer(p, f_p_ptr)
    f_p = f_p_ptr
  else
    ! p unset
  endif
  ! inout: f_c 0D_NOT_real
  if (assc(c)) then
    call c_f_pointer(c, f_c_ptr)
    f_c = f_c_ptr
  else
    ! c unset
  endif
  ! inout: f_s 0D_NOT_real
  if (assc(s)) then
    call c_f_pointer(s, f_s_ptr)
    f_s = f_s_ptr
  else
    ! s unset
  endif
  ! inout: f_err_tol 0D_NOT_real
  if (assc(err_tol)) then
    call c_f_pointer(err_tol, f_err_tol_ptr)
    f_err_tol = f_err_tol_ptr
  else
    ! err_tol unset
  endif
  if (assc(err_tol)) then

    f_value = gen_complete_elliptic(kc=f_kc, p=f_p, c=f_c, s=f_s, err_tol=f_err_tol)

  else
    f_value = gen_complete_elliptic(kc=f_kc, p=f_p, c=f_c, s=f_s)

  endif
  ! inout: f_kc 0D_NOT_real
  if (assc(kc)) then
    call c_f_pointer(kc, f_kc_ptr)
    f_kc_ptr = f_kc
  endif
  ! inout: f_p 0D_NOT_real
  if (assc(p)) then
    call c_f_pointer(p, f_p_ptr)
    f_p_ptr = f_p
  endif
  ! inout: f_c 0D_NOT_real
  if (assc(c)) then
    call c_f_pointer(c, f_c_ptr)
    f_c_ptr = f_c
  endif
  ! inout: f_s 0D_NOT_real
  if (assc(s)) then
    call c_f_pointer(s, f_s_ptr)
    f_s_ptr = f_s
  endif
  ! inout: f_err_tol 0D_NOT_real
  if (assc(err_tol)) then
    call c_f_pointer(err_tol, f_err_tol_ptr)
    f_err_tol_ptr = f_err_tol
  endif
  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_get_file_number (file_name, cnum_in, num_out, err_flag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096) :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  type(c_ptr), intent(in), value :: cnum_in
  character(len=4096) :: f_cnum_in
  character(kind=c_char), pointer :: f_cnum_in_ptr(:)
  type(c_ptr), intent(in), value :: num_out  ! 0D_NOT_integer
  integer(c_int) :: f_num_out
  integer(c_int), pointer :: f_num_out_ptr
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical :: f_err_flag_native
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! inout: f_file_name 0D_NOT_character
  if (.not. assc(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! inout: f_cnum_in 0D_NOT_character
  if (.not. assc(cnum_in)) return
  call c_f_pointer(cnum_in, f_cnum_in_ptr, [huge(0)])
  call to_f_str(f_cnum_in_ptr, f_cnum_in)
  ! inout: f_num_out 0D_NOT_integer
  if (assc(num_out)) then
    call c_f_pointer(num_out, f_num_out_ptr)
    f_num_out = f_num_out_ptr
  else
    ! num_out unset
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
  else
    ! err_flag unset
  endif
  call get_file_number(file_name=f_file_name, cnum_in=f_cnum_in, num_out=f_num_out, &
      err_flag=f_err_flag_native)

  ! inout: f_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_cnum_in 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_num_out 0D_NOT_integer
  if (assc(num_out)) then
    call c_f_pointer(num_out, f_num_out_ptr)
    f_num_out_ptr = f_num_out
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
end subroutine
subroutine fortran_get_file_time_stamp (file, time_stamp) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: file
  character(len=4096) :: f_file
  character(kind=c_char), pointer :: f_file_ptr(:)
  type(c_ptr), intent(in), value :: time_stamp
  character(len=4096) :: f_time_stamp
  character(kind=c_char), pointer :: f_time_stamp_ptr(:)
  ! ** End of parameters **
  ! inout: f_file 0D_NOT_character
  if (.not. assc(file)) return
  call c_f_pointer(file, f_file_ptr, [huge(0)])
  call to_f_str(f_file_ptr, f_file)
  ! inout: f_time_stamp 0D_NOT_character
  if (.not. assc(time_stamp)) return
  call c_f_pointer(time_stamp, f_time_stamp_ptr, [huge(0)])
  call to_f_str(f_time_stamp_ptr, f_time_stamp)
  call get_file_time_stamp(file=f_file, time_stamp=f_time_stamp)

  ! inout: f_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_time_stamp 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_i_bessel (m, arg, i_bes) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: i_bes  ! 0D_NOT_real
  real(rp) :: f_i_bes
  real(c_double), pointer :: f_i_bes_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: m  ! 0D_NOT_integer
  integer(c_int) :: f_m
  integer(c_int), pointer :: f_m_ptr
  type(c_ptr), intent(in), value :: arg  ! 0D_NOT_real
  real(c_double) :: f_arg
  real(c_double), pointer :: f_arg_ptr
  ! ** End of parameters **
  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m = f_m_ptr
  else
    ! m unset
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg = f_arg_ptr
  else
    ! arg unset
  endif
  f_i_bes = i_bessel(m=f_m, arg=f_arg)

  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m_ptr = f_m
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg_ptr = f_arg
  endif
  ! out: f_i_bes 0D_NOT_real
  call c_f_pointer(i_bes, f_i_bes_ptr)
  f_i_bes_ptr = f_i_bes
end subroutine
subroutine fortran_i_bessel_extended (m, arg, i_bes) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: i_bes  ! 0D_NOT_complex
  complex(rp) :: f_i_bes
  complex(c_double_complex), pointer :: f_i_bes_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: m  ! 0D_NOT_integer
  integer(c_int) :: f_m
  integer(c_int), pointer :: f_m_ptr
  type(c_ptr), intent(in), value :: arg  ! 0D_NOT_real
  real(c_double) :: f_arg
  real(c_double), pointer :: f_arg_ptr
  ! ** End of parameters **
  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m = f_m_ptr
  else
    ! m unset
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg = f_arg_ptr
  else
    ! arg unset
  endif
  f_i_bes = i_bessel_extended(m=f_m, arg=f_arg)

  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m_ptr = f_m
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg_ptr = f_arg
  endif
  ! out: f_i_bes 0D_NOT_complex
  call c_f_pointer(i_bes, f_i_bes_ptr)
  f_i_bes_ptr = f_i_bes
end subroutine
subroutine fortran_increment_file_number (file_name, digits, number, cnumber) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096) :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  type(c_ptr), intent(in), value :: digits  ! 0D_NOT_integer
  integer(c_int) :: f_digits
  integer(c_int), pointer :: f_digits_ptr
  type(c_ptr), intent(in), value :: number  ! 0D_NOT_integer
  integer(c_int) :: f_number
  integer(c_int), pointer :: f_number_ptr
  type(c_ptr), intent(in), value :: cnumber
  character(len=4096) :: f_cnumber
  character(kind=c_char), pointer :: f_cnumber_ptr(:)
  ! ** End of parameters **
  ! inout: f_file_name 0D_NOT_character
  if (.not. assc(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! inout: f_digits 0D_NOT_integer
  if (assc(digits)) then
    call c_f_pointer(digits, f_digits_ptr)
    f_digits = f_digits_ptr
  else
    ! digits unset
  endif
  ! inout: f_number 0D_NOT_integer
  if (assc(number)) then
    call c_f_pointer(number, f_number_ptr)
    f_number = f_number_ptr
  else
    ! number unset
  endif
  ! inout: f_cnumber 0D_NOT_character
  if (.not. assc(cnumber)) return
  call c_f_pointer(cnumber, f_cnumber_ptr, [huge(0)])
  call to_f_str(f_cnumber_ptr, f_cnumber)
  call increment_file_number(file_name=f_file_name, digits=f_digits, number=f_number, &
      cnumber=f_cnumber)

  ! inout: f_file_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_digits 0D_NOT_integer
  if (assc(digits)) then
    call c_f_pointer(digits, f_digits_ptr)
    f_digits_ptr = f_digits
  endif
  ! inout: f_number 0D_NOT_integer
  if (assc(number)) then
    call c_f_pointer(number, f_number_ptr)
    f_number_ptr = f_number
  endif
  ! inout: f_cnumber 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_index_nocase (string1, string2, indx) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: indx  ! 0D_NOT_integer
  integer :: f_indx
  integer(c_int), pointer :: f_indx_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string1
  character(len=4096) :: f_string1
  character(kind=c_char), pointer :: f_string1_ptr(:)
  type(c_ptr), intent(in), value :: string2
  character(len=4096) :: f_string2
  character(kind=c_char), pointer :: f_string2_ptr(:)
  ! ** End of parameters **
  ! inout: f_string1 0D_NOT_character
  if (.not. assc(string1)) return
  call c_f_pointer(string1, f_string1_ptr, [huge(0)])
  call to_f_str(f_string1_ptr, f_string1)
  ! inout: f_string2 0D_NOT_character
  if (.not. assc(string2)) return
  call c_f_pointer(string2, f_string2_ptr, [huge(0)])
  call to_f_str(f_string2_ptr, f_string2)
  f_indx = index_nocase(string1=f_string1, string2=f_string2)

  ! inout: f_string1 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_string2 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_indx 0D_NOT_integer
  call c_f_pointer(indx, f_indx_ptr)
  f_indx_ptr = f_indx
end subroutine
subroutine fortran_int_str (int_, width, str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: int_  ! 0D_NOT_integer
  integer(c_int) :: f_int
  integer(c_int), pointer :: f_int_ptr
  type(c_ptr), intent(in), value :: width  ! 0D_NOT_integer
  integer(c_int) :: f_width
  integer(c_int), pointer :: f_width_ptr
  ! ** End of parameters **
  ! inout: f_int 0D_NOT_integer
  if (assc(int_)) then
    call c_f_pointer(int_, f_int_ptr)
    f_int = f_int_ptr
  else
    ! int_ unset
  endif
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width = f_width_ptr
  else
    ! width unset
  endif
  if (assc(width)) then

    f_str = int_str(int=f_int, width=f_width)

  else
    f_str = int_str(int=f_int)

  endif
  ! inout: f_int 0D_NOT_integer
  if (assc(int_)) then
    call c_f_pointer(int_, f_int_ptr)
    f_int_ptr = f_int
  endif
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width_ptr = f_width
  endif
  ! out: f_str 0D_ALLOC_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1]) ! output-only string
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_is_alphabetic (string, valid_chars, is_alpha) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_alpha  ! 0D_NOT_logical
  logical :: f_is_alpha
  logical(c_bool), pointer :: f_is_alpha_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: valid_chars
  character(len=4096) :: f_valid_chars
  character(kind=c_char), pointer :: f_valid_chars_ptr(:)
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_valid_chars 0D_NOT_character
  if (assc(valid_chars)) then
    call c_f_pointer(valid_chars, f_valid_chars_ptr, [huge(0)])
    call to_f_str(f_valid_chars_ptr, f_valid_chars)
  endif
  if (assc(valid_chars)) then

    f_is_alpha = is_alphabetic(string=f_string, valid_chars=f_valid_chars)

  else
    f_is_alpha = is_alphabetic(string=f_string)

  endif
  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_valid_chars 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_is_alpha 0D_NOT_logical
  call c_f_pointer(is_alpha, f_is_alpha_ptr)
  f_is_alpha_ptr = f_is_alpha
end subroutine
subroutine fortran_is_decreasing_sequence (array, strict, is_decreasing) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: array
  type(real_container_alloc), pointer :: f_array
  type(c_ptr), intent(in), value :: strict  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_strict
  logical :: f_strict_native
  logical(c_bool), pointer :: f_strict_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_decreasing  ! 0D_NOT_logical
  logical :: f_is_decreasing
  logical(c_bool), pointer :: f_is_decreasing_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(array))   call c_f_pointer(array, f_array)
  ! in: f_strict 0D_NOT_logical
  if (assc(strict)) then
    call c_f_pointer(strict, f_strict_ptr)
    f_strict_native = f_strict_ptr
  else
    ! strict unset
  endif
  if (assc(strict)) then

    f_is_decreasing = is_decreasing_sequence(array=f_array%data, strict=f_strict_native)

  else
    f_is_decreasing = is_decreasing_sequence(array=f_array%data)

  endif
  ! out: f_is_decreasing 0D_NOT_logical
  call c_f_pointer(is_decreasing, f_is_decreasing_ptr)
  f_is_decreasing_ptr = f_is_decreasing
end subroutine
subroutine fortran_is_increasing_sequence (array, strict, is_increasing) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: array
  type(real_container_alloc), pointer :: f_array
  type(c_ptr), intent(in), value :: strict  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_strict
  logical :: f_strict_native
  logical(c_bool), pointer :: f_strict_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_increasing  ! 0D_NOT_logical
  logical :: f_is_increasing
  logical(c_bool), pointer :: f_is_increasing_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(array))   call c_f_pointer(array, f_array)
  ! in: f_strict 0D_NOT_logical
  if (assc(strict)) then
    call c_f_pointer(strict, f_strict_ptr)
    f_strict_native = f_strict_ptr
  else
    ! strict unset
  endif
  if (assc(strict)) then

    f_is_increasing = is_increasing_sequence(array=f_array%data, strict=f_strict_native)

  else
    f_is_increasing = is_increasing_sequence(array=f_array%data)

  endif
  ! out: f_is_increasing 0D_NOT_logical
  call c_f_pointer(is_increasing, f_is_increasing_ptr)
  f_is_increasing_ptr = f_is_increasing
end subroutine
subroutine fortran_is_integer (string, int_, delims, ix_word, valid) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: int_  ! 0D_NOT_integer
  integer(c_int) :: f_int
  integer(c_int), pointer :: f_int_ptr
  type(c_ptr), intent(in), value :: delims
  character(len=4096) :: f_delims
  character(kind=c_char), pointer :: f_delims_ptr(:)
  type(c_ptr), intent(in), value :: ix_word  ! 0D_NOT_integer
  integer(c_int) :: f_ix_word
  integer(c_int), pointer :: f_ix_word_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_int 0D_NOT_integer
  if (assc(int_)) then
    call c_f_pointer(int_, f_int_ptr)
    f_int = f_int_ptr
  else
    ! int_ unset
  endif
  ! inout: f_delims 0D_NOT_character
  if (assc(delims)) then
    call c_f_pointer(delims, f_delims_ptr, [huge(0)])
    call to_f_str(f_delims_ptr, f_delims)
  endif
  ! inout: f_ix_word 0D_NOT_integer
  if (assc(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
    f_ix_word = f_ix_word_ptr
  else
    ! ix_word unset
  endif
  if (assc(int_) .and. assc(delims) .and. assc(ix_word)) then

    f_valid = is_integer(string=f_string, int=f_int, delims=f_delims, ix_word=f_ix_word)

  elseif (assc(int_) .and. assc(delims)) then

    f_valid = is_integer(string=f_string, int=f_int, delims=f_delims)

  elseif (assc(int_) .and. assc(ix_word)) then

    f_valid = is_integer(string=f_string, int=f_int, ix_word=f_ix_word)

  elseif (assc(delims) .and. assc(ix_word)) then

    f_valid = is_integer(string=f_string, delims=f_delims, ix_word=f_ix_word)

  elseif (assc(int_)) then

    f_valid = is_integer(string=f_string, int=f_int)

  elseif (assc(delims)) then

    f_valid = is_integer(string=f_string, delims=f_delims)

  elseif (assc(ix_word)) then

    f_valid = is_integer(string=f_string, ix_word=f_ix_word)

  else
    f_valid = is_integer(string=f_string)

  endif
  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_int 0D_NOT_integer
  if (assc(int_)) then
    call c_f_pointer(int_, f_int_ptr)
    f_int_ptr = f_int
  endif
  ! inout: f_delims 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ix_word 0D_NOT_integer
  if (assc(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
    f_ix_word_ptr = f_ix_word
  endif
  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
end subroutine
subroutine fortran_is_logical (string, ignore, valid) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: ignore  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore
  logical :: f_ignore_native
  logical(c_bool), pointer :: f_ignore_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_ignore 0D_NOT_logical
  if (assc(ignore)) then
    call c_f_pointer(ignore, f_ignore_ptr)
    f_ignore_native = f_ignore_ptr
  else
    ! ignore unset
  endif
  if (assc(ignore)) then

    f_valid = is_logical(string=f_string, ignore=f_ignore_native)

  else
    f_valid = is_logical(string=f_string)

  endif
  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ignore 0D_NOT_logical
  if (assc(ignore)) then
    call c_f_pointer(ignore, f_ignore_ptr)
    f_ignore_ptr = f_ignore_native
  else
    ! f_ignore unset
  endif
  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
end subroutine
subroutine fortran_is_real (string, ignore, real_num, valid) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: ignore  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore
  logical :: f_ignore_native
  logical(c_bool), pointer :: f_ignore_ptr
  type(c_ptr), intent(in), value :: real_num  ! 0D_NOT_real
  real(c_double) :: f_real_num
  real(c_double), pointer :: f_real_num_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_ignore 0D_NOT_logical
  if (assc(ignore)) then
    call c_f_pointer(ignore, f_ignore_ptr)
    f_ignore_native = f_ignore_ptr
  else
    ! ignore unset
  endif
  ! inout: f_real_num 0D_NOT_real
  if (assc(real_num)) then
    call c_f_pointer(real_num, f_real_num_ptr)
    f_real_num = f_real_num_ptr
  else
    ! real_num unset
  endif
  if (assc(ignore) .and. assc(real_num)) then

    f_valid = is_real(string=f_string, ignore=f_ignore_native, real_num=f_real_num)

  elseif (assc(ignore)) then

    f_valid = is_real(string=f_string, ignore=f_ignore_native)

  elseif (assc(real_num)) then

    f_valid = is_real(string=f_string, real_num=f_real_num)

  else
    f_valid = is_real(string=f_string)

  endif
  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ignore 0D_NOT_logical
  if (assc(ignore)) then
    call c_f_pointer(ignore, f_ignore_ptr)
    f_ignore_ptr = f_ignore_native
  else
    ! f_ignore unset
  endif
  ! inout: f_real_num 0D_NOT_real
  if (assc(real_num)) then
    call c_f_pointer(real_num, f_real_num_ptr)
    f_real_num_ptr = f_real_num
  endif
  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
end subroutine
subroutine fortran_j_bessel (m, arg, j_bes) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: j_bes  ! 0D_NOT_real
  real(rp) :: f_j_bes
  real(c_double), pointer :: f_j_bes_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: m  ! 0D_NOT_integer
  integer(c_int) :: f_m
  integer(c_int), pointer :: f_m_ptr
  type(c_ptr), intent(in), value :: arg  ! 0D_NOT_real
  real(c_double) :: f_arg
  real(c_double), pointer :: f_arg_ptr
  ! ** End of parameters **
  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m = f_m_ptr
  else
    ! m unset
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg = f_arg_ptr
  else
    ! arg unset
  endif
  f_j_bes = j_bessel(m=f_m, arg=f_arg)

  ! inout: f_m 0D_NOT_integer
  if (assc(m)) then
    call c_f_pointer(m, f_m_ptr)
    f_m_ptr = f_m
  endif
  ! inout: f_arg 0D_NOT_real
  if (assc(arg)) then
    call c_f_pointer(arg, f_arg_ptr)
    f_arg_ptr = f_arg
  endif
  ! out: f_j_bes 0D_NOT_real
  call c_f_pointer(j_bes, f_j_bes_ptr)
  f_j_bes_ptr = f_j_bes
end subroutine
subroutine fortran_linear_fit (x, y, n_data, a, b, sig_a, sig_b) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x
  type(real_container_alloc), pointer :: f_x
  type(c_ptr), intent(in), value :: y
  type(real_container_alloc), pointer :: f_y
  type(c_ptr), intent(in), value :: n_data  ! 0D_NOT_integer
  integer(c_int) :: f_n_data
  integer(c_int), pointer :: f_n_data_ptr
  type(c_ptr), intent(in), value :: a  ! 0D_NOT_real
  real(c_double) :: f_a
  real(c_double), pointer :: f_a_ptr
  type(c_ptr), intent(in), value :: b  ! 0D_NOT_real
  real(c_double) :: f_b
  real(c_double), pointer :: f_b_ptr
  type(c_ptr), intent(in), value :: sig_a  ! 0D_NOT_real
  real(c_double) :: f_sig_a
  real(c_double), pointer :: f_sig_a_ptr
  type(c_ptr), intent(in), value :: sig_b  ! 0D_NOT_real
  real(c_double) :: f_sig_b
  real(c_double), pointer :: f_sig_b_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(x))   call c_f_pointer(x, f_x)
  !! container general array (1D_ALLOC_real)
  if (assc(y))   call c_f_pointer(y, f_y)
  ! inout: f_n_data 0D_NOT_integer
  if (assc(n_data)) then
    call c_f_pointer(n_data, f_n_data_ptr)
    f_n_data = f_n_data_ptr
  else
    ! n_data unset
  endif
  ! inout: f_a 0D_NOT_real
  if (assc(a)) then
    call c_f_pointer(a, f_a_ptr)
    f_a = f_a_ptr
  else
    ! a unset
  endif
  ! inout: f_b 0D_NOT_real
  if (assc(b)) then
    call c_f_pointer(b, f_b_ptr)
    f_b = f_b_ptr
  else
    ! b unset
  endif
  ! inout: f_sig_a 0D_NOT_real
  if (assc(sig_a)) then
    call c_f_pointer(sig_a, f_sig_a_ptr)
    f_sig_a = f_sig_a_ptr
  else
    ! sig_a unset
  endif
  ! inout: f_sig_b 0D_NOT_real
  if (assc(sig_b)) then
    call c_f_pointer(sig_b, f_sig_b_ptr)
    f_sig_b = f_sig_b_ptr
  else
    ! sig_b unset
  endif
  call linear_fit(x=f_x%data, y=f_y%data, n_data=f_n_data, a=f_a, b=f_b, sig_a=f_sig_a, &
      sig_b=f_sig_b)

  ! inout: f_n_data 0D_NOT_integer
  if (assc(n_data)) then
    call c_f_pointer(n_data, f_n_data_ptr)
    f_n_data_ptr = f_n_data
  endif
  ! inout: f_a 0D_NOT_real
  if (assc(a)) then
    call c_f_pointer(a, f_a_ptr)
    f_a_ptr = f_a
  endif
  ! inout: f_b 0D_NOT_real
  if (assc(b)) then
    call c_f_pointer(b, f_b_ptr)
    f_b_ptr = f_b
  endif
  ! inout: f_sig_a 0D_NOT_real
  if (assc(sig_a)) then
    call c_f_pointer(sig_a, f_sig_a_ptr)
    f_sig_a_ptr = f_sig_a
  endif
  ! inout: f_sig_b 0D_NOT_real
  if (assc(sig_b)) then
    call c_f_pointer(sig_b, f_sig_b_ptr)
    f_sig_b_ptr = f_sig_b
  endif
end subroutine
subroutine fortran_linear_fit_2d (x, y, z, coef) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: x
  type(real_container_alloc), pointer :: f_x
  type(c_ptr), intent(in), value :: y
  type(real_container_alloc), pointer :: f_y
  type(c_ptr), intent(in), value :: z
  type(real_container_alloc), pointer :: f_z
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: coef
  real(rp) :: f_coef(3)
  real(c_double), pointer :: f_coef_ptr(:)
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(x))   call c_f_pointer(x, f_x)
  !! container general array (1D_ALLOC_real)
  if (assc(y))   call c_f_pointer(y, f_y)
  !! container general array (1D_ALLOC_real)
  if (assc(z))   call c_f_pointer(z, f_z)
  call linear_fit_2d(x=f_x%data, y=f_y%data, z=f_z%data, coef=f_coef)

  ! out: f_coef 1D_NOT_real
  if (assc(coef)) f_coef_ptr = f_coef
end subroutine
subroutine fortran_logic_str (logic, str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: logic  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_logic
  logical :: f_logic_native
  logical(c_bool), pointer :: f_logic_ptr
  ! ** End of parameters **
  ! inout: f_logic 0D_NOT_logical
  if (assc(logic)) then
    call c_f_pointer(logic, f_logic_ptr)
    f_logic_native = f_logic_ptr
  else
    ! logic unset
  endif
  f_str = logic_str(logic=f_logic_native)

  ! inout: f_logic 0D_NOT_logical
  if (assc(logic)) then
    call c_f_pointer(logic, f_logic_ptr)
    f_logic_ptr = f_logic_native
  else
    ! f_logic unset
  endif
  ! out: f_str 0D_NOT_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1]) ! output-only string
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_lunget (func_retval__) bind(c)

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
subroutine fortran_match_reg (str, pat, is_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_match  ! 0D_NOT_logical
  logical :: f_is_match
  logical(c_bool), pointer :: f_is_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: pat
  character(len=4096) :: f_pat
  character(kind=c_char), pointer :: f_pat_ptr(:)
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. assc(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! inout: f_pat 0D_NOT_character
  if (.not. assc(pat)) return
  call c_f_pointer(pat, f_pat_ptr, [huge(0)])
  call to_f_str(f_pat_ptr, f_pat)
  f_is_match = match_reg(str=f_str, pat=f_pat)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_pat 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_is_match 0D_NOT_logical
  call c_f_pointer(is_match, f_is_match_ptr)
  f_is_match_ptr = f_is_match
end subroutine
subroutine fortran_milli_sleep (milli_sec) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: milli_sec  ! 0D_NOT_integer
  integer(c_int) :: f_milli_sec
  integer(c_int), pointer :: f_milli_sec_ptr
  ! ** End of parameters **
  ! inout: f_milli_sec 0D_NOT_integer
  if (assc(milli_sec)) then
    call c_f_pointer(milli_sec, f_milli_sec_ptr)
    f_milli_sec = f_milli_sec_ptr
  else
    ! milli_sec unset
  endif
  call milli_sleep(milli_sec=f_milli_sec)

  ! inout: f_milli_sec 0D_NOT_integer
  if (assc(milli_sec)) then
    call c_f_pointer(milli_sec, f_milli_sec_ptr)
    f_milli_sec_ptr = f_milli_sec
  endif
end subroutine
subroutine fortran_make_legal_comment (comment_in, comment_out) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: comment_in
  character(len=4096) :: f_comment_in
  character(kind=c_char), pointer :: f_comment_in_ptr(:)
  type(c_ptr), intent(in), value :: comment_out
  character(len=4096) :: f_comment_out
  character(kind=c_char), pointer :: f_comment_out_ptr(:)
  ! ** End of parameters **
  ! inout: f_comment_in 0D_NOT_character
  if (.not. assc(comment_in)) return
  call c_f_pointer(comment_in, f_comment_in_ptr, [huge(0)])
  call to_f_str(f_comment_in_ptr, f_comment_in)
  ! inout: f_comment_out 0D_NOT_character
  if (.not. assc(comment_out)) return
  call c_f_pointer(comment_out, f_comment_out_ptr, [huge(0)])
  call to_f_str(f_comment_out_ptr, f_comment_out)
  call make_legal_comment(comment_in=f_comment_in, comment_out=f_comment_out)

  ! inout: f_comment_in 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_comment_out 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_match_wild (string, template_, is_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_match  ! 0D_NOT_logical
  logical :: f_is_match
  logical(c_bool), pointer :: f_is_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: template_
  character(len=4096) :: f_template
  character(kind=c_char), pointer :: f_template_ptr(:)
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_template 0D_NOT_character
  if (.not. assc(template_)) return
  call c_f_pointer(template_, f_template_ptr, [huge(0)])
  call to_f_str(f_template_ptr, f_template)
  f_is_match = match_wild(string=f_string, template=f_template)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_template 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_is_match 0D_NOT_logical
  call c_f_pointer(is_match, f_is_match_ptr)
  f_is_match_ptr = f_is_match
end subroutine
subroutine fortran_n_choose_k (n, k, nck) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: nck  ! 0D_NOT_real
  real(rp) :: f_nck
  real(c_double), pointer :: f_nck_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: n  ! 0D_NOT_integer
  integer(c_int) :: f_n
  integer(c_int), pointer :: f_n_ptr
  type(c_ptr), intent(in), value :: k  ! 0D_NOT_integer
  integer(c_int) :: f_k
  integer(c_int), pointer :: f_k_ptr
  ! ** End of parameters **
  ! inout: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n = f_n_ptr
  else
    ! n unset
  endif
  ! inout: f_k 0D_NOT_integer
  if (assc(k)) then
    call c_f_pointer(k, f_k_ptr)
    f_k = f_k_ptr
  else
    ! k unset
  endif
  f_nck = n_choose_k(n=f_n, k=f_k)

  ! inout: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n_ptr = f_n
  endif
  ! inout: f_k 0D_NOT_integer
  if (assc(k)) then
    call c_f_pointer(k, f_k_ptr)
    f_k_ptr = f_k
  endif
  ! out: f_nck 0D_NOT_real
  call c_f_pointer(nck, f_nck_ptr)
  f_nck_ptr = f_nck
end subroutine
subroutine fortran_n_spline_create (deriv0, deriv1, x1, n_spline) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: deriv0
  type(real_container_alloc), pointer :: f_deriv0
  type(c_ptr), intent(in), value :: deriv1
  type(real_container_alloc), pointer :: f_deriv1
  real(c_double) :: x1  ! 0D_NOT_real
  real(rp) :: f_x1
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: n_spline
  type(real_container_alloc), pointer :: f_n_spline
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(deriv0))   call c_f_pointer(deriv0, f_deriv0)
  !! container general array (1D_ALLOC_real)
  if (assc(deriv1))   call c_f_pointer(deriv1, f_deriv1)
  ! in: f_x1 0D_NOT_real
  f_x1 = x1
  !! container general array (1D_ALLOC_real)
  if (assc(n_spline))   call c_f_pointer(n_spline, f_n_spline)
  call n_spline_create(deriv0=f_deriv0%data, deriv1=f_deriv1%data, x1=f_x1, &
      n_spline=f_n_spline%data)

end subroutine
subroutine fortran_ordinal_str (n, str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: n  ! 0D_NOT_integer
  integer(c_int) :: f_n
  integer(c_int), pointer :: f_n_ptr
  ! ** End of parameters **
  ! inout: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n = f_n_ptr
  else
    ! n unset
  endif
  f_str = ordinal_str(n=f_n)

  ! inout: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n_ptr = f_n
  endif
  ! out: f_str 0D_ALLOC_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1]) ! output-only string
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_parse_fortran_format (format_str, n_repeat, power, descrip, width, digits) &
    bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: format_str
  character(len=4096) :: f_format_str
  character(kind=c_char), pointer :: f_format_str_ptr(:)
  type(c_ptr), intent(in), value :: n_repeat  ! 0D_NOT_integer
  integer(c_int) :: f_n_repeat
  integer(c_int), pointer :: f_n_repeat_ptr
  type(c_ptr), intent(in), value :: power  ! 0D_NOT_integer
  integer(c_int) :: f_power
  integer(c_int), pointer :: f_power_ptr
  type(c_ptr), intent(in), value :: descrip
  character(len=4096) :: f_descrip
  character(kind=c_char), pointer :: f_descrip_ptr(:)
  type(c_ptr), intent(in), value :: width  ! 0D_NOT_integer
  integer(c_int) :: f_width
  integer(c_int), pointer :: f_width_ptr
  type(c_ptr), intent(in), value :: digits  ! 0D_NOT_integer
  integer(c_int) :: f_digits
  integer(c_int), pointer :: f_digits_ptr
  ! ** End of parameters **
  ! inout: f_format_str 0D_NOT_character
  if (.not. assc(format_str)) return
  call c_f_pointer(format_str, f_format_str_ptr, [huge(0)])
  call to_f_str(f_format_str_ptr, f_format_str)
  ! inout: f_n_repeat 0D_NOT_integer
  if (assc(n_repeat)) then
    call c_f_pointer(n_repeat, f_n_repeat_ptr)
    f_n_repeat = f_n_repeat_ptr
  else
    ! n_repeat unset
  endif
  ! inout: f_power 0D_NOT_integer
  if (assc(power)) then
    call c_f_pointer(power, f_power_ptr)
    f_power = f_power_ptr
  else
    ! power unset
  endif
  ! inout: f_descrip 0D_NOT_character
  if (.not. assc(descrip)) return
  call c_f_pointer(descrip, f_descrip_ptr, [huge(0)])
  call to_f_str(f_descrip_ptr, f_descrip)
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width = f_width_ptr
  else
    ! width unset
  endif
  ! inout: f_digits 0D_NOT_integer
  if (assc(digits)) then
    call c_f_pointer(digits, f_digits_ptr)
    f_digits = f_digits_ptr
  else
    ! digits unset
  endif
  call parse_fortran_format(format_str=f_format_str, n_repeat=f_n_repeat, power=f_power, &
      descrip=f_descrip, width=f_width, digits=f_digits)

  ! inout: f_format_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_n_repeat 0D_NOT_integer
  if (assc(n_repeat)) then
    call c_f_pointer(n_repeat, f_n_repeat_ptr)
    f_n_repeat_ptr = f_n_repeat
  endif
  ! inout: f_power 0D_NOT_integer
  if (assc(power)) then
    call c_f_pointer(power, f_power_ptr)
    f_power_ptr = f_power
  endif
  ! inout: f_descrip 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width_ptr = f_width
  endif
  ! inout: f_digits 0D_NOT_integer
  if (assc(digits)) then
    call c_f_pointer(digits, f_digits_ptr)
    f_digits_ptr = f_digits
  endif
end subroutine
subroutine fortran_poly_eval (poly, x, diff_coef, y) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: poly
  type(real_container_alloc), pointer :: f_poly
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: diff_coef  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_diff_coef
  logical :: f_diff_coef_native
  logical(c_bool), pointer :: f_diff_coef_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(poly))   call c_f_pointer(poly, f_poly)
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_diff_coef 0D_NOT_logical
  if (assc(diff_coef)) then
    call c_f_pointer(diff_coef, f_diff_coef_ptr)
    f_diff_coef_native = f_diff_coef_ptr
  else
    ! diff_coef unset
  endif
  if (assc(diff_coef)) then

    f_y = poly_eval(poly=f_poly%data, x=f_x, diff_coef=f_diff_coef_native)

  else
    f_y = poly_eval(poly=f_poly%data, x=f_x)

  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_probability_funct (x, prob) bind(c)

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
  f_prob = probability_funct(x=f_x)

  ! out: f_prob 0D_NOT_real
  call c_f_pointer(prob, f_prob_ptr)
  f_prob_ptr = f_prob
end subroutine
subroutine fortran_quadratic_roots (coefs, root) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: coefs
  real(rp) :: f_coefs(3)
  real(c_double), pointer :: f_coefs_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: root
  complex(rp) :: f_root(2)
  complex(c_double_complex), pointer :: f_root_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(coefs)) then
    call c_f_pointer(coefs, f_coefs_ptr, [3])
    f_coefs = f_coefs_ptr(:)
  endif
  f_root = quadratic_roots(coefs=f_coefs)

  ! out: f_root 1D_NOT_complex
  if (assc(root)) f_root_ptr = f_root
end subroutine
subroutine fortran_query_string (query_str, upcase, return_str, ix, ios) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: query_str
  character(len=4096) :: f_query_str
  character(kind=c_char), pointer :: f_query_str_ptr(:)
  type(c_ptr), intent(in), value :: upcase  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_upcase
  logical :: f_upcase_native
  logical(c_bool), pointer :: f_upcase_ptr
  type(c_ptr), intent(in), value :: return_str
  character(len=4096) :: f_return_str
  character(kind=c_char), pointer :: f_return_str_ptr(:)
  type(c_ptr), intent(in), value :: ix  ! 0D_NOT_integer
  integer(c_int) :: f_ix
  integer(c_int), pointer :: f_ix_ptr
  type(c_ptr), intent(in), value :: ios  ! 0D_NOT_integer
  integer(c_int) :: f_ios
  integer(c_int), pointer :: f_ios_ptr
  ! ** End of parameters **
  ! inout: f_query_str 0D_NOT_character
  if (.not. assc(query_str)) return
  call c_f_pointer(query_str, f_query_str_ptr, [huge(0)])
  call to_f_str(f_query_str_ptr, f_query_str)
  ! inout: f_upcase 0D_NOT_logical
  if (assc(upcase)) then
    call c_f_pointer(upcase, f_upcase_ptr)
    f_upcase_native = f_upcase_ptr
  else
    ! upcase unset
  endif
  ! inout: f_return_str 0D_NOT_character
  if (.not. assc(return_str)) return
  call c_f_pointer(return_str, f_return_str_ptr, [huge(0)])
  call to_f_str(f_return_str_ptr, f_return_str)
  ! inout: f_ix 0D_NOT_integer
  if (assc(ix)) then
    call c_f_pointer(ix, f_ix_ptr)
    f_ix = f_ix_ptr
  else
    ! ix unset
  endif
  ! inout: f_ios 0D_NOT_integer
  if (assc(ios)) then
    call c_f_pointer(ios, f_ios_ptr)
    f_ios = f_ios_ptr
  else
    ! ios unset
  endif
  call query_string(query_str=f_query_str, upcase=f_upcase_native, return_str=f_return_str, &
      ix=f_ix, ios=f_ios)

  ! inout: f_query_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_upcase 0D_NOT_logical
  if (assc(upcase)) then
    call c_f_pointer(upcase, f_upcase_ptr)
    f_upcase_ptr = f_upcase_native
  else
    ! f_upcase unset
  endif
  ! inout: f_return_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ix 0D_NOT_integer
  if (assc(ix)) then
    call c_f_pointer(ix, f_ix_ptr)
    f_ix_ptr = f_ix
  endif
  ! inout: f_ios 0D_NOT_integer
  if (assc(ios)) then
    call c_f_pointer(ios, f_ios_ptr)
    f_ios_ptr = f_ios
  endif
end subroutine
subroutine fortran_quote (str, q_str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_str
  character(len=4096) :: f_q_str
  character(kind=c_char), pointer :: f_q_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. assc(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  f_q_str = quote(str=f_str)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_q_str 0D_ALLOC_character
  call c_f_pointer(q_str, f_q_str_ptr, [len_trim(f_q_str) + 1]) ! output-only string
  call to_c_str(f_q_str, f_q_str_ptr)
end subroutine
subroutine fortran_real_to_string (real_num, width, n_signif, n_decimal, str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: real_num  ! 0D_NOT_real
  real(c_double) :: f_real_num
  real(c_double), pointer :: f_real_num_ptr
  type(c_ptr), intent(in), value :: width  ! 0D_NOT_integer
  integer(c_int) :: f_width
  integer(c_int), pointer :: f_width_ptr
  type(c_ptr), intent(in), value :: n_signif  ! 0D_NOT_integer
  integer(c_int) :: f_n_signif
  integer(c_int), pointer :: f_n_signif_ptr
  type(c_ptr), intent(in), value :: n_decimal  ! 0D_NOT_integer
  integer(c_int) :: f_n_decimal
  integer(c_int), pointer :: f_n_decimal_ptr
  ! ** End of parameters **
  ! inout: f_real_num 0D_NOT_real
  if (assc(real_num)) then
    call c_f_pointer(real_num, f_real_num_ptr)
    f_real_num = f_real_num_ptr
  else
    ! real_num unset
  endif
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width = f_width_ptr
  else
    ! width unset
  endif
  ! inout: f_n_signif 0D_NOT_integer
  if (assc(n_signif)) then
    call c_f_pointer(n_signif, f_n_signif_ptr)
    f_n_signif = f_n_signif_ptr
  else
    ! n_signif unset
  endif
  ! inout: f_n_decimal 0D_NOT_integer
  if (assc(n_decimal)) then
    call c_f_pointer(n_decimal, f_n_decimal_ptr)
    f_n_decimal = f_n_decimal_ptr
  else
    ! n_decimal unset
  endif
  if (assc(n_signif) .and. assc(n_decimal)) then

    f_str = real_to_string(real_num=f_real_num, width=f_width, n_signif=f_n_signif, &
        n_decimal=f_n_decimal)

  elseif (assc(n_signif)) then

    f_str = real_to_string(real_num=f_real_num, width=f_width, n_signif=f_n_signif)

  elseif (assc(n_decimal)) then

    f_str = real_to_string(real_num=f_real_num, width=f_width, n_decimal=f_n_decimal)

  else
    f_str = real_to_string(real_num=f_real_num, width=f_width)

  endif
  ! inout: f_real_num 0D_NOT_real
  if (assc(real_num)) then
    call c_f_pointer(real_num, f_real_num_ptr)
    f_real_num_ptr = f_real_num
  endif
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width_ptr = f_width
  endif
  ! inout: f_n_signif 0D_NOT_integer
  if (assc(n_signif)) then
    call c_f_pointer(n_signif, f_n_signif_ptr)
    f_n_signif_ptr = f_n_signif
  endif
  ! inout: f_n_decimal 0D_NOT_integer
  if (assc(n_decimal)) then
    call c_f_pointer(n_decimal, f_n_decimal_ptr)
    f_n_decimal_ptr = f_n_decimal
  endif
  ! out: f_str 0D_NOT_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1]) ! output-only string
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_real_num_fortran_format (number, width, n_blanks, fmt_str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: fmt_str
  character(len=4096) :: f_fmt_str
  character(kind=c_char), pointer :: f_fmt_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: number  ! 0D_NOT_real
  real(c_double) :: f_number
  real(c_double), pointer :: f_number_ptr
  type(c_ptr), intent(in), value :: width  ! 0D_NOT_integer
  integer(c_int) :: f_width
  integer(c_int), pointer :: f_width_ptr
  type(c_ptr), intent(in), value :: n_blanks  ! 0D_NOT_integer
  integer(c_int) :: f_n_blanks
  integer(c_int), pointer :: f_n_blanks_ptr
  ! ** End of parameters **
  ! inout: f_number 0D_NOT_real
  if (assc(number)) then
    call c_f_pointer(number, f_number_ptr)
    f_number = f_number_ptr
  else
    ! number unset
  endif
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width = f_width_ptr
  else
    ! width unset
  endif
  ! inout: f_n_blanks 0D_NOT_integer
  if (assc(n_blanks)) then
    call c_f_pointer(n_blanks, f_n_blanks_ptr)
    f_n_blanks = f_n_blanks_ptr
  else
    ! n_blanks unset
  endif
  if (assc(n_blanks)) then

    f_fmt_str = real_num_fortran_format(number=f_number, width=f_width, n_blanks=f_n_blanks)

  else
    f_fmt_str = real_num_fortran_format(number=f_number, width=f_width)

  endif
  ! inout: f_number 0D_NOT_real
  if (assc(number)) then
    call c_f_pointer(number, f_number_ptr)
    f_number_ptr = f_number
  endif
  ! inout: f_width 0D_NOT_integer
  if (assc(width)) then
    call c_f_pointer(width, f_width_ptr)
    f_width_ptr = f_width
  endif
  ! inout: f_n_blanks 0D_NOT_integer
  if (assc(n_blanks)) then
    call c_f_pointer(n_blanks, f_n_blanks_ptr)
    f_n_blanks_ptr = f_n_blanks
  endif
  ! out: f_fmt_str 0D_NOT_character
  call c_f_pointer(fmt_str, f_fmt_str_ptr, [len_trim(f_fmt_str) + 1]) ! output-only string
  call to_c_str(f_fmt_str, f_fmt_str_ptr)
end subroutine
subroutine fortran_str_count (str, match, num) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: num  ! 0D_NOT_integer
  integer :: f_num
  integer(c_int), pointer :: f_num_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: match
  character(len=4096) :: f_match
  character(kind=c_char), pointer :: f_match_ptr(:)
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. assc(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! inout: f_match 0D_NOT_character
  if (.not. assc(match)) return
  call c_f_pointer(match, f_match_ptr, [huge(0)])
  call to_f_str(f_match_ptr, f_match)
  f_num = str_count(str=f_str, match=f_match)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_match 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_num 0D_NOT_integer
  call c_f_pointer(num, f_num_ptr)
  f_num_ptr = f_num
end subroutine
subroutine fortran_real_path (path_in, path_out, is_ok) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_ok  ! 0D_NOT_logical
  logical :: f_is_ok
  logical(c_bool), pointer :: f_is_ok_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: path_in
  character(len=4096) :: f_path_in
  character(kind=c_char), pointer :: f_path_in_ptr(:)
  type(c_ptr), intent(in), value :: path_out
  character(len=4096) :: f_path_out
  character(kind=c_char), pointer :: f_path_out_ptr(:)
  ! ** End of parameters **
  ! inout: f_path_in 0D_NOT_character
  if (.not. assc(path_in)) return
  call c_f_pointer(path_in, f_path_in_ptr, [huge(0)])
  call to_f_str(f_path_in_ptr, f_path_in)
  ! inout: f_path_out 0D_NOT_character
  if (.not. assc(path_out)) return
  call c_f_pointer(path_out, f_path_out_ptr, [huge(0)])
  call to_f_str(f_path_out_ptr, f_path_out)
  f_is_ok = real_path(path_in=f_path_in, path_out=f_path_out)

  ! inout: f_path_in 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_path_out 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_is_ok 0D_NOT_logical
  call c_f_pointer(is_ok, f_is_ok_ptr)
  f_is_ok_ptr = f_is_ok
end subroutine
subroutine fortran_real_str (r_num, n_signif, n_decimal, str) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: r_num  ! 0D_NOT_real
  real(c_double) :: f_r_num
  real(c_double), pointer :: f_r_num_ptr
  type(c_ptr), intent(in), value :: n_signif  ! 0D_NOT_integer
  integer(c_int) :: f_n_signif
  integer(c_int), pointer :: f_n_signif_ptr
  type(c_ptr), intent(in), value :: n_decimal  ! 0D_NOT_integer
  integer(c_int) :: f_n_decimal
  integer(c_int), pointer :: f_n_decimal_ptr
  ! ** End of parameters **
  ! inout: f_r_num 0D_NOT_real
  if (assc(r_num)) then
    call c_f_pointer(r_num, f_r_num_ptr)
    f_r_num = f_r_num_ptr
  else
    ! r_num unset
  endif
  ! inout: f_n_signif 0D_NOT_integer
  if (assc(n_signif)) then
    call c_f_pointer(n_signif, f_n_signif_ptr)
    f_n_signif = f_n_signif_ptr
  else
    ! n_signif unset
  endif
  ! inout: f_n_decimal 0D_NOT_integer
  if (assc(n_decimal)) then
    call c_f_pointer(n_decimal, f_n_decimal_ptr)
    f_n_decimal = f_n_decimal_ptr
  else
    ! n_decimal unset
  endif
  if (assc(n_signif) .and. assc(n_decimal)) then

    f_str = real_str(r_num=f_r_num, n_signif=f_n_signif, n_decimal=f_n_decimal)

  elseif (assc(n_signif)) then

    f_str = real_str(r_num=f_r_num, n_signif=f_n_signif)

  elseif (assc(n_decimal)) then

    f_str = real_str(r_num=f_r_num, n_decimal=f_n_decimal)

  else
    f_str = real_str(r_num=f_r_num)

  endif
  ! inout: f_r_num 0D_NOT_real
  if (assc(r_num)) then
    call c_f_pointer(r_num, f_r_num_ptr)
    f_r_num_ptr = f_r_num
  endif
  ! inout: f_n_signif 0D_NOT_integer
  if (assc(n_signif)) then
    call c_f_pointer(n_signif, f_n_signif_ptr)
    f_n_signif_ptr = f_n_signif
  endif
  ! inout: f_n_decimal 0D_NOT_integer
  if (assc(n_decimal)) then
    call c_f_pointer(n_decimal, f_n_decimal_ptr)
    f_n_decimal_ptr = f_n_decimal
  endif
  ! out: f_str 0D_ALLOC_character
  call c_f_pointer(str, f_str_ptr, [len_trim(f_str) + 1]) ! output-only string
  call to_c_str(f_str, f_str_ptr)
end subroutine
subroutine fortran_rms_value (val_arr, good_val, ave_val, rms_val) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: val_arr
  type(real_container_alloc), pointer :: f_val_arr
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
  !! container general array (1D_ALLOC_real)
  if (assc(val_arr))   call c_f_pointer(val_arr, f_val_arr)
  !! container general array (1D_ALLOC_logical)
  if (assc(good_val))   call c_f_pointer(good_val, f_good_val)
  if (assc(good_val) .and. assc(ave_val)) then

    f_rms_val = rms_value(val_arr=f_val_arr%data, good_val=f_good_val%data, ave_val=f_ave_val)

  elseif (assc(good_val)) then

    f_rms_val = rms_value(val_arr=f_val_arr%data, good_val=f_good_val%data)

  elseif (assc(ave_val)) then

    f_rms_val = rms_value(val_arr=f_val_arr%data, ave_val=f_ave_val)

  else
    f_rms_val = rms_value(val_arr=f_val_arr%data)

  endif
  ! out: f_ave_val 0D_NOT_real
  if (assc(ave_val)) then
    call c_f_pointer(ave_val, f_ave_val_ptr)
    f_ave_val_ptr = f_ave_val
  endif
  ! out: f_rms_val 0D_NOT_real
  call c_f_pointer(rms_val, f_rms_val_ptr)
  f_rms_val_ptr = f_rms_val
end subroutine
subroutine fortran_rot_2d (vec_in, angle, vec_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: vec_in
  real(rp) :: f_vec_in(2)
  real(c_double), pointer :: f_vec_in_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: vec_out
  real(rp) :: f_vec_out(2)
  real(c_double), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(vec_in)) then
    call c_f_pointer(vec_in, f_vec_in_ptr, [2])
    f_vec_in = f_vec_in_ptr(:)
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  f_vec_out = rot_2d(vec_in=f_vec_in, angle=f_angle)

  ! out: f_vec_out 1D_NOT_real
  if (assc(vec_out)) f_vec_out_ptr = f_vec_out
end subroutine
subroutine fortran_run_timer (command, time, time0) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: command
  character(len=4096) :: f_command
  character(kind=c_char), pointer :: f_command_ptr(:)
  type(c_ptr), intent(in), value :: time  ! 0D_NOT_real
  real(c_double) :: f_time
  real(c_double), pointer :: f_time_ptr
  type(c_ptr), intent(in), value :: time0  ! 0D_NOT_real
  real(c_double) :: f_time0
  real(c_double), pointer :: f_time0_ptr
  ! ** End of parameters **
  ! inout: f_command 0D_NOT_character
  if (.not. assc(command)) return
  call c_f_pointer(command, f_command_ptr, [huge(0)])
  call to_f_str(f_command_ptr, f_command)
  ! inout: f_time 0D_NOT_real
  if (assc(time)) then
    call c_f_pointer(time, f_time_ptr)
    f_time = f_time_ptr
  else
    ! time unset
  endif
  ! inout: f_time0 0D_NOT_real
  if (assc(time0)) then
    call c_f_pointer(time0, f_time0_ptr)
    f_time0 = f_time0_ptr
  else
    ! time0 unset
  endif
  if (assc(time) .and. assc(time0)) then

    call run_timer(command=f_command, time=f_time, time0=f_time0)

  elseif (assc(time)) then

    call run_timer(command=f_command, time=f_time)

  elseif (assc(time0)) then

    call run_timer(command=f_command, time0=f_time0)

  else
    call run_timer(command=f_command)

  endif
  ! inout: f_command 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_time 0D_NOT_real
  if (assc(time)) then
    call c_f_pointer(time, f_time_ptr)
    f_time_ptr = f_time
  endif
  ! inout: f_time0 0D_NOT_real
  if (assc(time0)) then
    call c_f_pointer(time0, f_time0_ptr)
    f_time0_ptr = f_time0
  endif
end subroutine
subroutine fortran_sinc (x, nd, y) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** End of parameters **
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd = f_nd_ptr
  else
    ! nd unset
  endif
  if (assc(nd)) then

    f_y = sinc(x=f_x, nd=f_nd)

  else
    f_y = sinc(x=f_x)

  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd_ptr = f_nd
  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_sincc (x, nd, y) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** End of parameters **
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd = f_nd_ptr
  else
    ! nd unset
  endif
  if (assc(nd)) then

    f_y = sincc(x=f_x, nd=f_nd)

  else
    f_y = sincc(x=f_x)

  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd_ptr = f_nd
  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_sinhx_x (x, nd, y) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** End of parameters **
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd = f_nd_ptr
  else
    ! nd unset
  endif
  if (assc(nd)) then

    f_y = sinhx_x(x=f_x, nd=f_nd)

  else
    f_y = sinhx_x(x=f_x)

  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd_ptr = f_nd
  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_skip_header (ix_unit, error_flag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_unit  ! 0D_NOT_integer
  integer(c_int) :: f_ix_unit
  integer(c_int), pointer :: f_ix_unit_ptr
  type(c_ptr), intent(in), value :: error_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_error_flag
  logical :: f_error_flag_native
  logical(c_bool), pointer :: f_error_flag_ptr
  ! ** End of parameters **
  ! inout: f_ix_unit 0D_NOT_integer
  if (assc(ix_unit)) then
    call c_f_pointer(ix_unit, f_ix_unit_ptr)
    f_ix_unit = f_ix_unit_ptr
  else
    ! ix_unit unset
  endif
  ! inout: f_error_flag 0D_NOT_logical
  if (assc(error_flag)) then
    call c_f_pointer(error_flag, f_error_flag_ptr)
    f_error_flag_native = f_error_flag_ptr
  else
    ! error_flag unset
  endif
  call skip_header(ix_unit=f_ix_unit, error_flag=f_error_flag_native)

  ! inout: f_ix_unit 0D_NOT_integer
  if (assc(ix_unit)) then
    call c_f_pointer(ix_unit, f_ix_unit_ptr)
    f_ix_unit_ptr = f_ix_unit
  endif
  ! inout: f_error_flag 0D_NOT_logical
  if (assc(error_flag)) then
    call c_f_pointer(error_flag, f_error_flag_ptr)
    f_error_flag_ptr = f_error_flag_native
  else
    ! f_error_flag unset
  endif
end subroutine
subroutine fortran_sqrt_one (x, nd, ds1) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ds1  ! 0D_NOT_real
  real(rp) :: f_ds1
  real(c_double), pointer :: f_ds1_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  type(c_ptr), intent(in), value :: nd  ! 0D_NOT_integer
  integer(c_int) :: f_nd
  integer(c_int), pointer :: f_nd_ptr
  ! ** End of parameters **
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd = f_nd_ptr
  else
    ! nd unset
  endif
  if (assc(nd)) then

    f_ds1 = sqrt_one(x=f_x, nd=f_nd)

  else
    f_ds1 = sqrt_one(x=f_x)

  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! inout: f_nd 0D_NOT_integer
  if (assc(nd)) then
    call c_f_pointer(nd, f_nd_ptr)
    f_nd_ptr = f_nd
  endif
  ! out: f_ds1 0D_NOT_real
  call c_f_pointer(ds1, f_ds1_ptr)
  f_ds1_ptr = f_ds1
end subroutine
subroutine fortran_sqrt_alpha (alpha, x, y) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: alpha  ! 0D_NOT_real
  real(c_double) :: f_alpha
  real(c_double), pointer :: f_alpha_ptr
  type(c_ptr), intent(in), value :: x  ! 0D_NOT_real
  real(c_double) :: f_x
  real(c_double), pointer :: f_x_ptr
  ! ** End of parameters **
  ! inout: f_alpha 0D_NOT_real
  if (assc(alpha)) then
    call c_f_pointer(alpha, f_alpha_ptr)
    f_alpha = f_alpha_ptr
  else
    ! alpha unset
  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x = f_x_ptr
  else
    ! x unset
  endif
  f_y = sqrt_alpha(alpha=f_alpha, x=f_x)

  ! inout: f_alpha 0D_NOT_real
  if (assc(alpha)) then
    call c_f_pointer(alpha, f_alpha_ptr)
    f_alpha_ptr = f_alpha
  endif
  ! inout: f_x 0D_NOT_real
  if (assc(x)) then
    call c_f_pointer(x, f_x_ptr)
    f_x_ptr = f_x
  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_str_first_in_set (line, set, ignore_clauses, ix_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096) :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  type(c_ptr), intent(in), value :: ignore_clauses  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_clauses
  logical :: f_ignore_clauses_native
  logical(c_bool), pointer :: f_ignore_clauses_ptr
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_set 0D_NOT_character
  if (.not. assc(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  ! inout: f_ignore_clauses 0D_NOT_logical
  if (assc(ignore_clauses)) then
    call c_f_pointer(ignore_clauses, f_ignore_clauses_ptr)
    f_ignore_clauses_native = f_ignore_clauses_ptr
  else
    ! ignore_clauses unset
  endif
  if (assc(ignore_clauses)) then

    f_ix_match = str_first_in_set(line=f_line, set=f_set, &
        ignore_clauses=f_ignore_clauses_native)

  else
    f_ix_match = str_first_in_set(line=f_line, set=f_set)

  endif
  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_set 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ignore_clauses 0D_NOT_logical
  if (assc(ignore_clauses)) then
    call c_f_pointer(ignore_clauses, f_ignore_clauses_ptr)
    f_ignore_clauses_ptr = f_ignore_clauses_native
  else
    ! f_ignore_clauses unset
  endif
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_first_not_in_set (line, set, ix_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096) :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_set 0D_NOT_character
  if (.not. assc(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  f_ix_match = str_first_not_in_set(line=f_line, set=f_set)

  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_set 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_last_in_set (line, set, ix_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096) :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_set 0D_NOT_character
  if (.not. assc(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  f_ix_match = str_last_in_set(line=f_line, set=f_set)

  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_set 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_str_last_not_in_set (line, set, ix_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: set
  character(len=4096) :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_set 0D_NOT_character
  if (.not. assc(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  f_ix_match = str_last_not_in_set(line=f_line, set=f_set)

  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_set 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_string_to_int (line, default_, err_flag, err_print_flag, value) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_integer
  integer :: f_value
  integer(c_int), pointer :: f_value_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: default_  ! 0D_NOT_integer
  integer(c_int) :: f_default
  integer(c_int), pointer :: f_default_ptr
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical :: f_err_flag_native
  logical(c_bool), pointer :: f_err_flag_ptr
  type(c_ptr), intent(in), value :: err_print_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_print_flag
  logical :: f_err_print_flag_native
  logical(c_bool), pointer :: f_err_print_flag_ptr
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_default 0D_NOT_integer
  if (assc(default_)) then
    call c_f_pointer(default_, f_default_ptr)
    f_default = f_default_ptr
  else
    ! default_ unset
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
  else
    ! err_flag unset
  endif
  ! inout: f_err_print_flag 0D_NOT_logical
  if (assc(err_print_flag)) then
    call c_f_pointer(err_print_flag, f_err_print_flag_ptr)
    f_err_print_flag_native = f_err_print_flag_ptr
  else
    ! err_print_flag unset
  endif
  if (assc(err_print_flag)) then

    f_value = string_to_int(line=f_line, default=f_default, err_flag=f_err_flag_native, &
        err_print_flag=f_err_print_flag_native)

  else
    f_value = string_to_int(line=f_line, default=f_default, err_flag=f_err_flag_native)

  endif
  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_default 0D_NOT_integer
  if (assc(default_)) then
    call c_f_pointer(default_, f_default_ptr)
    f_default_ptr = f_default
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
  ! inout: f_err_print_flag 0D_NOT_logical
  if (assc(err_print_flag)) then
    call c_f_pointer(err_print_flag, f_err_print_flag_ptr)
    f_err_print_flag_ptr = f_err_print_flag_native
  else
    ! f_err_print_flag unset
  endif
  ! out: f_value 0D_NOT_integer
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_string_to_real (line, default_, err_flag, err_print_flag, value) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: default_  ! 0D_NOT_real
  real(c_double) :: f_default
  real(c_double), pointer :: f_default_ptr
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical :: f_err_flag_native
  logical(c_bool), pointer :: f_err_flag_ptr
  type(c_ptr), intent(in), value :: err_print_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_print_flag
  logical :: f_err_print_flag_native
  logical(c_bool), pointer :: f_err_print_flag_ptr
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_default 0D_NOT_real
  if (assc(default_)) then
    call c_f_pointer(default_, f_default_ptr)
    f_default = f_default_ptr
  else
    ! default_ unset
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
  else
    ! err_flag unset
  endif
  ! inout: f_err_print_flag 0D_NOT_logical
  if (assc(err_print_flag)) then
    call c_f_pointer(err_print_flag, f_err_print_flag_ptr)
    f_err_print_flag_native = f_err_print_flag_ptr
  else
    ! err_print_flag unset
  endif
  if (assc(err_print_flag)) then

    f_value = string_to_real(line=f_line, default=f_default, err_flag=f_err_flag_native, &
        err_print_flag=f_err_print_flag_native)

  else
    f_value = string_to_real(line=f_line, default=f_default, err_flag=f_err_flag_native)

  endif
  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_default 0D_NOT_real
  if (assc(default_)) then
    call c_f_pointer(default_, f_default_ptr)
    f_default_ptr = f_default
  endif
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
  ! inout: f_err_print_flag 0D_NOT_logical
  if (assc(err_print_flag)) then
    call c_f_pointer(err_print_flag, f_err_print_flag_ptr)
    f_err_print_flag_ptr = f_err_print_flag_native
  else
    ! f_err_print_flag unset
  endif
  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096) :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  type(c_ptr), intent(in), value :: delimitors
  character(len=4096) :: f_delimitors
  character(kind=c_char), pointer :: f_delimitors_ptr(:)
  type(c_ptr), intent(in), value :: out_str
  character(len=4096) :: f_out_str
  character(kind=c_char), pointer :: f_out_str_ptr(:)
  type(c_ptr), intent(in), value :: ix_word  ! 0D_NOT_integer
  integer(c_int) :: f_ix_word
  integer(c_int), pointer :: f_ix_word_ptr
  type(c_ptr), intent(in), value :: delim
  character(len=4096) :: f_delim
  character(kind=c_char), pointer :: f_delim_ptr(:)
  type(c_ptr), intent(in), value :: ix_next  ! 0D_NOT_integer
  integer(c_int) :: f_ix_next
  integer(c_int), pointer :: f_ix_next_ptr
  ! ** End of parameters **
  ! inout: f_in_str 0D_NOT_character
  if (.not. assc(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  ! inout: f_delimitors 0D_NOT_character
  if (.not. assc(delimitors)) return
  call c_f_pointer(delimitors, f_delimitors_ptr, [huge(0)])
  call to_f_str(f_delimitors_ptr, f_delimitors)
  ! inout: f_out_str 0D_NOT_character
  if (.not. assc(out_str)) return
  call c_f_pointer(out_str, f_out_str_ptr, [huge(0)])
  call to_f_str(f_out_str_ptr, f_out_str)
  ! inout: f_ix_word 0D_NOT_integer
  if (assc(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
    f_ix_word = f_ix_word_ptr
  else
    ! ix_word unset
  endif
  ! inout: f_delim 0D_NOT_character
  if (.not. assc(delim)) return
  call c_f_pointer(delim, f_delim_ptr, [huge(0)])
  call to_f_str(f_delim_ptr, f_delim)
  ! inout: f_ix_next 0D_NOT_integer
  if (assc(ix_next)) then
    call c_f_pointer(ix_next, f_ix_next_ptr)
    f_ix_next = f_ix_next_ptr
  else
    ! ix_next unset
  endif
  call string_trim2(in_str=f_in_str, delimitors=f_delimitors, out_str=f_out_str, &
      ix_word=f_ix_word, delim=f_delim, ix_next=f_ix_next)

  ! inout: f_in_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_delimitors 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_out_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ix_word 0D_NOT_integer
  if (assc(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
    f_ix_word_ptr = f_ix_word
  endif
  ! inout: f_delim 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ix_next 0D_NOT_integer
  if (assc(ix_next)) then
    call c_f_pointer(ix_next, f_ix_next_ptr)
    f_ix_next_ptr = f_ix_next
  endif
end subroutine
subroutine fortran_to_str (num, max_signif, string) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: num  ! 0D_NOT_real
  real(c_double) :: f_num
  real(c_double), pointer :: f_num_ptr
  type(c_ptr), intent(in), value :: max_signif  ! 0D_NOT_integer
  integer(c_int) :: f_max_signif
  integer(c_int), pointer :: f_max_signif_ptr
  ! ** End of parameters **
  ! inout: f_num 0D_NOT_real
  if (assc(num)) then
    call c_f_pointer(num, f_num_ptr)
    f_num = f_num_ptr
  else
    ! num unset
  endif
  ! inout: f_max_signif 0D_NOT_integer
  if (assc(max_signif)) then
    call c_f_pointer(max_signif, f_max_signif_ptr)
    f_max_signif = f_max_signif_ptr
  else
    ! max_signif unset
  endif
  if (assc(max_signif)) then

    f_string = to_str(num=f_num, max_signif=f_max_signif)

  else
    f_string = to_str(num=f_num)

  endif
  ! inout: f_num 0D_NOT_real
  if (assc(num)) then
    call c_f_pointer(num, f_num_ptr)
    f_num_ptr = f_num
  endif
  ! inout: f_max_signif 0D_NOT_integer
  if (assc(max_signif)) then
    call c_f_pointer(max_signif, f_max_signif_ptr)
    f_max_signif_ptr = f_max_signif
  endif
  ! out: f_string 0D_ALLOC_character
  call c_f_pointer(string, f_string_ptr, [len_trim(f_string) + 1]) ! output-only string
  call to_c_str(f_string, f_string_ptr)
end subroutine
subroutine fortran_type_this_file (filename) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: filename
  character(len=4096) :: f_filename
  character(kind=c_char), pointer :: f_filename_ptr(:)
  ! ** End of parameters **
  ! inout: f_filename 0D_NOT_character
  if (.not. assc(filename)) return
  call c_f_pointer(filename, f_filename_ptr, [huge(0)])
  call to_f_str(f_filename_ptr, f_filename)
  call type_this_file(filename=f_filename)

  ! inout: f_filename 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_upcase_string (string) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  call upcase_string(string=f_string)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_word_len (wording, wlen) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: wlen  ! 0D_NOT_integer
  integer :: f_wlen
  integer(c_int), pointer :: f_wlen_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: wording
  character(len=4096) :: f_wording
  character(kind=c_char), pointer :: f_wording_ptr(:)
  ! ** End of parameters **
  ! inout: f_wording 0D_NOT_character
  if (.not. assc(wording)) return
  call c_f_pointer(wording, f_wording_ptr, [huge(0)])
  call to_f_str(f_wording_ptr, f_wording)
  f_wlen = word_len(wording=f_wording)

  ! inout: f_wording 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_wlen 0D_NOT_integer
  call c_f_pointer(wlen, f_wlen_ptr)
  f_wlen_ptr = f_wlen
end subroutine
subroutine fortran_word_read (in_str, delim_list, word, ix_word, delim, delim_found, out_str, &
    ignore_interior) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096) :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  type(c_ptr), intent(in), value :: delim_list
  character(len=4096) :: f_delim_list
  character(kind=c_char), pointer :: f_delim_list_ptr(:)
  type(c_ptr), intent(in), value :: word
  character(len=4096) :: f_word
  character(kind=c_char), pointer :: f_word_ptr(:)
  type(c_ptr), intent(in), value :: ix_word  ! 0D_NOT_integer
  integer(c_int) :: f_ix_word
  integer(c_int), pointer :: f_ix_word_ptr
  type(c_ptr), intent(in), value :: delim
  character(len=4096) :: f_delim
  character(kind=c_char), pointer :: f_delim_ptr(:)
  type(c_ptr), intent(in), value :: delim_found  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_delim_found
  logical :: f_delim_found_native
  logical(c_bool), pointer :: f_delim_found_ptr
  type(c_ptr), intent(in), value :: out_str
  character(len=4096) :: f_out_str
  character(kind=c_char), pointer :: f_out_str_ptr(:)
  type(c_ptr), intent(in), value :: ignore_interior  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_interior
  logical :: f_ignore_interior_native
  logical(c_bool), pointer :: f_ignore_interior_ptr
  ! ** End of parameters **
  ! inout: f_in_str 0D_NOT_character
  if (.not. assc(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  ! inout: f_delim_list 0D_NOT_character
  if (.not. assc(delim_list)) return
  call c_f_pointer(delim_list, f_delim_list_ptr, [huge(0)])
  call to_f_str(f_delim_list_ptr, f_delim_list)
  ! inout: f_word 0D_NOT_character
  if (.not. assc(word)) return
  call c_f_pointer(word, f_word_ptr, [huge(0)])
  call to_f_str(f_word_ptr, f_word)
  ! inout: f_ix_word 0D_NOT_integer
  if (assc(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
    f_ix_word = f_ix_word_ptr
  else
    ! ix_word unset
  endif
  ! inout: f_delim 0D_NOT_character
  if (.not. assc(delim)) return
  call c_f_pointer(delim, f_delim_ptr, [huge(0)])
  call to_f_str(f_delim_ptr, f_delim)
  ! inout: f_delim_found 0D_NOT_logical
  if (assc(delim_found)) then
    call c_f_pointer(delim_found, f_delim_found_ptr)
    f_delim_found_native = f_delim_found_ptr
  else
    ! delim_found unset
  endif
  ! inout: f_out_str 0D_NOT_character
  if (.not. assc(out_str)) return
  call c_f_pointer(out_str, f_out_str_ptr, [huge(0)])
  call to_f_str(f_out_str_ptr, f_out_str)
  ! inout: f_ignore_interior 0D_NOT_logical
  if (assc(ignore_interior)) then
    call c_f_pointer(ignore_interior, f_ignore_interior_ptr)
    f_ignore_interior_native = f_ignore_interior_ptr
  else
    ! ignore_interior unset
  endif
  if (assc(ignore_interior)) then

    call word_read(in_str=f_in_str, delim_list=f_delim_list, word=f_word, ix_word=f_ix_word, &
        delim=f_delim, delim_found=f_delim_found_native, out_str=f_out_str, &
        ignore_interior=f_ignore_interior_native)

  else
    call word_read(in_str=f_in_str, delim_list=f_delim_list, word=f_word, ix_word=f_ix_word, &
        delim=f_delim, delim_found=f_delim_found_native, out_str=f_out_str)

  endif
  ! inout: f_in_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_delim_list 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_word 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ix_word 0D_NOT_integer
  if (assc(ix_word)) then
    call c_f_pointer(ix_word, f_ix_word_ptr)
    f_ix_word_ptr = f_ix_word
  endif
  ! inout: f_delim 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_delim_found 0D_NOT_logical
  if (assc(delim_found)) then
    call c_f_pointer(delim_found, f_delim_found_ptr)
    f_delim_found_ptr = f_delim_found_native
  else
    ! f_delim_found unset
  endif
  ! inout: f_out_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_ignore_interior 0D_NOT_logical
  if (assc(ignore_interior)) then
    call c_f_pointer(ignore_interior, f_ignore_interior_ptr)
    f_ignore_interior_ptr = f_ignore_interior_native
  else
    ! f_ignore_interior unset
  endif
end subroutine
subroutine fortran_str_substitute (string, str_match, str_replace, do_trim, ignore_escaped) &
    bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096) :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: str_match
  character(len=4096) :: f_str_match
  character(kind=c_char), pointer :: f_str_match_ptr(:)
  type(c_ptr), intent(in), value :: str_replace
  character(len=4096) :: f_str_replace
  character(kind=c_char), pointer :: f_str_replace_ptr(:)
  type(c_ptr), intent(in), value :: do_trim  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_do_trim
  logical :: f_do_trim_native
  logical(c_bool), pointer :: f_do_trim_ptr
  type(c_ptr), intent(in), value :: ignore_escaped  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_escaped
  logical :: f_ignore_escaped_native
  logical(c_bool), pointer :: f_ignore_escaped_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. assc(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_str_match 0D_NOT_character
  if (assc(str_match)) then
    call c_f_pointer(str_match, f_str_match_ptr, [huge(0)])
    call to_f_str(f_str_match_ptr, f_str_match)
  endif
  ! inout: f_str_replace 0D_NOT_character
  if (assc(str_replace)) then
    call c_f_pointer(str_replace, f_str_replace_ptr, [huge(0)])
    call to_f_str(f_str_replace_ptr, f_str_replace)
  endif
  ! inout: f_do_trim 0D_NOT_logical
  if (assc(do_trim)) then
    call c_f_pointer(do_trim, f_do_trim_ptr)
    f_do_trim_native = f_do_trim_ptr
  else
    ! do_trim unset
  endif
  ! inout: f_ignore_escaped 0D_NOT_logical
  if (assc(ignore_escaped)) then
    call c_f_pointer(ignore_escaped, f_ignore_escaped_ptr)
    f_ignore_escaped_native = f_ignore_escaped_ptr
  else
    ! ignore_escaped unset
  endif
  if (assc(str_match) .and. assc(str_replace) .and. assc(do_trim) .and. assc(ignore_escaped)) &
      then

    call str_substitute(string=f_string, str_match=f_str_match, str_replace=f_str_replace, &
        do_trim=f_do_trim_native, ignore_escaped=f_ignore_escaped_native)

  elseif (assc(str_match) .and. assc(str_replace) .and. assc(do_trim)) then

    call str_substitute(string=f_string, str_match=f_str_match, str_replace=f_str_replace, &
        do_trim=f_do_trim_native)

  elseif (assc(str_match) .and. assc(str_replace) .and. assc(ignore_escaped)) then

    call str_substitute(string=f_string, str_match=f_str_match, str_replace=f_str_replace, &
        ignore_escaped=f_ignore_escaped_native)

  elseif (assc(str_match) .and. assc(do_trim) .and. assc(ignore_escaped)) then

    call str_substitute(string=f_string, str_match=f_str_match, do_trim=f_do_trim_native, &
        ignore_escaped=f_ignore_escaped_native)

  elseif (assc(str_replace) .and. assc(do_trim) .and. assc(ignore_escaped)) then

    call str_substitute(string=f_string, str_replace=f_str_replace, do_trim=f_do_trim_native, &
        ignore_escaped=f_ignore_escaped_native)

  elseif (assc(str_match) .and. assc(str_replace)) then

    call str_substitute(string=f_string, str_match=f_str_match, str_replace=f_str_replace)

  elseif (assc(str_match) .and. assc(do_trim)) then

    call str_substitute(string=f_string, str_match=f_str_match, do_trim=f_do_trim_native)

  elseif (assc(str_match) .and. assc(ignore_escaped)) then

    call str_substitute(string=f_string, str_match=f_str_match, &
        ignore_escaped=f_ignore_escaped_native)

  elseif (assc(str_replace) .and. assc(do_trim)) then

    call str_substitute(string=f_string, str_replace=f_str_replace, do_trim=f_do_trim_native)

  elseif (assc(str_replace) .and. assc(ignore_escaped)) then

    call str_substitute(string=f_string, str_replace=f_str_replace, &
        ignore_escaped=f_ignore_escaped_native)

  elseif (assc(do_trim) .and. assc(ignore_escaped)) then

    call str_substitute(string=f_string, do_trim=f_do_trim_native, &
        ignore_escaped=f_ignore_escaped_native)

  elseif (assc(str_match)) then

    call str_substitute(string=f_string, str_match=f_str_match)

  elseif (assc(str_replace)) then

    call str_substitute(string=f_string, str_replace=f_str_replace)

  elseif (assc(do_trim)) then

    call str_substitute(string=f_string, do_trim=f_do_trim_native)

  elseif (assc(ignore_escaped)) then

    call str_substitute(string=f_string, ignore_escaped=f_ignore_escaped_native)

  else
    call str_substitute(string=f_string)

  endif
  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_str_match 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_str_replace 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_do_trim 0D_NOT_logical
  if (assc(do_trim)) then
    call c_f_pointer(do_trim, f_do_trim_ptr)
    f_do_trim_ptr = f_do_trim_native
  else
    ! f_do_trim unset
  endif
  ! inout: f_ignore_escaped 0D_NOT_logical
  if (assc(ignore_escaped)) then
    call c_f_pointer(ignore_escaped, f_ignore_escaped_ptr)
    f_ignore_escaped_ptr = f_ignore_escaped_native
  else
    ! f_ignore_escaped unset
  endif
end subroutine
subroutine fortran_str_match_wild (str, pat, a_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: a_match  ! 0D_NOT_logical
  logical :: f_a_match
  logical(c_bool), pointer :: f_a_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096) :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: pat
  character(len=4096) :: f_pat
  character(kind=c_char), pointer :: f_pat_ptr(:)
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. assc(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! inout: f_pat 0D_NOT_character
  if (.not. assc(pat)) return
  call c_f_pointer(pat, f_pat_ptr, [huge(0)])
  call to_f_str(f_pat_ptr, f_pat)
  f_a_match = str_match_wild(str=f_str, pat=f_pat)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_pat 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_a_match 0D_NOT_logical
  call c_f_pointer(a_match, f_a_match_ptr)
  f_a_match_ptr = f_a_match
end subroutine
subroutine fortran_str_upcase (dst, src) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: dst
  character(len=4096) :: f_dst
  character(kind=c_char), pointer :: f_dst_ptr(:)
  type(c_ptr), intent(in), value :: src
  character(len=4096) :: f_src
  character(kind=c_char), pointer :: f_src_ptr(:)
  ! ** End of parameters **
  ! inout: f_dst 0D_NOT_character
  if (.not. assc(dst)) return
  call c_f_pointer(dst, f_dst_ptr, [huge(0)])
  call to_f_str(f_dst_ptr, f_dst)
  ! inout: f_src 0D_NOT_character
  if (.not. assc(src)) return
  call c_f_pointer(src, f_src_ptr, [huge(0)])
  call to_f_str(f_src_ptr, f_src)
  call str_upcase(dst=f_dst, src=f_src)

  ! inout: f_dst 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_src 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_str_downcase (dst, src) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: dst
  character(len=4096) :: f_dst
  character(kind=c_char), pointer :: f_dst_ptr(:)
  type(c_ptr), intent(in), value :: src
  character(len=4096) :: f_src
  character(kind=c_char), pointer :: f_src_ptr(:)
  ! ** End of parameters **
  ! inout: f_dst 0D_NOT_character
  if (.not. assc(dst)) return
  call c_f_pointer(dst, f_dst_ptr, [huge(0)])
  call to_f_str(f_dst_ptr, f_dst)
  ! inout: f_src 0D_NOT_character
  if (.not. assc(src)) return
  call c_f_pointer(src, f_src_ptr, [huge(0)])
  call to_f_str(f_src_ptr, f_src)
  call str_downcase(dst=f_dst, src=f_src)

  ! inout: f_dst 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_src 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_system_command (line, err_flag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096) :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical :: f_err_flag_native
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. assc(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
  else
    ! err_flag unset
  endif
  if (assc(err_flag)) then

    call system_command(line=f_line, err_flag=f_err_flag_native)

  else
    call system_command(line=f_line)

  endif
  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_err_flag 0D_NOT_logical
  if (assc(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
end subroutine
subroutine fortran_string_trim (in_string, out_string, word_len) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: in_string
  character(len=4096) :: f_in_string
  character(kind=c_char), pointer :: f_in_string_ptr(:)
  type(c_ptr), intent(in), value :: out_string
  character(len=4096) :: f_out_string
  character(kind=c_char), pointer :: f_out_string_ptr(:)
  type(c_ptr), intent(in), value :: word_len  ! 0D_NOT_integer
  integer(c_int) :: f_word_len
  integer(c_int), pointer :: f_word_len_ptr
  ! ** End of parameters **
  ! inout: f_in_string 0D_NOT_character
  if (.not. assc(in_string)) return
  call c_f_pointer(in_string, f_in_string_ptr, [huge(0)])
  call to_f_str(f_in_string_ptr, f_in_string)
  ! inout: f_out_string 0D_NOT_character
  if (.not. assc(out_string)) return
  call c_f_pointer(out_string, f_out_string_ptr, [huge(0)])
  call to_f_str(f_out_string_ptr, f_out_string)
  ! inout: f_word_len 0D_NOT_integer
  if (assc(word_len)) then
    call c_f_pointer(word_len, f_word_len_ptr)
    f_word_len = f_word_len_ptr
  else
    ! word_len unset
  endif
  call string_trim(in_string=f_in_string, out_string=f_out_string, word_len=f_word_len)

  ! inout: f_in_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_out_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_word_len 0D_NOT_integer
  if (assc(word_len)) then
    call c_f_pointer(word_len, f_word_len_ptr)
    f_word_len_ptr = f_word_len
  endif
end subroutine
subroutine fortran_virtual_memory_usage (usage) bind(c)

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
subroutine fortran_find_location_real (arr, value, ix_match) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr
  type(real_container_alloc), pointer :: f_arr
  real(c_double) :: value  ! 0D_NOT_real
  real(rp) :: f_value
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(arr))   call c_f_pointer(arr, f_arr)
  ! in: f_value 0D_NOT_real
  f_value = value
  f_ix_match = find_location_real(arr=f_arr%data, value=f_value)

  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_find_location_int (arr, value, ix_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr
  type(integer_container_alloc), pointer :: f_arr
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_integer
  integer(c_int) :: f_value
  integer(c_int), pointer :: f_value_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_integer)
  if (assc(arr))   call c_f_pointer(arr, f_arr)
  ! inout: f_value 0D_NOT_integer
  if (assc(value)) then
    call c_f_pointer(value, f_value_ptr)
    f_value = f_value_ptr
  else
    ! value unset
  endif
  f_ix_match = find_location_int(arr=f_arr%data, value=f_value)

  ! inout: f_value 0D_NOT_integer
  if (assc(value)) then
    call c_f_pointer(value, f_value_ptr)
    f_value_ptr = f_value
  endif
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_find_location_logic (arr, value, ix_match) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_match  ! 0D_NOT_integer
  integer :: f_ix_match
  integer(c_int), pointer :: f_ix_match_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr
  type(logical_container_alloc), pointer :: f_arr
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_value
  logical :: f_value_native
  logical(c_bool), pointer :: f_value_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_logical)
  if (assc(arr))   call c_f_pointer(arr, f_arr)
  ! inout: f_value 0D_NOT_logical
  if (assc(value)) then
    call c_f_pointer(value, f_value_ptr)
    f_value_native = f_value_ptr
  else
    ! value unset
  endif
  f_ix_match = find_location_logic(arr=f_arr%data, value=f_value_native)

  ! inout: f_value 0D_NOT_logical
  if (assc(value)) then
    call c_f_pointer(value, f_value_ptr)
    f_value_ptr = f_value_native
  else
    ! f_value unset
  endif
  ! out: f_ix_match 0D_NOT_integer
  call c_f_pointer(ix_match, f_ix_match_ptr)
  f_ix_match_ptr = f_ix_match
end subroutine
subroutine fortran_coarse_frequency_estimate (data, error, frequency) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data
  type(real_container_alloc), pointer :: f_data
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: frequency  ! 0D_NOT_real
  real(rp) :: f_frequency
  real(c_double), pointer :: f_frequency_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_error
  logical :: f_error_native
  logical(c_bool), pointer :: f_error_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(data))   call c_f_pointer(data, f_data)
  ! inout: f_error 0D_NOT_logical
  if (assc(error)) then
    call c_f_pointer(error, f_error_ptr)
    f_error_native = f_error_ptr
  else
    ! error unset
  endif
  if (assc(error)) then

    f_frequency = coarse_frequency_estimate(data=f_data%data, error=f_error_native)

  else
    f_frequency = coarse_frequency_estimate(data=f_data%data)

  endif
  ! inout: f_error 0D_NOT_logical
  if (assc(error)) then
    call c_f_pointer(error, f_error_ptr)
    f_error_ptr = f_error_native
  else
    ! f_error unset
  endif
  ! out: f_frequency 0D_NOT_real
  call c_f_pointer(frequency, f_frequency_ptr)
  f_frequency_ptr = f_frequency
end subroutine
subroutine fortran_fine_frequency_estimate (data, frequency) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data
  type(real_container_alloc), pointer :: f_data
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: frequency  ! 0D_NOT_real
  real(rp) :: f_frequency
  real(c_double), pointer :: f_frequency_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(data))   call c_f_pointer(data, f_data)
  f_frequency = fine_frequency_estimate(data=f_data%data)

  ! out: f_frequency 0D_NOT_real
  call c_f_pointer(frequency, f_frequency_ptr)
  f_frequency_ptr = f_frequency
end subroutine
subroutine fortran_fourier_amplitude (data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp) &
    bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data
  type(real_container_alloc), pointer :: f_data
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
  !! container general array (1D_ALLOC_real)
  if (assc(data))   call c_f_pointer(data, f_data)
  ! in: f_frequency 0D_NOT_real
  f_frequency = frequency
  if (assc(dcos_amp) .and. assc(dsin_amp)) then

    call fourier_amplitude(data=f_data%data, frequency=f_frequency, cos_amp=f_cos_amp, &
        sin_amp=f_sin_amp, dcos_amp=f_dcos_amp, dsin_amp=f_dsin_amp)

  elseif (assc(dcos_amp)) then

    call fourier_amplitude(data=f_data%data, frequency=f_frequency, cos_amp=f_cos_amp, &
        sin_amp=f_sin_amp, dcos_amp=f_dcos_amp)

  elseif (assc(dsin_amp)) then

    call fourier_amplitude(data=f_data%data, frequency=f_frequency, cos_amp=f_cos_amp, &
        sin_amp=f_sin_amp, dsin_amp=f_dsin_amp)

  else
    call fourier_amplitude(data=f_data%data, frequency=f_frequency, cos_amp=f_cos_amp, &
        sin_amp=f_sin_amp)

  endif
  ! out: f_cos_amp 0D_NOT_real
  call c_f_pointer(cos_amp, f_cos_amp_ptr)
  f_cos_amp_ptr = f_cos_amp
  ! out: f_sin_amp 0D_NOT_real
  call c_f_pointer(sin_amp, f_sin_amp_ptr)
  f_sin_amp_ptr = f_sin_amp
  ! out: f_dcos_amp 0D_NOT_real
  if (assc(dcos_amp)) then
    call c_f_pointer(dcos_amp, f_dcos_amp_ptr)
    f_dcos_amp_ptr = f_dcos_amp
  endif
  ! out: f_dsin_amp 0D_NOT_real
  if (assc(dsin_amp)) then
    call c_f_pointer(dsin_amp, f_dsin_amp_ptr)
    f_dsin_amp_ptr = f_dsin_amp
  endif
end subroutine
subroutine fortran_w_mat_to_axis_angle (w_mat, axis, angle) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  type(c_ptr), intent(in), value :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  real(c_double), pointer :: f_angle_ptr
  ! ** End of parameters **
  !! general array (2D_NOT_real)
  if (assc(w_mat)) then
    call c_f_pointer(w_mat, f_w_mat_ptr, [3*3])
    call vec2mat(f_w_mat_ptr, f_w_mat)
  endif
  call w_mat_to_axis_angle(w_mat=f_w_mat, axis=f_axis, angle=f_angle)

  ! out: f_axis 1D_NOT_real
  if (assc(axis)) f_axis_ptr = f_axis
  ! out: f_angle 0D_NOT_real
  call c_f_pointer(angle, f_angle_ptr)
  f_angle_ptr = f_angle
end subroutine
subroutine fortran_w_mat_to_quat (w_mat, quat) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** End of parameters **
  !! general array (2D_NOT_real)
  if (assc(w_mat)) then
    call c_f_pointer(w_mat, f_w_mat_ptr, [3*3])
    call vec2mat(f_w_mat_ptr, f_w_mat)
  endif
  f_quat = w_mat_to_quat(w_mat=f_w_mat)

  ! out: f_quat 1D_NOT_real
  if (assc(quat)) f_quat_ptr = f_quat
end subroutine
subroutine fortran_quat_to_w_mat (quat, w_mat) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(quat)) then
    call c_f_pointer(quat, f_quat_ptr, [4])
    f_quat = f_quat_ptr(:)
  endif
  f_w_mat = quat_to_w_mat(quat=f_quat)

  ! out: f_w_mat 2D_NOT_real
! TODO: output 2D array: w_mat
end subroutine
subroutine fortran_axis_angle_to_w_mat (axis, angle, w_mat) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: w_mat
  real(rp) :: f_w_mat(3,3)
  real(c_double), pointer :: f_w_mat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(axis)) then
    call c_f_pointer(axis, f_axis_ptr, [3])
    f_axis = f_axis_ptr(:)
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  call axis_angle_to_w_mat(axis=f_axis, angle=f_angle, w_mat=f_w_mat)

  ! out: f_w_mat 2D_NOT_real
! TODO: output 2D array: w_mat
end subroutine
subroutine fortran_quat_to_omega (quat, omega) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: omega
  real(rp) :: f_omega(3)
  real(c_double), pointer :: f_omega_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(quat)) then
    call c_f_pointer(quat, f_quat_ptr, [4])
    f_quat = f_quat_ptr(:)
  endif
  f_omega = quat_to_omega(quat=f_quat)

  ! out: f_omega 1D_NOT_real
  if (assc(omega)) f_omega_ptr = f_omega
end subroutine
subroutine fortran_omega_to_quat (omega, quat) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: omega
  real(rp) :: f_omega(3)
  real(c_double), pointer :: f_omega_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(omega)) then
    call c_f_pointer(omega, f_omega_ptr, [3])
    f_omega = f_omega_ptr(:)
  endif
  f_quat = omega_to_quat(omega=f_omega)

  ! out: f_quat 1D_NOT_real
  if (assc(quat)) f_quat_ptr = f_quat
end subroutine
subroutine fortran_quat_to_axis_angle (quat, axis, angle) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  type(c_ptr), intent(in), value :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  real(c_double), pointer :: f_angle_ptr
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(quat)) then
    call c_f_pointer(quat, f_quat_ptr, [4])
    f_quat = f_quat_ptr(:)
  endif
  call quat_to_axis_angle(quat=f_quat, axis=f_axis, angle=f_angle)

  ! out: f_axis 1D_NOT_real
  if (assc(axis)) f_axis_ptr = f_axis
  ! out: f_angle 0D_NOT_real
  call c_f_pointer(angle, f_angle_ptr)
  f_angle_ptr = f_angle
end subroutine
subroutine fortran_axis_angle_to_quat (axis, angle, quat) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: axis
  real(rp) :: f_axis(3)
  real(c_double), pointer :: f_axis_ptr(:)
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(axis)) then
    call c_f_pointer(axis, f_axis_ptr, [3])
    f_axis = f_axis_ptr(:)
  endif
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  f_quat = axis_angle_to_quat(axis=f_axis, angle=f_angle)

  ! out: f_quat 1D_NOT_real
  if (assc(quat)) f_quat_ptr = f_quat
end subroutine
subroutine fortran_quat_conj_real (q_in, q_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: q_in
  real(rp) :: f_q_in(0:3)
  real(c_double), pointer :: f_q_in_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_out
  real(rp) :: f_q_out(0:3)
  real(c_double), pointer :: f_q_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(q_in)) then
    call c_f_pointer(q_in, f_q_in_ptr, [4])
    f_q_in = f_q_in_ptr(:)
  endif
  f_q_out = quat_conj_real(q_in=f_q_in)

  ! out: f_q_out 1D_NOT_real
  if (assc(q_out)) f_q_out_ptr = f_q_out
end subroutine
subroutine fortran_quat_conj_complex (q_in, q_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: q_in
  complex(rp) :: f_q_in(0:3)
  complex(c_double_complex), pointer :: f_q_in_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_out
  complex(rp) :: f_q_out(0:3)
  complex(c_double_complex), pointer :: f_q_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex)
  if (assc(q_in)) then
    call c_f_pointer(q_in, f_q_in_ptr, [4])
    f_q_in = f_q_in_ptr(:)
  endif
  f_q_out = quat_conj_complex(q_in=f_q_in)

  ! out: f_q_out 1D_NOT_complex
  if (assc(q_out)) f_q_out_ptr = f_q_out
end subroutine
subroutine fortran_quat_inverse (q_in, q_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: q_in
  real(rp) :: f_q_in(0:3)
  real(c_double), pointer :: f_q_in_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_out
  real(rp) :: f_q_out(0:3)
  real(c_double), pointer :: f_q_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(q_in)) then
    call c_f_pointer(q_in, f_q_in_ptr, [4])
    f_q_in = f_q_in_ptr(:)
  endif
  f_q_out = quat_inverse(q_in=f_q_in)

  ! out: f_q_out 1D_NOT_real
  if (assc(q_out)) f_q_out_ptr = f_q_out
end subroutine
subroutine fortran_quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9, q_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: q1
  real(rp) :: f_q1(0:3)
  real(c_double), pointer :: f_q1_ptr(:)
  type(c_ptr), intent(in), value :: q3
  real(rp) :: f_q3(0:3)
  real(c_double), pointer :: f_q3_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_out
  real(rp) :: f_q_out(0:3)
  real(c_double), pointer :: f_q_out_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: q2
  real(rp) :: f_q2(0:3)
  real(c_double), pointer :: f_q2_ptr(:)
  type(c_ptr), intent(in), value :: q4
  real(rp) :: f_q4(0:3)
  real(c_double), pointer :: f_q4_ptr(:)
  type(c_ptr), intent(in), value :: q5
  real(rp) :: f_q5(0:3)
  real(c_double), pointer :: f_q5_ptr(:)
  type(c_ptr), intent(in), value :: q6
  real(rp) :: f_q6(0:3)
  real(c_double), pointer :: f_q6_ptr(:)
  type(c_ptr), intent(in), value :: q7
  real(rp) :: f_q7(0:3)
  real(c_double), pointer :: f_q7_ptr(:)
  type(c_ptr), intent(in), value :: q8
  real(rp) :: f_q8(0:3)
  real(c_double), pointer :: f_q8_ptr(:)
  type(c_ptr), intent(in), value :: q9
  real(rp) :: f_q9(0:3)
  real(c_double), pointer :: f_q9_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(q1)) then
    call c_f_pointer(q1, f_q1_ptr, [4])
    f_q1 = f_q1_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q2)) then
    call c_f_pointer(q2, f_q2_ptr, [4])
    f_q2 = f_q2_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q3)) then
    call c_f_pointer(q3, f_q3_ptr, [4])
    f_q3 = f_q3_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q4)) then
    call c_f_pointer(q4, f_q4_ptr, [4])
    f_q4 = f_q4_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q5)) then
    call c_f_pointer(q5, f_q5_ptr, [4])
    f_q5 = f_q5_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q6)) then
    call c_f_pointer(q6, f_q6_ptr, [4])
    f_q6 = f_q6_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q7)) then
    call c_f_pointer(q7, f_q7_ptr, [4])
    f_q7 = f_q7_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q8)) then
    call c_f_pointer(q8, f_q8_ptr, [4])
    f_q8 = f_q8_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(q9)) then
    call c_f_pointer(q9, f_q9_ptr, [4])
    f_q9 = f_q9_ptr(:)
  endif
  if (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. &
      assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, &
        q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, &
        q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7)

  elseif (assc(q3) .and. assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7)

  elseif (assc(q4) .and. assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7)

  elseif (assc(q4) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q8=f_q8)

  elseif (assc(q4) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q9=f_q9)

  elseif (assc(q4) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4)

  elseif (assc(q3) .and. assc(q5)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5)

  elseif (assc(q3) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6)

  elseif (assc(q3) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7)

  elseif (assc(q3) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q8=f_q8)

  elseif (assc(q3) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3, q9=f_q9)

  elseif (assc(q4) .and. assc(q5)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5)

  elseif (assc(q4) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6)

  elseif (assc(q4) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7)

  elseif (assc(q4) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q8=f_q8)

  elseif (assc(q4) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4, q9=f_q9)

  elseif (assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6)

  elseif (assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7)

  elseif (assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q8=f_q8)

  elseif (assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5, q9=f_q9)

  elseif (assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7)

  elseif (assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q8=f_q8)

  elseif (assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6, q9=f_q9)

  elseif (assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q7=f_q7, q8=f_q8)

  elseif (assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q7=f_q7, q9=f_q9)

  elseif (assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q8=f_q8, q9=f_q9)

  elseif (assc(q3)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q3=f_q3)

  elseif (assc(q4)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q4=f_q4)

  elseif (assc(q5)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q5=f_q5)

  elseif (assc(q6)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q6=f_q6)

  elseif (assc(q7)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q7=f_q7)

  elseif (assc(q8)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q8=f_q8)

  elseif (assc(q9)) then

    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2, q9=f_q9)

  else
    f_q_out = quat_mul_real(q1=f_q1, q2=f_q2)

  endif
  ! out: f_q_out 1D_NOT_real
  if (assc(q_out)) f_q_out_ptr = f_q_out
end subroutine
subroutine fortran_quat_mul_complex (q1, q2, q3, q4, q5, q6, q7, q8, q9, q_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: q1
  complex(rp) :: f_q1(0:3)
  complex(c_double_complex), pointer :: f_q1_ptr(:)
  type(c_ptr), intent(in), value :: q3
  complex(rp) :: f_q3(0:3)
  complex(c_double_complex), pointer :: f_q3_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_out
  complex(rp) :: f_q_out(0:3)
  complex(c_double_complex), pointer :: f_q_out_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: q2
  complex(rp) :: f_q2(0:3)
  complex(c_double_complex), pointer :: f_q2_ptr(:)
  type(c_ptr), intent(in), value :: q4
  complex(rp) :: f_q4(0:3)
  complex(c_double_complex), pointer :: f_q4_ptr(:)
  type(c_ptr), intent(in), value :: q5
  complex(rp) :: f_q5(0:3)
  complex(c_double_complex), pointer :: f_q5_ptr(:)
  type(c_ptr), intent(in), value :: q6
  complex(rp) :: f_q6(0:3)
  complex(c_double_complex), pointer :: f_q6_ptr(:)
  type(c_ptr), intent(in), value :: q7
  complex(rp) :: f_q7(0:3)
  complex(c_double_complex), pointer :: f_q7_ptr(:)
  type(c_ptr), intent(in), value :: q8
  complex(rp) :: f_q8(0:3)
  complex(c_double_complex), pointer :: f_q8_ptr(:)
  type(c_ptr), intent(in), value :: q9
  complex(rp) :: f_q9(0:3)
  complex(c_double_complex), pointer :: f_q9_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex)
  if (assc(q1)) then
    call c_f_pointer(q1, f_q1_ptr, [4])
    f_q1 = f_q1_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q2)) then
    call c_f_pointer(q2, f_q2_ptr, [4])
    f_q2 = f_q2_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q3)) then
    call c_f_pointer(q3, f_q3_ptr, [4])
    f_q3 = f_q3_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q4)) then
    call c_f_pointer(q4, f_q4_ptr, [4])
    f_q4 = f_q4_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q5)) then
    call c_f_pointer(q5, f_q5_ptr, [4])
    f_q5 = f_q5_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q6)) then
    call c_f_pointer(q6, f_q6_ptr, [4])
    f_q6 = f_q6_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q7)) then
    call c_f_pointer(q7, f_q7_ptr, [4])
    f_q7 = f_q7_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q8)) then
    call c_f_pointer(q8, f_q8_ptr, [4])
    f_q8 = f_q8_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(q9)) then
    call c_f_pointer(q9, f_q9_ptr, [4])
    f_q9 = f_q9_ptr(:)
  endif
  if (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. &
      assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, &
        q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, &
        q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) &
      then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8, &
        q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q6=f_q6)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q6) .and. assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4) .and. assc(q5)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q5=f_q5)

  elseif (assc(q3) .and. assc(q4) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q6=f_q6)

  elseif (assc(q3) .and. assc(q4) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q7=f_q7)

  elseif (assc(q3) .and. assc(q4) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q8=f_q8)

  elseif (assc(q3) .and. assc(q4) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4, q9=f_q9)

  elseif (assc(q3) .and. assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q6=f_q6)

  elseif (assc(q3) .and. assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q7=f_q7)

  elseif (assc(q3) .and. assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q8=f_q8)

  elseif (assc(q3) .and. assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5, q9=f_q9)

  elseif (assc(q3) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q7=f_q7)

  elseif (assc(q3) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q8=f_q8)

  elseif (assc(q3) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6, q9=f_q9)

  elseif (assc(q3) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7, q8=f_q8)

  elseif (assc(q3) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7, q9=f_q9)

  elseif (assc(q3) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q8=f_q8, q9=f_q9)

  elseif (assc(q4) .and. assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q6=f_q6)

  elseif (assc(q4) .and. assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q7=f_q7)

  elseif (assc(q4) .and. assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q8=f_q8)

  elseif (assc(q4) .and. assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5, q9=f_q9)

  elseif (assc(q4) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q7=f_q7)

  elseif (assc(q4) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q8=f_q8)

  elseif (assc(q4) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6, q9=f_q9)

  elseif (assc(q4) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7, q8=f_q8)

  elseif (assc(q4) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7, q9=f_q9)

  elseif (assc(q4) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q8=f_q8, q9=f_q9)

  elseif (assc(q5) .and. assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q7=f_q7)

  elseif (assc(q5) .and. assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q8=f_q8)

  elseif (assc(q5) .and. assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6, q9=f_q9)

  elseif (assc(q5) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7, q8=f_q8)

  elseif (assc(q5) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7, q9=f_q9)

  elseif (assc(q5) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q8=f_q8, q9=f_q9)

  elseif (assc(q6) .and. assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7, q8=f_q8)

  elseif (assc(q6) .and. assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7, q9=f_q9)

  elseif (assc(q6) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q8=f_q8, q9=f_q9)

  elseif (assc(q7) .and. assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q7=f_q7, q8=f_q8, q9=f_q9)

  elseif (assc(q3) .and. assc(q4)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q4=f_q4)

  elseif (assc(q3) .and. assc(q5)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q5=f_q5)

  elseif (assc(q3) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q6=f_q6)

  elseif (assc(q3) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q7=f_q7)

  elseif (assc(q3) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q8=f_q8)

  elseif (assc(q3) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3, q9=f_q9)

  elseif (assc(q4) .and. assc(q5)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q5=f_q5)

  elseif (assc(q4) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q6=f_q6)

  elseif (assc(q4) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q7=f_q7)

  elseif (assc(q4) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q8=f_q8)

  elseif (assc(q4) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4, q9=f_q9)

  elseif (assc(q5) .and. assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q6=f_q6)

  elseif (assc(q5) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q7=f_q7)

  elseif (assc(q5) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q8=f_q8)

  elseif (assc(q5) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5, q9=f_q9)

  elseif (assc(q6) .and. assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q7=f_q7)

  elseif (assc(q6) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q8=f_q8)

  elseif (assc(q6) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6, q9=f_q9)

  elseif (assc(q7) .and. assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q7=f_q7, q8=f_q8)

  elseif (assc(q7) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q7=f_q7, q9=f_q9)

  elseif (assc(q8) .and. assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q8=f_q8, q9=f_q9)

  elseif (assc(q3)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q3=f_q3)

  elseif (assc(q4)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q4=f_q4)

  elseif (assc(q5)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q5=f_q5)

  elseif (assc(q6)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q6=f_q6)

  elseif (assc(q7)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q7=f_q7)

  elseif (assc(q8)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q8=f_q8)

  elseif (assc(q9)) then

    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2, q9=f_q9)

  else
    f_q_out = quat_mul_complex(q1=f_q1, q2=f_q2)

  endif
  ! out: f_q_out 1D_NOT_complex
  if (assc(q_out)) f_q_out_ptr = f_q_out
end subroutine
subroutine fortran_quat_rotate_real (quat, vec_in, vec_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: quat
  real(rp) :: f_quat(0:3)
  real(c_double), pointer :: f_quat_ptr(:)
  type(c_ptr), intent(in), value :: vec_in
  real(rp) :: f_vec_in(3)
  real(c_double), pointer :: f_vec_in_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: vec_out
  real(rp) :: f_vec_out(3)
  real(c_double), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(quat)) then
    call c_f_pointer(quat, f_quat_ptr, [4])
    f_quat = f_quat_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(vec_in)) then
    call c_f_pointer(vec_in, f_vec_in_ptr, [3])
    f_vec_in = f_vec_in_ptr(:)
  endif
  f_vec_out = quat_rotate_real(quat=f_quat, vec_in=f_vec_in)

  ! out: f_vec_out 1D_NOT_real
  if (assc(vec_out)) f_vec_out_ptr = f_vec_out
end subroutine
subroutine fortran_quat_rotate_complex (quat, vec_in, vec_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: quat
  complex(rp) :: f_quat(0:3)
  complex(c_double_complex), pointer :: f_quat_ptr(:)
  type(c_ptr), intent(in), value :: vec_in
  complex(rp) :: f_vec_in(3)
  complex(c_double_complex), pointer :: f_vec_in_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: vec_out
  complex(rp) :: f_vec_out(3)
  complex(c_double_complex), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_complex)
  if (assc(quat)) then
    call c_f_pointer(quat, f_quat_ptr, [4])
    f_quat = f_quat_ptr(:)
  endif
  !! general array (1D_NOT_complex)
  if (assc(vec_in)) then
    call c_f_pointer(vec_in, f_vec_in_ptr, [3])
    f_vec_in = f_vec_in_ptr(:)
  endif
  f_vec_out = quat_rotate_complex(quat=f_quat, vec_in=f_vec_in)

  ! out: f_vec_out 1D_NOT_complex
  if (assc(vec_out)) f_vec_out_ptr = f_vec_out
end subroutine
subroutine fortran_rotate_vec_given_axis_angle (vec_in, axis, angle, vec_out) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: vec_in
  real(rp) :: f_vec_in(3)
  real(c_double), pointer :: f_vec_in_ptr(:)
  type(c_ptr), intent(in), value :: axis
  type(real_container_alloc), pointer :: f_axis
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: vec_out
  real(rp) :: f_vec_out(3)
  real(c_double), pointer :: f_vec_out_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(vec_in)) then
    call c_f_pointer(vec_in, f_vec_in_ptr, [3])
    f_vec_in = f_vec_in_ptr(:)
  endif
  !! container general array (1D_ALLOC_real)
  if (assc(axis))   call c_f_pointer(axis, f_axis)
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  f_vec_out = rotate_vec_given_axis_angle(vec_in=f_vec_in, axis=f_axis%data, angle=f_angle)

  ! out: f_vec_out 1D_NOT_real
  if (assc(vec_out)) f_vec_out_ptr = f_vec_out
end subroutine
subroutine fortran_rotate_vec (vec, axis, angle) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: axis  ! 0D_NOT_integer
  integer :: f_axis
  real(c_double) :: angle  ! 0D_NOT_real
  real(rp) :: f_angle
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: vec
  type(real_container_alloc), pointer :: f_vec
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(vec))   call c_f_pointer(vec, f_vec)
  ! in: f_axis 0D_NOT_integer
  f_axis = axis
  ! in: f_angle 0D_NOT_real
  f_angle = angle
  call rotate_vec(vec=f_vec%data, axis=f_axis, angle=f_angle)

end subroutine
subroutine fortran_naff (cdata, freqs, amps, opt_dump_spectra, opt_zero_first) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: cdata
  type(complex_container_alloc), pointer :: f_cdata
  type(c_ptr), intent(in), value :: freqs
  type(real_container_alloc), pointer :: f_freqs
  type(c_ptr), intent(in), value :: amps
  type(complex_container_alloc), pointer :: f_amps
  type(c_ptr), intent(in), value :: opt_dump_spectra  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_spectra
  integer(c_int), pointer :: f_opt_dump_spectra_ptr
  type(c_ptr), intent(in), value :: opt_zero_first  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_opt_zero_first
  logical :: f_opt_zero_first_native
  logical(c_bool), pointer :: f_opt_zero_first_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_complex)
  if (assc(cdata))   call c_f_pointer(cdata, f_cdata)
  !! container general array (1D_ALLOC_real)
  if (assc(freqs))   call c_f_pointer(freqs, f_freqs)
  !! container general array (1D_ALLOC_complex)
  if (assc(amps))   call c_f_pointer(amps, f_amps)
  ! inout: f_opt_dump_spectra 0D_NOT_integer
  if (assc(opt_dump_spectra)) then
    call c_f_pointer(opt_dump_spectra, f_opt_dump_spectra_ptr)
    f_opt_dump_spectra = f_opt_dump_spectra_ptr
  else
    ! opt_dump_spectra unset
  endif
  ! inout: f_opt_zero_first 0D_NOT_logical
  if (assc(opt_zero_first)) then
    call c_f_pointer(opt_zero_first, f_opt_zero_first_ptr)
    f_opt_zero_first_native = f_opt_zero_first_ptr
  else
    ! opt_zero_first unset
  endif
  if (assc(opt_dump_spectra) .and. assc(opt_zero_first)) then

    call naff(cdata=f_cdata%data, freqs=f_freqs%data, amps=f_amps%data, &
        opt_dump_spectra=f_opt_dump_spectra, opt_zero_first=f_opt_zero_first_native)

  elseif (assc(opt_dump_spectra)) then

    call naff(cdata=f_cdata%data, freqs=f_freqs%data, amps=f_amps%data, &
        opt_dump_spectra=f_opt_dump_spectra)

  elseif (assc(opt_zero_first)) then

    call naff(cdata=f_cdata%data, freqs=f_freqs%data, amps=f_amps%data, &
        opt_zero_first=f_opt_zero_first_native)

  else
    call naff(cdata=f_cdata%data, freqs=f_freqs%data, amps=f_amps%data)

  endif
  ! inout: f_opt_dump_spectra 0D_NOT_integer
  if (assc(opt_dump_spectra)) then
    call c_f_pointer(opt_dump_spectra, f_opt_dump_spectra_ptr)
    f_opt_dump_spectra_ptr = f_opt_dump_spectra
  endif
  ! inout: f_opt_zero_first 0D_NOT_logical
  if (assc(opt_zero_first)) then
    call c_f_pointer(opt_zero_first, f_opt_zero_first_ptr)
    f_opt_zero_first_ptr = f_opt_zero_first_native
  else
    ! f_opt_zero_first unset
  endif
end subroutine
subroutine fortran_projdd (a, b, func_retval__) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: func_retval__  ! 0D_NOT_complex
  complex(rp) :: f_func_retval__
  complex(c_double_complex), pointer :: f_func_retval___ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: a
  type(complex_container_alloc), pointer :: f_a
  type(c_ptr), intent(in), value :: b
  type(complex_container_alloc), pointer :: f_b
  ! ** End of parameters **
  !! container general array (1D_ALLOC_complex)
  if (assc(a))   call c_f_pointer(a, f_a)
  !! container general array (1D_ALLOC_complex)
  if (assc(b))   call c_f_pointer(b, f_b)
  f_func_retval__ = projdd(a=f_a%data, b=f_b%data)

  ! out: f_func_retval__ 0D_NOT_complex
  call c_f_pointer(func_retval__, f_func_retval___ptr)
  f_func_retval___ptr = f_func_retval__
end subroutine
subroutine fortran_maximize_projection (seed, cdata, func_retval__) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: func_retval__  ! 0D_NOT_real
  real(rp) :: f_func_retval__
  real(c_double), pointer :: f_func_retval___ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: seed  ! 0D_NOT_real
  real(c_double) :: f_seed
  real(c_double), pointer :: f_seed_ptr
  type(c_ptr), intent(in), value :: cdata
  type(complex_container_alloc), pointer :: f_cdata
  ! ** End of parameters **
  ! inout: f_seed 0D_NOT_real
  if (assc(seed)) then
    call c_f_pointer(seed, f_seed_ptr)
    f_seed = f_seed_ptr
  else
    ! seed unset
  endif
  !! container general array (1D_ALLOC_complex)
  if (assc(cdata))   call c_f_pointer(cdata, f_cdata)
  f_func_retval__ = maximize_projection(seed=f_seed, cdata=f_cdata%data)

  ! inout: f_seed 0D_NOT_real
  if (assc(seed)) then
    call c_f_pointer(seed, f_seed_ptr)
    f_seed_ptr = f_seed
  endif
  ! out: f_func_retval__ 0D_NOT_real
  call c_f_pointer(func_retval__, f_func_retval___ptr)
  f_func_retval___ptr = f_func_retval__
end subroutine
subroutine fortran_interpolated_fft_gsl (cdata, calc_ok, opt_dump_spectrum, opt_dump_index, &
    this_fft) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_fft  ! 0D_NOT_real
  real(rp) :: f_this_fft
  real(c_double), pointer :: f_this_fft_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: cdata
  type(complex_container_alloc), pointer :: f_cdata
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_calc_ok
  logical :: f_calc_ok_native
  logical(c_bool), pointer :: f_calc_ok_ptr
  type(c_ptr), intent(in), value :: opt_dump_spectrum  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_spectrum
  integer(c_int), pointer :: f_opt_dump_spectrum_ptr
  type(c_ptr), intent(in), value :: opt_dump_index  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_index
  integer(c_int), pointer :: f_opt_dump_index_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_complex)
  if (assc(cdata))   call c_f_pointer(cdata, f_cdata)
  ! inout: f_calc_ok 0D_NOT_logical
  if (assc(calc_ok)) then
    call c_f_pointer(calc_ok, f_calc_ok_ptr)
    f_calc_ok_native = f_calc_ok_ptr
  else
    ! calc_ok unset
  endif
  ! inout: f_opt_dump_spectrum 0D_NOT_integer
  if (assc(opt_dump_spectrum)) then
    call c_f_pointer(opt_dump_spectrum, f_opt_dump_spectrum_ptr)
    f_opt_dump_spectrum = f_opt_dump_spectrum_ptr
  else
    ! opt_dump_spectrum unset
  endif
  ! inout: f_opt_dump_index 0D_NOT_integer
  if (assc(opt_dump_index)) then
    call c_f_pointer(opt_dump_index, f_opt_dump_index_ptr)
    f_opt_dump_index = f_opt_dump_index_ptr
  else
    ! opt_dump_index unset
  endif
  if (assc(opt_dump_spectrum) .and. assc(opt_dump_index)) then

    f_this_fft = interpolated_fft_gsl(cdata=f_cdata%data, calc_ok=f_calc_ok_native, &
        opt_dump_spectrum=f_opt_dump_spectrum, opt_dump_index=f_opt_dump_index)

  elseif (assc(opt_dump_spectrum)) then

    f_this_fft = interpolated_fft_gsl(cdata=f_cdata%data, calc_ok=f_calc_ok_native, &
        opt_dump_spectrum=f_opt_dump_spectrum)

  elseif (assc(opt_dump_index)) then

    f_this_fft = interpolated_fft_gsl(cdata=f_cdata%data, calc_ok=f_calc_ok_native, &
        opt_dump_index=f_opt_dump_index)

  else
    f_this_fft = interpolated_fft_gsl(cdata=f_cdata%data, calc_ok=f_calc_ok_native)

  endif
  ! inout: f_calc_ok 0D_NOT_logical
  if (assc(calc_ok)) then
    call c_f_pointer(calc_ok, f_calc_ok_ptr)
    f_calc_ok_ptr = f_calc_ok_native
  else
    ! f_calc_ok unset
  endif
  ! inout: f_opt_dump_spectrum 0D_NOT_integer
  if (assc(opt_dump_spectrum)) then
    call c_f_pointer(opt_dump_spectrum, f_opt_dump_spectrum_ptr)
    f_opt_dump_spectrum_ptr = f_opt_dump_spectrum
  endif
  ! inout: f_opt_dump_index 0D_NOT_integer
  if (assc(opt_dump_index)) then
    call c_f_pointer(opt_dump_index, f_opt_dump_index_ptr)
    f_opt_dump_index_ptr = f_opt_dump_index
  endif
  ! out: f_this_fft 0D_NOT_real
  call c_f_pointer(this_fft, f_this_fft_ptr)
  f_this_fft_ptr = f_this_fft
end subroutine
subroutine fortran_interpolated_fft (cdata, calc_ok, opt_dump_spectrum, opt_dump_index, &
    this_fft) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: this_fft  ! 0D_NOT_real
  real(rp) :: f_this_fft
  real(c_double), pointer :: f_this_fft_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: cdata
  type(complex_container_alloc), pointer :: f_cdata
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_calc_ok
  logical :: f_calc_ok_native
  logical(c_bool), pointer :: f_calc_ok_ptr
  type(c_ptr), intent(in), value :: opt_dump_spectrum  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_spectrum
  integer(c_int), pointer :: f_opt_dump_spectrum_ptr
  type(c_ptr), intent(in), value :: opt_dump_index  ! 0D_NOT_integer
  integer(c_int) :: f_opt_dump_index
  integer(c_int), pointer :: f_opt_dump_index_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_complex)
  if (assc(cdata))   call c_f_pointer(cdata, f_cdata)
  ! inout: f_calc_ok 0D_NOT_logical
  if (assc(calc_ok)) then
    call c_f_pointer(calc_ok, f_calc_ok_ptr)
    f_calc_ok_native = f_calc_ok_ptr
  else
    ! calc_ok unset
  endif
  ! inout: f_opt_dump_spectrum 0D_NOT_integer
  if (assc(opt_dump_spectrum)) then
    call c_f_pointer(opt_dump_spectrum, f_opt_dump_spectrum_ptr)
    f_opt_dump_spectrum = f_opt_dump_spectrum_ptr
  else
    ! opt_dump_spectrum unset
  endif
  ! inout: f_opt_dump_index 0D_NOT_integer
  if (assc(opt_dump_index)) then
    call c_f_pointer(opt_dump_index, f_opt_dump_index_ptr)
    f_opt_dump_index = f_opt_dump_index_ptr
  else
    ! opt_dump_index unset
  endif
  if (assc(opt_dump_spectrum) .and. assc(opt_dump_index)) then

    f_this_fft = interpolated_fft(cdata=f_cdata%data, calc_ok=f_calc_ok_native, &
        opt_dump_spectrum=f_opt_dump_spectrum, opt_dump_index=f_opt_dump_index)

  elseif (assc(opt_dump_spectrum)) then

    f_this_fft = interpolated_fft(cdata=f_cdata%data, calc_ok=f_calc_ok_native, &
        opt_dump_spectrum=f_opt_dump_spectrum)

  elseif (assc(opt_dump_index)) then

    f_this_fft = interpolated_fft(cdata=f_cdata%data, calc_ok=f_calc_ok_native, &
        opt_dump_index=f_opt_dump_index)

  else
    f_this_fft = interpolated_fft(cdata=f_cdata%data, calc_ok=f_calc_ok_native)

  endif
  ! inout: f_calc_ok 0D_NOT_logical
  if (assc(calc_ok)) then
    call c_f_pointer(calc_ok, f_calc_ok_ptr)
    f_calc_ok_ptr = f_calc_ok_native
  else
    ! f_calc_ok unset
  endif
  ! inout: f_opt_dump_spectrum 0D_NOT_integer
  if (assc(opt_dump_spectrum)) then
    call c_f_pointer(opt_dump_spectrum, f_opt_dump_spectrum_ptr)
    f_opt_dump_spectrum_ptr = f_opt_dump_spectrum
  endif
  ! inout: f_opt_dump_index 0D_NOT_integer
  if (assc(opt_dump_index)) then
    call c_f_pointer(opt_dump_index, f_opt_dump_index_ptr)
    f_opt_dump_index_ptr = f_opt_dump_index
  endif
  ! out: f_this_fft 0D_NOT_real
  call c_f_pointer(this_fft, f_this_fft_ptr)
  f_this_fft_ptr = f_this_fft
end subroutine
subroutine fortran_initfixedwindowls (N, dt, order, der, id) bind(c)

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
  f_id = initfixedwindowls(N=f_N, dt=f_dt, order=f_order, der=f_der)

  ! out: f_id 0D_NOT_integer
  call c_f_pointer(id, f_id_ptr)
  f_id_ptr = f_id
end subroutine
subroutine fortran_destfixedwindowls (id) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: id  ! 0D_NOT_integer
  integer(c_int) :: f_id
  integer(c_int), pointer :: f_id_ptr
  ! ** End of parameters **
  ! inout: f_id 0D_NOT_integer
  if (assc(id)) then
    call c_f_pointer(id, f_id_ptr)
    f_id = f_id_ptr
  else
    ! id unset
  endif
  call destfixedwindowls(id=f_id)

  ! inout: f_id 0D_NOT_integer
  if (assc(id)) then
    call c_f_pointer(id, f_id_ptr)
    f_id_ptr = f_id
  endif
end subroutine
subroutine fortran_fixedwindowls (ynew, id, z) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: z  ! 0D_NOT_real
  REAL(rp) :: f_z
  real(c_double), pointer :: f_z_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ynew  ! 0D_NOT_real
  real(c_double) :: f_ynew
  real(c_double), pointer :: f_ynew_ptr
  type(c_ptr), intent(in), value :: id  ! 0D_NOT_integer
  integer(c_int) :: f_id
  integer(c_int), pointer :: f_id_ptr
  ! ** End of parameters **
  ! inout: f_ynew 0D_NOT_real
  if (assc(ynew)) then
    call c_f_pointer(ynew, f_ynew_ptr)
    f_ynew = f_ynew_ptr
  else
    ! ynew unset
  endif
  ! inout: f_id 0D_NOT_integer
  if (assc(id)) then
    call c_f_pointer(id, f_id_ptr)
    f_id = f_id_ptr
  else
    ! id unset
  endif
  f_z = fixedwindowls(ynew=f_ynew, id=f_id)

  ! inout: f_ynew 0D_NOT_real
  if (assc(ynew)) then
    call c_f_pointer(ynew, f_ynew_ptr)
    f_ynew_ptr = f_ynew
  endif
  ! inout: f_id 0D_NOT_integer
  if (assc(id)) then
    call c_f_pointer(id, f_id_ptr)
    f_id_ptr = f_id
  endif
  ! out: f_z 0D_NOT_real
  call c_f_pointer(z, f_z_ptr)
  f_z_ptr = f_z
end subroutine
subroutine fortran_super_bicubic_interpolation (y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, &
    ansy, ansy1, ansy2) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: y
  real(rp) :: f_y(4)
  real(c_double), pointer :: f_y_ptr(:)
  type(c_ptr), intent(in), value :: y1
  real(rp) :: f_y1(4)
  real(c_double), pointer :: f_y1_ptr(:)
  type(c_ptr), intent(in), value :: y2
  real(rp) :: f_y2(4)
  real(c_double), pointer :: f_y2_ptr(:)
  type(c_ptr), intent(in), value :: y12
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
  !! general array (1D_NOT_real)
  if (assc(y)) then
    call c_f_pointer(y, f_y_ptr, [4])
    f_y = f_y_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(y1)) then
    call c_f_pointer(y1, f_y1_ptr, [4])
    f_y1 = f_y1_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(y2)) then
    call c_f_pointer(y2, f_y2_ptr, [4])
    f_y2 = f_y2_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(y12)) then
    call c_f_pointer(y12, f_y12_ptr, [4])
    f_y12 = f_y12_ptr(:)
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
  call super_bicubic_interpolation(y=f_y, y1=f_y1, y2=f_y2, y12=f_y12, x1l=f_x1l, x1u=f_x1u, &
      x2l=f_x2l, x2u=f_x2u, x1=f_x1, x2=f_x2, ansy=f_ansy, ansy1=f_ansy1, ansy2=f_ansy2)

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
subroutine fortran_super_bicubic_coef (y, y1, y2, y12, d1, d2, c) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: y
  real(dp) :: f_y(4)
  real(c_double), pointer :: f_y_ptr(:)
  type(c_ptr), intent(in), value :: y1
  real(dp) :: f_y1(4)
  real(c_double), pointer :: f_y1_ptr(:)
  type(c_ptr), intent(in), value :: y2
  real(dp) :: f_y2(4)
  real(c_double), pointer :: f_y2_ptr(:)
  type(c_ptr), intent(in), value :: y12
  real(dp) :: f_y12(4)
  real(c_double), pointer :: f_y12_ptr(:)
  real(c_double) :: d1  ! 0D_NOT_real
  real(dp) :: f_d1
  real(c_double) :: d2  ! 0D_NOT_real
  real(dp) :: f_d2
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: c
  real(dp) :: f_c(4,4)
  real(c_double), pointer :: f_c_ptr(:)
  ! ** End of parameters **
  !! general array (1D_NOT_real)
  if (assc(y)) then
    call c_f_pointer(y, f_y_ptr, [4])
    f_y = f_y_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(y1)) then
    call c_f_pointer(y1, f_y1_ptr, [4])
    f_y1 = f_y1_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(y2)) then
    call c_f_pointer(y2, f_y2_ptr, [4])
    f_y2 = f_y2_ptr(:)
  endif
  !! general array (1D_NOT_real)
  if (assc(y12)) then
    call c_f_pointer(y12, f_y12_ptr, [4])
    f_y12 = f_y12_ptr(:)
  endif
  ! in: f_d1 0D_NOT_real
  f_d1 = d1
  ! in: f_d2 0D_NOT_real
  f_d2 = d2
  call super_bicubic_coef(y=f_y, y1=f_y1, y2=f_y2, y12=f_y12, d1=f_d1, d2=f_d2, c=f_c)

  ! out: f_c 2D_NOT_real
! TODO: output 2D array: c
end subroutine
subroutine fortran_super_sort (arr) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr
  type(integer_container_alloc), pointer :: f_arr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_integer)
  if (assc(arr))   call c_f_pointer(arr, f_arr)
  call super_sort(arr=f_arr%data)

end subroutine
subroutine fortran_super_polint (xa, ya, x, y, dy) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: xa
  type(real_container_alloc), pointer :: f_xa
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y  ! 0D_NOT_real
  real(rp) :: f_y
  real(c_double), pointer :: f_y_ptr
  type(c_ptr), intent(in), value :: dy  ! 0D_NOT_real
  real(rp) :: f_dy
  real(c_double), pointer :: f_dy_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ya
  type(real_container_alloc), pointer :: f_ya
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(xa))   call c_f_pointer(xa, f_xa)
  !! container general array (1D_ALLOC_real)
  if (assc(ya))   call c_f_pointer(ya, f_ya)
  ! in: f_x 0D_NOT_real
  f_x = x
  call super_polint(xa=f_xa%data, ya=f_ya%data, x=f_x, y=f_y, dy=f_dy)

  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
  ! out: f_dy 0D_NOT_real
  call c_f_pointer(dy, f_dy_ptr)
  f_dy_ptr = f_dy
end subroutine
subroutine fortran_super_poly (x, coeffs, value) bind(c)

  implicit none
  ! ** In parameters **
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: coeffs
  type(real_container_alloc), pointer :: f_coeffs
  ! ** End of parameters **
  ! in: f_x 0D_NOT_real
  f_x = x
  !! container general array (1D_ALLOC_real)
  if (assc(coeffs))   call c_f_pointer(coeffs, f_coeffs)
  f_value = super_poly(x=f_x, coeffs=f_coeffs%data)

  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_ran_seed_put (seed, mpi_offset) bind(c)

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
  if (assc(mpi_offset)) then
    call c_f_pointer(mpi_offset, f_mpi_offset_ptr)
    f_mpi_offset = f_mpi_offset_ptr
  else
    ! mpi_offset unset
  endif
  if (assc(mpi_offset)) then

    call ran_seed_put(seed=f_seed, mpi_offset=f_mpi_offset)

  else
    call ran_seed_put(seed=f_seed)

  endif
end subroutine
subroutine fortran_ran_seed_get (seed) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: seed  ! 0D_NOT_integer
  integer :: f_seed
  integer(c_int), pointer :: f_seed_ptr
  ! ** End of parameters **
  call ran_seed_get(seed=f_seed)

  ! out: f_seed 0D_NOT_integer
  call c_f_pointer(seed, f_seed_ptr)
  f_seed_ptr = f_seed
end subroutine
subroutine fortran_allocate_thread_states () bind(c)

  implicit none
  ! ** End of parameters **
  call allocate_thread_states()

end subroutine
subroutine fortran_bicubic_cmplx_eval (x_norm, y_norm, bi_coef, df_dx, df_dy, f_val) bind(c)

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
  if (.not. assc(bi_coef)) return
  call c_f_pointer(bi_coef, f_bi_coef)
  if (assc(df_dx) .and. assc(df_dy)) then

    f_f_val = bicubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, bi_coef=f_bi_coef, &
        df_dx=f_df_dx, df_dy=f_df_dy)

  elseif (assc(df_dx)) then

    f_f_val = bicubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, bi_coef=f_bi_coef, &
        df_dx=f_df_dx)

  elseif (assc(df_dy)) then

    f_f_val = bicubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, bi_coef=f_bi_coef, &
        df_dy=f_df_dy)

  else
    f_f_val = bicubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, bi_coef=f_bi_coef)

  endif
  ! out: f_df_dx 0D_NOT_complex
  if (assc(df_dx)) then
    call c_f_pointer(df_dx, f_df_dx_ptr)
    f_df_dx_ptr = f_df_dx
  endif
  ! out: f_df_dy 0D_NOT_complex
  if (assc(df_dy)) then
    call c_f_pointer(df_dy, f_df_dy_ptr)
    f_df_dy_ptr = f_df_dy
  endif
  ! out: f_f_val 0D_NOT_complex
  call c_f_pointer(f_val, f_f_val_ptr)
  f_f_val_ptr = f_f_val
end subroutine
subroutine fortran_tricubic_cmplx_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz, &
    f_val) bind(c)

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
  if (.not. assc(tri_coef)) return
  call c_f_pointer(tri_coef, f_tri_coef)
  if (assc(df_dx) .and. assc(df_dy) .and. assc(df_dz)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dx=f_df_dx, df_dy=f_df_dy, df_dz=f_df_dz)

  elseif (assc(df_dx) .and. assc(df_dy)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dx=f_df_dx, df_dy=f_df_dy)

  elseif (assc(df_dx) .and. assc(df_dz)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dx=f_df_dx, df_dz=f_df_dz)

  elseif (assc(df_dy) .and. assc(df_dz)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dy=f_df_dy, df_dz=f_df_dz)

  elseif (assc(df_dx)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dx=f_df_dx)

  elseif (assc(df_dy)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dy=f_df_dy)

  elseif (assc(df_dz)) then

    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef, df_dz=f_df_dz)

  else
    f_f_val = tricubic_cmplx_eval(x_norm=f_x_norm, y_norm=f_y_norm, z_norm=f_z_norm, &
        tri_coef=f_tri_coef)

  endif
  ! out: f_df_dx 0D_NOT_complex
  if (assc(df_dx)) then
    call c_f_pointer(df_dx, f_df_dx_ptr)
    f_df_dx_ptr = f_df_dx
  endif
  ! out: f_df_dy 0D_NOT_complex
  if (assc(df_dy)) then
    call c_f_pointer(df_dy, f_df_dy_ptr)
    f_df_dy_ptr = f_df_dy
  endif
  ! out: f_df_dz 0D_NOT_complex
  if (assc(df_dz)) then
    call c_f_pointer(df_dz, f_df_dz_ptr)
    f_df_dz_ptr = f_df_dz
  endif
  ! out: f_f_val 0D_NOT_complex
  call c_f_pointer(f_val, f_f_val_ptr)
  f_f_val_ptr = f_f_val
end subroutine
subroutine fortran_bin_index (x, bin1_x_min, bin_delta, ix_bin) bind(c)

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
  f_ix_bin = bin_index(x=f_x, bin1_x_min=f_bin1_x_min, bin_delta=f_bin_delta)

  ! out: f_ix_bin 0D_NOT_integer
  call c_f_pointer(ix_bin, f_ix_bin_ptr)
  f_ix_bin_ptr = f_ix_bin
end subroutine
subroutine fortran_bin_x_center (ix_bin, bin1_x_min, bin_delta, x_center) bind(c)

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
  if (assc(ix_bin)) then
    call c_f_pointer(ix_bin, f_ix_bin_ptr)
    f_ix_bin = f_ix_bin_ptr
  else
    ! ix_bin unset
  endif
  ! in: f_bin1_x_min 0D_NOT_real
  f_bin1_x_min = bin1_x_min
  ! in: f_bin_delta 0D_NOT_real
  f_bin_delta = bin_delta
  f_x_center = bin_x_center(ix_bin=f_ix_bin, bin1_x_min=f_bin1_x_min, bin_delta=f_bin_delta)

  ! inout: f_ix_bin 0D_NOT_integer
  if (assc(ix_bin)) then
    call c_f_pointer(ix_bin, f_ix_bin_ptr)
    f_ix_bin_ptr = f_ix_bin
  endif
  ! out: f_x_center 0D_NOT_real
  call c_f_pointer(x_center, f_x_center_ptr)
  f_x_center_ptr = f_x_center
end subroutine
subroutine fortran_n_bins_automatic (n_data, n) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: n  ! 0D_NOT_integer
  integer :: f_n
  integer(c_int), pointer :: f_n_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: n_data  ! 0D_NOT_integer
  integer(c_int) :: f_n_data
  integer(c_int), pointer :: f_n_data_ptr
  ! ** End of parameters **
  ! inout: f_n_data 0D_NOT_integer
  if (assc(n_data)) then
    call c_f_pointer(n_data, f_n_data_ptr)
    f_n_data = f_n_data_ptr
  else
    ! n_data unset
  endif
  f_n = n_bins_automatic(n_data=f_n_data)

  ! inout: f_n_data 0D_NOT_integer
  if (assc(n_data)) then
    call c_f_pointer(n_data, f_n_data_ptr)
    f_n_data_ptr = f_n_data
  endif
  ! out: f_n 0D_NOT_integer
  call c_f_pointer(n, f_n_ptr)
  f_n_ptr = f_n
end subroutine
subroutine fortran_reallocate_spline (spline, n, n_min, exact) bind(c)

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
  logical :: f_exact_native
  logical(c_bool), pointer :: f_exact_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: spline
  type(spline_struct_container_alloc), pointer :: f_spline
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (assc(spline))   call c_f_pointer(spline, f_spline)
  ! in: f_n 0D_NOT_integer
  f_n = n
  ! in: f_n_min 0D_NOT_integer
  if (assc(n_min)) then
    call c_f_pointer(n_min, f_n_min_ptr)
    f_n_min = f_n_min_ptr
  else
    ! n_min unset
  endif
  ! in: f_exact 0D_NOT_logical
  if (assc(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_native = f_exact_ptr
  else
    ! exact unset
  endif
  if (assc(n_min) .and. assc(exact)) then

    call reallocate_spline(spline=f_spline%data, n=f_n, n_min=f_n_min, exact=f_exact_native)

  elseif (assc(n_min)) then

    call reallocate_spline(spline=f_spline%data, n=f_n, n_min=f_n_min)

  elseif (assc(exact)) then

    call reallocate_spline(spline=f_spline%data, n=f_n, exact=f_exact_native)

  else
    call reallocate_spline(spline=f_spline%data, n=f_n)

  endif
end subroutine
subroutine fortran_create_a_spline (r0, r1, slope0, slope1, spline) bind(c)

  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: r0
  type(real_container_alloc), pointer :: f_r0
  type(c_ptr), intent(in), value :: r1
  type(real_container_alloc), pointer :: f_r1
  real(c_double) :: slope0  ! 0D_NOT_real
  real(rp) :: f_slope0
  real(c_double) :: slope1  ! 0D_NOT_real
  real(rp) :: f_slope1
  ! ** Out parameters **
  type(c_ptr), value :: spline  ! 0D_NOT_type
  type(spline_struct), pointer :: f_spline
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(r0))   call c_f_pointer(r0, f_r0)
  !! container general array (1D_ALLOC_real)
  if (assc(r1))   call c_f_pointer(r1, f_r1)
  ! in: f_slope0 0D_NOT_real
  f_slope0 = slope0
  ! in: f_slope1 0D_NOT_real
  f_slope1 = slope1
  f_spline = create_a_spline(r0=f_r0%data, r1=f_r1%data, slope0=f_slope0, slope1=f_slope1)

  ! out: f_spline 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_spline_akima_interpolate (x_knot, y_knot, x, ok, y, dy) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: x_knot
  type(real_container_alloc), pointer :: f_x_knot
  type(c_ptr), intent(in), value :: y_knot
  type(real_container_alloc), pointer :: f_y_knot
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
  !! container general array (1D_ALLOC_real)
  if (assc(x_knot))   call c_f_pointer(x_knot, f_x_knot)
  !! container general array (1D_ALLOC_real)
  if (assc(y_knot))   call c_f_pointer(y_knot, f_y_knot)
  ! in: f_x 0D_NOT_real
  f_x = x
  if (assc(y) .and. assc(dy)) then

    call spline_akima_interpolate(x_knot=f_x_knot%data, y_knot=f_y_knot%data, x=f_x, ok=f_ok, &
        y=f_y, dy=f_dy)

  elseif (assc(y)) then

    call spline_akima_interpolate(x_knot=f_x_knot%data, y_knot=f_y_knot%data, x=f_x, ok=f_ok, &
        y=f_y)

  elseif (assc(dy)) then

    call spline_akima_interpolate(x_knot=f_x_knot%data, y_knot=f_y_knot%data, x=f_x, ok=f_ok, &
        dy=f_dy)

  else
    call spline_akima_interpolate(x_knot=f_x_knot%data, y_knot=f_y_knot%data, x=f_x, ok=f_ok)

  endif
  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
  ! out: f_y 0D_NOT_real
  if (assc(y)) then
    call c_f_pointer(y, f_y_ptr)
    f_y_ptr = f_y
  endif
  ! out: f_dy 0D_NOT_real
  if (assc(dy)) then
    call c_f_pointer(dy, f_dy_ptr)
    f_dy_ptr = f_dy
  endif
end subroutine
subroutine fortran_spline_evaluate (spline, x, ok, y, dy) bind(c)

  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: spline
  type(spline_struct_container_alloc), pointer :: f_spline
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
  !! container type array (1D_ALLOC_type)
  if (assc(spline))   call c_f_pointer(spline, f_spline)
  ! in: f_x 0D_NOT_real
  f_x = x
  if (assc(y) .and. assc(dy)) then

    call spline_evaluate(spline=f_spline%data, x=f_x, ok=f_ok, y=f_y, dy=f_dy)

  elseif (assc(y)) then

    call spline_evaluate(spline=f_spline%data, x=f_x, ok=f_ok, y=f_y)

  elseif (assc(dy)) then

    call spline_evaluate(spline=f_spline%data, x=f_x, ok=f_ok, dy=f_dy)

  else
    call spline_evaluate(spline=f_spline%data, x=f_x, ok=f_ok)

  endif
  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
  ! out: f_y 0D_NOT_real
  if (assc(y)) then
    call c_f_pointer(y, f_y_ptr)
    f_y_ptr = f_y
  endif
  ! out: f_dy 0D_NOT_real
  if (assc(dy)) then
    call c_f_pointer(dy, f_dy_ptr)
    f_dy_ptr = f_dy
  endif
end subroutine
subroutine fortran_bracket_index_for_spline (x_knot, x, ix0, strict, print_err, ok) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: x_knot
  type(real_container_alloc), pointer :: f_x_knot
  real(c_double) :: x  ! 0D_NOT_real
  real(rp) :: f_x
  type(c_ptr), intent(in), value :: strict  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_strict
  logical :: f_strict_native
  logical(c_bool), pointer :: f_strict_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical :: f_print_err_native
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix0  ! 0D_NOT_integer
  integer :: f_ix0
  integer(c_int), pointer :: f_ix0_ptr
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(x_knot))   call c_f_pointer(x_knot, f_x_knot)
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_strict 0D_NOT_logical
  if (assc(strict)) then
    call c_f_pointer(strict, f_strict_ptr)
    f_strict_native = f_strict_ptr
  else
    ! strict unset
  endif
  ! in: f_print_err 0D_NOT_logical
  if (assc(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
  else
    ! print_err unset
  endif
  if (assc(strict) .and. assc(print_err)) then

    f_ok = bracket_index_for_spline(x_knot=f_x_knot%data, x=f_x, ix0=f_ix0, &
        strict=f_strict_native, print_err=f_print_err_native)

  elseif (assc(strict)) then

    f_ok = bracket_index_for_spline(x_knot=f_x_knot%data, x=f_x, ix0=f_ix0, &
        strict=f_strict_native)

  elseif (assc(print_err)) then

    f_ok = bracket_index_for_spline(x_knot=f_x_knot%data, x=f_x, ix0=f_ix0, &
        print_err=f_print_err_native)

  else
    f_ok = bracket_index_for_spline(x_knot=f_x_knot%data, x=f_x, ix0=f_ix0)

  endif
  ! out: f_ix0 0D_NOT_integer
  call c_f_pointer(ix0, f_ix0_ptr)
  f_ix0_ptr = f_ix0
  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
end subroutine
subroutine fortran_spline1 (a_spline, x, n, y) bind(c)

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
  if (.not. assc(a_spline)) return
  call c_f_pointer(a_spline, f_a_spline)
  ! in: f_x 0D_NOT_real
  f_x = x
  ! in: f_n 0D_NOT_integer
  if (assc(n)) then
    call c_f_pointer(n, f_n_ptr)
    f_n = f_n_ptr
  else
    ! n unset
  endif
  if (assc(n)) then

    f_y = spline1(a_spline=f_a_spline, x=f_x, n=f_n)

  else
    f_y = spline1(a_spline=f_a_spline, x=f_x)

  endif
  ! out: f_y 0D_NOT_real
  call c_f_pointer(y, f_y_ptr)
  f_y_ptr = f_y
end subroutine
subroutine fortran_spline_akima (spline, ok) bind(c)

  use spline_mod, only: spline_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: spline
  type(spline_struct_container_alloc), pointer :: f_spline
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (assc(spline))   call c_f_pointer(spline, f_spline)
  call spline_akima(spline=f_spline%data, ok=f_ok)

  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
end subroutine
subroutine fortran_end_akima_spline_calc (spline, which_end) bind(c)

  use spline_mod, only: spline_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: which_end  ! 0D_NOT_integer
  integer :: f_which_end
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: spline
  type(spline_struct_container_alloc), pointer :: f_spline
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (assc(spline))   call c_f_pointer(spline, f_spline)
  ! in: f_which_end 0D_NOT_integer
  f_which_end = which_end
  call end_akima_spline_calc(spline=f_spline%data, which_end=f_which_end)

end subroutine
subroutine fortran_apfft_corr (rdata_in, bounds, window, phase, amp, freq, diag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: rdata_in
  type(real_container_alloc), pointer :: f_rdata_in
  type(c_ptr), intent(in), value :: bounds
  real(rp) :: f_bounds(2)
  real(c_double), pointer :: f_bounds_ptr(:)
  type(c_ptr), intent(in), value :: window
  character(len=4096) :: f_window
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
  !! container general array (1D_ALLOC_real)
  if (assc(rdata_in))   call c_f_pointer(rdata_in, f_rdata_in)
  !! general array (1D_NOT_real)
  if (assc(bounds)) then
    call c_f_pointer(bounds, f_bounds_ptr, [2])
    f_bounds = f_bounds_ptr(:)
  endif
  ! in: f_window 0D_NOT_character
  if (.not. assc(window)) return
  call c_f_pointer(window, f_window_ptr, [huge(0)])
  call to_f_str(f_window_ptr, f_window)
  ! in: f_diag 0D_NOT_integer
  if (assc(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
    f_diag = f_diag_ptr
  else
    ! diag unset
  endif
  if (assc(bounds) .and. assc(diag)) then

    call apfft_corr(rdata_in=f_rdata_in%data, window=f_window, phase=f_phase, amp=f_amp, &
        freq=f_freq, bounds=f_bounds, diag=f_diag)

  elseif (assc(bounds)) then

    call apfft_corr(rdata_in=f_rdata_in%data, window=f_window, phase=f_phase, amp=f_amp, &
        freq=f_freq, bounds=f_bounds)

  elseif (assc(diag)) then

    call apfft_corr(rdata_in=f_rdata_in%data, window=f_window, phase=f_phase, amp=f_amp, &
        freq=f_freq, diag=f_diag)

  else
    call apfft_corr(rdata_in=f_rdata_in%data, window=f_window, phase=f_phase, amp=f_amp, &
        freq=f_freq)

  endif
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
subroutine fortran_apfft (rdata_in, bounds, window, phase, diag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: rdata_in
  type(real_container_alloc), pointer :: f_rdata_in
  type(c_ptr), intent(in), value :: bounds
  real(rp) :: f_bounds(2)
  real(c_double), pointer :: f_bounds_ptr(:)
  type(c_ptr), intent(in), value :: window
  character(len=4096) :: f_window
  character(kind=c_char), pointer :: f_window_ptr(:)
  type(c_ptr), intent(in), value :: phase  ! 0D_NOT_real
  real(c_double) :: f_phase
  real(c_double), pointer :: f_phase_ptr
  type(c_ptr), intent(in), value :: diag  ! 0D_NOT_integer
  integer(c_int) :: f_diag
  integer(c_int), pointer :: f_diag_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(rdata_in))   call c_f_pointer(rdata_in, f_rdata_in)
  !! general array (1D_NOT_real)
  if (assc(bounds)) then
    call c_f_pointer(bounds, f_bounds_ptr, [2])
    f_bounds = f_bounds_ptr(:)
  endif
  ! inout: f_window 0D_NOT_character
  if (.not. assc(window)) return
  call c_f_pointer(window, f_window_ptr, [huge(0)])
  call to_f_str(f_window_ptr, f_window)
  ! inout: f_phase 0D_NOT_real
  if (assc(phase)) then
    call c_f_pointer(phase, f_phase_ptr)
    f_phase = f_phase_ptr
  else
    ! phase unset
  endif
  ! inout: f_diag 0D_NOT_integer
  if (assc(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
    f_diag = f_diag_ptr
  else
    ! diag unset
  endif
  if (assc(diag)) then

    call apfft(rdata_in=f_rdata_in%data, bounds=f_bounds, window=f_window, phase=f_phase, &
        diag=f_diag)

  else
    call apfft(rdata_in=f_rdata_in%data, bounds=f_bounds, window=f_window, phase=f_phase)

  endif
  ! inout: f_window 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_phase 0D_NOT_real
  if (assc(phase)) then
    call c_f_pointer(phase, f_phase_ptr)
    f_phase_ptr = f_phase
  endif
  ! inout: f_diag 0D_NOT_integer
  if (assc(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
    f_diag_ptr = f_diag
  endif
end subroutine
subroutine fortran_apfft_ext (rdata, bounds, window, phase, amp, freq, diag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: rdata
  type(real_container_alloc), pointer :: f_rdata
  type(c_ptr), intent(in), value :: bounds
  real(rp) :: f_bounds(2)
  real(c_double), pointer :: f_bounds_ptr(:)
  type(c_ptr), intent(in), value :: window
  character(len=4096) :: f_window
  character(kind=c_char), pointer :: f_window_ptr(:)
  type(c_ptr), intent(in), value :: phase  ! 0D_NOT_real
  real(c_double) :: f_phase
  real(c_double), pointer :: f_phase_ptr
  type(c_ptr), intent(in), value :: amp  ! 0D_NOT_real
  real(c_double) :: f_amp
  real(c_double), pointer :: f_amp_ptr
  type(c_ptr), intent(in), value :: freq  ! 0D_NOT_real
  real(c_double) :: f_freq
  real(c_double), pointer :: f_freq_ptr
  type(c_ptr), intent(in), value :: diag  ! 0D_NOT_integer
  integer(c_int) :: f_diag
  integer(c_int), pointer :: f_diag_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (assc(rdata))   call c_f_pointer(rdata, f_rdata)
  !! general array (1D_NOT_real)
  if (assc(bounds)) then
    call c_f_pointer(bounds, f_bounds_ptr, [2])
    f_bounds = f_bounds_ptr(:)
  endif
  ! inout: f_window 0D_NOT_character
  if (.not. assc(window)) return
  call c_f_pointer(window, f_window_ptr, [huge(0)])
  call to_f_str(f_window_ptr, f_window)
  ! inout: f_phase 0D_NOT_real
  if (assc(phase)) then
    call c_f_pointer(phase, f_phase_ptr)
    f_phase = f_phase_ptr
  else
    ! phase unset
  endif
  ! inout: f_amp 0D_NOT_real
  if (assc(amp)) then
    call c_f_pointer(amp, f_amp_ptr)
    f_amp = f_amp_ptr
  else
    ! amp unset
  endif
  ! inout: f_freq 0D_NOT_real
  if (assc(freq)) then
    call c_f_pointer(freq, f_freq_ptr)
    f_freq = f_freq_ptr
  else
    ! freq unset
  endif
  ! inout: f_diag 0D_NOT_integer
  if (assc(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
    f_diag = f_diag_ptr
  else
    ! diag unset
  endif
  if (assc(diag)) then

    call apfft_ext(rdata=f_rdata%data, bounds=f_bounds, window=f_window, phase=f_phase, &
        amp=f_amp, freq=f_freq, diag=f_diag)

  else
    call apfft_ext(rdata=f_rdata%data, bounds=f_bounds, window=f_window, phase=f_phase, &
        amp=f_amp, freq=f_freq)

  endif
  ! inout: f_window 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_phase 0D_NOT_real
  if (assc(phase)) then
    call c_f_pointer(phase, f_phase_ptr)
    f_phase_ptr = f_phase
  endif
  ! inout: f_amp 0D_NOT_real
  if (assc(amp)) then
    call c_f_pointer(amp, f_amp_ptr)
    f_amp_ptr = f_amp
  endif
  ! inout: f_freq 0D_NOT_real
  if (assc(freq)) then
    call c_f_pointer(freq, f_freq_ptr)
    f_freq_ptr = f_freq
  endif
  ! inout: f_diag 0D_NOT_integer
  if (assc(diag)) then
    call c_f_pointer(diag, f_diag_ptr)
    f_diag_ptr = f_diag
  endif
end subroutine
subroutine fortran_hanhan (N, hh) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: N  ! 0D_NOT_integer
  integer(c_int) :: f_N
  integer(c_int), pointer :: f_N_ptr
  type(c_ptr), intent(in), value :: hh
  type(real_container_alloc), pointer :: f_hh
  ! ** End of parameters **
  ! inout: f_N 0D_NOT_integer
  if (assc(N)) then
    call c_f_pointer(N, f_N_ptr)
    f_N = f_N_ptr
  else
    ! N unset
  endif
  !! container general array (1D_ALLOC_real)
  if (assc(hh))   call c_f_pointer(hh, f_hh)
  call hanhan(N=f_N, hh=f_hh%data)

  ! inout: f_N 0D_NOT_integer
  if (assc(N)) then
    call c_f_pointer(N, f_N_ptr)
    f_N_ptr = f_N
  endif
end subroutine
subroutine fortran_bit_set (word, pos, set_to_1) bind(c)

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
  if (assc(word)) then
    call c_f_pointer(word, f_word_ptr)
    f_word = f_word_ptr
  else
    ! word unset
  endif
  ! in: f_pos 0D_NOT_integer
  f_pos = pos
  ! in: f_set_to_1 0D_NOT_logical
  f_set_to_1 = set_to_1
  call bit_set(word=f_word, pos=f_pos, set_to_1=f_set_to_1)

  ! inout: f_word 0D_NOT_integer
  if (assc(word)) then
    call c_f_pointer(word, f_word_ptr)
    f_word_ptr = f_word
  endif
end subroutine

end module cppbmad_sim_utils_routines
