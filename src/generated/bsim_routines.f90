module cppbmad_bsim_routines

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

use bbu_track_mod, only: bbu_add_a_bunch, bbu_hom_voltage_calc, bbu_remove_head_bunch, &
    bbu_setup, bbu_track_a_stage, bbu_track_all, check_rf_freq, hom_voltage, logical_to_python, &
    rf_cav_names, write_bunch_by_bunch_info

use count_lines_in_file_mod, only: count_lines_in_file

use bsim_interface, only: insert_phase_trombone, set_tune_3d


use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

subroutine fortran_bbu_add_a_bunch (lat, bbu_beam, bbu_param, beam_init) bind(c)

  use bmad_struct, only: beam_init_struct, lat_struct
  use bbu_track_mod, only: bbu_beam_struct, bbu_param_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  type(c_ptr), value :: bbu_param  ! 0D_NOT_type
  type(bbu_param_struct), pointer :: f_bbu_param
  type(c_ptr), value :: beam_init  ! 0D_NOT_type
  type(beam_init_struct), pointer :: f_beam_init
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  ! inout: f_bbu_param 0D_NOT_type
  if (.not. c_associated(bbu_param)) return
  call c_f_pointer(bbu_param, f_bbu_param)
  ! inout: f_beam_init 0D_NOT_type
  if (.not. c_associated(beam_init)) return
  call c_f_pointer(beam_init, f_beam_init)
  call bbu_add_a_bunch(f_lat, f_bbu_beam, f_bbu_param, f_beam_init)

end subroutine
subroutine fortran_bbu_hom_voltage_calc (lat, bbu_beam, n_period, ix_stage_last_tracked) &
    bind(c)

  use bmad_struct, only: lat_struct
  use bbu_track_mod, only: bbu_beam_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  type(c_ptr), intent(in), value :: n_period  ! 0D_NOT_integer
  integer(c_int) :: f_n_period
  integer(c_int), pointer :: f_n_period_ptr
  type(c_ptr), intent(in), value :: ix_stage_last_tracked  ! 0D_NOT_integer
  integer(c_int) :: f_ix_stage_last_tracked
  integer(c_int), pointer :: f_ix_stage_last_tracked_ptr
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  ! inout: f_n_period 0D_NOT_integer
  if (c_associated(n_period)) then
    call c_f_pointer(n_period, f_n_period_ptr)
  else
    f_n_period_ptr => null()
  endif
  ! inout: f_ix_stage_last_tracked 0D_NOT_integer
  if (c_associated(ix_stage_last_tracked)) then
    call c_f_pointer(ix_stage_last_tracked, f_ix_stage_last_tracked_ptr)
  else
    f_ix_stage_last_tracked_ptr => null()
  endif
  call bbu_hom_voltage_calc(f_lat, f_bbu_beam, f_n_period_ptr, f_ix_stage_last_tracked_ptr)

  ! inout: f_n_period 0D_NOT_integer
  ! no output conversion for f_n_period
  ! inout: f_ix_stage_last_tracked 0D_NOT_integer
  ! no output conversion for f_ix_stage_last_tracked
end subroutine
subroutine fortran_bbu_remove_head_bunch (bbu_beam) bind(c)

  use bbu_track_mod, only: bbu_beam_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  ! ** End of parameters **
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  call bbu_remove_head_bunch(f_bbu_beam)

end subroutine
subroutine fortran_bbu_setup (lat, dt_bunch, bbu_param, bbu_beam) bind(c)

  use bmad_struct, only: lat_struct
  use bbu_track_mod, only: bbu_beam_struct, bbu_param_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), intent(in), value :: dt_bunch  ! 0D_NOT_real
  real(c_double) :: f_dt_bunch
  real(c_double), pointer :: f_dt_bunch_ptr
  type(c_ptr), value :: bbu_param  ! 0D_NOT_type
  type(bbu_param_struct), pointer :: f_bbu_param
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_dt_bunch 0D_NOT_real
  if (c_associated(dt_bunch)) then
    call c_f_pointer(dt_bunch, f_dt_bunch_ptr)
  else
    f_dt_bunch_ptr => null()
  endif
  ! inout: f_bbu_param 0D_NOT_type
  if (.not. c_associated(bbu_param)) return
  call c_f_pointer(bbu_param, f_bbu_param)
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  call bbu_setup(f_lat, f_dt_bunch_ptr, f_bbu_param, f_bbu_beam)

  ! inout: f_dt_bunch 0D_NOT_real
  ! no output conversion for f_dt_bunch
end subroutine
subroutine fortran_bbu_track_a_stage (lat, bbu_beam, bbu_param, lost, ix_stage_tracked) bind(c)

  use bmad_struct, only: lat_struct
  use bbu_track_mod, only: bbu_beam_struct, bbu_param_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  type(c_ptr), value :: bbu_param  ! 0D_NOT_type
  type(bbu_param_struct), pointer :: f_bbu_param
  type(c_ptr), intent(in), value :: lost  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_lost
  logical, target :: f_lost_native
  logical, pointer :: f_lost_native_ptr
  logical(c_bool), pointer :: f_lost_ptr
  type(c_ptr), intent(in), value :: ix_stage_tracked  ! 0D_NOT_integer
  integer(c_int) :: f_ix_stage_tracked
  integer(c_int), pointer :: f_ix_stage_tracked_ptr
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  ! inout: f_bbu_param 0D_NOT_type
  if (.not. c_associated(bbu_param)) return
  call c_f_pointer(bbu_param, f_bbu_param)
  ! inout: f_lost 0D_NOT_logical
  if (c_associated(lost)) then
    call c_f_pointer(lost, f_lost_ptr)
    f_lost_native = f_lost_ptr
    f_lost_native_ptr => f_lost_native
  else
    f_lost_native_ptr => null()
  endif
  ! inout: f_ix_stage_tracked 0D_NOT_integer
  if (c_associated(ix_stage_tracked)) then
    call c_f_pointer(ix_stage_tracked, f_ix_stage_tracked_ptr)
  else
    f_ix_stage_tracked_ptr => null()
  endif
  call bbu_track_a_stage(f_lat, f_bbu_beam, f_bbu_param, f_lost_native_ptr, &
      f_ix_stage_tracked_ptr)

  ! inout: f_lost 0D_NOT_logical
  if (c_associated(lost)) then
    call c_f_pointer(lost, f_lost_ptr)
    f_lost_ptr = f_lost_native
  else
    ! f_lost unset
  endif
  ! inout: f_ix_stage_tracked 0D_NOT_integer
  ! no output conversion for f_ix_stage_tracked
end subroutine
subroutine fortran_bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_normalized, &
    growth_rate, lost, irep) bind(c)

  use bmad_struct, only: beam_init_struct, lat_struct
  use bbu_track_mod, only: bbu_beam_struct, bbu_param_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  type(c_ptr), value :: bbu_param  ! 0D_NOT_type
  type(bbu_param_struct), pointer :: f_bbu_param
  type(c_ptr), value :: beam_init  ! 0D_NOT_type
  type(beam_init_struct), pointer :: f_beam_init
  type(c_ptr), intent(in), value :: hom_voltage_normalized  ! 0D_NOT_real
  real(c_double) :: f_hom_voltage_normalized
  real(c_double), pointer :: f_hom_voltage_normalized_ptr
  type(c_ptr), intent(in), value :: growth_rate  ! 0D_NOT_real
  real(c_double) :: f_growth_rate
  real(c_double), pointer :: f_growth_rate_ptr
  type(c_ptr), intent(in), value :: lost  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_lost
  logical, target :: f_lost_native
  logical, pointer :: f_lost_native_ptr
  logical(c_bool), pointer :: f_lost_ptr
  type(c_ptr), intent(in), value :: irep  ! 0D_NOT_integer
  integer(c_int) :: f_irep
  integer(c_int), pointer :: f_irep_ptr
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  ! inout: f_bbu_param 0D_NOT_type
  if (.not. c_associated(bbu_param)) return
  call c_f_pointer(bbu_param, f_bbu_param)
  ! inout: f_beam_init 0D_NOT_type
  if (.not. c_associated(beam_init)) return
  call c_f_pointer(beam_init, f_beam_init)
  ! inout: f_hom_voltage_normalized 0D_NOT_real
  if (c_associated(hom_voltage_normalized)) then
    call c_f_pointer(hom_voltage_normalized, f_hom_voltage_normalized_ptr)
  else
    f_hom_voltage_normalized_ptr => null()
  endif
  ! inout: f_growth_rate 0D_NOT_real
  if (c_associated(growth_rate)) then
    call c_f_pointer(growth_rate, f_growth_rate_ptr)
  else
    f_growth_rate_ptr => null()
  endif
  ! inout: f_lost 0D_NOT_logical
  if (c_associated(lost)) then
    call c_f_pointer(lost, f_lost_ptr)
    f_lost_native = f_lost_ptr
    f_lost_native_ptr => f_lost_native
  else
    f_lost_native_ptr => null()
  endif
  ! inout: f_irep 0D_NOT_integer
  if (c_associated(irep)) then
    call c_f_pointer(irep, f_irep_ptr)
  else
    f_irep_ptr => null()
  endif
  call bbu_track_all(f_lat, f_bbu_beam, f_bbu_param, f_beam_init, f_hom_voltage_normalized_ptr, &
      f_growth_rate_ptr, f_lost_native_ptr, f_irep_ptr)

  ! inout: f_hom_voltage_normalized 0D_NOT_real
  ! no output conversion for f_hom_voltage_normalized
  ! inout: f_growth_rate 0D_NOT_real
  ! no output conversion for f_growth_rate
  ! inout: f_lost 0D_NOT_logical
  if (c_associated(lost)) then
    call c_f_pointer(lost, f_lost_ptr)
    f_lost_ptr = f_lost_native
  else
    ! f_lost unset
  endif
  ! inout: f_irep 0D_NOT_integer
  ! no output conversion for f_irep
end subroutine
subroutine fortran_check_rf_freq (lat, fb) bind(c)

  use bmad_struct, only: lat_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), intent(in), value :: fb  ! 0D_NOT_real
  real(c_double) :: f_fb
  real(c_double), pointer :: f_fb_ptr
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_fb 0D_NOT_real
  if (c_associated(fb)) then
    call c_f_pointer(fb, f_fb_ptr)
  else
    f_fb_ptr => null()
  endif
  call check_rf_freq(f_lat, f_fb_ptr)

  ! inout: f_fb 0D_NOT_real
  ! no output conversion for f_fb
end subroutine
subroutine fortran_count_lines_in_file (file_name, lines) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: lines  ! 0D_NOT_integer
  INTEGER :: f_lines
  integer(c_int), pointer :: f_lines_ptr
  ! ** End of parameters **
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  call count_lines_in_file(f_file_name, f_lines)

  ! out: f_lines 0D_NOT_integer
  call c_f_pointer(lines, f_lines_ptr)
  f_lines_ptr = f_lines
end subroutine
subroutine fortran_hom_voltage (lr_wake, voltage) bind(c)

  use bmad_struct, only: wake_lr_mode_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: voltage  ! 0D_NOT_real
  real(rp) :: f_voltage
  real(c_double), pointer :: f_voltage_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: lr_wake  ! 0D_NOT_type
  type(wake_lr_mode_struct), pointer :: f_lr_wake
  ! ** End of parameters **
  ! inout: f_lr_wake 0D_NOT_type
  if (.not. c_associated(lr_wake)) return
  call c_f_pointer(lr_wake, f_lr_wake)
  f_voltage = hom_voltage(f_lr_wake)

  ! out: f_voltage 0D_NOT_real
  call c_f_pointer(voltage, f_voltage_ptr)
  f_voltage_ptr = f_voltage
end subroutine
subroutine fortran_insert_phase_trombone (branch) bind(c)

  use bmad_struct, only: branch_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  ! ** End of parameters **
  ! inout: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  call insert_phase_trombone(f_branch)

end subroutine
subroutine fortran_logical_to_python (logic, string) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: logic  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_logic
  logical, target :: f_logic_native
  logical, pointer :: f_logic_native_ptr
  logical(c_bool), pointer :: f_logic_ptr
  ! ** End of parameters **
  ! inout: f_logic 0D_NOT_logical
  if (c_associated(logic)) then
    call c_f_pointer(logic, f_logic_ptr)
    f_logic_native = f_logic_ptr
    f_logic_native_ptr => f_logic_native
  else
    f_logic_native_ptr => null()
  endif
  f_string = logical_to_python(f_logic_native_ptr)

  ! inout: f_logic 0D_NOT_logical
  if (c_associated(logic)) then
    call c_f_pointer(logic, f_logic_ptr)
    f_logic_ptr = f_logic_native
  else
    ! f_logic unset
  endif
  ! out: f_string 0D_NOT_character
  call c_f_pointer(string, f_string_ptr, [len_trim(f_string) + 1]) ! output-only string
  call to_c_str(f_string, f_string_ptr)
end subroutine
subroutine fortran_rf_cav_names (lat) bind(c)

  use bmad_struct, only: lat_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  call rf_cav_names(f_lat)

end subroutine
subroutine fortran_set_tune_3d (branch, target_tunes, mask, use_phase_trombone, z_tune_set, &
    group_knobs, print_err, everything_ok) bind(c)

  use bmad_struct, only: branch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: target_tunes
  real(rp) :: f_target_tunes(3)
  real(c_double), pointer :: f_target_tunes_ptr(:)
  type(c_ptr), intent(in), value :: use_phase_trombone  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_use_phase_trombone
  logical, target :: f_use_phase_trombone_native
  logical, pointer :: f_use_phase_trombone_native_ptr
  logical(c_bool), pointer :: f_use_phase_trombone_ptr
  type(c_ptr), intent(in), value :: z_tune_set  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_z_tune_set
  logical, target :: f_z_tune_set_native
  logical, pointer :: f_z_tune_set_native_ptr
  logical(c_bool), pointer :: f_z_tune_set_ptr
  type(c_ptr), intent(in), value :: group_knobs
  character(len=4096), target :: f_group_knobs
  character(kind=c_char), pointer :: f_group_knobs_ptr(:)
  character(len=4096), pointer :: f_group_knobs_call_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: everything_ok  ! 0D_NOT_logical
  logical :: f_everything_ok
  logical(c_bool), pointer :: f_everything_ok_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  type(c_ptr), intent(in), value :: mask
  character(len=4096), target :: f_mask
  character(kind=c_char), pointer :: f_mask_ptr(:)
  character(len=4096), pointer :: f_mask_call_ptr
  ! ** End of parameters **
  ! inout: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  !! general array (1D_NOT_real)
  if (c_associated(target_tunes)) then
    call c_f_pointer(target_tunes, f_target_tunes_ptr, [3])
    f_target_tunes = f_target_tunes_ptr(:)
  else
    f_target_tunes_ptr => null()
  endif
  ! inout: f_mask 0D_NOT_character
  if (c_associated(mask)) then
    call c_f_pointer(mask, f_mask_ptr, [huge(0)])
    call to_f_str(f_mask_ptr, f_mask)
    f_mask_call_ptr => f_mask
  else
    f_mask_call_ptr => null()
  endif
  ! in: f_use_phase_trombone 0D_NOT_logical
  if (c_associated(use_phase_trombone)) then
    call c_f_pointer(use_phase_trombone, f_use_phase_trombone_ptr)
    f_use_phase_trombone_native = f_use_phase_trombone_ptr
    f_use_phase_trombone_native_ptr => f_use_phase_trombone_native
  else
    f_use_phase_trombone_native_ptr => null()
  endif
  ! in: f_z_tune_set 0D_NOT_logical
  if (c_associated(z_tune_set)) then
    call c_f_pointer(z_tune_set, f_z_tune_set_ptr)
    f_z_tune_set_native = f_z_tune_set_ptr
    f_z_tune_set_native_ptr => f_z_tune_set_native
  else
    f_z_tune_set_native_ptr => null()
  endif
  ! in: f_group_knobs 1D_NOT_character
  if (c_associated(group_knobs)) then
    call c_f_pointer(group_knobs, f_group_knobs_ptr, [huge(0)])
    call to_f_str(f_group_knobs_ptr, f_group_knobs)
    f_group_knobs_call_ptr => f_group_knobs
  else
    f_group_knobs_call_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  f_everything_ok = set_tune_3d(f_branch, f_target_tunes, f_mask_call_ptr, &
      f_use_phase_trombone_native_ptr, f_z_tune_set_native_ptr, f_group_knobs_call_ptr, &
      f_print_err_native_ptr)

  ! inout: f_mask 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_everything_ok 0D_NOT_logical
  call c_f_pointer(everything_ok, f_everything_ok_ptr)
  f_everything_ok_ptr = f_everything_ok
end subroutine
subroutine fortran_write_bunch_by_bunch_info (lat, bbu_beam, bbu_param, this_stage) bind(c)

  use bmad_struct, only: lat_struct
  use bbu_track_mod, only: bbu_beam_struct, bbu_param_struct, bbu_stage_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), value :: bbu_beam  ! 0D_NOT_type
  type(bbu_beam_struct), pointer :: f_bbu_beam
  type(c_ptr), value :: bbu_param  ! 0D_NOT_type
  type(bbu_param_struct), pointer :: f_bbu_param
  type(c_ptr), value :: this_stage  ! 0D_PTR_type
  type(bbu_stage_struct), pointer :: f_this_stage
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_bbu_beam 0D_NOT_type
  if (.not. c_associated(bbu_beam)) return
  call c_f_pointer(bbu_beam, f_bbu_beam)
  ! inout: f_bbu_param 0D_NOT_type
  if (.not. c_associated(bbu_param)) return
  call c_f_pointer(bbu_param, f_bbu_param)
  ! inout: f_this_stage 0D_PTR_type
  if (.not. c_associated(this_stage)) return
  call c_f_pointer(this_stage, f_this_stage)
  call write_bunch_by_bunch_info(f_lat, f_bbu_beam, f_bbu_param, f_this_stage)

end subroutine

end module cppbmad_bsim_routines
