module cppbmad_tao_routines

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

use tao_data_and_eval_mod, only: integrate_max, integrate_min, tao_datum_integrate, &
    tao_datum_s_position, tao_do_wire_scan, tao_ele_geometry_with_misalignments, &
    tao_eval_floor_orbit, tao_evaluate_datum_at_s, tao_evaluate_lat_or_beam_data, &
    tao_expression_hash_substitute, tao_get_data, tao_load_this_datum, &
    tao_pointer_to_datum_ele, tao_scratch_values_calc, tao_to_int, &
    tao_to_phase_and_coupling_reading, tao_tracking_ele_index

use tao_c_interface_mod, only: re_allocate_c_double, tao_c_out_io_buffer_reset

use tao_interface, only: tao_abort_command_file, tao_alias_cmd, tao_beam_emit_calc, &
    tao_beam_track_endpoint, tao_branch_index, tao_chrom_calc_needed, tao_clear_cmd, &
    tao_clip_cmd, tao_close_command_file, tao_command, tao_constraint_type_name, &
    tao_control_tree_list, tao_count_strings, tao_curve_ele_ref, tao_curve_ix_uni, &
    tao_curve_name, tao_curve_rms_calc, tao_d2_d1_name, tao_data_check, tao_data_coupling_init, &
    tao_data_sanity_check, tao_datum_has_associated_ele, tao_datum_name, tao_de_optimizer, &
    tao_ele_shape_info, tao_evaluate_a_datum, tao_evaluate_tune, tao_find_plot_region, &
    tao_fixer, tao_floor_to_screen, tao_floor_to_screen_coords, tao_get_opt_vars, &
    tao_graph_name, tao_init, tao_init_find_elements, tao_init_lattice, tao_init_plotting, &
    tao_is_valid_name, tao_json_cmd, tao_key_info_to_str, tao_lat_bookkeeper, &
    tao_lat_emit_calc, tao_lat_sigma_calc_needed, tao_lattice_calc, tao_limit_calc, &
    tao_lmdif_optimizer, tao_locate_all_elements, tao_locate_elements, tao_mark_lattice_ele, &
    tao_merit, tao_one_turn_map_calc_needed, tao_open_file, tao_open_scratch_file, &
    tao_optimization_status, tao_oreint_building_wall_pt, tao_param_value_at_s, &
    tao_parse_command_args, tao_parse_element_param_str, tao_pause_cmd, tao_pick_universe, &
    tao_pipe_cmd, tao_place_cmd, tao_plot_cmd, tao_plot_setup, tao_plot_struct_transfer, &
    tao_pointer_to_building_wall_shape, tao_pointer_to_datum, tao_pointer_to_ele_shape, &
    tao_pointer_to_tao_lat, tao_pointer_to_universe, tao_pointer_to_universes, &
    tao_print_command_line_info, tao_ptc_normal_form, tao_python_cmd, tao_quiet_set, &
    tao_rad_int_calc_needed, tao_read_cmd, tao_read_phase_space_index, tao_regression_test, &
    tao_remove_blank_characters, tao_run_cmd, tao_scale_ping_data, tao_set_data_useit_opt, &
    tao_set_invalid, tao_set_opt_vars, tao_set_var_model_value, tao_set_var_useit_opt, &
    tao_setup_key_table, tao_shape_init, tao_show_cmd, tao_single_mode, &
    tao_spin_matrices_calc_needed, tao_spin_tracking_turn_on, tao_split_component, &
    tao_srdt_calc_needed, tao_subin_uni_number, tao_symbol_import_from_lat, tao_taper_cmd, &
    tao_to_real, tao_top_level, tao_turn_on_special_calcs_if_needed_for_plotting, &
    tao_uni_atsign_index, tao_universe_index, tao_use_data, tao_use_var, &
    tao_user_is_terminating_optimization, tao_var1_name, tao_var_attrib_name, tao_var_check, &
    tao_var_repoint, tao_var_target_calc, tao_var_useit_plot_calc, tao_write_cmd, &
    tao_x_axis_cmd

use tao_init_data_mod, only: tao_add_to_normal_mode_h_array, tao_allocate_data_array, &
    tao_d2_data_stuffit, tao_init_data, tao_init_data_end_stuff, tao_init_data_in_universe

use tao_init_variables_mod, only: tao_allocate_v1_var, tao_allocate_var_array, &
    tao_init_variables, tao_pointer_to_var_in_lattice, tao_pointer_to_var_in_lattice2

use tao_lattice_calc_mod, only: tao_beam_track, tao_inject_beam, tao_inject_particle, &
    tao_lat_sigma_track, tao_single_track, tao_too_many_particles_lost

use tao_graph_setup_mod, only: tao_calc_data_at_s_pts, tao_curve_beam_ellipse_setup, &
    tao_curve_check_universe, tao_curve_data_setup, tao_curve_datum_calc, &
    tao_data_type_substitute, tao_data_useit_plot_calc, tao_graph_controller_setup, &
    tao_graph_data_setup, tao_graph_data_slice_setup, tao_graph_dynamic_aperture_setup, &
    tao_graph_histogram_setup, tao_graph_phase_space_setup, tao_graph_s_min_max_calc, &
    tao_graph_setup, tao_particle_data_value, tao_phase_space_axis_index, tao_set_curve_invalid

use tao_wave_mod, only: tao_cbar_wave_anal, tao_orbit_beta_wave_anal, tao_phase_wave_anal, &
    tao_wave_analysis, tao_wave_cmd, tao_wave_fit

use tao_change_mod, only: tao_change_ele, tao_change_tune, tao_change_var, tao_change_z_tune, &
    tao_to_change_number

use tao_command_mod, only: tao_cmd_history_record, tao_next_word, tao_re_execute

use tao_plot_window_mod, only: tao_create_plot_window, tao_destroy_plot_window

use tao_struct, only: tao_deallocate_plot_cache, &
    tao_lattice_branches_equal_tao_lattice_branches, tao_lattice_equal_tao_lattice

use tao_dmerit_mod, only: tao_dmerit_calc, tao_dmodel_dvar_calc, tao_veto_vars_with_zero_dmodel

use tao_plot_mod, only: tao_draw_beam_chamber_wall, tao_draw_curve_data, &
    tao_draw_ele_for_floor_plan, tao_draw_floor_plan, tao_draw_graph_axes, &
    tao_draw_histogram_data, tao_draw_lat_layout, tao_draw_plots, tao_plot_data, &
    tao_plot_histogram, tao_plot_key_table, tao_plot_wave, tao_set_floor_plan_axis_label

use tao_geodesic_lm_optimizer_mod, only: tao_geodesic_lm_optimizer

use tao_get_user_input_mod, only: tao_get_user_input

use tao_init_mod, only: tao_init_beam_in_universe, tao_init_beams, tao_init_dynamic_aperture, &
    tao_init_global

use tao_lm_optimizer_mod, only: tao_lm_optimizer

use tao_scale_mod, only: tao_scale_cmd, tao_scale_graph, tao_scale_plot

use tao_set_mod, only: tao_set_beam_cmd, tao_set_beam_init_cmd, tao_set_bmad_com_cmd, &
    tao_set_branch_cmd, tao_set_calculate_cmd, tao_set_curve_cmd, tao_set_data_cmd, &
    tao_set_default_cmd, tao_set_drawing_cmd, tao_set_dynamic_aperture_cmd, &
    tao_set_elements_cmd, tao_set_geodesic_lm_cmd, tao_set_global_cmd, tao_set_graph_cmd, &
    tao_set_integer_value, tao_set_key_cmd, tao_set_lattice_cmd, tao_set_logical_value, &
    tao_set_openmp_n_threads, tao_set_opti_de_param_cmd, tao_set_particle_start_cmd, &
    tao_set_plot_cmd, tao_set_plot_page_cmd, tao_set_ptc_com_cmd, tao_set_qp_axis_struct, &
    tao_set_qp_point_struct, tao_set_qp_rect_struct, tao_set_ran_state_cmd, tao_set_real_value, &
    tao_set_region_cmd, tao_set_space_charge_com_cmd, tao_set_symbolic_number_cmd, &
    tao_set_tune_cmd, tao_set_universe_cmd, tao_set_var_cmd, tao_set_wave_cmd, &
    tao_set_z_tune_cmd

use tao_top10_mod, only: tao_show_constraints, tao_top10_derivative_print, &
    tao_top10_merit_categories_print, tao_var_write

use tao_svd_optimizer_mod, only: tao_svd_optimizer

use tao_x_scale_mod, only: tao_x_scale_cmd, tao_x_scale_graph, tao_x_scale_plot


use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

subroutine fortran_integrate_max (ix_start, ix_ele, datum_value, ix_m, branch, vec, datum) &
    bind(c)

  use bmad_struct, only: branch_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_start  ! 0D_NOT_integer
  integer(c_int) :: f_ix_start
  integer(c_int), pointer :: f_ix_start_ptr
  type(c_ptr), intent(in), value :: ix_ele  ! 0D_NOT_integer
  integer(c_int) :: f_ix_ele
  integer(c_int), pointer :: f_ix_ele_ptr
  type(c_ptr), intent(in), value :: datum_value  ! 0D_NOT_real
  real(c_double) :: f_datum_value
  real(c_double), pointer :: f_datum_value_ptr
  type(c_ptr), intent(in), value :: ix_m  ! 0D_NOT_integer
  integer(c_int) :: f_ix_m
  integer(c_int), pointer :: f_ix_m_ptr
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  type(c_ptr), intent(in), value :: vec
  type(real_container_alloc), pointer :: f_vec
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** End of parameters **
  ! inout: f_ix_start 0D_NOT_integer
  if (c_associated(ix_start)) then
    call c_f_pointer(ix_start, f_ix_start_ptr)
  else
    f_ix_start_ptr => null()
  endif
  ! inout: f_ix_ele 0D_NOT_integer
  if (c_associated(ix_ele)) then
    call c_f_pointer(ix_ele, f_ix_ele_ptr)
  else
    f_ix_ele_ptr => null()
  endif
  ! inout: f_datum_value 0D_NOT_real
  if (c_associated(datum_value)) then
    call c_f_pointer(datum_value, f_datum_value_ptr)
  else
    f_datum_value_ptr => null()
  endif
  ! inout: f_ix_m 0D_NOT_integer
  if (c_associated(ix_m)) then
    call c_f_pointer(ix_m, f_ix_m_ptr)
  else
    f_ix_m_ptr => null()
  endif
  ! inout: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  !! container general array (1D_ALLOC_real)
  if (c_associated(vec))   call c_f_pointer(vec, f_vec)
  ! inout: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  call integrate_max(f_ix_start_ptr, f_ix_ele_ptr, f_datum_value_ptr, f_ix_m_ptr, f_branch, &
      f_vec%data, f_datum)

  ! inout: f_ix_start 0D_NOT_integer
  ! no output conversion for f_ix_start
  ! inout: f_ix_ele 0D_NOT_integer
  ! no output conversion for f_ix_ele
  ! inout: f_datum_value 0D_NOT_real
  ! no output conversion for f_datum_value
  ! inout: f_ix_m 0D_NOT_integer
  ! no output conversion for f_ix_m
end subroutine
subroutine fortran_integrate_min (ix_start, ix_ele, datum_value, ix_m, branch, vec, datum) &
    bind(c)

  use bmad_struct, only: branch_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_start  ! 0D_NOT_integer
  integer(c_int) :: f_ix_start
  integer(c_int), pointer :: f_ix_start_ptr
  type(c_ptr), intent(in), value :: ix_ele  ! 0D_NOT_integer
  integer(c_int) :: f_ix_ele
  integer(c_int), pointer :: f_ix_ele_ptr
  type(c_ptr), intent(in), value :: datum_value  ! 0D_NOT_real
  real(c_double) :: f_datum_value
  real(c_double), pointer :: f_datum_value_ptr
  type(c_ptr), intent(in), value :: ix_m  ! 0D_NOT_integer
  integer(c_int) :: f_ix_m
  integer(c_int), pointer :: f_ix_m_ptr
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  type(c_ptr), intent(in), value :: vec
  type(real_container_alloc), pointer :: f_vec
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** End of parameters **
  ! inout: f_ix_start 0D_NOT_integer
  if (c_associated(ix_start)) then
    call c_f_pointer(ix_start, f_ix_start_ptr)
  else
    f_ix_start_ptr => null()
  endif
  ! inout: f_ix_ele 0D_NOT_integer
  if (c_associated(ix_ele)) then
    call c_f_pointer(ix_ele, f_ix_ele_ptr)
  else
    f_ix_ele_ptr => null()
  endif
  ! inout: f_datum_value 0D_NOT_real
  if (c_associated(datum_value)) then
    call c_f_pointer(datum_value, f_datum_value_ptr)
  else
    f_datum_value_ptr => null()
  endif
  ! inout: f_ix_m 0D_NOT_integer
  if (c_associated(ix_m)) then
    call c_f_pointer(ix_m, f_ix_m_ptr)
  else
    f_ix_m_ptr => null()
  endif
  ! inout: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  !! container general array (1D_ALLOC_real)
  if (c_associated(vec))   call c_f_pointer(vec, f_vec)
  ! inout: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  call integrate_min(f_ix_start_ptr, f_ix_ele_ptr, f_datum_value_ptr, f_ix_m_ptr, f_branch, &
      f_vec%data, f_datum)

  ! inout: f_ix_start 0D_NOT_integer
  ! no output conversion for f_ix_start
  ! inout: f_ix_ele 0D_NOT_integer
  ! no output conversion for f_ix_ele
  ! inout: f_datum_value 0D_NOT_real
  ! no output conversion for f_datum_value
  ! inout: f_ix_m 0D_NOT_integer
  ! no output conversion for f_ix_m
end subroutine
subroutine fortran_re_allocate_c_double (re, n, exact, init_val) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: n  ! 0D_NOT_integer
  integer :: f_n
  type(c_ptr), intent(in), value :: exact  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact
  logical, target :: f_exact_native
  logical, pointer :: f_exact_native_ptr
  logical(c_bool), pointer :: f_exact_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: re
  type(real_container_alloc), pointer :: f_re
  type(c_ptr), intent(in), value :: init_val  ! 0D_NOT_real
  real(c_double) :: f_init_val
  real(c_double), pointer :: f_init_val_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (c_associated(re))   call c_f_pointer(re, f_re)
  ! in: f_n 0D_NOT_integer
  f_n = n
  ! in: f_exact 0D_NOT_logical
  if (c_associated(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_native = f_exact_ptr
    f_exact_native_ptr => f_exact_native
  else
    f_exact_native_ptr => null()
  endif
  ! inout: f_init_val 0D_NOT_real
  if (c_associated(init_val)) then
    call c_f_pointer(init_val, f_init_val_ptr)
  else
    f_init_val_ptr => null()
  endif
  call re_allocate_c_double(f_re%data, f_n, f_exact_native_ptr, f_init_val_ptr)

  ! inout: f_init_val 0D_NOT_real
  ! no output conversion for f_init_val
end subroutine
subroutine fortran_tao_abort_command_file (force_abort) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: force_abort  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_force_abort
  logical, target :: f_force_abort_native
  logical, pointer :: f_force_abort_native_ptr
  logical(c_bool), pointer :: f_force_abort_ptr
  ! ** End of parameters **
  ! in: f_force_abort 0D_NOT_logical
  if (c_associated(force_abort)) then
    call c_f_pointer(force_abort, f_force_abort_ptr)
    f_force_abort_native = f_force_abort_ptr
    f_force_abort_native_ptr => f_force_abort_native
  else
    f_force_abort_native_ptr => null()
  endif
  call tao_abort_command_file(f_force_abort_native_ptr)

end subroutine
subroutine fortran_tao_add_to_normal_mode_h_array (h_str, h_array) bind(c)

  use bmad_struct, only: resonance_h_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: h_str
  character(len=4096), target :: f_h_str
  character(kind=c_char), pointer :: f_h_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: h_array
  type(resonance_h_struct_container_alloc), pointer :: f_h_array
  ! ** End of parameters **
  ! in: f_h_str 0D_NOT_character
  if (.not. c_associated(h_str)) return
  call c_f_pointer(h_str, f_h_str_ptr, [huge(0)])
  call to_f_str(f_h_str_ptr, f_h_str)
  !! container type array (1D_ALLOC_type)
  if (c_associated(h_array))   call c_f_pointer(h_array, f_h_array)
  call tao_add_to_normal_mode_h_array(f_h_str, f_h_array%data)

end subroutine
subroutine fortran_tao_alias_cmd (alias, string) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: alias
  character(len=4096), target :: f_alias
  character(kind=c_char), pointer :: f_alias_ptr(:)
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** End of parameters **
  ! in: f_alias 0D_NOT_character
  if (.not. c_associated(alias)) return
  call c_f_pointer(alias, f_alias_ptr, [huge(0)])
  call to_f_str(f_alias_ptr, f_alias)
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  call tao_alias_cmd(f_alias, f_string)

end subroutine
subroutine fortran_tao_allocate_data_array (u, n_data, exact) bind(c)

  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), intent(in), value :: n_data  ! 0D_NOT_integer
  integer(c_int) :: f_n_data
  integer(c_int), pointer :: f_n_data_ptr
  type(c_ptr), intent(in), value :: exact  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact
  logical, target :: f_exact_native
  logical, pointer :: f_exact_native_ptr
  logical(c_bool), pointer :: f_exact_ptr
  ! ** End of parameters **
  ! inout: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! inout: f_n_data 0D_NOT_integer
  if (c_associated(n_data)) then
    call c_f_pointer(n_data, f_n_data_ptr)
  else
    f_n_data_ptr => null()
  endif
  ! inout: f_exact 0D_NOT_logical
  if (c_associated(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_native = f_exact_ptr
    f_exact_native_ptr => f_exact_native
  else
    f_exact_native_ptr => null()
  endif
  call tao_allocate_data_array(f_u, f_n_data_ptr, f_exact_native_ptr)

  ! inout: f_n_data 0D_NOT_integer
  ! no output conversion for f_n_data
  ! inout: f_exact 0D_NOT_logical
  if (c_associated(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_ptr = f_exact_native
  else
    ! f_exact unset
  endif
end subroutine
subroutine fortran_tao_allocate_v1_var (n_v1, save_old) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: n_v1  ! 0D_NOT_integer
  integer(c_int) :: f_n_v1
  integer(c_int), pointer :: f_n_v1_ptr
  type(c_ptr), intent(in), value :: save_old  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_save_old
  logical, target :: f_save_old_native
  logical, pointer :: f_save_old_native_ptr
  logical(c_bool), pointer :: f_save_old_ptr
  ! ** End of parameters **
  ! inout: f_n_v1 0D_NOT_integer
  if (c_associated(n_v1)) then
    call c_f_pointer(n_v1, f_n_v1_ptr)
  else
    f_n_v1_ptr => null()
  endif
  ! inout: f_save_old 0D_NOT_logical
  if (c_associated(save_old)) then
    call c_f_pointer(save_old, f_save_old_ptr)
    f_save_old_native = f_save_old_ptr
    f_save_old_native_ptr => f_save_old_native
  else
    f_save_old_native_ptr => null()
  endif
  call tao_allocate_v1_var(f_n_v1_ptr, f_save_old_native_ptr)

  ! inout: f_n_v1 0D_NOT_integer
  ! no output conversion for f_n_v1
  ! inout: f_save_old 0D_NOT_logical
  if (c_associated(save_old)) then
    call c_f_pointer(save_old, f_save_old_ptr)
    f_save_old_ptr = f_save_old_native
  else
    ! f_save_old unset
  endif
end subroutine
subroutine fortran_tao_allocate_var_array (n_var, default_good_user) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: n_var  ! 0D_NOT_integer
  integer :: f_n_var
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: default_good_user  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_default_good_user
  logical, target :: f_default_good_user_native
  logical, pointer :: f_default_good_user_native_ptr
  logical(c_bool), pointer :: f_default_good_user_ptr
  ! ** End of parameters **
  ! in: f_n_var 0D_NOT_integer
  f_n_var = n_var
  ! inout: f_default_good_user 0D_NOT_logical
  if (c_associated(default_good_user)) then
    call c_f_pointer(default_good_user, f_default_good_user_ptr)
    f_default_good_user_native = f_default_good_user_ptr
    f_default_good_user_native_ptr => f_default_good_user_native
  else
    f_default_good_user_native_ptr => null()
  endif
  call tao_allocate_var_array(f_n_var, f_default_good_user_native_ptr)

  ! inout: f_default_good_user 0D_NOT_logical
  if (c_associated(default_good_user)) then
    call c_f_pointer(default_good_user, f_default_good_user_ptr)
    f_default_good_user_ptr = f_default_good_user_native
  else
    ! f_default_good_user unset
  endif
end subroutine
subroutine fortran_tao_beam_emit_calc (plane, emit_type, ele, bunch_params, emit) bind(c)

  use bmad_struct, only: bunch_params_struct, ele_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: plane  ! 0D_NOT_integer
  integer :: f_plane
  integer(c_int) :: emit_type  ! 0D_NOT_integer
  integer :: f_emit_type
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: bunch_params  ! 0D_NOT_type
  type(bunch_params_struct), pointer :: f_bunch_params
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: emit  ! 0D_NOT_real
  real(rp) :: f_emit
  real(c_double), pointer :: f_emit_ptr
  ! ** End of parameters **
  ! in: f_plane 0D_NOT_integer
  f_plane = plane
  ! in: f_emit_type 0D_NOT_integer
  f_emit_type = emit_type
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_bunch_params 0D_NOT_type
  if (.not. c_associated(bunch_params)) return
  call c_f_pointer(bunch_params, f_bunch_params)
  f_emit = tao_beam_emit_calc(f_plane, f_emit_type, f_ele, f_bunch_params)

  ! out: f_emit 0D_NOT_real
  call c_f_pointer(emit, f_emit_ptr)
  f_emit_ptr = f_emit
end subroutine
subroutine fortran_tao_beam_track (u, tao_lat, ix_branch, beam, calc_ok) bind(c)

  use tao_struct, only: tao_lattice_struct, tao_universe_struct
  use bmad_struct, only: beam_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  integer(c_int) :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  logical(c_bool), pointer :: f_calc_ok_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: beam  ! 0D_NOT_type
  type(beam_struct), pointer :: f_beam
  ! ** End of parameters **
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_ix_branch 0D_NOT_integer
  f_ix_branch = ix_branch
  ! inout: f_beam 0D_NOT_type
  if (.not. c_associated(beam)) return
  call c_f_pointer(beam, f_beam)
  call tao_beam_track(f_u, f_tao_lat, f_ix_branch, f_beam, f_calc_ok)

  ! out: f_calc_ok 0D_NOT_logical
  call c_f_pointer(calc_ok, f_calc_ok_ptr)
  f_calc_ok_ptr = f_calc_ok
end subroutine
subroutine fortran_tao_beam_track_endpoint (ele_id, lat, branch_str, where, u, ele) bind(c)

  use bmad_struct, only: ele_struct, lat_struct
  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: ele_id
  character(len=4096), target :: f_ele_id
  character(kind=c_char), pointer :: f_ele_id_ptr(:)
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  ! ** Out parameters **
  type(c_ptr), value :: ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele
  ! ** End of parameters **
  ! in: f_ele_id 0D_NOT_character
  if (.not. c_associated(ele_id)) return
  call c_f_pointer(ele_id, f_ele_id_ptr, [huge(0)])
  call to_f_str(f_ele_id_ptr, f_ele_id)
  ! in: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) return
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  f_ele = tao_beam_track_endpoint(f_ele_id, f_lat, f_branch_str, f_where, f_u)

  ! out: f_ele 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_branch_index (ix_branch, ix_this) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_this  ! 0D_NOT_integer
  integer :: f_ix_this
  integer(c_int), pointer :: f_ix_this_ptr
  ! ** End of parameters **
  ! in: f_ix_branch 0D_NOT_integer
  f_ix_branch = ix_branch
  f_ix_this = tao_branch_index(f_ix_branch)

  ! out: f_ix_this 0D_NOT_integer
  call c_f_pointer(ix_this, f_ix_this_ptr)
  f_ix_this_ptr = f_ix_this
end subroutine
subroutine fortran_tao_c_out_io_buffer_reset () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_c_out_io_buffer_reset()

end subroutine
subroutine fortran_tao_calc_data_at_s_pts (tao_lat, curve, comp_sign, good) bind(c)

  use tao_struct, only: tao_curve_struct, tao_lattice_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  type(c_ptr), intent(in), value :: comp_sign  ! 0D_NOT_real
  real(c_double) :: f_comp_sign
  real(c_double), pointer :: f_comp_sign_ptr
  type(c_ptr), intent(in), value :: good
  type(logical_container_alloc), pointer :: f_good
  ! ** End of parameters **
  ! inout: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! inout: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! inout: f_comp_sign 0D_NOT_real
  if (c_associated(comp_sign)) then
    call c_f_pointer(comp_sign, f_comp_sign_ptr)
  else
    f_comp_sign_ptr => null()
  endif
  !! container general array (1D_ALLOC_logical)
  if (c_associated(good))   call c_f_pointer(good, f_good)
  call tao_calc_data_at_s_pts(f_tao_lat, f_curve, f_comp_sign_ptr, f_good%data)

  ! inout: f_comp_sign 0D_NOT_real
  ! no output conversion for f_comp_sign
end subroutine
subroutine fortran_tao_cbar_wave_anal (plot) bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  call tao_cbar_wave_anal(f_plot)

end subroutine
subroutine fortran_tao_change_ele (ele_name, attrib_name, num_str, update, err_flag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: ele_name
  character(len=4096), target :: f_ele_name
  character(kind=c_char), pointer :: f_ele_name_ptr(:)
  type(c_ptr), intent(in), value :: attrib_name
  character(len=4096), target :: f_attrib_name
  character(kind=c_char), pointer :: f_attrib_name_ptr(:)
  type(c_ptr), intent(in), value :: num_str
  character(len=4096), target :: f_num_str
  character(kind=c_char), pointer :: f_num_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: update  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_update
  logical, target :: f_update_native
  logical, pointer :: f_update_native_ptr
  logical(c_bool), pointer :: f_update_ptr
  ! ** End of parameters **
  ! in: f_ele_name 0D_NOT_character
  if (.not. c_associated(ele_name)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(ele_name, f_ele_name_ptr, [huge(0)])
  call to_f_str(f_ele_name_ptr, f_ele_name)
  ! in: f_attrib_name 0D_NOT_character
  if (.not. c_associated(attrib_name)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(attrib_name, f_attrib_name_ptr, [huge(0)])
  call to_f_str(f_attrib_name_ptr, f_attrib_name)
  ! in: f_num_str 0D_NOT_character
  if (.not. c_associated(num_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(num_str, f_num_str_ptr, [huge(0)])
  call to_f_str(f_num_str_ptr, f_num_str)
  ! inout: f_update 0D_NOT_logical
  if (c_associated(update)) then
    call c_f_pointer(update, f_update_ptr)
    f_update_native = f_update_ptr
    f_update_native_ptr => f_update_native
  else
    f_update_native_ptr => null()
  endif
  call tao_change_ele(f_ele_name, f_attrib_name, f_num_str, f_update_native_ptr, f_err_flag)

  ! inout: f_update 0D_NOT_logical
  if (c_associated(update)) then
    call c_f_pointer(update, f_update_ptr)
    f_update_ptr = f_update_native
  else
    ! f_update unset
  endif
  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_change_tune (branch_str, mask_str, print_list, dqa_str, dqb_str, &
    err_flag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  type(c_ptr), intent(in), value :: mask_str
  character(len=4096), target :: f_mask_str
  character(kind=c_char), pointer :: f_mask_str_ptr(:)
  logical(c_bool) :: print_list  ! 0D_NOT_logical
  logical :: f_print_list
  type(c_ptr), intent(in), value :: dqa_str
  character(len=4096), target :: f_dqa_str
  character(kind=c_char), pointer :: f_dqa_str_ptr(:)
  type(c_ptr), intent(in), value :: dqb_str
  character(len=4096), target :: f_dqb_str
  character(kind=c_char), pointer :: f_dqb_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  ! in: f_mask_str 0D_NOT_character
  if (.not. c_associated(mask_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(mask_str, f_mask_str_ptr, [huge(0)])
  call to_f_str(f_mask_str_ptr, f_mask_str)
  ! in: f_print_list 0D_NOT_logical
  f_print_list = print_list
  ! in: f_dqa_str 0D_NOT_character
  if (.not. c_associated(dqa_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(dqa_str, f_dqa_str_ptr, [huge(0)])
  call to_f_str(f_dqa_str_ptr, f_dqa_str)
  ! in: f_dqb_str 0D_NOT_character
  if (.not. c_associated(dqb_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(dqb_str, f_dqb_str_ptr, [huge(0)])
  call to_f_str(f_dqb_str_ptr, f_dqb_str)
  call tao_change_tune(f_branch_str, f_mask_str, f_print_list, f_dqa_str, f_dqb_str, &
      f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_change_var (name, num_str, silent, err_flag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  type(c_ptr), intent(in), value :: num_str
  character(len=4096), target :: f_num_str
  character(kind=c_char), pointer :: f_num_str_ptr(:)
  logical(c_bool) :: silent  ! 0D_NOT_logical
  logical :: f_silent
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_num_str 0D_NOT_character
  if (.not. c_associated(num_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(num_str, f_num_str_ptr, [huge(0)])
  call to_f_str(f_num_str_ptr, f_num_str)
  ! in: f_silent 0D_NOT_logical
  f_silent = silent
  call tao_change_var(f_name, f_num_str, f_silent, f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_change_z_tune (branch_str, dq_str, err_flag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  type(c_ptr), intent(in), value :: dq_str
  character(len=4096), target :: f_dq_str
  character(kind=c_char), pointer :: f_dq_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  ! in: f_dq_str 0D_NOT_character
  if (.not. c_associated(dq_str)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(dq_str, f_dq_str_ptr, [huge(0)])
  call to_f_str(f_dq_str_ptr, f_dq_str)
  call tao_change_z_tune(f_branch_str, f_dq_str, f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_chrom_calc_needed (data_type, data_source, do_chrom) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: do_chrom  ! 0D_NOT_logical
  logical :: f_do_chrom
  logical(c_bool), pointer :: f_do_chrom_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: data_source
  character(len=4096), target :: f_data_source
  character(kind=c_char), pointer :: f_data_source_ptr(:)
  ! ** End of parameters **
  ! inout: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! inout: f_data_source 0D_NOT_character
  if (.not. c_associated(data_source)) return
  call c_f_pointer(data_source, f_data_source_ptr, [huge(0)])
  call to_f_str(f_data_source_ptr, f_data_source)
  f_do_chrom = tao_chrom_calc_needed(f_data_type, f_data_source)

  ! inout: f_data_type 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_data_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_do_chrom 0D_NOT_logical
  call c_f_pointer(do_chrom, f_do_chrom_ptr)
  f_do_chrom_ptr = f_do_chrom
end subroutine
subroutine fortran_tao_clear_cmd (cmd_line) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: cmd_line
  character(len=4096), target :: f_cmd_line
  character(kind=c_char), pointer :: f_cmd_line_ptr(:)
  ! ** End of parameters **
  ! in: f_cmd_line 0D_NOT_character
  if (.not. c_associated(cmd_line)) return
  call c_f_pointer(cmd_line, f_cmd_line_ptr, [huge(0)])
  call to_f_str(f_cmd_line_ptr, f_cmd_line)
  call tao_clear_cmd(f_cmd_line)

end subroutine
subroutine fortran_tao_clip_cmd (gang, where, value1, value2) bind(c)

  implicit none
  ! ** In parameters **
  logical(c_bool) :: gang  ! 0D_NOT_logical
  logical :: f_gang
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: value1  ! 0D_NOT_real
  real(c_double) :: f_value1
  real(c_double), pointer :: f_value1_ptr
  type(c_ptr), intent(in), value :: value2  ! 0D_NOT_real
  real(c_double) :: f_value2
  real(c_double), pointer :: f_value2_ptr
  ! ** End of parameters **
  ! in: f_gang 0D_NOT_logical
  f_gang = gang
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! inout: f_value1 0D_NOT_real
  if (c_associated(value1)) then
    call c_f_pointer(value1, f_value1_ptr)
  else
    f_value1_ptr => null()
  endif
  ! inout: f_value2 0D_NOT_real
  if (c_associated(value2)) then
    call c_f_pointer(value2, f_value2_ptr)
  else
    f_value2_ptr => null()
  endif
  call tao_clip_cmd(f_gang, f_where, f_value1_ptr, f_value2_ptr)

  ! inout: f_value1 0D_NOT_real
  ! no output conversion for f_value1
  ! inout: f_value2 0D_NOT_real
  ! no output conversion for f_value2
end subroutine
subroutine fortran_tao_close_command_file () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_close_command_file()

end subroutine
subroutine fortran_tao_cmd_history_record (cmd) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: cmd
  character(len=4096), target :: f_cmd
  character(kind=c_char), pointer :: f_cmd_ptr(:)
  ! ** End of parameters **
  ! inout: f_cmd 0D_NOT_character
  if (.not. c_associated(cmd)) return
  call c_f_pointer(cmd, f_cmd_ptr, [huge(0)])
  call to_f_str(f_cmd_ptr, f_cmd)
  call tao_cmd_history_record(f_cmd)

  ! inout: f_cmd 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_command (command_line, err, err_is_fatal) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: command_line
  character(len=4096), target :: f_command_line
  character(kind=c_char), pointer :: f_command_line_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_is_fatal  ! 0D_NOT_logical
  logical :: f_err_is_fatal
  logical(c_bool), pointer :: f_err_is_fatal_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err
  logical, target :: f_err_native
  logical, pointer :: f_err_native_ptr
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_command_line 0D_NOT_character
  if (.not. c_associated(command_line)) return
  call c_f_pointer(command_line, f_command_line_ptr, [huge(0)])
  call to_f_str(f_command_line_ptr, f_command_line)
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_native = f_err_ptr
    f_err_native_ptr => f_err_native
  else
    f_err_native_ptr => null()
  endif
  call tao_command(f_command_line, f_err_native_ptr, f_err_is_fatal)

  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_ptr = f_err_native
  else
    ! f_err unset
  endif
  ! out: f_err_is_fatal 0D_NOT_logical
  call c_f_pointer(err_is_fatal, f_err_is_fatal_ptr)
  f_err_is_fatal_ptr = f_err_is_fatal
end subroutine
subroutine fortran_tao_constraint_type_name (datum, datum_name) bind(c)

  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: datum_name
  character(len=4096), target :: f_datum_name
  character(kind=c_char), pointer :: f_datum_name_ptr(:)
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  f_datum_name = tao_constraint_type_name(f_datum)

  ! out: f_datum_name 0D_NOT_character
  call c_f_pointer(datum_name, f_datum_name_ptr, [len_trim(f_datum_name) + 1]) ! output-only string
  call to_c_str(f_datum_name, f_datum_name_ptr)
end subroutine
subroutine fortran_tao_control_tree_list (ele, tree) bind(c)

  use bmad_struct, only: ele_pointer_struct, ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), intent(in), value :: tree
  type(ele_pointer_struct_container_alloc), pointer :: f_tree
  ! ** End of parameters **
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  !! container type array (1D_ALLOC_type)
  if (c_associated(tree))   call c_f_pointer(tree, f_tree)
  call tao_control_tree_list(f_ele, f_tree%data)

end subroutine
subroutine fortran_tao_count_strings (string, pattern, num) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: pattern
  character(len=4096), target :: f_pattern
  character(kind=c_char), pointer :: f_pattern_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: num  ! 0D_NOT_integer
  integer :: f_num
  integer(c_int), pointer :: f_num_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! in: f_pattern 0D_NOT_character
  if (.not. c_associated(pattern)) return
  call c_f_pointer(pattern, f_pattern_ptr, [huge(0)])
  call to_f_str(f_pattern_ptr, f_pattern)
  call tao_count_strings(f_string, f_pattern, f_num)

  ! out: f_num 0D_NOT_integer
  call c_f_pointer(num, f_num_ptr)
  f_num_ptr = f_num
end subroutine
subroutine fortran_tao_create_plot_window () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_create_plot_window()

end subroutine
subroutine fortran_tao_curve_beam_ellipse_setup (curve) bind(c)

  use tao_struct, only: tao_curve_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** End of parameters **
  ! inout: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  call tao_curve_beam_ellipse_setup(f_curve)

end subroutine
subroutine fortran_tao_curve_check_universe (curve, uni, is_ok) bind(c)

  use tao_struct, only: tao_curve_struct, tao_universe_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: uni  ! 0D_PTR_type
  type(tao_universe_struct), pointer :: f_uni
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_ok  ! 0D_NOT_logical
  logical :: f_is_ok
  logical(c_bool), pointer :: f_is_ok_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** End of parameters **
  ! inout: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_uni 0D_PTR_type
  if (.not. c_associated(uni)) return
  call c_f_pointer(uni, f_uni)
  f_is_ok = tao_curve_check_universe(f_curve, f_uni)

  ! out: f_is_ok 0D_NOT_logical
  call c_f_pointer(is_ok, f_is_ok_ptr)
  f_is_ok_ptr = f_is_ok
end subroutine
subroutine fortran_tao_curve_data_setup (plot, graph, curve) bind(c)

  use tao_struct, only: tao_curve_struct, tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! inout: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  call tao_curve_data_setup(f_plot, f_graph, f_curve)

end subroutine
subroutine fortran_tao_curve_datum_calc (eles, plot, curve, who) bind(c)

  use bmad_struct, only: ele_pointer_struct
  use tao_struct, only: tao_curve_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: eles
  type(ele_pointer_struct_container_alloc), pointer :: f_eles
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(eles))   call c_f_pointer(eles, f_eles)
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  call tao_curve_datum_calc(f_eles%data, f_plot, f_curve, f_who)

end subroutine
subroutine fortran_tao_curve_ele_ref (curve, point_to_ele_ref, ele_track) bind(c)

  use tao_struct, only: tao_curve_struct
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** Out parameters **
  type(c_ptr), value :: ele_track  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele_track
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: point_to_ele_ref  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_point_to_ele_ref
  logical, target :: f_point_to_ele_ref_native
  logical, pointer :: f_point_to_ele_ref_native_ptr
  logical(c_bool), pointer :: f_point_to_ele_ref_ptr
  ! ** End of parameters **
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! inout: f_point_to_ele_ref 0D_NOT_logical
  if (c_associated(point_to_ele_ref)) then
    call c_f_pointer(point_to_ele_ref, f_point_to_ele_ref_ptr)
    f_point_to_ele_ref_native = f_point_to_ele_ref_ptr
    f_point_to_ele_ref_native_ptr => f_point_to_ele_ref_native
  else
    f_point_to_ele_ref_native_ptr => null()
  endif
  f_ele_track = tao_curve_ele_ref(f_curve, f_point_to_ele_ref_native_ptr)

  ! inout: f_point_to_ele_ref 0D_NOT_logical
  if (c_associated(point_to_ele_ref)) then
    call c_f_pointer(point_to_ele_ref, f_point_to_ele_ref_ptr)
    f_point_to_ele_ref_ptr = f_point_to_ele_ref_native
  else
    ! f_point_to_ele_ref unset
  endif
  ! out: f_ele_track 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_curve_ix_uni (curve, ix_uni) bind(c)

  use tao_struct, only: tao_curve_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  integer(c_int), pointer :: f_ix_uni_ptr
  ! ** End of parameters **
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  f_ix_uni = tao_curve_ix_uni(f_curve)

  ! out: f_ix_uni 0D_NOT_integer
  call c_f_pointer(ix_uni, f_ix_uni_ptr)
  f_ix_uni_ptr = f_ix_uni
end subroutine
subroutine fortran_tao_curve_name (curve, use_region, curve_name) bind(c)

  use tao_struct, only: tao_curve_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  type(c_ptr), intent(in), value :: use_region  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_use_region
  logical, target :: f_use_region_native
  logical, pointer :: f_use_region_native_ptr
  logical(c_bool), pointer :: f_use_region_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: curve_name
  character(len=4096), target :: f_curve_name
  character(kind=c_char), pointer :: f_curve_name_ptr(:)
  ! ** End of parameters **
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_use_region 0D_NOT_logical
  if (c_associated(use_region)) then
    call c_f_pointer(use_region, f_use_region_ptr)
    f_use_region_native = f_use_region_ptr
    f_use_region_native_ptr => f_use_region_native
  else
    f_use_region_native_ptr => null()
  endif
  f_curve_name = tao_curve_name(f_curve, f_use_region_native_ptr)

  ! out: f_curve_name 0D_NOT_character
  call c_f_pointer(curve_name, f_curve_name_ptr, [len_trim(f_curve_name) + 1]) ! output-only string
  call to_c_str(f_curve_name, f_curve_name_ptr)
end subroutine
subroutine fortran_tao_curve_rms_calc (curve, who, rms, mean) bind(c)

  use tao_struct, only: tao_curve_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: rms  ! 0D_NOT_real
  real(rp) :: f_rms
  real(c_double), pointer :: f_rms_ptr
  type(c_ptr), intent(in), value :: mean  ! 0D_NOT_real
  real(rp) :: f_mean
  real(c_double), pointer :: f_mean_ptr
  ! ** End of parameters **
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  call tao_curve_rms_calc(f_curve, f_who, f_rms, f_mean)

  ! out: f_rms 0D_NOT_real
  call c_f_pointer(rms, f_rms_ptr)
  f_rms_ptr = f_rms
  ! out: f_mean 0D_NOT_real
  call c_f_pointer(mean, f_mean_ptr)
  f_mean_ptr = f_mean
end subroutine
subroutine fortran_tao_d2_d1_name (d1, show_universe, d2_d1_name) bind(c)

  use tao_struct, only: tao_d1_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: d1  ! 0D_NOT_type
  type(tao_d1_data_struct), pointer :: f_d1
  type(c_ptr), intent(in), value :: show_universe  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_show_universe
  logical, target :: f_show_universe_native
  logical, pointer :: f_show_universe_native_ptr
  logical(c_bool), pointer :: f_show_universe_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: d2_d1_name
  character(len=4096), target :: f_d2_d1_name
  character(kind=c_char), pointer :: f_d2_d1_name_ptr(:)
  ! ** End of parameters **
  ! in: f_d1 0D_NOT_type
  if (.not. c_associated(d1)) return
  call c_f_pointer(d1, f_d1)
  ! in: f_show_universe 0D_NOT_logical
  if (c_associated(show_universe)) then
    call c_f_pointer(show_universe, f_show_universe_ptr)
    f_show_universe_native = f_show_universe_ptr
    f_show_universe_native_ptr => f_show_universe_native
  else
    f_show_universe_native_ptr => null()
  endif
  f_d2_d1_name = tao_d2_d1_name(f_d1, f_show_universe_native_ptr)

  ! out: f_d2_d1_name 0D_NOT_character
  call c_f_pointer(d2_d1_name, f_d2_d1_name_ptr, [len_trim(f_d2_d1_name) + 1]) ! output-only string
  call to_c_str(f_d2_d1_name, f_d2_d1_name_ptr)
end subroutine
subroutine fortran_tao_d2_data_stuffit (u, d2_name, n_d1_data) bind(c)

  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), intent(in), value :: d2_name
  character(len=4096), target :: f_d2_name
  character(kind=c_char), pointer :: f_d2_name_ptr(:)
  type(c_ptr), intent(in), value :: n_d1_data  ! 0D_NOT_integer
  integer(c_int) :: f_n_d1_data
  integer(c_int), pointer :: f_n_d1_data_ptr
  ! ** End of parameters **
  ! inout: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! inout: f_d2_name 0D_NOT_character
  if (.not. c_associated(d2_name)) return
  call c_f_pointer(d2_name, f_d2_name_ptr, [huge(0)])
  call to_f_str(f_d2_name_ptr, f_d2_name)
  ! inout: f_n_d1_data 0D_NOT_integer
  if (c_associated(n_d1_data)) then
    call c_f_pointer(n_d1_data, f_n_d1_data_ptr)
  else
    f_n_d1_data_ptr => null()
  endif
  call tao_d2_data_stuffit(f_u, f_d2_name, f_n_d1_data_ptr)

  ! inout: f_d2_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_n_d1_data 0D_NOT_integer
  ! no output conversion for f_n_d1_data
end subroutine
subroutine fortran_tao_data_check (err) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err
  logical, target :: f_err_native
  logical, pointer :: f_err_native_ptr
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_native = f_err_ptr
    f_err_native_ptr => f_err_native
  else
    f_err_native_ptr => null()
  endif
  call tao_data_check(f_err_native_ptr)

  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_ptr = f_err_native
  else
    ! f_err unset
  endif
end subroutine
subroutine fortran_tao_data_coupling_init (branch) bind(c)

  use bmad_struct, only: branch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  ! ** End of parameters **
  ! in: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  call tao_data_coupling_init(f_branch)

end subroutine
subroutine fortran_tao_data_sanity_check (datum, print_err, default_data_type, uni, is_valid) &
    bind(c)

  use tao_struct, only: tao_data_struct, tao_universe_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  logical(c_bool) :: print_err  ! 0D_NOT_logical
  logical :: f_print_err
  type(c_ptr), intent(in), value :: default_data_type
  character(len=4096), target :: f_default_data_type
  character(kind=c_char), pointer :: f_default_data_type_ptr(:)
  type(c_ptr), value :: uni  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_uni
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_valid  ! 0D_NOT_logical
  logical :: f_is_valid
  logical(c_bool), pointer :: f_is_valid_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_print_err 0D_NOT_logical
  f_print_err = print_err
  ! in: f_default_data_type 0D_NOT_character
  if (.not. c_associated(default_data_type)) return
  call c_f_pointer(default_data_type, f_default_data_type_ptr, [huge(0)])
  call to_f_str(f_default_data_type_ptr, f_default_data_type)
  ! in: f_uni 0D_NOT_type
  if (c_associated(uni))   call c_f_pointer(uni, f_uni)
  f_is_valid = tao_data_sanity_check(f_datum, f_print_err, f_default_data_type, f_uni)

  ! out: f_is_valid 0D_NOT_logical
  call c_f_pointer(is_valid, f_is_valid_ptr)
  f_is_valid_ptr = f_is_valid
end subroutine
subroutine fortran_tao_data_type_substitute (template_, str_out, curve, graph) bind(c)

  use tao_struct, only: tao_curve_struct, tao_graph_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: template_
  character(len=4096), target :: f_template
  character(kind=c_char), pointer :: f_template_ptr(:)
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: str_out
  character(len=4096), target :: f_str_out
  character(kind=c_char), pointer :: f_str_out_ptr(:)
  ! ** End of parameters **
  ! in: f_template 0D_NOT_character
  if (.not. c_associated(template_)) return
  call c_f_pointer(template_, f_template_ptr, [huge(0)])
  call to_f_str(f_template_ptr, f_template)
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_data_type_substitute(f_template, f_str_out, f_curve, f_graph)

  ! out: f_str_out 0D_NOT_character
  call c_f_pointer(str_out, f_str_out_ptr, [len_trim(f_str_out) + 1]) ! output-only string
  call to_c_str(f_str_out, f_str_out_ptr)
end subroutine
subroutine fortran_tao_data_useit_plot_calc (curve, graph, data, check_s_position, &
    most_invalid) bind(c)

  use tao_struct, only: tao_curve_struct, tao_data_struct, tao_graph_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  logical(c_bool) :: check_s_position  ! 0D_NOT_logical
  logical :: f_check_s_position
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: data
  type(tao_data_struct_container_alloc), pointer :: f_data
  type(c_ptr), intent(in), value :: most_invalid
  character(len=4096), target :: f_most_invalid
  character(kind=c_char), pointer :: f_most_invalid_ptr(:)
  ! ** End of parameters **
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  !! container type array (1D_ALLOC_type)
  if (c_associated(data))   call c_f_pointer(data, f_data)
  ! in: f_check_s_position 0D_NOT_logical
  f_check_s_position = check_s_position
  call tao_data_useit_plot_calc(f_curve, f_graph, f_data%data, f_check_s_position, &
      f_most_invalid)

  ! out: f_most_invalid 0D_NOT_character
  call c_f_pointer(most_invalid, f_most_invalid_ptr, [len_trim(f_most_invalid) + 1]) ! output-only string
  call to_c_str(f_most_invalid, f_most_invalid_ptr)
end subroutine
subroutine fortran_tao_datum_has_associated_ele (data_type, branch_geometry, &
    has_associated_ele) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: branch_geometry  ! 0D_NOT_integer
  integer(c_int) :: f_branch_geometry
  integer(c_int), pointer :: f_branch_geometry_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: has_associated_ele  ! 0D_NOT_integer
  integer :: f_has_associated_ele
  integer(c_int), pointer :: f_has_associated_ele_ptr
  ! ** End of parameters **
  ! in: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! in: f_branch_geometry 0D_NOT_integer
  if (c_associated(branch_geometry)) then
    call c_f_pointer(branch_geometry, f_branch_geometry_ptr)
  else
    f_branch_geometry_ptr => null()
  endif
  f_has_associated_ele = tao_datum_has_associated_ele(f_data_type, f_branch_geometry_ptr)

  ! out: f_has_associated_ele 0D_NOT_integer
  call c_f_pointer(has_associated_ele, f_has_associated_ele_ptr)
  f_has_associated_ele_ptr = f_has_associated_ele
end subroutine
subroutine fortran_tao_datum_integrate (datum, branch, s_pos, values, valid_value, why_invalid, &
    result) bind(c)

  use tao_struct, only: tao_data_struct
  use bmad_struct, only: branch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  type(c_ptr), intent(in), value :: s_pos
  type(real_container_alloc), pointer :: f_s_pos
  type(c_ptr), intent(in), value :: values
  type(real_container_alloc), pointer :: f_values
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid_value  ! 0D_NOT_logical
  logical :: f_valid_value
  logical(c_bool), pointer :: f_valid_value_ptr
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  type(c_ptr), intent(in), value :: result  ! 0D_NOT_real
  real(rp) :: f_result
  real(c_double), pointer :: f_result_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  !! container general array (1D_ALLOC_real)
  if (c_associated(s_pos))   call c_f_pointer(s_pos, f_s_pos)
  !! container general array (1D_ALLOC_real)
  if (c_associated(values))   call c_f_pointer(values, f_values)
  f_result = tao_datum_integrate(f_datum, f_branch, f_s_pos%data, f_values%data, f_valid_value, &
      f_why_invalid)

  ! out: f_valid_value 0D_NOT_logical
  call c_f_pointer(valid_value, f_valid_value_ptr)
  f_valid_value_ptr = f_valid_value
  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
  ! out: f_result 0D_NOT_real
  call c_f_pointer(result, f_result_ptr)
  f_result_ptr = f_result
end subroutine
subroutine fortran_tao_datum_name (datum, show_universe, datum_name) bind(c)

  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), intent(in), value :: show_universe  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_show_universe
  logical, target :: f_show_universe_native
  logical, pointer :: f_show_universe_native_ptr
  logical(c_bool), pointer :: f_show_universe_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: datum_name
  character(len=4096), target :: f_datum_name
  character(kind=c_char), pointer :: f_datum_name_ptr(:)
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_show_universe 0D_NOT_logical
  if (c_associated(show_universe)) then
    call c_f_pointer(show_universe, f_show_universe_ptr)
    f_show_universe_native = f_show_universe_ptr
    f_show_universe_native_ptr => f_show_universe_native
  else
    f_show_universe_native_ptr => null()
  endif
  f_datum_name = tao_datum_name(f_datum, f_show_universe_native_ptr)

  ! out: f_datum_name 0D_NOT_character
  call c_f_pointer(datum_name, f_datum_name_ptr, [len_trim(f_datum_name) + 1]) ! output-only string
  call to_c_str(f_datum_name, f_datum_name_ptr)
end subroutine
subroutine fortran_tao_datum_s_position (datum, ele, s_pos) bind(c)

  use tao_struct, only: tao_data_struct
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: s_pos  ! 0D_NOT_real
  real(rp) :: f_s_pos
  real(c_double), pointer :: f_s_pos_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  f_s_pos = tao_datum_s_position(f_datum, f_ele)

  ! out: f_s_pos 0D_NOT_real
  call c_f_pointer(s_pos, f_s_pos_ptr)
  f_s_pos_ptr = f_s_pos
end subroutine
subroutine fortran_tao_de_optimizer (abort) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: abort  ! 0D_NOT_logical
  logical :: f_abort
  logical(c_bool), pointer :: f_abort_ptr
  ! ** End of parameters **
  call tao_de_optimizer(f_abort)

  ! out: f_abort 0D_NOT_logical
  call c_f_pointer(abort, f_abort_ptr)
  f_abort_ptr = f_abort
end subroutine
subroutine fortran_tao_deallocate_plot_cache (plot_cache) bind(c)

  use tao_struct, only: tao_plot_cache_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: plot_cache
  type(tao_plot_cache_struct_container_alloc), pointer :: f_plot_cache
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(plot_cache))   call c_f_pointer(plot_cache, f_plot_cache)
  call tao_deallocate_plot_cache(f_plot_cache%data)

end subroutine
subroutine fortran_tao_destroy_plot_window () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_destroy_plot_window()

end subroutine
subroutine fortran_tao_dmerit_calc () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_dmerit_calc()

end subroutine
subroutine fortran_tao_dmodel_dvar_calc (force_calc, err_flag) bind(c)

  implicit none
  ! ** In parameters **
  logical(c_bool) :: force_calc  ! 0D_NOT_logical
  logical :: f_force_calc
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_force_calc 0D_NOT_logical
  f_force_calc = force_calc
  call tao_dmodel_dvar_calc(f_force_calc, f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_do_wire_scan (ele, theta, beam, moment) bind(c)

  use bmad_struct, only: beam_struct, ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  real(c_double) :: theta  ! 0D_NOT_real
  real(rp) :: f_theta
  type(c_ptr), value :: beam  ! 0D_NOT_type
  type(beam_struct), pointer :: f_beam
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: moment  ! 0D_NOT_real
  real(rp) :: f_moment
  real(c_double), pointer :: f_moment_ptr
  ! ** End of parameters **
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_theta 0D_NOT_real
  f_theta = theta
  ! in: f_beam 0D_NOT_type
  if (.not. c_associated(beam)) return
  call c_f_pointer(beam, f_beam)
  f_moment = tao_do_wire_scan(f_ele, f_theta, f_beam)

  ! out: f_moment 0D_NOT_real
  call c_f_pointer(moment, f_moment_ptr)
  f_moment_ptr = f_moment
end subroutine
subroutine fortran_tao_draw_beam_chamber_wall (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_draw_beam_chamber_wall(f_plot, f_graph)

end subroutine
subroutine fortran_tao_draw_curve_data (plot, graph, curve, have_data) bind(c)

  use tao_struct, only: tao_curve_struct, tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: have_data  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_have_data
  logical, target :: f_have_data_native
  logical, pointer :: f_have_data_native_ptr
  logical(c_bool), pointer :: f_have_data_ptr
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! inout: f_have_data 0D_NOT_logical
  if (c_associated(have_data)) then
    call c_f_pointer(have_data, f_have_data_ptr)
    f_have_data_native = f_have_data_ptr
    f_have_data_native_ptr => f_have_data_native
  else
    f_have_data_native_ptr => null()
  endif
  call tao_draw_curve_data(f_plot, f_graph, f_curve, f_have_data_native_ptr)

  ! inout: f_have_data 0D_NOT_logical
  if (c_associated(have_data)) then
    call c_f_pointer(have_data, f_have_data_ptr)
    f_have_data_ptr = f_have_data_native
  else
    ! f_have_data unset
  endif
end subroutine
subroutine fortran_tao_draw_ele_for_floor_plan (plot, graph, tao_lat, ele, ele_shape, &
    label_name, offset1, offset2) bind(c)

  use tao_struct, only: tao_ele_shape_struct, tao_graph_struct, tao_lattice_struct, tao_plot_struct
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: ele_shape  ! 0D_PTR_type
  type(tao_ele_shape_struct), pointer :: f_ele_shape
  type(c_ptr), intent(in), value :: label_name
  character(len=4096), target :: f_label_name
  character(kind=c_char), pointer :: f_label_name_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: offset1  ! 0D_NOT_real
  real(c_double) :: f_offset1
  real(c_double), pointer :: f_offset1_ptr
  type(c_ptr), intent(in), value :: offset2  ! 0D_NOT_real
  real(c_double) :: f_offset2
  real(c_double), pointer :: f_offset2_ptr
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_ele_shape 0D_PTR_type
  if (.not. c_associated(ele_shape)) return
  call c_f_pointer(ele_shape, f_ele_shape)
  ! in: f_label_name 0D_NOT_character
  if (.not. c_associated(label_name)) return
  call c_f_pointer(label_name, f_label_name_ptr, [huge(0)])
  call to_f_str(f_label_name_ptr, f_label_name)
  ! inout: f_offset1 0D_NOT_real
  if (c_associated(offset1)) then
    call c_f_pointer(offset1, f_offset1_ptr)
  else
    f_offset1_ptr => null()
  endif
  ! inout: f_offset2 0D_NOT_real
  if (c_associated(offset2)) then
    call c_f_pointer(offset2, f_offset2_ptr)
  else
    f_offset2_ptr => null()
  endif
  call tao_draw_ele_for_floor_plan(f_plot, f_graph, f_tao_lat, f_ele, f_ele_shape, &
      f_label_name, f_offset1_ptr, f_offset2_ptr)

  ! inout: f_offset1 0D_NOT_real
  ! no output conversion for f_offset1
  ! inout: f_offset2 0D_NOT_real
  ! no output conversion for f_offset2
end subroutine
subroutine fortran_tao_draw_floor_plan (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_draw_floor_plan(f_plot, f_graph)

end subroutine
subroutine fortran_tao_draw_graph_axes (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_draw_graph_axes(f_plot, f_graph)

end subroutine
subroutine fortran_tao_draw_histogram_data (plot, graph, curve, have_data) bind(c)

  use tao_struct, only: tao_curve_struct, tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: have_data  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_have_data
  logical, target :: f_have_data_native
  logical, pointer :: f_have_data_native_ptr
  logical(c_bool), pointer :: f_have_data_ptr
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! inout: f_have_data 0D_NOT_logical
  if (c_associated(have_data)) then
    call c_f_pointer(have_data, f_have_data_ptr)
    f_have_data_native = f_have_data_ptr
    f_have_data_native_ptr => f_have_data_native
  else
    f_have_data_native_ptr => null()
  endif
  call tao_draw_histogram_data(f_plot, f_graph, f_curve, f_have_data_native_ptr)

  ! inout: f_have_data 0D_NOT_logical
  if (c_associated(have_data)) then
    call c_f_pointer(have_data, f_have_data_ptr)
    f_have_data_ptr = f_have_data_native
  else
    ! f_have_data unset
  endif
end subroutine
subroutine fortran_tao_draw_lat_layout (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_draw_lat_layout(f_plot, f_graph)

end subroutine
subroutine fortran_tao_draw_plots (do_clear) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: do_clear  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_do_clear
  logical, target :: f_do_clear_native
  logical, pointer :: f_do_clear_native_ptr
  logical(c_bool), pointer :: f_do_clear_ptr
  ! ** End of parameters **
  ! in: f_do_clear 0D_NOT_logical
  if (c_associated(do_clear)) then
    call c_f_pointer(do_clear, f_do_clear_ptr)
    f_do_clear_native = f_do_clear_ptr
    f_do_clear_native_ptr => f_do_clear_native
  else
    f_do_clear_native_ptr => null()
  endif
  call tao_draw_plots(f_do_clear_native_ptr)

end subroutine
subroutine fortran_tao_ele_geometry_with_misalignments (datum, ele, valid_value, why_invalid, &
    value) bind(c)

  use tao_struct, only: tao_data_struct
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid_value  ! 0D_NOT_logical
  logical :: f_valid_value
  logical(c_bool), pointer :: f_valid_value_ptr
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  f_value = tao_ele_geometry_with_misalignments(f_datum, f_ele, f_valid_value, f_why_invalid)

  ! out: f_valid_value 0D_NOT_logical
  call c_f_pointer(valid_value, f_valid_value_ptr)
  f_valid_value_ptr = f_valid_value
  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_tao_ele_shape_info (ix_uni, ele, ele_shapes, e_shape, label_name, y1, y2, &
    ix_shape_min) bind(c)

  use bmad_struct, only: ele_struct
  use tao_struct, only: tao_ele_shape_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), intent(in), value :: ele_shapes
  type(tao_ele_shape_struct_container_alloc), pointer :: f_ele_shapes
  ! ** Out parameters **
  type(c_ptr), value :: e_shape  ! 0D_PTR_type
  type(tao_ele_shape_struct), pointer :: f_e_shape
  type(c_ptr), intent(in), value :: label_name
  character(len=4096), target :: f_label_name
  character(kind=c_char), pointer :: f_label_name_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: y1  ! 0D_NOT_real
  real(c_double) :: f_y1
  real(c_double), pointer :: f_y1_ptr
  type(c_ptr), intent(in), value :: y2  ! 0D_NOT_real
  real(c_double) :: f_y2
  real(c_double), pointer :: f_y2_ptr
  type(c_ptr), intent(in), value :: ix_shape_min  ! 0D_NOT_integer
  integer(c_int) :: f_ix_shape_min
  integer(c_int), pointer :: f_ix_shape_min_ptr
  ! ** End of parameters **
  ! in: f_ix_uni 0D_NOT_integer
  f_ix_uni = ix_uni
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  !! container type array (1D_ALLOC_type)
  if (c_associated(ele_shapes))   call c_f_pointer(ele_shapes, f_ele_shapes)
  ! out: f_e_shape 0D_PTR_type
  if (.not. c_associated(e_shape)) return
  call c_f_pointer(e_shape, f_e_shape)
  ! inout: f_y1 0D_NOT_real
  if (c_associated(y1)) then
    call c_f_pointer(y1, f_y1_ptr)
  else
    f_y1_ptr => null()
  endif
  ! inout: f_y2 0D_NOT_real
  if (c_associated(y2)) then
    call c_f_pointer(y2, f_y2_ptr)
  else
    f_y2_ptr => null()
  endif
  ! inout: f_ix_shape_min 0D_NOT_integer
  if (c_associated(ix_shape_min)) then
    call c_f_pointer(ix_shape_min, f_ix_shape_min_ptr)
  else
    f_ix_shape_min_ptr => null()
  endif
  call tao_ele_shape_info(f_ix_uni, f_ele, f_ele_shapes%data, f_e_shape, f_label_name, &
      f_y1_ptr, f_y2_ptr, f_ix_shape_min_ptr)

  ! out: f_e_shape 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
  ! out: f_label_name 0D_NOT_character
  call c_f_pointer(label_name, f_label_name_ptr, [len_trim(f_label_name) + 1]) ! output-only string
  call to_c_str(f_label_name, f_label_name_ptr)
  ! inout: f_y1 0D_NOT_real
  ! no output conversion for f_y1
  ! inout: f_y2 0D_NOT_real
  ! no output conversion for f_y2
  ! inout: f_ix_shape_min 0D_NOT_integer
  ! no output conversion for f_ix_shape_min
end subroutine
subroutine fortran_tao_eval_floor_orbit (datum, ele, orbit, bunch_params, valid_value, &
    why_invalid, value) bind(c)

  use tao_struct, only: tao_data_struct
  use bmad_struct, only: bunch_params_struct, coord_struct, ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: orbit  ! 0D_NOT_type
  type(coord_struct), pointer :: f_orbit
  type(c_ptr), value :: bunch_params  ! 0D_NOT_type
  type(bunch_params_struct), pointer :: f_bunch_params
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid_value  ! 0D_NOT_logical
  logical :: f_valid_value
  logical(c_bool), pointer :: f_valid_value_ptr
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_orbit 0D_NOT_type
  if (.not. c_associated(orbit)) return
  call c_f_pointer(orbit, f_orbit)
  ! in: f_bunch_params 0D_NOT_type
  if (.not. c_associated(bunch_params)) return
  call c_f_pointer(bunch_params, f_bunch_params)
  f_value = tao_eval_floor_orbit(f_datum, f_ele, f_orbit, f_bunch_params, f_valid_value, &
      f_why_invalid)

  ! out: f_valid_value 0D_NOT_logical
  call c_f_pointer(valid_value, f_valid_value_ptr)
  f_valid_value_ptr = f_valid_value
  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, &
    why_invalid, called_from_lat_calc, print_err) bind(c)

  use tao_struct, only: tao_data_struct, tao_lattice_struct, tao_universe_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  type(c_ptr), intent(in), value :: called_from_lat_calc  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_called_from_lat_calc
  logical, target :: f_called_from_lat_calc_native
  logical, pointer :: f_called_from_lat_calc_native_ptr
  logical(c_bool), pointer :: f_called_from_lat_calc_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: datum_value  ! 0D_NOT_real
  real(rp) :: f_datum_value
  real(c_double), pointer :: f_datum_value_ptr
  type(c_ptr), intent(in), value :: valid_value  ! 0D_NOT_logical
  logical :: f_valid_value
  logical(c_bool), pointer :: f_valid_value_ptr
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  character(len=4096), pointer :: f_why_invalid_call_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** End of parameters **
  ! inout: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_called_from_lat_calc 0D_NOT_logical
  if (c_associated(called_from_lat_calc)) then
    call c_f_pointer(called_from_lat_calc, f_called_from_lat_calc_ptr)
    f_called_from_lat_calc_native = f_called_from_lat_calc_ptr
    f_called_from_lat_calc_native_ptr => f_called_from_lat_calc_native
  else
    f_called_from_lat_calc_native_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call tao_evaluate_a_datum(f_datum, f_u, f_tao_lat, f_datum_value, f_valid_value, &
      f_why_invalid_call_ptr, f_called_from_lat_calc_native_ptr, f_print_err_native_ptr)

  ! out: f_datum_value 0D_NOT_real
  call c_f_pointer(datum_value, f_datum_value_ptr)
  f_datum_value_ptr = f_datum_value
  ! out: f_valid_value 0D_NOT_logical
  call c_f_pointer(valid_value, f_valid_value_ptr)
  f_valid_value_ptr = f_valid_value
  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
end subroutine
subroutine fortran_tao_evaluate_datum_at_s (datum, tao_lat, ele, ele_ref, valid_value, err_str, &
    bad_datum, value) bind(c)

  use tao_struct, only: tao_data_struct, tao_lattice_struct
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  type(c_ptr), value :: ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: ele_ref  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele_ref
  logical(c_bool) :: valid_value  ! 0D_NOT_logical
  logical :: f_valid_value
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_str
  character(len=4096), target :: f_err_str
  character(kind=c_char), pointer :: f_err_str_ptr(:)
  type(c_ptr), intent(in), value :: bad_datum  ! 0D_NOT_logical
  logical :: f_bad_datum
  logical(c_bool), pointer :: f_bad_datum_ptr
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_ele 0D_PTR_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_ele_ref 0D_PTR_type
  if (.not. c_associated(ele_ref)) return
  call c_f_pointer(ele_ref, f_ele_ref)
  ! in: f_valid_value 0D_NOT_logical
  f_valid_value = valid_value
  f_value = tao_evaluate_datum_at_s(f_datum, f_tao_lat, f_ele, f_ele_ref, f_valid_value, &
      f_err_str, f_bad_datum)

  ! out: f_err_str 0D_NOT_character
  call c_f_pointer(err_str, f_err_str_ptr, [len_trim(f_err_str) + 1]) ! output-only string
  call to_c_str(f_err_str, f_err_str_ptr)
  ! out: f_bad_datum 0D_NOT_logical
  call c_f_pointer(bad_datum, f_bad_datum_ptr)
  f_bad_datum_ptr = f_bad_datum
  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_tao_evaluate_lat_or_beam_data (err, data_name, values, print_err, &
    default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni, &
    dflt_eval_point, dflt_s_offset) bind(c)

  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data_name
  character(len=4096), target :: f_data_name
  character(kind=c_char), pointer :: f_data_name_ptr(:)
  logical(c_bool) :: print_err  ! 0D_NOT_logical
  logical :: f_print_err
  type(c_ptr), value :: dflt_ele_ref  ! 0D_PTR_type
  type(ele_struct), pointer :: f_dflt_ele_ref
  type(c_ptr), value :: dflt_ele_start  ! 0D_PTR_type
  type(ele_struct), pointer :: f_dflt_ele_start
  type(c_ptr), value :: dflt_ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_dflt_ele
  type(c_ptr), intent(in), value :: dflt_component
  character(len=4096), target :: f_dflt_component
  character(kind=c_char), pointer :: f_dflt_component_ptr(:)
  character(len=4096), pointer :: f_dflt_component_call_ptr
  type(c_ptr), intent(in), value :: dflt_uni  ! 0D_NOT_integer
  integer(c_int) :: f_dflt_uni
  integer(c_int), pointer :: f_dflt_uni_ptr
  type(c_ptr), intent(in), value :: dflt_eval_point  ! 0D_NOT_integer
  integer(c_int) :: f_dflt_eval_point
  integer(c_int), pointer :: f_dflt_eval_point_ptr
  type(c_ptr), intent(in), value :: dflt_s_offset  ! 0D_NOT_real
  real(c_double) :: f_dflt_s_offset
  real(c_double), pointer :: f_dflt_s_offset_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  type(c_ptr), intent(in), value :: values
  type(real_container_alloc), pointer :: f_values
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: default_source
  character(len=4096), target :: f_default_source
  character(kind=c_char), pointer :: f_default_source_ptr(:)
  ! ** End of parameters **
  ! in: f_data_name 0D_NOT_character
  if (.not. c_associated(data_name)) return
  call c_f_pointer(data_name, f_data_name_ptr, [huge(0)])
  call to_f_str(f_data_name_ptr, f_data_name)
  !! container general array (1D_ALLOC_real)
  if (c_associated(values))   call c_f_pointer(values, f_values)
  ! in: f_print_err 0D_NOT_logical
  f_print_err = print_err
  ! inout: f_default_source 0D_NOT_character
  if (.not. c_associated(default_source)) return
  call c_f_pointer(default_source, f_default_source_ptr, [huge(0)])
  call to_f_str(f_default_source_ptr, f_default_source)
  ! in: f_dflt_ele_ref 0D_PTR_type
  if (c_associated(dflt_ele_ref))   call c_f_pointer(dflt_ele_ref, f_dflt_ele_ref)
  ! in: f_dflt_ele_start 0D_PTR_type
  if (c_associated(dflt_ele_start))   call c_f_pointer(dflt_ele_start, f_dflt_ele_start)
  ! in: f_dflt_ele 0D_PTR_type
  if (c_associated(dflt_ele))   call c_f_pointer(dflt_ele, f_dflt_ele)
  ! in: f_dflt_component 0D_NOT_character
  if (c_associated(dflt_component)) then
    call c_f_pointer(dflt_component, f_dflt_component_ptr, [huge(0)])
    call to_f_str(f_dflt_component_ptr, f_dflt_component)
    f_dflt_component_call_ptr => f_dflt_component
  else
    f_dflt_component_call_ptr => null()
  endif
  ! in: f_dflt_uni 0D_NOT_integer
  if (c_associated(dflt_uni)) then
    call c_f_pointer(dflt_uni, f_dflt_uni_ptr)
  else
    f_dflt_uni_ptr => null()
  endif
  ! in: f_dflt_eval_point 0D_NOT_integer
  if (c_associated(dflt_eval_point)) then
    call c_f_pointer(dflt_eval_point, f_dflt_eval_point_ptr)
  else
    f_dflt_eval_point_ptr => null()
  endif
  ! in: f_dflt_s_offset 0D_NOT_real
  if (c_associated(dflt_s_offset)) then
    call c_f_pointer(dflt_s_offset, f_dflt_s_offset_ptr)
  else
    f_dflt_s_offset_ptr => null()
  endif
  call tao_evaluate_lat_or_beam_data(f_err, f_data_name, f_values%data, f_print_err, &
      f_default_source, f_dflt_ele_ref, f_dflt_ele_start, f_dflt_ele, &
      f_dflt_component_call_ptr, f_dflt_uni_ptr, f_dflt_eval_point_ptr, f_dflt_s_offset_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! inout: f_default_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_evaluate_tune (q_str, q0, delta_input, q_val) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: q_str
  character(len=4096), target :: f_q_str
  character(kind=c_char), pointer :: f_q_str_ptr(:)
  real(c_double) :: q0  ! 0D_NOT_real
  real(rp) :: f_q0
  logical(c_bool) :: delta_input  ! 0D_NOT_logical
  logical :: f_delta_input
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: q_val  ! 0D_NOT_real
  real(rp) :: f_q_val
  real(c_double), pointer :: f_q_val_ptr
  ! ** End of parameters **
  ! in: f_q_str 0D_NOT_character
  if (.not. c_associated(q_str)) return
  call c_f_pointer(q_str, f_q_str_ptr, [huge(0)])
  call to_f_str(f_q_str_ptr, f_q_str)
  ! in: f_q0 0D_NOT_real
  f_q0 = q0
  ! in: f_delta_input 0D_NOT_logical
  f_delta_input = delta_input
  f_q_val = tao_evaluate_tune(f_q_str, f_q0, f_delta_input)

  ! out: f_q_val 0D_NOT_real
  call c_f_pointer(q_val, f_q_val_ptr)
  f_q_val_ptr = f_q_val
end subroutine
subroutine fortran_tao_expression_hash_substitute (expression_in, expression_out, eval_ele) &
    bind(c)

  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: expression_in
  character(len=4096), target :: f_expression_in
  character(kind=c_char), pointer :: f_expression_in_ptr(:)
  type(c_ptr), value :: eval_ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_eval_ele
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: expression_out
  character(len=4096), target :: f_expression_out
  character(kind=c_char), pointer :: f_expression_out_ptr(:)
  ! ** End of parameters **
  ! in: f_expression_in 0D_NOT_character
  if (.not. c_associated(expression_in)) return
  call c_f_pointer(expression_in, f_expression_in_ptr, [huge(0)])
  call to_f_str(f_expression_in_ptr, f_expression_in)
  ! in: f_eval_ele 0D_PTR_type
  if (c_associated(eval_ele))   call c_f_pointer(eval_ele, f_eval_ele)
  call tao_expression_hash_substitute(f_expression_in, f_expression_out, f_eval_ele)

  ! out: f_expression_out 0D_NOT_character
  call c_f_pointer(expression_out, f_expression_out_ptr, [len_trim(f_expression_out) + 1]) ! output-only string
  call to_c_str(f_expression_out, f_expression_out_ptr)
end subroutine
subroutine fortran_tao_find_plot_region (err, where, region, print_flag) bind(c)

  use tao_struct, only: tao_plot_region_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  type(c_ptr), intent(in), value :: print_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_flag
  logical, target :: f_print_flag_native
  logical, pointer :: f_print_flag_native_ptr
  logical(c_bool), pointer :: f_print_flag_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  type(c_ptr), value :: region  ! 0D_PTR_type
  type(tao_plot_region_struct), pointer :: f_region
  ! ** End of parameters **
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! out: f_region 0D_PTR_type
  if (.not. c_associated(region)) return
  call c_f_pointer(region, f_region)
  ! in: f_print_flag 0D_NOT_logical
  if (c_associated(print_flag)) then
    call c_f_pointer(print_flag, f_print_flag_ptr)
    f_print_flag_native = f_print_flag_ptr
    f_print_flag_native_ptr => f_print_flag_native
  else
    f_print_flag_native_ptr => null()
  endif
  call tao_find_plot_region(f_err, f_where, f_region, f_print_flag_native_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! out: f_region 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_fixer (switch_, word1, word2) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: switch_
  character(len=4096), target :: f_switch
  character(kind=c_char), pointer :: f_switch_ptr(:)
  type(c_ptr), intent(in), value :: word1
  character(len=4096), target :: f_word1
  character(kind=c_char), pointer :: f_word1_ptr(:)
  type(c_ptr), intent(in), value :: word2
  character(len=4096), target :: f_word2
  character(kind=c_char), pointer :: f_word2_ptr(:)
  ! ** End of parameters **
  ! in: f_switch 0D_NOT_character
  if (.not. c_associated(switch_)) return
  call c_f_pointer(switch_, f_switch_ptr, [huge(0)])
  call to_f_str(f_switch_ptr, f_switch)
  ! in: f_word1 0D_NOT_character
  if (.not. c_associated(word1)) return
  call c_f_pointer(word1, f_word1_ptr, [huge(0)])
  call to_f_str(f_word1_ptr, f_word1)
  ! in: f_word2 0D_NOT_character
  if (.not. c_associated(word2)) return
  call c_f_pointer(word2, f_word2_ptr, [huge(0)])
  call to_f_str(f_word2_ptr, f_word2)
  call tao_fixer(f_switch, f_word1, f_word2)

end subroutine
subroutine fortran_tao_floor_to_screen (graph, r_floor, x_screen, y_screen) bind(c)

  use tao_struct, only: tao_graph_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: x_screen  ! 0D_NOT_real
  real(rp) :: f_x_screen
  real(c_double), pointer :: f_x_screen_ptr
  type(c_ptr), intent(in), value :: y_screen  ! 0D_NOT_real
  real(rp) :: f_y_screen
  real(c_double), pointer :: f_y_screen_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: r_floor
  real(rp) :: f_r_floor(3)
  real(c_double), pointer :: f_r_floor_ptr(:)
  ! ** End of parameters **
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  !! general array (1D_NOT_real)
  if (c_associated(r_floor)) then
    call c_f_pointer(r_floor, f_r_floor_ptr, [3])
    f_r_floor = f_r_floor_ptr(:)
  else
    f_r_floor_ptr => null()
  endif
  call tao_floor_to_screen(f_graph, f_r_floor, f_x_screen, f_y_screen)

  ! out: f_x_screen 0D_NOT_real
  call c_f_pointer(x_screen, f_x_screen_ptr)
  f_x_screen_ptr = f_x_screen
  ! out: f_y_screen 0D_NOT_real
  call c_f_pointer(y_screen, f_y_screen_ptr)
  f_y_screen_ptr = f_y_screen
end subroutine
subroutine fortran_tao_floor_to_screen_coords (graph, floor, screen) bind(c)

  use tao_struct, only: tao_graph_struct
  use bmad_struct, only: floor_position_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: floor  ! 0D_NOT_type
  type(floor_position_struct), pointer :: f_floor
  ! ** Out parameters **
  type(c_ptr), value :: screen  ! 0D_NOT_type
  type(floor_position_struct), pointer :: f_screen
  ! ** End of parameters **
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_floor 0D_NOT_type
  if (.not. c_associated(floor)) return
  call c_f_pointer(floor, f_floor)
  ! out: f_screen 0D_NOT_type
  if (.not. c_associated(screen)) return
  call c_f_pointer(screen, f_screen)
  call tao_floor_to_screen_coords(f_graph, f_floor, f_screen)

  ! out: f_screen 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_tao_geodesic_lm_optimizer (abort) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: abort  ! 0D_NOT_logical
  logical :: f_abort
  logical(c_bool), pointer :: f_abort_ptr
  ! ** End of parameters **
  call tao_geodesic_lm_optimizer(f_abort)

  ! out: f_abort 0D_NOT_logical
  call c_f_pointer(abort, f_abort_ptr)
  f_abort_ptr = f_abort
end subroutine
subroutine fortran_tao_get_data (data_value, data_weight, data_meas_value, data_ix_dModel) &
    bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: data_value
  type(real_container_alloc), pointer :: f_data_value
  type(c_ptr), intent(in), value :: data_weight
  type(real_container_alloc), pointer :: f_data_weight
  type(c_ptr), intent(in), value :: data_meas_value
  type(real_container_alloc), pointer :: f_data_meas_value
  type(c_ptr), intent(in), value :: data_ix_dModel
  type(integer_container_alloc), pointer :: f_data_ix_dModel
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (c_associated(data_value))   call c_f_pointer(data_value, f_data_value)
  !! container general array (1D_ALLOC_real)
  if (c_associated(data_weight))   call c_f_pointer(data_weight, f_data_weight)
  !! container general array (1D_ALLOC_real)
  if (c_associated(data_meas_value))   call c_f_pointer(data_meas_value, f_data_meas_value)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(data_ix_dModel))   call c_f_pointer(data_ix_dModel, f_data_ix_dModel)
  call tao_get_data(f_data_value%data, f_data_weight%data, f_data_meas_value%data, &
      f_data_ix_dModel%data)

end subroutine
subroutine fortran_tao_get_opt_vars (var_value, var_step, var_delta, var_weight, var_ix, &
    ignore_if_weight_is_zero, ignore_if_not_limited) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var_value
  type(real_container_alloc), pointer :: f_var_value
  type(c_ptr), intent(in), value :: var_step
  type(real_container_alloc), pointer :: f_var_step
  type(c_ptr), intent(in), value :: var_delta
  type(real_container_alloc), pointer :: f_var_delta
  type(c_ptr), intent(in), value :: var_weight
  type(real_container_alloc), pointer :: f_var_weight
  type(c_ptr), intent(in), value :: var_ix
  type(integer_container_alloc), pointer :: f_var_ix
  type(c_ptr), intent(in), value :: ignore_if_weight_is_zero  ! 0D_NOT_logical
  logical :: f_ignore_if_weight_is_zero
  logical(c_bool), pointer :: f_ignore_if_weight_is_zero_ptr
  type(c_ptr), intent(in), value :: ignore_if_not_limited  ! 0D_NOT_logical
  logical :: f_ignore_if_not_limited
  logical(c_bool), pointer :: f_ignore_if_not_limited_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (c_associated(var_value))   call c_f_pointer(var_value, f_var_value)
  !! container general array (1D_ALLOC_real)
  if (c_associated(var_step))   call c_f_pointer(var_step, f_var_step)
  !! container general array (1D_ALLOC_real)
  if (c_associated(var_delta))   call c_f_pointer(var_delta, f_var_delta)
  !! container general array (1D_ALLOC_real)
  if (c_associated(var_weight))   call c_f_pointer(var_weight, f_var_weight)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(var_ix))   call c_f_pointer(var_ix, f_var_ix)
  call tao_get_opt_vars(f_var_value%data, f_var_step%data, f_var_delta%data, f_var_weight%data, &
      f_var_ix%data, f_ignore_if_weight_is_zero, f_ignore_if_not_limited)

  ! out: f_ignore_if_weight_is_zero 0D_NOT_logical
  ! no output conversion for f_ignore_if_weight_is_zero
  ! out: f_ignore_if_not_limited 0D_NOT_logical
  ! no output conversion for f_ignore_if_not_limited
end subroutine
subroutine fortran_tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: prompt_str
  character(len=4096), target :: f_prompt_str
  character(kind=c_char), pointer :: f_prompt_str_ptr(:)
  character(len=4096), pointer :: f_prompt_str_call_ptr
  type(c_ptr), intent(in), value :: wait_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_wait_flag
  logical, target :: f_wait_flag_native
  logical, pointer :: f_wait_flag_native_ptr
  logical(c_bool), pointer :: f_wait_flag_ptr
  type(c_ptr), intent(in), value :: cmd_in
  character(len=4096), target :: f_cmd_in
  character(kind=c_char), pointer :: f_cmd_in_ptr(:)
  character(len=4096), pointer :: f_cmd_in_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: cmd_out
  character(len=4096), target :: f_cmd_out
  character(kind=c_char), pointer :: f_cmd_out_ptr(:)
  ! ** End of parameters **
  ! in: f_prompt_str 0D_NOT_character
  if (c_associated(prompt_str)) then
    call c_f_pointer(prompt_str, f_prompt_str_ptr, [huge(0)])
    call to_f_str(f_prompt_str_ptr, f_prompt_str)
    f_prompt_str_call_ptr => f_prompt_str
  else
    f_prompt_str_call_ptr => null()
  endif
  ! in: f_wait_flag 0D_NOT_logical
  if (c_associated(wait_flag)) then
    call c_f_pointer(wait_flag, f_wait_flag_ptr)
    f_wait_flag_native = f_wait_flag_ptr
    f_wait_flag_native_ptr => f_wait_flag_native
  else
    f_wait_flag_native_ptr => null()
  endif
  ! in: f_cmd_in 0D_NOT_character
  if (c_associated(cmd_in)) then
    call c_f_pointer(cmd_in, f_cmd_in_ptr, [huge(0)])
    call to_f_str(f_cmd_in_ptr, f_cmd_in)
    f_cmd_in_call_ptr => f_cmd_in
  else
    f_cmd_in_call_ptr => null()
  endif
  call tao_get_user_input(f_cmd_out, f_prompt_str_call_ptr, f_wait_flag_native_ptr, &
      f_cmd_in_call_ptr)

  ! out: f_cmd_out 0D_NOT_character
  call c_f_pointer(cmd_out, f_cmd_out_ptr, [len_trim(f_cmd_out) + 1]) ! output-only string
  call to_c_str(f_cmd_out, f_cmd_out_ptr)
end subroutine
subroutine fortran_tao_graph_controller_setup (graph) bind(c)

  use tao_struct, only: tao_graph_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_controller_setup(f_graph)

end subroutine
subroutine fortran_tao_graph_data_setup (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_data_setup(f_plot, f_graph)

end subroutine
subroutine fortran_tao_graph_data_slice_setup (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_data_slice_setup(f_plot, f_graph)

end subroutine
subroutine fortran_tao_graph_dynamic_aperture_setup (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_dynamic_aperture_setup(f_plot, f_graph)

end subroutine
subroutine fortran_tao_graph_histogram_setup (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_histogram_setup(f_plot, f_graph)

end subroutine
subroutine fortran_tao_graph_name (graph, use_region, graph_name) bind(c)

  use tao_struct, only: tao_graph_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), intent(in), value :: use_region  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_use_region
  logical, target :: f_use_region_native
  logical, pointer :: f_use_region_native_ptr
  logical(c_bool), pointer :: f_use_region_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: graph_name
  character(len=4096), target :: f_graph_name
  character(kind=c_char), pointer :: f_graph_name_ptr(:)
  ! ** End of parameters **
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_use_region 0D_NOT_logical
  if (c_associated(use_region)) then
    call c_f_pointer(use_region, f_use_region_ptr)
    f_use_region_native = f_use_region_ptr
    f_use_region_native_ptr => f_use_region_native
  else
    f_use_region_native_ptr => null()
  endif
  f_graph_name = tao_graph_name(f_graph, f_use_region_native_ptr)

  ! out: f_graph_name 0D_NOT_character
  call c_f_pointer(graph_name, f_graph_name_ptr, [len_trim(f_graph_name) + 1]) ! output-only string
  call to_c_str(f_graph_name, f_graph_name_ptr)
end subroutine
subroutine fortran_tao_graph_phase_space_setup (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_phase_space_setup(f_plot, f_graph)

end subroutine
subroutine fortran_tao_graph_s_min_max_calc (graph, branch, s_min, s_max) bind(c)

  use tao_struct, only: tao_graph_struct
  use bmad_struct, only: branch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: s_min  ! 0D_NOT_real
  real(rp) :: f_s_min
  real(c_double), pointer :: f_s_min_ptr
  type(c_ptr), intent(in), value :: s_max  ! 0D_NOT_real
  real(rp) :: f_s_max
  real(c_double), pointer :: f_s_max_ptr
  ! ** End of parameters **
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  call tao_graph_s_min_max_calc(f_graph, f_branch, f_s_min, f_s_max)

  ! out: f_s_min 0D_NOT_real
  call c_f_pointer(s_min, f_s_min_ptr)
  f_s_min_ptr = f_s_min
  ! out: f_s_max 0D_NOT_real
  call c_f_pointer(s_max, f_s_max_ptr)
  f_s_max_ptr = f_s_max
end subroutine
subroutine fortran_tao_graph_setup (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_graph_setup(f_plot, f_graph)

end subroutine
subroutine fortran_tao_init (err_flag) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  call tao_init(f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_init_beam_in_universe (u, beam_init, track_start, track_end, &
    comb_ds_save) bind(c)

  use tao_struct, only: tao_universe_struct
  use bmad_struct, only: beam_init_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), value :: beam_init  ! 0D_NOT_type
  type(beam_init_struct), pointer :: f_beam_init
  type(c_ptr), intent(in), value :: track_start
  character(len=4096), target :: f_track_start
  character(kind=c_char), pointer :: f_track_start_ptr(:)
  type(c_ptr), intent(in), value :: track_end
  character(len=4096), target :: f_track_end
  character(kind=c_char), pointer :: f_track_end_ptr(:)
  type(c_ptr), intent(in), value :: comb_ds_save  ! 0D_NOT_real
  real(c_double) :: f_comb_ds_save
  real(c_double), pointer :: f_comb_ds_save_ptr
  ! ** End of parameters **
  ! inout: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! inout: f_beam_init 0D_NOT_type
  if (.not. c_associated(beam_init)) return
  call c_f_pointer(beam_init, f_beam_init)
  ! inout: f_track_start 0D_NOT_character
  if (.not. c_associated(track_start)) return
  call c_f_pointer(track_start, f_track_start_ptr, [huge(0)])
  call to_f_str(f_track_start_ptr, f_track_start)
  ! inout: f_track_end 0D_NOT_character
  if (.not. c_associated(track_end)) return
  call c_f_pointer(track_end, f_track_end_ptr, [huge(0)])
  call to_f_str(f_track_end_ptr, f_track_end)
  ! inout: f_comb_ds_save 0D_NOT_real
  if (c_associated(comb_ds_save)) then
    call c_f_pointer(comb_ds_save, f_comb_ds_save_ptr)
  else
    f_comb_ds_save_ptr => null()
  endif
  call tao_init_beam_in_universe(f_u, f_beam_init, f_track_start, f_track_end, &
      f_comb_ds_save_ptr)

  ! inout: f_track_start 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_track_end 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_comb_ds_save 0D_NOT_real
  ! no output conversion for f_comb_ds_save
end subroutine
subroutine fortran_tao_init_beams (init_file) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: init_file
  character(len=4096), target :: f_init_file
  character(kind=c_char), pointer :: f_init_file_ptr(:)
  ! ** End of parameters **
  ! in: f_init_file 0D_NOT_character
  if (.not. c_associated(init_file)) return
  call c_f_pointer(init_file, f_init_file_ptr, [huge(0)])
  call to_f_str(f_init_file_ptr, f_init_file)
  call tao_init_beams(f_init_file)

end subroutine
subroutine fortran_tao_init_data (data_file) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data_file
  character(len=4096), target :: f_data_file
  character(kind=c_char), pointer :: f_data_file_ptr(:)
  ! ** End of parameters **
  ! in: f_data_file 0D_NOT_character
  if (.not. c_associated(data_file)) return
  call c_f_pointer(data_file, f_data_file_ptr, [huge(0)])
  call to_f_str(f_data_file_ptr, f_data_file)
  call tao_init_data(f_data_file)

end subroutine
subroutine fortran_tao_init_data_end_stuff () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_init_data_end_stuff()

end subroutine
subroutine fortran_tao_init_data_in_universe (u, n_d2_add, keep_existing_data) bind(c)

  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), intent(in), value :: n_d2_add  ! 0D_NOT_integer
  integer(c_int) :: f_n_d2_add
  integer(c_int), pointer :: f_n_d2_add_ptr
  type(c_ptr), intent(in), value :: keep_existing_data  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_keep_existing_data
  logical, target :: f_keep_existing_data_native
  logical, pointer :: f_keep_existing_data_native_ptr
  logical(c_bool), pointer :: f_keep_existing_data_ptr
  ! ** End of parameters **
  ! inout: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! inout: f_n_d2_add 0D_NOT_integer
  if (c_associated(n_d2_add)) then
    call c_f_pointer(n_d2_add, f_n_d2_add_ptr)
  else
    f_n_d2_add_ptr => null()
  endif
  ! inout: f_keep_existing_data 0D_NOT_logical
  if (c_associated(keep_existing_data)) then
    call c_f_pointer(keep_existing_data, f_keep_existing_data_ptr)
    f_keep_existing_data_native = f_keep_existing_data_ptr
    f_keep_existing_data_native_ptr => f_keep_existing_data_native
  else
    f_keep_existing_data_native_ptr => null()
  endif
  call tao_init_data_in_universe(f_u, f_n_d2_add_ptr, f_keep_existing_data_native_ptr)

  ! inout: f_n_d2_add 0D_NOT_integer
  ! no output conversion for f_n_d2_add
  ! inout: f_keep_existing_data 0D_NOT_logical
  if (c_associated(keep_existing_data)) then
    call c_f_pointer(keep_existing_data, f_keep_existing_data_ptr)
    f_keep_existing_data_ptr = f_keep_existing_data_native
  else
    ! f_keep_existing_data unset
  endif
end subroutine
subroutine fortran_tao_init_dynamic_aperture (init_file) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: init_file
  character(len=4096), target :: f_init_file
  character(kind=c_char), pointer :: f_init_file_ptr(:)
  ! ** End of parameters **
  ! in: f_init_file 0D_NOT_character
  if (.not. c_associated(init_file)) return
  call c_f_pointer(init_file, f_init_file_ptr, [huge(0)])
  call to_f_str(f_init_file_ptr, f_init_file)
  call tao_init_dynamic_aperture(f_init_file)

end subroutine
subroutine fortran_tao_init_find_elements (u, search_string, eles, attribute, found_one) &
    bind(c)

  use tao_struct, only: tao_universe_struct
  use bmad_struct, only: ele_pointer_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), intent(in), value :: search_string
  character(len=4096), target :: f_search_string
  character(kind=c_char), pointer :: f_search_string_ptr(:)
  type(c_ptr), intent(in), value :: attribute
  character(len=4096), target :: f_attribute
  character(kind=c_char), pointer :: f_attribute_ptr(:)
  character(len=4096), pointer :: f_attribute_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: eles
  type(ele_pointer_struct_container_alloc), pointer :: f_eles
  type(c_ptr), intent(in), value :: found_one  ! 0D_NOT_logical
  logical :: f_found_one
  logical(c_bool), pointer :: f_found_one_ptr
  ! ** End of parameters **
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! in: f_search_string 0D_NOT_character
  if (.not. c_associated(search_string)) return
  call c_f_pointer(search_string, f_search_string_ptr, [huge(0)])
  call to_f_str(f_search_string_ptr, f_search_string)
  !! container type array (1D_ALLOC_type)
  if (c_associated(eles))   call c_f_pointer(eles, f_eles)
  ! in: f_attribute 0D_NOT_character
  if (c_associated(attribute)) then
    call c_f_pointer(attribute, f_attribute_ptr, [huge(0)])
    call to_f_str(f_attribute_ptr, f_attribute)
    f_attribute_call_ptr => f_attribute
  else
    f_attribute_call_ptr => null()
  endif
  call tao_init_find_elements(f_u, f_search_string, f_eles%data, f_attribute_call_ptr, &
      f_found_one)

  ! out: f_found_one 0D_NOT_logical
  ! no output conversion for f_found_one
end subroutine
subroutine fortran_tao_init_global (init_file) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: init_file
  character(len=4096), target :: f_init_file
  character(kind=c_char), pointer :: f_init_file_ptr(:)
  ! ** End of parameters **
  ! in: f_init_file 0D_NOT_character
  if (.not. c_associated(init_file)) return
  call c_f_pointer(init_file, f_init_file_ptr, [huge(0)])
  call to_f_str(f_init_file_ptr, f_init_file)
  call tao_init_global(f_init_file)

end subroutine
subroutine fortran_tao_init_lattice (lat_file, err_flag) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: lat_file
  character(len=4096), target :: f_lat_file
  character(kind=c_char), pointer :: f_lat_file_ptr(:)
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical, target :: f_err_flag_native
  logical, pointer :: f_err_flag_native_ptr
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! inout: f_lat_file 0D_NOT_character
  if (.not. c_associated(lat_file)) return
  call c_f_pointer(lat_file, f_lat_file_ptr, [huge(0)])
  call to_f_str(f_lat_file_ptr, f_lat_file)
  ! inout: f_err_flag 0D_NOT_logical
  if (c_associated(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
    f_err_flag_native_ptr => f_err_flag_native
  else
    f_err_flag_native_ptr => null()
  endif
  call tao_init_lattice(f_lat_file, f_err_flag_native_ptr)

  ! inout: f_lat_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_err_flag 0D_NOT_logical
  if (c_associated(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
end subroutine
subroutine fortran_tao_init_plotting (plot_file) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: plot_file
  character(len=4096), target :: f_plot_file
  character(kind=c_char), pointer :: f_plot_file_ptr(:)
  ! ** End of parameters **
  ! inout: f_plot_file 0D_NOT_character
  if (.not. c_associated(plot_file)) return
  call c_f_pointer(plot_file, f_plot_file_ptr, [huge(0)])
  call to_f_str(f_plot_file_ptr, f_plot_file)
  call tao_init_plotting(f_plot_file)

  ! inout: f_plot_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_init_variables (var_file) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: var_file
  character(len=4096), target :: f_var_file
  character(kind=c_char), pointer :: f_var_file_ptr(:)
  ! ** End of parameters **
  ! in: f_var_file 0D_NOT_character
  if (.not. c_associated(var_file)) return
  call c_f_pointer(var_file, f_var_file_ptr, [huge(0)])
  call to_f_str(f_var_file_ptr, f_var_file)
  call tao_init_variables(f_var_file)

end subroutine
subroutine fortran_tao_inject_beam (u, model, ix_branch, beam, init_ok) bind(c)

  use tao_struct, only: tao_lattice_struct, tao_universe_struct
  use bmad_struct, only: beam_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), value :: model  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_model
  integer(c_int) :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  ! ** Out parameters **
  type(c_ptr), value :: beam  ! 0D_NOT_type
  type(beam_struct), pointer :: f_beam
  type(c_ptr), intent(in), value :: init_ok  ! 0D_NOT_logical
  logical :: f_init_ok
  logical(c_bool), pointer :: f_init_ok_ptr
  ! ** End of parameters **
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! in: f_model 0D_NOT_type
  if (.not. c_associated(model)) return
  call c_f_pointer(model, f_model)
  ! in: f_ix_branch 0D_NOT_integer
  f_ix_branch = ix_branch
  ! out: f_beam 0D_NOT_type
  if (.not. c_associated(beam)) return
  call c_f_pointer(beam, f_beam)
  call tao_inject_beam(f_u, f_model, f_ix_branch, f_beam, f_init_ok)

  ! out: f_beam 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
  ! out: f_init_ok 0D_NOT_logical
  call c_f_pointer(init_ok, f_init_ok_ptr)
  f_init_ok_ptr = f_init_ok
end subroutine
subroutine fortran_tao_inject_particle (u, model, ix_branch) bind(c)

  use tao_struct, only: tao_lattice_struct, tao_universe_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), value :: model  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_model
  type(c_ptr), intent(in), value :: ix_branch  ! 0D_NOT_integer
  integer(c_int) :: f_ix_branch
  integer(c_int), pointer :: f_ix_branch_ptr
  ! ** End of parameters **
  ! inout: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! inout: f_model 0D_NOT_type
  if (.not. c_associated(model)) return
  call c_f_pointer(model, f_model)
  ! inout: f_ix_branch 0D_NOT_integer
  if (c_associated(ix_branch)) then
    call c_f_pointer(ix_branch, f_ix_branch_ptr)
  else
    f_ix_branch_ptr => null()
  endif
  call tao_inject_particle(f_u, f_model, f_ix_branch_ptr)

  ! inout: f_ix_branch 0D_NOT_integer
  ! no output conversion for f_ix_branch
end subroutine
subroutine fortran_tao_is_valid_name (name, why_invalid, is_valid) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  type(c_ptr), intent(in), value :: is_valid  ! 0D_NOT_logical
  logical :: f_is_valid
  logical(c_bool), pointer :: f_is_valid_ptr
  ! ** End of parameters **
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  f_is_valid = tao_is_valid_name(f_name, f_why_invalid)

  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
  ! out: f_is_valid 0D_NOT_logical
  call c_f_pointer(is_valid, f_is_valid_ptr)
  f_is_valid_ptr = f_is_valid
end subroutine
subroutine fortran_tao_json_cmd (input_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: input_str
  character(len=4096), target :: f_input_str
  character(kind=c_char), pointer :: f_input_str_ptr(:)
  ! ** End of parameters **
  ! in: f_input_str 0D_NOT_character
  if (.not. c_associated(input_str)) return
  call c_f_pointer(input_str, f_input_str_ptr, [huge(0)])
  call to_f_str(f_input_str_ptr, f_input_str)
  call tao_json_cmd(f_input_str)

end subroutine
subroutine fortran_tao_key_info_to_str (ix_key, ix_min_key, ix_max_key, key_str, header_str) &
    bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_key  ! 0D_NOT_integer
  integer(c_int) :: f_ix_key
  integer(c_int), pointer :: f_ix_key_ptr
  type(c_ptr), intent(in), value :: ix_min_key  ! 0D_NOT_integer
  integer(c_int) :: f_ix_min_key
  integer(c_int), pointer :: f_ix_min_key_ptr
  type(c_ptr), intent(in), value :: ix_max_key  ! 0D_NOT_integer
  integer(c_int) :: f_ix_max_key
  integer(c_int), pointer :: f_ix_max_key_ptr
  type(c_ptr), intent(in), value :: key_str
  character(len=4096), target :: f_key_str
  character(kind=c_char), pointer :: f_key_str_ptr(:)
  type(c_ptr), intent(in), value :: header_str
  character(len=4096), target :: f_header_str
  character(kind=c_char), pointer :: f_header_str_ptr(:)
  ! ** End of parameters **
  ! inout: f_ix_key 0D_NOT_integer
  if (c_associated(ix_key)) then
    call c_f_pointer(ix_key, f_ix_key_ptr)
  else
    f_ix_key_ptr => null()
  endif
  ! inout: f_ix_min_key 0D_NOT_integer
  if (c_associated(ix_min_key)) then
    call c_f_pointer(ix_min_key, f_ix_min_key_ptr)
  else
    f_ix_min_key_ptr => null()
  endif
  ! inout: f_ix_max_key 0D_NOT_integer
  if (c_associated(ix_max_key)) then
    call c_f_pointer(ix_max_key, f_ix_max_key_ptr)
  else
    f_ix_max_key_ptr => null()
  endif
  ! inout: f_key_str 0D_NOT_character
  if (.not. c_associated(key_str)) return
  call c_f_pointer(key_str, f_key_str_ptr, [huge(0)])
  call to_f_str(f_key_str_ptr, f_key_str)
  ! inout: f_header_str 0D_NOT_character
  if (.not. c_associated(header_str)) return
  call c_f_pointer(header_str, f_header_str_ptr, [huge(0)])
  call to_f_str(f_header_str_ptr, f_header_str)
  call tao_key_info_to_str(f_ix_key_ptr, f_ix_min_key_ptr, f_ix_max_key_ptr, f_key_str, &
      f_header_str)

  ! inout: f_ix_key 0D_NOT_integer
  ! no output conversion for f_ix_key
  ! inout: f_ix_min_key 0D_NOT_integer
  ! no output conversion for f_ix_min_key
  ! inout: f_ix_max_key 0D_NOT_integer
  ! no output conversion for f_ix_max_key
  ! inout: f_key_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_header_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_lat_bookkeeper (u, tao_lat, err_flag) bind(c)

  use tao_struct, only: tao_lattice_struct, tao_universe_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(u, f_u)
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(tao_lat, f_tao_lat)
  call tao_lat_bookkeeper(f_u, f_tao_lat, f_err_flag)

  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_lat_emit_calc (plane, emit_type, ele, modes, emit) bind(c)

  use bmad_struct, only: ele_struct, normal_modes_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: plane  ! 0D_NOT_integer
  integer :: f_plane
  integer(c_int) :: emit_type  ! 0D_NOT_integer
  integer :: f_emit_type
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: modes  ! 0D_NOT_type
  type(normal_modes_struct), pointer :: f_modes
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: emit  ! 0D_NOT_real
  real(rp) :: f_emit
  real(c_double), pointer :: f_emit_ptr
  ! ** End of parameters **
  ! in: f_plane 0D_NOT_integer
  f_plane = plane
  ! in: f_emit_type 0D_NOT_integer
  f_emit_type = emit_type
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_modes 0D_NOT_type
  if (.not. c_associated(modes)) return
  call c_f_pointer(modes, f_modes)
  f_emit = tao_lat_emit_calc(f_plane, f_emit_type, f_ele, f_modes)

  ! out: f_emit 0D_NOT_real
  call c_f_pointer(emit, f_emit_ptr)
  f_emit_ptr = f_emit
end subroutine
subroutine fortran_tao_lat_sigma_calc_needed (data_type, data_source, do_lat_sigma) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: do_lat_sigma  ! 0D_NOT_logical
  logical :: f_do_lat_sigma
  logical(c_bool), pointer :: f_do_lat_sigma_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: data_source
  character(len=4096), target :: f_data_source
  character(kind=c_char), pointer :: f_data_source_ptr(:)
  ! ** End of parameters **
  ! inout: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! inout: f_data_source 0D_NOT_character
  if (.not. c_associated(data_source)) return
  call c_f_pointer(data_source, f_data_source_ptr, [huge(0)])
  call to_f_str(f_data_source_ptr, f_data_source)
  f_do_lat_sigma = tao_lat_sigma_calc_needed(f_data_type, f_data_source)

  ! inout: f_data_type 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_data_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_do_lat_sigma 0D_NOT_logical
  call c_f_pointer(do_lat_sigma, f_do_lat_sigma_ptr)
  f_do_lat_sigma_ptr = f_do_lat_sigma
end subroutine
subroutine fortran_tao_lat_sigma_track (tao_lat, calc_ok, ix_branch, print_err, force_calc) &
    bind(c)

  use tao_struct, only: tao_lattice_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  integer(c_int) :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  type(c_ptr), intent(in), value :: force_calc  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_force_calc
  logical, target :: f_force_calc_native
  logical, pointer :: f_force_calc_native_ptr
  logical(c_bool), pointer :: f_force_calc_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  logical(c_bool), pointer :: f_calc_ok_ptr
  ! ** End of parameters **
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_ix_branch 0D_NOT_integer
  f_ix_branch = ix_branch
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  ! in: f_force_calc 0D_NOT_logical
  if (c_associated(force_calc)) then
    call c_f_pointer(force_calc, f_force_calc_ptr)
    f_force_calc_native = f_force_calc_ptr
    f_force_calc_native_ptr => f_force_calc_native
  else
    f_force_calc_native_ptr => null()
  endif
  call tao_lat_sigma_track(f_tao_lat, f_calc_ok, f_ix_branch, f_print_err_native_ptr, &
      f_force_calc_native_ptr)

  ! out: f_calc_ok 0D_NOT_logical
  call c_f_pointer(calc_ok, f_calc_ok_ptr)
  f_calc_ok_ptr = f_calc_ok
end subroutine
subroutine fortran_tao_lattice_branches_equal_tao_lattice_branches (tlb1, tlb2) bind(c)

  use tao_struct, only: tao_lattice_branch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: tlb2
  type(tao_lattice_branch_struct_container_alloc), pointer :: f_tlb2
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: tlb1
  type(tao_lattice_branch_struct_container_alloc), pointer :: f_tlb1
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(tlb1))   call c_f_pointer(tlb1, f_tlb1)
  !! container type array (1D_ALLOC_type)
  if (c_associated(tlb2))   call c_f_pointer(tlb2, f_tlb2)
  call tao_lattice_branches_equal_tao_lattice_branches(f_tlb1%data, f_tlb2%data)

end subroutine
subroutine fortran_tao_lattice_calc (calc_ok, print_err) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  logical(c_bool), pointer :: f_calc_ok_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical :: f_print_err
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** End of parameters **
  call tao_lattice_calc(f_calc_ok, f_print_err)

  ! out: f_calc_ok 0D_NOT_logical
  call c_f_pointer(calc_ok, f_calc_ok_ptr)
  f_calc_ok_ptr = f_calc_ok
  ! out: f_print_err 0D_NOT_logical
  ! no output conversion for f_print_err
end subroutine
subroutine fortran_tao_lattice_equal_tao_lattice (lat1, lat2) bind(c)

  use tao_struct, only: tao_lattice_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: lat2  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_lat2
  ! ** Inout parameters **
  type(c_ptr), value :: lat1  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_lat1
  ! ** End of parameters **
  ! inout: f_lat1 0D_NOT_type
  if (.not. c_associated(lat1)) return
  call c_f_pointer(lat1, f_lat1)
  ! in: f_lat2 0D_NOT_type
  if (.not. c_associated(lat2)) return
  call c_f_pointer(lat2, f_lat2)
  call tao_lattice_equal_tao_lattice(f_lat1, f_lat2)

end subroutine
subroutine fortran_tao_limit_calc (limited) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: limited  ! 0D_NOT_logical
  logical :: f_limited
  logical(c_bool), pointer :: f_limited_ptr
  ! ** End of parameters **
  call tao_limit_calc(f_limited)

  ! out: f_limited 0D_NOT_logical
  call c_f_pointer(limited, f_limited_ptr)
  f_limited_ptr = f_limited
end subroutine
subroutine fortran_tao_lm_optimizer (abort) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: abort  ! 0D_NOT_logical
  logical :: f_abort
  logical(c_bool), pointer :: f_abort_ptr
  ! ** End of parameters **
  call tao_lm_optimizer(f_abort)

  ! out: f_abort 0D_NOT_logical
  call c_f_pointer(abort, f_abort_ptr)
  f_abort_ptr = f_abort
end subroutine
subroutine fortran_tao_lmdif_optimizer (abort) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: abort  ! 0D_NOT_logical
  logical :: f_abort
  logical(c_bool), pointer :: f_abort_ptr
  ! ** End of parameters **
  call tao_lmdif_optimizer(f_abort)

  ! out: f_abort 0D_NOT_logical
  call c_f_pointer(abort, f_abort_ptr)
  f_abort_ptr = f_abort
end subroutine
subroutine fortran_tao_load_this_datum (vec, ele_ref, ele_start, ele, datum_value, valid_value, &
    datum, branch, why_invalid, good) bind(c)

  use bmad_struct, only: branch_struct, ele_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: vec
  type(real_container_alloc), pointer :: f_vec
  type(c_ptr), value :: ele_ref  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele_ref
  type(c_ptr), value :: ele_start  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele_start
  type(c_ptr), value :: ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), intent(in), value :: datum_value  ! 0D_NOT_real
  real(c_double) :: f_datum_value
  real(c_double), pointer :: f_datum_value_ptr
  type(c_ptr), intent(in), value :: valid_value  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_valid_value
  logical, target :: f_valid_value_native
  logical, pointer :: f_valid_value_native_ptr
  logical(c_bool), pointer :: f_valid_value_ptr
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  character(len=4096), pointer :: f_why_invalid_call_ptr
  type(c_ptr), intent(in), value :: good
  type(logical_container_alloc), pointer :: f_good
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (c_associated(vec))   call c_f_pointer(vec, f_vec)
  ! inout: f_ele_ref 0D_PTR_type
  if (.not. c_associated(ele_ref)) return
  call c_f_pointer(ele_ref, f_ele_ref)
  ! inout: f_ele_start 0D_PTR_type
  if (.not. c_associated(ele_start)) return
  call c_f_pointer(ele_start, f_ele_start)
  ! inout: f_ele 0D_PTR_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! inout: f_datum_value 0D_NOT_real
  if (c_associated(datum_value)) then
    call c_f_pointer(datum_value, f_datum_value_ptr)
  else
    f_datum_value_ptr => null()
  endif
  ! inout: f_valid_value 0D_NOT_logical
  if (c_associated(valid_value)) then
    call c_f_pointer(valid_value, f_valid_value_ptr)
    f_valid_value_native = f_valid_value_ptr
    f_valid_value_native_ptr => f_valid_value_native
  else
    f_valid_value_native_ptr => null()
  endif
  ! inout: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! inout: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  ! inout: f_why_invalid 0D_NOT_character
  if (c_associated(why_invalid)) then
    call c_f_pointer(why_invalid, f_why_invalid_ptr, [huge(0)])
    call to_f_str(f_why_invalid_ptr, f_why_invalid)
    f_why_invalid_call_ptr => f_why_invalid
  else
    f_why_invalid_call_ptr => null()
  endif
  !! container general array (1D_ALLOC_logical)
  if (c_associated(good))   call c_f_pointer(good, f_good)
  call tao_load_this_datum(f_vec%data, f_ele_ref, f_ele_start, f_ele, f_datum_value_ptr, &
      f_valid_value_native_ptr, f_datum, f_branch, f_why_invalid_call_ptr, f_good%data)

  ! inout: f_datum_value 0D_NOT_real
  ! no output conversion for f_datum_value
  ! inout: f_valid_value 0D_NOT_logical
  if (c_associated(valid_value)) then
    call c_f_pointer(valid_value, f_valid_value_ptr)
    f_valid_value_ptr = f_valid_value_native
  else
    ! f_valid_value unset
  endif
  ! inout: f_why_invalid 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_locate_all_elements (ele_list, eles, err, ignore_blank) bind(c)

  use bmad_struct, only: ele_pointer_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: ele_list
  character(len=4096), target :: f_ele_list
  character(kind=c_char), pointer :: f_ele_list_ptr(:)
  type(c_ptr), intent(in), value :: ignore_blank  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_blank
  logical, target :: f_ignore_blank_native
  logical, pointer :: f_ignore_blank_native_ptr
  logical(c_bool), pointer :: f_ignore_blank_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: eles
  type(ele_pointer_struct_container_alloc), pointer :: f_eles
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_ele_list 0D_NOT_character
  if (.not. c_associated(ele_list)) return
  call c_f_pointer(ele_list, f_ele_list_ptr, [huge(0)])
  call to_f_str(f_ele_list_ptr, f_ele_list)
  !! container type array (1D_ALLOC_type)
  if (c_associated(eles))   call c_f_pointer(eles, f_eles)
  ! in: f_ignore_blank 0D_NOT_logical
  if (c_associated(ignore_blank)) then
    call c_f_pointer(ignore_blank, f_ignore_blank_ptr)
    f_ignore_blank_native = f_ignore_blank_ptr
    f_ignore_blank_native_ptr => f_ignore_blank_native
  else
    f_ignore_blank_native_ptr => null()
  endif
  call tao_locate_all_elements(f_ele_list, f_eles%data, f_err, f_ignore_blank_native_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, &
    ignore_blank, err_stat_level, above_ubound_is_err, ix_branch, multiple_eles_is_err) bind(c)

  use bmad_struct, only: ele_pointer_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: ele_list
  character(len=4096), target :: f_ele_list
  character(kind=c_char), pointer :: f_ele_list_ptr(:)
  integer(c_int) :: ix_universe  ! 0D_NOT_integer
  integer :: f_ix_universe
  type(c_ptr), intent(in), value :: lat_type  ! 0D_NOT_integer
  integer(c_int) :: f_lat_type
  integer(c_int), pointer :: f_lat_type_ptr
  type(c_ptr), intent(in), value :: ignore_blank  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_ignore_blank
  logical, target :: f_ignore_blank_native
  logical, pointer :: f_ignore_blank_native_ptr
  logical(c_bool), pointer :: f_ignore_blank_ptr
  type(c_ptr), intent(in), value :: err_stat_level  ! 0D_NOT_integer
  integer(c_int) :: f_err_stat_level
  integer(c_int), pointer :: f_err_stat_level_ptr
  type(c_ptr), intent(in), value :: ix_branch  ! 0D_NOT_integer
  integer(c_int) :: f_ix_branch
  integer(c_int), pointer :: f_ix_branch_ptr
  type(c_ptr), intent(in), value :: multiple_eles_is_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_multiple_eles_is_err
  logical, target :: f_multiple_eles_is_err_native
  logical, pointer :: f_multiple_eles_is_err_native_ptr
  logical(c_bool), pointer :: f_multiple_eles_is_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: eles
  type(ele_pointer_struct_container_alloc), pointer :: f_eles
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: above_ubound_is_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_above_ubound_is_err
  logical, target :: f_above_ubound_is_err_native
  logical, pointer :: f_above_ubound_is_err_native_ptr
  logical(c_bool), pointer :: f_above_ubound_is_err_ptr
  ! ** End of parameters **
  ! in: f_ele_list 0D_NOT_character
  if (.not. c_associated(ele_list)) return
  call c_f_pointer(ele_list, f_ele_list_ptr, [huge(0)])
  call to_f_str(f_ele_list_ptr, f_ele_list)
  ! in: f_ix_universe 0D_NOT_integer
  f_ix_universe = ix_universe
  !! container type array (1D_ALLOC_type)
  if (c_associated(eles))   call c_f_pointer(eles, f_eles)
  ! in: f_lat_type 0D_NOT_integer
  if (c_associated(lat_type)) then
    call c_f_pointer(lat_type, f_lat_type_ptr)
  else
    f_lat_type_ptr => null()
  endif
  ! in: f_ignore_blank 0D_NOT_logical
  if (c_associated(ignore_blank)) then
    call c_f_pointer(ignore_blank, f_ignore_blank_ptr)
    f_ignore_blank_native = f_ignore_blank_ptr
    f_ignore_blank_native_ptr => f_ignore_blank_native
  else
    f_ignore_blank_native_ptr => null()
  endif
  ! in: f_err_stat_level 0D_NOT_integer
  if (c_associated(err_stat_level)) then
    call c_f_pointer(err_stat_level, f_err_stat_level_ptr)
  else
    f_err_stat_level_ptr => null()
  endif
  ! inout: f_above_ubound_is_err 0D_NOT_logical
  if (c_associated(above_ubound_is_err)) then
    call c_f_pointer(above_ubound_is_err, f_above_ubound_is_err_ptr)
    f_above_ubound_is_err_native = f_above_ubound_is_err_ptr
    f_above_ubound_is_err_native_ptr => f_above_ubound_is_err_native
  else
    f_above_ubound_is_err_native_ptr => null()
  endif
  ! in: f_ix_branch 0D_NOT_integer
  if (c_associated(ix_branch)) then
    call c_f_pointer(ix_branch, f_ix_branch_ptr)
  else
    f_ix_branch_ptr => null()
  endif
  ! in: f_multiple_eles_is_err 0D_NOT_logical
  if (c_associated(multiple_eles_is_err)) then
    call c_f_pointer(multiple_eles_is_err, f_multiple_eles_is_err_ptr)
    f_multiple_eles_is_err_native = f_multiple_eles_is_err_ptr
    f_multiple_eles_is_err_native_ptr => f_multiple_eles_is_err_native
  else
    f_multiple_eles_is_err_native_ptr => null()
  endif
  call tao_locate_elements(f_ele_list, f_ix_universe, f_eles%data, f_err, f_lat_type_ptr, &
      f_ignore_blank_native_ptr, f_err_stat_level_ptr, f_above_ubound_is_err_native_ptr, &
      f_ix_branch_ptr, f_multiple_eles_is_err_native_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! inout: f_above_ubound_is_err 0D_NOT_logical
  if (c_associated(above_ubound_is_err)) then
    call c_f_pointer(above_ubound_is_err, f_above_ubound_is_err_ptr)
    f_above_ubound_is_err_ptr = f_above_ubound_is_err_native
  else
    ! f_above_ubound_is_err unset
  endif
end subroutine
subroutine fortran_tao_mark_lattice_ele (lat) bind(c)

  use bmad_struct, only: lat_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  call tao_mark_lattice_ele(f_lat)

end subroutine
subroutine fortran_tao_merit (calc_ok, this_merit) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  logical(c_bool), pointer :: f_calc_ok_ptr
  type(c_ptr), intent(in), value :: this_merit  ! 0D_NOT_real
  real(rp) :: f_this_merit
  real(c_double), pointer :: f_this_merit_ptr
  ! ** End of parameters **
  f_this_merit = tao_merit(f_calc_ok)

  ! out: f_calc_ok 0D_NOT_logical
  ! no output conversion for f_calc_ok
  ! out: f_this_merit 0D_NOT_real
  call c_f_pointer(this_merit, f_this_merit_ptr)
  f_this_merit_ptr = f_this_merit
end subroutine
subroutine fortran_tao_next_word (line, word) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: word
  character(len=4096), target :: f_word
  character(kind=c_char), pointer :: f_word_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: line
  character(len=4096), target :: f_line
  character(kind=c_char), pointer :: f_line_ptr(:)
  ! ** End of parameters **
  ! inout: f_line 0D_NOT_character
  if (.not. c_associated(line)) return
  call c_f_pointer(line, f_line_ptr, [huge(0)])
  call to_f_str(f_line_ptr, f_line)
  call tao_next_word(f_line, f_word)

  ! inout: f_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_word 0D_NOT_character
  call c_f_pointer(word, f_word_ptr, [len_trim(f_word) + 1]) ! output-only string
  call to_c_str(f_word, f_word_ptr)
end subroutine
subroutine fortran_tao_one_turn_map_calc_needed (data_type, data_source, do_one_turn_map) &
    bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: do_one_turn_map  ! 0D_NOT_logical
  logical :: f_do_one_turn_map
  logical(c_bool), pointer :: f_do_one_turn_map_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: data_source
  character(len=4096), target :: f_data_source
  character(kind=c_char), pointer :: f_data_source_ptr(:)
  ! ** End of parameters **
  ! inout: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! inout: f_data_source 0D_NOT_character
  if (.not. c_associated(data_source)) return
  call c_f_pointer(data_source, f_data_source_ptr, [huge(0)])
  call to_f_str(f_data_source_ptr, f_data_source)
  f_do_one_turn_map = tao_one_turn_map_calc_needed(f_data_type, f_data_source)

  ! inout: f_data_type 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_data_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_do_one_turn_map 0D_NOT_logical
  call c_f_pointer(do_one_turn_map, f_do_one_turn_map_ptr)
  f_do_one_turn_map_ptr = f_do_one_turn_map
end subroutine
subroutine fortran_tao_open_file (file, iunit, file_name, error_severity, binary) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: file_name
  character(len=4096), target :: f_file_name
  character(kind=c_char), pointer :: f_file_name_ptr(:)
  integer(c_int) :: error_severity  ! 0D_NOT_integer
  integer :: f_error_severity
  type(c_ptr), intent(in), value :: binary  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_binary
  logical, target :: f_binary_native
  logical, pointer :: f_binary_native_ptr
  logical(c_bool), pointer :: f_binary_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: iunit  ! 0D_NOT_integer
  integer :: f_iunit
  integer(c_int), pointer :: f_iunit_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: file
  character(len=4096), target :: f_file
  character(kind=c_char), pointer :: f_file_ptr(:)
  ! ** End of parameters **
  ! inout: f_file 0D_NOT_character
  if (.not. c_associated(file)) return
  call c_f_pointer(file, f_file_ptr, [huge(0)])
  call to_f_str(f_file_ptr, f_file)
  ! in: f_file_name 0D_NOT_character
  if (.not. c_associated(file_name)) return
  call c_f_pointer(file_name, f_file_name_ptr, [huge(0)])
  call to_f_str(f_file_name_ptr, f_file_name)
  ! in: f_error_severity 0D_NOT_integer
  f_error_severity = error_severity
  ! in: f_binary 0D_NOT_logical
  if (c_associated(binary)) then
    call c_f_pointer(binary, f_binary_ptr)
    f_binary_native = f_binary_ptr
    f_binary_native_ptr => f_binary_native
  else
    f_binary_native_ptr => null()
  endif
  call tao_open_file(f_file, f_iunit, f_file_name, f_error_severity, f_binary_native_ptr)

  ! inout: f_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_iunit 0D_NOT_integer
  call c_f_pointer(iunit, f_iunit_ptr)
  f_iunit_ptr = f_iunit
end subroutine
subroutine fortran_tao_open_scratch_file (err, iu) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  type(c_ptr), intent(in), value :: iu  ! 0D_NOT_integer
  integer :: f_iu
  integer(c_int), pointer :: f_iu_ptr
  ! ** End of parameters **
  f_iu = tao_open_scratch_file(f_err)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! out: f_iu 0D_NOT_integer
  call c_f_pointer(iu, f_iu_ptr)
  f_iu_ptr = f_iu
end subroutine
subroutine fortran_tao_optimization_status (datum, why_str) bind(c)

  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: why_str
  character(len=4096), target :: f_why_str
  character(kind=c_char), pointer :: f_why_str_ptr(:)
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  f_why_str = tao_optimization_status(f_datum)

  ! out: f_why_str 0D_NOT_character
  call c_f_pointer(why_str, f_why_str_ptr, [len_trim(f_why_str) + 1]) ! output-only string
  call to_c_str(f_why_str, f_why_str_ptr)
end subroutine
subroutine fortran_tao_orbit_beta_wave_anal (plot) bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  call tao_orbit_beta_wave_anal(f_plot)

end subroutine
subroutine fortran_tao_oreint_building_wall_pt (pt_in, pt_out) bind(c)

  use tao_struct, only: tao_building_wall_point_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: pt_in  ! 0D_NOT_type
  type(tao_building_wall_point_struct), pointer :: f_pt_in
  ! ** Out parameters **
  type(c_ptr), value :: pt_out  ! 0D_NOT_type
  type(tao_building_wall_point_struct), pointer :: f_pt_out
  ! ** End of parameters **
  ! in: f_pt_in 0D_NOT_type
  if (.not. c_associated(pt_in)) return
  call c_f_pointer(pt_in, f_pt_in)
  f_pt_out = tao_oreint_building_wall_pt(f_pt_in)

  ! out: f_pt_out 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_tao_param_value_at_s (dat_name, ele_to_s, ele_here, orbit, err_flag, &
    why_invalid, print_err, bad_datum, value) bind(c)

  use bmad_struct, only: coord_struct, ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ele_to_s  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele_to_s
  type(c_ptr), value :: ele_here  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele_here
  type(c_ptr), value :: orbit  ! 0D_NOT_type
  type(coord_struct), pointer :: f_orbit
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  character(len=4096), pointer :: f_why_invalid_call_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical :: f_print_err
  logical(c_bool), pointer :: f_print_err_ptr
  type(c_ptr), intent(in), value :: bad_datum  ! 0D_NOT_logical
  logical :: f_bad_datum
  logical(c_bool), pointer :: f_bad_datum_ptr
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: dat_name
  character(len=4096), target :: f_dat_name
  character(kind=c_char), pointer :: f_dat_name_ptr(:)
  ! ** End of parameters **
  ! inout: f_dat_name 0D_NOT_character
  if (.not. c_associated(dat_name)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(dat_name, f_dat_name_ptr, [huge(0)])
  call to_f_str(f_dat_name_ptr, f_dat_name)
  ! in: f_ele_to_s 0D_NOT_type
  if (.not. c_associated(ele_to_s)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(ele_to_s, f_ele_to_s)
  ! in: f_ele_here 0D_NOT_type
  if (.not. c_associated(ele_here)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(ele_here, f_ele_here)
  ! in: f_orbit 0D_NOT_type
  if (.not. c_associated(orbit)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(orbit, f_orbit)
  f_value = tao_param_value_at_s(f_dat_name, f_ele_to_s, f_ele_here, f_orbit, f_err_flag, &
      f_why_invalid_call_ptr, f_print_err, f_bad_datum)

  ! inout: f_dat_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
  ! out: f_print_err 0D_NOT_logical
  ! no output conversion for f_print_err
  ! out: f_bad_datum 0D_NOT_logical
  ! no output conversion for f_bad_datum
  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
end subroutine
subroutine fortran_tao_parse_command_args (error, cmd_line) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: cmd_line
  character(len=4096), target :: f_cmd_line
  character(kind=c_char), pointer :: f_cmd_line_ptr(:)
  character(len=4096), pointer :: f_cmd_line_call_ptr
  ! ** End of parameters **
  ! inout: f_cmd_line 0D_NOT_character
  if (c_associated(cmd_line)) then
    call c_f_pointer(cmd_line, f_cmd_line_ptr, [huge(0)])
    call to_f_str(f_cmd_line_ptr, f_cmd_line)
    f_cmd_line_call_ptr => f_cmd_line
  else
    f_cmd_line_call_ptr => null()
  endif
  call tao_parse_command_args(f_error, f_cmd_line_call_ptr)

  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
  ! inout: f_cmd_line 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_parse_element_param_str (err, in_str, uni, element, parameter, where, &
    component) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: in_str
  character(len=4096), target :: f_in_str
  character(kind=c_char), pointer :: f_in_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  type(c_ptr), intent(in), value :: uni
  character(len=4096), target :: f_uni
  character(kind=c_char), pointer :: f_uni_ptr(:)
  type(c_ptr), intent(in), value :: element
  character(len=4096), target :: f_element
  character(kind=c_char), pointer :: f_element_ptr(:)
  type(c_ptr), intent(in), value :: parameter
  character(len=4096), target :: f_parameter
  character(kind=c_char), pointer :: f_parameter_ptr(:)
  type(c_ptr), intent(in), value :: where  ! 0D_NOT_integer
  integer :: f_where
  integer(c_int), pointer :: f_where_ptr
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  ! ** End of parameters **
  ! in: f_in_str 0D_NOT_character
  if (.not. c_associated(in_str)) return
  call c_f_pointer(in_str, f_in_str_ptr, [huge(0)])
  call to_f_str(f_in_str_ptr, f_in_str)
  call tao_parse_element_param_str(f_err, f_in_str, f_uni, f_element, f_parameter, f_where, &
      f_component)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! out: f_uni 0D_NOT_character
  call c_f_pointer(uni, f_uni_ptr, [len_trim(f_uni) + 1]) ! output-only string
  call to_c_str(f_uni, f_uni_ptr)
  ! out: f_element 0D_NOT_character
  call c_f_pointer(element, f_element_ptr, [len_trim(f_element) + 1]) ! output-only string
  call to_c_str(f_element, f_element_ptr)
  ! out: f_parameter 0D_NOT_character
  call c_f_pointer(parameter, f_parameter_ptr, [len_trim(f_parameter) + 1]) ! output-only string
  call to_c_str(f_parameter, f_parameter_ptr)
  ! out: f_where 0D_NOT_integer
  call c_f_pointer(where, f_where_ptr)
  f_where_ptr = f_where
  ! out: f_component 0D_NOT_character
  call c_f_pointer(component, f_component_ptr, [len_trim(f_component) + 1]) ! output-only string
  call to_c_str(f_component, f_component_ptr)
end subroutine
subroutine fortran_tao_particle_data_value (data_type, p, value, err, ele, ix_bunch) bind(c)

  use bmad_struct, only: coord_struct, ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: p
  type(coord_struct_container_alloc), pointer :: f_p
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  integer(c_int) :: ix_bunch  ! 0D_NOT_integer
  integer :: f_ix_bunch
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value
  type(real_container_alloc), pointer :: f_value
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  !! container type array (1D_ALLOC_type)
  if (c_associated(p))   call c_f_pointer(p, f_p)
  !! container general array (1D_ALLOC_real)
  if (c_associated(value))   call c_f_pointer(value, f_value)
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_ix_bunch 0D_NOT_integer
  f_ix_bunch = ix_bunch
  call tao_particle_data_value(f_data_type, f_p%data, f_value%data, f_err, f_ele, f_ix_bunch)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_pause_cmd (time) bind(c)

  implicit none
  ! ** In parameters **
  real(c_double) :: time  ! 0D_NOT_real
  real(rp) :: f_time
  ! ** End of parameters **
  ! in: f_time 0D_NOT_real
  f_time = time
  call tao_pause_cmd(f_time)

end subroutine
subroutine fortran_tao_phase_space_axis_index (data_type, err, ix_axis) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  logical(c_bool) :: err  ! 0D_NOT_logical
  logical :: f_err
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_axis  ! 0D_NOT_integer
  integer :: f_ix_axis
  integer(c_int), pointer :: f_ix_axis_ptr
  ! ** End of parameters **
  ! in: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! in: f_err 0D_NOT_logical
  f_err = err
  f_ix_axis = tao_phase_space_axis_index(f_data_type, f_err)

  ! out: f_ix_axis 0D_NOT_integer
  call c_f_pointer(ix_axis, f_ix_axis_ptr)
  f_ix_axis_ptr = f_ix_axis
end subroutine
subroutine fortran_tao_phase_wave_anal (plot) bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  call tao_phase_wave_anal(f_plot)

end subroutine
subroutine fortran_tao_pick_universe (name_in, name_out, picked, err, ix_uni, explicit_uni, &
    dflt_uni, pure_uni) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name_in
  character(len=4096), target :: f_name_in
  character(kind=c_char), pointer :: f_name_in_ptr(:)
  type(c_ptr), intent(in), value :: dflt_uni  ! 0D_NOT_integer
  integer(c_int) :: f_dflt_uni
  integer(c_int), pointer :: f_dflt_uni_ptr
  type(c_ptr), intent(in), value :: pure_uni  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_pure_uni
  logical, target :: f_pure_uni_native
  logical, pointer :: f_pure_uni_native_ptr
  logical(c_bool), pointer :: f_pure_uni_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: name_out
  character(len=4096), target :: f_name_out
  character(kind=c_char), pointer :: f_name_out_ptr(:)
  type(c_ptr), intent(in), value :: picked
  type(logical_container_alloc), pointer :: f_picked
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  type(c_ptr), intent(in), value :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  integer(c_int), pointer :: f_ix_uni_ptr
  type(c_ptr), intent(in), value :: explicit_uni  ! 0D_NOT_logical
  logical :: f_explicit_uni
  logical(c_bool), pointer :: f_explicit_uni_ptr
  ! ** End of parameters **
  ! in: f_name_in 0D_NOT_character
  if (.not. c_associated(name_in)) return
  call c_f_pointer(name_in, f_name_in_ptr, [huge(0)])
  call to_f_str(f_name_in_ptr, f_name_in)
  !! container general array (1D_ALLOC_logical)
  if (c_associated(picked))   call c_f_pointer(picked, f_picked)
  ! in: f_dflt_uni 0D_NOT_integer
  if (c_associated(dflt_uni)) then
    call c_f_pointer(dflt_uni, f_dflt_uni_ptr)
  else
    f_dflt_uni_ptr => null()
  endif
  ! in: f_pure_uni 0D_NOT_logical
  if (c_associated(pure_uni)) then
    call c_f_pointer(pure_uni, f_pure_uni_ptr)
    f_pure_uni_native = f_pure_uni_ptr
    f_pure_uni_native_ptr => f_pure_uni_native
  else
    f_pure_uni_native_ptr => null()
  endif
  call tao_pick_universe(f_name_in, f_name_out, f_picked%data, f_err, f_ix_uni, f_explicit_uni, &
      f_dflt_uni_ptr, f_pure_uni_native_ptr)

  ! out: f_name_out 0D_NOT_character
  call c_f_pointer(name_out, f_name_out_ptr, [len_trim(f_name_out) + 1]) ! output-only string
  call to_c_str(f_name_out, f_name_out_ptr)
  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! out: f_ix_uni 0D_NOT_integer
  ! no output conversion for f_ix_uni
  ! out: f_explicit_uni 0D_NOT_logical
  ! no output conversion for f_explicit_uni
end subroutine
subroutine fortran_tao_pipe_cmd (input_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: input_str
  character(len=4096), target :: f_input_str
  character(kind=c_char), pointer :: f_input_str_ptr(:)
  ! ** End of parameters **
  ! in: f_input_str 0D_NOT_character
  if (.not. c_associated(input_str)) return
  call c_f_pointer(input_str, f_input_str_ptr, [huge(0)])
  call to_f_str(f_input_str_ptr, f_input_str)
  call tao_pipe_cmd(f_input_str)

end subroutine
subroutine fortran_tao_place_cmd (where, who, no_buffer) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: no_buffer  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_no_buffer
  logical, target :: f_no_buffer_native
  logical, pointer :: f_no_buffer_native_ptr
  logical(c_bool), pointer :: f_no_buffer_ptr
  ! ** End of parameters **
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_no_buffer 0D_NOT_logical
  if (c_associated(no_buffer)) then
    call c_f_pointer(no_buffer, f_no_buffer_ptr)
    f_no_buffer_native = f_no_buffer_ptr
    f_no_buffer_native_ptr => f_no_buffer_native
  else
    f_no_buffer_native_ptr => null()
  endif
  call tao_place_cmd(f_where, f_who, f_no_buffer_native_ptr)

end subroutine
subroutine fortran_tao_plot_cmd (where, component) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  ! ** End of parameters **
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  call tao_plot_cmd(f_where, f_component)

end subroutine
subroutine fortran_tao_plot_data (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_plot_data(f_plot, f_graph)

end subroutine
subroutine fortran_tao_plot_histogram (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_plot_histogram(f_plot, f_graph)

end subroutine
subroutine fortran_tao_plot_key_table (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_plot_key_table(f_plot, f_graph)

end subroutine
subroutine fortran_tao_plot_setup () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_plot_setup()

end subroutine
subroutine fortran_tao_plot_struct_transfer (plot_in, plot_out) bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot_in  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot_in
  ! ** Out parameters **
  type(c_ptr), value :: plot_out  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot_out
  ! ** End of parameters **
  ! in: f_plot_in 0D_NOT_type
  if (.not. c_associated(plot_in)) return
  call c_f_pointer(plot_in, f_plot_in)
  ! out: f_plot_out 0D_NOT_type
  if (.not. c_associated(plot_out)) return
  call c_f_pointer(plot_out, f_plot_out)
  call tao_plot_struct_transfer(f_plot_in, f_plot_out)

  ! out: f_plot_out 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
end subroutine
subroutine fortran_tao_plot_wave (plot, graph) bind(c)

  use tao_struct, only: tao_graph_struct, tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  call tao_plot_wave(f_plot, f_graph)

end subroutine
subroutine fortran_tao_pointer_to_building_wall_shape (wall_name, e_shape) bind(c)

  use tao_struct, only: tao_ele_shape_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: wall_name
  character(len=4096), target :: f_wall_name
  character(kind=c_char), pointer :: f_wall_name_ptr(:)
  ! ** Out parameters **
  type(c_ptr), value :: e_shape  ! 0D_PTR_type
  type(tao_ele_shape_struct), pointer :: f_e_shape
  ! ** End of parameters **
  ! in: f_wall_name 0D_NOT_character
  if (.not. c_associated(wall_name)) return
  call c_f_pointer(wall_name, f_wall_name_ptr, [huge(0)])
  call to_f_str(f_wall_name_ptr, f_wall_name)
  f_e_shape = tao_pointer_to_building_wall_shape(f_wall_name)

  ! out: f_e_shape 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_datum (d1, ele_name, datum_ptr) bind(c)

  use tao_struct, only: tao_d1_data_struct, tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: d1  ! 0D_NOT_type
  type(tao_d1_data_struct), pointer :: f_d1
  type(c_ptr), intent(in), value :: ele_name
  character(len=4096), target :: f_ele_name
  character(kind=c_char), pointer :: f_ele_name_ptr(:)
  ! ** Out parameters **
  type(c_ptr), value :: datum_ptr  ! 0D_PTR_type
  type(tao_data_struct), pointer :: f_datum_ptr
  ! ** End of parameters **
  ! in: f_d1 0D_NOT_type
  if (.not. c_associated(d1)) return
  call c_f_pointer(d1, f_d1)
  ! in: f_ele_name 0D_NOT_character
  if (.not. c_associated(ele_name)) return
  call c_f_pointer(ele_name, f_ele_name_ptr, [huge(0)])
  call to_f_str(f_ele_name_ptr, f_ele_name)
  f_datum_ptr = tao_pointer_to_datum(f_d1, f_ele_name)

  ! out: f_datum_ptr 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_datum_ele (lat, ele_name, ix_ele, datum, valid, why_invalid, &
    print_err, ele) bind(c)

  use bmad_struct, only: ele_struct, lat_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  integer(c_int) :: ix_ele  ! 0D_NOT_integer
  integer :: f_ix_ele
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: valid  ! 0D_NOT_logical
  logical :: f_valid
  logical(c_bool), pointer :: f_valid_ptr
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  character(len=4096), pointer :: f_why_invalid_call_ptr
  type(c_ptr), value :: ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ele_name
  character(len=4096), target :: f_ele_name
  character(kind=c_char), pointer :: f_ele_name_ptr(:)
  ! ** End of parameters **
  ! in: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  ! inout: f_ele_name 0D_NOT_character
  if (.not. c_associated(ele_name)) return
  call c_f_pointer(ele_name, f_ele_name_ptr, [huge(0)])
  call to_f_str(f_ele_name_ptr, f_ele_name)
  ! in: f_ix_ele 0D_NOT_integer
  f_ix_ele = ix_ele
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  f_ele = tao_pointer_to_datum_ele(f_lat, f_ele_name, f_ix_ele, f_datum, f_valid, &
      f_why_invalid_call_ptr, f_print_err_native_ptr)

  ! inout: f_ele_name 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_valid 0D_NOT_logical
  call c_f_pointer(valid, f_valid_ptr)
  f_valid_ptr = f_valid
  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
  ! out: f_ele 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, &
    dat_var_value, ix_shape_min, e_shape) bind(c)

  use bmad_struct, only: ele_struct
  use tao_struct, only: tao_ele_shape_struct
  implicit none
  ! ** In parameters **
  integer(c_int) :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), intent(in), value :: ele_shape
  type(tao_ele_shape_struct_container_alloc), pointer :: f_ele_shape
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: dat_var_name
  character(len=4096), target :: f_dat_var_name
  character(kind=c_char), pointer :: f_dat_var_name_ptr(:)
  character(len=4096), pointer :: f_dat_var_name_call_ptr
  type(c_ptr), intent(in), value :: dat_var_value  ! 0D_NOT_real
  real(rp) :: f_dat_var_value
  real(c_double), pointer :: f_dat_var_value_ptr
  type(c_ptr), value :: e_shape  ! 0D_PTR_type
  type(tao_ele_shape_struct), pointer :: f_e_shape
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_shape_min  ! 0D_NOT_integer
  integer(c_int) :: f_ix_shape_min
  integer(c_int), pointer :: f_ix_shape_min_ptr
  ! ** End of parameters **
  ! in: f_ix_uni 0D_NOT_integer
  f_ix_uni = ix_uni
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  !! container type array (1D_ALLOC_type)
  if (c_associated(ele_shape))   call c_f_pointer(ele_shape, f_ele_shape)
  ! inout: f_ix_shape_min 0D_NOT_integer
  if (c_associated(ix_shape_min)) then
    call c_f_pointer(ix_shape_min, f_ix_shape_min_ptr)
  else
    f_ix_shape_min_ptr => null()
  endif
  f_e_shape = tao_pointer_to_ele_shape(f_ix_uni, f_ele, f_ele_shape%data, &
      f_dat_var_name_call_ptr, f_dat_var_value, f_ix_shape_min_ptr)

  ! out: f_dat_var_name 0D_NOT_character
  call c_f_pointer(dat_var_name, f_dat_var_name_ptr, [len_trim(f_dat_var_name) + 1]) ! output-only string
  call to_c_str(f_dat_var_name, f_dat_var_name_ptr)
  ! out: f_dat_var_value 0D_NOT_real
  ! no output conversion for f_dat_var_value
  ! inout: f_ix_shape_min 0D_NOT_integer
  ! no output conversion for f_ix_shape_min
  ! out: f_e_shape 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_tao_lat (u, lat_type, tao_lat) bind(c)

  use tao_struct, only: tao_lattice_struct, tao_universe_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  type(c_ptr), intent(in), value :: lat_type  ! 0D_NOT_integer
  integer(c_int) :: f_lat_type
  integer(c_int), pointer :: f_lat_type_ptr
  ! ** Out parameters **
  type(c_ptr), value :: tao_lat  ! 0D_PTR_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  ! ** End of parameters **
  ! in: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  ! in: f_lat_type 0D_NOT_integer
  if (c_associated(lat_type)) then
    call c_f_pointer(lat_type, f_lat_type_ptr)
  else
    f_lat_type_ptr => null()
  endif
  f_tao_lat = tao_pointer_to_tao_lat(f_u, f_lat_type_ptr)

  ! out: f_tao_lat 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_universe_int (ix_uni, neg2_to_default, u) bind(c)

  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), value :: u  ! 0D_PTR_type
  type(tao_universe_struct), pointer :: f_u
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: ix_uni  ! 0D_NOT_integer
  integer(c_int) :: f_ix_uni
  integer(c_int), pointer :: f_ix_uni_ptr
  type(c_ptr), intent(in), value :: neg2_to_default  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_neg2_to_default
  logical, target :: f_neg2_to_default_native
  logical, pointer :: f_neg2_to_default_native_ptr
  logical(c_bool), pointer :: f_neg2_to_default_ptr
  ! ** End of parameters **
  ! inout: f_ix_uni 0D_NOT_integer
  if (c_associated(ix_uni)) then
    call c_f_pointer(ix_uni, f_ix_uni_ptr)
  else
    f_ix_uni_ptr => null()
  endif
  ! inout: f_neg2_to_default 0D_NOT_logical
  if (c_associated(neg2_to_default)) then
    call c_f_pointer(neg2_to_default, f_neg2_to_default_ptr)
    f_neg2_to_default_native = f_neg2_to_default_ptr
    f_neg2_to_default_native_ptr => f_neg2_to_default_native
  else
    f_neg2_to_default_native_ptr => null()
  endif
  f_u = tao_pointer_to_universe(f_ix_uni_ptr, f_neg2_to_default_native_ptr)

  ! inout: f_ix_uni 0D_NOT_integer
  ! no output conversion for f_ix_uni
  ! inout: f_neg2_to_default 0D_NOT_logical
  if (c_associated(neg2_to_default)) then
    call c_f_pointer(neg2_to_default, f_neg2_to_default_ptr)
    f_neg2_to_default_ptr = f_neg2_to_default_native
  else
    ! f_neg2_to_default unset
  endif
  ! out: f_u 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_universe_str (string, neg2_to_default, u) bind(c)

  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), value :: u  ! 0D_PTR_type
  type(tao_universe_struct), pointer :: f_u
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: neg2_to_default  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_neg2_to_default
  logical, target :: f_neg2_to_default_native
  logical, pointer :: f_neg2_to_default_native_ptr
  logical(c_bool), pointer :: f_neg2_to_default_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_neg2_to_default 0D_NOT_logical
  if (c_associated(neg2_to_default)) then
    call c_f_pointer(neg2_to_default, f_neg2_to_default_ptr)
    f_neg2_to_default_native = f_neg2_to_default_ptr
    f_neg2_to_default_native_ptr => f_neg2_to_default_native
  else
    f_neg2_to_default_native_ptr => null()
  endif
  f_u = tao_pointer_to_universe(f_string, f_neg2_to_default_native_ptr)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_neg2_to_default 0D_NOT_logical
  if (c_associated(neg2_to_default)) then
    call c_f_pointer(neg2_to_default, f_neg2_to_default_ptr)
    f_neg2_to_default_ptr = f_neg2_to_default_native
  else
    ! f_neg2_to_default unset
  endif
  ! out: f_u 0D_PTR_type
  ! TODO may require output conversion? 0D_PTR_type
end subroutine
subroutine fortran_tao_pointer_to_universes (name_in, unis, err, name_out, explicit_uni, &
    dflt_uni) bind(c)

  use tao_struct, only: tao_universe_pointer_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name_in
  character(len=4096), target :: f_name_in
  character(kind=c_char), pointer :: f_name_in_ptr(:)
  type(c_ptr), intent(in), value :: dflt_uni  ! 0D_NOT_integer
  integer(c_int) :: f_dflt_uni
  integer(c_int), pointer :: f_dflt_uni_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: unis
  type(tao_universe_pointer_struct_container_alloc), pointer :: f_unis
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  type(c_ptr), intent(in), value :: name_out
  character(len=4096), target :: f_name_out
  character(kind=c_char), pointer :: f_name_out_ptr(:)
  character(len=4096), pointer :: f_name_out_call_ptr
  type(c_ptr), intent(in), value :: explicit_uni  ! 0D_NOT_logical
  logical :: f_explicit_uni
  logical(c_bool), pointer :: f_explicit_uni_ptr
  ! ** End of parameters **
  ! in: f_name_in 0D_NOT_character
  if (.not. c_associated(name_in)) return
  call c_f_pointer(name_in, f_name_in_ptr, [huge(0)])
  call to_f_str(f_name_in_ptr, f_name_in)
  !! container type array (1D_ALLOC_type)
  if (c_associated(unis))   call c_f_pointer(unis, f_unis)
  ! in: f_dflt_uni 0D_NOT_integer
  if (c_associated(dflt_uni)) then
    call c_f_pointer(dflt_uni, f_dflt_uni_ptr)
  else
    f_dflt_uni_ptr => null()
  endif
  call tao_pointer_to_universes(f_name_in, f_unis%data, f_err, f_name_out_call_ptr, &
      f_explicit_uni, f_dflt_uni_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
  ! out: f_name_out 0D_NOT_character
  call c_f_pointer(name_out, f_name_out_ptr, [len_trim(f_name_out) + 1]) ! output-only string
  call to_c_str(f_name_out, f_name_out_ptr)
  ! out: f_explicit_uni 0D_NOT_logical
  ! no output conversion for f_explicit_uni
end subroutine
subroutine fortran_tao_pointer_to_var_in_lattice (var, ix_uni, ele, err) bind(c)

  use tao_struct, only: tao_var_struct
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: var  ! 0D_NOT_type
  type(tao_var_struct), pointer :: f_var
  integer(c_int) :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  ! ** End of parameters **
  ! in: f_var 0D_NOT_type
  if (.not. c_associated(var)) return
  call c_f_pointer(var, f_var)
  ! in: f_ix_uni 0D_NOT_integer
  f_ix_uni = ix_uni
  ! inout: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  call tao_pointer_to_var_in_lattice(f_var, f_ix_uni, f_ele, f_err)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_pointer_to_var_in_lattice2 (var, ix_uni, err) bind(c)

  use tao_struct, only: tao_var_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: var  ! 0D_NOT_type
  type(tao_var_struct), pointer :: f_var
  integer(c_int) :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_var 0D_NOT_type
  if (.not. c_associated(var)) return
  call c_f_pointer(var, f_var)
  ! in: f_ix_uni 0D_NOT_integer
  f_ix_uni = ix_uni
  call tao_pointer_to_var_in_lattice2(f_var, f_ix_uni, f_err)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_print_command_line_info () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_print_command_line_info()

end subroutine
subroutine fortran_tao_ptc_normal_form (do_calc, tao_lat, ix_branch, rf_on) bind(c)

  use tao_struct, only: tao_lattice_struct
  implicit none
  ! ** In parameters **
  logical(c_bool) :: do_calc  ! 0D_NOT_logical
  logical :: f_do_calc
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  integer(c_int) :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  type(c_ptr), intent(in), value :: rf_on  ! 0D_NOT_integer
  integer(c_int) :: f_rf_on
  integer(c_int), pointer :: f_rf_on_ptr
  ! ** End of parameters **
  ! in: f_do_calc 0D_NOT_logical
  f_do_calc = do_calc
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_ix_branch 0D_NOT_integer
  f_ix_branch = ix_branch
  ! in: f_rf_on 0D_NOT_integer
  if (c_associated(rf_on)) then
    call c_f_pointer(rf_on, f_rf_on_ptr)
  else
    f_rf_on_ptr => null()
  endif
  call tao_ptc_normal_form(f_do_calc, f_tao_lat, f_ix_branch, f_rf_on_ptr)

end subroutine
subroutine fortran_tao_python_cmd (input_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: input_str
  character(len=4096), target :: f_input_str
  character(kind=c_char), pointer :: f_input_str_ptr(:)
  ! ** End of parameters **
  ! in: f_input_str 0D_NOT_character
  if (.not. c_associated(input_str)) return
  call c_f_pointer(input_str, f_input_str_ptr, [huge(0)])
  call to_f_str(f_input_str_ptr, f_input_str)
  call tao_python_cmd(f_input_str)

end subroutine
subroutine fortran_tao_quiet_set (set) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: set
  character(len=4096), target :: f_set
  character(kind=c_char), pointer :: f_set_ptr(:)
  ! ** End of parameters **
  ! in: f_set 0D_NOT_character
  if (.not. c_associated(set)) return
  call c_f_pointer(set, f_set_ptr, [huge(0)])
  call to_f_str(f_set_ptr, f_set)
  call tao_quiet_set(f_set)

end subroutine
subroutine fortran_tao_rad_int_calc_needed (data_type, data_source, do_rad_int) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: do_rad_int  ! 0D_NOT_logical
  logical :: f_do_rad_int
  logical(c_bool), pointer :: f_do_rad_int_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: data_source
  character(len=4096), target :: f_data_source
  character(kind=c_char), pointer :: f_data_source_ptr(:)
  ! ** End of parameters **
  ! inout: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! inout: f_data_source 0D_NOT_character
  if (.not. c_associated(data_source)) return
  call c_f_pointer(data_source, f_data_source_ptr, [huge(0)])
  call to_f_str(f_data_source_ptr, f_data_source)
  f_do_rad_int = tao_rad_int_calc_needed(f_data_type, f_data_source)

  ! inout: f_data_type 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_data_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_do_rad_int 0D_NOT_logical
  call c_f_pointer(do_rad_int, f_do_rad_int_ptr)
  f_do_rad_int_ptr = f_do_rad_int
end subroutine
subroutine fortran_tao_re_execute (string, err) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err
  logical, target :: f_err_native
  logical, pointer :: f_err_native_ptr
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! inout: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_native = f_err_ptr
    f_err_native_ptr => f_err_native
  else
    f_err_native_ptr => null()
  endif
  call tao_re_execute(f_string, f_err_native_ptr)

  ! inout: f_string 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_ptr = f_err_native
  else
    ! f_err unset
  endif
end subroutine
subroutine fortran_tao_read_cmd (which, unis, file, silent) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: unis
  character(len=4096), target :: f_unis
  character(kind=c_char), pointer :: f_unis_ptr(:)
  logical(c_bool) :: silent  ! 0D_NOT_logical
  logical :: f_silent
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: which
  character(len=4096), target :: f_which
  character(kind=c_char), pointer :: f_which_ptr(:)
  type(c_ptr), intent(in), value :: file
  character(len=4096), target :: f_file
  character(kind=c_char), pointer :: f_file_ptr(:)
  ! ** End of parameters **
  ! inout: f_which 0D_NOT_character
  if (.not. c_associated(which)) return
  call c_f_pointer(which, f_which_ptr, [huge(0)])
  call to_f_str(f_which_ptr, f_which)
  ! in: f_unis 0D_NOT_character
  if (.not. c_associated(unis)) return
  call c_f_pointer(unis, f_unis_ptr, [huge(0)])
  call to_f_str(f_unis_ptr, f_unis)
  ! inout: f_file 0D_NOT_character
  if (.not. c_associated(file)) return
  call c_f_pointer(file, f_file_ptr, [huge(0)])
  call to_f_str(f_file_ptr, f_file)
  ! in: f_silent 0D_NOT_logical
  f_silent = silent
  call tao_read_cmd(f_which, f_unis, f_file, f_silent)

  ! inout: f_which 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_file 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_read_phase_space_index (name, ixc, print_err, ix_ps) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name
  character(len=4096), target :: f_name
  character(kind=c_char), pointer :: f_name_ptr(:)
  integer(c_int) :: ixc  ! 0D_NOT_integer
  integer :: f_ixc
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_ps  ! 0D_NOT_integer
  integer :: f_ix_ps
  integer(c_int), pointer :: f_ix_ps_ptr
  ! ** End of parameters **
  ! in: f_name 0D_NOT_character
  if (.not. c_associated(name)) return
  call c_f_pointer(name, f_name_ptr, [huge(0)])
  call to_f_str(f_name_ptr, f_name)
  ! in: f_ixc 0D_NOT_integer
  f_ixc = ixc
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  f_ix_ps = tao_read_phase_space_index(f_name, f_ixc, f_print_err_native_ptr)

  ! out: f_ix_ps 0D_NOT_integer
  call c_f_pointer(ix_ps, f_ix_ps_ptr)
  f_ix_ps_ptr = f_ix_ps
end subroutine
subroutine fortran_tao_regression_test () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_regression_test()

end subroutine
subroutine fortran_tao_remove_blank_characters (str) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  call tao_remove_blank_characters(f_str)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_run_cmd (which, abort) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: which
  character(len=4096), target :: f_which
  character(kind=c_char), pointer :: f_which_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: abort  ! 0D_NOT_logical
  logical :: f_abort
  logical(c_bool), pointer :: f_abort_ptr
  ! ** End of parameters **
  ! in: f_which 0D_NOT_character
  if (.not. c_associated(which)) return
  call c_f_pointer(which, f_which_ptr, [huge(0)])
  call to_f_str(f_which_ptr, f_which)
  call tao_run_cmd(f_which, f_abort)

  ! out: f_abort 0D_NOT_logical
  call c_f_pointer(abort, f_abort_ptr)
  f_abort_ptr = f_abort
end subroutine
subroutine fortran_tao_scale_cmd (where, y_min_in, y_max_in, axis, include_wall, gang, exact, &
    turn_autoscale_off) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  real(c_double) :: y_min_in  ! 0D_NOT_real
  real(rp) :: f_y_min_in
  real(c_double) :: y_max_in  ! 0D_NOT_real
  real(rp) :: f_y_max_in
  type(c_ptr), intent(in), value :: axis
  character(len=4096), target :: f_axis
  character(kind=c_char), pointer :: f_axis_ptr(:)
  character(len=4096), pointer :: f_axis_call_ptr
  type(c_ptr), intent(in), value :: include_wall  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_wall
  logical, target :: f_include_wall_native
  logical, pointer :: f_include_wall_native_ptr
  logical(c_bool), pointer :: f_include_wall_ptr
  type(c_ptr), intent(in), value :: gang
  character(len=4096), target :: f_gang
  character(kind=c_char), pointer :: f_gang_ptr(:)
  character(len=4096), pointer :: f_gang_call_ptr
  type(c_ptr), intent(in), value :: exact  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact
  logical, target :: f_exact_native
  logical, pointer :: f_exact_native_ptr
  logical(c_bool), pointer :: f_exact_ptr
  type(c_ptr), intent(in), value :: turn_autoscale_off  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_turn_autoscale_off
  logical, target :: f_turn_autoscale_off_native
  logical, pointer :: f_turn_autoscale_off_native_ptr
  logical(c_bool), pointer :: f_turn_autoscale_off_ptr
  ! ** End of parameters **
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! in: f_y_min_in 0D_NOT_real
  f_y_min_in = y_min_in
  ! in: f_y_max_in 0D_NOT_real
  f_y_max_in = y_max_in
  ! in: f_axis 0D_NOT_character
  if (c_associated(axis)) then
    call c_f_pointer(axis, f_axis_ptr, [huge(0)])
    call to_f_str(f_axis_ptr, f_axis)
    f_axis_call_ptr => f_axis
  else
    f_axis_call_ptr => null()
  endif
  ! in: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_native = f_include_wall_ptr
    f_include_wall_native_ptr => f_include_wall_native
  else
    f_include_wall_native_ptr => null()
  endif
  ! in: f_gang 0D_NOT_character
  if (c_associated(gang)) then
    call c_f_pointer(gang, f_gang_ptr, [huge(0)])
    call to_f_str(f_gang_ptr, f_gang)
    f_gang_call_ptr => f_gang
  else
    f_gang_call_ptr => null()
  endif
  ! in: f_exact 0D_NOT_logical
  if (c_associated(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_native = f_exact_ptr
    f_exact_native_ptr => f_exact_native
  else
    f_exact_native_ptr => null()
  endif
  ! in: f_turn_autoscale_off 0D_NOT_logical
  if (c_associated(turn_autoscale_off)) then
    call c_f_pointer(turn_autoscale_off, f_turn_autoscale_off_ptr)
    f_turn_autoscale_off_native = f_turn_autoscale_off_ptr
    f_turn_autoscale_off_native_ptr => f_turn_autoscale_off_native
  else
    f_turn_autoscale_off_native_ptr => null()
  endif
  call tao_scale_cmd(f_where, f_y_min_in, f_y_max_in, f_axis_call_ptr, &
      f_include_wall_native_ptr, f_gang_call_ptr, f_exact_native_ptr, &
      f_turn_autoscale_off_native_ptr)

end subroutine
subroutine fortran_tao_scale_graph (graph, y_min, y_max, axis, include_wall, y_range, y2_range) &
    bind(c)

  use tao_struct, only: tao_graph_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: y_min  ! 0D_NOT_real
  real(rp) :: f_y_min
  real(c_double) :: y_max  ! 0D_NOT_real
  real(rp) :: f_y_max
  type(c_ptr), intent(in), value :: axis
  character(len=4096), target :: f_axis
  character(kind=c_char), pointer :: f_axis_ptr(:)
  character(len=4096), pointer :: f_axis_call_ptr
  type(c_ptr), intent(in), value :: include_wall  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_wall
  logical, target :: f_include_wall_native
  logical, pointer :: f_include_wall_native_ptr
  logical(c_bool), pointer :: f_include_wall_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: y_range
  real(rp) :: f_y_range(2)
  real(c_double), pointer :: f_y_range_ptr(:)
  type(c_ptr), intent(in), value :: y2_range
  real(rp) :: f_y2_range(2)
  real(c_double), pointer :: f_y2_range_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! in: f_y_min 0D_NOT_real
  f_y_min = y_min
  ! in: f_y_max 0D_NOT_real
  f_y_max = y_max
  ! in: f_axis 0D_NOT_character
  if (c_associated(axis)) then
    call c_f_pointer(axis, f_axis_ptr, [huge(0)])
    call to_f_str(f_axis_ptr, f_axis)
    f_axis_call_ptr => f_axis
  else
    f_axis_call_ptr => null()
  endif
  ! in: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_native = f_include_wall_ptr
    f_include_wall_native_ptr => f_include_wall_native
  else
    f_include_wall_native_ptr => null()
  endif
  call tao_scale_graph(f_graph, f_y_min, f_y_max, f_axis_call_ptr, f_include_wall_native_ptr, &
      f_y_range, f_y2_range)

  ! out: f_y_range 1D_NOT_real
  if (c_associated(y_range)) then
    call c_f_pointer(y_range, f_y_range_ptr, [2])
    f_y_range_ptr = f_y_range(:)
  endif
  ! out: f_y2_range 1D_NOT_real
  if (c_associated(y2_range)) then
    call c_f_pointer(y2_range, f_y2_range_ptr, [2])
    f_y2_range_ptr = f_y2_range(:)
  endif
end subroutine
subroutine fortran_tao_scale_ping_data (u) bind(c)

  use tao_struct, only: tao_universe_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: u  ! 0D_NOT_type
  type(tao_universe_struct), pointer :: f_u
  ! ** End of parameters **
  ! inout: f_u 0D_NOT_type
  if (.not. c_associated(u)) return
  call c_f_pointer(u, f_u)
  call tao_scale_ping_data(f_u)

end subroutine
subroutine fortran_tao_scale_plot (plot, y_min_in, y_max_in, axis, include_wall, gang, &
    skip_lat_layout) bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: y_min_in  ! 0D_NOT_real
  real(rp) :: f_y_min_in
  real(c_double) :: y_max_in  ! 0D_NOT_real
  real(rp) :: f_y_max_in
  type(c_ptr), intent(in), value :: axis
  character(len=4096), target :: f_axis
  character(kind=c_char), pointer :: f_axis_ptr(:)
  character(len=4096), pointer :: f_axis_call_ptr
  type(c_ptr), intent(in), value :: include_wall  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_wall
  logical, target :: f_include_wall_native
  logical, pointer :: f_include_wall_native_ptr
  logical(c_bool), pointer :: f_include_wall_ptr
  type(c_ptr), intent(in), value :: gang
  character(len=4096), target :: f_gang
  character(kind=c_char), pointer :: f_gang_ptr(:)
  character(len=4096), pointer :: f_gang_call_ptr
  type(c_ptr), intent(in), value :: skip_lat_layout  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_skip_lat_layout
  logical, target :: f_skip_lat_layout_native
  logical, pointer :: f_skip_lat_layout_native_ptr
  logical(c_bool), pointer :: f_skip_lat_layout_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_y_min_in 0D_NOT_real
  f_y_min_in = y_min_in
  ! in: f_y_max_in 0D_NOT_real
  f_y_max_in = y_max_in
  ! in: f_axis 0D_NOT_character
  if (c_associated(axis)) then
    call c_f_pointer(axis, f_axis_ptr, [huge(0)])
    call to_f_str(f_axis_ptr, f_axis)
    f_axis_call_ptr => f_axis
  else
    f_axis_call_ptr => null()
  endif
  ! in: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_native = f_include_wall_ptr
    f_include_wall_native_ptr => f_include_wall_native
  else
    f_include_wall_native_ptr => null()
  endif
  ! in: f_gang 0D_NOT_character
  if (c_associated(gang)) then
    call c_f_pointer(gang, f_gang_ptr, [huge(0)])
    call to_f_str(f_gang_ptr, f_gang)
    f_gang_call_ptr => f_gang
  else
    f_gang_call_ptr => null()
  endif
  ! in: f_skip_lat_layout 0D_NOT_logical
  if (c_associated(skip_lat_layout)) then
    call c_f_pointer(skip_lat_layout, f_skip_lat_layout_ptr)
    f_skip_lat_layout_native = f_skip_lat_layout_ptr
    f_skip_lat_layout_native_ptr => f_skip_lat_layout_native
  else
    f_skip_lat_layout_native_ptr => null()
  endif
  call tao_scale_plot(f_plot, f_y_min_in, f_y_max_in, f_axis_call_ptr, &
      f_include_wall_native_ptr, f_gang_call_ptr, f_skip_lat_layout_native_ptr)

end subroutine
subroutine fortran_tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit) &
    bind(c)

  use bmad_struct, only: branch_struct, coord_struct, ele_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: ele_ref  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele_ref
  type(c_ptr), value :: ele_start  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele_start
  type(c_ptr), value :: ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), value :: branch  ! 0D_NOT_type
  type(branch_struct), pointer :: f_branch
  type(c_ptr), intent(in), value :: orbit
  type(coord_struct_container_alloc), pointer :: f_orbit
  ! ** End of parameters **
  ! inout: f_ele_ref 0D_PTR_type
  if (.not. c_associated(ele_ref)) return
  call c_f_pointer(ele_ref, f_ele_ref)
  ! inout: f_ele_start 0D_PTR_type
  if (.not. c_associated(ele_start)) return
  call c_f_pointer(ele_start, f_ele_start)
  ! inout: f_ele 0D_PTR_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! inout: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! inout: f_branch 0D_NOT_type
  if (.not. c_associated(branch)) return
  call c_f_pointer(branch, f_branch)
  !! container type array (1D_ALLOC_type)
  if (c_associated(orbit))   call c_f_pointer(orbit, f_orbit)
  call tao_scratch_values_calc(f_ele_ref, f_ele_start, f_ele, f_datum, f_branch, f_orbit%data)

end subroutine
subroutine fortran_tao_set_beam_cmd (who, value_str, branch_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) return
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  call tao_set_beam_cmd(f_who, f_value_str, f_branch_str)

end subroutine
subroutine fortran_tao_set_beam_init_cmd (who, value_str, branch_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) return
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  call tao_set_beam_init_cmd(f_who, f_value_str, f_branch_str)

end subroutine
subroutine fortran_tao_set_bmad_com_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_bmad_com_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_branch_cmd (branch_str, component_str, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  type(c_ptr), intent(in), value :: component_str
  character(len=4096), target :: f_component_str
  character(kind=c_char), pointer :: f_component_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) return
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  ! in: f_component_str 0D_NOT_character
  if (.not. c_associated(component_str)) return
  call c_f_pointer(component_str, f_component_str_ptr, [huge(0)])
  call to_f_str(f_component_str_ptr, f_component_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_branch_cmd(f_branch_str, f_component_str, f_value_str)

end subroutine
subroutine fortran_tao_set_calculate_cmd (switch_) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: switch_
  character(len=4096), target :: f_switch
  character(kind=c_char), pointer :: f_switch_ptr(:)
  character(len=4096), pointer :: f_switch_call_ptr
  ! ** End of parameters **
  ! inout: f_switch 0D_NOT_character
  if (c_associated(switch_)) then
    call c_f_pointer(switch_, f_switch_ptr, [huge(0)])
    call to_f_str(f_switch_ptr, f_switch)
    f_switch_call_ptr => f_switch
  else
    f_switch_call_ptr => null()
  endif
  call tao_set_calculate_cmd(f_switch_call_ptr)

  ! inout: f_switch 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_set_curve_cmd (curve_name, component, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: curve_name
  character(len=4096), target :: f_curve_name
  character(kind=c_char), pointer :: f_curve_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_curve_name 0D_NOT_character
  if (.not. c_associated(curve_name)) return
  call c_f_pointer(curve_name, f_curve_name_ptr, [huge(0)])
  call to_f_str(f_curve_name_ptr, f_curve_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_curve_cmd(f_curve_name, f_component, f_value_str)

end subroutine
subroutine fortran_tao_set_curve_invalid (curve, why_invalid, print_err) bind(c)

  use tao_struct, only: tao_curve_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  ! ** End of parameters **
  ! inout: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_why_invalid 0D_NOT_character
  if (.not. c_associated(why_invalid)) return
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [huge(0)])
  call to_f_str(f_why_invalid_ptr, f_why_invalid)
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call tao_set_curve_invalid(f_curve, f_why_invalid, f_print_err_native_ptr)

end subroutine
subroutine fortran_tao_set_data_cmd (who_str, value_str, silent) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who_str
  character(len=4096), target :: f_who_str
  character(kind=c_char), pointer :: f_who_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: silent  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_silent
  logical, target :: f_silent_native
  logical, pointer :: f_silent_native_ptr
  logical(c_bool), pointer :: f_silent_ptr
  ! ** End of parameters **
  ! in: f_who_str 0D_NOT_character
  if (.not. c_associated(who_str)) return
  call c_f_pointer(who_str, f_who_str_ptr, [huge(0)])
  call to_f_str(f_who_str_ptr, f_who_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  ! inout: f_silent 0D_NOT_logical
  if (c_associated(silent)) then
    call c_f_pointer(silent, f_silent_ptr)
    f_silent_native = f_silent_ptr
    f_silent_native_ptr => f_silent_native
  else
    f_silent_native_ptr => null()
  endif
  call tao_set_data_cmd(f_who_str, f_value_str, f_silent_native_ptr)

  ! inout: f_silent 0D_NOT_logical
  if (c_associated(silent)) then
    call c_f_pointer(silent, f_silent_ptr)
    f_silent_ptr = f_silent_native
  else
    ! f_silent unset
  endif
end subroutine
subroutine fortran_tao_set_data_useit_opt (data) bind(c)

  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: data
  type(tao_data_struct_container_alloc), pointer :: f_data
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(data))   call c_f_pointer(data, f_data)
  call tao_set_data_useit_opt(f_data%data)

end subroutine
subroutine fortran_tao_set_default_cmd (who_str, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who_str
  character(len=4096), target :: f_who_str
  character(kind=c_char), pointer :: f_who_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who_str 0D_NOT_character
  if (.not. c_associated(who_str)) return
  call c_f_pointer(who_str, f_who_str_ptr, [huge(0)])
  call to_f_str(f_who_str_ptr, f_who_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_default_cmd(f_who_str, f_value_str)

end subroutine
subroutine fortran_tao_set_drawing_cmd (drawing, component, value_str) bind(c)

  use tao_struct, only: tao_drawing_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: drawing  ! 0D_NOT_type
  type(tao_drawing_struct), pointer :: f_drawing
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_drawing 0D_NOT_type
  if (.not. c_associated(drawing)) return
  call c_f_pointer(drawing, f_drawing)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_drawing_cmd(f_drawing, f_component, f_value_str)

end subroutine
subroutine fortran_tao_set_dynamic_aperture_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_dynamic_aperture_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_elements_cmd (ele_list, attribute, value, update) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: ele_list
  character(len=4096), target :: f_ele_list
  character(kind=c_char), pointer :: f_ele_list_ptr(:)
  type(c_ptr), intent(in), value :: attribute
  character(len=4096), target :: f_attribute
  character(kind=c_char), pointer :: f_attribute_ptr(:)
  type(c_ptr), intent(in), value :: value
  character(len=4096), target :: f_value
  character(kind=c_char), pointer :: f_value_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: update  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_update
  logical, target :: f_update_native
  logical, pointer :: f_update_native_ptr
  logical(c_bool), pointer :: f_update_ptr
  ! ** End of parameters **
  ! in: f_ele_list 0D_NOT_character
  if (.not. c_associated(ele_list)) return
  call c_f_pointer(ele_list, f_ele_list_ptr, [huge(0)])
  call to_f_str(f_ele_list_ptr, f_ele_list)
  ! in: f_attribute 0D_NOT_character
  if (.not. c_associated(attribute)) return
  call c_f_pointer(attribute, f_attribute_ptr, [huge(0)])
  call to_f_str(f_attribute_ptr, f_attribute)
  ! in: f_value 0D_NOT_character
  if (.not. c_associated(value)) return
  call c_f_pointer(value, f_value_ptr, [huge(0)])
  call to_f_str(f_value_ptr, f_value)
  ! inout: f_update 0D_NOT_logical
  if (c_associated(update)) then
    call c_f_pointer(update, f_update_ptr)
    f_update_native = f_update_ptr
    f_update_native_ptr => f_update_native
  else
    f_update_native_ptr => null()
  endif
  call tao_set_elements_cmd(f_ele_list, f_attribute, f_value, f_update_native_ptr)

  ! inout: f_update 0D_NOT_logical
  if (c_associated(update)) then
    call c_f_pointer(update, f_update_ptr)
    f_update_ptr = f_update_native
  else
    ! f_update unset
  endif
end subroutine
subroutine fortran_tao_set_floor_plan_axis_label (graph, axis_in, axis_out, which) bind(c)

  use tao_struct, only: tao_graph_struct
  use quick_plot_struct, only: qp_axis_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), value :: axis_in  ! 0D_NOT_type
  type(qp_axis_struct), pointer :: f_axis_in
  type(c_ptr), value :: axis_out  ! 0D_NOT_type
  type(qp_axis_struct), pointer :: f_axis_out
  type(c_ptr), intent(in), value :: which
  character(len=4096), target :: f_which
  character(kind=c_char), pointer :: f_which_ptr(:)
  ! ** End of parameters **
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! inout: f_axis_in 0D_NOT_type
  if (.not. c_associated(axis_in)) return
  call c_f_pointer(axis_in, f_axis_in)
  ! inout: f_axis_out 0D_NOT_type
  if (.not. c_associated(axis_out)) return
  call c_f_pointer(axis_out, f_axis_out)
  ! inout: f_which 0D_NOT_character
  if (.not. c_associated(which)) return
  call c_f_pointer(which, f_which_ptr, [huge(0)])
  call to_f_str(f_which_ptr, f_which)
  call tao_set_floor_plan_axis_label(f_graph, f_axis_in, f_axis_out, f_which)

  ! inout: f_which 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_set_geodesic_lm_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_geodesic_lm_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_global_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_global_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_graph_cmd (graph_name, component, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: graph_name
  character(len=4096), target :: f_graph_name
  character(kind=c_char), pointer :: f_graph_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_graph_name 0D_NOT_character
  if (.not. c_associated(graph_name)) return
  call c_f_pointer(graph_name, f_graph_name_ptr, [huge(0)])
  call to_f_str(f_graph_name_ptr, f_graph_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_graph_cmd(f_graph_name, f_component, f_value_str)

end subroutine
subroutine fortran_tao_set_integer_value (var, var_str, value_str, error, min_val, max_val, &
    print_err) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: var_str
  character(len=4096), target :: f_var_str
  character(kind=c_char), pointer :: f_var_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  type(c_ptr), intent(in), value :: min_val  ! 0D_NOT_integer
  integer(c_int) :: f_min_val
  integer(c_int), pointer :: f_min_val_ptr
  type(c_ptr), intent(in), value :: max_val  ! 0D_NOT_integer
  integer(c_int) :: f_max_val
  integer(c_int), pointer :: f_max_val_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var  ! 0D_NOT_integer
  integer :: f_var
  integer(c_int), pointer :: f_var_ptr
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  ! ** End of parameters **
  ! in: f_var_str 0D_NOT_character
  if (.not. c_associated(var_str)) return
  call c_f_pointer(var_str, f_var_str_ptr, [huge(0)])
  call to_f_str(f_var_str_ptr, f_var_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  ! in: f_min_val 0D_NOT_integer
  if (c_associated(min_val)) then
    call c_f_pointer(min_val, f_min_val_ptr)
  else
    f_min_val_ptr => null()
  endif
  ! in: f_max_val 0D_NOT_integer
  if (c_associated(max_val)) then
    call c_f_pointer(max_val, f_max_val_ptr)
  else
    f_max_val_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call tao_set_integer_value(f_var, f_var_str, f_value_str, f_error, f_min_val_ptr, &
      f_max_val_ptr, f_print_err_native_ptr)

  ! out: f_var 0D_NOT_integer
  call c_f_pointer(var, f_var_ptr)
  f_var_ptr = f_var
  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
end subroutine
subroutine fortran_tao_set_invalid (datum, message, why_invalid, exterminate, err_level, &
    print_err) bind(c)

  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  type(c_ptr), intent(in), value :: message
  character(len=4096), target :: f_message
  character(kind=c_char), pointer :: f_message_ptr(:)
  type(c_ptr), intent(in), value :: exterminate  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exterminate
  logical, target :: f_exterminate_native
  logical, pointer :: f_exterminate_native_ptr
  logical(c_bool), pointer :: f_exterminate_ptr
  type(c_ptr), intent(in), value :: err_level  ! 0D_NOT_integer
  integer(c_int) :: f_err_level
  integer(c_int), pointer :: f_err_level_ptr
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  character(len=4096), pointer :: f_why_invalid_call_ptr
  ! ** End of parameters **
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  ! in: f_message 0D_NOT_character
  if (.not. c_associated(message)) return
  call c_f_pointer(message, f_message_ptr, [huge(0)])
  call to_f_str(f_message_ptr, f_message)
  ! in: f_exterminate 0D_NOT_logical
  if (c_associated(exterminate)) then
    call c_f_pointer(exterminate, f_exterminate_ptr)
    f_exterminate_native = f_exterminate_ptr
    f_exterminate_native_ptr => f_exterminate_native
  else
    f_exterminate_native_ptr => null()
  endif
  ! in: f_err_level 0D_NOT_integer
  if (c_associated(err_level)) then
    call c_f_pointer(err_level, f_err_level_ptr)
  else
    f_err_level_ptr => null()
  endif
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call tao_set_invalid(f_datum, f_message, f_why_invalid_call_ptr, f_exterminate_native_ptr, &
      f_err_level_ptr, f_print_err_native_ptr)

  ! out: f_why_invalid 0D_NOT_character
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [len_trim(f_why_invalid) + 1]) ! output-only string
  call to_c_str(f_why_invalid, f_why_invalid_ptr)
end subroutine
subroutine fortran_tao_set_key_cmd (key_str, cmd_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: key_str
  character(len=4096), target :: f_key_str
  character(kind=c_char), pointer :: f_key_str_ptr(:)
  type(c_ptr), intent(in), value :: cmd_str
  character(len=4096), target :: f_cmd_str
  character(kind=c_char), pointer :: f_cmd_str_ptr(:)
  ! ** End of parameters **
  ! in: f_key_str 0D_NOT_character
  if (.not. c_associated(key_str)) return
  call c_f_pointer(key_str, f_key_str_ptr, [huge(0)])
  call to_f_str(f_key_str_ptr, f_key_str)
  ! in: f_cmd_str 0D_NOT_character
  if (.not. c_associated(cmd_str)) return
  call c_f_pointer(cmd_str, f_cmd_str_ptr, [huge(0)])
  call to_f_str(f_cmd_str_ptr, f_cmd_str)
  call tao_set_key_cmd(f_key_str, f_cmd_str)

end subroutine
subroutine fortran_tao_set_lattice_cmd (dest_lat, source_lat) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: dest_lat
  character(len=4096), target :: f_dest_lat
  character(kind=c_char), pointer :: f_dest_lat_ptr(:)
  type(c_ptr), intent(in), value :: source_lat
  character(len=4096), target :: f_source_lat
  character(kind=c_char), pointer :: f_source_lat_ptr(:)
  ! ** End of parameters **
  ! in: f_dest_lat 0D_NOT_character
  if (.not. c_associated(dest_lat)) return
  call c_f_pointer(dest_lat, f_dest_lat_ptr, [huge(0)])
  call to_f_str(f_dest_lat_ptr, f_dest_lat)
  ! in: f_source_lat 0D_NOT_character
  if (.not. c_associated(source_lat)) return
  call c_f_pointer(source_lat, f_source_lat_ptr, [huge(0)])
  call to_f_str(f_source_lat_ptr, f_source_lat)
  call tao_set_lattice_cmd(f_dest_lat, f_source_lat)

end subroutine
subroutine fortran_tao_set_logical_value (var, var_str, value_str, error) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: var_str
  character(len=4096), target :: f_var_str
  character(kind=c_char), pointer :: f_var_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var  ! 0D_NOT_logical
  logical :: f_var
  logical(c_bool), pointer :: f_var_ptr
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  ! ** End of parameters **
  ! in: f_var_str 0D_NOT_character
  if (.not. c_associated(var_str)) return
  call c_f_pointer(var_str, f_var_str_ptr, [huge(0)])
  call to_f_str(f_var_str_ptr, f_var_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_logical_value(f_var, f_var_str, f_value_str, f_error)

  ! out: f_var 0D_NOT_logical
  call c_f_pointer(var, f_var_ptr)
  f_var_ptr = f_var
  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
end subroutine
subroutine fortran_tao_set_openmp_n_threads (n_threads) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: n_threads  ! 0D_NOT_integer
  integer :: f_n_threads
  ! ** End of parameters **
  ! in: f_n_threads 0D_NOT_integer
  f_n_threads = n_threads
  call tao_set_openmp_n_threads(f_n_threads)

end subroutine
subroutine fortran_tao_set_opt_vars (var_vec, print_limit_warning) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: var_vec
  type(real_container_alloc), pointer :: f_var_vec
  type(c_ptr), intent(in), value :: print_limit_warning  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_limit_warning
  logical, target :: f_print_limit_warning_native
  logical, pointer :: f_print_limit_warning_native_ptr
  logical(c_bool), pointer :: f_print_limit_warning_ptr
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (c_associated(var_vec))   call c_f_pointer(var_vec, f_var_vec)
  ! in: f_print_limit_warning 0D_NOT_logical
  if (c_associated(print_limit_warning)) then
    call c_f_pointer(print_limit_warning, f_print_limit_warning_ptr)
    f_print_limit_warning_native = f_print_limit_warning_ptr
    f_print_limit_warning_native_ptr => f_print_limit_warning_native
  else
    f_print_limit_warning_native_ptr => null()
  endif
  call tao_set_opt_vars(f_var_vec%data, f_print_limit_warning_native_ptr)

end subroutine
subroutine fortran_tao_set_opti_de_param_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_opti_de_param_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_particle_start_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_particle_start_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_plot_cmd (plot_name, component, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: plot_name
  character(len=4096), target :: f_plot_name
  character(kind=c_char), pointer :: f_plot_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_plot_name 0D_NOT_character
  if (.not. c_associated(plot_name)) return
  call c_f_pointer(plot_name, f_plot_name_ptr, [huge(0)])
  call to_f_str(f_plot_name_ptr, f_plot_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_plot_cmd(f_plot_name, f_component, f_value_str)

end subroutine
subroutine fortran_tao_set_plot_page_cmd (component, value_str, value_str2) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str2
  character(len=4096), target :: f_value_str2
  character(kind=c_char), pointer :: f_value_str2_ptr(:)
  character(len=4096), pointer :: f_value_str2_call_ptr
  ! ** End of parameters **
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  ! in: f_value_str2 0D_NOT_character
  if (c_associated(value_str2)) then
    call c_f_pointer(value_str2, f_value_str2_ptr, [huge(0)])
    call to_f_str(f_value_str2_ptr, f_value_str2)
    f_value_str2_call_ptr => f_value_str2
  else
    f_value_str2_call_ptr => null()
  endif
  call tao_set_plot_page_cmd(f_component, f_value_str, f_value_str2_call_ptr)

end subroutine
subroutine fortran_tao_set_ptc_com_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_ptc_com_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_qp_axis_struct (qp_axis_name, component, qp_axis, value, error, &
    ix_uni) bind(c)

  use quick_plot_struct, only: qp_axis_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: qp_axis_name
  character(len=4096), target :: f_qp_axis_name
  character(kind=c_char), pointer :: f_qp_axis_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value
  character(len=4096), target :: f_value
  character(kind=c_char), pointer :: f_value_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  type(c_ptr), intent(in), value :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  integer(c_int), pointer :: f_ix_uni_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: qp_axis  ! 0D_NOT_type
  type(qp_axis_struct), pointer :: f_qp_axis
  ! ** End of parameters **
  ! in: f_qp_axis_name 0D_NOT_character
  if (.not. c_associated(qp_axis_name)) return
  call c_f_pointer(qp_axis_name, f_qp_axis_name_ptr, [huge(0)])
  call to_f_str(f_qp_axis_name_ptr, f_qp_axis_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! inout: f_qp_axis 0D_NOT_type
  if (.not. c_associated(qp_axis)) return
  call c_f_pointer(qp_axis, f_qp_axis)
  ! in: f_value 0D_NOT_character
  if (.not. c_associated(value)) return
  call c_f_pointer(value, f_value_ptr, [huge(0)])
  call to_f_str(f_value_ptr, f_value)
  call tao_set_qp_axis_struct(f_qp_axis_name, f_component, f_qp_axis, f_value, f_error, &
      f_ix_uni)

  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
  ! out: f_ix_uni 0D_NOT_integer
  ! no output conversion for f_ix_uni
end subroutine
subroutine fortran_tao_set_qp_point_struct (qp_point_name, component, qp_point, value, error, &
    ix_uni) bind(c)

  use quick_plot_struct, only: qp_point_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: qp_point_name
  character(len=4096), target :: f_qp_point_name
  character(kind=c_char), pointer :: f_qp_point_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value
  character(len=4096), target :: f_value
  character(kind=c_char), pointer :: f_value_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  type(c_ptr), intent(in), value :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  integer(c_int), pointer :: f_ix_uni_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: qp_point  ! 0D_NOT_type
  type(qp_point_struct), pointer :: f_qp_point
  ! ** End of parameters **
  ! in: f_qp_point_name 0D_NOT_character
  if (.not. c_associated(qp_point_name)) return
  call c_f_pointer(qp_point_name, f_qp_point_name_ptr, [huge(0)])
  call to_f_str(f_qp_point_name_ptr, f_qp_point_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! inout: f_qp_point 0D_NOT_type
  if (.not. c_associated(qp_point)) return
  call c_f_pointer(qp_point, f_qp_point)
  ! in: f_value 0D_NOT_character
  if (.not. c_associated(value)) return
  call c_f_pointer(value, f_value_ptr, [huge(0)])
  call to_f_str(f_value_ptr, f_value)
  call tao_set_qp_point_struct(f_qp_point_name, f_component, f_qp_point, f_value, f_error, &
      f_ix_uni)

  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
  ! out: f_ix_uni 0D_NOT_integer
  ! no output conversion for f_ix_uni
end subroutine
subroutine fortran_tao_set_qp_rect_struct (qp_rect_name, component, qp_rect, value, error, &
    ix_uni) bind(c)

  use quick_plot_struct, only: qp_rect_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: qp_rect_name
  character(len=4096), target :: f_qp_rect_name
  character(kind=c_char), pointer :: f_qp_rect_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value
  character(len=4096), target :: f_value
  character(kind=c_char), pointer :: f_value_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  type(c_ptr), intent(in), value :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  integer(c_int), pointer :: f_ix_uni_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: qp_rect  ! 0D_NOT_type
  type(qp_rect_struct), pointer :: f_qp_rect
  ! ** End of parameters **
  ! in: f_qp_rect_name 0D_NOT_character
  if (.not. c_associated(qp_rect_name)) return
  call c_f_pointer(qp_rect_name, f_qp_rect_name_ptr, [huge(0)])
  call to_f_str(f_qp_rect_name_ptr, f_qp_rect_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! inout: f_qp_rect 0D_NOT_type
  if (.not. c_associated(qp_rect)) return
  call c_f_pointer(qp_rect, f_qp_rect)
  ! in: f_value 0D_NOT_character
  if (.not. c_associated(value)) return
  call c_f_pointer(value, f_value_ptr, [huge(0)])
  call to_f_str(f_value_ptr, f_value)
  call tao_set_qp_rect_struct(f_qp_rect_name, f_component, f_qp_rect, f_value, f_error, &
      f_ix_uni)

  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
  ! out: f_ix_uni 0D_NOT_integer
  ! no output conversion for f_ix_uni
end subroutine
subroutine fortran_tao_set_ran_state_cmd (state_string) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: state_string
  character(len=4096), target :: f_state_string
  character(kind=c_char), pointer :: f_state_string_ptr(:)
  ! ** End of parameters **
  ! in: f_state_string 0D_NOT_character
  if (.not. c_associated(state_string)) return
  call c_f_pointer(state_string, f_state_string_ptr, [huge(0)])
  call to_f_str(f_state_string_ptr, f_state_string)
  call tao_set_ran_state_cmd(f_state_string)

end subroutine
subroutine fortran_tao_set_real_value (var, var_str, value_str, error, min_val, max_val, &
    dflt_uni) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: var_str
  character(len=4096), target :: f_var_str
  character(kind=c_char), pointer :: f_var_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  type(c_ptr), intent(in), value :: min_val  ! 0D_NOT_real
  real(c_double) :: f_min_val
  real(c_double), pointer :: f_min_val_ptr
  type(c_ptr), intent(in), value :: max_val  ! 0D_NOT_real
  real(c_double) :: f_max_val
  real(c_double), pointer :: f_max_val_ptr
  type(c_ptr), intent(in), value :: dflt_uni  ! 0D_NOT_integer
  integer(c_int) :: f_dflt_uni
  integer(c_int), pointer :: f_dflt_uni_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var  ! 0D_NOT_real
  real(rp) :: f_var
  real(c_double), pointer :: f_var_ptr
  type(c_ptr), intent(in), value :: error  ! 0D_NOT_logical
  logical :: f_error
  logical(c_bool), pointer :: f_error_ptr
  ! ** End of parameters **
  ! in: f_var_str 0D_NOT_character
  if (.not. c_associated(var_str)) return
  call c_f_pointer(var_str, f_var_str_ptr, [huge(0)])
  call to_f_str(f_var_str_ptr, f_var_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  ! in: f_min_val 0D_NOT_real
  if (c_associated(min_val)) then
    call c_f_pointer(min_val, f_min_val_ptr)
  else
    f_min_val_ptr => null()
  endif
  ! in: f_max_val 0D_NOT_real
  if (c_associated(max_val)) then
    call c_f_pointer(max_val, f_max_val_ptr)
  else
    f_max_val_ptr => null()
  endif
  ! in: f_dflt_uni 0D_NOT_integer
  if (c_associated(dflt_uni)) then
    call c_f_pointer(dflt_uni, f_dflt_uni_ptr)
  else
    f_dflt_uni_ptr => null()
  endif
  call tao_set_real_value(f_var, f_var_str, f_value_str, f_error, f_min_val_ptr, f_max_val_ptr, &
      f_dflt_uni_ptr)

  ! out: f_var 0D_NOT_real
  call c_f_pointer(var, f_var_ptr)
  f_var_ptr = f_var
  ! out: f_error 0D_NOT_logical
  call c_f_pointer(error, f_error_ptr)
  f_error_ptr = f_error
end subroutine
subroutine fortran_tao_set_region_cmd (region_name, component, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: region_name
  character(len=4096), target :: f_region_name
  character(kind=c_char), pointer :: f_region_name_ptr(:)
  type(c_ptr), intent(in), value :: component
  character(len=4096), target :: f_component
  character(kind=c_char), pointer :: f_component_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_region_name 0D_NOT_character
  if (.not. c_associated(region_name)) return
  call c_f_pointer(region_name, f_region_name_ptr, [huge(0)])
  call to_f_str(f_region_name_ptr, f_region_name)
  ! in: f_component 0D_NOT_character
  if (.not. c_associated(component)) return
  call c_f_pointer(component, f_component_ptr, [huge(0)])
  call to_f_str(f_component_ptr, f_component)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_region_cmd(f_region_name, f_component, f_value_str)

end subroutine
subroutine fortran_tao_set_space_charge_com_cmd (who, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_space_charge_com_cmd(f_who, f_value_str)

end subroutine
subroutine fortran_tao_set_symbolic_number_cmd (sym_str, num_str, val) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: sym_str
  character(len=4096), target :: f_sym_str
  character(kind=c_char), pointer :: f_sym_str_ptr(:)
  type(c_ptr), intent(in), value :: num_str
  character(len=4096), target :: f_num_str
  character(kind=c_char), pointer :: f_num_str_ptr(:)
  character(len=4096), pointer :: f_num_str_call_ptr
  type(c_ptr), intent(in), value :: val  ! 0D_NOT_real
  real(c_double) :: f_val
  real(c_double), pointer :: f_val_ptr
  ! ** End of parameters **
  ! in: f_sym_str 0D_NOT_character
  if (.not. c_associated(sym_str)) return
  call c_f_pointer(sym_str, f_sym_str_ptr, [huge(0)])
  call to_f_str(f_sym_str_ptr, f_sym_str)
  ! in: f_num_str 0D_NOT_character
  if (c_associated(num_str)) then
    call c_f_pointer(num_str, f_num_str_ptr, [huge(0)])
    call to_f_str(f_num_str_ptr, f_num_str)
    f_num_str_call_ptr => f_num_str
  else
    f_num_str_call_ptr => null()
  endif
  ! in: f_val 0D_NOT_real
  if (c_associated(val)) then
    call c_f_pointer(val, f_val_ptr)
  else
    f_val_ptr => null()
  endif
  call tao_set_symbolic_number_cmd(f_sym_str, f_num_str_call_ptr, f_val_ptr)

end subroutine
subroutine fortran_tao_set_tune_cmd (branch_str, mask_str, print_list, qa_str, qb_str, &
    delta_input) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  type(c_ptr), intent(in), value :: mask_str
  character(len=4096), target :: f_mask_str
  character(kind=c_char), pointer :: f_mask_str_ptr(:)
  logical(c_bool) :: print_list  ! 0D_NOT_logical
  logical :: f_print_list
  type(c_ptr), intent(in), value :: qa_str
  character(len=4096), target :: f_qa_str
  character(kind=c_char), pointer :: f_qa_str_ptr(:)
  type(c_ptr), intent(in), value :: qb_str
  character(len=4096), target :: f_qb_str
  character(kind=c_char), pointer :: f_qb_str_ptr(:)
  logical(c_bool) :: delta_input  ! 0D_NOT_logical
  logical :: f_delta_input
  ! ** End of parameters **
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) return
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  ! in: f_mask_str 0D_NOT_character
  if (.not. c_associated(mask_str)) return
  call c_f_pointer(mask_str, f_mask_str_ptr, [huge(0)])
  call to_f_str(f_mask_str_ptr, f_mask_str)
  ! in: f_print_list 0D_NOT_logical
  f_print_list = print_list
  ! in: f_qa_str 0D_NOT_character
  if (.not. c_associated(qa_str)) return
  call c_f_pointer(qa_str, f_qa_str_ptr, [huge(0)])
  call to_f_str(f_qa_str_ptr, f_qa_str)
  ! in: f_qb_str 0D_NOT_character
  if (.not. c_associated(qb_str)) return
  call c_f_pointer(qb_str, f_qb_str_ptr, [huge(0)])
  call to_f_str(f_qb_str_ptr, f_qb_str)
  ! in: f_delta_input 0D_NOT_logical
  f_delta_input = delta_input
  call tao_set_tune_cmd(f_branch_str, f_mask_str, f_print_list, f_qa_str, f_qb_str, &
      f_delta_input)

end subroutine
subroutine fortran_tao_set_universe_cmd (uni, who, what) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: uni
  character(len=4096), target :: f_uni
  character(kind=c_char), pointer :: f_uni_ptr(:)
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: what
  character(len=4096), target :: f_what
  character(kind=c_char), pointer :: f_what_ptr(:)
  ! ** End of parameters **
  ! in: f_uni 0D_NOT_character
  if (.not. c_associated(uni)) return
  call c_f_pointer(uni, f_uni_ptr, [huge(0)])
  call to_f_str(f_uni_ptr, f_uni)
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_what 0D_NOT_character
  if (.not. c_associated(what)) return
  call c_f_pointer(what, f_what_ptr, [huge(0)])
  call to_f_str(f_what_ptr, f_what)
  call tao_set_universe_cmd(f_uni, f_who, f_what)

end subroutine
subroutine fortran_tao_set_var_cmd (var_str, value_str) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: var_str
  character(len=4096), target :: f_var_str
  character(kind=c_char), pointer :: f_var_str_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** End of parameters **
  ! in: f_var_str 0D_NOT_character
  if (.not. c_associated(var_str)) return
  call c_f_pointer(var_str, f_var_str_ptr, [huge(0)])
  call to_f_str(f_var_str_ptr, f_var_str)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_var_cmd(f_var_str, f_value_str)

end subroutine
subroutine fortran_tao_set_var_model_value (var, value, print_limit_warning) bind(c)

  use tao_struct, only: tao_var_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: var  ! 0D_NOT_type
  type(tao_var_struct), pointer :: f_var
  real(c_double) :: value  ! 0D_NOT_real
  real(rp) :: f_value
  type(c_ptr), intent(in), value :: print_limit_warning  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_limit_warning
  logical, target :: f_print_limit_warning_native
  logical, pointer :: f_print_limit_warning_native_ptr
  logical(c_bool), pointer :: f_print_limit_warning_ptr
  ! ** End of parameters **
  ! in: f_var 0D_NOT_type
  if (.not. c_associated(var)) return
  call c_f_pointer(var, f_var)
  ! in: f_value 0D_NOT_real
  f_value = value
  ! in: f_print_limit_warning 0D_NOT_logical
  if (c_associated(print_limit_warning)) then
    call c_f_pointer(print_limit_warning, f_print_limit_warning_ptr)
    f_print_limit_warning_native = f_print_limit_warning_ptr
    f_print_limit_warning_native_ptr => f_print_limit_warning_native
  else
    f_print_limit_warning_native_ptr => null()
  endif
  call tao_set_var_model_value(f_var, f_value, f_print_limit_warning_native_ptr)

end subroutine
subroutine fortran_tao_set_var_useit_opt () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_set_var_useit_opt()

end subroutine
subroutine fortran_tao_set_wave_cmd (who, value_str, err) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: who
  character(len=4096), target :: f_who
  character(kind=c_char), pointer :: f_who_ptr(:)
  type(c_ptr), intent(in), value :: value_str
  character(len=4096), target :: f_value_str
  character(kind=c_char), pointer :: f_value_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_who 0D_NOT_character
  if (.not. c_associated(who)) return
  call c_f_pointer(who, f_who_ptr, [huge(0)])
  call to_f_str(f_who_ptr, f_who)
  ! in: f_value_str 0D_NOT_character
  if (.not. c_associated(value_str)) return
  call c_f_pointer(value_str, f_value_str_ptr, [huge(0)])
  call to_f_str(f_value_str_ptr, f_value_str)
  call tao_set_wave_cmd(f_who, f_value_str, f_err)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_set_z_tune_cmd (branch_str, q_str, delta_input) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: branch_str
  character(len=4096), target :: f_branch_str
  character(kind=c_char), pointer :: f_branch_str_ptr(:)
  type(c_ptr), intent(in), value :: q_str
  character(len=4096), target :: f_q_str
  character(kind=c_char), pointer :: f_q_str_ptr(:)
  logical(c_bool) :: delta_input  ! 0D_NOT_logical
  logical :: f_delta_input
  ! ** End of parameters **
  ! in: f_branch_str 0D_NOT_character
  if (.not. c_associated(branch_str)) return
  call c_f_pointer(branch_str, f_branch_str_ptr, [huge(0)])
  call to_f_str(f_branch_str_ptr, f_branch_str)
  ! in: f_q_str 0D_NOT_character
  if (.not. c_associated(q_str)) return
  call c_f_pointer(q_str, f_q_str_ptr, [huge(0)])
  call to_f_str(f_q_str_ptr, f_q_str)
  ! in: f_delta_input 0D_NOT_logical
  f_delta_input = delta_input
  call tao_set_z_tune_cmd(f_branch_str, f_q_str, f_delta_input)

end subroutine
subroutine fortran_tao_setup_key_table () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_setup_key_table()

end subroutine
subroutine fortran_tao_shape_init (shape, err, print_err) bind(c)

  use tao_struct, only: tao_ele_shape_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: shape  ! 0D_NOT_type
  type(tao_ele_shape_struct), pointer :: f_shape
  ! ** End of parameters **
  ! inout: f_shape 0D_NOT_type
  if (.not. c_associated(shape)) return
  call c_f_pointer(shape, f_shape)
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call tao_shape_init(f_shape, f_err, f_print_err_native_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_show_cmd (what) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: what
  character(len=4096), target :: f_what
  character(kind=c_char), pointer :: f_what_ptr(:)
  ! ** End of parameters **
  ! in: f_what 0D_NOT_character
  if (.not. c_associated(what)) return
  call c_f_pointer(what, f_what_ptr, [huge(0)])
  call to_f_str(f_what_ptr, f_what)
  call tao_show_cmd(f_what)

end subroutine
subroutine fortran_tao_show_constraints (iunit, form) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: iunit  ! 0D_NOT_integer
  integer :: f_iunit
  type(c_ptr), intent(in), value :: form
  character(len=4096), target :: f_form
  character(kind=c_char), pointer :: f_form_ptr(:)
  ! ** End of parameters **
  ! in: f_iunit 0D_NOT_integer
  f_iunit = iunit
  ! in: f_form 0D_NOT_character
  if (.not. c_associated(form)) return
  call c_f_pointer(form, f_form_ptr, [huge(0)])
  call to_f_str(f_form_ptr, f_form)
  call tao_show_constraints(f_iunit, f_form)

end subroutine
subroutine fortran_tao_single_mode (char_) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: char_
  character(len=4096), target :: f_char
  character(kind=c_char), pointer :: f_char_ptr(:)
  ! ** End of parameters **
  ! in: f_char 0D_NOT_character
  if (.not. c_associated(char_)) return
  call c_f_pointer(char_, f_char_ptr, [huge(0)])
  call to_f_str(f_char_ptr, f_char)
  call tao_single_mode(f_char)

end subroutine
subroutine fortran_tao_single_track (tao_lat, calc_ok, ix_branch, print_err) bind(c)

  use tao_struct, only: tao_lattice_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: tao_lat  ! 0D_NOT_type
  type(tao_lattice_struct), pointer :: f_tao_lat
  integer(c_int) :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  type(c_ptr), intent(in), value :: print_err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_print_err
  logical, target :: f_print_err_native
  logical, pointer :: f_print_err_native_ptr
  logical(c_bool), pointer :: f_print_err_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: calc_ok  ! 0D_NOT_logical
  logical :: f_calc_ok
  logical(c_bool), pointer :: f_calc_ok_ptr
  ! ** End of parameters **
  ! in: f_tao_lat 0D_NOT_type
  if (.not. c_associated(tao_lat)) return
  call c_f_pointer(tao_lat, f_tao_lat)
  ! in: f_ix_branch 0D_NOT_integer
  f_ix_branch = ix_branch
  ! in: f_print_err 0D_NOT_logical
  if (c_associated(print_err)) then
    call c_f_pointer(print_err, f_print_err_ptr)
    f_print_err_native = f_print_err_ptr
    f_print_err_native_ptr => f_print_err_native
  else
    f_print_err_native_ptr => null()
  endif
  call tao_single_track(f_tao_lat, f_calc_ok, f_ix_branch, f_print_err_native_ptr)

  ! out: f_calc_ok 0D_NOT_logical
  call c_f_pointer(calc_ok, f_calc_ok_ptr)
  f_calc_ok_ptr = f_calc_ok
end subroutine
subroutine fortran_tao_spin_matrices_calc_needed (data_type, data_source, do_calc) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: do_calc  ! 0D_NOT_logical
  logical :: f_do_calc
  logical(c_bool), pointer :: f_do_calc_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: data_source
  character(len=4096), target :: f_data_source
  character(kind=c_char), pointer :: f_data_source_ptr(:)
  ! ** End of parameters **
  ! inout: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! inout: f_data_source 0D_NOT_character
  if (.not. c_associated(data_source)) return
  call c_f_pointer(data_source, f_data_source_ptr, [huge(0)])
  call to_f_str(f_data_source_ptr, f_data_source)
  f_do_calc = tao_spin_matrices_calc_needed(f_data_type, f_data_source)

  ! inout: f_data_type 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_data_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_do_calc 0D_NOT_logical
  call c_f_pointer(do_calc, f_do_calc_ptr)
  f_do_calc_ptr = f_do_calc
end subroutine
subroutine fortran_tao_spin_tracking_turn_on () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_spin_tracking_turn_on()

end subroutine
subroutine fortran_tao_split_component (comp_str, comp, err) bind(c)

  use tao_struct, only: tao_data_var_component_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: comp_str
  character(len=4096), target :: f_comp_str
  character(kind=c_char), pointer :: f_comp_str_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: comp
  type(tao_data_var_component_struct_container_alloc), pointer :: f_comp
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_comp_str 0D_NOT_character
  if (.not. c_associated(comp_str)) return
  call c_f_pointer(comp_str, f_comp_str_ptr, [huge(0)])
  call to_f_str(f_comp_str_ptr, f_comp_str)
  !! container type array (1D_ALLOC_type)
  if (c_associated(comp))   call c_f_pointer(comp, f_comp)
  call tao_split_component(f_comp_str, f_comp%data, f_err)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_srdt_calc_needed (data_type, data_source, do_srdt) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: do_srdt  ! 0D_NOT_integer
  integer :: f_do_srdt
  integer(c_int), pointer :: f_do_srdt_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: data_type
  character(len=4096), target :: f_data_type
  character(kind=c_char), pointer :: f_data_type_ptr(:)
  type(c_ptr), intent(in), value :: data_source
  character(len=4096), target :: f_data_source
  character(kind=c_char), pointer :: f_data_source_ptr(:)
  ! ** End of parameters **
  ! inout: f_data_type 0D_NOT_character
  if (.not. c_associated(data_type)) return
  call c_f_pointer(data_type, f_data_type_ptr, [huge(0)])
  call to_f_str(f_data_type_ptr, f_data_type)
  ! inout: f_data_source 0D_NOT_character
  if (.not. c_associated(data_source)) return
  call c_f_pointer(data_source, f_data_source_ptr, [huge(0)])
  call to_f_str(f_data_source_ptr, f_data_source)
  f_do_srdt = tao_srdt_calc_needed(f_data_type, f_data_source)

  ! inout: f_data_type 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_data_source 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_do_srdt 0D_NOT_integer
  call c_f_pointer(do_srdt, f_do_srdt_ptr)
  f_do_srdt_ptr = f_do_srdt
end subroutine
subroutine fortran_tao_subin_uni_number (name_in, ix_uni, name_out, ok) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: name_in
  character(len=4096), target :: f_name_in
  character(kind=c_char), pointer :: f_name_in_ptr(:)
  integer(c_int) :: ix_uni  ! 0D_NOT_integer
  integer :: f_ix_uni
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: name_out
  character(len=4096), target :: f_name_out
  character(kind=c_char), pointer :: f_name_out_ptr(:)
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  ! ** End of parameters **
  ! in: f_name_in 0D_NOT_character
  if (.not. c_associated(name_in)) return
  call c_f_pointer(name_in, f_name_in_ptr, [huge(0)])
  call to_f_str(f_name_in_ptr, f_name_in)
  ! in: f_ix_uni 0D_NOT_integer
  f_ix_uni = ix_uni
  f_ok = tao_subin_uni_number(f_name_in, f_ix_uni, f_name_out)

  ! out: f_name_out 0D_NOT_character
  call c_f_pointer(name_out, f_name_out_ptr, [len_trim(f_name_out) + 1]) ! output-only string
  call to_c_str(f_name_out, f_name_out_ptr)
  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
end subroutine
subroutine fortran_tao_svd_optimizer (abort) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: abort  ! 0D_NOT_logical
  logical :: f_abort
  logical(c_bool), pointer :: f_abort_ptr
  ! ** End of parameters **
  call tao_svd_optimizer(f_abort)

  ! out: f_abort 0D_NOT_logical
  call c_f_pointer(abort, f_abort_ptr)
  f_abort_ptr = f_abort
end subroutine
subroutine fortran_tao_symbol_import_from_lat (lat) bind(c)

  use bmad_struct, only: lat_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: lat  ! 0D_NOT_type
  type(lat_struct), pointer :: f_lat
  ! ** End of parameters **
  ! inout: f_lat 0D_NOT_type
  if (.not. c_associated(lat)) return
  call c_f_pointer(lat, f_lat)
  call tao_symbol_import_from_lat(f_lat)

end subroutine
subroutine fortran_tao_taper_cmd (except, uni_names) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: except
  character(len=4096), target :: f_except
  character(kind=c_char), pointer :: f_except_ptr(:)
  type(c_ptr), intent(in), value :: uni_names
  character(len=4096), target :: f_uni_names
  character(kind=c_char), pointer :: f_uni_names_ptr(:)
  ! ** End of parameters **
  ! in: f_except 0D_NOT_character
  if (.not. c_associated(except)) return
  call c_f_pointer(except, f_except_ptr, [huge(0)])
  call to_f_str(f_except_ptr, f_except)
  ! in: f_uni_names 0D_NOT_character
  if (.not. c_associated(uni_names)) return
  call c_f_pointer(uni_names, f_uni_names_ptr, [huge(0)])
  call to_f_str(f_uni_names_ptr, f_uni_names)
  call tao_taper_cmd(f_except, f_uni_names)

end subroutine
subroutine fortran_tao_to_change_number (num_str, n_size, change_number, abs_or_rel, err) &
    bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: num_str
  character(len=4096), target :: f_num_str
  character(kind=c_char), pointer :: f_num_str_ptr(:)
  type(c_ptr), intent(in), value :: n_size  ! 0D_NOT_integer
  integer(c_int) :: f_n_size
  integer(c_int), pointer :: f_n_size_ptr
  type(c_ptr), intent(in), value :: change_number
  type(real_container_alloc), pointer :: f_change_number
  type(c_ptr), intent(in), value :: abs_or_rel
  character(len=4096), target :: f_abs_or_rel
  character(kind=c_char), pointer :: f_abs_or_rel_ptr(:)
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err
  logical, target :: f_err_native
  logical, pointer :: f_err_native_ptr
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! inout: f_num_str 0D_NOT_character
  if (.not. c_associated(num_str)) return
  call c_f_pointer(num_str, f_num_str_ptr, [huge(0)])
  call to_f_str(f_num_str_ptr, f_num_str)
  ! inout: f_n_size 0D_NOT_integer
  if (c_associated(n_size)) then
    call c_f_pointer(n_size, f_n_size_ptr)
  else
    f_n_size_ptr => null()
  endif
  !! container general array (1D_ALLOC_real)
  if (c_associated(change_number))   call c_f_pointer(change_number, f_change_number)
  ! inout: f_abs_or_rel 0D_NOT_character
  if (.not. c_associated(abs_or_rel)) return
  call c_f_pointer(abs_or_rel, f_abs_or_rel_ptr, [huge(0)])
  call to_f_str(f_abs_or_rel_ptr, f_abs_or_rel)
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_native = f_err_ptr
    f_err_native_ptr => f_err_native
  else
    f_err_native_ptr => null()
  endif
  call tao_to_change_number(f_num_str, f_n_size_ptr, f_change_number%data, f_abs_or_rel, &
      f_err_native_ptr)

  ! inout: f_num_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_n_size 0D_NOT_integer
  ! no output conversion for f_n_size
  ! inout: f_abs_or_rel 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_ptr = f_err_native
  else
    ! f_err unset
  endif
end subroutine
subroutine fortran_tao_to_int (str, i_int, err) bind(c)

  implicit none
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: str
  character(len=4096), target :: f_str
  character(kind=c_char), pointer :: f_str_ptr(:)
  type(c_ptr), intent(in), value :: i_int  ! 0D_NOT_integer
  integer(c_int) :: f_i_int
  integer(c_int), pointer :: f_i_int_ptr
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err
  logical, target :: f_err_native
  logical, pointer :: f_err_native_ptr
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! inout: f_str 0D_NOT_character
  if (.not. c_associated(str)) return
  call c_f_pointer(str, f_str_ptr, [huge(0)])
  call to_f_str(f_str_ptr, f_str)
  ! inout: f_i_int 0D_NOT_integer
  if (c_associated(i_int)) then
    call c_f_pointer(i_int, f_i_int_ptr)
  else
    f_i_int_ptr => null()
  endif
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_native = f_err_ptr
    f_err_native_ptr => f_err_native
  else
    f_err_native_ptr => null()
  endif
  call tao_to_int(f_str, f_i_int_ptr, f_err_native_ptr)

  ! inout: f_str 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! inout: f_i_int 0D_NOT_integer
  ! no output conversion for f_i_int
  ! inout: f_err 0D_NOT_logical
  if (c_associated(err)) then
    call c_f_pointer(err, f_err_ptr)
    f_err_ptr = f_err_native
  else
    ! f_err unset
  endif
end subroutine
subroutine fortran_tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, &
    datum) bind(c)

  use bmad_struct, only: bpm_phase_coupling_struct, ele_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  ! ** Out parameters **
  type(c_ptr), value :: bpm_data  ! 0D_NOT_type
  type(bpm_phase_coupling_struct), pointer :: f_bpm_data
  type(c_ptr), intent(in), value :: valid_value  ! 0D_NOT_logical
  logical :: f_valid_value
  logical(c_bool), pointer :: f_valid_value_ptr
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: why_invalid
  character(len=4096), target :: f_why_invalid
  character(kind=c_char), pointer :: f_why_invalid_ptr(:)
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** End of parameters **
  ! in: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! out: f_bpm_data 0D_NOT_type
  if (.not. c_associated(bpm_data)) return
  call c_f_pointer(bpm_data, f_bpm_data)
  ! inout: f_why_invalid 0D_NOT_character
  if (.not. c_associated(why_invalid)) return
  call c_f_pointer(why_invalid, f_why_invalid_ptr, [huge(0)])
  call to_f_str(f_why_invalid_ptr, f_why_invalid)
  ! inout: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  call tao_to_phase_and_coupling_reading(f_ele, f_bpm_data, f_valid_value, f_why_invalid, &
      f_datum)

  ! out: f_bpm_data 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
  ! out: f_valid_value 0D_NOT_logical
  call c_f_pointer(valid_value, f_valid_value_ptr)
  f_valid_value_ptr = f_valid_value
  ! inout: f_why_invalid 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_tao_to_real (expression, value, err_flag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: expression
  character(len=4096), target :: f_expression
  character(kind=c_char), pointer :: f_expression_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: value  ! 0D_NOT_real
  real(rp) :: f_value
  real(c_double), pointer :: f_value_ptr
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical :: f_err_flag
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_expression 0D_NOT_character
  if (.not. c_associated(expression)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = .true.
    return
  endif
  call c_f_pointer(expression, f_expression_ptr, [huge(0)])
  call to_f_str(f_expression_ptr, f_expression)
  call tao_to_real(f_expression, f_value, f_err_flag)

  ! out: f_value 0D_NOT_real
  call c_f_pointer(value, f_value_ptr)
  f_value_ptr = f_value
  ! out: f_err_flag 0D_NOT_logical
  call c_f_pointer(err_flag, f_err_flag_ptr)
  f_err_flag_ptr = f_err_flag
end subroutine
subroutine fortran_tao_too_many_particles_lost (beam, no_beam) bind(c)

  use bmad_struct, only: beam_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: no_beam  ! 0D_NOT_logical
  logical :: f_no_beam
  logical(c_bool), pointer :: f_no_beam_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: beam  ! 0D_NOT_type
  type(beam_struct), pointer :: f_beam
  ! ** End of parameters **
  ! inout: f_beam 0D_NOT_type
  if (.not. c_associated(beam)) return
  call c_f_pointer(beam, f_beam)
  f_no_beam = tao_too_many_particles_lost(f_beam)

  ! out: f_no_beam 0D_NOT_logical
  call c_f_pointer(no_beam, f_no_beam_ptr)
  f_no_beam_ptr = f_no_beam
end subroutine
subroutine fortran_tao_top10_derivative_print () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_top10_derivative_print()

end subroutine
subroutine fortran_tao_top10_merit_categories_print (iunit) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: iunit  ! 0D_NOT_integer
  integer :: f_iunit
  ! ** End of parameters **
  ! in: f_iunit 0D_NOT_integer
  f_iunit = iunit
  call tao_top10_merit_categories_print(f_iunit)

end subroutine
subroutine fortran_tao_top_level (command, errcode) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: command
  character(len=4096), target :: f_command
  character(kind=c_char), pointer :: f_command_ptr(:)
  character(len=4096), pointer :: f_command_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: errcode  ! 0D_NOT_integer
  integer :: f_errcode
  integer(c_int), pointer :: f_errcode_ptr
  ! ** End of parameters **
  ! in: f_command 0D_NOT_character
  if (c_associated(command)) then
    call c_f_pointer(command, f_command_ptr, [huge(0)])
    call to_f_str(f_command_ptr, f_command)
    f_command_call_ptr => f_command
  else
    f_command_call_ptr => null()
  endif
  call tao_top_level(f_command_call_ptr, f_errcode)

  ! out: f_errcode 0D_NOT_integer
  ! no output conversion for f_errcode
end subroutine
subroutine fortran_tao_tracking_ele_index (ele, datum, ix_branch, ix_ele) bind(c)

  use bmad_struct, only: ele_struct
  use tao_struct, only: tao_data_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: ele  ! 0D_PTR_type
  type(ele_struct), pointer :: f_ele
  type(c_ptr), value :: datum  ! 0D_NOT_type
  type(tao_data_struct), pointer :: f_datum
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_branch  ! 0D_NOT_integer
  integer :: f_ix_branch
  integer(c_int), pointer :: f_ix_branch_ptr
  type(c_ptr), intent(in), value :: ix_ele  ! 0D_NOT_integer
  integer :: f_ix_ele
  integer(c_int), pointer :: f_ix_ele_ptr
  ! ** End of parameters **
  ! in: f_ele 0D_PTR_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_datum 0D_NOT_type
  if (.not. c_associated(datum)) return
  call c_f_pointer(datum, f_datum)
  f_ix_ele = tao_tracking_ele_index(f_ele, f_datum, f_ix_branch)

  ! out: f_ix_branch 0D_NOT_integer
  ! no output conversion for f_ix_branch
  ! out: f_ix_ele 0D_NOT_integer
  call c_f_pointer(ix_ele, f_ix_ele_ptr)
  f_ix_ele_ptr = f_ix_ele
end subroutine
subroutine fortran_tao_turn_on_special_calcs_if_needed_for_plotting () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_turn_on_special_calcs_if_needed_for_plotting()

end subroutine
subroutine fortran_tao_uni_atsign_index (string, ix_amp) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: string
  character(len=4096), target :: f_string
  character(kind=c_char), pointer :: f_string_ptr(:)
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ix_amp  ! 0D_NOT_integer
  integer :: f_ix_amp
  integer(c_int), pointer :: f_ix_amp_ptr
  ! ** End of parameters **
  ! in: f_string 0D_NOT_character
  if (.not. c_associated(string)) return
  call c_f_pointer(string, f_string_ptr, [huge(0)])
  call to_f_str(f_string_ptr, f_string)
  f_ix_amp = tao_uni_atsign_index(f_string)

  ! out: f_ix_amp 0D_NOT_integer
  call c_f_pointer(ix_amp, f_ix_amp_ptr)
  f_ix_amp_ptr = f_ix_amp
end subroutine
subroutine fortran_tao_universe_index (i_uni, neg2_to_default, i_this_uni) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: i_uni  ! 0D_NOT_integer
  integer :: f_i_uni
  type(c_ptr), intent(in), value :: neg2_to_default  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_neg2_to_default
  logical, target :: f_neg2_to_default_native
  logical, pointer :: f_neg2_to_default_native_ptr
  logical(c_bool), pointer :: f_neg2_to_default_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: i_this_uni  ! 0D_NOT_integer
  integer :: f_i_this_uni
  integer(c_int), pointer :: f_i_this_uni_ptr
  ! ** End of parameters **
  ! in: f_i_uni 0D_NOT_integer
  f_i_uni = i_uni
  ! in: f_neg2_to_default 0D_NOT_logical
  if (c_associated(neg2_to_default)) then
    call c_f_pointer(neg2_to_default, f_neg2_to_default_ptr)
    f_neg2_to_default_native = f_neg2_to_default_ptr
    f_neg2_to_default_native_ptr => f_neg2_to_default_native
  else
    f_neg2_to_default_native_ptr => null()
  endif
  f_i_this_uni = tao_universe_index(f_i_uni, f_neg2_to_default_native_ptr)

  ! out: f_i_this_uni 0D_NOT_integer
  call c_f_pointer(i_this_uni, f_i_this_uni_ptr)
  f_i_this_uni_ptr = f_i_this_uni
end subroutine
subroutine fortran_tao_use_data (action, data_name) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: action
  character(len=4096), target :: f_action
  character(kind=c_char), pointer :: f_action_ptr(:)
  type(c_ptr), intent(in), value :: data_name
  character(len=4096), target :: f_data_name
  character(kind=c_char), pointer :: f_data_name_ptr(:)
  ! ** End of parameters **
  ! in: f_action 0D_NOT_character
  if (.not. c_associated(action)) return
  call c_f_pointer(action, f_action_ptr, [huge(0)])
  call to_f_str(f_action_ptr, f_action)
  ! in: f_data_name 0D_NOT_character
  if (.not. c_associated(data_name)) return
  call c_f_pointer(data_name, f_data_name_ptr, [huge(0)])
  call to_f_str(f_data_name_ptr, f_data_name)
  call tao_use_data(f_action, f_data_name)

end subroutine
subroutine fortran_tao_use_var (action, var_name) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: action
  character(len=4096), target :: f_action
  character(kind=c_char), pointer :: f_action_ptr(:)
  type(c_ptr), intent(in), value :: var_name
  character(len=4096), target :: f_var_name
  character(kind=c_char), pointer :: f_var_name_ptr(:)
  ! ** End of parameters **
  ! in: f_action 0D_NOT_character
  if (.not. c_associated(action)) return
  call c_f_pointer(action, f_action_ptr, [huge(0)])
  call to_f_str(f_action_ptr, f_action)
  ! in: f_var_name 0D_NOT_character
  if (.not. c_associated(var_name)) return
  call c_f_pointer(var_name, f_var_name_ptr, [huge(0)])
  call to_f_str(f_var_name_ptr, f_var_name)
  call tao_use_var(f_action, f_var_name)

end subroutine
subroutine fortran_tao_user_is_terminating_optimization (is_terminating) bind(c)

  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: is_terminating  ! 0D_NOT_logical
  logical :: f_is_terminating
  logical(c_bool), pointer :: f_is_terminating_ptr
  ! ** End of parameters **
  f_is_terminating = tao_user_is_terminating_optimization()

  ! out: f_is_terminating 0D_NOT_logical
  call c_f_pointer(is_terminating, f_is_terminating_ptr)
  f_is_terminating_ptr = f_is_terminating
end subroutine
subroutine fortran_tao_var1_name (var, var1_name) bind(c)

  use tao_struct, only: tao_var_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: var  ! 0D_NOT_type
  type(tao_var_struct), pointer :: f_var
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var1_name
  character(len=4096), target :: f_var1_name
  character(kind=c_char), pointer :: f_var1_name_ptr(:)
  ! ** End of parameters **
  ! in: f_var 0D_NOT_type
  if (.not. c_associated(var)) return
  call c_f_pointer(var, f_var)
  f_var1_name = tao_var1_name(f_var)

  ! out: f_var1_name 0D_NOT_character
  call c_f_pointer(var1_name, f_var1_name_ptr, [len_trim(f_var1_name) + 1]) ! output-only string
  call to_c_str(f_var1_name, f_var1_name_ptr)
end subroutine
subroutine fortran_tao_var_attrib_name (var, var_attrib_name) bind(c)

  use tao_struct, only: tao_var_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: var  ! 0D_NOT_type
  type(tao_var_struct), pointer :: f_var
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var_attrib_name
  character(len=4096), target :: f_var_attrib_name
  character(kind=c_char), pointer :: f_var_attrib_name_ptr(:)
  ! ** End of parameters **
  ! in: f_var 0D_NOT_type
  if (.not. c_associated(var)) return
  call c_f_pointer(var, f_var)
  f_var_attrib_name = tao_var_attrib_name(f_var)

  ! out: f_var_attrib_name 0D_NOT_character
  call c_f_pointer(var_attrib_name, f_var_attrib_name_ptr, [len_trim(f_var_attrib_name) + 1]) ! output-only string
  call to_c_str(f_var_attrib_name, f_var_attrib_name_ptr)
end subroutine
subroutine fortran_tao_var_check (eles, attribute, silent) bind(c)

  use bmad_struct, only: ele_pointer_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: eles
  type(ele_pointer_struct_container_alloc), pointer :: f_eles
  type(c_ptr), intent(in), value :: attribute
  character(len=4096), target :: f_attribute
  character(kind=c_char), pointer :: f_attribute_ptr(:)
  logical(c_bool) :: silent  ! 0D_NOT_logical
  logical :: f_silent
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(eles))   call c_f_pointer(eles, f_eles)
  ! in: f_attribute 0D_NOT_character
  if (.not. c_associated(attribute)) return
  call c_f_pointer(attribute, f_attribute_ptr, [huge(0)])
  call to_f_str(f_attribute_ptr, f_attribute)
  ! in: f_silent 0D_NOT_logical
  f_silent = silent
  call tao_var_check(f_eles%data, f_attribute, f_silent)

end subroutine
subroutine fortran_tao_var_repoint () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_var_repoint()

end subroutine
subroutine fortran_tao_var_target_calc () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_var_target_calc()

end subroutine
subroutine fortran_tao_var_useit_plot_calc (graph, var) bind(c)

  use tao_struct, only: tao_graph_struct, tao_var_struct
  implicit none
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: var
  type(tao_var_struct_container_alloc), pointer :: f_var
  ! ** Inout parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  ! ** End of parameters **
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  !! container type array (1D_ALLOC_type)
  if (c_associated(var))   call c_f_pointer(var, f_var)
  call tao_var_useit_plot_calc(f_graph, f_var%data)

end subroutine
subroutine fortran_tao_var_write (out_file, show_good_opt_only, tao_format) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: out_file
  character(len=4096), target :: f_out_file
  character(kind=c_char), pointer :: f_out_file_ptr(:)
  type(c_ptr), intent(in), value :: show_good_opt_only  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_show_good_opt_only
  logical, target :: f_show_good_opt_only_native
  logical, pointer :: f_show_good_opt_only_native_ptr
  logical(c_bool), pointer :: f_show_good_opt_only_ptr
  type(c_ptr), intent(in), value :: tao_format  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_tao_format
  logical, target :: f_tao_format_native
  logical, pointer :: f_tao_format_native_ptr
  logical(c_bool), pointer :: f_tao_format_ptr
  ! ** End of parameters **
  ! in: f_out_file 0D_NOT_character
  if (.not. c_associated(out_file)) return
  call c_f_pointer(out_file, f_out_file_ptr, [huge(0)])
  call to_f_str(f_out_file_ptr, f_out_file)
  ! in: f_show_good_opt_only 0D_NOT_logical
  if (c_associated(show_good_opt_only)) then
    call c_f_pointer(show_good_opt_only, f_show_good_opt_only_ptr)
    f_show_good_opt_only_native = f_show_good_opt_only_ptr
    f_show_good_opt_only_native_ptr => f_show_good_opt_only_native
  else
    f_show_good_opt_only_native_ptr => null()
  endif
  ! in: f_tao_format 0D_NOT_logical
  if (c_associated(tao_format)) then
    call c_f_pointer(tao_format, f_tao_format_ptr)
    f_tao_format_native = f_tao_format_ptr
    f_tao_format_native_ptr => f_tao_format_native
  else
    f_tao_format_native_ptr => null()
  endif
  call tao_var_write(f_out_file, f_show_good_opt_only_native_ptr, f_tao_format_native_ptr)

end subroutine
subroutine fortran_tao_veto_vars_with_zero_dmodel () bind(c)

  implicit none
  ! ** End of parameters **
  call tao_veto_vars_with_zero_dmodel()

end subroutine
subroutine fortran_tao_wave_analysis (plot) bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  ! ** End of parameters **
  ! inout: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  call tao_wave_analysis(f_plot)

end subroutine
subroutine fortran_tao_wave_cmd (curve_name, plot_place, err_flag) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: curve_name
  character(len=4096), target :: f_curve_name
  character(kind=c_char), pointer :: f_curve_name_ptr(:)
  type(c_ptr), intent(in), value :: plot_place
  character(len=4096), target :: f_plot_place
  character(kind=c_char), pointer :: f_plot_place_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: err_flag  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_err_flag
  logical, target :: f_err_flag_native
  logical, pointer :: f_err_flag_native_ptr
  logical(c_bool), pointer :: f_err_flag_ptr
  ! ** End of parameters **
  ! in: f_curve_name 0D_NOT_character
  if (.not. c_associated(curve_name)) return
  call c_f_pointer(curve_name, f_curve_name_ptr, [huge(0)])
  call to_f_str(f_curve_name_ptr, f_curve_name)
  ! in: f_plot_place 0D_NOT_character
  if (.not. c_associated(plot_place)) return
  call c_f_pointer(plot_place, f_plot_place_ptr, [huge(0)])
  call to_f_str(f_plot_place_ptr, f_plot_place)
  ! inout: f_err_flag 0D_NOT_logical
  if (c_associated(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_native = f_err_flag_ptr
    f_err_flag_native_ptr => f_err_flag_native
  else
    f_err_flag_native_ptr => null()
  endif
  call tao_wave_cmd(f_curve_name, f_plot_place, f_err_flag_native_ptr)

  ! inout: f_err_flag 0D_NOT_logical
  if (c_associated(err_flag)) then
    call c_f_pointer(err_flag, f_err_flag_ptr)
    f_err_flag_ptr = f_err_flag_native
  else
    ! f_err_flag unset
  endif
end subroutine
subroutine fortran_tao_wave_fit (curve, ix1, n_dat, coef, rms, f1, f2, f3, f4) bind(c)

  use tao_struct, only: tao_curve_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: curve  ! 0D_NOT_type
  type(tao_curve_struct), pointer :: f_curve
  integer(c_int) :: ix1  ! 0D_NOT_integer
  integer :: f_ix1
  integer(c_int) :: n_dat  ! 0D_NOT_integer
  integer :: f_n_dat
  type(c_ptr), intent(in), value :: f1
  type(real_container_alloc), pointer :: f_f1
  type(c_ptr), intent(in), value :: f2
  type(real_container_alloc), pointer :: f_f2
  type(c_ptr), intent(in), value :: f3
  type(real_container_alloc), pointer :: f_f3
  type(c_ptr), intent(in), value :: f4
  type(real_container_alloc), pointer :: f_f4
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: coef
  type(real_container_alloc), pointer :: f_coef
  type(c_ptr), intent(in), value :: rms
  type(real_container_alloc), pointer :: f_rms
  ! ** End of parameters **
  ! in: f_curve 0D_NOT_type
  if (.not. c_associated(curve)) return
  call c_f_pointer(curve, f_curve)
  ! in: f_ix1 0D_NOT_integer
  f_ix1 = ix1
  ! in: f_n_dat 0D_NOT_integer
  f_n_dat = n_dat
  !! container general array (1D_ALLOC_real)
  if (c_associated(coef))   call c_f_pointer(coef, f_coef)
  !! container general array (1D_ALLOC_real)
  if (c_associated(rms))   call c_f_pointer(rms, f_rms)
  !! container general array (1D_ALLOC_real)
  if (c_associated(f1))   call c_f_pointer(f1, f_f1)
  !! container general array (1D_ALLOC_real)
  if (c_associated(f2))   call c_f_pointer(f2, f_f2)
  !! container general array (1D_ALLOC_real)
  if (c_associated(f3))   call c_f_pointer(f3, f_f3)
  !! container general array (1D_ALLOC_real)
  if (c_associated(f4))   call c_f_pointer(f4, f_f4)
  call tao_wave_fit(f_curve, f_ix1, f_n_dat, f_coef%data, f_rms%data, f_f1%data, f_f2%data, &
      f_f3%data, f_f4%data)

end subroutine
subroutine fortran_tao_write_cmd (what) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: what
  character(len=4096), target :: f_what
  character(kind=c_char), pointer :: f_what_ptr(:)
  ! ** End of parameters **
  ! in: f_what 0D_NOT_character
  if (.not. c_associated(what)) return
  call c_f_pointer(what, f_what_ptr, [huge(0)])
  call to_f_str(f_what_ptr, f_what)
  call tao_write_cmd(f_what)

end subroutine
subroutine fortran_tao_x_axis_cmd (where, what) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  type(c_ptr), intent(in), value :: what
  character(len=4096), target :: f_what
  character(kind=c_char), pointer :: f_what_ptr(:)
  ! ** End of parameters **
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! in: f_what 0D_NOT_character
  if (.not. c_associated(what)) return
  call c_f_pointer(what, f_what_ptr, [huge(0)])
  call to_f_str(f_what_ptr, f_what)
  call tao_x_axis_cmd(f_where, f_what)

end subroutine
subroutine fortran_tao_x_scale_cmd (where, x_min_in, x_max_in, err, include_wall, gang, exact, &
    turn_autoscale_off) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: where
  character(len=4096), target :: f_where
  character(kind=c_char), pointer :: f_where_ptr(:)
  real(c_double) :: x_min_in  ! 0D_NOT_real
  real(rp) :: f_x_min_in
  real(c_double) :: x_max_in  ! 0D_NOT_real
  real(rp) :: f_x_max_in
  type(c_ptr), intent(in), value :: include_wall  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_wall
  logical, target :: f_include_wall_native
  logical, pointer :: f_include_wall_native_ptr
  logical(c_bool), pointer :: f_include_wall_ptr
  type(c_ptr), intent(in), value :: gang
  character(len=4096), target :: f_gang
  character(kind=c_char), pointer :: f_gang_ptr(:)
  character(len=4096), pointer :: f_gang_call_ptr
  type(c_ptr), intent(in), value :: exact  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_exact
  logical, target :: f_exact_native
  logical, pointer :: f_exact_native_ptr
  logical(c_bool), pointer :: f_exact_ptr
  type(c_ptr), intent(in), value :: turn_autoscale_off  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_turn_autoscale_off
  logical, target :: f_turn_autoscale_off_native
  logical, pointer :: f_turn_autoscale_off_native_ptr
  logical(c_bool), pointer :: f_turn_autoscale_off_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: err  ! 0D_NOT_logical
  logical :: f_err
  logical(c_bool), pointer :: f_err_ptr
  ! ** End of parameters **
  ! in: f_where 0D_NOT_character
  if (.not. c_associated(where)) return
  call c_f_pointer(where, f_where_ptr, [huge(0)])
  call to_f_str(f_where_ptr, f_where)
  ! in: f_x_min_in 0D_NOT_real
  f_x_min_in = x_min_in
  ! in: f_x_max_in 0D_NOT_real
  f_x_max_in = x_max_in
  ! in: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_native = f_include_wall_ptr
    f_include_wall_native_ptr => f_include_wall_native
  else
    f_include_wall_native_ptr => null()
  endif
  ! in: f_gang 0D_NOT_character
  if (c_associated(gang)) then
    call c_f_pointer(gang, f_gang_ptr, [huge(0)])
    call to_f_str(f_gang_ptr, f_gang)
    f_gang_call_ptr => f_gang
  else
    f_gang_call_ptr => null()
  endif
  ! in: f_exact 0D_NOT_logical
  if (c_associated(exact)) then
    call c_f_pointer(exact, f_exact_ptr)
    f_exact_native = f_exact_ptr
    f_exact_native_ptr => f_exact_native
  else
    f_exact_native_ptr => null()
  endif
  ! in: f_turn_autoscale_off 0D_NOT_logical
  if (c_associated(turn_autoscale_off)) then
    call c_f_pointer(turn_autoscale_off, f_turn_autoscale_off_ptr)
    f_turn_autoscale_off_native = f_turn_autoscale_off_ptr
    f_turn_autoscale_off_native_ptr => f_turn_autoscale_off_native
  else
    f_turn_autoscale_off_native_ptr => null()
  endif
  call tao_x_scale_cmd(f_where, f_x_min_in, f_x_max_in, f_err, f_include_wall_native_ptr, &
      f_gang_call_ptr, f_exact_native_ptr, f_turn_autoscale_off_native_ptr)

  ! out: f_err 0D_NOT_logical
  call c_f_pointer(err, f_err_ptr)
  f_err_ptr = f_err
end subroutine
subroutine fortran_tao_x_scale_graph (graph, x_min, x_max, include_wall, have_scaled) bind(c)

  use tao_struct, only: tao_graph_struct
  implicit none
  ! ** Inout parameters **
  type(c_ptr), value :: graph  ! 0D_NOT_type
  type(tao_graph_struct), pointer :: f_graph
  type(c_ptr), intent(in), value :: x_min  ! 0D_NOT_real
  real(c_double) :: f_x_min
  real(c_double), pointer :: f_x_min_ptr
  type(c_ptr), intent(in), value :: x_max  ! 0D_NOT_real
  real(c_double) :: f_x_max
  real(c_double), pointer :: f_x_max_ptr
  type(c_ptr), intent(in), value :: include_wall  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_wall
  logical, target :: f_include_wall_native
  logical, pointer :: f_include_wall_native_ptr
  logical(c_bool), pointer :: f_include_wall_ptr
  type(c_ptr), intent(in), value :: have_scaled  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_have_scaled
  logical, target :: f_have_scaled_native
  logical, pointer :: f_have_scaled_native_ptr
  logical(c_bool), pointer :: f_have_scaled_ptr
  ! ** End of parameters **
  ! inout: f_graph 0D_NOT_type
  if (.not. c_associated(graph)) return
  call c_f_pointer(graph, f_graph)
  ! inout: f_x_min 0D_NOT_real
  if (c_associated(x_min)) then
    call c_f_pointer(x_min, f_x_min_ptr)
  else
    f_x_min_ptr => null()
  endif
  ! inout: f_x_max 0D_NOT_real
  if (c_associated(x_max)) then
    call c_f_pointer(x_max, f_x_max_ptr)
  else
    f_x_max_ptr => null()
  endif
  ! inout: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_native = f_include_wall_ptr
    f_include_wall_native_ptr => f_include_wall_native
  else
    f_include_wall_native_ptr => null()
  endif
  ! inout: f_have_scaled 0D_NOT_logical
  if (c_associated(have_scaled)) then
    call c_f_pointer(have_scaled, f_have_scaled_ptr)
    f_have_scaled_native = f_have_scaled_ptr
    f_have_scaled_native_ptr => f_have_scaled_native
  else
    f_have_scaled_native_ptr => null()
  endif
  call tao_x_scale_graph(f_graph, f_x_min_ptr, f_x_max_ptr, f_include_wall_native_ptr, &
      f_have_scaled_native_ptr)

  ! inout: f_x_min 0D_NOT_real
  ! no output conversion for f_x_min
  ! inout: f_x_max 0D_NOT_real
  ! no output conversion for f_x_max
  ! inout: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_ptr = f_include_wall_native
  else
    ! f_include_wall unset
  endif
  ! inout: f_have_scaled 0D_NOT_logical
  if (c_associated(have_scaled)) then
    call c_f_pointer(have_scaled, f_have_scaled_ptr)
    f_have_scaled_ptr = f_have_scaled_native
  else
    ! f_have_scaled unset
  endif
end subroutine
subroutine fortran_tao_x_scale_plot (plot, x_min_in, x_max_in, include_wall, gang, have_scaled) &
    bind(c)

  use tao_struct, only: tao_plot_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: plot  ! 0D_NOT_type
  type(tao_plot_struct), pointer :: f_plot
  real(c_double) :: x_min_in  ! 0D_NOT_real
  real(rp) :: f_x_min_in
  real(c_double) :: x_max_in  ! 0D_NOT_real
  real(rp) :: f_x_max_in
  type(c_ptr), intent(in), value :: include_wall  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_include_wall
  logical, target :: f_include_wall_native
  logical, pointer :: f_include_wall_native_ptr
  logical(c_bool), pointer :: f_include_wall_ptr
  type(c_ptr), intent(in), value :: gang
  character(len=4096), target :: f_gang
  character(kind=c_char), pointer :: f_gang_ptr(:)
  character(len=4096), pointer :: f_gang_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: have_scaled  ! 0D_NOT_logical
  logical :: f_have_scaled
  logical(c_bool), pointer :: f_have_scaled_ptr
  ! ** End of parameters **
  ! in: f_plot 0D_NOT_type
  if (.not. c_associated(plot)) return
  call c_f_pointer(plot, f_plot)
  ! in: f_x_min_in 0D_NOT_real
  f_x_min_in = x_min_in
  ! in: f_x_max_in 0D_NOT_real
  f_x_max_in = x_max_in
  ! in: f_include_wall 0D_NOT_logical
  if (c_associated(include_wall)) then
    call c_f_pointer(include_wall, f_include_wall_ptr)
    f_include_wall_native = f_include_wall_ptr
    f_include_wall_native_ptr => f_include_wall_native
  else
    f_include_wall_native_ptr => null()
  endif
  ! in: f_gang 0D_NOT_character
  if (c_associated(gang)) then
    call c_f_pointer(gang, f_gang_ptr, [huge(0)])
    call to_f_str(f_gang_ptr, f_gang)
    f_gang_call_ptr => f_gang
  else
    f_gang_call_ptr => null()
  endif
  call tao_x_scale_plot(f_plot, f_x_min_in, f_x_max_in, f_include_wall_native_ptr, &
      f_gang_call_ptr, f_have_scaled)

  ! out: f_have_scaled 0D_NOT_logical
  ! no output conversion for f_have_scaled
end subroutine

end module cppbmad_tao_routines
