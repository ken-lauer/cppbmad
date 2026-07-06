!------------------------------------------------------------------------------
! tao_c_hook_interface
!
! Hand-written Fortran trampolines that let a C (and, through cppbmad, Python)
! callable be installed into Tao's hook procedure pointers.
!
! For each supported Tao hook `tao_hook_<x>_def` this module provides:
!   * a module `type(c_funptr) :: cb_<x>` holding the C callback,
!   * a trampoline `tramp_<x>` matching the abstract interface `tao_hook_<x>_def`
!     that converts each argument to a C-callable form and forwards the call,
!   * a `bind(c)` registration entry point `tao_hook_register_<x>`.
!
! Shared C-ABI abstract interfaces and `cloc_*` address helpers live in
! hook_helpers. See plans/hook-plan.md.
!------------------------------------------------------------------------------

module tao_c_hook_interface

  use, intrinsic :: iso_c_binding
  use bmad
  use tao_interface
  use hook_helpers
  use hook_c_interfaces

  implicit none

  type(c_funptr) :: cb_lattice_calc = c_null_funptr
  type(c_funptr) :: cb_optimizer = c_null_funptr
  type(c_funptr) :: cb_merit_var = c_null_funptr
  type(c_funptr) :: cb_merit_data = c_null_funptr

contains

  subroutine tramp_lattice_calc(calc_ok)
    logical :: calc_ok
    procedure(ci_lattice_calc), pointer :: cproc
    integer(c_int) :: cok
    call c_f_procpointer(cb_lattice_calc, cproc)
    cok = merge(1_c_int, 0_c_int, calc_ok)
    call cproc(cok)
    calc_ok = (cok /= 0)
  end subroutine

  subroutine tao_hook_register_lattice_calc(cfun) bind(c, name='tao_hook_register_lattice_calc')
    type(c_funptr), value :: cfun
    cb_lattice_calc = cfun
    if (c_associated(cfun)) then
      tao_hook_lattice_calc_ptr => tramp_lattice_calc
    else
      tao_hook_lattice_calc_ptr => null()
    end if
  end subroutine

  subroutine tramp_optimizer(abort)
    logical :: abort
    procedure(ci_optimizer), pointer :: cproc
    integer(c_int) :: cabort
    call c_f_procpointer(cb_optimizer, cproc)
    ! abort is a hook output; Bmad passes it uninitialized, so default it to false --
    ! a callback returning None then leaves the optimizer running rather than
    ! inheriting a garbage abort flag.
    cabort = 0
    call cproc(cabort)
    abort = (cabort /= 0)
  end subroutine

  subroutine tao_hook_register_optimizer(cfun) bind(c, name='tao_hook_register_optimizer')
    type(c_funptr), value :: cfun
    cb_optimizer = cfun
    if (c_associated(cfun)) then
      tao_hook_optimizer_ptr => tramp_optimizer
    else
      tao_hook_optimizer_ptr => null()
    end if
  end subroutine

  subroutine tramp_merit_var(i_uni, j_var, var)
    integer, intent(in) :: i_uni, j_var
    type(tao_var_struct) :: var
    procedure(ci_merit_var), pointer :: cproc
    call c_f_procpointer(cb_merit_var, cproc)
    call cproc(int(i_uni, c_int), int(j_var, c_int), cloc_tao_var(var))
  end subroutine

  subroutine tao_hook_register_merit_var(cfun) bind(c, name='tao_hook_register_merit_var')
    type(c_funptr), value :: cfun
    cb_merit_var = cfun
    if (c_associated(cfun)) then
      tao_hook_merit_var_ptr => tramp_merit_var
    else
      tao_hook_merit_var_ptr => null()
    end if
  end subroutine

  subroutine tramp_merit_data(i_uni, j_data, data, valid_value_set)
    integer, intent(in) :: i_uni, j_data
    type(tao_data_struct) :: data
    logical :: valid_value_set
    procedure(ci_merit_data), pointer :: cproc
    integer(c_int) :: cvalid
    call c_f_procpointer(cb_merit_data, cproc)
    cvalid = merge(1_c_int, 0_c_int, valid_value_set)
    call cproc(int(i_uni, c_int), int(j_data, c_int), cloc_tao_data(data), cvalid)
    valid_value_set = (cvalid /= 0)
  end subroutine

  subroutine tao_hook_register_merit_data(cfun) bind(c, name='tao_hook_register_merit_data')
    type(c_funptr), value :: cfun
    cb_merit_data = cfun
    if (c_associated(cfun)) then
      tao_hook_merit_data_ptr => tramp_merit_data
    else
      tao_hook_merit_data_ptr => null()
    end if
  end subroutine

end module tao_c_hook_interface
