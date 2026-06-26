module cppbmad_cppbmad_extra_routines

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

use cppbmad_ele_misalignments_mod, only: set_ele_misalignments


use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

subroutine fortran_set_ele_misalignments (ele, x_offset, y_offset, z_offset, x_pitch, y_pitch, &
    tilt, check_free, ok) bind(c)

  use array_desc_mod
  use bmad_struct, only: ele_struct
  implicit none
  ! ** In parameters **
  real(c_double) :: x_offset  ! 0D_NOT_real
  real(rp) :: f_x_offset
  real(c_double) :: y_offset  ! 0D_NOT_real
  real(rp) :: f_y_offset
  real(c_double) :: z_offset  ! 0D_NOT_real
  real(rp) :: f_z_offset
  real(c_double) :: x_pitch  ! 0D_NOT_real
  real(rp) :: f_x_pitch
  real(c_double) :: y_pitch  ! 0D_NOT_real
  real(rp) :: f_y_pitch
  real(c_double) :: tilt  ! 0D_NOT_real
  real(rp) :: f_tilt
  type(c_ptr), intent(in), value :: check_free  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_check_free
  logical, target :: f_check_free_native
  logical, pointer :: f_check_free_native_ptr
  logical(c_bool), pointer :: f_check_free_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: ok  ! 0D_NOT_logical
  logical :: f_ok
  logical(c_bool), pointer :: f_ok_ptr
  ! ** Inout parameters **
  type(c_ptr), value :: ele  ! 0D_NOT_type
  type(ele_struct), pointer :: f_ele
  ! ** End of parameters **
  ! inout: f_ele 0D_NOT_type
  if (.not. c_associated(ele)) return
  call c_f_pointer(ele, f_ele)
  ! in: f_x_offset 0D_NOT_real
  f_x_offset = x_offset
  ! in: f_y_offset 0D_NOT_real
  f_y_offset = y_offset
  ! in: f_z_offset 0D_NOT_real
  f_z_offset = z_offset
  ! in: f_x_pitch 0D_NOT_real
  f_x_pitch = x_pitch
  ! in: f_y_pitch 0D_NOT_real
  f_y_pitch = y_pitch
  ! in: f_tilt 0D_NOT_real
  f_tilt = tilt
  ! in: f_check_free 0D_NOT_logical
  if (c_associated(check_free)) then
    call c_f_pointer(check_free, f_check_free_ptr)
    f_check_free_native = f_check_free_ptr
    f_check_free_native_ptr => f_check_free_native
  else
    f_check_free_native_ptr => null()
  endif
  f_ok = set_ele_misalignments(f_ele, f_x_offset, f_y_offset, f_z_offset, f_x_pitch, f_y_pitch, &
      f_tilt, f_check_free_native_ptr)

  ! out: f_ok 0D_NOT_logical
  call c_f_pointer(ok, f_ok_ptr)
  f_ok_ptr = f_ok
end subroutine

end module cppbmad_cppbmad_extra_routines
