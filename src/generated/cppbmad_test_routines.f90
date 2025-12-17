module cppbmad_cppbmad_test_routines

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

use cppbmad_test_mod, only: test_bunch_struct_array, test_bunch_struct_scalar, &
    test_character_scalar, test_complex_array, test_complex_scalar, test_integer8_array, &
    test_integer8_scalar, test_integer_array, test_integer_scalar, test_logical_array, &
    test_logical_scalar, test_real16_array, test_real16_scalar, test_real_array, &
    test_real_scalar


use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

subroutine fortran_test_integer_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int) :: val_in  ! 0D_NOT_integer
  integer :: f_val_in
  type(c_ptr), intent(in), value :: val_in_opt  ! 0D_NOT_integer
  integer(c_int) :: f_val_in_opt
  integer(c_int), pointer :: f_val_in_opt_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out  ! 0D_NOT_integer
  integer :: f_val_out
  integer(c_int), pointer :: f_val_out_ptr
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout  ! 0D_NOT_integer
  integer(c_int) :: f_val_inout
  integer(c_int), pointer :: f_val_inout_ptr
  type(c_ptr), intent(in), value :: val_inout_opt  ! 0D_NOT_integer
  integer(c_int) :: f_val_inout_opt
  integer(c_int), pointer :: f_val_inout_opt_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_integer
  f_val_in = val_in
  ! inout: f_val_inout 0D_NOT_integer
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
  else
    f_val_inout_ptr => null()
  endif
  ! in: f_val_in_opt 0D_NOT_integer
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr)
  else
    f_val_in_opt_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_integer
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
  else
    f_val_inout_opt_ptr => null()
  endif
  call test_integer_scalar(val_in=f_val_in, val_inout=f_val_inout_ptr, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt_ptr, val_inout_opt=f_val_inout_opt_ptr)

  ! inout: f_val_inout 0D_NOT_integer
  ! no output conversion for f_val_inout
  ! out: f_val_out 0D_NOT_integer
  call c_f_pointer(val_out, f_val_out_ptr)
  f_val_out_ptr = f_val_out
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_integer
  ! no output conversion for f_val_inout_opt
end subroutine
subroutine fortran_test_integer_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(integer_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(integer_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(integer_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(integer_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(integer_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container general array (1D_ALLOC_integer)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container general array (1D_ALLOC_integer)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_integer_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_integer8_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  integer(c_int64_t) :: val_in  ! 0D_NOT_integer8
  integer(8) :: f_val_in
  type(c_ptr), intent(in), value :: val_in_opt  ! 0D_NOT_integer8
  integer(c_int64_t) :: f_val_in_opt
  integer(c_int64_t), pointer :: f_val_in_opt_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out  ! 0D_NOT_integer8
  integer(8) :: f_val_out
  integer(c_int64_t), pointer :: f_val_out_ptr
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout  ! 0D_NOT_integer8
  integer(c_int64_t) :: f_val_inout
  integer(c_int64_t), pointer :: f_val_inout_ptr
  type(c_ptr), intent(in), value :: val_inout_opt  ! 0D_NOT_integer8
  integer(c_int64_t) :: f_val_inout_opt
  integer(c_int64_t), pointer :: f_val_inout_opt_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_integer8
  f_val_in = val_in
  ! inout: f_val_inout 0D_NOT_integer8
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
  else
    f_val_inout_ptr => null()
  endif
  ! in: f_val_in_opt 0D_NOT_integer8
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr)
  else
    f_val_in_opt_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_integer8
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
  else
    f_val_inout_opt_ptr => null()
  endif
  call test_integer8_scalar(val_in=f_val_in, val_inout=f_val_inout_ptr, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt_ptr, val_inout_opt=f_val_inout_opt_ptr)

  ! inout: f_val_inout 0D_NOT_integer8
  ! no output conversion for f_val_inout
  ! out: f_val_out 0D_NOT_integer8
  call c_f_pointer(val_out, f_val_out_ptr)
  f_val_out_ptr = f_val_out
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_integer8
  ! no output conversion for f_val_inout_opt
end subroutine
subroutine fortran_test_integer8_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(integer8_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(integer8_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(integer8_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(integer8_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(integer8_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container general array (1D_ALLOC_integer8)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container general array (1D_ALLOC_integer8)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container general array (1D_ALLOC_integer8)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container general array (1D_ALLOC_integer8)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container general array (1D_ALLOC_integer8)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_integer8_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_real_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  real(c_double) :: val_in  ! 0D_NOT_real
  real(rp) :: f_val_in
  type(c_ptr), intent(in), value :: val_in_opt  ! 0D_NOT_real
  real(c_double) :: f_val_in_opt
  real(c_double), pointer :: f_val_in_opt_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out  ! 0D_NOT_real
  real(rp) :: f_val_out
  real(c_double), pointer :: f_val_out_ptr
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout  ! 0D_NOT_real
  real(c_double) :: f_val_inout
  real(c_double), pointer :: f_val_inout_ptr
  type(c_ptr), intent(in), value :: val_inout_opt  ! 0D_NOT_real
  real(c_double) :: f_val_inout_opt
  real(c_double), pointer :: f_val_inout_opt_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_real
  f_val_in = val_in
  ! inout: f_val_inout 0D_NOT_real
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
  else
    f_val_inout_ptr => null()
  endif
  ! in: f_val_in_opt 0D_NOT_real
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr)
  else
    f_val_in_opt_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_real
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
  else
    f_val_inout_opt_ptr => null()
  endif
  call test_real_scalar(val_in=f_val_in, val_inout=f_val_inout_ptr, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt_ptr, val_inout_opt=f_val_inout_opt_ptr)

  ! inout: f_val_inout 0D_NOT_real
  ! no output conversion for f_val_inout
  ! out: f_val_out 0D_NOT_real
  call c_f_pointer(val_out, f_val_out_ptr)
  f_val_out_ptr = f_val_out
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_real
  ! no output conversion for f_val_inout_opt
end subroutine
subroutine fortran_test_real_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(real_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(real_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(real_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(real_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(real_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container general array (1D_ALLOC_real)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container general array (1D_ALLOC_real)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container general array (1D_ALLOC_real)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container general array (1D_ALLOC_real)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_real_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_real16_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  real(c_long_double) :: val_in  ! 0D_NOT_real16
  real(qp) :: f_val_in
  type(c_ptr), intent(in), value :: val_in_opt  ! 0D_NOT_real16
  real(c_long_double), pointer :: f_val_in_opt
  real(16), target :: f_val_in_opt_native
  real(16), pointer :: f_val_in_opt_native_ptr
  real(c_long_double), pointer :: f_val_in_opt_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out  ! 0D_NOT_real16
  real(qp) :: f_val_out
  real(c_long_double), pointer :: f_val_out_ptr
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout  ! 0D_NOT_real16
  real(c_long_double), pointer :: f_val_inout
  real(16), target :: f_val_inout_native
  real(16), pointer :: f_val_inout_native_ptr
  real(c_long_double), pointer :: f_val_inout_ptr
  type(c_ptr), intent(in), value :: val_inout_opt  ! 0D_NOT_real16
  real(c_long_double), pointer :: f_val_inout_opt
  real(16), target :: f_val_inout_opt_native
  real(16), pointer :: f_val_inout_opt_native_ptr
  real(c_long_double), pointer :: f_val_inout_opt_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_real16
  f_val_in = val_in
  ! inout: f_val_inout 0D_NOT_real16
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
    f_val_inout_native = f_val_inout_ptr
    f_val_inout_native_ptr => f_val_inout_native
  else
    f_val_inout_native_ptr => null()
  endif
  ! in: f_val_in_opt 0D_NOT_real16
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr)
    f_val_in_opt_native = f_val_in_opt_ptr
    f_val_in_opt_native_ptr => f_val_in_opt_native
  else
    f_val_in_opt_native_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_real16
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
    f_val_inout_opt_native = f_val_inout_opt_ptr
    f_val_inout_opt_native_ptr => f_val_inout_opt_native
  else
    f_val_inout_opt_native_ptr => null()
  endif
  call test_real16_scalar(val_in=f_val_in, val_inout=f_val_inout_native_ptr, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt_native_ptr, &
      val_inout_opt=f_val_inout_opt_native_ptr)

  ! inout: f_val_inout 0D_NOT_real16
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
    f_val_inout_ptr = f_val_inout_native
  else
    ! f_val_inout unset
  endif
  ! out: f_val_out 0D_NOT_real16
  call c_f_pointer(val_out, f_val_out_ptr)
  f_val_out_ptr = f_val_out
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_real16
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
    f_val_inout_opt_ptr = f_val_inout_opt_native
  else
    ! f_val_inout_opt unset
  endif
end subroutine
subroutine fortran_test_real16_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(real16_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(real16_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(real16_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(real16_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(real16_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container general array (1D_ALLOC_real16)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container general array (1D_ALLOC_real16)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container general array (1D_ALLOC_real16)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container general array (1D_ALLOC_real16)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container general array (1D_ALLOC_real16)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_real16_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_complex_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  complex(c_double_complex) :: val_in  ! 0D_NOT_complex
  complex(rp) :: f_val_in
  type(c_ptr), intent(in), value :: val_in_opt  ! 0D_NOT_complex
  complex(c_double_complex) :: f_val_in_opt
  complex(c_double_complex), pointer :: f_val_in_opt_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out  ! 0D_NOT_complex
  complex(rp) :: f_val_out
  complex(c_double_complex), pointer :: f_val_out_ptr
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout  ! 0D_NOT_complex
  complex(c_double_complex) :: f_val_inout
  complex(c_double_complex), pointer :: f_val_inout_ptr
  type(c_ptr), intent(in), value :: val_inout_opt  ! 0D_NOT_complex
  complex(c_double_complex) :: f_val_inout_opt
  complex(c_double_complex), pointer :: f_val_inout_opt_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_complex
  f_val_in = val_in
  ! inout: f_val_inout 0D_NOT_complex
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
  else
    f_val_inout_ptr => null()
  endif
  ! in: f_val_in_opt 0D_NOT_complex
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr)
  else
    f_val_in_opt_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_complex
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
  else
    f_val_inout_opt_ptr => null()
  endif
  call test_complex_scalar(val_in=f_val_in, val_inout=f_val_inout_ptr, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt_ptr, val_inout_opt=f_val_inout_opt_ptr)

  ! inout: f_val_inout 0D_NOT_complex
  ! no output conversion for f_val_inout
  ! out: f_val_out 0D_NOT_complex
  call c_f_pointer(val_out, f_val_out_ptr)
  f_val_out_ptr = f_val_out
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_complex
  ! no output conversion for f_val_inout_opt
end subroutine
subroutine fortran_test_complex_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(complex_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(complex_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(complex_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(complex_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(complex_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container general array (1D_ALLOC_complex)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container general array (1D_ALLOC_complex)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container general array (1D_ALLOC_complex)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container general array (1D_ALLOC_complex)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container general array (1D_ALLOC_complex)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_complex_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_logical_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  logical(c_bool) :: val_in  ! 0D_NOT_logical
  logical :: f_val_in
  type(c_ptr), intent(in), value :: val_in_opt  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_val_in_opt
  logical, target :: f_val_in_opt_native
  logical, pointer :: f_val_in_opt_native_ptr
  logical(c_bool), pointer :: f_val_in_opt_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out  ! 0D_NOT_logical
  logical :: f_val_out
  logical(c_bool), pointer :: f_val_out_ptr
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_val_inout
  logical, target :: f_val_inout_native
  logical, pointer :: f_val_inout_native_ptr
  logical(c_bool), pointer :: f_val_inout_ptr
  type(c_ptr), intent(in), value :: val_inout_opt  ! 0D_NOT_logical
  logical(c_bool), pointer :: f_val_inout_opt
  logical, target :: f_val_inout_opt_native
  logical, pointer :: f_val_inout_opt_native_ptr
  logical(c_bool), pointer :: f_val_inout_opt_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_logical
  f_val_in = val_in
  ! inout: f_val_inout 0D_NOT_logical
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
    f_val_inout_native = f_val_inout_ptr
    f_val_inout_native_ptr => f_val_inout_native
  else
    f_val_inout_native_ptr => null()
  endif
  ! in: f_val_in_opt 0D_NOT_logical
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr)
    f_val_in_opt_native = f_val_in_opt_ptr
    f_val_in_opt_native_ptr => f_val_in_opt_native
  else
    f_val_in_opt_native_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_logical
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
    f_val_inout_opt_native = f_val_inout_opt_ptr
    f_val_inout_opt_native_ptr => f_val_inout_opt_native
  else
    f_val_inout_opt_native_ptr => null()
  endif
  call test_logical_scalar(val_in=f_val_in, val_inout=f_val_inout_native_ptr, &
      val_out=f_val_out, opt_status=f_opt_status, val_in_opt=f_val_in_opt_native_ptr, &
      val_inout_opt=f_val_inout_opt_native_ptr)

  ! inout: f_val_inout 0D_NOT_logical
  if (c_associated(val_inout)) then
    call c_f_pointer(val_inout, f_val_inout_ptr)
    f_val_inout_ptr = f_val_inout_native
  else
    ! f_val_inout unset
  endif
  ! out: f_val_out 0D_NOT_logical
  call c_f_pointer(val_out, f_val_out_ptr)
  f_val_out_ptr = f_val_out
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_logical
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr)
    f_val_inout_opt_ptr = f_val_inout_opt_native
  else
    ! f_val_inout_opt unset
  endif
end subroutine
subroutine fortran_test_logical_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(logical_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(logical_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(logical_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(logical_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(logical_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container general array (1D_ALLOC_logical)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container general array (1D_ALLOC_logical)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container general array (1D_ALLOC_logical)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container general array (1D_ALLOC_logical)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container general array (1D_ALLOC_logical)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_logical_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_character_scalar (val_in, val_inout, val_out, opt_status, val_in_opt, &
    val_inout_opt) bind(c)

  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: val_in
  character(len=4096), target :: f_val_in
  character(kind=c_char), pointer :: f_val_in_ptr(:)
  type(c_ptr), intent(in), value :: val_in_opt
  character(len=4096), target :: f_val_in_opt
  character(kind=c_char), pointer :: f_val_in_opt_ptr(:)
  character(len=4096), pointer :: f_val_in_opt_call_ptr
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: val_out
  character(len=4096), target :: f_val_out
  character(kind=c_char), pointer :: f_val_out_ptr(:)
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: val_inout
  character(len=4096), target :: f_val_inout
  character(kind=c_char), pointer :: f_val_inout_ptr(:)
  type(c_ptr), intent(in), value :: val_inout_opt
  character(len=4096), target :: f_val_inout_opt
  character(kind=c_char), pointer :: f_val_inout_opt_ptr(:)
  character(len=4096), pointer :: f_val_inout_opt_call_ptr
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_character
  if (.not. c_associated(val_in)) return
  call c_f_pointer(val_in, f_val_in_ptr, [huge(0)])
  call to_f_str(f_val_in_ptr, f_val_in)
  ! inout: f_val_inout 0D_NOT_character
  if (.not. c_associated(val_inout)) return
  call c_f_pointer(val_inout, f_val_inout_ptr, [huge(0)])
  call to_f_str(f_val_inout_ptr, f_val_inout)
  ! in: f_val_in_opt 0D_NOT_character
  if (c_associated(val_in_opt)) then
    call c_f_pointer(val_in_opt, f_val_in_opt_ptr, [huge(0)])
    call to_f_str(f_val_in_opt_ptr, f_val_in_opt)
    f_val_in_opt_call_ptr => f_val_in_opt
  else
    f_val_in_opt_call_ptr => null()
  endif
  ! inout: f_val_inout_opt 0D_NOT_character
  if (c_associated(val_inout_opt)) then
    call c_f_pointer(val_inout_opt, f_val_inout_opt_ptr, [huge(0)])
    call to_f_str(f_val_inout_opt_ptr, f_val_inout_opt)
    f_val_inout_opt_call_ptr => f_val_inout_opt
  else
    f_val_inout_opt_call_ptr => null()
  endif
  call test_character_scalar(val_in=f_val_in, val_inout=f_val_inout, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt_call_ptr, &
      val_inout_opt=f_val_inout_opt_call_ptr)

  ! inout: f_val_inout 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
  ! out: f_val_out 0D_NOT_character
  call c_f_pointer(val_out, f_val_out_ptr, [len_trim(f_val_out) + 1]) ! output-only string
  call to_c_str(f_val_out, f_val_out_ptr)
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
  ! inout: f_val_inout_opt 0D_NOT_character
  ! TODO i/o string (max length issue; buffer overflow...)
end subroutine
subroutine fortran_test_bunch_struct_scalar (val_in, val_inout, val_out, opt_status, &
    val_in_opt, val_inout_opt) bind(c)

  use bmad_struct, only: bunch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), value :: val_in  ! 0D_NOT_type
  type(bunch_struct), pointer :: f_val_in
  type(c_ptr), value :: val_in_opt  ! 0D_NOT_type
  type(bunch_struct), pointer :: f_val_in_opt
  ! ** Out parameters **
  type(c_ptr), value :: val_out  ! 0D_NOT_type
  type(bunch_struct), pointer :: f_val_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), value :: val_inout  ! 0D_NOT_type
  type(bunch_struct), pointer :: f_val_inout
  type(c_ptr), value :: val_inout_opt  ! 0D_NOT_type
  type(bunch_struct), pointer :: f_val_inout_opt
  ! ** End of parameters **
  ! in: f_val_in 0D_NOT_type
  if (.not. c_associated(val_in)) return
  call c_f_pointer(val_in, f_val_in)
  ! inout: f_val_inout 0D_NOT_type
  if (.not. c_associated(val_inout)) return
  call c_f_pointer(val_inout, f_val_inout)
  ! out: f_val_out 0D_NOT_type
  if (.not. c_associated(val_out)) return
  call c_f_pointer(val_out, f_val_out)
  ! in: f_val_in_opt 0D_NOT_type
  if (c_associated(val_in_opt))   call c_f_pointer(val_in_opt, f_val_in_opt)
  ! inout: f_val_inout_opt 0D_NOT_type
  if (c_associated(val_inout_opt))   call c_f_pointer(val_inout_opt, f_val_inout_opt)
  call test_bunch_struct_scalar(val_in=f_val_in, val_inout=f_val_inout, val_out=f_val_out, &
      opt_status=f_opt_status, val_in_opt=f_val_in_opt, val_inout_opt=f_val_inout_opt)

  ! out: f_val_out 0D_NOT_type
  ! TODO may require output conversion? 0D_NOT_type
  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine
subroutine fortran_test_bunch_struct_array (arr_in, arr_inout, arr_out, opt_status, arr_in_opt, &
    arr_inout_opt) bind(c)

  use bmad_struct, only: bunch_struct
  implicit none
  ! ** In parameters **
  type(c_ptr), intent(in), value :: arr_in
  type(bunch_struct_container_alloc), pointer :: f_arr_in
  type(c_ptr), intent(in), value :: arr_in_opt
  type(bunch_struct_container_alloc), pointer :: f_arr_in_opt
  ! ** Out parameters **
  type(c_ptr), intent(in), value :: arr_out
  type(bunch_struct_container_alloc), pointer :: f_arr_out
  type(c_ptr), intent(in), value :: opt_status
  integer :: f_opt_status(2)
  integer(c_int), pointer :: f_opt_status_ptr(:)
  ! ** Inout parameters **
  type(c_ptr), intent(in), value :: arr_inout
  type(bunch_struct_container_alloc), pointer :: f_arr_inout
  type(c_ptr), intent(in), value :: arr_inout_opt
  type(bunch_struct_container_alloc), pointer :: f_arr_inout_opt
  ! ** End of parameters **
  !! container type array (1D_ALLOC_type)
  if (c_associated(arr_in))   call c_f_pointer(arr_in, f_arr_in)
  !! container type array (1D_ALLOC_type)
  if (c_associated(arr_inout))   call c_f_pointer(arr_inout, f_arr_inout)
  !! container type array (1D_ALLOC_type)
  if (c_associated(arr_out))   call c_f_pointer(arr_out, f_arr_out)
  !! container type array (1D_ALLOC_type)
  if (c_associated(arr_in_opt))   call c_f_pointer(arr_in_opt, f_arr_in_opt)
  !! container type array (1D_ALLOC_type)
  if (c_associated(arr_inout_opt))   call c_f_pointer(arr_inout_opt, f_arr_inout_opt)
  call test_bunch_struct_array(arr_in=f_arr_in%data, arr_inout=f_arr_inout%data, &
      arr_out=f_arr_out%data, opt_status=f_opt_status, arr_in_opt=f_arr_in_opt%data, &
      arr_inout_opt=f_arr_inout_opt%data)

  ! out: f_opt_status 1D_NOT_integer
  if (c_associated(opt_status)) then
    call c_f_pointer(opt_status, f_opt_status_ptr, [2])
    f_opt_status_ptr = f_opt_status(:)
  endif
end subroutine

end module cppbmad_cppbmad_test_routines
