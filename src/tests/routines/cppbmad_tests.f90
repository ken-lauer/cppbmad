module cppbmad_test_mod

  use iso_c_binding
  use bmad_struct, only: bunch_struct
  use precision_def, only: rp, qp
  implicit none

  integer, parameter :: i8 = 8 

contains

  ! ======================================================================
  ! INTEGER TESTS
  ! ======================================================================

  subroutine test_integer_scalar(val_in, val_inout, val_out, &
                                 opt_status, val_in_opt, val_inout_opt)
    integer, intent(in) :: val_in
    integer, intent(inout) :: val_inout
    integer, intent(out) :: val_out
    integer, intent(out) :: opt_status(2) 
    integer, optional, intent(in) :: val_in_opt
    integer, optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = val_inout + 1

    opt_status = 0
    if (present(val_in_opt)) then 
       opt_status(1) = 1
       val_out = val_out + val_in_opt
    end if

    if (present(val_inout_opt)) then
       opt_status(2) = 1
       val_inout_opt = val_inout_opt + 1
    end if

  end subroutine test_integer_scalar

  subroutine test_integer_array(arr_in, arr_inout, arr_out, &
                                opt_status, arr_in_opt, arr_inout_opt)
    integer, dimension(:), intent(in) :: arr_in
    integer, dimension(:), intent(inout) :: arr_inout
    integer, dimension(:), allocatable, intent(out) :: arr_out
    integer, intent(out) :: opt_status(2)
    integer, dimension(:), optional, intent(in) :: arr_in_opt
    integer, dimension(:), optional, intent(inout) :: arr_inout_opt

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    ! Ensure we have capacity before writing
    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if

    if (size(arr_inout) > 0) then
        arr_inout = arr_inout + 1
    end if

    opt_status = 0
    if (present(arr_in_opt)) then 
       opt_status(1) = 1
       if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
           arr_out = arr_out + arr_in_opt
       end if
    end if

    if (present(arr_inout_opt)) then
       opt_status(2) = 1
       if (size(arr_inout_opt) > 0) then
           arr_inout_opt = arr_inout_opt + 1
       end if
    end if
  end subroutine test_integer_array


  ! ======================================================================
  ! INTEGER8 TESTS
  ! ======================================================================

  subroutine test_integer8_scalar(val_in, val_inout, val_out, &
                                  opt_status, val_in_opt, val_inout_opt)
    integer(8), intent(in) :: val_in
    integer(8), intent(inout) :: val_inout
    integer(8), intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    integer(8), optional, intent(in) :: val_in_opt
    integer(8), optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = val_inout + 1_i8

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        val_out = val_out + val_in_opt
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt = val_inout_opt + 1_i8
    end if
  end subroutine test_integer8_scalar

  subroutine test_integer8_array(arr_in, arr_inout, arr_out, &
                                 opt_status, arr_in_opt, arr_inout_opt)
    integer(8), dimension(:), intent(in) :: arr_in
    integer(8), dimension(:), intent(inout) :: arr_inout
    integer(8), dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    integer(8), dimension(:), optional, intent(in) :: arr_in_opt
    integer(8), dimension(:), optional, intent(inout) :: arr_inout_opt

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if
    
    if (size(arr_inout) > 0) then
        arr_inout = arr_inout + 1_i8
    end if

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
        if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
            arr_out = arr_out + arr_in_opt
        end if
    end if
    if (present(arr_inout_opt)) then
        opt_status(2) = 1
         if (size(arr_inout_opt) > 0) then
            arr_inout_opt = arr_inout_opt + 1_i8
        end if
    end if
  end subroutine test_integer8_array


  ! ======================================================================
  ! REAL (Double) TESTS
  ! ======================================================================

  subroutine test_real_scalar(val_in, val_inout, val_out, &
                              opt_status, val_in_opt, val_inout_opt)
    real(rp), intent(in) :: val_in
    real(rp), intent(inout) :: val_inout
    real(rp), intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    real(rp), optional, intent(in) :: val_in_opt
    real(rp), optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = val_inout + 1.0_rp

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        val_out = val_out + val_in_opt
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt = val_inout_opt + 1.0_rp
    end if
  end subroutine test_real_scalar

  subroutine test_real_array(arr_in, arr_inout, arr_out, &
                                 opt_status, arr_in_opt, arr_inout_opt)
    real(rp), dimension(:), intent(in) :: arr_in
    real(rp), dimension(:), intent(inout) :: arr_inout
    real(rp), dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    real(rp), dimension(:), optional, intent(in) :: arr_in_opt
    real(rp), dimension(:), optional, intent(inout) :: arr_inout_opt

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if
    if (size(arr_inout) > 0) then
        arr_inout = arr_inout + 1.0_rp
    end if

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
        if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
            arr_out = arr_out + arr_in_opt
        end if
    end if
    if (present(arr_inout_opt)) then
        opt_status(2) = 1
        if (size(arr_inout_opt) > 0) then
            arr_inout_opt = arr_inout_opt + 1.0_rp
        end if
    end if
  end subroutine test_real_array


  ! ======================================================================
  ! REAL16 (Quad) TESTS
  ! ======================================================================

  subroutine test_real16_scalar(val_in, val_inout, val_out, &
                                opt_status, val_in_opt, val_inout_opt)
    real(qp), intent(in) :: val_in
    real(qp), intent(inout) :: val_inout
    real(qp), intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    real(qp), optional, intent(in) :: val_in_opt
    real(qp), optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = val_inout + 1.0_qp

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        val_out = val_out + val_in_opt
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt = val_inout_opt + 1.0_qp
    end if
  end subroutine test_real16_scalar

  subroutine test_real16_array(arr_in, arr_inout, arr_out, &
                               opt_status, arr_in_opt, arr_inout_opt)
    real(qp), dimension(:), intent(in) :: arr_in
    real(qp), dimension(:), intent(inout) :: arr_inout
    real(qp), dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    real(qp), dimension(:), optional, intent(in) :: arr_in_opt
    real(qp), dimension(:), optional, intent(inout) :: arr_inout_opt

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if
    if (size(arr_inout) > 0) then
        arr_inout = arr_inout + 1.0_qp
    end if

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
         if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
            arr_out = arr_out + arr_in_opt
        end if
    end if
    if (present(arr_inout_opt)) then
        opt_status(2) = 1
        if (size(arr_inout_opt) > 0) then
            arr_inout_opt = arr_inout_opt + 1.0_qp
        end if
    end if
  end subroutine test_real16_array


  ! ======================================================================
  ! COMPLEX TESTS
  ! ======================================================================

  subroutine test_complex_scalar(val_in, val_inout, val_out, &
                                 opt_status, val_in_opt, val_inout_opt)
    complex(rp), intent(in) :: val_in
    complex(rp), intent(inout) :: val_inout
    complex(rp), intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    complex(rp), optional, intent(in) :: val_in_opt
    complex(rp), optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = val_inout + cmplx(1.0_rp, 1.0_rp, rp)

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        val_out = val_out + val_in_opt
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt = val_inout_opt + cmplx(1.0_rp, 1.0_rp, rp)
    end if
  end subroutine test_complex_scalar

  subroutine test_complex_array(arr_in, arr_inout, arr_out, &
                                opt_status, arr_in_opt, arr_inout_opt)
    complex(rp), dimension(:), intent(in) :: arr_in
    complex(rp), dimension(:), intent(inout) :: arr_inout
    complex(rp), dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    complex(rp), dimension(:), optional, intent(in) :: arr_in_opt
    complex(rp), dimension(:), optional, intent(inout) :: arr_inout_opt

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if
    if (size(arr_inout) > 0) then
        arr_inout = arr_inout + cmplx(1.0_rp, 1.0_rp, rp)
    end if

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
        if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
            arr_out = arr_out + arr_in_opt
        end if
    end if
    if (present(arr_inout_opt)) then
        opt_status(2) = 1
         if (size(arr_inout_opt) > 0) then
            arr_inout_opt = arr_inout_opt + cmplx(1.0_rp, 1.0_rp, rp)
        end if
    end if
  end subroutine test_complex_array


  ! ======================================================================
  ! LOGICAL TESTS
  ! ======================================================================

  subroutine test_logical_scalar(val_in, val_inout, val_out, &
                                 opt_status, val_in_opt, val_inout_opt)
    logical, intent(in) :: val_in
    logical, intent(inout) :: val_inout
    logical, intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    logical, optional, intent(in) :: val_in_opt
    logical, optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = .not. val_inout

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        val_out = val_out .eqv. val_in_opt
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt = .not. val_inout_opt
    end if
  end subroutine test_logical_scalar

  subroutine test_logical_array(arr_in, arr_inout, arr_out, &
                                opt_status, arr_in_opt, arr_inout_opt)
    logical, dimension(:), intent(in) :: arr_in
    logical, dimension(:), intent(inout) :: arr_inout
    logical, dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    logical, dimension(:), optional, intent(in) :: arr_in_opt
    logical, dimension(:), optional, intent(inout) :: arr_inout_opt

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if
    if (size(arr_inout) > 0) then
        arr_inout = .not. arr_inout
    end if

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
        if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
            arr_out = arr_out .eqv. arr_in_opt
        end if
    end if
    if (present(arr_inout_opt)) then
        opt_status(2) = 1
        if (size(arr_inout_opt) > 0) then
            arr_inout_opt = .not. arr_inout_opt
        end if
    end if
  end subroutine test_logical_array


  ! ======================================================================
  ! CHARACTER TESTS
  ! ======================================================================

  subroutine test_character_scalar(val_in, val_inout, val_out, &
                                   opt_status, val_in_opt, val_inout_opt)
    character(*), intent(in) :: val_in
    character(*), intent(inout) :: val_inout
    character(*), intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    character(*), optional, intent(in) :: val_in_opt
    character(*), optional, intent(inout) :: val_inout_opt

    val_out = val_in
    val_inout = trim(val_inout) // "_mod"

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        val_out = trim(val_out) // trim(val_in_opt)
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt = trim(val_inout_opt) // "_mod"
    end if
  end subroutine test_character_scalar

  subroutine test_character_array(arr_in, arr_inout, arr_out, &
                                  opt_status, arr_in_opt, arr_inout_opt)
    character(*), dimension(:), intent(in) :: arr_in
    character(*), dimension(:), intent(inout) :: arr_inout
    character(*), dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    character(*), dimension(:), optional, intent(in) :: arr_in_opt
    character(*), dimension(:), optional, intent(inout) :: arr_inout_opt
    integer :: i

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        arr_out = arr_in
    end if

    do i = 1, size(arr_inout)
        arr_inout(i) = trim(arr_inout(i)) // "_mod"
    end do

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
         if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
            do i = 1, size(arr_out)
                arr_out(i) = trim(arr_out(i)) // trim(arr_in_opt(i))
            end do
        end if
    end if
    if (present(arr_inout_opt)) then
        opt_status(2) = 1
        do i = 1, size(arr_inout_opt)
            arr_inout_opt(i) = trim(arr_inout_opt(i)) // "_mod"
        end do
    end if
  end subroutine test_character_array


  ! ======================================================================
  ! STRUCT / TYPE TESTS (bunch_struct)
  ! ======================================================================

  subroutine test_bunch_struct_scalar(val_in, val_inout, val_out, &
                                       opt_status, val_in_opt, val_inout_opt)
    type(bunch_struct), intent(in) :: val_in
    type(bunch_struct), intent(inout) :: val_inout
    type(bunch_struct), intent(out) :: val_out
    integer, intent(out) :: opt_status(2)
    type(bunch_struct), optional, intent(in) :: val_in_opt
    type(bunch_struct), optional, intent(inout) :: val_inout_opt

    ! Intrinsic assignment for structs handles allocation of allocatable components 
    ! on the LHS (val_out) by copying the shape/data from val_in.
    val_out = val_in

    ! Modification: Increment ID, Add 1.0 to charge
    val_inout%ix_ele = val_inout%ix_ele + 1
    val_inout%charge_tot = val_inout%charge_tot + 1.0d0

    opt_status = 0
    if (present(val_in_opt)) then
        opt_status(1) = 1
        ! Modifying output based on opt input
        val_out%ix_ele = val_out%ix_ele + val_in_opt%ix_ele
        val_out%charge_tot = val_out%charge_tot + val_in_opt%charge_tot
    end if
    if (present(val_inout_opt)) then
        opt_status(2) = 1
        val_inout_opt%ix_ele = val_inout_opt%ix_ele + 1
        val_inout_opt%charge_tot = val_inout_opt%charge_tot + 1.0d0
    end if
  end subroutine test_bunch_struct_scalar

  subroutine test_bunch_struct_array(arr_in, arr_inout, arr_out, &
                                      opt_status, arr_in_opt, arr_inout_opt)
    type(bunch_struct), dimension(:), intent(in) :: arr_in
    type(bunch_struct), dimension(:), intent(inout) :: arr_inout
    type(bunch_struct), dimension(:), intent(out), allocatable :: arr_out
    integer, intent(out) :: opt_status(2)
    type(bunch_struct), dimension(:), optional, intent(in) :: arr_in_opt
    type(bunch_struct), dimension(:), optional, intent(inout) :: arr_inout_opt
    integer :: i

    if (allocated(arr_out)) deallocate(arr_out)
    allocate(arr_out(lbound(arr_in, 1):ubound(arr_in, 1)))

    ! Ensure calling side has allocated array shape correctly
    if (size(arr_out) == size(arr_in) .and. size(arr_in) > 0) then
        ! Fortran array assignment of derived types with allocatables 
        ! effectively does a deep copy element-by-element
        arr_out = arr_in
    end if
    
    do i = 1, size(arr_inout)
        arr_inout(i)%ix_ele = arr_inout(i)%ix_ele + 1
        arr_inout(i)%charge_tot = arr_inout(i)%charge_tot + 1.0d0
    end do

    opt_status = 0
    if (present(arr_in_opt)) then
        opt_status(1) = 1
         if (size(arr_in_opt) == size(arr_out) .and. size(arr_out) > 0) then
             do i = 1, size(arr_out)
                 arr_out(i)%ix_ele = arr_out(i)%ix_ele + arr_in_opt(i)%ix_ele
                 arr_out(i)%charge_tot = arr_out(i)%charge_tot + arr_in_opt(i)%charge_tot
             end do
         end if
    end if

    if (present(arr_inout_opt)) then
        opt_status(2) = 1
        do i = 1, size(arr_inout_opt)
            arr_inout_opt(i)%ix_ele = arr_inout_opt(i)%ix_ele + 1
            arr_inout_opt(i)%charge_tot = arr_inout_opt(i)%charge_tot + 1.0d0
        end do
    end if
  end subroutine test_bunch_struct_array

end module cppbmad_test_mod
