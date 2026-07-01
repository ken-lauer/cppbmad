!------------------------------------------------------------------------------
! hook_helpers
!
! Shared helpers for the hand-written Bmad/Tao hook trampolines
! (src/hooks/tao_hooks.f90, src/hooks/bmad_hooks.f90):
!
!   * `ci_*` : bind(c) abstract interfaces describing the C ABI of each hook's C
!     callback. Derived types cross as `type(c_ptr), value`; logical/integer
!     out-parameters cross by reference as `integer(c_int)`; `real(rp)` scalars
!     cross as `real(c_double)`; assumed-shape derived-type arrays cross as a base
!     pointer plus bounds and element size.
!
!   * `cloc_*` : take a derived-type argument into a `target` dummy and return its
!     address as a `type(c_ptr)`. Hook interfaces frequently do NOT declare their
!     dummies `target`, and gfortran enforces the TARGET attribute when
!     associating a trampoline with a hook procedure pointer. `c_loc` requires a
!     target, so trampolines route each derived-type argument through one of these
!     (passing a non-target actual to a target dummy is legal; the address is the
!     real storage address and remains valid after the call because scalar
!     derived-type arguments are passed by reference without a copy).
!
! Maintained by hand. See plans/hook-plan.md.
!------------------------------------------------------------------------------

module hook_helpers

  use, intrinsic :: iso_c_binding
  use bmad
  use tao_interface

  implicit none

  abstract interface

    ! logical out
    subroutine ci_i(iout) bind(c)
      import :: c_int
      integer(c_int) :: iout
    end subroutine

    ! int in, int in, derived
    subroutine ci_iip(i1, i2, p) bind(c)
      import :: c_int, c_ptr
      integer(c_int), value :: i1, i2
      type(c_ptr), value :: p
    end subroutine

    ! int in, int in, derived, logical out
    subroutine ci_iipi(i1, i2, p, iout) bind(c)
      import :: c_int, c_ptr
      integer(c_int), value :: i1, i2
      type(c_ptr), value :: p
      integer(c_int) :: iout
    end subroutine

    ! derived, derived, real in
    subroutine ci_pp_d(a, b, d) bind(c)
      import :: c_ptr, c_double
      type(c_ptr), value :: a, b
      real(c_double), value :: d
    end subroutine

    ! derived, derived, logical out
    subroutine ci_pp_i(a, b, iout) bind(c)
      import :: c_ptr, c_int
      type(c_ptr), value :: a, b
      integer(c_int) :: iout
    end subroutine

    ! four derived
    subroutine ci_p4(a, b, c, d) bind(c)
      import :: c_ptr
      type(c_ptr), value :: a, b, c, d
    end subroutine

    ! derived x3, real inout, int inout
    subroutine ci_ppp_di(a, b, c, d, iout) bind(c)
      import :: c_ptr, c_double, c_int
      type(c_ptr), value :: a, b, c
      real(c_double) :: d
      integer(c_int) :: iout
    end subroutine

    ! start_orb, ele, param, err_flag(out), finished(out), track(optional derived)
    subroutine ci_track1_custom(a, b, c, i1, i2, d) bind(c)
      import :: c_ptr, c_int
      type(c_ptr), value :: a, b, c
      integer(c_int) :: i1, i2
      type(c_ptr), value :: d
    end subroutine

    ! start_orb, ele, param, err_flag(out), finished(out), radiation_included(out), track(optional)
    subroutine ci_track1_preprocess(a, b, c, i1, i2, i3, d) bind(c)
      import :: c_ptr, c_int
      type(c_ptr), value :: a, b, c
      integer(c_int) :: i1, i2, i3
      type(c_ptr), value :: d
    end subroutine

    ! start_orb, ele, param, end_orb, err_flag(out), make_quaternion(optional logical, byref via ptr)
    subroutine ci_track1_spin(a, b, c, d, i1, e) bind(c)
      import :: c_ptr, c_int
      type(c_ptr), value :: a, b, c, d
      integer(c_int) :: i1
      type(c_ptr), value :: e
    end subroutine

    ! bunch, ele, err(out), centroid(optional array), direction(optional int ptr),
    ! finished(out), bunch_track(optional derived)
    subroutine ci_track1_bunch(bunch, ele, err, cdata, clb, cub, cesize, dir, finished, btrack) bind(c)
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), value :: bunch, ele
      integer(c_int) :: err
      type(c_ptr), value :: cdata
      integer(c_int), value :: clb, cub
      integer(c_size_t), value :: cesize
      type(c_ptr), value :: dir
      integer(c_int) :: finished
      type(c_ptr), value :: btrack
    end subroutine

    ! finished(out), lat, orbit(array), ix_start/end/direction(int in), ix_branch/track_state(optional int ptr)
    subroutine ci_track_many(finished, lat, odata, olb, oub, oesize, ix_start, ix_end, dir, ix_branch, track_state) &
        bind(c)
      import :: c_ptr, c_int, c_size_t
      integer(c_int) :: finished
      type(c_ptr), value :: lat
      type(c_ptr), value :: odata
      integer(c_int), value :: olb, oub
      integer(c_size_t), value :: oesize
      integer(c_int), value :: ix_start, ix_end, dir
      type(c_ptr), value :: ix_branch, track_state
    end subroutine

  end interface

contains

  function cloc_coord(x) result(p)
    type(coord_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_ele(x) result(p)
    type(ele_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_lat_param(x) result(p)
    type(lat_param_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_lat(x) result(p)
    type(lat_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_bunch(x) result(p)
    type(bunch_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_track(x) result(p)
    type(track_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_bunch_track(x) result(p)
    type(bunch_track_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_tao_var(x) result(p)
    type(tao_var_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  function cloc_tao_data(x) result(p)
    type(tao_data_struct), target :: x
    type(c_ptr) :: p
    p = c_loc(x)
  end function

  ! Byte size of one coord_struct element, for array element striding. Array
  ! trampolines take the base address via cloc_coord on the first element
  ! (contiguous actual assumed; the hook callers pass whole arrays).
  function coord_elem_size() result(n)
    integer(c_size_t) :: n
    type(coord_struct) :: tmp
    n = storage_size(tmp) / 8
  end function

end module hook_helpers
