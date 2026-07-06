!------------------------------------------------------------------------------
! hook_helpers
!
! Shared helpers for the hand-written Bmad/Tao hook trampolines
! (src/hooks/tao_hooks.f90, src/hooks/bmad_hooks.f90).
!
! The `ci_*` bind(c) abstract interfaces describing the C ABI of each hook's C
! callback live in the generated module `hook_c_interfaces`
! (src/generated/hook_interfaces.f90, from codegen/hooks.py); the trampolines
! `use` that module directly.
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
