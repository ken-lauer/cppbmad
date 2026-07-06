module tao_c_proxy_interface
  use bmad
  use tao_interface
  use fortran_cpp_utils, only: to_c_str

  implicit none

  ! C-compatible lattice type enumeration
  integer, parameter :: LATTICE_MODEL = 1
  integer, parameter :: LATTICE_DESIGN = 2  
  integer, parameter :: LATTICE_BASE = 3

contains

  function is_initialized() result(n) bind(c, name='tao_is_initialized')
    integer(c_int) :: n
    if (s%initialized) then
      n = 1
    else
      n = 0
    end if
  end function

  function c_get_space_charge_com() result(ptr) bind(c, name="bmad_get_space_charge_com")
    type(c_ptr) :: ptr
    ptr = c_loc(space_charge_com)
  end function

  function c_get_bmad_com() result(ptr) bind(c, name="bmad_get_bmad_com")
    type(c_ptr) :: ptr
    ptr = c_loc(bmad_com)
  end function

  function c_get_super_universe() result(ptr) bind(c, name="tao_get_super_universe_ptr")
    type(c_ptr) :: ptr
    ptr = c_loc(s)
  end function

end module tao_c_proxy_interface
