module ${module_name}

use bmad_interface
use bmad_struct
use fortran_cpp_utils
use precision_def ! , only: global_com, rp

use bmad_struct_proxy_mod

! ${use_lines}

use, intrinsic :: iso_c_binding

contains

! shorthand for c_associated since we're going to use it a lot here
elemental function assc(ptr) result(associated)
  type(c_ptr), intent(in) :: ptr
  logical :: associated
  
  associated = c_associated(ptr)
end function assc

! ${routines}

end module ${module_name}
