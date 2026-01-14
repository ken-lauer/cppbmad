module array_desc_mod
    use iso_c_binding
    implicit none

    integer, parameter :: MAX_ARRAY_RANK = 3

    type, bind(c) :: array_descriptor_t
        type(c_ptr) :: data_ptr
        integer(c_int) :: rank
        integer(c_int) :: dims(MAX_ARRAY_RANK)
        integer(c_int) :: strides(MAX_ARRAY_RANK)
    end type array_descriptor_t

! contains

end module array_desc_mod
