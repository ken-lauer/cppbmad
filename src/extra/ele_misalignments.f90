!+
! Module ele_misalignments_mod
!
! Fast batch setter for misalignment-class element attributes
! (x_offset, y_offset, z_offset, x_pitch, y_pitch, tilt).
!
! Usage:
!   if (.not. set_ele_misalignments(ele, dx, dy, dz, xp, yp, tlt)) cycle  ! skip
!   ...
!   call lattice_bookkeeper(lat)                                            ! once at end
!
! With check_free = .true. (default) the function rejects super_slaves,
! slice_slaves, and multipass_slaves (whose value() is overwritten by the
! next lattice_bookkeeper) — caller must route to the appropriate lord via
! pointer_to_lord() and retry on that. Girder slaves (minor_slave$) and free
! elements are accepted: a girder's contribution rides on *_tot and the
! slave's own x_offset / etc. survive the bookkeeper pass.
!-

module cppbmad_ele_misalignments_mod

use bmad

implicit none

private :: attr_ok, can_set_ele_misalignments
public :: set_ele_misalignments

contains

!+
! Function set_ele_misalignments (ele, x_offset, y_offset, z_offset, x_pitch, y_pitch, tilt, check_free) result (ok)
!
! Set the six misalignment attributes on ele in one call and mark bookkeeping
! flags once at the end. Equivalent to six set_ele_attribute calls minus the
! string parsing and minus the per-call attribute_set_bookkeeping work.
! Caller must invoke lattice_bookkeeper before tracking.
!
! The trailing underscore in the name marks this as a pybmad-side helper, not
! a routine from Bmad proper.
!
! When check_free is true (default), all six attributes are verified to
! resolve via attribute_index() to their expected slots AND to satisfy
! attribute_free(..., dependent_attribs_free = .true.). The check is
! all-or-nothing: if any attribute fails, the function returns .false. and ele
! is not modified. With check_free = .false. the writes proceed unconditionally;
! use only when the ele key is known to support all six attributes.
!
! Input:
!   ele         -- ele_struct: Element to modify.
!   x_offset    -- real(rp): the x-offset.
!   y_offset    -- real(rp): the y-offset.
!   z_offset    -- real(rp): the z-offset.
!   x_pitch     -- real(rp): the horizontal pitch.
!   y_pitch     -- real(rp): the vertical pitch.
!   tilt        -- real(rp): the tilt (or roll, for bends).
!   check_free  -- logical, optional: Default .true.. If true, validate slot
!                  layout and attribute freeness before any write.
!
! Output:
!   ele         -- ele_struct: Values stored and bookkeeping flags marked stale
!                  (attribute_group, mat6_group, control_group as appropriate);
!                  cached Taylor maps killed. Unchanged if ok == .false..
!   ok          -- logical: .true. on success, .false. if check_free caught a
!                  non-standard or non-free attribute.
!-

function set_ele_misalignments (ele, x_offset, y_offset, z_offset, x_pitch, y_pitch, tilt, check_free) result (ok)

type (ele_struct), target :: ele
real(rp), intent(in) :: x_offset, y_offset, z_offset
real(rp), intent(in) :: x_pitch, y_pitch, tilt
logical, intent(in), optional :: check_free
logical :: ok

!

if (logic_option(.true., check_free)) then
  if (.not. can_set_ele_misalignments(ele)) then
    ok = .false.
    return
  endif
endif

ele%value(x_offset$) = x_offset
ele%value(y_offset$) = y_offset
ele%value(z_offset$) = z_offset
ele%value(x_pitch$)  = x_pitch
ele%value(y_pitch$)  = y_pitch
ele%value(tilt$)     = tilt  ! tilt$ = roll$, importantly, per bmad_struct

! lattice_bookkeeper will run attribute_bookkeeper because attribute_group is
! stale, and that pass recomputes x_offset_tot etc.).

if (ele%key /= overlay$ .and. ele%key /= group$) then
  call set_ele_status_stale(ele, attribute_group$)
endif

if (ele%key /= overlay$ .and. ele%key /= group$ .and. &
    ele%lord_status /= multipass_lord$) then
  call set_ele_status_stale(ele, mat6_group$)
endif

if (ele%lord_status /= not_a_lord$) then
  call set_ele_status_stale(ele, control_group$)
endif

! Cached Taylor maps are invalidated by any geometry change.

if (ele%key /= taylor$ .and. associated(ele%taylor(1)%term)) then
  call kill_taylor(ele%taylor)
  call kill_taylor(ele%spin_taylor)
endif

ok = .true.

end function set_ele_misalignments

!--------------------------------------------------------------
! Private: check if the element misalignments can be set.
!
! Reject slave statuses whose value() is overwritten by lattice_bookkeeper —
! super_slave (makeup_super_slave1 re-derives misalignments from lord *_tot),
! multipass_slave (makeup_multipass_slave does wholesale slave%value =
! lord%value), and slice_slave (transfer_ele + makeup_super_slave1 from the
! parent). Caller must route to the appropriate lord. free$ and minor_slave$
! (girder slaves — girder's contribution rides on *_tot, the slave's own
! x_offset etc. survive) are OK to write directly.

function can_set_ele_misalignments (ele) result (ok)

type (ele_struct) :: ele
logical :: ok

!

select case (ele%slave_status)
case (super_slave$, multipass_slave$, slice_slave$)
  ok = .false.
  return
end select

ok = attr_ok(ele, 'X_OFFSET', x_offset$) .and. &
     attr_ok(ele, 'Y_OFFSET', y_offset$) .and. &
     attr_ok(ele, 'Z_OFFSET', z_offset$) .and. &
     attr_ok(ele, 'X_PITCH',  x_pitch$)  .and. &
     attr_ok(ele, 'Y_PITCH',  y_pitch$)

if (.not. ok) return

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  ok = attr_ok(ele, 'ROLL', roll$)
else
  ok = attr_ok(ele, 'TILT', tilt$)
endif

end function can_set_ele_misalignments

!--------------------------------------------------------------
! Private: name resolves to the expected slot AND is settable.

function attr_ok (ele, name, expected_ix) result (ok)

type (ele_struct) :: ele
character(*), intent(in) :: name
integer, intent(in) :: expected_ix
logical :: ok

integer :: ix

!

ix = attribute_index(ele, name)
if (ix /= expected_ix) then
  ok = .false.
  return
endif
ok = attribute_free(ele, name, .false., dependent_attribs_free = .true.)

end function attr_ok

end module cppbmad_ele_misalignments_mod
