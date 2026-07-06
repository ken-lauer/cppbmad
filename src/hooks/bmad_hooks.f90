!------------------------------------------------------------------------------
! bmad_c_hook_interface
!
! Hand-written Fortran trampolines that let a C (and, through cppbmad, Python)
! callable be installed into Bmad's tracking hook / custom procedure pointers
! declared in bmad_routine_interface.
!
! For each supported hook `<x>_def` this module provides:
!   * a module `type(c_funptr) :: cb_<x>` holding the C callback,
!   * a trampoline `tramp_<x>` matching the abstract interface `<x>_def`,
!   * a `bind(c)` registration entry point `bmad_hook_register_<x>`.
!
! Conventions (see plans/hook-plan.md):
!   * derived types cross as opaque `type(c_ptr)` via the shared `cloc_*` helpers,
!   * logical/integer out-parameters cross by reference as `integer(c_int)`,
!   * `real(rp)` scalars cross as `real(c_double)`,
!   * optional derived / integer / logical args cross as a `type(c_ptr)` that is
!     null when the argument is absent,
!   * assumed-shape `coord_struct` arrays cross as base pointer + bounds + element
!     size (contiguous actual assumed).
!------------------------------------------------------------------------------

module bmad_c_hook_interface

  use, intrinsic :: iso_c_binding
  use bmad_routine_interface
  use hook_helpers
  use hook_c_interfaces

  implicit none

  type(c_funptr) :: cb_time_rk = c_null_funptr
  type(c_funptr) :: cb_track1_bunch = c_null_funptr
  type(c_funptr) :: cb_track1_custom = c_null_funptr
  type(c_funptr) :: cb_track_many = c_null_funptr
  type(c_funptr) :: cb_track1_postprocess = c_null_funptr
  type(c_funptr) :: cb_track1_preprocess = c_null_funptr
  type(c_funptr) :: cb_track1_spin_custom = c_null_funptr
  type(c_funptr) :: cb_track1_wake = c_null_funptr
  type(c_funptr) :: cb_wall_hit = c_null_funptr

contains

  !----------------------------------------------------------------------------
  ! time_runge_kutta_periodic_kick_hook (orbit, ele, param, stop_time, init_needed)
  !----------------------------------------------------------------------------

  subroutine tramp_time_rk(orbit, ele, param, stop_time, init_needed)
    type(coord_struct) :: orbit
    type(ele_struct) :: ele
    type(lat_param_struct) :: param
    real(rp) :: stop_time
    integer :: init_needed
    procedure(ci_time_runge_kutta_periodic_kick), pointer :: cproc
    real(c_double) :: cst
    integer(c_int) :: cin
    call c_f_procpointer(cb_time_rk, cproc)
    cst = stop_time
    cin = init_needed
    call cproc(cloc_coord(orbit), cloc_ele(ele), cloc_lat_param(param), cst, cin)
    stop_time = cst
    init_needed = cin
  end subroutine

  subroutine bmad_hook_register_time_rk(cfun) &
      bind(c, name='bmad_hook_register_time_runge_kutta_periodic_kick')
    type(c_funptr), value :: cfun
    cb_time_rk = cfun
    if (c_associated(cfun)) then
      time_runge_kutta_periodic_kick_hook_ptr => tramp_time_rk
    else
      time_runge_kutta_periodic_kick_hook_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)
  !----------------------------------------------------------------------------

  subroutine tramp_track1_bunch(bunch, ele, err, centroid, direction, finished, bunch_track)
    type(bunch_struct), target :: bunch
    type(ele_struct), target :: ele
    logical :: err, finished
    type(coord_struct), optional :: centroid(0:)
    integer, optional :: direction
    type(bunch_track_struct), optional :: bunch_track
    procedure(ci_track1_bunch), pointer :: cproc
    integer(c_int) :: cerr, cfin, clb, cub
    integer(c_size_t) :: cesize
    type(c_ptr) :: cdata, dirp, btp
    integer(c_int), target :: cdir
    call c_f_procpointer(cb_track1_bunch, cproc)
    ! err/finished are hook outputs (Bmad may pass them uninitialized). Default them
    ! to false so a callback returning None (leave unchanged) observes safely rather
    ! than inheriting an uninitialized value that aborts/replaces the tracking.
    cerr = 0
    cfin = 0
    cdata = c_null_ptr
    clb = 0
    cub = -1
    if (present(centroid)) then
      if (size(centroid) > 0) then
        cdata = cloc_coord(centroid(lbound(centroid, 1)))
        clb = lbound(centroid, 1)
        cub = ubound(centroid, 1)
      end if
    end if
    cesize = coord_elem_size()
    if (present(direction)) then
      cdir = direction
      dirp = c_loc(cdir)
    else
      dirp = c_null_ptr
    end if
    if (present(bunch_track)) then
      btp = cloc_bunch_track(bunch_track)
    else
      btp = c_null_ptr
    end if
    call cproc(cloc_bunch(bunch), cloc_ele(ele), cerr, cdata, clb, cub, cesize, dirp, cfin, btp)
    err = (cerr /= 0)
    finished = (cfin /= 0)
    if (present(direction)) direction = cdir
  end subroutine

  subroutine bmad_hook_register_track1_bunch(cfun) bind(c, name='bmad_hook_register_track1_bunch')
    type(c_funptr), value :: cfun
    cb_track1_bunch = cfun
    if (c_associated(cfun)) then
      track1_bunch_hook_ptr => tramp_track1_bunch
    else
      track1_bunch_hook_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track1_custom (start_orb, ele, param, err_flag, finished, track)
  !----------------------------------------------------------------------------

  subroutine tramp_track1_custom(start_orb, ele, param, err_flag, finished, track)
    type(coord_struct) :: start_orb
    type(ele_struct) :: ele
    type(lat_param_struct) :: param
    logical :: err_flag, finished
    type(track_struct), optional :: track
    procedure(ci_track1_custom), pointer :: cproc
    integer(c_int) :: ce, cf
    type(c_ptr) :: tp
    call c_f_procpointer(cb_track1_custom, cproc)
    ! err_flag/finished are hook outputs; default them to false (see tramp_track1_bunch).
    ce = 0
    cf = 0
    if (present(track)) then
      tp = cloc_track(track)
    else
      tp = c_null_ptr
    end if
    call cproc(cloc_coord(start_orb), cloc_ele(ele), cloc_lat_param(param), ce, cf, tp)
    err_flag = (ce /= 0)
    finished = (cf /= 0)
  end subroutine

  subroutine bmad_hook_register_track1_custom(cfun) bind(c, name='bmad_hook_register_track1_custom')
    type(c_funptr), value :: cfun
    cb_track1_custom = cfun
    if (c_associated(cfun)) then
      track1_custom_ptr => tramp_track1_custom
    else
      track1_custom_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track_many_hook (finished, lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
  !----------------------------------------------------------------------------

  subroutine tramp_track_many(finished, lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
    logical :: finished
    type(lat_struct), target :: lat
    type(coord_struct) :: orbit(0:)
    integer :: ix_start, ix_end, direction
    integer, optional :: ix_branch, track_state
    procedure(ci_track_many), pointer :: cproc
    integer(c_int) :: cfin, olb, oub
    integer(c_size_t) :: oesize
    type(c_ptr) :: odata, ixbp, tsp
    integer(c_int), target :: cixb, cts
    call c_f_procpointer(cb_track_many, cproc)
    ! finished is a hook output; default it to false (see tramp_track1_bunch).
    cfin = 0
    olb = lbound(orbit, 1)
    oub = ubound(orbit, 1)
    if (size(orbit) > 0) then
      odata = cloc_coord(orbit(lbound(orbit, 1)))
    else
      odata = c_null_ptr
    end if
    oesize = coord_elem_size()
    if (present(ix_branch)) then
      cixb = ix_branch
      ixbp = c_loc(cixb)
    else
      ixbp = c_null_ptr
    end if
    if (present(track_state)) then
      cts = track_state
      tsp = c_loc(cts)
    else
      tsp = c_null_ptr
    end if
    call cproc(cfin, cloc_lat(lat), odata, olb, oub, oesize, &
        int(ix_start, c_int), int(ix_end, c_int), int(direction, c_int), ixbp, tsp)
    finished = (cfin /= 0)
    if (present(ix_branch)) ix_branch = cixb
    if (present(track_state)) track_state = cts
  end subroutine

  subroutine bmad_hook_register_track_many(cfun) bind(c, name='bmad_hook_register_track_many')
    type(c_funptr), value :: cfun
    cb_track_many = cfun
    if (c_associated(cfun)) then
      track_many_hook_ptr => tramp_track_many
    else
      track_many_hook_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track1_postprocess (start_orb, ele, param, end_orb)
  !----------------------------------------------------------------------------

  subroutine tramp_track1_postprocess(start_orb, ele, param, end_orb)
    type(coord_struct) :: start_orb, end_orb
    type(ele_struct) :: ele
    type(lat_param_struct) :: param
    procedure(ci_track1_postprocess), pointer :: cproc
    call c_f_procpointer(cb_track1_postprocess, cproc)
    call cproc(cloc_coord(start_orb), cloc_ele(ele), cloc_lat_param(param), cloc_coord(end_orb))
  end subroutine

  subroutine bmad_hook_register_track1_postprocess(cfun) &
      bind(c, name='bmad_hook_register_track1_postprocess')
    type(c_funptr), value :: cfun
    cb_track1_postprocess = cfun
    if (c_associated(cfun)) then
      track1_postprocess_ptr => tramp_track1_postprocess
    else
      track1_postprocess_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)
  !----------------------------------------------------------------------------

  subroutine tramp_track1_preprocess(start_orb, ele, param, err_flag, finished, radiation_included, track)
    type(coord_struct) :: start_orb
    type(ele_struct), target :: ele
    type(lat_param_struct) :: param
    logical :: err_flag, finished, radiation_included
    type(track_struct), optional :: track
    procedure(ci_track1_preprocess), pointer :: cproc
    integer(c_int) :: ce, cf, cr
    type(c_ptr) :: tp
    call c_f_procpointer(cb_track1_preprocess, cproc)
    ! err_flag/finished/radiation_included are hook outputs; Bmad passes err_flag
    ! uninitialized here, so default them to false (see tramp_track1_bunch). Without
    ! this, an observer returning None inherits a garbage err_flag and track1 aborts
    ! before tracking the element (and before track1_postprocess).
    ce = 0
    cf = 0
    cr = 0
    if (present(track)) then
      tp = cloc_track(track)
    else
      tp = c_null_ptr
    end if
    call cproc(cloc_coord(start_orb), cloc_ele(ele), cloc_lat_param(param), ce, cf, cr, tp)
    err_flag = (ce /= 0)
    finished = (cf /= 0)
    radiation_included = (cr /= 0)
  end subroutine

  subroutine bmad_hook_register_track1_preprocess(cfun) &
      bind(c, name='bmad_hook_register_track1_preprocess')
    type(c_funptr), value :: cfun
    cb_track1_preprocess = cfun
    if (c_associated(cfun)) then
      track1_preprocess_ptr => tramp_track1_preprocess
    else
      track1_preprocess_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track1_spin_custom (start_orb, ele, param, end_orb, err_flag, make_quaternion)
  !----------------------------------------------------------------------------

  subroutine tramp_track1_spin_custom(start_orb, ele, param, end_orb, err_flag, make_quaternion)
    type(coord_struct) :: start_orb, end_orb
    type(ele_struct) :: ele
    type(lat_param_struct) :: param
    logical :: err_flag
    logical, optional :: make_quaternion
    procedure(ci_track1_spin_custom), pointer :: cproc
    integer(c_int) :: ce
    integer(c_int), target :: cmq
    type(c_ptr) :: mqp
    call c_f_procpointer(cb_track1_spin_custom, cproc)
    ! err_flag is a hook output; default it to false (see tramp_track1_bunch).
    ! make_quaternion is an input the hook reads, so keep its incoming value.
    ce = 0
    if (present(make_quaternion)) then
      cmq = merge(1_c_int, 0_c_int, make_quaternion)
      mqp = c_loc(cmq)
    else
      mqp = c_null_ptr
    end if
    call cproc(cloc_coord(start_orb), cloc_ele(ele), cloc_lat_param(param), cloc_coord(end_orb), ce, mqp)
    err_flag = (ce /= 0)
    if (present(make_quaternion)) make_quaternion = (cmq /= 0)
  end subroutine

  subroutine bmad_hook_register_track1_spin_custom(cfun) &
      bind(c, name='bmad_hook_register_track1_spin_custom')
    type(c_funptr), value :: cfun
    cb_track1_spin_custom = cfun
    if (c_associated(cfun)) then
      track1_spin_custom_ptr => tramp_track1_spin_custom
    else
      track1_spin_custom_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! track1_wake_hook (bunch, ele, finished)
  !----------------------------------------------------------------------------

  subroutine tramp_track1_wake(bunch, ele, finished)
    type(bunch_struct) :: bunch
    type(ele_struct) :: ele
    logical :: finished
    procedure(ci_track1_wake), pointer :: cproc
    integer(c_int) :: cf
    call c_f_procpointer(cb_track1_wake, cproc)
    ! finished is a hook output; default it to false (see tramp_track1_bunch).
    cf = 0
    call cproc(cloc_bunch(bunch), cloc_ele(ele), cf)
    finished = (cf /= 0)
  end subroutine

  subroutine bmad_hook_register_track1_wake(cfun) bind(c, name='bmad_hook_register_track1_wake')
    type(c_funptr), value :: cfun
    cb_track1_wake = cfun
    if (c_associated(cfun)) then
      track1_wake_hook_ptr => tramp_track1_wake
    else
      track1_wake_hook_ptr => null()
    end if
  end subroutine

  !----------------------------------------------------------------------------
  ! wall_hit_handler_custom (orb, ele, s)
  !----------------------------------------------------------------------------

  subroutine tramp_wall_hit(orb, ele, s)
    type(coord_struct) :: orb
    type(ele_struct) :: ele
    real(rp) :: s
    procedure(ci_wall_hit_handler_custom), pointer :: cproc
    call c_f_procpointer(cb_wall_hit, cproc)
    call cproc(cloc_coord(orb), cloc_ele(ele), real(s, c_double))
  end subroutine

  subroutine bmad_hook_register_wall_hit(cfun) bind(c, name='bmad_hook_register_wall_hit_handler_custom')
    type(c_funptr), value :: cfun
    cb_wall_hit = cfun
    if (c_associated(cfun)) then
      wall_hit_handler_custom_ptr => tramp_wall_hit
    else
      wall_hit_handler_custom_ptr => null()
    end if
  end subroutine

end module bmad_c_hook_interface
