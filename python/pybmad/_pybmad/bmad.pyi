"""Bmad routines"""

from collections.abc import Sequence
from typing import overload

import _pybmad


class BmadHooks:
    """
    Registry of Bmad tracking-hook callbacks.

    Assign a callable to a property to install a hook, assign None to clear it, and read the property to get the current callback (or None). Proxy/array arguments are live, non-owning views valid only for the duration of the call; do not stash them. Exceptions raised in a callback are reported and swallowed (they never propagate into Fortran).
    """

    @property
    def time_runge_kutta_periodic_kick(self) -> object:
        """
        Called during time Runge-Kutta tracking to apply a periodic (e.g. RF) kick.

        Signature:
            fn(orbit: CoordStruct, ele: EleStruct, param: LatParamStruct,
               stop_time: float, init_needed: int) -> tuple[float, int] | None

        Return ``(stop_time, init_needed)`` to update them, or None to leave unchanged.
        """

    @time_runge_kutta_periodic_kick.setter
    def time_runge_kutta_periodic_kick(self, arg: object | None) -> None: ...

    @property
    def track1_bunch(self) -> object:
        """
        Called by track1_bunch before the standard single-element bunch tracking.
        Return finished=True to have the callback fully replace it.

        Signature:
            fn(bunch: BunchStruct, ele: EleStruct, err: bool,
               centroid: CoordStructArray1D | None, direction: int | None,
               finished: bool, bunch_track: BunchTrackStruct | None)
               -> tuple[bool, bool] | None

        Return ``(err, finished)``, or None to leave unchanged.
        """

    @track1_bunch.setter
    def track1_bunch(self, arg: object | None) -> None: ...

    @property
    def track1_custom(self) -> object:
        """
        Called by track1 to track an element whose tracking_method is `custom`
        (or a custom element). The callback performs the tracking: the first argument
        is in/out -- mutate it in place to the exit coordinates.

        Signature:
            fn(orbit: CoordStruct, ele: EleStruct, param: LatParamStruct,
               err_flag: bool, finished: bool, track: TrackStruct | None)
               -> tuple[bool, bool] | None

        Return ``(err_flag, finished)``, or None to leave unchanged.
        """

    @track1_custom.setter
    def track1_custom(self, arg: object | None) -> None: ...

    @property
    def track_many(self) -> object:
        """
        Called at the start of track_many (tracking through a range of elements).
        Return finished=True to have the callback fully replace the tracking.

        Signature:
            fn(finished: bool, lat: LatStruct, orbit: CoordStructArray1D,
               ix_start: int, ix_end: int, direction: int,
               ix_branch: int | None, track_state: int | None) -> bool | None

        `orbit` is a live array view. Return finished, or None to leave unchanged.
        """

    @track_many.setter
    def track_many(self, arg: object | None) -> None: ...

    @property
    def track1_postprocess(self) -> object:
        """
        Called after every track1, once an element has been tracked. `end_orb` is
        a live proxy -- mutate it to adjust the exit coordinates.

        Signature:
            fn(start_orb: CoordStruct, ele: EleStruct, param: LatParamStruct,
               end_orb: CoordStruct) -> None
        """

    @track1_postprocess.setter
    def track1_postprocess(self, arg: object | None) -> None: ...

    @property
    def track1_preprocess(self) -> object:
        """
        Called at the start of every track1, before the element is tracked.
        Return finished=True to have the callback fully replace the tracking.

        Signature:
            fn(start_orb: CoordStruct, ele: EleStruct, param: LatParamStruct,
               err_flag: bool, finished: bool, radiation_included: bool,
               track: TrackStruct | None) -> tuple[bool, bool, bool] | None

        Return ``(err_flag, finished, radiation_included)``, or None to leave unchanged.
        """

    @track1_preprocess.setter
    def track1_preprocess(self, arg: object | None) -> None: ...

    @property
    def track1_spin_custom(self) -> object:
        """
        Called by track1 to track spin when spin_tracking_method is `custom`.

        Signature:
            fn(start_orb: CoordStruct, ele: EleStruct, param: LatParamStruct,
               end_orb: CoordStruct, err_flag: bool, make_quaternion: bool | None)
               -> bool | tuple[bool, bool] | None

        Return err_flag, or ``(err_flag, make_quaternion)``, or None to leave unchanged.
        """

    @track1_spin_custom.setter
    def track1_spin_custom(self, arg: object | None) -> None: ...

    @property
    def track1_wake(self) -> object:
        """
        Called during bunch tracking to apply wakefields for an element.
        Return finished=True to have the callback fully replace the standard wake.

        Signature:
            fn(bunch: BunchStruct, ele: EleStruct, finished: bool) -> bool | None

        Return finished, or None to leave unchanged.
        """

    @track1_wake.setter
    def track1_wake(self, arg: object | None) -> None: ...

    @property
    def wall_hit_handler_custom(self) -> object:
        """
        Called during Runge-Kutta / time tracking when a particle hits the chamber
        wall, at longitudinal position `s`.

        Signature:
            fn(orb: CoordStruct, ele: EleStruct, s: float) -> None
        """

    @wall_hit_handler_custom.setter
    def wall_hit_handler_custom(self, arg: object | None) -> None: ...

hooks: BmadHooks = ...

class AbMultipoleKick:
    """ab_multipole_kick return type"""

    @property
    def kx(self) -> float: ...

    @property
    def ky(self) -> float: ...

    @property
    def dk(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ab_multipole_kick(a: float, b: float, n: int, ref_species: int, ele_orientation: int, coord: _pybmad.CoordStruct, pole_type: int | None = None, scale: float | None = None) -> AbMultipoleKick:
    """
    Subroutine to put in the kick due to an ab_multipole.

    Parameters
    ----------
    a : float
        Multipole skew component.

    b : float
        Multipole normal component.

    n : int
        Multipole order.

    ref_species : int
        Reference species.

    ele_orientation : int
        Element orientation +1 = normal, -1 = reversed, 0 = Ignore orientation and tracking species (used with
        pole_type = magnetic$).

    coord : CoordStruct
        Particle position and direction of travel.

    pole_type : int, optional
        Type of multipole. magnetic$ (default) or electric$.

    scale : float, optional
        Factor to scale the kicks. Default is 1. For pole_type = electric$, set scale to the longitudinal length
        of the field region.

    Returns
    -------
    kx : float
        X kick.

    ky : float
        Y kick.

    dk : 2D array of float (shape: 2,2), optional
        Kick derivative: dkick(x,y)/d(x,y).
    """

def ab_multipole_kicks(an: _pybmad.RealArray1D, bn: _pybmad.RealArray1D, ix_pole_max: int, ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, pole_type: int | None = None, scale: float | None = None, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Routine to put in the kick due to ab_multipole components in an element.
    The kick will be corrected for the orientation of the element and the particle direction of travel.
    Any difference between element p0c and orbit%p0c will be taken into account.

    Also see the multipole_kicks routine.

    Parameters
    ----------
    an : 1D array of float
        Skew multipole strengths.

    bn : 1D array of float
        Normal multipole strengths.

    ix_pole_max : int
        Maximum pole index.

    ele : EleStruct
        Lattice element containing the multipoles.

    orbit : CoordStruct
        Particle position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Kicked particle.

    pole_type : int, optional
        Type of multipole. magnetic$ (default) or electric$.

    scale : float, optional
        Factor to scale the kicks. Default is 1. For pole_type = electric$, set scale to the longitudinal length
        of the field region

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the multipole.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including multipole.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def absolute_photon_position(e_orb: _pybmad.CoordStruct, photon_orb: _pybmad.CoordStruct) -> None:
    """
    $OMP THREADPRIVATE(pinit_E_rel_target, pinit_va_E_rel, pinit_va_gamma, &
    $OMP                pinit_va_vert_angle, pinit_va_invert, pinit_alpha)

     Subroutine absolute_photon_position (e_orb, photon_orb)

     Routine to calculate the photon phase space coordinates given:
       1) The phase space coords of the emitting charged particle and
       2) The photon phase space coords relative to the emitting particle.
          The photon (x, y, z) position is ignored (it is assumed the photon is emitted at
          the charged particle position) and only the photon's (vx, vy, vz) velocity matters.

    Parameters
    ----------
    e_orb : CoordStruct
        charged particle position.

    photon_orb : CoordStruct
        Photon position relative to e_orb.
        This parameter is an input/output and is modified in-place.
        As an output, photon_orb: Absolute photon position.
    """

def absolute_time_tracking(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine absolute_time_tracking

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    Returns
    -------
    is_abs_time : bool
        True if absolute time tracking is needed.
    """

def ac_kicker_amp(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, true_time: float | None = None) -> float:
    """
    Wrapper for Fortran routine ac_kicker_amp

    Parameters
    ----------
    ele : EleStruct
        ac_kicker element.

    orbit : CoordStruct
        Contains the time to evaluate the amplitude at.

    true_time : float, optional
        The actual time. Normally this time is calculated using orbit.t or orbit.vec(5) but sometimes it is
        convenient to be able to override this. For example, time_runge_kutta uses this.

    Returns
    -------
    ac_amp : float
        Amplitude. Will be set to 1 if the element is not an ac_kicker.
    """

class ActionToXyz:
    """action_to_xyz return type"""

    @property
    def X(self) -> list[float]: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def action_to_xyz(ring: _pybmad.LatStruct, ix: int, J: Sequence[float]) -> ActionToXyz:
    """
    Given the normal mode invariants and phases J of a particle, returns the canonical coordinates.

    The J vector looks like:
    J = (sqrt(2Ja)cos(phia), -sqrt(2Ja)sin(phia), sqrt(2Jb)cos(phib), -sqrt(2Jb)sin(phib), sqrt(2Jc)cos(phic), -sqrt(2Jc)sin(phic))

    X is obtained from:
    X = N . J
    Where N is from the Eigen decomposition of the 1-turn transfer matrix.

    Parameters
    ----------
    ring : LatStruct
        lattice

    ix : int
        element index at which to calculate J

    J : 1D array of float (shape: 1:6)
        Vector containing normal mode invariants and phases

    Returns
    -------
    X : 1D array of float (shape: 1:6)
        canonical phase space coordinates of the particle

    err_flag : bool
        Set to true on error.  Often means Eigen decomposition failed.
    """

def add_lattice_control_structs(ele: _pybmad.EleStruct, n_add_slave: int | None = None, n_add_lord: int | None = None, n_add_slave_field: int | None = None, n_add_lord_field: int | None = None, add_at_end: bool | None = None) -> None:
    """
    Wrapper for Fortran routine add_lattice_control_structs

    Parameters
    ----------
    ele : EleStruct
        Lord or slave element that needs extra control elements.

    n_add_slave : int, optional
        Number of field slaves to add to lord. Default is zero.

    n_add_lord : int, optional
        Number of field lords to add to slave. Default is zero.

    n_add_slave_field : int, optional
        Number of field slaves to add to lord. Default is zero.

    n_add_lord_field : int, optional
        Number of field lords to add to slave. Default is zero.

    add_at_end : bool, optional
        Used when n_add_slave or n_add_slave_field is non-zero. If True then new space is added at the end of the
        array. If False then new space is added at the front of the array. Default is True.
    """

def add_ptc_layout_to_list(branch_ptc_info: _pybmad.PtcBranch1Struct, layout_end: _pybmad.Layout) -> None:
    """
    Routine to add a layout the a list of layouts.

    Parameters
    ----------
    branch_ptc_info : PtcBranch1Struct
        List of layouts
        This parameter is an input/output and is modified in-place.
        As an output, branch_ptc_info: Updated list.

    layout_end : Layout
        ptc layout
    """

class AddSuperimpose:
    """add_superimpose return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def super_ele_out(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def add_superimpose(lat: _pybmad.LatStruct, super_ele_in: _pybmad.EleStruct, ix_branch: int, save_null_drift: bool | None = None, create_jumbo_slave: bool | None = None, ix_insert: int | None = None, mangle_slave_names: bool | None = None, wrap: bool | None = None) -> AddSuperimpose:
    """
    save_null_drift, create_jumbo_slave, ix_insert, mangle_slave_names, wrap)

    Routine to superimpose an element. If the element can be inserted
    into the lat without making a super_lord element then this will be done.

    Note: This routine, since it handles only one superposition, is not sufficient for
      superposition in a multipass region. For historical reasons, the extra code needed
      is buried in the parser_add_superimpose code. If you need to do multipass superpositions
      please contact David Sagan and this situation will be rectified.

    Note: Bookkeeping like recalculating reference energies and recalculating transfer matrices
      is *not* done by this routine.

    Parameters
    ----------
    lat : LatStruct
        Lat to modify.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Modified lat.

    super_ele_in : EleStruct
        Element to superimpose.

    ix_branch : int
        Branch index to put element.

    save_null_drift : bool, optional
        Save a copy of a drift to be split as a null_ele? This is useful if further superpositions might use this
        drift as a reference element. After all superpositions are done, remove_eles_from_lat can be called to
        remove all null_eles. Default is False.

    create_jumbo_slave : bool, optional
        Default is False. If True then super_slaves that are created that have super_ele_in as their super_lord
        are em_field elements.

    ix_insert : int, optional
        If present and positive, and super_ele_in has zero length, use ix_insert as the index to insert
        super_ele_in at. ix_insert is useful when superposing next to another element that has zero or negative
        length (EG a patch) and you want to make sure that the superimposed element is on the correct side of the
        element.

    mangle_slave_names : bool, optional
        If True (default), adjust slave names appropriately. Name mangeling can take time so bmad_parser will do
        this all at once at the end.

    wrap : bool, optional
        If True (default), and if the superimposed element has an end that extends beyond the starting or ending
        edge of the lattice, wrap the element around the lattice so that the beginning portion of the element is
        at the lattice ending edge and the rest of the element is at the lattice start edge. If wrap = False, and
        the superimposed element has an end that extends beyound a lattice edge, extend the lattice to
        accommodate.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise

    super_ele_out : EleStruct, optional
        Pointer to the super element in the lattice.
    """

def add_this_multipass(lat: _pybmad.LatStruct, m_slaves: _pybmad.LatEleLocStructArray1D, lord_in: _pybmad.EleStruct | None = None) -> None:
    """
    Wrapper for Fortran routine add_this_multipass

    Parameters
    ----------
    lat : LatStruct

    m_slaves : 1D array of LatEleLocStruct

    lord_in : EleStruct, optional
    """

def add_this_name_to_list(ele: _pybmad.EleStruct, names: _pybmad.CharacterAlloc1D, an_indexx: _pybmad.IntAlloc1D, n_names: int, ix_match: int, has_been_added: bool, named_eles: _pybmad.ElePointerStructAlloc1D) -> None:
    """
    Wrapper for Fortran routine add_this_name_to_list

    Parameters
    ----------
    ele : EleStruct

    names : 1D array of str

    an_indexx : 1D array of int

    n_names : int

    ix_match : int

    has_been_added : bool

    named_eles : 1D array of ElePointerStruct
    """

def add_this_taylor_term(ele: _pybmad.EleStruct, i_out: int, coef: float, expn: Sequence[int]) -> None:
    """
    Subroutine used by bmad_parser and bmad_parser2 to parse the input file.
    This subroutine is not intended for general use.
    """

def adjust_super_slave_names(lat: _pybmad.LatStruct, ix1_lord: int, ix2_lord: int, first_time: bool | None = None) -> None:
    """
    Routine to adjust the names of the slaves.
    This routine is used by add_superimpose and is not meant for general use.
    """

def allocate_branch_array(lat: _pybmad.LatStruct, upper_bound: int) -> None:
    """
    Wrapper for Fortran routine allocate_branch_array

    Parameters
    ----------
    lat : LatStruct

    upper_bound : int
        Desired upper bound.
    """

def allocate_grid_field(g_field: _pybmad.GridFieldStructArray1D, n_gf: int) -> None:
    """
    Wrapper for Fortran routine allocate_grid_field

    Parameters
    ----------
    g_field : 1D array of GridFieldStruct

    n_gf : int
    """

def allocate_lat_ele_array(lat: _pybmad.LatStruct, upper_bound: int | None = None, ix_branch: int | None = None, do_ramper_slave_setup: bool | None = None) -> None:
    """
    Wrapper for Fortran routine allocate_lat_ele_array

    Parameters
    ----------
    lat : LatStruct
        Lattice with element array.

    upper_bound : int, optional
        Optional desired upper bound. Default: 1.3*ubound(ele(:)) or 10 if ele is not allocated.

    ix_branch : int, optional
        Branch index. Default is 0.

    do_ramper_slave_setup : bool, optional
        Default False. If true, setup ramper slaves. Generally this needs to be done if reallocating with a fully
        formed lattice.
    """

def allocate_plat(plat: _pybmad.ParserLatStruct, n_ele_max: int) -> None:
    """
    Subroutine to allocate allocatable array sizes.
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

def angle_between_polars(polar1: _pybmad.SpinPolarStruct, polar2: _pybmad.SpinPolarStruct) -> float:
    """
    Wrapper for Fortran routine angle_between_polars

    Parameters
    ----------
    polar1 : SpinPolarStruct
        (spin_polar_struct)

    polar2 : SpinPolarStruct
        (spin_polar_struct)

    Returns
    -------
    angle : float
        Angle between the polar vectors
    """

def angle_to_canonical_coords(orbit: _pybmad.CoordStruct, coord_type: str | None = None) -> None:
    r"""
    Wrapper for Fortran routine angle_to_canonical_coords

    Parameters
    ----------
    orbit : CoordStruct
        Orbit in angular coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit in canonical coordinates.

    coord_type : str, optional
        Angular coordinates type \'\' (default): (x, x' = dx/ds, y, y' = dy/ds, z, pz) 'ZGOUBI':     (x, x' = dx/ds,
        y, y' = dy/ds, dt = -z / (beta * c), pz)
    """

def aperture_bookkeeper(ele: _pybmad.EleStruct) -> None:
    """
    Routine to calculate aperture limits when ele%attribute_type is set to auto_aperture$

    Parameters
    ----------
    ele : EleStruct
        Element with aperture.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with apertures set.
    """

def apply_all_rampers(lat: _pybmad.LatStruct) -> bool:
    """
    Wrapper for Fortran routine apply_all_rampers

    Parameters
    ----------
    lat : LatStruct
        Lattice.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with rampers applied.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def apply_element_edge_kick(orb: _pybmad.CoordStruct, fringe_info: _pybmad.FringeFieldInfoStruct, track_ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, track_spin: bool, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None, rf_time: float | None = None, apply_sol_fringe: bool | None = None) -> None:
    """
    Wrapper for Fortran routine apply_element_edge_kick

    Parameters
    ----------
    orb : CoordStruct
        Starting coords in element reference frame.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Coords after application of the edge fringe field.

    fringe_info : FringeFieldInfoStruct
        Fringe information.

    track_ele : EleStruct
        Element being tracked through. Is different from fringe_info.hard_ele when there are superpositions and
        track_ele can be a super_slave of fringe_info.hard_ele.

    param : LatParamStruct
        lattice parameters.

    track_spin : bool
        Track the spin?

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before fringe.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including fringe.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    rf_time : float, optional
        RF clock time. If not present then the time will be calculated using the standard algorithm.

    apply_sol_fringe : bool, optional
        Apply the solenoid fringe kick? Default is True.
    """

def apply_energy_kick(dE: float, orbit: _pybmad.CoordStruct, ddE_dr: Sequence[float], mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine apply_energy_kick

    Parameters
    ----------
    dE : float
        Energy change

    orbit : CoordStruct
        Beginning coordinates
        This parameter is an input/output and is modified in-place.
        As an output, orbit: coordinates with added dE energy kick.

    ddE_dr : 1D array of float (shape: 2)
        real(rp), Derivatives of dE [ddE_dx, ddE_dy].

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before fringe.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including energy kick.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def apply_fft_3d_kicks(csr: _pybmad.CsrStruct, particle: _pybmad.CoordStructArray1D) -> None:
    """
    Routine to apply FFT-based 3D space charge kicks to particles.
    Deposits charge on a mesh, solves for the field, interpolates back, and applies kicks.

    Parameters
    ----------
    csr : CsrStruct
        Contains mesh, position arrays, and tracking parameters.

    particle : 1D array of CoordStruct
        Particles to kick.
        This parameter is an input/output and is modified in-place.
        As an output, particle: Particles with kicks applied.
    """

def apply_patch_to_ptc_fibre(ele: _pybmad.EleStruct) -> None:
    """
    Routine to take the patch parameters from a Bmad patch element and
    transfer them to the associated PTC fibre.

    Parameters
    ----------
    ele : EleStruct
        Patch element.
    """

def apply_rampers_to_slave(slave: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine apply_rampers_to_slave

    Parameters
    ----------
    slave : EleStruct
        Element to apply ramper elements to.

    Returns
    -------
    err_flag : bool
        Set true if there is an error. False otherwise.
    """

def array_re_str(arr: _pybmad.RealArray1D, parens_in: str | None = None) -> str:
    """
    Wrapper for Fortran routine array_re_str

    Parameters
    ----------
    arr : 1D array of float

    parens_in : str, optional

    Returns
    -------
    str_out : str
    """

def astra_max_field_reference(pt0: _pybmad.GridFieldPt1Struct, ele: _pybmad.EleStruct) -> float:
    """
    Wrapper for Fortran routine astra_max_field_reference

    Parameters
    ----------
    pt0 : GridFieldPt1Struct

    ele : EleStruct

    Returns
    -------
    field_value : float
    """

def at_this_ele_end(now_at: int, where_at: int) -> bool:
    """
    Wrapper for Fortran routine at_this_ele_end

    Parameters
    ----------
    now_at : int
        Which end is under consideration: entrance_end$, exit_end$, surface$, or in_between$.

    where_at : int
        Which ends have the aperture or fringe field: entrance_end$, exit_end$, continuous$, both_ends$,
        no_aperture$, surface$, wall_transition$.

    Returns
    -------
    is_at_this_end : bool
        True if at this end. False otherwise.
    """

def attribute_bookkeeper(ele: _pybmad.EleStruct, force_bookkeeping: bool | None = None) -> None:
    """
    Wrapper for Fortran routine attribute_bookkeeper

    Parameters
    ----------
    ele : EleStruct
        Element with attributes
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with self-consistant attributes.

    force_bookkeeping : bool, optional
        If present and True then force attribute bookkeeping to be done independent of the state of
        ele.bookkeeping_stat.attributes.
    """

class AttributeFree1:
    """attribute_free1 return type"""

    @property
    def why_not_free(self) -> int: ...

    @property
    def free(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

@overload
def attribute_free(ix_ele: int, attrib_name: str, lat: _pybmad.LatStruct, err_print_flag: bool | None = None, except_overlay: bool | None = None, dependent_attribs_free: bool | None = None) -> AttributeFree1:
    """
    Overloaded function for:
      Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)
      Function attribute_free2 (ele, attrib_name, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)
      Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)

    Routine to check if an attribute is free to vary.

    Attributes that cannot be changed directly include super_slave attributes (since
    these attributes are controlled by their super_lords) and attributes that
    are controlled by an overlay.

    Also dependent variables such as the angle of a bend cannot be
      freely variable.

    Parameters
    ----------
    ix_ele : int
        Index of element in element array.

    attrib_name : str
        Name of the attribute. Assumed upper case.

    lat : LatStruct
        Lattice structure.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message if attribute is not free.

    except_overlay : bool, optional
        If present and True then an attribute that is controlled by an overlay will be treated as free. This is
        used by, for example, the create_overlay routine.

    dependent_attribs_free : bool, optional
        If present and True then mark as free attributes that are dependent. For example, if ele.field_master = F,
        b1_field is dependent upon k1. Default is False. Also fixer Twiss/dispersion/orbit attributes are
        considered "dependent".

    Returns
    -------
    free : bool
        Set True if attribtute not found or attriubte cannot be changed directly.

    why_not_free : int, optional
        Possibilities are: field_master_dependent$  -> Dependent due to setting of ele.field_master. dependent$
        -> Not field_master_dependent$ but value is dependent upon the value of other attributes. does_not_exist$
        -> Attribute name is unrecognized or does not exist for the type of element. overlay_slave$           ->
        Attribute is controlled by an overlay lord. super_slave$             -> Attribute is controlled by
        element's super_lord. multipass_slave$         -> Attribute is controlled by element's multipass_lord.
    """

@overload
def attribute_free(ele: _pybmad.EleStruct, attrib_name: str, err_print_flag: bool | None = None, except_overlay: bool | None = None, dependent_attribs_free: bool | None = None) -> AttributeFree2:
    """
    Overloaded function for:
      Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)
      Function attribute_free2 (ele, attrib_name, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)
      Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)

    Routine to check if an attribute is free to vary.

    Attributes that cannot be changed directly include super_slave attributes (since
    these attributes are controlled by their super_lords) and attributes that
    are controlled by an overlay.

    Also dependent variables such as the angle of a bend cannot be
      freely variable.

    Parameters
    ----------
    ele : EleStruct
        Element containing the attribute

    attrib_name : str
        Name of the attribute. Assumed upper case.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message if attribute is not free.

    except_overlay : bool, optional
        If present and True then an attribute that is controlled by an overlay will be treated as free. This is
        used by, for example, the create_overlay routine.

    dependent_attribs_free : bool, optional
        If present and True then mark as free attributes that are dependent. For example, if ele.field_master = F,
        b1_field is dependent upon k1. Default is False. Also fixer Twiss/dispersion/orbit attributes are
        considered "dependent".

    Returns
    -------
    free : bool
        Set True if attribtute not found or attriubte cannot be changed directly.

    why_not_free : int, optional
        Possibilities are: field_master_dependent$  -> Dependent due to setting of ele.field_master. dependent$
        -> Not field_master_dependent$ but value is dependent upon the value of other attributes. does_not_exist$
        -> Attribute name is unrecognized or does not exist for the type of element. overlay_slave$           ->
        Attribute is controlled by an overlay lord. super_slave$             -> Attribute is controlled by
        element's super_lord. multipass_slave$         -> Attribute is controlled by element's multipass_lord.
    """

@overload
def attribute_free(ix_ele: int, ix_branch: int, attrib_name: str, lat: _pybmad.LatStruct, err_print_flag: bool | None = None, except_overlay: bool | None = None, dependent_attribs_free: bool | None = None) -> AttributeFree3:
    """
    Overloaded function for:
      Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)
      Function attribute_free2 (ele, attrib_name, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)
      Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, err_print_flag,
                                   except_overlay, dependent_attribs_free, why_not_free) result (free)

    Routine to check if an attribute is free to vary.

    Attributes that cannot be changed directly include super_slave attributes (since
    these attributes are controlled by their super_lords) and attributes that
    are controlled by an overlay.

    Also dependent variables such as the angle of a bend cannot be
      freely variable.

    Parameters
    ----------
    ix_ele : int
        Index of element in element array.

    ix_branch : int
        Branch index of element.

    attrib_name : str
        Name of the attribute. Assumed upper case.

    lat : LatStruct
        Lattice structure.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message if attribute is not free.

    except_overlay : bool, optional
        If present and True then an attribute that is controlled by an overlay will be treated as free. This is
        used by, for example, the create_overlay routine.

    dependent_attribs_free : bool, optional
        If present and True then mark as free attributes that are dependent. For example, if ele.field_master = F,
        b1_field is dependent upon k1. Default is False. Also fixer Twiss/dispersion/orbit attributes are
        considered "dependent".

    Returns
    -------
    free : bool
        Set True if attribtute not found or attriubte cannot be changed directly.

    why_not_free : int, optional
        Possibilities are: field_master_dependent$  -> Dependent due to setting of ele.field_master. dependent$
        -> Not field_master_dependent$ but value is dependent upon the value of other attributes. does_not_exist$
        -> Attribute name is unrecognized or does not exist for the type of element. overlay_slave$           ->
        Attribute is controlled by an overlay lord. super_slave$             -> Attribute is controlled by
        element's super_lord. multipass_slave$         -> Attribute is controlled by element's multipass_lord.
    """

class AttributeFree2:
    """attribute_free2 return type"""

    @property
    def why_not_free(self) -> int: ...

    @property
    def free(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

class AttributeFree3:
    """attribute_free3 return type"""

    @property
    def why_not_free(self) -> int: ...

    @property
    def free(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

class AttributeIndex1:
    """attribute_index1 return type"""

    @property
    def full_name(self) -> str: ...

    @property
    def attrib_index(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

@overload
def attribute_index(ele: _pybmad.EleStruct, name: str, can_abbreviate: bool | None = None, print_error: bool | None = None) -> AttributeIndex1:
    """
    Function to return the index of a attribute for a given BMAD element type
    and the name of the attribute. Abbreviations are by default permitted but must be at
    least 3 characters. Exception: overlay and group varialbe names may not
    be abbreviated.

    This routine is an overloaded name for:
      attribute_index1 (ele, name, full_name, can_abbreviate, print_error) result (attrib_index)
      attribute_index2 (key, name, full_name, can_abbreviate, print_error) result (attrib_index)

    Parameters
    ----------
    ele : EleStruct
        attribute_index will restrict the name search to valid attributes of the given element.

    name : str
        Attribute name. Must be uppercase.

    can_abbreviate : bool, optional
        Can abbreviate names? Default is True.

    print_error : bool, optional
        Default True. If false, do not print error message.

    Returns
    -------
    attrib_index : int
        Index of the attribute. If the attribute name is not appropriate then 0 will be returned. ix -> k1$

    full_name : str, optional
        Non-abbreviated name.

    Notes
    -----
    If ele%key or key = 0 -> Entire name table will be searched.
    See also:
    has_attribute
    attribute_info
    attribute_name
    """

@overload
def attribute_index(key: int, name: str, can_abbreviate: bool | None = None, print_error: bool | None = None) -> AttributeIndex2:
    """
    Function to return the index of a attribute for a given BMAD element type
    and the name of the attribute. Abbreviations are by default permitted but must be at
    least 3 characters. Exception: overlay and group varialbe names may not
    be abbreviated.

    This routine is an overloaded name for:
      attribute_index1 (ele, name, full_name, can_abbreviate, print_error) result (attrib_index)
      attribute_index2 (key, name, full_name, can_abbreviate, print_error) result (attrib_index)

    Parameters
    ----------
    key : int
        Equivalent to ele.key.

    name : str
        Attribute name. Must be uppercase.

    can_abbreviate : bool, optional
        Can abbreviate names? Default is True.

    print_error : bool, optional
        Default True. If false, do not print error message.

    Returns
    -------
    attrib_index : int
        Index of the attribute. If the attribute name is not appropriate then 0 will be returned. ix -> k1$

    full_name : str, optional
        Non-abbreviated name.

    Notes
    -----
    If ele%key or key = 0 -> Entire name table will be searched.
    See also:
    has_attribute
    attribute_info
    attribute_name
    """

class AttributeIndex2:
    """attribute_index2 return type"""

    @property
    def full_name(self) -> str: ...

    @property
    def attrib_index(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def attribute_info(ele: _pybmad.EleStruct, ix_att: int) -> _pybmad.EleAttributeStruct:
    """
    Function to return the info structure associated with an attribute for
    a particular type of BMAD element.

    Parameters
    ----------
    ele : EleStruct

    ix_att : int
        Index of attribute (e.g. k1$)

    Returns
    -------
    attrib_info : EleAttributeStruct
        Info on this attribute. Note: .value is not set since this info is contained in the ele argument. Use
        pointer_to_attribute to access the attribute.
    """

@overload
def attribute_name(key: int, ix_att: int, show_private: bool | None = None) -> str:
    """
    Function to return the name of an attribute for a particular type of
    Bmad element.

    This routine is an overloaded name for:
      attribute_name1 (ele, ix_att, show_private) result (attrib_name)
      attribute_name2 (key, ix_att, show_private) result (attrib_name)

    Note: attribute_name (key, ix_att) is not able to handle overlay/group control variables.
    Use attributge_name (ele, ix_att) is this is needed.

    Parameters
    ----------
    key : int
        Key name of element type (e.g. sbend$, etc.)

    ix_att : int
        Index of attribute (e.g. k1$)

    show_private : bool, optional
        If False (default) return null_name$ for private attributes.

    Returns
    -------
    attrib_name : str
        Name of attribute. First character is a "!" if there is a problem. Will always be upper case (even with
        private attributes). = "!BAD ELE KEY"           .key is invalid = "!BAD INDEX"             ix_att is
        invalid (out of range). = "!NULL" (null_name$)     ix_att does not correspond to an attribute or is
        private. name -> "K1"
    """

@overload
def attribute_name(ele: _pybmad.EleStruct, ix_att: int, show_private: bool | None = None) -> str:
    """
    Function to return the name of an attribute for a particular type of
    Bmad element.

    This routine is an overloaded name for:
      attribute_name1 (ele, ix_att, show_private) result (attrib_name)
      attribute_name2 (key, ix_att, show_private) result (attrib_name)

    Note: attribute_name (key, ix_att) is not able to handle overlay/group control variables.
    Use attributge_name (ele, ix_att) is this is needed.

    Parameters
    ----------
    ele : EleStruct

    ix_att : int
        Index of attribute (e.g. k1$)

    show_private : bool, optional
        If False (default) return null_name$ for private attributes.

    Returns
    -------
    attrib_name : str
        Name of attribute. First character is a "!" if there is a problem. Will always be upper case (even with
        private attributes). = "!BAD ELE KEY"           .key is invalid = "!BAD INDEX"             ix_att is
        invalid (out of range). = "!NULL" (null_name$)     ix_att does not correspond to an attribute or is
        private. name -> "K1"
    """

def attribute_set_bookkeeping(ele: _pybmad.EleStruct, attrib_name: str, err_flag: bool, attrib_ptr: _pybmad.AllPointerStruct | None = None) -> None:
    """
    Wrapper for Fortran routine attribute_set_bookkeeping

    Parameters
    ----------
    ele : EleStruct
        Element containing the attribute

    attrib_name : str
        Name of the attribute. Must be upper case.

    err_flag : bool

    attrib_ptr : AllPointerStruct, optional
        Pointer to the attribute. The presence of this argument saves a small amount of time.
    """

def attribute_type(attrib_name: str, ele: _pybmad.EleStruct | None = None) -> int:
    """
    Routine to return the logical type of an attribute.

    A "switch" attribute is an attribute whose value corresponds to some string.
    For example, the "COUPLER_AT" attirbute with value 1 corresponds to "ENTRANCE_END", etc.

    A "struct" attribute is an attribute that is the name for a "structure". For example,
    CARTESIAN_MAP is the name of the structure hoding a Cartesian map.

    If attrib_name corresponds to a switch attribute, The routine switch_attrib_value_name can
    be used to print the name corresponding to the attribute's value.

    Note: The "storage type" of an attribute is different from the "logical type" returned by
    this routine. For example, the logical type of attribute "n_slice" is integer. However, the
    value of "n_slice" is stored as a real number in the ele_struct [in ele%value(n_slice$)].

    Parameters
    ----------
    attrib_name : str
        Name of the attribute. Must be upper case.

    ele : EleStruct, optional
        Element associated with the attribute. Needed if attrib_name can correspond to an overlay or group
        variable.

    Returns
    -------
    attrib_type : int
        Attribute type: is_string$, is_logical$, is_integer$, is_real$, is_switch$, is_struct$, is_species$ or
        invalid_name$ Note: An overlay or group variable will be marked invalid_name$ if ele is missing.
    """

def attribute_units(attrib_name: str, unrecognized_units: str | None = None) -> str:
    """
    Routine to return the units associated with an attribute.
    Example: attrib_units('P0C') -> 'eV'

    Parameters
    ----------
    attrib_name : str
        Name of the attribute. Must be upper case.

    unrecognized_units : str, optional
        String to use if the attribute name is unrecognized. Note: Non-real attributes (EG: 'TRACKING_METHOD') are
        not recognized. Default is ""

    Returns
    -------
    attrib_units : str
        Units associated with the attribute.
    """

def autoscale_phase_and_amp(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, scale_phase: bool | None = None, scale_amp: bool | None = None, call_bookkeeper: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine autoscale_phase_and_amp

    Parameters
    ----------
    ele : EleStruct
        RF element or e_gun.
        This parameter is an input/output and is modified in-place.
        As an output, ele: element with phase and amplitude adjusted.

    param : LatParamStruct
        lattice parameters

    scale_phase : bool, optional
        Scale the phase? See above.

    scale_amp : bool, optional
        Scale the amplitude? See above.

    call_bookkeeper : bool, optional
        Call lattice_bookkeeper at end? Default is True.

    Returns
    -------
    err_flag : bool
        Logical, Set true if there is an error. False otherwise.
    """

def average_twiss(frac1: float, twiss1: _pybmad.TwissStruct, twiss2: _pybmad.TwissStruct) -> _pybmad.TwissStruct:
    """
    Wrapper for Fortran routine average_twiss

    Parameters
    ----------
    frac1 : float
        Fraction of twiss1 to use in the average.

    twiss1 : TwissStruct
        Twiss parameters to average.

    twiss2 : TwissStruct

    Returns
    -------
    ave_twiss : TwissStruct
        Average twiss.
    """

def bane1(ele: _pybmad.EleStruct, coulomb_log: float, n_part: float) -> _pybmad.IbsStruct:
    """
    This is an implementation of equations 10-15 from "Intrabeam
    scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
    It is a high energy approximation of the Bjorken-Mtingwa IBS
    formulation.
    """

class BbiKick:
    """bbi_kick return type"""

    @property
    def nk(self) -> list[float]: ...

    @property
    def dnk(self) -> list[list[float]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bbi_kick(x: float, y: float, sigma: Sequence[float], linear_kick: bool | None = None) -> BbiKick:
    """
    Wrapper for Fortran routine bbi_kick

    Parameters
    ----------
    x : float
        X coordinate.

    y : float
        Y coordinate.

    sigma : 1D array of float (shape: 2)
        Beam (x,y) sigmas.

    linear_kick : bool, optional
        Default False. If present and True, kick and dnk are computed using the extrapolated kick from the linear
        region.

    Returns
    -------
    nk : 1D array of float (shape: 2)
        Normalized, dimensionless kick component. In terms of the the actual kick: nk = [kick_x / (xi_x * sigma_x
        / beta_x), kick_y / (xi_y * sigma_y / beta_y) nk = -4 * pi * [x/sigma_x, y/sigma_y] in the linear region

    dnk : 2D array of float (shape: 2,2)
        derivatives of nk. EG: dnk(2,1) = nk(2)/dx
    """

def bbi_slice_calc(ele: _pybmad.EleStruct, n_slice: int, z_slice: _pybmad.RealArray1D) -> None:
    """
    Wrapper for Fortran routine bbi_slice_calc

    Parameters
    ----------
    ele : EleStruct
        beambeam element

    n_slice : int
        Number of slices

    z_slice : 1D array of float
        Array of slice positions 1:n_slice. zero padded for indexes greater than n_slice
    """

def beam_envelope_ibs(sigma_mat: Sequence[Sequence[float]], tail_cut: bool, tau: float, energy: float, n_part: float, species: int) -> list[list[float]]:
    """
    This is a sigma matrix based IBS calculation.
    It takes the beam sigma matrix and returns a matrix with changes to the 2nd order
    moments due to IBS.

    Use ibs_mat to change the sigma matrix like this:
    sigma_matrix_updated = sigma_matrix + ibs_mat*element_length
    See subroutine transport_with_sr_and_ibs in this module.

    Parameters
    ----------
    sigma_mat : 2D array of float (shape: 6,6)
        beam sigma_matrix at element entrance

    tail_cut : bool
        If true, then apply tail cut to coulomb logarithm.

    tau : float
        horizontal betatron damping rate.  Needed if tail_cut is true.

    energy : float
        beam energy in eV

    n_part : float
        number of particles in the bunch

    species : int
        Partical species.

    Returns
    -------
    ibs_mat : 2D array of float (shape: 6,6)
        changes in 2nd order moments due to IBS are ibs_mat*element_length
    """

def beam_equal_beam(beam1: _pybmad.BeamStruct, beam2: _pybmad.BeamStruct) -> None:
    """
    Wrapper for Fortran routine beam_equal_beam

    Parameters
    ----------
    beam1 : BeamStruct

    beam2 : BeamStruct
    """

class BeamInitSetup:
    """beam_init_setup return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def beam_init_set(self) -> _pybmad.BeamInitStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def beam_init_setup(beam_init_in: _pybmad.BeamInitStruct, ele: _pybmad.EleStruct, species: int, modes: _pybmad.NormalModesStruct | None = None) -> BeamInitSetup:
    """
    Wrapper for Fortran routine beam_init_setup

    Parameters
    ----------
    beam_init_in : BeamInitStruct
        Input parameters

    ele : EleStruct

    species : int
        Beam particle species.

    modes : NormalModesStruct, optional
        Normal mode parameters.

    Returns
    -------
    beam_init_set : BeamInitStruct
        See above.

    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

class BeamTilts:
    """beam_tilts return type"""

    @property
    def angle_xy(self) -> float: ...

    @property
    def angle_xz(self) -> float: ...

    @property
    def angle_yz(self) -> float: ...

    @property
    def angle_xpz(self) -> float: ...

    @property
    def angle_ypz(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def beam_tilts(S: Sequence[Sequence[float]]) -> BeamTilts:
    """
    Given a 6x6 matrix of second-order moments, this routine returns
    the beam tilts.

    angle_xy is obtained from the projection of the beam envelop into the
    xy plane.  The angle is that between the major axis of the projected
    beam envelope and the +x axis.  Positive angles are measured towards the
    +y axis.

    angle_xz is obtained from the projection of the beam envelop into the
    xy plane.  The angle is that between the major axis of the projected beam envelope
    and the +z axis.  Positive angles are measured towards the +x axis.

    angle_yz is obtained from the projection of the beam envelop into the
    yz plane.  The angle is that between the major axis of the projected beam envelope
    and the +z axis.  Positive angles are measured towards the +y axis.

    Parameters
    ----------
    S : 2D array of float (shape: 6,6)
        matrix of second order moments of beam envelope

    Returns
    -------
    angle_xy : float
        transverse tilt of beam envelope

    angle_xz : float
        horizontal crabbing of beam envelope

    angle_yz : float
        vertical crabbing of beam envelope

    angle_xpz : float
        x-pz coupling

    angle_ypz : float
        y-pz coupling
    """

def beambeam_fibre_setup(ele: _pybmad.EleStruct) -> _pybmad.Fibre:
    """
    Routine to setup a fibre to handle the beambeam interaction.

    Parameters
    ----------
    ele : EleStruct
        Bmad beambeam element.

    Returns
    -------
    ptc_fibre : Fibre
        Corresponding PTC fibre.
    """

def bend_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orb: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None, track_spin: bool | None = None) -> None:
    """
    Subroutine to track through the edge field of an sbend.
    This routine is called by apply_element_edge_kick only.

    Parameters
    ----------
    ele : EleStruct
        SBend element.

    param : LatParamStruct
        Rel charge.

    particle_at : int
        first_track_edge$, or second_track_edge$.

    orb : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Coords after tracking.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before fringe.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including fringe.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    track_spin : bool, optional
        If True then track the spin through the edge fields. Default: False.
    """

def bend_exact_multipole_field(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct, local_ref_frame: bool, calc_dfield: bool | None = None, calc_potential: bool | None = None) -> _pybmad.EmFieldStruct:
    """
    Wrapper for Fortran routine bend_exact_multipole_field

    Parameters
    ----------
    ele : EleStruct
        Bend element.

    param : LatParamStruct
        Lattice branch parameters.

    orbit : CoordStruct
        particle position.

    local_ref_frame : bool
        Is the particle position in the local element ref frame (as opposed to the lab frame)?

    calc_dfield : bool, optional
        If present and True then calculate the field derivatives.

    calc_potential : bool, optional
        Calc electric and magnetic potentials? Default is false.

    Returns
    -------
    field : EmFieldStruct
        Field
    """

def bend_length_has_been_set(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine bend_length_has_been_set

    Parameters
    ----------
    ele : EleStruct
        Element to be checked.

    Returns
    -------
    is_set : bool
        Note: will be set True for non-bend elements.
    """

def bend_photon_e_rel_init(r_in: float | None = None) -> float:
    """
    Routine to convert a random number in the interval [0,1] to a photon energy.
    The photon probability spectrum is:
      P(E_rel) = (3 / (5 * Pi)) * Integral_{E_rel}^{Infty} K_{5/3}(x) dx
    Where
      P(E_rel)) = Probability of finding a photon at relative energy E_rel.
      E_rel     = Relative photon energy: E / E_crit, E_crit = Critical energy.
      K_{5/3}   = Modified Bessel function.

    Notice that the P(E) is not the same as the distribution radiation energy since
    the photons must be energy weighted.

    There is a cut-off built into the calculation so that E_rel will be in the
    range [0, 31.4]. The error in neglecting photons with E_rel > 31.4 translates
    to neglecting one photon for every 10^15 generated.
    If r_in is present:
      r_in = 0 => E_rel = 0
      r_in = 1 => E_rel = 31.4

    Parameters
    ----------
    r_in : float, optional
        Integrated probability in the range [0,1]. If not present, a random number will be used.

    Returns
    -------
    E_rel : float
        Relative photon energy E/E_crit.
    """

def bend_photon_energy_integ_prob(E_photon: float, g_bend: float, gamma: float) -> float:
    """
    Routine to find the integrated probability corresponding to emitting a photon
    from a bend in the range [0, E_photon].

    Parameters
    ----------
    E_photon : float
        Photon energy.

    g_bend : float
        1/rho bending strength.

    gamma : float
        Relativistic gamma factor of generating charged particle.

    Returns
    -------
    integ_prob : float
        Integrated probability. Will be in the range [0, 1].
    """

def bend_photon_energy_normalized_probability(E_rel: float) -> float:
    """
    Routine to return the normalized probability that a photon will be emitted in a bend with energy
    E_rel relative to the critical energy. The probability is normalized such that
      Integral[0,Infinity] dE_rel P(E_rel) = 1

    Parameters
    ----------
    E_rel : float
        Photon energy relative to the critical energy.

    Returns
    -------
    prob : float
        Normalized probability.
    """

def bend_photon_init(g_bend_x: float, g_bend_y: float, gamma: float, E_min: float | None = None, E_max: float | None = None, E_integ_prob: float | None = None, vert_angle_min: float | None = None, vert_angle_max: float | None = None, vert_angle_symmetric: bool | None = None, emit_probability: float | None = None) -> _pybmad.CoordStruct:
    """
    vert_angle_min, vert_angle_max, vert_angle_symmetric, emit_probability)

    Routine to initalize a photon for dipole bends and wigglers (but not undulators).
    The photon is initialized using the standard formulas for bending radiation.

    The energy of the photon is calculated in one of two ways:

      1) If E_integ_prob is present and non-negative, the photon energy E will be such that the integrated
          probability  [E_min, E] relative to the integrated probability in the range [E_min, E_max] is E_integ_prob.
          That is, E_integ_prob can be used to to give a set of photon energies equally spaced in terms of the
          integrated probability distribution.

      2) If E_integ_prob is not present, or is negative, the photon energy is chosen at random in
          the range [E_min, E_max].

    An E_integ_prob of zero means that the generated photon will have energy E_min.
    An E_integ_prob of one means that the generated photon will have energy E_max.

    The photon's polarization, will have unit amplitude.

    This routine assumes that the emitting charged particle is on-axis and moving in
    the forward direction. To correct for the actual charged particle postion use the routine
      absolute_photon_position

    Parameters
    ----------
    g_bend_x : float
        Bending 1/rho component in horizontal plane.

    g_bend_y : float
        Bending 1/rho component in vertical plane.

    gamma : float
        Relativistic gamma factor of generating charged particle.

    E_min : float, optional
        Minimum photon energy. Default is zero. Ignored if negative.

    E_max : float, optional
        Maximum photon energy.  Default is Infinity. Ignored if negative. If non-positive then E_max will be taken
        to be Infinity.

    E_integ_prob : float, optional
        , optional :: integrated energy probability. See above. If E_integ_prob is non-negative, it must be in the
        range [0, 1].

    vert_angle_min : float, optional
        Minimum vertical angle to emit a photon. -pi/2 is used if argument not present or if argument is less than
        -pi/2.

    vert_angle_max : float, optional
        Maximum vertical angle to emit a photon. pi/2 is used if argument not present or if argument is greater
        than pi/2.

    vert_angle_symmetric : bool, optional
        Default is False. If True, photons will be emitted in the range [-vert_angle_max, -vert_angle_min] as well
        as the range [vert_angle_min, vert_angle_max]. In this case vert_angle_min/max must be positive.

    emit_probability : float, optional
        Probability of emitting a photon in the range [E_min, E_max] or in the vertical angular range given. The
        probability is normalized so that the probability of emitting if no ranges are given is 1.

    Returns
    -------
    orbit : CoordStruct
        Initialized photon.
    """

def bend_photon_polarization_init(g_bend_x: float, g_bend_y: float, E_rel: float, gamma_phi: float) -> _pybmad.CoordStruct:
    """
    Routine to set a photon's polarization.
    The photon's polarization will be either in the plane of the bend or out of the plane and
    the magnitude will be 1.

    Parameters
    ----------
    g_bend_x : float
        Bending 1/rho component in horizontal plane.

    g_bend_y : float
        Bending 1/rho component in vertical plane.

    E_rel : float
        Relative photon energy E/E_crit.

    gamma_phi : float
        gamma * phi where gamma is the beam relativistic factor and phi is the vertical photon angle (in radians).

    Returns
    -------
    orbit : CoordStruct
        Photon coords
    """

def bend_photon_vert_angle_init(E_rel: float, gamma: float, r_in: float | None = None, invert: bool | None = None) -> float:
    """
    Routine to convert an integrated probability to a vertical angle for emitting a photon from a bend.
    The integrated probability is in the range [0,1] with 0 corresponding to a phi = -pi/2 and
    integrated probability of 1 corresponding to phi = pi/2.

    Parameters
    ----------
    E_rel : float
        Relative photon energy E/E_crit.

    gamma : float
        beam relativistic factor

    r_in : float, optional
        Integrated probability in the range [0,1]. If not present, a random number will be used.

    invert : bool, optional
        If True then take r_in as the inverse integrated probability with inverted probability = 1 - probability.
        This is useful to avoid round-off errors when for looking at the tail of the distribution where the
        integrated prob is very close to 1 and small deviations can have large effects. Default is False.

    Returns
    -------
    phi : float
        The photon vertical emission angle (in radians). Note: phi is an increasing monotonic function of r_in.
    """

class BendShift:
    """bend_shift return type"""

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def position2(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bend_shift(position1: _pybmad.FloorPositionStruct, g: float, delta_s: float, ref_tilt: float | None = None) -> BendShift:
    """
    Wrapper for Fortran routine bend_shift

    Parameters
    ----------
    position1 : FloorPositionStruct
        Position of particle in inital coordinate frame.

    g : float
        Curvature (1/rho)

    delta_s : float
        S-position of final frame relative to the initial frame.

    ref_tilt : float, optional
        ref_tilt. Default: 0

    Returns
    -------
    position2 : FloorPositionStruct
        particle coordinates relative to the final frame.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix used in the transformation
    """

def bend_vert_angle_integ_prob(vert_angle: float, E_rel: float, gamma: float) -> float:
    """
    Routine to find the integrated probability corresponding to emitting a photon
    from a bend and with relative energy E_rel in the vertical angle range [-pi/2, vert_angle/2].

    Note: vert_angle is allowed to be out of the range [-pi/2, pi/2]. In this case, integ_prob
    will be set to 0 or 1 as appropriate.

    Parameters
    ----------
    vert_angle : float
        Vertical angle.

    E_rel : float
        Relative photon energy E/E_crit.

    gamma : float
        Relativistic gamma factor of generating charged particle.

    Returns
    -------
    integ_prob : float
        Integrated probability. Will be in the range [0, 1].
    """

def bjmt1(ele: _pybmad.EleStruct, coulomb_log: float, n_part: float) -> _pybmad.IbsStruct:
    """
    This is an implementation of equations 1-9 from "Intrabeam
    scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.

    This formulation is one of the slowest methods for calculating IBS rates.

    rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
    """

def bl_via_mat(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct, mode: _pybmad.NormalModesStruct, sig_z: float) -> None:
    """
    Calculates bunch length while taking PWD effects into account.
    PWD is approximated as a defocusing rf voltage.
    """

def bl_via_vlassov(current: float, alpha: float, Energy: float, sigma_p: float, Vrf: float, omega: float, U0: float, circ: float, R: float, L: float) -> float:
    """
    This is a frontend for get_bl_from_fwhm from longitudinal_profile_mod.
    See longitudinal_profile_mod for details.  In short, this implements a model of potential well distortion
    based on the Vlassov equation which uses an effective Resistive, Inductive, and Capacitive impedance.

    Parameters
    ----------
    current : float
        Beam current in amps

    alpha : float
        Momentum compaction

    Energy : float
        beam energy

    sigma_p : float
        energy spread

    Vrf : float
        total RF voltage in Volts

    omega : float
        rf frequency in radians/s

    U0 : float
        energy loss per turn (eV)

    circ : float
        circumpherence

    R : float
        Resistive part of effective impedance

    L : float
        Inductive part of effective impedance

    Returns
    -------
    sigma_z : float
        Bunch length. FWHM/TwoRootTwoLogTwo from bunch profile
    """

class BmadParser:
    """bmad_parser return type"""

    @property
    def lat(self) -> _pybmad.LatStruct: ...

    @property
    def digested_read_ok(self) -> bool: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def parse_lat(self) -> _pybmad.LatStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bmad_parser(lat_file: str, make_mats6: bool | None = None, use_line: str | None = None) -> BmadParser:
    """
    Wrapper for Fortran routine bmad_parser

    Parameters
    ----------
    lat_file : str
        Name of the input file.

    make_mats6 : bool, optional
        Compute the 6x6 transport matrices for the Elements? Default is True. Do not set False unless you know
        what you are doing.

    use_line : str, optional
        If present and not blank, override the use statement in the lattice file and use use_line instead.

    Returns
    -------
    lat : LatStruct
        Lat structure. See bmad_struct.f90 for more details.

    digested_read_ok : bool, optional
        Set True if the digested file was successfully read. False otherwise.

    err_flag : bool, optional
        Set true if there is an error, false otherwise. Note: err_flag does *not* include errors in lat_make_mat6
        since if there is a match element, there is an error raised since the Twiss parameters have not been set
        but this is expected.

    parse_lat : LatStruct, optional
        List of elements used to construct the lattice. Useful if bmad_parser2 will be called. See bmad_parser2
        documentation.
    """

def bmad_parser2(lat_file: str, lat: _pybmad.LatStruct, orbit: _pybmad.CoordStructArray1D | None = None, make_mats6: bool | None = None, err_flag: bool | None = None, parse_lat: _pybmad.LatStruct | None = None) -> None:
    """
    Wrapper for Fortran routine bmad_parser2

    Parameters
    ----------
    lat_file : str
        Input file name.

    lat : LatStruct
        lattice with existing layout.
        This parameter is an input/output and is modified in-place.
        As an output, lat: lattice with modifications.

    orbit : 1D array of CoordStruct, optional
        closed orbit for when bmad_parser2 calls lat_make_mat6

    make_mats6 : bool, optional
        Make the 6x6 transport matrices for then Elements? Default is True.

    err_flag : bool, optional

    parse_lat : LatStruct, optional
        Used by bmad_parser to pass to bmad_parser2 a list of elements that were defined in the lattice file but
        not used. This is useful in preventing errors being generated if group/overlay elements definded by
        lat_file refer to unused slaves in parse_lat.
    """

def bmad_parser_string_attribute_set(ele: _pybmad.EleStruct, attrib_name: str, delim: str, delim_found: bool, pele: _pybmad.ParserEleStruct | None = None, str_out: str | None = None) -> None:
    """
    Wrapper for Fortran routine bmad_parser_string_attribute_set

    Parameters
    ----------
    ele : EleStruct

    attrib_name : str

    delim : str

    delim_found : bool

    pele : ParserEleStruct, optional

    str_out : str, optional
    """

def bmad_patch_parameters_to_ptc(ang: Sequence[float], exi: Sequence[Sequence[float]]) -> None:
    """
    Wrapper for Fortran routine bmad_patch_parameters_to_ptc

    Parameters
    ----------
    ang : 1D array of float (shape: 3)

    exi : 2D array of float (shape: 3,3)
    """

def bp_set_ran_status() -> None:
    """Wrapper for Fortran routine bp_set_ran_status"""

def branch_equal_branch(branch1: _pybmad.BranchStruct, branch2: _pybmad.BranchStruct) -> None:
    """
    Wrapper for Fortran routine branch_equal_branch

    Parameters
    ----------
    branch1 : BranchStruct

    branch2 : BranchStruct
    """

def branch_name(branch: _pybmad.BranchStruct) -> str:
    """
    Wrapper for Fortran routine branch_name

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch

    Returns
    -------
    name : str
        Encoded name
    """

def branch_to_ptc_m_u(branch: _pybmad.BranchStruct) -> None:
    """
    Subroutine to create a PTC layout from a Bmad lattice branch.
    Note: If lat_to_ptc_layout has already been setup, you should first do a
              call kill_ptc_layouts(lat)
    This deallocates the pointers in PTC

    Note: If not already done, before you call this routine you need to first call:
       call set_ptc (...)
    [This is normally done in bmad_parser.]

    Note: If a Bmad element is using a hard edge model (EG: RFcavity element), there
    will be three corresponding PTC fibre elements: (drift, RF. drift) for example.
    In this case, ele%ptc_fibre will be set to point to the last PTC fibre. That is the
    exit end of ele will correspond to the exit end of ele%ptc_fibre.

    Parameters
    ----------
    branch : BranchStruct
        Input branch.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Pointers to generated layouts.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Pointer to PTC fibres
    """

def bunch_equal_bunch(bunch1: _pybmad.BunchStruct, bunch2: _pybmad.BunchStruct) -> None:
    """
    Wrapper for Fortran routine bunch_equal_bunch

    Parameters
    ----------
    bunch1 : BunchStruct

    bunch2 : BunchStruct
    """

def c_to_cbar(ele: _pybmad.EleStruct) -> list[list[float]]:
    """
    Wrapper for Fortran routine c_to_cbar

    Parameters
    ----------
    ele : EleStruct
        Element with C matrix and Twiss parameters.

    Returns
    -------
    cbar_mat : 2D array of float (shape: 2,2)
        Cbar matrix.
    """

class CalcBunchParams:
    """calc_bunch_params return type"""

    @property
    def bunch_params(self) -> _pybmad.BunchParamsStruct: ...

    @property
    def error(self) -> bool: ...

    @property
    def n_mat(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def calc_bunch_params(bunch: _pybmad.BunchStruct, print_err: bool | None = None, is_time_coords: bool | None = None, ele: _pybmad.EleStruct | None = None) -> CalcBunchParams:
    """
    Finds all bunch parameters defined in bunch_params_struct, both normal-mode
    and projected. Projected parameters are found purely from the geometrical
    distribution of the beam. Normal-Mode parameters are found using the method
    developed in:
      "Alternate approach to general coupled linear optics"
       A. Wolski, PRST AB 9, 024001 (2006)

    Note: If less than two particle remain then the various parameters will be
    set to zero.

    Parameters
    ----------
    bunch : BunchStruct
        Bunch_struct

    print_err : bool, optional
        If present and False then suppress "no eigen-system found" messages.

    is_time_coords : bool, optional
        Are particle coords using time coords. Default is False.

    ele : EleStruct, optional
        Element being tracked through. Must be present if is_time_coords = True.

    Returns
    -------
    bunch_params : BunchParamsStruct

    error : bool
        Set True if there is an error.

    n_mat : 2D array of float (shape: 6,6), optional
        N matrix defined in Wolski Eq 44 and used to convert from action-angle coords to lab coords (Wolski Eq
        51.).
    """

def calc_bunch_params_slice(bunch: _pybmad.BunchStruct, bunch_params: _pybmad.BunchParamsStruct, plane: int, slice_center: float, slice_spread: float, print_err: bool | None = None, is_time_coords: bool | None = None, ele: _pybmad.EleStruct | None = None) -> bool:
    """
    Finds bunch parameters for a slice of the beam.

    Parameters
    ----------
    bunch : BunchStruct
        bunch_struct

    plane : int
        plane to slice through (x$, px$, & etc...)

    slice_center : float
        Center to take slice about

    slice_spread : float
        +/- spread in slice about center.

    print_err : bool, optional
        If present and False then suppress "no eigen-system found" messages.

    is_time_coords : bool, optional
        Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
        to s-coords before bunch_params are computed.

    ele : EleStruct, optional
        Element being tracked through. Must be present if is_time_coords = True.

    Returns
    -------
    err : bool
        Set True if there is an error.
    """

def calc_bunch_params_z_slice(bunch: _pybmad.BunchStruct, bunch_params: _pybmad.BunchParamsStruct, slice_bounds: Sequence[float], print_err: bool | None = None, is_time_coords: bool | None = None, ele: _pybmad.EleStruct | None = None) -> bool:
    """
    Finds bunch parameters for a slice of the beam.

    The slice is specified in terms of percentage of particles ordered by z-position.
    For example, slice_bounds = [0.0, 0.5] specifies the trailing half of the bunch

    Parameters
    ----------
    bunch : BunchStruct
        bunch_struct

    slice_bounds : 1D array of float (shape: 2)
        Slice bounds in percentage of particles ordered by z-position. 0.0 is the back of the bunch and 1.0 is the
        front of the bunch.

    print_err : bool, optional
        If present and False then suppress "no eigen-system found" messages.

    is_time_coords : bool, optional
        Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
        to s-coords before bunch_params are computed.

    ele : EleStruct, optional
        Element being tracked through. Must be present if is_time_coords = True.

    Returns
    -------
    err : bool
        Set True if there is an error.
    """

def calc_bunch_sigma_matrix_etc(particle: _pybmad.CoordStructArray1D, charge: _pybmad.RealArray1D, is_time_coords: bool | None = None, ele: _pybmad.EleStruct | None = None) -> _pybmad.BunchParamsStruct:
    """
    Routine to find the sigma matrix elements of a particle distribution.

    Parameters
    ----------
    particle : 1D array of CoordStruct
        Array of particles.

    charge : 1D array of float
        Particle charge or photon intensity.

    Returns
    -------
    bunch_params : BunchParamsStruct
        Bunch parameters. .sigma(6,6) .centroid.vec(6) .centroid.p0c .rel_max(6) .rel_min(6)
    """

class CalcEmittancesAndTwissFromSigmaMatrix:
    """calc_emittances_and_twiss_from_sigma_matrix return type"""

    @property
    def bunch_params(self) -> _pybmad.BunchParamsStruct: ...

    @property
    def error(self) -> bool: ...

    @property
    def n_mat(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def calc_emittances_and_twiss_from_sigma_matrix(sigma_mat: Sequence[Sequence[float]], print_err: bool | None = None) -> CalcEmittancesAndTwissFromSigmaMatrix:
    """
    Routine to calc emittances and Twiss function from a beam sigma matrix.
    See: Andy Wolski "Alternative approach to general coupled linear optics".

    Parameters
    ----------
    sigma_mat : 2D array of float (shape: 6,6)
        Sigma matrix.

    print_err : bool, optional
        If present and False then suppress "no eigen-system found" messages.

    Returns
    -------
    bunch_params : BunchParamsStruct
        Holds Twiss and emittance info.

    error : bool
        Set True if there is an error. Can happen if the emittance of a mode is zero.

    n_mat : 2D array of float (shape: 6,6), optional
        N matrix defined in Wolski Eq 44 and used to convert from action-angle coords to lab coords (Wolski Eq
        51.).
    """

class CalcNextFringeEdge:
    """calc_next_fringe_edge return type"""

    @property
    def s_edge_body(self) -> float: ...

    @property
    def fringe_info(self) -> _pybmad.FringeFieldInfoStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def calc_next_fringe_edge(track_ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, init_needed: bool | None = None, time_tracking: bool | None = None) -> CalcNextFringeEdge:
    """
    Wrapper for Fortran routine calc_next_fringe_edge

    Parameters
    ----------
    track_ele : EleStruct
        Element being tracked through.

    orbit : CoordStruct
        Particle position

    init_needed : bool, optional
        If present and True then initialize.

    time_tracking : bool, optional
        If present and True then this routine is being called by the time Runge-Kutta tracker. Default is False.

    Returns
    -------
    s_edge_body : float
        S position of next hard edge in track_ele body frame. If there are no more hard edges then s_edge_body
        will be set to ele.value(l$) if orbit.direction*orbit.time_dir*ele.orientation = 1, and set to 0
        otherwise.

    fringe_info : FringeFieldInfoStruct
        Information on the next fringe to track through.
    """

def calc_spin_params(bunch: _pybmad.BunchStruct) -> _pybmad.BunchParamsStruct:
    """
    Rotine to calculate spin averages

    Parameters
    ----------
    bunch : BunchStruct
        Bunch of spins

    Returns
    -------
    bunch_params : BunchParamsStruct
        Structure holding average
    """

def calc_super_slave_key(lord1: _pybmad.EleStruct, lord2: _pybmad.EleStruct, create_jumbo_slave: bool | None = None) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine calc_super_slave_key

    Parameters
    ----------
    lord1 : EleStruct
        First slave. .key .sub_key

    lord2 : EleStruct
        Second slave. .key .sub_key

    create_jumbo_slave : bool, optional
        If True then slave.key will be set to em_field. Default is False.

    Returns
    -------
    slave : EleStruct
        Super_slave element.
    """

class CalcWallRadius:
    """calc_wall_radius return type"""

    @property
    def r_wall(self) -> float: ...

    @property
    def dr_dtheta(self) -> float: ...

    @property
    def ix_vertex(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def calc_wall_radius(v: _pybmad.Wall3DVertexStructArray1D, cos_ang: float, sin_ang: float) -> CalcWallRadius:
    """
    Routine to calculate the wall radius at a given angle for a given cross-section
    Additionally, the transverse directional derivative is calculated.

    Module needed:
      use wall3d_mod

    Parameters
    ----------
    v : 1D array of Wall3dVertexStruct
        Array of vertices that make up the cross-section.

    cos_ang : float
        cosine of the transverse photon position.

    sin_ang : float
        sine of the transverse photon position.

    Returns
    -------
    r_wall : float
        Wall radius at given angle.

    dr_dtheta : float
        derivative of r_wall.

    ix_vertex : int, optional
        Wall at given angle is between v(ix_vertex-1) and v(ix_vertex). If ix_vertex = 1 then Wall at given angle
        is between v(N) and v(1) where N = size(v).
    """

def calc_wiggler_g_params(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, s_rel: float, orb: _pybmad.CoordStruct, pt: _pybmad.RadIntTrackPointStruct, info: _pybmad.RadIntInfoStruct | None = None) -> None:
    """
    Wrapper for Fortran routine calc_wiggler_g_params

    Parameters
    ----------
    ele : EleStruct

    param : LatParamStruct

    s_rel : float

    orb : CoordStruct

    pt : RadIntTrackPointStruct

    info : RadIntInfoStruct, optional
    """

def calc_z_tune(branch: _pybmad.BranchStruct) -> None:
    """
    Wrapper for Fortran routine calc_z_tune

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch
        This parameter is an input/output and is modified in-place.
        As an output, branch: Synchrotron tune (radians). If unstable tune = 0.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Is the mode stable? If no rf then tune is zero but is stable.
        This parameter is an input/output and is modified in-place.
        As an output, branch: 6x6 1-turn matrix.
    """

def canonical_to_angle_coords(orbit: _pybmad.CoordStruct, coord_type: str | None = None) -> None:
    r"""
    Wrapper for Fortran routine canonical_to_angle_coords

    Parameters
    ----------
    orbit : CoordStruct
        Orbit in canonical coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit in angular coordinates.

    coord_type : str, optional
        Angular coordinates type \'\' (default): (x, x' = dx/ds, y, y' = dy/ds, z, pz) 'ZGOUBI':     (x, x' = dx/ds,
        y, y' = dy/ds, dt = -z / (beta * c), pz)
    """

def capillary_photon_hit_spot_calc(photon: _pybmad.PhotonTrackStruct, ele: _pybmad.EleStruct) -> None:
    """
    Routine to interpolate to where the photon has hit the capillary.

    Parameters
    ----------
    photon : PhotonTrackStruct
        Input coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, photon: Photon at capillary wall

    ele : EleStruct
        Capillary element
    """

def capillary_propagate_photon_a_step(photon: _pybmad.PhotonTrackStruct, ele: _pybmad.EleStruct, dlen: float) -> bool:
    """
    Routine to track a photon a step of a given length

    Parameters
    ----------
    photon : PhotonTrackStruct
        Input coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, photon: Output coordinates.

    ele : EleStruct
        Capillary element

    dlen : float
        Length to propagate a photon. The actual propagation length may be less if stop_at_boundary is True.

    Returns
    -------
    stop_at_boundary : bool
        If True then stop at cross-section boundary.
    """

def capillary_reflect_photon(photon: _pybmad.PhotonTrackStruct, ele: _pybmad.EleStruct) -> None:
    """
    Routine to reflect a photon from the capillary wall.

    Parameters
    ----------
    photon : PhotonTrackStruct
        Input coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, photon: Output coordinates.

    ele : EleStruct
        Capillary element
    """

def capillary_track_photon_to_wall(photon: _pybmad.PhotonTrackStruct, ele: _pybmad.EleStruct) -> None:
    """
    Routine to track through a capillary.

    Parameters
    ----------
    photon : PhotonTrackStruct
        Input coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, photon: Output coordinates.

    ele : EleStruct
        Capillary element
    """

def cbar_to_c(cbar_mat: Sequence[Sequence[float]], a: _pybmad.TwissStruct, b: _pybmad.TwissStruct) -> list[list[float]]:
    """
    Wrapper for Fortran routine cbar_to_c

    Parameters
    ----------
    cbar_mat : 2D array of float (shape: 2,2)
        Cbar matrix.

    a : TwissStruct
        a-mode Twiss parameters

    b : TwissStruct
        b-mode Twiss parameters

    Returns
    -------
    c_mat : 2D array of float (shape: 2,2)
        C matrix.
    """

def check_aperture_limit(orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, particle_at: int, param: _pybmad.LatParamStruct, old_orb: _pybmad.CoordStruct | None = None, check_momentum: bool | None = None) -> None:
    """
    Wrapper for Fortran routine check_aperture_limit

    Parameters
    ----------
    orb : CoordStruct
        coordinates of a particle.

    ele : EleStruct
        Element holding the aperture

    particle_at : int
        first_track_edge$, second_track_edge$, surface$, in_between$

    param : LatParamStruct
        Lattice global parameter structure.

    old_orb : CoordStruct, optional
        Old coordinates at last check. Needed if ele.aperture_at = wall_transition$. If not present then wall
        transitions will be ignored.

    check_momentum : bool, optional
        If present and false then checking of p_x and p_y will be disabled.
    """

def check_controller_controls(ele_key: int, contrl: _pybmad.ControlStructArray1D, name: str) -> bool:
    """
    Wrapper for Fortran routine check_controller_controls

    Parameters
    ----------
    ele_key : int
        Element type. overlay$, etc.

    contrl : 1D array of ControlStruct
        control info. 1 element for each slave.

    name : str
        Lord name. Used for error reporting.

    Returns
    -------
    err : bool
        Set true if there is a problem. False otherwise.
    """

def check_for_superimpose_problem(branch: _pybmad.BranchStruct, super_ele: _pybmad.EleStruct, err_flag: bool, wrap: bool, ref_ele: _pybmad.EleStruct | None = None) -> None:
    """
    Subroutine to check if there is a problem superimposing an element when there is multipass.
    In particular will check that:
      1) If the ref_ele is part of a multipass region then super_ele must be superimposed
         within the region.
    Or:
      2) If the ref_ele is not part of a multipass region then super_ele must also not
         be part of a multipass region.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

class CheckIfSInBounds:
    """check_if_s_in_bounds return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def translated_s(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def check_if_s_in_bounds(branch: _pybmad.BranchStruct, s: float, print_err: bool | None = None) -> CheckIfSInBounds:
    """
    Wrapper for Fortran routine check_if_s_in_bounds

    Parameters
    ----------
    branch : BranchStruct
        Branch

    s : float
        longitudinal position in the given branch.

    print_err : bool, optional
        Print error message if there is an error? Default is True.

    Returns
    -------
    err_flag : bool
        Set True if s position is out-of-bounds. False otherwise.

    translated_s : float, optional
        position translated to the range [0, branch_length]
    """

class ChooseQuadsForSetTune:
    """choose_quads_for_set_tune return type"""

    @property
    def dk1(self) -> _pybmad.RealAlloc1D: ...

    @property
    def eles(self) -> _pybmad.ElePointerStructAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def choose_quads_for_set_tune(branch: _pybmad.BranchStruct, mask: str | None = None) -> ChooseQuadsForSetTune:
    """
    Wrapper for Fortran routine choose_quads_for_set_tune

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch.

    mask : str, optional
        If present, assign weight of zero for all quads that do not match. That is, no variation for matching
        quads.

    Returns
    -------
    dk1 : 1D array of float
        Weights for the quadrupoles. All values will be +1 or -1.

    eles : 1D array of ElePointerStruct
        eles(i).ele points to element with dk1(i) weight.

    err_flag : bool, optional
        Set True if there is not one quad with positive dk1 and one quad with negative dk1.
    """

class ChromCalc:
    """chrom_calc return type"""

    @property
    def chrom_a(self) -> float: ...

    @property
    def chrom_b(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def low_E_lat(self) -> _pybmad.LatStruct: ...

    @property
    def high_E_lat(self) -> _pybmad.LatStruct: ...

    @property
    def low_E_orb(self) -> _pybmad.CoordStructAlloc1D: ...

    @property
    def high_E_orb(self) -> _pybmad.CoordStructAlloc1D: ...

    @property
    def delta_e(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def chrom_calc(lat: _pybmad.LatStruct, delta_e: float, pz: float | None = None, ix_branch: int | None = None, orb0: _pybmad.CoordStruct | None = None) -> ChromCalc:
    """
    Wrapper for Fortran routine chrom_calc

    Parameters
    ----------
    lat : LatStruct
        Lat

    delta_e : float
        +/- Delta energy used for the calculation. Notice that the energy difference between high and low is 2 *
        delta_e. If 0 then default of 1.0d-4 is used.
        This parameter is an input/output and is modified in-place.
        As an output, delta_e: Value used in computation. Set to 1.0d-4 if on input delta_e =< 0.

    pz : float, optional
        reference momentum about which to calculate. Default is 0.

    ix_branch : int, optional
        Index of the lattice branch to use. Default is 0.

    orb0 : CoordStruct, optional
        On-energy orbit at start (fixer point). Default is the branch.particle_start. Only needed if lattice
        branch has an open geometry.

    Returns
    -------
    delta_e : float
        +/- Delta energy used for the calculation. Notice that the energy difference between high and low is 2 *
        delta_e. If 0 then default of 1.0d-4 is used.
        This parameter is an input/output and is modified in-place.
        As an output, delta_e: Value used in computation. Set to 1.0d-4 if on input delta_e =< 0.

    chrom_a : float
        a-mode chromaticity.

    chrom_b : float
        b-mode chromaticity.

    err_flag : bool, optional
        Set true if there is an error. False otherwise.

    low_E_lat : LatStruct, optional
        Lattice with RF off and matrices computed at E_lat +pz - delta_e

    high_E_lat : LatStruct, optional
        Lattice with RF off and matrices computed at E_lat +pz + delta_e

    low_E_orb : 1D array of CoordStruct, optional
        Orbit computed at E_lat + pz - delta_e.

    high_E_orb : 1D array of CoordStruct, optional
        Orbit computed at E_lat + pz + delta_e.
    """

class ChromTune:
    """chrom_tune return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def delta_e(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def chrom_tune(lat: _pybmad.LatStruct, delta_e: float, target_x: float, target_y: float, err_tol: float) -> ChromTune:
    """
    Wrapper for Fortran routine chrom_tune

    Parameters
    ----------
    lat : LatStruct
        Lat to use,
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lat with sextupole set

    delta_e : float
        Delta energy used for the calculation. If 0 then default of 1.0d-4 is used.
        This parameter is an input/output and is modified in-place.
        As an output, delta_e: Set to 1.0d-4 if on input DELTA_E =< 0.

    target_x : float
        Target X Chromaticity

    target_y : float
        Target Y Chromaticity

    err_tol : float
        Max allowable Error: Error = | X_Target - X_Actual | + | Y_Target -Y_Actual | A good number is: err_tol =
        0.05_rp

    Returns
    -------
    delta_e : float
        Delta energy used for the calculation. If 0 then default of 1.0d-4 is used.
        This parameter is an input/output and is modified in-place.
        As an output, delta_e: Set to 1.0d-4 if on input DELTA_E =< 0.

    err_flag : bool
        .false. if match successful, .true. if failed Fails if takes longer than 100 iterations. If it fails the
        sextupoles are set to the last value calculated. Note: This subroutine assumes the Twiss parameters have
        been computed.
    """

def cimp1(ele: _pybmad.EleStruct, coulomb_log: float, n_part: float) -> _pybmad.IbsStruct:
    """
    This is an implementation of equations 34,38-40 from "Intrabeam
    scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
    It is a modified version of the Piwinski IBS formulation.
    The integral (34) is handled with a piecewise interpolation generated
    in mathematica.  The interpolation is accurate beyond 1% through it's
    effective range (.0001 - 3000).

    This is the quickest of the three IBS formuations in this module.

    rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
    """

def classical_radius(species: int) -> float:
    """
    Wrapper for Fortran routine classical_radius

    Parameters
    ----------
    species : int
        Species of particle.

    Returns
    -------
    radius : float
        Classical radius.
    """

def clear_lat_1turn_mats() -> _pybmad.LatStruct:
    """
    Wrapper for Fortran routine clear_lat_1turn_mats

    Returns
    -------
    lat : LatStruct
        Lat with 1-turn matrices cleared.
    """

def clear_taylor_maps_from_elements(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine clear_taylor_maps_from_elements

    Parameters
    ----------
    lat : LatStruct
        Lattice
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with all maps cleared
    """

def closed_orbit_calc(lat: _pybmad.LatStruct, closed_orb: _pybmad.CoordStructAlloc1D, i_dim: int | None = None, direction: int | None = None, ix_branch: int | None = None, print_err: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine closed_orbit_calc

    Parameters
    ----------
    lat : LatStruct
        Lat to track through.

    closed_orb : 1D array of CoordStruct
        closed_orb(nt) is the initial guess where nt = 0 for direction = 1 and nt = lat.n_ele_track for direction
        = -1. Additionally, if i_dim = 4, then closed_orb(nt).vec(6) is used as the energy around which the closed
        orbit is calculated.
        This parameter is an input/output and is modified in-place.
        As an output, closed_orb: Closed orbit. closed_orb(i)

    i_dim : int, optional
        Phase space dimensions to use: = 4  Transverse closed orbit at constant energy (RF off). (dE/E =
        closed_orb(0).vec(6)) = 5 Transverse closed orbit at constant energy (RF off) with the energy adjusted so
        that vec(5) is the same at the beginning and at the end. = 6 True closed orbit. Default: 4 if RF is off, 6
        if RF is on.

    direction : int, optional
        Direction of tracking.

    ix_branch : int, optional
        Lattice branch to find the closed orbit of. Default is 0 (main branch).

    print_err : bool, optional
        Print error message if calc does not converge? Default is True. Note: Condition messages like no RF
        voltage with i_dim = 6 will always be printed.

    Returns
    -------
    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

class ClosedOrbitFromTracking:
    """closed_orbit_from_tracking return type"""

    @property
    def closed_orb(self) -> _pybmad.CoordStructAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def closed_orbit_from_tracking(lat: _pybmad.LatStruct, i_dim: int, eps_rel: _pybmad.RealArray1D | None = None, eps_abs: _pybmad.RealArray1D | None = None, init_guess: _pybmad.CoordStruct | None = None) -> ClosedOrbitFromTracking:
    """
    Wrapper for Fortran routine closed_orbit_from_tracking

    Parameters
    ----------
    lat : LatStruct
        Lat to track through.

    i_dim : int
        = 2,4  Transverse closed orbit at constant energy. = 6    Full closed orbit using the entire transfer 6x6
        matrix.

    eps_rel : 1D array of float, optional
        Relative allowed error. Default is bmad_com.rel_tol_tracking

    eps_abs : 1D array of float, optional
        Absolute allowed error. Default is bmad_com.abs_tol_tracking

    init_guess : CoordStruct, optional
        Starting guess for the closed orbit at the start of the lattice. Set init_guess.vec(6) to the appropriate
        value of pz when calculating off-energy orbits. If not present then the origin will be used.

    Returns
    -------
    closed_orb : 1D array of CoordStruct
        closed orbit. This routine will allocate this array for you.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

def cmplx_re_str(cmp: complex) -> str:
    """
    Wrapper for Fortran routine cmplx_re_str

    Parameters
    ----------
    cmp : complex

    Returns
    -------
    str_out : str
    """

def combine_consecutive_elements(lat: _pybmad.LatStruct) -> bool:
    """
    Wrapper for Fortran routine combine_consecutive_elements

    Parameters
    ----------
    lat : LatStruct
        Lattice.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with elements combined.

    Returns
    -------
    error : bool
        Set True if there is an error. False otherwise.
    """

def complex_taylor_clean(complex_taylor: _pybmad.ComplexTaylorStruct) -> None:
    """
    Wrapper for Fortran routine complex_taylor_clean

    Parameters
    ----------
    complex_taylor : ComplexTaylorStruct
    """

@overload
def complex_taylor_coef(complex_taylor: _pybmad.ComplexTaylorStruct, exp: _pybmad.IntArray1D) -> complex: ...

@overload
def complex_taylor_coef(complex_taylor: _pybmad.ComplexTaylorStruct, i1: int | None = None, i2: int | None = None, i3: int | None = None, i4: int | None = None, i5: int | None = None, i6: int | None = None, i7: int | None = None, i8: int | None = None, i9: int | None = None) -> complex:
    """
    Function to return the coefficient for a particular complex_taylor term
    from a complex_taylor Series.

    Note: complex_taylor_coef is overloaded by:
      complex_taylor_coef1 (complex_taylor, exp)
      complex_taylor_coef2 (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
    Using the complex_taylor_coef2 form limits obtaining coefficients to 9th order
    or less. Also: complex_taylor_coef2 does not check that all i1, ..., i9 are between
    1 and 6.

    For example: To get the 2nd order term corresponding to
      y(out) = Coef * p_z(in)^2
    [This is somtimes refered to as the T_366 term]
    The call would be:
      type (complex_taylor_struct) complex_taylor(6)      ! complex_taylor Map
      ...
      coef = complex_taylor_coef (complex_taylor(3), 6, 6)  ! 1st possibility or ...
      coef = complex_taylor_coef (complex_taylor(3), [0, 0, 0, 0, 0, 2 ])

    Input (complex_taylor_coef1):
      complex_taylor -- complex_taylor_struct: complex_taylor series.
      exp(6)      -- Integer: Array of exponent indices.

    Input (complex_taylor_coef2):
      complex_taylor -- complex_taylor_struct: complex_taylor series.
      i1, ..., i9 -- Integer, optional: indexes (each between 1 and 6).
    """

def complex_taylor_equal_complex_taylor(complex_taylor1: _pybmad.ComplexTaylorStruct, complex_taylor2: _pybmad.ComplexTaylorStruct) -> None:
    """
    Wrapper for Fortran routine complex_taylor_equal_complex_taylor

    Parameters
    ----------
    complex_taylor1 : ComplexTaylorStruct

    complex_taylor2 : ComplexTaylorStruct
    """

def complex_taylor_exponent_index(expn: Sequence[int]) -> int:
    """
    Function to associate a unique number with a complex_taylor exponent.

    The number associated with a complex_taylor_term that is used for the sort is:
        number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
    where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.

    Parameters
    ----------
    expn : 1D array of int (shape: 6)
        complex_taylor exponent

    Returns
    -------
    index : int
        Sorted complex_taylor series.
    """

def complex_taylor_make_unit(complex_taylor: _pybmad.ComplexTaylorStructArray1D) -> None:
    """
    Subroutine to make the unit complex_taylor map:
          r(out) = Map * r(in) = r(in)
    """

class ComplexTaylorToMat6:
    """complex_taylor_to_mat6 return type"""

    @property
    def vec0(self) -> list[complex]: ...

    @property
    def mat6(self) -> list[list[complex]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def complex_taylor_to_mat6(a_complex_taylor: _pybmad.ComplexTaylorStructArray1D, r_in: _pybmad.ComplexArray1D, r_out: _pybmad.ComplexArray1D | None = None) -> ComplexTaylorToMat6:
    """
    Subroutine to calculate, from a complex_taylor map and about some trajectory:
      The 1st order (Jacobian) transfer matrix.

    Parameters
    ----------
    a_complex_taylor : 1D array of ComplexTaylorStruct (shape: 6)
        complex_taylor map.

    r_in : 1D array of complex
        Coordinates at the input.

    r_out : 1D array of complex, optional
        Coordinates at output.

    Returns
    -------
    vec0 : 1D array of complex (shape: 6)
        0th order tranfsfer map

    mat6 : 2D array of complex (shape: 6,6)
        1st order transfer map (6x6 matrix).
    """

def complex_taylors_equal_complex_taylors(complex_taylor1: _pybmad.ComplexTaylorStructArray1D, complex_taylor2: _pybmad.ComplexTaylorStructArray1D) -> None:
    """
    Wrapper for Fortran routine complex_taylors_equal_complex_taylors

    Parameters
    ----------
    complex_taylor1 : 1D array of ComplexTaylorStruct

    complex_taylor2 : 1D array of ComplexTaylorStruct
    """

def compute_slave_coupler(slave: _pybmad.EleStruct) -> None:
    """This routine is not meant for general use."""

def compute_super_lord_s(ref_ele: _pybmad.EleStruct, super_ele: _pybmad.EleStruct, pele: _pybmad.ParserEleStruct, ix_insert: int) -> None:
    """
    Wrapper for Fortran routine compute_super_lord_s

    Parameters
    ----------
    ref_ele : EleStruct

    super_ele : EleStruct

    pele : ParserEleStruct

    ix_insert : int
    """

def concat_ele_taylor(orb_taylor: _pybmad.TaylorStructArray1D, ele: _pybmad.EleStruct, spin_taylor: _pybmad.TaylorStructArray1D | None = None) -> bool:
    """
    Routine to concatinate an orbital taylor map and, optionally if present and
    bmad_com%spin_tracking_on = T, a spin taylor map.

    Transform:
      orb_taylor[x] -> ele_taylor(orb_taylor[x])
    If ele%taylor_map_includes_offsets = True:  ele_taylor == ele%taylor
    If ele%taylor_map_includes_offsets = False: ele_taylor == ele%taylor + offset corrections.

    Also see: concat_taylor

    Parameters
    ----------
    orb_taylor : 1D array of TaylorStruct
        Orbital Taylor map.
        This parameter is an input/output and is modified in-place.
        As an output, orb_taylor: Concatinated orbital map

    ele : EleStruct
        Element containing a Taylor map.

    spin_taylor : 1D array of TaylorStruct, optional
        Spin map to propagate
        This parameter is an input/output and is modified in-place.
        As an output, spin_taylor: Concatinated spin map.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def concat_taylor(taylor1: _pybmad.TaylorStructArray1D, taylor2: _pybmad.TaylorStructArray1D, taylor3: _pybmad.TaylorStructArray1D) -> None:
    """
    Subroutine to concatinate two taylor maps:
      taylor3[x] = taylor2(taylor1[x])

    Note: In general, if taylor2 is a component of an ele_struct, use
    concat_ele_taylor instead.

    Parameters
    ----------
    taylor1 : 1D array of TaylorStruct
        Taylor map.

    taylor2 : 1D array of TaylorStruct
        Taylor map.

    taylor3 : 1D array of TaylorStruct
        Concatinated map
    """

class ConcatTransferMat:
    """concat_transfer_mat return type"""

    @property
    def mat_out(self) -> list[list[float]]: ...

    @property
    def vec_out(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def concat_transfer_mat(mat_1: Sequence[Sequence[float]], vec_1: Sequence[float], mat_0: Sequence[Sequence[float]], vec_0: Sequence[float]) -> ConcatTransferMat:
    """
    Routine to concatinate two linear maps:
      mat_out = matmul(mat_1, mat_0)
      vec_out = matmul(mat_1, vec_0) + vec_1

    Parameters
    ----------
    mat_1 : 2D array of float (shape: 6,6)
        Map from s1 to s2

    vec_1 : 1D array of float (shape: 6)
        Map from s1 to s2

    mat_0 : 2D array of float (shape: 6,6)
        Map from s0 to s1

    vec_0 : 1D array of float (shape: 6)
        Map from s0 to s1

    Returns
    -------
    mat_out : 2D array of float (shape: 6,6)
        Map from s0 to s2

    vec_out : 1D array of float (shape: 6)
        Map from s0 to s2
    """

def control_bookkeeper(lat: _pybmad.LatStruct, ele: _pybmad.EleStruct | None = None, err_flag: bool | None = None) -> None:
    """
    Wrapper for Fortran routine control_bookkeeper

    Parameters
    ----------
    lat : LatStruct
        lattice to be used

    ele : EleStruct, optional
        Element whose attribute values have been changed. If not present bookkeeping will be done for all
        elements.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

def convert_bend_exact_multipole(g: float, out_type: int, an: Sequence[float], bn: Sequence[float]) -> None:
    """
    Wrapper for Fortran routine convert_bend_exact_multipole

    Parameters
    ----------
    g : float
        1/rho bending strength.

    out_type : int
        Output type: horizontally_pure$ or vertically_pure$.

    an : 1D array of float (shape: 0:n_pole_maxx)
        Skew multipoles.
        This parameter is an input/output and is modified in-place.
        As an output, an: Converted skew multipoles.

    bn : 1D array of float (shape: 0:n_pole_maxx)
        Non-skew multipoles.
        This parameter is an input/output and is modified in-place.
        As an output, bn: Converted Non-skew multipoles.
    """

class ConvertCoords:
    """convert_coords return type"""

    @property
    def out_type_str(self) -> str: ...

    @property
    def coord_out(self) -> _pybmad.CoordStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def convert_coords(in_type_str: str, coord_in: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> ConvertCoords:
    """
    Wrapper for Fortran routine convert_coords

    Parameters
    ----------
    in_type_str : str
        type of the input coords.

    coord_in : CoordStruct
        Input coordinates.

    ele : EleStruct
        Provides the Twiss parameters.

    Returns
    -------
    out_type_str : str
        type of the output coords.

    coord_out : CoordStruct
        Output coordinates.

    err_flag : bool, optional
        Set True if there is an error. False otherwise. in_type_str and out_type_str can be: 'LAB'
        {x, x', y, y', z, z'} 'MODE'               {a, a', b, b', z, z'} 'NORMALIZED'         {a_bar, a'_bar,
        b_bar, b'_bar, z_bar, z'_bar} 'ACTION-ANGLE'       {j_a, phi_a, j_b, phi_b, j_z,  phi_z} x_vec = V_mat *
        (a_vec + eta_vec * z') a_bar  =  sqrt(2*j_a) * cos(phi_a) a'_bar = -sqrt(2*j_a) * sin(phi_a)
    """

def convert_field_ele_to_lab(ele: _pybmad.EleStruct, s_here: float, forward_transform: bool, calc_dfield: bool | None = None, calc_potential: bool | None = None) -> _pybmad.EmFieldStruct:
    """
    Convert fields: ele to lab coords

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    s_here : float
        real(rp) S-position.

    forward_transform : bool
        Transform foward (to lab) or reverse.

    calc_dfield : bool, optional
        If present and True then calculate the field derivatives.

    calc_potential : bool, optional
        Calc electric and magnetic potentials? Default is false. This is experimental and only implemented for
        wigglers at present.

    Returns
    -------
    field : EmFieldStruct
        EM field.
    """

def convert_local_cartesian_to_local_curvilinear(x: float, z: float, g: float, xout: float, sout: float) -> None:
    """
    Wrapper for Fortran routine convert_local_cartesian_to_local_curvilinear

    Parameters
    ----------
    x : float

    z : float

    g : float

    xout : float

    sout : float
    """

def convert_local_curvilinear_to_local_cartesian(x: float, s: float, g: float, xout: float, zout: float) -> None:
    """
    Wrapper for Fortran routine convert_local_curvilinear_to_local_cartesian

    Parameters
    ----------
    x : float

    s : float

    g : float

    xout : float

    zout : float
    """

def convert_particle_coordinates_s_to_t(particle: _pybmad.CoordStruct, s_body: float, orientation: int) -> None:
    """
    Wrapper for Fortran routine convert_particle_coordinates_s_to_t

    Parameters
    ----------
    particle : CoordStruct
        Particle with .vec(:) in s-coords.

    s_body : float
        s-position in element body coords.

    orientation : int
        ele.orientation for vec(6).
    """

def convert_particle_coordinates_t_to_s(particle: _pybmad.CoordStruct, ele: _pybmad.EleStruct, use_downstream_p0c: bool | None = None) -> float:
    """
    Wrapper for Fortran routine convert_particle_coordinates_t_to_s

    Parameters
    ----------
    particle : CoordStruct
        Particle with .vec(:) in t-coords.

    ele : EleStruct
        Element particle is going through.

    use_downstream_p0c : bool, optional
        If True (the default), use ele.value(p0c$) as the reference momentum. If False, use ele.value(p0c_start$)
        as the reference.

    Returns
    -------
    s_body : float, optional
        s-position in element body coords.
    """

class ConvertPcTo:
    """convert_pc_to return type"""

    @property
    def E_tot(self) -> float: ...

    @property
    def gamma(self) -> float: ...

    @property
    def kinetic(self) -> float: ...

    @property
    def beta(self) -> float: ...

    @property
    def brho(self) -> float: ...

    @property
    def beta1(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def convert_pc_to(pc: float, particle: int) -> ConvertPcTo:
    """
    Wrapper for Fortran routine convert_pc_to

    Parameters
    ----------
    pc : float
        Particle momentum

    particle : int
        Type of particle. positron$, etc.

    Returns
    -------
    E_tot : float, optional
        Total energy of the particle.

    gamma : float, optional
        Gamma factor.

    kinetic : float, optional
        Kinetic energy

    beta : float, optional
        velocity / c_light

    brho : float, optional
        Nominal B_field*rho_bend

    beta1 : float, optional
        1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.

    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

class ConvertTotalEnergyTo:
    """convert_total_energy_to return type"""

    @property
    def gamma(self) -> float: ...

    @property
    def kinetic(self) -> float: ...

    @property
    def beta(self) -> float: ...

    @property
    def pc(self) -> float: ...

    @property
    def brho(self) -> float: ...

    @property
    def beta1(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def convert_total_energy_to(E_tot: float, particle: int, print_err: bool | None = None) -> ConvertTotalEnergyTo:
    """
    Wrapper for Fortran routine convert_total_energy_to

    Parameters
    ----------
    E_tot : float
        Total energy of the particle.

    particle : int
        Type of particle. positron$, etc.

    print_err : bool, optional
        Print error message if E_tot < particle mass? Default is True.

    Returns
    -------
    gamma : float, optional
        Gamma factor. Set to -1 for photons.

    kinetic : float, optional
        Kinetic energy

    beta : float, optional
        velocity / c_light

    pc : float, optional
        Particle momentum

    brho : float, optional
        Nominal B_field*rho_bend

    beta1 : float, optional
        1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.

    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

def converter_distribution_parser(ele: _pybmad.EleStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine converter_distribution_parser

    Parameters
    ----------
    ele : EleStruct

    delim : str

    delim_found : bool

    err_flag : bool
    """

def coord_equal_coord(coord2: _pybmad.CoordStruct) -> _pybmad.CoordStruct:
    """
    Subroutine that is used to set one coord equal to another.

     Note: This subroutine is called by the overloaded equal sign:
    		coord1 = coord2

    Parameters
    ----------
    coord2 : CoordStruct
        Input coord.

    Returns
    -------
    coord1 : CoordStruct
        Output coord.
    """

def coord_state_name(coord_state: int, one_word: bool | None = None) -> str:
    """
    Routine to return the string representation of a coord%state state.

    Parameters
    ----------
    coord_state : int
        coord.state value

    Returns
    -------
    state_str : str
        String representation.
    """

class CoordsBodyToLocal:
    """coords_body_to_local return type"""

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def local_position(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_body_to_local(body_position: _pybmad.FloorPositionStruct, ele: _pybmad.EleStruct, calculate_angles: bool | None = None) -> CoordsBodyToLocal:
    """
    Wrapper for Fortran routine coords_body_to_local

    Parameters
    ----------
    body_position : FloorPositionStruct
        Element body frame coordinates. .r(3)               [x, y, s] position with s = Position from entrance end
        of element.

    ele : EleStruct
        element that local_position coordinates are relative to.

    calculate_angles : bool, optional
        calculate angles for local_position Default: True. False returns local_position angles (.theta, .phi,
        .psi) = 0.

    Returns
    -------
    local_position : FloorPositionStruct
        Local laboratory coordinates. .r(3)               [x, y, s] position with s = Position from entrance end
        of element.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix at to transform vectors. v_local  = w_mat . v_body v_body   = transpose(w_mat) . v_local
    """

class CoordsBodyToRelExit:
    """coords_body_to_rel_exit return type"""

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def rel_exit(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_body_to_rel_exit(body_position: _pybmad.FloorPositionStruct, ele: _pybmad.EleStruct, calculate_angles: bool | None = None) -> CoordsBodyToRelExit:
    """
    Wrapper for Fortran routine coords_body_to_rel_exit

    Parameters
    ----------
    body_position : FloorPositionStruct
        Element body frame coordinates. .r                  [x, y, s] position with s = Position from entrance end
        of element .

    ele : EleStruct
        element that rel_exit coordinates are relative to.

    calculate_angles : bool, optional
        calculate angles for rel_exit Default: True. False returns rel_exit angles (.theta, .phi, .psi) = 0.

    Returns
    -------
    rel_exit : FloorPositionStruct
        Cartesian coordinates relative to exit of the element.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix at to transform vectors. v_rel_exit = w_mat . v_body v_body     = transpose(w_mat) . v_rel_exit
    """

class CoordsCurvilinearToFloor:
    """coords_curvilinear_to_floor return type"""

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_curvilinear_to_floor(xys: Sequence[float], branch: _pybmad.BranchStruct) -> CoordsCurvilinearToFloor:
    """
    Wrapper for Fortran routine coords_curvilinear_to_floor

    Parameters
    ----------
    xys : 1D array of float (shape: 3)
        (x, y, s) lab frame position vector.

    branch : BranchStruct
        Lattice branch that defines the local reference coordinates.

    Returns
    -------
    err_flag : bool
        Set True if global floor position cannot be computed.

    global : FloorPositionStruct
        Global floor position corresponding to (x, y, s) --    .w    -- W matrix to transform vectors: v_global =
        w_mat * v_local
    """

class CoordsFloorToCurvilinear:
    """coords_floor_to_curvilinear return type"""

    @property
    def ele1(self) -> _pybmad.EleStruct | None: ...

    @property
    def status(self) -> int: ...

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def local_coords(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_floor_to_curvilinear(floor_coords: _pybmad.FloorPositionStruct, ele0: _pybmad.EleStruct) -> CoordsFloorToCurvilinear:
    """
    Wrapper for Fortran routine coords_floor_to_curvilinear

    Parameters
    ----------
    floor_coords : FloorPositionStruct
        .r = [X, Y, Z] position in global coordinates

    ele0 : EleStruct
        Element to start the search at.

    Returns
    -------
    status : int
        ok$             -> Local_coords found. patch_problem$  -> No solution due to a patch element. outside$
        -> Outside of lattice ends (for open lattices).

    local_coords : FloorPositionStruct
        .r = [x, y, s] position in curvilinear coordinates with respect to ele1 with s relative to start the
        lattice branch.

    ele1 : EleStruct, optional
        Element that local_coords is with respect to.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix at s, to transform vectors from floor to local. w_mat will only be well defined if status = ok$
    """

class CoordsFloorToLocalCurvilinear:
    """coords_floor_to_local_curvilinear return type"""

    @property
    def status(self) -> int: ...

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def local_position(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_floor_to_local_curvilinear(global_position: _pybmad.FloorPositionStruct, ele: _pybmad.EleStruct, relative_to: int | None = None) -> CoordsFloorToLocalCurvilinear:
    """
    Wrapper for Fortran routine coords_floor_to_local_curvilinear

    Parameters
    ----------
    global_position : FloorPositionStruct
        .r = [X, Y, Z] position in global coordinates

    ele : EleStruct
        element to find local coordinates of.

    relative_to : int, optional
        not_set$ (default), upstream_end$, or downstream_end$. Force which end is used for z = 0. If
        upstream_end$, local_position.r(3) is relative to the upstream end which will not be the entrance end if
        ele.orientation = -1.

    Returns
    -------
    status : int
        longitudinal position: inside$: Inside the element. upstream_end$: At upstream end of element or beyound.
        downstream_end$: At downstream end of element or beyound.

    local_position : FloorPositionStruct
        .r = [x, y, z] position in local curvilinear coordinates.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix at s, to transform vectors. v_global = w_mat.v_local v_local = transpose(w_mat).v_global
    """

def coords_floor_to_relative(floor0: _pybmad.FloorPositionStruct, global_position: _pybmad.FloorPositionStruct, calculate_angles: bool | None = None, is_delta_position: bool | None = None) -> _pybmad.FloorPositionStruct:
    """
    Wrapper for Fortran routine coords_floor_to_relative

    Parameters
    ----------
    floor0 : FloorPositionStruct
        reference position

    global_position : FloorPositionStruct
        global position

    calculate_angles : bool, optional
        calculate angles for local_position Default: True. False returns local_position angles (.theta, .phi,
        .psi) = 0.

    is_delta_position : bool, optional
        If True then treat global_position.r as a difference position in global space and only rotate the position
        but not shift it. Default: False.

    Returns
    -------
    local_position : FloorPositionStruct
        position relative to floor0
    """

class CoordsLocalCurvilinearToBody:
    """coords_local_curvilinear_to_body return type"""

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def body_position(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_local_curvilinear_to_body(local_position: _pybmad.FloorPositionStruct, ele: _pybmad.EleStruct, calculate_angles: bool | None = None) -> CoordsLocalCurvilinearToBody:
    """
    Wrapper for Fortran routine coords_local_curvilinear_to_body

    Parameters
    ----------
    local_position : FloorPositionStruct
        local coordinates. .r(3)               [x, y, s] position with s = Position from entrance end of element.

    ele : EleStruct
        element that coordinates are relative to.

    calculate_angles : bool, optional
        calculate angles for body_position Default: True. False returns body_position angles (.theta, .phi, .psi)
        = 0.

    Returns
    -------
    body_position : FloorPositionStruct
        Element coordinates relative to exit of the element. .r(3)               [x, y, s] position with s =
        Position from entrance end of element.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix at to transform vectors. v_local  = w_mat . v_body v_body   = transpose(w_mat) . v_local
    """

class CoordsLocalCurvilinearToFloor:
    """coords_local_curvilinear_to_floor return type"""

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def global_position(self) -> _pybmad.FloorPositionStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def coords_local_curvilinear_to_floor(local_position: _pybmad.FloorPositionStruct, ele: _pybmad.EleStruct, in_body_frame: bool | None = None, calculate_angles: bool | None = None, end_origin: int | None = None, downstream_dir_ref: bool | None = None) -> CoordsLocalCurvilinearToFloor:
    """
    Wrapper for Fortran routine coords_local_curvilinear_to_floor

    Parameters
    ----------
    local_position : FloorPositionStruct
        Floor position in local curvilinear coordinates, with .r = [x, y, z_local] where z_local is wrt the
        entrance end of the element except if end_origin = downstream_end$. In this case, z_local is a distance
        -ele.value(l$) from the exit end (important for patch elements).

    ele : EleStruct
        element that local_position coordinates are relative to.

    in_body_frame : bool, optional
        True => local_position is in ele body frame and includes misalignments. Ignored if element is a patch.
        Default: False.

    calculate_angles : bool, optional
        calculate angles for global_position Default: True. False returns local_position angles (.theta, .phi,
        .psi) = 0.

    end_origin : int, optional
        not_set$ (default), upstream_end$, or downstream_end$. Force which end is used for z = 0. If
        upstream_end$, local_position.r(3) is relative to the upstream end which will not be the entrance end if
        ele.orientation = -1.

    downstream_dir_ref : bool, optional
        Default False. The output theta angle is calculated so that moduo 2pi this angle is near ele.floor.theta.
        If the element is reversed (ele.direction = -1), the element body coords point upstream which is not
        always wanted. If this arg is set True, ele.floor.theta+pi modulo to be in the range [-pi, pi] is the
        reference.

    Returns
    -------
    global_position : FloorPositionStruct
        Position in global coordinates.

    w_mat : 2D array of float (shape: 3,3), optional
        W matrix at z, to transform vectors. v_global     = w_mat . v_local/body v_local/body = transpose(w_mat) .
        v_global
    """

def coords_relative_to_floor(floor0: _pybmad.FloorPositionStruct, dr: Sequence[float], theta: float | None = None, phi: float | None = None, psi: float | None = None) -> _pybmad.FloorPositionStruct:
    """
    Wrapper for Fortran routine coords_relative_to_floor

    Parameters
    ----------
    floor0 : FloorPositionStruct
        Initial reference frame.

    dr : 1D array of float (shape: 3)
        (x, y, z) positional shift of the reference frame.

    theta : float, optional
        Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
        details. All angles must either be absent or present.

    phi : float, optional
        Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
        details. All angles must either be absent or present.

    psi : float, optional
        Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
        details. All angles must either be absent or present.

    Returns
    -------
    floor1 : FloorPositionStruct
        Shifted reference frame.
    """

def cos_phi(sigma: float, t: float, phi: float, d_param: _pybmad.DiffuseParamStruct) -> float:
    """
    computes  unnormalized cumulative distribution function in phi for a given x
     polar angles relative to surface normal
     azimuthal angle relative to plane of incidence (plane of incoming ray and surface normal)
     1/y suppressed

    Private routine to calculate integrated probability distribution in x = sin(graze_angle_out).
    """

def coulombfun(u: float, v: float, w: float, gam: float) -> float:
    """
    Wrapper for Fortran routine coulombfun

    Parameters
    ----------
    u : float

    v : float

    w : float

    gam : float

    Returns
    -------
    res : float
    """

def create_concatenated_wall3d(lat: _pybmad.LatStruct, err: bool) -> None:
    """
    Routine to concatinate lat%branch(i)ele(:)%wall3d%section(:) arrays into
    one lat%branch(i)%wall3d%section(:) array.

    Exceptions: capillary and aperture elements do not have their walls included.

    Module needed:
      use wall3d_mod

    Parameters
    ----------
    lat : LatStruct
        lattice
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice
    """

class CreateElementSlice:
    """create_element_slice return type"""

    @property
    def sliced_ele(self) -> _pybmad.EleStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def create_element_slice(ele_in: _pybmad.EleStruct, l_slice: float, offset: float, param: _pybmad.LatParamStruct, include_upstream_end: bool, include_downstream_end: bool, old_slice: _pybmad.EleStruct | None = None, orb_in: _pybmad.CoordStruct | None = None) -> CreateElementSlice:
    """
    Wrapper for Fortran routine create_element_slice

    Parameters
    ----------
    ele_in : EleStruct
        Original element to slice

    l_slice : float
        Length of the slice

    offset : float
        Offset of entrance end of sliced_ele from entrance end of ele_in.

    param : LatParamStruct
        lattice paramters.

    include_upstream_end : bool
        Sliced_ele contains the ele's entrance end?

    include_downstream_end : bool
        Sliced_ele contains the ele's exit end?

    old_slice : EleStruct, optional
        Previous slice or, if offset = 0, the previous element. If present this saves computation time of the
        reference energy and time at the start of the present slice. Also makes the ref energy continuous (there
        can be some small differences when using, say, runge_kutta tracking due to tracking tolerances).

    orb_in : CoordStruct, optional
        Incoming orbit if calling routine is doing tracking through the slice. This is used when old_slice is not
        present and there may be an adjustment needed to the orbit ref energy (EG space charge tracking does not
        keep track of ref energy through an lcavity).

    Returns
    -------
    sliced_ele : EleStruct
        Sliced_ele element with appropriate values set.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def create_feedback(lord: _pybmad.EleStruct, input: _pybmad.CharacterAlloc1D, output: _pybmad.CharacterAlloc1D, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine create_feedback

    Parameters
    ----------
    lord : EleStruct
        Feedback element.
        This parameter is an input/output and is modified in-place.
        As an output, lord: Modified feedback elment.

    input : 1D array of str
        Names of input slaves.

    output : 1D array of str
        Names of output slaves.

    err_flag : bool
        Set True if there is a problem.
    """

def create_field_overlap(lat: _pybmad.LatStruct, lord_name: str, slave_name: str) -> bool:
    """
    Wrapper for Fortran routine create_field_overlap

    Parameters
    ----------
    lat : LatStruct
        Lattice

    lord_name : str
        Name of the element with a field extending beyound it's bounds.

    slave_name : str
        Name of the element the lord's field overlaps.

    Returns
    -------
    err_flag : bool
        Set true if there is a problem (like no elements found).
    """

def create_girder(lat: _pybmad.LatStruct, ix_girder: int, contrl: _pybmad.ControlStructArray1D, girder_info: _pybmad.EleStruct, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine create_girder

    Parameters
    ----------
    lat : LatStruct
        Lat to modify.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Modified lattice.

    ix_girder : int
        Index of girder element.

    contrl : 1D array of ControlStruct
        Array of elements that are supported by the girder.

    girder_info : EleStruct
        Element containing attributes to be transfered to the Girder element: girder_info.name girder_info.alias
        girder_info.descrip girder_info.value(:)

    err_flag : bool
    """

def create_group(lord: _pybmad.EleStruct, contrl: _pybmad.ControlStructArray1D, err: bool) -> None:
    """
    Wrapper for Fortran routine create_group

    Parameters
    ----------
    lord : EleStruct
        Group element. .control.type
        This parameter is an input/output and is modified in-place.
        As an output, lord: Modified group elment

    contrl : 1D array of ControlStruct
        control info. 1 element for each slave.

    err : bool
        Set True if an attribute is not free to be controlled.
    """

def create_lat_ele_nametable(lat: _pybmad.LatStruct) -> _pybmad.NametableStruct:
    """
    Wrapper for Fortran routine create_lat_ele_nametable

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    Returns
    -------
    nametable : NametableStruct
        Nametable of the elment names
    """

def create_overlay(lord: _pybmad.EleStruct, contrl: _pybmad.ControlStructArray1D, err: bool) -> None:
    """
    Wrapper for Fortran routine create_overlay

    Parameters
    ----------
    lord : EleStruct
        Overlay element. .control.type
        This parameter is an input/output and is modified in-place.
        As an output, lord: Modified overlay elment

    contrl : 1D array of ControlStruct
        control info. 1 element for each slave.

    err : bool
        Set True if an attribute is not free to be controlled.
    """

class CreatePlanarWigglerModel:
    """create_planar_wiggler_model return type"""

    @property
    def lat(self) -> _pybmad.LatStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def create_planar_wiggler_model(wiggler_in: _pybmad.EleStruct, print_err: bool | None = None) -> CreatePlanarWigglerModel:
    """
    Routine to create series of bend and drift elements to serve as a replacement
    model for a planar wiggler.

    This routine is helpful for translating bmad lattices to a language that does not
    implement the Bmad wiggler model.

    This routine uses the mrqmin nonlinear optimizer to vary the parameters in the wiggler
    model to match:
      Integral g^2 (I_2 radiation integral)
      Integral g^3 (I_3 radiation integral)
      Transfer matrix.
    Also the endding horizontal transverse offset of the reference orbit (floor%r(1)) is
    matched to zero.

    Note: The resulting model does not have the vertical cubic nonlinearity that
    the actual wiggler has.

    Parameters
    ----------
    print_err : bool, optional
        If True (default) print an error message if there is an error.

    Returns
    -------
    lat : LatStruct
        Lattice containing the wiggler model

    err_flag : bool, optional
        Set True if there is an error.
    """

def create_ramper(lord: _pybmad.EleStruct, contrl: _pybmad.ControlStructArray1D, err: bool) -> None:
    """
    Wrapper for Fortran routine create_ramper

    Parameters
    ----------
    lord : EleStruct
        Ramper element. .control.type
        This parameter is an input/output and is modified in-place.
        As an output, lord: Modified ramper elment

    contrl : 1D array of ControlStruct
        control info. 1 element for each slave.

    err : bool
        Set True if an attribute is not free to be controlled.
    """

def create_sol_quad_model(sol_quad: _pybmad.EleStruct, lat: _pybmad.LatStruct) -> None:
    """
    Routine to create series of solenoid and quadrupole elements to serve as a replacement
    model for a sol_quad element.

    This routine is helpful for translating bmad lattices to a language that does not
    implement a combination solenoid/quadrupole.

    Not yet implemented!
    """

def create_unique_ele_names(lat: _pybmad.LatStruct, key: int, suffix: str) -> None:
    """
    Wrapper for Fortran routine create_unique_ele_names

    Parameters
    ----------
    lat : LatStruct
        Lattice holding the elements.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with names made unique.

    key : int
        Class key of elements to consider.

    suffix : str
        Suffix string. Must have a single "?" character.
    """

def create_wiggler_cartesian_map(ele: _pybmad.EleStruct) -> _pybmad.CartesianMapStruct:
    """
    Wrapper for Fortran routine create_wiggler_cartesian_map

    Parameters
    ----------
    ele : EleStruct
        Wiggler or undulator element.

    Returns
    -------
    cart_map : CartesianMapStruct
        Cartesian map.
    """

def crystal_attribute_bookkeeper(ele: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine crystal_attribute_bookkeeper

    Parameters
    ----------
    ele : EleStruct
        Crystal element.
    """

class CrystalDiffractionFieldCalc:
    """crystal_diffraction_field_calc return type"""

    @property
    def e_field(self) -> float: ...

    @property
    def e_phase(self) -> float: ...

    @property
    def orbit_state(self) -> int: ...

    @property
    def dr(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def crystal_diffraction_field_calc(cp: _pybmad.CrystalParamStruct, ele: _pybmad.EleStruct, thickness: float, param: _pybmad.LatParamStruct, p_factor: float, do_branch_calc: bool, follow_undiffracted: bool | None = None) -> CrystalDiffractionFieldCalc:
    """
    Routine to compute the photon field after reflection.

    Parameters
    ----------
    cp : CrystalParamStruct
        Crystal parameters.

    ele : EleStruct
        Crystal element.

    thickness : float
        Crystal thickness

    param : LatParamStruct

    p_factor : float
        Polarization factor.

    do_branch_calc : bool
        Calculate probability of branching to alpha or beta branches?

    follow_undiffracted : bool, optional
        Used with mosaic crystals to calcuate undefracted channel. Default is False.

    Returns
    -------
    e_field : float
        Output field amplitude assuming initial field is 1.

    e_phase : float
        Field phase advance.

    orbit_state : int
        Set to lost$ if crystal is to thick to transmit a photon.

    dr : 1D array of float (shape: 3)
        (x,y,z) orbit change.
    """

def crystal_h_misalign(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, h_vec: Sequence[float]) -> None:
    """
    Routine reorient the crystal H vector due to local imperfections in the crystal lattice.

    Parameters
    ----------
    ele : EleStruct
        Crystal element

    orbit : CoordStruct
        Photon position at crystal surface.

    h_vec : 1D array of float (shape: 3)
        H vector before misalignment.
        This parameter is an input/output and is modified in-place.
        As an output, h_vec: H vector after misalignment.
    """

def crystal_type_to_crystal_params(ele: _pybmad.EleStruct) -> bool:
    """
    Routine to set the crystal parameters based upon the crystal type.

    Crystal types are of the form:
      "ZZZ(ijk)"
    Where "ZZZ" is the atomic formula of the crystal material and "ijk" is the reciprical lattice
    vetor specifying the diffraction plans.

    Parameters
    ----------
    ele : EleStruct
        Crystal element.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Crystal element with computed parameter..

    Returns
    -------
    err_flag : bool
        Set True if crystal type is unrecognized. False otherwise.
    """

def csr_and_sc_apply_kicks(ele: _pybmad.EleStruct, csr: _pybmad.CsrStruct, particle: _pybmad.CoordStructArray1D) -> None:
    """
    Routine to calculate the longitudinal coherent synchrotron radiation kick.

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    csr : CsrStruct

    particle : 1D array of CoordStruct
        Particles to kick.
        This parameter is an input/output and is modified in-place.
        As an output, particle: Particles with kick applied.
    """

def csr_bin_kicks(ele: _pybmad.EleStruct, ds_kick_pt: float, csr: _pybmad.CsrStruct) -> bool:
    """
    Routine to cache intermediate values needed for the csr calculations.

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    ds_kick_pt : float
        Distance between the beginning of the element we are tracking through and the kick point (which is within
        this element).

    csr : CsrStruct

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise
    """

def csr_bin_particles(ele: _pybmad.EleStruct, particle: _pybmad.CoordStructArray1D, err_flag: bool) -> _pybmad.CsrStruct:
    """
    Routine to bin the particles longitudinally in s.

    To avoid noise in the calculation, every particle is considered to have a
    triangular distribution with a base length  given by
      space_charge_com%particle_bin_span * csr%dz_slice.
    That is, particles will, in general, overlap multiple bins.

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    particle : 1D array of CoordStruct
        Array of particles Other:

    Returns
    -------
    csr : CsrStruct
        The bin structure.
    """

class Cumulr:
    """cumulr return type"""

    @property
    def fn(self) -> float: ...

    @property
    def df(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def cumulr(phi: float, status: int) -> Cumulr:
    """
    Wrapper function passed to super_rtsafe by photon_diffuse_scattering.
    Made a module procedure (not nested) to avoid a stack trampoline.
    """

def custom_attribute_ubound_index(ele_class: int) -> int:
    """
    Routine to return, for a given element class, the upper bound index for the ele%custom(:)
    array which is needed to accomodate the registered custom attributes for that class.

    Parameters
    ----------
    ele_class : int
        Element class (key). EG: quadrupole$, etc.

    Returns
    -------
    ix_ubound : int
        Maximum index needed.
    """

class CustomEleAttribNameList:
    """custom_ele_attrib_name_list return type"""

    @property
    def index_list(self) -> _pybmad.IntAlloc1D: ...

    @property
    def name_list(self) -> _pybmad.CharacterAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def custom_ele_attrib_name_list() -> CustomEleAttribNameList:
    """
    Routine to create an array (index_list(i), name_list(i)) of custom element attribute names and indexes.
    Each name in the name_list is of the form:
      "{<class>::}<attribute_name>"
    where:
      <class>:: is an optional class prefix.
      <attribute_name> is the name of the attribute.

    Returns
    -------
    index_list : 1D array of int
        Index of the custom attribute.

    name_list : 1D array of str
        List of custom attributes.
    """

class DIntegral:
    """d_integral return type"""

    @property
    def fn(self) -> float: ...

    @property
    def df(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def d_integral(x: float, status: int) -> DIntegral:
    """
    Wrapper function passed to super_rtsafe by photon_diffuse_scattering.
    Made a module procedure (not nested) to avoid a stack trampoline.
    """

def damping_matrix_d(gamma: float, g_tot: float, B0: float, B1: float, delta: float, species: int) -> list[list[float]]:
    """
    Wrapper for Fortran routine damping_matrix_d

    Parameters
    ----------
    gamma : float

    g_tot : float

    B0 : float

    B1 : float

    delta : float

    species : int

    Returns
    -------
    mat : 2D array of float (shape: 6,6)
    """

def ddz_calc_csr(s_chord_source: float, status: int) -> float:
    """
    Routine to calculate the distance between the source particle and the
    kicked particle at constant time minus the target distance.
    Made a module procedure (not nested) to avoid a stack trampoline. State is passed
    from s_source_calc via the csr_*_ptr / csr_dr_match module variables.

    Parameters
    ----------
    s_chord_source : float
        Chord distance from start of element.

    Returns
    -------
    ddz_this : float
        Distance between source and kick particles: Calculated - Wanted.
    """

def deallocate_ele_pointers(ele: _pybmad.EleStruct, nullify_only: bool | None = None, nullify_branch: bool | None = None, dealloc_poles: bool | None = None) -> None:
    """
    Wrapper for Fortran routine deallocate_ele_pointers

    Parameters
    ----------
    ele : EleStruct
        Element with pointers.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with deallocated pointers.

    nullify_only : bool, optional
        If present and True: Nullify & do not deallocate.

    nullify_branch : bool, optional
        Nullify ele.branch? Default is True.

    dealloc_poles : bool, optional
        Dealloc ele.a/b_pole, ele.a/b_pole_elec? Default is True.
    """

def deallocate_expression_tree(tree: _pybmad.ExpressionTreeStruct) -> None:
    """
    Routine to deallocate an expression tree.

    Parameters
    ----------
    tree : ExpressionTreeStruct
        Tree to deallocate.
        This parameter is an input/output and is modified in-place.
        As an output, tree: Deallocated tree.
    """

def deallocate_lat_pointers(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine deallocate_lat_pointers

    Parameters
    ----------
    lat : LatStruct
        Lat with pointers.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lat with deallocated pointers.
    """

def default_tracking_species(param: _pybmad.LatParamStruct) -> int:
    """
    Wrapper for Fortran routine default_tracking_species

    Parameters
    ----------
    param : LatParamStruct
        Parameters for a lattice branch.

    Returns
    -------
    species : int
        Default species to be used for tracking.
    """

def deposit_particles(xa: _pybmad.RealArray1D, ya: _pybmad.RealArray1D, za: _pybmad.RealArray1D, qa: _pybmad.RealArray1D | None = None, total_charge: float | None = None, resize_mesh: bool | None = None, mesh_growth_factor: float | None = None, mesh_shrink_factor: float | None = None) -> _pybmad.Mesh3DStruct:
    """
    Deposits particle arrays onto mesh

    Parameters
    ----------
    xa : 1D array of float
        x coordinate array

    ya : 1D array of float
        y coordinate array

    za : 1D array of float
        z coordinate array

    qa : 1D array of float, optional
        charge coordinate array

    total_charge : float, optional
        total charge of particles, used only if qa is not present

    resize_mesh : bool, optional
        Set mesh bounds to fit bunch. default  : .true.

    mesh_growth_factor : float, optional
        Fractional padding when growing mesh (default: 0 = tight fit).

    mesh_shrink_factor : float, optional
        Fractional threshold for shrinking mesh (default: 0 = tight fit).

    Returns
    -------
    mesh3d : Mesh3dStruct
        .rho(:,:,:) .charge routine for charge deposition
    """

def detector_pixel_pt(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> list[int]:
    """
    Routine to return the pixel a particle is hitting.

    Parameters
    ----------
    orbit : CoordStruct
        Orbit at surface.

    ele : EleStruct
        Detector element.

    Returns
    -------
    ix_pix : 1D array of int (shape: 2)
        index of ele.photon.pixel.pt(:,:) the particle is in.
    """

def diffraction_plate_or_mask_hit_spot(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct) -> int:
    """
    Wrapper for Fortran routine diffraction_plate_or_mask_hit_spot

    Parameters
    ----------
    ele : EleStruct
        diffraction_plate or mask element.

    orbit : CoordStruct
        particle position.

    Returns
    -------
    ix_section : int
        integer, Set to index of the section where the particle is. Set to zero if the photon is outside all clear
        areas.
    """

def diffusion_matrix_b(gamma: float, g_tot: float, species: int) -> list[list[float]]:
    """
    Wrapper for Fortran routine diffusion_matrix_b

    Parameters
    ----------
    gamma : float

    g_tot : float

    species : int

    Returns
    -------
    mat : 2D array of float (shape: 6,6)
    """

class DistanceToAperture:
    """distance_to_aperture return type"""

    @property
    def no_aperture_here(self) -> bool: ...

    @property
    def dist(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def distance_to_aperture(orbit: _pybmad.CoordStruct, particle_at: int, ele: _pybmad.EleStruct) -> DistanceToAperture:
    """
    Wrapper for Fortran routine distance_to_aperture

    Parameters
    ----------
    orbit : CoordStruct
        Particle position.

    particle_at : int
        first_track_edge$, second_track_edge$, or in_between$

    ele : EleStruct
        Element containing aperture.

    Returns
    -------
    no_aperture_here : bool
        True if aperture does not exist at the longitudinal location of the particle.

    dist : float
        Normalized distance of the particle from the aperture.
    """

def do_mode_flip(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine do_mode_flip

    Parameters
    ----------
    ele : EleStruct
        Starting Element
        This parameter is an input/output and is modified in-place.
        As an output, ele: Flipped element

    Returns
    -------
    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

def dpc_given_de(pc_old: float, mass: float, dE: float) -> float:
    """
    Wrapper for Fortran routine dpc_given_de

    Parameters
    ----------
    pc_old : float

    mass : float

    dE : float

    Returns
    -------
    dpc : float
    """

def drift_and_pipe_track_methods_adjustment(lat: _pybmad.LatStruct) -> None:
    """
    Drift and pipe elements can be used in both photon and non-photon lines.
    A problem occures if, for example, a lattice file with both photon and
    non-photon lines contains a line like:
      drift::*[tracking_method] = taylor
    So this routine resets drift and pipe tracking_method and mat6_calc_method
    parameters in photon lines to bmad_standard if needed.

    Parameters
    ----------
    lat : LatStruct
        Lattice
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with tracking methods adjusted if needed.
    """

def drift_multipass_name_correction(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine drift_multipass_name_correction

    Parameters
    ----------
    lat : LatStruct
    """

def drift_orbit_time(orbit: _pybmad.CoordStruct, beta0: float, delta_s: float | None = None, delta_t: float | None = None) -> None:
    """
    Simple routine to drift a particle orbit in time-based coordinates by a distance delta_s
      or a time delta_t
      If the particle has zero longitudinal velocity, then the particle is not drifted
      and a warning is printed.

    Parameters
    ----------
    orbit : CoordStruct
        particle orbit in time-based coordinates.

    beta0 : float
        reference velocity v/c.

    delta_s : float, optional
        s-coordinate distance to drift particle.

    delta_t : float, optional
        -coordinate distancet to drift particle.
    """

def drift_particle_to_s(p: _pybmad.CoordStruct, s: float, branch: _pybmad.BranchStruct) -> None:
    """
    Drift a particle to a given s-coordinate

    Parameters
    ----------
    p : CoordStruct
        Init particle position.
        This parameter is an input/output and is modified in-place.
        As an output, p: Final particle position.

    s : float
        Target s coordinate.

    branch : BranchStruct
        Branch being tracked through.
    """

def drift_particle_to_t(p: _pybmad.CoordStruct, t: float, branch: _pybmad.BranchStruct) -> None:
    """
    Drift a particle to a given t-coordinate

    Parameters
    ----------
    p : CoordStruct
        Init particle position.
        This parameter is an input/output and is modified in-place.
        As an output, p: Final particle position.

    t : float
        Target t coordinate.

    branch : BranchStruct
        Lattice branch being tracked through.
    """

def dspline_len(s_chord0: float, s_chord1: float, spline: _pybmad.SplineStruct, dtheta_ref: float | None = None) -> float:
    """
    Routine to calculate the difference in length between the spline curve length and a referece line.
    Referece line is centroid chord (referece system of the spline) rotated by dtheta_ref.

    Parameters
    ----------
    s_chord0 : float
        Start position along centroid chord.

    s_chord1 : float
        Stop position along central_chord.

    spline : SplineStruct
        Spline of x-position as a function of s.

    dtheta_ref : float, optional
        angle to rotate the reference line from the centroid chord. Default is 0.

    Returns
    -------
    dlen : float
        L_spline - L_chord
    """

def dynamic_aperture_point(branch: _pybmad.BranchStruct, ele0: _pybmad.EleStruct, orb0: _pybmad.CoordStruct, theta_xy: float, ap_param: _pybmad.ApertureParamStruct, check_xy_init: bool | None = None) -> _pybmad.AperturePointStruct:
    """
    Subroutine to determine one dynamic aperture point by tracking.
    This routine works by determining where on a radial line y = const * x the aperture is.
    Here x and y are deviations from the reference orbit.

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch to track through.

    ele0 : EleStruct
        Lattice element at start of tracking

    orb0 : CoordStruct
        reference orbit at the start of tracking.

    theta_xy : float
        Angle of radial line (in radians) in x-y space. Angle is "normalized" by .x_init, .y_init.

    ap_param : ApertureParamStruct
        Structure holding the input data:

    check_xy_init : bool, optional
        If True, do not check that aperture_param.x_init and .y_init are non-zero. Default is True.

    Returns
    -------
    ap_point : AperturePointStruct
    """

def dynamic_aperture_scan(aperture_param: _pybmad.ApertureParamStruct, pz_start: _pybmad.RealArray1D, lat: _pybmad.LatStruct, print_timing: bool | None = None) -> _pybmad.ApertureScanStructAlloc1D:
    """
    Routine to do a set of dynamic aperture scans.

    Parameters
    ----------
    aperture_param : ApertureParamStruct
        Scan parameters

    pz_start : 1D array of float
        Starting phase space pz values.

    lat : LatStruct
        Lattice.

    print_timing : bool, optional
        If True print info on calculation time. Default is True.

    Returns
    -------
    aperture_scan : 1D array of ApertureScanStruct
        Set of scans. One for each pz_start(:).
    """

def e_accel_field(ele: _pybmad.EleStruct, voltage_or_gradient: int, bmad_standard_tracking: bool | None = None) -> float:
    """
    Wrapper for Fortran routine e_accel_field

    Parameters
    ----------
    ele : EleStruct
        Lcavity or rfcavity element.

    voltage_or_gradient : int
        voltage$ or gradient$

    bmad_standard_tracking : bool, optional
        Using bmad_standard tracking? Default is False.

    Returns
    -------
    field : float
        Cavity field or gradient.
    """

def e_crit_photon(gamma: float, g_bend: float) -> float:
    """
    Routine to calculate the photon critical energy in a bend.

    Parameters
    ----------
    gamma : float
        Gamma factor of charged particle emitting photon.

    g_bend : float
        1/radius bending strength.

    Returns
    -------
    E_crit : float
        Critical photon energy.
    """

class EigenDecomp6mat:
    """eigen_decomp_6mat return type"""

    @property
    def eval(self) -> list[complex]: ...

    @property
    def evec(self) -> list[list[complex]]: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def tunes(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def eigen_decomp_6mat(mat: Sequence[Sequence[float]]) -> EigenDecomp6mat:
    """
    Compute eigenvalues and eigenvectors of a real 6x6 matrix.
    The evals and evecs are in general complex.

    Parameters
    ----------
    mat : 2D array of float (shape: 6,6)
        6x6 real matrix.  Usually a transfer matrix or sigma matrix.

    Returns
    -------
    eval : 1D array of complex (shape: 6)
        complex eigenvalues.

    evec : 2D array of complex (shape: 6,6)
        complex eigenvectors arranged down columns.

    err_flag : bool
        set to true if an error has occured.

    tunes : 1D array of float (shape: 3), optional
        Mode tunes, in radians.
    """

def ele_compute_ref_energy_and_time(ele0: _pybmad.EleStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine ele_compute_ref_energy_and_time

    Parameters
    ----------
    ele0 : EleStruct
        Previous element in lattice with starting energy and time values.

    ele : EleStruct
        Lattice element
        This parameter is an input/output and is modified in-place.
        As an output, ele: Lattice element with reference energy and time.

    param : LatParamStruct
        Lattice parameters.

    err_flag : bool
        Set true if there is an error. False otherwise.
    """

def ele_equal_ele(ele_out: _pybmad.EleStruct, ele_in: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine ele_equal_ele

    Parameters
    ----------
    ele_out : EleStruct

    ele_in : EleStruct
    """

def ele_equals_ele(ele_in: _pybmad.EleStruct, update_nametable: bool) -> _pybmad.EleStruct:
    """
    Subroutine that is used to set an element equal to another.
    Note: Use ele_equal_ele instead unless you know what you are doing.

    Parameters
    ----------
    ele_in : EleStruct
        Input element.

    update_nametable : bool
        If true, update the nametable. If false, do not. Note: nametable updates can take time if this routine is
        called a many times. See remove_eles_from_lat as an example.

    Returns
    -------
    ele_out : EleStruct
        Output element.
    """

def ele_finalizer(ele: _pybmad.EleStruct) -> None:
    """
    Finalizer routine for ele_struct instances.
    NOTE: Not currently used.

    Parameters
    ----------
    ele : EleStruct
        Element to cleanup.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with pointers deallocated as needed.
    """

def ele_full_name(ele: _pybmad.EleStruct, template_: str | None = None) -> str:
    """
    Wrapper for Fortran routine ele_full_name

    Parameters
    ----------
    ele : EleStruct
        Element in a lattice

    Returns
    -------
    str : str
        : Name/location string.
    """

def ele_geometry(floor_start: _pybmad.FloorPositionStruct, ele: _pybmad.EleStruct, len_scale: float | None = None, ignore_patch_err: bool | None = None) -> _pybmad.FloorPositionStruct:
    """
    Wrapper for Fortran routine ele_geometry

    Parameters
    ----------
    floor_start : FloorPositionStruct
        Starting floor coordinates at upstream end. Not used for fiducial and girder elements.

    ele : EleStruct
        Element to propagate the geometry through.

    len_scale : float, optional
        factor to scale the length of the element. 1.0_rp => Output is geometry at end of element (default).
        0.5_rp => Output is geometry at center of element. -1.0_rp => Used to propagate geometry in reverse.

    ignore_patch_err : bool, optional
        If present and True, ignore flexible patch errors. This is used by ele_compute_ref_energy_and_time to
        suppress unnecessary messages.

    Returns
    -------
    floor_end : FloorPositionStruct, optional
        Output floor position. If not present then ele.floor will be used and ele.bookkeeping_state.floor_position
        will be set to ok$.
    """

def ele_geometry_with_misalignments(ele: _pybmad.EleStruct, len_scale: float | None = None) -> _pybmad.FloorPositionStruct:
    """
    Wrapper for Fortran routine ele_geometry_with_misalignments

    Parameters
    ----------
    ele : EleStruct
        Lattice element under consideration.

    len_scale : float, optional
        factor to scale the length of the element. 1.0_rp => Output is geometry at end of element (default).
        0.5_rp => Output is geometry at center of element. -1.0_rp => Used to propagate geometry in reverse.

    Returns
    -------
    floor : FloorPositionStruct
        Floor position with misalignments
    """

def ele_has_constant_ds_dt_ref(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine ele_has_constant_ds_dt_ref

    Parameters
    ----------
    ele : EleStruct
        Element.

    Returns
    -------
    is_const : bool
        True if reference velocity must be a constant.
    """

def ele_has_nonzero_kick(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine ele_has_nonzero_kick

    Parameters
    ----------
    ele : EleStruct
        Element with possible nonzero kicks.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with no kicks.

    Returns
    -------
    has_kick : bool
    """

def ele_has_nonzero_offset(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine ele_has_nonzero_offset

    Parameters
    ----------
    ele : EleStruct
        Element with possible nonzero offsets.

    Returns
    -------
    has_offset : bool
        Set true is element has a non-zero offset.
    """

def ele_is_monitor(ele: _pybmad.EleStruct, print_warning: bool | None = None) -> bool:
    """
    Routine to check that an element is either a detector, instrument, monitor, or marker.
    These are the elements where measurement errors can be defined.

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    print_warning : bool, optional
        If True print a warning message if the element not a monitor like element. Default is True.

    Returns
    -------
    is_monitor : bool
        Set True if the element is a monitor like element.
    """

def ele_loc(ele: _pybmad.EleStruct) -> _pybmad.LatEleLocStruct:
    """
    Wrapper for Fortran routine ele_loc

    Parameters
    ----------
    ele : EleStruct
        Element to be identified

    Returns
    -------
    loc : LatEleLocStruct
        Element identifier.
    """

def ele_loc_name(ele: _pybmad.EleStruct, show_branch0: bool | None = None, parens: str | None = None) -> str:
    """
    Wrapper for Fortran routine ele_loc_name

    Parameters
    ----------
    ele : EleStruct
        Element in a lattice

    show_branch0 : bool, optional
        Explicitly show branch for main lattice elements? Default is False.

    parens : str, optional
        If present, enclose location string using the two characters supplied. Typically parens will be set to
        "()" or "[]".

    Returns
    -------
    str : str
        Output string. Left justified.
    """

class EleMisalignmentLSCalc:
    """ele_misalignment_l_s_calc return type"""

    @property
    def L_mis(self) -> list[float]: ...

    @property
    def S_mis(self) -> list[list[float]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ele_misalignment_l_s_calc(ele: _pybmad.EleStruct) -> EleMisalignmentLSCalc:
    """
    Wrapper for Fortran routine ele_misalignment_l_s_calc

    Parameters
    ----------
    ele : EleStruct
        Element

    Returns
    -------
    L_mis : 1D array of float (shape: 3)
        Misalignment vector relative to center of element

    S_mis : 2D array of float (shape: 3,3)
        Misalignment matrix relative to center of element
    """

def ele_nametable_index(ele: _pybmad.EleStruct) -> int:
    """
    Wrapper for Fortran routine ele_nametable_index

    Parameters
    ----------
    ele : EleStruct
        Element in a lattice.

    Returns
    -------
    ix_nt : int
        Nametable index. lat.nametable.name(ix_nt) and lat.nametable.index(ix_nt) correspond with ele. Set to -1
        if ele is not a lattice element. For example, a slice_slave is not a lattice element.
    """

def ele_order_calc(lat: _pybmad.LatStruct) -> _pybmad.LatEleOrderStruct:
    """
    Wrapper for Fortran routine ele_order_calc

    Parameters
    ----------
    lat : LatStruct
        Lattice to analyze.

    Returns
    -------
    order : LatEleOrderStruct
        Structure holding the element order information.
    """

def ele_reference_energy_correction(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, particle_at: int, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine ele_reference_energy_correction

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    orbit : CoordStruct
        Coordinates to correct.

    particle_at : int
        first_track_edge$ (that is, entering the element), or second_track_edge$ (that is, leaving the element),
        or upstream_end$ (inherit ele.value(p0c_start$) ref), or downstream_end$ (inherit ele.value(p0c$)).
        inside$ (or anything else) -> Do nothing.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before correction.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including correction.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def ele_rf_step_index(E_ref: float, s_rel: float, ele: _pybmad.EleStruct) -> int:
    """
    Wrapper for Fortran routine ele_rf_step_index

    Parameters
    ----------
    E_ref : float
        Reference energy of step. If negative, ignore and use s_rel.

    s_rel : float
        S-position relative to the beginning of the element

    ele : EleStruct
        RF cavity.

    Returns
    -------
    ix_step : int
        Corresponding index in the ele.rf.steps(:) array.
    """

class EleToFibre:
    """ele_to_fibre return type"""

    @property
    def ptc_fibre(self) -> _pybmad.Fibre | None: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ele_to_fibre(ele: _pybmad.EleStruct, use_offsets: bool, integ_order: int | None = None, steps: int | None = None, for_layout: bool | None = None, ref_in: _pybmad.CoordStruct | None = None) -> EleToFibre:
    """
    Wrapper for Fortran routine ele_to_fibre

    Parameters
    ----------
    ele : EleStruct
        Bmad element.

    use_offsets : bool
        Does ptc_fibre include element offsets, pitches and tilt?

    integ_order : int, optional
        Order for the sympletic integrator. Possibilities are: 2, 4, or 6 Overrides ele.value(integrator_order$).
        default = 2 (if not set with set_ptc).

    steps : int, optional
        Number of integration steps. Overrides ele.value(ds_step$).

    for_layout : bool, optional
        If True then fibre will be put in the PTC layout. Default is False.

    ref_in : CoordStruct, optional
        Particle to be tracked. ref_particle$, electron$, etc. This argument should only be present when the fibre
        is not to be put in a layout.

    Returns
    -------
    err_flag : bool
        Set True if setup OK. False otherwise.

    ptc_fibre : Fibre, optional
        PTC fibre element.
    """

def ele_to_ptc_magnetic_bn_an(ele: _pybmad.EleStruct, bn: _pybmad.RealArray1D, an: _pybmad.RealArray1D) -> int:
    """
    Routine to compute the a(n) and b(n) magnetic multipole components of a magnet.
    This is used to interface between eles and PTC fibres

    Note: The multipole index uses the PTC convention of starting from 1 instead of zero.

    Note: On the PTC side bn(1) is error field when creating a fibre but
    is the total field when the fibre is being modified. This routine returns the error field.

    Parameters
    ----------
    ele : EleStruct
        Bmad Element.

    bn : 1D array of float
        Normal multipole component.

    an : 1D array of float
        Skew multipole component.

    Returns
    -------
    n_max : int, optional
        Maximum non-zero multipole component. Set to zero if there are no multipoles.
    """

def ele_to_spin_taylor(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orb0: _pybmad.CoordStruct) -> None:
    """
    Wrapper for Fortran routine ele_to_spin_taylor

    Parameters
    ----------
    ele : EleStruct
        Lattice element.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with spin map.

    param : LatParamStruct
        Branch parameters.

    orb0 : CoordStruct
        Starting ref coords.
    """

class EleToTaylor:
    """ele_to_taylor return type"""

    @property
    def orbital_taylor(self) -> _pybmad.TaylorStructArray1D: ...

    @property
    def spin_taylor(self) -> _pybmad.TaylorStructArray1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ele_to_taylor(ele: _pybmad.EleStruct, orb0: _pybmad.CoordStruct | None = None, taylor_map_includes_offsets: bool | None = None, include_damping: bool | None = None) -> EleToTaylor:
    """
    Wrapper for Fortran routine ele_to_taylor

    Parameters
    ----------
    ele : EleStruct
        Element to construct map for.

    orb0 : CoordStruct, optional
        Starting coords around which the Taylor map is evaluated. Default is the zero orbit.

    taylor_map_includes_offsets : bool, optional
        If present then value overrides ele.taylor_map_includes_offsets.

    include_damping : bool, optional
        Sets if radiation damping is included. Default is what is set in ptc_private.base_state.

    Returns
    -------
    orbital_taylor : 1D array of TaylorStruct (shape: 6), optional
        Orbital taylor map. If not present then the map is put in ele.taylor.

    spin_taylor : 1D array of TaylorStruct (shape: 0:3), optional
        Spin taylor map. If not present then the map is put in ele.spin_taylor.
    """

def ele_unique_name(ele: _pybmad.EleStruct, order: _pybmad.LatEleOrderStruct) -> str:
    """
    Wrapper for Fortran routine ele_unique_name

    Parameters
    ----------
    ele : EleStruct
        Element to construct a unique name for.

    order : LatEleOrderStruct
        Information on element ordering. Before calling this routine, use the routine ele_order_calc to compute
        this argument.

    Returns
    -------
    unique_name : str
        Unique name that can can be used to identify ele. The simplist name will be constructed. For example, if
        the element name is unique, unique_name will be set to the element name.
    """

def ele_value_has_changed(ele: _pybmad.EleStruct, list: _pybmad.IntArray1D, abs_tol: _pybmad.RealArray1D, set_old: bool) -> bool:
    """
    Wrapper for Fortran routine ele_value_has_changed

    Parameters
    ----------
    ele : EleStruct
        Element under consideration.
        This parameter is an input/output and is modified in-place.
        As an output, ele: ele.old_value may be set depending upon setting of set_old

    list : 1D array of int
        List of indexes of ele.value(:) array to check.

    abs_tol : 1D array of float
        List of values such that if the change in parameter value is less than this it is not considered to have
        changed significantly.

    set_old : bool
        If True then set ele.old_value(j) = ele.value(j) for j in list

    Returns
    -------
    has_changed : bool
        Set True if a value has changed significantly.
    """

def ele_vec_equal_ele_vec(ele1: _pybmad.EleStructArray1D, ele2: _pybmad.EleStructArray1D) -> None:
    """
    Wrapper for Fortran routine ele_vec_equal_ele_vec

    Parameters
    ----------
    ele1 : 1D array of EleStruct

    ele2 : 1D array of EleStruct
    """

class ElecMultipoleField:
    """elec_multipole_field return type"""

    @property
    def Ex(self) -> float: ...

    @property
    def Ey(self) -> float: ...

    @property
    def dE(self) -> list[list[float]] | None: ...

    @property
    def compute_dE(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def elec_multipole_field(a: float, b: float, n: int, coord: _pybmad.CoordStruct) -> ElecMultipoleField:
    """
    Wrapper for Fortran routine elec_multipole_field

    Parameters
    ----------
    a : float
        Multipole skew component.

    b : float
        Multipole normal component.

    n : int
        Multipole order.

    coord : CoordStruct

    Returns
    -------
    Ex : float
        X field component

    Ey : float
        Y field component.

    dE : 2D array of float (shape: 2,2), optional
        Field derivatives: dfield(x,y)/d(x,y).

    compute_dE : bool, optional
        If False, do not compute the field derivatives even if dE is present. Default is True.
    """

class ElementAtSBranch:
    """element_at_s_branch return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def s_eff(self) -> float: ...

    @property
    def position(self) -> _pybmad.CoordStruct: ...

    @property
    def ix_ele(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

@overload
def element_at_s(branch: _pybmad.BranchStruct, s: float, choose_max: bool, print_err: bool | None = None) -> ElementAtSBranch:
    """
    Function to return the index of the element at position s.

    element_at_s is an overloaded name for:
      function element_at_s_lat (lat, s, choose_max, ix_branch, err_flag, s_eff, position, print_err) result (ix_ele)
      function element_at_s_branch (branch, s, choose_max, err_flag, s_eff, position, print_err) result (ix_ele)

    The differnce between these two routine is that with element_at_s_lat, the branch is given by the lat
      and ix_ele arguments: branch = lat%branch(ix_ele). With element_at_s_branch, the branch is an argument.

    Also see: pointer_to_element_at_s

    ix_ele is choisen such that:
    If choose_max = True:
        If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
        Else: branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
    If choose_max = False:
        If s = branch%ele(0)%s: ix_ele = 0
        Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
    That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
        choose_max = True  => ix_ele = ix2
        choose_max = False => ix_ele = ix1

    The setting of choose_max only makes a difference when s corresponds to an element boundary.

    Note: For a circular lattice, s is evaluated at the effective s which
    is modulo the branch length:
        s_eff = s - branch_length * floor(s/branch_length)

    Note: If there are multiple elements that are at the given s position due to the presence of
    an element with a negative length, which of the possible elements is actually chosen is ill-defined.

    Parameters
    ----------
    branch : BranchStruct
        Branch to use

    s : float
        Longitudinal position.

    choose_max : bool
        See above

    print_err : bool, optional
        Print error message if there is an error? Default is True.

    Returns
    -------
    ix_ele : int
        Index of element at s.

    err_flag : bool, optional
        Set True if s is out of bounds. False otherwise.

    s_eff : float, optional
        Effective s. Equal to s with a open lattice. See above.

    position : CoordStruct, optional
        Positional information.
    """

@overload
def element_at_s(lat: _pybmad.LatStruct, s: float, choose_max: bool, ix_branch: int | None = None, print_err: bool | None = None) -> ElementAtSLat:
    """
    Function to return the index of the element at position s.

    element_at_s is an overloaded name for:
      function element_at_s_lat (lat, s, choose_max, ix_branch, err_flag, s_eff, position, print_err) result (ix_ele)
      function element_at_s_branch (branch, s, choose_max, err_flag, s_eff, position, print_err) result (ix_ele)

    The differnce between these two routine is that with element_at_s_lat, the branch is given by the lat
      and ix_ele arguments: branch = lat%branch(ix_ele). With element_at_s_branch, the branch is an argument.

    Also see: pointer_to_element_at_s

    ix_ele is choisen such that:
    If choose_max = True:
        If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
        Else: branch%ele(ix_ele)%s_start <= s < branch%ele(ix_ele)%s
    If choose_max = False:
        If s = branch%ele(0)%s: ix_ele = 0
        Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
    That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
        choose_max = True  => ix_ele = ix2
        choose_max = False => ix_ele = ix1

    The setting of choose_max only makes a difference when s corresponds to an element boundary.

    Note: For a circular lattice, s is evaluated at the effective s which
    is modulo the branch length:
        s_eff = s - branch_length * floor(s/branch_length)

    Note: If there are multiple elements that are at the given s position due to the presence of
    an element with a negative length, which of the possible elements is actually chosen is ill-defined.

    Parameters
    ----------
    lat : LatStruct
        Lattice of elements.

    s : float
        Longitudinal position.

    choose_max : bool
        See above

    ix_branch : int, optional
        Branch index. Default is 0.

    print_err : bool, optional
        Print error message if there is an error? Default is True.

    Returns
    -------
    ix_ele : int
        Index of element at s.

    err_flag : bool, optional
        Set True if s is out of bounds. False otherwise.

    s_eff : float, optional
        Effective s. Equal to s with a open lattice. See above.

    position : CoordStruct, optional
        Positional information.
    """

class ElementAtSLat:
    """element_at_s_lat return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def s_eff(self) -> float: ...

    @property
    def position(self) -> _pybmad.CoordStruct: ...

    @property
    def ix_ele(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def element_slice_iterator(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, i_slice: int, n_slice_tot: int, sliced_ele: _pybmad.EleStruct, s_start: float | None = None, s_end: float | None = None) -> None:
    """
    Wrapper for Fortran routine element_slice_iterator

    Parameters
    ----------
    ele : EleStruct
        Element to slice and dice.

    param : LatParamStruct
        Lattice parameters

    i_slice : int
        Slice index

    n_slice_tot : int
        Total number of slices.

    sliced_ele : EleStruct

    s_start : float, optional
        Starting edge of slice relative to beginning of element.

    s_end : float, optional
        Ending edge of slice relative to beginning of element.
    """

def ellipinc_test() -> None:
    """Wrapper for Fortran routine ellipinc_test"""

class EmFieldCalc:
    """em_field_calc return type"""

    @property
    def field(self) -> _pybmad.EmFieldStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def em_field_calc(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, s_pos: float, orbit: _pybmad.CoordStruct, local_ref_frame: bool, calc_dfield: bool | None = None, calc_potential: bool | None = None, use_overlap: bool | None = None, grid_allow_s_out_of_bounds: bool | None = None, rf_time: float | None = None, used_eles: _pybmad.ElePointerStructAlloc1D | None = None, print_err: bool | None = None, original_ele: _pybmad.EleStruct | None = None) -> EmFieldCalc:
    """
    Wrapper for Fortran routine em_field_calc

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    param : LatParamStruct
        Lattice parameters.

    s_pos : float
        Longitudinal position. If local_ref_frame = T: In Body coords relative to the entrance edge of the
        element. If local_ref_frame = F: In Lab coords relative to the upstream edge of the element.

    orbit : CoordStruct
        Transverse coordinates.

    local_ref_frame : bool
        Logical, If True then take the input coordinates and output fields as being with respect to the frame of
        referene of the element (ignore misalignments).

    calc_dfield : bool, optional
        If present and True then calculate the field derivatives.

    calc_potential : bool, optional
        Calc electric and magnetic potentials? Default is false. This is experimental and only implemented for
        wigglers at present.

    use_overlap : bool, optional
        Add in overlap fields from other elements? Default is True.

    grid_allow_s_out_of_bounds : bool, optional
        For grids, allow s-coordinate to be grossly out of bounds and return zero instead of an error? Default:
        False. Used internally for overlapping fields.

    rf_time : float, optional
        Set the time relative to the RF clock. Normally this time is calculated using orbit.t or orbit.vec(5) but
        sometimes it is convenient to be able to override this. For example, time_runge_kutta uses this.

    used_eles : 1D array of ElePointerStruct, optional
        For internal use only when this routine is called recursively. Used to prevent double counting when there
        is field overlap.

    print_err : bool, optional
        Print an error message? Default is True. For example, if the particle is out of bounds when the field is
        defined on a grid.

    original_ele : EleStruct, optional
        Used with recursive calls that pass the lord as the ele argument. In this case original_ele is the
        original ele argument.

    Returns
    -------
    field : EmFieldStruct
        E and B fields and derivatives.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

def em_field_derivatives(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, s_pos: float, orbit: _pybmad.CoordStruct, local_ref_frame: bool, grid_allow_s_out_of_bounds: bool | None = None, rf_time: float | None = None) -> _pybmad.EmFieldStruct:
    """
    Routine to calculate field derivatives.
    In theory this should be handled by em_filed_calc. In practice, em_field_calc is currently incomplete.

    Parameters
    ----------
    ele : EleStruct
        Element

    param : LatParamStruct
        Lattice parameters.

    s_pos : float
        Longitudinal position relative to the upstream edge of the element.

    orbit : CoordStruct
        Transverse coordinates.

    local_ref_frame : bool
        Logical, If True then take the input coordinates and output fields as being with respect to the frame of
        referene of the element (ignore misalignments).

    grid_allow_s_out_of_bounds : bool, optional
        For grids, allow s-coordinate to be grossly out of bounds and return zero instead of an error? Default:
        False. Used internally for overlapping fields.

    rf_time : float, optional
        RF clock time. If not present then the time will be calculated using the standard algorithm.

    Returns
    -------
    dfield : EmFieldStruct
        E and B field derivatives. dfield.E and dfield.B are not touched.
    """

def em_field_kick_vector_time(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, rf_time: float, orbit: _pybmad.CoordStruct, err_flag: bool, print_err: bool | None = None, extra_field: _pybmad.EmFieldStruct | None = None) -> list[float]:
    """
    Subroutine to convert particle coordinates from t-based to s-based system.

    Parameters
    ----------
    ele : EleStruct
        input particle

    param : LatParamStruct
        Reference momentum. The sign indicates direction of p_s.

    rf_time : float
        RF time.

    orbit : CoordStruct
        in t-based system

    err_flag : bool
        Set True if there is an error. False otherwise.

    print_err : bool, optional
        Passed to em_field_calc

    extra_field : EmFieldStruct, optional
        Static field to be added to the element field. Eg used with space charge.

    Returns
    -------
    dvec_dt : 1D array of float (shape: 10)
        Derivatives.
    """

def em_field_plus_em_field(field1: _pybmad.EmFieldStruct, field2: _pybmad.EmFieldStruct) -> _pybmad.EmFieldStruct:
    """
    Wrapper for Fortran routine em_field_plus_em_field

    Parameters
    ----------
    field1 : EmFieldStruct

    field2 : EmFieldStruct

    Returns
    -------
    field_tot : EmFieldStruct
    """

class Emit6d:
    """emit_6d return type"""

    @property
    def mode(self) -> _pybmad.NormalModesStruct: ...

    @property
    def sigma_mat(self) -> list[list[float]]: ...

    @property
    def rad_int_by_ele(self) -> _pybmad.RadIntAllEleStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def emit_6d(ele_ref: _pybmad.EleStruct, include_opening_angle: bool, closed_orbit: _pybmad.CoordStructArray1D | None = None) -> Emit6d:
    """
    Routine to calculate the three normal mode emittances, damping partition numbers, radiation integrals, etc.
    Since the emattances, etc. are only an invariant in the limit of zero damping, the calculated
    values will vary depending upon the reference element.

    If the lattice geometry is open, only the radiation integrals is computed.

    Parameters
    ----------
    ele_ref : EleStruct
        Origin of the 1-turn maps used to evaluate the emittances.

    include_opening_angle : bool
        If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
        comparing against other codes.

    closed_orbit : 1D array of CoordStruct, optional
        Closed orbit. If not present this routine will calculate it.

    Returns
    -------
    mode : NormalModesStruct
        Emittance and other info.

    sigma_mat : 2D array of float (shape: 6,6)
        Sigma matrix.

    rad_int_by_ele : RadIntAllEleStruct, optional
        Radiation integrals element-by-element.
    """

def energy_func(integ_prob: float, status: int) -> float:
    """
    Wrapper for Fortran routine energy_func

    Parameters
    ----------
    integ_prob : float

    status : int

    Returns
    -------
    de : float
    """

def entering_element(orbit: _pybmad.CoordStruct, particle_at: int) -> bool:
    """
    Wrapper for Fortran routine entering_element

    Parameters
    ----------
    orbit : CoordStruct
        Particle orbit.

    particle_at : int
        First_track_edge$ or second_track_edge$

    Returns
    -------
    is_entering : bool
        Set True if particle is going from outside to inside and vice versa.
    """

def envelope_radints(Lambda: Sequence[Sequence[complex]], Theta: Sequence[Sequence[complex]], Iota: Sequence[Sequence[complex]], alpha: Sequence[float], emit: Sequence[float]) -> None:
    """
    Calculates damping decrement and emittance of the three
    normal modes from the integrate diffusion, damping, and vertical
    excitation matrices names Lambda, Theta, and Iota, respectively.
    These three matrices are obtained from the subroutine integrated_mats.

    The damping times can obtained from alpha using:
       tau = lattice_length/c_light/alpha
    """

class EnvelopeRadintsIbs:
    """envelope_radints_ibs return type"""

    @property
    def alpha(self) -> list[float]: ...

    @property
    def emit(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def envelope_radints_ibs(Lambda: Sequence[Sequence[complex]], Theta: Sequence[Sequence[complex]], Iota: Sequence[Sequence[complex]], eles: _pybmad.EleStructArray1D, mode: _pybmad.NormalModesStruct, tail_cut: bool, npart: float, species: int) -> EnvelopeRadintsIbs:
    """
    Calculates damping decrement and emittance of the three
    normal modes by integrating the IBS, SR diffusion, and SR damping matrices.

    The IBS depends on the envelope, and so this routine iterates to
    locate the equilibrium beam envelope. This iterative process can fail to converge.

    The damping times can obtained from alpha using:
       tau = lattice_length/c_light/alpha

    alpha and emit are quantities for the three normal modes.
    alpha and emit are ordered by plane dominance.

    Only radiation from sbends and rbends is taken into account.
    The one-turn transfer matrix at each element (slice) is obtained
    by concatenating the individual element transfer matrices.

    Parameters
    ----------
    Lambda : 2D array of complex (shape: 6,6)
        Integrated damping matrix.

    Theta : 2D array of complex (shape: 6,6)
        Integrated diffusion matrix.

    Iota : 2D array of complex (shape: 6,6)
        Integrated vertical excitation matrix.

    eles : 1D array of EleStruct
        array of element structures representing ring.

    mode : NormalModesStruct
        normal_modes_struct

    tail_cut : bool
        apply tail cut.

    npart : float
        number of particles in bunch.

    species : int
        Particle species.

    Returns
    -------
    alpha : 1D array of float (shape: 3)
        Normal mode damping decrements.

    emit : 1D array of float (shape: 3)
        Normal mode emittances.
    """

def eq_ac_kicker(f1: _pybmad.AcKickerStruct, f2: _pybmad.AcKickerStruct) -> bool:
    """
    Wrapper for Fortran routine eq_ac_kicker

    Parameters
    ----------
    f1 : AcKickerStruct

    f2 : AcKickerStruct

    Returns
    -------
    is_eq : bool
    """

def eq_ac_kicker_freq(f1: _pybmad.AcKickerFreqStruct, f2: _pybmad.AcKickerFreqStruct) -> bool:
    """
    Wrapper for Fortran routine eq_ac_kicker_freq

    Parameters
    ----------
    f1 : AcKickerFreqStruct

    f2 : AcKickerFreqStruct

    Returns
    -------
    is_eq : bool
    """

def eq_ac_kicker_time(f1: _pybmad.AcKickerTimeStruct, f2: _pybmad.AcKickerTimeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_ac_kicker_time

    Parameters
    ----------
    f1 : AcKickerTimeStruct

    f2 : AcKickerTimeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_anormal_mode(f1: _pybmad.AnormalModeStruct, f2: _pybmad.AnormalModeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_anormal_mode

    Parameters
    ----------
    f1 : AnormalModeStruct

    f2 : AnormalModeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_aperture_param(f1: _pybmad.ApertureParamStruct, f2: _pybmad.ApertureParamStruct) -> bool:
    """
    Wrapper for Fortran routine eq_aperture_param

    Parameters
    ----------
    f1 : ApertureParamStruct

    f2 : ApertureParamStruct

    Returns
    -------
    is_eq : bool
    """

def eq_aperture_point(f1: _pybmad.AperturePointStruct, f2: _pybmad.AperturePointStruct) -> bool:
    """
    Wrapper for Fortran routine eq_aperture_point

    Parameters
    ----------
    f1 : AperturePointStruct

    f2 : AperturePointStruct

    Returns
    -------
    is_eq : bool
    """

def eq_aperture_scan(f1: _pybmad.ApertureScanStruct, f2: _pybmad.ApertureScanStruct) -> bool:
    """
    Wrapper for Fortran routine eq_aperture_scan

    Parameters
    ----------
    f1 : ApertureScanStruct

    f2 : ApertureScanStruct

    Returns
    -------
    is_eq : bool
    """

def eq_beam(f1: _pybmad.BeamStruct, f2: _pybmad.BeamStruct) -> bool:
    """
    Wrapper for Fortran routine eq_beam

    Parameters
    ----------
    f1 : BeamStruct

    f2 : BeamStruct

    Returns
    -------
    is_eq : bool
    """

def eq_beam_init(f1: _pybmad.BeamInitStruct, f2: _pybmad.BeamInitStruct) -> bool:
    """
    Wrapper for Fortran routine eq_beam_init

    Parameters
    ----------
    f1 : BeamInitStruct

    f2 : BeamInitStruct

    Returns
    -------
    is_eq : bool
    """

def eq_bmad_common(f1: _pybmad.BmadCommonStruct, f2: _pybmad.BmadCommonStruct) -> bool:
    """
    Wrapper for Fortran routine eq_bmad_common

    Parameters
    ----------
    f1 : BmadCommonStruct

    f2 : BmadCommonStruct

    Returns
    -------
    is_eq : bool
    """

def eq_bookkeeping_state(f1: _pybmad.BookkeepingStateStruct, f2: _pybmad.BookkeepingStateStruct) -> bool:
    """
    Wrapper for Fortran routine eq_bookkeeping_state

    Parameters
    ----------
    f1 : BookkeepingStateStruct

    f2 : BookkeepingStateStruct

    Returns
    -------
    is_eq : bool
    """

def eq_bpm_phase_coupling(f1: _pybmad.BpmPhaseCouplingStruct, f2: _pybmad.BpmPhaseCouplingStruct) -> bool:
    """
    Wrapper for Fortran routine eq_bpm_phase_coupling

    Parameters
    ----------
    f1 : BpmPhaseCouplingStruct

    f2 : BpmPhaseCouplingStruct

    Returns
    -------
    is_eq : bool
    """

def eq_branch(f1: _pybmad.BranchStruct, f2: _pybmad.BranchStruct) -> bool:
    """
    Wrapper for Fortran routine eq_branch

    Parameters
    ----------
    f1 : BranchStruct

    f2 : BranchStruct

    Returns
    -------
    is_eq : bool
    """

def eq_bunch(f1: _pybmad.BunchStruct, f2: _pybmad.BunchStruct) -> bool:
    """
    Wrapper for Fortran routine eq_bunch

    Parameters
    ----------
    f1 : BunchStruct

    f2 : BunchStruct

    Returns
    -------
    is_eq : bool
    """

def eq_bunch_params(f1: _pybmad.BunchParamsStruct, f2: _pybmad.BunchParamsStruct) -> bool:
    """
    Wrapper for Fortran routine eq_bunch_params

    Parameters
    ----------
    f1 : BunchParamsStruct

    f2 : BunchParamsStruct

    Returns
    -------
    is_eq : bool
    """

def eq_cartesian_map(f1: _pybmad.CartesianMapStruct, f2: _pybmad.CartesianMapStruct) -> bool:
    """
    Wrapper for Fortran routine eq_cartesian_map

    Parameters
    ----------
    f1 : CartesianMapStruct

    f2 : CartesianMapStruct

    Returns
    -------
    is_eq : bool
    """

def eq_cartesian_map_term(f1: _pybmad.CartesianMapTermStruct, f2: _pybmad.CartesianMapTermStruct) -> bool:
    """
    Wrapper for Fortran routine eq_cartesian_map_term

    Parameters
    ----------
    f1 : CartesianMapTermStruct

    f2 : CartesianMapTermStruct

    Returns
    -------
    is_eq : bool
    """

def eq_cartesian_map_term1(f1: _pybmad.CartesianMapTerm1Struct, f2: _pybmad.CartesianMapTerm1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_cartesian_map_term1

    Parameters
    ----------
    f1 : CartesianMapTerm1Struct

    f2 : CartesianMapTerm1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_complex_taylor(f1: _pybmad.ComplexTaylorStruct, f2: _pybmad.ComplexTaylorStruct) -> bool:
    """
    Wrapper for Fortran routine eq_complex_taylor

    Parameters
    ----------
    f1 : ComplexTaylorStruct

    f2 : ComplexTaylorStruct

    Returns
    -------
    is_eq : bool
    """

def eq_complex_taylor_term(f1: _pybmad.ComplexTaylorTermStruct, f2: _pybmad.ComplexTaylorTermStruct) -> bool:
    """
    Wrapper for Fortran routine eq_complex_taylor_term

    Parameters
    ----------
    f1 : ComplexTaylorTermStruct

    f2 : ComplexTaylorTermStruct

    Returns
    -------
    is_eq : bool
    """

def eq_control(f1: _pybmad.ControlStruct, f2: _pybmad.ControlStruct) -> bool:
    """
    Wrapper for Fortran routine eq_control

    Parameters
    ----------
    f1 : ControlStruct

    f2 : ControlStruct

    Returns
    -------
    is_eq : bool
    """

def eq_control_ramp1(f1: _pybmad.ControlRamp1Struct, f2: _pybmad.ControlRamp1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_control_ramp1

    Parameters
    ----------
    f1 : ControlRamp1Struct

    f2 : ControlRamp1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_control_var1(f1: _pybmad.ControlVar1Struct, f2: _pybmad.ControlVar1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_control_var1

    Parameters
    ----------
    f1 : ControlVar1Struct

    f2 : ControlVar1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_controller(f1: _pybmad.ControllerStruct, f2: _pybmad.ControllerStruct) -> bool:
    """
    Wrapper for Fortran routine eq_controller

    Parameters
    ----------
    f1 : ControllerStruct

    f2 : ControllerStruct

    Returns
    -------
    is_eq : bool
    """

def eq_coord(f1: _pybmad.CoordStruct, f2: _pybmad.CoordStruct) -> bool:
    """
    Wrapper for Fortran routine eq_coord

    Parameters
    ----------
    f1 : CoordStruct

    f2 : CoordStruct

    Returns
    -------
    is_eq : bool
    """

def eq_coord_array(f1: _pybmad.CoordArrayStruct, f2: _pybmad.CoordArrayStruct) -> bool:
    """
    Wrapper for Fortran routine eq_coord_array

    Parameters
    ----------
    f1 : CoordArrayStruct

    f2 : CoordArrayStruct

    Returns
    -------
    is_eq : bool
    """

def eq_cylindrical_map(f1: _pybmad.CylindricalMapStruct, f2: _pybmad.CylindricalMapStruct) -> bool:
    """
    Wrapper for Fortran routine eq_cylindrical_map

    Parameters
    ----------
    f1 : CylindricalMapStruct

    f2 : CylindricalMapStruct

    Returns
    -------
    is_eq : bool
    """

def eq_cylindrical_map_term(f1: _pybmad.CylindricalMapTermStruct, f2: _pybmad.CylindricalMapTermStruct) -> bool:
    """
    Wrapper for Fortran routine eq_cylindrical_map_term

    Parameters
    ----------
    f1 : CylindricalMapTermStruct

    f2 : CylindricalMapTermStruct

    Returns
    -------
    is_eq : bool
    """

def eq_cylindrical_map_term1(f1: _pybmad.CylindricalMapTerm1Struct, f2: _pybmad.CylindricalMapTerm1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_cylindrical_map_term1

    Parameters
    ----------
    f1 : CylindricalMapTerm1Struct

    f2 : CylindricalMapTerm1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_ele(f1: _pybmad.EleStruct, f2: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine eq_ele

    Parameters
    ----------
    f1 : EleStruct

    f2 : EleStruct

    Returns
    -------
    is_eq : bool
    """

def eq_ellipse_beam_init(f1: _pybmad.EllipseBeamInitStruct, f2: _pybmad.EllipseBeamInitStruct) -> bool:
    """
    Wrapper for Fortran routine eq_ellipse_beam_init

    Parameters
    ----------
    f1 : EllipseBeamInitStruct

    f2 : EllipseBeamInitStruct

    Returns
    -------
    is_eq : bool
    """

def eq_em_field(f1: _pybmad.EmFieldStruct, f2: _pybmad.EmFieldStruct) -> bool:
    """
    Wrapper for Fortran routine eq_em_field

    Parameters
    ----------
    f1 : EmFieldStruct

    f2 : EmFieldStruct

    Returns
    -------
    is_eq : bool
    """

def eq_expression_atom(f1: _pybmad.ExpressionAtomStruct, f2: _pybmad.ExpressionAtomStruct) -> bool:
    """
    Wrapper for Fortran routine eq_expression_atom

    Parameters
    ----------
    f1 : ExpressionAtomStruct

    f2 : ExpressionAtomStruct

    Returns
    -------
    is_eq : bool
    """

def eq_floor_position(f1: _pybmad.FloorPositionStruct, f2: _pybmad.FloorPositionStruct) -> bool:
    """
    Wrapper for Fortran routine eq_floor_position

    Parameters
    ----------
    f1 : FloorPositionStruct

    f2 : FloorPositionStruct

    Returns
    -------
    is_eq : bool
    """

def eq_gen_grad1(f1: _pybmad.GenGrad1Struct, f2: _pybmad.GenGrad1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_gen_grad1

    Parameters
    ----------
    f1 : GenGrad1Struct

    f2 : GenGrad1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_gen_grad_map(f1: _pybmad.GenGradMapStruct, f2: _pybmad.GenGradMapStruct) -> bool:
    """
    Wrapper for Fortran routine eq_gen_grad_map

    Parameters
    ----------
    f1 : GenGradMapStruct

    f2 : GenGradMapStruct

    Returns
    -------
    is_eq : bool
    """

def eq_gg_taylor(f1: _pybmad.GgTaylorStruct, f2: _pybmad.GgTaylorStruct) -> bool:
    """
    Wrapper for Fortran routine eq_gg_taylor

    Parameters
    ----------
    f1 : GgTaylorStruct

    f2 : GgTaylorStruct

    Returns
    -------
    is_eq : bool
    """

def eq_gg_taylor_term(f1: _pybmad.GgTaylorTermStruct, f2: _pybmad.GgTaylorTermStruct) -> bool:
    """
    Wrapper for Fortran routine eq_gg_taylor_term

    Parameters
    ----------
    f1 : GgTaylorTermStruct

    f2 : GgTaylorTermStruct

    Returns
    -------
    is_eq : bool
    """

def eq_grid_beam_init(f1: _pybmad.GridBeamInitStruct, f2: _pybmad.GridBeamInitStruct) -> bool:
    """
    Wrapper for Fortran routine eq_grid_beam_init

    Parameters
    ----------
    f1 : GridBeamInitStruct

    f2 : GridBeamInitStruct

    Returns
    -------
    is_eq : bool
    """

def eq_grid_field(f1: _pybmad.GridFieldStruct, f2: _pybmad.GridFieldStruct) -> bool:
    """
    Wrapper for Fortran routine eq_grid_field

    Parameters
    ----------
    f1 : GridFieldStruct

    f2 : GridFieldStruct

    Returns
    -------
    is_eq : bool
    """

def eq_grid_field_pt(f1: _pybmad.GridFieldPtStruct, f2: _pybmad.GridFieldPtStruct) -> bool:
    """
    Wrapper for Fortran routine eq_grid_field_pt

    Parameters
    ----------
    f1 : GridFieldPtStruct

    f2 : GridFieldPtStruct

    Returns
    -------
    is_eq : bool
    """

def eq_grid_field_pt1(f1: _pybmad.GridFieldPt1Struct, f2: _pybmad.GridFieldPt1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_grid_field_pt1

    Parameters
    ----------
    f1 : GridFieldPt1Struct

    f2 : GridFieldPt1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_high_energy_space_charge(f1: _pybmad.HighEnergySpaceChargeStruct, f2: _pybmad.HighEnergySpaceChargeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_high_energy_space_charge

    Parameters
    ----------
    f1 : HighEnergySpaceChargeStruct

    f2 : HighEnergySpaceChargeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_interval1_coef(f1: _pybmad.Interval1CoefStruct, f2: _pybmad.Interval1CoefStruct) -> bool:
    """
    Wrapper for Fortran routine eq_interval1_coef

    Parameters
    ----------
    f1 : Interval1CoefStruct

    f2 : Interval1CoefStruct

    Returns
    -------
    is_eq : bool
    """

def eq_kv_beam_init(f1: _pybmad.KvBeamInitStruct, f2: _pybmad.KvBeamInitStruct) -> bool:
    """
    Wrapper for Fortran routine eq_kv_beam_init

    Parameters
    ----------
    f1 : KvBeamInitStruct

    f2 : KvBeamInitStruct

    Returns
    -------
    is_eq : bool
    """

def eq_lat(f1: _pybmad.LatStruct, f2: _pybmad.LatStruct) -> bool:
    """
    Wrapper for Fortran routine eq_lat

    Parameters
    ----------
    f1 : LatStruct

    f2 : LatStruct

    Returns
    -------
    is_eq : bool
    """

def eq_lat_ele_loc(f1: _pybmad.LatEleLocStruct, f2: _pybmad.LatEleLocStruct) -> bool:
    """
    Wrapper for Fortran routine eq_lat_ele_loc

    Parameters
    ----------
    f1 : LatEleLocStruct

    f2 : LatEleLocStruct

    Returns
    -------
    is_eq : bool
    """

def eq_lat_param(f1: _pybmad.LatParamStruct, f2: _pybmad.LatParamStruct) -> bool:
    """
    Wrapper for Fortran routine eq_lat_param

    Parameters
    ----------
    f1 : LatParamStruct

    f2 : LatParamStruct

    Returns
    -------
    is_eq : bool
    """

def eq_linac_normal_mode(f1: _pybmad.LinacNormalModeStruct, f2: _pybmad.LinacNormalModeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_linac_normal_mode

    Parameters
    ----------
    f1 : LinacNormalModeStruct

    f2 : LinacNormalModeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_mode3(f1: _pybmad.Mode3Struct, f2: _pybmad.Mode3Struct) -> bool:
    """
    Wrapper for Fortran routine eq_mode3

    Parameters
    ----------
    f1 : Mode3Struct

    f2 : Mode3Struct

    Returns
    -------
    is_eq : bool
    """

def eq_mode_info(f1: _pybmad.ModeInfoStruct, f2: _pybmad.ModeInfoStruct) -> bool:
    """
    Wrapper for Fortran routine eq_mode_info

    Parameters
    ----------
    f1 : ModeInfoStruct

    f2 : ModeInfoStruct

    Returns
    -------
    is_eq : bool
    """

def eq_normal_modes(f1: _pybmad.NormalModesStruct, f2: _pybmad.NormalModesStruct) -> bool:
    """
    Wrapper for Fortran routine eq_normal_modes

    Parameters
    ----------
    f1 : NormalModesStruct

    f2 : NormalModesStruct

    Returns
    -------
    is_eq : bool
    """

def eq_photon_element(f1: _pybmad.PhotonElementStruct, f2: _pybmad.PhotonElementStruct) -> bool:
    """
    Wrapper for Fortran routine eq_photon_element

    Parameters
    ----------
    f1 : PhotonElementStruct

    f2 : PhotonElementStruct

    Returns
    -------
    is_eq : bool
    """

def eq_photon_material(f1: _pybmad.PhotonMaterialStruct, f2: _pybmad.PhotonMaterialStruct) -> bool:
    """
    Wrapper for Fortran routine eq_photon_material

    Parameters
    ----------
    f1 : PhotonMaterialStruct

    f2 : PhotonMaterialStruct

    Returns
    -------
    is_eq : bool
    """

def eq_photon_reflect_surface(f1: _pybmad.PhotonReflectSurfaceStruct, f2: _pybmad.PhotonReflectSurfaceStruct) -> bool:
    """
    Wrapper for Fortran routine eq_photon_reflect_surface

    Parameters
    ----------
    f1 : PhotonReflectSurfaceStruct

    f2 : PhotonReflectSurfaceStruct

    Returns
    -------
    is_eq : bool
    """

def eq_photon_reflect_table(f1: _pybmad.PhotonReflectTableStruct, f2: _pybmad.PhotonReflectTableStruct) -> bool:
    """
    Wrapper for Fortran routine eq_photon_reflect_table

    Parameters
    ----------
    f1 : PhotonReflectTableStruct

    f2 : PhotonReflectTableStruct

    Returns
    -------
    is_eq : bool
    """

def eq_photon_target(f1: _pybmad.PhotonTargetStruct, f2: _pybmad.PhotonTargetStruct) -> bool:
    """
    Wrapper for Fortran routine eq_photon_target

    Parameters
    ----------
    f1 : PhotonTargetStruct

    f2 : PhotonTargetStruct

    Returns
    -------
    is_eq : bool
    """

def eq_pixel_detec(f1: _pybmad.PixelDetecStruct, f2: _pybmad.PixelDetecStruct) -> bool:
    """
    Wrapper for Fortran routine eq_pixel_detec

    Parameters
    ----------
    f1 : PixelDetecStruct

    f2 : PixelDetecStruct

    Returns
    -------
    is_eq : bool
    """

def eq_pixel_pt(f1: _pybmad.PixelPtStruct, f2: _pybmad.PixelPtStruct) -> bool:
    """
    Wrapper for Fortran routine eq_pixel_pt

    Parameters
    ----------
    f1 : PixelPtStruct

    f2 : PixelPtStruct

    Returns
    -------
    is_eq : bool
    """

def eq_pre_tracker(f1: _pybmad.PreTrackerStruct, f2: _pybmad.PreTrackerStruct) -> bool:
    """
    Wrapper for Fortran routine eq_pre_tracker

    Parameters
    ----------
    f1 : PreTrackerStruct

    f2 : PreTrackerStruct

    Returns
    -------
    is_eq : bool
    """

def eq_rad_int1(f1: _pybmad.RadInt1Struct, f2: _pybmad.RadInt1Struct) -> bool:
    """
    Wrapper for Fortran routine eq_rad_int1

    Parameters
    ----------
    f1 : RadInt1Struct

    f2 : RadInt1Struct

    Returns
    -------
    is_eq : bool
    """

def eq_rad_int_all_ele(f1: _pybmad.RadIntAllEleStruct, f2: _pybmad.RadIntAllEleStruct) -> bool:
    """
    Wrapper for Fortran routine eq_rad_int_all_ele

    Parameters
    ----------
    f1 : RadIntAllEleStruct

    f2 : RadIntAllEleStruct

    Returns
    -------
    is_eq : bool
    """

def eq_rad_int_branch(f1: _pybmad.RadIntBranchStruct, f2: _pybmad.RadIntBranchStruct) -> bool:
    """
    Wrapper for Fortran routine eq_rad_int_branch

    Parameters
    ----------
    f1 : RadIntBranchStruct

    f2 : RadIntBranchStruct

    Returns
    -------
    is_eq : bool
    """

def eq_rad_map(f1: _pybmad.RadMapStruct, f2: _pybmad.RadMapStruct) -> bool:
    """
    Wrapper for Fortran routine eq_rad_map

    Parameters
    ----------
    f1 : RadMapStruct

    f2 : RadMapStruct

    Returns
    -------
    is_eq : bool
    """

def eq_rad_map_ele(f1: _pybmad.RadMapEleStruct, f2: _pybmad.RadMapEleStruct) -> bool:
    """
    Wrapper for Fortran routine eq_rad_map_ele

    Parameters
    ----------
    f1 : RadMapEleStruct

    f2 : RadMapEleStruct

    Returns
    -------
    is_eq : bool
    """

def eq_ramper_lord(f1: _pybmad.RamperLordStruct, f2: _pybmad.RamperLordStruct) -> bool:
    """
    Wrapper for Fortran routine eq_ramper_lord

    Parameters
    ----------
    f1 : RamperLordStruct

    f2 : RamperLordStruct

    Returns
    -------
    is_eq : bool
    """

def eq_space_charge_common(f1: _pybmad.SpaceChargeCommonStruct, f2: _pybmad.SpaceChargeCommonStruct) -> bool:
    """
    Wrapper for Fortran routine eq_space_charge_common

    Parameters
    ----------
    f1 : SpaceChargeCommonStruct

    f2 : SpaceChargeCommonStruct

    Returns
    -------
    is_eq : bool
    """

def eq_spin_polar(f1: _pybmad.SpinPolarStruct, f2: _pybmad.SpinPolarStruct) -> bool:
    """
    Wrapper for Fortran routine eq_spin_polar

    Parameters
    ----------
    f1 : SpinPolarStruct

    f2 : SpinPolarStruct

    Returns
    -------
    is_eq : bool
    """

def eq_spline(f1: _pybmad.SplineStruct, f2: _pybmad.SplineStruct) -> bool:
    """
    Wrapper for Fortran routine eq_spline

    Parameters
    ----------
    f1 : SplineStruct

    f2 : SplineStruct

    Returns
    -------
    is_eq : bool
    """

def eq_strong_beam(f1: _pybmad.StrongBeamStruct, f2: _pybmad.StrongBeamStruct) -> bool:
    """
    Wrapper for Fortran routine eq_strong_beam

    Parameters
    ----------
    f1 : StrongBeamStruct

    f2 : StrongBeamStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_curvature(f1: _pybmad.SurfaceCurvatureStruct, f2: _pybmad.SurfaceCurvatureStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_curvature

    Parameters
    ----------
    f1 : SurfaceCurvatureStruct

    f2 : SurfaceCurvatureStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_displacement(f1: _pybmad.SurfaceDisplacementStruct, f2: _pybmad.SurfaceDisplacementStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_displacement

    Parameters
    ----------
    f1 : SurfaceDisplacementStruct

    f2 : SurfaceDisplacementStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_displacement_pt(f1: _pybmad.SurfaceDisplacementPtStruct, f2: _pybmad.SurfaceDisplacementPtStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_displacement_pt

    Parameters
    ----------
    f1 : SurfaceDisplacementPtStruct

    f2 : SurfaceDisplacementPtStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_h_misalign(f1: _pybmad.SurfaceHMisalignStruct, f2: _pybmad.SurfaceHMisalignStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_h_misalign

    Parameters
    ----------
    f1 : SurfaceHMisalignStruct

    f2 : SurfaceHMisalignStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_h_misalign_pt(f1: _pybmad.SurfaceHMisalignPtStruct, f2: _pybmad.SurfaceHMisalignPtStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_h_misalign_pt

    Parameters
    ----------
    f1 : SurfaceHMisalignPtStruct

    f2 : SurfaceHMisalignPtStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_segmented(f1: _pybmad.SurfaceSegmentedStruct, f2: _pybmad.SurfaceSegmentedStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_segmented

    Parameters
    ----------
    f1 : SurfaceSegmentedStruct

    f2 : SurfaceSegmentedStruct

    Returns
    -------
    is_eq : bool
    """

def eq_surface_segmented_pt(f1: _pybmad.SurfaceSegmentedPtStruct, f2: _pybmad.SurfaceSegmentedPtStruct) -> bool:
    """
    Wrapper for Fortran routine eq_surface_segmented_pt

    Parameters
    ----------
    f1 : SurfaceSegmentedPtStruct

    f2 : SurfaceSegmentedPtStruct

    Returns
    -------
    is_eq : bool
    """

def eq_target_point(f1: _pybmad.TargetPointStruct, f2: _pybmad.TargetPointStruct) -> bool:
    """
    Wrapper for Fortran routine eq_target_point

    Parameters
    ----------
    f1 : TargetPointStruct

    f2 : TargetPointStruct

    Returns
    -------
    is_eq : bool
    """

def eq_taylor(f1: _pybmad.TaylorStruct, f2: _pybmad.TaylorStruct) -> bool:
    """
    Wrapper for Fortran routine eq_taylor

    Parameters
    ----------
    f1 : TaylorStruct

    f2 : TaylorStruct

    Returns
    -------
    is_eq : bool
    """

def eq_taylor_term(f1: _pybmad.TaylorTermStruct, f2: _pybmad.TaylorTermStruct) -> bool:
    """
    Wrapper for Fortran routine eq_taylor_term

    Parameters
    ----------
    f1 : TaylorTermStruct

    f2 : TaylorTermStruct

    Returns
    -------
    is_eq : bool
    """

def eq_track(f1: _pybmad.TrackStruct, f2: _pybmad.TrackStruct) -> bool:
    """
    Wrapper for Fortran routine eq_track

    Parameters
    ----------
    f1 : TrackStruct

    f2 : TrackStruct

    Returns
    -------
    is_eq : bool
    """

def eq_track_point(f1: _pybmad.TrackPointStruct, f2: _pybmad.TrackPointStruct) -> bool:
    """
    Wrapper for Fortran routine eq_track_point

    Parameters
    ----------
    f1 : TrackPointStruct

    f2 : TrackPointStruct

    Returns
    -------
    is_eq : bool
    """

def eq_twiss(f1: _pybmad.TwissStruct, f2: _pybmad.TwissStruct) -> bool:
    """
    Wrapper for Fortran routine eq_twiss

    Parameters
    ----------
    f1 : TwissStruct

    f2 : TwissStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wake(f1: _pybmad.WakeStruct, f2: _pybmad.WakeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wake

    Parameters
    ----------
    f1 : WakeStruct

    f2 : WakeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wake_lr(f1: _pybmad.WakeLrStruct, f2: _pybmad.WakeLrStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wake_lr

    Parameters
    ----------
    f1 : WakeLrStruct

    f2 : WakeLrStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wake_lr_mode(f1: _pybmad.WakeLrModeStruct, f2: _pybmad.WakeLrModeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wake_lr_mode

    Parameters
    ----------
    f1 : WakeLrModeStruct

    f2 : WakeLrModeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wake_sr(f1: _pybmad.WakeSrStruct, f2: _pybmad.WakeSrStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wake_sr

    Parameters
    ----------
    f1 : WakeSrStruct

    f2 : WakeSrStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wake_sr_mode(f1: _pybmad.WakeSrModeStruct, f2: _pybmad.WakeSrModeStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wake_sr_mode

    Parameters
    ----------
    f1 : WakeSrModeStruct

    f2 : WakeSrModeStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wake_sr_z_long(f1: _pybmad.WakeSrZLongStruct, f2: _pybmad.WakeSrZLongStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wake_sr_z_long

    Parameters
    ----------
    f1 : WakeSrZLongStruct

    f2 : WakeSrZLongStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wall3d(f1: _pybmad.Wall3DStruct, f2: _pybmad.Wall3DStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wall3d

    Parameters
    ----------
    f1 : Wall3dStruct

    f2 : Wall3dStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wall3d_section(f1: _pybmad.Wall3DSectionStruct, f2: _pybmad.Wall3DSectionStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wall3d_section

    Parameters
    ----------
    f1 : Wall3dSectionStruct

    f2 : Wall3dSectionStruct

    Returns
    -------
    is_eq : bool
    """

def eq_wall3d_vertex(f1: _pybmad.Wall3DVertexStruct, f2: _pybmad.Wall3DVertexStruct) -> bool:
    """
    Wrapper for Fortran routine eq_wall3d_vertex

    Parameters
    ----------
    f1 : Wall3dVertexStruct

    f2 : Wall3dVertexStruct

    Returns
    -------
    is_eq : bool
    """

def eq_xy_disp(f1: _pybmad.XyDispStruct, f2: _pybmad.XyDispStruct) -> bool:
    """
    Wrapper for Fortran routine eq_xy_disp

    Parameters
    ----------
    f1 : XyDispStruct

    f2 : XyDispStruct

    Returns
    -------
    is_eq : bool
    """

def equal_sign_here(ele: _pybmad.EleStruct, delim: str) -> bool:
    """
    Wrapper for Fortran routine equal_sign_here

    Parameters
    ----------
    ele : EleStruct

    delim : str

    Returns
    -------
    is_here : bool
    """

def equivalent_taylor_attributes(ele_taylor: _pybmad.EleStruct, ele2: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine equivalent_taylor_attributes

    Parameters
    ----------
    ele_taylor : EleStruct
        Element with a Taylor map

    ele2 : EleStruct
        Element that might receive the Taylor map from ele_taylor.

    Returns
    -------
    equiv : bool
        True if elements are equivalent.
    """

def etdiv(A: float, B: float, C: float, D: float, E: float, F: float) -> None:
    """
    Wrapper for Fortran routine etdiv

    Parameters
    ----------
    A : float

    B : float

    C : float

    D : float

    E : float

    F : float
    """

class EvaluateArrayIndex:
    """evaluate_array_index return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def word2(self) -> str: ...

    @property
    def delim2(self) -> str: ...

    @property
    def this_index(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def evaluate_array_index(delim_list1: str, delim_list2: str) -> EvaluateArrayIndex:
    """
    Function of evaluate the index of an array. Typically the text being parsed looks like:
         "5) = ..."         or
         "6).COMP = ..."

    Parameters
    ----------
    delim_list1 : str
        Delimitor after the integer. Normally ')'.

    delim_list2 : str
        Delimitor list to mark the end of word2. Normally '='.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.

    word2 : str
        Word found after delim1. Normally this should be blank.

    delim2 : str
        Actual delimitor found after word2.

    this_index : int
        Integer value
    """

class EvaluateLogical:
    """evaluate_logical return type"""

    @property
    def iostat(self) -> int: ...

    @property
    def this_logic(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def evaluate_logical(word: str) -> EvaluateLogical:
    """
    Function of convert a string into a logical value.
    Accepted possibilities are:
      .TRUE.  .FALSE.
       TRUE    FALSE
       T       F

    Parameters
    ----------
    word : str
        Input string.

    Returns
    -------
    iostat : int
        Status: Returns 0 if conversion successful.

    this_logic : bool
        Result.
    """

def exact_bend_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orb: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Subroutine to track through the edge field of an sbend.
    Uses routines adapted from PTC

    Parameters
    ----------
    ele : EleStruct
        SBend element.

    param : LatParamStruct

    particle_at : int
        first_track_edge$, or second_track_edge$.

    orb : CoordStruct
        Coords after tracking.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the edge.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix through the edge.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

def exp_bessi0(t: float, B1: float, B2: float) -> float:
    """
    This is essentially the Numercal Recipes bessi0 function multiplied by exp(-B1*t).

    This overcomes an issue where exp(B2*t) may be huge and exp(-B1*t) may be small.
    Evaluating exp(B2*t) may result in overflow, but exp((B2-B1)*t) has a moderate value.
    Simplifying the algebra of B2-B1 suggests that is should always have a moderate magnitude.

    Parameters
    ----------
    t : float
        Scalar agrument to evaluate function at.

    B1 : float
        Scalar value.  Eq. 33 from Piwinski's paper.

    B2 : float
        Scalar value.  Eq. 34 from Piwinski's paper.
    """

class ExpectOneOf:
    """expect_one_of return type"""

    @property
    def is_ok(self) -> bool: ...

    @property
    def delim(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def expect_one_of(delim_list: str, check_input_delim: bool, ele_name: str, delim: str, delim_found: bool) -> ExpectOneOf:
    """
    Routine to check either that the current delimitor or the next character in the parse stream is the
    expected delimitor.
    This routine is used for Bmad lattice file parsing and is not meant for general use.

    Also see: expect_this

    Parameters
    ----------
    delim_list : str
        List of expected (valid) delimitors. If list contains a space character then no delimitor (indicating the
        end of the command) is a valid possibility.

    check_input_delim : bool
        If True, then check if delim argument is in the delim_list. If False, check that the next character in the
        parse stream is an expected delimitor.

    ele_name : str
        Lattice element under construction. Used for error messages.

    delim : str
        Current delimitor that will be checked if check_input_delim = .true.
        This parameter is an input/output and is modified in-place.
        As an output, delim: Next delim if check_input_delim = False.

    Returns
    -------
    delim : str
        Current delimitor that will be checked if check_input_delim = .true.
        This parameter is an input/output and is modified in-place.
        As an output, delim: Next delim if check_input_delim = False.
    """

class ExpectThis:
    """expect_this return type"""

    @property
    def delim(self) -> str: ...

    @property
    def delim_found(self) -> bool: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def expect_this(expecting: str, check_delim: bool, call_check: bool, err_str: str, ele: _pybmad.EleStruct) -> ExpectThis:
    """
    Checks that the next character or characters in the parse stream corresponds to the
    characters in the expecting argument. For example, if expecting is ')={' these three characters
    should be the next non-blank characters in the parse stream.

    Also see: expect_one_of

    Parameters
    ----------
    expecting : str
        list of characters that are expected to be next in the parse stream.

    check_delim : bool
        If True then use delim argument as first token to check. A blank character indicates end of command is
        expected.

    call_check : bool
        If True then check for 'call::<filename>' construct.

    err_str : str
        String used for error messages.

    ele : EleStruct
        Element parameters being parsed.

    Returns
    -------
    delim : str
        Final delim

    delim_found : bool
        Is there a final delim (as opposed to end of command).
    """

def expression_stack_to_string(stack: _pybmad.ExpressionAtomStructArray1D, polish: bool | None = None) -> str:
    """
    Routine to convert an expression stack to a string

    Parameters
    ----------
    stack : 1D array of ExpressionAtomStruct
        arithmetic expression

    polish : bool, optional
        logical, optional, Construct expression in reverse polish? Default is False.

    Returns
    -------
    str : str
        : Expression in string form.
    """

class ExpressionStackValue:
    """expression_stack_value return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def err_str(self) -> str: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def expression_stack_value(stack: _pybmad.ExpressionAtomStructArray1D, var: _pybmad.ControlVar1StructArray1D | None = None, use_old: bool | None = None) -> ExpressionStackValue:
    """
    Routine to evaluate a mathematical expression represented by an "expression stack".
    Expression stacks are created by expression_string_to_stack.

    Note: Stack elements with stack(i)%type == variable$ need to be evalauated before
    calling this routine and the value placed in stack(i)%value.

    Also see:
      expression_value
      expression_string_to_stack

    Parameters
    ----------
    stack : 1D array of ExpressionAtomStruct
        Expression to evaluate.

    var : 1D array of ControlVar1Struct, optional
        Array of control variables. Used with Bmad controller elements.

    use_old : bool, optional
        Use var.old_value? Must be present if var(:) is present.

    Returns
    -------
    err_flag : bool
        True if there is an evaluation problem. False otherwise.

    err_str : str
        Error string explaining error if there is one.

    value : float
        Value of the expression.
    """

class ExpressionStringToStack:
    """expression_string_to_stack return type"""

    @property
    def stack(self) -> _pybmad.ExpressionAtomStructAlloc1D: ...

    @property
    def n_stack(self) -> int: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def err_str(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def expression_string_to_stack(string: str) -> ExpressionStringToStack:
    """
    This routine creates an expression stack array which can be used
    to evaluate an arithmethic expression.

    Stack end elements not used are marked stack(i)%type = end_stack$

    Stack elements with stack(i)%type = variable$ are elements that need
    to be evaluated before calling expression_stack_value.

    Also see:
      expression_value
      expression_stack_value

    Parameters
    ----------
    string : str
        Expression to be converted.

    Returns
    -------
    stack : 1D array of ExpressionAtomStruct
        Expression evaluation stack.

    n_stack : int
        number of "atoms" used by the expression

    err_flag : bool
        Set True if there is an error (EG divide by 0).

    err_str : str
        String describing the error.
    """

class ExpressionStringToTree:
    """expression_string_to_tree return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def err_str(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def expression_string_to_tree(string: str, root_tree: _pybmad.ExpressionTreeStruct) -> ExpressionStringToTree:
    """
    Routine to create an expression tree array which can be used
    to evaluate an arithmethic expression.

    Also see:
      expression_value
      expression_tree_value
      deallocate_expression_tree

    Important! trees use pointers as opposed to allocatable arrays due to the ifort compiler not being able to
    handle node%node(:) being an allocatable array. Thus deallocate_expression_tree must be called before
    any tree instance goes out of scope.

    Note types used:
      plus$, minus$, times$, divide$, power$, unary_minus$, unary_plus$
      constant$, numeric$, variable$, function$
      root$, parens$, func_parens$, square_brackets$, curly_brackets$
      arrow$, equal$, colon$, double_colon$, vertical_bar$, compound$

    An expression string will be split on:
      Two character operators: "->", "::"
      operators: + -  * / ^ = : &
      brackets: [] () {}
      comma: ,

    Root node name is "root" and is of type root$
    Brackets in the expression string must be matched.
    The corresponding tree node will have a name / type of:
      "[]" / square_brackets$,    "()" / parens$ or func_parens$,   "{}" / curley_brackets$

    The root node, equal nodes, and all bracket nodes, will have an array of child nodes all of which will be comma nodes.

    Parameters
    ----------
    string : str
        Expression to be converted.

    root_tree : ExpressionTreeStruct
        Only used when recursively called.

    Returns
    -------
    err_flag : bool
        Set True if there is an error (EG divide by 0).

    err_str : str
        String describing the error. Make length large to hold the expression.
    """

def expression_tree_to_string(tree: _pybmad.ExpressionTreeStruct, include_root: bool | None = None, n_node: int | None = None, parent: _pybmad.ExpressionTreeStruct | None = None) -> str:
    """
    Routine to convert an expression tree to a expression string.

    Parameters
    ----------
    tree : ExpressionTreeStruct
        Root of tree to print.

    include_root : bool, optional
        Default is True. If True, do not inculde in the output string the root node. Note: If the root node is of
        type root$, this node is always ignored.

    n_node : int, optional
        Node index. parent.node(n_node) === tree. Internal use only. Used with recursive calls.

    parent : ExpressionTreeStruct, optional
        Internal use only. Used with recusive calls.

    Returns
    -------
    str_out : str
        Expression string.
    """

class ExpressionValue:
    """expression_value return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def err_str(self) -> str: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def expression_value(expression: str, var: _pybmad.ControlVar1StructArray1D | None = None, use_old: bool | None = None) -> ExpressionValue:
    """
    Routine to evaluate a mathematical expression encoded in a string.

    Also see:
      expression_string_to_stack
      expression_stack_value

    Parameters
    ----------
    expression : str
        Expression string.

    var : 1D array of ControlVar1Struct, optional
        Array of control variables. Used with Bmad controller elements.

    use_old : bool, optional
        Use var.old_value? Must be present if var(:) is present.

    Returns
    -------
    err_flag : bool
        True if there is an evaluation problem. False otherwise.

    value : float
        Value of the expression.

    err_str : str, optional
        Error string explaining error if there is one.
    """

def fft1(a: _pybmad.RealArray1D, b: _pybmad.RealArray1D, n: int, isn: int) -> int:
    """
    Wrapper for Fortran routine fft1

    Parameters
    ----------
    a : 1D array of float

    b : 1D array of float

    n : int

    isn : int

    Returns
    -------
    ierr : int
    """

class FibreToEle:
    """fibre_to_ele return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def ix_ele(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def fibre_to_ele(ptc_fibre: _pybmad.Fibre, branch: _pybmad.BranchStruct, ix_ele: int, from_mad: bool | None = None) -> FibreToEle:
    """
    Wrapper for Fortran routine fibre_to_ele

    Parameters
    ----------
    ptc_fibre : Fibre
        PTC fibre.

    branch : BranchStruct
        branch containing elements.

    ix_ele : int
        Index in ele(:) array of element last used.
        This parameter is an input/output and is modified in-place.
        As an output, ix_ele: Index to element created (upper index if more than one created).

    from_mad : bool, optional
        If True, ignore PTC specific parameters like integrator_order. Default is False. True is used when the
        fibre has been created via MAD. In this case, the PTC specific parameters may not have good values.

    Returns
    -------
    ix_ele : int
        Index in ele(:) array of element last used.
        This parameter is an input/output and is modified in-place.
        As an output, ix_ele: Index to element created (upper index if more than one created).

    err_flag : bool
        Set true if there is an error. False otherwise. To do: lcavity energy change !? open or closed geometry?
        Energy patch
    """

def field_attribute_free(ele: _pybmad.EleStruct, attrib_name: str) -> bool:
    """
    Routine to check if a field attribute is free to vary.

    Field attributes are either normalized (EG K2 of a sextupole) or unnormalized (EG B2_GRADIENT of a sextupole).
    Whether normalized or unnormalized attributes are free to vary will depend on the setting  of ele%field_master.

    Generally, this routine should not be called directly. Use the routine attribute_free instead.

    Parameters
    ----------
    ele : EleStruct
        Element containing the attribute

    attrib_name : str
        Name of the field attribute. Assumed upper case.

    Returns
    -------
    free : bool
        Is the attribute free to vary? If the attribute is not recognized, free = True will be returned.
    """

def finalize_reflectivity_table(table: _pybmad.PhotonReflectTableStruct, in_degrees: bool) -> None:
    """
    Routine to finalize the construction of the reflectivity tables for a surface.

    Parameters
    ----------
    table : PhotonReflectTableStruct
        Surface tables to be finalized.
        This parameter is an input/output and is modified in-place.
        As an output, table: Finalized surface tables.

    in_degrees : bool
        Table angles in degrees?
    """

class FindElementEnds:
    """find_element_ends return type"""

    @property
    def ele1(self) -> _pybmad.EleStruct | None: ...

    @property
    def ele2(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def find_element_ends(ele: _pybmad.EleStruct, ix_multipass: int | None = None) -> FindElementEnds:
    """
    Wrapper for Fortran routine find_element_ends

    Parameters
    ----------
    ele : EleStruct
        Element to find the ends for.

    ix_multipass : int, optional
        Which multipass pass to follow. Default is 1. This is ignored if there is no multipass elements.

    Returns
    -------
    ele1 : EleStruct, optional
        Pointer to element just before ele.

    ele2 : EleStruct, optional
        Pointer to ele itself or the last sub-element within ele. Note: ele1 and ele2 will be nullified if ele is
        in the lord part of the lattice and does not have any slaves. Note: For an element in the tracking part of
        the lattice: ele1.ix_ele = ele.ix_ele - 1 ele2        => ele Exception: For Beginning element (index 0),
        ele1 => ele
    """

def find_fwhm(bound: float, args: Sequence[float]) -> float:
    """
    Finds the full width at half max of psi(t).  fwhm * c_light / TwoRtTwoLnTwo is taken as the bunch length.

    Steps followed:
      Find value for p(0) that normalizes the solution to dpsi/dt.
      Find max value of p(t) for the value of p(0) found in the previous step.
      Find find tlower, tlower < 0, such that p(tlower) = pmax/2.
      Find find tupper, tupper > 0, such that p(tupper) = pmax/2.
      fwhm is tupper-tlower

    Parameters
    ----------
    bound : float
        -bound and +bound is integration bound.

    args : 1D array of float (shape: 1:8)
        Parameters and constants of dpsi/dt.  See comments of psi_prime for details.

    Returns
    -------
    fwhm : float
        Full width at half max of psi(t)
    """

class FindMatchingFieldmap:
    """find_matching_fieldmap return type"""

    @property
    def match_ele(self) -> _pybmad.EleStruct | None: ...

    @property
    def ix_field(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def find_matching_fieldmap(file_name: str, ele: _pybmad.EleStruct, fm_type: int, ignore_slaves: bool | None = None) -> FindMatchingFieldmap:
    """
    Wrapper for Fortran routine find_matching_fieldmap

    Parameters
    ----------
    file_name : str
        File name associated with field to match to.

    ele : EleStruct
        Element holding the field to be matched.

    fm_type : int
        Type of fieldmap: cartesian_map$, cylindircal_map$, or gen_grad_map$, grid_field$

    ignore_slaves : bool, optional
        If True, ignore any multipass slaves. Default is False.

    Returns
    -------
    ix_field : int
        index of field. For example: matching field => match_ele.cartesian_map(ix_field) Set to -1 if no match
        found.

    match_ele : EleStruct, optional
        Pointer to element with matched field. Nullified if no match found.
    """

def find_normalization(bound: float, p0: float, args: Sequence[float]) -> float:
    """
    Finds value for boundary condition psi(0) that results in integral
    of psi(t) from -bound to +bound to be 1.0.  This is done with the secant method.
    Repeadedly calls integrate_psi with different values for psi(0).

    Parameters
    ----------
    bound : float
        -bound and +bound are integration boundaries

    p0 : float
        Boundary condition psi(0)

    args : 1D array of float (shape: 1:8)
        Parameters and constants of DEQ.  See psi_prime comments for details.

    Returns
    -------
    pnrml : float
        Value for psi(0) that results in integral of psi(t) from -bound to +bound being equal to 1.0
    """

class FloorAnglesToWMat:
    """floor_angles_to_w_mat return type"""

    @property
    def w_mat(self) -> list[list[float]] | None: ...

    @property
    def w_mat_inv(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def floor_angles_to_w_mat(theta: float, phi: float, psi: float) -> FloorAnglesToWMat:
    """
    Wrapper for Fortran routine floor_angles_to_w_mat

    Parameters
    ----------
    theta : float
        Azimuth angle.

    phi : float
        Pitch angle.

    psi : float
        Roll angle.

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3), optional
        Orientation matrix.

    w_mat_inv : 2D array of float (shape: 3,3), optional
        Inverse Orientation matrix.
    """

class FloorWMatToAngles:
    """floor_w_mat_to_angles return type"""

    @property
    def theta(self) -> float: ...

    @property
    def phi(self) -> float: ...

    @property
    def psi(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def floor_w_mat_to_angles(w_mat: Sequence[Sequence[float]], floor0: _pybmad.FloorPositionStruct | None = None) -> FloorWMatToAngles:
    """
    Wrapper for Fortran routine floor_w_mat_to_angles

    Parameters
    ----------
    w_mat : 2D array of float (shape: 3,3)
        Orientation matrix.

    floor0 : FloorPositionStruct, optional
        There are two solutions related by: [theta, phi, psi] & [pi+theta, pi-phi, pi+psi] If floor0 is present,
        choose the solution "nearest" the angles in floor0.

    Returns
    -------
    theta : float
        Azimuth angle.

    phi : float
        Pitch angle.

    psi : float
        Roll angle.
    """

def form_complex_taylor(re_taylor: _pybmad.TaylorStruct, im_taylor: _pybmad.TaylorStruct) -> _pybmad.ComplexTaylorStruct:
    """
    Subroutine to form a complex taylor from two taylor series representing
      the real and imaginary parts

    Parameters
    ----------
    re_taylor : TaylorStruct
        Real part

    im_taylor : TaylorStruct
        Imaginary part

    Returns
    -------
    complex_taylor : ComplexTaylorStruct
        combined complex taylor
    """

class FormDigestedBmadFileName:
    """form_digested_bmad_file_name return type"""

    @property
    def digested_file(self) -> str: ...

    @property
    def full_lat_file(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def form_digested_bmad_file_name(lat_file: str, use_line: str | None = None) -> FormDigestedBmadFileName:
    """
    Subroutine to form the standard name of the Bmad digested file.
    The standard digested file name has the suffix added to the file name:
        suffix = '.digested' + bmad_inc_version$
    Exception: If the use_line argument is present and not blank, the suffix will be:
        suffix = '.' + use_line + '.digested' + bmad_inc_version$

    Parameters
    ----------
    lat_file : str
        Input lattice file name.

    use_line : str, optional
        Line used for lattice expansion. If not present or blank, the line used is the one that was specified in
        the lattice file.

    Returns
    -------
    digested_file : str
        Name of the digested file.

    full_lat_file : str, optional
        Input lattice file name with full directory. Can be used for error messages.
    """

def fringe_here(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, particle_at: int) -> bool:
    """
    Wrapper for Fortran routine fringe_here

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    orbit : CoordStruct
        Particle position.

    particle_at : int
        Either first_track_edge$ or second_track_edge$.

    Returns
    -------
    is_here : bool
        True if there is a fringe. False if not.
    """

def g_bend_from_em_field(b: Sequence[float], e: Sequence[float], orbit: _pybmad.CoordStruct) -> list[float]:
    """
    Routine to calculate the bending strength (1/bending_radius) for a given particle for a given field.
    This will include the dipole bending field of an sbend.

    Parameters
    ----------
    B : 1D array of float (shape: 3)
        Magnetic field.

    E : 1D array of float (shape: 3)
        Electric field

    orbit : CoordStruct
        particle orbit

    Returns
    -------
    g_bend : 1D array of float (shape: 3)
        bending strength vector.
    """

class GBendingStrengthFromEmField:
    """g_bending_strength_from_em_field return type"""

    @property
    def g(self) -> list[float]: ...

    @property
    def dg(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def g_bending_strength_from_em_field(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, s_rel: float, orbit: _pybmad.CoordStruct, local_ref_frame: bool) -> GBendingStrengthFromEmField:
    """
    Wrapper for Fortran routine g_bending_strength_from_em_field

    Parameters
    ----------
    ele : EleStruct
        Element being tracked thorugh.

    param : LatParamStruct
        Lattice parameters.

    s_rel : float
        Distance from the start of the element to the particle.

    orbit : CoordStruct
        Particle position in lab (not element) frame.

    local_ref_frame : bool
        Logical, If True then take the input coordinates and output g as being with respect to the frame of
        referene of the element (ignore misalignments).

    Returns
    -------
    g : 1D array of float (shape: 3)
        g = (g_x, g_y, g_s) bending strength vector (|g| = 1/bend_radius).

    dg : 2D array of float (shape: 3,3), optional
        dg(:)/dr gradient. Takes into account dg_x/dx in a bend due to curvilinear coords.
    """

def g_integrals_calc(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine g_integrals_calc

    Parameters
    ----------
    lat : LatStruct
        Lattice to integrate through.
    """

def gamma_ref(ele: _pybmad.EleStruct) -> float:
    """
    Wrapper for Fortran routine gamma_ref

    Parameters
    ----------
    ele : EleStruct
        Element to evaluate at.

    Returns
    -------
    gamma : float
        Relativistic gamma factor Energy/mass*c^2.
    """

def gen_grad1_to_gg_taylor(ele: _pybmad.EleStruct, gen_grad: _pybmad.GenGradMapStruct, iz: int) -> _pybmad.GgTaylorStructArray1D:
    """
    Wrapper for Fortran routine gen_grad1_to_gg_taylor

    Parameters
    ----------
    ele : EleStruct
        Element containing the map.

    gen_grad : GenGradMapStruct
        Gen_grad map.

    iz : int
        z-plane index to evaluate.

    Returns
    -------
    gg_taylor : 1D array of GgTaylorStruct (shape: 3)
        Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
    """

def gen_grad_at_s_to_gg_taylor(ele: _pybmad.EleStruct, gen_grad: _pybmad.GenGradMapStruct, s_pos: float) -> _pybmad.GgTaylorStructArray1D:
    """
    Wrapper for Fortran routine gen_grad_at_s_to_gg_taylor

    Parameters
    ----------
    ele : EleStruct
        Element containing the map.

    gen_grad : GenGradMapStruct
        Gen_grad map.

    s_pos : float
        Position to evaluate gg_taylor at.

    Returns
    -------
    gg_taylor : 1D array of GgTaylorStruct (shape: 3)
        Map for (Bx, By, Bz) or (Ex, Ey, Ez) fields.
    """

def gen_grad_field(deriv: _pybmad.RealArray1D, gg: _pybmad.GenGrad1Struct, rho: float, theta: float) -> list[float]:
    """
    Wrapper for Fortran routine gen_grad_field

    Parameters
    ----------
    deriv : 1D array of float

    gg : GenGrad1Struct

    rho : float

    theta : float

    Returns
    -------
    field : 1D array of float (shape: 3)
    """

class GetAstraFieldgridNameAndScaling:
    """get_astra_fieldgrid_name_and_scaling return type"""

    @property
    def output_name(self) -> str: ...

    @property
    def field_scale(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def get_astra_fieldgrid_name_and_scaling(ele: _pybmad.EleStruct, name_indexx: _pybmad.StrIndexStruct, dimensions: int | None = None) -> GetAstraFieldgridNameAndScaling:
    """
    Subroutine to get a field grid filename and its scaling. Calls write_astra_field_grid_file.
      If the field grid file does not exist, it is written.

      Note: This is very similar to get_opal_fieldgrid_name_and_scaling

    Parameters
    ----------
    ele : EleStruct
        element to make map

    name_indexx : StrIndexStruct
        contains field grid filenames
        This parameter is an input/output and is modified in-place.
        As an output, name_indexx: updated if new name is added

    dimensions : int, optional
        1 or 3 dimensions. Default: 1

    Returns
    -------
    output_name : str
        output filename.

    field_scale : float
        the scaling of the field grid
    """

def get_bl_from_fwhm(bound: float, args: Sequence[float]) -> float:
    """
    Calculate bunch length as fwhm * c_light / TwoRtTwoLnTwo.
    Where fwhm is full width at half max of solution to dpsi/dt.

    Parameters
    ----------
    bound : float
        -bound and +bound are lower and upper integration bound.

    args : 1D array of float (shape: 1:8)
        Parameters and constants of dpsi/dt.  See comments of psi_prime for details.

    Returns
    -------
    sigma : float
        Bunch length
    """

def get_called_file(delim: str, call_file: str, err: bool) -> None:
    """
    Wrapper for Fortran routine get_called_file

    Parameters
    ----------
    delim : str

    call_file : str

    err : bool
    """

class GetEmitFromSigmaMat:
    """get_emit_from_sigma_mat return type"""

    @property
    def normal(self) -> list[float]: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def get_emit_from_sigma_mat(sigma_mat: Sequence[Sequence[float]], Nmat: Sequence[Sequence[float]] | None = None) -> GetEmitFromSigmaMat:
    r"""
    Given a beam envelop sigma matrix sigma_mat, this returns the 3 normal mode
    emittances.

    The normal mode emittance of the sigma matrix are the eigenvalues of
    sigma_mat . S

    If Nmat is present, then the modes are ordered such that the eigensystem most
    closely resembles Nmat.  If Nmat is not present, then the modes are ordered
    according to which plane they dominate.

        / 0  1  0  0  0  0 \
        |-1  0  0  0  0  0 |
    S = | 0  0  0  1  0  0 |
        | 0  0 -1  0  0  0 |
        | 0  0  0  0  0  1 |
        \ 0  0  0  0 -1  0 /

    Parameters
    ----------
    sigma_mat : 2D array of float (shape: 6,6)
        beam envelop sigma matrix

    Nmat : 2D array of float (shape: 6,6), optional
        If present, then the emittanced will be ordered such that the eigensystem most closely resembles Nmat.

    Returns
    -------
    normal : 1D array of float (shape: 3)
        normal mode emittances

    err_flag : bool
        Set to true if something went wrong.  Otherwise set to false.
    """

class GetGptFieldgridNameAndScaling:
    """get_gpt_fieldgrid_name_and_scaling return type"""

    @property
    def output_name(self) -> str: ...

    @property
    def field_scale(self) -> float: ...

    @property
    def ref_time(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def get_gpt_fieldgrid_name_and_scaling(ele: _pybmad.EleStruct, name_indexx: _pybmad.StrIndexStruct, dimensions: int | None = None) -> GetGptFieldgridNameAndScaling:
    """
    Subroutine to get a field grid filename and its scaling. Calls write_gpt_field_grid_file.
      If the field grid file does not exist, it is written.

      Note: This is very similar to get_opal_fieldgrid_name_and_scaling

    Parameters
    ----------
    ele : EleStruct
        element to make map

    name_indexx : StrIndexStruct
        contains field grid filenames
        This parameter is an input/output and is modified in-place.
        As an output, name_indexx: updated if new name is added

    dimensions : int, optional
        1, 2, or 3 dimensions.

    Returns
    -------
    output_name : str
        output filename.

    field_scale : float
        the scaling of the field grid

    ref_time : float
        time that the field was evaluated at
    """

def get_list_of_names(ele: _pybmad.EleStruct, err_str: str, name_list: _pybmad.CharacterAlloc1D, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

class GetNextWord:
    """get_next_word return type"""

    @property
    def ix_word(self) -> int: ...

    @property
    def delim(self) -> str: ...

    @property
    def delim_found(self) -> bool: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def get_next_word(word: str, delim_list: str, upper_case_word: bool | None = None, call_check: bool | None = None) -> GetNextWord:
    """
    Subroutine to get the next word from the input stream.
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Parameters
    ----------
    word : str
        Word returned

    delim_list : str
        List of valid delimiters

    upper_case_word : bool, optional
        if True then convert word to upper case. Default is True.

    call_check : bool, optional
        If present and True then check for 'call::<filename>' construct. Default is False.

    Returns
    -------
    ix_word : int
        length of word argument

    delim : str
        Actual delimiter found

    delim_found : bool
        Set true if a delimiter found. A delimiter may not be found if the end of the line is reached first.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

class GetOpalFieldgridNameAndScaling:
    """get_opal_fieldgrid_name_and_scaling return type"""

    @property
    def output_name(self) -> str: ...

    @property
    def field_scale(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def get_opal_fieldgrid_name_and_scaling(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, name_indexx: _pybmad.StrIndexStruct) -> GetOpalFieldgridNameAndScaling:
    """
    Subroutine to get a field grid filename and its scaling. Calls write_opal_field_grid_file.
      If the field grid file does not exist, it is written

    Parameters
    ----------
    ele : EleStruct
        element to make map

    param : LatParamStruct
        Contains lattice information

    name_indexx : StrIndexStruct
        contains field grid filenames
        This parameter is an input/output and is modified in-place.
        As an output, name_indexx: updated if new name is added

    Returns
    -------
    output_name : str
        output filename.

    field_scale : float
        the scaling of the field grid
    """

def get_overlay_group_names(ele: _pybmad.EleStruct, lat: _pybmad.LatStruct, pele: _pybmad.ParserEleStruct, delim: str, delim_found: bool, is_control_var_list: bool, err_flag: bool, names_out: _pybmad.CharacterAlloc1D | None = None) -> None:
    """
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Parameters
    ----------
    is_control_var_list : bool
        If True then parsing "var = {...}" list. If False then parsing "group/overlay/girder = {...}" list.
    """

def get_sequence_args(seq_name: str, arg_list: _pybmad.CharacterAlloc1D, delim: str, err_flag: bool) -> None:
    """
    Subroutine to get the argument list for a replacement_line or a list.
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

class GetSlaveList:
    """get_slave_list return type"""

    @property
    def slaves(self) -> _pybmad.ElePointerStructAlloc1D: ...

    @property
    def n_slave(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def get_slave_list(lord: _pybmad.EleStruct) -> GetSlaveList:
    """
    Wrapper for Fortran routine get_slave_list

    Parameters
    ----------
    lord : EleStruct
        The lord element.

    Returns
    -------
    slaves : 1D array of ElePointerStruct
        : Array of slaves.

    n_slave : int
        Number of slaves.
    """

def get_switch(name: str, name_list: _pybmad.CharacterAlloc1D, switch_: int, err: bool, ele: _pybmad.EleStruct, delim: str, delim_found: bool) -> None:
    """
    Wrapper for Fortran routine get_switch

    Parameters
    ----------
    name : str

    name_list : 1D array of str

    err : bool

    ele : EleStruct

    delim : str

    delim_found : bool
    """

def gg_taylor_equal_gg_taylor(gg_taylor1: _pybmad.GgTaylorStruct, gg_taylor2: _pybmad.GgTaylorStruct) -> None:
    """
    Wrapper for Fortran routine gg_taylor_equal_gg_taylor

    Parameters
    ----------
    gg_taylor1 : GgTaylorStruct

    gg_taylor2 : GgTaylorStruct
    """

def gg_taylors_equal_gg_taylors(gg_taylor1: _pybmad.GgTaylorStructArray1D, gg_taylor2: _pybmad.GgTaylorStructArray1D) -> None:
    """
    Wrapper for Fortran routine gg_taylors_equal_gg_taylors

    Parameters
    ----------
    gg_taylor1 : 1D array of GgTaylorStruct

    gg_taylor2 : 1D array of GgTaylorStruct
    """

def gpt_field_grid_scaling(ele: _pybmad.EleStruct, dimensions: int, field_scale: float, ref_time: float) -> None:
    """
    Wrapper for Fortran routine gpt_field_grid_scaling

    Parameters
    ----------
    ele : EleStruct

    dimensions : int

    field_scale : float

    ref_time : float
    """

def gpt_max_field_reference(pt0: _pybmad.GridFieldPt1Struct, ele: _pybmad.EleStruct) -> float:
    """
    Wrapper for Fortran routine gpt_max_field_reference

    Parameters
    ----------
    pt0 : GridFieldPt1Struct

    ele : EleStruct

    Returns
    -------
    field_value : float
    """

class GptToParticleBunch:
    """gpt_to_particle_bunch return type"""

    @property
    def bunch(self) -> _pybmad.BunchStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def gpt_to_particle_bunch(gpt_file: str, ele: _pybmad.EleStruct) -> GptToParticleBunch:
    """
    Routine to initialize a bunch of particles from a GPT screen file.

    Parameters
    ----------
    gpt_file : str
        Name of GPT data file.

    ele : EleStruct
        Lattice element whose downstream end coincident with the GPT screen.

    Returns
    -------
    bunch : BunchStruct
        Particle bunch

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def gradient_shift_sr_wake(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> float:
    """
    Wrapper for Fortran routine gradient_shift_sr_wake

    Parameters
    ----------
    ele : EleStruct
        Lcavity element.

    param : LatParamStruct
        Lattice parameters

    Returns
    -------
    grad_shift : float
        Shift in gradient
    """

def grid_field_interpolate(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, grid: _pybmad.GridFieldStruct, err_flag: bool, x1: float, x2: float | None = None, x3: float | None = None, allow_s_out_of_bounds: bool | None = None, print_err: bool | None = None) -> _pybmad.GridFieldPt1Struct:
    """
    allow_s_out_of_bounds, print_err)

    Subroutine to interpolate the E and B fields on a rectilinear grid.

    Parameters
    ----------
    ele : EleStruct
        Element containing the grid.

    orbit : CoordStruct
        Used for constructing an error message if the particle is out of bounds.

    grid : GridFieldStruct
        Grid to interpolate.

    err_flag : bool
        Set to true if there is an error. False otherwise.

    x1 : float
        dimension 1 interpolation point.

    x2 : float, optional
        dimension 2 interpolation point.

    x3 : float, optional
        dimension 3 interpolation point.

    allow_s_out_of_bounds : bool, optional
        allow s-coordinate grossly out of bounds to return zero field without an error. This is used when the
        field of one element overlaps the field of another. Default is False.

    print_err : bool, optional
        print an error message if the particle is out of bounds? Default is True.
    """

def hard_multipole_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Routine to track through the hard edge field of a multipole.
    The dipole component is ignored and only quadrupole and higher multipoles are included.

    This routine handles elements of type:
      sad_mult, sbend, quadrupole, sextupole

    For sad_mult elements, ele%a_pole and ele%b_pole ae used for the multipole values.
    For the other elements, k1 or k2 is used and it is assumed that we are in the element
    frame of reference so tilt = 0.

    Parameters
    ----------
    ele : EleStruct
        Element with fringe.

    param : LatParamStruct
        Tracking parameters.

    particle_at : int
        Either first_track_edge$ or second_track_edge$.

    orbit : CoordStruct
        Starting coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coordinates.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the fringe.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix including the fringe.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

def has_attribute(ele: _pybmad.EleStruct, attrib: str) -> bool:
    """
    Wrapper for Fortran routine has_attribute

    Parameters
    ----------
    ele : EleStruct

    attrib : str

    Returns
    -------
    has_it : bool
    """

def has_curvature(phot_ele: _pybmad.PhotonElementStruct) -> bool:
    """
    Routine to determine if a surface is potentially curved or is flat.

    Parameters
    ----------
    phot_ele : PhotonElementStruct
        From ele.photon

    Returns
    -------
    curved : bool
        Set True if phot_eleace is curved.
    """

def has_orientation_attributes(ele: _pybmad.EleStruct) -> bool:
    """
    Routine to determine whether an element has orientation attributes like x_offset, etc.
    Also see: has_attribute function.

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    Returns
    -------
    has_attribs : bool
        True if ele has orientation attributes. False otherwise.
    """

def hdf5_read_beam(file_name: str, beam: _pybmad.BeamStruct, error: bool, ele: _pybmad.EleStruct | None = None, pmd_header: _pybmad.PmdHeaderStruct | None = None, print_mom_shift_warning: bool | None = None, conserve_momentum: bool | None = None) -> None:
    """
    Wrapper for Fortran routine hdf5_read_beam

    Parameters
    ----------
    file_name : str

    beam : BeamStruct

    error : bool

    ele : EleStruct, optional

    pmd_header : PmdHeaderStruct, optional

    print_mom_shift_warning : bool, optional

    conserve_momentum : bool, optional
    """

def hdf5_read_grid_field(file_name: str, ele: _pybmad.EleStruct, g_field: _pybmad.GridFieldStructArray1D, err_flag: bool, pmd_header: _pybmad.PmdHeaderStruct | None = None, combine: bool | None = None) -> None:
    """
    Wrapper for Fortran routine hdf5_read_grid_field

    Parameters
    ----------
    file_name : str

    ele : EleStruct

    g_field : 1D array of GridFieldStruct

    err_flag : bool

    pmd_header : PmdHeaderStruct, optional

    combine : bool, optional
    """

def hdf5_write_beam(file_name: str, bunches: _pybmad.BunchStructArray1D, append: bool, error: bool, lat: _pybmad.LatStruct | None = None, alive_only: bool | None = None) -> None:
    """
    Wrapper for Fortran routine hdf5_write_beam

    Parameters
    ----------
    file_name : str

    bunches : 1D array of BunchStruct

    append : bool

    error : bool

    lat : LatStruct, optional

    alive_only : bool, optional
    """

def hdf5_write_grid_field(file_name: str, ele: _pybmad.EleStruct, g_field: _pybmad.GridFieldStructArray1D, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine hdf5_write_grid_field

    Parameters
    ----------
    file_name : str

    ele : EleStruct

    g_field : 1D array of GridFieldStruct

    err_flag : bool
    """

def hwang_bend_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orb: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Subroutine to track through the edge field of an sbend using a 2nd order map.
    Adapted from:
      Hwang and S. Y. Lee,
      "Dipole Fringe Field Thin Map for Compact Synchrotrons",
      Phys. Rev. ST Accel. Beams, 12, 122401, (2015).
    See the Bmad manual for details.

    Parameters
    ----------
    ele : EleStruct
        SBend element.

    param : LatParamStruct
        Rel charge.

    particle_at : int
        first_track_edge$, or second_track_edge$

    orb : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Coords after tracking.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the edge.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix including the edge.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

def i_csr(kick1: _pybmad.CsrKick1Struct, i_bin: int, csr: _pybmad.CsrStruct) -> None:
    """
    Routine to calculate the CSR kick integral (at y = 0)

    Parameters
    ----------
    kick1 : CsrKick1Struct

    i_bin : int
        Bin index.

    csr : CsrStruct
    """

def ibs1(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct, rates: _pybmad.IbsStruct, i: int | None = None, s: float | None = None) -> None:
    """
    Calculates IBS growth rates at some location in a lattice.
    The IBS rates are betatron growth rates.  That is, they are the rate of
    change in sigma_x, sigma_y, and sigma_p.  The emittance growth
    rate is twice the betatron growth rate.
    1/T_emit = 2/T_betatron.
    eg  emit(t) = emit_0 * exp(-2*t/T_betatron) because emit = sigma^2/beta

     Available IBS formulas (ibs_sim_params%formula):
       cimp - Completely Integrated Modified Piwinski
       bjmt - Bjorken-Mtingwa formulation general to bunched beams (time consuming)
       bane - Bane approximation of Bjorken-Mtingwa formulation
       mpzt - Modified Piwinski with Zotter's Integral
       mpxx - Modified Piwinski with a constant Coulomb log.
       kubo - Kubo and Oide's sigma matrix-based

    Either i or s, but not both, must be specified.
    """

def ibs_blowup1turn(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct) -> None:
    """
    Updates beam emittances with effect of IBS for
    one turn on the lattice.

    Parameters
    ----------
    lat : LatStruct
        lattice

    ibs_sim_params : IbsSimParamStruct
        Parameters for calculation of IBS rates
    """

class IbsDeltaCalc:
    """ibs_delta_calc return type"""

    @property
    def delta_sigma_energy(self) -> float: ...

    @property
    def delta_emit_a(self) -> float: ...

    @property
    def delta_emit_b(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ibs_delta_calc(lat: _pybmad.LatStruct, ix: int, ibs_sim_params: _pybmad.IbsSimParamStruct, sigma_mat: Sequence[Sequence[float]] | None = None) -> IbsDeltaCalc:
    """
    Calculates change in energy spread and emittances due to IBS for a single element.

    Parameters
    ----------
    lat : LatStruct
        lattice for tracking

    ix : int
        index of element to use: lat.ele(ix)

    ibs_sim_params : IbsSimParamStruct
        parameters for calculation of IBS rates.

    sigma_mat : 2D array of float (shape: 6,6), optional
        Beam's sigma matrix. Required for 'kubo' method.

    Returns
    -------
    delta_sigma_energy : float, optional
        change in energy spread in eV

    delta_emit_a : float, optional
        change in a-mode emittance (geometric)

    delta_emit_b : float, optional
        change in b-mode emittance (geometric)
    """

def ibs_equib_der(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct, inmode: _pybmad.NormalModesStruct, granularity: float) -> _pybmad.NormalModesStruct:
    """
    Computes equilibrium beam sizes by calculating emittance growth rates from IBS growth rates.
    Steps beam size through time till equilibrium is reached.

    Parameters
    ----------
    lat : LatStruct
        lattice for tracking

    ibs_sim_params : IbsSimParamStruct
        parameters for IBS calculation

    inmode : NormalModesStruct
        natural beam parameters

    granularity : float
        Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter. Set to -1 to calculate
        element-by-element.

    Returns
    -------
    ibsmode : NormalModesStruct
        beam parameters after IBS effects
    """

def ibs_equib_rlx(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct, inmode: _pybmad.NormalModesStruct, ratio: float, initial_blow_up: Sequence[float], granularity: float) -> _pybmad.NormalModesStruct:
    """
    $OMP THREADPRIVATE(bl_lat_ptr, bl_ibs_params_ptr, bl_mode_ptr)

     Subroutine ibs_equib_rlx(lat,ibs_sim_params,inmode,ibsmode,ratio,initial_blow_up,granularity)
     Iterates to equilibrium beam conditions using relaxation method

     This method requires that the initial beam size be larger than the equilibrium beam size.
     An initial_blow_up of 3 to 5 is a good place to start.

     See ibs_rates subroutine for available IBS rate formulas.

    Parameters
    ----------
    lat : LatStruct
        lattice for tracking

    ibs_sim_params : IbsSimParamStruct
        parameters for IBS calculation

    inmode : NormalModesStruct
        natural beam parameters

    ratio : float
        Ratio of vert_emit_coupling / vert_emit_total

    initial_blow_up : 1D array of float (shape: 3)
        Factor multiplied to all thre bunch dimensions prior to starting iteration.

    granularity : float
        Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter.

    Returns
    -------
    ibsmode : NormalModesStruct
        beam parameters after IBS effects
    """

def ibs_lifetime(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct, maxratio: _pybmad.IbsMaxratioStruct, granularity: float) -> _pybmad.IbsLifetimeStruct:
    """
    This module computes the beam lifetime due to
    the diffusion process according to equation 12
    from page 129 of The Handbook for Accelerator
    Physics and Engineering 2nd edition.

    Parameters
    ----------
    lat : LatStruct
        lattice for tracking.

    ibs_sim_params : IbsSimParamStruct
        parameters for calculation of IBS rates.

    maxratio : IbsMaxratioStruct
        Ax,y,p/sigma_x,y,p where Ax,y,p is the maximum sigma.  Note that this quantity is just the ratio, not the
        ratio squared.  For example, maxratio%Rx = 1.1 says that the maximum acceptable beamsize is 10% larger
        than the beamsize before IBS effects.

    granularity : float
        Step size when slicing lattice.  -1 for element-by-element.

    Returns
    -------
    lifetime : IbsLifetimeStruct
        structure returning IBS lifetimes
    """

def ibs_matrix_c(sigma_mat: Sequence[Sequence[float]], tail_cut: bool, tau: float, energy: float, n_part: float, species: int) -> list[list[float]]:
    """
    Wrapper for Fortran routine ibs_matrix_c

    Parameters
    ----------
    sigma_mat : 2D array of float (shape: 6,6)

    tail_cut : bool

    tau : float

    energy : float

    n_part : float

    species : int

    Returns
    -------
    ibs_mat : 2D array of float (shape: 6,6)
    """

def ibs_rates1turn(lat: _pybmad.LatStruct, ibs_sim_params: _pybmad.IbsSimParamStruct, granularity: float) -> _pybmad.IbsStruct:
    """
    Calculates IBS risetimes for given lat
    This is basically a front-end for the various formulas
    available in this module of calculating IBS rates.

    Parameters
    ----------
    lat : LatStruct
        lattice for tracking.

    ibs_sim_params : IbsSimParamStruct
        parameters for IBS calculation.

    granularity : float
        slice length.  -1 for element-by-element.

    Returns
    -------
    rates1turn : IbsStruct
        ibs rates for onr turn on the lattice.
    """

def igfcoulombfun(u: float, v: float, w: float, gam: float, dx: float, dy: float, dz: float) -> float:
    """
    Wrapper for Fortran routine igfcoulombfun

    Parameters
    ----------
    u : float

    v : float

    w : float

    gam : float

    dx : float

    dy : float

    dz : float

    Returns
    -------
    res : float
    """

def igfexfun(u: float, v: float, w: float, gam: float, dx: float, dy: float, dz: float) -> float:
    """
    Wrapper for Fortran routine igfexfun

    Parameters
    ----------
    u : float

    v : float

    w : float

    gam : float

    dx : float

    dy : float

    dz : float

    Returns
    -------
    res : float
    """

def igfeyfun(u: float, v: float, w: float, gam: float, dx: float, dy: float, dz: float) -> float:
    """
    Wrapper for Fortran routine igfeyfun

    Parameters
    ----------
    u : float

    v : float

    w : float

    gam : float

    dx : float

    dy : float

    dz : float

    Returns
    -------
    res : float
    """

def igfezfun(u: float, v: float, w: float, gam: float, dx: float, dy: float, dz: float) -> float:
    """
    Wrapper for Fortran routine igfezfun

    Parameters
    ----------
    u : float

    v : float

    w : float

    gam : float

    dx : float

    dy : float

    dz : float

    Returns
    -------
    res : float
    """

def image_charge_kick_calc(kick1: _pybmad.CsrKick1Struct, csr: _pybmad.CsrStruct) -> None:
    """
    Routine to calculate the image charge kick.

    Parameters
    ----------
    kick1 : CsrKick1Struct

    csr : CsrStruct
    """

class InitAttributeName1:
    """init_attribute_name1 return type"""

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def init_attribute_name1(is_ok: bool, ix_key: int, ix_attrib: int, name: str, attrib_state: int | None = None, override: bool | None = None) -> InitAttributeName1:
    """
    Routine to initialize a single name in the element attribute name table.

    Parameters
    ----------
    is_ok : bool
        Initial setting.
        This parameter is an input/output and is modified in-place.
        As an output, is_ok: Set False if there is a problem. Otherwise untouched.

    ix_key : int
        Key index.

    ix_attrib : int
        Attribute index.

    name : str
        Attribute name. Should be uppercase if attrib_state = is_free$. Should contain non-uppercase characters if
        attrib_state = private$.

    attrib_state : int, optional
        Class of attribute: does_not_exist$, is_free$, etc. Defaults to is_free$.

    override : bool, optional
        Normally this routine throws an error if the [ix_key, ix_attrib] has been set previously. If override =
        True then the set is done and no error is generated.

    Returns
    -------
    is_ok : bool
        Initial setting.
        This parameter is an input/output and is modified in-place.
        As an output, is_ok: Set False if there is a problem. Otherwise untouched.
    """

def init_attribute_name_array() -> None:
    """
    Private routine to initialize the attribute name array used by routines
    in attribute_mod. Not meant for general use.
    """

class InitBeamDistribution:
    """init_beam_distribution return type"""

    @property
    def beam(self) -> _pybmad.BeamStruct: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def beam_init_set(self) -> _pybmad.BeamInitStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def init_beam_distribution(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, beam_init: _pybmad.BeamInitStruct, modes: _pybmad.NormalModesStruct | None = None, print_p0c_shift_warning: bool | None = None, conserve_momentum: bool | None = None) -> InitBeamDistribution:
    """
    print_p0c_shift_warning, conserve_momentum)

    Subroutine to initialize a beam of particles.
    Initialization uses the downstream parameters of ele.

    Note: This routine sets the random number generator according to the settings
    in beam_int and at the end resets things to their initial state.

    For more information on individual bunch initialization, see the
    init_bunch_distribution routine.

    Note: The optional "modes" argument generally is used to pass in normal mode parameters as
    calculated from the lattice. If present, and if a parameter like beam_init%a_emit are
    set negative, then the corresponding parameter in the modes structure is used.
    If not present, a warning message is issued and the parameter is set to zero.
    This is only used for parameters that cannot be negative.

    Parameters
    ----------
    ele : EleStruct
        element to initialize distribution at (downstream end).

    param : LatParamStruct
        Lattice parameters

    beam_init : BeamInitStruct
        Use "getf beam_init_struct" for more details.

    modes : NormalModesStruct, optional
        Normal mode parameters. See above.

    print_p0c_shift_warning : bool, optional
        Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

    Returns
    -------
    beam : BeamStruct
        Structure with initialized particles.

    err_flag : bool, optional
        Set true if there is an error, false otherwise.

    beam_init_set : BeamInitStruct, optional
        Set to input beam_init with components like .a_emit set what is used in constructing the beam (which is
        different from beam_init.a_emit if this is set negative).
    """

def init_bmad() -> None:
    """Wrapper for Fortran routine init_bmad"""

def init_bmad_parser_common(lat: _pybmad.LatStruct | None = None) -> None:
    """
    Wrapper for Fortran routine init_bmad_parser_common

    Parameters
    ----------
    lat : LatStruct, optional
    """

class InitBunchDistribution:
    """init_bunch_distribution return type"""

    @property
    def bunch(self) -> _pybmad.BunchStruct: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def beam_init_used(self) -> _pybmad.BeamInitStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def init_bunch_distribution(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, beam_init: _pybmad.BeamInitStruct, ix_bunch: int, modes: _pybmad.NormalModesStruct | None = None, print_p0c_shift_warning: bool | None = None, conserve_momentum: bool | None = None) -> InitBunchDistribution:
    r"""
    print_p0c_shift_warning, conserve_momentum)

    Subroutine to initialize a distribution of particles of a bunch.
    Initialization uses the downstream parameters of ele.

    There are four distributions available:
      \'\', or 'ran_gauss' -- Random gaussian distribution.
      'ellipse'  -- concentric ellipses representing a Gaussian distribution
      'grid'     -- uniform rectangular grid
      'KV'       -- Kapchinsky-Vladimirsky distribution
    See the Bmad manual for more information.

    The distribution is matched to the Twiss parameters, centroid position, and Energy - z
    correlation as specified. Coupling in the element ele is incorporated into the distribution.

    Note: Except for the random number seed, the random number generator
    parameters used for this routine are set from the beam_init argument.
    That is, these parameters are independent of what is used everywhere else.

    Note: Make sure: |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z

    Note: The optional "modes" argument generally is used to pass in normal mode parameters as
    calculated from the lattice. If present, and if a parameter like beam_init%a_emit are
    set negative, then the corresponding parameter in the modes structure is used.
    If not present, a warning message is issued and the parameter is set to zero.
    This is only used for parameters that cannot be negative.

    Note: To get good results, It is important to make sure that for
    circular rings that beam_init%center is the correct closed orbit.
    The closed orbit will shift if, for example, radiation damping is turned on.

    Parameters
    ----------
    ele : EleStruct
        element to initialize distribution at (downstream end).

    param : LatParamStruct
        Lattice parameters

    beam_init : BeamInitStruct
        Use "getf beam_init_struct" for more details.

    ix_bunch : int
        Bunch index. 0 = bunch generated at time = 0.

    modes : NormalModesStruct, optional
        Normal mode parameters. See above.

    print_p0c_shift_warning : bool, optional
        Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

    Returns
    -------
    bunch : BunchStruct
        Structure with initialized particles.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.

    beam_init_used : BeamInitStruct, optional
        Set to input beam_init with components like .a_emit set what is used in constructing the beam (which can
        be different from beam_init.a_emit if this is set negative). If reading from a file, beam_init_used will
        equal beam_init.
    """

def init_complex_taylor_series(complex_taylor: _pybmad.ComplexTaylorStruct, n_term: int, save: bool | None = None) -> None:
    """
    Subroutine to initialize a Bmad complex_taylor series (6 of these series make
    a complex_taylor map). Note: This routine does not zero the structure. The calling
    routine is responsible for setting all values.

    Parameters
    ----------
    complex_taylor : ComplexTaylorStruct
        Old structure.
        This parameter is an input/output and is modified in-place.
        As an output, complex_taylor: Initalized structure.

    n_term : int
        Number of terms to allocate. n_term < 1 => complex_taylor.term pointer will be disassociated.

    save : bool, optional
        If True then save any old terms when complex_taylor is resized. Default is False.
    """

@overload
def init_coord(orb: _pybmad.CoordStruct, vec: Sequence[float], ele: _pybmad.EleStruct | None = None, element_end: int | None = None, particle: int | None = None, direction: int | None = None, E_photon: float | None = None, t_offset: float | None = None, shift_vec6: bool | None = None, spin: Sequence[float] | None = None, s_pos: float | None = None, random_on: bool | None = None) -> None:
    """
    Routine to initialize a coord_struct.

    This routine is an overloaded name for:
      Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
      Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
      Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

    Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice),
    or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6)
    corresponding to the beginning_ele's value of ele%value(p0c_start$).

    Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead
    of the standard which is to set orb%t from orb%vec(5).

    For photons:
      orb%vec(5) is set depending upon where the photon is relative to the element.
      If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle
          except if direction is set.

    Parameters
    ----------
    orb : CoordStruct
        Input orbit

    vec : 1D array of float (shape: 6)
        Coordinate vector. If not present then taken to be zero.

    ele : EleStruct, optional
        Particle is initialized to start at element_end of this ele.

    element_end : int, optional
        upstream_end$, downstream_end$, inside$, or start_end$. Must be present if ele argument is present.
        start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1. Default is
        upstream_end$. Note: If ele is the beginning element (index zero), the setting of element_end will not
        matter.

    particle : int, optional
        Particle type (electron$, etc.). If particle = not_set$ and orb_in is present, use orb_in.species instead.

    direction : int, optional
        +1 -> moving downstream +s direciton, -1 -> moving upstream. 0 -> Ignore. Default is to not change
        orb.direction except for photons which get set according to orb.vec(6).

    E_photon : float, optional
        Photon energy if particle is a photon. Ignored otherwise.

    t_offset : float, optional
        Offset of the reference time. This is non-zero when there are multiple bunches and the reference time for
        a particular particle is pegged to the time of the center of the bunch.

    shift_vec6 : bool, optional
        If present and False, prevent the shift of orb.vec(6).

    spin : 1D array of float (shape: 3), optional
        Particle spin. Taken to be zero if not present.

    s_pos : float, optional
        Particle s-position. Only relavent if element_end = inside$.

    random_on : bool, optional
        Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
        photon coords using a random number generator. If False, the photon coords will be centered within the
        distribution specified in the photon_init ele.
    """

@overload
def init_coord(orb_in: _pybmad.CoordStruct, ele: _pybmad.EleStruct | None = None, element_end: int | None = None, particle: int | None = None, direction: int | None = None, E_photon: float | None = None, t_offset: float | None = None, shift_vec6: bool | None = None, spin: Sequence[float] | None = None, s_pos: float | None = None, random_on: bool | None = None) -> _pybmad.CoordStruct:
    """
    Routine to initialize a coord_struct.

    This routine is an overloaded name for:
      Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
      Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
      Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

    Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice),
    or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6)
    corresponding to the beginning_ele's value of ele%value(p0c_start$).

    Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead
    of the standard which is to set orb%t from orb%vec(5).

    For photons:
      orb%vec(5) is set depending upon where the photon is relative to the element.
      If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle
          except if direction is set.

    Parameters
    ----------
    orb_in : CoordStruct
        Input orbit

    ele : EleStruct, optional
        Particle is initialized to start at element_end of this ele.

    element_end : int, optional
        upstream_end$, downstream_end$, inside$, or start_end$. Must be present if ele argument is present.
        start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1. Default is
        upstream_end$. Note: If ele is the beginning element (index zero), the setting of element_end will not
        matter.

    particle : int, optional
        Particle type (electron$, etc.). If particle = not_set$ and orb_in is present, use orb_in.species instead.

    direction : int, optional
        +1 -> moving downstream +s direciton, -1 -> moving upstream. 0 -> Ignore. Default is to not change
        orb.direction except for photons which get set according to orb.vec(6).

    E_photon : float, optional
        Photon energy if particle is a photon. Ignored otherwise.

    t_offset : float, optional
        Offset of the reference time. This is non-zero when there are multiple bunches and the reference time for
        a particular particle is pegged to the time of the center of the bunch.

    shift_vec6 : bool, optional
        If present and False, prevent the shift of orb.vec(6).

    spin : 1D array of float (shape: 3), optional
        Particle spin. Taken to be zero if not present.

    s_pos : float, optional
        Particle s-position. Only relavent if element_end = inside$.

    random_on : bool, optional
        Default is True. Used only for photons being initalized with a photon_init element. If True, vary the
        photon coords using a random number generator. If False, the photon coords will be centered within the
        distribution specified in the photon_init ele.

    Returns
    -------
    orb_out : CoordStruct
        Initialized coordinate
    """

@overload
def init_coord(orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct | None = None, element_end: int | None = None, particle: int | None = None, direction: int | None = None, E_photon: float | None = None, t_offset: float | None = None, shift_vec6: bool | None = None, spin: Sequence[float] | None = None) -> None:
    """
    Routine to initialize a coord_struct.

    This routine is an overloaded name for:
      Subroutine init_coord1 (orb, vec, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
      Subroutine init_coord2 (orb, orb_in, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)
      Subroutine init_coord3 (orb, ele, element_end, particle, direction, E_photon, t_offset, shift_vec6, spin, s_pos, random_on)

    Note: Unless shift_vec6 is set to False, if ele is a beginning_ele (IE, the element at the beginning of the lattice),
    or e_gun, orb%vec(6) is shifted so that a particle with orb%vec(6) = 0 will end up with a value of orb%vec(6)
    corresponding to the beginning_ele's value of ele%value(p0c_start$).

    Note: For non-photons, if orb_in%vec(5) is set to real_garbage$, orb_in%t will be used to set orb%vec(5) instead
    of the standard which is to set orb%t from orb%vec(5).

    For photons:
      orb%vec(5) is set depending upon where the photon is relative to the element.
      If orb is a photon, and orb_in is not a photon, photon is launched in same direciton as particle
          except if direction is set.

    Parameters
    ----------
    orb : CoordStruct
        Input orbit

    ele : EleStruct, optional
        Particle is initialized to start at element_end of this ele.

    element_end : int, optional
        upstream_end$, downstream_end$, inside$, or start_end$. Must be present if ele argument is present.
        start_end$ -> upstream_end$ if dir = 1 and start_end$ -> downstream_end$ if dir = -1. Default is
        upstream_end$. Note: If ele is the beginning element (index zero), the setting of element_end will not
        matter.

    particle : int, optional
        Particle type (electron$, etc.). If particle = not_set$ and orb_in is present, use orb_in.species instead.

    direction : int, optional
        +1 -> moving downstream +s direciton, -1 -> moving upstream. 0 -> Ignore. Default is to not change
        orb.direction except for photons which get set according to orb.vec(6).

    E_photon : float, optional
        Photon energy if particle is a photon. Ignored otherwise.

    t_offset : float, optional
        Offset of the reference time. This is non-zero when there are multiple bunches and the reference time for
        a particular particle is pegged to the time of the center of the bunch.

    shift_vec6 : bool, optional
        If present and False, prevent the shift of orb.vec(6).

    spin : 1D array of float (shape: 3), optional
        Particle spin. Taken to be zero if not present.
    """

def init_custom(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine init_custom

    Parameters
    ----------
    lat : LatStruct
    """

def init_ele(key: int | None = None, sub_key: int | None = None, ix_ele: int | None = None, branch: _pybmad.BranchStruct | None = None) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine init_ele

    Parameters
    ----------
    key : int, optional
        Key to initialize to. EG: quadrupole$, etc.

    sub_key : int, optional
        Sub-key to initialize to.

    ix_ele : int, optional
        ix_ele index to initalize to. Default = -1.

    branch : BranchStruct, optional
        Branch to point ele.branch and ele.ix_branch to. Otherwise ele.branch is nullified and ele.ix_branch = 0.

    Returns
    -------
    ele : EleStruct
        Initialized element.
    """

def init_fringe_info(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct | None = None, leng_sign: int | None = None) -> _pybmad.FringeFieldInfoStruct:
    """
    Wrapper for Fortran routine init_fringe_info

    Parameters
    ----------
    ele : EleStruct
        Lattice element associated with fringe_info.

    orbit : CoordStruct, optional
        Particle position. Must be present for a full init. If not full init only fringe_info.has_fringe will be
        set.

    leng_sign : int, optional
        Is element length positive (+1) or negative (-1)? Must be present if orbit is present.

    Returns
    -------
    fringe_info : FringeFieldInfoStruct
        Fringe information.
    """

def init_gg_taylor_series(gg_taylor: _pybmad.GgTaylorStruct, n_term: int, save_old: bool | None = None) -> None:
    """
    Subroutine to initialize a Bmad gg_taylor series (6 of these series make
    a gg_taylor map). Note: This routine does not zero the structure. The calling
    routine is responsible for setting all values.

    Parameters
    ----------
    gg_taylor : GgTaylorStruct
        Old structure.
        This parameter is an input/output and is modified in-place.
        As an output, gg_taylor: Initalized structure.

    n_term : int
        Number of terms to allocate. n_term < 0 => gg_taylor.term pointer will be disassociated.

    save_old : bool, optional
        If True then save any old terms when gg_taylor is resized. Default is False.
    """

def init_lat(n: int | None = None, init_beginning_ele: bool | None = None) -> _pybmad.LatStruct:
    """
    Wrapper for Fortran routine init_lat

    Parameters
    ----------
    n : int, optional
        Upper bound lat.ele(0:) array is initialized to. Default is 10.

    init_beginning_ele : bool, optional
        Init lat.ele(0)? Default is False.

    Returns
    -------
    lat : LatStruct
        Initialized lat.
    """

def init_multipole_cache(ele: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine init_multipole_cache

    Parameters
    ----------
    ele : EleStruct
        Element to init
        This parameter is an input/output and is modified in-place.
        As an output, ele: Initalized element.
    """

def init_photon_from_a_photon_init_ele(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct, random_on: bool | None = None) -> None:
    """
    Wrapper for Fortran routine init_photon_from_a_photon_init_ele

    Parameters
    ----------
    ele : EleStruct

    param : LatParamStruct

    orbit : CoordStruct

    random_on : bool, optional
    """

class InitPhotonIntegProb:
    """init_photon_integ_prob return type"""

    @property
    def E_photon(self) -> float: ...

    @property
    def integ_prob(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def init_photon_integ_prob(gamma: float, g: float, E_min: float, E_max: float, vert_angle_min: float | None = None, vert_angle_max: float | None = None, vert_angle_symmetric: bool | None = None, energy_integ_prob: float | None = None) -> InitPhotonIntegProb:
    """
    vert_angle_max, vert_angle_symmetric, energy_integ_prob, E_photon) result (integ_prob)

    Routine to calcuate the integrated probability of emitting a photon in a given vertical angle range
    and in a given energy range

    Parameters
    ----------
    gamma : float
        Gamma factor of charged particle emitting photon.

    g : float
        1/rho bending strength.

    E_min : float
        Minimum photon energy.

    E_max : float
        Maximum photon energy.

    vert_angle_min : float, optional
        Lower bound of vertical angle range.

    vert_angle_max : float, optional
        Upper bound of vertical angle range.

    vert_angle_symmetric : bool, optional
        Use two symmetric ranges [-vert_angle_max, -vert_angle_min] and [vert_angle_min, vert_angle_max] instead
        of just [vert_angle_min, vert_angle_max]?

    energy_integ_prob : float, optional
        If present, E_photon will be set to the photon energy such that the integrated probability of generating a
        photon in the given angle and energy range in the interval [E_min, E_photon] is energy_integ_prob. That
        is, energy_integ_prob = 0 => E_photon = E_min and energy_integ_prob = 1 => E_photon = E_max.

    Returns
    -------
    integ_prob : float
        Integrated probablility of emitting a photon in given angle and energy range.

    E_photon : float, optional
        See energy_integ_prob. E_photon must be present if energy_integ_prob is.
    """

def init_spin_distribution(beam_init: _pybmad.BeamInitStruct, ele: _pybmad.EleStruct) -> _pybmad.BunchStruct:
    """
    Initializes a spin distribution according to beam_init%spin.

    Parameters
    ----------
    beam_init : BeamInitStruct
        Initialization parameters

    Returns
    -------
    bunch : BunchStruct
        Bunch of particles. .particle(:).spin
    """

def init_surface_segment(phot: _pybmad.PhotonElementStruct, ix_pt: int, iy_pt: int) -> None:
    """
    Routine to init the componentes in ele%photon%segmented%pt(ix_pt,iy_pt) for use with segmented surface calculations.

    Parameters
    ----------
    phot : PhotonElementStruct
        Surface structure.

    ix_pt : int
        index of grid point to init.

    iy_pt : int
        index of grid point to init.
    """

def init_taylor_series(bmad_taylor: _pybmad.TaylorStruct, n_term: int, save_old: bool | None = None) -> None:
    """
    Wrapper for Fortran routine init_taylor_series

    Parameters
    ----------
    bmad_taylor : TaylorStruct
        Old structure.
        This parameter is an input/output and is modified in-place.
        As an output, bmad_taylor: Initalized structure.

    n_term : int
        Number of terms to allocate. n_term < 0 => bmad_taylor.term pointer will be disassociated.

    save_old : bool, optional
        If True then save any old terms and ref orbit when bmad_taylor is resized. If False zero the ref orbit.
        Default is False.
    """

def init_wake(n_sr_long: int, n_sr_trans: int, n_sr_z: int, n_lr_mode: int, always_allocate: bool | None = None) -> _pybmad.WakeStruct | None:
    """
    Wrapper for Fortran routine init_wake

    Parameters
    ----------
    n_sr_long : int
        Number of terms: wake.sr.long.

    n_sr_trans : int
        Number of terms: wake.sr.trans.

    n_sr_z : int
        Number of terms: wake.sr.z.

    n_lr_mode : int
        Number of terms: wake.lr.mode.

    always_allocate : bool, optional
        If present and True then allways allocate wake even if n_lr_mode, etc. are all 0. Default is False.

    Returns
    -------
    wake : WakeStruct, optional
        Initialized structure.
    """

def insert_element(lat: _pybmad.LatStruct, insert_ele: _pybmad.EleStruct, ix_ele: int, ix_branch: int | None = None, orbit: _pybmad.CoordStructAlloc1D | None = None) -> None:
    """
    Wrapper for Fortran routine insert_element

    Parameters
    ----------
    lat : LatStruct
        lattice that will be modified
        This parameter is an input/output and is modified in-place.
        As an output, lat: lattice with new element inserted

    insert_ele : EleStruct
        element to insert into the lat

    ix_ele : int
        branch.ele(:) index where the new element is inserted.

    ix_branch : int, optional
        : branch index for the insertion. Default = 0.

    orbit : 1D array of CoordStruct, optional
        orbit array to enlarge.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Enlarged orbit array.
    """

def integrand_base(t: float, args: _pybmad.RealArray1D) -> float:
    """
    This vectorized private function is the integrand in equation 31 of Piwinski's paper.

    This intetegrand has a sharp exponential decay, and so a change of variables from t to y where t=exp(y)
    is applied.  This COV makes the integrand more evenly distributed over the domain of integration,
    which makes it easier for qtrap to integrate.

    The change of variables is done using integrand_base_cov, which is then integrated
    using qtrap.

    Parameters
    ----------
    t : float
        Array of reals over which to evaluate the integrand.
    """

def integrate_psi(bound: float, p0: float, args: Sequence[float]) -> float:
    """
    Integrate psi(t) from -bound to +bound.  The integration is done in two parts.  First from 0 to -bound, then from
    0 to +bound.

    Parameters
    ----------
    bound : float
        integration bound

    p0 : float
        psi(0).  Boundary condition.

    args : 1D array of float (shape: 1:8)
        Parameters and constants of DEQ.  See psi_prime comments for details.

    Returns
    -------
    result : float
        Integral of psi from -bound to +bound.
    """

def integrated_mats(eles: _pybmad.EleStructArray1D, coos: _pybmad.CoordStructArray1D, Lambda: Sequence[Sequence[complex]], Theta: Sequence[Sequence[complex]], Iota: Sequence[Sequence[complex]], mode: _pybmad.NormalModesStruct) -> None:
    """No docstring available."""

@overload
def integration_timer(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, start: _pybmad.CoordStruct, orb_max: _pybmad.CoordStruct, tol: float) -> None:
    """
    Wrapper for Fortran routine integration_timer_ele

    Parameters
    ----------
    ele : EleStruct

    param : LatParamStruct

    start : CoordStruct

    orb_max : CoordStruct

    tol : float
    """

@overload
def integration_timer(a_fibre: _pybmad.Fibre, orbit: Sequence[float], orbit_max: Sequence[float], tol_dp: float) -> None:
    """
    Wrapper for Fortran routine integration_timer_fibre

    Parameters
    ----------
    a_fibre : Fibre

    orbit : 1D array of float (shape: 6)

    orbit_max : 1D array of float (shape: 6)

    tol_dp : float
    """

class InterpolateField:
    """interpolate_field return type"""

    @property
    def E(self) -> list[float]: ...

    @property
    def B(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def interpolate_field(x: float, y: float, z: float, mesh3d: _pybmad.Mesh3DStruct) -> InterpolateField:
    """
    Interpolate field on mesh

    Parameters
    ----------
    x : float
        coordinates to interpolate

    y : float
        coordinates to interpolate

    z : float
        coordinates to interpolate

    mesh3d : Mesh3dStruct
        contains efield, bfield

    Returns
    -------
    E : 1D array of float (shape: 3), optional
        interpolated electric field at x, y, z

    B : 1D array of float (shape: 3), optional
        interpolated magnetic field at x, y, z
    """

def ion_kick(orbit: _pybmad.CoordStruct, r_beam: Sequence[float], n_beam_part: float, a_twiss: _pybmad.TwissStruct, b_twiss: _pybmad.TwissStruct, sig_ee: float) -> list[float]:
    """
    Wrapper for Fortran routine ion_kick

    Parameters
    ----------
    orbit : CoordStruct
        Ion position.

    r_beam : 1D array of float (shape: 2)
        Beam (x, y) position.

    n_beam_part : float
        Number of beam particles.

    a_twiss : TwissStruct
        Horizontal like beam twiss parameters.

    b_twiss : TwissStruct
        vertical like beam twiss parameters.

    sig_ee : float
        Sigma_E/E beam energy spread.

    Returns
    -------
    kick : 1D array of float (shape: 3)
        (x, y, s) kick in m/sec.
    """

def is_attribute(ix_attrib: int, which: int) -> bool:
    """
    Routine to determine if an attribute index corresponds to a control variable for overlys/groups.

    Parameters
    ----------
    ix_attrib : int
        Attribute index.

    which : int
        control_var$, old_control_var$, all_control_var$, multipole$, elec_multipole$

    Returns
    -------
    is_attrib : bool
        True if a control variable
    """

def key_name_to_key_index(key_str: str, abbrev_allowed: bool | None = None) -> int:
    """
    Wrapper for Fortran routine key_name_to_key_index

    Parameters
    ----------
    key_str : str
        Name of the key. Result is case insensitive.

    abbrev_allowed : bool, optional
        Abbreviations (eg: "quad") allowed? Default is False. At least 3 characters are needed (except for
        rfcavity elements) if True.

    Returns
    -------
    key_index : int
        Index of the key. Set to -1 if key_name not recognized.
    """

class KickVectorCalc:
    """kick_vector_calc return type"""

    @property
    def dr_ds(self) -> list[float]: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def kick_vector_calc(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, s_body: float, orbit: _pybmad.CoordStruct, print_err: bool | None = None) -> KickVectorCalc:
    """
    Subroutine to calculate the dr/ds "kick vector" where
        r = [x, p_x, y, p_y, z, p_z, t, spin_x,y,z]

    Remember: In order to simplify the calculation, in the body of any element, P0 is taken to be
    the P0 at the exit end of the element.

      dr(1)/ds = dx/ds = dx/dt * dt/ds
      where:
        dx/dt = v_x = p_x / (1 + p_z)
        dt/ds = (1 + g*x) / v_s
        g = 1/rho, rho = bending radius (nonzero only in a dipole)

      dr(2)/ds = dp_x/ds = dP_x/dt * dt/ds / P0 + g_x * P_z
      where:
        dP_x/dt = EM_Force_x
        g_x = bending in x-plane.

      dr(3)/ds = dy/ds = dy/dt * dt/ds
      where:
        dy/dt = v_x

      dr(4)/ds = dp_y/ds = dP_y/dt * ds/dt / P0 + g_y * P_z
      where:
        dP_y/dt = EM_Force_y
        g_y = bending in y-plane.

      NOTE: dr(5)/ds IS IGNORED WHEN CALCULATING Z. SEE TRANSFER_THIS_ORBIT ABOVE.
      dr(5)/ds = dz/ds = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * c_light * [t(ref) - t]
                       = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * vec(5) / beta
      where:
        dt/ds(ref) = 1 / beta(ref)

      dr(6)/ds = dp_z/ds = d(EM_Force dot v_hat) * dt/ds / P0
      where:
         v_hat = velocity normalized to 1.

      dr(7)/ds = dt/ds

      dr(8:10)/ds = Spin omega vector

      dr(11)/ds = dt_ref/ds

    Parameters
    ----------
    ele : EleStruct
        Element being tracked thorugh.

    param : LatParamStruct
        Lattice parameters.

    orbit : CoordStruct
        Position of particle.

    Returns
    -------
    dr_ds : 1D array of float (shape: 11)
        Kick vector.

    err : bool
        Set True if there is an error.
    """

def kill_complex_taylor(complex_taylor: _pybmad.ComplexTaylorStructArray1D) -> None:
    """
    Subroutine to deallocate a Bmad complex_taylor map.

    Parameters
    ----------
    complex_taylor : 1D array of ComplexTaylorStruct
        complex_taylor to be deallocated. It is OK if complex_taylor has already been deallocated.
        This parameter is an input/output and is modified in-place.
        As an output, complex_taylor: deallocated complex_taylor structure.
    """

def kill_ptc_layouts(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine kill_ptc_layouts

    Parameters
    ----------
    lat : LatStruct
        Bmad lattice with associated layouts.
    """

def kill_taylor(bmad_taylor: _pybmad.TaylorStructArray1D) -> None:
    """
    Wrapper for Fortran routine kill_taylor

    Parameters
    ----------
    bmad_taylor : 1D array of TaylorStruct
        Taylor to be deallocated.
        This parameter is an input/output and is modified in-place.
        As an output, bmad_taylor: deallocated Taylor structure.
    """

def kind_name(this_kind: int) -> str:
    """
    function to return the name of a PTC kind.

    Parameters
    ----------
    this_kind : int
        PTC kind

    Returns
    -------
    kind_str : str
        String representation
    """

class KnotInterpolate:
    """knot_interpolate return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def y_pt(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def knot_interpolate(x_knot: _pybmad.RealArray1D, y_knot: _pybmad.RealArray1D, x_pt: float, interpolation: int) -> KnotInterpolate:
    """
    Wrapper for Fortran routine knot_interpolate

    Parameters
    ----------
    x_knot : 1D array of float
        Knot x-values.

    y_knot : 1D array of float
        Knot y-values.

    x_pt : float
        Point to evaluate at.

    interpolation : int
        Interpolation type. cubic$ or linear$.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.

    y_pt : float
        Interpolated y-value.
    """

def knots_to_string(x_knot: _pybmad.RealArray1D, y_knot: _pybmad.RealArray1D) -> str:
    """
    Wrapper for Fortran routine knots_to_string

    Parameters
    ----------
    x_knot : 1D array of float

    y_knot : 1D array of float

    Returns
    -------
    str : str
    """

def lafun(x: float, y: float, z: float) -> float:
    """
    Wrapper for Fortran routine lafun

    Parameters
    ----------
    x : float

    y : float

    z : float

    Returns
    -------
    res : float
    """

def lat_compute_ref_energy_and_time(lat: _pybmad.LatStruct) -> bool:
    """
    Wrapper for Fortran routine lat_compute_ref_energy_and_time

    Parameters
    ----------
    lat : LatStruct
        Input lattice.

    Returns
    -------
    err_flag : bool
        Set true if there is an error. False otherwise.
    """

class LatEleLocator:
    """lat_ele_locator return type"""

    @property
    def err(self) -> bool: ...

    @property
    def n_loc(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def lat_ele_locator(loc_str: str, lat: _pybmad.LatStruct, eles: _pybmad.ElePointerStructAlloc1D, n_loc: int, above_ubound_is_err: bool | None = None, ix_dflt_branch: int | None = None, order_by_index: bool | None = None, append_eles: bool | None = None) -> LatEleLocator:
    """
    Wrapper for Fortran routine lat_ele_locator

    Parameters
    ----------
    loc_str : str
        Element names or indexes. May be lower case.

    lat : LatStruct
        Lattice to search through.

    eles : 1D array of ElePointerStruct
        If append_eles is True, save existing elements.
        This parameter is an input/output and is modified in-place.
        As an output, eles: Array of matching elements.

    n_loc : int
        Number of existing elements. Used if append_eles is True.
        This parameter is an input/output and is modified in-place.
        As an output, n_loc: Number of locations found.

    above_ubound_is_err : bool, optional
        Default is True. If the upper bound "e2" on an "e1:e2" range construct is an integer and above the maximum
        element index then treat this as an error? If False, treat e2 as the maximum element index.

    ix_dflt_branch : int, optional
        If present and not -1 then restrict search to specified branch. If not present or -1: Search all branches.
        Exception: For elements specified using an integer index (EG: "43"), if ix_dflt_branch is not present or
        -1 use branch 0.

    order_by_index : bool, optional
        False is default. If True, order a component of loc_str like "quad::*" by element index instead of
        longitudinal s-position. Index ordering and s-position ordering are different when there are super lords
        and super slaves.

    append_eles : bool, optional
        Default is False. If True, found elements are appended to eles(:) array.

    Returns
    -------
    n_loc : int
        Number of existing elements. Used if append_eles is True.
        This parameter is an input/output and is modified in-place.
        As an output, n_loc: Number of locations found.

    err : bool, optional
        Set True if there is a decode error. Note: Not finding any matching element is not an error.
    """

def lat_equal_lat(lat_out: _pybmad.LatStruct, lat_in: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine lat_equal_lat

    Parameters
    ----------
    lat_out : LatStruct

    lat_in : LatStruct
    """

def lat_geometry(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine lat_geometry

    Parameters
    ----------
    lat : LatStruct
        The lattice.
    """

def lat_make_mat6(lat: _pybmad.LatStruct, ix_ele: int | None = None, ref_orb: _pybmad.CoordStructArray1D | None = None, ix_branch: int | None = None) -> bool:
    """
    Wrapper for Fortran routine lat_make_mat6

    Parameters
    ----------
    lat : LatStruct
        Lat containing the elements.

    ix_ele : int, optional
        Index of the element. If not present or negative, the matrices for all elements will be calculated.

    ref_orb : 1D array of CoordStruct, optional
        Coordinates of the reference orbit around which the matrix is calculated. If not present then the
        referemce is taken to be the origin.

    ix_branch : int, optional
        Branch index. Default is 0 (main lattice). -1 => All branches/all elements (ref_orb & ix_ele will be
        ignored).

    Returns
    -------
    err_flag : bool, optional
        True if there is an error. False otherwise.
    """

def lat_sanity_check(lat: _pybmad.LatStruct) -> bool:
    """
    Wrapper for Fortran routine lat_sanity_check

    Parameters
    ----------
    lat : LatStruct
        Lattice to check

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def lat_to_ptc_layout(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine lat_to_ptc_layout

    Parameters
    ----------
    lat : LatStruct
        Input lattice
    """

def lat_vec_equal_lat_vec(lat1: _pybmad.LatStructArray1D, lat2: _pybmad.LatStructArray1D) -> None:
    """
    Wrapper for Fortran routine lat_vec_equal_lat_vec

    Parameters
    ----------
    lat1 : 1D array of LatStruct

    lat2 : 1D array of LatStruct
    """

def lattice_bookkeeper(lat: _pybmad.LatStruct) -> bool:
    """
    Wrapper for Fortran routine lattice_bookkeeper

    Parameters
    ----------
    lat : LatStruct
        Lattice needing bookkeeping.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with bookkeeping done.

    Returns
    -------
    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

def lcavity_rf_step_setup(ele: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine lcavity_rf_step_setup

    Parameters
    ----------
    ele : EleStruct
        Lcavity element.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with ele.rf properly setup.
    """

def linear_bend_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orb: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Subroutine to track through the edge field of an sbend.
    Apply only the first order kick, which is edge focusing.

    Parameters
    ----------
    ele : EleStruct
        SBend element.

    param : LatParamStruct
        Rel charge.

    particle_at : int
        first_track_edge$, or second_track_edge$,

    orb : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Coords after tracking.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the edge.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix including the edge.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

class LinearCoef:
    """linear_coef return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def coef(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def linear_coef(stack: _pybmad.ExpressionAtomStructArray1D) -> LinearCoef:
    """
    Routine to return the linear coefficient of a linear expression.

    Parameters
    ----------
    stack : 1D array of ExpressionAtomStruct
        Expression stack.

    Returns
    -------
    err_flag : bool
        Set True if the expression is not linear

    coef : float
        Linear coefficient.
    """

def linear_to_spin_taylor(q_map: Sequence[Sequence[float]]) -> _pybmad.TaylorStructArray1D:
    """
    Wrapper for Fortran routine linear_to_spin_taylor

    Parameters
    ----------
    q_map : 2D array of float (shape: 0:3, 0:6)
        Linear quaternion map.

    Returns
    -------
    spin_taylor : 1D array of TaylorStruct (shape: 0:3)
        Taylor map
    """

class LoadParseLine:
    """load_parse_line return type"""

    @property
    def end_of_file(self) -> bool: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def load_parse_line(action: str, ix_start: int) -> LoadParseLine:
    """
    Subroutine to load characters from the input file.
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Parameters
    ----------
    action : str
        'continue', 'new_command', or 'init'

    ix_start : int
        Index in bp_com.parse_line string where to append stuff.

    Returns
    -------
    end_of_file : bool
        End of file reached?

    err_flag : bool, optional
        Set True if there is an error. False otherwise
    """

def lord_edge_aligned(slave: _pybmad.EleStruct, slave_edge: int, lord: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine lord_edge_aligned

    Parameters
    ----------
    slave : EleStruct
        Slave element.

    slave_edge : int
        End under consideration: entrance_end$, exit_end$, in_between$, etc.

    lord : EleStruct
        Lord element.

    Returns
    -------
    is_aligned : bool
        True if a lord edge is aligned with the slave edge. If slave_edge is not entrance_end$ nor exit_end$ then
        is_aligned is False.
    """

def low_energy_z_correction(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, ds: float, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> float:
    """
    Wrapper for Fortran routine low_energy_z_correction

    Parameters
    ----------
    orbit : CoordStruct
        Position before correction

    ele : EleStruct
        Element being tracked through.

    ds : float
        Longitudinal distance traveled by reference particle.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the multipole.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including multipole.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    dz : float
        Change in z.
    """

def lsc_kick_params_calc(ele: _pybmad.EleStruct, csr: _pybmad.CsrStruct) -> None:
    """
    Routine to cache intermediate values needed for the lsc calculation.
    This routine is not for image currents.

    Parameters
    ----------
    ele : EleStruct
        Element to set up cache for.

    csr : CsrStruct
        This parameter is an input/output and is modified in-place.
        As an output, csr: Binned particle averages.
    """

def mad_add_offsets_and_multipoles(ele: _pybmad.EleStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to add in the effect of element offsets and/or multipoles
    on the 2nd order transport map for the element.

    Parameters
    ----------
    ele : EleStruct
        Drift element.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_concat_map2(map1: _pybmad.MadMapStruct, map2: _pybmad.MadMapStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to concatinate two 2nd order transport maps.
        map3 = map2(map1)
    The equivalent MAD-8 routine is: TMCAT1

    Parameters
    ----------
    map1 : MadMapStruct
        First map in the beam line.

    map2 : MadMapStruct
        Second map in the beam line.

    Returns
    -------
    map3 : MadMapStruct
        Concatinated map.
    """

def mad_drift(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for a drift space.
    The equivalent MAD-8 routine is: TMDRF

    Parameters
    ----------
    ele : EleStruct
        Drift element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_elsep(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for an electric separator.
    The equivalent MAD-8 routine is: TMSEP

    Parameters
    ----------
    ele : EleStruct
        Electric seperator element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_map_to_taylor(map: _pybmad.MadMapStruct, energy: _pybmad.MadEnergyStruct, taylor: _pybmad.TaylorStructArray1D) -> None:
    """
    Subroutine to convert a MAD order 2 map to a Bmad taylor map.
    The conversion will also convert between MAD's (t, dE) and Bmad's (beta*t, dP) coords.

    Parameters
    ----------
    map : MadMapStruct
        Order 2 map.

    energy : MadEnergyStruct
        Energy numbers.

    taylor : 1D array of TaylorStruct
        Taylor map.
    """

def mad_quadrupole(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for an quadrupole element.
    The equivalent MAD-8 routine is: TMSEXT

    Parameters
    ----------
    ele : EleStruct
        Quadrupole element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_rfcavity(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for an rfcavity element.
    The equivalent MAD-8 routine is: TMRF

    Parameters
    ----------
    ele : EleStruct
        Rfcavity element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_sbend(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for a sector bend element.
    The equivalent MAD-8 routine is: TMBEND

    Parameters
    ----------
    ele : EleStruct
        Sbend element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_sbend_body(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for the body of a sector dipole.
    The equivalent MAD-8 routine is: TMSECT

    Parameters
    ----------
    ele : EleStruct
        Solenoid element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_sbend_fringe(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct, into: bool) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for the fringe field of a dipole.
    The equivalent MAD-8 routine is: TMFRNG

    Parameters
    ----------
    ele : EleStruct
        Solenoid element.

    energy : MadEnergyStruct
        particle energy structure.

    into : bool
        If True then map is for particle entering a dipole

    Returns
    -------
    map : MadMapStruct
        Fringe dipole map.
    """

def mad_sextupole(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for an sextupole.
    The equivalent MAD-8 routine is: TMSEXT

    Parameters
    ----------
    ele : EleStruct
        Sextupole element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

def mad_solenoid(ele: _pybmad.EleStruct, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to make a transport map for an solenoid.
    The equivalent MAD-8 routine is: TMSEXT

    Parameters
    ----------
    ele : EleStruct
        Solenoid element.

    energy : MadEnergyStruct
        particle energy structure.

    Returns
    -------
    map : MadMapStruct
        Structure holding the transfer map.
    """

class MadTmfoc:
    """mad_tmfoc return type"""

    @property
    def c(self) -> float: ...

    @property
    def s(self) -> float: ...

    @property
    def d(self) -> float: ...

    @property
    def f(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def mad_tmfoc(el: float, sk1: float) -> MadTmfoc:
    """
    Subroutine to compute the linear focussing functions.
    The equivalent MAD-8 routine is: TMFOC

    Parameters
    ----------
    el : float
        Length.

    sk1 : float
        Quadrupole strength.

    Returns
    -------
    c : float
        Cosine-like function.             c(k,l)

    s : float
        Sine-like function.               s(k,l)

    d : float
        Dispersion function.              d(k,l)

    f : float
        Integral of dispersion function.  f(k,l)
    """

def mad_tmsymm(te: Sequence[Sequence[Sequence[float]]]) -> None:
    """
    subroutine to symmertrize the 2nd order map t.
    The equivalent MAD-8 routine is: tmsymm

    Parameters
    ----------
    te : 3D array of float (shape: 6,6,6)
        array to be symmertrized.
        This parameter is an input/output and is modified in-place.
        As an output, te: symmetrized array.
    """

def mad_tmtilt(map: _pybmad.MadMapStruct, tilt: float) -> None:
    """
    Subroutine to apply a tilt to a transport map.
    The equivalent MAD-8 routine is: TMTILT

    Parameters
    ----------
    map : MadMapStruct
        Unrotated transport map.
        This parameter is an input/output and is modified in-place.
        As an output, map: Rotated transport map.

    tilt : float
        Tilt
    """

def mad_track1(c0: _pybmad.CoordStruct, map: _pybmad.MadMapStruct) -> _pybmad.CoordStruct:
    """
    Subroutine to track through a 2nd order transfer map.
    The equivalent MAD-8 routine is: TMTRAK

    Parameters
    ----------
    c0 : CoordStruct
        Starting coords.

    map : MadMapStruct
        2nd order map.

    Returns
    -------
    c1 : CoordStruct
        Ending coords.
    """

def make_g2_mats(twiss: _pybmad.TwissStruct, g2_mat: Sequence[Sequence[float]], g2_inv_mat: Sequence[Sequence[float]]) -> None:
    """
    Wrapper for Fortran routine make_g2_mats

    Parameters
    ----------
    twiss : TwissStruct
        Twiss parameters.

    g2_mat : 2D array of float (shape: 2,2)

    g2_inv_mat : 2D array of float (shape: 2,2)
    """

class MakeGMats:
    """make_g_mats return type"""

    @property
    def g_mat(self) -> list[list[float]]: ...

    @property
    def g_inv_mat(self) -> list[list[float]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_g_mats(ele: _pybmad.EleStruct) -> MakeGMats:
    """
    Wrapper for Fortran routine make_g_mats

    Parameters
    ----------
    ele : EleStruct
        Element

    Returns
    -------
    g_mat : 2D array of float (shape: 4,4)
        Normal mode to betaless coords

    g_inv_mat : 2D array of float (shape: 4,4)
        The inverse of G_MAT
    """

class MakeHvbp:
    """make_hvbp return type"""

    @property
    def B(self) -> list[list[float]]: ...

    @property
    def V(self) -> list[list[float]]: ...

    @property
    def H(self) -> list[list[float]]: ...

    @property
    def Vbar(self) -> list[list[float]] | None: ...

    @property
    def Hbar(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_hvbp(N: Sequence[Sequence[float]]) -> MakeHvbp:
    r"""
    Parameterizes the eigen-decomposition of the 6x6 transfer matrix into HVBP as defined in:
    "From the beam-envelop matrix to synchrotron-radiation integrals" by Ohmi, Hirata, and Oide.

    This routine takes N, which is usually made from make_N (also in this module), and decomposes
    it into H, V, B, and P.

    N is defined by:
    M = N.U.Inverse[N] where U is block diagonal and the blocks are 2x2 rotation matrices.
    and it is decomposed by this subroutine as,
    N = H.V.B.P
    P has the same free parameters as B
    B "Twiss matrix" has 6 free parameters (Twiss alphas and betas)
    B blocks have the form /     sqrt(beta)         0       \
                           \ -alpha/sqrt(beta) 1/sqrt(beta) /
    V "Teng matrix" has 4 free parameters (xy, xpy, ypx, and pxpy coupling)
    H "Dispersion matrix" has 8 free parameters (xz, xpz, pxz, pxpz, yz, ypz, pyz, pypz coupling)

    Parameters
    ----------
    N : 2D array of float (shape: 6,6)
        Matrix of eigenvectors prepared by make_N

    Returns
    -------
    B : 2D array of float (shape: 6,6)
        Block diagonal matrix of Twiss parameters

    V : 2D array of float (shape: 6,6)
        horizontal-vertical coupling information

    H : 2D array of float (shape: 6,6)
        horizontal-longitudinal and vertical-longitudinal coupling information

    Vbar : 2D array of float (shape: 6,6), optional
        mat_symp_conj(B).V.B

    Hbar : 2D array of float (shape: 6,6), optional
        mat_symp_conj(B).H.B
    """

def make_hybrid_lat(lat_in: _pybmad.LatStruct, use_taylor: bool | None = None, orb0_arr: _pybmad.CoordArrayStructArray1D | None = None) -> _pybmad.LatStruct:
    """
    Wrapper for Fortran routine make_hybrid_lat

    Parameters
    ----------
    lat_in : LatStruct
        Input lattice.

    use_taylor : bool, optional
        If present and True then the hybrid elements will have a taylor series instead of a simple linear matrix.
        If an element to be concatenated has a taylor series then this taylor series will be concatenated with the
        other elements in the hybrid element.

    orb0_arr : 1D array of CoordArrayStruct, optional
        Central orbit for taylor stuff. Each orb0_arr(i).orbit(:) holds the orbit for the i^th lattice branch

    Returns
    -------
    lat_out : LatStruct
        Lattice with hybrid elements. Note: Lat_out must not be the same actual argument as lat_in.
    """

class MakeMadMap:
    """make_mad_map return type"""

    @property
    def energy(self) -> _pybmad.MadEnergyStruct: ...

    @property
    def map(self) -> _pybmad.MadMapStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_mad_map(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> MakeMadMap:
    """
    Subroutine to make a 2nd order transport map a la MAD.

    Parameters
    ----------
    ele : EleStruct
        Element

    param : LatParamStruct
        particle id

    Returns
    -------
    energy : MadEnergyStruct
        Energy of the particle

    map : MadMapStruct
        Structure holding the transfer map.
    """

class MakeMat6:
    """make_mat6 return type"""

    @property
    def end_orb(self) -> _pybmad.CoordStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_mat6(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, start_orb: _pybmad.CoordStruct | None = None) -> MakeMat6:
    """
    Wrapper for Fortran routine make_mat6

    Parameters
    ----------
    ele : EleStruct
        Element holding the transfer matrix.

    param : LatParamStruct
        Lattice global parameters.

    start_orb : CoordStruct, optional
        Reference coordinates at the beginning of element. If not present, default is to use the zero orbit.

    Returns
    -------
    end_orb : CoordStruct, optional
        Reference coordinates at the end of element.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

class MakeMat6Bmad:
    """make_mat6_bmad return type"""

    @property
    def end_orb(self) -> _pybmad.CoordStruct: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_mat6_bmad(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, start_orb: _pybmad.CoordStruct) -> MakeMat6Bmad:
    """
    Wrapper for Fortran routine make_mat6_bmad

    Parameters
    ----------
    ele : EleStruct
        Element to track through.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with transfer matrix.

    param : LatParamStruct
        Parameters are needed for some elements.

    start_orb : CoordStruct
        Starting coords.

    Returns
    -------
    end_orb : CoordStruct
        Coordinates at the end of element.

    err : bool, optional
        Set True if there is an error. False otherwise.
    """

class MakeMat6BmadPhoton:
    """make_mat6_bmad_photon return type"""

    @property
    def end_orb(self) -> _pybmad.CoordStruct: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_mat6_bmad_photon(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, start_orb: _pybmad.CoordStruct) -> MakeMat6BmadPhoton:
    """
    Wrapper for Fortran routine make_mat6_bmad_photon

    Parameters
    ----------
    ele : EleStruct
        Element with transfer matrix
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with transfer matrix.

    param : LatParamStruct
        Parameters are needed for some elements.

    start_orb : CoordStruct
        Coordinates at the beginning of element.

    Returns
    -------
    end_orb : CoordStruct
        Coordinates at the end of element.

    err : bool, optional
        Set True if there is an error. False otherwise.
    """

def make_mat6_high_energy_space_charge(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> None:
    """
    Routine to add the ultra relativistic space charge kick to the element transfer matrix.
    The routine setup_space_charge_calc must be called
    initially before any tracking is done. This routine assumes a Gaussian
    bunch and is only valid with relativistic particles where the effect
    of the space charge is small.

    Parameters
    ----------
    ele : EleStruct
        Element tracked through.

    param : LatParamStruct
    """

def make_mat6_mad(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, c0: _pybmad.CoordStruct) -> _pybmad.CoordStruct:
    """
    Subroutine to make the 6x6 transfer matrix for an element from the
    2nd order MAD transport map. The map is stored in ele%taylor.
    If the map exists then it is simply used to calculate ele%mat6.
    If ele%taylor doesn't exist then calculate it.

    Parameters
    ----------
    ele : EleStruct
        Element with transfer matrix.

    param : LatParamStruct
        Lattice parameters.

    c0 : CoordStruct
        Coordinates at the beginning of element.

    Returns
    -------
    c1 : CoordStruct
        Coordinates at the end of element.
    """

def make_mat6_symp_lie_ptc(ele: _pybmad.EleStruct, start_orb: _pybmad.CoordStruct) -> _pybmad.CoordStruct:
    """
    Wrapper for Fortran routine make_mat6_symp_lie_ptc

    Parameters
    ----------
    ele : EleStruct
        Element with transfer matrix
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with transfer matrix.

    start_orb : CoordStruct
        Coordinates at the beginning of element.

    Returns
    -------
    end_orb : CoordStruct
        Coordinates at end of element.
    """

def make_mat6_taylor(ele: _pybmad.EleStruct, start_orb: _pybmad.CoordStruct, err_flag: bool | None = None) -> _pybmad.CoordStruct:
    """
    Wrapper for Fortran routine make_mat6_taylor

    Parameters
    ----------
    ele : EleStruct
        Element to track through.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with transfer matrix.

    start_orb : CoordStruct
        Starting coords.

    err_flag : bool, optional

    Returns
    -------
    end_orb : CoordStruct
        Coordinates at the end of element.
    """

class MakeMat6Tracking:
    """make_mat6_tracking return type"""

    @property
    def end_orb(self) -> _pybmad.CoordStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_mat6_tracking(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, start_orb: _pybmad.CoordStruct, spin_only: bool | None = None) -> MakeMat6Tracking:
    """
    Wrapper for Fortran routine make_mat6_tracking

    Parameters
    ----------
    ele : EleStruct
        Element with transfer matrix
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with transfer matrix.

    param : LatParamStruct
        Parameters are needed for some elements.

    start_orb : CoordStruct
        Coordinates at the beginning of element.

    spin_only : bool, optional
        Default False. If True, only calculate ele.spin_taylor.

    Returns
    -------
    end_orb : CoordStruct
        Coordinates at the end of element.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class MakeN:
    """make_n return type"""

    @property
    def N(self) -> list[list[float]]: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def tunes_out(self) -> list[float]: ...

    @property
    def U(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_n(t6: Sequence[Sequence[float]], abz_tunes: Sequence[float] | None = None) -> MakeN:
    """
    Given a 1-turn transfer matrix, this returns the matrix N.
    N converts between normal invarients and phases and canonical coordinates:
    X = N.J

    N is obtained from the Eigen decomposition of the 1-turn transfer matrix.
    It is obtained by applying certain normalizations to the matrix of Eigen vectors, then making
    the result real using Q.

    If abz_tunes is present, then the eigensystem is ordered by matching the tunes.
    If abz_tunes is not present, then the eigensystem is ordered by plane dominance.

    It is assumed that the synchrotron tune is less than pi.

    Parameters
    ----------
    t6 : 2D array of float (shape: 6,6)
        1-turn transfer matrix

    abz_tunes : 1D array of float (shape: 3), optional
        a-mode is abz_tunes(1), b-mode is abz_tunes(2), synch tune is abz_tunes(3)

    Returns
    -------
    N : 2D array of float (shape: 6,6)
        X = N.J

    err_flag : bool
        Set to true on error.  Often means Eigen decomposition failed.

    tunes_out : 1D array of float (shape: 3), optional
        Fractional tune (in radians) of the 3 normal modes of t6.

    U : 2D array of float (shape: 6,6), optional
        U = Inverse(N).t6.N.  Block diagonal matrix of 2x2 rotation matrices.
    """

class MakePbrh:
    """make_pbrh return type"""

    @property
    def P(self) -> list[list[complex]]: ...

    @property
    def Bp(self) -> list[list[complex]]: ...

    @property
    def R(self) -> list[list[complex]]: ...

    @property
    def H(self) -> list[list[complex]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_pbrh(M: Sequence[Sequence[float]], abz_tunes: Sequence[float]) -> MakePbrh:
    """
    Decomposes the 1-turn transfer matrix into normal mode twiss-like parameters,
    according to Sec. IIIB of Ohmi, Hirata, and Oide paper.

    Note:  The Twiss parameters generated by this function are identical to those delivered
           by mode3_mod.

    Parameters
    ----------
    M : 2D array of float (shape: 6,6)
        1-turn transfer matrix

    abz_tunes : 1D array of float (shape: 3)
        tunes for a,b, and c modes.  Used to identify which eigenvector is associated with which mode.

    Returns
    -------
    P : 2D array of complex (shape: 6,6)
        Eqn. 97.  Phase advances.

    Bp : 2D array of complex (shape: 6,6)
        Eqns. 89 & 101.  Beta functions.

    R : 2D array of complex (shape: 6,6)
        Eqn. 99.  Transverse coupling.

    H : 2D array of complex (shape: 6,6)
        Eqn. 100.  Longitudinal coupling.
    """

class MakeSmatFromAbc:
    """make_smat_from_abc return type"""

    @property
    def sigma_mat(self) -> list[list[float]]: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def Nout(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_smat_from_abc(t6: Sequence[Sequence[float]], mode: _pybmad.NormalModesStruct) -> MakeSmatFromAbc:
    """
    Given the 1-turn transfer matrix and a normal_modes_struct containing the normal mode
    emittances, this routine returns the beam envelop sigma matrix.

    sigma_mat = N.D.transpose(N)
    equivalent to: sigma_mat.S = N.D.mat_symp_conj(N)

    One way to populate mode%a%tune and mode%b%tune:
      mode%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi, twopi)
      mode%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi, twopi)

    Parameters
    ----------
    t6 : 2D array of float (shape: 6,6)
        1-turn transfer matrix

    mode : NormalModesStruct
        normal mode emittances

    Returns
    -------
    sigma_mat : 2D array of float (shape: 6,6)
        beam envelop sigma matrix

    err_flag : bool
        set to true if something goes wrong.  Usually means Eigen decomposition of the 1-turn matrix failed.

    Nout : 2D array of float (shape: 6,6), optional
        Contains the normalized eigenvectors that were used to make the sigma matrix.
    """

def make_unit_mad_map(map: _pybmad.MadMapStruct) -> None:
    """
    Subroutine to initialize a 2nd order transport map to unity.

    Parameters
    ----------
    map : MadMapStruct
        2nd order transport map.
        This parameter is an input/output and is modified in-place.
        As an output, map: Unity 2nd order map.
    """

def make_v(M: Sequence[Sequence[float]], V: Sequence[Sequence[complex]], abz_tunes: Sequence[float]) -> None:
    """
    For a one-turn transfer matrix M, this routine find the eigen matrix V.
    V is ordered such that the per turn phase advance of its column pairs agree with abz_tunes.
    It is normalized to be symplectic.
    """

class MakeVMats:
    """make_v_mats return type"""

    @property
    def v_mat(self) -> list[list[float]] | None: ...

    @property
    def v_inv_mat(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def make_v_mats(ele: _pybmad.EleStruct) -> MakeVMats:
    """
    Wrapper for Fortran routine make_v_mats

    Parameters
    ----------
    ele : EleStruct
        Element

    Returns
    -------
    v_mat : 2D array of float (shape: 4,4), optional
        Normal mode to X-Y coords transformation

    v_inv_mat : 2D array of float (shape: 4,4), optional
        X-Y coords to Normal mode transformation
    """

def makeup_control_slave(lat: _pybmad.LatStruct, slave: _pybmad.EleStruct, err_flag: bool) -> None:
    """This routine is not meant for general use."""

def makeup_group_lord(lat: _pybmad.LatStruct, lord: _pybmad.EleStruct, err_flag: bool) -> None:
    """
    Subroutine to calculate the attributes of group slave elements.
    This routine is private to bookkeeper_mod.
    """

def makeup_multipass_slave(lat: _pybmad.LatStruct, slave: _pybmad.EleStruct, err_flag: bool) -> None:
    """
    Subroutine to calcualte the attributes of multipass slave elements.
    This routine is not meant for guse.
    """

def makeup_super_slave(lat: _pybmad.LatStruct, slave: _pybmad.EleStruct, err_flag: bool) -> None:
    """
    Subroutine to calcualte the attributes of superposition slave elements.
    This routine is not meant for general use.
    """

def makeup_super_slave1(slave: _pybmad.EleStruct, lord: _pybmad.EleStruct, offset: float, param: _pybmad.LatParamStruct, include_upstream_end: bool, include_downstream_end: bool) -> bool:
    """
    Routine to construct a super_slave from a super_lord when the slave has only one lord.
    Note: Reference energy and times are not computed in this routine.

    Parameters
    ----------
    slave : EleStruct
        Slave element.
        This parameter is an input/output and is modified in-place.
        As an output, slave: Slave element with appropriate values set.

    lord : EleStruct
        Lord element.

    offset : float
        offset of entrance end of slave from entrance end of the lord.

    param : LatParamStruct
        lattice paramters.

    include_upstream_end : bool
        Slave contains the lord's entrance end?

    include_downstream_end : bool
        Slave contains the lord's exit end?

    Returns
    -------
    err_flag : bool
        Set true if there is an error. False otherwise.
    """

def map1_inverse(map1: _pybmad.SpinOrbitMap1Struct) -> _pybmad.SpinOrbitMap1Struct:
    """
    Wrapper for Fortran routine map1_inverse

    Parameters
    ----------
    map1 : SpinOrbitMap1Struct
        Input map.

    Returns
    -------
    inv_map1 : SpinOrbitMap1Struct
        Inverse map.
    """

def map1_make_unit() -> _pybmad.SpinOrbitMap1Struct:
    """
    Wrapper for Fortran routine map1_make_unit

    Returns
    -------
    map1 : SpinOrbitMap1Struct
        Unit map.
    """

def map1_times_map1(map2: _pybmad.SpinOrbitMap1Struct, map1: _pybmad.SpinOrbitMap1Struct) -> _pybmad.SpinOrbitMap1Struct:
    """
    Wrapper for Fortran routine map1_times_map1

    Parameters
    ----------
    map2 : SpinOrbitMap1Struct

    map1 : SpinOrbitMap1Struct

    Returns
    -------
    map_out : SpinOrbitMap1Struct
    """

def map_to_angle_coords(t_canon: _pybmad.TaylorStructArray1D) -> _pybmad.TaylorStructArray1D:
    """
    Wrapper for Fortran routine map_to_angle_coords

    Parameters
    ----------
    t_canon : 1D array of TaylorStruct (shape: 6)
        Taylor map in canonical coords.

    Returns
    -------
    t_angle : 1D array of TaylorStruct (shape: 6)
        Taylor map in angle coords.
    """

def mark_patch_regions(branch: _pybmad.BranchStruct) -> None:
    """
    Routine to mark which regions in a wall3d structure contain patch elements.
    This routine should be called by any routine that creates a beam chamber wall.

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch with .wall3d beam chamber wall.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Lattice branch with .wall3d.section(i).patch_in_region marked.
    """

def master_parameter_value(master_parameter: int, ele: _pybmad.EleStruct) -> float:
    """
    Wrapper for Fortran routine master_parameter_value

    Parameters
    ----------
    master_parameter : int
        Index of the master parameter.

    ele : EleStruct
        Element containing the fieldmap.

    Returns
    -------
    value : float
        Value of the master parameter.
    """

def mat4_multipole(knl: float, tilt: float, n: int, orbit: _pybmad.CoordStruct) -> list[list[float]]:
    """
    Wrapper for Fortran routine mat4_multipole

    Parameters
    ----------
    knl : float
        Strength of multipole

    tilt : float
        Tilt of multipole

    n : int

    orbit : CoordStruct
        coordinates of particle

    Returns
    -------
    kick_mat : 2D array of float (shape: 4,4)
        Kick matrix (Jacobian) at orbit.
    """

def mat6_add_offsets(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> None:
    """
    Wrapper for Fortran routine mat6_add_offsets

    Parameters
    ----------
    ele : EleStruct
        Element with given orientation.

    param : LatParamStruct
    """

def mat6_add_pitch(x_pitch_tot: float, y_pitch_tot: float, orientation: int, mat6: Sequence[Sequence[float]]) -> None:
    """
    Wrapper for Fortran routine mat6_add_pitch

    Parameters
    ----------
    x_pitch_tot : float
        Horizontal pitch

    y_pitch_tot : float
        Vertical pitch

    orientation : int
        Element longitudinal orientation. +1 or -1.

    mat6 : 2D array of float (shape: 6,6)
        1st order part of the transfer map (Jacobian).
        This parameter is an input/output and is modified in-place.
        As an output, mat6: 1st order xfer map with pitches.
    """

def mat6_to_complex_taylor(vec0: Sequence[complex], mat6: Sequence[Sequence[complex]]) -> _pybmad.ComplexTaylorStructArray1D:
    """
    Subroutine to form a first order complex_taylor map from the 6x6 transfer
    matrix and the 0th order transfer vector.

    Parameters
    ----------
    vec0 : 1D array of complex (shape: 6)
        0th order transfer vector.

    mat6 : 2D array of complex (shape: 6,6)
        6x6 transfer matrix.

    Returns
    -------
    complex_taylor : 1D array of ComplexTaylorStruct (shape: 6)
        first order complex_taylor map.
    """

class MatSympDecouple:
    """mat_symp_decouple return type"""

    @property
    def stat(self) -> int: ...

    @property
    def U(self) -> list[list[float]]: ...

    @property
    def V(self) -> list[list[float]]: ...

    @property
    def Ubar(self) -> list[list[float]]: ...

    @property
    def Vbar(self) -> list[list[float]]: ...

    @property
    def G(self) -> list[list[float]]: ...

    @property
    def twiss1(self) -> _pybmad.TwissStruct: ...

    @property
    def twiss2(self) -> _pybmad.TwissStruct: ...

    @property
    def gamma(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def mat_symp_decouple(t0: Sequence[Sequence[float]], type_out: bool) -> MatSympDecouple:
    """
    Wrapper for Fortran routine mat_symp_decouple

    Parameters
    ----------
    t0 : 2D array of float (shape: 4,4)
        Input matrix

    type_out : bool
        If .true. then an error message is typed out for a non ok$ STAT

    Returns
    -------
    stat : int
        status of results: ok$, in_stop_band$, or unstable$

    u : 2D array of float (shape: 4,4)
        See MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.

    v : 2D array of float (shape: 4,4)
        See MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.

    ubar : 2D array of float (shape: 4,4)
        See MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.

    vbar : 2D array of float (shape: 4,4)
        See MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.

    g : 2D array of float (shape: 4,4)
        See MGB CBN 85-2 and PPB/DLR PAC89 papers for more info.

    twiss1 : TwissStruct
        Twiss params for the "upper left" mode.

    twiss2 : TwissStruct
        Twiss params for the "lower right" mode.

    gamma : float
        gamma_c factor.
    """

class MatchEleToMat6:
    """match_ele_to_mat6 return type"""

    @property
    def mat6(self) -> list[list[float]]: ...

    @property
    def vec0(self) -> list[float]: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def match_ele_to_mat6(ele: _pybmad.EleStruct, start_orb: _pybmad.CoordStruct, include_delta_time: bool | None = None, set_trombone: bool | None = None) -> MatchEleToMat6:
    """
    Wrapper for Fortran routine match_ele_to_mat6

    Parameters
    ----------
    ele : EleStruct
        Match element.

    start_orb : CoordStruct
        Starting orbit.

    include_delta_time : bool, optional
        If False, ignore any finite ele.value(delta_time$). Default is True.

    set_trombone : bool, optional
        Default is False. If True, set the beginning and ending Twiss values in the element to create a phase
        trombone.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6)
        Transfer matrix (1st order part of xfer map).

    vec0 : 1D array of float (shape: 6)
        0th order part of the transfer map.

    err_flag : bool
        Set true if there is an error. False otherwise. Note: Currently err_flag is never set True.
    """

def mexp(x: float, m: int) -> float:
    """
    Wrapper for Fortran routine mexp

    Parameters
    ----------
    x : float
        Number.

    m : int
        Exponent.

    Returns
    -------
    this_exp : float
        Result.
    """

def mfft1(a: _pybmad.RealArray1D, b: _pybmad.RealArray1D, n: _pybmad.IntArray1D, ndim: int, isn: int) -> int:
    """
    Wrapper for Fortran routine mfft1

    Parameters
    ----------
    a : 1D array of float

    b : 1D array of float

    n : 1D array of int

    ndim : int

    isn : int

    Returns
    -------
    ierr : int
    """

def misalign_ptc_fibre(ele: _pybmad.EleStruct, use_offsets: bool, for_layout: bool) -> _pybmad.Fibre | None:
    """
    Routine to misalign a fibre associated with a Bmad element.

    Parameters
    ----------
    ele : EleStruct
        Bmad element with misalignments.

    use_offsets : bool
        Does ptc_fibre include element offsets, pitches and tilt? This argument is ignored if the element is a
        patch.

    for_layout : bool
        If True then fibre is being created as part of a layout as opposed to a stand-alone fibre

    Returns
    -------
    ptc_fibre : Fibre, optional
        PTC fibre element with misalignments.
    """

def momentum_compaction(branch: _pybmad.BranchStruct) -> float:
    """
    Wrapper for Fortran routine momentum_compaction

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch to calculate on.

    Returns
    -------
    mom_comp : float
        Momentum compaction.
    """

def mpxx1(ele: _pybmad.EleStruct, coulomb_log: float, n_part: float) -> _pybmad.IbsStruct:
    """
    Modified Piwinski, further modified to treat Coulomb Log
    in the same manner as Bjorken-Mtingwa, CIMP, Bane, Kubo & Oide, etc.
    This formula is derived in Section 2.8.4 of Michael Ehrlichman's Graduate Thesis.
    """

def mpzt1(ele: _pybmad.EleStruct, coulomb_log: float, n_part: float) -> _pybmad.IbsStruct:
    """
    Modified Piwinski with Zotter's integral.  This is Piwinski's original derivation,
    generalized to take the derivatives of the optics functions.  Also, Piwinski's
    original cumbersome triple integral is reaplaced by Zotter's single integral.  Zotter's
    integral is exact, and not an approximation.

    rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
    """

def multi_coulomb_log(ibs_sim_params: _pybmad.IbsSimParamStruct, ele: _pybmad.EleStruct, coulomb_log: float, n_part: float) -> None:
    """
    Calculates the value of the Coulomb log using various methods.

    ibs_sim_params%clog_to_use == 1   Classic coulomb log (pi/2 max scattering angle)
    ibs_sim_params%clog_to_use == 2   Integral based tail-cut prescribed by Raubenheimer.
    ibs_sim_params%clog_to_use == 3   Bane tail cut. 1 event/part/damping period.
    ibs_sim_params%clog_to_use == 4   Kubo and Oide tail cut. Used by CesrTA publications.
    """

class MultiTurnTrackingAnalysis:
    """multi_turn_tracking_analysis return type"""

    @property
    def track0(self) -> _pybmad.CoordStruct: ...

    @property
    def ele(self) -> _pybmad.EleStruct: ...

    @property
    def stable(self) -> bool: ...

    @property
    def growth_rate(self) -> float: ...

    @property
    def chi(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def multi_turn_tracking_analysis(track: _pybmad.CoordStructArray1D, i_dim: int) -> MultiTurnTrackingAnalysis:
    """
    Wrapper for Fortran routine multi_turn_tracking_analysis

    Parameters
    ----------
    track : 1D array of CoordStruct
        multi-turn tracking data to analyze. track(i) is the particle position at a given point in the lat on the
        i^th turn.

    i_dim : int
        number of dimensions used in the tracking: 2, or 4.

    Returns
    -------
    track0 : CoordStruct
        Closed orbit.

    ele : EleStruct
        structure holding the 1-turn matrix and Twiss parameters.

    stable : bool
        Is motion stable?

    growth_rate : float
        Unstable growth rate (= 0 if stable).

    chi : float
        How symplectic the computed 1-turn matrix is. See mat_symp_check for more details.

    err_flag : bool
        Set true if there is an error. False otherwise.
    """

def multilayer_type_to_multilayer_params(ele: _pybmad.EleStruct) -> bool:
    """
    Routine to set the multilayer parameters based upon the multilayer type.

    Multilayer types are of the form:
      "AAA:BBB"
    Where "AAA" is the atomic formula for the top layer crystal and "BBB" is the second layer atomic formula.

    Parameters
    ----------
    ele : EleStruct
        Multilayer element.

    Returns
    -------
    err_flag : bool
        Set True if multilayer type is unrecognized. False otherwise.
    """

def multipass_all_info(lat: _pybmad.LatStruct) -> _pybmad.MultipassAllInfoStruct:
    """
    Wrapper for Fortran routine multipass_all_info

    Parameters
    ----------
    lat : LatStruct
        Lattice

    Returns
    -------
    info : MultipassAllInfoStruct
        Multipass information.
    """

class MultipassChain:
    """multipass_chain return type"""

    @property
    def ix_pass(self) -> int: ...

    @property
    def n_links(self) -> int: ...

    @property
    def chain_ele(self) -> _pybmad.ElePointerStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def multipass_chain(ele: _pybmad.EleStruct, use_super_lord: bool | None = None) -> MultipassChain:
    """
    Wrapper for Fortran routine multipass_chain

    Parameters
    ----------
    ele : EleStruct
        Element in a multipass chain.

    use_super_lord : bool, optional
        If present and True and if ele is a super_slave, construct the chain_ele(:) array using the corresponding
        super_lords.

    Returns
    -------
    ix_pass : int
        Multipass pass number of the input element. Set to -1 if input element is not in a multipass section.

    n_links : int
        Number of times the physical element is passed through.

    chain_ele : 1D array of ElePointerStruct, optional
        pointers to the elements of the chain. Note: chain_ele(ix_pass).ele => ele
    """

def multipass_region_info(lat: _pybmad.LatStruct, mult_lat: _pybmad.MultipassRegionLatStruct, m_info: _pybmad.MultipassAllInfoStruct) -> None:
    """
    Wrapper for Fortran routine multipass_region_info

    Parameters
    ----------
    lat : LatStruct

    mult_lat : MultipassRegionLatStruct

    m_info : MultipassAllInfoStruct
    """

class Multipole1AbToKt:
    """multipole1_ab_to_kt return type"""

    @property
    def knl(self) -> float: ...

    @property
    def tn(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def multipole1_ab_to_kt(an: float, bn: float, n: int) -> Multipole1AbToKt:
    """
    Wrapper for Fortran routine multipole1_ab_to_kt

    Parameters
    ----------
    an : float
        Skew multipole component.

    bn : float
        Normal multipole component.

    n : int
        Order of multipole.

    Returns
    -------
    knl : float
        Multitude magnatude.

    tn : float
        Multipole angle.
    """

class Multipole1KtToAb:
    """multipole1_kt_to_ab return type"""

    @property
    def an(self) -> float: ...

    @property
    def bn(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def multipole1_kt_to_ab(knl: float, knsl: float, tn: float, n: int) -> Multipole1KtToAb:
    """
    Wrapper for Fortran routine multipole1_kt_to_ab

    Parameters
    ----------
    knl : float
        Normal multitude component.

    knsl : float
        Skew multitude component.

    tn : float
        Multipole angle.

    n : int
        Multipole order.

    Returns
    -------
    an : float
        Skew multipole component.

    bn : float
        Normal multipole component.
    """

def multipole_ab_to_kt(an: _pybmad.RealArray1D, bn: _pybmad.RealArray1D, knl: _pybmad.RealArray1D, tn: _pybmad.RealArray1D) -> None:
    """
    Wrapper for Fortran routine multipole_ab_to_kt

    Parameters
    ----------
    an : 1D array of float
        Skew multipole component.

    bn : 1D array of float
        Normal multipole component.

    knl : 1D array of float
        Multitude magnatude.

    tn : 1D array of float
        Multipole angle.
    """

class MultipoleEleToAb:
    """multipole_ele_to_ab return type"""

    @property
    def ix_pole_max(self) -> int: ...

    @property
    def a(self) -> list[float]: ...

    @property
    def b(self) -> list[float]: ...

    @property
    def b1(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def multipole_ele_to_ab(ele: _pybmad.EleStruct, use_ele_tilt: bool, pole_type: int | None = None, include_kicks: int | None = None, original: bool | None = None) -> MultipoleEleToAb:
    """
    Wrapper for Fortran routine multipole_ele_to_ab

    Parameters
    ----------
    ele : EleStruct
        Element.

    use_ele_tilt : bool
        If True then include ele.value(tilt_tot$) in calculations. use_ele_tilt is ignored in the case of
        multipole$ elements.

    pole_type : int, optional
        Type of multipole. magnetic$ (default) or electric$.

    include_kicks : int, optional
        Ignored for for pole_type == electric$ for non-elseparator elements. Possibilities are:

    original : bool, optional
        Default is false. If True, no scaling is applied.

    Returns
    -------
    ix_pole_max : int
        Index of largest nonzero a(:) or b(:) pole. Set to -1 if all multipoles are zero. ix_pole_max is set
        independent of a nonzero b1 (if present).

    a : 1D array of float (shape: 0:n_pole_maxx)
        Array of multipole values.

    b : 1D array of float (shape: 0:n_pole_maxx)
        Array of multipole values.

    b1 : float, optional
        If present, b1 is set to the value of the b(1) component of the b(:) array and b(1) is set to zero. Also
        ix_pole_max is ajusted as needed. This is used by routines that want to handle b(1) in a special way in
        tracking.
    """

def multipole_ele_to_kt(ele: _pybmad.EleStruct, use_ele_tilt: bool, knl: _pybmad.RealArray1D, tilt: _pybmad.RealArray1D, pole_type: int | None = None, include_kicks: int | None = None) -> int:
    """
    Wrapper for Fortran routine multipole_ele_to_kt

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    use_ele_tilt : bool
        If True then include ele.value(tilt_tot$) in calculations. use_ele_tilt is ignored in the case of
        multipole$ elements.

    knl : 1D array of float
        Vector of strengths, MAD units.

    tilt : 1D array of float
        Vector of tilts.

    pole_type : int, optional
        Type of multipole. magnetic$ (default) or electric$.

    include_kicks : int, optional
        Possibilities are:

    Returns
    -------
    ix_pole_max : int
        Index of largest nonzero pole.
    """

def multipole_init(who: int, zero: bool | None = None) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine multipole_init

    Parameters
    ----------
    who : int
        electric$, magnetic$, or all$

    zero : bool, optional
        If present and True then zero the arrays even if they already exist when this routine is called. Default
        is False which means that if the arrays already exist then this routine will do nothing.

    Returns
    -------
    ele : EleStruct
        Element holding the multipoles.
    """

def multipole_kick(knl: float, tilt: float, n: int, ref_species: int, ele_orientation: int, coord: _pybmad.CoordStruct, pole_type: int | None = None, ref_orb_offset: bool | None = None) -> None:
    """
    Subroutine to put in the kick due to a multipole.
    Note: The kick for an electric multipole does not include any energy change.

    Parameters
    ----------
    knl : float
        Multipole integrated strength.

    tilt : float
        Multipole tilt.

    n : int
        Multipole order.

    ref_species : int
        Reference species.

    ele_orientation : int
        Element orientation +1 = normal, -1 = reversed.

    coord : CoordStruct
        Particle position and direction of travel.

    pole_type : int, optional
        Type of multipole. magnetic$ (default) or electric$.

    ref_orb_offset : bool, optional
        If True and n = 0 then use the MAD convention and model the multipole as a zero length bend with bending
        angle knl. Default is False.
    """

def multipole_kick_mat(knl: _pybmad.RealArray1D, tilt: _pybmad.RealArray1D, ref_species: int, ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, factor: float) -> list[list[float]]:
    """
    Wrapper for Fortran routine multipole_kick_mat

    Parameters
    ----------
    knl : 1D array of float
        Strength of multipoles

    tilt : 1D array of float
        Tilt of multipoles

    ref_species : int
        Reference species.

    ele : EleStruct
        Lattice element containing multipoles.

    orbit : CoordStruct
        coordinates of particle around which the multipole kick matrix is computed.

    factor : float
        Factor to scale knl by.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6)
        matrix with kick values at mat6(2:4:2, 1:3:2). The rest of the matrix is untouched.
    """

def multipole_kicks(knl: _pybmad.RealArray1D, tilt: _pybmad.RealArray1D, ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, pole_type: int | None = None, ref_orb_offset: bool | None = None) -> None:
    """
    Subroutine to put in the kick due to a multipole element.
    Also see the ab_multipole_kicks routine.

    Parameters
    ----------
    knl : 1D array of float
        Multipole strengths.

    tilt : 1D array of float
        Multipole tilts.

    ele : EleStruct
        Lattice element containing the multipoles.

    orbit : CoordStruct
        Particle position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Kicked particle.

    pole_type : int, optional
        Type of multipole. magnetic$ (default) or electric$.

    ref_orb_offset : bool, optional
        If present and n = 0 then the multipole simulates a zero length bend with bending angle knl.
    """

def multipole_kt_to_ab(knl: _pybmad.RealArray1D, knsl: _pybmad.RealArray1D, tn: _pybmad.RealArray1D, an: _pybmad.RealArray1D, bn: _pybmad.RealArray1D) -> None:
    """
    Wrapper for Fortran routine multipole_kt_to_ab

    Parameters
    ----------
    knl : 1D array of float
        Normal multitude component.

    knsl : 1D array of float
        Skew multitude component.

    tn : 1D array of float
        Multipole angle.

    an : 1D array of float
        Skew multipole component.

    bn : 1D array of float
        Normal multipole component.
    """

def multipole_spin_tracking(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Wrapper for Fortran routine multipole_spin_tracking

    Parameters
    ----------
    ele : EleStruct
        Element

    param : LatParamStruct
        Lat_param_struct

    orbit : CoordStruct
        Particle coordinates.
    """

def mytan(y: float, x: float) -> float:
    """
    Wrapper for Fortran routine mytan

    Parameters
    ----------
    y : float

    x : float

    Returns
    -------
    arg : float
    """

def n_attrib_string_max_len() -> int:
    """
    Routine to return the the maximum number of characters in any attribute
    name known to bmad.

    Returns
    -------
    max_len : int
        Maximum number of characters in any attribute name.
    """

def new_control(lat: _pybmad.LatStruct, ele_name: str | None = None) -> int:
    """
    Wrapper for Fortran routine new_control

    Parameters
    ----------
    lat : LatStruct
        Lat used

    ele_name : str, optional
        Name of the new element.

    Returns
    -------
    ix_ele : int
        Index of the new control element
    """

def nint_chk(re_val: float) -> int:
    """
    Returns the nearest integer to re_val.
    Also does out-of-bounds error checking.
    Used with bmad parsing.

    Parameters
    ----------
    re_val : float
        Input real number.

    Returns
    -------
    int_val : int
        Output nearest integer.
    """

def normal_form_complex_taylors(one_turn_taylor: _pybmad.TaylorStructArray1D, rf_on: bool, F: _pybmad.ComplexTaylorStructArray1D | None = None, L: _pybmad.ComplexTaylorStructArray1D | None = None, A: _pybmad.TaylorStructArray1D | None = None, A_inverse: _pybmad.TaylorStructArray1D | None = None, order: int | None = None) -> None:
    """
    Wrapper for Fortran routine normal_form_complex_taylors

    Parameters
    ----------
    one_turn_taylor : 1D array of TaylorStruct (shape: 6)

    rf_on : bool

    F : 1D array of ComplexTaylorStruct (shape: 6), optional

    L : 1D array of ComplexTaylorStruct (shape: 6), optional

    A : 1D array of TaylorStruct (shape: 6), optional

    A_inverse : 1D array of TaylorStruct (shape: 6), optional

    order : int, optional
    """

class NormalFormTaylors:
    """normal_form_taylors return type"""

    @property
    def dhdj(self) -> _pybmad.TaylorStructArray1D: ...

    @property
    def A(self) -> _pybmad.TaylorStructArray1D: ...

    @property
    def A_inverse(self) -> _pybmad.TaylorStructArray1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def normal_form_taylors(one_turn_taylor: _pybmad.TaylorStructArray1D, rf_on: bool) -> NormalFormTaylors:
    """
    Do a normal form decomposition on a one-turn taylor map M:
      M = A o R o A_inverse
    where A maps Floquet (fully normalized) coordinates to lab coordinates.
    In Floquet coordinates, the amplitudes are defined as J_i = (1/2) (x_i^2 + p_i^2).
    The map R = exp(:h:) is a pure rotation with h = h(J) is a function of the amplitudes only.
    The angles (phase advances) are given by phi_i = 2pi*dh/dJ_i.
    The taylor terms of dhdj are therefore the tunes, chromaticities, amplitude dependent tune shifts, etc.

    The mapping procedure for one turn is:
     z_Floquet_in = A_inverse o z_Lab_in
     [phi_a, phi_b, phi_c] = 2 pi * dhdj o z_Floquet_in
     z_Floquet_out = RotationMatrix(phi_a, phi_b, phi_c) . z_Floquet_in
     z_Lab_out = A o z_Floquet_out

    Parameters
    ----------
    one_turn_taylor : 1D array of TaylorStruct (shape: 6)
        one turn taylor map

    rf_on : bool
        Was the map calculated with RF on?

    Returns
    -------
    dhdj : 1D array of TaylorStruct (shape: 6), optional
        Map from Floquet coordinates to phase advances

    A : 1D array of TaylorStruct (shape: 6), optional
        Map from Floquet coordinates to Lab coordinates

    A_inverse : 1D array of TaylorStruct (shape: 6), optional
        Map from Lab coordinates to Floquet coordinates
    """

class NormalMode3Calc:
    """normal_mode3_calc return type"""

    @property
    def tune(self) -> list[float]: ...

    @property
    def B(self) -> list[list[float]]: ...

    @property
    def HV(self) -> list[list[float]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def normal_mode3_calc(t6: Sequence[Sequence[float]], above_transition: bool | None = None, abz_tunes: Sequence[float] | None = None) -> NormalMode3Calc:
    """
    Does an Eigen decomposition of the 1-turn transfer matrix (mat) and generates
    B, V, H.

    If the above_transition argument is present and false, then the 3rd (z) mode is assumed
    to have a positive slip factor (z-mode rotates counter clockwise in phase space).
    Default is True ==> z-mode has a negative slip factor so the mode rotates clock-wise in phase space.

    Parameters
    ----------
    above_transition : bool, optional
        If present and false, then z-mode assumes positive slip factor. Else negative slip factor assumed.

    abz_tunes : 1D array of float (shape: 3), optional
        Tunes to order eigensystem by.

    Returns
    -------
    tune : 1D array of float (shape: 3)
        Tunes of the 3 normal modes (radians)

    B : 2D array of float (shape: 6,6)
        B is block diagonal and related to the normal mode Twiss parameters.

    HV : 2D array of float (shape: 6,6)
        Transforms from normal mode coordinates to canonical coordinates: x = H.V.a
    """

def normal_mode_dispersion(ele: _pybmad.EleStruct, reverse: bool | None = None) -> None:
    """
    Wrapper for Fortran routine normal_mode_dispersion

    Parameters
    ----------
    ele : EleStruct
        Element whose dispersions are to be adjusted.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with adjusted dispersions.

    reverse : bool, optional
        Default is False. If True, calculate the x,y dispersions from the normal mode ones.
    """

def normalize_evecs(evec: Sequence[Sequence[complex]]) -> bool:
    """
    Normalizes eigenvectors such that transpose(E).S.E = iS, where E = evec_r + i evec_i

    Parameters
    ----------
    evec : 2D array of complex (shape: 6,6)
        complex eigenvectors arranged down columns.
        This parameter is an input/output and is modified in-place.
        As an output, evec: Eigensystem normalized to be symplectic.

    Returns
    -------
    err_flag : bool
        Set true of normalization is not possible due to amplitude is zero.
    """

def num_field_eles(ele: _pybmad.EleStruct) -> int:
    """
    Wrapper for Fortran routine num_field_eles

    Parameters
    ----------
    ele : EleStruct
        Element with sum number of associated field elements.

    Returns
    -------
    n_field_ele : int
        Number of associated field elements.
    """

def num_lords(slave: _pybmad.EleStruct, lord_type: int) -> int:
    """
    Wrapper for Fortran routine num_lords

    Parameters
    ----------
    slave : EleStruct
        Slave element.

    lord_type : int
        Type of lord. super_lord$, multipass_lord$, girder_lord$, group_lord$, overlay_lord$, and governor$ (=
        group + overlay + control + girder)

    Returns
    -------
    num : int
        Number of lords of the given type.
    """

class OdeintBmad:
    """odeint_bmad return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def odeint_bmad(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, s1_body: float, s2_body: float, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> OdeintBmad:
    """
    $OMP THREADPRIVATE(rk_ele_ptr, rk_param_ptr, rk_orbit_ptr, rk_old_orbit_ptr, rk_s_body_ptr, &
    $OMP                rk_old_s, rk_err_flag_ptr)

     Subroutine odeint_bmad (orbit, ele, param, s1_body, s2_body, err_flag, track, mat6, make_matrix)

     Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
     Recipes.  See the NR book for more details.

     Notice that this routine has an two tolerances:
       bmad_com%rel_tol_adaptive_tracking
       bmad_com%abs_tol_adaptive_tracking

     Note: For elements where the reference energy is not constant (lcavity, etc.), and
     with elements where the reference particle does not follow the reference trajectory (wigglers for example),
     the calculation of z is "off" while the particle is inside the element. At the ends there is no problem.

    Parameters
    ----------
    orbit : CoordStruct
        Starting coords: (x, px, y, py, z, delta) in element body coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coords

    ele : EleStruct
        Element to track through.

    param : LatParamStruct
        Lattice parameters.

    s1_body : float
        Starting point relative to physical entrance.

    s2_body : float
        Ending point relative physical entrance.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix propagated through the element.

    make_matrix : bool, optional
        If True then make the 6x6 transfer matrix.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise. Note: a particle getting lost, for example hitting an
        aperture, is *not* an error.

    track : TrackStruct, optional
        Structure holding the track information.
    """

class OdeintBmadTime:
    """odeint_bmad_time return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def dt_step(self) -> float: ...

    @property
    def rf_time(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def odeint_bmad_time(orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, t_dir: int, rf_time: float, track: _pybmad.TrackStruct | None = None, t_end: float | None = None, extra_field: _pybmad.EmFieldStruct | None = None) -> OdeintBmadTime:
    """
    $OMP THREADPRIVATE(tt_ele_ptr, tt_param_ptr, tt_orb_ptr, tt_old_orb_ptr, tt_extra_field_ptr, &
    $OMP                tt_old_t_ptr, tt_rf_time_ptr, tt_s_fringe_ptr, tt_vec_err_ptr)

     Subroutine odeint_bmad_time (orb, ele, param, t_dir, rf_time, err_flag, track, t_end, dt_step, extra_field)

     Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
     Recipes.  See the NR book for more details.

     Tracking is done until the particle is lost or exits the element.

    Parameters
    ----------
    orb : CoordStruct
        Starting coords: (x, px, y, py, s, ps) [t-based]
        This parameter is an input/output and is modified in-place.
        As an output, orb: Ending coords

    ele : EleStruct
        Element to track through.

    param : LatParamStruct
        Beam parameters.

    t_dir : int
        Direction of time travel = +/-1. Can be negative for patches. Will be -1 if element has a negative length.

    rf_time : float
        Time relative to RF clock.
        This parameter is an input/output and is modified in-place.
        As an output, rf_time: Updated time.

    track : TrackStruct, optional
        Structure holding the track information.

    t_end : float, optional
        If present, maximum time to which the particle will be tracked. Used for tracking with given time steps.
        The time orb.t at which tracking stops may be less than this if the particle gets to the end of the
        element

    extra_field : EmFieldStruct, optional
        Static field to be added to the element field. Eg used with space charge.

    Returns
    -------
    rf_time : float
        Time relative to RF clock.
        This parameter is an input/output and is modified in-place.
        As an output, rf_time: Updated time.

    err_flag : bool
        Set True if there is an error. False otherwise.

    dt_step : float, optional
        Next RK time step that this tracker would take based on the error tolerance. Used by track_bunch_time.
    """

class OffsetParticle:
    """offset_particle return type"""

    @property
    def s_out(self) -> float: ...

    @property
    def spin_qrot(self) -> list[float]: ...

    @property
    def time(self) -> float | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def offset_particle(ele: _pybmad.EleStruct, set: bool, orbit: _pybmad.CoordStruct, set_tilt: bool | None = None, set_hvkicks: bool | None = None, drift_to_edge: int | None = None, s_pos: float | None = None, set_spin: bool | None = None, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None, time: float | None = None) -> OffsetParticle:
    """
    Wrapper for Fortran routine offset_particle

    Parameters
    ----------
    ele : EleStruct
        Element

    set : bool
        T (= set$)   -> Translate from lab coords to the local element coords. F (= unset$) -> Translate back from
        element to lab coords.

    orbit : CoordStruct
        Coordinates of the particle.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Coordinates of particle.

    set_tilt : bool, optional
        Default is True. T -> Rotate using ele.value(tilt$) and ele.value(roll$) for sbends. F -> Do not rotate

    set_hvkicks : bool, optional
        Default is True. T -> Apply 1/2 any hkick or vkick.

    drift_to_edge : int, optional
        no$             -> Do not propagate (drift) particle. no$ is default if s_pos is present. upstream_end$
        -> Propagate to upsteam edge. This is default if set = set$ and s_pos is not present. downstream_end$ ->
        Propagate to downsteam edge. This is default if set = unset$ and s_pos is not present. Note: "edge" is
        body edge if set = set$ and is laboratory (nominal non-misaligned) edge if set = unset$

    s_pos : float, optional
        Longitudinal particle position: If set = set$: Relative to upstream end (in lab coords). If set = unset$:
        Relative to entrance end (in body coords).

    set_spin : bool, optional
        Default if False. Rotate spin coordinates? Also bmad_com.spin_tracking_on must be T to rotate.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before off setting.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix after offsets applied.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    time : float, optional
        Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
        This parameter is an input/output and is modified in-place.
        As an output, time: Updated time.

    Returns
    -------
    s_out : float, optional
        Longitudinal particle position. If set = set$: Relative to entrance end (in body coords). If set = unset$:
        Relative to upstream end (in lab coords).

    spin_qrot : 1D array of float (shape: 0:3), optional
        Spin rotation quaternion

    time : float, optional
        Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
        This parameter is an input/output and is modified in-place.
        As an output, time: Updated time.
    """

def offset_photon(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, set: bool, offset_position_only: bool | None = None, rot_mat: Sequence[Sequence[float]] | None = None) -> None:
    """
    Wrapper for Fortran routine offset_photon

    Parameters
    ----------
    ele : EleStruct
        Element

    orbit : CoordStruct
        Coordinates of the particle.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Coordinates of particle.

    set : bool
        T (= set$)   -> Translate from lab coords to element coords. F (= unset$) -> Translate from element coords
        to lab coords.

    offset_position_only : bool, optional
        If present and True, only offset the position coordinates.

    rot_mat : 2D array of float (shape: 3,3), optional
        Rotation matrix from starting coords to ending coords.
    """

def one_turn_mat_at_ele(ele: _pybmad.EleStruct, phi_a: float, phi_b: float) -> list[list[float]]:
    """
    Wrapper for Fortran routine one_turn_mat_at_ele

    Parameters
    ----------
    ele : EleStruct
        Reference element.

    phi_a : float
        "a" mode tune in radians.

    phi_b : float
        "b" mode tune in radians.

    Returns
    -------
    mat4 : 2D array of float (shape: 4,4)
        1-Turn coupled matrix.
    """

class OpenBinaryFile:
    """open_binary_file return type"""

    @property
    def iu(self) -> int: ...

    @property
    def iver(self) -> int: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def open_binary_file(file_name: str, action: str, r_name: str) -> OpenBinaryFile:
    """
    Routine to open a binary file for reading or writing.

    Parameters
    ----------
    file_name : str
        File to create.

    action : str
        'READ' or 'WRITE'

    r_name : str
        Calling routine name for error messages.

    Returns
    -------
    iu : int
        Unit number of opened file.

    iver : int
        Version number if action = 'READ'

    is_ok : bool
        Open OK?
    """

class OrbitAmplitudeCalc:
    """orbit_amplitude_calc return type"""

    @property
    def amp_a(self) -> float: ...

    @property
    def amp_b(self) -> float: ...

    @property
    def amp_na(self) -> float: ...

    @property
    def amp_nb(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def orbit_amplitude_calc(ele: _pybmad.EleStruct, orb: _pybmad.CoordStruct) -> OrbitAmplitudeCalc:
    """
    Wrapper for Fortran routine orbit_amplitude_calc

    Parameters
    ----------
    ele : EleStruct
        Element holding the Twiss parameters, dispersion and coupling info.

    orb : CoordStruct
        Orbit coordinates at the exit end of ele.

    Returns
    -------
    amp_a : float, optional
        a-mode amplitude

    amp_b : float, optional
        b-mode amplitude

    amp_na : float, optional
        a-mode, energy normalized, amplitude.

    amp_nb : float, optional
        b-mode, energy normalized, amplitude.
    """

def orbit_reference_energy_correction(orbit: _pybmad.CoordStruct, p0c_new: float, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine orbit_reference_energy_correction

    Parameters
    ----------
    orbit : CoordStruct
        Coordinates to correct.

    p0c_new : float
        New reference momentum.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before correction.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix including correction.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def orbit_to_floor_phase_space(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> list[float]:
    """
    Wrapper for Fortran routine orbit_to_floor_phase_space

    Parameters
    ----------
    orbit : CoordStruct
        Particle orbit in local (not element) coordinates.

    ele : EleStruct
        Lattice element particle is in.

    Returns
    -------
    floor_phase_space : 1D array of float (shape: 6)
        Floor phase space
    """

def orbit_to_local_curvilinear(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, z_direction: int | None = None, relative_to: int | None = None) -> _pybmad.FloorPositionStruct:
    """
    Wrapper for Fortran routine orbit_to_local_curvilinear

    Parameters
    ----------
    orbit : CoordStruct
        Particle orbit in laboratory (not body) coordinates.

    ele : EleStruct
        Lattice element particle is in.

    z_direction : int, optional
        Set to +1 or -1.  Z-direction of particle velocity relative to element z-axis. Default is ele.orientation
        * orbit.direction.

    relative_to : int, optional
        not_set$ (default), upstream_end$, downstream_end$. If not_set$ then origin is at the entrance end.

    Returns
    -------
    local_position : FloorPositionStruct
        Position in local coordinates.
    """

class OrbitTooLarge:
    """orbit_too_large return type"""

    @property
    def param(self) -> _pybmad.LatParamStruct: ...

    @property
    def is_too_large(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def orbit_too_large(orbit: _pybmad.CoordStruct, check_momentum: bool | None = None) -> OrbitTooLarge:
    """
    Wrapper for Fortran routine orbit_too_large

    Parameters
    ----------
    orbit : CoordStruct
        Particle orbit.

    check_momentum : bool, optional
        If True (default) check the momentum.

    Returns
    -------
    is_too_large : bool
        True if orbit is too large. False otherwise.

    param : LatParamStruct, optional
    """

class OrderEvecsByNSimilarity:
    """order_evecs_by_n_similarity return type"""

    @property
    def evec(self) -> list[list[complex]]: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def order_evecs_by_n_similarity(eval: Sequence[complex], mat_tunes: Sequence[float], Nmat: Sequence[Sequence[float]]) -> OrderEvecsByNSimilarity:
    """
    This subroutine orderes the eigensystem such that Nmat.mat_symp_conj(N) is closest
    to the identity.  Nmat is supplied externally.

    Parameters
    ----------
    eval : 1D array of complex (shape: 6)
        complex eigenvalues.

    mat_tunes : 1D array of float (shape: 3)
        Three normal mode tunes, in radians.
        This parameter is an input/output and is modified in-place.
        As an output, mat_tunes: Ordered normal mode tunes, in radians.

    Nmat : 2D array of float (shape: 6,6)
        Normalized, real eigen matrix from make_N.

    Returns
    -------
    evec : 2D array of complex (shape: 6,6)
        complex eigenvectors arranged down columns.

    err_flag : bool
        Set True if there is an error. False otherwise
    """

def order_evecs_by_plane_dominance(evec: Sequence[Sequence[complex]], eval: Sequence[complex], mat_tunes: Sequence[float] | None = None) -> None:
    """
    This subroutine orderes the eigensystem according to which modes dominate the horizontal,
    vertical, and longitudinal planes.  This subroutine works well in machines
    that are not strongly coupled.  In machines with strong coupling, where the relation
    between the three eigenmodes a, b, c and the three lab coordinates x, y, z can change
    through the machine, this subroutine will not provide consistent ordering.

    Parameters
    ----------
    evec : 2D array of complex (shape: 6,6)
        complex eigenvectors arranged down columns.
        This parameter is an input/output and is modified in-place.
        As an output, evec: Ordered complex eigenvectors.

    eval : 1D array of complex (shape: 6)
        complex eigenvalues.
        This parameter is an input/output and is modified in-place.
        As an output, eval: Ordered complex eigenvalues.

    mat_tunes : 1D array of float (shape: 3), optional
        Three normal mode tunes, in radians.
        This parameter is an input/output and is modified in-place.
        As an output, mat_tunes: Reordered same as evecs.
    """

def order_evecs_by_tune(evec: Sequence[Sequence[complex]], eval: Sequence[complex], mat_tunes: Sequence[float], abz_tunes: Sequence[float]) -> bool:
    """
    This subroutine orders the eigensystem by matching the tunes of the eigensystem to
    externally supplied tunes abz_tunes.  abz_tunes is in radians.

    Parameters
    ----------
    evec : 2D array of complex (shape: 6,6)
        complex eigenvectors arranged down columns.
        This parameter is an input/output and is modified in-place.
        As an output, evec: Ordered eigenvectors.

    eval : 1D array of complex (shape: 6)
        complex eigenvalues.
        This parameter is an input/output and is modified in-place.
        As an output, eval: Ordered eigenvalues.

    mat_tunes : 1D array of float (shape: 3)
        Three normal mode tunes, in radians.

    abz_tunes : 1D array of float (shape: 3)
        Tunes to order eigensystem by.

    Returns
    -------
    err_flag : bool
        Set to true if an error occured.
    """

def order_particles_in_z(bunch: _pybmad.BunchStruct) -> None:
    """
    Routine to order the particles longitudinally in terms of decreasing %vec(5).
    That is from large z (head of bunch) to small z.
    Only live particles are ordered.

    Parameters
    ----------
    bunch : BunchStruct
        collection of particles.
    """

def order_super_lord_slaves(lat: _pybmad.LatStruct, ix_lord: int) -> None:
    """
    Wrapper for Fortran routine order_super_lord_slaves

    Parameters
    ----------
    lat : LatStruct
        Lat.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lat with fixed controls.

    ix_lord : int
        Index of lord element.
    """

def osc_alloc_freespace_array(nlo: Sequence[int], nhi: Sequence[int], npad: Sequence[int]) -> None:
    """
    Wrapper for Fortran routine osc_alloc_freespace_array

    Parameters
    ----------
    nlo : 1D array of int (shape: 3)

    nhi : 1D array of int (shape: 3)

    npad : 1D array of int (shape: 3)
    """

def osc_alloc_image_array(nlo: Sequence[int], nhi: Sequence[int], npad: Sequence[int]) -> None:
    """
    Wrapper for Fortran routine osc_alloc_image_array

    Parameters
    ----------
    nlo : 1D array of int (shape: 3)

    nhi : 1D array of int (shape: 3)

    npad : 1D array of int (shape: 3)
    """

def osc_alloc_rectpipe_arrays(nlo: Sequence[int], nhi: Sequence[int], npad: Sequence[int]) -> None:
    """
    Wrapper for Fortran routine osc_alloc_rectpipe_arrays

    Parameters
    ----------
    nlo : 1D array of int (shape: 3)

    nhi : 1D array of int (shape: 3)

    npad : 1D array of int (shape: 3)
    """

def osc_getgrnpipe(gam: float, a: float, b: float, delta: Sequence[float], umin: Sequence[float], npad: Sequence[int]) -> None:
    """
    Wrapper for Fortran routine osc_getgrnpipe

    Parameters
    ----------
    gam : float

    a : float

    b : float

    delta : 1D array of float (shape: 3)

    umin : 1D array of float (shape: 3)

    npad : 1D array of int (shape: 3)
    """

def osc_read_rectpipe_grn() -> None:
    """Wrapper for Fortran routine osc_read_rectpipe_grn"""

def osc_write_rectpipe_grn(apipe: float, bpipe: float, delta: Sequence[float], umin: Sequence[float], umax: Sequence[float], nlo: Sequence[int], nhi: Sequence[int], gamma: float) -> None:
    """
    Wrapper for Fortran routine osc_write_rectpipe_grn

    Parameters
    ----------
    apipe : float

    bpipe : float

    delta : 1D array of float (shape: 3)

    umin : 1D array of float (shape: 3)

    umax : 1D array of float (shape: 3)

    nlo : 1D array of int (shape: 3)

    nhi : 1D array of int (shape: 3)

    gamma : float
    """

def p_func(E_in: float) -> float:
    """
    Wrapper for Fortran routine p_func

    Parameters
    ----------
    E_in : float

    Returns
    -------
    rr1 : float
    """

def parse_cartesian_map(ct_map: _pybmad.CartesianMapStruct, ele: _pybmad.EleStruct, lat: _pybmad.LatStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    Subroutine to parse a "cartesian_map = {}" construct

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is private to bmad_parser_mod.
    This must read in:
    {type = ,
       dr = ,
       r0 = ,
       pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) )
       .
       .
       . ) },
    """

def parse_cylindrical_map(cl_map: _pybmad.CylindricalMapStruct, ele: _pybmad.EleStruct, lat: _pybmad.LatStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine parse_cylindrical_map

    Parameters
    ----------
    cl_map : CylindricalMapStruct

    ele : EleStruct

    lat : LatStruct

    delim : str

    delim_found : bool

    err_flag : bool
    """

def parse_gen_grad_map(gg_map: _pybmad.GenGradMapStruct, ele: _pybmad.EleStruct, lat: _pybmad.LatStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """Subroutine to parse a "gen_grad_map = {}" construct"""

def parse_grid_field(g_field: _pybmad.GridFieldStruct, ele: _pybmad.EleStruct, lat: _pybmad.LatStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine parse_grid_field

    Parameters
    ----------
    g_field : GridFieldStruct

    ele : EleStruct

    lat : LatStruct

    delim : str

    delim_found : bool

    err_flag : bool
    """

def parse_integer_list(err_str: str, lat: _pybmad.LatStruct, int_array: _pybmad.IntArray1D, exact_size: bool, delim: str, delim_found: bool, open_delim: str | None = None, separator: str | None = None, close_delim: str | None = None, default_value: int | None = None) -> bool:
    """
    separator, close_delim, default_value) result (is_ok)

    Routine to parse a list of integers of the form:
       open_delim integer_1 separator integer_2 . . . close_delim
    Example:   "(1.2, 2.3, 4.4, 8.5)"

    Similar to parse_integer_list2 except does not use allocatable array.
    See parse_integer_list2 for more details
    """

class ParseIntegerList2:
    """parse_integer_list2 return type"""

    @property
    def num_found(self) -> int: ...

    @property
    def delim(self) -> str: ...

    @property
    def delim_found(self) -> bool: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parse_integer_list2(err_str: str, lat: _pybmad.LatStruct, int_array: _pybmad.IntAlloc1D, num_expected: int | None = None, open_delim: str | None = None, separator: str | None = None, close_delim: str | None = None, default_value: int | None = None) -> ParseIntegerList2:
    """
    open_delim, separator, close_delim, default_value) result (is_ok)

    Routine to parse a list of integers of the form
       open_delim integer_1 separator integer_2 . . . close_delim
    Example:   (1, 2, 4, 8)

    Parameters
    ----------
    err_str : str
        Error string to print if there is an error.

    lat : LatStruct
        lattice

    int_array : 1D array of int
        the array to be read in Optional:
        This parameter is an input/output and is modified in-place.
        As an output, int_array: Array of values.

    Returns
    -------
    num_found : int
        number of elements.

    delim : str
        Delimiter found where the parsing of the input line stops.

    delim_found : bool
        Delimiter found? False if end of input command.

    is_ok : bool
        Set True if everything is ok.
    """

def parse_line_or_list(sequence: _pybmad.SeqStructAlloc1D, iseq_tot: int, lat: _pybmad.LatStruct, top_level: bool) -> None:
    """
    Subroutine to parse a sequence.
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

class ParseRealList:
    """parse_real_list return type"""

    @property
    def delim(self) -> str: ...

    @property
    def delim_found(self) -> bool: ...

    @property
    def num_found(self) -> int: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parse_real_list(lat: _pybmad.LatStruct, err_str: str, real_array: _pybmad.RealArray1D, exact_size: bool, open_delim: str | None = None, separator: str | None = None, close_delim: str | None = None, default_value: float | None = None) -> ParseRealList:
    """
    separator, close_delim, default_value, num_found) result (is_ok)

    Routine to parse a list of reals of the form:
       open_delim real_1 separator real_2 . . . close_delim
    Example:   "(1.2, 2.3, 4.4, 8.5)"

    Similar to parse_real_list2 except does not use allocatable array.
    Also see: parse_real_matrix.

    Parameters
    ----------
    lat : LatStruct
        Lattice

    err_str : str
        Error string to print if there is an error.

    real_array : 1D array of float

    exact_size : bool

    open_delim : str, optional

    separator : str, optional

    close_delim : str, optional

    default_value : float, optional

    Returns
    -------
    delim : str

    delim_found : bool

    num_found : int, optional
    """

class ParseRealList2:
    """parse_real_list2 return type"""

    @property
    def num_found(self) -> int: ...

    @property
    def delim(self) -> str: ...

    @property
    def delim_found(self) -> bool: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parse_real_list2(lat: _pybmad.LatStruct, err_str: str, real_array: _pybmad.RealAlloc1D, num_expected: int | None = None, open_brace: str | None = None, separator: str | None = None, close_brace: str | None = None, default_value: float | None = None, single_value: bool | None = None) -> ParseRealList2:
    """
    open_delim, separator, close_delim, default_value, single_value) result (is_ok)

    Routine to parse a list of reals of the form:
       open_brace real_1 separator real_2 . . . close_brace
    Example:   "(1.2, 2.3, 4.4, 8.5)"

    Also see:
      pase_real_list
      parse_real_matrix.

    Parameters
    ----------
    lat : LatStruct
        lattice

    err_str : str
        Error string to print if there is an error.

    real_array : 1D array of float
        the array to be read in
        This parameter is an input/output and is modified in-place.
        As an output, real_array: Array of values

    Returns
    -------
    num_found : int
        number of elements

    delim : str
        Delimiter found where the parsing of the input line stops.

    delim_found : bool
        Stopping delimiter found? False if end of input command.

    is_ok : bool
        Set True if everything is ok
    """

def parse_superimpose_command(lat: _pybmad.LatStruct, ele: _pybmad.EleStruct, pele: _pybmad.ParserEleStruct, delim: str) -> None:
    """No docstring available."""

def parser2_add_superimpose(lat: _pybmad.LatStruct, super_ele_in: _pybmad.EleStruct, pele: _pybmad.ParserEleStruct, in_lat: _pybmad.LatStruct | None = None) -> None:
    """
    Wrapper for Fortran routine parser2_add_superimpose

    Parameters
    ----------
    lat : LatStruct

    super_ele_in : EleStruct

    pele : ParserEleStruct

    in_lat : LatStruct, optional
    """

def parser_add_branch(fork_ele: _pybmad.EleStruct, lat: _pybmad.LatStruct, sequence: _pybmad.SeqStructAlloc1D, seq_name: _pybmad.CharacterAlloc1D, seq_indexx: _pybmad.IntAlloc1D, no_end_marker: bool, in_lat: _pybmad.LatStruct, plat: _pybmad.ParserLatStruct, created_new_branch: bool, new_branch_name: str | None = None) -> None:
    """
    seq_name, seq_indexx, no_end_marker, in_lat, plat, created_new_branch)

    Subroutine to do line expansion.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

def parser_add_constant(word: str, lat: _pybmad.LatStruct, redef_is_error: bool) -> None:
    """
    Wrapper for Fortran routine parser_add_constant

    Parameters
    ----------
    word : str

    lat : LatStruct

    redef_is_error : bool
    """

class ParserAddLords:
    """parser_add_lords return type"""

    @property
    def lat(self) -> _pybmad.LatStruct: ...

    @property
    def check_lat(self) -> _pybmad.LatStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parser_add_lords(lord_lat: _pybmad.LatStruct, n_ele_max: int, plat: _pybmad.ParserLatStruct) -> ParserAddLords:
    """
    Subroutine to add overlay, group, and girder lords to the lattice.
    For overlays and groups: If multiple elements have the same name then
    use all of them.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Parameters
    ----------
    lord_lat : LatStruct
        List of lord elements to add to lat.

    n_ele_max : int
        lord elements in lord_lat are in range [1:n_ele_max].

    plat : ParserLatStruct
        Extra info needed to place the lord elements

    Returns
    -------
    lat : LatStruct
        Lattice to add lord elements to.

    check_lat : LatStruct, optional
        If slave elements of a lord are not in lat but are in check_lat, do not issue error message about slave
        elements not found.
    """

def parser_add_superimpose(branch: _pybmad.BranchStruct, super_ele_in: _pybmad.EleStruct, pele: _pybmad.ParserEleStruct, in_lat: _pybmad.LatStruct, plat: _pybmad.ParserLatStruct) -> None:
    """
    Wrapper for Fortran routine parser_add_superimpose

    Parameters
    ----------
    branch : BranchStruct

    super_ele_in : EleStruct

    pele : ParserEleStruct

    in_lat : LatStruct

    plat : ParserLatStruct
    """

def parser_call_check(word: str, ix_word: int, delim: str, delim_found: bool, call_found: bool, err_flag: bool | None = None) -> None:
    """
    Routine to check if there is a "call::XXX" construct in the input stream.
    """

def parser_debug_print_info(lat: _pybmad.LatStruct, debug_line: str, sequence: _pybmad.SeqStructArray1D | None = None) -> None:
    """
    Subroutine to remove all null_ele elements.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

def parser_error(what1: str, what2: str | None = None, what3: str | None = None, what4: str | None = None, seq: _pybmad.SeqStruct | None = None, pele: _pybmad.ParserEleStruct | None = None, stop_here: bool | None = None, level: int | None = None, r_array: _pybmad.RealArray1D | None = None, i_array: _pybmad.IntArray1D | None = None) -> None:
    """
    Routine to print an error message generated when parsing a lattice.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Parameters
    ----------
    what1 : str
        First line in error message.

    what2 : str, optional
        Second line in error message.

    what3 : str, optional
        Third line in error message.

    what4 : str, optional
        Fourth line in error message.

    seq : SeqStruct, optional
        Used when error is generated during reading of a lattice file. Contains information such as file name, and
        line number where error was detected.

    pele : ParserEleStruct, optional
        Used when error is associated with a lattice element. Contains information on the lattice element.

    stop_here : bool, optional
        If present and True then immediately stop. Used with severe errors.

    level : int, optional
        Possibilities are:

    r_array : 1D array of float, optional
        Real numbers to be encoded in error message. See out_io doc.

    i_array : 1D array of int, optional
        Integer numbers to be encoded in error message. See out_io doc.
    """

class ParserExpandLine:
    """parser_expand_line return type"""

    @property
    def n_ele_expand(self) -> int: ...

    @property
    def expanded_line(self) -> _pybmad.BaseLineEleStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parser_expand_line(i_lev: int, line_name: str, sequence: _pybmad.SeqStructAlloc1D, seq_name: _pybmad.CharacterAlloc1D, seq_indexx: _pybmad.IntAlloc1D, no_end_marker: bool, lat: _pybmad.LatStruct | None = None, in_lat: _pybmad.LatStruct | None = None) -> ParserExpandLine:
    """
    seq_name, seq_indexx, no_end_marker, n_ele_expand, lat, in_lat, expanded_line)

    Subroutine to do line expansion.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Note: Either lat and in_lat must be present or expanded_line must be present.

    Parameters
    ----------
    i_lev : int
        Subsequence level. 1 => Root level.

    line_name : str
        Root line to expand.

    sequence : 1D array of SeqStruct
        Array of sequencies.

    seq_name : 1D array of str
        Array of sequence names.

    seq_indexx : 1D array of int
        Index array for the sequence names.

    no_end_marker : bool
        Put a marker named "end" at the end of the branch?

    lat : LatStruct, optional
        Lattice to put the expanded line
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with new line. Except if expanded_line is present.

    in_lat : LatStruct, optional
        Lattice with array of defined elements.

    Returns
    -------
    n_ele_expand : int
        Number of elements in the finished line.

    expanded_line : 1D array of BaseLineEleStruct, optional
        If present, lat argument will be ignored and the expanded line will be put into expanded_line.
    """

class ParserFastComplexRead:
    """parser_fast_complex_read return type"""

    @property
    def delim(self) -> str: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parser_fast_complex_read(cmplx_vec: _pybmad.ComplexArray1D, ele: _pybmad.EleStruct, err_str: str) -> ParserFastComplexRead:
    """
    Routine to read an array of complex numbers.

    This routine assumes that the array values are pure numbers in the form "<re>" or "(<re> <im>)"
    where <re> and <im> are real numbers (not expressions) and there are no commas except possibly
    at the end of the array.

    Parameters
    ----------
    cmplx_vec : 1D array of complex
        Complex vector.

    ele : EleStruct
        Lattice element associated with the array. Used for error messages.

    err_str : str
        String used when printing error messages identifying where in the lattice file the error is occuring.

    Returns
    -------
    delim : str
        Delimitor at end of array. Must be "," or "}"

    is_ok : bool
        True if everything OK. False otherwise.
    """

def parser_fast_integer_read(int_vec: _pybmad.IntArray1D, ele: _pybmad.EleStruct, delim_wanted: str, err_str: str) -> bool:
    """No docstring available."""

class ParserFastRealRead:
    """parser_fast_real_read return type"""

    @property
    def delim(self) -> str: ...

    @property
    def n_real(self) -> int: ...

    @property
    def is_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parser_fast_real_read(real_vec: _pybmad.RealArray1D, ele: _pybmad.EleStruct, end_delims: str, err_str: str, exact_size: bool | None = None) -> ParserFastRealRead:
    """
    Routine to read an array of real numbers.

    This routine assumes that the array values are pure numbers in the form "<re1> <re2> ...,"
    where <re1>, <re2>, etc. are real numbers (not expressions) and there are no commas except possibly,
    at the end of the array.

    Note: if end_delim is "," and next character is a delim but not ",", the next character is taken as the delim.

    Parameters
    ----------
    real_vec : 1D array of float
        Real vector.

    ele : EleStruct
        Lattice element associated with the array. Used for error messages.

    end_delims : str
        List of possible ending delimitors.

    err_str : str
        String used when printing error messages identifying where in the lattice file the error is occuring.

    exact_size : bool, optional
        If True (default), number of values must match real_vec size.

    Returns
    -------
    delim : str
        Delimitor at end of array.

    is_ok : bool
        True if everything OK. False otherwise.

    n_real : int, optional
        Number of elements found.
    """

def parser_file_stack(how: str, file_name_in: str | None = None, finished: bool | None = None, err: bool | None = None, open_file: bool | None = None, abort_on_open_error: bool | None = None) -> None:
    """
    Subroutine to keep track of the files that are opened for reading.
    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

def parser_get_integer(int_val: int, word: str, ix_word: int, delim: str, delim_found: bool, err: bool, str1: str | None = None, str2: str | None = None) -> None:
    """
    Wrapper for Fortran routine parser_get_integer

    Parameters
    ----------
    int_val : int

    word : str

    ix_word : int

    delim : str

    delim_found : bool

    err : bool

    str1 : str, optional

    str2 : str, optional
    """

def parser_get_logical(attrib_name: str, this_logic: bool, ele_name: str, delim: str, delim_found: bool, err: bool) -> None:
    """
    Wrapper for Fortran routine parser_get_logical

    Parameters
    ----------
    attrib_name : str

    this_logic : bool

    ele_name : str

    delim : str

    delim_found : bool

    err : bool
    """

def parser_identify_fork_to_element(lat: _pybmad.LatStruct) -> None:
    """
    Routine to identify the elements the forks in a lattice are branching to.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

def parser_init_custom_elements(lat: _pybmad.LatStruct) -> None:
    """No docstring available."""

def parser_print_line(lat: _pybmad.LatStruct, end_of_file: bool) -> None:
    """
    This routine is called when a print statement is found in the lattice file.
    """

def parser_read_lr_wake(ele: _pybmad.EleStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    Subroutine to read in a long-range wake field from an external file.
    This subroutine is used by bmad_parser and bmad_parser2.

    Parameters
    ----------
    ele : EleStruct
        Element containing wake structure.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with wake information.
    """

def parser_read_old_format_lr_wake(ele: _pybmad.EleStruct, lr_file_name: str) -> None:
    """
    Subroutine to read in a long-range wake field from an external file.
    This subroutine is used by bmad_parser and bmad_parser2.

    Parameters
    ----------
    ele : EleStruct
        Element containing wake structure.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with wake information.

    lr_file_name : str
        Name of long-range wake field file.
    """

def parser_read_old_format_sr_wake(ele: _pybmad.EleStruct, sr_file_name: str) -> None:
    """
    Subroutine to read in a short-range wake field from an external file.
    This subroutine is used by bmad_parser and bmad_parser2.

    Parameters
    ----------
    ele : EleStruct
        Element containing wake structure.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with wake information.

    sr_file_name : str
        Name of short-range wake field file.
    """

def parser_read_sr_wake(ele: _pybmad.EleStruct, delim: str, delim_found: bool, err_flag: bool) -> None:
    """
    Subroutine to read in a short-range wake field.
    This subroutine is used by bmad_parser and bmad_parser2.

    Parameters
    ----------
    ele : EleStruct
        Element containing wake structure.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with wake information.
    """

class ParserSetAttribute:
    """parser_set_attribute return type"""

    @property
    def delim(self) -> str: ...

    @property
    def delim_found(self) -> bool: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def pele(self) -> _pybmad.ParserEleStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def parser_set_attribute(how: int, ele: _pybmad.EleStruct, check_free: bool | None = None, heterogeneous_ele_list: bool | None = None, set_field_master: bool | None = None) -> ParserSetAttribute:
    """
    heterogeneous_ele_list, set_field_master)

    Subroutine used by bmad_parser and bmad_parser2 to get the value of
    an attribute from the input file and set the appropriate value in an element.

    This subroutine is not intended for general use.

    Parameters
    ----------
    how : int
        Either def$ if the element is being construct from scratch or redef$ if the element has already been
        formed and this is part of a "ele_name[attrib_name] = value" construct.

    ele : EleStruct
        Element whose attribute this is.

    check_free : bool, optional
        If present and True then an error will be generated if the attribute is not free to vary. Used by
        bmad_parser2.

    heterogeneous_ele_list : bool, optional
        If True (default = False), we are parsing something like something like "*[tracking_method] =
        runge_kutta". In this case, runge_kutta may not be valid for this ele but this is not an error.

    set_field_master : bool, optional
        If True (the default) set ele.field_master = T if the attribute to be set is something like DB_FIELD. If
        this routine is being called post lattice parsing, setting ele.field_master is *not* wanted.

    Returns
    -------
    delim : str
        Delimiter found where the parsing of the input line stops.

    delim_found : bool
        Delimiter found? False if end of input command.

    err_flag : bool
        Set True if there is a problem parsing the input.

    pele : ParserEleStruct, optional
        Structure to hold additional information that cannot be stored in the ele argument.
    """

def parser_transfer_control_struct(con_in: _pybmad.ControlStruct, lord: _pybmad.EleStruct, ix_var: int) -> _pybmad.ControlStruct:
    """
    Routine to transfer the information from an input control_struct (which stores
    the user input parameters) to a control_struct that will be stored in the lat%control
    or lord%control%ramp for a ramper.

    Parameters
    ----------
    con_in : ControlStruct
        Input control structure.

    lord : EleStruct
        Lord element associated with the control_struct.

    ix_var : int
        If an expression stack evaluates to a constant, this routine will modify the expression stack to evaluate
        to the value of: lord.control.var(ix_var) * constant

    Returns
    -------
    con_out : ControlStruct
        Output control structure.
    """

def particle_in_global_frame(orb: _pybmad.CoordStruct, branch: _pybmad.BranchStruct, in_time_coordinates: bool | None = None, in_body_frame: bool | None = None, w_mat_out: Sequence[Sequence[float]] | None = None) -> _pybmad.CoordStruct:
    """
    Returns the particle in global time coordinates given is coordinates orb in lattice lat.

    Parameters
    ----------
    orb : CoordStruct
        particle in s-coordinates

    branch : BranchStruct
        branch that contains branch.ele(orb.ix_ele)

    in_time_coordinates : bool, optional
        Default is false. If true, orb will taken as in time coordinates.

    in_body_frame : bool, optional
        Default is true. If false, ele offsets will be ignored.

    Returns
    -------
    particle : CoordStruct
        particle in global time coordinates
    """

def particle_is_moving_backwards(orbit: _pybmad.CoordStruct) -> bool:
    """
    Wrapper for Fortran routine particle_is_moving_backwards

    Parameters
    ----------
    orbit : CoordStruct
        Particle coordinates

    Returns
    -------
    is_moving_backwards : bool
        True if moving backward. False otherwise.
    """

def particle_is_moving_forward(orbit: _pybmad.CoordStruct, dir: int | None = None) -> bool:
    """
    Wrapper for Fortran routine particle_is_moving_forward

    Parameters
    ----------
    orbit : CoordStruct
        Particle coordinates

    dir : int, optional
        +1 if tracking forward(default) or -1 to return True if tracking backwards.

    Returns
    -------
    is_moving_forward : bool
        True if moving forward. False otherwise.
    """

def particle_rf_time(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, reference_active_edge: bool | None = None, s_rel: float | None = None, time_coords: bool | None = None, rf_freq: float | None = None, abs_time: bool | None = None) -> float:
    """
    Wrapper for Fortran routine particle_rf_time

    Parameters
    ----------
    orbit : CoordStruct
        Particle coordinates

    ele : EleStruct
        Element being tracked through.

    reference_active_edge : bool, optional
        If True, and ele is a rfcavity or lcavity, use the active edge (edge of the region with non-zero field) as
        the reference point.

    s_rel : float, optional
        Longitudinal position relative to the upstream edge of the element. Needed for relative time tracking when
        the particle is inside the element. Default is 0.

    time_coords : bool, optional
        Default False. If True then orbit is using time based phase space coordinates.

    rf_freq : float, optional
        If present, the returned time shifted by an integer multiple of 1/rf_freq to be in the range
        [-1/2*rf_freq, 1/2*rf_freq]. This is useful to avoid round-off errors.

    abs_time : bool, optional
        If False (default) use setting of bmad_com.absolute_time_tracking. If True, use absolute time instead of
        relative time.

    Returns
    -------
    time : float
        Current time.
    """

def patch_flips_propagation_direction(x_pitch: float, y_pitch: float) -> bool:
    """
    Wrapper for Fortran routine patch_flips_propagation_direction

    Parameters
    ----------
    x_pitch : float
        Rotaion around y-axis

    y_pitch : float
        Rotation around x-axis.

    Returns
    -------
    is_flip : bool
        True if patch does a flip
    """

def patch_length(patch: _pybmad.EleStruct, ref_coords: int | None = None) -> float:
    """
    Wrapper for Fortran routine patch_length

    Parameters
    ----------
    patch : EleStruct
        Patch element.

    ref_coords : int, optional
        Reference coords to use. entrance_end$, exit_end$ Default is nint(patch.value(ref_coords$)).

    Returns
    -------
    length : float
        Length of patch.
    """

class PhotonAbsorptionAndPhaseShift:
    """photon_absorption_and_phase_shift return type"""

    @property
    def absorption(self) -> float: ...

    @property
    def phase_shift(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def photon_absorption_and_phase_shift(material: str, Energy: float) -> PhotonAbsorptionAndPhaseShift:
    """
    Routine to calcualte the absorption and phase shift values for a photon with a given
    energy going through a particular material.

    Parameters
    ----------
    material : str
        Material name.

    Energy : float
        Photon energy (eV).

    Returns
    -------
    absorption : float
        E_field ~ Exp(-absorption * length)

    phase_shift : float
        E_field Phase shift (radians) per unit length relative to vacuum.

    err_flag : bool
        Set true if material not recognized.
    """

class PhotonAddToDetectorStatistics:
    """photon_add_to_detector_statistics return type"""

    @property
    def ix_pt(self) -> int: ...

    @property
    def iy_pt(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def photon_add_to_detector_statistics(orbit0: _pybmad.CoordStruct, orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, pixel_pt: _pybmad.PixelPtStruct | None = None) -> PhotonAddToDetectorStatistics:
    """
    Routine to add photon statistics to the appropriate pixel of a "detector" grid.

    It is assumed that track_to_surface has been called so that the photon is at the
    detector surface and that orbit%vec(1) and %vec(3) are in surface coords (needed for curved detectors).

    Parameters
    ----------
    orbit0 : CoordStruct
        Photon coords at beginning of lattice

    orbit : CoordStruct
        Photon coords at the detector.

    ele : EleStruct
        Element with grid.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with updatted grid.

    pixel_pt : PixelPtStruct, optional
        If present then use this grid point instead of the grid point determined by the (x, y) coords of the
        photon

    Returns
    -------
    ix_pt : int, optional
        Index of upgraded ele.photon.surface.grid.pt(:,:) point. These arguments are not set if the pixel_pt
        argument is present.

    iy_pt : int, optional
        Index of upgraded ele.photon.surface.grid.pt(:,:) point. These arguments are not set if the pixel_pt
        argument is present.
    """

class PhotonDiffuseScattering:
    """photon_diffuse_scattering return type"""

    @property
    def graze_angle_out(self) -> float: ...

    @property
    def phi_out(self) -> float: ...

    @property
    def diffuse_param(self) -> _pybmad.DiffuseParamStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def photon_diffuse_scattering(graze_angle_in: float, energy: float, surface: _pybmad.PhotonReflectSurfaceStruct) -> PhotonDiffuseScattering:
    """
    Routine to simulate the diffuse scattering of photons. The outgoing angles are
    choosen using the Dugan distribution.

    Also see: photon_reflection.
    Use photon_reflection_std_surface_init or read_surface_reflection_file to get surface info.

    Parameters
    ----------
    graze_angle_in : float
        Incident grazing (not polar) angle in radians.

    energy : float
        Photon energy in eV.

    surface : PhotonReflectSurfaceStruct
        surface info

    Returns
    -------
    graze_angle_out : float
        graze_angle in radians.

    phi_out : float
        Azimuthal angle in radians.

    diffuse_param : DiffuseParamStruct, optional
        Internal parameters used in the calculation. This is used for diagnostics and is not used in standard
        simulations.
    """

class PhotonHitFunc:
    """photon_hit_func return type"""

    @property
    def status(self) -> int: ...

    @property
    def d_radius(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def photon_hit_func(track_len: float) -> PhotonHitFunc:
    """
    Routine to be used as an argument in zbrent in the capillary_photon_hit_spot_calc.
    Made a module procedure (not nested) to avoid a stack trampoline.

    Parameters
    ----------
    track_len : float
        Place to position the photon.

    Returns
    -------
    status : int
        Not set.

    d_radius : float
        r_photon - r_wall
    """

def photon_read_spline(spline_dir: str) -> _pybmad.PhotonInitSplinesStruct:
    """
    Routine to initialize a photon using a set of spline fits.

    Parameters
    ----------
    spline_dir : str
        Root directory for the spline fits.

    Returns
    -------
    splines : PhotonInitSplinesStruct
        Spline structure
    """

class PhotonReflection:
    """photon_reflection return type"""

    @property
    def graze_angle_out(self) -> float: ...

    @property
    def phi_out(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def photon_reflection(graze_angle_in: float, energy: float, surface: _pybmad.PhotonReflectSurfaceStruct) -> PhotonReflection:
    """
    Routine to reflect a photon from a surface including both diffuse and specular reflections.

    Parameters
    ----------
    graze_angle_in : float
        Incident grazing (not polar) angle in radians.

    energy : float
        Photon energy in eV.

    surface : PhotonReflectSurfaceStruct
        surface info

    Returns
    -------
    graze_angle_out : float
        graze_angle in radians.

    phi_out : float
        Azimuthal angle in radians.
    """

def photon_reflection_std_surface_init() -> _pybmad.PhotonReflectSurfaceStruct:
    """
    $OMP THREADPRIVATE(dr_d_param_ptr, dr_surface_ptr, dr_old_integral, dr_tot_integral, &
    $OMP                dr_ran1, dr_ran2, dr_j)

     Subroutine photon_reflection_std_surface_init (surface)

     Routine to initialize the standard proton reflection probability tables.
     The standard tables are for 10 nm C film on Al substrate.
     The surface roughness for diffuse scattering is 200 nm and the
     the surface roughness correlation length is 5.5 um.

    Returns
    -------
    surface : PhotonReflectSurfaceStruct
        photon_reflect_surface_struct
    """

class PhotonReflectivity:
    """photon_reflectivity return type"""

    @property
    def p_reflect(self) -> float: ...

    @property
    def rel_p_specular(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def photon_reflectivity(angle: float, energy: float, surface: _pybmad.PhotonReflectSurfaceStruct) -> PhotonReflectivity:
    """
    Routine to evaluate the photon reflectivity.
      probability of absorption          = 1 - p_reflect
      probability of reflection          = p_reflect
      probability of specular reflection = p_reflect * rel_p_specular
      probability of diffuse reflection  = p_reflect * (1 - rel_p_specular)

    Use photon_reflection_std_surface_init or read_surface_reflection_file to get surface info.

    Parameters
    ----------
    angle : float
        Incident grazing angle in radians.

    energy : float
        Photon energy in eV.

    surface : PhotonReflectSurfaceStruct
        surface info

    Returns
    -------
    p_reflect : float
        Reflection probability.

    rel_p_specular : float
        Relative specular reflection probability.
    """

def photon_target_corner_calc(aperture_ele: _pybmad.EleStruct, x_lim: float, y_lim: float, z_lim: float, source_ele: _pybmad.EleStruct) -> _pybmad.TargetPointStruct:
    """
    Routine to calculate the corner coords in the source_ele ref frame.

    Parameters
    ----------
    aperture_ele : EleStruct
        Element containing the aperture

    x_lim : float
        Transverse corner points in aperture_ele coord frame.

    y_lim : float
        Transverse corner points in aperture_ele coord frame.

    source_ele : EleStruct
        Photon source element.

    Returns
    -------
    corner : TargetPointStruct
        Corner coords in source_ele ref frame.
    """

def photon_target_setup(ele: _pybmad.EleStruct) -> None:
    """
    Routine to calculate and store the parmeters needed for photon targeting.
    This routine is called by Bmad parsing routines and is not meant for general use.

    Photon initialization with targeting is done by the routine init_photon_from_a_photon_init_ele
    Which is called by init_coord.

    Parameters
    ----------
    ele : EleStruct
        Source element to setup. Element will be of type: sample, diffraction_plate or photon_init.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Source element with target parameters setup in ele.photon.target.
    """

def photon_type(ele: _pybmad.EleStruct) -> int:
    """
    Routine to return the type of photon to be tracked: coherent$ or incoherent$.

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    Returns
    -------
    e_type : int
        coherent$ or incoherent$
    """

def physical_ele_end(track_end: int, orbit: _pybmad.CoordStruct, ele_orientation: int, return_stream_end: bool | None = None) -> int:
    """
    Wrapper for Fortran routine physical_ele_end

    Parameters
    ----------
    track_end : int
        first_track_edge$, second_track_edge$, surface$, or in_between$

    orbit : CoordStruct
        Particle position.

    ele_orientation : int
        Either 1 = Normal or -1 = element reversed.

    return_stream_end : bool, optional
        If True return the stream end instead of the physical end. Default is False.

    Returns
    -------
    physical_end : int
        Return_stream_end ->  Possibilities False             ->  entrance_end$, exit_end$, surface$, or
        in_between$ True              ->  upstream_end$, downstream_end$
    """

def point_photon_emission(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct, direction: int, max_target_area: float, w_to_surface: Sequence[Sequence[float]] | None = None) -> None:
    """
    Routine to emit a photon from a point that may be on a surface.
    If there is a downstream target, the emission calc will take this into account.

    Parameters
    ----------
    ele : EleStruct
        Emitting element.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords of photon. --   Will be in curved surface coords if there is a curved surface.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Final phase-space coords

    direction : int
        +1 -> Emit in forward +z direction, -1 -> emit backwards.

    max_target_area : float
        Area of the solid angle photons may be emitted over. max_target_area is used for normalizing the photon
        field. generally will be equal to twopi or fourpi.

    w_to_surface : 2D array of float (shape: 3,3), optional
        Rotation matrix for curved surface.
    """

class PointerToAttribute:
    """pointer_to_attribute return type"""

    @property
    def a_ptr(self) -> _pybmad.AllPointerStruct: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def ix_attrib(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_attribute(ele: _pybmad.EleStruct, attrib_name: str, do_allocation: bool, err_print_flag: bool | None = None, do_unlink: bool | None = None) -> PointerToAttribute:
    """
    Wrapper for Fortran routine pointer_to_attribute

    Parameters
    ----------
    ele : EleStruct
        After this routine finishes Ptr_attrib will point to a variable within this element.

    attrib_name : str
        Name of attribute. Must be uppercase. For example: "HKICK".

    do_allocation : bool
        If True then do an allocation if needed. EG: The multipole An and Bn arrays need to be allocated before
        their use.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message on error.

    do_unlink : bool, optional
        Default False. If True and applicable, unlink the structure containing the attribute. See above for
        details.

    Returns
    -------
    a_ptr : AllPointerStruct
        Pointer to the attribute.

    err_flag : bool
        Set True if attribtute not found. False otherwise.

    ix_attrib : int, optional
        If applicable, this is the index to the attribute in the ele.value(:), ele.control.var(:), ele.a_pole(:)
        or ele.b_pole(:) arrays. Set to 0 if not in any of these arrays.
    """

@overload
def pointer_to_branch(ele: _pybmad.EleStruct) -> _pybmad.BranchStruct | None:
    """
    Routine to return a pointer to the lattice branch associated with a given name
    or a given element.

    This routine is an overloaded name for:
      pointer_to_branch_given_ele (ele) result (branch_ptr))
      pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

    The lattice branch *associated* with a given element is not necessarily the
    branch where the element is *located*. For example, all lords live in branch #0.
    But the branch associated with a super_lord element is the branch of its slaves.

    To get the branch where the element is located, simply use ele%ix_branch.

    Note: Result is ambiguous if ele argument is associated with multiple branches
    which can happen, for example, with overlay elements.

    Parameters
    ----------
    ele : EleStruct
        Element contained in the branch.

    Returns
    -------
    branch_ptr : BranchStruct, optional
        Pointer to the branch. Nullified if there is no associated branch.
    """

@overload
def pointer_to_branch(branch_name: str, lat: _pybmad.LatStruct, parameter_is_branch0: bool | None = None, blank_branch: int | None = None) -> _pybmad.BranchStruct | None:
    r"""
    Routine to return a pointer to the lattice branch associated with a given name
    or a given element.

    This routine is an overloaded name for:
      pointer_to_branch_given_ele (ele) result (branch_ptr))
      pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

    The lattice branch *associated* with a given element is not necessarily the
    branch where the element is *located*. For example, all lords live in branch #0.
    But the branch associated with a super_lord element is the branch of its slaves.

    To get the branch where the element is located, simply use ele%ix_branch.

    Note: Result is ambiguous if ele argument is associated with multiple branches
    which can happen, for example, with overlay elements.

    Parameters
    ----------
    branch_name : str
        May be a branch name or a branch index.

    lat : LatStruct
        Lattice to search.

    parameter_is_branch0 : bool, optional
        If True, 'PARAMETER' is taken to be an alternative name for branch(0). Default is False.

    blank_branch : int, optional
        Branch index if branch_name = \'\'. Default is blank is an error.

    Returns
    -------
    branch_ptr : BranchStruct, optional
        Pointer to the branch. Nullified if there is no associated branch.
    """

@overload
def pointer_to_ele(lat: _pybmad.LatStruct, ix_ele: int, ix_branch: int | None = None) -> _pybmad.EleStruct | None:
    """
    Routine to return a pointer to an element.
    pointer_to_ele is an overloaded name for:
        Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
        Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
        Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
        Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

    pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
    lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
    element in lat.

    Note that using ele_name to locate an element is potentially dangerous if there
    are multiple elements that have the same name. Better in this case is to use:
      lat_ele_locator

    Also see:
      pointer_to_slave
      pointer_to_lord

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    ix_ele : int
        Index of element in lat.branch(ix_branch).

    ix_branch : int, optional
        Index of the lat.branch(:) containing the element.

    Returns
    -------
    ele_ptr : EleStruct, optional
        Pointer to the element. Nullified if no match or error.
    """

@overload
def pointer_to_ele(lat: _pybmad.LatStruct, ele_loc: _pybmad.LatEleLocStruct) -> _pybmad.EleStruct | None:
    """
    Routine to return a pointer to an element.
    pointer_to_ele is an overloaded name for:
        Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
        Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
        Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
        Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

    pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
    lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
    element in lat.

    Note that using ele_name to locate an element is potentially dangerous if there
    are multiple elements that have the same name. Better in this case is to use:
      lat_ele_locator

    Also see:
      pointer_to_slave
      pointer_to_lord

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    ele_loc : LatEleLocStruct
        Location identification.

    Returns
    -------
    ele_ptr : EleStruct, optional
        Pointer to the element. Nullified if no match or error.
    """

@overload
def pointer_to_ele(lat: _pybmad.LatStruct, ele_name: str) -> _pybmad.EleStruct | None:
    """
    Routine to return a pointer to an element.
    pointer_to_ele is an overloaded name for:
        Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
        Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
        Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
        Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

    pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
    lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
    element in lat.

    Note that using ele_name to locate an element is potentially dangerous if there
    are multiple elements that have the same name. Better in this case is to use:
      lat_ele_locator

    Also see:
      pointer_to_slave
      pointer_to_lord

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    ele_name : str
        Name or index of element.

    Returns
    -------
    ele_ptr : EleStruct, optional
        Pointer to the element. Nullified if no match or error.
    """

@overload
def pointer_to_ele(lat: _pybmad.LatStruct, foreign_ele: _pybmad.EleStruct) -> _pybmad.EleStruct | None:
    """
    Routine to return a pointer to an element.
    pointer_to_ele is an overloaded name for:
        Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
        Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
        Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
        Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

    pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
    lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
    element in lat.

    Note that using ele_name to locate an element is potentially dangerous if there
    are multiple elements that have the same name. Better in this case is to use:
      lat_ele_locator

    Also see:
      pointer_to_slave
      pointer_to_lord

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    foreign_ele : EleStruct
        Lattice element in another lattice.

    Returns
    -------
    ele_ptr : EleStruct, optional
        Pointer to the element. Nullified if no match or error.
    """

class PointerToElementAtS:
    """pointer_to_element_at_s return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def s_eff(self) -> float: ...

    @property
    def position(self) -> _pybmad.CoordStruct: ...

    @property
    def ele(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_element_at_s(branch: _pybmad.BranchStruct, s: float, choose_max: bool, print_err: bool | None = None) -> PointerToElementAtS:
    """
    Function to return a pointer to the element at position s.
    That is, return ele => branch%ele(ix_ele) such that:
    If choose_max = True:
        If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
        Else: branch%ele(ix_ele)%s_strat <= s < branch%ele(ix_ele)%s
    If choose_max = False:
        If s = branch%ele(0): ix_ele = 0
        Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
    That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
        choose_max = True  => ix_ele = ix2
        choose_max = False => ix_ele = ix1

    Also see: element_at_s

    The setting of choose_max only makes a difference when s corresponds to an element boundary.

    Note: For a circular lattice, s is evaluated at the effective s which
    is modulo the branch length:
        s_eff = s - branch_length * floor(s/branch_length)

    Note: If there are multiple elements that are at the given s position due to the presence of
    an element with a negative length, which of the possible elements is actually chosen is ill-defined.

    Parameters
    ----------
    branch : BranchStruct
        Branch to use

    s : float
        Longitudinal position.

    choose_max : bool
        See above.

    print_err : bool, optional
        Print error message if there is an error? Default is True.

    Returns
    -------
    err_flag : bool, optional
        Set True if s is out of bounds. False otherwise.

    s_eff : float, optional
        Effective s. Equal to s with a open lattice. See above.

    position : CoordStruct, optional
        Positional information.

    ele : EleStruct, optional
        Pointer to element at s.
    """

def pointer_to_fibre(ele: _pybmad.EleStruct) -> _pybmad.Fibre | None:
    """
    Wrapper for Fortran routine pointer_to_fibre

    Parameters
    ----------
    ele : EleStruct
        Bmad element

    Returns
    -------
    assoc_fibre : Fibre, optional
        Pointer to the associated fibre.
    """

class PointerToFieldEle:
    """pointer_to_field_ele return type"""

    @property
    def dz_offset(self) -> float: ...

    @property
    def field_ele(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_field_ele(ele: _pybmad.EleStruct, ix_field_ele: int) -> PointerToFieldEle:
    """
    Wrapper for Fortran routine pointer_to_field_ele

    Parameters
    ----------
    ele : EleStruct
        Element with sum number of associated field elements.

    ix_field_ele : int
        Index of the field element to point to. This index runs from 1 to num_field_eles(ele).

    Returns
    -------
    dz_offset : float, optional
        Longitudinal offset of ele upstream edge from the field ele pointed to.

    field_ele : EleStruct, optional
        Pointer to the field element with index ix_field_ele. Will point to null if ix_field_ele is out of range.
    """

class PointerToGirder:
    """pointer_to_girder return type"""

    @property
    def ix_slave_back(self) -> int: ...

    @property
    def girder(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_girder(ele: _pybmad.EleStruct) -> PointerToGirder:
    """
    Wrapper for Fortran routine pointer_to_girder

    Parameters
    ----------
    ele : EleStruct
        Element to check.

    Returns
    -------
    ix_slave_back : int, optional
        Index back to ele. That is, pointer_to_slave(girder, ix_slave_back) will point back to ele. Set to -1 if
        no girder present

    girder : EleStruct, optional
        : Pointer to the girder. Null if ele is not girder supported.
    """

class PointerToIndexedAttribute:
    """pointer_to_indexed_attribute return type"""

    @property
    def a_ptr(self) -> _pybmad.AllPointerStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_indexed_attribute(ele: _pybmad.EleStruct, ix_attrib: int, do_allocation: bool, err_print_flag: bool | None = None) -> PointerToIndexedAttribute:
    """
    Wrapper for Fortran routine pointer_to_indexed_attribute

    Parameters
    ----------
    ele : EleStruct
        After this routine finishes A_ptr will point to a variable within this element.

    ix_attrib : int
        Integer, Attribute index.

    do_allocation : bool
        If True then do an allocation if needed. EG: The multipole An and Bn arrays need to be allocated before
        their use.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message on error.

    Returns
    -------
    a_ptr : AllPointerStruct
        Pointer to the attribute.

    err_flag : bool
        Set True if attribtute not found. False otherwise.
    """

class PointerToLord:
    """pointer_to_lord return type"""

    @property
    def control(self) -> _pybmad.ControlStruct | None: ...

    @property
    def ix_slave_back(self) -> int: ...

    @property
    def ix_control(self) -> int: ...

    @property
    def ix_ic(self) -> int: ...

    @property
    def lord_ptr(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_lord(slave: _pybmad.EleStruct, ix_lord: int, lord_type: int | None = None) -> PointerToLord:
    """
    Wrapper for Fortran routine pointer_to_lord

    Parameters
    ----------
    slave : EleStruct
        Slave element.

    ix_lord : int
        Index of the lord.

    lord_type : int, optional
        See above.

    Returns
    -------
    control : ControlStruct, optional
        Pointer to control info for this lord/slave relationship. Nullified if there is an error.

    ix_slave_back : int, optional
        Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) will point back to slave. Set
        to -1 if there is an error or the slave is a slice_slave.

    ix_control : int, optional
        Index in lat.control(:) array the control argument is at. For ramper lord elements, ix_control is index
        for the lord.control.ramper(:) array.

    ix_ic : int, optional
        Index of the lat.ic(:) element associated with the control argument.

    lord_ptr : EleStruct, optional
        Pointer to the lord. Nullified if there is an error.
    """

class PointerToMultipassLord:
    """pointer_to_multipass_lord return type"""

    @property
    def ix_pass(self) -> int: ...

    @property
    def super_lord(self) -> _pybmad.EleStruct | None: ...

    @property
    def multi_lord(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_multipass_lord(ele: _pybmad.EleStruct) -> PointerToMultipassLord:
    """
    Wrapper for Fortran routine pointer_to_multipass_lord

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    Returns
    -------
    ix_pass : int, optional
        Multipass turn number. Set to 0 if element is a multipass_lord. Set to -1 if element is not a
        multipass_slave.

    super_lord : EleStruct, optional
        super_lord of the element. Set to NULL if ele is not a super_slave or super_lord. Note: if ele is a
        multipass_lord there are multiple possible super_lord slaves.

    multi_lord : EleStruct, optional
        multipass_lord if there is one. Set to NULL if there is no multipass_lord.
    """

def pointer_to_next_ele(this_ele: _pybmad.EleStruct, offset: int | None = None, skip_beginning: bool | None = None, follow_fork: bool | None = None, ix_multipass: int | None = None) -> _pybmad.EleStruct | None:
    """
    Wrapper for Fortran routine pointer_to_next_ele

    Parameters
    ----------
    this_ele : EleStruct
        Starting element.

    offset : int, optional
        +1 -> return next element, +2 -> element after that, etc. Can be negative. Default = +1.

    skip_beginning : bool, optional
        If True then skip beginning element #0 when wrapping around. Default is False.

    follow_fork : bool, optional
        If True then fork at any fork element. Default is False.

    ix_multipass : int, optional
        Default = 1. Used to select the multipass branch if this_ele is a multipass_lord.

    Returns
    -------
    next_ele : EleStruct, optional
        Element after this_ele (if offset = 1). Nullified if there is an error. EG bad this_ele.
    """

class PointerToSlave:
    """pointer_to_slave return type"""

    @property
    def control(self) -> _pybmad.ControlStruct | None: ...

    @property
    def ix_lord_back(self) -> int: ...

    @property
    def ix_control(self) -> int: ...

    @property
    def ix_ic(self) -> int: ...

    @property
    def slave_ptr(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_slave(lord: _pybmad.EleStruct, ix_slave: int, slave_type: int | None = None) -> PointerToSlave:
    """
    Function to point to a slave of a lord.
    Note: Ramper lords do not have any associated slaves (slaves are assigned dynamically at run time).

    If slave_type = all$ (the default) the range for ix_slave is:
      1 to lord%n_slave                                 for "regular" slaves.
      lord%n_slave+1 to lord%n_slave+lord%n_slave_field for field overlap slaves.

    If slave_type = field_slave$, only the field overlap slaves may be accessed and the range for ix_slave is:
      1 to lord%n_slave_field

    Also see:
      pointer_to_lord
      pointer_to_super_lord
      pointer_to_ele
      num_lords

    Parameters
    ----------
    lord : EleStruct
        Lord element

    ix_slave : int
        Index of the slave in the list of slaves controled by the lord..

    slave_type : int, optional
        See above.

    Returns
    -------
    control : ControlStruct, optional
        Pointer to control info for this lord/slave relationship. Nullified if there is an error.

    ix_lord_back : int, optional
        Index back to the lord. That is, pointer_to_lord(slave_ptr, ix_lord_back) will point back to the lord. Set
        to -1 if there is an error.

    ix_control : int, optional
        Index in lat.control(:) array the control argument is at.

    ix_ic : int, optional
        Index of the lat.ic(:) element associated with the control argument.

    slave_ptr : EleStruct, optional
        Pointer to the slave. Nullified if there is an error.
    """

class PointerToSuperLord:
    """pointer_to_super_lord return type"""

    @property
    def control(self) -> _pybmad.ControlStruct | None: ...

    @property
    def ix_slave_back(self) -> int: ...

    @property
    def ix_control(self) -> int: ...

    @property
    def ix_ic(self) -> int: ...

    @property
    def lord_ptr(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_super_lord(slave: _pybmad.EleStruct, lord_type: int | None = None) -> PointerToSuperLord:
    """
    Wrapper for Fortran routine pointer_to_super_lord

    Parameters
    ----------
    slave : EleStruct
        Slave element.

    lord_type : int, optional
        If present, only return a super_lord of this type.

    Returns
    -------
    control : ControlStruct, optional
        Pointer to control info for this lord/slave relationship. Nullified if there is an error.

    ix_slave_back : int, optional
        Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) will point back to slave. Set
        to -1 if there is an error or the slave is a slice_slave.

    ix_control : int, optional
        Index in lat.control(:) array the control argument is at. For ramper lord elements, ix_control is index
        for the lord.control.ramper(:) array.

    ix_ic : int, optional
        Index of the lat.ic(:) element associated with the control argument.

    lord_ptr : EleStruct, optional
        Pointer to the lord.
    """

class PointerToSurfaceDisplacementPt:
    """pointer_to_surface_displacement_pt return type"""

    @property
    def ix(self) -> int: ...

    @property
    def iy(self) -> int: ...

    @property
    def xx(self) -> float: ...

    @property
    def yy(self) -> float: ...

    @property
    def pt(self) -> _pybmad.SurfaceDisplacementPtStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_surface_displacement_pt(ele: _pybmad.EleStruct, nearest: bool, x: float, y: float, extend_grid: bool | None = None) -> PointerToSurfaceDisplacementPt:
    """
    Routine to point to the grid point struct associated with point (x,y).

    Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.

    Parameters
    ----------
    ele : EleStruct
        Element containing the grid

    nearest : bool
        If True, return pointer to nearest grid point. If False, return pointer to the grid point lower and left
        of (x,y).

    x : float
        Photon position.

    y : float
        Photon position.

    extend_grid : bool, optional
        If (x,y) past grid pretend (x,y) is at grid boundary. Default is False.

    Returns
    -------
    ix : int, optional
        Grid point index.

    iy : int, optional
        Grid point index.

    xx : float, optional
        Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
        the nearest grid boundary point.

    yy : float, optional
        Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
        the nearest grid boundary point.

    pt : SurfaceDisplacementPtStruct, optional
        Pointer to grid point. Will not be associated if (x,y) outside the grid.
    """

class PointerToSurfaceSegmentedPt:
    """pointer_to_surface_segmented_pt return type"""

    @property
    def ix(self) -> int: ...

    @property
    def iy(self) -> int: ...

    @property
    def xx(self) -> float: ...

    @property
    def yy(self) -> float: ...

    @property
    def pt(self) -> _pybmad.SurfaceSegmentedPtStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_surface_segmented_pt(ele: _pybmad.EleStruct, nearest: bool, x: float, y: float, extend_grid: bool | None = None) -> PointerToSurfaceSegmentedPt:
    """
    Routine to point to the grid point struct associated with point (x,y).

    Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.

    Parameters
    ----------
    ele : EleStruct
        Element containing the grid

    nearest : bool
        If True, return pointer to nearest grid point. If False, return pointer to the grid point lower and left
        of (x,y).

    x : float
        Photon position.

    y : float
        Photon position.

    extend_grid : bool, optional
        If (x,y) past grid pretend (x,y) is at grid boundary. Default is False.

    Returns
    -------
    ix : int, optional
        Grid point index.

    iy : int, optional
        Grid point index.

    xx : float, optional
        Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
        the nearest grid boundary point.

    yy : float, optional
        Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set to be on
        the nearest grid boundary point.

    pt : SurfaceSegmentedPtStruct, optional
        Pointer to grid point. Will not be associated if (x,y) outside the grid.
    """

class PointerToWakeEle:
    """pointer_to_wake_ele return type"""

    @property
    def delta_s(self) -> float: ...

    @property
    def wake_ele(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_wake_ele(ele: _pybmad.EleStruct) -> PointerToWakeEle:
    """
    Wrapper for Fortran routine pointer_to_wake_ele

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    Returns
    -------
    delta_s : float, optional
        distance of wake locaiton from beginning of ele.

    wake_ele : EleStruct, optional
        Element having the associated wake. wake_ele will be nullified if there is no associated wake.
    """

class PointerToWall3d:
    """pointer_to_wall3d return type"""

    @property
    def ds_offset(self) -> float: ...

    @property
    def is_branch_wall(self) -> bool: ...

    @property
    def wall3d(self) -> _pybmad.Wall3DStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointer_to_wall3d(ele: _pybmad.EleStruct, ix_wall: int | None = None) -> PointerToWall3d:
    """
    Function to return a pointer to a wall3d structure associated
    with a given lattice element.

    Note: The wall associated with a the vacuum chamber is the branch%wall3d.

    Parameters
    ----------
    ele : EleStruct
        lattice element.

    ix_wall : int, optional
        index in wall3d(:) array. Default is 1.

    Returns
    -------
    ds_offset : float, optional
        Element offset: s(beginning of ele) - s(beginning of wall3d)

    is_branch_wall : bool, optional
        Set True if wall3d points to branch.wall3d.

    wall3d : Wall3dStruct, optional
        Pointer to the associated wall structure. Will be nullified if there is no associated wall.
    """

class PointersToAttribute:
    """pointers_to_attribute return type"""

    @property
    def ptr_array(self) -> _pybmad.AllPointerStructAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def eles(self) -> _pybmad.ElePointerStructAlloc1D: ...

    @property
    def ix_attrib(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def pointers_to_attribute(lat: _pybmad.LatStruct, ele_name: str, attrib_name: str, do_allocation: bool, err_print_flag: bool | None = None, do_unlink: bool | None = None) -> PointersToAttribute:
    """
    Wrapper for Fortran routine pointers_to_attribute

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    ele_name : str
        Element name. Must be uppercase

    attrib_name : str
        Attribute name. Must be uppercase. For example: "HKICK".

    do_allocation : bool
        If True then do an allocation if needed. EG: The multipole An and Bn arrays need to be allocated before
        their use.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message on error.

    do_unlink : bool, optional

    Returns
    -------
    ptr_array : 1D array of AllPointerStruct
        Pointer to the attribute. Size of ptr_array will be set to 0 if there is a problem.

    err_flag : bool
        Set True if attribtute not found.

    eles : 1D array of ElePointerStruct, optional
        Array of element pointers. size(eles) = size(ptr_array). If there are no associated lattice elements (EG
        if ele_name = 'PARTICLE_START'), eles(i).ele will be null.

    ix_attrib : int, optional
        If applicable then this is the index to the attribute in the ele.value(:) array.
    """

def polar_to_spinor(polar: _pybmad.SpinPolarStruct) -> list[complex]:
    """
    Wrapper for Fortran routine polar_to_spinor

    Parameters
    ----------
    polar : SpinPolarStruct
        includes polar phase

    Returns
    -------
    spinor : 1D array of complex (shape: 2)
        Spinor
    """

def polar_to_vec(polar: _pybmad.SpinPolarStruct) -> list[float]:
    """
    Wrapper for Fortran routine polar_to_vec

    Parameters
    ----------
    polar : SpinPolarStruct
        Spin_polar_struct

    Returns
    -------
    vec : 1D array of float (shape: 3)
        Real(3)
    """

def print_mesh3d(mesh3d: _pybmad.Mesh3DStruct) -> None:
    """
    Wrapper for Fortran routine print_mesh3d

    Parameters
    ----------
    mesh3d : Mesh3dStruct
    """

def prob_x_diffuse(x: float, d_param: _pybmad.DiffuseParamStruct, surface: _pybmad.PhotonReflectSurfaceStruct) -> float:
    """
    Contained routine to calculate integrated probability distribution in x = sin(graze_angle_out).

    Parameters
    ----------
    x : float
        sin(graze_angle_out)
    """

class ProjectEmitToXyz:
    """project_emit_to_xyz return type"""

    @property
    def sigma_x(self) -> float: ...

    @property
    def sigma_y(self) -> float: ...

    @property
    def sigma_z(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def project_emit_to_xyz(ring: _pybmad.LatStruct, ix: int, mode: _pybmad.NormalModesStruct) -> ProjectEmitToXyz:
    """
    Obtains the projected x, y, and z beamsizes by building the sigma matrix
    from the normal mode emittances and 1-turn transfer matrix.
    These projectes beamsize are what would be seen by instrumentation.

    This method of projecting takes into account transverse and longitudinal coupling.

    This method of obtaining the projected beam sizes is from "Alternitive approach to general
    coupled linear optics" by Andrzej Wolski.

    The normal mode emittances used to generate a beam envelop sigma matrix from the
    1-turn transfer matrix.  The projected sizes are from the 1, 1 3, 3 and 5, 5 elements of
    the sigma matrix.

    Parameters
    ----------
    ring : LatStruct
        the storage ring

    ix : int
        element at which to make the projection

    mode : NormalModesStruct
        normal mode emittances

    Returns
    -------
    sigma_x : float
        projected horizontal beamsize

    sigma_y : float
        projected vertical beamsize

    sigma_z : float
        projected longitudinal beamsize
    """

def propagate_part_way(orb_start: _pybmad.CoordStruct, param: _pybmad.LatParamStruct, pt: _pybmad.RadIntTrackPointStruct, info: _pybmad.RadIntInfoStruct, z_here: float, runt: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine propagate_part_way

    Parameters
    ----------
    orb_start : CoordStruct

    param : LatParamStruct

    pt : RadIntTrackPointStruct

    info : RadIntInfoStruct

    z_here : float

    runt : EleStruct
    """

def psi_prime_sca(t: float, p: float, args: Sequence[float]) -> float:
    """
    This wraps the array-valued psi_prime function as a scalar.

    See psi_prime comments for details.

    Parameters
    ----------
    t : float
        time relative to RF bucket

    p : float
        psi(t)

    args : 1D array of float (shape: 1:8)
        parameters and constants of DEQ

    Returns
    -------
    dpdt : float
        dpsi_dt
    """

def ptc_bookkeeper(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine ptc_bookkeeper

    Parameters
    ----------
    lat : LatStruct
        Bmad lattice.
    """

def ptc_calculate_tracking_step_size(ptc_layout: _pybmad.Layout, kl_max: float, ds_max: float | None = None, even_steps: _pybmad.BoolAlloc1D | None = None, r_typical: float | None = None, dx_tol_bend: float | None = None, use_2nd_order: bool | None = None, crossover: Sequence[int] | None = None, crossover_wiggler: Sequence[int] | None = None) -> None:
    """
    even_steps, r_typical, dx_tol_bend, use_2nd_order)

    Routine to calculate the optimum number of tracking steps and order
    of the integrator (2, 4, or 6) for each fibre in a layout.

    See the Bmad manual chapter on PTC for more details.

    Parameters
    ----------
    ptc_layout : Layout
        This parameter is an input/output and is modified in-place.
        As an output, ptc_layout: Lattice with the optimum number of tracking steps and integrator order.

    kl_max : float
        Maximum K1*L per tracking step.

    ds_max : float, optional
        Maximum ds for any step. Useful when including other physicas like space charge.

    even_steps : 1D array of bool (shape: 2), optional
        Always use an even number of steps for a fibre? Useful if need to evaluate at the center of fibres.

    r_typical : float, optional
        Typical transverse offset. Used for computing the effective contribution of K1*L due to sextupoles.

    dx_tol_bend : float, optional
        Tolerable residual orbit in a bend.

    use_2nd_order : bool, optional
        If present and True then force the use of 2nd order integrator.

    crossover : 1D array of int (shape: 2), optional
        crossover points between orders for all elements except wigglers. Default is [4, 18].

    crossover_wiggler : 1D array of int (shape: 2), optional
        crossover points for wigglers. Default is [30, 60].
    """

class PtcCheckForLostParticle:
    """ptc_check_for_lost_particle return type"""

    @property
    def state(self) -> int: ...

    @property
    def ptc_fibre(self) -> _pybmad.Fibre | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_check_for_lost_particle(do_reset: bool) -> PtcCheckForLostParticle:
    """
    Routine to check if a particle has been lost when tracking with PTC.

    Parameters
    ----------
    do_reset : bool
        If True then reset ptc flags.

    Returns
    -------
    state : int
        Same as coord_struct.state. alive$, lost$, lost_neg_x$, etc.

    ptc_fibre : Fibre, optional
        Pointer to fibre where particle lost. Nullified if particle alive.
    """

def ptc_closed_orbit_calc(branch: _pybmad.BranchStruct, radiation_damping_on: bool | None = None) -> _pybmad.CoordStructAlloc1D:
    """
    Routine to calculate the closed orbit of a lattice branch using PTC.
    This routine assumes the associated PTC layout has been crated
    with lat_to_ptc_layout.

    Parameters
    ----------
    branch : BranchStruct
        Branch of a lattice.

    radiation_damping_on : bool, optional
        If True, radiation dampling is included in the calculation. Default is the setting of
        bmad_com..radiation_damping_on.

    Returns
    -------
    closed_orbit : 1D array of CoordStruct
        closed_orbit
    """

class PtcEmitCalc:
    """ptc_emit_calc return type"""

    @property
    def norm_mode(self) -> _pybmad.NormalModesStruct: ...

    @property
    def closed_orb(self) -> _pybmad.CoordStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_emit_calc(ele: _pybmad.EleStruct, sigma_mat: Sequence[Sequence[float]]) -> PtcEmitCalc:
    """
    Routine to calculate emittances, etc.

    Note: This routine calls the PTC init_all routine.

    Parameters
    ----------
    ele : EleStruct
        Element at which to evaluate the parameters.

    Returns
    -------
    norm_mode : NormalModesStruct
        Normal_modes_struct %a%tune, %b%tune, %z%tune %a%alpha_damp, etc. %a%emittance, etc.

    closed_orb : CoordStruct
        Closed orbit at ele (Bmad coordinates). Notice: This closed orbit is the closed orbit with radiation on.
    """

def ptc_kill_map_with_radiation(rad_map: _pybmad.PtcRadMapStruct) -> None:
    """
    Routine to kill a binary file containing a ptc_rad_map_struct map

    Parameters
    ----------
    rad_map : PtcRadMapStruct
        Map with radiation included.
        This parameter is an input/output and is modified in-place.
        As an output, rad_map: Deallocated map.
    """

def ptc_layouts_resplit(dKL_max: float, l_max: float, l_max_drift_only: bool, bend_dorb: float, sex_dx: float, even: bool | None = None, crossover: Sequence[int] | None = None, crossover_wiggler: Sequence[int] | None = None) -> None:
    """
    even, crossover, crossover_wiggler)

    Routine to resplit (that is, recalculate the number of integration steps for an element)
    For the fibres in all layouts. After doing a resplit, the tune (and any other relavent
    "adjustable" parameters) should be adjusted to the correct values.

    Parameters
    ----------
    dKL_max : float
        Maximum K1 * L quadrupole strength allowed for an integration step. Reasonable value would be something
        like 0.04.

    l_max : float
        Maximum step length. Ignored if set to 0.

    l_max_drift_only : bool
        If True then l_max is only used for splitting drifts.

    bend_dorb : float
        Residual bend orbit error. With some integration methods a zero orbit at the start of the bend will not be
        zero at the end. In this case, bend_dorb sets a maximum allowable orbit deviation. If set to zero, this
        argument will be ignored. A resonable value is 10d-7. Note that the actual orbit deviation is not simply
        related to bend_dorb and can be larger. In any case, lowering bend_dorb (without making it zero) will
        lower the

    sex_dx : float
        To split sextupoles, sex_dx is used as the reference position about which the quadrupole strength is
        calculated. This quadrupole strength is then used with dKL_max to calculate the number of integration
        steps. Set to zero to ignore.

    even : bool, optional
        If True then each fibre  will have an even number of steps. If False then the number of steps will be odd.
        If not present then number of steps is not constrained to be even or odd.

    crossover : 1D array of int (shape: 2), optional
        crossover(1) sets the maximum number of 2nd order integration steps to use. If the number of steps would
        exceed crossover(1) then integration is switched to 4th order. crossover(2) sets the maximum number of 4th
        order integration steps. If this number is exceeded, 6th order integration is used. Currently the default
        in PTC is [4, 18].

    crossover_wiggler : 1D array of int (shape: 2), optional
        crossover for wiggler elements.
    """

def ptc_linear_isf_calc(branch: _pybmad.BranchStruct) -> _pybmad.LinearEleIsfStructAlloc1D:
    """
    Wrapper for Fortran routine ptc_linear_isf_calc

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch to analyze.

    Returns
    -------
    ele_isf : 1D array of LinearEleIsfStruct
        ISF at every element.
    """

def ptc_one_turn_mat_and_closed_orbit_calc(branch: _pybmad.BranchStruct, pz: float | None = None) -> None:
    """
    Routine to compute the transfer matrices for the individual elements and closed orbit
    for a lattice branch with closed geometry.

    Note: PTC itself does not compute Twiss parameters. Use twiss_from_mat6 to compute this.

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Lattice branch containing the matrices.

    pz : float, optional
        energy offset around which to calculate the matrices if there is no RF.
    """

def ptc_ran_seed_put(iseed: int) -> None:
    """
    Wrapper for Fortran routine ptc_ran_seed_put

    Parameters
    ----------
    iseed : int
        0 -> Use system clock.
    """

class PtcReadFlatFile:
    """ptc_read_flat_file return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def lat(self) -> _pybmad.LatStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_read_flat_file(flat_file: _pybmad.CharacterAlloc1D, create_end_marker: bool | None = None, from_mad: bool | None = None) -> PtcReadFlatFile:
    """
    Wrapper for Fortran routine ptc_read_flat_file

    Parameters
    ----------
    flat_file : 1D array of str
        Name(s) of PTC flat file(s).

    create_end_marker : bool, optional
        Put a marker element named END at the end of the lattice brances? Default is True.

    from_mad : bool, optional
        If True, ignore PTC specific parameters like integrator_order. Default is False. True is used when the
        fibre has been created via MAD. In this case, the PTC specific parameters may not have good values.

    Returns
    -------
    err_flag : bool
        Set True if there is a problem.

    lat : LatStruct, optional
        If present then setup a Bmad lattice.
    """

class PtcReadMapWithRadiation:
    """ptc_read_map_with_radiation return type"""

    @property
    def rad_map(self) -> _pybmad.PtcRadMapStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_read_map_with_radiation(file_name: str | None = None, file_unit: int | None = None) -> PtcReadMapWithRadiation:
    """
    Routine to read a binary file containing a ptc_rad_map_struct map

    Either file_name or file_unit must be present but not both.
    File_unit is used when there are multiple maps in a file.
    If file_unit is present, it is the responsibility of the calling routine to open the file beforehand
    and to close the file afterwards.

    Parameters
    ----------
    file_name : str, optional
        Name of binary file.

    file_unit : int, optional
        File unit number read from.

    Returns
    -------
    rad_map : PtcRadMapStruct
        Map with radiation included.

    err_flag : bool
        Set True if there is a read error.
    """

def ptc_set_rf_state_for_c_normal(nocavity: bool) -> None:
    """
    Wrapper for Fortran routine ptc_set_rf_state_for_c_normal

    Parameters
    ----------
    nocavity : bool
        True -> RF is off and vice versa.
    """

def ptc_set_taylor_order_if_needed() -> None:
    """
    Routine to see if the taylor_order for PTC needs to be set/changed.
    For example, for a change in bmad_com%taylor_order.
    """

class PtcSetupMapWithRadiation:
    """ptc_setup_map_with_radiation return type"""

    @property
    def rad_map(self) -> _pybmad.PtcRadMapStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_setup_map_with_radiation(ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct | None = None, map_order: int | None = None, include_damping: bool | None = None, create_symplectic_map: bool | None = None, orbit1: _pybmad.CoordStruct | None = None) -> PtcSetupMapWithRadiation:
    """
    create_symplectic_map, orbit1, err_flag)

    Routine to construct a map including radiation damping and excitation.
    Note: The setting of bmad_com%radiation_damping_on will determine if damping is included in the map.

    ele1/ele2 must have an associated PTC layout (which can be constructed by calling lat_to_ptc_layout).

    To track after calling this routine track by calling ptc_track_with_radiation.
    To cleanup memory after using, call ptc_kill_map_with_radiation.
    To save a map call ptc_write_map_with_radiation.
    To read a saved map call ptc_read_map_with_radiation.
    To set the random number seed call: ptc_ran_seed_put.

    Parameters
    ----------
    ele1 : EleStruct
        The map starts at the exit end of ele1.

    ele2 : EleStruct, optional
        The map ends at the exit end of ele2. If not present, the 1-turn map will be constructed.

    map_order : int, optional
        Order of the map. If not present or less than 1, the currently set order is used.

    include_damping : bool, optional
        If True (default), the map will be constructed with radiation damping included. If False, the map will not
        be constructed with radiation dampling included.

    create_symplectic_map : bool, optional
        If False (default), create a Taylor map. If True, create a partially inverted map which can be
        symplecitally tracked.

    orbit1 : CoordStruct, optional
        Orbit at ele1 about which the map is constructed. If not present then the orbit will be computed using PTC
        tracking.

    Returns
    -------
    rad_map : PtcRadMapStruct
        Transport map.

    err_flag : bool, optional
        Set True if there is an error such as not associated PTC layout.
    """

class PtcSpinCalc:
    """ptc_spin_calc return type"""

    @property
    def norm_mode(self) -> _pybmad.NormalModesStruct: ...

    @property
    def closed_orb(self) -> _pybmad.CoordStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_spin_calc(ele: _pybmad.EleStruct, sigma_mat: Sequence[Sequence[float]]) -> PtcSpinCalc:
    """
    Routine to equilibrium polarizations, etc.

    Parameters
    ----------
    ele : EleStruct
        Element at which to evaluate the parameters.

    Returns
    -------
    norm_mode : NormalModesStruct
        Normal_modes_struct %a%tune, %b%tune, %z%tune %a%alpha_damp, etc. %a%emittance, etc.

    closed_orb : CoordStruct
        Closed orbit at ele (Bmad coordinates). Notice: This closed orbit is the closed orbit with radiation on.
    """

def ptc_spin_matching_calc(branch: _pybmad.BranchStruct) -> _pybmad.SpinMatchingStructAlloc1D:
    """
    Wrapper for Fortran routine ptc_spin_matching_calc

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch to analyze.

    Returns
    -------
    match_info : 1D array of SpinMatchingStruct
        G-matrix and other stuff. The array will be allocated by this routine.
    """

class PtcTrackAll:
    """ptc_track_all return type"""

    @property
    def track_state(self) -> int: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ptc_track_all(branch: _pybmad.BranchStruct, orbit: _pybmad.CoordStructAlloc1D) -> PtcTrackAll:
    """
    Routine to track from the start to the end of a lattice branch.

    Parameters
    ----------
    branch : BranchStruct
        Lat to track through.

    orbit : 1D array of CoordStruct
        Coordinates at beginning of branch.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit array.

    Returns
    -------
    track_state : int, optional
        Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.

    err_flag : bool, optional
        Set true if particle lost or error. False otherwise
    """

def ptc_track_map_with_radiation(orbit: _pybmad.CoordStruct, rad_map: _pybmad.PtcRadMapStruct, rad_damp: bool | None = None, rad_fluct: bool | None = None) -> None:
    """
    Routine to track through a map that includes radiation.

    NOTE! Tracking without damping when the map was made with radiation (and vice versa)
    will not give good results. So avoid this situation unless testing.

    To construct the map, use the routine ptc_setup_map_with_radiation.
    To cleanup memory after using, call ptc_kill_map_with_radiation.
    To save a map call ptc_write_map_with_radiation.
    To read a saved map call ptc_read_map_with_radiation.
    To set the random number seed call: ptc_ran_seed_put.

    Parameters
    ----------
    orbit : CoordStruct
        Starting orbit.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending orbit after tracking through the map.

    rad_map : PtcRadMapStruct
        Map with radiation included.

    rad_damp : bool, optional
        Override the setting of bmad_com.radiation_damping_on.

    rad_fluct : bool, optional
        Override the setting of bmad_com.radiation_fluctuations_on
    """

def ptc_transfer_map_with_spin(branch: _pybmad.BranchStruct, t_map: _pybmad.TaylorStructArray1D, s_map: _pybmad.TaylorStructArray1D, orb0: _pybmad.CoordStruct, ix1: int | None = None, ix2: int | None = None, one_turn: bool | None = None, unit_start: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine ptc_transfer_map_with_spin

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch used in the calculation.

    t_map : 1D array of TaylorStruct (shape: 6)
        Initial orbital map (used when unit_start = False)
        This parameter is an input/output and is modified in-place.
        As an output, t_map: Orbital transfer map.

    s_map : 1D array of TaylorStruct (shape: 4)
        Initial spin map (used when unit_start = False)
        This parameter is an input/output and is modified in-place.
        As an output, s_map: Quaternion spin transfer map.

    orb0 : CoordStruct
        Initial orbit around which the map is made.

    ix1 : int, optional
        Element start index for the calculation. Default is 0.

    ix2 : int, optional
        Element end index for the calculation. Default is branch.n_ele_track.

    one_turn : bool, optional
        If present and True, and if ix1 = ix2, and the lattice is circular, then construct the one-turn map from
        ix1 back to ix1. Default = False.

    unit_start : bool, optional
        If present and False then t_map will be used as the starting map instead of the unit map. Default = True

    Returns
    -------
    err_flag : bool
        Set True if problem like number overflow, etc.
    """

def ptc_write_map_with_radiation(rad_map: _pybmad.PtcRadMapStruct, file_name: str | None = None, file_unit: int | None = None) -> None:
    """
    Routine to create or append to a binary file containing a ptc_rad_map_struct map.

    Either file_name or file_unit must be present but not both.
    If file_unit is present, it is the responsibility of the calling routine to open the file beforehand
    and to close the file afterwards.

    Parameters
    ----------
    rad_map : PtcRadMapStruct
        Map with radiation included.

    file_name : str, optional
        Name of binary file to create.

    file_unit : int, optional
        File unit number to append to.
    """

def ptwo(sigma: float, t: float, phi: float, d_param: _pybmad.DiffuseParamStruct) -> float:
    """
    unnormalized two-dimensional probability distribution in x and phi
    polar angles relative to surface normal
    azimuthal angle relative to plane of incidence (plane of incoming ray and surface normal)
    1/y suppressed

    Private routine.
    """

def pwd_mat(lat: _pybmad.LatStruct, t6: Sequence[Sequence[float]], inductance: float, sig_z: float) -> list[list[float]]:
    """
    Calculates potential well distortion as RF defocusing.  Calculates t6_pwd=t6.Mpwd,
    where Mpwd is identity with 65 element proportional to the inductance.

    Vpwd = -inductance * lat%param%n_part * e_charge * c_light**3 / SQRT(twopi) / sig_z**3 / omega_RF  !effective RF voltage from PWD
    Mpwd(6,5) = omega_RF * Vpwd / c_light / lat%ele(0)%value(E_TOT$) * branch%ele(i)%value(l$) / lat%param%total_length

    Parameters
    ----------
    lat : LatStruct
        TYPE(lat_struct)

    t6 : 2D array of float (shape: 6,6)
        1-turn transfer matrix

    inductance : float
        Longitudinal inductance in Henrys.  Something on the order of nH.

    sig_z : float
        Bunch length.

    Returns
    -------
    t6_pwd : 2D array of float (shape: 6,6)
        1-turn transfer matrix with PWD defocusing applied
    """

class Rad1DampAndStocMats:
    """rad1_damp_and_stoc_mats return type"""

    @property
    def rad_map(self) -> _pybmad.RadMapStruct: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def rad_int1(self) -> _pybmad.RadInt1Struct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def rad1_damp_and_stoc_mats(ele: _pybmad.EleStruct, include_opening_angle: bool, orb_in: _pybmad.CoordStruct, orb_out: _pybmad.CoordStruct, g2_tol: float, g3_tol: float, ele0: _pybmad.EleStruct | None = None) -> Rad1DampAndStocMats:
    """
    Routine to calculate the damping and stochastic matrices for a given lattice element.

    Parameters
    ----------
    ele : EleStruct
        Element under consideration.

    include_opening_angle : bool
        If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
        comparing against other codes.

    orb_in : CoordStruct
        Entrance orbit about which to compute the matrices.

    orb_out : CoordStruct
        Exit orbit.

    g2_tol : float
        Tollerance on g^2 per unit length (damping tolerance).

    g3_tol : float
        Tollerance on g^3 per unit length (stocastic tolerance).

    ele0 : EleStruct, optional
        Element before `ele`. Needed if and only if rad_int1 is present

    Returns
    -------
    rad_map : RadMapStruct
        Damping and stochastic matrices.

    err_flag : bool
        Set true if there is an error. False otherwise.

    rad_int1 : RadInt1Struct, optional
        Radiation integrals
    """

class RadDampAndStocMats:
    """rad_damp_and_stoc_mats return type"""

    @property
    def rmap(self) -> _pybmad.RadMapStruct: ...

    @property
    def mode(self) -> _pybmad.NormalModesStruct: ...

    @property
    def xfer_nodamp_mat(self) -> list[list[float]]: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def rad_int_branch(self) -> _pybmad.RadIntBranchStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def rad_damp_and_stoc_mats(ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct, include_opening_angle: bool, closed_orbit: _pybmad.CoordStructArray1D | None = None) -> RadDampAndStocMats:
    """
    Routine to calculate the damping and stochastic variance matrices from exit end of ele1
    to the exit end of ele2. Use ele1 = ele2 to get 1-turn matrices.

    If ele2 is before ele1 the integration range if from ele1 to the branch end plus
    from the beginning to ele2.

    Note: The ele%mat6 matrices will be remade. By convention, these matrices
    do not include damping.

    Parameters
    ----------
    ele1 : EleStruct
        Start element of integration range.

    ele2 : EleStruct
        End element of integration range.

    include_opening_angle : bool
        If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
        comparing against other codes.

    closed_orbit : 1D array of CoordStruct, optional
        Closed orbit. If not present this routine will calculate it.

    Returns
    -------
    rmap : RadMapStruct
        Damping and stochastic mats

    mode : NormalModesStruct

    xfer_nodamp_mat : 2D array of float (shape: 6,6)
        Transfer matrix without damping.

    err_flag : bool
        Set true if there is a problem.

    rad_int_branch : RadIntBranchStruct, optional
        Array of element-by-element radiation integrals.
    """

class RadGIntegrals:
    """rad_g_integrals return type"""

    @property
    def int_g(self) -> list[float]: ...

    @property
    def int_g3(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def rad_g_integrals(ele: _pybmad.EleStruct, where: int, orb_in: _pybmad.CoordStruct, orb_out: _pybmad.CoordStruct, int_g2: float, g_tol: float, g2_tol: float, g3_tol: float) -> RadGIntegrals:
    """
    Routine to calculate bending strength integrals (g(s) = 1/trajectory_bending_radius(s)) in
    laboratory coords.

    Parameters
    ----------
    ele : EleStruct
        Element under consideration.

    where : int
        What part of ele to integrate over. upstream$ -> 1st half of element, downsteam$ -> 2nd half, all$ ->
        everything.

    orb_in : CoordStruct
        Entrance orbit about which to compute the matrices.

    orb_out : CoordStruct
        Exit orbit.

    g_tol : float
        Tollerance on |g| per unit length.

    g2_tol : float
        Tollerance on g^2 per unit length.

    g3_tol : float
        Tollerance on g^3 per unit length.

    Returns
    -------
    int_g : 1D array of float (shape: 2)
        Integrals of (gx,gy) vector.

    int_g3 : float
        integrals of |g|^2 and |g|^3.
    """

class RadiationIntegrals:
    """radiation_integrals return type"""

    @property
    def mode(self) -> _pybmad.NormalModesStruct: ...

    @property
    def rad_int_by_ele(self) -> _pybmad.RadIntAllEleStruct: ...

    @property
    def ix_cache(self) -> int | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def radiation_integrals(lat: _pybmad.LatStruct, orbit: _pybmad.CoordStructArray1D, ix_cache: int | None = None, ix_branch: int | None = None) -> RadiationIntegrals:
    """
    Wrapper for Fortran routine radiation_integrals

    Parameters
    ----------
    lat : LatStruct
        Lattice to use. The calculation assumes that the Twiss parameters have been calculated.

    orbit : 1D array of CoordStruct
        Closed orbit for the branch.

    ix_cache : int, optional
        Cache pointer.
        This parameter is an input/output and is modified in-place.
        As an output, ix_cache: Cache pointer. If ix_cache = 0 at input then

    ix_branch : int, optional
        Lattice branch index. Default is 0.

    Returns
    -------
    mode : NormalModesStruct
        Parameters for the ("horizontal like") a-mode, ("vertical like") b-mode, and the z-mode

    ix_cache : int, optional
        Cache pointer.
        This parameter is an input/output and is modified in-place.
        As an output, ix_cache: Cache pointer. If ix_cache = 0 at input then

    rad_int_by_ele : RadIntAllEleStruct, optional
        Radiation integrals element by element.
    """

def radiation_map_setup(ele: _pybmad.EleStruct, ref_orbit_in: _pybmad.CoordStruct | None = None) -> bool:
    """
    Routine to calculate the radiation kick for a lattice element.

    Parameters
    ----------
    ele : EleStruct
        Element whose map is to be setup.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with map calculated.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def ramper_slave_setup(lat: _pybmad.LatStruct, force_setup: bool | None = None) -> None:
    """
    Wrapper for Fortran routine ramper_slave_setup

    Parameters
    ----------
    lat : LatStruct
        Lattice to be setup.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with ramper slaves setup.

    force_setup : bool, optional
        Default False. If True, do the setup even if lat.ramper_slave_bookkeeping = ok$. But the setup will never
        be done if lat.ramper_slave_bookkeeping = super_ok$.
    """

class RamperValue:
    """ramper_value return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ramper_value(ramper: _pybmad.EleStruct, r1: _pybmad.ControlRamp1Struct) -> RamperValue:
    """
    Wrapper for Fortran routine ramper_value

    Parameters
    ----------
    ramper : EleStruct
        Ramper lord.

    r1 : ControlRamp1Struct
        Slave function.

    Returns
    -------
    err_flag : bool
        Set True if there is an error, False otherwise.

    value : float
        Value of the slave function.
    """

def randomize_lr_wake_frequencies(ele: _pybmad.EleStruct) -> bool:
    """
    Routine to randomize the frequencies of the lr wake HOMs according to:
      freq = freq_in * (1 + lr_freq_spread) * rr)
    where rr is a Gaussian distributed random number with unit variance.

    Parameters
    ----------
    ele : EleStruct
        Element with wake. If no wake then nothing is done.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with wake frequencies set.

    Returns
    -------
    set_done : bool, optional
        Set True if there where lr wakes to be set. False otherwise.
    """

def rchomp(rel: float, plc: int) -> str:
    """
    Wrapper for Fortran routine rchomp

    Parameters
    ----------
    rel : float

    plc : int

    Returns
    -------
    out : str
    """

def re_allocate_eles(eles: _pybmad.ElePointerStructAlloc1D, n: int, save_old: bool | None = None, exact: bool | None = None) -> None:
    """
    Wrapper for Fortran routine re_allocate_eles

    Parameters
    ----------
    eles : 1D array of ElePointerStruct
        Array of element pointers with possible old data.
        This parameter is an input/output and is modified in-place.
        As an output, eles: Array of element pointers.

    n : int
        Array size to set.

    save_old : bool, optional
        If present and True then save the old data.

    exact : bool, optional
        If present and True then eles will have size = n If False (default), reallcation will not be done if eles
        is already large enough
    """

@overload
def re_allocate(section: _pybmad.Wall3DSectionStructAlloc1D, n: int, exact: bool | None = None) -> None:
    """
    Wrapper for Fortran routine re_allocate_wall3d_section_array

    Parameters
    ----------
    section : 1D array of Wall3dSectionStruct

    n : int

    exact : bool, optional
    """

@overload
def re_allocate(v: _pybmad.Wall3DVertexStructAlloc1D, n: int, exact: bool | None = None) -> None:
    """
    Wrapper for Fortran routine re_allocate_wall3d_vertex_array

    Parameters
    ----------
    v : 1D array of Wall3dVertexStruct

    n : int

    exact : bool, optional
    """

def re_associate_node_array(tree: _pybmad.ExpressionTreeStruct, n: int, exact: bool | None = None) -> None:
    """
    Routine to resize the tree%node(:) array.

    Note: The data of the array is preserved but data at the end of the
    array will be lost if n is less than the original size of the array

    Parameters
    ----------
    tree : ExpressionTreeStruct

    n : int
        Size wanted.

    exact : bool, optional
        Default is False. If False, the size of the output array is permitted to be larger than n.
    """

@overload
def re_str(rel: float) -> str:
    """
    Wrapper for Fortran routine re_str_qp

    Parameters
    ----------
    rel : float

    Returns
    -------
    str_out : str
    """

@overload
def re_str(rel: float) -> str:
    """
    Wrapper for Fortran routine re_str_rp

    Parameters
    ----------
    rel : float

    Returns
    -------
    str_out : str
    """

class ReadBeamAscii:
    """read_beam_ascii return type"""

    @property
    def beam(self) -> _pybmad.BeamStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def read_beam_ascii(file_name: str, beam_init: _pybmad.BeamInitStruct) -> ReadBeamAscii:
    """
    Subroutine to read in a beam definition file.
    If non_zero, the following components of beam_init are used to rescale the beam:
        %n_bunch
        %n_particle
        %charge_tot

    If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
    Otherwise the file is assumed to be ASCII.

    Parameters
    ----------
    file_name : str
        Name of beam file.

    beam_init : BeamInitStruct
        See above.

    Returns
    -------
    beam : BeamStruct
        Structure holding the beam information.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class ReadBeamFile:
    """read_beam_file return type"""

    @property
    def beam(self) -> _pybmad.BeamStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def read_beam_file(file_name: str, beam_init: _pybmad.BeamInitStruct, ele: _pybmad.EleStruct | None = None, print_mom_shift_warning: bool | None = None, conserve_momentum: bool | None = None) -> ReadBeamFile:
    """
    Subroutine to read in a beam definition file.
    If non_zero, the following components of beam_init are used to rescale the beam:
        %n_bunch
        %n_particle
        %bunch_charge -> charge_tot
        %species

    If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
    Otherwise the file is assumed to be ASCII.

    Parameters
    ----------
    file_name : str
        Name of beam file.

    beam_init : BeamInitStruct
        See above.

    ele : EleStruct, optional
        Element with reference energy, etc.

    print_mom_shift_warning : bool, optional
        Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

    Returns
    -------
    beam : BeamStruct
        Structure holding the beam information.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class ReadBinaryCartesianMap:
    """read_binary_cartesian_map return type"""

    @property
    def cart_map(self) -> _pybmad.CartesianMapStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def read_binary_cartesian_map(file_name: str, ele: _pybmad.EleStruct) -> ReadBinaryCartesianMap:
    """
    Routine to read a binary cartesian_map structure.

    Parameters
    ----------
    file_name : str
        File to create.

    ele : EleStruct
        Element associated with the map.

    Returns
    -------
    cart_map : CartesianMapStruct
        cartesian_map_struct, cartesian map.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class ReadBinaryCylindricalMap:
    """read_binary_cylindrical_map return type"""

    @property
    def cl_map(self) -> _pybmad.CylindricalMapStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def read_binary_cylindrical_map(file_name: str, ele: _pybmad.EleStruct) -> ReadBinaryCylindricalMap:
    """
    Routine to read a binary cylindrical_map structure.

    Parameters
    ----------
    file_name : str
        File to create.

    ele : EleStruct
        Element associated with the map.

    Returns
    -------
    cl_map : CylindricalMapStruct
        cylindrical_map_struct, cylindrical map.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class ReadBinaryGridField:
    """read_binary_grid_field return type"""

    @property
    def g_field(self) -> _pybmad.GridFieldStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def read_binary_grid_field(file_name: str, ele: _pybmad.EleStruct) -> ReadBinaryGridField:
    """
    Routine to read a binary grid_field structure.

    Parameters
    ----------
    file_name : str
        File to create.

    ele : EleStruct
        Element associated with the map.

    Returns
    -------
    g_field : GridFieldStruct
        grid_field_struct, cylindrical map.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class ReadDigestedBmadFile:
    """read_digested_bmad_file return type"""

    @property
    def lat(self) -> _pybmad.LatStruct: ...

    @property
    def inc_version(self) -> int: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def parser_calling(self) -> bool: ...

    @property
    def lat_files(self) -> _pybmad.CharacterAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def read_digested_bmad_file(digested_file: str) -> ReadDigestedBmadFile:
    """
    Wrapper for Fortran routine read_digested_bmad_file

    Parameters
    ----------
    digested_file : str
        Name of the digested file.

    Returns
    -------
    lat : LatStruct
        Output lattice structure

    inc_version : int
        bmad version number stored in the lattice file. If the file is current this number should be the same as
        the global parameter bmad_inc_version$. Set to -1 if there is a read error.

    err_flag : bool, optional
        Set True if there is an error. False otherwise.

    parser_calling : bool, optional
        Is this routine being called from a parser routine (like bmad_parser)? Default is False. This argument
        determines what are considered errors. For example, a moved digested file is considered an error if this
        routine is called from a parser but not otherwise. The reason for this dichotomy is that a parser is able
        to reread the original lattice file.

    lat_files : 1D array of str, optional
        List of Bmad lattice files that defined this lattice.
    """

def read_surface_reflection_file(file_name: str) -> _pybmad.PhotonReflectSurfaceStruct:
    """
    Routine to read the reflection probability data for a given type of surface from a file.

    Parameters
    ----------
    file_name : str
        Name of the file.

    Returns
    -------
    surface : PhotonReflectSurfaceStruct
        Surface info.
    """

def reallocate_beam(beam: _pybmad.BeamStruct, n_bunch: int, n_particle: int | None = None, extend: bool | None = None) -> None:
    """
    Wrapper for Fortran routine reallocate_beam

    Parameters
    ----------
    beam : BeamStruct
        Beam bunches are saved if save = True.
        This parameter is an input/output and is modified in-place.
        As an output, beam: Allocated beam_struct structure.

    n_bunch : int
        Number of bunches.

    n_particle : int, optional
        Number of particles. Must be non-negative. If save = True then the number of particles in existing bunches
        will not be touched. If not present, beam.bunch(i).particle(:) will be in an undefined state.

    extend : bool, optional
    """

def reallocate_bp_com_const() -> None:
    """Wrapper for Fortran routine reallocate_bp_com_const"""

def reallocate_bunch(n_particle: int, save: bool | None = None) -> _pybmad.BunchStruct:
    """
    Wrapper for Fortran routine reallocate_bunch

    Parameters
    ----------
    n_particle : int
        Number of particles. Must be non-negative.

    save : bool, optional
        If present and True then save the old bunch info.

    Returns
    -------
    bunch : BunchStruct
        Allocated bunch_struct structure.
    """

def reallocate_control(lat: _pybmad.LatStruct, n: int) -> None:
    """
    Wrapper for Fortran routine reallocate_control

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    n : int
        Array size for lat.control(:) and lat.ic(:).
    """

@overload
def reallocate_coord(coord_array: _pybmad.CoordArrayStructAlloc1D, lat: _pybmad.LatStruct) -> None:
    """
    Routine to allocate or reallocate at allocatable coord_struct array.
    reallocate_coord is an overloaded name for:
      reallocate_coord_n (coord, n_coord)
      reallocate_coord_lat (coord, lat, ix_branch)

    Subroutine to allocate an allocatable coord_struct array to at least:
        coord(0:n_coord)                            if n_coord arg is used.
        coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

    The old coordinates are saved
    If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
    In any case, coord(n)%vec for n > 0 is set to zero.

    Parameters
    ----------
    lat : LatStruct
        Lattice
    """

@overload
def reallocate_coord(coord: _pybmad.CoordStructAlloc1D, lat: _pybmad.LatStruct, ix_branch: int | None = None) -> None:
    """
    Routine to allocate or reallocate at allocatable coord_struct array.
    reallocate_coord is an overloaded name for:
      reallocate_coord_n (coord, n_coord)
      reallocate_coord_lat (coord, lat, ix_branch)

    Subroutine to allocate an allocatable coord_struct array to at least:
        coord(0:n_coord)                            if n_coord arg is used.
        coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

    The old coordinates are saved
    If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
    In any case, coord(n)%vec for n > 0 is set to zero.

    Parameters
    ----------
    coord : 1D array of CoordStruct
        Allocatable array.
        This parameter is an input/output and is modified in-place.
        As an output, coord: Allocated array.

    lat : LatStruct
        Lattice

    ix_branch : int, optional
        Branch to use. Default is 0 (main branch).
    """

@overload
def reallocate_coord(coord: _pybmad.CoordStructAlloc1D, n_coord: int) -> None:
    """
    Routine to allocate or reallocate at allocatable coord_struct array.
    reallocate_coord is an overloaded name for:
      reallocate_coord_n (coord, n_coord)
      reallocate_coord_lat (coord, lat, ix_branch)

    Subroutine to allocate an allocatable coord_struct array to at least:
        coord(0:n_coord)                            if n_coord arg is used.
        coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

    The old coordinates are saved
    If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
    In any case, coord(n)%vec for n > 0 is set to zero.

    Parameters
    ----------
    coord : 1D array of CoordStruct
        Allocatable array.
        This parameter is an input/output and is modified in-place.
        As an output, coord: Allocated array.

    n_coord : int
        Minimum array upper bound wanted.
    """

def reallocate_expression_stack(stack: _pybmad.ExpressionAtomStructAlloc1D, n: int, exact: bool | None = None) -> None:
    """
    Wrapper for Fortran routine reallocate_expression_stack

    Parameters
    ----------
    stack : 1D array of ExpressionAtomStruct
        Existing stack array.
        This parameter is an input/output and is modified in-place.
        As an output, stack: Resized stack.

    n : int
        Array size needed.

    exact : bool, optional
        If present and False then the size of the output array is permitted to be larger than n. Default is True.
    """

def reallocate_sequence(sequence: _pybmad.SeqStructAlloc1D, n_seq: int) -> None:
    """No docstring available."""

def rel_tracking_charge_to_mass(orbit: _pybmad.CoordStruct, ref_species: int) -> float:
    """
    Wrapper for Fortran routine rel_tracking_charge_to_mass

    Parameters
    ----------
    orbit : CoordStruct
        Particle position structure.

    ref_species : int
        Reference species

    Returns
    -------
    rel_charge : float
        Relative charge/mass
    """

def relative_mode_flip(ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine relative_mode_flip

    Parameters
    ----------
    ele1 : EleStruct
        Elements to compare.

    ele2 : EleStruct
        Elements to compare.
    """

class ReleaseRadIntCache:
    """release_rad_int_cache return type"""

    @property
    def ix_cache(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def release_rad_int_cache(ix_cache: int) -> ReleaseRadIntCache:
    """
    Subroutine to release the memory associated with caching wiggler values.
    See the radiation_integrals routine for further details.

    Parameters
    ----------
    ix_cache : int
        Cache number.
        This parameter is an input/output and is modified in-place.
        As an output, ix_cache: Cache number set to 0,

    Returns
    -------
    ix_cache : int
        Cache number.
        This parameter is an input/output and is modified in-place.
        As an output, ix_cache: Cache number set to 0,
    """

def remove_constant_taylor(taylor_in: _pybmad.TaylorStructArray1D, taylor_out: _pybmad.TaylorStructArray1D, c0: _pybmad.RealArray1D, remove_higher_order_terms: bool) -> None:
    """
    Subroutine to remove the constant part of a taylor map.
    Optionally terms that are higher order than bmad_com%taylor_order can
    be removed.

    Note: It is assumed that taylor_out has been deallocated before the call to
    this routine. Calling this routine with the first two actual arguments the
    same is prohibited.

    Parameters
    ----------
    taylor_in : 1D array of TaylorStruct
        Input taylor map.

    taylor_out : 1D array of TaylorStruct
        Taylor with constant terms removed.

    c0 : 1D array of float
        The constant part of the taylor map

    remove_higher_order_terms : bool
        If True then terms that are higher order than bmad_com.taylor_order are removed.
    """

def remove_dead_from_bunch(bunch_in: _pybmad.BunchStruct) -> _pybmad.BunchStruct:
    """
    Wrapper for Fortran routine remove_dead_from_bunch

    Parameters
    ----------
    bunch_in : BunchStruct
        Input bunch with alive and dead particles.

    Returns
    -------
    bunch_out : BunchStruct
        Output bunch with only alive and pre_born particles. Note: bunch_out can be the same actual argument as
        bunch_in.
    """

def remove_eles_from_lat(lat: _pybmad.LatStruct, check_sanity: bool | None = None) -> None:
    """
    Wrapper for Fortran routine remove_eles_from_lat

    Parameters
    ----------
    lat : LatStruct
        Lattice to compress.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Compressed lattice.

    check_sanity : bool, optional
        If True (default) then call lat_sanity_check
    """

def remove_lord_slave_link(lord: _pybmad.EleStruct, slave: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine remove_lord_slave_link

    Parameters
    ----------
    lord : EleStruct
        Lord element
        This parameter is an input/output and is modified in-place.
        As an output, lord: Lord element with link info removed

    slave : EleStruct
        Slave element
        This parameter is an input/output and is modified in-place.
        As an output, slave: Slave element with link info removed
    """

def residual_pwd_sig_z(zz: float, status: int | None = None) -> float:
    """
    Wrapper for Fortran routine residual_pwd_sig_z

    Parameters
    ----------
    zz : float

    status : int, optional
    """

def reverse_lat(lat_in: _pybmad.LatStruct, track_antiparticle: bool | None = None) -> _pybmad.LatStruct:
    """
    Wrapper for Fortran routine reverse_lat

    Parameters
    ----------
    lat_in : LatStruct
        Input lattice to reverse.

    track_antiparticle : bool, optional
        Set the particle species of the reversed lat to the anti-particle of lat_in? Default is True.

    Returns
    -------
    lat_rev : LatStruct
        Reversed lattice.
    """

def rf_coupler_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, phase: float, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    No longer in the codebase
    function rf_clock_setup (branch, n_rf_included, n_rf_excluded) result (ok)
      import
      implicit none
      type (branch_struct), target :: branch
      integer n_rf_included, n_rf_excluded
      logical ok
    end function

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through

    param : LatParamStruct
        branch parameters.

    particle_at : int
        first_track_edge$, or second_track_edge$.

    phase : float
        phase of cavity

    orbit : CoordStruct
        Position before kick.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Position after kick.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def rf_is_on(branch: _pybmad.BranchStruct, ix_ele1: int | None = None, ix_ele2: int | None = None) -> bool:
    """
    Wrapper for Fortran routine rf_is_on

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch to check.

    ix_ele1 : int, optional
        Start of range of elements to check. Default is 0.

    ix_ele2 : int, optional
        End of range of elements to check. Default is branch.n_ele_track.

    Returns
    -------
    is_on : bool
        True if any rfcavity is powered. False otherwise.
    """

def rf_ref_time_offset(ele: _pybmad.EleStruct, ds: float | None = None) -> float:
    """
    Wrapper for Fortran routine rf_ref_time_offset

    Parameters
    ----------
    ele : EleStruct
        RF Element being tracked through.

    ds : float, optional
        Distance of particle from start edge. Default is zero.

    Returns
    -------
    time : float
        Offset time.
    """

def rfun(u: float, v: float, w: float, gam: float, a: float, b: float, hz: float, i: int, j: int) -> float:
    """
    Wrapper for Fortran routine rfun

    Parameters
    ----------
    u : float

    v : float

    w : float

    gam : float

    a : float

    b : float

    hz : float

    i : int

    j : int

    Returns
    -------
    res : float
    """

def rk_adaptive_time_step(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orb: _pybmad.CoordStruct, t_dir: int, rf_time: float, dt_try: float, dt_did: float, dt_next: float, err_flag: bool, extra_field: _pybmad.EmFieldStruct | None = None) -> None:
    """
    Wrapper for Fortran routine rk_adaptive_time_step

    Parameters
    ----------
    ele : EleStruct

    param : LatParamStruct

    orb : CoordStruct

    t_dir : int

    rf_time : float

    dt_try : float

    dt_did : float

    dt_next : float

    err_flag : bool

    extra_field : EmFieldStruct, optional
    """

def rk_time_step1(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, rf_time: float, orb: _pybmad.CoordStruct, dt: float, new_orb: _pybmad.CoordStruct, err_flag: bool, dr_dt: Sequence[float] | None = None, print_err: bool | None = None, extra_field: _pybmad.EmFieldStruct | None = None) -> list[float]:
    """
    Wrapper for Fortran routine rk_time_step1

    Parameters
    ----------
    ele : EleStruct

    param : LatParamStruct

    rf_time : float

    orb : CoordStruct

    dt : float

    new_orb : CoordStruct

    err_flag : bool

    dr_dt : 1D array of float (shape: 10), optional

    print_err : bool, optional

    extra_field : EmFieldStruct, optional

    Returns
    -------
    r_err : 1D array of float (shape: 10)
    """

def rotate3(vec: Sequence[float], angle: float) -> list[float]:
    """
    Wrapper for Fortran routine rotate3

    Parameters
    ----------
    vec : 1D array of float (shape: 3)

    angle : float

    Returns
    -------
    rvec : 1D array of float (shape: 3)
    """

def rotate_em_field(field: _pybmad.EmFieldStruct, w_mat: Sequence[Sequence[float]], w_inv: Sequence[Sequence[float]], calc_dfield: bool | None = None, calc_potential: bool | None = None) -> None:
    """
    Routine to transform the fields using the given rotation matrices.

    Parameters
    ----------
    field : EmFieldStruct
        E and B fields and derivatives.

    w_mat : 2D array of float (shape: 3,3)
        rotation matrix.

    w_inv : 2D array of float (shape: 3,3)
        rotation matrix inverse = transpose(w_mat)

    calc_dfield : bool, optional
        If present and True then rotate the field derivatives.

    calc_potential : bool, optional
        Rotate the magnetic vector potential? Default is false.
    """

def rotate_field_zx(field: _pybmad.EmFieldStruct, theta: float) -> None:
    """
    Wrapper for Fortran routine rotate_field_zx

    Parameters
    ----------
    field : EmFieldStruct

    theta : float
    """

def rotate_for_curved_surface(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, set: bool, rot_mat: Sequence[Sequence[float]]) -> None:
    """
    Wrapper for Fortran routine rotate_for_curved_surface

    Parameters
    ----------
    ele : EleStruct
        reflecting element

    orbit : CoordStruct
        Photon position.

    set : bool
        True -> Transform body coords to local curved body coords. False -> Transform local curved body to body
        coords.

    rot_mat : 2D array of float (shape: 3,3)
        When set = False, rotation matrix calculated from previous call with set = True.
        This parameter is an input/output and is modified in-place.
        As an output, rot_mat: When set = True, calculated rotation matrix.
    """

def rotate_spin(rot_vec: Sequence[float], spin: Sequence[float]) -> list[float]:
    """
    Wrapper for Fortran routine rotate_spin

    Parameters
    ----------
    rot_vec : 1D array of float (shape: 3)
        Rotation axis. Magnitude of rot_vec is the rotation angle.

    spin : 1D array of float (shape: 3)
        Initial coords.
        This parameter is an input/output and is modified in-place.
        As an output, spin: Final coords.

    Returns
    -------
    qrot : 1D array of float (shape: 0:3), optional
        : rotation quaternion.
    """

def rotate_spin_a_step(orbit: _pybmad.CoordStruct, field: _pybmad.EmFieldStruct, ele: _pybmad.EleStruct, ds: float) -> None:
    """
    Wrapper for Fortran routine rotate_spin_a_step

    Parameters
    ----------
    orbit : CoordStruct
        Initial orbit.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit with rotated spin

    field : EmFieldStruct
        EM Field

    ele : EleStruct
        ele_struct, Element being tracked through.

    ds : float
        Longitudinal step in element body frame.
    """

def rotate_spin_given_field(orbit: _pybmad.CoordStruct, sign_z_vel: int, BL: Sequence[float] | None = None, EL: Sequence[float] | None = None, qrot: Sequence[float] | None = None) -> None:
    """
    Wrapper for Fortran routine rotate_spin_given_field

    Parameters
    ----------
    orbit : CoordStruct
        Initial orbit.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit with rotated spin

    sign_z_vel : int
        +/- 1. Sign of direction of travel relative to the element.

    BL : 1D array of float (shape: 3), optional
        Integrated field strength. Assumed zero if not present.

    EL : 1D array of float (shape: 3), optional
        Integrated field strength. Assumed zero if not present.

    qrot : 1D array of float (shape: 0:3), optional
        Initial rotation quaternion.
        This parameter is an input/output and is modified in-place.
        As an output, qrot: Rotation quaternion with rotation due to the field added in.
    """

def s_body_calc(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> float:
    """
    Wrapper for Fortran routine s_body_calc

    Parameters
    ----------
    orbit : CoordStruct
        Particle coordinates.

    ele : EleStruct
        Lattice element

    Returns
    -------
    s_body : float
        Body postion.
    """

def s_calc(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine s_calc

    Parameters
    ----------
    lat : LatStruct
    """

def s_ref_to_s_chord(s_ref: float, eleinfo: _pybmad.CsrEleInfoStruct) -> float:
    """
    Routine to calculate s_chord given s_ref.

    Parameters
    ----------
    s_ref : float
        s-position along element ref coords.

    eleinfo : CsrEleInfoStruct
        Element info

    Returns
    -------
    s_chord : float
        s-position along centroid chord.
    """

class SSourceCalc:
    """s_source_calc return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def s_source(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def s_source_calc(kick1: _pybmad.CsrKick1Struct, csr: _pybmad.CsrStruct, dr_match: Sequence[float]) -> SSourceCalc:
    """
    Routine to calculate the distance between source and kick points.

    Parameters
    ----------
    kick1 : CsrKick1Struct

    csr : CsrStruct

    dr_match : 1D array of float (shape: 3)
        Discontinuity factor if there is a match element between source and kick elements.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. Untouched otherwise.

    s_source : float
        source s-position.
    """

def sad_mult_hard_bend_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Routine to track through the hard edge bend fringe field for a bend or sad_mult element.
    Only the bend field is taken into account here. Higher order multipolse must be handled elsewhere.

    This routine assumes that the particle coordinates are with respect to the actual magnet face.
    Thus finite e1/e2 must be taken into account by other routines.

    SAD calls this the "linear" fringe even though it is nonlinear.

    Parameters
    ----------
    ele : EleStruct
        Element with fringe.

    param : LatParamStruct
        Tracking parameters.

    particle_at : int
        Either first_track_edge$ or second_track_edge$.

    orbit : CoordStruct
        Starting coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coordinates.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the fringe.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix including the fringe.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

def sad_soft_bend_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orb: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Subroutine to track through the ("linear") bend soft edge field of an sbend or sad_mult.

    Parameters
    ----------
    ele : EleStruct
        SBend or sad_mult element.

    param : LatParamStruct

    particle_at : int
        first_track_edge$, or second_track_edge$.

    orb : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Coords after tracking.

    mat6 : 2D array of float (shape: 6,6), optional
        Starting matrix
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix after fringe field

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

def save_a_beam_step(ele: _pybmad.EleStruct, beam: _pybmad.BeamStruct, bunch_tracks: _pybmad.BunchTrackStructArray1D | None = None, s_body: float | None = None, is_time_coords: bool | None = None) -> None:
    """
    Wrapper for Fortran routine save_a_beam_step

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    beam : BeamStruct
        Bunches in the beam whose parameters are to be saved.

    bunch_tracks : 1D array of BunchTrackStruct, optional
        Track up to now. If bunch_tracks.n_pt < 0, the structure will be reinitialized.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_tracks: Track with current bunch info appended on. This routine does nothing

    s_body : float, optional
        Body s-position from beginning of element.

    is_time_coords : bool, optional
        Default is False. If True, input beam is using time coordinates in which case there will be a conversion
        to s-coords before bunch_params are computed.
    """

def save_a_bunch_step(ele: _pybmad.EleStruct, bunch: _pybmad.BunchStruct, bunch_track: _pybmad.BunchTrackStruct | None = None, s_body: float | None = None, is_time_coords: bool | None = None) -> None:
    """
    Wrapper for Fortran routine save_a_bunch_step

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    bunch : BunchStruct
        Bunch whose parameters are to be saved.

    bunch_track : BunchTrackStruct, optional
        Track up to now. If bunch_track.n_pt < 0, the structure will be reinitialized.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: Track with current bunch info appended on. This routine does nothing

    s_body : float, optional
        Body s-position from beginning of element.

    is_time_coords : bool, optional
        Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
        to s-coords before bunch_params are computed.
    """

def save_a_step(track: _pybmad.TrackStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, local_ref_frame: bool, orb: _pybmad.CoordStruct, s_rel: float, save_field: bool | None = None, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None, rf_time: float | None = None, strong_beam: _pybmad.StrongBeamStruct | None = None) -> None:
    """
    Wrapper for Fortran routine save_a_step

    Parameters
    ----------
    track : TrackStruct
        Track up to now. If track.n_pt < 0, the structure will be reinitialized.
        This parameter is an input/output and is modified in-place.
        As an output, track: Track with current position appended on.

    ele : EleStruct
        Element being tracked through.

    param : LatParamStruct
        Lattice parameters.

    local_ref_frame : bool
        If True then input orb is with respect to body coordinates.

    orb : CoordStruct
        trajectory at s with respect to element coordinates.

    s_rel : float
        Longitudinal position wrt the element. If local_ref_frame = F: Lab coords. If local_ref_frame = T: body
        coords.

    save_field : bool, optional
        Save electric and magnetic field values? Default is False.

    mat6 : 2D array of float (shape: 6,6), optional
        Matrix to store.

    make_matrix : bool, optional
        Is mat6 a valid matrix? Default is False.

    rf_time : float, optional
        RF clock time used for calculating the field.. If not present then the time will be calculated using the
        standard algorithm. This is only needed if save_field = True.

    strong_beam : StrongBeamStruct, optional
        Strong beam info if tracking through a beambeam element.
    """

def sbend_body_with_k1_map(ele: _pybmad.EleStruct, dg: float, b1: float, param: _pybmad.LatParamStruct, n_step: int, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine sbend_body_with_k1_map

    Parameters
    ----------
    ele : EleStruct
        Sbend element.

    dg : float
        Field error.

    b1 : float
        b1 quadrupole strength * rel_charge_dir

    param : LatParamStruct
        Branch parameters.

    n_step : int
        Number of steps to divide the bend into. Only one step is taken by this routine.

    orbit : CoordStruct
        Orbit at beginning of the bend.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coordinates.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix with body added in.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

class ScAdaptiveStep:
    """sc_adaptive_step return type"""

    @property
    def dt_next(self) -> float: ...

    @property
    def include_image(self) -> bool: ...

    @property
    def dt_step(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def sc_adaptive_step(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, include_image: bool, t_now: float, dt_step: float, sc_field: _pybmad.EmFieldStructArray1D) -> ScAdaptiveStep:
    """
    Routine to track a bunch of particles with space charge for one step using
    adaptive step size control and determine appropriate step size for the next step

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position in t-based coordinates
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position in t-based coordinates.

    ele : EleStruct
        Nominal lattice element being tracked through.

    include_image : bool
        Include image charge forces?
        This parameter is an input/output and is modified in-place.
        As an output, include_image: Set False if image charge calc no longer needed (Note

    t_now : float
        Current time at the beginning of tracking

    dt_step : float
        Initial SC time step to take
        This parameter is an input/output and is modified in-place.
        As an output, dt_step: Step done.

    sc_field : 1D array of EmFieldStruct
        : Array to hold space charge fields. Its length should be the number of particles.

    Returns
    -------
    include_image : bool
        Include image charge forces?
        This parameter is an input/output and is modified in-place.
        As an output, include_image: Set False if image charge calc no longer needed (Note

    dt_step : float
        Initial SC time step to take
        This parameter is an input/output and is modified in-place.
        As an output, dt_step: Step done.

    dt_next : float
        Next SC time step the tracker would take based on the error tolerance
    """

class ScStep:
    """sc_step return type"""

    @property
    def n_emit(self) -> int: ...

    @property
    def include_image(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def sc_step(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, include_image: bool, t_end: float, sc_field: _pybmad.EmFieldStructArray1D) -> ScStep:
    """
    Subroutine to track a bunch through a given time step with space charge

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position in t-based coordinates
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position in t-based coordinates after space charge kick.

    ele : EleStruct
        Nominal element being tracked through.

    include_image : bool
        Include image charge forces?
        This parameter is an input/output and is modified in-place.
        As an output, include_image: Set False if image charge calc no longer needed (Note

    t_end : float
        Time at which the tracking ends.

    sc_field : 1D array of EmFieldStruct
        : Array to hold space charge fields. Its length should be the number of particles.

    Returns
    -------
    include_image : bool
        Include image charge forces?
        This parameter is an input/output and is modified in-place.
        As an output, include_image: Set False if image charge calc no longer needed (Note

    n_emit : int, optional
        The number of particles emitted in this step.
    """

def set_active_fixer(fixer: _pybmad.EleStruct, turn_on: bool | None = None) -> _pybmad.CoordStruct:
    """
    Set the acvitive fixer element.
    All other fixer/beginning_ele elements in the branch will be deactivated.

    If turn_on is True (default), the fixer argument becomes the active fixer.
    If turn_on is False, and fixer%is_on is also False, there is nothing to be done.
    If turn_on is False, and fixer%is_on is True, turn this fixer off and turn on the beginning element.

    Parameters
    ----------
    fixer : EleStruct
        Fixer element to make active.
        This parameter is an input/output and is modified in-place.
        As an output, fixer: Element is now active.

    turn_on : bool, optional
        If True (default), make this fixer the active element. If False, make the beginning element active.

    Returns
    -------
    orbit : CoordStruct, optional
        Load with stored fixer phase space and spin values.
    """

class SetBranchAndEleForOmp:
    """set_branch_and_ele_for_omp return type"""

    @property
    def branch(self) -> _pybmad.BranchStruct | None: ...

    @property
    def ele0(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def set_branch_and_ele_for_omp(ix_lat: int, ele0_loc: _pybmad.LatEleLocStruct, lats: _pybmad.LatPointerStructArray1D, lat: _pybmad.LatStruct) -> SetBranchAndEleForOmp:
    """
    Routine to set branch and ele for the lattice to be used to track through.

    Different lattices are used to avoid problems when there is ramping since ramping
    will modify lattice element parameters.

    Parameters
    ----------
    ix_lat : int
        Lattice index.

    ele0_loc : LatEleLocStruct
        ele0 location in lattice.

    lats : 1D array of LatPointerStruct
        Pointers to lattices to track through.
        This parameter is an input/output and is modified in-place.
        As an output, lats: Properly initialized if needed.

    lat : LatStruct
        Original lattice.

    Returns
    -------
    branch : BranchStruct, optional
        Branch to track through

    ele0 : EleStruct, optional
        starting element for tracking.
    """

def set_custom_attribute_name(custom_name: str, custom_index: int | None = None) -> bool:
    """
    Routine to add custom element attributes to the element attribute name table.

    Parameters
    ----------
    custom_name : str
        Name of the custom attribute. If prefixed by "<class>::" then the custom name will be set only for that
        element class. Example: "quadrupole::error" will set the alias custom namefor quadrupoles.

    custom_index : int, optional
        Index used in assigning where in the ele_struct the custom attribute is put. If not present or 0 then the
        next unused slot is used.

    Returns
    -------
    err_flag : bool
        Set True if an error. False otherwise.
    """

class SetEleAttribute:
    """set_ele_attribute return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def err_id(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def set_ele_attribute(ele: _pybmad.EleStruct, set_string: str, err_print_flag: bool | None = None, set_lords: bool | None = None) -> SetEleAttribute:
    """
    Wrapper for Fortran routine set_ele_attribute

    Parameters
    ----------
    ele : EleStruct
        Element with attribute to set.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with attribute set.

    set_string : str
        Attribute and value for set.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message if attribute is, for example, not free.

    set_lords : bool, optional
        Default False. If True, set the super_lord(s) if the element is a super_slave.

    Returns
    -------
    err_flag : bool
        Set True if there is an error, False otherwise.

    err_id : int, optional
        Set to an integer which identifies the error type. 0 = no error. The higher the error the further along
        the error was encountered.
    """

def set_ele_defaults(ele: _pybmad.EleStruct, do_allocate: bool | None = None) -> None:
    """
    Wrapper for Fortran routine set_ele_defaults

    Parameters
    ----------
    ele : EleStruct
        Element to init.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Initialized element.

    do_allocate : bool, optional
        Do default allocation of element components? Default is True.
    """

def set_ele_name(ele: _pybmad.EleStruct, name: str) -> None:
    """
    Wrapper for Fortran routine set_ele_name

    Parameters
    ----------
    ele : EleStruct
        Element whose name is to be set.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with name set.

    name : str
        Name to set.
    """

def set_ele_real_attribute(ele: _pybmad.EleStruct, attrib_name: str, value: float, err_print_flag: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine set_ele_real_attribute

    Parameters
    ----------
    ele : EleStruct
        Element with attribute to set.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with attribute set.

    attrib_name : str
        Attribute name.

    value : float
        value to set to.

    err_print_flag : bool, optional
        If present and False then suppress printing of an error message if attribute is, for example, not free.

    Returns
    -------
    err_flag : bool
        Set True if there is an error, False otherwise.
    """

def set_ele_status_stale(ele: _pybmad.EleStruct, status_group: int, set_slaves: bool | None = None, old_eles: _pybmad.ElePointerStructAlloc1D | None = None) -> None:
    """
    Wrapper for Fortran routine set_ele_status_stale

    Parameters
    ----------
    ele : EleStruct
        Element to set.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element.

    status_group : int
        Which flag groups to set. Possibilities are: attribute_group$, control_group$, floor_position_group$,
        s_position_group$, s_and_floor_position_group$, ref_energy_group$, or mat6_group$, all_groups$

    set_slaves : bool, optional
        If present and False then do not set the status for any slaves. Default is True.

    old_eles : 1D array of ElePointerStruct, optional
        List of elements already set. This argument is only used when this routine is called recursively.
    """

@overload
def set_flags_for_changed_attribute(ele: _pybmad.EleStruct, all_attrib: _pybmad.AllPointerStruct, set_dependent: bool | None = None) -> None:
    """
    Routine to mark an element or lattice as modified.
    Also will do some dependent variable bookkeeping when a particular attribute has
    been altered.

    This routine should be called after the attribute has been set.

    set_flags_for_changed_attribute is an overloaded name for:
      set_flags_for_changed_lat_attribute (lat, set_dependent)
      set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
      set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
      set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
      set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

    The set_flags_for_changed_lat_attribute (lat) routine is used when one
    does not know what has changed and wants a complete bookkeeping done.

    NOTE: The attribute argument MUST be the component that was changed. For example:
        ele%value(x_offset$) = off_value
        call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
    And NOT:
        call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

    Parameters
    ----------
    ele : EleStruct
        ele_struct, Element being modified.

    all_attrib : AllPointerStruct
        Pointer to attribute.

    set_dependent : bool, optional
        If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
        when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
        doing.
    """

@overload
def set_flags_for_changed_attribute(ele: _pybmad.EleStruct, attrib: int, set_dependent: bool | None = None) -> None: ...

@overload
def set_flags_for_changed_attribute(lat: _pybmad.LatStruct, set_dependent: bool | None = None) -> None:
    """
    Routine to mark an element or lattice as modified.
    Also will do some dependent variable bookkeeping when a particular attribute has
    been altered.

    This routine should be called after the attribute has been set.

    set_flags_for_changed_attribute is an overloaded name for:
      set_flags_for_changed_lat_attribute (lat, set_dependent)
      set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
      set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
      set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
      set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

    The set_flags_for_changed_lat_attribute (lat) routine is used when one
    does not know what has changed and wants a complete bookkeeping done.

    NOTE: The attribute argument MUST be the component that was changed. For example:
        ele%value(x_offset$) = off_value
        call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
    And NOT:
        call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

    Parameters
    ----------
    lat : LatStruct
        Lattice being modified.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with appropriate changes.

    set_dependent : bool, optional
        If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
        when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
        doing.
    """

@overload
def set_flags_for_changed_attribute(ele: _pybmad.EleStruct, attrib: bool, set_dependent: bool | None = None) -> None: ...

@overload
def set_flags_for_changed_attribute(ele: _pybmad.EleStruct, attrib: float | None = None, set_dependent: bool | None = None) -> None:
    """
    Routine to mark an element or lattice as modified.
    Also will do some dependent variable bookkeeping when a particular attribute has
    been altered.

    This routine should be called after the attribute has been set.

    set_flags_for_changed_attribute is an overloaded name for:
      set_flags_for_changed_lat_attribute (lat, set_dependent)
      set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
      set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
      set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
      set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

    The set_flags_for_changed_lat_attribute (lat) routine is used when one
    does not know what has changed and wants a complete bookkeeping done.

    NOTE: The attribute argument MUST be the component that was changed. For example:
        ele%value(x_offset$) = off_value
        call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
    And NOT:
        call set_flags_for_changed_attribute (ele, off_value)  ! WRONG

    Parameters
    ----------
    ele : EleStruct
        ele_struct, Element being modified.

    set_dependent : bool, optional
        If False then dependent parameter bookkeeping will not be done. False is used, for example, during parsing
        when dependent bookkeepin is not wanted. Default is True. Do not set False unless you know what you are
        doing.
    """

class SetFringeOnOff:
    """set_fringe_on_off return type"""

    @property
    def fringe_at(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def set_fringe_on_off(fringe_at: float, ele_end: int, on_or_off: int) -> SetFringeOnOff:
    """
    Wrapper for Fortran routine set_fringe_on_off

    Parameters
    ----------
    fringe_at : float
        Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
        This parameter is an input/output and is modified in-place.
        As an output, fringe_at: Modified fringe setting.

    ele_end : int
        Element edge: entrance_end$ or exit_end$

    on_or_off : int
        Turn on$ or off$

    Returns
    -------
    fringe_at : float
        Present fringe_at setting. entrance_end$, exit_end$, both_ends$, or no_end$
        This parameter is an input/output and is modified in-place.
        As an output, fringe_at: Modified fringe setting.
    """

def set_lords_status_stale(ele: _pybmad.EleStruct, stat_group: int, control_bookkeeping: bool | None = None, flag: int | None = None) -> None:
    """
    Wrapper for Fortran routine set_lords_status_stale

    Parameters
    ----------
    ele : EleStruct
        Element

    stat_group : int
        which status group to set. floor_position_group$, etc. See set_ele_status_stale for more details.

    control_bookkeeping : bool, optional
        Call control_bookkeeper for each lord if needed? Default if False.

    flag : int, optional
        Do not use. For coordinating recursion.
    """

def set_on_off(key: int, lat: _pybmad.LatStruct, switch_: int, orb: _pybmad.CoordStructArray1D | None = None, use_ref_orb: bool | None = None, ix_branch: int | None = None, saved_values: _pybmad.RealAlloc1D | None = None, attribute: str | None = None, set_val: int | None = None) -> None:
    """
    Wrapper for Fortran routine set_on_off

    Parameters
    ----------
    key : int
        Class name of elements to be turned on or off. [quadrupole$, etc.]

    lat : LatStruct
        lattice structure holding the elements.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Modified lattice.

    orb : 1D array of CoordStruct, optional
        Needed for lat_make_mat6

    use_ref_orb : bool, optional
        If present and true then use ele.map_ref_orb for the reference orbit for calculating .mat6. Default is
        false.

    ix_branch : int, optional
        If present then only set for this lattice branch.

    saved_values : 1D array of float, optional
        Element-by element saved values of the component. Must be present if needed (EG if switch =
        restore_state$, etc.).
        This parameter is an input/output and is modified in-place.
        As an output, saved_values: Saved values of the component.

    attribute : str, optional
        Attribute to turn on/off. Eg: 'K2', 'MULTIPOLE_ON', etc. Default is 'IS_ON'. Must be upper case.

    set_val : int, optional
        Value to set to. Overrides normal set value.
    """

def set_orbit_to_zero(orbit: _pybmad.CoordStructArray1D, n1: int, n2: int, ix_noset: int | None = None) -> None:
    """
    Wrapper for Fortran routine set_orbit_to_zero

    Parameters
    ----------
    orbit : 1D array of CoordStruct
        Array with particle positions in the range orbit(n1:n2) set to zero except for orbit(ix_noset).

    n1 : int
        Lower bound of orbit(:) array subset.

    n2 : int
        Upper bound of orbit(:) array subset.

    ix_noset : int, optional
        If present then orbit(ix_noset) will not be zeroed.
    """

def set_ptc(e_tot: float | None = None, particle: int | None = None, taylor_order: int | None = None, integ_order: int | None = None, n_step: int | None = None, no_cavity: bool | None = None, force_init: bool | None = None) -> None:
    """
    Wrapper for Fortran routine set_ptc

    Parameters
    ----------
    e_tot : float, optional
        Energy in eV.

    particle : int, optional
        Type of particle: electron$, proton$, etc.

    taylor_order : int, optional
        Maximum order of the taylor polynomials. 0 => Use default.

    integ_order : int, optional
        Default Order for the drift-kick-drift sympletic integrator. Possibilities are: 2, 4, or 6 Default = 2

    n_step : int, optional
        Default Number of integration steps. Default = 1

    no_cavity : bool, optional
        No RF Cavity exists? Default = False. Corresponds to the nocavity option of the PTC init routine.
        no_cavity = .true. will turn any cavity into a drift.

    force_init : bool, optional
        If present and True then force a PTC init.
    """

def set_ptc_base_state(component: str, set_val: bool) -> bool:
    """
    Wrapper for Fortran routine set_ptc_base_state

    Parameters
    ----------
    component : str
        Name of component. "TOTALPATH", "SPIN", "NOCAVITY", "TIME", etc. See the PTC internal_state structure for
        component names.

    set_val : bool
        Value to set to. For TOTALPATH, True => 1, False => 0.

    Returns
    -------
    old_val : bool, optional
        Old value.
    """

def set_ptc_com_pointers() -> None:
    """Routine to set ptc_com pointers to PTC global variables."""

class SetPtcQuiet:
    """set_ptc_quiet return type"""

    @property
    def old_val(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def set_ptc_quiet(channel: int, set: bool, old_val: int) -> SetPtcQuiet:
    """
    Routine to set the lielib_print(:) array or c_verbose logical to suppress informational messages
    that can clutter the output from a program using PTC.

    Note: Only suppress printing if ptc_com%print_info_messages = F.

    Parameters
    ----------
    channel : int
        Index in the lielib_print(:) array to set. 0 => c_verbose.

    set : bool
        If set$ then set lielib_print(:). If unset$ then undo a previous set$.

    old_val : int
        Old value needed for set = unset$.
        This parameter is an input/output and is modified in-place.
        As an output, old_val: Saved value for set = set$.

    Returns
    -------
    old_val : int
        Old value needed for set = unset$.
        This parameter is an input/output and is modified in-place.
        As an output, old_val: Saved value for set = set$.
    """

def set_ptc_verbose(on: bool) -> None:
    """
    Wrapper for Fortran routine set_ptc_verbose

    Parameters
    ----------
    on : bool
    """

def set_pwd_ele(lat: _pybmad.LatStruct, mode0: _pybmad.NormalModesStruct, inductance: float) -> None:
    """
    Simulates the effect of potential well distortion by adjusting lat%ele(ix_pwd)%taylor(6)%term(2)%coef for an
    element in the lattice.  This element will apply a pz kick based on the z coordinate.
    Element is assumed to be at lat%ele(1).  The ibs_ring driver program
    inserts a taylor element into lat%ele(1) if set to perform pwd calculations.

    Parameters
    ----------
    lat : LatStruct
        lattice

    mode0 : NormalModesStruct
        .sig_z and .z.sige_e should be populated before calling this subroutine.

    inductance : float
        An inductance-like parameter describing the distortion of the potential well.
    """

def set_status_flags(stat: int) -> _pybmad.BookkeepingStateStruct:
    """
    Wrapper for Fortran routine set_status_flags

    Parameters
    ----------
    stat : int
        bookkeeping status. ok$, stale$, etc.

    Returns
    -------
    bookkeeping_state : BookkeepingStateStruct
    """

def set_tune(phi_a_set: float, phi_b_set: float, dk1: _pybmad.RealArray1D, eles: _pybmad.ElePointerStructArray1D, branch: _pybmad.BranchStruct, orb: _pybmad.CoordStructAlloc1D, print_err: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine set_tune

    Parameters
    ----------
    phi_a_set : float
        Horizontal set tune (radians)

    phi_b_set : float
        Vertical set tune (radians)

    dk1 : 1D array of float
        Relative amount to vary a quad in tuning. The variation will be proportional to dk1. Those quads with a
        positive dk1(i) will be varied as one group and the quads with negative dk1(i) will be varied as another
        group. The routine choose_quads_for_set_tune can be used to calculate values for dk1.

    eles : 1D array of ElePointerStruct
        eles(i).ele points to quadrupole corresponding to dk1(i).

    branch : BranchStruct
        Lattice branch to tune.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Q_tuned lattice branch

    orb : 1D array of CoordStruct
        If RF is off: Energy dE/E at which the tune is computed.
        This parameter is an input/output and is modified in-place.
        As an output, orb: New closed orbit.

    print_err : bool, optional
        Print error message if there is a problem? Default is True.

    Returns
    -------
    ok : bool
        Set True if everything is ok. False otherwise.
    """

def set_tune_via_group_knobs(phi_set: Sequence[float], branch: _pybmad.BranchStruct, group_knobs: Sequence[str], orb: _pybmad.CoordStructAlloc1D, print_err: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine set_tune_via_group_knobs

    Parameters
    ----------
    phi_set : 1D array of float (shape: 2)
        Set tunes (radians).

    branch : BranchStruct
        Lattice branch to tune.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Q_tuned lattice branch

    group_knobs : 1D array of str (shape: 2)
        Names of group knobs to vary.

    orb : 1D array of CoordStruct
        If RF is off: Energy dE/E at which the tune is computed.
        This parameter is an input/output and is modified in-place.
        As an output, orb: New closed orbit.

    print_err : bool, optional
        Print error message if there is a problem? Default is True.

    Returns
    -------
    ok : bool
        Set True if everything is ok. False otherwise.
    """

def set_twiss(branch: _pybmad.BranchStruct, twiss_ele: _pybmad.EleStruct, ix_ele: int, match_deta_ds: bool, err_flag: bool, print_err: bool | None = None) -> None:
    """
    Wrapper for Fortran routine set_twiss

    Parameters
    ----------
    branch : BranchStruct
        Branch to modify.

    twiss_ele : EleStruct
        Element with desired Twiss parameters.

    ix_ele : int
        Match branch.ele(ix_ele) Twiss to twiss_ele.

    match_deta_ds : bool
        If True, match deta_ds. If False, match etap.

    err_flag : bool
        Set True if there is an error. False otherwise.

    print_err : bool, optional
        Print an error message if there is an error? Default is True.
    """

def set_z_tune(branch: _pybmad.BranchStruct, z_tune: float, print_err: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine set_z_tune

    Parameters
    ----------
    branch : BranchStruct

    z_tune : float
        Longitudinal tune in radians (must be negative above transition).

    print_err : bool, optional
        Default is True. If False, suppress error messages

    Returns
    -------
    ok : bool, optional
        If present, returns true or false if set was successful. If not present, set_z_tune will bomb if tune
        could not be set.
    """

def settable_dep_var_bookkeeping(ele: _pybmad.EleStruct) -> None:
    """
    Subroutine to initialize dependent variables in an element.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.
    """

def setup_high_energy_space_charge_calc(calc_on: bool, branch: _pybmad.BranchStruct, n_part: float, mode: _pybmad.NormalModesStruct, beam_init: _pybmad.BeamInitStruct | None = None, closed_orb: _pybmad.CoordStructArray1D | None = None) -> None:
    """
    Routine to initialize constants needed by the ultra relativistic space charge
    tracking routine track1_high_energy_space_charge. This setup routine must be called if
    the lattice or any of the other input parameters are changed.

    Parameters used:
        a-mode emittance
        b-mode emittance
        sig_z bunch length
        sig_pz relative energy spread

    Parameters
    ----------
    calc_on : bool
        Turns on or off the space charge calculation.

    branch : BranchStruct
        Lattice for tracking.

    n_part : float
        Number of actual particles in a bunch. Used to compute the bunch charge.

    mode : NormalModesStruct
        Structure holding the beam info. Will be combined with info in beam_init.

    beam_init : BeamInitStruct, optional
        Structure holding beam info. Will be combined with info in mode.

    closed_orb : 1D array of CoordStruct, optional
        Closed orbit. If not present the closed orbit is taken to be zero.
    """

def sigma_mat_ptc_to_bmad(sigma_mat_ptc: Sequence[Sequence[float]], beta0: float) -> list[list[float]]:
    """
    Routine to convert a PTC sigma matrix to a Bmad sigma matrix.
    The conversion includes the conversion between Bmad and PTC time coordinate systems.

    Since PTC uses delta_E/P0c and Bmad uses delta_P/P0c coordinates, and since
    the relationship between delta_E and delta_P is nonlinear, this routine
    simplifies the calculation and assumes that the particle beta is constant
    over the range of particle energies.

    Parameters
    ----------
    sigma_mat_ptc : 2D array of float (shape: 6,6)
        PTC sigma matrix.

    beta0 : float
        Reference particle velocity

    Returns
    -------
    sigma_mat_bmad : 2D array of float (shape: 6,6)
        Bmad sigma matrix.
    """

def significant_difference(value1: float, value2: float, abs_tol: float | None = None, rel_tol: float | None = None) -> bool:
    """
    Wrapper for Fortran routine significant_difference

    Parameters
    ----------
    value1 : float
        First value.

    value2 : float
        Second value.

    abs_tol : float, optional
        Absolute tolerance. Default is 0.

    rel_tol : float, optional
        Relative tolerance. Default is 0.

    Returns
    -------
    is_different : bool
        Set True if the difference is significant. False otherwise.
    """

def skip_ele_blender(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine skip_ele_blender

    Parameters
    ----------
    ele : EleStruct

    Returns
    -------
    skip : bool
    """

def slice_lattice(lat: _pybmad.LatStruct, ele_list: str, do_bookkeeping: bool | None = None, set_phase_zero: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine slice_lattice

    Parameters
    ----------
    lat : LatStruct
        Lattice to slice.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with unwanted elements sliced out.

    ele_list : str
        List of elements to retain. See the documentation for the lat_ele_locator routine for the syntax of the
        list.

    do_bookkeeping : bool, optional
        Default is True. If false, the calling routine is responsible for: * Modifying lat.particle_start if
        needed. * Calculating Twiss functions.

    set_phase_zero : bool, optional
        Default is True. Set betatron phase to zero?

    Returns
    -------
    error : bool
        Set True if there is an error Set False if not.
    """

def soft_quadrupole_edge_kick(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, particle_at: int, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Routine to add the SAD "linear" soft edge (for finite f1 or f2).
    This routine assumes that the particle orbit has been rotated to the element reference frame.
    This routine is called with sad_mult and quadrupole elements.

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through

    param : LatParamStruct
        Tracking parameters.

    particle_at : int
        first_track_edge$, or second_track_edge$.

    orbit : CoordStruct
        Position before kick.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Position after kick.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the edge.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix with edge kick added on.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.
    """

def sol_quad_mat6_calc(ks_in: float, k1_in: float, tilt: float, length: float, ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine sol_quad_mat6_calc

    Parameters
    ----------
    ks_in : float

    k1_in : float

    tilt : float
        quadrupole tilt.

    length : float
        Sol_quad length.

    ele : EleStruct
        Sol_quad element.

    orbit : CoordStruct
        Orbit at beginning of the sol_quad.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the sol_quad.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix includeing the sol_quad.

    make_matrix : bool, optional
        Extend the matrix?
    """

def solve_psi_adaptive(t0: float, t1: float, p0: float, args: Sequence[float]) -> float:
    """
    Solve dpsi/dt for psi(t1) using adaptive steps and method:
      "Implicit Bulirsch-Stoer method of Bader and Deuflhard."

    The boundary condition p0 is psi(t0)

    Parameters
    ----------
    t0 : float
        initial time

    t1 : float
        final time

    p0 : float
        Boundary condition psi(t0)

    args : 1D array of float (shape: 1:8)
        Parameters.  See psi_prime comments for details.

    Returns
    -------
    p1 : float
        psi(t1)
    """

def solve_psi_fixed_steps(t0: float, t1: float, p0: float, args: Sequence[float], t: _pybmad.RealArray1D, p: _pybmad.RealArray1D) -> None:
    """
    Solve dpsi/dt for psi(t1) using fixed steps and method:
      "Implicit Bulirsch-Stoer method of Bader and Deuflhard."

    The boundary condition p0 is psi(t0).

    Number of steps is determined by SIZE(p).

    Parameters
    ----------
    t0 : float
        initial time

    t1 : float
        final time

    p0 : float
        Boundary condition psi(t0)

    args : 1D array of float (shape: 1:8)
        Parameters.  See psi_prime comments for details.

    t : 1D array of float
        Array of times from t0 to t1

    p : 1D array of float
        Array of psi evaluated at t(:)
    """

def sort_complex_taylor_terms(complex_taylor_in: _pybmad.ComplexTaylorStruct) -> _pybmad.ComplexTaylorStruct:
    """
    Subroutine to sort the complex_taylor terms from "lowest" to "highest" of
    a complex_taylor series.
    This subroutine is needed because what comes out of PTC is not sorted.

    Uses function complex_taylor_exponent_index to sort.

    Note: complex_taylor_sorted needs to have been initialized.
    Note: complex_taylor_sorted cannot be complex_taylor_in. That is it is not legal to write:
              call sort_complex_taylor_terms (this_complex_taylor, this_complex_taylor)

    Parameters
    ----------
    complex_taylor_in : ComplexTaylorStruct
        Unsorted complex_taylor series.

    Returns
    -------
    complex_taylor_sorted : ComplexTaylorStruct
        Sorted complex_taylor series.
    """

def space_charge_cathodeimages(mesh3d: _pybmad.Mesh3DStruct, direct_field_calc: bool | None = None, integrated_green_function: bool | None = None, image_method: int | None = None) -> None:
    """
    Wrapper for Fortran routine space_charge_cathodeimages

    Parameters
    ----------
    mesh3d : Mesh3dStruct

    direct_field_calc : bool, optional

    integrated_green_function : bool, optional

    image_method : int, optional
    """

def space_charge_freespace(mesh3d: _pybmad.Mesh3DStruct, direct_field_calc: bool | None = None, integrated_green_function: bool | None = None) -> None:
    """
    Wrapper for Fortran routine space_charge_freespace

    Parameters
    ----------
    mesh3d : Mesh3dStruct

    direct_field_calc : bool, optional

    integrated_green_function : bool, optional
    """

def space_charge_rectpipe(mesh3d: _pybmad.Mesh3DStruct, apipe: float, bpipe: float, direct_field_calc: bool | None = None, integrated_green_function: bool | None = None) -> None:
    """
    Wrapper for Fortran routine space_charge_rectpipe

    Parameters
    ----------
    mesh3d : Mesh3dStruct

    apipe : float

    bpipe : float

    direct_field_calc : bool, optional

    integrated_green_function : bool, optional
    """

class SpinConcatLinearMaps:
    """spin_concat_linear_maps return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def map1(self) -> _pybmad.SpinOrbitMap1Struct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spin_concat_linear_maps(branch: _pybmad.BranchStruct, n1: int, n2: int, map1_ele: _pybmad.SpinOrbitMap1StructArray1D | None = None, orbit: _pybmad.CoordStructArray1D | None = None, excite_zero: Sequence[str] | None = None) -> SpinConcatLinearMaps:
    """
    Wrapper for Fortran routine spin_concat_linear_maps

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch.

    n1 : int
        Starting element index. Start at element downstream end.

    n2 : int
        Ending element index. End at element downstream end

    map1_ele : 1D array of SpinOrbitMap1Struct, optional
        Individual spin/orbit maps.

    orbit : 1D array of CoordStruct, optional
        Reference orbit used if maps must be created.

    excite_zero : 1D array of str (shape: 3), optional
        Three lists of elements where ds_vec/dr_vec terms are zeroed.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.

    map1 : SpinOrbitMap1Struct
        Map with element spin/orbit maps concatenated.
    """

def spin_depolarization_rate(branch: _pybmad.BranchStruct, match_info: _pybmad.SpinMatchingStructArray1D, rad_int_by_ele: _pybmad.RadIntAllEleStruct) -> float:
    """
    Wrapper for Fortran routine spin_depolarization_rate

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch the beam is going through.

    match_info : 1D array of SpinMatchingStruct

    rad_int_by_ele : RadIntAllEleStruct
        Element-by-element radiation integrals.

    Returns
    -------
    depol_rate : float
        Depolarization rate (1/sec). Will be positive.
    """

class SpinDnDpzFromMat8:
    """spin_dn_dpz_from_mat8 return type"""

    @property
    def error(self) -> bool: ...

    @property
    def dn_dpz(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spin_dn_dpz_from_mat8(mat_1turn: Sequence[Sequence[float]], dn_dpz_partial: Sequence[Sequence[float]] | None = None) -> SpinDnDpzFromMat8:
    """
    Wrapper for Fortran routine spin_dn_dpz_from_mat8

    Parameters
    ----------
    mat_1turn : 2D array of float (shape: 8,8)
        Spin-orbital matrix.

    dn_dpz_partial : 2D array of float (shape: 3,3), optional
        dn_dpz_partial(i,:) is dn_dpz with only one osccilation mode "excited". So dn_dpz_partial(1,:) represents
        a-mode excitation, etc.

    Returns
    -------
    error : bool
        Set True if there is an error. False otherwise.

    dn_dpz : 1D array of float (shape: 3)
        dn_dpz (l,n,m) coordinates.
    """

class SpinDnDpzFromQmap:
    """spin_dn_dpz_from_qmap return type"""

    @property
    def error(self) -> bool: ...

    @property
    def dn_dpz(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spin_dn_dpz_from_qmap(orb_mat: Sequence[Sequence[float]], q_map: Sequence[Sequence[float]], dn_dpz_partial: Sequence[Sequence[float]], dn_dpz_partial2: Sequence[Sequence[float]], n0: Sequence[float] | None = None) -> SpinDnDpzFromQmap:
    """
    Wrapper for Fortran routine spin_dn_dpz_from_qmap

    Parameters
    ----------
    orb_mat : 2D array of float (shape: 6,6)
        1-turn orbital matrix.

    q_map : 2D array of float (shape: 0:3,0:6)
        1-turn spin linear quaternion map.

    dn_dpz_partial : 2D array of float (shape: 3,3)
        ) is dn_dpz with only one osccilation mode "excited". So dn_dpz_partial(1,:) represents a-mode excitation,
        etc.

    dn_dpz_partial2 : 2D array of float (shape: 3,3)
        ) is dn_dpz with only two osccilation modes "excited". So dn_dpz_partial(1,:) represents b-mode and c-mode
        excitation without the a-mode, etc.

    n0 : 1D array of float (shape: 3), optional
        3,0).

    Returns
    -------
    error : bool
        Set True if there is an error. False otherwise.

    dn_dpz : 1D array of float (shape: 3)
        dn_dpz.
    """

def spin_map1_normalize(spin1: Sequence[Sequence[float]]) -> None:
    """
    Wrapper for Fortran routine spin_map1_normalize

    Parameters
    ----------
    spin1 : 2D array of float (shape: 0:3,0:6)
        Unnormalized spin map.
        This parameter is an input/output and is modified in-place.
        As an output, spin1: Normalized spin map.
    """

class SpinMat8ResonanceStrengths:
    """spin_mat8_resonance_strengths return type"""

    @property
    def xi_sum(self) -> float: ...

    @property
    def xi_diff(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spin_mat8_resonance_strengths(orb_evec: Sequence[complex], mat8: Sequence[Sequence[float]]) -> SpinMat8ResonanceStrengths:
    """
    Wrapper for Fortran routine spin_mat8_resonance_strengths

    Parameters
    ----------
    orb_evec : 1D array of complex (shape: 6)
        Orbital eigenvector.

    mat8 : 2D array of float (shape: 6,6)
        Spin/orbital matrix.

    Returns
    -------
    xi_sum : float
        Sum resonance strength.

    xi_diff : float
        Difference resonance strength.
    """

class SpinMatToEigen:
    """spin_mat_to_eigen return type"""

    @property
    def orb_eval(self) -> list[complex]: ...

    @property
    def orb_evec(self) -> list[list[complex]]: ...

    @property
    def n0(self) -> list[float]: ...

    @property
    def spin_evec(self) -> list[list[complex]]: ...

    @property
    def error(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spin_mat_to_eigen(orb_mat: Sequence[Sequence[float]], spin_map: Sequence[Sequence[float]]) -> SpinMatToEigen:
    """
    Wrapper for Fortran routine spin_mat_to_eigen

    Parameters
    ----------
    orb_mat : 2D array of float (shape: 6,6)
        Orbital matrix.

    spin_map : 2D array of float (shape: 0:3,0:6)
        Quaternion 0th & 1st order map.

    Returns
    -------
    orb_eval : 1D array of complex (shape: 6)
        Eigenvalues.

    orb_evec : 2D array of complex (shape: 6,6)
        Orbital eigenvectors. orb_evec(j,:) is the j^th vector.

    n0 : 1D array of float (shape: 3)
        n_0 invariant spin

    spin_evec : 2D array of complex (shape: 6,3)
        Spin eigenvectors. spin_evec(j,:) is the j^th vector.

    error : bool
        Set true if there is an error. False otherwise.
    """

def spin_omega(field: _pybmad.EmFieldStruct, orbit: _pybmad.CoordStruct, sign_z_vel: int, phase_space_coords: bool | None = None) -> list[float]:
    """
    Wrapper for Fortran routine spin_omega

    Parameters
    ----------
    field : EmFieldStruct

    orbit : CoordStruct

    sign_z_vel : int

    phase_space_coords : bool, optional

    Returns
    -------
    omega : 1D array of float (shape: 3)
    """

class SpinQuatResonanceStrengths:
    """spin_quat_resonance_strengths return type"""

    @property
    def xi_sum(self) -> float: ...

    @property
    def xi_diff(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spin_quat_resonance_strengths(orb_evec: Sequence[complex], spin_q: Sequence[Sequence[float]]) -> SpinQuatResonanceStrengths:
    """
    Wrapper for Fortran routine spin_quat_resonance_strengths

    Parameters
    ----------
    orb_evec : 1D array of complex (shape: 6)
        Orbital eigenvector.

    spin_q : 2D array of float (shape: 0:3,0:6)
        First order spin map.

    Returns
    -------
    xi_sum : float
        Sum resonance strength.

    xi_diff : float
        Difference resonance strength.
    """

def spin_taylor_to_linear(spin_taylor: _pybmad.TaylorStructArray1D, normalize: bool, dref_orb: Sequence[float], is_on: bool) -> list[list[float]]:
    """
    Wrapper for Fortran routine spin_taylor_to_linear

    Parameters
    ----------
    spin_taylor : 1D array of TaylorStruct (shape: 0:3)
        Taylor spin map.

    normalize : bool
        If True, normalize the linear map.

    dref_orb : 1D array of float (shape: 6)
        Change in Reference orbit: output_map1_ref - input_taylor_ref.

    is_on : bool
        Is map turned on? If not spin_map1 will be the unit map.

    Returns
    -------
    spin_map1 : 2D array of float (shape: 0:3,0:6)
        First order spin map.
    """

def spinor_to_polar(spinor: Sequence[complex]) -> _pybmad.SpinPolarStruct:
    """
    Wrapper for Fortran routine spinor_to_polar

    Parameters
    ----------
    spinor : 1D array of complex (shape: 2)
        Spinor

    Returns
    -------
    polar : SpinPolarStruct
        The resultant Unitary Vector in polar coordinates
    """

def spinor_to_vec(spinor: Sequence[complex]) -> list[float]:
    """
    Wrapper for Fortran routine spinor_to_vec

    Parameters
    ----------
    spinor : 1D array of complex (shape: 2)
        Spinor

    Returns
    -------
    vec : 1D array of float (shape: 3)
        spin vector in cartesian coordinates
    """

def spline_fit_orbit(start_orb: _pybmad.CoordStruct, end_orb: _pybmad.CoordStruct, spline_x: Sequence[float], spline_y: Sequence[float]) -> None:
    """
    Wrapper for Fortran routine spline_fit_orbit

    Parameters
    ----------
    start_orb : CoordStruct
        Starting coords.

    end_orb : CoordStruct
        Ending coords.

    spline_x : 1D array of float (shape: 0:3)
        Spline coefs for the horizontal trajectory.

    spline_y : 1D array of float (shape: 0:3)
        Spline coefs for vertical trajectory.
    """

def split_expression_string(expr: str, width: int, indent: int, break_str: str | None = None) -> _pybmad.CharacterAlloc1D:
    """
    Routine to break an expression into a number of lines for a nicer display.
    Used when printing expressions.

    Parameters
    ----------
    expr : str
        String containing the expression.

    width : int
        Maximum width of split expression.

    indent : int
        If positive: Number of spaces to indent for every line after the first. If negative: No indentation but
        first line is shortened by |indent|.

    break_str : str, optional
        If present, only break lines at places where this string is.

    Returns
    -------
    lines : 1D array of str
        Split expression.
    """

class SplitLat:
    """split_lat return type"""

    @property
    def ix_split(self) -> int: ...

    @property
    def split_done(self) -> bool: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def split_lat(lat: _pybmad.LatStruct, s_split: float, ix_branch: int, add_suffix: bool | None = None, check_sanity: bool | None = None, save_null_drift: bool | None = None, choose_max: bool | None = None, ix_insert: int | None = None) -> SplitLat:
    """
    Wrapper for Fortran routine split_lat

    Parameters
    ----------
    lat : LatStruct
        Original lat structure.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Modified lat structure.

    s_split : float
        Position at which lat.branch(ix_branch) is to be split.

    ix_branch : int
        Index of lat.branch(:) to use.

    add_suffix : bool, optional
        If True (default) add '#1' and '#2" suffixes to the split elements.

    check_sanity : bool, optional
        If True (default) then call lat_sanity_check after the split to make sure everything is ok.

    save_null_drift : bool, optional
        Save a copy of a drift to be split as a null_ele? This is useful when superpositions are done. See
        add_superimpose for more info. Default is False.

    choose_max : bool, optional
        If no splitting of an element is needed, that is, s_split is at an element boundary, there can be multiple
        possible values for ix_split if there exist zero length elements at the split point. If choose_max = True,
        ix_split will be chosen to be the maximum possible index and if choose_max = False ix_split will be chosen
        to be the minimal possible index. If s_split is not at an element boundary, the setting of choose_max is
        immaterial. If ix_insert is present, the default value of choose_max is set to give the closest element to
        ix_insert. If ix_insert is not present, the default value of choose_max is False.

    ix_insert : int, optional
        Element index near the point to be split. ix_insert is useful in the case where there is a patch with a
        negative length which can create an ambiguity as to where to do the split In this case ix_insert will
        remove the ambiguity. Also useful to ensure where to split if there are elements with zero length nearby.
        Ignored if negative.

    Returns
    -------
    ix_split : int
        Index of element just before the s = s_split point.

    split_done : bool
        True if lat was split.

    err_flag : bool, optional
        Set true if there is an error, false otherwise.
    """

def sprint_spin_taylor_map(ele: _pybmad.EleStruct, start_orbit: Sequence[float] | None = None) -> None:
    """
    Wrapper for Fortran routine sprint_spin_taylor_map

    Parameters
    ----------
    ele : EleStruct
        Element to form map for.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with map.

    start_orbit : 1D array of float (shape: 6), optional
        Reference orbit for the map. Default is zero orbit.
    """

def sr_longitudinal_wake_particle(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to apply the short-range wake longitudinal component kick to a particle and then add
    to the existing longitudinal wake the contribution from the particle.

    Parameters
    ----------
    ele : EleStruct
        Element with wakes.

    orbit : CoordStruct
        Particle coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: coords after the kick.
    """

def sr_transverse_wake_particle(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Subroutine to apply the short-range wake transverse component of the kick to a particle and then add
    to the existing transverse wake the contribution from the particle.

    Parameters
    ----------
    ele : EleStruct
        Element with wakes.

    orbit : CoordStruct
        Starting particle coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending particle coords.
    """

def sr_z_long_wake(ele: _pybmad.EleStruct, bunch: _pybmad.BunchStruct, z_ave: float) -> None:
    """
    Subroutine to apply the short-range z-wake kick to a particle.

    Parameters
    ----------
    ele : EleStruct
        Element with wake.

    bunch : BunchStruct
        Bunch before wake applied.

    z_ave : float
        Average z-position of all live particles.
    """

def srdt_calc(lat: _pybmad.LatStruct, order: int, n_slices_gen_opt: int | None = None, n_slices_sxt_opt: int | None = None, per_ele_out: _pybmad.SummationRdtStructAlloc1D | None = None) -> _pybmad.SummationRdtStruct:
    """
    Calculate summation RDT terms up to order=1 or order=2 while slicing sextupoles
    n_slices_sxt_opt times and all other elements n_slices_gen_opt times.

    These formulas are documented in "The Sextupole Scheme for the Swiss Light Source (SLS): An Analytic Approach"
    by Johan Bengtsson.  SLS Note 9/97.

    The 2nd order formulas are documented in "Second-order driving terms due to sextupoles and
    chromatic effects of quadrupoles" by Chun-xi Wang.  AOP-TN-2009-020.

    Parameters
    ----------
    lat : LatStruct
        lattice with Twiss parameters calculated.

    order : int
        1 to calculate only first order terms.  2 to also calculate 2nd order terms.

    n_slices_gen_opt : int, optional
        number of times to slice elements other than sextupoles.  Default is 10.

    n_slices_sxt_opt : int, optional
        nubmer of times to slice sextupoles.  Default is 20.

    Returns
    -------
    srdt_sums : SummationRdtStruct
        contains complex RDT strengths.
    """

def srdt_lsq_solution(lat: _pybmad.LatStruct, var_indexes: _pybmad.IntArray1D, n_slices_gen_opt: int | None = None, n_slices_sxt_opt: int | None = None, chrom_set_x_opt: float | None = None, chrom_set_y_opt: float | None = None, weight_in: Sequence[float] | None = None) -> _pybmad.RealAlloc1D:
    """
    chrom_set_x_opt, chrom_set_y_opt, weight_in)

    Given lat, finds K2 moments that set the chromaticity and zeros-out the real
    and complex parts of the first order driving terms, that minimizes the sum of the squares
    of the K2 moments.  i.e. the weakest sextupole scheme that sets chromaticity
    and zeros out the first order terms.

    Note:  This subroutine does not, in its present form, work well with knobs, overlays, or in lattices where
           multiple elements have the same name.

    This subroutine assumes that Nsext > 18.

    Parameters
    ----------
    lat : LatStruct
        lattice with Twiss parameters calculated.

    var_indexes : 1D array of int
        indexes in lat.ele that are K2 variables.  Must be sorted smallest index to largest index.

    n_slices_gen_opt : int, optional
        number of times to slice elements other than sextupoles.  Default is 10.

    n_slices_sxt_opt : int, optional
        nubmer of times to slice sextupoles.  Default is 20.

    chrom_set_x_opt : float, optional
        what to set x chromaticity to.  Default zero.

    chrom_set_y_opt : float, optional
        what to set y chromaticity to.  Default zero.

    weight_in : 1D array of float (shape: 10), optional
        moment weights. Terms are: [wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002, wgt_h21000,
        wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200, If present, any terms equal to zero are given default
        values which is 1.0e4 for wgt_chrom_x and wgt_chrom_y and is 1.0 for everything else.

    Returns
    -------
    ls_soln : 1D array of float
        contains K2 for the indexes in var_indexes
    """

def start_branch_at(lat: _pybmad.LatStruct, ele_start: str, move_end_marker: bool) -> bool:
    """
    Wrapper for Fortran routine start_branch_at

    Parameters
    ----------
    lat : LatStruct
        Lattice to modify.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Modified lattice.

    ele_start : str
        Start element. Ele_start will identify the lattice branch to modify.

    move_end_marker : bool
        If True then the end marker (if it is present) will be shifted like any other element. False means that
        the end marker will stay at the end.

    Returns
    -------
    error : bool
        Set True if there is an error Set False if not.
    """

def stream_ele_end(physical_end: int, ele_orientation: int) -> int:
    """
    Wrapper for Fortran routine stream_ele_end

    Parameters
    ----------
    physical_end : int
        entrance_end$, exit_end$, surface$, etc.

    ele_orientation : int
        Either 1 = Normal or -1 = element reversed.

    Returns
    -------
    stream_end : int
        upstream_end$, downstream_end$, or set equal to physical_end if physical_end is neither entrance_end$ nor
        exit_end$
    """

def string_attrib(attrib_name: str, ele: _pybmad.EleStruct) -> str:
    """
    Routine to return the value of a string attribute of a lattice element.
    This routine is useful when attrib_name is specified by the program user.

    For example:
      call string_attrib ('NAME', ele, attrib_value)  ! Will return attrib_value = ele%name

    Parameters
    ----------
    attrib_name : str
        Name of the type of element attribute.

    ele : EleStruct
        Lattice element.

    Returns
    -------
    attrib_value : str
        The string associated with the attribute.
    """

class StrongBeamSigmaCalc:
    """strong_beam_sigma_calc return type"""

    @property
    def sigma(self) -> list[float]: ...

    @property
    def bbi_const(self) -> float: ...

    @property
    def dsigma_ds(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def strong_beam_sigma_calc(ele: _pybmad.EleStruct, s_pos: float) -> StrongBeamSigmaCalc:
    """
    Wrapper for Fortran routine strong_beam_sigma_calc

    Parameters
    ----------
    ele : EleStruct
        Beambeam element.

    s_pos : float
        Longitudinal position in lab coords of slice (used with hourglass effect correction).

    Returns
    -------
    sigma : 1D array of float (shape: 2)
        Strong beam x,y sigmas.

    bbi_const : float
        BBI kick scale factor.

    dsigma_ds : 1D array of float (shape: 2)
        sig_x and sig_y longitudinal derivatives.
    """

def strong_beam_strength(ele: _pybmad.EleStruct) -> float:
    """
    Wrapper for Fortran routine strong_beam_strength

    Parameters
    ----------
    ele : EleStruct
        Beambeam element.

    Returns
    -------
    strength : float
        Strong beam strength.
    """

class SurfaceGridDisplacement:
    """surface_grid_displacement return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def z(self) -> float: ...

    @property
    def dz_dxy(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def surface_grid_displacement(ele: _pybmad.EleStruct, x: float, y: float, extend_grid: bool | None = None) -> SurfaceGridDisplacement:
    """
    Routine to add in the z displacement defined by the grid

    Parameters
    ----------
    ele : EleStruct
        Element containing the grid

    x : float
        Photon coords at surface.

    y : float
        Photon coords at surface.

    extend_grid : bool, optional
        If (x,y) past grid pretend (x,y) is at grid boundary. Default is False.

    Returns
    -------
    err_flag : bool
        Set True if there is a problem.

    z : float
        surface height at (x, y).

    dz_dxy : 1D array of float (shape: 2), optional
        Surface slope at (x, y).
    """

class SwitchAttribValueName:
    """switch_attrib_value_name return type"""

    @property
    def is_default(self) -> bool: ...

    @property
    def name_list(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def attrib_val_name(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def switch_attrib_value_name(attrib_name: str, attrib_value: float, ele: _pybmad.EleStruct) -> SwitchAttribValueName:
    """
    is_default, name_list) result (attrib_val_name)

    Routine to return the name corresponding to the value of a given switch attribute.

    This routine is for "switch" and "species" attributes. For example, the "aperture_type" attribute
    can have value names of "Entrance_End", "Exit_End", etc.

    Optionally, this routine can determine if the attribute value corresponds
    to the default value. That is, the value that the attribute would have if
    not specified in the lattice file.

    Use the routine attribute_type to first test if the type of the attribute
    corresponds to is_switch$.

    Parameters
    ----------
    attrib_name : str
        Name of the attribute. Must be upper case.

    attrib_value : float
        Value of the attribute.

    ele : EleStruct
        Lattice element that the attribute is contained in. Generally only needed to determine the default value.

    Returns
    -------
    attrib_val_name : str
        Name corresponding to the value. Set to null_name$ if there is a problem.

    is_default : bool, optional
        If True then the value of the attiribute corresponds to the default value. If this argument is present,
        the ele argument must also be present.

    name_list : 1D array of str, optional
        List of names the switch can take. Deallocated if there is an error.
    """

def symp_lie_bmad(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None, offset_ele: bool | None = None) -> _pybmad.TrackStruct:
    """
    Wrapper for Fortran routine symp_lie_bmad

    Parameters
    ----------
    ele : EleStruct
        Element with transfer matrix
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with transfer matrix.

    param : LatParamStruct
        Parameters are needed for some elements.

    orbit : CoordStruct
        Coordinates at the beginning of element.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Coordinates at the end of element.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix propagated through the element.

    make_matrix : bool, optional
        If True then make the 6x6 transfer matrix.

    offset_ele : bool, optional
        Offset the element using ele.value(x_offset$), etc. Default is True.

    Returns
    -------
    track : TrackStruct, optional
        Structure holding the track information. When tracking through multiple elements, the trajectory in an
        element is appended to the existing trajectory. To reset: Set track.n_pt = -1.
    """

class T6ToB123:
    """t6_to_b123 return type"""

    @property
    def B1(self) -> list[list[float]]: ...

    @property
    def B2(self) -> list[list[float]]: ...

    @property
    def B3(self) -> list[list[float]]: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def t6_to_b123(t6: Sequence[Sequence[float]], abz_tunes: Sequence[float]) -> T6ToB123:
    """
    This decomposes the one-turn matrix according to Equation 56 from
    "Alternative approach to general coupled linear optics" by A. Wolski. PRSTAB.

    Note that a sigma matrix can be assembeled from:  sigma = B1*emit_a + B2*emit_b + B3*emit_c

    Parameters
    ----------
    t6 : 2D array of float (shape: 6,6)
        1-turn transfer matrix.  RF assumed to be on.

    abz_tunes : 1D array of float (shape: 3)
        a-mode and b-mode tunes.  Used to order eigensystem.

    Returns
    -------
    B1 : 2D array of float (shape: 6,6)
        Beta matrix associated with a-mode.

    B2 : 2D array of float (shape: 6,6)
        Beta matrix associated with b-mode.

    B3 : 2D array of float (shape: 6,6)
        Beta matrix associated with c-mode.

    err_flag : bool
        Set True if there is an error. False otherwise
    """

def taper_mag_strengths(lat: _pybmad.LatStruct, ref_lat: _pybmad.LatStruct | None = None, except_: str | None = None, err_flag: bool | None = None) -> None:
    """
    Wrapper for Fortran routine taper_mag_strengths

    Parameters
    ----------
    lat : LatStruct
        Lattice to vary.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with magnet strengths varied.

    ref_lat : LatStruct, optional
        Reference lattice. If not present, lat will be used as the ref.

    except : str, optional
        List of elements not to vary.

    err_flag : bool, optional
    """

class TargetMinMaxCalc:
    """target_min_max_calc return type"""

    @property
    def y_min(self) -> float: ...

    @property
    def y_max(self) -> float: ...

    @property
    def phi_min(self) -> float: ...

    @property
    def phi_max(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def target_min_max_calc(r_corner1: Sequence[float], r_corner2: Sequence[float], y_min: float, y_max: float, phi_min: float, phi_max: float, initial: bool | None = None) -> TargetMinMaxCalc:
    """
    Routine to calculate the min/max values for (y, phi).
    min/max values are cumulative.

    Parameters
    ----------
    r_corner1 : 1D array of float (shape: 3)
        In target coords: A corner of the target. Must be normalized to 1.

    r_corner2 : 1D array of float (shape: 3)
        In target coords: Adjacent corner of the target. Must be normalized to 1.

    y_min : float
        min/max values. Only needed if initial = False.

    y_max : float
        min/max values. Only needed if initial = False.

    phi_min : float
        min/max values. Only needed if initial = False.

    phi_max : float
        min/max values. Only needed if initial = False.

    initial : bool, optional
        If present and True then this is the first edge for computation.

    Returns
    -------
    y_min : float
        min/max values. Only needed if initial = False.

    y_max : float
        min/max values. Only needed if initial = False.

    phi_min : float
        min/max values. Only needed if initial = False.

    phi_max : float
        min/max values. Only needed if initial = False.
    """

class TargetRotMats:
    """target_rot_mats return type"""

    @property
    def w_to_target(self) -> list[list[float]]: ...

    @property
    def w_to_ele(self) -> list[list[float]]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def target_rot_mats(r_center: Sequence[float]) -> TargetRotMats:
    """
    Routine to calculate the rotation matrices between ele coords and "target" coords.
    By definition, in target coords r_center = [0, 0, 1].

    Parameters
    ----------
    r_center : 1D array of float (shape: 3)
        In lab coords: Center of target relative to phton emission point.

    Returns
    -------
    w_to_target : 2D array of float (shape: 3,3)
        Rotation matrix from ele to target coords.

    w_to_ele : 2D array of float (shape: 3,3)
        Rotation matrix from target to ele coords.
    """

def taylor_equal_taylor(taylor1: _pybmad.TaylorStruct, taylor2: _pybmad.TaylorStruct) -> None:
    """
    Wrapper for Fortran routine taylor_equal_taylor

    Parameters
    ----------
    taylor1 : TaylorStruct

    taylor2 : TaylorStruct
    """

def taylor_inverse(taylor_in: _pybmad.TaylorStructArray1D, taylor_inv: _pybmad.TaylorStructArray1D) -> bool:
    """
    Subroutine to invert a taylor map. Since the inverse map is truncated, it is not exact.

    Parameters
    ----------
    taylor_in : 1D array of TaylorStruct
        Input taylor map.

    taylor_inv : 1D array of TaylorStruct
        Inverted taylor map.

    Returns
    -------
    err : bool, optional
        Set True if there is no inverse. If not present then print an error message.
    """

def taylor_propagate1(orb_taylor: _pybmad.TaylorStructArray1D, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, ref_in: _pybmad.CoordStruct | None = None, spin_taylor: _pybmad.TaylorStructArray1D | None = None) -> bool:
    """
    Subroutine to track (symplectic integration) a orbital map, and optionally a spin map, through an element.
    The spin tracking is only done if spin_taylor is present and bmad_com%spin_tracking_on = T.
    The alternative routine, if ele has a taylor map, is concat_taylor.

    This routine will fail if there is no corresponding ptc fibre for this
    element. In general, the transfer_map_calc routine should be used instead.

    Parameters
    ----------
    orb_taylor : 1D array of TaylorStruct
        Map to be tracked
        This parameter is an input/output and is modified in-place.
        As an output, orb_taylor: Map through element.

    ele : EleStruct
        Element to track through

    param : LatParamStruct

    ref_in : CoordStruct, optional
        Particle to be tracked. Must be present if the particle to be tracked is not the reference particle or if
        the direction of propagation is backwards.

    spin_taylor : 1D array of TaylorStruct, optional
        Spin map to be tracked
        This parameter is an input/output and is modified in-place.
        As an output, spin_taylor: Tracked spin map.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def taylor_to_mad_map(taylor: _pybmad.TaylorStructArray1D, energy: _pybmad.MadEnergyStruct) -> _pybmad.MadMapStruct:
    """
    Subroutine to convert a Taylor map to a mad order 2 map.
    If any of the Taylor terms have order greater than 2 they are ignored.

    Parameters
    ----------
    taylor : 1D array of TaylorStruct
        Taylor map.

    energy : MadEnergyStruct
        Energy numbers.

    Returns
    -------
    map : MadMapStruct
        Order 2 map.
    """

def taylors_equal_taylors(taylor1: _pybmad.TaylorStructArray1D, taylor2: _pybmad.TaylorStructArray1D) -> None:
    """
    Wrapper for Fortran routine taylors_equal_taylors

    Parameters
    ----------
    taylor1 : 1D array of TaylorStruct

    taylor2 : 1D array of TaylorStruct
    """

def tilt_coords(tilt_val: float, coord: _pybmad.RealArray1D, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tilt_coords

    Parameters
    ----------
    tilt_val : float
        Tilt value (could be the roll value for a bend)

    coord : 1D array of float
        Coordinates of particle before rotation.
        This parameter is an input/output and is modified in-place.
        As an output, coord: Coordinates of particle after rotation.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before tilt.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix transfer matrix after tilt applied.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def tilt_coords_photon(tilt_val: float, coord: _pybmad.RealArray1D, w_mat: Sequence[Sequence[float]] | None = None) -> None:
    """
    Wrapper for Fortran routine tilt_coords_photon

    Parameters
    ----------
    tilt_val : float
        Tilt value (could be the roll value for a bend)

    coord : 1D array of float
        Coordinates of particle before rotation.
        This parameter is an input/output and is modified in-place.
        As an output, coord: Coordinates of particle after rotation.

    w_mat : 2D array of float (shape: 3,3), optional
        Rotation matrix before tilt.
        This parameter is an input/output and is modified in-place.
        As an output, w_mat: Rotation matrix after tilt.
    """

def tilt_mat6(mat6: Sequence[Sequence[float]], tilt: float) -> None:
    """
    Wrapper for Fortran routine tilt_mat6

    Parameters
    ----------
    mat6 : 2D array of float (shape: 6,6)
        Untilted matrix.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Tilted matrix.

    tilt : float
        Tilt angle.
    """

class ToEtaReading:
    """to_eta_reading return type"""

    @property
    def reading(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def to_eta_reading(eta_actual: _pybmad.RealArray1D, ele: _pybmad.EleStruct, axis: int, add_noise: bool) -> ToEtaReading:
    """
    Compute the measured dispersion reading given the true dispersion and the
    monitor offsets, noise, etc.

    This routine will only give a nonzero reading for Bmad markers,
    monitors, and instruments.

    Parameters
    ----------
    eta_actual : 1D array of float
        Actual (eta_x, eta_y) dispersion.

    ele : EleStruct
        Element where the orbit is measured.

    axis : int
        x_plane$ or y_plane$

    add_noise : bool
        If True add noise to the reading

    Returns
    -------
    reading : float
        BPM reading

    err : bool
        Set True if there is an error. False otherwise.
    """

class ToFieldmapCoords:
    """to_fieldmap_coords return type"""

    @property
    def x(self) -> float: ...

    @property
    def y(self) -> float: ...

    @property
    def z(self) -> float: ...

    @property
    def cos_ang(self) -> float: ...

    @property
    def sin_ang(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def to_fieldmap_coords(ele: _pybmad.EleStruct, local_orb: _pybmad.CoordStruct, s_body: float, ele_anchor_pt: int, r0: Sequence[float], curved_ref_frame: bool) -> ToFieldmapCoords:
    """
    x, y, z, cos_ang, sin_ang, err_flag)

    Routine to return the (x,y,s) position relative to a field map.

    Parameters
    ----------
    ele : EleStruct
        Element being tracked through.

    local_orb : CoordStruct
        Particle orbit. Must be in local element coordinates.

    s_body : float
        Longitudinal position relative to the entrance end of the element.

    ele_anchor_pt : int
        anchor point of the field map (anchor_beginning$, anchor_center$, or anchor_end$).

    r0 : 1D array of float (shape: 3)
        origin point of the fieldmap.

    curved_ref_frame : bool
        If the element is a bend: Does the field map follow the bend reference coords?

    Returns
    -------
    x : float
        Coords relative to the field map.

    y : float
        Coords relative to the field map.

    z : float
        Coords relative to the field map.

    cos_ang : float
        cos and sin of coordinate rotation angle.

    sin_ang : float
        cos and sin of coordinate rotation angle.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

class ToOrbitReading:
    """to_orbit_reading return type"""

    @property
    def reading(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def to_orbit_reading(orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, axis: int, add_noise: bool) -> ToOrbitReading:
    """
    Calculate the measured reading on a bpm given the actual orbit and the
    BPM's offsets, noise, etc.

    This routine will only give a nonzero reading for Bmad markers,
    monitors, and instruments.

    Parameters
    ----------
    orb : CoordStruct
        Orbit position at BPM.

    ele : EleStruct
        Element where the orbit is measured.

    axis : int
        x_plane$ or y_plane$

    add_noise : bool
        If True add noise to the reading

    Returns
    -------
    reading : float
        BPM reading

    err : bool
        Set True if there is no valid reading. For example, if ele.is_on = False.
    """

class ToPhaseAndCouplingReading:
    """to_phase_and_coupling_reading return type"""

    @property
    def reading(self) -> _pybmad.BpmPhaseCouplingStruct: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def to_phase_and_coupling_reading(ele: _pybmad.EleStruct, add_noise: bool) -> ToPhaseAndCouplingReading:
    """
    Find the measured coupling values given the actual ones

    This routine will only give a nonzero reading for Bmad markers,
    monitors, and instruments.

    Parameters
    ----------
    ele : EleStruct
        Element where phase is measured.

    add_noise : bool
        If True add noise to the reading

    Returns
    -------
    reading : BpmPhaseCouplingStruct
        K and Cbar coupling parameters

    err : bool
        Set True if there is an error. False otherwise.
    """

def to_photon_angle_coords(orb_in: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> _pybmad.CoordStruct:
    """
    Routine to convert from standard photon coords to "angle" coords defined as:
          x, angle_x, y, angle_y, z, E-E_ref

    Parameters
    ----------
    orb_in : CoordStruct
        orbit in standard photon coords.

    ele : EleStruct
        Reference element (generally the detector element.)

    Returns
    -------
    orb_out : CoordStruct
        Transformed coordinates.
    """

def to_surface_coords(lab_orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> _pybmad.CoordStruct:
    """
    Wrapper for Fortran routine to_surface_coords

    Parameters
    ----------
    lab_orbit : CoordStruct
        Photon position in laboratory coords.

    ele : EleStruct
        Detector element.

    Returns
    -------
    surface_orbit : CoordStruct
        Photon position in element body coordinates.
    """

def touschek_lifetime(mode: _pybmad.NormalModesStruct, lat: _pybmad.LatStruct) -> float:
    """
    Calculates the touschek lifetime for a lattice by calling touschek_rate1
    for each element.
    The loss rate at each element is averaged over one turn to obtain the lifetime.

    This function assumes that the twiss parameters and closed orbit have
    been calculated, and that mode has been populated.

    This subroutine assumes a fixed momentum aperture.  The loss rate at each element
    uses the same momentum aperture, mode%pz_aperture.

    A common way to call this function is to first populate mode using
    radiation integrals.  If an ideal lattice is used, the vertical
    emittance must also be set to a reasonable value.  If the vertical
    emittance is due only to quantum excitation, then it will likely be
    several orders of magnitude smaller than any real physical situation, in which
    case the integral in this function will have problems converging.

    In addition to setting mode, also set lat%param%n_part to the number of particles
    per bunch.

    Parameters
    ----------
    mode : NormalModesStruct
        beam properties

    lat : LatStruct
        Accelerator Lattice

    Returns
    -------
    Tl : float
        Touschek lifetime in seconds
    """

def touschek_lifetime_ele_by_ele(mode: _pybmad.NormalModesStruct, lat: _pybmad.LatStruct, momentum_aperture: _pybmad.MomentumApertureStructArray1D) -> float:
    """
    Calculates the touschek lifetime for a lattice by calling touschek_rate1
    for each element the momentum_aperture array of momentum_aperture_structs.
    This calculation is based on Piwinski 1998 "The Touschek Effect In
    Strong Focusing Storage Rings".  This is the most general case, equation 31.
    42.

    A common way to call this function is to first populate mode using
    radiation integrals.  If an ideal lattice is used, the vertical
    emittance must also be set to a reasonable value.  If the vertical
    emittance is due only to quantum excitation, then it will likely be
    several orders of magnitude smaller than any real physical situation, in which
    case the integral in this function will have problems converging.

    In addition to setting mode, also set lat%param%n_part to the number of particles
    per bunch.

    This function assumes that the twiss parameters
    been calculated, and that mode has been populated with emittance and bunch length.

    Parameters
    ----------
    mode : NormalModesStruct
        beam properties

    lat : LatStruct
        Lattice

    momentum_aperture : 1D array of MomentumApertureStruct
        ele-by-ele unsymmatric apertures

    Returns
    -------
    Tl : float
        Touschek lifetime in seconds
    """

def touschek_lifetime_with_aperture(mode: _pybmad.NormalModesStruct, lat: _pybmad.LatStruct, momentum_aperture: _pybmad.MomentumApertureStructArray1D) -> float:
    """
    Calculates the touschek lifetime for a lattice by calling touschek_rate1
    for each s-coordinate in the momentum_aperture array of momentum_aperture_structs.
    This calculation is based on Piwinski 1998 "The Touschek Effect In
    Strong Focusing Storage Rings".  This is the most general case, equation 31.
    42.

    A common way to call this function is to first populate mode using
    radiation integrals.  If an ideal lattice is used, the vertical
    emittance must also be set to a reasonable value.  If the vertical
    emittance is due only to quantum excitation, then it will likely be
    several orders of magnitude smaller than any real physical situation, in which
    case the integral in this function will have problems converging.

    In addition to setting mode, also set lat%param%n_part to the number of particles
    per bunch.

    This function assumes that the twiss parameters
    been calculated, and that mode has been populated with emittance and bunch length.

    This function assumes that momentum_aperture(0)%s==0 and momentum_aperture(last)%s==lat%param%total_length.

    Parameters
    ----------
    mode : NormalModesStruct
        beam properties

    lat : LatStruct
        Lattice

    momentum_aperture : 1D array of MomentumApertureStruct
        loc-by-loc unsymmatric apertures

    Returns
    -------
    Tl : float
        Touschek lifetime in seconds
    """

def touschek_rate1(mode: _pybmad.NormalModesStruct, lat: _pybmad.LatStruct, ix: int | None = None, s: float | None = None) -> float:
    """
    Calculates the touschek rate at the location specified by s or ix
    This calculation is based on Piwinski 1998 "The Touschek Effect In
    Strong Focusing Storage Rings".  This is the most general case, equation
    31.

    This function uses twiss_and_track_at_s to determine the Twiss parameters
    at the location s or element index ix.

    A common way to call this function is to first populate mode using
    radiation integrals.  If an ideal lattice is used, the vertical
    emittance must also be set to a reasonable value.  If the vertical
    emittance is due only to quantum excitation, then it will likely be
    several orders of magnitude smaller than any real physical situation, in which
    case the integral in this function will have problems converging.
    Additionally, mode%pz_aperture needs to be set to the momentum aperture.

    In addition to setting mode, also set lat%param%n_part to the number of particles
    per bunch.

    IMPORTANT NOTE: If the lattice type is a circular lattice, then
                    mode%a%emittance and mode%b%emittance are assumed to
                    contain the normalized emittences.  If lattice geometry is
                    open, the emittances are assumed to be
                    unnormalized.

    IMPORTANT NOTE: The output of this subroutine is the loss rate assuming
                    that two particles are lost per collision, one with too
                    much energy, and one with too little energy.  This agrees
                    with Piwinski's original derivation, which assumes that the
                    positive energy aperture is equal in magnitude to the
                    negative energy aperture.  If you are studying an
                    accelerator with a non-symmetric energy aperture, then
                    this subroutine should be called twice, once with the positive
                    aperture, and once with the negative aperture, and rate from
                    each call should be halved and summed.

    Parameters
    ----------
    mode : NormalModesStruct
        beam properties

    lat : LatStruct
        Lattice

    ix : int, optional
        element index (either s or ix must be specified)

    s : float, optional
        location in meters (either s or ix must be specified)

    Returns
    -------
    rate : float
        Touschek rate, in units particle per second, assuming two particles per event.
    """

def touschek_rate1_zap(mode: _pybmad.NormalModesStruct, rate: float, lat: _pybmad.LatStruct, ix: int | None = None, s: float | None = None) -> None:
    """
    Wrapper for Fortran routine touschek_rate1_zap

    Parameters
    ----------
    mode : NormalModesStruct

    rate : float

    lat : LatStruct

    ix : int, optional

    s : float, optional
    """

class Track1:
    """track1 return type"""

    @property
    def end_orb(self) -> _pybmad.CoordStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track1(start_orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, track: _pybmad.TrackStruct | None = None, ignore_radiation: bool | None = None, make_map1: bool | None = None, init_to_edge: bool | None = None) -> Track1:
    """
    Wrapper for Fortran routine track1

    Parameters
    ----------
    start_orb : CoordStruct
        Starting position.

    ele : EleStruct
        Element to track through.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Modified if make_map1 is True.

    param : LatParamStruct
        Reference particle info.

    track : TrackStruct, optional
        Structure holding existing track.
        This parameter is an input/output and is modified in-place.
        As an output, track: Structure holding the track information if the

    ignore_radiation : bool, optional
        If present and True then do not include radiation effects along with space charge effects.

    make_map1 : bool, optional
        Make ele.mat6 and ele.spin_q components? Default is false.

    init_to_edge : bool, optional
        Default is True. If True then force the tracked particle to begin at the element's edge. See above. Do not
        use this argument unless you know what you are doing.

    Returns
    -------
    end_orb : CoordStruct
        End position.

    err_flag : bool, optional
        Set true if there is an error. False otherwise. Note: The particle getting lost (EG hitting an aperture)
        is *not* an error. An error is something like start_orb not being properly initialized.
    """

def track1_beam(beam: _pybmad.BeamStruct, ele: _pybmad.EleStruct, centroid: _pybmad.CoordStructArray1D | None = None, direction: int | None = None) -> bool:
    """
    Subroutine to track a beam of particles through an element.

    Parameters
    ----------
    beam : BeamStruct
        Starting beam position.
        This parameter is an input/output and is modified in-place.
        As an output, beam: Ending beam position.

    ele : EleStruct
        element to track through.

    centroid : 1D array of CoordStruct, optional
        Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before beam tracking by
        tracking a single particle.

    direction : int, optional
        +1 (default) -> Track forward, -1 -> Track backwards.

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost for a CSR calc.
    """

class Track1Bmad:
    """track1_bmad return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track1_bmad(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> Track1Bmad:
    """
    Wrapper for Fortran routine track1_bmad

    Parameters
    ----------
    orbit : CoordStruct
        Starting position
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Element

    param : LatParamStruct

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix propagated through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    err_flag : bool, optional
        Set true if there is an error. False otherwise.

    track : TrackStruct, optional
        Structure holding the track information if the lattice element does tracking step-by-step. See track1 for
        more details.
    """

def track1_bmad_photon(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> bool:
    """
    Wrapper for Fortran routine track1_bmad_photon

    Parameters
    ----------
    orbit : CoordStruct
        Starting position
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position

    ele : EleStruct
        Element

    param : LatParamStruct

    Returns
    -------
    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

def track1_bunch(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, centroid: _pybmad.CoordStructArray1D | None = None, direction: int | None = None, bunch_track: _pybmad.BunchTrackStruct | None = None) -> bool:
    """
    Subroutine to track a bunch of particles through an element.

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position.

    ele : EleStruct
        element to track through.

    centroid : 1D array of CoordStruct, optional
        Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before beam tracking by
        tracking a single particle.

    direction : int, optional
        +1 (default) -> Track forward, -1 -> Track backwards.

    bunch_track : BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: Track information appended to track.

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost for a CSR calc.
    """

def track1_bunch_csr(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, centroid: _pybmad.CoordStructArray1D, s_start: float | None = None, s_end: float | None = None, bunch_track: _pybmad.BunchTrackStruct | None = None) -> bool:
    """
    $OMP THREADPRIVATE(csr_kick1_ptr, csr_csr_ptr, csr_einfo_s_ptr, csr_dr_match_ptr)

     Subroutine track1_bunch_csr (bunch, ele, centroid, err, s_start, s_end, bunch_track)

     Routine to track a bunch of particles through an element with csr radiation effects.

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position.

    ele : EleStruct
        The element to track through. Must be part of a lattice.

    centroid : 1D array of CoordStruct
        coord_struct, Approximate beam centroid orbit for the lattice branch. Calculate this before beam tracking
        by tracking a single particle.

    s_start : float, optional
        Starting position relative to ele. Default = 0

    s_end : float, optional
        Ending position. Default is ele length.

    bunch_track : BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: track information if the tracking method does

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost.
    """

def track1_bunch_csr3d(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, centroid: _pybmad.CoordStructArray1D, s_start: float | None = None, s_end: float | None = None, bunch_track: _pybmad.BunchTrackStruct | None = None) -> bool:
    """
    EXPERIMENTAL. NOT CURRENTLY OPERATIONAL!

    Routine to track a bunch of particles through an element using
    steady-state 3D CSR.

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position.

    ele : EleStruct
        The element to track through. Must be part of a lattice.

    centroid : 1D array of CoordStruct
        coord_struct, Approximate beam centroid orbit for the lattice branch. Calculate this before beam tracking
        by tracking a single particle.

    s_start : float, optional
        Starting position relative to ele. Default = 0

    s_end : float, optional
        Ending position. Default is ele length.

    bunch_track : BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: track information if the tracking method does

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost.

    Notes
    -----
    The core routines are from the OpenCSR package developed at:
    https://github.com/ChristopherMayes/OpenCSR
    """

def track1_bunch_hom(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, direction: int | None = None, bunch_track: _pybmad.BunchTrackStruct | None = None) -> None:
    """
    Subroutine to track a bunch of particles through an element including wakefields.

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position.

    ele : EleStruct
        The element to track through.

    direction : int, optional
        +1 (default) -> Track forward, -1 -> Track backwards.

    bunch_track : BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: Track information appended to track.
    """

def track1_bunch_space_charge(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, track_to_same_s: bool | None = None, bunch_track: _pybmad.BunchTrackStruct | None = None) -> bool:
    """
    Wrapper for Fortran routine track1_bunch_space_charge

    Parameters
    ----------
    bunch : BunchStruct
        Starting bunch position.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Ending bunch position.

    ele : EleStruct
        Element to track through. Must be part of a lattice.

    track_to_same_s : bool, optional
        Default is True. If True, drift particles to all have the same s-position.

    bunch_track : BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: track information if the tracking method does

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost for a CSR calc.
    """

def track1_crystal(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track diffraction from a crystal.

    Parameters
    ----------
    ele : EleStruct
        Element tracking through.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_diffraction_plate_or_mask(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track through diffraction plate and mask elements.

    Parameters
    ----------
    ele : EleStruct
        Diffraction plate or mask element.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_high_energy_space_charge(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to apply the ultra-relative space charge kick to a particle at the end of an element.
    The routine setup_high_energy_space_charge_calc must be called initially before any tracking is done.
    This routine assumes a Gaussian bunch and is only valid with relativistic particles where the
    effect of the space charge is small.

    Parameters
    ----------
    ele : EleStruct
        Element tracked through.

    param : LatParamStruct

    orbit : CoordStruct
        Starting position
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position
    """

def track1_lens(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track through a lens.

    Parameters
    ----------
    ele : EleStruct
        Element tracking through.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_linear(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> None:
    """
    Wrapper for Fortran routine track1_linear

    Parameters
    ----------
    orbit : CoordStruct
        Starting position
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position

    ele : EleStruct
        Element

    param : LatParamStruct
    """

def track1_lr_wake(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct) -> None:
    """
    Subroutine to put in the long-range wakes for particle tracking.

    Parameters
    ----------
    bunch : BunchStruct
        Bunch to track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Kicked bunch.

    ele : EleStruct
        Element with wakes.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with updated wake amplitudes.
    """

def track1_mad(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> None:
    """
    Subroutine to track through an element using a 2nd order transfer map.
    Note: If map does not exist then one will be created.

    Parameters
    ----------
    orbit : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coords.

    ele : EleStruct
        Element to track through.

    param : LatParamStruct
        Lattice parameters.
    """

def track1_mirror(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track reflection from a mirror.

    Parameters
    ----------
    ele : EleStruct
        Element tracking through.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_mosaic_crystal(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track diffraction from a crystal.

    Parameters
    ----------
    ele : EleStruct
        Element tracking through.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_multilayer_mirror(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track reflection from a multilayer_mirror.
    Basic equations are from Kohn, "On the Theory of Reflectivity of an X-Ray Multilayer Mirror".

    Parameters
    ----------
    ele : EleStruct
        Element tracking through.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_radiation(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, edge: int) -> None:
    """
    Subroutine to apply a kick to a particle to account for radiation dampling and/or fluctuations.

    For tracking through a given element, this routine should be called initially when
    the particle is at the entrance end and at the end when the particle is at the exit end, when
    the orbit is with respect to laboratory (not element body) coordinates.
    That is, each time this routine is called it applies half the radiation kick for the entire element.

    Note: This routine is called by track1.

    Parameters
    ----------
    orbit : CoordStruct
        Particle position before radiation applied.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Particle position after radiation has been applied.

    ele : EleStruct
        Element generating radiation.

    edge : int
        Where the particle is: start_edge$ or end_edge$.
    """

def track1_radiation_center(orbit: _pybmad.CoordStruct, ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct, rad_damp: bool | None = None, rad_fluct: bool | None = None) -> None:
    """
    Used for elements that have been split in half: This routine applies a kick to a particle
    to account for radiation dampling and/or fluctuations.

    Also see: track1_radiation.

    Parameters
    ----------
    orbit : CoordStruct
        Particle at center of element before radiation applied.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Particle position after radiation has been applied.

    ele1 : EleStruct
        First half of the split element.

    ele2 : EleStruct
        Second half of the split element.

    rad_damp : bool, optional
        If present, override setting of bmad_com.radiation_damping_on.

    rad_fluct : bool, optional
        If present, override setting of bmad_com.radiation_fluctuations_on.
    """

class Track1RungeKutta:
    """track1_runge_kutta return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track1_runge_kutta(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> Track1RungeKutta:
    """
    Wrapper for Fortran routine track1_runge_kutta

    Parameters
    ----------
    orbit : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coords.

    ele : EleStruct
        Ele_struct

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix propagated through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.

    track : TrackStruct, optional
        Structure holding the track information.
    """

def track1_sample(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, orbit: _pybmad.CoordStruct) -> None:
    """
    Routine to track reflection from a sample element.

    Parameters
    ----------
    ele : EleStruct
        Element tracking through.

    param : LatParamStruct
        lattice parameters.

    orbit : CoordStruct
        phase-space coords to be transformed
        This parameter is an input/output and is modified in-place.
        As an output, orbit: final phase-space coords
    """

def track1_spin(start_orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, end_orb: _pybmad.CoordStruct, make_quaternion: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track1_spin

    Parameters
    ----------
    start_orb : CoordStruct
        Starting coords.

    ele : EleStruct
        Element to track through.

    param : LatParamStruct
        Beam parameters.

    end_orb : CoordStruct
        Ending coords.

    make_quaternion : bool, optional
        If present and true then calculate the 1st order spin map which is represented as a quaternion.
    """

def track1_spin_integration(start_orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, end_orb: _pybmad.CoordStruct) -> None:
    """
    Wrapper for Fortran routine track1_spin_integration

    Parameters
    ----------
    start_orb : CoordStruct
        Starting coords.

    ele : EleStruct
        Element to track through.

    param : LatParamStruct
        Beam parameters.

    end_orb : CoordStruct
        Ending coords.
    """

def track1_spin_taylor(start_orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> _pybmad.CoordStruct:
    """
    Wrapper for Fortran routine track1_spin_taylor

    Parameters
    ----------
    start_orb : CoordStruct
        Starting coords.

    ele : EleStruct
        Element to track through.

    param : LatParamStruct
        Beam parameters.

    Returns
    -------
    end_orb : CoordStruct
    """

def track1_sr_wake(bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct) -> None:
    """
    Subroutine to apply the short range wake fields to a bunch.

    Parameters
    ----------
    bunch : BunchStruct
        Bunch of particles.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Bunch with wakefields applied to the particles.

    ele : EleStruct
        Element with wakefields.
    """

def track1_symp_lie_ptc(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> _pybmad.TrackStruct:
    """
    Wrapper for Fortran routine track1_symp_lie_ptc

    Parameters
    ----------
    orbit : CoordStruct
        Starting position
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position

    ele : EleStruct
        Element

    param : LatParamStruct

    Returns
    -------
    track : TrackStruct, optional
        Structure holding the track information.
    """

def track1_taylor(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, taylor: _pybmad.TaylorStructArray1D | None = None, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track1_taylor

    Parameters
    ----------
    orbit : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coords.

    ele : EleStruct
        Element to track through.

    taylor : 1D array of TaylorStruct (shape: 6), optional
        Alternative map to use instead of ele.taylor.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

class Track1TimeRungeKutta:
    """track1_time_runge_kutta return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    @property
    def dt_step(self) -> float | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track1_time_runge_kutta(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, t_end: float | None = None, dt_step: float | None = None) -> Track1TimeRungeKutta:
    """
    Wrapper for Fortran routine track1_time_runge_kutta

    Parameters
    ----------
    orbit : CoordStruct
        starting position, z-based coords
        This parameter is an input/output and is modified in-place.
        As an output, orbit: end position, z-based coords

    ele : EleStruct
        element

    param : LatParamStruct
        lattice parameters

    t_end : float, optional
        If present, maximum time to which the particle will be tracked. Used for tracking with given time steps.
        The time orb.t at which tracking stops may be less than this if the particle gets to the end of the
        element

    dt_step : float, optional
        If positive, next RK time step to take. This overrides bmad_com.init_ds_adaptive_tracking. Used by
        track_bunch_time.
        This parameter is an input/output and is modified in-place.
        As an output, dt_step: Next RK time step that this tracker would take based on the error tolerance.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise

    track : TrackStruct, optional
        Contains array of the step-by-step particle trajectory along with the field at these positions. When
        tracking through multiple elements, the trajectory in an element is appended to the existing trajectory.
        To reset: Set track.n_pt = -1.

    dt_step : float, optional
        If positive, next RK time step to take. This overrides bmad_com.init_ds_adaptive_tracking. Used by
        track_bunch_time.
        This parameter is an input/output and is modified in-place.
        As an output, dt_step: Next RK time step that this tracker would take based on the error tolerance.
    """

class TrackABeambeam:
    """track_a_beambeam return type"""

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    @property
    def mat6(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_a_beambeam(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> TrackABeambeam:
    """
    Wrapper for Fortran routine track_a_beambeam

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Beambeam element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    track : TrackStruct, optional
        Structure holding the track information if the lattice element does tracking step-by-step. See track1 for
        more details.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_bend(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track_a_bend

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Bend element.

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix to the element end.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def track_a_bend_photon(orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct, length: float) -> None:
    """
    Routine to track a photon through a dipole bend.
    The photon is traveling in a straight line but the reference frame
    is curved in a circular shape.

    Parameters
    ----------
    orb : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orb: End position.

    ele : EleStruct
        Bend element.

    length : float
        length to track.
    """

def track_a_capillary(orb: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> None:
    """
    $OMP THREADPRIVATE(cap_ele_ptr, cap_photon_ptr, cap_photon1_ptr)

     Subroutine track_a_capillary (orb, ele)

     Routine to track through a capillary.

    Parameters
    ----------
    orb : CoordStruct
        Input photon coordinates.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Output photon coordinates.

    ele : EleStruct
        Capillary element
    """

def track_a_converter(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_converter

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        converter element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_crab_cavity(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_crab_cavity

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        crab_cavity element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

class TrackADrift:
    """track_a_drift return type"""

    @property
    def time(self) -> float | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_a_drift(orb: _pybmad.CoordStruct, length: float, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None, ele_orientation: int | None = None, include_ref_motion: bool | None = None, time: float | None = None) -> TrackADrift:
    """
    Wrapper for Fortran routine track_a_drift

    Parameters
    ----------
    orb : CoordStruct
        Orbit at start of the drift.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Orbit at end of the drift.

    length : float
        Length to drift through in body coordinates. --    If orb.direction = 1, positive length is in +z
        direction and vice versa.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the drift.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix including the drift.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    ele_orientation : int, optional
        Element orientation. Default is orb.direction.

    include_ref_motion : bool, optional
        Include effect of the motion of the reference particle? Default is True. False is basically only used by
        offset_particle. Additionally, if False, orb.s is not changed.

    time : float, optional
        Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
        This parameter is an input/output and is modified in-place.
        As an output, time: Updated time.

    Returns
    -------
    time : float, optional
        Particle time before drifting. Typically this is an RF clock time which may not be equal to orb.t
        This parameter is an input/output and is modified in-place.
        As an output, time: Updated time.
    """

def track_a_drift_photon(orb: _pybmad.CoordStruct, length: float, phase_relative_to_ref: bool) -> None:
    """
    Wrapper for Fortran routine track_a_drift_photon

    Parameters
    ----------
    orb : CoordStruct
        Orbit at start of the drift.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Orbit at end of the drift

    length : float
        Longitudinal length to drift through.

    phase_relative_to_ref : bool
        If true then E field phase shift is relative to ref particle.
    """

def track_a_foil(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_foil

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        foil element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is False.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_gkicker(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track_a_gkicker

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Gkicker

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def track_a_lcavity(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track_a_lcavity

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Thick multipole element.

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def track_a_lcavity_old(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track_a_lcavity_old

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Thick multipole element.

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def track_a_mask(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_mask

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Mask element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_match(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, err_flag: bool | None = None, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_match

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Match element.

    param : LatParamStruct
        Lattice parameters.

    err_flag : bool, optional

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

class TrackAPatch:
    """track_a_patch return type"""

    @property
    def s_ent(self) -> float: ...

    @property
    def ds_ref(self) -> float: ...

    @property
    def mat6(self) -> list[list[float]] | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_a_patch(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, drift_to_exit: bool | None = None, track_spin: bool | None = None, make_matrix: bool | None = None) -> TrackAPatch:
    """
    Wrapper for Fortran routine track_a_patch

    Parameters
    ----------
    ele : EleStruct
        patch element.

    orbit : CoordStruct
        Starting phase space coords
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Coords after applying a patch transformation.

    drift_to_exit : bool, optional
        If False then do not drift the particle from beginning to end face. Also do not correct for a reference
        energy shift. Default is True.

    track_spin : bool, optional
        If True rotate the spin vector appropriately. If ele.spin_tracking_method = symp_lie_ptc -> default =
        True. Else default = False.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    s_ent : float, optional
        Longitudinal coordinate of the initial particle position in the frame of reference of the face where the
        particle exits. For a patch with positive z_offset and all other attributes zero, s_ent = -z_offset.

    ds_ref : float, optional
        Distance reference particle travels from entrance to exit.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_patch_photon(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, drift_to_exit: bool | None = None, use_z_pos: bool | None = None) -> None:
    """
    Routine to track through a patch element with a photon.
    The steps for tracking are:
      1) Transform from entrance to exit coordinates.
      2) Drift particle from the entrance to the exit coordinants.

    Parameters
    ----------
    ele : EleStruct
        patch element.

    orbit : CoordStruct
        Starting phase space coords
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Coords after applying a patch transformation.

    drift_to_exit : bool, optional
        If False then do not drift the particle from start to ending faces. Default is True.

    use_z_pos : bool, optional
        If present and True, use orbit.vec(5) as the true z-position relative to the start of the element instead
        of assuming that the particle is at the patch edge.
    """

def track_a_pickup(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, err_flag: bool | None = None, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_pickup

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Pickup element.

    param : LatParamStruct
        Lattice parameters.

    err_flag : bool, optional

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_quadrupole(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_quadrupole

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Quadrupole element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_rfcavity(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_rfcavity

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        rfcavity element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_sad_mult(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track_a_sad_mult

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Sad_mult element.

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix up to the sad_mult.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def track_a_sol_quad(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_sol_quad

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Sol_quad or solenoid element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

def track_a_thick_multipole(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, mat6: Sequence[Sequence[float]] | None = None, make_matrix: bool | None = None) -> None:
    """
    Wrapper for Fortran routine track_a_thick_multipole

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Thick multipole element.

    param : LatParamStruct
        Lattice parameters.

    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix before the element.
        This parameter is an input/output and is modified in-place.
        As an output, mat6: Transfer matrix through the element.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.
    """

def track_a_wiggler(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, make_matrix: bool | None = None) -> list[list[float]] | None:
    """
    Wrapper for Fortran routine track_a_wiggler

    Parameters
    ----------
    orbit : CoordStruct
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: End position.

    ele : EleStruct
        Wiggler element.

    param : LatParamStruct
        Lattice parameters.

    make_matrix : bool, optional
        Propagate the transfer matrix? Default is false.

    Returns
    -------
    mat6 : 2D array of float (shape: 6,6), optional
        Transfer matrix through the element.
    """

class TrackAZeroLengthElement:
    """track_a_zero_length_element return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_a_zero_length_element(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> TrackAZeroLengthElement:
    """
    Wrapper for Fortran routine track_a_zero_length_element

    Parameters
    ----------
    orbit : CoordStruct
        Starting coords.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Ending coords.

    ele : EleStruct
        Element tracked through.

    param : LatParamStruct
        Lattice parameters.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.

    track : TrackStruct, optional
        Structure holding the track information.
    """

class TrackAll:
    """track_all return type"""

    @property
    def track_state(self) -> int: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def orbit0(self) -> _pybmad.CoordStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_all(lat: _pybmad.LatStruct, orbit: _pybmad.CoordStructAlloc1D, ix_branch: int | None = None, init_lost: bool | None = None) -> TrackAll:
    """
    Wrapper for Fortran routine track_all

    Parameters
    ----------
    lat : LatStruct
        Lat to track through.

    orbit : 1D array of CoordStruct
        orbit(0) is the starting coordinates for tracking. If not allocated, the zero orbit will be used.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit array.

    ix_branch : int, optional
        Index of branch to track. Default is 0 (main branch).

    init_lost : bool, optional
        Default if False. If True, initialize orbit(N) terms that are not tracked through due to particle loss.

    Returns
    -------
    track_state : int, optional
        Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.

    err_flag : bool, optional
        Set true if particle lost or error. False otherwise

    orbit0 : 1D array of CoordStruct, optional
        Orbit array for branch 0. Used to fill in the orbit at lord elemenets. Only needed when orbit(:) is not
        the orbit for branch 0.
    """

def track_beam(lat: _pybmad.LatStruct, beam: _pybmad.BeamStruct, ele1: _pybmad.EleStruct | None = None, ele2: _pybmad.EleStruct | None = None, centroid: _pybmad.CoordStructArray1D | None = None, direction: int | None = None, bunch_tracks: _pybmad.BunchTrackStructArray1D | None = None) -> bool:
    """
    Subroutine to track a beam of particles from the end of
    ele1 Through to the end of ele2. Both must be in the same lattice branch.

    Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.

    Parameters
    ----------
    lat : LatStruct
        Lattice to track through.

    beam : BeamStruct
        Beam at end of element ix1.
        This parameter is an input/output and is modified in-place.
        As an output, beam: Beam at end of element ix2.

    ele1 : EleStruct, optional
        Starting element (this element is NOT tracked through). Default is lat.ele(0).

    ele2 : EleStruct, optional
        Ending element. Default is lat.ele(lat.n_ele_track).

    centroid : 1D array of CoordStruct, optional
        Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before beam tracking by
        tracking a single particle.

    direction : int, optional
        +1 (default) -> Track forward, -1 -> Track backwards.

    bunch_tracks : 1D array of BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_tracks: track information if the tracking method does

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost for a CSR calc.
    """

def track_bunch(lat: _pybmad.LatStruct, bunch: _pybmad.BunchStruct, ele1: _pybmad.EleStruct | None = None, ele2: _pybmad.EleStruct | None = None, centroid: _pybmad.CoordStructArray1D | None = None, direction: int | None = None, bunch_track: _pybmad.BunchTrackStruct | None = None) -> bool:
    """
    Subroutine to track a particle bunch from the end of ele1 Through to the end of ele2.
    Both must be in the same lattice branch.
    With forward tracking, if ele2 is at or before ele1, the tracking will "wrap" around
    the ends of the lattice.

    Note: To zero wakes between runs, zero_lr_wakes_in_lat needs to be called.

    Parameters
    ----------
    lat : LatStruct
        Lattice to track through.

    bunch : BunchStruct
        Bunch at end of element ix1.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Bunch at end of element ix2.

    ele1 : EleStruct, optional
        Starting element (this element is NOT tracked through). Default is lat.ele(0).

    ele2 : EleStruct, optional
        Ending element. Default is lat.ele(lat.n_ele_track).

    centroid : 1D array of CoordStruct, optional
        Approximate centroid orbit. Only needed if CSR is on. Hint: Calculate this before bunch tracking by
        tracking a single particle.

    direction : int, optional
        +1 (default) -> Track forward, -1 -> Track backwards.

    bunch_track : BunchTrackStruct, optional
        Existing tracks. If bunch_track.n_pt = -1 then Overwrite any existing track.
        This parameter is an input/output and is modified in-place.
        As an output, bunch_track: track information if the tracking method does

    Returns
    -------
    err : bool
        Set true if there is an error. EG: Too many particles lost for a CSR calc.
    """

def track_bunch_time(bunch: _pybmad.BunchStruct, branch: _pybmad.BranchStruct, t_end: float, s_end: float, dt_step: _pybmad.RealArray1D | None = None, extra_field: _pybmad.EmFieldStructArray1D | None = None) -> None:
    """
    Wrapper for Fortran routine track_bunch_time

    Parameters
    ----------
    bunch : BunchStruct
        Coordinates must be time-coords in element body frame.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Coordinates will be time-coords in element body frame.

    branch : BranchStruct
        Lattice branch being tracked through.

    t_end : float
        Ending time.

    s_end : float
        Ending s-position.

    dt_step : 1D array of float, optional
        Initial step to take for each particle. Overrides bmad_com.init_ds_adaptive_tracking.
        This parameter is an input/output and is modified in-place.
        As an output, dt_step: Next RK time step that this tracker would take based on the error tolerance.

    extra_field : 1D array of EmFieldStruct, optional
        Per particle static field to be added to the lattice element field. Eg used with space charge.
    """

def track_bunch_to_s(bunch: _pybmad.BunchStruct, s: float, branch: _pybmad.BranchStruct) -> None:
    """
    Drift a bunch of particles to the same s coordinate

    Parameters
    ----------
    bunch : BunchStruct
        Input bunch position in s-based coordinate.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Output bunch position in s-based coordinate. Particles will be at the same s
        coordinate

    s : float
        Target coordinate.

    branch : BranchStruct
        Branch being tracked through.
    """

def track_bunch_to_t(bunch: _pybmad.BunchStruct, t_target: float, branch: _pybmad.BranchStruct) -> None:
    """
    Drift a bunch of particles to the same t coordinate

    Parameters
    ----------
    bunch : BunchStruct
        Input bunch position in s-based coordinate.
        This parameter is an input/output and is modified in-place.
        As an output, bunch: Output bunch position in s-based coordinate. Particles will be at the same t
        coordinate

    t_target : float
        Target t coordinate.

    branch : BranchStruct
        Lattice branch being tracked through.
    """

def track_complex_taylor(start_orb: _pybmad.ComplexArray1D, complex_taylor: _pybmad.ComplexTaylorStructArray1D, end_orb: _pybmad.ComplexArray1D) -> None:
    """
    Subroutine to track using a complex_taylor map.

    Parameters
    ----------
    start_orb : 1D array of complex
        Starting coords.

    complex_taylor : 1D array of ComplexTaylorStruct
        complex_taylor map.

    end_orb : 1D array of complex
        Ending coords.
    """

class TrackFromSToS:
    """track_from_s_to_s return type"""

    @property
    def orbit_end(self) -> _pybmad.CoordStruct: ...

    @property
    def all_orb(self) -> _pybmad.CoordStructAlloc1D: ...

    @property
    def track_state(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_from_s_to_s(lat: _pybmad.LatStruct, s_start: float, s_end: float, orbit_start: _pybmad.CoordStruct, ix_branch: int | None = None, ix_ele_end: int | None = None) -> TrackFromSToS:
    """
    Wrapper for Fortran routine track_from_s_to_s

    Parameters
    ----------
    lat : LatStruct
        Lattice to track through

    s_start : float
        Starting s-position.

    s_end : float
        Ending s-position. If <= s_start then will wrap

    orbit_start : CoordStruct
        Starting coordinates.

    ix_branch : int, optional
        Lattice branch index. Default is 0 (main branch).

    ix_ele_end : int, optional
        If present, ignore s_end and track to in between ix_ele_end and ix_ele_end+1

    Returns
    -------
    orbit_end : CoordStruct
        Ending coordinates.

    all_orb : 1D array of CoordStruct, optional
        If present then the orbit at the exit ends of the elements tracked through will be recorded in this
        structure.

    track_state : int, optional
        Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.
    """

def track_func(s_target: float, status: int) -> float:
    """
    Wrapper for Fortran routine track_func

    Parameters
    ----------
    s_target : float

    status : int

    Returns
    -------
    dt : float
    """

def track_many(lat: _pybmad.LatStruct, orbit: _pybmad.CoordStructArray1D, ix_start: int, ix_end: int, direction: int, ix_branch: int | None = None) -> int:
    """
    Wrapper for Fortran routine track_many

    Parameters
    ----------
    lat : LatStruct
        Lat to track through.

    orbit : 1D array of CoordStruct
        Coordinates at start of tracking.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Orbit.

    ix_start : int
        Start index (See Note).

    ix_end : int
        End index (See Note).

    direction : int
        Direction to track. = +1 -> Track forward (+s) = -1 -> Track backward (-s)

    ix_branch : int, optional
        Branch to track. Default is 0 (main lattice).

    Returns
    -------
    track_state : int, optional
        Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.
    """

def track_to_surface(ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, param: _pybmad.LatParamStruct) -> list[list[float]]:
    """
    Wrapper for Fortran routine track_to_surface

    Parameters
    ----------
    ele : EleStruct
        Element

    orbit : CoordStruct
        Coordinates in the element coordinate frame
        This parameter is an input/output and is modified in-place.
        As an output, orbit: At surface in local surface coordinate frame

    param : LatParamStruct
        Branch parameters.

    Returns
    -------
    w_surface : 2D array of float (shape: 3,3)
        real(rp), rotation matrix to transform to surface coords.
    """

class TrackUntilDead:
    """track_until_dead return type"""

    @property
    def end_orb(self) -> _pybmad.CoordStruct: ...

    @property
    def track(self) -> _pybmad.TrackStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def track_until_dead(start_orb: _pybmad.CoordStruct, lat: _pybmad.LatStruct) -> TrackUntilDead:
    """
    Subroutine to track a particle arbitrarily through a lattice, forwards or backwards,
      until it is lost or exits the lattice.

      The starting element is located using start_orb%s.

    Parameters
    ----------
    start_orb : CoordStruct
        Starting coords.

    lat : LatStruct
        lattice that contains and element at start_orb.s

    Returns
    -------
    end_orb : CoordStruct
        final coords

    track : TrackStruct, optional
        (optional)
    """

class TrackingRadMapSetup:
    """tracking_rad_map_setup return type"""

    @property
    def rad_map(self) -> _pybmad.RadMapStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tracking_rad_map_setup(ele: _pybmad.EleStruct, tollerance: float, ref_edge: int) -> TrackingRadMapSetup:
    """
    Wrapper for Fortran routine tracking_rad_map_setup

    Parameters
    ----------
    ele : EleStruct
        Element to setup. Matrices will be with respect to the map reference orbit.

    tollerance : float
        Tolerance used for the computation.

    ref_edge : int
        Edge that the matrices are referenced to. upstream_end$ or downstream_end$.

    Returns
    -------
    rad_map : RadMapStruct
        Structure holding the matrices.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def transfer_ac_kick(ac_in: _pybmad.AcKickerStruct) -> _pybmad.AcKickerStruct | None:
    """
    Wrapper for Fortran routine transfer_ac_kick

    Parameters
    ----------
    ac_in : AcKickerStruct
        Input

    Returns
    -------
    ac_out : AcKickerStruct, optional
        Gets set equal to ac_in
    """

def transfer_branch(branch1: _pybmad.BranchStruct) -> _pybmad.BranchStruct:
    """
    Wrapper for Fortran routine transfer_branch

    Parameters
    ----------
    branch1 : BranchStruct

    Returns
    -------
    branch2 : BranchStruct
    """

def transfer_branch_parameters(branch_in: _pybmad.BranchStruct) -> _pybmad.BranchStruct:
    """
    Wrapper for Fortran routine transfer_branch_parameters

    Parameters
    ----------
    branch_in : BranchStruct
        Input branch.

    Returns
    -------
    branch_out : BranchStruct
        Output branch with parameters set.
    """

def transfer_branches(branch1: _pybmad.BranchStructArray1D, branch2: _pybmad.BranchStructArray1D) -> None:
    """
    Wrapper for Fortran routine transfer_branches

    Parameters
    ----------
    branch1 : 1D array of BranchStruct

    branch2 : 1D array of BranchStruct
    """

def transfer_ele(ele1: _pybmad.EleStruct, nullify_pointers: bool | None = None) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine transfer_ele

    Parameters
    ----------
    ele1 : EleStruct

    nullify_pointers : bool, optional
        If present and True then nullify the pointers in ele2 except for the ele2.lat and ele2.lord pointers. This
        gives a "bare bones" copy where one does not have to worry about deallocating allocated structure
        components later.

    Returns
    -------
    ele2 : EleStruct
    """

def transfer_ele_taylor(ele_in: _pybmad.EleStruct, taylor_order: int | None = None) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine transfer_ele_taylor

    Parameters
    ----------
    ele_in : EleStruct
        Element with the Taylor map.

    taylor_order : int, optional
        Order to truncate the Taylor map at.

    Returns
    -------
    ele_out : EleStruct
        Element receiving the Taylor map truncated to order taylor_order.
    """

def transfer_eles(ele1: _pybmad.EleStructArray1D, ele2: _pybmad.EleStructArray1D) -> None:
    """
    Wrapper for Fortran routine transfer_eles

    Parameters
    ----------
    ele1 : 1D array of EleStruct

    ele2 : 1D array of EleStruct
    """

def transfer_fieldmap(ele_in: _pybmad.EleStruct, who: int) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine transfer_fieldmap

    Parameters
    ----------
    ele_in : EleStruct
        Input element.

    who : int
        Possibilities are: all$, cartesian_map$, cylindrical_map$, or grid_field$

    Returns
    -------
    ele_out : EleStruct
        Output element.
    """

def transfer_fixer_params(fixer: _pybmad.EleStruct, to_stored: bool, orbit: _pybmad.CoordStruct | None = None, who: str | None = None) -> bool:
    """
    Set parameters of fixer.

    Parameters
    ----------
    fixer : EleStruct
        Fixer element to set.

    to_stored : bool
        If False, set real Twiss from stored. If True, set stored Twiss from real.

    orbit : CoordStruct, optional
        Used for 'phase_space' transfers. Used for input if to_stored = True.
        This parameter is an input/output and is modified in-place.
        As an output, orbit: Used for 'phase_space' transfers. Used for output if to_stored = False.

    who : str, optional
        Who to set. Possibilities are: Groups: 'all', ' ' (default and same as 'all') Note: This excludes all
        'start' sets., 'twiss', 'a_twiss', 'b_twiss', 'cmat', 'x_dispersion', 'y_dispersion', 'dispersion',
        'chromatic', 'orbit', 'phase_space', 'spin', 'x_plane', 'y_plane', 'z_plane', 'start', 'start_spin',
        'start_phase_space', Individula Parameters: 'x', 'px', 'cmat_11', etc.

    Returns
    -------
    is_ok : bool
        logical
    """

def transfer_lat(lat1: _pybmad.LatStruct) -> _pybmad.LatStruct:
    """
    Wrapper for Fortran routine transfer_lat

    Parameters
    ----------
    lat1 : LatStruct

    Returns
    -------
    lat2 : LatStruct
    """

def transfer_lat_parameters(lat_in: _pybmad.LatStruct) -> _pybmad.LatStruct:
    """
    Wrapper for Fortran routine transfer_lat_parameters

    Parameters
    ----------
    lat_in : LatStruct
        Input lat.

    Returns
    -------
    lat_out : LatStruct
        Output lat with parameters set.
    """

def transfer_map_calc(lat: _pybmad.LatStruct, orb_map: _pybmad.TaylorStructArray1D, ix1: int | None = None, ix2: int | None = None, ref_orb: _pybmad.CoordStruct | None = None, ix_branch: int | None = None, one_turn: bool | None = None, unit_start: bool | None = None, concat_if_possible: bool | None = None, spin_map: _pybmad.TaylorStructArray1D | None = None) -> bool:
    """
    Wrapper for Fortran routine transfer_map_calc

    Parameters
    ----------
    lat : LatStruct
        Lattice used in the calculation.

    orb_map : 1D array of TaylorStruct
        Initial map (used when unit_start = False)
        This parameter is an input/output and is modified in-place.
        As an output, orb_map: Transfer map.

    ix1 : int, optional
        Element start index for the calculation. Default is 0.

    ix2 : int, optional
        Element end index for the calculation. Default is lat.n_ele_track.

    ref_orb : CoordStruct, optional
        Reference orbit/particle at s1 around which the map is made. This arg is needed if: unit_start = True or
        particle is not the same as the reference particle of the lattice.

    ix_branch : int, optional
        Lattice branch index. Default is 0.

    one_turn : bool, optional
        If present and True, and if ix1 = ix2, and the lattice is circular, then construct the one-turn map from
        ix1 back to ix1. Default = False.

    unit_start : bool, optional
        If present and False then orb_map will be used as the starting map instead of the unit map. Default = True

    concat_if_possible : bool, optional
        If present and True then use map concatenation rather than tracking if a map is present for a given
        lattice element. See above. Default is False.

    spin_map : 1D array of TaylorStruct, optional
        Input quaternion spin map. Output only computed if bmad_com.spin_tracking_on = T
        This parameter is an input/output and is modified in-place.
        As an output, spin_map: Quaternion spin map.

    Returns
    -------
    err_flag : bool
        Set True if problem like number overflow, etc.
    """

class TransferMapFromSToS:
    """transfer_map_from_s_to_s return type"""

    @property
    def ref_orb_out(self) -> _pybmad.CoordStruct: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def transfer_map_from_s_to_s(lat: _pybmad.LatStruct, t_map: _pybmad.TaylorStructArray1D, s1: float | None = None, s2: float | None = None, ref_orb_in: _pybmad.CoordStruct | None = None, ix_branch: int | None = None, one_turn: bool | None = None, unit_start: bool | None = None, concat_if_possible: bool | None = None, spin_map: _pybmad.TaylorStructArray1D | None = None) -> TransferMapFromSToS:
    """
    one_turn, unit_start, err_flag, concat_if_possible, spin_map)

    Subroutine to calculate the transfer map between longitudinal positions s1 to s2.

    If s2 < s1 and lat%param%geometry is closed$ then the
    calculation will 'wrap around' the lattice end.
    For example, if s1 = 900 and s2 = 10 then the t_map is the map from
    element 900 to the lattice end plus from 0 through 10.

    If s2 < s1 and lat%param%geometry is open$ then the inverse of the forward map of s2 -> s1 is computed.

    If s2 = s1 then you get the unit map except if one_turn = True and the lattice is circular.

    Parameters
    ----------
    lat : LatStruct
        Lattice used in the calculation.

    t_map : 1D array of TaylorStruct
        Initial map (used when unit_start = False)
        This parameter is an input/output and is modified in-place.
        As an output, t_map: Transfer map.

    s1 : float, optional
        Element start position for the calculation. Default is 0.

    s2 : float, optional
        Element end position for the calculation. Default is lat.param.total_length.

    ref_orb_in : CoordStruct, optional
        Reference orbit/particle at s1 around which the map is made. This arg is needed if: unit_start = True or
        particle is not the same as the reference particle of the lattice.

    ix_branch : int, optional
        Lattice branch index. Default is 0 (main branch).

    one_turn : bool, optional
        If present and True, and s1 = s2, and the lattice is circular: Construct the one-turn map from s1 back to
        s1. Otherwise t_map is unchanged or the unit map if unit_start = T. Default = False.

    unit_start : bool, optional
        If present and False then t_map will be used as the starting map instead of the unit map. Default = True

    concat_if_possible : bool, optional
        If present and True then use map concatenation rather than tracking if a map is present for a given
        lattice element. See above. Default is False.

    spin_map : 1D array of TaylorStruct, optional
        Initial spin map.
        This parameter is an input/output and is modified in-place.
        As an output, spin_map: Final spin map. Only computed if bmad_com.spin_tracking_on = T.

    Returns
    -------
    ref_orb_out : CoordStruct, optional
        Ending coordinates of the reference orbit. This is also the actual orbit of particle

    err_flag : bool, optional
        Set true if there is an error. False otherwise.
    """

def transfer_mat2_from_twiss(twiss1: _pybmad.TwissStruct, twiss2: _pybmad.TwissStruct) -> list[list[float]]:
    """
    Wrapper for Fortran routine transfer_mat2_from_twiss

    Parameters
    ----------
    twiss1 : TwissStruct
        Twiss parameters at the initial point.

    twiss2 : TwissStruct
        Twiss parameters at the end point.

    Returns
    -------
    mat : 2D array of float (shape: 2,2)
        Transfer matrix between the two points.
    """

def transfer_mat_from_twiss(ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct, orb1: Sequence[float], orb2: Sequence[float]) -> list[list[float]]:
    """
    Wrapper for Fortran routine transfer_mat_from_twiss

    Parameters
    ----------
    ele1 : EleStruct
        Element with twiss parameters for the starting point.

    ele2 : EleStruct
        Element with twiss parameters for the ending point.

    orb1 : 1D array of float (shape: 6)
        Reference orbit at ele1 (affects m(i,6) dispersion terms).

    orb2 : 1D array of float (shape: 6)
        Reference orbit at ele2 (affects m(i,6) dispersion terms).

    Returns
    -------
    m : 2D array of float (shape: 6,6)
        Transfer matrix between the two points.
    """

def transfer_matrix_calc(lat: _pybmad.LatStruct, xfer_mat: Sequence[Sequence[float]], xfer_vec: Sequence[float] | None = None, ix1: int | None = None, ix2: int | None = None, ix_branch: int | None = None, one_turn: bool | None = None) -> None:
    """
    Wrapper for Fortran routine transfer_matrix_calc

    Parameters
    ----------
    lat : LatStruct
        Lattice used in the calculation.

    xfer_mat : 2D array of float (shape: 6,6)

    xfer_vec : 1D array of float (shape: 6), optional

    ix1 : int, optional
        Element start index for the calculation. Default is 0.

    ix2 : int, optional
        Element end index for the calculation. Defaults: If ix1 is not present: ix2 = lat.n_ele_track If ix1 is
        present and lattice is closed: Calculate the one-turn matrix from ix1 back to ix1.

    ix_branch : int, optional
        Branch index. Default is 0.

    one_turn : bool, optional
        If present and True, and ix1 = ix2, and the lattice is closed: Construct the one-turn matrix from ix1 back
        to ix1. If False, (the default), and ix1 = ix2, mat6 is the unit matrix.
    """

def transfer_twiss(ele_in: _pybmad.EleStruct, reverse: bool | None = None) -> _pybmad.EleStruct:
    """
    Wrapper for Fortran routine transfer_twiss

    Parameters
    ----------
    ele_in : EleStruct
        Element with existing Twiss parameters.

    reverse : bool, optional
        Reverse alpha and coupling as if particle is going in the reversed direction? Default is False.

    Returns
    -------
    ele_out : EleStruct
        Element receiving the Twiss parameters.
    """

def transfer_wake(wake_in: _pybmad.WakeStruct) -> _pybmad.WakeStruct | None:
    """
    Wrapper for Fortran routine transfer_wake

    Parameters
    ----------
    wake_in : WakeStruct
        Input wake.

    Returns
    -------
    wake_out : WakeStruct, optional
        Output wake.
    """

def truncate_complex_taylor_to_order(complex_taylor_in: _pybmad.ComplexTaylorStructArray1D, order: int, complex_taylor_out: _pybmad.ComplexTaylorStructArray1D) -> None:
    """
    Subroutine to throw out all terms in a complex_taylor map that are above a certain order.

    Parameters
    ----------
    complex_taylor_in : 1D array of ComplexTaylorStruct
        Input complex_taylor map.

    order : int
        Order above which terms are dropped.

    complex_taylor_out : 1D array of ComplexTaylorStruct
        Truncated complex_taylor map.
    """

class Twiss1Propagate:
    """twiss1_propagate return type"""

    @property
    def twiss2(self) -> _pybmad.TwissStruct: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def twiss1_propagate(twiss1: _pybmad.TwissStruct, mat2: Sequence[Sequence[float]], ele_key: int, length: float) -> Twiss1Propagate:
    """
    Wrapper for Fortran routine twiss1_propagate

    Parameters
    ----------
    twiss1 : TwissStruct
        Input Twiss parameters.

    mat2 : 2D array of float (shape: 2,2)
        The transfer matrix.

    ele_key : int
        quadrupole$, etc.

    length : float
        Determines whether the phase is increasing or decreasing.

    Returns
    -------
    twiss2 : TwissStruct
        Output Twiss parameters.

    err : bool
        Set True if there is an error, false otherwise.
    """

def twiss3_at_start(lat: _pybmad.LatStruct, err_flag: bool, ix_branch: int | None = None) -> list[float]:
    """
    Subroutine to calculate the 3D twiss parameters of the three modes of the full 6D 1-turn transfer matrix.
    This routine is for lattices with closed geometries. For open lattices see: twiss3_from_twiss2.

    Note: The rf must be on for this calculation.

    Parameters
    ----------
    lat : LatStruct
        Lattice with

    ix_branch : int, optional
        Branch index. 0 = default.

    Returns
    -------
    tune3 : 1D array of float (shape: 3), optional
        Normal mode tunes
    """

def twiss3_from_twiss2(ele: _pybmad.EleStruct) -> None:
    """
    Routine to calculate the 3D Twiss parameters given the 2D transverse Twiss parameters and some
    longitudinal parameters.
    Also see: twiss3_at_start

    Parameters
    ----------
    ele : EleStruct
        Lattice element at which the calculation is made.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element
    """

def twiss3_propagate1(ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct, err_flag: bool) -> None:
    """
    Subroutine to propagate the twiss parameters using all three normal modes.
    Subroutine from original mode3_mod.
    """

def twiss3_propagate_all(lat: _pybmad.LatStruct, ix_branch: int | None = None) -> None:
    """
    Subroutine to propagate the twiss parameters using all three normal modes.
    Subroutine from original mode3_mod.

    Parameters
    ----------
    lat : LatStruct
        Lattice

    ix_branch : int, optional
        : Branch index. 0 = default.
    """

@overload
def twiss_and_track(lat: _pybmad.LatStruct, orb_array: _pybmad.CoordArrayStructAlloc1D, print_err: bool | None = None, calc_chrom: bool | None = None, use_particle_start: bool | None = None) -> int:
    """
    This routine is an overloaded name for:
      Subroutine twiss_and_track_branch (lat, orb, status, ix_branch, print_err, calc_chrom, orb_start, use_particle_start)
      Subroutine twiss_and_track_all (lat, orb_array, status, print_err, calc_chrom, use_particle_start)

    Routine to calculate the twiss parameters, transport matrices and orbit.

    The essential difference between these two procedures is that
    twiss_and_track_branch only does the main branch while twiss_and_track_all
    does everything but the photon_fork elements.

    Note: This is not necessarily the fastest way to do things since this
    routine does the entire calculation from scratch.

    For a circular ring: If the RF is on, the computed orbit will be the 6D closed orbit.
    If the RF is off, the 4D transverse closed orbit using orbi(0)%vec(6) is computed.

    For an open lattice, the orbit will be computed using orb(0) as
    starting conditions.

    If there is a problem the status argument settings are: in_stop_band$,
    unstable$, non_symplectic$, in_stop_band$, non_symplectic$, xfer_mat_clac_failure$,
    twiss_propagate_failure$, no_complete_orbit$, or no_closed_orbit$. Note: in_stop_band$, unstable$,
    and non_symplectic$ refer to the 1-turn matrix which is computed with closed lattices.
    For an open geometry branch, status = no_complete_orbit$ is for
    where the particle is lost in tracking. A negative sign is used to differentiate an
    error occuring in the first call to twiss_at_start from the second call to twiss_at_start.

    If there is a problem in an open geometry branch, status argument setting is -N where N is the element
    where the particle was lost in tracking (negative numbers are used here to avoid confusion with ok$
    which is mapped to 1.

    Parameters
    ----------
    lat : LatStruct
        lattice.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lat with computed twiss parameters.

    orb_array : 1D array of CoordArrayStruct
        Array of orbit arrays.
        As an output, orb_array: Used as the starting point for an open lattice.

    print_err : bool, optional
        Default is True. If False, suppress error messages.

    calc_chrom : bool, optional
        Default is False. If True, calculate the chromatic functions.

    use_particle_start : bool, optional
        Default is False. If True use branch.particle_start for the starting orbit. Do not use both this and
        orb_start.

    Returns
    -------
    status : int, optional
        Set ok$ if everything is OK and set to something else otherwise. See above for more details.
    """

@overload
def twiss_and_track(lat: _pybmad.LatStruct, orb: _pybmad.CoordStructAlloc1D, ix_branch: int | None = None, print_err: bool | None = None, calc_chrom: bool | None = None, orb_start: _pybmad.CoordStruct | None = None, use_particle_start: bool | None = None) -> int:
    """
    This routine is an overloaded name for:
      Subroutine twiss_and_track_branch (lat, orb, status, ix_branch, print_err, calc_chrom, orb_start, use_particle_start)
      Subroutine twiss_and_track_all (lat, orb_array, status, print_err, calc_chrom, use_particle_start)

    Routine to calculate the twiss parameters, transport matrices and orbit.

    The essential difference between these two procedures is that
    twiss_and_track_branch only does the main branch while twiss_and_track_all
    does everything but the photon_fork elements.

    Note: This is not necessarily the fastest way to do things since this
    routine does the entire calculation from scratch.

    For a circular ring: If the RF is on, the computed orbit will be the 6D closed orbit.
    If the RF is off, the 4D transverse closed orbit using orbi(0)%vec(6) is computed.

    For an open lattice, the orbit will be computed using orb(0) as
    starting conditions.

    If there is a problem the status argument settings are: in_stop_band$,
    unstable$, non_symplectic$, in_stop_band$, non_symplectic$, xfer_mat_clac_failure$,
    twiss_propagate_failure$, no_complete_orbit$, or no_closed_orbit$. Note: in_stop_band$, unstable$,
    and non_symplectic$ refer to the 1-turn matrix which is computed with closed lattices.
    For an open geometry branch, status = no_complete_orbit$ is for
    where the particle is lost in tracking. A negative sign is used to differentiate an
    error occuring in the first call to twiss_at_start from the second call to twiss_at_start.

    If there is a problem in an open geometry branch, status argument setting is -N where N is the element
    where the particle was lost in tracking (negative numbers are used here to avoid confusion with ok$
    which is mapped to 1.

    Parameters
    ----------
    lat : LatStruct
        lattice.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lat with computed twiss parameters.

    orb : 1D array of CoordStruct
        Orbit to be computed
        As an output, orb: Initial conditions to be used for an open geometry lattices.
        As an output, orb: Energy at which the closed orbit is computed.
        This parameter is an input/output and is modified in-place.
        As an output, orb: Computed orbit.

    ix_branch : int, optional
        Branch to track.

    print_err : bool, optional
        Default is True. If False, suppress error messages.

    calc_chrom : bool, optional
        Default is False. If True, calculate the chromatic functions.

    orb_start : CoordStruct, optional
        If present, use this as the starting orbit.

    use_particle_start : bool, optional
        Default is False. If True use branch.particle_start for the starting orbit. Do not use both this and
        orb_start.

    Returns
    -------
    status : int, optional
        Set ok$ if everything is OK and set to something else otherwise. See above for more details.
    """

def twiss_and_track_at_s(lat: _pybmad.LatStruct, s: float, ele_at_s: _pybmad.EleStruct | None = None, orb: _pybmad.CoordStructArray1D | None = None, orb_at_s: _pybmad.CoordStruct | None = None, ix_branch: int | None = None, use_last: bool | None = None, compute_floor_coords: bool | None = None) -> bool:
    """
    Subroutine to return the twiss parameters and particle orbit at a
    given longitudinal position.

    When calculating the Twiss parameters, this routine assumes
    that the lattice elements already contain the Twiss parameters calculated
    for the ends of the elements.

    Additionally, the orbit at the ends of the elements (contained in orb(:)) must be
    precomputed when orb_at_s is present.

    Precomputation of Twiss and orbit at the element ends may be done with the twiss_and_track routine.

    See also:
      twiss_and_track_from_s_to_s
      twiss_and_track_intra_ele

    Parameters
    ----------
    lat : LatStruct
        Lattice.

    s : float
        Longitudinal position. If s is negative the the position is taken to be lat.param.total_length - s.

    ele_at_s : EleStruct, optional
        If the use_last argument is True, ele_at_s is taken to contain valid Twiss parameters stored from a
        previous call to this routine.
        This parameter is an input/output and is modified in-place.
        As an output, ele_at_s: Element structure holding the Twiss parameters.

    orb : 1D array of CoordStruct, optional
        Orbit through the Lattice.

    orb_at_s : CoordStruct, optional
        If the use_last argument is True, orb_at_s is taken to contain the valid orbit stored from a previous
        call.
        This parameter is an input/output and is modified in-place.
        As an output, orb_at_s: Particle position at the position s.

    ix_branch : int, optional
        Branch index, Default is 0 (main lattice).

    use_last : bool, optional
        If present and True, and if ele_at_s.s < s, then use ele_at_s and orb_at_s as the starting point for the
        present calculation. This can speed things up when the present s-position is in the middle of a long
        complicated element and the tracking (EG: Runge-Kutta) is slow.

    compute_floor_coords : bool, optional
        If present and True then the global "floor" coordinates (without misalignments) will be calculated and put
        in ele_at_s.floor.

    Returns
    -------
    err : bool, optional
        Set True if there is a problem in the calculation, False otherwise.
    """

class TwissAndTrackFromSToS:
    """twiss_and_track_from_s_to_s return type"""

    @property
    def orbit_end(self) -> _pybmad.CoordStruct: ...

    @property
    def ele_end(self) -> _pybmad.EleStruct: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def twiss_and_track_from_s_to_s(branch: _pybmad.BranchStruct, orbit_start: _pybmad.CoordStruct, s_end: float, ele_start: _pybmad.EleStruct | None = None, compute_floor_coords: bool | None = None, compute_twiss: bool | None = None) -> TwissAndTrackFromSToS:
    """
    Wrapper for Fortran routine twiss_and_track_from_s_to_s

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch to track through.

    orbit_start : CoordStruct
        Starting phase space coordinates at s_start.

    s_end : float
        Ending position.

    ele_start : EleStruct, optional
        Holds the starting parameters at s_start.

    compute_floor_coords : bool, optional
        If present and True then the global "floor" coordinates will be calculated and put in ele_end.floor.

    compute_twiss : bool, optional
        Default True. If False, to save a little time, do not compute Twiss parameters.

    Returns
    -------
    orbit_end : CoordStruct
        End phase space coordinates.

    ele_end : EleStruct, optional
        Holds the ending Twiss parameters and the transfer matrix. If present then the ele_start argument must
        also be present.

    err : bool, optional
        Set True if there is a problem like the particle gets lost in tracking
    """

class TwissAndTrackIntraEle:
    """twiss_and_track_intra_ele return type"""

    @property
    def orbit_end(self) -> _pybmad.CoordStruct: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def twiss_and_track_intra_ele(ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct, l_start: float, l_end: float, track_upstream_end: bool, track_downstream_end: bool, orbit_start: _pybmad.CoordStruct | None = None, ele_start: _pybmad.EleStruct | None = None, ele_end: _pybmad.EleStruct | None = None, compute_floor_coords: bool | None = None, compute_twiss: bool | None = None, reuse_ele_end: bool | None = None) -> TwissAndTrackIntraEle:
    """
    Wrapper for Fortran routine twiss_and_track_intra_ele

    Parameters
    ----------
    ele : EleStruct
        Element to track through.

    param : LatParamStruct

    l_start : float
        Start position measured from the beginning of the element.

    l_end : float
        Stop position measured from the beginning of the element.

    track_upstream_end : bool
        If True then entrance effects are included in the tracking. But only if l_start = 0 and
        orbit_start.location /= inside$.

    track_downstream_end : bool
        If True then exit effects are included in the tracking but only if l_end = ele.value(l$) (within
        bmad_com.significant_length tol)

    orbit_start : CoordStruct, optional
        Starting phase space coordinates at l_start.

    ele_start : EleStruct, optional
        Holds the starting Twiss parameters at l_start.

    ele_end : EleStruct, optional
        If reuse_ele_end is set True then reuse ele_end from trancking instead of recomputing ele_end from
        scratch. This can save time.
        This parameter is an input/output and is modified in-place.
        As an output, ele_end: Holds the ending Twiss parameters at l_end (except for photons).

    compute_floor_coords : bool, optional
        If present and True then the global "floor" coordinates (without misalignments) will be calculated and put
        in ele_end.floor.

    compute_twiss : bool, optional
        Default True. If False, to save a little time, do not compute Twiss parameters. Also if ele_start is not
        present, no Twiss parameters are computed.

    reuse_ele_end : bool, optional
        If present and True, and if ele_end has the correct lonigitudianal length and key type, reuse ele_end from
        trancking instead of recomputing ele_end from scratch. This can save time.

    Returns
    -------
    orbit_end : CoordStruct, optional
        End phase space coordinates. If present then the orbit_start argument must also be present.

    err : bool, optional
        Set True if there is a problem like the particle gets lost in tracking
    """

class TwissAtElement:
    """twiss_at_element return type"""

    @property
    def start(self) -> _pybmad.EleStruct: ...

    @property
    def end(self) -> _pybmad.EleStruct: ...

    @property
    def average(self) -> _pybmad.EleStruct: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def twiss_at_element(ele: _pybmad.EleStruct) -> TwissAtElement:
    """
    Wrapper for Fortran routine twiss_at_element

    Parameters
    ----------
    ele : EleStruct
        Element to be averaged

    Returns
    -------
    start : EleStruct, optional
        Twiss and s at start of element.

    end : EleStruct, optional
        Twiss and s at end of element.

    average : EleStruct, optional
        Average Twiss and s of element.
    """

def twiss_at_start(lat: _pybmad.LatStruct, ix_branch: int | None = None, type_out: bool | None = None) -> int:
    """
    Wrapper for Fortran routine twiss_at_start

    Parameters
    ----------
    lat : LatStruct
        Lat
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with twiss parameters computed.

    ix_branch : int, optional
        Branch to use. Default is 0 (main branch).

    type_out : bool, optional
        If True (the default), print an error message If the 1-turn matrix is unstable.

    Returns
    -------
    status : int, optional
        Calculation status: ok$, in_stop_band$, unstable$, or non_symplectic$
    """

class TwissFromTracking:
    """twiss_from_tracking return type"""

    @property
    def symp_err(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def twiss_from_tracking(lat: _pybmad.LatStruct, ref_orb0: _pybmad.CoordStruct, d_orb: _pybmad.RealArray1D | None = None) -> TwissFromTracking:
    """
    Wrapper for Fortran routine twiss_from_tracking

    Parameters
    ----------
    lat : LatStruct
        Lat to track through.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Structure holding the Twiss parameters.

    ref_orb0 : CoordStruct
        Reference orbit at lat.ele(0).

    d_orb : 1D array of float, optional
        Vector of offsets to use. If not present or zero bmad_com.d_orb(:) will be used.

    Returns
    -------
    symp_err : float
        A measure of how symplectic the constructed matrices were before symplecitification. mat_symp_check for
        more details.

    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def twiss_propagate1(ele1: _pybmad.EleStruct, ele2: _pybmad.EleStruct, forward: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine twiss_propagate1

    Parameters
    ----------
    ele1 : EleStruct
        Element holding the starting Twiss parameters for forwards propagation.
        This parameter is an input/output and is modified in-place.
        As an output, ele1: Element for the ending Twiss parameters for backwards propagation.

    ele2 : EleStruct
        Element holding the transfer matrix and, if backwards propagation, the starting Twiss.
        This parameter is an input/output and is modified in-place.
        As an output, ele2: Element for the ending Twiss parameters for forward propagation.

    forward : bool, optional
        Default is True. If false, propagate the Twiss backwards.

    Returns
    -------
    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

def twiss_propagate_all(lat: _pybmad.LatStruct, ix_branch: int | None = None, ie_start: int | None = None, ie_end: int | None = None) -> bool:
    """
    Wrapper for Fortran routine twiss_propagate_all

    Parameters
    ----------
    lat : LatStruct
        lattice.
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with parameters computed for the branch.

    ix_branch : int, optional
        Branch index. Default is 0 (main lattice).

    ie_start : int, optional
        Starting element index. Default is 0. Note: The first element at which the Twiss parameters are calculated
        is ie_start+1.

    ie_end : int, optional
        Ending element index, Default is branch.n_ele_track.

    Returns
    -------
    err_flag : bool, optional
        Set True if there is an error. False otherwise.
    """

def twiss_to_1_turn_mat(twiss: _pybmad.TwissStruct, phi: float) -> list[list[float]]:
    """
    Wrapper for Fortran routine twiss_to_1_turn_mat

    Parameters
    ----------
    twiss : TwissStruct
        Structure holding the Twiss parameters. .beta .alpha

    phi : float
        Tune in radians.

    Returns
    -------
    mat2 : 2D array of float (shape: 2,2)
        1-turn matrix.
    """

class TypeComplexTaylors:
    """type_complex_taylors return type"""

    @property
    def lines(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def n_lines(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def type_complex_taylors(complex_taylor: _pybmad.ComplexTaylorStructArray1D, max_order: int | None = None, file_id: int | None = None, out_type: str | None = None, clean: bool | None = None) -> TypeComplexTaylors:
    r"""
    Subroutine to print or put in a string array a Bmad taylor map.
    If the lines(:) argument is not present, the element information is printed to the terminal.

    Parameters
    ----------
    complex_taylor : 1D array of ComplexTaylorStruct
        Array of taylors.

    max_order : int, optional
        Maximum order to print.

    file_id : int, optional
        If present, write output to a file with handle file_id.

    out_type : str, optional
        Determins the string to be used for the output type column. \'\' (default)  -> '1', '2', '3', etc. 'PHASE'
        -> 'X', 'Px, 'Y', 'Py', 'Z', 'Pz' 'SPIN'        -> 'S1', 'Sx', 'Sy', 'Sz'  ! If size = 4: quaternion
        representation 'SPIN'        -> 'Sx', 'Sy', 'Sz'  ! If size = 3: 'NONE'        -> No out column Anything
        else -> Use this for the output column.

    clean : bool, optional
        If True then do not include terms whose coefficients are negligible. Default is false

    Returns
    -------
    lines : 1D array of str, optional
        : Character array to hold the output. If not present, the information is printed to the terminal. Char
        width should be 120 or above for out_type = 'PHASE' but can be less for other out_types.

    n_lines : int, optional
        Number of lines in lines(:) that hold valid output. n_lines must be present if lines(:) is.
    """

def type_coord(coord: _pybmad.CoordStruct) -> None:
    """
    Wrapper for Fortran routine type_coord

    Parameters
    ----------
    coord : CoordStruct
        Coordinate
    """

class TypeEle:
    """type_ele return type"""

    @property
    def lines(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def n_lines(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def type_ele(ele: _pybmad.EleStruct, type_zero_attrib: bool | None = None, type_mat6: int | None = None, type_taylor: bool | None = None, twiss_out: int | None = None, type_control: int | None = None, type_wake: bool | None = None, type_floor_coords: bool | None = None, type_field: int | None = None, type_wall: bool | None = None, type_rad_kick: bool | None = None, type_internal: bool | None = None) -> TypeEle:
    """
    Wrapper for Fortran routine type_ele

    Parameters
    ----------
    ele : EleStruct
        Element

    type_zero_attrib : bool, optional
        If False then surpress printing of real attributes whose value is 0 or switch attributes that have their
        default value. Default is False.

    type_mat6 : int, optional
        = 0   => Do not type ele.mat6 = 4   => Type 4X4 xy submatrix = 6   => Type full 6x6 matrix (Default)

    type_taylor : bool, optional
        Print out taylor map terms? If ele.taylor is not allocated then this is ignored. Default is False.

    twiss_out : int, optional
        Print the Twiss parameters at the element end? = 0         => Do not print the Twiss parameters = radians$
        => Print Twiss, phi in radians (Default). = degrees$  => Print Twiss, phi in degrees. = cycles$   => Print
        Twiss, phi in radians/2pi.

    type_control : int, optional
        Print control status? If ele.branch.lat is not associated cannot print status info. = no$      => One line
        of info. = short$   => Almost all info except long knot point lists are truncated (default). = all$     =>
        Everything.

    type_wake : bool, optional
        If True then print the long-range and short-range wakes information. If False then just print how many
        terms the wake has. Default is True. If ele.wake is not allocated then this is ignored.

    type_floor_coords : bool, optional
        Default is False. If True then print the global ("floor") coordinates at the exit end of the element.

    type_field : int, optional
        Print field maps? = no$      => One line of info (default). = short$   => Header info. No tables. = all$
        => Everything.

    type_wall : bool, optional
        Default is False. If True, print wall info.

    type_rad_kick : bool, optional
        Default is False. If True, print synch rad kick info.

    type_internal : bool, optional
        Default is False. If True, print some internal parameters.

    Returns
    -------
    lines : 1D array of str, optional
        Character array to hold the output. If not present, the information is printed to the terminal.

    n_lines : int, optional
        Number of lines in lines(:) that hold valid output. n_lines must be present if lines(:) is.
    """

def type_end_stuff(li: _pybmad.CharacterAlloc1D, nl: int, lines: _pybmad.CharacterAlloc1D | None = None, n_lines: int | None = None) -> None:
    """
    Wrapper for Fortran routine type_end_stuff

    Parameters
    ----------
    li : 1D array of str

    nl : int

    lines : 1D array of str, optional

    n_lines : int, optional
    """

def type_expression_tree(tree: _pybmad.ExpressionTreeStruct, indent: int | None = None) -> None:
    """
    Routine to print an expression tree in tree form.
    Good for debugging.

    Parameters
    ----------
    tree : ExpressionTreeStruct
        Tree to print.

    indent : int, optional
        Initial indent. Default is zero.
    """

class TypePtcFibre:
    """type_ptc_fibre return type"""

    @property
    def lines(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def n_lines(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def type_ptc_fibre(ptc_fibre: _pybmad.Fibre, print_coords: bool | None = None) -> TypePtcFibre:
    """
    Routine to put information on a PTC fibre element into a string array.
    If "lines" is not present, the information will be printed to the screen.

    Parameters
    ----------
    ptc_fibre : Fibre
        Fibre to type info of.

    print_coords : bool, optional
        If True then print coordinate and patch information. Default is True.

    Returns
    -------
    lines : 1D array of str, optional
        Character array to hold the output.

    n_lines : int, optional
        Number of lines used in lines(:)
    """

def type_ptc_layout(lay: _pybmad.Layout) -> None:
    """Subroutine to print the global information in a layout"""

class TypeTaylors:
    """type_taylors return type"""

    @property
    def n_lines(self) -> int | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def type_taylors(bmad_taylor: _pybmad.TaylorStructArray1D, max_order: int | None = None, lines: _pybmad.CharacterAlloc1D | None = None, n_lines: int | None = None, file_id: int | None = None, out_style: str | None = None, clean: bool | None = None, out_var_suffix: str | None = None, append: bool | None = None) -> TypeTaylors:
    r"""
    Wrapper for Fortran routine type_taylors

    Parameters
    ----------
    bmad_taylor : 1D array of TaylorStruct
        Array of taylors.

    max_order : int, optional
        Maximum order to print.

    lines : 1D array of str, optional
        Used with append = True. Output will start at n_lines+1
        This parameter is an input/output and is modified in-place.
        As an output, lines: Array to hold the output.

    n_lines : int, optional
        Used with append = True. Output will start at n_lines+1.
        This parameter is an input/output and is modified in-place.
        As an output, n_lines: Number of lines in lines(

    file_id : int, optional
        If present, write output to a file with handle file_id.

    out_style : str, optional
        Determins the string to be used for the output type column. \'\' (default) -> 'X', 'Px, 'Y', 'Py', 'Z', 'Pz'
        If size(bmad_taylor) = 6 -> 'S1', 'Sx', 'Sy', 'Sz' (Spin quaternions) If size(bmad_taylor) = 4 -> '1',
        '2', '3', etc. Otherwiase 'NUMBER'     -> '1', '2', '3', etc. 'BMAD'       -> Bmad lattice file format.

    clean : bool, optional
        If True then do not include terms whose coefficients are negligible. Default is false.

    out_var_suffix : str, optional
        If out_style = 'SCIBMAD', out_var_suffix is used as the suffix of the variable holding the taylor map.
        Default is "z". For example, if "z" is the suffix then: Descriptor = "d_z", orbital map name = "v_z", ref
        orbit name = v0_z, and spin map name = "q_z".

    append : bool, optional
        Default is False. If True, n_lines on input is the number of existing lines in lines(:) to save.

    Returns
    -------
    n_lines : int, optional
        Used with append = True. Output will start at n_lines+1.
        This parameter is an input/output and is modified in-place.
        As an output, n_lines: Number of lines in lines(
    """

def type_twiss(ele: _pybmad.EleStruct, frequency_units: int | None = None, compact_format: bool | None = None, lines: _pybmad.CharacterAlloc1D | None = None) -> int:
    """
    Wrapper for Fortran routine type_twiss

    Parameters
    ----------
    ele : EleStruct
        Element containing the Twiss parameters.

    frequency_units : int, optional
        Units for phi: = radians$  => Type Twiss, use radians for phi (Default). = degrees$  => Type Twiss, use
        degrees for phi. = cycles$   => Type Twiss, use cycles (1 = 2pi) units.

    compact_format : bool, optional
        If present and True then use a compact output form.

    lines : 1D array of str, optional
        : Character array to hold the output. The string length should be at least 120 characters. 13 lines are
        needed for the verbose form. If not present, the information is printed to the terminal.

    Returns
    -------
    n_lines : int, optional
        Number of lines in lines(:) that hold valid output. n_lines must be present if lines(:) is.
    """

def update_ele_from_fibre(ele: _pybmad.EleStruct) -> None:
    """
    Routine to update a bmad lattice element when the associated PTC fibre has been modified.
    Remember to call lattice_bookkeeper after calling this routine.

    Parameters
    ----------
    ele : EleStruct
        Element with corresponding ele.ptc_fibre fibre.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Modified element.
    """

def update_fibre_from_ele(ele: _pybmad.EleStruct) -> bool:
    """
    Wrapper for Fortran routine update_fibre_from_ele

    Parameters
    ----------
    ele : EleStruct
        Element with corresponding PTC fibre.

    Returns
    -------
    survey_needed : bool
        Set True if a call to survey will be needed. Calling survey is avoided in this routine to save time if
        multiple elements are being updated.
    """

def update_floor_angles(floor: _pybmad.FloorPositionStruct, floor0: _pybmad.FloorPositionStruct | None = None) -> None:
    """
    Wrapper for Fortran routine update_floor_angles

    Parameters
    ----------
    floor : FloorPositionStruct
        Position with input w matrix.
        This parameter is an input/output and is modified in-place.
        As an output, floor: Position with output angles.

    floor0 : FloorPositionStruct, optional
        Reference position. There are two solutions related by: [theta, phi, psi] & [pi+theta, pi-phi, pi+psi] If
        floor0 is present, choose the solution "nearest" the angles in floor0.
    """

def valid_field_calc(ele: _pybmad.EleStruct, field_calc: int) -> bool:
    """
    Wrapper for Fortran routine valid_field_calc

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    field_calc : int
        bmad_standard$, etc.

    Returns
    -------
    is_valid : bool
        True if a valid method. False otherwise.
    """

def valid_fringe_type(ele: _pybmad.EleStruct, fringe_type: int) -> bool:
    """
    Wrapper for Fortran routine valid_fringe_type

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    fringe_type : int
        bmad_standard$, etc.

    Returns
    -------
    is_valid : bool
        True if a valid method. False otherwise.
    """

def valid_mat6_calc_method(ele: _pybmad.EleStruct, species: int, mat6_calc_method: int) -> bool:
    """
    Wrapper for Fortran routine valid_mat6_calc_method

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    species : int
        Type of particle being tracked. electron$, etc. or not_set$

    mat6_calc_method : int
        bmad_standard$, etc.

    Returns
    -------
    is_valid : bool
        True if a valid method. False otherwise.
    """

def valid_spin_tracking_method(ele: _pybmad.EleStruct, spin_tracking_method: int) -> bool:
    """
    Wrapper for Fortran routine valid_spin_tracking_method

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    spin_tracking_method : int
        bmad_standard$, etc.

    Returns
    -------
    is_valid : bool
        True if a valid method. False otherwise.
    """

def valid_tracking_method(ele: _pybmad.EleStruct, species: int, tracking_method: int) -> bool:
    """
    Wrapper for Fortran routine valid_tracking_method

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    species : int
        Type of particle being tracked. electron$, etc. or not_set$

    tracking_method : int
        bmad_standard$, etc.

    Returns
    -------
    is_valid : bool
        True if a valid method. False otherwise.
    """

class ValueOfAttribute:
    """value_of_attribute return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def value_of_attribute(ele: _pybmad.EleStruct, attrib_name: str, err_print_flag: bool | None = None, err_value: float | None = None) -> ValueOfAttribute:
    """
    Wrapper for Fortran routine value_of_attribute

    Parameters
    ----------
    ele : EleStruct
        After this routine finishes Ptr_attrib will point to a variable within this element.

    attrib_name : str
        Name of attribute. Must be uppercase. For example: "HKICK".

    err_print_flag : bool, optional
        If present and True then print an error message if there is an  error.

    err_value : float, optional
        Value to set value argument if there is an error. Default is 0.

    Returns
    -------
    value : float
        Value of the attribute. Set to err_value if not found.

    err_flag : bool, optional
        Set True if attribtute not found. False otherwise.
    """

def value_to_line(line: str, value: float, str: str, typ: str, ignore_if_zero: bool | None = None, use_comma: bool | None = None) -> None:
    """
    Wrapper for Fortran routine value_to_line

    Parameters
    ----------
    line : str

    value : float

    str : str

    typ : str

    ignore_if_zero : bool, optional

    use_comma : bool, optional
    """

def vec_to_polar(vec: Sequence[float], phase: float | None = None) -> _pybmad.SpinPolarStruct:
    """
    Wrapper for Fortran routine vec_to_polar

    Parameters
    ----------
    vec : 1D array of float (shape: 3)
        unitary spin vector

    phase : float, optional
        Phase of the spinor, if not given then set to zero

    Returns
    -------
    polar : SpinPolarStruct
    """

def vec_to_spinor(vec: Sequence[float], phase: float | None = None) -> list[complex]:
    """
    Wrapper for Fortran routine vec_to_spinor

    Parameters
    ----------
    vec : 1D array of float (shape: 3)
        Spin vector in cartesian coordinates

    phase : float, optional
        Phase of the spinor, if not given then set to zero

    Returns
    -------
    spinor : 1D array of complex (shape: 2)
        Spinor.
    """

def verify_valid_name(name: str, ix_name: int, pure_name: bool | None = None, include_wild: bool | None = None) -> bool:
    r"""
    Routine to check if a name is well formed. Examples:
      "0>>Q0"                           -- Invalid (will only be valid after lattice expansion).
      "Q1##1"                           -- Invalid (double hash not accepted).
      "Q2A_C.\7#"                       -- Pure name (no "[", "]", "(", ")", "%" characters present).
      "Q3[GRID_FIELD(1)%FIELD_SCALE]"   -- Valid but not a pure name.
      "RFCAVITY::*"                     -- Valid if include_wild = True.

    This subroutine is used by bmad_parser and bmad_parser2.
    This subroutine is not intended for general use.

    Parameters
    ----------
    name : str
        Name(1:ix_name) is the string to check.

    ix_name : int
        Number of characters in the name.

    pure_name : bool, optional
        If True, reject names that contain "[", "]", "(", ")", "." characters. Default is False.

    include_wild : bool, optional
        Name can include wild card characters and additionally type prefixes like "QUAD::". Default is False.

    Returns
    -------
    is_valid : bool
        True if name is well formed. False otherwise.
    """

def vert_angle_func(integ_prob: float, status: int) -> float:
    """
    Wrapper for Fortran routine vert_angle_func

    Parameters
    ----------
    integ_prob : float

    status : int

    Returns
    -------
    d_angle : float
    """

def w_mat_for_bend_angle(angle: float, ref_tilt: float, r_vec: Sequence[float] | None = None) -> list[list[float]]:
    """
    Wrapper for Fortran routine w_mat_for_bend_angle

    Parameters
    ----------
    angle : float
        Bending angle.

    ref_tilt : float
        Reference tilt.

    r_vec : 1D array of float (shape: 3), optional
        Starting position.
        This parameter is an input/output and is modified in-place.
        As an output, r_vec: position with ref_tilt transformation

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3)
        W matrix
    """

def w_mat_for_tilt(tilt: float, return_inverse: bool | None = None) -> list[list[float]]:
    """
    Wrapper for Fortran routine w_mat_for_tilt

    Parameters
    ----------
    tilt : float
        pitch angle

    return_inverse : bool, optional
        If True, return the inverse matrix. Default is False.

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3)
        Transformation matrix.
    """

def w_mat_for_x_pitch(x_pitch: float, return_inverse: bool | None = None) -> list[list[float]]:
    """
    Wrapper for Fortran routine w_mat_for_x_pitch

    Parameters
    ----------
    x_pitch : float
        pitch angle

    return_inverse : bool, optional
        If True, return the inverse matrix. Default is False.

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3)
        Transformation matrix.
    """

def w_mat_for_y_pitch(y_pitch: float, return_inverse: bool | None = None) -> list[list[float]]:
    """
    Wrapper for Fortran routine w_mat_for_y_pitch

    Parameters
    ----------
    y_pitch : float
        pitch angle

    return_inverse : bool, optional
        If True, return the inverse matrix. Default is False.

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3)
        Transformation matrix.
    """

class Wall3dDRadius:
    """wall3d_d_radius return type"""

    @property
    def perp(self) -> list[float]: ...

    @property
    def ix_section(self) -> int: ...

    @property
    def no_wall_here(self) -> bool: ...

    @property
    def origin(self) -> list[float]: ...

    @property
    def radius_wall(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def d_radius(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def wall3d_d_radius(position: _pybmad.RealArray1D, ele: _pybmad.EleStruct, ix_wall: int | None = None) -> Wall3dDRadius:
    """
    no_wall_here, origin, radius_wall, err_flag) result (d_radius)

    Routine to calculate the difference radius = particle_radius - wall_radius.
    Radiuses are measured along a line from the wall origin with the line passing through
    the particle position.
    The wall origin itself lies on a line connecting the centers of the bounding sections.

    Module needed:
      use wall3d_mod

    Parameters
    ----------
    position : 1D array of float
        Particle position in element coordinates. In a patch, with respect to entrance coords. [position(1),
        position(3)] = [x, y] transverse coords. position(5)                = Longitudinal position relative to
        beginning of element. position(6)                = Longitudinal velocity (only +/- sign matters).

    ele : EleStruct
        Element with wall

    ix_wall : int, optional
        Index of wall in .wall3d(:) array. Default is 1.

    Returns
    -------
    d_radius : float
        r_particle - r_wall

    perp : 1D array of float (shape: 3), optional
        Perpendicular normal to the wall.

    ix_section : int, optional
        Set to wall slice section particle is in. That is between ix_section and ix_section+1.

    no_wall_here : bool, optional
        True if the sub-chamber under consideration does not exist at the longitudinal location of the particle.

    origin : 1D array of float (shape: 3), optional
        (x, y, s) origin with respect to the radius is measured. Uses the same coords as position.

    radius_wall : float, optional
        Radius of the wall.

    err_flag : bool, optional
        Set True if error. (EG noassociated .wall3d), false otherwise.
    """

def wall3d_initializer(wall3d: _pybmad.Wall3DStruct) -> bool:
    """
    Routine to initialize a wall3d_struct
      1) Add vertex points if there is symmetry.
      2) Compute circular and elliptical centers.
      3) Compute spline coefficients, etc.

    Parameters
    ----------
    wall3d : Wall3dStruct
        Wall.
        This parameter is an input/output and is modified in-place.
        As an output, wall3d: Initialized wall.

    Returns
    -------
    err : bool
        Set true if there is a problem.
    """

def wall3d_section_initializer(section: _pybmad.Wall3DSectionStruct) -> bool:
    """
    Routine to initialize a wall3d_section_struct:
      1) Add vertex points if there is symmetry.
      2) Compute circular and elliptical centers.

    Parameters
    ----------
    section : Wall3dSectionStruct
        Wall3d section.
        This parameter is an input/output and is modified in-place.
        As an output, section: Initialized section-section.

    Returns
    -------
    err : bool
        Set true if there is a problem.
    """

def wall3d_to_position(orbit: _pybmad.CoordStruct, ele: _pybmad.EleStruct) -> list[float]:
    """
    Routine to return the suitable postion to be used in calling wall3d_d_radius

    This routine assumes that if in a patch the coordinates of orbit are with respect
    to the downstream end if orbit%direction*orbit%time_dir = 1 and vice versa.

    Parameters
    ----------
    orbit : CoordStruct
        Particle position.

    ele : EleStruct
        Element particle is in.

    Returns
    -------
    position : 1D array of float (shape: 6)
        Position used in wall3d_d_radius call.
    """

def word_to_value(word: str, lat: _pybmad.LatStruct, value: float, err_flag: bool, ele: _pybmad.EleStruct | None = None) -> None:
    """
    Wrapper for Fortran routine word_to_value

    Parameters
    ----------
    word : str

    lat : LatStruct

    value : float

    err_flag : bool

    ele : EleStruct, optional
    """

def write_ascii_beam_file(file_name: str, beam: _pybmad.BeamStruct, new_file: bool | None = None, alive_only: bool | None = None) -> None:
    """
    Routine to write a beam file in ASCII format (version 4).

    Parameters
    ----------
    file_name : str
        Name of file.

    beam : BeamStruct
        Beam to write

    new_file : bool, optional
        New file or append? Default = True.

    alive_only : bool, optional
        Only write live (includes pre_born) particles to the file? Default is False.
    """

def write_astra_bend(iu: int, strength: float, id: int, d1: Sequence[float], d2: Sequence[float], d3: Sequence[float], d4: Sequence[float]) -> None:
    """
    Wrapper for Fortran routine write_astra_bend

    Parameters
    ----------
    iu : int

    strength : float

    id : int

    d1 : 1D array of float (shape: 2)

    d2 : 1D array of float (shape: 2)

    d3 : 1D array of float (shape: 2)

    d4 : 1D array of float (shape: 2)
    """

def write_astra_ele(iu: int, ele: _pybmad.EleStruct, id: int, fieldgrid_names: _pybmad.StrIndexStruct | None = None, dimensions: int | None = None) -> None:
    """
    Wrapper for Fortran routine write_astra_ele

    Parameters
    ----------
    iu : int

    ele : EleStruct

    id : int

    fieldgrid_names : StrIndexStruct, optional

    dimensions : int, optional
    """

class WriteAstraFieldGridFile:
    """write_astra_field_grid_file return type"""

    @property
    def maxfield(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_astra_field_grid_file(astra_file_unit: int, ele: _pybmad.EleStruct, dz: float | None = None) -> WriteAstraFieldGridFile:
    """
    Write 1-D field map files for Astra. The format is:
    z field
    ...

    Note: Simplified from write_opal_field_grid_file

    Parameters
    ----------
    astra_file_unit : int
        unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

    ele : EleStruct
        element to make map

    dz : float, optional
        z step size in m. Default: 0.001 m

    Returns
    -------
    maxfield : float
        absolute maximum found for element field scaling

    err : bool, optional
        Set True if, say a file could not be opened.
    """

class WriteAstraFieldGridFile3d:
    """write_astra_field_grid_file_3d return type"""

    @property
    def maxfield(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_astra_field_grid_file_3d(base_filename: str, ele: _pybmad.EleStruct, dz: float | None = None) -> WriteAstraFieldGridFile3d:
    r"""
    Writes 3-D field map files for Astra. The format is:
    Nx x[1] x[2] ....... x[Nx-1] x[Nx]
    Ny y[1] y[2] ....... y[Ny-1] y[Ny]
    Nz z[1] z[2] ....... z[Nz-1] z[Nz]
    <field values>
    where field values are produced from a loop as in:
    do iz = 1, Nz
      do iy = 1, Ny
        write single line: field(:, iy, iz)

    Note: similar to write_astra_field_grid_file

    Parameters
    ----------
    base_filename : str
        Base filename. Files will be written as: base_filename.ex, .ey, .ez, .bx, .by, .bz If set to \'\', no files
        will be written

    ele : EleStruct
        element to make map

    dz : float, optional
        z step size in m. Default: 0.001 m

    Returns
    -------
    maxfield : float
        absolute maximum on-axis field found for element field scaling

    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_astra_lattice_file(astra_file_unit: int, lat: _pybmad.LatStruct, astra_lattice_param: _pybmad.AstraLatticeParamStruct) -> bool:
    """
    Subroutine to write an Astra lattice file using the information in a lat_struct.

    Parameters
    ----------
    astra_file_unit : int
        unit number to write to

    lat : LatStruct
        Holds the lattice information.

    Returns
    -------
    err : bool
        Set True if, say a file could not be opened.
    """

def write_beam_file(file_name: str, beam: _pybmad.BeamStruct, new_file: bool | None = None, file_format: int | None = None, lat: _pybmad.LatStruct | None = None, alive_only: bool | None = None) -> None:
    """
    Routine to write a beam file.

    A '.h5' suffix will be appended to the created file if hdf5$ format is used and file_name does not
    already have a '.h5' or '.hdf5' suffix.

    Parameters
    ----------
    file_name : str
        Name of file.

    beam : BeamStruct
        Beam to write

    new_file : bool, optional
        New file or append? Default = True.

    file_format : int, optional
        ascii$, or hdf5$ (default). old_ascii$ (deprecated) is still accepted.

    lat : LatStruct, optional
        If present, lattice info will be writen to hdf5 files.

    alive_only : bool, optional
        Only write live (includes pre_born) particles to the file? Default is False.
    """

def write_beam_floor_positions(file_name: str, beam: _pybmad.BeamStruct, ele: _pybmad.EleStruct, new_file: bool | None = None) -> None:
    """
    Wrapper for Fortran routine write_beam_floor_positions

    Parameters
    ----------
    file_name : str
        Name of file.

    beam : BeamStruct
        Beam to write

    ele : EleStruct
        Element that the beam is at.

    new_file : bool, optional
        New file or append? Default = True.
    """

def write_binary_cartesian_map(file_name: str, ele: _pybmad.EleStruct, cart_map: _pybmad.CartesianMapStruct) -> bool:
    """
    Routine to write a binary cartesian_map structure.
    Note: The file name should have a ".bin" suffix.

    Parameters
    ----------
    file_name : str
        File to create.

    ele : EleStruct
        Element associated with the map.

    cart_map : CartesianMapStruct
        Cartesian map.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def write_binary_cylindrical_map(file_name: str, ele: _pybmad.EleStruct, cl_map: _pybmad.CylindricalMapStruct) -> bool:
    """
    Routine to write a binary cylindrical_map structure.
    Note: The file name should have a ".bin" suffix.

    Parameters
    ----------
    file_name : str
        File to create.

    ele : EleStruct
        Element associated with the map.

    cl_map : CylindricalMapStruct
        Cylindrical map.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def write_binary_grid_field(file_name: str, ele: _pybmad.EleStruct, g_field: _pybmad.GridFieldStruct) -> bool:
    """
    Routine to write a binary grid_field structure.
    Note: The file name should have a ".bin" suffix.

    Parameters
    ----------
    file_name : str
        File to create.

    ele : EleStruct
        Element associated with the map.

    g_field : GridFieldStruct
        Cylindrical map.

    Returns
    -------
    err_flag : bool
        Set True if there is an error. False otherwise.
    """

def write_blender_ele(iu: int, ele: _pybmad.EleStruct, old_format: bool | None = None) -> None:
    """
    Wrapper for Fortran routine write_blender_ele

    Parameters
    ----------
    iu : int

    ele : EleStruct

    old_format : bool, optional
    """

def write_blender_lat_layout(file_name: str, lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine write_blender_lat_layout

    Parameters
    ----------
    file_name : str

    lat : LatStruct
    """

def write_bmad_lattice_file(bmad_file: str, lat: _pybmad.LatStruct, output_form: int | None = None, orbit0: _pybmad.CoordStruct | None = None) -> bool:
    """
    Wrapper for Fortran routine write_bmad_lattice_file

    Parameters
    ----------
    bmad_file : str
        Name of the output lattice file.

    lat : LatStruct
        Holds the lattice information.

    output_form : int, optional
        binary$   -> Write grid_field info in binary hdf5 form in separate files. Default. All other fields are
        writen in separate files in ASCII ascii$    -> Fields will be put in separate ASCII files. one_file$ ->
        Everything in one file.

    orbit0 : CoordStruct, optional
        Initial orbit. Used to write the inital orbit if the lattice geometry is closed.

    Returns
    -------
    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_digested_bmad_file(digested_name: str, lat: _pybmad.LatStruct, n_files: int | None = None, file_names: _pybmad.CharacterAlloc1D | None = None, extra: _pybmad.ExtraParsingInfoStruct | None = None) -> bool:
    """
    Wrapper for Fortran routine write_digested_bmad_file

    Parameters
    ----------
    digested_name : str
        Name for the digested file.

    lat : LatStruct
        Input lat structure.

    n_files : int, optional
        Number of original files

    file_names : 1D array of str, optional
        Names of the original files used to create the lat structure.

    extra : ExtraParsingInfoStruct, optional
        Extra info that can be stored in the digested file.

    Returns
    -------
    err_flag : bool, optional
        Set True if there is a problem. EG: No write permission. Set False if everything is OK.
    """

def write_gpt_ele(iu: int, ele: _pybmad.EleStruct, name: str, dimensions: int, fieldgrid_names: _pybmad.StrIndexStruct | None = None, only_phasing: bool | None = None) -> None:
    """
    Wrapper for Fortran routine write_gpt_ele

    Parameters
    ----------
    iu : int

    ele : EleStruct

    name : str

    dimensions : int

    fieldgrid_names : StrIndexStruct, optional

    only_phasing : bool, optional
    """

class WriteGptFieldGridFile1d:
    """write_gpt_field_grid_file_1d return type"""

    @property
    def maxfield(self) -> float: ...

    @property
    def ref_time(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_gpt_field_grid_file_1d(gpt_file_unit: int, ele: _pybmad.EleStruct, dz: float | None = None) -> WriteGptFieldGridFile1d:
    """
    Write 1-D field map files for gpt. The format is:
    z field
    ...

    Note: Simplified from write_opal_field_grid_file

    Parameters
    ----------
    gpt_file_unit : int
        unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

    ele : EleStruct
        element to make map

    dz : float, optional
        z step size in m. Default: 0.001 m

    Returns
    -------
    maxfield : float
        absolute maximum found for element field scaling

    ref_time : float
        time that the field was evaluated at

    err : bool, optional
        Set True if, say a file could not be opened.
    """

class WriteGptFieldGridFile2d:
    """write_gpt_field_grid_file_2d return type"""

    @property
    def maxfield(self) -> float: ...

    @property
    def ref_time(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_gpt_field_grid_file_2d(gpt_file_unit: int, ele: _pybmad.EleStruct, dr: float | None = None, dz: float | None = None, r_max: float | None = None) -> WriteGptFieldGridFile2d:
    """
    Subroutine to write an GPT lattice file using the information in
    a lat_struct. Optionally only part of the lattice can be generated.

    Parameters
    ----------
    gpt_file_unit : int
        unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

    ele : EleStruct
        element to make map

    dr : float, optional
        r step size in m. Default: 0.001 m

    dz : float, optional
        z step size in m. Default: 0.001 m

    r_max : float, optional
        maximum radius in m. Default: 0.02 m

    Returns
    -------
    maxfield : float
        absolute maximum found for element field scaling

    ref_time : float
        time that the field was evaluated at

    err : bool, optional
        Set True if, say a file could not be opened.
    """

class WriteGptFieldGridFile3d:
    """write_gpt_field_grid_file_3d return type"""

    @property
    def maxfield(self) -> float: ...

    @property
    def ref_time(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_gpt_field_grid_file_3d(base_filename: str, ele: _pybmad.EleStruct, dz: float | None = None) -> WriteGptFieldGridFile3d:
    r"""
    Writes 3-D field map files for gpt. The format is:

    E-fields:
    'x', 'y', 'z', 'ExRe', 'EyRe', 'EzRe', 'ExIm ', 'EyIm ', 'EzIm '
    H-fields
    'x', 'y', 'z', 'HxRe', 'HyRe', 'HzRe', 'HxIm ', 'HyIm ', 'HzIm '

    where the fields oscillate as exp(+i \omega t)

    Note: similar to write_gpt_field_grid_file

    Parameters
    ----------
    base_filename : str
        Base filename. Files will be written as: base_filename_E_ASCII.gpt, _H_ASCII.gpt If set to \'\', no files
        will be written

    ele : EleStruct
        element to make map

    dz : float, optional
        z step size in m. Default: 0.001 m

    Returns
    -------
    maxfield : float
        absolute maximum on-axis field found for element field scaling

    ref_time : float
        time that the field was evaluated at

    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_gpt_lattice_file(lat: _pybmad.LatStruct, gpt_lat_param: _pybmad.GptLatParamStruct) -> bool:
    """
    Subroutine to write a gpt lattice file using the information in a Bmad lattice.

    Parameters
    ----------
    lat : LatStruct
        Holds the lattice information.

    gpt_lat_param : GptLatParamStruct
        Parameters for constructing the lattice

    Returns
    -------
    err : bool
        Set True if there is an error
    """

class WriteLatLine:
    """write_lat_line return type"""

    @property
    def line(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_lat_line(line: str, iu: int, end_is_neigh: bool, do_split: bool | None = None, ampersand_at_ends: bool | None = None) -> WriteLatLine:
    """
    Routine to write strings to a lattice file.
    This routine will break the string up into multiple lines
    if the string is too long and add a continuation character if needed.

    If the "line" arg does not represent a full "sentence" (end_is_neigh = False),
    then only part of the line may be written and the part not written will be returned.

    Parameters
    ----------
    line : str
        String of text.
        This parameter is an input/output and is modified in-place.
        As an output, line: part of the string not written.

    iu : int
        Unit number to write to.

    end_is_neigh : bool
        If true then write out everything. Otherwise wait for a full line of max_char characters or so.

    do_split : bool, optional
        Split line if overlength? Default is True. False is used when line has already been split for expressions
        since the expression splitting routine does a much better job of it.

    ampersand_at_ends : bool, optional
        Default True. If False then do not include "&" line continuation

    Returns
    -------
    line : str
        String of text.
        This parameter is an input/output and is modified in-place.
        As an output, line: part of the string not written.
    """

def write_lattice_elegant_format(out_file_name: str, lat: _pybmad.LatStruct, ref_orbit: _pybmad.CoordStructAlloc1D | None = None, use_matrix_model: bool | None = None, include_apertures: bool | None = None, dr12_drift_max: float | None = None, ix_branch: int | None = None) -> bool:
    """
    Wrapper for Fortran routine write_lattice_elegant_format

    Parameters
    ----------
    out_file_name : str
        Name of the mad output lattice file.

    lat : LatStruct
        Holds the lattice information.

    ref_orbit : 1D array of CoordStruct, optional
        Referece orbit for sad_mult and patch elements. This argument must be present if the lattice has sad_mult
        or patch elements and is being translated to MAD-8 or SAD.

    use_matrix_model : bool, optional
        Use a drift-matrix_drift model for wigglers/undulators? [A MAD "matrix" is a 2nd order Taylor map.] This
        switch is ignored for SAD conversion. Default is False -> Use a bend-drift-bend model. Note: sol_quad
        elements always use a drift-matrix-drift model.

    include_apertures : bool, optional
        If True (the default), add to the output lattice a zero length collimator element next to any non-
        collimator element that has an aperture. Note: MADX translations for non-drift elements can handle non-
        collimator elements with an aperture so in this case this argument is ignored.

    dr12_drift_max : float, optional
        Max deviation for drifts allowed before a correction matrix element is added. Default value is 1d-5. A
        negative number means use default.

    ix_branch : int, optional
        Index of lattice branch to use. Default = 0.

    Returns
    -------
    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_lattice_foreign_format(out_type: str, out_file_name: str, lat: _pybmad.LatStruct, ref_orbit: _pybmad.CoordStructAlloc1D | None = None, use_matrix_model: bool | None = None, include_apertures: bool | None = None, dr12_drift_max: float | None = None, ix_branch: int | None = None) -> bool:
    """
    Wrapper for Fortran routine write_lattice_foreign_format

    Parameters
    ----------
    out_type : str
        Either 'ELEGANT', 'MAD-8', 'MAD-X', 'SAD', 'OPAL-T', 'PALS', or 'SCIBMAD'.

    out_file_name : str
        Name of the mad output lattice file.

    lat : LatStruct
        Holds the lattice information.

    ref_orbit : 1D array of CoordStruct, optional
        Referece orbit for sad_mult and patch elements. This argument must be present if the lattice has sad_mult
        or patch elements and is being translated to MAD-8 or SAD.

    use_matrix_model : bool, optional
        Use a drift-matrix_drift model for wigglers/undulators? [A MAD "matrix" is a 2nd order Taylor map.] This
        switch is ignored for SAD conversion. Default is False -> Use a bend-drift-bend model. Note: sol_quad
        elements always use a drift-matrix-drift model.

    include_apertures : bool, optional
        If True (the default), add to the output lattice a zero length collimator element next to any non-
        collimator element that has an aperture. Note: MADX translations for non-drift elements can handle non-
        collimator elements with an aperture so in this case this argument is ignored.

    dr12_drift_max : float, optional
        Max deviation for drifts allowed before a correction matrix element is added. Default value is 1d-5. A
        negative number means use default.

    ix_branch : int, optional
        Index of lattice branch to use. Default = 0.

    Returns
    -------
    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_lattice_mad_format(out_type: str, out_file_name: str, lat: _pybmad.LatStruct, ref_orbit: _pybmad.CoordStructAlloc1D | None = None, use_matrix_model: bool | None = None, include_apertures: bool | None = None, dr12_drift_max: float | None = None, ix_branch: int | None = None) -> bool:
    """
    Wrapper for Fortran routine write_lattice_mad_format

    Parameters
    ----------
    out_type : str
        Either 'MAD-8', or 'MAD-X'

    out_file_name : str
        Name of the mad output lattice file.

    lat : LatStruct
        Holds the lattice information.

    ref_orbit : 1D array of CoordStruct, optional
        Referece orbit for sad_mult and patch elements. This argument must be present if the lattice has sad_mult
        or patch elements and is being translated to MAD-8 or SAD.

    use_matrix_model : bool, optional
        Use a drift-matrix_drift model for wigglers/undulators? [A MAD "matrix" is a 2nd order Taylor map.] This
        switch is ignored for SAD conversion. Default is False -> Use a bend-drift-bend model. Note: sol_quad
        elements always use a drift-matrix-drift model.

    include_apertures : bool, optional
        If True (the default), add to the output lattice a zero length collimator element next to any non-
        collimator element that has an aperture. Note: MADX translations for non-drift elements can handle non-
        collimator elements with an aperture so in this case this argument is ignored.

    dr12_drift_max : float, optional
        Max deviation for drifts allowed before a correction matrix element is added. Default value is 1d-5. A
        negative number means use default.

    ix_branch : int, optional
        Index of lattice branch to use. Default = 0.

    Returns
    -------
    err : bool, optional
        Set True if, say a file could not be opened.
    """

class WriteLatticePalsFormat:
    """write_lattice_pals_format return type"""

    @property
    def pals_file(self) -> str: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_lattice_pals_format(lat: _pybmad.LatStruct) -> WriteLatticePalsFormat:
    """
    Wrapper for Fortran routine write_lattice_pals_format

    Parameters
    ----------
    lat : LatStruct
        Lattice

    Returns
    -------
    pals_file : str
        Pals lattice file name.

    err_flag : bool, optional
        Error flag
    """

def write_lattice_sad_format(out_file_name: str, lat: _pybmad.LatStruct, include_apertures: bool | None = None, ix_branch: int | None = None, err: bool | None = None) -> None:
    """
    Wrapper for Fortran routine write_lattice_sad_format

    Parameters
    ----------
    out_file_name : str

    lat : LatStruct

    include_apertures : bool, optional

    ix_branch : int, optional

    err : bool, optional
    """

class WriteLatticeScibmadFormat:
    """write_lattice_scibmad_format return type"""

    @property
    def scibmad_file(self) -> str: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_lattice_scibmad_format(lat: _pybmad.LatStruct) -> WriteLatticeScibmadFormat:
    """
    Wrapper for Fortran routine write_lattice_scibmad_format

    Parameters
    ----------
    lat : LatStruct
        Lattice

    Returns
    -------
    scibmad_file : str
        SciBmad lattice file name.

    err_flag : bool, optional
        Error flag
    """

def write_line_element(line: str, iu: int, ele: _pybmad.EleStruct, lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine write_line_element

    Parameters
    ----------
    line : str

    iu : int

    ele : EleStruct

    lat : LatStruct
    """

class WriteOpalFieldGridFile:
    """write_opal_field_grid_file return type"""

    @property
    def maxfield(self) -> float: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def write_opal_field_grid_file(opal_file_unit: int, ele: _pybmad.EleStruct, param: _pybmad.LatParamStruct) -> WriteOpalFieldGridFile:
    """
    Subroutine to write an OPAL lattice file using the information in
    a lat_struct. Optionally only part of the lattice can be generated.

    Parameters
    ----------
    opal_file_unit : int
        unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

    ele : EleStruct
        element to make map

    param : LatParamStruct
        Contains lattice information

    Returns
    -------
    maxfield : float
        absolute maximum found for element field scaling

    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_opal_lattice_file(opal_file_unit: int, lat: _pybmad.LatStruct) -> bool:
    """
    Subroutine to write an OPAL lattice file using the information in
    a lat_struct. Optionally only part of the lattice can be generated.

    Parameters
    ----------
    opal_file_unit : int
        unit number to write to

    lat : LatStruct
        Holds the lattice information.

    Returns
    -------
    err : bool, optional
        Set True if, say a file could not be opened.
    """

def write_time_particle_distribution(time_file_unit: int, bunch: _pybmad.BunchStruct, ele: _pybmad.EleStruct, style: str | None = None, branch: _pybmad.BranchStruct | None = None, format: str | None = None) -> bool:
    r"""
    Subroutine to write a time-based bunch from a standard Bmad bunch

    Note: 'BMAD' style (absolute curvilinear coordinates):
          n_particles_alive
          x/m  m*c^2 \beta_x*\gamma/eV y/m m*c^2\beta_y*\gamma/eV s/m m*c^2\beta_z*\gamma/eV time/s charge/C

          'OPAL' style (absolute curvilinear coordinates):
          n_particles_alive
          x/m  \beta_x*\gamma  y/m \beta_y*\gamma s/m \beta_s*\gamma

          'ASTRA' style (global Cartesian coordinates, first line is the reference particle used for z, pz, and t calculation):
          x/m y/m  z/m  m*c^2 \beta_x*\gamma/eV m*c^2 \beta_y*\gamma/eV m*c^2 \beta_z*\gamma/eV time/ns charge/nC species status

          'GPT' style (global Cartesian coordinates, with header labeling the columns)
          x/m y/m z/m \beta_x*\gamma \beta_y*\gamma \beta_z*\gamma t/s elementary_charge/C charge/elementary_charge

    Parameters
    ----------
    time_file_unit : int
        unit number to write to, if > 0

    bunch : BunchStruct
        bunch to be written. Particles are drifted to bmad_bunch.t_center for output

    ele : EleStruct
        Element being tracked through.

    style : str, optional
        Style of output file: 'BMAD' (default), 'OPAL', 'ASTRA', 'GPT'

    branch : BranchStruct, optional
        Required for 'ASTRA' style

    format : str, optional
        format for numerical output. default: 'es15.7'

    Returns
    -------
    err : bool, optional
        Set True if, say a file could not be opened.
    """

def xlafun(x: float, y: float, z: float) -> float:
    """
    Wrapper for Fortran routine xlafun

    Parameters
    ----------
    x : float

    y : float

    z : float

    Returns
    -------
    res : float
    """

def xraylib_nist_compound(name: str) -> int:
    """
    Routine to return the xraylib index for a given NIST compound.
    Taken from file xraylib/include/xraylib-nist_compounds.h

    Parameters
    ----------
    name : str
        Name of compound

    Returns
    -------
    indx : int
        Compound index. -1 if not found.
    """

def ylafun(x: float, y: float, z: float) -> float:
    """
    Wrapper for Fortran routine ylafun

    Parameters
    ----------
    x : float

    y : float

    z : float

    Returns
    -------
    res : float
    """

class ZAtSurface:
    """z_at_surface return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def dz_dxy(self) -> list[float]: ...

    @property
    def z(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def z_at_surface(ele: _pybmad.EleStruct, x: float, y: float, extend_grid: bool | None = None) -> ZAtSurface:
    """
    Routine return the height (z) of the surface for a particular (x,y) position.
    Remember: +z points into the element.

    Parameters
    ----------
    ele : EleStruct
        Element

    x : float
        Photon coordinates on surface.

    y : float
        Photon coordinates on surface.

    extend_grid : bool, optional
        If a grid is involved and (x, y) is outside of the grid, and extend_grid = True: Pretend (x, y) is at
        edge. Default is False.

    Returns
    -------
    err_flag : bool
        Set True if cannot compute z due to, say, point being outside of ellipseoid or grid bounds.

    z : float
        z coordinate.

    dz_dxy : 1D array of float (shape: 2), optional
        Surface slope at (x, y).
    """

def zero_ele_kicks(ele: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine zero_ele_kicks

    Parameters
    ----------
    ele : EleStruct
        Element with possible nonzero kicks.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with no kicks.
    """

def zero_ele_offsets(ele: _pybmad.EleStruct) -> None:
    """
    Wrapper for Fortran routine zero_ele_offsets

    Parameters
    ----------
    ele : EleStruct
        Element with possible nonzero offsets, etc.
        This parameter is an input/output and is modified in-place.
        As an output, ele: Element with no (mis)orientation.
    """

def zero_lr_wakes_in_lat(lat: _pybmad.LatStruct) -> None:
    """
    Routine to zero the long range wake amplitudes for the elements that have
    long range wakes in a lattice.

    Parameters
    ----------
    lat : LatStruct
        Lattice
    """

def zlafun(x: float, y: float, z: float) -> float:
    """
    Wrapper for Fortran routine zlafun

    Parameters
    ----------
    x : float

    y : float

    z : float

    Returns
    -------
    res : float
    """
