"""SimUtils routines"""

from collections.abc import Sequence
from typing import overload

import _pybmad


def all_pointer_to_string(a_ptr: _pybmad.AllPointerStruct, err: bool | None = None) -> str:
    """
    Wrapper for Fortran routine all_pointer_to_string

    Parameters
    ----------
    a_ptr : AllPointerStruct

    err : bool, optional

    Returns
    -------
    str : str
    """

def allocate_thread_states() -> None:
    """
    Routine to allocate random number state structures when openMP is used.
    """

def anomalous_moment_of(species: int) -> float:
    """
    Routine to return the anomolous moment for subatomic species type. Otherwise returns 0.

    Parameters
    ----------
    species : int
        Species ID.

    Returns
    -------
    moment : float
        Anomalous moment.
    """

def antiparticle(species: int) -> int:
    """
    Routine to return the antiparticle ID given the particle ID.
    For a molecule the anti-species is just the molecude with the charge reversed.

    Parameters
    ----------
    species : int
        Particle ID.

    Returns
    -------
    anti_species : int
        Antiparticle ID.
    """

def apfft(rdata_in: _pybmad.RealArray1D, bounds: Sequence[float], window: str, phase: float, diag: int | None = None) -> None:
    """
    Implements the All Phase FFT method for obtaining accurate phase from signal data.

    The signal data is truncated to an odd length, and the phase is relative to the central point.
    """

class ApfftCorr:
    """apfft_corr return type"""

    @property
    def phase(self) -> float: ...

    @property
    def amp(self) -> float: ...

    @property
    def freq(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def apfft_corr(rdata_in: _pybmad.RealArray1D, window: str, bounds: Sequence[float] | None = None, diag: int | None = None) -> ApfftCorr:
    """
    For real signal rdata_in, computes phase, frequency, and amplitude
    of peak found within bounds.  Algorithm is corrected all-phase FFT and should.

    This routine finds only one peak:  the largest amplitude within the bound.  Signals with multiple
    components can be investigated by varying bounds appropriately.

    Parameters
    ----------
    rdata_in : 1D array of float
        signal data.

    window : str
        'rec' or 'han' for rectangular or Hann window.

    bounds : 1D array of float (shape: 2), optional
        range within which to search for peak.

    diag : int, optional
        causes low-level routine apfft_ext to produce a fort.X file where X=9000+fid containing diag information.

    Returns
    -------
    phase : float
        phase of peak found in signal.

    amp : float
        amplitude of peak

    freq : float
        frequency of peak
    """

def apfft_ext(rdata: _pybmad.RealArray1D, bounds: Sequence[float], window: str, phase: float, amp: float, freq: float, diag: int | None = None) -> None:
    """
    Implements the All Phase FFT method for obtaining accurate phase from signal data.

    This "extended" apfft subroutine returns the amplitudes and frequency as well, for use
    by the corrected apfft subroutine in this module.
    """

def asinc(x: float, nd: int | None = None) -> float:
    """
    Wrapper for Fortran routine asinc

    Parameters
    ----------
    x : float

    nd : int, optional
        Derivative order. nd = 0 (default) -> compute arcsin(x) / x NOTE: Currently only nd = 0 and nd = 1 are
        implemented.

    Returns
    -------
    y : float
        nd^th derivative. of arcsin(x)/x
    """

def assert_equal(int_arr: _pybmad.IntArray1D, err_str: str) -> int:
    """
    Wrapper for Fortran routine assert_equal

    Parameters
    ----------
    int_arr : 1D array of int

    err_str : str

    Returns
    -------
    ival : int
    """

def atomic_number(species: int) -> int:
    """
    Routine to return the atomic number Z if species argument corresponds to an atomic particle  or is a proton.
    Set to the charge for atomic particles.
    Set to zero for molecules.

    Parameters
    ----------
    species : int
        Spicies ID.

    Returns
    -------
    atomic_num : int
        Atomic index. Set to zero if a molecule
    """

def atomic_species_id(charge: int, is_anti: bool, atomic_num: int, n_nuc: int) -> int:
    """
    Routine to return the species ID for an atom

    Parameters
    ----------
    charge : int
        Charge of the atom.

    is_anti : bool
        Is an anti-atom.

    atomic_num : int
        Atomic number.

    n_nuc : int
        Number of nucleons.

    Returns
    -------
    species_id : int
        Species ID number.
    """

def axis_angle_to_quat(axis: Sequence[float], angle: float) -> list[float]:
    """
    Routine to convert from axis + angle representation to a quaternion.

    Parameters
    ----------
    axis : 1D array of float (shape: 3)
        Axis of rotation.

    angle : float
        angle of rotation.

    Returns
    -------
    quat : 1D array of float (shape: 0:3)
        Rotation quaternion.
    """

def axis_angle_to_w_mat(axis: Sequence[float], angle: float) -> list[list[float]]:
    """
    Routine to construct the 3D rotation matrix w_mat given an axis of rotation
    and a rotation angle.

    Parameters
    ----------
    axis : 1D array of float (shape: 3)
        Rotation axis. Does not have to be normalized.

    angle : float
        Rotation angle in the range [-pi, pi].

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3)
        Rotation matrix
    """

class BicubicCmplxEval:
    """bicubic_cmplx_eval return type"""

    @property
    def df_dx(self) -> complex: ...

    @property
    def df_dy(self) -> complex: ...

    @property
    def f_val(self) -> complex: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bicubic_cmplx_eval(x_norm: float, y_norm: float, bi_coef: _pybmad.BicubicCmplxCoefStruct) -> BicubicCmplxEval:
    """
    Routine to evaluate a bicubic interpolating complex function.

    Use the routine bicubic_interpolation_cmplx_coefs to generate bi_coef.

    Note: In the equations below, the four points of the grid box being interpolated range
    from (x0, y0) to (x0+dx, y0+dy).

    Parameters
    ----------
    x_norm : float
        x_norm = (x - x0) / dx

    y_norm : float
        y_norm = (y - y0) / dy

    bi_coef : BicubicCmplxCoefStruct
        Coefficients.

    Returns
    -------
    f_val : complex
        Value of f.

    df_dx : complex, optional
        Normalized first derivative: True df/dx = df_dx * dx

    df_dy : complex, optional
        Normalized first derivative: True df/dy = df_dy * dy
    """

def bin_index(x: float, bin1_x_min: float, bin_delta: float) -> int:
    """
    Helper function to locate the appropriate histogram bin index.

    Parameters
    ----------
    x : float
        Input value to bin.

    bin1_x_min : float
        Minimum value of bin with index 1.

    bin_delta : float
        Bin width.

    Returns
    -------
    ix_bin : int
        Index of bin x is in.
    """

class BinXCenter:
    """bin_x_center return type"""

    @property
    def x_center(self) -> float: ...

    @property
    def ix_bin(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bin_x_center(ix_bin: int, bin1_x_min: float, bin_delta: float) -> BinXCenter:
    """
    Helper function to locate the center of a histogram bin.

    Parameters
    ----------
    ix_bin : int
        Index of bin under question.

    bin1_x_min : float
        Minimum value of bin with index 1.

    bin_delta : float
        Bin width.

    Returns
    -------
    ix_bin : int
        Index of bin under question.
    """

class BitSet:
    """bit_set return type"""

    @property
    def word(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bit_set(word: int, pos: int, set_to_1: bool) -> BitSet:
    """
    Routine to set a bit in a word.

    Parameters
    ----------
    word : int
        Input word
        This parameter is an input/output and is modified in-place.
        As an output, word: Word with bit set.

    pos : int
        position to set.

    set_to_1 : bool
        If True then bit is set to 1. If False bit is set to 0.

    Returns
    -------
    word : int
        Input word
        This parameter is an input/output and is modified in-place.
        As an output, word: Word with bit set.
    """

class BracketIndexForSpline:
    """bracket_index_for_spline return type"""

    @property
    def ix0(self) -> int: ...

    @property
    def ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bracket_index_for_spline(x_knot: _pybmad.RealArray1D, x: float, strict: bool | None = None, print_err: bool | None = None) -> BracketIndexForSpline:
    """
    Routine to find which interval to use for evaluating a spline.
    If strict = False (default), x is in range if
          x_knot(1) - (x_knot(2) - x_knot(1)) < x < x_knot(n) + (x_knot(n) - x_knot(n-1))
    If stric = True, x is in range if
          x_knot(1) <= x <= x_knot(n)
    where n = size(x_knot)

    Parameters
    ----------
    x_knot : 1D array of float
        Array of x values.

    x : float
        Evaluation point.

    strict : bool, optional
        Default is False. Determines acceptible range.

    print_err : bool, optional
        Default is True. Print error message if out of range?

    Returns
    -------
    ix0 : int
        If ok = True, x is in the interval [x_knot(ix0), x_knot(ix0+1)]

    ok : bool
        True if x is in range. False otherwise.
    """

def calc_file_number(file_name: str, num_in: int, num_out: int, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine calc_file_number

    Parameters
    ----------
    file_name : str

    num_in : int

    num_out : int

    err_flag : bool
    """

class Celbd:
    """celbd return type"""

    @property
    def elb(self) -> float: ...

    @property
    def eld(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def celbd(mc: float) -> Celbd:
    """
    Wrapper for Fortran routine celbd

    Parameters
    ----------
    mc : float

    Returns
    -------
    elb : float

    eld : float
    """

def cesr_getarg(i_arg: int) -> str:
    """
    Platform independent function to return the i'th command line argument.
    Use this with cesr_iargc.

    Note: The difference between this routine and the Fortran instrinsic
    get_command_argument is that for i_arg = 0, this routine returns the
    command line with the name of the executable removed from the beginning of the line.
    get_command_argument, on the other hand returns the name of the
    executable when the argument is 0.

    Parameters
    ----------
    i_arg : int
        Index of argument to return. i_arg = 0 => Entire line minus the executable string. i_arg = 1 => First
        argument.

    Returns
    -------
    arg : str
        i'th command line argument. If i_arg > number_of_args then arg is a blank string.
    """

def cesr_iargc() -> int:
    """
    Note: Use the Fortran intrinsic command_argument_count instead

    Platform independent function to return the number of command line arguments.
    Use this with cesr_getarg.
    """

def change_file_number(file_name: str, change: int) -> None:
    """
    Wrapper for Fortran routine change_file_number

    Parameters
    ----------
    file_name : str

    change : int
    """

def charge_of(species: int, default_: int | None = None) -> int:
    """
    Routine to return the charge, in units of e+, of a particle.

    Parameters
    ----------
    species : int
        Species ID.

    Returns
    -------
    charge : int
        particle charge.
    """

def charge_to_mass_of(species: int) -> float:
    """
    Routine to return the charge (in units of e+) to mass (in units of eV) ratio of a particle.

    Parameters
    ----------
    species : int
        Species ID.

    Returns
    -------
    charge_mass_ratio : float
        particle charge to mass ratio. (1/eV)
    """

def coarse_frequency_estimate(data: _pybmad.RealArray1D, error: bool | None = None) -> float:
    """
    Simple function to take periodic data and estimate
    the most dominant frequency by FFT.

    Parameters
    ----------
    data : 1D array of float
        data to analyze. Preferably size(data) is a power of 2 Otherwise the data is padded with zeros.

    Returns
    -------
    frequency : float
        Frequency corresponding to the largest FFT amplitude
    """

def complex_error_function(wr: float, wi: float, zr: float, zi: float) -> None:
    """
    Wrapper for Fortran routine complex_error_function

    Parameters
    ----------
    wr : float

    wi : float

    zr : float

    zi : float
    """

def cos_one(angle: float) -> float:
    """
    Wrapper for Fortran routine cos_one

    Parameters
    ----------
    angle : float
        Angle.

    Returns
    -------
    cos1 : float
        Result.
    """

def cosc(x: float, nd: int | None = None) -> float:
    """
    Wrapper for Fortran routine cosc

    Parameters
    ----------
    x : float

    nd : int, optional
        Derivative order. nd = 0 (default) -> compute (1 - cos(x)) / x^2 NOTE: Currently only nd = 0 and nd = 1
        are implemented.

    Returns
    -------
    y : float
        nd^th derivative of (1 - cos(x)) / x^2
    """

def create_a_spline(r0: _pybmad.RealArray1D, r1: _pybmad.RealArray1D, slope0: float, slope1: float) -> _pybmad.SplineStruct:
    """
    Routine to create a single spline given end point positions and slopes.
    The spline will pass through the data points and have the given slopes
    at these points.

    Modules used:
      use spline_mod

    Parameters
    ----------
    r0 : 1D array of float
        Start (x, y) point.

    r1 : 1D array of float
        End (x, y) point.

    slope0 : float
        Starting slope.

    slope1 : float
        End slope.

    Returns
    -------
    spline : SplineStruct
        Spline.
    """

def cross_product(a: _pybmad.RealArray1D, b: _pybmad.RealArray1D) -> list[float]:
    """
    Wrapper for Fortran routine cross_product

    Parameters
    ----------
    a : 1D array of float
        Input vectors.

    b : 1D array of float
        Input vectors.

    Returns
    -------
    c : 1D array of float (shape: 3)
        Cross product: a X b.
    """

def date_and_time_stamp(string: str, numeric_month: bool | None = None, include_zone: bool | None = None) -> None:
    """
    Wrapper for Fortran routine date_and_time_stamp

    Parameters
    ----------
    string : str

    numeric_month : bool, optional

    include_zone : bool, optional
    """

def destfixedwindowls(id: int) -> None:
    """
    Wrapper for Fortran routine destfixedwindowls

    Parameters
    ----------
    id : int
    """

def detab(str: str) -> None:
    """
    Wrapper for Fortran routine detab

    Parameters
    ----------
    str : str
    """

def display_size_and_resolution(ix_screen: int, x_size: float, y_size: float, x_res: float, y_res: float) -> None:
    """
    Wrapper for Fortran routine display_size_and_resolution

    Parameters
    ----------
    ix_screen : int

    x_size : float

    y_size : float

    x_res : float

    y_res : float
    """

def dj_bessel(m: int, arg: float) -> float:
    """
    Wrapper for Fortran routine dj_bessel

    Parameters
    ----------
    m : int

    arg : float
        Bessel argument.

    Returns
    -------
    dj_bes : float
        Bessel value.
    """

def djb_hash(str: str, old_hash: int | None = None) -> int:
    """
    Wrapper for Fortran routine djb_hash

    Parameters
    ----------
    str : str

    old_hash : int, optional

    Returns
    -------
    hash : int
    """

def djb_str_hash(in_str: str) -> str:
    """
    Wrapper for Fortran routine djb_str_hash

    Parameters
    ----------
    in_str : str

    Returns
    -------
    hash_str : str
    """

def downcase_string(string: str) -> None:
    """
    Wrapper for Fortran routine downcase_string

    Parameters
    ----------
    string : str
    """

class Elbd:
    """elbd return type"""

    @property
    def b(self) -> float: ...

    @property
    def d(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def elbd(phi: float, phic: float, mc: float) -> Elbd:
    """
    Wrapper for Fortran routine elbd

    Parameters
    ----------
    phi : float

    phic : float

    mc : float

    Returns
    -------
    b : float

    d : float
    """

class Elcbd:
    """elcbd return type"""

    @property
    def b(self) -> float: ...

    @property
    def dx(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def elcbd(c0: float, mc: float) -> Elcbd:
    """
    Wrapper for Fortran routine elcbd

    Parameters
    ----------
    c0 : float

    mc : float

    Returns
    -------
    b : float

    dx : float
    """

class Ellipinc:
    """ellipinc return type"""

    @property
    def ellipkinc(self) -> float: ...

    @property
    def ellipeinc(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ellipinc(phi: float, m: float) -> Ellipinc:
    """
    Calculates the first and second incomplete elliptic integrals,
    using methods from T. Fukushima, (2011, 2018)

    Uses classical transformations to handle negative m.
    This package needs a function for the third kind to use the new 2018 transformations.
    """

class Elsbd:
    """elsbd return type"""

    @property
    def b(self) -> float: ...

    @property
    def d(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def elsbd(s0: float, mc: float) -> Elsbd:
    """
    Wrapper for Fortran routine elsbd

    Parameters
    ----------
    s0 : float

    mc : float

    Returns
    -------
    b : float

    d : float
    """

def end_akima_spline_calc(spline: _pybmad.SplineStructArray1D, which_end: int) -> None:
    """
    Routine to calculate the slopes at the ends of a spline array

    Parameters
    ----------
    spline : 1D array of SplineStruct
        Array of splines.
        This parameter is an input/output and is modified in-place.
        As an output, spline: Array with slopes at end calculated.

    which_end : int
        0 => calculate slopes for the start end of the array. 1 => calculate slopes for the end end of the array.
    """

def err_exit(err_str: str | None = None) -> None:
    """
    Wrapper for Fortran routine err_exit

    Parameters
    ----------
    err_str : str, optional
    """

def factorial(n: int) -> float:
    """
    Wrapper for Fortran routine factorial

    Parameters
    ----------
    n : int
        Must be non-negative

    Returns
    -------
    fact : float
        n!. Will return negative number if there is an error.
    """

def faddeeva_function(z: Sequence[float], w: Sequence[float], dw: Sequence[Sequence[float]]) -> None:
    """
    Wrapper for Fortran routine faddeeva_function

    Parameters
    ----------
    z : 1D array of float (shape: 2)

    w : 1D array of float (shape: 2)

    dw : 2D array of float (shape: 2,2)
    """

def fft_1d(arr: _pybmad.ComplexArray1D, isign: int) -> None:
    """
    no longer exists
    subroutine fff_sub(line, error)
      implicit none
      character(*) line
      logical error
    end subroutine

    Parameters
    ----------
    arr : 1D array of complex
        Input array.
        This parameter is an input/output and is modified in-place.
        As an output, arr: FFT of array.

    isign : int
        -1 => "Forward" transform, +1 => "Backwards" transform.
    """

def file_directorizer(in_file: str, out_file: str, directory: str, add_switch: bool) -> None:
    """
    Wrapper for Fortran routine file_directorizer

    Parameters
    ----------
    in_file : str

    out_file : str

    directory : str

    add_switch : bool
    """

def file_get(string: str, dflt_file_name: str, file_name: str) -> None:
    """
    Wrapper for Fortran routine file_get

    Parameters
    ----------
    string : str

    dflt_file_name : str

    file_name : str
    """

def file_get_open(string: str, dflt_file_name: str, file_name: str, file_unit: int, readonly: bool) -> None:
    """
    Wrapper for Fortran routine file_get_open

    Parameters
    ----------
    string : str

    dflt_file_name : str

    file_name : str

    file_unit : int

    readonly : bool
    """

def file_suffixer(in_file_name: str, out_file_name: str, suffix: str, add_switch: bool) -> None:
    """
    Wrapper for Fortran routine file_suffixer

    Parameters
    ----------
    in_file_name : str

    out_file_name : str

    suffix : str

    add_switch : bool
    """

@overload
def find_location(arr: _pybmad.IntArray1D, value: int) -> int:
    """
    Wrapper for Fortran routine find_location_int

    Parameters
    ----------
    arr : 1D array of int

    value : int

    Returns
    -------
    ix_match : int
    """

@overload
def find_location(arr: _pybmad.BoolAlloc1D, value: bool) -> int:
    """
    Wrapper for Fortran routine find_location_logic

    Parameters
    ----------
    arr : 1D array of bool

    value : bool

    Returns
    -------
    ix_match : int
    """

@overload
def find_location(arr: _pybmad.RealArray1D, value: float) -> int:
    """
    Wrapper for Fortran routine find_location_real

    Parameters
    ----------
    arr : 1D array of float
        real(rp), logical, or integer

    value : float
        :).

    Returns
    -------
    ix_match : int
        Index of match. Zero if no match found.
    """

@overload
def find_location(arr: _pybmad.CharacterAlloc1D, value: str) -> int:
    """
    Wrapper for Fortran routine find_location_str

    Parameters
    ----------
    arr : 1D array of str

    value : str

    Returns
    -------
    ix_match : int
    """

def fine_frequency_estimate(data: _pybmad.RealArray1D) -> float:
    """
    Uses Laskar's method to accurately find the most dominant frequency
    A coarse estimate is first made by FFT.

    Parameters
    ----------
    data : 1D array of float
        data to analyze

    Returns
    -------
    frequency : float
        Frequency corresponding to the largest FFT amplitude
    """

def fixedwindowls(ynew: float, id: int) -> float:
    """
    Main function of the windowLS modult.  Each call to this function adds a data point to the fit
    and returns the derivative evaluated at the end of the window.  It is assumed that all data points
    are separeted by the same interval.
    This module is initialized with zeros for all data points, and so the results are unreliable until
    a number of data points equal to N has been entered.

    initFixedWindowLS must be called prior to calling this function.  destFixedWindowLS should be
    called when the instance is no longer needed.

    Parameters
    ----------
    """

class FourierAmplitude:
    """fourier_amplitude return type"""

    @property
    def cos_amp(self) -> float: ...

    @property
    def sin_amp(self) -> float: ...

    @property
    def dcos_amp(self) -> float: ...

    @property
    def dsin_amp(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def fourier_amplitude(data: _pybmad.RealArray1D, frequency: float) -> FourierAmplitude:
    """
    Computes cos_amp = (1/N) * sum_n=0^{N-1} data(n-1) cos(twopi*frequency*n)
        and  sin_amp = (1/N) * sum_n=0^{N-1} data(n-1) sin(twopi*frequency*n)
        and optionally dcos_amp = d/dfrequency cos_amp
                       dsin_amp = d/dfrequency sin_amp

    Parameters
    ----------
    data : 1D array of float
        data to analyze

    frequency : float
        frequency

    Returns
    -------
    cos_amp : float
        cosine amplitude

    sin_amp : float
        sine amplitude

    dcos_amp : float, optional
        cosine amplitude derivative

    dsin_amp : float, optional
        sine amplitude derivative
    """

class Gelbd:
    """gelbd return type"""

    @property
    def elb(self) -> float: ...

    @property
    def eld(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def gelbd(phi: float, mc: float) -> Gelbd:
    """
    Wrapper for Fortran routine gelbd

    Parameters
    ----------
    phi : float

    mc : float

    Returns
    -------
    elb : float

    eld : float
    """

def gen_complete_elliptic(kc: float, p: float, c: float, s: float, err_tol: float | None = None) -> float:
    """
    Wrapper for Fortran routine gen_complete_elliptic

    Parameters
    ----------
    kc : float
        Fuction input values.

    p : float
        Fuction input values.

    c : float
        Fuction input values.

    s : float
        Fuction input values.

    err_tol : float, optional
        Relative error tolerance. Default = 1d-12

    Returns
    -------
    value : float
        Output value.
    """

def get_a_char(wait: bool, ignore_this: _pybmad.CharacterAlloc1D | None = None) -> str:
    """
    Subroutine for getting a single character from the terminal.
    Also see: get_tty_char

    System Libraries that need to be linked to:
      readline curses

    Parameters
    ----------
    wait : bool
        If True then routine will wait until a keystroke has occured. If False and no keystroke is in the buffer
        then achar(0) will be returned as this_char.

    ignore_this : 1D array of str, optional
        List of characters to ignore. If a keystroke matches a character on this list the keystroke is ignored.

    Returns
    -------
    this_char : str
        Character returned
    """

def get_file_number(file_name: str, cnum_in: str, num_out: int, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine get_file_number

    Parameters
    ----------
    file_name : str

    cnum_in : str

    num_out : int

    err_flag : bool
    """

def get_file_time_stamp(file: str, time_stamp: str) -> None:
    """
    no longer exists
    subroutine get_next_number (filein, cnum, digits)
      implicit none
      character(*) filein
      character(*) cnum
      integer digits
    end subroutine
    """

def get_tty_char(wait: bool, flush: bool) -> str:
    """
    Subroutine for getting a single character from the terminal.
    Also see: get_a_char

    System Libraries that need to be linked to:
      readline curses

    Parameters
    ----------
    wait : bool
        If True then routine will wait until a keystroke has occured. If False and no keystroke is in the buffer
        then achar(0) will be returned as this_char.

    flush : bool
        If True then the keystroke buffer will be cleared first before any processing.

    Returns
    -------
    this_char : str
        Character returned
    """

def hanhan(N: int, hh: _pybmad.RealArray1D) -> None:
    """
    Wrapper for Fortran routine hanhan

    Parameters
    ----------
    N : int

    hh : 1D array of float
    """

def i_bessel(m: int, arg: float) -> float:
    """
    Wrapper for Fortran routine i_bessel

    Parameters
    ----------
    m : int
        Bessel order.

    arg : float
        Bessel argument.

    Returns
    -------
    i_bes : float
        Bessel value.
    """

def i_bessel_extended(m: int, arg: float) -> complex:
    """
    Wrapper for Fortran routine i_bessel_extended

    Parameters
    ----------
    m : int
        Bessel order.

    arg : float
        Bessel argument.

    Returns
    -------
    i_bes : complex
        Bessel value.
    """

def increment_file_number(file_name: str, digits: int, number: int, cnumber: str) -> None:
    """
    Wrapper for Fortran routine increment_file_number

    Parameters
    ----------
    file_name : str

    digits : int

    number : int

    cnumber : str
    """

def index_nocase(string1: str, string2: str) -> int:
    """
    Wrapper for Fortran routine index_nocase

    Parameters
    ----------
    string1 : str

    string2 : str

    Returns
    -------
    indx : int
    """

def initfixedwindowls(N: int, dt: float, order: int, der: int) -> int:
    """
    Initializes an instance of the fixed window least squares module.
    See module documentation (getf windowLS_mod) for use details.
    Any instance of windowLS created with this module should be destroyed with destFixedWindowLS.

    Parameters
    ----------
    N : int
        Number of data points to fit over. aka window size.

    dt : float
        Time interval between data points. It is assumed that the data is separated by fixed time intervals.

    order : int
        Order of fit polynomial.  Must be greater than or equal to der.

    der : int
        Order of derivative to be returned. Set der=0 to obtain the fit.
    """

def initial_lmdif() -> None:
    """Wrapper for Fortran routine initial_lmdif"""

def int_str(int_: int, width: int | None = None) -> str:
    """
    Wrapper for Fortran routine int_str

    Parameters
    ----------
    width : int, optional

    Returns
    -------
    str : str
    """

def interpolated_fft(cdata: _pybmad.ComplexArray1D, calc_ok: bool, opt_dump_spectrum: int | None = None, opt_dump_index: int | None = None) -> float:
    """
    Windows the complex data and used Numerical Recipes four1 to find the peak in the spectrum.
    The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
    available.
    """

def interpolated_fft_gsl(cdata: _pybmad.ComplexArray1D, calc_ok: bool, opt_dump_spectrum: int | None = None, opt_dump_index: int | None = None) -> float:
    """
    Windows the complex data and uses a mixed-radix GSL routine to find the peak in the spectrum.
    The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
    available.
    """

def is_alphabetic(string: str, valid_chars: str | None = None) -> bool:
    """
    no longer exists
    function inverse_prob (val) result (prob)
      import
      implicit none
      real(rp) prob
      real(rp) val
    end function
    """

def is_decreasing_sequence(array: _pybmad.RealArray1D, strict: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine is_decreasing_sequence

    Parameters
    ----------
    array : 1D array of float
        Sequence.

    strict : bool, optional
        If True (default) sequence must be strictly decreasing.

    Returns
    -------
    is_decreasing : bool
        Set True if sequence is decreasing.
    """

def is_false(param: float) -> bool:
    """
    Routine to translate from a real number to a boolian True or False.
    Translation: 0 = False, nonzero = True

    Also see: is_true and int_logic

    The typical use of this routine is for parameters in ele_struct%value(:) which
    is a real array. Some of the elements in the %value array are used to specify
    boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).

    Parameters
    ----------
    param : float
        Real number to be translated

    Returns
    -------
    this_false : bool
        Set True if param is zero. False otherwise.
    """

def is_increasing_sequence(array: _pybmad.RealArray1D, strict: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine is_increasing_sequence

    Parameters
    ----------
    array : 1D array of float
        Sequence.

    strict : bool, optional
        If True (default) sequence must be strictly increasing.

    Returns
    -------
    is_increasing : bool
        Set True if sequence is increasing.
    """

def is_integer(string: str, int_: int | None = None, delims: str | None = None, ix_word: int | None = None) -> bool:
    """
    Wrapper for Fortran routine is_integer

    Parameters
    ----------
    string : str

    delims : str, optional

    ix_word : int, optional

    Returns
    -------
    valid : bool
    """

def is_logical(string: str, ignore: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine is_logical

    Parameters
    ----------
    string : str

    ignore : bool, optional

    Returns
    -------
    valid : bool
    """

def is_real(string: str, ignore: bool | None = None, real_num: float | None = None) -> bool:
    """
    Wrapper for Fortran routine is_real

    Parameters
    ----------
    string : str

    ignore : bool, optional

    real_num : float, optional

    Returns
    -------
    valid : bool
    """

def is_subatomic_species(species: int) -> bool:
    """
    Routine to return True if species argument corresponds to a subatomic particle.

    Parameters
    ----------
    species : int
        Spicies ID.

    Returns
    -------
    is_subatomic : bool
        Set True if species corresponds to a subatomic particle.
    """

def is_true(param: float) -> bool:
    """
    Routine to translate from a real number to a boolian True or False.
    Translation: 0 = False, nonzero = True

    Also see: is_false and int_logic

    The typical use of this routine is for parameters in ele_struct%value(:) which
    is a real array. Some of the elements in the %value array are used to specify
    boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).

    Parameters
    ----------
    param : float
        Real number to be translated

    Returns
    -------
    this_true : bool
        Set False if param is zero. True otherwise.
    """

def j_bessel(m: int, arg: float) -> float:
    """
    Wrapper for Fortran routine j_bessel

    Parameters
    ----------
    m : int

    arg : float
        Bessel argument.

    Returns
    -------
    j_bes : float
        Bessel value.
    """

def linear_fit(x: _pybmad.RealArray1D, y: _pybmad.RealArray1D, n_data: int, a: float, b: float, sig_a: float, sig_b: float) -> None:
    """
    Wrapper for Fortran routine linear_fit

    Parameters
    ----------
    x : 1D array of float

    y : 1D array of float

    n_data : int

    a : float

    b : float

    sig_a : float

    sig_b : float
    """

def linear_fit_2d(x: _pybmad.RealArray1D, y: _pybmad.RealArray1D, z: _pybmad.RealArray1D) -> list[float]:
    """
    Wrapper for Fortran routine linear_fit_2d

    Parameters
    ----------
    x : 1D array of float
        Array of x-values.

    y : 1D array of float
        Array of y-values.

    z : 1D array of float
        Array of z-values

    Returns
    -------
    coef : 1D array of float (shape: 3)
        Coefficients of the linear fit
    """

def logic_str(logic: bool) -> str:
    """
    Wrapper for Fortran routine logic_str

    Parameters
    ----------
    logic : bool

    Returns
    -------
    str : str
    """

def lunget() -> int:
    """
    Wrapper for Fortran routine lunget

    Parameters
    ----------
    """

def make_legal_comment(comment_in: str, comment_out: str) -> None:
    """
    Wrapper for Fortran routine make_legal_comment

    Parameters
    ----------
    comment_in : str

    comment_out : str
    """

def mass_of(species: int) -> float:
    """
    Routine to return the mass, in units of eV/c^2, of a particle.
    To convert to AMU divide mass_of value by the constant atomic_mass_unit.

    Note: For atoms where the isotopic number is given, the mass is calculated using the neutral atomic mass
    adjusted by the weight of any added or missing electrons. The calculated mass is off very slightly due to
    binding energy effects. Exception: For #1H+ (proton) and #2H+ (deuteron) the exact mass is used since it is known.

    Parameters
    ----------
    species : int
        Species ID.

    Returns
    -------
    mass : float
        particle mass. Set to real_garbage$ if species value is invalid.
    """

def match_reg(str: str, pat: str) -> bool:
    """
    Wrapper for Fortran routine match_reg

    Parameters
    ----------
    str : str

    pat : str

    Returns
    -------
    is_match : bool
    """

def match_wild(string: str, template_: str) -> bool:
    """
    Wrapper for Fortran routine match_wild

    Parameters
    ----------
    string : str

    Returns
    -------
    is_match : bool
    """

def match_word(string: str, names: _pybmad.CharacterAlloc1D, ix: int, exact_case: bool | None = None, can_abbreviate: bool | None = None, matched_name: str | None = None) -> None:
    """
    Wrapper for Fortran routine match_word

    Parameters
    ----------
    string : str

    names : 1D array of str

    ix : int

    exact_case : bool, optional

    can_abbreviate : bool, optional

    matched_name : str, optional
    """

def maximize_projection(seed: float, cdata: _pybmad.ComplexArray1D) -> float:
    """
    Optimizer that uses Numerical Recipes brent to find a local maximum,
    which is the frequency that maximizes the projection.
    """

def milli_sleep(milli_sec: int) -> None:
    """
    Wrapper for Fortran routine milli_sleep

    Parameters
    ----------
    milli_sec : int
    """

def modulo2_dp(x: float, amp: float) -> float:
    """
    Function to return
        mod2 = x + 2 * n * amp
    where n is an integer chosen such that
       -amp <= mod2 < amp

    Parameters
    ----------
    x : float
        Real(sp), Real(rp), or Integer

    amp : float
        Must be positive.

    Returns
    -------
    mod2 : float
        Result
    """

def modulo2_int(x: int, amp: int) -> int:
    """
    Function to return
        mod2 = x + 2 * n * amp
    where n is an integer chosen such that
       -amp <= mod2 < amp

    Parameters
    ----------
    x : int
        Real(sp), Real(rp), or Integer

    amp : int
        Must be positive.

    Returns
    -------
    mod2 : int
        Result
    """

def modulo2_qp(x: float, amp: float) -> float:
    """
    Function to return
        mod2 = x + 2 * n * amp
    where n is an integer chosen such that
       -amp <= mod2 < amp

    Parameters
    ----------
    x : float
        Real(sp), Real(rp), or Integer

    amp : float
        Must be positive.

    Returns
    -------
    mod2 : float
        Result
    """

def modulo2_sp(x: float, amp: float) -> float:
    """
    Function to return
        mod2 = x + 2 * n * amp
    where n is an integer chosen such that
       -amp <= mod2 < amp

    Parameters
    ----------
    x : float
        Real(sp), Real(rp), or Integer

    amp : float
        Must be positive.

    Returns
    -------
    mod2 : float
        Result
    """

def n_bins_automatic(n_data: int) -> int:
    """Function to automatically select the number of bins"""

def n_choose_k(n: int, k: int) -> float:
    """
    Wrapper for Fortran routine n_choose_k

    Parameters
    ----------
    n : int
        Must be non-negative with n >= k.

    k : int
        Must be non-negative with n >= k.

    Returns
    -------
    nck : float
        N choose K will return negative number if there is an error.
    """

def n_spline_create(deriv0: _pybmad.RealArray1D, deriv1: _pybmad.RealArray1D, x1: float, n_spline: _pybmad.RealArray1D) -> None:
    """
    Wrapper for Fortran routine n_spline_create

    Parameters
    ----------
    deriv0 : 1D array of float
        Derivative vector from order 0 to some order n at x = 0.

    deriv1 : 1D array of float
        Derivative vector from order 0 to some order n at x = x1.

    x1 : float
        Location where deriv1 derivatives have been evaluated.

    n_spline : 1D array of float
        real(rp), Derivative vector from order 0 to order 2*n+1 of the interpolation spline.
    """

def naff(cdata: _pybmad.ComplexArray1D, freqs: _pybmad.RealArray1D, amps: _pybmad.ComplexArray1D, opt_dump_spectra: int | None = None, opt_zero_first: bool | None = None) -> None:
    """
    This subroutine implements the NAFF algorithm for calculating the spectra
    of periodic data.

    See naff_mod documentation for details.

    Frequencies returned are in units of 2pi. That is, freqs ranges from 0 to 1.

    freqs and amps must be allocated before hand.  This subroutine will repeat the
    decomposition loop until all elements of freqs and amps are populated.
    """

def nametable_add(nametable: _pybmad.NametableStruct, name: str, ix_name: int) -> None:
    """
    Wrapper for Fortran routine nametable_add

    Parameters
    ----------
    nametable : NametableStruct

    name : str

    ix_name : int
    """

def nametable_bracket_indexx(nametable: _pybmad.NametableStruct, name: str, n_match: int | None = None) -> int:
    """
    Wrapper for Fortran routine nametable_bracket_indexx

    Parameters
    ----------
    nametable : NametableStruct

    name : str

    n_match : int, optional

    Returns
    -------
    ix_max : int
    """

def nametable_change1(nametable: _pybmad.NametableStruct, name: str, ix_name: int) -> None:
    """
    Wrapper for Fortran routine nametable_change1

    Parameters
    ----------
    nametable : NametableStruct

    name : str

    ix_name : int
    """

def nametable_init(nametable: _pybmad.NametableStruct, n_min: int | None = None, n_max: int | None = None) -> None:
    """
    Wrapper for Fortran routine nametable_init

    Parameters
    ----------
    nametable : NametableStruct

    n_min : int, optional

    n_max : int, optional
    """

def nametable_remove(nametable: _pybmad.NametableStruct, ix_name: int) -> None:
    """
    Wrapper for Fortran routine nametable_remove

    Parameters
    ----------
    nametable : NametableStruct

    ix_name : int
    """

def negative_ampsquared(frequency: float, status: int | None = None) -> float:
    """
    Wrapper for Fortran routine negative_ampsquared

    Parameters
    ----------
    frequency : float

    status : int, optional

    Returns
    -------
    amp : float
    """

def negative_dampsquared(frequency: float, status: int | None = None) -> float:
    """
    Wrapper for Fortran routine negative_dampsquared

    Parameters
    ----------
    frequency : float

    status : int, optional

    Returns
    -------
    damp : float
    """

def omega_to_quat(omega: Sequence[float]) -> list[float]:
    """
    Routine to convert from omega + angle representation to a quaternion.

    Parameters
    ----------
    omega : 1D array of float (shape: 3)
        Axis of rotation + magnitude = rotation angle.

    Returns
    -------
    quat : 1D array of float (shape: 0:3)
        Rotation quaternion.
    """

def openpmd_species_name(species: int) -> str:
    """
    Routine to return the openPMD name of a particle species given the Bmad species ID.
    Note: the pmd_name does not include the particle charge. For example, if species
    corresponds to He+ then the pmd_name will be "He".

    Parameters
    ----------
    species : int
        Bmad species ID number.

    Returns
    -------
    pmd_name : str
        Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
    """

def ordinal_str(n: int) -> str:
    """
    Wrapper for Fortran routine ordinal_str

    Parameters
    ----------
    n : int

    Returns
    -------
    str : str
    """

def out_io_buffer_get_line(ix_line: int) -> str:
    """
    Routine to return the nuber of lines in the internal buffer.
    See the output_direct documentation for more details.
    """

def out_io_buffer_num_lines() -> int:
    """
    Routine to return the nuber of lines in the internal buffer.
    See the output_direct documentation for more details.
    """

def out_io_buffer_reset() -> None:
    """Routine to initialize the buffer used for capturing output."""

@overload
def out_io(level: int, routine_name: str, line: str, i_num: int, insert_tag_line: bool | None = None) -> None:
    """
    Wrapper for Fortran routine out_io_int

    Parameters
    ----------
    level : int

    routine_name : str

    line : str

    i_num : int

    insert_tag_line : bool, optional
    """

@overload
def out_io(level: int, routine_name: str, line1: str, line2: str | None = None, line3: str | None = None, line4: str | None = None, line5: str | None = None, line6: str | None = None, line7: str | None = None, line8: str | None = None, line9: str | None = None, line10: str | None = None, line11: str | None = None, line12: str | None = None, r_array: _pybmad.RealArray1D | None = None, i_array: _pybmad.IntArray1D | None = None, l_array: _pybmad.BoolAlloc1D | None = None, insert_tag_line: bool | None = None) -> None:
    """
    Wrapper for Fortran routine out_io_line12

    Parameters
    ----------
    level : int

    routine_name : str

    line1 : str

    line2 : str, optional

    line3 : str, optional

    line4 : str, optional

    line5 : str, optional

    line6 : str, optional

    line7 : str, optional

    line8 : str, optional

    line9 : str, optional

    line10 : str, optional

    line11 : str, optional

    line12 : str, optional

    r_array : 1D array of float, optional

    i_array : 1D array of int, optional

    l_array : 1D array of bool, optional

    insert_tag_line : bool, optional
    """

@overload
def out_io(level: int, routine_name: str, lines: _pybmad.CharacterAlloc1D, r_array: _pybmad.RealArray1D | None = None, i_array: _pybmad.IntArray1D | None = None, l_array: _pybmad.BoolAlloc1D | None = None, insert_tag_line: bool | None = None) -> None:
    """
    Wrapper for Fortran routine out_io_lines

    Parameters
    ----------
    level : int

    routine_name : str

    lines : 1D array of str

    r_array : 1D array of float, optional

    i_array : 1D array of int, optional

    l_array : 1D array of bool, optional

    insert_tag_line : bool, optional
    """

@overload
def out_io(level: int, routine_name: str, line: str, l_num: bool, insert_tag_line: bool | None = None) -> None:
    """
    Wrapper for Fortran routine out_io_logical

    Parameters
    ----------
    level : int

    routine_name : str

    line : str

    l_num : bool

    insert_tag_line : bool, optional
    """

@overload
def out_io(level: int, routine_name: str, line: str, r_num: float, insert_tag_line: bool | None = None) -> None:
    """
    Wrapper for Fortran routine out_io_real

    Parameters
    ----------
    level : int

    routine_name : str

    line : str

    r_num : float

    insert_tag_line : bool, optional
    """

def out_io_print_and_capture_setup(print_on: bool | None = None, capture_state: str | None = None, capture_add_null: bool | None = None) -> None:
    """
    Set whether a message from a call to out_io is sent to the terminal for printing and/or captured for program use.

    Capture may be desired, for example, to display the output in a separate window or captured output could be passed
    to a python process for processing.

    The procedure for how a message is handled is as follows:
      First: When out_io is called, the message level is used to determine if anything is to be printed or captured at all.
        When a program is started, everything will pass this test for printing and/or capturing.
        This behavior can be modified by calls to the output_direct routine.
      Second: If a message is to be printed and/or captured (passes the first step), then the internal print_on flag is used
        to determine if printing to the terminal and the internal capture_state flag is used to determine if capture is
        to be done. The initial setting of these flags is print_on = True and capture_state = 'OFF'.
        These internal flags can be set using the print_on and capture_state arguments of this routine.

    Notice that whether a message is also written to a file is independent of print and capture logic (see output_direct for more details).

    There are two capture modes. buffered (blocked) and unbuffered (unblocked) output.
    If a message is to be captured as outlined above, one and only one capture mode is used

    Unbuffered output is used when running multithreaded so that the program does not have to wait for output. For example, with a GUI.
    With unbuffered output, out_io calls three routines:
      out_io_called(level, routine_name)  ! Called at the start of a message.
      out_io_line(line)                   ! Called for each line of a message.
      out_io_end()                        ! Called at end of a message.
    The versions of these routines in the sim_utils library are just dummies.
    The idea is that modified versions of these routines can be used to capture the output.

    Buffered output uses an internal buffer to store the output.
    Output that has been buffered is retrieved by using the routines:
      out_io_buffer_reset
      out_io_buffer_num_lines and
      out_io_buffer_get_line

    Parameters
    ----------
    print_on : bool, optional
        If present, set the internal print_on flag to the value of this argument.

    capture_state : str, optional
        If present, set the internal capture_state to the value of this argument. Possible values:

    capture_add_null : bool, optional
        Is captured output null terminated (for interfacing with C/C++)?
    """

def parse_fortran_format(format_str: str, n_repeat: int, power: int, descrip: str, width: int, digits: int) -> None:
    """
    Wrapper for Fortran routine parse_fortran_format

    Parameters
    ----------
    format_str : str

    n_repeat : int

    power : int

    descrip : str

    width : int

    digits : int
    """

def pointer_to_locations(string: str, array: _pybmad.IntAlloc1D, num: int, ix_min: int, ix_max: int, names: _pybmad.CharacterAlloc1D | None = None, exact_case: bool | None = None, print_err: bool | None = None) -> None:
    """
    Wrapper for Fortran routine pointer_to_locations

    Parameters
    ----------
    string : str

    array : 1D array of int

    num : int

    ix_min : int

    ix_max : int

    names : 1D array of str, optional

    exact_case : bool, optional

    print_err : bool, optional
    """

def pointer_to_ran_state(ran_state: _pybmad.RandomStateStruct | None = None, ix_thread: int | None = None) -> _pybmad.RandomStateStruct | None:
    """
    Routine to point to the appropriate state structure for generating random numbers

    Parameters
    ----------
    ran_state : RandomStateStruct, optional
        Point to this if present. Otherwise point to the global saved state.

    ix_thread : int, optional
        Thread index.

    Returns
    -------
    ran_state_ptr : RandomStateStruct, optional
        Pointer to the appropriate state.
    """

def poly_eval(poly: _pybmad.RealArray1D, x: float, diff_coef: bool | None = None) -> float:
    """
    Wrapper for Fortran routine poly_eval

    Parameters
    ----------
    poly : 1D array of float
        Polynomial

    x : float
        Point to evaluate at.

    diff_coef : bool, optional
        poly(:) array are differentials? Default is False.

    Returns
    -------
    y : float
        Value of polynomial.
    """

def probability_funct(x: float) -> float:
    """
    Wrapper for Fortran routine probability_funct

    Parameters
    ----------
    x : float
        Function argument.

    Returns
    -------
    prob : float
    """

def projdd(a: _pybmad.ComplexArray1D, b: _pybmad.ComplexArray1D) -> complex:
    """
    Wrapper for Fortran routine projdd

    Parameters
    ----------
    a : 1D array of complex

    b : 1D array of complex
    """

def quadratic_roots(coefs: Sequence[float]) -> list[complex]:
    """
    Wrapper for Fortran routine quadratic_roots

    Parameters
    ----------
    coefs : 1D array of float (shape: 3)
        Coefficients of the quadratic equation with 0 = coefs(1) + coefs(2) * x + coefs(3) * x^2

    Returns
    -------
    root : 1D array of complex (shape: 2)
        Complex roots.
    """

@overload
def quat_conj(q_in: Sequence[complex]) -> list[complex]:
    """
    Overloaded name to create the conjugate of a quaternian.
    Overloaded functions are:
      Function quat_conj_real (q_in) result (q_out)
      Function quat_conj_complex (q_in) result (q_out)

    Parameters
    ----------
    q_in : 1D array of complex (shape: 0:3)
        Quaternion input.

    Returns
    -------
    q_out : 1D array of complex (shape: 0:3)
        Conjugate quaternion.
    """

@overload
def quat_conj(q_in: Sequence[float]) -> list[float]:
    """
    Overloaded name to create the conjugate of a quaternian.
    Overloaded functions are:
      Function quat_conj_real (q_in) result (q_out)
      Function quat_conj_complex (q_in) result (q_out)

    Parameters
    ----------
    q_in : 1D array of float (shape: 0:3)
        Quaternion input.

    Returns
    -------
    q_out : 1D array of float (shape: 0:3)
        Conjugate quaternion.
    """

def quat_inverse(q_in: Sequence[float]) -> list[float]:
    """
    Routine to create the inverse of a quaternian.

    Parameters
    ----------
    q_in : 1D array of float (shape: 0:3)
        Quaternion input.

    Returns
    -------
    q_out : 1D array of float (shape: 0:3)
        Inverse quaternion.
    """

@overload
def quat_mul(q1: Sequence[complex], q2: Sequence[complex], q3: Sequence[complex] | None = None, q4: Sequence[complex] | None = None, q5: Sequence[complex] | None = None, q6: Sequence[complex] | None = None, q7: Sequence[complex] | None = None, q8: Sequence[complex] | None = None, q9: Sequence[complex] | None = None) -> list[complex]:
    """
    Overloaded name to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
    Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
    Overloaded functions are:
      Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
      Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

    Parameters
    ----------
    q1 : 1D array of complex (shape: 0:3)
        Quaternions.

    q2 : 1D array of complex (shape: 0:3)
        Quaternions.

    q3 : 1D array of complex (shape: 0:3), optional
        More quaternions.

    q9 : 1D array of complex (shape: 0:3), optional
        More quaternions.

    Returns
    -------
    q_out : 1D array of complex (shape: 0:3)
        Resultant q1 * q2
    """

@overload
def quat_mul(q1: Sequence[float], q2: Sequence[float], q3: Sequence[float] | None = None, q4: Sequence[float] | None = None, q5: Sequence[float] | None = None, q6: Sequence[float] | None = None, q7: Sequence[float] | None = None, q8: Sequence[float] | None = None, q9: Sequence[float] | None = None) -> list[float]:
    """
    Overloaded name to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
    Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
    Overloaded functions are:
      Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
      Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

    Parameters
    ----------
    q1 : 1D array of float (shape: 0:3)
        Quaternions.

    q2 : 1D array of float (shape: 0:3)
        Quaternions.

    q3 : 1D array of float (shape: 0:3), optional
        More quaternions.

    q9 : 1D array of float (shape: 0:3), optional
        More quaternions.

    Returns
    -------
    q_out : 1D array of float (shape: 0:3)
        Resultant q1 * q2
    """

@overload
def quat_rotate(quat: Sequence[complex], vec_in: Sequence[complex]) -> list[complex]:
    """
    Overloaded name to rotate a vector using a quaternion..
    Overloaded functions are:
      Function quat_rotate_real (quat, vec_in) result (vec_out)
      Function quat_rotate_complex (quat, vec_in) result (vec_out)

    Parameters
    ----------
    quat : 1D array of complex (shape: 0:3)
        Quaternion to rotate with. Does not have to be normalized.

    vec_in : 1D array of complex (shape: 3)
        Initial vector.

    Returns
    -------
    vec_out : 1D array of complex (shape: 3)
        Final vector.
    """

@overload
def quat_rotate(quat: Sequence[float], vec_in: Sequence[float]) -> list[float]:
    """
    Overloaded name to rotate a vector using a quaternion..
    Overloaded functions are:
      Function quat_rotate_real (quat, vec_in) result (vec_out)
      Function quat_rotate_complex (quat, vec_in) result (vec_out)

    Parameters
    ----------
    quat : 1D array of float (shape: 0:3)
        Quaternion to rotate with. Does not have to be normalized.

    vec_in : 1D array of float (shape: 3)
        Initial vector.

    Returns
    -------
    vec_out : 1D array of float (shape: 3)
        Final vector.
    """

class QuatToAxisAngle:
    """quat_to_axis_angle return type"""

    @property
    def axis(self) -> list[float]: ...

    @property
    def angle(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def quat_to_axis_angle(quat: Sequence[float]) -> QuatToAxisAngle:
    """
    Routine to convert from quaternion to axis + angle representation.
    The angle will be in the range 0 <= angle <= pi.

    Parameters
    ----------
    quat : 1D array of float (shape: 0:3)
        Rotation quaternion. Assumed normalized.

    Returns
    -------
    axis : 1D array of float (shape: 3)
        Axis of rotation.

    angle : float
        angle of rotation in range [0, pi].
    """

def quat_to_omega(quat: Sequence[float]) -> list[float]:
    """
    Routine to convert rotation from quaternion representation to omega (axis + angle).

    Parameters
    ----------
    quat : 1D array of float (shape: 0:3)
        Rotation quaternion. Assumed normalized.

    Returns
    -------
    omega : 1D array of float (shape: 3)
        Axis of rotation + magnitude = rotation angle.
    """

def quat_to_w_mat(quat: Sequence[float]) -> list[list[float]]:
    """
    Routine to construct the 3D rotation matrix w_mat given a rotation quaternion

    Parameters
    ----------
    quat : 1D array of float (shape: 0:3)
        Quaternion.

    Returns
    -------
    w_mat : 2D array of float (shape: 3,3)
        Rotation matrix
    """

def query_string(query_str: str, upcase: bool, return_str: str, ix: int, ios: int) -> None:
    """
    Wrapper for Fortran routine query_string

    Parameters
    ----------
    query_str : str

    upcase : bool

    return_str : str

    ix : int

    ios : int
    """

def quote(str: str) -> str:
    """
    Wrapper for Fortran routine quote

    Parameters
    ----------
    str : str

    Returns
    -------
    q_str : str
    """

def quoten(str: _pybmad.CharacterAlloc1D, delim: str | None = None) -> str:
    """
    Wrapper for Fortran routine quoten

    Parameters
    ----------
    str : 1D array of str

    delim : str, optional

    Returns
    -------
    q_str : str
    """

def ran_default_state(set_state: _pybmad.RandomStateStruct | None = None) -> _pybmad.RandomStateStruct:
    """
    Routine to set or get the state of the default random number generator.
    See the ran_seed_put documentation for more details

    Parameters
    ----------
    set_state : RandomStateStruct, optional
        State to set the default generator to.

    Returns
    -------
    get_state : RandomStateStruct, optional
        Returns the state of the default generator.
    """

def ran_engine(set: str | None = None, get: str | None = None, ran_state: _pybmad.RandomStateStruct | None = None) -> None:
    r"""
    Routine to set what random number generator algorithm is used.
    If this routine is never called then pseudo_random$ is used.
    With sobseq quasi-random numbers the maximum dimension is 6.

    Parameters
    ----------
    set : str, optional
        Set the random number engine. Possibilities are: 'pseudo' -> Uses ran from Numerical Recipies (F90).
        'quasi'  -> Uses sobseq from Numerical Recipes. \'\'       -> Do nothing.

    get : str, optional
        Get the current (before any set) random number engine.

    ran_state : RandomStateStruct, optional
        Internal state. See the ran_seed_put documentation for more details.
    """

class RanGaussConverter:
    """ran_gauss_converter return type"""

    @property
    def get(self) -> str: ...

    @property
    def get_sigma_cut(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def ran_gauss_converter(set: str | None = None, set_sigma_cut: float | None = None, ran_state: _pybmad.RandomStateStruct | None = None) -> RanGaussConverter:
    r"""
    Routine to set what conversion routine is used for converting
    uniformly distributed random numbers to Gaussian distributed random numbers.

    If this routine is not called then exact_gaussian$ is used.

    exact_gaussian$ is a straight forward converter as explained in Numerical recipes.

    quick_gaussian$ is a quick a dirty approximation with a cutoff so that no
    numbers will be generated beyound what is set for sigma_cut.

    A negative sigma_cut means that the exact_gaussian$ will not be limited
    and the quick_gaussian$ will use a default of 10.0

    Note: Because of technical issues, when using the quasi_random$ number generator
    (see the ran_engine routine), the quick_gaussian$ method will automatically be
    used independent of what was set with this routine.

    Parameters
    ----------
    set : str, optional
        Set the random number engine. Possibilities are: 'exact' 'quick'  ! Old deprecated: 'limited' 'ziggurat'
        \'\'       ! Do nothing

    set_sigma_cut : float, optional
        Sigma cutoff. Initially: sigma_cut = -1.

    ran_state : RandomStateStruct, optional
        Internal state. See the ran_seed_put documentation for more details.

    Returns
    -------
    get : str, optional
        Get the current (before any set) gaussian converter.

    get_sigma_cut : float, optional
        Get the current (before any set) sigma cutoff.
    """

def ran_gauss_scalar(ran_state: _pybmad.RandomStateStruct | None = None, sigma_cut: float | None = None, index_quasi: int | None = None) -> float:
    """
    Routine to return a gaussian distributed random number with unit sigma.
    This routine uses the same algorithm as gasdev from Numerical Recipes.

    Note: ran_gauss is an overloaded name for:
        ran_gauss_scalar   ! harvest is a scalar
        ran_gauss_vector   ! harvest is a 1-D array.

    Note: Use ran_seed_put for initialization.
    Note: Use ran_engine to set which random number generator to use.
    Note: Use ran_gauss_converter to set which conversion routine to use.

    Parameters
    ----------
    ran_state : RandomStateStruct, optional
        Internal state. See the ran_seed_put documentation for more details.

    sigma_cut : float, optional
        If present and positive will override setting of ran_state.gauss_sigma_cut.

    Returns
    -------
    harvest : float
        Random number. Or
        As an output, harvest: Random number array.
    """

def ran_gauss_vector(harvest: _pybmad.RealArray1D, ran_state: _pybmad.RandomStateStruct | None = None, sigma_cut: float | None = None) -> None:
    """
    Routine to return a gaussian distributed random number with unit sigma.
    This routine uses the same algorithm as gasdev from Numerical Recipes.

    Note: ran_gauss is an overloaded name for:
        ran_gauss_scalar   ! harvest is a scalar
        ran_gauss_vector   ! harvest is a 1-D array.

    Note: Use ran_seed_put for initialization.
    Note: Use ran_engine to set which random number generator to use.
    Note: Use ran_gauss_converter to set which conversion routine to use.

    Parameters
    ----------
    harvest : 1D array of float
        Random number. Or
        As an output, harvest: Random number array.

    ran_state : RandomStateStruct, optional
        Internal state. See the ran_seed_put documentation for more details.

    sigma_cut : float, optional
        If present and positive will override setting of ran_state.gauss_sigma_cut.
    """

def ran_seed_get() -> int:
    """
    Routine to return the seed used for the random number generator.

    Parameters
    ----------

    Returns
    -------
    seed : int
        Random number seed used.
    """

def ran_seed_put(seed: int, mpi_offset: int | None = None) -> None:
    """
    Routine to seed a random number generator.

    If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
    the system clock will be used to generate the seed.

    Note: The seed is only used with the pseudo_random$ engine.
    Note: Use the subroutine ran_seed_get(seed) to get the seed used.
    Note: Use pointer_to_ran_state() to access the ran state directly.

    Parameters
    ----------
    seed : int
        Seed number. If seed = 0 then a seed will be choosen based upon the system clock.

    mpi_offset : int, optional
        Offset added to seed. Default is zero. Used with MPI processes ensure different threads use different
        random numbers.
    """

@overload
def ran_uniform(ran_state: _pybmad.RandomStateStruct | None = None, index_quasi: int | None = None) -> float:
    """
    Routine to return a random number uniformly distributed in the
    interval [0, 1]. This routine uses the same algorithm as ran or sobseq
    from Numberical Recipes in Fortran90.
    See ran_engine.

    Note: ran_uniform is an overloaded name for:
        ran_uniform_scalar   ! harvest is a scalar
        ran_uniform_vector   ! harvest is a 1-D array.

    Note: Use ran_seed_put for initialization.
    Note: Use ran_engine to set which random number generator to use.

    Parameters
    ----------
    ran_state : RandomStateStruct, optional
        Internal state. See the ran_seed_put documentation for more details.

    Returns
    -------
    harvest : float
        Random number. Or
        As an output, harvest: Random number array.
    """

@overload
def ran_uniform(harvest: _pybmad.RealArray1D, ran_state: _pybmad.RandomStateStruct | None = None) -> None:
    """
    Routine to return a random number uniformly distributed in the
    interval [0, 1]. This routine uses the same algorithm as ran or sobseq
    from Numberical Recipes in Fortran90.
    See ran_engine.

    Note: ran_uniform is an overloaded name for:
        ran_uniform_scalar   ! harvest is a scalar
        ran_uniform_vector   ! harvest is a 1-D array.

    Note: Use ran_seed_put for initialization.
    Note: Use ran_engine to set which random number generator to use.

    Parameters
    ----------
    harvest : 1D array of float
        Random number. Or
        As an output, harvest: Random number array.

    ran_state : RandomStateStruct, optional
        Internal state. See the ran_seed_put documentation for more details.
    """

def rcelbd(mc: float, elb: float, eld: float) -> None:
    """
    Wrapper for Fortran routine rcelbd

    Parameters
    ----------
    mc : float

    elb : float

    eld : float
    """

def read_a_line(prompt: str, trim_prompt: bool | None = None, prompt_color: str | None = None, prompt_bold: bool | None = None, history_file: str | None = None) -> str:
    """
    Subroutine to read a line of input from the terminal.
    The line is also add to the history buffer so that the up-arrow
    and down-arrow keys can be used to recall past commands.

    Also see:
      readline_read_history
      readline_write_history

    System Libraries that need to be linked to:
      readline curses

    Parameters
    ----------
    prompt : str
        Prompt string to use.

    trim_prompt : bool, optional
        If present and True then trim the prompt string and add a single blank before printing the prompt string.
        Default is True.

    prompt_color : str, optional
        Color of the prompt. Possibilities are: 'BLACK', 'RED', 'GREEN', 'YELLOW', 'BLUE', 'MAGENTA', 'CYAN',
        'GRAY', 'DEFAULT'. The 'DEFAULT' setting (the default) does not set the prompt color.

    prompt_bold : bool, optional
        If present and True then the prompt will be printed in bold.

    history_file : str, optional
        If present, add line_out to a file whose name is given by history_file. History files are useful for
        saving the command history in between when a program is run multiple times.

    Returns
    -------
    line_out : str
        Line typed by the user. Note: If cntl-D is pressed, line_out = achar(24).
    """

def readline_read_history(history_file: str) -> int:
    """
    Routine to add the contents of a file to the readline history list.
    Use this routine with the read_a_line routine.

    Parameters
    ----------
    history_file : str
        Name of the history file. EG: '~/.my_history'

    Returns
    -------
    status : int
        0 = Success, otherwise failure.
    """

def readline_write_history(history_file: str) -> int:
    """
    Routine to write the contents of the readline history list to a file.
    Use this routine with the read_a_line routine.

    Parameters
    ----------
    history_file : str
        Name of the history file. EG: '~/.my_history'

    Returns
    -------
    status : int
        0 = Success, otherwise failure.
    """

def real_num_fortran_format(number: float, width: int, n_blanks: int | None = None) -> str:
    """
    Wrapper for Fortran routine real_num_fortran_format

    Parameters
    ----------
    number : float

    width : int

    n_blanks : int, optional

    Returns
    -------
    fmt_str : str
    """

def real_path(path_in: str, path_out: str) -> bool:
    """
    Wrapper for Fortran routine real_path

    Parameters
    ----------
    path_in : str

    path_out : str

    Returns
    -------
    is_ok : bool
    """

def real_str(r_num: float, n_signif: int | None = None, n_decimal: int | None = None) -> str:
    """
    Wrapper for Fortran routine real_str

    Parameters
    ----------
    r_num : float

    n_signif : int, optional

    n_decimal : int, optional

    Returns
    -------
    str : str
    """

def real_to_string(real_num: float, width: int, n_signif: int | None = None, n_decimal: int | None = None) -> str:
    """
    Wrapper for Fortran routine real_to_string

    Parameters
    ----------
    real_num : float

    width : int

    n_signif : int, optional

    n_decimal : int, optional

    Returns
    -------
    str : str
    """

def reallocate_spline(spline: _pybmad.SplineStructAlloc1D, n: int, n_min: int | None = None, exact: bool | None = None) -> None:
    """
    Subroutine to allocate an allocatable spline_struct array.
    The data of the array is preserved but data at the end of the
    array will be lost if n is less than the original size of the array

    Parameters
    ----------
    spline : 1D array of SplineStruct
        Spline to reallocate.
        This parameter is an input/output and is modified in-place.
        As an output, spline: Allocated spline.

    n : int
        Upper bound needed for 1-dimensional arrays.

    n_min : int, optional
        Lower bound of spline array. Default is 1.

    exact : bool, optional
        If present and False then the size of the output array is permitted to be larger than n. Default is True.
    """

def relbd(phi: float, phic: float, mc: float, b: float, d: float) -> None:
    """
    Wrapper for Fortran routine relbd

    Parameters
    ----------
    phi : float

    phic : float

    mc : float

    b : float

    d : float
    """

def relcbd(c0: float, mc: float, b: float, dx: float) -> None:
    """
    Wrapper for Fortran routine relcbd

    Parameters
    ----------
    c0 : float

    mc : float

    b : float

    dx : float
    """

def relsbd(s0: float, mc: float, b: float, d: float) -> None:
    """
    Wrapper for Fortran routine relsbd

    Parameters
    ----------
    s0 : float

    mc : float

    b : float

    d : float
    """

def rgelbd(phi: float, mc: float, elb: float, eld: float) -> None:
    """
    Wrapper for Fortran routine rgelbd

    Parameters
    ----------
    phi : float

    mc : float

    elb : float

    eld : float
    """

class RmsValue:
    """rms_value return type"""

    @property
    def ave_val(self) -> float: ...

    @property
    def rms_val(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def rms_value(val_arr: _pybmad.RealArray1D, good_val: _pybmad.BoolAlloc1D | None = None) -> RmsValue:
    """
    Wrapper for Fortran routine rms_value

    Parameters
    ----------
    val_arr : 1D array of float
        Array of reals.

    good_val : 1D array of bool, optional
        If present, only calculate RMS where good_val(i) = True.

    Returns
    -------
    rms_val : float
        RMS value. Set to real_garbage$ if there is a problem.

    ave_val : float, optional
        average value.
    """

def rot_2d(vec_in: Sequence[float], angle: float) -> list[float]:
    """
    Wrapper for Fortran routine rot_2d

    Parameters
    ----------
    vec_in : 1D array of float (shape: 2)
        Init vec

    angle : float
        angle in radians.

    Returns
    -------
    vec_out : 1D array of float (shape: 2)
        Rotated vec.
    """

def rotate_vec(vec: _pybmad.RealArray1D, axis: int, angle: float) -> None:
    """
    Basic routine to rotate vector components around the x, y, or z axis.

    Parameters
    ----------
    vec : 1D array of float
        vector
        This parameter is an input/output and is modified in-place.
        As an output, vec: Rotated vector.

    axis : int
        x_axis$, y_axis$, or z_axis$

    angle : float
        angle to rotate.
    """

def rotate_vec_given_axis_angle(vec_in: Sequence[float], axis: _pybmad.RealArray1D, angle: float) -> list[float]:
    """
    Routine to rotate a vector.

    Parameters
    ----------
    vec_in : 1D array of float (shape: 3)
        Initial vector.

    axis : 1D array of float
        Axis of rotation. Must be normalized to 1.

    angle : float
        Angle to rotate by

    Returns
    -------
    vec_out : 1D array of float (shape: 3)
        Final vector.
    """

def rp8(int_in: int) -> float:
    """
    Routine to convert from integer to real of type rp.
    This routine is used to avoid the implicit integer to single precision that happens when
    multiplying int*real(rp).

    Parameters
    ----------
    int_in : int
        Input integer.

    Returns
    -------
    re_out : float
        Equiv real.
    """

def rserbd(y: float, m: float, b: float, d: float) -> None:
    """
    Wrapper for Fortran routine rserbd

    Parameters
    ----------
    y : float

    m : float

    b : float

    d : float
    """

def run_timer(command: str, time: float | None = None, time0: float | None = None) -> None:
    """
    Wrapper for Fortran routine run_timer

    Parameters
    ----------
    command : str

    time : float, optional

    time0 : float, optional
    """

class Serbd:
    """serbd return type"""

    @property
    def b(self) -> float: ...

    @property
    def d(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def serbd(y: float, m: float) -> Serbd:
    """
    Wrapper for Fortran routine serbd

    Parameters
    ----------
    y : float

    m : float

    Returns
    -------
    b : float

    d : float
    """

def set_all_ptr(a_ptr: _pybmad.AllPointerStruct, value: float, delta: bool | None = None, value_set: float | None = None) -> None:
    """
    Wrapper for Fortran routine set_all_ptr

    Parameters
    ----------
    a_ptr : AllPointerStruct

    value : float

    delta : bool, optional

    value_set : float, optional
    """

def set_env(env_name: str, env_value: str, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine set_env

    Parameters
    ----------
    env_name : str

    env_value : str

    err_flag : bool
    """

@overload
def set_parameter(param_val: int, set_val: int) -> int:
    """
    Wrapper for Fortran routine set_parameter_int

    Parameters
    ----------
    param_val : int

    set_val : int

    Returns
    -------
    save_val : int
    """

@overload
def set_parameter(param_val: bool, set_val: bool) -> bool:
    """
    Wrapper for Fortran routine set_parameter_logic

    Parameters
    ----------
    param_val : bool

    set_val : bool

    Returns
    -------
    save_val : bool
    """

@overload
def set_parameter(param_val: float, set_val: float) -> float:
    """
    Wrapper for Fortran routine set_parameter_real

    Parameters
    ----------
    param_val : float

    set_val : float

    Returns
    -------
    save_val : float
    """

def set_species_charge(species_in: int, charge: int) -> int:
    """
    Routine to return the ID for a particle of the same type as species_in but with a different charge.
    Exception: If species_in corresponds to a subatomic particle, the charge argument is ignored and
    species_charged will be set equal to species_in.

    Parameters
    ----------
    species_in : int
        Input species.

    charge : int
        Charge to set species_charged to.

    Returns
    -------
    species_charged : int
        Species of the same type as species_in but with different charge.
    """

@overload
def sign_of(num: int, zero_is_zero: bool | None = None) -> int:
    """
    Routine to return the sign of a number.
    Note: Fortran instrinsic sign function is similar to sign_of with zero_is_zero = False.

    Parameters
    ----------
    num : int
        Input number

    zero_is_zero : bool, optional
        If True (default), num = 0 gives num_sign = 0. If False, num = 0 gives num_sign = 1.

    Returns
    -------
    num_sign : int
        +1 if num is positive, -1 if num is negative, and 0 or +1 if num is zero depending upon setting of
        zero_is_zero.
    """

@overload
def sign_of(num: float, zero_is_zero: bool | None = None) -> int:
    """
    Routine to return the sign of a number.
    Note: Fortran instrinsic sign function is similar to sign_of with zero_is_zero = False.

    Parameters
    ----------
    num : float
        Input number

    zero_is_zero : bool, optional
        If True (default), num = 0 gives num_sign = 0. If False, num = 0 gives num_sign = 1.

    Returns
    -------
    num_sign : int
        +1 if num is positive, -1 if num is negative, and 0 or +1 if num is zero depending upon setting of
        zero_is_zero.
    """

def sinc(x: float, nd: int | None = None) -> float:
    """
    Wrapper for Fortran routine sinc

    Parameters
    ----------
    x : float
        Number.

    nd : int, optional
        Derivative order. nd = 0 (default) -> compute sin(x) / x NOTE: Currently only nd = 0 and nd = 1 are
        implemented.

    Returns
    -------
    y : float
        nd^th derivative of sin(x) / x
    """

def sincc(x: float, nd: int | None = None) -> float:
    """
    Wrapper for Fortran routine sincc

    Parameters
    ----------
    x : float
        Number.

    nd : int, optional
        Derivative order. nd = 0 (default) -> compute (x - sin(x)) / x^3 NOTE: Currently only nd = 0 and nd = 1
        are implemented.

    Returns
    -------
    y : float
        nd^th derivative of (x - sin(x)) / x^3
    """

def sinhx_x(x: float, nd: int | None = None) -> float:
    """
    Wrapper for Fortran routine sinhx_x

    Parameters
    ----------
    x : float
        Number.

    nd : int, optional
        Derivative order. nd = 0 (default) -> compute sinh(x) / x NOTE: Currently only nd = 0 and nd = 1 are
        implemented.

    Returns
    -------
    y : float
        nd^th derivative of sinh(x) / x.
    """

def skip_header(ix_unit: int, error_flag: bool) -> None:
    """
    Wrapper for Fortran routine skip_header

    Parameters
    ----------
    ix_unit : int

    error_flag : bool
    """

def special_projection(f: float, status: int | None = None) -> float:
    """
    Calculates <cdata | exp(i theta)>

    Used only by maximize projection. Uses data global to the function to accomodate stock NR routine.
    """

def species_id(name: str, default_: int | None = None, print_err: bool | None = None) -> int:
    """
    Routine to return the integer ID index of a particle species given the name.

    For subatomic particles, the case does not matter.
    For all other types of particles, the case does matter.

    Parameters
    ----------
    name : str
        Name of the species.

    print_err : bool, optional
        Print error message? Default is True. If False, return species = invalid$,

    Returns
    -------
    species : int
        Species ID. Will return invalid$ if name is not valid. Will return not_set$ if name is blank
    """

def species_id_from_openpmd(pmd_name: str, charge: int) -> int:
    """
    Routine to return the Bmad species ID given the openPMD species name and given particle charge.
    Note: If pmd_name corresponds to a subatomic particle, the charge argument is ignored.

    Parameters
    ----------
    pmd_name : str
        OpenPMD species name.

    charge : int
        Species charge. Ignored for subatomic particles.

    Returns
    -------
    species : int
        Bmad spicies ID number.
    """

def species_name(species: int) -> str:
    """
    Routine to return the name of a particle species given the integer index.

    Parameters
    ----------
    species : int
        Species ID.

    Returns
    -------
    name : str
        Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
    """

def species_of(mass: float, charge: int) -> int:
    """
    Routine to return the integer ID index of a particle species given the mass and charge.
    Note: Currently this routine only works for subatomic particles and is used for decoding PTC flat files.

    Parameters
    ----------
    mass : float
        Mass of the particle

    charge : int
        Charge of the particle.

    Returns
    -------
    species : int
        Species ID. Will return invalid$ if name is not valid.
    """

def spin_of(species: int, non_subatomic_default: float | None = None) -> float:
    """
    Routine to return the spin, in units of hbar, of a particle.
    This routine is only valid for subatomic particles.
    For all other particles, the returned spin value will be the value of non_subatomic_default.

    Parameters
    ----------
    species : int
        Species ID.

    non_subatomic_default : float, optional
        Default value to be used for non-subatomic species. Default value of this argument is zero.

    Returns
    -------
    spin : float
        Particle spin.
    """

def spline1(a_spline: _pybmad.SplineStruct, x: float, n: int | None = None) -> float:
    """
    Function for spline evaluation using a single spline (instead of a spline array).
    Also see:
      spline_evaluate
      spline_akima_interpolate

    Modules used:
      use spline_mod

    Parameters
    ----------
    a_spline : SplineStruct
        Single spline structure.

    x : float
        Point for evaluation.

    n : int, optional
        Output derivative order. May be -1, 0, 1, 2, or 3. Default is 0. n = -1 => output is integral of y from
        a_spline.x0 to x. n = 1 => output is dy/dx, n = 2 => output is d^2y/dx^2, etc.

    Returns
    -------
    y : float
        Interpolated spline value or derivative.
    """

def spline_akima(spline: _pybmad.SplineStructArray1D) -> bool:
    """
    Given a set of (x,y) points we want to interpolate between the points.
    This subroutine computes the semi-hermite cubic spline developed by
    Hiroshi Akima. The spline goes thorugh all the points (that is, it is
    not a smoothing spline). For interpolation use:
      spline_evaluate
      spline_akima_interpolate ! You do not need to call spline_akima if you use this routine.

    Reference:
      H Akima, "A New Method of Interpolation and Smooth Curve Fitting Based
      on Local Procedures", J. Assoc. Comp. Mach., Vol 17(4), 589-602 (1970).

    Modules used:
      use spline_mod

    Parameters
    ----------
    spline : 1D array of SplineStruct

    Returns
    -------
    ok : bool
        Set .false. if something is wrong (like less than 2 points used).
    """

class SplineAkimaInterpolate:
    """spline_akima_interpolate return type"""

    @property
    def ok(self) -> bool: ...

    @property
    def y(self) -> float: ...

    @property
    def dy(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spline_akima_interpolate(x_knot: _pybmad.RealArray1D, y_knot: _pybmad.RealArray1D, x: float) -> SplineAkimaInterpolate:
    """
    Routine to interpolate using an akima spline.

    When evaluating at enough points, this routine is slower than calling spline_akima to
    first evaluate the spline coefficients and then repeatedly calling spline_evaluate.

    The advantage of this routine is that only the (x, y) knot points need to be stored
    and it will be faster if the number of evaluations is small.

    This routine will extrapolate past the range of x_knot(:) up to a distance equal to the
    length between an end point and the point just inside the end point.

    Parameters
    ----------
    x_knot : 1D array of float
        Array of x values for the knot points. Must have more than 2 points and be in asending order.

    y_knot : 1D array of float
        Array of y values for the knot points. Must be same size as x_knot(:).

    x : float
        Point to evaluate at.

    Returns
    -------
    ok : bool
        Set .true. if everything ok, That is, x is within the spline range.

    y : float, optional
        Spline interpolation.

    dy : float, optional
        Spline derivative interpolation.
    """

class SplineEvaluate:
    """spline_evaluate return type"""

    @property
    def ok(self) -> bool: ...

    @property
    def y(self) -> float: ...

    @property
    def dy(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def spline_evaluate(spline: _pybmad.SplineStructArray1D, x: float) -> SplineEvaluate:
    """
    Subroutine to evalueate a spline at a set of points.

    A point outside of the range of knot points is an error.
    Also see:
      spline1
      spline_akima_interpolate

    A spline may be generated using, for example, the spline_akima routine.

    Modules used:
      use spline_mod

    Parameters
    ----------
    spline : 1D array of SplineStruct
        Spline structure.

    x : float
        point for evaluation.

    Returns
    -------
    ok : bool
        Set .true. if everything ok. That is, x is within the spline range.

    y : float, optional
        Spline interpolation.

    dy : float, optional
        Spline derivative interpolation.
    """

def sqrt_alpha(alpha: float, x: float) -> float:
    """
    Wrapper for Fortran routine sqrt_alpha

    Parameters
    ----------
    alpha : float
        Number

    x : float
        Number

    Returns
    -------
    y : float
        Result.
    """

def sqrt_one(x: float, nd: int | None = None) -> float:
    """
    Wrapper for Fortran routine sqrt_one

    Parameters
    ----------
    x : float
        Number

    nd : int, optional
        Derivative order. nd = 0 (default) -> compute Sqrt[1+x] - 1. NOTE: Currently only nd = 0 and nd = 1 are
        implemented.

    Returns
    -------
    ds1 : float
    """

def str_count(str: str, match_: str) -> int:
    """
    Wrapper for Fortran routine str_count

    Parameters
    ----------
    str : str

    match : str

    Returns
    -------
    num : int
    """

def str_downcase(src: str) -> str:
    """
    Wrapper for Fortran routine str_downcase

    Parameters
    ----------
    src : str

    Returns
    -------
    dst : str
    """

def str_first_in_set(line: str, set: str, ignore_clauses: bool | None = None) -> int:
    """
    Wrapper for Fortran routine str_first_in_set

    Parameters
    ----------
    line : str

    set : str

    ignore_clauses : bool, optional

    Returns
    -------
    ix_match : int
    """

def str_first_not_in_set(line: str, set: str) -> int:
    """
    Wrapper for Fortran routine str_first_not_in_set

    Parameters
    ----------
    line : str

    set : str

    Returns
    -------
    ix_match : int
    """

def str_last_in_set(line: str, set: str) -> int:
    """
    Wrapper for Fortran routine str_last_in_set

    Parameters
    ----------
    line : str

    set : str

    Returns
    -------
    ix_match : int
    """

def str_last_not_in_set(line: str, set: str) -> int:
    """
    Wrapper for Fortran routine str_last_not_in_set

    Parameters
    ----------
    line : str

    set : str

    Returns
    -------
    ix_match : int
    """

def str_match_wild(str: str, pat: str) -> bool:
    """
    Wrapper for Fortran routine str_match_wild

    Parameters
    ----------
    str : str

    pat : str

    Returns
    -------
    a_match : bool
    """

def str_substitute(string: str, str_match: str | None = None, str_replace: str | None = None, do_trim: bool | None = None, ignore_escaped: bool | None = None) -> None:
    """
    Wrapper for Fortran routine str_substitute

    Parameters
    ----------
    string : str

    str_match : str, optional

    str_replace : str, optional

    do_trim : bool, optional

    ignore_escaped : bool, optional
    """

def str_upcase(src: str) -> str:
    """
    Wrapper for Fortran routine str_upcase

    Parameters
    ----------
    src : str

    Returns
    -------
    dst : str
    """

def string_to_int(line: str, default_: int, err_flag: bool, err_print_flag: bool | None = None) -> int:
    """
    Wrapper for Fortran routine string_to_int

    Parameters
    ----------
    line : str

    err_flag : bool

    err_print_flag : bool, optional

    Returns
    -------
    value : int
    """

def string_to_real(line: str, default_: float, err_flag: bool, err_print_flag: bool | None = None) -> float:
    """
    Wrapper for Fortran routine string_to_real

    Parameters
    ----------
    line : str

    err_flag : bool

    err_print_flag : bool, optional

    Returns
    -------
    value : float
    """

def string_trim(in_string: str, out_string: str, word_len: int) -> None:
    """
    Wrapper for Fortran routine string_trim

    Parameters
    ----------
    in_string : str

    out_string : str

    word_len : int
    """

def string_trim2(in_str: str, delimitors: str, out_str: str, ix_word: int, delim: str, ix_next: int) -> None:
    """
    Wrapper for Fortran routine string_trim2

    Parameters
    ----------
    in_str : str

    delimitors : str

    out_str : str

    ix_word : int

    delim : str

    ix_next : int
    """

def suggest_lmdif(XV: _pybmad.RealArray1D, FV: _pybmad.RealArray1D, eps: float, itermx: int, reset_flag: bool | None = None) -> bool:
    """
    Reverse communication subroutine.
    It suggests values for your input variables based on
    the previous value of your merit function.

    Use initial_lmdif to initialize internal variables

    Parameters
    ----------
    xv : 1D array of float
        Array of variables
        This parameter is an input/output and is modified in-place.
        As an output, xv: Suggested new values

    fv : 1D array of float
        Array of function value/s that should be optimized to zero
        This parameter is an input/output and is modified in-place.
        As an output, fv: After the last optimization this returns the best values ever.

    eps : float
        Desired accuracy with which the optimum should be found.

    itermx : int
        Max number of iterations

    reset_flag : bool, optional
        Optional. Used by initial_lmdif to clear previous saved values

    Returns
    -------
    at_end : bool
        Set to False if more optimization is recommended. If set to True then xv(:) will be the minimum found.
    """

def super_bicubic_coef(y: Sequence[float], y1: Sequence[float], y2: Sequence[float], y12: Sequence[float], d1: float, d2: float) -> list[list[float]]:
    """
    Routine to compute coefficients for bicubic interpolation.
    This is from NR bcucof.

    Parameters
    ----------
    y : 1D array of float (shape: 4)
        Function values at grid points.

    y1 : 1D array of float (shape: 4)
        dy/dx1 derivatives.

    y2 : 1D array of float (shape: 4)
        dy/dx2 derivatives.

    y12 : 1D array of float (shape: 4)
        d2y/dx1*dx2 second derivatives.

    d1 : float
        Grid width in 1-direction.

    d2 : float
        Grid width in 2-direction.

    Returns
    -------
    c : 2D array of float (shape: 4,4)
        Coefficients.
    """

class SuperBicubicInterpolation:
    """super_bicubic_interpolation return type"""

    @property
    def ansy(self) -> float: ...

    @property
    def ansy1(self) -> float: ...

    @property
    def ansy2(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def super_bicubic_interpolation(y: Sequence[float], y1: Sequence[float], y2: Sequence[float], y12: Sequence[float], x1l: float, x1u: float, x2l: float, x2u: float, x1: float, x2: float) -> SuperBicubicInterpolation:
    """
    Routine to do bicubic interpolation.
    This is from NR bcuint.

    Note! The four grid points are arrayed in counter-clockwise order beginning from the lower left.
    So, for example, y = [y_ll, y_lu, y_uu, y_ul] where "l" = lower, "u" = upper index.

    Parameters
    ----------
    y : 1D array of float (shape: 4)
        Function values at grid points.

    y1 : 1D array of float (shape: 4)
        dy/dx1 derivatives.

    y2 : 1D array of float (shape: 4)
        dy/dx2 derivatives.

    y12 : 1D array of float (shape: 4)
        d2y/dx1*dx2 second derivatives.

    x1l : float
        1-direction coordinate at lower points.

    x1u : float
        1-direction coordinate at upper points

    x2l : float
        2-direction coordinate at lower points.

    x2u : float
        2-direction coordinate at upper points

    x1 : float
        1-direction coordinate at point to evaluate.

    x2 : float
        2-direction coordinate at point to evaluate.

    Returns
    -------
    ansy : float
        Interpolation value.

    ansy1 : float
        1-direction derivative at interpolation point.

    ansy2 : float
        2-direction derivative at interpolation point.
    """

class SuperPolint:
    """super_polint return type"""

    @property
    def y(self) -> float: ...

    @property
    def dy(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def super_polint(xa: _pybmad.RealArray1D, ya: _pybmad.RealArray1D, x: float) -> SuperPolint:
    """
    This is essentially polint from Numerical Recipes.

    Parameters
    ----------
    xa : 1D array of float

    ya : 1D array of float

    x : float

    Returns
    -------
    y : float

    dy : float
    """

def super_poly(x: float, coeffs: _pybmad.RealArray1D) -> float:
    """
    Routine to compute Sum: coef(i)*x^i

    Parameters
    ----------
    x : float
        Variable.

    Returns
    -------
    value : float
        Polynomial value.
    """

def super_sobseq(x: _pybmad.RealArray1D, ran_state: _pybmad.RandomStateStruct | None = None) -> None:
    """
    Routine patterened after sobseq in Numerical Recipes.
    Difference is that this version has an argument for the internal state.

    Parameters
    ----------
    x : 1D array of float
        Random vector.

    ran_state : RandomStateStruct, optional
        Generator state. See the ran_seed_put documentation for more details.
    """

def super_sort(arr: _pybmad.IntArray1D) -> None:
    """
    Routine to sort an integer array in place.
    This is the NR routine sort modified to sort integers.

    Parameters
    ----------
    arr : 1D array of int
        Array of integers.
        This parameter is an input/output and is modified in-place.
        As an output, arr: Sorted array.
    """

def system_command(line: str) -> bool:
    """
    Wrapper for Fortran routine system_command

    Parameters
    ----------
    line : str
        Command to pass to the system shell.

    Returns
    -------
    err_flag : bool, optional
        Set True if there is an error (bad command, etc.).
    """

def test_xgelbd() -> None:
    """Wrapper for Fortran routine test_xgelbd"""

def to_str(num: float, max_signif: int | None = None) -> str:
    """
    no longer exists
    subroutine test_tune_tracker_lock (tracker_locked)
      implicit none
      logical tracker_locked(2)
    end subroutine
    """

class TricubicCmplxEval:
    """tricubic_cmplx_eval return type"""

    @property
    def df_dx(self) -> complex: ...

    @property
    def df_dy(self) -> complex: ...

    @property
    def df_dz(self) -> complex: ...

    @property
    def f_val(self) -> complex: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tricubic_cmplx_eval(x_norm: float, y_norm: float, z_norm: float, tri_coef: _pybmad.TricubicCmplxCoefStruct) -> TricubicCmplxEval:
    """
    Routine to evaluate a tricubic interpolating complex function.

    Use the routine tricubic_interpolation_cmplx_coefs to generate tri_coef.

    Note: In the equations below, the eight points of the grid box being interpolated range
    from (x0, y0, z0) to (x0+dx, y0+dy, z0+dz).

    Parameters
    ----------
    x_norm : float
        x_norm = (x - x0) / dx

    y_norm : float
        y_norm = (y - y0) / dy

    z_norm : float
        z_norm = (z - z0) / dz

    tri_coef : TricubicCmplxCoefStruct
        Coefficients.

    Returns
    -------
    f_val : complex
        Value of f.

    df_dx : complex, optional
        Normalized first derivative: True df/dx = df_dx * dx

    df_dy : complex, optional
        Normalized first derivative: True df/dy = df_dy * dy

    df_dz : complex, optional
        Normalized first derivative: True df/dz = df_dz * dz
    """

def type_this_file(filename: str) -> None:
    """
    Wrapper for Fortran routine type_this_file

    Parameters
    ----------
    filename : str
    """

def upcase_string(string: str) -> None:
    """
    Wrapper for Fortran routine upcase_string

    Parameters
    ----------
    string : str
    """

def value_of_all_ptr(a_ptr: _pybmad.AllPointerStruct) -> float:
    """
    Wrapper for Fortran routine value_of_all_ptr

    Parameters
    ----------
    a_ptr : AllPointerStruct

    Returns
    -------
    value : float
    """

def virtual_memory_usage() -> int:
    """
    Wrapper for Fortran routine virtual_memory_usage

    Parameters
    ----------

    Returns
    -------
    usage : int
    """

class WMatToAxisAngle:
    """w_mat_to_axis_angle return type"""

    @property
    def axis(self) -> list[float]: ...

    @property
    def angle(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def w_mat_to_axis_angle(w_mat: Sequence[Sequence[float]]) -> WMatToAxisAngle:
    """
    Routine to find the rotation axis and rotation angle corresponding to a given
    3D rotation matrix.

    The rotation angle is chosen in the range [0, pi].

    Parameters
    ----------
    w_mat : 2D array of float (shape: 3,3)
        Rotation matrix.

    Returns
    -------
    axis : 1D array of float (shape: 3)
        Rotation axis. Normalized to 1.

    angle : float
        Rotation angle in the range [0, pi].
    """

def w_mat_to_quat(w_mat: Sequence[Sequence[float]]) -> list[float]:
    """
    Routine to find the quaternion corresponding to a given 3D rotation matrix.

    Parameters
    ----------
    w_mat : 2D array of float (shape: 3,3)
        Rotation matrix

    Returns
    -------
    quat : 1D array of float (shape: 0:3)
        Quaternion.
    """

def word_len(wording: str) -> int:
    """
    Wrapper for Fortran routine word_len

    Parameters
    ----------
    wording : str

    Returns
    -------
    wlen : int
    """

def word_read(in_str: str, delim_list: str, word: str, ix_word: int, delim: str, delim_found: bool, out_str: str, ignore_interior: bool | None = None) -> None:
    """
    Wrapper for Fortran routine word_read

    Parameters
    ----------
    in_str : str

    delim_list : str

    word : str

    ix_word : int

    delim : str

    delim_found : bool

    out_str : str

    ignore_interior : bool, optional
    """

def x0_radiation_length(species: int) -> float:
    """
    Routine to return the X0 raidation length for atomes.

    Parameters
    ----------
    Species : int
        Species ID.

    Returns
    -------
    x0 : float
        Radiation length in kg/m^2. Set to real_garbage$ if species is not atomic or has atomic index greater than
        92.
    """

def zig_table_init() -> None:
    """
    Private routine to initialize the Ziggurat lookup tables on first use.
    Based on Marsaglia & Tsang (2000), "The Ziggurat Method for Generating
    Random Variables", Journal of Statistical Software 5(8).
    """
