"""Tao routines"""

from collections.abc import Sequence
from typing import overload

import _pybmad


def integrate_max(ix_start: int, ix_ele: int, datum_value: float, ix_m: int, branch: _pybmad.BranchStruct, vec: _pybmad.RealArray1D, datum: _pybmad.TaoDataStruct) -> None:
    """
    Wrapper for Fortran routine integrate_max

    Parameters
    ----------
    ix_start : int

    ix_ele : int

    datum_value : float

    ix_m : int

    branch : BranchStruct

    vec : 1D array of float

    datum : TaoDataStruct
    """

def integrate_min(ix_start: int, ix_ele: int, datum_value: float, ix_m: int, branch: _pybmad.BranchStruct, vec: _pybmad.RealArray1D, datum: _pybmad.TaoDataStruct) -> None:
    """
    Wrapper for Fortran routine integrate_min

    Parameters
    ----------
    ix_start : int

    ix_ele : int

    datum_value : float

    ix_m : int

    branch : BranchStruct

    vec : 1D array of float

    datum : TaoDataStruct
    """

def tao_abort_command_file(force_abort: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_abort_command_file

    Parameters
    ----------
    force_abort : bool, optional
        : If present and True, ignore s.global.cmd_file_abort_on_error and abort any open command files.
    """

def tao_add_to_normal_mode_h_array(h_str: str) -> _pybmad.ResonanceHStructAlloc1D:
    """
    Routine to add on to the "h(:)" array holding the list of normal form
    resonance driving terms to calculate.
    If h_str is already in the h_array(:) list, nothing is done.

    Parameters
    ----------
    h_str : str
        Resonance driving term ID. EG: "110000"

    Returns
    -------
    h_array : 1D array of ResonanceHStruct
        Array of resonance driving terms.
    """

def tao_alias_cmd(alias: str, string: str) -> None:
    """
    Wrapper for Fortran routine tao_alias_cmd

    Parameters
    ----------
    alias : str
        Name of the tao command file.

    string : str
        Command file arguments.
    """

def tao_allocate_data_array(u: _pybmad.TaoUniverseStruct, n_data: int, exact: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_allocate_data_array

    Parameters
    ----------
    u : TaoUniverseStruct

    n_data : int

    exact : bool, optional
    """

def tao_allocate_v1_var(n_v1: int, save_old: bool) -> None:
    """
    Wrapper for Fortran routine tao_allocate_v1_var

    Parameters
    ----------
    n_v1 : int

    save_old : bool
    """

def tao_allocate_var_array(n_var: int, default_good_user: bool) -> None:
    """
    Routine to increase the s%var(:) array size.

    Parameters
    ----------
    n_var : int
        Size of s.var(:) wanted.
    """

def tao_beam_emit_calc(plane: int, emit_type: int, ele: _pybmad.EleStruct, bunch_params: _pybmad.BunchParamsStruct) -> float:
    """
    Wrapper for Fortran routine tao_beam_emit_calc

    Parameters
    ----------
    plane : int
        x_plane$ or y_plane$.

    emit_type : int
        Either projected_emit$ or apparent_emit$

    ele : EleStruct
        Element.

    bunch_params : BunchParamsStruct
        Bunch sigma matrix

    Returns
    -------
    emit : float
        emittance.
    """

def tao_beam_track(u: _pybmad.TaoUniverseStruct, tao_lat: _pybmad.TaoLatticeStruct, ix_branch: int, beam: _pybmad.BeamStruct) -> bool:
    """
    Routine to track a a beam of particles.

    Parameters
    ----------
    u : TaoUniverseStruct
        Universe to track through.

    tao_lat : TaoLatticeStruct
        Structure containing the lattice.

    ix_branch : int
        Branch index to track through.

    beam : BeamStruct
        Initial beam distribution
        This parameter is an input/output and is modified in-place.
        As an output, beam: Final beam distribution.

    Returns
    -------
    calc_ok : bool
        Set True if there were no problems, False otherwise.
    """

def tao_beam_track_endpoint(ele_id: str, lat: _pybmad.LatStruct, branch_str: str, where: str, u: _pybmad.TaoUniverseStruct) -> _pybmad.EleStruct | None:
    r"""
    Wrapper for Fortran routine tao_beam_track_endpoint

    Parameters
    ----------
    ele_id : str
        Name or index of the element.

    lat : LatStruct
        Lattice.

    branch_str : str
        Branch where the tracking is done. \'\' => Branch not specified.

    where : str
        'TRACK_END', 'TRACK_START', etc.. Used for error messages.

    u : TaoUniverseStruct
        Universe beam is being tracked in.

    Returns
    -------
    ele : EleStruct, optional
        Pointer to the track endpoint element. Nullified if error. Note: blank ele_id is handled if "where"
        contains 'END' or 'START'
    """

def tao_branch_index(ix_branch: int) -> int:
    """
    Wrapper for Fortran routine tao_branch_index

    Parameters
    ----------
    ix_branch : int
        Nominal branch number.

    Returns
    -------
    ix_this : int
        Branch number.
    """

def tao_calc_data_at_s_pts(tao_lat: _pybmad.TaoLatticeStruct, curve: _pybmad.TaoCurveStruct, comp_sign: float, good: _pybmad.BoolAlloc1D) -> None:
    """
    Wrapper for Fortran routine tao_calc_data_at_s_pts

    Parameters
    ----------
    tao_lat : TaoLatticeStruct

    curve : TaoCurveStruct

    comp_sign : float

    good : 1D array of bool
    """

def tao_call_cmd(file_name: str, cmd_arg: _pybmad.CharacterAlloc1D | None = None) -> None:
    """
    Wrapper for Fortran routine tao_call_cmd

    Parameters
    ----------
    file_name : str
        Name of the tao command file.

    cmd_arg : 1D array of str, optional
        Command file arguments.
    """

def tao_cbar_wave_anal(plot: _pybmad.TaoPlotStruct) -> None:
    """No docstring available."""

def tao_change_ele(ele_name: str, attrib_name: str, num_str: str, update: bool) -> bool:
    """
    Routine to change a variable in the model lattice.

    Parameters
    ----------
    ele_name : str
        Name of variable or element.

    attrib_name : str
        Attribute name of element.

    num_str : str
        Change in value. A '@' signifies a absolute set. A 'd' signifies a set relative design.

    Returns
    -------
    err_flag : bool
        logical, Set true if there is an error, false otherwise.
    """

def tao_change_tune(branch_str: str, mask_str: str, print_list: bool, dqa_str: str, dqb_str: str) -> bool:
    """
    No docstring available.

    Parameters
    ----------
    branch_str : str
        List of branches to apply tune set to.

    mask_str : str
        List of quadrupoles to veto.

    print_list : bool
        If True, print a list of elements varied and coefficients.

    dqa_str : str
        Expression for dQa tune.

    dqb_str : str
        Expression for dQb tune.

    Returns
    -------
    err_flag : bool
        logical, Set true if there is an error, false otherwise.
    """

def tao_change_var(name: str, num_str: str, silent: bool) -> bool:
    """
    Routine to change a variable in the model lattice.

    Parameters
    ----------
    name : str
        Name of variable or element.

    num_str : str
        Change in value. A '@' signifies a absolute set. A 'd' signifies a set relative design.

    silent : bool
        If True then do not print any info.

    Returns
    -------
    err_flag : bool
        logical, Set true if there is an error, false otherwise.
    """

def tao_change_z_tune(branch_str: str, dq_str: str) -> bool:
    """
    No docstring available.

    Parameters
    ----------
    branch_str : str
        List of branches to apply tune set to.

    dq_str : str
        Expression for dQc tune.

    Returns
    -------
    err_flag : bool
        logical, Set true if there is an error, false otherwise.
    """

def tao_chrom_calc_needed(data_type: str, data_source: str) -> bool:
    """
    Wrapper for Fortran routine tao_chrom_calc_needed

    Parameters
    ----------
    data_type : str

    data_source : str

    Returns
    -------
    do_chrom : bool
    """

def tao_clear_cmd(cmd_line: str) -> None:
    """
    Wrapper for Fortran routine tao_clear_cmd

    Parameters
    ----------
    cmd_line : str
        Should be set to 'maps'.
    """

def tao_clip_cmd(gang: bool, where: str, value1: float, value2: float) -> None:
    """
    Wrapper for Fortran routine tao_clip_cmd

    Parameters
    ----------
    gang : bool
        Gang all data d1 arrays together.

    where : str
        Graph() to clip. Eg: 'top:x'

    value1 : float

    value2 : float
    """

def tao_close_command_file() -> None:
    """Wrapper for Fortran routine tao_close_command_file"""

def tao_cmd_history_record(cmd: str) -> None:
    """Subroutine to record a cmd in the command history stack"""

def tao_cmd_split(cmd_line: str, n_word: int, cmd_word: _pybmad.CharacterAlloc1D, extra_words_is_error: bool, separator: str | None = None) -> bool:
    """
    This routine splits the command line into "words" (everything between separators).

    Parameters
    ----------
    cmd_line : str
        The command line.

    n_word : int
        Maximum number of words to split command line into.

    cmd_word : 1D array of str
        The individual words.

    extra_words_is_error : bool
        are extra words allowed at the end? If True then err argument is set True. If False then cmd_word(n_word)
        will contain everything after the n_word-1 word.

    separator : str, optional
        a list of characters that, besides a blank space, signify a word boundary.

    Returns
    -------
    err : bool
        error in splitting words For example: separator = '-+' cmd_line = 'model-design' cmd_word(1) = 'model'
        cmd_word(2) = '-' cmd_word(3) = 'design'

    Notes
    -----
    Anything between single or double quotes is treated as a single word.
    Quoted words have quote marks removed.
    Whitespace or a separator inside of "{}", "()", or "[]" is ignored.
    Whitespace after or before a comma is ignored.
    """

def tao_command(command_line: str, err: bool) -> bool:
    """
    Wrapper for Fortran routine tao_command

    Parameters
    ----------
    command_line : str
        command line

    err : bool

    Returns
    -------
    err_is_fatal : bool
        Set True on non-recoverable error. False otherwise
    """

def tao_constraint_type_name(datum: _pybmad.TaoDataStruct) -> str:
    """
    Wrapper for Fortran routine tao_constraint_type_name

    Parameters
    ----------
    datum : TaoDataStruct
        Datum

    Returns
    -------
    datum_name : str
        Appropriate name.
    """

def tao_control_tree_list(ele: _pybmad.EleStruct) -> _pybmad.ElePointerStructAlloc1D:
    """
    Wrapper for Fortran routine tao_control_tree_list

    Parameters
    ----------
    ele : EleStruct
        Lattice element to start at.

    Returns
    -------
    tree : 1D array of ElePointerStruct
        Array of elements.
    """

def tao_count_strings(string: str, pattern: str) -> int:
    """
    Wrapper for Fortran routine tao_count_strings

    Parameters
    ----------
    string : str
        the string to look at

    pattern : str
        the search pattern

    Returns
    -------
    num : int
        number of occurances
    """

def tao_create_plot_window() -> None:
    """
    Subroutine to create the plot window.
    This soubroutine knows not to create a second window if one already exists.
    """

def tao_curve_beam_ellipse_setup(curve: _pybmad.TaoCurveStruct) -> None:
    """
    Wrapper for Fortran routine tao_curve_beam_ellipse_setup

    Parameters
    ----------
    curve : TaoCurveStruct
    """

def tao_curve_check_universe(curve: _pybmad.TaoCurveStruct, uni: _pybmad.TaoUniverseStruct) -> bool:
    """
    Routine to check if the universe associated with a curve exists and is on.

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve to check.
        This parameter is an input/output and is modified in-place.
        As an output, curve: Curve.valid set to False if needed.

    uni : TaoUniverseStruct
        Associated universe

    Returns
    -------
    is_ok : bool
        Set True if associated universe exists and is on.
    """

def tao_curve_data_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct, curve: _pybmad.TaoCurveStruct) -> None:
    """
    Wrapper for Fortran routine tao_curve_data_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct

    curve : TaoCurveStruct
    """

def tao_curve_datum_calc(eles: _pybmad.ElePointerStructAlloc1D, plot: _pybmad.TaoPlotStruct, curve: _pybmad.TaoCurveStruct, who: str) -> None:
    """
    Routine to calculate datum values.
    The values are calculated at the end of each eles(:)%ele element.

    Parameters
    ----------
    eles : 1D array of ElePointerStruct
        Array of elements.

    plot : TaoPlotStruct

    curve : TaoCurveStruct
        This parameter is an input/output and is modified in-place.
        As an output, curve: Structure holding the datum values

    who : str
        Where to put the data. Either: "SYMBOL" or "LINE".
    """

def tao_curve_ele_ref(curve: _pybmad.TaoCurveStruct, point_to_ele_ref: bool) -> _pybmad.EleStruct | None:
    """
    Wrapper for Fortran routine tao_curve_ele_ref

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve with ref ele.

    point_to_ele_ref : bool

    Returns
    -------
    ele_track : EleStruct
    """

def tao_curve_ix_uni(curve: _pybmad.TaoCurveStruct) -> int:
    """
    Wrapper for Fortran routine tao_curve_ix_uni

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve.

    Returns
    -------
    ix_uni : int
        Universe index.
    """

def tao_curve_name(curve: _pybmad.TaoCurveStruct, use_region: bool | None = None) -> str:
    """
    Wrapper for Fortran routine tao_curve_name

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve

    use_region : bool, optional
        If present and True then use the region name instead of the plot name. Region name is 'NULL_REGION' if
        there is no assocaited region.

    Returns
    -------
    curve_name : str
        Appropriate name.
    """

class TaoCurveRmsCalc:
    """tao_curve_rms_calc return type"""

    @property
    def rms(self) -> float: ...

    @property
    def mean(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_curve_rms_calc(curve: _pybmad.TaoCurveStruct, who: str) -> TaoCurveRmsCalc:
    """
    Wrapper for Fortran routine tao_curve_rms_calc

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve to analyze.

    who : str
        "LINE" or "SYMBOL".

    Returns
    -------
    rms : float
        RMS. -1 => Curve has no data.

    mean : float
        Mean.
    """

def tao_d2_d1_name(d1: _pybmad.TaoD1DataStruct, show_universe: bool | None = None) -> str:
    """
    Wrapper for Fortran routine tao_d2_d1_name

    Parameters
    ----------
    d1 : TaoD1DataStruct
        Data array.

    show_universe : bool, optional
        Show the datum's universe. Default is True.

    Returns
    -------
    d2_d1_name : str
        Appropriate name.
    """

def tao_d2_data_stuffit(u: _pybmad.TaoUniverseStruct, d2_name: str, n_d1_data: int) -> None:
    """
    Wrapper for Fortran routine tao_d2_data_stuffit

    Parameters
    ----------
    u : TaoUniverseStruct

    d2_name : str

    n_d1_data : int
    """

def tao_data_check(err: bool) -> None:
    """
    Wrapper for Fortran routine tao_data_check

    Parameters
    ----------
    err : bool
    """

def tao_data_coupling_init(branch: _pybmad.BranchStruct) -> None:
    """
    Wrapper for Fortran routine tao_data_coupling_init

    Parameters
    ----------
    branch : BranchStruct
        New lattice branch.
    """

def tao_data_sanity_check(datum: _pybmad.TaoDataStruct, print_err: bool, default_data_type: str, uni: _pybmad.TaoUniverseStruct | None = None) -> bool:
    """
    Wrapper for Fortran routine tao_data_sanity_check

    Parameters
    ----------
    datum : TaoDataStruct
        Datum to check.

    print_err : bool
        Print error message if data is not valid?

    default_data_type : str
        Default data type associated with the datum's d2 structure.

    uni : TaoUniverseStruct, optional
        Universe to use instead of datum.d1.d2.ix_universe

    Returns
    -------
    is_valid : bool
        True if internally consistent.
    """

def tao_data_show_use(d2_data: _pybmad.TaoD2DataStruct, lines: _pybmad.CharacterAlloc1D | None = None, nl: int | None = None) -> None:
    """
    Wrapper for Fortran routine tao_data_show_use

    Parameters
    ----------
    d2_data : TaoD2DataStruct

    lines : 1D array of str, optional

    nl : int, optional
    """

def tao_data_type_substitute(template_: str, curve: _pybmad.TaoCurveStruct, graph: _pybmad.TaoGraphStruct) -> str:
    """
    Routine substitute the appropriate data type string for instances of "#ref" and
    "#comp" in template.

    Additionally, if template does not have a "|" character,
    the string "|" + component will be added at the end of str_out.

    Parameters
    ----------
    curve : TaoCurveStruct
        curve.ele_ref_name is substituted for all instances of "#ref".

    graph : TaoGraphStruct

    Returns
    -------
    str_out : str
        String with substitutions.
    """

def tao_data_useit_plot_calc(curve: _pybmad.TaoCurveStruct, graph: _pybmad.TaoGraphStruct, data: _pybmad.TaoDataStructArray1D, check_s_position: bool) -> str:
    """
    Routine to set the data for plotting.

    Parameters
    ----------
    curve : TaoCurveStruct
        tao_curve_struct

    graph : TaoGraphStruct
        tao_graph_struct

    data : 1D array of TaoDataStruct

    check_s_position : bool
        If present and True then veto data that does not have an s-position.

    Returns
    -------
    most_invalid : str
        String documenting biggest invalid data problem.
    """

def tao_datum_has_associated_ele(data_type: str, branch_geometry: int | None = None) -> int:
    """
    Wrapper for Fortran routine tao_datum_has_associated_ele

    Parameters
    ----------
    data_type : str
        Type of data.

    branch_geometry : int, optional
        Geometry of the associated lattice branch. open$ or closed$.

    Returns
    -------
    has_associated_ele : int
    """

class TaoDatumIntegrate:
    """tao_datum_integrate return type"""

    @property
    def valid_value(self) -> bool: ...

    @property
    def why_invalid(self) -> str: ...

    @property
    def result(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_datum_integrate(datum: _pybmad.TaoDataStruct, branch: _pybmad.BranchStruct, s_pos: _pybmad.RealArray1D, values: _pybmad.RealArray1D) -> TaoDatumIntegrate:
    """
    Routine to calculate the integral, rms, or average of an array of values associated with a datum.

    Parameters
    ----------
    datum : TaoDataStruct
        Datum under consideration.

    branch : BranchStruct
        Associated lattice branch.

    s_pos : 1D array of float
        Array of s-positions of the values.

    values : 1D array of float
        Array of values.

    Returns
    -------
    valid_value : bool
        Set false if, for example, all s_pos(:) are the same.

    why_invalid : str
        Information string if there is a problem.

    result : float
        Integral, rms, or average depending upon datum.merit_type.
    """

def tao_datum_name(datum: _pybmad.TaoDataStruct, show_universe: bool | None = None) -> str:
    """
    Wrapper for Fortran routine tao_datum_name

    Parameters
    ----------
    datum : TaoDataStruct
        Datum

    show_universe : bool, optional
        Show the datum's universe. Default is True.

    Returns
    -------
    datum_name : str
        Appropriate name.
    """

def tao_datum_s_position(datum: _pybmad.TaoDataStruct, ele: _pybmad.EleStruct) -> float:
    """
    Routine to calculate the longitudinal position associated with a datum.

    Parameters
    ----------
    datum : TaoDataStruct
        Datum under conideration.

    ele : EleStruct
        Associated lattice element.

    Returns
    -------
    s_pos : float
        Associated longitudinal position.
    """

def tao_de_optimizer() -> bool:
    """
    Wrapper for Fortran routine tao_de_optimizer

    Returns
    -------
    abort : bool
        Set True if an user stop signal detected.
    """

def tao_deallocate_plot_cache(plot_cache: _pybmad.TaoPlotCacheStructAlloc1D) -> None:
    """
    Wrapper for Fortran routine tao_deallocate_plot_cache

    Parameters
    ----------
    plot_cache : 1D array of TaoPlotCacheStruct
    """

def tao_deallocate_tree(tree: _pybmad.TaoEvalNodeStruct) -> None:
    """
    Routine to deallocate tree%node(:) and everything below it

    Parameters
    ----------
    tree : TaoEvalNodeStruct
        Root of tree to deallocate.
        This parameter is an input/output and is modified in-place.
        As an output, tree: Deallocated tree.
    """

def tao_destroy_plot_window() -> None:
    """Wrapper for Fortran routine tao_destroy_plot_window"""

def tao_dmerit_calc() -> None:
    """No docstring available."""

def tao_dmodel_dvar_calc(force_calc: bool) -> bool:
    """
    Subroutine to calculate the dModel_dVar derivative matrix.

    Parameters
    ----------
    force_calc : bool
        If true then force recalculation of the matrix. If False then only calculate matrix if it doesn't exist.

    Returns
    -------
    err_flag : bool
        Set true if there is an error. False otherwise.
    """

def tao_do_wire_scan(ele: _pybmad.EleStruct, theta: float, beam: _pybmad.BeamStruct) -> float:
    """
    Returns the beam's second moment using the wire along the specified angle.
    Keep in mind that the actual correlation axis is 90 degrees off of the
    wire angle

    This simulates a fast wire scanner that performs the scan over only one
    bunch. Obviously, this isn't realistic. Any dynamic effects will not be
    accounted for!

    Parameters
    ----------
    ele : EleStruct

    theta : float
        wire angle wrt x axis (in degrees)

    beam : BeamStruct
        contains the beam distribution

    Returns
    -------
    moment : float
        second moment along axis specified by angle.
    """

def tao_draw_beam_chamber_wall(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    NOTE: THIS ROUTINE IS NOT CURRENTLY ACITVE (NOT CALLED BY ANY OTHER ROUTINE).

    Subroutine tao_draw_beam_chamber_wall (plot, graph)

    Routine to draw the beam chamber wall.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

class TaoDrawCurveData:
    """tao_draw_curve_data return type"""

    @property
    def have_data(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_draw_curve_data(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct, curve: _pybmad.TaoCurveStruct, have_data: bool) -> TaoDrawCurveData:
    """
    Routine to draw a graph with data and/or variable curves.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph containing the curve.

    curve : TaoCurveStruct
        Curve to draw.

    have_data : bool
        Intitial state.
        This parameter is an input/output and is modified in-place.
        As an output, have_data: Is there any data to plot? Set True if so.

    Returns
    -------
    have_data : bool
        Intitial state.
        This parameter is an input/output and is modified in-place.
        As an output, have_data: Is there any data to plot? Set True if so.
    """

def tao_draw_ele_for_floor_plan(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct, tao_lat: _pybmad.TaoLatticeStruct, ele: _pybmad.EleStruct, ele_shape: _pybmad.TaoEleShapeStruct, label_name: str, offset1: float, offset2: float) -> None:
    """
    Routine to draw one lattice element or one datum location for the floor plan graph.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.

    tao_lat : TaoLatticeStruct
        Lattice containing the element.

    ele : EleStruct
        Element to draw.

    ele_shape : TaoEleShapeStruct
        Shape to draw from s.plot_page.floor_plan.ele_shape(:) array. Will be NULL if no associated shape for this
        element.

    label_name : str
        Shape label.

    offset1 : float
        Transverse distances used to scale the drawing of the element shape.

    offset2 : float
        Transverse distances used to scale the drawing of the element shape.
    """

def tao_draw_floor_plan(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw a floor plan graph.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

def tao_draw_graph_axes(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw a just the graph part of a data graph.
    The calling routine takes care of drawing any curves.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

class TaoDrawHistogramData:
    """tao_draw_histogram_data return type"""

    @property
    def have_data(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_draw_histogram_data(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct, curve: _pybmad.TaoCurveStruct, have_data: bool) -> TaoDrawHistogramData:
    """
    Routine to draw a graph with data and/or variable histograms.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph containing the histogram.

    curve : TaoCurveStruct
        Histogram to draw.

    have_data : bool
        Intitial state.
        This parameter is an input/output and is modified in-place.
        As an output, have_data: Is there any data to plot? Set True if so.

    Returns
    -------
    have_data : bool
        Intitial state.
        This parameter is an input/output and is modified in-place.
        As an output, have_data: Is there any data to plot? Set True if so.
    """

def tao_draw_lat_layout(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw a lattice layout graph.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

def tao_draw_plots(do_clear: bool | None = None) -> None:
    """
    Subroutine to draw the plots on the plot window.

    Parameters
    ----------
    do_clear : bool, optional
        If present and False then call qp_clear_page. This argument is used when drawing PS or GIF.
    """

class TaoEleGeometryWithMisalignments:
    """tao_ele_geometry_with_misalignments return type"""

    @property
    def valid_value(self) -> bool: ...

    @property
    def why_invalid(self) -> str: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_ele_geometry_with_misalignments(datum: _pybmad.TaoDataStruct, ele: _pybmad.EleStruct) -> TaoEleGeometryWithMisalignments:
    """
    Routine to evaluate a floor position with misalignments at a given element.
    This routine is private and not for general use.

    Parameters
    ----------
    datum : TaoDataStruct
        Datum info

    ele : EleStruct
        Lattice element to evaluate at.

    Returns
    -------
    valid_value : bool
        Was able to evalute the datum?

    why_invalid : str
        If not valid, why not.

    value : float
        Datum value.
    """

class TaoEleShapeInfo:
    """tao_ele_shape_info return type"""

    @property
    def e_shape(self) -> _pybmad.TaoEleShapeStruct | None: ...

    @property
    def label_name(self) -> str: ...

    @property
    def y1(self) -> float: ...

    @property
    def y2(self) -> float: ...

    @property
    def ix_shape_min(self) -> int | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_ele_shape_info(ix_uni: int, ele: _pybmad.EleStruct, ele_shapes: _pybmad.TaoEleShapeStructArray1D, ix_shape_min: int | None = None) -> TaoEleShapeInfo:
    """
    Wrapper for Fortran routine tao_ele_shape_info

    Parameters
    ----------
    ix_uni : int
        Universe index.

    ele : EleStruct
        Lattice element.

    ele_shapes : 1D array of TaoEleShapeStruct
        Array of shapes to search.

    ix_shape_min : int, optional
        Index of minimum ele_shape(:) index to start search from. Default is 1.
        This parameter is an input/output and is modified in-place.
        As an output, ix_shape_min: Ele_shape(

    Returns
    -------
    label_name : str
        Label name.

    y1 : float
        shape transverse sizes.

    y2 : float
        shape transverse sizes.

    e_shape : TaoEleShapeStruct, optional
        element shape. Will be nullified if no associated shape.

    ix_shape_min : int, optional
        Index of minimum ele_shape(:) index to start search from. Default is 1.
        This parameter is an input/output and is modified in-place.
        As an output, ix_shape_min: Ele_shape(
    """

def tao_ele_shape_input_to_struct(shape_input: _pybmad.TaoEleShapeInput, namelist_name: str | None = None) -> _pybmad.TaoEleShapeStruct:
    """
    Wrapper for Fortran routine tao_ele_shape_input_to_struct

    Parameters
    ----------
    shape_input : TaoEleShapeInput

    namelist_name : str, optional

    Returns
    -------
    shape_struct : TaoEleShapeStruct
    """

def tao_ele_shape_struct_to_input(shape_struct: _pybmad.TaoEleShapeStruct) -> _pybmad.TaoEleShapeInput:
    """
    Wrapper for Fortran routine tao_ele_shape_struct_to_input

    Parameters
    ----------
    shape_struct : TaoEleShapeStruct

    Returns
    -------
    shape_input : TaoEleShapeInput
    """

class TaoEvalFloorOrbit:
    """tao_eval_floor_orbit return type"""

    @property
    def valid_value(self) -> bool: ...

    @property
    def why_invalid(self) -> str: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_eval_floor_orbit(datum: _pybmad.TaoDataStruct, ele: _pybmad.EleStruct, orbit: _pybmad.CoordStruct, bunch_params: _pybmad.BunchParamsStruct) -> TaoEvalFloorOrbit:
    """
    Routine to evaluate a floor_orbit datum at a given element.
    This routine is private and not for general use.

    Parameters
    ----------
    datum : TaoDataStruct
        Datum info

    ele : EleStruct
        Lattice element to evaluate at.

    orbit : CoordStruct
        Particle orbit at element.

    bunch_params : BunchParamsStruct
        Bunch parameters at element.

    Returns
    -------
    valid_value : bool
        Was able to evalute the datum?

    why_invalid : str
        If not valid, why not.

    value : float
        Datum value.
    """

class TaoEvaluateADatum:
    """tao_evaluate_a_datum return type"""

    @property
    def datum_value(self) -> float: ...

    @property
    def valid_value(self) -> bool: ...

    @property
    def why_invalid(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_a_datum(datum: _pybmad.TaoDataStruct, u: _pybmad.TaoUniverseStruct, tao_lat: _pybmad.TaoLatticeStruct, called_from_lat_calc: bool | None = None, print_err: bool | None = None) -> TaoEvaluateADatum:
    """
    Wrapper for Fortran routine tao_evaluate_a_datum

    Parameters
    ----------
    datum : TaoDataStruct
        What type of datum

    u : TaoUniverseStruct
        Which universe to use.

    tao_lat : TaoLatticeStruct
        Lattice to use.

    called_from_lat_calc : bool, optional
        Default is False. If true, prevents infinite loop of this routine calling tao_lattice_calc

    print_err : bool, optional
        Default is True. If False, do not print an error message.

    Returns
    -------
    datum_value : float
        Value of the datum.

    valid_value : bool
        Set false when there is a problem. Set true otherwise.

    why_invalid : str, optional
        Tells why datum value is invalid.
    """

class TaoEvaluateDatumAtS:
    """tao_evaluate_datum_at_s return type"""

    @property
    def err_str(self) -> str: ...

    @property
    def bad_datum(self) -> bool: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_datum_at_s(datum: _pybmad.TaoDataStruct, tao_lat: _pybmad.TaoLatticeStruct, ele: _pybmad.EleStruct, ele_ref: _pybmad.EleStruct, valid_value: bool) -> TaoEvaluateDatumAtS:
    """
    Routine to evaluate a datum at a given s-position in the lattice

    Parameters
    ----------
    datum : TaoDataStruct
        Datum to evaluate.

    tao_lat : TaoLatticeStruct

    ele : EleStruct
        Evaluation element.

    ele_ref : EleStruct
        Reference element.

    valid_value : bool
        True if evaluation was sucessful. False if not.

    Returns
    -------
    err_str : str
        Error string for printing an error message.

    bad_datum : bool
        True -> datum is malformed. False -> Could evaluate or evaluation problem was not due to the datum itself
        (EG: the lattice was unstable).

    value : float
        Datum value.
    """

class TaoEvaluateElementParameters:
    """tao_evaluate_element_parameters return type"""

    @property
    def err(self) -> bool: ...

    @property
    def values(self) -> _pybmad.RealAlloc1D: ...

    @property
    def info(self) -> _pybmad.TaoExpressionInfoStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_element_parameters(param_name: str, print_err: bool, dflt_source: str, dflt_ele: _pybmad.EleStruct | None = None, dflt_component: str | None = None, dflt_uni: int | None = None, eval_point: int | None = None) -> TaoEvaluateElementParameters:
    """
    Wrapper for Fortran routine tao_evaluate_element_parameters

    Parameters
    ----------
    param_name : str
        parameter name.

    print_err : bool
        Print error message?

    dflt_source : str
        Default source

    dflt_ele : EleStruct, optional
        Default element if not specified by param_name.

    dflt_component : str, optional
        Default component

    dflt_uni : int, optional
        Default universe to use.

    eval_point : int, optional

    Returns
    -------
    err : bool
        True if there is an error in syntax. False otherwise

    values : 1D array of float
        Array of datum values.

    info : 1D array of TaoExpressionInfoStruct, optional
    """

class TaoEvaluateExpression:
    """tao_evaluate_expression return type"""

    @property
    def value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def info(self) -> _pybmad.TaoExpressionInfoStructAlloc1D: ...

    @property
    def stack(self) -> _pybmad.TaoEvalNodeStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_expression(expression: str, n_size: int, use_good_user: bool, print_err: bool | None = None, dflt_component: str | None = None, dflt_source: str | None = None, dflt_ele_ref: _pybmad.EleStruct | None = None, dflt_ele_start: _pybmad.EleStruct | None = None, dflt_ele: _pybmad.EleStruct | None = None, dflt_dat_or_var_index: str | None = None, dflt_uni: int | None = None, dflt_eval_point: int | None = None, dflt_s_offset: float | None = None, dflt_orbit: _pybmad.CoordStruct | None = None, datum: _pybmad.TaoDataStruct | None = None) -> TaoEvaluateExpression:
    r"""
    Wrapper for Fortran routine tao_evaluate_expression

    Parameters
    ----------
    expression : str
        Arithmetic expression.

    n_size : int
        Size of the value array. If the expression evaluates to a a scalar, each value in the value array will get
        this value. If n_size = 0, the natural size is determined by the expression itself.

    use_good_user : bool
        Use the good_user logical in evaluating good(:)

    print_err : bool, optional
        If False then supress evaluation error messages. This does not affect syntax error messages. Default is
        True.

    dflt_component : str, optional
        Component to use if not specified in the expression. 'model' (default), 'base', or 'design'.

    dflt_source : str, optional
        Default source ('lat', 'data', etc.). Default is \'\'.

    dflt_ele_ref : EleStruct, optional
        Default reference element.

    dflt_ele_start : EleStruct, optional
        Default start element for ranges.

    dflt_ele : EleStruct, optional
        Default element to evaluate at.

    dflt_dat_or_var_index : str, optional
        Default datum or variable index to use.

    dflt_uni : int, optional
        Default universe to use. If 0 or not present, use viewed universe.

    dflt_eval_point : int, optional
        Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

    dflt_s_offset : float, optional
        Default offset of eval_point. Default = 0.

    dflt_orbit : CoordStruct, optional
        Default orbit to evaluate at.

    datum : TaoDataStruct, optional
        If present, check to see that the expression does not depend upon a datum that will be evaluated after
        this datum. If so, this is an error.

    Returns
    -------
    value : 1D array of float
        Value of arithmetic expression.

    err_flag : bool
        True on an error. EG: Invalid expression. A divide by zero is not an error but good(:) will be set to
        False.

    info : 1D array of TaoExpressionInfoStruct, optional
        Is the value valid?, etc. Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
        orbit.x[23]|good_user is False.

    stack : 1D array of TaoEvalNodeStruct, optional
        Array of nodes of variable names. This is useful to check what datums or variables are used in the
        expression.
    """

class TaoEvaluateExpressionNew:
    """tao_evaluate_expression_new return type"""

    @property
    def value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def info(self) -> _pybmad.TaoExpressionInfoStructAlloc1D: ...

    @property
    def stack(self) -> _pybmad.TaoEvalNodeStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_expression_new(expression: str, n_size: int, use_good_user: bool, print_err: bool | None = None, dflt_component: str | None = None, dflt_source: str | None = None, dflt_ele_ref: _pybmad.EleStruct | None = None, dflt_ele_start: _pybmad.EleStruct | None = None, dflt_ele: _pybmad.EleStruct | None = None, dflt_dat_or_var_index: str | None = None, dflt_uni: int | None = None, dflt_eval_point: int | None = None, dflt_s_offset: float | None = None, dflt_orbit: _pybmad.CoordStruct | None = None, datum: _pybmad.TaoDataStruct | None = None) -> TaoEvaluateExpressionNew:
    r"""
    Wrapper for Fortran routine tao_evaluate_expression_new

    Parameters
    ----------
    expression : str
        Arithmetic expression.

    n_size : int
        Size of the value array. If the expression evaluates to a a scalar, each value in the value array will get
        this value. If n_size = 0, the natural size is determined by the expression itself.

    use_good_user : bool
        Use the good_user logical in evaluating good(:)

    print_err : bool, optional
        If False then supress evaluation error messages. This does not affect syntax error messages. Default is
        True.

    dflt_component : str, optional
        Component to use if not specified in the expression. 'model' (default), 'base', or 'design'.

    dflt_source : str, optional
        Default source ('lat', 'data', etc.). Default is \'\'.

    dflt_ele_ref : EleStruct, optional
        Default reference element.

    dflt_ele_start : EleStruct, optional
        Default start element for ranges.

    dflt_ele : EleStruct, optional
        Default element to evaluate at.

    dflt_dat_or_var_index : str, optional
        Default datum or variable index to use.

    dflt_uni : int, optional
        Default universe to use. If 0 or not present, use viewed universe.

    dflt_eval_point : int, optional
        Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

    dflt_s_offset : float, optional
        Default offset of eval_point. Default = 0.

    dflt_orbit : CoordStruct, optional
        Default orbit to evaluate at.

    datum : TaoDataStruct, optional
        If present, check to see that the expression does not depend upon a datum that will be evaluated after
        this datum. If so, this is an error.

    Returns
    -------
    value : 1D array of float
        Value of arithmetic expression.

    err_flag : bool
        True on an error. EG: Invalid expression. A divide by zero is not an error but good(:) will be set to
        False.

    info : 1D array of TaoExpressionInfoStruct, optional
        Is the value valid?, etc. Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
        orbit.x[23]|good_user is False.

    stack : 1D array of TaoEvalNodeStruct, optional
        Array of nodes of variable names. This is useful to check what datums or variables are used in the
        expression.
    """

class TaoEvaluateExpressionOld:
    """tao_evaluate_expression_old return type"""

    @property
    def value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    @property
    def info(self) -> _pybmad.TaoExpressionInfoStructAlloc1D: ...

    @property
    def stack(self) -> _pybmad.TaoEvalNodeStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_expression_old(expression: str, n_size: int, use_good_user: bool, print_err: bool | None = None, dflt_component: str | None = None, dflt_source: str | None = None, dflt_ele_ref: _pybmad.EleStruct | None = None, dflt_ele_start: _pybmad.EleStruct | None = None, dflt_ele: _pybmad.EleStruct | None = None, dflt_dat_or_var_index: str | None = None, dflt_uni: int | None = None, dflt_eval_point: int | None = None, dflt_s_offset: float | None = None, dflt_orbit: _pybmad.CoordStruct | None = None, datum: _pybmad.TaoDataStruct | None = None) -> TaoEvaluateExpressionOld:
    r"""
    Wrapper for Fortran routine tao_evaluate_expression_old

    Parameters
    ----------
    expression : str
        Arithmetic expression.

    n_size : int
        Size of the value array. If the expression evaluates to a a scalar, each value in the value array will get
        this value. If n_size = 0, the natural size is determined by the expression itself.

    use_good_user : bool
        Use the good_user logical in evaluating good(:)

    print_err : bool, optional
        If False then supress evaluation error messages. This does not affect syntax error messages. Default is
        True.

    dflt_component : str, optional
        Component to use if not specified in the expression. 'model' (default), 'base', or 'design'.

    dflt_source : str, optional
        Default source ('lat', 'data', etc.). Default is \'\'.

    dflt_ele_ref : EleStruct, optional
        Default reference element.

    dflt_ele_start : EleStruct, optional
        Default start element for ranges.

    dflt_ele : EleStruct, optional
        Default element to evaluate at.

    dflt_dat_or_var_index : str, optional
        Default datum or variable index to use.

    dflt_uni : int, optional
        Default universe to use. If 0 or not present, use viewed universe.

    dflt_eval_point : int, optional
        Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

    dflt_s_offset : float, optional
        Default offset of eval_point. Default = 0.

    dflt_orbit : CoordStruct, optional
        Default orbit to evaluate at.

    datum : TaoDataStruct, optional
        If present, check to see that the expression does not depend upon a datum that will be evaluated after
        this datum. If so, this is an error.

    Returns
    -------
    value : 1D array of float
        Value of arithmetic expression.

    err_flag : bool
        True on an error. EG: Invalid expression. A divide by zero is not an error but good(:) will be set to
        False.

    info : 1D array of TaoExpressionInfoStruct, optional
        Is the value valid?, etc. Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
        orbit.x[23]|good_user is False.

    stack : 1D array of TaoEvalNodeStruct, optional
        Array of nodes of variable names. This is useful to check what datums or variables are used in the
        expression.
    """

class TaoEvaluateLatOrBeamData:
    """tao_evaluate_lat_or_beam_data return type"""

    @property
    def err(self) -> bool: ...

    @property
    def values(self) -> _pybmad.RealAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_lat_or_beam_data(data_name: str, print_err: bool, default_source: str, dflt_ele_ref: _pybmad.EleStruct | None = None, dflt_ele_start: _pybmad.EleStruct | None = None, dflt_ele: _pybmad.EleStruct | None = None, dflt_component: str | None = None, dflt_uni: int | None = None, dflt_eval_point: int | None = None, dflt_s_offset: float | None = None) -> TaoEvaluateLatOrBeamData:
    """
    ! private tao_scratch_values_calc, tao_eval_floor_orbit, tao_ele_geometry_with_misalignments

     Subroutine tao_evaluate_lat_or_beam_data (err, data_name, values, print_err, default_source, default_source,
                   dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni, dflt_eval_point, dflt_s_offset)

     Routine to evaluate data with a lat or beam source of the form:
         <universe>@lat::<data_type>[<ix_ele_start>&<ix_ele>]|<component>

    Parameters
    ----------
    data_name : str
        data name.

    print_err : bool
        Print error message?

    dflt_ele_ref : EleStruct, optional
        Default reference element.

    dflt_ele_start : EleStruct, optional
        Default start element.

    dflt_ele : EleStruct, optional
        Default element to evaluate at.

    dflt_component : str, optional
        Default component: 'model' (default), 'base', or 'design'.

    dflt_uni : int, optional
        Default universe to use.

    dflt_eval_point : int, optional
        Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.

    dflt_s_offset : float, optional
        Default offset of eval_point. Default = 0.

    Returns
    -------
    err : bool
        True if there is an error. False otherwise.

    values : 1D array of float
        Array of datum valuse.
    """

class TaoEvaluateStackOld:
    """tao_evaluate_stack_old return type"""

    @property
    def value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_stack_old(stack: _pybmad.TaoEvalNodeStructArray1D, n_size_in: int, use_good_user: bool, print_err: bool, expression: str, info_in: _pybmad.TaoExpressionInfoStructAlloc1D | None = None) -> TaoEvaluateStackOld:
    """
    Routine to evaluate an expression stack.

    Parameters
    ----------
    stack : 1D array of TaoEvalNodeStruct
        Expression stack

    n_size_in : int
        Desired array size. If the expression evaluates to a a scalar, each value in the value array will get this
        value. If n_size = 0, the natural size is determined by the expression itself.

    use_good_user : bool
        Use the good_user logical in evaluating good(:)

    print_err : bool
        If False then supress evaluation error messages. This does not affect syntax error messages. Default is
        True.

    expression : str
        Original expression. Used for error messages.

    Returns
    -------
    value : 1D array of float
        Value of arithmetic expression.

    err_flag : bool
        True on error. False otherwise
    """

class TaoEvaluateTree:
    """tao_evaluate_tree return type"""

    @property
    def value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_evaluate_tree(tao_tree: _pybmad.TaoEvalNodeStruct, n_size: int, use_good_user: bool, print_err: bool, expression: str, info_in: _pybmad.TaoExpressionInfoStructAlloc1D | None = None) -> TaoEvaluateTree:
    """
    Wrapper for Fortran routine tao_evaluate_tree

    Parameters
    ----------
    tao_tree : TaoEvalNodeStruct
        Expression tree

    n_size : int

    use_good_user : bool
        Use the good_user logical in evaluating good(:)

    print_err : bool
        If False then supress evaluation error messages. This does not affect syntax error messages. Default is
        True.

    expression : str
        Original expression. Used for error messages.

    info_in : 1D array of TaoExpressionInfoStruct, optional

    Returns
    -------
    value : 1D array of float
        Value(s) of the arithmetic expression.

    err_flag : bool
        True on error. False otherwise
    """

def tao_evaluate_tune(q_str: str, q0: float, delta_input: bool) -> float:
    """
    Wrapper for Fortran routine tao_evaluate_tune

    Parameters
    ----------
    q_str : str
        String expression.

    q0 : float
        Default to use if q_str evaluates to zero. Also used to set the integer part of the tune.

    delta_input : bool
        If true then qa_str and qb_str are deltas from present tune.

    Returns
    -------
    q_val : float
        Tune value. Set zero if there is an error.
    """

def tao_expression_hash_substitute(expression_in: str, eval_ele: _pybmad.EleStruct | None = None) -> str:
    """
    Routine to, in the expression, substitute the evaluation lattice element name in place
    of hash ("#") characters. Care is taken to only do this where it makes sense.
    For example, "Q1##3" where here "##3" means the third instance of Q1, does not qualify.

    Specifically, a substitution will be done if the character before the hash and the
    character after are one of:
      [,]-*+/:|@<>, or a blank character, or the beginning or end of the expression

    Parameters
    ----------
    expression_in : str
        Expression.

    eval_ele : EleStruct, optional
        Evaluation element name to substitute in. If not present, expression will not be modified.

    Returns
    -------
    expression_out : str
        Expression with substitutions made.
    """

def tao_expression_tree_to_string(tree: _pybmad.TaoEvalNodeStruct, include_root: bool | None = None, n_node: int | None = None, parent: _pybmad.TaoEvalNodeStruct | None = None) -> str:
    """
    Routine to convert an expression tree to a expression string.

    Parameters
    ----------
    tree : TaoEvalNodeStruct
        Tree to print.

    include_root : bool, optional
        Default is True. If True, do not inculde in the output string the root node. Note: If the root node is of
        type root$, this node is always ignored.

    n_node : int, optional
        Internal use only. Used with recursive calls.

    parent : TaoEvalNodeStruct, optional
        Internal use only. Used with recusive calls.

    Returns
    -------
    str_out : str
        Expression string.
    """

class TaoFindData:
    """tao_find_data return type"""

    @property
    def err(self) -> bool: ...

    @property
    def d2_array(self) -> _pybmad.TaoD2DataArrayStructAlloc1D: ...

    @property
    def d1_array(self) -> _pybmad.TaoD1DataArrayStructAlloc1D: ...

    @property
    def d_array(self) -> _pybmad.TaoDataArrayStructAlloc1D: ...

    @property
    def re_array(self) -> _pybmad.TaoRealPointerStructAlloc1D: ...

    @property
    def log_array(self) -> _pybmad.TaoLogicalArrayStructAlloc1D: ...

    @property
    def str_array(self) -> _pybmad.TaoStringArrayStructAlloc1D: ...

    @property
    def int_array(self) -> _pybmad.TaoIntegerArrayStructAlloc1D: ...

    @property
    def component(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_find_data(data_name: str, ix_uni: int | None = None, dflt_index: str | None = None, print_err: bool | None = None) -> TaoFindData:
    """
    Wrapper for Fortran routine tao_find_data

    Parameters
    ----------
    data_name : str
        The data name type. Eg: "3@orbit.x[2:5,10]|meas"

    ix_uni : int, optional
        Index of default universe to use. If ix_uni = 0 then "viewed" universe will be used. Also, if not present
        then the "viewed" universe will be used.

    dflt_index : str, optional
        If present and non-negative, and if no index is specified by the data_name argument, this index is used in
        the evaluation.

    print_err : bool, optional
        Print error message if data is not found? Default is True.

    Returns
    -------
    err : bool
        Err condition

    d2_array : 1D array of TaoD2DataArrayStruct, optional
        Array of pointers to all the matching d2_data structure. Size(d2_array) = 0 if no structures found.

    d1_array : 1D array of TaoD1DataArrayStruct, optional
        Array of pointers to all the matching d1_data structures. Size(d1_array) = 0 if no structures found.

    d_array : 1D array of TaoDataArrayStruct, optional
        Array of pointers to all the matching tao_data_structs.  Size(d_array) = 0 if no structures found.

    re_array : 1D array of TaoRealPointerStruct, optional
        Array of pointers to real component values.  Size(re_array) = 0 if no structures found.

    log_array : 1D array of TaoLogicalArrayStruct, optional
        Array of pointers to logical component values.  Size(log_array) = 0 if no structures found.

    str_array : 1D array of TaoStringArrayStruct, optional
        Array of pointers to character component values.  Size(str_array) = 0 if no structures found.

    int_array : 1D array of TaoIntegerArrayStruct, optional
        Array of pointers to integer component values.  Size(int_array) = 0 if no structures found.

    component : str, optional
        Name of the component. E.G: 'good_user' set to ' ' if no component present.
    """

class TaoFindPlotRegion:
    """tao_find_plot_region return type"""

    @property
    def err(self) -> bool: ...

    @property
    def region(self) -> _pybmad.TaoPlotRegionStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_find_plot_region(where: str, print_flag: bool | None = None) -> TaoFindPlotRegion:
    """
    Wrapper for Fortran routine tao_find_plot_region

    Parameters
    ----------
    where : str
        Region name.

    print_flag : bool, optional
        If present and False then surpress error messages. Default is True.

    Returns
    -------
    err : bool
        Set True on error. False otherwise.

    region : TaoPlotRegionStruct, optional
        Region found.
    """

class TaoFindPlots:
    """tao_find_plots return type"""

    @property
    def err(self) -> bool: ...

    @property
    def plot(self) -> _pybmad.TaoPlotArrayStructAlloc1D: ...

    @property
    def graph(self) -> _pybmad.TaoGraphArrayStructAlloc1D: ...

    @property
    def curve(self) -> _pybmad.TaoCurveArrayStructAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_find_plots(name: str, where: str, print_flag: bool | None = None, blank_means_all: bool | None = None, only_visible: bool | None = None) -> TaoFindPlots:
    """
    Wrapper for Fortran routine tao_find_plots

    Parameters
    ----------
    name : str
        Name of plot or region.

    where : str
        Where to look: 'TEMPLATE', 'REGION', 'BOTH' For where = 'BOTH', if something is found in a plot region,
        then the templates will not be searched

    print_flag : bool, optional
        If present and False then surpress error messages. Default is True.

    blank_means_all : bool, optional
        If present and True then blank graph or curve fields get  interpreted as "*".

    only_visible : bool, optional
        Default is True. If True and s.global.plot_on = True then only include visible plots. If False then plot
        visible setting is ignored.

    Returns
    -------
    err : bool
        Set True on error. False otherwise.

    plot : 1D array of TaoPlotArrayStruct, optional
        Array of plots. If error => size set to 0.

    graph : 1D array of TaoGraphArrayStruct, optional
        Array of graphs. If error => size set to 0.

    curve : 1D array of TaoCurveArrayStruct, optional
        Array of curves. If error => size set to 0.
    """

class TaoFindVar:
    """tao_find_var return type"""

    @property
    def err(self) -> bool: ...

    @property
    def v1_array(self) -> _pybmad.TaoV1VarArrayStructAlloc1D: ...

    @property
    def v_array(self) -> _pybmad.TaoVarArrayStructAlloc1D: ...

    @property
    def re_array(self) -> _pybmad.TaoRealPointerStructAlloc1D: ...

    @property
    def log_array(self) -> _pybmad.TaoLogicalArrayStructAlloc1D: ...

    @property
    def str_array(self) -> _pybmad.TaoStringArrayStructAlloc1D: ...

    @property
    def component(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_find_var(var_name: str, print_err: bool | None = None, dflt_var_index: str | None = None) -> TaoFindVar:
    """
    Wrapper for Fortran routine tao_find_var

    Parameters
    ----------
    var_name : str
        Name of the variable.

    print_err : bool, optional
        Print error message if data is not found? Default is True.

    dflt_var_index : str, optional
        If present and "[...]" var selection substring is not present, then dflt_var_index will be used. [Do not
        include the brackets in this string.]

    Returns
    -------
    err : bool
        err condition

    v1_array : 1D array of TaoV1VarArrayStruct, optional
        Array of pointers to all the v1_var structures.

    v_array : 1D array of TaoVarArrayStruct, optional
        Array of pointers to the variable data point.

    re_array : 1D array of TaoRealPointerStruct, optional
        Array of pointers to the real component values.

    log_array : 1D array of TaoLogicalArrayStruct, optional
        Array of pointers to logical component values.

    str_array : 1D array of TaoStringArrayStruct, optional
        Array of pointers to character component values.

    component : str, optional
        Name of the component. E.G: 'good_user' set to ' ' if no component present.
    """

def tao_fixer(switch_: str, word1: str, word2: str) -> None:
    """
    Wrapper for Fortran routine tao_fixer

    Parameters
    ----------
    word1 : str
        First word of command.

    word2 : str
        Secton word of command.
    """

class TaoFloorToScreen:
    """tao_floor_to_screen return type"""

    @property
    def x_screen(self) -> float: ...

    @property
    def y_screen(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_floor_to_screen(graph: _pybmad.TaoGraphStruct, r_floor: Sequence[float]) -> TaoFloorToScreen:
    """
    Wrapper for Fortran routine tao_floor_to_screen

    Parameters
    ----------
    graph : TaoGraphStruct
        Graph defining the projection plane.

    r_floor : 1D array of float (shape: 3)

    Returns
    -------
    x_screen : float
        x-coordinate of projected point.

    y_screen : float
        y-coordinate of projected point.
    """

def tao_floor_to_screen_coords(graph: _pybmad.TaoGraphStruct, floor: _pybmad.FloorPositionStruct) -> _pybmad.FloorPositionStruct:
    """
    Wrapper for Fortran routine tao_floor_to_screen_coords

    Parameters
    ----------
    graph : TaoGraphStruct
        Graph defining the projection plane.

    floor : FloorPositionStruct
        3D coordinate.

    Returns
    -------
    screen : FloorPositionStruct
        Projected point
    """

def tao_geodesic_lm_optimizer() -> bool:
    """
    Routine to minimize the merit function by varying variables until
    the "data" as calculated from the model matches the measured data.

    This subroutine is a wrapper for the "geodesic"
    Levenburg - Marquardt method.

    Returns
    -------
    abort : bool
        Set True if an user stop signal detected.
    """

class TaoGetData:
    """tao_get_data return type"""

    @property
    def data_value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def data_weight(self) -> _pybmad.RealAlloc1D: ...

    @property
    def data_meas_value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def data_ix_dModel(self) -> _pybmad.IntAlloc1D: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_get_data() -> TaoGetData:
    """
    Subroutine to get the values of the data used in optimization and put them
    in an array. The data is ordered starting with the first universe

    Returns
    -------
    data_value : 1D array of float, optional
        Data model values.

    data_weight : 1D array of float, optional
        Data weights in the merit function.

    data_meas_value : 1D array of float, optional
        Data values when the data was taken.

    data_ix_dModel : 1D array of int, optional
        Data ix_dModel indices
    """

class TaoGetOptVars:
    """tao_get_opt_vars return type"""

    @property
    def var_value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def var_step(self) -> _pybmad.RealAlloc1D: ...

    @property
    def var_delta(self) -> _pybmad.RealAlloc1D: ...

    @property
    def var_weight(self) -> _pybmad.RealAlloc1D: ...

    @property
    def var_ix(self) -> _pybmad.IntAlloc1D: ...

    @property
    def ignore_if_weight_is_zero(self) -> bool: ...

    @property
    def ignore_if_not_limited(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_get_opt_vars() -> TaoGetOptVars:
    """
    Wrapper for Fortran routine tao_get_opt_vars

    Returns
    -------
    var_value : 1D array of float, optional
        Variable model values.

    var_step : 1D array of float, optional
        Variable step sizes.

    var_delta : 1D array of float, optional
        Variable Merit deltas.

    var_weight : 1D array of float, optional
        Variable weights in the merit function.

    var_ix : 1D array of int, optional
        Variable s.var(:) indexes

    ignore_if_weight_is_zero : bool, optional
        If present and True then ignore all variables whose merit weight is zero.

    ignore_if_not_limited : bool, optional
        If present and True then ignore all variables with limit constraint that are not limited.
    """

def tao_get_user_input(prompt_str: str | None = None, wait_flag: bool | None = None, cmd_in: str | None = None) -> str:
    """
    Subroutine to get the next Tao command. In order of precedence, input may come from:
      1) s%com%cmd string (if s%com%use_cmd_here is set to True).
         Used for recalling commands from the history stack.
      3) A command file.
      4) The cmd_in argument (if present). Used, for example, when interfacing with Python.
      5) The terminal.

    Note: A saved command string is present if a prior input string contained multiple commands.
    For example, the following string is read from a command file or terminal or passed via cmd_in:
            "show ele 1; set opti de; run"
    Then cmd_out would be "show ele 1" and "set opti de; run" would be saved for the next call to this routine.

    Note: In single character mode, the input precedence order is ignored and input is taken from the terminal.

    Parameters
    ----------
    prompt_str : str, optional
        Primpt string to print at terminal. If not present then s.global.prompt_string will be used.

    wait_flag : bool, optional
        Used for single mode: Wait state for get_a_char call.

    cmd_in : str, optional
        Command to be used in place getting user input.

    Returns
    -------
    cmd_out : str
        Command from the user.
    """

def tao_graph_controller_setup(graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_controller_setup

    Parameters
    ----------
    graph : TaoGraphStruct
    """

def tao_graph_data_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_data_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct
    """

def tao_graph_data_slice_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_data_slice_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct
    """

def tao_graph_dynamic_aperture_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_dynamic_aperture_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct
    """

def tao_graph_histogram_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_histogram_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct
    """

def tao_graph_name(graph: _pybmad.TaoGraphStruct, use_region: bool | None = None) -> str:
    """
    Wrapper for Fortran routine tao_graph_name

    Parameters
    ----------
    graph : TaoGraphStruct
        Graph

    use_region : bool, optional
        If present and True then use the region name instead of the plot name. Region name is 'NULL_REGION' if
        there is no assocaited region.

    Returns
    -------
    graph_name : str
        Appropriate name.
    """

def tao_graph_phase_space_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_phase_space_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct
    """

class TaoGraphSMinMaxCalc:
    """tao_graph_s_min_max_calc return type"""

    @property
    def s_min(self) -> float: ...

    @property
    def s_max(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_graph_s_min_max_calc(graph: _pybmad.TaoGraphStruct, branch: _pybmad.BranchStruct) -> TaoGraphSMinMaxCalc:
    """
    Routine to calculate min and max for a graph when plot%x_axis_type is set to "s".

    Parameters
    ----------
    graph : TaoGraphStruct
        Graph to calculate for.

    branch : BranchStruct
        Associated lattice branch.

    Returns
    -------
    s_min : float
        Graph min. May be negative with graph.allow_wrap_around = T.

    s_max : float
        Graph max.
    """

def tao_graph_setup(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Wrapper for Fortran routine tao_graph_setup

    Parameters
    ----------
    plot : TaoPlotStruct

    graph : TaoGraphStruct
    """

class TaoHelp:
    """tao_help return type"""

    @property
    def lines(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def n_lines(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_help(what1: str, what2: str) -> TaoHelp:
    """
    Wrapper for Fortran routine tao_help

    Parameters
    ----------
    what1 : str
        command to query. EG: "show".

    what2 : str
        subcommand to query. EG: "element".

    Returns
    -------
    lines : 1D array of str, optional
        If present then the output will be put in this string array instead of printing to the terminal.

    n_lines : int, optional
        Must be present if lines is present. Number of lines used in the lines(:) array.
    """

def tao_init() -> bool:
    """
    Wrapper for Fortran routine tao_init

    Returns
    -------
    err_flag : bool
        Set Treu if there is an error. False otherwise.
    """

def tao_init_beam_in_universe(u: _pybmad.TaoUniverseStruct, beam_init: _pybmad.BeamInitStruct, track_start: str, track_end: str, comb_ds_save: float) -> None:
    """
    Wrapper for Fortran routine tao_init_beam_in_universe

    Parameters
    ----------
    u : TaoUniverseStruct

    beam_init : BeamInitStruct

    track_start : str

    track_end : str

    comb_ds_save : float
    """

def tao_init_beams(init_file: str) -> None:
    """
    Subroutine to initialize beam stuff.

    Parameters
    ----------
    init_file : str
        Tao initialization file. If blank, there is no file so just use the defaults.
    """

def tao_init_data(data_file: str) -> None:
    """
    Subroutine to initialize the tao data structures.

    Parameters
    ----------
    data_file : str
        Tao data initialization file. If blank, there is no file so just use the defaults.
    """

def tao_init_data_end_stuff() -> None:
    """Wrapper for Fortran routine tao_init_data_end_stuff"""

def tao_init_data_in_universe(u: _pybmad.TaoUniverseStruct, n_d2_add: int, keep_existing_data: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_init_data_in_universe

    Parameters
    ----------
    u : TaoUniverseStruct

    n_d2_add : int

    keep_existing_data : bool, optional
    """

def tao_init_dynamic_aperture(init_file: str) -> None:
    """
    Routine to initalize dynamic aperture simulations.

    Parameters
    ----------
    init_file : str
        File setting dynamic_aperture parameters.
    """

class TaoInitFindElements:
    """tao_init_find_elements return type"""

    @property
    def eles(self) -> _pybmad.ElePointerStructAlloc1D: ...

    @property
    def found_one(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_init_find_elements(u: _pybmad.TaoUniverseStruct, search_string: str, attribute: str | None = None) -> TaoInitFindElements:
    """
    Wrapper for Fortran routine tao_init_find_elements

    Parameters
    ----------
    u : TaoUniverseStruct
        Universe to search

    search_string : str
        What to search for

    attribute : str, optional
        Check that attribute of element is free to vary.

    Returns
    -------
    eles : 1D array of ElePointerStruct
        List of matching elements. Size is zero if no elements found.

    found_one : bool, optional
        Set True if a matching element is found. However: Not set if no matching element found.
    """

def tao_init_global(init_file: str) -> None:
    """
    Subroutine to initialize the tao global structures.

    Parameters
    ----------
    init_file : str
        Tao initialization file. If blank, there is no file so just use the defaults.
    """

def tao_init_lattice(lat_file: str, err_flag: bool) -> None:
    """
    Wrapper for Fortran routine tao_init_lattice

    Parameters
    ----------
    lat_file : str

    err_flag : bool
    """

def tao_init_plotting(plot_file: str) -> None:
    """
    Wrapper for Fortran routine tao_init_plotting

    Parameters
    ----------
    plot_file : str
    """

def tao_init_variables(var_file: str) -> None:
    """
    Subroutine to initialize the tao variable structures.

    Parameters
    ----------
    var_file : str
        Tao variable initialization file. If blank, there is no file so just use the defaults.
    """

class TaoInjectBeam:
    """tao_inject_beam return type"""

    @property
    def beam(self) -> _pybmad.BeamStruct: ...

    @property
    def init_ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_inject_beam(u: _pybmad.TaoUniverseStruct, model: _pybmad.TaoLatticeStruct, ix_branch: int) -> TaoInjectBeam:
    """
    This will initialize the beam for a given lattice branch.

    Trying to inject a beam of one species into a branch with a different ref species
    (example: electron bunch into photon branch) is problematical. To avoid problems, Tao
    will set not inject (init_ok = False) if there is a mismatch.

    Parameters
    ----------
    u : TaoUniverseStruct
        Universe containing the lattice.

    model : TaoLatticeStruct
        Universe parameters.

    ix_branch : int
        Lattice branch index to inject into.

    Returns
    -------
    beam : BeamStruct
        Initial beam.

    init_ok : bool
        Set False if there are problems. True otherwise.
    """

def tao_inject_particle(u: _pybmad.TaoUniverseStruct, model: _pybmad.TaoLatticeStruct, ix_branch: int) -> None:
    """
    Wrapper for Fortran routine tao_inject_particle

    Parameters
    ----------
    u : TaoUniverseStruct

    model : TaoLatticeStruct

    ix_branch : int
    """

class TaoIsValidName:
    """tao_is_valid_name return type"""

    @property
    def why_invalid(self) -> str: ...

    @property
    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_is_valid_name(name: str) -> TaoIsValidName:
    """
    Wrapper for Fortran routine tao_is_valid_name

    Parameters
    ----------
    name : str
        Name to be checked.

    Returns
    -------
    why_invalid : str
        Why invalid description.

    is_valid : bool
        True if valid. False otherwise.
    """

def tao_json_cmd(input_str: str) -> None:
    """
    Wrapper for Fortran routine tao_json_cmd

    Parameters
    ----------
    input_str : str
        What to show.
    """

def tao_key_info_to_str(ix_key: int, ix_min_key: int, ix_max_key: int, key_str: str, header_str: str) -> None:
    """
    Wrapper for Fortran routine tao_key_info_to_str

    Parameters
    ----------
    ix_key : int

    ix_min_key : int

    ix_max_key : int

    key_str : str

    header_str : str
    """

def tao_lat_bookkeeper(u: _pybmad.TaoUniverseStruct, tao_lat: _pybmad.TaoLatticeStruct) -> bool:
    """
    Wrapper for Fortran routine tao_lat_bookkeeper

    Parameters
    ----------
    u : TaoUniverseStruct

    tao_lat : TaoLatticeStruct

    Returns
    -------
    err_flag : bool
        Set True if there is a problem. False otherwise.
    """

def tao_lat_emit_calc(plane: int, emit_type: int, ele: _pybmad.EleStruct, modes: _pybmad.NormalModesStruct) -> float:
    """
    Wrapper for Fortran routine tao_lat_emit_calc

    Parameters
    ----------
    plane : int
        x_plane$ or y_plane$.

    emit_type : int
        Either projected_emit$ or apparent_emit$

    ele : EleStruct
        Element holding the Twiss and coupling parameters.

    modes : NormalModesStruct
        Structure holding the emittances

    Returns
    -------
    emit : float
        emittance.
    """

def tao_lat_sigma_calc_needed(data_type: str, data_source: str) -> bool:
    """
    Wrapper for Fortran routine tao_lat_sigma_calc_needed

    Parameters
    ----------
    data_type : str

    data_source : str

    Returns
    -------
    do_lat_sigma : bool
    """

def tao_lat_sigma_track(tao_lat: _pybmad.TaoLatticeStruct, ix_branch: int, print_err: bool | None = None, force_calc: bool | None = None) -> bool:
    """
    Routine to track the 6x6 sigma matrix through the lattice using the lattice linear transfer matrices.

    Parameters
    ----------
    tao_lat : TaoLatticeStruct
        Structure containing the lattice.

    ix_branch : int
        Branch index to track through.

    print_err : bool, optional
        Default is False. Print error messages if, eg, lattice is unstable?

    force_calc : bool, optional
        Default is False. If True, force the calculation to be done.

    Returns
    -------
    calc_ok : bool
        Set True if there were no problems, False otherwise.
    """

def tao_lattice_branches_equal_tao_lattice_branches(tlb1: _pybmad.TaoLatticeBranchStructArray1D, tlb2: _pybmad.TaoLatticeBranchStructArray1D) -> None:
    """
    Wrapper for Fortran routine tao_lattice_branches_equal_tao_lattice_branches

    Parameters
    ----------
    tlb1 : 1D array of TaoLatticeBranchStruct

    tlb2 : 1D array of TaoLatticeBranchStruct
    """

class TaoLatticeCalc:
    """tao_lattice_calc return type"""

    @property
    def calc_ok(self) -> bool: ...

    @property
    def print_err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_lattice_calc() -> TaoLatticeCalc:
    """
    Wrapper for Fortran routine tao_lattice_calc

    Returns
    -------
    calc_ok : bool
        Set False if there was an error in the calculation like a particle was lost or a lat is unstable.

    print_err : bool, optional
        Default True. If False, do not print error messages if, for example, the lattice is unstable.
    """

def tao_lattice_equal_tao_lattice(lat1: _pybmad.TaoLatticeStruct, lat2: _pybmad.TaoLatticeStruct) -> None:
    """
    Wrapper for Fortran routine tao_lattice_equal_tao_lattice

    Parameters
    ----------
    lat1 : TaoLatticeStruct

    lat2 : TaoLatticeStruct
    """

def tao_limit_calc() -> bool:
    """
    Wrapper for Fortran routine tao_limit_calc

    Returns
    -------
    limited : bool
        Set True if a variable is past a limit.
    """

def tao_lm_optimizer() -> bool:
    """
    Routine to minimize the merit function by varying variables until
    the "data" as calculated from the model matches the measured data.

    This subroutine is a wrapper for the mrqmin routine of Numerical Recipes.
    See the Numerical Recipes writeup for more details.
    'lm' stands for Levenburg - Marquardt. Otherwise known as LMDIF.

    Returns
    -------
    abort : bool
        Set True if an user stop signal detected.
    """

def tao_lmdif_optimizer() -> bool:
    """
    Wrapper for Fortran routine tao_lmdif_optimizer

    Returns
    -------
    abort : bool
        Set True if an user stop signal detected or there is a problem with calculating the merit function.
    """

def tao_load_this_datum(vec: _pybmad.RealArray1D, ele_ref: _pybmad.EleStruct, ele_start: _pybmad.EleStruct, ele: _pybmad.EleStruct, datum_value: float, valid_value: bool, datum: _pybmad.TaoDataStruct, branch: _pybmad.BranchStruct, why_invalid: str | None = None, good: _pybmad.BoolAlloc1D | None = None) -> None:
    """
    Wrapper for Fortran routine tao_load_this_datum

    Parameters
    ----------
    vec : 1D array of float

    ele_ref : EleStruct

    ele_start : EleStruct

    ele : EleStruct

    datum_value : float

    valid_value : bool

    datum : TaoDataStruct

    branch : BranchStruct

    why_invalid : str, optional

    good : 1D array of bool, optional
    """

class TaoLocateAllElements:
    """tao_locate_all_elements return type"""

    @property
    def eles(self) -> _pybmad.ElePointerStructAlloc1D: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_locate_all_elements(ele_list: str, ignore_blank: bool | None = None) -> TaoLocateAllElements:
    """
    Wrapper for Fortran routine tao_locate_all_elements

    Parameters
    ----------
    ele_list : str
        String with element names using element list format.

    ignore_blank : bool, optional
        If present and true then do nothing if ele_list is blank. otherwise a blank is treated as an error.

    Returns
    -------
    eles : 1D array of ElePointerStruct
        : Array of elements in the model lat.

    err : bool
        Set true on error.
    """

class TaoLocateElements:
    """tao_locate_elements return type"""

    @property
    def eles(self) -> _pybmad.ElePointerStructAlloc1D: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_locate_elements(ele_list: str, ix_universe: int, lat_type: int | None = None, ignore_blank: bool | None = None, err_stat_level: int | None = None, above_ubound_is_err: bool | None = None, ix_branch: int | None = None, multiple_eles_is_err: bool | None = None) -> TaoLocateElements:
    """
    Wrapper for Fortran routine tao_locate_elements

    Parameters
    ----------
    ele_list : str
        String with element names using element list format.

    ix_universe : int
        Universe to search. -1 => search s.global.default_universe. -2 (all unis) => error. ix_universe is ignored
        if ele_list starts with a universe specifier "N@".

    lat_type : int, optional
        model$ (default), design$, or base$.

    ignore_blank : bool, optional
        If present and true then do nothing if ele_list is blank. otherwise treated as an error.

    err_stat_level : int, optional
        Status level for error messages. If not present, print with level s_error$. Use s_nooutput$ to prevent
        printing.

    above_ubound_is_err : bool, optional

    ix_branch : int, optional
        If present and non-negative then use this as the branch index for elements specified using an integer
        index (EG: "43"). If -1 use the default branch, search all branches.

    multiple_eles_is_err : bool, optional
        If present and True then matching to more than one element is an error.

    Returns
    -------
    eles : 1D array of ElePointerStruct
        : Array of elements in the model lat.

    err : bool
        Set true on error.
    """

def tao_mark_lattice_ele(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine tao_mark_lattice_ele

    Parameters
    ----------
    lat : LatStruct
        Input lattice
        This parameter is an input/output and is modified in-place.
        As an output, lat: Lattice with elements marked.
    """

class TaoMerit:
    """tao_merit return type"""

    @property
    def calc_ok(self) -> bool: ...

    @property
    def this_merit(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_merit() -> TaoMerit:
    """
    Wrapper for Fortran routine tao_merit

    Returns
    -------
    this_merit : float
        Merit value.

    calc_ok : bool, optional
        Set False if there was an error in the calculation like a particle was lost or a lat is unstable.
    """

class TaoNextSwitch:
    """tao_next_switch return type"""

    @property
    def err(self) -> bool: ...

    @property
    def line(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_next_switch(line: str, switch_list: _pybmad.CharacterAlloc1D, return_next_word: bool, neg_num_not_switch: bool | None = None, print_err: bool | None = None) -> TaoNextSwitch:
    r"""
    Subroutine look at the next word on the command line and match this word to a list of "switches"
    given by the switch_list argument.

    If switch_list(1) starts with a "-" or "#" character, switches are assumed to start with this character.
    If switch_list(1) starts with any other character, everything is considered to be a switch.

    Switch abbreviations are permitted.

    If return_next_word = True then, when a non-switch word is encountered, the switch argument
    will be set to that word and that word will be removed from the line argument.

    If return_next_word = False then, when a non-switch word is encountered, the switch argument
    will be set to \'\' and the non-switch word will be left on the line argument.

    If the first non-blank character in line is a single or double quote. The word returned will be the
    substring from the initial quote mark to the next matching quote mark. The quote marks will be removed
    from the returned switch argument.

    Parameters
    ----------
    line : str
        Command line
        This parameter is an input/output and is modified in-place.
        As an output, line: Command line with first word removed if

    switch_list : 1D array of str
        List of valid switches.

    return_next_word : bool
        See above.

    neg_num_not_switch : bool, optional
        If present and True then a word like "-34" will be treated as a non-switch.

    print_err : bool, optional
        Default is True. If False, do not print unknown switch error.

    Returns
    -------
    line : str
        Command line
        This parameter is an input/output and is modified in-place.
        As an output, line: Command line with first word removed if

    err : bool
        Set True if the next word begins with '-' but there is no match to anything in switch_list.
    """

class TaoNextWord:
    """tao_next_word return type"""

    @property
    def word(self) -> str: ...

    @property
    def line(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_next_word(line: str) -> TaoNextWord:
    """
    Routine to return the next word in a line.

    Words are delimited by a space character except if the space is within quotes.
    Additionally, spaces within brackets "(...)", "{...}", and "[...]" are ignored.
    Outer quote marks will be removed in the returned word.

    Parameters
    ----------
    line : str
        String to parse.
        This parameter is an input/output and is modified in-place.
        As an output, line: String with first word removed.

    Returns
    -------
    line : str
        String to parse.
        This parameter is an input/output and is modified in-place.
        As an output, line: String with first word removed.

    word : str
        First word of line.
    """

def tao_one_turn_map_calc_needed(data_type: str, data_source: str) -> bool:
    """
    Wrapper for Fortran routine tao_one_turn_map_calc_needed

    Parameters
    ----------
    data_type : str

    data_source : str

    Returns
    -------
    do_one_turn_map : bool
    """

def tao_open_file(file: str, file_name: str, error_severity: int, binary: bool | None = None) -> int:
    """
    Wrapper for Fortran routine tao_open_file

    Parameters
    ----------
    file : str

    file_name : str
        File name.

    error_severity : int
        Severity level used in the error message. Possibilities are s_fatal$, etc. See out_io doc for more
        details. Use -1 to not print a message if file cannot be opened.

    binary : bool, optional
        If present and True then open a binary file, Defaut is False.

    Returns
    -------
    iunit : int
        Logical unit number. Set to 0 if file not openable.
    """

class TaoOpenScratchFile:
    """tao_open_scratch_file return type"""

    @property
    def err(self) -> bool: ...

    @property
    def iu(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_open_scratch_file() -> TaoOpenScratchFile:
    """
    Wrapper for Fortran routine tao_open_scratch_file

    Returns
    -------
    err : bool
        Set True if there is an error. False otherwise.

    iu : int
        File handle unit number.
    """

def tao_optimization_status(datum: _pybmad.TaoDataStruct) -> str:
    """
    Wrapper for Fortran routine tao_optimization_status

    Parameters
    ----------
    datum : TaoDataStruct
        Datum to evaluate.

    Returns
    -------
    why_str : str
        Optimization status of the datum.
    """

def tao_orbit_beta_wave_anal(plot: _pybmad.TaoPlotStruct) -> None:
    """No docstring available."""

def tao_oreint_building_wall_pt(pt_in: _pybmad.TaoBuildingWallPointStruct) -> _pybmad.TaoBuildingWallPointStruct:
    """
    Wrapper for Fortran routine tao_oreint_building_wall_pt

    Parameters
    ----------
    pt_in : TaoBuildingWallPointStruct
        Building wall point.

    Returns
    -------
    pt_out : TaoBuildingWallPointStruct
        Building wall point with orientation params applied.
    """

class TaoParamValueAtS:
    """tao_param_value_at_s return type"""

    @property
    def err_flag(self) -> bool: ...

    @property
    def why_invalid(self) -> str: ...

    @property
    def print_err(self) -> bool: ...

    @property
    def bad_datum(self) -> bool: ...

    @property
    def value(self) -> float: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_param_value_at_s(dat_name: str, ele_to_s: _pybmad.EleStruct, ele_here: _pybmad.EleStruct, orbit: _pybmad.CoordStruct) -> TaoParamValueAtS:
    """
    Wrapper for Fortran routine tao_param_value_at_s

    Parameters
    ----------
    dat_name : str

    ele_to_s : EleStruct
        Element whose exit end is at the evaluation s-position.

    ele_here : EleStruct
        Lattice element that overlaps the s-position ele.s.

    orbit : CoordStruct
        Orbit at the evaluation s-position.

    Returns
    -------
    err_flag : bool
        Set true if parameter cannot be evaluated.

    value : float
        Parameter value.

    why_invalid : str, optional
        Set if  err_flag = True to document why is there a problem.

    print_err : bool, optional
        Print error message on error? Default is True.

    bad_datum : bool, optional
        Data_type is malformed.
    """

def tao_param_value_routine(str: str, use_good_user: bool, saved_prefix: str, stack: _pybmad.TaoEvalNodeStruct, err_flag: bool, print_err: bool, dflt_component: str | None = None, dflt_source: str | None = None, dflt_ele_ref: _pybmad.EleStruct | None = None, dflt_ele_start: _pybmad.EleStruct | None = None, dflt_ele: _pybmad.EleStruct | None = None, dflt_dat_or_var_index: str | None = None, dflt_uni: int | None = None, dflt_eval_point: int | None = None, dflt_s_offset: float | None = None, dflt_orbit: _pybmad.CoordStruct | None = None, datum: _pybmad.TaoDataStruct | None = None) -> None:
    """
    Wrapper for Fortran routine tao_param_value_routine

    Parameters
    ----------
    str : str

    use_good_user : bool

    saved_prefix : str

    stack : TaoEvalNodeStruct

    err_flag : bool

    print_err : bool

    dflt_component : str, optional

    dflt_source : str, optional

    dflt_ele_ref : EleStruct, optional

    dflt_ele_start : EleStruct, optional

    dflt_ele : EleStruct, optional

    dflt_dat_or_var_index : str, optional

    dflt_uni : int, optional

    dflt_eval_point : int, optional

    dflt_s_offset : float, optional

    dflt_orbit : CoordStruct, optional

    datum : TaoDataStruct, optional
    """

def tao_parse_command_args(cmd_line: str | None = None) -> bool:
    """
    Wrapper for Fortran routine tao_parse_command_args

    Parameters
    ----------
    cmd_line : str, optional

    Returns
    -------
    error : bool
        Set True if there is an error. False otherwise.
    """

class TaoParseElementParamStr:
    """tao_parse_element_param_str return type"""

    @property
    def err(self) -> bool: ...

    @property
    def uni(self) -> str: ...

    @property
    def element(self) -> str: ...

    @property
    def parameter(self) -> str: ...

    @property
    def where(self) -> int: ...

    @property
    def component(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_parse_element_param_str(in_str: str) -> TaoParseElementParamStr:
    """
    Wrapper for Fortran routine tao_parse_element_param_str

    Parameters
    ----------
    in_str : str
        String specifying a parameter of an element or elements.

    Returns
    -------
    err : bool
        Set True if there is a parse error. False otherwise.

    uni : str
        Universe substring.

    element : str
        Element name.

    parameter : str
        Element parameter name.

    where : int
        One of not_set$, anchor_beginning$, anchor_center$, or anchor_end$.

    component : str
        One of "model", "design", or "base".
    """

class TaoParticleDataValue:
    """tao_particle_data_value return type"""

    @property
    def value(self) -> _pybmad.RealAlloc1D: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_particle_data_value(data_type: str, p: _pybmad.CoordStructArray1D, ele: _pybmad.EleStruct, ix_bunch: int) -> TaoParticleDataValue:
    """
    Routine to calculate the value array of a data_type for an array of particles.

    Parameters
    ----------
    data_type : str
        Type of data.

    p : 1D array of CoordStruct
        coord_struct, Array of particles containing the data.

    ele : EleStruct
        Needed for "Ja" evaluation.

    ix_bunch : int
        Bunch index.

    Returns
    -------
    value : 1D array of float
        Array of values.

    err : bool
        Set True if there is an error. False otherwise.
    """

def tao_pause_cmd(time: float) -> None:
    """
    Wrapper for Fortran routine tao_pause_cmd

    Parameters
    ----------
    time : float
        Time to pause in seconds.
    """

def tao_phase_space_axis_index(data_type: str, err: bool) -> int:
    """
    Routine to calculate the phase space axis index for a given data type.

    Parameters
    ----------
    data_type : str
        Type of data.

    err : bool
        Set True if there is an error.

    Returns
    -------
    ix_axis : int
        Axis index.
    """

def tao_phase_wave_anal(plot: _pybmad.TaoPlotStruct) -> None:
    """No docstring available."""

class TaoPickUniverse:
    """tao_pick_universe return type"""

    @property
    def name_out(self) -> str: ...

    @property
    def picked(self) -> _pybmad.BoolAlloc1D: ...

    @property
    def err(self) -> bool: ...

    @property
    def ix_uni(self) -> int: ...

    @property
    def explicit_uni(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_pick_universe(name_in: str, dflt_uni: int | None = None, pure_uni: bool | None = None) -> TaoPickUniverse:
    """
    Wrapper for Fortran routine tao_pick_universe

    Parameters
    ----------
    name_in : str
        data name with possible universe spec.

    dflt_uni : int, optional
        Default universe to use. Set to -1 if explicit universe is required.

    pure_uni : bool, optional
        Default is False. See above

    Returns
    -------
    name_out : str
        name_in without any "n@" beginning.

    picked : 1D array of bool
        Array showing picked universes. The array will be resized if necessary.

    err : bool
        Set True if an error is detected.

    ix_uni : int, optional
        Set to the picked universe with the highest index.

    explicit_uni : bool, optional
        Set True if name_in has explicit universe "n@" specification.
    """

def tao_pipe_cmd(input_str: str) -> None:
    """
    Wrapper for Fortran routine tao_pipe_cmd

    Parameters
    ----------
    input_str : str
        What to show.
    """

def tao_place_cmd(where: str, who: str, no_buffer: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_place_cmd

    Parameters
    ----------
    where : str
        Region where the plot goes. Eg: 'top'.

    who : str
        Type of plot. Eg: 'orbit'.

    no_buffer : bool, optional
        If present and True then prevents buffering in the case when s.global.external_plotting = T
    """

def tao_plot_cmd(where: str, component: str) -> None:
    """
    Wrapper for Fortran routine tao_plot_cmd

    Parameters
    ----------
    where : str
        Region name to identify the plot to set.

    component : str
        Who to plot. EG: 'meas - design'
    """

def tao_plot_data(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw a graph with data and/or variable curves.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

def tao_plot_histogram(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw one graph for the histogram analysis plot.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

def tao_plot_key_table(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw a key table graph.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

def tao_plot_setup() -> None:
    """Wrapper for Fortran routine tao_plot_setup"""

def tao_plot_struct_transfer(plot_in: _pybmad.TaoPlotStruct) -> _pybmad.TaoPlotStruct:
    """
    Wrapper for Fortran routine tao_plot_struct_transfer

    Parameters
    ----------
    plot_in : TaoPlotStruct
        Input structure.

    Returns
    -------
    plot_out : TaoPlotStruct
        Output struture.
    """

def tao_plot_wave(plot: _pybmad.TaoPlotStruct, graph: _pybmad.TaoGraphStruct) -> None:
    """
    Routine to draw one graph for the wave analysis plot.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot containing the graph.

    graph : TaoGraphStruct
        Graph to plot.
    """

class TaoPointerToBranches:
    """tao_pointer_to_branches return type"""

    @property
    def branches(self) -> _pybmad.BranchPointerStructAlloc1D: ...

    @property
    def unis(self) -> _pybmad.TaoUniversePointerStructAlloc1D: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_pointer_to_branches(branch_str: str) -> TaoPointerToBranches:
    """
    Wrapper for Fortran routine tao_pointer_to_branches

    Parameters
    ----------
    branch_str : str
        String specifying what branches to use.

    Returns
    -------
    branches : 1D array of BranchPointerStruct
        Array of pointers to branches.

    unis : 1D array of TaoUniversePointerStruct
        Array of corresponding universes.

    err : bool
        Set True if there is an error.
    """

def tao_pointer_to_building_wall_shape(wall_name: str) -> _pybmad.TaoEleShapeStruct | None:
    """
    Wrapper for Fortran routine tao_pointer_to_building_wall_shape

    Parameters
    ----------
    wall_name : str
        Name of the wall.

    Returns
    -------
    e_shape : TaoEleShapeStruct, optional
        Associated shape. Nullified if there is no associated shape.
    """

def tao_pointer_to_datum(d1: _pybmad.TaoD1DataStruct, ele_name: str) -> _pybmad.TaoDataStruct | None:
    """
    Wrapper for Fortran routine tao_pointer_to_datum

    Parameters
    ----------
    d1 : TaoD1DataStruct
        D1 data struct to search.

    ele_name : str
        Name of lattice element to match to.

    Returns
    -------
    datum_ptr : TaoDataStruct, optional
        Pointer to the matched datum. Will be null if no match found.
    """

class TaoPointerToDatumEle:
    """tao_pointer_to_datum_ele return type"""

    @property
    def valid(self) -> bool: ...

    @property
    def why_invalid(self) -> str: ...

    @property
    def ele(self) -> _pybmad.EleStruct | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_pointer_to_datum_ele(lat: _pybmad.LatStruct, ele_name: str, ix_ele: int, datum: _pybmad.TaoDataStruct, print_err: bool | None = None) -> TaoPointerToDatumEle:
    """
    Routine to see if an element index corresponds to an element with a definite
    location such as an overlay or multipass element.

    If the element is a super_lord then the super_slave element at the exit end
    of the lord will be returned. Otherwise ix_loc will be set to ix_ele.

    Parameters
    ----------
    lat : LatStruct
        Lattice

    ix_ele : int
        Index of element.

    datum : TaoDataStruct
        Used for error messages and gives branch index.

    print_err : bool, optional
        Default is True. If False, do not print an error message.

    Returns
    -------
    valid : bool
        Set False if element does not have a definite location. Set True otherwise

    why_invalid : str, optional
        Tells why datum value is invalid.

    ele : EleStruct, optional
        : Pointer to the element. Set to NULL if not valid or no associated element.
    """

class TaoPointerToEleShape:
    """tao_pointer_to_ele_shape return type"""

    @property
    def dat_var_name(self) -> str: ...

    @property
    def dat_var_value(self) -> float: ...

    @property
    def e_shape(self) -> _pybmad.TaoEleShapeStruct | None: ...

    @property
    def ix_shape_min(self) -> int | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_pointer_to_ele_shape(ix_uni: int, ele: _pybmad.EleStruct, ele_shape: _pybmad.TaoEleShapeStructArray1D, ix_shape_min: int | None = None) -> TaoPointerToEleShape:
    """
    Wrapper for Fortran routine tao_pointer_to_ele_shape

    Parameters
    ----------
    ix_uni : int
        Universe index.

    ele : EleStruct
        Lattice element.

    ele_shape : 1D array of TaoEleShapeStruct
        Array of shapes to search.

    ix_shape_min : int, optional
        Index of minimum ele_shape(:) index to start search from. Default is 1.
        This parameter is an input/output and is modified in-place.
        As an output, ix_shape_min: Ele_shape(

    Returns
    -------
    dat_var_name : str, optional
        Name of datum or variable associated with e_shape. Will be set to "" if there is no associated datum or
        variable.

    dat_var_value : float, optional
        Value of datum or variable associated with e_shape. Will be set to zero if there is no associated datum or
        variable.

    ix_shape_min : int, optional
        Index of minimum ele_shape(:) index to start search from. Default is 1.
        This parameter is an input/output and is modified in-place.
        As an output, ix_shape_min: Ele_shape(

    e_shape : TaoEleShapeStruct, optional
        Associated shape. Nullified if there is no associated shape.
    """

def tao_pointer_to_tao_lat(u: _pybmad.TaoUniverseStruct, lat_type: int | None = None) -> _pybmad.TaoLatticeStruct | None:
    """
    Wrapper for Fortran routine tao_pointer_to_tao_lat

    Parameters
    ----------
    u : TaoUniverseStruct
        Universe to work with

    lat_type : int, optional
        model$ (default), design$, or base$.

    Returns
    -------
    tao_lat : TaoLatticeStruct, optional
        Tao_lat pointer. Points to u.model, u.design, or u.base
    """

@overload
def tao_pointer_to_universe(ix_uni: int, neg2_to_default: bool | None = None) -> _pybmad.TaoUniverseStruct | None:
    """
    Routine to set a pointer to a universe.

    This is an overloaded routine for the:
     tao_pointer_to_universe_int (ix_uni, neg2_to_default) result (u)
     tao_pointer_to_universe_str (string, neg2_to_default) result (u)

    Note: With a string argument, this routine can only handle single universe picks.
    That is, it cannot handlle something like "[1,3,4]@...". To handle multiple universe picks, use:
      tao_pointer_to_universes

    Parameters
    ----------
    ix_uni : int
        Index to the s.u(:) array If ix_uni is -1 -> u(s.global.default_universe) will be used.

    neg2_to_default : bool, optional
        i_uni = -2 (all universes) maps to the default uni? Default if False.

    Returns
    -------
    u : TaoUniverseStruct, optional
        Universe pointer. u will be nullified if there is an error and an error message will be printed.
    """

@overload
def tao_pointer_to_universe(string: str, neg2_to_default: bool | None = None) -> TaoPointerToUniverseStr:
    """
    Routine to set a pointer to a universe.

    This is an overloaded routine for the:
     tao_pointer_to_universe_int (ix_uni, neg2_to_default) result (u)
     tao_pointer_to_universe_str (string, neg2_to_default) result (u)

    Note: With a string argument, this routine can only handle single universe picks.
    That is, it cannot handlle something like "[1,3,4]@...". To handle multiple universe picks, use:
      tao_pointer_to_universes

    Parameters
    ----------
    string : str
        String in the form "<ix_uni>@..." or, if no "@" is present, u will point to the default universe.
        This parameter is an input/output and is modified in-place.
        As an output, string: String with universe prefix stripped off.

    neg2_to_default : bool, optional
        i_uni = -2 (all universes) maps to the default uni? Default if False.

    Returns
    -------
    string : str
        String in the form "<ix_uni>@..." or, if no "@" is present, u will point to the default universe.
        This parameter is an input/output and is modified in-place.
        As an output, string: String with universe prefix stripped off.

    u : TaoUniverseStruct, optional
        Universe pointer. u will be nullified if there is an error and an error message will be printed.
    """

class TaoPointerToUniverseStr:
    """tao_pointer_to_universe_str return type"""

    @property
    def u(self) -> _pybmad.TaoUniverseStruct | None: ...

    @property
    def string(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

class TaoPointerToUniverses:
    """tao_pointer_to_universes return type"""

    @property
    def unis(self) -> _pybmad.TaoUniversePointerStructAlloc1D: ...

    @property
    def err(self) -> bool: ...

    @property
    def name_out(self) -> str: ...

    @property
    def explicit_uni(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_pointer_to_universes(name_in: str, dflt_uni: int | None = None) -> TaoPointerToUniverses:
    """
    Wrapper for Fortran routine tao_pointer_to_universes

    Parameters
    ----------
    name_in : str
        data name with possible universe spec.

    dflt_uni : int, optional
        Default universe to use. Set to -1 if explicit universe is required.

    Returns
    -------
    unis : 1D array of TaoUniversePointerStruct
        Array of pointers to picked universes. The array will be resized if necessary.

    err : bool
        Set True if an error is detected.

    name_out : str, optional
        name_in without any "n@" beginning.

    explicit_uni : bool, optional
        Set True if name_in has explicit universe "n@" specification.
    """

def tao_pointer_to_var_in_lattice(var: _pybmad.TaoVarStruct, ix_uni: int, ele: _pybmad.EleStruct) -> bool:
    """
    Routine to add a pointer from a given Tao variable
    to the appropriate variable in a lattice.

    Parameters
    ----------
    var : TaoVarStruct
        Structure has the info of where to point.

    ix_uni : int
        the universe to use

    Returns
    -------
    err : bool
        Set True if there is an error. False otherwise.
    """

def tao_pointer_to_var_in_lattice2(var: _pybmad.TaoVarStruct, ix_uni: int) -> bool:
    """
    Routine to add a pointer from a given Tao variable
    to the appropriate variable in a lattice.

    Parameters
    ----------
    var : TaoVarStruct
        Structure has the info of where to point.

    ix_uni : int
        the universe to use

    Returns
    -------
    err : bool
        Set True if there is an error. False otherwise.
    """

def tao_print_command_line_info() -> None:
    """Wrapper for Fortran routine tao_print_command_line_info"""

def tao_print_vars(iu: int, ix_uni: int, show_good_opt_only: bool | None = None, tao_format: bool | None = None, v_array: _pybmad.TaoVarArrayStructArray1D | None = None) -> None:
    """
    Routine to print a list of set statements for the Bmad parameters controlled by the Tao variables.
    The set statements are Bmad lattice format compatible.

    When tao_format = True, the output is in the form "set variable <name> = <value>"
    so the file can be used as a Tao command file. If tao_format = False, the format
    is suitable for inclusion in a Bmad lattice file.

    Parameters
    ----------
    iu : int
        File unit number. 0 => print to the terminal.

    ix_uni : int
        Universe index. If zero print slave parameters for all universes. If non-zero, only print set statements
        for slave parameters of this universe. Ignored if tao_format = True.

    show_good_opt_only : bool, optional
        If True, only show slave parameters of variables used in optimization.

    tao_format : bool, optional
        Output format. Default False. See above.

    v_array : 1D array of TaoVarArrayStruct, optional
        Variable array. If present, restrict printing to parameters of these variables.
    """

def tao_ptc_normal_form(do_calc: bool, tao_lat: _pybmad.TaoLatticeStruct, ix_branch: int, rf_on: int | None = None) -> None:
    """
    Wrapper for Fortran routine tao_ptc_normal_form

    Parameters
    ----------
    do_calc : bool
        Set True to do the calculation.

    tao_lat : TaoLatticeStruct
        Lattice to work on.

    ix_branch : int
        Branch of lattice to work on.

    rf_on : int, optional
        RF state for calculation. yes$, no$, or maybe$ (default) maybe$ means that RF state in branch is used.
    """

def tao_python_cmd(input_str: str) -> None:
    """
    Wrapper for Fortran routine tao_python_cmd

    Parameters
    ----------
    input_str : str
        What to show.
    """

def tao_quiet_set(set: str) -> None:
    """
    Wrapper for Fortran routine tao_quiet_set

    Parameters
    ----------
    set : str
        True is silent running is wanted.
    """

def tao_rad_int_calc_needed(data_type: str, data_source: str) -> bool:
    """
    Wrapper for Fortran routine tao_rad_int_calc_needed

    Parameters
    ----------
    data_type : str

    data_source : str

    Returns
    -------
    do_rad_int : bool
    """

def tao_re_allocate_expression_info(info: _pybmad.TaoExpressionInfoStructAlloc1D, n: int, exact: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_re_allocate_expression_info

    Parameters
    ----------
    info : 1D array of TaoExpressionInfoStruct
        This parameter is an input/output and is modified in-place.
        As an output, info: Allocated array with size(re) >= n.

    n : int
        Size wanted.

    exact : bool, optional
        If present and False then the size of the output array is permitted to be larger than n. Default is True.
    """

def tao_re_associate_node_array(tree: _pybmad.TaoEvalNodeStruct, n: int, exact: bool | None = None) -> None:
    """
    Routine to resize the tree%node(:) array.

    Note: The data of the array is preserved but data at the end of the
    array will be lost if n is less than the original size of the array

    Parameters
    ----------
    tree : TaoEvalNodeStruct

    n : int
        Size wanted.

    exact : bool, optional
        Default is False. If False, the size of the output array is permitted to be larger than n.
    """

def tao_re_execute(string: str, err: bool) -> None:
    """Subroutine to execute a previous command."""

def tao_read_cmd(which: str, unis: str, file: str, silent: bool) -> None:
    """
    Wrapper for Fortran routine tao_read_cmd

    Parameters
    ----------
    which : str

    unis : str
        Universes to apply to

    file : str

    silent : bool
        Silent
    """

def tao_read_phase_space_index(name: str, ixc: int, print_err: bool | None = None) -> int:
    """
    Wrapper for Fortran routine tao_read_phase_space_index

    Parameters
    ----------
    name : str
        character array holding the index. Must be in the range 1-6.

    ixc : int
        location within <name> to evaluate index.

    print_err : bool, optional
        If present and False then do not print an error message

    Returns
    -------
    ix_ps : int
        Index at <name>(<ixc>:<ixc>). Returns 0 if bad index.
    """

def tao_regression_test(cmd_str: str) -> None:
    """
    Wrapper for Fortran routine tao_regression_test

    Parameters
    ----------
    cmd_str : str
    """

class TaoRemoveBlankCharacters:
    """tao_remove_blank_characters return type"""

    @property
    def str(self) -> str: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_remove_blank_characters(str: str) -> TaoRemoveBlankCharacters:
    """
    Wrapper for Fortran routine tao_remove_blank_characters

    Parameters
    ----------
    str : str
        Input string.
        This parameter is an input/output and is modified in-place.
        As an output, str: String with blank characters removed.

    Returns
    -------
    str : str
        Input string.
        This parameter is an input/output and is modified in-place.
        As an output, str: String with blank characters removed.
    """

def tao_run_cmd(which: str) -> bool:
    """
    Wrapper for Fortran routine tao_run_cmd

    Parameters
    ----------
    which : str
        which optimizer to use.

    Returns
    -------
    abort : bool
        Set True if the run was aborted by the user, an at minimum condition, a singular matrix condition, etc..
        False otherwise.
    """

def tao_scale_cmd(where: str, y_min_in: float, y_max_in: float, axis: str | None = None, include_wall: bool | None = None, gang: str | None = None, exact: bool | None = None, turn_autoscale_off: bool | None = None) -> None:
    r"""
    Routine to scale a plot.
    If y_min = y_max, the scales will be chosen to show all the data.

    Parameters
    ----------
    where : str
        Region to scale. Eg: "top:x"

    y_min_in : float
        Plot y-axis min value.

    y_max_in : float
        Plot y-axis max value.

    axis : str, optional
        'y', 'y2', or \'\' (both). Default = \'\'.

    include_wall : bool, optional
        Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
        include the building wall position will be included in determining the the scale.

    gang : str, optional
        'gang', 'nogang', \'\'. Default = \'\'.

    exact : bool, optional
        Exact plot y_max, y_min to correspond to y_min_in, y_max_in? Default is False. Only relavent when y_min_in
        /= y_max_in.

    turn_autoscale_off : bool, optional
        If True (default) then turn off plot.autoscale_y logical for all plots that are scaled.
    """

class TaoScaleGraph:
    """tao_scale_graph return type"""

    @property
    def y_range(self) -> list[float]: ...

    @property
    def y2_range(self) -> list[float]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_scale_graph(graph: _pybmad.TaoGraphStruct, y_min: float, y_max: float, axis: str | None = None, include_wall: bool | None = None) -> TaoScaleGraph:
    r"""
    Routine to scale the y-axis and/or y2-axis of a graph
    If y_min = y_max then autoscaling will be done and the particular value of y_min and y_max is ignored.
    Note: y_min/y_max is ignored if scaling y2-axis and graph%y2_mirrors_y = T.

    Parameters
    ----------
    graph : TaoGraphStruct
        Graph with axis/axes to be scaled.
        This parameter is an input/output and is modified in-place.
        As an output, graph: Graph with scaled axis/axes.

    y_min : float
        Axis [min, max] must cover [y_min, y_max] if not autoscaling.

    y_max : float
        Axis [min, max] must cover [y_min, y_max] if not autoscaling.

    axis : str, optional
        Axis to scale. \'\'   -> scale y and y2 (default). 'y'  -> scale y-axis. 'y2' -> scale y2-axis

    include_wall : bool, optional
        Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
        include the building wall position will be included in determining the the scale.

    Returns
    -------
    y_range : 1D array of float (shape: 2), optional
        Only used by tao_scale_plot when ganging graphs.

    y2_range : 1D array of float (shape: 2), optional
        Only used by tao_scale_plot when ganging graphs.
    """

def tao_scale_ping_data(u: _pybmad.TaoUniverseStruct) -> None:
    """
    Wrapper for Fortran routine tao_scale_ping_data

    Parameters
    ----------
    u : TaoUniverseStruct
    """

def tao_scale_plot(plot: _pybmad.TaoPlotStruct, y_min_in: float, y_max_in: float, axis: str | None = None, include_wall: bool | None = None, gang: str | None = None, skip_lat_layout: bool | None = None) -> None:
    r"""
    Routine to scale the y-axis and/or y2-axis of the graphs of the plot.
    If y_min_in = y_max_in then autoscaling will be done and the particular value
    of y_min_in and y_max_in is ignored.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot with graphs to be scaled.
        This parameter is an input/output and is modified in-place.
        As an output, plot: Plot with scaled graphs.

    y_min_in : float
        Axis [min, max] must cover [y_min_in, y_max_in] if not autoscaling.

    y_max_in : float
        Axis [min, max] must cover [y_min_in, y_max_in] if not autoscaling.

    axis : str, optional
        Axis to scale. \'\'   -> scale y and y2 (default). 'y'  -> scale y-axis. 'y2' -> scale y2-axis

    include_wall : bool, optional
        Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
        include the building wall position will be included in determining the the scale.

    gang : str, optional
        If autoscale then make all graph y-axes the same and/or make all y2-axes the same? \'\'        -> (default)
        Use setting of plot.autoscale_gang_y 'gang'    -> Gang graphs. 'nogang'  -> Do not gang graphs.

    skip_lat_layout : bool, optional
        If True, skip scaling any lat_layout graphs. Default is false.
    """

def tao_scratch_values_calc(ele_ref: _pybmad.EleStruct, ele_start: _pybmad.EleStruct, ele: _pybmad.EleStruct, datum: _pybmad.TaoDataStruct, branch: _pybmad.BranchStruct, orbit: _pybmad.CoordStructArray1D) -> None:
    """
    Wrapper for Fortran routine tao_scratch_values_calc

    Parameters
    ----------
    ele_ref : EleStruct

    ele_start : EleStruct

    ele : EleStruct

    datum : TaoDataStruct

    branch : BranchStruct

    orbit : 1D array of CoordStruct
    """

def tao_set_beam_cmd(who: str, value_str: str, branch_str: str) -> None:
    r"""
    Routine to set various beam parameters.

    Parameters
    ----------
    who : str
        which parameter to set.

    value_str : str
        Value to set to.

    branch_str : str
        Branch to use. \'\' => branch 0.
    """

def tao_set_beam_init_cmd(who: str, value_str: str, branch_str: str) -> None:
    r"""
    Routine to set beam_init variables

    Parameters
    ----------
    who : str
        which beam_init variable to set

    value_str : str
        Value to set to.

    branch_str : str
        Branch to use. \'\' => branch 0
    """

def tao_set_bmad_com_cmd(who: str, value_str: str) -> None:
    """
    Routine to set bmad_com variables

    Parameters
    ----------
    who : str
        which bmad_com variable to set

    value_str : str
        Value to set to.
    """

def tao_set_branch_cmd(branch_str: str, component_str: str, value_str: str) -> None:
    """
    Routine to set lattice branch values.

    Parameters
    ----------
    branch_str : str
        Which branch to set.

    component_str : str
        Which branch parameter to set.

    value_str : str
        What value to set it to.
    """

def tao_set_calculate_cmd(switch_: str | None = None) -> None:
    """Toggles off lattice calc and plotting."""

def tao_set_curve_cmd(curve_name: str, component: str, value_str: str) -> None:
    """
    Routine to set var values.

    Parameters
    ----------
    curve_name : str
        Which curve to set.

    component : str
        Which component to set.

    value_str : str
        What value to set it to.
    """

def tao_set_curve_invalid(curve: _pybmad.TaoCurveStruct, why_invalid: str, print_err: bool | None = None) -> None:
    """
    Routine to set curve%valid to False.

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve to set.
        This parameter is an input/output and is modified in-place.
        As an output, curve: Curve properly set.

    why_invalid : str
        Invalid information.

    print_err : bool, optional
        If present and True then also print an error message.
    """

def tao_set_data_cmd(who_str: str, value_str: str, silent: bool | None = None) -> None:
    """
    Routine to set data values.

    Parameters
    ----------
    who_str : str
        Which data component(s) to set.

    value_str : str
        What value to set it to.
    """

def tao_set_data_useit_opt(data: _pybmad.TaoDataStructArray1D | None = None) -> None:
    """
    Wrapper for Fortran routine tao_set_data_useit_opt

    Parameters
    ----------
    data : 1D array of TaoDataStruct, optional
        Data to work on. Default is all data in all universes.
    """

def tao_set_default_cmd(who_str: str, value_str: str) -> None:
    """
    Routine to set default values.

    Parameters
    ----------
    who_str : str
        Which default component(s) to set.

    value_str : str
        What value to set it to.
    """

def tao_set_drawing_cmd(drawing: _pybmad.TaoDrawingStruct, component: str, value_str: str) -> None:
    """
    Routine to set floor_plan and lat_layout parameters.

    Parameters
    ----------
    drawing : TaoDrawingStruct
        s.plot_page.floor_plan or s.plot_page.lat_layout.

    component : str
        Which shape component to set.

    value_str : str
        Value to set to.
    """

def tao_set_dynamic_aperture_cmd(who: str, value_str: str) -> None:
    """
    Sets dynamic aperture parameters.

    Parameters
    ----------
    who : str
        which parameter to set.

    value_str : str
        Value to set to.
    """

def tao_set_elements_cmd(ele_list: str, attribute: str, value: str, update: bool) -> None:
    """
    Sets element parameters.

    Parameters
    ----------
    ele_list : str
        which elements.

    attribute : str
        Attribute to set.

    value : str
        Value to set.
    """

def tao_set_flags_for_changed_attribute(u: _pybmad.TaoUniverseStruct, ele_name: str, ele_ptr: _pybmad.EleStruct | None = None, val_ptr: _pybmad.AllPointerStruct | None = None, who: str | None = None) -> None:
    """
    Wrapper for Fortran routine tao_set_flags_for_changed_attribute

    Parameters
    ----------
    u : TaoUniverseStruct
        Universe containing the lattice.

    ele_name : str
        Associated "element" of the changed parameter.

    ele_ptr : EleStruct, optional
        Pointer to the element. May be null, for example, if ele_name = "PARTICLE_START".

    val_ptr : AllPointerStruct, optional
        optional: Pointer to the attribute that was changed. Must be present if ele_ptr is present.

    who : str, optional
        Name of changed attribute. Only used with PARTICLE_START.
    """

def tao_set_floor_plan_axis_label(graph: _pybmad.TaoGraphStruct, axis_in: _pybmad.QpAxisStruct, axis_out: _pybmad.QpAxisStruct, which: str) -> None:
    """
    Wrapper for Fortran routine tao_set_floor_plan_axis_label

    Parameters
    ----------
    graph : TaoGraphStruct

    axis_in : QpAxisStruct

    axis_out : QpAxisStruct

    which : str
    """

def tao_set_geodesic_lm_cmd(who: str, value_str: str) -> None:
    """
    Routine to set geodesic_lm variables

    Parameters
    ----------
    who : str
        which geodesic_lm variable to set

    value_str : str
        Value to set to.
    """

def tao_set_global_cmd(who: str, value_str: str) -> None:
    """
    Routine to set global variables

    Parameters
    ----------
    who : str
        which global variable to set

    value_str : str
        Value to set to.
    """

def tao_set_graph_cmd(graph_name: str, component: str, value_str: str) -> None:
    """
    Routine to set var values.

    Parameters
    ----------
    graph_name : str
        Which graph to set.

    component : str
        Which component to set.

    value_str : str
        What value to set it to.
    """

class TaoSetIntegerValue:
    """tao_set_integer_value return type"""

    @property
    def var(self) -> int: ...

    @property
    def error(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_set_integer_value(var_str: str, value_str: str, min_val: int | None = None, max_val: int | None = None, print_err: bool | None = None) -> TaoSetIntegerValue:
    """
    Subroutine to read and set the value of an integer varialbe.

    If the value is out of the range [min_val, max_val] then an error message will
    be generated and the variable will not be set.

    Parameters
    ----------
    var_str : str
        Used for error messages.

    value_str : str
        String with encoded value.

    min_val : int, optional
        Minimum value.

    max_val : int, optional
        Maximum value.

    print_err : bool, optional
        If True, print error message. Default is true

    Returns
    -------
    var : int
        Variable to set.

    error : bool
        Set True on an error. False otherwise.
    """

def tao_set_invalid(datum: _pybmad.TaoDataStruct, message: str, exterminate: bool | None = None, err_level: int | None = None, print_err: bool | None = None) -> str:
    """
    Wrapper for Fortran routine tao_set_invalid

    Parameters
    ----------
    datum : TaoDataStruct
        Bad datum.

    message : str
        Error message.

    exterminate : bool, optional
        Default is False. If True, set datum.exists to False so that Tao will ignore this datum from now on.

    err_level : int, optional
        s_error$ (default), s_warn$, etc.

    print_err : bool, optional
        Default is True. If False, do not print an error message.

    Returns
    -------
    why_invalid : str, optional
        Set to message if present.
    """

def tao_set_key_cmd(key_str: str, cmd_str: str) -> None:
    """
    Associates a command with a key press for single mode.

    Parameters
    ----------
    key_str : str
        keyboard key.

    cmd_str : str
        Command associated with key.
    """

def tao_set_lattice_cmd(dest_lat: str, source_lat: str) -> None:
    """
    Sets a lattice equal to another. This will also update the data structs

    Parameters
    ----------
    dest_lat : str
        Maybe: 'model', 'design', or 'base' with optional '@n' at beginning to indicate the universe

    source_lat : str
        Maybe: 'model', 'design', or 'base'
    """

class TaoSetLogicalValue:
    """tao_set_logical_value return type"""

    @property
    def var(self) -> bool: ...

    @property
    def error(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_set_logical_value(var_str: str, value_str: str) -> TaoSetLogicalValue:
    """
    Subroutine to read and set the value of an logical varialbe.

    If the value is out of the range [min_val, max_val] then an error message will
    be generated and the variable will not be set.

    Parameters
    ----------
    var_str : str
        Used for error messages.

    value_str : str
        String with encoded value.

    Returns
    -------
    var : bool
        Variable to set.

    error : bool
        Set True on an error. False otherwise.
    """

def tao_set_openmp_n_threads(n_threads: int) -> None:
    """
    Routine to set OpenMP thread count.  Errors if OpenMP is not available.

    Parameters
    ----------
    n_threads : int
        Number of threads.
    """

def tao_set_opt_vars(var_vec: _pybmad.RealArray1D, print_limit_warning: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_set_opt_vars

    Parameters
    ----------
    var_vec : 1D array of float
        Vector of variables.

    print_limit_warning : bool, optional
        Print a warning if the value is past the variable's limits. Default is True.
    """

def tao_set_opti_de_param_cmd(who: str, value_str: str) -> None:
    """
    Routine to set opti_de_param variables

    Parameters
    ----------
    who : str
        which opti_de_param variable to set

    value_str : str
        Value to set to.
    """

def tao_set_particle_start_cmd(who: str, value_str: str) -> None:
    """
    Routine to set particle_start variables.

    Parameters
    ----------
    who : str
        which particle_start variable to set

    value_str : str
        Value to set to.
    """

def tao_set_plot_cmd(plot_name: str, component: str, value_str: str) -> None:
    """
    Routine to set plot parameters.

    Parameters
    ----------
    plot_name : str
        Which plot to set.

    component : str
        Which component to set.

    value_str : str
        What value to set it to.
    """

def tao_set_plot_page_cmd(component: str, value_str: str, value_str2: str | None = None) -> None:
    """
    Set various aspects of the plotting window

    Parameters
    ----------
    component : str
        Which component to set.

    value_str : str
        What value to set to.

    value_str2 : str, optional
        2nd value if component is an array.
    """

def tao_set_plotting(plot_input: _pybmad.TaoPlotPageInput, plot_page: _pybmad.TaoPlotPageStruct, use_cmd_line_geom: bool, reverse: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_set_plotting

    Parameters
    ----------
    plot_input : TaoPlotPageInput

    plot_page : TaoPlotPageStruct

    use_cmd_line_geom : bool

    reverse : bool, optional
    """

def tao_set_ptc_com_cmd(who: str, value_str: str) -> None:
    """
    Routine to set ptc_com variables

    Parameters
    ----------
    who : str
        which ptc_com variable to set

    value_str : str
        Value to set to.
    """

class TaoSetQpAxisStruct:
    """tao_set_qp_axis_struct return type"""

    @property
    def error(self) -> bool: ...

    @property
    def ix_uni(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_set_qp_axis_struct(qp_axis_name: str, component: str, qp_axis: _pybmad.QpAxisStruct, value: str) -> TaoSetQpAxisStruct:
    """
    Routine to set qp_axis_names of a qp_axis_struct.

    Parameters
    ----------
    qp_axis_name : str
        qp_axis name. Used for error messages.

    component : str
        qp_axis component name.

    qp_axis : QpAxisStruct
        qp_axis_struct with component to modify
        This parameter is an input/output and is modified in-place.
        As an output, qp_axis: qp_axis_struct with changed component value.

    value : str
        Component value.

    Returns
    -------
    error : bool
        Set true if there is an error. False otherwise.

    ix_uni : int, optional
        Tao universe number in case the value depends upon a parameter of a particular universe.
    """

class TaoSetQpPointStruct:
    """tao_set_qp_point_struct return type"""

    @property
    def error(self) -> bool: ...

    @property
    def ix_uni(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_set_qp_point_struct(qp_point_name: str, component: str, qp_point: _pybmad.QpPointStruct, value: str) -> TaoSetQpPointStruct:
    """
    Routine to set qp_point_names of a qp_point_struct.

    Parameters
    ----------
    qp_point_name : str
        qp_point name. Used for error messages.

    component : str
        qp_point component name.

    qp_point : QpPointStruct
        qp_point_struct with component to modify
        This parameter is an input/output and is modified in-place.
        As an output, qp_point: qp_point_struct with changed component value.

    value : str
        Component value.

    Returns
    -------
    error : bool
        Set true if there is an error. False otherwise.

    ix_uni : int, optional
        Tao universe number in case the value depends upon a parameter of a particular universe.
    """

class TaoSetQpRectStruct:
    """tao_set_qp_rect_struct return type"""

    @property
    def error(self) -> bool: ...

    @property
    def ix_uni(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_set_qp_rect_struct(qp_rect_name: str, component: str, qp_rect: _pybmad.QpRectStruct, value: str) -> TaoSetQpRectStruct:
    """
    Routine to set qp_rect_names of a qp_rect_struct.

    Parameters
    ----------
    qp_rect_name : str
        qp_rect name. Used for error messages.

    component : str
        qp_rect component name.

    qp_rect : QpRectStruct
        qp_rect_struct with component to modify
        This parameter is an input/output and is modified in-place.
        As an output, qp_rect: qp_rect_struct with changed component value.

    value : str
        Component value.

    Returns
    -------
    error : bool
        Set true if there is an error. False otherwise.

    ix_uni : int, optional
        Tao universe number in case the value depends upon a parameter of a particular universe.
    """

def tao_set_ran_state_cmd(state_string: str) -> None:
    """
    Sets the random number generator state.

    Parameters
    ----------
    state_string : str
        Encoded random number state.
    """

class TaoSetRealValue:
    """tao_set_real_value return type"""

    @property
    def var(self) -> float: ...

    @property
    def error(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_set_real_value(var_str: str, value_str: str, min_val: float | None = None, max_val: float | None = None, dflt_uni: int | None = None) -> TaoSetRealValue:
    """
    Subroutine to read and set the value of a real variable.

    If the value is out of the range [min_val, max_val] then an error message will
    be generated and the variable will not be set.

    Parameters
    ----------
    var_str : str
        Used for error messages.

    value_str : str
        String with encoded value.

    min_val : float, optional
        Minimum value.

    max_val : float, optional
        Maximum value.

    dflt_uni : int, optional
        Default universe used to evaluate parameters.

    Returns
    -------
    var : float
        Variable to set.

    error : bool
        Set True on an error. False otherwise.
    """

def tao_set_region_cmd(region_name: str, component: str, value_str: str) -> None:
    """
    Routine to set region parameters.

    Parameters
    ----------
    region_name : str
        Which region to set.

    component : str
        Which component to set.

    value_str : str
        What value to set it to.
    """

def tao_set_space_charge_com_cmd(who: str, value_str: str) -> None:
    """
    Routine to set space_charge_com variables

    Parameters
    ----------
    who : str
        which space_charge_com variable to set

    value_str : str
        Value to set to.
    """

def tao_set_symbolic_number_cmd(sym_str: str, num_str: str | None = None, val: float | None = None) -> None:
    """
    Associates a given symbol with a given number.
    Note: Either num_str or val argument must be present.

    Parameters
    ----------
    sym_str : str
        Symbol.

    num_str : str, optional
        Symbol value expression.

    val : float, optional
        Value of symbol
    """

def tao_set_tune_cmd(branch_str: str, mask_str: str, print_list: bool, qa_str: str, qb_str: str, delta_input: bool) -> None:
    """
    Routine to set the transverse tunes.

    Parameters
    ----------
    branch_str : str
        List of branches to apply tune set to.

    mask_str : str
        List of quadrupoles to veto.

    print_list : bool
        If True, print a list of elements varied and coefficients.

    qa_str : str
        Expression for Qa tune.

    qb_str : str
        Expression for Qb tune.

    delta_input : bool
        If true then qa_str and qb_str are deltas from present tune.
    """

def tao_set_universe_cmd(uni: str, who: str, what: str) -> None:
    """
    Sets a universe on or off, or sets the recalculate or twiss_calc logicals, etc.

    Parameters
    ----------
    uni : str
        which universe; 0 => current viewed universe

    who : str
        "on", "off", "recalculate", "dynamic_aperture_calc", "one_turn_map_calc", or "twiss_calc"

    what : str
        "on" or "off" for who = "dynamic_aperture_calc", "one_turn_map_calc" or "twiss_calc".
    """

def tao_set_var_cmd(var_str: str, value_str: str) -> None:
    """
    Routine to set var values.

    Parameters
    ----------
    var_str : str
        Which var name to set.

    value_str : str
        What value to set it to.
    """

def tao_set_var_model_value(var: _pybmad.TaoVarStruct, value: float, print_limit_warning: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_set_var_model_value

    Parameters
    ----------
    var : TaoVarStruct
        Variable to set

    value : float
        Value to set to

    print_limit_warning : bool, optional
        Print a warning if the value is past the variable's limits. Default is True.
    """

def tao_set_var_useit_opt() -> None:
    """Wrapper for Fortran routine tao_set_var_useit_opt"""

def tao_set_wave_cmd(who: str, value_str: str) -> bool:
    """
    Routine to set wave variables

    Parameters
    ----------
    who : str
        which wave variable to set

    value_str : str
        Value to set to.

    Returns
    -------
    err : bool
        Set True if there is an error. False otherwise.
    """

def tao_set_z_tune_cmd(branch_str: str, q_str: str, delta_input: bool) -> None:
    """
    Routine to set the z-tune.

    Parameters
    ----------
    branch_str : str
        List of branches to apply tune set to.

    q_str : str
        Expression for Qc tune.

    delta_input : bool
        If true then qa_str and qb_str are deltas from present tune.
    """

def tao_setup_key_table() -> None:
    """Wrapper for Fortran routine tao_setup_key_table"""

def tao_shape_init(shape: _pybmad.TaoEleShapeStruct, print_err: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine tao_shape_init

    Parameters
    ----------
    shape : TaoEleShapeStruct
        Shape

    print_err : bool, optional
        If True then print an error message if there is a problem. Default is True.

    Returns
    -------
    err : bool
        Set true if there is a problem translating the element class.
    """

def tao_show_cmd(what: str) -> None:
    """
    Wrapper for Fortran routine tao_show_cmd

    Parameters
    ----------
    what : str
        What to show.
    """

def tao_show_constraints(iunit: int, form: str) -> None:
    """
    Routine to show a list of datums and variables and how they contribute to the merit function.

    Parameters
    ----------
    iunit : int
        File unit to write to. 0 => print to the terminal.

    form : str
        What to output: 'ALL'   -> All datums and variables. 'TOP10' -> Top datums and variables that contribute
        to the merit function.
    """

class TaoShowThis:
    """tao_show_this return type"""

    @property
    def result_id(self) -> str: ...

    @property
    def lines(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def nl(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_show_this(what: str) -> TaoShowThis:
    """
    Wrapper for Fortran routine tao_show_this

    Parameters
    ----------
    what : str
        What to show.

    Returns
    -------
    result_id : str
        ID string idendifying what was shown. Used by the Tao GUI.

    lines : 1D array of str
        Output.

    nl : int
        Number of lines(:).
    """

def tao_single_mode(char_: str) -> None:
    """
    Wrapper for Fortran routine tao_single_mode

    Parameters
    ----------
    """

def tao_single_track(tao_lat: _pybmad.TaoLatticeStruct, ix_branch: int, print_err: bool | None = None) -> bool:
    """
    Routine to track a single particle and calculate lattice functions through a lattice.

    Parameters
    ----------
    tao_lat : TaoLatticeStruct
        Structure containing the lattice.

    ix_branch : int
        Branch index to track through.

    print_err : bool, optional
        Default False. Print error messages if, eg, lattice is unstable?

    Returns
    -------
    calc_ok : bool
        Set True if there were no problems, False otherwise.
    """

def tao_spin_matrices_calc_needed(data_type: str, data_source: str) -> bool:
    """
    Wrapper for Fortran routine tao_spin_matrices_calc_needed

    Parameters
    ----------
    data_type : str

    data_source : str

    Returns
    -------
    do_calc : bool
    """

def tao_spin_tracking_turn_on() -> None:
    """Wrapper for Fortran routine tao_spin_tracking_turn_on"""

class TaoSplitComponent:
    """tao_split_component return type"""

    @property
    def comp(self) -> _pybmad.TaoDataVarComponentStructAlloc1D: ...

    @property
    def err(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_split_component(comp_str: str) -> TaoSplitComponent:
    """
    Wrapper for Fortran routine tao_split_component

    Parameters
    ----------
    comp_str : str
        Components. EG: 'meas - design'

    Returns
    -------
    comp : 1D array of TaoDataVarComponentStruct
        Array of individual components.

    err : bool
        Set True if there is an error, False otherwise.
    """

def tao_srdt_calc_needed(data_type: str, data_source: str) -> int:
    """
    Wrapper for Fortran routine tao_srdt_calc_needed

    Parameters
    ----------
    data_type : str

    data_source : str

    Returns
    -------
    do_srdt : int
    """

class TaoSubinUniNumber:
    """tao_subin_uni_number return type"""

    @property
    def name_out(self) -> str: ...

    @property
    def ok(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_subin_uni_number(name_in: str, ix_uni: int) -> TaoSubinUniNumber:
    """
    Wrapper for Fortran routine tao_subin_uni_number

    Parameters
    ----------
    name_in : str
        Input name with "#" character

    ix_uni : int
        Universe index.

    Returns
    -------
    name_out : str
        Output name.

    ok : bool
        False if multiple universes and no "#" in name_in. True otherwise.
    """

def tao_svd_optimizer() -> bool:
    """
    Routine to minimize the merit function using svd.

    Returns
    -------
    abort : bool
        Set True if svd step increases the merit function.
    """

def tao_symbol_import_from_lat(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine tao_symbol_import_from_lat

    Parameters
    ----------
    lat : LatStruct
    """

def tao_taper_cmd(except_: str, uni_names: str) -> None:
    """
    Wrapper for Fortran routine tao_taper_cmd

    Parameters
    ----------
    except : str
        List of elements not to vary.

    uni_names : str
        Universes to taper.
    """

def tao_to_change_number(num_str: str, n_size: int, change_number: _pybmad.RealAlloc1D, abs_or_rel: str, err: bool) -> None:
    """
    Wrapper for Fortran routine tao_to_change_number

    Parameters
    ----------
    num_str : str

    n_size : int

    change_number : 1D array of float

    abs_or_rel : str

    err : bool
    """

def tao_to_int(str: str, i_int: int, err: bool) -> None:
    """
    Converts a string to an integer

    If the string str is blank then i_int = 0
    """

class TaoToPhaseAndCouplingReading:
    """tao_to_phase_and_coupling_reading return type"""

    @property
    def bpm_data(self) -> _pybmad.BpmPhaseCouplingStruct: ...

    @property
    def valid_value(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_to_phase_and_coupling_reading(ele: _pybmad.EleStruct, why_invalid: str, datum: _pybmad.TaoDataStruct) -> TaoToPhaseAndCouplingReading:
    """
    Buffer routine for to_phase_and_coupling_reading.

    Parameters
    ----------
    ele : EleStruct
        The monitor.

    Returns
    -------
    bpm_data : BpmPhaseCouplingStruct
        Monitor values

    valid_value : bool
        Valid data value?
    """

class TaoToReal:
    """tao_to_real return type"""

    @property
    def value(self) -> float: ...

    @property
    def err_flag(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_to_real(expression: str) -> TaoToReal:
    """
    Wrapper for Fortran routine tao_to_real

    Parameters
    ----------
    expression : str
        arithmetic expression

    Returns
    -------
    value : float
        Value of arithmetic expression.

    err_flag : bool
        TRUE on error.
    """

def tao_to_top10(top10: _pybmad.TaoTop10StructArray1D, value: float, name: str, c_index: int, order: str) -> None:
    """
    Routine to order the largest contributors to the merit function in
    a list. Call this routine for each contributor.

    Note: Before first calling this routine set:
      top10(:)%valid = .false.

    Parameters
    ----------
    top10 : 1D array of TaoTop10Struct
        List of top contributors. Note that the list is not limited to 10 entries.

    value : float
        value of the contributor.

    name : str
        Name of the contributor..

    c_index : int
        Index of the contributor.

    order : str
        Ordering of the list. Possibilities are:
    """

def tao_too_many_particles_lost(beam: _pybmad.BeamStruct) -> bool:
    """
    Wrapper for Fortran routine tao_too_many_particles_lost

    Parameters
    ----------
    beam : BeamStruct

    Returns
    -------
    no_beam : bool
    """

def tao_top10_derivative_print() -> None:
    """Routine to print out the top10 contributors to the merit function."""

def tao_top10_merit_categories_print(iunit: int) -> None:
    """
    Routine to print the top data and variable categories that contribute to
    the merit function.

    Parameters
    ----------
    iunit : int
        File unit to write to. 0 => print to the terminal.
    """

def tao_top_level(command: str | None = None) -> int:
    """
    Wrapper for Fortran routine tao_top_level

    Parameters
    ----------
    command : str, optional
        Tao command string. If present, getting user input from the terminal is bypassed. This is used when
        interfacing to Python.

    Returns
    -------
    errcode : int, optional
        Return error code: 0 => OK, Not 0 => Err.
    """

class TaoTrackingEleIndex:
    """tao_tracking_ele_index return type"""

    @property
    def ix_branch(self) -> int: ...

    @property
    def ix_ele(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def tao_tracking_ele_index(ele: _pybmad.EleStruct, datum: _pybmad.TaoDataStruct) -> TaoTrackingEleIndex:
    """
    Routine to return the index in the tracking part of a lattice that corresponds to ele.

    Parameters
    ----------
    ele : EleStruct
        Lattice element.

    datum : TaoDataStruct
        Datum

    Returns
    -------
    ix_ele : int
        Element index associated with ele.

    ix_branch : int, optional
        Lattice branch associated with element
    """

def tao_turn_on_special_calcs_if_needed_for_plotting() -> None:
    """
    Wrapper for Fortran routine tao_turn_on_special_calcs_if_needed_for_plotting
    """

def tao_type_expression_tree(tree: _pybmad.TaoEvalNodeStruct, indent: int | None = None) -> None:
    """
    Routine to print an expression tree in tree form.
    Good for debugging.

    Parameters
    ----------
    tree : TaoEvalNodeStruct
        Tree to print.

    indent : int, optional
        Initial indent. Default is zero.
    """

def tao_uni_atsign_index(string: str) -> int:
    """
    Routine to return the index of an atsign ("@") character in a string if the atsign is
    being used as a separator between a universe spec and the rest of the string.

    For example:
      string = "[1:3]@orbit.x[5] => ix_amp = 6
      string = "orbit.x[5@0.2]   => ix_amp = 0 (no universe "@" present)

    Parameters
    ----------
    string : str
        String to parse

    Returns
    -------
    ix_amp : int
        Index of universe "@". Set to zero if no universe "@" found.
    """

def tao_universe_index(i_uni: int, neg2_to_default: bool | None = None) -> int:
    """
    Wrapper for Fortran routine tao_universe_index

    Parameters
    ----------
    i_uni : int
        Nominal universe number.

    neg2_to_default : bool, optional
        i_uni = -2 (all universes) maps to the default uni? Default if False.

    Returns
    -------
    i_this_uni : int
        Universe number.
    """

def tao_use_data(action: str, data_name: str) -> None:
    """
    Wrapper for Fortran routine tao_use_data

    Parameters
    ----------
    action : str
        veto, use or restore

    data_name : str
        the selected data
    """

def tao_use_var(action: str, var_name: str) -> None:
    """
    Wrapper for Fortran routine tao_use_var

    Parameters
    ----------
    action : str
        'use', 'veto', or 'restore'

    var_name : str
        the selected variable name or all
    """

def tao_user_is_terminating_optimization() -> bool:
    """
    Wrapper for Fortran routine tao_user_is_terminating_optimization

    Returns
    -------
    is_terminating : bool
        Set True of '.' is detected. False otherwise.
    """

def tao_var1_name(var: _pybmad.TaoVarStruct) -> str:
    """
    Wrapper for Fortran routine tao_var1_name

    Parameters
    ----------
    var : TaoVarStruct
        Variable

    Returns
    -------
    var1_name : str
        Appropriate name.
    """

def tao_var_attrib_name(var: _pybmad.TaoVarStruct) -> str:
    """
    Wrapper for Fortran routine tao_var_attrib_name

    Parameters
    ----------
    var : TaoVarStruct
        Variable

    Returns
    -------
    var_attrib_name : str
        Attribute list.
    """

def tao_var_check(eles: _pybmad.ElePointerStructAlloc1D, attribute: str, silent: bool) -> None:
    """
    Wrapper for Fortran routine tao_var_check

    Parameters
    ----------
    eles : 1D array of ElePointerStruct
        Array of elements which have a changed attribute.

    attribute : str
        Name of attribute changed.

    silent : bool
        If True and the problem can be fixed, do not issue an error message.
    """

def tao_var_repoint() -> None:
    """Wrapper for Fortran routine tao_var_repoint"""

def tao_var_show_use(v1_var: _pybmad.TaoV1VarStruct, lines: _pybmad.CharacterAlloc1D | None = None, nl: int | None = None) -> None:
    """
    Wrapper for Fortran routine tao_var_show_use

    Parameters
    ----------
    v1_var : TaoV1VarStruct
        tao_v1_var_struct

    lines : 1D array of str, optional

    nl : int, optional
    """

def tao_var_target_calc() -> None:
    """Wrapper for Fortran routine tao_var_target_calc"""

def tao_var_useit_plot_calc(graph: _pybmad.TaoGraphStruct, var: _pybmad.TaoVarStructArray1D) -> None:
    """
    Wrapper for Fortran routine tao_var_useit_plot_calc

    Parameters
    ----------
    graph : TaoGraphStruct

    var : 1D array of TaoVarStruct
    """

def tao_var_write(out_file: str, show_good_opt_only: bool | None = None, tao_format: bool | None = None) -> None:
    r"""
    Routine to write the optimized variables. One file will be created for each universe.
    The created file will have three sections:
      1) The variable values
      2) The list of constraints.
      3) A list of the top 10 constraints.
    If out_file = \'\' the information will be dumped to the terminal.
    In this case, only the variable values will be printed.

    When tao_format = True, the output is in the form "set variable <name> = <value>"
    so the file can be used as a Tao command file. If tao_format = False, the format
    is suitable for inclusion in a Bmad lattice file.

    Parameters
    ----------
    out_file : str
        Name of output file. If blank. Ouput to the terminal.

    show_good_opt_only : bool, optional
        Write only the variables used in the optimization? Default is False.

    tao_format : bool, optional
        Output format. Default False. See above.
    """

def tao_veto_vars_with_zero_dmodel() -> None:
    """
    Routine to veto all variables with zero effect on data used in the merit function.
    """

def tao_wave_analysis(plot: _pybmad.TaoPlotStruct) -> None:
    """
    Routine to do a wave anaysis.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot region setup by tao_wave_cmd.
        This parameter is an input/output and is modified in-place.
        As an output, plot: Plot with wave analysis curves.
    """

def tao_wave_cmd(curve_name: str, plot_place: str, err_flag: bool) -> None:
    """
    Routine to do the initial setup for wave plotting.
    The wave analysis is done by the routine tao_wave_analysis.

    Parameters
    ----------
    curve_name : str
        Character(*) curve for wave analysis.

    plot_place : str
        Character(*) place on plot page to put the wave plot.
    """

def tao_wave_fit(curve: _pybmad.TaoCurveStruct, ix1: int, n_dat: int, coef: _pybmad.RealArray1D, rms: _pybmad.RealArray1D, f1: _pybmad.RealArray1D, f2: _pybmad.RealArray1D | None = None, f3: _pybmad.RealArray1D | None = None, f4: _pybmad.RealArray1D | None = None) -> None:
    """
    Routine for fitting the curve data to up to four functions using a least squares
    SVD fit.

    Parameters
    ----------
    curve : TaoCurveStruct
        Curve containing the data.

    ix1 : int
        Index of first point in the data array.

    n_dat : int
        Number of data points.

    coef : 1D array of float
        Fit coefficients.

    rms : 1D array of float
        Variances with rms(n_func+1) = sqrt(chi^2/n_dat).

    f1 : 1D array of float
        First fit function.

    f2 : 1D array of float, optional
        Second fit function.

    f3 : 1D array of float, optional
        third fit function.

    f4 : 1D array of float, optional
        fourth fit function.
    """

def tao_write_cmd(what: str) -> None:
    """
    Wrapper for Fortran routine tao_write_cmd

    Parameters
    ----------
    what : str
        What to output. See the code for more details.
    """

def tao_write_lines(iunit: int, line: _pybmad.CharacterAlloc1D) -> None:
    """
    Subroutine to write out a series of lines to a file or to the terminal.
    It is assumed that any file has already been opened.

    Parameters
    ----------
    iunit : int
        File unit to write to. 0 => print to the terminal.

    line : 1D array of str
        A series of lines.
    """

def tao_x_axis_cmd(where: str, what: str) -> None:
    """
    Wrapper for Fortran routine tao_x_axis_cmd

    Parameters
    ----------
    where : str
        Region to axis. Eg: "top"

    what : str
        "s" or "index"
    """

def tao_x_scale_cmd(where: str, x_min_in: float, x_max_in: float, include_wall: bool | None = None, gang: str | None = None, exact: bool | None = None, turn_autoscale_off: bool | None = None) -> bool:
    r"""
    Routine to scale a plot. If x_min = x_max
    Then the scales will be chosen to show all the data.

    Parameters
    ----------
    where : str
        Region to scale. Eg: "top"

    x_min_in : float
        Plot x-axis min value.

    x_max_in : float
        Plot x-axis max value.

    include_wall : bool, optional
        Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
        include the building wall position will be included in determining the the scale.

    gang : str, optional
        'gang', 'nogang', \'\'. Default = \'\'.

    exact : bool, optional
        Exact plot y_max, y_min to correspond to y_min_in, y_max_in? Default is False. Only relavent when y_min_in
        /= y_max_in.

    turn_autoscale_off : bool, optional
        If True (default) then turn off plot.autoscale_x logical for all plots that are scaled.

    Returns
    -------
    err : bool
        Set to True if the plot cannot be found. False otherwise.
    """

def tao_x_scale_graph(graph: _pybmad.TaoGraphStruct, x_min: float, x_max: float, include_wall: bool | None = None, have_scaled: bool | None = None) -> None:
    """
    Wrapper for Fortran routine tao_x_scale_graph

    Parameters
    ----------
    graph : TaoGraphStruct

    x_min : float

    x_max : float

    include_wall : bool, optional

    have_scaled : bool, optional
    """

def tao_x_scale_plot(plot: _pybmad.TaoPlotStruct, x_min_in: float, x_max_in: float, include_wall: bool | None = None, gang: str | None = None) -> bool:
    r"""
    Routine to scale a plot. If x_min = x_max
    Then the scales will be chosen to show all the data.

    Parameters
    ----------
    plot : TaoPlotStruct
        Plot to scale. Eg: "top"

    x_min_in : float
        Plot x-axis min value.

    x_max_in : float
        Plot x-axis max value.

    include_wall : bool, optional
        Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in. If present and True
        include the building wall position will be included in determining the the scale.

    gang : str, optional
        'gang', 'nogang', \'\'. Default = \'\'.

    Returns
    -------
    have_scaled : bool, optional
        Has a graph been scaled?
    """
