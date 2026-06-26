"""bsim routines"""

from collections.abc import Sequence

import _pybmad


def bbu_add_a_bunch(lat: _pybmad.LatStruct, bbu_beam: _pybmad.BbuBeamStruct, bbu_param: _pybmad.BbuParamStruct, beam_init: _pybmad.BeamInitStruct) -> None:
    """
    Wrapper for Fortran routine bbu_add_a_bunch

    Parameters
    ----------
    lat : LatStruct

    bbu_beam : BbuBeamStruct

    bbu_param : BbuParamStruct

    beam_init : BeamInitStruct
    """

def bbu_hom_voltage_calc(lat: _pybmad.LatStruct, bbu_beam: _pybmad.BbuBeamStruct, n_period: int, ix_stage_last_tracked: int) -> None:
    """
    Wrapper for Fortran routine bbu_hom_voltage_calc

    Parameters
    ----------
    lat : LatStruct

    bbu_beam : BbuBeamStruct

    n_period : int

    ix_stage_last_tracked : int
    """

def bbu_remove_head_bunch(bbu_beam: _pybmad.BbuBeamStruct) -> None:
    """
    Wrapper for Fortran routine bbu_remove_head_bunch

    Parameters
    ----------
    bbu_beam : BbuBeamStruct
    """

def bbu_setup(lat: _pybmad.LatStruct, dt_bunch: float, bbu_param: _pybmad.BbuParamStruct, bbu_beam: _pybmad.BbuBeamStruct) -> None:
    """
    Wrapper for Fortran routine bbu_setup

    Parameters
    ----------
    lat : LatStruct

    dt_bunch : float

    bbu_param : BbuParamStruct

    bbu_beam : BbuBeamStruct
    """

def bbu_track_a_stage(lat: _pybmad.LatStruct, bbu_beam: _pybmad.BbuBeamStruct, bbu_param: _pybmad.BbuParamStruct, lost: bool, ix_stage_tracked: int) -> None:
    """
    Wrapper for Fortran routine bbu_track_a_stage

    Parameters
    ----------
    lat : LatStruct

    bbu_beam : BbuBeamStruct

    bbu_param : BbuParamStruct

    lost : bool

    ix_stage_tracked : int
    """

class BbuTrackAll:
    """bbu_track_all return type"""

    @property
    def hom_voltage_normalized(self) -> float: ...

    @property
    def growth_rate(self) -> float: ...

    @property
    def lost(self) -> bool: ...

    @property
    def irep(self) -> int: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def bbu_track_all(lat: _pybmad.LatStruct, bbu_beam: _pybmad.BbuBeamStruct, bbu_param: _pybmad.BbuParamStruct, beam_init: _pybmad.BeamInitStruct) -> BbuTrackAll:
    """
    Wrapper for Fortran routine bbu_track_all

    Parameters
    ----------
    lat : LatStruct

    bbu_beam : BbuBeamStruct

    bbu_param : BbuParamStruct

    beam_init : BeamInitStruct

    Returns
    -------
    hom_voltage_normalized : float
        HOM voltage normalized

    growth_rate : float
        Growth rate

    lost : bool
        Lost

    irep : int
        irep
    """

def check_rf_freq(lat: _pybmad.LatStruct, fb: float) -> None:
    """
    Wrapper for Fortran routine check_rf_freq

    Parameters
    ----------
    lat : LatStruct

    fb : float
    """

def count_lines_in_file(file_name: str) -> int:
    """
    Wrapper for Fortran routine count_lines_in_file

    Parameters
    ----------
    file_name : str

    Returns
    -------
    lines : int
    """

def hom_voltage(lr_wake: _pybmad.WakeLrModeStruct) -> float:
    """
    Wrapper for Fortran routine hom_voltage

    Parameters
    ----------
    lr_wake : WakeLrModeStruct

    Returns
    -------
    voltage : float
    """

def insert_phase_trombone(branch: _pybmad.BranchStruct) -> None:
    """
    Wrapper for Fortran routine insert_phase_trombone

    Parameters
    ----------
    branch : BranchStruct
        Lattice branch.
        This parameter is an input/output and is modified in-place.
        As an output, branch: Lattice branch with trumbone at branch.ele(1).
    """

def logical_to_python(logic: bool) -> str:
    """
    Wrapper for Fortran routine logical_to_python

    Parameters
    ----------
    logic : bool

    Returns
    -------
    string : str
    """

def rf_cav_names(lat: _pybmad.LatStruct) -> None:
    """
    Wrapper for Fortran routine rf_cav_names

    Parameters
    ----------
    lat : LatStruct
    """

def set_tune_3d(branch: _pybmad.BranchStruct, target_tunes: Sequence[float], mask: str | None = None, use_phase_trombone: bool | None = None, z_tune_set: bool | None = None, group_knobs: Sequence[str] | None = None, print_err: bool | None = None) -> bool:
    """
    Wrapper for Fortran routine set_tune_3d

    Parameters
    ----------
    branch : BranchStruct
        This parameter is an input/output and is modified in-place.
        As an output, branch: with adjusted quads and RF to match desired tunes.

    target_tunes : 1D array of float (shape: 3)
        tunes for a, b, z modes (rad/2pi). Must include integer part.

    mask : str, optional

    use_phase_trombone : bool, optional
        Default False. If true, use a match element in phase trombone mode to adjust the tunes. The match element
        must be the first element in the lattice. Use insert_phase_trombone to insert one.

    z_tune_set : bool, optional
        Default True. If false, do not try to set the synch tune.

    group_knobs : 1D array of str (shape: 2), optional
        If set non-blank, use these group elements for tuning.

    print_err : bool, optional
        Print error message if there is a problem? Default is True.

    Returns
    -------
    everything_ok : bool
        Returns true or false if set was successful.
    """

def write_bunch_by_bunch_info(lat: _pybmad.LatStruct, bbu_beam: _pybmad.BbuBeamStruct, bbu_param: _pybmad.BbuParamStruct, this_stage: _pybmad.BbuStageStruct) -> None:
    """
    Wrapper for Fortran routine write_bunch_by_bunch_info

    Parameters
    ----------
    lat : LatStruct

    bbu_beam : BbuBeamStruct

    bbu_param : BbuParamStruct

    this_stage : BbuStageStruct
    """
