from __future__ import annotations

import numpy as np

import pybmad


def get_info_example(_tao):
    """
    This is meant to be used by tao-subproc.py.

    SubprocessTao may only call importable functions (i.e., those outside of
    __main__).
    """
    s = pybmad.get_super_universe()
    univ = s.u[0]
    lattice = univ.model.lat.use_name
    ele1_name = univ.model.lat.ele[1].name
    orbit = univ.design.tao_branch[0].orbit[0].vec

    # Note that we are also somewhat limited on supported return types here.
    # Native types: float, int, str, bytes, booll
    # Native containers: list, tuple, set, dict
    # And np.ndarray, which is recommended for passing back large arrays.
    return {
        "lattice": lattice,
        "ele1_name": ele1_name,
        "orbit": np.asarray(orbit),
    }
