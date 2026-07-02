"""Tests for the hand-written Bmad tracking hooks (``pybmad.bmad.hooks``).

Only the hooks that can be triggered from a plain single-particle / single-bunch
track through a simple lattice are exercised here. ``wall_hit_handler_custom``
(needs a Runge-Kutta wall intersection), ``time_runge_kutta_periodic_kick`` and
``track1_wake`` (need specific tracking methods / wakefields) are not covered.
"""

from __future__ import annotations

import logging
import os
import pathlib
import subprocess
import sys

import pybmad
import pytest
from conftest import TESTS_ROOT
from pybmad import bmad

logger = logging.getLogger(__name__)

LAT_FILE = TESTS_ROOT / "data" / "wall_generator" / "lat.bmad"

ALL_HOOKS = (
    "time_runge_kutta_periodic_kick",
    "track1_bunch",
    "track1_custom",
    "track_many",
    "track1_postprocess",
    "track1_preprocess",
    "track1_spin_custom",
    "track1_wake",
    "wall_hit_handler_custom",
)


@pytest.fixture(autouse=True)
def clear_hooks():
    """Restore global state (hooks + ``bmad_com``) after each test.

    Hooks and ``bmad_com`` are process-global singletons, so a test that installs
    a hook or flips e.g. ``spin_tracking_on`` would otherwise leak into every
    later test.
    """
    spin_tracking_on = pybmad.get_bmad_com().spin_tracking_on
    yield
    for name in ALL_HOOKS:
        setattr(bmad.hooks, name, None)
    pybmad.get_bmad_com().spin_tracking_on = spin_tracking_on


@pytest.fixture
def lat() -> pybmad.LatStruct:
    res = pybmad.bmad_parser(str(LAT_FILE))
    assert not res.err_flag
    return res.lat


def make_orbit(lat: pybmad.LatStruct, x: float = 1e-4) -> pybmad.CoordStructAlloc1D:
    """Build a fresh, initialized orbit array for ``lat``'s root branch."""
    branch = lat.branch[0]
    orbit = pybmad.CoordStruct.new_array1d(0)
    pybmad.reallocate_coord(orbit, lat, 0)
    pybmad.init_coord(orbit[0], [x, 0, 0, 0, 0, 0], branch.ele[0], pybmad.DOWNSTREAM_END)
    pybmad.lattice_bookkeeper(lat)
    return orbit


def make_bunch(lat: pybmad.LatStruct, x: float = 1e-4) -> pybmad.BunchStruct:
    branch = lat.branch[0]
    bunch = pybmad.BunchStruct()
    bunch.particle.resize(1)
    pybmad.init_coord(bunch.particle[0], [x, 0, 0, 0, 0, 0], branch.ele[0], pybmad.DOWNSTREAM_END)
    bunch.n_live = 1
    return bunch


def test_track1_postprocess_fires_with_live_proxies(lat):
    branch = lat.branch[0]
    orbit = make_orbit(lat)
    seen = []

    def hook(start_orb, ele, param, end_orb):
        logger.info(f"\nHook: {start_orb=} {ele=} {param=} {end_orb=}")
        seen.append((ele.name, float(end_orb.vec[0])))

    bmad.hooks.track1_postprocess = hook
    result = pybmad.track1(orbit[0], branch.ele[1], branch.param)

    assert len(seen) == 1
    assert seen[0][0] == branch.ele[1].name
    assert seen[0][1] == pytest.approx(float(result.end_orb.vec[0]))


def test_assignment_readback_and_none_clears(lat):
    """The property reads back the current callable and clears when set to None."""
    branch = lat.branch[0]
    orbit = make_orbit(lat)
    calls = []

    def postprocess_hook(*_):
        calls.append(1)

    assert bmad.hooks.track1_postprocess is None
    bmad.hooks.track1_postprocess = postprocess_hook
    assert bmad.hooks.track1_postprocess is postprocess_hook

    pybmad.track1(orbit[0], branch.ele[1], branch.param)
    assert len(calls) == 1

    bmad.hooks.track1_postprocess = None
    assert bmad.hooks.track1_postprocess is None
    pybmad.track1(orbit[0], branch.ele[1], branch.param)
    assert len(calls) == 1


@pytest.mark.filterwarnings("ignore::pytest.PytestUnraisableExceptionWarning")
def test_exception_in_hook_does_not_crash(lat):
    """A raising callback is reported, not propagated into Fortran."""
    branch = lat.branch[0]
    orbit = make_orbit(lat)
    calls = []

    def hook(start_orb, ele, param, end_orb):
        logger.info("\nHook: %s %s %s %s", start_orb, ele, param, end_orb)
        calls.append(1)
        raise ValueError("intentional")

    bmad.hooks.track1_postprocess = hook
    result = pybmad.track1(orbit[0], branch.ele[1], branch.param)  # must not crash
    assert len(calls) == 1
    assert result.end_orb is not None


def test_track1_preprocess_fires_and_optional_track_is_none(lat):
    branch = lat.branch[0]
    orbit = make_orbit(lat)
    seen = []

    def hook(start_orb, ele, param, err_flag, finished, radiation_included, track):
        logger.info(
            f"\nHook: {start_orb=} {ele=} {param=} {err_flag=} {finished=} {radiation_included=} {track=}"
        )
        seen.append((ele.name, track))

    bmad.hooks.track1_preprocess = hook
    pybmad.track1(orbit[0], branch.ele[1], branch.param)

    assert len(seen) == 1
    assert seen[0][1] is None  # `track` argument was absent


def test_track1_custom_replaces_tracking_and_writes_back(lat):
    """A custom-tracking element routes track1 through the hook; the orbit proxy
    is inout, so mutating it propagates back into the track1 result."""
    branch = lat.branch[0]
    ele = branch.ele[1]
    ele.tracking_method = pybmad.CUSTOM
    orbit = make_orbit(lat)

    def hook(orb, e, param, err_flag, finished, track):
        logger.info(f"\nHook: {orb=} {e=} {param=} {err_flag=} {finished=} {track=}")
        orb.vec[0] = 0.00123  # this element "tracked" the particle to here
        return (False, True)  # err_flag=False, finished=True (hook did the tracking)

    bmad.hooks.track1_custom = hook
    result = pybmad.track1(orbit[0], ele, branch.param)

    assert float(result.end_orb.vec[0]) == pytest.approx(0.00123)


def test_track1_spin_custom_fires(lat):
    branch = lat.branch[0]
    ele = branch.ele[1]
    pybmad.get_bmad_com().spin_tracking_on = True
    ele.spin_tracking_method = pybmad.CUSTOM
    orbit = make_orbit(lat)
    seen = []

    def hook(start_orb, e, param, end_orb, err_flag, make_quaternion):
        logger.info(f"\nHook: {start_orb=} {e=} {param=} {end_orb=} {err_flag=} {make_quaternion=}")
        seen.append((e.name, make_quaternion))
        return False  # err_flag

    bmad.hooks.track1_spin_custom = hook
    pybmad.track1(orbit[0], ele, branch.param)

    assert len(seen) == 1
    assert seen[0][0] == ele.name


def test_track_many_receives_live_coord_array_view(lat):
    branch = lat.branch[0]
    orbit = make_orbit(lat)
    captured = {}

    def hook(finished, lat_, orbit_, ix_start, ix_end, direction, ix_branch, track_state):
        logger.info(
            f"\nHook: {finished=} {lat_=} {orbit_=} {ix_start=} {ix_end=} {direction=} {ix_branch=} {track_state}"
        )
        captured["n"] = len(orbit_)
        captured["start_x"] = float(orbit_[0].vec[0])
        captured["use_name"] = lat_.use_name
        captured["ix_branch"] = ix_branch
        return False  # not finished: let Bmad do the actual tracking

    bmad.hooks.track_many = hook
    pybmad.track_many(lat, orbit.view(), 0, branch.n_ele_track, 1)

    assert captured["n"] == branch.n_ele_track + 1
    assert captured["start_x"] == pytest.approx(1e-4)
    assert captured["use_name"] == lat.use_name
    assert captured["ix_branch"] is None  # optional argument was absent


def test_track1_bunch_fires(lat):
    branch = lat.branch[0]
    bunch = make_bunch(lat)
    seen = []

    def hook(b, ele, err, centroid, direction, finished, bunch_track):
        logger.info(f"\nHook: {b=} {ele=} {err=} {centroid=} {direction=} {finished=} {bunch_track}")
        seen.append((ele.name, centroid, direction, bunch_track, len(b.particle)))
        return (False, True)  # err=False, finished=True

    bmad.hooks.track1_bunch = hook
    pybmad.track1_bunch(bunch, branch.ele[1])

    assert len(seen) == 1
    name, centroid, _direction, bunch_track, n_particle = seen[0]
    assert name == branch.ele[1].name
    assert centroid is None
    assert bunch_track is None  # optional arguments absent
    assert n_particle == 1


def test_installed_hook_survives_interpreter_shutdown():
    """Regression: a hook left installed at exit must not crash the process during
    interpreter finalization.

    The Python callable is captured by an owning ``nb::object`` in a C++ file-static
    ``std::function``; if not released before ``Py_Finalize`` it segfaults at exit.
    Installing a hook is enough to arm this -- no tracking is needed.
    """
    code = (
        "import pybmad\n"
        "from pybmad import bmad, tao\n"
        "bmad.hooks.track1_postprocess = lambda *a: None\n"
        "tao.hooks.lattice_calc = lambda ok: None\n"
    )
    pybmad_parent = str(pathlib.Path(pybmad.__file__).resolve().parents[1])
    proc = subprocess.run(
        [sys.executable, "-c", code],
        env={**os.environ, "PYTHONPATH": pybmad_parent},
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr.decode()
