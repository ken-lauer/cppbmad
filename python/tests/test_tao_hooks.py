"""Tests for the hand-written Tao hooks (``pybmad.tao.hooks``).

The registration/readback tests exercise the property wiring directly. The
firing test drives a real Tao session (via ``pytao``) to confirm the
Fortran -> C -> Python hook path works end to end; it is skipped when pytao or
the Tao example lattice is unavailable.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pybmad
import pytest
from pybmad import tao

if TYPE_CHECKING:
    from pytao import Tao

logger = logging.getLogger(__name__)

TAO_HOOKS = ("lattice_calc", "optimizer", "merit_var", "merit_data")

TAO_INIT_FILE = "$ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init"


@pytest.fixture(autouse=True)
def clear_tao_hooks():
    """Tao hooks are process-global singletons; clear them after each test."""
    yield
    for name in TAO_HOOKS:
        setattr(tao.hooks, name, None)


@pytest.mark.parametrize("name", TAO_HOOKS)
def test_hook_registration_readback_and_clear(name):
    """Each property starts None, reads back the assigned callable, and clears."""

    def hook(*_args):
        return None

    assert getattr(tao.hooks, name) is None
    setattr(tao.hooks, name, hook)
    assert getattr(tao.hooks, name) is hook
    setattr(tao.hooks, name, None)
    assert getattr(tao.hooks, name) is None


@pytest.fixture
def tao_session():
    Tao = pytest.importorskip("pytao").Tao
    try:
        return Tao(init_file=TAO_INIT_FILE, noplot=True)
    except Exception as exc:  # ACC_ROOT_DIR unset, example lattice missing, ...
        pytest.skip(f"Could not initialize Tao: {exc}")


def test_lattice_calc_hook_fires_with_current_status(tao_session):
    """Installing ``lattice_calc`` routes Tao's lattice calculation through the
    hook, which receives the current ``calc_ok`` status."""
    seen = []

    def hook(calc_ok):
        seen.append(calc_ok)
        # implicit return None - leave calc_ok unchanged

    tao.hooks.lattice_calc = hook
    # Changing a variable forces a lattice recalculation, firing the hook.
    tao_session.cmd("set ele Q1 k1 = 0.30001")

    assert seen, "lattice_calc hook did not fire"
    assert all(v is True for v in seen)


def _custom_merit_datum():
    """Route the first datum of universe 1 to the ``merit_data`` hook.

    Tao only calls ``tao_hook_merit_data`` for data whose ``merit_type`` is not
    one it recognizes, and its ``set`` command rejects unknown merit types -- so
    set it directly on the datum struct instead. Returns the live datum proxy.
    """
    datum = pybmad.get_super_universe().u[0].data[1]
    datum.merit_type = "custom_hook"
    return datum


def _run_optimizer(tao_session: Tao):
    """Force merit evaluations (which fire the merit hooks) via a short optimizer
    run. A custom ``merit_type`` makes some commands raise ``TaoCommandError``
    (e.g. once a datum is marked non-existent); that is expected here and ignored
    per command -- the observable is the datum/var state, not the command status."""
    tao_session.cmds(
        [
            "use var *",
            "set global optimizer = lmdif",
            "set global n_opti_loops = 1",
            "run",
        ],
        raises=False,
    )


def test_merit_var_hook_fires_for_custom_merit_type(tao_session):
    """A variable with a custom ``merit_type`` routes its merit contribution
    through the ``merit_var`` hook, which receives a live ``TaoVarStruct``."""
    var = pybmad.get_super_universe().var[1]
    var.merit_type = "custom_hook"
    seen = []

    def hook(i_uni, j_var, var):
        logger.info(f"\nHook: {i_uni=} {j_var=} {var=}")
        seen.append((i_uni, j_var, var.merit_type))
        var.delta_merit = 0.0  # the hook supplies the contribution

    tao.hooks.merit_var = hook
    _run_optimizer(tao_session)

    assert seen, "merit_var hook did not fire"
    assert seen[0][2] == "custom_hook"


def test_merit_data_none_return_preserves_valid_value_set(tao_session):
    """Regression for the None-return contract: returning None from ``merit_data``
    must leave ``valid_value_set`` unchanged (it arrives True), not coerce it to
    False -- which Tao treats as "merit_type not recognized" and marks the datum
    non-existent."""
    datum = _custom_merit_datum()
    received = []

    def hook(i_uni, j_data, datum, valid_value_set):
        logger.info(f"\nHook: {i_uni=} {j_data=} {datum=} {valid_value_set=}")
        received.append(valid_value_set)
        # implicit return None -> leave valid_value_set unchanged

    tao.hooks.merit_data = hook
    _run_optimizer(tao_session)

    assert received
    assert received[0] is True
    assert datum.exists is True


def test_merit_data_false_return_is_honored(tao_session):
    """Positive control for the above: an explicit False return is written back
    (marking the datum non-existent), confirming None is distinct from False."""
    datum = _custom_merit_datum()

    def hook(i_uni, j_data, datum, valid_value_set):
        logger.info(f"\nHook: {i_uni=} {j_data=} {datum=} {valid_value_set=}")
        return False

    tao.hooks.merit_data = hook
    _run_optimizer(tao_session)

    assert datum.exists is False
