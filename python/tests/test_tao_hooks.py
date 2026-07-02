"""Tests for the hand-written Tao hooks (``pybmad.tao.hooks``).

The registration/readback tests exercise the property wiring directly. The
firing test drives a real Tao session (via ``pytao``) to confirm the
Fortran -> C -> Python hook path works end to end; it is skipped when pytao or
the Tao example lattice is unavailable.
"""

from __future__ import annotations

import logging

import pytest
from pybmad import tao

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
