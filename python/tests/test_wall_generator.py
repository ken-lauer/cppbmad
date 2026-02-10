from __future__ import annotations

import sys

import numpy as np
import pytest
from conftest import TESTS_ROOT
from examples.wall_generator import get_wall, get_wall_contour

import pybmad
from pybmad import BranchStruct

WALL_ROOT = TESTS_ROOT / "data" / "wall_generator"


@pytest.fixture
def lattice() -> pybmad.LatStruct:
    res = pybmad.bmad_parser(str(WALL_ROOT / "lat.bmad"))
    assert not res.err_flag
    return res.lat


@pytest.fixture
def branch(lattice: pybmad.LatStruct) -> pybmad.BranchStruct:
    branch = lattice.branch[0]
    assert branch.wall3d
    return branch


def test_get_wall(branch: BranchStruct):
    # generated via: $ACC_ROOT_DIR/debug/bin/wall_generator lat.bmad 8 0.1
    expected = np.loadtxt(WALL_ROOT / "wall.out")

    data = []
    for section, i_angle, ele, orb_out in get_wall(branch, n_angles=8):
        data.append([*orb_out.vec, ele.ix_ele, i_angle, section.s])

    data = np.asarray(data)
    np.testing.assert_allclose(data, expected)


def test_get_wall_contour(branch: BranchStruct):
    # generated via: $ACC_ROOT_DIR/debug/bin/wall_generator lat.bmad 8 0.1
    expected = np.loadtxt(WALL_ROOT / "wall_contour.out")

    data = []
    for point_s, i_angle, ele, orb_out in get_wall_contour(branch, n_angles=8, ds=0.1):
        data.append([*orb_out.vec, ele.ix_ele, i_angle, point_s])

    data = np.asarray(data)
    np.testing.assert_allclose(data, expected)


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", *sys.argv[1:], __file__]))
