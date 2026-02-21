from __future__ import annotations

import gc
import sys

import pytest
from conftest import TESTS_ROOT

import pybmad

WALL_ROOT = TESTS_ROOT / "data" / "wall_generator"


def _parse_lattice():
    res = pybmad.bmad_parser(str(WALL_ROOT / "lat.bmad"))
    assert not res.err_flag
    return res


class TestPropertyKeepAlive:
    def test_branch_property_keeps_lattice_alive(self):
        res = _parse_lattice()
        lattice = res.lat
        rc_before = sys.getrefcount(lattice)

        branch_alloc = lattice.branch

        rc_after = sys.getrefcount(lattice)
        assert rc_after > rc_before, (
            f"lattice refcount did not increase after .branch access "
            f"({rc_before} -> {rc_after}); keep_alive on def_property_readonly is broken"
        )
        del branch_alloc

    def test_getitem_keeps_alloc_alive(self):
        res = _parse_lattice()
        lattice = res.lat
        branch_alloc = lattice.branch
        rc_before = sys.getrefcount(branch_alloc)

        branch = branch_alloc[0]

        rc_after = sys.getrefcount(branch_alloc)
        assert rc_after > rc_before, (
            f"branch_alloc refcount did not increase after [0] access "
            f"({rc_before} -> {rc_after}); keep_alive on __getitem__ is broken"
        )
        del branch

    def test_lat_property_keeps_res_alive(self):
        res = _parse_lattice()
        rc_before = sys.getrefcount(res)

        lattice = res.lat

        rc_after = sys.getrefcount(res)
        assert rc_after > rc_before, (
            f"res refcount did not increase after .lat access "
            f"({rc_before} -> {rc_after}); reference_internal on def_readonly is broken"
        )
        del lattice


class TestLifetimeChainAfterDeletion:
    @pytest.fixture(autouse=True, params=list(range(5)))
    def _iter(self, request: pytest.FixtureRequest) -> int:
        return request.param

    def test_branch_survives_fixture_pattern(self):
        def fixture_branch():
            res = _parse_lattice()
            lattice = res.lat
            branch = lattice.branch[0]
            assert branch.wall3d
            return branch

        branch = fixture_branch()
        gc.collect()

        # branch.wall3d must still be accessible
        assert branch.wall3d, "wall3d is falsy after intermediates were GC'd"
        wall3d = branch.wall3d[0]
        assert len(wall3d.section) > 0, "deep wall3d.section access failed"

    def test_branch_survives_explicit_deletion(self):
        res = _parse_lattice()
        lattice = res.lat
        branch_alloc = lattice.branch
        branch = branch_alloc[0]
        assert branch.wall3d

        del branch_alloc
        gc.collect()
        assert branch.wall3d, "wall3d invalid after del branch_alloc"

        del lattice
        gc.collect()
        assert branch.wall3d, "wall3d invalid after del lattice"

        del res
        gc.collect()
        assert branch.wall3d, "wall3d invalid after del res"

        # Deep access should still work
        wall3d = branch.wall3d[0]
        assert len(wall3d.section) > 0

    def test_oneliner_branch_access(self):
        res = _parse_lattice()
        lattice = res.lat
        branch = lattice.branch[0]  # alloc is a temporary

        del res
        del lattice
        gc.collect()

        assert branch.wall3d, "wall3d invalid after oneliner pattern"


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", *sys.argv[1:], __file__]))
