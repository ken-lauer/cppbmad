from __future__ import annotations

import pathlib

import pytest

import pybmad

TESTS_ROOT = pathlib.Path(__file__).resolve().parent
CPPBMAD_REPO_ROOT = TESTS_ROOT.parents[1]


@pytest.fixture(autouse=True, scope="session")
def do_not_print():
    pybmad.out_io_print_and_capture_setup(
        print_on=False,
        capture_state="BUFFERED",
    )
