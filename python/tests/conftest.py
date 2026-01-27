from __future__ import annotations

import pytest

import pybmad


@pytest.fixture(autouse=True, scope="session")
def do_not_print():
    pybmad.out_io_print_and_capture_setup(
        print_on=False,
        capture_state="BUFFERED",
    )
