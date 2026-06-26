from __future__ import annotations

from pybmad import AllEncompassingStruct, TestSubStruct


def test_string_member_default_is_empty():
    assert TestSubStruct().sr.file == ""


def test_string_member_roundtrip():
    s = TestSubStruct()
    s.sr.file = "lattice.bmad"
    assert s.sr.file == "lattice.bmad"


def test_string_member_trailing_whitespace_trimmed():
    # Fortran character(N) is space-padded; the getter trims trailing spaces.
    s = TestSubStruct()
    s.sr.file = "abc"
    assert s.sr.file == "abc"
    assert len(s.sr.file) == 3


def test_string_member_via_nested_struct():
    aes = AllEncompassingStruct()
    aes.type_0d.sr.file = "nested.bmad"
    assert aes.type_0d.sr.file == "nested.bmad"


def test_string_member_truncated_to_capacity():
    s = TestSubStruct()
    s.sr.file = "x" * 250
    assert len(s.sr.file) == 200  # character(200)
