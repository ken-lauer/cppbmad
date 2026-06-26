from __future__ import annotations

from pybmad import AllEncompassingStruct, TestSubStruct


def test_same_underlying_struct_is_equal():
    aes = AllEncompassingStruct()
    a = aes.type_0d
    b = aes.type_0d
    assert a == b
    assert hash(a) == hash(b)


def test_distinct_structs_are_not_equal():
    s1 = TestSubStruct()
    s2 = TestSubStruct()
    assert s1 != s2
    assert hash(s1) != hash(s2)


def test_same_underlying_struct_collapses_in_set():
    aes = AllEncompassingStruct()
    assert len({aes.type_0d, aes.type_0d}) == 1


def test_struct_usable_as_dict_key():
    s1 = TestSubStruct()
    s2 = TestSubStruct()
    d = {s1: "one", s2: "two"}
    assert len(d) == 2
    assert d[s1] == "one"


def test_eq_with_non_struct_is_false():
    s1 = TestSubStruct()
    assert (s1 == 12345) is False
    assert s1 != 12345
