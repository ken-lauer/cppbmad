from __future__ import annotations

import enum

import pybmad
from pybmad import (
    ALIVE,
    AUTO_APERTURE,
    BMAD_STANDARD,
    DRIFT,
    ELLIPTICAL,
    ENTRANCE_END,
    EXIT_END,
    FREE,
    GROUP_LORD,
    KICKER,
    LINEAR,
    LOST,
    LOST_NEG_X,
    LOST_NEG_Y,
    LOST_POS_X,
    LOST_POS_Y,
    LOST_PZ,
    LOST_Z,
    MARKER,
    OFF,
    OLD_CONTROL_VAR_OFFSET,
    ON,
    PRE_BORN,
    QUADRUPOLE,
    RECTANGULAR,
    RFCAVITY,
    RUNGE_KUTTA,
    SBEND,
    SEXTUPOLE,
    SOLENOID,
    TAYLOR_OFFSET,
    VAR_OFFSET,
    WALL3D,
    X_PLANE,
    Y_PLANE,
    Z_PLANE,
    EleAttribute,
    EleKey,
)


class TestEnumConstants:
    """Test that enum constants are accessible and have correct relationships."""

    def test_element_types_are_distinct(self):
        elements = {DRIFT, SBEND, QUADRUPOLE, SEXTUPOLE, RFCAVITY, MARKER, SOLENOID, KICKER}
        # All should be positive integers
        for elem in elements:
            assert isinstance(elem, int)
            assert elem > 0
        # All should be distinct
        assert len(elements) == 8

    def test_tracking_methods_are_distinct(self):
        methods = {BMAD_STANDARD, RUNGE_KUTTA, LINEAR}
        for m in methods:
            assert isinstance(m, int)
            assert m > 0
        assert len(methods) == 3

    def test_particle_states(self):
        assert PRE_BORN == 0  # OpenPMD convention
        assert ALIVE > 0
        assert LOST > 0
        assert ALIVE != LOST

        # Lost in different planes should be distinct
        lost_states = {LOST_NEG_X, LOST_POS_X, LOST_NEG_Y, LOST_POS_Y, LOST_Z, LOST_PZ}
        assert len(lost_states) == 6

    def test_on_off(self):
        assert OFF != ON
        assert isinstance(OFF, int)
        assert isinstance(ON, int)

    def test_offset_relationships(self):
        assert VAR_OFFSET > OLD_CONTROL_VAR_OFFSET
        assert TAYLOR_OFFSET > VAR_OFFSET

    def test_aperture_types_are_distinct(self):
        apertures = {AUTO_APERTURE, RECTANGULAR, ELLIPTICAL, WALL3D}
        assert len(apertures) == 4

    def test_geometry_constants(self):
        assert ENTRANCE_END != EXIT_END
        planes = {X_PLANE, Y_PLANE, Z_PLANE}
        assert len(planes) == 3

    def test_control_types_are_distinct(self):
        assert FREE != GROUP_LORD


class TestEleKeyEnum:
    """Test the EleKey IntEnum class."""

    def test_is_int_enum(self):
        assert issubclass(EleKey, enum.IntEnum)

    def test_element_lookup_by_name(self):
        assert EleKey["DRIFT"] == EleKey.DRIFT
        assert EleKey["QUADRUPOLE"] == EleKey.QUADRUPOLE

    def test_matches_module_constants(self):
        """EleKey enum values should match the top-level constants."""
        assert EleKey.DRIFT == DRIFT
        assert EleKey.SBEND == SBEND
        assert EleKey.QUADRUPOLE == QUADRUPOLE
        assert EleKey.SEXTUPOLE == SEXTUPOLE
        assert EleKey.RFCAVITY == RFCAVITY
        assert EleKey.SOLENOID == SOLENOID
        assert EleKey.MARKER == MARKER

    def test_usable_as_int(self):
        """IntEnum should be usable wherever an int is expected."""
        assert EleKey.DRIFT + 0 == DRIFT
        assert int(EleKey.DRIFT) == DRIFT

    def test_members_are_positive(self):
        for member in EleKey:
            assert member.value > 0, f"{member.name} has non-positive value {member.value}"


class TestEleAttributeEnum:
    """Test the EleAttribute IntEnum class."""

    def test_is_int_enum(self):
        assert issubclass(EleAttribute, enum.IntEnum)

    def test_l_is_first(self):
        """L (length) is always index 0 (Fortran 1 - 1)."""
        assert EleAttribute.L == 0

    def test_tilt_roll_alias(self):
        """TILT and ROLL are documented aliases for the same attribute."""
        assert EleAttribute.TILT == EleAttribute.ROLL

    def test_k1_k2_distinct(self):
        assert EleAttribute.K1 != EleAttribute.K2


class TestEnumAccessFromModule:
    """Test that enums are accessible via the pybmad module directly."""

    def test_constants_on_module(self):
        assert hasattr(pybmad, "DRIFT")
        assert hasattr(pybmad, "ALIVE")
        assert hasattr(pybmad, "BMAD_STANDARD")

    def test_enum_classes_on_module(self):
        assert hasattr(pybmad, "EleKey")
        assert hasattr(pybmad, "EleAttribute")
