#include <bmad.hpp>

#include "doctest.h"

using namespace Bmad;

TEST_CASE("Enum constants - element types are distinct") {
  // Core element types should all be distinct positive integers
  CHECK(DRIFT > 0);
  CHECK(SBEND > 0);
  CHECK(QUADRUPOLE > 0);
  CHECK(RFCAVITY > 0);
  CHECK(MARKER > 0);
  CHECK(SOLENOID > 0);

  CHECK(DRIFT != SBEND);
  CHECK(DRIFT != QUADRUPOLE);
  CHECK(SBEND != QUADRUPOLE);
  CHECK(RFCAVITY != DRIFT);
  CHECK(MARKER != DRIFT);
}

TEST_CASE("Enum constants - tracking methods are distinct") {
  CHECK(BMAD_STANDARD > 0);
  CHECK(RUNGE_KUTTA > 0);
  CHECK(LINEAR > 0);
  CHECK(TIME_RUNGE_KUTTA > 0);

  CHECK(BMAD_STANDARD != RUNGE_KUTTA);
  CHECK(BMAD_STANDARD != LINEAR);
  CHECK(RUNGE_KUTTA != LINEAR);
  CHECK(RUNGE_KUTTA != TIME_RUNGE_KUTTA);
}

TEST_CASE("Enum constants - particle states") {
  // PRE_BORN is 0 by convention (conforms to OpenPMD standard)
  CHECK(PRE_BORN == 0);
  // ALIVE and LOST must be distinct and positive
  CHECK(ALIVE > 0);
  CHECK(LOST > 0);
  CHECK(ALIVE != LOST);

  // Lost states in different planes should all be distinct
  CHECK(LOST_NEG_X != LOST_POS_X);
  CHECK(LOST_NEG_Y != LOST_POS_Y);
  CHECK(LOST_NEG_X != LOST_NEG_Y);

  // Old names should match new names
  CHECK(LOST_NEG_X == LOST_NEG_X_APERTURE);
  CHECK(LOST_POS_X == LOST_POS_X_APERTURE);
  CHECK(LOST_NEG_Y == LOST_NEG_Y_APERTURE);
  CHECK(LOST_POS_Y == LOST_POS_Y_APERTURE);
  CHECK(LOST_Z == LOST_Z_APERTURE);
  CHECK(LOST_PZ == LOST_PZ_APERTURE);
}

TEST_CASE("Enum constants - on/off ordering") {
  CHECK(OFF != ON);
  CHECK(OFF > 0);
  CHECK(ON > 0);
}

TEST_CASE("Enum constants - offset relationships") {
  // VAR_OFFSET must be greater than OLD_CONTROL_VAR_OFFSET (documented invariant)
  CHECK(VAR_OFFSET > OLD_CONTROL_VAR_OFFSET);
  // TAYLOR_OFFSET should be much larger than other offsets
  CHECK(TAYLOR_OFFSET > VAR_OFFSET);
  CHECK(TAYLOR_OFFSET > OLD_CONTROL_VAR_OFFSET);
  CHECK(N_VAR_MAX > 0);
}

TEST_CASE("Enum constants - aperture types are distinct") {
  CHECK(AUTO_APERTURE != RECTANGULAR);
  CHECK(RECTANGULAR != ELLIPTICAL);
  CHECK(ELLIPTICAL != WALL3D);
  CHECK(AUTO_APERTURE > 0);
}

TEST_CASE("Enum constants - geometry") {
  CHECK(ENTRANCE_END != EXIT_END);
  CHECK(X_PLANE != Y_PLANE);
  CHECK(Y_PLANE != Z_PLANE);
}

TEST_CASE("EleKey enum class") {
  // EleKey values should match the corresponding const int values
  CHECK(static_cast<size_t>(EleKey::DRIFT) == 1);
}

TEST_CASE("EleAttribute enum class") {
  // EleAttribute uses 0-based indexing (Fortran value - 1)
  // L is always the first attribute (index 0)
  CHECK(static_cast<size_t>(EleAttribute::L) == Bmad::L - 1);

  // TILT and ROLL are documented aliases for the same attribute
  CHECK(EleAttribute::TILT == EleAttribute::ROLL);

  // K1 and K2 should be distinct
  CHECK(EleAttribute::K1 != EleAttribute::K2);
}
