#include <bmad.hpp>

#include "doctest.h"

using namespace Bmad;

// ---------------------------------------------------------------------------
// Exception type hierarchy
// ---------------------------------------------------------------------------

TEST_CASE("BmadException inherits from std::runtime_error") {
  BmadException ex("test message");
  // Should be catchable as runtime_error
  CHECK_THROWS_AS(throw ex, std::runtime_error);
  CHECK(std::string(ex.what()) == "test message");
}

TEST_CASE("InvalidIndexException") {
  SUBCASE("message format") {
    InvalidIndexException ex("array", 5, 10);
    std::string msg = ex.what();
    CHECK(msg.find("Invalid") != std::string::npos);
    CHECK(msg.find("array") != std::string::npos);
    CHECK(msg.find("5") != std::string::npos);
    // max_value=10 means valid range 0-9
    CHECK(msg.find("0-9") != std::string::npos);
  }

  SUBCASE("is a BmadException") {
    InvalidIndexException ex("test", 0, 0);
    CHECK_THROWS_AS(throw ex, BmadException);
  }

  SUBCASE("is a runtime_error") {
    InvalidIndexException ex("test", 0, 0);
    CHECK_THROWS_AS(throw ex, std::runtime_error);
  }
}

TEST_CASE("NullPointerException") {
  SUBCASE("message format") {
    NullPointerException ex("spline_struct");
    std::string msg = ex.what();
    CHECK(msg.find("Null pointer") != std::string::npos);
    CHECK(msg.find("spline_struct") != std::string::npos);
  }

  SUBCASE("is a BmadException") {
    NullPointerException ex("test");
    CHECK_THROWS_AS(throw ex, BmadException);
  }
}

// ---------------------------------------------------------------------------
// FortranProxy null pointer rejection
// ---------------------------------------------------------------------------

TEST_CASE("FortranProxy rejects null pointer") {
  // Constructing a proxy from nullptr should throw NullPointerException
  CHECK_THROWS_AS(SplineStruct(nullptr, false), NullPointerException);
  CHECK_THROWS_AS(CoordStruct(nullptr, false), NullPointerException);
  CHECK_THROWS_AS(TestSubStruct(nullptr, false), NullPointerException);
}

// ---------------------------------------------------------------------------
// FortranProxy ownership and validity
// ---------------------------------------------------------------------------

TEST_CASE("FortranProxy default construction") {
  // Default construction should allocate and own memory
  auto spline = SplineStruct();
  CHECK(spline.is_valid());
  CHECK(spline.owns_memory());
  CHECK(spline.get_fortran_ptr() != nullptr);
}

TEST_CASE("FortranProxy move semantics") {
  auto original = SplineStruct();
  void *original_ptr = original.get_fortran_ptr();
  CHECK(original.is_valid());

  // Move should transfer ownership
  auto moved = std::move(original);
  CHECK(moved.is_valid());
  CHECK(moved.get_fortran_ptr() == original_ptr);
  CHECK(moved.owns_memory());
}

TEST_CASE("FortranProxy clone independence") {
  auto original = SplineStruct();
  original.set_x0(42.0);

  auto cloned = original.clone();
  CHECK(cloned.is_valid());
  CHECK(cloned.owns_memory());
  // Clone should be at a different address
  CHECK(cloned.get_fortran_ptr() != original.get_fortran_ptr());

  // Modifying clone should not affect original
  cloned.set_x0(99.0);
  CHECK(original.x0() == doctest::Approx(42.0));
  CHECK(cloned.x0() == doctest::Approx(99.0));
}

TEST_CASE("FortranProxy copy constructor independence") {
  auto original = CoordStruct();
  original.vec()[0] = 1.0;

  CoordStruct copy(original);
  CHECK(copy.is_valid());
  CHECK(copy.owns_memory());
  CHECK(copy.get_fortran_ptr() != original.get_fortran_ptr());

  // Modifying copy should not affect original
  copy.vec()[0] = 99.0;
  CHECK(original.vec()[0] == doctest::Approx(1.0));
  CHECK(copy.vec()[0] == doctest::Approx(99.0));
}
