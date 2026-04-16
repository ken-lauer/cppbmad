#include <bmad.hpp>

#include "doctest.h"

using namespace Bmad;

// ---------------------------------------------------------------------------
// Basic type overloads
// ---------------------------------------------------------------------------

TEST_CASE("to_string - complex") {
  SUBCASE("positive imaginary") {
    auto s = to_string(std::complex<double>(1.5, 2.0));
    CHECK(s.find("1.5") != std::string::npos);
    CHECK(s.find("2") != std::string::npos);
    CHECK(s.find("i") != std::string::npos);
    CHECK(s.find("+") != std::string::npos);
  }

  SUBCASE("negative imaginary") {
    auto s = to_string(std::complex<double>(1.0, -3.0));
    CHECK(s.find("1") != std::string::npos);
    CHECK(s.find("-3") != std::string::npos);
    CHECK(s.find("i") != std::string::npos);
    // Should not have '+' before negative imaginary
    CHECK(s.find("+-") == std::string::npos);
  }

  SUBCASE("zero") {
    auto s = to_string(std::complex<double>(0.0, 0.0));
    CHECK(s.find("i") != std::string::npos);
  }
}

TEST_CASE("to_string - optional") {
  SUBCASE("has value") {
    std::optional<double> opt(3.14);
    auto s = to_string(opt);
    CHECK(s.find("3.14") != std::string::npos);
  }

  SUBCASE("nullopt") {
    std::optional<double> opt;
    CHECK(to_string(opt) == "nullopt");
  }
}

TEST_CASE("to_string - void pointer") {
  SUBCASE("nullptr") {
    void *ptr = nullptr;
    CHECK(to_string(ptr) == "nullptr");
  }

  SUBCASE("non-null") {
    int x = 42;
    void *ptr = &x;
    auto s = to_string(ptr);
    CHECK(s != "nullptr");
  }
}

TEST_CASE("to_string - typed pointer") {
  SUBCASE("nullptr") {
    double *ptr = nullptr;
    CHECK(to_string(ptr) == "nullptr");
  }

  SUBCASE("non-null") {
    double x = 2.5;
    double *ptr = &x;
    auto s = to_string(ptr);
    CHECK(s.find("2.5") != std::string::npos);
  }
}

// ---------------------------------------------------------------------------
// STL containers
// ---------------------------------------------------------------------------

TEST_CASE("to_string - vector") {
  SUBCASE("empty") {
    std::vector<int> v;
    CHECK(to_string(v) == "[]");
  }

  SUBCASE("single element") {
    std::vector<int> v = {42};
    CHECK(to_string(v) == "[42]");
  }

  SUBCASE("multiple elements") {
    std::vector<int> v = {1, 2, 3};
    auto s = to_string(v);
    CHECK(s.front() == '[');
    CHECK(s.back() == ']');
    CHECK(s.find("1") != std::string::npos);
    CHECK(s.find("2") != std::string::npos);
    CHECK(s.find("3") != std::string::npos);
  }
}

TEST_CASE("to_string - array") {
  std::array<int, 3> arr = {10, 20, 30};
  auto s = to_string(arr);
  CHECK(s.front() == '[');
  CHECK(s.back() == ']');
  CHECK(s.find("10") != std::string::npos);
  CHECK(s.find("20") != std::string::npos);
  CHECK(s.find("30") != std::string::npos);
  CHECK(s.find("30") > s.find("20"));
  CHECK(s.find("20") > s.find("10"));
}

// ---------------------------------------------------------------------------
// repr() helper
// ---------------------------------------------------------------------------

TEST_CASE("repr helper") {
  SUBCASE("with fields") {
    int dummy = 0;
    auto s = repr(&dummy, "MyClass", {{"x", "1"}, {"y", "hello"}});
    CHECK(s.find("MyClass(0x") != std::string::npos);
    CHECK(s.find("x=1") != std::string::npos);
    CHECK(s.find("y=hello") != std::string::npos);
    CHECK(s.back() == ')');
  }

  SUBCASE("no fields") {
    int dummy = 0;
    auto s = repr(&dummy, "Empty", {});
    CHECK(s.find("Empty(0x") != std::string::npos);
    CHECK(s.back() == ')');
  }

  SUBCASE("null pointer") {
    auto s = repr(nullptr, "NullObj", {{"a", "1"}});
    CHECK(s.find("NullObj(0x0") != std::string::npos);
  }
}

// ---------------------------------------------------------------------------
// Struct to_string (generated)
// ---------------------------------------------------------------------------

TEST_CASE("to_string - SplineStruct") {
  auto spline = SplineStruct();
  spline.set_x0(1.0);
  spline.set_y0(2.0);

  auto s = to_string(spline);
  CHECK(s.find("SplineStruct(0x") != std::string::npos);
  CHECK(s.find("x0=") != std::string::npos);
  CHECK(s.find("y0=") != std::string::npos);
  CHECK(s.back() == ')');
}

TEST_CASE("to_string - CoordStruct") {
  auto coord = CoordStruct();
  auto s = to_string(coord);
  CHECK(s.find("CoordStruct(0x") != std::string::npos);
  CHECK(s.find("vec=") != std::string::npos);
  CHECK(s.back() == ')');
}

TEST_CASE("to_string - AllEncompassingStruct") {
  auto aes = AllEncompassingStruct();
  aes.set_real_dp_0d(3.14);
  aes.set_int_0d(42);

  auto s = to_string(aes);
  CHECK(s.find("AllEncompassingStruct(0x") != std::string::npos);
  CHECK(s.find("real_dp_0d=") != std::string::npos);
  CHECK(s.find("int_0d=") != std::string::npos);
  CHECK(s.back() == ')');
}

TEST_CASE("to_string - TestSubStruct with nesting") {
  auto sub = TestSubStruct();
  sub.sr().set_bbb(99);

  auto s = to_string(sub);
  CHECK(s.find("TestSubStruct(0x") != std::string::npos);
  // Should contain nested TestSubSubStruct repr
  CHECK(s.find("sr=") != std::string::npos);
}

// ---------------------------------------------------------------------------
// Fortran array to_string
// ---------------------------------------------------------------------------

TEST_CASE("to_string - FArray1D via proxy") {
  auto aes = AllEncompassingStruct();
  auto arr = aes.real_dp_1d();
  arr[0] = 1.5;
  arr[1] = 2.5;
  arr[2] = 3.5;

  auto s = to_string(arr);
  CHECK(s.front() == '[');
  CHECK(s.back() == ']');
  CHECK(s.find("1.5") != std::string::npos);
  CHECK(s.find("2.5") != std::string::npos);
  CHECK(s.find("3.5") != std::string::npos);
}

TEST_CASE("to_string - FAlloc1D via container") {
  auto container = RealAlloc1D();
  container.resize(3);
  auto view = container.view();
  view[0] = 10.0;
  view[1] = 20.0;
  view[2] = 30.0;

  auto s = to_string(container);
  CHECK(s.front() == '[');
  CHECK(s.back() == ']');
  CHECK(s.find("10") != std::string::npos);
  CHECK(s.find("20") != std::string::npos);
  CHECK(s.find("30") != std::string::npos);
}
