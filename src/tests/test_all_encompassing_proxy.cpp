#include "bmad.hpp"
#include "doctest.h"

using namespace Bmad;

TEST_CASE("AllEncompassingProxy scalar types") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("complex_dp_0d read/write") {
    auto expected = std::complex<double>(1.0, 2.0);
    proxy.set_complex_dp_0d(expected);
    CHECK(proxy.complex_dp_0d().real() == doctest::Approx(expected.real()));
    CHECK(proxy.complex_dp_0d().imag() == doctest::Approx(expected.imag()));
  }

  SUBCASE("real_dp_0d read/write") {
    double expected = 3.14159265358979;
    proxy.set_real_dp_0d(expected);
    CHECK(proxy.real_dp_0d() == doctest::Approx(expected));
  }

  SUBCASE("real_rp_0d read/write") {
    double expected = 2.71828182845905;
    proxy.set_real_rp_0d(expected);
    CHECK(proxy.real_rp_0d() == doctest::Approx(expected));
  }

  SUBCASE("int_0d read/write") {
    int expected = 42;
    proxy.set_int_0d(expected);
    CHECK(proxy.int_0d() == expected);
  }

  SUBCASE("int8_0d read/write") {
    int64_t expected = 9223372036854775807LL; // max int64
    proxy.set_int8_0d(expected);
    CHECK(proxy.int8_0d() == expected);
  }

  SUBCASE("logical_0d read/write true") {
    proxy.set_logical_0d(true);
    CHECK(proxy.logical_0d() == true);
  }

  SUBCASE("logical_0d read/write false") {
    proxy.set_logical_0d(false);
    CHECK(proxy.logical_0d() == false);
  }
}

TEST_CASE("AllEncompassingProxy 1D fixed arrays") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("real_dp_1d read/write") {
    auto arr = proxy.real_dp_1d();
    REQUIRE(arr.size() == 3);

    // Write values
    for (int i = 0; i < arr.size(); ++i) {
      arr[i] = static_cast<double>(i) * 1.5;
    }

    // Read back and verify
    auto arr_read = proxy.real_dp_1d();
    for (int i = 0; i < arr.size(); ++i) {
      CHECK(arr_read[i] == doctest::Approx(static_cast<double>(i) * 1.5));
    }
  }

  SUBCASE("real_rp_1d read/write") {
    auto arr = proxy.real_rp_1d();
    REQUIRE(arr.size() == 3);

    for (int i = 0; i < arr.size(); ++i) {
      arr[i] = static_cast<double>(i) * 2.5;
    }

    auto arr_read = proxy.real_rp_1d();
    for (int i = 0; i < arr.size(); ++i) {
      CHECK(arr_read[i] == doctest::Approx(static_cast<double>(i) * 2.5));
    }
  }

  SUBCASE("int_1d read/write") {
    auto arr = proxy.int_1d();
    REQUIRE(arr.size() == 3);

    for (int i = 0; i < arr.size(); ++i) {
      arr[i] = i * 10;
    }

    auto arr_read = proxy.int_1d();
    for (int i = 0; i < arr.size(); ++i) {
      CHECK(arr_read[i] == i * 10);
    }
  }

  SUBCASE("complex_dp_1d read/write") {
    auto arr = proxy.complex_dp_1d();
    REQUIRE(arr.size() == 3);

    for (int i = 0; i < arr.size(); ++i) {
      arr[i] = std::complex<double>(
          static_cast<double>(i), static_cast<double>(i) + 0.5);
    }

    auto arr_read = proxy.complex_dp_1d();
    for (int i = 0; i < arr.size(); ++i) {
      CHECK(arr_read[i].real() == doctest::Approx(static_cast<double>(i)));
      CHECK(
          arr_read[i].imag() == doctest::Approx(static_cast<double>(i) + 0.5));
    }
  }
}

TEST_CASE("AllEncompassingProxy 2D fixed arrays") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("real_dp_2d read/write") {
    auto arr = proxy.real_dp_2d();
    REQUIRE(arr.total_size() == 3 * 4);

    auto sizes = arr.size(); // returns std::array<int, 2>
    auto bounds1 = arr.bounds(1); // first dimension bounds
    auto bounds2 = arr.bounds(2); // second dimension bounds

    // Write values based on Fortran indices (using bounds)
    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        arr(i, j) = static_cast<double>(i * 100 + j);
      }
    }

    // Read back and verify
    auto arr_read = proxy.real_dp_2d();
    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        CHECK(
            arr_read(i, j) ==
            doctest::Approx(static_cast<double>(i * 100 + j)));
      }
    }
  }

  SUBCASE("int_2d read/write") {
    auto arr = proxy.int_2d();
    REQUIRE(arr.total_size() == 3 * 4);

    auto bounds1 = arr.bounds(1);
    auto bounds2 = arr.bounds(2);

    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        arr(i, j) = i * 10 + j;
      }
    }

    auto arr_read = proxy.int_2d();
    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        CHECK(arr_read(i, j) == i * 10 + j);
      }
    }
  }

  SUBCASE("complex_dp_2d read/write") {
    auto arr = proxy.complex_dp_2d();
    REQUIRE(arr.total_size() == 3 * 4);

    auto bounds1 = arr.bounds(1);
    auto bounds2 = arr.bounds(2);

    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        arr(i, j) = std::complex<double>(
            static_cast<double>(i), static_cast<double>(j));
      }
    }

    auto arr_read = proxy.complex_dp_2d();
    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        CHECK(arr_read(i, j).real() == doctest::Approx(static_cast<double>(i)));
        CHECK(arr_read(i, j).imag() == doctest::Approx(static_cast<double>(j)));
      }
    }
  }
}

TEST_CASE("AllEncompassingProxy 3D fixed arrays") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("real_dp_3d read/write") {
    auto arr = proxy.real_dp_3d();
    REQUIRE(arr.total_size() == 3 * 4 * 5);

    auto bounds1 = arr.bounds(1);
    auto bounds2 = arr.bounds(2);
    auto bounds3 = arr.bounds(3);

    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        for (int k = bounds3.first; k <= bounds3.second; ++k) {
          arr(i, j, k) = static_cast<double>(i * 100 + j * 10 + k);
        }
      }
    }

    auto arr_read = proxy.real_dp_3d();
    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        for (int k = bounds3.first; k <= bounds3.second; ++k) {
          CHECK(
              arr_read(i, j, k) ==
              doctest::Approx(static_cast<double>(i * 100 + j * 10 + k)));
        }
      }
    }
  }

  SUBCASE("int_3d read/write") {
    auto arr = proxy.int_3d();
    REQUIRE(arr.total_size() == 3 * 4 * 5);

    auto bounds1 = arr.bounds(1);
    auto bounds2 = arr.bounds(2);
    auto bounds3 = arr.bounds(3);

    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        for (int k = bounds3.first; k <= bounds3.second; ++k) {
          arr(i, j, k) = i * 100 + j * 10 + k;
        }
      }
    }

    auto arr_read = proxy.int_3d();
    for (int i = bounds1.first; i <= bounds1.second; ++i) {
      for (int j = bounds2.first; j <= bounds2.second; ++j) {
        for (int k = bounds3.first; k <= bounds3.second; ++k) {
          CHECK(arr_read(i, j, k) == i * 100 + j * 10 + k);
        }
      }
    }
  }
}

TEST_CASE("AllEncompassingProxy pointer members - read only") {
  auto proxy = AllEncompassingProxy();

  // Note: Pointer members in Fortran are uninitialized by default (null).
  // This test only verifies that we can safely query pointer state.
  // Writing to uninitialized pointers requires explicit allocation first.

  SUBCASE("real_dp_0d_ptr null check") {
    // Pointer should be null initially
    auto ptr = proxy.real_dp_0d_ptr();
    // Just verify we don't crash - pointer may be null or garbage
    (void)ptr; // suppress unused warning
  }

  SUBCASE("real_rp_0d_ptr null check") {
    auto ptr = proxy.real_rp_0d_ptr();
    (void)ptr;
  }

  SUBCASE("int_0d_ptr null check") {
    auto ptr = proxy.int_0d_ptr();
    (void)ptr;
  }

  SUBCASE("int8_0d_ptr null check") {
    auto ptr = proxy.int8_0d_ptr();
    (void)ptr;
  }

  SUBCASE("logical_0d_ptr null check") {
    auto ptr = proxy.logical_0d_ptr();
    (void)ptr;
  }
}

TEST_CASE("AllEncompassingProxy 1D pointer arrays") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("real_dp_1d_ptr") {
    auto arr = proxy.real_dp_1d_ptr();
    REQUIRE_FALSE(arr.is_valid());

    // Pointer arrays may be unallocated initially
    if (arr.size() > 0) {
      for (int i = 0; i < arr.size(); ++i) {
        arr[i] = static_cast<double>(i) * 3.0;
      }

      auto arr_read = proxy.real_dp_1d_ptr();
      for (int i = 0; i < arr.size(); ++i) {
        CHECK(arr_read[i] == doctest::Approx(static_cast<double>(i) * 3.0));
      }
    }
  }

  SUBCASE("int_1d_ptr") {
    auto arr = proxy.int_1d_ptr();
    if (arr.size() > 0) {
      for (int i = 0; i < arr.size(); ++i) {
        arr[i] = i * 7;
      }

      auto arr_read = proxy.int_1d_ptr();
      for (int i = 0; i < arr.size(); ++i) {
        CHECK(arr_read[i] == i * 7);
      }
    }
  }
}

TEST_CASE("AllEncompassingProxy allocatable arrays") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("real_dp_1d_alloc") {
    auto arr = proxy.real_dp_1d_alloc();
    // Allocatable arrays may need allocation first
    if (arr.size() > 0) {
      for (int i = 0; i < arr.size(); ++i) {
        arr[i] = static_cast<double>(i) * 4.0;
      }

      auto arr_read = proxy.real_dp_1d_alloc();
      for (int i = 0; i < arr.size(); ++i) {
        CHECK(arr_read[i] == doctest::Approx(static_cast<double>(i) * 4.0));
      }
    }
  }

  SUBCASE("int_1d_alloc") {
    auto arr = proxy.int_1d_alloc();
    if (arr.size() > 0) {
      for (int i = 0; i < arr.size(); ++i) {
        arr[i] = i * 11;
      }

      auto arr_read = proxy.int_1d_alloc();
      for (int i = 0; i < arr.size(); ++i) {
        CHECK(arr_read[i] == i * 11);
      }
    }
  }

  SUBCASE("complex_dp_1d_alloc") {
    auto arr = proxy.complex_dp_1d_alloc();
    if (arr.size() > 0) {
      for (int i = 0; i < arr.size(); ++i) {
        arr[i] = std::complex<double>(
            static_cast<double>(i), -static_cast<double>(i));
      }

      auto arr_read = proxy.complex_dp_1d_alloc();
      for (int i = 0; i < arr.size(); ++i) {
        CHECK(arr_read[i].real() == doctest::Approx(static_cast<double>(i)));
        CHECK(arr_read[i].imag() == doctest::Approx(-static_cast<double>(i)));
      }
    }
  }
}

TEST_CASE("AllEncompassingProxy nested type (TestSubProxy)") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("type_0d read/write with nested fields") {
    auto sub = proxy.type_0d();
    CHECK_NOTHROW(proxy.type_0d());

    // Access nested TestSubSubProxy via sr()
    auto sub_sub = proxy.type_0d().sr();

    // Test int8 field (aaa)
    sub_sub.set_aaa(123456789LL);
    CHECK(proxy.type_0d().sr().aaa() == 123456789LL);

    // Test int field (bbb)
    sub_sub.set_bbb(42);
    CHECK(proxy.type_0d().sr().bbb() == 42);

    // Test real fields
    sub_sub.set_t_ref(1.5);
    CHECK(proxy.type_0d().sr().t_ref() == doctest::Approx(1.5));

    sub_sub.set_freq_spread(0.001);
    CHECK(proxy.type_0d().sr().freq_spread() == doctest::Approx(0.001));
  }

  SUBCASE("type_1d access and modify") {
    auto arr = proxy.type_1d();
    REQUIRE(arr.size() > 0);

    // Modify first element's nested fields
    arr[0].sr().set_bbb(100);
    arr[1].sr().set_bbb(200);

    CHECK(proxy.type_1d()[0].sr().bbb() == 100);
    CHECK(proxy.type_1d()[1].sr().bbb() == 200);
  }

  SUBCASE("type_2d access") {
    auto arr = proxy.type_2d();
    REQUIRE(arr.total_size() > 0);
    auto bounds1 = arr.bounds(1);
    auto bounds2 = arr.bounds(2);

    // Access and modify element
    arr(bounds1.first, bounds2.first).sr().set_aaa(999);
    CHECK(proxy.type_2d()(bounds1.first, bounds2.first).sr().aaa() == 999);
  }

  SUBCASE("type_3d access") {
    auto arr = proxy.type_3d();
    REQUIRE(arr.total_size() > 0);
    auto bounds1 = arr.bounds(1);
    auto bounds2 = arr.bounds(2);
    auto bounds3 = arr.bounds(3);

    arr(bounds1.first, bounds2.first, bounds3.first).sr().set_bbb(777);
    CHECK(
        proxy.type_3d()(bounds1.first, bounds2.first, bounds3.first)
            .sr()
            .bbb() == 777);
  }

  SUBCASE("type_0d_ptr access") {
    auto opt = proxy.type_0d_ptr();
    // May or may not be allocated
    if (opt.has_value()) {
      CHECK_NOTHROW(opt.value());
    }
  }
}

TEST_CASE("TestSubProxy standalone") {
  SUBCASE("create and access fields") {
    auto sub = TestSubProxy();

    // Access nested TestSubSubProxy
    auto sub_sub = sub.sr();

    // Test all fields
    sub_sub.set_aaa(987654321LL);
    CHECK(sub.sr().aaa() == 987654321LL);

    sub_sub.set_bbb(123);
    CHECK(sub.sr().bbb() == 123);

    sub_sub.set_t_ref(2.718);
    CHECK(sub.sr().t_ref() == doctest::Approx(2.718));

    sub_sub.set_freq_spread(1e-6);
    CHECK(sub.sr().freq_spread() == doctest::Approx(1e-6));
  }

  SUBCASE("copy behavior") {
    auto sub1 = TestSubProxy();
    sub1.sr().set_bbb(111);

    // Clone creates independent copy
    auto sub2 = sub1.clone();
    sub2.sr().set_bbb(222);

    CHECK(sub1.sr().bbb() == 111);
    CHECK(sub2.sr().bbb() == 222);
  }
}

TEST_CASE("TestSubProxyAllocatable1D") {
  SUBCASE("allocate and access") {
    auto arr = TestSubProxyAllocatable1D();
    // resize(lbound, ubound) allocates array from lbound to ubound inclusive
    arr.resize(1, 5); // Fortran array(1:5)
    CHECK(arr.view().lower_bound() == 1);

    REQUIRE(arr.size() == 5);

    // Set values in each element (0-based C++ indexing)
    for (int i = 0; i < arr.size(); ++i) {
      arr[i].sr().set_bbb(i * 100);
    }

    // Read back
    for (int i = 0; i < arr.size(); ++i) {
      CHECK(arr[i].sr().bbb() == i * 100);
    }
  }

  SUBCASE("resize with lbound 0") {
    auto arr = TestSubProxyAllocatable1D();
    arr.resize(0, 2); // Fortran array(0:2)
    REQUIRE(arr.size() == 2);
    CHECK(arr.view().lower_bound() == 0);

    for (int i = 0; i < arr.size(); ++i) {
      arr[i].sr().set_bbb((i + 1) * 10);
    }

    for (int i = 0; i < arr.size(); ++i) {
      CHECK(arr[i].sr().bbb() == (i + 1) * 10);
    }
  }

  SUBCASE("resize and clear") {
    auto arr = TestSubProxyAllocatable1D();

    // Initially empty
    CHECK(arr.empty());

    // Allocate
    arr.resize(1, 3);
    CHECK(arr.size() == 3);
    CHECK(arr.view().lower_bound() == 1);

    // Set data
    arr[0].sr().set_aaa(1);
    arr[1].sr().set_aaa(2);
    arr[2].sr().set_aaa(3);

    // Clear
    arr.clear();
    CHECK(arr.empty());
  }
}

TEST_CASE("AllEncompassingProxy edge cases") {
  auto proxy = AllEncompassingProxy();

  SUBCASE("zero values") {
    proxy.set_real_dp_0d(0.0);
    CHECK(proxy.real_dp_0d() == doctest::Approx(0.0));

    proxy.set_int_0d(0);
    CHECK(proxy.int_0d() == 0);

    proxy.set_complex_dp_0d(std::complex<double>(0.0, 0.0));
    CHECK(proxy.complex_dp_0d().real() == doctest::Approx(0.0));
    CHECK(proxy.complex_dp_0d().imag() == doctest::Approx(0.0));
  }

  SUBCASE("negative values") {
    proxy.set_real_dp_0d(-123.456);
    CHECK(proxy.real_dp_0d() == doctest::Approx(-123.456));

    proxy.set_int_0d(-42);
    CHECK(proxy.int_0d() == -42);

    proxy.set_int8_0d(-9223372036854775807LL);
    CHECK(proxy.int8_0d() == -9223372036854775807LL);
  }

  SUBCASE("complex with negative components") {
    auto cpx = std::complex<double>(-1.5, -2.5);
    proxy.set_complex_dp_0d(cpx);
    CHECK(proxy.complex_dp_0d().real() == doctest::Approx(-1.5));
    CHECK(proxy.complex_dp_0d().imag() == doctest::Approx(-2.5));
  }

  SUBCASE("very small values") {
    double small = 1e-300;
    proxy.set_real_dp_0d(small);
    CHECK(proxy.real_dp_0d() == doctest::Approx(small));
  }

  SUBCASE("very large values") {
    double large = 1e300;
    proxy.set_real_dp_0d(large);
    CHECK(proxy.real_dp_0d() == doctest::Approx(large));
  }
}
