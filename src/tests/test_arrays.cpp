#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "bmad/fortran_arrays.hpp" // Assumes the header is here
#include "doctest.h"

#include <numeric>
#include <vector>

// =============================================================================
// MOCKING INFRASTRUCTURE
// These simulate the Fortran runtime allocation/descriptor behavior
// =============================================================================

namespace Mocks {

// Global counters to verify RAII / Deletion logic
static int g_alloc_count = 0;
static int g_dealloc_count = 0;

void reset_counters() {
  g_alloc_count = 0;
  g_dealloc_count = 0;
}

// -----------------------------------------------------------------------------
// Mock: FAlloc1D Handle (Simulating a Fortran Descriptor)
// -----------------------------------------------------------------------------
struct IntDescriptor {
  std::vector<int> buffer;
  int cur_lbound = 1;

  void* data_ptr() {
    return buffer.empty() ? nullptr : buffer.data();
  }
  int size() const {
    return static_cast<int>(buffer.size());
  }
};

void* Alloc_IntDesc() {
  g_alloc_count++;
  return new IntDescriptor();
}

void Dealloc_IntDesc(void* ptr) {
  if (ptr) {
    g_dealloc_count++;
    delete static_cast<IntDescriptor*>(ptr);
  }
}

void Realloc_IntDesc(void* ptr, int lb, size_t n) {
  auto* desc = static_cast<IntDescriptor*>(ptr);
  desc->cur_lbound = lb;
  desc->buffer.resize(n);
  // Fill with predictable data for testing: i * 10
  for (size_t i = 0; i < n; ++i)
    desc->buffer[i] = static_cast<int>(i) * 10;
}

void Access_IntDesc(
    void* ptr,
    void** d,
    int* lb,
    int* sz,
    size_t* es,
    bool* alloc) {
  auto* desc = static_cast<IntDescriptor*>(ptr);
  *d = desc->data_ptr();
  *lb = desc->cur_lbound;
  *sz = desc->size();
  *es = sizeof(int);
  *alloc = (*d != nullptr);
}

// -----------------------------------------------------------------------------
// Mock: FTypeArray1D Lifecycle (Simulating Derived Type allocators)
// -----------------------------------------------------------------------------

struct MyType {
  int val;
  void* raw_addr; // Store access address to verify pointer arithmetic
};

// A "Proxy" class that users would write to wrap the raw void*
class MyTypeProxy {
  void* ptr_;

 public:
  explicit MyTypeProxy(void* p) : ptr_(p) {}

  // Getter simply treats the void* as an int* for this test
  int get_val() const {
    return *static_cast<int*>(ptr_);
  }
  void set_val(int v) {
    *static_cast<int*>(ptr_) = v;
  }
  void* raw() const {
    return ptr_;
  }
};

void* Alloc_Type(int n, size_t* elem_size) {
  g_alloc_count++;
  *elem_size = sizeof(int);
  return new int[n]; // Simple int array to simulate a structure
}

void Dealloc_Type(void* ptr, int /*size_context*/) {
  if (ptr) {
    g_dealloc_count++;
    delete[] static_cast<int*>(ptr);
  }
}

} // namespace Mocks

using namespace Bmad;

// =============================================================================
// TEST CASES
// =============================================================================

TEST_CASE("FArrayND<Primitive, 1> : Explicit 1D Specialization") {
  // Setup raw data: [100, 101, 102, 103, 104]
  std::vector<int> raw_data(5);
  std::iota(raw_data.begin(), raw_data.end(), 100);

  // Scenario: Array runs from index 10 to 14 (size 5)
  FArray1D<int> arr(raw_data.data(), 5, 10, 14, true);

  SUBCASE("Metadata is correct") {
    CHECK(arr.is_valid());
    CHECK(arr.size() == 5);
    CHECK(arr.lower_bound() == 10);
    CHECK(arr.upper_bound() == 14);
    CHECK(arr.rank() == 1);
    CHECK(!arr.empty());
  }

  SUBCASE("Access via Fortran Indexing (Curve-based)") {
    CHECK(arr(10) == 100); // lower bound
    CHECK(arr(14) == 104); // upper bound
    CHECK_THROWS_AS(arr(9), std::out_of_range);
    CHECK_THROWS_AS(arr(15), std::out_of_range);
  }

  SUBCASE("Access via C Indexing (0-based)") {
    CHECK(arr[0] == 100);
    CHECK(arr[4] == 104);
    CHECK(arr.at(0) == 100);
    CHECK_THROWS_AS(arr[-1], std::out_of_range);
    CHECK_THROWS_AS(arr[5], std::out_of_range);
  }

  SUBCASE("Iteration") {
    int sum = 0;
    for (int val : arr)
      sum += val;
    CHECK(sum == 510);
  }

  SUBCASE("Vector conversion") {
    auto vec = arr.to_vector();
    CHECK(vec.size() == 5);
    CHECK(vec[0] == 100);
  }
}

TEST_CASE("FArrayND<T, 2> : Multidimensional Array") {
  // 2x3 Array (Row-major linear storage usually, but let's define manual strides)
  // Fortran is Column-Major.
  // Let's pretend we have a 2x3 array.
  // (1,1) (1,2) (1,3)
  // (2,1) (2,2) (2,3)
  //
  // Storing lineary as 1,2,3,4,5,6
  std::vector<int> raw = {1, 2, 3, 4, 5, 6};

  // Configuring strides for standard C layout (Row Major) to make math easy in head
  // Dims: 2, 3.
  // Stride for dim 2 (cols) = 1.
  // Stride for dim 1 (rows) = 3.

  // Bounds: 1:2, 1:3
  FArray2D<int> arr(
      raw.data(), // ptr
      2,
      1,
      2, // Dim 1: size=2, lb=1, ub=2
      3,
      1,
      3, // Dim 2: size=3, lb=1, ub=3
      3,
      1, // Strides: dim1_stride=3, dim2_stride=1 (Row Major simulation)
      true // valid
  );

  SUBCASE("Metadata") {
    CHECK(arr.size(1) == 2);
    CHECK(arr.size(2) == 3);
    CHECK(arr.total_size() == 6);
  }

  SUBCASE("Coordinate access") {
    // (1,1) -> index 0 -> 1
    CHECK(arr(1, 1) == 1);

    // (1,2) -> row 1 is offset 0, col 2 is offset 1 -> index 1 -> 2
    CHECK(arr(1, 2) == 2);

    // (2,1) -> row 2 is offset 3, col 1 is offset 0 -> index 3 -> 4
    CHECK(arr(2, 1) == 4);

    // Access via C reference
    CHECK(arr.at(0, 0) == 1); // C index 0,0
    CHECK(arr.at(1, 0) == 4); // C index 1,0
  }

  SUBCASE("Out of bounds") {
    CHECK_THROWS(arr(3, 1)); // Row too high
    CHECK_THROWS(arr(1, 4)); // Col too high
  }
}

TEST_CASE("FCharArray1D : String helper") {
  // 3 strings of length 5 packed together
  // "Hello", "War  ", "World" (packed tight for test)
  std::string buffer = "HelloWar  World";
  // Note: buffer is const char, cast const away for the view, but we won't write to const memory
  std::vector<char> writable_buf(buffer.begin(), buffer.end());

  int n_strings = 3;
  int str_len = 5;

  FCharArray1D strings(
      writable_buf.data(),
      n_strings,
      1,
      n_strings, // Bounds 1:3
      str_len,
      true);

  SUBCASE("Reading strings trims whitespace") {
    CHECK(strings(1).str() == "Hello");
    CHECK(strings(2).str() == "War"); // "War  " -> "War"
    CHECK(strings(3).str() == "World");
  }

  SUBCASE("Implicit string conversion") {
    std::string s = strings(1);
    CHECK(s == "Hello");
  }

  SUBCASE("Writing strings handles padding") {
    strings(1) = "Hi"; // Should become "Hi   " internaly
    CHECK(writable_buf[0] == 'H');
    CHECK(writable_buf[1] == 'i');
    CHECK(writable_buf[2] == ' ');
    CHECK(strings(1).str() == "Hi");
  }

  SUBCASE("Iteration") {
    auto it = strings.begin();
    CHECK((*it).str() == "Hello");
    ++it;
    CHECK((*it).str() == "War");
  }
}

TEST_CASE("FTypeArrayND<T, 1> : Ownership and Lifecycle") {
  Mocks::reset_counters();

  // 1. Test "Owned" array (created via .allocate)
  {
    // Allocate array of 10 items, lower bound 5
    auto arr = FTypeArray1D<
        Mocks::MyTypeProxy,
        Mocks::Alloc_Type,
        Mocks::Dealloc_Type>::allocate(10, 5);

    CHECK(Mocks::g_alloc_count == 1);
    CHECK(arr.size() == 10);
    CHECK(arr.lower_bound() == 5);
    CHECK(arr.upper_bound() == 14);

    // Access/Modify
    arr(5).set_val(999);
    CHECK(arr(5).get_val() == 999);

    // C-style access
    arr[0].set_val(111);
    CHECK(arr[0].get_val() == 111); // Same memory loc as arr(5)

    // Address consistency
    CHECK(arr.at(0).raw() == arr(5).raw());
  }
  // Destructor should have successfully fired DeallocFunc here (or shared_ptr deleter)
  CHECK(Mocks::g_dealloc_count == 1);

  // 2. Test "View" array (Non-owning)
  int raw_int = 42;
  {
    // Create view manually
    FTypeArray1D<Mocks::MyTypeProxy, nullptr, nullptr> view(
        &raw_int, 1, 1, 1, true, sizeof(int));
    CHECK(view(1).get_val() == 42);
  }
  // View destruction should NOT trigger mock dealloc
  CHECK(Mocks::g_dealloc_count == 1); // logic check: still 1 from previous test
}

TEST_CASE("FAlloc1D : Allocatable Container Refactoring Support") {
  Mocks::reset_counters();

  // Instantiate container. Constructor calls AllocFunc immediately to create handle.
  Bmad::FAlloc1D<
      int,
      Mocks::Alloc_IntDesc,
      Mocks::Dealloc_IntDesc,
      Mocks::Realloc_IntDesc,
      Mocks::Access_IntDesc>
      container;

  CHECK(Mocks::g_alloc_count == 1); // Handle created
  CHECK(container.empty());

  // Resize (Allocation)
  container.resize(1, 5);
  CHECK(!container.empty());
  CHECK(container.size() == 5);

  // Check internal View generation
  // logic in Mock: fill with i*10 -> 0, 10, 20...
  CHECK(container(1) == 0);
  CHECK(container(3) == 20); // 0, 10, 20
  CHECK(container[2] == 20);

  // Modify data
  container(1) = 999;
  CHECK(container(1) == 999);

  // Clear
  container.clear();
  CHECK(container.empty()); // View size should be 0 or !valid

  // Scope exit: destructor should clean up handle
}
// g_dealloc_count check outside scope
TEST_CASE("FAlloc1D : RAII Check") {
  // Relying on previous test running, or simpler:
  Mocks::reset_counters();
  {
    Bmad::FAlloc1D<
        int,
        Mocks::Alloc_IntDesc,
        Mocks::Dealloc_IntDesc,
        Mocks::Realloc_IntDesc,
        Mocks::Access_IntDesc>
        temp;
    // Created
  }
  // Destroyed
  CHECK(Mocks::g_dealloc_count == 1);
}
