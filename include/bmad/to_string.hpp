#pragma once
#include <array>
#include <complex>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "bmad/fortran_arrays.hpp"
namespace Bmad {
using std::complex;
using std::string;
using std::to_string; // Import standard to_string for basic types

// -----------------------------------------------------------------------------
// Forward declarations for STL containers to allow recursive nesting
// -----------------------------------------------------------------------------
template <typename T>
string to_string(const std::vector<T>& vec);

template <typename T, std::size_t N>
string to_string(const std::array<T, N>& arr);

// -----------------------------------------------------------------------------
// Existing Arithmetic / Optional Overloads
// -----------------------------------------------------------------------------

template <typename T>
string to_string(const complex<T>&);

template <typename T>
string to_string(const complex<T>& c) {
  std::ostringstream oss;
  oss << c.real();
  if (c.imag() >= 0) {
    oss << "+";
  }
  oss << c.imag() << "i";
  return oss.str();
}

template <typename T>
std::string to_string(const std::optional<T>& opt) {
  if (opt.has_value()) {
    return to_string(opt.value());
  }
  return "nullopt";
}

template <typename T>
std::string to_string(const T* opt) {
  if (opt) {
    return to_string(*opt);
  }
  return "nullptr";
}

template string to_string(const complex<double>&);

std::string repr(
    const void* obj,
    const std::string& class_name,
    const std::initializer_list<std::pair<std::string, std::string>>& fields);

// Explicit implementation for std::vector (VariableArrayXD)
template <typename T>
string to_string(const std::vector<T>& vec) {
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i > 0)
      oss << ", ";
    // This recursion works because to_string is overloaded for the contained type T
    oss << to_string(vec[i]);
  }
  oss << "]";
  return oss.str();
}

// Explicit implementation for std::array (FixedArrayXD)
template <typename T, std::size_t N>
string to_string(const std::array<T, N>& arr) {
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < N; ++i) {
    if (i > 0)
      oss << ", ";
    oss << to_string(arr[i]);
  }
  oss << "]";
  return oss.str();
}

// =============================================================================
// to_string overloads for Fortran arrays
// =============================================================================

std::string to_string(const FCharArray1D& arr);
namespace detail {
// Helper to unpack std::array indices into variadic arguments for FArrayND::at()
template <typename T, std::size_t N, std::size_t... Is>
const T& farray_get_packed(
    const FArrayND<T, N>& arr,
    const std::array<int, N>& idxs,
    std::index_sequence<Is...>) {
  return arr.at(idxs[Is]...);
}
} // namespace detail

template <typename T, std::size_t N>
std::string to_string(const FArrayND<T, N>& arr) {
  if (!arr.is_valid())
    return "[]";

  std::ostringstream oss;

  // Optimized path for 1D arrays (most common usage)
  if constexpr (N == 1) {
    oss << "[";
    // arr.size() returns std::array<int, 1>, we need the first element
    int len = arr.size();
    for (int i = 0; i < len; ++i) {
      if (i > 0)
        oss << ", ";
      // Use .at(i) because operator[] is not available on FArrayND
      oss << to_string(arr.at(i));
    }
    oss << "]";
  }
  // General path for N-dimensional arrays
  else {
    // Current state of indices for recursion
    std::array<int, N> indices{};

    // Recursive lambda to handle printing nested dimensions
    // Uses auto&& self to allow recursion within the lambda
    auto print_recursive = [&](auto&& self, size_t current_dim) -> void {
      int dim_size = arr.size()[current_dim];
      oss << "[";

      for (int i = 0; i < dim_size; ++i) {
        if (i > 0)
          oss << ", ";

        indices[current_dim] = i;

        if (current_dim == N - 1) {
          // Base case: We are at the deepest level (Is Leaf).
          // Unpack the 'indices' array into arguments for arr.at(i, j, k...)
          oss << to_string(
              detail::farray_get_packed(
                  arr, indices, std::make_index_sequence<N>{}));
        } else {
          // Recursive step: Go deeper
          self(self, current_dim + 1);
        }
      }
      oss << "]";
    };

    // Start recursion at dimension 0
    print_recursive(print_recursive, 0);
  }

  return oss.str();
}

template <typename ProxyType, std::size_t N, auto AllocFunc, auto DeallocFunc>
std::string to_string(
    const FTypeArrayND<ProxyType, N, AllocFunc, DeallocFunc>& arr) {
  std::ostringstream oss;
  oss << "[";
  if constexpr (N == 1) {
    for (int i = 0; i < arr.size(); ++i) {
      if (i > 0)
        oss << ", ";
      oss << to_string(arr[i]);
    }
  } else {
    oss << "(N-dim type array, size=" << arr.total_size() << ")";
  }
  oss << "]";
  return oss.str();
}

template std::string to_string(const FArray1D<std::complex<double>>&);

// Template overload for FTypeAlloc1D (Allocatable Derived Types)
template <
    typename ViewType,
    auto AllocFunc,
    auto DeallocFunc,
    auto ReallocFunc,
    auto AccessFunc>
std::string to_string(
    const FTypeAlloc1D<
        ViewType,
        AllocFunc,
        DeallocFunc,
        ReallocFunc,
        AccessFunc>& arr) {
  return to_string(arr.view());
}

// Template overload for FAlloc1D (Allocatable Native/Primitive Types)
template <
    typename T,
    auto AllocFunc,
    auto DeallocFunc,
    auto ReallocFunc,
    auto AccessFunc>
std::string to_string(
    const FAlloc1D<T, AllocFunc, DeallocFunc, ReallocFunc, AccessFunc>& arr) {
  return to_string(arr.view());
}
} // namespace Bmad
