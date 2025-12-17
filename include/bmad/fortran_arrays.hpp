#pragma once

#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits> // std::enable_if / std::is_same
#include <utility>
#include <vector>

// #include "to_string.hpp"

namespace Bmad {

// Forward declarations
template <typename T, std::size_t N>
class FArrayND;
template <typename T>
class FArrayND<T, 1>; // 1D partial specialization

template <typename T>
using FArray1D = FArrayND<T, 1>;
template <typename T>
using FArray2D = FArrayND<T, 2>;
template <typename T>
using FArray3D = FArrayND<T, 3>;

template <
    typename ProxyType,
    std::size_t N,
    void* (*AllocFunc)(int, size_t*) = nullptr,
    void (*DeallocFunc)(void*, int) = nullptr>
class FTypeArrayND;

template <
    typename ProxyType,
    void* (*AllocFunc)(int, size_t*),
    void (*DeallocFunc)(void*, int)>
class FTypeArrayND<
    ProxyType,
    1,
    AllocFunc,
    DeallocFunc>; // 1D partial specialization

template <
    typename ProxyType,
    void* (*AllocFunc)(int, size_t*) = nullptr,
    void (*DeallocFunc)(void*, int) = nullptr>
using FTypeArray1D = FTypeArrayND<ProxyType, 1, AllocFunc, DeallocFunc>;

template <typename ProxyType>
using FTypeArray2D = FTypeArrayND<ProxyType, 2>;
template <typename ProxyType>
using FTypeArray3D = FTypeArrayND<ProxyType, 3>;

class FCharArray1D;

template <
    typename ViewType,
    void* (*AllocFunc)(),
    void (*DeallocFunc)(void*),
    void (*ReallocFunc)(void*, int, size_t),
    void (*AccessFunc)(void*, void**, int*, int*, size_t*, bool*)>
class FTypeAlloc1D;

template <
    typename T,
    void* (*AllocFunc)(),
    void (*DeallocFunc)(void*),
    void (*ReallocFunc)(void*, int, size_t),
    void (*AccessFunc)(void*, void**, int*, int*, size_t*, bool*)>
class FAlloc1D;

using std::to_string;

// =============================================================================
// FArrayND<T, N> - Primary template for N-dimensional primitive arrays
// =============================================================================

template <typename T, std::size_t N>
class FArrayND {
 public:
  // Standard container typedefs
  using value_type = T;

 private:
  T* data_;
  std::array<int, N> sizes_;
  std::array<int, N> lower_bounds_;
  std::array<int, N> upper_bounds_;
  std::array<int, N> strides_;
  bool valid_;

  template <typename... Indices>
  int linear_index_fortran(Indices... indices) const {
    std::array<int, N> idx{static_cast<int>(indices)...};
    int index = 0;
    for (std::size_t i = 0; i < N; ++i)
      index += (idx[i] - lower_bounds_[i]) * strides_[i];
    return index;
  }

  template <typename... Indices>
  int linear_index_c(Indices... indices) const {
    std::array<int, N> idx{static_cast<int>(indices)...};
    int index = 0;
    for (std::size_t i = 0; i < N; ++i)
      index += idx[i] * strides_[i];
    return index;
  }

  template <typename... Indices>
  void check_fortran_bounds(Indices... indices) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    std::array<int, N> idx{static_cast<int>(indices)...};
    for (std::size_t i = 0; i < N; ++i) {
      if (idx[i] < lower_bounds_[i] || idx[i] > upper_bounds_[i])
        throw std::out_of_range("Index out of bounds");
    }
  }

  template <typename... Indices>
  void check_c_bounds(Indices... indices) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    std::array<int, N> idx{static_cast<int>(indices)...};
    for (std::size_t i = 0; i < N; ++i) {
      if (idx[i] < 0 || idx[i] >= sizes_[i])
        throw std::out_of_range("Index out of bounds");
    }
  }

 public:
  FArrayND() : data_(nullptr), valid_(false) {
    sizes_.fill(0);
    lower_bounds_.fill(0);
    upper_bounds_.fill(-1);
    strides_.fill(0);
  }

  FArrayND(
      T* data,
      const std::array<int, N>& sizes,
      const std::array<int, N>& lower_bounds,
      const std::array<int, N>& upper_bounds,
      const std::array<int, N>& strides,
      bool valid)
      : data_(data),
        sizes_(sizes),
        lower_bounds_(lower_bounds),
        upper_bounds_(upper_bounds),
        strides_(strides),
        valid_(valid) {}

  // Variadic constructor: size1, lb1, ub1, size2, lb2, ub2, ..., stride1, stride2, ..., valid
  template <
      typename... Args,
      typename = std::enable_if_t<sizeof...(Args) == 4 * N + 1>>
  FArrayND(T* data, Args... args) : data_(data) {
    std::array<int, 4 * N + 1> a{static_cast<int>(args)...};
    for (std::size_t i = 0; i < N; ++i) {
      sizes_[i] = a[i * 3];
      lower_bounds_[i] = a[i * 3 + 1];
      upper_bounds_[i] = a[i * 3 + 2];
    }
    for (std::size_t i = 0; i < N; ++i)
      strides_[i] = a[3 * N + i];
    valid_ = static_cast<bool>(a[4 * N]);
  }

  // Fortran-style indexing (uses bounds)
  template <typename... Indices>
  T& operator()(Indices... indices) {
    static_assert(sizeof...(indices) == N, "Wrong number of indices");
    check_fortran_bounds(indices...);
    return data_[linear_index_fortran(indices...)];
  }
  template <typename... Indices>
  const T& operator()(Indices... indices) const {
    static_assert(sizeof...(indices) == N, "Wrong number of indices");
    check_fortran_bounds(indices...);
    return data_[linear_index_fortran(indices...)];
  }

  // C-style indexing (0-based)
  template <typename... Indices>
  T& at(Indices... indices) {
    static_assert(sizeof...(indices) == N, "Wrong number of indices");
    check_c_bounds(indices...);
    return data_[linear_index_c(indices...)];
  }
  template <typename... Indices>
  const T& at(Indices... indices) const {
    static_assert(sizeof...(indices) == N, "Wrong number of indices");
    check_c_bounds(indices...);
    return data_[linear_index_c(indices...)];
  }

  template <typename... Indices>
  T& at_fortran(Indices... i) {
    return operator()(i...);
  }
  template <typename... Indices>
  const T& at_fortran(Indices... i) const {
    return operator()(i...);
  }

  bool is_valid() const {
    return valid_;
  }
  const std::array<int, N>& size() const {
    return sizes_;
  }
  int size(int dim) const {
    return sizes_[dim - 1];
  }
  int total_size() const {
    if (!valid_)
      return 0;
    int t = 1;
    for (auto s : sizes_)
      t *= s;
    return t;
  }

  std::array<std::pair<int, int>, N> bounds() const {
    std::array<std::pair<int, int>, N> r;
    for (std::size_t i = 0; i < N; ++i)
      r[i] = {lower_bounds_[i], upper_bounds_[i]};
    return r;
  }
  std::pair<int, int> bounds(int dim) const {
    return {lower_bounds_[dim - 1], upper_bounds_[dim - 1]};
  }
  int lower_bound(int dim) const {
    return lower_bounds_[dim - 1];
  }
  int upper_bound(int dim) const {
    return upper_bounds_[dim - 1];
  }
  const std::array<int, N>& strides() const {
    return strides_;
  }
  constexpr std::size_t rank() const {
    return N;
  }

  T* data() {
    return valid_ ? data_ : nullptr;
  }
  const T* data() const {
    return valid_ ? data_ : nullptr;
  }

  std::vector<T> to_flat_vector() const {
    if (!valid_)
      return {};
    std::vector<T> result;
    result.reserve(total_size());
    for (int i = 0; i < total_size(); ++i)
      result.push_back(data_[i]);
    return result;
  }

  bool empty() const {
    return !valid_ || total_size() == 0;
  }
};

// =============================================================================
// FArrayND<T, 1> - Explicit 1D specialization
// =============================================================================

template <typename T>
class FArrayND<T, 1> {
 public:
  // Standard container typedefs for binding introspection
  using value_type = T;
  using iterator = T*;
  using const_iterator = const T*;
  using reference = T&;
  using const_reference = const T&;

 private:
  T* data_;
  int size_;
  int lower_bound_;
  int upper_bound_;
  bool valid_;

 public:
  // Constructors
  FArrayND()
      : data_(nullptr),
        size_(0),
        lower_bound_(0),
        upper_bound_(-1),
        valid_(false) {}

  FArrayND(T* data, int size, int lower, int upper, bool valid)
      : data_(data),
        size_(size),
        lower_bound_(lower),
        upper_bound_(upper),
        valid_(valid) {}

  // Fortran-style indexing (uses bounds)
  T& operator()(int i) {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < lower_bound_ || i > upper_bound_)
      throw std::out_of_range("Index " + std::to_string(i) + " out of bounds");
    return data_[i - lower_bound_];
  }
  const T& operator()(int i) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < lower_bound_ || i > upper_bound_)
      throw std::out_of_range("Index " + std::to_string(i) + " out of bounds");
    return data_[i - lower_bound_];
  }

  // C-style indexing (0-based)
  T& operator[](int i) {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < 0 || i >= size_)
      throw std::out_of_range("Index out of bounds");
    return data_[i];
  }
  const T& operator[](int i) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < 0 || i >= size_)
      throw std::out_of_range("Index out of bounds");
    return data_[i];
  }

  T& at(int i) {
    return operator[](i);
  }
  const T& at(int i) const {
    return operator[](i);
  }
  T& at_fortran(int i) {
    return operator()(i);
  }
  const T& at_fortran(int i) const {
    return operator()(i);
  }

  // Properties
  bool is_valid() const {
    return valid_;
  }
  int size() const {
    return size_;
  }
  int total_size() const {
    return size_;
  }
  std::pair<int, int> bounds() const {
    return {lower_bound_, upper_bound_};
  }
  int lower_bound() const {
    return lower_bound_;
  }
  int upper_bound() const {
    return upper_bound_;
  }
  constexpr std::size_t rank() const {
    return 1;
  }

  // Data access
  T* data() {
    return valid_ ? data_ : nullptr;
  }
  const T* data() const {
    return valid_ ? data_ : nullptr;
  }

  // Iterators
  T* begin() {
    return valid_ ? data_ : nullptr;
  }
  T* end() {
    return valid_ ? data_ + size_ : nullptr;
  }
  const T* begin() const {
    return valid_ ? data_ : nullptr;
  }
  const T* end() const {
    return valid_ ? data_ + size_ : nullptr;
  }

  // Conversion
  std::vector<T> to_vector() const {
    return valid_ ? std::vector<T>(data_, data_ + size_) : std::vector<T>();
  }
  std::vector<T> to_flat_vector() const {
    return to_vector();
  }
  bool empty() const {
    return !valid_ || size_ == 0;
  }
};

// =============================================================================
// FArrayND<bool, 1> - Specialization for Fortran LOGICAL arrays
// =============================================================================

template <>
class FArrayND<bool, 1> {
 public:
  // Fortran Logical is typically stored as a 4-byte integer (C_INT).
  // We manage the storage as int*, but present the interface as bool.
  using storage_type = int;
  using value_type = bool;

  // Custom proxy to handle conversion between int storage and bool value
  class ReferenceProxy {
    storage_type* ptr_;

   public:
    explicit ReferenceProxy(storage_type* p) : ptr_(p) {}

    // Read: C++ bool from Fortran int (non-zero is true)
    operator bool() const {
      return *ptr_ != 0;
    }

    // Write: Fortran int from C++ bool (1 is true, 0 is false)
    ReferenceProxy& operator=(bool v) {
      *ptr_ = v ? 1 : 0;
      return *this;
    }
    ReferenceProxy& operator=(const ReferenceProxy& other) {
      *ptr_ = other.operator bool() ? 1 : 0;
      return *this;
    }
  };

  class ConstReferenceProxy {
    const storage_type* ptr_;

   public:
    explicit ConstReferenceProxy(const storage_type* p) : ptr_(p) {}
    operator bool() const {
      return *ptr_ != 0;
    }
  };

  // Iterator support
  class iterator {
    storage_type* ptr_;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = bool;
    using pointer = void; // No actual bool* exists
    using reference = ReferenceProxy;
    using iterator_category = std::random_access_iterator_tag;

    iterator(storage_type* p) : ptr_(p) {}
    reference operator*() const {
      return ReferenceProxy(ptr_);
    }
    iterator& operator++() {
      ++ptr_;
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++ptr_;
      return tmp;
    }
    iterator& operator--() {
      --ptr_;
      return *this;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      --ptr_;
      return tmp;
    }
    iterator operator+(int n) const {
      return iterator(ptr_ + n);
    }
    iterator operator-(int n) const {
      return iterator(ptr_ - n);
    }
    difference_type operator-(const iterator& other) const {
      return ptr_ - other.ptr_;
    }
    bool operator==(const iterator& other) const {
      return ptr_ == other.ptr_;
    }
    bool operator!=(const iterator& other) const {
      return ptr_ != other.ptr_;
    }
  };

 private:
  storage_type* data_;
  int size_;
  int lower_bound_;
  int upper_bound_;
  bool valid_;

 public:
  FArrayND()
      : data_(nullptr),
        size_(0),
        lower_bound_(0),
        upper_bound_(-1),
        valid_(false) {}

  // The Constructor receives bool* because FAlloc1D logic casts void* to T* (bool*).
  // We immediately reinterpret_cast it back to the correct storage type (int*).
  FArrayND(bool* data, int size, int lower, int upper, bool valid)
      : data_(reinterpret_cast<storage_type*>(data)),
        size_(size),
        lower_bound_(lower),
        upper_bound_(upper),
        valid_(valid) {}

  // Fortran-style indexing
  ReferenceProxy operator()(int i) {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < lower_bound_ || i > upper_bound_)
      throw std::out_of_range("Index " + std::to_string(i) + " out of bounds");
    return ReferenceProxy(data_ + (i - lower_bound_));
  }

  ConstReferenceProxy operator()(int i) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < lower_bound_ || i > upper_bound_)
      throw std::out_of_range("Index " + std::to_string(i) + " out of bounds");
    return ConstReferenceProxy(data_ + (i - lower_bound_));
  }

  // C-style indexing
  ReferenceProxy operator[](int i) {
    // Omitting validity check for speed in []
    return ReferenceProxy(data_ + i);
  }
  ConstReferenceProxy operator[](int i) const {
    return ConstReferenceProxy(data_ + i);
  }

  ReferenceProxy at(int i) {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < 0 || i >= size_)
      throw std::out_of_range("Index " + std::to_string(i) + " out of bounds");
    return operator[](i);
  }

  ConstReferenceProxy at(int i) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < 0 || i >= size_)
      throw std::out_of_range("Index " + std::to_string(i) + " out of bounds");
    return operator[](i);
  }

  bool is_valid() const {
    return valid_;
  }
  int size() const {
    return size_;
  }
  int total_size() const {
    return size_;
  }
  std::pair<int, int> bounds() const {
    return {lower_bound_, upper_bound_};
  }
  int lower_bound() const {
    return lower_bound_;
  }
  int upper_bound() const {
    return upper_bound_;
  }

  // Important: We cannot return bool* because the memory is int*.
  // We explicitly delete this or return nullptr to avoid dangerous pointer math logic.
  bool* data() = delete;

  // Provide access to the raw storage if strictly necessary
  storage_type* get_storage_ptr() {
    return valid_ ? data_ : nullptr;
  }
  const storage_type* get_storage_ptr() const {
    return valid_ ? data_ : nullptr;
  }

  iterator begin() {
    return iterator(valid_ ? data_ : nullptr);
  }
  iterator end() {
    return iterator(valid_ ? data_ + size_ : nullptr);
  }

  std::vector<bool> to_vector() const {
    if (!valid_)
      return {};
    std::vector<bool> vec;
    vec.reserve(size_);
    for (int i = 0; i < size_; ++i)
      vec.push_back((data_[i] != 0));
    return vec;
  }

  bool empty() const {
    return !valid_ || size_ == 0;
  }
};
// =============================================================================
// FTypeArrayND<ProxyType, N> - Primary template for N-D derived type arrays
// =============================================================================

template <
    typename ProxyType,
    std::size_t N,
    void* (*AllocFunc)(int, size_t*),
    void (*DeallocFunc)(void*, int)>
class FTypeArrayND {
 public:
  using value_type = ProxyType;

 private:
  void* data_;
  std::array<int, N> sizes_;
  std::array<int, N> lower_bounds_;
  std::array<int, N> upper_bounds_;
  std::array<size_t, N> strides_;
  bool valid_;
  size_t element_size_;

  template <typename... Indices>
  size_t linear_index_fortran(Indices... indices) const {
    std::array<int, N> idx{static_cast<int>(indices)...};
    size_t lin = 0;
    for (std::size_t d = 0; d < N; ++d)
      lin += static_cast<size_t>(idx[d] - lower_bounds_[d]) * strides_[d];
    return lin;
  }

  template <typename... Indices>
  size_t linear_index_c(Indices... indices) const {
    std::array<int, N> idx{static_cast<int>(indices)...};
    size_t lin = 0;
    for (std::size_t d = 0; d < N; ++d)
      lin += static_cast<size_t>(idx[d]) * strides_[d];
    return lin;
  }

  void* element_ptr(size_t lin) const {
    return static_cast<char*>(data_) + lin * element_size_;
  }

 public:
  FTypeArrayND() : data_(nullptr), valid_(false), element_size_(0) {
    sizes_.fill(0);
    lower_bounds_.fill(0);
    upper_bounds_.fill(-1);
    strides_.fill(0);
  }

  FTypeArrayND(
      void* data,
      const std::array<int, N>& sizes,
      const std::array<int, N>& lower_bounds,
      const std::array<int, N>& upper_bounds,
      const std::array<size_t, N>& strides,
      bool valid,
      size_t element_size)
      : data_(data),
        sizes_(sizes),
        lower_bounds_(lower_bounds),
        upper_bounds_(upper_bounds),
        strides_(strides),
        valid_(valid),
        element_size_(element_size) {}

  template <typename... Indices>
  ProxyType operator()(Indices... indices) {
    static_assert(sizeof...(indices) == N, "Wrong number of indices");
    return ProxyType(element_ptr(linear_index_fortran(indices...)));
  }
  template <typename... Indices>
  const ProxyType operator()(Indices... indices) const {
    static_assert(sizeof...(indices) == N, "Wrong number of indices");
    return ProxyType(element_ptr(linear_index_fortran(indices...)));
  }

  template <typename... Indices>
  ProxyType at(Indices... indices) {
    return ProxyType(element_ptr(linear_index_c(indices...)));
  }
  template <typename... Indices>
  const ProxyType at(Indices... indices) const {
    return ProxyType(element_ptr(linear_index_c(indices...)));
  }

  bool is_valid() const {
    return valid_;
  }
  size_t element_size() const {
    return element_size_;
  }
  size_t total_size() const {
    if (!valid_)
      return 0;
    size_t t = 1;
    for (auto s : sizes_)
      t *= s;
    return t;
  }
  const std::array<int, N>& sizes() const {
    return sizes_;
  }
  int size(int dim) const {
    return sizes_[dim - 1];
  }
  std::array<std::pair<int, int>, N> bounds() const {
    std::array<std::pair<int, int>, N> b;
    for (std::size_t i = 0; i < N; ++i)
      b[i] = {lower_bounds_[i], upper_bounds_[i]};
    return b;
  }
  std::pair<int, int> bounds(int dim) const {
    return {lower_bounds_[dim - 1], upper_bounds_[dim - 1]};
  }
  int lower_bound(int dim) const {
    return lower_bounds_[dim - 1];
  }
  int upper_bound(int dim) const {
    return upper_bounds_[dim - 1];
  }
  const std::array<size_t, N>& strides() const {
    return strides_;
  }

  void* data() {
    return valid_ ? data_ : nullptr;
  }
  const void* data() const {
    return valid_ ? data_ : nullptr;
  }
  void* get_fortran_ptr() const {
    return data_;
  }
  bool empty() const {
    return !valid_ || total_size() == 0;
  }

  std::vector<ProxyType> to_vector() const {
    if (!valid_)
      return {};
    std::vector<ProxyType> r;
    r.reserve(total_size());
    for (size_t i = 0; i < total_size(); ++i)
      r.emplace_back(element_ptr(i));
    return r;
  }
  std::vector<ProxyType> to_flat_vector() const {
    return to_vector();
  }

  class iterator {
    const FTypeArrayND* arr_;
    size_t idx_;

   public:
    iterator(const FTypeArrayND* a, size_t i) : arr_(a), idx_(i) {}
    ProxyType operator*() {
      return ProxyType(arr_->element_ptr(idx_));
    }
    iterator& operator++() {
      ++idx_;
      return *this;
    }
    bool operator==(const iterator& other) const {
      return idx_ == other.idx_ && arr_ == other.arr_;
    }

    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }
  };
  iterator begin() {
    return iterator(this, 0);
  }
  iterator end() {
    return iterator(this, total_size());
  }
};

// =============================================================================
// FTypeArrayND<ProxyType, 1> - Explicit 1D specialization with ownership
// =============================================================================

template <
    typename ProxyType,
    void* (*AllocFunc)(int, size_t*),
    void (*DeallocFunc)(void*, int)>
class FTypeArrayND<ProxyType, 1, AllocFunc, DeallocFunc> {
 public:
  using value_type = ProxyType;

 private:
  std::shared_ptr<void> data_;
  int size_;
  int lower_bound_;
  int upper_bound_;
  bool valid_;
  size_t element_size_;

  void* element_ptr(int i) const {
    return static_cast<char*>(data_.get()) + i * element_size_;
  }

 public:
  // Constructors
  FTypeArrayND()
      : data_(nullptr),
        size_(0),
        lower_bound_(0),
        upper_bound_(-1),
        valid_(false),
        element_size_(0) {}

  FTypeArrayND(
      void* data,
      int size,
      int lower,
      int upper,
      bool valid,
      size_t elem_size)
      : size_(size),
        lower_bound_(lower),
        upper_bound_(upper),
        valid_(valid),
        element_size_(elem_size) {
    if (data) {
      // Non-owning view: use no-op deleter
      data_ = std::shared_ptr<void>(data, [](void*) {});
    }
  }

  // Static factory for owned arrays
  // NOTE: depends on 'U' to defer evaluation until the function is actually called.
  template <typename U = void>
  static std::enable_if_t<
      std::is_same<U, void>::value && (AllocFunc != nullptr) &&
          (DeallocFunc != nullptr),
      FTypeArrayND>
  allocate(int size, int lower_bound = 1) {
    if (size < 0)
      throw std::invalid_argument("Size must be non-negative");

    size_t elem_size = 0;
    void* raw_ptr = AllocFunc(size, &elem_size);
    if (!raw_ptr)
      throw std::runtime_error("Failed to allocate");

    FTypeArrayND arr;
    arr.size_ = size;
    // Helper lambda to bridge C++ shared_ptr deleter to Fortran DeallocFunc
    arr.data_ = std::shared_ptr<void>(
        raw_ptr, [size](void* p) { DeallocFunc(p, size); });
    arr.lower_bound_ = lower_bound;
    arr.upper_bound_ = lower_bound + size - 1;
    arr.valid_ = true;
    arr.element_size_ = elem_size;
    return arr;
  }

  // Destructor: Rule of Zero (handled by shared_ptr)

  // Copy: Rule of Zero (shared ownership enabled via shared_ptr)

  // Move: Rule of Zero (handled by shared_ptr)

  FTypeArrayND as_view() const {
    return *this; // shared_ptr enables implicit view semantics
  }

  // is_owned() is conceptually removed as shared_ptr handles lifetime abstractly

  // Fortran-style indexing
  ProxyType operator()(int i) {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < lower_bound_ || i > upper_bound_)
      throw std::out_of_range("Index out of bounds");
    return ProxyType(element_ptr(i - lower_bound_));
  }
  const ProxyType operator()(int i) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < lower_bound_ || i > upper_bound_)
      throw std::out_of_range("Index out of bounds");
    return ProxyType(element_ptr(i - lower_bound_));
  }

  // C-style indexing
  ProxyType operator[](int i) {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < 0 || i >= size_)
      throw std::out_of_range("Index out of bounds");
    return ProxyType(element_ptr(i));
  }
  const ProxyType operator[](int i) const {
    if (!valid_)
      throw std::runtime_error("Array not allocated");
    if (i < 0 || i >= size_)
      throw std::out_of_range("Index out of bounds");
    return ProxyType(element_ptr(i));
  }

  ProxyType at(int i) {
    return operator[](i);
  }
  const ProxyType at(int i) const {
    return operator[](i);
  }
  ProxyType at_fortran(int i) {
    return operator()(i);
  }
  const ProxyType at_fortran(int i) const {
    return operator()(i);
  }

  // Properties
  bool is_valid() const {
    return valid_;
  }
  int size() const {
    return size_;
  }
  size_t total_size() const {
    return static_cast<size_t>(size_);
  }
  size_t element_size() const {
    return element_size_;
  }
  std::pair<int, int> bounds() const {
    return {lower_bound_, upper_bound_};
  }
  int lower_bound() const {
    return lower_bound_;
  }
  int upper_bound() const {
    return upper_bound_;
  }

  void* data() {
    return valid_ ? data_.get() : nullptr;
  }
  const void* data() const {
    return valid_ ? data_.get() : nullptr;
  }
  void* get_fortran_ptr() const {
    return data_.get();
  }
  bool empty() const {
    return !valid_ || size_ == 0;
  }

  std::vector<ProxyType> to_vector() const {
    if (!valid_)
      return {};
    std::vector<ProxyType> r;
    r.reserve(size_);
    for (int i = 0; i < size_; ++i)
      r.emplace_back(element_ptr(i));
    return r;
  }
  std::vector<ProxyType> to_flat_vector() const {
    return to_vector();
  }

  // Iterators
  class iterator {
    const FTypeArrayND* arr_;
    int idx_;

   public:
    iterator(const FTypeArrayND* a, int i) : arr_(a), idx_(i) {}
    ProxyType operator*() {
      return ProxyType(arr_->element_ptr(idx_));
    }
    iterator& operator++() {
      ++idx_;
      return *this;
    }
    bool operator==(const iterator& o) const {
      return idx_ == o.idx_ && arr_ == o.arr_;
    }
    bool operator!=(const iterator& o) const {
      return !(*this == o);
    }
  };
  iterator begin() {
    return iterator(this, 0);
  }
  iterator end() {
    return iterator(this, size_);
  }
  iterator begin() const {
    return iterator(this, 0);
  }
  iterator end() const {
    return iterator(this, size_);
  }
};

// =============================================================================
// FCharArray1D - Character string arrays (special handling)
// =============================================================================

class FCharArray1D {
 private:
  char* data_;
  int size_;
  int lower_bound_;
  int upper_bound_;
  int str_len_;
  bool valid_;

  std::string trim(const char* str, int len) const {
    int end = len - 1;
    while (end >= 0 && (str[end] == ' ' || str[end] == '\0'))
      --end;
    return std::string(str, end + 1);
  }

  void pad(char* dest, const std::string& src, int len) const {
    int n = std::min(static_cast<int>(src.length()), len);
    std::memcpy(dest, src.c_str(), n);
    for (int i = n; i < len; ++i)
      dest[i] = ' ';
  }

 public:
  FCharArray1D()
      : data_(nullptr),
        size_(0),
        lower_bound_(0),
        upper_bound_(-1),
        str_len_(0),
        valid_(false) {}

  FCharArray1D(
      char* data,
      int size,
      int lower,
      int upper,
      int str_len,
      bool valid)
      : data_(data),
        size_(size),
        lower_bound_(lower),
        upper_bound_(upper),
        str_len_(str_len),
        valid_(valid) {}

  class StringProxy {
    char* data_;
    int len_;
    FCharArray1D* parent_;

   public:
    StringProxy(char* d, int l, FCharArray1D* p)
        : data_(d), len_(l), parent_(p) {}
    operator std::string() const {
      return parent_->trim(data_, len_);
    }
    StringProxy& operator=(const std::string& s) {
      parent_->pad(data_, s, len_);
      return *this;
    }
    std::string str() const {
      return parent_->trim(data_, len_);
    }
    char* data() {
      return data_;
    }
    int length() const {
      return len_;
    }
  };

  class ConstStringProxy {
    const char* data_;
    int len_;
    const FCharArray1D* parent_;

   public:
    ConstStringProxy(const char* d, int l, const FCharArray1D* p)
        : data_(d), len_(l), parent_(p) {}
    operator std::string() const {
      return parent_->trim(data_, len_);
    }
    std::string str() const {
      return parent_->trim(data_, len_);
    }
    const char* data() const {
      return data_;
    }
    int length() const {
      return len_;
    }
  };

  StringProxy operator()(int i) {
    return StringProxy(data_ + (i - lower_bound_) * str_len_, str_len_, this);
  }
  ConstStringProxy operator()(int i) const {
    return ConstStringProxy(
        data_ + (i - lower_bound_) * str_len_, str_len_, this);
  }
  StringProxy operator[](int i) {
    return StringProxy(data_ + i * str_len_, str_len_, this);
  }
  ConstStringProxy operator[](int i) const {
    return ConstStringProxy(data_ + i * str_len_, str_len_, this);
  }

  std::string get_string(int i) const {
    return trim(data_ + (i - lower_bound_) * str_len_, str_len_);
  }
  void set_string(int i, const std::string& s) {
    pad(data_ + (i - lower_bound_) * str_len_, s, str_len_);
  }

  bool is_valid() const {
    return valid_;
  }
  int size() const {
    return size_;
  }
  int string_length() const {
    return str_len_;
  }
  std::pair<int, int> bounds() const {
    return {lower_bound_, upper_bound_};
  }
  int lower_bound() const {
    return lower_bound_;
  }
  int upper_bound() const {
    return upper_bound_;
  }
  char* data() {
    return valid_ ? data_ : nullptr;
  }
  const char* data() const {
    return valid_ ? data_ : nullptr;
  }
  bool empty() const {
    return !valid_ || size_ == 0;
  }

  std::vector<std::string> to_vector() const {
    if (!valid_)
      return {};
    std::vector<std::string> r;
    r.reserve(size_);
    for (int i = 0; i < size_; ++i)
      r.push_back(trim(data_ + i * str_len_, str_len_));
    return r;
  }

  class Iterator {
    FCharArray1D* p_;
    int i_;

    friend bool operator==(const Iterator& a, const Iterator& b) {
      return a.p_ == b.p_ && a.i_ == b.i_;
    }
    friend bool operator!=(const Iterator& a, const Iterator& b) {
      return !(a == b);
    }

   public:
    Iterator(FCharArray1D* p, int i) : p_(p), i_(i) {}
    StringProxy operator*() {
      return (*p_)[i_];
    }
    Iterator& operator++() {
      ++i_;
      return *this;
    }
    bool operator!=(const Iterator& o) const {
      return i_ != o.i_;
    }
  };
  Iterator begin() {
    return Iterator(this, 0);
  }
  Iterator end() {
    return Iterator(this, size_);
  }
};

// =============================================================================
// FTypeAlloc1D - Container for Fortran allocatable arrays
// =============================================================================

template <
    typename ViewType,
    void* (*AllocFunc)(),
    void (*DeallocFunc)(void*),
    void (*ReallocFunc)(void*, int, size_t),
    void (*AccessFunc)(void*, void**, int*, int*, size_t*, bool*)>
class FTypeAlloc1D {
 public:
  using view_type = ViewType;

 private:
  std::shared_ptr<void> handle_;
  mutable ViewType view_;
  mutable bool stale_ = true;

  void refresh() const {
    void* data = nullptr;
    int lbound = 1, size = 0;
    size_t elem_size = 0;
    bool alloc = false;
    AccessFunc(handle_.get(), &data, &lbound, &size, &elem_size, &alloc);
    view_ = (alloc && data && size > 0)
        ? ViewType(data, size, lbound, lbound + size - 1, true, elem_size)
        : ViewType();
    stale_ = false;
  }

 public:
  FTypeAlloc1D() {
    handle_ = std::shared_ptr<void>(AllocFunc(), DeallocFunc);
  }
  explicit FTypeAlloc1D(int lbound, int n) : FTypeAlloc1D() {
    if (n > 0)
      resize(lbound, n);
  }
  // Destructor: Rule of Zero (handled by shared_ptr + DeallocFunc)

  void resize(int lbound, int n) {
    ReallocFunc(handle_.get(), lbound, n);
    stale_ = true;
  }
  void clear() {
    ReallocFunc(handle_.get(), 0, 0);
    stale_ = true;
  }

  void* get_fortran_ptr() const {
    return handle_.get();
  }
  ViewType& view() const {
    refresh();
    return view_;
  }
  int size() const {
    return view().size();
  }
  bool empty() const {
    return size() == 0;
  }

  auto begin() {
    return view().begin();
  }
  auto end() {
    return view().end();
  }
  decltype(auto) operator[](int i) {
    return view()[i];
  }
  decltype(auto) operator[](int i) const {
    return view()[i];
  }
};

// =============================================================================
// FAlloc1D - Container for Fortran allocatable arrays (native types)
// =============================================================================

template <
    typename T,
    void* (*AllocFunc)(),
    void (*DeallocFunc)(void*),
    void (*ReallocFunc)(void*, int, size_t),
    void (*AccessFunc)(void*, void**, int*, int*, size_t*, bool*)>
class FAlloc1D {
 public:
  using view_type = FArray1D<T>;

 private:
  std::shared_ptr<void> handle_;
  mutable FArray1D<T> view_;
  mutable bool stale_ = true;

  void refresh() const {
    void* data_ptr = nullptr;
    int lbound = 1, size = 0;
    size_t elem_size = 0;
    bool alloc = false;
    AccessFunc(handle_.get(), &data_ptr, &lbound, &size, &elem_size, &alloc);
    T* data = static_cast<T*>(data_ptr);
    view_ = (alloc && data && size > 0)
        ? FArray1D<T>(data, size, lbound, lbound + size - 1, true)
        : FArray1D<T>();
    stale_ = false;
  }

 public:
  FAlloc1D() {
    handle_ = std::shared_ptr<void>(AllocFunc(), DeallocFunc);
  }
  explicit FAlloc1D(int lbound, int n) : FAlloc1D() {
    if (n > 0)
      resize(lbound, n);
  }
  // Destructor: Rule of Zero (handled by shared_ptr + DeallocFunc)

  void resize(int lbound, int n) {
    ReallocFunc(handle_.get(), lbound, n);
    stale_ = true;
  }
  void clear() {
    ReallocFunc(handle_.get(), 0, 0);
    stale_ = true;
  }

  void* get_fortran_ptr() const {
    return handle_.get();
  }
  FArray1D<T>& view() const {
    refresh();
    return view_;
  }
  FArray1D<T>& get() const {
    return view();
  }
  int size() const {
    return view().size();
  }
  bool empty() const {
    return size() == 0;
  }

  auto begin() {
    return view().begin();
  }
  auto end() {
    return view().end();
  }
  decltype(auto) operator[](int i) {
    return view()[i];
  }
  decltype(auto) operator[](int i) const {
    return view()[i];
  }
  decltype(auto) operator()(int i) {
    return view()(i);
  }
  decltype(auto) operator()(int i) const {
    return view()(i);
  }
};

// =============================================================================
// ProxyHelpers - Helper functions for proxy generation
// =============================================================================

class ProxyHelpers {
 public:
  template <typename T, typename Func>
  static FArray1D<T> get_array_1d(const void* ptr, Func f) {
    T* data = nullptr;
    int bounds[2];
    bool alloc = false;
    f(ptr, &data, bounds, &alloc);
    int size = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    return FArray1D<T>(data, size, bounds[0], bounds[1], alloc);
  }

  template <typename T, typename Func>
  static FArray2D<T> get_array_2d(const void* ptr, Func f) {
    T* data = nullptr;
    int bounds[4], strides[2];
    bool alloc = false;
    f(ptr, &data, bounds, strides, &alloc);
    int s1 = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    int s2 = alloc ? (bounds[3] - bounds[2] + 1) : 0;
    return FArray2D<T>(
        data,
        s1,
        bounds[0],
        bounds[1],
        s2,
        bounds[2],
        bounds[3],
        strides[0],
        strides[1],
        alloc);
  }

  template <typename T, typename Func>
  static FArray3D<T> get_array_3d(const void* ptr, Func f) {
    T* data = nullptr;
    int bounds[6], strides[3];
    bool alloc = false;
    f(ptr, &data, bounds, strides, &alloc);
    int s1 = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    int s2 = alloc ? (bounds[3] - bounds[2] + 1) : 0;
    int s3 = alloc ? (bounds[5] - bounds[4] + 1) : 0;
    return FArray3D<T>(
        data,
        s1,
        bounds[0],
        bounds[1],
        s2,
        bounds[2],
        bounds[3],
        s3,
        bounds[4],
        bounds[5],
        strides[0],
        strides[1],
        strides[2],
        alloc);
  }

  template <typename Func>
  static FCharArray1D get_char_array_1d(const void* ptr, Func f) {
    char* data = nullptr;
    int bounds[2], str_len = 0;
    bool alloc = false;
    f(ptr, &data, bounds, &str_len, &alloc);
    int size = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    return FCharArray1D(data, size, bounds[0], bounds[1], str_len, alloc);
  }

  template <typename Func>
  static std::string get_string(const void* ptr, Func f) {
    char* data = nullptr;
    int len = 0;
    bool alloc = false;
    f(ptr, &data, &len, &alloc);
    return (alloc && data && len > 0) ? std::string(data, len) : std::string();
  }

  template <typename ArrayT, typename Func>
  static ArrayT get_type_array_1d(const void* ptr, Func f) {
    void* data = nullptr;
    int bounds[2] = {0, 0};
    bool alloc = false;
    size_t elem_size = 0;
    f(ptr, &data, bounds, &alloc, &elem_size);
    int size = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    return ArrayT(data, size, bounds[0], bounds[1], alloc, elem_size);
  }

  template <typename ArrayT, typename Func>
  static ArrayT get_type_array_2d(const void* ptr, Func f) {
    void* data = nullptr;
    int bounds[4], strides_in[2];
    bool alloc = false;
    size_t elem_size = 0;
    f(ptr, &data, bounds, strides_in, &alloc, &elem_size);
    int s1 = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    int s2 = alloc ? (bounds[3] - bounds[2] + 1) : 0;
    std::array<int, 2> sizes = {s1, s2};
    std::array<int, 2> lowers = {bounds[0], bounds[2]};
    std::array<int, 2> uppers = {bounds[1], bounds[3]};
    std::array<size_t, 2> strides = {
        static_cast<size_t>(strides_in[0]), static_cast<size_t>(strides_in[1])};
    return ArrayT(data, sizes, lowers, uppers, strides, alloc, elem_size);
  }

  template <typename ArrayT, typename Func>
  static ArrayT get_type_array_3d(const void* ptr, Func f) {
    void* data = nullptr;
    int bounds[6], strides_in[3];
    bool alloc = false;
    size_t elem_size = 0;
    f(ptr, &data, bounds, strides_in, &alloc, &elem_size);
    int s1 = alloc ? (bounds[1] - bounds[0] + 1) : 0;
    int s2 = alloc ? (bounds[3] - bounds[2] + 1) : 0;
    int s3 = alloc ? (bounds[5] - bounds[4] + 1) : 0;
    std::array<int, 3> sizes = {s1, s2, s3};
    std::array<int, 3> lowers = {bounds[0], bounds[2], bounds[4]};
    std::array<int, 3> uppers = {bounds[1], bounds[3], bounds[5]};
    std::array<size_t, 3> strides = {
        static_cast<size_t>(strides_in[0]),
        static_cast<size_t>(strides_in[1]),
        static_cast<size_t>(strides_in[2])};
    return ArrayT(data, sizes, lowers, uppers, strides, alloc, elem_size);
  }
};

} // namespace Bmad
