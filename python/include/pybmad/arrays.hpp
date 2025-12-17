#pragma once

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/generated/to_string.hpp"

namespace py = pybind11;
using namespace Bmad;

namespace Pybmad {

inline int normalize_index(int i, int size);

// =============================================================================
// 1. Primitive Array Views (FArray1D<T>)
// =============================================================================
// Binds generic FArray1D types (int, double, complex, etc.)
template <typename T>
void bind_FArray1D(py::module& m, const std::string& name) {
  using Class = FArray1D<T>;

  py::class_<Class>(m, name.c_str(), py::buffer_protocol())
      .def(py::init<>())
      // Buffer protocol allows zero-copy numpy access!
      .def_buffer([](Class& a) -> py::buffer_info {
        if (!a.is_valid()) {
          throw std::runtime_error(
              "Attempted to access invalid/unallocated Fortran array");
        }
        return py::buffer_info(
            a.data(), /* Pointer to buffer */
            sizeof(T), /* Size of one scalar */
            py::format_descriptor<
                T>::format(), /* Python struct-style format descriptor */
            1, /* Number of dimensions */
            {(size_t)a.size()}, /* Buffer dimensions */
            {sizeof(T)} /* Strides (in bytes) */
        );
      })
      .def("__len__", &Class::size)
      // C-style indexing (0-based) for Python consistency
      .def(
          "__getitem__",
          [](Class& self, int i) -> T {
            // handle negative indexing commonly used in Python
            if (i < 0)
              i += self.size();
            // relying on C++ at() checks, but capturing standard exception if needed
            return self.at(i);
          })
      .def(
          "__setitem__",
          [](Class& self, int i, T val) {
            if (i < 0)
              i += self.size();
            self.at(i) = val;
          })
      .def("__str__", [](const Class& self) { return Bmad::to_string(self); })
      .def("__repr__", [](const Class& self) { return Bmad::to_string(self); })
      .def("is_valid", &Class::is_valid)
      .def("to_list", &Class::to_flat_vector)
      .def_property_readonly("lower_bound", &Class::lower_bound)
      .def_property_readonly("upper_bound", &Class::upper_bound);
}

// =============================================================================
// Specialization: bind_FArray1D<bool>
// =============================================================================

template <>
inline void bind_FArray1D<bool>(py::module& m, const std::string& name) {
  using Class = FArray1D<bool>;

  // Note: py::buffer_protocol() is not supported here.
  py::class_<Class>(m, name.c_str())
      .def(py::init<>())
      .def("__len__", &Class::size)
      .def(
          "__getitem__",
          [](Class& self, int i) -> bool {
            if (i < 0)
              i += self.size();
            return self.at(i);
          })
      .def(
          "__setitem__",
          [](Class& self, int i, bool val) {
            if (i < 0)
              i += self.size();
            self.at(i) = val;
          })
      .def("__str__", [](const Class& self) { return Bmad::to_string(self); })
      .def("__repr__", [](const Class& self) { return Bmad::to_string(self); })
      .def("is_valid", &Class::is_valid)
      .def("to_list", &Class::to_vector)
      .def_property_readonly("lower_bound", &Class::lower_bound)
      .def_property_readonly("upper_bound", &Class::upper_bound);
}

// =============================================================================
// ND array templates
// =============================================================================

template <typename T, size_t Rank>
void bind_FArrayND(py::module_& m, const std::string& name) {
  using Class = FArrayND<T, Rank>;

  auto cls =
      py::class_<Class>(m, name.c_str(), py::buffer_protocol())
          .def(py::init<>())
          // ---------------------------------------------------------------------
          // Buffer Protocol: Exposes memory directly to NumPy
          // ---------------------------------------------------------------------
          .def_buffer([](Class& a) -> py::buffer_info {
            if (!a.is_valid()) {
              throw std::runtime_error("Attempted to access invalid FArray");
            }

            // Prepare shape and strides vectors
            std::vector<size_t> shape;
            std::vector<size_t> strides_bytes;

            // Access raw internal arrays for dimensions
            // Note: Bmad sizes/strides are std::array internally accessible via public getters
            auto bmad_sizes = a.size();
            auto bmad_strides = a.strides();

            for (size_t i = 0; i < Rank; ++i) {
              shape.push_back(bmad_sizes[i]);
              // Bmad stores element strides; NumPy wants byte strides
              strides_bytes.push_back(bmad_strides[i] * sizeof(T));
            }

            return py::buffer_info(
                a.data(), /* Pointer to buffer */
                sizeof(T), /* Size of one scalar */
                py::format_descriptor<
                    T>::format(), /* Python struct-style format descriptor */
                Rank, /* Number of dimensions */
                shape, /* Buffer dimensions */
                strides_bytes /* Strides (in bytes) */
            );
          })
          // ---------------------------------------------------------------------
          // Standard Methods
          // ---------------------------------------------------------------------
          .def("is_valid", &Class::is_valid)
          .def("__len__", &Class::total_size)
          .def(
              "__str__",
              [](const Class& self) { return Bmad::to_string(self); })
          .def(
              "__repr__",
              [](const Class& self) { return Bmad::to_string(self); })
          .def("to_list", &Class::to_flat_vector)
          .def_property_readonly("total_size", &Class::total_size);

  // -------------------------------------------------------------------------
  // Dimensionality Specific Indexing (SFINAE / if constexpr usually,
  // but explicit specialization is cleaner for pybind lambdas)
  // -------------------------------------------------------------------------

  if constexpr (Rank == 1) {
    // 1D Indexing: arr[i]
    cls.def("__getitem__", [](Class& self, int i) -> T {
      return self.at(normalize_index(i, self.size(1)));
    });
    cls.def("__setitem__", [](Class& self, int i, T val) {
      self.at(normalize_index(i, self.size(1))) = val;
    });
  } else if constexpr (Rank == 2) {
    // 2D Indexing: arr[i, j] -> Python passes a tuple
    cls.def("__getitem__", [](Class& self, py::tuple idx) -> T {
      if (idx.size() != 2)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      return self.at(i, j);
    });
    cls.def("__setitem__", [](Class& self, py::tuple idx, T val) {
      if (idx.size() != 2)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      self.at(i, j) = val;
    });
  } else if constexpr (Rank == 3) {
    // 3D Indexing: arr[i, j, k]
    cls.def("__getitem__", [](Class& self, py::tuple idx) -> T {
      if (idx.size() != 3)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      int k = normalize_index(idx[2].cast<int>(), self.size(3));
      return self.at(i, j, k);
    });
    cls.def("__setitem__", [](Class& self, py::tuple idx, T val) {
      if (idx.size() != 3)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      int k = normalize_index(idx[2].cast<int>(), self.size(3));
      self.at(i, j, k) = val;
    });
  }
}
// =============================================================================
// 2. Primitive Allocator Containers (FAlloc1D<T, ...>)
// =============================================================================
// Helper to bind the wrapper containers (RealAlloc1D, etc.)
template <typename AllocClass>
void bind_FAlloc1D(py::module& m, const std::string& name) {
  py::class_<AllocClass>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int, int>(), py::arg("lbound"), py::arg("n"))
      .def("resize", &AllocClass::resize, py::arg("lbound"), py::arg("n"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)
      .def(
          "__getitem__",
          [](AllocClass& self, int i) -> auto {
            // Forward to the inner view's checked access
            // Requires calling view() ensures refresh
            if (i < 0)
              i += self.size();
            return self.view().at(i);
          })
      .def(
          "__setitem__",
          [](AllocClass& self,
             int i,
             typename AllocClass::view_type::value_type val) {
            if (i < 0)
              i += self.size();
            self.view().at(i) = val;
          })
      // Allow treating it like an iterable
      .def(
          "__iter__",
          [](AllocClass& self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      // Expose the view object if they want the numpy-compatible object
      .def("view", [](AllocClass& self) { return self.view(); });
}

// BoolAlloc1D specialization; avoid wrapping ReferenceProxy

template <>
inline void bind_FAlloc1D<BoolAlloc1D>(py::module& m, const std::string& name) {
  // 1. Alias the concrete type for cleaner code
  using AllocClass = BoolAlloc1D;

  py::class_<AllocClass>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int, int>(), py::arg("lbound"), py::arg("n"))
      // Standard container methods
      .def("resize", &AllocClass::resize, py::arg("lbound"), py::arg("n"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)
      // Expose the inner view (returns FArrayND<bool,1>, which must be bound separately!)
      .def("view", [](AllocClass& self) { return self.view(); })

      .def(
          "__getitem__",
          [](AllocClass& self, int i) -> bool {
            if (i < 0)
              i += self.size();
            return self.view().at(i);
          })

      // 3. Safe Item Assignment
      .def(
          "__setitem__",
          [](AllocClass& self, int i, bool val) {
            if (i < 0)
              i += self.size();
            self.view().at(i) = val;
          })

      // 4. Custom Iterator (Yields bools, bypasses ReferenceProxy)
      .def(
          "__iter__",
          [](AllocClass& self) {
            // Helper struct to iterate by index and return simple bools
            struct BoolIndexIterator {
              int index;
              AllocClass*
                  container; // Keep unsafe ptr, reliant on keep_alive below

              using iterator_category = std::forward_iterator_tag;
              using value_type = bool;
              using difference_type = std::ptrdiff_t;
              using pointer = void;
              using reference = bool;

              bool operator==(const BoolIndexIterator& other) const {
                return index == other.index;
              }
              bool operator!=(const BoolIndexIterator& other) const {
                return index != other.index;
              }
              BoolIndexIterator& operator++() {
                ++index;
                return *this;
              }
              // DEREFERENCE: Use .at() to get Proxy, convert to bool
              bool operator*() const {
                return container->view().at(index);
              }
            };

            return py::make_iterator(
                BoolIndexIterator{0, &self},
                BoolIndexIterator{self.size(), &self});
          },
          // IMPORTANT: Keep the container alive while the iterator exists
          py::keep_alive<0, 1>());
}
// =============================================================================
// 3. Derived Type (Proxy) Views (FTypeArrayND<ProxyType, N...>)
// =============================================================================

// Helper for N-Dimensional Proxy Arrays (Works for 1, 2, 3...)
// We pass the full View template instantiation as ArrayType
// e.g. ArrayType = SplineProxyArray1D
template <typename ArrayType>
void bind_FTypeArrayND(py::module& m, const std::string& name) {
  py::class_<ArrayType>(m, name.c_str())
      .def(
          py::init<>()) // Usually these are created C++ side and returned, but init useful for tests
      .def("__len__", [](const ArrayType& a) { return a.total_size(); })
      .def("is_valid", &ArrayType::is_valid)

      // __getitem__ variadic depending on N is hard in generic lambda.
      // We will default to flat C-style indexing or 1D indexing for simplicity here,
      // unless you want to specialize for N=2, N=3 specifically for python syntax arr[x,y].
      // Pybind generally maps [x,y] to a tuple argument for __getitem__.

      .def(
          "__getitem__",
          [](ArrayType& self, int i) {
            // Accessing the proxy creates a temporary object attached to that memory
            // which Pybind will cast to Python.
            // Note: using linear_index_c logic via 'at'
            if (i < 0)
              i += static_cast<int>(self.total_size());
            return self.at(i);
          })
      // We don't usually Set the proxy itself, we get the proxy and set its fields

      .def(
          "__iter__",
          [](ArrayType& self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())

      .def("__str__", [](const ArrayType& self) {
        return Bmad::to_string(self);
      });
}

// Special check: FTypeArrayND<..., 1, ...> has "allocate" static method and ownership
// We can use SFINAE or just overload names if we know specific types.
// For the 1D proxies with allocators (like SplineProxyArray1D), we might want access to .resize()
// essentially provided by the FTypeArrayND 1D logic if it was owned,
// but usually in Bmad, FTypeAlloc1D manages the memory.
// We stick to the View binding above for pure arrays.

// =============================================================================
// 4. Derived Type (Proxy) Allocators (FTypeAlloc1D<ViewType...>)
// =============================================================================

template <typename AllocClass>
void bind_FTypeAlloc1D(py::module& m, const std::string& name) {
  py::class_<AllocClass>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int, int>(), py::arg("lbound"), py::arg("n"))
      .def("resize", &AllocClass::resize, py::arg("lbound"), py::arg("n"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)

      // Item access returns the Proxy object
      .def(
          "__getitem__",
          [](AllocClass& self, int i) {
            if (i < 0)
              i += self.size();
            // .view() refreshes valid pointers, .at() does bounds check
            return self.view().at(i);
          })
      .def(
          "__iter__",
          [](AllocClass& self) {
            return py::make_iterator(self.view().begin(), self.view().end());
          },
          py::keep_alive<0, 1>())

      .def("view", [](AllocClass& self) { return self.view(); });
}

// =============================================================================
// 5. Char Array Binding (Strings)
// =============================================================================

void bind_FCharArray1D(py::module& m);

void bind_standard_arrays(py::module& m);
}; // namespace Pybmad
