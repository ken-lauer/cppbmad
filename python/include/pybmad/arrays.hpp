#pragma once

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"

namespace py = pybind11;
using namespace Bmad;

namespace Pybmad {

inline int normalize_index(int i, int size);

// =============================================================================
// Primitive Array Views (FArray1D<T>)
// =============================================================================
// Binds generic FArray1D types (int, double, complex, etc.)
// Special handling for bool via constexpr.
template <typename T, typename Allocator = void>
void bind_FArray1D(py::module &m, const std::string &name) {
  using Class = FArray1D<T>;
  constexpr bool is_bool = std::is_same_v<T, bool>;

  // 1. Create the class definition
  // Use buffer protocol for standard types, but standard only for bools
  auto cls = [&]() {
    if constexpr (is_bool) {
      return py::class_<Class>(m, name.c_str());
    } else {
      return py::class_<Class>(m, name.c_str(), py::buffer_protocol());
    }
  }();

  cls.def(py::init<>())
      .def("__len__", &Class::size)
      .def("is_valid", &Class::is_valid)
      .def("__str__", [](const Class &self) { return Bmad::to_string(self); })
      .def("__repr__", [](const Class &self) { return Bmad::to_string(self); })
      .def_property_readonly("lower_bound", &Class::lower_bound)
      .def_property_readonly("upper_bound", &Class::upper_bound);

  // 2. Buffer Protocol (Non-Bool only)
  if constexpr (!is_bool) {
    cls.def_buffer([](Class &a) -> py::buffer_info {
      if (!a.is_valid()) {
        throw std::runtime_error("Attempted to access invalid/unallocated Fortran array");
      }
      return py::buffer_info(
          a.data(),
          sizeof(T),
          py::format_descriptor<T>::format(),
          1,
          {(size_t)a.size()},
          {sizeof(T)}
      );
    });
    cls.def("to_list", &Class::to_flat_vector);
  } else {
    cls.def("to_list", &Class::to_vector);
  }

  cls.def(
         "__getitem__",
         [](Class &self, int i) -> T {
           if (i < 0)
             i += self.size();
           return self.at(i);
         }
  ).def("__setitem__", [](Class &self, int i, T val) {
    if (i < 0)
      i += self.size();
    self.at(i) = val; // Works for both reference and Proxy
  });

  // 4. Implicit Conversion from Allocator
  if constexpr (!std::is_void_v<Allocator>) {
    // Define a constructor that takes the Allocator and Views it.
    // keep_alive<1, 2>: The new View (1) ensures Allocator (2) stays alive.
    cls.def(py::init([](Allocator &a) { return a.view(); }), py::keep_alive<1, 2>());

    // Register the implicit conversion handling in pybind11
    py::implicitly_convertible<Allocator, Class>();
  }
}

// =============================================================================
// ND array templates
// =============================================================================

template <typename T, size_t Rank>
void bind_FArrayND(py::module_ &m, const std::string &name) {
  using Class = FArrayND<T, Rank>;

  auto cls =
      py::class_<Class>(m, name.c_str(), py::buffer_protocol())
          .def(py::init<>())
          .def_buffer([](Class &a) -> py::buffer_info {
            if (!a.is_valid()) {
              throw std::runtime_error("Attempted to access invalid FArray");
            }

            std::vector<size_t> shape;
            std::vector<size_t> strides_bytes;

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
                py::format_descriptor<T>::format(), /* Python struct-style format descriptor */
                Rank, /* Number of dimensions */
                shape, /* Buffer dimensions */
                strides_bytes /* Strides (in bytes) */
            );
          })
          .def("is_valid", &Class::is_valid)
          .def("__len__", &Class::total_size)
          .def("__str__", [](const Class &self) { return Bmad::to_string(self); })
          .def("__repr__", [](const Class &self) { return Bmad::to_string(self); })
          .def("to_list", &Class::to_flat_vector)
          .def_property_readonly("total_size", &Class::total_size);

  if constexpr (Rank == 1) {
    // 1D Indexing: arr[i]
    cls.def("__getitem__", [](Class &self, int i) -> T {
      return self.at(normalize_index(i, self.size(1)));
    });
    cls.def("__setitem__", [](Class &self, int i, T val) {
      self.at(normalize_index(i, self.size(1))) = val;
    });
  } else if constexpr (Rank == 2) {
    // 2D Indexing: arr[i, j] -> Python passes a tuple
    cls.def("__getitem__", [](Class &self, py::tuple idx) -> T {
      if (idx.size() != 2)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      return self.at(i, j);
    });
    cls.def("__setitem__", [](Class &self, py::tuple idx, T val) {
      if (idx.size() != 2)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      self.at(i, j) = val;
    });
  } else if constexpr (Rank == 3) {
    // 3D Indexing: arr[i, j, k]
    cls.def("__getitem__", [](Class &self, py::tuple idx) -> T {
      if (idx.size() != 3)
        throw py::index_error("Index tuple size mismatch");
      int i = normalize_index(idx[0].cast<int>(), self.size(1));
      int j = normalize_index(idx[1].cast<int>(), self.size(2));
      int k = normalize_index(idx[2].cast<int>(), self.size(3));
      return self.at(i, j, k);
    });
    cls.def("__setitem__", [](Class &self, py::tuple idx, T val) {
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
// Primitive Allocator Containers (FAlloc1D<T, ...>)
// =============================================================================
// Helper to bind the wrapper containers (RealAlloc1D, etc.)
template <typename AllocClass>
void bind_FAlloc1D(py::module &m, const std::string &name) {
  py::class_<AllocClass>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int>(), py::arg("n"))
      .def("resize", &AllocClass::resize, py::arg("n"))
      .def("resize_bounds", &AllocClass::resize_bounds, py::arg("lbound"), py::arg("ubound"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)
      .def(
          "__getitem__",
          [](AllocClass &self, int i) -> auto {
            if (i < 0)
              i += self.size();
            return self.view().at(i);
          }
      )
      .def(
          "__setitem__",
          [](AllocClass &self, int i, typename AllocClass::view_type::value_type val) {
            if (i < 0)
              i += self.size();
            self.view().at(i) = val;
          }
      )
      .def(
          "__iter__",
          [](AllocClass &self) { return py::make_iterator(self.begin(), self.end()); },
          py::keep_alive<0, 1>()
      )
      .def("view", [](AllocClass &self) { return self.view(); });
}

// BoolAlloc1D specialization; avoid wrapping ReferenceProxy
// Specialization must match the typedef in bmad/fortran_arrays.hpp
// using BoolAlloc1D = FAlloc1D<bool, ...>;
// We cannot verify the exact full type easily, so we use a templated helper if possible
// or just rely on the fact that if we use the EXACT typedef, it works.

template <>
inline void bind_FAlloc1D<BoolAlloc1D>(py::module &m, const std::string &name) {
  using AllocClass = BoolAlloc1D;

  py::class_<AllocClass>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int>(), py::arg("n"))
      .def("resize", &AllocClass::resize, py::arg("n"))
      .def("resize_bounds", &AllocClass::resize_bounds, py::arg("lbound"), py::arg("ubound"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)
      .def("view", [](AllocClass &self) { return self.view(); })

      .def(
          "__getitem__",
          [](AllocClass &self, int i) -> bool {
            if (i < 0)
              i += self.size();
            // .view() is FArray1D<bool>, .at(i) returns ReferenceProxy
            return static_cast<bool>(self.view().at(i));
          }
      )
      .def(
          "__setitem__",
          [](AllocClass &self, int i, bool val) {
            if (i < 0)
              i += self.size();
            self.view().at(i) = val;
          }
      )
      .def(
          "__iter__",
          [](AllocClass &self) {
            // Safe iterator wrapper
            struct BoolAllocIterator {
              int index;
              AllocClass *container;

              using iterator_category = std::forward_iterator_tag;
              using value_type = bool;
              using difference_type = std::ptrdiff_t;
              using pointer = void;
              using reference = bool;

              bool operator==(const BoolAllocIterator &other) const {
                return index == other.index && container == other.container;
              }
              bool operator!=(const BoolAllocIterator &other) const { return !(*this == other); }
              BoolAllocIterator &operator++() {
                ++index;
                return *this;
              }
              bool operator*() const { return static_cast<bool>(container->view().at(index)); }
            };

            return py::make_iterator(
                BoolAllocIterator{0, &self},
                BoolAllocIterator{self.size(), &self}
            );
          },
          py::keep_alive<0, 1>()
      );
}
// =============================================================================
// Derived Type (Proxy) Views (FTypeArrayND<ProxyType, N...>)
// =============================================================================

template <typename ArrayType>
void bind_FTypeArrayND(py::module &m, const std::string &name) {
  py::class_<ArrayType>(m, name.c_str())
      .def(py::init<>())
      .def("__len__", [](const ArrayType &a) { return a.total_size(); })
      .def("is_valid", &ArrayType::is_valid)
      .def(
          "__getitem__",
          [](ArrayType &self, int i) {
            if (i < 0)
              i += static_cast<int>(self.total_size());
            if (i < 0 || i >= static_cast<int>(self.total_size()))
              throw py::index_error();
            return self.at(i);
          }
      )
      .def(
          "__getitem__",
          [](ArrayType &self, py::slice slice) {
            size_t start, stop, step, slice_length;
            if (!slice.compute(self.total_size(), &start, &stop, &step, &slice_length))
              throw py::error_already_set();

            py::list list;
            for (size_t i = 0; i < slice_length; ++i) {
              list.append(self.at(start));
              start += step;
            }
            return list;
          }
      )
      .def(
          "__setitem__",
          [](ArrayType &self, int i, typename ArrayType::value_type &other) {
            if (i < 0)
              i += static_cast<int>(self.total_size());
            if (i < 0 || i >= static_cast<int>(self.total_size()))
              throw py::index_error();

            FortranTraits<typename ArrayType::value_type>::copy(
                other.get_fortran_ptr(),
                self.at(i).get_fortran_ptr()
            );
          }
      )

      .def(
          "__iter__",
          [](ArrayType &self) { return py::make_iterator(self.begin(), self.end()); },
          py::keep_alive<0, 1>()
      )

      .def("__str__", [](const ArrayType &self) { return Bmad::to_string(self); });
}

// =============================================================================
// Derived Type (Proxy) Allocators (FTypeAlloc1D<ViewType...>)
// =============================================================================

template <typename AllocClass>
void bind_FTypeAlloc1D(py::module &m, const std::string &name) {
  py::class_<AllocClass>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int>(), py::arg("n"))
      .def("resize", &AllocClass::resize, py::arg("n"))
      .def("resize_bounds", &AllocClass::resize_bounds, py::arg("lbound"), py::arg("ubound"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)

      .def(
          "__getitem__",
          [](AllocClass &self, int i) {
            if (i < 0)
              i += self.size();
            // .view() refreshes valid pointers, .at() does bounds check
            return self.view().at(i);
          }
      )
      .def(
          "__getitem__",
          [](AllocClass &self, py::slice slice) {
            size_t start, stop, step, slice_length;
            auto &view = self.view();
            if (!slice.compute(view.total_size(), &start, &stop, &step, &slice_length))
              throw py::error_already_set();

            py::list list;
            for (size_t i = 0; i < slice_length; ++i) {
              list.append(view.at(start));
              start += step;
            }
            return list;
          }
      )
      .def(
          "__setitem__",
          [](AllocClass &self, int i, typename AllocClass::view_type::value_type &other) {
            if (i < 0)
              i += static_cast<int>(self.size());

            // .view() refreshes valid pointers, .at() does bounds check
            FortranTraits<typename AllocClass::view_type::value_type>::copy(
                other.get_fortran_ptr(),
                self.view().at(i).get_fortran_ptr()
            );
          }
      )
      .def(
          "__iter__",
          [](AllocClass &self) {
            return py::make_iterator(self.view().begin(), self.view().end());
          },
          py::keep_alive<0, 1>()
      )

      .def("view", [](AllocClass &self) { return self.view(); });

  py::implicitly_convertible<AllocClass, typename AllocClass::view_type>();
}

// =============================================================================
// Char Array Binding (Strings)
// =============================================================================

void bind_FCharArray1D(py::module &m);

void bind_standard_arrays(py::module &m);
}; // namespace Pybmad
