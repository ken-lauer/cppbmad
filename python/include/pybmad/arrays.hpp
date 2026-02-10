#pragma once

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"

namespace py = pybind11;
using namespace Bmad;

namespace Pybmad {

inline int normalize_index(int i, int size);

// =============================================================================
// Unified Binder: View (FArray1D) + Allocator (FAlloc1D) + std::vector
// =============================================================================

template <typename T, typename AllocClass>
void bind_1D_array_pair(
    py::module &m,
    const std::string &view_name,
    const std::string &alloc_name
) {
  using ViewClass = FArray1D<T>;
  constexpr bool is_bool = std::is_same_v<T, bool>;

  auto view_cls = [&]() {
    if constexpr (is_bool) {
      return py::class_<ViewClass>(m, view_name.c_str());
    } else {
      return py::class_<ViewClass>(m, view_name.c_str(), py::buffer_protocol());
    }
  }();

  view_cls.def(py::init<>())
      .def("__len__", &ViewClass::size)
      .def("is_valid", &ViewClass::is_valid)
      .def("__str__", [](const ViewClass &self) { return Bmad::to_string(self); })
      .def("__repr__", [](const ViewClass &self) { return Bmad::to_string(self); })
      .def_property_readonly("lower_bound", &ViewClass::lower_bound)
      .def_property_readonly("upper_bound", &ViewClass::upper_bound);

  // View Get/Set Items
  view_cls
      .def(
          "__getitem__",
          [](ViewClass &self, int i) -> T {
            if (i < 0)
              i += self.size();
            return self.at(i);
          }
      )
      .def("__setitem__", [](ViewClass &self, int i, T val) {
        if (i < 0)
          i += self.size();
        self.at(i) = val;
      });

  view_cls.def("to_list", &ViewClass::to_vector);

  if constexpr (!is_bool) {
    view_cls.def_buffer([](ViewClass &a) -> py::buffer_info {
      if (!a.is_valid())
        throw std::runtime_error("Invalid Fortran array access");
      return py::buffer_info(
          a.data(),
          sizeof(T),
          py::format_descriptor<T>::format(),
          1,
          {(size_t)a.size()},
          {sizeof(T)}
      );
    });
  }

  // --------------------------------------------------------
  // 2. Definition of Allocator Class (RealAlloc1D, etc.)
  // --------------------------------------------------------
  auto alloc_cls =
      py::class_<AllocClass>(m, alloc_name.c_str())
          .def(py::init<>())
          .def(py::init<int>(), py::arg("n"))
          .def("resize", &AllocClass::resize, py::arg("n"))
          // resize_bounds missing in some subtypes? If universal, keep it.
          .def("resize_bounds", &AllocClass::resize_bounds, py::arg("lbound"), py::arg("ubound"))
          .def("clear", &AllocClass::clear)
          .def("__len__", &AllocClass::size)
          .def("view", [](AllocClass &self) { return self.view(); }, py::keep_alive<0, 1>());

  if constexpr (!is_bool) {
    alloc_cls.def(
        "__iter__",
        [](AllocClass &self) { return py::make_iterator(self.begin(), self.end()); },
        py::keep_alive<0, 1>()
    );
  } else {
    alloc_cls.def("__iter__", [](AllocClass &self) {
      py::list result;
      auto &view = self.view();
      for (int i = 0; i < view.size(); ++i) {
        result.append(static_cast<bool>(view.at(i)));
      }
      return result.attr("__iter__")();
    });
  }

  // Allocator Proxy Get/Set (Passes through to View)
  // We use T for return type to force evaluation of Proxy->Value
  alloc_cls
      .def(
          "__getitem__",
          [](AllocClass &self, int i) -> T {
            if (i < 0)
              i += self.size();
            return self.view().at(i);
          }
      )
      .def("__setitem__", [](AllocClass &self, int i, T val) {
        if (i < 0)
          i += self.size();
        self.view().at(i) = val;
      });

  // --------------------------------------------------------
  // 3. Implicit Conversations & Interop
  // --------------------------------------------------------

  // A. View constructed from Allocator (View <--- Alloc)
  //    Essential for functions returning Alloc to Python as Views,
  //    or implicit conversion in args.
  view_cls.def(py::init([](AllocClass &a) { return a.view(); }), py::keep_alive<1, 2>());
  py::implicitly_convertible<AllocClass, ViewClass>();

  // B. Allocator constructed from View (Alloc <--- View) (Deep Copy)
  alloc_cls.def(py::init([](const ViewClass &v) {
    if (!v.is_valid())
      throw std::runtime_error("Cannot copy from invalid view");
    AllocClass a;
    a.resize(v.size());
    for (int i = 0; i < v.size(); i++)
      a[i] = v[i];
    return a;
  }));
  py::implicitly_convertible<ViewClass, AllocClass>();

  // C. Allocator constructed from std::vector/List (Alloc <--- List) (Deep Copy)
  //    This logic merges the std::vector equivalent into the Allocator system.
  alloc_cls.def(py::init([](const std::vector<T> &vec) {
    AllocClass a;
    a.resize(vec.size());
    auto &v = a.view();
    for (size_t i = 0; i < vec.size(); i++)
      v[i] = vec[i];
    return a;
  }));
  py::implicitly_convertible<std::vector<T>, AllocClass>();
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
}
// =============================================================================
// Derived Type (Proxy) Alloc + View (FTypeArrayND<ProxyType, N...>)
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
          },
          py::keep_alive<0, 1>()
      )
      .def(
          "__getitem__",
          [](py::object self_py, py::slice slice) {
            auto &self = self_py.cast<ArrayType &>();
            size_t start, stop, step, slice_length;
            if (!slice.compute(self.total_size(), &start, &stop, &step, &slice_length))
              throw py::error_already_set();

            py::list list;
            for (size_t i = 0; i < slice_length; ++i) {
              auto item = self.at(start);

              py::object py_item =
                  py::cast(item, py::return_value_policy::reference_internal, self_py);

              list.append(py_item);
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

template <typename ArrayType, typename AllocClass>
void bind_1d_type_array_pair(
    py::module &m,
    const std::string &view_name,
    const std::string &alloc_name
) {
  using ProxyType = typename ArrayType::value_type;

  // --------------------------------------------------------
  // 1. Definition of View Class (FTypeArrayND)
  // --------------------------------------------------------
  py::class_<ArrayType> view_cls(m, view_name.c_str());

  view_cls.def(py::init<>())
      .def("__len__", [](const ArrayType &a) { return a.total_size(); })
      .def("is_valid", &ArrayType::is_valid)
      .def("__str__", [](const ArrayType &self) { return Bmad::to_string(self); });

  // View: Single Item Access
  view_cls
      .def(
          "__getitem__",
          [](ArrayType &self, int i) {
            if (i < 0)
              i += static_cast<int>(self.total_size());
            if (i < 0 || i >= static_cast<int>(self.total_size()))
              throw py::index_error();
            return self.at(i);
          },
          py::keep_alive<0, 1>()
      )
      .def("__setitem__", [](ArrayType &self, int i, ProxyType &other) {
        if (i < 0)
          i += static_cast<int>(self.total_size());
        if (i < 0 || i >= static_cast<int>(self.total_size()))
          throw py::index_error();

        // Deep copy via FortranTraits
        // The Proxy assignment operator usually only copies the pointer wrapper,
        // not the underlying data contents.
        FortranTraits<ProxyType>::copy(other.get_fortran_ptr(), self.at(i).get_fortran_ptr());
      });

  // View: Slice Access
  view_cls.def(
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
      },
      py::keep_alive<0, 1>()
  );

  // View: Iteration
  view_cls.def(
      "__iter__",
      [](ArrayType &self) { return py::make_iterator(self.begin(), self.end()); },
      py::keep_alive<0, 1>()
  );

  // View: Export to standard vector
  view_cls.def("to_list", &ArrayType::to_vector, py::keep_alive<0, 1>());

  // --------------------------------------------------------
  // 2. Definition of Allocator Class (FTypeAlloc1D)
  // --------------------------------------------------------
  py::class_<AllocClass> alloc_cls(m, alloc_name.c_str());

  alloc_cls.def(py::init<>())
      .def(py::init<int>(), py::arg("n"))
      .def("resize", &AllocClass::resize, py::arg("n"))
      .def("resize_bounds", &AllocClass::resize_bounds, py::arg("lbound"), py::arg("ubound"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)
      .def("view", [](AllocClass &self) { return self.view(); }, py::keep_alive<0, 1>());

  // Allocator: Single Item Access (Delegates to View logic)
  alloc_cls
      .def(
          "__getitem__",
          [](AllocClass &self, int i) {
            if (i < 0)
              i += self.size();
            // .view() refreshes pointers if underlying realloc happened
            return self.view().at(i);
          },
          py::keep_alive<0, 1>()
      )
      .def("__setitem__", [](AllocClass &self, int i, ProxyType &other) {
        auto &v = self.view();
        int size = static_cast<int>(v.total_size());

        if (i < 0)
          i += size;
        if (i < 0 || i >= size)
          throw py::index_error();

        FortranTraits<ProxyType>::copy(other.get_fortran_ptr(), v.at(i).get_fortran_ptr());
      });

  // Allocator: Slice Access
  alloc_cls.def("__getitem__", [](py::object self_py, py::slice slice) {
    auto &self = self_py.cast<AllocClass &>();

    auto &view = self.view();
    size_t start, stop, step, slice_length;
    if (!slice.compute(view.total_size(), &start, &stop, &step, &slice_length))
      throw py::error_already_set();

    py::list list;
    for (size_t i = 0; i < slice_length; ++i) {
      auto item = view.at(start);
      py::object py_item = py::cast(item, py::return_value_policy::reference_internal, self_py);
      list.append(py_item);
      start += step;
    }
    return list;
  });

  // Allocator: Iteration
  alloc_cls.def(
      "__iter__",
      [](AllocClass &self) { return py::make_iterator(self.view().begin(), self.view().end()); },
      py::keep_alive<0, 1>()
  );

  // --------------------------------------------------------
  // 3. Implicit Conversions & Interop
  // --------------------------------------------------------

  // A. View constructed from Allocator (View <--- Alloc)
  // Allows functions returning AllocClass references to be caught by py::view_class bindings
  view_cls.def(py::init([](AllocClass &a) { return a.view(); }), py::keep_alive<1, 2>());
  py::implicitly_convertible<AllocClass, ArrayType>();

  // B. Allocator constructed from View (Alloc <--- View) (Deep Copy)
  alloc_cls.def(py::init([](const ArrayType &v) {
    if (!v.is_valid())
      throw std::runtime_error("Cannot copy from invalid view");

    AllocClass a;
    a.resize(static_cast<int>(v.total_size()));

    auto &target_view = a.view();

    for (size_t i = 0; i < v.total_size(); i++) {
      // Explicit Deep Copy for derived types
      FortranTraits<ProxyType>::copy(
          v.at(i).get_fortran_ptr(),
          target_view.at(i).get_fortran_ptr()
      );
    }
    return a;
  }));
  py::implicitly_convertible<ArrayType, AllocClass>();

  // C. Allocator constructed from std::vector/List (Alloc <--- List) (Deep Copy)
  // Allows passing a Python list of proxies to a function expecting AllocClass
  alloc_cls.def(py::init([](const std::vector<ProxyType> &vec) {
    AllocClass a;
    a.resize(static_cast<int>(vec.size()));

    auto &target_view = a.view();

    for (size_t i = 0; i < vec.size(); i++) {
      // vec[i] is a Proxy wrapper pointing to source data,
      // target_view[i] is wrapper pointing to alloc's new memory.
      FortranTraits<ProxyType>::copy(vec[i].get_fortran_ptr(), target_view.at(i).get_fortran_ptr());
    }
    return a;
  }));
  py::implicitly_convertible<std::vector<ProxyType>, AllocClass>();
}

// =============================================================================
// Char Array Binding (Strings)
// =============================================================================

void bind_FCharArray1D(py::module &m);

void bind_standard_arrays(py::module &m);
}; // namespace Pybmad
