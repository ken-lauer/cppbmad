#pragma once

#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstring>
#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"

namespace nb = nanobind;
using namespace Bmad;

namespace Pybmad {

inline int normalize_index(int i, int size);

// numpy 2.x calls __array__(dtype=None, copy=None). Returns true only when a
// guaranteed copy was explicitly requested (e.g. np.array(x), copy=True). For
// copy=False/None we may hand back a zero-copy view into the Fortran memory.
inline bool array_copy_requested(const nb::kwargs &kwargs) {
  if (kwargs.contains("copy")) {
    nb::object c = kwargs["copy"];
    if (!c.is_none())
      return nb::cast<bool>(c);
  }
  return false;
}

// Build a numpy array that owns its own copy of the data, so mutations do not
// propagate back into the Fortran-managed buffer. ArrT is the (ranked) ndarray
// return type, so the copy and zero-copy paths share one type.
template <typename ArrT, typename T>
ArrT make_numpy_copy(const T *src, size_t ndim, const size_t *shape, char order = 'C') {
  size_t total = 1;
  for (size_t i = 0; i < ndim; ++i)
    total *= shape[i];
  T *buf = new T[total];
  std::memcpy(buf, src, total * sizeof(T));
  nb::capsule owner(buf, [](void *p) noexcept { delete[] static_cast<T *>(p); });
  return ArrT(buf, ndim, shape, owner, nullptr, nb::dtype<T>(), nb::device::cpu::value, 0, order);
}

// =============================================================================
// Unified Binder: View (FArray1D) + Allocator (FAlloc1D) + std::vector
// =============================================================================

template <typename T, typename AllocClass>
void bind_1D_array_pair(
    nb::module_ &m,
    const std::string &view_name,
    const std::string &alloc_name
) {
  using ViewClass = FArray1D<T>;
  constexpr bool is_bool = std::is_same_v<T, bool>;

  auto view_cls = nb::class_<ViewClass>(m, view_name.c_str());

  view_cls.def(nb::init<>())
      .def("__len__", &ViewClass::size)
      .def("is_valid", &ViewClass::is_valid)
      .def(
          "__str__",
          [](const ViewClass &self) {
            return self.is_valid() ? Bmad::to_string(self) : "<invalid>";
          }
      )
      .def(
          "__repr__",
          [](const ViewClass &self) {
            return self.is_valid() ? Bmad::to_string(self) : "<invalid>";
          }
      )
      .def_prop_ro("lower_bound", &ViewClass::lower_bound)
      .def_prop_ro("upper_bound", &ViewClass::upper_bound);

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
      .def(
          "__setitem__",
          [](ViewClass &self, int i, T val) {
            if (i < 0)
              i += self.size();
            self.at(i) = val;
          }
      )
      .def("__setitem__", [](ViewClass &self, nb::slice slice, nb::sequence val) {
        auto [start, stop, step, slicelength] = slice.compute(self.size());
        if (nb::len(val) != slicelength)
          throw nb::value_error(
              "Left and right hand size of slice assignment have different sizes"
          );
        for (size_t i = 0; i < slicelength; ++i) {
          self.at(start + i * step) = nb::cast<T>(val[i]);
        }
      });

  // Slice read: zero-copy strided numpy view for numeric dtypes (matches numpy
  // slicing semantics). Fortran logical is 4 bytes/elem, so bool falls back to
  // a list copy.
  if constexpr (!is_bool) {
    view_cls.def("__getitem__", [](nb::handle_t<ViewClass> self_py, nb::slice slice) {
      auto &self = nb::cast<ViewClass &>(self_py);
      auto [start, stop, step, slicelength] = slice.compute(self.size());
      size_t shape[1] = {(size_t)slicelength};
      int64_t strides[1] = {(int64_t)step};
      return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
          self.data() + start,
          1,
          shape,
          nb::handle(self_py),
          strides
      );
    });
  } else {
    view_cls.def("__getitem__", [](ViewClass &self, nb::slice slice) {
      auto [start, stop, step, slicelength] = slice.compute(self.size());
      nb::list ret;
      for (size_t i = 0; i < slicelength; ++i)
        ret.append(static_cast<bool>(self.at(start + i * step)));
      return ret;
    });
  }

  view_cls.def("to_list", &ViewClass::to_vector);

  if constexpr (!is_bool) {
    view_cls.def(
        "__iter__",
        [](ViewClass &self) {
          return nb::make_iterator(nb::type<ViewClass>(), "ViewIterator", self.begin(), self.end());
        },
        nb::keep_alive<0, 1>()
    );
  } else {
    view_cls.def("__iter__", [](ViewClass &self) {
      nb::list result;
      for (int i = 0; i < self.size(); ++i) {
        result.append(static_cast<bool>(self.at(i)));
      }
      return result.attr("__iter__")();
    });
  }

  if constexpr (!is_bool) {
    view_cls.def("__array__", [](nb::handle_t<ViewClass> self, nb::kwargs kwargs) {
      using Ret = nb::ndarray<nb::numpy, T, nb::ndim<1>>;
      auto &a = nb::cast<ViewClass &>(self);
      if (!a.is_valid())
        throw std::runtime_error("Invalid Fortran array access");
      size_t shape[1] = {(size_t)a.size()};
      if (array_copy_requested(kwargs))
        return make_numpy_copy<Ret>(a.data(), 1, shape);
      return Ret(a.data(), 1, shape, nb::handle(self));
    });
  }

  // --------------------------------------------------------
  // 2. Definition of Allocator Class (RealAlloc1D, etc.)
  // --------------------------------------------------------
  auto alloc_cls =
      nb::class_<AllocClass>(m, alloc_name.c_str())
          .def(nb::init<>())
          .def(nb::init<int>(), nb::arg("n"))
          .def("resize", &AllocClass::resize, nb::arg("n"))
          // resize_bounds missing in some subtypes? If universal, keep it.
          .def("resize_bounds", &AllocClass::resize_bounds, nb::arg("lbound"), nb::arg("ubound"))
          .def("clear", &AllocClass::clear)
          .def("__len__", &AllocClass::size)
          .def("view", [](AllocClass &self) { return self.view(); }, nb::keep_alive<0, 1>());

  if constexpr (!is_bool) {
    alloc_cls.def(
        "__iter__",
        [](AllocClass &self) {
          return nb::make_iterator(
              nb::type<AllocClass>(),
              "AllocIterator",
              self.begin(),
              self.end()
          );
        },
        nb::keep_alive<0, 1>()
    );

    alloc_cls.def("__array__", [](nb::handle_t<AllocClass> self, nb::kwargs kwargs) {
      using Ret = nb::ndarray<nb::numpy, T, nb::ndim<1>>;
      auto &a = nb::cast<AllocClass &>(self);
      auto &v = a.view();
      if (!v.is_valid())
        throw std::runtime_error("Invalid Fortran array access");
      size_t shape[1] = {(size_t)v.size()};
      if (array_copy_requested(kwargs))
        return make_numpy_copy<Ret>(v.data(), 1, shape);
      return Ret(v.data(), 1, shape, nb::handle(self));
    });
  } else {
    alloc_cls.def("__iter__", [](AllocClass &self) {
      nb::list result;
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
            auto &v = self.view();
            if (i < 0)
              i += v.size();
            return v.at(i);
          }
      )
      .def(
          "__setitem__",
          [](AllocClass &self, int i, T val) {
            auto &v = self.view();
            if (i < 0)
              i += v.size();
            v.at(i) = val;
          }
      )
      .def("__setitem__", [](AllocClass &self, nb::slice slice, nb::sequence val) {
        auto &v = self.view();
        auto [start, stop, step, slicelength] = slice.compute(v.size());
        if (nb::len(val) != slicelength)
          throw nb::value_error(
              "Left and right hand size of slice assignment have different sizes"
          );
        for (size_t i = 0; i < slicelength; ++i) {
          v.at(start + i * step) = nb::cast<T>(val[i]);
        }
      });

  // Slice read: zero-copy strided numpy view for numeric dtypes, list for bool.
  if constexpr (!is_bool) {
    alloc_cls.def("__getitem__", [](nb::handle_t<AllocClass> self_py, nb::slice slice) {
      auto &self = nb::cast<AllocClass &>(self_py);
      auto &v = self.view();
      auto [start, stop, step, slicelength] = slice.compute(v.size());
      size_t shape[1] = {(size_t)slicelength};
      int64_t strides[1] = {(int64_t)step};
      return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
          v.data() + start,
          1,
          shape,
          nb::handle(self_py),
          strides
      );
    });
  } else {
    alloc_cls.def("__getitem__", [](AllocClass &self, nb::slice slice) {
      auto &v = self.view();
      auto [start, stop, step, slicelength] = slice.compute(v.size());
      nb::list ret;
      for (size_t i = 0; i < slicelength; ++i)
        ret.append(static_cast<bool>(v.at(start + i * step)));
      return ret;
    });
  }

  // --------------------------------------------------------
  // 3. Implicit Conversations & Interop
  // --------------------------------------------------------

  // A. View constructed from Allocator (View <--- Alloc)
  //    Essential for functions returning Alloc to Python as Views,
  //    or implicit conversion in args.
  view_cls.def(
      "__init__",
      [](ViewClass *self, AllocClass &a) { new (self) ViewClass(a.view()); },
      nb::keep_alive<1, 2>()
  );
  nb::implicitly_convertible<AllocClass, ViewClass>();

  // B. Allocator constructed from View (Alloc <--- View) (Deep Copy)
  alloc_cls.def("__init__", [](AllocClass *self, const ViewClass &v) {
    if (!v.is_valid())
      throw std::runtime_error("Cannot copy from invalid view");
    new (self) AllocClass();
    self->resize(v.size());
    for (int i = 0; i < v.size(); i++)
      (*self)[i] = v[i];
  });
  nb::implicitly_convertible<ViewClass, AllocClass>();

  // C. Allocator constructed from std::vector/List (Alloc <--- List) (Deep Copy)
  //    This logic merges the std::vector equivalent into the Allocator system.
  alloc_cls.def("__init__", [](AllocClass *self, const std::vector<T> &vec) {
    new (self) AllocClass();
    self->resize(vec.size());
    auto &v = self->view();
    for (size_t i = 0; i < vec.size(); i++)
      v[i] = vec[i];
  });
  nb::implicitly_convertible<std::vector<T>, AllocClass>();
}

// =============================================================================
// ND array templates
// =============================================================================

template <typename T, size_t Rank>
void bind_FArrayND(nb::module_ &m, const std::string &name) {
  using Class = FArrayND<T, Rank>;

  auto cls = nb::class_<Class>(m, name.c_str())
                 .def(nb::init<>())
                 .def("is_valid", &Class::is_valid)
                 .def("__len__", &Class::total_size)
                 .def(
                     "__str__",
                     [](const Class &self) {
                       return self.is_valid() ? Bmad::to_string(self) : "<invalid>";
                     }
                 )
                 .def(
                     "__repr__",
                     [](const Class &self) {
                       return self.is_valid() ? Bmad::to_string(self) : "<invalid>";
                     }
                 )
                 .def("to_list", &Class::to_flat_vector)
                 .def_prop_ro("total_size", &Class::total_size);

  // Fortran logical uses 4 bytes per element; numpy bool uses 1 byte.
  // Zero-copy __array__ is only safe for types with matching element sizes.
  if constexpr (!std::is_same_v<T, bool>) {
    cls.def("__array__", [](nb::handle_t<Class> self, nb::kwargs kwargs) {
      // Fortran arrays are column-major and contiguous (non-contiguous arrays
      // are reported as invalid by the access layer), hence ndim<Rank>/f_contig.
      using Ret = nb::ndarray<nb::numpy, T, nb::ndim<Rank>, nb::f_contig>;
      auto &a = nb::cast<Class &>(self);
      if (!a.is_valid()) {
        throw std::runtime_error("Attempted to access invalid FArray");
      }

      size_t shape[Rank];
      auto bmad_sizes = a.size();
      for (size_t i = 0; i < Rank; ++i) {
        shape[i] = bmad_sizes[i];
      }

      if (array_copy_requested(kwargs))
        return make_numpy_copy<Ret>(a.data(), Rank, shape, 'F');
      return Ret(
          a.data(),
          Rank,
          shape,
          nb::handle(self),
          /* strides */ nullptr,
          nanobind::dtype<T>(),
          nanobind::device::cpu::value,
          0,
          'F'
      );
    });
  }
}
// =============================================================================
// Derived Type (Proxy) Alloc + View (FTypeArrayND<ProxyType, N...>)
// =============================================================================

template <typename ArrayType>
void bind_FTypeArrayND(nb::module_ &m, const std::string &name) {
  nb::class_<ArrayType>(m, name.c_str())
      .def(nb::init<>())
      .def("__len__", [](const ArrayType &a) { return a.total_size(); })
      .def("is_valid", &ArrayType::is_valid)
      .def(
          "__getitem__",
          [](ArrayType &self, int i) {
            if (i < 0)
              i += static_cast<int>(self.total_size());
            if (i < 0 || i >= static_cast<int>(self.total_size()))
              throw nb::index_error();
            return self.at(i);
          },
          nb::keep_alive<0, 1>()
      )
      .def(
          "__getitem__",
          [](nb::object self_py, nb::slice slice) {
            auto &self = nb::cast<ArrayType &>(self_py);
            auto [start, stop, step, slice_length] = slice.compute(self.total_size());

            nb::list list;
            for (size_t i = 0; i < slice_length; ++i) {
              auto item = self.at(start);

              nb::object py_item =
                  nb::cast(item, nb::rv_policy::reference_internal, nb::handle(self_py));

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
              throw nb::index_error();

            FortranTraits<typename ArrayType::value_type>::copy(
                other.get_fortran_ptr(),
                self.at(i).get_fortran_ptr()
            );
          }
      )

      .def(
          "__iter__",
          [](ArrayType &self) {
            return nb::make_iterator(
                nb::type<ArrayType>(),
                "FTypeArrayNDIterator",
                self.begin(),
                self.end()
            );
          },
          nb::keep_alive<0, 1>()
      )

      .def("__str__", [](const ArrayType &self) {
        return self.is_valid() ? Bmad::to_string(self) : "<invalid>";
      });
}

template <typename ArrayType, typename AllocClass>
void bind_1d_type_array_pair(
    nb::module_ &m,
    const std::string &view_name,
    const std::string &alloc_name
) {
  using ProxyType = typename ArrayType::value_type;

  // --------------------------------------------------------
  // 1. Definition of View Class (FTypeArrayND)
  // --------------------------------------------------------
  nb::class_<ArrayType> view_cls(m, view_name.c_str());

  view_cls.def(nb::init<>())
      .def("__len__", [](const ArrayType &a) { return a.total_size(); })
      .def("is_valid", &ArrayType::is_valid)
      .def("__str__", [](const ArrayType &self) {
        return self.is_valid() ? Bmad::to_string(self) : "<invalid>";
      });

  // View: Single Item Access
  view_cls
      .def(
          "__getitem__",
          [](ArrayType &self, int i) {
            if (i < 0)
              i += static_cast<int>(self.total_size());
            if (i < 0 || i >= static_cast<int>(self.total_size()))
              throw nb::index_error();
            return self.at(i);
          },
          nb::keep_alive<0, 1>()
      )
      .def("__setitem__", [](ArrayType &self, int i, ProxyType &other) {
        if (i < 0)
          i += static_cast<int>(self.total_size());
        if (i < 0 || i >= static_cast<int>(self.total_size()))
          throw nb::index_error();

        // Deep copy via FortranTraits
        // The Proxy assignment operator usually only copies the pointer wrapper,
        // not the underlying data contents.
        FortranTraits<ProxyType>::copy(other.get_fortran_ptr(), self.at(i).get_fortran_ptr());
      });

  // View: Slice Access
  view_cls.def(
      "__getitem__",
      [](ArrayType &self, nb::slice slice) {
        auto [start, stop, step, slice_length] = slice.compute(self.total_size());

        nb::list list;
        for (size_t i = 0; i < slice_length; ++i) {
          list.append(self.at(start));
          start += step;
        }
        return list;
      },
      nb::keep_alive<0, 1>()
  );

  // View: Iteration
  view_cls.def(
      "__iter__",
      [](ArrayType &self) {
        return nb::make_iterator(
            nb::type<ArrayType>(),
            "TypeArray1DIterator",
            self.begin(),
            self.end()
        );
      },
      nb::keep_alive<0, 1>()
  );

  // View: Export to standard vector
  view_cls.def("to_list", &ArrayType::to_vector, nb::keep_alive<0, 1>());

  // --------------------------------------------------------
  // 2. Definition of Allocator Class (FTypeAlloc1D)
  // --------------------------------------------------------
  nb::class_<AllocClass> alloc_cls(m, alloc_name.c_str());

  alloc_cls.def(nb::init<>())
      .def(nb::init<int>(), nb::arg("n"))
      .def("resize", &AllocClass::resize, nb::arg("n"))
      .def("resize_bounds", &AllocClass::resize_bounds, nb::arg("lbound"), nb::arg("ubound"))
      .def("clear", &AllocClass::clear)
      .def("__len__", &AllocClass::size)
      .def("view", [](AllocClass &self) { return self.view(); }, nb::keep_alive<0, 1>());

  // Allocator: Single Item Access (Delegates to View logic)
  alloc_cls
      .def(
          "__getitem__",
          [](AllocClass &self, int i) {
            auto &v = self.view();
            if (i < 0)
              i += v.size();
            return v.at(i);
          },
          nb::keep_alive<0, 1>()
      )
      .def("__setitem__", [](AllocClass &self, int i, ProxyType &other) {
        auto &v = self.view();
        int size = static_cast<int>(v.total_size());

        if (i < 0)
          i += size;
        if (i < 0 || i >= size)
          throw nb::index_error();

        FortranTraits<ProxyType>::copy(other.get_fortran_ptr(), v.at(i).get_fortran_ptr());
      });

  // Allocator: Slice Access
  alloc_cls.def("__getitem__", [](nb::object self_py, nb::slice slice) {
    auto &self = nb::cast<AllocClass &>(self_py);

    auto &view = self.view();
    auto [start, stop, step, slice_length] = slice.compute(view.total_size());

    nb::list list;
    for (size_t i = 0; i < slice_length; ++i) {
      nb::object py_item = nb::cast(view.at(start), nb::rv_policy::move);
      // tie the item back to AllocClass's lifetime
      nb::detail::keep_alive(py_item.ptr(), self_py.ptr());
      list.append(py_item);
      start += step;
    }
    return list;
  });

  // Allocator: Iteration
  alloc_cls.def(
      "__iter__",
      [](AllocClass &self) {
        auto &v = self.view();
        return nb::make_iterator(nb::type<AllocClass>(), "TypeAlloc1DIterator", v.begin(), v.end());
      },
      nb::keep_alive<0, 1>()
  );

  // --------------------------------------------------------
  // 3. Implicit Conversions & Interop
  // --------------------------------------------------------

  // A. View constructed from Allocator (View <--- Alloc)
  // Allows functions returning AllocClass references to be caught by view_class bindings
  view_cls.def(
      "__init__",
      [](ArrayType *self, AllocClass &a) { new (self) ArrayType(a.view()); },
      nb::keep_alive<1, 2>()
  );
  nb::implicitly_convertible<AllocClass, ArrayType>();

  // B. Allocator constructed from View (Alloc <--- View) (Deep Copy)
  alloc_cls.def("__init__", [](AllocClass *self, const ArrayType &v) {
    if (!v.is_valid())
      throw std::runtime_error("Cannot copy from invalid view");

    new (self) AllocClass();
    self->resize(static_cast<int>(v.total_size()));

    auto &target_view = self->view();

    for (size_t i = 0; i < v.total_size(); i++) {
      // Explicit Deep Copy for derived types
      FortranTraits<ProxyType>::copy(
          v.at(i).get_fortran_ptr(),
          target_view.at(i).get_fortran_ptr()
      );
    }
  });
  nb::implicitly_convertible<ArrayType, AllocClass>();

  // C. Allocator constructed from std::vector/List (Alloc <--- List) (Deep Copy)
  // Allows passing a Python list of proxies to a function expecting AllocClass
  alloc_cls.def("__init__", [](AllocClass *self, const std::vector<ProxyType> &vec) {
    new (self) AllocClass();
    self->resize(static_cast<int>(vec.size()));

    auto &target_view = self->view();

    for (size_t i = 0; i < vec.size(); i++) {
      // vec[i] is a Proxy wrapper pointing to source data,
      // target_view[i] is wrapper pointing to alloc's new memory.
      FortranTraits<ProxyType>::copy(vec[i].get_fortran_ptr(), target_view.at(i).get_fortran_ptr());
    }
  });
  nb::implicitly_convertible<std::vector<ProxyType>, AllocClass>();
}

// =============================================================================
// Char Array Binding (Strings)
// =============================================================================

void bind_FCharArray1D(nb::module_ &m);
void bind_FCharAlloc1D(nb::module_ &m);

void bind_standard_arrays(nb::module_ &m);
}; // namespace Pybmad
