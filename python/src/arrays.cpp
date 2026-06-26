#include "pybmad/arrays.hpp"

namespace Pybmad {

inline int normalize_index(int i, int size) {
  if (i < 0)
    i += size;
  if (i < 0 || i >= size)
    throw nb::index_error();
  return i;
}
void bind_FCharAlloc1D(nb::module_ &m) {
  nb::class_<CharacterAlloc1D>(m, "CharacterAlloc1D")
      .def(nb::init<>())
      .def(nb::init<int, int>(), nb::arg("n"), nb::arg("str_len") = 200)
      .def("resize", &CharacterAlloc1D::resize, nb::arg("n"), nb::arg("str_len"))
      .def("clear", &CharacterAlloc1D::clear)
      .def("__len__", &CharacterAlloc1D::size)
      .def("string_length", &CharacterAlloc1D::string_length)
      .def("to_list", &CharacterAlloc1D::to_vector)
      .def(
          "__getitem__",
          [](CharacterAlloc1D &self, int i) -> std::string {
            auto &v = self.view();
            int n = v.size();
            if (i < 0)
              i += n;
            if (i < 0 || i >= n)
              throw nb::index_error();
            return v[i];
          }
      )
      .def(
          "__setitem__",
          [](CharacterAlloc1D &self, int i, const std::string &val) {
            auto &v = self.view();
            int n = v.size();
            if (i < 0)
              i += n;
            if (i < 0 || i >= n)
              throw nb::index_error();
            v[i] = val;
          }
      )
      .def(
          "__iter__",
          [](CharacterAlloc1D &self) {
            // Convert to Python list and iterate (avoids StringProxy)
            nb::list result;
            auto &view = self.view();
            for (int i = 0; i < view.size(); ++i)
              result.append(view.get_string(view.lower_bound() + i));
            return result.attr("__iter__")();
          }
      )
      .def(
          "view",
          [](CharacterAlloc1D &self) { return self.view(); },
          nb::keep_alive<0, 1>()
      )
      .def(
          "__repr__",
          [](CharacterAlloc1D &self) {
            auto vec = self.to_vector();
            std::string result = "CharacterAlloc1D([";
            for (size_t i = 0; i < vec.size(); ++i) {
              if (i > 0)
                result += ", ";
              result += "'" + vec[i] + "'";
            }
            result += "])";
            return result;
          }
      )
      .def("__str__", [](CharacterAlloc1D &self) {
        auto vec = self.to_vector();
        std::string result = "[";
        for (size_t i = 0; i < vec.size(); ++i) {
          if (i > 0)
            result += ", ";
          result += "'" + vec[i] + "'";
        }
        result += "]";
        return result;
      });

  // Implicit conversion from list[str]
  nb::implicitly_convertible<std::vector<std::string>, CharacterAlloc1D>();
}

void bind_FCharArray1D(nb::module_ &m) {
  nb::class_<FCharArray1D>(m, "FCharArray1D")
      .def(nb::init<>())
      .def("is_valid", &FCharArray1D::is_valid)
      .def("__len__", &FCharArray1D::size)
      // Return std::string directly for Python convenience
      .def(
          "__getitem__",
          [](FCharArray1D &self, int i) {
            if (i < 0)
              i += self.size();
            if (i < 0 || i >= self.size())
              throw nb::index_error();
            return self[i];
          }
      )
      .def(
          "__setitem__",
          [](FCharArray1D &self, int i, const std::string &val) {
            if (i < 0)
              i += self.size();
            if (i < 0 || i >= self.size())
              throw nb::index_error();
            // FCharArray1D::operator[] returns proxy which implements operator=(string)
            self[i] = val;
          }
      )
      .def(
          "__iter__",
          [](FCharArray1D &self) {
            return nb::make_iterator(
                nb::type<FCharArray1D>(),
                "FCharArray1DIterator",
                self.begin(),
                self.end()
            );
          },
          nb::keep_alive<0, 1>()
      )
      .def("to_list", &FCharArray1D::to_vector)
      .def("__str__", [](const FCharArray1D &t) { return Bmad::to_string(t); });
}
void bind_standard_arrays(nb::module_ &m) {
  // 1. Primitive Arrays (View + Allocator + Vector Interop)
  //    This generates: RealArray1D, RealAlloc1D, and all permutations of conversion
  //    between them and python lists (std::vector).
  bind_1D_array_pair<double, RealAlloc1D>(m, "RealArray1D", "RealAlloc1D");

#ifdef __FLOAT128__
  bind_1D_array_pair<__float128, Real16Alloc1D>(m, "Real16Array1D", "Real16Alloc1D");
#else
  bind_1D_array_pair<long double, Real16Alloc1D>(m, "Real16Array1D", "Real16Alloc1D");
#endif

  bind_1D_array_pair<int, IntAlloc1D>(m, "IntArray1D", "IntAlloc1D");
  bind_1D_array_pair<int64_t, Int8Alloc1D>(m, "Int8Array1D", "Int8Alloc1D");
  bind_1D_array_pair<bool, BoolAlloc1D>(m, "BoolArray1D", "BoolAlloc1D");
  bind_1D_array_pair<std::complex<double>, ComplexAlloc1D>(m, "ComplexArray1D", "ComplexAlloc1D");

  bind_FArrayND<double, 2>(m, "RealArray2D");
  bind_FArrayND<int, 2>(m, "IntArray2D");
  bind_FArrayND<int64_t, 2>(m, "Int8Array2D");
  bind_FArrayND<bool, 2>(m, "BoolArray2D");
  bind_FArrayND<std::complex<double>, 2>(m, "ComplexArray2D");

  bind_FArrayND<double, 3>(m, "RealArray3D");
  bind_FArrayND<int, 3>(m, "IntArray3D");
  bind_FArrayND<int64_t, 3>(m, "Int8Array3D");
  bind_FArrayND<bool, 3>(m, "BoolArray3D");
  bind_FArrayND<std::complex<double>, 3>(m, "ComplexArray3D");

  // 3. String Arrays
  bind_FCharArray1D(m);
  bind_FCharAlloc1D(m);
}
} // namespace Pybmad
