#include "pybmad/arrays.hpp"

namespace Pybmad {

inline int normalize_index(int i, int size) {
  if (i < 0)
    i += size;
  if (i < 0 || i >= size)
    throw py::index_error();
  return i;
}
void bind_FCharArray1D(py::module& m) {
  py::class_<FCharArray1D>(m, "FCharArray1D")
      .def(py::init<>())
      .def("is_valid", &FCharArray1D::is_valid)
      .def("__len__", &FCharArray1D::size)
      // Return std::string directly for Python convenience
      .def(
          "__getitem__",
          [](FCharArray1D& self, int i) {
            if (i < 0)
              i += self.size();
            if (i < 0 || i >= self.size())
              throw py::index_error();
            return self[i];
          })
      .def(
          "__setitem__",
          [](FCharArray1D& self, int i, const std::string& val) {
            if (i < 0)
              i += self.size();
            if (i < 0 || i >= self.size())
              throw py::index_error();
            // FCharArray1D::operator[] returns proxy which implements operator=(string)
            self[i] = val;
          })
      .def(
          "__iter__",
          [](FCharArray1D& self) {
            // Helper to make iterator that yields strings instead of StringProxy objects
            // to make Python life easier
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def("to_list", &FCharArray1D::to_vector)
      .def("__str__", [](const FCharArray1D& t) { return Bmad::to_string(t); });
}
void bind_standard_arrays(py::module& m) {
  bind_FArray1D<double>(m, "RealArray1D");
#ifdef __FLOAT128__
  bind_FArray1D<__float128>(m, "Real16Array1D");
#else
  bind_FArray1D<long double>(m, "Real16Array1D");
#endif
  bind_FArray1D<int>(m, "IntArray1D");
  bind_FArray1D<int64_t>(m, "Int8Array1D");
  bind_FArray1D<bool>(m, "BoolArray1D");
  bind_FArray1D<std::complex<double>>(m, "ComplexArray1D");

  bind_FArrayND<double, 2>(m, "RealArray2D");
  bind_FArrayND<int, 2>(m, "IntArray2D");
  bind_FArrayND<bool, 2>(m, "BoolArray2D");
  bind_FArrayND<std::complex<double>, 2>(m, "ComplexArray2D");

  bind_FArrayND<double, 3>(m, "RealArray3D");
  bind_FArrayND<int, 3>(m, "IntArray3D");
  bind_FArrayND<bool, 3>(m, "BoolArray3D");
  bind_FArrayND<std::complex<double>, 3>(m, "ComplexArray3D");

  bind_FAlloc1D<RealAlloc1D>(m, "RealAlloc1D");
  bind_FAlloc1D<IntAlloc1D>(m, "IntAlloc1D");
  bind_FAlloc1D<Int8Alloc1D>(m, "Int8Alloc1D");
  //!! this is the type, duh
  bind_FAlloc1D<BoolAlloc1D>(m, "BoolAlloc1D");
  bind_FAlloc1D<ComplexAlloc1D>(m, "ComplexAlloc1D");
  bind_FAlloc1D<Real16Alloc1D>(m, "Real16Alloc1D");

  // 3. String Arrays
  bind_FCharArray1D(m);
}
} // namespace Pybmad
