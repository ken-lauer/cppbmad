#include "pybmad/generated/SimUtils_routines_h.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyHanhan {
  int N;
};
PyHanhan python_hanhan(int N, RealAlloc1D& hh) {
  SimUtils::hanhan(N, hh);
  auto py_result{PyHanhan{N}};
  return py_result;
}

void init_SimUtils_routines_h(py::module& m) {
  m.def(
      "hanhan",
      &python_hanhan,
      py::arg("N"),
      py::arg("hh"),
      R"""(Parameters
  ----------
  N : 
  hh : 
  )""");
  py::class_<PyHanhan, std::unique_ptr<PyHanhan>>(
      m, "Hanhan", "Fortran routine hanhan return value")
      .def_readonly("N", &PyHanhan::N)
      .def("__len__", [](const PyHanhan&) { return 1; })
      .def("__getitem__", [](const PyHanhan& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.N);
        throw py::index_error();
      });
}
