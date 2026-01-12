#include "pybmad/generated/SimUtils_routines_h.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyHanhan python_hanhan(int N, RealAlloc1D& hh) {
  SimUtils::hanhan(N, hh);
  auto py_result{PyHanhan{N}};
  return py_result;
}

void init_SimUtils_routines_h(py::module& m) {
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
}
