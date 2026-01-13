#include "pybmad/generated/SimUtils_routines_j.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyJBessel python_j_bessel(int m, double arg, double j_bes) {
  SimUtils::j_bessel(m, arg, j_bes);
  auto py_result{PyJBessel{m, arg, j_bes}};
  return py_result;
}

void init_SimUtils_routines_j(py::module& m) {
  py::class_<PyJBessel, std::unique_ptr<PyJBessel>>(
      m, "JBessel", "j_bessel return type")
      .def_readonly("m", &PyJBessel::m)
      .def_readonly("arg", &PyJBessel::arg)
      .def_readonly("j_bes", &PyJBessel::j_bes)
      .def("__len__", [](const PyJBessel&) { return 3; })
      .def("__getitem__", [](const PyJBessel& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.m);
        if (i == 1)
          return py::cast(s.arg);
        if (i == 2)
          return py::cast(s.j_bes);
        throw py::index_error();
      });
  m.def(
      "j_bessel",
      &python_j_bessel,
      py::arg("m"),
      py::arg("arg"),
      py::arg("j_bes"),
      R"""(Parameters
  ----------
  m : 
  arg : 
  j_bes : 
  )""");
}
