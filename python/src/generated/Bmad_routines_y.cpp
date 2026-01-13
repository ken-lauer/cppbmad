#include "pybmad/generated/Bmad_routines_y.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyYlafun python_ylafun(double x, double y, double z, double res) {
  Bmad::ylafun(x, y, z, res);
  auto py_result{PyYlafun{x, y, z, res}};
  return py_result;
}

void init_Bmad_routines_y(py::module& m) {
  py::class_<PyYlafun, std::unique_ptr<PyYlafun>>(
      m, "Ylafun", "ylafun return type")
      .def_readonly("x", &PyYlafun::x)
      .def_readonly("y", &PyYlafun::y)
      .def_readonly("z", &PyYlafun::z)
      .def_readonly("res", &PyYlafun::res)
      .def("__len__", [](const PyYlafun&) { return 4; })
      .def("__getitem__", [](const PyYlafun& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.x);
        if (i == 1)
          return py::cast(s.y);
        if (i == 2)
          return py::cast(s.z);
        if (i == 3)
          return py::cast(s.res);
        throw py::index_error();
      });
  m.def(
      "ylafun",
      &python_ylafun,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::arg("res"),
      R"""(Parameters
  ----------
  x : 
  y : 
  z : 
  res : 
  )""");
}
