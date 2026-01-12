#include "pybmad/generated/Bmad_routines_x.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyXlafun python_xlafun(double x, double y, double z, double res) {
  Bmad::xlafun(x, y, z, res);
  auto py_result{PyXlafun{x, y, z, res}};
  return py_result;
}

void init_Bmad_routines_x(py::module& m) {
  py::class_<PyXlafun, std::unique_ptr<PyXlafun>>(
      m, "Xlafun", "Fortran routine xlafun return value")
      .def_readonly("x", &PyXlafun::x)
      .def_readonly("y", &PyXlafun::y)
      .def_readonly("z", &PyXlafun::z)
      .def_readonly("res", &PyXlafun::res)
      .def("__len__", [](const PyXlafun&) { return 4; })
      .def("__getitem__", [](const PyXlafun& s, int i) -> py::object {
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
      "xlafun",
      &python_xlafun,
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
  m.def(
      "xraylib_nist_compound",
      &Bmad::xraylib_nist_compound,
      py::arg("name"),
      R"""(Function xraylib_nist_compound (name) result (indx)

  Routine to return the xraylib index for a given NIST compound.
  Taken from file xraylib/include/xraylib-nist_compounds.h

  Parameters
  ----------
  name : unknown
      Name of compound

  Returns
  -------
  indx : int
      Compound index. -1 if not found.
  )""");
}
