#include "pybmad/generated/SimUtils_routines_u.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyUpcaseString {
  std::string string;
};
PyUpcaseString python_upcase_string(std::string string) {
  SimUtils::upcase_string(string);
  auto py_result{PyUpcaseString{string}};
  return py_result;
}

void init_SimUtils_routines_u(py::module& m) {
  m.def(
      "upcase_string",
      &python_upcase_string,
      py::arg("string"),
      R"""(Parameters
  ----------
  string : 
  )""");
  py::class_<PyUpcaseString, std::unique_ptr<PyUpcaseString>>(
      m, "UpcaseString", "Fortran routine upcase_string return value")
      .def_readonly("string", &PyUpcaseString::string)
      .def("__len__", [](const PyUpcaseString&) { return 1; })
      .def("__getitem__", [](const PyUpcaseString& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.string);
        throw py::index_error();
      });
}
