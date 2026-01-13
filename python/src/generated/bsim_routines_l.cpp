#include "pybmad/generated/bsim_routines_l.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyLogicalToPython python_logical_to_python(bool logic, std::string string) {
  bsim::logical_to_python(logic, string);
  auto py_result{PyLogicalToPython{logic, string}};
  return py_result;
}

void init_bsim_routines_l(py::module& m) {
  py::class_<PyLogicalToPython, std::unique_ptr<PyLogicalToPython>>(
      m, "LogicalToPython", "logical_to_python return type")
      .def_readonly("logic", &PyLogicalToPython::logic)
      .def_readonly("string", &PyLogicalToPython::string)
      .def("__len__", [](const PyLogicalToPython&) { return 2; })
      .def("__getitem__", [](const PyLogicalToPython& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.logic);
        if (i == 1)
          return py::cast(s.string);
        throw py::index_error();
      });
  m.def(
      "logical_to_python",
      &python_logical_to_python,
      py::arg("logic"),
      py::arg("string"),
      R"""(Parameters
  ----------
  logic : 
  string : 
  )""");
}
