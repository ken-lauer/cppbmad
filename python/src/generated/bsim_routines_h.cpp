#include "pybmad/generated/bsim_routines_h.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyHomVoltage {
  double voltage;
};
PyHomVoltage python_hom_voltage(WakeLrModeProxy& lr_wake, double voltage) {
  bsim::hom_voltage(lr_wake, voltage);
  auto py_result{PyHomVoltage{voltage}};
  return py_result;
}

void init_bsim_routines_h(py::module& m) {
  m.def(
      "hom_voltage",
      &python_hom_voltage,
      py::arg("lr_wake"),
      py::arg("voltage"),
      R"""(Parameters
  ----------
  lr_wake : 
  voltage : 
  )""");
  py::class_<PyHomVoltage, std::unique_ptr<PyHomVoltage>>(
      m, "HomVoltage", "Fortran routine hom_voltage return value")
      .def_readonly("voltage", &PyHomVoltage::voltage)
      .def("__len__", [](const PyHomVoltage&) { return 1; })
      .def("__getitem__", [](const PyHomVoltage& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.voltage);
        throw py::index_error();
      });
}
