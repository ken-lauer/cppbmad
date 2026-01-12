#include "pybmad/generated/Tao_routines_i.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyIntegrateMax python_integrate_max(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum) {
  Tao::integrate_max(ix_start, ix_ele, datum_value, ix_m, branch, vec, datum);
  auto py_result{PyIntegrateMax{ix_start, ix_ele, datum_value, ix_m}};
  return py_result;
}
PyIntegrateMin python_integrate_min(
    int ix_start,
    int ix_ele,
    double datum_value,
    int ix_m,
    BranchProxy& branch,
    RealAlloc1D& vec,
    TaoDataProxy& datum) {
  Tao::integrate_min(ix_start, ix_ele, datum_value, ix_m, branch, vec, datum);
  auto py_result{PyIntegrateMin{ix_start, ix_ele, datum_value, ix_m}};
  return py_result;
}

void init_Tao_routines_i(py::module& m) {
  py::class_<PyIntegrateMax, std::unique_ptr<PyIntegrateMax>>(
      m, "IntegrateMax", "Fortran routine integrate_max return value")
      .def_readonly("ix_start", &PyIntegrateMax::ix_start)
      .def_readonly("ix_ele", &PyIntegrateMax::ix_ele)
      .def_readonly("datum_value", &PyIntegrateMax::datum_value)
      .def_readonly("ix_m", &PyIntegrateMax::ix_m)
      .def("__len__", [](const PyIntegrateMax&) { return 4; })
      .def("__getitem__", [](const PyIntegrateMax& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.ix_start);
        if (i == 1)
          return py::cast(s.ix_ele);
        if (i == 2)
          return py::cast(s.datum_value);
        if (i == 3)
          return py::cast(s.ix_m);
        throw py::index_error();
      });
  m.def(
      "integrate_max",
      &python_integrate_max,
      py::arg("ix_start"),
      py::arg("ix_ele"),
      py::arg("datum_value"),
      py::arg("ix_m"),
      py::arg("branch"),
      py::arg("vec"),
      py::arg("datum"),
      R"""(Parameters
  ----------
  ix_start : 
  ix_ele : 
  datum_value : 
  ix_m : 
  branch : 
  vec : 
  datum : 
  )""");
  py::class_<PyIntegrateMin, std::unique_ptr<PyIntegrateMin>>(
      m, "IntegrateMin", "Fortran routine integrate_min return value")
      .def_readonly("ix_start", &PyIntegrateMin::ix_start)
      .def_readonly("ix_ele", &PyIntegrateMin::ix_ele)
      .def_readonly("datum_value", &PyIntegrateMin::datum_value)
      .def_readonly("ix_m", &PyIntegrateMin::ix_m)
      .def("__len__", [](const PyIntegrateMin&) { return 4; })
      .def("__getitem__", [](const PyIntegrateMin& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.ix_start);
        if (i == 1)
          return py::cast(s.ix_ele);
        if (i == 2)
          return py::cast(s.datum_value);
        if (i == 3)
          return py::cast(s.ix_m);
        throw py::index_error();
      });
  m.def(
      "integrate_min",
      &python_integrate_min,
      py::arg("ix_start"),
      py::arg("ix_ele"),
      py::arg("datum_value"),
      py::arg("ix_m"),
      py::arg("branch"),
      py::arg("vec"),
      py::arg("datum"),
      R"""(Parameters
  ----------
  ix_start : 
  ix_ele : 
  datum_value : 
  ix_m : 
  branch : 
  vec : 
  datum : 
  )""");
}
