#include "pybmad/generated/SimUtils_routines_l.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyLinearFit {
  int n_data;
  double a;
  double b;
  double sig_a;
  double sig_b;
};
PyLinearFit python_linear_fit(
    RealAlloc1D& x,
    RealAlloc1D& y,
    int n_data,
    double a,
    double b,
    double sig_a,
    double sig_b) {
  SimUtils::linear_fit(x, y, n_data, a, b, sig_a, sig_b);
  auto py_result{PyLinearFit{n_data, a, b, sig_a, sig_b}};
  return py_result;
}
struct PyLogicStr {
  bool logic;
  std::string str;
};
PyLogicStr python_logic_str(bool logic, std::string str) {
  SimUtils::logic_str(logic, str);
  auto py_result{PyLogicStr{logic, str}};
  return py_result;
}

void init_SimUtils_routines_l(py::module& m) {
  m.def(
      "linear_fit",
      &python_linear_fit,
      py::arg("x"),
      py::arg("y"),
      py::arg("n_data"),
      py::arg("a"),
      py::arg("b"),
      py::arg("sig_a"),
      py::arg("sig_b"),
      R"""(Parameters
  ----------
  x : 
  y : 
  n_data : 
  a : 
  b : 
  sig_a : 
  sig_b : 
  )""");
  py::class_<PyLinearFit, std::unique_ptr<PyLinearFit>>(
      m, "LinearFit", "Fortran routine linear_fit return value")
      .def_readonly("n_data", &PyLinearFit::n_data)
      .def_readonly("a", &PyLinearFit::a)
      .def_readonly("b", &PyLinearFit::b)
      .def_readonly("sig_a", &PyLinearFit::sig_a)
      .def_readonly("sig_b", &PyLinearFit::sig_b)
      .def("__len__", [](const PyLinearFit&) { return 5; })
      .def("__getitem__", [](const PyLinearFit& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.n_data);
        if (i == 1)
          return py::cast(s.a);
        if (i == 2)
          return py::cast(s.b);
        if (i == 3)
          return py::cast(s.sig_a);
        if (i == 4)
          return py::cast(s.sig_b);
        throw py::index_error();
      });
  m.def(
      "linear_fit_2d",
      &SimUtils::linear_fit_2d,
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      R"""(Parameters
  ----------
  x : float
      Array of x-values.
  y : float
      Array of y-values.
  z : float
      Array of z-values
  coef : float
      Coefficients of the linear fit
  )""");
  m.def(
      "logic_str",
      &python_logic_str,
      py::arg("logic"),
      py::arg("str"),
      R"""(Parameters
  ----------
  logic : 
  str : 
  )""");
  py::class_<PyLogicStr, std::unique_ptr<PyLogicStr>>(
      m, "LogicStr", "Fortran routine logic_str return value")
      .def_readonly("logic", &PyLogicStr::logic)
      .def_readonly("str", &PyLogicStr::str)
      .def("__len__", [](const PyLogicStr&) { return 2; })
      .def("__getitem__", [](const PyLogicStr& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.logic);
        if (i == 1)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "lunget",
      &SimUtils::lunget,
      R"""(Parameters
  ----------
  lunget : 
  )""");
}
