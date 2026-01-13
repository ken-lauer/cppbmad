#include "pybmad/generated/SimUtils_routines_p.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyParseFortranFormat python_parse_fortran_format(
    std::string format_str,
    int n_repeat,
    int power,
    std::string descrip,
    int width,
    int digits) {
  SimUtils::parse_fortran_format(
      format_str, n_repeat, power, descrip, width, digits);
  auto py_result{PyParseFortranFormat{
      format_str, n_repeat, power, descrip, width, digits}};
  return py_result;
}
PyPolyEval python_poly_eval(
    RealAlloc1D& poly,
    double x,
    std::optional<bool> diff_coef,
    double y) {
  SimUtils::poly_eval(poly, x, diff_coef, y);
  auto py_result{PyPolyEval{y}};
  return py_result;
}
PyProbabilityFunct python_probability_funct(double x, double prob) {
  SimUtils::probability_funct(x, prob);
  auto py_result{PyProbabilityFunct{prob}};
  return py_result;
}
PyProjdd python_projdd(
    ComplexAlloc1D& a,
    ComplexAlloc1D& b,
    std::complex<double> func_retval__) {
  SimUtils::projdd(a, b, func_retval__);
  auto py_result{PyProjdd{func_retval__}};
  return py_result;
}

void init_SimUtils_routines_p(py::module& m) {
  py::class_<PyParseFortranFormat, std::unique_ptr<PyParseFortranFormat>>(
      m, "ParseFortranFormat", "parse_fortran_format return type")
      .def_readonly("format_str", &PyParseFortranFormat::format_str)
      .def_readonly("n_repeat", &PyParseFortranFormat::n_repeat)
      .def_readonly("power", &PyParseFortranFormat::power)
      .def_readonly("descrip", &PyParseFortranFormat::descrip)
      .def_readonly("width", &PyParseFortranFormat::width)
      .def_readonly("digits", &PyParseFortranFormat::digits)
      .def("__len__", [](const PyParseFortranFormat&) { return 6; })
      .def(
          "__getitem__",
          [](const PyParseFortranFormat& s, int i) -> py::object {
            if (i < 0)
              i += 6;
            if (i == 0)
              return py::cast(s.format_str);
            if (i == 1)
              return py::cast(s.n_repeat);
            if (i == 2)
              return py::cast(s.power);
            if (i == 3)
              return py::cast(s.descrip);
            if (i == 4)
              return py::cast(s.width);
            if (i == 5)
              return py::cast(s.digits);
            throw py::index_error();
          });
  m.def(
      "parse_fortran_format",
      &python_parse_fortran_format,
      py::arg("format_str"),
      py::arg("n_repeat"),
      py::arg("power"),
      py::arg("descrip"),
      py::arg("width"),
      py::arg("digits"),
      R"""(Parameters
  ----------
  format_str : 
  n_repeat : 
  power : 
  descrip : 
  width : 
  digits : 
  )""");
  m.def(
      "pointer_to_ran_state",
      &SimUtils::pointer_to_ran_state,
      py::arg("ran_state") = py::none(),
      py::arg("ix_thread") = py::none(),
      R"""(Function pointer_to_ran_state(ran_state, ix_thread) result (ran_state_ptr)

  Routine to point to the appropriate state structure for generating random numbers

  Parameters
  ----------
  ran_state : RandomStateStruct, optional
      Point to this if present. Otherwise point to the global saved state.
  ix_thread : int, optional
      Thread index.

  Returns
  -------
  ran_state_ptr : RandomStateStruct
      Pointer to the appropriate state.
  )""");
  py::class_<PyPolyEval, std::unique_ptr<PyPolyEval>>(
      m, "PolyEval", "poly_eval return type")
      .def_readonly("y", &PyPolyEval::y)
      .def("__len__", [](const PyPolyEval&) { return 1; })
      .def("__getitem__", [](const PyPolyEval& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.y);
        throw py::index_error();
      });
  m.def(
      "poly_eval",
      &python_poly_eval,
      py::arg("poly"),
      py::arg("x"),
      py::arg("diff_coef") = py::none(),
      py::arg("y"),
      R"""(Parameters
  ----------
  poly : float
      Polynomial
  x : float
      Point to evaluate at.
  diff_coef : bool, optional
      poly(:) array are differentials? Default is False.
  y : 
  )""");
  py::class_<PyProbabilityFunct, std::unique_ptr<PyProbabilityFunct>>(
      m, "ProbabilityFunct", "probability_funct return type")
      .def_readonly("prob", &PyProbabilityFunct::prob)
      .def("__len__", [](const PyProbabilityFunct&) { return 1; })
      .def("__getitem__", [](const PyProbabilityFunct& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.prob);
        throw py::index_error();
      });
  m.def(
      "probability_funct",
      &python_probability_funct,
      py::arg("x"),
      py::arg("prob"),
      R"""(Parameters
  ----------
  x : float
      Function argument.
  prob : 
  )""");
  py::class_<PyProjdd, std::unique_ptr<PyProjdd>>(
      m, "Projdd", "projdd return type")
      .def_readonly("func_retval__", &PyProjdd::func_retval__)
      .def("__len__", [](const PyProjdd&) { return 1; })
      .def("__getitem__", [](const PyProjdd& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.func_retval__);
        throw py::index_error();
      });
  m.def(
      "projdd",
      &python_projdd,
      py::arg("a"),
      py::arg("b"),
      py::arg("func_retval__"),
      R"""(Parameters
  ----------
  a : 
  b : 
  projdd : 
  )""");
}
