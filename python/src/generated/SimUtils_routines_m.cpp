#include "pybmad/generated/SimUtils_routines_m.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyMakeLegalComment python_make_legal_comment(
    std::string comment_in,
    std::string comment_out) {
  SimUtils::make_legal_comment(comment_in, comment_out);
  auto py_result{PyMakeLegalComment{comment_in, comment_out}};
  return py_result;
}
PyMatchReg python_match_reg(std::string str, std::string pat, bool is_match) {
  SimUtils::match_reg(str, pat, is_match);
  auto py_result{PyMatchReg{str, pat, is_match}};
  return py_result;
}
PyMatchWild python_match_wild(
    std::string string,
    std::string template_,
    bool is_match) {
  SimUtils::match_wild(string, template_, is_match);
  auto py_result{PyMatchWild{string, template_, is_match}};
  return py_result;
}
PyMaximizeProjection python_maximize_projection(
    double seed,
    ComplexAlloc1D& cdata,
    double func_retval__) {
  SimUtils::maximize_projection(seed, cdata, func_retval__);
  auto py_result{PyMaximizeProjection{seed, func_retval__}};
  return py_result;
}
PyMilliSleep python_milli_sleep(int milli_sec) {
  SimUtils::milli_sleep(milli_sec);
  auto py_result{PyMilliSleep{milli_sec}};
  return py_result;
}

void init_SimUtils_routines_m(py::module& m) {
  py::class_<PyMakeLegalComment, std::unique_ptr<PyMakeLegalComment>>(
      m, "MakeLegalComment", "Fortran routine make_legal_comment return value")
      .def_readonly("comment_in", &PyMakeLegalComment::comment_in)
      .def_readonly("comment_out", &PyMakeLegalComment::comment_out)
      .def("__len__", [](const PyMakeLegalComment&) { return 2; })
      .def("__getitem__", [](const PyMakeLegalComment& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.comment_in);
        if (i == 1)
          return py::cast(s.comment_out);
        throw py::index_error();
      });
  m.def(
      "make_legal_comment",
      &python_make_legal_comment,
      py::arg("comment_in"),
      py::arg("comment_out"),
      R"""(Parameters
  ----------
  comment_in : 
  comment_out : 
  )""");
  m.def(
      "mass_of",
      &SimUtils::mass_of,
      py::arg("species"),
      R"""(Function mass_of (species) result (mass)

  Routine to return the mass, in units of eV/c^2, of a particle.
  To convert to AMU divide mass_of value by the constant atomic_mass_unit.

  Note: For atoms where the isotopic number is given, the mass is calculated using the neutral atomic mass
  adjusted by the weight of any added or missing electrons. The calculated mass is off very slightly due to
  binding energy effects. Exception: For #1H+ (proton) and #2H+ (deuteron) the exact mass is used since it is known.

  Parameters
  ----------
  species : int
      Species ID.

  Returns
  -------
  mass : float
      particle mass. Set to real_garbage$ if species value is invalid.
  )""");
  py::class_<PyMatchReg, std::unique_ptr<PyMatchReg>>(
      m, "MatchReg", "Fortran routine match_reg return value")
      .def_readonly("str", &PyMatchReg::str)
      .def_readonly("pat", &PyMatchReg::pat)
      .def_readonly("is_match", &PyMatchReg::is_match)
      .def("__len__", [](const PyMatchReg&) { return 3; })
      .def("__getitem__", [](const PyMatchReg& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.str);
        if (i == 1)
          return py::cast(s.pat);
        if (i == 2)
          return py::cast(s.is_match);
        throw py::index_error();
      });
  m.def(
      "match_reg",
      &python_match_reg,
      py::arg("str"),
      py::arg("pat"),
      py::arg("is_match"),
      R"""(Parameters
  ----------
  str : 
  pat : 
  is_match : 
  )""");
  py::class_<PyMatchWild, std::unique_ptr<PyMatchWild>>(
      m, "MatchWild", "Fortran routine match_wild return value")
      .def_readonly("string", &PyMatchWild::string)
      .def_readonly("template_", &PyMatchWild::template_)
      .def_readonly("is_match", &PyMatchWild::is_match)
      .def("__len__", [](const PyMatchWild&) { return 3; })
      .def("__getitem__", [](const PyMatchWild& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.template_);
        if (i == 2)
          return py::cast(s.is_match);
        throw py::index_error();
      });
  m.def(
      "match_wild",
      &python_match_wild,
      py::arg("string"),
      py::arg("template_"),
      py::arg("is_match"),
      R"""(Parameters
  ----------
  string : 
  template : 
  is_match : 
  )""");
  py::class_<PyMaximizeProjection, std::unique_ptr<PyMaximizeProjection>>(
      m,
      "MaximizeProjection",
      "Fortran routine maximize_projection return value")
      .def_readonly("seed", &PyMaximizeProjection::seed)
      .def_readonly("func_retval__", &PyMaximizeProjection::func_retval__)
      .def("__len__", [](const PyMaximizeProjection&) { return 2; })
      .def(
          "__getitem__",
          [](const PyMaximizeProjection& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.seed);
            if (i == 1)
              return py::cast(s.func_retval__);
            throw py::index_error();
          });
  m.def(
      "maximize_projection",
      &python_maximize_projection,
      py::arg("seed"),
      py::arg("cdata"),
      py::arg("func_retval__"),
      R"""(function maximize_projection

  Optimizer that uses Numerical Recipes brent to find a local maximum,
  which is the frequency that maximizes the projection.


  )""");
  py::class_<PyMilliSleep, std::unique_ptr<PyMilliSleep>>(
      m, "MilliSleep", "Fortran routine milli_sleep return value")
      .def_readonly("milli_sec", &PyMilliSleep::milli_sec)
      .def("__len__", [](const PyMilliSleep&) { return 1; })
      .def("__getitem__", [](const PyMilliSleep& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.milli_sec);
        throw py::index_error();
      });
  m.def(
      "milli_sleep",
      &python_milli_sleep,
      py::arg("milli_sec"),
      R"""(Parameters
  ----------
  milli_sec : 
  )""");
}
