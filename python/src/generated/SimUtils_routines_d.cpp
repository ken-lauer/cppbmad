#include "pybmad/generated/SimUtils_routines_d.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyDateAndTimeStamp python_date_and_time_stamp(
    std::string string,
    std::optional<bool> numeric_month = std::nullopt,
    std::optional<bool> include_zone = std::nullopt) {
  SimUtils::date_and_time_stamp(
      string, make_opt_ref(numeric_month), make_opt_ref(include_zone));
  auto py_result{PyDateAndTimeStamp{string, numeric_month, include_zone}};
  return py_result;
}
PyDetab python_detab(std::string str) {
  SimUtils::detab(str);
  auto py_result{PyDetab{str}};
  return py_result;
}
PyDisplaySizeAndResolution python_display_size_and_resolution(
    int ix_screen,
    double x_size,
    double y_size,
    double x_res,
    double y_res) {
  SimUtils::display_size_and_resolution(
      ix_screen, x_size, y_size, x_res, y_res);
  auto py_result{
      PyDisplaySizeAndResolution{ix_screen, x_size, y_size, x_res, y_res}};
  return py_result;
}
PyDjBessel python_dj_bessel(int m, double arg, double dj_bes) {
  SimUtils::dj_bessel(m, arg, dj_bes);
  auto py_result{PyDjBessel{m, arg, dj_bes}};
  return py_result;
}
PyDjbHash python_djb_hash(
    std::string str,
    std::optional<int> old_hash,
    int hash) {
  SimUtils::djb_hash(str, make_opt_ref(old_hash), hash);
  auto py_result{PyDjbHash{str, old_hash, hash}};
  return py_result;
}
PyDjbStrHash python_djb_str_hash(std::string in_str, std::string hash_str) {
  SimUtils::djb_str_hash(in_str, hash_str);
  auto py_result{PyDjbStrHash{in_str, hash_str}};
  return py_result;
}
PyDowncaseString python_downcase_string(std::string string) {
  SimUtils::downcase_string(string);
  auto py_result{PyDowncaseString{string}};
  return py_result;
}

void init_SimUtils_routines_d(py::module& m) {
  py::class_<PyDateAndTimeStamp, std::unique_ptr<PyDateAndTimeStamp>>(
      m, "DateAndTimeStamp", "Fortran routine date_and_time_stamp return value")
      .def_readonly("string", &PyDateAndTimeStamp::string)
      .def_readonly("numeric_month", &PyDateAndTimeStamp::numeric_month)
      .def_readonly("include_zone", &PyDateAndTimeStamp::include_zone)
      .def("__len__", [](const PyDateAndTimeStamp&) { return 3; })
      .def("__getitem__", [](const PyDateAndTimeStamp& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.string);
        if (i == 1)
          return py::cast(s.numeric_month);
        if (i == 2)
          return py::cast(s.include_zone);
        throw py::index_error();
      });
  m.def(
      "date_and_time_stamp",
      &python_date_and_time_stamp,
      py::arg("string"),
      py::arg("numeric_month") = py::none(),
      py::arg("include_zone") = py::none(),
      R"""(Parameters
  ----------
  string : 
  numeric_month : 
  include_zone : 
  )""");
  m.def(
      "destfixedwindowls",
      &SimUtils::destfixedwindowls,
      py::arg("id"),
      R"""(Parameters
  ----------
  id : 
  )""");
  py::class_<PyDetab, std::unique_ptr<PyDetab>>(
      m, "Detab", "Fortran routine detab return value")
      .def_readonly("str", &PyDetab::str)
      .def("__len__", [](const PyDetab&) { return 1; })
      .def("__getitem__", [](const PyDetab& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "detab",
      &python_detab,
      py::arg("str"),
      R"""(Parameters
  ----------
  str : 
  )""");
  py::class_<
      PyDisplaySizeAndResolution,
      std::unique_ptr<PyDisplaySizeAndResolution>>(
      m,
      "DisplaySizeAndResolution",
      "Fortran routine display_size_and_resolution return value")
      .def_readonly("ix_screen", &PyDisplaySizeAndResolution::ix_screen)
      .def_readonly("x_size", &PyDisplaySizeAndResolution::x_size)
      .def_readonly("y_size", &PyDisplaySizeAndResolution::y_size)
      .def_readonly("x_res", &PyDisplaySizeAndResolution::x_res)
      .def_readonly("y_res", &PyDisplaySizeAndResolution::y_res)
      .def("__len__", [](const PyDisplaySizeAndResolution&) { return 5; })
      .def(
          "__getitem__",
          [](const PyDisplaySizeAndResolution& s, int i) -> py::object {
            if (i < 0)
              i += 5;
            if (i == 0)
              return py::cast(s.ix_screen);
            if (i == 1)
              return py::cast(s.x_size);
            if (i == 2)
              return py::cast(s.y_size);
            if (i == 3)
              return py::cast(s.x_res);
            if (i == 4)
              return py::cast(s.y_res);
            throw py::index_error();
          });
  m.def(
      "display_size_and_resolution",
      &python_display_size_and_resolution,
      py::arg("ix_screen"),
      py::arg("x_size"),
      py::arg("y_size"),
      py::arg("x_res"),
      py::arg("y_res"),
      R"""(Parameters
  ----------
  ix_screen : 
  x_size : 
  y_size : 
  x_res : 
  y_res : 
  )""");
  py::class_<PyDjBessel, std::unique_ptr<PyDjBessel>>(
      m, "DjBessel", "Fortran routine dj_bessel return value")
      .def_readonly("m", &PyDjBessel::m)
      .def_readonly("arg", &PyDjBessel::arg)
      .def_readonly("dj_bes", &PyDjBessel::dj_bes)
      .def("__len__", [](const PyDjBessel&) { return 3; })
      .def("__getitem__", [](const PyDjBessel& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.m);
        if (i == 1)
          return py::cast(s.arg);
        if (i == 2)
          return py::cast(s.dj_bes);
        throw py::index_error();
      });
  m.def(
      "dj_bessel",
      &python_dj_bessel,
      py::arg("m"),
      py::arg("arg"),
      py::arg("dj_bes"),
      R"""(Parameters
  ----------
  m : 
  arg : 
  dj_bes : 
  )""");
  py::class_<PyDjbHash, std::unique_ptr<PyDjbHash>>(
      m, "DjbHash", "Fortran routine djb_hash return value")
      .def_readonly("str", &PyDjbHash::str)
      .def_readonly("old_hash", &PyDjbHash::old_hash)
      .def_readonly("hash", &PyDjbHash::hash)
      .def("__len__", [](const PyDjbHash&) { return 3; })
      .def("__getitem__", [](const PyDjbHash& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.str);
        if (i == 1)
          return py::cast(s.old_hash);
        if (i == 2)
          return py::cast(s.hash);
        throw py::index_error();
      });
  m.def(
      "djb_hash",
      &python_djb_hash,
      py::arg("str"),
      py::arg("old_hash") = py::none(),
      py::arg("hash"),
      R"""(Parameters
  ----------
  str : 
  old_hash : 
  hash : 
  )""");
  py::class_<PyDjbStrHash, std::unique_ptr<PyDjbStrHash>>(
      m, "DjbStrHash", "Fortran routine djb_str_hash return value")
      .def_readonly("in_str", &PyDjbStrHash::in_str)
      .def_readonly("hash_str", &PyDjbStrHash::hash_str)
      .def("__len__", [](const PyDjbStrHash&) { return 2; })
      .def("__getitem__", [](const PyDjbStrHash& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.in_str);
        if (i == 1)
          return py::cast(s.hash_str);
        throw py::index_error();
      });
  m.def(
      "djb_str_hash",
      &python_djb_str_hash,
      py::arg("in_str"),
      py::arg("hash_str"),
      R"""(Parameters
  ----------
  in_str : 
  hash_str : 
  )""");
  py::class_<PyDowncaseString, std::unique_ptr<PyDowncaseString>>(
      m, "DowncaseString", "Fortran routine downcase_string return value")
      .def_readonly("string", &PyDowncaseString::string)
      .def("__len__", [](const PyDowncaseString&) { return 1; })
      .def("__getitem__", [](const PyDowncaseString& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.string);
        throw py::index_error();
      });
  m.def(
      "downcase_string",
      &python_downcase_string,
      py::arg("string"),
      R"""(Parameters
  ----------
  string : 
  )""");
}
