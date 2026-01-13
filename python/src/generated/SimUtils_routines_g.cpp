#include "pybmad/generated/SimUtils_routines_g.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyGenCompleteElliptic python_gen_complete_elliptic(
    double kc,
    double p,
    double c,
    double s,
    std::optional<double> err_tol,
    double value) {
  SimUtils::gen_complete_elliptic(kc, p, c, s, make_opt_ref(err_tol), value);
  auto py_result{PyGenCompleteElliptic{kc, p, c, s, err_tol, value}};
  return py_result;
}
PyGetFileNumber python_get_file_number(
    std::string file_name,
    std::string cnum_in,
    int num_out,
    bool err_flag) {
  SimUtils::get_file_number(file_name, cnum_in, num_out, err_flag);
  auto py_result{PyGetFileNumber{file_name, cnum_in, num_out, err_flag}};
  return py_result;
}
PyGetFileTimeStamp python_get_file_time_stamp(
    std::string file,
    std::string time_stamp) {
  SimUtils::get_file_time_stamp(file, time_stamp);
  auto py_result{PyGetFileTimeStamp{file, time_stamp}};
  return py_result;
}

void init_SimUtils_routines_g(py::module& m) {
  py::class_<PyGenCompleteElliptic, std::unique_ptr<PyGenCompleteElliptic>>(
      m, "GenCompleteElliptic", "gen_complete_elliptic return type")
      .def_readonly("kc", &PyGenCompleteElliptic::kc)
      .def_readonly("p", &PyGenCompleteElliptic::p)
      .def_readonly("c", &PyGenCompleteElliptic::c)
      .def_readonly("s", &PyGenCompleteElliptic::s)
      .def_readonly("err_tol", &PyGenCompleteElliptic::err_tol)
      .def_readonly("value", &PyGenCompleteElliptic::value)
      .def("__len__", [](const PyGenCompleteElliptic&) { return 6; })
      .def(
          "__getitem__",
          [](const PyGenCompleteElliptic& s, int i) -> py::object {
            if (i < 0)
              i += 6;
            if (i == 0)
              return py::cast(s.kc);
            if (i == 1)
              return py::cast(s.p);
            if (i == 2)
              return py::cast(s.c);
            if (i == 3)
              return py::cast(s.s);
            if (i == 4)
              return py::cast(s.err_tol);
            if (i == 5)
              return py::cast(s.value);
            throw py::index_error();
          });
  m.def(
      "gen_complete_elliptic",
      &python_gen_complete_elliptic,
      py::arg("kc"),
      py::arg("p"),
      py::arg("c"),
      py::arg("s"),
      py::arg("err_tol") = py::none(),
      py::arg("value"),
      R"""(Parameters
  ----------
  kc : 
  p : 
  c : 
  s : 
  err_tol : 
  value : 
  )""");
  py::class_<PyGetFileNumber, std::unique_ptr<PyGetFileNumber>>(
      m, "GetFileNumber", "get_file_number return type")
      .def_readonly("file_name", &PyGetFileNumber::file_name)
      .def_readonly("cnum_in", &PyGetFileNumber::cnum_in)
      .def_readonly("num_out", &PyGetFileNumber::num_out)
      .def_readonly("err_flag", &PyGetFileNumber::err_flag)
      .def("__len__", [](const PyGetFileNumber&) { return 4; })
      .def("__getitem__", [](const PyGetFileNumber& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.file_name);
        if (i == 1)
          return py::cast(s.cnum_in);
        if (i == 2)
          return py::cast(s.num_out);
        if (i == 3)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "get_file_number",
      &python_get_file_number,
      py::arg("file_name"),
      py::arg("cnum_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      R"""(Parameters
  ----------
  file_name : 
  cnum_in : 
  num_out : 
  err_flag : 
  )""");
  py::class_<PyGetFileTimeStamp, std::unique_ptr<PyGetFileTimeStamp>>(
      m, "GetFileTimeStamp", "get_file_time_stamp return type")
      .def_readonly("file", &PyGetFileTimeStamp::file)
      .def_readonly("time_stamp", &PyGetFileTimeStamp::time_stamp)
      .def("__len__", [](const PyGetFileTimeStamp&) { return 2; })
      .def("__getitem__", [](const PyGetFileTimeStamp& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.file);
        if (i == 1)
          return py::cast(s.time_stamp);
        throw py::index_error();
      });
  m.def(
      "get_file_time_stamp",
      &python_get_file_time_stamp,
      py::arg("file"),
      py::arg("time_stamp"),
      R"""(no longer exists
  subroutine get_next_number (filein, cnum, digits)
    implicit none
    character(*) filein
    character(*) cnum
    integer digits
  end subroutine

  )""");
}
