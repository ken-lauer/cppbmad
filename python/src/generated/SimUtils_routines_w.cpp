#include "pybmad/generated/SimUtils_routines_w.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyWordLen {
  std::string wording;
  int wlen;
};
PyWordLen python_word_len(std::string wording, int wlen) {
  SimUtils::word_len(wording, wlen);
  auto py_result{PyWordLen{wording, wlen}};
  return py_result;
}
struct PyWordRead {
  std::string in_str;
  std::string delim_list;
  std::string word;
  int ix_word;
  std::string delim;
  bool delim_found;
  std::string out_str;
  std::optional<bool> ignore_interior;
};
PyWordRead python_word_read(
    std::string in_str,
    std::string delim_list,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    std::string out_str,
    std::optional<bool> ignore_interior = std::nullopt) {
  SimUtils::word_read(
      in_str,
      delim_list,
      word,
      ix_word,
      delim,
      delim_found,
      out_str,
      make_opt_ref(ignore_interior));
  auto py_result{PyWordRead{
      in_str,
      delim_list,
      word,
      ix_word,
      delim,
      delim_found,
      out_str,
      ignore_interior}};
  return py_result;
}

void init_SimUtils_routines_w(py::module& m) {
  m.def(
      "w_mat_to_axis_angle",
      &SimUtils::w_mat_to_axis_angle,
      py::arg("w_mat"),
      R"""(Subroutine w_mat_to_axis_angle (w_mat, axis, angle)

  Routine to find the rotation axis and rotation angle corresponding to a given
  3D rotation matrix.

  The rotation angle is chosen in the range [0, pi].

  Parameters
  ----------
  w_mat : float
      Rotation matrix.

  Returns
  -------
  axis : float
      Rotation axis. Normalized to 1.
  angle : float
      Rotation angle in the range [0, pi].
  )""");
  py::class_<
      SimUtils::WMatToAxisAngle,
      std::unique_ptr<SimUtils::WMatToAxisAngle>>(
      m, "WMatToAxisAngle", "Fortran routine w_mat_to_axis_angle return value")
      .def_readonly("axis", &SimUtils::WMatToAxisAngle::axis)
      .def_readonly("angle", &SimUtils::WMatToAxisAngle::angle)
      .def("__len__", [](const SimUtils::WMatToAxisAngle&) { return 2; })
      .def(
          "__getitem__",
          [](const SimUtils::WMatToAxisAngle& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.axis);
            if (i == 1)
              return py::cast(s.angle);
            throw py::index_error();
          });
  m.def(
      "w_mat_to_quat",
      &SimUtils::w_mat_to_quat,
      py::arg("w_mat"),
      R"""(Function w_mat_to_quat (w_mat) result (quat)

  Routine to find the quaternion corresponding to a given 3D rotation matrix.

  Parameters
  ----------
  w_mat : float
      Rotation matrix

  Returns
  -------
  quat : float
      Quaternion.
  )""");
  m.def(
      "word_len",
      &python_word_len,
      py::arg("wording"),
      py::arg("wlen"),
      R"""(Parameters
  ----------
  wording : 
  wlen : 
  )""");
  py::class_<PyWordLen, std::unique_ptr<PyWordLen>>(
      m, "WordLen", "Fortran routine word_len return value")
      .def_readonly("wording", &PyWordLen::wording)
      .def_readonly("wlen", &PyWordLen::wlen)
      .def("__len__", [](const PyWordLen&) { return 2; })
      .def("__getitem__", [](const PyWordLen& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.wording);
        if (i == 1)
          return py::cast(s.wlen);
        throw py::index_error();
      });
  m.def(
      "word_read",
      &python_word_read,
      py::arg("in_str"),
      py::arg("delim_list"),
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("out_str"),
      py::arg("ignore_interior") = py::none(),
      R"""(Parameters
  ----------
  in_str : 
  delim_list : 
  word : 
  ix_word : 
  delim : 
  delim_found : 
  out_str : 
  ignore_interior : 
  )""");
  py::class_<PyWordRead, std::unique_ptr<PyWordRead>>(
      m, "WordRead", "Fortran routine word_read return value")
      .def_readonly("in_str", &PyWordRead::in_str)
      .def_readonly("delim_list", &PyWordRead::delim_list)
      .def_readonly("word", &PyWordRead::word)
      .def_readonly("ix_word", &PyWordRead::ix_word)
      .def_readonly("delim", &PyWordRead::delim)
      .def_readonly("delim_found", &PyWordRead::delim_found)
      .def_readonly("out_str", &PyWordRead::out_str)
      .def_readonly("ignore_interior", &PyWordRead::ignore_interior)
      .def("__len__", [](const PyWordRead&) { return 8; })
      .def("__getitem__", [](const PyWordRead& s, int i) -> py::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return py::cast(s.in_str);
        if (i == 1)
          return py::cast(s.delim_list);
        if (i == 2)
          return py::cast(s.word);
        if (i == 3)
          return py::cast(s.ix_word);
        if (i == 4)
          return py::cast(s.delim);
        if (i == 5)
          return py::cast(s.delim_found);
        if (i == 6)
          return py::cast(s.out_str);
        if (i == 7)
          return py::cast(s.ignore_interior);
        throw py::index_error();
      });
}
