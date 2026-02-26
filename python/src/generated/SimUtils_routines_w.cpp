#include "pybmad/generated/SimUtils_routines_w.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_w(py::module &m) {
  py::class_<SimUtils::WMatToAxisAngle, std::unique_ptr<SimUtils::WMatToAxisAngle>>(
      m,
      "WMatToAxisAngle",
      "w_mat_to_axis_angle return type"
  )
      .def_readonly("axis", &SimUtils::WMatToAxisAngle::axis)
      .def_readonly("angle", &SimUtils::WMatToAxisAngle::angle)
      .def("__len__", [](const SimUtils::WMatToAxisAngle &) { return 2; })
      .def("__getitem__", [](const SimUtils::WMatToAxisAngle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.axis);
        if (i == 1)
          return py::cast(s.angle);
        throw py::index_error();
      });
  m.def(
      "w_mat_to_axis_angle",
      &SimUtils::w_mat_to_axis_angle,
      py::arg("w_mat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine w_mat_to_axis_angle (w_mat, axis, angle)

Routine to find the rotation axis and rotation angle corresponding to a given
3D rotation matrix.

The rotation angle is chosen in the range [0, pi].

Parameters
----------
w_mat : 2D array of float (shape: 3,3)
    Rotation matrix.

Returns
-------
axis : 1D array of float (shape: 3)
    Rotation axis. Normalized to 1.

angle : float
    Rotation angle in the range [0, pi].
)"""
  );
  m.def(
      "w_mat_to_quat",
      &SimUtils::w_mat_to_quat,
      py::arg("w_mat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function w_mat_to_quat (w_mat) result (quat)

Routine to find the quaternion corresponding to a given 3D rotation matrix.

Parameters
----------
w_mat : 2D array of float (shape: 3,3)
    Rotation matrix

Returns
-------
quat : 1D array of float (shape: 0:3)
    Quaternion.
)"""
  );
  m.def(
      "word_len",
      &SimUtils::word_len,
      py::arg("wording"),
      py::arg("wlen"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine word_len

Parameters
----------
wording : str

wlen : int
)"""
  );
  m.def(
      "word_read",
      &SimUtils::word_read,
      py::arg("in_str"),
      py::arg("delim_list"),
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("out_str"),
      py::arg("ignore_interior") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine word_read

Parameters
----------
in_str : str

delim_list : str

word : str

ix_word : int

delim : str

delim_found : bool

out_str : str

ignore_interior : bool, optional
)"""
  );
}
