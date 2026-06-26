#include "pybmad/generated/SimUtils_routines_w.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_w(nb::module_ &m) {
  nb::class_<SimUtils::WMatToAxisAngle>(m, "WMatToAxisAngle", "w_mat_to_axis_angle return type")
      .def_ro("axis", &SimUtils::WMatToAxisAngle::axis)
      .def_ro("angle", &SimUtils::WMatToAxisAngle::angle)
      .def("__len__", [](const SimUtils::WMatToAxisAngle &) { return 2; })
      .def("__getitem__", [](const SimUtils::WMatToAxisAngle &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.axis);
        if (i == 1)
          return nb::cast(s.angle);
        throw nb::index_error();
      });
  m.def(
      "w_mat_to_axis_angle",
      &SimUtils::w_mat_to_axis_angle,
      nb::arg("w_mat"),
      R"""(Routine to find the rotation axis and rotation angle corresponding to a given
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
      nb::arg("w_mat"),
      R"""(Routine to find the quaternion corresponding to a given 3D rotation matrix.

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
      nb::arg("wording"),
      R"""(Wrapper for Fortran routine word_len

Parameters
----------
wording : str

Returns
-------
wlen : int
)"""
  );
  m.def(
      "word_read",
      &SimUtils::word_read,
      nb::arg("in_str"),
      nb::arg("delim_list"),
      nb::arg("word"),
      nb::arg("ix_word"),
      nb::arg("delim"),
      nb::arg("delim_found"),
      nb::arg("out_str"),
      nb::arg("ignore_interior") = nb::none(),
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
