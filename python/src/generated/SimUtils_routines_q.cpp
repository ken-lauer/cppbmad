#include "pybmad/generated/SimUtils_routines_q.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_q(py::module &m) {
  m.def(
      "quadratic_roots",
      &SimUtils::quadratic_roots,
      py::arg("coefs"),
      py::arg("root"),
      R"""(Parameters
----------
coefs : float
    Coefficients of the quadratic equation with 0 = coefs(1) + coefs(2) * x + coefs(3) * x^2
root : 
)"""
  );
  m.def(
      "quat_conj",
      py::overload_cast<FixedArray1D<Complex, 4>>(&SimUtils::quat_conj),
      py::arg("q_in"),
      R"""(Function quat_conj (q_in) result (q_out)

Overloaded name to create the conjugate of a quaternian.
Overloaded functions are:
  Function quat_conj_real (q_in) result (q_out)
  Function quat_conj_complex (q_in) result (q_out)

Parameters
----------
q_in : float
    Quaternion input.

Returns
-------
q_out : float
    Conjugate quaternion.
)"""
  );
  m.def(
      "quat_conj",
      py::overload_cast<FixedArray1D<Real, 4>>(&SimUtils::quat_conj),
      py::arg("q_in"),
      R"""(Function quat_conj (q_in) result (q_out)

Overloaded name to create the conjugate of a quaternian.
Overloaded functions are:
  Function quat_conj_real (q_in) result (q_out)
  Function quat_conj_complex (q_in) result (q_out)

Parameters
----------
q_in : float
    Quaternion input.

Returns
-------
q_out : float
    Conjugate quaternion.
)"""
  );
  m.def(
      "quat_inverse",
      &SimUtils::quat_inverse,
      py::arg("q_in"),
      R"""(Function quat_inverse (q_in) result (q_out)

Routine to create the inverse of a quaternian.

Parameters
----------
q_in : float
    Quaternion input.

Returns
-------
q_out : float
    Inverse quaternion.
)"""
  );
  m.def(
      "quat_mul",
      py::overload_cast<
          FixedArray1D<Complex, 4>,
          FixedArray1D<Complex, 4>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>>(&SimUtils::quat_mul),
      py::arg("q1"),
      py::arg("q2"),
      py::arg("q3") = py::none(),
      py::arg("q4") = py::none(),
      py::arg("q5") = py::none(),
      py::arg("q6") = py::none(),
      py::arg("q7") = py::none(),
      py::arg("q8") = py::none(),
      py::arg("q9") = py::none(),
      R"""(Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

Overloaded name to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
Overloaded functions are:
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

Parameters
----------
q1 : float
    Quaternions.
q3 : float, optional
    More quaternions.

Returns
-------
q_out : float
    Resultant q1 * q2
)"""
  );
  m.def(
      "quat_mul",
      py::overload_cast<
          FixedArray1D<Real, 4>,
          FixedArray1D<Real, 4>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>>(&SimUtils::quat_mul),
      py::arg("q1"),
      py::arg("q2"),
      py::arg("q3") = py::none(),
      py::arg("q4") = py::none(),
      py::arg("q5") = py::none(),
      py::arg("q6") = py::none(),
      py::arg("q7") = py::none(),
      py::arg("q8") = py::none(),
      py::arg("q9") = py::none(),
      R"""(Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

Overloaded name to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
Overloaded functions are:
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

Parameters
----------
q1 : float
    Quaternions.
q3 : float, optional
    More quaternions.

Returns
-------
q_out : float
    Resultant q1 * q2
)"""
  );
  m.def(
      "quat_rotate",
      py::overload_cast<FixedArray1D<Complex, 4>, FixedArray1D<Complex, 3>>(&SimUtils::quat_rotate),
      py::arg("quat"),
      py::arg("vec_in"),
      R"""(Function quat_rotate (quat, vec_in) result (vec_out)

Overloaded name to rotate a vector using a quaternion..
Overloaded functions are:
  Function quat_rotate_real (quat, vec_in) result (vec_out)
  Function quat_rotate_complex (quat, vec_in) result (vec_out)

Parameters
----------
quat : float
    Quaternion to rotate with. Does not have to be normalized.
vec_in : float
    Initial vector.

Returns
-------
vec_out : float
    Final vector.
)"""
  );
  m.def(
      "quat_rotate",
      py::overload_cast<FixedArray1D<Real, 4>, FixedArray1D<Real, 3>>(&SimUtils::quat_rotate),
      py::arg("quat"),
      py::arg("vec_in"),
      R"""(Function quat_rotate (quat, vec_in) result (vec_out)

Overloaded name to rotate a vector using a quaternion..
Overloaded functions are:
  Function quat_rotate_real (quat, vec_in) result (vec_out)
  Function quat_rotate_complex (quat, vec_in) result (vec_out)

Parameters
----------
quat : float
    Quaternion to rotate with. Does not have to be normalized.
vec_in : float
    Initial vector.

Returns
-------
vec_out : float
    Final vector.
)"""
  );
  py::class_<SimUtils::QuatToAxisAngle, std::unique_ptr<SimUtils::QuatToAxisAngle>>(
      m,
      "QuatToAxisAngle",
      "quat_to_axis_angle return type"
  )
      .def_readonly("axis", &SimUtils::QuatToAxisAngle::axis)
      .def_readonly("angle", &SimUtils::QuatToAxisAngle::angle)
      .def("__len__", [](const SimUtils::QuatToAxisAngle &) { return 2; })
      .def("__getitem__", [](const SimUtils::QuatToAxisAngle &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.axis);
        if (i == 1)
          return py::cast(s.angle);
        throw py::index_error();
      });
  m.def(
      "quat_to_axis_angle",
      &SimUtils::quat_to_axis_angle,
      py::arg("quat"),
      R"""(Subroutine quat_to_axis_angle (quat, axis, angle)

Routine to convert from quaternion to axis + angle representation.
The angle will be in the range 0 <= angle <= pi.

Parameters
----------
quat : float
    Rotation quaternion. Assumed normalized.

Returns
-------
axis : float
    Axis of rotation.
angle : float
    angle of rotation in range [0, pi].
)"""
  );
  m.def(
      "quat_to_omega",
      &SimUtils::quat_to_omega,
      py::arg("quat"),
      R"""(Function quat_to_omega (quat) result (omega)

Routine to convert rotation from quaternion representation to omega (axis + angle).

Parameters
----------
quat : float
    Rotation quaternion. Assumed normalized.

Returns
-------
omega : float
    Axis of rotation + magnitude = rotation angle.
)"""
  );
  m.def(
      "quat_to_w_mat",
      &SimUtils::quat_to_w_mat,
      py::arg("quat"),
      R"""(Function quat_to_w_mat (quat) result (w_mat)

Routine to construct the 3D rotation matrix w_mat given a rotation quaternion

Parameters
----------
quat : float
    Quaternion.

Returns
-------
w_mat : float
    Rotation matrix
)"""
  );
  m.def(
      "query_string",
      &SimUtils::query_string,
      py::arg("query_str"),
      py::arg("upcase"),
      py::arg("return_str"),
      py::arg("ix"),
      py::arg("ios"),
      R"""(Parameters
----------
query_str : 
upcase : 
return_str : 
ix : 
ios : 
)"""
  );
  m.def(
      "quote",
      &SimUtils::quote,
      py::arg("str"),
      py::arg("q_str"),
      R"""(Parameters
----------
str : 
q_str : 
)"""
  );
}
