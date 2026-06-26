#include "pybmad/generated/SimUtils_routines_q.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_q(nb::module_ &m) {
  m.def(
      "quadratic_roots",
      &SimUtils::quadratic_roots,
      nb::arg("coefs"),
      R"""(Wrapper for Fortran routine quadratic_roots

Parameters
----------
coefs : 1D array of float (shape: 3)
    Coefficients of the quadratic equation with 0 = coefs(1) + coefs(2) * x + coefs(3) * x^2

Returns
-------
root : 1D array of complex (shape: 2)
    Complex roots.
)"""
  );
  m.def(
      "quat_conj",
      nb::overload_cast<FixedArray1D<Complex, 4>>(&SimUtils::quat_conj),
      nb::arg("q_in"),
      R"""(Overloaded name to create the conjugate of a quaternian.
Overloaded functions are:
  Function quat_conj_real (q_in) result (q_out)
  Function quat_conj_complex (q_in) result (q_out)

Parameters
----------
q_in : 1D array of complex (shape: 0:3)
    Quaternion input.

Returns
-------
q_out : 1D array of complex (shape: 0:3)
    Conjugate quaternion.
)"""
  );
  m.def(
      "quat_conj",
      nb::overload_cast<FixedArray1D<Real, 4>>(&SimUtils::quat_conj),
      nb::arg("q_in"),
      R"""(Overloaded name to create the conjugate of a quaternian.
Overloaded functions are:
  Function quat_conj_real (q_in) result (q_out)
  Function quat_conj_complex (q_in) result (q_out)

Parameters
----------
q_in : 1D array of float (shape: 0:3)
    Quaternion input.

Returns
-------
q_out : 1D array of float (shape: 0:3)
    Conjugate quaternion.
)"""
  );
  m.def(
      "quat_inverse",
      &SimUtils::quat_inverse,
      nb::arg("q_in"),
      R"""(Routine to create the inverse of a quaternian.

Parameters
----------
q_in : 1D array of float (shape: 0:3)
    Quaternion input.

Returns
-------
q_out : 1D array of float (shape: 0:3)
    Inverse quaternion.
)"""
  );
  m.def(
      "quat_mul",
      nb::overload_cast<
          FixedArray1D<Complex, 4>,
          FixedArray1D<Complex, 4>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>,
          std::optional<FixedArray1D<Complex, 4>>>(&SimUtils::quat_mul),
      nb::arg("q1"),
      nb::arg("q2"),
      nb::arg("q3") = nb::none(),
      nb::arg("q4") = nb::none(),
      nb::arg("q5") = nb::none(),
      nb::arg("q6") = nb::none(),
      nb::arg("q7") = nb::none(),
      nb::arg("q8") = nb::none(),
      nb::arg("q9") = nb::none(),
      R"""(Overloaded name to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
Overloaded functions are:
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

Parameters
----------
q1 : 1D array of complex (shape: 0:3)
    Quaternions.

q2 : 1D array of complex (shape: 0:3)
    Quaternions.

q3 : 1D array of complex (shape: 0:3), optional
    More quaternions.

q9 : 1D array of complex (shape: 0:3), optional
    More quaternions.

Returns
-------
q_out : 1D array of complex (shape: 0:3)
    Resultant q1 * q2
)"""
  );
  m.def(
      "quat_mul",
      nb::overload_cast<
          FixedArray1D<Real, 4>,
          FixedArray1D<Real, 4>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>,
          std::optional<FixedArray1D<Real, 4>>>(&SimUtils::quat_mul),
      nb::arg("q1"),
      nb::arg("q2"),
      nb::arg("q3") = nb::none(),
      nb::arg("q4") = nb::none(),
      nb::arg("q5") = nb::none(),
      nb::arg("q6") = nb::none(),
      nb::arg("q7") = nb::none(),
      nb::arg("q8") = nb::none(),
      nb::arg("q9") = nb::none(),
      R"""(Overloaded name to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
Overloaded functions are:
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
  Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

Parameters
----------
q1 : 1D array of float (shape: 0:3)
    Quaternions.

q2 : 1D array of float (shape: 0:3)
    Quaternions.

q3 : 1D array of float (shape: 0:3), optional
    More quaternions.

q9 : 1D array of float (shape: 0:3), optional
    More quaternions.

Returns
-------
q_out : 1D array of float (shape: 0:3)
    Resultant q1 * q2
)"""
  );
  m.def(
      "quat_rotate",
      nb::overload_cast<FixedArray1D<Complex, 4>, FixedArray1D<Complex, 3>>(&SimUtils::quat_rotate),
      nb::arg("quat"),
      nb::arg("vec_in"),
      R"""(Overloaded name to rotate a vector using a quaternion..
Overloaded functions are:
  Function quat_rotate_real (quat, vec_in) result (vec_out)
  Function quat_rotate_complex (quat, vec_in) result (vec_out)

Parameters
----------
quat : 1D array of complex (shape: 0:3)
    Quaternion to rotate with. Does not have to be normalized.

vec_in : 1D array of complex (shape: 3)
    Initial vector.

Returns
-------
vec_out : 1D array of complex (shape: 3)
    Final vector.
)"""
  );
  m.def(
      "quat_rotate",
      nb::overload_cast<FixedArray1D<Real, 4>, FixedArray1D<Real, 3>>(&SimUtils::quat_rotate),
      nb::arg("quat"),
      nb::arg("vec_in"),
      R"""(Overloaded name to rotate a vector using a quaternion..
Overloaded functions are:
  Function quat_rotate_real (quat, vec_in) result (vec_out)
  Function quat_rotate_complex (quat, vec_in) result (vec_out)

Parameters
----------
quat : 1D array of float (shape: 0:3)
    Quaternion to rotate with. Does not have to be normalized.

vec_in : 1D array of float (shape: 3)
    Initial vector.

Returns
-------
vec_out : 1D array of float (shape: 3)
    Final vector.
)"""
  );
  nb::class_<SimUtils::QuatToAxisAngle>(m, "QuatToAxisAngle", "quat_to_axis_angle return type")
      .def_ro("axis", &SimUtils::QuatToAxisAngle::axis)
      .def_ro("angle", &SimUtils::QuatToAxisAngle::angle)
      .def("__len__", [](const SimUtils::QuatToAxisAngle &) { return 2; })
      .def("__getitem__", [](const SimUtils::QuatToAxisAngle &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.axis);
        if (i == 1)
          return nb::cast(s.angle);
        throw nb::index_error();
      });
  m.def(
      "quat_to_axis_angle",
      &SimUtils::quat_to_axis_angle,
      nb::arg("quat"),
      R"""(Routine to convert from quaternion to axis + angle representation.
The angle will be in the range 0 <= angle <= pi.

Parameters
----------
quat : 1D array of float (shape: 0:3)
    Rotation quaternion. Assumed normalized.

Returns
-------
axis : 1D array of float (shape: 3)
    Axis of rotation.

angle : float
    angle of rotation in range [0, pi].
)"""
  );
  m.def(
      "quat_to_omega",
      &SimUtils::quat_to_omega,
      nb::arg("quat"),
      R"""(Routine to convert rotation from quaternion representation to omega (axis + angle).

Parameters
----------
quat : 1D array of float (shape: 0:3)
    Rotation quaternion. Assumed normalized.

Returns
-------
omega : 1D array of float (shape: 3)
    Axis of rotation + magnitude = rotation angle.
)"""
  );
  m.def(
      "quat_to_w_mat",
      &SimUtils::quat_to_w_mat,
      nb::arg("quat"),
      R"""(Routine to construct the 3D rotation matrix w_mat given a rotation quaternion

Parameters
----------
quat : 1D array of float (shape: 0:3)
    Quaternion.

Returns
-------
w_mat : 2D array of float (shape: 3,3)
    Rotation matrix
)"""
  );
  m.def(
      "query_string",
      &SimUtils::query_string,
      nb::arg("query_str"),
      nb::arg("upcase"),
      nb::arg("return_str"),
      nb::arg("ix"),
      nb::arg("ios"),
      R"""(Wrapper for Fortran routine query_string

Parameters
----------
query_str : str

upcase : bool

return_str : str

ix : int

ios : int
)"""
  );
  m.def(
      "quote",
      &SimUtils::quote,
      nb::arg("str"),
      R"""(Wrapper for Fortran routine quote

Parameters
----------
str : str

Returns
-------
q_str : str
)"""
  );
  m.def(
      "quoten",
      &SimUtils::quoten,
      nb::arg("str"),
      nb::arg("delim") = nb::none(),
      R"""(Wrapper for Fortran routine quoten

Parameters
----------
str : 1D array of str

delim : str, optional

Returns
-------
q_str : str
)"""
  );
}
