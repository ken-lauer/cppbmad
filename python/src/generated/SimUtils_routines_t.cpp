#include "pybmad/generated/SimUtils_routines_t.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_t(nb::module_ &m) {
  m.def(
      "test_xgelbd",
      &SimUtils::test_xgelbd,
      R"""(Wrapper for Fortran routine test_xgelbd
)"""
  );
  m.def(
      "to_str",
      &SimUtils::to_str,
      nb::arg("num"),
      nb::arg("max_signif") = nb::none(),
      R"""(no longer exists
subroutine test_tune_tracker_lock (tracker_locked)
  implicit none
  logical tracker_locked(2)
end subroutine
)"""
  );
  nb::class_<SimUtils::TricubicCmplxEval>(m, "TricubicCmplxEval", "tricubic_cmplx_eval return type")
      .def_ro("df_dx", &SimUtils::TricubicCmplxEval::df_dx)
      .def_ro("df_dy", &SimUtils::TricubicCmplxEval::df_dy)
      .def_ro("df_dz", &SimUtils::TricubicCmplxEval::df_dz)
      .def_ro("f_val", &SimUtils::TricubicCmplxEval::f_val)
      .def("__len__", [](const SimUtils::TricubicCmplxEval &) { return 4; })
      .def("__getitem__", [](const SimUtils::TricubicCmplxEval &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.df_dx);
        if (i == 1)
          return nb::cast(s.df_dy);
        if (i == 2)
          return nb::cast(s.df_dz);
        if (i == 3)
          return nb::cast(s.f_val);
        throw nb::index_error();
      });
  m.def(
      "tricubic_cmplx_eval",
      &SimUtils::tricubic_cmplx_eval,
      nb::arg("x_norm"),
      nb::arg("y_norm"),
      nb::arg("z_norm"),
      nb::arg("tri_coef"),
      R"""(Routine to evaluate a tricubic interpolating complex function.

Use the routine tricubic_interpolation_cmplx_coefs to generate tri_coef.

Note: In the equations below, the eight points of the grid box being interpolated range
from (x0, y0, z0) to (x0+dx, y0+dy, z0+dz).

Parameters
----------
x_norm : float
    x_norm = (x - x0) / dx

y_norm : float
    y_norm = (y - y0) / dy

z_norm : float
    z_norm = (z - z0) / dz

tri_coef : TricubicCmplxCoefStruct
    Coefficients.

Returns
-------
f_val : complex
    Value of f.

df_dx : complex, optional
    Normalized first derivative: True df/dx = df_dx * dx

df_dy : complex, optional
    Normalized first derivative: True df/dy = df_dy * dy

df_dz : complex, optional
    Normalized first derivative: True df/dz = df_dz * dz
)"""
  );
  nb::class_<SimUtils::TricubicEval>(m, "TricubicEval", "tricubic_eval return type")
      .def_ro("df_dx", &SimUtils::TricubicEval::df_dx)
      .def_ro("df_dy", &SimUtils::TricubicEval::df_dy)
      .def_ro("df_dz", &SimUtils::TricubicEval::df_dz)
      .def_ro("f_val", &SimUtils::TricubicEval::f_val)
      .def("__len__", [](const SimUtils::TricubicEval &) { return 4; })
      .def("__getitem__", [](const SimUtils::TricubicEval &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.df_dx);
        if (i == 1)
          return nb::cast(s.df_dy);
        if (i == 2)
          return nb::cast(s.df_dz);
        if (i == 3)
          return nb::cast(s.f_val);
        throw nb::index_error();
      });
  m.def(
      "tricubic_eval",
      &SimUtils::tricubic_eval,
      nb::arg("x_norm"),
      nb::arg("y_norm"),
      nb::arg("z_norm"),
      nb::arg("tri_coef"),
      R"""(Routine to evaluate a tricubic interpolating function.

Use the routine tricubic_interpolation_coefs to generate tri_coef.

Note: In the equations below, the eight points of the grid box being interpolated range
from (x0, y0, z0) to (x0+dx, y0+dy, z0+dz).

Parameters
----------
x_norm : float
    x_norm = (x - x0) / dx

y_norm : float
    y_norm = (y - y0) / dy

z_norm : float
    z_norm = (z - z0) / dz

tri_coef : TricubicCoefStruct
    Coefficients.

Returns
-------
f_val : float
    Value of f.

df_dx : float, optional
    Normalized first derivative: True df/dx = df_dx * dx

df_dy : float, optional
    Normalized first derivative: True df/dy = df_dy * dy

df_dz : float, optional
    Normalized first derivative: True df/dz = df_dz * dz
)"""
  );
  m.def(
      "tricubic_interpolation_cmplx_coefs",
      &SimUtils::tricubic_interpolation_cmplx_coefs,
      nb::arg("field_at_box"),
      R"""(Routine to compute the complex coefficients for tricubic interpolation.

Use the routine tricubic_cmplx_eval to evaluate the interpolation function.

Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy

Parameters
----------
field_at_box : CmplxFieldAt3dBoxStruct
    Field and normalized derivatives at the 4 grid points.

Returns
-------
tri_coef : TricubicCmplxCoefStruct
    Coefficients.
)"""
  );
  m.def(
      "tricubic_interpolation_coefs",
      &SimUtils::tricubic_interpolation_coefs,
      nb::arg("field_at_box"),
      R"""(Routine to compute the coefficients for tricubic interpolation.

Use the routine tricubic_eval to evaluate the interpolation function.

Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy

Parameters
----------
field_at_box : FieldAt3dBoxStruct
    Field and normalized derivatives at the 4 grid points.

Returns
-------
tri_coef : TricubicCoefStruct
    Coefficients.
)"""
  );
  m.def(
      "type_this_file",
      &SimUtils::type_this_file,
      nb::arg("filename"),
      R"""(Wrapper for Fortran routine type_this_file

Parameters
----------
filename : str
)"""
  );
}
