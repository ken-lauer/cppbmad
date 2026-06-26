#include "pybmad/generated/Bmad_routines_k.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Bmad_routines_k(nb::module_ &m) {
  m.def(
      "key_name_to_key_index",
      &Bmad::key_name_to_key_index,
      nb::arg("key_str"),
      nb::arg("abbrev_allowed") = nb::none(),
      R"""(Wrapper for Fortran routine key_name_to_key_index

Parameters
----------
key_str : str
    Name of the key. Result is case insensitive.

abbrev_allowed : bool, optional
    Abbreviations (eg: "quad") allowed? Default is False. At least 3 characters are needed (except for
    rfcavity elements) if True.

Returns
-------
key_index : int
    Index of the key. Set to -1 if key_name not recognized.
)"""
  );
  nb::class_<Bmad::KickVectorCalc>(m, "KickVectorCalc", "kick_vector_calc return type")
      .def_ro("dr_ds", &Bmad::KickVectorCalc::dr_ds)
      .def_ro("err", &Bmad::KickVectorCalc::err)
      .def("__len__", [](const Bmad::KickVectorCalc &) { return 2; })
      .def("__getitem__", [](const Bmad::KickVectorCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.dr_ds);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "kick_vector_calc",
      &Bmad::kick_vector_calc,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("s_body"),
      nb::arg("orbit"),
      nb::arg("print_err") = nb::none(),
      R"""(Subroutine to calculate the dr/ds "kick vector" where
    r = [x, p_x, y, p_y, z, p_z, t, spin_x,y,z]

Remember: In order to simplify the calculation, in the body of any element, P0 is taken to be
the P0 at the exit end of the element.

  dr(1)/ds = dx/ds = dx/dt * dt/ds
  where:
    dx/dt = v_x = p_x / (1 + p_z)
    dt/ds = (1 + g*x) / v_s
    g = 1/rho, rho = bending radius (nonzero only in a dipole)

  dr(2)/ds = dp_x/ds = dP_x/dt * dt/ds / P0 + g_x * P_z
  where:
    dP_x/dt = EM_Force_x
    g_x = bending in x-plane.

  dr(3)/ds = dy/ds = dy/dt * dt/ds
  where:
    dy/dt = v_x

  dr(4)/ds = dp_y/ds = dP_y/dt * ds/dt / P0 + g_y * P_z
  where:
    dP_y/dt = EM_Force_y
    g_y = bending in y-plane.

  NOTE: dr(5)/ds IS IGNORED WHEN CALCULATING Z. SEE TRANSFER_THIS_ORBIT ABOVE.
  dr(5)/ds = dz/ds = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * c_light * [t(ref) - t]
                   = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * vec(5) / beta
  where:
    dt/ds(ref) = 1 / beta(ref)

  dr(6)/ds = dp_z/ds = d(EM_Force dot v_hat) * dt/ds / P0
  where:
     v_hat = velocity normalized to 1.

  dr(7)/ds = dt/ds

  dr(8:10)/ds = Spin omega vector

  dr(11)/ds = dt_ref/ds

Parameters
----------
ele : EleStruct
    Element being tracked thorugh.

param : LatParamStruct
    Lattice parameters.

orbit : CoordStruct
    Position of particle.

Returns
-------
dr_ds : 1D array of float (shape: 11)
    Kick vector.

err : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "kill_complex_taylor",
      &Bmad::kill_complex_taylor,
      nb::arg("complex_taylor"),
      R"""(Subroutine to deallocate a Bmad complex_taylor map.

Parameters
----------
complex_taylor : 1D array of ComplexTaylorStruct
    complex_taylor to be deallocated. It is OK if complex_taylor has already been deallocated.
    This parameter is an input/output and is modified in-place.
    As an output, complex_taylor: deallocated complex_taylor structure.
)"""
  );
  m.def(
      "kill_ptc_layouts",
      &Bmad::kill_ptc_layouts,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine kill_ptc_layouts

Parameters
----------
lat : LatStruct
    Bmad lattice with associated layouts.
)"""
  );
  m.def(
      "kill_taylor",
      &Bmad::kill_taylor,
      nb::arg("bmad_taylor"),
      R"""(Wrapper for Fortran routine kill_taylor

Parameters
----------
bmad_taylor : 1D array of TaylorStruct
    Taylor to be deallocated.
    This parameter is an input/output and is modified in-place.
    As an output, bmad_taylor: deallocated Taylor structure.
)"""
  );
  m.def(
      "kind_name",
      &Bmad::kind_name,
      nb::arg("this_kind"),
      R"""(function to return the name of a PTC kind.

Parameters
----------
this_kind : int
    PTC kind

Returns
-------
kind_str : str
    String representation
)"""
  );
  nb::class_<Bmad::KnotInterpolate>(m, "KnotInterpolate", "knot_interpolate return type")
      .def_ro("err_flag", &Bmad::KnotInterpolate::err_flag)
      .def_ro("y_pt", &Bmad::KnotInterpolate::y_pt)
      .def("__len__", [](const Bmad::KnotInterpolate &) { return 2; })
      .def("__getitem__", [](const Bmad::KnotInterpolate &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.y_pt);
        throw nb::index_error();
      });
  m.def(
      "knot_interpolate",
      &Bmad::knot_interpolate,
      nb::arg("x_knot"),
      nb::arg("y_knot"),
      nb::arg("x_pt"),
      nb::arg("interpolation"),
      R"""(Wrapper for Fortran routine knot_interpolate

Parameters
----------
x_knot : 1D array of float
    Knot x-values.

y_knot : 1D array of float
    Knot y-values.

x_pt : float
    Point to evaluate at.

interpolation : int
    Interpolation type. cubic$ or linear$.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.

y_pt : float
    Interpolated y-value.
)"""
  );
  m.def(
      "knots_to_string",
      &Bmad::knots_to_string,
      nb::arg("x_knot"),
      nb::arg("y_knot"),
      R"""(Wrapper for Fortran routine knots_to_string

Parameters
----------
x_knot : 1D array of float

y_knot : 1D array of float

Returns
-------
str : str
)"""
  );
}
