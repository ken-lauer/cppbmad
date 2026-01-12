#include "pybmad/generated/Bmad_routines_k.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyKeyNameToKeyIndex python_key_name_to_key_index(
    std::string key_str,
    std::optional<bool> abbrev_allowed,
    int key_index) {
  Bmad::key_name_to_key_index(key_str, abbrev_allowed, key_index);
  auto py_result{PyKeyNameToKeyIndex{key_index}};
  return py_result;
}
PyKickVectorCalc python_kick_vector_calc(
    EleProxy& ele,
    LatParamProxy& param,
    double s_body,
    CoordProxy& orbit,
    std::optional<bool> print_err = std::nullopt) {
  auto _result = Bmad::kick_vector_calc(
      ele, param, s_body, orbit, make_opt_ref(print_err));
  auto py_result{PyKickVectorCalc{_result, print_err}};
  return py_result;
}
PyKnotInterpolate python_knot_interpolate(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    double x_pt,
    int interpolation,
    double y_pt) {
  auto _result =
      Bmad::knot_interpolate(x_knot, y_knot, x_pt, interpolation, y_pt);
  auto py_result{PyKnotInterpolate{_result, y_pt}};
  return py_result;
}
PyKnotsToString python_knots_to_string(
    RealAlloc1D& x_knot,
    RealAlloc1D& y_knot,
    std::string str) {
  Bmad::knots_to_string(x_knot, y_knot, str);
  auto py_result{PyKnotsToString{str}};
  return py_result;
}

void init_Bmad_routines_k(py::module& m) {
  py::class_<PyKeyNameToKeyIndex, std::unique_ptr<PyKeyNameToKeyIndex>>(
      m,
      "KeyNameToKeyIndex",
      "Fortran routine key_name_to_key_index return value")
      .def_readonly("key_index", &PyKeyNameToKeyIndex::key_index)
      .def("__len__", [](const PyKeyNameToKeyIndex&) { return 1; })
      .def(
          "__getitem__", [](const PyKeyNameToKeyIndex& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.key_index);
            throw py::index_error();
          });
  m.def(
      "key_name_to_key_index",
      &python_key_name_to_key_index,
      py::arg("key_str"),
      py::arg("abbrev_allowed") = py::none(),
      py::arg("key_index"),
      R"""(Parameters
  ----------
  key_str : unknown
      Name of the key. Result is case insensitive.
  abbrev_allowed : bool, optional
      Abbreviations (eg: "quad") allowed? Default is False. At least 3 characters are needed (except for
      rfcavity elements) if True.
  key_index : 
  )""");
  py::class_<PyKickVectorCalc, std::unique_ptr<PyKickVectorCalc>>(
      m, "KickVectorCalc", "Fortran routine kick_vector_calc return value")
      .def_readonly("dr_ds", &PyKickVectorCalc::dr_ds)
      .def_readonly("err", &PyKickVectorCalc::err)
      .def_readonly("print_err", &PyKickVectorCalc::print_err)
      .def("__len__", [](const PyKickVectorCalc&) { return 3; })
      .def("__getitem__", [](const PyKickVectorCalc& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.dr_ds);
        if (i == 1)
          return py::cast(s.err);
        if (i == 2)
          return py::cast(s.print_err);
        throw py::index_error();
      });
  m.def(
      "kick_vector_calc",
      &python_kick_vector_calc,
      py::arg("ele"),
      py::arg("param"),
      py::arg("s_body"),
      py::arg("orbit"),
      py::arg("print_err") = py::none(),
      R"""(Subroutine kick_vector_calc (ele, param, s_rel, orbit, dr_ds, field, err, print_err)

  Subroutine to calculate the dr/ds "kick vector" where
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
  s_rel : float
      Distance from the start of the element to the particle.
  orbit : CoordStruct
      Position of particle.
  local_ref_frame : !
      Logical, If True then take the input coordinates -- Logical, If True then take the input coordinates as
      being with respect to the frame of referene of the element.

  Returns
  -------
  dr_ds : float
      Kick vector.
  field : EmFieldStruct
      Local field.
  err : bool
      Set True if there is an error.

  Notes
  -----
  Remember: In order to simplify the calculation, in the body of any element, P0 is taken to be
  )""");
  m.def(
      "kill_complex_taylor",
      &Bmad::kill_complex_taylor,
      py::arg("complex_taylor"),
      R"""(Subroutine kill_complex_taylor (complex_taylor)

  Subroutine to deallocate a Bmad complex_taylor map.

  Parameters
  ----------
  complex_taylor : ComplexTaylorStruct
      complex_taylor to be deallocated. It is OK if complex_taylor has already been deallocated.
      This parameter is an input/output and is modified in-place. As an output: deallocated complex_taylor
      structure.
  )""");
  m.def(
      "kill_ptc_layouts",
      &Bmad::kill_ptc_layouts,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Bmad lattice with associated layouts.
  )""");
  m.def(
      "kill_taylor",
      &Bmad::kill_taylor,
      py::arg("bmad_taylor"),
      R"""(Parameters
  ----------
  bmad_taylor : TaylorStruct
      Taylor to be deallocated.
      This parameter is an input/output and is modified in-place. As an output: deallocated Taylor structure.
  )""");
  m.def(
      "kind_name",
      &Bmad::kind_name,
      py::arg("this_kind"),
      R"""(Function kind_name (this_kind) result (kind_str)

  function to return the name of a PTC kind.

  Parameters
  ----------
  this_kind : int
      PTC kind

  Returns
  -------
  kind_str : unknown
      String representation
  )""");
  py::class_<PyKnotInterpolate, std::unique_ptr<PyKnotInterpolate>>(
      m, "KnotInterpolate", "Fortran routine knot_interpolate return value")
      .def_readonly("err_flag", &PyKnotInterpolate::err_flag)
      .def_readonly("y_pt", &PyKnotInterpolate::y_pt)
      .def("__len__", [](const PyKnotInterpolate&) { return 2; })
      .def("__getitem__", [](const PyKnotInterpolate& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.y_pt);
        throw py::index_error();
      });
  m.def(
      "knot_interpolate",
      &python_knot_interpolate,
      py::arg("x_knot"),
      py::arg("y_knot"),
      py::arg("x_pt"),
      py::arg("interpolation"),
      py::arg("y_pt"),
      R"""(Parameters
  ----------
  x_knot : float
      Knot x-values.
  y_knot : float
      Knot y-values.
  x_pt : float
      Point to evaluate at.
  interpolation : int
      Interpolation type. cubic$ or linear$.
  err_flag : bool
      Set True if there is an error. False otherwise.
  y_pt : 
  )""");
  py::class_<PyKnotsToString, std::unique_ptr<PyKnotsToString>>(
      m, "KnotsToString", "Fortran routine knots_to_string return value")
      .def_readonly("str", &PyKnotsToString::str)
      .def("__len__", [](const PyKnotsToString&) { return 1; })
      .def("__getitem__", [](const PyKnotsToString& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.str);
        throw py::index_error();
      });
  m.def(
      "knots_to_string",
      &python_knots_to_string,
      py::arg("x_knot"),
      py::arg("y_knot"),
      py::arg("str"),
      R"""(Parameters
  ----------
  x_knot : 
  y_knot : 
  str : 
  )""");
}
