#include "pybmad/generated/SimUtils_routines_s.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_s(py::module &m) {
  py::class_<SimUtils::Serbd, std::unique_ptr<SimUtils::Serbd>>(m, "Serbd", "serbd return type")
      .def_readonly("b", &SimUtils::Serbd::b)
      .def_readonly("d", &SimUtils::Serbd::d)
      .def("__len__", [](const SimUtils::Serbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Serbd &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.b);
        if (i == 1)
          return py::cast(s.d);
        throw py::index_error();
      });
  m.def(
      "serbd",
      &SimUtils::serbd,
      py::arg("y"),
      py::arg("m"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine serbd

Parameters
----------
y : float

m : float

Returns
-------
b : float

d : float
)"""
  );
  m.def(
      "set_env",
      &SimUtils::set_env,
      py::arg("env_name"),
      py::arg("env_value"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_env

Parameters
----------
env_name : str

env_value : str

err_flag : bool
)"""
  );
  m.def(
      "set_parameter",
      py::overload_cast<int, int, int>(&SimUtils::set_parameter),
      py::arg("param_val"),
      py::arg("set_val"),
      py::arg("save_val"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_parameter_int

Parameters
----------
param_val : int

set_val : int

save_val : int
)"""
  );
  m.def(
      "set_parameter",
      py::overload_cast<bool, bool, bool>(&SimUtils::set_parameter),
      py::arg("param_val"),
      py::arg("set_val"),
      py::arg("save_val"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_parameter_logic

Parameters
----------
param_val : bool

set_val : bool

save_val : bool
)"""
  );
  m.def(
      "set_parameter",
      py::overload_cast<double, double, double>(&SimUtils::set_parameter),
      py::arg("param_val"),
      py::arg("set_val"),
      py::arg("save_val"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine set_parameter_real

Parameters
----------
param_val : float

set_val : float

save_val : float
)"""
  );
  m.def(
      "set_species_charge",
      &SimUtils::set_species_charge,
      py::arg("species_in"),
      py::arg("charge"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function set_species_charge(species_in, charge) result(species_charged)

Routine to return the ID for a particle of the same type as species_in but with a different charge.
Exception: If species_in corresponds to a subatomic particle, the charge argument is ignored and
species_charged will be set equal to species_in.

Parameters
----------
species_in : int
    Input species.

charge : int
    Charge to set species_charged to.

Returns
-------
species_charged : int
    Species of the same type as species_in but with different charge.
)"""
  );
  m.def(
      "sign_of",
      py::overload_cast<int, std::optional<bool>>(&SimUtils::sign_of),
      py::arg("num"),
      py::arg("zero_is_zero") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function sign_of (num, zero_is_zero) result (num_sign)

Routine to return the sign of a number.
Note: Fortran instrinsic sign function is similar to sign_of with zero_is_zero = False.

Parameters
----------
num : int
    Input number

zero_is_zero : bool, optional
    If True (default), num = 0 gives num_sign = 0. If False, num = 0 gives num_sign = 1.

Returns
-------
num_sign : int
    +1 if num is positive, -1 if num is negative, and 0 or +1 if num is zero depending upon setting of
    zero_is_zero.
)"""
  );
  m.def(
      "sign_of",
      py::overload_cast<double, std::optional<bool>>(&SimUtils::sign_of),
      py::arg("num"),
      py::arg("zero_is_zero") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function sign_of (num, zero_is_zero) result (num_sign)

Routine to return the sign of a number.
Note: Fortran instrinsic sign function is similar to sign_of with zero_is_zero = False.

Parameters
----------
num : float
    Input number

zero_is_zero : bool, optional
    If True (default), num = 0 gives num_sign = 0. If False, num = 0 gives num_sign = 1.

Returns
-------
num_sign : int
    +1 if num is positive, -1 if num is negative, and 0 or +1 if num is zero depending upon setting of
    zero_is_zero.
)"""
  );
  m.def(
      "sinc",
      &SimUtils::sinc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sinc

Parameters
----------
x : float
    Number.

nd : int, optional
    Derivative order. nd = 0 (default) -> compute sin(x) / x NOTE: Currently only nd = 0 and nd = 1 are
    implemented.

Returns
-------
y : float
    nd^th derivative of sin(x) / x
)"""
  );
  m.def(
      "sincc",
      &SimUtils::sincc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sincc

Parameters
----------
x : float
    Number.

nd : int, optional
    Derivative order. nd = 0 (default) -> compute (x - sin(x)) / x^3 NOTE: Currently only nd = 0 and nd = 1
    are implemented.

Returns
-------
y : float
    nd^th derivative of (x - sin(x)) / x^3
)"""
  );
  m.def(
      "sinhx_x",
      &SimUtils::sinhx_x,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sinhx_x

Parameters
----------
x : float
    Number.

nd : int, optional
    Derivative order. nd = 0 (default) -> compute sinh(x) / x NOTE: Currently only nd = 0 and nd = 1 are
    implemented.

Returns
-------
y : float
    nd^th derivative of sinh(x) / x.
)"""
  );
  m.def(
      "skip_header",
      &SimUtils::skip_header,
      py::arg("ix_unit"),
      py::arg("error_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine skip_header

Parameters
----------
ix_unit : int

error_flag : bool
)"""
  );
  m.def(
      "special_projection",
      &SimUtils::special_projection,
      py::arg("f"),
      py::arg("func_retval__"),
      py::arg("status") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(function special_projection

Calculates <cdata | exp(i theta)>

Used only by maximize projection. Uses data global to the function to accomodate stock NR routine.
)"""
  );
  m.def(
      "species_id",
      &SimUtils::species_id,
      py::arg("name"),
      py::arg("default_") = py::none(),
      py::arg("print_err") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function species_id (name, default, print_err) result(species)

Routine to return the integer ID index of a particle species given the name.

For subatomic particles, the case does not matter.
For all other types of particles, the case does matter.

Parameters
----------
name : str
    Name of the species.

print_err : bool, optional
    Print error message? Default is True. If False, return species = invalid$,

Returns
-------
species : int
    Species ID. Will return invalid$ if name is not valid. Will return not_set$ if name is blank
)"""
  );
  m.def(
      "species_id_from_openpmd",
      &SimUtils::species_id_from_openpmd,
      py::arg("pmd_name"),
      py::arg("charge"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function species_id_from_openpmd (pmd_name, charge) result(species)

Routine to return the Bmad species ID given the openPMD species name and given particle charge.
Note: If pmd_name corresponds to a subatomic particle, the charge argument is ignored.

Parameters
----------
pmd_name : str
    OpenPMD species name.

charge : int
    Species charge. Ignored for subatomic particles.

Returns
-------
species : int
    Bmad spicies ID number.
)"""
  );
  m.def(
      "species_name",
      &SimUtils::species_name,
      py::arg("species"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function species_name (species) result(name)

Routine to return the name of a particle species given the integer index.

Parameters
----------
species : int
    Species ID.

Returns
-------
name : str
    Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
)"""
  );
  m.def(
      "species_of",
      &SimUtils::species_of,
      py::arg("mass"),
      py::arg("charge"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function species_of (mass, charge) result (species)

Routine to return the integer ID index of a particle species given the mass and charge.
Note: Currently this routine only works for subatomic particles and is used for decoding PTC flat files.

Parameters
----------
mass : float
    Mass of the particle

charge : int
    Charge of the particle.

Returns
-------
species : int
    Species ID. Will return invalid$ if name is not valid.
)"""
  );
  m.def(
      "spin_of",
      &SimUtils::spin_of,
      py::arg("species"),
      py::arg("non_subatomic_default") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function spin_of (species, non_subatomic_default) result (spin)

Routine to return the spin, in units of hbar, of a particle.
This routine is only valid for subatomic particles.
For all other particles, the returned spin value will be the value of non_subatomic_default.

Parameters
----------
species : int
    Species ID.

non_subatomic_default : float, optional
    Default value to be used for non-subatomic species. Default value of this argument is zero.

Returns
-------
spin : float
    Particle spin.
)"""
  );
  m.def(
      "spline1",
      &SimUtils::spline1,
      py::arg("a_spline"),
      py::arg("x"),
      py::arg("n") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function spline1 (a_spline, x, n) result (y)

Function for spline evaluation using a single spline (instead of a spline array).
Also see:
  spline_evaluate
  spline_akima_interpolate

Modules used:
  use spline_mod

Parameters
----------
a_spline : SplineStruct
    Single spline structure.

x : float
    Point for evaluation.

n : int, optional
    Output derivative order. May be -1, 0, 1, 2, or 3. Default is 0. n = -1 => output is integral of y from
    a_spline.x0 to x. n = 1 => output is dy/dx, n = 2 => output is d^2y/dx^2, etc.

Returns
-------
y : float
    Interpolated spline value or derivative.
)"""
  );
  m.def(
      "spline_akima",
      &SimUtils::spline_akima,
      py::arg("spline"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine spline_akima (spline, ok)

Given a set of (x,y) points we want to interpolate between the points.
This subroutine computes the semi-hermite cubic spline developed by
Hiroshi Akima. The spline goes thorugh all the points (that is, it is
not a smoothing spline). For interpolation use:
  spline_evaluate
  spline_akima_interpolate ! You do not need to call spline_akima if you use this routine.

Reference:
  H Akima, "A New Method of Interpolation and Smooth Curve Fitting Based
  on Local Procedures", J. Assoc. Comp. Mach., Vol 17(4), 589-602 (1970).

Modules used:
  use spline_mod

Parameters
----------
spline : 1D array of SplineStruct

Returns
-------
ok : bool
    Set .false. if something is wrong (like less than 2 points used).
)"""
  );
  py::class_<SimUtils::SplineAkimaInterpolate, std::unique_ptr<SimUtils::SplineAkimaInterpolate>>(
      m,
      "SplineAkimaInterpolate",
      "spline_akima_interpolate return type"
  )
      .def_readonly("ok", &SimUtils::SplineAkimaInterpolate::ok)
      .def_readonly("y", &SimUtils::SplineAkimaInterpolate::y)
      .def_readonly("dy", &SimUtils::SplineAkimaInterpolate::dy)
      .def("__len__", [](const SimUtils::SplineAkimaInterpolate &) { return 3; })
      .def("__getitem__", [](const SimUtils::SplineAkimaInterpolate &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ok);
        if (i == 1)
          return py::cast(s.y);
        if (i == 2)
          return py::cast(s.dy);
        throw py::index_error();
      });
  m.def(
      "spline_akima_interpolate",
      &SimUtils::spline_akima_interpolate,
      py::arg("x_knot"),
      py::arg("y_knot"),
      py::arg("x"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine spline_akima_interpolate (x_knot, y_knot, x, ok, y, dy)

Routine to interpolate using an akima spline.

When evaluating at enough points, this routine is slower than calling spline_akima to
first evaluate the spline coefficients and then repeatedly calling spline_evaluate.

The advantage of this routine is that only the (x, y) knot points need to be stored
and it will be faster if the number of evaluations is small.

This routine will extrapolate past the range of x_knot(:) up to a distance equal to the
length between an end point and the point just inside the end point.

Parameters
----------
x_knot : 1D array of float
    Array of x values for the knot points. Must have more than 2 points and be in asending order.

y_knot : 1D array of float
    Array of y values for the knot points. Must be same size as x_knot(:).

x : float
    Point to evaluate at.

Returns
-------
ok : bool
    Set .true. if everything ok, That is, x is within the spline range.

y : float, optional
    Spline interpolation.

dy : float, optional
    Spline derivative interpolation.
)"""
  );
  py::class_<SimUtils::SplineEvaluate, std::unique_ptr<SimUtils::SplineEvaluate>>(
      m,
      "SplineEvaluate",
      "spline_evaluate return type"
  )
      .def_readonly("ok", &SimUtils::SplineEvaluate::ok)
      .def_readonly("y", &SimUtils::SplineEvaluate::y)
      .def_readonly("dy", &SimUtils::SplineEvaluate::dy)
      .def("__len__", [](const SimUtils::SplineEvaluate &) { return 3; })
      .def("__getitem__", [](const SimUtils::SplineEvaluate &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ok);
        if (i == 1)
          return py::cast(s.y);
        if (i == 2)
          return py::cast(s.dy);
        throw py::index_error();
      });
  m.def(
      "spline_evaluate",
      &SimUtils::spline_evaluate,
      py::arg("spline"),
      py::arg("x"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine spline_evaluate (spline, x, ok, y, dy)

Subroutine to evalueate a spline at a set of points.

A point outside of the range of knot points is an error.
Also see:
  spline1
  spline_akima_interpolate

A spline may be generated using, for example, the spline_akima routine.

Modules used:
  use spline_mod

Parameters
----------
spline : 1D array of SplineStruct
    Spline structure.

x : float
    point for evaluation.

Returns
-------
ok : bool
    Set .true. if everything ok. That is, x is within the spline range.

y : float, optional
    Spline interpolation.

dy : float, optional
    Spline derivative interpolation.
)"""
  );
  m.def(
      "sqrt_alpha",
      &SimUtils::sqrt_alpha,
      py::arg("alpha"),
      py::arg("x"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sqrt_alpha

Parameters
----------
alpha : float
    Number

x : float
    Number

Returns
-------
y : float
    Result.
)"""
  );
  m.def(
      "sqrt_one",
      &SimUtils::sqrt_one,
      py::arg("x"),
      py::arg("ds1"),
      py::arg("nd") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine sqrt_one

Parameters
----------
x : float
    Number

ds1 : float

nd : int, optional
    Derivative order. nd = 0 (default) -> compute Sqrt[1+x] - 1. NOTE: Currently only nd = 0 and nd = 1 are
    implemented.
)"""
  );
  m.def(
      "str_count",
      &SimUtils::str_count,
      py::arg("str"),
      py::arg("match_"),
      py::arg("num"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_count

Parameters
----------
str : str

match : str

num : int
)"""
  );
  m.def(
      "str_downcase",
      &SimUtils::str_downcase,
      py::arg("src"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_downcase

Parameters
----------
src : str

Returns
-------
dst : str
)"""
  );
  m.def(
      "str_first_in_set",
      &SimUtils::str_first_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      py::arg("ignore_clauses") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_first_in_set

Parameters
----------
line : str

set : str

ix_match : int

ignore_clauses : bool, optional
)"""
  );
  m.def(
      "str_first_not_in_set",
      &SimUtils::str_first_not_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_first_not_in_set

Parameters
----------
line : str

set : str

ix_match : int
)"""
  );
  m.def(
      "str_last_in_set",
      &SimUtils::str_last_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_last_in_set

Parameters
----------
line : str

set : str

ix_match : int
)"""
  );
  m.def(
      "str_last_not_in_set",
      &SimUtils::str_last_not_in_set,
      py::arg("line"),
      py::arg("set"),
      py::arg("ix_match"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_last_not_in_set

Parameters
----------
line : str

set : str

ix_match : int
)"""
  );
  m.def(
      "str_match_wild",
      &SimUtils::str_match_wild,
      py::arg("str"),
      py::arg("pat"),
      py::arg("a_match"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_match_wild

Parameters
----------
str : str

pat : str

a_match : bool
)"""
  );
  m.def(
      "str_substitute",
      &SimUtils::str_substitute,
      py::arg("string"),
      py::arg("str_match") = py::none(),
      py::arg("str_replace") = py::none(),
      py::arg("do_trim") = py::none(),
      py::arg("ignore_escaped") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_substitute

Parameters
----------
string : str

str_match : str, optional

str_replace : str, optional

do_trim : bool, optional

ignore_escaped : bool, optional
)"""
  );
  m.def(
      "str_upcase",
      &SimUtils::str_upcase,
      py::arg("src"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine str_upcase

Parameters
----------
src : str

Returns
-------
dst : str
)"""
  );
  m.def(
      "string_to_int",
      &SimUtils::string_to_int,
      py::arg("line"),
      py::arg("default_"),
      py::arg("err_flag"),
      py::arg("value"),
      py::arg("err_print_flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine string_to_int

Parameters
----------
line : str

err_flag : bool

value : int

err_print_flag : bool, optional
)"""
  );
  m.def(
      "string_to_real",
      &SimUtils::string_to_real,
      py::arg("line"),
      py::arg("default_"),
      py::arg("err_flag"),
      py::arg("value"),
      py::arg("err_print_flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine string_to_real

Parameters
----------
line : str

err_flag : bool

value : float

err_print_flag : bool, optional
)"""
  );
  m.def(
      "string_trim",
      &SimUtils::string_trim,
      py::arg("in_string"),
      py::arg("out_string"),
      py::arg("word_len"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine string_trim

Parameters
----------
in_string : str

out_string : str

word_len : int
)"""
  );
  m.def(
      "string_trim2",
      &SimUtils::string_trim2,
      py::arg("in_str"),
      py::arg("delimitors"),
      py::arg("out_str"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("ix_next"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine string_trim2

Parameters
----------
in_str : str

delimitors : str

out_str : str

ix_word : int

delim : str

ix_next : int
)"""
  );
  m.def(
      "suggest_lmdif",
      &SimUtils::suggest_lmdif,
      py::arg("XV"),
      py::arg("FV"),
      py::arg("eps"),
      py::arg("itermx"),
      py::arg("reset_flag") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(subroutine suggest_lmdif (xv, fv, eps, itermx, at_end, reset_flag)

Reverse communication subroutine.
It suggests values for your input variables based on
the previous value of your merit function.

Use initial_lmdif to initialize internal variables

Parameters
----------
xv : 1D array of float
    Array of variables
    This parameter is an input/output and is modified in-place.
    As an output, xv: Suggested new values

fv : 1D array of float
    Array of function value/s that should be optimized to zero
    This parameter is an input/output and is modified in-place.
    As an output, fv: After the last optimization this returns the best values ever.

eps : float
    Desired accuracy with which the optimum should be found.

itermx : int
    Max number of iterations

reset_flag : bool, optional
    Optional. Used by initial_lmdif to clear previous saved values

Returns
-------
at_end : bool
    Set to False if more optimization is recommended. If set to True then xv(:) will be the minimum found.
)"""
  );
  m.def(
      "super_bicubic_coef",
      &SimUtils::super_bicubic_coef,
      py::arg("y"),
      py::arg("y1"),
      py::arg("y2"),
      py::arg("y12"),
      py::arg("d1"),
      py::arg("d2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine super_bicubic_coef(y, y1, y2, y12, d1, d2, c)

Routine to compute coefficients for bicubic interpolation.
This is from NR bcucof.

Parameters
----------
y : 1D array of float (shape: 4)
    Function values at grid points.

y1 : 1D array of float (shape: 4)
    dy/dx1 derivatives.

y2 : 1D array of float (shape: 4)
    dy/dx2 derivatives.

y12 : 1D array of float (shape: 4)
    d2y/dx1*dx2 second derivatives.

d1 : float
    Grid width in 1-direction.

d2 : float
    Grid width in 2-direction.

Returns
-------
c : 2D array of float (shape: 4,4)
    Coefficients.
)"""
  );
  py::class_<
      SimUtils::SuperBicubicInterpolation,
      std::unique_ptr<SimUtils::SuperBicubicInterpolation>>(
      m,
      "SuperBicubicInterpolation",
      "super_bicubic_interpolation return type"
  )
      .def_readonly("ansy", &SimUtils::SuperBicubicInterpolation::ansy)
      .def_readonly("ansy1", &SimUtils::SuperBicubicInterpolation::ansy1)
      .def_readonly("ansy2", &SimUtils::SuperBicubicInterpolation::ansy2)
      .def("__len__", [](const SimUtils::SuperBicubicInterpolation &) { return 3; })
      .def("__getitem__", [](const SimUtils::SuperBicubicInterpolation &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.ansy);
        if (i == 1)
          return py::cast(s.ansy1);
        if (i == 2)
          return py::cast(s.ansy2);
        throw py::index_error();
      });
  m.def(
      "super_bicubic_interpolation",
      &SimUtils::super_bicubic_interpolation,
      py::arg("y"),
      py::arg("y1"),
      py::arg("y2"),
      py::arg("y12"),
      py::arg("x1l"),
      py::arg("x1u"),
      py::arg("x2l"),
      py::arg("x2u"),
      py::arg("x1"),
      py::arg("x2"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine super_bicubic_interpolation(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, ansy, ansy1, ansy2)

Routine to do bicubic interpolation.
This is from NR bcuint.

Note! The four grid points are arrayed in counter-clockwise order beginning from the lower left.
So, for example, y = [y_ll, y_lu, y_uu, y_ul] where "l" = lower, "u" = upper index.

Parameters
----------
y : 1D array of float (shape: 4)
    Function values at grid points.

y1 : 1D array of float (shape: 4)
    dy/dx1 derivatives.

y2 : 1D array of float (shape: 4)
    dy/dx2 derivatives.

y12 : 1D array of float (shape: 4)
    d2y/dx1*dx2 second derivatives.

x1l : float
    1-direction coordinate at lower points.

x1u : float
    1-direction coordinate at upper points

x2l : float
    2-direction coordinate at lower points.

x2u : float
    2-direction coordinate at upper points

x1 : float
    1-direction coordinate at point to evaluate.

x2 : float
    2-direction coordinate at point to evaluate.

Returns
-------
ansy : float
    Interpolation value.

ansy1 : float
    1-direction derivative at interpolation point.

ansy2 : float
    2-direction derivative at interpolation point.
)"""
  );
  py::class_<SimUtils::SuperPolint, std::unique_ptr<SimUtils::SuperPolint>>(
      m,
      "SuperPolint",
      "super_polint return type"
  )
      .def_readonly("y", &SimUtils::SuperPolint::y)
      .def_readonly("dy", &SimUtils::SuperPolint::dy)
      .def("__len__", [](const SimUtils::SuperPolint &) { return 2; })
      .def("__getitem__", [](const SimUtils::SuperPolint &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.y);
        if (i == 1)
          return py::cast(s.dy);
        throw py::index_error();
      });
  m.def(
      "super_polint",
      &SimUtils::super_polint,
      py::arg("xa"),
      py::arg("ya"),
      py::arg("x"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function super_polint (xa, ya, x, y, dy)

This is essentially polint from Numerical Recipes.

Parameters
----------
xa : 1D array of float

ya : 1D array of float

x : float

Returns
-------
y : float

dy : float
)"""
  );
  m.def(
      "super_poly",
      &SimUtils::super_poly,
      py::arg("x"),
      py::arg("coeffs"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function super_poly (x, coef) result (value)

Routine to compute Sum: coef(i)*x^i

Parameters
----------
x : float
    Variable.

Returns
-------
value : float
    Polynomial value.
)"""
  );
  m.def(
      "super_sobseq",
      &SimUtils::super_sobseq,
      py::arg("x"),
      py::arg("ran_state") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine super_sobseq (x, ran_state)

Routine patterened after sobseq in Numerical Recipes.
Difference is that this version has an argument for the internal state.

Parameters
----------
x : 1D array of float
    Random vector.

ran_state : RandomStateStruct, optional
    Generator state. See the ran_seed_put documentation for more details.
)"""
  );
  m.def(
      "super_sort",
      &SimUtils::super_sort,
      py::arg("arr"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine super_sort(arr)

Routine to sort an integer array in place.
This is the NR routine sort modified to sort integers.

Parameters
----------
arr : 1D array of int
    Array of integers.
    This parameter is an input/output and is modified in-place.
    As an output, arr: Sorted array.
)"""
  );
  m.def(
      "system_command",
      &SimUtils::system_command,
      py::arg("line"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine system_command

Parameters
----------
line : str
    Command to pass to the system shell.

Returns
-------
err_flag : bool, optional
    Set True if there is an error (bad command, etc.).
)"""
  );
}
