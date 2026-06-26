#include "pybmad/generated/SimUtils_routines_s.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_s(nb::module_ &m) {
  nb::class_<SimUtils::Serbd>(m, "Serbd", "serbd return type")
      .def_ro("b", &SimUtils::Serbd::b)
      .def_ro("d", &SimUtils::Serbd::d)
      .def("__len__", [](const SimUtils::Serbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Serbd &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.b);
        if (i == 1)
          return nb::cast(s.d);
        throw nb::index_error();
      });
  m.def(
      "serbd",
      &SimUtils::serbd,
      nb::arg("y"),
      nb::arg("m"),
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
      "set_all_ptr",
      &SimUtils::set_all_ptr,
      nb::arg("a_ptr"),
      nb::arg("value"),
      nb::arg("delta") = nb::none(),
      nb::arg("value_set") = nb::none(),
      R"""(Wrapper for Fortran routine set_all_ptr

Parameters
----------
a_ptr : AllPointerStruct

value : float

delta : bool, optional

value_set : float, optional
)"""
  );
  m.def(
      "set_env",
      &SimUtils::set_env,
      nb::arg("env_name"),
      nb::arg("env_value"),
      nb::arg("err_flag"),
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
      nb::overload_cast<int, int>(&SimUtils::set_parameter),
      nb::arg("param_val"),
      nb::arg("set_val"),
      R"""(Wrapper for Fortran routine set_parameter_int

Parameters
----------
param_val : int

set_val : int

Returns
-------
save_val : int
)"""
  );
  m.def(
      "set_parameter",
      nb::overload_cast<bool, bool>(&SimUtils::set_parameter),
      nb::arg("param_val"),
      nb::arg("set_val"),
      R"""(Wrapper for Fortran routine set_parameter_logic

Parameters
----------
param_val : bool

set_val : bool

Returns
-------
save_val : bool
)"""
  );
  m.def(
      "set_parameter",
      nb::overload_cast<double, double>(&SimUtils::set_parameter),
      nb::arg("param_val"),
      nb::arg("set_val"),
      R"""(Wrapper for Fortran routine set_parameter_real

Parameters
----------
param_val : float

set_val : float

Returns
-------
save_val : float
)"""
  );
  m.def(
      "set_species_charge",
      &SimUtils::set_species_charge,
      nb::arg("species_in"),
      nb::arg("charge"),
      R"""(Routine to return the ID for a particle of the same type as species_in but with a different charge.
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
      nb::overload_cast<int, std::optional<bool>>(&SimUtils::sign_of),
      nb::arg("num"),
      nb::arg("zero_is_zero") = nb::none(),
      R"""(Routine to return the sign of a number.
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
      nb::overload_cast<double, std::optional<bool>>(&SimUtils::sign_of),
      nb::arg("num"),
      nb::arg("zero_is_zero") = nb::none(),
      R"""(Routine to return the sign of a number.
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
      nb::arg("x"),
      nb::arg("nd") = nb::none(),
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
      nb::arg("x"),
      nb::arg("nd") = nb::none(),
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
      nb::arg("x"),
      nb::arg("nd") = nb::none(),
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
      nb::arg("ix_unit"),
      nb::arg("error_flag"),
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
      nb::arg("f"),
      nb::arg("status") = nb::none(),
      R"""(Calculates <cdata | exp(i theta)>

Used only by maximize projection. Uses data global to the function to accomodate stock NR routine.
)"""
  );
  m.def(
      "species_id",
      &SimUtils::species_id,
      nb::arg("name"),
      nb::arg("default_") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to return the integer ID index of a particle species given the name.

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
      nb::arg("pmd_name"),
      nb::arg("charge"),
      R"""(Routine to return the Bmad species ID given the openPMD species name and given particle charge.
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
      nb::arg("species"),
      R"""(Routine to return the name of a particle species given the integer index.

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
      nb::arg("mass"),
      nb::arg("charge"),
      R"""(Routine to return the integer ID index of a particle species given the mass and charge.
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
      nb::arg("species"),
      nb::arg("non_subatomic_default") = nb::none(),
      R"""(Routine to return the spin, in units of hbar, of a particle.
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
      nb::arg("a_spline"),
      nb::arg("x"),
      nb::arg("n") = nb::none(),
      R"""(Function for spline evaluation using a single spline (instead of a spline array).
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
      nb::arg("spline"),
      R"""(Given a set of (x,y) points we want to interpolate between the points.
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
  nb::class_<SimUtils::SplineAkimaInterpolate>(
      m,
      "SplineAkimaInterpolate",
      "spline_akima_interpolate return type"
  )
      .def_ro("ok", &SimUtils::SplineAkimaInterpolate::ok)
      .def_ro("y", &SimUtils::SplineAkimaInterpolate::y)
      .def_ro("dy", &SimUtils::SplineAkimaInterpolate::dy)
      .def("__len__", [](const SimUtils::SplineAkimaInterpolate &) { return 3; })
      .def("__getitem__", [](const SimUtils::SplineAkimaInterpolate &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.ok);
        if (i == 1)
          return nb::cast(s.y);
        if (i == 2)
          return nb::cast(s.dy);
        throw nb::index_error();
      });
  m.def(
      "spline_akima_interpolate",
      &SimUtils::spline_akima_interpolate,
      nb::arg("x_knot"),
      nb::arg("y_knot"),
      nb::arg("x"),
      R"""(Routine to interpolate using an akima spline.

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
  nb::class_<SimUtils::SplineEvaluate>(m, "SplineEvaluate", "spline_evaluate return type")
      .def_ro("ok", &SimUtils::SplineEvaluate::ok)
      .def_ro("y", &SimUtils::SplineEvaluate::y)
      .def_ro("dy", &SimUtils::SplineEvaluate::dy)
      .def("__len__", [](const SimUtils::SplineEvaluate &) { return 3; })
      .def("__getitem__", [](const SimUtils::SplineEvaluate &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.ok);
        if (i == 1)
          return nb::cast(s.y);
        if (i == 2)
          return nb::cast(s.dy);
        throw nb::index_error();
      });
  m.def(
      "spline_evaluate",
      &SimUtils::spline_evaluate,
      nb::arg("spline"),
      nb::arg("x"),
      R"""(Subroutine to evalueate a spline at a set of points.

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
      nb::arg("alpha"),
      nb::arg("x"),
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
      nb::arg("x"),
      nb::arg("nd") = nb::none(),
      R"""(Wrapper for Fortran routine sqrt_one

Parameters
----------
x : float
    Number

nd : int, optional
    Derivative order. nd = 0 (default) -> compute Sqrt[1+x] - 1. NOTE: Currently only nd = 0 and nd = 1 are
    implemented.

Returns
-------
ds1 : float
)"""
  );
  m.def(
      "str_count",
      &SimUtils::str_count,
      nb::arg("str"),
      nb::arg("match_"),
      R"""(Wrapper for Fortran routine str_count

Parameters
----------
str : str

match : str

Returns
-------
num : int
)"""
  );
  m.def(
      "str_downcase",
      &SimUtils::str_downcase,
      nb::arg("src"),
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
      nb::arg("line"),
      nb::arg("set"),
      nb::arg("ignore_clauses") = nb::none(),
      R"""(Wrapper for Fortran routine str_first_in_set

Parameters
----------
line : str

set : str

ignore_clauses : bool, optional

Returns
-------
ix_match : int
)"""
  );
  m.def(
      "str_first_not_in_set",
      &SimUtils::str_first_not_in_set,
      nb::arg("line"),
      nb::arg("set"),
      R"""(Wrapper for Fortran routine str_first_not_in_set

Parameters
----------
line : str

set : str

Returns
-------
ix_match : int
)"""
  );
  m.def(
      "str_last_in_set",
      &SimUtils::str_last_in_set,
      nb::arg("line"),
      nb::arg("set"),
      R"""(Wrapper for Fortran routine str_last_in_set

Parameters
----------
line : str

set : str

Returns
-------
ix_match : int
)"""
  );
  m.def(
      "str_last_not_in_set",
      &SimUtils::str_last_not_in_set,
      nb::arg("line"),
      nb::arg("set"),
      R"""(Wrapper for Fortran routine str_last_not_in_set

Parameters
----------
line : str

set : str

Returns
-------
ix_match : int
)"""
  );
  m.def(
      "str_match_wild",
      &SimUtils::str_match_wild,
      nb::arg("str"),
      nb::arg("pat"),
      R"""(Wrapper for Fortran routine str_match_wild

Parameters
----------
str : str

pat : str

Returns
-------
a_match : bool
)"""
  );
  m.def(
      "str_substitute",
      &SimUtils::str_substitute,
      nb::arg("string"),
      nb::arg("str_match") = nb::none(),
      nb::arg("str_replace") = nb::none(),
      nb::arg("do_trim") = nb::none(),
      nb::arg("ignore_escaped") = nb::none(),
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
      nb::arg("src"),
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
      nb::arg("line"),
      nb::arg("default_"),
      nb::arg("err_flag"),
      nb::arg("err_print_flag") = nb::none(),
      R"""(Wrapper for Fortran routine string_to_int

Parameters
----------
line : str

err_flag : bool

err_print_flag : bool, optional

Returns
-------
value : int
)"""
  );
  m.def(
      "string_to_real",
      &SimUtils::string_to_real,
      nb::arg("line"),
      nb::arg("default_"),
      nb::arg("err_flag"),
      nb::arg("err_print_flag") = nb::none(),
      R"""(Wrapper for Fortran routine string_to_real

Parameters
----------
line : str

err_flag : bool

err_print_flag : bool, optional

Returns
-------
value : float
)"""
  );
  m.def(
      "string_trim",
      &SimUtils::string_trim,
      nb::arg("in_string"),
      nb::arg("out_string"),
      nb::arg("word_len"),
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
      nb::arg("in_str"),
      nb::arg("delimitors"),
      nb::arg("out_str"),
      nb::arg("ix_word"),
      nb::arg("delim"),
      nb::arg("ix_next"),
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
      nb::arg("XV"),
      nb::arg("FV"),
      nb::arg("eps"),
      nb::arg("itermx"),
      nb::arg("reset_flag") = nb::none(),
      R"""(Reverse communication subroutine.
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
      nb::arg("y"),
      nb::arg("y1"),
      nb::arg("y2"),
      nb::arg("y12"),
      nb::arg("d1"),
      nb::arg("d2"),
      R"""(Routine to compute coefficients for bicubic interpolation.
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
  nb::class_<SimUtils::SuperBicubicInterpolation>(
      m,
      "SuperBicubicInterpolation",
      "super_bicubic_interpolation return type"
  )
      .def_ro("ansy", &SimUtils::SuperBicubicInterpolation::ansy)
      .def_ro("ansy1", &SimUtils::SuperBicubicInterpolation::ansy1)
      .def_ro("ansy2", &SimUtils::SuperBicubicInterpolation::ansy2)
      .def("__len__", [](const SimUtils::SuperBicubicInterpolation &) { return 3; })
      .def("__getitem__", [](const SimUtils::SuperBicubicInterpolation &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.ansy);
        if (i == 1)
          return nb::cast(s.ansy1);
        if (i == 2)
          return nb::cast(s.ansy2);
        throw nb::index_error();
      });
  m.def(
      "super_bicubic_interpolation",
      &SimUtils::super_bicubic_interpolation,
      nb::arg("y"),
      nb::arg("y1"),
      nb::arg("y2"),
      nb::arg("y12"),
      nb::arg("x1l"),
      nb::arg("x1u"),
      nb::arg("x2l"),
      nb::arg("x2u"),
      nb::arg("x1"),
      nb::arg("x2"),
      R"""(Routine to do bicubic interpolation.
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
  nb::class_<SimUtils::SuperPolint>(m, "SuperPolint", "super_polint return type")
      .def_ro("y", &SimUtils::SuperPolint::y)
      .def_ro("dy", &SimUtils::SuperPolint::dy)
      .def("__len__", [](const SimUtils::SuperPolint &) { return 2; })
      .def("__getitem__", [](const SimUtils::SuperPolint &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.y);
        if (i == 1)
          return nb::cast(s.dy);
        throw nb::index_error();
      });
  m.def(
      "super_polint",
      &SimUtils::super_polint,
      nb::arg("xa"),
      nb::arg("ya"),
      nb::arg("x"),
      R"""(This is essentially polint from Numerical Recipes.

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
      nb::arg("x"),
      nb::arg("coeffs"),
      R"""(Routine to compute Sum: coef(i)*x^i

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
      [](FArray1D<Real> &x, RandomStateStruct *ran_state) {
        auto fn = static_cast<void (*)(FArray1D<Real> &, optional_ref<RandomStateStruct>)>(
            &SimUtils::super_sobseq
        );
        return fn(x, ptr_to_opt_ref(ran_state));
      },
      nb::arg("x"),
      nb::arg("ran_state") = nb::none(),
      R"""(Routine patterened after sobseq in Numerical Recipes.
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
      nb::arg("arr"),
      R"""(Routine to sort an integer array in place.
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
      nb::arg("line"),
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
