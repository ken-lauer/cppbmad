#include "pybmad/generated/SimUtils_routines_m.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_m(py::module &m) {
  m.def(
      "make_legal_comment",
      &SimUtils::make_legal_comment,
      py::arg("comment_in"),
      py::arg("comment_out"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine make_legal_comment

Parameters
----------
comment_in : str

comment_out : str
)"""
  );
  m.def(
      "mass_of",
      &SimUtils::mass_of,
      py::arg("species"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function mass_of (species) result (mass)

Routine to return the mass, in units of eV/c^2, of a particle.
To convert to AMU divide mass_of value by the constant atomic_mass_unit.

Note: For atoms where the isotopic number is given, the mass is calculated using the neutral atomic mass
adjusted by the weight of any added or missing electrons. The calculated mass is off very slightly due to
binding energy effects. Exception: For #1H+ (proton) and #2H+ (deuteron) the exact mass is used since it is known.

Parameters
----------
species : int
    Species ID.

Returns
-------
mass : float
    particle mass. Set to real_garbage$ if species value is invalid.
)"""
  );
  m.def(
      "match_reg",
      &SimUtils::match_reg,
      py::arg("str"),
      py::arg("pat"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine match_reg

Parameters
----------
str : str

pat : str

Returns
-------
is_match : bool
)"""
  );
  m.def(
      "match_wild",
      &SimUtils::match_wild,
      py::arg("string"),
      py::arg("template_"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine match_wild

Parameters
----------
string : str

Returns
-------
is_match : bool
)"""
  );
  m.def(
      "match_word",
      &SimUtils::match_word,
      py::arg("string"),
      py::arg("names"),
      py::arg("ix"),
      py::arg("exact_case") = py::none(),
      py::arg("can_abbreviate") = py::none(),
      py::arg("matched_name") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine match_word

Parameters
----------
string : str

names : 1D array of str

ix : int

exact_case : bool, optional

can_abbreviate : bool, optional

matched_name : str, optional
)"""
  );
  m.def(
      "maximize_projection",
      &SimUtils::maximize_projection,
      py::arg("seed"),
      py::arg("cdata"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(function maximize_projection

Optimizer that uses Numerical Recipes brent to find a local maximum,
which is the frequency that maximizes the projection.
)"""
  );
  m.def(
      "milli_sleep",
      &SimUtils::milli_sleep,
      py::arg("milli_sec"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine milli_sleep

Parameters
----------
milli_sec : int
)"""
  );
  m.def(
      "modulo2_dp",
      &SimUtils::modulo2_dp,
      py::arg("x"),
      py::arg("amp"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function modulo2 (x, amp) result (mod2)

Function to return
    mod2 = x + 2 * n * amp
where n is an integer chosen such that
   -amp <= mod2 < amp

Parameters
----------
x : float
    Real(sp), Real(rp), or Integer

amp : float
    Must be positive.

Returns
-------
mod2 : float
    Result
)"""
  );
  m.def(
      "modulo2_int",
      &SimUtils::modulo2_int,
      py::arg("x"),
      py::arg("amp"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function modulo2 (x, amp) result (mod2)

Function to return
    mod2 = x + 2 * n * amp
where n is an integer chosen such that
   -amp <= mod2 < amp

Parameters
----------
x : int
    Real(sp), Real(rp), or Integer

amp : int
    Must be positive.

Returns
-------
mod2 : int
    Result
)"""
  );
  m.def(
      "modulo2_qp",
      &SimUtils::modulo2_qp,
      py::arg("x"),
      py::arg("amp"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function modulo2 (x, amp) result (mod2)

Function to return
    mod2 = x + 2 * n * amp
where n is an integer chosen such that
   -amp <= mod2 < amp

Parameters
----------
x : float
    Real(sp), Real(rp), or Integer

amp : float
    Must be positive.

Returns
-------
mod2 : float
    Result
)"""
  );
  m.def(
      "modulo2_sp",
      &SimUtils::modulo2_sp,
      py::arg("x"),
      py::arg("amp"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function modulo2 (x, amp) result (mod2)

Function to return
    mod2 = x + 2 * n * amp
where n is an integer chosen such that
   -amp <= mod2 < amp

Parameters
----------
x : float
    Real(sp), Real(rp), or Integer

amp : float
    Must be positive.

Returns
-------
mod2 : float
    Result
)"""
  );
}
