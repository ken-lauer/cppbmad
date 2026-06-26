#include "pybmad/generated/SimUtils_routines_m.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_m(nb::module_ &m) {
  m.def(
      "make_legal_comment",
      &SimUtils::make_legal_comment,
      nb::arg("comment_in"),
      nb::arg("comment_out"),
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
      nb::arg("species"),
      R"""(Routine to return the mass, in units of eV/c^2, of a particle.
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
      nb::arg("str"),
      nb::arg("pat"),
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
      nb::arg("string"),
      nb::arg("template_"),
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
      nb::arg("string"),
      nb::arg("names"),
      nb::arg("ix"),
      nb::arg("exact_case") = nb::none(),
      nb::arg("can_abbreviate") = nb::none(),
      nb::arg("matched_name") = nb::none(),
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
      nb::arg("seed"),
      nb::arg("cdata"),
      R"""(Optimizer that uses Numerical Recipes brent to find a local maximum,
which is the frequency that maximizes the projection.
)"""
  );
  m.def(
      "milli_sleep",
      &SimUtils::milli_sleep,
      nb::arg("milli_sec"),
      R"""(Wrapper for Fortran routine milli_sleep

Parameters
----------
milli_sec : int
)"""
  );
  m.def(
      "modulo2_dp",
      &SimUtils::modulo2_dp,
      nb::arg("x"),
      nb::arg("amp"),
      R"""(Function to return
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
      nb::arg("x"),
      nb::arg("amp"),
      R"""(Function to return
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
      nb::arg("x"),
      nb::arg("amp"),
      R"""(Function to return
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
      nb::arg("x"),
      nb::arg("amp"),
      R"""(Function to return
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
      "molecular_components",
      &SimUtils::molecular_components,
      nb::arg("molecule"),
      R"""(Routine to decompose a molecule into its components.
For example: molecule = 'H2O' => component(1) = ('H', 2), component(2) = ('O', 1)

Parameters
----------
molecule : str
    Molecular name.

Returns
-------
component : 1D array of MolecularComponentStruct
    Array of components.
)"""
  );
}
