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
      R"""(Parameters
----------
comment_in : 
comment_out : 
)"""
  );
  m.def(
      "mass_of",
      &SimUtils::mass_of,
      py::arg("species"),
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
      py::arg("is_match"),
      R"""(Parameters
----------
str : 
pat : 
is_match : 
)"""
  );
  m.def(
      "match_wild",
      &SimUtils::match_wild,
      py::arg("string"),
      py::arg("template_"),
      py::arg("is_match"),
      R"""(Parameters
----------
string : 
template : 
is_match : 
)"""
  );
  m.def(
      "maximize_projection",
      &SimUtils::maximize_projection,
      py::arg("seed"),
      py::arg("cdata"),
      py::arg("func_retval__"),
      R"""(function maximize_projection

Optimizer that uses Numerical Recipes brent to find a local maximum,
which is the frequency that maximizes the projection.


)"""
  );
  m.def(
      "milli_sleep",
      &SimUtils::milli_sleep,
      py::arg("milli_sec"),
      R"""(Parameters
----------
milli_sec : 
)"""
  );
}
