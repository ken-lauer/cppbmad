#include "pybmad/generated/SimUtils_routines_p.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_p(py::module &m) {
  m.def(
      "parse_fortran_format",
      &SimUtils::parse_fortran_format,
      py::arg("format_str"),
      py::arg("n_repeat"),
      py::arg("power"),
      py::arg("descrip"),
      py::arg("width"),
      py::arg("digits"),
      R"""(Wrapper for Fortran routine parse_fortran_format

Parameters
----------
format_str : character

n_repeat : int

power : int

descrip : character

width : int

digits : int
)"""
  );
  m.def(
      "pointer_to_ran_state",
      &SimUtils::pointer_to_ran_state,
      py::arg("ran_state") = py::none(),
      py::arg("ix_thread") = py::none(),
      R"""(Function pointer_to_ran_state(ran_state, ix_thread) result (ran_state_ptr)

Routine to point to the appropriate state structure for generating random numbers

Parameters
----------
ran_state : RandomStateStruct, optional
    Point to this if present. Otherwise point to the global saved state.

ix_thread : int, optional
    Thread index.

Returns
-------
ran_state_ptr : RandomStateStruct, optional
    Pointer to the appropriate state.
)"""
  );
  m.def(
      "poly_eval",
      &SimUtils::poly_eval,
      py::arg("poly"),
      py::arg("x"),
      py::arg("diff_coef") = py::none(),
      R"""(Wrapper for Fortran routine poly_eval

Parameters
----------
poly : 1D array of float
    Polynomial

x : float
    Point to evaluate at.

diff_coef : bool, optional
    poly(:) array are differentials? Default is False.

Returns
-------
y : float
    Value of polynomial.
)"""
  );
  m.def(
      "probability_funct",
      &SimUtils::probability_funct,
      py::arg("x"),
      py::arg("prob"),
      R"""(Wrapper for Fortran routine probability_funct

Parameters
----------
x : float
    Function argument.

prob : float
)"""
  );
  m.def(
      "projdd",
      &SimUtils::projdd,
      py::arg("a"),
      py::arg("b"),
      py::arg("func_retval__"),
      R"""(Wrapper for Fortran routine projdd

Parameters
----------
a : 1D array of complex

b : 1D array of complex
)"""
  );
}
