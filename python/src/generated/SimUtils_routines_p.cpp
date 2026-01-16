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
      R"""(Parameters
----------
format_str : 
n_repeat : 
power : 
descrip : 
width : 
digits : 
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
ran_state_ptr : RandomStateStruct
    Pointer to the appropriate state.
)"""
  );
  m.def(
      "poly_eval",
      &SimUtils::poly_eval,
      py::arg("poly"),
      py::arg("x"),
      py::arg("diff_coef") = py::none(),
      R"""(Parameters
----------
poly : float
    Polynomial
x : float
    Point to evaluate at.
diff_coef : bool, optional
    poly(:) array are differentials? Default is False.
y : float
    Value of polynomial.
)"""
  );
  m.def(
      "probability_funct",
      &SimUtils::probability_funct,
      py::arg("x"),
      py::arg("prob"),
      R"""(Parameters
----------
x : float
    Function argument.
prob : 
)"""
  );
  m.def(
      "projdd",
      &SimUtils::projdd,
      py::arg("a"),
      py::arg("b"),
      py::arg("func_retval__"),
      R"""(Parameters
----------
a : 
b : 
projdd : 
)"""
  );
}
