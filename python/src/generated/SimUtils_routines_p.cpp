#include "pybmad/generated/SimUtils_routines_p.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_p(nb::module_ &m) {
  m.def(
      "parse_fortran_format",
      &SimUtils::parse_fortran_format,
      nb::arg("format_str"),
      nb::arg("n_repeat"),
      nb::arg("power"),
      nb::arg("descrip"),
      nb::arg("width"),
      nb::arg("digits"),
      R"""(Wrapper for Fortran routine parse_fortran_format

Parameters
----------
format_str : str

n_repeat : int

power : int

descrip : str

width : int

digits : int
)"""
  );
  m.def(
      "pointer_to_locations",
      [](std::string string,
         IntAlloc1D &array,
         int num,
         int ix_min,
         int ix_max,
         CharacterAlloc1D *names,
         std::optional<bool> exact_case,
         std::optional<bool> print_err) {
        auto fn = static_cast<void (*)(
            std::string,
            IntAlloc1D &,
            int,
            int,
            int,
            optional_ref<CharacterAlloc1D>,
            std::optional<bool>,
            std::optional<bool>
        )>(&SimUtils::pointer_to_locations);
        return fn(string, array, num, ix_min, ix_max, ptr_to_opt_ref(names), exact_case, print_err);
      },
      nb::arg("string"),
      nb::arg("array"),
      nb::arg("num"),
      nb::arg("ix_min"),
      nb::arg("ix_max"),
      nb::arg("names") = nb::none(),
      nb::arg("exact_case") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Wrapper for Fortran routine pointer_to_locations

Parameters
----------
string : str

array : 1D array of int

num : int

ix_min : int

ix_max : int

names : 1D array of str, optional

exact_case : bool, optional

print_err : bool, optional
)"""
  );
  m.def(
      "pointer_to_ran_state",
      [](RandomStateStruct *ran_state, std::optional<int> ix_thread) {
        auto fn = static_cast<std::optional<
            RandomStateStruct> (*)(optional_ref<RandomStateStruct>, std::optional<int>)>(
            &SimUtils::pointer_to_ran_state
        );
        return fn(ptr_to_opt_ref(ran_state), ix_thread);
      },
      nb::arg("ran_state") = nb::none(),
      nb::arg("ix_thread") = nb::none(),
      R"""(Routine to point to the appropriate state structure for generating random numbers

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
      nb::arg("poly"),
      nb::arg("x"),
      nb::arg("diff_coef") = nb::none(),
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
      nb::arg("x"),
      R"""(Wrapper for Fortran routine probability_funct

Parameters
----------
x : float
    Function argument.

Returns
-------
prob : float
)"""
  );
  m.def(
      "projdd",
      &SimUtils::projdd,
      nb::arg("a"),
      nb::arg("b"),
      R"""(Wrapper for Fortran routine projdd

Parameters
----------
a : 1D array of complex

b : 1D array of complex
)"""
  );
}
