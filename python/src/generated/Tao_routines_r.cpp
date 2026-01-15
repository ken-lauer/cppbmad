#include "pybmad/generated/Tao_routines_r.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Tao_routines_r(py::module &m) {
  m.def(
      "re_allocate_c_double",
      &Tao::re_allocate_c_double,
      py::arg("re"),
      py::arg("n"),
      py::arg("exact") = py::none(),
      py::arg("init_val") = py::none(),
      R"""(Subroutine re_allocate_c_double (re, n, exact, init_val)

Routine to reallocate an array of c_double reals.
This is modeled after the reallocate functions in Numerical Recipes.
Note: The data of the array is preserved but data at the end of the
array will be lost if n is less than the original size of the array

Parameters
----------
re : float
    Real array.
    This parameter is an input/output and is modified in-place. As an output: Allocated array with size(re) >=
    n.
n : int
    Size wanted.
exact : bool, optional
    If present and False then the size of the output array is permitted to be larger than n. Default is True.
)"""
  );
}
