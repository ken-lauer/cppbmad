#include "pybmad/generated/Tao_routines_r.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyReAllocateCDouble {
  std::optional<double> init_val;
};
PyReAllocateCDouble python_re_allocate_c_double(
    RealAlloc1D& re,
    int n,
    std::optional<bool> exact = std::nullopt,
    std::optional<double> init_val = std::nullopt) {
  Tao::re_allocate_c_double(re, n, exact, make_opt_ref(init_val));
  auto py_result{PyReAllocateCDouble{init_val}};
  return py_result;
}

void init_Tao_routines_r(py::module& m) {
  m.def(
      "re_allocate_c_double",
      &python_re_allocate_c_double,
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
  )""");
  py::class_<PyReAllocateCDouble, std::unique_ptr<PyReAllocateCDouble>>(
      m,
      "ReAllocateCDouble",
      "Fortran routine re_allocate_c_double return value")
      .def_readonly("init_val", &PyReAllocateCDouble::init_val)
      .def("__len__", [](const PyReAllocateCDouble&) { return 1; })
      .def(
          "__getitem__", [](const PyReAllocateCDouble& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.init_val);
            throw py::index_error();
          });
}
