#include "pybmad/generated/SimUtils_routines_e.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyErrExit python_err_exit(std::optional<std::string> err_str = std::nullopt) {
  SimUtils::err_exit(make_opt_ref(err_str));
  auto py_result{PyErrExit{err_str}};
  return py_result;
}

void init_SimUtils_routines_e(py::module& m) {
  m.def(
      "end_akima_spline_calc",
      &SimUtils::end_akima_spline_calc,
      py::arg("spline"),
      py::arg("which_end"),
      R"""(Subroutine end_akima_spline_calc (spline, which_end)

  Routine to calculate the slopes at the ends of a spline array

  Parameters
  ----------
  spline : SplineStruct
      Array of splines.
      This parameter is an input/output and is modified in-place. As an output: Array with slopes at end
      calculated.
  which_end : int
      0 => calculate slopes for the start end of the array. 1 => calculate slopes for the end end of the array.
  )""");
  py::class_<PyErrExit, std::unique_ptr<PyErrExit>>(
      m, "ErrExit", "Fortran routine err_exit return value")
      .def_readonly("err_str", &PyErrExit::err_str)
      .def("__len__", [](const PyErrExit&) { return 1; })
      .def("__getitem__", [](const PyErrExit& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.err_str);
        throw py::index_error();
      });
  m.def(
      "err_exit",
      &python_err_exit,
      py::arg("err_str") = py::none(),
      R"""(Parameters
  ----------
  err_str : 
  )""");
}
