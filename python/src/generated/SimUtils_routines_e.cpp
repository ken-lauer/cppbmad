#include "pybmad/generated/SimUtils_routines_e.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_e(py::module &m) {
  m.def(
      "end_akima_spline_calc",
      &SimUtils::end_akima_spline_calc,
      py::arg("spline"),
      py::arg("which_end"),
      R"""(Subroutine end_akima_spline_calc (spline, which_end)

Routine to calculate the slopes at the ends of a spline array

Parameters
----------
spline : 1D array of SplineStruct
    Array of splines.
    This parameter is an input/output and is modified in-place.
    As an output, spline: Array with slopes at end calculated.

which_end : int
    0 => calculate slopes for the start end of the array. 1 => calculate slopes for the end end of the array.

Returns
-------
spline : 1D array of SplineStruct
    Array of splines.
    This parameter is an input/output and is modified in-place.
    As an output, spline: Array with slopes at end calculated.
)"""
  );
  m.def(
      "err_exit",
      &SimUtils::err_exit,
      py::arg("err_str") = py::none(),
      R"""(Wrapper for Fortran routine err_exit

Parameters
----------
err_str : character, optional

Returns
-------
err_str : character, optional
)"""
  );
}
