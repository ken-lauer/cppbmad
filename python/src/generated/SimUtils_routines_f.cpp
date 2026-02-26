#include "pybmad/generated/SimUtils_routines_f.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_f(py::module &m) {
  m.def(
      "factorial",
      &SimUtils::factorial,
      py::arg("n"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine factorial

Parameters
----------
n : int
    Must be non-negative

Returns
-------
fact : float
    n!. Will return negative number if there is an error.
)"""
  );
  m.def(
      "faddeeva_function",
      &SimUtils::faddeeva_function,
      py::arg("z"),
      py::arg("w"),
      py::arg("dw"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine faddeeva_function

Parameters
----------
z : 1D array of float (shape: 2)

w : 1D array of float (shape: 2)

dw : 2D array of float (shape: 2,2)
)"""
  );
  m.def(
      "fft_1d",
      &SimUtils::fft_1d,
      py::arg("arr"),
      py::arg("isign"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(no longer exists
subroutine fff_sub(line, error)
  implicit none
  character(*) line
  logical error
end subroutine

Parameters
----------
arr : 1D array of complex
    Input array.
    This parameter is an input/output and is modified in-place.
    As an output, arr: FFT of array.

isign : int
    -1 => "Forward" transform, +1 => "Backwards" transform.
)"""
  );
  m.def(
      "file_directorizer",
      &SimUtils::file_directorizer,
      py::arg("in_file"),
      py::arg("out_file"),
      py::arg("directory"),
      py::arg("add_switch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine file_directorizer

Parameters
----------
in_file : str

out_file : str

directory : str

add_switch : bool
)"""
  );
  m.def(
      "file_get",
      &SimUtils::file_get,
      py::arg("string"),
      py::arg("dflt_file_name"),
      py::arg("file_name"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine file_get

Parameters
----------
string : str

dflt_file_name : str

file_name : str
)"""
  );
  m.def(
      "file_get_open",
      &SimUtils::file_get_open,
      py::arg("string"),
      py::arg("dflt_file_name"),
      py::arg("file_name"),
      py::arg("file_unit"),
      py::arg("readonly"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine file_get_open

Parameters
----------
string : str

dflt_file_name : str

file_name : str

file_unit : int

readonly : bool
)"""
  );
  m.def(
      "file_suffixer",
      &SimUtils::file_suffixer,
      py::arg("in_file_name"),
      py::arg("out_file_name"),
      py::arg("suffix"),
      py::arg("add_switch"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine file_suffixer

Parameters
----------
in_file_name : str

out_file_name : str

suffix : str

add_switch : bool
)"""
  );
  m.def(
      "find_location",
      py::overload_cast<FArray1D<Int> &, int, int>(&SimUtils::find_location),
      py::arg("arr"),
      py::arg("value"),
      py::arg("ix_match"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine find_location_int

Parameters
----------
arr : 1D array of int

value : int

ix_match : int
)"""
  );
  m.def(
      "find_location",
      py::overload_cast<BoolAlloc1D &, bool, int>(&SimUtils::find_location),
      py::arg("arr"),
      py::arg("value"),
      py::arg("ix_match"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine find_location_logic

Parameters
----------
arr : 1D array of bool

value : bool

ix_match : int
)"""
  );
  m.def(
      "find_location",
      py::overload_cast<FArray1D<Real> &, double>(&SimUtils::find_location),
      py::arg("arr"),
      py::arg("value"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine find_location_real

Parameters
----------
arr : 1D array of float
    real(rp), logical, or integer

value : float
    :).

Returns
-------
ix_match : int
    Index of match. Zero if no match found.
)"""
  );
  m.def(
      "fine_frequency_estimate",
      &SimUtils::fine_frequency_estimate,
      py::arg("data"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function fine_frequency_estimate(data) result(frequency)

Uses Laskar's method to accurately find the most dominant frequency
A coarse estimate is first made by FFT.

Parameters
----------
data : 1D array of float
    data to analyze

Returns
-------
frequency : float
    Frequency corresponding to the largest FFT amplitude
)"""
  );
  m.def(
      "fixedwindowls",
      &SimUtils::fixedwindowls,
      py::arg("ynew"),
      py::arg("id"),
      py::arg("z"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function fixedWindowLS

Main function of the windowLS modult.  Each call to this function adds a data point to the fit
and returns the derivative evaluated at the end of the window.  It is assumed that all data points
are separeted by the same interval.
This module is initialized with zeros for all data points, and so the results are unreliable until
a number of data points equal to N has been entered.

initFixedWindowLS must be called prior to calling this function.  destFixedWindowLS should be
called when the instance is no longer needed.

Parameters
----------
)"""
  );
  py::class_<SimUtils::FourierAmplitude, std::unique_ptr<SimUtils::FourierAmplitude>>(
      m,
      "FourierAmplitude",
      "fourier_amplitude return type"
  )
      .def_readonly("cos_amp", &SimUtils::FourierAmplitude::cos_amp)
      .def_readonly("sin_amp", &SimUtils::FourierAmplitude::sin_amp)
      .def_readonly("dcos_amp", &SimUtils::FourierAmplitude::dcos_amp)
      .def_readonly("dsin_amp", &SimUtils::FourierAmplitude::dsin_amp)
      .def("__len__", [](const SimUtils::FourierAmplitude &) { return 4; })
      .def("__getitem__", [](const SimUtils::FourierAmplitude &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.cos_amp);
        if (i == 1)
          return py::cast(s.sin_amp);
        if (i == 2)
          return py::cast(s.dcos_amp);
        if (i == 3)
          return py::cast(s.dsin_amp);
        throw py::index_error();
      });
  m.def(
      "fourier_amplitude",
      &SimUtils::fourier_amplitude,
      py::arg("data"),
      py::arg("frequency"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine fourier_amplitude(data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp)

Computes cos_amp = (1/N) * sum_n=0^{N-1} data(n-1) cos(twopi*frequency*n)
    and  sin_amp = (1/N) * sum_n=0^{N-1} data(n-1) sin(twopi*frequency*n)
    and optionally dcos_amp = d/dfrequency cos_amp
                   dsin_amp = d/dfrequency sin_amp

Parameters
----------
data : 1D array of float
    data to analyze

frequency : float
    frequency

Returns
-------
cos_amp : float
    cosine amplitude

sin_amp : float
    sine amplitude

dcos_amp : float, optional
    cosine amplitude derivative

dsin_amp : float, optional
    sine amplitude derivative
)"""
  );
}
