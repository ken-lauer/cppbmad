#include "pybmad/generated/SimUtils_routines_c.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_c(py::module &m) {
  m.def(
      "calc_file_number",
      &SimUtils::calc_file_number,
      py::arg("file_name"),
      py::arg("num_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine calc_file_number

Parameters
----------
file_name : str

num_in : int

num_out : int

err_flag : bool
)"""
  );
  py::class_<SimUtils::Celbd, std::unique_ptr<SimUtils::Celbd>>(m, "Celbd", "celbd return type")
      .def_readonly("elb", &SimUtils::Celbd::elb)
      .def_readonly("eld", &SimUtils::Celbd::eld)
      .def("__len__", [](const SimUtils::Celbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Celbd &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.elb);
        if (i == 1)
          return py::cast(s.eld);
        throw py::index_error();
      });
  m.def(
      "celbd",
      &SimUtils::celbd,
      py::arg("mc"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine celbd

Parameters
----------
mc : float

Returns
-------
elb : float

eld : float
)"""
  );
  m.def(
      "cesr_getarg",
      &SimUtils::cesr_getarg,
      py::arg("i_arg"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine cesr_getarg (i_arg, arg)

Platform independent function to return the i'th command line argument.
Use this with cesr_iargc.

Note: The difference between this routine and the Fortran instrinsic
get_command_argument is that for i_arg = 0, this routine returns the
command line with the name of the executable removed from the beginning of the line.
get_command_argument, on the other hand returns the name of the
executable when the argument is 0.

Parameters
----------
i_arg : int
    Index of argument to return. i_arg = 0 => Entire line minus the executable string. i_arg = 1 => First
    argument.

Returns
-------
arg : str
    i'th command line argument. If i_arg > number_of_args then arg is a blank string.
)"""
  );
  m.def(
      "cesr_iargc",
      &SimUtils::cesr_iargc,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function cesr_iargc ()

Note: Use the Fortran intrinsic command_argument_count instead

Platform independent function to return the number of command line arguments.
Use this with cesr_getarg.
)"""
  );
  m.def(
      "change_file_number",
      &SimUtils::change_file_number,
      py::arg("file_name"),
      py::arg("change"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine change_file_number

Parameters
----------
file_name : str

change : int
)"""
  );
  m.def(
      "charge_of",
      &SimUtils::charge_of,
      py::arg("species"),
      py::arg("default_") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function charge_of (species, default) result (charge)

Routine to return the charge, in units of e+, of a particle.

Parameters
----------
species : int
    Species ID.

Returns
-------
charge : int
    particle charge.
)"""
  );
  m.def(
      "charge_to_mass_of",
      &SimUtils::charge_to_mass_of,
      py::arg("species"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function charge_to_mass_of (species) result (charge_mass_ratio)

Routine to return the charge (in units of e+) to mass (in units of eV) ratio of a particle.

Parameters
----------
species : int
    Species ID.

Returns
-------
charge_mass_ratio : float
    particle charge to mass ratio. (1/eV)
)"""
  );
  m.def(
      "coarse_frequency_estimate",
      &SimUtils::coarse_frequency_estimate,
      py::arg("data"),
      py::arg("error") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function coarse_frequency_estimate(data, error) result(frequency)

Simple function to take periodic data and estimate
the most dominant frequency by FFT.

Parameters
----------
data : 1D array of float
    data to analyze. Preferably size(data) is a power of 2 Otherwise the data is padded with zeros.

Returns
-------
frequency : float
    Frequency corresponding to the largest FFT amplitude
)"""
  );
  m.def(
      "complex_error_function",
      &SimUtils::complex_error_function,
      py::arg("wr"),
      py::arg("wi"),
      py::arg("zr"),
      py::arg("zi"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine complex_error_function

Parameters
----------
wr : float

wi : float

zr : float

zi : float
)"""
  );
  m.def(
      "cos_one",
      &SimUtils::cos_one,
      py::arg("angle"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine cos_one

Parameters
----------
angle : float
    Angle.

Returns
-------
cos1 : float
    Result.
)"""
  );
  m.def(
      "cosc",
      &SimUtils::cosc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine cosc

Parameters
----------
x : float

nd : int, optional
    Derivative order. nd = 0 (default) -> compute (1 - cos(x)) / x^2 NOTE: Currently only nd = 0 and nd = 1
    are implemented.

Returns
-------
y : float
    nd^th derivative of (1 - cos(x)) / x^2
)"""
  );
  m.def(
      "create_a_spline",
      &SimUtils::create_a_spline,
      py::arg("r0"),
      py::arg("r1"),
      py::arg("slope0"),
      py::arg("slope1"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function create_a_spline (r0, r1, slope0, slope1) result (spline)

Routine to create a single spline given end point positions and slopes.
The spline will pass through the data points and have the given slopes
at these points.

Modules used:
  use spline_mod

Parameters
----------
r0 : 1D array of float
    Start (x, y) point.

r1 : 1D array of float
    End (x, y) point.

slope0 : float
    Starting slope.

slope1 : float
    End slope.

Returns
-------
spline : SplineStruct
    Spline.
)"""
  );
  m.def(
      "cross_product",
      &SimUtils::cross_product,
      py::arg("a"),
      py::arg("b"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine cross_product

Parameters
----------
a : 1D array of float
    Input vectors.

b : 1D array of float
    Input vectors.

Returns
-------
c : 1D array of float (shape: 3)
    Cross product: a X b.
)"""
  );
}
