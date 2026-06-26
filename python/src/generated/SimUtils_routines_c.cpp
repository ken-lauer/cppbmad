#include "pybmad/generated/SimUtils_routines_c.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_c(nb::module_ &m) {
  m.def(
      "calc_file_number",
      &SimUtils::calc_file_number,
      nb::arg("file_name"),
      nb::arg("num_in"),
      nb::arg("num_out"),
      nb::arg("err_flag"),
      R"""(Wrapper for Fortran routine calc_file_number

Parameters
----------
file_name : str

num_in : int

num_out : int

err_flag : bool
)"""
  );
  nb::class_<SimUtils::Celbd>(m, "Celbd", "celbd return type")
      .def_ro("elb", &SimUtils::Celbd::elb)
      .def_ro("eld", &SimUtils::Celbd::eld)
      .def("__len__", [](const SimUtils::Celbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Celbd &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.elb);
        if (i == 1)
          return nb::cast(s.eld);
        throw nb::index_error();
      });
  m.def(
      "celbd",
      &SimUtils::celbd,
      nb::arg("mc"),
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
      nb::arg("i_arg"),
      R"""(Platform independent function to return the i'th command line argument.
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
      R"""(Note: Use the Fortran intrinsic command_argument_count instead

Platform independent function to return the number of command line arguments.
Use this with cesr_getarg.
)"""
  );
  m.def(
      "change_file_number",
      &SimUtils::change_file_number,
      nb::arg("file_name"),
      nb::arg("change"),
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
      nb::arg("species"),
      nb::arg("default_") = nb::none(),
      R"""(Routine to return the charge, in units of e+, of a particle.

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
      nb::arg("species"),
      R"""(Routine to return the charge (in units of e+) to mass (in units of eV) ratio of a particle.

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
      nb::arg("data"),
      nb::arg("error") = nb::none(),
      R"""(Simple function to take periodic data and estimate
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
      nb::arg("wr"),
      nb::arg("wi"),
      nb::arg("zr"),
      nb::arg("zi"),
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
      nb::arg("angle"),
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
      nb::arg("x"),
      nb::arg("nd") = nb::none(),
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
      nb::arg("r0"),
      nb::arg("r1"),
      nb::arg("slope0"),
      nb::arg("slope1"),
      R"""(Routine to create a single spline given end point positions and slopes.
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
      nb::arg("a"),
      nb::arg("b"),
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
