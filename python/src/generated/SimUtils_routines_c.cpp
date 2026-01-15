#include "pybmad/generated/SimUtils_routines_c.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_c(py::module& m) {
  m.def(
      "calc_file_number",
      &SimUtils::calc_file_number,
      py::arg("file_name"),
      py::arg("num_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      R"""(Parameters
----------
file_name : 
num_in : 
num_out : 
err_flag : 
)""");
  m.def(
      "change_file_number",
      &SimUtils::change_file_number,
      py::arg("file_name"),
      py::arg("change"),
      R"""(Parameters
----------
file_name : 
change : 
)""");
  m.def(
      "charge_of",
      &SimUtils::charge_of,
      py::arg("species"),
      py::arg("default_") = py::none(),
      R"""(Function charge_of (species, default) result (charge)

Routine to return the charge, in units of e+, of a particle.

Parameters
----------
species : int
    Species ID.
default : int, optional
    If present then use default value if species = not_set$.

Returns
-------
charge : int
    particle charge.
)""");
  m.def(
      "charge_to_mass_of",
      &SimUtils::charge_to_mass_of,
      py::arg("species"),
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
)""");
  m.def(
      "coarse_frequency_estimate",
      &SimUtils::coarse_frequency_estimate,
      py::arg("data"),
      py::arg("error") = py::none(),
      R"""(Function coarse_frequency_estimate(data, error) result(frequency)

Simple function to take periodic data and estimate
the most dominant frequency by FFT.

Parameters
----------
data : float
    data to analyze. Preferably size(data) is a power of 2 Otherwise the data is padded with zeros.

Returns
-------
frequency : float
    Frequency corresponding to the largest FFT amplitude
err : bool
    Error: not enough data. Frequency is near 0 or 0.5
)""");
  m.def(
      "complex_error_function",
      &SimUtils::complex_error_function,
      py::arg("wr"),
      py::arg("wi"),
      py::arg("zr"),
      py::arg("zi"),
      R"""(Parameters
----------
wr : 
wi : 
zr : 
zi : 
)""");
  m.def(
      "cos_one",
      &SimUtils::cos_one,
      py::arg("angle"),
      py::arg("cos1"),
      R"""(Parameters
----------
angle : 
cos1 : 
)""");
  m.def(
      "cosc",
      &SimUtils::cosc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("y"),
      R"""(Parameters
----------
x : 
nd : 
y : 
)""");
  m.def(
      "create_a_spline",
      &SimUtils::create_a_spline,
      py::arg("r0"),
      py::arg("r1"),
      py::arg("slope0"),
      py::arg("slope1"),
      R"""(Function create_a_spline (r0, r1, slope0, slope1) result (spline)

Routine to create a single spline given end point positions and slopes.
The spline will pass through the data points and have the given slopes
at these points.

Modules used:
  use spline_mod

Parameters
----------
r0 : float
    Start (x, y) point.
r1 : float
    End (x, y) point.
slope0 : float
    Starting slope.
slope1 : float
    End slope.

Returns
-------
spline : SplineStruct
    Spline.
)""");
  m.def(
      "cross_product",
      &SimUtils::cross_product,
      py::arg("a"),
      py::arg("b"),
      py::arg("c"),
      R"""(Parameters
----------
a : float
    Input vectors.
b : 
c : 
)""");
}
