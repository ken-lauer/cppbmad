#include "pybmad/generated/SimUtils_routines_i.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_i(py::module &m) {
  m.def(
      "i_bessel",
      &SimUtils::i_bessel,
      py::arg("m"),
      py::arg("arg"),
      py::arg("i_bes"),
      R"""(Parameters
----------
m : 
arg : 
i_bes : 
)"""
  );
  m.def(
      "i_bessel_extended",
      &SimUtils::i_bessel_extended,
      py::arg("m"),
      py::arg("arg"),
      py::arg("i_bes"),
      R"""(Parameters
----------
m : 
arg : 
i_bes : 
)"""
  );
  m.def(
      "increment_file_number",
      &SimUtils::increment_file_number,
      py::arg("file_name"),
      py::arg("digits"),
      py::arg("number"),
      py::arg("cnumber"),
      R"""(Parameters
----------
file_name : 
digits : 
number : 
cnumber : 
)"""
  );
  m.def(
      "index_nocase",
      &SimUtils::index_nocase,
      py::arg("string1"),
      py::arg("string2"),
      py::arg("indx"),
      R"""(Parameters
----------
string1 : 
string2 : 
indx : 
)"""
  );
  m.def(
      "initfixedwindowls",
      &SimUtils::initfixedwindowls,
      py::arg("N"),
      py::arg("dt"),
      py::arg("order"),
      py::arg("der"),
      py::arg("id"),
      R"""(Function initFixedWindowLS

Initializes an instance of the fixed window least squares module.
See module documentation (getf windowLS_mod) for use details.
Any instance of windowLS created with this module should be destroyed with destFixedWindowLS.

Parameters
----------
N : int
    Number of data points to fit over. aka window size.
dt : float
    Time interval between data points. It is assumed that the data is
separated by fixed time intervals. : 
order : int
    Order of fit polynomial.  Must be greater than or equal to der.
der : int
    Order of derivative to be returned. Set der=0 to obtain the fit.

Returns
-------
<return value> : int
    id of windowLS instance created.
)"""
  );
  m.def(
      "int_str",
      &SimUtils::int_str,
      py::arg("int_"),
      py::arg("width") = py::none(),
      py::arg("str"),
      R"""(Parameters
----------
int : 
width : 
str : 
)"""
  );
  m.def(
      "interpolated_fft",
      &SimUtils::interpolated_fft,
      py::arg("cdata"),
      py::arg("calc_ok"),
      py::arg("opt_dump_spectrum") = py::none(),
      py::arg("opt_dump_index") = py::none(),
      py::arg("this_fft"),
      R"""(Function interpolated_fft (cdata, calc_ok, opt_dump_spectrum, opt_dump_index) result (this_fft)

Windows the complex data and used Numerical Recipes four1 to find the peak in the spectrum.
The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
available.


Returns
-------
this_fft
)"""
  );
  m.def(
      "interpolated_fft_gsl",
      &SimUtils::interpolated_fft_gsl,
      py::arg("cdata"),
      py::arg("calc_ok"),
      py::arg("opt_dump_spectrum") = py::none(),
      py::arg("opt_dump_index") = py::none(),
      py::arg("this_fft"),
      R"""(function interpolated_fft_gsl

Windows the complex data and uses a mixed-radix GSL routine to find the peak in the spectrum.
The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
available.

)"""
  );
  m.def(
      "is_alphabetic",
      &SimUtils::is_alphabetic,
      py::arg("string"),
      py::arg("valid_chars") = py::none(),
      py::arg("is_alpha"),
      R"""(no longer exists
function inverse_prob (val) result (prob)
  import
  implicit none
  real(rp) prob
  real(rp) val
end function


Returns
-------
prob
)"""
  );
  m.def(
      "is_decreasing_sequence",
      &SimUtils::is_decreasing_sequence,
      py::arg("array"),
      py::arg("strict") = py::none(),
      R"""(Parameters
----------
array : float
    Sequence.
strict : bool, optional
    If True (default) sequence must be strictly decreasing.
is_decreasing : bool
    Set True if sequence is decreasing.
)"""
  );
  m.def(
      "is_false",
      &SimUtils::is_false,
      py::arg("param"),
      R"""(Function is_false (param) result (this_false)

Routine to translate from a real number to a boolian True or False.
Translation: 0 = False, nonzero = True

Also see: is_true and int_logic

The typical use of this routine is for parameters in ele_struct%value(:) which
is a real array. Some of the elements in the %value array are used to specify
boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).

Parameters
----------
param : float
    Real number to be translated

Returns
-------
this_false : bool
    Set True if param is zero. False otherwise.

Notes
-----
Related routines:
is_true int_logic ) which is a real array. Some of the elements in the %value array are used to specify
boolian attributes. For example quadrupoles use ele%value(scale_multipoles$).
)"""
  );
  m.def(
      "is_increasing_sequence",
      &SimUtils::is_increasing_sequence,
      py::arg("array"),
      py::arg("strict") = py::none(),
      R"""(Parameters
----------
array : float
    Sequence.
strict : bool, optional
    If True (default) sequence must be strictly increasing.
is_increasing : bool
    Set True if sequence is increasing.
)"""
  );
  m.def(
      "is_integer",
      &SimUtils::is_integer,
      py::arg("string"),
      py::arg("int_") = py::none(),
      py::arg("delims") = py::none(),
      py::arg("ix_word") = py::none(),
      py::arg("valid"),
      R"""(Parameters
----------
string : 
int : 
delims : 
ix_word : 
valid : 
)"""
  );
  m.def(
      "is_logical",
      &SimUtils::is_logical,
      py::arg("string"),
      py::arg("ignore") = py::none(),
      py::arg("valid"),
      R"""(Parameters
----------
string : 
ignore : 
valid : 
)"""
  );
  m.def(
      "is_real",
      &SimUtils::is_real,
      py::arg("string"),
      py::arg("ignore") = py::none(),
      py::arg("real_num") = py::none(),
      py::arg("valid"),
      R"""(Parameters
----------
string : 
ignore : 
real_num : 
valid : 
)"""
  );
  m.def(
      "is_subatomic_species",
      &SimUtils::is_subatomic_species,
      py::arg("species"),
      R"""(Function is_subatomic_species(species) result (is_subatomic)

Routine to return True if species argument corresponds to a subatomic particle.

Parameters
----------
species : int
    Spicies ID.

Returns
-------
is_subatomic : bool
    Set True if species corresponds to a subatomic particle.
)"""
  );
  m.def(
      "is_true",
      &SimUtils::is_true,
      py::arg("param"),
      R"""(Function is_true (param) result (this_true)

Routine to translate from a real number to a boolian True or False.
Translation: 0 = False, nonzero = True

Also see: is_false and int_logic

The typical use of this routine is for parameters in ele_struct%value(:) which
is a real array. Some of the elements in the %value array are used to specify
boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).

Parameters
----------
param : float
    Real number to be translated

Returns
-------
this_true : bool
    Set False if param is zero. True otherwise.

Notes
-----
Related routines:
is_false int_logic ) which is a real array. Some of the elements in the %value array are used to specify
boolian attributes. For example quadrupoles use ele%value(scale_multipoles$).
)"""
  );
}
