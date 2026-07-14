#include "pybmad/generated/SimUtils_routines_i.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_i(nb::module_ &m) {
  m.def(
      "i_bessel",
      &SimUtils::i_bessel,
      nb::arg("m"),
      nb::arg("arg"),
      R"""(Wrapper for Fortran routine i_bessel

Parameters
----------
m : int
    Bessel order.

arg : float
    Bessel argument.

Returns
-------
i_bes : float
    Bessel value.
)"""
  );
  m.def(
      "i_bessel_extended",
      &SimUtils::i_bessel_extended,
      nb::arg("m"),
      nb::arg("arg"),
      R"""(Wrapper for Fortran routine i_bessel_extended

Parameters
----------
m : int
    Bessel order.

arg : float
    Bessel argument.

Returns
-------
i_bes : complex
    Bessel value.
)"""
  );
  m.def(
      "increment_file_number",
      &SimUtils::increment_file_number,
      nb::arg("file_name"),
      nb::arg("digits"),
      nb::arg("number"),
      nb::arg("cnumber"),
      R"""(Wrapper for Fortran routine increment_file_number

Parameters
----------
file_name : str

digits : int

number : int

cnumber : str
)"""
  );
  m.def(
      "index_nocase",
      &SimUtils::index_nocase,
      nb::arg("string1"),
      nb::arg("string2"),
      R"""(Wrapper for Fortran routine index_nocase

Parameters
----------
string1 : str

string2 : str

Returns
-------
indx : int
)"""
  );
  m.def(
      "initfixedwindowls",
      &SimUtils::initfixedwindowls,
      nb::arg("N"),
      nb::arg("dt"),
      nb::arg("order"),
      nb::arg("der"),
      R"""(Initializes an instance of the fixed window least squares module.
See module documentation (getf windowLS_mod) for use details.
Any instance of windowLS created with this module should be destroyed with destFixedWindowLS.

Parameters
----------
N : int
    Number of data points to fit over. aka window size.

dt : float
    Time interval between data points. It is assumed that the data is separated by fixed time intervals.

order : int
    Order of fit polynomial.  Must be greater than or equal to der.

der : int
    Order of derivative to be returned. Set der=0 to obtain the fit.
)"""
  );
  m.def(
      "initial_lmdif",
      &SimUtils::initial_lmdif,
      R"""(Wrapper for Fortran routine initial_lmdif
)"""
  );
  m.def(
      "int_str",
      &SimUtils::int_str,
      nb::arg("int_"),
      nb::arg("width") = nb::none(),
      R"""(Wrapper for Fortran routine int_str

Parameters
----------
width : int, optional

Returns
-------
str : str
)"""
  );
  m.def(
      "interpolated_fft",
      &SimUtils::interpolated_fft,
      nb::arg("cdata"),
      nb::arg("calc_ok"),
      nb::arg("opt_dump_spectrum") = nb::none(),
      nb::arg("opt_dump_index") = nb::none(),
      R"""(Windows the complex data and used Numerical Recipes four1 to find the peak in the spectrum.
The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
available.
)"""
  );
  m.def(
      "interpolated_fft_gsl",
      &SimUtils::interpolated_fft_gsl,
      nb::arg("cdata"),
      nb::arg("calc_ok"),
      nb::arg("opt_dump_spectrum") = nb::none(),
      nb::arg("opt_dump_index") = nb::none(),
      R"""(Windows the complex data and uses a mixed-radix GSL routine to find the peak in the spectrum.
The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
available.
)"""
  );
  m.def(
      "inverse_prob",
      &SimUtils::inverse_prob,
      nb::arg("val"),
      R"""(Wrapper for Fortran routine inverse_prob

Parameters
----------
val : float

Returns
-------
prob : float
)"""
  );
  m.def(
      "is_alphabetic",
      &SimUtils::is_alphabetic,
      nb::arg("string"),
      nb::arg("valid_chars") = nb::none(),
      R"""(Wrapper for Fortran routine is_alphabetic

Parameters
----------
string : str

valid_chars : str, optional

Returns
-------
is_alpha : bool
)"""
  );
  m.def(
      "is_decreasing_sequence",
      &SimUtils::is_decreasing_sequence,
      nb::arg("array"),
      nb::arg("strict") = nb::none(),
      R"""(Wrapper for Fortran routine is_decreasing_sequence

Parameters
----------
array : 1D array of float
    Sequence.

strict : bool, optional
    If True (default) sequence must be strictly decreasing.

Returns
-------
is_decreasing : bool
    Set True if sequence is decreasing.
)"""
  );
  m.def(
      "is_false",
      &SimUtils::is_false,
      nb::arg("param"),
      R"""(Routine to translate from a real number to a boolian True or False.
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
)"""
  );
  m.def(
      "is_increasing_sequence",
      &SimUtils::is_increasing_sequence,
      nb::arg("array"),
      nb::arg("strict") = nb::none(),
      R"""(Wrapper for Fortran routine is_increasing_sequence

Parameters
----------
array : 1D array of float
    Sequence.

strict : bool, optional
    If True (default) sequence must be strictly increasing.

Returns
-------
is_increasing : bool
    Set True if sequence is increasing.
)"""
  );
  m.def(
      "is_integer",
      &SimUtils::is_integer,
      nb::arg("string"),
      nb::arg("int_") = nb::none(),
      nb::arg("delims") = nb::none(),
      nb::arg("ix_word") = nb::none(),
      R"""(Wrapper for Fortran routine is_integer

Parameters
----------
string : str

delims : str, optional

ix_word : int, optional

Returns
-------
valid : bool
)"""
  );
  m.def(
      "is_logical",
      &SimUtils::is_logical,
      nb::arg("string"),
      nb::arg("ignore") = nb::none(),
      R"""(Wrapper for Fortran routine is_logical

Parameters
----------
string : str

ignore : bool, optional

Returns
-------
valid : bool
)"""
  );
  m.def(
      "is_real",
      &SimUtils::is_real,
      nb::arg("string"),
      nb::arg("ignore") = nb::none(),
      nb::arg("real_num") = nb::none(),
      R"""(Wrapper for Fortran routine is_real

Parameters
----------
string : str

ignore : bool, optional

real_num : float, optional

Returns
-------
valid : bool
)"""
  );
  m.def(
      "is_subatomic_species",
      &SimUtils::is_subatomic_species,
      nb::arg("species"),
      R"""(Routine to return True if species argument corresponds to a subatomic particle.

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
      nb::arg("param"),
      R"""(Routine to translate from a real number to a boolian True or False.
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
)"""
  );
}
