#include "pybmad/generated/SimUtils_routines_n.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_n(nb::module_ &m) {
  m.def(
      "n_bins_automatic",
      &SimUtils::n_bins_automatic,
      nb::arg("n_data"),
      R"""(Function to automatically select the number of bins
)"""
  );
  m.def(
      "n_choose_k",
      &SimUtils::n_choose_k,
      nb::arg("n"),
      nb::arg("k"),
      R"""(Wrapper for Fortran routine n_choose_k

Parameters
----------
n : int
    Must be non-negative with n >= k.

k : int
    Must be non-negative with n >= k.

Returns
-------
nck : float
    N choose K will return negative number if there is an error.
)"""
  );
  m.def(
      "n_spline_create",
      &SimUtils::n_spline_create,
      nb::arg("deriv0"),
      nb::arg("deriv1"),
      nb::arg("x1"),
      nb::arg("n_spline"),
      R"""(Wrapper for Fortran routine n_spline_create

Parameters
----------
deriv0 : 1D array of float
    Derivative vector from order 0 to some order n at x = 0.

deriv1 : 1D array of float
    Derivative vector from order 0 to some order n at x = x1.

x1 : float
    Location where deriv1 derivatives have been evaluated.

n_spline : 1D array of float
    real(rp), Derivative vector from order 0 to order 2*n+1 of the interpolation spline.
)"""
  );
  m.def(
      "naff",
      &SimUtils::naff,
      nb::arg("cdata"),
      nb::arg("freqs"),
      nb::arg("amps"),
      nb::arg("opt_dump_spectra") = nb::none(),
      nb::arg("opt_zero_first") = nb::none(),
      R"""(This subroutine implements the NAFF algorithm for calculating the spectra
of periodic data.

See naff_mod documentation for details.

Frequencies returned are in units of 2pi. That is, freqs ranges from 0 to 1.

freqs and amps must be allocated before hand.  This subroutine will repeat the
decomposition loop until all elements of freqs and amps are populated.
)"""
  );
  m.def(
      "nametable_add",
      &SimUtils::nametable_add,
      nb::arg("nametable"),
      nb::arg("name"),
      nb::arg("ix_name"),
      R"""(Wrapper for Fortran routine nametable_add

Parameters
----------
nametable : NametableStruct

name : str

ix_name : int
)"""
  );
  m.def(
      "nametable_bracket_indexx",
      &SimUtils::nametable_bracket_indexx,
      nb::arg("nametable"),
      nb::arg("name"),
      nb::arg("n_match") = nb::none(),
      R"""(Wrapper for Fortran routine nametable_bracket_indexx

Parameters
----------
nametable : NametableStruct

name : str

n_match : int, optional

Returns
-------
ix_max : int
)"""
  );
  m.def(
      "nametable_change1",
      &SimUtils::nametable_change1,
      nb::arg("nametable"),
      nb::arg("name"),
      nb::arg("ix_name"),
      R"""(Wrapper for Fortran routine nametable_change1

Parameters
----------
nametable : NametableStruct

name : str

ix_name : int
)"""
  );
  m.def(
      "nametable_init",
      &SimUtils::nametable_init,
      nb::arg("nametable"),
      nb::arg("n_min") = nb::none(),
      nb::arg("n_max") = nb::none(),
      R"""(Wrapper for Fortran routine nametable_init

Parameters
----------
nametable : NametableStruct

n_min : int, optional

n_max : int, optional
)"""
  );
  m.def(
      "nametable_remove",
      &SimUtils::nametable_remove,
      nb::arg("nametable"),
      nb::arg("ix_name"),
      R"""(Wrapper for Fortran routine nametable_remove

Parameters
----------
nametable : NametableStruct

ix_name : int
)"""
  );
  m.def(
      "negative_ampsquared",
      &SimUtils::negative_ampsquared,
      nb::arg("frequency"),
      nb::arg("status") = nb::none(),
      R"""(Wrapper for Fortran routine negative_ampsquared

Parameters
----------
frequency : float

status : int, optional

Returns
-------
amp : float
)"""
  );
  m.def(
      "negative_dampsquared",
      &SimUtils::negative_dampsquared,
      nb::arg("frequency"),
      nb::arg("status") = nb::none(),
      R"""(Wrapper for Fortran routine negative_dampsquared

Parameters
----------
frequency : float

status : int, optional

Returns
-------
damp : float
)"""
  );
}
