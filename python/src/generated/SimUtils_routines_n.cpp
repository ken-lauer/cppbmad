#include "pybmad/generated/SimUtils_routines_n.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_n(py::module& m) {
  m.def(
      "n_bins_automatic",
      &SimUtils::n_bins_automatic,
      py::arg("n_data"),
      py::arg("n"),
      R"""(Function to automatically select the number of bins

)""");
  m.def(
      "n_choose_k",
      &SimUtils::n_choose_k,
      py::arg("n"),
      py::arg("k"),
      py::arg("nck"),
      R"""(Parameters
----------
n : 
k : 
nck : 
)""");
  m.def(
      "n_spline_create",
      &SimUtils::n_spline_create,
      py::arg("deriv0"),
      py::arg("deriv1"),
      py::arg("x1"),
      py::arg("n_spline"),
      R"""(Parameters
----------
deriv0 : float
    Derivative vector from order 0 to some order n at x = 0.
deriv1 : float
    Derivative vector from order 0 to some order n at x = x1.
x1 : float
    Location where deriv1 derivatives have been evaluated.
n_spline : 
    real(rp), Derivative vector from order 0 to order 2*n+1 of the interpolation spline.
)""");
  m.def(
      "naff",
      &SimUtils::naff,
      py::arg("cdata"),
      py::arg("freqs"),
      py::arg("amps"),
      py::arg("opt_dump_spectra") = py::none(),
      py::arg("opt_zero_first") = py::none(),
      R"""(subroutine naff(cdata,freqs,amps,opt_dump_spectra,opt_zero_first)

This subroutine implements the NAFF algorithm for calculating the spectra
of periodic data.

See naff_mod documentation for details.

Frequencies returned are in units of 2pi. That is, freqs ranges from 0 to 1.

freqs and amps must be allocated before hand.  This subroutine will repeat the
decomposition loop until all elements of freqs and amps are populated.

)""");
  m.def(
      "nametable_add",
      &SimUtils::nametable_add,
      py::arg("nametable"),
      py::arg("name"),
      py::arg("ix_name"),
      R"""(Parameters
----------
nametable : 
name : 
ix_name : 
)""");
  m.def(
      "nametable_bracket_indexx",
      &SimUtils::nametable_bracket_indexx,
      py::arg("nametable"),
      py::arg("name"),
      py::arg("n_match") = py::none(),
      py::arg("ix_max"),
      R"""(Parameters
----------
nametable : 
name : 
n_match : 
ix_max : 
)""");
  m.def(
      "nametable_change1",
      &SimUtils::nametable_change1,
      py::arg("nametable"),
      py::arg("name"),
      py::arg("ix_name"),
      R"""(Parameters
----------
nametable : 
name : 
ix_name : 
)""");
  m.def(
      "nametable_init",
      &SimUtils::nametable_init,
      py::arg("nametable"),
      py::arg("n_min") = py::none(),
      py::arg("n_max") = py::none(),
      R"""(Parameters
----------
nametable : 
n_min : 
n_max : 
)""");
  m.def(
      "nametable_remove",
      &SimUtils::nametable_remove,
      py::arg("nametable"),
      py::arg("ix_name"),
      R"""(Parameters
----------
nametable : 
ix_name : 
)""");
}
