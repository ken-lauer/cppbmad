#include "pybmad/generated/SimUtils_routines_d.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_d(nb::module_ &m) {
  m.def(
      "date_and_time_stamp",
      &SimUtils::date_and_time_stamp,
      nb::arg("string"),
      nb::arg("numeric_month") = nb::none(),
      nb::arg("include_zone") = nb::none(),
      R"""(Wrapper for Fortran routine date_and_time_stamp

Parameters
----------
string : str

numeric_month : bool, optional

include_zone : bool, optional
)"""
  );
  m.def(
      "destfixedwindowls",
      &SimUtils::destfixedwindowls,
      nb::arg("id"),
      R"""(Wrapper for Fortran routine destfixedwindowls

Parameters
----------
id : int
)"""
  );
  m.def(
      "detab",
      &SimUtils::detab,
      nb::arg("str"),
      R"""(Wrapper for Fortran routine detab

Parameters
----------
str : str
)"""
  );
  m.def(
      "display_size_and_resolution",
      &SimUtils::display_size_and_resolution,
      nb::arg("ix_screen"),
      nb::arg("x_size"),
      nb::arg("y_size"),
      nb::arg("x_res"),
      nb::arg("y_res"),
      R"""(Wrapper for Fortran routine display_size_and_resolution

Parameters
----------
ix_screen : int

x_size : float

y_size : float

x_res : float

y_res : float
)"""
  );
  m.def(
      "dj_bessel",
      &SimUtils::dj_bessel,
      nb::arg("m"),
      nb::arg("arg"),
      R"""(Wrapper for Fortran routine dj_bessel

Parameters
----------
m : int

arg : float
    Bessel argument.

Returns
-------
dj_bes : float
    Bessel value.
)"""
  );
  m.def(
      "djb_hash",
      &SimUtils::djb_hash,
      nb::arg("str"),
      nb::arg("old_hash") = nb::none(),
      R"""(Wrapper for Fortran routine djb_hash

Parameters
----------
str : str

old_hash : int, optional

Returns
-------
hash : int
)"""
  );
  m.def(
      "djb_str_hash",
      &SimUtils::djb_str_hash,
      nb::arg("in_str"),
      R"""(Wrapper for Fortran routine djb_str_hash

Parameters
----------
in_str : str

Returns
-------
hash_str : str
)"""
  );
  m.def(
      "downcase_string",
      &SimUtils::downcase_string,
      nb::arg("string"),
      R"""(Wrapper for Fortran routine downcase_string

Parameters
----------
string : str
)"""
  );
}
