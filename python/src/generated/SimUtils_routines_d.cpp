#include "pybmad/generated/SimUtils_routines_d.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_d(py::module &m) {
  m.def(
      "date_and_time_stamp",
      &SimUtils::date_and_time_stamp,
      py::arg("string"),
      py::arg("numeric_month") = py::none(),
      py::arg("include_zone") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("id"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine destfixedwindowls

Parameters
----------
id : int
)"""
  );
  m.def(
      "detab",
      &SimUtils::detab,
      py::arg("str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine detab

Parameters
----------
str : str
)"""
  );
  m.def(
      "display_size_and_resolution",
      &SimUtils::display_size_and_resolution,
      py::arg("ix_screen"),
      py::arg("x_size"),
      py::arg("y_size"),
      py::arg("x_res"),
      py::arg("y_res"),
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("m"),
      py::arg("arg"),
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("str"),
      py::arg("hash"),
      py::arg("old_hash") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine djb_hash

Parameters
----------
str : str

hash : int

old_hash : int, optional
)"""
  );
  m.def(
      "djb_str_hash",
      &SimUtils::djb_str_hash,
      py::arg("in_str"),
      py::arg("hash_str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine djb_str_hash

Parameters
----------
in_str : str

hash_str : str
)"""
  );
  m.def(
      "downcase_string",
      &SimUtils::downcase_string,
      py::arg("string"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine downcase_string

Parameters
----------
string : str
)"""
  );
}
