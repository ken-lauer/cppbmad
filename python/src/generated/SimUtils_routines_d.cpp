#include "pybmad/generated/SimUtils_routines_d.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_d(py::module& m) {
  m.def(
      "date_and_time_stamp",
      &SimUtils::date_and_time_stamp,
      py::arg("string"),
      py::arg("numeric_month") = py::none(),
      py::arg("include_zone") = py::none(),
      R"""(Parameters
----------
string : 
numeric_month : 
include_zone : 
)""");
  m.def(
      "destfixedwindowls",
      &SimUtils::destfixedwindowls,
      py::arg("id"),
      R"""(Parameters
----------
id : 
)""");
  m.def(
      "detab",
      &SimUtils::detab,
      py::arg("str"),
      R"""(Parameters
----------
str : 
)""");
  m.def(
      "display_size_and_resolution",
      &SimUtils::display_size_and_resolution,
      py::arg("ix_screen"),
      py::arg("x_size"),
      py::arg("y_size"),
      py::arg("x_res"),
      py::arg("y_res"),
      R"""(Parameters
----------
ix_screen : 
x_size : 
y_size : 
x_res : 
y_res : 
)""");
  m.def(
      "dj_bessel",
      &SimUtils::dj_bessel,
      py::arg("m"),
      py::arg("arg"),
      py::arg("dj_bes"),
      R"""(Parameters
----------
m : 
arg : 
dj_bes : 
)""");
  m.def(
      "djb_hash",
      &SimUtils::djb_hash,
      py::arg("str"),
      py::arg("old_hash") = py::none(),
      py::arg("hash"),
      R"""(Parameters
----------
str : 
old_hash : 
hash : 
)""");
  m.def(
      "djb_str_hash",
      &SimUtils::djb_str_hash,
      py::arg("in_str"),
      py::arg("hash_str"),
      R"""(Parameters
----------
in_str : 
hash_str : 
)""");
  m.def(
      "downcase_string",
      &SimUtils::downcase_string,
      py::arg("string"),
      R"""(Parameters
----------
string : 
)""");
}
