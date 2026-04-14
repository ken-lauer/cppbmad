#include "pybmad/generated/SimUtils_routines_g.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_g(py::module &m) {
  py::class_<SimUtils::Gelbd, std::unique_ptr<SimUtils::Gelbd>>(m, "Gelbd", "gelbd return type")
      .def_readonly("elb", &SimUtils::Gelbd::elb)
      .def_readonly("eld", &SimUtils::Gelbd::eld)
      .def("__len__", [](const SimUtils::Gelbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Gelbd &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.elb);
        if (i == 1)
          return py::cast(s.eld);
        throw py::index_error();
      });
  m.def(
      "gelbd",
      &SimUtils::gelbd,
      py::arg("phi"),
      py::arg("mc"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine gelbd

Parameters
----------
phi : float

mc : float

Returns
-------
elb : float

eld : float
)"""
  );
  m.def(
      "gen_complete_elliptic",
      &SimUtils::gen_complete_elliptic,
      py::arg("kc"),
      py::arg("p"),
      py::arg("c"),
      py::arg("s"),
      py::arg("err_tol") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine gen_complete_elliptic

Parameters
----------
kc : float
    Fuction input values.

p : float
    Fuction input values.

c : float
    Fuction input values.

s : float
    Fuction input values.

err_tol : float, optional
    Relative error tolerance. Default = 1d-12

Returns
-------
value : float
    Output value.
)"""
  );
  m.def(
      "get_a_char",
      &SimUtils::get_a_char,
      py::arg("wait"),
      py::arg("ignore_this") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine get_a_char (this_char, wait, ignore_this)

Subroutine for getting a single character from the terminal.
Also see: get_tty_char

System Libraries that need to be linked to:
  readline curses

Parameters
----------
wait : bool
    If True then routine will wait until a keystroke has occured. If False and no keystroke is in the buffer
    then achar(0) will be returned as this_char.

ignore_this : 1D array of str, optional
    List of characters to ignore. If a keystroke matches a character on this list the keystroke is ignored.

Returns
-------
this_char : str
    Character returned
)"""
  );
  m.def(
      "get_file_number",
      &SimUtils::get_file_number,
      py::arg("file_name"),
      py::arg("cnum_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine get_file_number

Parameters
----------
file_name : str

cnum_in : str

num_out : int

err_flag : bool
)"""
  );
  m.def(
      "get_file_time_stamp",
      &SimUtils::get_file_time_stamp,
      py::arg("file"),
      py::arg("time_stamp"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(no longer exists
subroutine get_next_number (filein, cnum, digits)
  implicit none
  character(*) filein
  character(*) cnum
  integer digits
end subroutine
)"""
  );
  m.def(
      "get_tty_char",
      &SimUtils::get_tty_char,
      py::arg("wait"),
      py::arg("flush"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine get_tty_char (this_char, wait, flush)

Subroutine for getting a single character from the terminal.
Also see: get_a_char

System Libraries that need to be linked to:
  readline curses

Parameters
----------
wait : bool
    If True then routine will wait until a keystroke has occured. If False and no keystroke is in the buffer
    then achar(0) will be returned as this_char.

flush : bool
    If True then the keystroke buffer will be cleared first before any processing.

Returns
-------
this_char : str
    Character returned
)"""
  );
}
