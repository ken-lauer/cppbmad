#include "pybmad/generated/SimUtils_routines_g.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_g(nb::module_ &m) {
  nb::class_<SimUtils::Gelbd>(m, "Gelbd", "gelbd return type")
      .def_ro("elb", &SimUtils::Gelbd::elb)
      .def_ro("eld", &SimUtils::Gelbd::eld)
      .def("__len__", [](const SimUtils::Gelbd &) { return 2; })
      .def("__getitem__", [](const SimUtils::Gelbd &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.elb);
        if (i == 1)
          return nb::cast(s.eld);
        throw nb::index_error();
      });
  m.def(
      "gelbd",
      &SimUtils::gelbd,
      nb::arg("phi"),
      nb::arg("mc"),
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
      nb::arg("kc"),
      nb::arg("p"),
      nb::arg("c"),
      nb::arg("s"),
      nb::arg("err_tol") = nb::none(),
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
      "general_bin_count",
      &SimUtils::general_bin_count,
      nb::arg("bin_data"),
      nb::arg("ix1"),
      nb::arg("ix2") = nb::none(),
      nb::arg("ix3") = nb::none(),
      R"""(Function for getting the count at a general index. Count will be 0 if out of bounds
)"""
  );
  m.def(
      "general_bin_index",
      &SimUtils::general_bin_index,
      nb::arg("bin_data"),
      nb::arg("ix1"),
      nb::arg("ix2") = nb::none(),
      nb::arg("ix3") = nb::none(),
      R"""(Function for looking up an index in the 1D count array
)"""
  );
  m.def(
      "general_bin_index_in_bounds",
      &SimUtils::general_bin_index_in_bounds,
      nb::arg("bin_data"),
      nb::arg("ix1"),
      nb::arg("ix2") = nb::none(),
      nb::arg("ix3") = nb::none(),
      R"""(Function for checking bounds
)"""
  );
  m.def(
      "get_a_char",
      [](bool wait, CharacterAlloc1D *ignore_this) {
        auto fn =
            static_cast<std::string (*)(bool, optional_ref<CharacterAlloc1D>)>(&SimUtils::get_a_char
            );
        return fn(wait, ptr_to_opt_ref(ignore_this));
      },
      nb::arg("wait"),
      nb::arg("ignore_this") = nb::none(),
      R"""(Subroutine for getting a single character from the terminal.
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
      nb::arg("file_name"),
      nb::arg("cnum_in"),
      nb::arg("num_out"),
      nb::arg("err_flag"),
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
      nb::arg("file"),
      nb::arg("time_stamp"),
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
      nb::arg("wait"),
      nb::arg("flush"),
      R"""(Subroutine for getting a single character from the terminal.
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
