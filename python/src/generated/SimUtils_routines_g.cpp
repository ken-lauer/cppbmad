#include "pybmad/generated/SimUtils_routines_g.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_g(py::module &m) {
  m.def(
      "gen_complete_elliptic",
      &SimUtils::gen_complete_elliptic,
      py::arg("kc"),
      py::arg("p"),
      py::arg("c"),
      py::arg("s"),
      py::arg("err_tol") = py::none(),
      py::arg("value"),
      R"""(Wrapper for Fortran routine gen_complete_elliptic

Parameters
----------
kc : float

p : float

c : float

s : float

err_tol : float, optional

value : float

Returns
-------
kc : float

p : float

c : float

s : float

err_tol : float, optional

value : float
)"""
  );
  m.def(
      "get_file_number",
      &SimUtils::get_file_number,
      py::arg("file_name"),
      py::arg("cnum_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      R"""(Wrapper for Fortran routine get_file_number

Parameters
----------
file_name : character

cnum_in : character

num_out : int

err_flag : bool

Returns
-------
file_name : character

cnum_in : character

num_out : int

err_flag : bool
)"""
  );
  m.def(
      "get_file_time_stamp",
      &SimUtils::get_file_time_stamp,
      py::arg("file"),
      py::arg("time_stamp"),
      R"""(no longer exists
subroutine get_next_number (filein, cnum, digits)
  implicit none
  character(*) filein
  character(*) cnum
  integer digits
end subroutine
)"""
  );
}
