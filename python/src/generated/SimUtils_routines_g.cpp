#include "pybmad/generated/SimUtils_routines_g.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_g(py::module& m) {
  m.def(
      "gen_complete_elliptic",
      &SimUtils::gen_complete_elliptic,
      py::arg("kc"),
      py::arg("p"),
      py::arg("c"),
      py::arg("s"),
      py::arg("err_tol") = py::none(),
      py::arg("value"),
      R"""(Parameters
----------
kc : 
p : 
c : 
s : 
err_tol : 
value : 
)""");
  m.def(
      "get_file_number",
      &SimUtils::get_file_number,
      py::arg("file_name"),
      py::arg("cnum_in"),
      py::arg("num_out"),
      py::arg("err_flag"),
      R"""(Parameters
----------
file_name : 
cnum_in : 
num_out : 
err_flag : 
)""");
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

)""");
}
